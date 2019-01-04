﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.Analyzers.Interfaces;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Discretization.Providers;
using ISAAR.MSolve.FEM.Providers;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.Solvers.Assemblers;
using ISAAR.MSolve.Solvers.Commons;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.Problems
{
    public class ProblemPorous_v2 : IImplicitIntegrationProvider_v2, IStaticProvider_v2, INonLinearProvider_v2
    {
        private bool providersInitialized = false;
        private double scalingCoefficient;
        private Dictionary<int, IMatrix> ms, cs, ks;
        private Dictionary<int, CsrMatrix> qs;
        private readonly Model_v2 model;
        private readonly ISolver_v2 solver;
        private readonly IReadOnlyDictionary<int, ILinearSystem_v2> linearSystems;
        private ElementPoreStiffnessProvider_v2 stiffnessProvider;
        private ElementPoreDampingProvider_v2 dampingProvider;
        private ElementPoreMassProvider_v2 massProvider;

        public ProblemPorous_v2(Model_v2 model, ISolver_v2 solver)
        {
            this.model = model;
            this.solver = solver;
            this.linearSystems = solver.LinearSystems;
        }

        private IDictionary<int, IMatrix> Ms
        {
            get
            {
                if (ms == null) BuildMs();
                return ms;
            }
        }

        private IDictionary<int, IMatrix> Cs
        {
            get
            {
                if (cs == null) BuildCs();
                return cs;
            }
        }

        private IDictionary<int, IMatrix> Ks
        {
            get
            {
                if (ks == null)
                    BuildKs();
                else
                    RebuildKs();
                return ks;
            }
        }

        private void InitializeProvidersAndBuildQs(ImplicitIntegrationCoefficients c)
        {
            //if (providersInitialized) return;

            //providersInitialized = true;
            massProvider = new ElementPoreMassProvider_v2(new ElementStructuralMassProvider_v2(), c.Damping);
            dampingProvider = new ElementPoreDampingProvider_v2(new ElementStructuralDampingProvider_v2(), c.Damping);
            stiffnessProvider = new ElementPoreStiffnessProvider_v2(new ElementStructuralStiffnessProvider_v2(), c.Damping);
            BuildQs();
        }

        private void BuildKs()
        {
            ks = new Dictionary<int, IMatrix>(model.Subdomains.Count);
            foreach (ISubdomain_v2 subdomain in model.Subdomains)
            {
                ks.Add(subdomain.ID, solver.BuildGlobalMatrix(subdomain, stiffnessProvider));
            }
        }

        private void RebuildKs()
        {
            foreach (ISubdomain_v2 subdomain in model.Subdomains)
            {
                if (subdomain.MaterialsModified)
                {
                    ks[subdomain.ID] = solver.BuildGlobalMatrix(subdomain, stiffnessProvider);
                    subdomain.ResetMaterialsModifiedProperty();
                }
            }
        }

        private void BuildMs()
        {
            ms = new Dictionary<int, IMatrix>(model.Subdomains.Count);
            foreach (ISubdomain_v2 subdomain in model.Subdomains)
            {
                ms.Add(subdomain.ID, solver.BuildGlobalMatrix(subdomain, massProvider));
            }
        }

        private void BuildCs()
        {
            cs = new Dictionary<int, IMatrix>(model.Subdomains.Count);
            foreach (ISubdomain_v2 subdomain in model.Subdomains)
            {
                cs.Add(subdomain.ID, solver.BuildGlobalMatrix(subdomain, dampingProvider));
            }
        }

        private void BuildQs()
        {
            qs = new Dictionary<int, CsrMatrix>(model.SubdomainsDictionary.Count);
            foreach (Subdomain_v2 subdomain in model.SubdomainsDictionary.Values)
            {
                qs.Add(subdomain.ID, BuildQFromSubdomain(subdomain));
            }
        }

        //TODO: this should be done by an assembler class
        //TODO: make sure this is called whenever the ordering changes
        private CsrMatrix BuildQFromSubdomain(Subdomain_v2 subdomain) 
        {
            int numFreeDofs = subdomain.DofOrdering.NumFreeDofs;
            var qSubdomain = DokRowMajor.CreateEmpty(numFreeDofs, numFreeDofs);
            DofTable allDofs = subdomain.DofOrdering.FreeDofs;
            foreach (Element_v2 element in subdomain.Elements)
            {
                if (!(element.ElementType is IPorousFiniteElement_v2)) continue;

                var e = (IPorousFiniteElement_v2)element.ElementType;
                IMatrix q = e.CouplingMatrix(element);

                int iElementMatrixRow = 0;
                for (int i = 0; i < element.ElementType.DofEnumerator.GetDOFTypes(element).Count; i++)
                {
                    Node_v2 nodeRow = element.Nodes[i];
                    foreach (DOFType dofTypeRow in element.ElementType.DofEnumerator.GetDOFTypes(element)[i])
                    {
                        if (dofTypeRow != DOFType.Pore) continue;

                        int dofRow = allDofs[nodeRow, dofTypeRow];
                        int iElementMatrixColumn = 0;

                        for (int j = 0; j < element.ElementType.DofEnumerator.GetDOFTypes(element).Count; j++)
                        {
                            Node_v2 nodeColumn = element.Nodes[j];
                            foreach (DOFType dofTypeColumn in element.ElementType.DofEnumerator.GetDOFTypes(element)[j])
                            {
                                if (dofTypeColumn == DOFType.Pore) continue;

                                int dofColumn = allDofs[nodeColumn, dofTypeColumn];
                                qSubdomain.AddToEntry(dofColumn, dofRow, q[iElementMatrixRow, iElementMatrixColumn]);
                                iElementMatrixColumn++;
                            }
                        }
                        iElementMatrixRow++;
                    }
                }
            }

            return qSubdomain.BuildCsrMatrix(true);
        }

        private void ScaleSubdomainSolidVector(ILinearSystem_v2 linearSystem, IVector vector)
        {
            int id = linearSystem.Subdomain.ID;
            foreach ((INode node, DOFType dofType, int dofIdx) in model.SubdomainsDictionary[id].DofOrdering.FreeDofs)
            {
                if (dofType!= DOFType.Pore) vector.Set(dofIdx, vector[dofIdx] * this.scalingCoefficient);
            }
        }

        private void CalculateEffectiveMatrixInternal(ILinearSystem_v2 linearSystem, ImplicitIntegrationCoefficients coefficients)
        {
            int id = linearSystem.Subdomain.ID;
            IMatrix matrix = this.Ks[id];
            matrix.LinearCombinationIntoThis(coefficients.Stiffness, Ms[id], coefficients.Mass);
            matrix.AxpyIntoThis(Cs[id], coefficients.Damping);
            linearSystem.SetMatrix(this.Ks[id]);
        }

        #region IAnalyzerProvider Members
        public void Reset()
        {
            foreach (Subdomain_v2 subdomain in model.SubdomainsDictionary.Values)
                foreach (Element_v2 element in subdomain.Elements)
                    element.ElementType.ClearMaterialState();

            cs = null;
            ks = null;
            ms = null;
        }
        #endregion

        #region IImplicitIntegrationProvider Members

        public void CalculateEffectiveMatrix(ILinearSystem_v2 linearSystem, ImplicitIntegrationCoefficients coefficients)
        {
            InitializeProvidersAndBuildQs(coefficients);
            scalingCoefficient = coefficients.Damping;
            CalculateEffectiveMatrixInternal(linearSystem, coefficients);
        }

        public void ProcessRhs(ILinearSystem_v2 linearSystem, ImplicitIntegrationCoefficients coefficients)
        {
            scalingCoefficient = coefficients.Damping;
            ScaleSubdomainSolidVector(linearSystem, linearSystem.RhsVector);
        }

        public void GetRhsFromHistoryLoad(int timeStep)
        {
            foreach (Subdomain_v2 subdomain in model.Subdomains) subdomain.Forces.Clear();

            model.AssignLoads();
            model.AssignMassAccelerationHistoryLoads(timeStep);
        }

        public IDictionary<int, IVector> GetAccelerationsOfTimeStep(int timeStep)
        {
            var d = new Dictionary<int, IVector>();
            foreach (ILinearSystem_v2 linearSystem in linearSystems.Values)
            {
                d.Add(linearSystem.Subdomain.ID, linearSystem.CreateZeroVector());
            }

            if (model.MassAccelerationHistoryLoads.Count > 0)
            {
                List<MassAccelerationLoad> m = new List<MassAccelerationLoad>(model.MassAccelerationHistoryLoads.Count);
                foreach (IMassAccelerationHistoryLoad l in model.MassAccelerationHistoryLoads)
                    m.Add(new MassAccelerationLoad() { Amount = l[timeStep], DOF = l.DOF });

                foreach (ISubdomain_v2 subdomain in model.Subdomains)
                {
                    int[] subdomainToGlobalDofs = model.GlobalDofOrdering.MapFreeDofsSubdomainToGlobal(subdomain);
                    foreach ((INode node, DOFType dofType, int subdomainDofIdx) in subdomain.DofOrdering.FreeDofs)
                    {
                        int globalDofIdx = subdomainToGlobalDofs[subdomainDofIdx];
                        foreach (var l in m)
                        {
                            if (dofType == l.DOF) d[subdomain.ID].Set(globalDofIdx, l.Amount);
                        }
                    }
                }
            }

            //foreach (ElementMassAccelerationHistoryLoad load in model.ElementMassAccelerationHistoryLoads)
            //{
            //    MassAccelerationLoad hl = new MassAccelerationLoad() { Amount = load.HistoryLoad[timeStep] * 564000000, DOF = load.HistoryLoad.DOF };
            //    load.Element.Subdomain.AddLocalVectorToGlobal(load.Element,
            //        load.Element.ElementType.CalculateAccelerationForces(load.Element, (new MassAccelerationLoad[] { hl }).ToList()),
            //        load.Element.Subdomain.Forces);
            //}

            return d;
        }

        public IDictionary<int, IVector> GetVelocitiesOfTimeStep(int timeStep)
        {
            var d = new Dictionary<int, IVector>();

            foreach (ILinearSystem_v2 linearSystem in linearSystems.Values)
            {
                d.Add(linearSystem.Subdomain.ID, linearSystem.CreateZeroVector());
            }

            return d;
        }

        public IVector MassMatrixVectorProduct(ILinearSystem_v2 linearSystem, IVectorView lhsVector)
            => this.Ms[linearSystem.Subdomain.ID].Multiply(lhsVector);

        public IVector DampingMatrixVectorProduct(ILinearSystem_v2 linearSystem, IVectorView lhsVector)
        {
            int id = linearSystem.Subdomain.ID;
            IVector result = this.Cs[id].Multiply(lhsVector);
            result.SubtractIntoThis(qs[id].Multiply(lhsVector));
            return result;
        }

        #endregion

        #region IStaticProvider Members

        public void CalculateMatrix(ILinearSystem_v2 linearSystem)
        {
            throw new NotImplementedException();
        }

        //public void CalculateMatrices()
        //{
        //    throw new NotImplementedException();
        //}

        #endregion

        #region INonLinearProvider Members

        public double CalculateRhsNorm(IVectorView rhs)
        {
            //TODO: cache the relevant indices.
            //return (new Vector<double>(rhs)).Norm;

            double norm = 0;
            foreach ((INode node, DOFType dofType, int dofIdx) in model.GlobalDofOrdering.GlobalFreeDofs)
            {
                if (dofType != DOFType.Pore) norm += rhs[dofIdx] * rhs[dofIdx];
            }
            return Math.Sqrt(norm);
        }

        public void ProcessInternalRhs(ILinearSystem_v2 linearSystem, IVector rhs, IVectorView solution)
        {
            //ScaleSubdomainSolidVector(subdomain, new Vector<double>(rhs));
            //return;

            rhs.AddIntoThis(qs[linearSystem.Subdomain.ID].Multiply(solution));
            ScaleSubdomainSolidVector(linearSystem, rhs);
        }

        #endregion
    }
}
