﻿using System;
using System.Collections.Generic;
using System.Diagnostics;
using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.Discretization.Integration;
using ISAAR.MSolve.Discretization.Integration.Quadratures;
using ISAAR.MSolve.Discretization.Mesh;
using ISAAR.MSolve.FEM.Interpolation;
using ISAAR.MSolve.FEM.Interpolation.GaussPointExtrapolation;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.Enrichments.Items;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.FreedomDegrees;
using ISAAR.MSolve.XFEM.Integration;
using ISAAR.MSolve.XFEM.Materials;
using ISAAR.MSolve.XFEM.Utilities;

namespace ISAAR.MSolve.XFEM.Elements
{
    /// <summary>
    /// TODO: Uses the same interpolation and nodes as the underlying std element! This must change
    /// TODO: Enumerating artificial dofs may be needed to be done by this class. (e.g if structural FE introduce more 
    ///     artificial dofs than continuum FE)
    /// TODO: The calculation of Kss uses the same gauss points as the calculation of Kes, Kee. 
    ///     Pros: only need to track one set of Gauss points, which simplifies non linear analysis. 
    ///     Cons: calculating Kss with the Gauss points of an enriched element is much more expensive
    /// </summary>
    public class XContinuumElement2D : ICell<XNode>, IFiniteElement2DView<XNode, IIsoparametricInterpolation2D>
    {
        public XContinuumElement2D(IReadOnlyList<XNode> nodes, IIsoparametricInterpolation2D interpolation,
            IGaussPointExtrapolation2D gaussPointExtrapolation, IQuadrature2D standardQuadrature, 
            IIntegrationStrategy2D<XContinuumElement2D> integrationStrategy, 
            IIntegrationStrategy2D<XContinuumElement2D> jIntegralStrategy, IMaterialField2D material)
        {
            this.Nodes = nodes;
            this.Interpolation = interpolation;
            this.GaussPointExtrapolation = gaussPointExtrapolation;
            this.StandardQuadrature = standardQuadrature;
            this.IntegrationStrategy = integrationStrategy;
            this.JintegralStrategy = jIntegralStrategy;
            this.Material = material;
            this.EnrichmentItems = new List<IEnrichmentItem2D>();
            this.NumStandardDofs = 2 * nodes.Count;
        }

        public CellType CellType => Interpolation.CellType;

        /// <summary>
        /// ERROR: elements should not be enriched explicitly. 
        /// Instead the enrichment items should store which elements they interact with. 
        /// If the element needs to access the enrichment items it should do so through its nodes.
        /// Ok, but how would the integration strategy access the enrichment item?
        /// </summary>
        public List<IEnrichmentItem2D> EnrichmentItems { get; }

        public IGaussPointExtrapolation2D GaussPointExtrapolation { get; }
        internal IIntegrationStrategy2D<XContinuumElement2D> IntegrationStrategy { get; }

        /// <summary>
        /// Common interpolation for standard and enriched nodes.
        /// </summary>
        public IIsoparametricInterpolation2D Interpolation { get; }

        internal IIntegrationStrategy2D<XContinuumElement2D> JintegralStrategy { get; }
        internal IMaterialField2D Material { get; }

        /// <summary>
        /// All nodes are enriched for now.
        /// </summary>
        public IReadOnlyList<XNode> Nodes { get; }
        public int NumStandardDofs { get; }
        internal IQuadrature2D StandardQuadrature { get; }


        // TODO: return a symmetric matrix
        public Matrix BuildStandardStiffnessMatrix()
        {
            var stiffness = Matrix.CreateZero(NumStandardDofs, NumStandardDofs);
            IReadOnlyList<EvalInterpolation2D> evaluatedInterpolations = 
                Interpolation.EvaluateAllAtGaussPoints(Nodes, StandardQuadrature);

            for (int gp = 0; gp < StandardQuadrature.IntegrationPoints.Count; ++gp)
            {
                GaussPoint gaussPoint = StandardQuadrature.IntegrationPoints[gp];

                // Material properties
                Matrix constitutive = Material.CalculateConstitutiveMatrixAt(gaussPoint, evaluatedInterpolations[gp]);
                //TODO: The thickness is constant per element in FEM, but what about XFEM? Different materials within the same element are possible.
                double thickness = Material.GetThicknessAt(gaussPoint, evaluatedInterpolations[gp]);

                // Calculate the necessary quantities for the integration
                var jacobian = evaluatedInterpolations[gp].Jacobian;
                Matrix deformation = CalculateStandardDeformationMatrix(evaluatedInterpolations[gp].ShapeGradientsCartesian);

                // Contribution of this gauss point to the element stiffness matrix
                Matrix partial = deformation.ThisTransposeTimesOtherTimesThis(constitutive);
                Debug.Assert(partial.NumRows == NumStandardDofs);
                Debug.Assert(partial.NumColumns == NumStandardDofs);
                double dV = jacobian.DirectDeterminant * gaussPoint.Weight * thickness;
                stiffness.AxpyIntoThis(partial, dV);
            }
            return stiffness;
        }

        public (Matrix stiffnessEnrichedEnriched, Matrix stiffnessEnrichedStandard) BuildEnrichedStiffnessMatricesLower()
        {
            int numStandardDofs = NumStandardDofs;
            int numEnrichedDofs = CountEnrichedDofs();
            var stiffnessEnrichedStandard = Matrix.CreateZero(numEnrichedDofs, numStandardDofs);
            var stiffnessEnrichedEnriched = Matrix.CreateZero(numEnrichedDofs, numEnrichedDofs);

            foreach (GaussPoint gaussPoint in IntegrationStrategy.GenerateIntegrationPoints(this))
            {
                // Calculate the necessary quantities for the integration
                EvalInterpolation2D evaluatedInterpolation = Interpolation.EvaluateAllAt(Nodes, gaussPoint);
                double thickness = Material.GetThicknessAt(gaussPoint, evaluatedInterpolation);
                Matrix constitutive = Material.CalculateConstitutiveMatrixAt(gaussPoint, evaluatedInterpolation);
                Matrix Bstd = CalculateStandardDeformationMatrix(evaluatedInterpolation.ShapeGradientsCartesian);
                Matrix Benr = CalculateEnrichedDeformationMatrix(numEnrichedDofs,
                    gaussPoint, evaluatedInterpolation);

                // Contributions of this gauss point to the element stiffness matrices. 
                // Kee = SUM(Benr^T * E * Benr * dV), Kes = SUM(Benr^T * E * Bstd * dV)
                // TODO: Scale only the smallest matrix (constitutive) before the multiplications. Probably requires a copy of the constitutive matrix.
                double dVolume = thickness * evaluatedInterpolation.Jacobian.DirectDeterminant * gaussPoint.Weight;
                Matrix transposeBenrTimesConstitutive = Benr.MultiplyRight(constitutive, true, false); // cache the result

                Matrix Kes = transposeBenrTimesConstitutive * Bstd;  // enriched-standard part
                stiffnessEnrichedStandard.AxpyIntoThis(Kes, dVolume);

                Matrix Kee = transposeBenrTimesConstitutive * Benr;  // enriched-enriched part
                stiffnessEnrichedEnriched.AxpyIntoThis(Kee, dVolume);
            }

            return (stiffnessEnrichedEnriched, stiffnessEnrichedStandard);
        }

        public (Matrix stiffnessEnrichedEnriched, Matrix stiffnessStandardEnriched) BuildEnrichedStiffnessMatricesUpper()
        {
            int numStandardDofs = NumStandardDofs;
            int numEnrichedDofs = CountEnrichedDofs();
            var stiffnessStandardEnriched = Matrix.CreateZero(numStandardDofs, numEnrichedDofs);
            var stiffnessEnrichedEnriched = Matrix.CreateZero(numEnrichedDofs, numEnrichedDofs);

            foreach (GaussPoint gaussPoint in IntegrationStrategy.GenerateIntegrationPoints(this))
            {
                // Calculate the necessary quantities for the integration
                EvalInterpolation2D evaluatedInterpolation = Interpolation.EvaluateAllAt(Nodes, gaussPoint);
                double thickness = Material.GetThicknessAt(gaussPoint, evaluatedInterpolation);
                Matrix constitutive = Material.CalculateConstitutiveMatrixAt(gaussPoint, evaluatedInterpolation);
                Matrix Bstd = CalculateStandardDeformationMatrix(evaluatedInterpolation.ShapeGradientsCartesian);
                Matrix Benr = CalculateEnrichedDeformationMatrix(numEnrichedDofs,
                    gaussPoint, evaluatedInterpolation);

                // Contributions of this gauss point to the element stiffness matrices. 
                // Kee = SUM(Benr^T * E * Benr * dV), Kse = SUM(Bstd^T * E * Benr * dV)
                // TODO: Scale only the smallest matrix (constitutive) before the multiplications. Probably requires a copy of the constitutive matrix.
                double dVolume = thickness * evaluatedInterpolation.Jacobian.DirectDeterminant * gaussPoint.Weight;
                Matrix constitutiveTimesBenr = constitutive * Benr; // cache the result

                Matrix Kse = Bstd.MultiplyRight(constitutiveTimesBenr, true, false);  // enriched-standard part
                stiffnessStandardEnriched.AxpyIntoThis(Kse, dVolume);

                Matrix Kee = Benr.MultiplyRight(constitutiveTimesBenr, true, false);  // enriched-enriched part
                stiffnessEnrichedEnriched.AxpyIntoThis(Kee, dVolume);
            }

            return (stiffnessEnrichedEnriched, stiffnessStandardEnriched);
        }

        /// <summary>
        /// Calculates the deformation matrix B. Dimensions = 3x8.
        /// B is a linear transformation FROM the nodal values of the displacement field TO the the derivatives of
        /// the displacement field in respect to the cartesian axes (i.e. the stresses): {dU/dX} = [B] * {d} => 
        /// {u,x v,y u,y, v,x} = [... Bk ...] * {u1 v1 u2 v2 u3 v3 u4 v4}, where k = 1, ... nodesCount is a node and
        /// Bk = [dNk/dx 0; 0 dNk/dY; dNk/dy dNk/dx] (3x2)
        /// </summary>
        /// <param name="evaluatedInterpolation">The shape function derivatives calculated at a specific 
        ///     integration point</param>
        /// <returns></returns>
        internal Matrix CalculateStandardDeformationMatrix(Matrix shapeGradientsCartesian)
        {
            var deformation = Matrix.CreateZero(3, NumStandardDofs);
            for (int nodeIdx = 0; nodeIdx < Nodes.Count; ++nodeIdx)
            {
                int col0 = 2 * nodeIdx;
                int col1 = 2 * nodeIdx + 1;

                deformation[0, col0] = shapeGradientsCartesian[nodeIdx, 0];
                deformation[1, col1] = shapeGradientsCartesian[nodeIdx, 1];
                deformation[2, col0] = shapeGradientsCartesian[nodeIdx, 1];
                deformation[2, col1] = shapeGradientsCartesian[nodeIdx, 0];
            }
            return deformation;
        }

        // TODO: the argument asrtificialDofsCount was added when this method was private and only called by 
        // BuildStiffnessMatrix() that already counted the dofs. Since it is now used by other modules 
        // (J-integral, output), it would be better to obscure it, at the cost of recounting the dofs in some cases.
        internal Matrix CalculateEnrichedDeformationMatrix(int artificialDofsCount,
            NaturalPoint gaussPoint, EvalInterpolation2D evaluatedInterpolation)
        {
            //CartesianPoint cartesianPoint = evaluatedInterpolation.TransformPointNaturalToGlobalCartesian(gaussPoint);
            var uniqueEnrichments = new Dictionary<IEnrichmentItem2D, EvaluatedFunction2D[]>();

            var deformationMatrix = Matrix.CreateZero(3, artificialDofsCount);
            int currentColumn = 0;
            for (int nodeIdx = 0; nodeIdx < Nodes.Count; ++nodeIdx)
            {
                double N = evaluatedInterpolation.ShapeFunctions[nodeIdx];
                double dNdx = evaluatedInterpolation.ShapeGradientsCartesian[nodeIdx, 0];
                double dNdy = evaluatedInterpolation.ShapeGradientsCartesian[nodeIdx, 1];

                foreach (var enrichment in Nodes[nodeIdx].EnrichmentItems)
                {
                    IEnrichmentItem2D enrichmentItem = enrichment.Key;
                    double[] nodalEnrichmentValues = enrichment.Value;

                    // The enrichment function probably has been evaluated when processing a previous node. Avoid reevaluation.
                    EvaluatedFunction2D[] evaluatedEnrichments;
                    if (!(uniqueEnrichments.TryGetValue(enrichmentItem, out evaluatedEnrichments)))
                    {
                        evaluatedEnrichments = enrichmentItem.EvaluateAllAt(gaussPoint, this, evaluatedInterpolation);
                        uniqueEnrichments[enrichmentItem] = evaluatedEnrichments;
                    }

                    for (int i = 0; i < evaluatedEnrichments.Length; ++i)
                    {
                        // For each node and with all derivatives w.r.t. cartesian coordinates, the enrichment derivatives 
                        // are: Bx = enrN,x = N,x(x,y) * [H(x,y) - H(node)] + N(x,y) * H,x(x,y), where H is the enrichment 
                        // function
                        double Bx = dNdx * (evaluatedEnrichments[i].Value - nodalEnrichmentValues[i])
                            + N * evaluatedEnrichments[i].CartesianDerivatives[0];
                        double By = dNdy * (evaluatedEnrichments[i].Value - nodalEnrichmentValues[i])
                            + N * evaluatedEnrichments[i].CartesianDerivatives[1];

                        // This depends on the convention: node major or enrichment major. The following is node major.
                        int col1 = currentColumn++;
                        int col2 = currentColumn++;

                        deformationMatrix[0, col1] = Bx;
                        deformationMatrix[1, col2] = By;
                        deformationMatrix[2, col1] = By;
                        deformationMatrix[2, col2] = Bx;
                    }
                }
            }
            Debug.Assert(currentColumn == artificialDofsCount);
            return deformationMatrix;
        }

        /// <summary>
        /// This only works for points that do not lie on the crack interface. As such it is safe to pass GPs only
        /// </summary>
        /// <param name="gaussPoint"></param>
        /// <param name="evaluatedInterpolation"></param>
        /// <param name="standardNodalDisplacements"></param>
        /// <param name="enrichedNodalDisplacements"></param>
        /// <returns></returns>
        internal Vector2 CalculateDisplacementField(NaturalPoint gaussPoint, EvalInterpolation2D evaluatedInterpolation,
            Vector standardNodalDisplacements, Vector enrichedNodalDisplacements)
        {
            #region debug
            double tol = 1e-6;
            if ((Math.Abs(gaussPoint.Xi) <= tol) || (Math.Abs(gaussPoint.Eta) <= tol))
            {
                Console.WriteLine("Found an intersection point that isn't a node.");
            }
            #endregion
            var displacements = new double[2];

            // Standard contributions
            for (int nodeIdx = 0; nodeIdx < Nodes.Count; ++nodeIdx)
            {
                double shapeFunction = evaluatedInterpolation.ShapeFunctions[nodeIdx];
                displacements[0] += shapeFunction * standardNodalDisplacements[2 * nodeIdx];
                displacements[1] += shapeFunction * standardNodalDisplacements[2 * nodeIdx + 1];
            }

            // Enriched contributions
            //TODO: this should be taken as input, so that it is only computed once if it is needed for displacements, strains, 
            //      stresses, just like the evaluated interpolation.
            IReadOnlyDictionary<IEnrichmentItem2D, EvaluatedFunction2D[]> evalEnrichments =
                EvaluateEnrichments(gaussPoint, evaluatedInterpolation);
            int dof = 0;
            for (int nodeIdx = 0; nodeIdx < Nodes.Count; ++nodeIdx)
            {
                double shapeFunction = evaluatedInterpolation.ShapeFunctions[nodeIdx];
                foreach (var nodalEnrichment in Nodes[nodeIdx].EnrichmentItems)
                {
                    EvaluatedFunction2D[] currentEvalEnrichments = evalEnrichments[nodalEnrichment.Key];
                    for (int e = 0; e < currentEvalEnrichments.Length; ++e)
                    {
                        double basisFunction = shapeFunction * (currentEvalEnrichments[e].Value - nodalEnrichment.Value[e]);
                        //double basisFunction = shapeFunction * currentEvalEnrichments[e].Value; // for debugging
                        displacements[0] += basisFunction * enrichedNodalDisplacements[dof++];
                        displacements[1] += basisFunction * enrichedNodalDisplacements[dof++];
                    }
                }
            }

            return Vector2.CreateFromArray(displacements);
        }

        /// <summary>
        /// The displacement field derivatives are a 2x2 matrix: gradientU[i,j] = dui/dj where i is the vector component 
        /// and j is the coordinate, w.r.t which the differentiation is done. The differentation coordinates and the
        /// vector components refer to the global cartesian system. 
        /// </summary>
        /// <param name="evaluatedInterpolation"></param>
        /// <param name="nodalDisplacementsX"></param>
        /// <param name="nodalDisplacementsY"></param>
        /// <returns></returns>
        internal Matrix2by2 CalculateDisplacementFieldGradient(NaturalPoint gaussPoint,
            EvalInterpolation2D evaluatedInterpolation, Vector standardNodalDisplacements,
            Vector enrichedNodalDisplacements) //TODO: this must only allow evaluations at Gauss points. It doesn't work for points on the crack interface
        {
            var displacementGradient = Matrix2by2.CreateZero();

            // Standard contributions
            for (int nodeIdx = 0; nodeIdx < Nodes.Count; ++nodeIdx)
            {
                double displacementX = standardNodalDisplacements[2 * nodeIdx];
                double displacementY = standardNodalDisplacements[2 * nodeIdx + 1];

                double dNdx = evaluatedInterpolation.ShapeGradientsCartesian[nodeIdx, 0];
                double dNdy = evaluatedInterpolation.ShapeGradientsCartesian[nodeIdx, 1];
                displacementGradient[0, 0] += dNdx * displacementX;
                displacementGradient[0, 1] += dNdy * displacementX;
                displacementGradient[1, 0] += dNdx * displacementY;
                displacementGradient[1, 1] += dNdy * displacementY;
            }

            // Enriched contributions. TODO: Extract the common steps with building B into a separate method 
            IReadOnlyDictionary<IEnrichmentItem2D, EvaluatedFunction2D[]> evalEnrichments =
                EvaluateEnrichments(gaussPoint, evaluatedInterpolation);
            int dof = 0;
            for (int nodeIdx = 0; nodeIdx < Nodes.Count; ++nodeIdx)
            {
                double N = evaluatedInterpolation.ShapeFunctions[nodeIdx];
                double dNdx = evaluatedInterpolation.ShapeGradientsCartesian[nodeIdx, 0];
                double dNdy = evaluatedInterpolation.ShapeGradientsCartesian[nodeIdx, 1];

                foreach (var nodalEnrichment in Nodes[nodeIdx].EnrichmentItems)
                {
                    EvaluatedFunction2D[] currentEvalEnrichments = evalEnrichments[nodalEnrichment.Key];
                    for (int e = 0; e < currentEvalEnrichments.Length; ++e)
                    {
                        double psi = currentEvalEnrichments[e].Value;
                        Vector2 gradPsi = currentEvalEnrichments[e].CartesianDerivatives;
                        double deltaPsi = psi - nodalEnrichment.Value[e];

                        double Bx = dNdx * deltaPsi + N * gradPsi[0];
                        double By = dNdy * deltaPsi + N * gradPsi[1];

                        double enrDisplacementX = enrichedNodalDisplacements[dof++];
                        double enrDisplacementY = enrichedNodalDisplacements[dof++];

                        displacementGradient[0, 0] += Bx * enrDisplacementX;
                        displacementGradient[0, 1] += By * enrDisplacementX;
                        displacementGradient[1, 0] += Bx * enrDisplacementY;
                        displacementGradient[1, 1] += By * enrDisplacementY;
                    }
                }
            }

            return displacementGradient;
        }

        // In a non linear problem I would also have to pass the new displacements or I would have to update the
        // material state elsewhere.
        public Tensor2D CalculateStressTensor(Matrix2by2 displacementFieldGradient, Matrix constitutive)
        {
            double strainXX = displacementFieldGradient[0, 0];
            double strainYY = displacementFieldGradient[1, 1];
            double strainXYtimes2 = displacementFieldGradient[0, 1] + displacementFieldGradient[1, 0];

            // Should constitutive also be a tensor? Or  should I use matrices and vectors instead of tensors?
            double stressXX = constitutive[0, 0] * strainXX + constitutive[0, 1] * strainYY;
            double stressYY = constitutive[1, 0] * strainXX + constitutive[1, 1] * strainYY;
            double stressXY = constitutive[2, 2] * strainXYtimes2;

            return new Tensor2D(stressXX, stressYY, stressXY);
        }

        #region Dofs (perhaps all these should be delegated to element specific std and enr DofOrderers)
        public int CountEnrichedDofs()
        {
            int count = 0;
            foreach (XNode node in Nodes) count += node.EnrichedDofsCount; // in all nodes or in enriched interpolation nodes?
            return count;
        }


        internal FreedomDegrees.Ordering.DofTable<EnrichedDof> GetEnrichedDofs()
        {
            var elementDofs = new FreedomDegrees.Ordering.DofTable<EnrichedDof>();
            int dofCounter = 0;
            foreach (XNode node in Nodes)
            {
                foreach (var enrichment in node.EnrichmentItems.Keys)
                {
                    foreach (var enrichedDof in enrichment.Dofs) // there are different dofs for x and y axes
                    {
                        elementDofs[node, enrichedDof] = dofCounter++;
                    }
                }
            }
            return elementDofs;
        }

        /// <summary>
        /// TODO: Perhaps this should be saved as a DofOrderer object (the dofs themselves would be created on  
        /// demand though). XElement will have a mutable one, while others will get a view. I could still use a  
        /// DofOrderer even if I do not save it. Transfering most of the code to the Enumerator class, also reduces  
        /// code duplication with the standard ContinuumElement2D
        /// </summary>
        /// <returns></returns>
        internal FreedomDegrees.Ordering.DofTable<DisplacementDof> GetStandardDofs()
        {
            var elementDofs = new FreedomDegrees.Ordering.DofTable<DisplacementDof>();
            int dofCounter = 0;
            foreach (XNode node in Nodes)
            {
                elementDofs[node, DisplacementDof.X] = dofCounter++;
                elementDofs[node, DisplacementDof.Y] = dofCounter++;
            }
            return elementDofs;
        }
        #endregion

        private IReadOnlyDictionary<IEnrichmentItem2D, EvaluatedFunction2D[]> EvaluateEnrichments(
            NaturalPoint gaussPoint, EvalInterpolation2D evaluatedInterpolation)
        {
            var cachedEvalEnrichments = new Dictionary<IEnrichmentItem2D, EvaluatedFunction2D[]>();
            foreach (XNode node in Nodes)
            {
                foreach (var enrichment in node.EnrichmentItems)
                {
                    IEnrichmentItem2D enrichmentItem = enrichment.Key;
                    double[] nodalEnrichmentValues = enrichment.Value;

                    // The enrichment function probably has been evaluated when processing a previous node. Avoid reevaluation.
                    if (!(cachedEvalEnrichments.TryGetValue(enrichmentItem, out EvaluatedFunction2D[] evaluatedEnrichments)))
                    {
                        evaluatedEnrichments = enrichmentItem.EvaluateAllAt(gaussPoint, this, evaluatedInterpolation);
                        cachedEvalEnrichments[enrichmentItem] = evaluatedEnrichments;
                    }
                }
            }
            return cachedEvalEnrichments;
        }
    }
}