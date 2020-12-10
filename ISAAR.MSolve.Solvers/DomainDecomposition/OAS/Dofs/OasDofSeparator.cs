using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.IGA.Entities;
using ISAAR.MSolve.IGA.Problems.SupportiveClasses;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.DomainDecomposition.OAS.Decomposition;
using ISAAR.MSolve.Solvers.DomainDecomposition.OAS.Mappings;
using ISAAR.MSolve.Solvers.DomainDecomposition.Overlapping.Schwarz.Additive;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.OAS.Dofs
{
    public class OasDofSeparator:IOasDofSeparator
    {
        private readonly IStructuralAsymmetricModel model;
        private readonly int numberOfSubdomainsX;
        private readonly int numberOfSubdomainsY;
        private readonly int overlapping;
        private readonly Dictionary<int, DofTable> subdomainDofOrderingsLocal=new Dictionary<int, DofTable>();
        private readonly Dictionary<int ,int[]> subdomainDofsLocalToFree=new Dictionary<int, int[]>();
        private readonly Dictionary<int, BooleanMatrixRowsToColumns> subdomainToFreeDofMappings=new Dictionary<int, BooleanMatrixRowsToColumns>();
        // private CsrMatrix globalToCoarseMappingMatrix;
        private CscMatrix globalToCoarseMappingMatrix;

        public OasDofSeparator(IStructuralAsymmetricModel model, int numberOfSubdomainsX, int numberOfSubdomainsY, int overlapping=1)
        {
            this.model = model;
            this.numberOfSubdomainsX = numberOfSubdomainsX;
            this.numberOfSubdomainsY = numberOfSubdomainsY;
            this.overlapping = overlapping;
        }

        public int NumberOfSubdomains => numberOfSubdomainsX * numberOfSubdomainsY;

        public BooleanMatrixRowsToColumns GetDofMappingBoundaryClusterToSubdomain(int subdomainID) =>
            subdomainToFreeDofMappings[subdomainID];

        // public CsrMatrix GetGlobalToCoarseMapping => globalToCoarseMappingMatrix;
        public CscMatrix GetGlobalToCoarseMapping => globalToCoarseMappingMatrix;

        public int[] GetDofsSubdomainToFree(int subdomainId) => subdomainDofsLocalToFree[subdomainId];

        public int GetNumFreeDofs(int subdomainID) => subdomainDofOrderingsLocal[subdomainID].EntryCount;

        public void SeparateSubdomainAndCoarseProblemDofs()
        {
            var ksiDecomposition = new AxisDecomposition(model.Subdomains[0].NumberOfCpKsi, numberOfSubdomainsX, overlapping);
            var hetaDecomposition = new AxisDecomposition(model.Subdomains[0].NumberOfCpHeta, numberOfSubdomainsY, overlapping);
            var overlappingDecomposition = new OverlappingDecomposition2D(ksiDecomposition, hetaDecomposition);

            var freeDofs = model.Subdomains[0].FreeDofRowOrdering.FreeDofs;
            
            for (int i = 0; i < overlappingDecomposition.GetNumberOfSubdomains; i++)
            {
                var (subdomainDofOrdering, localToGlobal) = OrderSubdomainDof(overlappingDecomposition, i, freeDofs);

                subdomainDofOrderingsLocal[i] = subdomainDofOrdering;
                subdomainDofsLocalToFree[i] = localToGlobal.ToArray();
            }

            for (int i = 0; i < overlappingDecomposition.GetNumberOfSubdomains; i++)
            {
                var numberOfFreeDofs = freeDofs.EntryCount;
                var numberOfSubdomainDofs = subdomainDofOrderingsLocal[i].EntryCount;
                subdomainToFreeDofMappings[i] = new BooleanMatrixRowsToColumns(numberOfSubdomainDofs, numberOfFreeDofs,
                    subdomainDofsLocalToFree[i]);
            }

            CalculateGlobalToCoarseMapping2(ksiDecomposition, hetaDecomposition);
        }

        private void CalculateGlobalToCoarseMapping2(AxisDecomposition ksiDecomposition, AxisDecomposition hetaDecomposition)
        {
            var cpPositionsKsi =
                CalculateControlPointAxisPositions(model.Subdomains[0].KnotValueVectorKsi, model.Subdomains[0].DegreeKsi);
            var cpPositionsHeta =
                CalculateControlPointAxisPositions(model.Subdomains[0].KnotValueVectorHeta, model.Subdomains[0].DegreeHeta);

            var coarseKnotsKsi = CalculateCoarseAxis(ksiDecomposition, cpPositionsKsi);
            var coarseKnotsHeta = CalculateCoarseAxis(hetaDecomposition, cpPositionsHeta);

            Vector coarseKvvKsi = CalculateKnotValueVectorFromKnots(coarseKnotsKsi, model.Subdomains[0].DegreeKsi);
            Vector coarseKvvHeta = CalculateKnotValueVectorFromKnots(coarseKnotsHeta, model.Subdomains[0].DegreeHeta);

            var cpCoarsePositionsKsi = CalculateControlPointAxisPositions(coarseKvvKsi, model.Subdomains[0].DegreeKsi);
            var cpCoarsePositionsHeta = CalculateControlPointAxisPositions(coarseKvvHeta, model.Subdomains[0].DegreeHeta);

            var coarseCp= new List<IWeightedPoint>();
            var id = 0;
            foreach (var ksi in cpCoarsePositionsKsi)
            {
                foreach (var heta in cpCoarsePositionsHeta)
                {
                    coarseCp.Add(new ControlPoint_Solve()
                    {
                        ID = id++, Ksi = ksi, Heta = heta, WeightFactor = 1.0,
                        X = ksi, Y = heta,
                    });
                }
            }

            var freeModelDof = (model.Subdomains[0].NumberOfCpKsi - 1) * model.Subdomains[0].NumberOfCpHeta * 2;
            var freeCoarseProblemDof = (cpCoarsePositionsKsi.Length - 1) * cpCoarsePositionsHeta.Length * 2;

            var finePoints = new List<NaturalPoint>();
            foreach (var x in cpPositionsKsi)
                finePoints.AddRange(cpPositionsHeta.Select(y => new NaturalPoint(x, y, 0.0)));

            // CSR format
            // var globalToCoarseMatrix = DokRowMajor.CreateEmpty(freeCoarseProblemDof, freeModelDof);
            // for (int i = cpPositionsHeta.Length; i < finePoints.Count; i++)
            // {
            //     var pointShapeFunctions = new NURBS2D_Solve(model.Subdomains[0].DegreeKsi, model.Subdomains[0].DegreeHeta,
            //         coarseKvvKsi, coarseKvvHeta, finePoints[i], coarseCp, true);
            //
            //     var fineIdxFree = i - cpPositionsHeta.Length;
            //     for (int j = cpCoarsePositionsHeta.Length; j < pointShapeFunctions.NurbsValues.NumRows; j++)
            //     {
            //         var coarseIdxFree = j - cpCoarsePositionsHeta.Length;
            //
            //         if (!(Math.Abs(pointShapeFunctions.NurbsValues[j, 0]) > 1e-16)) continue;
            //         globalToCoarseMatrix[2 * coarseIdxFree, 2 * fineIdxFree] = pointShapeFunctions.NurbsValues[j, 0];
            //         globalToCoarseMatrix[2 * coarseIdxFree + 1, 2 * fineIdxFree + 1] = pointShapeFunctions.NurbsValues[j, 0];
            //     }
            // }
            //
            // globalToCoarseMappingMatrix = globalToCoarseMatrix.BuildCsrMatrix(true);
            
            // CSC format
            var globalToCoarseMatrix = DokColMajor.CreateEmpty(freeCoarseProblemDof, freeModelDof);
            for (int i = cpPositionsHeta.Length; i < finePoints.Count; i++)
            {
                var pointShapeFunctions = new NURBS2D_Solve(model.Subdomains[0].DegreeKsi, model.Subdomains[0].DegreeHeta,
                    coarseKvvKsi, coarseKvvHeta, finePoints[i], coarseCp, true);
            
                var fineIdxFree = i - cpPositionsHeta.Length;
                for (int j = cpCoarsePositionsHeta.Length; j < pointShapeFunctions.NurbsValues.NumRows; j++)
                {
                    var coarseIdxFree = j - cpCoarsePositionsHeta.Length;
            
                    if (!(Math.Abs(pointShapeFunctions.NurbsValues[j, 0]) > 1e-16)) continue;
                    globalToCoarseMatrix[2 * coarseIdxFree, 2 * fineIdxFree] = pointShapeFunctions.NurbsValues[j, 0];
                    globalToCoarseMatrix[2 * coarseIdxFree + 1, 2 * fineIdxFree + 1] = pointShapeFunctions.NurbsValues[j, 0];
                }
            }
            
            globalToCoarseMappingMatrix = globalToCoarseMatrix.BuildCscMatrix(true);
        }

        private Vector CalculateKnotValueVectorFromKnots(double[] knots, int degree)
        {
            var knotValueVector = Vector.CreateZero(knots.Length + 2 * degree);

            for (int i = 0; i < degree; i++)
            {
                knotValueVector[i] = knots[0];
            }

            for (int i = 0; i < knots.Length; i++)
            {
                knotValueVector[i + degree] = knots[i];
            }

            for (int i = knots.Length+degree; i < knotValueVector.Length; i++)
            {
                knotValueVector[i] = knots.Last();
            }

            return knotValueVector;
        }

        // private void CalculateGlobalToCoarseMapping1(AxisDecomposition ksiDecomposition, AxisDecomposition hetaDecomposition)
        // {
        //     var cpPositionsKsi =
        //         CalculateControlPointAxisPositions(model.Subdomains[0].KnotValueVectorKsi, model.Subdomains[0].DegreeKsi);
        //     var cpPositionsHeta =
        //         CalculateControlPointAxisPositions(model.Subdomains[0].KnotValueVectorHeta, model.Subdomains[0].DegreeHeta);
        //
        //     var coarseCpKsi = CalculateCoarseAxis(ksiDecomposition, cpPositionsKsi);
        //     var coarseCpHeta = CalculateCoarseAxis(hetaDecomposition, cpPositionsHeta);
        //
        //     var freeModelDof = (model.Subdomains[0].NumberOfCpKsi * model.Subdomains[0].NumberOfCpHeta -
        //                         model.Subdomains[0].NumberOfCpHeta) * 2;
        //     var freeCoarseProblemDof = (coarseCpKsi.Length * coarseCpHeta.Length - coarseCpHeta.Length) * 2;
        //
        //     var coarsePoints = new List<NaturalPoint>();
        //     foreach (var x in coarseCpKsi)
        //         coarsePoints.AddRange(coarseCpHeta.Select(y => new NaturalPoint(x, y, 0.0)));
        //
        //     var modelCps = model.Nodes.Select(x => x as IWeightedPoint).ToList();
        //
        //     var globalToCoarseMatrix = DokRowMajor.CreateEmpty(freeModelDof, freeCoarseProblemDof);
        //     for (int i = coarseCpHeta.Length; i < coarsePoints.Count; i++)
        //     {
        //         var pointShapeFunctions = new NURBS2D_Solve(model.Subdomains[0].DegreeKsi, model.Subdomains[0].DegreeHeta,
        //             model.Subdomains[0].KnotValueVectorKsi,
        //             model.Subdomains[0].KnotValueVectorHeta, coarsePoints[i], modelCps, true);
        //         var coarseIdxFree = i - coarseCpHeta.Length;
        //         for (int j = cpPositionsHeta.Length; j < pointShapeFunctions.NurbsValues.NumRows; j++)
        //         {
        //             var freeIdxGlobal = j - cpPositionsHeta.Length;
        //
        //             if (!(Math.Abs(pointShapeFunctions.NurbsValues[j, 0]) > 1e-16)) continue;
        //             globalToCoarseMatrix[2 * freeIdxGlobal, 2 * coarseIdxFree] =
        //                 pointShapeFunctions.NurbsValues[j, 0];
        //             globalToCoarseMatrix[2 * freeIdxGlobal + 1, 2 * coarseIdxFree + 1] = pointShapeFunctions.NurbsValues[j, 0];
        //         }
        //     }
        //
        //     globalToCoarseMappingMatrix = globalToCoarseMatrix.BuildCsrMatrix(true);
        // }

        private double[] CalculateCoarseAxis(AxisDecomposition axisDecomposition, double[] cpPositions)
        {
            var coarseIndices = new int[axisDecomposition.GetNumberOfSubdomains + 1];
            coarseIndices[0] = 0;
            for (int i = 0; i < axisDecomposition.GetNumberOfSubdomains; i++)
            {
                coarseIndices[i + 1] = axisDecomposition.GetIndicesOfSubdomainCp(i).Last();
            }

            var coarseAxisCoordinates = new double[coarseIndices.Length];
            for (int i = 0; i < coarseIndices.Length; i++)
            {
                coarseAxisCoordinates[i] = cpPositions[coarseIndices[i]];
            }

            return coarseAxisCoordinates;
        }

        private double[] CalculateControlPointAxisPositions(Vector knotValueVector, int degree)
        {
            var numberOfCp = knotValueVector.Length - degree - 1;
            var cpPositions = new double[numberOfCp];
            
            for (int i = 0; i < numberOfCp; i++)
            {
                var sum = 0.0;
                for (int j = 0; j < degree; j++)
                {
                    sum += knotValueVector[i + j+1];
                }

                cpPositions[i] = sum / degree;
            }


            return cpPositions;
        }

        private (DofTable subdomainDofOrdering, List<int> localToGlobal) OrderSubdomainDof(
            OverlappingDecomposition2D overlappingDecomposition, int i, DofTable freeDofs)
        {
            var subdomainDofIdx = 0;
            var subdomainDofOrdering = new DofTable();
            var localToGlobal = new List<int>();
            var cpIndices = overlappingDecomposition.GetControlPointIndicesOfSubdomain(i);
            // var freeControlPoints = freeDofs.GetRows().ToArray();
            // var freeControlPoints = freeDofs.GetRows().ToDictionary(x => x.ID);
            var freeControlPoints = freeDofs.GetRows().ToDictionary(x => x.ID);

            for (int j = 0; j < cpIndices.Length; j++)
            {
                // var node=freeControlPoints.FirstOrDefault(x => x.ID == cpIndices[j]);
                // if (node==null) continue;
                if (!freeControlPoints.ContainsKey(cpIndices[j])) continue;
                var node=freeControlPoints[cpIndices[j]];
                
                var dofsOfNode = freeDofs.GetDataOfRow(node);
                foreach (var dofTypeIdxPair in dofsOfNode)
                {
                    subdomainDofOrdering[node, dofTypeIdxPair.Key] = subdomainDofIdx++;
                    localToGlobal.Add(dofTypeIdxPair.Value);
                }
            }

            return (subdomainDofOrdering, localToGlobal);
        }
    }
}
