using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.DomainDecomposition.Overlapping.Schwarz.Additive.Interfaces;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Overlapping.Schwarz.Additive
{
    public class ModelOverlappingDecomposition2D : IModelOverlappingDecomposition
    {
        private int[][] _connectivity;
        private IAxisOverlappingDecomposition _ksiDecomposition;
        private IAxisOverlappingDecomposition _hetaDecomposition;
        private int _numberOfCpHeta;
        private List<IWeightedPoint> _patchControlPoints;
        private CsrMatrix _coarseSpaceInterpolation;

        public ModelOverlappingDecomposition2D(IAxisOverlappingDecomposition ksiDecomposition,
            IAxisOverlappingDecomposition hetaDecomposition, int numberOfCpHeta, List<IWeightedPoint> patchControlPoints)
        {
            _ksiDecomposition = ksiDecomposition;
            _hetaDecomposition = hetaDecomposition;
            _numberOfCpHeta = numberOfCpHeta;
            _patchControlPoints = patchControlPoints;
        }

        public int NumberOfSubdomains => _connectivity.GetLength(0);

        public CsrMatrix CoarseSpaceInterpolation => _coarseSpaceInterpolation;

        public void DecomposeMatrix()
        {
            _ksiDecomposition.DecomposeAxis();
            _hetaDecomposition.DecomposeAxis();

            CalculateLocalProblems();
            CreateCoarseSpaceInterpolation();
        }

        private void CreateCoarseSpaceInterpolation()
        {
            var degreKsi = _ksiDecomposition.Degree;
            var knotValueVectorKsi = _ksiDecomposition.KnotValueVector;
            var coarsePointsKsi = _ksiDecomposition.GetAxisCoarsePoints();

            var degreeHeta = _hetaDecomposition.Degree;
            var knotValueVectorHeta = _hetaDecomposition.KnotValueVector;
            var coarsePointsHeta = _hetaDecomposition.GetAxisCoarsePoints();
            var coarsePoints2D =
                (from ksi in coarsePointsKsi from heta in coarsePointsHeta select new NaturalPoint(ksi, heta)).ToList();

            var R0 = DokRowMajor.CreateEmpty(_patchControlPoints.Count * 2, coarsePoints2D.Count * 2);
            for (int i = 0; i < coarsePoints2D.Count; i++)
            {
                var NurbsValues = CalculateNURBS2D(degreKsi, degreeHeta, knotValueVectorKsi,
                    knotValueVectorHeta, coarsePoints2D[i], _patchControlPoints);

                for (int j = 0; j < NurbsValues.NumRows; j++)
                {
                    if (Math.Abs(NurbsValues[j, 0]) < 10e-9) continue;
                    R0[2 * j, 2 * i] = NurbsValues[j, 0];
                    R0[2 * j + 1, 2 * i + 1] = NurbsValues[j, 0];
                }
            }

            _coarseSpaceInterpolation = R0.BuildCsrMatrix(true);
        }

        private Matrix CalculateNURBS2D(int degreeKsi, int degreeHeta, IVector knotValueVectorKsi,
            IVector knotValueVectorHeta, NaturalPoint naturalPoint, List<IWeightedPoint> patchControlPoints)
        {
            var numberOfCPKsi = knotValueVectorKsi.Length - degreeKsi - 1;
            var numberOfCPHeta = knotValueVectorHeta.Length - degreeHeta - 1;

            var spanKsi = FindSpan(numberOfCPKsi, degreeKsi, naturalPoint.Xi, knotValueVectorKsi);
            var spanHeta = FindSpan(numberOfCPHeta, degreeHeta, naturalPoint.Eta, knotValueVectorHeta);

            var pointFunctionsKsi = BasisFunctions(spanKsi, naturalPoint.Xi, degreeKsi, knotValueVectorKsi);
            var pointFunctionsHeta = BasisFunctions(spanHeta, naturalPoint.Eta, degreeHeta, knotValueVectorHeta);

            int numberOfElementControlPoints = patchControlPoints.Count;

            var nurbsValues = Matrix.CreateZero(numberOfElementControlPoints, 1);

            double sumKsiHeta = 0;

            for (int k = 0; k < numberOfElementControlPoints; k++)
            {
                int indexKsi = patchControlPoints[k].ID / numberOfCPHeta;
                int indexHeta = patchControlPoints[k].ID % numberOfCPHeta;
                sumKsiHeta += pointFunctionsKsi[indexKsi] *
                              pointFunctionsHeta[indexHeta] *
                              patchControlPoints[k].WeightFactor;
            }

            for (int k = 0; k < numberOfElementControlPoints; k++)
            {
                int indexKsi = patchControlPoints[k].ID / numberOfCPHeta;
                int indexHeta = patchControlPoints[k].ID % numberOfCPHeta;

                nurbsValues[k, 0] =
                    pointFunctionsKsi[indexKsi] *
                    pointFunctionsHeta[indexHeta] *
                    patchControlPoints[k].WeightFactor / sumKsiHeta;
            }

            return nurbsValues;
        }

        public static int FindSpan(int numberOfBasisFunctions, int degree, double pointCoordinate, IVector knotValueVector)
        {
            if (pointCoordinate == knotValueVector[numberOfBasisFunctions + 1]) return numberOfBasisFunctions;
            int minimum = degree;
            int maximum = numberOfBasisFunctions + 1;
            int mid = (minimum + maximum) / 2;
            while (pointCoordinate < knotValueVector[mid] || pointCoordinate >= knotValueVector[mid + 1])
            {
                if (pointCoordinate < knotValueVector[mid])
                    maximum = mid;
                else
                    minimum = mid;
                mid = (minimum + maximum) / 2;
            }

            return mid;
        }

        public static Vector BasisFunctions(int spanId, double pointCoordinate, int degree, IVector knotValueVector)
        {
            var basisFunctions = Vector.CreateZero(degree + 1);
            var left = Vector.CreateZero(degree + 1);
            var right = Vector.CreateZero(degree + 1);
            basisFunctions[0] = 1;
            for (int j = 1; j <= degree; j++)
            {
                left[j] = pointCoordinate - knotValueVector[spanId + 1 - j];
                right[j] = knotValueVector[spanId + j] - pointCoordinate;
                var saved = 0.0;
                for (int r = 0; r < j; r++)
                {
                    var temp = basisFunctions[r] / (right[r + 1] + left[j - r]);
                    basisFunctions[r] = saved + right[r + 1] * temp;
                    saved = left[j - r] * temp;
                }

                basisFunctions[j] = saved;
            }

            return basisFunctions;
        }

        private void CalculateLocalProblems()
        {
            var numberOfSubdomainKsi = _ksiDecomposition.NumberOfAxisSubdomains;
            var numberOfSubdomainHeta = _hetaDecomposition.NumberOfAxisSubdomains;
            _connectivity = new int[numberOfSubdomainKsi * numberOfSubdomainHeta][];

            var indexSubdomain = -1;
            for (int i = 0; i < numberOfSubdomainKsi; i++)
            {
                var subdomainKsiConnectivity = _ksiDecomposition.GetAxisIndicesOfSubdomain(i);
                for (int j = 0; j < numberOfSubdomainHeta; j++)
                {
                    var subdomainHetaConnectivity = _hetaDecomposition.GetAxisIndicesOfSubdomain(j);
                    indexSubdomain++;
                    _connectivity[indexSubdomain] = new int[subdomainKsiConnectivity.Length * subdomainHetaConnectivity.Length];

                    var indexCP = 0;
                    foreach (var indexKsi in subdomainKsiConnectivity)
                    {
                        foreach (var indexHeta in subdomainHetaConnectivity)
                        {
                            var globalIndex = indexKsi * _numberOfCpHeta + indexHeta;
                            _connectivity[indexSubdomain][2 * indexCP] = 2 * globalIndex;
                            _connectivity[indexSubdomain][2 * indexCP + 1] = 2 * globalIndex + 1;
                            indexCP++;
                        }
                    }
                }
            }
        }

        public int[] GetConnectivityOfSubdomain(int indexSubdomain) => _connectivity[indexSubdomain];

    }
}