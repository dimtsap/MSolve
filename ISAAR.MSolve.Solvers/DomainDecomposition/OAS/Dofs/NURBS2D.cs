using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.LinearAlgebra.Interpolation;
using TriangleNet.Topology.DCEL;

namespace ISAAR.MSolve.IGA.Problems.SupportiveClasses
{
    public class NURBS2D_Solve
    {
		public Matrix NurbsValues { get; private set; }
        public Matrix NurbsDerivativeValuesKsi { get; private set; }
        public Matrix NurbsDerivativeValuesHeta { get; private set; }
	    public Matrix NurbsSecondDerivativeValueKsi { get; private set; }
		public Matrix NurbsSecondDerivativeValueHeta { get; private set; }
		public Matrix NurbsSecondDerivativeValueKsiHeta { get; private set; }


        public NURBS2D_Solve(int degreeKsi, int degreeHeta, Vector knotValueVectorKsi, Vector knotValueVectorHeta, NaturalPoint collocationPoint, IList<IWeightedPoint> controlPoints, bool calculateAllFunctions)
        {
            var numberOfControlPointsHeta = knotValueVectorHeta.Length - degreeHeta - 1;

            BSPLines1D bsplinesKsi = new BSPLines1D(degreeKsi, knotValueVectorKsi, Vector.CreateFromArray(new double[] { collocationPoint.Xi }));
            BSPLines1D bsplinesHeta = new BSPLines1D(degreeHeta, knotValueVectorHeta, Vector.CreateFromArray(new double[] { collocationPoint.Eta }));
            bsplinesKsi.calculateBSPLinesAndDerivatives();
            bsplinesHeta.calculateBSPLinesAndDerivatives();


            int numberOfGPKsi = 1;
            int numberOfGPHeta = 1;
            int numberOfElementControlPoints = controlPoints.Count;

            NurbsValues = Matrix.CreateZero(numberOfElementControlPoints, 1);
            NurbsDerivativeValuesKsi = Matrix.CreateZero(numberOfElementControlPoints, 1);
            NurbsDerivativeValuesHeta = Matrix.CreateZero(numberOfElementControlPoints, 1);
            NurbsSecondDerivativeValueKsi = Matrix.CreateZero(numberOfElementControlPoints, 1);
            NurbsSecondDerivativeValueHeta = Matrix.CreateZero(numberOfElementControlPoints, 1);
            NurbsSecondDerivativeValueKsiHeta = Matrix.CreateZero(numberOfElementControlPoints, 1);

            for (int i = 0; i < numberOfGPKsi; i++)
            {
                for (int j = 0; j < numberOfGPHeta; j++)
                {
                    double sumKsiHeta = 0;
                    double sumdKsiHeta = 0;
                    double sumKsidHeta = 0;
                    double sumdKsidKsi = 0;
                    double sumdHetadHeta = 0;
                    double sumdKsidHeta = 0;

                    for (int k = 0; k < numberOfElementControlPoints; k++)
                    {
                        // Why type casting is needed.?

                        int indexKsi = controlPoints[k].ID / numberOfControlPointsHeta;
                        int indexHeta = controlPoints[k].ID % numberOfControlPointsHeta;
                        sumKsiHeta += bsplinesKsi.BSPLineValues[indexKsi, i] *
                                      bsplinesHeta.BSPLineValues[indexHeta, j] *
                                      controlPoints[k].WeightFactor;
                        sumdKsiHeta += bsplinesKsi.BSPLineDerivativeValues[indexKsi, i] *
                                       bsplinesHeta.BSPLineValues[indexHeta, j] *
                                       controlPoints[k].WeightFactor;
                        sumKsidHeta += bsplinesKsi.BSPLineValues[indexKsi, i] *
                                       bsplinesHeta.BSPLineDerivativeValues[indexHeta, j] *
                                       controlPoints[k].WeightFactor;
                        sumdKsidKsi += bsplinesKsi.BSPLineSecondDerivativeValues[indexKsi, i] *
                                       bsplinesHeta.BSPLineValues[indexHeta, j] *
                                       controlPoints[k].WeightFactor;
                        sumdHetadHeta += bsplinesKsi.BSPLineValues[indexKsi, i] *
                                         bsplinesHeta.BSPLineSecondDerivativeValues[indexHeta, j] *
                                         controlPoints[k].WeightFactor;
                        sumdKsidHeta += bsplinesKsi.BSPLineDerivativeValues[indexKsi, i] *
                                        bsplinesHeta.BSPLineDerivativeValues[indexHeta, j] *
                                        controlPoints[k].WeightFactor;
                    }
                    for (int k = 0; k < numberOfElementControlPoints; k++)
                    {
                        int indexKsi = controlPoints[k].ID / numberOfControlPointsHeta;
                        int indexHeta = controlPoints[k].ID % numberOfControlPointsHeta;

                        NurbsValues[k, i * numberOfGPHeta + j] =
                            bsplinesKsi.BSPLineValues[indexKsi, i] *
                            bsplinesHeta.BSPLineValues[indexHeta, j] *
                            controlPoints[k].WeightFactor / sumKsiHeta;

                        NurbsDerivativeValuesKsi[k, i * numberOfGPHeta + j] =
                            bsplinesHeta.BSPLineValues[indexHeta, j] * controlPoints[k].WeightFactor *
                            (bsplinesKsi.BSPLineDerivativeValues[indexKsi, i] * sumKsiHeta -
                            bsplinesKsi.BSPLineValues[indexKsi, i] * sumdKsiHeta) / Math.Pow(sumKsiHeta, 2);

                        NurbsDerivativeValuesHeta[k, i * numberOfGPHeta + j] =
                            bsplinesKsi.BSPLineValues[indexKsi, i] * controlPoints[k].WeightFactor *
                            (bsplinesHeta.BSPLineDerivativeValues[indexHeta, j] * sumKsiHeta -
                            bsplinesHeta.BSPLineValues[indexHeta, j] * sumKsidHeta) / Math.Pow(sumKsiHeta, 2);

                        NurbsSecondDerivativeValueKsi[k, i * numberOfGPHeta + j] =
                            bsplinesHeta.BSPLineValues[indexHeta, j] * controlPoints[k].WeightFactor *
                            (bsplinesKsi.BSPLineSecondDerivativeValues[indexKsi, i] / sumKsiHeta -
                             2 * bsplinesKsi.BSPLineDerivativeValues[indexKsi, i] * sumdKsiHeta / Math.Pow(sumKsiHeta, 2) -
                             bsplinesKsi.BSPLineValues[indexKsi, i] * sumdKsidKsi / Math.Pow(sumKsiHeta, 2) +
                             2 * bsplinesKsi.BSPLineValues[indexKsi, i] * Math.Pow(sumdKsiHeta, 2) / Math.Pow(sumKsiHeta, 3));

                        NurbsSecondDerivativeValueHeta[k, i * numberOfGPHeta + j] =
                            bsplinesKsi.BSPLineValues[indexKsi, i] * controlPoints[k].WeightFactor *
                            (bsplinesHeta.BSPLineSecondDerivativeValues[indexHeta, j] / sumKsiHeta -
                             2 * bsplinesHeta.BSPLineDerivativeValues[indexHeta, j] * sumKsidHeta / Math.Pow(sumKsiHeta, 2) -
                             bsplinesHeta.BSPLineValues[indexHeta, j] * sumdHetadHeta / Math.Pow(sumKsiHeta, 2) +
                             2 * bsplinesHeta.BSPLineValues[indexHeta, j] * Math.Pow(sumKsidHeta, 2) / Math.Pow(sumKsiHeta, 3));

                        NurbsSecondDerivativeValueKsiHeta[k, i * numberOfGPHeta + j] =
                            controlPoints[k].WeightFactor *
                            (bsplinesKsi.BSPLineDerivativeValues[indexKsi, i] *
                             bsplinesHeta.BSPLineDerivativeValues[indexHeta, j] / sumKsiHeta -
                             bsplinesKsi.BSPLineDerivativeValues[indexKsi, i] * bsplinesHeta.BSPLineValues[indexHeta, j] *
                             sumKsidHeta / Math.Pow(sumKsiHeta, 2) -
                             bsplinesKsi.BSPLineValues[indexKsi, i] * bsplinesHeta.BSPLineDerivativeValues[indexHeta, j] *
                             sumdKsiHeta / Math.Pow(sumKsiHeta, 2) -
                             bsplinesKsi.BSPLineValues[indexKsi, i] * bsplinesHeta.BSPLineValues[indexHeta, j] *
                             sumdKsidHeta / Math.Pow(sumKsiHeta, 2) +
                             2 * bsplinesKsi.BSPLineValues[indexKsi, i] * bsplinesHeta.BSPLineValues[indexHeta, j] *
                             sumdKsiHeta * sumKsidHeta / Math.Pow(sumKsiHeta, 3));
                    }
                }
            }
        }


	}
}
