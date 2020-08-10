using System.Collections.Generic;
using ISAAR.MSolve.Discretization.Integration.Quadratures;
using ISAAR.MSolve.Geometry.Coordinates;

namespace ISAAR.MSolve.IGA.SupportiveClasses.Interfaces
{
    public interface IShapeFunction1D
    {
        double[,] Values { get; }
        double[,] DerivativeValues { get; }
        double[,] SecondDerivativeValues { get; }

        IReadOnlyList<double[]> EvaluateFunctionsAtGaussPoints(IQuadrature1D quadrature);
        IReadOnlyList<double[,]> EvaluateNaturalDerivativesAtGaussPoints(IQuadrature1D quadrature);
        IReadOnlyList<double[,]> EvaluateNaturalSecondDerivativesAtGaussPoints(IQuadrature1D quadrature);

        double[] EvaluateFunctionsAt(NaturalPoint naturalPoint);
        double[,] EvaluateNaturalDerivativesAt(NaturalPoint naturalPoint);
        double[,] EvaluateNaturalSecondDerivativesAt(NaturalPoint naturalPoint);
    }
}