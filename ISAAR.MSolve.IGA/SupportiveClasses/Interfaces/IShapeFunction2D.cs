using System.Collections.Generic;
using ISAAR.MSolve.Discretization.Integration.Quadratures;
using ISAAR.MSolve.Geometry.Coordinates;

namespace ISAAR.MSolve.IGA.SupportiveClasses.Interfaces
{
    public interface IShapeFunction2D
    {
        double[,] Values { get; }
        double[,] DerivativeValuesHeta { get;  }
        double[,] DerivativeValuesKsi { get; }
        double[,] SecondDerivativeValuesHeta { get; }
        double[,] SecondDerivativeValuesKsi { get; }
        double[,] SecondDerivativeValuesKsiHeta { get; }
        
        double[] EvaluateFunctionsAt(NaturalPoint naturalPoint);
        double[,] EvaluateNaturalDerivativesAt(NaturalPoint naturalPoint);
        double[,] EvaluateNaturalSecondDerivativesAt(NaturalPoint naturalPoint);

        IReadOnlyList<double[]> EvaluateFunctionsAtGaussPoints(IQuadrature2D quadrature);
        IReadOnlyList<double[,]> EvaluateNaturalDerivativesAtGaussPoints(IQuadrature2D quadrature);
        IReadOnlyList<double[,]> EvaluateNaturalSecondDerivativesAtGaussPoints(IQuadrature2D quadrature);
    }
}