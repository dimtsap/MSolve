using System.Collections.Generic;
using ISAAR.MSolve.Discretization.Integration.Quadratures;
using ISAAR.MSolve.Geometry.Coordinates;

namespace ISAAR.MSolve.IGA.SupportiveClasses.Interfaces
{
    public interface IShapeFunction3D
    {
        double[,] Values { get; }

        double[,] DerivativeValuesKsi { get; }
        double[,] DerivativeValuesHeta { get; }
        double[,] DerivativeValuesZeta { get; }

        double[,] SecondDerivativeValuesKsi { get; }
        double[,] SecondDerivativeValuesHeta { get; }
        double[,] SecondDerivativeValuesZeta { get; }

        double[,] SecondDerivativeValuesKsiHeta { get; }
        double[,] SecondDerivativeValuesHetaZeta { get; }
        double[,] SecondDerivativeValuesKsiZeta { get; }

        double[] EvaluateFunctionsAt(NaturalPoint naturalPoint);
        double[,] EvaluateNaturalDerivativesAt(NaturalPoint naturalPoint);
        double[,] EvaluateNaturalSecondDerivativesAt(NaturalPoint naturalPoint);

        IReadOnlyList<double[]> EvaluateFunctionsAtGaussPoints(IQuadrature3D quadrature);
        IReadOnlyList<double[,]> EvaluateNaturalDerivativesAtGaussPoints(IQuadrature3D quadrature);
        IReadOnlyList<double[,]> EvaluateNaturalSecondDerivativesAtGaussPoints(IQuadrature3D quadrature);
        
    }
}