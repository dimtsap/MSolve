﻿using ISAAR.MSolve.Geometry.Coordinates;

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

        double[,] CalculateShapeFunctionsAt(NaturalPoint point);

        double[,] CalculateShapeFunctionsAt(NaturalPoint[] points);
    }
}