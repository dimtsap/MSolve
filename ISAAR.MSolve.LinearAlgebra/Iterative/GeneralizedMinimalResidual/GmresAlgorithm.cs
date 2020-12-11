﻿using System;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Iterative.Preconditioning;
using ISAAR.MSolve.LinearAlgebra.Iterative.Termination;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.LinearAlgebra.Iterative.GeneralizedMinimalResidual
{
    /// <summary>
    /// Based on Restarted GMRES implementation provided in https://people.sc.fsu.edu/~jburkardt/m_src/mgmres/mgmres.html
    /// </summary>
    public class GmresAlgorithm
    {
        private const string name = "Restarted Generalized minimal residual method";

        private double absoluteTolerance;
        private double relativeTolerance;
        private int maximumIterations;
        private IMaxIterationsProvider innerIterationsProvider;
        protected IVector residual;

        public GmresAlgorithm(double absoluteTolerance, double relativeTolerance, int maximumIterations,
            IMaxIterationsProvider innerIterationsProvider)
        {
            this.absoluteTolerance = absoluteTolerance;
            this.relativeTolerance = relativeTolerance;
            this.maximumIterations = maximumIterations;
            this.innerIterationsProvider = innerIterationsProvider;
        }

        public IterativeStatistics Solve(IMatrixView matrix, IPreconditioner preconditioner,IVectorView rhs, IVector solution,
            bool initialGuessIsZero, Func<IVector> zeroVectorInitializer)
        {
            return Solve(new ExplicitMatrixTransformation(matrix), preconditioner,rhs, solution, initialGuessIsZero,
                zeroVectorInitializer);
        }


        public IterativeStatistics Solve(ILinearTransformation matrix, IPreconditioner preconditioner, IVectorView rhs, IVector solution,
            bool initialGuessIsZero, Func<IVector> zeroVectorInitializer)
        {
            Preconditions.CheckMultiplicationDimensions(matrix.NumColumns, solution.Length);
            Preconditions.CheckSystemSolutionDimensions(matrix.NumRows, rhs.Length);

            var innerIterations = innerIterationsProvider.GetMaxIterations(matrix.NumRows);
            IVector[] v =new Vector[innerIterations+1];
            var y = Vector.CreateZero(innerIterations + 1);
            var c= Vector.CreateZero(innerIterations + 1);
            var s = Vector.CreateZero(innerIterations + 1);
            var delta = 0.001;
            double residualNorm = double.MaxValue;
            var usedIterations = 0;

            if (initialGuessIsZero) residual = rhs.Copy();
            else residual = ExactResidual.Calculate(matrix, rhs, solution);

            for ( var iteration = 0; iteration < maximumIterations; iteration++)
            {
                var preconditionedResidual = Vector.CreateZero(residual.Length);
                preconditioner.SolveLinearSystem(residual, preconditionedResidual);
                residual = preconditionedResidual;
                //preconditioner.SolveLinearSystem(residual, residual);
                //var residual = ExactResidual.Calculate(matrix, rhs, solution);

                residualNorm = residual.Norm2();

                double residualTolerance;
                if (iteration == 0)
                    residualTolerance = residualNorm * relativeTolerance;

                v[0] = residual.Scale(1 / residualNorm);

                var g = Vector.CreateZero(innerIterations+1);
                g[0] = residualNorm;
                var hessenbergMatrix = DokColMajor.CreateEmpty(innerIterations + 1, innerIterations);
                //var hessenbergMatrix = Matrix.CreateZero(innerIterations + 1, innerIterations);

                var indexIteration = 0;
                for (int innerIteration = 0; innerIteration < innerIterations; innerIteration++)
                {
                    indexIteration = innerIteration;
                    v[innerIteration + 1] = Vector.CreateZero(v[innerIteration].Length);

                    matrix.Multiply(v[innerIteration], v[innerIteration + 1]);

                    var preconditionedV = Vector.CreateZero(residual.Length);
                    preconditioner.SolveLinearSystem(v[innerIteration + 1], preconditionedV);
                    v[innerIteration + 1] = preconditionedV;

                    //preconditioner.SolveLinearSystem(v[innerIteration + 1], v[innerIteration + 1]);

                    var av = v[innerIteration + 1].Norm2();

                    for (var j = 0; j <= innerIteration; j++)
                    {
                        hessenbergMatrix[j, innerIteration] = v[j].DotProduct(v[innerIteration + 1]);
                        v[innerIteration + 1] = v[innerIteration + 1].Subtract(v[j].Scale(hessenbergMatrix[j, innerIteration]));
                    }

                    hessenbergMatrix[innerIteration + 1, innerIteration] = v[innerIteration + 1].Norm2();


                    if (Math.Abs(av + delta * hessenbergMatrix[innerIteration + 1, innerIteration] - av) < 10e-9)
                    {
                        for (int j = 0; j <= innerIteration; j++)
                        {
                            var htmp = v[j].DotProduct(v[innerIteration + 1]);
                            hessenbergMatrix[j, innerIteration] += htmp;
                            v[innerIteration+1].LinearCombinationIntoThis(1.0,v[j],-htmp);
                        }

                        hessenbergMatrix[innerIteration + 1, innerIteration] = v[innerIteration + 1].Norm2();
                    }

                    if (Math.Abs(hessenbergMatrix[innerIteration+1,innerIteration]) > 10e-17)
                        v[innerIteration+1].ScaleIntoThis(1/hessenbergMatrix[innerIteration+1,innerIteration]);

                    if (innerIteration > 0)
                    {

                        //y = hessenbergMatrix.GetColumn(innerIteration).GetSubvector(0,innerIteration+2);
                        y = Vector.CreateZero(innerIteration + 2);
                        for (int i = 0; i < innerIteration+2; i++)
                        {
                            y[i] = hessenbergMatrix[i, innerIteration];
                        }


                        for (int i = 0; i <= innerIteration-1; i++)
                        {
                            y = CalculateGivensRotation(c[i], s[i], i, y);
                        }

                        //hessenbergMatrix.SetSubcolumn(innerIteration, y);
                        for (int i = 0; i < y.Length; i++)
                        {
                            hessenbergMatrix[i, innerIteration] = y[i];
                        }
                    }

                    var mu = Math.Sqrt(hessenbergMatrix[innerIteration, innerIteration] *
                                       hessenbergMatrix[innerIteration, innerIteration] +
                                       hessenbergMatrix[innerIteration + 1, innerIteration] *
                                       hessenbergMatrix[innerIteration + 1, innerIteration]);
                    c[innerIteration] = hessenbergMatrix[innerIteration, innerIteration] / mu;
                    s[innerIteration] = -hessenbergMatrix[innerIteration + 1, innerIteration] / mu;

                    hessenbergMatrix[innerIteration, innerIteration] =
                        c[innerIteration] * hessenbergMatrix[innerIteration, innerIteration] -
                        s[innerIteration] * hessenbergMatrix[innerIteration + 1, innerIteration];

                    hessenbergMatrix[innerIteration + 1, innerIteration] = 0.0;
                    g = CalculateGivensRotation(c[innerIteration], s[innerIteration], innerIteration, g);

                    residualNorm = Math.Abs(g[innerIteration + 1]);
                    usedIterations++;

                    if (residualNorm <= relativeTolerance && residualNorm <= absoluteTolerance)
                        break;
                }

                indexIteration = indexIteration - 1;

                y[indexIteration + 1] =g[indexIteration + 1] / hessenbergMatrix[indexIteration + 1, indexIteration + 1];
                for (int i = indexIteration; i >= 0; i--)
                {
                    var aux = Vector.CreateZero(indexIteration + 2 - (i + 1));
                    for (int j = i+1; j < indexIteration+2; j++)
                    {
                        aux[j-(i+1)] = hessenbergMatrix[i, j];
                    }

                    //y[i] = (g[i] - (hessenbergMatrix.GetRow(i).GetSubvector(i+1, indexIteration + 2)
                    //                    .DotProduct(y.GetSubvector(i + 1, indexIteration + 2)) ))/ hessenbergMatrix[i, i];
                    y[i] = (g[i] - (aux.DotProduct(y.GetSubvector(i + 1, indexIteration + 2)) ))/ hessenbergMatrix[i, i];
                }


                for (int i = 0; i < matrix.NumRows; i++)
                {
                    var subV = Vector.CreateZero(indexIteration + 2);
                    for (int j = 0; j < indexIteration + 2; j++)
                    {
                        subV[j] = v[j][i];
                    }
                    
                    solution.Set(i,solution[i]+ subV.DotProduct(y.GetSubvector(0,indexIteration+2)));
                }

                if (residualNorm <= relativeTolerance && residualNorm <= absoluteTolerance)
                    break;
            }
            return new IterativeStatistics()
            {
                HasConverged = residualNorm <= relativeTolerance && residualNorm <= absoluteTolerance,
                AlgorithmName = name,
                NumIterationsRequired = usedIterations,
                ResidualNormRatioEstimation = residualNorm
            };
        }

        private Vector CalculateGivensRotation(double c, double s, int k, Vector g)
        {
            var g1 = c * g[k] - s * g[k + 1];
            var g2 = s * g[k] + c * g[k + 1];

            g.Set(k, g1);
            g.Set(k+1, g2);

            return g;
        }

        public class Builder
        {
            public int MaximumIterations { get; set; } = 100000;

            public IMaxIterationsProvider InnerIterationsProvider { get; set; } =
                new PercentageMaxIterationsProvider(1.0);

            /// <summary>
            /// Absolute tolerance of the current residual.
            /// </summary>
            public double AbsoluteTolerance { get; set; } = 10e-6;

            /// <summary>
            /// Relative tolerance of the current residual compared to the initial residual.
            /// </summary>
            public double RelativeTolerance { get; set; } = 10e-9;

            public GmresAlgorithm Build() => new GmresAlgorithm(AbsoluteTolerance, RelativeTolerance, MaximumIterations,
                InnerIterationsProvider);
        }
    }
}