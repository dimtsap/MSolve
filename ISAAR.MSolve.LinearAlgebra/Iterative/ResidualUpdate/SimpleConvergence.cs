﻿using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.LinearAlgebra.Iterative.ResidualUpdate
{
    /// <summary>
    /// Checks if the residual vector has reached the required tolerance, without performing corrections.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class SimpleConvergence : IResidualConvergence
    {
        private double limitDotResidual;

        /// <summary>
        /// See 
        /// <see cref="IResidualConvergence.HasConverged(IVectorView, IVector, ref double, bool, Func{IVectorView, double})"/>.
        /// </summary>
        public bool HasConverged(IVectorView solutionVector, IVector residualVector, ref double residualDotProduct,
            bool isResidualCorrected, Func<IVectorView, double> residualDotCalculation)
            => residualDotProduct <= limitDotResidual;

        /// <summary>
        /// See <see cref="IResidualConvergence.Initialize(IMatrixView, IVectorView, double, double)"/>.
        /// </summary>
        public void Initialize(IMatrixView matrix, IVectorView rhsVector, double residualTolerance, 
            double initialResidualDotProduct)
        {
            // ε^2 * δ0
            // This is more efficient than normalizing and computing square roots. However the order is important, since 
            // ε^2 could be very small and risk precision loss
            this.limitDotResidual = residualTolerance * (residualTolerance * initialResidualDotProduct);
        }
    }
}
