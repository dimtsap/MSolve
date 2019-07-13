using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Iterative.Preconditioning;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.DomainDecomposition.Overlapping.Schwarz.Additive.CoarseProblem;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Overlapping.Schwarz.Additive
{
    public class OverlappingAdditiveSchwarzPreconditioner : IPreconditioner
    {
        private readonly IMatrix _preconditioner;
        public OverlappingAdditiveSchwarzPreconditioner(ICoarseProblemFactory coarseProblemFactory,
            LocalProblemsFactory localProblemsFactory, IOverlappingDecomposition overlappingDecomposition)
        {
            coarseProblemFactory.GenerateMatrices(overlappingDecomposition);
            localProblemsFactory.GenerateMatrices(overlappingDecomposition);
            var coarseProblemContribution = coarseProblemFactory.RetrievePreconditionerContribution();
            var localProblemsContribution = localProblemsFactory.RetrievePreconditionerContribution();
            _preconditioner = coarseProblemContribution.Add(localProblemsContribution);
            Order = _preconditioner.NumRows;
        }

        /// <summary>
        /// The number of rows/columns of the preconditioner and the original matrix
        /// </summary>
        public int Order { get; }

        public void SolveLinearSystem(IVectorView rhsVector, IVector lhsVector)
        {
            Preconditions.CheckSystemSolutionDimensions(Order, rhsVector.Length);
            lhsVector = _preconditioner.Multiply(rhsVector);
        }
    }

    public class Factory : IPreconditionerFactory
    {
        public readonly IOverlappingDecomposition OverlappingDecomposition;

        public Factory(IOverlappingDecomposition overlappingDecomposition) =>
            OverlappingDecomposition = overlappingDecomposition;

        public ICoarseProblemFactory CoarseProblemFactory { get; set; } = new NestedInterpolationCoarseProblemFactory();

        public LocalProblemsFactory LocalProblemsFactory { get; set; } = new LocalProblemsFactory();

        public IPreconditioner CreatePreconditionerFor(IMatrixView matrix) =>
            new OverlappingAdditiveSchwarzPreconditioner(CoarseProblemFactory, LocalProblemsFactory,
                OverlappingDecomposition);
    }
}