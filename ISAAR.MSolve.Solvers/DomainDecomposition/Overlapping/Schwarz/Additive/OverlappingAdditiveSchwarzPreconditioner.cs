using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Iterative.Preconditioning;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.DomainDecomposition.Overlapping.Schwarz.Additive.CoarseProblem;
using ISAAR.MSolve.Solvers.DomainDecomposition.Overlapping.Schwarz.Additive.Interfaces;
using ISAAR.MSolve.Solvers.Ordering;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Overlapping.Schwarz.Additive
{
    public class OverlappingAdditiveSchwarzPreconditioner : IPreconditioner
    {
        private readonly IMatrix _preconditioner;

        public OverlappingAdditiveSchwarzPreconditioner(IMatrixView matrix, ICoarseProblemFactory coarseProblemFactory,
            LocalProblemsFactory localProblemsFactory, IModelOverlappingDecomposition modelOverlappingDecomposition,
            IStructuralAsymmetricModel model)
        {
            modelOverlappingDecomposition.DecomposeMatrix();
            coarseProblemFactory.GenerateMatrices(matrix, modelOverlappingDecomposition);
            localProblemsFactory.GenerateMatrices(matrix, modelOverlappingDecomposition);
            var coarseProblemContribution = coarseProblemFactory.RetrievePreconditionerContribution();
            var localProblemsContribution = localProblemsFactory.RetrievePreconditionerContribution();

            var freeSubdomainDofs= model.GlobalRowDofOrdering.MapFreeDofsSubdomainToGlobal(model.Subdomains[0]);
            _preconditioner = coarseProblemContribution.Add(localProblemsContribution).GetSubmatrix(freeSubdomainDofs, freeSubdomainDofs);
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
        private readonly IStructuralAsymmetricModel _model;
        public readonly IModelOverlappingDecomposition ModelOverlappingDecomposition;
        private readonly List<IWeightedPoint> _patchControlPoints;

        public Factory(IModelOverlappingDecomposition modelOverlappingDecomposition, IStructuralAsymmetricModel model)
        {
            _model = model;
            ModelOverlappingDecomposition = modelOverlappingDecomposition;
        }

        public ICoarseProblemFactory CoarseProblemFactory { get; set; } = new NestedInterpolationCoarseProblemFactory();

        public LocalProblemsFactory LocalProblemsFactory { get; set; } = new LocalProblemsFactory();

        public IPreconditioner CreatePreconditionerFor(IMatrixView matrix) =>
            new OverlappingAdditiveSchwarzPreconditioner(matrix, CoarseProblemFactory, LocalProblemsFactory,
                ModelOverlappingDecomposition, _model);
    }
}