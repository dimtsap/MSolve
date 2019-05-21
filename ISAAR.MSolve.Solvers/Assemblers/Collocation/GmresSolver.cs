﻿using System;
using System.Collections.Generic;
using System.Net.Http.Headers;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Iterative;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.Solvers.Commons;
using ISAAR.MSolve.Solvers.Ordering;
using ISAAR.MSolve.Solvers.Ordering.Reordering;

namespace ISAAR.MSolve.Solvers.Assemblers.Collocation
{
    public class GmresSolver : SingleSubdomainRectangularSolverBase<CsrMatrix>
    {
        private readonly GmresAlgorithm gmresAlgorithm;

        public GmresSolver(IStructuralAsymmetricModel model,  AsymmetricDofOrderer dofRowOrderer, IDofOrderer dofColOrderer)
            : base(model, dofRowOrderer, dofColOrderer, new CsrRectangularAssembler(true), "GmresSolver")
        {
            this.gmresAlgorithm= new GmresAlgorithm.Builder().Build();
        }

        public override void Initialize() { }

        public override void HandleMatrixWillBeSet()
        {
        }

        public override void PreventFromOverwrittingSystemMatrices()
        {
        }

        public override void Solve()
        {
            if (linearSystem == null) linearSystem.Solution = linearSystem.CreateZeroVector();
            else linearSystem.Solution.Clear();

            IterativeStatistics stats = gmresAlgorithm.Solve(linearSystem.Matrix, linearSystem.RhsVector,
                linearSystem.Solution, true, () => linearSystem.CreateZeroVector());

            if (!stats.HasConverged)
                throw new IterativeSolverNotConvergedException("Gmres did not converge");
        }

        protected override Matrix InverseSystemMatrixTimesOtherMatrix(IMatrixView otherMatrix)
        {
            throw new NotImplementedException();
        }

        public class Builder : ISolverBuilder
        {
            public AsymmetricDofOrderer RowDofOrderer { get; set; } = new AsymmetricDofOrderer(new RowDofOrderingStrategy());

            public IDofOrderer ColumnDofOrderer { get; set; } = new DofOrderer(new NodeMajorDofOrderingStrategy(), new NullReordering());

            public ISolver BuildSolver(IStructuralModel model)
            {
                if (!(model is IStructuralAsymmetricModel asymmetricModel))
                    throw new ArgumentException("Gmres solver builder can be used only with asymmetric models.");

                return new GmresSolver(asymmetricModel, RowDofOrderer, ColumnDofOrderer);
            }
        }
    }
}