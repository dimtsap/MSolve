using System;
using System.Diagnostics;
using System.IO;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Iterative;
using ISAAR.MSolve.LinearAlgebra.Iterative.GeneralizedMinimalResidual;
using ISAAR.MSolve.LinearAlgebra.Iterative.Preconditioning;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.Solvers.Assemblers;
using ISAAR.MSolve.Solvers.Commons;
using ISAAR.MSolve.Solvers.Ordering;
using ISAAR.MSolve.Solvers.Ordering.Reordering;

namespace ISAAR.MSolve.Solvers.Iterative
{
    public class GmresSolver : SingleSubdomainNonSymmetricSolverBase<CsrMatrix>
    {
        private readonly GmresAlgorithm gmresAlgorithm;
        private readonly IPreconditionerFactory preconditionerFactory;

        private bool mustUpdatePreconditioner = true;
        private IPreconditioner preconditioner;

        public GmresSolver(IStructuralAsymmetricModel model, GmresAlgorithm gmresAlgorithm,
            IPreconditionerFactory preconditionerFactory, AsymmetricDofOrderer dofRowOrderer, IDofOrderer dofColOrderer)
            : base(model, dofRowOrderer, dofColOrderer, new CsrNonSymmetricAssembler(true), "GmresSolver")
        {
            this.gmresAlgorithm = new GmresAlgorithm.Builder().Build();
            this.preconditionerFactory = preconditionerFactory;
        }

        public override void Initialize()
        {
        }

        public override void HandleMatrixWillBeSet()
        {
            mustUpdatePreconditioner = true;
            preconditioner = null;
        }

        public override void PreventFromOverwrittingSystemMatrices()
        {
        }

        public override void Solve()
        {
            var watch = new Stopwatch();
            if (linearSystem.SolutionConcrete == null) linearSystem.SolutionConcrete = linearSystem.CreateZeroVectorConcrete();
            else linearSystem.SolutionConcrete.Clear();

            // Preconditioning
            if (mustUpdatePreconditioner)
            {
                watch.Start();
                preconditioner = preconditionerFactory.CreatePreconditionerFor(linearSystem.Matrix);
                watch.Stop();
                Logger.LogTaskDuration("Calculating preconditioner", watch.ElapsedMilliseconds);
                watch.Reset();
                mustUpdatePreconditioner = false;
            }

            watch.Start();
            IterativeStatistics stats = gmresAlgorithm.Solve(linearSystem.Matrix, preconditioner, linearSystem.RhsConcrete,
                linearSystem.SolutionConcrete, true, () => linearSystem.CreateZeroVector());
            if (!stats.HasConverged)
                throw new IterativeSolverNotConvergedException("Gmres did not converge");
            watch.Stop();
            using (StreamWriter writer = new StreamWriter(Path.Combine(Directory.GetCurrentDirectory(),"OasGmres.txt")))
            {
                writer.WriteLine($"Converged: {stats.HasConverged}");
                writer.WriteLine($"Gmres iterations: {stats.NumIterationsRequired}");
                writer.WriteLine($"Gmres residual norm: {stats.ResidualNormRatioEstimation}");
            }
            Logger.LogTaskDuration("Iterative algorithm", watch.ElapsedMilliseconds);
            Logger.LogIterativeAlgorithm(stats.NumIterationsRequired, stats.ResidualNormRatioEstimation);
            Logger.IncrementAnalysisStep();
        }

        protected override Matrix InverseSystemMatrixTimesOtherMatrix(IMatrixView otherMatrix)
        {
            throw new NotImplementedException();
        }

        public class Builder : ISolverBuilder
        {
            public AsymmetricDofOrderer RowDofOrderer { get; set; } =
                new AsymmetricDofOrderer(new RowDofOrderingStrategy());

            public IDofOrderer ColumnDofOrderer { get; set; } =
                new DofOrderer(new NodeMajorDofOrderingStrategy(), new NullReordering());

            public GmresAlgorithm GmresAlgorithm { get; set; } = (new GmresAlgorithm.Builder()).Build();

            public IPreconditionerFactory PreconditionerFactory { get; set; } = new IdentityPreconditioner.Factory();

            public ISolver BuildSolver(IStructuralModel model)
            {
                if (!(model is IStructuralAsymmetricModel asymmetricModel))
                    throw new ArgumentException("Gmres solver builder can be used only with asymmetric models.");

                return new GmresSolver(asymmetricModel, GmresAlgorithm,PreconditionerFactory, RowDofOrderer, ColumnDofOrderer);
            }
        }
    }
}