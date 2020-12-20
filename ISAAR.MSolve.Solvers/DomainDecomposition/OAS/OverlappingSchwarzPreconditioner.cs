using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Iterative.Preconditioning;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.DomainDecomposition.OAS.Dofs;
using ISAAR.MSolve.Solvers.DomainDecomposition.OAS.StiffnessMatrices;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.OAS
{
    public class OverlappingSchwarzPreconditioner:IPreconditioner
    {
        private readonly OasDofSeparator _dofSeparator;
        private readonly OasMatrixManager _matrixManager;

        public OverlappingSchwarzPreconditioner(OasDofSeparator dofSeparator, OasMatrixManager matrixManager)
        {
            _dofSeparator = dofSeparator;
            _matrixManager = matrixManager;
        }

        public void SolveLinearSystem(IVectorView rhsVector, IVector lhsVector)
        {
            var numberOfSubdomains = _dofSeparator.NumberOfSubdomains;

            for (int i = 0; i < numberOfSubdomains; i++)
            {
                var Rkl = _dofSeparator.GetDofMappingBoundaryClusterToSubdomain(i);
                var r1 = Rkl.Multiply((Vector)rhsVector, false);
                var r2 = _matrixManager.SolveSubdomainLinearSystem(r1, i);
                var r3 = Rkl.Multiply(r2, true);
                lhsVector.AddIntoThis(r3);
            }

            //var R0 = _dofSeparator.GetGlobalToCoarseMapping;
            //var r01 = R0.Multiply((Vector) rhsVector, true);
            //var r02 = _matrixManager.SolveCoarseProblemLinearSystem(r01);
            //var r03 = R0.Multiply(r02, false);

            //var R0 = _dofSeparator.GetGlobalToCoarseMapping;
            //var r01 = R0.Multiply((Vector) rhsVector, false);
            //var r02 = _matrixManager.SolveCoarseProblemLinearSystem(r01);
            //var r03 = R0.Multiply(r02, true);


            //lhsVector.AddIntoThis(r03);
        }

        public class Factory: IPreconditionerFactory
        {
            private readonly IStructuralAsymmetricModel model;

            private readonly OasDofSeparator DofSeparator;

            private readonly OasMatrixManager MatrixManager;

            /// <summary>
            /// Initializes a new instance of <see cref="IdentityPreconditioner.Factory"/>.
            /// </summary>
            public Factory(IStructuralAsymmetricModel model, int numberOfSubdomainsX, int numberOfSubdomainsY, int overlapping=1)
            {
                this.model = model;
                this.DofSeparator=new OasDofSeparator(model,numberOfSubdomainsX, numberOfSubdomainsY, overlapping);
                this.MatrixManager= new OasMatrixManager(DofSeparator);
            }

            /// <summary>
            /// See <see cref="IPreconditionerFactory.CreatePreconditionerFor(IMatrixView)"/>.
            /// </summary>
            public IPreconditioner CreatePreconditionerFor(IMatrixView matrix)
            {
                DofSeparator.SeparateSubdomainAndCoarseProblemDofs();
                MatrixManager.ExtractSubdomainMatrices(matrix);
                var preconditioner= new OverlappingSchwarzPreconditioner(DofSeparator, MatrixManager);
                return preconditioner;
            }
        }
    }
}
