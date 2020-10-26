using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Triangulation;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.DomainDecomposition.OAS.Dofs;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.OAS.StiffnessMatrices
{
    public class OasMatrixManager
    {
        private readonly OasDofSeparator dofSeparator;
        private Dictionary<int, CscMatrix> Akl = new Dictionary<int, CscMatrix>();
        private Dictionary<int, LUCSparseNet> invAkl = new Dictionary<int, LUCSparseNet>();
        //private Dictionary<int, LUFactorization> invAkl = new Dictionary<int, LUFactorization>(); 

        public OasMatrixManager(OasDofSeparator dofSeparator)
        {
            this.dofSeparator = dofSeparator;
        }


        public void ExtractSubdomainMatrices(IMatrixView matrix)
        {
            var numberOfSubdomains = dofSeparator.NumberOfSubdomains;
            var Kff = (CsrMatrix) matrix;
            var globalFreeDofs = Enumerable.Range(0, Kff.NumRows).ToArray();

            for (int i = 0; i < numberOfSubdomains; i++)
            {
                var subdomainDofs = dofSeparator.GetDofsSubdomainToFree(i);
                var subdomainDofsCount = subdomainDofs.Length;

                var remainingDofs = globalFreeDofs.Except(subdomainDofs).ToArray();
                var submatrixExtractor= new SubmatrixExtractorFullCsrCsc();

                submatrixExtractor.ExtractSubmatrices(Kff,remainingDofs,subdomainDofs);
                //TODO: This might need better handling
                var Akl = submatrixExtractor.Submatrix11;

                invAkl[i] = LUCSparseNet.Factorize(Akl);
            }
        }

        public Vector SolveSubdomainLinearSystem(Vector rhsVector, int subdomainId)
        {
            var solution = Vector.CreateZero(rhsVector.Length);
            invAkl[subdomainId].SolveLinearSystem(rhsVector, solution);
            return solution;
        }

    }
}
