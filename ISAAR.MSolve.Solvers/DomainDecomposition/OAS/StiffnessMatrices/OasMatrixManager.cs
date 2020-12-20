using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.LinearAlgebra;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;
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
        private LUCSparseNet invA0;

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
            //var A0 = ExtractCoarseProblemMatrix(Kff);
            //invA0 = LUCSparseNet.Factorize(A0);
        }

        //private CscMatrix ExtractCoarseProblemMatrix(CsrMatrix Kff)
        //{
        //    // Global To Local Interpolation
        //    //var R0 = dofSeparator.GetGlobalToCoarseMapping;
        //    //var aux1 = Kff.MultiplyLeft(R0, false, true);
        //    //var A0 = aux1.MultiplyRight(R0, false);

        //    // Local to Global Interpolation

            
        //     var R0 = dofSeparator.GetGlobalToCoarseMapping;

        //     LibrarySettings.LinearAlgebraProviders = LinearAlgebraProviderChoice.MKL;
        //     var A0=MultiplyR0TimesKffTimesR0Transpose(R0, Kff);
        //     LibrarySettings.LinearAlgebraProviders = LinearAlgebraProviderChoice.Managed;
             
             
        //     // var A0=R0.ThisTimesOtherTimesThisTranspose(Kff);
             
        //     //var aux1 = R0.MultiplyRight(Kff,false, false);
        //    // var aux1 = Kff.MultiplyLeft(R0, false,false); // this is the last working one 
        //    //
        //    // var A0 = R0.MultiplyLeft(aux1,true); // this is the last working one
        //    //var A0 = aux1.MultiplyRight(R0, false, true);

        //    var dokCol = DokColMajor.CreateEmpty(A0.NumRows, A0.NumColumns);
        //    for (int i = 0; i < A0.NumRows; i++)
        //    {
        //        for (int j = 0; j < A0.NumColumns; j++)
        //        {
        //            if (Math.Abs(A0[i, j]) < 1e-16) continue;
        //            dokCol[i, j] = A0[i, j];
        //        }
        //    }

        //    return dokCol.BuildCscMatrix(true);
        //}

        private Matrix MultiplyR0TimesKffTimesR0Transpose(CscMatrix r0, CsrMatrix kff)
        {
            var coarseProblemSize = r0.NumRows;
            var globalProblemSize = r0.NumColumns;
            var columnsOfAux = new Vector[coarseProblemSize];
            for (int i = 0; i < coarseProblemSize; i++)
            {
                var columnOfTranspose=r0.GetRow(i);
                columnsOfAux[i]=kff.Multiply(columnOfTranspose);
            }

            var matrix = Matrix.CreateZero(coarseProblemSize, coarseProblemSize);
            for (int i = 0; i < coarseProblemSize; i++)
            {
               matrix.SetSubcolumn(i,r0.Multiply(columnsOfAux[i])); 
            }

            return matrix;
        }
        
        private Matrix MultiplyR0TimesKffTimesR0Transpose(CsrMatrix r0, CsrMatrix kff)
        {
            var coarseProblemSize = r0.NumRows;
            var globalProblemSize = r0.NumColumns;
            var columnsOfAux = new Vector[coarseProblemSize];
            for (int i = 0; i < coarseProblemSize; i++)
            {
                var columnOfTranspose=r0.GetRow(i);
                columnsOfAux[i]=kff.Multiply(columnOfTranspose);
            }

            var matrix = Matrix.CreateZero(coarseProblemSize, coarseProblemSize);
            for (int i = 0; i < coarseProblemSize; i++)
            {
                matrix.SetSubcolumn(i,r0.Multiply(columnsOfAux[i])); 
            }

            return matrix;
        }

        public Vector SolveCoarseProblemLinearSystem(Vector rhsVector)
        {
            var solution = Vector.CreateZero(rhsVector.Length);
            invA0.SolveLinearSystem(rhsVector, solution);
            return solution;
        }



        public Vector SolveSubdomainLinearSystem(Vector rhsVector, int subdomainId)
        {
            var solution = Vector.CreateZero(rhsVector.Length);
            invAkl[subdomainId].SolveLinearSystem(rhsVector, solution);
            return solution;
        }

    }
}
