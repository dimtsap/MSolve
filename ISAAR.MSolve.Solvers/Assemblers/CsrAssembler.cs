﻿using System.Collections.Generic;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.Solvers.Commons;
using ISAAR.MSolve.Solvers.Ordering;

//TODO: this should work with MatrixProviders instead of asking the elements directly.
namespace ISAAR.MSolve.Solvers.Assemblers
{
    /// <summary>
    /// Builds the global matrix of the linear system that will be solved. This matrix is in CSR format.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class CsrAssembler
    {
        private readonly bool sortColsOfEachRow;

        /// <summary>
        /// 
        /// </summary>
        /// <param name="sortColsOfEachRow">Sorting the columns of each row in the CSR storage format may increase performance 
        ///     of the matrix vector multiplications. It is recommended to set it to true, especially for iterative linear 
        ///     system solvers.</param>
        public CsrAssembler(bool sortColsOfEachRow = true)
        {
            this.sortColsOfEachRow = sortColsOfEachRow;
        }

        public (CsrMatrix Kff, DokRowMajor Kfc) BuildGlobalMatrices(IEnumerable<ContinuumElement2D> elements,
            AllDofOrderer dofOrderer)
        {
            int numConstrainedDofs = dofOrderer.NumConstrainedDofs;
            int numFreeDofs = dofOrderer.NumFreeDofs;
            var Kff = DokRowMajor.CreateEmpty(numFreeDofs, numFreeDofs);
            var Kfc = DokRowMajor.CreateEmpty(numFreeDofs, numConstrainedDofs);

            foreach (ContinuumElement2D element in elements)
            {
                // TODO: perhaps that could be done and cached during the dof enumeration to avoid iterating over the dofs twice
                (IReadOnlyDictionary<int, int> mapStandard, IReadOnlyDictionary<int, int> mapConstrained) =
                    dofOrderer.MapDofsElementToGlobal(element);
                Matrix k = Conversions.MatrixOldToNew(element.BuildStiffnessMatrix());
                Kff.AddSubmatrixSymmetric(k, mapStandard);
                Kfc.AddSubmatrix(k, mapStandard, mapConstrained);
            }

            //TODO: perhaps I should filter the matrices in the concrete class before returning (e.g. dok.Build())
            return (Kff.BuildCsrMatrix(sortColsOfEachRow), Kfc);
        }

        public CsrMatrix BuildGlobalMatrix(IEnumerable<ContinuumElement2D> elements, FreeDofOrderer dofOrderer)
        {
            int numFreeDofs = dofOrderer.NumFreeDofs;
            var Kff = DokRowMajor.CreateEmpty(numFreeDofs, numFreeDofs);

            foreach (ContinuumElement2D element in elements)
            {
                // TODO: perhaps that could be done and cached during the dof enumeration to avoid iterating over the dofs twice
                IReadOnlyDictionary<int, int> mapStandard = dofOrderer.MapFreeDofsElementToGlobal(element);
                Matrix k = Conversions.MatrixOldToNew(element.BuildStiffnessMatrix());
                Kff.AddSubmatrixSymmetric(k, mapStandard);
            }

            return Kff.BuildCsrMatrix(sortColsOfEachRow);
        }
    }
}
