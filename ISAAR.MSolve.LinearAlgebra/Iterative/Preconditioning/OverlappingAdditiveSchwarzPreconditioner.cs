using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.LinearAlgebra.Iterative.Preconditioning
{
    /// <summary>
    /// Implements the Overlapping Additive Schwarz for a square matrix.
    /// Authors: Dimitris Tsapetis
    /// </summary>
    public class OverlappingAdditiveSchwarzPreconditioner:IPreconditioner
    {
        public void SolveLinearSystem(IVectorView rhsVector, IVector lhsVector)
        {
            throw new NotImplementedException();
        }
    }


    public class Factory
    {

    }
}
