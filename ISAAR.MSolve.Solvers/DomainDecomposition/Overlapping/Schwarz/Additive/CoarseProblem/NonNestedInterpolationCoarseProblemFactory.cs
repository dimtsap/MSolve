using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Matrices;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Overlapping.Schwarz.Additive.CoarseProblem
{
    public class NonNestedInterpolationCoarseProblemFactory : ICoarseProblemFactory
    {
        public void GenerateMatrices(IOverlappingDecomposition overlappingDecomposition)
        {
            throw new NotImplementedException();
        }

        public CsrMatrix RetrievePreconditionerContribution()
        {
            throw new NotImplementedException();
        }
    }
}
