using System.Collections.Generic;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.Solvers.DomainDecomposition.Overlapping.Schwarz.Additive.Interfaces;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Overlapping.Schwarz.Additive.CoarseProblem
{
    public class NestedStiffnessCoarseProblemFactory : ICoarseProblemFactory
    {
        public void GenerateMatrices(IMatrixView matrix, IModelOverlappingDecomposition modelOverlappingDecomposition)
        {
            throw new System.NotImplementedException();
        }

        public Matrix RetrievePreconditionerContribution()
        {
            throw new System.NotImplementedException();
        }
    }
}
