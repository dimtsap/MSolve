using ISAAR.MSolve.LinearAlgebra.Matrices;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Overlapping.Schwarz.Additive.CoarseProblem
{
    public interface ICoarseProblemFactory
    {
        void GenerateMatrices(IOverlappingDecomposition overlappingDecomposition);
        CsrMatrix RetrievePreconditionerContribution();
    }
}