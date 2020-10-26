using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.OAS.StiffnessMatrices
{
    public interface IOasMatrixManager
    {
        Vector SolveWithSubdomainMatrix(int subdomainId, Vector rhsVector);
    }
}