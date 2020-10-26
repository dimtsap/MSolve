using ISAAR.MSolve.Solvers.DomainDecomposition.OAS.Mappings;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.OAS.Dofs
{
    public interface IOasDofSeparator
    {
        int NumberOfSubdomains { get; }

        int GetNumFreeDofs(int subdomainID);

        BooleanMatrixRowsToColumns GetDofMappingBoundaryClusterToSubdomain(int subdomainID);

        void SeparateSubdomainAndCoarseProblemDofs();
    }
}