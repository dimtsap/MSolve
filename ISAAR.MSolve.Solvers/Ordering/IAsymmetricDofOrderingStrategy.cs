using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;

namespace ISAAR.MSolve.Solvers.Ordering
{
	public interface IAsymmetricDofOrderingStrategy
	{
		(int numGlobalFreeDofs, DofTable globalFreeDofs) OrderGlobalDofs(IStructuralAsymmetricModel model);

		(int numSubdomainFreeDofs, DofTable subdomainFreeDofs) OrderSubdomainDofs(IAsymmetricSubdomain subdomain);
	}
}