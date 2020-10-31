using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.Discretization.Interfaces
{
	public interface IAsymmetricSubdomain:ISubdomain
	{
		ISubdomainFreeDofOrdering FreeDofRowOrdering { get; set; }
		ISubdomainFreeDofOrdering FreeDofColOrdering { get; set; }

        ISubdomainConstrainedDofOrdering ConstrainedDofRowOrdering { get; set; }
        ISubdomainConstrainedDofOrdering ConstrainedDofColOrdering { get; set; }

		int NumberOfCpKsi { get; }
        int NumberOfCpHeta { get; }

        Vector KnotValueVectorKsi { get; set; }
        Vector KnotValueVectorHeta { get; set; }

        int DegreeKsi { get; set; }
        int DegreeHeta { get; set; }
    }
}