using ISAAR.MSolve.LinearAlgebra.Matrices;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Overlapping.Schwarz.Additive.Interfaces
{
    public interface IModelOverlappingDecomposition
    {
        int NumberOfSubdomains { get; }
        void DecomposeMatrix();
        int[] GetConnectivityOfSubdomain(int indexSubdomain);
        CsrMatrix CoarseSpaceInterpolation { get; }
    }
}
