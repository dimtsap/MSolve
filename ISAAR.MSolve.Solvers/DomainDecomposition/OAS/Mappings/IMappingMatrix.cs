using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.OAS.Mappings
{
    public interface IMappingMatrix
    {
        int NumColumns { get; }

        int NumRows { get; }

        Matrix CopyToFullMatrix();

        Vector Multiply(Vector vector, bool transpose);
    }
}