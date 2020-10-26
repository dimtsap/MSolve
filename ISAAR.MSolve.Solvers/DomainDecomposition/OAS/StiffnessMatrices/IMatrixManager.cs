using System.Collections.Generic;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.LinearSystems;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.OAS.StiffnessMatrices
{
    public interface IMatrixManager
    {
        ILinearSystem GetLinearSystem(int subdomainID);

        IMatrix BuildKff(int subdomainID, ISubdomainFreeDofOrdering dofOrdering, IEnumerable<IElement> elements,
            IElementMatrixProvider matrixProvider);

        //TODO: Refactor this
        void SetSolution(int subdomainID, Vector solution);
    }
}