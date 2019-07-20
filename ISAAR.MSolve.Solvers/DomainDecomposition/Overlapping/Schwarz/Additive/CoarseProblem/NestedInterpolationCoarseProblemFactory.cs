using System.Collections.Generic;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.Solvers.DomainDecomposition.Overlapping.Schwarz.Additive.Interfaces;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Overlapping.Schwarz.Additive.CoarseProblem
{
    public class NestedInterpolationCoarseProblemFactory : ICoarseProblemFactory
    {
        private Matrix inverseCoarseMatrix;
        private IModelOverlappingDecomposition _modelOverlappingDecomposition;

        public void GenerateMatrices(IMatrixView matrix, IModelOverlappingDecomposition modelOverlappingDecomposition)
        {
            _modelOverlappingDecomposition = modelOverlappingDecomposition;
            var interpolationMatrix = modelOverlappingDecomposition.CoarseSpaceInterpolation.Transpose();
            var coarseSpaceMatrix = interpolationMatrix.ThisTimesOtherTimesThisTranspose(matrix);
            inverseCoarseMatrix = coarseSpaceMatrix.Invert();
        }

        public Matrix RetrievePreconditionerContribution()
        {
            var interpolationMatrix = _modelOverlappingDecomposition.CoarseSpaceInterpolation;
            return interpolationMatrix.ThisTimesOtherTimesThisTranspose(inverseCoarseMatrix);
        }
    }
}