using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.DomainDecomposition.Overlapping.Schwarz.Additive;
using Xunit;

namespace ISAAR.MSolve.Solvers.Tests
{
    public class OverlappingSchwarzTests
    {
        [Fact]
        public void ParametricAxisDecompositionTest()
        {
            var knotValueVector = Vector.CreateFromArray(new double[]
                {0, 0, 0, 0, 1 / 6.0, 1 / 3.0, 1 / 2.0, 2 / 3.0, 5 / 6.0, 1, 1, 1, 1});

            var d = new ParametricAxisDecomposition(knotValueVector, 3, 2);
        }
    }
}
