using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Overlapping.Schwarz.Additive
{
    public class ParametricAxisDecomposition:IOverlappingDecomposition
    {
        private List<double> _decompositionIndices;
        public ParametricAxisDecomposition(IVector knotValueVector,int degree, int numberOfProcessors)
        {
            var numberOfShapeFunctions = knotValueVector.Length - degree - 1;
            var shapeFunctionSupports= new List<FunctionSupport>();
            var indexFunction = 0;
            for (int i = 0; i < numberOfShapeFunctions; i++)
            {
                shapeFunctionSupports.Add(new FunctionSupport
                {
                    ID = indexFunction++,
                    Start = knotValueVector[i],
                    End = knotValueVector[i+degree+1]
                });
            }

            var knots = knotValueVector.CopyToArray().Distinct().ToList();
            var axisElements = knots.Count() - 1;
            var subdomainElements = (axisElements % numberOfProcessors == 0)
                ? axisElements / numberOfProcessors
                : axisElements / numberOfProcessors + 1;

            _decompositionIndices = new List<double>();
            for (int i = 0; i < knots.Count(); i+=subdomainElements)
                _decompositionIndices.Add(knots[i]);
            _decompositionIndices.Add(knots[knots.Count-1]);

            var shapeFunctionsIds=new List<int>();
            shapeFunctionsIds.Add(0);
            foreach (var index in _decompositionIndices)
            {
                var id = shapeFunctionSupports.First(s => s.Start < index && s.End > index).ID;
                shapeFunctionsIds.Add(id);
            }

            var subdomainIndices = new int[shapeFunctionsIds.Count-1][];
            for (int i = 0; i < shapeFunctionsIds.Count-2; i++)
            {
                var indexCount = shapeFunctionsIds[i + 1] - shapeFunctionsIds[i];
                subdomainIndices[i] = Enumerable.Range(shapeFunctionsIds[i], indexCount).ToArray();
            }
            
        }
    }

    public class FunctionSupport
    {
        public int ID { get; set; }
        public double Start { get; set; }
        public double End { get; set; }
    }
}
