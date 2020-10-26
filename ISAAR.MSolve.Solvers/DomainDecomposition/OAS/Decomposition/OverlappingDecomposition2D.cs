using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.OAS.Decomposition
{
    public class OverlappingDecomposition2D
    {
        private readonly AxisDecomposition _ksiDecomposition;
        private readonly AxisDecomposition _hetaDecomposition;
        private List<int[]> subdomainControlPointIndices2D=new List<int[]>();
        public OverlappingDecomposition2D(AxisDecomposition ksiDecomposition, AxisDecomposition hetaDecomposition)
        {
            _ksiDecomposition = ksiDecomposition;
            _hetaDecomposition = hetaDecomposition;

            var numberOfCpHeta = _hetaDecomposition.GetNumberOfControlPoints;
            for (int i = 0; i < _ksiDecomposition.GetNumberOfSubdomains; i++)
            {
                var axisKsiIndices = _ksiDecomposition.GetIndicesOfSubdomainCp(i);
                var numberOfIndicesKsi = axisKsiIndices.Length;
                for (int j = 0; j < _hetaDecomposition.GetNumberOfSubdomains; j++)
                {
                    var axisHetaIndices = _hetaDecomposition.GetIndicesOfSubdomainCp(j);
                    var numberOfIndicesHeta = axisKsiIndices.Length;

                    var indices2D = new int[numberOfIndicesKsi*numberOfIndicesHeta];
                    var counter = 0;

                    for (var indexKsi = 0; indexKsi < axisKsiIndices.Length; indexKsi++)
                    {
                        for (int indexHeta = 0; indexHeta < axisHetaIndices.Length; indexHeta++)
                        {
                            indices2D[counter++] = axisKsiIndices[indexKsi] * numberOfCpHeta + axisHetaIndices[indexHeta];
                        }
                    }
                    subdomainControlPointIndices2D.Add(indices2D);
                }
            }
        }

        public int GetNumberOfSubdomains => subdomainControlPointIndices2D.Count;

        public int[] GetControlPointIndicesOfSubdomain(int subdomainIndex) =>
            subdomainControlPointIndices2D[subdomainIndex];
    }
}
