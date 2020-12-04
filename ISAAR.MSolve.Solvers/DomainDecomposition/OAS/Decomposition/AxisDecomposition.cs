using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.OAS.Decomposition
{
    public class AxisDecomposition
    {
        private readonly int numberOfControlPoints;
        private int cpPerSubdomain;
        private List<int[]> subdomainCpIndices=new List<int[]>();

        public AxisDecomposition(int numberOfControlPoints, int numberOfSubdomains, int overlapping=1)
        {
            this.numberOfControlPoints = numberOfControlPoints;
            if ((numberOfControlPoints + (numberOfSubdomains - 1) * overlapping) % numberOfSubdomains == 0)
            {
                cpPerSubdomain=
                    (numberOfControlPoints + (numberOfSubdomains - 1) * overlapping) / numberOfSubdomains;
            }
            else
            {
                cpPerSubdomain = (numberOfControlPoints + (numberOfSubdomains - 1) * overlapping) / numberOfSubdomains + 1;
            }

            for (int i = 0; i < numberOfSubdomains; i++)
            {
                var subdomainStartingCpId = i * cpPerSubdomain - i* overlapping;
                var subdomainEndingCpId = subdomainStartingCpId + cpPerSubdomain - 1;

                if (subdomainStartingCpId > numberOfControlPoints-1)
                {
                    subdomainStartingCpId = numberOfControlPoints-1;
                }
                if (subdomainEndingCpId > numberOfControlPoints-1)
                {
                    subdomainEndingCpId = numberOfControlPoints-1;
                }

                if (!(subdomainStartingCpId == numberOfControlPoints && subdomainEndingCpId == numberOfControlPoints))
                {
                    subdomainCpIndices.Add(Enumerable.Range(subdomainStartingCpId, subdomainEndingCpId - subdomainStartingCpId+1)
                        .ToArray());
                }
            }
        }

        public int GetNumberOfSubdomains => subdomainCpIndices.Count;

        public int[] GetIndicesOfSubdomainCp(int subdomainNumber) => subdomainCpIndices[subdomainNumber];

        public int GetNumberOfControlPoints => numberOfControlPoints;
    }
}
