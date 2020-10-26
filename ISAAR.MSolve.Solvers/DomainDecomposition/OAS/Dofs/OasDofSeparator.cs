using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Solvers.DomainDecomposition.OAS.Decomposition;
using ISAAR.MSolve.Solvers.DomainDecomposition.OAS.Mappings;
using ISAAR.MSolve.Solvers.DomainDecomposition.Overlapping.Schwarz.Additive;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.OAS.Dofs
{
    public class OasDofSeparator:IOasDofSeparator
    {
        private readonly IStructuralAsymmetricModel model;
        private readonly int numberOfSubdomainsX;
        private readonly int numberOfSubdomainsY;
        private readonly Dictionary<int, DofTable> subdomainDofOrderingsLocal=new Dictionary<int, DofTable>();
        private readonly Dictionary<int ,int[]> subdomainDofsLocalToFree=new Dictionary<int, int[]>();
        private readonly Dictionary<int, BooleanMatrixRowsToColumns> subdomainToFreeDofMappings=new Dictionary<int, BooleanMatrixRowsToColumns>();

        public OasDofSeparator(IStructuralAsymmetricModel model, int numberOfSubdomainsX, int numberOfSubdomainsY)
        {
            this.model = model;
            this.numberOfSubdomainsX = numberOfSubdomainsX;
            this.numberOfSubdomainsY = numberOfSubdomainsY;
        }

        public int NumberOfSubdomains => numberOfSubdomainsX * numberOfSubdomainsY;

        public BooleanMatrixRowsToColumns GetDofMappingBoundaryClusterToSubdomain(int subdomainID) =>
            subdomainToFreeDofMappings[subdomainID];

        public int[] GetDofsSubdomainToFree(int subdomainId) => subdomainDofsLocalToFree[subdomainId];

        public int GetNumFreeDofs(int subdomainID) => subdomainDofOrderingsLocal[subdomainID].EntryCount;

        public void SeparateSubdomainAndCoarseProblemDofs()
        {
            var ksiDecomposition = new AxisDecomposition(model.Subdomains[0].NumberOfCpKsi, numberOfSubdomainsX);
            var hetaDecomposition = new AxisDecomposition(model.Subdomains[0].NumberOfCpHeta, numberOfSubdomainsY);
            var overlappingDecomposition = new OverlappingDecomposition2D(ksiDecomposition, hetaDecomposition);

            var freeDofs = model.Subdomains[0].FreeDofRowOrdering.FreeDofs;
            
            for (int i = 0; i < overlappingDecomposition.GetNumberOfSubdomains; i++)
            {
                var (subdomainDofOrdering, localToGlobal) = OrderSubdomainDof(overlappingDecomposition, i, freeDofs);

                subdomainDofOrderingsLocal[i] = subdomainDofOrdering;
                subdomainDofsLocalToFree[i] = localToGlobal.ToArray();
            }

            for (int i = 0; i < overlappingDecomposition.GetNumberOfSubdomains; i++)
            {
                var numberOfFreeDofs = freeDofs.EntryCount;
                var numberOfSubdomainDofs = subdomainDofOrderingsLocal[i].EntryCount;
                subdomainToFreeDofMappings[i] = new BooleanMatrixRowsToColumns(numberOfSubdomainDofs, numberOfFreeDofs,
                    subdomainDofsLocalToFree[i]);
            }
        }

        private (DofTable subdomainDofOrdering, List<int> localToGlobal) OrderSubdomainDof(
            OverlappingDecomposition2D overlappingDecomposition, int i, DofTable freeDofs)
        {
            var subdomainDofIdx = 0;
            var subdomainDofOrdering = new DofTable();
            var localToGlobal = new List<int>();
            var cpIndices = overlappingDecomposition.GetControlPointIndicesOfSubdomain(i);
            var freeControlPoints = freeDofs.GetRows();

            for (int j = 0; j < cpIndices.Length; j++)
            {
                var node=freeControlPoints.SingleOrDefault(x => x.ID == cpIndices[j]);
                if (node==null) continue;
                var dofsOfNode = freeDofs.GetDataOfRow(node);
                foreach (var dofTypeIdxPair in dofsOfNode)
                {
                    subdomainDofOrdering[node, dofTypeIdxPair.Key] = subdomainDofIdx++;
                    localToGlobal.Add(dofTypeIdxPair.Value);
                }
            }

            return (subdomainDofOrdering, localToGlobal);
        }
    }
}
