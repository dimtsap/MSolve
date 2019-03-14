﻿using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.LinearAlgebra.Matrices;

//TODO: This should only calculate them. Another object should manage them.
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Feti
{
    /// <summary>
    /// Calculates the signed boolean matrices of the equations that enforce continuity between the multiple instances of 
    /// boundary dofs.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class ContinuityEquationsCalculator
    {
        private readonly ICrosspointStrategy crosspointStrategy;

        internal ContinuityEquationsCalculator(ICrosspointStrategy crosspointStrategy)
        {
            this.crosspointStrategy = crosspointStrategy;
        }

        public Dictionary<int, SignedBooleanMatrix> BooleanMatrices { get; private set; }
        public int NumContinuityEquations { get; private set; }
        public DiagonalMatrix LagrangeMultiplierMultiplicity { get; private set; } //TODO: so far this is useless

        public void CreateBooleanMatrices(IStructuralModel_v2 model)
        {
            // Find boundary nodes and continuity equations
            var boundaryNodes = 
                new List<(INode node, DOFType[] dofs, ISubdomain_v2[] subdomainsPlus, ISubdomain_v2[] subdomainsMinus)>();
            NumContinuityEquations = 0;
            foreach (INode node in model.Nodes) //TODO: this probably doesn't work if there are embedded nodes.
            {
                int nodeMultiplicity = node.SubdomainsDictionary.Count;
                if (nodeMultiplicity > 1)
                {
                    DOFType[] dofsOfNode = model.GlobalDofOrdering.GlobalFreeDofs.GetColumnsOfRow(node).ToArray(); //TODO: interacting with the table can be optimized.
                    
                    // If all dofs of this node are constrained, then it is not considered boundary.
                    if (dofsOfNode.Length == 0) continue;
                    
                    IEnumerable<ISubdomain_v2> nodeSubdomains = node.SubdomainsDictionary.Values;
                    (ISubdomain_v2[] subdomainsPlus, ISubdomain_v2[] subdomainsMinus) =
                        crosspointStrategy.FindSubdomainCombinations(nodeMultiplicity, nodeSubdomains);
                    
                    boundaryNodes.Add((node, dofsOfNode, subdomainsPlus, subdomainsMinus));
                    NumContinuityEquations += dofsOfNode.Length * subdomainsPlus.Length;
                }
            }

            // Create the signed boolean matrices.
            BooleanMatrices = new Dictionary<int, SignedBooleanMatrix>();
            foreach (ISubdomain_v2 subdomain in model.Subdomains)
            {
                BooleanMatrices.Add(subdomain.ID, 
                    new SignedBooleanMatrix(NumContinuityEquations, subdomain.FreeDofOrdering.NumFreeDofs));
            }

            var lagrangeMultiplicities = new double[NumContinuityEquations];

            // Fill the boolean matrices: node major, subdomain medium, dof minor. TODO: not sure about this order.
            int equationCounter = 0;
            foreach ((INode node, DOFType[] dofs, ISubdomain_v2[] subdomainsPlus, ISubdomain_v2[] subdomainsMinus) 
                in boundaryNodes)
            {
                double nodeMultiplicity = node.SubdomainsDictionary.Count;
                int numSubdomainCombos = subdomainsPlus.Length;
                for (int c = 0; c < numSubdomainCombos; ++c)
                {
                    //TODO: each subdomain appears in many combinations. It would be faster to cache its indices for the specific dofs.
                    SignedBooleanMatrix booleanPlus = BooleanMatrices[subdomainsPlus[c].ID];
                    SignedBooleanMatrix booleanMinus = BooleanMatrices[subdomainsMinus[c].ID];
                    bool nodeIsOrdered = subdomainsPlus[c].FreeDofOrdering.FreeDofs.TryGetDataOfRow(node,
                        out IReadOnlyDictionary<DOFType, int> dofsPlus);
                    Debug.Assert(nodeIsOrdered);
                    nodeIsOrdered = subdomainsMinus[c].FreeDofOrdering.FreeDofs.TryGetDataOfRow(node,
                        out IReadOnlyDictionary<DOFType, int> dofsMinus);
                    Debug.Assert(nodeIsOrdered);

                    foreach (DOFType dof in dofs)
                    {
                        booleanPlus.AddEntry(equationCounter, dofsPlus[dof], true);
                        booleanMinus.AddEntry(equationCounter, dofsMinus[dof], false);
                        lagrangeMultiplicities[equationCounter] = nodeMultiplicity;
                        ++equationCounter;
                    }
                }
            }

            LagrangeMultiplierMultiplicity = DiagonalMatrix.CreateFromArray(lagrangeMultiplicities);
        }
    }
}
