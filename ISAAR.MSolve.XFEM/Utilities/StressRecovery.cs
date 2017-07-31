﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Interpolation;
using ISAAR.MSolve.XFEM.LinearAlgebra;
using ISAAR.MSolve.XFEM.Tensors;

namespace ISAAR.MSolve.XFEM.Utilities
{
    class StressRecovery
    {
        private readonly Model2D model;

        public StressRecovery(Model2D model)
        {
            this.model = model;
        }

        public IReadOnlyList<Tensor2D> ComputeSmoothedNodalStresses(double[] solution)
        {
            var stressesFromAllElements = new Dictionary<XNode2D, List<Tensor2D>>();
            foreach (var node in model.Nodes) stressesFromAllElements[node] = new List<Tensor2D>();
            double[] constrainedDisplacements = model.CalculateConstrainedDisplacements();

            foreach (var element in model.Elements)
            {
                IReadOnlyDictionary<XNode2D, Tensor2D> elementStresses =
                    ComputeNodalStressesOfElement(element, solution, constrainedDisplacements);
                foreach (var nodeStressPair in elementStresses)
                {
                    stressesFromAllElements[nodeStressPair.Key].Add(nodeStressPair.Value);
                }
            }

            // Average with equal weights for all elements. TODO: perhaps vary the weights depending on the element type/area
            var nodalStresses = new Tensor2D[model.Nodes.Count];
            for (int i = 0; i < model.Nodes.Count; ++i)
            {
                XNode2D node = model.Nodes[i];
                double stressXX = 0.0, stressYY = 0.0, stressXY = 0.0;
                foreach (var tensor in stressesFromAllElements[node])
                {
                    stressXX += tensor.XX;
                    stressYY += tensor.YY;
                    stressXY += tensor.XY;
                }
                int contributingElementsCount = stressesFromAllElements[node].Count;
                stressXX /= contributingElementsCount;
                stressYY /= contributingElementsCount;
                stressXY /= contributingElementsCount;
                nodalStresses[i] = new Tensor2D(stressXX, stressYY, stressXY);
            }

            return nodalStresses;
        }

        // Computes stresses directly at the nodes. The other approach is to compute them at Gauss points and then extrapolate
        private IReadOnlyDictionary<XNode2D, Tensor2D> ComputeNodalStressesOfElement(XContinuumElement2D element, 
            double[] freeDisplacements, double[] constrainedDisplacements)
        {
            double[] standardDisplacements = model.DofEnumerator.ExtractDisplacementVectorOfElementFromGlobal(element, 
                freeDisplacements, constrainedDisplacements);
            double[] enrichedDisplacements = 
                model.DofEnumerator.ExtractEnrichedDisplacementsOfElementFromGlobal(element, freeDisplacements);

            IReadOnlyList<INaturalPoint2D> naturalNodes = element.ElementType.NaturalCoordinatesOfNodes;
            var nodalStresses = new Dictionary<XNode2D, Tensor2D>();
            for (int i = 0; i < element.Nodes.Count; ++i)
            {
                EvaluatedInterpolation2D evaluatedInterpolation =
                    element.Interpolation.EvaluateAt(element.Nodes, naturalNodes[i]);
                DenseMatrix displacementGradient = element.CalculateDisplacementFieldGradient(
                    naturalNodes[i], evaluatedInterpolation, standardDisplacements, enrichedDisplacements);
                Matrix2D constitutive = 
                    element.Material.CalculateConstitutiveMatrixAt(naturalNodes[i], evaluatedInterpolation);
                nodalStresses[element.Nodes[i]] = element.CalculateStressTensor(displacementGradient, constitutive);
            }

            return nodalStresses;
        }
    }
}
