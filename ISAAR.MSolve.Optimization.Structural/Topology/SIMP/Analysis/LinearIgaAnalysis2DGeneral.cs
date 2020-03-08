using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Discretization.Providers;
using ISAAR.MSolve.IGA.Elements.Continuum;
using ISAAR.MSolve.IGA.Entities;
using ISAAR.MSolve.IGA.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers;

namespace ISAAR.MSolve.Optimization.Structural.Topology.SIMP.Analysis
{
    public class LinearIgaAnalysis2DGeneral:ILinearAnalysis
    {
        private readonly Model model;
        private readonly Patch subdomain;
        private readonly ISolver solver;
        private readonly IList<Element> modelElements;
        private readonly IElement[] elementWrappers;
        private readonly ScalingDofEnumerator[] penalizers;
        private readonly IElementMatrixProvider problemMatrixProvider = new ElementStructuralStiffnessProvider();
        private bool isFirstAnalysis = true;
        private ProblemStructural problem;
        private readonly ContinuumElement2D[] elementTypes;
        private LinearAnalyzer linearAnalyzer;
        private StaticAnalyzer staticAnalyzer;
        private Vector elementVolumes;
        private IVectorView globalDisplacements;

        public LinearIgaAnalysis2DGeneral(Model model, ISolver solver)
        {
            if (model.PatchesDictionary.Count != 1) throw new NotImplementedException(
                "Cannot perform topology optimization with multiple subdomains yet.");
            this.model = model;
            subdomain = model.PatchesDictionary.First().Value;
            this.solver = solver;

            // Elements
            modelElements = model.Elements;
            NumElements = modelElements.Count;
            elementTypes = new ContinuumElement2D[NumElements];
            elementWrappers = new IElement[NumElements];
            penalizers = new ScalingDofEnumerator[NumElements];
            for (int e = 0; e < NumElements; ++e)
            {
                if (modelElements[e].ElementType is ContinuumElement2D continuum) elementTypes[e] = continuum;
                else throw new ArgumentException("2D topology optimization only works with 2D continuum elements,"
                                                 + $" but the element with ID = {modelElements[e].ID} was not.");
                elementWrappers[e] = modelElements[e];
                penalizers[e] = new ScalingDofEnumerator();
            }
            NumLoadCases = 1;
        }

        public int NumElements { get; }

        public int NumLoadCases { get; }

        public double CalculateTotalVolume(IVectorView densities) => elementVolumes.DotProduct(densities);

        public Vector GetElementDisplacements(int elementIdx, int loadCaseIdx)
            => Vector.CreateFromArray(subdomain.FreeDofOrdering.ExtractVectorElementFromSubdomain(
                elementWrappers[elementIdx], globalDisplacements));

        public IMatrixView GetElementStiffnessForCurrentMaterial(int elementIdx)
            => elementTypes[elementIdx].StiffnessMatrix(modelElements[elementIdx]);

        public IMatrixView GetElementStiffnessForUnitMaterial(int elementIdx)
        {
            //TODO: It would be nice to have a GetOriginalMatrix() method in the dofEnumerator class and a way to call it.
            double currentModulus = penalizers[elementIdx].ScaleFactor;
            penalizers[elementIdx].ScaleFactor = 1.0;
            IMatrix unitStiffness = problemMatrixProvider.Matrix(elementWrappers[elementIdx]);
            penalizers[elementIdx].ScaleFactor = currentModulus;
            return unitStiffness;
        }

        public void Initialize()
        {
            elementVolumes = CalcElementVolumes();

            //TODO: does this work with embedding? What if someone else overwrites the element's property?
            for (int e = 0; e < NumElements; ++e)
            {
                elementTypes[e].DofEnumerator = penalizers[e];
            }

            problem = new ProblemStructural(model, solver);
            linearAnalyzer = new LinearAnalyzer(model, solver, problem);
            staticAnalyzer = new StaticAnalyzer(model, solver, problem, linearAnalyzer);
            
        }

        public void UpdateModelAndAnalyze(IVectorView youngModuli)
        {
            for (int e = 0; e < NumElements; ++e) penalizers[e].ScaleFactor = youngModuli[e];

            // This is necessary in order to rebuild the stiffness matrix. However it also resets each element's state, which is
            // wasteful here. 
            //TODO: If I mess with analyzer operations, I should implement a new analyzer instead of using StaticAnalyzer.
            //TODO: I could avoid this, if the analyzer decided when the matrix will be rebuilt, instead of the problem provider. 
            problem.Reset();

            if (isFirstAnalysis)
            {
                staticAnalyzer.Initialize(true);
                isFirstAnalysis = false;
            }
            else staticAnalyzer.Initialize(false);
            staticAnalyzer.Solve();
            globalDisplacements = solver.LinearSystems[0].Solution;
        }

        private Vector CalcElementVolumes()
        {
            var volumes = new double[NumElements];
            for (int e = 0; e < NumElements; ++e)
            {
                volumes[e] = elementTypes[e].Thickness * elementTypes[e].CalculateArea(modelElements[e]);
            }
            return Vector.CreateFromArray(volumes);
        }
    }
}
