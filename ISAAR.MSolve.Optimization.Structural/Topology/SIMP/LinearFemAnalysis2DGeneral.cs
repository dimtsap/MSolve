﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers;

namespace ISAAR.MSolve.Optimization.Structural.Topology.SIMP
{
    public class LinearFemAnalysis2DGeneral : ILinearFemAnalysis
    {
        private readonly Model_v2 model;
        private readonly Subdomain_v2 subdomain;
        private readonly ISolver_v2 solver;
        private readonly ContinuumElement2D[] elementTypes;
        private readonly IElement_v2[] elementWrappers;
        private readonly PenaltyDofEnumerator[] penalizers;
        private ProblemStructural_v2 problem;
        private LinearAnalyzer_v2 linearAnalyzer;
        private StaticAnalyzer_v2 staticAnalyzer;
        private Vector elementVolumes;
        private IVectorView globalDisplacements;

        public LinearFemAnalysis2DGeneral(Model_v2 model, ISolver_v2 solver)
        {
            if (model.SubdomainsDictionary.Count != 1) throw new NotImplementedException(
                "Cannot perform topology optimization with multiple subdomains yet.");
            this.model = model;
            subdomain = model.SubdomainsDictionary.First().Value;
            this.solver = solver;

            // Elements
            IList<Element_v2> modelElements = model.Elements;
            NumElements = modelElements.Count;
            elementTypes = new ContinuumElement2D[NumElements];
            elementWrappers = new IElement_v2[NumElements];
            penalizers = new PenaltyDofEnumerator[NumElements];
            for (int e = 0; e < NumElements; ++e)
            {
                if (modelElements[e].ElementType is ContinuumElement2D continuum) elementTypes[e] = continuum;
                else throw new ArgumentException("2D topology optimization only works with 2D continuum elements,"
                    + $" but the element with ID = {elementTypes[e].ID} was not.");
                elementWrappers[e] = modelElements[e];
                penalizers[e] = new PenaltyDofEnumerator();
            }
            NumLoadCases = 1; //TODO: implement multiple load cases in Model, analyzers and solvers
        }

        public int NumElements { get; }

        public int NumLoadCases { get; }

        public void AnalyzeModelWithDensities(IVectorView densities)
        {
            for (int e = 0; e < NumElements; ++e) penalizers[e].Penalty = densities[e];

            // This is necessary in order to rebuild the stiffness matrix. However it also resets each element's state, which is
            // wasteful here. 
            //TODO: If I mess with analyzer operations, I should implement a new analyzer instead of using StaticAnalyzer.
            //TODO: I could avoid this, if the analyzer decided when the matrix will be rebuilt, instead of the problem provider. 
            problem.Reset();

            staticAnalyzer.Solve();
            globalDisplacements = solver.LinearSystems[0].Solution;
        }

        public double CalculateTotalVolume(IVectorView densities) => elementVolumes.DotProduct(densities);

        public IMatrixView GetPenalizedStiffnessOfElement(int elementIdx, IVectorView densities)
            => elementTypes[elementIdx].BuildStiffnessMatrix();

        public Vector GetElementDisplacements(int elementIdx, int loadCaseIdx)
            => Vector.CreateFromArray(subdomain.DofOrdering.ExtractVectorElementFromSubdomain(
                elementWrappers[elementIdx], globalDisplacements));

        public void Initialize()
        {
            elementVolumes = CalcElementVolumes();

            //TODO: does this work with embedding? What if someone else overwrites the element's property?
            for (int e = 0; e < NumElements; ++e)
            {
                elementTypes[e].DofEnumerator = penalizers[e];
            }

            problem = new ProblemStructural_v2(model, solver);
            linearAnalyzer = new LinearAnalyzer_v2(solver);
            staticAnalyzer = new StaticAnalyzer_v2(model, solver, problem, linearAnalyzer);
            staticAnalyzer.Initialize(true);
        }

        private Vector CalcElementVolumes()
        {
            var volumes = new double[NumElements];
            for (int e = 0; e < NumElements; ++e)
            {
                volumes[e] = elementTypes[e].Thickness * elementTypes[e].CalculateArea();
            }
            return Vector.CreateFromArray(volumes);
        }
    }
}
