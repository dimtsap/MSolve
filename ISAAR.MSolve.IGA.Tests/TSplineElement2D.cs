using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.IGA.Entities;
using ISAAR.MSolve.IGA.Postprocessing;
using ISAAR.MSolve.IGA.Readers;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers.Direct;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.Solvers.Ordering;
using ISAAR.MSolve.Solvers.Ordering.Reordering;
using Xunit;
using VectorExtensions = ISAAR.MSolve.Numerical.LinearAlgebra.VectorExtensions;

namespace ISAAR.MSolve.IGA.Tests
{
	public class TSplineElement2D
	{
		[Fact]
		public void SquareTSpline2D()
		{
			VectorExtensions.AssignTotalAffinityCount();
			Model model = new Model();
			var filename = "square_structured";
			string filepath = $"..\\..\\..\\InputFiles\\{filename}.iga";
			IGAFileReader modelReader= new IGAFileReader(model, filepath);
			modelReader.CreateTSplineShellsModelFromFile();

			model.PatchesDictionary[0].Material = new ElasticMaterial2D(StressState2D.PlaneStrain)
			{
				PoissonRatio = 0.3,
				YoungModulus = 250
			};
			model.PatchesDictionary[0].Thickness = 1;

			foreach (var controlPoint in model.ControlPointsDictionary.Values.Where(cp => cp.X < 1e-8))
				model.ControlPointsDictionary[controlPoint.ID].Constrains.Add(new Constraint() { DOF = DOFType.X });

			foreach (var controlPoint in model.ControlPointsDictionary.Values.Where(cp => cp.Y < 1e-8&&cp.X<1-8))
				model.ControlPointsDictionary[controlPoint.ID].Constrains.Add(new Constraint() { DOF = DOFType.Y });

			foreach (var controlPoint in model.ControlPointsDictionary.Values.Where(cp => Math.Abs(cp.X-1) < 1e-8))
				model.Loads.Add(new Load()
				{
					Amount = 10,
					ControlPoint = model.ControlPointsDictionary[controlPoint.ID],
					DOF = DOFType.X
				});

			var solverBuilder = new SuiteSparseSolver.Builder();
			solverBuilder.DofOrderer = new DofOrderer(
				new NodeMajorDofOrderingStrategy(), new NullReordering());
			ISolver_v2 solver = solverBuilder.BuildSolver(model);

			// Structural problem provider
			var provider = new ProblemStructural_v2(model, solver);

			// Linear static analysis
			var childAnalyzer = new LinearAnalyzer_v2(solver);
			var parentAnalyzer = new StaticAnalyzer_v2(model, solver, provider, childAnalyzer);

			// Run the analysis
			parentAnalyzer.Initialize();
			parentAnalyzer.Solve();


			var paraview = new ParaviewTsplines2D(model, solver.LinearSystems[0].Solution, filename);
			paraview.CreateParaviewFile();
		}
	}
}
