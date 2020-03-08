using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.IGA.Entities;
using ISAAR.MSolve.IGA.Readers;
using ISAAR.MSolve.LinearAlgebra.Input;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Optimization.Algorithms.GradientBased.OC;
using ISAAR.MSolve.Optimization.Algorithms.GradientBased.OC.Bisection;
using ISAAR.MSolve.Optimization.Algorithms.GradientBased.OC.Convergence;
using ISAAR.MSolve.Optimization.Logging;
using ISAAR.MSolve.Optimization.Structural.Topology.SIMP;
using ISAAR.MSolve.Optimization.Structural.Topology.SIMP.Analysis;
using ISAAR.MSolve.Optimization.Structural.Topology.SIMP.Filtering;
using ISAAR.MSolve.Optimization.Structural.Topology.SIMP.MaterialInterpolation;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers;
using ISAAR.MSolve.Solvers.Direct;
using ISAAR.MSolve.Solvers.Ordering;
using ISAAR.MSolve.Solvers.Ordering.Reordering;
using Xunit;

namespace ISAAR.MSolve.Optimization.Structural.Tests.Topology.SIMP
{
    public class TopologySimpIGATests
    {

        [Fact]
        public void TestCantileverBeam()
        {
            string filename = Path.Combine(Directory.GetCurrentDirectory(),"InputFiles","Cantilever2D.txt");
            var youngModulus = 1000.0;
            var material=new ElasticMaterial2D(StressState2D.PlaneStress)
            {
                YoungModulus = youngModulus,
                PoissonRatio = 0.3
            };
            var modelReader = new IsogeometricReader(filename,material2D:material);
            var model=modelReader.CreateModelFromFile();

            // Forces and Boundary Conditions
            foreach (ControlPoint controlPoint in model.PatchesDictionary[0].EdgesDictionary[1].ControlPointsDictionary
                .Values)
                model.Loads.Add(new Load()
                {
                    Amount = -1,
                    Node = model.ControlPoints.ToList()[controlPoint.ID],
                    DOF = StructuralDof.TranslationY
                });

            // Boundary Conditions - Dirichlet
            foreach (ControlPoint controlPoint in model.PatchesDictionary[0].EdgesDictionary[0].ControlPointsDictionary
                .Values)
            {
                model.ControlPointsDictionary[controlPoint.ID].Constraints.Add(new Constraint() {DOF = StructuralDof.TranslationX});
                model.ControlPointsDictionary[controlPoint.ID].Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationY });
            }

            var solverBuilder = new SuiteSparseSolver.Builder();
            var solver = solverBuilder.BuildSolver(model);

            double filterAreaRadius = 1.2;
            var iga = new LinearIgaAnalysis2DGeneral(model, solver);
            var filter = new ProximityDensityFilterIGA2D(model, filterAreaRadius);

            double volumeFraction = 0.4;
            var penalty = 3.0;

            // Define the optimization
            var materialInterpolation = new PowerLawMaterialInterpolation(youngModulus, penalty, 1E-3);

            var optimAlgorithmBuilder = new OptimalityCriteriaBuilder()
            {
                DampingCoeff = 0.5,
                MoveLimit = 0.2,
                InitialBisectionLimitLower = 0.0,
                InitialBisectionLimitUpper = 1E5,
                BisectionConvergence = new BisectionConvergenceAbsoluteChange(1E-4),
                OptimalityCriteriaConvergence = new DesignVariableChangeConvergence(1E-2)
            };

            var simp = new TopologySimpLinear2D(iga, optimAlgorithmBuilder, filter, materialInterpolation, volumeFraction);
            var logger = new ObjectiveFunctionLogger();
            simp.Logger = logger;

            // Run the optimization
            simp.Initialize();
            (double compliance, Vector densities) = simp.Optimize();

        }

    }
}
