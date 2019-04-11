﻿using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Analyzers.Interfaces;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Elements.SupportiveClasses;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Materials;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.Solvers.Skyline;
using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Interfaces;
using Xunit;
using ISAAR.MSolve.Solvers.Direct;
using ISAAR.MSolve.LinearAlgebra.Output;
using ISAAR.MSolve.Logging;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Analyzers.NonLinear;
using ISAAR.MSolve.Solvers;

namespace ISAAR.MSolve.Tests
{
    public class Beam3DQuaternionNonlinearTest
    {
        private const string output = @"C:\Users\Serafeim\Desktop\output.txt";

        [Fact]
        public void SolveNLBeam()
        {
            var m = new Model();
            m.NodesDictionary.Add(1, new Node() { ID = 1, X = 0, Y = 0, Z = 0 });
            m.NodesDictionary.Add(2, new Node() { ID = 2, X = 5, Y = 0, Z = 0 });
            m.NodesDictionary[1].Constraints.Add(new Constraint { DOF = DOFType.X });
            m.NodesDictionary[1].Constraints.Add(new Constraint { DOF = DOFType.Y });
            m.NodesDictionary[1].Constraints.Add(new Constraint { DOF = DOFType.Z });
            m.NodesDictionary[1].Constraints.Add(new Constraint { DOF = DOFType.RotX });
            m.NodesDictionary[1].Constraints.Add(new Constraint { DOF = DOFType.RotY });
            m.NodesDictionary[1].Constraints.Add(new Constraint { DOF = DOFType.RotZ });
            m.ElementsDictionary.Add(1, new Element()
            {
                ID = 1,
                ElementType = new Beam3DCorotationalQuaternion(m.Nodes, new ElasticMaterial3D() { YoungModulus = 2.1e6, PoissonRatio = 0.2 }, 1,
                new BeamSection3D(0.06, 0.0002, 0.00045, 0.000818, 0.05, 0.05))
            });
            m.ElementsDictionary[1].AddNodes(m.Nodes);
            m.SubdomainsDictionary.Add(1, new Subdomain() { ID = 1 });
            m.SubdomainsDictionary[1].ElementsDictionary.Add(1, m.ElementsDictionary[1]);
            m.Loads.Add(new Load() { Node = m.NodesDictionary[2], Amount = 100, DOF = DOFType.Y });

            m.ConnectDataStructures();
            VectorExtensions.AssignTotalAffinityCount();

            var linearSystems = new Dictionary<int, ILinearSystem>();
            linearSystems[1] = new SkylineLinearSystem(1, m.Subdomains[0].Forces);
            var linearSystemsArray = new[] { linearSystems[1] };
            var subdomainUpdaters = new[] { new NonLinearSubdomainUpdater(m.Subdomains[0]) };
            var subdomainMappers = new[] { new SubdomainGlobalMapping(m.Subdomains[0]) };
            int increments = 10;
            int totalDOFs = m.TotalDOFs;
            SolverSkyline solver = new SolverSkyline(linearSystems[1]);
            ProblemStructural provider = new ProblemStructural(m, linearSystems);
            NewtonRaphsonNonLinearAnalyzer childAnalyzer = new NewtonRaphsonNonLinearAnalyzer(solver, linearSystemsArray, subdomainUpdaters, subdomainMappers,
                provider, increments, totalDOFs);

            // Choose parent analyzer -> Parent: Static
            StaticAnalyzer parentAnalyzer = new StaticAnalyzer(provider, childAnalyzer, linearSystems);
            parentAnalyzer.BuildMatrices();
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();
        }

        public void SolveNLBeam_v2()
        {
            var m = new Model_v2();
            m.NodesDictionary.Add(1, new Node_v2() { ID = 1, X = 0, Y = 0, Z = 0 });
            m.NodesDictionary.Add(2, new Node_v2() { ID = 2, X = 5, Y = 0, Z = 0 });
            m.NodesDictionary[1].Constraints.Add(new Constraint { DOF = DOFType.X });
            m.NodesDictionary[1].Constraints.Add(new Constraint { DOF = DOFType.Y });
            m.NodesDictionary[1].Constraints.Add(new Constraint { DOF = DOFType.Z });
            m.NodesDictionary[1].Constraints.Add(new Constraint { DOF = DOFType.RotX });
            m.NodesDictionary[1].Constraints.Add(new Constraint { DOF = DOFType.RotY });
            m.NodesDictionary[1].Constraints.Add(new Constraint { DOF = DOFType.RotZ });
            m.ElementsDictionary.Add(1, new Element_v2()
            {
                ID = 1,
                ElementType = new Beam3DCorotationalQuaternion_v2(m.Nodes, new ElasticMaterial3D_v2() { YoungModulus = 2.1e6, PoissonRatio = 0.2 }, 1,
                new BeamSection3D(0.06, 0.0002, 0.00045, 0.000818, 0.05, 0.05))
            });
            m.ElementsDictionary[1].AddNodes(m.Nodes);
            m.SubdomainsDictionary.Add(1, new Subdomain_v2(1));
            m.SubdomainsDictionary[1].Elements.Add(m.ElementsDictionary[1]);
            m.Loads.Add(new Load_v2() { Node = m.NodesDictionary[2], Amount = 100, DOF = DOFType.Y });

            // Solver
            var solverBuilder = new SkylineSolver.Builder();
            ISolver_v2 solver = solverBuilder.BuildSolver(m);

            // Problem type
            var provider = new ProblemStructural_v2(m, solver);

            // Analyzers
            int increments = 10;
            var childAnalyzerBuilder = new LoadControlAnalyzer_v2.Builder(m, solver, provider, increments);
            //childAnalyzerBuilder.SubdomainUpdaters = new[] { new NonLinearSubdomainUpdater_v2(model.SubdomainsDictionary[subdomainID]) }; // This is the default
            LoadControlAnalyzer_v2 childAnalyzer = childAnalyzerBuilder.Build();
            var parentAnalyzer = new StaticAnalyzer_v2(m, solver, provider, childAnalyzer);

            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();
        }

        [Fact]
        public void CantileverYBeam3DQuaternionNonlinearTest()
        {
            VectorExtensions.AssignTotalAffinityCount();
            double youngModulus = 21000.0;
            double poissonRatio = 0.3;
            double nodalLoad = 20000.0;
            double area = 91.04;
            double inertiaY = 2843.0;
            double inertiaZ = 8091.0;
            double torsionalInertia = 76.57;
            double effectiveAreaY = 91.04;
            double effectiveAreaZ = 91.04;
            int nNodes = 3;
            int nElems = 2;
            int monitorNode = 3;

            // Create new 3D material
            ElasticMaterial3D material = new ElasticMaterial3D
            {
                YoungModulus = youngModulus,
                PoissonRatio = poissonRatio,
            };

            // Node creation
            IList<Node> nodes = new List<Node>();
            Node node1 = new Node { ID = 1, X = 0.0,   Y = 0.0, Z = 0.0 };
            Node node2 = new Node { ID = 2, X = 100.0, Y = 0.0, Z = 0.0 };
            Node node3 = new Node { ID = 3, X = 200.0, Y = 0.0, Z = 0.0 };
            
            nodes.Add(node1);
            nodes.Add(node2);
            nodes.Add(node3);
            
            // Model creation
            Model model = new Model();

            // Add a single subdomain to the model
            model.SubdomainsDictionary.Add(1, new Subdomain() { ID = 1 });

            // Add nodes to the nodes dictonary of the model
            for (int i = 0; i < nodes.Count; ++i)
            {
                model.NodesDictionary.Add(i + 1, nodes[i]);
            }

            // Constrain bottom nodes of the model
            model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = DOFType.X });
            model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = DOFType.Y });
            model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = DOFType.Z });
            model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = DOFType.RotX });
            model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = DOFType.RotY });
            model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = DOFType.RotZ });  

            // Generate elements of the structure
            int iNode = 1;
            for (int iElem = 0; iElem < nElems; iElem++)
            {
                // element nodes
                IList<Node> elementNodes = new List<Node>();
                elementNodes.Add(model.NodesDictionary[iNode]);
                elementNodes.Add(model.NodesDictionary[iNode + 1]);

                // Create new Beam3D section and element
                var beamSection = new BeamSection3D(area, inertiaY, inertiaZ, torsionalInertia, effectiveAreaY, effectiveAreaZ);
                var beam = new Beam3DCorotationalQuaternion(elementNodes, material, 7.85, beamSection);

                // Create elements
                var element = new Element()
                {
                    ID = iElem + 1,
                    ElementType = beam
                };

                // Add nodes to the created element
                element.AddNode(model.NodesDictionary[iNode]);
                element.AddNode(model.NodesDictionary[iNode + 1]);

                var a = beam.StiffnessMatrix(element);
                //a.PrintMatrix2D(output);

                // Add beam element to the element and subdomains dictionary of the model
                model.ElementsDictionary.Add(element.ID, element);
                model.SubdomainsDictionary[1].ElementsDictionary.Add(element.ID, element);
                iNode++;
            }    

            // Add nodal load values at the top nodes of the model
            model.Loads.Add(new Load() { Amount = nodalLoad, Node = model.NodesDictionary[monitorNode], DOF = DOFType.Y });

            // Needed in order to make all the required data structures
            model.ConnectDataStructures();

            // Choose linear equation system solver
            var linearSystems = new Dictionary<int, ILinearSystem>();
            linearSystems[1] = new SkylineLinearSystem(1, model.Subdomains[0].Forces);
            SolverSkyline solver = new SolverSkyline(linearSystems[1]);

            // Choose the provider of the problem -> here a structural problem
            ProblemStructural provider = new ProblemStructural(model, linearSystems);

            // Choose child analyzer -> Child: NewtonRaphsonNonLinearAnalyzer
            var linearSystemsArray = new[] { linearSystems[1] };
            var subdomainUpdaters = new[] { new NonLinearSubdomainUpdater(model.Subdomains[0]) };
            var subdomainMappers = new[] { new SubdomainGlobalMapping(model.Subdomains[0]) };
            int increments = 10;
            int totalDOFs = model.TotalDOFs;
            NewtonRaphsonNonLinearAnalyzer childAnalyzer = new NewtonRaphsonNonLinearAnalyzer(solver, linearSystemsArray, subdomainUpdaters, subdomainMappers,
            provider, increments, totalDOFs);

            // Choose parent analyzer -> Parent: Static
            StaticAnalyzer parentAnalyzer = new StaticAnalyzer(provider, childAnalyzer, linearSystems);

            parentAnalyzer.BuildMatrices();
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            Assert.Equal(148.936792350562, linearSystems[1].Solution[7], 2);
        }

        [Fact]
        public void CantileverYBeam3DQuaternionNonlinearTest_v2()
        {
            double youngModulus = 21000.0;
            double poissonRatio = 0.3;
            double nodalLoad = 20000.0;
            double area = 91.04;
            double inertiaY = 2843.0;
            double inertiaZ = 8091.0;
            double torsionalInertia = 76.57;
            double effectiveAreaY = 91.04;
            double effectiveAreaZ = 91.04;
            int nNodes = 3;
            int nElems = 2;
            int monitorNode = 3;

            // Create new 3D material
            var material = new ElasticMaterial3D_v2
            {
                YoungModulus = youngModulus,
                PoissonRatio = poissonRatio,
            };

            // Node creation
            IList<Node_v2> nodes = new List<Node_v2>();
            Node_v2 node1 = new Node_v2 { ID = 1, X = 0.0, Y = 0.0, Z = 0.0 };
            Node_v2 node2 = new Node_v2 { ID = 2, X = 100.0, Y = 0.0, Z = 0.0 };
            Node_v2 node3 = new Node_v2 { ID = 3, X = 200.0, Y = 0.0, Z = 0.0 };

            nodes.Add(node1);
            nodes.Add(node2);
            nodes.Add(node3);

            // Model creation
            Model_v2 model = new Model_v2();

            // Add a single subdomain to the model
            model.SubdomainsDictionary.Add(1, new Subdomain_v2(1));

            // Add nodes to the nodes dictonary of the model
            for (int i = 0; i < nodes.Count; ++i)
            {
                model.NodesDictionary.Add(i + 1, nodes[i]);
            }

            // Constrain bottom nodes of the model
            model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = DOFType.X });
            model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = DOFType.Y });
            model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = DOFType.Z });
            model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = DOFType.RotX });
            model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = DOFType.RotY });
            model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = DOFType.RotZ });

            // Generate elements of the structure
            int iNode = 1;
            for (int iElem = 0; iElem < nElems; iElem++)
            {
                // element nodes
                IList<Node_v2> elementNodes = new List<Node_v2>();
                elementNodes.Add(model.NodesDictionary[iNode]);
                elementNodes.Add(model.NodesDictionary[iNode + 1]);

                // Create new Beam3D section and element
                var beamSection = new BeamSection3D(area, inertiaY, inertiaZ, torsionalInertia, effectiveAreaY, effectiveAreaZ);
                var beam = new Beam3DCorotationalQuaternion_v2(elementNodes, material, 7.85, beamSection);

                // Create elements
                var element = new Element_v2()
                {
                    ID = iElem + 1,
                    ElementType = beam
                };

                // Add nodes to the created element
                element.AddNode(model.NodesDictionary[iNode]);
                element.AddNode(model.NodesDictionary[iNode + 1]);

                var a = beam.StiffnessMatrix(element);
                //var writer = new FullMatrixWriter();
                //writer.WriteToFile(a, output);


                // Add beam element to the element and subdomains dictionary of the model
                model.ElementsDictionary.Add(element.ID, element);
                model.SubdomainsDictionary[1].Elements.Add(element);
                iNode++;
            }

            // Add nodal load values at the top nodes of the model
            model.Loads.Add(new Load_v2() { Amount = nodalLoad, Node = model.NodesDictionary[monitorNode], DOF = DOFType.Y });

            // Solver
            var solverBuilder = new SkylineSolver.Builder();
            ISolver_v2 solver = solverBuilder.BuildSolver(model);

            // Problem type
            var provider = new ProblemStructural_v2(model, solver);

            // Analyzers
            int increments = 10;
            var childAnalyzerBuilder = new LoadControlAnalyzer_v2.Builder(model, solver, provider, increments);
            childAnalyzerBuilder.ResidualTolerance = 1E-3;
            //childAnalyzerBuilder.SubdomainUpdaters = new[] { new NonLinearSubdomainUpdater_v2(model.SubdomainsDictionary[subdomainID]) }; // This is the default
            LoadControlAnalyzer_v2 childAnalyzer = childAnalyzerBuilder.Build();
            var parentAnalyzer = new StaticAnalyzer_v2(model, solver, provider, childAnalyzer);

            // Request output
            childAnalyzer.LogFactories[1] = new LinearAnalyzerLogFactory_v2(new int[] { 7 });

            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            // Check output
            DOFSLog_v2 log = (DOFSLog_v2)childAnalyzer.Logs[1][0]; //There is a list of logs for each subdomain and we want the first one
            Assert.Equal(148.936792350562, log.DOFValues[7], 2);
        }

        [Fact]
        public void PlaneFrameTest()
        {
            VectorExtensions.AssignTotalAffinityCount();
            double youngModulus = 21000.0;
            double poissonRatio = 0.3;
            double nodalLoad = 500000.0;
            double area = 91.04;
            double inertiaY = 2843.0;
            double inertiaZ = 8091.0;
            double torsionalInertia = 76.57;
            double effectiveAreaY = 91.04;
            double effectiveAreaZ = 91.04;
            int nNodes = 4;
            int nElems = 3;
            int monitorNode = 2;

            // Create new 3D material
            ElasticMaterial3D material = new ElasticMaterial3D
            {
                YoungModulus = youngModulus,
                PoissonRatio = poissonRatio,
            };

            // Node creation
            IList<Node> nodes = new List<Node>();
            Node node1 = new Node { ID = 1, X = 0.0, Y = 0.0, Z = 0.0 };
            Node node2 = new Node { ID = 2, X = 0.0, Y = 100.0, Z = 0.0 };
            Node node3 = new Node { ID = 3, X = 100.0, Y = 100.0, Z = 0.0 };
            Node node4 = new Node { ID = 4, X = 100.0, Y = 0.0, Z = 0.0 };

            nodes.Add(node1);
            nodes.Add(node2);
            nodes.Add(node3);
            nodes.Add(node4);

            // Model creation
            Model model = new Model();

            // Add a single subdomain to the model
            model.SubdomainsDictionary.Add(1, new Subdomain() { ID = 1 });

            // Add nodes to the nodes dictonary of the model
            for (int i = 0; i < nodes.Count; ++i)
            {
                model.NodesDictionary.Add(i + 1, nodes[i]);
            }

            // Constrain first and last nodes of the model
            model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = DOFType.X });
            model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = DOFType.Y });
            model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = DOFType.Z });
            model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = DOFType.RotX });
            model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = DOFType.RotY });
            model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = DOFType.RotZ });
            model.NodesDictionary[4].Constraints.Add(new Constraint { DOF = DOFType.X });
            model.NodesDictionary[4].Constraints.Add(new Constraint { DOF = DOFType.Y });
            model.NodesDictionary[4].Constraints.Add(new Constraint { DOF = DOFType.Z });
            model.NodesDictionary[4].Constraints.Add(new Constraint { DOF = DOFType.RotX });
            model.NodesDictionary[4].Constraints.Add(new Constraint { DOF = DOFType.RotY });
            model.NodesDictionary[4].Constraints.Add(new Constraint { DOF = DOFType.RotZ });

            // Generate elements of the structure
            int iNode = 1;
            for (int iElem = 0; iElem < nElems; iElem++)
            {
                // element nodes
                IList<Node> elementNodes = new List<Node>();
                elementNodes.Add(model.NodesDictionary[iNode]);
                elementNodes.Add(model.NodesDictionary[iNode + 1]);

                // Create new Beam3D section and element
                var beamSection = new BeamSection3D(area, inertiaY, inertiaZ, torsionalInertia, effectiveAreaY, effectiveAreaZ);
                var beam = new Beam3DCorotationalQuaternion(elementNodes, material, 7.85, beamSection);

                // Create elements
                var element = new Element()
                {
                    ID = iElem + 1,
                    ElementType = beam
                };                

                var a = beam.StiffnessMatrix(element);

                // Add nodes to the created element
                element.AddNode(model.NodesDictionary[iNode]);
                element.AddNode(model.NodesDictionary[iNode + 1]);

                // Add beam element to the element and subdomains dictionary of the model
                model.ElementsDictionary.Add(element.ID, element);
                model.SubdomainsDictionary[1].ElementsDictionary.Add(element.ID, element);
                iNode++;
            }

            // Add nodal load values at the top nodes of the model
            model.Loads.Add(new Load() { Amount = nodalLoad, Node = model.NodesDictionary[monitorNode], DOF = DOFType.X });

            // Needed in order to make all the required data structures
            model.ConnectDataStructures();

            // Choose linear equation system solver
            var linearSystems = new Dictionary<int, ILinearSystem>();
            linearSystems[1] = new SkylineLinearSystem(1, model.Subdomains[0].Forces);
            SolverSkyline solver = new SolverSkyline(linearSystems[1]);

            // Choose the provider of the problem -> here a structural problem
            ProblemStructural provider = new ProblemStructural(model, linearSystems);

            // Choose child analyzer -> Child: NewtonRaphsonNonLinearAnalyzer
            var linearSystemsArray = new[] { linearSystems[1] };
            var subdomainUpdaters = new[] { new NonLinearSubdomainUpdater(model.Subdomains[0]) };
            var subdomainMappers = new[] { new SubdomainGlobalMapping(model.Subdomains[0]) };
            int increments = 10;
            int totalDOFs = model.TotalDOFs;
            NewtonRaphsonNonLinearAnalyzer childAnalyzer = new NewtonRaphsonNonLinearAnalyzer(solver, linearSystemsArray, subdomainUpdaters, subdomainMappers,
            provider, increments, totalDOFs);

            // Choose parent analyzer -> Parent: Static
            StaticAnalyzer parentAnalyzer = new StaticAnalyzer(provider, childAnalyzer, linearSystems);

            parentAnalyzer.BuildMatrices();
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            Assert.Equal(120.1108698752, linearSystems[1].Solution[0], 2);
        }

        [Fact]
        public void PlaneFrameTest_v2()
        {
            double youngModulus = 21000.0;
            double poissonRatio = 0.3;
            double nodalLoad = 500000.0;
            double area = 91.04;
            double inertiaY = 2843.0;
            double inertiaZ = 8091.0;
            double torsionalInertia = 76.57;
            double effectiveAreaY = 91.04;
            double effectiveAreaZ = 91.04;
            int nNodes = 4;
            int nElems = 3;
            int monitorNode = 2;

            // Create new 3D material
            var material = new ElasticMaterial3D_v2
            {
                YoungModulus = youngModulus,
                PoissonRatio = poissonRatio,
            };

            // Node creation
            IList<Node_v2> nodes = new List<Node_v2>();
            Node_v2 node1 = new Node_v2 { ID = 1, X = 0.0, Y = 0.0, Z = 0.0 };
            Node_v2 node2 = new Node_v2 { ID = 2, X = 0.0, Y = 100.0, Z = 0.0 };
            Node_v2 node3 = new Node_v2 { ID = 3, X = 100.0, Y = 100.0, Z = 0.0 };
            Node_v2 node4 = new Node_v2 { ID = 4, X = 100.0, Y = 0.0, Z = 0.0 };

            nodes.Add(node1);
            nodes.Add(node2);
            nodes.Add(node3);
            nodes.Add(node4);

            // Model creation
            Model_v2 model = new Model_v2();

            // Add a single subdomain to the model
            model.SubdomainsDictionary.Add(1, new Subdomain_v2(1));

            // Add nodes to the nodes dictonary of the model
            for (int i = 0; i < nodes.Count; ++i)
            {
                model.NodesDictionary.Add(i + 1, nodes[i]);
            }

            // Constrain first and last nodes of the model
            model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = DOFType.X });
            model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = DOFType.Y });
            model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = DOFType.Z });
            model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = DOFType.RotX });
            model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = DOFType.RotY });
            model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = DOFType.RotZ });
            model.NodesDictionary[4].Constraints.Add(new Constraint { DOF = DOFType.X });
            model.NodesDictionary[4].Constraints.Add(new Constraint { DOF = DOFType.Y });
            model.NodesDictionary[4].Constraints.Add(new Constraint { DOF = DOFType.Z });
            model.NodesDictionary[4].Constraints.Add(new Constraint { DOF = DOFType.RotX });
            model.NodesDictionary[4].Constraints.Add(new Constraint { DOF = DOFType.RotY });
            model.NodesDictionary[4].Constraints.Add(new Constraint { DOF = DOFType.RotZ });

            // Generate elements of the structure
            int iNode = 1;
            for (int iElem = 0; iElem < nElems; iElem++)
            {
                // element nodes
                IList<Node_v2> elementNodes = new List<Node_v2>();
                elementNodes.Add(model.NodesDictionary[iNode]);
                elementNodes.Add(model.NodesDictionary[iNode + 1]);

                // Create new Beam3D section and element
                var beamSection = new BeamSection3D(area, inertiaY, inertiaZ, torsionalInertia, effectiveAreaY, effectiveAreaZ);
                var beam = new Beam3DCorotationalQuaternion_v2(elementNodes, material, 7.85, beamSection);

                // Create elements
                var element = new Element_v2()
                {
                    ID = iElem + 1,
                    ElementType = beam
                };

                var a = beam.StiffnessMatrix(element);

                // Add nodes to the created element
                element.AddNode(model.NodesDictionary[iNode]);
                element.AddNode(model.NodesDictionary[iNode + 1]);

                // Add beam element to the element and subdomains dictionary of the model
                model.ElementsDictionary.Add(element.ID, element);
                model.SubdomainsDictionary[1].Elements.Add(element);
                iNode++;
            }

            // Add nodal load values at the top nodes of the model
            model.Loads.Add(new Load_v2() { Amount = nodalLoad, Node = model.NodesDictionary[monitorNode], DOF = DOFType.X });

            // Solver
            var solverBuilder = new SkylineSolver.Builder();
            ISolver_v2 solver = solverBuilder.BuildSolver(model);

            // Problem type
            var provider = new ProblemStructural_v2(model, solver);

            // Analyzers
            int increments = 10;
            var childAnalyzerBuilder = new LoadControlAnalyzer_v2.Builder(model, solver, provider, increments);
            childAnalyzerBuilder.ResidualTolerance = 1E-3;
            //childAnalyzerBuilder.SubdomainUpdaters = new[] { new NonLinearSubdomainUpdater_v2(model.SubdomainsDictionary[subdomainID]) }; // This is the default
            LoadControlAnalyzer_v2 childAnalyzer = childAnalyzerBuilder.Build();
            var parentAnalyzer = new StaticAnalyzer_v2(model, solver, provider, childAnalyzer);

            // Request output
            childAnalyzer.LogFactories[1] = new LinearAnalyzerLogFactory_v2(new int[] { 0 });

            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();


            // Check output
            DOFSLog_v2 log = (DOFSLog_v2)childAnalyzer.Logs[1][0]; //There is a list of logs for each subdomain and we want the first one
            Assert.Equal(120.1108698752, log.DOFValues[0], 2);
        }
    }
}
