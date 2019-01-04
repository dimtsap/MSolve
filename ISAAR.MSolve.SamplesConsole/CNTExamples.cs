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
using ISAAR.MSolve.Discretization.Interfaces;
using Xunit;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.FEM;
using ISAAR.MSolve.Discretization.Providers;
using System.IO;
using System.Linq;
using ISAAR.MSolve.Discretization;

namespace ISAAR.MSolve.SamplesConsole
{
    public class CNTExamples
    {
        public static void CNT_4_4_DisplacementControl()
        {
            VectorExtensions.AssignTotalAffinityCount();
            double youngModulus = 16710.0;
            double poissonRatio = 0.034;
            double nodalDisplacement = 23.73;
            double area = 5.594673861218848e-003;
            double inertiaY = 2.490804749753243e-006;
            double inertiaZ = 2.490804749753243e-006;
            double torsionalInertia = inertiaY / 2.0;
            double effectiveAreaY = area;
            double effectiveAreaZ = area;

            string workingDirectory = @"..\..\..\Resources\Beam3DInputFiles";
            string geometryFileName = "CNT-4-4-L=25-Geometry.inp";
            string connectivityFileName = "CNT-4-4-L=25-ConnMatr.inp";
            int increments = 20;

            int nNodes = File.ReadLines(workingDirectory + '\\' + geometryFileName).Count();
            int nElems = File.ReadLines(workingDirectory + '\\' + connectivityFileName).Count();
            int monitorNode_1 = nNodes - 1;
            int monitorNode_2 = nNodes;

            // Create new 3D material
            ElasticMaterial3D material_1 = new ElasticMaterial3D
            {
                YoungModulus = youngModulus,
                PoissonRatio = poissonRatio,
            };

            ElasticMaterial3D material_2 = new ElasticMaterial3D
            {
                YoungModulus = 100.0 * youngModulus,
                PoissonRatio = poissonRatio,
            };

            // Model creation
            Model model = new Model();

            // Subdomains
            model.SubdomainsDictionary.Add(1, new Subdomain() { ID = 1 });

            // Node creation
            IList<Node> nodes = new List<Node>();
            using (TextReader reader = File.OpenText(workingDirectory + '\\' + geometryFileName))
            {
                for (int i = 0; i < nNodes; i++)
                {
                    string text = reader.ReadLine();
                    string[] bits = text.Split(',');
                    int nodeID = int.Parse(bits[0]);
                    double nodeX = double.Parse(bits[1]);
                    double nodeY = double.Parse(bits[2]);
                    double nodeZ = double.Parse(bits[3]);
                    nodes.Add(new Node { ID = nodeID, X = nodeX, Y = nodeY, Z = nodeZ });
                    model.NodesDictionary.Add(nodeID, nodes[i]);
                }
            }

            // Constraints
            IList<Node> constraintsNodes = new List<Node>();
            constraintsNodes.Add(nodes[1 - 1]);
            constraintsNodes.Add(nodes[2 - 1]);
            constraintsNodes.Add(nodes[617 - 1]);
            constraintsNodes.Add(nodes[618 - 1]);
            constraintsNodes.Add(nodes[1029 - 1]);
            constraintsNodes.Add(nodes[1030 - 1]);
            constraintsNodes.Add(nodes[1441 - 1]);
            constraintsNodes.Add(nodes[1442 - 1]);
            constraintsNodes.Add(nodes[1649 - 1]);

            for (int i = 0; i < 9; i++)
            {
                int iNode = constraintsNodes[i].ID;
                model.NodesDictionary[iNode].Constraints.Add(new Constraint { DOF = DOFType.X });
                model.NodesDictionary[iNode].Constraints.Add(new Constraint { DOF = DOFType.Y });
                model.NodesDictionary[iNode].Constraints.Add(new Constraint { DOF = DOFType.Z });
                model.NodesDictionary[iNode].Constraints.Add(new Constraint { DOF = DOFType.RotX });
                model.NodesDictionary[iNode].Constraints.Add(new Constraint { DOF = DOFType.RotY });
                model.NodesDictionary[iNode].Constraints.Add(new Constraint { DOF = DOFType.RotZ });
            }

            // Applied displacement
            model.NodesDictionary[monitorNode_2].Constraints.Add(new Constraint { DOF = DOFType.Y, Amount = nodalDisplacement });

            // Create new Beam3D section and element
            var beamSection = new BeamSection3D(area, inertiaY, inertiaZ, torsionalInertia, effectiveAreaY, effectiveAreaZ);

            // Generate elements
            using (TextReader reader = File.OpenText(workingDirectory + '\\' + connectivityFileName))
            {
                for (int i = 0; i < (nElems - 16); i++)
                {
                    string text = reader.ReadLine();
                    string[] bits = text.Split(',');
                    int elementID = int.Parse(bits[0]);
                    int node1 = int.Parse(bits[1]);
                    int node2 = int.Parse(bits[2]);
                    // element nodes
                    IList<Node> elementNodes = new List<Node>();
                    elementNodes.Add(model.NodesDictionary[node1]);
                    elementNodes.Add(model.NodesDictionary[node2]);
                    // create element
                    var beam_1 = new Beam3DCorotationalQuaternion(elementNodes, material_1, 7.85, beamSection);
                    var element = new Element { ID = elementID, ElementType = beam_1 };
                    // Add nodes to the created element
                    element.AddNode(model.NodesDictionary[node1]);
                    element.AddNode(model.NodesDictionary[node2]);
                    // beam stiffness matrix
                    var a = beam_1.StiffnessMatrix(element);
                    // Add beam element to the element and subdomains dictionary of the model
                    model.ElementsDictionary.Add(element.ID, element);
                    model.SubdomainsDictionary[1].ElementsDictionary.Add(element.ID, element);
                }
                for (int i = (nElems - 16); i < nElems; i++)
                {
                    string text = reader.ReadLine();
                    string[] bits = text.Split(',');
                    int elementID = int.Parse(bits[0]);
                    int node1 = int.Parse(bits[1]);
                    int node2 = int.Parse(bits[2]);
                    // element nodes
                    IList<Node> elementNodes = new List<Node>();
                    elementNodes.Add(model.NodesDictionary[node1]);
                    elementNodes.Add(model.NodesDictionary[node2]);
                    // create element
                    var beam_2 = new Beam3DCorotationalQuaternion(elementNodes, material_2, 7.85, beamSection);
                    var element = new Element { ID = elementID, ElementType = beam_2 };
                    // Add nodes to the created element
                    element.AddNode(model.NodesDictionary[node1]);
                    element.AddNode(model.NodesDictionary[node2]);
                    // beam stiffness matrix
                    var a = beam_2.StiffnessMatrix(element);
                    // Add beam element to the element and subdomains dictionary of the model
                    model.ElementsDictionary.Add(element.ID, element);
                    model.SubdomainsDictionary[1].ElementsDictionary.Add(element.ID, element);
                }
            }

            // Needed in order to make all the required data structures
            model.ConnectDataStructures();

            // Choose linear equation system solver
            var linearSystems = new Dictionary<int, ILinearSystem>();
            linearSystems[1] = new SkylineLinearSystem(1, model.Subdomains[0].Forces);

            // Skyline Solver
            SolverSkyline solver = new SolverSkyline(linearSystems[1]);

            // Choose the provider of the problem -> here a structural problem
            ProblemStructural provider = new ProblemStructural(model, linearSystems);

            var linearSystemsArray = new[] { linearSystems[1] };
            var subdomainUpdaters = new[] { new NonLinearSubdomainUpdater(model.Subdomains[0]) };
            var subdomainMappers = new[] { new SubdomainGlobalMapping(model.Subdomains[0]) };
            var equivalentLoadsAssemblers = new[] { new EquivalentLoadsAssembler(model.Subdomains[0], new ElementStructuralStiffnessProvider()) };
            int totalDOFs = model.TotalDOFs;
            DisplacementControlNonLinearAnalyzer childAnalyzer = new DisplacementControlNonLinearAnalyzer(solver, linearSystemsArray, subdomainUpdaters, subdomainMappers,
            equivalentLoadsAssemblers, provider, increments, totalDOFs);

            // Choose parent analyzer -> Parent: Static
            StaticAnalyzer parentAnalyzer = new StaticAnalyzer(provider, childAnalyzer, linearSystems);

            parentAnalyzer.BuildMatrices();
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            int monDOF = linearSystems[1].Solution.Length - 5;  //(6 * monitorNode_2 - 4);

            double monitorDisplacement = linearSystems[1].Solution[9841];

            Assert.Equal(-21.2328445476855, monitorDisplacement, 8);
        }

        public static void CNT_4_4_NewtonRaphson()
        {
            VectorExtensions.AssignTotalAffinityCount();
            double youngModulus = 16710.0;
            double poissonRatio = 0.034;
            double nodalLoad = 1.2;
            double area = 5.594673861218848e-003;
            double inertiaY = 2.490804749753243e-006;
            double inertiaZ = 2.490804749753243e-006;
            double torsionalInertia = inertiaY / 2.0;
            double effectiveAreaY = area;
            double effectiveAreaZ = area;
            string workingDirectory = @"..\..\..\Resources\Beam3DInputFiles";
            string geometryFileName = "CNT-4-4-L=25-Geometry.inp";
            string connectivityFileName = "CNT-4-4-L=25-ConnMatr.inp";
            int increments = 100;

            //Read number of nodes and number of elements from input files
            int nNodes = File.ReadLines(workingDirectory + '\\' + geometryFileName).Count();
            int nElems = File.ReadLines(workingDirectory + '\\' + connectivityFileName).Count();
            int monitorNode_1 = nNodes - 1;
            int monitorNode_2 = nNodes;

            // Create new 3D material
            ElasticMaterial3D material_1 = new ElasticMaterial3D
            {
                YoungModulus = youngModulus,
                PoissonRatio = poissonRatio,
            };

            ElasticMaterial3D material_2 = new ElasticMaterial3D
            {
                YoungModulus = 100.0 * youngModulus,
                PoissonRatio = poissonRatio,
            };

            // Model creation
            Model model = new Model();

            // Subdomains
            model.SubdomainsDictionary.Add(1, new Subdomain() { ID = 1 });

            // Node creation
            IList<Node> nodes = new List<Node>();
            using (TextReader reader = File.OpenText(workingDirectory + '\\' + geometryFileName))
            {
                for (int i = 0; i < nNodes; i++)
                {
                    string text = reader.ReadLine();
                    string[] bits = text.Split(',');
                    int nodeID = int.Parse(bits[0]);
                    double nodeX = double.Parse(bits[1]);
                    double nodeY = double.Parse(bits[2]);
                    double nodeZ = double.Parse(bits[3]);
                    nodes.Add(new Node { ID = nodeID, X = nodeX, Y = nodeY, Z = nodeZ });
                    model.NodesDictionary.Add(nodeID, nodes[i]);
                }
            }

            // Constraints
            IList<Node> constraintsNodes = new List<Node>();
            constraintsNodes.Add(nodes[1 - 1]);
            constraintsNodes.Add(nodes[2 - 1]);
            constraintsNodes.Add(nodes[617 - 1]);
            constraintsNodes.Add(nodes[618 - 1]);
            constraintsNodes.Add(nodes[1029 - 1]);
            constraintsNodes.Add(nodes[1030 - 1]);
            constraintsNodes.Add(nodes[1441 - 1]);
            constraintsNodes.Add(nodes[1442 - 1]);
            constraintsNodes.Add(nodes[1649 - 1]);

            for (int i = 0; i < 9; i++)
            {
                int iNode = constraintsNodes[i].ID;
                model.NodesDictionary[iNode].Constraints.Add(new Constraint { DOF = DOFType.X });
                model.NodesDictionary[iNode].Constraints.Add(new Constraint { DOF = DOFType.Y });
                model.NodesDictionary[iNode].Constraints.Add(new Constraint { DOF = DOFType.Z });
                model.NodesDictionary[iNode].Constraints.Add(new Constraint { DOF = DOFType.RotX });
                model.NodesDictionary[iNode].Constraints.Add(new Constraint { DOF = DOFType.RotY });
                model.NodesDictionary[iNode].Constraints.Add(new Constraint { DOF = DOFType.RotZ });
            }

            // Create new Beam3D section and element
            var beamSection = new BeamSection3D(area, inertiaY, inertiaZ, torsionalInertia, effectiveAreaY, effectiveAreaZ);

            // Generate elements
            using (TextReader reader = File.OpenText(workingDirectory + '\\' + connectivityFileName))
            {
                for (int i = 0; i < (nElems - 16); i++)
                {
                    string text = reader.ReadLine();
                    string[] bits = text.Split(',');
                    int elementID = int.Parse(bits[0]);
                    int node1 = int.Parse(bits[1]);
                    int node2 = int.Parse(bits[2]);
                    // element nodes
                    IList<Node> elementNodes = new List<Node>();
                    elementNodes.Add(model.NodesDictionary[node1]);
                    elementNodes.Add(model.NodesDictionary[node2]);
                    // create element
                    var beam_1 = new Beam3DCorotationalQuaternion(elementNodes, material_1, 7.85, beamSection);
                    var element = new Element { ID = elementID, ElementType = beam_1 };
                    // Add nodes to the created element
                    element.AddNode(model.NodesDictionary[node1]);
                    element.AddNode(model.NodesDictionary[node2]);
                    // beam stiffness matrix
                    var a = beam_1.StiffnessMatrix(element);
                    // Add beam element to the element and subdomains dictionary of the model
                    model.ElementsDictionary.Add(element.ID, element);
                    model.SubdomainsDictionary[1].ElementsDictionary.Add(element.ID, element);
                }
                for (int i = (nElems - 16); i < nElems; i++)
                {
                    string text = reader.ReadLine();
                    string[] bits = text.Split(',');
                    int elementID = int.Parse(bits[0]);
                    int node1 = int.Parse(bits[1]);
                    int node2 = int.Parse(bits[2]);
                    // element nodes
                    IList<Node> elementNodes = new List<Node>();
                    elementNodes.Add(model.NodesDictionary[node1]);
                    elementNodes.Add(model.NodesDictionary[node2]);
                    // create element
                    var beam_2 = new Beam3DCorotationalQuaternion(elementNodes, material_2, 7.85, beamSection);
                    var element = new Element { ID = elementID, ElementType = beam_2 };
                    // Add nodes to the created element
                    element.AddNode(model.NodesDictionary[node1]);
                    element.AddNode(model.NodesDictionary[node2]);
                    // beam stiffness matrix
                    var a = beam_2.StiffnessMatrix(element);
                    // Add beam element to the element and subdomains dictionary of the model
                    model.ElementsDictionary.Add(element.ID, element);
                    model.SubdomainsDictionary[1].ElementsDictionary.Add(element.ID, element);
                }
            }

            // Add nodal load values at the top nodes of the model
            model.Loads.Add(new Load() { Amount = nodalLoad, Node = model.NodesDictionary[monitorNode_2], DOF = DOFType.Y });

            // Needed in order to make all the required data structures
            model.ConnectDataStructures();

            // Choose linear equation system solver
            var linearSystems = new Dictionary<int, ILinearSystem>();
            linearSystems[1] = new SkylineLinearSystem(1, model.Subdomains[0].Forces);

            // Skyline Solver
            SolverSkyline solver = new SolverSkyline(linearSystems[1]);

            // Choose the provider of the problem -> here a structural problem
            ProblemStructural provider = new ProblemStructural(model, linearSystems);

            // Choose child analyzer -> Child: NewtonRaphsonNonLinearAnalyzer
            var linearSystemsArray = new[] { linearSystems[1] };
            var subdomainUpdaters = new[] { new NonLinearSubdomainUpdater(model.Subdomains[0]) };
            var subdomainMappers = new[] { new SubdomainGlobalMapping(model.Subdomains[0]) };
            int totalDOFs = model.TotalDOFs;
            NewtonRaphsonNonLinearAnalyzer childAnalyzer = new NewtonRaphsonNonLinearAnalyzer(solver, linearSystemsArray, subdomainUpdaters, subdomainMappers,
            provider, increments, totalDOFs);

            // Choose parent analyzer -> Parent: Static
            StaticAnalyzer parentAnalyzer = new StaticAnalyzer(provider, childAnalyzer, linearSystems);

            parentAnalyzer.BuildMatrices();
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            int monDOF = linearSystems[1].Solution.Length - 5;  //(6 * monitorNode_2 - 4);

            double monitorDisplacement = linearSystems[1].Solution[monDOF];

            Assert.Equal(23.7373042863345, monitorDisplacement, 8);
        }
    }
}
