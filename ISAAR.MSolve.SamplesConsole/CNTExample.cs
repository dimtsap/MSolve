﻿using System;
using System.Collections.Generic;
using System.Linq;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Meshes;
using ISAAR.MSolve.FEM.Meshes.Custom;
using ISAAR.MSolve.FEM.Meshes.GMSH;
using ISAAR.MSolve.Logging.VTK;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Materials.VonMisesStress;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.Solvers.Skyline;
using System.IO;
using ISAAR.MSolve.FEM.Materials;
using ISAAR.MSolve.FEM.Elements.SupportiveClasses;

namespace ISAAR.MSolve.SamplesConsole
{
    public class CNTExample
    {
        public static void Run()
        {
            VectorExtensions.AssignTotalAffinityCount();
            double youngModulus = 16710.0;
            double poissonRatio = 0.034;
            double nodalLoad = -1.0;
            double area = 5.594673861218848e-003;
            double inertiaY = 2.490804749753243e-006;
            double inertiaZ = 2.490804749753243e-006;
            double torsionalInertia = inertiaY/2.0;
            double effectiveAreaY = area;
            double effectiveAreaZ = area;
            string workingDirectory = @"E:\GEORGE_DATA\DESKTOP\CNTExample";
            string geometryFileName = "CNT-8-8-L=100-h=3-Geometry.inp";
            string connectivityFileName = "CNT-8-8-L=100-h=3-ConnMatr.inp";

            int nNodes = File.ReadLines(workingDirectory + '\\'  + geometryFileName).Count();
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
                YoungModulus = 100.0*youngModulus,
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
            model.NodesDictionary[monitorNode_1].Constraints.Add(DOFType.X);
            model.NodesDictionary[monitorNode_1].Constraints.Add(DOFType.Y);
            model.NodesDictionary[monitorNode_1].Constraints.Add(DOFType.Z);
            model.NodesDictionary[monitorNode_2].Constraints.Add(DOFType.X);
            model.NodesDictionary[monitorNode_2].Constraints.Add(DOFType.Y);

            // Create new Beam3D section and element
            var beamSection = new BeamSection3D(area, inertiaY, inertiaZ, torsionalInertia, effectiveAreaY, effectiveAreaZ);
            
            // element nodes
            IList<Node> elementNodes = new List<Node>();

            // Generate elements
            using (TextReader reader = File.OpenText(workingDirectory + '\\' + connectivityFileName))
            {
                for (int i = 0; i < (nElems - 32); i++)
                {
                    string text = reader.ReadLine();
                    string[] bits = text.Split(',');
                    int elementID = int.Parse(bits[0]);
                    int node1 = int.Parse(bits[1]);
                    int node2 = int.Parse(bits[2]);
                    // element nodes
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
                for (int i = (nElems - 31); i < nElems; i++)
                {
                    string text = reader.ReadLine();
                    string[] bits = text.Split(',');
                    int elementID = int.Parse(bits[0]);
                    int node1 = int.Parse(bits[1]);
                    int node2 = int.Parse(bits[2]);
                    // element nodes
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
            model.Loads.Add(new Load() { Amount = nodalLoad, Node = model.NodesDictionary[monitorNode_2], DOF = DOFType.Z });

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
            //int maximumIteration = 120;
            //int iterationStepsForMatrixRebuild = 500;
            NewtonRaphsonNonLinearAnalyzer childAnalyzer = new NewtonRaphsonNonLinearAnalyzer(solver, linearSystemsArray, subdomainUpdaters, subdomainMappers,
            provider, increments, totalDOFs);

            // Choose parent analyzer -> Parent: Static
            StaticAnalyzer parentAnalyzer = new StaticAnalyzer(provider, childAnalyzer, linearSystems);

            parentAnalyzer.BuildMatrices();
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            int monDOF = (6 * monitorNode_2 - 3);
            double monitorDisplacement = linearSystems[1].Solution[monDOF - 1];
        }

    }
}
