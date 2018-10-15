﻿using System;
using System.Collections.Generic;
using System.IO;
using System.Text;
using Ansys.ACT.Interfaces.Common;
using Ansys.ACT.Interfaces.Mechanical;
using Ansys.ACT.Interfaces.Mesh;
using Ansys.ACT.Interfaces.UserObject;
using ISAAR.MSolve.FEM.Entities;
using System.Linq;
using Ansys.ACT.Automation.Mechanical;
using Ansys.ACT.Core.Callbacks;
using Ansys.ACT.Interfaces.Geometry;
using Ansys.ACT.Interfaces.Graphics.Entities;
using Ansys.ACT.Interfaces.Project;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.Geometry.Shapes;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.Solvers.Skyline;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Discretization.Interfaces;
using Newtonsoft.Json;
using SimpleExpressionEvaluator;
using INode = Ansys.ACT.Interfaces.Mesh.INode;
using Model = ISAAR.MSolve.FEM.Entities.Model;

namespace ISAAR.MSolve.ANSYS
{
	public class MainMSolve
	{
		private Dictionary<ElementTypeEnum, Dictionary<int, int[]>> ElementTypesNodesDictionary =
			new Dictionary<ElementTypeEnum, Dictionary<int, int[]>>()
			{
				{
					ElementTypeEnum.kTri3,
					new Dictionary<int, int[]>
					{
						{0, new int[] {2, 1}},
						{1, new int[] {0, 2}},
						{2, new int[] {1, 0}}
					}
				},
				{
					ElementTypeEnum.kTri6,
					new Dictionary<int, int[]>
					{
						{0, new int[] {2, 1}},
						{1, new int[] {0, 2}},
						{2, new int[] {1, 0}},
						{3, new int[] {0, 1}},
						{4, new int[] {1, 2}},
						{5, new int[] {2, 3}}
					}
				},
				{
					ElementTypeEnum.kQuad4,
					new Dictionary<int, int[]>
					{
						{0, new int[] {3, 1}},
						{1, new int[] {0, 2}},
						{2, new int[] {1, 3}},
						{3, new int[] {2, 0}},
					}
				},
				{
					ElementTypeEnum.kQuad8,
					new Dictionary<int, int[]>
					{
						{0, new int[] {3, 1}},
						{1, new int[] {0, 2}},
						{2, new int[] {1, 3}},
						{3, new int[] {2, 0}},
						{4, new int[] {0, 1}},
						{5, new int[] {1, 2}},
						{6, new int[] {2, 3}},
						{7, new int[] {3, 0}},
					}
				},
				{
					ElementTypeEnum.kTet4,
					new Dictionary<int, int[]>
					{
						{0, new int[] {2, 1, 3}},
						{1, new int[] {0, 2, 3}},
						{2, new int[] {1, 0, 3}},
						{3, new int[] {0, 1, 2}}
					}
				},
				{
					ElementTypeEnum.kTet10,
					new Dictionary<int, int[]>
					{
						{0, new int[] {2, 1, 3}},
						{1, new int[] {0, 2, 3}},
						{2, new int[] {1, 0, 3}},
						{3, new int[] {0, 1, 2}},
						{4, new int[] {0, 1}},
						{5, new int[] {1, 2}},
						{6, new int[] {2, 0}},
						{7, new int[] {0, 3}},
						{8, new int[] {1, 3}},
						{9, new int[] {2, 3}},
					}
				},
				{
					ElementTypeEnum.kPyramid5,
					new Dictionary<int, int[]>
					{
						{0, new int[] {3, 1, 4}},
						{1, new int[] {0, 2, 4}},
						{2, new int[] {1, 3, 4}},
						{3, new int[] {2, 0, 4}},
						{4, new int[] {0, 1, 2, 3}},
					}
				},
				{
					ElementTypeEnum.kPyramid13,
					new Dictionary<int, int[]>
					{
						{0, new int[] {3, 1, 4}},
						{1, new int[] {0, 2, 4}},
						{2, new int[] {1, 3, 4}},
						{3, new int[] {2, 0, 4}},
						{4, new int[] {0, 1, 2, 3}},
						{5, new int[] {0, 1}},
						{6, new int[] {1, 2}},
						{7, new int[] {2, 3}},
						{8, new int[] {3, 0}},
						{9, new int[] {0, 4}},
						{10, new int[] {1, 4}},
						{11, new int[] {2, 4}},
						{12, new int[] {3, 4}},
					}
				},
				{
					ElementTypeEnum.kWedge6,
					new Dictionary<int, int[]>
					{
						{0, new int[] {2, 1, 3}},
						{1, new int[] {0, 2, 4}},
						{2, new int[] {1, 0, 5}},
						{3, new int[] {5, 4, 0}},
						{4, new int[] {3, 5, 1}},
						{5, new int[] {4, 3, 2}},
					}
				},
				{
					ElementTypeEnum.kWedge15,
					new Dictionary<int, int[]>
					{
						{0, new int[] {2, 1, 3}},
						{1, new int[] {0, 2, 4}},
						{2, new int[] {1, 0, 5}},
						{3, new int[] {5, 4, 0}},
						{4, new int[] {3, 5, 1}},
						{5, new int[] {4, 3, 2}},
						{6, new int[] {0, 1}},
						{7, new int[] {1, 2}},
						{8, new int[] {2, 0}},
						{9, new int[] {3, 4}},
						{10, new int[] {4, 5}},
						{11, new int[] {5, 3}},
						{12, new int[] {0, 3}},
						{13, new int[] {1, 4}},
						{14, new int[] {2, 5}},
					}
				},
				{
					ElementTypeEnum.kHex8,
					new Dictionary<int, int[]>
					{
						{0, new int[] {3, 1, 4}},
						{1, new int[] {0, 2, 5}},
						{2, new int[] {1, 3, 6}},
						{3, new int[] {2, 0, 7}},
						{4, new int[] {5, 0, 7}},
						{5, new int[] {1, 4, 6}},
						{6, new int[] {5, 7, 2}},
						{7, new int[] {6, 3, 4}}
					}
				},
				{
					ElementTypeEnum.kHex20,
					new Dictionary<int, int[]>
					{
						{0, new int[] {3, 1, 4}},
						{1, new int[] {0, 2, 5}},
						{2, new int[] {1, 3, 6}},
						{3, new int[] {2, 0, 7}},
						{4, new int[] {5, 0, 7}},
						{5, new int[] {1, 4, 6}},
						{6, new int[] {5, 7, 2}},
						{7, new int[] {6, 3, 4}},
						{8, new int[] {0, 1}},
						{9, new int[] {1, 2}},
						{10, new int[] {2, 3}},
						{11, new int[] {3, 0}},
						{12, new int[] {4, 5}},
						{13, new int[] {5, 6}},
						{14, new int[] {6, 7}},
						{15, new int[] {7, 4}},
						{16, new int[] {0, 4}},
						{17, new int[] {1, 5}},
						{18, new int[] {2, 6}},
						{19, new int[] {3, 7}},
					}
				},
			};

		private Dictionary<ElementTypeEnum, CellType3D> AnsysMSolveElementDictionary =
			new Dictionary<ElementTypeEnum, CellType3D>
			{
				{ElementTypeEnum.kHex8, CellType3D.Hexa8},
				{ElementTypeEnum.kHex20, CellType3D.Hexa20},
				{ElementTypeEnum.kTet4, CellType3D.Tet4},
				{ElementTypeEnum.kTet10, CellType3D.Tet10},
				{ElementTypeEnum.kWedge6, CellType3D.Wedge6},
				{ElementTypeEnum.kWedge15, CellType3D.Wedge15},
				{ElementTypeEnum.kPyramid5, CellType3D.Pyra5},
				{ElementTypeEnum.kPyramid13, CellType3D.Pyra13},
			};


		private readonly IMechanicalExtAPI _api;

		public MainMSolve(IExtAPI api)
		{
			_api = api as IMechanicalExtAPI;
		}

		public bool SolveStatic(IUserSolver userSolver, Func<int, string, bool> func)
		{
			var solver = userSolver as IMechanicalUserSolver;
			try
			{
				var model = new Model();
				model.SubdomainsDictionary.Add(0, new Subdomain() {ID = 0});

				solver.Analysis.MeshData.Nodes.ToList()
					.ForEach(n => model.NodesDictionary.Add(n.Id, new Node3D(n.Id, n.X, n.Y, n.Z)));
				var mat = (_api.DataModel.GeoData.Assemblies[0].Parts[0].Bodies[0] as IGeoBody).Material;
				System.IO.File.WriteAllText(@"C:\Users\Dimitris\Desktop\ANSYS Models\materials.json", JsonConvert.SerializeObject(mat));

				var elementsList = new List<Element>();
				var factory = new ContinuumElement3DFactory(null, null);

				System.IO.File.WriteAllText(@"C:\Users\Dimitris\Desktop\ANSYS Models\elements.json", JsonConvert.SerializeObject(solver.Analysis.MeshData.Elements));
				
				//foreach (var ansysElement in solver.Analysis.MeshData.Elements)
				//{
				//	var elementNodes = new List<Node3D>();
				//	ansysElement.NodeIds.ToList().ForEach(id => elementNodes.Add((Node3D) model.NodesDictionary[id]));
				//	var element = factory.CreateElement(AnsysMSolveElementDictionary[ansysElement.Type], elementNodes);
				//	var elementWrapper = new Element() {ID = ansysElement.Id, ElementType = element};
				//	foreach (Node node in element.Nodes) elementWrapper.AddNode(node);
				//	model.ElementsDictionary.Add(ansysElement.Id, elementWrapper);
				//	model.SubdomainsDictionary[0].ElementsDictionary.Add(ansysElement.Id, elementWrapper);
				//}

				//var staticStructural = _api.DataModel.Project.Model.Analyses[0];

				//model.ConnectDataStructures();
				//var linearSystems = new Dictionary<int, ILinearSystem>();
				//linearSystems[0] = new SkylineLinearSystem(0, model.Subdomains[0].Forces);
				//SolverSkyline solverSkyline = new SolverSkyline(linearSystems[0]);
				//ProblemStructural provider = new ProblemStructural(model, linearSystems);
				//LinearAnalyzer childAnalyzer = new LinearAnalyzer(solverSkyline, linearSystems);
				//StaticAnalyzer parentAnalyzer = new StaticAnalyzer(provider, childAnalyzer, linearSystems);
				//parentAnalyzer.BuildMatrices();
				//parentAnalyzer.Initialize();
				//parentAnalyzer.Solve();
			}
			catch (Exception e)
			{
				System.IO.File.WriteAllText(@"C:\Users\Dimitris\Desktop\ANSYS Models\error.json", JsonConvert.SerializeObject(e));
				return false;
			}

			return true;
		}

		public virtual IEnumerable<double> NodeLoadValues(IUserLoad load, IEnumerable<int> nodeIds)
		{
			IMeshData mesh = ((IMechanicalUserLoad) load).Analysis.MeshData;
			var valueX = Convert.ToDouble(load.Properties["ValueX"]);
			var valueY = Convert.ToDouble(load.Properties["ValueY"]);
			var valueZ = Convert.ToDouble(load.Properties["ValueZ"]);

			var a= nodeIds.Select(nodeId => mesh.NodeById(nodeId))
				.Select(node => Math.Sqrt(valueX * valueX + valueY * valueY + valueZ * valueZ)).ToList();

			System.IO.File.WriteAllText(@"C:\Users\Dimitris\Desktop\ANSYS Models\loadNodes.json", JsonConvert.SerializeObject(a));
			return a;
		}

		private List<TempLoad> _loads = new List<TempLoad>();

		public virtual void WriteInitialLoadValues(IUserLoad load, StringWriter filename)
		{
			var res = new List<double>();
			IMeshData mesh = ((IMechanicalUserLoad) load).Analysis.MeshData;
			var valueX = Convert.ToDouble(load.Properties["ValueX"]);
			var valueY = Convert.ToDouble(load.Properties["ValueY"]);
			var valueZ = Convert.ToDouble(load.Properties["ValueZ"]);

			var propGeo = load.Properties["Geometry"];
			var refIds = (propGeo.Value as Surface).Ids;

			foreach (var refId in refIds)
			{
				var meshRegion = mesh.MeshRegionById(refId);
				var nodeIds = meshRegion.NodeIds;
				foreach (var nodeId in nodeIds)
				{
					INode node = mesh.NodeById(nodeId);
					_loads.Add(new TempLoad {Amount = valueX, DofType = DOFType.X, nodeId = nodeId});
					_loads.Add(new TempLoad {Amount = valueY, DofType = DOFType.Y, nodeId = nodeId});
					_loads.Add(new TempLoad {Amount = valueZ, DofType = DOFType.Z, nodeId = nodeId});
				}
			}
		}

		public virtual IEnumerable<double> NodeBoundaryValues(IUserLoad load, IEnumerable<int> nodeIds)
		{
			IMeshData mesh = ((IMechanicalUserLoad) load).Analysis.MeshData;
			var b= nodeIds.Select(nodeId => mesh.NodeById(nodeId)).Select(node => 0.0).ToList();

			System.IO.File.WriteAllText(@"C:\Users\Dimitris\Desktop\ANSYS Models\boundaryNodes.json", JsonConvert.SerializeObject(b));

			return b;
		}

		private Dictionary<int, List<DOFType>> _nodeConstrains = new Dictionary<int, List<DOFType>>();

		public virtual void WriteInitialBoundaryValues(IUserLoad load, StringWriter filename)
		{
			var res = new List<double>();
			IMeshData mesh = ((IMechanicalUserLoad) load).Analysis.MeshData;
			var valueX = Convert.ToInt32(load.Properties["ConstraintX"]);
			var valueY = Convert.ToInt32(load.Properties["ConstraintX"]);
			var valueZ = Convert.ToInt32(load.Properties["ConstraintZ"]);

			var propGeo = load.Properties["Geometry"];
			var refIds = (propGeo.Value as Surface).Ids;

			foreach (var refId in refIds)
			{
				var meshRegion = mesh.MeshRegionById(refId);
				var nodeIds = meshRegion.NodeIds;
				foreach (var nodeId in nodeIds)
				{
					if (!_nodeConstrains.ContainsKey(nodeId))
						_nodeConstrains.Add(nodeId, new List<DOFType>());
					if (valueX == 1) _nodeConstrains[nodeId].Add(DOFType.X);
					if (valueY == 1) _nodeConstrains[nodeId].Add(DOFType.Y);
					if (valueZ == 1) _nodeConstrains[nodeId].Add(DOFType.Z);
					_nodeConstrains[nodeId].Distinct();
				}
			}
		}
	}

	public class TempLoad
	{
		public double Amount { get; set; }
		public DOFType DofType { get; set; }
		public int nodeId { get; set; }
	}
}