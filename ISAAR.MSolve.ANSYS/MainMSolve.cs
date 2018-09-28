using System;
using System.Collections.Generic;
using System.Text;
using Ansys.ACT.Interfaces.Common;
using Ansys.ACT.Interfaces.Mechanical;
using Ansys.ACT.Interfaces.Mesh;
using Ansys.ACT.Interfaces.UserObject;
using ISAAR.MSolve.FEM.Entities;
using System.Linq;
using Ansys.ACT.Interfaces.Geometry;
using Ansys.ACT.Interfaces.Project;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.Geometry.Shapes;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.Solvers.Skyline;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Analyzers;

namespace ISAAR.MSolve.ANSYS
{
	public class MainMSolve
	{
		private Dictionary<ElementTypeEnum, Dictionary<int, int[]>> ElementTypesNodesDictionary =new Dictionary<ElementTypeEnum, Dictionary<int, int[]>>()
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
		private Dictionary<ElementTypeEnum, CellType3D> AnsysMSolveElementDictionary= new Dictionary<ElementTypeEnum, CellType3D>
		{
			{ElementTypeEnum.kHex8, CellType3D.Hexa8 },
			{ElementTypeEnum.kHex20, CellType3D.Hexa20 },
			{ElementTypeEnum.kTet4, CellType3D.Tet4 },
			{ElementTypeEnum.kTet10, CellType3D.Tet10 },
			{ElementTypeEnum.kWedge6, CellType3D.Wedge6 },
			{ElementTypeEnum.kWedge15, CellType3D.Wedge15 },
			{ElementTypeEnum.kPyramid5, CellType3D.Pyra5 },
			{ElementTypeEnum.kPyramid13, CellType3D.Pyra13 },
		};


		private readonly IMechanicalExtAPI _api;
		private readonly IMechanicalUserSolver _solver;

		public MainMSolve(IExtAPI api, IUserSolver solver)
		{
			_api = api as IMechanicalExtAPI;
			_solver = solver as IMechanicalUserSolver;
		}

		public bool Solve(IUserResult solver, Func<int, string, bool> func)
		{
			try
			{
				var model = new Model();
				model.SubdomainsDictionary.Add(0, new Subdomain() { ID = 0 });

				_solver.Analysis.MeshData.Nodes.ToList().ForEach(n => model.NodesDictionary.Add(n.Id,new Node3D(n.Id, n.X, n.Y, n.Z)));
				//TODO:Materials
				var mat = (_api.DataModel.GeoData.Assemblies[0].Parts[0].Bodies[0] as IGeoBody).Material.ToString();
				Console.WriteLine(mat);

				var elementsList= new List<Element>();
				var factory = new ContinuumElement3DFactory(null, null);
				foreach (var ansysElement in _solver.Analysis.MeshData.Elements)
				{
					var elementNodes= new List<Node3D>();
					ansysElement.NodeIds.ToList().ForEach(id=>elementNodes.Add((Node3D)model.NodesDictionary[id]));
					var element=factory.CreateElement(AnsysMSolveElementDictionary[ansysElement.Type],elementNodes);
					var elementWrapper = new Element() { ID = ansysElement.Id, ElementType = element };
					foreach (Node node in element.Nodes) elementWrapper.AddNode(node);
					model.ElementsDictionary.Add(ansysElement.Id, elementWrapper);
					model.SubdomainsDictionary[0].ElementsDictionary.Add(ansysElement.Id, elementWrapper);
				}

				var staticStructural = _api.DataModel.Project.Model.Analyses[0];

				//TODO: find load and constrains
				var loads= staticStructural.GetLoadObjects("ExtensionName").ToList();
				foreach (var load in loads)
				{
					var userLoad=load as IMechanicalUserLoad;
					if (userLoad.IsLoad)
					{
					}
					else if (userLoad.IsSupport)
					{

					}
				}


				model.ConnectDataStructures();
				var linearSystems = new Dictionary<int, ILinearSystem>();
				linearSystems[0] = new SkylineLinearSystem(0, model.Subdomains[0].Forces);
				SolverSkyline solverSkyline = new SolverSkyline(linearSystems[0]);
				ProblemStructural provider = new ProblemStructural(model, linearSystems);
				LinearAnalyzer childAnalyzer = new LinearAnalyzer(solverSkyline, linearSystems);
				StaticAnalyzer parentAnalyzer = new StaticAnalyzer(provider, childAnalyzer, linearSystems);
				parentAnalyzer.BuildMatrices();
				parentAnalyzer.Initialize();
				parentAnalyzer.Solve();

			}
			catch (Exception e)
			{
				Console.WriteLine("Error:");
				Console.WriteLine(e.Message);
				return false;
			}

			return true;
		}
	}
}