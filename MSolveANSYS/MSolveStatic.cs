﻿using Ansys.ACT.Interfaces.Mechanical;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;
using System.Text;
using Ansys.ACT.Automation.Mechanical;
using Ansys.ACT.Automation.Mechanical.BoundaryConditions;
using Ansys.ACT.Interfaces.Common;
using Ansys.ACT.Interfaces.Geometry;
using Ansys.ACT.Interfaces.Mesh;
using Ansys.ACT.Interfaces.Post;
using Ansys.ACT.Interfaces.UserObject;
using Ansys.Core.Units;
using Ansys.EngineeringData.Material;
using ISAAR.MSolve.Discretization.Interfaces;
using Newtonsoft.Json;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Materials;
using ISAAR.MSolve.Geometry.Shapes;
using Model = ISAAR.MSolve.FEM.Entities.Model;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.Solvers.Skyline;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;

namespace MSolveANSYS
{
	public class MSolveStatic
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

		private Dictionary<ElementTypeEnum, int[]> AnsysMSolveElementLocalCoordinates =
			new Dictionary<ElementTypeEnum, int[]>()
			{
				{ElementTypeEnum.kHex8, new[] {4,5,1,0,7,6,2,3}}
			};

		private readonly IMechanicalExtAPI _api;
		private readonly IMechanicalUserSolver _solver;

		public MSolveStatic(IExtAPI api,IUserSolver solver)
		{
			_api = api as IMechanicalExtAPI;
			_solver= solver as IMechanicalUserSolver;
		}


		public Model model;
		public IVector SolutionVector;
		public bool onsolve(IUserSolver userSolver, Func<int, string, bool> func)
		{
			var solver = userSolver as IMechanicalUserSolver;
			try
			{
				VectorExtensions.AssignTotalAffinityCount();
				model = new Model();
				model.SubdomainsDictionary.Add(0, new Subdomain() { ID = 0 });

				solver.Analysis.MeshData.Nodes.ToList()
					.ForEach(n => model.NodesDictionary.Add(n.Id, new Node3D(n.Id, n.X, n.Y, n.Z)));
				
				var mat = (_api.DataModel.GeoData.Assemblies[0].Parts[0].Bodies[0] as IGeoBody).Material as MaterialClass;
				var dataDictionary=GetMaterialPropertyByName(mat, "Elasticity");
				var material = new ElasticMaterial3D()
				{
					YoungModulus = dataDictionary["Young's Modulus"].Value,
					PoissonRatio = dataDictionary["Poisson's Ratio"].Value
				};

				var elementsList = new List<Element>();
				var factory = new ContinuumElement3DFactory(material, null);
				

				foreach (var ansysElement in solver.Analysis.MeshData.Elements)
				{
					var ansysNodes = new List<Node3D>();
					ansysElement.NodeIds.ToList().ForEach(id => ansysNodes.Add((Node3D)model.NodesDictionary[id]));
					var msolveNodes = RenumberNodesFromAnsysToMSolve(ansysNodes, ansysElement.Type);
					var element = factory.CreateElement(AnsysMSolveElementDictionary[ansysElement.Type], msolveNodes);
					var elementWrapper = new Element() { ID = ansysElement.Id, ElementType = element };
					foreach (var node in element.Nodes) elementWrapper.AddNode(node);
					model.ElementsDictionary.Add(ansysElement.Id, elementWrapper);
					model.SubdomainsDictionary[0].ElementsDictionary.Add(ansysElement.Id, elementWrapper);
				}

				var force = _api.Application.InvokeUIThread(()=>_api.DataModel.Project.Model.Analyses[0].Children[1]) as Force;

				var forceLocation =_api.Application.InvokeUIThread(() => force.Location) as ISelectionInfo;
				var xValue = _api.Application.InvokeUIThread(() => force.XComponent.Output.DiscreteValues[0]) as Quantity;
				var yValue = _api.Application.InvokeUIThread(() => force.YComponent.Output.DiscreteValues[0]) as Quantity;
				var zValue = _api.Application.InvokeUIThread(() => force.ZComponent.Output.DiscreteValues[0]) as Quantity;
				var magnitude = _api.Application.InvokeUIThread(() => force.Magnitude.Output.DiscreteValues[1]) as Quantity;
				var forceSurfaceId = forceLocation.Ids[0];
				var forceNodes=solver.Analysis.MeshData.MeshRegionById(forceSurfaceId).Nodes;

				foreach (var node in forceNodes)
				{
					model.Loads.Add(new Load()
					{
						Amount = xValue.Value * magnitude.Value / forceNodes.Count,
						Node = model.NodesDictionary[node.Id],
						DOF = DOFType.X
					});
					model.Loads.Add(new Load()
					{
						Amount = yValue.Value * magnitude.Value / forceNodes.Count,
						Node = model.NodesDictionary[node.Id],
						DOF = DOFType.Y
					});
					model.Loads.Add(new Load()
					{
						Amount = zValue.Value * magnitude.Value / forceNodes.Count,
						Node = model.NodesDictionary[node.Id],
						DOF = DOFType.Z
					});
				}


				var fixedSupport = _api.Application.InvokeUIThread(() => _api.DataModel.Project.Model.Analyses[0].Children[2]) as FixedSupport;
				var fixedLocation = _api.Application.InvokeUIThread(() => fixedSupport.Location) as ISelectionInfo;
				var fixedSurfaceId = fixedLocation.Ids[0];
				var fixedNodes = solver.Analysis.MeshData.MeshRegionById(fixedSurfaceId).Nodes;

				foreach (var node in fixedNodes)
				{
					model.NodesDictionary[node.Id].Constraints.Add(DOFType.X);
					model.NodesDictionary[node.Id].Constraints.Add(DOFType.Y);
					model.NodesDictionary[node.Id].Constraints.Add(DOFType.Z);
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

				SolutionVector = linearSystems[0].Solution;
			}
			catch (Exception e)
			{
				System.IO.File.WriteAllText(@"C:\Users\Dimitris\Desktop\ANSYS Models\error.json", JsonConvert.SerializeObject(e));
				return false;
			}

			return true;
		}

		private List<Node3D> RenumberNodesFromAnsysToMSolve(IReadOnlyList<Node3D> ansysNodes, ElementTypeEnum elementType)
		{
			var connectivity = AnsysMSolveElementLocalCoordinates[elementType];
			return ansysNodes.Select((t, i) => ansysNodes[connectivity[i]]).ToList();
		}

		/// <summary>
		/// Return the list of material property names for the given material.
		/// </summary>
		/// <param name="material"></param>
		/// <returns></returns>
		private List<string> GetListMaterialProperties(MaterialClass material)
		{
			var propertyNames = new List<string>();
			for (int i = 0; i < material.MaterialProperties.Count; i++)
			{
				propertyNames.Add(material.MaterialProperties.Get(i).TypeName);
			}

			return propertyNames;
		}

		private Dictionary<string, (string Unit, double Value)> GetMaterialPropertyByName(MaterialClass material, string propertyName)
		{
			var properties = new List<IMaterialProperty>();
			for (var i = 0; i < material.MaterialProperties.Count; i++)
			{
				properties.Add(material.MaterialProperties.Get(i));
			}

			var property= properties.Single(p => p.TypeName == propertyName);
			
			var variablesDictionary=new Dictionary<string, (string unit ,double value)>();
			for (int i = 0; i < property.MaterialPropertyDatas.Get(0).Variables.Count; i++)
			{
				var variable = property.MaterialPropertyDatas.Get(0).Variables.Get(i);
				variablesDictionary.Add(variable.TypeName, (variable.Unit,variable.DatumColl.Get(0).Value));
			}

			return variablesDictionary;
		}

		public IEnumerable<string> getreader(IUserSolver userSolver)
		{
			return new[]
			{
				"MSolveANSYS.MSolveStaticReader" ,
				JsonConvert.SerializeObject(model.NodalDOFsDictionary),
				JsonConvert.SerializeObject(SolutionVector),
			};
		}
	}
}
