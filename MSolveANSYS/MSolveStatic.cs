using Ansys.ACT.Interfaces.Mechanical;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;
using System.Text;
using Ansys.ACT.Automation.Mechanical;
using Ansys.ACT.Interfaces.Common;
using Ansys.ACT.Interfaces.Geometry;
using Ansys.ACT.Interfaces.UserObject;
using Ansys.EngineeringData.Material;
using Newtonsoft.Json;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Entities;
using Model = ISAAR.MSolve.FEM.Entities.Model;

namespace MSolveANSYS
{
	public class MSolveStatic
	{
		private readonly IMechanicalExtAPI _api;
		private readonly IMechanicalUserSolver _solver;

		public MSolveStatic(IExtAPI api,IUserSolver solver)
		{
			_api = api as IMechanicalExtAPI;
			_solver= solver as IMechanicalUserSolver;
		}

		public bool onsolve(IUserSolver userSolver, Func<int, string, bool> func)
		{
			var solver = userSolver as IMechanicalUserSolver;
			try
			{
				var model = new Model();
				model.SubdomainsDictionary.Add(0, new Subdomain() { ID = 0 });

				solver.Analysis.MeshData.Nodes.ToList()
					.ForEach(n => model.NodesDictionary.Add(n.Id, new Node3D(n.Id, n.X, n.Y, n.Z)));
				var mat = (_api.DataModel.GeoData.Assemblies[0].Parts[0].Bodies[0] as IGeoBody).Material as MaterialClass;
				
				var materialProperties = GetListMaterialProperties(mat);

				System.IO.File.WriteAllText(@"C:\Users\Dimitris\Desktop\ANSYS Models\materialProperties1.json",
					JsonConvert.SerializeObject(materialProperties));

				//System.IO.File.WriteAllText(@"C:\Users\Dimitris\Desktop\ANSYS Models\materials.json", mat.ToString());

				//var elementsList = new List<Element>();
				//var factory = new ContinuumElement3DFactory(null, null);

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

		private List<string> GetListMaterialProperties(MaterialClass material)
		{

			var a = material.MaterialProperties;
			
			return material.GetType().GetProperties().Select(p=>p.Name).ToList();
		}

		private object GetMaterialPropertyByName(object material, string propertyName)
		{
			return material.GetType().GetProperties().Single(p => p.Name == propertyName).GetValue(material, null);
		}

	}
}
