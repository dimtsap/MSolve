using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Ansys.ACT.Automation.Mechanical;
using Ansys.ACT.Automation.Mechanical.AnalysisSettings;
using Ansys.ACT.Automation.Mechanical.BoundaryConditions;
using Ansys.ACT.Interfaces.Common;
using Ansys.ACT.Interfaces.Geometry;
using Ansys.ACT.Interfaces.Mechanical;
using Ansys.ACT.Interfaces.Mesh;
using Ansys.ACT.Interfaces.UserObject;
using Ansys.EngineeringData.Material;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Materials;
using ISAAR.MSolve.Geometry.Shapes;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.Solvers.Skyline;
using Newtonsoft.Json;
using Model = ISAAR.MSolve.FEM.Entities.Model;

namespace MSolveANSYS
{
	public class MSolveDynamic
	{
		private readonly Dictionary<ElementTypeEnum, CellType3D> _ansysMSolveElementDictionary =
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

		private readonly Dictionary<ElementTypeEnum, int[]> _ansysMSolveElementLocalCoordinates =
			new Dictionary<ElementTypeEnum, int[]>()
			{
				{ElementTypeEnum.kHex8, new[] {4,5,1,0,7,6,2,3}},
				{ElementTypeEnum.kHex20, new int[20] {4,5,1,0,7,6,2,3,12,16,15,17,13,8,9,11,14,19,18,10}},
				{ElementTypeEnum.kTet4, new int[4]{0,3,1,2} },
				{ElementTypeEnum.kTet10, new int[10]{0,3,1,2,7,8,4,6,5,9} },
				{ElementTypeEnum.kWedge6, new int[6]{5,4,3,2,1,0} },
				{ElementTypeEnum.kWedge15, new int[15]{5,4,3,2,1,0,10,11,14,9,13,12,7,8,6} }
			};

		private readonly IMechanicalExtAPI _api;
		private readonly IMechanicalUserSolver _solver;

		public MSolveDynamic(IExtAPI api, IUserSolver solver)
		{
			_api = api as IMechanicalExtAPI;
			_solver = solver as IMechanicalUserSolver;
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

				var dynamicMaterial= new DynamicMaterial(1,0.05,0.05);

				var material = AnsysUtilities.GetAnsysMaterial(_api);
				AnsysUtilities.GenerateElements(material, dynamicMaterial, solver, model);
				AnsysUtilities.CalculateAcceleration(_api,solver,model);
				//AnsysUtilities.CalculateForces(_api, solver, model);
				AnsysUtilities.ImposeFixedSupport(_api, solver, model);

				model.ConnectDataStructures();

				var a=_api.Application.InvokeUIThread(()=>_api.DataModel.Project.Model.Analyses[0].PropertyNames) as List<string>;

				var linearSystems = new Dictionary<int, ILinearSystem>();
				linearSystems[0] = new SkylineLinearSystem(0, model.Subdomains[0].Forces);
				SolverSkyline solverSkyline = new SolverSkyline(linearSystems[0]);
				ProblemStructural provider = new ProblemStructural(model, linearSystems);
				LinearAnalyzer childAnalyzer = new LinearAnalyzer(solverSkyline, linearSystems);
				NewmarkDynamicAnalyzer parentAnalyzer = new NewmarkDynamicAnalyzer(provider, childAnalyzer, linearSystems, 0.6, 1, 0.02, 53.74);
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
