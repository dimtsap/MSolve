using System;
using System.Collections.Generic;
using System.Linq;
using Ansys.ACT.Automation.Mechanical;
using Ansys.ACT.Automation.Mechanical.BoundaryConditions;
using Ansys.ACT.Interfaces.Common;
using Ansys.ACT.Interfaces.Geometry;
using Ansys.ACT.Interfaces.Mechanical;
using Ansys.ACT.Interfaces.Mesh;
using Ansys.Core.Units;
using Ansys.EngineeringData.Material;
using Ansys.Mechanical.DataModel.Enums;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Mesh;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Materials;
using Model = ISAAR.MSolve.FEM.Entities.Model;

namespace AnsysMSolve
{
	public static class AnsysUtilities
	{
		private static readonly Dictionary<ElementTypeEnum, CellType> _ansysMSolveElementDictionary =
			new Dictionary<ElementTypeEnum, CellType>
			{
				{ElementTypeEnum.kHex8, CellType.Hexa8},
				{ElementTypeEnum.kHex20, CellType.Hexa20},
				{ElementTypeEnum.kTet4, CellType.Tet4},
				{ElementTypeEnum.kTet10, CellType.Tet10},
				{ElementTypeEnum.kWedge6, CellType.Wedge6},
				{ElementTypeEnum.kWedge15, CellType.Wedge15},
				{ElementTypeEnum.kPyramid5, CellType.Pyra5},
				{ElementTypeEnum.kPyramid13, CellType.Pyra13},
			};

		private static readonly Dictionary<ElementTypeEnum, int[]> _ansysMSolveElementLocalCoordinates =
			new Dictionary<ElementTypeEnum, int[]>()
			{
				{ElementTypeEnum.kHex8, new[] {4,5,1,0,7,6,2,3}},
				{ElementTypeEnum.kHex20, new int[20] {4,5,1,0,7,6,2,3,12,16,15,17,13,8,9,11,14,19,18,10}},
				{ElementTypeEnum.kTet4, new int[4]{0,3,1,2} },
				{ElementTypeEnum.kTet10, new int[10]{0,3,1,2,7,8,4,6,5,9} },
				{ElementTypeEnum.kWedge6, new int[6]{5,4,3,2,1,0} },
				{ElementTypeEnum.kWedge15, new int[15]{5,4,3,2,1,0,10,11,14,9,13,12,7,8,6} }
			};

        private static readonly Table<ElementTypeEnum, int, int[]> _elementFacesNodes =
            new Table<ElementTypeEnum, int, int[]>()
            {
                [ElementTypeEnum.kTet4, 1] = new int[] {0, 1, 3},
                [ElementTypeEnum.kTet4, 2] = new int[] {1, 2, 3},
                [ElementTypeEnum.kTet4, 3] = new int[] {2, 0, 3},
                [ElementTypeEnum.kTet4, 4] = new int[] {0, 2, 1},

                [ElementTypeEnum.kTet10, 1] = new int[] {0, 1, 3, 4, 8, 7},
                [ElementTypeEnum.kTet10, 2] = new int[] {1, 2, 3, 5, 9, 8},
                [ElementTypeEnum.kTet10, 3] = new int[] {2, 0, 3, 6, 7, 9},
                [ElementTypeEnum.kTet10, 4] = new int[] {0, 2, 1, 6, 5, 4},

                [ElementTypeEnum.kHex8, 1] = new int[] {0, 1, 5, 4},
                [ElementTypeEnum.kHex8, 2] = new int[] {1, 2, 6, 5},
                [ElementTypeEnum.kHex8, 3] = new int[] {2, 3, 7, 6},
                [ElementTypeEnum.kHex8, 4] = new int[] {3, 0, 4, 7},
                [ElementTypeEnum.kHex8, 5] = new int[] {0, 3, 2, 1},
                [ElementTypeEnum.kHex8, 6] = new int[] {4, 5, 6, 7},

                [ElementTypeEnum.kHex20, 1] = new int[] {0, 1, 5, 4, 8, 17, 12, 16},
                [ElementTypeEnum.kHex20, 2] = new int[] {1, 2, 6, 5, 9, 18, 13, 17},
                [ElementTypeEnum.kHex20, 3] = new int[] {2, 3, 7, 6, 10, 19, 14, 18},
                [ElementTypeEnum.kHex20, 4] = new int[] {3, 0, 4, 7, 11, 16, 15, 19},
                [ElementTypeEnum.kHex20, 5] = new int[] {0, 3, 2, 1, 11, 10, 9, 8},
                [ElementTypeEnum.kHex20, 6] = new int[] {4, 5, 6, 7, 12, 13, 14, 15},

                [ElementTypeEnum.kWedge6, 1] = new int[] {0, 2, 1},
                [ElementTypeEnum.kWedge6, 2] = new int[] {3, 4, 5},
                [ElementTypeEnum.kWedge6, 3] = new int[] {0, 1, 4, 3},
                [ElementTypeEnum.kWedge6, 4] = new int[] {1, 2, 5, 4},
                [ElementTypeEnum.kWedge6, 5] = new int[] {0, 3, 5, 2},

                [ElementTypeEnum.kWedge15, 1] = new int[] {0, 2, 1, 8, 7, 6},
                [ElementTypeEnum.kWedge15, 2] = new int[] {3, 4, 5, 9, 10, 11},
                [ElementTypeEnum.kWedge15, 3] = new int[] {0, 1, 4, 3, 6, 13, 9, 12},
                [ElementTypeEnum.kWedge15, 4] = new int[] {1, 2, 5, 4, 7, 14, 10, 13},
                [ElementTypeEnum.kWedge15, 5] = new int[] {0, 3, 5, 2, 12, 11, 14, 8},

                [ElementTypeEnum.kPyramid5, 1] = new int[] {0, 3, 2, 1},
                [ElementTypeEnum.kPyramid5, 2] = new int[] {0, 1, 4},
                [ElementTypeEnum.kPyramid5, 3] = new int[] {1, 2, 4},
                [ElementTypeEnum.kPyramid5, 4] = new int[] {3, 4, 2},
                [ElementTypeEnum.kPyramid5, 5] = new int[] {0, 4, 3},

                [ElementTypeEnum.kPyramid13, 1] = new int[] {0, 3, 2, 1, 8, 7, 6, 5},
                [ElementTypeEnum.kPyramid13, 1] = new int[] {0, 1, 4, 5, 10, 9},
                [ElementTypeEnum.kPyramid13, 1] = new int[] {1, 2, 4, 6, 11, 10},
                [ElementTypeEnum.kPyramid13, 1] = new int[] {3, 4, 2, 12, 11, 7},
                [ElementTypeEnum.kPyramid13, 1] = new int[] {0, 4, 3, 9, 12, 8},
            };

        

		public static ElasticMaterial3D GetAnsysMaterial(IMechanicalExtAPI _api)
		{
			var mat = (_api.DataModel.GeoData.Assemblies[0].Parts[0].Bodies[0] as IGeoBody).Material as MaterialClass;
			var dataDictionary = GetMaterialPropertyByName(mat, "Elasticity");
			var material = new ElasticMaterial3D()
			{
				YoungModulus = dataDictionary["Young's Modulus"].Value,
				PoissonRatio = dataDictionary["Poisson's Ratio"].Value
			};
			return material;
		}

		private static List<string> GetListMaterialProperties(MaterialClass material)
		{
			var propertyNames = new List<string>();
			for (int i = 0; i < material.MaterialProperties.Count; i++)
			{
				propertyNames.Add(material.MaterialProperties.Get(i).TypeName);
			}

			return propertyNames;
		}

		private static Dictionary<string, (string Unit, double Value)> GetMaterialPropertyByName(MaterialClass material, string propertyName)
		{
			var properties = new List<IMaterialProperty>();
			for (var i = 0; i < material.MaterialProperties.Count; i++)
			{
				properties.Add(material.MaterialProperties.Get(i));
			}

			var property = properties.Single(p => p.TypeName == propertyName);

			var variablesDictionary = new Dictionary<string, (string unit, double value)>();
			for (int i = 0; i < property.MaterialPropertyDatas.Get(0).Variables.Count; i++)
			{
				var variable = property.MaterialPropertyDatas.Get(0).Variables.Get(i);
				variablesDictionary.Add(variable.TypeName, (variable.Unit, variable.DatumColl.Get(0).Value));
			}

			return variablesDictionary;
		}

		private static List<Node> RenumberNodesFromAnsysToMSolve(IReadOnlyList<Node> ansysNodes, ElementTypeEnum elementType)
		{
			var connectivity = _ansysMSolveElementLocalCoordinates[elementType];
			return ansysNodes.Select((t, i) => ansysNodes[connectivity[i]]).ToList();
		}


		public static void ImposeFixedSupport(IMechanicalExtAPI _api,IMechanicalUserSolver solver, Model model)
		{
			var ansysFixedSupports =
				_api.Application.InvokeUIThread(() => _api.DataModel.Project.Model.Analyses[0].Children
						.Where(c => c.GetType() == typeof(FixedSupport)).ToList()) as List<DataModelObject>;

			foreach (var ansysFixedSupport in ansysFixedSupports)
			{
				var fixedSupport = ansysFixedSupport as FixedSupport;
				var fixedLocation = _api.Application.InvokeUIThread(() => fixedSupport.Location) as ISelectionInfo;
				var fixedSurfaceId = fixedLocation.Ids[0];
				var fixedNodes = solver.Analysis.MeshData.MeshRegionById(fixedSurfaceId).Nodes;

				foreach (var node in fixedNodes)
				{
					model.NodesDictionary[node.Id].Constraints.Add(new Constraint(){DOF = StructuralDof.TranslationX});
                    model.NodesDictionary[node.Id].Constraints.Add(new Constraint(){DOF = StructuralDof.TranslationY});
                    model.NodesDictionary[node.Id].Constraints.Add(new Constraint(){DOF = StructuralDof.TranslationZ});
				}
			}
		}

		public static void GenerateElements(ElasticMaterial3D material, DynamicMaterial dynamicMaterial, IMechanicalUserSolver solver, Model model)
		{
			var factory = new ContinuumElement3DFactory(material, dynamicMaterial);
			foreach (var ansysElement in solver.Analysis.MeshData.Elements)
			{
				var ansysNodes = new List<Node>();
				ansysElement.NodeIds.ToList().ForEach(id => ansysNodes.Add(model.NodesDictionary[id]));
				var msolveNodes = RenumberNodesFromAnsysToMSolve(ansysNodes, ansysElement.Type);
				var element = factory.CreateElement(_ansysMSolveElementDictionary[ansysElement.Type], msolveNodes);
				var elementWrapper = new Element() { ID = ansysElement.Id, ElementType = element };
				foreach (var node in element.Nodes) elementWrapper.AddNode(node);
				model.ElementsDictionary.Add(ansysElement.Id, elementWrapper);
				model.SubdomainsDictionary[0].Elements.Add(elementWrapper);
			}
		}
		
		public static void CalculateForces(IMechanicalExtAPI _api,IMechanicalUserSolver solver, Model model)
		{
			var forces =
				_api.Application.InvokeUIThread(() => _api.DataModel.Project.Model.Analyses[0].Children
					.Where(c => c.GetType() == typeof(Force)).ToList()) as List<DataModelObject>;

			foreach (var ansysForce in forces)
			{
				var force = ansysForce as Force;
				var isParsed = Enum.TryParse<LoadDefineBy>(_api.Application.InvokeUIThread(() => force.DefineBy).ToString(), out var defineBy);
				if (defineBy == LoadDefineBy.Vector)
					CalculateVectorForce(_api,solver,model, force);
				else
					CalculateComponentsForce(_api,solver,model, force);
			}
		}
		
        public static void CalculatePressure(IMechanicalExtAPI _api,IMechanicalUserSolver solver, Model model)
        {
            solver.Analysis.MeshData.GetQuad4ExteriorFaces(out int[] quad4Elements);
            solver.Analysis.MeshData.GetQuad8ExteriorFaces(out int[] quad8Elements);
            solver.Analysis.MeshData.GetTri3ExteriorFaces(out int[] tri3Elements);
            solver.Analysis.MeshData.GetTri6ExteriorFaces(out int[] tri6Elements);

            List<ElementFace> quad4ExteriorFaces=ConvertToElementFacePair(quad4Elements);
            List<ElementFace> quad8ExteriorFaces=ConvertToElementFacePair(quad8Elements);
            List<ElementFace> tri3ExteriorFaces=ConvertToElementFacePair(tri3Elements);
            List<ElementFace> tri6ExteriorFaces=ConvertToElementFacePair(tri6Elements);

            var pressures =
                _api.Application.InvokeUIThread(() => _api.DataModel.Project.Model.Analyses[0].Children
                    .Where(c => c.GetType() == typeof(Pressure)).ToList()) as List<DataModelObject>;

            foreach (var ansysPressures in pressures)
            {
                var pressure = ansysPressures as Pressure;
                var isParsed = Enum.TryParse<LoadDefineBy>(_api.Application.InvokeUIThread(() => pressure.DefineBy).ToString(), out var defineBy);
				var pressureLocation = _api.Application.InvokeUIThread(() => pressure.Location) as ISelectionInfo;
                var pressureSurfaceId = pressureLocation.Ids[0];
                var elements = solver.Analysis.MeshData.MeshRegionById(pressureSurfaceId).Elements;

                switch (defineBy)
                {
                    case LoadDefineBy.NormalToOrTangential:
                    case LoadDefineBy.Vector:
                        CalculateVectorPressure(_api, model, elements, pressure, quad4ExteriorFaces, quad8ExteriorFaces,
                            tri3ExteriorFaces, tri6ExteriorFaces);
                        break;
					case LoadDefineBy.Components:
                        CalculateComponentsPressure(_api,model, elements, pressure, quad4ExteriorFaces, quad8ExteriorFaces,
                            tri3ExteriorFaces, tri6ExteriorFaces);
                        break;
                }
            }
        }

		private static void CalculateVectorPressure(IMechanicalExtAPI _api, Model model, IList<IElement> elements,
            Pressure pressure, List<ElementFace> quad4ExteriorFaces, List<ElementFace> quad8ExteriorFaces,
            List<ElementFace> tri3ExteriorFaces, List<ElementFace> tri6ExteriorFaces)
		{
            var magnitude = _api.Application.InvokeUIThread(() => pressure.Magnitude);

            var magnitudeValue = _api.Application.InvokeUIThread(() => pressure.Magnitude.Output.DiscreteValues[1]) as Quantity;

            foreach (var element in elements)
            {
                var a = pressure.Magnitude;
                var quads4 = quad4ExteriorFaces.Where(x => x.ElementId == element.Id).ToList();
                var quads8 = quad8ExteriorFaces.Where(x => x.ElementId == element.Id).ToList();
                var tri3 = tri3ExteriorFaces.Where(x => x.ElementId == element.Id).ToList();
                var tri6 = tri6ExteriorFaces.Where(x => x.ElementId == element.Id).ToList();

                foreach (var face in quads4)
                    CreateGeometryElement(CellType.Quad4, element, model);

                foreach (var face in quads8)
                    CreateGeometryElement(CellType.Quad8, element, model);

                foreach (var face in tri3)
                    CreateGeometryElement(CellType.Tri3, element, model);

                foreach (var face in tri6)
                    CreateGeometryElement(CellType.Tri6, element, model);
            }
		}

		private static void CalculateComponentsPressure(IMechanicalExtAPI _api,Model model, IList<IElement> elements, Pressure pressure, List<ElementFace> quad4ExteriorFaces,
            List<ElementFace> quad8ExteriorFaces, List<ElementFace> tri3ExteriorFaces, List<ElementFace> tri6ExteriorFaces)
        {
            var magnitude = _api.Application.InvokeUIThread(() => pressure.Magnitude);

            var magnitudeValue = _api.Application.InvokeUIThread(() => pressure.Magnitude.Output.DiscreteValues[1]) as Quantity;
            foreach (var element in elements)
            {
                var a = pressure.Magnitude;
                var quads4 = quad4ExteriorFaces.Where(x => x.ElementId == element.Id).ToList();
                var quads8 = quad8ExteriorFaces.Where(x => x.ElementId == element.Id).ToList();
                var tri3 = tri3ExteriorFaces.Where(x => x.ElementId == element.Id).ToList();
                var tri6 = tri6ExteriorFaces.Where(x => x.ElementId == element.Id).ToList();

                foreach (var face in quads4)
                    CreateGeometryElement(CellType.Quad4, element, model);

                foreach (var face in quads8)
                    CreateGeometryElement(CellType.Quad8, element, model);

                foreach (var face in tri3)
                    CreateGeometryElement(CellType.Tri3, element, model);

                foreach (var face in tri6)
                    CreateGeometryElement(CellType.Tri6, element, model);
            }
        }

        private static void CreateGeometryElement(CellType quad4, IElement element, Model model)
		{
		}

		private static List<ElementFace> ConvertToElementFacePair(int[] elements)
		{
			var list= new List<ElementFace>();
            if (elements == null) return list;
            for (int i = 0; i < elements.Length; i+=2)
            {
                list.Add(new ElementFace()
                {
					ElementId = elements[i],
					FaceId =elements[i+1]
                });
            }
            return list;
        }

		//public static void CalculateAcceleration(IMechanicalExtAPI _api, IMechanicalUserSolver solver, Model model)
		//{
		//	var accelerations = _api.Application.InvokeUIThread(() => _api.DataModel.Project.Model.Analyses[0].Children
		//		.Where(c => c.GetType() == typeof(Acceleration)).ToList()) as List<DataModelObject>;

		//	foreach (var ansysAcceleration in accelerations)
		//	{
		//		var acceleration = ansysAcceleration as Acceleration;
		//		var accelerationLocation =
		//			_api.Application.InvokeUIThread(() => acceleration.Location) as ISelectionInfo;

		//		var xValues = _api.Application.InvokeUIThread(() => acceleration.XComponent.Output.DiscreteValues) as List<Quantity>;
		//		var yValues = _api.Application.InvokeUIThread(() => acceleration.YComponent.Output.DiscreteValues) as List<Quantity>;
		//		var zValues = _api.Application.InvokeUIThread(() => acceleration.ZComponent.Output.DiscreteValues) as List<Quantity>;

		//		model.MassAccelerationHistoryLoads.Add(
		//			new MassAccelerationHistoryLoad(xValues.Select(v => v.Value).ToList()) {DOF =  StructuralDof.TranslationX});
		//		model.MassAccelerationHistoryLoads.Add(
		//			new MassAccelerationHistoryLoad(yValues.Select(v => v.Value).ToList()) { DOF =  StructuralDof.TranslationY });
		//		model.MassAccelerationHistoryLoads.Add(
		//			new MassAccelerationHistoryLoad(zValues.Select(v => v.Value).ToList()) { DOF =  StructuralDof.TranslationZ });
		//	}
		//}


		public static void CalculateComponentsForce(IMechanicalExtAPI _api,IMechanicalUserSolver solver,Model model, Force force)
		{
			var forceLocation = _api.Application.InvokeUIThread(() => force.Location) as ISelectionInfo;
			var xValue = _api.Application.InvokeUIThread(() => force.XComponent.Output.DiscreteValues[1]) as Quantity;
			var yValue = _api.Application.InvokeUIThread(() => force.YComponent.Output.DiscreteValues[1]) as Quantity;
			var zValue = _api.Application.InvokeUIThread(() => force.ZComponent.Output.DiscreteValues[1]) as Quantity;
			var forceSurfaceId = forceLocation.Ids[0];
			var forceNodes = solver.Analysis.MeshData.MeshRegionById(forceSurfaceId).Nodes;

			foreach (var node in forceNodes)
			{
				model.Loads.Add(new Load()
				{
					Amount = xValue.Value / forceNodes.Count,
					Node = model.NodesDictionary[node.Id],
					DOF =StructuralDof.TranslationX
				});
				model.Loads.Add(new Load()
				{
					Amount = yValue.Value / forceNodes.Count,
					Node = model.NodesDictionary[node.Id],
					DOF = StructuralDof.TranslationY
				});
				model.Loads.Add(new Load()
				{
					Amount = zValue.Value / forceNodes.Count,
					Node = model.NodesDictionary[node.Id],
					DOF = StructuralDof.TranslationZ
				});
			}
		}

		private static void CalculateVectorForce(IMechanicalExtAPI _api,IMechanicalUserSolver solver, Model model, Force force)
		{
			var forceLocation = _api.Application.InvokeUIThread(() => force.Location) as ISelectionInfo;
			var xValue = _api.Application.InvokeUIThread(() => force.XComponent.Output.DiscreteValues[0]) as Quantity;
			var yValue = _api.Application.InvokeUIThread(() => force.YComponent.Output.DiscreteValues[0]) as Quantity;
			var zValue = _api.Application.InvokeUIThread(() => force.ZComponent.Output.DiscreteValues[0]) as Quantity;
			var magnitude = _api.Application.InvokeUIThread(() => force.Magnitude.Output.DiscreteValues[1]) as Quantity;
			var forceSurfaceId = forceLocation.Ids[0];
			var forceNodes = solver.Analysis.MeshData.MeshRegionById(forceSurfaceId).Nodes;
            var a = solver.Analysis.MeshData.MeshRegionById(forceSurfaceId);

            a.Mesh.GetQuad4ExteriorFaces(out int[] quad4Elements);
            var element = solver.Analysis.MeshData.ElementById(quad4Elements[0]);
			
			foreach (var node in forceNodes)
			{
				model.Loads.Add(new Load()
				{
					Amount = xValue.Value * magnitude.Value / forceNodes.Count,
					Node = model.NodesDictionary[node.Id],
					DOF = StructuralDof.TranslationX
				});
				model.Loads.Add(new Load()
				{
					Amount = yValue.Value * magnitude.Value / forceNodes.Count,
					Node = model.NodesDictionary[node.Id],
					DOF = StructuralDof.TranslationY
				});
				model.Loads.Add(new Load()
				{
					Amount = zValue.Value * magnitude.Value / forceNodes.Count,
					Node = model.NodesDictionary[node.Id],
					DOF = StructuralDof.TranslationZ
				});
			}
		}


	}

    public class ElementFace
    {
		public int ElementId { get; set; }
		public int FaceId { get; set; }
    }
}
