﻿using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Discretization.Loads;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.IGA.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Numerical.Commons;

//TODO: find what is going on with the dynamic loads and refactor them. That 564000000 in AssignMassAccelerationHistoryLoads()
//      cannot be correct.
namespace ISAAR.MSolve.IGA.Entities
{
	public class Model : IStructuralModel_v2
	{
		private int numberOfPatches = 0;
		private int numberOfInterfaces = 0;

		private readonly Dictionary<int, ControlPoint> controlPointsDictionary = new Dictionary<int, ControlPoint>();
		private readonly Dictionary<int, Element> elementsDictionary = new Dictionary<int, Element>();
		private readonly Dictionary<int, Patch> patchesDictionary = new Dictionary<int, Patch>();

		private readonly IList<Load> loads = new List<Load>();

		private readonly IList<IMassAccelerationHistoryLoad> massAccelerationHistoryLoads =
			new List<IMassAccelerationHistoryLoad>();
		public IList<MassAccelerationLoad> MassAccelerationLoads { get; } = new List<MassAccelerationLoad>();

		private IGlobalFreeDofOrdering globalDofOrdering;

		//public IList<EmbeddedNode> EmbeddedNodes { get; } = new List<EmbeddedNode>();

		public IList<Cluster> Clusters => ClustersDictionary.Values.ToList();
		public Dictionary<int, Cluster> ClustersDictionary { get; } = new Dictionary<int, Cluster>();

		IReadOnlyList<IElement> IStructuralModel_v2.Elements => ElementsDictionary.Values.ToList();
		public IList<Element> Elements => ElementsDictionary.Values.ToList();
		public Dictionary<int, Element> ElementsDictionary => elementsDictionary;

		public IList<Load> Loads { get; } = new List<Load>();

		public IList<IMassAccelerationHistoryLoad> MassAccelerationHistoryLoads { get; } =
			new List<IMassAccelerationHistoryLoad>();

		public IList<ElementMassAccelerationHistoryLoad> ElementMassAccelerationHistoryLoads { get; } =
			new List<ElementMassAccelerationHistoryLoad>();

		public IList<ControlPoint> ControlPoints => controlPointsDictionary.Values.ToList();
		IReadOnlyList<INode> IStructuralModel_v2.Nodes => controlPointsDictionary.Values.ToList();

		public Dictionary<int, ControlPoint> ControlPointsDictionary
		{
			get => controlPointsDictionary;
		}

		public int NumberOfPatches
		{
			get { return numberOfPatches; }
			set { numberOfPatches = value; }
		}

		public int NumberOfInterfaces
		{
			get { return numberOfInterfaces; }
			set { numberOfInterfaces = value; }
		}

		IReadOnlyList<ISubdomain_v2> IStructuralModel_v2.Subdomains => patchesDictionary.Values.ToList();
		public IList<Patch> Patches => patchesDictionary.Values.ToList();

		public Dictionary<int, Patch> PatchesDictionary
		{
			get => patchesDictionary;
		}

		public Table<INode, DOFType, double> Constraints { get; private set; } =
			new Table<INode, DOFType, double>(); //TODOMaria: maybe it's useless in model class

		public IGlobalFreeDofOrdering GlobalDofOrdering
		{
			get => globalDofOrdering;
			set
			{
				globalDofOrdering = value;
				foreach (Patch patch in Patches)
				{
					patch.DofOrdering = GlobalDofOrdering.SubdomainDofOrderings[patch];
					patch.Forces = Vector.CreateZero(patch.DofOrdering.NumFreeDofs);
				}

				//EnumerateSubdomainLagranges();
				//EnumerateDOFMultiplicity();
			}
		}

		public void AssignLoads()
		{
			foreach (Patch patch in PatchesDictionary.Values) patch.Forces.Clear();
			AssignControlPointLoads();
			AssignBoundaryLoads();
			AssignMassAccelerationLoads();
		}

		private void AssignMassAccelerationLoads()
		{
			if (MassAccelerationLoads.Count<1) return;

			foreach (var patch in PatchesDictionary.Values)
			{
				foreach (var element in patch.Elements)
				{
					patch.DofOrdering.AddVectorElementToSubdomain(element,
						Vector.CreateFromArray(element.ElementType.CalculateAccelerationForces(element,MassAccelerationLoads)),
						patch.Forces);
				}
			}
		}

		public void AssignMassAccelerationHistoryLoads(int timeStep)
		{
			if (MassAccelerationHistoryLoads.Count > 0)
			{
				List<MassAccelerationLoad> m=new List<MassAccelerationLoad>(MassAccelerationHistoryLoads.Count);
				foreach (var load in MassAccelerationHistoryLoads)
				{
					m.Add(new MassAccelerationLoad() {Amount = load[timeStep], DOF = load.DOF});
				}

				foreach (var patch in PatchesDictionary.Values)
				{
					foreach (var element in patch.Elements)
					{
							patch.DofOrdering.AddVectorElementToSubdomain(element,
								Vector.CreateFromArray(element.ElementType.CalculateAccelerationForces(element,m)),
								patch.Forces);
					}
				}
			}

			foreach (var load in ElementMassAccelerationHistoryLoads)
			{
				MassAccelerationLoad hl = new MassAccelerationLoad() { Amount = load.HistoryLoad[timeStep] * 564000000, DOF = load.HistoryLoad.DOF };
				Element element = load.Element as Element;
				ISubdomain_v2 subdomain = element.Patch;
				var accelerationForces = Vector.CreateFromArray(element.ElementType.CalculateAccelerationForces(
					((Element)load.Element), (new MassAccelerationLoad[] { hl }).ToList()));
				globalDofOrdering.SubdomainDofOrderings[subdomain].AddVectorElementToSubdomain(element, accelerationForces,
					subdomain.Forces);
			}
		}

		private void AssignBoundaryLoads()
		{
			foreach (Patch patch in patchesDictionary.Values)
			{
				foreach (Edge edge in patch.EdgesDictionary.Values)
				{
					Dictionary<int, double> edgeLoadDictionary = edge.CalculateLoads();
					foreach (int dof in edgeLoadDictionary.Keys)
						if (dof != -1)
							patch.Forces[dof] += edgeLoadDictionary[dof];
				}

				foreach (Face face in patch.FacesDictionary.Values)
				{
					Dictionary<int, double> faceLoadDictionary = face.CalculateLoads();
					foreach (int dof in faceLoadDictionary.Keys)
						if (dof != -1)
							patch.Forces[dof] += faceLoadDictionary[dof];
				}
			}
		}

		private void AssignControlPointLoads()
		{
			foreach (Patch patch in PatchesDictionary.Values)
			{
				patch.ControlPointLoads = new Table<ControlPoint, DOFType, double>();
			}

			foreach (Load load in Loads)
			{
				var cp = ((ControlPoint) load.ControlPoint);
				double amountPerPatch = load.Amount / cp.PatchesDictionary.Count;
				foreach (Patch patch in cp.PatchesDictionary.Values)
				{
					bool wasNotContained = patch.ControlPointLoads.TryAdd(cp, load.DOF, amountPerPatch);
				}
			}

			//TODO: this should be done by the subdomain when the analyzer decides.
			foreach (Patch patch in PatchesDictionary.Values)
			{
				foreach ((ControlPoint node, DOFType dofType, double amount) in patch.ControlPointLoads)
				{
					if (!patch.DofOrdering.FreeDofs.Contains(node, dofType)) continue;
					int patchDofIdx = patch.DofOrdering.FreeDofs[node, dofType];
					patch.Forces[patchDofIdx] = amount;
				}
			}
		}


		//What is the purpose of this method? If someone wanted to clear the Model, they could just create a new one.
		public void Clear()
		{
			Loads.Clear();
			ClustersDictionary.Clear();
			PatchesDictionary.Clear();
			ElementsDictionary.Clear();
			ControlPointsDictionary.Clear();
			globalDofOrdering = null;
			Constraints.Clear();
			MassAccelerationHistoryLoads.Clear();
		}

		public void ConnectDataStructures()
		{
			BuildInterconnectionData();
			AssignConstraints();
		}

		//TODO: constraints should not be saved inside the nodes. As it is right now (22/11/2018) the same constraint 
		//      is saved in the node, the model constraints table and the subdomain constraints table. Furthermore,
		//      displacement control analyzer updates the subdomain constraints table only (another bad design decision).  
		//      It is too easy to access the wrong instance of the constraint. 
		private void AssignConstraints()
		{
			foreach (ControlPoint controlPoint in ControlPointsDictionary.Values)
			{
				if (controlPoint.Constraints == null) continue;
				foreach (Constraint constraint in controlPoint.Constraints)
					Constraints[controlPoint, constraint.DOF] = constraint.Amount;
			}

			foreach (Patch patch in PatchesDictionary.Values) patch.ExtractConstraintsFromGlobal(Constraints);
		}

		private void BuildElementDictionaryOfEachControlPoint()
		{
			foreach (Element element in elementsDictionary.Values)
			foreach (ControlPoint controlPoint in element.ControlPoints)
				controlPoint.ElementsDictionary.Add(element.ID, element);
		}

		private void
			BuildInterconnectionData() //TODOMaria: maybe I have to generate the constraints dictionary for each subdomain here
		{
			BuildPatchOfEachElement();
			BuildElementDictionaryOfEachControlPoint();
			foreach (ControlPoint controlPoint in ControlPointsDictionary.Values) controlPoint.BuildPatchesDictionary();

			//foreach (Patch patch in PatchesDictionary.Values) patch.DefineControlPointsFromElements();
		}

		private void BuildPatchOfEachElement()
		{
			foreach (Patch patch in patchesDictionary.Values)
			foreach (Element element in patch.Elements)
			{
				element.Patch = patch;
				element.Model = this;
			}
		}
	}
}