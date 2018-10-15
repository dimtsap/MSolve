using System;
using System.Collections.Generic;
using System.Text;
using Ansys.ACT.Interfaces.Mechanical;
using Ansys.ACT.Interfaces;
using Ansys.ACT.Interfaces.UserObject;
using Ansys.ACT.Interfaces.Common;
using Ansys.ACT.Interfaces.Mesh;

namespace ISAAR.MSolve.ANSYS
{
	public class Load
	{
		private readonly IMechanicalExtAPI _api;
		private readonly IMechanicalUserLoad _load;

		public Load(IExtAPI api, IUserLoad load)
		{
			_api = (IMechanicalExtAPI)api;
			_load = (IMechanicalUserLoad) load;
		}

		public virtual IEnumerable<double> getnodalvaluesfordisplay(IUserLoad load, IEnumerable<int> nodeIds)
		{
			var res = new List<double>();
			IMeshData mesh = _load.Analysis.MeshData;
			foreach (var nodeId in nodeIds)
			{
				INode node = mesh.NodeById(nodeId);
				res.Add(Math.Sqrt(node.X * node.X + node.Y * node.Y + node.Z * node.Z));
			}

			return res;
		}
	}
}
