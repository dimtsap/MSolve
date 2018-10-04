using System;
using System.Collections.Generic;
using System.Data;
using System.IO;
using System.Text;
using Ansys.ACT.Interfaces.Common;
using Ansys.ACT.Interfaces.Mechanical;
using Ansys.ACT.Interfaces.Mesh;
using Ansys.ACT.Interfaces.UserObject;
using SimpleExpressionEvaluator;

namespace ISAAR.MSolve.ANSYS
{
	public class FaceLoad
	{
		private readonly IMechanicalExtAPI _api;
		private readonly IMechanicalUserLoad _load;

		public FaceLoad(IExtAPI api, IUserLoad load)
		{
			_api = (IMechanicalExtAPI) api;
			_load = (IMechanicalUserLoad) load;
		}

		
	}
}
