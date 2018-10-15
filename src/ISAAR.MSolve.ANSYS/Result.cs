using System;
using System.Collections.Generic;
using System.Text;
using Ansys.ACT.Interfaces.Mechanical;
using Ansys.ACT.Interfaces.Post;
using Ansys.ACT.Interfaces.UserObject;

namespace ISAAR.MSolve.ANSYS
{
	public class Result
	{
		internal double[] res = new double[1];

		public Result(IMechanicalExtAPI extApi, IUserResult result)
		{
			
		}

		public void evaluate(IUserResult entity, IStepInfo stepInfo, IResultCollector collector)
		{
			foreach (var id in collector.Ids)
			{
				res[0] = id;
				collector.SetValues(id, res);
			}
		} 
	}
}
