using System;
using System.Collections.Generic;
using System.Text;
using Ansys.ACT.Interfaces.Post;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using Newtonsoft.Json;

namespace MSolveANSYS
{
	public class MSolveStaticReader : ICustomResultReader
	{
		private string[] _infos;
		private Dictionary<int, Dictionary<DOFType, int>> NodalDOFsDictionary;
		private SolutionVector SolutionVector;

		public MSolveStaticReader(string[] infos)
		{
			_infos = infos;
			NodalDOFsDictionary = JsonConvert.DeserializeObject<Dictionary<int, Dictionary<DOFType, int>>>(_infos[0]);
			SolutionVector = JsonConvert.DeserializeObject<SolutionVector>(_infos[1]);
		}

		public IEnumerable<string> GetComponentNames(string resultName)
		{
			return new[] { "X", "Y", "Z" };
		}

		public string GetComponentUnit(string resultName, string componentName)
		{
			return "Length";
		}

		public ResultLocationEnum GetResultLocation(string resultName)
		{
			return ResultLocationEnum.Node;
		}

		public IEnumerable<string> GetResultNames()
		{
			return new[] { "U" };
		}

		public ResultTypeEnum GetResultType(string resultName)
		{
			return ResultTypeEnum.Vector;
		}

		public IEnumerable<double> GetStepValues()
		{
			return new[] { 1.0 };
		}

		public void GetValues(string resultName, IResultCollector collector)
		{
			if (resultName != "U") return;
			foreach (var id in collector.Ids)
			{
				var msolveId = id;
				collector.Values[3* (msolveId-1)] = (NodalDOFsDictionary[msolveId][DOFType.X] == -1) ? 0
					: SolutionVector.Data[NodalDOFsDictionary[msolveId][DOFType.X]];
				collector.Values[3 * (msolveId - 1) + 1] = (NodalDOFsDictionary[msolveId][DOFType.Y] == -1) ? 0
					: SolutionVector.Data[NodalDOFsDictionary[msolveId][DOFType.Y]];
				collector.Values[3 * (msolveId - 1) + 2] = (NodalDOFsDictionary[msolveId][DOFType.Z] == -1) ? 0
					: SolutionVector.Data[NodalDOFsDictionary[msolveId][DOFType.Z]];
			}
		}

		public void SetCurrentStep(IStepInfo stepInfo)
		{
		}
	}
}
