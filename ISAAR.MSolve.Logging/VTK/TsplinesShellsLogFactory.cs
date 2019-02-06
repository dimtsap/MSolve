using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.IGA.Entities;
using ISAAR.MSolve.IGA.Postprocessing;
using ISAAR.MSolve.Logging.Interfaces;

namespace ISAAR.MSolve.Logging.VTK
{
	public class TsplinesShellsLogFactory : ILogFactory_v2
	{
		private readonly string _directory;
		private readonly Model _model;
		private readonly string _filename;

		public TsplinesShellsLogFactory(Model model, string filename)
		{
			_filename = filename;
			_model = model;
		}

		public IAnalyzerLog_v2[] CreateLogs()
		{
			var logs = new List<IAnalyzerLog_v2>(1);

			return new IAnalyzerLog_v2[]
			{
				new TsplineShellsLog(_model, _filename)
			};
		}
	}
}
