using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.IGA.Entities;
using ISAAR.MSolve.IGA.Postprocessing;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Logging.Interfaces;

namespace ISAAR.MSolve.Logging.VTK
{
	public class TsplineShellsLog : IAnalyzerLog_v2
	{
		private static int iteration=0;
		private Model _model;
		private string _filename;
		public TsplineShellsLog(Model model, string filename)
		{
			_model = model;
			_filename = filename;
		}

		public void StoreResults(DateTime startTime, DateTime endTime, IVectorView solution)
		{
			var paraview = new ParaviewTsplineShells(_model);
			paraview.CreateParaviewFile(solution, $"{_filename}_{iteration}");
			iteration++;
		}
	}
}
