using ISAAR.MSolve.Analyzers.Interfaces;
using ISAAR.MSolve.IGA.Entities;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using IVector = ISAAR.MSolve.LinearAlgebra.Vectors.IVector;

namespace ISAAR.MSolve.Analyzers
{
    public class NonLinearPatchUpdater_v2 : INonLinearSubdomainUpdater_v2
    {
		private readonly ISubdomain_v2 patch;

		public NonLinearPatchUpdater_v2(ISubdomain_v2 patch)
		{
			this.patch = patch;
		}

		public void ScaleConstraints(double scalingFactor)
		{
			this.patch.ScaleConstraints(scalingFactor);
		}

		public IVector GetRhsFromSolution(IVectorView solution, IVectorView dSolution)
		{
			return patch.GetRhsFromSolution(solution, dSolution);
		}

		public void ResetState()
		{
			this.patch.ClearMaterialStresses();
		}

		public void UpdateState()
		{
			this.patch.SaveMaterialState();
		}
	}
}

