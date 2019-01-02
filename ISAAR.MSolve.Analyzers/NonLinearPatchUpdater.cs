using ISAAR.MSolve.Analyzers.Interfaces;
using ISAAR.MSolve.IGA.Entities;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Text;

namespace ISAAR.MSolve.Analyzers
{
    public class NonLinearPatchUpdater : INonLinearSubdomainUpdater
    {
        private readonly Patch patch;

        public NonLinearPatchUpdater(Patch patch)
        {
            this.patch = patch;
        }

        public void ScaleConstraints(double scalingFactor)
        {
	        throw new NotImplementedException();
	        //this.patch.ScaleConstraints(scalingFactor);
		}

        public IVector GetRHSFromSolution(IVector solution, IVector dSolution) //TODO leave 
        {
            return this.patch.GetRHSFromSolution(solution, dSolution);
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

