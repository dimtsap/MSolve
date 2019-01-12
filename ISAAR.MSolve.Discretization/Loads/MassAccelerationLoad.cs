using ISAAR.MSolve.Discretization.Interfaces;

namespace ISAAR.MSolve.Discretization.Loads
{
    public class MassAccelerationLoad
    {
        public DOFType DOF { get; set; }
        public double Amount { get; set; }
    }
}
