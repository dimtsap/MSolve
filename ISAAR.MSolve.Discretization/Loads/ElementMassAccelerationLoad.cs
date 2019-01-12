using ISAAR.MSolve.Discretization.Interfaces;

namespace ISAAR.MSolve.Discretization.Loads
{
    public class ElementMassAccelerationLoad
    {
        public IElement Element { get; set; }
        public DOFType DOF { get; set; }
        public double Amount { get; set; }
    }
}
