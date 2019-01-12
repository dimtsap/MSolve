using ISAAR.MSolve.Discretization.Interfaces;

namespace ISAAR.MSolve.Discretization.Loads
{
    public class ElementMassAccelerationHistoryLoad
    {
        public IElement Element { get; set; }
        public MassAccelerationHistoryLoad HistoryLoad { get; set; }
    }
}
