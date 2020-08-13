using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.IGA.Readers.NurbsMesh
{
    public class ControlPointDto
    {
        public int ID { get; set; }
        public double X { get; set; }
        public double Y { get; set; }
        public double Z { get; set; }
        public double Weight { get; set; }
    }
}
