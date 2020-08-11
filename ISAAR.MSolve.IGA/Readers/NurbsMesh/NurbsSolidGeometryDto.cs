using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.IGA.Readers.NurbsMesh
{
    public class NurbsSolidGeometryDto
    {
        public int ID { get; set; }

        public int NumberOfCpKsi { get; set; }

        public int NumberOfCpHeta { get; set; }

        public int NumberOfCpZeta { get; set; }

        public int DegreeKsi { get; set; }

        public int DegreeHeta { get; set; }

        public int DegreeZeta { get; set; }

        public double[] KnotValueVectorKsi { get; set; }

        public double[] KnotValueVectorHeta { get; set; }

        public double[] KnotValueVectorZeta { get; set; }

        public int[] ControlPointIDs { get; set; }
    }
}
