using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.IGA.Readers.NurbsMesh
{
    public class ModelDto
    {
        public List<NurbsSurfaceGeometryDto> NurbsSurfacePatches { get; set; }
        public List<NurbsSolidGeometryDto> NurbsSolidPatches { get; set; }
        public List<ControlPointDto> ControlPoints { get; set; }
    }
}
