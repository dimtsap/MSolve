﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Discretization.Mesh;
using ISAAR.MSolve.XFEM.Elements;

namespace ISAAR.MSolve.XFEM.Output.VTK
{
    class VTKCell
    {
        public static readonly IDictionary<CellType, int> CellTypeCodes =
            new Dictionary<CellType, int>()
            {
                { CellType.Quad4, 9 }
            };

        public XContinuumElement2D OriginalElement { get; }
        public int TypeCode { get; }
        public IReadOnlyList<VtkPoint2D> Vertices { get; }

        public VTKCell(XContinuumElement2D originalElement, int cellType, IReadOnlyList<VtkPoint2D> vertices)
        {
            this.OriginalElement = originalElement;
            this.TypeCode = cellType;
            this.Vertices = vertices;
        }
    }
}
