﻿using ISAAR.MSolve.Discretization.Mesh;
using ISAAR.MSolve.XFEM.Elements;
using System;
using System.Collections.Generic;
using System.Text;

//TODO: perhaps I should store the point ids directly, to save memory. Then I could just use LinkedLists for both points and cells
namespace ISAAR.MSolve.XFEM.Output.VTK
{
    class VtkCell2D
    {
        public static readonly IDictionary<CellType, int> CellTypeCodes =
            new Dictionary<CellType, int>()
            {
                { CellType.Quad4, 9 },
                { CellType.Tri3, 5 }
            };

        public int TypeCode { get; }
        public IReadOnlyList<VtkPoint2D> Vertices { get; } 

        public VtkCell2D(int cellType, IReadOnlyList<VtkPoint2D> vertices)
        {
            this.TypeCode = cellType;
            this.Vertices = vertices;
        }
    }
}
