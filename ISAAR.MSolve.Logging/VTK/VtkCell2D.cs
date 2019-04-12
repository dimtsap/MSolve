﻿using System.Collections.Generic;
using ISAAR.MSolve.Geometry.Shapes;

namespace ISAAR.MSolve.Logging.VTK
{
    /// <summary>
    /// Cell used to represent VTK grids.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class VtkCell2D
    {
        public static readonly IReadOnlyDictionary<CellType, int> cellTypeCodes = 
            new Dictionary<CellType, int>
            {
                                            // 3 ---- 2
                                            // |      |
                                            // |      |
                { CellType.Quad4, 9 },    // 0 ---- 1

                                            // 3 -- 6 -- 2
                                            // |         |
                                            // 7         5
                                            // |         |
                { CellType.Quad8, 23 },   // 0 -- 4 -- 1

                                            // 3 -- 6 -- 2
                                            // |    |    |
                                            // 7 -- 8 -- 5
                                            // |    |    |
                { CellType.Quad9, 28 },   // 0 -- 4 -- 1

                                            //    2
                                            //   /  \
                                            //  /    \
                { CellType.Tri3, 5 },     // 0 ---  1

                                            //     2
                                            //    /  \
                                            //   5    4
                                            //  /       \
                { CellType.Tri6, 22 }     // 0 -- 3 -- 1
            };

        public VtkCell2D(int code, IReadOnlyList<VtkPoint2D> vertices)
        {
            this.Code = code;
            this.Vertices = vertices;
        }

        public int Code { get; }
        public IReadOnlyList<VtkPoint2D> Vertices { get; }
    }
}
