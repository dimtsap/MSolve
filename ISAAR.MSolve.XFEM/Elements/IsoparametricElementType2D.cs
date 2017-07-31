﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Integration.Quadratures;
using ISAAR.MSolve.XFEM.Interpolation;

namespace ISAAR.MSolve.XFEM.Elements
{
    abstract class IsoparametricElementType2D
    {
        public static readonly IsoparametricElementType2D Quad4 = new IsoparametricQuad4();
        public static readonly IsoparametricElementType2D Quad9 = new IsoparametricQuad9();

        private IsoparametricElementType2D() // prevents derived classes, except for nested ones
        { }

        /// <summary>
        /// TODO: Create a standard integration rule interface that guarantees gauss points that are
        /// i) immutable, ii) precached for fast generation, iii) stored globally for all elements
        /// TODO: Should this be stored as a field? It is a static property of the concrete (aka derived) class...
        /// </summary>
        public abstract IStandardQuadrature2D StandardQuadrature { get; }
        public abstract IsoparametricInterpolation2D Interpolation { get; }
        public abstract void CheckNodes(IReadOnlyList<Node2D> nodes);
        public abstract IReadOnlyList<INaturalPoint2D> NaturalCoordinatesOfNodes { get; }

        private class IsoparametricQuad4: IsoparametricElementType2D
        {
            public override IsoparametricInterpolation2D Interpolation
            {
                get { return IsoparametricInterpolation2D.Quad4; }
            }

            public override IStandardQuadrature2D StandardQuadrature
            {
                get { return GaussLegendre2D.Order2x2; }
            }

            public override void CheckNodes(IReadOnlyList<Node2D> nodes)
            {
                if (nodes.Count != 4) throw new ArgumentException("A Quad4 finite element has 4 nodes, but "
                + nodes.Count + " nodes were provided.");
                // TODO: Perhaps I can check the order of the nodes too or even the shape
            }

            public override IReadOnlyList<INaturalPoint2D> NaturalCoordinatesOfNodes
            {
                get
                {
                    return new INaturalPoint2D[]
                    {
                        new NaturalPoint2D(-1.0, -1.0),
                        new NaturalPoint2D(1.0, -1.0),
                        new NaturalPoint2D(1.0, 1.0),
                        new NaturalPoint2D(-1.0, 1.0)
                    };
                }
            }
        }

        private class IsoparametricQuad9 : IsoparametricElementType2D
        {
            public override IsoparametricInterpolation2D Interpolation
            {
                get { return IsoparametricInterpolation2D.Quad9; }
            }

            public override IStandardQuadrature2D StandardQuadrature
            {
                get { return GaussLegendre2D.Order3x3; }
            }

            public override void CheckNodes(IReadOnlyList<Node2D> nodes)
            {
                if (nodes.Count != 9) throw new ArgumentException("A Quad9 finite element has 9 nodes, but "
                + nodes.Count + " nodes were provided.");
                // TODO: Perhaps I can check the order of the nodes too or even the shape
            }

            public override IReadOnlyList<INaturalPoint2D> NaturalCoordinatesOfNodes
            {
                get
                {
                    return new INaturalPoint2D[]
                    {
                        new NaturalPoint2D(-1.0, -1.0),
                        new NaturalPoint2D(1.0, -1.0),
                        new NaturalPoint2D(1.0, 1.0),
                        new NaturalPoint2D(-1.0, 1.0),

                        new NaturalPoint2D(0.0, -1.0),
                        new NaturalPoint2D(1.0, 0.0),
                        new NaturalPoint2D(0.0, 1.0),
                        new NaturalPoint2D(-1.0, 0.0),
                        new NaturalPoint2D(0.0, 0.0),
                    };
                }
            }
        }
    }
}
