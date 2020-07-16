using System.IO;
using System.Linq;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.IGA.Elements;
using ISAAR.MSolve.IGA.Entities;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.IGA.Postprocessing
{
    /// <summary>
	/// Paraview file  generator for NURBS Shells geometries.
	/// </summary>
	public class ParaviewNurbsShells
	{
		private readonly string _filename;
		private readonly Model _model;
		private readonly IVectorView _solution;

		/// <summary>
		/// Defines a Paraview FileWriter.
		/// </summary>
		/// <param name="model">An isogeometric <see cref="Model"/>.</param>
		/// <param name="solution">An <see cref="IVectorView"/> containing the solution of the linear system.</param>
		/// <param name="filename">The name of the paraview file to be generated.</param>
		public ParaviewNurbsShells(Model model, IVectorView solution, string filename)
		{
			_model = model;
			_solution = solution;
			_filename = filename;
		}

        #region Unstructured
        public void CreateParaviewFile()
        {
            var projectiveControlPoints = CalculateProjectiveControlPoints();
            var numberOfPointsPerElement = 4;
            var nodes = new double[_model.Elements.Count * numberOfPointsPerElement, 3];
            var pointIndex = 0;
            foreach (var element in _model.Elements.Where(x => !(x.ElementType is PenaltyDofPair)))
            {
                var tsplineElement = element as NurbsKirchhoffLoveShellElementNL;
                var elementPoints = tsplineElement.CalculatePointsForPostProcessing(tsplineElement);

                for (int i = 0; i < elementPoints.GetLength(0); i++)
                {
                    nodes[pointIndex, 0] = elementPoints[i, 0];
                    nodes[pointIndex, 1] = elementPoints[i, 1];
                    nodes[pointIndex++, 2] = elementPoints[i, 2];
                }
            }

            var elementConnectivity = CreateTsplineConnectivity();

            var pointDisplacements = new double[nodes.GetLength(0), 3];
            foreach (var element in _model.Elements.Where(x => !(x.ElementType is PenaltyDofPair)))
            {
                var localDisplacements = Matrix.CreateZero(element.ControlPointsDictionary.Count, 3);
                var counterCP = 0;
                foreach (var controlPoint in element.ControlPoints)
                {
                    localDisplacements[counterCP, 0] =
                        (!_model.GlobalDofOrdering.GlobalFreeDofs.Contains(controlPoint, StructuralDof.TranslationX))
                            ? 0.0
                            : _solution[_model.GlobalDofOrdering.GlobalFreeDofs[controlPoint, StructuralDof.TranslationX]];
                    localDisplacements[counterCP, 1] =
                        (!_model.GlobalDofOrdering.GlobalFreeDofs.Contains(controlPoint, StructuralDof.TranslationY))
                            ? 0.0
                            : _solution[_model.GlobalDofOrdering.GlobalFreeDofs[controlPoint, StructuralDof.TranslationY]];
                    localDisplacements[counterCP++, 2] =
                        (!_model.GlobalDofOrdering.GlobalFreeDofs.Contains(controlPoint, StructuralDof.TranslationZ))
                            ? 0.0
                            : _solution[_model.GlobalDofOrdering.GlobalFreeDofs[controlPoint, StructuralDof.TranslationZ]];
                }
                var elementKnotDisplacements = element.ElementType.CalculateDisplacementsForPostProcessing(element, localDisplacements);
                for (int i = 0; i < elementConnectivity.GetLength(1); i++)
                {
                    var knotConnectivity = elementConnectivity[element.ID, i];
                    pointDisplacements[knotConnectivity, 0] = elementKnotDisplacements[i, 0];
                    pointDisplacements[knotConnectivity, 1] = elementKnotDisplacements[i, 1];
                    pointDisplacements[knotConnectivity, 2] = elementKnotDisplacements[i, 2];
                }
            }

            WriteTSplineShellsFile(nodes, elementConnectivity, pointDisplacements);
        }

        private void WriteTSplineShellsFile(double[,] nodeCoordinates, int[,] elementConnectivity, double[,] displacements)
		{
			var numberOfPoints = nodeCoordinates.GetLength(0);
			var numberOfCells = elementConnectivity.GetLength(0);

			const int numberOfVerticesPerCell = 4;
			const int paraviewCellCode = 9;

			using (StreamWriter outputFile = new StreamWriter(Path.Combine(Directory.GetCurrentDirectory(),$"{_filename}Paraview.vtu")))
			{
				outputFile.WriteLine("<?xml version=\"1.0\"?>");
				outputFile.WriteLine("<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">");
				outputFile.WriteLine("<UnstructuredGrid>");

				outputFile.WriteLine($"<Piece NumberOfPoints=\"{numberOfPoints}\" NumberOfCells=\"{numberOfCells}\">");
				outputFile.WriteLine($"<PointData Vectors=\"U\">");
				outputFile.WriteLine($"<DataArray type=\"Float32\" Name=\"U\" format=\"ascii\" NumberOfComponents=\"3\">");

				for (int i = 0; i < numberOfPoints; i++)
					outputFile.WriteLine($"{displacements[i, 0]} {displacements[i, 1]} {displacements[i, 2]}");

				outputFile.WriteLine("</DataArray>");
				outputFile.WriteLine("</PointData>");
				outputFile.WriteLine("<Points>");
				outputFile.WriteLine("<DataArray type=\"Float32\" NumberOfComponents=\"3\">");

				for (int i = 0; i < numberOfPoints; i++)
					outputFile.WriteLine($"{nodeCoordinates[i, 0]} {nodeCoordinates[i, 1]} {nodeCoordinates[i, 2]}");

				outputFile.WriteLine("</DataArray>");
				outputFile.WriteLine("</Points>");
				outputFile.WriteLine("<Cells>");
				outputFile.WriteLine("<DataArray type=\"Int32\" Name=\"connectivity\">");

				for (int i = 0; i < numberOfCells; i++)
				{
					for (int j = 0; j < elementConnectivity.GetLength(1); j++)
						outputFile.Write($"{elementConnectivity[i, j]} ");
					outputFile.WriteLine("");
				}

				outputFile.WriteLine("</DataArray>");
				outputFile.WriteLine("<DataArray type=\"Int32\" Name=\"offsets\">");

				var offset = 0;
				for (int i = 0; i < numberOfCells; i++)
				{
					offset += numberOfVerticesPerCell;
					outputFile.WriteLine(offset);
				}

				outputFile.WriteLine("</DataArray>");
				outputFile.WriteLine("<DataArray type=\"Int32\" Name=\"types\">");

				for (int i = 0; i < numberOfCells; i++)
					outputFile.WriteLine(paraviewCellCode);

				outputFile.WriteLine("</DataArray>");
				outputFile.WriteLine("</Cells>");
				outputFile.WriteLine("</Piece>");
				outputFile.WriteLine("</UnstructuredGrid>");
				outputFile.WriteLine("</VTKFile>");
			}
		}

        private double[,] CalculateProjectiveControlPoints()
        {
            var projectiveCPs = new double[_model.PatchesDictionary[0].ControlPoints.Count, 4];
            foreach (var controlPoint in _model.PatchesDictionary[0].ControlPoints)
            {
                var weight = controlPoint.WeightFactor;
                projectiveCPs[controlPoint.ID, 0] = controlPoint.X * weight;
                projectiveCPs[controlPoint.ID, 1] = controlPoint.Y * weight;
                projectiveCPs[controlPoint.ID, 2] = controlPoint.Z * weight;
                projectiveCPs[controlPoint.ID, 3] = weight;
            }

            return projectiveCPs;
        }

        private int[,] CreateTsplineConnectivity()
		{
			var elementConnectivity = new int[_model.Elements.Count, 4];
			var pointCounter = 0;
			for (int i = 0; i < _model.Elements.Count; i++)
			{
				elementConnectivity[i, 0] = pointCounter++;
				elementConnectivity[i, 1] = pointCounter++;
				elementConnectivity[i, 2] = pointCounter++;
				elementConnectivity[i, 3] = pointCounter++;
			}

			return elementConnectivity;
		}



        #endregion
        ///// <summary>
        ///// Creates Paraview File of the NURBS Shells geometry.
        ///// </summary>
        //public void CreateParaviewFile()
        //{
        //    var uniqueKnotsKsi = _model.PatchesDictionary[0].KnotValueVectorKsi.RemoveDuplicatesFindMultiplicity();
        //    var uniqueKnotsHeta = _model.PatchesDictionary[0].KnotValueVectorHeta.RemoveDuplicatesFindMultiplicity();

        //    var numberOfKnotsKsi = uniqueKnotsKsi[0].Length;
        //    var numberOfKnotsHeta = uniqueKnotsHeta[0].Length;

        //    var knots = new double[numberOfKnotsKsi * numberOfKnotsHeta, 3];
        //    var count = 0;
        //    var patch = _model.PatchesDictionary[0];

        //    var projectiveControlPoints = CalculateProjectiveControlPoints();

        //    for (var knotKsiIndex = 0; knotKsiIndex < numberOfKnotsKsi; knotKsiIndex++)
        //    {
        //        for (var knotHetaIndex = 0; knotHetaIndex < numberOfKnotsHeta; knotHetaIndex++)
        //        {
        //            var hetaCoordinate = uniqueKnotsHeta[0][knotHetaIndex];
        //            var ksiCoordinate = uniqueKnotsKsi[0][knotKsiIndex];
        //            var point3D = ParaviewNurbs2D.SurfacePoint2D(patch.NumberOfControlPointsKsi - 1, patch.DegreeKsi,
        //                patch.KnotValueVectorKsi, patch.NumberOfControlPointsHeta - 1, patch.DegreeHeta,
        //                patch.KnotValueVectorHeta, projectiveControlPoints, ksiCoordinate, hetaCoordinate); //TODO: Find what goes wrong here
        //            knots[count, 0] = point3D[0] / point3D[3];
        //            knots[count, 1] = point3D[1] / point3D[3];
        //            knots[count, 2] = point3D[2] / point3D[3];
        //            count++;
        //        }
        //    }

        //    var incrementKsi = numberOfKnotsHeta;
        //    var incrementHeta = 1;
        //    var nodePattern = new int[] { 0, incrementKsi, incrementKsi + 1, 1 };
        //    var elementConnectivity = CreateElement2DConnectivity(nodePattern, uniqueKnotsKsi[0].Length - 1,
        //        uniqueKnotsHeta[0].Length - 1, incrementKsi, incrementHeta);
        //    var knotDisplacements = new double[knots.GetLength(0), 3];

        //    foreach (var element in _model.Elements.Where(x => !(x.ElementType is PenaltyDofPair)))
        //    {
        //        var localDisplacements = Matrix.CreateZero(element.ControlPointsDictionary.Count, 3);
        //        var counterCP = 0;
        //        foreach (var controlPoint in element.ControlPoints)
        //        {
        //            localDisplacements[counterCP, 0] =
        //                (!_model.GlobalDofOrdering.GlobalFreeDofs.Contains(controlPoint, StructuralDof.TranslationX))
        //                    ? 0.0
        //                    : _solution[_model.GlobalDofOrdering.GlobalFreeDofs[controlPoint, StructuralDof.TranslationX]];
        //            localDisplacements[counterCP, 1] =
        //                (!_model.GlobalDofOrdering.GlobalFreeDofs.Contains(controlPoint, StructuralDof.TranslationY))
        //                    ? 0.0
        //                    : _solution[_model.GlobalDofOrdering.GlobalFreeDofs[controlPoint, StructuralDof.TranslationY]];
        //            localDisplacements[counterCP, 2] =
        //                (!_model.GlobalDofOrdering.GlobalFreeDofs.Contains(controlPoint, StructuralDof.TranslationZ))
        //                    ? 0.0
        //                    : _solution[_model.GlobalDofOrdering.GlobalFreeDofs[controlPoint, StructuralDof.TranslationZ]];
        //            counterCP++;

        //        }
        //        var elementKnotDisplacements = element.ElementType.CalculateDisplacementsForPostProcessing(element, localDisplacements);
        //        for (int i = 0; i < elementConnectivity.GetLength(1); i++)
        //        {
        //            var knotConnectivity = elementConnectivity[element.ID, i];
        //            knotDisplacements[knotConnectivity, 0] = elementKnotDisplacements[i, 0];
        //            knotDisplacements[knotConnectivity, 1] = elementKnotDisplacements[i, 1];
        //            knotDisplacements[knotConnectivity, 2] = elementKnotDisplacements[i, 2];
        //        }
        //    }

        //    Write2DNurbsFile(knots, elementConnectivity, "Quad4", knotDisplacements);
        //}

        //private void Write2DNurbsFile(double[,] nodeCoordinates, int[,] elementConnectivity, string elementType, double[,] displacements)
        //{
        //    var numberOfNodes = nodeCoordinates.GetLength(0);
        //    var numberOfCells = elementConnectivity.GetLength(0);

        //    int numberOfVerticesPerCell = 0;
        //    int paraviewCellCode = 0;

        //    if (elementType == "Quad4")
        //    {
        //        numberOfVerticesPerCell = 4;
        //        paraviewCellCode = 9;
        //    }

        //    using (StreamWriter outputFile = new StreamWriter(Path.Combine(Directory.GetCurrentDirectory(), $"{_filename}Paraview.vtu")))
        //    {
        //        outputFile.WriteLine("<VTKFile type=\"UnstructuredGrid\"  version=\"0.1\"   >");
        //        outputFile.WriteLine("<UnstructuredGrid>");
        //        outputFile.WriteLine($"<Piece  NumberOfPoints=\"{numberOfNodes}\" NumberOfCells=\"{numberOfCells}\">");

        //        outputFile.WriteLine("<Points>");
        //        outputFile.WriteLine("<DataArray  type=\"Float64\"  NumberOfComponents=\"3\"  format=\"ascii\" >");
        //        for (int i = 0; i < numberOfNodes; i++)
        //            outputFile.WriteLine($"{nodeCoordinates[i, 0]} {nodeCoordinates[i, 1]} {nodeCoordinates[i, 2]}");

        //        outputFile.WriteLine("</DataArray>");
        //        outputFile.WriteLine("</Points>");

        //        outputFile.WriteLine("<Cells>");
        //        outputFile.WriteLine("<DataArray  type=\"Int32\"  Name=\"connectivity\"  format=\"ascii\">");
        //        for (int i = 0; i < numberOfCells; i++)
        //        {
        //            for (int j = 0; j < elementConnectivity.GetLength(1); j++)
        //                outputFile.Write($"{elementConnectivity[i, j]} ");
        //            outputFile.WriteLine("");
        //        }

        //        outputFile.WriteLine("</DataArray>");
        //        outputFile.WriteLine("<DataArray  type=\"Int32\"  Name=\"offsets\"  format=\"ascii\">");

        //        var offset = 0;
        //        for (int i = 0; i < numberOfCells; i++)
        //        {
        //            offset += numberOfVerticesPerCell;
        //            outputFile.WriteLine(offset);
        //        }

        //        outputFile.WriteLine("</DataArray>");
        //        outputFile.WriteLine("<DataArray  type=\"UInt8\"  Name=\"types\"  format=\"ascii\">");
        //        for (int i = 0; i < numberOfCells; i++)
        //            outputFile.WriteLine(paraviewCellCode);

        //        outputFile.WriteLine("</DataArray>");
        //        outputFile.WriteLine("</Cells>");

        //        outputFile.WriteLine("<PointData  Vectors=\"U\">");

        //        outputFile.WriteLine("<DataArray  type=\"Float64\"  Name=\"U\" NumberOfComponents=\"3\" format=\"ascii\">");

        //        for (int i = 0; i < numberOfNodes; i++)
        //            outputFile.WriteLine($"{displacements[i, 0]} {displacements[i, 1]} {displacements[i, 2]}");

        //        outputFile.WriteLine("</DataArray>");
        //        outputFile.WriteLine("</PointData>");
        //        outputFile.WriteLine("</Piece>");
        //        outputFile.WriteLine("</UnstructuredGrid>");
        //        outputFile.WriteLine("</VTKFile>");
        //    }
        //}

        ///// <summary>
        ///// Creates a control point coordinates matrix in projective coordinates
        ///// </summary>
        ///// <returns></returns>
        //private double[,] CalculateProjectiveControlPoints()
        //{
        //    var projectiveCPs = new double[_model.PatchesDictionary[0].ControlPoints.Count, 4];
        //    foreach (var controlPoint in _model.PatchesDictionary[0].ControlPoints)
        //    {
        //        var weight = controlPoint.WeightFactor;
        //        projectiveCPs[controlPoint.ID, 0] = controlPoint.X * weight;
        //        projectiveCPs[controlPoint.ID, 1] = controlPoint.Y * weight;
        //        projectiveCPs[controlPoint.ID, 2] = controlPoint.Z * weight;
        //        projectiveCPs[controlPoint.ID, 3] = weight;
        //    }

        //    return projectiveCPs;
        //}

        //private int[,] CreateElement2DConnectivity(int[] nodePattern, int numberOfElementsKsi, int numberOfElementsHeta,
        //            int incrementKsi, int incrementHeta)
        //{
        //    var increment = 0;
        //    var elementConnectivity = new int[numberOfElementsKsi * numberOfElementsHeta, nodePattern.Length];
        //    var elementCounter = 0;
        //    for (int elementKsi = 0; elementKsi < numberOfElementsKsi; elementKsi++)
        //    {
        //        increment = elementKsi * incrementKsi;
        //        for (int elementHeta = 0; elementHeta < numberOfElementsHeta; elementHeta++)
        //        {
        //            for (int i = 0; i < nodePattern.Length; i++)
        //                elementConnectivity[elementCounter, i] = nodePattern[i] + increment;
        //            increment += incrementHeta;
        //            elementCounter++;
        //        }
        //    }

        //    return elementConnectivity;
        //}
    }
}
