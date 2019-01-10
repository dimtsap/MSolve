using System;
using System.Collections.Generic;
using System.IO;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.IGA.Elements;
using ISAAR.MSolve.IGA.Entities;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Solvers.Interfaces;

namespace ISAAR.MSolve.IGA.Postprocessing
{
	public class ParaviewTsplines2D
	{
		private Model _model;
		private IVectorView _solution;
		private string _filename;

		public ParaviewTsplines2D(Model model, IVectorView solution, string filename)
		{
			_model = model;
			_solution = solution;
			_filename = filename;
		}

		public void CreateParaviewFile()
		{
			var projectiveControlPoints = CalculateProjectiveControlPoints();
			var numberOfPointsPerElement = 4;
			var nodes = new double[_model.Elements.Count * numberOfPointsPerElement, 2];
			var pointIndex = 0;
			foreach (var element in _model.Elements)
			{
				var tsplineElement = element as TSplineElement2D;
				var elementPoints = tsplineElement.CalculatePointsForPostProcessing(tsplineElement);
				for (int i = 0; i < elementPoints.GetLength(0); i++)
				{
					nodes[pointIndex, 0] = elementPoints[i, 0];
					nodes[pointIndex++, 1] = elementPoints[i, 1];
				}
			}

			var elementConnectivity = CreateTsplineConnectivity();

			var pointDisplacements = new double[nodes.GetLength(0), 2];
			var pointStrains = new double[nodes.GetLength(0), 3];
			var pointStresses = new double[nodes.GetLength(0), 3];
			foreach (var element in _model.Elements)
			{
				var localDisplacements = new double[element.ControlPoints.Count, 2];
				var counterCP = 0;
				foreach (var controlPoint in element.ControlPoints)
				{
					localDisplacements[counterCP, 0] =
						(!_model.GlobalDofOrdering.GlobalFreeDofs.Contains(controlPoint, DOFType.X))
							? 0.0
							: _solution[_model.GlobalDofOrdering.GlobalFreeDofs[controlPoint, DOFType.X]];
					localDisplacements[counterCP, 1] =
						(!_model.GlobalDofOrdering.GlobalFreeDofs.Contains(controlPoint, DOFType.Y))
							? 0.0
							: _solution[_model.GlobalDofOrdering.GlobalFreeDofs[controlPoint, DOFType.Y]];
				}

				var elementKnotDisplacements =
					element.ElementType.CalculateDisplacementsForPostProcessing(element, localDisplacements);
				var (knotStrains, knotStresses) =
					element.ElementType.CalculateStressesForPostProcessing(element, localDisplacements);
				for (int i = 0; i < elementConnectivity.GetLength(1); i++)
				{
					var knotConnectivity = elementConnectivity[element.ID, i];
					pointDisplacements[knotConnectivity, 0] = elementKnotDisplacements[i, 0];
					pointDisplacements[knotConnectivity, 1] = elementKnotDisplacements[i, 1];

					pointStrains[knotConnectivity, 0] = knotStrains[i, 0];
					pointStrains[knotConnectivity, 1] = knotStrains[i, 1];
					pointStrains[knotConnectivity, 2] = knotStrains[i, 2];

					pointStresses[knotConnectivity, 0] = knotStresses[i, 0];
					pointStresses[knotConnectivity, 1] = knotStresses[i, 1];
					pointStresses[knotConnectivity, 2] = knotStresses[i, 2];
				}
			}

			Write2DTSplinesFile(nodes, elementConnectivity, pointDisplacements,pointStrains,pointStresses);
		}

		public void Write2DTSplinesFile(double[,] nodeCoordinates, int[,] elementConnectivity, double[,] displacements, double[,] strains, double[,] stresses)
		{
			var numberOfPoints = nodeCoordinates.GetLength(0);
			var numberOfCells = elementConnectivity.GetLength(0);

			int numberOfVerticesPerCell = 0;
			int paraviewCellCode = 0;

			numberOfVerticesPerCell = 4;
			paraviewCellCode = 9;

			using (StreamWriter outputFile = new StreamWriter($"..\\..\\..\\OutputFiles\\{_filename}Paraview.vtu"))
			{
				outputFile.WriteLine("<?xml version=\"1.0\"?>");
				outputFile.WriteLine("<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">");
				outputFile.WriteLine("<UnstructuredGrid>");

				outputFile.WriteLine($"<Piece NumberOfPoints=\"{numberOfPoints}\" NumberOfCells=\"{numberOfCells}\">");
				outputFile.WriteLine($"<PointData Vectors=\"U\">");
				outputFile.WriteLine($"<DataArray type=\"Float32\" Name=\"U\" format=\"ascii\" NumberOfComponents=\"3\">");

				for (int i = 0; i < numberOfPoints; i++)
					outputFile.WriteLine($"{displacements[i, 0]} {displacements[i, 1]} 0.0");

				outputFile.WriteLine("</DataArray>");
				outputFile.WriteLine($"<DataArray type=\"Float32\" Name=\"Strain\" format=\"ascii\" NumberOfComponents=\"3\">");

				for (int i = 0; i < numberOfPoints; i++)
					outputFile.WriteLine($"{strains[i, 0]} {strains[i, 1]} {strains[i, 2]}");


				outputFile.WriteLine("</DataArray>");
				outputFile.WriteLine($"<DataArray type=\"Float32\" Name=\"Stress\" format=\"ascii\" NumberOfComponents=\"3\">");

				for (int i = 0; i < numberOfPoints; i++)
					outputFile.WriteLine($"{stresses[i, 0]} {stresses[i, 1]} {stresses[i, 2]}");


				outputFile.WriteLine("</DataArray>");
				outputFile.WriteLine("</PointData>");
				outputFile.WriteLine("<Points>");
				outputFile.WriteLine("<DataArray type=\"Float32\" NumberOfComponents=\"3\">");

				for (int i = 0; i < numberOfPoints; i++)
					outputFile.WriteLine($"{nodeCoordinates[i, 0]} {nodeCoordinates[i, 1]} 0.0");
				

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

	}
}
