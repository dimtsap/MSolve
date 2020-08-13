using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using ISAAR.MSolve.IGA.Elements;
using ISAAR.MSolve.IGA.Elements.Structural;
using ISAAR.MSolve.IGA.Entities;
using ISAAR.MSolve.IGA.SupportiveClasses;
using ISAAR.MSolve.IGA.SupportiveClasses.Interpolation;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Materials.Interfaces;

namespace ISAAR.MSolve.IGA.Readers
{
    /// <summary>
	/// Read .iga files exported from Rhino.
	/// </summary>
	public class IgaFileReader
	{
		private readonly string _filename;

		private readonly Model _model;

		private int controlPointIDcounter = 0;

		private int elementIDCounter = 0;

		private int numberOfDimensions;

		/// <summary>
		/// Create an .iga file reader.
		/// </summary>
		/// <param name="model">An isogeometric <see cref="Model"/>.</param>
		/// <param name="filename">The name of the file to be read.</param>
		public IgaFileReader(Model model, string filename)
		{
			_model = model;
			_filename = filename;
		}

		private enum Attributes
		{
			type, noden, elemn, node, belem, set
		}

		private enum Types
		{
			plane, surface
		}

        public enum TSplineShellType
        {
            Linear, Section, Thickness
        }

        /// <summary>
        /// Create Model from reading an .iga file.
        /// </summary>
        /// <param name="shellType"><see cref="Enum"/> that specifies the type of the shells that will be generated.</param>
        /// <param name="shellMaterial">The material of the shell.</param>
        /// <param name="thickness">The thickness of the shell.</param>
        public void CreateTSplineShellsModelFromFile(TSplineShellType shellType = TSplineShellType.Linear, 
            IShellMaterial shellMaterial = null, IShellSectionMaterial sectionMaterial=null, double thickness = 1)
		{
			char[] delimeters = { ' ', '=', '\t' };
			Attributes? name = null;

			String[] text = System.IO.File.ReadAllLines(_filename);

			_model.PatchesDictionary.Add(0, new Patch());
			for (int i = 0; i < text.Length; i++)
			{
				var line = text[i].Split(delimeters, StringSplitOptions.RemoveEmptyEntries);
				if (line.Length == 0) continue;
				try
				{
					name = (Attributes)Enum.Parse(typeof(Attributes), line[0].ToLower());
				}
				catch (Exception exception)
				{
					throw new KeyNotFoundException($"Variable name {line[0]} is not found. {exception.Message}");
				}
				switch (name)
				{
					case Attributes.type:
						Types type;
						try
						{
							type = (Types)Enum.Parse(typeof(Types), line[1].ToLower());
						}
						catch (Exception exception)
						{
							throw new KeyNotFoundException($"Variable name {line[0]} is not found. {exception.Message}");
						}
						numberOfDimensions = type == Types.plane ? 2 : 3;
						break;

					case Attributes.noden:
						break;

					case Attributes.elemn:
						var numberOfElements = int.Parse(line[1]);
						break;

					case Attributes.node:
						var controlPoint = new ControlPoint
						{
							ID = controlPointIDcounter,
							X = double.Parse(line[1], CultureInfo.InvariantCulture),
							Y = double.Parse(line[2], CultureInfo.InvariantCulture),
							Z = double.Parse(line[3], CultureInfo.InvariantCulture),
							WeightFactor = double.Parse(line[4], CultureInfo.InvariantCulture)
						};
						_model.ControlPointsDictionary.Add(controlPointIDcounter, controlPoint);
						((List<ControlPoint>)_model.PatchesDictionary[0].ControlPoints).Add(controlPoint);
						controlPointIDcounter++;
						break;

					case Attributes.belem:
						var numberOfElementNodes = int.Parse(line[1]);
						var elementDegreeKsi = int.Parse(line[2]);
						var elementDegreeHeta = int.Parse(line[3]);
						i++;
						line = text[i].Split(delimeters);
						int[] connectivity = new int[numberOfElementNodes];
						for (int j = 0; j < numberOfElementNodes; j++)
							connectivity[j] = Int32.Parse(line[j]);

						var extractionOperator = Matrix.CreateZero(numberOfElementNodes,
							(elementDegreeKsi + 1) * (elementDegreeHeta + 1));
						for (int j = 0; j < numberOfElementNodes; j++)
						{
							line = text[++i].Split(delimeters);
							for (int k = 0; k < (elementDegreeKsi + 1) * (elementDegreeHeta + 1); k++)
							{
								extractionOperator[j, k] = double.Parse(line[k]);
							}
						}

						if (numberOfDimensions == 2)
						{
							throw new NotImplementedException("TSpline2D not yet implemented");
						}
						else
						{
							switch (shellType)
							{
								case TSplineShellType.Linear:
									CreateLinearShell(elementDegreeKsi, elementDegreeHeta, extractionOperator, connectivity, sectionMaterial, thickness);
									break;

								case TSplineShellType.Thickness:
									CreateThicknessShell(elementDegreeKsi, elementDegreeHeta, extractionOperator, connectivity, shellMaterial, thickness);
									break;
							}
						}
						break;

					case Attributes.set:
						break;
				}
			}
		}


		public void CreateTSplineShellsModelFromFile(TSplineShellType shellType = TSplineShellType.Linear, 
            List<IShellMaterial> shellMaterials=null, IShellSectionMaterial sectionMaterial=null, double thickness = 1)
		{
			char[] delimeters = { ' ', '=', '\t' };
			Attributes? name = null;

			String[] text = System.IO.File.ReadAllLines(_filename);

			_model.PatchesDictionary.Add(0, new Patch());
			for (int i = 0; i < text.Length; i++)
			{
				var line = text[i].Split(delimeters, StringSplitOptions.RemoveEmptyEntries);
				if (line.Length == 0) continue;
				try
				{
					name = (Attributes)Enum.Parse(typeof(Attributes), line[0].ToLower());
				}
				catch (Exception exception)
				{
					throw new KeyNotFoundException($"Variable name {line[0]} is not found. {exception.Message}");
				}
				switch (name)
				{
					case Attributes.type:
						Types type;
						try
						{
							type = (Types)Enum.Parse(typeof(Types), line[1].ToLower());
						}
						catch (Exception exception)
						{
							throw new KeyNotFoundException($"Variable name {line[0]} is not found. {exception.Message}");
						}
						numberOfDimensions = type == Types.plane ? 2 : 3;
						break;

					case Attributes.noden:
						break;

					case Attributes.elemn:
						var numberOfElements = int.Parse(line[1]);
						break;

					case Attributes.node:
						var controlPoint = new ControlPoint
						{
							ID = controlPointIDcounter,
							X = double.Parse(line[1], CultureInfo.InvariantCulture),
							Y = double.Parse(line[2], CultureInfo.InvariantCulture),
							Z = double.Parse(line[3], CultureInfo.InvariantCulture),
							WeightFactor = double.Parse(line[4], CultureInfo.InvariantCulture)
						};
						_model.ControlPointsDictionary.Add(controlPointIDcounter, controlPoint);
						((List<ControlPoint>)_model.PatchesDictionary[0].ControlPoints).Add(controlPoint);
						controlPointIDcounter++;
						break;

					case Attributes.belem:
						var numberOfElementNodes = int.Parse(line[1]);
						var elementDegreeKsi = int.Parse(line[2]);
						var elementDegreeHeta = int.Parse(line[3]);
						i++;
						line = text[i].Split(delimeters);
						int[] connectivity = new int[numberOfElementNodes];
						for (int j = 0; j < numberOfElementNodes; j++)
							connectivity[j] = Int32.Parse(line[j]);

						var extractionOperator = Matrix.CreateZero(numberOfElementNodes,
							(elementDegreeKsi + 1) * (elementDegreeHeta + 1));
						for (int j = 0; j < numberOfElementNodes; j++)
						{
							line = text[++i].Split(delimeters);
							for (int k = 0; k < (elementDegreeKsi + 1) * (elementDegreeHeta + 1); k++)
							{
								extractionOperator[j, k] = double.Parse(line[k]);
							}
						}

						if (numberOfDimensions == 2)
						{
							throw new NotImplementedException("TSpline2D not yet implemented");
						}
						else
						{
							switch (shellType)
							{
								case TSplineShellType.Linear:
									CreateLinearShell(elementDegreeKsi, elementDegreeHeta, extractionOperator, connectivity, sectionMaterial, thickness);
									break;

								case TSplineShellType.Thickness:
									CreateThicknessShell(elementDegreeKsi, elementDegreeHeta, extractionOperator, connectivity, shellMaterials, thickness);
									break;
							}
						}
						break;

					case Attributes.set:
						break;
				}
			}
		}


		private void CreateLinearShell(int elementDegreeKsi, int elementDegreeHeta, Matrix extractionOperator,
			int[] connectivity, IShellSectionMaterial material, double thickness)
		{
            var elementControlPoints= connectivity.Select(t => _model.ControlPointsDictionary[t]).ToArray();
            var tsplines = new ShapeTSplines2DFromBezierExtraction(elementDegreeKsi, elementDegreeHeta,
                extractionOperator, elementControlPoints);
			var gauss= new GaussQuadrature();
			var gaussPoints= gauss.CalculateElementGaussPoints(elementDegreeKsi, elementDegreeHeta, new List<Knot>
					{
						new Knot(){ID=0,Ksi=-1,Heta = -1,Zeta = 0},
						new Knot(){ID=1,Ksi=-1,Heta = 1,Zeta = 0},
						new Knot(){ID=2,Ksi=1,Heta = -1,Zeta = 0},
						new Knot(){ID=3,Ksi=1,Heta = 1,Zeta = 0}
					}).ToArray();

			var element = new Element
			{
				ID = elementIDCounter,
				Patch = _model.PatchesDictionary[0],
				ElementType = new KirchhoffLoveShell(material, tsplines, gaussPoints,thickness),
			};
            element.AddControlPoints(elementControlPoints);
			_model.ElementsDictionary.Add(elementIDCounter++, element);
			_model.PatchesDictionary[0].Elements.Add(element);
		}

		private void CreateThicknessShell(int elementDegreeKsi, int elementDegreeHeta, Matrix extractionOperator,
			int[] connectivity, IShellMaterial shellMaterial, double thickness)
		{
            var elementControlPoints = connectivity.Select(t => _model.ControlPointsDictionary[t]).ToArray();
            var tsplines = new ShapeTSplines2DFromBezierExtraction(elementDegreeKsi, elementDegreeHeta,
                extractionOperator, elementControlPoints);
            Element element = new Element
            {
                ID = elementIDCounter,
                Patch = _model.PatchesDictionary[0],
				ElementType = new KirchhoffLoveShellNL(shellMaterial, new List<Knot>
                    {
                        new Knot() {ID = 0, Ksi = -1, Heta = -1, Zeta = 0},
                        new Knot() {ID = 1, Ksi = -1, Heta = 1, Zeta = 0},
                        new Knot() {ID = 2, Ksi = 1, Heta = -1, Zeta = 0},
                        new Knot() {ID = 3, Ksi = 1, Heta = 1, Zeta = 0},
                    }, tsplines, elementControlPoints,
					_model.PatchesDictionary[0], thickness, elementDegreeKsi, elementDegreeHeta)
			};
                
            element.AddControlPoints(elementControlPoints);
			_model.ElementsDictionary.Add(elementIDCounter++, element);
			_model.PatchesDictionary[0].Elements.Add(element);
		}

        private void CreateThicknessShell(int elementDegreeKsi, int elementDegreeHeta, Matrix extractionOperator,
            int[] connectivity, List<IShellMaterial> shellMaterials, double thickness)
        {
            var elementControlPoints = connectivity.Select(t => _model.ControlPointsDictionary[t]).ToArray();
            var tsplines = new ShapeTSplines2DFromBezierExtraction(elementDegreeKsi, elementDegreeHeta,
                extractionOperator, elementControlPoints);
            Element element = new Element
            {
                ID = elementIDCounter,
                Patch = _model.PatchesDictionary[0],
                ElementType = new KirchhoffLoveShellNL(shellMaterials, new List<Knot>
                    {
                        new Knot() {ID = 0, Ksi = -1, Heta = -1, Zeta = 0},
                        new Knot() {ID = 1, Ksi = -1, Heta = 1, Zeta = 0},
                        new Knot() {ID = 2, Ksi = 1, Heta = -1, Zeta = 0},
                        new Knot() {ID = 3, Ksi = 1, Heta = 1, Zeta = 0},
                    }, tsplines, elementControlPoints,
                    _model.PatchesDictionary[0], thickness, elementDegreeKsi, elementDegreeHeta)
            };
                
            element.AddControlPoints(elementControlPoints);
            _model.ElementsDictionary.Add(elementIDCounter++, element);
            _model.PatchesDictionary[0].Elements.Add(element);
        }
	}
}
