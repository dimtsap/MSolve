﻿using System;
using System.Collections.Generic;
using System.Linq;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.IGA.Entities;
using ISAAR.MSolve.IGA.Entities.Loads;
using ISAAR.MSolve.IGA.Interfaces;
using ISAAR.MSolve.IGA.Problems.SupportiveClasses;
using ISAAR.MSolve.IGA.SupportiveClasses;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Materials.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;

namespace ISAAR.MSolve.IGA.Elements
{
	public class TSplineKirchhoffLoveShellElementMaterial : Element, IStructuralIsogeometricElement
	{
		public Matrix2D ExtractionOperator { get; set; }
		public int DegreeKsi { get; set; }
		public int DegreeHeta { get; set; }
		protected readonly static DOFType[] controlPointDOFTypes = new DOFType[] { DOFType.X, DOFType.Y, DOFType.Z };
		protected DOFType[][] dofTypes;
		protected IElementDOFEnumerator dofEnumerator = new GenericDOFEnumerator();
		private DynamicMaterial dynamicProperties;
		private IReadOnlyList<IShellMaterial> materialsAtGaussPoints;

		private Dictionary<GaussLegendrePoint3D, Dictionary<GaussLegendrePoint3D, IShellMaterial>>
			materialsAtThicknessGP =
				new Dictionary<GaussLegendrePoint3D, Dictionary<GaussLegendrePoint3D, IShellMaterial>>();


		public TSplineKirchhoffLoveShellElementMaterial(int id, Patch patch, int degreeKsi, int degreeHeta,
			double thickness, Matrix2D extractionOperator, IShellMaterial shellMaterial)
		{
			this.ID = id;
			this.Patch = patch;
			this.DegreeKsi = degreeKsi;
			this.DegreeHeta = degreeHeta;
			this.Thickness = thickness;
			this.ExtractionOperator = extractionOperator;

			CreateElementGaussPoints(this);
			foreach (var medianSurfaceGP in thicknessIntegrationPoints.Keys)
			{
				materialsAtThicknessGP.Add(medianSurfaceGP, new Dictionary<GaussLegendrePoint3D, IShellMaterial>());
				foreach (var point in thicknessIntegrationPoints[medianSurfaceGP])
				{
					materialsAtThicknessGP[medianSurfaceGP].Add(point, shellMaterial.Clone());
				}
			}
		}

		public IElementDOFEnumerator DOFEnumerator
		{
			get
			{
				return dofEnumerator;
			}

			set
			{
				this.dofEnumerator = value;
			}
		}

		public ElementDimensions ElementDimensions
		{
			get
			{
				return ElementDimensions.ThreeD;
			}
		}

		public double Thickness { get; set; }

		public bool MaterialModified
		{
			get
			{
				for (int j = 0; j < materialsAtThicknessGP.Keys.Count; j++)
				{
					foreach (var keyValuePair in materialsAtThicknessGP[materialsAtThicknessGP.Keys.ToList()[j]])
					{
						var material = keyValuePair.Value;
						if (material.Modified)
						{
							return true;
						}
					}
				}
				return false;
			}
		}

		public double[] CalculateForces(Element element, double[] localDisplacements, double[] localdDisplacements)
		{
			var shellElement = (TSplineKirchhoffLoveShellElementMaterial)element;
			var gaussPoints = materialsAtThicknessGP.Keys.ToList();
			Vector ElementNodalForces = new Vector(shellElement.ControlPointsDictionary.Count * 3);
			ShapeTSplines2DFromBezierExtraction tsplines = new ShapeTSplines2DFromBezierExtraction(shellElement, shellElement.ControlPoints);


			for (int j = 0; j < gaussPoints.Count; j++)
			{
				var jacobianMatrix = CalculateJacobian(shellElement, tsplines, j);

				var hessianMatrix = CalculateHessian(shellElement, tsplines, j);

				var surfaceBasisVector1 = CalculateSurfaceBasisVector1(jacobianMatrix, 0);

				var surfaceBasisVector2 = CalculateSurfaceBasisVector1(jacobianMatrix, 1);

				var surfaceBasisVector3 = new[]
				{
					surfaceBasisVector1[1] * surfaceBasisVector2[2] - surfaceBasisVector1[2] * surfaceBasisVector2[1],
					surfaceBasisVector1[2] * surfaceBasisVector2[0] - surfaceBasisVector1[0] * surfaceBasisVector2[2],
					surfaceBasisVector1[0] * surfaceBasisVector2[1] - surfaceBasisVector1[1] * surfaceBasisVector2[0]
				};

				double norm = 0;
				for (int i = 0; i < surfaceBasisVector3.Length; i++)
					norm += surfaceBasisVector3[i] * surfaceBasisVector3[i];
				var J1 = Math.Sqrt(norm);

				for (int i = 0; i < surfaceBasisVector3.Length; i++)
					surfaceBasisVector3[i] = surfaceBasisVector3[i] / J1;

				var surfaceBasisVectorDerivative1 = CalculateSurfaceBasisVector1(hessianMatrix, 0);
				var surfaceBasisVectorDerivative2 = CalculateSurfaceBasisVector1(hessianMatrix, 1);
				var surfaceBasisVectorDerivative12 = CalculateSurfaceBasisVector1(hessianMatrix, 2);

				var Bmembrane = CalculateMembraneDeformationMatrix(tsplines, j, surfaceBasisVector1, surfaceBasisVector2, shellElement);
				var Bbending = CalculateBendingDeformationMatrix(surfaceBasisVector3, tsplines, j, surfaceBasisVector2, surfaceBasisVectorDerivative1, surfaceBasisVector1, J1, surfaceBasisVectorDerivative2, surfaceBasisVectorDerivative12, shellElement);

				var (MembraneForces, BendingMoments) =
					IntegratedStressesOverThickness(gaussPoints, j);

				var wfactor = J1 * gaussPoints[j].WeightFactor;
				for (int i = 0; i < Bmembrane.GetLength(1); i++)
				{
					for (int k = 0; k < Bmembrane.GetLength(0); k++)
					{
						ElementNodalForces[i] += Bmembrane[k, i] * MembraneForces[k] * wfactor +
												 Bbending[k, i] * BendingMoments[k] * wfactor;
					}
				}
			}

			return ElementNodalForces.Data;
		}

		public double[] CalculateForcesForLogging(Element element, double[] localDisplacements)
		{
			throw new NotImplementedException();
		}


		public Dictionary<int, double> CalculateLoadingCondition(Element element, Edge edge, NeumannBoundaryCondition neumann)
		{
			throw new NotImplementedException();
		}

		public Dictionary<int, double> CalculateLoadingCondition(Element element, Face face, NeumannBoundaryCondition neumann)
		{
			throw new NotImplementedException();
		}

		public Dictionary<int, double> CalculateLoadingCondition(Element element, Edge edge, PressureBoundaryCondition pressure)
		{
			throw new NotImplementedException();
		}

		public Dictionary<int, double> CalculateLoadingCondition(Element element, Face face, PressureBoundaryCondition pressure)
		{
			throw new NotImplementedException();
		}

		public Tuple<double[], double[]> CalculateStresses(Element element, double[] localDisplacements, double[] localdDisplacements)
		{
			var shellElement = (TSplineKirchhoffLoveShellElementMaterial)element;
			//Matrix2D stiffnessMatrixElement = new Matrix2D(shellElement.ControlPointsDictionary.Count * 3, shellElement.ControlPointsDictionary.Count * 3);

			ShapeTSplines2DFromBezierExtraction tsplines = new ShapeTSplines2DFromBezierExtraction(shellElement, shellElement.ControlPoints);

			for (int j = 0; j < materialsAtThicknessGP.Keys.Count; j++)
			{
				var jacobianMatrix = CalculateJacobian(shellElement, tsplines, j);

				var hessianMatrix = CalculateHessian(shellElement, tsplines, j);

				var surfaceBasisVector1 = CalculateSurfaceBasisVector1(jacobianMatrix, 0);

				var surfaceBasisVector2 = CalculateSurfaceBasisVector1(jacobianMatrix, 1);

				var surfaceBasisVector3 = new[]
				{
					surfaceBasisVector1[1] * surfaceBasisVector2[2] - surfaceBasisVector1[2] * surfaceBasisVector2[1],
					surfaceBasisVector1[2] * surfaceBasisVector2[0] - surfaceBasisVector1[0] * surfaceBasisVector2[2],
					surfaceBasisVector1[0] * surfaceBasisVector2[1] - surfaceBasisVector1[1] * surfaceBasisVector2[0]
				};

				double norm = 0;
				for (int i = 0; i < surfaceBasisVector3.Length; i++)
					norm += surfaceBasisVector3[i] * surfaceBasisVector3[i];
				var J1 = Math.Sqrt(norm);

				var surfaceBasisVectorDerivative1 = CalculateSurfaceBasisVector1(hessianMatrix, 0);
				var surfaceBasisVectorDerivative2 = CalculateSurfaceBasisVector1(hessianMatrix, 1);
				var surfaceBasisVectorDerivative12 = CalculateSurfaceBasisVector1(hessianMatrix, 2);

				var Bmembrane = CalculateMembraneDeformationMatrix(tsplines, j, surfaceBasisVector1, surfaceBasisVector2, shellElement);
				var Bbending = CalculateBendingDeformationMatrix(surfaceBasisVector3, tsplines, j, surfaceBasisVector2,
					surfaceBasisVectorDerivative1, surfaceBasisVector1, J1, surfaceBasisVectorDerivative2,
					surfaceBasisVectorDerivative12, shellElement);

				var membraneStrain = new double[Bmembrane.GetLength(0)];
				var bendingStrain = new double[Bmembrane.GetLength(0)];
				for (int i = 0; i < Bmembrane.GetLength(0); i++)
				{
					for (int k = 0; k < Bmembrane.GetLength(1); k++)
					{
						membraneStrain[i] += Bmembrane[i,k] * localDisplacements[k];
						bendingStrain[i] += Bbending[i,k] * localDisplacements[k];
					}
				}

				foreach (var keyValuePair in materialsAtThicknessGP[materialsAtThicknessGP.Keys.ToList()[j]])
				{
					var thicknessPoint = keyValuePair.Key;
					var material = keyValuePair.Value;
					var gpStrain = new double[bendingStrain.Length];
					var z = thicknessPoint.Zeta;
					for (int i = 0; i < bendingStrain.Length; i++)
					{
						gpStrain[i] += membraneStrain[i] + bendingStrain[i] * z;
					}
					material.UpdateMaterial(gpStrain);
				}


			}
			return new Tuple<double[], double[]>(new double[0], new double[0]);

		}

		public void ClearMaterialState()
		{
			throw new NotImplementedException();
		}

		public void SaveMaterialState()
		{
			for (int j = 0; j < materialsAtThicknessGP.Keys.Count; j++)
			{
				foreach (var keyValuePair in materialsAtThicknessGP[materialsAtThicknessGP.Keys.ToList()[j]])
				{
					var thicknessPoint = keyValuePair.Key;
					var material = keyValuePair.Value;
					material.SaveState();
				}
			}
		}

        public Dictionary<GaussLegendrePoint3D, double[]> GetStressesAtLowerSurfaceForLogging()
        {
            Dictionary<GaussLegendrePoint3D, double[]> StressesAtLowerSurfaceGpS = new Dictionary<GaussLegendrePoint3D, double[]>();

            for (int j = 0; j < materialsAtThicknessGP.Keys.Count; j++)
            {
                double[] lowerSurfaceStresses = materialsAtThicknessGP[materialsAtThicknessGP.Keys.ToList()[j]].ElementAt(0).Value.Stresses;
                StressesAtLowerSurfaceGpS.Add(materialsAtThicknessGP.Keys.ToList()[j], lowerSurfaceStresses);
            }

            return StressesAtLowerSurfaceGpS;
        }

        public Dictionary<GaussLegendrePoint3D, double[]> GetStressesAtMidSurfaceForLogging()
        {
            Dictionary<GaussLegendrePoint3D, double[]> StressesAtLowerSurfaceGpS = new Dictionary<GaussLegendrePoint3D, double[]>();

            for (int j = 0; j < materialsAtThicknessGP.Keys.Count; j++)
            {
                double[] lowerSurfaceStresses = materialsAtThicknessGP[materialsAtThicknessGP.Keys.ToList()[j]].ElementAt(1).Value.Stresses;
                StressesAtLowerSurfaceGpS.Add(materialsAtThicknessGP.Keys.ToList()[j], lowerSurfaceStresses);
            }

            return StressesAtLowerSurfaceGpS;
        }


        public void ClearMaterialStresses()
		{
			//TODO: the next throws an exception. Investigate. Possible changes in Analyzers may be the cause.
			//for (int j = 0; j < materialsAtThicknessGP.Keys.Count; j++)
			//{
			//    foreach (var keyValuePair in materialsAtThicknessGP[materialsAtThicknessGP.Keys.ToList()[j]])
			//    {
			//        var thicknessPoint = keyValuePair.Key;
			//        var material = keyValuePair.Value;
			//        material.ClearState();
			//    }
			//}
		}

		public IMatrix2D DampingMatrix(IElement element)
		{
			throw new NotImplementedException();
		}

		public IList<IList<DOFType>> GetElementDOFTypes(IElement element)
		{
			var nurbsElement = (TSplineKirchhoffLoveShellElementMaterial)element;
			dofTypes = new DOFType[nurbsElement.ControlPoints.Count][];
			for (int i = 0; i < nurbsElement.ControlPoints.Count; i++)
			{
				dofTypes[i] = controlPointDOFTypes;
			}
			return dofTypes;
		}

		public IMatrix2D MassMatrix(IElement element)
		{
			throw new NotImplementedException();
		}

		public void ResetMaterialModified()
		{
			throw new NotImplementedException();
		}
		

		public IMatrix2D StiffnessMatrix(IElement element)
		{
			var shellElement = (TSplineKirchhoffLoveShellElementMaterial)element;
			var gaussPoints = materialsAtThicknessGP.Keys.ToArray();
			ShapeTSplines2DFromBezierExtraction tsplines =
				new ShapeTSplines2DFromBezierExtraction(shellElement, shellElement.ControlPoints);


			var bRows = 3;
			var bCols = shellElement.ControlPoints.Count * 3;
			var stiffnessMatrix = new double[bCols, bCols];
			var BmTranspose = new double[bCols, bRows];
			var BbTranspose = new double[bCols, bRows];

			var BmTransposeMultStiffness = new double[bCols, bRows];
			var BbTransposeMultStiffness = new double[bCols, bRows];
			var BmbTransposeMultStiffness = new double[bCols, bRows];
			var BbmTransposeMultStiffness = new double[bCols, bRows];

			for (int j = 0; j < gaussPoints.Length; j++)
			{
				var jacobianMatrix = CalculateJacobian(shellElement, tsplines, j);

				var hessianMatrix = CalculateHessian(shellElement, tsplines, j);
				var surfaceBasisVector1 = CalculateSurfaceBasisVector1(jacobianMatrix, 0);

				var surfaceBasisVector2 = CalculateSurfaceBasisVector1(jacobianMatrix, 1);

				var surfaceBasisVector3 = new[]
				{
					surfaceBasisVector1[1] * surfaceBasisVector2[2] - surfaceBasisVector1[2] * surfaceBasisVector2[1],
					surfaceBasisVector1[2] * surfaceBasisVector2[0] - surfaceBasisVector1[0] * surfaceBasisVector2[2],
					surfaceBasisVector1[0] * surfaceBasisVector2[1] - surfaceBasisVector1[1] * surfaceBasisVector2[0]
				};

				foreach (var integrationPointMaterial in materialsAtThicknessGP[gaussPoints[j]].Values)
				{
					integrationPointMaterial.TangentVectorV1 = surfaceBasisVector1;
					integrationPointMaterial.TangentVectorV2 = surfaceBasisVector2;
					integrationPointMaterial.NormalVectorV3 = surfaceBasisVector3;
				}

				double norm = 0;
				for (int i = 0; i < surfaceBasisVector3.Length; i++)
					norm += surfaceBasisVector3[i] * surfaceBasisVector3[i];
				var J1 = Math.Sqrt(norm);

				for (int i = 0; i < surfaceBasisVector3.Length; i++)
					surfaceBasisVector3[i] = surfaceBasisVector3[i] / J1;

				var surfaceBasisVectorDerivative1 = CalculateSurfaceBasisVector1(hessianMatrix, 0);
				var surfaceBasisVectorDerivative2 = CalculateSurfaceBasisVector1(hessianMatrix, 1);
				var surfaceBasisVectorDerivative12 = CalculateSurfaceBasisVector1(hessianMatrix, 2);

				var Bmembrane = CalculateMembraneDeformationMatrix(tsplines, j, surfaceBasisVector1, surfaceBasisVector2, shellElement);
				var Bbending = CalculateBendingDeformationMatrix(surfaceBasisVector3, tsplines, j, surfaceBasisVector2,
					surfaceBasisVectorDerivative1, surfaceBasisVector1, J1, surfaceBasisVectorDerivative2,
					surfaceBasisVectorDerivative12, shellElement);

				var (MembraneConstitutiveMatrix, BendingConstitutiveMatrix, CouplingConstitutiveMatrix) =
					IntegratedConstitutiveOverThickness(gaussPoints, j);

				double wFactor = J1 * gaussPoints[j].WeightFactor;
				double tempb = 0;
				double tempm = 0;
				Array.Clear(BmTranspose, 0, bRows * bCols);
				Array.Clear(BbTranspose, 0, bRows * bCols);
				for (int i = 0; i < bRows; i++)
				{
					for (int k = 0; k < bCols; k++)
					{
						BmTranspose[k, i] = Bmembrane[i, k] * wFactor;
						BbTranspose[k, i] = Bbending[i, k] * wFactor;
					}
				}

				double tempcm = 0;
				double tempcb = 0;
				double tempcc = 0;
				Array.Clear(BmTransposeMultStiffness, 0, bRows * bCols);
				Array.Clear(BbTransposeMultStiffness, 0, bRows * bCols);
				Array.Clear(BmbTransposeMultStiffness, 0, bRows * bCols);
				Array.Clear(BbmTransposeMultStiffness, 0, bRows * bCols);
				for (int i = 0; i < bCols; i++)
				{
					for (int k = 0; k < bRows; k++)
					{
						tempm = BmTranspose[i, k];
						tempb = BbTranspose[i, k];
						for (int m = 0; m < bRows; m++)
						{
							tempcm = MembraneConstitutiveMatrix[k, m];
							tempcb = BendingConstitutiveMatrix[k, m];
							tempcc = CouplingConstitutiveMatrix[k, m];

							BmTransposeMultStiffness[i, m] += tempm * tempcm;
							BbTransposeMultStiffness[i, m] += tempb * tempcb;
							BmbTransposeMultStiffness[i, m] += tempm * tempcc;
							BbmTransposeMultStiffness[i, m] += tempb * tempcc;
						}
					}
				}

				double tempmb = 0;
				double tempbm = 0;
				double mem = 0;
				double ben = 0;
				for (int i = 0; i < bCols; i++)
				{
					for (int k = 0; k < bRows; k++)
					{
						tempm = BmTransposeMultStiffness[i, k];
						tempb = BbTransposeMultStiffness[i, k];
						tempmb = BmbTransposeMultStiffness[i, k];
						tempbm = BbmTransposeMultStiffness[i, k];

						for (int m = 0; m < bCols; m++)
						{
							mem = Bmembrane[k, m];
							ben = Bbending[k, m];
							stiffnessMatrix[i, m] += tempm * mem + tempb * ben + tempmb * ben + tempbm * mem;
						}
					}
				}
			}
			return new Matrix2D(stiffnessMatrix);
		}

		public (double[,] MembraneConstitutiveMatrix, double[,] BendingConstitutiveMatrix, double[,] CouplingConstitutiveMatrix) IntegratedConstitutiveOverThickness(IList<GaussLegendrePoint3D> gaussPoints, int j)
		{
			var MembraneConstitutiveMatrix = new double[3, 3];
			var BendingConstitutiveMatrix = new double[3, 3];
			var CouplingConstitutiveMatrix = new double[3, 3];

			foreach (var keyValuePair in materialsAtThicknessGP[gaussPoints[j]])
			{
				var thicknessPoint = keyValuePair.Key;
				var material = keyValuePair.Value;
				var constitutiveMatrixM = ((Matrix2D)material.ConstitutiveMatrix).Data;
				double tempc = 0;
				double w = thicknessPoint.WeightFactor;
				double z = thicknessPoint.Zeta;
				for (int i = 0; i < 3; i++)
				{
					for (int k = 0; k < 3; k++)
					{
						tempc = constitutiveMatrixM[i, k];
						MembraneConstitutiveMatrix[i, k] += tempc * w;
						CouplingConstitutiveMatrix[i, k] += tempc * w * z;
						BendingConstitutiveMatrix[i, k] += tempc * w * z * z;
					}
				}
			}

			return (MembraneConstitutiveMatrix, BendingConstitutiveMatrix, CouplingConstitutiveMatrix);
		}

		private double[,] CopyConstitutiveMatrix(double[,] f)
		{
			var g = new double[f.GetLength(0), f.GetLength(1)];
			Array.Copy(f, 0, g, 0, f.Length);
			return g;
		}

		public (Vector MembraneForces, Vector BendingMoments) IntegratedStressesOverThickness(
			IList<GaussLegendrePoint3D> gaussPoints, int j)
		{
			var MembraneForces = new Vector(3);
			var BendingMoments = new Vector(3);

			foreach (var keyValuePair in materialsAtThicknessGP[gaussPoints[j]])
			{
				var thicknessPoint = keyValuePair.Key;
				var material = keyValuePair.Value;
				var gpForces = new Vector(material.Stresses.GetLength(0));
				material.Stresses.CopyTo(gpForces.Data, 0);
				gpForces.Scale(thicknessPoint.WeightFactor);
				MembraneForces.Add(gpForces);

				var gpMoments = new Vector(material.Stresses.GetLength(0));
				material.Stresses.CopyTo(gpMoments.Data, 0);
				gpMoments.Scale(thicknessPoint.WeightFactor *
				                Math.Pow(thicknessPoint.Zeta, 1));
				BendingMoments.Add(gpMoments);

			}

			return (MembraneForces, BendingMoments);
		}


		private double[,] CalculateBendingDeformationMatrix(double[] surfaceBasisVector3, ShapeTSplines2DFromBezierExtraction tsplines, int j,
			double[] surfaceBasisVector2, double[] surfaceBasisVectorDerivative1, double[] surfaceBasisVector1, double J1,
			double[] surfaceBasisVectorDerivative2, double[] surfaceBasisVectorDerivative12, TSplineKirchhoffLoveShellElementMaterial element)
		{
			var Bbending = new double[3, element.ControlPoints.Count * 3];
			var s1 = new Vector(surfaceBasisVector1);
			var s2 = new Vector(surfaceBasisVector2);
			var s3 = new Vector(surfaceBasisVector3);
			var s11 = new Vector(surfaceBasisVectorDerivative1);
			var s22 = new Vector(surfaceBasisVectorDerivative2);
			var s12 = new Vector(surfaceBasisVectorDerivative12);
			for (int column = 0; column < element.ControlPoints.Count * 3; column += 3)
			{
				#region BI1

				var BI1 = s3 ^ s3;
				BI1.Multiply(tsplines.TSplineDerivativeValuesHeta[column / 3, j]);
				var auxVector = s2 ^ s3;
				auxVector.Multiply(tsplines.TSplineDerivativeValuesKsi[column / 3, j]);
				BI1.Add(auxVector);
				BI1.Multiply(s3.DotProduct(s11));
				auxVector = s1 ^ s11;
				auxVector.Multiply(tsplines.TSplineDerivativeValuesHeta[column / 3, j]);
				BI1.Add(auxVector);
				BI1.Multiply(1 / J1);
				auxVector[0] = surfaceBasisVector3[0];
				auxVector[1] = surfaceBasisVector3[1];
				auxVector[2] = surfaceBasisVector3[2];
				auxVector.Multiply(-tsplines.TSplineSecondDerivativesValueKsi[column / 3, j]);
				BI1.Add(auxVector);

				#endregion

				#region BI2

				Vector BI2 = s3 ^ s3;
				BI2.Multiply(tsplines.TSplineDerivativeValuesHeta[column / 3, j]);
				auxVector = s2 ^ s3;
				auxVector.Multiply(tsplines.TSplineDerivativeValuesKsi[column / 3, j]);
				BI2.Add(auxVector);
				BI2.Multiply(s3.DotProduct(s22));
				auxVector = s1 ^ s22;
				auxVector.Multiply(tsplines.TSplineDerivativeValuesHeta[column / 3, j]);
				BI2.Add(auxVector);
				auxVector = s22 ^ s2;
				auxVector.Multiply(tsplines.TSplineDerivativeValuesKsi[column / 3, j]);
				BI2.Add(auxVector);
				BI2.Multiply(1 / J1);
				auxVector[0] = surfaceBasisVector3[0];
				auxVector[1] = surfaceBasisVector3[1];
				auxVector[2] = surfaceBasisVector3[2];
				auxVector.Multiply(-tsplines.TSplineSecondDerivativesValueHeta[column / 3, j]);
				BI2.Add(auxVector);

				#endregion

				#region BI3

				Vector BI3 = s3 ^ s3;
				BI3.Multiply(tsplines.TSplineDerivativeValuesHeta[column / 3, j]);
				auxVector = s2 ^ s3;
				auxVector.Multiply(tsplines.TSplineDerivativeValuesKsi[column / 3, j]);
				BI3.Add(auxVector);
				BI3.Multiply(s3.DotProduct(s12));
				auxVector = s1 ^ s12;
				auxVector.Multiply(tsplines.TSplineDerivativeValuesHeta[column / 3, j]);
				BI3.Add(auxVector);
				auxVector = s22 ^ s2;
				auxVector.Multiply(tsplines.TSplineDerivativeValuesKsi[column / 3, j]);
				BI3.Add(auxVector);
				BI3.Multiply(1 / J1);
				auxVector[0] = surfaceBasisVector3[0];
				auxVector[1] = surfaceBasisVector3[1];
				auxVector[2] = surfaceBasisVector3[2];
				auxVector.Multiply(-tsplines.TSplineSecondDerivativesValueKsiHeta[column / 3, j]);
				BI3.Add(auxVector);

				#endregion

				Bbending[0, column] = BI1[0];
				Bbending[0, column + 1] = BI1[1];
				Bbending[0, column + 2] = BI1[2];

				Bbending[1, column] = BI2[0];
				Bbending[1, column + 1] = BI2[1];
				Bbending[1, column + 2] = BI2[2];

				Bbending[2, column] = 2 * BI3[0];
				Bbending[2, column + 1] = 2 * BI3[1];
				Bbending[2, column + 2] = 2 * BI3[2];
			}

			return Bbending;
		}

		private double[,] CalculateMembraneDeformationMatrix(ShapeTSplines2DFromBezierExtraction tsplines, int j, double[] surfaceBasisVector1,
			double[] surfaceBasisVector2, TSplineKirchhoffLoveShellElementMaterial element)
		{
			var dRIa = new double[3, element.ControlPoints.Count * 3];
			for (int i = 0; i < element.ControlPoints.Count; i++)
			{
				for (int m = 0; m < 3; m++)
				{
					dRIa[m, i] = tsplines.TSplineDerivativeValuesHeta[i, j] * surfaceBasisVector1[m] +
								 tsplines.TSplineDerivativeValuesKsi[i, j] * surfaceBasisVector2[m];
				}
			}

			var Bmembrane = new double[3, element.ControlPoints.Count * 3];
			for (int column = 0; column < element.ControlPoints.Count * 3; column += 3)
			{
				Bmembrane[0, column] = tsplines.TSplineDerivativeValuesKsi[column / 3, j] * surfaceBasisVector1[0];
				Bmembrane[0, column + 1] = tsplines.TSplineDerivativeValuesKsi[column / 3, j] * surfaceBasisVector1[1];
				Bmembrane[0, column + 2] = tsplines.TSplineDerivativeValuesKsi[column / 3, j] * surfaceBasisVector1[2];

				Bmembrane[1, column] = tsplines.TSplineDerivativeValuesHeta[column / 3, j] * surfaceBasisVector2[0];
				Bmembrane[1, column + 1] = tsplines.TSplineDerivativeValuesHeta[column / 3, j] * surfaceBasisVector2[1];
				Bmembrane[1, column + 2] = tsplines.TSplineDerivativeValuesHeta[column / 3, j] * surfaceBasisVector2[2];

				Bmembrane[2, column] = dRIa[0, column / 3];
				Bmembrane[2, column + 1] = dRIa[1, column / 3];
				Bmembrane[2, column + 2] = dRIa[2, column / 3];
			}

			return Bmembrane;
		}

		private static double[] CalculateSurfaceBasisVector1(double[,] Matrix, int row)
		{
			var surfaceBasisVector1 = new double[3];
			surfaceBasisVector1[0] = Matrix[row, 0];
			surfaceBasisVector1[1] = Matrix[row, 1];
			surfaceBasisVector1[2] = Matrix[row, 2];
			return surfaceBasisVector1;
		}

		private static double[,] CalculateHessian(TSplineKirchhoffLoveShellElementMaterial shellElement, ShapeTSplines2DFromBezierExtraction tsplines, int j)
		{
			double[,] hessianMatrix = new double[3, 3];
			for (int k = 0; k < shellElement.ControlPoints.Count; k++)
			{
				hessianMatrix[0, 0] += tsplines.TSplineSecondDerivativesValueKsi[k, j] * shellElement.ControlPoints[k].X;
				hessianMatrix[0, 1] += tsplines.TSplineSecondDerivativesValueKsi[k, j] * shellElement.ControlPoints[k].Y;
				hessianMatrix[0, 2] += tsplines.TSplineSecondDerivativesValueKsi[k, j] * shellElement.ControlPoints[k].Z;
				hessianMatrix[1, 0] += tsplines.TSplineSecondDerivativesValueHeta[k, j] * shellElement.ControlPoints[k].X;
				hessianMatrix[1, 1] += tsplines.TSplineSecondDerivativesValueHeta[k, j] * shellElement.ControlPoints[k].Y;
				hessianMatrix[1, 2] += tsplines.TSplineSecondDerivativesValueHeta[k, j] * shellElement.ControlPoints[k].Z;
				hessianMatrix[2, 0] += tsplines.TSplineSecondDerivativesValueKsiHeta[k, j] * shellElement.ControlPoints[k].X;
				hessianMatrix[2, 1] += tsplines.TSplineSecondDerivativesValueKsiHeta[k, j] * shellElement.ControlPoints[k].Y;
				hessianMatrix[2, 2] += tsplines.TSplineSecondDerivativesValueKsiHeta[k, j] * shellElement.ControlPoints[k].Z;
			}

			return hessianMatrix;
		}

		private static double[,] CalculateJacobian(TSplineKirchhoffLoveShellElementMaterial shellElement, ShapeTSplines2DFromBezierExtraction tsplines, int j)
		{
			double[,] jacobianMatrix = new double[2, 3];
			for (int k = 0; k < shellElement.ControlPoints.Count; k++)
			{
				jacobianMatrix[0, 0] += tsplines.TSplineDerivativeValuesKsi[k, j] * shellElement.ControlPoints[k].X;
				jacobianMatrix[0, 1] += tsplines.TSplineDerivativeValuesKsi[k, j] * shellElement.ControlPoints[k].Y;
				jacobianMatrix[0, 2] += tsplines.TSplineDerivativeValuesKsi[k, j] * shellElement.ControlPoints[k].Z;
				jacobianMatrix[1, 0] += tsplines.TSplineDerivativeValuesHeta[k, j] * shellElement.ControlPoints[k].X;
				jacobianMatrix[1, 1] += tsplines.TSplineDerivativeValuesHeta[k, j] * shellElement.ControlPoints[k].Y;
				jacobianMatrix[1, 2] += tsplines.TSplineDerivativeValuesHeta[k, j] * shellElement.ControlPoints[k].Z;
			}

			return jacobianMatrix;
		}


		const int thicknessIntegrationDegree = 2;
		Dictionary<GaussLegendrePoint3D, List<GaussLegendrePoint3D>> thicknessIntegrationPoints = new Dictionary<GaussLegendrePoint3D, List<GaussLegendrePoint3D>>();
		private IList<GaussLegendrePoint3D> CreateElementGaussPoints(TSplineKirchhoffLoveShellElementMaterial element)
		{
			GaussQuadrature gauss = new GaussQuadrature();
			var medianSurfaceGP = gauss.CalculateElementGaussPoints(element.DegreeKsi, element.DegreeHeta, new List<Knot>
			{
				new Knot() {ID = 0, Ksi = -1, Heta = -1, Zeta = 0},
				new Knot() {ID = 1, Ksi = -1, Heta = 1, Zeta = 0},
				new Knot() {ID = 2, Ksi = 1, Heta = -1, Zeta = 0},
				new Knot() {ID = 3, Ksi = 1, Heta = 1, Zeta = 0},
			});

			foreach (var point in medianSurfaceGP)
			{
				thicknessIntegrationPoints.Add(point, gauss.CalculateElementGaussPoints(thicknessIntegrationDegree,
					new List<Knot>
					{
						new Knot() {ID = 0, Ksi = -Thickness / 2, Heta = point.Heta },
						new Knot() {ID = 1, Ksi = Thickness / 2, Heta = point.Heta},
					}).Select(gp => new GaussLegendrePoint3D(point.Ksi, point.Heta, gp.Ksi, null, gp.WeightFactor)).ToList());
			}

			//new Knot() { ID = 0, Ksi = -1, Heta = -1, Zeta = -Thickness / 2 },
			//new Knot() { ID = 1, Ksi = -1, Heta = -1, Zeta = Thickness / 2 },
			//new Knot() { ID = 2, Ksi = -1, Heta = 1, Zeta = -Thickness / 2 },
			//new Knot() { ID = 3, Ksi = -1, Heta = 1, Zeta = Thickness / 2 },
			//new Knot() { ID = 4, Ksi = 1, Heta = -1, Zeta = -Thickness / 2 },
			//new Knot() { ID = 5, Ksi = 1, Heta = -1, Zeta = Thickness / 2 },
			//new Knot() { ID = 6, Ksi = 1, Heta = 1, Zeta = -Thickness / 2 },
			//new Knot() { ID = 7, Ksi = 1, Heta = 1, Zeta = Thickness / 2 }


			return medianSurfaceGP;
		}


		public double[,] CalculateDisplacementsForPostProcessing(Element element, double[,] localDisplacements)
		{
			var tsplineElement = (TSplineKirchhoffLoveShellElementMaterial)element;
			var knotParametricCoordinatesKsi = new Vector(new double[] { -1, 1 });
			var knotParametricCoordinatesHeta = new Vector(new double[] { -1, 1 });

			ShapeTSplines2DFromBezierExtraction tsplines = new ShapeTSplines2DFromBezierExtraction(tsplineElement, tsplineElement.ControlPoints, knotParametricCoordinatesKsi, knotParametricCoordinatesHeta);

			var knotDisplacements = new double[4, 3];
			var paraviewKnotRenumbering = new int[] { 0, 3, 1, 2 };
			for (int j = 0; j < knotDisplacements.GetLength(0); j++)
			{
				for (int i = 0; i < element.ControlPoints.Count; i++)
				{
					knotDisplacements[paraviewKnotRenumbering[j], 0] += tsplines.TSplineValues[i, j] * localDisplacements[i, 0];
					knotDisplacements[paraviewKnotRenumbering[j], 1] += tsplines.TSplineValues[i, j] * localDisplacements[i, 1];
					knotDisplacements[paraviewKnotRenumbering[j], 2] += tsplines.TSplineValues[i, j] * localDisplacements[i, 2];
				}
			}

			return knotDisplacements;
		}

		public (double[,] knotStrains, double[,] knotStresses) CalculateStressesForPostProcessing(Element element, double[,] localDisplacements)
		{
			var tsplineElement = (TSplineKirchhoffLoveShellElementMaterial)element;
			var knotParametricCoordinatesKsi = new Vector(new double[] { -1, 1 });
			var knotParametricCoordinatesHeta = new Vector(new double[] { -1, 1 });

			ShapeTSplines2DFromBezierExtraction tsplines = new ShapeTSplines2DFromBezierExtraction(tsplineElement, tsplineElement.ControlPoints, knotParametricCoordinatesKsi, knotParametricCoordinatesHeta);

			var knotStresses = new double[4, 3];
			var knotStrains = new double[4, 3];

			for (int j = 0; j < knotStresses.GetLength(0); j++)
			{
				var jacobianMatrix = CalculateJacobian(tsplineElement, tsplines, j);

				var hessianMatrix = CalculateHessian(tsplineElement, tsplines, j);

				var surfaceBasisVector1 = CalculateSurfaceBasisVector1(jacobianMatrix, 0);

				var surfaceBasisVector2 = CalculateSurfaceBasisVector1(jacobianMatrix, 1);

				var surfaceBasisVector3 = new[]
				{
					surfaceBasisVector1[1] * surfaceBasisVector2[2] - surfaceBasisVector1[2] * surfaceBasisVector2[1],
					surfaceBasisVector1[2] * surfaceBasisVector2[0] - surfaceBasisVector1[0] * surfaceBasisVector2[2],
					surfaceBasisVector1[0] * surfaceBasisVector2[1] - surfaceBasisVector1[1] * surfaceBasisVector2[0]
				};

				double norm = 0;
				for (int i = 0; i < surfaceBasisVector3.Length; i++)
					norm += surfaceBasisVector3[i] * surfaceBasisVector3[i];
				var J1 = Math.Sqrt(norm);

				var surfaceBasisVectorDerivative1 = CalculateSurfaceBasisVector1(hessianMatrix, 0);
				var surfaceBasisVectorDerivative2 = CalculateSurfaceBasisVector1(hessianMatrix, 1);
				var surfaceBasisVectorDerivative12 = CalculateSurfaceBasisVector1(hessianMatrix, 2);

				var Bmembrane = CalculateMembraneDeformationMatrix(tsplines, j, surfaceBasisVector1, surfaceBasisVector2, tsplineElement);
				var Bbending = CalculateBendingDeformationMatrix(surfaceBasisVector3, tsplines, j, surfaceBasisVector2,
					surfaceBasisVectorDerivative1, surfaceBasisVector1, J1, surfaceBasisVectorDerivative2,
					surfaceBasisVectorDerivative12, tsplineElement);

				var membraneStrain = new double[Bmembrane.GetLength(0)];
				var bendingStrain = new double[Bmembrane.GetLength(0)];
				var ldisp = new double[localDisplacements.Length];
				var index = 0;
				for (int i = 0; i < localDisplacements.GetLength(0); i++)
				for (int k = 0; k < localDisplacements.GetLength(1); k++)
					ldisp[index++] = localDisplacements[i, k];

				for (int i = 0; i < Bmembrane.GetLength(0); i++)
				{
					for (int k = 0; k < Bmembrane.GetLength(1); k++)
					{
						membraneStrain[i] += Bmembrane[i, k] * ldisp[k];
						bendingStrain[i] += Bbending[i, k] * ldisp[k];
					}
				}

				var material = (tsplineElement.materialsAtThicknessGP.Values).ToList()[0].Values.ToList()[0].Clone();
				material.TangentVectorV1 = surfaceBasisVector1;
				material.TangentVectorV2 = surfaceBasisVector2;
				material.NormalVectorV3 = surfaceBasisVector3;
				var gpStrain = new double[bendingStrain.Length];
				var z = -Thickness/2;
				for (int i = 0; i < bendingStrain.Length; i++)
				{
					gpStrain[i] += membraneStrain[i] + bendingStrain[i] * z;
					knotStrains[j, i] += membraneStrain[i] + bendingStrain[i] * z;
				}

				material.UpdateMaterial(gpStrain);
				var stresses = material.Stresses;
				for (int i = 0; i < stresses.Length; i++)
				{
					knotStresses[j, i] = stresses[i];
				}
			}
			return (knotStrains, knotStresses);

		}

		public double[,] CalculatePointsForPostProcessing(TSplineKirchhoffLoveShellElementMaterial element)
		{
			var localCoordinates = new double[4, 2]
			{
				{-1, -1},
				{-1, 1},
				{1, -1},
				{1, 1}
			};

			var knotParametricCoordinatesKsi = new Vector(new double[] { -1, 1 });
			var knotParametricCoordinatesHeta = new Vector(new double[] { -1, 1 });

			ShapeTSplines2DFromBezierExtraction tsplines = new ShapeTSplines2DFromBezierExtraction(element, element.ControlPoints, knotParametricCoordinatesKsi, knotParametricCoordinatesHeta);

			var knotDisplacements = new double[4, 3];
			var paraviewKnotRenumbering = new int[] { 0, 3, 1, 2 };
			for (int j = 0; j < localCoordinates.GetLength(0); j++)
			{
				for (int i = 0; i < element.ControlPoints.Count; i++)
				{
					knotDisplacements[paraviewKnotRenumbering[j], 0] += tsplines.TSplineValues[i, j] * element.ControlPoints[i].X;
					knotDisplacements[paraviewKnotRenumbering[j], 1] += tsplines.TSplineValues[i, j] * element.ControlPoints[i].Y;
					knotDisplacements[paraviewKnotRenumbering[j], 2] += tsplines.TSplineValues[i, j] * element.ControlPoints[i].Z;
				}
			}

			return knotDisplacements;

		}

	}
}
