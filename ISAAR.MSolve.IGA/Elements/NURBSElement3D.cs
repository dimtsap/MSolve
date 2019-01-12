using System;
using System.Collections.Generic;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Discretization.Loads;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.IGA.Entities;
using ISAAR.MSolve.IGA.Entities.Loads;
using ISAAR.MSolve.IGA.Interfaces;
using ISAAR.MSolve.IGA.Problems.SupportiveClasses;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Materials.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using Element = ISAAR.MSolve.IGA.Entities.Element;

namespace ISAAR.MSolve.IGA.Elements
{
	public class NURBSElement3D : Element, IStructuralIsogeometricElement
	{
		protected readonly static DOFType[] controlPointDOFTypes = new DOFType[] { DOFType.X, DOFType.Y, DOFType.Z };
		protected DOFType[][] dofTypes;
		protected IElementDOFEnumerator dofEnumerator = new GenericDOFEnumerator();
		private DynamicMaterial dynamicProperties;


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

		public IList<IList<DOFType>> GetElementDOFTypes(IElement element)
		{
			var nurbsElement = (NURBSElement3D)element;
			dofTypes = new DOFType[nurbsElement.ControlPoints.Count][];
			for (int i = 0; i < nurbsElement.ControlPoints.Count; i++)
			{
				dofTypes[i] = controlPointDOFTypes;
			}
			return dofTypes;
		}

		public bool MaterialModified
		{
			get
			{
				throw new NotImplementedException();
			}
		}

		public void ResetMaterialModified()
		{
			throw new NotImplementedException();
		}

		public double[] CalculateForces(Element element, double[] localDisplacements, double[] localdDisplacements)
		{
			throw new NotImplementedException();
		}

		public double[] CalculateForcesForLogging(Element element, double[] localDisplacements)
		{
			throw new NotImplementedException();
		}

		public double[] CalculateAccelerationForces(FEM.Entities.Element element, IList<MassAccelerationLoad> loads)
		{
			throw new NotImplementedException();
		}

		public Tuple<double[], double[]> CalculateStresses(Element element, double[] localDisplacements, double[] localdDisplacements)
		{
			throw new NotImplementedException();
		}

		public void SaveMaterialState()
		{
			throw new NotImplementedException();
		}

		public void ClearMaterialState()
		{
			throw new NotImplementedException();
		}

		public void ClearMaterialStresses()
		{
			throw new NotImplementedException();
		}

		public IMatrix2D DampingMatrix(IElement element)
		{
			throw new NotImplementedException();
		}

		public IMatrix2D MassMatrix(IElement element)
		{
			throw new NotImplementedException();
		}

		public IMatrix2D StiffnessMatrix(IElement element)
		{
			var nurbsElement = (NURBSElement3D)element;
			IList<GaussLegendrePoint3D> gaussPoints = CreateElementGaussPoints(nurbsElement);
			Matrix2D stiffnessMatrixElement = new Matrix2D(nurbsElement.ControlPointsDictionary.Count * 3, nurbsElement.ControlPointsDictionary.Count * 3);

			NURBS3D nurbs = new NURBS3D(nurbsElement, nurbsElement.ControlPoints);

			for (int j = 0; j < gaussPoints.Count; j++)
			{
				Matrix2D jacobianMatrix = CalculateJacobian(nurbsElement.ControlPoints, nurbs, j);

				double jacdet = CalculateJacobianDeterminant(jacobianMatrix);

				Matrix2D inverseJacobian = CalculateInverseJacobian(jacobianMatrix, jacdet);

				Matrix2D B1 = CalculateDeformationMatrix1(inverseJacobian);

				Matrix2D B2 = CalculateDeformationMatrix2(nurbsElement.ControlPoints, nurbs, j);

				Matrix2D B = B1 * B2;

				Matrix2D stiffnessMatrixGaussPoint = B.Transpose() * ((IContinuumMaterial3D)nurbsElement.Patch.Material).ConstitutiveMatrix * B * jacdet * gaussPoints[j].WeightFactor;

				for (int m = 0; m < nurbsElement.ControlPoints.Count * 3; m++)
				{
					for (int n = 0; n < nurbsElement.ControlPoints.Count * 3; n++)
					{
						stiffnessMatrixElement[m, n] += stiffnessMatrixGaussPoint[m, n];
					}
				}
			}
			return stiffnessMatrixElement;

		}

		private static Matrix2D CalculateDeformationMatrix2(IList<ControlPoint> elementControlPoints, NURBS3D nurbs, int j)
		{
			Matrix2D B2 = new Matrix2D(9, 3 * elementControlPoints.Count);
			for (int column = 0; column < 3 * elementControlPoints.Count; column += 3)
			{
				B2[0, column] += nurbs.NurbsDerivativeValuesKsi[column / 3, j];
				B2[1, column] += nurbs.NurbsDerivativeValuesHeta[column / 3, j];
				B2[2, column] += nurbs.NurbsDerivativeValuesZeta[column / 3, j];

				B2[3, column + 1] += nurbs.NurbsDerivativeValuesKsi[column / 3, j];
				B2[4, column + 1] += nurbs.NurbsDerivativeValuesHeta[column / 3, j];
				B2[5, column + 1] += nurbs.NurbsDerivativeValuesZeta[column / 3, j];

				B2[6, column + 2] += nurbs.NurbsDerivativeValuesKsi[column / 3, j];
				B2[7, column + 2] += nurbs.NurbsDerivativeValuesHeta[column / 3, j];
				B2[8, column + 2] += nurbs.NurbsDerivativeValuesZeta[column / 3, j];
			}

			return B2;
		}

		private static Matrix2D CalculateDeformationMatrix1(Matrix2D inverseJacobian)
		{
			Matrix2D B1 = new Matrix2D(6, 9);

			B1[0, 0] += inverseJacobian[0, 0];
			B1[0, 1] += inverseJacobian[0, 1];
			B1[0, 2] += inverseJacobian[0, 2];

			B1[1, 3] += inverseJacobian[1, 0];
			B1[1, 4] += inverseJacobian[1, 1];
			B1[1, 5] += inverseJacobian[1, 2];

			B1[2, 6] += inverseJacobian[2, 0];
			B1[2, 7] += inverseJacobian[2, 1];
			B1[2, 8] += inverseJacobian[2, 2];

			B1[3, 0] += inverseJacobian[1, 0];
			B1[3, 1] += inverseJacobian[1, 1];
			B1[3, 2] += inverseJacobian[1, 2];
			B1[3, 3] += inverseJacobian[0, 0];
			B1[3, 4] += inverseJacobian[0, 1];
			B1[3, 5] += inverseJacobian[0, 2];

			B1[4, 3] += inverseJacobian[2, 0];
			B1[4, 4] += inverseJacobian[2, 1];
			B1[4, 5] += inverseJacobian[2, 2];
			B1[4, 6] += inverseJacobian[1, 0];
			B1[4, 7] += inverseJacobian[1, 1];
			B1[4, 8] += inverseJacobian[1, 2];

			B1[5, 0] += inverseJacobian[2, 0];
			B1[5, 1] += inverseJacobian[2, 1];
			B1[5, 2] += inverseJacobian[2, 2];
			B1[5, 6] += inverseJacobian[0, 0];
			B1[5, 7] += inverseJacobian[0, 1];
			B1[5, 8] += inverseJacobian[0, 2];
			return B1;
		}

		private static Matrix2D CalculateInverseJacobian(Matrix2D jacobianMatrix, double jacdet)
		{
			Matrix2D inverseJacobian = new Matrix2D(3, 3);

			inverseJacobian[0, 0] = jacobianMatrix[1, 1] * jacobianMatrix[2, 2] - jacobianMatrix[1, 2] * jacobianMatrix[2, 1];
			inverseJacobian[0, 1] = jacobianMatrix[0, 2] * jacobianMatrix[2, 1] - jacobianMatrix[0, 1] * jacobianMatrix[2, 2];
			inverseJacobian[0, 2] = jacobianMatrix[0, 1] * jacobianMatrix[1, 2] - jacobianMatrix[0, 2] * jacobianMatrix[1, 1];

			inverseJacobian[1, 0] = jacobianMatrix[1, 2] * jacobianMatrix[2, 0] - jacobianMatrix[1, 0] * jacobianMatrix[2, 2];
			inverseJacobian[1, 1] = jacobianMatrix[0, 0] * jacobianMatrix[2, 2] - jacobianMatrix[0, 2] * jacobianMatrix[2, 0];
			inverseJacobian[1, 2] = jacobianMatrix[0, 2] * jacobianMatrix[1, 0] - jacobianMatrix[0, 0] * jacobianMatrix[1, 2];

			inverseJacobian[2, 0] = jacobianMatrix[1, 0] * jacobianMatrix[2, 1] - jacobianMatrix[1, 1] * jacobianMatrix[2, 0];
			inverseJacobian[2, 1] = jacobianMatrix[0, 1] * jacobianMatrix[2, 0] - jacobianMatrix[0, 0] * jacobianMatrix[2, 1];
			inverseJacobian[2, 2] = jacobianMatrix[0, 0] * jacobianMatrix[1, 1] - jacobianMatrix[0, 1] * jacobianMatrix[1, 0];

			inverseJacobian = inverseJacobian * (1 / jacdet);
			return inverseJacobian;
		}

		private static double CalculateJacobianDeterminant(Matrix2D jacobianMatrix)
		{
			return jacobianMatrix[0, 0] * (jacobianMatrix[1, 1] * jacobianMatrix[2, 2] - jacobianMatrix[2, 1] * jacobianMatrix[1, 2])
								- jacobianMatrix[0, 1] * (jacobianMatrix[1, 0] * jacobianMatrix[2, 2] - jacobianMatrix[2, 0] * jacobianMatrix[1, 2])
								+ jacobianMatrix[0, 2] * (jacobianMatrix[1, 0] * jacobianMatrix[2, 1] - jacobianMatrix[2, 0] * jacobianMatrix[1, 1]);
		}

		private static Matrix2D CalculateJacobian(IList<ControlPoint> elementControlPoints, NURBS3D nurbs, int j)
		{
			Matrix2D jacobianMatrix = new Matrix2D(3, 3);

			for (int k = 0; k < elementControlPoints.Count; k++)
			{
				jacobianMatrix[0, 0] += nurbs.NurbsDerivativeValuesKsi[k, j] * elementControlPoints[k].X;
				jacobianMatrix[0, 1] += nurbs.NurbsDerivativeValuesKsi[k, j] * elementControlPoints[k].Y;
				jacobianMatrix[0, 2] += nurbs.NurbsDerivativeValuesKsi[k, j] * elementControlPoints[k].Z;
				jacobianMatrix[1, 0] += nurbs.NurbsDerivativeValuesHeta[k, j] * elementControlPoints[k].X;
				jacobianMatrix[1, 1] += nurbs.NurbsDerivativeValuesHeta[k, j] * elementControlPoints[k].Y;
				jacobianMatrix[1, 2] += nurbs.NurbsDerivativeValuesHeta[k, j] * elementControlPoints[k].Z;
				jacobianMatrix[2, 0] += nurbs.NurbsDerivativeValuesZeta[k, j] * elementControlPoints[k].X;
				jacobianMatrix[2, 1] += nurbs.NurbsDerivativeValuesZeta[k, j] * elementControlPoints[k].Y;
				jacobianMatrix[2, 2] += nurbs.NurbsDerivativeValuesZeta[k, j] * elementControlPoints[k].Z;
			}

			return jacobianMatrix;
		}

		private IList<GaussLegendrePoint3D> CreateElementGaussPoints(Element element)
		{
			GaussQuadrature gauss = new GaussQuadrature();
			return gauss.CalculateElementGaussPoints(element.Patch.DegreeKsi, element.Patch.DegreeHeta, element.Patch.DegreeZeta, element.Knots);
		}

		public Dictionary<int, double> CalculateLoadingCondition(Element element, Edge edge, NeumannBoundaryCondition neumann)
		{
			throw new NotSupportedException();
		}

		public Dictionary<int, double> CalculateLoadingCondition(Element element, Face face, NeumannBoundaryCondition neumann)
		{
			throw new NotSupportedException();
		}

		public Dictionary<int, double> CalculateLoadingCondition(Element element, Edge edge, PressureBoundaryCondition pressure)
		{
			throw new NotSupportedException();
		}

		public Dictionary<int, double> CalculateLoadingCondition(Element element, Face face, PressureBoundaryCondition pressure)
		{
			throw new NotSupportedException();
		}

		public double[,] CalculateDisplacementsForPostProcessing(Element element, double[,] localDisplacements)
		{
			var nurbsElement = (NURBSElement3D)element;
			var knotParametricCoordinatesKsi = new Vector(new double[] { element.Knots[0].Ksi, element.Knots[4].Ksi });
			var knotParametricCoordinatesHeta = new Vector(new double[] { element.Knots[0].Heta, element.Knots[2].Heta });
			var knotParametricCoordinates�eta = new Vector(new double[] { element.Knots[0].Zeta, element.Knots[1].Zeta });
			NURBS3D nurbs = new NURBS3D(nurbsElement, nurbsElement.ControlPoints, knotParametricCoordinatesKsi,
				knotParametricCoordinatesHeta, knotParametricCoordinates�eta);
			var knotDisplacements = new double[8, 3];
			var paraviewKnotRenumbering = new int[] { 0, 4, 2, 6, 1, 5, 3, 7 };
			for (int j = 0; j < element.Knots.Count; j++)
			{
				for (int i = 0; i < element.ControlPoints.Count; i++)
				{
					knotDisplacements[paraviewKnotRenumbering[j], 0] += nurbs.NurbsValues[i, j] * localDisplacements[i, 0];
					knotDisplacements[paraviewKnotRenumbering[j], 1] += nurbs.NurbsValues[i, j] * localDisplacements[i, 1];
					knotDisplacements[paraviewKnotRenumbering[j], 2] += nurbs.NurbsValues[i, j] * localDisplacements[i, 2];
				}
			}

			return knotDisplacements;
		}

		public (double[,] knotStrains, double[,] knotStresses) CalculateStressesForPostProcessing(Element element, double[,] localDisplacements)
		{
			throw new NotImplementedException();
		}

		public double[] CalculateAccelerationForces(Element element, IList<MassAccelerationLoad> loads)
		{
			throw new NotImplementedException();
		}
	}
}