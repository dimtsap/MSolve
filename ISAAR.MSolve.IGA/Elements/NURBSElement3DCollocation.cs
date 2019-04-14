﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Security.Cryptography.X509Certificates;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.IGA.Entities;
using ISAAR.MSolve.IGA.Entities.Loads;
using ISAAR.MSolve.IGA.Interfaces;
using ISAAR.MSolve.IGA.Problems.SupportiveClasses;
using ISAAR.MSolve.LinearAlgebra.Factorizations;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Solvers.LinearSystems;
using Vector = ISAAR.MSolve.LinearAlgebra.Vectors.Vector;

namespace ISAAR.MSolve.IGA.Elements
{
	/// <summary>
	/// Three dimensional collocation point.
	/// Calculated according to "Improving the Computational Performance of Isogeometric Analysis" of Gkritzalis, Christos.
	/// Authors: Dimitris Tsapetis.
	/// </summary>
	public class NURBSElement3DCollocation:Element, IStructuralIsogeometricElement
	{
		public NaturalPoint3D CollocationPoint;
		public ElementDimensions ElementDimensions => ElementDimensions.ThreeD;
		public IElementDofEnumerator_v2 DofEnumerator { get; set; }
		public bool MaterialModified { get; }

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

		public void ResetMaterialModified()
		{
			throw new NotImplementedException();
		}

		public Tuple<double[], double[]> CalculateStresses(Element element, double[] localDisplacements, double[] localdDisplacements)
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

		public double[,] CalculateDisplacementsForPostProcessing(Element element, double[,] localDisplacements)
		{
			throw new NotImplementedException();
		}

		public void ClearMaterialState()
		{
			throw new NotImplementedException();
		}

		public IMatrix StiffnessMatrix(IElement_v2 element)
		{
			var elementCollocation = (NURBSElement3DCollocation)element;

			var nurbs = new NURBS3D(elementCollocation.Patch.NumberOfControlPointsKsi,
				elementCollocation.Patch.NumberOfControlPointsHeta, elementCollocation.Patch.NumberOfControlPointsZeta,
				elementCollocation.Patch.DegreeKsi, elementCollocation.Patch.DegreeHeta,
				elementCollocation.Patch.DegreeZeta, elementCollocation.Patch.KnotValueVectorKsi,
				elementCollocation.Patch.KnotValueVectorHeta, elementCollocation.Patch.KnotValueVectorZeta,
				elementCollocation.ControlPoints.ToArray(), elementCollocation.CollocationPoint);

			var jacobianMatrix = CalculateJacobian(elementCollocation.ControlPoints.ToArray(), nurbs, 0);

			var hessianMatrix=CalculateHessian(elementCollocation.ControlPoints, nurbs);

			var squareDerivatives = CalculateSquareDerivatives(jacobianMatrix);

			var inverseJacobian = jacobianMatrix.Invert();

			var dR = CalculateNaturalDerivatives(nurbs, inverseJacobian);

			var ddR = CalculateNaturalSecondDerivatives(nurbs, hessianMatrix, dR, squareDerivatives);

			return CollocationPointStiffness(elementCollocation, ddR);
		}

		private static Matrix CollocationPointStiffness(NURBSElement3DCollocation elementCollocation, double[,] ddR)
		{
			var collocationPointStiffness = Matrix.CreateZero(3, elementCollocation.ControlPoints.Count * 3);

			var E = elementCollocation.Patch.Material.YoungModulus;
			var nu = elementCollocation.Patch.Material.PoissonRatio;

			var lambda = nu * E / (1 + nu) / (1 - 2 * nu);
			var m = E / 2 / (1 + nu);

			for (int i = 0; i < elementCollocation.ControlPoints.Count * 3; i += 3)
			{
				var index = i / 3;
				collocationPointStiffness[0, i] = (lambda + 2 * m) * ddR[0, index] + m * ddR[1, index] + m * ddR[2, index];
				collocationPointStiffness[1, i] = (lambda + m) * ddR[3, index];
				collocationPointStiffness[2, i] = (lambda + m) * ddR[4, index];

				collocationPointStiffness[0, i + 1] = (lambda + m) * ddR[3, index];
				collocationPointStiffness[1, i + 1] = (lambda + 2 * m) * ddR[1, index] + m * ddR[0, index] + m * ddR[2, index];
				collocationPointStiffness[2, i + 1] = (lambda + m) * ddR[5, index];

				collocationPointStiffness[0, i + 2] = (lambda + m) * ddR[4, index];
				collocationPointStiffness[1, i + 2] = (lambda + m) * ddR[5, index];
				collocationPointStiffness[2, i + 2] = (lambda + 2 * m) * ddR[2, 0] + m * ddR[0, index] + m * ddR[1, index];
			}

			return collocationPointStiffness;
		}

		private double[,] CalculateNaturalSecondDerivatives(NURBS3D nurbs, Matrix hessianMatrix, Matrix dR, Matrix squareDerivatives)
		{
			var ddR2 = hessianMatrix * dR;

			var ddR3 = new double[3, nurbs.NurbsSecondDerivativeValueKsi.GetLength(0)];
			for (int i = 0; i < ddR3.GetLength(1); i++)
			{
				ddR3[0, i] = nurbs.NurbsSecondDerivativeValueKsi[i, 0] - ddR2[0, i];
				ddR3[1, i] = nurbs.NurbsSecondDerivativeValueHeta[i, 0] - ddR2[1, i];
				ddR3[2, i] = nurbs.NurbsSecondDerivativeValueZeta[i, 0] - ddR2[2, i];

				ddR3[3, i] = nurbs.NurbsSecondDerivativeValueKsiHeta[i, 0] - ddR2[3, i];
				ddR3[4, i] = nurbs.NurbsSecondDerivativeValueKsiZeta[i, 0] - ddR2[4, i];
				ddR3[5, i] = nurbs.NurbsSecondDerivativeValueHetaZeta[i, 0] - ddR2[5, i];
			}

			var ddR = new double[3, nurbs.NurbsSecondDerivativeValueKsi.GetLength(0)];

			var LU=squareDerivatives.FactorLU();

			for (int i = 0; i < ddR.GetLength(1); i++)
			{
				var solution = Vector.CreateZero(6);
				var rhs = Vector.CreateFromArray(new[]
					{ddR3[0, i], ddR3[1, i], ddR3[2, i], ddR3[3, i], ddR3[4, i], ddR3[5, i]});
				LU.SolveLinearSystem(rhs, solution);
				ddR[0, i] = solution[0];
				ddR[1, i] = solution[1];
				ddR[2, i] = solution[2];
				ddR[3, i] = solution[3];
				ddR[4, i] = solution[4];
				ddR[5, i] = solution[5];
			}

			return ddR;
		}

		private Matrix CalculateNaturalDerivatives(NURBS3D nurbs, Matrix3by3 inverseJacobian)
		{
			var dR = Matrix.CreateZero(3, nurbs.NurbsDerivativeValuesKsi.GetLength(0));
			for (int i = 0; i < dR.NumColumns; i++)
			{
				var dKsi = nurbs.NurbsDerivativeValuesKsi[i, 0];
				var dHeta = nurbs.NurbsDerivativeValuesHeta[i, 0];
				var dZeta = nurbs.NurbsDerivativeValuesZeta[i, 0];

				dR[0, i] = inverseJacobian[0, 0] * dKsi + inverseJacobian[0, 1] * dHeta + inverseJacobian[0, 2] * dZeta;
				dR[1, i] = inverseJacobian[1, 0] * dKsi + inverseJacobian[1, 1] * dHeta + inverseJacobian[1, 2] * dZeta;
				dR[2, i] = inverseJacobian[2, 0] * dKsi + inverseJacobian[2, 1] * dHeta + inverseJacobian[2, 2] * dZeta;
			}

			return dR;
		}

		private Matrix CalculateSquareDerivatives(Matrix3by3 jacobianMatrix)
		{
			var J11 = jacobianMatrix[0, 0];
			var J12 = jacobianMatrix[0, 1];
			var J13 = jacobianMatrix[0, 2];
			var J21 = jacobianMatrix[1, 0];
			var J22 = jacobianMatrix[1, 1];
			var J23 = jacobianMatrix[1, 2];
			var J31 = jacobianMatrix[2, 0];
			var J32 = jacobianMatrix[2, 1];
			var J33 = jacobianMatrix[2, 2];


			var squareDerivatives = Matrix.CreateZero(6, 6);
			squareDerivatives[0, 0] = J11*J11;
			squareDerivatives[0, 1] = J12*J12;
			squareDerivatives[0, 2] = J13*J13;
			squareDerivatives[0, 3] = 2*J11*J12;
			squareDerivatives[0, 4] = 2*J11*J13;
			squareDerivatives[0, 5] = 2*J12*J13;

			squareDerivatives[1, 0] = J21*J21;
			squareDerivatives[1, 1] = J22*J22;
			squareDerivatives[1, 2] = J23*J23;
			squareDerivatives[1, 3] = 2*J21*J22;
			squareDerivatives[1, 4] = 2*J21*J23;
			squareDerivatives[1, 5] = 2*J22*J23;

			squareDerivatives[2, 0] = J31*J31;
			squareDerivatives[2, 1] = J32*J32;
			squareDerivatives[2, 2] = J33*J33;
			squareDerivatives[2, 3] = 2*J31*J32;
			squareDerivatives[2, 4] = 2*J31*J33;
			squareDerivatives[2, 5] = 2*J32*J33;

			squareDerivatives[3, 0] = J11*J21;
			squareDerivatives[3, 1] = J12*J22;
			squareDerivatives[3, 2] = J13*J23;
			squareDerivatives[3, 3] = J11*J22+J21*J12;
			squareDerivatives[3, 4] = J11*J23+J21*J13;
			squareDerivatives[3, 5] = J12*J23+J22*J13;

			squareDerivatives[4, 0] = J11*J31;
			squareDerivatives[4, 1] = J12*J32;
			squareDerivatives[4, 2] = J13*J33;
			squareDerivatives[4, 3] = J11*J32+J31*J12;
			squareDerivatives[4, 4] = J11*J33+J31*J13;
			squareDerivatives[4, 5] = J12*J33+J32*J13;

			squareDerivatives[5, 0] = J21*J31;
			squareDerivatives[5, 1] = J22*J32;
			squareDerivatives[5, 2] = J23*J33;
			squareDerivatives[5, 3] = J21*J32+J31*J22;
			squareDerivatives[5, 4] = J21*J33+J31*J23;
			squareDerivatives[5, 5] = J22*J33+J32*J23;

			return squareDerivatives;
		}

		private Matrix CalculateHessian(IList<ControlPoint> controlPoints, NURBS3D nurbs)
		{
			var hessianMatrix = Matrix.CreateZero(6,3);
			for (int k = 0; k < controlPoints.Count; k++)
			{
				hessianMatrix[0, 0] += nurbs.NurbsSecondDerivativeValueKsi[k, 0] * ControlPoints[k].X;
				hessianMatrix[0, 1] += nurbs.NurbsSecondDerivativeValueKsi[k, 0] * ControlPoints[k].Y;
				hessianMatrix[0, 2] += nurbs.NurbsSecondDerivativeValueKsi[k, 0] * ControlPoints[k].Z;

				hessianMatrix[1, 0] += nurbs.NurbsSecondDerivativeValueHeta[k, 0] * ControlPoints[k].X;
				hessianMatrix[1, 1] += nurbs.NurbsSecondDerivativeValueHeta[k, 0] * ControlPoints[k].Y;
				hessianMatrix[1, 2] += nurbs.NurbsSecondDerivativeValueHeta[k, 0] * ControlPoints[k].Z;

				hessianMatrix[2, 0] += nurbs.NurbsSecondDerivativeValueZeta[k, 0] * ControlPoints[k].X;
				hessianMatrix[2, 1] += nurbs.NurbsSecondDerivativeValueZeta[k, 0] * ControlPoints[k].Y;
				hessianMatrix[2, 2] += nurbs.NurbsSecondDerivativeValueZeta[k, 0] * ControlPoints[k].Z;

				hessianMatrix[3, 0] += nurbs.NurbsSecondDerivativeValueKsiHeta[k, 0] * ControlPoints[k].X;
				hessianMatrix[3, 1] += nurbs.NurbsSecondDerivativeValueKsiHeta[k, 0] * ControlPoints[k].Y;
				hessianMatrix[3, 2] += nurbs.NurbsSecondDerivativeValueKsiHeta[k, 0] * ControlPoints[k].Z;

				hessianMatrix[4, 0] += nurbs.NurbsSecondDerivativeValueKsiZeta[k, 0] * ControlPoints[k].X;
				hessianMatrix[4, 1] += nurbs.NurbsSecondDerivativeValueKsiZeta[k, 0] * ControlPoints[k].Y;
				hessianMatrix[4, 2] += nurbs.NurbsSecondDerivativeValueKsiZeta[k, 0] * ControlPoints[k].Z;

				hessianMatrix[5, 0] += nurbs.NurbsSecondDerivativeValueHetaZeta[k, 0] * ControlPoints[k].X;
				hessianMatrix[5, 1] += nurbs.NurbsSecondDerivativeValueHetaZeta[k, 0] * ControlPoints[k].Y;
				hessianMatrix[5, 2] += nurbs.NurbsSecondDerivativeValueHetaZeta[k, 0] * ControlPoints[k].Z;
			}

			return hessianMatrix;
		}

		private Matrix3by3 CalculateJacobian(ControlPoint[] elementControlPoints, NURBS3D nurbs, int j)
		{
			var jacobianMatrix = Matrix3by3.CreateZero();

			for (int k = 0; k < elementControlPoints.Length; k++)
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

		public IMatrix MassMatrix(IElement_v2 element)
		{
			throw new NotImplementedException();
		}

		public IMatrix DampingMatrix(IElement_v2 element)
		{
			throw new NotImplementedException();
		}

		public IList<IList<DOFType>> GetElementDOFTypes(IElement_v2 element)
		{
			throw new NotImplementedException();
		}
	}
}