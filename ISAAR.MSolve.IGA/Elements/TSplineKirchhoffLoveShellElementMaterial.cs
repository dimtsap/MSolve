using System;
using System.Collections.Generic;
using System.Linq;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Discretization.Mesh;
using ISAAR.MSolve.IGA.Entities;
using ISAAR.MSolve.IGA.Entities.Loads;
using ISAAR.MSolve.IGA.Interfaces;
using ISAAR.MSolve.IGA.Problems.SupportiveClasses;
using ISAAR.MSolve.IGA.SupportiveClasses;
using ISAAR.MSolve.LinearAlgebra.Input;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Materials.Interfaces;

namespace ISAAR.MSolve.IGA.Elements
{
    /// <summary>
    /// Thesis Kiendl
    /// </summary>
    public class TSplineKirchhoffLoveShellElementMaterial : Element, IStructuralIsogeometricElement
    {
        public Matrix ExtractionOperator { get; set; }
        public int DegreeKsi { get; set; }
        public int DegreeHeta { get; set; }

        protected readonly static IDofType[] controlPointDOFTypes = new IDofType[]
            {StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ};

        protected IDofType[][] dofTypes;
        protected IElementDofEnumerator dofEnumerator = new GenericDofEnumerator();
        private DynamicMaterial dynamicProperties;
        private IReadOnlyList<IShellMaterial> materialsAtGaussPoints;

        private Dictionary<GaussLegendrePoint3D, Dictionary<GaussLegendrePoint3D, IShellMaterial>>
            materialsAtThicknessGP =
                new Dictionary<GaussLegendrePoint3D, Dictionary<GaussLegendrePoint3D, IShellMaterial>>();

        public CellType CellType { get; } = CellType.Unknown;

        public TSplineKirchhoffLoveShellElementMaterial(int id, Patch patch, int degreeKsi, int degreeHeta,
            double thickness, Matrix extractionOperator, IShellMaterial shellMaterial)
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

        public IElementDofEnumerator DofEnumerator
        {
            get { return dofEnumerator; }

            set { this.dofEnumerator = value; }
        }

        public ElementDimensions ElementDimensions
        {
            get { return ElementDimensions.ThreeD; }
        }

        public double Thickness { get; set; }

        public bool MaterialModified => throw new NotImplementedException();

        public double[] CalculateForces(Element element, double[] localDisplacements, double[] localdDisplacements)
        {
            var shellElement = (TSplineKirchhoffLoveShellElementMaterial) element;
            IList<GaussLegendrePoint3D> gaussPoints = CreateElementGaussPoints(shellElement);
            var ElementNodalForces = new double[shellElement.ControlPointsDictionary.Count * 3];
            var tsplines = new ShapeTSplines2DFromBezierExtraction(shellElement, shellElement.ControlPoints);


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

                var Bmembrane = CalculateMembraneDeformationMatrix(tsplines, j, surfaceBasisVector1,
                    surfaceBasisVector2, shellElement);
                var Bbending = CalculateBendingDeformationMatrix(surfaceBasisVector3, tsplines, j, surfaceBasisVector2,
                    surfaceBasisVectorDerivative1, surfaceBasisVector1, J1, surfaceBasisVectorDerivative2,
                    surfaceBasisVectorDerivative12, shellElement);

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

            return ElementNodalForces;
        }

        public double[] CalculateForcesForLogging(Element element, double[] localDisplacements)
        {
            throw new NotImplementedException();
        }

        public void SaveMaterialState()
        {
            throw new NotImplementedException();
        }

        public Dictionary<int, double> CalculateLoadingCondition(Element element, Edge edge,
            NeumannBoundaryCondition neumann)
        {
            throw new NotImplementedException();
        }

        public Dictionary<int, double> CalculateLoadingCondition(Element element, Face face,
            NeumannBoundaryCondition neumann)
        {
            throw new NotImplementedException();
        }

        public Dictionary<int, double> CalculateLoadingCondition(Element element, Edge edge,
            PressureBoundaryCondition pressure)
        {
            throw new NotImplementedException();
        }

        public Dictionary<int, double> CalculateLoadingCondition(Element element, Face face,
            PressureBoundaryCondition pressure)
        {
            throw new NotImplementedException();
        }

        public Tuple<double[], double[]> CalculateStresses(Element element, double[] localDisplacements,
            double[] localdDisplacements)
        {
            var shellElement = (TSplineKirchhoffLoveShellElementMaterial) element;
            //Matrix stiffnessMatrixElement = new Matrix(shellElement.ControlPointsDictionary.Count * 3, shellElement.ControlPointsDictionary.Count * 3);

            ShapeTSplines2DFromBezierExtraction tsplines =
                new ShapeTSplines2DFromBezierExtraction(shellElement, shellElement.ControlPoints);

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

                var Bmembrane = CalculateMembraneDeformationMatrix(tsplines, j, surfaceBasisVector1,
                    surfaceBasisVector2, shellElement);
                var Bbending = CalculateBendingDeformationMatrix(surfaceBasisVector3, tsplines, j, surfaceBasisVector2,
                    surfaceBasisVectorDerivative1, surfaceBasisVector1, J1, surfaceBasisVectorDerivative2,
                    surfaceBasisVectorDerivative12, shellElement);

                var membraneStrain = new double[Bmembrane.GetLength(0)];
                var bendingStrain = new double[Bmembrane.GetLength(0)];
                for (int i = 0; i < Bmembrane.GetLength(1); i++)
                {
                    for (int k = 0; k < Bmembrane.GetLength(0); k++)
                    {
                        membraneStrain[i] += Bmembrane[k, i] * localDisplacements[k];
                        bendingStrain[i] += Bbending[k, i] * localDisplacements[k];
                    }
                }

                foreach (var keyValuePair in materialsAtThicknessGP[materialsAtThicknessGP.Keys.ToList()[j]])
                {
                    var thicknessPoint = keyValuePair.Key;
                    var material = keyValuePair.Value;
                    var gpStrain = new double[bendingStrain.Length];
                    var z = -thicknessPoint.Zeta;
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

        public void ClearMaterialStresses()
        {
            throw new NotImplementedException();
        }

        public IMatrix DampingMatrix(IElement element)
        {
            throw new NotImplementedException();
        }

        public IReadOnlyList<IReadOnlyList<IDofType>> GetElementDofTypes(IElement element)
        {
            var nurbsElement = (TSplineKirchhoffLoveShellElementMaterial) element;
            dofTypes = new IDofType[nurbsElement.ControlPoints.Count][];
            for (int i = 0; i < nurbsElement.ControlPoints.Count; i++)
            {
                dofTypes[i] = controlPointDOFTypes;
            }

            return dofTypes;
        }

        public IMatrix MassMatrix(IElement element)
        {
            throw new NotImplementedException();
        }

        public void ResetMaterialModified()
        {
            throw new NotImplementedException();
        }

        //public IMatrix StiffnessMatrix(IElement element)
        //{
        //	var shellElement = (TSplineKirchhoffLoveShellElementMaterial)element;
        //	Matrix stiffnessMatrixElement = new Matrix(shellElement.ControlPointsDictionary.Count * 3, shellElement.ControlPointsDictionary.Count * 3);
        //	var gaussPoints = materialsAtThicknessGP.Keys.ToArray();
        //	ShapeTSplines2DFromBezierExtraction tsplines = new ShapeTSplines2DFromBezierExtraction(shellElement, shellElement.ControlPoints);


        //	var bRows = 3;
        //	var bCols = shellElement.ControlPoints.Count * 3;
        //	var stiffnessMatrix = new double[bCols, bCols];
        //	var BmTranspose = new double[bCols, bRows];
        //	var BbTranspose = new double[bCols, bRows];
        //	var BmbTranspose = new double[bCols, bRows];
        //	var BbmTranspose = new double[bCols, bRows];

        //	var BmTransposeMultStiffness = new double[bCols, bRows];
        //	var BbTransposeMultStiffness = new double[bCols, bRows];
        //	var BmbTransposeMultStiffness = new double[bCols, bRows];
        //	var BbmTransposeMultStiffness = new double[bCols, bRows];

        //	for (int j = 0; j < gaussPoints.Length; j++)
        //	{
        //		var jacobianMatrix = CalculateJacobian(shellElement, tsplines, j);

        //		var hessianMatrix = CalculateHessian(shellElement, tsplines, j);
        //		var surfaceBasisVector1 = CalculateSurfaceBasisVector1(jacobianMatrix, 0);

        //		var surfaceBasisVector2 = CalculateSurfaceBasisVector1(jacobianMatrix, 1);

        //		var surfaceBasisVector3 = surfaceBasisVector1^surfaceBasisVector2;

        //		foreach (var integrationPointMaterial in materialsAtThicknessGP[gaussPoints[j]].Values)
        //		{
        //			integrationPointMaterial.TangentVectorV1 = surfaceBasisVector1.Data;
        //			integrationPointMaterial.TangentVectorV2 = surfaceBasisVector2.Data;
        //			integrationPointMaterial.NormalVectorV3 = surfaceBasisVector3.Data;
        //		}
        //		var J1 = surfaceBasisVector3.Norm;
        //		surfaceBasisVector3.ScaleIntoThis(1 / J1);

        //		var surfaceBasisVectorDerivative1 = CalculateSurfaceBasisVector1(hessianMatrix, 0);
        //		var surfaceBasisVectorDerivative2 = CalculateSurfaceBasisVector1(hessianMatrix, 1);
        //		var surfaceBasisVectorDerivative12 = CalculateSurfaceBasisVector1(hessianMatrix, 2);

        //		var Bmembrane = CalculateMembraneDeformationMatrix(tsplines, j, surfaceBasisVector1, surfaceBasisVector2, shellElement);
        //		var Bbending = CalculateBendingDeformationMatrix(surfaceBasisVector3, tsplines, j, surfaceBasisVector2,
        //			surfaceBasisVectorDerivative1, surfaceBasisVector1, J1, surfaceBasisVectorDerivative2,
        //			surfaceBasisVectorDerivative12, shellElement);

        //		var (MembraneConstitutiveMatrix, BendingConstitutiveMatrix, CouplingConstitutiveMatrix) =
        //			IntegratedConstitutiveOverThickness(gaussPoints, j);

        //		var Kmembrane = Bmembrane.Transpose() * MembraneConstitutiveMatrix * Bmembrane  * J1 *
        //		                gaussPoints[j].WeightFactor;
        //		var Kbending = Bbending.Transpose() * BendingConstitutiveMatrix * Bbending  * J1 *
        //					   gaussPoints[j].WeightFactor;

        //		var KMembraneBending = Bmembrane.Transpose() * CouplingConstitutiveMatrix * Bbending * J1 *
        //		                       gaussPoints[j].WeightFactor;

        //		var KBendingMembrane = Bbending.Transpose() * CouplingConstitutiveMatrix * Bmembrane * J1 *
        //		                       gaussPoints[j].WeightFactor;


        //		stiffnessMatrixElement.AddIntoThis(Kmembrane);
        //		stiffnessMatrixElement.AddIntoThis(Kbending);
        //		stiffnessMatrixElement.AddIntoThis(KMembraneBending);
        //		stiffnessMatrixElement.AddIntoThis(KBendingMembrane);
        //	}
        //	return stiffnessMatrixElement;
        //}

        private double[][] initialSurfaceBasisVectors1;
        private double[][] initialSurfaceBasisVectors2;
        private double[][] initialSurfaceBasisVectors3;

        private double[][] initialSurfaceBasisVectorDerivative1;
        private double[][] initialSurfaceBasisVectorDerivative2;
        private double[][] initialSurfaceBasisVectorDerivative12;

        private bool isInitialized;

        public IMatrix StiffnessMatrix(IElement element)
        {
            var shellElement = (TSplineKirchhoffLoveShellElementMaterial) element;
            var gaussPoints = materialsAtThicknessGP.Keys.ToArray();
            ShapeTSplines2DFromBezierExtraction tsplines =
                new ShapeTSplines2DFromBezierExtraction(shellElement, shellElement.ControlPoints);

            if (!isInitialized)
            {
                CalculateInitialConfigurationData(shellElement, tsplines, gaussPoints);
                isInitialized = true;
            }


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

                var Bmembrane = CalculateMembraneDeformationMatrix(tsplines, j, surfaceBasisVector1,
                    surfaceBasisVector2, shellElement);
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

                var (MembraneForces, BendingMoments) = IntegratedStressesOverThickness(gaussPoints, j);


                var KmembraneNL = CalculateKmembraneNL(shellElement, MembraneForces, tsplines, j);
                // TODO: add to element stiffness Matrix


                var KbendingNL = CalculateKbendingNL(shellElement, BendingMoments, tsplines,
                    Vector.CreateFromArray(surfaceBasisVector1), Vector.CreateFromArray(surfaceBasisVector2),
                    Vector.CreateFromArray(surfaceBasisVector3),
                    Vector.CreateFromArray(surfaceBasisVectorDerivative1),
                    Vector.CreateFromArray(surfaceBasisVectorDerivative2),
                    Vector.CreateFromArray(surfaceBasisVectorDerivative12), J1, j);
            }

            return Matrix.CreateFromArray(stiffnessMatrix);
        }

        private IMatrix CalculateKbendingNL(TSplineKirchhoffLoveShellElementMaterial shellElement,
            double[] bendingMoments, ShapeTSplines2DFromBezierExtraction tsplines, Vector surfaceBasisVector1,
            Vector surfaceBasisVector2, Vector surfaceBasisVector3, Vector surfaceBasisVectorDerivative1, Vector surfaceBasisVectorDerivative2,
            Vector surfaceBasisVectorDerivative12,double J1, int j)
        {
            var KbendingNL =
                Matrix.CreateZero(shellElement.ControlPoints.Count * 3, shellElement.ControlPoints.Count * 3);


            for (int i = 0; i < shellElement.ControlPoints.Count; i++)
            {
                var a1r = Matrix3by3.CreateIdentity().Scale(tsplines.TSplineDerivativeValuesKsi[i, j]);
                var a2r = Matrix3by3.CreateIdentity().Scale(tsplines.TSplineDerivativeValuesHeta[i, j]);

                var a11r = Matrix3by3.CreateIdentity().Scale(tsplines.TSplineSecondDerivativesValueKsi[i, j]);
                var a22r = Matrix3by3.CreateIdentity().Scale(tsplines.TSplineSecondDerivativesValueHeta[i, j]);
                var a12r = Matrix3by3.CreateIdentity().Scale(tsplines.TSplineSecondDerivativesValueKsiHeta[i, j]);
                for (int k = 0; k < shellElement.ControlPoints.Count; k++)
                {
                    var a11s = Matrix3by3.CreateIdentity().Scale(tsplines.TSplineSecondDerivativesValueKsi[k, j]);
                    var a22s = Matrix3by3.CreateIdentity().Scale(tsplines.TSplineSecondDerivativesValueHeta[k, j]);
                    var a12s = Matrix3by3.CreateIdentity().Scale(tsplines.TSplineSecondDerivativesValueKsiHeta[k, j]);

                    var a3r = CalculateA3r(tsplines, i, j, surfaceBasisVector2, surfaceBasisVector1);
                    var a3s = CalculateA3r(tsplines, k, j, surfaceBasisVector2, surfaceBasisVector1);

                    var a1s = Matrix3by3.CreateIdentity().Scale(tsplines.TSplineDerivativeValuesKsi[k, j]);
                    var a2s = Matrix3by3.CreateIdentity().Scale(tsplines.TSplineDerivativeValuesHeta[k, j]);

                    #region B
                    var term1_532 = new Vector[3, 3];
                    for (int m = 0; m < 3; m++)
                    {
                        for (int n = 0; n < 3; n++)
                        {
                            var temp = a1r.GetRow(m).CrossProduct(a2s.GetRow(n));
                            temp.ScaleIntoThis(J1);
                            term1_532[m, n] = temp;
                        }
                    }


                    var term2_532 = new Vector[3, 3];
                    for (int m = 0; m < 3; m++)
                    {
                        var a3r_dashed = a1r.GetRow(m).CrossProduct(surfaceBasisVector2) +
                                         surfaceBasisVector1.CrossProduct(a2r.GetRow(m));
                        for (int n = 0; n < 3; n++)
                        {
                            //TODO: a3s_dashed, a3r_dashed calculated out of the loop for all cp 
                            var a3s_dashed = a1s.GetRow(n).CrossProduct(surfaceBasisVector2) +
                                             surfaceBasisVector1.CrossProduct(a2s.GetRow(n));
                            var term_525 = surfaceBasisVector3 * a3s_dashed;
                            term2_532[m, n] = a3r_dashed.Scale(-term_525 / J1 / J1);
                        }
                    }


                    var term3_532 = new Vector[3, 3];
                    for (int m = 0; m < 3; m++)
                    {
                        var a3r_dashed = a1r.GetRow(m).CrossProduct(surfaceBasisVector2) +
                                         surfaceBasisVector1.CrossProduct(a2r.GetRow(m));
                        for (int n = 0; n < 3; n++)
                        {
                            var a3s_dashed = a1s.GetRow(n).CrossProduct(surfaceBasisVector2) +
                                             surfaceBasisVector1.CrossProduct(a2s.GetRow(n));
                            var term_525 = surfaceBasisVector3 * a3r_dashed;
                            term3_532[m, n] = a3s_dashed.Scale(-term_525 / J1 / J1);
                        }
                    }


                    var term4_532 = new Vector[3, 3];
                    for (int m = 0; m < 3; m++)
                    {
                        var a3r_dashed = a1r.GetRow(m).CrossProduct(surfaceBasisVector2) +
                                         surfaceBasisVector1.CrossProduct(a2r.GetRow(m));
                        for (int n = 0; n < 3; n++)
                        {
                            var a3s_dashed = a1s.GetRow(n).CrossProduct(surfaceBasisVector2) +
                                             surfaceBasisVector1.CrossProduct(a2s.GetRow(n));
                            // term 5_31
                            var a3_rs = term1_532[m, n] * surfaceBasisVector3 * J1 + a3r_dashed * a3s_dashed / J1 -
                                        (a3r_dashed * surfaceBasisVector3) * (a3s_dashed * surfaceBasisVector3) / J1;
                            term4_532[m, n] = surfaceBasisVector3.Scale(-a3_rs / J1);
                        }
                    }

                    var term5_532 = new Vector[3, 3];
                    for (int m = 0; m < 3; m++)
                    {
                        var a3r_dashed = a1r.GetRow(m).CrossProduct(surfaceBasisVector2) +
                                         surfaceBasisVector1.CrossProduct(a2r.GetRow(m));
                        var term_525_r = surfaceBasisVector3 * a3r_dashed;
                        for (int n = 0; n < 3; n++)
                        {
                            var a3s_dashed = a1s.GetRow(n).CrossProduct(surfaceBasisVector2) +
                                             surfaceBasisVector1.CrossProduct(a2s.GetRow(n));
                            var term_525_s = surfaceBasisVector3 * a3s_dashed;
                            term5_532[m, n] = surfaceBasisVector3.Scale(2 / J1 / J1 * term_525_r * term_525_s);
                        }
                    }

                    var a3rs = new Vector[3, 3];
                    for (int m = 0; m < 3; m++)
                    {
                        for (int n = 0; n < 3; n++)
                        {
                            a3rs[m, n] = term1_532[m, n] + term2_532[m, n] + term3_532[m, n] + term4_532[m, n] +
                                         term5_532[m, n];
                        }
                    }


                    #endregion

                    var termA = bendingMoments[0] * (a11r * a3s + a11s * a3r) +
                                bendingMoments[2] * (a22r * a3s + a22s * a3r) +
                                bendingMoments[3] * (a12r * a3s + a12s * a3r) * 2;

                    var aux = bendingMoments[0] * surfaceBasisVectorDerivative1 +
                               bendingMoments[1] * surfaceBasisVectorDerivative2 +
                               2 * bendingMoments[2] * surfaceBasisVectorDerivative12;
                    var termB = Matrix3by3.CreateZero();
                    for (int m = 0; m < 3; m++)
                    {
                        for (int n = 0; n < 3; n++)
                        {
                            termB[m, n] = aux * a3rs[m, n];
                        }
                    }

                    var gaussPointStiffness = termA + termB;
                    for (int l = i * 3; l < i * 3 + 3; l++)
                    for (int m = k * 3; m < k * 3 + 3; m++)
                        KbendingNL[l, m] += gaussPointStiffness[l / 3, m / 3];
                }
            }

            return KbendingNL;
        }


        private Matrix3by3 CalculateA3r(ShapeTSplines2DFromBezierExtraction tsplines, int i, int j,
            Vector surfaceBasisVector2, Vector surfaceBasisVector1)
        {
            var aux1 = Vector.CreateFromArray(new double[] {tsplines.TSplineDerivativeValuesKsi[i, j], 0, 0})
                .CrossProduct(surfaceBasisVector2);
            var aux2 = Vector.CreateFromArray(new double[] {0, tsplines.TSplineDerivativeValuesKsi[i, j], 0})
                .CrossProduct(surfaceBasisVector2);
            var aux3 = Vector.CreateFromArray(new double[] {0, 0, tsplines.TSplineDerivativeValuesKsi[i, j]})
                .CrossProduct(surfaceBasisVector2);

            var aux4 = surfaceBasisVector1.CrossProduct(Vector.CreateFromArray(new double[]
                {tsplines.TSplineDerivativeValuesHeta[i, j], 0, 0}));
            var aux5 = surfaceBasisVector1.CrossProduct(Vector.CreateFromArray(new double[]
                {0, tsplines.TSplineDerivativeValuesHeta[i, j], 0}));
            var aux6 = surfaceBasisVector1.CrossProduct(Vector.CreateFromArray(new double[]
                {0, 0, tsplines.TSplineDerivativeValuesHeta[i, j]}));

            var a3r = Matrix3by3.CreateZero();
            a3r[0, 0] = aux1[0] + aux4[0];
            a3r[0, 0] = aux1[1] + aux4[1];
            a3r[0, 0] = aux1[2] + aux4[2];

            a3r[0, 0] = aux2[0] + aux5[0];
            a3r[0, 0] = aux2[1] + aux5[1];
            a3r[0, 0] = aux2[2] + aux5[2];

            a3r[0, 0] = aux3[0] + aux6[0];
            a3r[0, 0] = aux3[1] + aux6[1];
            a3r[0, 0] = aux3[2] + aux6[2];
            return a3r;
        }

        private void CalculateInitialConfigurationData(TSplineKirchhoffLoveShellElementMaterial shellElement,
            ShapeTSplines2DFromBezierExtraction tsplines, IList<GaussLegendrePoint3D> gaussPoints)
        {
            var numberOfGP = gaussPoints.Count;
            initialSurfaceBasisVectors1 = new double[numberOfGP][];
            initialSurfaceBasisVectors2 = new double[numberOfGP][];
            initialSurfaceBasisVectors3 = new double[numberOfGP][];
            initialSurfaceBasisVectorDerivative1 = new double[numberOfGP][];
            initialSurfaceBasisVectorDerivative2 = new double[numberOfGP][];
            initialSurfaceBasisVectorDerivative12 = new double[numberOfGP][];

            for (int j = 0; j < gaussPoints.Count; j++)
            {
                var jacobianMatrix = CalculateJacobian(shellElement, tsplines, j);

                var hessianMatrix = CalculateHessian(shellElement, tsplines, j);
                initialSurfaceBasisVectors1[j] = CalculateSurfaceBasisVector1(jacobianMatrix, 0);
                initialSurfaceBasisVectors2[j] = CalculateSurfaceBasisVector1(jacobianMatrix, 1);
                initialSurfaceBasisVectors3[j] =
                    CalculateCrossProduct(initialSurfaceBasisVectors1[j], initialSurfaceBasisVectors2[j]);

                initialSurfaceBasisVectorDerivative1[j] = CalculateSurfaceBasisVector1(hessianMatrix, 0);
                initialSurfaceBasisVectorDerivative2[j] = CalculateSurfaceBasisVector1(hessianMatrix, 1);
                initialSurfaceBasisVectorDerivative12[j] = CalculateSurfaceBasisVector1(hessianMatrix, 2);
            }
        }

        private double[] CalculateCrossProduct(double[] vector1, double[] vector2)
        {
            return new[]
            {
                vector1[1] * vector2[2] - vector1[2] * vector2[1],
                vector1[2] * vector2[0] - vector1[0] * vector2[2],
                vector1[0] * vector2[1] - vector1[1] * vector2[0]
            };
        }


        private Matrix CalculateKmembraneNL(TSplineKirchhoffLoveShellElementMaterial shellElement,
            double[] membraneForces, ShapeTSplines2DFromBezierExtraction tSplines, int j)
        {
            var KmembraneNL =
                Matrix.CreateZero(shellElement.ControlPoints.Count * 3, shellElement.ControlPoints.Count * 3);

            for (int i = 0; i < shellElement.ControlPoints.Count; i++)
            {
                var a1r = Matrix3by3.CreateIdentity().Scale(tSplines.TSplineDerivativeValuesKsi[i, j]);
                var a2r = Matrix3by3.CreateIdentity().Scale(tSplines.TSplineDerivativeValuesHeta[i, j]);
                for (int k = 0; k < shellElement.ControlPoints.Count; k++)
                {
                    var a1s = Matrix3by3.CreateIdentity().Scale(tSplines.TSplineDerivativeValuesKsi[k, j]);
                    var a2s = Matrix3by3.CreateIdentity().Scale(tSplines.TSplineDerivativeValuesHeta[k, j]);

                    var klocal = membraneForces[0] * a1r * a1s + membraneForces[1] * a2r * a2s +
                                 membraneForces[2] * (a1r * a2s + a1s * a2r);

                    for (int l = i * 3; l < i * 3 + 3; l++)
                    for (int m = k * 3; m < k * 3 + 3; m++)
                        KmembraneNL[l, m] += klocal[l / 3, m / 3];
                }
            }

            return KmembraneNL;
        }


        public (double[,] MembraneConstitutiveMatrix, double[,] BendingConstitutiveMatrix, double[,]
            CouplingConstitutiveMatrix) IntegratedConstitutiveOverThickness(IList<GaussLegendrePoint3D> gaussPoints,
                int j)
        {
            var MembraneConstitutiveMatrix = new double[3, 3];
            var BendingConstitutiveMatrix = new double[3, 3];
            var CouplingConstitutiveMatrix = new double[3, 3];

            foreach (var keyValuePair in materialsAtThicknessGP[gaussPoints[j]])
            {
                var thicknessPoint = keyValuePair.Key;
                var material = keyValuePair.Value;
                var constitutiveMatrixM = material.ConstitutiveMatrix;
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

        public (double[] MembraneForces, Double[] BendingMoments) IntegratedStressesOverThickness(
            IList<GaussLegendrePoint3D> gaussPoints, int j)
        {
            var MembraneForces = new double[3];
            var BendingMoments = new double[3];

            foreach (var keyValuePair in materialsAtThicknessGP[gaussPoints[j]])
            {
                var thicknessPoint = keyValuePair.Key;
                var material = keyValuePair.Value;
                var w = thicknessPoint.WeightFactor;
                var z = thicknessPoint.Zeta;
                for (int i = 0; i < 3; i++)
                {
                    MembraneForces[i] += material.Stresses[i] * w * Thickness / 2;
                    BendingMoments[i] += material.Stresses[i] * w * z * z * Thickness / 2;
                }
            }

            return (MembraneForces, BendingMoments);
        }


        private double[,] CalculateBendingDeformationMatrix(double[] surfaceBasisVector3,
            ShapeTSplines2DFromBezierExtraction tsplines, int j,
            double[] surfaceBasisVector2, double[] surfaceBasisVectorDerivative1, double[] surfaceBasisVector1,
            double J1,
            double[] surfaceBasisVectorDerivative2, double[] surfaceBasisVectorDerivative12,
            TSplineKirchhoffLoveShellElementMaterial element)
        {
            var Bbending = new double[3, element.ControlPoints.Count * 3];
            var s1 = Vector.CreateFromArray(surfaceBasisVector1);
            var s2 = Vector.CreateFromArray(surfaceBasisVector2);
            var s3 = Vector.CreateFromArray(surfaceBasisVector3);
            var s11 = Vector.CreateFromArray(surfaceBasisVectorDerivative1);
            var s22 = Vector.CreateFromArray(surfaceBasisVectorDerivative2);
            var s12 = Vector.CreateFromArray(surfaceBasisVectorDerivative12);
            for (int column = 0; column < element.ControlPoints.Count * 3; column += 3)
            {
                #region BI1

                var BI1 = s3.CrossProduct(s3);
                BI1.ScaleIntoThis(tsplines.TSplineDerivativeValuesHeta[column / 3, j]);
                var auxVector = s2.CrossProduct(s3);
                auxVector.ScaleIntoThis(tsplines.TSplineDerivativeValuesKsi[column / 3, j]);
                BI1.AddIntoThis(auxVector);
                BI1.ScaleIntoThis(s3.DotProduct(s11));
                auxVector = s1.CrossProduct(s11);
                auxVector.ScaleIntoThis(tsplines.TSplineDerivativeValuesHeta[column / 3, j]);
                BI1.AddIntoThis(auxVector);
                BI1.ScaleIntoThis(1 / J1);
                auxVector[0] = surfaceBasisVector3[0];
                auxVector[1] = surfaceBasisVector3[1];
                auxVector[2] = surfaceBasisVector3[2];
                auxVector.ScaleIntoThis(-tsplines.TSplineSecondDerivativesValueKsi[column / 3, j]);
                BI1.AddIntoThis(auxVector);

                #endregion

                #region BI2

                IVector BI2 = s3.CrossProduct(s3);
                BI2.ScaleIntoThis(tsplines.TSplineDerivativeValuesHeta[column / 3, j]);
                auxVector = s2.CrossProduct(s3);
                auxVector.ScaleIntoThis(tsplines.TSplineDerivativeValuesKsi[column / 3, j]);
                BI2.AddIntoThis(auxVector);
                BI2.ScaleIntoThis(s3.DotProduct(s22));
                auxVector = s1.CrossProduct(s22);
                auxVector.ScaleIntoThis(tsplines.TSplineDerivativeValuesHeta[column / 3, j]);
                BI2.AddIntoThis(auxVector);
                auxVector = s22.CrossProduct(s2);
                auxVector.ScaleIntoThis(tsplines.TSplineDerivativeValuesKsi[column / 3, j]);
                BI2.AddIntoThis(auxVector);
                BI2.ScaleIntoThis(1 / J1);
                auxVector[0] = surfaceBasisVector3[0];
                auxVector[1] = surfaceBasisVector3[1];
                auxVector[2] = surfaceBasisVector3[2];
                auxVector.ScaleIntoThis(-tsplines.TSplineSecondDerivativesValueHeta[column / 3, j]);
                BI2.AddIntoThis(auxVector);

                #endregion

                #region BI3

                Vector BI3 = s3.CrossProduct(s3);
                BI3.ScaleIntoThis(tsplines.TSplineDerivativeValuesHeta[column / 3, j]);
                auxVector = s2.CrossProduct(s3);
                auxVector.ScaleIntoThis(tsplines.TSplineDerivativeValuesKsi[column / 3, j]);
                BI3.AddIntoThis(auxVector);
                BI3.ScaleIntoThis(s3.DotProduct(s12));
                auxVector = s1.CrossProduct(s12);
                auxVector.ScaleIntoThis(tsplines.TSplineDerivativeValuesHeta[column / 3, j]);
                BI3.AddIntoThis(auxVector);
                auxVector = s22.CrossProduct(s2);
                auxVector.ScaleIntoThis(tsplines.TSplineDerivativeValuesKsi[column / 3, j]);
                BI3.AddIntoThis(auxVector);
                BI3.ScaleIntoThis(1 / J1);
                auxVector[0] = surfaceBasisVector3[0];
                auxVector[1] = surfaceBasisVector3[1];
                auxVector[2] = surfaceBasisVector3[2];
                auxVector.ScaleIntoThis(-tsplines.TSplineSecondDerivativesValueKsiHeta[column / 3, j]);
                BI3.AddIntoThis(auxVector);

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

        private double[,] CalculateMembraneDeformationMatrix(ShapeTSplines2DFromBezierExtraction tsplines, int j,
            double[] surfaceBasisVector1,
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

        private static double[,] CalculateHessian(TSplineKirchhoffLoveShellElementMaterial shellElement,
            ShapeTSplines2DFromBezierExtraction tsplines, int j)
        {
            double[,] hessianMatrix = new double[3, 3];
            for (int k = 0; k < shellElement.ControlPoints.Count; k++)
            {
                hessianMatrix[0, 0] +=
                    tsplines.TSplineSecondDerivativesValueKsi[k, j] * shellElement.ControlPoints[k].X;
                hessianMatrix[0, 1] +=
                    tsplines.TSplineSecondDerivativesValueKsi[k, j] * shellElement.ControlPoints[k].Y;
                hessianMatrix[0, 2] +=
                    tsplines.TSplineSecondDerivativesValueKsi[k, j] * shellElement.ControlPoints[k].Z;
                hessianMatrix[1, 0] +=
                    tsplines.TSplineSecondDerivativesValueHeta[k, j] * shellElement.ControlPoints[k].X;
                hessianMatrix[1, 1] +=
                    tsplines.TSplineSecondDerivativesValueHeta[k, j] * shellElement.ControlPoints[k].Y;
                hessianMatrix[1, 2] +=
                    tsplines.TSplineSecondDerivativesValueHeta[k, j] * shellElement.ControlPoints[k].Z;
                hessianMatrix[2, 0] +=
                    tsplines.TSplineSecondDerivativesValueKsiHeta[k, j] * shellElement.ControlPoints[k].X;
                hessianMatrix[2, 1] +=
                    tsplines.TSplineSecondDerivativesValueKsiHeta[k, j] * shellElement.ControlPoints[k].Y;
                hessianMatrix[2, 2] +=
                    tsplines.TSplineSecondDerivativesValueKsiHeta[k, j] * shellElement.ControlPoints[k].Z;
            }

            return hessianMatrix;
        }

        private static double[,] CalculateJacobian(TSplineKirchhoffLoveShellElementMaterial shellElement,
            ShapeTSplines2DFromBezierExtraction tsplines, int j)
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

        Dictionary<GaussLegendrePoint3D, List<GaussLegendrePoint3D>> thicknessIntegrationPoints =
            new Dictionary<GaussLegendrePoint3D, List<GaussLegendrePoint3D>>();

        private IList<GaussLegendrePoint3D> CreateElementGaussPoints(TSplineKirchhoffLoveShellElementMaterial element)
        {
            GaussQuadrature gauss = new GaussQuadrature();
            var medianSurfaceGP = gauss.CalculateElementGaussPoints(element.DegreeKsi, element.DegreeHeta,
                new List<Knot>
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
                            new Knot() {ID = 0, Ksi = -Thickness / 2, Heta = point.Heta},
                            new Knot() {ID = 1, Ksi = Thickness / 2, Heta = point.Heta},
                        }).Select(gp => new GaussLegendrePoint3D(point.Ksi, point.Heta, gp.Ksi, null, gp.WeightFactor))
                    .ToList());
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
            var tsplineElement = (TSplineKirchhoffLoveShellElementMaterial) element;
            var knotParametricCoordinatesKsi = Vector.CreateFromArray(new double[] {-1, 1});
            var knotParametricCoordinatesHeta = Vector.CreateFromArray(new double[] {-1, 1});

            ShapeTSplines2DFromBezierExtraction tsplines = new ShapeTSplines2DFromBezierExtraction(tsplineElement,
                tsplineElement.ControlPoints, knotParametricCoordinatesKsi, knotParametricCoordinatesHeta);

            var knotDisplacements = new double[4, 3];
            var paraviewKnotRenumbering = new int[] {0, 3, 1, 2};
            for (int j = 0; j < knotDisplacements.GetLength(0); j++)
            {
                for (int i = 0; i < element.ControlPoints.Count; i++)
                {
                    knotDisplacements[paraviewKnotRenumbering[j], 0] +=
                        tsplines.TSplineValues[i, j] * localDisplacements[i, 0];
                    knotDisplacements[paraviewKnotRenumbering[j], 1] +=
                        tsplines.TSplineValues[i, j] * localDisplacements[i, 1];
                    knotDisplacements[paraviewKnotRenumbering[j], 2] +=
                        tsplines.TSplineValues[i, j] * localDisplacements[i, 2];
                }
            }

            return knotDisplacements;
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

            var knotParametricCoordinatesKsi = Vector.CreateFromArray(new double[] {-1, 1});
            var knotParametricCoordinatesHeta = Vector.CreateFromArray(new double[] {-1, 1});

            ShapeTSplines2DFromBezierExtraction tsplines = new ShapeTSplines2DFromBezierExtraction(element,
                element.ControlPoints, knotParametricCoordinatesKsi, knotParametricCoordinatesHeta);

            var knotDisplacements = new double[4, 3];
            var paraviewKnotRenumbering = new int[] {0, 3, 1, 2};
            for (int j = 0; j < localCoordinates.GetLength(0); j++)
            {
                for (int i = 0; i < element.ControlPoints.Count; i++)
                {
                    knotDisplacements[paraviewKnotRenumbering[j], 0] +=
                        tsplines.TSplineValues[i, j] * element.ControlPoints[i].X;
                    knotDisplacements[paraviewKnotRenumbering[j], 1] +=
                        tsplines.TSplineValues[i, j] * element.ControlPoints[i].Y;
                    knotDisplacements[paraviewKnotRenumbering[j], 2] +=
                        tsplines.TSplineValues[i, j] * element.ControlPoints[i].Z;
                }
            }

            return knotDisplacements;
        }
    }
}