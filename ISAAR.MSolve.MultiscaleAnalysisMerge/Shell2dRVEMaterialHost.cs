using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Materials;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.MultiscaleAnalysis.Interfaces;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.FEM.Providers;
using ISAAR.MSolve.MultiscaleAnalysis.SupportiveClasses;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Discretization.Providers;
using ISAAR.MSolve.Materials.Interfaces;
using System;
using System.Collections.Generic;
using ISAAR.MSolve.MultiscaleAnalysis;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra;

namespace ISAAR.MSolve.MultiscaleAnalysisMerge
{
    public class Shell2dRVEMaterialHost : IShellMaterial
    {
        public double[] NormalVectorV3 { get; set; }
        public double[] TangentVectorV1 { get; set; }
        public double[] TangentVectorV2 { get; set; }
        private IMatrix2D consMat;

        private int rveDatabaseSize;
        private int gpsPerRve;
        private int gpCounter;
        private IdegenerateRVEbuilder rveBuilderToClone;

        IShellMaterial coreMaterial;

        private Dictionary<int, IShellMaterial> coreMaterials;

        public Shell2dRVEMaterialHost(int rveDatabaseSize, int gpsPerRve, int gpCounter, IdegenerateRVEbuilder rVEbuilder)
        {
            this.rveDatabaseSize = rveDatabaseSize;
            this.gpsPerRve = gpsPerRve;
            this.gpCounter = gpCounter;
            coreMaterials = new Dictionary<int, IShellMaterial>();
            this.rveBuilderToClone = rVEbuilder;
        }

        private Shell2dRVEMaterialHost(IShellMaterial coreMaterial)
        {
            this.coreMaterial = coreMaterial;
        }

        public IShellMaterial Clone()
        {
            gpCounter += 1;

            if ((gpCounter/ (rveDatabaseSize*gpsPerRve)==1)&&(gpCounter% (rveDatabaseSize * gpsPerRve)==1))
            {
                gpCounter += -(rveDatabaseSize * gpsPerRve);
            }

            int rve_id = ((gpCounter-1) / gpsPerRve) + 1;

            if(!coreMaterials.ContainsKey(rve_id))
            {
                coreMaterials[rve_id] = new
                Microstructure3DevelopMultipleSubdomainsUseBaseSmallStrainsShelltransformationSimu((IdegenerateRVEbuilder)rveBuilderToClone.Clone(rve_id), true);
                coreMaterials[rve_id].TangentVectorV1 = new double[3] { 1, 0, 0 };
                coreMaterials[rve_id].TangentVectorV2 = new double[3] { 0, 1, 0 };
                coreMaterials[rve_id].NormalVectorV3 = new double[3] { 0, 0, 1 };

            }

            var sharedRVEmat = new Shell2dRVEMaterialHost(coreMaterials[rve_id]);
            return sharedRVEmat;
        }

        object ICloneable.Clone()
        {
            return this.Clone();
        }

        public IMatrix2D ConstitutiveMatrix
        {
            get
            {
                if (consMat == null) BuildConstitutive();
                return consMat;
            }
        }

        private void BuildConstitutive()
        {
            var transformationMatrix=this.CalculateTransformationMatrix(new Vector(TangentVectorV1), new Vector(TangentVectorV2));
            var coreConstCopy = new double[3, 3];
            for (int i1 = 0; i1 < 3; i1++) { for (int i2 = 0; i2 < 3; i2++) { coreConstCopy[i1, i2] = coreMaterial.ConstitutiveMatrix[i1, i2]; } }
            consMat = (transformationMatrix.Transpose() * (new Matrix2D(coreConstCopy)) * transformationMatrix);
        }

        private Matrix2D CalculateTransformationMatrix(Vector surfaceBasisVector1, Vector surfaceBasisVector2)
        {
            var auxMatrix1 = new Matrix2D(2, 2);  //auxMatrix: covariant metric coefficients gab
            auxMatrix1[0, 0] = surfaceBasisVector1.DotProduct(surfaceBasisVector1);
            auxMatrix1[0, 1] = surfaceBasisVector1.DotProduct(surfaceBasisVector2);
            auxMatrix1[1, 0] = surfaceBasisVector2.DotProduct(surfaceBasisVector1);
            auxMatrix1[1, 1] = surfaceBasisVector2.DotProduct(surfaceBasisVector2);
            (Matrix2D inverse, double det) = auxMatrix1.Invert2x2AndDeterminant(1e-20); //inverse: contravariant metric coefficients g_ab (ekthetis ta a,b)

            //Contravariant base vectors
            double[][] G_i = new double[2][];
            for (int i1 = 0; i1 < 2; i1++)
            {
                G_i[i1] = new double[3];
                for (int i2 = 0; i2 < 3; i2++)
                {
                    G_i[i1][i2] = inverse[i1, 0] * surfaceBasisVector1[i2] + inverse[i1, 1] * surfaceBasisVector2[i2];
                }
            }

            //Normalised covariant base vectors
            double[][] Ei = new double[2][];// to trito den xreiazetai

            Ei[0] = new double[3]; Array.Copy(surfaceBasisVector1.Data, Ei[0], 3);
            double G1_norm = surfaceBasisVector1.Norm;
            for (int i1 = 0; i1 < 3; i1++) { Ei[0][i1] = Ei[0][i1] / G1_norm; }

            double G2_dot_E1 = 0;
            for (int i1 = 0; i1 < 3; i1++) { G2_dot_E1 += surfaceBasisVector2[i1] * Ei[0][i1]; }

            double[] projection = new double[3];
            for (int i1 = 0; i1 < 3; i1++) { projection[i1] = G2_dot_E1 * Ei[0][i1]; }

            Ei[1] = new double[3];
            for (int i1 = 0; i1 < 3; i1++) { Ei[1][i1] = surfaceBasisVector2[i1] - projection[i1]; }
            double norm1 = new Vector(Ei[1]).Norm;
            for (int i1 = 0; i1 < 3; i1++) { Ei[1][i1] = Ei[1][i1] / norm1; }

            double[,] EiDOTG_j = new double[2, 2];

            for (int i1 = 0; i1 < 2; i1++)
            {
                for (int i2 = 0; i2 < 2; i2++)
                {
                    EiDOTG_j[i1, i2] = new Vector(Ei[i1]).DotProduct(new Vector(G_i[i2]));
                }
            }

            var transformationMatrix = new Matrix2D(new double[3, 3] { {EiDOTG_j[0,0]*EiDOTG_j[0,0],EiDOTG_j[0,1]*EiDOTG_j[0,1],EiDOTG_j[0,0]*EiDOTG_j[0,1]  },
                 {EiDOTG_j[1,0]*EiDOTG_j[1,0],EiDOTG_j[1,1]*EiDOTG_j[1,1],EiDOTG_j[1,0]*EiDOTG_j[1,1]  },
                {2*EiDOTG_j[1,0]*EiDOTG_j[0,0],2*EiDOTG_j[1,1]*EiDOTG_j[0,1],EiDOTG_j[1,0]*EiDOTG_j[0,1]+EiDOTG_j[1,1]*EiDOTG_j[0,0]   } });

            return transformationMatrix;
        }

        public double[] Stresses => throw new NotSupportedException();

        public int ID => throw new NotImplementedException();

        public bool Modified => throw new NotImplementedException();

        public double[] Coordinates { get => throw new NotImplementedException(); set => throw new NotImplementedException(); }

        public double YoungModulus => throw new NotSupportedException();

        public double PoissonRatio => throw new NotSupportedException();

        
        public void UpdateMaterial(double[] GLvec)
        {
            throw new NotSupportedException();
        }

        public void ResetModified()
        {
            throw new NotSupportedException();
        }

        public void SaveState()
        {
            throw new NotSupportedException();
        }

        public void ClearState()
        {
            throw new NotSupportedException();
        }

        public void ClearStresses()
        {
            throw new NotSupportedException();
        }
    }
}
