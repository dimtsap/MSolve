﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.LinearAlgebra.ArrayManipulations
{
    // TODO: Have a separate conversions class for testing and use MKL (BLAS) routines.
    public static class Conversions
    {
        public static double[,] Array2DLowerToSymmetric(double[,] array2D)
        {
            if (array2D.GetLength(0) != array2D.GetLength(1))
            {
                throw new ArgumentException("The provided matrix is not square");
            }
            int n = array2D.GetLength(0);
            double[,] symm = new double[n, n];
            Array.Copy(array2D, symm, n * n);
            for (int i = 0; i < n; ++i)
            {
                for (int j = 0; j < i; ++j)
                {
                    symm[j, i] = symm[i, j];
                }
            }
            return symm;
        }

        public static double[,] Array2DUpperToSymmetric(double[,] array2D)
        {
            if (array2D.GetLength(0) != array2D.GetLength(1))
            {
                throw new ArgumentException("The provided matrix is not square");
            }
            int n = array2D.GetLength(0);
            double[,] symm = new double[n, n];
            Array.Copy(array2D, symm, n * n);
            for (int i = 0; i < n; ++i)
            {
                for (int j = 0; j < i; ++j)
                {
                    symm[i, j] = symm[j, i];
                }
            }
            return symm;
        }

        public static double[] Array2DToFullColMajor(double[,] array2D)
        {
            int numRows = array2D.GetLength(0);
            int numColumns = array2D.GetLength(1);
            double[] array1D = new double[numRows * numColumns];
            int idxCounter = -1;
            for (int j = 0; j < numColumns; ++j) //The order of loops is important
            {
                for (int i = 0; i < numRows; ++i)
                {
                    array1D[++idxCounter] = array2D[i, j];
                }
            }
            return array1D;
        }

        public static double[] Array2DToFullRowMajor(double[,] array2D)
        {
            int numRows = array2D.GetLength(0);
            int numColumns = array2D.GetLength(1);
            double[] array1D = new double[numRows * numColumns];
            for (int i = 0; i < numRows; ++i)
            {
                for (int j = 0; j < numColumns; ++j)
                {
                    array1D[i * numColumns + j] = array2D[i, j];
                }
            }
            return array1D;
        }

        public static double[] Array2DToPackedLowerColMajor(double[,] array2D)
        {
            if (array2D.GetLength(0) != array2D.GetLength(1))
            {
                throw new ArgumentException("The provided matrix is not square");
            }
            int n = array2D.GetLength(0);
            double[] array1D = new double[(n * (n + 1)) / 2];
            int counter = 0; // Simplifies indexing but the outer and inner loops cannot be interchanged
            for (int j = 0; j < n; ++j)
            {
                for (int i = j; i < n; ++i)
                {
                    array1D[counter] = array2D[i, j];
                    ++counter; // Clearer than post-incrementing during indexing.
                }
            }
            return array1D;
        }

        public static double[] Array2DToPackedLowerRowMajor(double[,] array2D)
        {
            if (array2D.GetLength(0) != array2D.GetLength(1))
            {
                throw new ArgumentException("The provided matrix is not square");
            }
            int n = array2D.GetLength(0);
            double[] array1D = new double[(n * (n + 1)) / 2];
            int counter = 0; // Simplifies indexing but the outer and inner loops cannot be interchanged
            for (int i = 0; i < n; ++i)
            {
                for (int j = 0; j <= i; ++j)
                {
                    array1D[counter] = array2D[i, j];
                    ++counter; // Clearer than post-incrementing during indexing.
                }
            }
            return array1D;
        }

        public static double[] Array2DToPackedUpperColMajor(double[,] array2D)
        {
            if (array2D.GetLength(0) != array2D.GetLength(1))
            {
                throw new ArgumentException("The provided matrix is not square");
            }
            int n = array2D.GetLength(0);
            double[] array1D = new double[(n * (n + 1)) / 2];
            int counter = 0; // Simplifies indexing but the outer and inner loops cannot be interchanged
            for (int j = 0; j < n; ++j)
            {
                for (int i = 0; i <= j; ++i)
                {
                    array1D[counter] = array2D[i, j];
                    ++counter; // Clearer than post-incrementing during indexing.
                }
            }
            return array1D;
        }

        public static double[] Array2DToPackedUpperRowMajor(double[,] array2D)
        {
            if (array2D.GetLength(0) != array2D.GetLength(1))
            {
                throw new ArgumentException("The provided matrix is not square");
            }
            int n = array2D.GetLength(0);
            double[] array1D = new double[(n * (n + 1)) / 2];
            int counter = 0; // Simplifies indexing but the outer and inner loops cannot be interchanged
            for (int i = 0; i < n; ++i)
            {
                for (int j = i; j < n; ++j)
                {
                    array1D[counter] = array2D[i, j];
                    ++counter; // Clearer than post-incrementing during indexing.
                }
            }
            return array1D;
        }

        public static double[] ColumnMajorToRowMajor(double[] colMajor, int numRows, int numCols)
        {
            double[] rowMajor = new double[numRows * numCols];
            int idxCounter = -1;
            for (int i = 0; i < numRows; ++i) // The order of loops is important
            {
                for (int j = 0; j < numCols; ++j)
                {
                    rowMajor[++idxCounter] = colMajor[j * numRows + i];
                }
            }
            return rowMajor;
        }

        // TODO: This is no longer a conversion per se. I need to organize these methods
        public static void CopyUpperToLowerColMajor(double[] squareFullMatrix, int order)
        {
            for (int j = 0; j < order; ++j)
            {
                for (int i = 0; i < j; ++i)
                {
                    squareFullMatrix[i * order + j] = squareFullMatrix[j * order + i];
                }
            }
        }

        public static double[,] FullColMajorToArray2D(double[] array1D, int numRows, int numColumns)
        {
            double[,] array2D = new double[numRows, numColumns];
            for (int j = 0; j < numColumns; ++j)
            {
                for (int i = 0; i < numRows; ++i)
                {
                    array2D[i, j] = array1D[j * numRows + i];
                }
            }
            return array2D;
        }

        /// <summary>
        /// Extract the lower part only. Both I/O are column major with full storage.
        /// </summary>
        /// <param name="full"></param>
        /// <param name="unitDiagonal"></param>
        /// <returns></returns>
        public static double[] FullColMajorToFullLowerRowMajor(double[] full, bool unitDiagonal)
        {
            int n = FullLengthToOrder(full.Length);
            double[] lower = new double[n * n];
            for (int j = 0; j < n; ++j)
            {
                for (int i = j + 1; i < n; ++i)
                {
                    lower[i * n + j] = full[j * n + i]; // lower is row major, full is col major
                }
            }
            if (unitDiagonal)
            {
                for (int d = 0; d < n; ++d) lower[d * n + d] = 1.0;
            }
            else
            {
                for (int d = 0; d < n; ++d) lower[d * n + d] = full[d * n + d];
            }
            return lower;
        }

        /// <summary>
        /// Extract the lower part only. Both I/O are column major with full storage.
        /// </summary>
        /// <param name="full"></param>
        /// <param name="unitDiagonal"></param>
        /// <returns></returns>
        public static double[] FullColMajorToFullLowerColMajor(double[] full, bool unitDiagonal)
        {
            int n = FullLengthToOrder(full.Length);
            double[] lower = new double[n * n];
            for (int j = 0; j < n; ++j)
            {
                for (int i = j + 1; i < n; ++i)
                {
                    lower[j * n + i] = full[j * n + i];
                }
            }
            if (unitDiagonal)
            {
                for (int d = 0; d < n; ++d) lower[d * n + d] = 1.0;
            }
            else
            {
                for (int d = 0; d < n; ++d) lower[d * n + d] = full[d * n + d];
            }
            return lower;
        }

        /// <summary>
        /// Extract the lower trapezoid (subdiagonal) part only. Both I/O are column major with full storage. Appropriate for 
        /// rectangular matrices.
        /// </summary>
        /// <param name="full"></param>
        /// <returns></returns>
        public static double[] FullColMajorToFullLowerColMajorRect(int numRows, int numCols, double[] full)
        {
            if (numRows > numCols) throw new NotImplementedException("For now: numRows <= numCols");
            double[] lower = new double[full.Length];
            for (int j = 0; j < numRows; ++j) // numRows < numCols, so this will not scan past the square submatrix that contains the diagonal
            {
                for (int i = j; i < numRows; ++i) // this will scan all rows of the relevant columns
                {
                    lower[j * numRows + i] = full[j * numRows + i];
                }
            }
            return lower;
        }

        /// <summary>
        /// Extract the upper part only. Both I/O are column major with full storage.
        /// </summary>
        /// <param name="full"></param>
        /// <param name="unitDiagonal"></param>
        /// <returns></returns>
        public static double[] FullColMajorToFullUpperColMajor(double[] full, bool unitDiagonal)
        {
            int n = FullLengthToOrder(full.Length);
            double[] upper = new double[n * n];
            for (int j = 0; j < n; ++j)
            {
                for (int i = 0; i < j; ++i)
                {
                    upper[j * n + i] = full[j * n + i];
                }
            }
            if (unitDiagonal)
            {
                for (int d = 0; d < n; ++d) upper[d * n + d] = 1.0;
            }
            else
            {
                for (int d = 0; d < n; ++d) upper[d * n + d] = full[d * n + d];
            }
            return upper;
        }

        /// <summary>
        /// Extract the upper trapezoid (superdiagonal) part only. Both I/O are column major with full storage. Appropriate for 
        /// rectangular matrices.
        /// </summary>
        /// <param name="full"></param>
        /// <returns></returns>
        public static double[] FullColMajorToFullUpperColMajorRect(int numRows, int numCols, double[] full)
        {
            if (numCols > numRows) throw new NotImplementedException("For now: numRows >= numCols");
            double[] upper = new double[full.Length];
            for (int j = 0; j < numCols; ++j) // numCols < numRows, so this will scan all columns
            {
                for (int i = 0; i <= j; ++i) // numCols < numRows, so this will not scan past the square submatrix that contains the diagonal
                {
                    upper[j * numRows + i] = full[j * numRows + i];
                }
            }
            return upper;
        }

        /// <summary>
        /// Extract the upper part only. Both I/O are column major with full storage.
        /// </summary>
        /// <param name="full"></param>
        /// <param name="order"></param>
        /// <returns></returns>
        public static double[] FullColMajorToPackedUpperColMajor(double[] full, int order)
        {
            //int n = FullLengthToOrder(full.Length);
            int n = order;
            double[] upper = new double[(n * (n + 1)) / 2];
            int counter = 0; // Simplifies indexing but the outer and inner loops cannot be interchanged
            for (int j = 0; j < n; ++j)
            {
                for (int i = 0; i <= j; ++i)
                {
                    upper[counter] = full[j * n + i];
                    ++counter; // Clearer than post-incrementing during indexing.
                }
            }
            return upper;
        }

        public static double[,] FullRowMajorToArray2D(double[] array1D, int numRows, int numColumns)
        {
            double[,] array2D = new double[numRows, numColumns];
            for (int i = 0; i < numRows; ++i)
            {
                for (int j = 0; j < numColumns; ++j)
                {
                    array2D[i, j] = array1D[i * numColumns + j];
                }
            }
            return array2D;
        }

        public static double[,] FullLowerColMajorToArray2D(double[] array1D, bool unitDiagonal)
        {
            int n = FullLengthToOrder(array1D.Length);
            double[,] array2D = new double[n, n];
            for (int j = 0; j < n; ++j)
            {
                for (int i = j + 1; i < n; ++i)
                {
                    array2D[i, j] = array1D[j * n + i];
                }
            }
            if (unitDiagonal)
            {
                for (int d = 0; d < n; ++d) array2D[d, d] = 1.0;
            }
            else
            {
                for (int d = 0; d < n; ++d) array2D[d, d] = array1D[d * n + d];
            }
            return array2D;
        }

        public static double[,] FullUpperColMajorToArray2D(double[] array1D, bool unitDiagonal)
        {
            int n = FullLengthToOrder(array1D.Length);
            double[,] array2D = new double[n, n];
            for (int j = 0; j < n; ++j)
            {
                for (int i = 0; i < j; ++i)
                {
                    array2D[i, j] = array1D[j * n + i];
                }
            }
            if (unitDiagonal)
            {
                for (int d = 0; d < n; ++d) array2D[d, d] = 1.0;
            }
            else
            {
                for (int d = 0; d < n; ++d) array2D[d, d] = array1D[d * n + d];
            }
            return array2D;
        }

        public static double[,] PackedLowerColMajorToArray2D(double[] array1D)
        {
            int n = PackedLengthToOrder(array1D.Length);
            double[,] array2D = new double[n, n];
            int counter = 0; // Simplifies indexing but the outer and inner loops cannot be interchanged
            for (int j = 0; j < n; ++j)
            {
                for (int i = j; i < n; ++i)
                {
                    array2D[i, j] = array1D[counter];
                    ++counter; // Clearer than post-incrementing during indexing.
                }
            }
            return array2D;
        }

        public static double[,] PackedLowerRowMajorToArray2D(double[] array1D)
        {
            int n = PackedLengthToOrder(array1D.Length);
            double[,] array2D = new double[n, n];
            int counter = 0; // Simplifies indexing but the outer and inner loops cannot be interchanged
            for (int i = 0; i < n; ++i)
            {
                for (int j = 0; j <= i; ++j)
                {
                    array2D[i, j] = array1D[counter];
                    ++counter; // Clearer than post-incrementing during indexing.
                }
            }
            return array2D;
        }

        public static double[,] PackedUpperColMajorToArray2D(double[] array1D, int order = 0)
        {
            int n = (order == 0) ? PackedLengthToOrder(array1D.Length) : order;
            double[,] array2D = new double[n, n];
            int counter = 0; // Simplifies indexing but the outer and inner loops cannot be interchanged
            for (int j = 0; j < n; ++j)
            {
                for (int i = 0; i <= j; ++i)
                {
                    array2D[i, j] = array1D[counter];
                    ++counter; // Clearer than post-incrementing during indexing.
                }
            }
            return array2D;
        }

        public static double[,] PackedUpperColMajorToArray2DSymm(double[] packed, int order = 0)
        {
            int n = (order == 0) ? PackedLengthToOrder(packed.Length) : order;
            double[, ] array2D = new double[n, n];
            for (int j = 0; j < n; ++j)
            {
                for (int i = 0; i < j; ++i)
                {
                    double val = packed[i + (j * (j + 1)) / 2]; ;
                    array2D[i, j] = val;
                    array2D[j, i] = val;
                }
            }
            for (int i = 0; i < n; ++i) // The diagonal entries need to be proceessed only once.
            {
                array2D[i, i] = packed[i + (i * (i + 1)) / 2];
            }
            return array2D;
        }

        public static double[] PackedUpperColMajorToFullColMajor(double[] packed, int order = 0)
        {
            int n = (order == 0) ? PackedLengthToOrder(packed.Length) : order;
            double[] full = new double[n * n];
            for (int j = 0; j < n; ++j)
            {
                for (int i = 0; i <= j; ++i)
                {
                    double val = packed[i + (j * (j + 1)) / 2]; ;
                    full[j * n + i] = val;
                }
            }
            return full;
        }

        /// <summary>
        /// Converts a symmetric matrix from packed column major format to full column major format.
        /// </summary>
        /// <param name="packed">The matrix entries in packed format. Its length must be n*(n+1)/2, where n is the order of the 
        ///     matrix.</param>
        /// <param name="order">The order of the matrix. It must be positive and match the length of <see cref="packed"/>. If a 
        ///     value is provided, these will not be checked. If no value is provided, the order will be calculated from 
        ///     <see cref="packed"/> instead.</param>
        /// <returns></returns>
        public static double[] PackedUpperColMajorToFullSymmColMajor(double[] packed, int order = 0)
        {
            int n = (order == 0) ? PackedLengthToOrder(packed.Length) : order;
            double[] full = new double[n * n];
            for (int j = 0; j < n; ++j)
            {
                for (int i = 0; i < j; ++i)
                {
                    double val = packed[i + (j * (j + 1)) / 2]; ;
                    full[j * n + i] = val;
                    full[i * n + j] = val;
                }
            }
            for (int i = 0; i < n; ++i) // The diagonal entries need to be proceessed only once.
            {
                full[i * n + i] = packed[i + (i * (i + 1)) / 2];
            }
            return full;
        }

        public static double[,] PackedUpperRowMajorToArray2D(double[] array1D)
        {
            int n = PackedLengthToOrder(array1D.Length);
            double[,] array2D = new double[n, n];
            int counter = 0; // Simplifies indexing but the outer and inner loops cannot be interchanged
            for (int i = 0; i < n; ++i)
            {
                for (int j = i; j < n; ++j)
                {
                    array2D[i, j] = array1D[counter];
                    ++counter; // Clearer than post-incrementing during indexing.
                }
            }
            return array2D;
        }

        public static double[] PackedLowerRowMajorToFullColMajor(double[] packed, int order = 0)
        {
            int n = (order == 0) ? PackedLengthToOrder(packed.Length) : order;
            double[] full = new double[n * n];
            for (int j = 0; j < n; ++j)
            {
                for (int i = 0; i <= j; ++i)
                {
                    double val = packed[j + ((i + 1) * i) / 2]; ;
                    full[j * n + i] = val;
                }
            }
            return full;
        }

        public static int FullLengthToOrder(int length)
        {
            // length = n^2 => n = sqrt(length)
            double n = Math.Sqrt(length);
            int order = (int)Math.Round(n);
            if (order * order != length)
            {
                throw new ArgumentException("The length of the 1D array must be an integer L such that L=n*n");
            }
            return order;
        }

        public static int PackedLengthToOrder(int length)
        {
            // length = n*(n+1)/2 => n = ( -1+sqrt(1+8*length) )/2
            double n = (-1.0 + Math.Sqrt(1 + 8 * length)) / 2;
            int order = (int)Math.Round(n);
            if (order * (order + 1) / 2 != length)
            {
                throw new ArgumentException("The length of the 1D array must be an integer L such that L=n*(n+1)/2");
            }
            return order;
        }
    }
}
