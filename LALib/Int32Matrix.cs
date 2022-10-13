//
// MatrixInt32.cs
//
// Author:
//       Tiago C. Teixeira <>
//
// Copyright (c) 2022 
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.


using System;

namespace LALib
{

    /// <summary>
    /// Matrix of int32 elements.
    /// </summary>
    public class Int32Matrix : Matrix<Int32>
    {

        public Int32Matrix()
            : base( new IntMathOperations() )
        {
        }


        /// <summary>
        /// Clones an int matrix to a double matrix.
        /// </summary>
        /// <returns>Returns a double matrix.</returns>
        /// <param name="A">Source matrix.</param>
        private static double[][] CloneToDoubleMatrix(int[][] A)
        {
            int numR = GetNumRows(A);
            int numC = GetNumCols(A);
            double[][] FResult = Matrix<double>.CreateJaggedArray(numR, numC);

            for (int i = 0; i < numR; i++)
            {
                for (int j = 0; j < numC; j++)
                {
                    FResult[i][j] = A[i][j];
                }
            }

            return FResult;
        }


        /// <summary>
        /// Computes the determinant for double matrices using Bareiss alg.
        /// </summary>
        /// <returns>Returns the determinant of matrix.</returns>
        /// <param name="m1">The original matrix.</param>
        /// <param name="n">Matrix<T> dimension.</param>
        /// <remarks>
        ///     Note: The determinant applies only to a square matrix.
        /// </remarks>
        public static double DET_BareissAlg(int[][] m1, int n)
        {
            double FResult;
            double[][] mat = CloneToDoubleMatrix(m1);

            mat[0][0] = 1.0d;

            for (int k = 1; k < n; k++)
            {
                for (int i = k + 1; i < n; i++)
                {
                    for (int j = k + 1; j < n; j++)
                    {
                        mat[i][j] = ( mat[i][j] * mat[k][k] - mat[i][k] * mat[k][j] ) / mat[k - 1][k - 1];
                    }
                }
            }

            FResult = mat[n - 1][n - 1];
            return FResult;
        }


        /// <summary>
        /// Generates a matrix filled with random values.
        /// </summary>
        /// <returns>Returns a matrix filled with random values.</returns>
        /// <param name="rows">Number of matrix rows.</param>
        /// <param name="cols">Number of matrix cols.</param>
        /// <param name="minVal">Minimum random value.</param>
        /// <param name="maxVal">Maximum random value.</param>
        /// <param name="seed">Seed value for random generator.</param>
        public static int[][] Random(int rows, int cols,
                        int minVal, int maxVal, int seed)
        {
            // return a matrix with random values
            Random ran = new Random(seed);
            int[][] result = Matrix<int>.CreateJaggedArray(rows, cols);

            for (int i = 0; i < rows; ++i)
                for (int j = 0; j < cols; ++j)
                    result[i][j] = ran.Next(minVal, maxVal);

            return result;
        }


    }

}
