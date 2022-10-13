//
// MatrixDouble.cs
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
    /// Matrix of double elements.
    /// </summary>
    public class DoubleMatrix : Matrix<double>
    {

        public DoubleMatrix() 
            : base( new DoubleMathOperations() )
        {
        }


        /// <summary>
        /// Computes the determinant for double matrices using Bareiss alg.
        /// </summary>
        /// <returns>Returns the determinant of matrix.</returns>
        /// <param name="m1">The original matrix.</param>
        /// <param name="n">Matrix<T> dimension.</param>
        /// <param name="overwrite">If <c>true</c> original matrix will be changed.</param>
        /// <remarks>
        ///     Note: The determinant applies only to a square matrix.
        /// </remarks>
        public double DET_BareissAlg( double[][] m1, int n, bool overwrite )
        {
            double FResult = base.DET_BareissAlg_Base_T(m1, n, overwrite);
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
        public static double[][] Random(int rows, int cols,
                        double minVal, double maxVal, int seed)
        {
            // return a matrix with random values
            Random ran = new Random(seed);
            double[][] result = Matrix<double>.CreateJaggedArray(rows, cols);

            for (int i = 0; i < rows; ++i)
                for (int j = 0; j < cols; ++j)
                    result[i][j] = (maxVal - minVal) *
                        ran.NextDouble() + minVal;

            return result;
        }


    }

}
