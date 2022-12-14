/*****************************************************************
 * Copyright 2022, Tiago C. Teixeira
 *
 * MIT license (Massachusetts Institute of Technology)
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy 
 * of this software and associated documentation files (the "Software"), to deal 
 * in the Software without restriction, including without limitation the rights 
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies 
 * of the Software, and to permit persons to whom the Software is furnished to 
 * do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all 
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING 
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS 
 * IN THE SOFTWARE.
 *  
 * 
 * This code was adapted from James D. McCaffrey blog
 * https://jamesmccaffrey.wordpress.com/2015/03/06/inverting-a-matrix-using-c/
 * 
 * Changes were made to suport several numeric types 
 * (like int, complex numbers, double, etc) instead of 
 * only double type matrices
 * 
 * 
 *****************************************************************/

using System;
using System.Threading.Tasks;

namespace LALib
{

    /// <summary>
    /// Helps calculating the inverse of a matrix.
    /// </summary>
    public class InverseHelper<T> //where T : IComparable
    {

        private IBasicMathOperations<T> Op { get; set; }

        public InverseHelper(IBasicMathOperations<T> _op)
        {
            this.Op = _op;
        }


        /// <summary>
        /// Helpers to solve linear equations of LU decomposed matrix.
        /// </summary>
        /// <returns>Returns the result for inverse calculation.</returns>
        /// <param name="luMatrix">LU decomposed matrix.</param>
        /// <param name="b">The second member terms.</param>
        /// <typeparam name="T">The 1st type parameter.</typeparam>
        /// <remarks>
        /// IMPORTANT NOTE: 
        ///     Before calling this helper, permute b using the perm array
        ///     from MatrixDecompose that generated luMatrix
        /// </remarks>
        private T[] HelperSolve(T[][] luMatrix, T[] b)
        {
            // before calling this helper, permute b using the perm array
            // from MatrixDecompose that generated luMatrix
            int n = luMatrix.Length;
            T[] x = new T[n];
            b.CopyTo(x, 0);
            T sum = this.Op.Zero;// 0.0;
            //dynamic d1 = 0.0;
            //dynamic d2 = 0.0;

            for (int i = 1; i < n; ++i)
            {
                sum = x[i];
                for (int j = 0; j < i; ++j)
                {
                    sum = this.Op.Sub(sum, this.Op.Mult(luMatrix[i][j], x[j]));

                    //d1 = luMatrix[i][j]; d2 = x[j];
                    //sum -= d1 * d2;
                    ////sum -= luMatrix[i][j] * x[j];
                }

                x[i] = sum;
            }

            //d1 = luMatrix[n - 1][n - 1];
            //d2 = x[n - 1];
            //d2 /= d1;
            //x[n - 1] = d2;

            x[n - 1] = this.Op.Div( x[n - 1], luMatrix[n - 1][n - 1] );


            //x[n - 1] /= luMatrix[n - 1][n - 1];
            for (int i = n - 2; i >= 0; --i)
            {

                sum = x[i];
                for (int j = i + 1; j < n; ++j)
                {
                    //d1 = luMatrix[i][j]; //d2 = x[j];
                    //sum -= d1 * x[j];

                    sum = this.Op.Sub(sum, this.Op.Mult(luMatrix[i][j], x[j]));

                    //sum -= luMatrix[i][j] * x[j];
                }

                x[i] = this.Op.Div( sum, luMatrix[i][i] );
            }

            return x;
        }


        /// <summary>
        /// Computes the inverse a specified original matrix usind his LUP decomposed matrix.
        /// </summary>
        /// <returns>Returns the inverse matrix if succeeded.</returns>
        /// <param name="LUP">The LU decomposed matrix of the original matrix.</param>
        /// <param name="perm">The permutation vector from decomposition process.</param>
        /// <typeparam name="T">The 1st type parameter.</typeparam>
        /// <remarks>
        /// 
        ///     Note: Not all matrices are invertible.
        /// 
        ///     Non-square matrices (m-by-n matrices for which m ≠ n) do not have 
        ///     an inverse. However, in some cases such a matrix may have a left 
        ///     inverse or right inverse.
        /// 
        ///     In linear algebra, an n-by-n square matrix A is called invertible 
        ///     (also nonsingular or nondegenerate), if there exists an n-by-n 
        ///     square matrix B such that AxB = BxA = In  
        ///     where In denotes the n-by-n identity matrix and the multiplication 
        ///     used is ordinary matrix multiplication.
        ///
        ///     A square matrix that is not invertible is called singular or degenerate. 
        ///     A square matrix is singular if and only if its determinant is zero.
        /// 
        /// </remarks>
        public T[][] Inverse(T[][] LUP, int[] perm)
        {
            T[][] FResult = null;

            if (LUP == null)
                throw new Exception("Unable to compute inverse.");

            int n = LUP.Length;
            FResult = Matrix<T>.CreateJaggedArray(n, n); // initialize result arrays

            T[] b = new T[n];
            T val;
            for (int i = 0; i < n; ++i)
            {
                for (int j = 0; j < n; ++j)
                {
                    if (i == perm[j])
                        b[j] = (val = this.Op.One);
                    else
                        b[j] = (val = this.Op.Zero);
                }

                T[] x = this.HelperSolve(LUP, b);

                for (int j = 0; j < n; ++j)
                    FResult[j][i] = x[j];
            }

            return FResult;
        }


        /// <summary>
        /// Computes the inverse a specified original matrix usind his LUP decomposed matrix
        /// using parallelization.
        /// </summary>
        /// <returns>Returns the inverse matrix if succeeded.</returns>
        /// <param name="LUP">The LU decomposed matrix of the original matrix.</param>
        /// <param name="perm">The permutation vector from decomposition process.</param>
        /// <typeparam name="T">The 1st type parameter.</typeparam>
        /// <remarks>
        /// 
        ///     Note: Not all matrices are invertible.
        /// 
        ///     Non-square matrices (m-by-n matrices for which m ≠ n) do not have 
        ///     an inverse. However, in some cases such a matrix may have a left 
        ///     inverse or right inverse.
        /// 
        ///     In linear algebra, an n-by-n square matrix A is called invertible 
        ///     (also nonsingular or nondegenerate), if there exists an n-by-n 
        ///     square matrix B such that AxB = BxA = In  
        ///     where In denotes the n-by-n identity matrix and the multiplication 
        ///     used is ordinary matrix multiplication.
        ///
        ///     A square matrix that is not invertible is called singular or degenerate. 
        ///     A square matrix is singular if and only if its determinant is zero.
        /// 
        /// </remarks>
        public T[][] InversePar(T[][] LUP, int[] perm)
        {
            T[][] FResult = null;

            if (LUP == null)
                throw new Exception("Unable to compute inverse.");

            int n = LUP.Length;
            FResult = Matrix<T>.CreateJaggedArrayPar(n, n); // initialize result arrays

            T[] b = new T[n];
            //T val;

            for (int i = 0; i < n; ++i) 
            {
                Parallel.For(0, n, delegate(int j)
                    {
                        if (i == perm[j])
                            b[j] = this.Op.One;     //(val = this.Op.One);
                        else
                            b[j] = this.Op.Zero;    //(val = this.Op.Zero);
                    }
                );

                T[] x = this.HelperSolve(LUP, b);

                Parallel.For(0, n, delegate (int j)
                    {
                        FResult[j][i] = x[j];
                    }
                );
            }

            return FResult;
        }


    }

}
