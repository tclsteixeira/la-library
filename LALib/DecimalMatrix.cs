//********************************************************************************
// Copyright 2022, Tiago C. Teixeira
//
// MIT license (Massachusetts Institute of Technology)
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal 
// in the Software without restriction, including without limitation the rights 
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies 
// of the Software, and to permit persons to whom the Software is furnished to 
// do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all 
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING 
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS 
// IN THE SOFTWARE.
// 
//*******************************************************************************/
//
using System;
using LALib.DecMath;

namespace LALib
{

    /// <summary>
    /// Represents a matrix with decimal (128bits) precision floating point values elements.
    /// </summary>
    public class DecimalMatrix : Matrix<Decimal>
    {
        public DecimalMatrix() : base(new DecimalMathOperations())
        {
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
        public static decimal[][] Random(int rows, int cols,
                        decimal minVal, decimal maxVal, int seed)
        {
            // return a matrix with random values
            Random ran = new Random(seed);
            decimal[][] result = Matrix<decimal>.CreateJaggedArray(rows, cols);

            for (int i = 0; i < rows; ++i)
                for (int j = 0; j < cols; ++j)
                    result[i][j] = 
                        NextDecimal(ran, minVal, maxVal);

            return result;
        }


        /// <summary>
        /// Computes the Frobenius norm of a given matrix <paramref name="A"/>
        /// with 128 bits floating point number precision.
        /// </summary>
        /// <returns>Returns the Frobenius norm of the given matrix.</returns>
        /// <param name="A">Source matrix.</param>
        /// <remarks>
        /// In mathematics, a matrix norm is a vector norm in a vector space 
        /// whose elements (vectors) are matrices (of given dimensions). 
        /// </remarks>
        public static decimal FrobeniusNorm128bit(decimal[][] A)
        {
            decimal FResult = decimal.Zero;
            int numRA = Matrix<decimal>.GetNumRows(A);
            int numCA = Matrix<decimal>.GetNumCols(A);

            decimal sumSq = decimal.Zero;// this.Op.Zero; // 0;
            for (int i = 0; i < numRA; i++)
            {
                for (int j = 0; j < numCA; j++)
                {
                    //sumSq +=  A[i][j] * A[i][j];
                    sumSq = decimal.Add( sumSq, decimal.Multiply( A[i][j], A[i][j] ));  // this.Op.Add(sumSq, this.Op.Mult(A[i][j], A[i][j]));
                }
            }

            FResult = DecimalMath.DecimalSqrt(sumSq); // this.Op.SqrtDbl(sumSq);
            return FResult;
        }


        #region Matrix induced infinity Norm


        /// <summary>
        /// Calculates the induced infinity norm of this matrix with 128 bits precision.
        /// </summary>
        /// <returns>Returns the maximum absolute row sum of the matrix.</returns>
        /// <param name="A">Source matrix.</param>
        public static decimal InfinityNorm128bits(decimal[][] A)
        {
            decimal FResult = decimal.MinValue;
            decimal partialRes = decimal.Zero;
            int numR = Matrix<decimal>.GetNumRows(A);
            int numC = Matrix<decimal>.GetNumCols(A);

            for (var i = 0; i < numR; i++)
            {
                decimal sum = decimal.Zero;
                for (var j = 0; j < numC; j++)
                {
                    //sum += Math.Abs( A[i][j] );
                    sum = decimal.Add(sum, Math.Abs(A[i][j]));
                }

                partialRes = (partialRes > sum) ? partialRes : sum;
            }

            return FResult = partialRes;
        }


        #endregion Matrix induced infinity Norm



        #region Matrix P-norms by rows and columns


        /// <summary>
        /// Calculates the p-norms of all row vectors of matrix <paramref name="A"/>
        /// with 128 bits precision.
        /// </summary>
        /// <returns>Returns the p-norms of all row vectors.</returns>
        /// <param name="A">Source matrix.</param>
        /// <param name="p">Norm type (1, 2, infinity).</param>
        /// <remarks>
        /// Typical values for p are 1.0 (L1, Manhattan norm), 2.0 
        /// (L2, Euclidean norm) and positive infinity (infinity norm).
        /// </remarks>
        public decimal[] P_NormByRows128bits(decimal[][] A, double p)
        {
            decimal[] FResult = null;
            int numR = Matrix<decimal>.GetNumRows(A);
            int numC = Matrix<decimal>.GetNumCols(A);

            if (p <= 0.0)
            {
                throw new ArgumentOutOfRangeException(nameof(p), "Value must be positive.");
            }

            var ret = new decimal[numR];
            decimal s;

            if (Utils.AlmostEqual(p, 2.0))
            {
                for (int i = 0; i < numR; i++)
                {
                    s = decimal.Zero;
                    for (int j = 0; j < numC; j++)
                    {
                        s = decimal.Add(s, decimal.Multiply(A[i][j], A[i][j]));
                    }

                    ret[i] = DecimalMath.DecimalSqrt(s);// Math.Sqrt(s);
                }
            }
            else if (Utils.AlmostEqual(p, 1.0))
            {
                // which is simply the sum of each row absolute values of the matrix; 
                for (int i = 0; i < numR; i++)
                {
                    s = 0;
                    for (int j = 0; j < numC; j++)
                    {
                        s = decimal.Add(s, Math.Abs(A[i][j]));
                        //s += this.Op.ToDouble(this.Op.Abs(A[i][j]));
                    }

                    ret[i] = s;
                }
            }
            else if (double.IsPositiveInfinity(p))
            {
                // which is simply the maximum absolute row value of the matrix; 
                for (int i = 0; i < numR; i++)
                {
                    s = 0;
                    for (int j = 0; j < numC; j++)
                    {
                        s = Math.Max(s, Math.Abs(A[i][j]));
                        //s = Math.Max(s, this.Op.ToDouble(this.Op.Abs(A[i][j])));
                    }

                    ret[i] = s;
                }
            }
            else
            {
                decimal decP = (decimal)p;
                decimal invnorm = decimal.Divide(decimal.One, decP);  // 1.0 / p;

                for (int i = 0; i < numR; i++)
                {
                    s = 0;
                    for (int j = 0; j < numC; j++)
                    {
                        s = Decimal.Add(s, DecimalMath.DecimalPow(Math.Abs(A[i][j]), decP));
                        //s += Math.Pow(this.Op.ToDouble(this.Op.Abs(A[i][j])), p);
                    }

                    ret[i] = DecimalMath.DecimalPow(s, invnorm);   //Math.Pow(s, invnorm);
                }
            }

            FResult = ret;
            return FResult;
        }


        /// <summary>
        /// Calculates the p-norms of all column vectors of matrix <paramref name="A"/>
        /// with 128 bits precision.
        /// </summary>
        /// <returns>Returns the p-norms of all column vectors.</returns>
        /// <param name="A">Source matrix.</param>
        /// <param name="p">Norm type (1, 2, infinity).</param>
        /// <remarks>
        /// Typical values for p are 1.0 (L1, Manhattan norm), 2.0 
        /// (L2, Euclidean norm) and positive infinity (infinity norm).
        /// </remarks>
        public decimal[] P_NormByCols128bits(decimal[][] A, double p)
        {
            decimal[] FResult = null;
            int numR = Matrix<decimal>.GetNumRows(A);
            int numC = Matrix<decimal>.GetNumCols(A);

            if (p <= 0.0)
            {
                throw new ArgumentOutOfRangeException(nameof(p), "Value must be positive (1, 2, +infinity).");
            }

            var ret = new decimal[numC];
            decimal s;

            if (Utils.AlmostEqual(p, 2.0))
            {
                for (int j = 0; j < numC; j++)
                {
                    s = decimal.Zero;   //0;
                    for (int i = 0; i < numR; i++)
                    {
                        s = decimal.Add(s, Decimal.Multiply( A[i][j], A[i][j]) );
                        //s += Math.Pow(this.Op.ToDouble(this.Op.Abs(A[i][j])), 2);
                    }

                    ret[j] = DecimalMath.DecimalSqrt(s);   //Math.Sqrt(s);
                }
            }
            else if (Utils.AlmostEqual(p, 1.0))
            {
                // which is simply the sum of each column absolute values of the matrix; 
                for (int j = 0; j < numC; j++)
                {
                    s = decimal.Zero;   // 0;
                    for (int i = 0; i < numR; i++)
                    {
                        s = Decimal.Add(s, Math.Abs(A[i][j]));
                        //s += this.Op.ToDouble(this.Op.Abs(A[i][j]));
                    }

                    ret[j] = s;
                }
            }
            else if (double.IsPositiveInfinity(p))
            {
                // which is simply the maximum absolute column value of the matrix; 
                for (int j = 0; j < numC; j++)
                {
                    s = decimal.Zero;   // 0;
                    for (int i = 0; i < numR; i++)
                    {
                        s = Math.Max(s, Math.Abs(A[i][j]));
                        //s = Math.Max(s, this.Op.ToDouble(this.Op.Abs(A[i][j])));
                    }

                    ret[j] = s;
                }
            }
            else
            {
                decimal decP = (decimal)p;
                decimal invnorm = decimal.Divide(decimal.One, decP);  //1.0 / p;

                for (int j = 0; j < numC; j++)
                {
                    s = decimal.Zero;   //0;
                    for (int i = 0; i < numR; i++)
                    {
                        s = decimal.Add(s, DecimalMath.DecimalPow(Math.Abs(A[i][j]), decP));
                        //s += Math.Pow(this.Op.ToDouble(this.Op.Abs(A[i][j])), p);
                    }

                    ret[j] = DecimalMath.DecimalPow(s, invnorm);  //Math.Pow(s, invnorm);
                }
            }

            FResult = ret;
            return FResult;
        }


        #endregion Matrix P-norms by rows and columns



        /// <summary>
        /// Computes the determinant for decimal (128 bits FPN) matrices using Bareiss alg.
        /// </summary>
        /// <returns>Returns the determinant of matrix.</returns>
        /// <param name="m1">The original matrix.</param>
        /// <param name="n">Matrix dimension.</param>
        /// <param name="overwrite">If <c>true</c> original matrix will be changed.</param>
        /// <remarks>
        ///     Note: The determinant applies only to a square matrix.
        /// </remarks>
        public decimal DET_BareissAlg(decimal[][] m1, int n, bool overwrite)
        {
            decimal FResult = base.DET_BareissAlg_Base_T(m1, n, overwrite);
            return FResult;
        }



        #region Auxiliary methods

        /// <summary>
        /// Computes a random value between 0..1 with 128 bits precision.
        /// </summary>
        /// <returns>Returns the random value.</returns>
        /// <param name="rng">The random generator instance.</param>
        private static decimal NextDecimal(Random rng)
        {
            double RandH, RandL;
            do
            {
                RandH = rng.NextDouble();
                RandL = rng.NextDouble();
            } while ((RandH > 0.99999999999999d) || (RandL > 0.99999999999999d));
            return (decimal)RandH + (decimal)RandL / 1E14m;
        }


        /// <summary>
        /// Computes a random value in a given range with 128 bits precision.
        /// </summary>
        /// <returns>Returns the random value.</returns>
        /// <param name="rng">The random generator instance.</param>
        /// <param name="minValue">The lower range limit.</param>
        /// <param name="maxValue">The upper range limit.</param>
        private static decimal NextDecimal(Random rng, decimal minValue, decimal maxValue)
        {
            return NextDecimal(rng) * (maxValue - minValue) + minValue;
        }

        #endregion Auxiliary methods


    }

}
