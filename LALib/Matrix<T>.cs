/********************************************************************************
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
 * Parts of the code were extracted from the following sources: 
 * 
 *      https://jamesmccaffrey.wordpress.com (James D. McCaffrey blog)
 *      https://en.wikipedia.org/
 *      https://www.geeksforgeeks.org
 *      https://numerics.mathdotnet.com/
 * 
 ********************************************************************************/


using System;
using System.Collections.Generic;
using System.Threading.Tasks;

// using System.Numerics?
#if SYSNUMERICS
using cpx = System.Numerics.Complex;

// using ComplexN library?
#elif COMPLEXN
    using cpx as ComplexN.ComplexNumber;

#else
    throw new exception("Must use a complex number library, i.e., System.Numerics, or ComplexN or other compatible with System.Numerics library.");

#endif

namespace LALib
{

    /// <summary>
    /// Implements the base generic class to execute operations with matrices.
    /// </summary>
    /// <remarks>
    /// This class should be derived from other classes. Do not be used directly, except
    /// to execute static methods.
    /// 
    /// Use 'Int32Matrix' class for matrix of integer elements,
    /// use 'DoubleMatrix' class for matrix of double elements,
    /// use 'ComplexMatrix' class for matrix of complex number elements,
    /// 
    /// 
    ///                         *** IMPORTANT ***
    ///     Matrices here are represented by two dimensional jagged arrays 
    ///     (ex: int[][], double[][] or complex[][])
    /// 
    /// Matrices of integer elements should be used only for limited purposes,
    /// because they only support integer division operations. Decimal values are truncated.
    /// 
    /// Use double or complex number matrices instead for the full range of funtionalities.
    /// 
    /// </remarks>
    public class Matrix<T>
    {

        private readonly string C_MSG_FROB_INNER_PROD = "In order to compute the Frobenius inner product, both " +
                    "matrices must have the same size.";



        #region Constructor

        protected Matrix(IBasicMathOperations<T> op)
        {
            this.Op = op;
        }

        #endregion Constructor



        #region Properties

        /// <summary>
        /// Gets or sets the instance to execute basic math operations like addition,
        /// subtractions, multiplication, ..., of elements of type <typeparamref name="T"/>.
        /// </summary>
        /// <value>The operations instance.</value>
        /// <remarks>Should be set in constructor of derived class.</remarks>
        protected IBasicMathOperations<T> Op { get; set; }


        #endregion Properties



        #region Methods


        /// <summary>
        /// Utility function to initialize a jagged array of any type.
        /// </summary>
        /// <returns>The jagged array.</returns>
        /// <param name="numRows">Number rows.</param>
        /// <param name="numCols">Number columns.</param>
        /// <typeparam name="T">The array base type parameter.</typeparam>
        public static T[][] CreateJaggedArray(int numRows, int numCols)
        {
            T[][] FResult = null;

            FResult = new T[numRows][];
            for (int i = 0; i < numRows; i++)
            {
                FResult[i] = new T[numCols];
            }

            return FResult;
        }


        /// <summary>
        /// Utility function to initialize a jagged array of any type
        /// using parallelization.
        /// </summary>
        /// <returns>The jagged array.</returns>
        /// <param name="numRows">Number rows.</param>
        /// <param name="numCols">Number columns.</param>
        /// <typeparam name="T">The array base type parameter.</typeparam>
        public static T[][] CreateJaggedArrayPar(int numRows, int numCols)
        {
            T[][] FResult = null;
            FResult = new T[numRows][];
            Parallel.For(0, numRows, delegate (int i) 
                {
                    FResult[i] = new T[numCols];
                }
            );

            return FResult;
        }


        /// <summary>
        /// Gets the number rows of rectangular jagged array <paramref name="A"/>.
        /// </summary>
        /// <returns>Returns the number rows.</returns>
        /// <param name="A">The rectangular jagged array.</param>
        public static int GetNumRows(T[][] A)
        {
            return A.GetLength(0);
        }


        /// <summary>
        /// Gets the number columns of rectangular jagged array <paramref name="A"/>.
        /// </summary>
        /// <returns>Returns the number columns.</returns>
        /// <param name="A">The rectangular jagged array.</param>
        public static int GetNumCols(T[][] A)
        {
            if (A.GetLength(0) == 0)
                return 0;

            return A[0].GetLength(0);
        }


        /// <summary>
        /// Clones the matrix (hard copy).
        /// </summary>
        /// <returns>The new matri.</returns>
        /// <param name="m1">The source matrix to be cloned.</param>
        public static T[][] Clone(T[][] m1)
        {
            T[][] FResult = null;
            int rows = m1.GetLength(0);
            int cols = 0;
            if (rows > 0)
            {
                cols = m1[0].GetLength(0);
                FResult = new T[rows][];

                for (int i = 0; i < rows; i++)
                {
                    FResult[i] = new T[cols];
                    for (int j = 0; j < cols; j++)
                    {
                        FResult[i][j] = m1[i][j];
                    }
                }
            }

            return FResult;
        }


        /// <summary>
        /// Clones the matrix using paralellism (hard copy).
        /// </summary>
        /// <returns>The new matrix.</returns>
        /// <param name="m1">The source matrix to be cloned.</param>
        public static T[][] ClonePar(T[][] m1)
        {
            T[][] FResult = null;
            int rows = m1.GetLength(0);
            int cols = 0;
            if (rows > 0)
            {
                cols = m1[0].GetLength(0);
                FResult = new T[rows][];

                Parallel.For(0, rows, delegate(int i) 
                    {
                        FResult[i] = new T[cols];
                        for (int j = 0; j < cols; j++)
                        {
                            FResult[i][j] = m1[i][j];
                        }
                    }
                );
            }

            return FResult;
        }


        /// <summary>
        /// Computes the values sum of each matrix column.
        /// </summary>
        /// <returns>Returns array with sum of columns values.</returns>
        /// <param name="A">The source matrix.</param>
        public T[] SumCols(T[][] A)
        {
            int R = Matrix<T>.GetNumRows(A);
            int C = Matrix<T>.GetNumCols(A);
            T[] FResult = new T[C];

            for (int i = 0; i < C; i++)
            {
                for (int j = 0; j < R; j++)
                {
                    FResult[i] = this.Op.Add( FResult[i], A[j][i] );
                }
            }

            return FResult;
        }


        /// <summary>
        /// Computes the values sum of each matrix row.
        /// </summary>
        /// <returns>Returns array with sum of rows values.</returns>
        /// <param name="A">The source matrix.</param>
        public T[] SumRows(T[][] A)
        {
            int R = Matrix<T>.GetNumRows(A);
            int C = Matrix<T>.GetNumCols(A);
            T[] FResult = new T[R];

            for (int i = 0; i < R; i++)
            {
                for (int j = 0; j < C; j++)
                {
                    FResult[i] = this.Op.Add( FResult[i], A[i][j]);
                }
            }

            return FResult;
        }


        ///// <summary>
        ///// Computes the values sum of each matrix column.
        ///// </summary>
        ///// <returns>Returns array with sum of columns values.</returns>
        ///// <param name="A">The source matrix.</param>
        //public static double[] SumCols(double[][] A)
        //{
        //    int R = Matrix<T>.GetNumRows<double>(A);
        //    int C = Matrix<T>.GetNumCols<double>(A);
        //    double[] FResult = new double[C];

        //    for (int i = 0; i < C; i++)
        //    {
        //        for (int j = 0; j < R; j++)
        //        {
        //            FResult[i] += A[j][i];
        //        }
        //    }

        //    return FResult;
        //}


        ///// <summary>
        ///// Computes the values sum of each matrix row.
        ///// </summary>
        ///// <returns>Returns array with sum of rows values.</returns>
        ///// <param name="A">The source matrix.</param>
        //public static double[] SumRows(double[][] A)
        //{
        //    int R = Matrix<T>.GetNumRows<double>(A);
        //    int C = Matrix<T>.GetNumCols<double>(A);
        //    double[] FResult = new double[R];

        //    for (int i = 0; i < R; i++)
        //    {
        //        for (int j = 0; j < C; j++)
        //        {
        //            FResult[i] += A[i][j];
        //        }
        //    }

        //    return FResult;
        //}


        ///// <summary>
        ///// Computes the values sum of each matrix column.
        ///// </summary>
        ///// <returns>Returns array with sum of columns values.</returns>
        ///// <param name="A">The source matrix.</param>
        //public static cpx[] SumCols(cpx[][] A)
        //{
        //    int R = Matrix<T>.GetNumRows<cpx>(A);
        //    int C = Matrix<T>.GetNumCols<cpx>(A);
        //    cpx[] FResult = new cpx[C];

        //    for (int i = 0; i < C; i++)
        //    {
        //        for (int j = 0; j < R; j++)
        //        {
        //            FResult[i] += A[j][i];
        //        }
        //    }

        //    return FResult;
        //}


        ///// <summary>
        ///// Computes the values sum of each matrix row.
        ///// </summary>
        ///// <returns>Returns array with sum of rows values.</returns>
        ///// <param name="A">The source matrix.</param>
        //public static cpx[] SumRows(cpx[][] A)
        //{
        //    int R = Matrix<T>.GetNumRows<cpx>(A);
        //    int C = Matrix<T>.GetNumCols<cpx>(A);
        //    cpx[] FResult = new cpx[R];

        //    for (int i = 0; i < R; i++)
        //    {
        //        for (int j = 0; j < C; j++)
        //        {
        //            FResult[i] += A[i][j];
        //        }
        //    }

        //    return FResult;
        //}


        /// <summary>
        /// Computes the values sum of each matrix column using parallelism.
        /// </summary>
        /// <returns>Returns array with sum of columns values.</returns>
        /// <param name="A">The source matrix.</param>
        public T[] SumColsPar(T[][] A)
        {
            int R = Matrix<T>.GetNumRows(A);
            int C = Matrix<T>.GetNumCols(A);
            T[] FResult = new T[C];

            Parallel.For(0, C, delegate(int i) {
                    for (int j = 0; j < R; j++)
                    {
                        FResult[i] = this.Op.Add( FResult[i], A[j][i] );
                    }
                }
            );

            return FResult;
        }


        /// <summary>
        /// Computes the values sum of each matrix row using parallelism.
        /// </summary>
        /// <returns>Returns array with sum of rows values.</returns>
        /// <param name="A">The source matrix.</param>
        public T[] SumRowsPar(T[][] A)
        {
            int R = Matrix<T>.GetNumRows(A);
            int C = Matrix<T>.GetNumCols(A);
            T[] FResult = new T[R];

            Parallel.For(0, R, delegate (int i) {
                    for (int j = 0; j < C; j++)
                    {
                        FResult[i] = this.Op.Add( FResult[i], A[i][j] );
                    }
                }
            );

            return FResult;
        }


        ///// <summary>
        ///// Computes the values sum of each matrix column using parallelism.
        ///// </summary>
        ///// <returns>Returns array with sum of columns values.</returns>
        ///// <param name="A">The source matrix.</param>
        //public static double[] SumColsPar(double[][] A)
        //{
        //    int R = Matrix<T>.GetNumRows<double>(A);
        //    int C = Matrix<T>.GetNumCols<double>(A);
        //    double[] FResult = new double[C];

        //    Parallel.For(0, C, delegate (int i) {
        //            for (int j = 0; j < R; j++)
        //            {
        //                FResult[i] += A[j][i];
        //            }
        //        }
        //    );

        //    return FResult;
        //}


        ///// <summary>
        ///// Computes the values sum of each matrix row using parallelism.
        ///// </summary>
        ///// <returns>Returns array with sum of rows values.</returns>
        ///// <param name="A">The source matrix.</param>
        //public static double[] SumRowsPar(double[][] A)
        //{
        //    int R = Matrix<T>.GetNumRows<double>(A);
        //    int C = Matrix<T>.GetNumCols<double>(A);
        //    double[] FResult = new double[R];

        //    Parallel.For(0, R, delegate (int i) {
        //            for (int j = 0; j < C; j++)
        //            {
        //                FResult[i] += A[i][j];
        //            }
        //        }
        //    );

        //    return FResult;
        //}


        ///// <summary>
        ///// Computes the values sum of each matrix column using parallelism.
        ///// </summary>
        ///// <returns>Returns array with sum of columns values.</returns>
        ///// <param name="A">The source matrix.</param>
        //public static cpx[] SumColsPar(cpx[][] A)
        //{
        //    int R = Matrix<T>.GetNumRows<cpx>(A);
        //    int C = Matrix<T>.GetNumCols<cpx>(A);
        //    cpx[] FResult = new cpx[C];

        //    Parallel.For(0, C, delegate (int i) {
        //            for (int j = 0; j < R; j++)
        //            {
        //                FResult[i] += A[j][i];
        //            }
        //        }
        //    );

        //    return FResult;
        //}


        ///// <summary>
        ///// Computes the values sum of each matrix row using parallelism.
        ///// </summary>
        ///// <returns>Returns array with sum of rows values.</returns>
        ///// <param name="A">The source matrix.</param>
        //public static cpx[] SumRowsPar(cpx[][] A)
        //{
        //    int R = Matrix<T>.GetNumRows<cpx>(A);
        //    int C = Matrix<T>.GetNumCols<cpx>(A);
        //    cpx[] FResult = new cpx[R];

        //    Parallel.For(0, R, delegate (int i) {
        //            for (int j = 0; j < C; j++)
        //            {
        //                FResult[i] += A[i][j];
        //            }
        //        }
        //    );

        //    return FResult;
        //}


        /// <summary>
        /// Checks if the matrix product between <paramref name="m1"/> and <paramref name="m2"/> is defined.
        /// </summary>
        /// <returns>Returns <c>true</c>, if matrix product is defined, <c>false</c> otherwise.</returns>
        /// <param name="m1">Left matrix.</param>
        /// <param name="m2">Right matrix.</param>
        /// <typeparam name="T">The 1st type parameter.</typeparam>
        /// <remarks>
        /// In order for the matrix product to be defined,
        /// the number of columns of the left matrix must match the number of rows of right matrix.
        /// Note: In matrix product the commutative property is not applicable.
        /// </remarks>
        public static bool IsMatrixMultDefined(T[][] m1, T[][] m2)
        {
            return (GetNumCols(m1) == GetNumRows(m2));
        }


        ////public abstract void Set(T a, T value);
        //public abstract T Div(T a, T b);
        //public abstract T Mult(T a, T b);
        //public abstract T Sub(T a, T b);
        //public abstract T Add(T a, T b);
        //public abstract T Conjugate(T a);



        /// <summary>
        /// Multiplies two matrices using single core.
        /// </summary>
        /// <param name="numR1">Number of rows of first matrix.</param>
        /// <param name="numC1">Number of columns of first matrix.</param>
        /// <param name="m1">First matrix.</param>
        /// <param name="numR2">Number of rows of second matrix.</param>
        /// <param name="numC2">Number of columns of second matrix.</param>
        /// <param name="m2">Second matrix.</param>
        /// <remarks>
        /// In order for matrix multiplication to be defined, the number of columns in the 
        /// first matrix must be equal to the number of rows in the second matrix.
        /// Result matrix will have number of rows of first matrix and number of columns of second matris.
        /// i.e. let A be (nxm) ans B (mxp), result A.B will be C (nxp). 
        /// </remarks>
        protected T[][] Mult(int numR1, int numC1, T[][] m1, int numR2, int numC2, T[][] m2)
        {
            T[][] FResult = null;

            // check sizes
            if ((numC1 > 0) && (numR1 > 0) && (numC2 > 0) && (numR2 > 0))
            {
                if (numC1 != numR2)
                    throw new ArgumentException("In order for matrix multiplication to be defined, the number of columns in the first matrix must be equal to the number of rows in the second matrix.");
                else
                {
                    // initialize result matrix
                    FResult = Matrix<T>.CreateJaggedArray(numR1, numC2);

                    for (int i = 0; i < numR1; i++)
                    {
                        for (int j = 0; j < numC2; j++)
                        {
                            //FResult[i][j] = 0;
                            for (int k = 0; k < numC1; k++)
                            {
                               FResult[i][j] = this.Op.Add( FResult[i][j], this.Op.Mult( m1[i][k], m2[k][j] ) );
                            }
                        }
                    }
                }
            }
            else
            {
                throw new ArgumentException("In matrix multiplication, the number of rows and columns must be greater than zero in both matrices.");
            }

            return FResult;
        }



        public T[][] Mult(T[][] m1, T[][] m2)
        {
            return Mult(GetNumRows(m1), GetNumCols(m1), m1, GetNumRows(m2), GetNumCols(m2), m2);
        }


        /// <summary>
        /// Multiplies two matrices using parallelism (multiple cores).
        /// </summary>
        /// <param name="numR1">Number of rows of first matrix.</param>
        /// <param name="numC1">Number of columns of first matrix.</param>
        /// <param name="m1">First matrix.</param>
        /// <param name="numR2">Number of rows of second matrix.</param>
        /// <param name="numC2">Number of columns of second matrix.</param>
        /// <param name="m2">Second matrix.</param>
        /// <remarks>
        /// In order for matrix multiplication to be defined, the number of columns in the 
        /// first matrix must be equal to the number of rows in the second matrix.
        /// Result matrix will have number of rows of first matrix and number of columns of second matris.
        /// i.e. let A be (nxm) ans B (mxp), result A.B will be C (nxp). 
        /// </remarks>
        private T[][] MultPar(int numR1, int numC1, T[][] m1, int numR2, int numC2, T[][] m2)
        {
            T[][] FResult = null;

            // check sizes
            if ((numC1 > 0) && (numR1 > 0) && (numC2 > 0) && (numR2 > 0))
            {
                if (numC1 != numR2)
                    throw new ArgumentException("In order for matrix multiplication to be defined, the number of columns in the first matrix must be equal to the number of rows in the second matrix.");
                else
                {
                    // initialize result matrix
                    FResult = Matrix<T>.CreateJaggedArray(numR1, numC2);

                    Parallel.For(0, numR1, delegate(int i)
                        {
                            for (int j = 0; j < numC2; j++)
                            {
                                //FResult[i][j] = 0;
                                for (int k = 0; k < numC1; k++)
                                {
                                    FResult[i][j] = this.Op.Add( FResult[i][j], this.Op.Mult( m1[i][k], m2[k][j]) );
                                }
                            }
                        }
                    );
                }
            }
            else
            {
                throw new ArgumentException("In matrix multiplication, the number of rows and columns must be greater than zero in both matrices.");
            }

            return FResult;
        }


        public T[][] MultPar(T[][] m1, T[][] m2)
        {
            return MultPar(GetNumRows(m1), GetNumCols(m1), m1, GetNumRows(m2), GetNumCols(m2), m2);
        }


        ///// <summary>
        ///// Multiplies two matrices using parallelism (multiple cores).
        ///// </summary>
        ///// <param name="numR1">Number of rows of first matrix.</param>
        ///// <param name="numC1">Number of columns of first matrix.</param>
        ///// <param name="m1">First matrix.</param>
        ///// <param name="numR2">Number of rows of second matrix.</param>
        ///// <param name="numC2">Number of columns of second matrix.</param>
        ///// <param name="m2">Second matrix.</param>
        ///// <remarks>
        ///// In order for matrix multiplication to be defined, the number of columns in the 
        ///// first matrix must be equal to the number of rows in the second matrix.
        ///// Result matrix will have number of rows of first matrix and number of columns of second matris.
        ///// i.e. let A be (nxm) ans B (mxp), result A.B will be C (nxp). 
        ///// </remarks>
        //public static double[][] MultPar(int numR1, int numC1, double[][] m1, int numR2, int numC2, double[][] m2)
        //{
        //    double[][] FResult = null;

        //    // check sizes
        //    if ((numC1 > 0) && (numR1 > 0) && (numC2 > 0) && (numR2 > 0))
        //    {
        //        if (numC1 != numR2)
        //            throw new ArgumentException("In order for matrix multiplication to be defined, the number of columns in the first matrix must be equal to the number of rows in the second matrix.");
        //        else
        //        {
        //            // initialize result matrix
        //            FResult = Matrix<T>.CreateJaggedArray<double>(numR1, numC2);

        //            Parallel.For(0, numR1, delegate (int i)
        //            {
        //                for (int j = 0; j < numC2; j++)
        //                {
        //                    FResult[i][j] = 0;
        //                    for (int k = 0; k < numC1; k++)
        //                    {
        //                        FResult[i][j] += m1[i][k] * m2[k][j];
        //                    }
        //                }
        //            }
        //                );
        //        }
        //    }
        //    else
        //    {
        //        throw new ArgumentException("In matrix multiplication, the number of rows and columns must be greater than zero in both matrices.");
        //    }

        //    return FResult;
        //}


        ///// <summary>
        ///// Multiplies two matrices using single core.
        ///// </summary>
        ///// <param name="numR1">Number of rows of first matrix.</param>
        ///// <param name="numC1">Number of columns of first matrix.</param>
        ///// <param name="m1">First matrix.</param>
        ///// <param name="numR2">Number of rows of second matrix.</param>
        ///// <param name="numC2">Number of columns of second matrix.</param>
        ///// <param name="m2">Second matrix.</param>
        ///// <remarks>
        ///// In order for matrix multiplication to be defined, the number of columns in the 
        ///// first matrix must be equal to the number of rows in the second matrix.
        ///// Result matrix will have number of rows of first matrix and number of columns of second matris.
        ///// i.e. let A be (nxm) ans B (mxp), result A.B will be C (nxp). 
        ///// </remarks>
        //public static double[][] Mult(int numR1, int numC1, double[][] m1, int numR2, int numC2, double[][] m2)
        //{
        //    double[][] FResult = null;

        //    // check sizes
        //    if ((numC1 > 0) && (numR1 > 0) && (numC2 > 0) && (numR2 > 0))
        //    {
        //        if (numC1 != numR2)
        //            throw new ArgumentException("In order for matrix multiplication to be defined, the number of columns in the first matrix must be equal to the number of rows in the second matrix.");
        //        else
        //        {
        //            // initialize result matrix
        //            FResult = Matrix<T>.CreateJaggedArray<double>(numR1, numC2);

        //            for (int i = 0; i < numR1; i++)
        //            {
        //                for (int j = 0; j < numC2; j++)
        //                {
        //                    FResult[i][j] = 0;
        //                    for (int k = 0; k < numC1; k++)
        //                    {
        //                        FResult[i][j] += m1[i][k] * m2[k][j];
        //                    }
        //                }
        //            }
        //        }
        //    }
        //    else
        //    {
        //        throw new ArgumentException("In matrix multiplication, the number of rows and columns must be greater than zero in both matrices.");
        //    }

        //    return FResult;
        //}


        ///// <summary>
        ///// Multiplies two matrices using parallelism (multiple cores).
        ///// </summary>
        ///// <param name="numR1">Number of rows of first matrix.</param>
        ///// <param name="numC1">Number of columns of first matrix.</param>
        ///// <param name="m1">First matrix.</param>
        ///// <param name="numR2">Number of rows of second matrix.</param>
        ///// <param name="numC2">Number of columns of second matrix.</param>
        ///// <param name="m2">Second matrix.</param>
        ///// <remarks>
        ///// In order for matrix multiplication to be defined, the number of columns in the 
        ///// first matrix must be equal to the number of rows in the second matrix.
        ///// Result matrix will have number of rows of first matrix and number of columns of second matris.
        ///// i.e. let A be (nxm) ans B (mxp), result A.B will be C (nxp). 
        ///// </remarks>
        //public static cpx[][] MultPar(int numR1, int numC1, cpx[][] m1, int numR2, int numC2, cpx[][] m2)
        //{
        //    cpx[][] FResult = null;

        //    // check sizes
        //    if ((numC1 > 0) && (numR1 > 0) && (numC2 > 0) && (numR2 > 0))
        //    {
        //        if (numC1 != numR2)
        //            throw new ArgumentException("In order for matrix multiplication to be defined, the number of columns in the first matrix must be equal to the number of rows in the second matrix.");
        //        else
        //        {
        //            // initialize result matrix
        //            FResult = Matrix<T>.CreateJaggedArray<cpx>(numR1, numC2);

        //            Parallel.For(0, numR1, delegate (int i)
        //            {
        //                for (int j = 0; j < numC2; j++)
        //                {
        //                    FResult[i][j] = 0;
        //                    for (int k = 0; k < numC1; k++)
        //                    {
        //                        FResult[i][j] += m1[i][k] * m2[k][j];
        //                    }
        //                }
        //            }
        //                );
        //        }
        //    }
        //    else
        //    {
        //        throw new ArgumentException("In matrix multiplication, the number of rows and columns must be greater than zero in both matrices.");
        //    }

        //    return FResult;
        //}


        ///// <summary>
        ///// Multiplies two matrices using single core.
        ///// </summary>
        ///// <param name="numR1">Number of rows of first matrix.</param>
        ///// <param name="numC1">Number of columns of first matrix.</param>
        ///// <param name="m1">First matrix.</param>
        ///// <param name="numR2">Number of rows of second matrix.</param>
        ///// <param name="numC2">Number of columns of second matrix.</param>
        ///// <param name="m2">Second matrix.</param>
        ///// <remarks>
        ///// In order for matrix multiplication to be defined, the number of columns in the 
        ///// first matrix must be equal to the number of rows in the second matrix.
        ///// Result matrix will have number of rows of first matrix and number of columns of second matris.
        ///// i.e. let A be (nxm) ans B (mxp), result A.B will be C (nxp). 
        ///// </remarks>
        //public static cpx[][] Mult(int numR1, int numC1, cpx[][] m1, int numR2, int numC2, cpx[][] m2)
        //{
        //    cpx[][] FResult = null;

        //    // check sizes
        //    if ((numC1 > 0) && (numR1 > 0) && (numC2 > 0) && (numR2 > 0))
        //    {
        //        if (numC1 != numR2)
        //            throw new ArgumentException("In order for matrix multiplication to be defined, the number of columns in the first matrix must be equal to the number of rows in the second matrix.");
        //        else
        //        {
        //            // initialize result matrix
        //            FResult = Matrix<T>.CreateJaggedArray<cpx>(numR1, numC2);

        //            for (int i = 0; i < numR1; i++)
        //            {
        //                for (int j = 0; j < numC2; j++)
        //                {
        //                    FResult[i][j] = 0;
        //                    for (int k = 0; k < numC1; k++)
        //                    {
        //                        FResult[i][j] += m1[i][k] * m2[k][j];
        //                    }
        //                }
        //            }
        //        }
        //    }
        //    else
        //    {
        //        throw new ArgumentException("In matrix multiplication, the number of rows and columns must be greater than zero in both matrices.");
        //    }

        //    return FResult;
        //}


        /// <summary>
        /// Transpose the specified matrix <paramref name="m1"/> and returns result in new matrix.
        /// </summary>
        /// <returns>The resulting transposed matrix.</returns>
        /// <param name="numR">Number rows of source matrix.</param>
        /// <param name="numC">Number columns of source matrix.</param>
        /// <param name="m1">Source matrix.</param>
        /// <typeparam name="T">The matrix base element type.</typeparam>
        public static T[][] TransposePar(int numR, int numC, T[][] m1)
        {
            T[][] FResult = null;
            FResult = Matrix<T>.CreateJaggedArray(numC, numR);

            Parallel.For(0, numR, delegate (int i)
                {
                    for (int j = 0; j < numC; j++)
                    {
                        FResult[j][i] = m1[i][j];
                    }
                }
            );

            return FResult;
        }


        /// <summary>
        /// Transpose the specified matrix <paramref name="m1"/> and returns result in new matrix.
        /// </summary>
        /// <returns>The resulting transposed matrix.</returns>
        /// <param name="numR">Number rows of source matrix.</param>
        /// <param name="numC">Number columns of source matrix.</param>
        /// <param name="m1">Source matrix.</param>
        /// <typeparam name="T">The matrix base element type.</typeparam>
        public static T[][] Transpose(int numR, int numC, T[][] m1)
        {
            T[][] FResult = null;
            FResult = Matrix<T>.CreateJaggedArray(numC, numR);

            for(int i = 0; i < numR; i++)
            {
                for (int j = 0; j < numC; j++)
                {
                    FResult[j][i] = m1[i][j];
                }
            }

            return FResult;
        }


        /// <summary>
        /// Computes the conjugate transpose of the specified matrix <paramref name="m1"/> and returns result in new matrix.
        /// </summary>
        /// <returns>Returns the resulting conjugate transposed matrix.</returns>
        /// <param name="numR">Number rows of source matrix.</param>
        /// <param name="numC">Number columns of source matrix.</param>
        /// <param name="m1">Source matrix.</param>
        /// <typeparam name="T">The matrix base element type.</typeparam>
        public T[][] ConjugateTranspose(int numR, int numC, T[][] m1)
        {
            T[][] FResult = null;
            FResult = Matrix<T>.CreateJaggedArray(numC, numR);

            for (int i = 0; i < numR; i++)
            {
                for (int j = 0; j < numC; j++)
                {
                    FResult[j][i] = this.Op.Conjugate( m1[i][j] );
                }
            }

            return FResult;
        }

        //public static double DET(int[][] a, int n)
        //{
        //    int i, j, k;
        //    double det = 0;
        //    for (i = 0; i < n - 1; i++)
        //    {
        //        for (j = i + 1; j < n; j++)
        //        {
        //            det = a[j][i] / a[i][i];
        //            for (k = i; k < n; k++)
        //                a[j][k] = a[j][k] - det * a[i][k]; // Here's exception
        //        }
        //    }
        //    det = 1;
        //    for (i = 0; i < n; i++)
        //        det = det * a[i][i];
        //    return det;
        //}


        //public static double DETInPlace(double[][] a, int n)
        //{
        //    int i, j, k;
        //    double det = 0;
        //    for (i = 0; i < n - 1; i++)
        //    {
        //        for (j = i + 1; j < n; j++)
        //        {
        //            det = a[j][i] / a[i][i];
        //            for (k = i; k < n; k++)
        //                a[j][k] = a[j][k] - det * a[i][k]; // Here's exception
        //        }
        //    }
        //    det = 1;
        //    for (i = 0; i < n; i++)
        //        det = det * a[i][i];
        //    return det;
        //}


        /// <summary>
        /// Computes the determinant for double matrices using Bareiss alg.
        /// </summary>
        /// <returns>Returns the determinant of matrix.</returns>
        /// <param name="m1">The original matrix.</param>
        /// <param name="n">Matrix<T> dimension.</param>
        /// <param name="overwrite">If <c>true</c> original matrix will be changed.</param>
        /// <remarks>
        ///     Note: The determinant applies only to a square matrix.
        /// 
        ///     --- This method should only be called from a derived class. ---
        /// 
        /// </remarks>
        protected T DET_BareissAlg_Base_T(T[][] m1, int n, bool overwrite)
        {
            T FResult;
            T[][] mat = m1;

            if (!overwrite)
                mat = Matrix<T>.Clone(m1);

            mat[0][0] = this.Op.One;

            for (int k = 1; k < n; k++)
            {
                for (int i = k + 1; i < n; i++)
                {
                    for (int j = k + 1; j < n; j++)
                    {
                        mat[i][j] = this.Op.Div( this.Op.Sub( this.Op.Mult( mat[i][j], mat[k][k] ), this.Op.Mult( mat[i][k], mat[k][j] ) ), mat[k - 1][k - 1] );
                    }
                }
            }

            FResult = mat[n-1][n-1];
            return FResult;
        }


        ///// <summary>
        ///// Computes the determinant for integer matrices using Bareiss alg.
        ///// </summary>
        ///// <returns>Returns the determinant of matrix.</returns>
        ///// <param name="m1">The original matrix.</param>
        ///// <param name="n">Matrix<T> dimension.</param>
        ///// <param name="overwrite">If <c>true</c> original matrix will be changed.</param>
        ///// <remarks>
        /////     Note: The determinant applies only to a square matrix.
        ///// </remarks>
        //public static double DET_BareissAlg(int[][] m1, int n, bool overwrite)
        //{
        //    double FResult = double.NaN;
        //    int[][] mat = m1;

        //    if (!overwrite)
        //        mat = Matrix<T>.Clone<int>(m1);

        //    mat[0][0] = 1;

        //    for (int k = 1; k < n; k++)
        //    {
        //        for (int i = k + 1; i < n; i++)
        //        {
        //            for (int j = k + 1; j < n; j++)
        //            {
        //                mat[i][j] = (mat[i][j] * mat[k][k] - mat[i][k] * mat[k][j]) / mat[k - 1][k - 1];
        //            }
        //        }
        //    }

        //    FResult = mat[n - 1][n - 1];
        //    return FResult;
        //}


        ///// <summary>
        ///// Computes the determinant for complex number matrices using Bareiss alg.
        ///// </summary>
        ///// <returns>Returns the determinant of matrix.</returns>
        ///// <param name="m1">The original matrix.</param>
        ///// <param name="n">Matrix<T> dimension.</param>
        ///// <param name="overwrite">If <c>true</c> original matrix will be changed.</param>
        ///// <remarks>
        /////     Note: The determinant applies only to a square matrix.
        ///// </remarks>
        //public static cpx DET_BareissAlg(cpx[][] m1, int n, bool overwrite)
        //{
        //    cpx FResult = double.NaN;
        //    cpx[][] mat = m1;

        //    if (!overwrite)
        //        mat = Matrix<T>.Clone<cpx>(m1);

        //    mat[0][0] = 1;

        //    for (int k = 1; k < n; k++)
        //    {
        //        for (int i = k + 1; i < n; i++)
        //        {
        //            for (int j = k + 1; j < n; j++)
        //            {
        //                mat[i][j] = (mat[i][j] * mat[k][k] - mat[i][k] * mat[k][j]) / mat[k - 1][k - 1];
        //            }
        //        }
        //    }

        //    FResult = mat[n - 1][n - 1];
        //    return FResult;
        //}


        /// <summary>
        /// Computes determinant from given decomposed matrix <paramref name="LUP"/>.
        /// </summary>
        /// <returns>Returns the determinant from given decomposed matrix <paramref name="LUP"/>.</returns>
        /// <param name="LUP">The source decomposed matrix.</param>
        /// <typeparam name="T">The 1st type parameter.</typeparam>
        /// <remarks>
        ///     Note: Applies to integer and double matrices only.
        ///           The determinant applies only to a square matrix.
        /// 
        /// </remarks>
        public T DeterminanFromLU(T[][] LUP)
        {
            T FResult = (new LUSolver<T>(this.Op)).DeterminantFromLU(LUP);
            return FResult;
        }


        ///// <summary>
        ///// Computes determinant from given complex number decomposed matrix <paramref name="LUP"/>.
        ///// </summary>
        ///// <returns>Returns the determinant from given decomposed matrix <paramref name="LUP"/>.</returns>
        ///// <param name="LUP">The source complex number decomposed matrix.</param>
        ///// <typeparam name="T">The 1st type parameter.</typeparam>
        ///// <remarks>
        /////     Note: Applies to complex number matrices only.
        /////           The determinant applies only to a square matrix.
        ///// 
        ///// </remarks>
        //public static cpx DeterminanFromLU(cpx[][] LUP)
        //{
        //    cpx FResult = LUSolver<cpx>.DeterminantFromLU(LUP);
        //    return FResult;
        //}

        //    int NRows = numR - fromR + 1;

        //    // Note: must be a square matrix
        //    if (numR != numC)
        //        throw new ArgumentException("In order to calculate determinant the matrix must be square.");

        //    if (NRows == 0)
        //        throw new ArgumentException("To calculate determinant the matrix must have at least one row and column.");
        //    else if (NRows == 1)
        //        return m1[0][0];
        //    else if (NRows == 2)
        //        return m1[fromR][fromC] * m1[fromR + 1][fromC + 1] - m1[fromR][fromC + 1] * m1[fromR + 1][fromC];
        //    else if (NRows == 3) // Leibniz formula 3x3
        //        return m1[fromR][fromC] * m1[fromR + 1][fromC + 1] * m1[fromR + 2][fromC + 2] + m1[fromR][fromC + 1] * m1[fromR + 1][fromC + 2] * m1[fromR][fromC + 2]
        //                    + m1[fromR][fromC + 2] * m1[fromR + 1][fromC] * m1[fromR + 2][fromC + 1]
        //                    - m1[fromR][fromC + 2] * m1[fromR + 1][fromC + 1] * m1[fromR + 2][fromC]
        //                    - m1[fromR][fromC + 1] * m1[fromR + 1][fromC] * m1[fromR + 2][fromC + 2]
        //                    - m1[fromR][fromC] * m1[fromR + 1][fromC + 2] * m1[fromR + 2][fromC + 1];
        //    else
        //    {
        //        int row = 1;
        //        Int64 Ai = 0;
        //        int nCols = numC - fromC + 1;

        //        for (int i = fromC; i < nCols; i++)
        //        {
        //            Ai = Determinant(numR, numC, m1, fromR + 1, fromC + 1);

        //            (m1[row][i] * Math.Pow(-1, row + i)) * v(rowR * (i + 1))
        //        }
        //    }
        //}


        //public static Int64 Determinant(int numR, int numC, int[][] m1, int fromR, int fromC)
        //{
        //    Int64 FResult = 1;
        //    int NRows = numR - fromR + 1;

        //    // Note: must be a square matrix
        //    if (numR != numC)
        //        throw new ArgumentException("In order to calculate determinant the matrix must be square.");

        //    if (NRows == 0)
        //        throw new ArgumentException("To calculate determinant the matrix must have at least one row and column.");
        //    else if (NRows == 1)
        //        return m1[0][0];
        //    else if (NRows == 2)
        //        return m1[fromR][fromC] * m1[fromR + 1][fromC + 1] - m1[fromR][fromC + 1] * m1[fromR + 1][fromC];
        //    else if (NRows== 3) // Leibniz formula 3x3
        //        return m1[fromR][fromC] * m1[fromR + 1][fromC + 1] * m1[fromR + 2][fromC + 2] + m1[fromR][fromC + 1] * m1[fromR + 1][fromC + 2] * m1[fromR][fromC + 2]
        //                    + m1[fromR][fromC + 2] * m1[fromR + 1][fromC] * m1[fromR + 2][fromC + 1]
        //                    - m1[fromR][fromC + 2] * m1[fromR + 1][fromC + 1] * m1[fromR + 2][fromC]
        //                    - m1[fromR][fromC + 1] * m1[fromR + 1][fromC] * m1[fromR + 2][fromC + 2]
        //                    - m1[fromR][fromC] * m1[fromR + 1][fromC + 2] * m1[fromR + 2][fromC + 1];
        //    else
        //    {
        //        int row = 1;
        //        Int64 Ai = 0;
        //        int nCols = numC - fromC + 1;

        //        for (int i = fromC; i < nCols; i++)
        //        {
        //            Ai = Determinant(numR, numC, m1, fromR + 1, fromC + 1);

        //             (m1[row][i] * Math.Pow(-1, row + i)) *   v (rowR*(i+1))
        //        }
        //    }
        //}


        #region Matrix exponentiation


        // Exponentiation only applies to square matrices

        /// <summary>
        /// Computes a double matrix raised to a given integer power <paramref name="exp"/>.
        /// </summary>
        /// <returns>Returns matrix raised to a given integer power.</returns>
        /// <param name="mat">Source matrix.</param>
        /// <param name="exp">Integer exponent.</param>
        /// <remarks>
        ///   Only applies to square matrices.
        /// </remarks>
        public T[][] Pow(T[][] mat, int exp)
        {
            int numR = mat.GetLength(0);
            if (numR != mat[0].GetLength(0)) throw new ArgumentException("Matrix<T> must be square. Exponentiation can only be applied to square matrices.");
            T[][] accumulator = this.Identity(numR);
            for (int i = 0; i < exp; i++)
            {
                accumulator = this.Mult(numR, numR, accumulator, numR, numR, mat);
            }
            return accumulator;
        }


        ///// <summary>
        ///// Computes an integer matrix raised to a given integer power <paramref name="exp"/>.
        ///// </summary>
        ///// <returns>Returns matrix raised to a given integer power.</returns>
        ///// <param name="mat">Source matrix.</param>
        ///// <param name="exp">Integer exponent.</param>
        ///// <remarks>
        /////   Only applies to square matrices.
        ///// </remarks>
        //public static int[][] Pow(int[][] mat, int exp)
        //{
        //    int numR = mat.GetLength(0);
        //    if (numR != mat[0].GetLength(0)) throw new ArgumentException("Matrix<T> must be square. Exponentiation can only be applied to square matrices.");
        //    int[][] accumulator = Matrix<T>.IdentityInt(numR);
        //    for (int i = 0; i < exp; i++)
        //    {
        //        accumulator = Matrix<T>.Mult (numR, numR, accumulator, numR, numR, mat);
        //    }
        //    return accumulator;
        //}


        ///// <summary>
        ///// Computes a complex number matrix raised to a given integer power <paramref name="exp"/>.
        ///// </summary>
        ///// <returns>Returns matrix raised to a given integer power.</returns>
        ///// <param name="mat">Source matrix.</param>
        ///// <param name="exp">Integer exponent.</param>
        ///// <remarks>
        /////   Only applies to square matrices.
        ///// </remarks>
        //public static cpx[][] Pow(cpx[][] mat, int exp)
        //{
        //    int numR = mat.GetLength(0);
        //    if (numR != mat[0].GetLength(0)) throw new ArgumentException("Matrix<T> must be square. Exponentiation can only be applied to square matrices.");
        //    cpx[][] accumulator = Matrix<T>.IdentityCpx(numR);
        //    for (int i = 0; i < exp; i++)
        //    {
        //        accumulator = Matrix<T>.Mult(numR, numR, accumulator, numR, numR, mat);
        //    }
        //    return accumulator;
        //}


        #endregion Matrix exponentiation



        #region Complex conjugate matrix


        /// <summary>
        /// Computes the conjugate the specified matrix <paramref name="A"/>.
        /// </summary>
        /// <returns>Returns the conjugate matrix.</returns>
        /// <param name="A">The source matrix.</param>
        /// <param name="overwrite">
        /// If set to <c>true</c> overwrites matrix <paramref name="A"/>, 
        /// otherwise 'A' remains unchanged.
        /// 
        /// For non complex matrices, result should be equal to source matrix.
        /// </param>
        public T[][] Conjugate(T[][] A, bool overwrite)
        {
            int numR = A.GetLength(0);
            int numC = A[0].GetLength(0);
            T[][] FResult = null;

            if (overwrite)
                FResult = A;
            else
                FResult = Matrix<T>.CreateJaggedArray(numR, numC);

            for (int i = 0; i < numR; i++)
            {
                for (int j = 0; j < numC; j++)
                {
                    FResult[i][j] = this.Op.Conjugate( A[i][j] );
                }
            }

            return FResult;
        }


        /// <summary>
        /// Computes the conjugate the specified matrix <paramref name="A"/> using parallelism.
        /// </summary>
        /// <returns>Returns the conjugate matrix.</returns>
        /// <param name="A">The source matrix.</param>
        /// <param name="overwrite">
        /// If set to <c>true</c> overwrites matrix <paramref name="A"/>, 
        /// otherwise 'A' remains unchanged.
        /// 
        /// For non complex matrices, result should be equal to source matrix.
        /// </param>
        public T[][] ConjugatePar(T[][] A, bool overwrite)
        {
            int numR = A.GetLength(0);
            int numC = A[0].GetLength(0);
            T[][] FResult = null;

            if (overwrite)
                FResult = A;
            else
                FResult = Matrix<T>.CreateJaggedArray(numR, numC);

            Parallel.For(0, numR, delegate(int i) {
                    for (int j = 0; j < numC; j++)
                    {
                        FResult[i][j] = this.Op.Conjugate(A[i][j]);
                    }
                }
            );

            return FResult;
        }



        ///// <summary>
        ///// Computes the conjugate the specified matrix <paramref name="A"/>.
        ///// </summary>
        ///// <returns>Returns the conjugate matrix.</returns>
        ///// <param name="A">The source matrix.</param>
        ///// <param name="overwrite">
        ///// If set to <c>true</c> overwrites matrix <paramref name="A"/>, 
        ///// otherwise 'A' remains unchanged.
        ///// </param>
        //public static cpx[][] Conjugate(cpx[][] A, bool overwrite)
        //{
        //    int numR = A.GetLength(0);
        //    int numC = A[0].GetLength(0);
        //    cpx[][] FResult = null;

        //    if (overwrite)
        //        FResult = A;
        //    else
        //        FResult = Matrix<T>.CreateJaggedArray<cpx>(numR, numC);

        //    for (int i = 0; i < numR; i++)
        //    {
        //        for (int j = 0; j < numC; j++)
        //        {
        //            FResult[i][j] = cpx.Conjugate(A[i][j]);
        //        }
        //    }

        //    return FResult;
        //}


        #endregion Complex conjugate matrix



        #region Identity matrix


        /// <summary>
        /// Creates the identity matrix of size <paramref name="size"/>.
        /// </summary>
        /// <returns>Returns the identity matrix.</returns>
        /// <param name="size">Size of matrix (number of rows and columns).</param>
        /// It is "square" (has same number of rows as columns),
        /// It has 1s on the diagonal and 0s everywhere else.
        /// Its symbol is the capital letter I.
        public T[][] Identity(int size)
        {
            T[][] FResult = Matrix<T>.CreateJaggedArray(size, size);
            for (int i = 0; i < size; i++) FResult[i][i] = this.Op.One;
            return FResult;
        }



        ///// <summary>
        ///// Creates the identity matrix of size <paramref name="size"/>.
        ///// </summary>
        ///// <returns>Returns the identity matrix.</returns>
        ///// <param name="size">Size of matrix (number of rows and columns).</param>
        ///// It is "square" (has same number of rows as columns),
        ///// It has 1s on the diagonal and 0s everywhere else.
        ///// Its symbol is the capital letter I.
        //public static cpx[][] IdentityCpx(int size)
        //{
        //    cpx[][] FResult = CreateJaggedArray<cpx>(size, size);
        //    for (int i = 0; i < size; i++) FResult[i][i] = 1;
        //    return FResult;
        //}


        ///// <summary>
        ///// Creates the identity matrix of size <paramref name="size"/>.
        ///// </summary>
        ///// <returns>Returns the identity matrix.</returns>
        ///// <param name="size">Size of matrix (number of rows and columns).</param>
        ///// <remarks>
        ///// It is "square" (has same number of rows as columns),
        ///// It has 1s on the diagonal and 0s everywhere else.
        ///// Its symbol is the capital letter I.
        ///// </remarks>
        //public static int[][] IdentityInt(int size)
        //{
        //    int[][] FResult = CreateJaggedArray<int>(size, size);
        //    for (int i = 0; i < size; i++) FResult[i][i] = 1;
        //    return FResult;
        //}


        ///// <summary>
        ///// Creates the identity matrix of size <paramref name="size"/>.
        ///// </summary>
        ///// <returns>Returns the identity matrix.</returns>
        ///// <param name="size">Size of matrix (number of rows and columns).</param>
        ///// <remarks>
        ///// It is "square" (has same number of rows as columns),
        ///// It has 1s on the diagonal and 0s everywhere else.
        ///// Its symbol is the capital letter I.
        ///// </remarks>
        //public static double[][] IdentityDbl(int size)
        //{
        //    double[][] FResult = CreateJaggedArray<double>(size, size);
        //    for (int i = 0; i < size; i++) FResult[i][i] = 1;
        //    return FResult;
        //}


        /// <summary>
        /// Creates the identity matrix of size <paramref name="size"/> using parallelism.
        /// </summary>
        /// <returns>Returns the identity matrix.</returns>
        /// <param name="size">Size of matrix (number of rows and columns).</param>
        /// <remarks>
        /// It is "square" (has same number of rows as columns),
        /// It has 1s on the diagonal and 0s everywhere else.
        /// Its symbol is the capital letter I.
        /// </remarks>
        public T[][] IdentityPar(int size)
        {
            T[][] FResult = Matrix<T>.CreateJaggedArray(size, size);
            Parallel.For(0, size, delegate (int i)
                {
                    FResult[i][i] = this.Op.One;
                }
            );
            return FResult;
        }


        ///// <summary>
        ///// Creates the identity matrix of size <paramref name="size"/> using parallelism.
        ///// </summary>
        ///// <returns>Returns the identity matrix.</returns>
        ///// <param name="size">Size of matrix (number of rows and columns).</param>
        ///// <remarks>
        ///// It is "square" (has same number of rows as columns),
        ///// It has 1s on the diagonal and 0s everywhere else.
        ///// Its symbol is the capital letter I.
        ///// </remarks>
        //public static int[][] IdentityIntPar(int size)
        //{
        //    int[][] FResult = CreateJaggedArray<int>(size, size);
        //    Parallel.For(0, size, delegate (int i)
        //        {
        //            FResult[i][i] = 1;
        //        }
        //    );
        //    return FResult;
        //}


        ///// <summary>
        ///// Creates the identity matrix of size <paramref name="size"/> using parallelism.
        ///// </summary>
        ///// <returns>Returns the identity matrix.</returns>
        ///// <param name="size">Size of matrix (number of rows and columns).</param>
        ///// <remarks>
        ///// It is "square" (has same number of rows as columns),
        ///// It has 1s on the diagonal and 0s everywhere else.
        ///// Its symbol is the capital letter I.
        ///// </remarks>
        //public static cpx[][] IdentityCpxPar(int size)
        //{
        //    cpx[][] FResult = CreateJaggedArray<cpx>(size, size);
        //    Parallel.For(0, size, delegate (int i)
        //        {
        //            FResult[i][i] = 1;
        //        }
        //    );
        //    return FResult;
        //}


        #endregion Identity matrix



        #region Trace


        /// <summary>
        /// Computes the trace of a double matrix.
        /// </summary>
        /// <returns>Returns the trace of a double matrix.</returns>
        /// <param name="A">The source matrix.</param>
        /// <param name="N">The matrix dimension.</param>
        /// <remarks>
        /// In linear algebra, the trace of a square matrix A, denoted tr(A),
        ///  is defined to be the sum of elements on the main diagonal 
        /// (from the upper left to the lower right) of A. The trace is only 
        /// defined for a square matrix (n × n). 
        /// </remarks>
        public T Trace(T[][] A, int N)
        {
            T FResult = this.Op.Zero;// 0.0;
            for (int i = 0; i < N; i++)
            {
                FResult = this.Op.Add( FResult, A[i][i] );
            }

            return FResult;
        }


        ///// <summary>
        ///// Computes the trace of an integer matrix.
        ///// </summary>
        ///// <returns>Returns the trace of an integer matrix.</returns>
        ///// <param name="A">The source matrix.</param>
        ///// <param name="N">The matrix dimension.</param>
        ///// <remarks>
        ///// In linear algebra, the trace of a square matrix A, denoted tr(A),
        /////  is defined to be the sum of elements on the main diagonal 
        ///// (from the upper left to the lower right) of A. The trace is only 
        ///// defined for a square matrix (n × n). 
        ///// </remarks>
        //public static int Trace(int[][] A, int N)
        //{
        //    int FResult = 0;
        //    for (int i = 0; i < N; i++)
        //    {
        //        FResult += A[i][i];
        //    }

        //    return FResult;
        //}


        ///// <summary>
        ///// Computes the trace of a complex matrix.
        ///// </summary>
        ///// <returns>Returns the trace of a complex matrix.</returns>
        ///// <param name="A">The source matrix.</param>
        ///// <param name="N">The matrix dimension.</param>
        ///// <remarks>
        ///// In linear algebra, the trace of a square matrix A, denoted tr(A),
        /////  is defined to be the sum of elements on the main diagonal 
        ///// (from the upper left to the lower right) of A. The trace is only 
        ///// defined for a square matrix (n × n). 
        ///// </remarks>
        //public static cpx Trace(cpx[][] A, int N)
        //{
        //    cpx FResult = 0;
        //    for (int i = 0; i < N; i++)
        //    {
        //        FResult += A[i][i];
        //    }

        //    return FResult;
        //}


        #endregion Trace



        #region Symmetric, skew-Symmetric and Hermitian and skew-Hermitian check


        /// <summary>
        /// Checks if a matrix is symmetric or Hermitian.
        /// </summary>
        /// <returns>Returns <c>true</c>, if is hermitian, <c>false</c> otherwise.</returns>
        /// <param name="A">The source matrix.</param>
        /// <param name="N">The source matrix dimension.</param>
        /// <remarks>
        /// In linear algebra, a symmetric matrix is a square matrix that is 
        /// equal to its transpose. 
        /// </remarks>
        public static bool IsSymmetric(T[][] A, int N)
        {
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    if (!(A[i][j].Equals( A[j][i] ) ) )
                        return false;
                }
            }

            return true;
        }


        ///// <summary>
        ///// Checks if a matrix is symmetric or Hermitian.
        ///// </summary>
        ///// <returns>Returns <c>true</c>, if is hermitian, <c>false</c> otherwise.</returns>
        ///// <param name="A">The source matrix.</param>
        ///// <param name="N">The source matrix dimension.</param>
        ///// <remarks>
        ///// In linear algebra, a symmetric matrix is a square matrix that is 
        ///// equal to its transpose. 
        ///// </remarks>
        //public static bool IsSymmetric(int[][] A, int N)
        //{
        //    for (int i = 0; i < N; i++)
        //    {
        //        for (int j = 0; j < N; j++)
        //        {
        //            if (A[i][j] != A[j][i])
        //                return false;
        //        }
        //    }

        //    return true;
        //}


        /// <summary>
        /// Checks if a matrix is skew-symmetric.
        /// </summary>
        /// <returns>Returns <c>true</c>, if is skew-symmetric, <c>false</c> otherwise.</returns>
        /// <param name="A">The source matrix.</param>
        /// <param name="N">The source matrix dimension.</param>
        /// <remarks>
        /// In mathematics, particularly in linear algebra, a skew-symmetric 
        /// (or antisymmetric or antimetric) matrix is a square matrix whose 
        /// transpose equals its negative. 
        /// 
        ///     That is: A[i][j] = -A[j][i]
        /// 
        /// </remarks>
        public bool IsSkewSymmetric(T[][] A, int N)
        {
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    if (!(A[i][j].Equals( this.Op.Minus( A[j][i] ) ) ) )
                        return false;
                }
            }

            return true;
        }


        /// <summary>
        /// Checks if a matrix is symmetric or Hermitian.
        /// </summary>
        /// <returns>Returns <c>true</c>, if is hermitian, <c>false</c> otherwise.</returns>
        /// <param name="A">The source matrix.</param>
        /// <param name="N">The source matrix dimension.</param>
        /// <remarks>
        /// In linear algebra, a symmetric matrix is a square matrix that is 
        /// equal to its transpose. 
        /// </remarks>
        public static bool IsSymmetric(double[][] A, int N, double tol = Utils.DBL_EPSILON)
        {
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    if (Math.Abs(A[i][j] - A[j][i]) > tol)
                        return false;
                }
            }

            return true;
        }


        /// <summary>
        /// Checks if a double matrix is skew-symmetric.
        /// </summary>
        /// <returns>Returns <c>true</c>, if is skew-symmetric, <c>false</c> otherwise.</returns>
        /// <param name="A">The source matrix.</param>
        /// <param name="N">The source matrix dimension.</param>
        /// <remarks>
        /// In linear algebra, a skew-symmetric matrix is a square matrix that 
        /// its transpose is equals its negative. 
        /// </remarks>
        public static bool IsSkewSymmetric(double[][] A, int N, double tol = Utils.DBL_EPSILON)
        {
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    if (Math.Abs(-A[i][j] - A[j][i]) > tol)
                        return false;
                }
            }

            return true;
        }



        /// <summary>
        /// Checks if a complex matrix is hermitian.
        /// </summary>
        /// <returns>Returns <c>true</c>, if is hermitian, <c>false</c> otherwise.</returns>
        /// <param name="A">The source matrix.</param>
        /// <param name="N">The source matrix dimension.</param>
        /// <remarks>
        /// In mathematics, a Hermitian matrix (or self-adjoint matrix) is a 
        /// complex square matrix that is equal to its own conjugate 
        /// transpose—that is, the element in the i-th row and j-th column is 
        /// equal to the complex conjugate of the element in the j-th row and i-th column, for all indices i and j.
        ///                        _______
        ///     That is: A[i][j] = A[j][i]
        /// 
        /// </remarks>
        public static bool IsHermitian(cpx[][] A, int N)
        {
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    if (A[i][j] != cpx.Conjugate(A[j][i]) )
                        return false;
                }
            }

            return true;
        }


        /// <summary>
        /// Checks if a complex matrix is skew-Hermitian.
        /// </summary>
        /// <returns>Returns <c>true</c>, if is skew-Hermitian, <c>false</c> otherwise.</returns>
        /// <param name="A">The source matrix.</param>
        /// <param name="N">The source matrix dimension.</param>
        /// <remarks>
        /// In linear algebra, a square matrix with complex entries is said to 
        /// be skew-Hermitian or anti-Hermitian if its conjugate transpose is 
        /// the negative of the original matrix.
        ///                         _______
        ///     That is: A[i][j] = -A[j][i]
        /// 
        /// Skew-Hermitian matrices can be understood as the complex versions 
        /// of real skew-symmetric matrices, or as the matrix analogue of the 
        /// purely imaginary numbers.
        /// 
        /// </remarks>
        public static bool IsSkewHermitian(cpx[][] A, int N)
        {
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    if (A[i][j] != -cpx.Conjugate(A[j][i]))
                        return false;
                }
            }

            return true;
        }


        /// <summary>
        /// Checks if a square matrix is a Hankel matrix.
        /// </summary>
        /// <returns>Returns <c>true</c>, if it's a Hankel matrix, <c>false</c> otherwise.</returns>
        /// <param name="A">The source matrix.</param>
        /// <param name="N">The source matrix dimension.</param>
        /// <remarks>
        /// 
        /// In linear algebra, a Hankel matrix (or catalecticant matrix), named 
        /// after Hermann Hankel, is a square matrix in which each ascending 
        /// skew-diagonal from left to right is constant, e.g.: 
        /// 
        /// 
        ///     |a b c d e|
        ///     |b c d e f|
        ///     |c d e f g|
        ///     |d e f g h|
        ///     |e f g h i|
        /// 
        ///     That is: A[i][j] = A[i + k][j − k] for all k = 0, ..., j−i
        ///
        /// Notes:  The Hankel matrix is a symmetric matrix.
        ///         The Hilbert matrix is an example of a Hankel matrix.
        ///                         
        /// </remarks>
        public static bool IsHankelMatrix(T[][] A, int N)
        {
            // for each row
            for (int i = 0; i < N; i++)
            {
                // for each column
                for (int j = 0; j < N; j++)
                {

                    // checking if i + j is less
                    // than n
                    if (i + j < N)
                    {
                        // checking if the element
                        // is equal to the
                        // corresponding diagonal
                        // constant
                        if (!(A[i][j].Equals(A[i + j][0])))
                            return false;
                    }
                    else
                    {
                        // checking if the element
                        // is equal to the
                        // corresponding diagonal
                        // constant
                        if (!(A[i][j].Equals(
                           A[i + j - N + 1][N - 1] ) ) )
                            return false;
                    }
                }
            }

            return true;
        }


        ///// <summary>
        ///// Checks if an integer square matrix is a Hankel matrix.
        ///// </summary>
        ///// <returns>Returns <c>true</c>, if it's a Hankel matrix, <c>false</c> otherwise.</returns>
        ///// <param name="A">The source matrix.</param>
        ///// <param name="N">The source matrix dimension.</param>
        ///// <remarks>
        ///// 
        ///// In linear algebra, a Hankel matrix (or catalecticant matrix), named 
        ///// after Hermann Hankel, is a square matrix in which each ascending 
        ///// skew-diagonal from left to right is constant, e.g.: 
        ///// 
        ///// 
        /////     |a b c d e|
        /////     |b c d e f|
        /////     |c d e f g|
        /////     |d e f g h|
        /////     |e f g h i|
        ///// 
        /////     That is: A[i][j] = A[i + k][j − k] for all k = 0, ..., j−i
        /////
        ///// Notes:  The Hankel matrix is a symmetric matrix.
        /////         The Hilbert matrix is an example of a Hankel matrix.
        /////                         
        ///// </remarks>
        //public static bool IsHankelMatrix(int[][] A, int N)
        //{
        //    // for each row
        //    for (int i = 0; i < N; i++)
        //    {
        //        // for each column
        //        for (int j = 0; j < N; j++)
        //        {

        //            // checking if i + j is less
        //            // than n
        //            if (i + j < N)
        //            {
        //                // checking if the element
        //                // is equal to the
        //                // corresponding diagonal
        //                // constant
        //                if (A[i][j] != A[i + j][0])
        //                    return false;
        //            }
        //            else
        //            {
        //                // checking if the element
        //                // is equal to the
        //                // corresponding diagonal
        //                // constant
        //                if (A[i][j] !=
        //                   A[i + j - N + 1][N - 1])
        //                    return false;
        //            }
        //        }
        //    }

        //    return true;
        //}


        /// <summary>
        /// Checks if a double square matrix is a Hankel matrix.
        /// </summary>
        /// <returns>Returns <c>true</c>, if it's a Hankel matrix, <c>false</c> otherwise.</returns>
        /// <param name="A">The source matrix.</param>
        /// <param name="N">The source matrix dimension.</param>
        /// <param name="tol">Tolerance for floating point values comparison (default is C DBL_EPSILON).</param>
        /// <remarks>
        /// 
        /// In linear algebra, a Hankel matrix (or catalecticant matrix), named 
        /// after Hermann Hankel, is a square matrix in which each ascending 
        /// skew-diagonal from left to right is constant, e.g.: 
        /// 
        /// 
        ///     |a b c d e|
        ///     |b c d e f|
        ///     |c d e f g|
        ///     |d e f g h|
        ///     |e f g h i|
        /// 
        ///     That is: A[i][j] = A[i + k][j − k] for all k = 0, ..., j−i
        ///
        /// Notes:  The Hankel matrix is a symmetric matrix.
        ///         The Hilbert matrix is an example of a Hankel matrix.
        ///                         
        /// </remarks>
        public static bool IsHankelMatrix(double[][] A, int N, double tol = Utils.DBL_EPSILON)
        {
            // for each row
            for (int i = 0; i < N; i++)
            {
                // for each column
                for (int j = 0; j < N; j++)
                {

                    // checking if i + j is less
                    // than n
                    if (i + j < N)
                    {
                        // checking if the element
                        // is equal to the
                        // corresponding diagonal
                        // constant
                        if (Math.Abs( A[i][j] - A[i+j][0] ) > tol)  // if (A[i][j] != A[i + j][0])
                            return false;
                    }
                    else
                    {
                        // checking if the element
                        // is equal to the
                        // corresponding diagonal
                        // constant
                        if (Math.Abs(A[i][j] - A[i+j-N+1][N-1]) > tol) // if (A[i][j] != A[i + j - N + 1][N - 1])
                            return false;
                    }
                }
            }

            return true;
        }


        ///// <summary>
        ///// Checks if a complex square matrix is a Hankel matrix.
        ///// </summary>
        ///// <returns>Returns <c>true</c>, if it's a Hankel matrix, <c>false</c> otherwise.</returns>
        ///// <param name="A">The source matrix.</param>
        ///// <param name="N">The source matrix dimension.</param>
        ///// <remarks>
        ///// 
        ///// In linear algebra, a Hankel matrix (or catalecticant matrix), named 
        ///// after Hermann Hankel, is a square matrix in which each ascending 
        ///// skew-diagonal from left to right is constant, e.g.: 
        ///// 
        ///// 
        /////     |a b c d e|
        /////     |b c d e f|
        /////     |c d e f g|
        /////     |d e f g h|
        /////     |e f g h i|
        ///// 
        /////     That is: A[i][j] = A[i + k][j − k] for all k = 0, ..., j−i
        /////
        ///// Notes:  The Hankel matrix is a symmetric matrix.
        /////         The Hilbert matrix is an example of a Hankel matrix.
        /////                         
        ///// </remarks>
        //public static bool IsHankelMatrix(cpx[][] A, int N)
        //{
        //    // for each row
        //    for (int i = 0; i < N; i++)
        //    {
        //        // for each column
        //        for (int j = 0; j < N; j++)
        //        {

        //            // checking if i + j is less
        //            // than n
        //            if (i + j < N)
        //            {
        //                // checking if the element
        //                // is equal to the
        //                // corresponding diagonal
        //                // constant
        //                if (A[i][j] != A[i + j][0])
        //                    return false;
        //            }
        //            else
        //            {
        //                // checking if the element
        //                // is equal to the
        //                // corresponding diagonal
        //                // constant
        //                if (A[i][j] != A[i+j-N+1][N-1])
        //                    return false;
        //            }
        //        }
        //    }

        //    return true;
        //}


        #endregion Symmetric, skew-Symmetric and Hermitian and skew-Hermitian check



        /// <summary>
        /// Creates a Hilbert matrix of size <paramref name="N"/>.
        /// </summary>
        /// <returns>Returns a Hilbert matrix.</returns>
        /// <param name="N">The square matrix size.</param>
        /// <remarks>
        /// 
        /// In linear algebra, a Hilbert matrix, introduced by Hilbert(1894), is 
        /// a square matrix with entries being the unit fractions.
        /// 
        /// 
        ///     |1   1/2 1/3 1/4 1/5|
        ///     |1/2 1/3 1/4 1/5 1/6|
        ///     |1/3 1/4 1/5 1/6 1/7|
        ///     |1/4 1/5 1/6 1/7 1/8|
        ///     |1/5 1/6 1/7 1/8 1/9|
        /// 
        ///     That is: A[i][j] = 1 / ( (i+1)+(j+1)-1) )
        ///
        /// Notes:  The Hilbert matrix is symmetric and positive definite. 
        ///         The Hilbert matrix is also totally positive (meaning that the 
        ///         determinant of every submatrix is positive).
        ///
        ///         The Hilbert matrix is an example of a Hankel matrix.It is also a specific example of a Cauchy matrix. 
        /// 
        /// Appliations:
        /// 
        ///     The method of moments applied to polynomial distributions results
        ///     in a Hankel matrix, which in the special case of approximating a 
        ///     probability distribution on the interval [0, 1] results in a 
        ///     Hilbert matrix. This matrix needs to be inverted to obtain the 
        ///     weight parameters of the polynomial distribution approximation.
        ///                         
        /// </remarks>
        public static double[][] CreateHilbertMatrix(int N)
        {
            double[][] FResult = Matrix<double>.CreateJaggedArray(N, N);

            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    // using the formula
                    // to generate
                    // hilbert matrix
                    FResult[i][j] = (double)1.0 /
                            ((i + 1) + (j + 1) -
                            (double)1.0);
                }
            }

            return FResult;
        }


        /// <summary>
        /// Creates a Cauchy matrix of size m x n.
        /// </summary>
        /// <returns>Returns a Cauchy matrix .</returns>
        /// <param name="m">Number of matrix rows.</param>
        /// <param name="n">Number of matrix columns.</param>
        /// <param name="X">First list of values that are injective sequences (they contain distinct elements). Must have at least <paramref name="m"/> elements.</param>
        /// <param name="Y">Second list of values that are injective sequences (they contain distinct elements). Must have at least <paramref name="n"/> elements.</param>
        /// <remarks>
        /// 
        /// In mathematics, a Cauchy matrix, named after Augustin-Louis Cauchy, 
        /// is an m×n matrix with elements a[i][j] in the form
        /// 
        /// a[i][j] = 1 / (X[i] - Y[j]);   X[i]-Y[j] != 0  
        /// 
        /// where X[i] and Y[j] are elements of a field F, and (X[i]) and (Y[j]) are injective sequences (they contain distinct elements).
        ///
        /// The Hilbert matrix is a special case of the Cauchy matrix, where
        ///    X[i] − Y[j] = i + j − 1. 
        ///
        /// Every submatrix of a Cauchy matrix is itself a Cauchy matrix.
        /// 
        /// 
        /// </remarks>
        public static double[][] CreateCauchyMatrix(int m, int n, int[] X, int[] Y)
        {
            const string C_MSG_INVALID_ARGS = "Invalid argument values for Cauchy matrix";

            if (m > X.Length)
                throw new Exception($"Argument '{nameof(X)}' must have at least {m} elements.");

            if (n > Y.Length)
                throw new Exception($"Argument '{nameof(Y)}' must have at least {n} elements.");
                
            double[][] FResult = Matrix<double>.CreateJaggedArray(m, n);

            // condition: (X[i] - Y[j]) != 0

            for (int i = 0; i < m; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    if ((X[i] - Y[j]) != 0)
                        // using the formula
                        // to generate
                        // Cauchy matrix
                        FResult[i][j] = (double)1.0 /
                            (X[i] - Y[j]);
                    else
                        throw new Exception(C_MSG_INVALID_ARGS + ". Division by zero.");
                }
            }

            return FResult;
        }


        /// <summary>
        /// Checks if a matrix is a Toeplitz matrix.
        /// </summary>
        /// <returns>Returns <c>true</c>, if it's a Toeplitz matrix, <c>false</c> otherwise.</returns>
        /// <param name="A">The source matrix.</param>
        /// <remarks>
        /// 
        /// In linear algebra, a Toeplitz matrix or diagonal-constant matrix, 
        /// named after Otto Toeplitz, is a matrix in which each descending 
        /// diagonal from left to right is constant. For instance, the following 
        /// matrix is a Toeplitz matrix: 
        /// 
        /// 
        ///     |a b c d e|
        ///     |f a b c d|
        ///     |g f a b c|
        ///     |h g f a b|
        ///     |i h g f a|
        /// 
        ///     That is: A[i][j] = A[i + 1][j + 1] 
        ///
        /// Notes:  A Toeplitz matrix is not necessarily square. 
        /// 
        ///         Toeplitz matrices are also closely connected with Fourier series
        ///         , because the multiplication operator by a trigonometric polynomial
        ///         , compressed to a finite-dimensional space, can be represented 
        ///         by such a matrix. Similarly, one can represent linear 
        ///         convolution as multiplication by a Toeplitz matrix.
        ///        
        /// Time Complexity: O(mn), where m is number of rows and n is number of columns.
        /// Space Complexity: O(m+n), because at worst case, if a matrix is 
        ///                   Toeplitz, we have store exactly(m+n-1) key, value 
        ///                   pairs. (In first row we have n distinct keys and then for 
        ///                   next each m-1 rows, we keep adding one unique key to the map.
        /// 
        /// Credits to:  Aarti_Rathi and Aditya Goel                 
        /// </remarks>
        public static bool IsToeplitz(T[][] A)
        {
            // row = number of rows
            // col = number of columns
            int row = A.Length;
            int col = A[0].Length;

            // HashMap to store key,value pairs
            IDictionary<int, T> map
                = new Dictionary<int, T>();

            for (int i = 0; i < row; i++)
            {
                for (int j = 0; j < col; j++)
                {
                    int key = i - j;
                    // if key value exists in the hashmap,
                    if (map.ContainsKey(key))
                    {
                        // we check whether the current value
                        // stored in this key matches to element
                        // at current index or not. If not,
                        // return false
                        if (!EqualityComparer<T>.Default.Equals(map[key], A[i][j]))//;)  (map[key] != A[i][j])
                            return false;
                    }
                    // else we put key,value pair in hashmap
                    else
                    {
                        map.Add(key, A[i][j]);
                    }
                }
            }

            return true;
        }


        /// <summary>
        /// Adjusts matrix cell to zero if number is less than given <paramref name="tolerance"/>.
        /// </summary>
        /// <returns>Return the adjusted matrix.</returns>
        /// <param name="A">The source matrix.</param>
        /// <param name="tolerance">The tolerance value.</param>
        public static double[][] AdjustToValue(double value, double[][] A, double tolerance = Utils.DBL_TOLERANCE)
        {
            int numR = A.GetLength(0);
            int numC = A[0].GetLength(0);

            for (int i = 0; i < numR; i++)
            {
                for (int j = 0; j < numC; j++)
                {
                    if (Math.Abs(value - A[i][j]) <= tolerance)
                        A[i][j] = value; 
                }
            }

            return A;
        }


        /// <summary>
        /// Adjusts matrix cell to zero if number is less than given <paramref name="tolerance"/>
        /// using parallelism.
        /// </summary>
        /// <returns>Return the adjusted matrix.</returns>
        /// <param name="A">The source matrix.</param>
        /// <param name="tolerance">The tolerance value.</param>
        public static double[][] AdjustToValuePar(double value, double[][] A, double tolerance = Utils.DBL_TOLERANCE)
        {
            int numR = A.GetLength(0);
            int numC = A[0].GetLength(0);
            if (numR > numC)
            {
                Parallel.For(0, numR, delegate (int i) 
                    {
                        for (int j = 0; j < numC; j++)
                        {
                            if (Math.Abs(value - A[i][j]) <= tolerance)
                                A[i][j] = value;
                        }
                    }
                );
            }
            else
                for (int i = 0; i < numR; i++)
                {
                    Parallel.For(0, numC, delegate (int j)
                        {
                            if (Math.Abs(value - A[i][j]) <= tolerance)
                                A[i][j] = value;
                        }
                    );
                }

            return A;
        }


        #region Random matrix



        ///// <summary>
        ///// Generates a matrix filled with random values.
        ///// </summary>
        ///// <returns>Returns a matrix filled with random values.</returns>
        ///// <param name="rows">Number of matrix rows.</param>
        ///// <param name="cols">Number of matrix cols.</param>
        ///// <param name="minVal">Minimum random value.</param>
        ///// <param name="maxVal">Maximum random value.</param>
        ///// <param name="seed">Seed value for random generator.</param>
        //public static int[][] MatrixRandom(int rows, int cols,
        //                int minVal, int maxVal, int seed)
        //{
        //    // return a matrix with random values
        //    Random ran = new Random(seed);
        //    int[][] result = Matrix<int>.CreateJaggedArray(rows, cols);

        //    for (int i = 0; i < rows; ++i)
        //        for (int j = 0; j < cols; ++j)
        //            result[i][j] = ran.Next(minVal, maxVal);

        //    return result;
        //}





        #endregion Random matrix




        /// <summary>
        /// Decomposes matrix <paramref name="A"/> using LU method.
        /// </summary>
        /// <returns>Returns the decomposed matrix and the permutation vector <paramref name="P1"/> if succeeded, otherwise <c>null</c> or exception..</returns>
        /// <param name="A">The original square matrix.</param>
        /// <param name="P1">The permutation vector as out parameter.</param>
        public T[][] LUPDecompose(T[][] A, out int[] P1)
        {
            T[][] FResult = null;
            try
            {
                LUSolver<T> solv = new LUSolver<T>(this.Op);
                FResult = solv.LUPDecompose(A, out P1);
            }
            catch (Exception ex)
            {
                throw ex;
            }

            return FResult;
        }


        /// <summary>
        /// Solves linear equations system for given coeficients matrix <paramref name="A"/> and second member terms <paramref name="b"/>.
        /// </summary>
        /// <returns>Returns solutions vector X..</returns>
        /// <param name="A">Source matrix with equations coeficients.</param>
        /// <param name="b">Equations second member terms vector.</param>
        /// <param name="overwrite">If <c>true</c> source matrix will be overwritten.</param>
        /// <typeparam name="T">The 1st type parameter.</typeparam>
        public T[] LUSolve(T[][] A, T[] b, bool overwrite)
        {
            T[] FResult = null;
            LUSolver<T> solv = new LUSolver<T>(this.Op);
            FResult = solv.Solve(A, b, overwrite);
            return FResult;
        }


        /// <summary>
        /// Solves linear equations system for given coeficients matrix <paramref name="A"/> and second member terms <paramref name="b"/>
        /// using parallelization.
        /// </summary>
        /// <returns>Returns solutions vector X..</returns>
        /// <param name="A">Source matrix with equations coeficients.</param>
        /// <param name="b">Equations second member terms vector.</param>
        /// <param name="overwrite">If <c>true</c> source matrix will be overwritten.</param>
        /// <typeparam name="T">The 1st type parameter.</typeparam>
        public T[] LUSolvePar(T[][] A, T[] b, bool overwrite)
        {
            T[] FResult = null;
            LUSolver<T> solv = new LUSolver<T>(this.Op);
            FResult = solv.SolvePar(A, b, overwrite);
            return FResult;
        }


        /// <summary>
        /// Computes the lower triangular matrix from LU decomposition.
        /// </summary>
        /// <returns>Returns the lower triangular matrix.</returns>
        /// <param name="A">Source matrix.</param>
        /// <typeparam name="T">The 1st type parameter.</typeparam>
        public T[][] LULower(T[][] A)
        {
            int N = Matrix<T>.GetNumRows(A);
            LUSolver<T> solv = new LUSolver<T>(this.Op);
            solv.lUUpperLower(A, N, out T[][] lower, out T[][] upper);
            return lower;
        }


        /// <summary>
        /// Computes the lower triangular matrix from LU decomposition using
        /// parallelization.
        /// </summary>
        /// <returns>Returns the lower triangular matrix.</returns>
        /// <param name="A">Source matrix.</param>
        /// <typeparam name="T">The 1st type parameter.</typeparam>
        public T[][] LULowerPar(T[][] A)
        {
            int N = Matrix<T>.GetNumRows(A);
            LUSolver<T> solv = new LUSolver<T>(this.Op);
            solv.LUUpperLowerPar(A, N, out T[][] lower, out T[][] upper);
            return lower;
        }


        /// <summary>
        /// Computes the upper triangular matrix from LU decomposition.
        /// </summary>
        /// <returns>Returns the upper triangular matrix.</returns>
        /// <param name="A">Source matrix.</param>
        /// <typeparam name="T">The 1st type parameter.</typeparam>
        public T[][] LUUpper(T[][] A)
        {
            int N = Matrix<T>.GetNumRows(A);
            LUSolver<T> solv = new LUSolver<T>(this.Op);
            solv.lUUpperLower(A, N, out T[][] lower, out T[][] upper);
            return upper;
        }


        /// <summary>
        /// Computes the upper triangular matrix from LU decomposition
        /// using parallelization.
        /// </summary>
        /// <returns>Returns the upper triangular matrix.</returns>
        /// <param name="A">Source matrix.</param>
        /// <typeparam name="T">The 1st type parameter.</typeparam>
        public T[][] LUUpperPar(T[][] A)
        {
            int N = Matrix<T>.GetNumRows(A);
            LUSolver<T> solv = new LUSolver<T>(this.Op);
            solv.LUUpperLowerPar(A, N, out T[][] lower, out T[][] upper);
            return upper;
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
            T[][] FResult = (new InverseHelper<T>(this.Op)).Inverse(LUP, perm);
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
            T[][] FResult = (new InverseHelper<T>(this.Op)).InversePar(LUP, perm);
            return FResult;
        }


        /// <summary>
        /// Computes the inverse of the specified matrix <paramref name="A"/>.
        /// </summary>
        /// <returns>Returns the inverse matrix of <paramref name="A"/> if succeeded.</returns>
        /// <param name="A">The original matrix.</param>
        /// <typeparam name="T">The 1st type parameter.</typeparam>
        public T[][] Inverse(T[][] A) 
        {
            T[][] FResult = null;
            int n = A.Length;
            T[][] lum = (new LUSolver<T>(this.Op)).LUPDecompose(A, out int[] perm);

            if (lum == null)
                throw new Exception("Unable to compute inverse.");

            FResult = (new InverseHelper<T>(this.Op)).Inverse(lum, perm);
            return FResult;
        }


        /// <summary>
        /// Computes the inverse of the specified matrix <paramref name="A"/>
        /// using Parallelization.
        /// </summary>
        /// <returns>Returns the inverse matrix of <paramref name="A"/> if succeeded.</returns>
        /// <param name="A">The original matrix.</param>
        /// <typeparam name="T">The 1st type parameter.</typeparam>
        public T[][] InversePar(T[][] A)
        {
            T[][] FResult = null;
            int n = A.Length;
            T[][] lum = (new LUSolver<T>(this.Op)).LUPDecomposePar(A, out int[] perm);

            if (lum == null)
                throw new Exception("Unable to compute inverse.");

            FResult = (new InverseHelper<T>(this.Op)).InversePar(lum, perm);
            return FResult;
        }


        /// <summary>
        /// Computes the rank of the specified matrix <paramref name="A"/>.
        /// </summary>
        /// <returns>Returns the rank of matrix if succeeded, throws exception otherwise..</returns>
        /// <param name="A">A.</param>
        /// <typeparam name="T">The 1st type parameter.</typeparam>
        public int Rank(T[][] A)
        {
            int FResult = 0;
            FResult = (new RankHelper<T>(this.Op)).RankOfMatrix(A);
            return FResult;
        }


        /// <summary>
        /// Computes the rank of the specified matrix <paramref name="A"/> using parallelism.
        /// </summary>
        /// <returns>Returns the rank of matrix if succeeded, throws exception otherwise..</returns>
        /// <param name="A">A.</param>
        /// <typeparam name="T">The 1st type parameter.</typeparam>
        public int RankPar(T[][] A)
        {
            int FResult = 0;
            FResult = (new RankHelper<T>(this.Op)).RankOfMatrixPar(A);
            return FResult;
        }


        /// <summary>
        /// Creates a matrix from an array of rows (each element represents a 
        /// row in the matrix).
        /// </summary>
        /// <returns>Returns the matrix generated from vector of rows.</returns>
        /// <param name="v">The vector (array) of rows.</param>
        /// <typeparam name="T">The 1st type parameter.</typeparam>
        /// <example>
        /// 
        ///  {1, 2, 3} => { 
        ///                 {1}, 
        ///                 {2}, 
        ///                 {3} 
        ///               }
        /// 
        /// </example>
        public static T[][] MatrixFromRowsVector(T[] v)
        {
            int N = v.Length;
            T[][] FResult = new T[N][];

            for (int i = 0; i < N; i++)
                FResult[i] = new T[1] { v[i] };
                
            return FResult;
        }


        /// <summary>
        /// Creates a matrix from an array of columns (each element represents a 
        /// column in the matrix).
        /// </summary>
        /// <returns>Returns the matrix generated from vector of columns.</returns>
        /// <param name="v">The vector (array) of columns.</param>
        /// <typeparam name="T">The 1st type parameter.</typeparam>
        /// <example>
        /// 
        ///  {1, 2, 3} => { 
        ///                 {1, 2, 3} 
        ///               }
        /// 
        /// </example>
        public static T[][] MatrixFromColumnsVector(T[] v)
        {
            int N = v.Length;
            T[][] FResult = new T[1][];
            FResult[0] = new T[N];

            for (int i = 0; i < N; i++)
                FResult[0][i] = v[i];

            return FResult;
        }



        #region Matrix Frobenius inner product


        /// <summary>
        /// Computes the inner product of two matrices <paramref name="A"/> and <paramref name="B"/>.
        /// </summary>
        /// <returns>Returns the inner product of the two matrices.</returns>
        /// <param name="A">First matrix.</param>
        /// <param name="B">Second matrix.</param>
        /// <remarks>
        /// The two matrices must have the same dimension - same number of rows 
        /// and columns, but are not restricted to be square matrices.
        /// 
        ///               _______
        /// result = sum( A[i][j] * B[i][j] )
        /// 
        /// </remarks>
        public T FrobeniusInnerProd(T[][] A, T[][] B)
        {
            T FResult = this.Op.Zero;

            int numRA = Matrix<T>.GetNumRows(A);
            int numRB = Matrix<T>.GetNumRows(B);

            int numCA = Matrix<T>.GetNumCols(A);
            int numCB = Matrix<T>.GetNumCols(B);

            // A nd B must have same size
            if (!((numRA == numRB) && (numCA == numCB)))
                throw new ArgumentException(C_MSG_FROB_INNER_PROD);

            T sum = this.Op.Zero;
            for (int i = 0; i < numRA; i++)
            {
                for (int j = 0; j < numCA; j++)
                {
                    //       _______
                    //sum += A[i][j] * B[i][j];
                    sum = this.Op.Add( sum, this.Op.Mult( this.Op.Conjugate(A[i][j]), B[i][j] ) );
                }
            }

            return FResult = sum;
        }


        ///// <summary>
        ///// Computes the inner product of two matrices <paramref name="A"/> and <paramref name="B"/>.
        ///// </summary>
        ///// <returns>Returns the inner product of the two matrices.</returns>
        ///// <param name="A">First matrix.</param>
        ///// <param name="B">Second matrix.</param>
        ///// <remarks>
        ///// 
        ///// The two matrices must have the same dimension - same number of rows
        ///// and columns, but are not restricted to be square matrices.
        ///// 
        ///// result = sum( A[i][j] * B[i][j] )
        ///// 
        ///// </remarks>
        //public static double FrobeniusInnerProd(double[][] A, double[][] B)
        //{
        //    double FResult = 0;

        //    int numRA = Matrix<T>.GetNumRows<double>(A);
        //    int numRB = Matrix<T>.GetNumRows<double>(B);

        //    int numCA = Matrix<T>.GetNumCols<double>(A);
        //    int numCB = Matrix<T>.GetNumCols<double>(B);

        //    // A nd B must have same size
        //    if (!((numRA == numRB) && (numCA == numCB)))
        //        throw new ArgumentException(C_MSG_FROB_INNER_PROD);

        //    double sum = 0;
        //    for (int i = 0; i < numRA; i++)
        //    {
        //        for (int j = 0; j < numCA; j++)
        //        {
        //            sum += A[i][j] * B[i][j];
        //        }
        //    }

        //    return FResult;
        //}


        ///// <summary>
        ///// Computes the inner product of two matrices <paramref name="A"/> and <paramref name="B"/>.
        ///// </summary>
        ///// <returns>Returns the inner product of the two matrices.</returns>
        ///// <param name="A">First matrix.</param>
        ///// <param name="B">Second matrix.</param>
        ///// <remarks>
        ///// 
        ///// For complex number matrices,
        /////               _______
        ///// result = sum( A[i][j] * B[i][j] )
        ///// 
        ///// </remarks>
        //public static cpx FrobeniusInnerProd(cpx[][] A, cpx[][] B)
        //{
        //    cpx FResult = 0;

        //    int numRA = Matrix<T>.GetNumRows<cpx>(A);
        //    int numRB = Matrix<T>.GetNumRows<cpx>(B);

        //    int numCA = Matrix<T>.GetNumCols<cpx>(A);
        //    int numCB = Matrix<T>.GetNumCols<cpx>(B);

        //    // A nd B must have same size
        //    if (!((numRA == numRB) && (numCA == numCB)))
        //        throw new ArgumentException(C_MSG_FROB_INNER_PROD);

        //    cpx sum = 0;
        //    for (int i = 0; i < numRA; i++)
        //    {
        //        for (int j = 0; j < numCA; j++)
        //        {
        //            sum += cpx.Conjugate( A[i][j] ) * B[i][j];
        //        }
        //    }

        //    return FResult;
        //}


        #endregion Matrix Frobenius inner product



        #region Matrix Norm


        #region Matrix Frobenius Norm


        /// <summary>
        /// Computes the Frobenius norm of a given matrix <paramref name="A"/>.
        /// </summary>
        /// <returns>Returns the Frobenius norm of the given matrix.</returns>
        /// <param name="A">Source matrix.</param>
        /// <remarks>
        /// In mathematics, a matrix norm is a vector norm in a vector space 
        /// whose elements (vectors) are matrices (of given dimensions). 
        /// </remarks>
        public virtual double FrobeniusNorm(T[][] A)
        {
            double FResult = 0.0;

            int numRA = Matrix<T>.GetNumRows(A);
            int numCA = Matrix<T>.GetNumCols(A);

            T sumSq = this.Op.Zero; // 0;
            for (int i = 0; i < numRA; i++)
            {
                for (int j = 0; j < numCA; j++)
                {
                    //sumSq +=  A[i][j] * A[i][j];
                    sumSq = this.Op.Add(sumSq, this.Op.Mult( A[i][j], A[i][j] ) );

                }
            }

            FResult = this.Op.SqrtDbl(sumSq);
            return FResult;
        }


        ///// <summary>
        ///// Computes the Frobenius norm of a given matrix <paramref name="A"/>.
        ///// </summary>
        ///// <returns>Returns the Frobenius norm of the given matrix.</returns>
        ///// <param name="A">Source matrix.</param>
        //public static double FrobeniusNorm(double[][] A)
        //{
        //    double FResult = 0;

        //    int numRA = Matrix<T>.GetNumRows<double>(A);
        //    int numCA = Matrix<T>.GetNumCols<double>(A);

        //    double sumSq = 0;
        //    for (int i = 0; i < numRA; i++)
        //    {
        //        for (int j = 0; j < numCA; j++)
        //        {
        //            sumSq += A[i][j] * A[i][j];
        //        }
        //    }

        //    FResult = Math.Sqrt(sumSq);
        //    return FResult;
        //}


        


        #endregion Matrix Frobenius Norm



        #region Matrix L1 induced norm

        ///// <summary>
        ///// Calculates the induced L1 norm of given matrix <paramref name="A"/>.
        ///// </summary>
        ///// <param name="A">Source matrix.</param>
        ///// <returns>
        ///// Returns the maximum absolute column sum of the matrix.
        ///// </returns>
        //public static double L1Norm(int[][] A)
        //{
        //    var norm = 0d;
        //    int numR = Matrix<T>.GetNumRows<int>(A);
        //    int numC = Matrix<T>.GetNumCols<int>(A);

        //    for (var j = 0; j < numC; j++)
        //    {
        //        var s = 0d;
        //        for (var i = 0; i < numR; i++)
        //        {
        //            s += Math.Abs(A[i][j]);
        //        }
        //        norm = Math.Max(norm, s);
        //    }
        //    return norm;
        //}


        /// <summary>
        /// Calculates the induced L1 norm of given matrix <paramref name="A"/>.
        /// </summary>
        /// <param name="A">Source matrix.</param>
        /// <returns>
        /// Returns the maximum absolute column sum of the matrix.
        /// </returns>
        public double L1Norm(T[][] A)
        {
            T norm = this.Op.Zero;  // 0d;
            int numR = Matrix<T>.GetNumRows(A);
            int numC = Matrix<T>.GetNumCols(A);

            for (var j = 0; j < numC; j++)
            {
                T s = this.Op.Zero;// 0d;
                for (var i = 0; i < numR; i++)
                {
                    s = this.Op.Add( s, this.Op.Abs( A[i][j] ) );  // Abs is same as Magnitude foe complex numbers
                }

                norm = this.Op.IsGreaterThan(norm, s) ? norm : s;// ;  Math.Max(norm, s);
            }

            return this.Op.ToDouble( norm );
        }


        ///// <summary>
        ///// Calculates the induced L1 norm of given matrix <paramref name="A"/>.
        ///// </summary>
        ///// <param name="A">Source matrix.</param>
        ///// <returns>
        ///// Returns the maximum absolute column sum of the matrix.
        ///// </returns>
        //public static double L1Norm(cpx[][] A)
        //{
        //    var norm = 0d;
        //    int numR = Matrix<T>.GetNumRows<cpx>(A);
        //    int numC = Matrix<T>.GetNumCols<cpx>(A);

        //    for (var j = 0; j < numC; j++)
        //    {
        //        var s = 0d;
        //        for (var i = 0; i < numR; i++)
        //        {
        //            s += A[i][j].Magnitude;
        //        }
        //        norm = Math.Max(norm, s);
        //    }
        //    return norm;
        //}


        #endregion Matrix L1 induced norm



        #region Matrix induced infinity Norm
            

        /// <summary>
        /// Calculates the induced infinity norm of this matrix.
        /// </summary>
        /// <returns>Returns the maximum absolute row sum of the matrix.</returns>
        /// <param name="A">Source matrix.</param>
        public double InfinityNorm(T[][] A)
        {
            double FResult = int.MinValue;
            T partialRes = this.Op.Zero;
            int numR = Matrix<T>.GetNumRows(A);
            int numC = Matrix<T>.GetNumCols(A);

            for (var i = 0; i < numR; i++)
            {
                T sum = this.Op.Zero;   // 0.0d;
                for (var j = 0; j < numC; j++)
                {
                    //sum += Math.Abs( A[i][j] );
                    sum = this.Op.Add(sum, this.Op.Abs( A[i][j] ));
                }

                partialRes = this.Op.IsGreaterThan(partialRes, sum) ? partialRes : sum; //
                //(int)Math.Max(FResult, sum);
            }

            return FResult = this.Op.ToDouble(partialRes);
        }


        ///// <summary>
        ///// Calculates the induced infinity norm of this matrix.
        ///// </summary>
        ///// <returns>Returns the maximum absolute row sum of the matrix.</returns>
        ///// <param name="A">Source matrix.</param>
        //public static double IninityNorm(double[][] A)
        //{
        //    double FResult = double.MinValue;
        //    int numR = Matrix<T>.GetNumRows<double>(A);
        //    int numC = Matrix<T>.GetNumCols<double>(A);

        //    for (var i = 0; i < numR; i++)
        //    {
        //        var sum = 0.0d;
        //        for (var j = 0; j < numC; j++)
        //        {
        //            sum += Math.Abs(A[i][j]);
        //        }

        //        FResult = Math.Max(FResult, sum);
        //    }

        //    return FResult;
        //}


        ///// <summary>
        ///// Calculates the induced infinity norm of this matrix.
        ///// </summary>
        ///// <returns>Returns the maximum absolute row sum of the matrix.</returns>
        ///// <param name="A">Source matrix.</param>
        //public static double IninityNorm(cpx[][] A)
        //{
        //    double FResult = double.MinValue;
        //    int numR = Matrix<T>.GetNumRows<cpx>(A);
        //    int numC = Matrix<T>.GetNumCols<cpx>(A);

        //    for (var i = 0; i < numR; i++)
        //    {
        //        double sum = 0.0d;
        //        for (var j = 0; j < numC; j++)
        //        {
        //            sum += A[i][j].Magnitude;
        //        }

        //        FResult = Math.Max(FResult, sum);
        //    }

        //    return FResult;
        //}


        #endregion Matrix induced infinity Norm



        #region Matrix P-norms by rows and columns


        /// <summary>
        /// Calculates the p-norms of all row vectors of matrix <paramref name="A"/>.
        /// </summary>
        /// <returns>Returns the p-norms of all row vectors.</returns>
        /// <param name="A">Source matrix.</param>
        /// <param name="p">Norm type (1, 2, infinity).</param>
        /// <remarks>
        /// Typical values for p are 1.0 (L1, Manhattan norm), 2.0 
        /// (L2, Euclidean norm) and positive infinity (infinity norm).
        /// </remarks>
        public double[] P_NormByRows(T[][] A, double p)
        {
            double[] FResult = null;
            int numR = Matrix<T>.GetNumRows(A);
            int numC = Matrix<T>.GetNumCols(A);

            if (p <= 0.0)
            {
                throw new ArgumentOutOfRangeException(nameof(p), "Value must be positive.");
            }

            var ret = new double[numR];
            double s;

            if ( Utils.AlmostEqual(p, 2.0) )
            {
                for (int i = 0; i < numR; i++)
                {
                    s = 0;
                    for (int j = 0; j < numC; j++)
                    {
                        s += Math.Pow( this.Op.ToDouble( this.Op.Abs( A[i][j] ) ), 2 );
                    }

                    ret[i] = Math.Sqrt(s);
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
                        s += this.Op.ToDouble( this.Op.Abs( A[i][j] ) );
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
                        s = Math.Max( s, this.Op.ToDouble( this.Op.Abs( A[i][j] ) ) );
                    }

                    ret[i] = s;
                }
            }
            else
            {
                double invnorm = 1.0 / p;

                for (int i = 0; i < numR; i++)
                {
                    s = 0;
                    for (int j = 0; j < numC; j++)
                    {
                        s += Math.Pow( this.Op.ToDouble( this.Op.Abs( A[i][j] ) ), p );
                    }

                    ret[i] = Math.Pow( s, invnorm );
                }
            }

            FResult = ret;
            return FResult;
        }


        /// <summary>
        /// Calculates the p-norms of all column vectors of matrix <paramref name="A"/>.
        /// </summary>
        /// <returns>Returns the p-norms of all column vectors.</returns>
        /// <param name="A">Source matrix.</param>
        /// <param name="p">Norm type (1, 2, infinity).</param>
        /// <remarks>
        /// Typical values for p are 1.0 (L1, Manhattan norm), 2.0 
        /// (L2, Euclidean norm) and positive infinity (infinity norm).
        /// </remarks>
        public double[] P_NormByCols(T[][] A, double p)
        {
            double[] FResult = null;
            int numR = Matrix<T>.GetNumRows(A);
            int numC = Matrix<T>.GetNumCols(A);

            if (p <= 0.0)
            {
                throw new ArgumentOutOfRangeException(nameof(p), "Value must be positive.");
            }

            var ret = new double[numC];
            double s;

            if (Utils.AlmostEqual(p, 2.0))
            {
                for (int j = 0; j < numC; j++)
                {
                    s = 0;
                    for (int i = 0; i < numR; i++)
                    {
                        s += Math.Pow( this.Op.ToDouble( this.Op.Abs( A[i][j] ) ), 2);
                    }

                    ret[j] = Math.Sqrt(s);
                }
            }
            else if (Utils.AlmostEqual(p, 1.0))
            {
                // which is simply the sum of each column absolute values of the matrix; 
                for (int j = 0; j < numC; j++)
                {
                    s = 0;
                    for (int i = 0; i < numR; i++)
                    {
                        s += this.Op.ToDouble( this.Op.Abs( A[i][j] ) );
                    }

                    ret[j] = s;
                }
            }
            else if (double.IsPositiveInfinity(p))
            {
                // which is simply the maximum absolute column value of the matrix; 
                for (int j = 0; j < numC; j++)
                {
                    s = 0;
                    for (int i = 0; i < numR; i++)
                    {
                        s = Math.Max(s, this.Op.ToDouble( this.Op.Abs( A[i][j] ) ) );
                    }

                    ret[j] = s;
                }
            }
            else
            {
                double invnorm = 1.0 / p;

                for (int j = 0; j < numC; j++)
                {
                    s = 0;
                    for (int i = 0; i < numR; i++)
                    {
                        s += Math.Pow(this.Op.ToDouble( this.Op.Abs( A[i][j] ) ), p);
                    }

                    ret[j] = Math.Pow(s, invnorm);
                }
            }

            FResult = ret;
            return FResult;
        }


        ///// <summary>
        ///// Calculates the p-norms of all row vectors of matrix <paramref name="A"/>.
        ///// </summary>
        ///// <returns>Returns the p-norms of all row vectors.</returns>
        ///// <param name="A">Source matrix.</param>
        ///// <param name="p">Norm type (1, 2, infinity).</param>
        ///// <remarks>
        ///// Typical values for p are 1.0 (L1, Manhattan norm), 2.0 
        ///// (L2, Euclidean norm) and positive infinity (infinity norm).
        ///// </remarks>
        //public static double[] P_NormByRows(double[][] A, double p)
        //{
        //    double[] FResult = null;
        //    int numR = Matrix<T>.GetNumRows<double>(A);
        //    int numC = Matrix<T>.GetNumCols<double>(A);

        //    if (p <= 0.0)
        //    {
        //        throw new ArgumentOutOfRangeException(nameof(p), "Value must be positive.");
        //    }

        //    var ret = new double[numR];
        //    double s;

        //    if (Utils.AlmostEqual(p, 2.0))
        //    {
        //        for (int i = 0; i < numR; i++)
        //        {
        //            s = 0;
        //            for (int j = 0; j < numC; j++)
        //            {
        //                s += Math.Pow(Math.Abs(A[i][j]), 2);
        //            }

        //            ret[i] = Math.Sqrt(s);
        //        }
        //    }
        //    else if (Utils.AlmostEqual(p, 1.0))
        //    {
        //        // which is simply the sum of each row absolute values of the matrix; 
        //        for (int i = 0; i < numR; i++)
        //        {
        //            s = 0;
        //            for (int j = 0; j < numC; j++)
        //            {
        //                s += Math.Abs(A[i][j]);
        //            }

        //            ret[i] = s;
        //        }
        //    }
        //    else if (double.IsPositiveInfinity(p))
        //    {
        //        // which is simply the maximum absolute row value of the matrix; 
        //        for (int i = 0; i < numR; i++)
        //        {
        //            s = 0;
        //            for (int j = 0; j < numC; j++)
        //            {
        //                s = Math.Max(s, Math.Abs(A[i][j]));
        //            }

        //            ret[i] = s;
        //        }
        //    }
        //    else
        //    {
        //        double invnorm = 1.0 / p;

        //        for (int i = 0; i < numR; i++)
        //        {
        //            s = 0;
        //            for (int j = 0; j < numC; j++)
        //            {
        //                s += Math.Pow(Math.Abs(A[i][j]), p);
        //            }

        //            ret[i] = Math.Pow(s, invnorm);
        //        }
        //    }

        //    FResult = ret;
        //    return FResult;
        //}


        ///// <summary>
        ///// Calculates the p-norms of all column vectors of matrix <paramref name="A"/>.
        ///// </summary>
        ///// <returns>Returns the p-norms of all column vectors.</returns>
        ///// <param name="A">Source matrix.</param>
        ///// <param name="p">Norm type (1, 2, infinity).</param>
        ///// <remarks>
        ///// Typical values for p are 1.0 (L1, Manhattan norm), 2.0 
        ///// (L2, Euclidean norm) and positive infinity (infinity norm).
        ///// </remarks>
        //public static double[] P_NormByCols(double[][] A, double p)
        //{
        //    double[] FResult = null;
        //    int numR = Matrix<T>.GetNumRows<double>(A);
        //    int numC = Matrix<T>.GetNumCols<double>(A);

        //    if (p <= 0.0)
        //    {
        //        throw new ArgumentOutOfRangeException(nameof(p), "Value must be positive.");
        //    }

        //    var ret = new double[numC];
        //    double s;

        //    if (Utils.AlmostEqual(p, 2.0))
        //    {
        //        for (int j = 0; j < numC; j++)
        //        {
        //            s = 0;
        //            for (int i = 0; i < numR; i++)
        //            {
        //                s += Math.Pow(Math.Abs(A[i][j]), 2);
        //            }

        //            ret[j] = Math.Sqrt(s);
        //        }
        //    }
        //    else if (Utils.AlmostEqual(p, 1.0))
        //    {
        //        // which is simply the sum of each column absolute values of the matrix; 
        //        for (int j = 0; j < numC; j++)
        //        {
        //            s = 0;
        //            for (int i = 0; i < numR; i++)
        //            {
        //                s += Math.Abs(A[i][j]);
        //            }

        //            ret[j] = s;
        //        }
        //    }
        //    else if (double.IsPositiveInfinity(p))
        //    {
        //        // which is simply the maximum absolute column value of the matrix; 
        //        for (int j = 0; j < numC; j++)
        //        {
        //            s = 0;
        //            for (int i = 0; i < numR; i++)
        //            {
        //                s = Math.Max(s, Math.Abs(A[i][j]));
        //            }

        //            ret[j] = s;
        //        }
        //    }
        //    else
        //    {
        //        double invnorm = 1.0 / p;

        //        for (int j = 0; j < numC; j++)
        //        {
        //            s = 0;
        //            for (int i = 0; i < numR; i++)
        //            {
        //                s += Math.Pow(Math.Abs(A[i][j]), p);
        //            }

        //            ret[j] = Math.Pow(s, invnorm);
        //        }
        //    }

        //    FResult = ret;
        //    return FResult;
        //}


        ///// <summary>
        ///// Calculates the p-norms of all row vectors of matrix <paramref name="A"/>.
        ///// </summary>
        ///// <returns>Returns the p-norms of all row vectors.</returns>
        ///// <param name="A">Source matrix.</param>
        ///// <param name="p">Norm type (1, 2, infinity).</param>
        ///// <remarks>
        ///// Typical values for p are 1.0 (L1, Manhattan norm), 2.0 
        ///// (L2, Euclidean norm) and positive infinity (infinity norm).
        ///// </remarks>
        //public static double[] P_NormByRows(cpx[][] A, double p)
        //{
        //    double[] FResult = null;
        //    int numR = Matrix<T>.GetNumRows<cpx>(A);
        //    int numC = Matrix<T>.GetNumCols<cpx>(A);

        //    if (p <= 0.0)
        //    {
        //        throw new ArgumentOutOfRangeException(nameof(p), "Value must be positive.");
        //    }

        //    var ret = new double[numR];
        //    double s;

        //    if (Utils.AlmostEqual(p, 2.0))
        //    {
        //        for (int i = 0; i < numR; i++)
        //        {
        //            s = 0;
        //            for (int j = 0; j < numC; j++)
        //            {
        //                s += Math.Pow(A[i][j].Magnitude, 2);
        //            }

        //            ret[i] = Math.Sqrt(s);
        //        }
        //    }
        //    else if (Utils.AlmostEqual(p, 1.0))
        //    {
        //        // which is simply the sum of each row absolute values of the matrix; 
        //        for (int i = 0; i < numR; i++)
        //        {
        //            s = 0;
        //            for (int j = 0; j < numC; j++)
        //            {
        //                s += A[i][j].Magnitude;
        //            }

        //            ret[i] = s;
        //        }
        //    }
        //    else if (double.IsPositiveInfinity(p))
        //    {
        //        // which is simply the maximum absolute row value of the matrix; 
        //        for (int i = 0; i < numR; i++)
        //        {
        //            s = 0;
        //            for (int j = 0; j < numC; j++)
        //            {
        //                s = Math.Max(s, A[i][j].Magnitude);
        //            }

        //            ret[i] = s;
        //        }
        //    }
        //    else
        //    {
        //        double invnorm = 1.0 / p;

        //        for (int i = 0; i < numR; i++)
        //        {
        //            s = 0;
        //            for (int j = 0; j < numC; j++)
        //            {
        //                s += Math.Pow(A[i][j].Magnitude, p);
        //            }

        //            ret[i] = Math.Pow(s, invnorm);
        //        }
        //    }

        //    FResult = ret;
        //    return FResult;
        //}


        ///// <summary>
        ///// Calculates the p-norms of all column vectors of matrix <paramref name="A"/>.
        ///// </summary>
        ///// <returns>Returns the p-norms of all column vectors.</returns>
        ///// <param name="A">Source matrix.</param>
        ///// <param name="p">Norm type (1, 2, infinity).</param>
        ///// <remarks>
        ///// Typical values for p are 1.0 (L1, Manhattan norm), 2.0 
        ///// (L2, Euclidean norm) and positive infinity (infinity norm).
        ///// </remarks>
        //public static double[] P_NormByCols(cpx[][] A, double p)
        //{
        //    double[] FResult = null;
        //    int numR = Matrix<T>.GetNumRows<cpx>(A);
        //    int numC = Matrix<T>.GetNumCols<cpx>(A);

        //    if (p <= 0.0)
        //    {
        //        throw new ArgumentOutOfRangeException(nameof(p), "Value must be positive.");
        //    }

        //    var ret = new double[numC];
        //    double s;

        //    if (Utils.AlmostEqual(p, 2.0))
        //    {
        //        for (int j = 0; j < numC; j++)
        //        {
        //            s = 0;
        //            for (int i = 0; i < numR; i++)
        //            {
        //                s += Math.Pow(A[i][j].Magnitude, 2);
        //            }

        //            ret[j] = Math.Sqrt(s);
        //        }
        //    }
        //    else if (Utils.AlmostEqual(p, 1.0))
        //    {
        //        // which is simply the sum of each column absolute values of the matrix; 
        //        for (int j = 0; j < numC; j++)
        //        {
        //            s = 0;
        //            for (int i = 0; i < numR; i++)
        //            {
        //                s += A[i][j].Magnitude;
        //            }

        //            ret[j] = s;
        //        }
        //    }
        //    else if (double.IsPositiveInfinity(p))
        //    {
        //        // which is simply the maximum absolute column value of the matrix; 
        //        for (int j = 0; j < numC; j++)
        //        {
        //            s = 0;
        //            for (int i = 0; i < numR; i++)
        //            {
        //                s = Math.Max(s, A[i][j].Magnitude);
        //            }

        //            ret[j] = s;
        //        }
        //    }
        //    else
        //    {
        //        double invnorm = 1.0 / p;

        //        for (int j = 0; j < numC; j++)
        //        {
        //            s = 0;
        //            for (int i = 0; i < numR; i++)
        //            {
        //                s += Math.Pow(A[i][j].Magnitude, p);
        //            }

        //            ret[j] = Math.Pow(s, invnorm);
        //        }
        //    }

        //    FResult = ret;
        //    return FResult;
        //}


        #endregion Matrix P-norms by rows and columns



        //#region Row vectors p-norms



        //#endregion Row vectors p-norms




        #endregion Matrix Norm


        /// <summary>
        /// Computes the QR decomposition matrices for a given matrix <paramref name="A"/> using 
        /// Gram–Schmidt Algorithm.
        /// The Gram-Schmidt process is inherently numerically unstable. While the application of 
        /// the projections has an appealing geometric analogy to orthogonalization, 
        /// the orthogonalization itself is prone to numerical error. A significant advantage 
        /// however is the ease of implementation, which makes this a useful algorithm to use for 
        /// prototyping if a pre-built linear algebra library is unavailable.
        /// </summary>
        /// <returns>Returns zero and the QR decomposition matrices if succeeded.</returns>
        /// <param name="A">Input matrix A.</param>
        /// <param name="Q">Output matrix Q (orthogonal matrix).</param>
        /// <param name="R">Output matrix R (Upper/right triangular matrix).</param>
        public int QRDecomp(T[][] A, out T[][] Q, out T[][] R)
        {
            int FResult = -1;
            QRDecomp<T> qr = new QRDecomp<T>(this.Op);
            FResult = qr.MatDecompQR(A, out Q, out R);
            return FResult;
        }


        /// <summary>
        /// Computes the QR decomposition matrices for a given matrix <paramref name="A"/> using 
        /// Gram–Schmidt Algorithm using parallelism.
        /// The Gram-Schmidt process is inherently numerically unstable. While the application of 
        /// the projections has an appealing geometric analogy to orthogonalization, 
        /// the orthogonalization itself is prone to numerical error. A significant advantage 
        /// however is the ease of implementation, which makes this a useful algorithm to use for 
        /// prototyping if a pre-built linear algebra library is unavailable.
        /// </summary>
        /// <returns>Returns zero and the QR decomposition matrices if succeeded.</returns>
        /// <param name="A">Input matrix A.</param>
        /// <param name="Q">Output matrix Q (orthogonal matrix).</param>
        /// <param name="R">Output matrix R (Upper/right triangular matrix).</param>
        public int QRDecompPar(T[][] A, out T[][] Q, out T[][] R)
        {
            int FResult = -1;
            QRDecomp<T> qr = new QRDecomp<T>(this.Op);
            FResult = qr.MatDecompQRPar(A, out Q, out R);
            return FResult;
        }


        /// <summary>
        /// Constructs a square diagonal matrix from a vector with diagonal values.
        /// </summary>
        /// <returns>Returns a square diagonal matrix.</returns>
        /// <param name="a">The input vector with diagonal values.</param>
        /// <remarks>
        /// In linear algebra, a diagonal matrix is a matrix in which the entries outside the 
        /// main diagonal are all zero; the term usually refers to square matrices. 
        /// Elements of the main diagonal can either be zero or nonzero.
        /// </remarks>
        public static T[][] DiagonalFromVector(T[] a)
        {
            int size = a.GetLength(0);
            T[][] FResult = Matrix<T>.CreateJaggedArray(size, size);
            for (int i = 0; i < size; i++)
            {
                FResult[i][i] = a[i];
            }

            return FResult;
        }


        /// <summary>
        ///  From a given matrix, constructs a vector of its diagonal entries. .
        /// </summary>
        /// <returns>Returns a vector with matrix diagonal entries.</returns>
        /// <remarks>
        /// </remarks>
        public static T[] DiagonalToVector(T[][] D)
        {
            int size = Math.Min(Matrix<T>.GetNumRows(D), Matrix<T>.GetNumCols(D));
            T[] FResult = new T[size];
            for (int i = 0; i < size; i++)
            {
                FResult[i] = D[i][i];
            }

            return FResult;
        }


        #endregion Methods


    }


}
