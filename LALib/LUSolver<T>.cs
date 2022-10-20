/*********************************************************************************
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
 ********************************************************************************/

using System;
using System.Threading.Tasks;

#if SYSNUMERICS
#else
using ComplexN;
#endif

namespace LALib
{


    /// <summary>
    /// Implements a class to solve linear equations systems by LU matrix factorization using rectangular arrays.
    /// </summary>
    /// <remarks>
    /// A = L*U
    /// A -> original matrix
    /// L -> lower triangular matrix
    /// U -> upper triangular matrix
    /// 
    /// Not all square matrices have an LU decomposition, and it may be necessary 
    /// to permute the rows of a matrix before obtaining its LU factorization. 
    /// 
    /// </remarks>
    public class LUSolver<T> //where T : IComparable
    {
        private IBasicMathOperations<T> Op { get; set; }

        public LUSolver(IBasicMathOperations<T> _op)
        {
            this.Op = _op;
        }


        #region methods


        #region Public methods


        /// <summary>
        /// Solve X for the given LU decomposed matrix <paramref name="LUP"/> and the permutation vector <paramref name="perm"/>" and right side vector <paramref name="b"/> (LUP*X=b).
        /// </summary>
        /// <param name="LUP">LU decomposed matrix from equations coeficients matrix.</param>
        /// <param name="perm">The permutation vector obtained from LU decomposition process.</param>
        /// <param name="b">The right side vector for LUPX=b (second member terms).</param>
        /// <returns>Returns the X vector results from LUPX=b.</returns>
        public T[] Solve(T[][] LUP, int[] perm, T[] b)
        {
            T[] FResult = null;
            int numRows = Matrix<T>.GetNumRows(LUP);
            int numCols = Matrix<T>.GetNumCols(LUP);

            // note: must be a square matrix
            if (numRows != numCols)
                throw new ArgumentException($"Equations first member coeficients must be a square matrix.");

            int N = LUP.GetLength(0);  // number of rows
            FResult = new T[N];

            // solve X
            LUPSolve(LUP, perm, b, N, ref FResult);

            return FResult;
        }


        /// <summary>
        /// Solve X for the specified coeficients matrix <paramref name="a"/> and right side vector <paramref name="b"/> (AX=b).
        /// </summary>
        /// <param name="a">Linear equations first member coeficients matrix.</param>
        /// <param name="b">The right side vector for AX=b (second member terms).</param>
        /// <param name="overwrite">If <c>true</c> original matrix  <paramref name="a"/> will be changed.</param>
        /// <returns>Returns the X vector results from AX=b.</returns>
        public T[] Solve(T[][] a, T[] b, bool overwrite)
        {
            T[] FResult = null;
            int numRows = Matrix<T>.GetNumRows(a);
            int numCols = Matrix<T>.GetNumCols(a);

            // note: must be a square matrix
            if (numRows != numCols)
                throw new ArgumentException($"Equations first member coeficients '{nameof(a)}' must be a square matrix.");

            // create permutation vector
            int[] Pv = new int[numRows + 1];
            T[][] LEU = null;

            if (overwrite)
                LEU = a;
            else
                LEU = Matrix<T>.Clone(a);

            // compute LU decomposition (PA = LU)
            if (LUPDecomposeInPlace(ref LEU, numRows, this.Op.Zero, ref Pv))
            {
                // LEU is changed, it contains a copy of both matrices L-E and U as LEU=(L-E)+U such that P*A=L*U.
                FResult = this.Solve(LEU, Pv, b);
            }
            else
            {
                throw new ApplicationException($"Input matrix can not be decomposed by LU method.");
            }

            return FResult;
        }


        /// <summary>
        /// Solve X for the specified coeficients matrix <paramref name="a"/> and right side vector <paramref name="b"/> (AX=b)
        /// using parallelization.
        /// </summary>
        /// <param name="a">Linear equations first member coeficients matrix.</param>
        /// <param name="b">The right side vector for AX=b (second member terms).</param>
        /// <param name="overwrite">If <c>true</c> original matrix  <paramref name="a"/> will be changed.</param>
        /// <returns>Returns the X vector results from AX=b.</returns>
        public T[] SolvePar(T[][] a, T[] b, bool overwrite)
        {
            T[] FResult = null;
            int numRows = Matrix<T>.GetNumRows(a);
            int numCols = Matrix<T>.GetNumCols(a);

            // note: must be a square matrix
            if (numRows != numCols)
                throw new ArgumentException($"Equations first member coeficients '{nameof(a)}' must be a square matrix.");

            // create permutation vector
            int[] Pv = new int[numRows + 1];
            T[][] LEU = null;

            if (overwrite)
                LEU = a;
            else
                LEU = Matrix<T>.ClonePar(a);

            // compute LU decomposition (PA = LU)
            if (LUPDecomposeInPlacePar(ref LEU, numRows, this.Op.Zero, ref Pv))
            {
                // LEU is changed, it contains a copy of both matrices L-E and U as LEU=(L-E)+U such that P*A=L*U.
                FResult = this.Solve(LEU, Pv, b);
            }
            else
            {
                throw new ApplicationException($"Input matrix can not be decomposed by LU method.");
            }

            return FResult;
        }


        /// <summary>
        /// Decomposes matrix <paramref name="A"/> using LU algorithm.
        /// </summary>
        /// <returns>Returns the decomposed matrix and the permutation vector <paramref name="P1"/> if succeeded, otherwise <c>null</c> or exception..</returns>
        /// <param name="A">The original square matrix.</param>
        /// <param name="P1">The permutation vector as out parameter.</param>
        public T[][] LUPDecompose(T[][] A, out int[] P1)
        {
            T[][] FResult = null;
            try
            {
                FResult = Matrix<T>.Clone(A);  // create copy because A will be chenged
                int N = A.Length;
                P1 = new int[N + 1];

                if (LUPDecomposeInPlace(ref FResult, N, this.Op.Zero, ref P1))
                {
                }
            }
            catch (Exception ex)
            {
                throw ex;
            }

            return FResult;
        }


        /// <summary>
        /// Decomposes matrix <paramref name="A"/> using LU algorithm using parallelization.
        /// </summary>
        /// <returns>Returns the decomposed matrix and the permutation vector <paramref name="P1"/> if succeeded, otherwise <c>null</c> or exception..</returns>
        /// <param name="A">The original square matrix.</param>
        /// <param name="P1">The permutation vector as out parameter.</param>
        public T[][] LUPDecomposePar(T[][] A, out int[] P1)
        {
            T[][] FResult = null;
            try
            {
                FResult = Matrix<T>.ClonePar(A);  // create copy because A will be chenged
                int N = A.Length;
                P1 = new int[N + 1];

                if (LUPDecomposeInPlacePar(ref FResult, N, this.Op.Zero, ref P1))
                {
                }
            }
            catch (Exception ex)
            {
                throw ex;
            }

            return FResult;
        }


        /// <summary>
        /// Computes the upper and lower triangular matrices from LU decomposition.
        /// </summary>
        /// <param name="mat">The original square matrix.</param>
        /// <param name="N">The matrix size.</param>
        /// <param name="lower">The computed lower triangular matrix.</param>
        /// <param name="upper">The computed upper triangular matrix.</param>
        public void LUUpperLower(T[][] mat, int N, out T[][] lower, out T[][] upper)
        {
            lower = Matrix<T>.CreateJaggedArray(N, N);
            upper = Matrix<T>.CreateJaggedArray(N, N);

            //dynamic d1 = 0.0;
            //dynamic d2 = 0.0;
            T sum = this.Op.Zero;

            // Decomposing matrix into Upper and Lower
            // triangular matrix
            for (int i = 0; i < N; i++)
            {
                // Upper Triangular
                for (int k = i; k < N; k++)
                {
                    // Summation of L(i, j) * U(j, k)
                    sum = this.Op.Zero;
                    for (int j = 0; j < i; j++)
                    {
                        //d1 = lower[i][j];
                        //d2 = upper[j][k];

                        sum = this.Op.Add( sum, this.Op.Mult( lower[i][j], upper[j][k] ) );
                    }

                    // Evaluating U(i, k)
                    //d1 = mat[i][k];
                    upper[i][k] = this.Op.Sub( mat[i][k], sum );
                }

                // Lower Triangular
                for (int k = i; k < N; k++)
                {
                    if (i == k)
                        lower[i][i] = this.Op.One; // Diagonal as 1
                    else
                    {
                        // Summation of L(k, j) * U(j, i)
                        sum = this.Op.Zero; //sum = 0;
                        for (int j = 0; j < i; j++)
                        {
                            //d1 = lower[k][j];
                            //d2 = upper[j][i];
                            sum = this.Op.Add( sum, this.Op.Mult( lower[k][j], upper[j][i] ) );
                        }

                        // Evaluating L(k, i)
                        //d1 = (d1 = mat[k][i]) - sum;
                        //d2 = upper[i][i];
                        lower[k][i] = this.Op.Div( this.Op.Sub( mat[k][i], sum ), upper[i][i] );
                    }
                }
            }
        }


        /// <summary>
        /// Computes the upper and lower triangular matrices from LU decomposition
        /// using parallelization.
        /// </summary>
        /// <param name="mat">The original square matrix.</param>
        /// <param name="N">The matrix size.</param>
        /// <param name="lower">The computed lower triangular matrix.</param>
        /// <param name="upper">The computed upper triangular matrix.</param>
        public void LUUpperLowerPar(T[][] mat, int N, out T[][] lower, out T[][] upper)
        {
            lower = Matrix<T>.CreateJaggedArray(N, N);
            upper = Matrix<T>.CreateJaggedArray(N, N);

            //dynamic d1 = 0.0;
            //dynamic d2 = 0.0;
            T sum = this.Op.Zero;

            T[][] l = lower;
            T[][] u = upper;

            // Decomposing matrix into Upper and Lower
            // triangular matrix
            for (int i = 0; i < N; i++)
            {
                // Upper Triangular
                Parallel.For(i, N, delegate (int k)
                    {
                        // Summation of L(i, j) * U(j, k)
                        sum = this.Op.Zero;
                        for (int j = 0; j < i; j++)
                        {
                            //d1 = lower[i][j];
                            //d2 = upper[j][k];

                            sum = this.Op.Add(sum, this.Op.Mult(l[i][j], u[j][k]));
                        }

                        // Evaluating U(i, k)
                        //d1 = mat[i][k];
                        u[i][k] = this.Op.Sub(mat[i][k], sum);
                    }
                );

                // Lower Triangular
                Parallel.For(i, N, delegate (int k)
                    {
                        if (i == k)
                            l[i][i] = this.Op.One; // Diagonal as 1
                        else
                        {
                            // Summation of L(k, j) * U(j, i)
                            sum = this.Op.Zero; //sum = 0;
                            for (int j = 0; j < i; j++)
                            {
                                //d1 = lower[k][j];
                                //d2 = upper[j][i];
                                sum = this.Op.Add(sum, this.Op.Mult(l[k][j], u[j][i]));
                            }

                            // Evaluating L(k, i)
                            //d1 = (d1 = mat[k][i]) - sum;
                            //d2 = upper[i][i];
                            l[k][i] = this.Op.Div(this.Op.Sub(mat[k][i], sum), u[i][i]);
                        }
                    }
                );
            }
        }



        #region Determinant calculation


        /// <summary>
        /// Computes the matrix determinant from decomposed matrix <paramref name="A"/>.
        /// </summary>
        /// <returns>Returns the determinant from a decomposed matrix.</returns>
        /// <param name="LUP">The decomposed matrix in the form A = LU * P.</param>
        /// <remarks>
        /// 
        ///     Note: this function applies to int, double and complex matrices.
        /// 
        /// </remarks>
        public T DeterminantFromLU(T[][] LUP)
        {
            if (LUP == null)
                throw new Exception("Unable to compute matrix determinant.");

            T FResult = this.Op.One;
            T val = this.Op.Zero; // 0.0;

            for (int i = 0; i < LUP.Length; ++i)
            {
                val = LUP[i][i];
                FResult = this.Op.Mult(FResult, val);
            }

            return FResult;
        }


        /// <summary>
        /// Computes the matrix determinant from decomposed matrix <paramref name="A"/>.
        /// </summary>
        /// <returns>Returns the determinant from a decomposed matrix.</returns>
        /// <param name="A">The original matrix.</param>
        /// <remarks>
        ///     Determinant = Product of decomposed matrix diagonal values.
        /// 
        ///     Note: this function applies to integer, double and complex matrices.
        /// </remarks>
        public T Determinant(T[][] A)
        {
            T FResult = this.Op.Zero;// 0.0;

            T[][] LUP = this.LUPDecompose(A, out int[] perm);

            //T[][] LUP = (new Matrix<T>(this.Op)).LUPDecompose(A, out int[] perm);

            if (LUP == null)
                throw new Exception("Unable to compute matrix determinant.");

            FResult = this.DeterminantFromLU(LUP);
            return FResult;
        }


        #endregion Determinant calculation



        #endregion Public methods


        /// <summary>
        /// Computes LU decomposition for a given matrix <paramref name="A"/>.
        /// </summary>
        /// <returns>Returns <c>true</c>, if decomposition was succeeded, <c>false</c> otherwise.</returns>
        /// <param name="A">Original square matrix to decompose (will be modified).</param>
        /// <param name="N">Dimension of square matrix (nº rows).</param>
        /// <param name="Tol">Small tolerance number to detect failure when the matrix is near degenerate.</param>
        /// <param name="P1">Permutation vector of size N+1. Will contain column indexes where the 
        ///                  permutation matrix has "1".
        /// </param>
        /// <remarks>
        ///     Output ->
        ///       Matrix A is changed, it contains a copy of both matrices L-E and U as A=(L-E)+U such that P*A=L*U.
        ///       The permutation matrix is not stored as a matrix, but in an integer vector P of size N+1 
        ///       containing column indexes where the permutation matrix has "1". The last element P[N] = S + N,
        ///       where S is the number of row exchanges needed for determinant computation, det(P)=(-1)^S
        ///
        /// </remarks>
        private bool LUPDecomposeInPlace(ref T[][] A, int N, T Tol, ref int[] P1)
        {
            int i, j, k, imax;
            T maxA, absA;
            //dynamic maxA, absA;
            T[] ptr;    // matrix row (used for changing rows)

            if (P1.Length != (N + 1))   // alloc memory if necessary
                P1 = new int[N + 1];

            for (i = 0; i <= N; i++)
                P1[i] = i; //Unit permutation matrix, P[N] initialized with N

            for (i = 0; i < N; i++)
            {
                maxA = this.Op.Zero;// 0.0;
                imax = i;

                for (k = i; k < N; k++)
                    if ( this.Op.IsGreaterThan( ( absA = this.Op.Abs( A[k][i] ) ), maxA) )
                    {
                        maxA = absA;
                        imax = k;
                    }

                if (this.Op.IsLesserThan( maxA, Tol) ) return false; //failure, matrix is degenerate

                if (imax != i)
                {
                    //pivoting P
                    j = P1[i];
                    P1[i] = P1[imax];
                    P1[imax] = j;

                    //pivoting rows of A
                    ptr = A[i];
                    A[i] = A[imax];
                    A[imax] = ptr;

                    //counting pivots starting from N (for determinant)
                    P1[N]++;
                }

                for (j = i + 1; j < N; j++)
                {
                    T a = A[i][i];
                    A[j][i] = this.Op.Div(A[j][i], a);
                    a = A[j][i];

                    for (k = i + 1; k < N; k++)
                    {
                        T b = A[i][k];
                        A[j][k] = this.Op.Sub( A[j][k], this.Op.Mult( a, b ) );
                    }
                }
            }

            return true;  //decomposition done 
        }


        /// <summary>
        /// Computes LU decomposition for a given matrix <paramref name="A"/> using parallelization.
        /// </summary>
        /// <returns>Returns <c>true</c>, if decomposition was succeeded, <c>false</c> otherwise.</returns>
        /// <param name="A">Original square matrix to decompose (will be modified).</param>
        /// <param name="N">Dimension of square matrix (nº rows).</param>
        /// <param name="Tol">Small tolerance number to detect failure when the matrix is near degenerate.</param>
        /// <param name="P1">Permutation vector of size N+1. Will contain column indexes where the 
        ///                  permutation matrix has "1".
        /// </param>
        /// <remarks>
        ///     Output ->
        ///       Matrix A is changed, it contains a copy of both matrices L-E and U as A=(L-E)+U such that P*A=L*U.
        ///       The permutation matrix is not stored as a matrix, but in an integer vector P of size N+1 
        ///       containing column indexes where the permutation matrix has "1". The last element P[N] = S + N,
        ///       where S is the number of row exchanges needed for determinant computation, det(P)=(-1)^S
        ///
        /// </remarks>
        private bool LUPDecomposeInPlacePar(ref T[][] A, int N, T Tol, ref int[] P1)
        {
            int i, j, k, imax;
            T maxA, absA;
            //dynamic maxA, absA;
            T[] ptr;    // matrix row (used for changing rows)

            if (P1.Length != (N + 1))   // alloc memory if necessary
                P1 = new int[N + 1];

            for (i = 0; i <= N; i++)
                P1[i] = i; //Unit permutation matrix, P[N] initialized with N

            for (i = 0; i < N; i++)
            {
                maxA = this.Op.Zero;// 0.0;
                imax = i;

                for (k = i; k < N; k++)
                    if (this.Op.IsGreaterThan((absA = this.Op.Abs(A[k][i])), maxA))
                    {
                        maxA = absA;
                        imax = k;
                    }

                if (this.Op.IsLesserThan(maxA, Tol)) return false; //failure, matrix is degenerate

                if (imax != i)
                {
                    //pivoting P
                    j = P1[i];
                    P1[i] = P1[imax];
                    P1[imax] = j;

                    //pivoting rows of A
                    ptr = A[i];
                    A[i] = A[imax];
                    A[imax] = ptr;

                    //counting pivots starting from N (for determinant)
                    P1[N]++;
                }

                for (j = i + 1; j < N; j++)
                {
                    T a = A[i][i];
                    A[j][i] = this.Op.Div(A[j][i], a);
                    a = A[j][i];

                    T[][] A1 = A;
                    Parallel.For(i + 1, N, delegate (int k1)
                        {
                            T b = A1[i][k1];
                            A1[j][k1] = this.Op.Sub(A1[j][k1], this.Op.Mult(a, b));
                        }
                    );
                    //for (k = i + 1; k < N; k++)
                    //{

                    //}
                }
            }

            return true;  //decomposition done 
        }


        /// <summary>
        /// Solves a linear equations system using LU decomposition method.
        /// </summary>
        /// <param name="A">Output matrix filled in LUPDecompose() method <see cref="LUPDecomposeInPlace(ref T[][], int, double, ref int[])"/>.
        ///                 "A" contains a copy of both matrices L-E and U as A=(L-E)+U such that P*A=L*U.
        /// </param>
        /// <param name="P1">Output Permutation vector filled in LUPDecompose() method <see cref="LUPDecomposeInPlace(ref T[][], int, double, ref int[])"/>.</param>
        /// <param name="b">Vector filled with second member terms of linear equation system.</param>
        /// <param name="N">Matrix dimension of A (nº rows).</param>
        /// <param name="x">Returns the solution vector of A*x=b.</param>
        /// <remarks>
        ///     OUTPUT -> x - solution vector of A*x=b
        /// </remarks>
        private void LUPSolve(T[][] A, int[] P1, T[] b, int N, ref T[] x)
        {
            for (int i = 0; i < N; i++)
            {
                x[i] = b[P1[i]];

                T a1;
                T a2;
                for (int k = 0; k < i; k++)
                {
                    a1 = A[i][k];
                    a2 = x[k];
                    x[i] = this.Op.Sub( x[i], this.Op.Mult( a1, a2 ) );
                    //x[i] -= A[i][k] * x[k];
                }
            }

            for (int i = N - 1; i >= 0; i--)
            {
                //dynamic a1;
                //dynamic a2;
                for (int k = i + 1; k < N; k++)
                {
                    //a1 = A[i][k];
                    //a2 = x[k];
                    //x[i] -= a1 * a2;
                    x[i] = this.Op.Sub( x[i], this.Op.Mult( A[i][k], x[k] ) );
                }

                //a1 = x[i];
                //a2 = A[i][i];
                //x[i] = a1 / a2;
                x[i] = this.Op.Div( x[i], A[i][i] );

                //x[i] = x[i] / A[i][i];
            }
        }


        #endregion methods


    }

}
