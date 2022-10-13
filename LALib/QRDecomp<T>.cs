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
 * Credits and source code based on following sources:
 * 
 *  https://jamesmccaffrey.wordpress.com/2020/04/20/implementing-matrix-qr-decomposition-from-scratch-using-the-gram-schmidt-algorithm-with-c/
 *  https://en.wikipedia.org/wiki/QR_decomposition
 *
 ********************************************************************************/


using System;
using System.Threading.Tasks;

namespace LALib
{

    /// <summary>
    /// Implements a class to compute the QR decomposition matrices for a given matrix <paramref name="A"/> 
    /// using Gram–Schmidt Algorithm.
    /// The Gram-Schmidt process is inherently numerically unstable. While the application of 
    /// the projections has an appealing geometric analogy to orthogonalization, 
    /// the orthogonalization itself is prone to numerical error. A significant advantage 
    /// however is the ease of implementation, which makes this a useful algorithm to use for 
    /// prototyping if a pre-built linear algebra library is unavailable.
    /// </summary>
    /// <returns>Returns the QR decomposition matrices if succeeded.</returns>
    /// <param name="A">Input matrix A.</param>
    /// <param name="Q">Output matrix Q (orthogonal matrix).</param>
    /// <param name="R">Output matrix R (Upper/right triangular matrix).</param>
    /// <remarks>
    /// 
    /// This implementation should work well with double type values matrices.
    /// I do not know if this implemention is right for matrices of complex numbers 
    /// (sorry for my ignorance).
    /// 
    /// Do not use integer type matrices, because in the computations, decimal values will be 
    /// truncated to integers.
    /// 
    /// </remarks>
    public class QRDecomp<T>
    {

        private IBasicMathOperations<T> Op { get; set; }


        public QRDecomp(IBasicMathOperations<T> _op)
        {
            this.Op = _op;
        }


        /// <summary>
        /// Computes the norm of a vector (one dimensional array).
        /// </summary>
        /// <returns>Returns the norm of a given vector.</returns>
        /// <param name="vec">Source vector.</param>
        public T VecNorm(T[] vec)
        {
            T sum = this.Op.Zero;   // 0.0;
            for (int i = 0; i < vec.Length; ++i)
                sum = this.Op.Add(sum,  this.Op.Mult(vec[i], vec[i]) );

            return this.Op.Sqrt( sum );
        }


        /// <summary>
        /// Computes the dot product of two vectors (one dimensional arrays).
        /// </summary>
        /// <returns>Returns the dot product of two vectors.</returns>
        /// <param name="u">First vector.</param>
        /// <param name="v">Second vector.</param>
        private T VecDotProd(T[] u, T[] v)
        {
            T result = this.Op.Zero;
            for (int i = 0; i < u.Length; ++i)
                result = this.Op.Add(result, this.Op.Mult( this.Op.Conjugate( u[i] ), v[i] ) );
            return result;
        }


        /// <summary>
        /// Computes the projection of vector <paramref name="u"/> over vector <paramref name="a"/>.
        /// </summary>
        /// <returns>Returns the projection of vector <paramref name="u"/> over vector <paramref name="a"/>.</returns>
        /// <param name="u">First vector (1 dim array).</param>
        /// <param name="a">Second vector (1 dim array).</param>
        private T[] VecProjection(T[] u, T[] a)
        {
            // proj(u, a) = (inner(u,a) / inner(u, u)) * u
            // u cannot be all 0s
            int n = u.Length;
            T dotUA = VecDotProd(u, a);
            T dotUU = VecDotProd(u, u);
            T[] result = new T[n];

            for (int i = 0; i < n; ++i)
                result[i] = this.Op.Mult( this.Op.Div(dotUA, dotUU), u[i] );

            return result;
        }


        /// <summary>
        /// Computes the projection of vector <paramref name="u"/> over vector <paramref name="a"/>
        /// using parallelism.
        /// </summary>
        /// <returns>Returns the projection of vector <paramref name="u"/> over vector <paramref name="a"/>.</returns>
        /// <param name="u">First vector (1 dim array).</param>
        /// <param name="a">Second vector (1 dim array).</param>
        private T[] VecProjectionPar(T[] u, T[] a)
        {
            // proj(u, a) = (inner(u,a) / inner(u, u)) * u
            // u cannot be all 0s
            int n = u.Length;
            T dotUA = VecDotProd(u, a);
            T dotUU = VecDotProd(u, u);
            T[] result = new T[n];

            Parallel.For(0, n, delegate (int i)
               {
                   result[i] = this.Op.Mult(this.Op.Div(dotUA, dotUU), u[i]);
               }
            );
            return result;
        }


        /// <summary>
        /// Multiplies two matrices, matrix <paramref name="matA"/> and matrix <paramref name="matB"/>.
        /// </summary>
        /// <returns>Returns the product of two matrices.</returns>
        /// <param name="matA">First matrix.</param>
        /// <param name="matB">second matrix.</param>
        private T[][] MatProduct(T[][] matA,
                                 T[][] matB)
        {
            int aRows = matA.Length;
            int aCols = matA[0].Length;
            int bRows = matB.Length;
            int bCols = matB[0].Length;
            if (aCols != bRows)
                throw new Exception("Non-conformable matrices");

            T[][] result = Matrix<T>.CreateJaggedArray(aRows, bCols); // MatCreate(aRows, bCols);

            for (int i = 0; i < aRows; ++i) 
                for (int j = 0; j < bCols; ++j)
                    for (int k = 0; k < aCols; ++k) 
                        result[i][j] = this.Op.Add(result[i][j], this.Op.Mult( matA[i][k], matB[k][j] ) );

            return result;
        }


        /// <summary>
        /// Multiplies two matrices, matrix <paramref name="matA"/> and matrix <paramref name="matB"/>
        /// using parallelism.
        /// </summary>
        /// <returns>Returns the product of two matrices.</returns>
        /// <param name="matA">First matrix.</param>
        /// <param name="matB">second matrix.</param>
        private T[][] MatProductPar(T[][] matA,
                                 T[][] matB)
        {
            int aRows = matA.Length;
            int aCols = matA[0].Length;
            int bRows = matB.Length;
            int bCols = matB[0].Length;
            if (aCols != bRows)
                throw new Exception("Non-conformable matrices");

            T[][] result = Matrix<T>.CreateJaggedArray(aRows, bCols); // MatCreate(aRows, bCols);

            Parallel.For(0, aRows, delegate(int i)
                {
                    for (int j = 0; j < bCols; ++j)
                        for (int k = 0; k < aCols; ++k)
                            result[i][j] = this.Op.Add(result[i][j], this.Op.Mult(matA[i][k], matB[k][j]));
                }
            );

            return result;
        }


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
        public int MatDecompQR( T[][] A,
                                       out T[][] Q, out T[][] R )
        {
            // work with rows of the transpose
            // return another transpose at end
            int rows = Matrix<T>.GetNumRows(A);
            int cols = Matrix<T>.GetNumCols(A);

            T[][] a = Matrix<T>.Transpose(rows, cols, A); //MatTranspose(A);
            T[][] u = Matrix<T>.Clone(a);  //MatDuplicate(a);
            //int rows = a.Length;  // of the transpose
            //int cols = a[0].Length;

            Q = Matrix<T>.CreateJaggedArray(cols, rows);  //MatCreate(cols, rows);
            R = Matrix<T>.CreateJaggedArray(cols, rows);  //MatCreate(cols, rows);

            // first row of a (first col of M)
            for (int j = 0; j < cols; ++j)
                u[0][j] = a[0][j];

            T[] accum = new T[cols];
            // remaining rows of a
            for (int i = 1; i < rows; ++i)  
            {
                for (int j = 0; j < cols; ++j)
                {
                    // accumulate projections
                    accum = new T[cols];
                    for (int t = 0; t < i; ++t)
                    {
                        T[] proj = VecProjection(u[t], a[i]);
                        for (int k = 0; k < cols; ++k)
                            accum[k] = this.Op.Add(accum[k], proj[k]);
                    }
                }
                for (int k = 0; k < cols; ++k)
                    u[i][k] = this.Op.Sub( a[i][k], accum[k] );
            }

            for (int i = 0; i < rows; ++i)
            {
                T norm = VecNorm(u[i]);
                for (int j = 0; j < cols; ++j)
                    u[i][j] = this.Op.Div( u[i][j], norm );
            }
            // at this point u is Q(trans)

            T[][] q = Matrix<T>.Transpose( Matrix<T>.GetNumRows(u), Matrix<T>.GetNumCols(u), u);  // MatTranspose(u);
            for (int i = 0; i < q.Length; ++i)
                for (int j = 0; j < q[0].Length; ++j)
                    Q[i][j] = q[i][j];

            T[][] r =  this.MatProduct(u, A);
            for (int i = 0; i < r.Length; ++i)
                for (int j = 0; j < r[0].Length; ++j)
                    R[i][j] = r[i][j];

            return 0;
        } // 


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
        public int MatDecompQRPar(T[][] A,
                                       out T[][] Q, out T[][] R)
        {
            // work with rows of the transpose
            // return another transpose at end
            int rows = Matrix<T>.GetNumRows(A);
            int cols = Matrix<T>.GetNumCols(A);

            T[][] a = Matrix<T>.TransposePar(rows, cols, A); //MatTranspose(A);
            T[][] u = Matrix<T>.ClonePar(a);  //MatDuplicate(a);
            //int rows = a.Length;  // of the transpose
            //int cols = a[0].Length;

            Q = Matrix<T>.CreateJaggedArrayPar(cols, rows);  //MatCreate(cols, rows);
            R = Matrix<T>.CreateJaggedArrayPar(cols, rows);  //MatCreate(cols, rows);

            // first row of a (first col of M)
            Parallel.For(0, cols, delegate (int j)
                {
                    u[0][j] = a[0][j];
                }
            );

            T[] accum = new T[cols];
            // remaining rows of a
            Parallel.For(1, rows, delegate (int i)
            //for (int i = 1; i < rows; ++i)
                {
                    for (int j = 0; j < cols; ++j)
                    {
                        // accumulate projections
                        accum = new T[cols];
                        for (int t = 0; t < i; ++t)
                        {
                            T[] proj = VecProjectionPar(u[t], a[i]);
                            for (int k = 0; k < cols; ++k)
                                accum[k] = this.Op.Add(accum[k], proj[k]);
                        }
                    }
                    for (int k = 0; k < cols; ++k)
                        u[i][k] = this.Op.Sub(a[i][k], accum[k]);
                }
            );

            //for (int i = 0; i < rows; ++i)
            Parallel.For(0, rows, delegate (int i) 
                {
                    T norm = VecNorm(u[i]);
                    for (int j = 0; j < cols; ++j)
                        u[i][j] = this.Op.Div(u[i][j], norm);
                }
            );
               
            // at this point u is Q(trans)
            T[][] q = Matrix<T>.TransposePar(Matrix<T>.GetNumRows(u), Matrix<T>.GetNumCols(u), u);  // MatTranspose(u);
            T[][] Q1 = Q;
            Parallel.For(0, q.Length, delegate (int i) 
                {
                    for (int j = 0; j < q[0].Length; ++j)
                        Q1[i][j] = q[i][j];
                }
            );
                
            T[][] r = this.MatProductPar(u, A);
            T[][] R1 = R;
            Parallel.For(0, r.Length, delegate(int i) 
                {
                    for (int j = 0; j < r[0].Length; ++j)
                        R1[i][j] = r[i][j];
                }
            );
                
            return 0;
        } // 



    }

}
