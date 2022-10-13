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
 * 
 * 
 ********************************************************************************/


using System;
using System.Threading.Tasks;

namespace LALib
{

    /// <summary>
    /// Computes the rank of a matrix.
    /// </summary>
    public class RankHelper<T>
    {

        private IBasicMathOperations<T> Op { get; set; }


        public RankHelper(IBasicMathOperations<T> _op)
        {
            this.Op = _op; 
        }


        /// <summary>
        /// Function for exchanging two rows of a matrix.
        /// </summary>
        /// <param name="mat">Original matrix.</param>
        /// <param name="row1">First row to exchange.</param>
        /// <param name="row2">Second row to exchange.</param>
        private static void Swap(T[][] mat,
              int row1, int row2)
        {
            T[] temp = mat[row1];
            mat[row1] = mat[row2];
            mat[row2] = temp;
        }


        /// <summary>
        /// Function for finding rank of matrix.
        /// </summary>
        /// <returns>Returns the rank of given matrix <paramref name="mat"/>.</returns>
        /// <param name="mat">Source matrix to find its rank.</param>
        /// <typeparam name="T">The 1st type parameter.</typeparam>
        /// <remarks>
        /// 
        ///     In linear algebra, the rank of a matrix A is the dimension of the vector
        ///     space generated (or spanned) by its columns. This corresponds to the 
        ///     maximal number of linearly independent columns of A. This, in turn,
        ///     is identical to the dimension of the vector space spanned by its 
        ///     rows. Rank is thus a measure of the "nondegenerateness" of the 
        ///     system of linear equations and linear transformation encoded by A.
        ///     There are multiple equivalent definitions of rank. A matrix's 
        ///     rank is one of its most fundamental characteristics.
        /// 
        ///     The rank is commonly denoted by rank(A) or rk(A); sometimes the 
        ///     parentheses are not written, as in rank A.
        /// 
        /// </remarks>
        public int RankOfMatrix(T[][] mat)
        {
            int R = mat.Length;
            int C = mat[0].Length;

            // min of R x C
            int rank =  (C < R ) ? C : R;
            T val = this.Op.Zero;
            T mult = this.Op.Zero; //0.0;
            //dynamic v1 = 0;
            //dynamic v2 = 0;

            for (int row = 0; row < rank; row++)
            {
                // Before we visit current row
                // 'row', we make sure that
                // mat[row][0],....mat[row][row-1]
                // are 0.

                val = mat[row][row];

                // Diagonal element is not zero
                if (!(val.Equals(this.Op.Zero)))  // != 0)
                {
                    for (int col = 0; col < R; col++)
                    {
                        if (col != row)
                        {
                            // This makes all entries
                            // of current column
                            // as 0 except entry
                            // 'mat[row][row]'
                            //v1 = mat[col][row]; v2 = mat[row][row];
                            //mult = v1 / v2;

                            mult = this.Op.Div( mat[col][row], mat[row][row] );

                            for (int i = 0; i < rank; i++)
                            {
                                //v1 = mat[row][i];
                                //mat[col][i] -= mult
                                             //* v1;

                                mat[col][i] = this.Op.Sub(mat[col][i], this.Op.Mult(mult, mat[row][i]) );

                            }
                        }
                    }
                }

                // Diagonal element is already zero.
                // Two cases arise:
                // 1) If there is a row below it
                // with non-zero entry, then swap
                // this row with that row and process
                // that row
                // 2) If all elements in current
                // column below mat[r][row] are 0,
                // then remove this column by
                // swapping it with last column and
                // reducing number of columns by 1.
                else
                {
                    bool reduce = true;

                    // Find the non-zero element
                    // in current column
                    for (int i = row + 1; i < R; i++)
                    {
                        val = mat[i][row];
                        // Swap the row with non-zero
                        // element with this row.
                        if (!(val.Equals(this.Op.Zero) ) )  // != 0)
                        {
                            RankHelper<T>.Swap(mat, row, i);//rank);
                            reduce = false;
                            break;
                        }
                    }

                    // If we did not find any row with
                    // non-zero element in current
                    // column, then all values in
                    // this column are 0.
                    if (reduce)
                    {
                        // Reduce number of columns
                        rank--;

                        // Copy the last column here
                        for (int i = 0; i < R; i++)
                            mat[i][row] = mat[i][rank];
                    }

                    // Process this row again
                    row--;
                }

                // Uncomment these lines to see
                // intermediate results display(mat, R, C);
                // printf("\n");
            }

            return rank;
        }


        /// <summary>
        /// Function for finding rank of matrix using parallelism.
        /// </summary>
        /// <returns>Returns the rank of given matrix <paramref name="mat"/>.</returns>
        /// <param name="mat">Source matrix to find its rank.</param>
        /// <typeparam name="T">The 1st type parameter.</typeparam>
        /// <remarks>
        /// 
        ///     In linear algebra, the rank of a matrix A is the dimension of the vector
        ///     space generated (or spanned) by its columns. This corresponds to the 
        ///     maximal number of linearly independent columns of A. This, in turn,
        ///     is identical to the dimension of the vector space spanned by its 
        ///     rows. Rank is thus a measure of the "nondegenerateness" of the 
        ///     system of linear equations and linear transformation encoded by A.
        ///     There are multiple equivalent definitions of rank. A matrix's 
        ///     rank is one of its most fundamental characteristics.
        /// 
        ///     The rank is commonly denoted by rank(A) or rk(A); sometimes the 
        ///     parentheses are not written, as in rank A.
        /// 
        /// </remarks>
        public int RankOfMatrixPar(T[][] mat)
        {
            int R = mat.Length;
            int C = mat[0].Length;

            // min of R x C
            int rank = (C < R) ? C : R;
            T val = this.Op.Zero;
            T mult = this.Op.Zero; //0.0;
            //dynamic v1 = 0;
            //dynamic v2 = 0;

            for (int row = 0; row < rank; row++)
            {
                // Before we visit current row
                // 'row', we make sure that
                // mat[row][0],....mat[row][row-1]
                // are 0.

                val = mat[row][row];

                // Diagonal element is not zero
                if (!(val.Equals(this.Op.Zero)))  // != 0)
                {
                    Parallel.For(0, R, delegate( int col )
                        {
                            if (col != row)
                            {
                                // This makes all entries
                                // of current column
                                // as 0 except entry
                                // 'mat[row][row]'
                                //v1 = mat[col][row]; v2 = mat[row][row];
                                //mult = v1 / v2;

                                mult = this.Op.Div(mat[col][row], mat[row][row]);

                                for (int i = 0; i < rank; i++)
                                {
                                    //v1 = mat[row][i];
                                    //mat[col][i] -= mult
                                    //* v1;

                                    mat[col][i] = this.Op.Sub(mat[col][i], this.Op.Mult(mult, mat[row][i]));

                                }
                            }
                        }
                    );
                    //for (int col = 0; col < R; col++)
                }

                // Diagonal element is already zero.
                // Two cases arise:
                // 1) If there is a row below it
                // with non-zero entry, then swap
                // this row with that row and process
                // that row
                // 2) If all elements in current
                // column below mat[r][row] are 0,
                // then remove this column by
                // swapping it with last column and
                // reducing number of columns by 1.
                else
                {
                    bool reduce = true;

                    // Find the non-zero element
                    // in current column
                    for (int i = row + 1; i < R; i++)
                    {
                        val = mat[i][row];
                        // Swap the row with non-zero
                        // element with this row.
                        if (!(val.Equals(this.Op.Zero)))  // != 0)
                        {
                            RankHelper<T>.Swap(mat, row, i);//rank);
                            reduce = false;
                            break;
                        }
                    }

                    // If we did not find any row with
                    // non-zero element in current
                    // column, then all values in
                    // this column are 0.
                    if (reduce)
                    {
                        // Reduce number of columns
                        rank--;

                        // Copy the last column here
                        for (int i = 0; i < R; i++)
                            mat[i][row] = mat[i][rank];
                    }

                    // Process this row again
                    row--;
                }

                // Uncomment these lines to see
                // intermediate results display(mat, R, C);
                // printf("\n");
            }

            return rank;
        }



    }


}
