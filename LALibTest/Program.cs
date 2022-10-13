using System;
using System.Numerics;
using LALib;

namespace LALibTest
{

    class MainClass
    {
        /// <summary>
        /// The entry point of the program, where the program control starts and ends.
        /// Performs a test to several matrix functions of LALib library.
        /// </summary>
        /// <param name="args">The command-line arguments.</param>
        public static void Main(string[] args)
        {
            Console.WriteLine("LALib tests:");

            Console.WriteLine("Let A (2x3) matrix and B (3x2) matrix,");
            Console.WriteLine();

            Int32Matrix i32AR = new Int32Matrix();
            ComplexMatrix cpxMat = new ComplexMatrix();

            int[][] A = Matrix<int>.CreateJaggedArray(2, 3);
            int[][] B = Matrix<int>.CreateJaggedArray(3, 2);

            int n = 2; // A.GetLength(0);     // num rows A
            int m = 3; // A[0].GetLength(0);     // num cols A
            int p = 2; //B[0].GetLength(0);     // Num cols B (num rows B must be equal to num cols A)

            Console.WriteLine("Create randow matrix A:");
            A = Int32Matrix.Random(n, m, -10, 10, 1254);
            //FillRandom(A);

            Console.WriteLine("Matrix A:");
            PrintMatrix<int>(A);

            Console.WriteLine();
            Console.WriteLine("Create randow matrix B:");
            B = Int32Matrix.Random(m, p, -10, 10, 1124);

            Console.WriteLine("Matrix B:");
            PrintMatrix<int>(B);

            Console.WriteLine();
            Console.WriteLine();

            Console.WriteLine("Matrix C=AxB:");

            int[][] C = i32AR.MultPar( A, B );

            PrintMatrix(C);

            Console.WriteLine();

            // until aprox. 50000 cells single core execution seems faster than parallel.

            n = 500; m = 500;
            p = 500;

            //A = Matrix<int>.CreateJaggedArray(n, m);
            //B = Matrix<int>.CreateJaggedArray(m, p);

            Console.WriteLine("Creating random matrix A and B 500x500:");

            A = Int32Matrix.Random(n, m, -10, 10, 1222 );
            B = Int32Matrix.Random(m, p, -10, 10, 1054);

            Console.WriteLine($"Multiply A ({n},{m}) by B ({m},{p}) using parallel multi-core execution:");

            DateTimeOffset dtIni = DateTimeOffset.Now;
            C = i32AR.MultPar(A, B);
            TimeSpan ts = TimeSpan.FromTicks(DateTimeOffset.Now.Ticks - dtIni.Ticks);

            Console.WriteLine($"Elapsed time: {ts.ToString()}.");

            Console.WriteLine();

            Console.WriteLine($"Multiply A ({n},{m}) by B ({m},{p}) using only single core execution:");

            DateTimeOffset dtIni2 = DateTimeOffset.Now;
            C = i32AR.Mult( A, B );
            ts = TimeSpan.FromTicks(DateTimeOffset.Now.Ticks - dtIni2.Ticks);

            Console.WriteLine($"Elapsed time: {ts.ToString()}.");

            Console.WriteLine();
            Console.WriteLine();

            dtIni = DateTimeOffset.Now;

            Console.WriteLine($"Clone Matrix A ({n}x{m}) (no parallelism):");

            // Clone A
            int[][] AClone = Matrix<int>.Clone(A);

            ts = TimeSpan.FromTicks(DateTimeOffset.Now.Ticks - dtIni.Ticks);
            Console.WriteLine($"Elapsed time: {ts.ToString()}.");

            Console.WriteLine();

            dtIni = DateTimeOffset.Now;

            Console.WriteLine($"Clone Matrix B ({m}x{p}) (with parallelism):");

            // Clone A
            int[][] BClone = Matrix<int>.Clone(B);

            ts = TimeSpan.FromTicks(DateTimeOffset.Now.Ticks - dtIni.Ticks);
            Console.WriteLine($"Elapsed time: {ts.ToString()}.");

            Console.WriteLine();

            // Transpose
            //A = Matrix<int>.CreateJaggedArray(2, 3);
            A = Int32Matrix.Random(2, 3, -10, 10, 1264);
            B = Matrix<int>.Transpose(2, 3, A);

            Console.WriteLine("Transpose A (2x3) and save result in B");

            Console.WriteLine("Matrix A:");
            PrintMatrix<int>(A);

            Console.WriteLine();
            Console.WriteLine("B = Transpose(A):");
            PrintMatrix<int>(B);


            Console.WriteLine();
            Console.WriteLine("A = random double matrix (3x3) - Compute determinant of A = det(A):");

            //double[][] Adbl = Matrix<double>.CreateJaggedArray(3, 3);
            double[][] Adbl = DoubleMatrix.Random(3, 3, 0, 1, 1334);

            PrintMatrix<double>(Adbl);

            Console.WriteLine();

            DoubleMatrix dblmat = new DoubleMatrix();

            //double detA = Matrix.DETInPlace(Adbl, 3);
            double detA = dblmat.DET_BareissAlg(Adbl, 3, false);

            Console.WriteLine($"Det(A) using Bareiss alg = {detA}");

            //PrintMatrix<double>(Adbl);

            Console.WriteLine();
            Console.WriteLine("Identity matrix size 10");
            PrintMatrix<double>(dblmat.IdentityPar(10));

            // Now test random matrices

            Console.WriteLine();
            Console.WriteLine("Integer random 5 x 5 matrix with range from -100 to 100");
            PrintMatrix<int>(Int32Matrix.Random( 5, 5, -100, 100, 1245 ));

            // Linear equation solver
            Console.WriteLine();

            double[][] coefs = Matrix<double>.CreateJaggedArray(4, 4);
            coefs[0][0] = 1; coefs[0][1] = 1; coefs[0][2] = -3; coefs[0][3] = 1;
            coefs[1][0] = -5; coefs[1][1] = 3; coefs[1][2] = -4; coefs[1][3] = 1;
            coefs[2][0] = 1; coefs[2][1] = 0; coefs[2][2] = 2; coefs[2][3] = -1;
            coefs[3][0] = 1; coefs[3][1] = 2; coefs[3][2] = 0; coefs[3][3] = 0;

            double[] b = new double[4] { 2, 0, 1, 12 };

            // Solve it ...
            double[] X = dblmat.LUSolve( coefs, b, false );

            Console.WriteLine("-------------------------------------");
            Console.WriteLine("Linear equation solver using LU method");
            Console.WriteLine();
            Console.WriteLine("Coeficients matrix:");
            PrintMatrix<double>(coefs);

            Console.WriteLine();
            Console.WriteLine("Second member terms:");
            PrintVector<double>(b);

            Console.WriteLine();
            Console.WriteLine("Solutions X vector:");
            PrintVector<double>(X);

            // Using parallelization
            Console.WriteLine();
            Console.WriteLine("-------------------------------------");
            Console.WriteLine("Linear equation solver with parallelization using LU method");
            Console.WriteLine();
            Console.WriteLine("Coeficients matrix:");
            PrintMatrix<double>(coefs);

            Console.WriteLine();
            Console.WriteLine("Second member terms:");
            PrintVector<double>(b);

            // Solve it with parallelization ----------
            double[] Xpar = dblmat.LUSolvePar(coefs, b, false);

            Console.WriteLine();
            Console.WriteLine("Solutions X vector:");
            PrintVector<double>(Xpar);
            //------------------------------------------


            Console.WriteLine();
            Console.WriteLine("--------------------------");
            Console.WriteLine("LU lower triangular matrix:");
            double[][] lower = dblmat.LULower( coefs );
            PrintMatrix<double>(lower);

            Console.WriteLine();
            Console.WriteLine("LU upper triangular matrix:");
            double[][] upper = dblmat.LUUpper( coefs );
            PrintMatrix<double>(upper);

            Console.WriteLine();
            Console.WriteLine("Prof -> L*U matrix:");
            double[][] LU = dblmat.Mult(lower, upper);
            PrintMatrix<double>(LU);

            Console.WriteLine();
            Console.WriteLine("--------------------------");
            Console.WriteLine("LU lower triangular matrix with parallelization:");
            double[][] lowerPar = dblmat.LULowerPar(coefs);
            PrintMatrix<double>(lowerPar);

            Console.WriteLine();
            Console.WriteLine("LU upper triangular matrix using Parallelization:");
            double[][] upperPar = dblmat.LUUpperPar(coefs);
            PrintMatrix<double>(upperPar);

            Console.WriteLine();
            Console.WriteLine("Prof -> L*U matrix with parallelization:");
            double[][] LUPar = dblmat.MultPar(lowerPar, upperPar);
            PrintMatrix<double>(LUPar);

            Console.WriteLine();
            Console.WriteLine("--------------------------");
            Console.WriteLine("Sum of rows of coefs matrix:");
            double[] sum = dblmat.SumRows(coefs);
            PrintVector<double>(sum);

            Console.WriteLine();
            Console.WriteLine("Sum of columns of coefs matrix:");
            sum = dblmat.SumCols(coefs);
            PrintVector<double>(sum);

            double[][] invCoefs = dblmat.Inverse(coefs);

            Console.WriteLine();
            Console.WriteLine("-------------------------------------");
            Console.WriteLine("Inverse of coeficients matrix");
            Console.WriteLine();
            Console.WriteLine("Inverse matrix:");
            PrintMatrix<double>(invCoefs);

            Console.WriteLine();
            Console.WriteLine("Coeficients matrix (A) times his inverse (A*A^(-1):");

            // must be equal to 4x4 identity matrix
            double[][] invProf = dblmat.Mult( coefs, invCoefs );
            invProf = Matrix<double>.AdjustToValue( 0, invProf );

            PrintMatrix<double>(invProf);

            Console.WriteLine();
            Console.WriteLine("-------------------------------------");
            Console.WriteLine("Inverse of coeficients matrix using parallelism");
            Console.WriteLine();
            Console.WriteLine("Inverse matrix:");
            double[][] invCoefsPar = dblmat.InversePar(coefs);
            PrintMatrix<double>(invCoefsPar);

            Console.WriteLine();
            Console.WriteLine("Coeficients matrix (A) times his inverse (A*A^(-1):");

            // must be equal to 4x4 identity matrix
            double[][] invProfPar = dblmat.MultPar(coefs, invCoefsPar);
            invProfPar = Matrix<double>.AdjustToValuePar(0, invProfPar);
            PrintMatrix<double>(invProfPar);

            Console.WriteLine();
            Console.WriteLine("-------------------------------------");
            Console.WriteLine("Rank of matrix 'r':");
            int[][] r = Int32Matrix.Random(5, 5, -3, 3, 9234);// new double[2][] { new double[] { -1, -3 }, new double[] { 1, 3 } };
            int[][] r2 = Matrix<int>.Clone(r);// Int32Matrix.Random(10, 10, -5, 5, 234);// new double[2][] { new double[] { -1, -3 }, new double[] { 1, 3 } };

            //double[][] r = new double[2][] { new double[] { -1, -3 }, new double[] { 1, 3 } };
            PrintMatrix<int>(r);
            Console.WriteLine();
            Console.WriteLine($"Rank = {i32AR.Rank(r)}");

            Console.WriteLine();
            Console.WriteLine("-------------------------------------");
            Console.WriteLine("Rank of matrix 'r' using parallelism:");
            //double[][] rpar = new double[2][] { new double[] { -1, -3 }, new double[] { 1, 3 } };
            PrintMatrix<int>(r2);
            Console.WriteLine();
            Console.WriteLine($"Rank = {i32AR.RankPar(r2)}");

            Console.WriteLine();
            Console.WriteLine("-------------------------------------");
            Console.WriteLine("Create Hilbert matrix 5x5:");
            double[][] hilb = Matrix<double>.CreateHilbertMatrix(5);
            PrintMatrix<double>(hilb);

            Console.WriteLine();
            Console.WriteLine("-------------------------------------");
            Console.WriteLine("Create matrix from array whose elements represents rows:");
            int[] elRows = new int[3] { 1, 2, 3 };
            PrintArray<int>(elRows);
            int[][] rMat = Matrix<int>.MatrixFromRowsVector(elRows);
            Console.WriteLine("Result matrix:");
            PrintMatrix<int>(rMat);

            Console.WriteLine();
            Console.WriteLine("-------------------------------------");
            Console.WriteLine("Create matrix from array whose elements represents columns:");
            int[] elCols = new int[3] { 4, 5, 6};
            PrintArray<int>(elCols);
            int[][] cMat = Matrix<int>.MatrixFromColumnsVector(elCols);
            Console.WriteLine("Result matrix:");
            PrintMatrix<int>(cMat);


            Console.WriteLine();
            Console.WriteLine("-------------------------------------");
            Console.WriteLine("Matrix exponentiation - Matrix A^3:");
            A = Int32Matrix.Random(3, 3, -10, 10, 2321);
            Console.WriteLine("Matrix A:");
            PrintMatrix<int>(A);
            Console.WriteLine();
            Console.WriteLine("Matrix A^3:");
            B = i32AR.Pow(A, 3);
            PrintMatrix<int>(B);

            Console.WriteLine();
            Console.WriteLine("-------------------------------------");
            Console.WriteLine("Complex numbers - Matrix M^3:");
            Complex[][] M = ComplexMatrix.Random(3, 3, -10, 10, 2621);
            Console.WriteLine("Matrix M:");
            PrintMatrix<Complex>(M);
            Console.WriteLine();
            Console.WriteLine("Matrix M^3:");
            Complex[][] M2 = cpxMat.Pow(M, 3);
            PrintMatrix<Complex>(M2);

            Console.WriteLine();
            Console.WriteLine("-- NORMS -----------------------------");
            Console.WriteLine("Complex matrix 'N1' 3x3:");
            Complex[][] N1 = Matrix<Complex>.CreateJaggedArray( 3, 3 );
            N1[0] = new Complex[] { new Complex(2,1), new Complex(7.1, 0.4), new Complex(-2.5, 3) };
            N1[1] = new Complex[] { new Complex(-2, -1), new Complex(9.1, -5.4), new Complex(-2.9, -3) };
            N1[2] = new Complex[] { new Complex(1, -2), new Complex(-3.7, 0.9), new Complex(2.5, -5) };
            Console.WriteLine("Matrix N1:");
            PrintMatrix<Complex>(N1);
            Console.WriteLine();
            Console.WriteLine($"Induced L1 norm: {cpxMat.L1Norm(N1)}");
            Console.WriteLine($"Frobenius norm: {cpxMat.FrobeniusNorm(N1)}");
            Console.WriteLine($"Infinity norm: {cpxMat.InfinityNorm(N1)}");
            Console.WriteLine("--------_-----------------------------");

            Console.WriteLine();
            Console.WriteLine("-- INNER (DOT) PRODUCT ---------------");
            Console.WriteLine("Complex matrix 'N2' 3x3:");
            Complex[][] N2 = Matrix<Complex>.CreateJaggedArray(3, 3);
            N2[0] = new Complex[] { new Complex(-4, 2), new Complex(2.1, 9.4), new Complex(-0.5, 8) };
            N2[1] = new Complex[] { new Complex(2, -1.6), new Complex(8, -4), new Complex(-9, -3) };
            N2[2] = new Complex[] { new Complex(1.3, -2.6), new Complex(-3, -0.9), new Complex(2.1, 5) };
            Console.WriteLine("Matrix N2:");
            PrintMatrix<Complex>(N2);
            Console.WriteLine();
            Console.WriteLine($"Frobenius inner (dot) product (N1.N2): {cpxMat.FrobeniusInnerProd(N1, N2)}");
            Console.WriteLine("-------------------------------------");

            Console.WriteLine();
            Console.WriteLine("-- QR Decomposition ---------------");
            Console.WriteLine("-------------------------------------");
            Console.WriteLine("Double matrix 'D' 3x3:");
            double[][] D = Matrix<double>.CreateJaggedArray(3, 3);
            D[0] = new double[] { 12, -51, 4 };
            D[1] = new double[] { 6,  167,  -68};
            D[2] = new double[] { -4, 24, -41};
            Console.WriteLine("Matrix D:");
            PrintMatrix<double>(D);

            if (dblmat.QRDecomp(D, out double[][] Q, out double[][] R) == 0)
            {
                Console.WriteLine();
                Console.WriteLine("Matrix Q:");
                Q = Matrix<double>.AdjustToValue(0, Q);
                PrintMatrix<double>(Q);
                Console.WriteLine();
                Console.WriteLine("Matrix R:");
                R = Matrix<double>.AdjustToValue(0, R, 1E-13);
                PrintMatrix<double>(R);

                Console.WriteLine();
                Console.WriteLine("Prof: D = QxR:");
                Console.WriteLine("Matrix QxR:");
                PrintMatrix<double>(dblmat.Mult(Q, R));
            }
            else
                Console.WriteLine("Failed to compute QR decomposition on matrix 'D'.");

            Console.WriteLine();
            Console.WriteLine("-- QR Decomposition using parallelism ---------------");
            Console.WriteLine("-------------------------------------");
            Console.WriteLine("Double matrix 'D' 3x3:");
            double[][] Dpar = Matrix<double>.CreateJaggedArrayPar(3, 3);
            Dpar[0] = new double[] { 12, -51, 4 };
            Dpar[1] = new double[] { 6, 167, -68 };
            Dpar[2] = new double[] { -4, 24, -41 };
            Console.WriteLine("Matrix D:");
            PrintMatrix<double>(Dpar);

            if (dblmat.QRDecompPar(Dpar, out double[][] Qpar, out double[][] Rpar) == 0)
            {
                Console.WriteLine();
                Console.WriteLine("Matrix Q:");
                Qpar = Matrix<double>.AdjustToValuePar(0, Qpar);
                PrintMatrix<double>(Qpar);
                Console.WriteLine();
                Console.WriteLine("Matrix R:");
                Rpar = Matrix<double>.AdjustToValue(0, Rpar, 1E-13);
                PrintMatrix<double>(Rpar);

                Console.WriteLine();
                Console.WriteLine("Prof: D = QxR:");
                Console.WriteLine("Matrix QxR:");
                PrintMatrix<double>(dblmat.MultPar(Qpar, Rpar));
            }
            else
                Console.WriteLine("Failed to compute QR decomposition using parallelism on matrix 'D'.");
               
            Console.WriteLine();
            Console.WriteLine("------------------------------------");
            Console.WriteLine("Create diagonal matrix from vector:");
            int[] v1 = new int[3] { 1, 2, 3 };
            Console.WriteLine("Vector a:");
            PrintVector<int>(v1);
            Console.WriteLine();
            int[][] diag = Matrix<int>.DiagonalFromVector(v1);
            Console.WriteLine("Diagonal matrix D:");
            PrintMatrix<int>(diag);

            Console.WriteLine();
            Console.WriteLine("------------------------------------");
            Console.WriteLine("Create vector from diagonal matrix elements:");
            int[][] D2 = Int32Matrix.Random(4, 3, 0, 10, 64578);
            Console.WriteLine("Matrix D2:");
            PrintMatrix<int>(D2);
            Console.WriteLine();
            Console.WriteLine("Vector with 'D2' matrix diagonal elements:");
            v1 = Matrix<int>.DiagonalToVector(D2);
            PrintVector<int>(v1);

            //Console.WriteLine("Press any key to continue ...");
            //Console.ReadKey();
        }


        ///// <summary>
        ///// Fills a matrix of integers with random values.
        ///// </summary>
        ///// <param name="m1">The source matrix.</param>
        //private static void FillRandom(int[][] m1)
        //{
        //    int n = m1.GetLength(0);
        //    int m = m1[0].GetLength(0);
        //    Random rnd = new Random();

        //    for (int i = 0; i < n; i++)
        //    {
        //        for (int j = 0; j < m; j++)
        //        {
        //            m1[i][j] = rnd.Next(-10, 10);
        //        }
        //    }
        //}


        ///// <summary>
        ///// Fills a matrix of doubles with random values.
        ///// </summary>
        ///// <param name="m1">The source matrix.</param>
        //private static void FillRandom(double[][] m1)
        //{
        //    int n = m1.GetLength(0);
        //    int m = m1[0].GetLength(0);
        //    Random rnd = new Random();

        //    for (int i = 0; i < n; i++)
        //    {
        //        for (int j = 0; j < m; j++)
        //        {
        //            m1[i][j] = rnd.NextDouble();
        //        }
        //    }
        //}


        /// <summary>
        /// Prints the matrix.
        /// </summary>
        /// <param name="m1">The matrix instance.</param>
        /// <typeparam name="T">The matrix cell element type.</typeparam>
        private static void PrintMatrix<T>(T[][] m1)
        {
            if (m1 != null)
            {
                int spaces = 3;
                int n = m1.GetLength(0);
                int m = m1[0].GetLength(0);

                for (int i = 0; i < n; i++)
                {
                    for (int j = 0; j < m; j++)
                    {
                        Console.Write( m1[i][j].ToString());
                        Console.Write("".PadLeft(spaces));
                    }

                    Console.WriteLine();
                }
            }
        }


        /// <summary>
        /// Prints an array vector.
        /// </summary>
        /// <param name="m1">The array instance.</param>
        /// <typeparam name="T">The array cell element type.</typeparam>
        private static void PrintVector<T>(T[] m1)
        {
            if (m1 != null)
            {
                int spaces = 3;
                int n = m1.Length;

                for (int i = 0; i < n; i++)
                {
                    Console.Write(m1[i].ToString());
                    Console.Write("".PadLeft(spaces));

                    Console.WriteLine();
                }
            }
        }


        /// <summary>
        /// Prints an array.
        /// </summary>
        /// <param name="arr">The array instance.</param>
        /// <typeparam name="T">The array cell element type.</typeparam>
        private static void PrintArray<T>(T[] arr)
        {
            if (arr != null)
            {
                Console.Write("Array: ");
                Console.Write("{");
                string sep = ", ";
                int n = arr.Length;

                for (int i = 0; i < n; i++)
                {
                    Console.Write(arr[i].ToString());

                    if (i < (n-1))
                        Console.Write(sep);
                }

                Console.WriteLine("}");
            }
        }


    }
}
