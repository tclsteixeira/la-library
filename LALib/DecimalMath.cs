// /*********************************************************************************
//  * Copyright 2022, Tiago C. Teixeira
//  *
//  * MIT license (Massachusetts Institute of Technology)
//  *
//  * Permission is hereby granted, free of charge, to any person obtaining a copy
//  * of this software and associated documentation files (the "Software"), to deal 
//  * in the Software without restriction, including without limitation the rights 
//  * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies 
//  * of the Software, and to permit persons to whom the Software is furnished to 
//  * do so, subject to the following conditions:
//  *
//  * The above copyright notice and this permission notice shall be included in all 
//  * copies or substantial portions of the Software.
//  *
//  * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
//  * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
//  * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
//  * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
//  * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING 
//  * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS 
//  * IN THE SOFTWARE.
//  * 
//  ********************************************************************************/
//
using System;
namespace LALib.DecMath
{

    /// <summary>
    /// Math functions for double QUAD (128 bits) floating point numbers (decimal type).
    /// </summary>
    public static class DecimalMath
    {

        // Math constant 'e' with sixty decimal numbers precision
        public const decimal C_E_CONSTANT = 2.718281828459045235360287471352662497757247093699959574966967m;


        //public DecimalMath()
        //{
        //}


        /// <summary>
        /// Calculates the power of a base raised to a specified exponent for QUAD decimal numbers (128 bits).
        /// </summary>
        /// <param name="_base">The base to be raised.</param>
        /// <param name="_power">The power to raise the base.</param>
        /// <returns>Returns the power of a base raised to a given power if succeeded, throws an exception otherwise.</returns>
        public static decimal DecimalPow(decimal _base, decimal _power)
        {
            decimal FResult = decimal.Zero;
            try
            {
                FResult = DecimalExp(_power * DecimalLn(_base));
            }
            catch (Exception ex)
            {
                throw ex;
            }

            return FResult;
        }


        /// <summary>
        /// Calculates the exponential number (constant e) raised to a specifiec QUAD precision decimal number (128 bits).
        /// </summary>
        /// <param name="Exponent">The exponent.</param>
        /// <remarks>
        ///     Result = e^Exponent
        /// </remarks>
        /// <returns>Returns the exponential number (constant e) raised to a specifiec decimal number.</returns>
        public static decimal DecimalExp(decimal Exponent)
        {
            decimal FResult = decimal.Zero;
            try
            {
                decimal X, P, Frac, I, L;
                X = Exponent;
                Frac = X;
                P = (1.0m + X);
                I = 1.0m;

                do
                {
                    I++;
                    Frac *= (X / I);
                    L = P;
                    P += Frac;

                    //System.Diagnostics.Debug.WriteLine("L = " + L);
                    //System.Diagnostics.Debug.WriteLine("P = " + P);
                    //System.Diagnostics.Debug.WriteLine("L-P = " + (L - P));

                    //if (L-P == 0)
                    //    System.Diagnostics.Debug.WriteLine("L-P = " + (L - P));

                } while (L != P);

                FResult = P;
            }
            catch (Exception ex)
            {
                throw ex;
            }

            return FResult;
        }


        /// <summary>
        /// Calculates the natural logarithm of a specified QUAD precision decimal number (128 bits).
        /// </summary>
        /// <param name="arg">The function argument.</param>
        /// <remarks>
        ///     Result = Ln(arg)
        /// </remarks>
        /// <returns>Returns the natural logarithm.</returns>
        public static decimal DecimalLn(decimal arg)
        {
            decimal N, P, L, R, A, E;
            try
            {
                if (arg < 0)
                    throw new NotSupportedException(string.Format("Negative argument is not supported for {0} function.", "ln"));
                else
                {
                    E = C_E_CONSTANT;   // 2.71828182845905....;
                    P = arg;
                    N = 0.0m;

                    // This speeds up the convergence by calculating the integral
                    while (P >= E)
                    {
                        P /= E;
                        N++;
                    }

                    N += (P / E);
                    P = arg;

                    // Note: Sometimes the error can not converge to the exact zero because of precision roundings.
                    //       So we need to check if the error remains unchanged betwwen iterations to exit loop.
                    decimal error = decimal.Zero;
                    decimal newError = decimal.Zero;
                    do
                    {
                        error = newError;
                        A = N;
                        L = (P / (DecimalExp(N - 1.0m)));
                        R = ((N - 1.0m) * E);
                        N = ((L + R) / E);

                        //System.Diagnostics.Debug.WriteLine("A = " + A);
                        //System.Diagnostics.Debug.WriteLine("N = " + N);
                        //System.Diagnostics.Debug.WriteLine("A-N = " + (A - N));

                        if (A < N)
                            newError = N - A;
                        else
                            newError = A - N;

                        //if (A - N == 0)
                        //    System.Diagnostics.Debug.WriteLine("A-N = " + (A - N));

                        // If error remains equal no more convergence is possible.
                        if (newError == error)
                            break;

                    } while (N != A);
                }
            }
            catch (Exception ex)
            {
                throw ex;
            }

            return N;
        }


        /// <summary>
        /// Compute the square root for decimal type.
        /// </summary>
        /// <returns>Returns the square root.</returns>
        /// <param name="x">The square root argument.</param>
        public static decimal DecimalSqrt(decimal x)
        {
            decimal FResult = decimal.Zero;
            if (x > 0)
            {
                decimal root = x / 3;
                int i;
                for (i = 0; i < 32; i++)
                    root = (root + x / root) / 2;

                FResult = root;

                //Console.WriteLine("Square root: {0}", root);
                //Console.WriteLine("Math.Sqrt: {0}", Math.Sqrt((double)x));
            }
            else
            {
                throw new ArgumentException("Argument for square root must be a positive non zero value.");
            }

            return FResult;
        }


    }
}
