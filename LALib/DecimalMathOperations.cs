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
using LALib.DecMath;

namespace LALib
{
    /// <summary>
    /// Basic decimal math operations.
    /// </summary>
    public class DecimalMathOperations : IBasicMathOperations<decimal>
    {
        public DecimalMathOperations()
        {
        }

        public decimal Round(decimal a, int numDec)
        {
            return Math.Round(a, numDec);
        }

        public decimal Minus(decimal a)
        {
            return -a;
        }

        public bool IsLesserThan(decimal a, decimal b)
        {
            return (a < b);
        }

        public bool IsGreaterThan(decimal a, decimal b)
        {
            return (a > b);
        }

        public decimal Abs(decimal value)
        {
            return Math.Abs(value);
        }

        public decimal Zero
        {
            get { return Decimal.Zero; }
        }

        public decimal One
        {
            get { return Decimal.One; }
        }

        public decimal Add(decimal a, decimal b)
        {
            return decimal.Add( a, b );
        }

        public decimal Div(decimal a, decimal b)
        {
            return a / b;
        }

        public decimal Mult(decimal a, decimal b)
        {
            return a * b;
        }

        public decimal Sub(decimal a, decimal b)
        {
            return a - b;
        }

        public decimal Conjugate(decimal a)
        {
            return a;
        }

        public double SqrtDbl(decimal a)
        {
            return Decimal.ToDouble( DecimalMath.DecimalSqrt(a) ); //Math.Sqrt(Decimal.ToDouble(a));
        }


        public decimal Sqrt(decimal a)
        {
            return DecimalMath.DecimalSqrt(a);
        }


        public double ToDouble(decimal a)
        {
            return Decimal.ToDouble(a);
        }


    }

}
