/******************************************************************************
* DoubelMathOperations.cs
*
* Author:
*       Tiago C. Teixeira <>
*
* Copyright (c) 2022 
*
* Permission is hereby granted, free of charge, to any person obtaining a copy
* of this software and associated documentation files (the "Software"), to deal
* in the Software without restriction, including without limitation the rights
* to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the Software is
* furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included in
* all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
* OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
* THE SOFTWARE.
******************************************************************************/

using System;
namespace LALib
{
    public class DoubleMathOperations : IBasicMathOperations<double>
    {

        public DoubleMathOperations()
        {
        }


        public double Round(double a, int numDec)
        {
            return Math.Round(a, numDec);
        }

        public double Minus(double a)
        {
            return -a;
        }

        public bool IsLesserThan(double a, double b)
        {
            return (a < b);
        }

        public bool IsGreaterThan(double a, double b)
        {
            return (a > b);
        }

        public double Abs(double value)
        {
            return Math.Abs(value);
        }

        public double Zero
        {
            get { return 0.0d; }
        }

        public double One
        {
            get { return 1.0d; }
        }

        public double Add(double a, double b)
        {
            return a + b;
        }

        public double Div(double a, double b)
        {
            return a / b;
        }

        public double Mult(double a, double b)
        {
            return a * b;
        }

        public double Sub(double a, double b)
        {
            return a - b;
        }

        public double Conjugate(double a)
        {
            return a;
        }

        public double SqrtDbl(double a)
        {
            return Math.Sqrt(a);
        }

        public double Sqrt(double a)
        {
            return Math.Sqrt(a);
        }

        public double ToDouble(double a)
        {
            return a;
        }

    }

}
