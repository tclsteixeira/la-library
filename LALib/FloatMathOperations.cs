/********************************************************************************
 *  Copyright 2022, Tiago C. Teixeira
 * 
 *  MIT license (Massachusetts Institute of Technology)
 * 
 *  Permission is hereby granted, free of charge, to any person obtaining a copy
 *  of this software and associated documentation files (the "Software"), to deal 
 *  in the Software without restriction, including without limitation the rights 
 *  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies 
 *  of the Software, and to permit persons to whom the Software is furnished to 
 *  do so, subject to the following conditions:
 * 
 *  The above copyright notice and this permission notice shall be included in all 
 *  copies or substantial portions of the Software.
 * 
 *  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
 *  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
 *  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
 *  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
 *  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING 
 *  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS 
 *  IN THE SOFTWARE.
 *
 ********************************************************************************/

using System;
namespace LALib
{
    public class FloatMathOperations : IBasicMathOperations<float>
    {

        public FloatMathOperations()
        {
        }


        public float Round(float a, int numDec)
        {
            return (float)Math.Round(a, numDec);
        }

        public float Minus(float a)
        {
            return -a;
        }

        public bool IsLesserThan(float a, float b)
        {
            return (a < b);
        }

        public bool IsGreaterThan(float a, float b)
        {
            return (a > b);
        }

        public float Abs(float value)
        {
            return Math.Abs(value);
        }

        public float Zero
        {
            get { return 0.0f; }
        }

        public float One
        {
            get { return 1.0f; }
        }

        public float Add(float a, float b)
        {
            return a + b;
        }

        public float Div(float a, float b)
        {
            return a / b;
        }

        public float Mult(float a, float b)
        {
            return a * b;
        }

        public float Sub(float a, float b)
        {
            return a - b;
        }

        public float Conjugate(float a)
        {
            return a;
        }

        public double SqrtDbl(float a)
        {
            return Math.Sqrt(a);
        }

        public float Sqrt(float a)
        {
            return (float)Math.Sqrt(a);
        }

        public double ToDouble(float a)
        {
            return a;
        }

    }
}

