//
// IntMathOperations.cs
//
// Author:
//       Tiago C. Teixeira <>
//
// Copyright (c) 2022 
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

using System;
namespace LALib
{
    public class IntMathOperations : IBasicMathOperations<Int32>
    {
        public IntMathOperations()
        {
        }

        //public bool IsEqual(int a, int b)
        //{
        //    return (a == b);
        //}

        public int Minus(int a)
        {
            return -a;
        }

        public bool IsLesserThan(int a, int b)
        {
            return (a < b);
        }

        public bool IsGreaterThan(int a, int b)
        {
            return (a > b);
        }

        public int Abs(int value)
        {
            return Math.Abs(value);
        }

        public int Zero
        {
            get { return 0; }
        }

        public int One
        {
            get { return 1; }
        }

        public int Add(int a, int b)
        {
            return a + b;
        }

        public int Div(int a, int b)
        {
            Math.DivRem(a, b, out int result);
            return result;
        }

        public int Mult(int a, int b)
        {
            return a * b;
        }

        public int Sub(int a, int b)
        {
            return a - b;
        }

        public int Conjugate(int a)
        {
            return a;
        }

        public double SqrtDbl(int a)
        {
            return Math.Sqrt(a);
        }


        public int Sqrt(int a)
        {
            return (int)Math.Sqrt(a);
        }


        public double ToDouble(int a)
        {
            return (double)a;
        }

    }
}
