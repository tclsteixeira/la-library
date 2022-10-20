/******************************************************************************
 * ComplexMathOperations.cs
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
 *
 ******************************************************************************/

using System;

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
    public class ComplexMathOperations : IBasicMathOperations<cpx>
    {
        public ComplexMathOperations()
        {
        }


        public cpx Round(cpx a, int numDec)
        {
            return new cpx( Math.Round(a.Real, numDec), Math.Round(a.Imaginary, numDec) );
        }

        public cpx Minus(cpx a)
        {
            return -a;
        }

        public bool IsLesserThan(cpx a, cpx b)
        {
            return (a.Magnitude < b.Magnitude);
        }

        public bool IsGreaterThan(cpx a, cpx b)
        {
            return (a.Magnitude > b.Magnitude);
        }

        public cpx Abs(cpx value)
        {
            return cpx.Abs(value);
        }

        public cpx Zero
        {
            get { return cpx.Zero; }
        }

        public cpx One
        {
            get { return cpx.One; }
        }

        public cpx Add(cpx a, cpx b)
        {
            return a + b;
        }

        public cpx Div(cpx a, cpx b)
        {
            return a / b;
        }

        public cpx Mult(cpx a, cpx b)
        {
            return a * b;
        }

        public cpx Sub(cpx a, cpx b)
        {
            return a - b;
        }

        public cpx Conjugate(cpx a)
        {
            return a;
        }

        public double SqrtDbl(cpx a)
        {
            throw new NotSupportedException("'cpx SqrtDbl(cpx a)' is not supported in current context.");
        }

        public cpx Sqrt(cpx a)
        {
            return cpx.Sqrt(a);
        }

        public double ToDouble(cpx a)
        {
            return a.Real;
        }

    }
}
