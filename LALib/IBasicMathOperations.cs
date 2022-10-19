//
// IBasicMathOperations.cs
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

    public interface IBasicMathOperations<T> //where T : IComparable
    {
        //bool IsEqual(T a, T b);
        bool IsLesserThan(T a, T b);
        bool IsGreaterThan(T a, T b);
        T Minus (T value);
        T Abs(T value); // same as Magnitude for complex numbers
        T Zero { get; }
        T One { get; }
        T Div(T a, T b);
        T Mult(T a, T b);
        T Sub(T a, T b);
        T Add(T a, T b);
        T Conjugate(T a);
        double SqrtDbl(T a);
        T Sqrt(T a);
        T Round(T a, int numDec);

        //double Magnitude(T a); // same as Abs
        double ToDouble(T a);

    }

}

