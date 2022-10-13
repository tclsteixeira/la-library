using System;
#if SYSNUMERICS
using cpx = System.Numerics.Complex;
#else
using ComplexN; // use another complex number library
#endif


namespace LALib
{

    /// <summary>
    /// Provides constants and auxiliary functions for using
    /// with matrix compurations.
    /// </summary>
    public static class Utils
    {

        // C dbl_epsilon
        public const double DBL_EPSILON = 2.2204460492503131e-16;
        public const double DBL_TOLERANCE = 1e-15;


        /// <summary>
        /// Checks if two floating point values <paramref name="a"/> and <paramref name="b"/> are almost equals.
        /// </summary>
        /// <returns>Returns <c>true</c>, if two floating point values are almost equals for the given tolerance <paramref name="tol"/>, <c>false</c> otherwise.</returns>
        /// <param name="a">The first value.</param>
        /// <param name="b">The second value.</param>
        /// <param name="tol">The tolerance for comparison (C DBL_EPSILON by default).</param>
        public static bool AlmostEqual(double a, double b, double tol = DBL_EPSILON)
        {
            if (Math.Abs(a - b) < tol)
                return true;
            else
                return false;
        }

    }
}
