using System;
using System.Collections.Generic;
using System.Numerics;
using System.Text;
using System.Runtime.InteropServices;
namespace BlazorServer
{
    public static class SpanExtension
    {
        public static void Assign<T, T2>(this Span<T> res, ReadOnlySpan<T> from, T2 coef) where T : struct where T2 : struct
        {
            if (typeof(T) == typeof(double))
            {
                if (coef is double dblc)
                {
                    var dblfrom = MemoryMarshal.Cast<T, double>(from);
                    var dblres = MemoryMarshal.Cast<T, double>(res);
                    for (int i = 0; i < dblfrom.Length; i++) dblres[i] = dblfrom[i] * dblc;
                    return;
                }
            }
            if (typeof(T) == typeof(Complex))
            {
                var cplxfrom = MemoryMarshal.Cast<T, Complex>(from);
                var cplxres = MemoryMarshal.Cast<T, Complex>(res);
                if (coef is double dblcoef)
                {
                    for (int i = 0; i < cplxfrom.Length; i++) cplxres[i] = cplxfrom[i] * dblcoef;
                    return;
                }
                if (coef is Complex cplxcoef)
                {
                    for (int i = 0; i < cplxfrom.Length; i++) cplxres[i] = cplxfrom[i] * cplxcoef;
                    return;
                }
            }
            throw new ArgumentOutOfRangeException("Only complex and double values available");
        }
        public static void Mult<T>(this Span<T> vec, double val) where T : struct
        {
            if (typeof(T) == typeof(double))
            {
                var dblvec = MemoryMarshal.Cast<T, double>(vec);
                for (int i = 0; i < dblvec.Length; i++) dblvec[i] *= val;
                return;
            }
            if (typeof(T) == typeof(Complex))
            {
                var cplxvec = MemoryMarshal.Cast<T, Complex>(vec);
                for (int i = 0; i < cplxvec.Length; i++) cplxvec[i] *=val;
                return;
            }
            throw new ArgumentOutOfRangeException("Only complex and double elements available");
        }

        public static T Mult<T>(this Span<T> vec1, ReadOnlySpan<T> vec2) where T : struct => Mult((ReadOnlySpan<T>)vec1, vec2);
        public static T Mult<T>(this ReadOnlySpan<T> vec1, ReadOnlySpan<T> vec2) where T : struct
        {
            if (typeof(T) == typeof(double))
            {
                var dblvec1 = MemoryMarshal.Cast<T, double>(vec1);
                var dblvec2 = MemoryMarshal.Cast<T, double>(vec2);
                double sum = 0;
                for (int i = 0; i < dblvec1.Length; i++) sum += dblvec1[i] * dblvec2[i];
                return (T)(ValueType)sum;
            }
            if (typeof(T) == typeof(Complex))
            {
                var cplxvec1 = MemoryMarshal.Cast<T, Complex>(vec1);
                var cplxvec2 = MemoryMarshal.Cast<T, Complex>(vec2);
                Complex res = 0;
                for (int i = 0; i < cplxvec1.Length; i++) res += cplxvec1[i] * cplxvec2[i];
                return (T)(ValueType)res;
            }
            throw new ArgumentOutOfRangeException("Only complex and double elements available");
        }
        public static void MultAll<T>(this Span<T> vec1, ReadOnlySpan<T> vec2) where T : struct
        {
            if (typeof(T) == typeof(double))
            {
                var dblvec1 = MemoryMarshal.Cast<T, double>(vec1);
                var dblvec2 = MemoryMarshal.Cast<T, double>(vec2);
                for (int i = 0; i < dblvec1.Length; i++) dblvec1[i] *= dblvec2[i];
                return;
            }
            if (typeof(T) == typeof(Complex))
            {
                var cplxvec1 = MemoryMarshal.Cast<T, Complex>(vec1);
                var cplxvec2 = MemoryMarshal.Cast<T, Complex>(vec2);
                for (int i = 0; i < cplxvec1.Length; i++) cplxvec1[i] *= cplxvec2[i];
                return;
            }
            throw new ArgumentOutOfRangeException("Only complex and double elements available");
        }
        public static void MultAll<T>(this Span<T> vec1, ReadOnlySpan<double> vec2) where T : struct
        {
            if (typeof(T) == typeof(double))
            {
                var dblvec1 = MemoryMarshal.Cast<T, double>(vec1);
                for (int i = 0; i < dblvec1.Length; i++) dblvec1[i] *= vec2[i];
                return;
            }
            if (typeof(T) == typeof(Complex))
            {
                var cplxvec1 = MemoryMarshal.Cast<T, Complex>(vec1);
                for (int i = 0; i < cplxvec1.Length; i++) cplxvec1[i] *= vec2[i];
                return;
            }
            throw new ArgumentOutOfRangeException("Only complex and double elements available");
        }
        public static void DevideAll<T>(this Span<T> vec1, ReadOnlySpan<T> vec2) where T : struct
        {
            if (typeof(T) == typeof(double))
            {
                var dblvec1 = MemoryMarshal.Cast<T, double>(vec1);
                var dblvec2 = MemoryMarshal.Cast<T, double>(vec2);
                for (int i = 0; i < dblvec1.Length; i++) dblvec1[i] /= dblvec2[i];
                return;
            }
            if (typeof(T) == typeof(Complex))
            {
                var cplxvec1 = MemoryMarshal.Cast<T, Complex>(vec1);
                var cplxvec2 = MemoryMarshal.Cast<T, Complex>(vec2);
                for (int i = 0; i < cplxvec1.Length; i++) cplxvec1[i] /= cplxvec2[i];
                return;
            }
            throw new ArgumentOutOfRangeException("Only complex and double elements available");
        }
        public static void DevideAll<T>(this Span<T> vec1, T value) where T : struct
        {
            if (value is double dval)
            {
                var dblvec1 = MemoryMarshal.Cast<T, double>(vec1);
                for (int i = 0; i < dblvec1.Length; i++) dblvec1[i] /= dval;
                return;
            }
            if (value is Complex cval)
            {
                var cplxvec1 = MemoryMarshal.Cast<T, Complex>(vec1);
                for (int i = 0; i < cplxvec1.Length; i++) cplxvec1[i] /= cval;
                return;
            }
            throw new ArgumentOutOfRangeException("Only complex and double elements available");
        }
        public static T Dot<T>(this Span<T> vec1, ReadOnlySpan<T> vec2) where T : struct => Dot((ReadOnlySpan<T>)(vec1), vec2);
        public static T Dot<T>(this ReadOnlySpan<T> vec1, ReadOnlySpan<T> vec2) where T : struct
        {
            if (typeof(T) == typeof(double))
            {
                var dblvec1 = MemoryMarshal.Cast<T, double>(vec1);
                var dblvec2 = MemoryMarshal.Cast<T, double>(vec2);
                double sum = 0;
                for (int i = 0; i < dblvec1.Length; i++) sum += dblvec1[i] * dblvec2[i];
                return (T)(ValueType)sum;
            }
            if (typeof(T) == typeof(Complex))
            {
                var cplxvec1 = MemoryMarshal.Cast<T, Complex>(vec1);
                var cplxvec2 = MemoryMarshal.Cast<T, Complex>(vec2);
                Complex res = 0;
                for (int i = 0; i < cplxvec1.Length; i++)
                {
                    Complex cmp = new Complex(cplxvec2[i].Real, -cplxvec2[i].Imaginary);
                    res += cplxvec1[i] * cmp;
                }
                return (T)(ValueType)res;
            }
            throw new ArgumentOutOfRangeException("Only complex and double elements available");
        }
        public static void Assign<T>(this Span<T> res, ReadOnlySpan<T> from) where T : struct
        {
            for (int i = 0; i < from.Length; i++) res[i] = from[i];
        }
        public static void Add<T>(this Span<T> res, ReadOnlySpan<T> from, double coef) where T : struct
        {
            if (typeof(T) == typeof(double))
            {
                var dblfrom = MemoryMarshal.Cast<T, double>(from);
                var dblres  = MemoryMarshal.Cast<T, double>(res);
                for (int i = 0; i < dblres.Length; i++) dblres[i] += dblfrom[i] * coef;
                return;
            }
            if (typeof(T) == typeof(Complex))
            {
                var cplxfrom = MemoryMarshal.Cast<T, Complex>(from);
                var cplxres  = MemoryMarshal.Cast<T, Complex>(res);
                for (int i = 0; i < cplxres.Length; i++) cplxres[i] += cplxfrom[i] * coef;
                return;
            }
            throw new ArgumentOutOfRangeException("Only complex and double elements available");
        }
        public static void Add<T>(this Span<T> res, ReadOnlySpan<T> from, T coef) where T : struct
        {
            if (typeof(T) == typeof(double) && coef is double coeff)  
            {
                var dblfrom = MemoryMarshal.Cast<T, double>(from);
                var dblres  = MemoryMarshal.Cast<T, double>(res);
                for (int i = 0; i < dblres.Length; i++) dblres[i] += dblfrom[i] * coeff;
                return;
            }
            if (typeof(T) == typeof(Complex) && coef is Complex coefcplx)
            {
                var cplxfrom = MemoryMarshal.Cast<T, Complex>(from);
                var cplxres  = MemoryMarshal.Cast<T, Complex>(res);
                for (int i = 0; i < cplxres.Length; i++) cplxres[i] += cplxfrom[i] * coefcplx;
                return;
            }
            throw new ArgumentOutOfRangeException("Only complex and double elements available");
        }
        public static void Add<T>(this Span<T> res, ReadOnlySpan<T> from) where T : struct
        {
            if (typeof(T) == typeof(double))
            {
                var dblfrom = MemoryMarshal.Cast<T, double>(from);
                var dblres = MemoryMarshal.Cast<T, double>(res);
                for (int i = 0; i < dblres.Length; i++) dblres[i] += dblfrom[i];
                return;
            }
            if (typeof(T) == typeof(Complex))
            {
                var cplxfrom = MemoryMarshal.Cast<T, Complex>(from);
                var cplxres  = MemoryMarshal.Cast<T, Complex>(res);
                for (int i = 0; i < cplxres.Length; i++) cplxres[i] += cplxfrom[i];
                return;
            }
            throw new ArgumentOutOfRangeException("Only complex and double elements available");
        }
        public static void SqrtAll<T>(this Span<T> vec) where T : struct
        {
            if (typeof(T) == typeof(double))
            {
                var dbl = MemoryMarshal.Cast<T, double>(vec);
                for (int i = 0; i < dbl.Length; i++) dbl[i]=Math.Sqrt(dbl[i]);
                return;
            }
            if (typeof(T) == typeof(Complex))
            {
                var cplx = MemoryMarshal.Cast<T, Complex>(vec);
                for (int i = 0; i < cplx.Length; i++)
                {
                    cplx[i] = Complex.Sqrt(cplx[i]);
                }
                return;
            }
            throw new ArgumentOutOfRangeException("Only complex and double elements available");
        }
        /// <summary>
        /// Квадрат нормы
        /// </summary>
        public static double NormSqr<T>(this Span<T> vec) where T : struct => NormSqr((ReadOnlySpan<T>)vec);

        public static double NormSqr<T>(this ReadOnlySpan<T> vec) where T : struct
        {
            if (typeof(T) == typeof(double))
            {
                var dbl = MemoryMarshal.Cast<T, double>(vec);
                double res = 0;
                for (int i = 0; i < dbl.Length; i++) res += dbl[i] * dbl[i];
                return res;
            }
            if (typeof(T) == typeof(Complex))
            {
                var cplx = MemoryMarshal.Cast<T, Complex>(vec);
                double s = 0;
                for (int i = 0; i < cplx.Length; i++)
                {
                    double re = cplx[i].Real;
                    double im = cplx[i].Imaginary;
                    s += re * re + im * im;
                }
                return s;
            }
            throw new ArgumentOutOfRangeException("Only complex and double elements available");
        }
        public static double Norm<T>(this Span<T> vec) where T : struct => Norm((ReadOnlySpan<T>)vec);
        public static double Norm<T>(this ReadOnlySpan<T> vec) where T : struct
        {
            return Math.Sqrt(vec.NormSqr());
        }

    }
}
