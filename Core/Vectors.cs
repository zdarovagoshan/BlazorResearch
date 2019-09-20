using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Runtime.CompilerServices;

namespace BlazorServer.Core
{
    public struct Vector2D : IEquatable<Vector2D>
    {
        public double X { get; }
        public double Y { get; }
        public double[] AsArray() => new[] { X, Y };
        public Vector2D(double x, double y) : this()
        {
            X = x;
            Y = y;
        }
        public double this[int i]
        {
            get
            {
                return i switch
                {
                    0 => X,
                    1 => Y,
                    _ => throw new ArgumentOutOfRangeException(),
                };
            }
        }
        public double Norm => Math.Sqrt(X * X + Y * Y);
        public static double Distance(Vector2D a, Vector2D b) => (a - b).Norm;
        public static double SqrDistance(Vector2D a, Vector2D b)
        {
            Vector2D diff = a - b;
            return diff * diff;
        }
        public double Distance(Vector2D b) => Distance(this, b);
        public Vector2D Normalize() => this / Norm;
        public override string ToString() => $"Vec({X}, {Y})";

        public bool Equals([AllowNull] Vector2D a) => a.X == X && a.Y == Y;
        public override bool Equals(object obj) => obj is Vector2D && Equals((Vector2D)obj);
        public override int GetHashCode() => X.GetHashCode() ^ Y.GetHashCode();

        #region Static operators
        public static Vector2D operator -(Vector2D a) => new Vector2D(-a.X, -a.Y);
        public static Vector2D operator +(Vector2D a, Vector2D b) => new Vector2D(a.X + b.X, a.Y + b.Y);
        public static Vector2D operator -(Vector2D a, Vector2D b) => new Vector2D(a.X - b.X, a.Y - b.Y);
        public static Vector2D operator /(Vector2D a, double b) => new Vector2D(a.X / b, a.Y / b);
        public static Vector2D operator *(Vector2D a, double b) => new Vector2D(a.X * b, a.Y * b);
        public static Vector2D operator *(double b, Vector2D a) => new Vector2D(b * a.X, b * a.Y);
        public static double operator *(Vector2D a, Vector2D b) => a.X * b.X + a.Y * b.Y;
        public static bool operator ==(Vector2D a, Vector2D b) => a.X == b.X && a.Y == b.Y;
        public static bool operator !=(Vector2D a, Vector2D b) => a.X != b.X || a.Y != b.Y;
        public static Vector2D Cross(Vector2D a) => new Vector2D(a.Y, -a.X);
        public static double Mixed(Vector2D a, Vector2D b) => a.Y * a.X - a.X * b.Y;
        public static Vector2D Sum(Vector2D a, Vector2D b) => new Vector2D(a.X + b.X, a.Y + b.Y);

        #endregion
    }
    public struct Vector3D : IEquatable<Vector3D>
    {
        public static readonly Vector3D Zero = new Vector3D(0, 0, 0);
        public static readonly Vector3D XAxis = new Vector3D(1, 0, 0);
        public static readonly Vector3D YAxis = new Vector3D(0, 1, 0);
        public static readonly Vector3D ZAxis = new Vector3D(0, 0, 1);
        public static readonly Vector3D[] Axes = { XAxis, YAxis, ZAxis };

        public double X { get; }
        public double Y { get; }
        public double Z { get; }

        public Vector3D(double x, double y, double z)
        {
            X = x;
            Y = y;
            Z = z;
        }
        public Vector3D(Vector2D vec, double z)
        {
            X = vec.X;
            Y = vec.Y;
            Z = z;
        }
        public Vector3D(ReadOnlySpan<double> arr)
        {
#if DEBUG
            if (arr.Length != 3) throw new ArgumentException();
#endif
            X = arr[0];
            Y = arr[1];
            Z = arr[2];
        }
        public double Distance(Vector3D b) => Distance(this, b);

        public double SqrDistance(Vector3D b) => SqrDistance(this, b);
        public static double Distance(Vector3D a, Vector3D b) => (a - b).Norm;

        public static double SqrDistance(Vector3D a, Vector3D b)
        {
            var diff = a - b;
            return diff * diff;
        }
        public double this[int i]
        {
            get
            {
                switch (i)
                {
                    case 0: return X;
                    case 1: return Y;
                    case 2: return Z;
                    default: throw new IndexOutOfRangeException();
                }
            }
        }


        public double[] AsArray() => new[] { X, Y, Z };
        public Vector2D As2D() => new Vector2D(X, Y);

        public double Norm => Math.Sqrt(X * X + Y * Y + Z * Z);

        public double MaxNorm => Math.Max(Math.Abs(X), Math.Max(Math.Abs(Y), Math.Abs(Z)));

        public Vector3D Projection(Vector3D p) => (this * p) * p;

        public Vector3D Normalize() => this / Norm;

        public Vector3D Round(int digits) => new Vector3D(Math.Round(X, digits), Math.Round(Y, digits), Math.Round(Z, digits));

        public override string ToString() => $"Vec({X}, {Y}, {Z})";

        public override bool Equals(object obj) => obj is Vector3D && Equals((Vector3D)obj);

        public override int GetHashCode() => X.GetHashCode() ^ Y.GetHashCode() ^ Z.GetHashCode();

        public bool Equals(Vector3D a) => a.X == X && a.Y == Y && a.Z == Z;

        public static bool TryParse(string line, out Vector3D res)
        {
            double x, y, z;
            var words = line.Split(new[] { ' ', '\t', ',', '(', ')', '<', '>' }, StringSplitOptions.RemoveEmptyEntries);
            if (words[0] == "Vec")
            {
                if (words.Length != 4 || !double.TryParse(words[1], out x) || !double.TryParse(words[2], out y)
                    || !double.TryParse(words[3], out z))
                {
                    res = Vector3D.Zero;
                    return false;
                }
                res = new Vector3D(x, y, z);
                return true;
            }
            if (words.Length != 3 || !double.TryParse(words[0], out x) || !double.TryParse(words[1], out y)
                || !double.TryParse(words[2], out z))
            {
                res = Vector3D.Zero;
                return false;
            }

            res = new Vector3D(x, y, z);
            return true;
        }
        public static Vector3D Vec(double x, double y, double z) => new Vector3D(x, y, z);

        public static Vector3D Parse(string line)
        {
            Vector3D res;
            if (!TryParse(line, out res))
                throw new FormatException("Can't parse Vector3D!");
            return res;
        }

        #region Static operators

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Vector3D operator -(Vector3D a) => new Vector3D(-a.X, -a.Y, -a.Z);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Vector3D operator +(Vector3D a) => a;

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double operator *(Vector3D a, Vector3D b) => a.X * b.X + a.Y * b.Y + a.Z * b.Z;

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Vector3D operator *(double a, Vector3D b) => new Vector3D(a * b.X, a * b.Y, a * b.Z);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Vector3D operator *(Vector3D b, double a) => new Vector3D(a * b.X, a * b.Y, a * b.Z);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Vector3D operator /(Vector3D a, double v) => new Vector3D(a.X / v, a.Y / v, a.Z / v);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Vector3D operator +(Vector3D a, Vector3D b) => new Vector3D(a.X + b.X, a.Y + b.Y, a.Z + b.Z);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Vector3D operator -(Vector3D a, Vector3D b) => new Vector3D(a.X - b.X, a.Y - b.Y, a.Z - b.Z);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static bool operator ==(Vector3D a, Vector3D b) => a.X == b.X && a.Y == b.Y && a.Z == b.Z;

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static bool operator !=(Vector3D a, Vector3D b) => a.X != b.X || a.Y != b.Y || a.Z != b.Z;

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Vector3D Cross(Vector3D v1, Vector3D v2) =>
            new Vector3D(v1.Y * v2.Z - v2.Y * v1.Z, v1.Z * v2.X - v1.X * v2.Z, v1.X * v2.Y - v1.Y * v2.X);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double Mixed(Vector3D v1, Vector3D v2, Vector3D v3) =>
            (v1.Y * v2.Z - v2.Y * v1.Z) * v3.X + (v1.Z * v2.X - v1.X * v2.Z) * v3.Y + (v1.X * v2.Y - v1.Y * v2.X) * v3.Z;

        public static Vector3D Sum(Vector3D a, Vector3D b) => a + b;

        public static Vector3D Min(Vector3D a, Vector3D b) =>
            new Vector3D(Math.Min(a.X, b.X), Math.Min(a.Y, b.Y), Math.Min(a.Z, b.Z));

        public static Vector3D Max(Vector3D a, Vector3D b) =>
            new Vector3D(Math.Max(a.X, b.X), Math.Max(a.Y, b.Y), Math.Max(a.Z, b.Z));

        #endregion

        #region EqualityComparer

        private class EqualityComparer : IEqualityComparer<Vector3D>
        {
            public int Digits { get; set; }

            public bool Equals(Vector3D v1, Vector3D v2)
            {
                return v1.Round(Digits) == v2.Round(Digits);
            }

            public int GetHashCode(Vector3D obj)
            {
                return obj.Round(Digits).GetHashCode();
            }
        }

        public static IEqualityComparer<Vector3D> CreateComparer(int digits = 7)
        {
            return new EqualityComparer { Digits = digits };
        }

        internal object GetCoords()
        {
            throw new NotImplementedException();
        }

        #endregion
    }
}