using System;
using System.Collections.Generic;
using System.IO;
using System.Numerics;
using System.Runtime.InteropServices;
using System.Threading;
using TelmaQuasar.Core;
using TelmaQuasar.Core.MKL;

namespace BlazorServer.Matrix
{
    public class SimpleSparsePortrait : ISparsePortrait
    {
        Dictionary<int, SortedSet<int>> Rows = new Dictionary<int, SortedSet<int>>();
        Dictionary<int, SortedSet<int>> Columns = new Dictionary<int, SortedSet<int>>();
        public int RowSize => Rows.Count;

        public int ColumnSize => Columns.Count;

        public void Add(IEnumerable<int> rows, IEnumerable<int> columns)
        {
            foreach(var row in rows)
            {
                if (Rows.TryGetValue(row, out var Row)) Row.UnionWith(columns);
                else Rows.Add(row, new SortedSet<int>(columns));
            }
            foreach (var column in columns)
            {
                if (Columns.TryGetValue(column, out var Column)) Column.UnionWith(rows);
                else Columns.Add(column, new SortedSet<int>(rows));
            }
        }
        public void Add(ReadOnlySpan<int> rows, ReadOnlySpan<int> columns)
        {
            foreach (var row in rows)
            {
                if (!Rows.TryGetValue(row, out var Row))
                {
                    Row = new SortedSet<int>();
                    Rows.Add(row, Row);
                }
                foreach (var c in columns) Row.Add(c);
            }
            foreach (var column in columns)
            {
                if (!Columns.TryGetValue(column, out var Column))
                {
                    Column = new SortedSet<int>();
                    Columns.Add(column, Column);
                }
                foreach (var r in rows) Column.Add(r);
            }
        }

        public IEnumerable<int> Column(int column) => Columns[column];

        public IEnumerable<int> Row(int row) => Rows[row];
    }
    public static class MatrixAssistant
    {
        public static bool CanWorkAsLUFactorizer(this FactorizerType factorizer)
        {
            switch (factorizer)
            {
                case FactorizerType.LUSQ:return true;
                case FactorizerType.LU:return true;
                case FactorizerType.DI:return true;
                case FactorizerType.NONE:return true;
                default:return false;
            }
        }
        public static bool CanUseUxMult(this FactorizerType factorizer)
        {
            switch (factorizer)
            {
                case FactorizerType.LUSQ: return true;
                case FactorizerType.LU: return true;
                case FactorizerType.DI: return true;
                case FactorizerType.NONE: return true;
                default: return false;
            }
        }
        public static bool CanUse(this SolverType solver,FactorizerType type)
        {
            switch (solver)
            {
                case SolverType.LOS:
                case SolverType.CG:return type.CanWorkAsLUFactorizer();
                case SolverType.GMRES:return type.CanWorkAsLUFactorizer() && type.CanUseUxMult();
                case SolverType.BiCG:return type.CanWorkAsLUFactorizer() && type.CanUseUxMult();
                case SolverType.BiCGStab: return type.CanWorkAsLUFactorizer() && type.CanUseUxMult();
                case SolverType.PARDISO:return false;
                default:return false;
            }
        }
        public static bool IsEnabledFor(this SolverType type, SparseMatrixType matrix)
        {
            switch (matrix)
            {
                case SparseMatrixType.ROWCOL_SYM_SPARSE:
                    switch (type)
                    {
                        case SolverType.CG: return true;
                        case SolverType.LOS: return true;
                        case SolverType.GMRES: return true;
                        case SolverType.BiCG: return false;
                        case SolverType.BiCGStab: return true;
                        case SolverType.PARDISO: return false;
                        default: return false;
                    }
                case SparseMatrixType.ROWCOL_NONSYM_SPARSE:
                    switch (type)
                    {
                        case SolverType.CG: return false;
                        case SolverType.LOS: return true;
                       case SolverType.GMRES: return true;
                        case SolverType.BiCG: return true;
                        case SolverType.BiCGStab: return true;
                        case SolverType.PARDISO: return false;
                        default: return false;
                    }
                case SparseMatrixType.ROW_SPARSE_FOR_PARDISO:
                    switch (type)
                    {
                        case SolverType.CG: return false;
                        case SolverType.LOS: return true;
                        case SolverType.GMRES: return true;
                        case SolverType.BiCG: return true;
                        case SolverType.BiCGStab: return true;
                        case SolverType.PARDISO: return true;
                        default: return false;
                    }
                case SparseMatrixType.ROW_SYM_SPARSE_FOR_PARDISO:
                    switch (type)
                    {
                        case SolverType.CG: return true;
                        case SolverType.LOS: return true;
                        case SolverType.GMRES: return true;
                        case SolverType.BiCG: return false;
                        case SolverType.BiCGStab: return true;
                        case SolverType.PARDISO: return true;
                        default: return false;
                    }
                case SparseMatrixType.ROWCOL_SYM_SPARSE_Harmonic:
                    switch (type)
                    {
                        case SolverType.CG: return true;
                        case SolverType.LOS: return true;
                        case SolverType.GMRES: return true;
                        case SolverType.BiCG: return false;
                        case SolverType.BiCGStab: return true;
                        case SolverType.PARDISO: return false;
                        default: return false;
                    }
                case SparseMatrixType.ROWCOL_NONSYM_SPARSE_Harmonic:
                    switch (type)
                    {
                        case SolverType.CG: return false;
                        case SolverType.LOS: return true;
                        case SolverType.GMRES: return true;
                        case SolverType.BiCG: return true;
                        case SolverType.BiCGStab: return true;
                        case SolverType.PARDISO: return false;
                        default: return false;
                    }
                case SparseMatrixType.ROW_SPARSE_HARMONIC_FOR_PARDISO:
                    switch (type)
                    {
                        case SolverType.CG: return false;
                        case SolverType.LOS: return true;
                        case SolverType.GMRES: return true;
                        case SolverType.BiCG: return true;
                        case SolverType.BiCGStab: return true;
                        case SolverType.PARDISO: return true;
                        default: return false;
                    }

                default:
                    return false;
            }
        }
        public static void CalcResidual<T>(this IRectangularMatrix<T> matrix, ReadOnlySpan<T> pr, ReadOnlySpan<T> u, Span<T> r)where T:struct
        {
            matrix.MultVect(u, r);
            r.Add(pr, -1.0);
            r.Mult(-1.0);
        }

        static void SimpleMult(IRectangularMatrix<double> matr, ReadOnlySpan<double> to, Span<double> res)
        {
            res.Fill(0);
            foreach(var (i,j,m) in matr.AllNonZeroElements) res[i] += to[j] * m;
        }
        static void SimpleTMult(IRectangularMatrix<double> matr, ReadOnlySpan<double> to, Span<double> res)
        {
            res.Fill(0);
            foreach (var (i, j, m) in matr.AllNonZeroElements) res[j] += to[i] * m;
        }
        static void SimpleMult(this IRectangularMatrix<Complex> matr, ReadOnlySpan<Complex> to, Span<Complex> res)
        {
            res.Fill(0);
            foreach (var (i, j, m) in matr.AllNonZeroElements) res[i] += to[j] * m;
        }
        static void SimpleTMult(this IRectangularMatrix<Complex> matr, ReadOnlySpan<Complex> to, Span<Complex> res)
        {
            res.Fill(0);
            foreach (var (i, j, m) in matr.AllNonZeroElements) res[j] += to[i] * m;
        }
        public static void SimpleMultRelease<T>(this IRectangularMatrix<T> matr, ReadOnlySpan<T> to, Span<T> res) where T : struct
        {
            switch (matr)
            {
                case IRectangularMatrix<double> dblmatr:
                    SimpleMult(dblmatr, MemoryMarshal.Cast<T, double>(to), MemoryMarshal.Cast<T, double>(res));
                    return;
                case IRectangularMatrix<Complex> cplxmatr:
                    SimpleMult(cplxmatr, MemoryMarshal.Cast<T, Complex>(to), MemoryMarshal.Cast<T, Complex>(res));
                    return;
                default: throw new ArgumentOutOfRangeException("Only complex and double matrix elements available");
            }
        }
        public static void SimpleTMultRelease<T>(this IRectangularMatrix<T> matr, ReadOnlySpan<T> to, Span<T> res) where T : struct
        {
            switch (matr)
            {
                case IRectangularMatrix<double> dblmatr:
                    SimpleTMult(dblmatr, MemoryMarshal.Cast<T, double>(to), MemoryMarshal.Cast<T, double>(res));
                    return;
                case IRectangularMatrix<Complex> cplxmatr:
                    SimpleTMult(cplxmatr, MemoryMarshal.Cast<T, Complex>(to), MemoryMarshal.Cast<T, Complex>(res));
                    return;
                default: throw new ArgumentOutOfRangeException("Only complex and double matrix elements available");
            }
        }
        public static string DisplayName(this FactorizerType type)
        {
            switch (type)
            {
                case FactorizerType.DI: return "Diagonal";
                case FactorizerType.LU: return "LU";
                case FactorizerType.LUSQ: return "LUsq";
                case FactorizerType.NONE: return "None";
                default: return "unknown";
            }
        }
        public static string DisplayName(this SolverType type)
        {

            switch (type)
            {
                case SolverType.BiCG: return "BiCG";
                case SolverType.BiCGStab: return "BiCGStub";
                case SolverType.CG: return "CG";
                case SolverType.GMRES: return "GMRES";
                case SolverType.LOS: return "LOS";
                case SolverType.PARDISO: return "Pardiso";
                default: return "unknown";
            }
        }
        public static string DisplayName(this IMatrixFactorizerBase fact)
        {
            return DisplayName(fact.DisplayCode);
        }
        public static string DisplayName(this ISlaeSolverBase fact)
        {
            return DisplayName(fact.DisplayCode);
        }
        public static ISparseMatrixBase CreateMatrix(SparseMatrixType type)
        {
            switch (type)
            {
                case SparseMatrixType.ROWCOL_SYM_SPARSE: return new RowSparseSymmetricSlaeOfDouble();
                case SparseMatrixType.ROWCOL_NONSYM_SPARSE: return new RowSparseNonSymSlaeOfDouble();
                // case SparseMatrixType.ROW_SPARSE_FOR_PARDISO: return new CSRPardisoMatrix();

                case SparseMatrixType.ROWCOL_SYM_SPARSE_Harmonic: return new RowSparseSymmetricMatrixOfComplex();
                case SparseMatrixType.ROWCOL_NONSYM_SPARSE_Harmonic: return new RowSparseNonSymSlaeOfComplex();
                // case SparseMatrixType.ROW_SPARSE_HARMONIC_FOR_PARDISO: return new RowOnlySparseNonSymSlaeOfComplex();
                // case SparseMatrixType.ROW_SYM_SPARSE_FOR_PARDISO:return new SymmCSRPardisoMatrix();
                default:
                    throw new ArgumentOutOfRangeException();
            }
        }
        public static IMatrixFactorizerBase CreateFactorizer(FactorizerType code, ISparseMatrixBase sl)
        {
            if (sl is SparseMatrix<double> sld)
            {
                switch (code)
                {
                    case FactorizerType.NONE: return null;
                    case FactorizerType.DI:
                        return new DiagonalFactorizer<double>(sld);
                    case FactorizerType.LUSQ:
                        {
                            if (sl is RowSparseNonSymMatrix<double> sldns)
                                return new LuSQFactorizerForRowColNonSymOfDouble(sldns);
                            if (sl is RowSparseSymmetricMatrix<double> slds) return new LuSQFactorizerForRowColSymOfDouble(slds);
                        }
                        break;

                    case FactorizerType.LU:
                        {
                            if (sl is RowSparseNonSymMatrix<double>sldns) return new LuFactorizerForRowColNonSymOfDouble(sldns);
                            break;
                        }
                    default:
                        break;
                }
            }
            if (sl is SparseMatrix<Complex> slc)
            {
                switch (code)
                {
                    case FactorizerType.NONE: return null;
                    case FactorizerType.DI:
                        return new DiagonalFactorizer<Complex>(slc);
                    case FactorizerType.LUSQ:
                        {
                            if (sl is RowSparseNonSymMatrix<Complex>slcns)
                                return new LuSQFactorizerForRowColNonSymOfComplex(slcns);
                            if (sl is RowSparseSymmetricMatrix<Complex> slcs) return new LuSQFactorizerForRowColSymOfComplex(slcs);
                        }
                        break;

                    case FactorizerType.LU:
                        {
                            if (sl is RowSparseNonSymMatrix<Complex> slcns) return new LuFactorizerForRowColNonSymOfComplex(slcns);
                            break; 
                        }
                    default:
                        break;
                }
            } 

            throw new NotImplementedException();
        }
            public static ISlaeSolver<T> CreateSolver<T>(SolverType code) where T : struct
        {
            switch (code)
            {
                case SolverType.CG: return new CGSolver<T>();
                case SolverType.LOS: return new LOSSolver<T>();
                case SolverType.GMRES: return (ISlaeSolver<T>)new GMRESSolver();
                case SolverType.PARDISO: return (ISlaeSolver<T>)new PardisoSolver();
                case SolverType.BiCGStab: return new BiCGStabSolver<T>();
                default:throw new ArgumentOutOfRangeException();
            }
        }
        public static void ThreadSafeAddDouble(ref double to, double value)
        {
            double s, sum;
            do
            {
                s = to;
                sum = s + value;
            } while (s != Interlocked.CompareExchange(ref to, sum, s) && !double.IsNaN(value));
        }

        public static void ThreadSafeAddSpanElement(Span<double> to, int toind, double value)
        {
            ThreadSafeAddDouble(ref to[toind], value);
        }
        public static void ThreadSafeAddSpanElement(Span<Complex> to, int toind, Complex value)
        {
            var dbl = MemoryMarshal.Cast<Complex,double>(to.Slice(toind)).Slice(0, 2);
            ThreadSafeAddDouble(ref dbl[0], value.Real);
            ThreadSafeAddDouble(ref dbl[1], value.Imaginary);
        }
        public static void AddToSetWithoutCheck(IList<ISet<int>> sets, int ind1, int ind2)
        {
            sets[ind1].Add(ind2);
        }

        public static void AddToSet(IList<ISet<int>> sets, int ind1, int ind2)
        {
            if (ind1 < ind2)
                AddToSetWithoutCheck(sets, ind2, ind1);
            else
                AddToSetWithoutCheck(sets, ind1, ind2);
        }

        public static int QuickFindInOrderedIgJg(ReadOnlySpan<int> ig, ReadOnlySpan<int> jg, int ind1, int ind2)//! Поиск адреса элемента ind1, ind2 в матрице
        {

            return jg.Slice(ig[ind1], ig[ind1 + 1] - ig[ind1]).BinarySearch(ind2) + ig[ind1];
            //	for(int ind=edgeig[ind1];ind<edgeig[ind1+1];ind++)if(edgejg[ind] == ind2)break;
            //        if(jg[ind] != ind2 || ind<ig[ind1 + 1] || ind >= ig[ind1])throw 
        }

        public static double SparseMult(ReadOnlySpan<double> gg, ReadOnlySpan<int> jg, ReadOnlySpan<double> vec)
        {
            double sum = 0;
            for (int i = 0; i < jg.Length; i++) sum += vec[jg[i]] * gg[i];
            return sum;
        }
        public static void SparseAdd(Span<double> w, ReadOnlySpan<int> jg, ReadOnlySpan<double> gg, double val)
        {
            for (int i = 0; i < jg.Length; i++) w[jg[i]] += gg[i] * val;
        }
        //public static void SparseAdd<T>(Span<T> w, ReadOnlySpan<int> jg, ReadOnlySpan<T> gg, T val) where T:struct
        //{
        //    if (val is double dbl) SparseAdd(w.NonPortableCast<T, double>(), jg, gg.NonPortableCast<T,double>(), dbl);
        //    if (val is Complex cplx) SparseAdd(w.NonPortableCast<T, Complex>(), jg, gg.NonPortableCast<T, Complex>(), cplx);
        //    throw new NotImplementedException();
        //}
        public static Complex SparseMult(ReadOnlySpan<Complex> gg, ReadOnlySpan<int> jg, ReadOnlySpan<Complex> vec)
        {
            Complex sum = 0;
            for (int i = 0; i < jg.Length; i++) sum += vec[jg[i]] * gg[i];
            return sum;
        }
        public static void SparseAdd(Span<Complex> w, ReadOnlySpan<int> jg, ReadOnlySpan<Complex> gg, Complex val)
        {
            for (int i = 0; i < jg.Length; i++) w[jg[i]] += gg[i] * val;
        }
        public static void WriteBinary<T>(string path, ReadOnlySpan<T> array) where T : struct
        {
            using (var writer = new System.IO.BinaryWriter(File.OpenWrite(path)))
            {
                var span = MemoryMarshal.AsBytes(array);
                for (int i = 0; i < span.Length; i++) writer.Write(span[i]);
            }
        }
        public static void ReadBinary<T>(string path, Span<T> array) where T : struct
        {
            using (var reader = new System.IO.BinaryReader(File.OpenRead(path)))
            {
                var span = MemoryMarshal.AsBytes(array);
                for (int i = 0; i < span.Length; i++) span[i] = reader.ReadByte();
            }
        }
    }
    /*Класс локальных прямоугольных матриц*/

    public class SimpleRectangularMatrix<T> : IRectangularMatrix<T> where T : struct
    {
        protected int row, col;
        protected T[] matr;
        public SimpleRectangularMatrix(int r = 0, int c = 0)
        {
            row = r;
            col = c == 0 ? r : c;
            matr = new T[row * col];
        }
        public void MultVect(ReadOnlySpan<T> to, Span<T> res)
        {
            this.SimpleMultRelease(to, res);
        }

        public void TMultVect(ReadOnlySpan<T> to, Span<T> res)
        {
            this.SimpleTMultRelease(to, res);
        }
        public IEnumerable<(int, int, T)> AllNonZeroElements
        {
            get
            {
                for (int i = 0; i < row; i++)
                    for (int j = 0; j < col; j++) yield return (i, j, matr[i * col + j]);
            }
        }

        public int CSize => col;

        public int RSize => row;


        public T this[int i, int j]
        {
            get => matr[i * col + j];
            set { matr[i * col + j] = value; }
        }
        public void Nullify()
        {
            for (int i = 0; i < matr.Length; i++) matr[i] = default;
        }

        public void SymmetryScale(ReadOnlySpan<double> Multiplier)
        {
            var arithmetic = TemplateArithmetic.Create<T>();
            for (int i = 0; i < CSize; i++)
                for (int j = 0; j < RSize; j++)
                    this[i, j] = arithmetic.Multiply(this[i, j], Multiplier[i] * Multiplier[j]);
        }
    }
    public class SimpleRectangularMatrixOfDouble:SimpleRectangularMatrix<double>
    {
        public SimpleRectangularMatrixOfDouble(int r = 0, int c = 0) : base(r, c) { }
        public void ThreadSafeAdd(int i,int j,double v)
        {
            MatrixAssistant.ThreadSafeAddSpanElement(matr.AsSpan(), i * col + j, v);
        }
    }
    public class SimpleSquareSymmetryMatrixOfDouble : IInversableMatrix<double>
    {
        int sz;
        double[] matr;
        double[][] lmatr;
        int[] lmatrrownum;

        public int Size => sz;

        public int RSize => sz;

        public int CSize => sz;

        double this[int i, int j] => matr[(i * (i + 1)) / 2 + j];

        void PrepareGauss()
        {
            lmatr = new double[sz][];
            lmatrrownum = new int[sz];
            for (int i = 0; i < sz; i++)
            {
                lmatr[i] = new double[sz];
                for (int j = 0; j <= i; j++)
                {
                    lmatr[j][i] = lmatr[i][j] = this[i, j];
                }
            }
            for (int i = 0; i < sz; i++)
            {
                int maxj = i;
                for (int j = i + 1; j < sz; j++)
                {
                    if (Math.Abs(lmatr[j][i]) > Math.Abs(lmatr[maxj][i])) maxj = j;
                }
                lmatrrownum[i] = maxj;
                if (maxj != i)
                {
                    for (int j = i; j < sz; j++)
                    {
                        var d = lmatr[i][j];
                        lmatr[i][j] = lmatr[maxj][j];
                        lmatr[maxj][j] = d;
                    }
                }
                for (int j = i + 1; j < sz; j++) lmatr[i][j] /= lmatr[i][i];
                for (int j = i + 1; j < sz; j++)
                {
                    for (int k = i + 1; k < sz; k++)
                        lmatr[j][k] -= lmatr[i][k] * lmatr[j][i];
                }
            }
        }
        void ProcessGaussRightPart(Span<double> rp)
        {
            for (int i = 0; i < sz; i++)
            {
                if (lmatrrownum[i] != i)
                {
                    var d = rp[i];
                    rp[i] = rp[lmatrrownum[i]];
                    rp[lmatrrownum[i]] = d;
                }
                if (Math.Abs(lmatr[i][i]) < 1e-30) throw new DivideByZeroException();
                rp[i] /= lmatr[i][i];
                for (int j = i + 1; j < sz; j++)
                {
                    rp[j] -= rp[i] * lmatr[j][i];
                }
            }
        }

        void UMult(Span<double> to) // решение СЛАУ с верхней треугольной матрицей разложения
        {
            for (int i = sz - 1; i >= 0; i--)
            {
                for (int j = 0; j < i; j++) to[j] -= lmatr[j][i] * to[i];
            }
        }

        public void InverseMult(ReadOnlySpan<double> to, Span<double> res)
        {
            if (lmatr == null)
            {
                PrepareGauss();
            }
            for (int i = 0; i < sz; i++) res[i] = to[i];
            ProcessGaussRightPart(res);
            UMult(res);
        }

        public void MultVect(ReadOnlySpan<double> to, Span<double> res)
        {
            this.SimpleMultRelease(to, res);
        }

        public void TMultVect(ReadOnlySpan<double> to, Span<double> res)
        {
            this.SimpleTMultRelease(to, res);
        }

        public IEnumerable<(int, int, double)> AllNonZeroElements
        {
            get
            {
                for (int i = 0; i < sz; i++)
                {
                    var start = (i * (i + 1)) / 2;
                    yield return (i, i, matr[start + i]);
                    for (int j = 0; j < i; j++)
                    {
                        var d = matr[start + j];
                        yield return (i, j, d);
                        yield return (j, i, d);
                    }
                }
            }
        }
        public void SymmetryScale(ReadOnlySpan<double> Multiplier)
        {
            for (int i = 0; i < CSize; i++)
                for (int j = 0; j <= i; j++)
                    matr[(i * (i + 1)) / 2 + j] *= Multiplier[i] * Multiplier[j];
        }

        public SimpleSquareSymmetryMatrixOfDouble(int size = 0) { Resize(size); }
        void Resize(int _size)
        {
            sz = _size;
            matr = new double[(sz * (sz + 1)) / 2];
            lmatr = null;
            lmatrrownum = null;
        }
        public void Set(int i,int j,double val)
        {
            if(j<=i)matr[(i * (i + 1)) / 2 + j]=val;
        }
        public void Add(int i, int j, double val)
        {
            if (j <= i) matr[(i * (i + 1)) / 2 + j] += val;
        }
    }
    public class LocalRectangularMatrixWithZeros<T> : IRectangularMatrix<T> where T : struct
    {
        int row, col;// текущее кол-во столбцов строк
        Dictionary<int, Dictionary<int, T>> NotZero;

        public LocalRectangularMatrixWithZeros(int rows = 0, int columns = 0)
        {
            Resize(rows, columns);
        }
        void Resize(int r, int c)//Установить размеры прямоугольной матрицы
        {
            row = r;
            col = c;
            NotZero = new Dictionary<int, Dictionary<int, T>>();
        }
        public int RSize => row;
        public int CSize => col;

        public void MultVect(ReadOnlySpan<T> to, Span<T> res)
        {
            this.SimpleMultRelease(to, res);
        }
        public void TMultVect(ReadOnlySpan<T> to, Span<T> res)
        {
            this.SimpleTMultRelease(to, res);
        }
        public void SymmetryScale(ReadOnlySpan<double> Multiplier)
        {
            var old = NotZero;
            NotZero = new Dictionary<int, Dictionary<int, T>>();
            var arithmetic = TemplateArithmetic.Create<T>();
                    foreach (var i in old)
                        foreach (var j in i.Value)
                        {
                            this[i.Key, j.Key]= arithmetic.Multiply(j.Value, Multiplier[i.Key] * Multiplier[j.Key]); 
                        }
        }

        public IEnumerable<(int, int, T)> AllNonZeroElements
        {
            get
            {
                foreach (var i in NotZero)
                    foreach (var j in i.Value)
                    {
                        yield return (i.Key, j.Key, j.Value);
                    }
            }
        }
        public T this[int i, int j]
        {
            get
            {
                if (NotZero.TryGetValue(i, out var dict))
                {
                    if (dict.TryGetValue(j, out T val)) return val;
                }
                return default;
            }
            set
            {
                if (NotZero.TryGetValue(i, out var dict))
                    dict[j] = value;
                else
                {
                    var newdict = new Dictionary<int, T>();
                    newdict[j] = value;
                    NotZero[i] = newdict;
                }
            }
        }
    }
    public class IdentityMatrix : IInversableMatrix<double>
    {
        int sz;
        public IdentityMatrix(int size) { sz = size; }
        public double Det => 1;
        public void InverseMult(ReadOnlySpan<double> to, Span<double> res)
        {
            for (int i = 0; i < sz; i++) res[i] = to[i];
        }

        public void MultVect(ReadOnlySpan<double> to, Span<double> res)
        {
            for (int i = 0; i < sz; i++) res[i] = to[i];
        }

        public void TMultVect(ReadOnlySpan<double> to, Span<double> res)
        {
            for (int i = 0; i < sz; i++) res[i] = to[i];
        }

        public void SymmetryScale(ReadOnlySpan<double> Multiplier)
        {
            throw new NotImplementedException();
        }

        public IEnumerable<(int, int, double)> AllNonZeroElements
        {
            get
            {
                for (int i = 0; i < sz; i++) yield return (i, i, 1);
            }
        }
        public int Size => sz;
        public int RSize => sz;
        public int CSize => sz;
    }

  
    public class SpecialPsiPhiTMatrix<T> : IRectangularMatrix<T> where T : struct
    {
        T[] psi; T[] phi;
        public void Set(ReadOnlySpan<T> Psi, ReadOnlySpan<T> Phi)
        {
            psi = new T[Psi.Length];
            phi = new T[Phi.Length];
            psi.AsSpan().Assign(Psi);
            phi.AsSpan().Assign(Phi);
        }
        public void Clear() { psi = null; phi = null; }
        public int CSize => phi == null ? 0 : phi.Length;
        public int RSize => psi == null ? 0 : psi.Length;
        public void MultVect(ReadOnlySpan<T> to, Span<T> res)
        {
            T d = to.Mult(phi);
            res.Assign(psi, d);
        }

        public void TMultVect(ReadOnlySpan<T> to, Span<T> res)
        {
            T d = to.Mult(psi);
            res.Assign(phi, d);
        }
        public void SymmetryScale(ReadOnlySpan<double> Multiplier)
        {
            psi.AsSpan().MultAll(Multiplier);
            phi.AsSpan().MultAll(Multiplier);
        }

        public IEnumerable<(int, int, T)> AllNonZeroElements
        {
            get
            {
                if (psi is double[] dblpsi && phi is double[] dblphi)
                {
                    for (int i = 0; i < dblpsi.Length; i++)
                        for (int j = 0; j < dblphi.Length; j++)
                        {
                            yield return (i, j, (T)(ValueType)(dblpsi[i] * dblphi[j]));
                        }
                }
                else if (psi is Complex[] cplxpsi && phi is Complex[] cplxphi)
                {
                    for (int i = 0; i < cplxpsi.Length; i++)
                        for (int j = 0; j < cplxphi.Length; j++)
                        {
                            yield return (i, j, (T)(ValueType)(cplxpsi[i] * cplxphi[j]));
                        }
                }
                else throw new NotImplementedException();
            }
        }
    }
    public class Matrix2x2 : IInversableMatrix<double>
    {
        double[,] matr = new double[2, 2];
        double[,] invmatr = new double[2, 2];
        double det = 0;
        bool inversed = false;

        void CalcInverse()
        {
            det = matr[1, 1] * matr[0, 0] - matr[0, 1] * matr[1, 0];
            invmatr[0, 0] = matr[1, 1] / det;
            invmatr[0, 1] = -matr[0, 1] / det;
            invmatr[1, 0] = -matr[1, 0] / det;
            invmatr[1, 1] = matr[0, 0] / det;
            inversed = true;
        }

        public int Size => 2;
        public void MultVect(ReadOnlySpan<double> to, Span<double> res)
        {
            res[0] = matr[0, 0] * to[0] + matr[0, 1] * to[1];
            res[1] = matr[1, 0] * to[0] + matr[1, 1] * to[1];
        }

        public void TMultVect(ReadOnlySpan<double> to, Span<double> res)
        {
            res[0] = matr[0, 0] * to[0] + matr[1, 0] * to[1];
            res[1] = matr[0, 1] * to[0] + matr[1, 1] * to[1];
        }

        public IEnumerable<(int, int, double)> AllNonZeroElements
            {
                get
                {
                    yield return (0, 0, matr[0, 0]);
                    yield return (0, 1, matr[0, 1]);
                    yield return (1, 0, matr[1, 0]);
                    yield return (1, 1, matr[1, 1]);
                }
            }
        void SetColumns(ReadOnlySpan<double> p1, ReadOnlySpan<double> p2)
        {
            matr[0, 0] = p1[0]; matr[0, 1] = p2[0]; matr[1, 0] = p1[1]; matr[1, 1] = p2[1];
            inversed = false;
        }
        public double Det
        {
            get
            {
                if (!inversed) CalcInverse();
                return det;
            }
        }

        public int RSize => 2;

        public int CSize => 2;

        double[] Solve(ReadOnlySpan<double> rp)
        {
            if (!inversed) CalcInverse();
            return new double[] { invmatr[0, 0] * rp[0] + invmatr[0, 1] * rp[1], invmatr[1, 0] * rp[0] + invmatr[1, 1] * rp[1] };
        }
        public double this[int i, int j]
        {
            get { return matr[i, j]; }
            set { matr[i, j] = value; }
        }

        public void InverseMult(ReadOnlySpan<double> to, Span<double> res)
        {
            res.Assign(Solve(to));
        }
        public void InverseTMult(ReadOnlySpan<double> to, Span<double> res)
        {
            if (!inversed) CalcInverse();
            res[0] = invmatr[0, 0] * to[0] + invmatr[1, 0] * to[1];
            res[1] = invmatr[0, 1] * to[0] + invmatr[1, 1] * to[1];
        }

        public void SymmetryScale(ReadOnlySpan<double> Multiplier)
        {
            throw new NotImplementedException();
        }
    }
    public class Matrix3x3 : IInversableMatrix<double>
    {
        protected readonly double [] matr = new double[9];
        protected readonly double[] inv = new double[9];
        protected double det = 0;
        protected bool hasinverse = false;
        protected void Inverse()
        {
            inv[0] = matr[4] * matr[8] - matr[7] * matr[5];
            inv[1] = matr[2] * matr[7] - matr[1] * matr[8];
            inv[2] = matr[1] * matr[5] - matr[2] * matr[4];
            inv[3] = matr[6] * matr[5] - matr[3] * matr[8];
            inv[4] = matr[0] * matr[8] - matr[6] * matr[2];
            inv[5] = matr[2] * matr[3] - matr[0] * matr[5];
            inv[6] = matr[3] * matr[7] - matr[6] * matr[4];
            inv[7] = matr[6] * matr[1] - matr[0] * matr[7];
            inv[8] = matr[0] * matr[4] - matr[3] * matr[1];
            //inv[0, 0] = matr[1, 1] * matr[2, 2] - matr[2, 1] * matr[1, 2];
            //inv[0, 1] = matr[0, 2] * matr[2, 1] - matr[0, 1] * matr[2, 2];
            //inv[0, 2] = matr[0, 1] * matr[1, 2] - matr[0, 2] * matr[1, 1];
            //inv[1, 0] = matr[2, 0] * matr[1, 2] - matr[1, 0] * matr[2, 2];
            //inv[1, 1] = matr[0, 0] * matr[2, 2] - matr[2, 0] * matr[0, 2];
            //inv[1, 2] = matr[0, 2] * matr[1, 0] - matr[0, 0] * matr[1, 2];
            //inv[2, 0] = matr[1, 0] * matr[2, 1] - matr[2, 0] * matr[1, 1];
            //inv[2, 1] = matr[2, 0] * matr[0, 1] - matr[0, 0] * matr[2, 1];
            //inv[2, 2] = matr[0, 0] * matr[1, 1] - matr[1, 0] * matr[0, 1];

            det = matr[0] * inv[0] + matr[1] * inv[3] + matr[2] * inv[6];
            hasinverse = true;
        }
        public int Size => 3;
        public void MultVect(ReadOnlySpan<double> to, Span<double> res)
        {
            res[0] = matr[0] * to[0] + matr[1] * to[1] + matr[2] * to[2];
            res[1] = matr[3] * to[0] + matr[4] * to[1] + matr[5] * to[2];
            res[2] = matr[6] * to[0] + matr[7] * to[1] + matr[8] * to[2];
        }
        public void TMultVect(ReadOnlySpan<double> to, Span<double> res)
        {
            res[0] = matr[0] * to[0] + matr[3] * to[1] + matr[6] * to[2];
            res[1] = matr[1] * to[0] + matr[4] * to[1] + matr[7] * to[2];
            res[2] = matr[2] * to[0] + matr[5] * to[1] + matr[7] * to[2];
        }

        public IEnumerable<(int, int, double)> AllNonZeroElements
            {
                get
                {
                    for (int i = 0; i < 3; i++)
                        for (int j = 0; j < 3; j++) yield return (i, j, matr[i*3+ j]);
                }
            }
        public void SetColumns(ReadOnlySpan<double> p1, ReadOnlySpan<double> p2, ReadOnlySpan<double> p3)
        {
            hasinverse = false;
            for (int i = 0; i < 3; i++)
            {
                matr[i*3+ 0] = p1[i];
                matr[i*3+ 1] = p2[i];
                matr[i*3+ 2] = p3[i];
            }
        }
        public void Nullify()
        {
            hasinverse = false;
            for (int i = 0; i < 9; i++) matr[i] = 0;
        }
        public double[] Solve(ReadOnlySpan<double> rp)
        {
            if (!hasinverse) Inverse();
            return new double[] {
                (inv[0] * rp[0] + inv[1] * rp[1] + inv[2] * rp[2])/det,
                (inv[3] * rp[0] + inv[4] * rp[1] + inv[5] * rp[2])/det,
                (inv[6] * rp[0] + inv[7] * rp[1] + inv[8] * rp[2])/det};
        }
        public double this[int i, int j]
        {
            get => matr[i*3+j]; set { matr[i*3+j] = value; }
        }

        public double Det
        {
            get
            {
                if (!hasinverse) Inverse();
                return det;
            }
        }

        public int RSize => 3;

        public int CSize => 3;

        public void InverseMult(ReadOnlySpan<double> rp, Span<double> res)
        {
            if (!hasinverse) Inverse();
            res[0] = (inv[0] * rp[0] + inv[1] * rp[1] + inv[2] * rp[2]) / det;
            res[1] = (inv[3] * rp[0] + inv[4] * rp[1] + inv[5] * rp[2]) / det;
            res[2] = (inv[6] * rp[0] + inv[7] * rp[1] + inv[8] * rp[2]) / det;
        }
        void InverseTMult(ReadOnlySpan<double> to, Span<double> res)
        {
            if (!hasinverse) Inverse();
            res[0] = (inv[ 0] * to[0] + inv[3] * to[1] + inv[6] * to[2]) / det;
            res[1] = (inv[ 1] * to[0] + inv[4] * to[1] + inv[7] * to[2]) / det;
            res[2] = (inv[ 2] * to[0] + inv[5] * to[1] + inv[7] * to[2]) / det;
        }

        public void SymmetryScale(ReadOnlySpan<double> Multiplier)
        {
            throw new NotImplementedException();
        }
    }

}