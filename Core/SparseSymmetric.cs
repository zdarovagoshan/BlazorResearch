using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;
using System.Threading;
namespace Telma.Matrix
{
    public abstract class RowSparseMatrix<T> : SparseMatrix<T> where T : struct
    {
        protected int[] ig;
        protected int[] jg;
        protected T[] di;
        public ReadOnlySpan<int> Ig => ig;
        public ReadOnlySpan<int> Jg => jg;
        public ReadOnlySpan<T> Di => di;
        public override int Size => di == null?0:di.Length;
        public override void SetPortrait(ISparsePortrait portrait)
        {
            var NUnknowns = portrait.RowSize;
            ig = new int[NUnknowns + 1];
            ig[0] = 0;
            for (int i = 0; i < NUnknowns; i++)
            {
                ig[i + 1] = ig[i];
                foreach (var l in portrait.Row(i))
                {
                    if (l < i) ig[i + 1]++;
                }
            }
            jg = new int[ig[NUnknowns]];
            for (int i = 0; i < NUnknowns; i++)
            {
                long jgaddr = ig[i];
                foreach (var j in portrait.Row(i).OrderBy(m => m))
                {
                    if(j<i) jg[jgaddr++] = j;
                }
            }
            AllocateArrays();
        }
        public override T GetDi(int i)
        {
            return di[i];
        }
        protected override void SetDi(int i, T val)
        {
            di[i] = val;
        }

        protected int FindElem(int i, int j)// возвращает индекс в массиве jg
        {
            if (i < j) return MatrixAssistant.QuickFindInOrderedIgJg(ig, jg, j, i);
            else return MatrixAssistant.QuickFindInOrderedIgJg(ig, jg, i, j);
        }
        public int BlockSize { get { return ig.Length - 1; } } //!< Для неблочной матрицы совпадает с обычным размером
    }

    //!симметричная матрица СЛАУ в разрежено-строчном формате
    public abstract class RowSparseSymmetricMatrix<T> : RowSparseMatrix<T> where T : struct
    {
        protected T[] gg;
        public ReadOnlySpan<T> Gg => gg;
        public override void ClearRow(int row)
        {
            throw new NotImplementedException();
        }
        public override bool IsMatrixValid
        {
            get
            {
                return BlockSize > 0 && di.Length == BlockSize &&
                    jg.Length > 0 && gg.Length == jg.Length;
            }
        }
        public override bool IsSymmetry => true;
        public override void PrintFormatted(TextWriter writer, int row = -1)
        {
            throw new NotImplementedException();
        }
        public override IEnumerable<(int, int, T)> AllNonZeroElements
        {
            get
            {
                var sz = BlockSize;
                for (int i = 0; i < sz; i++) yield return (i, i, di[i]);
                for (int i = 0; i < sz; i++)
                {
                    var k1 = ig[i];
                    var k2 = ig[i + 1];
                    for (var j = k1; j < k2; j++)
                    {
                        yield return (i, jg[j], gg[j]);
                        yield return (jg[j], i, gg[j]);
                    }
                }
            }
        }
        public override void Nullify()
        {
            di.AsSpan().Fill(default);
            gg.AsSpan().Fill(default);
            if(Solver!=null)
                SetSolverAndFactorizationType(Solver.DisplayCode, Solver.Factorizer == null ? FactorizerType.NONE : Solver.Factorizer.DisplayCode);
        }
        public override void Assign(ISparseMatrix<T> matr)
        {
            if (matr is RowSparseSymmetricMatrix<T> sl)
            {
                ig = new int[sl.ig.Length];
                ig.AsSpan().Assign(sl.ig);
                jg = new int[sl.jg.Length];
                jg.AsSpan().Assign(sl.jg);
                AllocateArrays();
                di.AsSpan().Assign(sl.di);
                gg.AsSpan().Assign(sl.gg);
                if (Solver != null)
                {
                    var stype = Solver.DisplayCode;
                    var ftype = Solver.Factorizer.DisplayCode;
                    SetSolverAndFactorizationType(stype, ftype);
                    Solver.Dispose();
                }
            }
            else throw new ArgumentOutOfRangeException();
        }
        public override void Read(string path)
        {
            int kuslau = int.Parse(File.ReadAllLines(path + "\\kuslau").FirstOrDefault());
            ig = new int[kuslau + 1];
            MatrixAssistant.ReadBinary(path + "\\ig.dat", ig.AsSpan());
            jg = new int[ig[kuslau]];
            MatrixAssistant.ReadBinary(path + "\\jg.dat", jg.AsSpan());
            AllocateArrays();
            MatrixAssistant.ReadBinary(path + "\\di.dat", di.AsSpan());
            MatrixAssistant.ReadBinary(path + "\\gg.dat", gg.AsSpan());
        }
        public override void Write(string path)
        {
            MatrixAssistant.WriteBinary<int>(path+"\\ig.dat",ig);
            MatrixAssistant.WriteBinary<int>(path+"\\jg.dat", jg);
            MatrixAssistant.WriteBinary<T>(path + "\\di.dat", di);
            MatrixAssistant.WriteBinary<T>(path + "\\gg.dat", gg);
            File.WriteAllText(path + "\\kuslau", $"{BlockSize}\n");
        }
        public override void TMultVect(ReadOnlySpan<T> to, Span<T> res)
        {
            MultVect(to, res);
        }
        protected override void AllocateArrays()
        {
            di = new T[ig.Length - 1];
            gg = new T[jg.Length];
        }
    }
    public class RowSparseSymmetricSlaeOfDouble : RowSparseSymmetricMatrix<double>
    {
        public override void CorrectDiagonal()
        {
            for (int i = 0; i < di.Length; i++) if (di[i] == 0.0) di[i] = 1;
        }
        public override SparseMatrixType MatrixType => SparseMatrixType.ROWCOL_SYM_SPARSE;

        public override int SolverDepth { get => throw new NotImplementedException(); set => throw new NotImplementedException(); }
        override public void AddElement(int i, int j, double d)
        {
            if (i == j)
                MatrixAssistant.ThreadSafeAddDouble(ref di[i], d);
            else if(i>j)
            {
                var k = FindElem(i, j);
                MatrixAssistant.ThreadSafeAddDouble(ref gg[k], d);
            }
        }
        public override void MultVect(ReadOnlySpan<double> to, Span<double> result)
        {
            /*умножение матрицы на вектор*/
            int sz = BlockSize;
            for (int i = 0; i < sz; i++)
            {
                result[i] = di[i] * to[i];
            }
            for (int i=0;i<sz;i++)
            {
                var k1 = ig[i];
                var k2 = ig[i + 1];
                result[i] += MatrixAssistant.SparseMult(gg.AsSpan().Slice(k1, k2 - k1), jg.AsSpan().Slice(k1, k2 - k1), to);
            }
            for (int i = 0; i<sz; i++)
			{
				var k1 = ig[i];
                var k2 = ig[i + 1];
                MatrixAssistant.SparseAdd(result, jg.AsSpan().Slice(k1, k2 - k1), gg.AsSpan().Slice(k1,k2-k1), to[i]);
			}
        }


        protected override void Add(double mult, SparseMatrix<double> matr)
        {
            if (matr is RowSparseSymmetricSlaeOfDouble sl)
            {
                for (int i = 0; i < sl.di.Length; i++) MatrixAssistant.ThreadSafeAddDouble(ref di[i], mult * sl.di[i]);
                for (int i = 0; i < sl.gg.Length; i++) MatrixAssistant.ThreadSafeAddDouble(ref gg[i], mult * sl.gg[i]);
            }
            else throw new ArgumentOutOfRangeException("BadMatrixType");
        }
        public override void Diagonalize()
        {
            var sz = BlockSize;
            if (ig[BlockSize] == 0) return;
            for (int i = 0; i < sz; i++)
            {
                var k1 = ig[i];
                var k2 = ig[i + 1];
                for (var j = k1; j < k2; j++)
                {
                    di[i] += gg[j];
                    di[jg[j]] += gg[j];
                }
            }
            gg = null;
            jg = null;
            for (int i = 0; i < ig.Length; i++) ig[i] = 0;
        }
    }
    public class RowSparseSymmetricMatrixOfComplex : RowSparseSymmetricMatrix<Complex>
    {
        public override void CorrectDiagonal()
        {
            for (int i = 0; i < di.Length; i++) if (di[i] == 0.0) di[i] = 1;
        }
        public override SparseMatrixType MatrixType => SparseMatrixType.ROWCOL_SYM_SPARSE_Harmonic;

        public override int SolverDepth { get => throw new NotImplementedException(); set => throw new NotImplementedException(); }

        override public void AddElement(int i, int j, Complex d)
        {
            if (i == j)
                MatrixAssistant.ThreadSafeAddSpanElement(di,i, d);
            else if(i>j)
            {
                var k = FindElem(i, j);
                MatrixAssistant.ThreadSafeAddSpanElement(gg,k, d);
            }
        }
        public override void MultVect(ReadOnlySpan<Complex> to, Span<Complex> result)
        {
            /*умножение матрицы на вектор*/
            int sz = BlockSize;
            for (int i = 0; i < sz; i++)
            {
                result[i] = di[i] * to[i];
            }
            for (int i = 0; i < sz; i++)
            {
                var k1 = ig[i];
                var k2 = ig[i + 1];
                result[i] += MatrixAssistant.SparseMult(gg.AsSpan().Slice(k1, k2 - k1), jg.AsSpan().Slice(k1, k2 - k1), to);
            }
            for (int i = 0; i < sz; i++)
            {
                var k1 = ig[i];
                var k2 = ig[i + 1];
                MatrixAssistant.SparseAdd(result, jg.AsSpan().Slice(k1, k2 - k1), gg.AsSpan().Slice(k1, k2 - k1), to[i]);
            }
        }


        protected override void Add(double mult, SparseMatrix<Complex> matr)
        {
            if (matr is RowSparseSymmetricMatrixOfComplex sl)
            {
                for (int i = 0; i < sl.di.Length; i++) MatrixAssistant.ThreadSafeAddSpanElement(di,i, mult * sl.di[i]);
                for (int i = 0; i < sl.gg.Length; i++) MatrixAssistant.ThreadSafeAddSpanElement(gg,i, mult * sl.gg[i]);
            }
            else throw new ArgumentOutOfRangeException("BadMatrixType");
        }
        public override void Diagonalize()
        {
            var sz = BlockSize;
            if (ig[BlockSize] == 0) return;
            for (int i = 0; i < sz; i++)
            {
                var k1 = ig[i];
                var k2 = ig[i + 1];
                for (var j = k1; j < k2; j++)
                {
                    di[i] += gg[j];
                    di[jg[j]] += gg[j];
                }
            }
            gg = null;
            jg = null;
            for (int i = 0; i < ig.Length; i++) ig[i] = 0;
        }
    }
}