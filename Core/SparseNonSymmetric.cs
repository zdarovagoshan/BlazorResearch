using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;
using System.Threading;

namespace BlazorServer.Matrix
{
    public abstract class RowSparseNonSymMatrix<T> : RowSparseMatrix<T> where T : struct
    {

        protected T[] ggl;
        protected T[] ggu;
        public ReadOnlyMemory<T> Ggl => ggl;
        public ReadOnlyMemory<T> Ggu => ggu;
        public override bool IsMatrixValid
        {
            get
            {
                var sz = BlockSize;
                return sz > 0 && di.Length == sz && ggu.Length == jg.Length && ggl.Length == jg.Length && ig.Length == sz + 1 && jg.Length == ig[sz];
            }
        }
        protected override void AllocateArrays()
        {
            var sz = BlockSize;
            di = new T[sz];
            ggl = new T[jg.Length];
            ggu = new T[jg.Length];
        }
        public override void Nullify()
        {
            di.AsSpan().Fill(default);
            ggl.AsSpan().Fill(default);
            ggu.AsSpan().Fill(default);
            SetSolverAndFactorizationType(Solver.DisplayCode, Solver.Factorizer == null?FactorizerType.NONE:Solver.Factorizer.DisplayCode);
        }
        public override bool IsSymmetry => false;
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
                    for (int j = k1; j < k2; j++)
                    {
                        yield return (i, jg[j], ggl[j]);
                        yield return (jg[j], i, ggu[j]);

                    }
                }
            }
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
            MatrixAssistant.ReadBinary(path + "\\ggl.dat", ggl.AsSpan());
            MatrixAssistant.ReadBinary(path + "\\ggu.dat", ggu.AsSpan());
        }
        public override void Write(string path)
        {
            MatrixAssistant.WriteBinary<int>(path + "\\ig.dat", ig.AsSpan());
            MatrixAssistant.WriteBinary<int>(path + "\\jg.dat", jg.AsSpan());
            MatrixAssistant.WriteBinary<T>(path + "\\di.dat", di.AsSpan());
            MatrixAssistant.WriteBinary<T>(path + "\\ggl.dat", ggl.AsSpan());
            MatrixAssistant.WriteBinary<T>(path + "\\ggu.dat", ggu.AsSpan());
            File.WriteAllText(path + "\\kuslau", $"{BlockSize}\n");
        }
        public override void ClearRow(int trow)
        {
            var sz = BlockSize;
            for (var i = trow + 1; i < sz; i++)
                for (var j = ig[i]; j < ig[i + 1]; j++)
                    if (jg[j] == trow) ggu[j] = default;
            for (var i = ig[trow]; i < ig[trow + 1]; i++) ggl[i] = default;
            di[trow] = default;
        }
        public override void PrintFormatted(TextWriter writer, int row = -1)
        {
            throw new NotImplementedException();
        }
        public override void Assign(ISparseMatrix<T> matr)
        {
            if (matr is RowSparseNonSymMatrix<T> sl)
            {
                ig = new int[sl.ig.Length];
                ig.AsSpan().Assign(sl.ig);
                jg = new int[sl.jg.Length];
                jg.AsSpan().Assign(sl.jg);
                AllocateArrays();
                di.AsSpan().Assign(sl.di);
                ggl.AsSpan().Assign(sl.ggl);
                ggu.AsSpan().Assign(sl.ggu);
            }
            else if (matr is RowSparseSymmetricMatrix<T> sls)
            {
                ig = new int[sls.BlockSize + 1];
                ig.AsSpan().Assign(sls.Ig);
                jg = new int[sls.Jg.Length];
                jg.AsSpan().Assign(sls.Jg);
                AllocateArrays();
                di.AsSpan().Assign(sls.Di);
                ggl.AsSpan().Assign(sls.Gg);
                ggu.AsSpan().Assign(sls.Gg);
            }
            else throw new ArgumentOutOfRangeException();
        }
    }
    public class RowSparseNonSymSlaeOfDouble : RowSparseNonSymMatrix<double>
    {
        public override void AddLocalMatrix(ReadOnlySpan<int> indrows, ReadOnlySpan<int> indcolumns, IRectangularMatrix<double> matr)
        {
            throw new NotImplementedException();
        }
        public override void AddElement(int tdi, int tdj, double val)
        {
            if (tdi == tdj)
                MatrixAssistant.ThreadSafeAddDouble(ref di[tdi], val);
            else
            {
                var addr = FindElem(tdi, tdj);
                if (tdi < tdj)
                {
                    MatrixAssistant.ThreadSafeAddDouble(ref ggu[addr], val);
                }
                else
                {
                    MatrixAssistant.ThreadSafeAddDouble(ref ggl[addr], val);
                }
            }
        }
        void AddElem(int i, int j, double dIJ, double dJI)
        {
            if (i == j)
            {
                MatrixAssistant.ThreadSafeAddDouble(ref di[i], dIJ);
            }
            else
            {
                var addr = FindElem(i, j);
                if (i < j)
                {
                    MatrixAssistant.ThreadSafeAddDouble(ref ggu[addr], dIJ);
                    MatrixAssistant.ThreadSafeAddDouble(ref ggl[addr], dJI);
                }
                else
                {
                    MatrixAssistant.ThreadSafeAddDouble(ref ggl[addr], dIJ);
                    MatrixAssistant.ThreadSafeAddDouble(ref ggu[addr], dJI);
                }
            }
        }
        void AddDiagElem(int i, double d)
        {
            MatrixAssistant.ThreadSafeAddDouble(ref di[i], d);
        }
        public override void MultVect(ReadOnlySpan<double> vec, Span<double> res)
        {
            var sz = BlockSize;
            for (int i = 0; i < sz; i++) res[i] = di[i] * vec[i];
            for (int i = 0; i < sz; i++)
            {
                var k1 = ig[i];
                var k2 = ig[i + 1];
                res[i] += MatrixAssistant.SparseMult(ggl.AsSpan().Slice(k1, k2 - k1), jg.AsSpan().Slice(k1, k2 - k1), vec);
            }
            for (int i = 0; i < sz; i++)
            {
                var k1 = ig[i];
                var k2 = ig[i + 1];
                MatrixAssistant.SparseAdd(res, jg.AsSpan().Slice(k1, k2 - k1), ggu.AsSpan().Slice(k1, k2 - k1), vec[i]);
            }
        }
        public override void TMultVect(ReadOnlySpan<double> vec, Span<double> res)
        {
            var sz = BlockSize;
            for (int i = 0; i < sz; i++) res[i] = di[i] * vec[i];
            for (int i = 0; i < sz; i++)
            {
                var k1 = ig[i];
                var k2 = ig[i + 1];
                res[i] += MatrixAssistant.SparseMult(ggu.AsSpan().Slice(k1, k2 - k1), jg.AsSpan().Slice(k1, k2 - k1), vec);
            }
            for (int i = 0; i < sz; i++)
            {
                var k1 = ig[i];
                var k2 = ig[i + 1];
                MatrixAssistant.SparseAdd(res, jg.AsSpan().Slice(k1, k2 - k1), ggl.AsSpan().Slice(k1, k2 - k1), vec[i]);
            }
        }
        protected override void Add(double mult, SparseMatrix<double> matr)
        {
            if (matr is RowSparseNonSymSlaeOfDouble sl)
            {
                for (int i = 0; i < di.Length; i++) MatrixAssistant.ThreadSafeAddDouble(ref di[i], sl.di[i] * mult);
                for (int i = 0; i < ggl.Length; i++) MatrixAssistant.ThreadSafeAddDouble(ref ggl[i], sl.ggl[i] * mult);
                for (int i = 0; i < ggu.Length; i++) MatrixAssistant.ThreadSafeAddDouble(ref ggu[i], sl.ggu[i] * mult);
            }
            else if (matr is RowSparseSymmetricSlaeOfDouble sls)
            {
                for (int i = 0; i < di.Length; i++) MatrixAssistant.ThreadSafeAddDouble(ref di[i], sls.Di[i] * mult);
                for (int i = 0; i < ggl.Length; i++)
                {
                    var d = sls.Gg[i] * mult;
                    MatrixAssistant.ThreadSafeAddDouble(ref ggl[i], d);
                    MatrixAssistant.ThreadSafeAddDouble(ref ggu[i], d);
                }
            }
            else throw new ArgumentOutOfRangeException();
        }
        public override void Diagonalize()
        {
            if (ig[BlockSize] == 0) return;
            var sz = BlockSize;
            for (int i = 0; i < sz; i++)
            {
                var k1 = ig[i];
                var k2 = ig[i + 1];
                for (var j = k1; j < k2; j++)
                {
                    di[i] += ggl[j];
                    di[jg[j]] += ggu[j];
                }
            }
            ggl = null;
            ggu = null;
            for (int i = 0; i < ig.Length; i++) ig[i] = 0;
            jg = null;
        }

        public override void CorrectDiagonal()
        {
            for (int i = 0; i < di.Length; i++) if (di[i] == 0) di[i] = 1;
        }

        public override SparseMatrixType MatrixType => SparseMatrixType.ROWCOL_NONSYM_SPARSE;

        public override int SolverDepth { get; set; }
    }
    public class RowSparseNonSymSlaeOfComplex : RowSparseNonSymMatrix<Complex>
    {
        public override void AddElement(int tdi, int tdj, Complex val)
        {
            if (tdi == tdj)
                MatrixAssistant.ThreadSafeAddSpanElement(di, tdi, val);
            else
            {
                var addr = FindElem(tdi, tdj);
                if (tdi < tdj)
                {
                    MatrixAssistant.ThreadSafeAddSpanElement(ggu, addr, val);
                }
                else
                {
                    MatrixAssistant.ThreadSafeAddSpanElement(ggl, addr, val);
                }
            }
        }
        void AddElem(int i, int j, Complex dIJ, Complex dJI)
        {
            if (i == j)
            {
                MatrixAssistant.ThreadSafeAddSpanElement(di, i, dIJ);
            }
            else
            {
                var addr = FindElem(i, j);
                if (i < j)
                {
                    MatrixAssistant.ThreadSafeAddSpanElement(ggu, addr, dIJ);
                    MatrixAssistant.ThreadSafeAddSpanElement(ggl, addr, dJI);
                }
                else
                {
                    MatrixAssistant.ThreadSafeAddSpanElement(ggl, addr, dIJ);
                    MatrixAssistant.ThreadSafeAddSpanElement(ggu, addr, dJI);
                }
            }
        }
        void AddDiagElem(int i, Complex d)
        {
            MatrixAssistant.ThreadSafeAddSpanElement(di, i, d);
        }
        public override void MultVect(ReadOnlySpan<Complex> vec, Span<Complex> res)
        {
            var sz = BlockSize;
            for (int i = 0; i < sz; i++) res[i] = di[i] * vec[i];
            for (int i = 0; i < sz; i++)
            {
                var k1 = ig[i];
                var k2 = ig[i + 1];
                res[i] += MatrixAssistant.SparseMult(ggl.AsSpan().Slice(k1, k2 - k1), jg.AsSpan().Slice(k1, k2 - k1), vec);
            }
            for (int i = 0; i < sz; i++)
            {
                var k1 = ig[i];
                var k2 = ig[i + 1];
                MatrixAssistant.SparseAdd(res, jg.AsSpan().Slice(k1, k2 - k1), ggu.AsSpan().Slice(k1, k2 - k1), vec[i]);
            }
        }
        public override void TMultVect(ReadOnlySpan<Complex> vec, Span<Complex> res)
        {
            var sz = BlockSize;
            for (int i = 0; i < sz; i++) res[i] = di[i] * vec[i];
            for (int i = 0; i < sz; i++)
            {
                var k1 = ig[i];
                var k2 = ig[i + 1];
                res[i] += MatrixAssistant.SparseMult(ggu.AsSpan().Slice(k1, k2 - k1), jg.AsSpan().Slice(k1, k2 - k1), vec);
            }
            for (int i = 0; i < sz; i++)
            {
                var k1 = ig[i];
                var k2 = ig[i + 1];
                MatrixAssistant.SparseAdd(res, jg.AsSpan().Slice(k1, k2 - k1), ggl.AsSpan().Slice(k1, k2 - k1), vec[i]);
            }
        }
        protected override void Add(double mult, SparseMatrix<Complex> matr)
        {
            if (matr is RowSparseNonSymSlaeOfComplex sl)
            {
                for (int i = 0; i < di.Length; i++) MatrixAssistant.ThreadSafeAddSpanElement(di, i, sl.di[i] * mult);
                for (int i = 0; i < ggl.Length; i++) MatrixAssistant.ThreadSafeAddSpanElement(ggl, i, sl.ggl[i] * mult);
                for (int i = 0; i < ggu.Length; i++) MatrixAssistant.ThreadSafeAddSpanElement(ggu, i, sl.ggu[i] * mult);
            }
            else if (matr is RowSparseSymmetricMatrixOfComplex sls)
            {
                for (int i = 0; i < di.Length; i++) MatrixAssistant.ThreadSafeAddSpanElement(di, i, sls.Di[i] * mult);
                for (int i = 0; i < ggl.Length; i++)
                {
                    var d = sls.Gg[i] * mult;
                    MatrixAssistant.ThreadSafeAddSpanElement(ggl, i, d);
                    MatrixAssistant.ThreadSafeAddSpanElement(ggu, i, d);
                }
            }
            else throw new ArgumentOutOfRangeException();
        }
        public override void Diagonalize()
        {
            if (ig[BlockSize] == 0) return;
            var sz = BlockSize;
            for (int i = 0; i < sz; i++)
            {
                var k1 = ig[i];
                var k2 = ig[i + 1];
                for (var j = k1; j < k2; j++)
                {
                    di[i] += ggl[j];
                    di[jg[j]] += ggu[j];
                }
            }
            ggl = null;
            ggu = null;
            for (int i = 0; i < ig.Length; i++) ig[i] = 0;
            jg = null;
        }

        public override void CorrectDiagonal()
        {
            for (int i = 0; i < di.Length; i++) if (di[i] == 0) di[i] = 1;
        }

        public override SparseMatrixType MatrixType => SparseMatrixType.ROWCOL_NONSYM_SPARSE;

        public override int SolverDepth { get; set; }
    }
    public abstract class RowOnlySparseNonSymSlae<T> : SparseMatrix<T> where T : struct
    {
        protected int[] ig;
        protected int[] jg;
        protected int[] diaddr;
        protected T[] gg;
        public override void Nullify()
        {
            gg.AsSpan().Fill(default);
            if(Solver != null)SetSolverAndFactorizationType(Solver.DisplayCode, Solver.Factorizer == null ? FactorizerType.NONE : Solver.Factorizer.DisplayCode);
        }
        public override void SetPortrait(ISparsePortrait portrait)
        {
            var NUnknowns = portrait.RowSize;
            diaddr = new int[NUnknowns];
            ig = new int[NUnknowns + 1];
            ig[0] = 0;
            for (int i = 0; i < NUnknowns; i++)
            {
                ig[i + 1] = ig[i] + portrait.Row(i).Count();
            }
            jg = new int[ig[NUnknowns]];
            for (int i = 0; i < NUnknowns; i++)
            {
                var jgaddr = ig[i];
                foreach (var j in portrait.Row(i).OrderBy(m => m))
                {
                    if (j == i) diaddr[i] = jgaddr;
                    jg[jgaddr++] = j;
                }
            }
            gg = new T[jg.Length];
            AllocateArrays();
        }
        public override T GetDi(int i) => gg[diaddr[i]];
        protected override void SetDi(int i, T val)
        {
            gg[diaddr[i]] = val;
        }
        protected int FindElem(int i, int j)// возвращает индекс в массиве jg
        {
            return MatrixAssistant.QuickFindInOrderedIgJg(ig, jg, i, j);
        }
        public override int Size => BlockSize;
        public int BlockSize => ig.Length - 1; // Для неблочной матрицы совпадает с обычным размером
        public override bool IsSymmetry => false;
        public override void ClearRow(int row)
        {

            for (var i = ig[row]; i < ig[row + 1]; i++) gg[i] = default;
        }
        public override bool IsMatrixValid => BlockSize > 0 && diaddr.Length == BlockSize && jg.Length == ig[BlockSize] && gg.Length == jg.Length;
        protected override void AllocateArrays()
        {
            var sz = BlockSize;
            gg = new T[jg.Length];
        }
        public override void PrintFormatted(TextWriter writer, int row = -1)
        {
            throw new NotImplementedException();
        }
        public override void Read(string path)
        {
            throw new NotImplementedException();
        }
        public override void Write(string path)
        {
            throw new NotImplementedException();
        }
        public override int SolverDepth { get; set; }
        public override IEnumerable<(int, int, T)> AllNonZeroElements
        {
            get
            {
                var sz = BlockSize;
                for (int i = 0; i < sz; i++)
                {
                    var k1 = ig[i];
                    var k2 = ig[i + 1];
                    for (var j = k1; j < k2; j++) yield return (i, jg[j], gg[j]);
                }
            }
        }
        public override void Assign(ISparseMatrix<T> matr)
        {
            if (matr is RowOnlySparseNonSymSlae<T> sl)
            {
                ig = (int[])sl.ig.Clone();
                jg = (int[])sl.jg.Clone();

                diaddr = (int[])sl.diaddr.Clone();

                gg = (T[])sl.gg.Clone();


            }
            throw new ArgumentOutOfRangeException("BadMatrixType");
        }
    }
    public class RowOnlySparseNonSymSlaeOfDouble : RowOnlySparseNonSymSlae<double>
    {
        void SetElement(int tdi, int tdj, double val)
        {
            if (tdi == tdj)
                MatrixAssistant.ThreadSafeAddDouble(ref gg[diaddr[tdi]], val);
            else
            {
                var addr = FindElem(tdi, tdj);
                MatrixAssistant.ThreadSafeAddDouble(ref gg[addr], val);
            }
        }
        void AddElem(int i, int j, double dIJ, double dJI)
        {
            AddElement(i, j, dIJ);
            AddElement(j, i, dJI);
        }
        void AddDiagElem(int i, double d)
        {
            MatrixAssistant.ThreadSafeAddDouble(ref gg[diaddr[i]], d);
        }
        public override void AddElement(int tdi, int tdj, double val)
        {
            if (tdi == tdj)
                MatrixAssistant.ThreadSafeAddDouble(ref gg[diaddr[tdi]], val);
            else
            {
                var addr = FindElem(tdi, tdj);
                MatrixAssistant.ThreadSafeAddDouble(ref gg[addr], val);
            }
        }
        public override void MultVect(ReadOnlySpan<double> vec, Span<double> res)
        {
            var sz = BlockSize;
            for (int i = 0; i < sz; i++) res[i] = 0;
            for (int i = 0; i < sz; i++)
            {
                var k1 = ig[i] ;
                var k2 = ig[i + 1] ;
                //                            res[i]+=gg.SparseMult(k1,k2-k1,jg,vec);
                for (var j = k1; j < k2; j++) res[i] += gg[j] * vec[jg[j] ];
            }
        }
        public override void TMultVect(ReadOnlySpan<double> vec, Span<double> res)
        {
            var sz = BlockSize;
            for (int i = 0; i < sz; i++) res[i] = 0;
            for (int i = 0; i < sz; i++)
            {
                var k1 = ig[i] ;
                var k2 = ig[i + 1] ;
                //                            res[i]+=gg.SparseMult(k1,k2-k1,jg,vec);
                for (var j = k1; j < k2; j++) res[jg[j]] += gg[j] * vec[i];
            }
        }
        protected override void Add(double mult, SparseMatrix<double> matr)
        {
            if (matr is RowOnlySparseNonSymSlaeOfDouble sl)
            {
                for (int i = 0; i < sl.gg.Length; i++) MatrixAssistant.ThreadSafeAddDouble(ref gg[i], mult * sl.gg[i]);
                return;
            }
            throw new ArgumentOutOfRangeException("BadMatrixType");
        }
        public override void Diagonalize()
        {
            throw new NotImplementedException();
        }
        public override void CorrectDiagonal()
        {
            for (int i = 0; i < BlockSize; i++) if (GetDi(i) == 0) SetDi(i, 1);
        }

        public override SparseMatrixType MatrixType => SparseMatrixType.ROW_SPARSE_FOR_PARDISO;
    }
    public class RowOnlySparseNonSymSlaeOfComplex: RowOnlySparseNonSymSlae<Complex>
    {
        void SetElement(int tdi, int tdj, Complex val)
        {
            if (tdi == tdj)
                MatrixAssistant.ThreadSafeAddSpanElement(gg,diaddr[tdi], val);
            else
            {
                var addr = FindElem(tdi, tdj);
                MatrixAssistant.ThreadSafeAddSpanElement(gg,addr, val);
            }
        }
        void AddElem(int i, int j, Complex dIJ, Complex dJI)
        {
            AddElement(i, j, dIJ);
            AddElement(j, i, dJI);
        }
        void AddDiagElem(int i, Complex d)
        {
            MatrixAssistant.ThreadSafeAddSpanElement(gg,diaddr[i], d);
        }
        public override void AddElement(int tdi, int tdj, Complex val)
        {
            if (tdi == tdj)
                MatrixAssistant.ThreadSafeAddSpanElement(gg,diaddr[tdi], val);
            else
            {
                var addr = FindElem(tdi, tdj);
                MatrixAssistant.ThreadSafeAddSpanElement(gg,addr, val);
            }
        }
        public override void MultVect(ReadOnlySpan<Complex> vect, Span<Complex> res)
        {
            var sz = BlockSize;
            res.Fill(0);
            for (int i = 0; i < sz; i++)
            {
                var k1 = ig[i] ;
                var k2 = ig[i + 1] ;
                //                            res[i]+=gg.SparseMult(k1,k2-k1,jg,vec);
                for (var j = k1; j < k2; j++) res[i] += gg[j] * vect[jg[j] ];
            }
        }
        public override void TMultVect(ReadOnlySpan<Complex> vec, Span<Complex> res)
        { 
            var sz = BlockSize;
            res.Fill(0);
            for (int i = 0; i < sz; i++)
            {
                var k1 = ig[i] ;
                var k2 = ig[i + 1] ;
                //                            res[i]+=gg.SparseMult(k1,k2-k1,jg,vec);
                for (var j = k1; j < k2; j++) res[jg[j] ] += gg[j] * vec[i];
            }
        }
        protected override void Add(double mult, SparseMatrix<Complex> matr)
        {
            if (matr is RowOnlySparseNonSymSlaeOfComplex sl)
            {
                for (int i = 0; i < sl.gg.Length; i++) MatrixAssistant.ThreadSafeAddSpanElement(gg,i, mult * sl.gg[i]);
                return;
            }
            throw new ArgumentOutOfRangeException("BadMatrixType");
        }
        public override void Diagonalize()
        {
            throw new NotImplementedException();
        }
        public override void CorrectDiagonal()
        {
            for (int i = 0; i < BlockSize; i++) if (GetDi(i) == 0) SetDi(i, 1);
        }

        public override SparseMatrixType MatrixType => SparseMatrixType.ROW_SPARSE_FOR_PARDISO;
    }
}