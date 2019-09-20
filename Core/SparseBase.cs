using System;
using System.Collections.Generic;
using System.IO;
using System.Numerics;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading;
using System.Threading.Tasks;

namespace BlazorServer.Matrix
{
    public interface ILinearOperator<T> : ISquareMatrix<T>, ISparseMatrix<T> where T : struct
    {

    }
    public interface ISquareSparseMatrix<T>: ILinearOperator<T> where T : struct
    {
        void AddElement(int i, int j, T elem);
        void Assign(ISparseMatrix<T> matr);
    }
    public abstract class SparseMatrix<T> : ISquareSparseMatrix<T> where T : struct
    {
        protected ISlaeSolver<T> solver;

        public abstract void PrintFormatted(TextWriter writer, int row = -1);
        //        public abstract void AllocateData();
        //        public abstract void ClearMatrix();
        public abstract void Write(string path);
        public abstract void Read(string path);
        public abstract void CorrectDiagonal();

        public ISlaeSolver<T> Solver => solver;
        public abstract SparseMatrixType MatrixType { get; }
        public abstract bool IsSymmetry { get; }
        public abstract bool IsMatrixValid { get; }
        public abstract int SolverDepth { get; set; }
        public int RSize => Size;
        public int CSize => Size;
        public abstract int Size { get; }

        ISlaeSolverBase ISparseMatrixBase.Solver => solver;

        public abstract IEnumerable<(int, int, T)> AllNonZeroElements { get; }

        public abstract void Assign(ISparseMatrix<T> matr);
        protected abstract void Add(double mult, SparseMatrix<T> matr);
        //		protected void SwapRightPart(vector<T>& newrp) { pr.swap(newrp); }
        public void SetSolverAndFactorizationType(SolverType stype, FactorizerType ftype)
        {
            solver?.Dispose();
            solver = MatrixAssistant.CreateSolver<T>(stype);
            if (!stype.IsEnabledFor(MatrixType)) throw new Exception("Bad solver " + solver.DisplayName());
            solver.Matrix = this;
            if (ftype != FactorizerType.NONE)
            {
                var spfact = (IMatrixFactorizer<T>)MatrixAssistant.CreateFactorizer(ftype, this);
                if (spfact == null) throw new Exception("BadSlaeFactorizer " + MatrixAssistant.DisplayName(ftype));
                if (!stype.CanUse(ftype)) throw new Exception(solver.DisplayName() + " can't use factorizer " + spfact.DisplayName());
                solver.Factorizer = spfact;
            }
        }
        /// <summary>
        /// Считает невяку. Возвращает квадрат нормы невязки и квадрат нормы правой части
        /// </summary>        
        //public void CalcResidual(TMatrixStruct TMatr, ReadOnlySpan<T> vect, Span<T> residual,out double f1, out double f2)
        //{
        //    var u = Solver.ResVector.Span;
        //    var v = vect.Span;
        //    for (int i = 0; i < u.Length; i++) u[i] = v[TMatr.Reorder.UsedUnk(i)];
        //    CalcInternalResidual(Solver.ResVector, out f1, out f2);
        //}
        public void CalcInternalResidual(ReadOnlySpan<T> Pr, ReadOnlySpan<T> vect, Span<T> residual, out double f1, out double f2)  //! Возвращает квадрат нормы невязки и квадрат нормы правой части
        {
            MultVect(vect, residual);
            // посчитать невязку
            f2 = Pr.NormSqr();
            long sz = Size;
            residual.Add(Pr, -1);
            residual.Assign(residual, -1.0);
            f1 = residual.NormSqr();
        }

        public abstract T GetDi(int i);//!< вернуть i элемент диагонали
        protected abstract void SetDi(int i, T val);
        public void DiagonaleSolve(Span<T> u)
        {
            if (typeof(T) == typeof(double))
            {
                var dbl = MemoryMarshal.Cast<T, double>(u);
                for (int i = 0; i < Size; i++)
                {
                    dbl[i] /= (double)(ValueType)GetDi(i);
                }
            }
            else if (typeof(T) == typeof(Complex))
            {
                var cplx = MemoryMarshal.Cast<T, Complex>(u);
                for (int i = 0; i < Size; i++)
                {
                    cplx[i] /= (Complex)(ValueType)GetDi(i);
                }
            }
            else throw new NotImplementedException();
        }
        public abstract void Diagonalize(); //! Диагонализация матрицы массы для явной схемы
        public abstract void ClearRow(int row);
        public abstract void AddElement(int i, int j, T elem);
        public virtual void AddLocalMatrix(ReadOnlySpan<int> indrows, ReadOnlySpan<int> indcolumns, IRectangularMatrix<T> matr)//!< добавить локальную матрицу в СЛАУ
        {
            foreach(var elem in matr.AllNonZeroElements)
            {
                AddElement(indrows[elem.Item1], indcolumns[elem.Item2], elem.Item3);
            }
        }
        //public abstract void AddLocalMatrix(TMatrixStruct TMatr, ReadOnlySpan<int> ind, IRectangularMatrix<T> matr);//!< добавить локальную матрицу в СЛАУ

        public abstract void SetPortrait(ISparsePortrait portrait);
        public abstract void MultVect(ReadOnlySpan<T> to, Span<T> res);
        public abstract void TMultVect(ReadOnlySpan<T> to, Span<T> res);
        protected abstract void AllocateArrays();

        public void SymmetryScale(ReadOnlySpan<double> Multiplier)
        {
            throw new NotImplementedException();
        }

        public abstract void Nullify();
    }
}
