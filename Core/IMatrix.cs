using System;
using System.Collections.Generic;
using System.IO;
using System.Numerics;
using System.Runtime.InteropServices;
using System.Threading;
using System.Threading.Tasks;

namespace BlazorServer.Matrix
{
    public interface IMatrixBase
    {
        /// <summary>
        /// Число строк
        /// </summary>
        int RSize { get; }
        /// <summary>
        /// число столбцов
        /// </summary>
        int CSize { get; }
    }
    public interface IRectangularMatrix<T> : IMatrixBase where T : struct
    {
        void MultVect(ReadOnlySpan<T> to, Span<T> res); //! Умножение на вектор
        void TMultVect(ReadOnlySpan<T> to, Span<T> res);//! Умножение на вектор транспонированной матрицы
        /// <summary>
        /// Возвращает все ненулевые элементы (например, для занесения в глобальную). 
        /// Для симметричной марицы должна пройти по двум половинкам !!!
        /// выбрасывает исключение NotImplementedException если процедура невозможна
        /// </summary>
        IEnumerable<(int, int, T)> AllNonZeroElements { get; }

        void SymmetryScale(ReadOnlySpan<double> Multiplier);//! Умножение строк и столбцов на веса
    }

    public interface ISquareMatrix<T> : IRectangularMatrix<T> where T : struct
    {
        int Size { get; } //! Размерность для квадратной матрицы
    }
    public interface IInversableMatrix<T> : ISquareMatrix<T> where T : struct
    {
        void InverseMult(ReadOnlySpan<T> to, Span<T> res);
    }
    public enum FactorizerType { LUSQ, LU, DI, NONE };

    public interface IMatrixFactorizerBase : IDisposable
    {
        FactorizerType DisplayCode { get; }
        void Factorize(/* IProgressBar progress = null,  */CancellationToken? cancellationToken = null);
    }
    public enum SolverType { CG, LOS, GMRES, BiCG, BiCGStab, PARDISO };
    public interface ISlaeSolverBase : IDisposable
    {
        IMatrixFactorizerBase Factorizer { get; }
        //bool IsEnabledFor(ISparseMatrixBase slae);
        //bool CanUse(IMatrixFactorizerBase factorizer);
        SolverType DisplayCode { get; }
    };
    public enum SparseMatrixType { ROWCOL_SYM_SPARSE, ROWCOL_NONSYM_SPARSE, ROW_SPARSE_FOR_PARDISO, ROW_SYM_SPARSE_FOR_PARDISO, ROWCOL_SYM_SPARSE_Harmonic, ROWCOL_NONSYM_SPARSE_Harmonic, ROW_SPARSE_HARMONIC_FOR_PARDISO };

    public interface ISparsePortrait //всегда полный (оба треугольника и диагональ) !!!
    {
        void Add(IEnumerable<int> rows, IEnumerable<int> columns);
        void Add(ReadOnlySpan<int> rows, ReadOnlySpan<int> columns);
        IEnumerable<int> Row(int row);
        IEnumerable<int> Column(int column);
        int RowSize { get; }
        int ColumnSize { get; }
    }
    public interface ISparseMatrixBase : IMatrixBase
    {
        ISlaeSolverBase Solver { get; }
        /// <summary>
        /// форматный вывод матрицы
        /// </summary>
        /// <param name="writer">куда выводить</param>
        /// <param name="row">строка, -1 означает все</param>
        void PrintFormatted(System.IO.TextWriter writer, int row = -1);
        SparseMatrixType MatrixType { get; }
        void SetPortrait(ISparsePortrait portrait);
        bool IsSymmetry { get; }
        bool IsMatrixValid { get; }
        /// <summary>
        /// только для отладки. Пишет файлы с фиксированными именами по указанному пути
        /// </summary>
        void Write(string path);
        void Read(string path);
        /// <summary>
        /// корректировка главной дигонали СЛАУ (ставит единички в нулевых строках)
        /// </summary>
        void CorrectDiagonal();
        void SetSolverAndFactorizationType(SolverType solver, FactorizerType factorizer);
    }
    public interface IMatrixFactorizer<T> : IMatrixFactorizerBase where T : struct
    {
        /// <summary>
        /// //! умножение на L^-1
        /// </summary>
        /// <param name="vector">результат в него запишется</param>
        void LMult(Span<T> vector);
        /// <summary>
        /// умножение на U^-1
        /// </summary>
        /// <param name="vector">результат в него запишется</param>
		void UMult(Span<T> vector);
        /// <summary>
        /// //! умножение на U
        /// </summary>
		void UxMult(ReadOnlySpan<T> to, Span<T> res);
    };
    public interface ISparseMatrix<T> : ISparseMatrixBase where T : struct
    {
        new ISlaeSolver<T> Solver { get; }
        T GetDi(int i);//!< вернуть i элемент диагонали
        void DiagonaleSolve(Span<T> u);
        void Diagonalize(); //! Диагонализация матрицы массы для явной схемы
        void Nullify();
    }
    public interface ISlaeSolverParameters
    {
        double Epsru0 { get; }//выход по относительной невязке
        double Epsrvh { get; }//выход по относительной невязке от входной
        int Maxiter { get; }
        int GMResDepth { get; }
    }
    public interface ISlaeSolver<T> : ISlaeSolverBase where T : struct
    {
        // ILinearOperator<T> Matrix { get; set; }
        new IMatrixFactorizer<T> Factorizer { get; set; }
        void UseNeumannCorrection(ReadOnlySpan<T> vector);
        /// <summary>
        /// начальное приближение и решение в нумерации СЛАУ
        /// </summary>
        /// <summary>
        ////Решение СЛАУ
        /// </summary>
        /// <param name="res">возвращаемый результат, содержит начальное приближение
        /// <param name="progress">информационная строка, текущая длина, максимальная длина</param>
        /// <param name="ct"></param>
        /// <returns>результат, невязка, число итераций</returns>
        (double residual, int iter) Solve(ReadOnlySpan<T> RightPart, Span<T> res, ISlaeSolverParameters parameters = null, /* IProgressBar progress = null, */
                               CancellationToken? ct = null);
    }

}
