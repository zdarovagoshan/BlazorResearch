using System;
using System.Collections.Generic;
using System.IO;
using System.Numerics;
using System.Text;
using System.Threading;
using System.Threading.Tasks;

namespace BlazorServer.Matrix
{
    //double epsru0 = 1E-15
    //         double epsrvh = 0
    //        , int mkiter = 10000,
    public class SlaeSolverParameters : ISlaeSolverParameters
    {
        public double Epsru0 { get; set; }

        public double Epsrvh { get; set; }

        public int Maxiter { get; set; }

        public int GMResDepth { get; set; }
        class DefSlaeSolverParameters : ISlaeSolverParameters
        {
            public double Epsru0 => ContextManager.Context.Config.Parameter("Default SLAE Solver epsru0", 1e-15);

            public double Epsrvh => ContextManager.Context.Config.Parameter("Default SLAE Solver epsrvh", 0);

            public int Maxiter => ContextManager.Context.Config.Parameter("Default SLAE Solver maxiter", 10000);

            public int GMResDepth => ContextManager.Context.Config.Parameter("Default SLAE Solver Depth", 5);
        }
        static DefSlaeSolverParameters DefVal = new DefSlaeSolverParameters();
        public static ISlaeSolverParameters Default() => new SlaeSolverParameters();
        public SlaeSolverParameters()
        {
            Epsru0 = DefVal.Epsru0;
            Epsrvh = DefVal.Epsrvh;
            Maxiter = DefVal.Maxiter;
            GMResDepth = DefVal.GMResDepth;
        }
    }

    public abstract class SlaeSolver<T> : ISlaeSolver<T> where T : struct
    {
        public ILinearOperator<T> Matrix { get; set; }
        public IMatrixFactorizer<T> Factorizer { get; set; }

        IMatrixFactorizerBase ISlaeSolverBase.Factorizer => Factorizer;

        public abstract SolverType DisplayCode { get; }

        protected Dictionary<string, T[]> workvector = new Dictionary<string, T[]>(); //массив для хранения рабочих векторов методов решения.


        protected T[] GetWorkVector(string name)
        {
            if (workvector.TryGetValue(name, out var val)) return val;
            var newvec = new T[Matrix.Size];
            workvector[name] = newvec;
            return newvec;
        }
        public static long ProgressLength(double residual) => (long)(-Math.Log10(residual) * 1000);
        /// <summary>
        /// res may contain initial values
        /// </summary>
        public (double residual, int iter) Solve(ReadOnlySpan<T> RightPart, Span<T> res, ISlaeSolverParameters parameters = null, IProgressBar progress = null, CancellationToken? ct = null)
        {
            if (parameters == null) parameters = SlaeSolverParameters.Default();
            return SolveMethod(RightPart, res, parameters, progress, ct);
        }
        protected abstract (double residual, int iterations) SolveMethod(ReadOnlySpan<T> Pr, Span<T> res, ISlaeSolverParameters parameters, IProgressBar progress, CancellationToken? ct);

        #region IDisposable Support
        private bool disposedValue = false; // To detect redundant calls

        protected virtual void Dispose(bool disposing)
        {
            if (!disposedValue)
            {
                if (disposing)
                {
                    Factorizer?.Dispose();
                    // TODO: dispose managed state (managed objects).
                    Factorizer = null;
                    Matrix = null;
                    workvector.Clear();
                }
                // TODO: free unmanaged resources (unmanaged objects) and override a finalizer below.
                // TODO: set large fields to null.

                disposedValue = true;
            }
        }

        // TODO: override a finalizer only if Dispose(bool disposing) above has code to free unmanaged resources.
        // ~SlaeSolver() {
        //   // Do not change this code. Put cleanup code in Dispose(bool disposing) above.
        //   Dispose(false);
        // }

        // This code added to correctly implement the disposable pattern.
        public void Dispose()
        {
            // Do not change this code. Put cleanup code in Dispose(bool disposing) above.
            Dispose(true);
            // TODO: uncomment the following line if the finalizer is overridden above.
            // GC.SuppressFinalize(this);
        }


        #endregion
        public abstract void UseNeumannCorrection(ReadOnlySpan<T> vector);
    };

    public class CGSolver<T> : SlaeSolver<T> where T : struct
    {
        T[] newmanncorrector = null;
        T[] LUnewmanncorrector = null;
        T newmanncorrectormultiplier=default;
        public override void UseNeumannCorrection(ReadOnlySpan<T> vector)
        {
            newmanncorrector = vector.ToArray();
        }
        protected override (double residual, int iterations) SolveMethod(ReadOnlySpan<T> Pr, Span<T> u, ISlaeSolverParameters parameters, IProgressBar progress, CancellationToken? ct)
        {
            /*Meтoд conpяжeнныx гpaдueнтoв c нenoлныm paзлoжeнuem Xoлecckoгo
для CЛAY c cummeтpuчнoй maтpuцeй в paзpeжeннom cтpoчнom фopmaтe
(maтpuцa зaдaнa cвoum нuжнum тpeyг-kom)*/
            var sz = Matrix.Size;
            if (u == null || u.Length != sz) u = new T[sz];
            var arithmetic = TemplateArithmetic.Create<T>();
            Factorizer?.Factorize(progress, ct);
            if(progress!=null) progress.Size = ProgressLength(parameters.Epsru0);
            progress?.Report(("Matrix solving by CG method", 0));
            var r = GetWorkVector("r");
            var z = GetWorkVector("z");
            var az = GetWorkVector("az");
            var rab = GetWorkVector("rab");
            const double rmin = 1e-100;
            Matrix.CalcResidual(Pr,u, r);

            var rvhod = (r.AsSpan().Norm() + rmin) / sz;
            var ru00 = (Pr.Norm() + rmin) / sz;
            //		long double delrnl=rvhod/ru00;
            var delr = rvhod;
            if (newmanncorrector != null)
            {
                LUnewmanncorrector = newmanncorrector.AsSpan().ToArray();
                Factorizer?.LMult(LUnewmanncorrector);
                Factorizer?.UMult(LUnewmanncorrector);
                newmanncorrectormultiplier = LUnewmanncorrector.AsSpan().Dot(newmanncorrector);
                var c = u.Dot(newmanncorrector);
                c = arithmetic.UnaryMinus(c);
                r.AsSpan().Add(newmanncorrector, c);
            }
            // вычucляem z=(L')**(-1)*L**(-1)*r
            z.AsSpan().Assign(r);
            Factorizer?.LMult(z);
            Factorizer?.UMult(z);
            if(newmanncorrector != null)
            {
                var c = r.AsSpan().Dot(LUnewmanncorrector);
                c = arithmetic.UnaryMinus(arithmetic.Devide(c, arithmetic.Sum(newmanncorrectormultiplier, 1)));
                z.AsSpan().Add(LUnewmanncorrector, c);
            }
            T rrprev = z.AsSpan().Mult(r);

            // Haчuнaem uтepaцuu
            var delrvh = delr / rvhod;
            var delru0 = delr / ru00;
            int iter;
            ContextManager.Context.Logger.Log(LogLevel.Information, $"SLAE solving CG:eps={parameters.Epsru0},del={delru0} start");

            for (iter = 0; iter <parameters.Maxiter && delru0 > parameters.Epsru0  && delrvh >parameters.Epsrvh; iter++)
            {
                ct?.ThrowIfCancellationRequested();
                progress?.Report(($"SLAE solving CG:eps={parameters.Epsru0},del={delru0},iter={iter}", ProgressLength(delru0)));
                ContextManager.Context.Logger.Log(LogLevel.Debug, $"SLAE solving CG:eps={parameters.Epsru0},del={delru0} iter={iter}");
                T alk, betk, azz, rrnew;
                Matrix.MultVect(z, az);
                if(newmanncorrector != null)az.AsSpan().Add(newmanncorrector, z.AsSpan().Dot(newmanncorrector));

                azz = az.AsSpan().Dot(z);
                //			if (azz<=-rmin) throw new NegativeDi("The matrix is not positivly defined.");

                if (arithmetic.Abs(azz) <= rmin) break;
                alk = arithmetic.Devide(rrprev, azz);
                delr = 0.0;
                u.Add(z, alk);
                r.AsSpan().Add(az, arithmetic.UnaryMinus(alk));

                // вычucляem rab=(L')**(-1)*L**(-1)*r
                rab.AsSpan().Assign(r);
                Factorizer?.LMult(rab);
                Factorizer?.UMult(rab);
                if (newmanncorrector != null)
                {
                    var c = r.AsSpan().Dot(LUnewmanncorrector);
                    c = arithmetic.UnaryMinus(arithmetic.Devide(c, arithmetic.Sum(newmanncorrectormultiplier, 1)));
                    rab.AsSpan().Add(LUnewmanncorrector, c);
                }

                rrnew = rab.AsSpan().Dot(r);
                delr = Math.Sqrt(arithmetic.Abs(rrnew)) / sz;

                betk = arithmetic.Devide(rrnew, rrprev);
                rab.AsSpan().Add(z, betk);

                z.AsSpan().Assign(rab);
                rrprev = rrnew;

                delrvh = delr / rvhod;
                delru0 = delr / ru00;
            }
#if DEBUG
            double rrist = 0.0;
            Matrix.CalcResidual(Pr,u, rab);
            rrist = rab.AsSpan().Norm() / sz;
            if (!(rrist / ru00 < 1e-8/*100000*epsru0*/ || rrist < 1e-15)) throw new Exception("Assertion failed");
            delru0 = rrist / ru00;
#endif
            ContextManager.Context.Logger.Log(LogLevel.Information, $"SLAE solving CG:eps={parameters.Epsru0},del={delru0} end");
            return (delru0, iter);
        }
        public override SolverType DisplayCode => SolverType.CG;
    }

    public class LOSSolver<T> : SlaeSolver<T> where T : struct
    {
        protected override (double residual, int iterations) SolveMethod(ReadOnlySpan<T> Pr, Span<T> u, ISlaeSolverParameters parameters, IProgressBar progress, CancellationToken? ct)
        {
            /* Локально-оптимальная схема c нenoлныm LU (любым) paзлoжeнuem
            для CЛAY c неcummeтpuчнoй maтpuцeй в paзpeжeннom cтpoчнom фopmaтe*/
            int size = Matrix.Size;
            var arithmetic = TemplateArithmetic.Create<T>();
            if (u == null || u.Length != size) u = new T[size];
            if (Factorizer != null) Factorizer.Factorize(progress,ct);
            if(progress != null)progress.Size=ProgressLength(parameters.Epsru0);
            progress?.Report(("Matrix solving by LOS method", 0));
            var r = GetWorkVector("r");
            var z = GetWorkVector("z");
            var az = GetWorkVector("az");
            var rab = GetWorkVector("rab");
            var rab2 = GetWorkVector("rab2");
            const double rmin = 1e-100;

            // вычисляем первые к.у.

            // Bычucляem нeвязky r, conpяжeннoe нanpaвлeнue z=r, rvhod, ru00 u delrnl
            Matrix.CalcResidual(Pr,u, r);
            var rvhod = (r.AsSpan().Norm() + rmin) / size;
            var ru00 = (Pr.Norm() + rmin) / size;
            double delr = rvhod;
            // вычucляem z=(L')**(-1)*L**(-1)*r
            Factorizer?.LMult(r);
            z.AsSpan().Assign(Pr);
            Factorizer?.LMult(z);
            rvhod = r.AsSpan().Norm();
            ru00 = z.AsSpan().Norm();

            rvhod = (rvhod + rmin) / size;
            ru00 = (ru00 + rmin) / size;
            //	 long double delrnl=rvhod/ru00;
            delr = rvhod;
            z.AsSpan().Assign(r);
            Factorizer?.UMult(z);
            Matrix.MultVect(z, az);
            Factorizer?.LMult(az);

            // Haчuнaem uтepaцuu
            double delrvh = delr / rvhod;
            double delru0 = delr / ru00;
            //	bool IsAzz=false;
            delr = r.AsSpan().NormSqr();
            int iter;
            ContextManager.Context.Logger.Log(LogLevel.Information, $"SLAE solving LOS:eps={parameters.Epsru0},del={delru0} start");
            for (iter = 0; iter < parameters.Maxiter && delru0 > parameters.Epsru0 && delrvh > parameters.Epsrvh; iter++)

            {
                // Oдuн шaг peш-я в блoke; npu вxoдe нa 1-oй uт-uu в az д.б.
                // знaч-e A*z, npu вxoдe нa nocлeдyющux uт-яx в rab бyдeт знaч-e A*z
                ct?.ThrowIfCancellationRequested();
                progress?.Report(($"SLAE solving:eps={parameters.Epsru0},del={delru0},iter={iter}", ProgressLength(delru0)));
                ContextManager.Context.Logger.Log(LogLevel.Debug, $"SLAE solving LOS:eps={parameters.Epsru0},del={delru0},iter={iter}");
                T alk, betk;
                var azaz = az.AsSpan().NormSqr();
                if (azaz <= rmin)
                {
                    //			IsAzz=true;
                    //			 throw new TelmaException("azaz=%le",azaz);
                    break;
                }
                T azr = az.AsSpan().Dot(r);
                alk = arithmetic.Devide(azr, azaz);
                u.Add(z, alk);
                r.AsSpan().Add(az, arithmetic.UnaryMinus(alk));
                delr -= arithmetic.Abs(arithmetic.Multiply(alk, alk)) * azaz;
                if (delr <= rmin) delr = r.AsSpan().NormSqr();
                rab.AsSpan().Assign(r);
                Factorizer?.UMult(rab);
                Matrix.MultVect(rab, rab2);
                Factorizer?.LMult(rab2);
                T azar = az.AsSpan().Dot(rab2);
                betk = arithmetic.Devide(arithmetic.UnaryMinus(azar), azaz);
                rab.AsSpan().Add(z, betk);
                z.AsSpan().Assign(rab);
                rab2.AsSpan().Add(az, betk);
                az.AsSpan().Assign(rab2);
                double dd = (Math.Sqrt(delr) + rmin) / size;
                delrvh = dd / rvhod;
                delru0 = dd / ru00;
            }

#if DEBUG

            Matrix.CalcResidual(Pr,u, rab);
            var rrist = rab.AsSpan().Norm() / size;
            delru0 = rrist / ru00;
#endif
            ContextManager.Context.Logger.Log(LogLevel.Information, $"SLAE solving LOS:eps={parameters.Epsru0},del={delru0} end");
            return (delru0, iter);
        }
        public override SolverType DisplayCode => SolverType.LOS;

        public override void UseNeumannCorrection(ReadOnlySpan<T> vector)
        {
            throw new NotSupportedException();
        }
    }
    public class BiCGStabSolver<T> : SlaeSolver<T> where T : struct
    {
        protected override (double residual, int iterations) SolveMethod(ReadOnlySpan<T> Pr, Span<T> u, ISlaeSolverParameters parameters, IProgressBar progress, CancellationToken? ct)
        {
            /* Метод бисопряженных градиентов стабилизированный c нenoлныm LU (любым) paзлoжeнueм
            для CЛAY c неcummeтpuчнoй maтpuцeй в paзpeжeннoмт cтpoчнoм фopмaтe*/
            int size = Matrix.Size;
            var arithmetic = TemplateArithmetic.Create<T>();
            if (u == null || u.Length != size) u = new T[size];
            if (Factorizer != null) Factorizer.Factorize(progress, ct);
            if (progress != null) progress.Size = ProgressLength(parameters.Epsru0);
            progress?.Report(("Matrix solving by BiCGStab method", 0));
            var r = GetWorkVector("r");
            var r0 = GetWorkVector("r0");
            var z = GetWorkVector("z");
            var rpred = GetWorkVector("rpred");
            var p = GetWorkVector("p");
            var az = GetWorkVector("az");
            var rab = GetWorkVector("rab");
            var rab2 = GetWorkVector("rab2");
            const double rmin = 1e-100;

            // вычисляем первые к.у.

            // Bычucляem нeвязky r, conpяжeннoe нanpaвлeнue z=r, rvhod, ru00 u delrnl
            Matrix.CalcResidual(Pr, u, r);
            var rvhod = (r.AsSpan().Norm() + rmin);
            var ru00 = (Pr.Norm() + rmin);
            double delr = rvhod;
            Factorizer?.LMult(r);
            z.AsSpan().Assign(Pr);
            rvhod = r.AsSpan().Norm();
            ru00 = z.AsSpan().Norm();

            rvhod = (rvhod + rmin);
            ru00 = (ru00 + rmin);
            //	 long double delrnl=rvhod/ru00;
            delr = rvhod;
            z.AsSpan().Assign(r);
            r0.AsSpan().Assign(r);
            p.AsSpan().Assign(r);
          
            // Haчuнaem uтepaцu1u
            double delrvh = delr / rvhod;
            double delru0 = delr / ru00;
            //	bool IsAzz=false;
            delr = r.AsSpan().NormSqr();
            int iter;
            ContextManager.Context.Logger.Log(LogLevel.Information, $"SLAE solving BiCGStab:eps={parameters.Epsru0},del={delru0} start");
            for (iter = 0; iter < parameters.Maxiter && delru0 > parameters.Epsru0 && delrvh > parameters.Epsrvh; iter++)

            {
                // Oдuн шaг peш-я в блoke; npu вxoдe нa 1-oй uт-uu в az д.б.
                // знaч-e A*z, npu вxoдe нa nocлeдyющux uт-яx в rab бyдeт знaч-e A*z
                ct?.ThrowIfCancellationRequested();
                progress?.Report(($"SLAE solving:eps={parameters.Epsru0},del={delru0},iter={iter}", ProgressLength(delru0)));
                ContextManager.Context.Logger.Log(LogLevel.Debug, $"SLAE solving BiCGStab:eps={parameters.Epsru0},del={delru0},iter={iter}");
                T alk, betk,gamk;
                var pnorm = p.AsSpan().NormSqr();
                if (pnorm <= rmin)
                {
                    //			IsAzz=true;
                    //			 throw new TelmaException("azaz=%le",azaz);
                    break;
                }
                T arr0 = r.AsSpan().Dot(r0);
                rab.AsSpan().Assign(z);
                Factorizer?.UMult(rab);
                Matrix.MultVect(rab, az);
                Factorizer?.LMult(az);
                T ar0z = r0.AsSpan().Dot(az);
                alk = arithmetic.Devide(arr0, ar0z);
                delr -= arithmetic.Abs(arithmetic.Multiply(alk, alk)) * pnorm;
                if (delr <= rmin) delr = r.AsSpan().NormSqr();
                p.AsSpan().Assign(r);
                p.AsSpan().Add(az, arithmetic.UnaryMinus(alk));
                rab.AsSpan().Assign(p);
                Factorizer?.UMult(rab);
                Matrix.MultVect(rab, rab2);
                Factorizer?.LMult(rab2);
                T pap = p.AsSpan().Dot(rab2);
                var apap= rab2.AsSpan().NormSqr();
                gamk = arithmetic.Devide(pap, apap);
                u.Add(z, alk);
                u.Add(p, gamk);
                p.AsSpan().Add(rab2, arithmetic.UnaryMinus(gamk));
                rpred.AsSpan().Assign(r);
                r.AsSpan().Assign(p);
                T rr0;
                rr0= r.AsSpan().Dot(r0);
                betk= arithmetic.Devide(arithmetic.Multiply(rr0,alk), arithmetic.Multiply(arr0, gamk));
                z.AsSpan().Assign(r);
                z.AsSpan().Add(az, arithmetic.UnaryMinus(arithmetic.Multiply(betk, gamk)));
                z.AsSpan().Add(rpred, betk);
                
                double dd = (Math.Sqrt(delr) + rmin);
                delrvh = dd / rvhod;
                delru0 = dd / ru00;
            }
            Factorizer?.UMult(u);
#if DEBUG

            Matrix.CalcResidual(Pr, u, z);
            var rrist = z.AsSpan().Norm();
            delru0 = rrist / ru00;
#endif
            ContextManager.Context.Logger.Log(LogLevel.Information, $"SLAE solving BiCGStab:eps={parameters.Epsru0},del={delru0} end");
            return (delru0, iter);
        }
        public override SolverType DisplayCode => SolverType.BiCGStab;

        public override void UseNeumannCorrection(ReadOnlySpan<T> vector)
        {
            throw new NotSupportedException();
        }
    }
    public class DiagonalFactorizer<T> : IMatrixFactorizer<T> where T : struct, IEquatable<T>
    {
        protected SparseMatrix<T> Matrix { get; private set; }
        protected T[] cdi;
        public DiagonalFactorizer(SparseMatrix<T> slae)
        {
            Matrix = slae;
        }

        public FactorizerType DisplayCode => FactorizerType.DI;

        public void UxMult(ReadOnlySpan<T> to, Span<T> res)
        {
            // вычucляem rabnew=(U)*rabold
            res.Assign(to);
            res.MultAll(cdi);
        }

        public void Factorize(IProgressBar progress=null, CancellationToken? cancellationToken=null)
        {
            var sz = Matrix.Size;
            if (cdi != null && cdi.Length == sz) return;
            cdi = new T[sz];

            for (int i = 0; i < sz; i++)
            {
                progress?.Report(($"Factorize {i} row from {sz}",i));
                cancellationToken?.ThrowIfCancellationRequested();
                cdi[i] = Matrix.GetDi(i);
                if (cdi[i].Equals(default)) throw new NotSupportedException("Zero diagonal element");
            }
            cdi.AsSpan().SqrtAll();
        }
        public bool IsEnabledFor(ISparseMatrixBase slae) => true;
        public void LMult(Span<T> vector) //! умножение на L^-1
        {
            vector.DevideAll(cdi);
        }
        public void UMult(Span<T> vector)//! умножение на U^-1
        {
            vector.DevideAll(cdi);
        }

        #region IDisposable Support
        private bool disposedValue = false; // To detect redundant calls

        protected virtual void Dispose(bool disposing)
        {
            if (!disposedValue)
            {
                if (disposing)
                {
                    cdi = null;
                    Matrix = null;
                    // TODO: dispose managed state (managed objects).
                }

                // TODO: free unmanaged resources (unmanaged objects) and override a finalizer below.
                // TODO: set large fields to null.

                disposedValue = true;
            }
        }

        // TODO: override a finalizer only if Dispose(bool disposing) above has code to free unmanaged resources.
        // ~DiagonalFactorizer() {
        //   // Do not change this code. Put cleanup code in Dispose(bool disposing) above.
        //   Dispose(false);
        // }

        // This code added to correctly implement the disposable pattern.
        public void Dispose()
        {
            // Do not change this code. Put cleanup code in Dispose(bool disposing) above.
            Dispose(true);
            // TODO: uncomment the following line if the finalizer is overridden above.
            // GC.SuppressFinalize(this);
        }
        #endregion
    }
    public abstract class LuSQFactorizerForRowColNonSym<T> : IMatrixFactorizer<T> where T : struct, IEquatable<T>
    {
        protected RowSparseNonSymMatrix<T> Matrix { get; private set; }

        public bool CanWorkAsLUFactorizer => true;

        public bool CanUseUxMult => true;

        public FactorizerType DisplayCode => FactorizerType.LUSQ;

        protected T[] cdi;
        protected T[] cggl;
        protected T[] cggu;
        public LuSQFactorizerForRowColNonSym(RowSparseNonSymMatrix<T> _slae) { Matrix = _slae; }
        void Clear() //!Удаляет данные факторизации
        {
            cdi = null;
            cggl = null;
            cggu = null;
        }
        public void Factorize(IProgressBar progress, CancellationToken? cancellationToken)
        {
            if (cdi != null && cdi.Length == Matrix.BlockSize) return; //AllocateData обнуляет !!!
            var ggl = Matrix.Ggl;
            var ggu = Matrix.Ggu;
            cdi = new T[Matrix.BlockSize];
            cggl = new T[ggl.Length];
            cggu = new T[ggu.Length];
            // Henoлнoe LUsq paзлoжeнue  для CЛAУ с неcummeтpuчнoй maтpuцeй
            // в paзpeжeннom cтpoчнom фopmaтe
            const double rumen = 1.0e-20;
            int nerazl = 0;
            //nerazl - k-вo нepaзлoжeнныx эл-тoв CЛAУ
            T r = default;
            int nerazinstr = 0;
//            bool badstringexist = false;
            var sz = Matrix.BlockSize;
            var arithmetic = TemplateArithmetic.Create<T>();
            if (progress != null) progress.Size = sz;
            for (int i = 0; i < sz; i++)
            {
                progress?.Report(($"Factorize {i} row from {sz}",i));
                nerazinstr = 0;
                r = Matrix.GetDi(i);
                do
                {
                    T sumdi = default;
                    if (nerazinstr > 0)
                    {
                        cggl[Matrix.Ig[i] + nerazinstr - 1] = default;
                        cggu[Matrix.Ig[i] + nerazinstr - 1] = default;
  //                      badstringexist = true;
                    }
                    for (int ii = Matrix.Ig[i] + nerazinstr; ii < Matrix.Ig[i + 1]; ii++)
                    {
                        // вычucляem cggl(i,j) и cggu(j,i)
                        var j = Matrix.Jg[ii];
                        cggl[ii] = ggl.Span[ii];
                        cggu[ii] = ggu.Span[ii];
                        //  ymнoжaem i-ю блок-cтpoky нa j-ю
                        var jbeg = Matrix.Ig[j];
                        var jend = Matrix.Ig[j + 1] - 1;
                        if (jbeg <= jend)
                        {
                            var nggi = Matrix.Ig[i];
                            // nggi - нomep тekyщeгo эл-тa в maccuвe gg в i-й cтpoke
                            var nggj = jbeg;
                            // nggj - нomep тekyщeгo эл-тa в maccuвe gg в j-й cтpoke
                            while (true)
                            {
                                if (Matrix.Jg[nggi] <= Matrix.Jg[nggj])
                                {
                                    if (Matrix.Jg[nggi] == Matrix.Jg[nggj])
                                    {
                                        cggl[ii] = arithmetic.Subtruct(cggl[ii], arithmetic.Multiply(cggu[nggj], cggl[nggi]));
                                        cggu[ii] = arithmetic.Subtruct(cggu[ii], arithmetic.Multiply(cggl[nggj], cggu[nggi]));
                                    }
                                    //							if(++nggi>=ig[i+1])break;
                                    if (++nggi >= ii) break;
                                }
                                else
                                {
                                    if (++nggj > jend) break;
                                }
                            }
                        }
                        cggl[ii] = arithmetic.Devide(cggl[ii], cdi[j]);
                        cggu[ii] = arithmetic.Devide(cggu[ii], cdi[j]);
                        sumdi = arithmetic.Sum(sumdi, arithmetic.Multiply(cggl[ii], cggu[ii]));
                    }
                    r = arithmetic.Subtruct(Matrix.GetDi(i), sumdi);
                } while (arithmetic.IsLess(arithmetic.Devide(r, rumen), Matrix.GetDi(i)) && Matrix.Ig[i] + (++nerazinstr) <= Matrix.Ig[i + 1]);
                nerazl += nerazinstr;
                if (nerazl > 0) break;
                cdi[i] = arithmetic.Sqrt(r);
            }
            //    _ASSERTE(!nerazl);//message("Bad precondition: nerazl=%d",nerazl);
            if (nerazl > 0)
            {
                cggl = null;
                cggu = null;
                for (int i = 0; i < sz; i++)
                {
                    cdi[i] = Matrix.GetDi(i);
                    if (arithmetic.IsLess(cdi[i], default(T))) throw new ArgumentOutOfRangeException($"Negative di[{i}]={cdi[i]}");
                }
                cdi.AsSpan().SqrtAll();
            }
            //    if (badstringexist) TelmaStd::logonlyerrormessage(u"LUsq error");
        }

        public abstract void LMult(Span<T> vector);
        public abstract void UMult(Span<T> vector);
        public abstract void UxMult(ReadOnlySpan<T> to, Span<T> res);

        #region IDisposable Support
        private bool disposedValue = false; // To detect redundant calls

        protected virtual void Dispose(bool disposing)
        {
            if (!disposedValue)
            {
                if (disposing)
                {
                    Clear();
                    Matrix = null;
                    // TODO: dispose managed state (managed objects).
                }

                // TODO: free unmanaged resources (unmanaged objects) and override a finalizer below.
                // TODO: set large fields to null.

                disposedValue = true;
            }
        }

        // TODO: override a finalizer only if Dispose(bool disposing) above has code to free unmanaged resources.
        // ~LuSQFactorizerForRowColNonSym() {
        //   // Do not change this code. Put cleanup code in Dispose(bool disposing) above.
        //   Dispose(false);
        // }

        // This code added to correctly implement the disposable pattern.
        public void Dispose()
        {
            // Do not change this code. Put cleanup code in Dispose(bool disposing) above.
            Dispose(true);
            // TODO: uncomment the following line if the finalizer is overridden above.
            // GC.SuppressFinalize(this);
        }
        #endregion
    }
    public class LuSQFactorizerForRowColNonSymOfDouble : LuSQFactorizerForRowColNonSym<double>
    {
        public LuSQFactorizerForRowColNonSymOfDouble(RowSparseNonSymMatrix<double> _slae) : base(_slae)        {        }

        public override void UxMult(ReadOnlySpan<double> to, Span<double> res)
        {// вычucляem rabnew=(U)*rabold
            var sz = cdi.Length;
            res.Assign(to);
            res.MultAll(cdi);
            if (cggu != null && cggu.Length == Matrix.Jg.Length)
            {
                for (int ii = 0; ii < sz; ii++)
                {
                    var start = Matrix.Ig[ii];
                    var len = Matrix.Ig[ii + 1] - start;
                    MatrixAssistant.SparseAdd(res, Matrix.Jg.Slice(start, len), cggu.AsSpan().Slice(start, len), to[ii]);
                }
            }
        }
        public override void LMult(Span<double> vector)
        {
            // вычucляem r=(L)**(-1)**r
            var sz = cdi.Length;
            for (int ii = 0; ii < sz; ii++)
            {
                double sum = default;
                if (cggl != null)
                {
                    var start = Matrix.Ig[ii];
                    var len = Matrix.Ig[ii + 1] - Matrix.Ig[ii];
                    sum = MatrixAssistant.SparseMult(cggl.AsSpan().Slice(start, len), Matrix.Jg.Slice(start, len), vector);
                }
                vector[ii] = (vector[ii] - sum) / cdi[ii];
            }
        }
        public override void UMult(Span<double> vector)
        {
            var sz = cdi.Length;
            for (int ii = sz - 1; ii >= 0; ii--)
            {
                double rii = vector[ii];
                rii /= cdi[ii];
                vector[ii] = rii;
                if (cggu != null)
                {
                    var start = Matrix.Ig[ii];
                    var len = Matrix.Ig[ii + 1] - Matrix.Ig[ii];
                    MatrixAssistant.SparseAdd(vector, Matrix.Jg.Slice(start, len), cggu.AsSpan().Slice(start, len), -rii);
                }
            }
        }
    }
    public class LuSQFactorizerForRowColNonSymOfComplex : LuSQFactorizerForRowColNonSym<Complex>
    {
        public LuSQFactorizerForRowColNonSymOfComplex(RowSparseNonSymMatrix<Complex> _slae) : base(_slae)
        {
        }

        public override void UxMult(ReadOnlySpan<Complex> to, Span<Complex> res)
        {// вычucляem rabnew=(U)*rabold
            var sz = cdi.Length;
            res.Assign(to);
            res.MultAll(cdi);
            if (cggu != null && cggu.Length == Matrix.Jg.Length)
            {
                for (int ii = 0; ii < sz; ii++)
                {
                    var start = Matrix.Ig[ii];
                    var len = Matrix.Ig[ii + 1] - start;
                    MatrixAssistant.SparseAdd(res, Matrix.Jg.Slice(start, len), cggu.AsSpan().Slice(start, len), to[ii]);
                }
            }
        }
        public override void LMult(Span<Complex> vector)
        {
            // вычucляem r=(L)**(-1)**r
            var sz = cdi.Length;
            for (int ii = 0; ii < sz; ii++)
            {
                Complex sum = default;
                if (cggl != null)
                {
                    var start = Matrix.Ig[ii];
                    var len = Matrix.Ig[ii + 1] - Matrix.Ig[ii];
                    sum = MatrixAssistant.SparseMult(cggl.AsSpan().Slice(start, len), Matrix.Jg.Slice(start, len), vector);
                }
                vector[ii] = (vector[ii] - sum) / cdi[ii];
            }
        }
        public override void UMult(Span<Complex> vector)
        {
            var sz = cdi.Length;
            for (int ii = sz - 1; ii >= 0; ii--)
            {
                Complex rii = vector[ii];
                rii /= cdi[ii];
                vector[ii] = rii;
                if (cggu != null)
                {
                    var start = Matrix.Ig[ii];
                    var len = Matrix.Ig[ii + 1] - Matrix.Ig[ii];
                    MatrixAssistant.SparseAdd(vector, Matrix.Jg.Slice(start, len), cggu.AsSpan().Slice(start, len), -rii);
                }
            }
        }
    }
    public abstract class LuSQFactorizerForRowColSym<T> : IMatrixFactorizer<T> where T : struct
    {
        protected RowSparseSymmetricMatrix<T> Matrix { get; private set; }

        public bool CanWorkAsLUFactorizer => true;

        public bool CanUseUxMult => true;

        public FactorizerType DisplayCode => FactorizerType.LUSQ;

        protected T[] cdi;
        protected T[] cgg;
        public LuSQFactorizerForRowColSym(RowSparseSymmetricMatrix<T> _slae) { Matrix = _slae; }
        void Clear() //!Удаляет данные факторизации
        {
            cdi = null;
            cgg = null;
        }
        public void Factorize(IProgressBar progress, CancellationToken? cancellationToken)
        {
            if (cdi != null) return; //AllocateData обнуляет !!!
            Matrix.CorrectDiagonal();
            cdi = new T[Matrix.BlockSize];
            cgg = new T[Matrix.Gg.Length];
            // Henoлнoe paзлoжeнue Xoлecckoгo для CЛAУ cummeтpuчнoй maтpuцeй
            // в paзpeжeннom cтpoчнom фopmaтe (maтpuцa зaдaнa cвoum нuжнum тpeyг-kom)
            const double rumen = 1e-5;
            int nerazl = 0;
//            bool badstringexist = false;
            //nerazl - k-вo нepaзлoжeнныx эл-тoв CЛAУ
            T r = default;
            int nerazinstr = 0;
            var sz = Matrix.BlockSize;
            var arithmetic = TemplateArithmetic.Create<T>();
            for (int ii = 0; ii < sz; ii++)
            {
                nerazinstr = 0;
                r = Matrix.GetDi(ii);
                do
                {
                    T sumdi = default;
                    if (nerazinstr > 0)
                    {
                        cgg[Matrix.Ig[ii] + nerazinstr - 1] = default;
//                        badstringexist = true;
                    }
                    for (int i = Matrix.Ig[ii] + nerazinstr; i < Matrix.Ig[ii + 1]; i++)
                    {
                        // вычucляem cgg(i)
                        T sum = default;
                        var jj = Matrix.Jg[i];
                        if (i != Matrix.Ig[ii] + nerazinstr)
                        {
                            var najj = Matrix.Ig[jj];
                            var kojj = Matrix.Ig[jj + 1] - 1;
                            if (najj <= kojj)
                            {
                                //  ymнoжaem ii-ю cтpoky нa jj-ю
                                var nggi = Matrix.Ig[ii] + nerazinstr;
                                // nggi - нomep тekyщeгo эл-тa в maccuвe gg в ii-й cтpoke
                                var nggj = najj;
                                // nggj - нomep тekyщeгo эл-тa в maccuвe gg в jj-й cтpoke
                                while (true)
                                {
                                    if (Matrix.Jg[nggi] <= Matrix.Jg[nggj])
                                    {
                                        if (Matrix.Jg[nggi] == Matrix.Jg[nggj]) sum = arithmetic.Sum(sum, arithmetic.Multiply(cgg[nggi], cgg[nggj]));
                                        if (++nggi >= i) break;
                                    }
                                    else
                                    {
                                        if (++nggj > kojj) break;
                                    }
                                }
                            }
                        }
                        cgg[i] = arithmetic.Devide(arithmetic.Subtruct(Matrix.Gg[i], sum), cdi[jj]);
                    }
                    //  ymнoжaem ii-ю cтpoky нa ii-ю

                    sumdi = arithmetic.Sum(sumdi, cgg.AsSpan().Slice(Matrix.Ig[ii] + nerazinstr, Matrix.Ig[ii + 1] - Matrix.Ig[ii] - nerazinstr).NormSqr());
                    //for(i=ig[ii]+nerazinstr;i<ig[ii+1];i++)
                    //{
                    //	long double rr=cgg[i];
                    //	sumdi+=rr*rr;
                    //}
                    r = arithmetic.Subtruct(Matrix.Di[ii], sumdi);
                }
                while (arithmetic.IsLess(arithmetic.Devide(r, rumen), Matrix.Di[ii]) && Matrix.Ig[ii] + (++nerazinstr) <= Matrix.Ig[ii + 1]);
                nerazl += nerazinstr;
                cdi[ii] = arithmetic.Sqrt(r);
            }
            //		_ASSERTE(!nerazl);//message("Bad precondition: nerazl=%d",nerazl);
            if (nerazl > 0)
            {
                cgg = null;
                for (int i = 0; i < sz; i++)
                {
                    if (arithmetic.IsLess(Matrix.GetDi(i), default)) throw new ArgumentOutOfRangeException($"Negative di[{i}]={cdi[i]}");
                    cdi[i] = arithmetic.Sqrt(Matrix.GetDi(i));
                }
            }
            //    if (badstringexist) TelmaStd::logonlyerrormessage(u"LUsq error");
        }
        public abstract void LMult(Span<T> vector);
        public abstract void UMult(Span<T> vector);
        public abstract void UxMult(ReadOnlySpan<T> to, Span<T> res);

        #region IDisposable Support
        private bool disposedValue = false; // To detect redundant calls

        protected virtual void Dispose(bool disposing)
        {
            if (!disposedValue)
            {
                if (disposing)
                {
                    Clear();
                    Matrix = null;
                    // TODO: dispose managed state (managed objects).
                }

                // TODO: free unmanaged resources (unmanaged objects) and override a finalizer below.
                // TODO: set large fields to null.

                disposedValue = true;
            }
        }

        // TODO: override a finalizer only if Dispose(bool disposing) above has code to free unmanaged resources.
        // ~LuSQFactorizerForRowColSym() {
        //   // Do not change this code. Put cleanup code in Dispose(bool disposing) above.
        //   Dispose(false);
        // }

        // This code added to correctly implement the disposable pattern.
        public void Dispose()
        {
            // Do not change this code. Put cleanup code in Dispose(bool disposing) above.
            Dispose(true);
            // TODO: uncomment the following line if the finalizer is overridden above.
            // GC.SuppressFinalize(this);
        }
        #endregion
    }
    public class LuSQFactorizerForRowColSymOfDouble : LuSQFactorizerForRowColSym<double>
    {
        public LuSQFactorizerForRowColSymOfDouble(RowSparseSymmetricMatrix<double> _slae) : base(_slae)
        {
        }

        public override void LMult(Span<double> f)
        {
            // вычucляem r=(L)**(-1)**r
            var siz = Matrix.BlockSize;

            for (int ii = 0; ii < siz; ii++)
            {
                double sum = 0;

                if (cgg != null)
                {
                    var start = Matrix.Ig[ii];
                    var len = Matrix.Ig[ii + 1] - Matrix.Ig[ii];
                    sum = MatrixAssistant.SparseMult(cgg.AsSpan().Slice(start, len), Matrix.Jg.Slice(start, len), f);
                }
                f[ii] = (f[ii] - sum) / cdi[ii];
            }
        }
        public override void UMult(Span<double> rab)
        {
            // вычucляem rabnew=(L')**(-1)*rabold (т.e. peшaem L'*rabnew=rabold)
            var siz = Matrix.BlockSize;
            for (var ii = siz - 1; ii >= 0; ii--)
            {
                var rabii = rab[ii] / cdi[ii];
                rab[ii] = rabii;
                if (cgg != null)
                {
                    var start = Matrix.Ig[ii];
                    var len = Matrix.Ig[ii + 1] - Matrix.Ig[ii];
                    MatrixAssistant.SparseAdd(rab, Matrix.Jg.Slice(start, len), cgg.AsSpan().Slice(start, len), -rabii);
                }
            }
        }
        public override void UxMult(ReadOnlySpan<double> to, Span<double> res)
        {
            // вычucляem rabnew=(U)*rabold
            var sz = Matrix.BlockSize;
            res.Assign(to);
            res.MultAll(cdi);
            for (int ii = 0; ii < sz; ii++)
            {
                for (int jj = Matrix.Ig[ii]; jj < Matrix.Ig[ii + 1]; jj++)
                {
                    res[Matrix.Jg[jj]] += cgg[jj] * to[ii];
                }
            }
        }
    }
    public class LuSQFactorizerForRowColSymOfComplex : LuSQFactorizerForRowColSym<Complex>
    {
        public LuSQFactorizerForRowColSymOfComplex(RowSparseSymmetricMatrix<Complex> _slae) : base(_slae)
        {
        }

        public override void LMult(Span<Complex> f)
        {
            // вычucляem r=(L)**(-1)**r
            var siz = Matrix.BlockSize;

            for (int ii = 0; ii < siz; ii++)
            {
                Complex sum = 0;

                if (cgg != null)
                {
                    var start = Matrix.Ig[ii];
                    var len = Matrix.Ig[ii + 1] - Matrix.Ig[ii];
                    sum = MatrixAssistant.SparseMult(cgg.AsSpan().Slice(start, len), Matrix.Jg.Slice(start, len), f);
                }
                f[ii] = (f[ii] - sum) / cdi[ii];
            }
        }
        public override void UMult(Span<Complex> rab)
        {
            // вычucляem rabnew=(L')**(-1)*rabold (т.e. peшaem L'*rabnew=rabold)
            var siz = Matrix.BlockSize;
            for (var ii = siz - 1; ii >= 0; ii--)
            {
                var rabii = rab[ii] / cdi[ii];
                rab[ii] = rabii;
                if (cgg != null)
                {
                    var start = Matrix.Ig[ii];
                    var len = Matrix.Ig[ii + 1] - Matrix.Ig[ii];
                    MatrixAssistant.SparseAdd(rab, Matrix.Jg.Slice(start, len), cgg.AsSpan().Slice(start, len), -rabii);
                }
            }
        }
        public override void UxMult(ReadOnlySpan<Complex> to, Span<Complex> res)
        {
            // вычucляem rabnew=(U)*rabold
            var sz = Matrix.BlockSize;
            res.Assign(to);
            res.MultAll(cdi);
            for (int ii = 0; ii < sz; ii++)
            {
                for (int jj = Matrix.Ig[ii]; jj < Matrix.Ig[ii + 1]; jj++)
                {
                    res[Matrix.Jg[jj]] += cgg[jj] * to[ii];
                }
            }
        }
    }
    public abstract class LuFactorizerForRowColNonSym<T> : IMatrixFactorizer<T> where T : struct
    {
        protected RowSparseNonSymMatrix<T> Matrix { get; private set; }

        public bool CanWorkAsLUFactorizer => true;

        public bool CanUseUxMult => true;

        public FactorizerType DisplayCode => FactorizerType.LUSQ;

        protected T[] cdi;
        protected T[] cggl;
        protected T[] cggu;
        public LuFactorizerForRowColNonSym(RowSparseNonSymMatrix<T> _slae) { Matrix = _slae; }
        void Clear() //!Удаляет данные факторизации
        {
            cdi = null;
            cggl = null;
            cggu = null;
        }
        public void Factorize(IProgressBar progress, CancellationToken? cancellationToken)
        {
            if (cdi != null && cdi.Length == Matrix.BlockSize) return; //AllocateData обнуляет !!!
            var ggl = Matrix.Ggl;
            var ggu = Matrix.Ggu;
            cdi = new T[Matrix.BlockSize];
            cggl = new T[ggl.Length];
            cggu = new T[ggu.Length];
            // Henoлнoe LUsq paзлoжeнue  для CЛAУ с неcummeтpuчнoй maтpuцeй
            // в paзpeжeннom cтpoчнom фopmaтe
            //const double rumen = 1.0e-20;
            //int nerazl = 0;
            //nerazl - k-вo нepaзлoжeнныx эл-тoв CЛAУ
            T r = default;
            var sz = Matrix.BlockSize;
            var arithmetic = TemplateArithmetic.Create<T>();
            if (progress != null) progress.Size = sz;
            for (int i = 0; i < sz; i++)
            {
                progress?.Report(($"LU factorization {i} row from {sz}", i));
                r = Matrix.GetDi(i);
                T sumdi = default;
                for (int ii = Matrix.Ig[i]; ii < Matrix.Ig[i + 1]; ii++)
                {
                    // вычucляem cggl(i,j) и cggu(j,i)
                    var j = Matrix.Jg[ii];
                    cggl[ii] = ggl.Span[ii];
                    cggu[ii] = ggu.Span[ii];
                    //  ymнoжaem i-ю блок-cтpoky нa j-ю
                    var jbeg = Matrix.Ig[j];
                    var jend = Matrix.Ig[j + 1] - 1;
                    if (jbeg <= jend)
                    {
                        var nggi = Matrix.Ig[i];
                        // nggi - нomep тekyщeгo эл-тa в maccuвe gg в i-й cтpoke
                        var nggj = jbeg;
                        // nggj - нomep тekyщeгo эл-тa в maccuвe gg в j-й cтpoke
                        while (true)
                        {
                            if (Matrix.Jg[nggi] <= Matrix.Jg[nggj])
                            {
                                if (Matrix.Jg[nggi] == Matrix.Jg[nggj])
                                {
                                    cggl[ii] = arithmetic.Subtruct(cggl[ii], arithmetic.Multiply(cggu[nggj], cggl[nggi]));
                                    cggu[ii] = arithmetic.Subtruct(cggu[ii], arithmetic.Multiply(cggl[nggj], cggu[nggi]));
                                }
                                //							if(++nggi>=ig[i+1])break;
                                if (++nggi >= ii) break;
                            }
                            else
                            {
                                if (++nggj > jend) break;
                            }
                        }
                    }
                    cggl[ii] = arithmetic.Devide(cggl[ii], cdi[j]);
                    sumdi = arithmetic.Sum(sumdi, arithmetic.Multiply(cggl[ii], cggu[ii]));
                }
                r = arithmetic.Subtruct(Matrix.GetDi(i), sumdi);
                cdi[i] = r;
            }
            //    _ASSERTE(!nerazl);//message("Bad precondition: nerazl=%d",nerazl);
            //    if (badstringexist) TelmaStd::logonlyerrormessage(u"LUsq error");
        }

        public abstract void LMult(Span<T> vector);
        public abstract void UMult(Span<T> vector);
        public abstract void UxMult(ReadOnlySpan<T> to, Span<T> res);

        #region IDisposable Support
        private bool disposedValue = false; // To detect redundant calls

        protected virtual void Dispose(bool disposing)
        {
            if (!disposedValue)
            {
                if (disposing)
                {
                    Clear();
                    Matrix = null;
                    // TODO: dispose managed state (managed objects).
                }

                // TODO: free unmanaged resources (unmanaged objects) and override a finalizer below.
                // TODO: set large fields to null.

                disposedValue = true;
            }
        }

        // TODO: override a finalizer only if Dispose(bool disposing) above has code to free unmanaged resources.
        // ~LuFactorizerForRowColNonSym() {
        //   // Do not change this code. Put cleanup code in Dispose(bool disposing) above.
        //   Dispose(false);
        // }

        // This code added to correctly implement the disposable pattern.
        public void Dispose()
        {
            // Do not change this code. Put cleanup code in Dispose(bool disposing) above.
            Dispose(true);
            // TODO: uncomment the following line if the finalizer is overridden above.
            // GC.SuppressFinalize(this);
        }
        #endregion
    }
    public class LuFactorizerForRowColNonSymOfDouble : LuFactorizerForRowColNonSym<double>
    {
        public LuFactorizerForRowColNonSymOfDouble(RowSparseNonSymMatrix<double> _slae) : base(_slae)
        {
        }

        public override void UxMult(ReadOnlySpan<double> to, Span<double> res)
        {
            res.Assign(to);
            res.MultAll(cdi);
            if (cggu != null)
            {
                var sz = Matrix.Size;
                for (int ii = 0; ii < sz; ii++)
                {
                    var start = Matrix.Ig[ii];
                    var len = Matrix.Ig[ii + 1] - start;
                    MatrixAssistant.SparseAdd(res, Matrix.Jg.Slice(start, len), cggu.AsSpan().Slice(start, len), to[ii]);
                }
            }
        }
        public override void LMult(Span<double> rab)
        {
            // вычucляem r=(L)**(-1)**r
            var siz = Matrix.BlockSize;
            for (int ii = 0; ii < siz; ii++)
            {
                double sum = 0;
                if (cggl != null)
                {
                    var start = Matrix.Ig[ii];
                    var len = Matrix.Ig[ii + 1] - start;
                    sum = MatrixAssistant.SparseMult(cggl.AsSpan().Slice(start, len), Matrix.Jg.Slice(start, len), rab);
                }
                rab[ii] = (rab[ii] - sum);
            }
        }
        public override void UMult(Span<double> rab)
        {
            var siz = Matrix.BlockSize;
            for (var ii = siz - 1; ii >= 0; ii--)
            {
                var rabii = rab[ii];
                rabii /= cdi[ii];
                rab[ii] = rabii;
                if (cggu != null)
                {
                    var start = Matrix.Ig[ii];
                    var len = Matrix.Ig[ii + 1] - start;
                    MatrixAssistant.SparseAdd(rab, Matrix.Jg.Slice(start, len), cggu.AsSpan().Slice(start, len), -rabii);
                }
            }
        }
    }
    public class LuFactorizerForRowColNonSymOfComplex : LuFactorizerForRowColNonSym<Complex>
    {
        public LuFactorizerForRowColNonSymOfComplex(RowSparseNonSymMatrix<Complex> _slae) : base(_slae)
        {
        }

        public override void UxMult(ReadOnlySpan<Complex> to, Span<Complex> res)
        {
            res.Assign(to);
            res.MultAll(cdi);
            if (cggu != null)
            {
                var sz = Matrix.Size;
                for (int ii = 0; ii < sz; ii++)
                {
                    var start = Matrix.Ig[ii];
                    var len = Matrix.Ig[ii + 1] - start;
                    MatrixAssistant.SparseAdd(res, Matrix.Jg.Slice(start, len), cggu.AsSpan().Slice(start, len), to[ii]);
                }
            }
        }
        public override void LMult(Span<Complex> rab)
        {
            // вычucляem r=(L)**(-1)**r
            var siz = Matrix.BlockSize;
            for (int ii = 0; ii < siz; ii++)
            {
                Complex sum = 0;
                if (cggl != null)
                {
                    var start = Matrix.Ig[ii];
                    var len = Matrix.Ig[ii + 1] - start;
                    sum = MatrixAssistant.SparseMult(cggl.AsSpan().Slice(start, len), Matrix.Jg.Slice(start, len), rab);
                }
                rab[ii] = (rab[ii] - sum);
            }
        }
        public override void UMult(Span<Complex> rab)
        {
            var siz = Matrix.BlockSize;
            for (var ii = siz - 1; ii >= 0; ii--)
            {
                var rabii = rab[ii];
                rabii /= cdi[ii];
                rab[ii] = rabii;
                if (cggu != null)
                {
                    var start = Matrix.Ig[ii];
                    var len = Matrix.Ig[ii + 1] - start;
                    MatrixAssistant.SparseAdd(rab, Matrix.Jg.Slice(start, len), cggu.AsSpan().Slice(start, len), -rabii);
                }
            }
        }
    }
    public class GMRESSolver : SlaeSolver<double>
    {
        public int SolverDepth { get; set; }
        double[][] HHt;
        double[][] DepthVector;
        double[] Hessenberg;
        double[] D;
        double[] Z;
        protected override void Dispose(bool disposing)
        {
            base.Dispose(disposing);
            if (disposing)
            {
                HHt = null;
                DepthVector = null;
                Hessenberg = null;
                D = null;
                Z = null;
            }
        }
        void HessenbergSolve()// построение базиса пространства
        {
            int sz = Matrix.Size;
            int depth = SolverDepth;
            if (HHt == null || HHt.Length != depth) HHt = new double[depth][];
            for (int i = 0; i < depth; i++)
            {
                if (HHt[i] == null) HHt[i] = new double[depth+1];
                else for (int j = 0; j < depth; j++) HHt[i][j] = 0;
            }
            for (int i = 0; i < depth; i++)
            {
                for (int j = 0; j < depth; j++)
                {
                    for (int k = 0; k <= depth; k++)
                    {
                        HHt[i][j] += Hessenberg[k * (depth + 1) + i] * Hessenberg[k * (depth + 1) + j];
                    }
                }
                for (int k = 0; k < depth; k++)
                    HHt[i][depth] += Hessenberg[k * (depth + 1) + i] * D[k];
            }
            //Решение СЛАУ методом Гаусса c выбором по столбцу
            for (int i = 0; i < depth; i++)
            {
                int maxdiag = i;
                for (int j = i + 1; j < depth; j++) if (Math.Abs(HHt[j][i]) > Math.Abs(HHt[maxdiag][i])) maxdiag = j;
                if (maxdiag != i)
                {
                    var hh = HHt[i];
                    HHt[i] = HHt[maxdiag];
                    HHt[maxdiag] = hh;
                }
                for (int j = i + 1; j <= depth; j++) HHt[i][j] /= HHt[i][i];
                for (int j = i + 1; j < depth; j++)
                {
                    for (int k = i + 1; k <= depth; k++) HHt[j][k] -= HHt[i][k] * HHt[j][i];
                }
            }
            if (Z == null || Z.Length != depth) Z = new double[depth];
            for (int i = depth - 1; i >= 0; i--)
            {
                Z[i] = HHt[i][depth];
                for (int j = i + 1; j < depth; j++) Z[i] -= HHt[i][j] * Z[j];
            }
        }
        protected override (double residual, int iterations) SolveMethod(ReadOnlySpan<double> Pr, Span<double> u, ISlaeSolverParameters parameters, IProgressBar progress, CancellationToken? ct)
        {
            /*GMRES c нenoлныm LU(sq) paзлoжeнuem для CЛAY c неcummeтpuчнoй maтpuцeй в paзpeжeннom cтpoчнom фopmaтe*/
            var size = Matrix.Size;
            if (SolverDepth == 0) SolverDepth = 5;
             var depth = SolverDepth;
            if (DepthVector == null || DepthVector.Length != depth) DepthVector = new double[depth + 1][];
            if (Hessenberg == null || Hessenberg.Length != (depth + 1) * (depth + 1)) Hessenberg = new double[(depth + 1) * (depth + 1)];
            if (D == null || D.Length != depth) D = new double[depth];
            //HHt.resize(this->depth);
            //Z.resize(this->depth);
            if (u == null || u.Length != size) u = new double[size];
            if (Factorizer != null) Factorizer.Factorize(progress, ct);
            var r = GetWorkVector("r");
            var z = GetWorkVector("z");
            var az = GetWorkVector("az");
            var rab = GetWorkVector("rab");
            var rab2 = GetWorkVector("rab2");
            const double rmin = 1e-100;
            // вычисляем первые к.у.
            // Bычucляem нeвязky r, conpяжeннoe нanpaвлeнue z=r, rvhod, ru00 u delrnl
            double rvhodist = 0;
            double ru00ist = 0;
            Matrix.CalcResidual(Pr,u, r);
            rvhodist = (Math.Sqrt(rvhodist) + rmin) / size;
            ru00ist = (Math.Sqrt(ru00ist) + rmin) / size;
            //	long double delrnlist=rvhodist/ru00ist;
            double delr = rvhodist;
            Factorizer.LMult(r);
            z.AsSpan().Assign(Pr);
            Factorizer.LMult(z);
            double rvhod = r.AsSpan().NormSqr();
            double ru00 = z.AsSpan().NormSqr();
            rvhod = (Math.Sqrt(rvhod) + rmin) / size;
            ru00 = (Math.Sqrt(ru00) + rmin) / size;
            //	 long double delrnl=rvhod/ru00;
            delr = rvhod;
            rab.AsSpan().Assign(u);
            Factorizer.UxMult(rab, u);
            // Haчuнaem uтepaцuu
            double delrvh = delr / rvhod;
            double delru0 = delr / ru00;
            bool IsAzz = false;
            var resres = r.AsSpan().NormSqr();
            delr = resres;
            int iter = 0;
            if(progress != null)progress.Size = ProgressLength(parameters.Epsru0);
            if (delr > rmin)
                for (iter = 0; iter < parameters.Maxiter && delru0 > parameters.Epsru0 && delrvh > parameters.Epsrvh; iter++)
                {
                    progress?.Report(($"SLAE solving GMRES({depth}):eps={parameters.Epsru0},del={delru0},iter={iter}", ProgressLength(delru0)));

                    var normnev = Math.Sqrt(resres);
                    if (DepthVector[0] == null || DepthVector[0].Length != r.Length) DepthVector[0] = new double[r.Length];
                    DepthVector[0].AsSpan().Assign(r, 1 / normnev);
                    for (int IterDepth = 0; IterDepth < depth; IterDepth++)
                    {
                        rab.AsSpan().Assign(DepthVector[IterDepth]);
                        Factorizer.UMult(rab);
                        Matrix.MultVect(rab, rab2);
                        Factorizer.LMult(rab2);
                        if (DepthVector[IterDepth + 1] == null || DepthVector[IterDepth + 1].Length != size) DepthVector[IterDepth + 1] = new double[size];
                        DepthVector[IterDepth + 1].AsSpan().Assign(rab2);
                        for (int l = 0; l <= IterDepth; l++)
                        {
                            var hss = DepthVector[l].AsSpan().Dot(rab2);
                            Hessenberg[l * (depth + 1) + IterDepth] = hss;
                            DepthVector[IterDepth + 1].AsSpan().Add(DepthVector[l], -hss);
                        }
                        //double hs=sqrt(Scal(DepthVector[IterDepth+1],DepthVector[IterDepth+1]));
                        var hs = DepthVector[IterDepth + 1].AsSpan().Norm();
                        if (hs < rmin)
                        {
                            IsAzz = true;
                            break;
                        }
                        DepthVector[IterDepth + 1].AsSpan().DevideAll(hs);
                        Hessenberg[(IterDepth + 1) * (depth + 1) + IterDepth] = hs;
                    }
                    if (IsAzz) break;
                    for (int i = 0; i < D.Length; i++) D[i] = 0;
                    D[0] = normnev;
                    HessenbergSolve();
                    for (int kk = 0; kk < size; kk++)
                    {
                        for (int l = 0; l < depth; l++) u[kk] += DepthVector[l][kk] * Z[l];
                    }
                    rab.AsSpan().Assign(u);
                    Factorizer.UMult(rab);
                    Matrix.CalcResidual(Pr,rab, r);
                    delr = r.AsSpan().Norm();
                    Factorizer.LMult(r);
                    resres = r.AsSpan().NormSqr();
                    double dd = (delr + rmin) / size;
                    ru00 = (Pr.Norm() + rmin) / size;
                    delrvh = dd / rvhod;
                    delru0 = dd / ru00;
                }
            Factorizer.UMult(u);
            var rrist = delru0 * ru00ist;
#if DEBUG
            Matrix.CalcResidual(Pr,u, r);
            rrist = r.AsSpan().Norm() / size;
            if (!(rrist / ru00ist < 100000 * parameters.Epsru0 || rrist < 1e-15)) throw new Exception("Assertion failed");
#endif
            return (rrist / ru00ist, iter);
        }
        public override SolverType DisplayCode => SolverType.GMRES;

        public override void UseNeumannCorrection(ReadOnlySpan<double> vector)
        {
            throw new NotSupportedException();
        }
    }
}