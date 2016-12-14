using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Diagnostics;


namespace CMAESharp
{
    [Serializable]
    public class CMAEvolutionStrategy
    {

        private static readonly long serialVersionUID = 2918241407634253526L;
        public readonly string versionNumber = "0.99.40";

        //stop condition
        int index = 0;
        string[] messages = new string[0];
        double lastcounteval;

        public bool isTrue()
        {
            return test() > 0;
        }

        public bool isFalse()
        {
            return !isTrue();
        }

        public int getNumber()
        {
            return test();
        }

        public string[] getMessages()
        {
            return messages;
        }

        public void clear()
        {
            messages = new string[0];
            index = 0;
        }

        private void appendMessage(string s)
        {
            string[] mold = messages;
            messages = new string[index + 1];

            for (int i = 0; i < index; i++)
            {
                messages[i] = mold[i];
            }
            messages[index++] = s + " (iter=" + countiter + ",eval=" + counteval + ")";
        }

        int test()
        {
            if (state < 0)
                return 0;  // not yet initialized
            if (index > 0 && (counteval == lastcounteval
                                || counteval == lastcounteval + 1)) // one evaluation for xmean is ignored
                return index;  // termination criterion already met

            lastcounteval = counteval;

            /* FUNCTION VALUE */
            //IRWAN for maximize+minimize
            if (taskType == "Maximize")
            {
                if ((countiter > 1 || state >= 3) && bestever_fit >= options.stopFitness)
                    appendMessage("Fitness: Objective function value dropped below the target function value " +
                            options.stopFitness);
            }
            else if (taskType == "Minimize")
            {
                if ((countiter > 1 || state >= 3) && bestever_fit <= options.stopFitness)
                    appendMessage("Fitness: Objective function value dropped below the target function value " +
                            options.stopFitness);
            }

            /* #Fevals */
            if (counteval >= options.stopMaxFunEvals)
                appendMessage("MaxFunEvals: maximum number of function evaluations " + options.stopMaxFunEvals + " reached");

            /* #iterations */
            if (countiter >= options.stopMaxIter)
                appendMessage("MaxIter: maximum number of iterations reached");

            /* TOLFUN */
            if ((countiter > 1 || state >= 3)
                && Math.Max(math.max(fit.history), fit.fitness[fit.fitness.Length - 1].val)
                        - Math.Min(math.min(fit.history), fit.fitness[0].val) <= options.stopTolFun)
                appendMessage("TolFun: function value changes below stopTolFun=" + options.stopTolFun);

            /* TOLFUNHIST */
            if (options.stopTolFunHist >= 0 && countiter > fit.history.Length)
            {
                if (math.max(fit.history) - math.min(fit.history) <= options.stopTolFunHist)
                    appendMessage("TolFunHist: history of function value changes below stopTolFunHist=" + options.stopTolFunHist);
            }

            /* TOLX */
            double tolx = Math.Max(options.stopTolX, options.stopTolXFactor * minstartsigma);
            if (sigma * maxsqrtdiagC < tolx
                    && sigma * math.max(math.abs(pc)) < tolx)
                appendMessage("TolX or TolXFactor: standard deviation below " + tolx);

            /* TOLXUP */
            if (sigma * maxsqrtdiagC > options.stopTolUpXFactor * maxstartsigma)
                appendMessage("TolUpX: standard deviation increased by more than stopTolUpXFactor=" +
                        options.stopTolUpXFactor +
                        ", larger initial standard deviation recommended");

            /* STOPNOW */
            if (options.stopnow)
                appendMessage("Manual: flag Options.stopnow set or stop now in .properties file");

            /* Internal (numerical) stopping termination criteria */

            /* Test each principal axis i, whether x == x + 0.1 * sigma * rgD[i] * B[i] */
            for (int iAchse = 0; iAchse < N; ++iAchse)
            {
                int iKoo;
                int l = flgdiag ? iAchse : 0;
                int u = flgdiag ? iAchse + 1 : N;
                double fac = 0.1 * sigma * diagD[iAchse];
                for (iKoo = l; iKoo < u; ++iKoo)
                {
                    if (xmean[iKoo] != xmean[iKoo] + fac * B[iKoo][iAchse])
                        break; // is OK for this iAchse
                }
                if (iKoo == u) // no break, therefore no change for axis iAchse
                    appendMessage("NoEffectAxis: Mutation " + 0.1 * sigma * diagD[iAchse] +
                            " in a principal axis " + iAchse + " has no effect");
            } /* for iAchse */

            /* Test whether one component of xmean is stuck */
            for (int iKoo = 0; iKoo < N; ++iKoo)
            {
                if (xmean[iKoo] == xmean[iKoo] + 0.2 * sigma * Math.Sqrt(C[iKoo][iKoo]))
                    appendMessage("NoEffectCoordinate: Mutation of size " +
                            0.2 * sigma * Math.Sqrt(C[iKoo][iKoo]) +
                            " in coordinate " + iKoo + " has no effect");
            } /* for iKoo */

            /* Condition number */
            if (math.min(diagD) <= 0)
                appendMessage("ConditionNumber: smallest eigenvalue smaller or equal zero");
            else if (math.max(diagD) / math.min(diagD) > 1e7)
                appendMessage("ConditionNumber: condition number of the covariance matrix exceeds 1e14");
            return index; // call to appendMessage increments index
        }

        void testAndCorrectNumerics()
        { // not much left here

            /* Flat Fitness, Test if function values are identical */
            if (getCountIter() > 1 || (getCountIter() == 1 && state >= 3))
                if (fit.fitness[0].val == fit.fitness[Math.Min(sp.getLambda() - 1, sp.getLambda() / 2 + 1) - 1].val)
                {
                    warning("flat fitness landscape, consider reformulation of fitness, step-size increased");
                    sigma *= Math.Exp(0.2 + sp.getCs() / sp.getDamps());
                }

            /* Align (renormalize) scale C (and consequently sigma) */
            /* e.g. for infinite stationary state simulations (noise
             * handling needs to be introduced for that) */
            double fac = 1.0;
            if (math.max(diagD) < 1e-6)
                fac = 1.0 / math.max(diagD);
            else if (math.min(diagD) > 1e4)
                fac = 1.0 / math.min(diagD);

            if (fac != 1.0)
            {
                sigma /= fac;
                for (int i = 0; i < N; ++i)
                {
                    pc[i] *= fac;
                    diagD[i] *= fac;
                    for (int j = 0; j <= i; ++j)
                        C[i][j] *= fac * fac;
                }
            }
        }

        public CMAOptions options = new CMAOptions();

        private CMAParameters sp = new CMAParameters();
        public CMAParameters parameters
        {
            get { return sp; }
            set { sp = value; }
        }

        int N;
        long seed = CurrentTimeMillis();
        Random rand = new Random((int)CurrentTimeMillis());

        readonly MyMath math = new MyMath();
        double axisratio;
        long counteval;
        long countiter;

        long bestever_eval;
        double[] bestever_x;
        double bestever_fit = double.NaN;

        double sigma = 0.0;
        double[] typicalX;
        double[] initialX;
        public double[] LBound;
        public double[] UBound;
        double[] xmean;
        double xmean_fit = double.NaN;
        double[] pc;
        double[] ps;
        double[][] C;
        double maxsqrtdiagC;
        double minsqrtdiagC;
        double[][] B;
        double[] diagD;
        bool flgdiag;

        double[] startsigma;
        double maxstartsigma;
        double minstartsigma;

        bool iniphase;

        double state = -1;
        long citerlastwritten = 0;
        long countwritten = 0;
        int lockDimension = 0;
        int mode = 0;
        readonly int SINGLE_MODE = 1;
        readonly int PARALLEL_MODE = 2;

        long countCupdatesSinceEigenupdate;

        readonly FitnessCollector fit = new FitnessCollector();

        double recentFunctionValue;
        double recentMaxFunctionValue;
        double recentMinFunctionValue;
        int idxRecentOffspring;

        double[][] arx;
        public double[][] population;
        double[] xold;

        double[] BDz;
        double[] artmp;

        string propertiesFileName = "CMAEvolutionStrategy.properties";

        public CMAEvolutionStrategy(string TaskType)
        {
            state = -1;
            this.taskType = TaskType;
        }

        //not implemented properties

        public CMAEvolutionStrategy(int dimension)
        {
            setDimension(dimension);
            state = -1;
        }

        public double[] init(int dimension, double[] initialX, double[] initialStandardDeviations)
        {
            setInitialX(initialX);
            setInitialStandardDeviations(initialStandardDeviations);
            return init(dimension);
        }

        private double[] getArrayOf(double x, int dim)
        {
            double[] res = new double[dim];
            for (int i = 0; i < dim; ++i)
                res[i] = x;
            return res;
        }

        private double[] expandToDimension(double[] x, int dim)
        {
            if (x == null)
                return null;
            if (x.Length == dim)
                return x;
            if (x.Length != 1)
                error("x must have length one or length dimension");

            return getArrayOf(x[0], dim);
        }

        public double[] init(int dimension)
        {
            setDimension(dimension);
            return init();
        }

        public double[] init()
        {
            int i;
            if (N <= 0)
                error("dimension needs to be determined, use eg. setDimension() or setInitialX()");
            if (state >= 0)
                error("init() cannot be called twice");
            if (state == 0) // less save variant 
                return new double[sp.getLambda()];
            if (state > 0)
                error("init() cannot be called after the first population was sampled");

            sp = parameters; /* just in case the user assigned parameters */
            if (sp.supplemented == 0) // a bit a hack
                sp.supplementRemainders(N, options);
            sp.locked = 1; // lambda cannot be changed anymore

            diagD = new double[N];
            for (i = 0; i < N; ++i)
                diagD[i] = 1;

            /* expand Boundaries */
            LBound = expandToDimension(LBound, N);
            if (LBound == null)
            {
                LBound = new double[N];
                for (i = 0; i < N; ++i)
                    LBound[i] = Double.NegativeInfinity;
            }

            UBound = expandToDimension(UBound, N);
            if (UBound == null)
            {
                UBound = new double[N];
                for (i = 0; i < N; ++i)
                    UBound[i] = Double.PositiveInfinity;
            }

            /* Initialization of sigmas */
            if (startsigma != null)
            { // 
                if (startsigma.Length == 1)
                {
                    sigma = startsigma[0];
                }
                else if (startsigma.Length == N)
                {
                    sigma = math.max(startsigma);
                    if (sigma <= 0)
                        error("initial standard deviation sigma must be positive");
                    for (i = 0; i < N; ++i)
                    {
                        diagD[i] = startsigma[i] / sigma;
                    }
                }
                else
                    Debug.Assert(false);
            }
            else
            {
                // we might use boundaries here to find startsigma, but I prefer to have stddevs mandatory 
                error("no initial standard deviation specified, use setInitialStandardDeviations()");
                sigma = 0.5;
            }

            if (sigma <= 0 || math.min(diagD) <= 0)
            {
                error("initial standard deviations not specified or non-positive, " +
                "use setInitialStandarddeviations()");
                sigma = 1;
            }
            /* save initial standard deviation */
            if (startsigma == null || startsigma.Length == 1)
            {
                startsigma = new double[N];
                for (i = 0; i < N; ++i)
                {
                    startsigma[i] = sigma * diagD[i];
                }
            }
            maxstartsigma = math.max(startsigma);
            minstartsigma = math.min(startsigma);
            axisratio = maxstartsigma / minstartsigma; // axis parallel distribution

            /* expand typicalX, might still be null afterwards */
            typicalX = expandToDimension(typicalX, N);

            /* Initialization of xmean */
            xmean = expandToDimension(xmean, N);
            if (xmean == null)
            {
                /* set via typicalX */
                if (typicalX != null)
                {
                    xmean = (double[])typicalX.Clone();
                    for (i = 0; i < N; ++i)
                        xmean[i] += sigma * diagD[i] * nextGaussian(rand);
                    /* set via boundaries, is depriciated */
                }
                else if (math.max(UBound) < Double.MaxValue
                      && math.min(LBound) > -Double.MaxValue)
                {
                    error("no initial search point (solution) X or typical X specified");
                    xmean = new double[N];
                    for (i = 0; i < N; ++i)
                    { /* TODO: reconsider this algorithm to set X0 */
                        double offset = sigma * diagD[i];
                        double range = (UBound[i] - LBound[i] - 2 * sigma * diagD[i]);
                        if (offset > 0.4 * (UBound[i] - LBound[i]))
                        {
                            offset = 0.4 * (UBound[i] - LBound[i]);
                            range = 0.2 * (UBound[i] - LBound[i]);
                        }
                        xmean[i] = LBound[i] + offset + nextDouble(rand) * range;
                    }
                }
                else
                {
                    error("no initial search point (solution) X or typical X specified");
                    xmean = new double[N];
                    for (i = 0; i < N; ++i)
                        xmean[i] = nextDouble(rand);
                }
            }

            Debug.Assert(xmean != null);
            Debug.Assert(sigma > 0);

            /* interpret missing option value */
            if (options.diagonalCovarianceMatrix < 0) // necessary for hello world message
                options.diagonalCovarianceMatrix = 1 * 150 * N / sp.lambda; // cave: duplication below

            /* non-settable parameters */
            pc = new double[N];
            ps = new double[N];
            B = new double[N][];
            for (int j = 0; j < B.Length; j++)
            {
                B[j] = new double[N];
            }

            C = new double[N][];
            for (int j = 0; j < C.Length; j++)
            {
                C[j] = new double[N];
            }

            xold = new double[N];
            BDz = new double[N];
            bestever_x = (double[])xmean.Clone();
            // bestever = new CMASolution(xmean);
            artmp = new double[N];


            fit.deltaFitHist = new double[5];
            fit.idxDeltaFitHist = -1;
            for (i = 0; i < fit.deltaFitHist.Length; ++i)
                fit.deltaFitHist[i] = 1.0;

            // code to be duplicated in reSizeLambda
            fit.fitness = new IntDouble[sp.getLambda()];   // including penalties, used yet
            fit.raw = new IntDouble[sp.getLambda()];       // raw function values
            fit.history = new double[10 + 30 * N / sp.getLambda()];

            arx = new double[sp.getLambda()][];
            for (int j = 0; j < arx.Length; j++)
            {
                arx[j] = new double[N];
            }

            population = new double[sp.getLambda()][];
            for (int j = 0; j < population.Length; j++)
            {
                population[j] = new double[N];
            }

            for (i = 0; i < sp.getLambda(); ++i)
            {
                fit.fitness[i] = new IntDouble();
                fit.raw[i] = new IntDouble();
            }

            // initialization
            for (i = 0; i < N; ++i)
            {
                pc[i] = 0;
                ps[i] = 0;
                for (int j = 0; j < N; ++j)
                {
                    B[i][j] = 0;
                }
                for (int j = 0; j < i; ++j)
                {
                    C[i][j] = 0;
                }
                B[i][i] = 1;
                C[i][i] = diagD[i] * diagD[i];
            }
            maxsqrtdiagC = Math.Sqrt(math.max(math.diag(C)));
            minsqrtdiagC = Math.Sqrt(math.min(math.diag(C)));
            countCupdatesSinceEigenupdate = 0;
            iniphase = false; // obsolete

            /* Some consistency check */
            for (i = 0; i < N; ++i)
            {
                if (LBound[i] > UBound[i])
                    error("lower bound is greater than upper bound");
                if (typicalX != null)
                {
                    if (LBound[i] > typicalX[i])
                        error("lower bound '" + LBound[i] + "'is greater than typicalX" + typicalX[i]);
                    if (UBound[i] < typicalX[i])
                        error("upper bound '" + UBound[i] + "' is smaller than typicalX " + typicalX[i]);
                }
            }
            String[] s = getMessages();
            
            //IRWAN DEBUG
            if (s.Length != 0)
            {
                if (!s[0].Equals(""))
                    warning("termination condition satisfied at initialization: \n  " + s[0]);
            }
            initialX = (double[])xmean.Clone(); // keep finally chosen initialX

            timings.start = CurrentTimeMillis();
            timings.starteigen = CurrentTimeMillis();

            state = 0;
            if (options.verbosity > -1)
                printlnHelloWorld();

            return new double[sp.getLambda()];

        }

        public CMAParameters getParameterDefaults()
        {
            return sp.getDefaults(N);
        }

        public CMAParameters getParameterDefaults(int N)
        {
            return sp.getDefaults(N);
        }

        //not implemented properties

        private void warning(String s)
        {
            println(" CMA-ES warning: " + s);
        }

        public void error(String s)
        {
            println(" CMA-ES error: " + s);
            throw new Exception(" CMA-ES error: " + s);
        }

        Timing timings = new Timing();

        void eigendecomposition(int flgforce)
        {
            /* Update B and D, calculate eigendecomposition */
            int i, j;

            if (countCupdatesSinceEigenupdate == 0 && flgforce < 2)
                return;

            //           20% is usually better in terms of running *time* (only on fast to evaluate functions)
            if (!flgdiag && flgforce <= 0 &&
                    (timings.eigendecomposition > 1000 + options.maxTimeFractionForEigendecomposition
                            * (CurrentTimeMillis() - timings.starteigen)
                            || countCupdatesSinceEigenupdate < 1.0 / sp.getCcov() / N / 5.0))
                return;

            if (flgdiag)
            {
                for (i = 0; i < N; ++i)
                {
                    diagD[i] = Math.Sqrt(C[i][i]);
                }
                countCupdatesSinceEigenupdate = 0;
                timings.starteigen = CurrentTimeMillis(); // reset starting time
                timings.eigendecomposition = 0;             // not really necessary
            }
            else
            {
                // set B <- C
                for (i = 0; i < N; ++i)
                    for (j = 0; j <= i; ++j)
                        B[i][j] = B[j][i] = C[i][j];

                // eigendecomposition
                double[] offdiag = new double[N];
                long firsttime = CurrentTimeMillis();
                tred2(N, B, diagD, offdiag);
                tql2(N, diagD, offdiag, B);
                timings.eigendecomposition += CurrentTimeMillis() - firsttime;

                if (options.checkEigenSystem > 0)
                    checkEigenSystem(N, C, diagD, B); // for debugging 

                // assign diagD to eigenvalue square roots
                for (i = 0; i < N; ++i)
                {
                    if (diagD[i] < 0) // numerical problem?
                        error("an eigenvalue has become negative");
                    diagD[i] = Math.Sqrt(diagD[i]);
                }
                countCupdatesSinceEigenupdate = 0;
            } // end Update B and D
            if (math.min(diagD) == 0) // error management is done elsewhere
                axisratio = Double.PositiveInfinity;
            else
                axisratio = math.max(diagD) / math.min(diagD);
        }

        int checkEigenSystem(int N, double[][] C, double[] diag, double[][] Q)
        {
            /* compute Q diag Q^T and Q Q^T to check */
            int i, j, k, res = 0;
            double cc, dd;
            String s;

            for (i = 0; i < N; ++i)
                for (j = 0; j < N; ++j)
                {
                    for (cc = 0.0, dd = 0.0, k = 0; k < N; ++k)
                    {
                        cc += diag[k] * Q[i][k] * Q[j][k];
                        dd += Q[i][k] * Q[j][k];
                    }
                    /* check here, is the normalization the right one? */
                    if (Math.Abs(cc - C[i > j ? i : j][i > j ? j : i]) / Math.Sqrt(C[i][i] * C[j][j]) > 1e-10
                            && Math.Abs(cc - C[i > j ? i : j][i > j ? j : i]) > 1e-9)
                    { /* quite large */
                        s = " " + i + " " + j + " " + cc + " " + C[i > j ? i : j][i > j ? j : i] + " " + (cc - C[i > j ? i : j][i > j ? j : i]);
                        warning("cmaes_t:Eigen(): imprecise result detected " + s);
                        ++res;
                    }
                    if (Math.Abs(dd - (i == j ? 1 : 0)) > 1e-10)
                    {
                        s = i + " " + j + " " + dd;
                        warning("cmaes_t:Eigen(): imprecise result detected (Q not orthog.) " + s);
                        ++res;
                    }
                }
            return res;
        }

        private void tred2(int n, double[][] V, double[] d, double[] e)
        {

            //  This is derived from the Algol procedures tred2 by
            //  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
            //  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
            //  Fortran subroutine in EISPACK.

            for (int j = 0; j < n; j++)
            {
                d[j] = V[n - 1][j];
            }

            // Householder reduction to tridiagonal form.

            for (int i = n - 1; i > 0; i--)
            {

                // Scale to avoid under/overflow.

                double scale = 0.0;
                double h = 0.0;
                for (int k = 0; k < i; k++)
                {
                    scale = scale + Math.Abs(d[k]);
                }
                if (scale == 0.0)
                {
                    e[i] = d[i - 1];
                    for (int j = 0; j < i; j++)
                    {
                        d[j] = V[i - 1][j];
                        V[i][j] = 0.0;
                        V[j][i] = 0.0;
                    }
                }
                else
                {

                    // Generate Householder vector.

                    for (int k = 0; k < i; k++)
                    {
                        d[k] /= scale;
                        h += d[k] * d[k];
                    }
                    double f = d[i - 1];
                    double g = Math.Sqrt(h);
                    if (f > 0)
                    {
                        g = -g;
                    }
                    e[i] = scale * g;
                    h = h - f * g;
                    d[i - 1] = f - g;
                    for (int j = 0; j < i; j++)
                    {
                        e[j] = 0.0;
                    }

                    // Apply similarity transformation to remaining columns.

                    for (int j = 0; j < i; j++)
                    {
                        f = d[j];
                        V[j][i] = f;
                        g = e[j] + V[j][j] * f;
                        for (int k = j + 1; k <= i - 1; k++)
                        {
                            g += V[k][j] * d[k];
                            e[k] += V[k][j] * f;
                        }
                        e[j] = g;
                    }
                    f = 0.0;
                    for (int j = 0; j < i; j++)
                    {
                        e[j] /= h;
                        f += e[j] * d[j];
                    }
                    double hh = f / (h + h);
                    for (int j = 0; j < i; j++)
                    {
                        e[j] -= hh * d[j];
                    }
                    for (int j = 0; j < i; j++)
                    {
                        f = d[j];
                        g = e[j];
                        for (int k = j; k <= i - 1; k++)
                        {
                            V[k][j] -= (f * e[k] + g * d[k]);
                        }
                        d[j] = V[i - 1][j];
                        V[i][j] = 0.0;
                    }
                }
                d[i] = h;
            }

            // Accumulate transformations.

            for (int i = 0; i < n - 1; i++)
            {
                V[n - 1][i] = V[i][i];
                V[i][i] = 1.0;
                double h = d[i + 1];
                if (h != 0.0)
                {
                    for (int k = 0; k <= i; k++)
                    {
                        d[k] = V[k][i + 1] / h;
                    }
                    for (int j = 0; j <= i; j++)
                    {
                        double g = 0.0;
                        for (int k = 0; k <= i; k++)
                        {
                            g += V[k][i + 1] * V[k][j];
                        }
                        for (int k = 0; k <= i; k++)
                        {
                            V[k][j] -= g * d[k];
                        }
                    }
                }
                for (int k = 0; k <= i; k++)
                {
                    V[k][i + 1] = 0.0;
                }
            }
            for (int j = 0; j < n; j++)
            {
                d[j] = V[n - 1][j];
                V[n - 1][j] = 0.0;
            }
            V[n - 1][n - 1] = 1.0;
            e[0] = 0.0;
        }

        private void tql2(int n, double[] d, double[] e, double[][] V)
        {

            //  This is derived from the Algol procedures tql2, by
            //  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
            //  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
            //  Fortran subroutine in EISPACK.

            for (int i = 1; i < n; i++)
            {
                e[i - 1] = e[i];
            }
            e[n - 1] = 0.0;

            double f = 0.0;
            double tst1 = 0.0;
            double eps = Math.Pow(2.0, -52.0);
            for (int l = 0; l < n; l++)
            {

                // Find small subdiagonal element

                tst1 = Math.Max(tst1, Math.Abs(d[l]) + Math.Abs(e[l]));
                int m = l;
                while (m < n)
                {
                    if (Math.Abs(e[m]) <= eps * tst1)
                    {
                        break;
                    }
                    m++;
                }

                // If m == l, d[l] is an eigenvalue,
                // otherwise, iterate.

                if (m > l)
                {
                    int iter = 0;
                    do
                    {
                        iter = iter + 1;  // (Could check iteration count here.)

                        // Compute implicit shift

                        double g = d[l];
                        double p = (d[l + 1] - g) / (2.0 * e[l]);
                        double r = math.hypot(p, 1.0);
                        if (p < 0)
                        {
                            r = -r;
                        }
                        d[l] = e[l] / (p + r);
                        d[l + 1] = e[l] * (p + r);
                        double dl1 = d[l + 1];
                        double h = g - d[l];
                        for (int i = l + 2; i < n; i++)
                        {
                            d[i] -= h;
                        }
                        f = f + h;

                        // Implicit QL transformation.

                        p = d[m];
                        double c = 1.0;
                        double c2 = c;
                        double c3 = c;
                        double el1 = e[l + 1];
                        double s = 0.0;
                        double s2 = 0.0;
                        for (int i = m - 1; i >= l; i--)
                        {
                            c3 = c2;
                            c2 = c;
                            s2 = s;
                            g = c * e[i];
                            h = c * p;
                            r = math.hypot(p, e[i]);
                            e[i + 1] = s * r;
                            s = e[i] / r;
                            c = p / r;
                            p = c * d[i] - s * g;
                            d[i + 1] = h + s * (c * g + s * d[i]);

                            // Accumulate transformation.

                            for (int k = 0; k < n; k++)
                            {
                                h = V[k][i + 1];
                                V[k][i + 1] = s * V[k][i] + c * h;
                                V[k][i] = c * V[k][i] - s * h;
                            }
                        }
                        p = -s * s2 * c3 * el1 * e[l] / dl1;
                        e[l] = s * p;
                        d[l] = c * p;

                        // Check for convergence.

                    } while (Math.Abs(e[l]) > eps * tst1);
                }
                d[l] = d[l] + f;
                e[l] = 0.0;
            }

            // Sort eigenvalues and corresponding vectors.

            for (int i = 0; i < n - 1; i++)
            {
                int k = i;
                double p = d[i];
                for (int j = i + 1; j < n; j++)
                {
                    if (d[j] < p)
                    { // NH find smallest k>i
                        k = j;
                        p = d[j];
                    }
                }
                if (k != i)
                {
                    d[k] = d[i]; // swap k and i 
                    d[i] = p;
                    for (int j = 0; j < n; j++)
                    {
                        p = V[j][i];
                        V[j][i] = V[j][k];
                        V[j][k] = p;
                    }
                }
            }
        } // tql2

        double[][] genoPhenoTransformation(double[][] popx, double[][] popy)
        {
            if (popy == null || popy == popx || popy.Length != popx.Length)
                popy = new double[popx.Length][];

            for (int i = 0; i < popy.Length; ++i)
                popy[i] = genoPhenoTransformation(popx[i], popy[i]);

            return popy;
        }

        double[][] phenoGenoTransformation(double[][] popx, double[][] popy)
        {
            if (popy == null || popy == popx || popy.Length != popx.Length)
                popy = new double[popx.Length][];

            for (int i = 0; i < popy.Length; ++i)
                popy[i] = phenoGenoTransformation(popx[i], popy[i]);

            return popy;
        }

        double[] genoPhenoTransformation(double[] x, double[] y)
        {
            if (y == null || y == x || y.Length != x.Length)
            {
                y = (double[])x.Clone();
                return y; // for now return an identical copy
            }
            for (int i = 0; i < N; ++i)
                y[i] = x[i];
            return y;
        }

        double[] phenoGenoTransformation(double[] x, double[] y)
        {
            if (y == null || y == x || y.Length != x.Length)
            {
                y = (double[])x.Clone();
                return y; // for now return an identical copy
            }
            for (int i = 0; i < N; ++i)
                y[i] = x[i];
            return y;
        }

        public double[][] samplePopulation()
        {
            int i, j, iNk;
            double sum;

            if (state < 0)
                init();
            else if (state < 3 && state > 2)
                error("mixing of calls to updateSingle() and samplePopulation() is not possible");
            else
                eigendecomposition(0); // latest possibility to generate B and diagD

            if (state != 1)
                ++countiter;
            state = 1; // can be repeatedly called without problem
            idxRecentOffspring = sp.getLambda() - 1; // not really necessary at the moment


            // ensure maximal and minimal standard deviations
            if (options.lowerStandardDeviations != null && options.lowerStandardDeviations.Length > 0)
                for (i = 0; i < N; ++i)
                {
                    double d = options.lowerStandardDeviations[Math.Min(i, options.lowerStandardDeviations.Length - 1)];
                    if (d > sigma * minsqrtdiagC)
                        sigma = d / minsqrtdiagC;
                }
            if (options.upperStandardDeviations != null && options.upperStandardDeviations.Length > 0)
                for (i = 0; i < N; ++i)
                {
                    double d = options.upperStandardDeviations[Math.Min(i, options.upperStandardDeviations.Length - 1)];
                    if (d < sigma * maxsqrtdiagC)
                        sigma = d / maxsqrtdiagC;
                }

            testAndCorrectNumerics();

            /* sample the distribution */
            for (iNk = 0; iNk < sp.getLambda(); ++iNk)
            { /*
            * generate scaled
            * random vector (D * z)
            */

                // code duplication from resampleSingle because of possible future resampling before GenoPheno
                /* generate scaled random vector (D * z) */
                if (flgdiag)
                    for (i = 0; i < N; ++i)
                        arx[iNk][i] = xmean[i] + sigma * diagD[i] * nextGaussian(rand);
                else
                {
                    for (i = 0; i < N; ++i)
                        artmp[i] = diagD[i] * nextGaussian(rand);

                    /* add mutation (sigma * B * (D*z)) */
                    for (i = 0; i < N; ++i)
                    {
                        for (j = 0, sum = 0; j < N; ++j)
                            sum += B[i][j] * artmp[j];
                        arx[iNk][i] = xmean[i] + sigma * sum;
                    }
                }
                // redo this while isOutOfBounds(arx[iNk])
            }

            // I am desperately missing a const/readonly/visible qualifier. 
            return population = genoPhenoTransformation(arx, population);

        }

        public double[] resampleSingle(int index)
        {
            int i, j;
            double sum;
            if (state != 1)
                error("call samplePopulation before calling resampleSingle(int index)");

            /* sample the distribution */
            /* generate scaled random vector (D * z) */
            if (flgdiag)
                for (i = 0; i < N; ++i)
                    arx[index][i] = xmean[i] + sigma * diagD[i] * nextGaussian(rand);
            else
            {
                for (i = 0; i < N; ++i)
                    artmp[i] = diagD[i] * nextGaussian(rand);

                /* add mutation (sigma * B * (D*z)) */
                for (i = 0; i < N; ++i)
                {
                    for (j = 0, sum = 0; j < N; ++j)
                        sum += B[i][j] * artmp[j];
                    arx[index][i] = xmean[i] + sigma * sum;
                }
            }
            return population[index] = genoPhenoTransformation(arx[index], population[index]);
        }

        public double mahalanobisNorm(double[] x, double[] mean)
        {
            double yi, snorm = 0;
            int i, j;
            // snorm = (x-mean)' Cinverse (x-mean) = (x-mean)' (BD^2B')^-1 (x-mean)
            //       = (x-mean)' B'^-1 D^-2 B^-1 (x-mean) 
            //       = (x-mean)' B D^-1 D^-1 B' (x-mean)
            //       = (D^-1 B' (x-mean))' * (D^-1 B' (x-mean))
            /* calculate z := D^(-1) * B^(-1) * BDz into artmp, we could have stored z instead */
            for (i = 0; i < N; ++i)
            {
                for (j = 0, yi = 0.0; j < N; ++j)
                    yi += B[j][i] * (x[j] - mean[j]);
                // yi = i-th component of B' (x-mean)
                snorm += yi * yi / diagD[i] / diagD[i];
            }
            return Math.Sqrt(snorm) / sigma;
        }

        public void updateDistribution(double[][] population, double[] functionValues)
        {
            updateDistribution(population, functionValues, 0);
        }

        public void updateDistribution(double[][] population, double[] functionValues, int nInjected)
        {
            // TODO: Needs to be tested yet for nInjected > 0
            // pass first input argument
            arx = phenoGenoTransformation(population, null); // TODO should still be tested
            for (int i = 0; i < nInjected; ++i)
            {
                warning("TODO: checking of injected solution has not yet been tested");
                // if (mahalanobisNorm(arx[0], xmean) > Math.sqrt(N) + 2) // testing: seems fine
                //     System.out.println(mahalanobisNorm(arx[i], xmean)/Math.sqrt(N));
                double upperLength = Math.Sqrt(N) + 2.0 * N / (N + 2.0);  // should become an interfaced parameter? 
                double fac = upperLength / mahalanobisNorm(arx[i], xmean);
                if (fac < 1)
                    for (int j = 0; j < N; ++j)
                        arx[i][j] = xmean[j] + fac * (arx[i][j] - xmean[j]);
            }
            updateDistribution(functionValues);
        }

        public void updateDistribution(double[] functionValues)
        {
            if (state == 3)
            {
                error("updateDistribution() was already called");
            }
            if (functionValues.Length != sp.getLambda())
                error("argument double[] funcionValues.length=" + functionValues.Length
                        + "!=" + "lambda=" + sp.getLambda());

            /* pass input argument */
            for (int i = 0; i < sp.getLambda(); ++i)
            {
                fit.raw[i].val = functionValues[i];
                fit.raw[i].i = i;
            }

            counteval += sp.getLambda();
            recentFunctionValue = math.min(fit.raw).val;
            recentMaxFunctionValue = math.max(fit.raw).val;
            recentMinFunctionValue = math.min(fit.raw).val;
            updateDistribution();
        }

        private void updateDistribution()
        {

            int i, j, k, iNk, hsig;
            double sum;
            double psxps;

            if (state == 3)
            {
                error("updateDistribution() was already called");
            }

            /* sort function values */
            
            //irwan to maximize
            Array.Sort(fit.raw, fit.raw[0]);
            if(taskType == "Maximize")
            {
                Array.Reverse(fit.raw);
            }

            for (iNk = 0; iNk < sp.getLambda(); ++iNk)
            {
                fit.fitness[iNk].val = fit.raw[iNk].val; // superfluous at time
                fit.fitness[iNk].i = fit.raw[iNk].i;
            }

            /* update fitness history */
            for (i = fit.history.Length - 1; i > 0; --i)
                fit.history[i] = fit.history[i - 1];
            fit.history[0] = fit.raw[0].val;

            /* save/update bestever-value */
            updateBestEver(arx[fit.raw[0].i], fit.raw[0].val,
                    counteval - sp.getLambda() + fit.raw[0].i + 1);

            /* re-calculate diagonal flag */
            flgdiag = (options.diagonalCovarianceMatrix == 1 || options.diagonalCovarianceMatrix >= countiter);
            if (options.diagonalCovarianceMatrix == -1) // options might have been re-read
                flgdiag = (countiter <= 1 * 150 * N / sp.lambda);  // CAVE: duplication of "default"

            /* calculate xmean and BDz~N(0,C) */
            for (i = 0; i < N; ++i)
            {
                xold[i] = xmean[i];
                xmean[i] = 0.0;
                for (iNk = 0; iNk < sp.getMu(); ++iNk)
                    xmean[i] += sp.getWeights()[iNk] * arx[fit.fitness[iNk].i][i];
                BDz[i] = Math.Sqrt(sp.getMueff()) * (xmean[i] - xold[i]) / sigma;
            }

            /* cumulation for sigma (ps) using B*z */
            if (flgdiag)
            {
                /* given B=I we have B*z = z = D^-1 BDz  */
                for (i = 0; i < N; ++i)
                {
                    ps[i] = (1.0 - sp.getCs()) * ps[i]
                                                   + Math.Sqrt(sp.getCs() * (2.0 - sp.getCs()))
                                                   * BDz[i] / diagD[i];
                }
            }
            else
            {
                /* calculate z := D^(-1) * B^(-1) * BDz into artmp, we could have stored z instead */
                for (i = 0; i < N; ++i)
                {
                    for (j = 0, sum = 0.0; j < N; ++j)
                        sum += B[j][i] * BDz[j];
                    artmp[i] = sum / diagD[i];
                }
                /* cumulation for sigma (ps) using B*z */
                for (i = 0; i < N; ++i)
                {
                    for (j = 0, sum = 0.0; j < N; ++j)
                        sum += B[i][j] * artmp[j];
                    ps[i] = (1.0 - sp.getCs()) * ps[i]
                                                   + Math.Sqrt(sp.getCs() * (2.0 - sp.getCs())) * sum;
                }
            }

            /* calculate norm(ps)^2 */
            psxps = 0;
            for (i = 0; i < N; ++i)
                psxps += ps[i] * ps[i];

            /* cumulation for covariance matrix (pc) using B*D*z~N(0,C) */
            hsig = 0;
            if (Math.Sqrt(psxps)
                    / Math.Sqrt(1.0 - Math.Pow(1.0 - sp.getCs(), 2.0 * countiter))
                    / sp.chiN < 1.4 + 2.0 / (N + 1.0))
            {
                hsig = 1;
            }
            for (i = 0; i < N; ++i)
            {
                pc[i] = (1.0 - sp.getCc()) * pc[i] + hsig
                * Math.Sqrt(sp.getCc() * (2.0 - sp.getCc())) * BDz[i];
            }

            /* stop initial phase, not in use anymore as hsig does the job */
            if (iniphase
                    && countiter > Math.Min(1 / sp.getCs(), 1 + N / sp.getMucov()))
                if (psxps / sp.getDamps()
                        / (1.0 - Math.Pow((1.0 - sp.getCs()), countiter)) < N * 1.05)
                    iniphase = false;
            /* update of C */
            if (sp.getCcov() > 0 && iniphase == false)
            {

                ++countCupdatesSinceEigenupdate;

                /* update covariance matrix */
                for (i = 0; i < N; ++i)
                    for (j = (flgdiag ? i : 0);
                         j <= i; ++j)
                    {
                        C[i][j] = (1 - sp.getCcov(flgdiag))
                        * C[i][j]
                               + sp.getCcov()
                               * (1.0 / sp.getMucov())
                               * (pc[i] * pc[j] + (1 - hsig) * sp.getCc()
                                       * (2.0 - sp.getCc()) * C[i][j]);
                        for (k = 0; k < sp.getMu(); ++k)
                        { /*
                    * additional rank mu
                    * update
                    */
                            C[i][j] += sp.getCcov() * (1 - 1.0 / sp.getMucov())
                            * sp.getWeights()[k]
                                              * (arx[fit.fitness[k].i][i] - xold[i])
                                              * (arx[fit.fitness[k].i][j] - xold[j]) / sigma
                                              / sigma;
                        }
                    }
                maxsqrtdiagC = Math.Sqrt(math.max(math.diag(C)));
                minsqrtdiagC = Math.Sqrt(math.min(math.diag(C)));
            } // update of C

            /* update of sigma */
            sigma *= Math.Exp(((Math.Sqrt(psxps) / sp.chiN) - 1) * sp.getCs()
                    / sp.getDamps());

            state = 3;
        }

        double[] assignNew(double[] rhs, double[] lhs)
        {
            Debug.Assert(rhs != null);
            if (lhs != null && lhs != rhs && lhs.Length == rhs.Length)
                for (int i = 0; i < lhs.Length; ++i)
                    lhs[i] = rhs[i];
            else
                lhs = (double[])rhs.Clone();
            return lhs;
        }

        void updateBestEver(double[] x, double fitness, long eval)
        {
            //IRWAN for maximize+maximize
            if (taskType == "Maximize")
            {
                if (fitness > bestever_fit || Double.IsNaN(bestever_fit))
                {  // countiter == 1 not needed anymore
                    bestever_fit = fitness;
                    bestever_eval = eval;
                    bestever_x = assignNew(x, bestever_x); // save (hopefully) efficient assignment
                }
            }
            else if (taskType == "Minimize")
            {
                if (fitness < bestever_fit || Double.IsNaN(bestever_fit))
                {  // countiter == 1 not needed anymore
                    bestever_fit = fitness;
                    bestever_eval = eval;
                    bestever_x = assignNew(x, bestever_x); // save (hopefully) efficient assignment
                }
            }
        }

        public double getAxisRatio()
        {
            return axisratio;
        }

        public CMASolution getBestSolution()
        {
            return new CMASolution(bestever_x, bestever_fit, bestever_eval);
        }

        public CMASolution setFitnessOfMeanX(double fitness)
        {
            xmean_fit = fitness;
            ++counteval;
            updateBestEver(xmean, fitness, counteval);
            return new CMASolution(bestever_x, bestever_fit, bestever_eval);
        }

        public double[] getBestX()
        {
            if (state < 0)
                return null;
            return (double[])bestever_x.Clone();
        }

        public double getBestFunctionValue()
        {
            if (state < 0)
                return Double.NaN;
            return bestever_fit;
        }

        public long getBestEvaluationNumber()
        {
            return bestever_eval;
        }

        public ISolutionPoint getBestRecentSolution()
        {
            return new CMASolution(genoPhenoTransformation(arx[fit.raw[0].i], null),
                    fit.raw[0].val,
                    counteval - sp.getLambda() + fit.raw[0].i + 1);
        }

        public double[] getBestRecentX()
        {
            return genoPhenoTransformation(arx[fit.raw[0].i], null);
        }

        public double getBestRecentFunctionValue()
        {
            return recentMinFunctionValue;
        }

        public double getWorstRecentFunctionValue()
        {
            return recentMaxFunctionValue;
        }

        public double[] getMeanX()
        {
            return (double[])xmean.Clone();
        }

        public int getDimension()
        {
            return N;
        }

        public long getCountEval()
        {
            return counteval;
        }

        public long getCountIter()
        {
            return countiter;
        }

        public double[] getInitialX()
        {
            if (state < 0)
                error("initiaX not yet available, init() must be called first");
            return (double[])initialX.Clone();
        }

        public Random getRand()
        {
            return rand;
        }

        //not implemented properties

        public long getSeed()
        {
            return seed;
        }

        public long setCountEval(long c)
        {
            return counteval = c;
        }

        public void setDimension(int n)
        {
            if ((lockDimension > 0 || state >= 0) && N != n)
                error("dimension cannot be changed anymore or contradicts to initialX");
            N = n;
        }

        public void setTypicalX(double x)
        {
            if (state >= 0)
                error("typical x cannot be set anymore");
            typicalX = new double[] { x }; // allows "late binding" of dimension
        }

        public void setTypicalX(double[] x)
        {
            if (state >= 0)
                error("typical x cannot be set anymore");
            if (x.Length == 1)
            { // to make properties work
                setTypicalX(x[0]);
                return;
            }
            if (N < 1)
                setDimension(x.Length);
            if (N != x.Length)
                error("dimensions N=" + N + " and input x.length=" + x.Length + "do not agree");
            typicalX = new double[N];
            for (int i = 0; i < N; ++i)
                typicalX[i] = x[i];
            lockDimension = 1;
        }

        public void setInitialStandardDeviation(double startsigma)
        {
            if (state >= 0)
                error("standard deviations cannot be set anymore");
            this.startsigma = new double[] { startsigma };
        }

        public void setInitialStandardDeviations(double[] startsigma)
        {
            // assert startsigma != null; // assert should not be used for public arg check
            if (state >= 0)
                error("standard deviations cannot be set anymore");
            if (startsigma.Length == 1)
            { // to make properties work
                setInitialStandardDeviation(startsigma[0]);
                return;
            }
            if (N > 0 && N != startsigma.Length)
                error("dimensions N=" + N + " and input startsigma.length="
                        + startsigma.Length + "do not agree");
            if (N == 0)
                setDimension(startsigma.Length);

            Debug.Assert(N == startsigma.Length);

            this.startsigma = (double[])startsigma.Clone();
            lockDimension = 1;
        }

        public void setInitialX(double x)
        {
            if (state >= 0)
                error("initial x cannot be set anymore");
            xmean = new double[] { x }; // allows "late binding" of dimension N
        }

        public void setInitialX(double l, double u)
        {
            if (state >= 0)
                error("initial x cannot be set anymore");
            if (N < 1)
                error("dimension must have been specified before");
            xmean = new double[N];
            for (int i = 0; i < xmean.Length; ++i)
                xmean[i] = l + (u - l) * nextDouble(rand);
            lockDimension = 1;
        }

        public void setInitialX(double[] l, double[] u)
        {
            if (state >= 0)
                error("initial x cannot be set anymore");
            if (l.Length != u.Length)
                error("length of lower and upper values disagree");
            setDimension(l.Length);
            xmean = new double[N];
            for (int i = 0; i < xmean.Length; ++i)
                xmean[i] = l[i] + (u[i] - l[i]) * nextDouble(rand);
            lockDimension = 1;
        }

        public void setInitialX(double[] x)
        {
            if (state >= 0)
                error("initial x cannot be set anymore");
            if (x.Length == 1)
            { // to make properties work
                setInitialX(x[0]);
                return;
            }
            if (N > 0 && N != x.Length)
                error("dimensions do not match");
            if (N == 0)
                setDimension(x.Length);

            Debug.Assert(N == x.Length);

            xmean = new double[N];
            for (int i = 0; i < N; ++i)
                xmean[i] = x[i];
            lockDimension = 1; // because xmean is set up
        }

        public void setRand(Random rand)
        {
            this.rand = rand;
        }

        public void setSeed(long seed)
        {
            if (state >= 0)
                throw new Exception("setting seed has no effect at this point");
            else
            {
                if (seed <= 0)
                    seed = CurrentTimeMillis();
                this.seed = seed;
                rand = new Random((int)seed);
            }
        }

        public String getPrintLine()
        {
            String s;
            if (state < 0)
            {
                s = string.Format("\n{0}", countiter) +
                        string.Format("\n{0}", 0) +
                        string.Format("\n{0}", counteval);
            }
            else
            {
                s = string.Format("\ncountiter = {0}", countiter) +
                    string.Format("\nidxRecentOffspring = {0}", idxRecentOffspring) +
                    string.Format("\nrecentFunctionValue = {0}", recentFunctionValue) +
                    string.Format("\ngetBestFunctionValue() - recentFunctionValue = {0}", getBestFunctionValue() - recentFunctionValue) +
                    string.Format("\nrecentMaxFunctionValue - recentFunctionValue = {0}", recentMaxFunctionValue - recentFunctionValue) +
                    string.Format("\nmath.maxidx(math.diag(C)) = {0}", math.maxidx(math.diag(C))) +
                    string.Format("\nsigma * maxsqrtdiagC = {0}", sigma * maxsqrtdiagC) +
                    string.Format("\nmath.minidx(math.diag(C)) = {0}", math.minidx(math.diag(C))) +
                    string.Format("\nsigma * minsqrtdiagC = {0}", sigma * minsqrtdiagC) +
                    string.Format("\nsigma * math.min(diagD) = {0}", sigma * math.min(diagD)) +
                    string.Format("\nsigma = {0}", sigma) +
                    string.Format("\naxisratio = {0}", axisratio) +
                    string.Format("\ntimespan = {0}", (CurrentTimeMillis() - timings.start) / 1000) +
                    string.Format("\neigendecomposetime = {0}", (timings.eigendecomposition / 1000));
            }

            return s;
        }

        public String getPrintAnnotation()
        {
            String s = "Iteration,#Fevals: rb Function Value Delta( best ,worst) |idx: Max SD idx: Min SD  | minsigD  sigma Axisratio | time, in eig";
            return s;
        }

        public String helloWorld()
        {
            DateTime now = DateTime.Now;
            string nowstring = now.ToString();
            String s;

            s = "(" + sp.getMu() + "," + sp.getLambda()
            + ")-CMA-ES(mu_eff=" + Math.Round(10.0 * sp.getMueff()) / 10.0 + "), Ver=\""
            + versionNumber
            + "\", dimension=" + N
            + ", " + options.diagonalCovarianceMatrix + " diagonal iter."
            + ", randomSeed=" + seed
            + " (" + nowstring + ")";

            return s;
        }

        public void println(String s)
        {
            Console.WriteLine(s);
            if (options.writeDisplayToFile > 0)
            { }
            // IRWAN writeToFile(options.outputFileNamesPrefix + "disp" + ".dat", s, 1);
        }

        public void println()
        {
            println(getPrintLine());
        }

        public void printlnAnnotation()
        {
            println(getPrintAnnotation());
        }

        public void printlnHelloWorld()
        {
            println(helloWorld());
        }

        public String getDataRowFitness()
        {
            String s;
            s = countiter + " " + counteval + " " + sigma + " " + axisratio + " "
            + bestever_fit + " ";
            if (mode == SINGLE_MODE)
                s += recentFunctionValue + " ";
            else
            {
                s += fit.raw[0].val + " ";
                s += fit.raw[sp.getLambda() / 2].val + " ";
                s += fit.raw[sp.getLambda() - 1].val + " ";
                s += math.min(diagD) + " "
                    + (math.maxidx(math.diag(C)) + 1) + " " + sigma * maxsqrtdiagC + " "
                    + (math.minidx(math.diag(C)) + 1) + " " + sigma * minsqrtdiagC;
                //for (int i = 0; i < sp.getLambda(); ++i) {
                //    s += fit.funValues[i].d + " ";
                //}
            }
            return s;
        }

        public String getDataRowXRecentBest()
        {
            int idx = 0;
            if (mode == SINGLE_MODE)
                idx = idxRecentOffspring;
            String s;
            s = countiter + " " + counteval + " " + sigma + " 0 "
                + (state == 1 ? Double.NaN : fit.raw[idx].val) + " ";
            for (int i = 0; i < N; ++i)
            {
                s += arx[fit.raw[idx].i][i] + " ";
            }
            return s;
        }

        public String getDataRowXMean()
        {
            String s;
            s = countiter + " " + counteval + " " + sigma + " 0 0 ";
            for (int i = 0; i < N; ++i)
            {
                s += xmean[i] + " ";
            }
            return s;
        }

        public String getDataRowAxlen()
        {
            String s; ;
            s = countiter + " " + counteval + " " + sigma + " " + axisratio + " "
               + maxsqrtdiagC / minsqrtdiagC + " ";
            double[] tmp = (double[])diagD.Clone();
            Array.Sort(tmp);
            for (int i = 0; i < N; ++i)
            {
                s += tmp[i] + " ";
            }
            return s;
        }

        public String getDataRowStddev()
        {
            String s;
            s = countiter + " " + counteval + " " + sigma + " "
            + (1 + math.maxidx(math.diag(C))) + " " + (1 + math.minidx(math.diag(C))) + " ";
            for (int i = 0; i < N; ++i)
            {
                s += sigma * Math.Sqrt(C[i][i]) + " ";
            }
            return s;
        }

        public String getDataC()
        {
            int i, j;
            String s;
            s = "%# " + countiter + " " + counteval + " " + sigma + "\n";
            for (i = 0; i < N; ++i)
            {
                for (j = 0; j < i; ++j) // ouput correlation in the lower half
                    s += C[i][j] / Math.Sqrt(C[i][i] * C[j][j]) + " ";
                for (j = i; j < N; ++j)
                    s += sigma * sigma * C[i][j] + " ";
                s += "\n";
            }
            return s;
        }

        private String[] fileswritten = new String[] { "" };

        //writetofile

        //writetodefaultfile

        private static readonly DateTime Jan1st2016 = new DateTime(2016, 1, 1, 0, 0, 0, DateTimeKind.Utc);

        public static long CurrentTimeMillis()
        {
            return (long)(DateTime.UtcNow - Jan1st2016).TotalMilliseconds;
        }

        private static readonly object syncLock = new object();

        public static double nextDouble(Random rand)
        {
            lock (syncLock)
            {
                return rand.NextDouble();
            }
        }

        public static double nextGaussian(Random rand)
        {
            double randStdNormal = Math.Sqrt(-2.0 * Math.Log(nextDouble(rand))) *
                         Math.Sin(2.0 * Math.PI * nextDouble(rand)); //random normal(0,1)
            return randStdNormal;
        }

        public string taskType;
    }

}
