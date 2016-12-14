using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace cmaes
{
    public class FunctionCollector : AbstractObjectiveFunction
    {

        public FunctionCollector(double function_number,
        int flgRotate,
        double axisratio)
        {

            actFun = (int)(function_number);
            rotate = flgRotate;
            scaling = axisratio == 0 ? 1.0 : axisratio;

            if (actFun > maxFuncNumber)
                actFun = 1; /* sphere */

            // assign all functions by number here
            funs[0] = new RandFun();
            funs[10] = new Sphere();

            // convex-quadratic
            funs[30] = new Cigar(axisratio == 0 ? 1e3 : scaling);
            funs[40] = new Tablet(axisratio == 0 ? 1e3 : scaling);
            funs[50] = new Elli(axisratio == 0 ? 1e3 : scaling);
            funs[60] = new CigTab(axisratio == 0 ? 1e4 : scaling);
            funs[70] = new TwoAxes(axisratio == 0 ? 1e3 : scaling);

            // uni-modal, well, essentially 
            funs[80] = new Rosen();
            funs[90] = new DiffPow();
            funs[91] = new ssDiffPow();

            // multi-modal
            funs[150] = new Rastrigin(scaling, 10);
            funs[160] = new Ackley(scaling);

            //      funs[999]  = new Experimental();
            //      funs[]  = new ();
            //      funs[]  = new ();

        }

        static readonly int maxFuncNumber = 999;
        IObjectiveFunction[] funs = new IObjectiveFunction[maxFuncNumber + 1];
        int actFun = 0;
        int rotate = 0;
        double scaling = 1;
        Basis B = new Basis();


        public override double valueOf(double[] x)
        {
            x = (double[])x.Clone(); // regard input as imutable, not really Java philosophy
            if (rotate > 0)     // rotate
                x = B.Rotate(x);
            if (scaling != 1)
            { // scale 
                for (int i = 0; i < x.Length; ++i)
                    x[i] = Math.Pow(10, i / (x.Length - 1.0)) * x[i];
            }
            return funs[actFun] == null ? funs[0].valueOf(x) : funs[actFun].valueOf(x);
        }
        public bool isFeasible(double[] x)
        { // unfortunate code duplication
            //int i;
            //for (i = 0; i < x.length; ++i)
            //	if (x[i] < 0.01)
            //		return false;
            //return true;
            return funs[actFun].isFeasible(x);
        }

    }



    class RandFun : AbstractObjectiveFunction
    {
        Random rand = new Random();
        public override double valueOf(double[] x)
        {
            double res = CMAEvolutionStrategy.nextDouble(rand);
            return res;
        }
    }

    class Sphere : AbstractObjectiveFunction
    {
        public override double valueOf(double[] x)
        {
            double res = 0;
            for (int i = 0; i < x.Length; ++i)
                res += x[i] * x[i];
            return res;
        }
        /*
        public bool isFeasible(double[] x)
        {
            return true;
        }*/
    }

    public class Cigar : AbstractObjectiveFunction
    {
        public Cigar() : this(1e3) { }

        public Cigar(double axisratio)
        {
            factor = axisratio * axisratio;
        }
        public double factor = 1e6;

        public override double valueOf(double[] x)
        {
            double res = x[0] * x[0];
            for (int i = 1; i < x.Length; ++i)
                res += factor * x[i] * x[i];
            return res;
        }
    }

    public class Tablet : AbstractObjectiveFunction
    {
        public Tablet() : this(1e3) { }

        public Tablet(double axisratio)
        {
            factor = axisratio * axisratio;
        }
        public double factor = 1e6;

        public override double valueOf(double[] x)
        {
            double res = factor * x[0] * x[0];
            for (int i = 1; i < x.Length; ++i)
                res += x[i] * x[i];
            return res;
        }
    }

    public class CigTab : AbstractObjectiveFunction
    {
        public CigTab() : this(1e4) { }
        public CigTab(double axisratio)
        {
            factor = axisratio;
        }
        public double factor = 1e6;

        public override double valueOf(double[] x)
        {
            int end = x.Length - 1;
            double res = x[0] * x[0] / factor + factor * x[end] * x[end];
            for (int i = 1; i < end; ++i)
                res += x[i] * x[i];
            return res;
        }
    }

    public class TwoAxes : AbstractObjectiveFunction
    {
        public double factor = 1e6;
        public TwoAxes()
        {
        }
        public TwoAxes(double axisratio)
        {
            factor = axisratio * axisratio;
        }

        public override double valueOf(double[] x)
        {
            double res = 0;
            for (int i = 0; i < x.Length; ++i)
                res += (i < x.Length / 2 ? factor : 1) * x[i] * x[i];
            return res;
        }
    }

    public class ElliRotated : AbstractObjectiveFunction
    {
        public ElliRotated() : this(1e3) { }

        public ElliRotated(double axisratio)
        {
            factor = axisratio * axisratio;
        }
        public Basis B = new Basis();
        public double factor = 1e6;

        public override double valueOf(double[] x)
        {
            x = B.Rotate(x);
            double res = 0;
            for (int i = 0; i < x.Length; ++i)
                res += Math.Pow(factor, i / (x.Length - 1.0)) * x[i] * x[i];
            return res;
        }
    }

    public class Elli : AbstractObjectiveFunction
    {
        public Elli() : this(1e3) { }
        public Elli(double axisratio)
        {
            factor = axisratio * axisratio;
        }
        public double factor = 1e6;

        public override double valueOf(double[] x)
        {
            double res = 0;
            for (int i = 0; i < x.Length; ++i)
                res += Math.Pow(factor, i / (x.Length - 1.0)) * x[i] * x[i];
            return res;
        }
    }

    public class DiffPow : AbstractObjectiveFunction
    {
        public override double valueOf(double[] x)
        {
            double res = 0;
            for (int i = 0; i < x.Length; ++i)
                res += Math.Pow(Math.Abs(x[i]), 2.0 + 10 * (double)i / (x.Length - 1.0));
            return res;
        }
    }

    public class ssDiffPow : AbstractObjectiveFunction
    {
        public override double valueOf(double[] x)
        {
            return Math.Pow(new DiffPow().valueOf(x), 0.25);
        }
    }

    public class Rosen : AbstractObjectiveFunction
    {
        public override double valueOf(double[] x)
        {
            double res = 0;
            for (int i = 0; i < x.Length - 1; ++i)
                res += 1e2 * (x[i] * x[i] - x[i + 1]) * (x[i] * x[i] - x[i + 1]) +
                (x[i] - 1.0) * (x[i] - 1.0);
            return res;
        }
    }

    public class Ackley : AbstractObjectiveFunction
    {
        double axisratio = 1.0;
        public Ackley(double axra)
        {
            axisratio = axra;
        }
        public Ackley()
        {
        }

        public override double valueOf(double[] x)
        {
            double res = 0;
            double res2 = 0;
            double fac = 0;
            for (int i = 0; i < x.Length; ++i)
            {
                fac = Math.Pow(axisratio, (i - 1.0) / (x.Length - 1.0));
                res += fac * fac * x[i] * x[i];
                res2 += Math.Cos(2.0 * Math.PI * fac * x[i]);
            }
            return (20.0 - 20.0 * Math.Exp(-0.2 * Math.Sqrt(res / x.Length))
                    + Math.Exp(1.0) - Math.Exp(res2 / x.Length));
        }
    }

    public class Rastrigin : AbstractObjectiveFunction
    {
        public Rastrigin()
            : this(1, 10)
        {
        }
        public Rastrigin(double axisratio, double amplitude)
        {
            this.axisratio = axisratio;
            this.amplitude = amplitude;
        }
        public double axisratio = 1;
        public double amplitude = 10;

        public override double valueOf(double[] x)
        {
            double fac;
            double res = 0;
            for (int i = 0; i < x.Length; ++i)
            {
                fac = Math.Pow(axisratio, (i - 1.0) / (x.Length - 1.0));
                if (i == 0 && x[i] < 0)
                    fac *= 1.0;
                res += fac * fac * x[i] * x[i]
                   + amplitude * (1.0 - Math.Cos(2.0 * Math.PI * fac * x[i]));
            }
            return res;
        }
    }

    public class Basis
    {
        double[][] B; // usually field names should be lower case
        Random rand = new Random(2); // use not always the same basis

        public double[] Rotate(double[] x)
        {
            GenBasis(x.Length);
            double[] y = new double[x.Length];
            for (int i = 0; i < x.Length; ++i)
            {
                y[i] = 0;
                for (int j = 0; j < x.Length; ++j)
                    y[i] += B[i][j] * x[j];
            }
            return y;
        }
        double[][] Rotate(double[][] pop)
        {
            double[][] y = new double[pop.Length][];
            for (int i = 0; i < pop.Length; ++i)
            {
                y[i] = Rotate(pop[i]);
            }
            return y;
        }

        void GenBasis(int DIM)
        {
            if (B != null ? B.Length == DIM : false)
                return;

            double sp;
            int i, j, k;

            /* generate orthogonal basis */
            B = new double[DIM][];
            for (int l = 0; l < B.Length; l++)
            {
                B[l] = new double[DIM];
            }


            for (i = 0; i < DIM; ++i)
            {
                /* sample components gaussian */
                for (j = 0; j < DIM; ++j)
                    B[i][j] = CMAEvolutionStrategy.nextGaussian(rand);
                /* substract projection of previous vectors */
                for (j = i - 1; j >= 0; --j)
                {
                    for (sp = 0.0, k = 0; k < DIM; ++k)
                        sp += B[i][k] * B[j][k]; /* scalar product */
                    for (k = 0; k < DIM; ++k)
                        B[i][k] -= sp * B[j][k]; /* substract */
                }
                /* normalize */
                for (sp = 0.0, k = 0; k < DIM; ++k)
                    sp += B[i][k] * B[i][k]; /* squared norm */
                for (k = 0; k < DIM; ++k)
                    B[i][k] /= Math.Sqrt(sp);
            }
        }
    }
}
