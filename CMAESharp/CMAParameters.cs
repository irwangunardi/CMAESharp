using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Diagnostics;

namespace CMAESharp
{
    [Serializable]
    public class CMAParameters
    {
        private static readonly long serialVersionUID = -1305062342816588003L;
        public int supplemented;
        public int locked;
        public int lambda;
        int mu;
        public double mucov;
        public double mueff;
        public double[] weights;
        public double damps;
        public double cs;
        public double cc;
        public double ccov;
        public double ccovsep;

        public double chiN;

        public CMAParameters()
        {
            mucov = -1;
            ccov = -1;
        }

        public string check()
        {
            if (lambda <= 1)
                return "offspring population size lambda must be greater than onem is " + lambda;
            if (mu < 1)
                return "parent number mu must be greater or equal to one, is " + mu;
            if (mu > lambda)
                return "parent number mu " + mu + " must be smaller or equal to offspring population size lambda " + lambda;
            if (weights.Length != mu)
                return "number of recombination weights " + weights.Length + " disagrees with parent number mu " + mu;
            if (cs <= 0 || cs > 1)
                return "0 < cs <= 1 must hold for step-size cumulation parameter cs, is " + cs;
            if (damps <= 0)
                return "step-size damping parameter damps must be greater than zero, is " + damps;
            if (cc <= 0 || cc > 1)
                return "0 < cc <= 1 must hold for cumulation parameter cc, is " + cc;
            if (mucov < 0)
                return "mucov >= 0 must hold, is " + mucov;
            if (ccov < 0)
                return "learning parameter ccov >= 0 must hold, is " + ccov;
            return "";
        }

        public CMAParameters getDefaults(int N)
        {
            if (N == 0)
                error("default parameters needs dimension been set");

            CMAParameters p = new CMAParameters();
            p.supplementRemainders(N, new CMAOptions());
            return p;
        }

        public CMAParameters getDefaults(int N, int lambda)
        {
            CMAParameters p = new CMAParameters();
            p.setLambda(lambda);
            p.supplementRemainders(N, new CMAOptions());
            return p;
        }

        public void supplementRemainders(int N, CMAOptions opts)
        {
            if (supplemented > 0)
                error("defaults cannot be supplemented twice");
            if (N == 0)
                error("dimension must be greater than zero");

            supplemented = 1;
            locked = 1;

            chiN = Math.Sqrt(N)
            * (1.0 - 1.0 / (4.0 * N) + 1.0 / (21.0 * N * N));

            if (lambda <= 0)
                lambda = (int)(4.0 + 3.0 * Math.Log(N));
            if (mu <= 0)
                mu = (int)Math.Floor(lambda / 2.0);

            if (weights == null)
                setWeights(mu, recombinationType);
            else if (weights.Length == 0)
                setWeights(mu, recombinationType);

            if (cs <= 0)
                cs = (mueff + 2) / (N + mueff + 3);

            if (damps <= 0)
                damps =
                    (1 + 2 * Math.Max(0, Math.Sqrt((mueff - 1.0) / (N + 1.0)) - 1))
                    * Math.Max(0.3, 1 -
                            N / (1e-6 + Math.Min(opts.stopMaxIter,
                                    opts.stopMaxFunEvals / lambda)))
                                    + cs;

            if (cc <= 0)
                cc = 4.0 / (N + 4.0);

            if (mucov < 0)
                mucov = mueff;

            if (ccov < 0)
            {
                ccov = 2.0 / (N + 1.41) / (N + 1.41) / mucov
                + (1 - (1.0 / mucov))
                * Math.Min(1, (2 * mueff - 1) / (mueff + (N + 2) * (N + 2)));
                ccovsep = Math.Min(1, ccov * (N + 1.5) / 3.0);
            }

            String s = check();
            if (s == null)
                ;
            else if (s.Equals(""))
                ;
            else
                error(s);
        }

        public int getMu()
        {
            return mu;
        }

        public void setMu(int mu)
        {
            if (locked != 0)
                error("parameters are locked");
            this.mu = mu;
        }

        public int getLambda()
        {
            return lambda;
        }

        int flgLambdaChanged = 0;

        void setLambda(int lambda)
        {
            if (locked != 0)
                error("parameters cannot be set anymore");
            this.lambda = lambda;
        }

        public int getPopulationSize()
        {
            return getLambda();
        }

        public void setPopulationSize(int lambda)
        {
            setLambda(lambda);
        }

        public enum RecombinationType { superlinear, linear, equal };
        RecombinationType recombinationType = RecombinationType.superlinear;

        public double[] getWeights()
        {
            return this.weights;
        }

        public void setRecombinationWeights(RecombinationType recombinationType)
        {
            if (locked != 0)
                error("parameters cannot be set anymore");
            this.recombinationType = recombinationType;
        }

        public void setRecombination(int mu, RecombinationType recombinationType)
        {
            if (locked != 0)
                error("parameters are locked");
            this.mu = mu;
            this.recombinationType = recombinationType;
        }

        private void setWeights(int mu, RecombinationType recombinationType)
        {
            double[] w = new double[mu];
            if (recombinationType == RecombinationType.equal)
                for (int i = 0; i < mu; ++i)
                    w[i] = 1;
            else if (recombinationType == RecombinationType.linear)
                for (int i = 0; i < mu; ++i)
                    w[i] = mu - i;
            else
                for (int i = 0; i < mu; ++i)
                    w[i] = (Math.Log(mu + 1) - Math.Log(i + 1));

            setWeights(w);
        }

        protected void setWeights(double[] weights)
        {
            //DEBUG IRWAN Debug.Assert(locked == 0);
            double sum = 0;
            for (int i = 0; i < weights.Length; ++i)
                sum += weights[i];
            for (int i = 0; i < weights.Length; ++i)
                weights[i] /= sum;
            this.weights = weights;
            double sum1 = 0;
            double sum2 = 0;
            for (int i = 0; i < mu; ++i)
            {
                sum1 += weights[i];
                sum2 += weights[i] * weights[i];
            }
            this.mueff = sum1 * sum1 / sum2;
        }

        public double getMueff()
        {
            return mueff;
        }

        public double getMucov()
        {
            return mucov;
        }

        public void setMucov(double mucov)
        {
            if (locked != 0)
                error("parameters cannot be set anymore");
            this.mucov = mucov;
        }

        public double getCcov(bool flgdiag)
        {
            if (flgdiag)
                return ccovsep;
            return ccov;
        }

        public double getCcov()
        {
            return ccov;
        }

        public void setCcov(double ccov)
        {
            this.ccov = ccov; // can be set anytime, cave: switching from diagonal to full cov
        }

        public double getDamps()
        {
            return damps;
        }

        public void setDamps(double damps)
        {
            if (locked != 0) // not really necessary!?
                error("parameters cannot be set anymore");
            this.damps = damps;
        }

        public double getCc()
        {
            return cc;
        }

        public void setCc(double cc)
        {
            this.cc = cc;
        }

        public double getCs()
        {
            return cs;
        }

        public void setCs(double cs)
        {
            if (locked != 0)
                error("parameters cannot be set anymore");
            this.cs = cs;
        }

        private void error(String s)
        {
            Console.WriteLine(" CMA-ES error: " + s);
            throw new Exception(" CMA-ES error: " + s);
        }
    }
}
