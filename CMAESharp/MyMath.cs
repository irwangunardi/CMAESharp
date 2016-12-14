using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace CMAESharp
{
    public class MyMath
    {
        public int itest;

        double square(double d)
        {
            return d * d;
        }
        double prod(double[] ar)
        {
            double res = 1.0;
            for (int i = 0; i < ar.Length; ++i)
                res *= ar[i];
            return res;
        }

        public double median(double[] ar)
        {
            // need a copy of ar
            double[] ar2 = new double[ar.Length];
            for (int i = 0; i < ar.Length; ++i)
                ar2[i] = ar[i];
            Array.Sort(ar2);
            if (ar2.Length % 2 == 0)
                return (ar2[ar.Length / 2] + ar2[ar.Length / 2 - 1]) / 2.0;
            else
                return ar2[ar.Length / 2];
        }

        /** @return Maximum value of 1-D double array */
        public double max(double[] ar)
        {
            int i;
            double m;
            m = ar[0];
            for (i = 1; i < ar.Length; ++i)
            {
                if (m < ar[i])
                    m = ar[i];
            }
            return m;
        }

        /** sqrt(a^2 + b^2) without under/overflow. **/
        public double hypot(double a, double b)
        {
            double r = 0;
            if (Math.Abs(a) > Math.Abs(b))
            {
                r = b / a;
                r = Math.Abs(a) * Math.Sqrt(1 + r * r);
            }
            else if (b != 0)
            {
                r = a / b;
                r = Math.Abs(b) * Math.Sqrt(1 + r * r);
            }
            return r;
        }
        /** @return index of minium value of 1-D double array */
        public int minidx(double[] ar)
        {
            return minidx(ar, ar.Length - 1);
        }

        /** @return index of minium value of 1-D double 
         *   array between index 0 and maxidx 
         * @param ar double[] 
         * @param maxidx last index to be considered */
        public int minidx(double[] ar, int maxidx)
        {
            int i, idx;
            idx = 0;
            for (i = 1; i < maxidx; ++i)
            {
                if (ar[idx] > ar[i])
                    idx = i;
            }
            return idx;
        }

        /** @return index of minium value of 1-D double 
         *   array between index 0 and maxidx 
         * @param ar double[] 
         * @param maxidx last index to be considered */
        protected int minidx(IntDouble[] ar, int maxidx)
        {
            int i, idx;
            idx = 0;
            for (i = 1; i < maxidx; ++i)
            {
                if (ar[idx].val > ar[i].val)
                    idx = i;
            }
            return idx;
        }

        /** @return index of maximum value of 1-D double array */
        public int maxidx(double[] ar)
        {
            int i, idx;
            idx = 0;
            for (i = 1; i < ar.Length; ++i)
            {
                if (ar[idx] < ar[i])
                    idx = i;
            }
            return idx;
        }
        /** @return Minimum value of 1-D double array */
        public double min(double[] ar)
        {
            int i;
            double m;
            m = ar[0];
            for (i = 1; i < ar.Length; ++i)
            {
                if (m > ar[i])
                    m = ar[i];
            }
            return m;
        }

        /** @return Maximum value of 1-D Object array where the object implements Comparator 
         *    Example: max(Double arx, arx[0]) */
        public Double max(Double[] ar, IComparer<Double> c)
        {
            int i;
            Double m;
            m = ar[0];
            for (i = 1; i < ar.Length; ++i)
            {
                if (c.Compare(m, ar[i]) > 0)
                    m = ar[i];
            }
            return m;
        }

        /** @return Maximum value of 1-D IntDouble array */
        public IntDouble max(IntDouble[] ar)
        {
            int i;
            IntDouble m;
            m = ar[0];
            for (i = 1; i < ar.Length; ++i)
            {
                if (m.Compare(m, ar[i]) < 0)
                    m = ar[i];
            }
            return m;
        }

        /** @return Minimum value of 1-D IntDouble array */
        public IntDouble min(IntDouble[] ar)
        {
            int i;
            IntDouble m;
            m = ar[0];
            for (i = 1; i < ar.Length; ++i)
            {
                if (m.Compare(m, ar[i]) > 0)
                    m = ar[i];
            }
            return m;
        }

        /** @return Minimum value of 1-D Object array defining a Comparator */
        public Double min(Double[] ar, IComparer<Double> c)
        {
            int i;
            Double m;
            m = ar[0];
            for (i = 1; i < ar.Length; ++i)
            {
                if (c.Compare(m, ar[i]) < 0)
                    m = ar[i];
            }
            return m;
        }

        /**
         * @return Diagonal of an 2-D double array
         */
        public double[] diag(double[][] ar)
        {
            int i;
            double[] diag = new double[ar.Length];
            for (i = 0; i < ar.Length && i < ar[i].Length; ++i)
                diag[i] = ar[i][i];
            return diag;
        }

        /**
         * @return 1-D double array of absolute values of an 1-D double array
         */
        public double[] abs(double[] v)
        {
            double[] res = new double[v.Length];
            for (int i = 0; i < v.Length; ++i)
                res[i] = Math.Abs(v[i]);
            return res;
        }
    }
}
