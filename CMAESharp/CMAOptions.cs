using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Text.RegularExpressions;

namespace CMAESharp
{
    [Serializable]
    public class CMAOptions
    {
        private static readonly long serialVersionUID = 2255162105325585121L;
        public long diagonalCovarianceMatrix = 0;
        public double[] lowerStandardDeviations;
        public double[] upperStandardDeviations;
        //IRWAN for maximize
        public double stopFitness = double.NaN;
        public double stopTolFun = 1e-12;
        public double stopTolFunHist = 1e-13;
        public double stopTolX = 0.0;
        public double stopTolXFactor = 1e-11;
        public double stopTolUpXFactor = 1e3;
        public long stopMaxFunEvals = long.MaxValue;
        public long stopMaxIter = long.MaxValue;
        public bool stopnow = false;
        public int verbosity = 1;
        public string outputFileNamesPrefix = "outcmaes";
        public int writeDisplayToFile = 1;
        public double maxTimeFractionForEigendecomposition = 0.2;
        public double maxTimeFractionForWriteToDefaultFiles = 0.1;
        public int checkEigenSystem = 0;

        void setOptions()
        {
            //unimplemented
        }

        public double getFirstToken(string s, double def)
        {
            if (s == null)
            {
                return def;
            }
            string[] ar = Regex.Split(s, "\\s+");
            if (ar[0].Equals(""))
            {
                return def;
            }
            return double.Parse(ar[0]);
        }

        public string getFirstToken(string s)
        {
            if (s == null)
            {
                return "";
            }
            string[] ar = Regex.Split(s, "\\s+");
            return ar[0];
        }

        public int getFirstToken(string s, int def)
        {
            if (s == null)
            {
                return def;
            }
            string[] ar = Regex.Split(s, "\\s+");
            if (ar[0].Equals(""))
            {
                return def;
            }
            return int.Parse(ar[0]);
        }

        private string removeComments(string s)
        {
            int i;
            { i = s.IndexOf("#"); }
            if (i >= 0)
            { s = s.Substring(0, i); }
            i = s.IndexOf("!");
            if (i >= 0)
            { s = s.Substring(0, i); }
            i = s.IndexOf("%");
            if (i >= 0)
            { s = s.Substring(0, i); }
            i = s.IndexOf("//");
            if (i >= 0)
            { s = s.Substring(0, i); }
            return s;
        }

        private long getFirstToken(String s, long def)
        {
            if (s == null)
                return def;
            String[] ar = Regex.Split(removeComments(s),"\\s+");
            if (ar[0].Equals(""))
                return def;
            return long.Parse(ar[0]);
        }

        string[] getAllToken(String s)
        {
            return Regex.Split(removeComments(s), "\\s+");
        }

        double[] parseDouble(String[] ars)
        {
            double[] ard = new double[ars.Length];
            for (int i = 0; i < ars.Length; ++i)
            {
                ard[i] = Double.Parse(ars[i]);
            }
            return ard;
        }
    }
}
