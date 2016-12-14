using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace CMAESharp
{
    public class FitnessCollector
    {
        public double[] history;
        public IntDouble[] fitness; //not sure
        public IntDouble[] raw;
        public double[] deltaFitHist = new double[5];
        public int idxDeltaFitHist = 0;
    }
}
