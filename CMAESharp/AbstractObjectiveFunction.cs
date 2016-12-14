using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace CMAESharp
{
    public abstract class AbstractObjectiveFunction : IObjectiveFunction, IObjectiveFunctionParallel
    {
        abstract public double valueOf(double[] x);

        public double[] valuesOf(double[][] pop)
        {
            double[] res = new double[pop.Length];
            for (int i = 0; i < pop.Length; ++i)
                res[i] = valueOf(pop[i]);
            return res;
        }

        public bool isFeasible(double[] x)
        {
            return true;
        }
    }
}
