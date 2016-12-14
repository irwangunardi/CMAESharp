using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace CMAESharp
{
    public interface IObjectiveFunctionParallel
    {
        double[] valuesOf(double[][] pop);
    }
}
