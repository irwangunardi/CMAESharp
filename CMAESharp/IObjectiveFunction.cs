using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace CMAESharp
{
    public interface IObjectiveFunction
    {
        /** @param x  a point (candidate solution) in the pre-image of the objective function 
        @return  objective function value of the input search point  
        */
        double valueOf(double[] x);
        bool isFeasible(double[] x);
    }
}
