using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace CMAESharp
{
    public interface ISolutionPoint
    {
        /** objective function value (fitness) of the search point x */
        double getFitness();
        /** count at what evaluation number the search point x was evaluated */
        long getEvaluationNumber();
        /** value of the point in search space, that is in the 
         * preimage of the objective function to be optimized */
        double[] getX();

        /** objective function value (fitness) of the search point x */
        void setFitness(double fitness); // TODO better FunctionValue than Fitness ? 
        /** count at what evaluation number the search point x was evaluated */
        void setEvaluationNumber(long evaluation);
        /** value of the solution point in search space, the 
         * preimage of the objective function to be optimized */
        void setX(double[] x);
    }
}
