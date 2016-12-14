using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace CMAESharp
{
    public class Timing
    {
        public Timing()
        {
            birth = CMAEvolutionStrategy.CurrentTimeMillis();
            start = birth; // on the save side 
        }
        public long birth; // time at construction, not really in use
        public long start; // time at end of init()
        public long starteigen; // time after flgdiag was turned off, ie when calls to eigen() start
        public long eigendecomposition = 0; // spent time in eigendecomposition
        public long writedefaultfiles = 0;        // spent time in writeToDefaultFiles

    }
}
