using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace CMAESharp
{
    public class IntDouble : IComparer<IntDouble>
    {
        public int i;    // unique integer value, useful after sorting
        public double val; // double value
        public IntDouble(double d, int i)
        {
            this.val = d;
            this.i = i;
        }
        public IntDouble(double d)
        {
            this.val = d;
        }
        public IntDouble()
        {
        }

        public int Compare(IntDouble o1, IntDouble o2)
        {
            if (o1.val < o2.val)
                return -1;
            if (o1.val > o2.val)
                return 1;
            if (o1.i < o2.i)
                return -1;
            if (o1.i > o2.i)
                return 1;
            return 0;
        }

        public bool Equals(IntDouble o1, IntDouble o2)
        {
            if (o1.Compare(o1, o2) == 0) // && o1.hashCode() == o2.hashCode()
                return true;
            return false;
        }
    }
}
