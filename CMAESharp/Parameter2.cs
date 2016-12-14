using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace cmaes
{
    public class Parameter
    {
        public Parameter(double min, double max, double precision)
        {
            double candidate = min;
            while (true)
            {
                _CandidateValue.Add(candidate);
                candidate += precision;
                if (candidate > max)
                {
                    break;
                }
            }

        }

        public Parameter()
        {

        }

        public ParameterType Type { get; set; }
        public ParameterSource Source { get; set; }

        //if type is double/integer
        private List<double> _CandidateValue = new List<double>();
        public List<double> CandidateValue
        {
            get
            {
                List<double> SortedCandidateValue = new List<double>(_CandidateValue);
                SortedCandidateValue.Sort();
                return SortedCandidateValue;
            }
            set
            {
                _CandidateValue = value;
            }
        }

        //if type is text
        private List<NumericalTextPair> _NumTextPairCandidateValue = new List<NumericalTextPair>();
        public List<NumericalTextPair> NumTextPairCandidateValue
        {
            get { return _NumTextPairCandidateValue; }
            set { _NumTextPairCandidateValue = value; }
        }

        public double MaxCandidateValue()//add some gap to overcome round down
        {
            if (this.Type != ParameterType.Text)
            {
                List<double> SortedCandidateValue = new List<double>(_CandidateValue);
                SortedCandidateValue.Sort();
                return SortedCandidateValue.Last() + (SortedCandidateValue.Last() - SortedCandidateValue[SortedCandidateValue.Count - 2]);
            }
            else
            {
                return _NumTextPairCandidateValue.Count();
            }
        }

        public double MinCandidateValue()
        {
            if (this.Type != ParameterType.Text)
            {
                List<double> SortedCandidateValue = new List<double>(_CandidateValue);
                SortedCandidateValue.Sort();
                return SortedCandidateValue.First();
            }
            else
            {
                return 0;
            }
        }

    }

    public struct NumericalTextPair
    {
        public double NumericalValue { get; set; }
        public string TextValue { get; set; }
    }

    public enum ParameterType
    {
        Real = 0,
        Integer = 1,
        Text = 2,
    }

    public enum ParameterSource
    {
        CandidateValue = 0,
        Formula = 1
    }

}
