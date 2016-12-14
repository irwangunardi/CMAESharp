using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace cmaes
{
    public class SimulationFilter : IObjectiveFunction
    {
        private double worstValue = 0;
        private bool worstValueFilled = false;

        private List<double[]> inputDatabase = new List<double[]>();
        private List<double> fitnessDatabase = new List<double>();

        public List<Parameter> Parameters;
        public Simulator Simulator;

        public double valueOf(double[] x)
        {
            double res = worstValue;
            int outOfRange = 0;

            for (int i = 0; i < x.Length; i++)
            {
                //penalty for out of range
                if (x[i] > Parameters[i].MaxCandidateValue())
                {
                    outOfRange++;
                    res = res - Math.Abs(x[i] - Parameters[i].MaxCandidateValue());
                }
                else if (x[i] < Parameters[i].MinCandidateValue())
                {
                    outOfRange++;
                    res = res - Math.Abs(Parameters[i].MinCandidateValue() - x[i]);
                }
            }

            if (outOfRange != 0)//if its out of range, return trash function
            {
                return res;
            }
            else
            {
                //round down to nearest candidate value
                for (int i = 0; i < x.Length; i++)
                {
                    if (Parameters[i].Type != ParameterType.Text)
                    {
                        for (int j = 0; j < Parameters[i].CandidateValue.Count; j++)
                        {

                            if (Parameters[i].CandidateValue[j] > x[i])
                            {
                                x[i] = Parameters[i].CandidateValue[j - 1];
                                break;
                            }

                            if (j == Parameters[i].CandidateValue.Count - 1) //if its bigger than any candidate value, round down to the biggest candidate value
                            {
                                x[i] = Parameters[i].CandidateValue.Last();
                            }
                        }
                    }
                    else //if it is a text, round down by flooring as they dont have candidate values
                    {
                        x[i] = Math.Floor(x[i]);
                    }
                }

                //check if it already simulated
                for (int i = 0; i < inputDatabase.Count; i++)
                {
                    int match = 0;
                    for (int j = 0; j < x.Length; j++)
                    {
                        if (inputDatabase[i][j] == x[j])
                        {
                            match++;
                        }
                    }
                    if (match == x.Length)
                    {
                        return fitnessDatabase[i];
                    }
                }

                //run simulator
                res = Simulator.Run(ToSimulator(x));

                //update database
                inputDatabase.Add((double[])x.Clone());
                fitnessDatabase.Add(res);

                //update worst value
                if (res < worstValue || worstValueFilled == false)
                {
                    worstValue = res;
                    worstValueFilled = true;
                }

                //write to console
                Console.WriteLine("Job ID : {0}, Fitness : {1}",fitnessDatabase.Count(),res);

                return res;
            }

        }

        public bool isFeasible(double[] x) { return true; }

        private List<string> ToSimulator(double[] x)
        {
            List<string> temp = new List<string>();
            for (int i = 0; i < x.Length; i++)
            {
                if (Parameters[i].Type == ParameterType.Text)
                {
                    temp.Add(Parameters[i].NumTextPairCandidateValue[(int)x[i]].TextValue + "@" + Parameters[i].NumTextPairCandidateValue[(int)x[i]].NumericalValue.ToString());
                }
                else if (Parameters[i].Type == ParameterType.Integer)
                {
                    temp.Add(((int)x[i]).ToString());
                }
                else if (Parameters[i].Type == ParameterType.Real)
                {
                    temp.Add(x[i].ToString());
                }
            }

            temp.Add((((int)double.Parse(temp[4])) - 2).ToString());
            temp.Add((((int)double.Parse(temp[4])) - 1).ToString());
            temp.Add((((int)double.Parse(temp[10])) - 2).ToString());
            temp.Add((((int)double.Parse(temp[10])) - 1).ToString());
            temp.Add((((int)double.Parse(temp[16])) - 2).ToString());
            temp.Add((((int)double.Parse(temp[16])) - 1).ToString());
            temp.Add((((int)double.Parse(temp[22])) - 2).ToString());
            temp.Add((((int)double.Parse(temp[22])) - 1).ToString());
            temp.Add((((int)double.Parse(temp[28])) - 2).ToString());
            temp.Add((((int)double.Parse(temp[28])) - 1).ToString());

            if (((int)double.Parse(temp[5])) >= 3) { temp.Add("OPEN"); } else { temp.Add("CLOSED"); }
            if (((int)double.Parse(temp[11])) >= 3) { temp.Add("OPEN"); } else { temp.Add("CLOSED"); }
            if (((int)double.Parse(temp[17])) >= 3) { temp.Add("OPEN"); } else { temp.Add("CLOSED"); }
            if (((int)double.Parse(temp[23])) >= 3) { temp.Add("OPEN"); } else { temp.Add("CLOSED"); }
            if (((int)double.Parse(temp[29])) >= 3) { temp.Add("OPEN"); } else { temp.Add("CLOSED"); }
            return temp;
        }
    }
}
