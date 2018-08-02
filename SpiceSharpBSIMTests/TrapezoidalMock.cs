using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using SpiceSharp.Algebra;

namespace SpiceSharp.IntegrationMethods
{
    /// <summary>
    /// Trapezoidal method for testing
    /// </summary>
    public class TrapezoidalMock : Trapezoidal
    {
        /// <summary>
        /// Expose time
        /// </summary>
        /// <param name="time"></param>
        public void SetTime(double time) => Time = time;

        /// <summary>
        /// Set the old timesteps
        /// </summary>
        /// <param name="delta"></param>
        public void SetTimesteps(IEnumerable<double> delta)
        {
            using (var it = delta.GetEnumerator())
            {
                DeltaOld.Clear(index =>
                {
                    it.MoveNext();
                    return it.Current;
                });
            }
        }

        /// <summary>
        /// Initialize all solutions
        /// </summary>
        /// <param name="size"></param>
        public void InitializeSolutions(int size)
        {
            Solutions.Clear(index => new DenseVector<double>(size));
        }

        /// <summary>
        /// Expose solutions
        /// </summary>
        public void SetSolution(int index, IEnumerable<double> solution)
        {
            var i = 0;
            foreach (var value in solution)
            {
                Solutions[index][i] = value;
                i++;
            }
        }
    }
}
