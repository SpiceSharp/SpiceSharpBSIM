using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using SpiceSharp.IntegrationMethods;

namespace SpiceSharp.Simulations
{
    /// <summary>
    /// Transient mock
    /// </summary>
    public class TransientMock : Transient
    {
        /// <summary>
        /// Private variables
        /// </summary>
        public TrapezoidalMock MockMethod;

        private Dictionary<int, int> StateMap { get; } = new Dictionary<int, int>();

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="name">Name</param>
        public TransientMock(Identifier name) : base(name)
        {
            var config = ParameterSets.Get<TimeConfiguration>();
            MockMethod = new TrapezoidalMock();
            config.Method = MockMethod;
        }

        /// <summary>
        /// Setup the transient analysis
        /// </summary>
        public new void Setup(Circuit ckt)
        {
            base.Setup(ckt);
        }

        /// <summary>
        /// Unsetup the transient analysis
        /// </summary>
        /// <param name="ckt"></param>
        public void Unsetup(Circuit ckt)
        {
            base.Unsetup();
        }

        /// <summary>
        /// Calculate the solution
        /// </summary>
        public void Calculate()
        {
            if (!TimeIterate(TimeConfiguration.TranMaxIterations))
                throw new CircuitException("Failed iterations");
        }

        /// <summary>
        /// Set the solution
        /// </summary>
        /// <param name="index"></param>
        /// <param name="solution"></param>
        public void SetSolution(int index, IEnumerable<double> solution)
        {
            var rstate = States.Get<RealState>();

            // Solution 0 is the reference solution
            if (index == 0)
            {
                var i = 0;
                foreach (var v in solution)
                    rstate.Solution[i++] = v;
            }

            // Also update the solution
            MockMethod.SetSolution(index, solution);
        }

        /// <summary>
        /// Set the states in history
        /// </summary>
        /// <param name="index">Index</param>
        /// <param name="states">States</param>
        public void SetState(int index, IEnumerable<double> states)
        {
            if (states == null)
                return;
            var i = 0;

            // Set the current state
            if (index == 0)
            {
                i = 0;
                foreach (var value in states)
                {
                    if (StateMap.TryGetValue(i, out var mapped))
                        StatePool.History[index][mapped] = value;
                    else
                        StatePool.History[index][i] = value;
                    i++;
                }
            }

            i = 0;
            foreach (var value in states)
            {
                if (StateMap.TryGetValue(i, out var mapped))
                    StatePool.History[index + 1][mapped] = value;
                else
                    StatePool.History[index + 1][i] = value;
                i = 0;
            }
        }
    }
}
