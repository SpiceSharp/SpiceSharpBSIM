using System;
using SpiceSharp.Algebra;
using SpiceSharp.Behaviors;
using SpiceSharp.IntegrationMethods;
using SpiceSharp.Simulations;
namespace SpiceSharp.Components.BSIM1Behaviors
{

    /// <summary>
    /// Transient behavior for a <see cref="BSIM1"/>
    /// </summary>
    public class TransientBehavior : BaseTransientBehavior
    {

        /// <summary>
        /// Necessary behaviors and parameters
        /// </summary>
        private LoadBehavior _load;

        /// <summary>
        /// States
        /// </summary>
        public StateDerivative Qb { get; private set; }
        public StateDerivative Qg { get; private set; }
        public StateDerivative Qd { get; private set; }

        /// <summary>
        /// Constructor
        /// </summary>
        public TransientBehavior(Identifier name) : base(name)
        {
        }

        /// <summary>
        /// Setup the behavior
        /// </summary>
        public override void Setup(SetupDataProvider provider)
        {
            if (provider == null)
                throw new ArgumentNullException(nameof(provider));
            _load = provider.GetBehavior<LoadBehavior>("entity");
            _load.TranBehavior = this;
        }

        /// <summary>
        /// Create states
        /// </summary>
        public override void CreateStates(StatePool states)
        {
            Qb = states.CreateDerivative();
            Qg = states.CreateDerivative();
            Qd = states.CreateDerivative();
        }

        /// <summary>
        /// Get equation pointers
        /// </summary>
        public override void GetEquationPointers(Solver<double> solver)
        {
        }

        /// <summary>
        /// Transient behavior
        /// </summary>
        public override void Transient(TimeSimulation simulation)
        {
            // Do nothing
        }

        /// <summary>
        /// Truncate
        /// </summary>
        /// <returns></returns>
        public override double Truncate()
        {
            var timetmp = double.PositiveInfinity;
            timetmp = Math.Min(timetmp, Qb.LocalTruncationError());
            timetmp = Math.Min(timetmp, Qg.LocalTruncationError());
            timetmp = Math.Min(timetmp, Qd.LocalTruncationError());
            return timetmp;
        }
    }
}