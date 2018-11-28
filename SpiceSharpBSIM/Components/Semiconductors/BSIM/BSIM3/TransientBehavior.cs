using System;
using SpiceSharp.Algebra;
using SpiceSharp.Behaviors;
using SpiceSharp.IntegrationMethods;
using SpiceSharp.Simulations;
using SpiceSharp.Simulations.Behaviors;

namespace SpiceSharp.Components.BSIM3Behaviors
{
    /// <summary>
    /// Transient behavior for a <see cref="BSIM3"/>
    /// </summary>
    public class TransientBehavior : ExportingBehavior, ITimeBehavior
    {
        /// <summary>
        /// Necessary behaviors and parameters
        /// </summary>
        private BiasingBehavior _load;

        /// <summary>
        /// Properties
        /// </summary>
        public StateDerivative Qb { get; private set; }
        public StateDerivative Qg { get; private set; }
        public StateDerivative Qd { get; private set; }
        public StateDerivative Qcheq { get; private set; }
        public StateDerivative Qcdump { get; private set; }

        /// <summary>
        /// Constructor
        /// </summary>
        public TransientBehavior(string name) : base(name)
        {
        }

        /// <summary>
        /// Setup the behavior
        /// </summary>
        public override void Setup(Simulation simulation, SetupDataProvider provider)
        {
            if (provider == null)
                throw new ArgumentNullException(nameof(provider));
            this._load = provider.GetBehavior<BiasingBehavior>("instance");
            _load.TranBehavior = this;
        }

        /// <summary>
        /// Create states
        /// </summary>
        public void CreateStates(IntegrationMethod method)
        {
            Qb = method.CreateDerivative();
            Qg = method.CreateDerivative();
            Qd = method.CreateDerivative();
            Qcheq = method.CreateDerivative(false);
            Qcdump = method.CreateDerivative(false);
        }

        /// <summary>
        /// Gets the state of the dc.
        /// </summary>
        /// <param name="simulation">The simulation.</param>
        public void GetDcState(TimeSimulation simulation)
        {
        }

        /// <summary>
        /// Get equation pointers
        /// </summary>
        public void GetEquationPointers(Solver<double> solver)
        {
        }

        /// <summary>
        /// Transient behavior
        /// </summary>
        public void Transient(TimeSimulation simulation)
        {
        }
    }
}