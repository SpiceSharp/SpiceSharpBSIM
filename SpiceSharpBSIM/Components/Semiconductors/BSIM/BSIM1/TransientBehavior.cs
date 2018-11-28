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
        private BiasingBehavior _load;

        /// <summary>
        /// States
        /// </summary>
        public StateDerivative Qb { get; private set; }
        public StateDerivative Qg { get; private set; }
        public StateDerivative Qd { get; private set; }

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
            _load = provider.GetBehavior<BiasingBehavior>();
            _load.TranBehavior = this;
        }

        /// <summary>
        /// Create states
        /// </summary>
        public override void CreateStates(IntegrationMethod method)
        {
            if (method == null)
                throw new ArgumentNullException(nameof(method));
            Qb = method.CreateDerivative();
            Qg = method.CreateDerivative();
            Qd = method.CreateDerivative();
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
    }
}