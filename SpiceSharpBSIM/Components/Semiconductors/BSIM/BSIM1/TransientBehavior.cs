using System;
using SpiceSharp.Behaviors;
using SpiceSharp.IntegrationMethods;
using SpiceSharp.Simulations;

namespace SpiceSharp.Components.BSIM1Behaviors
{

    /// <summary>
    /// Transient behavior for a <see cref="BSIM1"/>
    /// </summary>
    public class TransientBehavior : Behavior, ITimeBehavior
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
        public override void Bind(Simulation simulation, BindingContext context)
        {
            if (context == null)
                throw new ArgumentNullException(nameof(context));
            _load = context.GetBehavior<BiasingBehavior>();
            _load.TranBehavior = this;

            var method = ((TimeSimulation)simulation).Method;
            Qb = method.CreateDerivative();
            Qg = method.CreateDerivative();
            Qd = method.CreateDerivative();
        }

        /// <summary>
        /// Gets the state of the dc.
        /// </summary>
        void ITimeBehavior.InitializeStates()
        {
        }

        /// <summary>
        /// Transient behavior
        /// </summary>
        void ITimeBehavior.Load()
        {
        }
    }
}