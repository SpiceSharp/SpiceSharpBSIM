using System;
using SpiceSharp.Algebra;
using SpiceSharp.Behaviors;
using SpiceSharp.IntegrationMethods;
using SpiceSharp.Simulations;

namespace SpiceSharp.Components.BSIM3Behaviors
{
    /// <summary>
    /// Transient behavior for a <see cref="BSIM3"/>
    /// </summary>
    public class TransientBehavior : Behavior, ITimeBehavior
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
        public override void Bind(Simulation simulation, BindingContext context)
        {
            base.Bind(simulation, context);

            this._load = context.GetBehavior<BiasingBehavior>("instance");
            _load.TranBehavior = this;

            var method = ((TimeSimulation)simulation).Method;
            Qb = method.CreateDerivative();
            Qg = method.CreateDerivative();
            Qd = method.CreateDerivative();
            Qcheq = method.CreateDerivative(false);
            Qcdump = method.CreateDerivative(false);
        }
        
        /// <summary>
        /// Initialize states.
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