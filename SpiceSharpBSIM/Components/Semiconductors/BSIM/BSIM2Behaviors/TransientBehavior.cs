using SpiceSharp.Behaviors;
using SpiceSharp.IntegrationMethods;
using SpiceSharp.Simulations;

namespace SpiceSharp.Components.BSIM2Behaviors
{	
	/// <summary>
	/// Transient behavior for a <see cref="BSIM2"/>
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
        /// <param name="simulation">The simulation.</param>
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