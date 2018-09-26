using System;
using SpiceSharp.Algebra;
using SpiceSharp.Behaviors;
using SpiceSharp.IntegrationMethods;
using SpiceSharp.Simulations;
namespace SpiceSharp.Components.BSIM2Behaviors
{
	
	/// <summary>
	/// Transient behavior for a <see cref="BSIM2"/>
	/// </summary>
	public class TransientBehavior : BaseTransientBehavior
	{
		
		/// <summary>
		/// Necessary behaviors and parameters
		/// </summary>
		private LoadBehavior _load;
		
		/// <summary>
		/// Properties
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
		public override void Setup(Simulation simulation, SetupDataProvider provider)
		{
			if (provider == null)
				throw new ArgumentNullException(nameof(provider));
			_load = provider.GetBehavior<LoadBehavior>();
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
			
		}
	}
}