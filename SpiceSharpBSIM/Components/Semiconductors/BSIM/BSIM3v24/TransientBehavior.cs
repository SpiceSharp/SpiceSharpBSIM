using System;
using SpiceSharp.Algebra;
using SpiceSharp.Behaviors;
using SpiceSharp.IntegrationMethods;
using SpiceSharp.Simulations;
namespace SpiceSharp.Components.BSIM3v24Behaviors
{
	
	/// <summary>
	/// Transient behavior for a <see cref="BSIM3"/>
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
		public StateDerivative Qcheq { get; private set; }
		public StateDerivative Qcdump { get; private set; }
		
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
			_load = provider.GetBehavior<LoadBehavior>("instance");
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
			Qcheq = states.CreateDerivative();
			Qcdump = states.CreateDerivative();
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