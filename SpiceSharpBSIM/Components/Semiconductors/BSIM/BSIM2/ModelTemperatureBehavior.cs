using System;
using System.Collections.Generic;
using SpiceSharp.Behaviors;
using SpiceSharp.Simulations;

namespace SpiceSharp.Components.BSIM2Behaviors
{
	/// <summary>
	/// Temperature behavior for a <see cref="BSIM2Model" />
	/// </summary>
	public class ModelTemperatureBehavior : Behavior, ITemperatureBehavior
	{
		/// <summary>
		/// Necessary behaviors and parameters
		/// </summary>
		protected ModelBaseParameters ModelParameters { get; private set; }
		
		/// <summary>
		/// Properties
		/// </summary>
		public double Vtm { get; private set; }

        /// <summary>
        /// Size-dependent parameters for this model
        /// </summary>
        public Dictionary<Tuple<double, double>, BSIM2SizeDependParams> Params { get; } = new Dictionary<Tuple<double, double>, BSIM2SizeDependParams>();

	    /// <summary>
	    /// Constructor
	    /// </summary>
	    public ModelTemperatureBehavior(string name) : base(name)
	    {
	    }
		
		/// <summary>
		/// Setup the behavior
		/// </summary>
		public override void Bind(Simulation simulation, BindingContext context)
		{
            base.Bind(simulation, context);

            // Get parameter sets
			ModelParameters = context.GetParameterSet<ModelBaseParameters>();
		}
		
		/// <summary>
		/// Temperature behavior
		/// </summary>
		void ITemperatureBehavior.Temperature()
		{
			Vtm = 8.625e-5 * (ModelParameters.Temp + 273.0);
            Params.Clear();
		}
	}
}