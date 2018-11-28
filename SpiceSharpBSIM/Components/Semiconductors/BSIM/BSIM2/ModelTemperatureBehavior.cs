using System;
using System.Collections.Generic;
using SpiceSharp.Behaviors;
using SpiceSharp.Simulations;
using SpiceSharp.Simulations.Behaviors;

namespace SpiceSharp.Components.BSIM2Behaviors
{
	/// <summary>
	/// Temperature behavior for a <see cref="BSIM2Model" />
	/// </summary>
	public class ModelTemperatureBehavior : ExportingBehavior, ITemperatureBehavior
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
		public override void Setup(Simulation simulation, SetupDataProvider provider)
		{
			if (provider == null)
				throw new ArgumentNullException(nameof(provider));

            // Get parameter sets
			ModelParameters = provider.GetParameterSet<ModelBaseParameters>();
		}
		
		/// <summary>
		/// Temperature behavior
		/// </summary>
		public void Temperature(BaseSimulation simulation)
		{
			Vtm = 8.625e-5 * (ModelParameters.Temp + 273.0);
            Params.Clear();
		}
	}
}