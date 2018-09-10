using System;
using System.Collections.Generic;
using SpiceSharp.Behaviors;
using SpiceSharp.Simulations;
namespace SpiceSharp.Components.BSIM2Behaviors
{
	
	/// <summary>
	/// Temperature behavior for a <see cref="BSIM2Model" />
	/// </summary>
	public class ModelTemperatureBehavior : BaseTemperatureBehavior
	{
		
		/// <summary>
		/// Necessary behaviors and parameters
		/// </summary>
		private ModelBaseParameters _mbp;
		
		/// <summary>
		/// Properties
		/// </summary>
		public double Cox { get; private set; }
		public double Vdd2 { get; private set; }
		public double Vgg2 { get; private set; }
		public double Vbb2 { get; private set; }
		public double Vtm { get; private set; }

        /// <summary>
        /// Size-dependent parameters for this model
        /// </summary>
        public Dictionary<Tuple<double, double>, BSIM2SizeDependParams> Params { get; } = new Dictionary<Tuple<double, double>, BSIM2SizeDependParams>();

	    /// <summary>
	    /// Constructor
	    /// </summary>
	    public ModelTemperatureBehavior(Identifier name) : base(name)
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
			_mbp = provider.GetParameterSet<ModelBaseParameters>("entity");
		}
		
		/// <summary>
		/// Temperature behavior
		/// </summary>
		public override void Temperature(BaseSimulation simulation)
		{
			if (_mbp.BulkJctPotential < 0.1)
			{
				_mbp.BulkJctPotential.Value = 0.1;
			}
			if (_mbp.SidewallJctPotential < 0.1)
			{
				_mbp.SidewallJctPotential.Value = 0.1;
			}
			Cox = 3.453e-13 / (_mbp.Tox * 1.0e-4);
			Vdd2 = 2.0 * _mbp.Vdd;
			Vgg2 = 2.0 * _mbp.Vgg;
			Vbb2 = 2.0 * _mbp.Vbb;
			Vtm = 8.625e-5 * (_mbp.Temp + 273.0);
            Params.Clear();
		}
	}
}