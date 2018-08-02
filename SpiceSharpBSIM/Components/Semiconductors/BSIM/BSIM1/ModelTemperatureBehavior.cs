using System;
using SpiceSharp.Behaviors;
using SpiceSharp.Simulations;
namespace SpiceSharp.Components.BSIM1Behaviors
{

    /// <summary>
    /// Temperature behavior for a <see cref="BSIM1Model" />
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

        /// <summary>
        /// Constructor
        /// </summary>
        public ModelTemperatureBehavior(Identifier name) : base(name)
        {

        }

        /// <summary>
        /// Setup the behavior
        /// </summary>
        public override void Setup(SetupDataProvider provider)
        {
            if (provider == null)
                throw new ArgumentNullException(nameof(provider));
            _mbp = provider.GetParameterSet<ModelBaseParameters>("entity");
        }

        /// <summary>
        /// Temperature behavior
        /// </summary>
        public override void Temperature(BaseSimulation simulation)
        {
            double cox;
            if (_mbp.BulkJctPotential < 0.1)
            {
                _mbp.BulkJctPotential.Value = 0.1;
            }
            if (_mbp.SidewallJctPotential < 0.1)
            {
                _mbp.SidewallJctPotential.Value = 0.1;
            }
            cox = 3.453e-13 / (_mbp.OxideThickness * 1.0e-4);
            Cox = cox;
        }
    }
}