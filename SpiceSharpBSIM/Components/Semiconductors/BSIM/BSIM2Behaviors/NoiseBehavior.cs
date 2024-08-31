using SpiceSharp.Attributes;
using SpiceSharp.Behaviors;
using SpiceSharp.Components;
using SpiceSharp.Components.NoiseSources;
using SpiceSharp.Simulations;
using System;

namespace SpiceSharpBSIM.Components.Semiconductors.BSIM.BSIM2Behaviors
{
    /// <summary>
    /// The noise behavior for a <see cref="BSIM2"/> transistor.
    /// </summary>
    [BehaviorFor(typeof(BSIM2)), AddBehaviorIfNo(typeof(INoiseBehavior))]
    [GeneratedParameters]
    public partial class NoiseBehavior : FrequencyBehavior, INoiseBehavior
    {
        private readonly ITemperatureSimulationState _temperature;
        private readonly INoiseSimulationState _state;
        private readonly NoiseThermal _rd, _rs, _id;
        private readonly NoiseGain _flicker;

        /// <inheritdoc/>
        [ParameterName("noise"), ParameterInfo("The total output noise density")]
        public double OutputNoiseDensity => _rd.OutputNoiseDensity + _rs.OutputNoiseDensity + _id.OutputNoiseDensity + _flicker.OutputNoiseDensity;

        /// <inheritdoc/>
        [ParameterName("onoise"), ParameterInfo("The total integrated output noise")]
        public double TotalOutputNoise => _rd.TotalOutputNoise + _rs.TotalOutputNoise + _id.TotalOutputNoise + _flicker.TotalOutputNoise;

        /// <inheritdoc/>
        [ParameterName("inoise"), ParameterInfo("The total integrated input noise")]
        public double TotalInputNoise => _rd.TotalInputNoise + _rs.TotalInputNoise + _id.TotalInputNoise + _flicker.TotalInputNoise;

        /// <summary>
        /// Gets the thermal noise of the drain resistor.
        /// </summary>
        [ParameterName("rd"), ParameterInfo("The thermal noise of the drain resistor")]
        public INoiseSource ThermalDrain => _rd;

        /// <summary>
        /// Gets the thermal noise of the source resistor.
        /// </summary>
        [ParameterName("rs"), ParameterInfo("The thermal noise of the source resistor")]
        public INoiseSource ThermalSource => _rs;

        /// <summary>
        /// Gets the shot noise of the drain current.
        /// </summary>
        [ParameterName("id"), ParameterInfo("The shot noise of the drain current")]
        public INoiseSource ShotDrainCurrent => _id;

        /// <summary>
        /// Gets the flicker noise.
        /// </summary>
        [ParameterName("flicker"), ParameterInfo("The flicker noise")]
        public INoiseSource Flicker => _flicker;

        /// <summary>
        /// Creates a new <see cref="NoiseBehavior"/>.
        /// </summary>
        /// <param name="context">The binding context.</param>
        public NoiseBehavior(ComponentBindingContext context)
            : base(context)
        {
            _temperature = context.GetState<ITemperatureSimulationState>();
            _state = context.GetState<INoiseSimulationState>();
            
            _rd = new NoiseThermal("rd", _drain, _drainPrime);
            _rs = new NoiseThermal("rs", _source, _sourcePrime);
            _id = new NoiseThermal("id", _drainPrime, _sourcePrime);
            _flicker = new NoiseGain("flicker", _drainPrime, _sourcePrime);
        }

        /// <inheritdoc />
        void INoiseSource.Initialize()
        {
            _rd.Initialize();
            _rs.Initialize();
            _id.Initialize();
            _flicker.Initialize();
        }

        void INoiseBehavior.Compute()
        {
            double m = Parameters.Multiplier;
            _rd.Compute(DrainConductance * m, _temperature.Temperature);
            _rs.Compute(SourceConductance * m, _temperature.Temperature);
            _id.Compute(2.0 / 3.0 * Math.Abs(Gm * m), _temperature.Temperature);
            _flicker.Compute(ModelParameters.FNcoef * m *
                 Math.Exp(ModelParameters.FNexp *
                 Math.Log(Math.Max(Math.Abs(Cd), 1e-38))) /
                 (_state.Point.Value.Frequency *
                 (Parameters.Width - ModelParameters.DeltaW * 1e-6) *
                 (Parameters.Length - ModelParameters.DeltaL * 1e-6) *
                 ModelTemperature.Cox * ModelTemperature.Cox));
        }
    }
}
