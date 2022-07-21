using SpiceSharp;
using SpiceSharp.Attributes;
using SpiceSharp.Behaviors;
using SpiceSharp.Components;
using SpiceSharp.Components.NoiseSources;
using SpiceSharp.Simulations;
using System;

namespace SpiceSharpBSIM.Components.Semiconductors.BSIM.BSIM3Behaviors
{
    /// <summary>
    /// Noise behavior for a <see cref="BSIM3"/> transistor.
    /// </summary>
    [BehaviorFor(typeof(BSIM3)), AddBehaviorIfNo(typeof(INoiseBehavior))]
    public class NoiseBehavior : FrequencyBehavior, INoiseBehavior
    {
        private readonly ITemperatureSimulationState _temperature;
        private readonly INoiseSimulationState _state;
        private readonly NoiseThermal _rd, _rs, _id;
        private readonly NoiseGain _flicker;

        /*
         Channel thermal and flicker noises are calculated based on the value
         of model->BSIM3noiMod.
         If model->BSIM3noiMod = 1,
            Channel thermal noise = SPICE2 model
            Flicker noise         = SPICE2 model
         If model->BSIM3noiMod = 2,
            Channel thermal noise = BSIM3 model
            Flicker noise         = BSIM3 model
         If model->BSIM3noiMod = 3,
            Channel thermal noise = SPICE2 model
            Flicker noise         = BSIM3 model
         If model->BSIM3noiMod = 4,
            Channel thermal noise = BSIM3 model
            Flicker noise         = SPICE2 model
         If model->BSIM3noiMod = 5,
            Channel thermal noise = SPICE2 model with linear/sat fix
            Flicker noise         = SPICE2 model
         If model->BSIM3noiMod = 6,
            Channel thermal noise = SPICE2 model with linear/sat fix
            Flicker noise         = BSIM3 model
         */

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
        /// <param name="context">The context.</param>
        public NoiseBehavior(ComponentBindingContext context)
            : base(context)
        {
            _temperature = context.GetState<ITemperatureSimulationState>();
            _state = context.GetState<INoiseSimulationState>();

            _rd = new NoiseThermal("rd", _drainPrime, _drain);
            _rs = new NoiseThermal("rs", _sourcePrime, _source);
            _id = new NoiseThermal("id", _drainPrime, _sourcePrime);
            _flicker = new NoiseGain("1overf", _drainPrime, _sourcePrime);
        }

        /// <inheritdoc />
        void INoiseSource.Initialize()
        {
            _rd.Initialize();
            _rs.Initialize();
            _id.Initialize();
            _flicker.Initialize();
        }

        /// <inheritdoc />
        void INoiseBehavior.Compute()
        {
            double vds;

            double T1, T10, T11;
            double Ssi, Swi;

            double m = Parameters.M;
            _rd.Compute(DrainConductance * m, _temperature.Temperature);
            _rs.Compute(SourceConductance * m, _temperature.Temperature);
            switch (ModelParameters.NoiMod.Value)
            {
                case 1:
                case 3:
                    _id.Compute(2.0 * Math.Abs(this._gm + this._gds + this._gmbs) / 3.0 * m,
                        _temperature.Temperature);
                    break;
                case 5:
                case 6:
                    vds = Math.Min(this._vds, this._vdsat);
                    _id.Compute((3.0 - vds / this._vdsat)
                             * Math.Abs(this._gm + this._gds + this._gmbs) / 3.0 * m,
                             _temperature.Temperature);
                    break;
                case 2:
                case 4:
                    _id.Compute((m * this._ueff * Math.Abs(this._qinv) /
                        (Param.BSIM3leff * Param.BSIM3leff + this._ueff * Math.Abs(this._qinv) * this._rds)),
                        _temperature.Temperature);    /* bugfix */
                    break;
            }
            switch (ModelParameters.NoiMod.Value)
            {
                case 1:
                case 4:
                case 5:
                    _flicker.Compute(m * ModelParameters.Kf * 
                        Math.Exp(ModelParameters.Af * Math.Log(Math.Max(Math.Abs(this._cd), 1e-38))) /
                        (Math.Pow(_state.Point.Value.Frequency, ModelParameters.Ef) * Param.BSIM3leff * Param.BSIM3leff * ModelTemperature.Cox));
                    break;
                case 2:
                case 3:
                case 6:
                    vds = this._vds;
                    if (vds < 0.0)
                        vds = -vds;
                    Ssi = StrongInversionNoiseEval(vds, _state.Point.Value.Frequency, _temperature.Temperature);
                    T10 = ModelParameters.OxideTrapDensityA
                    * 8.62e-5 * _temperature.Temperature;
                    T11 = Param.BSIM3weff * Param.BSIM3leff * Math.Pow(_state.Point.Value.Frequency, ModelParameters.Ef) * 4.0e36;
                    Swi = T10 / T11 * this._cd * this._cd;
                    T1 = Swi + Ssi;
                    if (T1 > 0.0)
                        _flicker.Compute(m * (Ssi * Swi) / T1);
                    else
                        _flicker.Compute(0);
                    break;
            }
        }

        private double StrongInversionNoiseEval(double Vds, double freq, double temp)
        {
            double cd, esat, DelClm, EffFreq, N0, Nl, Leff, Leffsq;
            double T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, Ssi;

            cd = Math.Abs(this._cd);
            Leff = Param.BSIM3leff - 2.0 * ModelParameters.Lintnoi;
            Leffsq = Leff * Leff;
            esat = 2.0 * Param.BSIM3vsattemp / this._ueff;
            if (ModelParameters.Em <= 0.0) DelClm = 0.0;
            else
            {
                T0 = ((((Vds - this._vdseff) / Param.BSIM3litl)
                    + ModelParameters.Em) / esat);
                DelClm = Param.BSIM3litl * Math.Log(Math.Max(T0, 1e-38));
                if (DelClm < 0.0) DelClm = 0.0;  /* bugfix */
            }
            EffFreq = Math.Pow(freq, ModelParameters.Ef);
            T1 = Constants.Charge * Constants.Charge * 8.62e-5 * cd * temp * this._ueff;
            T2 = 1.0e8 * EffFreq * this._abulk * ModelTemperature.Cox * Leffsq;
            N0 = ModelTemperature.Cox * this._vgsteff / Constants.Charge;
            Nl = ModelTemperature.Cox * this._vgsteff
                  * (1.0 - this._abovVgst2Vtm * this._vdseff) / Constants.Charge;

            T3 = ModelParameters.OxideTrapDensityA
               * Math.Log(Math.Max(((N0 + 2.0e14) / (Nl + 2.0e14)), 1e-38));
            T4 = ModelParameters.OxideTrapDensityB * (N0 - Nl);
            T5 = ModelParameters.OxideTrapDensityC * 0.5 * (N0 * N0 - Nl * Nl);

            T6 = 8.62e-5 * temp * cd * cd;
            T7 = 1.0e8 * EffFreq * Leffsq * Param.BSIM3weff;
            T8 = ModelParameters.OxideTrapDensityA + ModelParameters.OxideTrapDensityB * Nl
               + ModelParameters.OxideTrapDensityC * Nl * Nl;
            T9 = (Nl + 2.0e14) * (Nl + 2.0e14);

            Ssi = T1 / T2 * (T3 + T4 + T5) + T6 / T7 * DelClm * T8 / T9;
            return Ssi;
        }
    }
}
