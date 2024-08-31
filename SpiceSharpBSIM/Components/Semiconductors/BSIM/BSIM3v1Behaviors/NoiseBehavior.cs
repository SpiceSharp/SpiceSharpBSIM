using SpiceSharp;
using SpiceSharp.Attributes;
using SpiceSharp.Behaviors;
using SpiceSharp.Components;
using SpiceSharp.Components.NoiseSources;
using SpiceSharp.Simulations;
using System;

namespace SpiceSharpBSIM.Components.Semiconductors.BSIM.BSIM3v1Behaviors
{
    /// <summary>
    /// Noise behavior for a <see cref="BSIM3v1"/>.
    /// </summary>
    [BehaviorFor(typeof(BSIM3v1)), AddBehaviorIfNo(typeof(INoiseBehavior))]
    public class NoiseBehavior : FrequencyBehavior, INoiseBehavior
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
        void INoiseBehavior.Load() { }

        /// <inheritdoc />
        void INoiseBehavior.Compute()
        {
            double vgs, vds, Slimit;
            double T1, T10, T11;
            double Ssi, Swi;

            _rd.Compute(this._drainConductance * Parameters.M, _temperature.Temperature);
            _rs.Compute(this._sourceConductance * Parameters.M, _temperature.Temperature);

            switch (ModelParameters.NoiMod.Value)
            {
                case 1:
                case 3:
                    _id.Compute(2.0 / 3.0 * Math.Abs(this._gm + this._gds + this._gmbs) * Parameters.M, _temperature.Temperature);
                    break;
                case 2:
                case 4:
                    _id.Compute(this._ueff * Math.Abs((this._qinv * Parameters.M) / (Param.BSIM3v1leff * Param.BSIM3v1leff)), _temperature.Temperature);
                    break;
            }

            switch (ModelParameters.NoiMod.Value)
            {
                case 1:
                case 4:
                    _flicker.Compute(ModelParameters.Kf * Math.Exp(ModelParameters.Af
                  * Math.Log(Math.Max(Math.Abs(this._cd * Parameters.M), 1e-38))) /
                  (Math.Pow(_state.Point.Value.Frequency, ModelParameters.Ef) * Param.BSIM3v1leff * Param.BSIM3v1leff * ModelTemperature.Cox));
                    break;
                case 2:
                case 3:
                    vgs = this._vgs;
                    vds = this._vds;
                    if (vds < 0.0)
                    {
                        vds = -vds;
                        vgs = vgs + vds;
                    }
                    if (vgs >= this._von + 0.1)
                    {
                        Ssi = StrongInversionNoiseEval_b3(vgs, vds, _state.Point.Value.Frequency, _temperature.Temperature);
                        _flicker.Compute(Ssi);
                    }
                    else
                    {
                        T10 = ModelParameters.OxideTrapDensityA
                            * 8.62e-5 * _temperature.Temperature;
                        T11 = Param.BSIM3v1weff * Parameters.M * Param.BSIM3v1leff
                            * Math.Pow(_state.Point.Value.Frequency, ModelParameters.Ef)
                            * 4.0e36;
                        Swi = T10 / T11 * this._cd * Parameters.M
                            * this._cd * Parameters.M;
                        Slimit = StrongInversionNoiseEval_b3(this._von + 0.1, vds, _state.Point.Value.Frequency, _temperature.Temperature);
                        T1 = Swi + Slimit;
                        if (T1 > 0.0)
                            _flicker.Compute((Slimit * Swi) / T1);
                        else
                            _flicker.Compute(0.0);
                    }
                    break;
            }
        }

        private double StrongInversionNoiseEval_b3(double vgs, double vds, double freq, double temp)
        {
            double cd, esat, DelClm, EffFreq, N0, Nl, Vgst;
            double T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, Ssi;

            cd = Math.Abs(this._cd) * Parameters.M;
            if (vds > this._vdsat)
            {
                esat = 2.0 * Param.BSIM3v1vsattemp / this._ueff;
                T0 = ((((vds - this._vdsat) / Param.BSIM3v1litl) + ModelParameters.Em)
                   / esat);
                DelClm = Param.BSIM3v1litl * Math.Log(Math.Max(T0, 1e-38));
            }
            else
                DelClm = 0.0;
            EffFreq = Math.Pow(freq, ModelParameters.Ef);
            T1 = Constants.Charge * Constants.Charge * 8.62e-5 * cd * temp * this._ueff;
            T2 = 1.0e8 * EffFreq * ModelTemperature.Cox
               * Param.BSIM3v1leff * Param.BSIM3v1leff;
            Vgst = vgs - this._von;
            N0 = ModelTemperature.Cox * Vgst / Constants.Charge;
            if (N0 < 0.0)
                N0 = 0.0;
            Nl = ModelTemperature.Cox * (Vgst - Math.Min(vds, this._vdsat)) / Constants.Charge;
            if (Nl < 0.0)
                Nl = 0.0;

            T3 = ModelParameters.OxideTrapDensityA
               * Math.Log(Math.Max(((N0 + 2.0e14) / (Nl + 2.0e14)), 1e-38));
            T4 = ModelParameters.OxideTrapDensityB * (N0 - Nl);
            T5 = ModelParameters.OxideTrapDensityC * 0.5 * (N0 * N0 - Nl * Nl);

            T6 = 8.62e-5 * temp * cd * cd;
            T7 = 1.0e8 * EffFreq * Param.BSIM3v1leff
               * Param.BSIM3v1leff * Param.BSIM3v1weff * Parameters.M;
            T8 = ModelParameters.OxideTrapDensityA + ModelParameters.OxideTrapDensityB * Nl
               + ModelParameters.OxideTrapDensityC * Nl * Nl;
            T9 = (Nl + 2.0e14) * (Nl + 2.0e14);

            Ssi = T1 / T2 * (T3 + T4 + T5) + T6 / T7 * DelClm * T8 / T9;
            return Ssi;
        }
    }
}
