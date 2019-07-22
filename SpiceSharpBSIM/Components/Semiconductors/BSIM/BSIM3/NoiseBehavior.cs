using System;
using SpiceSharp.Behaviors;
using SpiceSharp.Components.NoiseSources;
using SpiceSharp.Simulations;
using SpiceSharp.Simulations.Behaviors;

namespace SpiceSharp.Components.BSIM3Behaviors
{
    /// <summary>
    /// Noise behavior for a <see cref="BSIM3"/>
    /// </summary>
    public class NoiseBehavior : ExportingBehavior, INoiseBehavior, IConnectedBehavior
    {
        private const double N_MINLOG = 1e-38;

        /// <summary>
        /// Necessary parameters and behaviors
        /// </summary>
        private BaseParameters _bp;
        private ModelBaseParameters _mbp;
        private BiasingBehavior _load;
        private TemperatureBehavior _temp;

        /// <summary>
        /// Nodes
        /// </summary>
        private int _drainNode, _gateNode, _sourceNode, _bulkNode, _sourceNodePrime, _drainNodePrime;

        /// <summary>
        /// Get the transistor noise sources
        /// </summary>
        public ComponentNoise TransistorNoise { get; } = new ComponentNoise(
            new NoiseThermal("rd", 0, 4),
            new NoiseThermal("rs", 1, 5),
            new NoiseThermal("id", 4, 5),
            new NoiseGain("1overf", 4, 5)
        );

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="name">Name</param>
        public NoiseBehavior(string name) : base(name)
        {
        }

        /// <summary>
        /// Setup behavior
        /// </summary>
        /// <param name="provider"></param>
        public override void Setup(Simulation simulation, SetupDataProvider provider)
        {
            if (provider == null)
                throw new ArgumentNullException(nameof(provider));

            // Get parameter sets
            _bp = provider.GetParameterSet<BaseParameters>();
            _mbp = provider.GetParameterSet<ModelBaseParameters>();

            // Get behaviors
            _temp = provider.GetBehavior<TemperatureBehavior>();
            _load = provider.GetBehavior<BiasingBehavior>();
        }

        /// <summary>
        /// Connect the behavior
        /// </summary>
        /// <param name="pins">Pins</param>
        public void Connect(params int[] pins)
        {
            _drainNode = pins[0];
            _gateNode = pins[1];
            _sourceNode = pins[2];
            _bulkNode = pins[3];
        }

        /// <summary>
        /// Connect the noise sources
        /// </summary>
        public void ConnectNoise()
        {
            _drainNodePrime = _load.DrainNodePrime;
            _sourceNodePrime = _load.SourceNodePrime;
            TransistorNoise.Setup(_drainNode, _gateNode, _sourceNode, _bulkNode, _drainNodePrime, _sourceNodePrime);
        }

        /// <summary>
        /// Noise behavior
        /// </summary>
        /// <param name="simulation"></param>
        public void Noise(Noise simulation)
        {
            double vds;
            var state = simulation.NoiseState;
            var pParam = _temp.Param;
            var m = _bp.Multiplier;

            TransistorNoise.Generators[0].SetCoefficients(_temp.DrainConductance);
            TransistorNoise.Generators[1].SetCoefficients(_temp.SourceConductance);

            switch (_mbp.NoiMod)
            {
                case 1:
                case 3:
                    TransistorNoise.Generators[2]
                        .SetCoefficients(2.0 * Math.Abs(_load.Gm + _load.Gds + _load.Gmbs) / 3.0 * m);
                    break;
                case 5:
                case 6:
                    vds = Math.Min(_load.Vds, _load.Vdsat);
                    TransistorNoise.Generators[2]
                        .SetCoefficients((3.0 - vds / _load.Vdsat) * Math.Abs(_load.Gm + _load.Gds + _load.Gmbs) / 3.0 *
                                         m);
                    break;
                case 2:
                case 4:
                    TransistorNoise.Generators[2]
                        .SetCoefficients(_load.Ueff * Math.Abs(_load.Qinv) /
                                         (pParam.BSIM3leff * pParam.BSIM3leff +
                                          _load.Ueff * Math.Abs(_load.Qinv) * _load.Rds) * m);
                    break;
                default:
                    CircuitWarning.Warning(this, "{0}: Invalid noise mode".FormatString(Name));
                    break;
            }

            switch (_mbp.NoiMod)
            {
                case 1:
                case 4:
                case 5:
                    TransistorNoise.Generators[3]
                        .SetCoefficients(
                            _mbp.Kf * Math.Exp(_mbp.Af * Math.Log(Math.Max(Math.Abs(_load.Cd), N_MINLOG))) /
                            (Math.Pow(state.Frequency, _mbp.Ef) * pParam.BSIM3leff * pParam.BSIM3leff * _mbp.Cox) * m);
                    break;
                case 2:
                case 3:
                case 6:
                    vds = _load.Vds;
                    if (vds < 0.0)
                        vds = -vds;
                    var ssi = StrongInversionNoiseEval(vds, state.Frequency, simulation.RealState.Temperature);
                    var t10 = _mbp.OxideTrapDensityA * 8.62e-5 * simulation.RealState.Temperature;
                    var t11 = pParam.BSIM3weff * pParam.BSIM3leff * Math.Pow(state.Frequency, _mbp.Ef) * 4.0e36;
                    var swi = t10 / t11 * _load.Cd * _load.Cd;
                    var t1 = swi + ssi;
                    if (t1 > 0.0)
                        TransistorNoise.Generators[3].SetCoefficients(ssi * swi / t1 * m);
                    else
                        TransistorNoise.Generators[3].SetCoefficients(0.0);
                    break;
            }

            // Evaluate noise sources
            TransistorNoise.Evaluate(simulation);
        }

        /// <summary>
        /// Strong inversion noise evaluation
        /// </summary>
        /// <param name="vds">Vds</param>
        /// <param name="freq">Frequency</param>
        /// <param name="temp">Temperature</param>
        /// <returns></returns>
        private double StrongInversionNoiseEval(double vds, double freq, double temp)
        {
            double delClm;

            var pParam = _temp.Param;
            var cd = Math.Abs(_load.Cd);
            var leff = pParam.BSIM3leff - 2.0 * _mbp.Lintnoi;
            var leffsq = leff * leff;
            var esat = 2.0 * pParam.BSIM3vsattemp / _load.Ueff;
            if (_mbp.Em <= 0.0)
                delClm = 0.0;
            else
            {
                var t0 = ((vds - _load.Vdseff) / pParam.BSIM3litl + _mbp.Em) / esat;
                delClm = pParam.BSIM3litl * Math.Log(Math.Max(t0, N_MINLOG));
                if (delClm < 0.0) delClm = 0.0; /* bugfix */
            }

            var effFreq = Math.Pow(freq, _mbp.Ef);
            var t1 = Constants.Charge * Constants.Charge * 8.62e-5 * cd * temp * _load.Ueff;
            var t2 = 1.0e8 * effFreq * _load.Abulk * _mbp.Cox * leffsq;
            var n0 = _mbp.Cox * _load.Vgsteff / Constants.Charge;
            var nl = _mbp.Cox * _load.Vgsteff * (1.0 - _load.AbovVgst2Vtm * _load.Vdseff) / Constants.Charge;

            var t3 = _mbp.OxideTrapDensityA * Math.Log(Math.Max((n0 + 2.0e14) / (nl + 2.0e14), N_MINLOG));
            var t4 = _mbp.OxideTrapDensityB * (n0 - nl);
            var t5 = _mbp.OxideTrapDensityC * 0.5 * (n0 * n0 - nl * nl);

            var t6 = 8.62e-5 * temp * cd * cd;
            var t7 = 1.0e8 * effFreq * leffsq * pParam.BSIM3weff;
            var t8 = _mbp.OxideTrapDensityA + _mbp.OxideTrapDensityB * nl + _mbp.OxideTrapDensityC * nl * nl;
            var t9 = (nl + 2.0e14) * (nl + 2.0e14);

            var ssi = t1 / t2 * (t3 + t4 + t5) + t6 / t7 * delClm * t8 / t9;
            return ssi;
        }
    }
}