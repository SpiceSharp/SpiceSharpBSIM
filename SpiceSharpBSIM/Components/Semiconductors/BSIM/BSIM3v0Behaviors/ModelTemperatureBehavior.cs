using SpiceSharp;
using SpiceSharp.Attributes;
using SpiceSharp.Behaviors;
using SpiceSharp.Components;
using SpiceSharp.Entities;
using SpiceSharp.ParameterSets;
using SpiceSharp.Simulations;
using System;
using System.Collections.Generic;

namespace SpiceSharpBSIM.Components.Semiconductors.BSIM.BSIM3v0Behaviors
{
    /// <summary>
    /// Temperature behavior for a <see cref="BSIM3v0Model"/>.
    /// </summary>
    [BehaviorFor(typeof(BSIM3v0Model)), AddBehaviorIfNo(typeof(ITemperatureBehavior))]
    public class ModelTemperatureBehavior : Behavior, ITemperatureBehavior, IParameterized<ModelParameters>
    {
        private readonly ITemperatureSimulationState _temperature;

        public const double Kb = 1.3806226e-23;
        public const double KboQ = 8.617087e-5;  /* Kb / q  where q = 1.60219e-19 */
        public const double EPSOX = 3.453133e-11;
        public const double EPSSI = 1.03594e-10;
        public const double PI = 3.141592654;
        public const double MAX_EXP = 5.834617425e14;
        public const double MIN_EXP = 1.713908431e-15;
        public const double EXP_THRESHOLD = 34.0;
        public const double Charge_q = 1.60219e-19;

        /// <inheritdoc />
        public ModelParameters Parameters { get; }
        public double Cox { get; private set; }

        /// <summary>
        /// Gets the size-dependent properties.
        /// </summary>
        public Dictionary<Tuple<double, double>, SizeDependentProperties> SizeDependentProperties { get; } = new Dictionary<Tuple<double, double>, SizeDependentProperties>();
        public double Vcrit { get; private set; }
        public double Factor1 { get; private set; }
        public double Vtm { get; private set; }

        /// <summary>
        /// Creates a new <see cref="ModelTemperatureBehavior"/>.
        /// </summary>
        /// <param name="context">The context.</param>
        public ModelTemperatureBehavior(BindingContext context)
            : base(context)
        {
            _temperature = context.GetState<ITemperatureSimulationState>();
            Parameters = context.GetParameterSet<ModelParameters>();
            Setup();
        }

        private void Setup()
        {
            this.Cox = new GivenParameter<double>(3.453133e-11 / Parameters.Tox, false);
            if (!Parameters.Dsub.Given)
                Parameters.Dsub = new GivenParameter<double>(Parameters.Drout, false);
            if (!Parameters.Vth0.Given)
                Parameters.Vth0 = new GivenParameter<double>((Parameters.Type > 0) ? 0.7 : -0.7, false);
            if (!Parameters.Uc.Given)
                Parameters.Uc = new GivenParameter<double>((Parameters.MobMod == 3) ? -0.0465 : -0.0465e-9, false);
            if (!Parameters.Uc1.Given)
                Parameters.Uc1 = new GivenParameter<double>((Parameters.MobMod == 3) ? -0.056 : -0.056e-9, false);
            if (!Parameters.U0.Given)
                Parameters.U0 = new GivenParameter<double>((Parameters.Type > 0) ? 0.067 : 0.025, false);
            else if (Parameters.U0 > 1.0)
                Parameters.U0 = new GivenParameter<double>(Parameters.U0 / 1.0e4, false);
            if (!Parameters.Lu0.Given)
                Parameters.Lu0 = new GivenParameter<double>(0.0, false);
            else if (Parameters.U0 > 1.0)
                Parameters.Lu0 = new GivenParameter<double>(Parameters.Lu0 / 1.0e4, false);
            if (!Parameters.Wu0.Given)
                Parameters.Wu0 = new GivenParameter<double>(0.0, false);
            else if (Parameters.U0 > 1.0)
                Parameters.Wu0 = new GivenParameter<double>(Parameters.Wu0 / 1.0e4, false);
            if (!Parameters.Pu0.Given)
                Parameters.Pu0 = new GivenParameter<double>(0.0, false);
            else if (Parameters.U0 > 1.0)
                Parameters.Pu0 = new GivenParameter<double>(Parameters.Pu0 / 1.0e4, false);
            if (!Parameters.Tnom.Given)
                Parameters.Tnom = new GivenParameter<double>(_temperature.NominalTemperature, false);
            if (!Parameters.Dwc.Given)
                Parameters.Dwc = new GivenParameter<double>(Parameters.Wint, false);
            if (!Parameters.Dlc.Given)
                Parameters.Dlc = new GivenParameter<double>(Parameters.Lint, false);
            if (!Parameters.Cf.Given)
                Parameters.Cf = new GivenParameter<double>(2.0 * EPSOX / Math.PI
                   * Math.Log(1.0 + 0.4e-6 / Parameters.Tox), false);
            if (!Parameters.Cgdo.Given)
            {
                if (Parameters.Dlc.Given && (Parameters.Dlc > 0.0))
                {
                    Parameters.Cgdo = new GivenParameter<double>(Parameters.Dlc * this.Cox
                         - Parameters.Cgdl, false);
                    if (Parameters.Cgdo < 0.0)
                        Parameters.Cgdo = new GivenParameter<double>(0.0, false);
                }
                else
                    Parameters.Cgdo = new GivenParameter<double>(0.6 * Parameters.Xj * this.Cox, false);
            }
            if (!Parameters.Cgso.Given)
            {
                if (Parameters.Dlc.Given && (Parameters.Dlc > 0.0))
                {
                    Parameters.Cgso = new GivenParameter<double>(Parameters.Dlc * this.Cox
                         - Parameters.Cgsl, false);
                    if (Parameters.Cgso < 0.0)
                        Parameters.Cgso = new GivenParameter<double>(0.0, false);
                }
                else
                    Parameters.Cgso = new GivenParameter<double>(0.6 * Parameters.Xj * this.Cox, false);
            }
            if (!Parameters.OxideTrapDensityA.Given)
            {
                if (Parameters.Type > 0)
                    Parameters.OxideTrapDensityA = new GivenParameter<double>(1e20, false);
                else
                    Parameters.OxideTrapDensityA = new GivenParameter<double>(9.9e18, false);
            }
            if (!Parameters.OxideTrapDensityB.Given)
            {
                if (Parameters.Type > 0)
                    Parameters.OxideTrapDensityB = new GivenParameter<double>(5e4, false);
                else
                    Parameters.OxideTrapDensityB = new GivenParameter<double>(2.4e3, false);
            }
            if (!Parameters.OxideTrapDensityC.Given)
            {
                if (Parameters.Type > 0)
                    Parameters.OxideTrapDensityC = new GivenParameter<double>(-1.4e-12, false);
                else
                    Parameters.OxideTrapDensityC = new GivenParameter<double>(1.4e-12, false);
            }
        }

        /// <inheritdoc />
        void ITemperatureBehavior.Temperature()
        {
            SizeDependentProperties.Clear();
            this.Vcrit = Constants.Vt0 * Math.Log(Constants.Vt0
                  / (Constants.Root2 * 1.0e-14));
            this.Factor1 = Math.Sqrt(EPSSI / EPSOX * Parameters.Tox);
            this.Vtm = KboQ * _temperature.Temperature;
        }
    }
}
