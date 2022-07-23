using SpiceSharp.Behaviors;
using SpiceSharp.Entities;
using SpiceSharp.ParameterSets;
using SpiceSharp.Simulations;
using SpiceSharp;
using System;
using System.Collections.Generic;
using SpiceSharp.Attributes;
using SpiceSharp.Components;

namespace SpiceSharpBSIM.Components.Semiconductors.BSIM.BSIM3v2Behaviors
{
    /// <summary>
    /// Temperature behavior for a <see cref="BSIM3v2Model"/>.
    /// </summary>
    [BehaviorFor(typeof(BSIM3v2Model)), AddBehaviorIfNo(typeof(ITemperatureBehavior))]
    public class ModelTemperatureBehavior : Behavior, ITemperatureBehavior, IParameterized<ModelParameters>
    {
        private readonly ITemperatureSimulationState _temperature;

        public const double KboQ = 8.617087e-5;
        public const double MAX_EXP = 5.834617425e14;
        public const double MIN_EXP = 1.713908431e-15;
        public const double EXP_THRESHOLD = 34.0;
        public const double SMOOTHFACTOR = 0.1;
        public const double EPSOX = 3.453133e-11;
        public const double EPSSI = 1.03594e-10;
        public const double Meter2Micron = 1.0e6;

        /// <inheritdoc />
        public ModelParameters Parameters { get; }

        /// <summary>
        /// Gets the size-dependent parameters.
        /// </summary>
        public Dictionary<Tuple<double, double>, SizeDependParams> SizeDependParams { get; } = new Dictionary<Tuple<double, double>, SizeDependParams>();

        public Version IntVersion { get; private set; }
        public double Cox { get; private set; }
        public double TRatio { get; private set; }
        public double Vcrit { get; private set; }
        public double Factor1 { get; private set; }
        public double Vtm0 { get; private set; }
        public double Ni { get; private set; }
        public double Vtm { get; private set; }
        public double JctTempSatCurDensity { get; private set; }
        public double JctSidewallTempSatCurDensity { get; private set; }
        public double UnitAreaTempJctCap { get; private set; }
        public double UnitLengthSidewallTempJctCap { get; private set; }
        public double UnitLengthGateSidewallTempJctCap { get; private set; }
        public double PhiB { get; private set; }
        public double PhiBSW { get; private set; }
        public double PhiBSWG { get; private set; }
        public double T0 { get; private set; }

        /// <summary>
        /// Creates a new <see cref="ModelTemperatureBehavior"/>.
        /// </summary>
        /// <param name="context">The context.</param>
        public ModelTemperatureBehavior(BindingContext context)
            : base(context)
        {
            _temperature = context.GetState<ITemperatureSimulationState>();
            Parameters = context.GetParameterSet<ModelParameters>();
        }

        /// <inheritdoc />
        void ITemperatureBehavior.Temperature()
        {
            double Eg, Eg0, T0, T1;
            double delTemp, Temp, Tnom;
            SizeDependParams.Clear();

            /* Default value Processing for BSIM3v32 MOSFET Models */
            if (!Parameters.NqsMod.Given)
                Parameters.NqsMod = 0;
            else if ((Parameters.NqsMod != 0) && (Parameters.NqsMod != 1))
            {
                Parameters.NqsMod = 0;
                SpiceSharpWarning.Warning(this, "Warning: nqsMod has been set to its default value: 0.");
            }
            this.Cox = 3.453133e-11 / Parameters.Tox;
            if (!Parameters.Toxm.Given)
                Parameters.Toxm = new GivenParameter<double>(Parameters.Tox, false);

            if (!Parameters.Dsub.Given)
                Parameters.Dsub = new GivenParameter<double>(Parameters.Drout, false);
            if (!Parameters.Vth0.Given)
                Parameters.Vth0 = new GivenParameter<double>(Parameters.Type > 0 ? 0.7 : -0.7, false);

            if (!Parameters.Uc.Given)
                Parameters.Uc = new GivenParameter<double>((Parameters.MobMod.Value == 3) ? -0.0465 : -0.0465e-9, false);
            if (!Parameters.Uc1.Given)
                Parameters.Uc1 = new GivenParameter<double>((Parameters.MobMod.Value == 3) ? -0.056 : -0.056e-9, false);
            if (!Parameters.U0.Given)
                Parameters.U0 = new GivenParameter<double>((Parameters.Type > 0) ? 0.067 : 0.025, false);

            if (!Parameters.Tnom.Given)
                Parameters.Tnom = new GivenParameter<double>(_temperature.NominalTemperature, false);

            if (!Parameters.Llc.Given)
                Parameters.Llc = new GivenParameter<double>(Parameters.Ll, false);

            if (!Parameters.Lwc.Given)
                Parameters.Lwc = new GivenParameter<double>(Parameters.Lw, false);
            if (!Parameters.Lwlc.Given)
                Parameters.Lwlc = new GivenParameter<double>(Parameters.Lwl, false);

            if (!Parameters.Wlc.Given)
                Parameters.Wlc = new GivenParameter<double>(Parameters.Wl, false);
            if (!Parameters.Wwc.Given)
                Parameters.Wwc = new GivenParameter<double>(Parameters.Ww, false);
            if (!Parameters.Wwlc.Given)
                Parameters.Wwlc = new GivenParameter<double>(Parameters.Wwl, false);
            if (!Parameters.Dwc.Given)
                Parameters.Dwc = new GivenParameter<double>(Parameters.Wint, false);
            if (!Parameters.Dlc.Given)
                Parameters.Dlc = new GivenParameter<double>(Parameters.Lint, false);

            if (!Parameters.Cf.Given)
                Parameters.Cf = new GivenParameter<double>(2.0 * EPSOX / Math.PI
                                   * Math.Log(1.0 + 0.4e-6 / Parameters.Tox.Value), false);

            /* If the user does not provide the model revision,
             * we always choose the most recent.
             */

            /* I have added below the code that translate model string
             * into an integer. This trick is meant to speed up the
             * revision testing instruction, since comparing integer
             * is faster than comparing strings.
             * Paolo Nenzi 2002
             */
            if (Parameters.Version == "3.2.4" || Parameters.Version == "3.24")
                this.IntVersion = Version.BSIM3v32V324;
            else if (Parameters.Version == "3.2.3" || Parameters.Version == "3.23")
                this.IntVersion = Version.BSIM3v32V323;
            else if (Parameters.Version == "3.2.2" || Parameters.Version == "3.22")
                this.IntVersion = Version.BSIM3v32V322;
            else if (Parameters.Version == "3.2" || Parameters.Version == "3.20")
                this.IntVersion = Version.BSIM3v32V32;
            else
                this.IntVersion = Version.BSIM3v32V3OLD;

            /* BSIM3v32V3OLD is a placeholder for pre 3.2 revision
             * This model should not be used for pre 3.2 models.
             */

            if (!Parameters.Cgdo.Given)
            {
                if (Parameters.Dlc.Given && (Parameters.Dlc > 0.0))
                {
                    Parameters.Cgdo = Parameters.Dlc * this.Cox
                                     - Parameters.Cgdl;
                }
                else
                    Parameters.Cgdo = 0.6 * Parameters.Xj * this.Cox;
            }
            if (!Parameters.Cgso.Given)
            {
                if (Parameters.Dlc.Given && (Parameters.Dlc > 0.0))
                {
                    Parameters.Cgso = Parameters.Dlc * this.Cox
                                     - Parameters.Cgsl;
                }
                else
                    Parameters.Cgso = 0.6 * Parameters.Xj * this.Cox;
            }

            if (!Parameters.Cgbo.Given)
            {
                Parameters.Cgbo = 2.0 * Parameters.Dwc * this.Cox;
            }

            if (!Parameters.UnitLengthGateSidewallJctCap.Given)
                Parameters.UnitLengthGateSidewallJctCap = new GivenParameter<double>(Parameters.UnitLengthSidewallJctCap, false);
            if (!Parameters.GatesidewallJctPotential.Given)
                Parameters.GatesidewallJctPotential = new GivenParameter<double>(Parameters.SidewallJctPotential, false);
            if (!Parameters.BulkJctGateSideGradingCoeff.Given)
                Parameters.BulkJctGateSideGradingCoeff = new GivenParameter<double>(Parameters.BulkJctSideGradingCoeff, false);

            if (!Parameters.OxideTrapDensityA.Given)
            {
                if (Parameters.Type > 0)
                    Parameters.OxideTrapDensityA = 1e20;
                else
                    Parameters.OxideTrapDensityA = 9.9e18;
            }
            if (!Parameters.OxideTrapDensityB.Given)
            {
                if (Parameters.Type > 0)
                    Parameters.OxideTrapDensityB = 5e4;
                else
                    Parameters.OxideTrapDensityB = 2.4e3;
            }
            if (!Parameters.OxideTrapDensityC.Given)
            {
                if (Parameters.Type > 0)
                    Parameters.OxideTrapDensityC = -1.4e-12;
                else
                    Parameters.OxideTrapDensityC = 1.4e-12;
            }

            Temp = _temperature.Temperature;
            Tnom = Parameters.Tnom;
            this.TRatio = Temp / Tnom;

            this.Vcrit = Constants.Vt0 * Math.Log(Constants.Vt0 / (Constants.Root2 * 1.0e-14));
            this.Factor1 = Math.Sqrt(EPSSI / EPSOX * Parameters.Tox);

            this.Vtm0 = KboQ * Tnom;
            Eg0 = 1.16 - 7.02e-4 * Tnom * Tnom / (Tnom + 1108.0);
            this.Ni = 1.45e10 * (Tnom / 300.15) * Math.Sqrt(Tnom / 300.15)
               * Math.Exp(21.5565981 - Eg0 / (2.0 * this.Vtm0));

            this.Vtm = KboQ * Temp;
            Eg = 1.16 - 7.02e-4 * Temp * Temp / (Temp + 1108.0);
            if (Temp != Tnom)
            {
                T0 = Eg0 / this.Vtm0 - Eg / this.Vtm + Parameters.JctTempExponent
                   * Math.Log(Temp / Tnom);
                T1 = Math.Exp(T0 / Parameters.JctEmissionCoeff);
                this.JctTempSatCurDensity = Parameters.JctSatCurDensity
                                                 * T1;
                this.JctSidewallTempSatCurDensity
                                     = Parameters.JctSidewallSatCurDensity * T1;
            }
            else
            {
                this.JctTempSatCurDensity = Parameters.JctSatCurDensity;
                this.JctSidewallTempSatCurDensity
                           = Parameters.JctSidewallSatCurDensity;
            }

            if (this.JctTempSatCurDensity < 0.0)
                this.JctTempSatCurDensity = 0.0;
            if (this.JctSidewallTempSatCurDensity < 0.0)
                this.JctSidewallTempSatCurDensity = 0.0;

            /* Temperature dependence of D/B and S/B diode capacitance begins */
            delTemp = _temperature.Temperature - Parameters.Tnom;
            T0 = Parameters.Tcj * delTemp;
            if (T0 >= -1.0)
            {
                /* Added revision dependent code */
                switch (this.IntVersion)
                {
                    case Version.BSIM3v32V324:
                    case Version.BSIM3v32V323:
                        this.UnitAreaTempJctCap =
                                Parameters.UnitAreaJctCap * (1.0 + T0);
                        break;
                    case Version.BSIM3v32V322:
                    case Version.BSIM3v32V32:
                    default:
                        Parameters.UnitAreaJctCap *= 1.0 + T0;
                        break;
                }
            }
            else if (Parameters.UnitAreaJctCap > 0.0)
            {
                /* Added revision dependent code */
                switch (this.IntVersion)
                {
                    case Version.BSIM3v32V324:
                    case Version.BSIM3v32V323:
                        this.UnitAreaTempJctCap = 0.0;
                        break;
                    case Version.BSIM3v32V322:
                    case Version.BSIM3v32V32:
                    default:
                        Parameters.UnitAreaJctCap = 0.0;
                        break;
                }
                SpiceSharpWarning.Warning(this, "Temperature effect has caused cj to be negative. Cj is clamped to zero.");
            }
            T0 = Parameters.Tcjsw * delTemp;
            if (T0 >= -1.0)
            {
                /* Added revision dependent code */
                switch (this.IntVersion)
                {
                    case Version.BSIM3v32V324:
                    case Version.BSIM3v32V323:
                        this.UnitLengthSidewallTempJctCap =
                                Parameters.UnitLengthSidewallJctCap * (1.0 + T0);
                        break;
                    case Version.BSIM3v32V322:
                    case Version.BSIM3v32V32:
                    default:
                        Parameters.UnitLengthSidewallJctCap *= 1.0 + T0;
                        break;
                }
            }
            else if (Parameters.UnitLengthSidewallJctCap > 0.0)
            {
                /* Added revision dependent code */
                switch (this.IntVersion)
                {
                    case Version.BSIM3v32V324:
                    case Version.BSIM3v32V323:
                        this.UnitLengthSidewallTempJctCap = 0.0;
                        break;
                    case Version.BSIM3v32V322:
                    case Version.BSIM3v32V32:
                    default:
                        Parameters.UnitLengthSidewallJctCap = 0.0;
                        break;
                }
                SpiceSharpWarning.Warning(this, "Temperature effect has caused cjsw to be negative. Cjsw is clamped to zero.");
            }
            T0 = Parameters.Tcjswg * delTemp;
            if (T0 >= -1.0)
            {
                /* Added revision dependent code */
                switch (this.IntVersion)
                {
                    case Version.BSIM3v32V324:
                    case Version.BSIM3v32V323:
                        this.UnitLengthGateSidewallTempJctCap =
                                Parameters.UnitLengthGateSidewallJctCap * (1.0 + T0);
                        break;
                    case Version.BSIM3v32V322:
                    case Version.BSIM3v32V32:
                    default:
                        Parameters.UnitLengthGateSidewallJctCap *= 1.0 + T0;
                        break;
                }
            }
            else if (Parameters.UnitLengthGateSidewallJctCap > 0.0)
            {
                /* Added revision dependent code */
                switch (this.IntVersion)
                {
                    case Version.BSIM3v32V324:
                    case Version.BSIM3v32V323:
                        this.UnitLengthGateSidewallTempJctCap = 0.0;
                        break;
                    case Version.BSIM3v32V322:
                    case Version.BSIM3v32V32:
                    default:
                        Parameters.UnitLengthGateSidewallJctCap = 0.0;
                        break;
                }
                SpiceSharpWarning.Warning(this, "Temperature effect has caused cjswg to be negative. Cjswg is clamped to zero.");
            }

            this.PhiB = Parameters.BulkJctPotential
                             - Parameters.Tpb * delTemp;
            if (this.PhiB < 0.01)
            {
                this.PhiB = 0.01;
                SpiceSharpWarning.Warning(this, "Temperature effect has caused pb to be less than 0.01. Pb is clamped to 0.01.");
            }
            this.PhiBSW = Parameters.SidewallJctPotential
                               - Parameters.Tpbsw * delTemp;
            if (this.PhiBSW <= 0.01)
            {
                this.PhiBSW = 0.01;
                SpiceSharpWarning.Warning(this, "Temperature effect has caused pbsw to be less than 0.01. Pbsw is clamped to 0.01.");
            }
            this.PhiBSWG = Parameters.GatesidewallJctPotential
                                - Parameters.Tpbswg * delTemp;
            if (this.PhiBSWG <= 0.01)
            {
                this.PhiBSWG = 0.01;
                SpiceSharpWarning.Warning(this, "Temperature effect has caused pbswg to be less than 0.01. Pbswg is clamped to 0.01.");
            }
            /* End of junction capacitance */

            // Why is this even necessary...
            this.T0 = T0;
        }
    }
}
