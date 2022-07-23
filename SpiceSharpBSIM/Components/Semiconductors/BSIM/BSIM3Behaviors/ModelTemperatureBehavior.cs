using System;
using System.Collections.Generic;
using SpiceSharp;
using SpiceSharp.Attributes;
using SpiceSharp.Behaviors;
using SpiceSharp.Components;
using SpiceSharp.Entities;
using SpiceSharp.ParameterSets;
using SpiceSharp.Simulations;

namespace SpiceSharpBSIM.Components.Semiconductors.BSIM.BSIM3Behaviors
{
    /// <summary>
    /// Temperature behavior for a <see cref="BSIM3Model" />
    /// </summary>
    [BehaviorFor(typeof(BSIM3Model)), AddBehaviorIfNo(typeof(ITemperatureBehavior))]
    public class ModelTemperatureBehavior : Behavior, ITemperatureBehavior, IParameterized<ModelParameters>
    {
        private readonly ITemperatureSimulationState _temperature;

        public const double Kb = 1.3806226e-23;
        public const double KboQ = 8.617087e-5;
        public const double EPSSI = 1.03594e-10;
        public const double EPSOX = 3.453133e-11;

        /// <inheritdoc />
        public ModelParameters Parameters { get; }

        /// <summary>
        /// Size-dependent parameters cache
        /// </summary>
        public Dictionary<Tuple<double, double>, SizeDependParams> SizeDependParams { get; } = new Dictionary<Tuple<double, double>, SizeDependParams>();

        /// <summary>
        /// Properties
        /// </summary>
        public double Vcrit { get; private set; }
        public double Factor1 { get; private set; }
        public double Vtm { get; private set; }
        public double Vtm0 { get; private set; }
        public double JctTempSatCurDensity { get; private set; }
        public double JctSidewallTempSatCurDensity { get; private set; }
        public double UnitAreaTempJctCap { get; private set; }
        public double UnitLengthSidewallTempJctCap { get; private set; }
        public double UnitLengthGateSidewallTempJctCap { get; private set; }
        public double PhiB { get; private set; }
        public double PhiBSW { get; private set; }
        public double PhiBSWG { get; private set; }
        public double Ni { get; private set; }
        public double TRatio { get; private set; }
        public double Cox { get; private set; }
        public double T0 { get; private set; }

        /// <summary>
        /// Constructor
        /// </summary>
        public ModelTemperatureBehavior(BindingContext context)
            : base(context)
        {
            _temperature = context.GetState<ITemperatureSimulationState>();
            Parameters = context.GetParameterSet<ModelParameters>();
        }

        /// <summary>
        /// Temperature behavior
        /// </summary>
        void ITemperatureBehavior.Temperature()
        {
            double T1, Tnom, delTemp, Eg, Eg0;
            SizeDependParams.Clear();

            // Things done during setup that we need to put here
            if (Parameters.NqsMod.Given)
            {
                if (Parameters.NqsMod != 0 && Parameters.NqsMod != 1)
                {
                    Parameters.NqsMod = 0;
                    SpiceSharpWarning.Warning(this, "Warning: nqsMod has been set to its default value: 0.");
                }
            }
            if (Parameters.AcnqsMod.Given)
            {
                if ((Parameters.AcnqsMod != 0) && (Parameters.AcnqsMod != 1))
                {
                    Parameters.AcnqsMod = 0;
                    SpiceSharpWarning.Warning(this, "Warning: acnqsMod has been set to its default value: 0.");
                }
            }
            Cox = 3.453133e-11 / Parameters.Tox;
            if (!Parameters.Toxm.Given)
                Parameters.Toxm = new GivenParameter<double>(Parameters.Tox, false);

            if (!Parameters.Dsub.Given)
                Parameters.Dsub = new GivenParameter<double>(Parameters.Drout, false);
            if (!Parameters.Vth0.Given)
                Parameters.Vth0 = new GivenParameter<double>(Parameters.B3Type > 0 ? 0.7 : -0.7, false);

            if (!Parameters.Uc.Given)
                Parameters.Uc = new GivenParameter<double>((Parameters.MobMod.Value == 3) ? -0.0465 : -0.0465e-9, false);
            if (!Parameters.Uc1.Given)
                Parameters.Uc1 = new GivenParameter<double>((Parameters.MobMod.Value == 3) ? -0.056 : -0.056e-9, false);
            if (!Parameters.U0.Given)
                Parameters.U0 = new GivenParameter<double>((Parameters.B3Type > 0) ? 0.067 : 0.025, false);

            /* unit degree celcius */
            if (!Parameters.Tnom.Given)
                Parameters.Tnom = new GivenParameter<double>(_temperature.NominalTemperature, false);
            /*        else
                        Parameters.tnom = Parameters.tnom + 273.15; we make this transform in b3mpar.c in the first run */

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
                               * Math.Log(1.0 + 0.4e-6 / Parameters.Tox), false);
            if (!Parameters.Cgdo.Given)
            {
                if (Parameters.Dlc.Given && (Parameters.Dlc > 0.0))
                {
                    Parameters.Cgdo = Parameters.Dlc * Cox
                                     - Parameters.Cgdl;
                }
                else
                    Parameters.Cgdo = 0.6 * Parameters.Xj * Cox;
            }
            if (!Parameters.Cgso.Given)
            {
                if (Parameters.Dlc.Given && (Parameters.Dlc > 0.0))
                {
                    Parameters.Cgso = Parameters.Dlc * Cox
                                     - Parameters.Cgsl;
                }
                else
                    Parameters.Cgso = 0.6 * Parameters.Xj * Cox;
            }

            if (!Parameters.Cgbo.Given)
            {
                Parameters.Cgbo = 2.0 * Parameters.Dwc * Cox;
            }
            
            if (!Parameters.UnitLengthGateSidewallJctCap.Given)
                Parameters.UnitLengthGateSidewallJctCap = new GivenParameter<double>(Parameters.UnitLengthSidewallJctCap, false);
            if (!Parameters.GatesidewallJctPotential.Given)
                Parameters.GatesidewallJctPotential = new GivenParameter<double>(Parameters.SidewallJctPotential, false);
            if (!Parameters.BulkJctGateSideGradingCoeff.Given)
                Parameters.BulkJctGateSideGradingCoeff = new GivenParameter<double>(Parameters.BulkJctSideGradingCoeff, false);

            if (!Parameters.OxideTrapDensityA.Given)
            {
                if (Parameters.B3Type > 0)
                    Parameters.OxideTrapDensityA = 1e20;
                else
                    Parameters.OxideTrapDensityA = 9.9e18;
            }
            if (!Parameters.OxideTrapDensityB.Given)
            {
                if (Parameters.B3Type > 0)
                    Parameters.OxideTrapDensityB = 5e4;
                else
                    Parameters.OxideTrapDensityB = 2.4e3;
            }
            if (!Parameters.OxideTrapDensityC.Given)
            {
                if (Parameters.B3Type > 0)
                    Parameters.OxideTrapDensityC = -1.4e-12;
                else
                    Parameters.OxideTrapDensityC = 1.4e-12;

            }

            // Original temperature behavior
            Tnom = Parameters.Tnom;
            TRatio = _temperature.Temperature / Tnom;

            Vcrit = Constants.Vt0 * Math.Log(Constants.Vt0 / (Constants.Root2 * 1.0e-14));
            Factor1 = Math.Sqrt(EPSSI / EPSOX * Parameters.Tox);

            Vtm0 = KboQ * Tnom;
            Eg0 = 1.16 - 7.02e-4 * Tnom * Tnom / (Tnom + 1108.0);
            Ni = 1.45e10 * (Tnom / 300.15) * Math.Sqrt(Tnom / 300.15)
               * Math.Exp(21.5565981 - Eg0 / (2.0 * Vtm0));

            Vtm = KboQ * _temperature.Temperature;
            Eg = 1.16 - 7.02e-4 * _temperature.Temperature * _temperature.Temperature / (_temperature.Temperature + 1108.0);
            if (_temperature.Temperature != Tnom)
            {
                T0 = Eg0 / Vtm0 - Eg / Vtm + Parameters.JctTempExponent
                   * Math.Log(_temperature.Temperature / Tnom);
                T1 = Math.Exp(T0 / Parameters.JctEmissionCoeff);
                JctTempSatCurDensity = Parameters.JctSatCurDensity
                                                 * T1;
                JctSidewallTempSatCurDensity
                            = Parameters.JctSidewallSatCurDensity * T1;
            }
            else
            {
                JctTempSatCurDensity = Parameters.JctSatCurDensity;
                JctSidewallTempSatCurDensity
                           = Parameters.JctSidewallSatCurDensity;
            }

            if (JctTempSatCurDensity < 0.0)
                JctTempSatCurDensity = 0.0;
            if (JctSidewallTempSatCurDensity < 0.0)
                JctSidewallTempSatCurDensity = 0.0;

            /* Temperature dependence of D/B and S/B diode capacitance begins */
            delTemp = _temperature.Temperature - Parameters.Tnom;
            T0 = Parameters.Tcj * delTemp;
            if (T0 >= -1.0)
            {
                UnitAreaTempJctCap = Parameters.UnitAreaJctCap * (1.0 + T0);
            }
            else if (Parameters.UnitAreaJctCap > 0.0)
            {
                UnitAreaTempJctCap = 0.0;
                SpiceSharpWarning.Warning(this, "Temperature effect has caused cj to be negative.Cj is clamped to zero.");
            }
            T0 = Parameters.Tcjsw * delTemp;
            if (T0 >= -1.0)
            {
                UnitLengthSidewallTempJctCap = Parameters.UnitLengthSidewallJctCap * (1.0 + T0);
            }
            else if (Parameters.UnitLengthSidewallJctCap > 0.0)
            {
                UnitLengthSidewallTempJctCap = 0.0;
                SpiceSharpWarning.Warning(this, "Temperature effect has caused cjsw to be negative. Cjsw is clamped to zero.");
            }
            T0 = Parameters.Tcjswg * delTemp;
            if (T0 >= -1.0)
            {
                UnitLengthGateSidewallTempJctCap = Parameters.UnitLengthGateSidewallJctCap * (1.0 + T0);
            }
            else if (Parameters.UnitLengthGateSidewallJctCap > 0.0)
            {
                UnitLengthGateSidewallTempJctCap = 0.0;
                SpiceSharpWarning.Warning(this, "Temperature effect has caused cjswg to be negative. Cjswg is clamped to zero.");
            }

            PhiB = Parameters.BulkJctPotential
                             - Parameters.Tpb * delTemp;
            if (PhiB < 0.01)
            {
                PhiB = 0.01;
                SpiceSharpWarning.Warning(this, "Temperature effect has caused pb to be less than 0.01. Pb is clamped to 0.01.");
            }
            PhiBSW = Parameters.SidewallJctPotential
                               - Parameters.Tpbsw * delTemp;
            if (PhiBSW <= 0.01)
            {
                PhiBSW = 0.01;
                SpiceSharpWarning.Warning(this, "Temperature effect has caused pbsw to be less than 0.01. Pbsw is clamped to 0.01.");
            }
            PhiBSWG = Parameters.GatesidewallJctPotential
                                - Parameters.Tpbswg * delTemp;
            if (PhiBSWG <= 0.01)
            {
                PhiBSWG = 0.01;
                SpiceSharpWarning.Warning(this, "Temperature effect has caused pbswg to be less than 0.01. Pbswg is clamped to 0.01.");
            }
            /* End of junction capacitance */
        }
    }
}