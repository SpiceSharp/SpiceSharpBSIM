using SpiceSharp;
using SpiceSharp.Attributes;
using SpiceSharp.Behaviors;
using SpiceSharp.Components;
using SpiceSharp.Entities;
using SpiceSharp.ParameterSets;
using SpiceSharp.Simulations;
using System;
using System.Collections.Generic;

namespace SpiceSharpBSIM.Components.Semiconductors.BSIM.BSIM3v1Behaviors
{
    /// <summary>
    /// Temperature behavior for a <see cref="BSIM3v1Model"/>.
    /// </summary>
    [BehaviorFor(typeof(BSIM3v1Model)), AddBehaviorIfNo(typeof(ITemperatureBehavior))]
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
        public Dictionary<Tuple<double, double>, SizeDependentProperties> SizeDependentProperties { get; } = new Dictionary<Tuple<double, double>, SizeDependentProperties>();
        public double Vcrit { get; private set; }
        public double Factor1 { get; private set; }
        public double Ni { get; private set; }
        public double Vtm { get; private set; }
        public double JctTempSatCurDensity { get; private set; }
        public double JctSidewallTempSatCurDensity { get; private set; }
        public double Vtm0 { get; private set; }

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
            this.Cox = 3.453133e-11 / Parameters.Tox;
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
            if (!Parameters.Tnom.Given)
                Parameters.Tnom = new GivenParameter<double>(_temperature.NominalTemperature, false);
            if (!Parameters.Dwc.Given)
                Parameters.Dwc = new GivenParameter<double>(Parameters.Wint, false);
            if (!Parameters.Dlc.Given)
                Parameters.Dlc = new GivenParameter<double>(Parameters.Lint, false);
            if (!Parameters.Cf.Given)
                Parameters.Cf = new GivenParameter<double>(2.0 * EPSOX / PI
                   * Math.Log(1.0 + 0.4e-6 / Parameters.Tox), false);
            if (!Parameters.Cgdo.Given)
            {
                if (Parameters.Dlc.Given && (Parameters.Dlc > 0.0))
                {
                    Parameters.Cgdo = new GivenParameter<double>(Parameters.Dlc * this.Cox
                         - Parameters.Cgdl, false);
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
                }
                else
                    Parameters.Cgso = new GivenParameter<double>(0.6 * Parameters.Xj * this.Cox, false);
            }
            if (!Parameters.Cgbo.Given)
            {
                Parameters.Cgbo = new GivenParameter<double>(2.0 * Parameters.Dwc * this.Cox, false);
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
            double Eg, Eg0, T0, T1;
            double Temp, TRatio, Vtm0, Tnom;

            Temp = _temperature.Temperature;
            SizeDependentProperties.Clear();

            Tnom = Parameters.Tnom;
            TRatio = Temp / Tnom;

            this.Vcrit = Constants.Vt0 * Math.Log(Constants.Vt0 / (Constants.Root2 * 1.0e-14));
            this.Factor1 = Math.Sqrt(EPSSI / EPSOX * Parameters.Tox);

            Vtm0 = KboQ * Tnom;
            Eg0 = 1.16 - 7.02e-4 * Tnom * Tnom / (Tnom + 1108.0);
            this.Ni = 1.45e10 * (Tnom / 300.15) * Math.Sqrt(Tnom / 300.15)
               * Math.Exp(21.5565981 - Eg0 / (2.0 * Vtm0));

            this.Vtm = KboQ * Temp;
            Eg = 1.16 - 7.02e-4 * Temp * Temp / (Temp + 1108.0);
            if (Temp != Tnom)
            {
                T0 = Eg0 / Vtm0 - Eg / this.Vtm + Parameters.JctTempExponent
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

            this.Vtm0 = Vtm0;
        }
    }
}
