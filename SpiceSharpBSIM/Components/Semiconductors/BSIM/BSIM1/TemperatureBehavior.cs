﻿using System;
using SpiceSharp.Behaviors;
using SpiceSharp.Simulations;
using SpiceSharp.Simulations.Behaviors;

namespace SpiceSharp.Components.BSIM1Behaviors
{
    /// <summary>
    /// Temperature behavior for a <see cref="BSIM1" />
    /// </summary>
    public class TemperatureBehavior : ExportingBehavior, ITemperatureBehavior
    {
        /// <summary>
        /// Necessary behaviors and parameters
        /// </summary>
        protected BaseParameters BaseParameters { get; private set; }
        protected ModelBaseParameters ModelParameters { get; private set; }

        /// <summary>
        /// Properties
        /// </summary>
        public double GDoverlapCap { get; private set; }
        public double GSoverlapCap { get; private set; }
        public double GBoverlapCap { get; private set; }
        public double DrainConductance { get; private set; }
        public double SourceConductance { get; private set; }
        public double Vfb { get; private set; }
        public double Phi { get; private set; }
        public double K1 { get; private set; }
        public double K2 { get; private set; }
        public double Eta { get; private set; }
        public double EtaB { get; private set; }
        public double EtaD { get; private set; }
        public double BetaZero { get; private set; }
        public double BetaZeroB { get; private set; }
        public double Ugs { get; private set; }
        public double UgsB { get; private set; }
        public double Uds { get; private set; }
        public double UdsB { get; private set; }
        public double UdsD { get; private set; }
        public double BetaVdd { get; private set; }
        public double BetaVddB { get; private set; }
        public double BetaVddD { get; private set; }
        public double SubthSlope { get; private set; }
        public double SubthSlopeB { get; private set; }
        public double SubthSlopeD { get; private set; }
        public double Vt0 { get; private set; }
        public double Von { get; set; }

        /// <summary>
        /// Constructor
        /// </summary>
        public TemperatureBehavior(string name) : base(name)
        {

        }

        /// <summary>
        /// Setup the behavior
        /// </summary>
        public override void Setup(Simulation simulation, SetupDataProvider provider)
        {
            if (provider == null)
                throw new ArgumentNullException(nameof(provider));

            // Get parameters
            BaseParameters = provider.GetParameterSet<BaseParameters>();
            ModelParameters = provider.GetParameterSet<ModelBaseParameters>("model");
        }

        /// <summary>
        /// Temperature behavior
        /// </summary>
        public void Temperature(BaseSimulation simulation)
        {
            double effChanLength, effChanWidth, coxWoverL, leff, weff;
            if ((effChanLength = BaseParameters.Length - ModelParameters.DeltaL * 1e-6) <= 0)
            {
                throw new CircuitException("B1: mosfet {0}: Effective channel length <= 0".FormatString(Name));
            }
            if ((effChanWidth = BaseParameters.Width - ModelParameters.DeltaW * 1e-6) <= 0)
            {
                throw new CircuitException("B1: mosfet {0}: Effective channel width <= 0".FormatString(Name));
            }
            GDoverlapCap = effChanWidth * ModelParameters.GateDrainOverlapCap;
            GSoverlapCap = effChanWidth * ModelParameters.GateSourceOverlapCap;
            GBoverlapCap = BaseParameters.Length * ModelParameters.GateBulkOverlapCap;
            if ((DrainConductance = ModelParameters.SheetResistance * BaseParameters.DrainSquares) != 0.0)
            {
                DrainConductance = 1.0 / DrainConductance;
            }
            if ((SourceConductance = ModelParameters.SheetResistance * BaseParameters.SourceSquares) != 0.0)
            {
                SourceConductance = 1.0 / SourceConductance;
            }
            leff = effChanLength * 1.0e6;
            weff = effChanWidth * 1.0e6;
            coxWoverL = ModelParameters.Cox * weff / leff;
            Vfb = ModelParameters.Vfb0 + ModelParameters.VfbL / leff + ModelParameters.VfbW / weff;
            Phi = ModelParameters.Phi0 + ModelParameters.PhiL / leff + ModelParameters.PhiW / weff;
            K1 = ModelParameters.K10 + ModelParameters.K1L / leff + ModelParameters.K1W / weff;
            K2 = ModelParameters.K20 + ModelParameters.K2L / leff + ModelParameters.K2W / weff;
            Eta = ModelParameters.Eta0 + ModelParameters.EtaL / leff + ModelParameters.EtaW / weff;
            EtaB = ModelParameters.EtaB0 + ModelParameters.EtaBl / leff + ModelParameters.EtaBw / weff;
            EtaD = ModelParameters.EtaD0 + ModelParameters.EtaDl / leff + ModelParameters.EtaDw / weff;
            BetaZero = ModelParameters.MobZero;
            BetaZeroB = ModelParameters.MobZeroB0 + ModelParameters.MobZeroBl / leff + ModelParameters.MobZeroBw / weff;
            Ugs = ModelParameters.Ugs0 + ModelParameters.UgsL / leff + ModelParameters.UgsW / weff;
            UgsB = ModelParameters.UgsB0 + ModelParameters.UgsBL / leff + ModelParameters.UgsBW / weff;
            Uds = ModelParameters.Uds0 + ModelParameters.UdsL / leff + ModelParameters.UdsW / weff;
            UdsB = ModelParameters.UdsB0 + ModelParameters.UdsBL / leff + ModelParameters.UdsBW / weff;
            UdsD = ModelParameters.UdsD0 + ModelParameters.UdsDL / leff + ModelParameters.UdsDW / weff;
            BetaVdd = ModelParameters.MobVdd0 + ModelParameters.MobVddl / leff + ModelParameters.MobVddw / weff;
            BetaVddB = ModelParameters.MobVddB0 + ModelParameters.MobVddBl / leff + ModelParameters.MobVddBw / weff;
            BetaVddD = ModelParameters.MobVddD0 + ModelParameters.MobVddDl / leff + ModelParameters.MobVddDw / weff;
            SubthSlope = ModelParameters.SubthSlope0 + ModelParameters.SubthSlopeL / leff + ModelParameters.SubthSlopeW / weff;
            SubthSlopeB = ModelParameters.SubthSlopeB0 + ModelParameters.SubthSlopeBL / leff + ModelParameters.SubthSlopeBW / weff;
            SubthSlopeD = ModelParameters.SubthSlopeD0 + ModelParameters.SubthSlopeDL / leff + ModelParameters.SubthSlopeDW / weff;
            if (Phi < 0.1)
                Phi = 0.1;
            if (K1 < 0.0)
                K1 = 0.0;
            if (K2 < 0.0)
                K2 = 0.0;
            Vt0 = Vfb + Phi + K1 * Math.Sqrt(Phi) - K2 * Phi;
            Von = Vt0;
            BetaZero = BetaZero * coxWoverL;
            BetaZeroB = BetaZeroB * coxWoverL;
            BetaVdd = BetaVdd * coxWoverL;
            BetaVddB = BetaVddB * coxWoverL;
            BetaVddD = Math.Max(BetaVddD * coxWoverL, 0.0);
        }
    }
}