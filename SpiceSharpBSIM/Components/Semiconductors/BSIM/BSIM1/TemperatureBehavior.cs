using System;
using SpiceSharp.Behaviors;
using SpiceSharp.Simulations;
namespace SpiceSharp.Components.BSIM1Behaviors
{

    /// <summary>
    /// Temperature behavior for a <see cref="BSIM1" />
    /// </summary>
    public class TemperatureBehavior : BaseTemperatureBehavior
    {

        /// <summary>
        /// Necessary behaviors and parameters
        /// </summary>
        private BaseParameters _bp;
        private ModelBaseParameters _mbp;
        private ModelTemperatureBehavior _modelTemp;

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
        public TemperatureBehavior(Identifier name) : base(name)
        {

        }

        /// <summary>
        /// Setup the behavior
        /// </summary>
        public override void Setup(Simulation simulation, SetupDataProvider provider)
        {
            if (provider == null)
                throw new ArgumentNullException(nameof(provider));
            _bp = provider.GetParameterSet<BaseParameters>("entity");
            _mbp = provider.GetParameterSet<ModelBaseParameters>("model");
            _modelTemp = provider.GetBehavior<ModelTemperatureBehavior>("model");
        }

        /// <summary>
        /// Temperature behavior
        /// </summary>
        public override void Temperature(BaseSimulation simulation)
        {
            double effChanLength, effChanWidth, coxWoverL, leff, weff;
            if ((effChanLength = _bp.Length - _mbp.DeltaL * 1e-6) <= 0)
            {
                throw new CircuitException("B1: mosfet {0}, model {1}: Effective channel length <=0".FormatString(_modelTemp.Name, Name));
            }
            if ((effChanWidth = _bp.Width - _mbp.DeltaW * 1e-6) <= 0)
            {
                throw new CircuitException("B1: mosfet {0}, model {1}: Effective channel width <=0".FormatString(_modelTemp.Name, Name));
            }
            GDoverlapCap = effChanWidth * _mbp.GateDrainOverlapCap;
            GSoverlapCap = effChanWidth * _mbp.GateSourceOverlapCap;
            GBoverlapCap = _bp.Length * _mbp.GateBulkOverlapCap;
            if ((DrainConductance = _mbp.SheetResistance * _bp.DrainSquares) != 0.0)
            {
                DrainConductance = 1.0 / DrainConductance;
            }
            if ((SourceConductance = _mbp.SheetResistance * _bp.SourceSquares) != 0.0)
            {
                SourceConductance = 1.0 / SourceConductance;
            }
            leff = effChanLength * 1.0e6;
            weff = effChanWidth * 1.0e6;
            coxWoverL = _modelTemp.Cox * weff / leff;
            Vfb = _mbp.Vfb0 + _mbp.VfbL / leff + _mbp.VfbW / weff;
            Phi = _mbp.Phi0 + _mbp.PhiL / leff + _mbp.PhiW / weff;
            K1 = _mbp.K10 + _mbp.K1L / leff + _mbp.K1W / weff;
            K2 = _mbp.K20 + _mbp.K2L / leff + _mbp.K2W / weff;
            Eta = _mbp.Eta0 + _mbp.EtaL / leff + _mbp.EtaW / weff;
            EtaB = _mbp.EtaB0 + _mbp.EtaBl / leff + _mbp.EtaBw / weff;
            EtaD = _mbp.EtaD0 + _mbp.EtaDl / leff + _mbp.EtaDw / weff;
            BetaZero = _mbp.MobZero;
            BetaZeroB = _mbp.MobZeroB0 + _mbp.MobZeroBl / leff + _mbp.MobZeroBw / weff;
            Ugs = _mbp.Ugs0 + _mbp.UgsL / leff + _mbp.UgsW / weff;
            UgsB = _mbp.UgsB0 + _mbp.UgsBL / leff + _mbp.UgsBW / weff;
            Uds = _mbp.Uds0 + _mbp.UdsL / leff + _mbp.UdsW / weff;
            UdsB = _mbp.UdsB0 + _mbp.UdsBL / leff + _mbp.UdsBW / weff;
            UdsD = _mbp.UdsD0 + _mbp.UdsDL / leff + _mbp.UdsDW / weff;
            BetaVdd = _mbp.MobVdd0 + _mbp.MobVddl / leff + _mbp.MobVddw / weff;
            BetaVddB = _mbp.MobVddB0 + _mbp.MobVddBl / leff + _mbp.MobVddBw / weff;
            BetaVddD = _mbp.MobVddD0 + _mbp.MobVddDl / leff + _mbp.MobVddDw / weff;
            SubthSlope = _mbp.SubthSlope0 + _mbp.SubthSlopeL / leff + _mbp.SubthSlopeW / weff;
            SubthSlopeB = _mbp.SubthSlopeB0 + _mbp.SubthSlopeBL / leff + _mbp.SubthSlopeBW / weff;
            SubthSlopeD = _mbp.SubthSlopeD0 + _mbp.SubthSlopeDL / leff + _mbp.SubthSlopeDW / weff;
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