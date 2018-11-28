using SpiceSharp.Attributes;

namespace SpiceSharp.Components.BSIM1Behaviors
{
    /// <summary>
    /// Base parameters for a <see cref="BSIM1Model" />
    /// </summary>
    public class ModelBaseParameters : ParameterSet
    {
        /*
         * Please note that Spice 3f5 did not use default parameters. We will add default values as described here:
         * http://literature.cdn.keysight.com/litweb/pdf/ads2001/spicent/snt0814.html
         */

        /// <summary>
        /// Properties
        /// </summary>
        public double Type { get; private set; } = 1;
        [ParameterName("vfb"), ParameterInfo("Flat band voltage")]
        public GivenParameter<double> Vfb0 { get; } = new GivenParameter<double>(-0.3);
        [ParameterName("lvfb"), ParameterInfo("Length dependence of vfb")]
        public GivenParameter<double> VfbL { get; } = new GivenParameter<double>();
        [ParameterName("wvfb"), ParameterInfo("Width dependence of vfb")]
        public GivenParameter<double> VfbW { get; } = new GivenParameter<double>();
        [ParameterName("phi"), ParameterInfo("Strong inversion surface potential ")]
        public GivenParameter<double> Phi0 { get; } = new GivenParameter<double>(0.6);
        [ParameterName("lphi"), ParameterInfo("Length dependence of phi")]
        public GivenParameter<double> PhiL { get; } = new GivenParameter<double>();
        [ParameterName("wphi"), ParameterInfo("Width dependence of phi")]
        public GivenParameter<double> PhiW { get; } = new GivenParameter<double>();
        [ParameterName("k1"), ParameterInfo("Bulk effect coefficient 1")]
        public GivenParameter<double> K10 { get; } = new GivenParameter<double>(0.5);
        [ParameterName("lk1"), ParameterInfo("Length dependence of k1")]
        public GivenParameter<double> K1L { get; } = new GivenParameter<double>();
        [ParameterName("wk1"), ParameterInfo("Width dependence of k1")]
        public GivenParameter<double> K1W { get; } = new GivenParameter<double>();
        [ParameterName("k2"), ParameterInfo("Bulk effect coefficient 2")]
        public GivenParameter<double> K20 { get; } = new GivenParameter<double>();
        [ParameterName("lk2"), ParameterInfo("Length dependence of k2")]
        public GivenParameter<double> K2L { get; } = new GivenParameter<double>();
        [ParameterName("wk2"), ParameterInfo("Width dependence of k2")]
        public GivenParameter<double> K2W { get; } = new GivenParameter<double>();
        [ParameterName("eta"), ParameterInfo("VDS dependence of threshold voltage")]
        public GivenParameter<double> Eta0 { get; } = new GivenParameter<double>();
        [ParameterName("leta"), ParameterInfo("Length dependence of eta")]
        public GivenParameter<double> EtaL { get; } = new GivenParameter<double>();
        [ParameterName("weta"), ParameterInfo("Width dependence of eta")]
        public GivenParameter<double> EtaW { get; } = new GivenParameter<double>();
        [ParameterName("x2e"), ParameterInfo("VBS dependence of eta")]
        public GivenParameter<double> EtaB0 { get; } = new GivenParameter<double>(-0.07);
        [ParameterName("lx2e"), ParameterInfo("Length dependence of x2e")]
        public GivenParameter<double> EtaBl { get; } = new GivenParameter<double>();
        [ParameterName("wx2e"), ParameterInfo("Width dependence of x2e")]
        public GivenParameter<double> EtaBw { get; } = new GivenParameter<double>();
        [ParameterName("x3e"), ParameterInfo("VDS dependence of eta")]
        public GivenParameter<double> EtaD0 { get; } = new GivenParameter<double>();
        [ParameterName("lx3e"), ParameterInfo("Length dependence of x3e")]
        public GivenParameter<double> EtaDl { get; } = new GivenParameter<double>();
        [ParameterName("wx3e"), ParameterInfo("Width dependence of x3e")]
        public GivenParameter<double> EtaDw { get; } = new GivenParameter<double>();
        [ParameterName("dl"), ParameterInfo("Channel length reduction in um")]
        public GivenParameter<double> DeltaL { get; } = new GivenParameter<double>();
        [ParameterName("dw"), ParameterInfo("Channel width reduction in um")]
        public GivenParameter<double> DeltaW { get; } = new GivenParameter<double>();
        [ParameterName("muz"), ParameterInfo("Zero field mobility at VDS=0 VGS=VTH")]
        public GivenParameter<double> MobZero { get; } = new GivenParameter<double>(600);
        [ParameterName("x2mz"), ParameterInfo("VBS dependence of muz")]
        public GivenParameter<double> MobZeroB0 { get; } = new GivenParameter<double>();
        [ParameterName("lx2mz"), ParameterInfo("Length dependence of x2mz")]
        public GivenParameter<double> MobZeroBl { get; } = new GivenParameter<double>();
        [ParameterName("wx2mz"), ParameterInfo("Width dependence of x2mz")]
        public GivenParameter<double> MobZeroBw { get; } = new GivenParameter<double>();
        [ParameterName("mus"), ParameterInfo("Mobility at VDS=VDD VGS=VTH, channel length modulation")]
        public GivenParameter<double> MobVdd0 { get; } = new GivenParameter<double>(1082);
        [ParameterName("lmus"), ParameterInfo("Length dependence of mus")]
        public GivenParameter<double> MobVddl { get; } = new GivenParameter<double>();
        [ParameterName("wmus"), ParameterInfo("Width dependence of mus")]
        public GivenParameter<double> MobVddw { get; } = new GivenParameter<double>();
        [ParameterName("x2ms"), ParameterInfo("VBS dependence of mus")]
        public GivenParameter<double> MobVddB0 { get; } = new GivenParameter<double>();
        [ParameterName("lx2ms"), ParameterInfo("Length dependence of x2ms")]
        public GivenParameter<double> MobVddBl { get; } = new GivenParameter<double>();
        [ParameterName("wx2ms"), ParameterInfo("Width dependence of x2ms")]
        public GivenParameter<double> MobVddBw { get; } = new GivenParameter<double>();
        [ParameterName("x3ms"), ParameterInfo("VDS dependence of mus")]
        public GivenParameter<double> MobVddD0 { get; } = new GivenParameter<double>();
        [ParameterName("lx3ms"), ParameterInfo("Length dependence of x3ms")]
        public GivenParameter<double> MobVddDl { get; } = new GivenParameter<double>();
        [ParameterName("wx3ms"), ParameterInfo("Width dependence of x3ms")]
        public GivenParameter<double> MobVddDw { get; } = new GivenParameter<double>();
        [ParameterName("u0"), ParameterInfo("VGS dependence of mobility")]
        public GivenParameter<double> Ugs0 { get; } = new GivenParameter<double>(670.0);
        [ParameterName("lu0"), ParameterInfo("Length dependence of u0")]
        public GivenParameter<double> UgsL { get; } = new GivenParameter<double>();
        [ParameterName("wu0"), ParameterInfo("Width dependence of u0")]
        public GivenParameter<double> UgsW { get; } = new GivenParameter<double>();
        [ParameterName("x2u0"), ParameterInfo("VBS dependence of u0")]
        public GivenParameter<double> UgsB0 { get; } = new GivenParameter<double>();
        [ParameterName("lx2u0"), ParameterInfo("Length dependence of x2u0")]
        public GivenParameter<double> UgsBL { get; } = new GivenParameter<double>();
        [ParameterName("wx2u0"), ParameterInfo("Width dependence of x2u0")]
        public GivenParameter<double> UgsBW { get; } = new GivenParameter<double>();
        [ParameterName("u1"), ParameterInfo("VDS depence of mobility, velocity saturation")]
        public GivenParameter<double> Uds0 { get; } = new GivenParameter<double>();
        [ParameterName("lu1"), ParameterInfo("Length dependence of u1")]
        public GivenParameter<double> UdsL { get; } = new GivenParameter<double>();
        [ParameterName("wu1"), ParameterInfo("Width dependence of u1")]
        public GivenParameter<double> UdsW { get; } = new GivenParameter<double>();
        [ParameterName("x2u1"), ParameterInfo("VBS depence of u1")]
        public GivenParameter<double> UdsB0 { get; } = new GivenParameter<double>();
        [ParameterName("lx2u1"), ParameterInfo("Length depence of x2u1")]
        public GivenParameter<double> UdsBL { get; } = new GivenParameter<double>();
        [ParameterName("wx2u1"), ParameterInfo("Width depence of x2u1")]
        public GivenParameter<double> UdsBW { get; } = new GivenParameter<double>();
        [ParameterName("x3u1"), ParameterInfo("VDS depence of u1")]
        public GivenParameter<double> UdsD0 { get; } = new GivenParameter<double>();
        [ParameterName("lx3u1"), ParameterInfo("Length dependence of x3u1")]
        public GivenParameter<double> UdsDL { get; } = new GivenParameter<double>();
        [ParameterName("wx3u1"), ParameterInfo("Width depence of x3u1")]
        public GivenParameter<double> UdsDW { get; } = new GivenParameter<double>();
        [ParameterName("n0"), ParameterInfo("Subthreshold slope")]
        public GivenParameter<double> SubthSlope0 { get; } = new GivenParameter<double>(0.5);
        [ParameterName("ln0"), ParameterInfo("Length dependence of n0")]
        public GivenParameter<double> SubthSlopeL { get; } = new GivenParameter<double>();
        [ParameterName("wn0"), ParameterInfo("Width dependence of n0")]
        public GivenParameter<double> SubthSlopeW { get; } = new GivenParameter<double>();
        [ParameterName("nb"), ParameterInfo("VBS dependence of subthreshold slope")]
        public GivenParameter<double> SubthSlopeB0 { get; } = new GivenParameter<double>();
        [ParameterName("lnb"), ParameterInfo("Length dependence of nb")]
        public GivenParameter<double> SubthSlopeBL { get; } = new GivenParameter<double>();
        [ParameterName("wnb"), ParameterInfo("Width dependence of nb")]
        public GivenParameter<double> SubthSlopeBW { get; } = new GivenParameter<double>();
        [ParameterName("nd"), ParameterInfo("VDS dependence of subthreshold slope")]
        public GivenParameter<double> SubthSlopeD0 { get; } = new GivenParameter<double>();
        [ParameterName("lnd"), ParameterInfo("Length dependence of nd")]
        public GivenParameter<double> SubthSlopeDL { get; } = new GivenParameter<double>();
        [ParameterName("wnd"), ParameterInfo("Width dependence of nd")]
        public GivenParameter<double> SubthSlopeDW { get; } = new GivenParameter<double>();
        [ParameterName("tox"), ParameterInfo("Gate oxide thickness in um")]
        public GivenParameter<double> OxideThickness { get; } = new GivenParameter<double>(1e-7);
        [ParameterName("temp"), ParameterInfo("Temperature in degree Celcius")]
        public GivenParameter<double> Temp { get; } = new GivenParameter<double>(25);
        [ParameterName("vdd"), ParameterInfo("Supply voltage to specify mus")]
        public GivenParameter<double> Vdd { get; } = new GivenParameter<double>(5.0);
        [ParameterName("cgso"), ParameterInfo("Gate source overlap capacitance per unit channel width(m)")]
        public GivenParameter<double> GateSourceOverlapCap { get; } = new GivenParameter<double>();
        [ParameterName("cgdo"), ParameterInfo("Gate drain overlap capacitance per unit channel width(m)")]
        public GivenParameter<double> GateDrainOverlapCap { get; } = new GivenParameter<double>();
        [ParameterName("cgbo"), ParameterInfo("Gate bulk overlap capacitance per unit channel length(m)")]
        public GivenParameter<double> GateBulkOverlapCap { get; } = new GivenParameter<double>();
        [ParameterName("xpart"), ParameterInfo("Flag for channel charge partitioning")]
        public GivenParameter<double> ChannelChargePartitionFlag { get; } = new GivenParameter<double>(1.0);
        [ParameterName("rsh"), ParameterInfo("Source drain diffusion sheet resistance in ohm per square")]
        public GivenParameter<double> SheetResistance { get; } = new GivenParameter<double>();
        [ParameterName("js"), ParameterInfo("Source drain junction saturation current per unit area")]
        public GivenParameter<double> JctSatCurDensity { get; } = new GivenParameter<double>();
        [ParameterName("pb"), ParameterInfo("Source drain junction built in potential")]
        public GivenParameter<double> BulkJctPotential { get; } = new GivenParameter<double>(0.8);
        [ParameterName("mj"), ParameterInfo("Source drain bottom junction capacitance grading coefficient")]
        public GivenParameter<double> BulkJctBotGradingCoeff { get; } = new GivenParameter<double>(0.5);
        [ParameterName("pbsw"), ParameterInfo("Source drain side junction capacitance built in potential")]
        public GivenParameter<double> SidewallJctPotential { get; } = new GivenParameter<double>(1.0);
        [ParameterName("mjsw"), ParameterInfo("Source drain side junction capacitance grading coefficient")]
        public GivenParameter<double> BulkJctSideGradingCoeff { get; } = new GivenParameter<double>(0.33);
        [ParameterName("cj"), ParameterInfo("Source drain bottom junction capacitance per unit area")]
        public GivenParameter<double> UnitAreaJctCap { get; } = new GivenParameter<double>();
        [ParameterName("cjsw"), ParameterInfo("Source drain side junction capacitance per unit area")]
        public GivenParameter<double> UnitLengthSidewallJctCap { get; } = new GivenParameter<double>();
        [ParameterName("wdf"), ParameterInfo("Default width of source drain diffusion in um")]
        public GivenParameter<double> DefaultWidth { get; } = new GivenParameter<double>();
        [ParameterName("dell"), ParameterInfo("Length reduction of source drain diffusion")]
        public GivenParameter<double> DeltaLength { get; } = new GivenParameter<double>();
        [ParameterName("nmos"), ParameterInfo("Flag to indicate NMOS")]
        public void SetNMOS(bool flag = true)
        {
            if (flag)
                Type = 1;
        }
        [ParameterName("pmos"), ParameterInfo("Flag to indicate PMOS")]
        public void SetPMOS(bool flag = true)
        {
            if (flag)
                Type = -1;
        }

        /// <summary>
        /// Gets the oxide capacitance.
        /// </summary>
        /// <value>
        /// The oxide capacitance.
        /// </value>
        public double Cox { get; private set; }

        /// <summary>
        /// Deep cloning
        /// </summary>
        /// <returns></returns>
        public override ParameterSet DeepClone()
        {
            var cloned = (ModelBaseParameters) base.DeepClone();
            cloned.Type = Type;
            return cloned;
        }

        /// <summary>
        /// Calculates the defaults.
        /// </summary>
        public override void CalculateDefaults()
        {
            if (BulkJctPotential < 0.1)
                BulkJctPotential.Value = 0.1;
            if (SidewallJctPotential < 0.1)
                SidewallJctPotential.Value = 0.1;
            Cox = 3.453e-13 / (OxideThickness * 1.0e-4);
        }
    }
}