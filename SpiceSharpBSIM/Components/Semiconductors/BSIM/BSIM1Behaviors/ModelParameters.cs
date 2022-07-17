using SpiceSharp;
using SpiceSharp.Attributes;
using SpiceSharp.ParameterSets;

namespace SpiceSharpBSIM.Components.Semiconductors.BSIM.BSIM1Behaviors
{
    /// <summary>
    /// Base parameters for a <see cref="BSIM1Model" />
    /// </summary>
    [GeneratedParameters]
    public partial class ModelParameters : ParameterSet<ModelParameters>
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
        [Finite]
        private GivenParameter<double> _vfb0 = new GivenParameter<double>(-0.3);

        [ParameterName("lvfb"), ParameterInfo("Length dependence of vfb")]
        [Finite]
        private GivenParameter<double> _vfbL = new GivenParameter<double>();

        [ParameterName("wvfb"), ParameterInfo("Width dependence of vfb")]
        [Finite]
        private GivenParameter<double> _vfbW = new GivenParameter<double>();

        [ParameterName("phi"), ParameterInfo("Strong inversion surface potential ")]
        [Finite]
        private GivenParameter<double> _phi0 = new GivenParameter<double>(0.6);

        [ParameterName("lphi"), ParameterInfo("Length dependence of phi")]
        [Finite]
        private GivenParameter<double> _phiL = new GivenParameter<double>();

        [ParameterName("wphi"), ParameterInfo("Width dependence of phi")]
        [Finite]
        private GivenParameter<double> _phiW = new GivenParameter<double>();

        [ParameterName("k1"), ParameterInfo("Bulk effect coefficient 1")]
        [Finite]
        private GivenParameter<double> _k10 = new GivenParameter<double>(0.5);

        [ParameterName("lk1"), ParameterInfo("Length dependence of k1")]
        [Finite]
        private GivenParameter<double> _k1L = new GivenParameter<double>();

        [ParameterName("wk1"), ParameterInfo("Width dependence of k1")]
        [Finite]
        private GivenParameter<double> _k1W = new GivenParameter<double>();

        [ParameterName("k2"), ParameterInfo("Bulk effect coefficient 2")]
        [Finite]
        private GivenParameter<double> _k20 = new GivenParameter<double>();

        [ParameterName("lk2"), ParameterInfo("Length dependence of k2")]
        [Finite]
        private GivenParameter<double> _k2L = new GivenParameter<double>();

        [ParameterName("wk2"), ParameterInfo("Width dependence of k2")]
        [Finite]
        private GivenParameter<double> _k2W = new GivenParameter<double>();

        [ParameterName("eta"), ParameterInfo("VDS dependence of threshold voltage")]
        [Finite]
        private GivenParameter<double> _eta0 = new GivenParameter<double>();

        [ParameterName("leta"), ParameterInfo("Length dependence of eta")]
        [Finite]
        private GivenParameter<double> _etaL = new GivenParameter<double>();

        [ParameterName("weta"), ParameterInfo("Width dependence of eta")]
        [Finite]
        private GivenParameter<double> _etaW = new GivenParameter<double>();

        [ParameterName("x2e"), ParameterInfo("VBS dependence of eta")]
        [Finite]
        private GivenParameter<double> _etaB0 = new GivenParameter<double>(-0.07);

        [ParameterName("lx2e"), ParameterInfo("Length dependence of x2e")]
        [Finite]
        private GivenParameter<double> _etaBl = new GivenParameter<double>();

        [ParameterName("wx2e"), ParameterInfo("Width dependence of x2e")]
        [Finite]
        private GivenParameter<double> _etaBw = new GivenParameter<double>();

        [ParameterName("x3e"), ParameterInfo("VDS dependence of eta")]
        [Finite]
        private GivenParameter<double> _etaD0 = new GivenParameter<double>();

        [ParameterName("lx3e"), ParameterInfo("Length dependence of x3e")]
        [Finite]
        private GivenParameter<double> _etaDl = new GivenParameter<double>();

        [ParameterName("wx3e"), ParameterInfo("Width dependence of x3e")]
        [Finite]
        private GivenParameter<double> _etaDw = new GivenParameter<double>();

        [ParameterName("dl"), ParameterInfo("Channel length reduction in um")]
        [Finite]
        private GivenParameter<double> _deltaL = new GivenParameter<double>();

        [ParameterName("dw"), ParameterInfo("Channel width reduction in um")]
        [Finite]
        private GivenParameter<double> _deltaW = new GivenParameter<double>();

        [ParameterName("muz"), ParameterInfo("Zero field mobility at VDS=0 VGS=VTH")]
        [Finite, GreaterThan(0)]
        private GivenParameter<double> _mobZero = new GivenParameter<double>(600);

        [ParameterName("x2mz"), ParameterInfo("VBS dependence of muz")]
        [Finite]
        private GivenParameter<double> _mobZeroB0 = new GivenParameter<double>();

        [ParameterName("lx2mz"), ParameterInfo("Length dependence of x2mz")]
        [Finite]
        private GivenParameter<double> _mobZeroBl = new GivenParameter<double>();

        [ParameterName("wx2mz"), ParameterInfo("Width dependence of x2mz")]
        [Finite]
        private GivenParameter<double> _mobZeroBw = new GivenParameter<double>();

        [ParameterName("mus"), ParameterInfo("Mobility at VDS=VDD VGS=VTH, channel length modulation")]
        [Finite, GreaterThan(0)]
        private GivenParameter<double> _mobVdd0 = new GivenParameter<double>(1082);

        [ParameterName("lmus"), ParameterInfo("Length dependence of mus")]
        [Finite]
        private GivenParameter<double> _mobVddl = new GivenParameter<double>();

        [ParameterName("wmus"), ParameterInfo("Width dependence of mus")]
        [Finite]
        private GivenParameter<double> _mobVddw = new GivenParameter<double>();

        [ParameterName("x2ms"), ParameterInfo("VBS dependence of mus")]
        [Finite]
        private GivenParameter<double> _mobVddB0 = new GivenParameter<double>();

        [ParameterName("lx2ms"), ParameterInfo("Length dependence of x2ms")]
        [Finite]
        private GivenParameter<double> _mobVddBl = new GivenParameter<double>();

        [ParameterName("wx2ms"), ParameterInfo("Width dependence of x2ms")]
        [Finite]
        private GivenParameter<double> _mobVddBw = new GivenParameter<double>();

        [ParameterName("x3ms"), ParameterInfo("VDS dependence of mus")]
        [Finite]
        private GivenParameter<double> _mobVddD0 = new GivenParameter<double>();

        [ParameterName("lx3ms"), ParameterInfo("Length dependence of x3ms")]
        [Finite]
        private GivenParameter<double> _mobVddDl = new GivenParameter<double>();

        [ParameterName("wx3ms"), ParameterInfo("Width dependence of x3ms")]
        [Finite]
        private GivenParameter<double> _mobVddDw = new GivenParameter<double>();

        [ParameterName("u0"), ParameterInfo("VGS dependence of mobility")]
        [Finite]
        private GivenParameter<double> _ugs0 = new GivenParameter<double>(670.0);

        [ParameterName("lu0"), ParameterInfo("Length dependence of u0")]
        [Finite]
        private GivenParameter<double> _ugsL = new GivenParameter<double>();

        [ParameterName("wu0"), ParameterInfo("Width dependence of u0")]
        [Finite]
        private GivenParameter<double> _ugsW = new GivenParameter<double>();

        [ParameterName("x2u0"), ParameterInfo("VBS dependence of u0")]
        [Finite]
        private GivenParameter<double> _ugsB0 = new GivenParameter<double>();

        [ParameterName("lx2u0"), ParameterInfo("Length dependence of x2u0")]
        [Finite]
        private GivenParameter<double> _ugsBL = new GivenParameter<double>();

        [ParameterName("wx2u0"), ParameterInfo("Width dependence of x2u0")]
        [Finite]
        private GivenParameter<double> _ugsBW = new GivenParameter<double>();

        [ParameterName("u1"), ParameterInfo("VDS depence of mobility, velocity saturation")]
        [Finite]
        private GivenParameter<double> _uds0 = new GivenParameter<double>();

        [ParameterName("lu1"), ParameterInfo("Length dependence of u1")]
        [Finite]
        private GivenParameter<double> _udsL = new GivenParameter<double>();

        [ParameterName("wu1"), ParameterInfo("Width dependence of u1")]
        [Finite]
        private GivenParameter<double> _udsW = new GivenParameter<double>();

        [ParameterName("x2u1"), ParameterInfo("VBS depence of u1")]
        [Finite]
        private GivenParameter<double> _udsB0 = new GivenParameter<double>();

        [ParameterName("lx2u1"), ParameterInfo("Length depence of x2u1")]
        [Finite]
        private GivenParameter<double> _udsBL = new GivenParameter<double>();

        [ParameterName("wx2u1"), ParameterInfo("Width depence of x2u1")]
        [Finite]
        private GivenParameter<double> _udsBW = new GivenParameter<double>();

        [ParameterName("x3u1"), ParameterInfo("VDS depence of u1")]
        [Finite]
        private GivenParameter<double> _udsD0 = new GivenParameter<double>();

        [ParameterName("lx3u1"), ParameterInfo("Length dependence of x3u1")]
        [Finite]
        private GivenParameter<double> _udsDL = new GivenParameter<double>();

        [ParameterName("wx3u1"), ParameterInfo("Width depence of x3u1")]
        [Finite]
        private GivenParameter<double> _udsDW = new GivenParameter<double>();

        [ParameterName("n0"), ParameterInfo("Subthreshold slope")]
        [Finite]
        private GivenParameter<double> _subthSlope0 = new GivenParameter<double>(0.5);

        [ParameterName("ln0"), ParameterInfo("Length dependence of n0")]
        [Finite]
        private GivenParameter<double> _subthSlopeL = new GivenParameter<double>();

        [ParameterName("wn0"), ParameterInfo("Width dependence of n0")]
        [Finite]
        private GivenParameter<double> _subthSlopeW = new GivenParameter<double>();

        [ParameterName("nb"), ParameterInfo("VBS dependence of subthreshold slope")]
        [Finite]
        private GivenParameter<double> _subthSlopeB0 = new GivenParameter<double>();

        [ParameterName("lnb"), ParameterInfo("Length dependence of nb")]
        [Finite]
        private GivenParameter<double> _subthSlopeBL = new GivenParameter<double>();

        [ParameterName("wnb"), ParameterInfo("Width dependence of nb")]
        [Finite]
        private GivenParameter<double> _subthSlopeBW = new GivenParameter<double>();

        [ParameterName("nd"), ParameterInfo("VDS dependence of subthreshold slope")]
        [Finite]
        private GivenParameter<double> _subthSlopeD0 = new GivenParameter<double>();

        [ParameterName("lnd"), ParameterInfo("Length dependence of nd")]
        [Finite]
        private GivenParameter<double> _subthSlopeDL = new GivenParameter<double>();

        [ParameterName("wnd"), ParameterInfo("Width dependence of nd")]
        [Finite]
        private GivenParameter<double> _subthSlopeDW = new GivenParameter<double>();

        [ParameterName("tox"), ParameterInfo("Gate oxide thickness in um")]
        [Finite, GreaterThan(0)]
        private GivenParameter<double> _oxideThickness = new GivenParameter<double>(1e-7);

        [ParameterName("temp"), ParameterInfo("Temperature in degree Celcius")]
        [Finite, GreaterThan(-Constants.CelsiusKelvin)]
        private GivenParameter<double> _temp = new GivenParameter<double>(25);

        [ParameterName("vdd"), ParameterInfo("Supply voltage to specify mus")]
        [Finite]
        private GivenParameter<double> _vdd = new GivenParameter<double>(5.0);

        [ParameterName("cgso"), ParameterInfo("Gate source overlap capacitance per unit channel width(m)")]
        [Finite, GreaterThanOrEquals(0)]
        private GivenParameter<double> _gateSourceOverlapCap = new GivenParameter<double>();

        [ParameterName("cgdo"), ParameterInfo("Gate drain overlap capacitance per unit channel width(m)")]
        [Finite, GreaterThanOrEquals(0)]
        private GivenParameter<double> _gateDrainOverlapCap = new GivenParameter<double>();

        [ParameterName("cgbo"), ParameterInfo("Gate bulk overlap capacitance per unit channel length(m)")]
        [Finite, GreaterThanOrEquals(0)]
        private GivenParameter<double> _gateBulkOverlapCap = new GivenParameter<double>();

        [ParameterName("xpart"), ParameterInfo("Flag for channel charge partitioning")]
        [Finite, LowerLimit(0), UpperLimit(1)]
        private GivenParameter<double> _channelChargePartitionFlag = new GivenParameter<double>(1.0);

        [ParameterName("rsh"), ParameterInfo("Source drain diffusion sheet resistance in ohm per square")]
        [Finite, GreaterThanOrEquals(0)]
        private GivenParameter<double> _sheetResistance = new GivenParameter<double>();

        [ParameterName("js"), ParameterInfo("Source drain junction saturation current per unit area")]
        [Finite, GreaterThanOrEquals(0)]
        private GivenParameter<double> _jctSatCurDensity = new GivenParameter<double>();

        [ParameterName("pb"), ParameterInfo("Source drain junction built in potential")]
        [Finite, LowerLimit(0.1)]
        private GivenParameter<double> _bulkJctPotential = new GivenParameter<double>(0.8);

        [ParameterName("mj"), ParameterInfo("Source drain bottom junction capacitance grading coefficient")]
        [Finite]
        private GivenParameter<double> _bulkJctBotGradingCoeff = new GivenParameter<double>(0.5);

        [ParameterName("pbsw"), ParameterInfo("Source drain side junction capacitance built in potential")]
        [Finite, LowerLimit(0.1)]
        private GivenParameter<double> _sidewallJctPotential = new GivenParameter<double>(1.0);

        [ParameterName("mjsw"), ParameterInfo("Source drain side junction capacitance grading coefficient")]
        [Finite]
        private GivenParameter<double> _bulkJctSideGradingCoeff = new GivenParameter<double>(0.33);

        [ParameterName("cj"), ParameterInfo("Source drain bottom junction capacitance per unit area")]
        [Finite, GreaterThanOrEquals(0)]
        private GivenParameter<double> _unitAreaJctCap = new GivenParameter<double>();

        [ParameterName("cjsw"), ParameterInfo("Source drain side junction capacitance per unit area")]
        [Finite, GreaterThanOrEquals(0)]
        private GivenParameter<double> _unitLengthSidewallJctCap = new GivenParameter<double>();

        [ParameterName("wdf"), ParameterInfo("Default width of source drain diffusion in um")]
        [Finite, GreaterThan(0)]
        private GivenParameter<double> _defaultWidth = new GivenParameter<double>();

        [ParameterName("dell"), ParameterInfo("Length reduction of source drain diffusion")]
        [Finite, GreaterThan(0)]
        private GivenParameter<double> _deltaLength = new GivenParameter<double>();

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
        /// Gets the name of the type.
        /// </summary>
        /// <value>
        /// The name of the type.
        /// </value>
        [ParameterName("type"), ParameterInfo("N-channel or P-channel MOS")]
        public string TypeName
        {
            get
            {
                if (Type > 0)
                    return "nmos";
                return "pmos";
            }
        }
    }
}