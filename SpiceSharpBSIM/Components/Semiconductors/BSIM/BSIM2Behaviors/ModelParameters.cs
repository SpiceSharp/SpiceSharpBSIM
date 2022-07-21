using SpiceSharp;
using SpiceSharp.Attributes;
using SpiceSharp.ParameterSets;

namespace SpiceSharpBSIM.Components.Semiconductors.BSIM.BSIM2Behaviors
{
    /// <summary>
    /// Base parameters for a <see cref="BSIM2Model" />
    /// </summary>
    [GeneratedParameters]
    public partial class ModelParameters : ParameterSet<ModelParameters>
    {
        public double Type { get; private set; } = 1;

        [ParameterName("vfb"), ParameterInfo("Flat band voltage")]
        [Finite]
        private GivenParameter<double> _vfb0 = new GivenParameter<double>(-1);

        [ParameterName("lvfb"), ParameterInfo("Length dependence of vfb")]
        [Finite]
        private GivenParameter<double> _vfbL = new GivenParameter<double>();

        [ParameterName("wvfb"), ParameterInfo("Width dependence of vfb")]
        [Finite]
        private GivenParameter<double> _vfbW = new GivenParameter<double>();

        [ParameterName("phi"), ParameterInfo("Strong inversion surface potential")]
        [Finite]
        private GivenParameter<double> _phi0 = new GivenParameter<double>(0.75);

        [ParameterName("lphi"), ParameterInfo("Length dependence of phi")]
        [Finite]
        private GivenParameter<double> _phiL = new GivenParameter<double>();

        [ParameterName("wphi"), ParameterInfo("Width dependence of phi")]
        [Finite]
        private GivenParameter<double> _phiW = new GivenParameter<double>();

        [ParameterName("k1"), ParameterInfo("Bulk effect coefficient 1")]
        [Finite]
        private GivenParameter<double> _k10 = new GivenParameter<double>(0.8);

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

        [ParameterName("eta0"), ParameterInfo("VDS dependence of threshold voltage at VDD=0")]
        [Finite]
        private GivenParameter<double> _eta00 = new GivenParameter<double>();

        [ParameterName("leta0"), ParameterInfo("Length dependence of eta0")]
        [Finite]
        private GivenParameter<double> _eta0L = new GivenParameter<double>();

        [ParameterName("weta0"), ParameterInfo("Width dependence of eta0")]
        [Finite]
        private GivenParameter<double> _eta0W = new GivenParameter<double>();

        [ParameterName("etab"), ParameterInfo("VBS dependence of eta")]
        [Finite]
        private GivenParameter<double> _etaB0 = new GivenParameter<double>();

        [ParameterName("letab"), ParameterInfo("Length dependence of etab")]
        [Finite]
        private GivenParameter<double> _etaBL = new GivenParameter<double>();

        [ParameterName("wetab"), ParameterInfo("Width dependence of etab")]
        [Finite]
        private GivenParameter<double> _etaBW = new GivenParameter<double>();

        [ParameterName("dl"), ParameterInfo("Channel length reduction in um")]
        [Finite]
        private GivenParameter<double> _deltaL = new GivenParameter<double>();

        [ParameterName("dw"), ParameterInfo("Channel width reduction in um")]
        [Finite]
        private GivenParameter<double> _deltaW = new GivenParameter<double>();

        [ParameterName("mu0"), ParameterInfo("Low-field mobility, at VDS=0 VGS=VTH")]
        [Finite]
        private GivenParameter<double> _mob00 = new GivenParameter<double>(400);

        [ParameterName("mu0b"), ParameterInfo("VBS dependence of low-field mobility")]
        [Finite]
        private GivenParameter<double> _mob0B0 = new GivenParameter<double>();

        [ParameterName("lmu0b"), ParameterInfo("Length dependence of mu0b")]
        [Finite]
        private GivenParameter<double> _mob0BL = new GivenParameter<double>();

        [ParameterName("wmu0b"), ParameterInfo("Width dependence of mu0b")]
        [Finite]
        private GivenParameter<double> _mob0BW = new GivenParameter<double>();

        [ParameterName("mus0"), ParameterInfo("Mobility at VDS=VDD VGS=VTH")]
        [Finite, GreaterThan(0)]
        private GivenParameter<double> _mobs00 = new GivenParameter<double>(500);

        [ParameterName("lmus0"), ParameterInfo("Length dependence of mus0")]
        [Finite]
        private GivenParameter<double> _mobs0L = new GivenParameter<double>();

        [ParameterName("wmus0"), ParameterInfo("Width dependence of mus")]
        [Finite]
        private GivenParameter<double> _mobs0W = new GivenParameter<double>();

        [ParameterName("musb"), ParameterInfo("VBS dependence of mus")]
        [Finite]
        private GivenParameter<double> _mobsB0 = new GivenParameter<double>();

        [ParameterName("lmusb"), ParameterInfo("Length dependence of musb")]
        [Finite]
        private GivenParameter<double> _mobsBL = new GivenParameter<double>();

        [ParameterName("wmusb"), ParameterInfo("Width dependence of musb")]
        [Finite]
        private GivenParameter<double> _mobsBW = new GivenParameter<double>();

        [ParameterName("mu20"), ParameterInfo("VDS dependence of mu in tanh term")]
        [Finite]
        private GivenParameter<double> _mob200 = new GivenParameter<double>(1.5);

        [ParameterName("lmu20"), ParameterInfo("Length dependence of mu20")]
        [Finite]
        private GivenParameter<double> _mob20L = new GivenParameter<double>();

        [ParameterName("wmu20"), ParameterInfo("Width dependence of mu20")]
        [Finite]
        private GivenParameter<double> _mob20W = new GivenParameter<double>();

        [ParameterName("mu2b"), ParameterInfo("VBS dependence of mu2")]
        [Finite]
        private GivenParameter<double> _mob2B0 = new GivenParameter<double>();

        [ParameterName("lmu2b"), ParameterInfo("Length dependence of mu2b")]
        [Finite]
        private GivenParameter<double> _mob2BL = new GivenParameter<double>();

        [ParameterName("wmu2b"), ParameterInfo("Width dependence of mu2b")]
        [Finite]
        private GivenParameter<double> _mob2BW = new GivenParameter<double>();

        [ParameterName("mu2g"), ParameterInfo("VGS dependence of mu2")]
        [Finite]
        private GivenParameter<double> _mob2G0 = new GivenParameter<double>();

        [ParameterName("lmu2g"), ParameterInfo("Length dependence of mu2g")]
        [Finite]
        private GivenParameter<double> _mob2GL = new GivenParameter<double>();

        [ParameterName("wmu2g"), ParameterInfo("Width dependence of mu2g")]
        [Finite]
        private GivenParameter<double> _mob2GW = new GivenParameter<double>();

        [ParameterName("mu30"), ParameterInfo("VDS dependence of mu in linear term")]
        [Finite]
        private GivenParameter<double> _mob300 = new GivenParameter<double>(10);

        [ParameterName("lmu30"), ParameterInfo("Length dependence of mu30")]
        [Finite]
        private GivenParameter<double> _mob30L = new GivenParameter<double>();

        [ParameterName("wmu30"), ParameterInfo("Width dependence of mu30")]
        [Finite]
        private GivenParameter<double> _mob30W = new GivenParameter<double>();

        [ParameterName("mu3b"), ParameterInfo("VBS dependence of mu3")]
        [Finite]
        private GivenParameter<double> _mob3B0 = new GivenParameter<double>();

        [ParameterName("lmu3b"), ParameterInfo("Length dependence of mu3b")]
        [Finite]
        private GivenParameter<double> _mob3BL = new GivenParameter<double>();

        [ParameterName("wmu3b"), ParameterInfo("Width dependence of mu3b")]
        [Finite]
        private GivenParameter<double> _mob3BW = new GivenParameter<double>();

        [ParameterName("mu3g"), ParameterInfo("VGS dependence of mu3")]
        [Finite]
        private GivenParameter<double> _mob3G0 = new GivenParameter<double>();

        [ParameterName("lmu3g"), ParameterInfo("Length dependence of mu3g")]
        [Finite]
        private GivenParameter<double> _mob3GL = new GivenParameter<double>();

        [ParameterName("wmu3g"), ParameterInfo("Width dependence of mu3g")]
        [Finite]
        private GivenParameter<double> _mob3GW = new GivenParameter<double>();

        [ParameterName("mu40"), ParameterInfo("VDS dependence of mu in linear term")]
        [Finite]
        private GivenParameter<double> _mob400 = new GivenParameter<double>();

        [ParameterName("lmu40"), ParameterInfo("Length dependence of mu40")]
        [Finite]
        private GivenParameter<double> _mob40L = new GivenParameter<double>();

        [ParameterName("wmu40"), ParameterInfo("Width dependence of mu40")]
        [Finite]
        private GivenParameter<double> _mob40W = new GivenParameter<double>();

        [ParameterName("mu4b"), ParameterInfo("VBS dependence of mu4")]
        [Finite]
        private GivenParameter<double> _mob4B0 = new GivenParameter<double>();

        [ParameterName("lmu4b"), ParameterInfo("Length dependence of mu4b")]
        [Finite]
        private GivenParameter<double> _mob4BL = new GivenParameter<double>();

        [ParameterName("wmu4b"), ParameterInfo("Width dependence of mu4b")]
        [Finite]
        private GivenParameter<double> _mob4BW = new GivenParameter<double>();

        [ParameterName("mu4g"), ParameterInfo("VGS dependence of mu4")]
        [Finite]
        private GivenParameter<double> _mob4G0 = new GivenParameter<double>();

        [ParameterName("lmu4g"), ParameterInfo("Length dependence of mu4g")]
        [Finite]
        private GivenParameter<double> _mob4GL = new GivenParameter<double>();

        [ParameterName("wmu4g"), ParameterInfo("Width dependence of mu4g")]
        [Finite]
        private GivenParameter<double> _mob4GW = new GivenParameter<double>();

        [ParameterName("ua0"), ParameterInfo("Linear VGS dependence of mobility")]
        [Finite]
        private GivenParameter<double> _ua00 = new GivenParameter<double>(0.2);

        [ParameterName("lua0"), ParameterInfo("Length dependence of ua0")]
        [Finite]
        private GivenParameter<double> _ua0L = new GivenParameter<double>();

        [ParameterName("wua0"), ParameterInfo("Width dependence of ua0")]
        [Finite]
        private GivenParameter<double> _ua0W = new GivenParameter<double>();

        [ParameterName("uab"), ParameterInfo("VBS dependence of ua")]
        [Finite]
        private GivenParameter<double> _uaB0 = new GivenParameter<double>();

        [ParameterName("luab"), ParameterInfo("Length dependence of uab")]
        [Finite]
        private GivenParameter<double> _uaBL = new GivenParameter<double>();

        [ParameterName("wuab"), ParameterInfo("Width dependence of uab")]
        [Finite]
        private GivenParameter<double> _uaBW = new GivenParameter<double>();

        [ParameterName("ub0"), ParameterInfo("Quadratic VGS dependence of mobility")]
        [Finite]
        private GivenParameter<double> _ub00 = new GivenParameter<double>();

        [ParameterName("lub0"), ParameterInfo("Length dependence of ub0")]
        [Finite]
        private GivenParameter<double> _ub0L = new GivenParameter<double>();

        [ParameterName("wub0"), ParameterInfo("Width dependence of ub0")]
        [Finite]
        private GivenParameter<double> _ub0W = new GivenParameter<double>();

        [ParameterName("ubb"), ParameterInfo("VBS dependence of ub")]
        [Finite]
        private GivenParameter<double> _ubB0 = new GivenParameter<double>();

        [ParameterName("lubb"), ParameterInfo("Length dependence of ubb")]
        [Finite]
        private GivenParameter<double> _ubBL = new GivenParameter<double>();

        [ParameterName("wubb"), ParameterInfo("Width dependence of ubb")]
        [Finite]
        private GivenParameter<double> _ubBW = new GivenParameter<double>();

        [ParameterName("u10"), ParameterInfo("VDS depence of mobility")]
        [Finite]
        private GivenParameter<double> _u100 = new GivenParameter<double>(0.1);

        [ParameterName("lu10"), ParameterInfo("Length dependence of u10")]
        [Finite]
        private GivenParameter<double> _u10L = new GivenParameter<double>();

        [ParameterName("wu10"), ParameterInfo("Width dependence of u10")]
        [Finite]
        private GivenParameter<double> _u10W = new GivenParameter<double>();

        [ParameterName("u1b"), ParameterInfo("VBS depence of u1")]
        [Finite]
        private GivenParameter<double> _u1B0 = new GivenParameter<double>();

        [ParameterName("lu1b"), ParameterInfo("Length depence of u1b")]
        [Finite]
        private GivenParameter<double> _u1BL = new GivenParameter<double>();

        [ParameterName("wu1b"), ParameterInfo("Width depence of u1b")]
        [Finite]
        private GivenParameter<double> _u1BW = new GivenParameter<double>();

        [ParameterName("u1d"), ParameterInfo("VDS depence of u1")]
        [Finite]
        private GivenParameter<double> _u1D0 = new GivenParameter<double>();

        [ParameterName("lu1d"), ParameterInfo("Length depence of u1d")]
        [Finite]
        private GivenParameter<double> _u1DL = new GivenParameter<double>();

        [ParameterName("wu1d"), ParameterInfo("Width depence of u1d")]
        [Finite]
        private GivenParameter<double> _u1DW = new GivenParameter<double>();

        [ParameterName("n0"), ParameterInfo("Subthreshold slope at VDS=0 VBS=0")]
        [Finite]
        private GivenParameter<double> _n00 = new GivenParameter<double>(1.4);

        [ParameterName("ln0"), ParameterInfo("Length dependence of n0")]
        [Finite]
        private GivenParameter<double> _n0L = new GivenParameter<double>();

        [ParameterName("wn0"), ParameterInfo("Width dependence of n0")]
        [Finite]
        private GivenParameter<double> _n0W = new GivenParameter<double>();

        [ParameterName("nb"), ParameterInfo("VBS dependence of n")]
        [Finite]
        private GivenParameter<double> _nB0 = new GivenParameter<double>(0.5);

        [ParameterName("lnb"), ParameterInfo("Length dependence of nb")]
        [Finite]
        private GivenParameter<double> _nBL = new GivenParameter<double>();

        [ParameterName("wnb"), ParameterInfo("Width dependence of nb")]
        [Finite]
        private GivenParameter<double> _nBW = new GivenParameter<double>();

        [ParameterName("nd"), ParameterInfo("VDS dependence of n")]
        [Finite]
        private GivenParameter<double> _nD0 = new GivenParameter<double>();

        [ParameterName("lnd"), ParameterInfo("Length dependence of nd")]
        [Finite]
        private GivenParameter<double> _nDL = new GivenParameter<double>();

        [ParameterName("wnd"), ParameterInfo("Width dependence of nd")]
        [Finite]
        private GivenParameter<double> _nDW = new GivenParameter<double>();

        [ParameterName("vof0"), ParameterInfo("Threshold voltage offset AT VDS=0 VBS=0")]
        [Finite]
        private GivenParameter<double> _vof00 = new GivenParameter<double>(1.8);

        [ParameterName("lvof0"), ParameterInfo("Length dependence of vof0")]
        [Finite]
        private GivenParameter<double> _vof0L = new GivenParameter<double>();

        [ParameterName("wvof0"), ParameterInfo("Width dependence of vof0")]
        [Finite]
        private GivenParameter<double> _vof0W = new GivenParameter<double>();

        [ParameterName("vofb"), ParameterInfo("VBS dependence of vof")]
        [Finite]
        private GivenParameter<double> _vofB0 = new GivenParameter<double>();

        [ParameterName("lvofb"), ParameterInfo("Length dependence of vofb")]
        [Finite]
        private GivenParameter<double> _vofBL = new GivenParameter<double>();

        [ParameterName("wvofb"), ParameterInfo("Width dependence of vofb")]
        [Finite]
        private GivenParameter<double> _vofBW = new GivenParameter<double>();

        [ParameterName("vofd"), ParameterInfo("VDS dependence of vof")]
        [Finite]
        private GivenParameter<double> _vofD0 = new GivenParameter<double>();

        [ParameterName("lvofd"), ParameterInfo("Length dependence of vofd")]
        [Finite]
        private GivenParameter<double> _vofDL = new GivenParameter<double>();

        [ParameterName("wvofd"), ParameterInfo("Width dependence of vofd")]
        [Finite]
        private GivenParameter<double> _vofDW = new GivenParameter<double>();

        [ParameterName("ai0"), ParameterInfo("Pre-factor of hot-electron effect.")]
        [Finite]
        private GivenParameter<double> _ai00 = new GivenParameter<double>();

        [ParameterName("lai0"), ParameterInfo("Length dependence of ai0")]
        [Finite]
        private GivenParameter<double> _ai0L = new GivenParameter<double>();

        [ParameterName("wai0"), ParameterInfo("Width dependence of ai0")]
        [Finite]
        private GivenParameter<double> _ai0W = new GivenParameter<double>();

        [ParameterName("aib"), ParameterInfo("VBS dependence of ai")]
        [Finite]
        private GivenParameter<double> _aiB0 = new GivenParameter<double>();

        [ParameterName("laib"), ParameterInfo("Length dependence of aib")]
        [Finite]
        private GivenParameter<double> _aiBL = new GivenParameter<double>();

        [ParameterName("waib"), ParameterInfo("Width dependence of aib")]
        [Finite]
        private GivenParameter<double> _aiBW = new GivenParameter<double>();

        [ParameterName("bi0"), ParameterInfo("Exponential factor of hot-electron effect.")]
        [Finite]
        private GivenParameter<double> _bi00 = new GivenParameter<double>();

        [ParameterName("lbi0"), ParameterInfo("Length dependence of bi0")]
        [Finite]
        private GivenParameter<double> _bi0L = new GivenParameter<double>();

        [ParameterName("wbi0"), ParameterInfo("Width dependence of bi0")]
        [Finite]
        private GivenParameter<double> _bi0W = new GivenParameter<double>();

        [ParameterName("bib"), ParameterInfo("VBS dependence of bi")]
        [Finite]
        private GivenParameter<double> _biB0 = new GivenParameter<double>();

        [ParameterName("lbib"), ParameterInfo("Length dependence of bib")]
        [Finite]
        private GivenParameter<double> _biBL = new GivenParameter<double>();

        [ParameterName("wbib"), ParameterInfo("Width dependence of bib")]
        [Finite]
        private GivenParameter<double> _biBW = new GivenParameter<double>();

        [ParameterName("vghigh"), ParameterInfo("Upper bound of the cubic spline function.")]
        [Finite]
        private GivenParameter<double> _vghigh0 = new GivenParameter<double>(0.2);

        [ParameterName("lvghigh"), ParameterInfo("Length dependence of vghigh")]
        [Finite]
        private GivenParameter<double> _vghighL = new GivenParameter<double>();

        [ParameterName("wvghigh"), ParameterInfo("Width dependence of vghigh")]
        [Finite]
        private GivenParameter<double> _vghighW = new GivenParameter<double>();

        [ParameterName("vglow"), ParameterInfo("Lower bound of the cubic spline function.")]
        [Finite]
        private GivenParameter<double> _vglow0 = new GivenParameter<double>(-0.15);

        [ParameterName("lvglow"), ParameterInfo("Length dependence of vglow")]
        [Finite]
        private GivenParameter<double> _vglowL = new GivenParameter<double>();

        [ParameterName("wvglow"), ParameterInfo("Width dependence of vglow")]
        [Finite]
        private GivenParameter<double> _vglowW = new GivenParameter<double>();

        [ParameterName("tox"), ParameterInfo("Gate oxide thickness in um")]
        [Finite, GreaterThan(0)]
        private GivenParameter<double> _tox = new GivenParameter<double>(0.03);

        [ParameterName("temp"), ParameterInfo("Temperature in degree Celcius")]
        [Finite, GreaterThan(-Constants.CelsiusKelvin)]
        private GivenParameter<double> _temp = new GivenParameter<double>(27);

        [ParameterName("vdd"), ParameterInfo("Maximum Vds")]
        [Finite]
        private GivenParameter<double> _vdd = new GivenParameter<double>(5);

        [ParameterName("vgg"), ParameterInfo("Maximum Vgs")]
        [Finite]
        private GivenParameter<double> _vgg = new GivenParameter<double>(5);

        [ParameterName("vbb"), ParameterInfo("Maximum Vbs")]
        [Finite]
        private GivenParameter<double> _vbb = new GivenParameter<double>(5);

        [ParameterName("cgso"), ParameterInfo("Gate source overlap capacitance per unit channel width (m)")]
        [Finite]
        private GivenParameter<double> _gateSourceOverlapCap = new GivenParameter<double>();

        [ParameterName("cgdo"), ParameterInfo("Gate drain overlap capacitance per unit channel width (m)")]
        [Finite]
        private GivenParameter<double> _gateDrainOverlapCap = new GivenParameter<double>();

        [ParameterName("cgbo"), ParameterInfo("Gate bulk overlap capacitance per unit channel length (m)")]
        [Finite]
        private GivenParameter<double> _gateBulkOverlapCap = new GivenParameter<double>();

        [ParameterName("xpart"), ParameterInfo("Flag for channel charge partitioning")]
        [Finite]
        private GivenParameter<double> _channelChargePartitionFlag = new GivenParameter<double>();

        [ParameterName("rsh"), ParameterInfo("Source drain diffusion sheet resistance in ohm per square")]
        [Finite, GreaterThanOrEquals(0)]
        private GivenParameter<double> _sheetResistance = new GivenParameter<double>();

        [ParameterName("js"), ParameterInfo("Source drain junction saturation current per unit area")]
        [Finite, GreaterThanOrEquals(0)]
        private GivenParameter<double> _jctSatCurDensity = new GivenParameter<double>();

        [ParameterName("pb"), ParameterInfo("Source drain junction built in potential")]
        [Finite, LowerLimit(0.1)]
        private GivenParameter<double> _bulkJctPotential = new GivenParameter<double>();

        [ParameterName("mj"), ParameterInfo("Source drain bottom junction capacitance grading coefficient")]
        [Finite]
        private GivenParameter<double> _bulkJctBotGradingCoeff = new GivenParameter<double>();

        [ParameterName("pbsw"), ParameterInfo("Source drain side junction capacitance built in potential")]
        [Finite, LowerLimit(0.1)]
        private GivenParameter<double> _sidewallJctPotential = new GivenParameter<double>();

        [ParameterName("mjsw"), ParameterInfo("Source drain side junction capacitance grading coefficient")]
        [Finite]
        private GivenParameter<double> _bulkJctSideGradingCoeff = new GivenParameter<double>();

        [ParameterName("cj"), ParameterInfo("Source drain bottom junction capacitance per unit area")]
        [Finite, GreaterThanOrEquals(0)]
        private GivenParameter<double> _unitAreaJctCap = new GivenParameter<double>();

        [ParameterName("cjsw"), ParameterInfo("Source drain side junction capacitance per unit area")]
        [Finite, GreaterThanOrEquals(0)]
        private GivenParameter<double> _unitLengthSidewallJctCap = new GivenParameter<double>();

        [ParameterName("wdf"), ParameterInfo("Default width of source drain diffusion in um")]
        [Finite, GreaterThan(0)]
        private GivenParameter<double> _defaultWidth = new GivenParameter<double>(10);

        [ParameterName("dell"), ParameterInfo("Length reduction of source drain diffusion")]
        [Finite]
        private GivenParameter<double> _deltaLength = new GivenParameter<double>();

        [ParameterName("kf"), ParameterInfo("Flicker noise coefficient")]
        [Finite]
        private GivenParameter<double> _fNcoef = new GivenParameter<double>();

        [ParameterName("af"), ParameterInfo("Flicker noise exponent")]
        [Finite]
        private GivenParameter<double> _fNexp = new GivenParameter<double>(1.0);

        [ParameterName("nmos"), ParameterInfo("Flag to indicate NMOS")]
        public void SetNmos(bool flag = true)
        {
            if (flag)
                Type = 1;
        }

        [ParameterName("pmos"), ParameterInfo("Flag to indicate PMOS")]
        public void SetPmos(bool flag = true)
        {
            if (flag)
                Type = -1;
        }
    }
}