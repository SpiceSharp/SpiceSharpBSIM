﻿using SpiceSharp;
using SpiceSharp.Attributes;
using SpiceSharp.ParameterSets;

namespace SpiceSharpBSIM.Components.Semiconductors.BSIM.BSIM3v2Behaviors
{
    /// <summary>
    /// Model parameters for a <see cref="BSIM3v2Model"/>.
    /// </summary>
    [GeneratedParameters]
    public partial class ModelParameters : ParameterSet<ModelParameters>
    {
        [ParameterName("path"), ParameterInfo("Filename for used to log parameter checks")]
        public string CheckPath { get; set; } = "b3v2check.log";

        [ParameterName("capmod"), ParameterInfo("Capacitance model selector")]
        private GivenParameter<int> _capMod = new GivenParameter<int>(3);

        [ParameterName("mobmod"), ParameterInfo("Mobility model selector")]
        private GivenParameter<int> _mobMod = new GivenParameter<int>(1);

        [ParameterName("noimod"), ParameterInfo("Noise model selector")]
        private GivenParameter<int> _noiMod = new GivenParameter<int>(1);

        [ParameterName("nqsmod"), ParameterInfo("Non-quasi-static model selector")]
        private GivenParameter<int> _nqsMod = new GivenParameter<int>();

        [ParameterName("acm"), ParameterInfo("Area calculation method selector")]
        private GivenParameter<int> _acmMod = new GivenParameter<int>(0);

        [ParameterName("calcacm"), ParameterInfo("Area calculation method ACM=12")]
        private GivenParameter<int> _calcacm = new GivenParameter<int>(0);

        [ParameterName("paramchk"), ParameterInfo("Model parameter checking selector")]
        private GivenParameter<int> _paramChk = new GivenParameter<int>(0);

        [ParameterName("binunit"), ParameterInfo("Bin  unit  selector")]
        private GivenParameter<int> _binUnit = new GivenParameter<int>(1);

        [ParameterName("version"), ParameterInfo(" parameter for model version")]
        private GivenParameter<string> _version = new GivenParameter<string>("3.2.4");

        [ParameterName("tox"), ParameterInfo("Gate oxide thickness in meters")]
        [Finite]
        private GivenParameter<double> _tox = new GivenParameter<double>(150.0e-10);

        [ParameterName("toxm"), ParameterInfo("Gate oxide thickness used in extraction")]
        [Finite]
        private GivenParameter<double> _toxm = new GivenParameter<double>();

        [ParameterName("cdsc"), ParameterInfo("Drain/Source and channel coupling capacitance")]
        [Finite]
        private GivenParameter<double> _cdsc = new GivenParameter<double>(2.4e-4);

        [ParameterName("cdscb"), ParameterInfo("Body-bias dependence of cdsc")]
        [Finite]
        private GivenParameter<double> _cdscb = new GivenParameter<double>(0.0);

        [ParameterName("cdscd"), ParameterInfo("Drain-bias dependence of cdsc")]
        [Finite]
        private GivenParameter<double> _cdscd = new GivenParameter<double>(0.0);

        [ParameterName("cit"), ParameterInfo("Interface state capacitance")]
        [Finite]
        private GivenParameter<double> _cit = new GivenParameter<double>(0.0);

        [ParameterName("nfactor"), ParameterInfo("Subthreshold swing Coefficient")]
        [Finite]
        private GivenParameter<double> _nfactor = new GivenParameter<double>(1);

        [ParameterName("xj"), ParameterInfo("Junction depth in meters")]
        [Finite]
        private GivenParameter<double> _xj = new GivenParameter<double>(.15e-6);

        [ParameterName("vsat"), ParameterInfo("Saturation velocity at tnom")]
        [Finite]
        private GivenParameter<double> _vsat = new GivenParameter<double>(8.0e4);

        [ParameterName("at"), ParameterInfo("Temperature coefficient of vsat")]
        [Finite]
        private GivenParameter<double> _at = new GivenParameter<double>(3.3e4);

        [ParameterName("a0"), ParameterInfo("Non-uniform depletion width effect coefficient.")]
        [Finite]
        private GivenParameter<double> _a0 = new GivenParameter<double>(1.0);

        [ParameterName("ags"), ParameterInfo("Gate bias  coefficient of Abulk.")]
        [Finite]
        private GivenParameter<double> _ags = new GivenParameter<double>(0.0);

        [ParameterName("a1"), ParameterInfo("Non-saturation effect coefficient")]
        [Finite]
        private GivenParameter<double> _a1 = new GivenParameter<double>(0.0);

        [ParameterName("a2"), ParameterInfo("Non-saturation effect coefficient")]
        [Finite]
        private GivenParameter<double> _a2 = new GivenParameter<double>(1.0);

        [ParameterName("keta"), ParameterInfo("Body-bias coefficient of non-uniform depletion width effect.")]
        [Finite]
        private GivenParameter<double> _keta = new GivenParameter<double>(-0.047);

        [ParameterName("nsub"), ParameterInfo("Substrate doping concentration")]
        [Finite]
        private GivenParameter<double> _nsub = new GivenParameter<double>(6.0e16);

        [ParameterName("nch"), ParameterInfo("Channel doping concentration", Units = "1/cm^3")]
        [Finite]
        public GivenParameter<double> Npeak
        {
            get => _npeak;
            set
            {
                Utility.Finite(value, nameof(Npeak));
                if (value > 1e20)
                    _npeak = value * 1e-6;
                else
                    _npeak = value;
            }
        }
        private GivenParameter<double> _npeak = new GivenParameter<double>(1.7e+17);

        [ParameterName("ngate"), ParameterInfo("Poly-gate doping concentration", Units = "1/cm^3")]
        [Finite]
        public GivenParameter<double> Ngate
        {
            get => _ngate;
            set
            {
                Utility.Finite(value, nameof(Ngate));
                if (value > 1.000001e24)
                    _ngate = value * 1e-6;
                else
                    _ngate = value;
            }
        }
        private GivenParameter<double> _ngate = new GivenParameter<double>();

        [ParameterName("gamma1"), ParameterInfo("Vth body coefficient")]
        [Finite]
        private GivenParameter<double> _gamma1 = new GivenParameter<double>();

        [ParameterName("gamma2"), ParameterInfo("Vth body coefficient")]
        [Finite]
        private GivenParameter<double> _gamma2 = new GivenParameter<double>();

        [ParameterName("vbx"), ParameterInfo("Vth transition body Voltage")]
        [Finite]
        private GivenParameter<double> _vbx = new GivenParameter<double>();

        [ParameterName("vbm"), ParameterInfo("Maximum body voltage")]
        [Finite]
        private GivenParameter<double> _vbm = new GivenParameter<double>(-3.0);

        [ParameterName("xt"), ParameterInfo("Doping depth")]
        [Finite]
        private GivenParameter<double> _xt = new GivenParameter<double>(1.55e-7);

        [ParameterName("k1"), ParameterInfo("Bulk effect coefficient 1")]
        [Finite]
        private GivenParameter<double> _k1 = new GivenParameter<double>();

        [ParameterName("kt1"), ParameterInfo("Temperature coefficient of Vth")]
        [Finite]
        private GivenParameter<double> _kt1 = new GivenParameter<double>(-0.11);

        [ParameterName("kt1l"), ParameterInfo("Temperature coefficient of Vth")]
        [Finite]
        private GivenParameter<double> _kt1l = new GivenParameter<double>(0.0);

        [ParameterName("kt2"), ParameterInfo("Body-coefficient of kt1")]
        [Finite]
        private GivenParameter<double> _kt2 = new GivenParameter<double>(0.022);

        [ParameterName("k2"), ParameterInfo("Bulk effect coefficient 2")]
        [Finite]
        private GivenParameter<double> _k2 = new GivenParameter<double>();

        [ParameterName("k3"), ParameterInfo("Narrow width effect coefficient")]
        [Finite]
        private GivenParameter<double> _k3 = new GivenParameter<double>(80.0);

        [ParameterName("k3b"), ParameterInfo("Body effect coefficient of k3")]
        [Finite]
        private GivenParameter<double> _k3b = new GivenParameter<double>(0.0);

        [ParameterName("w0"), ParameterInfo("Narrow width effect parameter")]
        [Finite]
        private GivenParameter<double> _w0 = new GivenParameter<double>(2.5e-6);

        [ParameterName("nlx"), ParameterInfo("Lateral non-uniform doping effect")]
        [Finite]
        private GivenParameter<double> _nlx = new GivenParameter<double>(1.74e-7);

        [ParameterName("dvt0"), ParameterInfo("Short channel effect coeff. 0")]
        [Finite]
        private GivenParameter<double> _dvt0 = new GivenParameter<double>(2.2);

        [ParameterName("dvt1"), ParameterInfo("Short channel effect coeff. 1")]
        [Finite]
        private GivenParameter<double> _dvt1 = new GivenParameter<double>(0.53);

        [ParameterName("dvt2"), ParameterInfo("Short channel effect coeff. 2")]
        [Finite]
        private GivenParameter<double> _dvt2 = new GivenParameter<double>(-0.032);

        [ParameterName("dvt0w"), ParameterInfo("Narrow Width coeff. 0")]
        [Finite]
        private GivenParameter<double> _dvt0w = new GivenParameter<double>(0.0);

        [ParameterName("dvt1w"), ParameterInfo("Narrow Width effect coeff. 1")]
        [Finite]
        private GivenParameter<double> _dvt1w = new GivenParameter<double>(5.3e6);

        [ParameterName("dvt2w"), ParameterInfo("Narrow Width effect coeff. 2")]
        [Finite]
        private GivenParameter<double> _dvt2w = new GivenParameter<double>(-0.032);

        [ParameterName("drout"), ParameterInfo("DIBL coefficient of output resistance")]
        [Finite]
        private GivenParameter<double> _drout = new GivenParameter<double>(0.56);

        [ParameterName("dsub"), ParameterInfo("DIBL coefficient in the subthreshold region")]
        [Finite]
        private GivenParameter<double> _dsub = new GivenParameter<double>();

        [ParameterName("vtho"), ParameterName("vth0"), ParameterInfo("Threshold voltage")]
        [Finite]
        private GivenParameter<double> _vth0 = new GivenParameter<double>();

        [ParameterName("ua"), ParameterInfo("Linear gate dependence of mobility")]
        [Finite]
        private GivenParameter<double> _ua = new GivenParameter<double>(2.25e-9);

        [ParameterName("ua1"), ParameterInfo("Temperature coefficient of ua")]
        [Finite]
        private GivenParameter<double> _ua1 = new GivenParameter<double>(4.31e-9);

        [ParameterName("ub"), ParameterInfo("Quadratic gate dependence of mobility")]
        [Finite]
        private GivenParameter<double> _ub = new GivenParameter<double>(5.87e-19);

        [ParameterName("ub1"), ParameterInfo("Temperature coefficient of ub")]
        [Finite]
        private GivenParameter<double> _ub1 = new GivenParameter<double>(-7.61e-18);

        [ParameterName("uc"), ParameterInfo("Body-bias dependence of mobility")]
        [Finite]
        private GivenParameter<double> _uc = new GivenParameter<double>();

        [ParameterName("uc1"), ParameterInfo("Temperature coefficient of uc")]
        [Finite]
        private GivenParameter<double> _uc1 = new GivenParameter<double>();

        [ParameterName("u0"), ParameterInfo("Low-field mobility at Tnom")]
        [Finite]
        private GivenParameter<double> _u0 = new GivenParameter<double>();

        [ParameterName("ute"), ParameterInfo("Temperature coefficient of mobility")]
        [Finite]
        private GivenParameter<double> _ute = new GivenParameter<double>(-1.5);

        [ParameterName("voff"), ParameterInfo("Threshold voltage offset")]
        [Finite]
        private GivenParameter<double> _voff = new GivenParameter<double>(-0.08);

        [ParameterName("cgso"), ParameterInfo("Gate-source overlap capacitance per width")]
        [Finite]
        private GivenParameter<double> _cgso = new GivenParameter<double>();

        [ParameterName("cgdo"), ParameterInfo("Gate-drain overlap capacitance per width")]
        [Finite]
        private GivenParameter<double> _cgdo = new GivenParameter<double>();

        [ParameterName("cgbo"), ParameterInfo("Gate-bulk overlap capacitance per length")]
        [Finite]
        private GivenParameter<double> _cgbo = new GivenParameter<double>();

        [ParameterName("xpart"), ParameterInfo("Channel charge partitioning")]
        [Finite]
        private GivenParameter<double> _xpart = new GivenParameter<double>(0.0);

        [ParameterName("elm"), ParameterInfo("Non-quasi-static Elmore Constant Parameter")]
        [Finite]
        private GivenParameter<double> _elm = new GivenParameter<double>(5.0);

        [ParameterName("delta"), ParameterInfo("Effective Vds parameter")]
        [Finite]
        private GivenParameter<double> _delta = new GivenParameter<double>(0.01);

        [ParameterName("rsh"), ParameterInfo("Source-drain sheet resistance")]
        [Finite]
        private GivenParameter<double> _sheetResistance = new GivenParameter<double>(0.0);

        [ParameterName("rdsw"), ParameterInfo("Source-drain resistance per width")]
        [Finite]
        private GivenParameter<double> _rdsw = new GivenParameter<double>(0);

        [ParameterName("prwg"), ParameterInfo("Gate-bias effect on parasitic resistance ")]
        [Finite]
        private GivenParameter<double> _prwg = new GivenParameter<double>(0.0);

        [ParameterName("prwb"), ParameterInfo("Body-effect on parasitic resistance ")]
        [Finite]
        private GivenParameter<double> _prwb = new GivenParameter<double>(0.0);

        [ParameterName("prt"), ParameterInfo("Temperature coefficient of parasitic resistance ")]
        [Finite]
        private GivenParameter<double> _prt = new GivenParameter<double>(0.0);

        [ParameterName("eta0"), ParameterInfo("Subthreshold region DIBL coefficient")]
        [Finite]
        private GivenParameter<double> _eta0 = new GivenParameter<double>(0.08);

        [ParameterName("etab"), ParameterInfo("Subthreshold region DIBL coefficient")]
        [Finite]
        private GivenParameter<double> _etab = new GivenParameter<double>(-0.07);

        [ParameterName("pclm"), ParameterInfo("Channel length modulation Coefficient")]
        [Finite]
        private GivenParameter<double> _pclm = new GivenParameter<double>(1.3);

        [ParameterName("pdiblc1"), ParameterInfo("Drain-induced barrier lowering coefficient")]
        [Finite]
        private GivenParameter<double> _pdibl1 = new GivenParameter<double>(.39);

        [ParameterName("pdiblc2"), ParameterInfo("Drain-induced barrier lowering coefficient")]
        [Finite]
        private GivenParameter<double> _pdibl2 = new GivenParameter<double>(0.0086);

        [ParameterName("pdiblcb"), ParameterInfo("Body-effect on drain-induced barrier lowering")]
        [Finite]
        private GivenParameter<double> _pdiblb = new GivenParameter<double>(0.0);

        [ParameterName("pscbe1"), ParameterInfo("Substrate current body-effect coefficient")]
        [Finite]
        private GivenParameter<double> _pscbe1 = new GivenParameter<double>(4.24e8);

        [ParameterName("pscbe2"), ParameterInfo("Substrate current body-effect coefficient")]
        [Finite]
        private GivenParameter<double> _pscbe2 = new GivenParameter<double>(1.0e-5);

        [ParameterName("pvag"), ParameterInfo("Gate dependence of output resistance parameter")]
        [Finite]
        private GivenParameter<double> _pvag = new GivenParameter<double>(0.0);

        [ParameterName("js"), ParameterInfo("Source/drain junction reverse saturation current density")]
        [Finite]
        private GivenParameter<double> _jctSatCurDensity = new GivenParameter<double>(1.0E-4);

        [ParameterName("jsw"), ParameterInfo("Sidewall junction reverse saturation current density")]
        [Finite]
        private GivenParameter<double> _jctSidewallSatCurDensity = new GivenParameter<double>(0.0);

        [ParameterName("pb"), ParameterInfo("Source/drain junction built-in potential")]
        [Finite, LowerLimit(0.1)]
        private GivenParameter<double> _bulkJctPotential = new GivenParameter<double>(1.0);

        [ParameterName("nj"), ParameterInfo("Source/drain junction emission coefficient")]
        [Finite]
        private GivenParameter<double> _jctEmissionCoeff = new GivenParameter<double>(1.0);

        [ParameterName("xti"), ParameterInfo("Junction current temperature exponent")]
        [Finite]
        private GivenParameter<double> _jctTempExponent = new GivenParameter<double>(3.0);

        [ParameterName("mj"), ParameterInfo("Source/drain bottom junction capacitance grading coefficient")]
        [Finite]
        private GivenParameter<double> _bulkJctBotGradingCoeff = new GivenParameter<double>(0.5);

        [ParameterName("pbsw"), ParameterInfo("Source/drain sidewall junction capacitance built in potential")]
        [Finite, LowerLimit(0.1)]
        private GivenParameter<double> _sidewallJctPotential = new GivenParameter<double>(1.0);

        [ParameterName("mjsw"), ParameterInfo("Source/drain sidewall junction capacitance grading coefficient")]
        [Finite]
        private GivenParameter<double> _bulkJctSideGradingCoeff = new GivenParameter<double>(0.33);

        [ParameterName("pbswg"), ParameterInfo("Source/drain (gate side) sidewall junction capacitance built in potential")]
        [Finite, LowerLimit(0.1)]
        private GivenParameter<double> _gatesidewallJctPotential = new GivenParameter<double>();

        [ParameterName("mjswg"), ParameterInfo("Source/drain (gate side) sidewall junction capacitance grading coefficient")]
        [Finite]
        private GivenParameter<double> _bulkJctGateSideGradingCoeff = new GivenParameter<double>();

        [ParameterName("cj"), ParameterInfo("Source/drain bottom junction capacitance per unit area")]
        [Finite]
        private GivenParameter<double> _unitAreaJctCap = new GivenParameter<double>(5.0E-4);

        [ParameterName("vfbcv"), ParameterInfo("Flat Band Voltage parameter for capmod=0 only")]
        [Finite]
        private GivenParameter<double> _vfbcv = new GivenParameter<double>(-1.0);

        [ParameterName("vfb"), ParameterInfo("Flat Band Voltage")]
        [Finite]
        private GivenParameter<double> _vfb = new GivenParameter<double>();

        [ParameterName("cjsw"), ParameterInfo("Source/drain sidewall junction capacitance per unit periphery")]
        [Finite]
        private GivenParameter<double> _unitLengthSidewallJctCap = new GivenParameter<double>(5.0E-10);

        [ParameterName("cjswg"), ParameterInfo("Source/drain (gate side) sidewall junction capacitance per unit width")]
        [Finite]
        private GivenParameter<double> _unitLengthGateSidewallJctCap = new GivenParameter<double>();

        [ParameterName("tpb"), ParameterInfo("Temperature coefficient of pb")]
        [Finite]
        private GivenParameter<double> _tpb = new GivenParameter<double>(0.0);

        [ParameterName("tcj"), ParameterInfo("Temperature coefficient of cj")]
        [Finite]
        private GivenParameter<double> _tcj = new GivenParameter<double>(0.0);

        [ParameterName("tpbsw"), ParameterInfo("Temperature coefficient of pbsw")]
        [Finite]
        private GivenParameter<double> _tpbsw = new GivenParameter<double>(0.0);

        [ParameterName("tcjsw"), ParameterInfo("Temperature coefficient of cjsw")]
        [Finite]
        private GivenParameter<double> _tcjsw = new GivenParameter<double>(0.0);

        [ParameterName("tpbswg"), ParameterInfo("Temperature coefficient of pbswg")]
        [Finite]
        private GivenParameter<double> _tpbswg = new GivenParameter<double>(0.0);

        [ParameterName("tcjswg"), ParameterInfo("Temperature coefficient of cjswg")]
        [Finite]
        private GivenParameter<double> _tcjswg = new GivenParameter<double>(0.0);

        [ParameterName("acde"), ParameterInfo("Exponential coefficient for finite charge thickness")]
        [Finite]
        private GivenParameter<double> _acde = new GivenParameter<double>(1.0);

        [ParameterName("moin"), ParameterInfo("Coefficient for gate-bias dependent surface potential")]
        [Finite]
        private GivenParameter<double> _moin = new GivenParameter<double>(15.0);

        [ParameterName("noff"), ParameterInfo("C-V turn-on/off parameter")]
        [Finite]
        private GivenParameter<double> _noff = new GivenParameter<double>(1.0);

        [ParameterName("voffcv"), ParameterInfo("C-V lateral-shift parameter")]
        [Finite]
        private GivenParameter<double> _voffcv = new GivenParameter<double>(0.0);

        [ParameterName("lint"), ParameterInfo("Length reduction parameter")]
        [Finite]
        private GivenParameter<double> _lint = new GivenParameter<double>(0.0);

        [ParameterName("ll"), ParameterInfo("Length reduction parameter")]
        [Finite]
        private GivenParameter<double> _ll = new GivenParameter<double>(0.0);

        [ParameterName("llc"), ParameterInfo("Length reduction parameter for CV")]
        [Finite]
        private GivenParameter<double> _llc = new GivenParameter<double>();

        [ParameterName("lln"), ParameterInfo("Length reduction parameter")]
        [Finite]
        private GivenParameter<double> _lln = new GivenParameter<double>(1.0);

        [ParameterName("lw"), ParameterInfo("Length reduction parameter")]
        [Finite]
        private GivenParameter<double> _lw = new GivenParameter<double>(0.0);

        [ParameterName("lwc"), ParameterInfo("Length reduction parameter for CV")]
        [Finite]
        private GivenParameter<double> _lwc = new GivenParameter<double>();

        [ParameterName("lwn"), ParameterInfo("Length reduction parameter")]
        [Finite]
        private GivenParameter<double> _lwn = new GivenParameter<double>(1.0);

        [ParameterName("lwl"), ParameterInfo("Length reduction parameter")]
        [Finite]
        private GivenParameter<double> _lwl = new GivenParameter<double>(0.0);

        [ParameterName("lwlc"), ParameterInfo("Length reduction parameter for CV")]
        [Finite]
        private GivenParameter<double> _lwlc = new GivenParameter<double>();

        [ParameterName("lmin"), ParameterInfo("Minimum length for the model")]
        [Finite]
        private GivenParameter<double> _lmin = new GivenParameter<double>(0.0);

        [ParameterName("lmax"), ParameterInfo("Maximum length for the model")]
        [Finite]
        private GivenParameter<double> _lmax = new GivenParameter<double>(1.0);

        [ParameterName("xl"), ParameterInfo("Length correction parameter")]
        [Finite]
        private GivenParameter<double> _xl = new GivenParameter<double>(0.0);

        [ParameterName("xw"), ParameterInfo("Width correction parameter")]
        [Finite]
        private GivenParameter<double> _xw = new GivenParameter<double>(0.0);

        [ParameterName("wr"), ParameterInfo("Width dependence of rds")]
        [Finite]
        private GivenParameter<double> _wr = new GivenParameter<double>(1.0);

        [ParameterName("wint"), ParameterInfo("Width reduction parameter")]
        [Finite]
        private GivenParameter<double> _wint = new GivenParameter<double>(0.0);

        [ParameterName("dwg"), ParameterInfo("Width reduction parameter")]
        [Finite]
        private GivenParameter<double> _dwg = new GivenParameter<double>(0.0);

        [ParameterName("dwb"), ParameterInfo("Width reduction parameter")]
        [Finite]
        private GivenParameter<double> _dwb = new GivenParameter<double>(0.0);

        [ParameterName("wl"), ParameterInfo("Width reduction parameter")]
        [Finite]
        private GivenParameter<double> _wl = new GivenParameter<double>(0.0);

        [ParameterName("wlc"), ParameterInfo("Width reduction parameter for CV")]
        [Finite]
        private GivenParameter<double> _wlc = new GivenParameter<double>();

        [ParameterName("wln"), ParameterInfo("Width reduction parameter")]
        [Finite]
        private GivenParameter<double> _wln = new GivenParameter<double>(1.0);

        [ParameterName("ww"), ParameterInfo("Width reduction parameter")]
        [Finite]
        private GivenParameter<double> _ww = new GivenParameter<double>(0.0);

        [ParameterName("wwc"), ParameterInfo("Width reduction parameter for CV")]
        [Finite]
        private GivenParameter<double> _wwc = new GivenParameter<double>();

        [ParameterName("wwn"), ParameterInfo("Width reduction parameter")]
        [Finite]
        private GivenParameter<double> _wwn = new GivenParameter<double>(1.0);

        [ParameterName("wwl"), ParameterInfo("Width reduction parameter")]
        [Finite]
        private GivenParameter<double> _wwl = new GivenParameter<double>(0.0);

        [ParameterName("wwlc"), ParameterInfo("Width reduction parameter for CV")]
        [Finite]
        private GivenParameter<double> _wwlc = new GivenParameter<double>();

        [ParameterName("wmin"), ParameterInfo("Minimum width for the model")]
        [Finite]
        private GivenParameter<double> _wmin = new GivenParameter<double>(0.0);

        [ParameterName("wmax"), ParameterInfo("Maximum width for the model")]
        [Finite]
        private GivenParameter<double> _wmax = new GivenParameter<double>(1.0);

        [ParameterName("b0"), ParameterInfo("Abulk narrow width parameter")]
        [Finite]
        private GivenParameter<double> _b0 = new GivenParameter<double>(0.0);

        [ParameterName("b1"), ParameterInfo("Abulk narrow width parameter")]
        [Finite]
        private GivenParameter<double> _b1 = new GivenParameter<double>(0.0);

        [ParameterName("cgsl"), ParameterInfo("New C-V model parameter")]
        [Finite]
        private GivenParameter<double> _cgsl = new GivenParameter<double>(0.0);

        [ParameterName("cgdl"), ParameterInfo("New C-V model parameter")]
        [Finite]
        private GivenParameter<double> _cgdl = new GivenParameter<double>(0.0);

        [ParameterName("ckappa"), ParameterInfo("New C-V model parameter")]
        [Finite]
        private GivenParameter<double> _ckappa = new GivenParameter<double>(0.6);

        [ParameterName("cf"), ParameterInfo("Fringe capacitance parameter")]
        [Finite]
        private GivenParameter<double> _cf = new GivenParameter<double>();

        [ParameterName("clc"), ParameterInfo("Vdsat parameter for C-V model")]
        [Finite]
        private GivenParameter<double> _clc = new GivenParameter<double>(0.1e-6);

        [ParameterName("cle"), ParameterInfo("Vdsat parameter for C-V model")]
        [Finite]
        private GivenParameter<double> _cle = new GivenParameter<double>(0.6);

        [ParameterName("dwc"), ParameterInfo("Delta W for C-V model")]
        [Finite]
        private GivenParameter<double> _dwc = new GivenParameter<double>();

        [ParameterName("dlc"), ParameterInfo("Delta L for C-V model")]
        [Finite]
        private GivenParameter<double> _dlc = new GivenParameter<double>();

        [ParameterName("hdif"), ParameterInfo("ACM Parameter: Distance Gate - contact")]
        [Finite]
        private GivenParameter<double> _hdif = new GivenParameter<double>(0.0);

        [ParameterName("ldif"), ParameterInfo("ACM Parameter: Length of LDD Gate-Source/Drain")]
        [Finite]
        private GivenParameter<double> _ldif = new GivenParameter<double>(0.0);

        [ParameterName("ld"), ParameterInfo("ACM Parameter: Length of LDD under Gate")]
        [Finite]
        private GivenParameter<double> _ld = new GivenParameter<double>(0.0);

        [ParameterName("rd"), ParameterInfo("ACM Parameter: Resistance of LDD drain side")]
        [Finite]
        private GivenParameter<double> _rd = new GivenParameter<double>(0.0);

        [ParameterName("rs"), ParameterInfo("ACM Parameter: Resistance of LDD source side")]
        [Finite]
        private GivenParameter<double> _rs = new GivenParameter<double>(0.0);

        [ParameterName("rdc"), ParameterInfo("ACM Parameter: Resistance contact drain side")]
        [Finite]
        private GivenParameter<double> _rdc = new GivenParameter<double>(0.0);

        [ParameterName("rsc"), ParameterInfo("ACM Parameter: Resistance contact source side")]
        [Finite]
        private GivenParameter<double> _rsc = new GivenParameter<double>(0.0);

        [ParameterName("wmlt"), ParameterInfo("ACM Parameter: Width shrink factor")]
        [Finite]
        private GivenParameter<double> _wmlt = new GivenParameter<double>(1.0);

        [ParameterName("lmlt"), ParameterInfo("Channel length shrink factor")]
        [Finite]
        private GivenParameter<double> _lmlt = new GivenParameter<double>(1.0);

        [ParameterName("alpha0"), ParameterInfo("substrate current model parameter")]
        [Finite]
        private GivenParameter<double> _alpha0 = new GivenParameter<double>(0.0);

        [ParameterName("alpha1"), ParameterInfo("substrate current model parameter")]
        [Finite]
        private GivenParameter<double> _alpha1 = new GivenParameter<double>(0.0);

        [ParameterName("beta0"), ParameterInfo("substrate current model parameter")]
        [Finite]
        private GivenParameter<double> _beta0 = new GivenParameter<double>(30.0);

        [ParameterName("ijth"), ParameterInfo("Diode limiting current")]
        [Finite]
        private GivenParameter<double> _ijth = new GivenParameter<double>(0.1);

        [ParameterName("lcdsc"), ParameterInfo("Length dependence of cdsc")]
        [Finite]
        private GivenParameter<double> _lcdsc = new GivenParameter<double>(0.0);

        [ParameterName("lcdscb"), ParameterInfo("Length dependence of cdscb")]
        [Finite]
        private GivenParameter<double> _lcdscb = new GivenParameter<double>(0.0);

        [ParameterName("lcdscd"), ParameterInfo("Length dependence of cdscd")]
        [Finite]
        private GivenParameter<double> _lcdscd = new GivenParameter<double>(0.0);

        [ParameterName("lcit"), ParameterInfo("Length dependence of cit")]
        [Finite]
        private GivenParameter<double> _lcit = new GivenParameter<double>(0.0);

        [ParameterName("lnfactor"), ParameterInfo("Length dependence of nfactor")]
        [Finite]
        private GivenParameter<double> _lnfactor = new GivenParameter<double>(0.0);

        [ParameterName("lxj"), ParameterInfo("Length dependence of xj")]
        [Finite]
        private GivenParameter<double> _lxj = new GivenParameter<double>(0.0);

        [ParameterName("lvsat"), ParameterInfo("Length dependence of vsat")]
        [Finite]
        private GivenParameter<double> _lvsat = new GivenParameter<double>(0.0);

        [ParameterName("lat"), ParameterInfo("Length dependence of at")]
        [Finite]
        private GivenParameter<double> _lat = new GivenParameter<double>(0.0);

        [ParameterName("la0"), ParameterInfo("Length dependence of a0")]
        [Finite]
        private GivenParameter<double> _la0 = new GivenParameter<double>(0.0);

        [ParameterName("lags"), ParameterInfo("Length dependence of ags")]
        [Finite]
        private GivenParameter<double> _lags = new GivenParameter<double>(0.0);

        [ParameterName("la1"), ParameterInfo("Length dependence of a1")]
        [Finite]
        private GivenParameter<double> _la1 = new GivenParameter<double>(0.0);

        [ParameterName("la2"), ParameterInfo("Length dependence of a2")]
        [Finite]
        private GivenParameter<double> _la2 = new GivenParameter<double>(0.0);

        [ParameterName("lketa"), ParameterInfo("Length dependence of keta")]
        [Finite]
        private GivenParameter<double> _lketa = new GivenParameter<double>(0.0);

        [ParameterName("lnsub"), ParameterInfo("Length dependence of nsub")]
        [Finite]
        private GivenParameter<double> _lnsub = new GivenParameter<double>(0.0);

        [ParameterName("lnch"), ParameterInfo("Length dependence of nch")]
        [Finite]
        public GivenParameter<double> Lnpeak
        {
            get => _lnpeak;
            set
            {
                Utility.Finite(value, nameof(Lnpeak));
                if (value > 1e20)
                    _lnpeak = value * 1e-6;
                else
                    _lnpeak = value;
            }
        }
        private GivenParameter<double> _lnpeak = new GivenParameter<double>();

        [ParameterName("lngate"), ParameterInfo("Length dependence of ngate")]
        [Finite]
        public GivenParameter<double> Lngate
        {
            get => _lngate;
            set
            {
                Utility.Finite(value, nameof(Lngate));
                if (value > 1.0e23)
                    _lngate = value * 1e-6;
                else
                    _lngate = value;
            }
        }
        private GivenParameter<double> _lngate = new GivenParameter<double>();

        [ParameterName("lgamma1"), ParameterInfo("Length dependence of gamma1")]
        [Finite]
        private GivenParameter<double> _lgamma1 = new GivenParameter<double>();

        [ParameterName("lgamma2"), ParameterInfo("Length dependence of gamma2")]
        [Finite]
        private GivenParameter<double> _lgamma2 = new GivenParameter<double>();

        [ParameterName("lvbx"), ParameterInfo("Length dependence of vbx")]
        [Finite]
        private GivenParameter<double> _lvbx = new GivenParameter<double>();

        [ParameterName("lvbm"), ParameterInfo("Length dependence of vbm")]
        [Finite]
        private GivenParameter<double> _lvbm = new GivenParameter<double>(0.0);

        [ParameterName("lxt"), ParameterInfo("Length dependence of xt")]
        [Finite]
        private GivenParameter<double> _lxt = new GivenParameter<double>(0.0);

        [ParameterName("lk1"), ParameterInfo("Length dependence of k1")]
        [Finite]
        private GivenParameter<double> _lk1 = new GivenParameter<double>();

        [ParameterName("lkt1"), ParameterInfo("Length dependence of kt1")]
        [Finite]
        private GivenParameter<double> _lkt1 = new GivenParameter<double>(0.0);

        [ParameterName("lkt1l"), ParameterInfo("Length dependence of kt1l")]
        [Finite]
        private GivenParameter<double> _lkt1l = new GivenParameter<double>(0.0);

        [ParameterName("lkt2"), ParameterInfo("Length dependence of kt2")]
        [Finite]
        private GivenParameter<double> _lkt2 = new GivenParameter<double>(0.0);

        [ParameterName("lk2"), ParameterInfo("Length dependence of k2")]
        [Finite]
        private GivenParameter<double> _lk2 = new GivenParameter<double>();

        [ParameterName("lk3"), ParameterInfo("Length dependence of k3")]
        [Finite]
        private GivenParameter<double> _lk3 = new GivenParameter<double>(0.0);

        [ParameterName("lk3b"), ParameterInfo("Length dependence of k3b")]
        [Finite]
        private GivenParameter<double> _lk3b = new GivenParameter<double>(0.0);

        [ParameterName("lw0"), ParameterInfo("Length dependence of w0")]
        [Finite]
        private GivenParameter<double> _lw0 = new GivenParameter<double>(0.0);

        [ParameterName("lnlx"), ParameterInfo("Length dependence of nlx")]
        [Finite]
        private GivenParameter<double> _lnlx = new GivenParameter<double>(0.0);

        [ParameterName("ldvt0"), ParameterInfo("Length dependence of dvt0")]
        [Finite]
        private GivenParameter<double> _ldvt0 = new GivenParameter<double>(0.0);

        [ParameterName("ldvt1"), ParameterInfo("Length dependence of dvt1")]
        [Finite]
        private GivenParameter<double> _ldvt1 = new GivenParameter<double>(0.0);

        [ParameterName("ldvt2"), ParameterInfo("Length dependence of dvt2")]
        [Finite]
        private GivenParameter<double> _ldvt2 = new GivenParameter<double>(0.0);

        [ParameterName("ldvt0w"), ParameterInfo("Length dependence of dvt0w")]
        [Finite]
        private GivenParameter<double> _ldvt0w = new GivenParameter<double>(0.0);

        [ParameterName("ldvt1w"), ParameterInfo("Length dependence of dvt1w")]
        [Finite]
        private GivenParameter<double> _ldvt1w = new GivenParameter<double>(0.0);

        [ParameterName("ldvt2w"), ParameterInfo("Length dependence of dvt2w")]
        [Finite]
        private GivenParameter<double> _ldvt2w = new GivenParameter<double>(0.0);

        [ParameterName("ldrout"), ParameterInfo("Length dependence of drout")]
        [Finite]
        private GivenParameter<double> _ldrout = new GivenParameter<double>(0.0);

        [ParameterName("ldsub"), ParameterInfo("Length dependence of dsub")]
        [Finite]
        private GivenParameter<double> _ldsub = new GivenParameter<double>(0.0);

        [ParameterName("lvtho"), ParameterInfo("Length dependence of vtho")]
        [Finite]
        private GivenParameter<double> _lvth0 = new GivenParameter<double>(0.0);

        [ParameterName("lua"), ParameterInfo("Length dependence of ua")]
        [Finite]
        private GivenParameter<double> _lua = new GivenParameter<double>(0.0);

        [ParameterName("lua1"), ParameterInfo("Length dependence of ua1")]
        [Finite]
        private GivenParameter<double> _lua1 = new GivenParameter<double>(0.0);

        [ParameterName("lub"), ParameterInfo("Length dependence of ub")]
        [Finite]
        private GivenParameter<double> _lub = new GivenParameter<double>(0.0);

        [ParameterName("lub1"), ParameterInfo("Length dependence of ub1")]
        [Finite]
        private GivenParameter<double> _lub1 = new GivenParameter<double>(0.0);

        [ParameterName("luc"), ParameterInfo("Length dependence of uc")]
        [Finite]
        private GivenParameter<double> _luc = new GivenParameter<double>(0.0);

        [ParameterName("luc1"), ParameterInfo("Length dependence of uc1")]
        [Finite]
        private GivenParameter<double> _luc1 = new GivenParameter<double>(0.0);

        [ParameterName("lu0"), ParameterInfo("Length dependence of u0")]
        [Finite]
        private GivenParameter<double> _lu0 = new GivenParameter<double>(0.0);

        [ParameterName("lute"), ParameterInfo("Length dependence of ute")]
        [Finite]
        private GivenParameter<double> _lute = new GivenParameter<double>(0.0);

        [ParameterName("lvoff"), ParameterInfo("Length dependence of voff")]
        [Finite]
        private GivenParameter<double> _lvoff = new GivenParameter<double>(0.0);

        [ParameterName("lelm"), ParameterInfo("Length dependence of elm")]
        [Finite]
        private GivenParameter<double> _lelm = new GivenParameter<double>(0.0);

        [ParameterName("ldelta"), ParameterInfo("Length dependence of delta")]
        [Finite]
        private GivenParameter<double> _ldelta = new GivenParameter<double>(0.0);

        [ParameterName("lrdsw"), ParameterInfo("Length dependence of rdsw ")]
        [Finite]
        private GivenParameter<double> _lrdsw = new GivenParameter<double>(0.0);

        [ParameterName("lprwg"), ParameterInfo("Length dependence of prwg ")]
        [Finite]
        private GivenParameter<double> _lprwg = new GivenParameter<double>(0.0);

        [ParameterName("lprwb"), ParameterInfo("Length dependence of prwb ")]
        [Finite]
        private GivenParameter<double> _lprwb = new GivenParameter<double>(0.0);

        [ParameterName("lprt"), ParameterInfo("Length dependence of prt ")]
        [Finite]
        private GivenParameter<double> _lprt = new GivenParameter<double>(0.0);

        [ParameterName("leta0"), ParameterInfo("Length dependence of eta0")]
        [Finite]
        private GivenParameter<double> _leta0 = new GivenParameter<double>(0.0);

        [ParameterName("letab"), ParameterInfo("Length dependence of etab")]
        [Finite]
        private GivenParameter<double> _letab = new GivenParameter<double>(-0.0);

        [ParameterName("lpclm"), ParameterInfo("Length dependence of pclm")]
        [Finite]
        private GivenParameter<double> _lpclm = new GivenParameter<double>(0.0);

        [ParameterName("lpdiblc1"), ParameterInfo("Length dependence of pdiblc1")]
        [Finite]
        private GivenParameter<double> _lpdibl1 = new GivenParameter<double>(0.0);

        [ParameterName("lpdiblc2"), ParameterInfo("Length dependence of pdiblc2")]
        [Finite]
        private GivenParameter<double> _lpdibl2 = new GivenParameter<double>(0.0);

        [ParameterName("lpdiblcb"), ParameterInfo("Length dependence of pdiblcb")]
        [Finite]
        private GivenParameter<double> _lpdiblb = new GivenParameter<double>(0.0);

        [ParameterName("lpscbe1"), ParameterInfo("Length dependence of pscbe1")]
        [Finite]
        private GivenParameter<double> _lpscbe1 = new GivenParameter<double>(0.0);

        [ParameterName("lpscbe2"), ParameterInfo("Length dependence of pscbe2")]
        [Finite]
        private GivenParameter<double> _lpscbe2 = new GivenParameter<double>(0.0);

        [ParameterName("lpvag"), ParameterInfo("Length dependence of pvag")]
        [Finite]
        private GivenParameter<double> _lpvag = new GivenParameter<double>(0.0);

        [ParameterName("lwr"), ParameterInfo("Length dependence of wr")]
        [Finite]
        private GivenParameter<double> _lwr = new GivenParameter<double>(0.0);

        [ParameterName("ldwg"), ParameterInfo("Length dependence of dwg")]
        [Finite]
        private GivenParameter<double> _ldwg = new GivenParameter<double>(0.0);

        [ParameterName("ldwb"), ParameterInfo("Length dependence of dwb")]
        [Finite]
        private GivenParameter<double> _ldwb = new GivenParameter<double>(0.0);

        [ParameterName("lb0"), ParameterInfo("Length dependence of b0")]
        [Finite]
        private GivenParameter<double> _lb0 = new GivenParameter<double>(0.0);

        [ParameterName("lb1"), ParameterInfo("Length dependence of b1")]
        [Finite]
        private GivenParameter<double> _lb1 = new GivenParameter<double>(0.0);

        [ParameterName("lcgsl"), ParameterInfo("Length dependence of cgsl")]
        [Finite]
        private GivenParameter<double> _lcgsl = new GivenParameter<double>(0.0);

        [ParameterName("lcgdl"), ParameterInfo("Length dependence of cgdl")]
        [Finite]
        private GivenParameter<double> _lcgdl = new GivenParameter<double>(0.0);

        [ParameterName("lckappa"), ParameterInfo("Length dependence of ckappa")]
        [Finite]
        private GivenParameter<double> _lckappa = new GivenParameter<double>(0.0);

        [ParameterName("lcf"), ParameterInfo("Length dependence of cf")]
        [Finite]
        private GivenParameter<double> _lcf = new GivenParameter<double>(0.0);

        [ParameterName("lclc"), ParameterInfo("Length dependence of clc")]
        [Finite]
        private GivenParameter<double> _lclc = new GivenParameter<double>(0.0);

        [ParameterName("lcle"), ParameterInfo("Length dependence of cle")]
        [Finite]
        private GivenParameter<double> _lcle = new GivenParameter<double>(0.0);

        [ParameterName("lalpha0"), ParameterInfo("Length dependence of alpha0")]
        [Finite]
        private GivenParameter<double> _lalpha0 = new GivenParameter<double>(0.0);

        [ParameterName("lalpha1"), ParameterInfo("Length dependence of alpha1")]
        [Finite]
        private GivenParameter<double> _lalpha1 = new GivenParameter<double>(0.0);

        [ParameterName("lbeta0"), ParameterInfo("Length dependence of beta0")]
        [Finite]
        private GivenParameter<double> _lbeta0 = new GivenParameter<double>(0.0);

        [ParameterName("lvfbcv"), ParameterInfo("Length dependence of vfbcv")]
        [Finite]
        private GivenParameter<double> _lvfbcv = new GivenParameter<double>(0.0);

        [ParameterName("lvfb"), ParameterInfo("Length dependence of vfb")]
        [Finite]
        private GivenParameter<double> _lvfb = new GivenParameter<double>(0.0);

        [ParameterName("lacde"), ParameterInfo("Length dependence of acde")]
        [Finite]
        private GivenParameter<double> _lacde = new GivenParameter<double>(0.0);

        [ParameterName("lmoin"), ParameterInfo("Length dependence of moin")]
        [Finite]
        private GivenParameter<double> _lmoin = new GivenParameter<double>(0.0);

        [ParameterName("lnoff"), ParameterInfo("Length dependence of noff")]
        [Finite]
        private GivenParameter<double> _lnoff = new GivenParameter<double>(0.0);

        [ParameterName("lvoffcv"), ParameterInfo("Length dependence of voffcv")]
        [Finite]
        private GivenParameter<double> _lvoffcv = new GivenParameter<double>(0.0);

        [ParameterName("wcdsc"), ParameterInfo("Width dependence of cdsc")]
        [Finite]
        private GivenParameter<double> _wcdsc = new GivenParameter<double>(0.0);

        [ParameterName("wcdscb"), ParameterInfo("Width dependence of cdscb")]
        [Finite]
        private GivenParameter<double> _wcdscb = new GivenParameter<double>(0.0);

        [ParameterName("wcdscd"), ParameterInfo("Width dependence of cdscd")]
        [Finite]
        private GivenParameter<double> _wcdscd = new GivenParameter<double>(0.0);

        [ParameterName("wcit"), ParameterInfo("Width dependence of cit")]
        [Finite]
        private GivenParameter<double> _wcit = new GivenParameter<double>(0.0);

        [ParameterName("wnfactor"), ParameterInfo("Width dependence of nfactor")]
        [Finite]
        private GivenParameter<double> _wnfactor = new GivenParameter<double>(0.0);

        [ParameterName("wxj"), ParameterInfo("Width dependence of xj")]
        [Finite]
        private GivenParameter<double> _wxj = new GivenParameter<double>(0.0);

        [ParameterName("wvsat"), ParameterInfo("Width dependence of vsat")]
        [Finite]
        private GivenParameter<double> _wvsat = new GivenParameter<double>(0.0);

        [ParameterName("wat"), ParameterInfo("Width dependence of at")]
        [Finite]
        private GivenParameter<double> _wat = new GivenParameter<double>(0.0);

        [ParameterName("wa0"), ParameterInfo("Width dependence of a0")]
        [Finite]
        private GivenParameter<double> _wa0 = new GivenParameter<double>(0.0);

        [ParameterName("wags"), ParameterInfo("Width dependence of ags")]
        [Finite]
        private GivenParameter<double> _wags = new GivenParameter<double>(0.0);

        [ParameterName("wa1"), ParameterInfo("Width dependence of a1")]
        [Finite]
        private GivenParameter<double> _wa1 = new GivenParameter<double>(0.0);

        [ParameterName("wa2"), ParameterInfo("Width dependence of a2")]
        [Finite]
        private GivenParameter<double> _wa2 = new GivenParameter<double>(0.0);

        [ParameterName("wketa"), ParameterInfo("Width dependence of keta")]
        [Finite]
        private GivenParameter<double> _wketa = new GivenParameter<double>(0.0);

        [ParameterName("wnsub"), ParameterInfo("Width dependence of nsub")]
        [Finite]
        private GivenParameter<double> _wnsub = new GivenParameter<double>(0.0);

        [ParameterName("wnch"), ParameterInfo("Width dependence of nch")]
        [Finite]
        public GivenParameter<double> Wnpeak
        {
            get => _wnpeak;
            set
            {
                Utility.Finite(value, nameof(Wnpeak));
                if (value > 1e20)
                    _wnpeak = value * 1e-6;
                else
                    _wnpeak = value;
            }
        }
        private GivenParameter<double> _wnpeak = new GivenParameter<double>();

        [ParameterName("wngate"), ParameterInfo("Width dependence of ngate")]
        [Finite]
        public GivenParameter<double> Wngate
        {
            get => _wngate;
            set
            {
                Utility.Finite(value, nameof(Wngate));
                if (value > 1e23)
                    _wngate = value * 1e-6;
                else
                    _wngate = value;
            }
        }
        private GivenParameter<double> _wngate = new GivenParameter<double>();

        [ParameterName("wgamma1"), ParameterInfo("Width dependence of gamma1")]
        [Finite]
        private GivenParameter<double> _wgamma1 = new GivenParameter<double>();

        [ParameterName("wgamma2"), ParameterInfo("Width dependence of gamma2")]
        [Finite]
        private GivenParameter<double> _wgamma2 = new GivenParameter<double>();

        [ParameterName("wvbx"), ParameterInfo("Width dependence of vbx")]
        [Finite]
        private GivenParameter<double> _wvbx = new GivenParameter<double>();

        [ParameterName("wvbm"), ParameterInfo("Width dependence of vbm")]
        [Finite]
        private GivenParameter<double> _wvbm = new GivenParameter<double>(0.0);

        [ParameterName("wxt"), ParameterInfo("Width dependence of xt")]
        [Finite]
        private GivenParameter<double> _wxt = new GivenParameter<double>(0.0);

        [ParameterName("wk1"), ParameterInfo("Width dependence of k1")]
        [Finite]
        private GivenParameter<double> _wk1 = new GivenParameter<double>();

        [ParameterName("wkt1"), ParameterInfo("Width dependence of kt1")]
        [Finite]
        private GivenParameter<double> _wkt1 = new GivenParameter<double>(0.0);

        [ParameterName("wkt1l"), ParameterInfo("Width dependence of kt1l")]
        [Finite]
        private GivenParameter<double> _wkt1l = new GivenParameter<double>(0.0);

        [ParameterName("wkt2"), ParameterInfo("Width dependence of kt2")]
        [Finite]
        private GivenParameter<double> _wkt2 = new GivenParameter<double>(0.0);

        [ParameterName("wk2"), ParameterInfo("Width dependence of k2")]
        [Finite]
        private GivenParameter<double> _wk2 = new GivenParameter<double>();

        [ParameterName("wk3"), ParameterInfo("Width dependence of k3")]
        [Finite]
        private GivenParameter<double> _wk3 = new GivenParameter<double>(0.0);

        [ParameterName("wk3b"), ParameterInfo("Width dependence of k3b")]
        [Finite]
        private GivenParameter<double> _wk3b = new GivenParameter<double>(0.0);

        [ParameterName("ww0"), ParameterInfo("Width dependence of w0")]
        [Finite]
        private GivenParameter<double> _ww0 = new GivenParameter<double>(0.0);

        [ParameterName("wnlx"), ParameterInfo("Width dependence of nlx")]
        [Finite]
        private GivenParameter<double> _wnlx = new GivenParameter<double>(0.0);

        [ParameterName("wdvt0"), ParameterInfo("Width dependence of dvt0")]
        [Finite]
        private GivenParameter<double> _wdvt0 = new GivenParameter<double>(0.0);

        [ParameterName("wdvt1"), ParameterInfo("Width dependence of dvt1")]
        [Finite]
        private GivenParameter<double> _wdvt1 = new GivenParameter<double>(0.0);

        [ParameterName("wdvt2"), ParameterInfo("Width dependence of dvt2")]
        [Finite]
        private GivenParameter<double> _wdvt2 = new GivenParameter<double>(0.0);

        [ParameterName("wdvt0w"), ParameterInfo("Width dependence of dvt0w")]
        [Finite]
        private GivenParameter<double> _wdvt0w = new GivenParameter<double>(0.0);

        [ParameterName("wdvt1w"), ParameterInfo("Width dependence of dvt1w")]
        [Finite]
        private GivenParameter<double> _wdvt1w = new GivenParameter<double>(0.0);

        [ParameterName("wdvt2w"), ParameterInfo("Width dependence of dvt2w")]
        [Finite]
        private GivenParameter<double> _wdvt2w = new GivenParameter<double>(0.0);

        [ParameterName("wdrout"), ParameterInfo("Width dependence of drout")]
        [Finite]
        private GivenParameter<double> _wdrout = new GivenParameter<double>(0.0);

        [ParameterName("wdsub"), ParameterInfo("Width dependence of dsub")]
        [Finite]
        private GivenParameter<double> _wdsub = new GivenParameter<double>(0.0);

        [ParameterName("wvtho"), ParameterInfo("Width dependence of vtho")]
        [Finite]
        private GivenParameter<double> _wvth0 = new GivenParameter<double>(0.0);

        [ParameterName("wua"), ParameterInfo("Width dependence of ua")]
        [Finite]
        private GivenParameter<double> _wua = new GivenParameter<double>(0.0);

        [ParameterName("wua1"), ParameterInfo("Width dependence of ua1")]
        [Finite]
        private GivenParameter<double> _wua1 = new GivenParameter<double>(0.0);

        [ParameterName("wub"), ParameterInfo("Width dependence of ub")]
        [Finite]
        private GivenParameter<double> _wub = new GivenParameter<double>(0.0);

        [ParameterName("wub1"), ParameterInfo("Width dependence of ub1")]
        [Finite]
        private GivenParameter<double> _wub1 = new GivenParameter<double>(0.0);

        [ParameterName("wuc"), ParameterInfo("Width dependence of uc")]
        [Finite]
        private GivenParameter<double> _wuc = new GivenParameter<double>(0.0);

        [ParameterName("wuc1"), ParameterInfo("Width dependence of uc1")]
        [Finite]
        private GivenParameter<double> _wuc1 = new GivenParameter<double>(0.0);

        [ParameterName("wu0"), ParameterInfo("Width dependence of u0")]
        [Finite]
        private GivenParameter<double> _wu0 = new GivenParameter<double>(0.0);

        [ParameterName("wute"), ParameterInfo("Width dependence of ute")]
        [Finite]
        private GivenParameter<double> _wute = new GivenParameter<double>(0.0);

        [ParameterName("wvoff"), ParameterInfo("Width dependence of voff")]
        [Finite]
        private GivenParameter<double> _wvoff = new GivenParameter<double>(0.0);

        [ParameterName("welm"), ParameterInfo("Width dependence of elm")]
        [Finite]
        private GivenParameter<double> _welm = new GivenParameter<double>(0.0);

        [ParameterName("wdelta"), ParameterInfo("Width dependence of delta")]
        [Finite]
        private GivenParameter<double> _wdelta = new GivenParameter<double>(0.0);

        [ParameterName("wrdsw"), ParameterInfo("Width dependence of rdsw ")]
        [Finite]
        private GivenParameter<double> _wrdsw = new GivenParameter<double>(0.0);

        [ParameterName("wprwg"), ParameterInfo("Width dependence of prwg ")]
        [Finite]
        private GivenParameter<double> _wprwg = new GivenParameter<double>(0.0);

        [ParameterName("wprwb"), ParameterInfo("Width dependence of prwb ")]
        [Finite]
        private GivenParameter<double> _wprwb = new GivenParameter<double>(0.0);

        [ParameterName("wprt"), ParameterInfo("Width dependence of prt")]
        [Finite]
        private GivenParameter<double> _wprt = new GivenParameter<double>(0.0);

        [ParameterName("weta0"), ParameterInfo("Width dependence of eta0")]
        [Finite]
        private GivenParameter<double> _weta0 = new GivenParameter<double>(0.0);

        [ParameterName("wetab"), ParameterInfo("Width dependence of etab")]
        [Finite]
        private GivenParameter<double> _wetab = new GivenParameter<double>(0.0);

        [ParameterName("wpclm"), ParameterInfo("Width dependence of pclm")]
        [Finite]
        private GivenParameter<double> _wpclm = new GivenParameter<double>(0.0);

        [ParameterName("wpdiblc1"), ParameterInfo("Width dependence of pdiblc1")]
        [Finite]
        private GivenParameter<double> _wpdibl1 = new GivenParameter<double>(0.0);

        [ParameterName("wpdiblc2"), ParameterInfo("Width dependence of pdiblc2")]
        [Finite]
        private GivenParameter<double> _wpdibl2 = new GivenParameter<double>(0.0);

        [ParameterName("wpdiblcb"), ParameterInfo("Width dependence of pdiblcb")]
        [Finite]
        private GivenParameter<double> _wpdiblb = new GivenParameter<double>(0.0);

        [ParameterName("wpscbe1"), ParameterInfo("Width dependence of pscbe1")]
        [Finite]
        private GivenParameter<double> _wpscbe1 = new GivenParameter<double>(0.0);

        [ParameterName("wpscbe2"), ParameterInfo("Width dependence of pscbe2")]
        [Finite]
        private GivenParameter<double> _wpscbe2 = new GivenParameter<double>(0.0);

        [ParameterName("wpvag"), ParameterInfo("Width dependence of pvag")]
        [Finite]
        private GivenParameter<double> _wpvag = new GivenParameter<double>(0.0);

        [ParameterName("wwr"), ParameterInfo("Width dependence of wr")]
        [Finite]
        private GivenParameter<double> _wwr = new GivenParameter<double>(0.0);

        [ParameterName("wdwg"), ParameterInfo("Width dependence of dwg")]
        [Finite]
        private GivenParameter<double> _wdwg = new GivenParameter<double>(0.0);

        [ParameterName("wdwb"), ParameterInfo("Width dependence of dwb")]
        [Finite]
        private GivenParameter<double> _wdwb = new GivenParameter<double>(0.0);

        [ParameterName("wb0"), ParameterInfo("Width dependence of b0")]
        [Finite]
        private GivenParameter<double> _wb0 = new GivenParameter<double>(0.0);

        [ParameterName("wb1"), ParameterInfo("Width dependence of b1")]
        [Finite]
        private GivenParameter<double> _wb1 = new GivenParameter<double>(0.0);

        [ParameterName("wcgsl"), ParameterInfo("Width dependence of cgsl")]
        [Finite]
        private GivenParameter<double> _wcgsl = new GivenParameter<double>(0.0);

        [ParameterName("wcgdl"), ParameterInfo("Width dependence of cgdl")]
        [Finite]
        private GivenParameter<double> _wcgdl = new GivenParameter<double>(0.0);

        [ParameterName("wckappa"), ParameterInfo("Width dependence of ckappa")]
        [Finite]
        private GivenParameter<double> _wckappa = new GivenParameter<double>(0.0);

        [ParameterName("wcf"), ParameterInfo("Width dependence of cf")]
        [Finite]
        private GivenParameter<double> _wcf = new GivenParameter<double>(0.0);

        [ParameterName("wclc"), ParameterInfo("Width dependence of clc")]
        [Finite]
        private GivenParameter<double> _wclc = new GivenParameter<double>(0.0);

        [ParameterName("wcle"), ParameterInfo("Width dependence of cle")]
        [Finite]
        private GivenParameter<double> _wcle = new GivenParameter<double>(0.0);

        [ParameterName("walpha0"), ParameterInfo("Width dependence of alpha0")]
        [Finite]
        private GivenParameter<double> _walpha0 = new GivenParameter<double>(0.0);

        [ParameterName("walpha1"), ParameterInfo("Width dependence of alpha1")]
        [Finite]
        private GivenParameter<double> _walpha1 = new GivenParameter<double>(0.0);

        [ParameterName("wbeta0"), ParameterInfo("Width dependence of beta0")]
        [Finite]
        private GivenParameter<double> _wbeta0 = new GivenParameter<double>(0.0);

        [ParameterName("wvfbcv"), ParameterInfo("Width dependence of vfbcv")]
        [Finite]
        private GivenParameter<double> _wvfbcv = new GivenParameter<double>(0.0);

        [ParameterName("wvfb"), ParameterInfo("Width dependence of vfb")]
        [Finite]
        private GivenParameter<double> _wvfb = new GivenParameter<double>(0.0);

        [ParameterName("wacde"), ParameterInfo("Width dependence of acde")]
        [Finite]
        private GivenParameter<double> _wacde = new GivenParameter<double>(0.0);

        [ParameterName("wmoin"), ParameterInfo("Width dependence of moin")]
        [Finite]
        private GivenParameter<double> _wmoin = new GivenParameter<double>(0.0);

        [ParameterName("wnoff"), ParameterInfo("Width dependence of noff")]
        [Finite]
        private GivenParameter<double> _wnoff = new GivenParameter<double>(0.0);

        [ParameterName("wvoffcv"), ParameterInfo("Width dependence of voffcv")]
        [Finite]
        private GivenParameter<double> _wvoffcv = new GivenParameter<double>(0.0);

        [ParameterName("pcdsc"), ParameterInfo("Cross-term dependence of cdsc")]
        [Finite]
        private GivenParameter<double> _pcdsc = new GivenParameter<double>(0.0);

        [ParameterName("pcdscb"), ParameterInfo("Cross-term dependence of cdscb")]
        [Finite]
        private GivenParameter<double> _pcdscb = new GivenParameter<double>(0.0);

        [ParameterName("pcdscd"), ParameterInfo("Cross-term dependence of cdscd")]
        [Finite]
        private GivenParameter<double> _pcdscd = new GivenParameter<double>(0.0);

        [ParameterName("pcit"), ParameterInfo("Cross-term dependence of cit")]
        [Finite]
        private GivenParameter<double> _pcit = new GivenParameter<double>(0.0);

        [ParameterName("pnfactor"), ParameterInfo("Cross-term dependence of nfactor")]
        [Finite]
        private GivenParameter<double> _pnfactor = new GivenParameter<double>(0.0);

        [ParameterName("pxj"), ParameterInfo("Cross-term dependence of xj")]
        [Finite]
        private GivenParameter<double> _pxj = new GivenParameter<double>(0.0);

        [ParameterName("pvsat"), ParameterInfo("Cross-term dependence of vsat")]
        [Finite]
        private GivenParameter<double> _pvsat = new GivenParameter<double>(0.0);

        [ParameterName("pat"), ParameterInfo("Cross-term dependence of at")]
        [Finite]
        private GivenParameter<double> _pat = new GivenParameter<double>(0.0);

        [ParameterName("pa0"), ParameterInfo("Cross-term dependence of a0")]
        [Finite]
        private GivenParameter<double> _pa0 = new GivenParameter<double>(0.0);

        [ParameterName("pags"), ParameterInfo("Cross-term dependence of ags")]
        [Finite]
        private GivenParameter<double> _pags = new GivenParameter<double>(0.0);

        [ParameterName("pa1"), ParameterInfo("Cross-term dependence of a1")]
        [Finite]
        private GivenParameter<double> _pa1 = new GivenParameter<double>(0.0);

        [ParameterName("pa2"), ParameterInfo("Cross-term dependence of a2")]
        [Finite]
        private GivenParameter<double> _pa2 = new GivenParameter<double>(0.0);

        [ParameterName("pketa"), ParameterInfo("Cross-term dependence of keta")]
        [Finite]
        private GivenParameter<double> _pketa = new GivenParameter<double>(0.0);

        [ParameterName("pnsub"), ParameterInfo("Cross-term dependence of nsub")]
        [Finite]
        private GivenParameter<double> _pnsub = new GivenParameter<double>(0.0);

        [ParameterName("pnch"), ParameterInfo("Cross-term dependence of nch")]
        [Finite]
        public GivenParameter<double> Pnpeak
        {
            get => _pnpeak;
            set
            {
                Utility.Finite(value, nameof(Pnpeak));
                if (value > 1e20)
                    _pnpeak = value * 1e-6;
                else
                    _pnpeak = value;
            }
        }
        private GivenParameter<double> _pnpeak = new GivenParameter<double>();

        [ParameterName("pngate"), ParameterInfo("Cross-term dependence of ngate")]
        [Finite]
        public GivenParameter<double> Pngate
        {
            get => _pngate;
            set
            {
                Utility.Finite(value, nameof(Pngate));
                if (value > 1e23)
                    _pngate = value * 1e-6;
                else
                    _pngate = value;
            }
        }
        private GivenParameter<double> _pngate = new GivenParameter<double>();

        [ParameterName("pgamma1"), ParameterInfo("Cross-term dependence of gamma1")]
        [Finite]
        private GivenParameter<double> _pgamma1 = new GivenParameter<double>();

        [ParameterName("pgamma2"), ParameterInfo("Cross-term dependence of gamma2")]
        [Finite]
        private GivenParameter<double> _pgamma2 = new GivenParameter<double>();

        [ParameterName("pvbx"), ParameterInfo("Cross-term dependence of vbx")]
        [Finite]
        private GivenParameter<double> _pvbx = new GivenParameter<double>();

        [ParameterName("pvbm"), ParameterInfo("Cross-term dependence of vbm")]
        [Finite]
        private GivenParameter<double> _pvbm = new GivenParameter<double>(0.0);

        [ParameterName("pxt"), ParameterInfo("Cross-term dependence of xt")]
        [Finite]
        private GivenParameter<double> _pxt = new GivenParameter<double>(0.0);

        [ParameterName("pk1"), ParameterInfo("Cross-term dependence of k1")]
        [Finite]
        private GivenParameter<double> _pk1 = new GivenParameter<double>();

        [ParameterName("pkt1"), ParameterInfo("Cross-term dependence of kt1")]
        [Finite]
        private GivenParameter<double> _pkt1 = new GivenParameter<double>(0.0);

        [ParameterName("pkt1l"), ParameterInfo("Cross-term dependence of kt1l")]
        [Finite]
        private GivenParameter<double> _pkt1l = new GivenParameter<double>(0.0);

        [ParameterName("pkt2"), ParameterInfo("Cross-term dependence of kt2")]
        [Finite]
        private GivenParameter<double> _pkt2 = new GivenParameter<double>(0.0);

        [ParameterName("pk2"), ParameterInfo("Cross-term dependence of k2")]
        [Finite]
        private GivenParameter<double> _pk2 = new GivenParameter<double>();

        [ParameterName("pk3"), ParameterInfo("Cross-term dependence of k3")]
        [Finite]
        private GivenParameter<double> _pk3 = new GivenParameter<double>(0.0);

        [ParameterName("pk3b"), ParameterInfo("Cross-term dependence of k3b")]
        [Finite]
        private GivenParameter<double> _pk3b = new GivenParameter<double>(0.0);

        [ParameterName("pw0"), ParameterInfo("Cross-term dependence of w0")]
        [Finite]
        private GivenParameter<double> _pw0 = new GivenParameter<double>(0.0);

        [ParameterName("pnlx"), ParameterInfo("Cross-term dependence of nlx")]
        [Finite]
        private GivenParameter<double> _pnlx = new GivenParameter<double>(0.0);

        [ParameterName("pdvt0"), ParameterInfo("Cross-term dependence of dvt0")]
        [Finite]
        private GivenParameter<double> _pdvt0 = new GivenParameter<double>(0.0);

        [ParameterName("pdvt1"), ParameterInfo("Cross-term dependence of dvt1")]
        [Finite]
        private GivenParameter<double> _pdvt1 = new GivenParameter<double>(0.0);

        [ParameterName("pdvt2"), ParameterInfo("Cross-term dependence of dvt2")]
        [Finite]
        private GivenParameter<double> _pdvt2 = new GivenParameter<double>(0.0);

        [ParameterName("pdvt0w"), ParameterInfo("Cross-term dependence of dvt0w")]
        [Finite]
        private GivenParameter<double> _pdvt0w = new GivenParameter<double>(0.0);

        [ParameterName("pdvt1w"), ParameterInfo("Cross-term dependence of dvt1w")]
        [Finite]
        private GivenParameter<double> _pdvt1w = new GivenParameter<double>(0.0);

        [ParameterName("pdvt2w"), ParameterInfo("Cross-term dependence of dvt2w")]
        [Finite]
        private GivenParameter<double> _pdvt2w = new GivenParameter<double>(0.0);

        [ParameterName("pdrout"), ParameterInfo("Cross-term dependence of drout")]
        [Finite]
        private GivenParameter<double> _pdrout = new GivenParameter<double>(0.0);

        [ParameterName("pdsub"), ParameterInfo("Cross-term dependence of dsub")]
        [Finite]
        private GivenParameter<double> _pdsub = new GivenParameter<double>(0.0);

        [ParameterName("pvtho"), ParameterInfo("Cross-term dependence of vtho")]
        [Finite]
        private GivenParameter<double> _pvth0 = new GivenParameter<double>(0.0);

        [ParameterName("pua"), ParameterInfo("Cross-term dependence of ua")]
        [Finite]
        private GivenParameter<double> _pua = new GivenParameter<double>(0.0);

        [ParameterName("pua1"), ParameterInfo("Cross-term dependence of ua1")]
        [Finite]
        private GivenParameter<double> _pua1 = new GivenParameter<double>(0.0);

        [ParameterName("pub"), ParameterInfo("Cross-term dependence of ub")]
        [Finite]
        private GivenParameter<double> _pub = new GivenParameter<double>(0.0);

        [ParameterName("pub1"), ParameterInfo("Cross-term dependence of ub1")]
        [Finite]
        private GivenParameter<double> _pub1 = new GivenParameter<double>(0.0);

        [ParameterName("puc"), ParameterInfo("Cross-term dependence of uc")]
        [Finite]
        private GivenParameter<double> _puc = new GivenParameter<double>(0.0);

        [ParameterName("puc1"), ParameterInfo("Cross-term dependence of uc1")]
        [Finite]
        private GivenParameter<double> _puc1 = new GivenParameter<double>(0.0);

        [ParameterName("pu0"), ParameterInfo("Cross-term dependence of u0")]
        [Finite]
        private GivenParameter<double> _pu0 = new GivenParameter<double>(0.0);

        [ParameterName("pute"), ParameterInfo("Cross-term dependence of ute")]
        [Finite]
        private GivenParameter<double> _pute = new GivenParameter<double>(0.0);

        [ParameterName("pvoff"), ParameterInfo("Cross-term dependence of voff")]
        [Finite]
        private GivenParameter<double> _pvoff = new GivenParameter<double>(0.0);

        [ParameterName("pelm"), ParameterInfo("Cross-term dependence of elm")]
        [Finite]
        private GivenParameter<double> _pelm = new GivenParameter<double>(0.0);

        [ParameterName("pdelta"), ParameterInfo("Cross-term dependence of delta")]
        [Finite]
        private GivenParameter<double> _pdelta = new GivenParameter<double>(0.0);

        [ParameterName("prdsw"), ParameterInfo("Cross-term dependence of rdsw ")]
        [Finite]
        private GivenParameter<double> _prdsw = new GivenParameter<double>(0.0);

        [ParameterName("pprwg"), ParameterInfo("Cross-term dependence of prwg ")]
        [Finite]
        private GivenParameter<double> _pprwg = new GivenParameter<double>(0.0);

        [ParameterName("pprwb"), ParameterInfo("Cross-term dependence of prwb ")]
        [Finite]
        private GivenParameter<double> _pprwb = new GivenParameter<double>(0.0);

        [ParameterName("pprt"), ParameterInfo("Cross-term dependence of prt ")]
        [Finite]
        private GivenParameter<double> _pprt = new GivenParameter<double>(0.0);

        [ParameterName("peta0"), ParameterInfo("Cross-term dependence of eta0")]
        [Finite]
        private GivenParameter<double> _peta0 = new GivenParameter<double>(0.0);

        [ParameterName("petab"), ParameterInfo("Cross-term dependence of etab")]
        [Finite]
        private GivenParameter<double> _petab = new GivenParameter<double>(0.0);

        [ParameterName("ppclm"), ParameterInfo("Cross-term dependence of pclm")]
        [Finite]
        private GivenParameter<double> _ppclm = new GivenParameter<double>(0.0);

        [ParameterName("ppdiblc1"), ParameterInfo("Cross-term dependence of pdiblc1")]
        [Finite]
        private GivenParameter<double> _ppdibl1 = new GivenParameter<double>(0.0);

        [ParameterName("ppdiblc2"), ParameterInfo("Cross-term dependence of pdiblc2")]
        [Finite]
        private GivenParameter<double> _ppdibl2 = new GivenParameter<double>(0.0);

        [ParameterName("ppdiblcb"), ParameterInfo("Cross-term dependence of pdiblcb")]
        [Finite]
        private GivenParameter<double> _ppdiblb = new GivenParameter<double>(0.0);

        [ParameterName("ppscbe1"), ParameterInfo("Cross-term dependence of pscbe1")]
        [Finite]
        private GivenParameter<double> _ppscbe1 = new GivenParameter<double>(0.0);

        [ParameterName("ppscbe2"), ParameterInfo("Cross-term dependence of pscbe2")]
        [Finite]
        private GivenParameter<double> _ppscbe2 = new GivenParameter<double>(0.0);

        [ParameterName("ppvag"), ParameterInfo("Cross-term dependence of pvag")]
        [Finite]
        private GivenParameter<double> _ppvag = new GivenParameter<double>(0.0);

        [ParameterName("pwr"), ParameterInfo("Cross-term dependence of wr")]
        [Finite]
        private GivenParameter<double> _pwr = new GivenParameter<double>(0.0);

        [ParameterName("pdwg"), ParameterInfo("Cross-term dependence of dwg")]
        [Finite]
        private GivenParameter<double> _pdwg = new GivenParameter<double>(0.0);

        [ParameterName("pdwb"), ParameterInfo("Cross-term dependence of dwb")]
        [Finite]
        private GivenParameter<double> _pdwb = new GivenParameter<double>(0.0);

        [ParameterName("pb0"), ParameterInfo("Cross-term dependence of b0")]
        [Finite]
        private GivenParameter<double> _pb0 = new GivenParameter<double>(0.0);

        [ParameterName("pb1"), ParameterInfo("Cross-term dependence of b1")]
        [Finite]
        private GivenParameter<double> _pb1 = new GivenParameter<double>(0.0);

        [ParameterName("pcgsl"), ParameterInfo("Cross-term dependence of cgsl")]
        [Finite]
        private GivenParameter<double> _pcgsl = new GivenParameter<double>(0.0);

        [ParameterName("pcgdl"), ParameterInfo("Cross-term dependence of cgdl")]
        [Finite]
        private GivenParameter<double> _pcgdl = new GivenParameter<double>(0.0);

        [ParameterName("pckappa"), ParameterInfo("Cross-term dependence of ckappa")]
        [Finite]
        private GivenParameter<double> _pckappa = new GivenParameter<double>(0.0);

        [ParameterName("pcf"), ParameterInfo("Cross-term dependence of cf")]
        [Finite]
        private GivenParameter<double> _pcf = new GivenParameter<double>(0.0);

        [ParameterName("pclc"), ParameterInfo("Cross-term dependence of clc")]
        [Finite]
        private GivenParameter<double> _pclc = new GivenParameter<double>(0.0);

        [ParameterName("pcle"), ParameterInfo("Cross-term dependence of cle")]
        [Finite]
        private GivenParameter<double> _pcle = new GivenParameter<double>(0.0);

        [ParameterName("palpha0"), ParameterInfo("Cross-term dependence of alpha0")]
        [Finite]
        private GivenParameter<double> _palpha0 = new GivenParameter<double>(0.0);

        [ParameterName("palpha1"), ParameterInfo("Cross-term dependence of alpha1")]
        [Finite]
        private GivenParameter<double> _palpha1 = new GivenParameter<double>(0.0);

        [ParameterName("pbeta0"), ParameterInfo("Cross-term dependence of beta0")]
        [Finite]
        private GivenParameter<double> _pbeta0 = new GivenParameter<double>(0.0);

        [ParameterName("pvfbcv"), ParameterInfo("Cross-term dependence of vfbcv")]
        [Finite]
        private GivenParameter<double> _pvfbcv = new GivenParameter<double>(0.0);

        [ParameterName("pvfb"), ParameterInfo("Cross-term dependence of vfb")]
        [Finite]
        private GivenParameter<double> _pvfb = new GivenParameter<double>(0.0);

        [ParameterName("pacde"), ParameterInfo("Cross-term dependence of acde")]
        [Finite]
        private GivenParameter<double> _pacde = new GivenParameter<double>(0.0);

        [ParameterName("pmoin"), ParameterInfo("Cross-term dependence of moin")]
        [Finite]
        private GivenParameter<double> _pmoin = new GivenParameter<double>(0.0);

        [ParameterName("pnoff"), ParameterInfo("Cross-term dependence of noff")]
        [Finite]
        private GivenParameter<double> _pnoff = new GivenParameter<double>(0.0);

        [ParameterName("pvoffcv"), ParameterInfo("Cross-term dependence of voffcv")]
        [Finite]
        private GivenParameter<double> _pvoffcv = new GivenParameter<double>(0.0);

        [ParameterName("noia"), ParameterInfo("Flicker noise parameter")]
        [Finite]
        private GivenParameter<double> _oxideTrapDensityA = new GivenParameter<double>();

        [ParameterName("noib"), ParameterInfo("Flicker noise parameter")]
        [Finite]
        private GivenParameter<double> _oxideTrapDensityB = new GivenParameter<double>();

        [ParameterName("noic"), ParameterInfo("Flicker noise parameter")]
        [Finite]
        private GivenParameter<double> _oxideTrapDensityC = new GivenParameter<double>();

        [ParameterName("em"), ParameterInfo("Flicker noise parameter")]
        [Finite]
        private GivenParameter<double> _em = new GivenParameter<double>(4.1e7);

        [ParameterName("ef"), ParameterInfo("Flicker noise frequency exponent")]
        [Finite]
        private GivenParameter<double> _ef = new GivenParameter<double>(1.0);

        [ParameterName("af"), ParameterInfo("Flicker noise exponent")]
        [Finite]
        private GivenParameter<double> _af = new GivenParameter<double>(1.0);

        [ParameterName("kf"), ParameterInfo("Flicker noise coefficient")]
        [Finite]
        private GivenParameter<double> _kf = new GivenParameter<double>(0.0);

        [ParameterName("vgs_max"), ParameterInfo("maximum voltage G-S branch")]
        [Finite]
        private GivenParameter<double> _vgsMax = new GivenParameter<double>(1e99);

        [ParameterName("vgd_max"), ParameterInfo("maximum voltage G-D branch")]
        [Finite]
        private GivenParameter<double> _vgdMax = new GivenParameter<double>(1e99);

        [ParameterName("vgb_max"), ParameterInfo("maximum voltage G-B branch")]
        [Finite]
        private GivenParameter<double> _vgbMax = new GivenParameter<double>(1e99);

        [ParameterName("vds_max"), ParameterInfo("maximum voltage D-S branch")]
        [Finite]
        private GivenParameter<double> _vdsMax = new GivenParameter<double>(1e99);

        [ParameterName("vbs_max"), ParameterInfo("maximum voltage B-S branch")]
        [Finite]
        private GivenParameter<double> _vbsMax = new GivenParameter<double>(1e99);

        [ParameterName("vbd_max"), ParameterInfo("maximum voltage B-D branch")]
        [Finite]
        private GivenParameter<double> _vbdMax = new GivenParameter<double>(1e99);

        [ParameterName("vgsr_max"), ParameterInfo("maximum voltage G-S branch")]
        [Finite]
        private GivenParameter<double> _vgsrMax = new GivenParameter<double>(1e99);

        [ParameterName("vgdr_max"), ParameterInfo("maximum voltage G-D branch")]
        [Finite]
        private GivenParameter<double> _vgdrMax = new GivenParameter<double>(1e99);

        [ParameterName("vgbr_max"), ParameterInfo("maximum voltage G-B branch")]
        [Finite]
        private GivenParameter<double> _vgbrMax = new GivenParameter<double>(1e99);

        [ParameterName("vbsr_max"), ParameterInfo("maximum voltage B-S branch")]
        [Finite]
        private GivenParameter<double> _vbsrMax = new GivenParameter<double>(1e99);

        [ParameterName("vbdr_max"), ParameterInfo("maximum voltage B-D branch")]
        [Finite]
        private GivenParameter<double> _vbdrMax = new GivenParameter<double>(1e99);

        [ParameterName("tnom"), ParameterInfo("Parameter measurement temperature")]
        [DerivedProperty, Finite, GreaterThan(-Constants.CelsiusKelvin)]
        public double TnomCelsius
        {
            get => Tnom - Constants.CelsiusKelvin;
            set => Tnom = value + Constants.CelsiusKelvin;
        }

        [GreaterThan(0), Finite]
        private GivenParameter<double> _tnom = new GivenParameter<double>(300.15);

        /// <summary>
        /// Gets the mosfet type (1 = NMOS, -1 = PMOS).
        /// </summary>
        public double Type { get; private set; } = 1.0;

        [ParameterName("nmos"), ParameterInfo("Flag to indicate NMOS")]
        public void SetNmos(bool flag = true)
        {
            if (flag)
                Type = 1.0;
        }

        [ParameterName("pmos"), ParameterInfo("Flag to indicate PMOS")]
        public void SetPmos(bool flag = true)
        {
            if (flag)
                Type = -1.0;
        }
    }
}
