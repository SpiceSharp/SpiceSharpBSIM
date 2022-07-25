using SpiceSharp.Attributes;
using SpiceSharp.ParameterSets;
using System;
using System.Collections.Generic;
using System.Text;
using SpiceSharp.Components;
using SpiceSharp;

namespace SpiceSharpBSIM.Components.Semiconductors.BSIM.BSIM4Behaviors
{
    /// <summary>
    /// Parameters for a <see cref="BSIM4Model"/>
    /// </summary>
    [GeneratedParameters]
    public partial class ModelParameters : ParameterSet<ModelParameters>
    {
        [ParameterName("path"), ParameterInfo("Filename for used to log parameter checks")]
        public string CheckPath { get; set; } = "b4v8check.log";

        [ParameterName("cvchargemod"), ParameterInfo("Capacitance Charge model selector")]
        private GivenParameter<int> _cvchargeMod = new GivenParameter<int>(0);

        [ParameterName("capmod"), ParameterInfo("Capacitance model selector")]
        private GivenParameter<int> _capMod = new GivenParameter<int>();

        [ParameterName("diomod"), ParameterInfo("Diode IV model selector")]
        private GivenParameter<int> _dioMod = new GivenParameter<int>();

        [ParameterName("rdsmod"), ParameterInfo("Bias-dependent S/D resistance model selector")]
        private GivenParameter<int> _rdsMod = new GivenParameter<int>();

        [ParameterName("trnqsmod"), ParameterInfo("Transient NQS model selector")]
        private GivenParameter<int> _trnqsMod = new GivenParameter<int>();

        [ParameterName("acnqsmod"), ParameterInfo("AC NQS model selector")]
        private GivenParameter<int> _acnqsMod = new GivenParameter<int>();

        [ParameterName("mobmod"), ParameterInfo("Mobility model selector")]
        private GivenParameter<int> _mobMod = new GivenParameter<int>();

        [ParameterName("rbodymod"), ParameterInfo("Distributed body R model selector")]
        private GivenParameter<int> _rbodyMod = new GivenParameter<int>();

        [ParameterName("rgatemod"), ParameterInfo("Gate R model selector")]
        private GivenParameter<int> _rgateMod = new GivenParameter<int>();

        [ParameterName("permod"), ParameterInfo("Pd and Ps model selector")]
        private GivenParameter<int> _perMod = new GivenParameter<int>();

        [ParameterName("geomod"), ParameterInfo("Geometry dependent parasitics model selector")]
        private GivenParameter<int> _geoMod = new GivenParameter<int>(0);

        [ParameterName("rgeomod"), ParameterInfo("S/D resistance and contact model selector")]
        private GivenParameter<int> _rgeoMod = new GivenParameter<int>();

        [ParameterName("fnoimod"), ParameterInfo("Flicker noise model selector")]
        private GivenParameter<int> _fnoiMod = new GivenParameter<int>();

        [ParameterName("tnoimod"), ParameterInfo("Thermal noise model selector")]
        private GivenParameter<int> _tnoiMod = new GivenParameter<int>();

        [ParameterName("mtrlmod"), ParameterInfo("parameter for non-silicon substrate or metal gate selector")]
        private GivenParameter<int> _mtrlMod = new GivenParameter<int>();

        [ParameterName("mtrlcompatmod"), ParameterInfo("New Material Mod backward compatibility selector")]
        private GivenParameter<int> _mtrlCompatMod = new GivenParameter<int>();

        [ParameterName("igcmod"), ParameterInfo("Gate-to-channel Ig model selector")]
        private GivenParameter<int> _igcMod = new GivenParameter<int>();

        [ParameterName("igbmod"), ParameterInfo("Gate-to-body Ig model selector")]
        private GivenParameter<int> _igbMod = new GivenParameter<int>();

        [ParameterName("tempmod"), ParameterInfo("Temperature model selector")]
        private GivenParameter<int> _tempMod = new GivenParameter<int>();

        [ParameterName("gidlmod"), ParameterInfo("parameter for GIDL selector")]
        private GivenParameter<int> _gidlMod = new GivenParameter<int>(0);

        [ParameterName("paramchk"), ParameterInfo("Model parameter checking selector")]
        private GivenParameter<int> _paramChk = new GivenParameter<int>(1);

        [ParameterName("binunit"), ParameterInfo("Bin  unit  selector")]
        private GivenParameter<int> _binUnit = new GivenParameter<int>(1);

        [ParameterName("version"), ParameterInfo("parameter for model version")]
        private GivenParameter<string> _version = new GivenParameter<string>();

        [ParameterName("eot"), ParameterInfo("Equivalent gate oxide thickness in meters")]
        [Finite]
        private GivenParameter<double> _eot = new GivenParameter<double>(15.0e-10);

        [ParameterName("vddeot"), ParameterInfo("Voltage for extraction of Equivalent gate oxide thickness")]
        [Finite]
        private GivenParameter<double> _vddeot = new GivenParameter<double>();

        [ParameterName("tempeot"), ParameterInfo(" Temperature for extraction of EOT")]
        [Finite]
        private GivenParameter<double> _tempeot = new GivenParameter<double>(300.15);

        [ParameterName("leffeot"), ParameterInfo(" Effective length for extraction of EOT")]
        [Finite]
        private GivenParameter<double> _leffeot = new GivenParameter<double>(1);

        [ParameterName("weffeot"), ParameterInfo("Effective width for extraction of EOT")]
        [Finite]
        private GivenParameter<double> _weffeot = new GivenParameter<double>(10);

        [ParameterName("ados"), ParameterInfo("Charge centroid parameter")]
        [Finite]
        private GivenParameter<double> _ados = new GivenParameter<double>(1.0);

        [ParameterName("bdos"), ParameterInfo("Charge centroid parameter")]
        [Finite]
        private GivenParameter<double> _bdos = new GivenParameter<double>(1.0);

        [ParameterName("toxe"), ParameterInfo("Electrical gate oxide thickness in meters")]
        [Finite]
        private GivenParameter<double> _toxe = new GivenParameter<double>(30.0e-10);

        [ParameterName("toxp"), ParameterInfo("Physical gate oxide thickness in meters")]
        [Finite]
        private GivenParameter<double> _toxp = new GivenParameter<double>();

        [ParameterName("toxm"), ParameterInfo("Gate oxide thickness at which parameters are extracted")]
        [Finite]
        private GivenParameter<double> _toxm = new GivenParameter<double>();

        [ParameterName("toxref"), ParameterInfo("Target tox value")]
        [Finite]
        private GivenParameter<double> _toxref = new GivenParameter<double>(30.0e-10);

        [ParameterName("dtox"), ParameterInfo("Defined as (toxe - toxp) ")]
        [Finite]
        private GivenParameter<double> _dtox = new GivenParameter<double>(0.0);

        [ParameterName("epsrox"), ParameterInfo("Dielectric constant of the gate oxide relative to vacuum")]
        [Finite]
        private GivenParameter<double> _epsrox = new GivenParameter<double>(3.9);

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
        private GivenParameter<double> _nfactor = new GivenParameter<double>(1.0);

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

        [ParameterName("phig"), ParameterInfo("Work function of gate")]
        [Finite]
        private GivenParameter<double> _phig = new GivenParameter<double>(4.05);

        [ParameterName("epsrgate"), ParameterInfo("Dielectric constant of gate relative to vacuum")]
        [Finite]
        private GivenParameter<double> _epsrgate = new GivenParameter<double>(11.7);

        [ParameterName("easub"), ParameterInfo("Electron affinity of substrate")]
        [Finite]
        private GivenParameter<double> _easub = new GivenParameter<double>(4.05);

        [ParameterName("epsrsub"), ParameterInfo("Dielectric constant of substrate relative to vacuum")]
        [Finite]
        private GivenParameter<double> _epsrsub = new GivenParameter<double>(11.7);

        [ParameterName("ni0sub"), ParameterInfo("Intrinsic carrier concentration of substrate at 300.15K")]
        [Finite]
        private GivenParameter<double> _ni0sub = new GivenParameter<double>(1.45e10);

        [ParameterName("bg0sub"), ParameterInfo("Band-gap of substrate at T=0K")]
        [Finite]
        private GivenParameter<double> _bg0sub = new GivenParameter<double>(1.16);

        [ParameterName("tbgasub"), ParameterInfo("First parameter of band-gap change due to temperature")]
        [Finite]
        private GivenParameter<double> _tbgasub = new GivenParameter<double>(7.02e-4);

        [ParameterName("tbgbsub"), ParameterInfo("Second parameter of band-gap change due to temperature")]
        [Finite]
        private GivenParameter<double> _tbgbsub = new GivenParameter<double>(1108.0);

        [ParameterName("nsub"), ParameterInfo("Substrate doping concentration")]
        [Finite]
        private GivenParameter<double> _nsub = new GivenParameter<double>(6.0e16);

        [ParameterName("ndep"), ParameterInfo("Channel doping concentration at the depletion edge")]
        [Finite]
        public GivenParameter<double> Ndep
        {
            get => _ndep;
            set
            {
                Utility.Finite(value, nameof(Ndep));
                if (value > 1.0e20)
                    _ndep = value * 1e-6;
                else
                    _ndep = value;
            }
        }
        private GivenParameter<double> _ndep = new GivenParameter<double>(1.7e17);

        [ParameterName("nsd"), ParameterInfo("S/D doping concentration")]
        [Finite]
        public GivenParameter<double> Nsd
        {
            get => _nsd;
            set
            {
                Utility.Finite(value, nameof(Nsd));
                if (value > 1.000001e24)
                    _nsd = value * 1e-6;
                else
                    _nsd = value;
            }
        }
        private GivenParameter<double> _nsd = new GivenParameter<double>(1.0e20);

        [ParameterName("phin"), ParameterInfo("Adjusting parameter for surface potential due to non-uniform vertical doping")]
        [Finite]
        private GivenParameter<double> _phin = new GivenParameter<double>(0.0);

        [ParameterName("ngate"), ParameterInfo("Poly-gate doping concentration")]
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
        private GivenParameter<double> _ngate = new GivenParameter<double>(0);

        [ParameterName("gamma1"), ParameterInfo("Vth body coefficient")]
        [Finite]
        private GivenParameter<double> _gamma1 = new GivenParameter<double>(0.0);

        [ParameterName("gamma2"), ParameterInfo("Vth body coefficient")]
        [Finite]
        private GivenParameter<double> _gamma2 = new GivenParameter<double>(0.0);

        [ParameterName("vbx"), ParameterInfo("Vth transition body Voltage")]
        [Finite]
        private GivenParameter<double> _vbx = new GivenParameter<double>(0.0);

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

        [ParameterName("dvtp0"), ParameterInfo("First parameter for Vth shift due to pocket")]
        [Finite]
        private GivenParameter<double> _dvtp0 = new GivenParameter<double>(0.0);

        [ParameterName("dvtp1"), ParameterInfo("Second parameter for Vth shift due to pocket")]
        [Finite]
        private GivenParameter<double> _dvtp1 = new GivenParameter<double>(0.0);

        [ParameterName("dvtp2"), ParameterInfo("3rd parameter for Vth shift due to pocket")]
        [Finite]
        private GivenParameter<double> _dvtp2 = new GivenParameter<double>(0.0);

        [ParameterName("dvtp3"), ParameterInfo("4th parameter for Vth shift due to pocket")]
        [Finite]
        private GivenParameter<double> _dvtp3 = new GivenParameter<double>(0.0);

        [ParameterName("dvtp4"), ParameterInfo("5th parameter for Vth shift due to pocket")]
        [Finite]
        private GivenParameter<double> _dvtp4 = new GivenParameter<double>(0.0);

        [ParameterName("dvtp5"), ParameterInfo("6th parameter for Vth shift due to pocket")]
        [Finite]
        private GivenParameter<double> _dvtp5 = new GivenParameter<double>(0.0);

        [ParameterName("lpe0"), ParameterInfo("Equivalent length of pocket region at zero bias")]
        [Finite]
        private GivenParameter<double> _lpe0 = new GivenParameter<double>(1.74e-7);

        [ParameterName("lpeb"), ParameterInfo("Equivalent length of pocket region accounting for body bias")]
        [Finite]
        private GivenParameter<double> _lpeb = new GivenParameter<double>(0.0);

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

        [ParameterName("vth0"), ParameterName("vtho"), ParameterInfo("Threshold voltage")]
        [Finite]
        private GivenParameter<double> _vth0 = new GivenParameter<double>();

        [ParameterName("ua"), ParameterInfo("Linear gate dependence of mobility")]
        [Finite]
        private GivenParameter<double> _ua = new GivenParameter<double>();

        [ParameterName("ua1"), ParameterInfo("Temperature coefficient of ua")]
        [Finite]
        private GivenParameter<double> _ua1 = new GivenParameter<double>(1.0e-9);

        [ParameterName("ub"), ParameterInfo("Quadratic gate dependence of mobility")]
        [Finite]
        private GivenParameter<double> _ub = new GivenParameter<double>(1.0e-19);

        [ParameterName("ub1"), ParameterInfo("Temperature coefficient of ub")]
        [Finite]
        private GivenParameter<double> _ub1 = new GivenParameter<double>(-1.0e-18);

        [ParameterName("uc"), ParameterInfo("Body-bias dependence of mobility")]
        [Finite]
        private GivenParameter<double> _uc = new GivenParameter<double>();

        [ParameterName("uc1"), ParameterInfo("Temperature coefficient of uc")]
        [Finite]
        private GivenParameter<double> _uc1 = new GivenParameter<double>();

        [ParameterName("ud"), ParameterInfo("Coulomb scattering factor of mobility")]
        [Finite]
        private GivenParameter<double> _ud = new GivenParameter<double>(0.0);

        [ParameterName("ud1"), ParameterInfo("Temperature coefficient of ud")]
        [Finite]
        private GivenParameter<double> _ud1 = new GivenParameter<double>(0.0);

        [ParameterName("up"), ParameterInfo("Channel length linear factor of mobility")]
        [Finite]
        private GivenParameter<double> _up = new GivenParameter<double>(0.0);

        [ParameterName("lp"), ParameterInfo("Channel length exponential factor of mobility")]
        [Finite]
        private GivenParameter<double> _lp = new GivenParameter<double>(1.0e-8);

        [ParameterName("u0"), ParameterInfo("Low-field mobility at Tnom")]
        [Finite]
        private GivenParameter<double> _u0 = new GivenParameter<double>();

        [ParameterName("eu"), ParameterInfo("Mobility exponent")]
        [Finite]
        private GivenParameter<double> _eu = new GivenParameter<double>();

        [ParameterName("ucs"), ParameterInfo("Colombic scattering exponent")]
        [Finite]
        private GivenParameter<double> _ucs = new GivenParameter<double>();

        [ParameterName("ute"), ParameterInfo("Temperature coefficient of mobility")]
        [Finite]
        private GivenParameter<double> _ute = new GivenParameter<double>(-1.5);

        [ParameterName("ucste"), ParameterInfo("Temperature coefficient of colombic mobility")]
        [Finite]
        private GivenParameter<double> _ucste = new GivenParameter<double>(-4.775e-3);

        [ParameterName("voff"), ParameterInfo("Threshold voltage offset")]
        [Finite]
        private GivenParameter<double> _voff = new GivenParameter<double>(-0.08);

        [ParameterName("minv"), ParameterInfo("Fitting parameter for moderate inversion in Vgsteff")]
        [Finite]
        private GivenParameter<double> _minv = new GivenParameter<double>(0.0);

        [ParameterName("minvcv"), ParameterInfo("Fitting parameter for moderate inversion in Vgsteffcv")]
        [Finite]
        private GivenParameter<double> _minvcv = new GivenParameter<double>(0.0);

        [ParameterName("voffl"), ParameterInfo("Length dependence parameter for Vth offset")]
        [Finite]
        private GivenParameter<double> _voffl = new GivenParameter<double>(0.0);

        [ParameterName("voffcvl"), ParameterInfo("Length dependence parameter for Vth offset in CV")]
        [Finite]
        private GivenParameter<double> _voffcvl = new GivenParameter<double>(0.0);

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

        [ParameterName("delta"), ParameterInfo("Effective Vds parameter")]
        [Finite]
        private GivenParameter<double> _delta = new GivenParameter<double>(0.01);

        [ParameterName("rsh"), ParameterInfo("Source-drain sheet resistance")]
        [Finite]
        private GivenParameter<double> _sheetResistance = new GivenParameter<double>(0.0);

        [ParameterName("rdsw"), ParameterInfo("Source-drain resistance per width")]
        [Finite]
        private GivenParameter<double> _rdsw = new GivenParameter<double>(200.0);

        [ParameterName("rdswmin"), ParameterInfo("Source-drain resistance per width at high Vg")]
        [Finite]
        private GivenParameter<double> _rdswmin = new GivenParameter<double>(0.0);

        [ParameterName("rsw"), ParameterInfo("Source resistance per width")]
        [Finite]
        private GivenParameter<double> _rsw = new GivenParameter<double>(100.0);

        [ParameterName("rdw"), ParameterInfo("Drain resistance per width")]
        [Finite]
        private GivenParameter<double> _rdw = new GivenParameter<double>(100.0);

        [ParameterName("rdwmin"), ParameterInfo("Drain resistance per width at high Vg")]
        [Finite]
        private GivenParameter<double> _rdwmin = new GivenParameter<double>(0.0);

        [ParameterName("rswmin"), ParameterInfo("Source resistance per width at high Vg")]
        [Finite]
        private GivenParameter<double> _rswmin = new GivenParameter<double>(0.0);

        [ParameterName("prwg"), ParameterInfo("Gate-bias effect on parasitic resistance ")]
        [Finite]
        private GivenParameter<double> _prwg = new GivenParameter<double>(1.0);

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
        private GivenParameter<double> _pdibl1 = new GivenParameter<double>(0.39);

        [ParameterName("pdiblc2"), ParameterInfo("Drain-induced barrier lowering coefficient")]
        [Finite]
        private GivenParameter<double> _pdibl2 = new GivenParameter<double>(0.0086);

        [ParameterName("pdiblcb"), ParameterInfo("Body-effect on drain-induced barrier lowering")]
        [Finite]
        private GivenParameter<double> _pdiblb = new GivenParameter<double>(0.0);

        [ParameterName("fprout"), ParameterInfo("Rout degradation coefficient for pocket devices")]
        [Finite]
        private GivenParameter<double> _fprout = new GivenParameter<double>(0.0);

        [ParameterName("pdits"), ParameterInfo("Coefficient for drain-induced Vth shifts")]
        [Finite]
        private GivenParameter<double> _pdits = new GivenParameter<double>(0.0);

        [ParameterName("pditsl"), ParameterInfo("Length dependence of drain-induced Vth shifts")]
        [Finite]
        private GivenParameter<double> _pditsl = new GivenParameter<double>(0.0);

        [ParameterName("pditsd"), ParameterInfo("Vds dependence of drain-induced Vth shifts")]
        [Finite]
        private GivenParameter<double> _pditsd = new GivenParameter<double>(0.0);

        [ParameterName("pscbe1"), ParameterInfo("Substrate current body-effect coefficient")]
        [Finite]
        private GivenParameter<double> _pscbe1 = new GivenParameter<double>(4.24e8);

        [ParameterName("pscbe2"), ParameterInfo("Substrate current body-effect coefficient")]
        [Finite]
        private GivenParameter<double> _pscbe2 = new GivenParameter<double>(1.0e-5);

        [ParameterName("pvag"), ParameterInfo("Gate dependence of output resistance parameter")]
        [Finite]
        private GivenParameter<double> _pvag = new GivenParameter<double>(0.0);

        [ParameterName("jss"), ParameterInfo("Bottom source junction reverse saturation current density")]
        [Finite]
        private GivenParameter<double> _sjctSatCurDensity = new GivenParameter<double>(1.0E-4);

        [ParameterName("jsws"), ParameterInfo("Isolation edge sidewall source junction reverse saturation current density")]
        [Finite]
        private GivenParameter<double> _sjctSidewallSatCurDensity = new GivenParameter<double>(0.0);

        [ParameterName("jswgs"), ParameterInfo("Gate edge source junction reverse saturation current density")]
        [Finite]
        private GivenParameter<double> _sjctGateSidewallSatCurDensity = new GivenParameter<double>(0.0);

        [ParameterName("pbs"), ParameterInfo("Source junction built-in potential")]
        [Finite]
        private GivenParameter<double> _sbulkJctPotential = new GivenParameter<double>(1.0);

        [ParameterName("njs"), ParameterInfo("Source junction emission coefficient")]
        [Finite]
        private GivenParameter<double> _sjctEmissionCoeff = new GivenParameter<double>(1.0);

        [ParameterName("xtis"), ParameterInfo("Source junction current temperature exponent")]
        [Finite]
        private GivenParameter<double> _sjctTempExponent = new GivenParameter<double>(3.0);

        [ParameterName("mjs"), ParameterInfo("Source bottom junction capacitance grading coefficient")]
        [Finite]
        private GivenParameter<double> _sbulkJctBotGradingCoeff = new GivenParameter<double>(0.5);

        [ParameterName("pbsws"), ParameterInfo("Source sidewall junction capacitance built in potential")]
        [Finite]
        private GivenParameter<double> _ssidewallJctPotential = new GivenParameter<double>(1.0);

        [ParameterName("mjsws"), ParameterInfo("Source sidewall junction capacitance grading coefficient")]
        [Finite]
        private GivenParameter<double> _sbulkJctSideGradingCoeff = new GivenParameter<double>(0.33);

        [ParameterName("pbswgs"), ParameterInfo("Source (gate side) sidewall junction capacitance built in potential")]
        [Finite]
        private GivenParameter<double> _sGatesidewallJctPotential = new GivenParameter<double>();

        [ParameterName("mjswgs"), ParameterInfo("Source (gate side) sidewall junction capacitance grading coefficient")]
        [Finite]
        private GivenParameter<double> _sbulkJctGateSideGradingCoeff = new GivenParameter<double>();

        [ParameterName("cjs"), ParameterInfo("Source bottom junction capacitance per unit area")]
        [Finite]
        private GivenParameter<double> _sunitAreaJctCap = new GivenParameter<double>(5.0E-4);

        [ParameterName("cjsws"), ParameterInfo("Source sidewall junction capacitance per unit periphery")]
        [Finite]
        private GivenParameter<double> _sunitLengthSidewallJctCap = new GivenParameter<double>(5.0E-10);

        [ParameterName("cjswgs"), ParameterInfo("Source (gate side) sidewall junction capacitance per unit width")]
        [Finite]
        private GivenParameter<double> _sunitLengthGateSidewallJctCap = new GivenParameter<double>();

        [ParameterName("jsd"), ParameterInfo("Bottom drain junction reverse saturation current density")]
        [Finite]
        private GivenParameter<double> _djctSatCurDensity = new GivenParameter<double>();

        [ParameterName("jswd"), ParameterInfo("Isolation edge sidewall drain junction reverse saturation current density")]
        [Finite]
        private GivenParameter<double> _djctSidewallSatCurDensity = new GivenParameter<double>();

        [ParameterName("jswgd"), ParameterInfo("Gate edge drain junction reverse saturation current density")]
        [Finite]
        private GivenParameter<double> _djctGateSidewallSatCurDensity = new GivenParameter<double>();

        [ParameterName("pbd"), ParameterInfo("Drain junction built-in potential")]
        [Finite]
        private GivenParameter<double> _dbulkJctPotential = new GivenParameter<double>();

        [ParameterName("njd"), ParameterInfo("Drain junction emission coefficient")]
        [Finite]
        private GivenParameter<double> _djctEmissionCoeff = new GivenParameter<double>();

        [ParameterName("xtid"), ParameterInfo("Drainjunction current temperature exponent")]
        [Finite]
        private GivenParameter<double> _djctTempExponent = new GivenParameter<double>();

        [ParameterName("mjd"), ParameterInfo("Drain bottom junction capacitance grading coefficient")]
        [Finite]
        private GivenParameter<double> _dbulkJctBotGradingCoeff = new GivenParameter<double>();

        [ParameterName("pbswd"), ParameterInfo("Drain sidewall junction capacitance built in potential")]
        [Finite]
        private GivenParameter<double> _dsidewallJctPotential = new GivenParameter<double>();

        [ParameterName("mjswd"), ParameterInfo("Drain sidewall junction capacitance grading coefficient")]
        [Finite]
        private GivenParameter<double> _dbulkJctSideGradingCoeff = new GivenParameter<double>();

        [ParameterName("pbswgd"), ParameterInfo("Drain (gate side) sidewall junction capacitance built in potential")]
        [Finite]
        private GivenParameter<double> _dGatesidewallJctPotential = new GivenParameter<double>();

        [ParameterName("mjswgd"), ParameterInfo("Drain (gate side) sidewall junction capacitance grading coefficient")]
        [Finite]
        private GivenParameter<double> _dbulkJctGateSideGradingCoeff = new GivenParameter<double>();

        [ParameterName("cjd"), ParameterInfo("Drain bottom junction capacitance per unit area")]
        [Finite]
        private GivenParameter<double> _dunitAreaJctCap = new GivenParameter<double>();

        [ParameterName("cjswd"), ParameterInfo("Drain sidewall junction capacitance per unit periphery")]
        [Finite]
        private GivenParameter<double> _dunitLengthSidewallJctCap = new GivenParameter<double>();

        [ParameterName("cjswgd"), ParameterInfo("Drain (gate side) sidewall junction capacitance per unit width")]
        [Finite]
        private GivenParameter<double> _dunitLengthGateSidewallJctCap = new GivenParameter<double>();

        [ParameterName("vfbcv"), ParameterInfo("Flat Band Voltage parameter for capmod=0 only")]
        [Finite]
        private GivenParameter<double> _vfbcv = new GivenParameter<double>(-1.0);

        [ParameterName("vfb"), ParameterInfo("Flat Band Voltage")]
        [Finite]
        private GivenParameter<double> _vfb = new GivenParameter<double>(-1.0);

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

        [ParameterName("dmcg"), ParameterInfo("Distance of Mid-Contact to Gate edge")]
        [Finite]
        private GivenParameter<double> _dmcg = new GivenParameter<double>(0.0);

        [ParameterName("dmci"), ParameterInfo("Distance of Mid-Contact to Isolation")]
        [Finite]
        private GivenParameter<double> _dmci = new GivenParameter<double>();

        [ParameterName("dmdg"), ParameterInfo("Distance of Mid-Diffusion to Gate edge")]
        [Finite]
        private GivenParameter<double> _dmdg = new GivenParameter<double>(0.0);

        [ParameterName("dmcgt"), ParameterInfo("Distance of Mid-Contact to Gate edge in Test structures")]
        [Finite]
        private GivenParameter<double> _dmcgt = new GivenParameter<double>(0.0);

        [ParameterName("xgw"), ParameterInfo("Distance from gate contact center to device edge")]
        [Finite]
        private GivenParameter<double> _xgw = new GivenParameter<double>(0.0);

        [ParameterName("xgl"), ParameterInfo("Variation in Ldrawn")]
        [Finite]
        private GivenParameter<double> _xgl = new GivenParameter<double>(0.0);

        [ParameterName("rshg"), ParameterInfo("Gate sheet resistance")]
        [Finite]
        private GivenParameter<double> _rshg = new GivenParameter<double>(0.1);

        [ParameterName("ngcon"), ParameterInfo("Number of gate contacts")]
        [Finite]
        private GivenParameter<double> _ngcon = new GivenParameter<double>(1.0);

        [ParameterName("xrcrg1"), ParameterInfo("First fitting parameter the bias-dependent Rg")]
        [Finite]
        private GivenParameter<double> _xrcrg1 = new GivenParameter<double>(12.0);

        [ParameterName("xrcrg2"), ParameterInfo("Second fitting parameter the bias-dependent Rg")]
        [Finite]
        private GivenParameter<double> _xrcrg2 = new GivenParameter<double>(1.0);

        [ParameterName("lambda"), ParameterInfo(" Velocity overshoot parameter")]
        [Finite]
        private GivenParameter<double> _lambda = new GivenParameter<double>(0.0);

        [ParameterName("vtl"), ParameterInfo(" thermal velocity")]
        [Finite]
        private GivenParameter<double> _vtl = new GivenParameter<double>(2.0e5);

        [ParameterName("lc"), ParameterInfo(" back scattering parameter")]
        [Finite]
        private GivenParameter<double> _lc = new GivenParameter<double>(5.0e-9);

        [ParameterName("xn"), ParameterInfo(" back scattering parameter")]
        [Finite]
        private GivenParameter<double> _xn = new GivenParameter<double>(3.0);

        [ParameterName("vfbsdoff"), ParameterInfo("S/D flatband voltage offset")]
        [Finite]
        private GivenParameter<double> _vfbsdoff = new GivenParameter<double>(0.0);

        [ParameterName("tvfbsdoff"), ParameterInfo("Temperature parameter for vfbsdoff")]
        [Finite]
        private GivenParameter<double> _tvfbsdoff = new GivenParameter<double>(0.0);

        [ParameterName("tvoff"), ParameterInfo("Temperature parameter for voff")]
        [Finite]
        private GivenParameter<double> _tvoff = new GivenParameter<double>(0.0);

        [ParameterName("tnfactor"), ParameterInfo("Temperature parameter for nfactor")]
        [Finite]
        private GivenParameter<double> _tnfactor = new GivenParameter<double>(0.0);

        [ParameterName("teta0"), ParameterInfo("Temperature parameter for eta0")]
        [Finite]
        private GivenParameter<double> _teta0 = new GivenParameter<double>(0.0);

        [ParameterName("tvoffcv"), ParameterInfo("Temperature parameter for tvoffcv")]
        [Finite]
        private GivenParameter<double> _tvoffcv = new GivenParameter<double>(0.0);

        [ParameterName("lintnoi"), ParameterInfo("lint offset for noise calculation")]
        [Finite]
        private GivenParameter<double> _lintnoi = new GivenParameter<double>(0.0);

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

        [ParameterName("ckappas"), ParameterInfo("S/G overlap C-V parameter ")]
        [Finite]
        private GivenParameter<double> _ckappas = new GivenParameter<double>(0.6);

        [ParameterName("ckappad"), ParameterInfo("D/G overlap C-V parameter")]
        [Finite]
        private GivenParameter<double> _ckappad = new GivenParameter<double>();

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

        [ParameterName("xw"), ParameterInfo("W offset for channel width due to mask/etch effect")]
        [Finite]
        private GivenParameter<double> _xw = new GivenParameter<double>(0.0);

        [ParameterName("xl"), ParameterInfo("L offset for channel length due to mask/etch effect")]
        [Finite]
        private GivenParameter<double> _xl = new GivenParameter<double>(0.0);

        [ParameterName("dlcig"), ParameterInfo("Delta L for Ig model")]
        [Finite]
        private GivenParameter<double> _dlcig = new GivenParameter<double>();

        [ParameterName("dlcigd"), ParameterInfo("Delta L for Ig model drain side")]
        [Finite]
        private GivenParameter<double> _dlcigd = new GivenParameter<double>();

        [ParameterName("dwj"), ParameterInfo("Delta W for S/D junctions")]
        [Finite]
        private GivenParameter<double> _dwj = new GivenParameter<double>();

        [ParameterName("alpha0"), ParameterInfo("substrate current model parameter")]
        [Finite]
        private GivenParameter<double> _alpha0 = new GivenParameter<double>(0.0);

        [ParameterName("alpha1"), ParameterInfo("substrate current model parameter")]
        [Finite]
        private GivenParameter<double> _alpha1 = new GivenParameter<double>(0.0);

        [ParameterName("beta0"), ParameterInfo("substrate current model parameter")]
        [Finite]
        private GivenParameter<double> _beta0 = new GivenParameter<double>(0.0);

        [ParameterName("agidl"), ParameterInfo("Pre-exponential constant for GIDL")]
        [Finite]
        private GivenParameter<double> _agidl = new GivenParameter<double>(0.0);

        [ParameterName("bgidl"), ParameterInfo("Exponential constant for GIDL")]
        [Finite]
        private GivenParameter<double> _bgidl = new GivenParameter<double>(2.3e9);

        [ParameterName("cgidl"), ParameterInfo("Parameter for body-bias dependence of GIDL")]
        [Finite]
        private GivenParameter<double> _cgidl = new GivenParameter<double>(0.5);

        [ParameterName("rgidl"), ParameterInfo("GIDL vg parameter")]
        [Finite]
        private GivenParameter<double> _rgidl = new GivenParameter<double>(1.0);

        [ParameterName("kgidl"), ParameterInfo("GIDL vb parameter")]
        [Finite]
        private GivenParameter<double> _kgidl = new GivenParameter<double>(0.0);

        [ParameterName("fgidl"), ParameterInfo("GIDL vb parameter")]
        [Finite]
        private GivenParameter<double> _fgidl = new GivenParameter<double>(1.0);

        [ParameterName("egidl"), ParameterInfo("Fitting parameter for Bandbending")]
        [Finite]
        private GivenParameter<double> _egidl = new GivenParameter<double>(0.8);

        [ParameterName("agisl"), ParameterInfo("Pre-exponential constant for GISL")]
        [Finite]
        private GivenParameter<double> _agisl = new GivenParameter<double>();

        [ParameterName("bgisl"), ParameterInfo("Exponential constant for GISL")]
        [Finite]
        private GivenParameter<double> _bgisl = new GivenParameter<double>();

        [ParameterName("cgisl"), ParameterInfo("Parameter for body-bias dependence of GISL")]
        [Finite]
        private GivenParameter<double> _cgisl = new GivenParameter<double>();

        [ParameterName("rgisl"), ParameterInfo("GISL vg parameter")]
        [Finite]
        private GivenParameter<double> _rgisl = new GivenParameter<double>();

        [ParameterName("kgisl"), ParameterInfo("GISL vb parameter")]
        [Finite]
        private GivenParameter<double> _kgisl = new GivenParameter<double>();

        [ParameterName("fgisl"), ParameterInfo("GISL vb parameter")]
        [Finite]
        private GivenParameter<double> _fgisl = new GivenParameter<double>();

        [ParameterName("egisl"), ParameterInfo("Fitting parameter for Bandbending")]
        [Finite]
        private GivenParameter<double> _egisl = new GivenParameter<double>();

        [ParameterName("aigc"), ParameterInfo("Parameter for Igc")]
        [Finite]
        private GivenParameter<double> _aigc = new GivenParameter<double>();

        [ParameterName("bigc"), ParameterInfo("Parameter for Igc")]
        [Finite]
        private GivenParameter<double> _bigc = new GivenParameter<double>();

        [ParameterName("cigc"), ParameterInfo("Parameter for Igc")]
        [Finite]
        private GivenParameter<double> _cigc = new GivenParameter<double>();

        [ParameterName("aigsd"), ParameterInfo("Parameter for Igs,d")]
        [Finite]
        private GivenParameter<double> _aigsd = new GivenParameter<double>();

        [ParameterName("bigsd"), ParameterInfo("Parameter for Igs,d")]
        [Finite]
        private GivenParameter<double> _bigsd = new GivenParameter<double>();

        [ParameterName("cigsd"), ParameterInfo("Parameter for Igs,d")]
        [Finite]
        private GivenParameter<double> _cigsd = new GivenParameter<double>();

        [ParameterName("aigs"), ParameterInfo("Parameter for Igs")]
        [Finite]
        private GivenParameter<double> _aigs = new GivenParameter<double>();

        [ParameterName("bigs"), ParameterInfo("Parameter for Igs")]
        [Finite]
        private GivenParameter<double> _bigs = new GivenParameter<double>();

        [ParameterName("cigs"), ParameterInfo("Parameter for Igs")]
        [Finite]
        private GivenParameter<double> _cigs = new GivenParameter<double>();

        [ParameterName("aigd"), ParameterInfo("Parameter for Igd")]
        [Finite]
        private GivenParameter<double> _aigd = new GivenParameter<double>();

        [ParameterName("bigd"), ParameterInfo("Parameter for Igd")]
        [Finite]
        private GivenParameter<double> _bigd = new GivenParameter<double>();

        [ParameterName("cigd"), ParameterInfo("Parameter for Igd")]
        [Finite]
        private GivenParameter<double> _cigd = new GivenParameter<double>();

        [ParameterName("aigbacc"), ParameterInfo("Parameter for Igb")]
        [Finite]
        private GivenParameter<double> _aigbacc = new GivenParameter<double>(1.36e-2);

        [ParameterName("bigbacc"), ParameterInfo("Parameter for Igb")]
        [Finite]
        private GivenParameter<double> _bigbacc = new GivenParameter<double>(1.71e-3);

        [ParameterName("cigbacc"), ParameterInfo("Parameter for Igb")]
        [Finite]
        private GivenParameter<double> _cigbacc = new GivenParameter<double>(0.075);

        [ParameterName("aigbinv"), ParameterInfo("Parameter for Igb")]
        [Finite]
        private GivenParameter<double> _aigbinv = new GivenParameter<double>(1.11e-2);

        [ParameterName("bigbinv"), ParameterInfo("Parameter for Igb")]
        [Finite]
        private GivenParameter<double> _bigbinv = new GivenParameter<double>(9.49e-4);

        [ParameterName("cigbinv"), ParameterInfo("Parameter for Igb")]
        [Finite]
        private GivenParameter<double> _cigbinv = new GivenParameter<double>(0.006);

        [ParameterName("nigc"), ParameterInfo("Parameter for Igc slope")]
        [Finite]
        private GivenParameter<double> _nigc = new GivenParameter<double>(1.0);

        [ParameterName("nigbinv"), ParameterInfo("Parameter for Igbinv slope")]
        [Finite]
        private GivenParameter<double> _nigbinv = new GivenParameter<double>(3.0);

        [ParameterName("nigbacc"), ParameterInfo("Parameter for Igbacc slope")]
        [Finite]
        private GivenParameter<double> _nigbacc = new GivenParameter<double>(1.0);

        [ParameterName("ntox"), ParameterInfo("Exponent for Tox ratio")]
        [Finite]
        private GivenParameter<double> _ntox = new GivenParameter<double>(1.0);

        [ParameterName("eigbinv"), ParameterInfo("Parameter for the Si bandgap for Igbinv")]
        [Finite]
        private GivenParameter<double> _eigbinv = new GivenParameter<double>(1.1);

        [ParameterName("pigcd"), ParameterInfo("Parameter for Igc partition")]
        [Finite]
        private GivenParameter<double> _pigcd = new GivenParameter<double>(1.0);

        [ParameterName("poxedge"), ParameterInfo("Factor for the gate edge Tox")]
        [Finite]
        private GivenParameter<double> _poxedge = new GivenParameter<double>(1.0);

        [ParameterName("ijthdfwd"), ParameterInfo("Forward drain diode forward limiting current")]
        [Finite]
        private GivenParameter<double> _ijthdfwd = new GivenParameter<double>();

        [ParameterName("ijthsfwd"), ParameterInfo("Forward source diode forward limiting current")]
        [Finite]
        private GivenParameter<double> _ijthsfwd = new GivenParameter<double>(0.1);

        [ParameterName("ijthdrev"), ParameterInfo("Reverse drain diode forward limiting current")]
        [Finite]
        private GivenParameter<double> _ijthdrev = new GivenParameter<double>();

        [ParameterName("ijthsrev"), ParameterInfo("Reverse source diode forward limiting current")]
        [Finite]
        private GivenParameter<double> _ijthsrev = new GivenParameter<double>(0.1);

        [ParameterName("xjbvd"), ParameterInfo("Fitting parameter for drain diode breakdown current")]
        [Finite]
        private GivenParameter<double> _xjbvd = new GivenParameter<double>();

        [ParameterName("xjbvs"), ParameterInfo("Fitting parameter for source diode breakdown current")]
        [Finite]
        private GivenParameter<double> _xjbvs = new GivenParameter<double>(1.0);

        [ParameterName("bvd"), ParameterInfo("Drain diode breakdown voltage")]
        [Finite]
        private GivenParameter<double> _bvd = new GivenParameter<double>();

        [ParameterName("bvs"), ParameterInfo("Source diode breakdown voltage")]
        [Finite]
        private GivenParameter<double> _bvs = new GivenParameter<double>(10.0);

        [ParameterName("jtss"), ParameterInfo("Source bottom trap-assisted saturation current density")]
        [Finite]
        private GivenParameter<double> _jtss = new GivenParameter<double>(0.0);

        [ParameterName("jtsd"), ParameterInfo("Drain bottom trap-assisted saturation current density")]
        [Finite]
        private GivenParameter<double> _jtsd = new GivenParameter<double>();

        [ParameterName("jtssws"), ParameterInfo("Source STI sidewall trap-assisted saturation current density")]
        [Finite]
        private GivenParameter<double> _jtssws = new GivenParameter<double>(0.0);

        [ParameterName("jtsswd"), ParameterInfo("Drain STI sidewall trap-assisted saturation current density")]
        [Finite]
        private GivenParameter<double> _jtsswd = new GivenParameter<double>();

        [ParameterName("jtsswgs"), ParameterInfo("Source gate-edge sidewall trap-assisted saturation current density")]
        [Finite]
        private GivenParameter<double> _jtsswgs = new GivenParameter<double>(0.0);

        [ParameterName("jtsswgd"), ParameterInfo("Drain gate-edge sidewall trap-assisted saturation current density")]
        [Finite]
        private GivenParameter<double> _jtsswgd = new GivenParameter<double>();

        [ParameterName("jtweff"), ParameterInfo("TAT current width dependence")]
        [Finite]
        private GivenParameter<double> _jtweff = new GivenParameter<double>(0.0);

        [ParameterName("njts"), ParameterInfo("Non-ideality factor for bottom junction")]
        [Finite]
        private GivenParameter<double> _njts = new GivenParameter<double>(20.0);

        [ParameterName("njtssw"), ParameterInfo("Non-ideality factor for STI sidewall junction")]
        [Finite]
        private GivenParameter<double> _njtssw = new GivenParameter<double>(20.0);

        [ParameterName("njtsswg"), ParameterInfo("Non-ideality factor for gate-edge sidewall junction")]
        [Finite]
        private GivenParameter<double> _njtsswg = new GivenParameter<double>(20.0);

        [ParameterName("njtsd"), ParameterInfo("Non-ideality factor for bottom junction drain side")]
        [Finite]
        private GivenParameter<double> _njtsd = new GivenParameter<double>();

        [ParameterName("njtsswd"), ParameterInfo("Non-ideality factor for STI sidewall junction drain side")]
        [Finite]
        private GivenParameter<double> _njtsswd = new GivenParameter<double>();

        [ParameterName("njtsswgd"), ParameterInfo("Non-ideality factor for gate-edge sidewall junction drain side")]
        [Finite]
        private GivenParameter<double> _njtsswgd = new GivenParameter<double>();

        [ParameterName("xtss"), ParameterInfo("Power dependence of JTSS on temperature")]
        [Finite]
        private GivenParameter<double> _xtss = new GivenParameter<double>(0.02);

        [ParameterName("xtsd"), ParameterInfo("Power dependence of JTSD on temperature")]
        [Finite]
        private GivenParameter<double> _xtsd = new GivenParameter<double>();

        [ParameterName("xtssws"), ParameterInfo("Power dependence of JTSSWS on temperature")]
        [Finite]
        private GivenParameter<double> _xtssws = new GivenParameter<double>(0.02);

        [ParameterName("xtsswd"), ParameterInfo("Power dependence of JTSSWD on temperature")]
        [Finite]
        private GivenParameter<double> _xtsswd = new GivenParameter<double>();

        [ParameterName("xtsswgs"), ParameterInfo("Power dependence of JTSSWGS on temperature")]
        [Finite]
        private GivenParameter<double> _xtsswgs = new GivenParameter<double>(0.02);

        [ParameterName("xtsswgd"), ParameterInfo("Power dependence of JTSSWGD on temperature")]
        [Finite]
        private GivenParameter<double> _xtsswgd = new GivenParameter<double>();

        [ParameterName("tnjts"), ParameterInfo("Temperature coefficient for NJTS")]
        [Finite]
        private GivenParameter<double> _tnjts = new GivenParameter<double>(0.0);

        [ParameterName("tnjtssw"), ParameterInfo("Temperature coefficient for NJTSSW")]
        [Finite]
        private GivenParameter<double> _tnjtssw = new GivenParameter<double>(0.0);

        [ParameterName("tnjtsswg"), ParameterInfo("Temperature coefficient for NJTSSWG")]
        [Finite]
        private GivenParameter<double> _tnjtsswg = new GivenParameter<double>(0.0);

        [ParameterName("tnjtsd"), ParameterInfo("Temperature coefficient for NJTSD")]
        [Finite]
        private GivenParameter<double> _tnjtsd = new GivenParameter<double>();

        [ParameterName("tnjtsswd"), ParameterInfo("Temperature coefficient for NJTSSWD")]
        [Finite]
        private GivenParameter<double> _tnjtsswd = new GivenParameter<double>();

        [ParameterName("tnjtsswgd"), ParameterInfo("Temperature coefficient for NJTSSWGD")]
        [Finite]
        private GivenParameter<double> _tnjtsswgd = new GivenParameter<double>();

        [ParameterName("vtss"), ParameterInfo("Source bottom trap-assisted voltage dependent parameter")]
        [Finite]
        private GivenParameter<double> _vtss = new GivenParameter<double>(10.0);

        [ParameterName("vtsd"), ParameterInfo("Drain bottom trap-assisted voltage dependent parameter")]
        [Finite]
        private GivenParameter<double> _vtsd = new GivenParameter<double>();

        [ParameterName("vtssws"), ParameterInfo("Source STI sidewall trap-assisted voltage dependent parameter")]
        [Finite]
        private GivenParameter<double> _vtssws = new GivenParameter<double>(10.0);

        [ParameterName("vtsswd"), ParameterInfo("Drain STI sidewall trap-assisted voltage dependent parameter")]
        [Finite]
        private GivenParameter<double> _vtsswd = new GivenParameter<double>();

        [ParameterName("vtsswgs"), ParameterInfo("Source gate-edge sidewall trap-assisted voltage dependent parameter")]
        [Finite]
        private GivenParameter<double> _vtsswgs = new GivenParameter<double>(10.0);

        [ParameterName("vtsswgd"), ParameterInfo("Drain gate-edge sidewall trap-assisted voltage dependent parameter")]
        [Finite]
        private GivenParameter<double> _vtsswgd = new GivenParameter<double>();

        [ParameterName("gbmin"), ParameterInfo("Minimum body conductance")]
        [Finite]
        private GivenParameter<double> _gbmin = new GivenParameter<double>(1.0e-12);

        [ParameterName("rbdb"), ParameterInfo("Resistance between bNode and dbNode")]
        [Finite]
        private GivenParameter<double> _rbdb = new GivenParameter<double>(50.0);

        [ParameterName("rbpb"), ParameterInfo("Resistance between bNodePrime and bNode")]
        [Finite]
        private GivenParameter<double> _rbpb = new GivenParameter<double>(50.0);

        [ParameterName("rbsb"), ParameterInfo("Resistance between bNode and sbNode")]
        [Finite]
        private GivenParameter<double> _rbsb = new GivenParameter<double>(50.0);

        [ParameterName("rbps"), ParameterInfo("Resistance between bNodePrime and sbNode")]
        [Finite]
        private GivenParameter<double> _rbps = new GivenParameter<double>(50.0);

        [ParameterName("rbpd"), ParameterInfo("Resistance between bNodePrime and bNode")]
        [Finite]
        private GivenParameter<double> _rbpd = new GivenParameter<double>(50.0);

        [ParameterName("rbps0"), ParameterInfo("Body resistance RBPS scaling")]
        [Finite]
        private GivenParameter<double> _rbps0 = new GivenParameter<double>(50.0);

        [ParameterName("rbpsl"), ParameterInfo("Body resistance RBPS L scaling")]
        [Finite]
        private GivenParameter<double> _rbpsl = new GivenParameter<double>(0.0);

        [ParameterName("rbpsw"), ParameterInfo("Body resistance RBPS W scaling")]
        [Finite]
        private GivenParameter<double> _rbpsw = new GivenParameter<double>(0.0);

        [ParameterName("rbpsnf"), ParameterInfo("Body resistance RBPS NF scaling")]
        [Finite]
        private GivenParameter<double> _rbpsnf = new GivenParameter<double>(0.0);

        [ParameterName("rbpd0"), ParameterInfo("Body resistance RBPD scaling")]
        [Finite]
        private GivenParameter<double> _rbpd0 = new GivenParameter<double>(50.0);

        [ParameterName("rbpdl"), ParameterInfo("Body resistance RBPD L scaling")]
        [Finite]
        private GivenParameter<double> _rbpdl = new GivenParameter<double>(0.0);

        [ParameterName("rbpdw"), ParameterInfo("Body resistance RBPD W scaling")]
        [Finite]
        private GivenParameter<double> _rbpdw = new GivenParameter<double>(0.0);

        [ParameterName("rbpdnf"), ParameterInfo("Body resistance RBPD NF scaling")]
        [Finite]
        private GivenParameter<double> _rbpdnf = new GivenParameter<double>(0.0);

        [ParameterName("rbpbx0"), ParameterInfo("Body resistance RBPBX  scaling")]
        [Finite]
        private GivenParameter<double> _rbpbx0 = new GivenParameter<double>(100.0);

        [ParameterName("rbpbxl"), ParameterInfo("Body resistance RBPBX L scaling")]
        [Finite]
        private GivenParameter<double> _rbpbxl = new GivenParameter<double>(0.0);

        [ParameterName("rbpbxw"), ParameterInfo("Body resistance RBPBX W scaling")]
        [Finite]
        private GivenParameter<double> _rbpbxw = new GivenParameter<double>(0.0);

        [ParameterName("rbpbxnf"), ParameterInfo("Body resistance RBPBX NF scaling")]
        [Finite]
        private GivenParameter<double> _rbpbxnf = new GivenParameter<double>(0.0);

        [ParameterName("rbpby0"), ParameterInfo("Body resistance RBPBY  scaling")]
        [Finite]
        private GivenParameter<double> _rbpby0 = new GivenParameter<double>(100.0);

        [ParameterName("rbpbyl"), ParameterInfo("Body resistance RBPBY L scaling")]
        [Finite]
        private GivenParameter<double> _rbpbyl = new GivenParameter<double>(0.0);

        [ParameterName("rbpbyw"), ParameterInfo("Body resistance RBPBY W scaling")]
        [Finite]
        private GivenParameter<double> _rbpbyw = new GivenParameter<double>(0.0);

        [ParameterName("rbpbynf"), ParameterInfo("Body resistance RBPBY NF scaling")]
        [Finite]
        private GivenParameter<double> _rbpbynf = new GivenParameter<double>(0.0);

        [ParameterName("rbsbx0"), ParameterInfo("Body resistance RBSBX  scaling")]
        [Finite]
        private GivenParameter<double> _rbsbx0 = new GivenParameter<double>(100.0);

        [ParameterName("rbsby0"), ParameterInfo("Body resistance RBSBY  scaling")]
        [Finite]
        private GivenParameter<double> _rbsby0 = new GivenParameter<double>(100.0);

        [ParameterName("rbdbx0"), ParameterInfo("Body resistance RBDBX  scaling")]
        [Finite]
        private GivenParameter<double> _rbdbx0 = new GivenParameter<double>(100.0);

        [ParameterName("rbdby0"), ParameterInfo("Body resistance RBDBY  scaling")]
        [Finite]
        private GivenParameter<double> _rbdby0 = new GivenParameter<double>(100.0);

        [ParameterName("rbsdbxl"), ParameterInfo("Body resistance RBSDBX L scaling")]
        [Finite]
        private GivenParameter<double> _rbsdbxl = new GivenParameter<double>(0.0);

        [ParameterName("rbsdbxw"), ParameterInfo("Body resistance RBSDBX W scaling")]
        [Finite]
        private GivenParameter<double> _rbsdbxw = new GivenParameter<double>(0.0);

        [ParameterName("rbsdbxnf"), ParameterInfo("Body resistance RBSDBX NF scaling")]
        [Finite]
        private GivenParameter<double> _rbsdbxnf = new GivenParameter<double>(0.0);

        [ParameterName("rbsdbyl"), ParameterInfo("Body resistance RBSDBY L scaling")]
        [Finite]
        private GivenParameter<double> _rbsdbyl = new GivenParameter<double>(0.0);

        [ParameterName("rbsdbyw"), ParameterInfo("Body resistance RBSDBY W scaling")]
        [Finite]
        private GivenParameter<double> _rbsdbyw = new GivenParameter<double>(0.0);

        [ParameterName("rbsdbynf"), ParameterInfo("Body resistance RBSDBY NF scaling")]
        [Finite]
        private GivenParameter<double> _rbsdbynf = new GivenParameter<double>(0.0);

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

        [ParameterName("lndep"), ParameterInfo("Length dependence of ndep")]
        [Finite]
        public GivenParameter<double> Lndep
        {
            get => _lndep;
            set
            {
                Utility.Finite(value, nameof(Lndep));
                if (value > 1.0e20)
                    _lndep = value * 1e-6;
                else
                    _lndep = value;
            }
        }
        private GivenParameter<double> _lndep = new GivenParameter<double>(0.0);

        [ParameterName("lnsd"), ParameterInfo("Length dependence of nsd")]
        [Finite]
        public GivenParameter<double> Lnsd
        {
            get => _lnsd;
            set
            {
                Utility.Finite(value, nameof(Lnsd));
                if (value > 1.0e23)
                    _lnsd = value * 1e-6;
                else
                    _lnsd = value;
            }
        }
        private GivenParameter<double> _lnsd = new GivenParameter<double>(0.0);

        [ParameterName("lphin"), ParameterInfo("Length dependence of phin")]
        [Finite]
        private GivenParameter<double> _lphin = new GivenParameter<double>(0.0);

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
        private GivenParameter<double> _lngate = new GivenParameter<double>(0.0);

        [ParameterName("lgamma1"), ParameterInfo("Length dependence of gamma1")]
        [Finite]
        private GivenParameter<double> _lgamma1 = new GivenParameter<double>(0.0);

        [ParameterName("lgamma2"), ParameterInfo("Length dependence of gamma2")]
        [Finite]
        private GivenParameter<double> _lgamma2 = new GivenParameter<double>(0.0);

        [ParameterName("lvbx"), ParameterInfo("Length dependence of vbx")]
        [Finite]
        private GivenParameter<double> _lvbx = new GivenParameter<double>(0.0);

        [ParameterName("lvbm"), ParameterInfo("Length dependence of vbm")]
        [Finite]
        private GivenParameter<double> _lvbm = new GivenParameter<double>(0.0);

        [ParameterName("lxt"), ParameterInfo("Length dependence of xt")]
        [Finite]
        private GivenParameter<double> _lxt = new GivenParameter<double>(0.0);

        [ParameterName("lk1"), ParameterInfo("Length dependence of k1")]
        [Finite]
        private GivenParameter<double> _lk1 = new GivenParameter<double>(0.0);

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
        private GivenParameter<double> _lk2 = new GivenParameter<double>(0.0);

        [ParameterName("lk3"), ParameterInfo("Length dependence of k3")]
        [Finite]
        private GivenParameter<double> _lk3 = new GivenParameter<double>(0.0);

        [ParameterName("lk3b"), ParameterInfo("Length dependence of k3b")]
        [Finite]
        private GivenParameter<double> _lk3b = new GivenParameter<double>(0.0);

        [ParameterName("lw0"), ParameterInfo("Length dependence of w0")]
        [Finite]
        private GivenParameter<double> _lw0 = new GivenParameter<double>(0.0);

        [ParameterName("ldvtp0"), ParameterInfo("Length dependence of dvtp0")]
        [Finite]
        private GivenParameter<double> _ldvtp0 = new GivenParameter<double>(0.0);

        [ParameterName("ldvtp1"), ParameterInfo("Length dependence of dvtp1")]
        [Finite]
        private GivenParameter<double> _ldvtp1 = new GivenParameter<double>(0.0);

        [ParameterName("ldvtp2"), ParameterInfo("Length dependence of dvtp2")]
        [Finite]
        private GivenParameter<double> _ldvtp2 = new GivenParameter<double>(0.0);

        [ParameterName("ldvtp3"), ParameterInfo("Length dependence of dvtp3")]
        [Finite]
        private GivenParameter<double> _ldvtp3 = new GivenParameter<double>(0.0);

        [ParameterName("ldvtp4"), ParameterInfo("Length dependence of dvtp4")]
        [Finite]
        private GivenParameter<double> _ldvtp4 = new GivenParameter<double>(0.0);

        [ParameterName("ldvtp5"), ParameterInfo("Length dependence of dvtp5")]
        [Finite]
        private GivenParameter<double> _ldvtp5 = new GivenParameter<double>(0.0);

        [ParameterName("llpe0"), ParameterInfo("Length dependence of lpe0")]
        [Finite]
        private GivenParameter<double> _llpe0 = new GivenParameter<double>(0.0);

        [ParameterName("llpeb"), ParameterInfo("Length dependence of lpeb")]
        [Finite]
        private GivenParameter<double> _llpeb = new GivenParameter<double>(0.0);

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

        [ParameterName("lvth0"), ParameterName("lvtho"), ParameterInfo("Length dependence of vth0")]
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

        [ParameterName("lud"), ParameterInfo("Length dependence of ud")]
        [Finite]
        private GivenParameter<double> _lud = new GivenParameter<double>(0.0);

        [ParameterName("lud1"), ParameterInfo("Length dependence of ud1")]
        [Finite]
        private GivenParameter<double> _lud1 = new GivenParameter<double>(0.0);

        [ParameterName("lup"), ParameterInfo("Length dependence of up")]
        [Finite]
        private GivenParameter<double> _lup = new GivenParameter<double>(0.0);

        [ParameterName("llp"), ParameterInfo("Length dependence of lp")]
        [Finite]
        private GivenParameter<double> _llp = new GivenParameter<double>(0.0);

        [ParameterName("lu0"), ParameterInfo("Length dependence of u0")]
        [Finite]
        private GivenParameter<double> _lu0 = new GivenParameter<double>(0.0);

        [ParameterName("lute"), ParameterInfo("Length dependence of ute")]
        [Finite]
        private GivenParameter<double> _lute = new GivenParameter<double>(0.0);

        [ParameterName("lucste"), ParameterInfo("Length dependence of ucste")]
        [Finite]
        private GivenParameter<double> _lucste = new GivenParameter<double>(0.0);

        [ParameterName("lvoff"), ParameterInfo("Length dependence of voff")]
        [Finite]
        private GivenParameter<double> _lvoff = new GivenParameter<double>(0.0);

        [ParameterName("lminv"), ParameterInfo("Length dependence of minv")]
        [Finite]
        private GivenParameter<double> _lminv = new GivenParameter<double>(0.0);

        [ParameterName("lminvcv"), ParameterInfo("Length dependence of minvcv")]
        [Finite]
        private GivenParameter<double> _lminvcv = new GivenParameter<double>(0.0);

        [ParameterName("ldelta"), ParameterInfo("Length dependence of delta")]
        [Finite]
        private GivenParameter<double> _ldelta = new GivenParameter<double>(0.0);

        [ParameterName("lrdsw"), ParameterInfo("Length dependence of rdsw ")]
        [Finite]
        private GivenParameter<double> _lrdsw = new GivenParameter<double>(0.0);

        [ParameterName("lrsw"), ParameterInfo("Length dependence of rsw")]
        [Finite]
        private GivenParameter<double> _lrsw = new GivenParameter<double>(0.0);

        [ParameterName("lrdw"), ParameterInfo("Length dependence of rdw")]
        [Finite]
        private GivenParameter<double> _lrdw = new GivenParameter<double>(0.0);

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

        [ParameterName("lfprout"), ParameterInfo("Length dependence of pdiblcb")]
        [Finite]
        private GivenParameter<double> _lfprout = new GivenParameter<double>(0.0);

        [ParameterName("lpdits"), ParameterInfo("Length dependence of pdits")]
        [Finite]
        private GivenParameter<double> _lpdits = new GivenParameter<double>(0.0);

        [ParameterName("lpditsd"), ParameterInfo("Length dependence of pditsd")]
        [Finite]
        private GivenParameter<double> _lpditsd = new GivenParameter<double>(0.0);

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

        [ParameterName("lckappas"), ParameterInfo("Length dependence of ckappas")]
        [Finite]
        private GivenParameter<double> _lckappas = new GivenParameter<double>(0.0);

        [ParameterName("lckappad"), ParameterInfo("Length dependence of ckappad")]
        [Finite]
        private GivenParameter<double> _lckappad = new GivenParameter<double>(0.0);

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

        [ParameterName("lagidl"), ParameterInfo("Length dependence of agidl")]
        [Finite]
        private GivenParameter<double> _lagidl = new GivenParameter<double>(0.0);

        [ParameterName("lbgidl"), ParameterInfo("Length dependence of bgidl")]
        [Finite]
        private GivenParameter<double> _lbgidl = new GivenParameter<double>(0.0);

        [ParameterName("lcgidl"), ParameterInfo("Length dependence of cgidl")]
        [Finite]
        private GivenParameter<double> _lcgidl = new GivenParameter<double>(0.0);

        [ParameterName("lrgidl"), ParameterInfo("Length dependence of rgidl")]
        [Finite]
        private GivenParameter<double> _lrgidl = new GivenParameter<double>(0.0);

        [ParameterName("lkgidl"), ParameterInfo("Length dependence of kgidl")]
        [Finite]
        private GivenParameter<double> _lkgidl = new GivenParameter<double>(0.0);

        [ParameterName("lfgidl"), ParameterInfo("Length dependence of fgidl")]
        [Finite]
        private GivenParameter<double> _lfgidl = new GivenParameter<double>(0.0);

        [ParameterName("legidl"), ParameterInfo("Length dependence of egidl")]
        [Finite]
        private GivenParameter<double> _legidl = new GivenParameter<double>(0.0);

        [ParameterName("lagisl"), ParameterInfo("Length dependence of agisl")]
        [Finite]
        private GivenParameter<double> _lagisl = new GivenParameter<double>();

        [ParameterName("lbgisl"), ParameterInfo("Length dependence of bgisl")]
        [Finite]
        private GivenParameter<double> _lbgisl = new GivenParameter<double>();

        [ParameterName("lcgisl"), ParameterInfo("Length dependence of cgisl")]
        [Finite]
        private GivenParameter<double> _lcgisl = new GivenParameter<double>();

        [ParameterName("lrgisl"), ParameterInfo("Length dependence of rgisl")]
        [Finite]
        private GivenParameter<double> _lrgisl = new GivenParameter<double>();

        [ParameterName("lkgisl"), ParameterInfo("Length dependence of kgisl")]
        [Finite]
        private GivenParameter<double> _lkgisl = new GivenParameter<double>();

        [ParameterName("lfgisl"), ParameterInfo("Length dependence of fgisl")]
        [Finite]
        private GivenParameter<double> _lfgisl = new GivenParameter<double>();

        [ParameterName("legisl"), ParameterInfo("Length dependence of egisl")]
        [Finite]
        private GivenParameter<double> _legisl = new GivenParameter<double>();

        [ParameterName("laigc"), ParameterInfo("Length dependence of aigc")]
        [Finite]
        private GivenParameter<double> _laigc = new GivenParameter<double>(0.0);

        [ParameterName("lbigc"), ParameterInfo("Length dependence of bigc")]
        [Finite]
        private GivenParameter<double> _lbigc = new GivenParameter<double>(0.0);

        [ParameterName("lcigc"), ParameterInfo("Length dependence of cigc")]
        [Finite]
        private GivenParameter<double> _lcigc = new GivenParameter<double>(0.0);

        [ParameterName("laigsd"), ParameterInfo("Length dependence of aigsd")]
        [Finite]
        private GivenParameter<double> _laigsd = new GivenParameter<double>();

        [ParameterName("lbigsd"), ParameterInfo("Length dependence of bigsd")]
        [Finite]
        private GivenParameter<double> _lbigsd = new GivenParameter<double>();

        [ParameterName("lcigsd"), ParameterInfo("Length dependence of cigsd")]
        [Finite]
        private GivenParameter<double> _lcigsd = new GivenParameter<double>();

        [ParameterName("laigs"), ParameterInfo("Length dependence of aigs")]
        [Finite]
        private GivenParameter<double> _laigs = new GivenParameter<double>();

        [ParameterName("lbigs"), ParameterInfo("Length dependence of bigs")]
        [Finite]
        private GivenParameter<double> _lbigs = new GivenParameter<double>();

        [ParameterName("lcigs"), ParameterInfo("Length dependence of cigs")]
        [Finite]
        private GivenParameter<double> _lcigs = new GivenParameter<double>();

        [ParameterName("laigd"), ParameterInfo("Length dependence of aigd")]
        [Finite]
        private GivenParameter<double> _laigd = new GivenParameter<double>(0.0);

        [ParameterName("lbigd"), ParameterInfo("Length dependence of bigd")]
        [Finite]
        private GivenParameter<double> _lbigd = new GivenParameter<double>(0.0);

        [ParameterName("lcigd"), ParameterInfo("Length dependence of cigd")]
        [Finite]
        private GivenParameter<double> _lcigd = new GivenParameter<double>(0.0);

        [ParameterName("laigbacc"), ParameterInfo("Length dependence of aigbacc")]
        [Finite]
        private GivenParameter<double> _laigbacc = new GivenParameter<double>(0.0);

        [ParameterName("lbigbacc"), ParameterInfo("Length dependence of bigbacc")]
        [Finite]
        private GivenParameter<double> _lbigbacc = new GivenParameter<double>(0.0);

        [ParameterName("lcigbacc"), ParameterInfo("Length dependence of cigbacc")]
        [Finite]
        private GivenParameter<double> _lcigbacc = new GivenParameter<double>(0.0);

        [ParameterName("laigbinv"), ParameterInfo("Length dependence of aigbinv")]
        [Finite]
        private GivenParameter<double> _laigbinv = new GivenParameter<double>(0.0);

        [ParameterName("lbigbinv"), ParameterInfo("Length dependence of bigbinv")]
        [Finite]
        private GivenParameter<double> _lbigbinv = new GivenParameter<double>(0.0);

        [ParameterName("lcigbinv"), ParameterInfo("Length dependence of cigbinv")]
        [Finite]
        private GivenParameter<double> _lcigbinv = new GivenParameter<double>(0.0);

        [ParameterName("lnigc"), ParameterInfo("Length dependence of nigc")]
        [Finite]
        private GivenParameter<double> _lnigc = new GivenParameter<double>(0.0);

        [ParameterName("lnigbinv"), ParameterInfo("Length dependence of nigbinv")]
        [Finite]
        private GivenParameter<double> _lnigbinv = new GivenParameter<double>(0.0);

        [ParameterName("lnigbacc"), ParameterInfo("Length dependence of nigbacc")]
        [Finite]
        private GivenParameter<double> _lnigbacc = new GivenParameter<double>(0.0);

        [ParameterName("lntox"), ParameterInfo("Length dependence of ntox")]
        [Finite]
        private GivenParameter<double> _lntox = new GivenParameter<double>(0.0);

        [ParameterName("leigbinv"), ParameterInfo("Length dependence for eigbinv")]
        [Finite]
        private GivenParameter<double> _leigbinv = new GivenParameter<double>(0.0);

        [ParameterName("lpigcd"), ParameterInfo("Length dependence for pigcd")]
        [Finite]
        private GivenParameter<double> _lpigcd = new GivenParameter<double>(0.0);

        [ParameterName("lpoxedge"), ParameterInfo("Length dependence for poxedge")]
        [Finite]
        private GivenParameter<double> _lpoxedge = new GivenParameter<double>(0.0);

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

        [ParameterName("lxrcrg1"), ParameterInfo("Length dependence of xrcrg1")]
        [Finite]
        private GivenParameter<double> _lxrcrg1 = new GivenParameter<double>(0.0);

        [ParameterName("lxrcrg2"), ParameterInfo("Length dependence of xrcrg2")]
        [Finite]
        private GivenParameter<double> _lxrcrg2 = new GivenParameter<double>(0.0);

        [ParameterName("llambda"), ParameterInfo("Length dependence of lambda")]
        [Finite]
        private GivenParameter<double> _llambda = new GivenParameter<double>(0.0);

        [ParameterName("lvtl"), ParameterInfo(" Length dependence of vtl")]
        [Finite]
        private GivenParameter<double> _lvtl = new GivenParameter<double>(0.0);

        [ParameterName("lxn"), ParameterInfo(" Length dependence of xn")]
        [Finite]
        private GivenParameter<double> _lxn = new GivenParameter<double>(0.0);

        [ParameterName("leu"), ParameterInfo(" Length dependence of eu")]
        [Finite]
        private GivenParameter<double> _leu = new GivenParameter<double>(0.0);

        [ParameterName("lucs"), ParameterInfo("Length dependence of lucs")]
        [Finite]
        private GivenParameter<double> _lucs = new GivenParameter<double>(0.0);

        [ParameterName("lvfbsdoff"), ParameterInfo("Length dependence of vfbsdoff")]
        [Finite]
        private GivenParameter<double> _lvfbsdoff = new GivenParameter<double>(0.0);

        [ParameterName("ltvfbsdoff"), ParameterInfo("Length dependence of tvfbsdoff")]
        [Finite]
        private GivenParameter<double> _ltvfbsdoff = new GivenParameter<double>(0.0);

        [ParameterName("ltvoff"), ParameterInfo("Length dependence of tvoff")]
        [Finite]
        private GivenParameter<double> _ltvoff = new GivenParameter<double>(0.0);

        [ParameterName("ltnfactor"), ParameterInfo("Length dependence of tnfactor")]
        [Finite]
        private GivenParameter<double> _ltnfactor = new GivenParameter<double>(0.0);

        [ParameterName("lteta0"), ParameterInfo("Length dependence of teta0")]
        [Finite]
        private GivenParameter<double> _lteta0 = new GivenParameter<double>(0.0);

        [ParameterName("ltvoffcv"), ParameterInfo("Length dependence of tvoffcv")]
        [Finite]
        private GivenParameter<double> _ltvoffcv = new GivenParameter<double>(0.0);

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

        [ParameterName("wndep"), ParameterInfo("Width dependence of ndep")]
        [Finite]
        public GivenParameter<double> Wndep
        {
            get => _wndep;
            set
            {
                Utility.Finite(value, nameof(Wndep));
                if (value > 1.0e20)
                    _wndep = value * 1e-6;
                else
                    _wndep = value;
            }
        }
        private GivenParameter<double> _wndep = new GivenParameter<double>(0.0);

        [ParameterName("wnsd"), ParameterInfo("Width dependence of nsd")]
        [Finite]
        public GivenParameter<double> Wnsd
        {
            get => _wnsd ;
            set
            {
                Utility.Finite(value, nameof(Wnsd));
                if (value > 1.0e23)
                    _wnsd = value * 1e-6;
                else
                    _wnsd = value;
            }
        }
        private GivenParameter<double> _wnsd = new GivenParameter<double>(0.0);

        [ParameterName("wphin"), ParameterInfo("Width dependence of phin")]
        [Finite]
        private GivenParameter<double> _wphin = new GivenParameter<double>(0.0);

        [ParameterName("wngate"), ParameterInfo("Width dependence of ngate")]
        [Finite]
        public GivenParameter<double> Wngate
        {
            get => _wngate ;
            set
            {
                Utility.Finite(value, nameof(Wngate));
                if (value > 1.0e23)
                    _wngate = value * 1e-6;
                else
                    _wngate = value;
            }
        }
        private GivenParameter<double> _wngate = new GivenParameter<double>(0.0);

        [ParameterName("wgamma1"), ParameterInfo("Width dependence of gamma1")]
        [Finite]
        private GivenParameter<double> _wgamma1 = new GivenParameter<double>(0.0);

        [ParameterName("wgamma2"), ParameterInfo("Width dependence of gamma2")]
        [Finite]
        private GivenParameter<double> _wgamma2 = new GivenParameter<double>(0.0);

        [ParameterName("wvbx"), ParameterInfo("Width dependence of vbx")]
        [Finite]
        private GivenParameter<double> _wvbx = new GivenParameter<double>(0.0);

        [ParameterName("wvbm"), ParameterInfo("Width dependence of vbm")]
        [Finite]
        private GivenParameter<double> _wvbm = new GivenParameter<double>(0.0);

        [ParameterName("wxt"), ParameterInfo("Width dependence of xt")]
        [Finite]
        private GivenParameter<double> _wxt = new GivenParameter<double>(0.0);

        [ParameterName("wk1"), ParameterInfo("Width dependence of k1")]
        [Finite]
        private GivenParameter<double> _wk1 = new GivenParameter<double>(0.0);

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
        private GivenParameter<double> _wk2 = new GivenParameter<double>(0.0);

        [ParameterName("wk3"), ParameterInfo("Width dependence of k3")]
        [Finite]
        private GivenParameter<double> _wk3 = new GivenParameter<double>(0.0);

        [ParameterName("wk3b"), ParameterInfo("Width dependence of k3b")]
        [Finite]
        private GivenParameter<double> _wk3b = new GivenParameter<double>(0.0);

        [ParameterName("ww0"), ParameterInfo("Width dependence of w0")]
        [Finite]
        private GivenParameter<double> _ww0 = new GivenParameter<double>(0.0);

        [ParameterName("wdvtp0"), ParameterInfo("Width dependence of dvtp0")]
        [Finite]
        private GivenParameter<double> _wdvtp0 = new GivenParameter<double>(0.0);

        [ParameterName("wdvtp1"), ParameterInfo("Width dependence of dvtp1")]
        [Finite]
        private GivenParameter<double> _wdvtp1 = new GivenParameter<double>(0.0);

        [ParameterName("wdvtp2"), ParameterInfo("Width dependence of dvtp2")]
        [Finite]
        private GivenParameter<double> _wdvtp2 = new GivenParameter<double>(0.0);

        [ParameterName("wdvtp3"), ParameterInfo("Width dependence of dvtp3")]
        [Finite]
        private GivenParameter<double> _wdvtp3 = new GivenParameter<double>(0.0);

        [ParameterName("wdvtp4"), ParameterInfo("Width dependence of dvtp4")]
        [Finite]
        private GivenParameter<double> _wdvtp4 = new GivenParameter<double>(0.0);

        [ParameterName("wdvtp5"), ParameterInfo("Width dependence of dvtp5")]
        [Finite]
        private GivenParameter<double> _wdvtp5 = new GivenParameter<double>(0.0);

        [ParameterName("wlpe0"), ParameterInfo("Width dependence of lpe0")]
        [Finite]
        private GivenParameter<double> _wlpe0 = new GivenParameter<double>(0.0);

        [ParameterName("wlpeb"), ParameterInfo("Width dependence of lpeb")]
        [Finite]
        private GivenParameter<double> _wlpeb = new GivenParameter<double>(0.0);

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

        [ParameterName("wvth0"), ParameterName("wvtho"), ParameterInfo("Width dependence of vth0")]
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

        [ParameterName("wud"), ParameterInfo("Width dependence of ud")]
        [Finite]
        private GivenParameter<double> _wud = new GivenParameter<double>(0.0);

        [ParameterName("wud1"), ParameterInfo("Width dependence of ud1")]
        [Finite]
        private GivenParameter<double> _wud1 = new GivenParameter<double>(0.0);

        [ParameterName("wup"), ParameterInfo("Width dependence of up")]
        [Finite]
        private GivenParameter<double> _wup = new GivenParameter<double>(0.0);

        [ParameterName("wlp"), ParameterInfo("Width dependence of lp")]
        [Finite]
        private GivenParameter<double> _wlp = new GivenParameter<double>(0.0);

        [ParameterName("wu0"), ParameterInfo("Width dependence of u0")]
        [Finite]
        private GivenParameter<double> _wu0 = new GivenParameter<double>(0.0);

        [ParameterName("wute"), ParameterInfo("Width dependence of ute")]
        [Finite]
        private GivenParameter<double> _wute = new GivenParameter<double>(0.0);

        [ParameterName("wucste"), ParameterInfo("Width dependence of ucste")]
        [Finite]
        private GivenParameter<double> _wucste = new GivenParameter<double>(0.0);

        [ParameterName("wvoff"), ParameterInfo("Width dependence of voff")]
        [Finite]
        private GivenParameter<double> _wvoff = new GivenParameter<double>(0.0);

        [ParameterName("wminv"), ParameterInfo("Width dependence of minv")]
        [Finite]
        private GivenParameter<double> _wminv = new GivenParameter<double>(0.0);

        [ParameterName("wminvcv"), ParameterInfo("Width dependence of minvcv")]
        [Finite]
        private GivenParameter<double> _wminvcv = new GivenParameter<double>(0.0);

        [ParameterName("wdelta"), ParameterInfo("Width dependence of delta")]
        [Finite]
        private GivenParameter<double> _wdelta = new GivenParameter<double>(0.0);

        [ParameterName("wrdsw"), ParameterInfo("Width dependence of rdsw ")]
        [Finite]
        private GivenParameter<double> _wrdsw = new GivenParameter<double>(0.0);

        [ParameterName("wrsw"), ParameterInfo("Width dependence of rsw")]
        [Finite]
        private GivenParameter<double> _wrsw = new GivenParameter<double>(0.0);

        [ParameterName("wrdw"), ParameterInfo("Width dependence of rdw")]
        [Finite]
        private GivenParameter<double> _wrdw = new GivenParameter<double>(0.0);

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

        [ParameterName("wfprout"), ParameterInfo("Width dependence of pdiblcb")]
        [Finite]
        private GivenParameter<double> _wfprout = new GivenParameter<double>(0.0);

        [ParameterName("wpdits"), ParameterInfo("Width dependence of pdits")]
        [Finite]
        private GivenParameter<double> _wpdits = new GivenParameter<double>(0.0);

        [ParameterName("wpditsd"), ParameterInfo("Width dependence of pditsd")]
        [Finite]
        private GivenParameter<double> _wpditsd = new GivenParameter<double>(0.0);

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

        [ParameterName("wckappas"), ParameterInfo("Width dependence of ckappas")]
        [Finite]
        private GivenParameter<double> _wckappas = new GivenParameter<double>(0.0);

        [ParameterName("wckappad"), ParameterInfo("Width dependence of ckappad")]
        [Finite]
        private GivenParameter<double> _wckappad = new GivenParameter<double>(0.0);

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

        [ParameterName("wagidl"), ParameterInfo("Width dependence of agidl")]
        [Finite]
        private GivenParameter<double> _wagidl = new GivenParameter<double>(0.0);

        [ParameterName("wbgidl"), ParameterInfo("Width dependence of bgidl")]
        [Finite]
        private GivenParameter<double> _wbgidl = new GivenParameter<double>(0.0);

        [ParameterName("wcgidl"), ParameterInfo("Width dependence of cgidl")]
        [Finite]
        private GivenParameter<double> _wcgidl = new GivenParameter<double>(0.0);

        [ParameterName("wrgidl"), ParameterInfo("Width dependence of rgidl")]
        [Finite]
        private GivenParameter<double> _wrgidl = new GivenParameter<double>(0.0);

        [ParameterName("wkgidl"), ParameterInfo("Width dependence of kgidl")]
        [Finite]
        private GivenParameter<double> _wkgidl = new GivenParameter<double>(0.0);

        [ParameterName("wfgidl"), ParameterInfo("Width dependence of fgidl")]
        [Finite]
        private GivenParameter<double> _wfgidl = new GivenParameter<double>(0.0);

        [ParameterName("wegidl"), ParameterInfo("Width dependence of egidl")]
        [Finite]
        private GivenParameter<double> _wegidl = new GivenParameter<double>(0.0);

        [ParameterName("wagisl"), ParameterInfo("Width dependence of agisl")]
        [Finite]
        private GivenParameter<double> _wagisl = new GivenParameter<double>();

        [ParameterName("wbgisl"), ParameterInfo("Width dependence of bgisl")]
        [Finite]
        private GivenParameter<double> _wbgisl = new GivenParameter<double>();

        [ParameterName("wcgisl"), ParameterInfo("Width dependence of cgisl")]
        [Finite]
        private GivenParameter<double> _wcgisl = new GivenParameter<double>();

        [ParameterName("wrgisl"), ParameterInfo("Width dependence of rgisl")]
        [Finite]
        private GivenParameter<double> _wrgisl = new GivenParameter<double>();

        [ParameterName("wkgisl"), ParameterInfo("Width dependence of kgisl")]
        [Finite]
        private GivenParameter<double> _wkgisl = new GivenParameter<double>();

        [ParameterName("wfgisl"), ParameterInfo("Width dependence of fgisl")]
        [Finite]
        private GivenParameter<double> _wfgisl = new GivenParameter<double>();

        [ParameterName("wegisl"), ParameterInfo("Width dependence of egisl")]
        [Finite]
        private GivenParameter<double> _wegisl = new GivenParameter<double>();

        [ParameterName("waigc"), ParameterInfo("Width dependence of aigc")]
        [Finite]
        private GivenParameter<double> _waigc = new GivenParameter<double>(0.0);

        [ParameterName("wbigc"), ParameterInfo("Width dependence of bigc")]
        [Finite]
        private GivenParameter<double> _wbigc = new GivenParameter<double>(0.0);

        [ParameterName("wcigc"), ParameterInfo("Width dependence of cigc")]
        [Finite]
        private GivenParameter<double> _wcigc = new GivenParameter<double>(0.0);

        [ParameterName("waigsd"), ParameterInfo("Width dependence of aigsd")]
        [Finite]
        private GivenParameter<double> _waigsd = new GivenParameter<double>();

        [ParameterName("wbigsd"), ParameterInfo("Width dependence of bigsd")]
        [Finite]
        private GivenParameter<double> _wbigsd = new GivenParameter<double>();

        [ParameterName("wcigsd"), ParameterInfo("Width dependence of cigsd")]
        [Finite]
        private GivenParameter<double> _wcigsd = new GivenParameter<double>();

        [ParameterName("waigs"), ParameterInfo("Width dependence of aigs")]
        [Finite]
        private GivenParameter<double> _waigs = new GivenParameter<double>();

        [ParameterName("wbigs"), ParameterInfo("Width dependence of bigs")]
        [Finite]
        private GivenParameter<double> _wbigs = new GivenParameter<double>();

        [ParameterName("wcigs"), ParameterInfo("Width dependence of cigs")]
        [Finite]
        private GivenParameter<double> _wcigs = new GivenParameter<double>();

        [ParameterName("waigd"), ParameterInfo("Width dependence of aigd")]
        [Finite]
        private GivenParameter<double> _waigd = new GivenParameter<double>(0.0);

        [ParameterName("wbigd"), ParameterInfo("Width dependence of bigd")]
        [Finite]
        private GivenParameter<double> _wbigd = new GivenParameter<double>(0.0);

        [ParameterName("wcigd"), ParameterInfo("Width dependence of cigd")]
        [Finite]
        private GivenParameter<double> _wcigd = new GivenParameter<double>(0.0);

        [ParameterName("waigbacc"), ParameterInfo("Width dependence of aigbacc")]
        [Finite]
        private GivenParameter<double> _waigbacc = new GivenParameter<double>(0.0);

        [ParameterName("wbigbacc"), ParameterInfo("Width dependence of bigbacc")]
        [Finite]
        private GivenParameter<double> _wbigbacc = new GivenParameter<double>(0.0);

        [ParameterName("wcigbacc"), ParameterInfo("Width dependence of cigbacc")]
        [Finite]
        private GivenParameter<double> _wcigbacc = new GivenParameter<double>(0.0);

        [ParameterName("waigbinv"), ParameterInfo("Width dependence of aigbinv")]
        [Finite]
        private GivenParameter<double> _waigbinv = new GivenParameter<double>(0.0);

        [ParameterName("wbigbinv"), ParameterInfo("Width dependence of bigbinv")]
        [Finite]
        private GivenParameter<double> _wbigbinv = new GivenParameter<double>(0.0);

        [ParameterName("wcigbinv"), ParameterInfo("Width dependence of cigbinv")]
        [Finite]
        private GivenParameter<double> _wcigbinv = new GivenParameter<double>(0.0);

        [ParameterName("wnigc"), ParameterInfo("Width dependence of nigc")]
        [Finite]
        private GivenParameter<double> _wnigc = new GivenParameter<double>(0.0);

        [ParameterName("wnigbinv"), ParameterInfo("Width dependence of nigbinv")]
        [Finite]
        private GivenParameter<double> _wnigbinv = new GivenParameter<double>(0.0);

        [ParameterName("wnigbacc"), ParameterInfo("Width dependence of nigbacc")]
        [Finite]
        private GivenParameter<double> _wnigbacc = new GivenParameter<double>(0.0);

        [ParameterName("wntox"), ParameterInfo("Width dependence of ntox")]
        [Finite]
        private GivenParameter<double> _wntox = new GivenParameter<double>(0.0);

        [ParameterName("weigbinv"), ParameterInfo("Width dependence for eigbinv")]
        [Finite]
        private GivenParameter<double> _weigbinv = new GivenParameter<double>(0.0);

        [ParameterName("wpigcd"), ParameterInfo("Width dependence for pigcd")]
        [Finite]
        private GivenParameter<double> _wpigcd = new GivenParameter<double>(0.0);

        [ParameterName("wpoxedge"), ParameterInfo("Width dependence for poxedge")]
        [Finite]
        private GivenParameter<double> _wpoxedge = new GivenParameter<double>(0.0);

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

        [ParameterName("wxrcrg1"), ParameterInfo("Width dependence of xrcrg1")]
        [Finite]
        private GivenParameter<double> _wxrcrg1 = new GivenParameter<double>(0.0);

        [ParameterName("wxrcrg2"), ParameterInfo("Width dependence of xrcrg2")]
        [Finite]
        private GivenParameter<double> _wxrcrg2 = new GivenParameter<double>(0.0);

        [ParameterName("wlambda"), ParameterInfo("Width dependence of lambda")]
        [Finite]
        private GivenParameter<double> _wlambda = new GivenParameter<double>(0.0);

        [ParameterName("wvtl"), ParameterInfo("Width dependence of vtl")]
        [Finite]
        private GivenParameter<double> _wvtl = new GivenParameter<double>(0.0);

        [ParameterName("wxn"), ParameterInfo("Width dependence of xn")]
        [Finite]
        private GivenParameter<double> _wxn = new GivenParameter<double>(0.0);

        [ParameterName("weu"), ParameterInfo("Width dependence of eu")]
        [Finite]
        private GivenParameter<double> _weu = new GivenParameter<double>(0.0);

        [ParameterName("wucs"), ParameterInfo("Width dependence of ucs")]
        [Finite]
        private GivenParameter<double> _wucs = new GivenParameter<double>(0.0);

        [ParameterName("wvfbsdoff"), ParameterInfo("Width dependence of vfbsdoff")]
        [Finite]
        private GivenParameter<double> _wvfbsdoff = new GivenParameter<double>(0.0);

        [ParameterName("wtvfbsdoff"), ParameterInfo("Width dependence of tvfbsdoff")]
        [Finite]
        private GivenParameter<double> _wtvfbsdoff = new GivenParameter<double>(0.0);

        [ParameterName("wtvoff"), ParameterInfo("Width dependence of tvoff")]
        [Finite]
        private GivenParameter<double> _wtvoff = new GivenParameter<double>(0.0);

        [ParameterName("wtnfactor"), ParameterInfo("Width dependence of tnfactor")]
        [Finite]
        private GivenParameter<double> _wtnfactor = new GivenParameter<double>(0.0);

        [ParameterName("wteta0"), ParameterInfo("Width dependence of teta0")]
        [Finite]
        private GivenParameter<double> _wteta0 = new GivenParameter<double>(0.0);

        [ParameterName("wtvoffcv"), ParameterInfo("Width dependence of tvoffcv")]
        [Finite]
        private GivenParameter<double> _wtvoffcv = new GivenParameter<double>(0.0);

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

        [ParameterName("pndep"), ParameterInfo("Cross-term dependence of ndep")]
        [Finite]
        public GivenParameter<double> Pndep
        {
            get => _pndep;
            set
            {
                Utility.Finite(value, nameof(Pndep));
                if (value > 1.0e20)
                    _pndep = value * 1e-6;
                else
                    _pndep = value;
            }
        }
        private GivenParameter<double> _pndep = new GivenParameter<double>(0.0);

        [ParameterName("pnsd"), ParameterInfo("Cross-term dependence of nsd")]
        [Finite]
        public GivenParameter<double> Pnsd
        {
            get => _pnsd;
            set
            {
                Utility.Finite(value, nameof(Pnsd));
                if (value > 1.0e23)
                    _pnsd = value * 1e-6;
                else
                    _pnsd = value;
            }
        }
        private GivenParameter<double> _pnsd = new GivenParameter<double>(0.0);

        [ParameterName("pphin"), ParameterInfo("Cross-term dependence of phin")]
        [Finite]
        private GivenParameter<double> _pphin = new GivenParameter<double>(0.0);

        [ParameterName("pngate"), ParameterInfo("Cross-term dependence of ngate")]
        [Finite]
        public GivenParameter<double> Pngate
        {
            get => _pngate;
            set
            {
                Utility.Finite(value, nameof(Pngate));
                if (value > 1.0e23)
                    _pngate = value * 1e-6;
                else
                    _pngate = value;
            }
        }
        private GivenParameter<double> _pngate = new GivenParameter<double>(0.0);

        [ParameterName("pgamma1"), ParameterInfo("Cross-term dependence of gamma1")]
        [Finite]
        private GivenParameter<double> _pgamma1 = new GivenParameter<double>(0.0);

        [ParameterName("pgamma2"), ParameterInfo("Cross-term dependence of gamma2")]
        [Finite]
        private GivenParameter<double> _pgamma2 = new GivenParameter<double>(0.0);

        [ParameterName("pvbx"), ParameterInfo("Cross-term dependence of vbx")]
        [Finite]
        private GivenParameter<double> _pvbx = new GivenParameter<double>(0.0);

        [ParameterName("pvbm"), ParameterInfo("Cross-term dependence of vbm")]
        [Finite]
        private GivenParameter<double> _pvbm = new GivenParameter<double>(0.0);

        [ParameterName("pxt"), ParameterInfo("Cross-term dependence of xt")]
        [Finite]
        private GivenParameter<double> _pxt = new GivenParameter<double>(0.0);

        [ParameterName("pk1"), ParameterInfo("Cross-term dependence of k1")]
        [Finite]
        private GivenParameter<double> _pk1 = new GivenParameter<double>(0.0);

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
        private GivenParameter<double> _pk2 = new GivenParameter<double>(0.0);

        [ParameterName("pk3"), ParameterInfo("Cross-term dependence of k3")]
        [Finite]
        private GivenParameter<double> _pk3 = new GivenParameter<double>(0.0);

        [ParameterName("pk3b"), ParameterInfo("Cross-term dependence of k3b")]
        [Finite]
        private GivenParameter<double> _pk3b = new GivenParameter<double>(0.0);

        [ParameterName("pw0"), ParameterInfo("Cross-term dependence of w0")]
        [Finite]
        private GivenParameter<double> _pw0 = new GivenParameter<double>(0.0);

        [ParameterName("pdvtp0"), ParameterInfo("Cross-term dependence of dvtp0")]
        [Finite]
        private GivenParameter<double> _pdvtp0 = new GivenParameter<double>(0.0);

        [ParameterName("pdvtp1"), ParameterInfo("Cross-term dependence of dvtp1")]
        [Finite]
        private GivenParameter<double> _pdvtp1 = new GivenParameter<double>(0.0);

        [ParameterName("pdvtp2"), ParameterInfo("Cross-term dependence of dvtp2")]
        [Finite]
        private GivenParameter<double> _pdvtp2 = new GivenParameter<double>(0.0);

        [ParameterName("pdvtp3"), ParameterInfo("Cross-term dependence of dvtp3")]
        [Finite]
        private GivenParameter<double> _pdvtp3 = new GivenParameter<double>(0.0);

        [ParameterName("pdvtp4"), ParameterInfo("Cross-term dependence of dvtp4")]
        [Finite]
        private GivenParameter<double> _pdvtp4 = new GivenParameter<double>(0.0);

        [ParameterName("pdvtp5"), ParameterInfo("Cross-term dependence of dvtp5")]
        [Finite]
        private GivenParameter<double> _pdvtp5 = new GivenParameter<double>(0.0);

        [ParameterName("plpe0"), ParameterInfo("Cross-term dependence of lpe0")]
        [Finite]
        private GivenParameter<double> _plpe0 = new GivenParameter<double>(0.0);

        [ParameterName("plpeb"), ParameterInfo("Cross-term dependence of lpeb")]
        [Finite]
        private GivenParameter<double> _plpeb = new GivenParameter<double>(0.0);

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

        [ParameterName("pvth0"), ParameterName("pvtho"), ParameterInfo("Cross-term dependence of vth0")]
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

        [ParameterName("pud"), ParameterInfo("Cross-term dependence of ud")]
        [Finite]
        private GivenParameter<double> _pud = new GivenParameter<double>(0.0);

        [ParameterName("pud1"), ParameterInfo("Cross-term dependence of ud1")]
        [Finite]
        private GivenParameter<double> _pud1 = new GivenParameter<double>(0.0);

        [ParameterName("pup"), ParameterInfo("Cross-term dependence of up")]
        [Finite]
        private GivenParameter<double> _pup = new GivenParameter<double>(0.0);

        [ParameterName("plp"), ParameterInfo("Cross-term dependence of lp")]
        [Finite]
        private GivenParameter<double> _plp = new GivenParameter<double>(0.0);

        [ParameterName("pu0"), ParameterInfo("Cross-term dependence of u0")]
        [Finite]
        private GivenParameter<double> _pu0 = new GivenParameter<double>(0.0);

        [ParameterName("pute"), ParameterInfo("Cross-term dependence of ute")]
        [Finite]
        private GivenParameter<double> _pute = new GivenParameter<double>(0.0);

        [ParameterName("pucste"), ParameterInfo("Cross-term dependence of ucste")]
        [Finite]
        private GivenParameter<double> _pucste = new GivenParameter<double>(0.0);

        [ParameterName("pvoff"), ParameterInfo("Cross-term dependence of voff")]
        [Finite]
        private GivenParameter<double> _pvoff = new GivenParameter<double>(0.0);

        [ParameterName("pminv"), ParameterInfo("Cross-term dependence of minv")]
        [Finite]
        private GivenParameter<double> _pminv = new GivenParameter<double>(0.0);

        [ParameterName("pminvcv"), ParameterInfo("Cross-term dependence of minvcv")]
        [Finite]
        private GivenParameter<double> _pminvcv = new GivenParameter<double>(0.0);

        [ParameterName("pdelta"), ParameterInfo("Cross-term dependence of delta")]
        [Finite]
        private GivenParameter<double> _pdelta = new GivenParameter<double>(0.0);

        [ParameterName("prdsw"), ParameterInfo("Cross-term dependence of rdsw ")]
        [Finite]
        private GivenParameter<double> _prdsw = new GivenParameter<double>(0.0);

        [ParameterName("prsw"), ParameterInfo("Cross-term dependence of rsw")]
        [Finite]
        private GivenParameter<double> _prsw = new GivenParameter<double>(0.0);

        [ParameterName("prdw"), ParameterInfo("Cross-term dependence of rdw")]
        [Finite]
        private GivenParameter<double> _prdw = new GivenParameter<double>(0.0);

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

        [ParameterName("pfprout"), ParameterInfo("Cross-term dependence of pdiblcb")]
        [Finite]
        private GivenParameter<double> _pfprout = new GivenParameter<double>(0.0);

        [ParameterName("ppdits"), ParameterInfo("Cross-term dependence of pdits")]
        [Finite]
        private GivenParameter<double> _ppdits = new GivenParameter<double>(0.0);

        [ParameterName("ppditsd"), ParameterInfo("Cross-term dependence of pditsd")]
        [Finite]
        private GivenParameter<double> _ppditsd = new GivenParameter<double>(0.0);

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

        [ParameterName("pckappas"), ParameterInfo("Cross-term dependence of ckappas")]
        [Finite]
        private GivenParameter<double> _pckappas = new GivenParameter<double>(0.0);

        [ParameterName("pckappad"), ParameterInfo("Cross-term dependence of ckappad")]
        [Finite]
        private GivenParameter<double> _pckappad = new GivenParameter<double>(0.0);

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

        [ParameterName("pagidl"), ParameterInfo("Cross-term dependence of agidl")]
        [Finite]
        private GivenParameter<double> _pagidl = new GivenParameter<double>(0.0);

        [ParameterName("pbgidl"), ParameterInfo("Cross-term dependence of bgidl")]
        [Finite]
        private GivenParameter<double> _pbgidl = new GivenParameter<double>(0.0);

        [ParameterName("pcgidl"), ParameterInfo("Cross-term dependence of cgidl")]
        [Finite]
        private GivenParameter<double> _pcgidl = new GivenParameter<double>(0.0);

        [ParameterName("prgidl"), ParameterInfo("Cross-term dependence of rgidl")]
        [Finite]
        private GivenParameter<double> _prgidl = new GivenParameter<double>(0.0);

        [ParameterName("pkgidl"), ParameterInfo("Cross-term dependence of kgidl")]
        [Finite]
        private GivenParameter<double> _pkgidl = new GivenParameter<double>(0.0);

        [ParameterName("pfgidl"), ParameterInfo("Cross-term dependence of fgidl")]
        [Finite]
        private GivenParameter<double> _pfgidl = new GivenParameter<double>(0.0);

        [ParameterName("pegidl"), ParameterInfo("Cross-term dependence of egidl")]
        [Finite]
        private GivenParameter<double> _pegidl = new GivenParameter<double>(0.0);

        [ParameterName("pagisl"), ParameterInfo("Cross-term dependence of agisl")]
        [Finite]
        private GivenParameter<double> _pagisl = new GivenParameter<double>();

        [ParameterName("pbgisl"), ParameterInfo("Cross-term dependence of bgisl")]
        [Finite]
        private GivenParameter<double> _pbgisl = new GivenParameter<double>();

        [ParameterName("pcgisl"), ParameterInfo("Cross-term dependence of cgisl")]
        [Finite]
        private GivenParameter<double> _pcgisl = new GivenParameter<double>();

        [ParameterName("pegisl"), ParameterInfo("Cross-term dependence of egisl")]
        [Finite]
        private GivenParameter<double> _pegisl = new GivenParameter<double>();

        [ParameterName("prgisl"), ParameterInfo("Cross-term dependence of rgisl")]
        [Finite]
        private GivenParameter<double> _prgisl = new GivenParameter<double>();

        [ParameterName("pkgisl"), ParameterInfo("Cross-term dependence of kgisl")]
        [Finite]
        private GivenParameter<double> _pkgisl = new GivenParameter<double>();

        [ParameterName("pfgisl"), ParameterInfo("Cross-term dependence of fgisl")]
        [Finite]
        private GivenParameter<double> _pfgisl = new GivenParameter<double>();

        [ParameterName("paigc"), ParameterInfo("Cross-term dependence of aigc")]
        [Finite]
        private GivenParameter<double> _paigc = new GivenParameter<double>(0.0);

        [ParameterName("pbigc"), ParameterInfo("Cross-term dependence of bigc")]
        [Finite]
        private GivenParameter<double> _pbigc = new GivenParameter<double>(0.0);

        [ParameterName("pcigc"), ParameterInfo("Cross-term dependence of cigc")]
        [Finite]
        private GivenParameter<double> _pcigc = new GivenParameter<double>(0.0);

        [ParameterName("paigsd"), ParameterInfo("Cross-term dependence of aigsd")]
        [Finite]
        private GivenParameter<double> _paigsd = new GivenParameter<double>();

        [ParameterName("pbigsd"), ParameterInfo("Cross-term dependence of bigsd")]
        [Finite]
        private GivenParameter<double> _pbigsd = new GivenParameter<double>();

        [ParameterName("pcigsd"), ParameterInfo("Cross-term dependence of cigsd")]
        [Finite]
        private GivenParameter<double> _pcigsd = new GivenParameter<double>();

        [ParameterName("paigs"), ParameterInfo("Cross-term dependence of aigs")]
        [Finite]
        private GivenParameter<double> _paigs = new GivenParameter<double>();

        [ParameterName("pbigs"), ParameterInfo("Cross-term dependence of bigs")]
        [Finite]
        private GivenParameter<double> _pbigs = new GivenParameter<double>();

        [ParameterName("pcigs"), ParameterInfo("Cross-term dependence of cigs")]
        [Finite]
        private GivenParameter<double> _pcigs = new GivenParameter<double>();

        [ParameterName("paigd"), ParameterInfo("Cross-term dependence of aigd")]
        [Finite]
        private GivenParameter<double> _paigd = new GivenParameter<double>(0.0);

        [ParameterName("pbigd"), ParameterInfo("Cross-term dependence of bigd")]
        [Finite]
        private GivenParameter<double> _pbigd = new GivenParameter<double>(0.0);

        [ParameterName("pcigd"), ParameterInfo("Cross-term dependence of cigd")]
        [Finite]
        private GivenParameter<double> _pcigd = new GivenParameter<double>(0.0);

        [ParameterName("paigbacc"), ParameterInfo("Cross-term dependence of aigbacc")]
        [Finite]
        private GivenParameter<double> _paigbacc = new GivenParameter<double>(0.0);

        [ParameterName("pbigbacc"), ParameterInfo("Cross-term dependence of bigbacc")]
        [Finite]
        private GivenParameter<double> _pbigbacc = new GivenParameter<double>(0.0);

        [ParameterName("pcigbacc"), ParameterInfo("Cross-term dependence of cigbacc")]
        [Finite]
        private GivenParameter<double> _pcigbacc = new GivenParameter<double>(0.0);

        [ParameterName("paigbinv"), ParameterInfo("Cross-term dependence of aigbinv")]
        [Finite]
        private GivenParameter<double> _paigbinv = new GivenParameter<double>(0.0);

        [ParameterName("pbigbinv"), ParameterInfo("Cross-term dependence of bigbinv")]
        [Finite]
        private GivenParameter<double> _pbigbinv = new GivenParameter<double>(0.0);

        [ParameterName("pcigbinv"), ParameterInfo("Cross-term dependence of cigbinv")]
        [Finite]
        private GivenParameter<double> _pcigbinv = new GivenParameter<double>(0.0);

        [ParameterName("pnigc"), ParameterInfo("Cross-term dependence of nigc")]
        [Finite]
        private GivenParameter<double> _pnigc = new GivenParameter<double>(0.0);

        [ParameterName("pnigbinv"), ParameterInfo("Cross-term dependence of nigbinv")]
        [Finite]
        private GivenParameter<double> _pnigbinv = new GivenParameter<double>(0.0);

        [ParameterName("pnigbacc"), ParameterInfo("Cross-term dependence of nigbacc")]
        [Finite]
        private GivenParameter<double> _pnigbacc = new GivenParameter<double>(0.0);

        [ParameterName("pntox"), ParameterInfo("Cross-term dependence of ntox")]
        [Finite]
        private GivenParameter<double> _pntox = new GivenParameter<double>(0.0);

        [ParameterName("peigbinv"), ParameterInfo("Cross-term dependence for eigbinv")]
        [Finite]
        private GivenParameter<double> _peigbinv = new GivenParameter<double>(0.0);

        [ParameterName("ppigcd"), ParameterInfo("Cross-term dependence for pigcd")]
        [Finite]
        private GivenParameter<double> _ppigcd = new GivenParameter<double>(0.0);

        [ParameterName("ppoxedge"), ParameterInfo("Cross-term dependence for poxedge")]
        [Finite]
        private GivenParameter<double> _ppoxedge = new GivenParameter<double>(0.0);

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

        [ParameterName("pxrcrg1"), ParameterInfo("Cross-term dependence of xrcrg1")]
        [Finite]
        private GivenParameter<double> _pxrcrg1 = new GivenParameter<double>(0.0);

        [ParameterName("pxrcrg2"), ParameterInfo("Cross-term dependence of xrcrg2")]
        [Finite]
        private GivenParameter<double> _pxrcrg2 = new GivenParameter<double>(0.0);

        [ParameterName("plambda"), ParameterInfo("Cross-term dependence of lambda")]
        [Finite]
        private GivenParameter<double> _plambda = new GivenParameter<double>(0.0);

        [ParameterName("pvtl"), ParameterInfo("Cross-term dependence of vtl")]
        [Finite]
        private GivenParameter<double> _pvtl = new GivenParameter<double>(0.0);

        [ParameterName("pxn"), ParameterInfo("Cross-term dependence of xn")]
        [Finite]
        private GivenParameter<double> _pxn = new GivenParameter<double>(0.0);

        [ParameterName("peu"), ParameterInfo("Cross-term dependence of eu")]
        [Finite]
        private GivenParameter<double> _peu = new GivenParameter<double>(0.0);

        [ParameterName("pucs"), ParameterInfo("Cross-term dependence of ucs")]
        [Finite]
        private GivenParameter<double> _pucs = new GivenParameter<double>(0.0);

        [ParameterName("pvfbsdoff"), ParameterInfo("Cross-term dependence of vfbsdoff")]
        [Finite]
        private GivenParameter<double> _pvfbsdoff = new GivenParameter<double>(0.0);

        [ParameterName("ptvfbsdoff"), ParameterInfo("Cross-term dependence of tvfbsdoff")]
        [Finite]
        private GivenParameter<double> _ptvfbsdoff = new GivenParameter<double>(0.0);

        [ParameterName("ptvoff"), ParameterInfo("Cross-term dependence of tvoff")]
        [Finite]
        private GivenParameter<double> _ptvoff = new GivenParameter<double>(0.0);

        [ParameterName("ptnfactor"), ParameterInfo("Cross-term dependence of tnfactor")]
        [Finite]
        private GivenParameter<double> _ptnfactor = new GivenParameter<double>(0.0);

        [ParameterName("pteta0"), ParameterInfo("Cross-term dependence of teta0")]
        [Finite]
        private GivenParameter<double> _pteta0 = new GivenParameter<double>(0.0);

        [ParameterName("ptvoffcv"), ParameterInfo("Cross-term dependence of tvoffcv")]
        [Finite]
        private GivenParameter<double> _ptvoffcv = new GivenParameter<double>(0.0);

        [ParameterName("saref"), ParameterInfo("Reference distance between OD edge to poly of one side")]
        [Finite]
        private GivenParameter<double> _saref = new GivenParameter<double>(1e-6);

        [ParameterName("sbref"), ParameterInfo("Reference distance between OD edge to poly of the other side")]
        [Finite]
        private GivenParameter<double> _sbref = new GivenParameter<double>(1e-6);

        [ParameterName("wlod"), ParameterInfo("Width parameter for stress effect")]
        [Finite]
        private GivenParameter<double> _wlod = new GivenParameter<double>(0);

        [ParameterName("ku0"), ParameterInfo("Mobility degradation/enhancement coefficient for LOD")]
        [Finite]
        private GivenParameter<double> _ku0 = new GivenParameter<double>(0);

        [ParameterName("kvsat"), ParameterInfo("Saturation velocity degradation/enhancement parameter for LOD")]
        [Finite]
        private GivenParameter<double> _kvsat = new GivenParameter<double>(0);

        [ParameterName("kvth0"), ParameterInfo("Threshold degradation/enhancement parameter for LOD")]
        [Finite]
        private GivenParameter<double> _kvth0 = new GivenParameter<double>(0);

        [ParameterName("tku0"), ParameterInfo("Temperature coefficient of KU0")]
        [Finite]
        private GivenParameter<double> _tku0 = new GivenParameter<double>(0);

        [ParameterName("llodku0"), ParameterInfo("Length parameter for u0 LOD effect")]
        [Finite]
        private GivenParameter<double> _llodku0 = new GivenParameter<double>(0);

        [ParameterName("wlodku0"), ParameterInfo("Width parameter for u0 LOD effect")]
        [Finite]
        private GivenParameter<double> _wlodku0 = new GivenParameter<double>(0);

        [ParameterName("llodvth"), ParameterInfo("Length parameter for vth LOD effect")]
        [Finite]
        private GivenParameter<double> _llodvth = new GivenParameter<double>(0);

        [ParameterName("wlodvth"), ParameterInfo("Width parameter for vth LOD effect")]
        [Finite]
        private GivenParameter<double> _wlodvth = new GivenParameter<double>(0);

        [ParameterName("lku0"), ParameterInfo("Length dependence of ku0")]
        [Finite]
        private GivenParameter<double> _lku0 = new GivenParameter<double>(0);

        [ParameterName("wku0"), ParameterInfo("Width dependence of ku0")]
        [Finite]
        private GivenParameter<double> _wku0 = new GivenParameter<double>(0);

        [ParameterName("pku0"), ParameterInfo("Cross-term dependence of ku0")]
        [Finite]
        private GivenParameter<double> _pku0 = new GivenParameter<double>(0);

        [ParameterName("lkvth0"), ParameterInfo("Length dependence of kvth0")]
        [Finite]
        private GivenParameter<double> _lkvth0 = new GivenParameter<double>(0);

        [ParameterName("wkvth0"), ParameterInfo("Width dependence of kvth0")]
        [Finite]
        private GivenParameter<double> _wkvth0 = new GivenParameter<double>(0);

        [ParameterName("pkvth0"), ParameterInfo("Cross-term dependence of kvth0")]
        [Finite]
        private GivenParameter<double> _pkvth0 = new GivenParameter<double>(0);

        [ParameterName("stk2"), ParameterInfo("K2 shift factor related to stress effect on vth")]
        [Finite]
        private GivenParameter<double> _stk2 = new GivenParameter<double>(0);

        [ParameterName("lodk2"), ParameterInfo("K2 shift modification factor for stress effect")]
        [Finite]
        private GivenParameter<double> _lodk2 = new GivenParameter<double>(1.0);

        [ParameterName("steta0"), ParameterInfo("eta0 shift factor related to stress effect on vth")]
        [Finite]
        private GivenParameter<double> _steta0 = new GivenParameter<double>(0);

        [ParameterName("lodeta0"), ParameterInfo("eta0 shift modification factor for stress effect")]
        [Finite]
        private GivenParameter<double> _lodeta0 = new GivenParameter<double>(1.0);

        [ParameterName("web"), ParameterInfo("Coefficient for SCB")]
        [Finite]
        private GivenParameter<double> _web = new GivenParameter<double>(0.0);

        [ParameterName("wec"), ParameterInfo("Coefficient for SCC")]
        [Finite]
        private GivenParameter<double> _wec = new GivenParameter<double>(0.0);

        [ParameterName("kvth0we"), ParameterInfo("Threshold shift factor for well proximity effect")]
        [Finite]
        private GivenParameter<double> _kvth0we = new GivenParameter<double>(0.0);

        [ParameterName("k2we"), ParameterInfo(" K2 shift factor for well proximity effect ")]
        [Finite]
        private GivenParameter<double> _k2we = new GivenParameter<double>(0.0);

        [ParameterName("ku0we"), ParameterInfo(" Mobility degradation factor for well proximity effect ")]
        [Finite]
        private GivenParameter<double> _ku0we = new GivenParameter<double>(0.0);

        [ParameterName("scref"), ParameterInfo(" Reference distance to calculate SCA, SCB and SCC")]
        [Finite]
        private GivenParameter<double> _scref = new GivenParameter<double>(1.0E-6);

        [ParameterName("wpemod"), ParameterInfo(" Flag for WPE model (WPEMOD=1 to activate this model) ")]
        [Finite]
        private GivenParameter<double> _wpemod = new GivenParameter<double>();

        [ParameterName("lkvth0we"), ParameterInfo("Length dependence of kvth0we")]
        [Finite]
        private GivenParameter<double> _lkvth0we = new GivenParameter<double>(0);

        [ParameterName("lk2we"), ParameterInfo(" Length dependence of k2we ")]
        [Finite]
        private GivenParameter<double> _lk2we = new GivenParameter<double>(0);

        [ParameterName("lku0we"), ParameterInfo(" Length dependence of ku0we ")]
        [Finite]
        private GivenParameter<double> _lku0we = new GivenParameter<double>(0);

        [ParameterName("wkvth0we"), ParameterInfo("Width dependence of kvth0we")]
        [Finite]
        private GivenParameter<double> _wkvth0we = new GivenParameter<double>(0);

        [ParameterName("wk2we"), ParameterInfo(" Width dependence of k2we ")]
        [Finite]
        private GivenParameter<double> _wk2we = new GivenParameter<double>(0);

        [ParameterName("wku0we"), ParameterInfo(" Width dependence of ku0we ")]
        [Finite]
        private GivenParameter<double> _wku0we = new GivenParameter<double>(0);

        [ParameterName("pkvth0we"), ParameterInfo("Cross-term dependence of kvth0we")]
        [Finite]
        private GivenParameter<double> _pkvth0we = new GivenParameter<double>(0);

        [ParameterName("pk2we"), ParameterInfo(" Cross-term dependence of k2we ")]
        [Finite]
        private GivenParameter<double> _pk2we = new GivenParameter<double>(0);

        [ParameterName("pku0we"), ParameterInfo(" Cross-term dependence of ku0we ")]
        [Finite]
        private GivenParameter<double> _pku0we = new GivenParameter<double>(0);

        [ParameterName("noia"), ParameterInfo("Flicker noise parameter")]
        [Finite]
        private GivenParameter<double> _oxideTrapDensityA = new GivenParameter<double>();

        [ParameterName("noib"), ParameterInfo("Flicker noise parameter")]
        [Finite]
        private GivenParameter<double> _oxideTrapDensityB = new GivenParameter<double>();

        [ParameterName("noic"), ParameterInfo("Flicker noise parameter")]
        [Finite]
        private GivenParameter<double> _oxideTrapDensityC = new GivenParameter<double>(8.75e9);

        [ParameterName("tnoia"), ParameterInfo("Thermal noise parameter")]
        [Finite]
        private GivenParameter<double> _tnoia = new GivenParameter<double>(1.5);

        [ParameterName("tnoib"), ParameterInfo("Thermal noise parameter")]
        [Finite]
        private GivenParameter<double> _tnoib = new GivenParameter<double>(3.5);

        [ParameterName("tnoic"), ParameterInfo("Thermal noise parameter")]
        [Finite]
        private GivenParameter<double> _tnoic = new GivenParameter<double>(0.0);

        [ParameterName("rnoia"), ParameterInfo("Thermal noise coefficient")]
        [Finite]
        private GivenParameter<double> _rnoia = new GivenParameter<double>(0.577);

        [ParameterName("rnoib"), ParameterInfo("Thermal noise coefficient")]
        [Finite]
        private GivenParameter<double> _rnoib = new GivenParameter<double>(0.5164);

        [ParameterName("rnoic"), ParameterInfo("Thermal noise coefficient")]
        [Finite]
        private GivenParameter<double> _rnoic = new GivenParameter<double>(0.395);

        [ParameterName("ntnoi"), ParameterInfo("Thermal noise parameter")]
        [Finite]
        private GivenParameter<double> _ntnoi = new GivenParameter<double>(1.0);

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

        /// <summary>
        /// Gets the mosfet type (1 for NMOS, -1 for PMOS).
        /// </summary>
        public double Type { get; private set; } = 1.0;

        [ParameterName("tnom"), ParameterInfo("Parameter measurement temperature")]
        [DerivedProperty, Finite, GreaterThan(-Constants.CelsiusKelvin)]
        public double TnomCelsius
        {
            get => Tnom - Constants.CelsiusKelvin;
            set => Tnom = value + Constants.CelsiusKelvin;
        }

        [GreaterThan(0), Finite]
        private GivenParameter<double> _tnom = new GivenParameter<double>(300.15);

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
