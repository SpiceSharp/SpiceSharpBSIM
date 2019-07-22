using System;
using SpiceSharp.Attributes;
namespace SpiceSharp.Components.BSIM3v24Behaviors
{
	
	/// <summary>
	/// Base parameters for a <see cref="BSIM3v24Model" />
	/// </summary>
	public class ModelBaseParameters : ParameterSet
	{
	    /// <summary>
	    /// Gets or sets the path for parameter checking
	    /// </summary>
	    [ParameterName("path"), ParameterInfo("File path for used to log parameter checks")]
	    public string CheckPath { get; set; } = "bsim3check.log";
		
		/// <summary>
		/// Properties
		/// </summary>
		[ParameterName("mobmod"), ParameterInfo("Mobility model selector")]
		public GivenParameter<int> MobMod { get; } = new GivenParameter<int>(1);
		[ParameterName("binunit"), ParameterInfo("Bin  unit  selector")]
		public GivenParameter<int> BinUnit { get; } = new GivenParameter<int>(1);
		[ParameterName("paramchk"), ParameterInfo("Model parameter checking selector")]
		public GivenParameter<int> ParamChk { get; } = new GivenParameter<int>();
		[ParameterName("capmod"), ParameterInfo("Capacitance model selector")]
		public GivenParameter<int> CapMod { get; } = new GivenParameter<int>(3);
		[ParameterName("noimod"), ParameterInfo("Noise model selector")]
		public GivenParameter<int> NoiMod { get; } = new GivenParameter<int>(1);
	    [ParameterName("version"), ParameterInfo(" parameter for model version")]
	    public string Version { get; set; } = "3.2.4";
		[ParameterName("tox"), ParameterInfo("Gate oxide thickness in meters")]
		public GivenParameter<double> Tox { get; } = new GivenParameter<double>(1.5e-08);
		[ParameterName("toxm"), ParameterInfo("Gate oxide thickness used in extraction")]
		public GivenParameter<double> Toxm { get; } = new GivenParameter<double>();
		[ParameterName("cdsc"), ParameterInfo("Drain/Source and channel coupling capacitance")]
		public GivenParameter<double> Cdsc { get; } = new GivenParameter<double>(0.00024);
		[ParameterName("cdscb"), ParameterInfo("Body-bias dependence of cdsc")]
		public GivenParameter<double> Cdscb { get; } = new GivenParameter<double>();
		[ParameterName("cdscd"), ParameterInfo("Drain-bias dependence of cdsc")]
		public GivenParameter<double> Cdscd { get; } = new GivenParameter<double>();
		[ParameterName("cit"), ParameterInfo("Interface state capacitance")]
		public GivenParameter<double> Cit { get; } = new GivenParameter<double>();
		[ParameterName("nfactor"), ParameterInfo("Subthreshold swing Coefficient")]
		public GivenParameter<double> Nfactor { get; } = new GivenParameter<double>(1);
		[ParameterName("xj"), ParameterInfo("Junction depth in meters")]
		public GivenParameter<double> Xj { get; } = new GivenParameter<double>(1.5e-07);
		[ParameterName("vsat"), ParameterInfo("Saturation velocity at tnom")]
		public GivenParameter<double> Vsat { get; } = new GivenParameter<double>(80000);
		[ParameterName("a0"), ParameterInfo("Non-uniform depletion width effect coefficient.")]
		public GivenParameter<double> A0 { get; } = new GivenParameter<double>(1);
		[ParameterName("ags"), ParameterInfo("Gate bias  coefficient of Abulk.")]
		public GivenParameter<double> Ags { get; } = new GivenParameter<double>();
		[ParameterName("a1"), ParameterInfo("Non-saturation effect coefficient")]
		public GivenParameter<double> A1 { get; } = new GivenParameter<double>();
		[ParameterName("a2"), ParameterInfo("Non-saturation effect coefficient")]
		public GivenParameter<double> A2 { get; } = new GivenParameter<double>(1);
		[ParameterName("at"), ParameterInfo("Temperature coefficient of vsat")]
		public GivenParameter<double> At { get; } = new GivenParameter<double>(33000);
		[ParameterName("keta"), ParameterInfo("Body-bias coefficient of non-uniform depletion width effect.")]
		public GivenParameter<double> Keta { get; } = new GivenParameter<double>(-0.047);
		[ParameterName("nsub"), ParameterInfo("Substrate doping concentration")]
		public GivenParameter<double> Nsub { get; } = new GivenParameter<double>(6e+16);
		[ParameterName("nch"), ParameterInfo("Channel doping concentration")]
		public GivenParameter<double> Npeak { get; } = new GivenParameter<double>(1.7e+17);
		[ParameterName("ngate"), ParameterInfo("Poly-gate doping concentration")]
		public GivenParameter<double> Ngate { get; } = new GivenParameter<double>();
		[ParameterName("gamma1"), ParameterInfo("Vth body coefficient")]
		public GivenParameter<double> Gamma1 { get; } = new GivenParameter<double>();
		[ParameterName("gamma2"), ParameterInfo("Vth body coefficient")]
		public GivenParameter<double> Gamma2 { get; } = new GivenParameter<double>();
		[ParameterName("vbx"), ParameterInfo("Vth transition body Voltage")]
		public GivenParameter<double> Vbx { get; } = new GivenParameter<double>();
		[ParameterName("vbm"), ParameterInfo("Maximum body voltage")]
		public GivenParameter<double> Vbm { get; } = new GivenParameter<double>(-3);
		[ParameterName("xt"), ParameterInfo("Doping depth")]
		public GivenParameter<double> Xt { get; } = new GivenParameter<double>(1.55e-07);
		[ParameterName("k1"), ParameterInfo("Bulk effect coefficient 1")]
		public GivenParameter<double> K1 { get; } = new GivenParameter<double>();
		[ParameterName("kt1"), ParameterInfo("Temperature coefficient of Vth")]
		public GivenParameter<double> Kt1 { get; } = new GivenParameter<double>(-0.11);
		[ParameterName("kt1l"), ParameterInfo("Temperature coefficient of Vth")]
		public GivenParameter<double> Kt1l { get; } = new GivenParameter<double>();
		[ParameterName("kt2"), ParameterInfo("Body-coefficient of kt1")]
		public GivenParameter<double> Kt2 { get; } = new GivenParameter<double>(0.022);
		[ParameterName("k2"), ParameterInfo("Bulk effect coefficient 2")]
		public GivenParameter<double> K2 { get; } = new GivenParameter<double>();
		[ParameterName("k3"), ParameterInfo("Narrow width effect coefficient")]
		public GivenParameter<double> K3 { get; } = new GivenParameter<double>(80);
		[ParameterName("k3b"), ParameterInfo("Body effect coefficient of k3")]
		public GivenParameter<double> K3b { get; } = new GivenParameter<double>();
		[ParameterName("nlx"), ParameterInfo("Lateral non-uniform doping effect")]
		public GivenParameter<double> Nlx { get; } = new GivenParameter<double>(1.74e-07);
		[ParameterName("w0"), ParameterInfo("Narrow width effect parameter")]
		public GivenParameter<double> W0 { get; } = new GivenParameter<double>(2.5e-06);
		[ParameterName("dvt0"), ParameterInfo("Short channel effect coeff. 0")]
		public GivenParameter<double> Dvt0 { get; } = new GivenParameter<double>(2.2);
		[ParameterName("dvt1"), ParameterInfo("Short channel effect coeff. 1")]
		public GivenParameter<double> Dvt1 { get; } = new GivenParameter<double>(0.53);
		[ParameterName("dvt2"), ParameterInfo("Short channel effect coeff. 2")]
		public GivenParameter<double> Dvt2 { get; } = new GivenParameter<double>(-0.032);
		[ParameterName("dvt0w"), ParameterInfo("Narrow Width coeff. 0")]
		public GivenParameter<double> Dvt0w { get; } = new GivenParameter<double>();
		[ParameterName("dvt1w"), ParameterInfo("Narrow Width effect coeff. 1")]
		public GivenParameter<double> Dvt1w { get; } = new GivenParameter<double>(5300000);
		[ParameterName("dvt2w"), ParameterInfo("Narrow Width effect coeff. 2")]
		public GivenParameter<double> Dvt2w { get; } = new GivenParameter<double>(-0.032);
		[ParameterName("drout"), ParameterInfo("DIBL coefficient of output resistance")]
		public GivenParameter<double> Drout { get; } = new GivenParameter<double>(0.56);
		[ParameterName("dsub"), ParameterInfo("DIBL coefficient in the subthreshold region")]
		public GivenParameter<double> Dsub { get; } = new GivenParameter<double>();
		[ParameterName("vth0"), ParameterName("vtho"), ParameterInfo("Threshold voltage")]
		public GivenParameter<double> Vth0 { get; } = new GivenParameter<double>();
		[ParameterName("ua"), ParameterInfo("Linear gate dependence of mobility")]
		public GivenParameter<double> Ua { get; } = new GivenParameter<double>(2.25e-09);
		[ParameterName("ua1"), ParameterInfo("Temperature coefficient of ua")]
		public GivenParameter<double> Ua1 { get; } = new GivenParameter<double>(4.31e-09);
		[ParameterName("ub"), ParameterInfo("Quadratic gate dependence of mobility")]
		public GivenParameter<double> Ub { get; } = new GivenParameter<double>(5.87e-19);
		[ParameterName("ub1"), ParameterInfo("Temperature coefficient of ub")]
		public GivenParameter<double> Ub1 { get; } = new GivenParameter<double>(-7.61e-18);
		[ParameterName("uc"), ParameterInfo("Body-bias dependence of mobility")]
		public GivenParameter<double> Uc { get; } = new GivenParameter<double>();
		[ParameterName("uc1"), ParameterInfo("Temperature coefficient of uc")]
		public GivenParameter<double> Uc1 { get; } = new GivenParameter<double>();
		[ParameterName("u0"), ParameterInfo("Low-field mobility at Tnom")]
		public GivenParameter<double> U0 { get; } = new GivenParameter<double>();
		[ParameterName("ute"), ParameterInfo("Temperature coefficient of mobility")]
		public GivenParameter<double> Ute { get; } = new GivenParameter<double>(-1.5);
		[ParameterName("voff"), ParameterInfo("Threshold voltage offset")]
		public GivenParameter<double> Voff { get; } = new GivenParameter<double>(-0.08);
		[ParameterName("delta"), ParameterInfo("Effective Vds parameter")]
		public GivenParameter<double> Delta { get; } = new GivenParameter<double>(0.01);
		[ParameterName("rdsw"), ParameterInfo("Source-drain resistance per width")]
		public GivenParameter<double> Rdsw { get; } = new GivenParameter<double>();
		[ParameterName("prwg"), ParameterInfo("Gate-bias effect on parasitic resistance ")]
		public GivenParameter<double> Prwg { get; } = new GivenParameter<double>();
		[ParameterName("prwb"), ParameterInfo("Body-effect on parasitic resistance ")]
		public GivenParameter<double> Prwb { get; } = new GivenParameter<double>();
		[ParameterName("prt"), ParameterInfo("Temperature coefficient of parasitic resistance ")]
		public GivenParameter<double> Prt { get; } = new GivenParameter<double>();
		[ParameterName("eta0"), ParameterInfo("Subthreshold region DIBL coefficient")]
		public GivenParameter<double> Eta0 { get; } = new GivenParameter<double>(0.08);
		[ParameterName("etab"), ParameterInfo("Subthreshold region DIBL coefficient")]
		public GivenParameter<double> Etab { get; } = new GivenParameter<double>(-0.07);
		[ParameterName("pclm"), ParameterInfo("Channel length modulation Coefficient")]
		public GivenParameter<double> Pclm { get; } = new GivenParameter<double>(1.3);
		[ParameterName("pdiblc1"), ParameterInfo("Drain-induced barrier lowering coefficient")]
		public GivenParameter<double> Pdibl1 { get; } = new GivenParameter<double>(0.39);
		[ParameterName("pdiblc2"), ParameterInfo("Drain-induced barrier lowering coefficient")]
		public GivenParameter<double> Pdibl2 { get; } = new GivenParameter<double>(0.0086);
		[ParameterName("pdiblcb"), ParameterInfo("Body-effect on drain-induced barrier lowering")]
		public GivenParameter<double> Pdiblb { get; } = new GivenParameter<double>();
		[ParameterName("pscbe1"), ParameterInfo("Substrate current body-effect coefficient")]
		public GivenParameter<double> Pscbe1 { get; } = new GivenParameter<double>(424000000);
		[ParameterName("pscbe2"), ParameterInfo("Substrate current body-effect coefficient")]
		public GivenParameter<double> Pscbe2 { get; } = new GivenParameter<double>(1e-05);
		[ParameterName("pvag"), ParameterInfo("Gate dependence of output resistance parameter")]
		public GivenParameter<double> Pvag { get; } = new GivenParameter<double>();
		[ParameterName("wr"), ParameterInfo("Width dependence of rds")]
		public GivenParameter<double> Wr { get; } = new GivenParameter<double>(1);
		[ParameterName("dwg"), ParameterInfo("Width reduction parameter")]
		public GivenParameter<double> Dwg { get; } = new GivenParameter<double>();
		[ParameterName("dwb"), ParameterInfo("Width reduction parameter")]
		public GivenParameter<double> Dwb { get; } = new GivenParameter<double>();
		[ParameterName("b0"), ParameterInfo("Abulk narrow width parameter")]
		public GivenParameter<double> B0 { get; } = new GivenParameter<double>();
		[ParameterName("b1"), ParameterInfo("Abulk narrow width parameter")]
		public GivenParameter<double> B1 { get; } = new GivenParameter<double>();
		[ParameterName("alpha0"), ParameterInfo("substrate current model parameter")]
		public GivenParameter<double> Alpha0 { get; } = new GivenParameter<double>();
		[ParameterName("alpha1"), ParameterInfo("substrate current model parameter")]
		public GivenParameter<double> Alpha1 { get; } = new GivenParameter<double>();
		[ParameterName("beta0"), ParameterInfo("substrate current model parameter")]
		public GivenParameter<double> Beta0 { get; } = new GivenParameter<double>(30);
		[ParameterName("ijth"), ParameterInfo("Diode limiting current")]
		public GivenParameter<double> Ijth { get; } = new GivenParameter<double>(0.1);
		[ParameterName("vfb"), ParameterInfo("Flat Band Voltage")]
		public GivenParameter<double> Vfb { get; } = new GivenParameter<double>();
		[ParameterName("elm"), ParameterInfo("Non-quasi-static Elmore Constant Parameter")]
		public GivenParameter<double> Elm { get; } = new GivenParameter<double>(5);
		[ParameterName("cgsl"), ParameterInfo("New C-V model parameter")]
		public GivenParameter<double> Cgsl { get; } = new GivenParameter<double>();
		[ParameterName("cgdl"), ParameterInfo("New C-V model parameter")]
		public GivenParameter<double> Cgdl { get; } = new GivenParameter<double>();
		[ParameterName("ckappa"), ParameterInfo("New C-V model parameter")]
		public GivenParameter<double> Ckappa { get; } = new GivenParameter<double>(0.6);
		[ParameterName("cf"), ParameterInfo("Fringe capacitance parameter")]
		public GivenParameter<double> Cf { get; } = new GivenParameter<double>();
		[ParameterName("clc"), ParameterInfo("Vdsat parameter for C-V model")]
		public GivenParameter<double> Clc { get; } = new GivenParameter<double>(1e-07);
		[ParameterName("cle"), ParameterInfo("Vdsat parameter for C-V model")]
		public GivenParameter<double> Cle { get; } = new GivenParameter<double>(0.6);
		[ParameterName("dwc"), ParameterInfo("Delta W for C-V model")]
		public GivenParameter<double> Dwc { get; } = new GivenParameter<double>();
		[ParameterName("dlc"), ParameterInfo("Delta L for C-V model")]
		public GivenParameter<double> Dlc { get; } = new GivenParameter<double>();
		[ParameterName("vfbcv"), ParameterInfo("Flat Band Voltage parameter for capmod=0 only")]
		public GivenParameter<double> Vfbcv { get; } = new GivenParameter<double>(-1);
		[ParameterName("acde"), ParameterInfo("Exponential coefficient for finite charge thickness")]
		public GivenParameter<double> Acde { get; } = new GivenParameter<double>(1);
		[ParameterName("moin"), ParameterInfo("Coefficient for gate-bias dependent surface potential")]
		public GivenParameter<double> Moin { get; } = new GivenParameter<double>(15);
		[ParameterName("noff"), ParameterInfo("C-V turn-on/off parameter")]
		public GivenParameter<double> Noff { get; } = new GivenParameter<double>(1);
		[ParameterName("voffcv"), ParameterInfo("C-V lateral-shift parameter")]
		public GivenParameter<double> Voffcv { get; } = new GivenParameter<double>();
		[ParameterName("tcj"), ParameterInfo("Temperature coefficient of cj")]
		public GivenParameter<double> Tcj { get; } = new GivenParameter<double>();
		[ParameterName("tpb"), ParameterInfo("Temperature coefficient of pb")]
		public GivenParameter<double> Tpb { get; } = new GivenParameter<double>();
		[ParameterName("tcjsw"), ParameterInfo("Temperature coefficient of cjsw")]
		public GivenParameter<double> Tcjsw { get; } = new GivenParameter<double>();
		[ParameterName("tpbsw"), ParameterInfo("Temperature coefficient of pbsw")]
		public GivenParameter<double> Tpbsw { get; } = new GivenParameter<double>();
		[ParameterName("tcjswg"), ParameterInfo("Temperature coefficient of cjswg")]
		public GivenParameter<double> Tcjswg { get; } = new GivenParameter<double>();
		[ParameterName("tpbswg"), ParameterInfo("Temperature coefficient of pbswg")]
		public GivenParameter<double> Tpbswg { get; } = new GivenParameter<double>();
		[ParameterName("lcdsc"), ParameterInfo("Length dependence of cdsc")]
		public GivenParameter<double> Lcdsc { get; } = new GivenParameter<double>();
		[ParameterName("lcdscb"), ParameterInfo("Length dependence of cdscb")]
		public GivenParameter<double> Lcdscb { get; } = new GivenParameter<double>();
		[ParameterName("lcdscd"), ParameterInfo("Length dependence of cdscd")]
		public GivenParameter<double> Lcdscd { get; } = new GivenParameter<double>();
		[ParameterName("lcit"), ParameterInfo("Length dependence of cit")]
		public GivenParameter<double> Lcit { get; } = new GivenParameter<double>();
		[ParameterName("lnfactor"), ParameterInfo("Length dependence of nfactor")]
		public GivenParameter<double> Lnfactor { get; } = new GivenParameter<double>();
		[ParameterName("lxj"), ParameterInfo("Length dependence of xj")]
		public GivenParameter<double> Lxj { get; } = new GivenParameter<double>();
		[ParameterName("lvsat"), ParameterInfo("Length dependence of vsat")]
		public GivenParameter<double> Lvsat { get; } = new GivenParameter<double>();
		[ParameterName("la0"), ParameterInfo("Length dependence of a0")]
		public GivenParameter<double> La0 { get; } = new GivenParameter<double>();
		[ParameterName("lags"), ParameterInfo("Length dependence of ags")]
		public GivenParameter<double> Lags { get; } = new GivenParameter<double>();
		[ParameterName("la1"), ParameterInfo("Length dependence of a1")]
		public GivenParameter<double> La1 { get; } = new GivenParameter<double>();
		[ParameterName("la2"), ParameterInfo("Length dependence of a2")]
		public GivenParameter<double> La2 { get; } = new GivenParameter<double>();
		[ParameterName("lat"), ParameterInfo("Length dependence of at")]
		public GivenParameter<double> Lat { get; } = new GivenParameter<double>();
		[ParameterName("lketa"), ParameterInfo("Length dependence of keta")]
		public GivenParameter<double> Lketa { get; } = new GivenParameter<double>();
		[ParameterName("lnsub"), ParameterInfo("Length dependence of nsub")]
		public GivenParameter<double> Lnsub { get; } = new GivenParameter<double>();
		[ParameterName("lnch"), ParameterInfo("Length dependence of nch")]
		public GivenParameter<double> Lnpeak { get; } = new GivenParameter<double>();
		[ParameterName("lngate"), ParameterInfo("Length dependence of ngate")]
		public GivenParameter<double> Lngate { get; } = new GivenParameter<double>();
		[ParameterName("lgamma1"), ParameterInfo("Length dependence of gamma1")]
		public GivenParameter<double> Lgamma1 { get; } = new GivenParameter<double>();
		[ParameterName("lgamma2"), ParameterInfo("Length dependence of gamma2")]
		public GivenParameter<double> Lgamma2 { get; } = new GivenParameter<double>();
		[ParameterName("lvbx"), ParameterInfo("Length dependence of vbx")]
		public GivenParameter<double> Lvbx { get; } = new GivenParameter<double>();
		[ParameterName("lvbm"), ParameterInfo("Length dependence of vbm")]
		public GivenParameter<double> Lvbm { get; } = new GivenParameter<double>();
		[ParameterName("lxt"), ParameterInfo("Length dependence of xt")]
		public GivenParameter<double> Lxt { get; } = new GivenParameter<double>();
		[ParameterName("lk1"), ParameterInfo("Length dependence of k1")]
		public GivenParameter<double> Lk1 { get; } = new GivenParameter<double>();
		[ParameterName("lkt1"), ParameterInfo("Length dependence of kt1")]
		public GivenParameter<double> Lkt1 { get; } = new GivenParameter<double>();
		[ParameterName("lkt1l"), ParameterInfo("Length dependence of kt1l")]
		public GivenParameter<double> Lkt1l { get; } = new GivenParameter<double>();
		[ParameterName("lkt2"), ParameterInfo("Length dependence of kt2")]
		public GivenParameter<double> Lkt2 { get; } = new GivenParameter<double>();
		[ParameterName("lk2"), ParameterInfo("Length dependence of k2")]
		public GivenParameter<double> Lk2 { get; } = new GivenParameter<double>();
		[ParameterName("lk3"), ParameterInfo("Length dependence of k3")]
		public GivenParameter<double> Lk3 { get; } = new GivenParameter<double>();
		[ParameterName("lk3b"), ParameterInfo("Length dependence of k3b")]
		public GivenParameter<double> Lk3b { get; } = new GivenParameter<double>();
		[ParameterName("lnlx"), ParameterInfo("Length dependence of nlx")]
		public GivenParameter<double> Lnlx { get; } = new GivenParameter<double>();
		[ParameterName("lw0"), ParameterInfo("Length dependence of w0")]
		public GivenParameter<double> Lw0 { get; } = new GivenParameter<double>();
		[ParameterName("ldvt0"), ParameterInfo("Length dependence of dvt0")]
		public GivenParameter<double> Ldvt0 { get; } = new GivenParameter<double>();
		[ParameterName("ldvt1"), ParameterInfo("Length dependence of dvt1")]
		public GivenParameter<double> Ldvt1 { get; } = new GivenParameter<double>();
		[ParameterName("ldvt2"), ParameterInfo("Length dependence of dvt2")]
		public GivenParameter<double> Ldvt2 { get; } = new GivenParameter<double>();
		[ParameterName("ldvt0w"), ParameterInfo("Length dependence of dvt0w")]
		public GivenParameter<double> Ldvt0w { get; } = new GivenParameter<double>();
		[ParameterName("ldvt1w"), ParameterInfo("Length dependence of dvt1w")]
		public GivenParameter<double> Ldvt1w { get; } = new GivenParameter<double>();
		[ParameterName("ldvt2w"), ParameterInfo("Length dependence of dvt2w")]
		public GivenParameter<double> Ldvt2w { get; } = new GivenParameter<double>();
		[ParameterName("ldrout"), ParameterInfo("Length dependence of drout")]
		public GivenParameter<double> Ldrout { get; } = new GivenParameter<double>();
		[ParameterName("ldsub"), ParameterInfo("Length dependence of dsub")]
		public GivenParameter<double> Ldsub { get; } = new GivenParameter<double>();
		[ParameterName("lvth0"), ParameterName("lvtho"), ParameterInfo("Length dependence of vto")]
		public GivenParameter<double> Lvth0 { get; } = new GivenParameter<double>();
		[ParameterName("lua"), ParameterInfo("Length dependence of ua")]
		public GivenParameter<double> Lua { get; } = new GivenParameter<double>();
		[ParameterName("lua1"), ParameterInfo("Length dependence of ua1")]
		public GivenParameter<double> Lua1 { get; } = new GivenParameter<double>();
		[ParameterName("lub"), ParameterInfo("Length dependence of ub")]
		public GivenParameter<double> Lub { get; } = new GivenParameter<double>();
		[ParameterName("lub1"), ParameterInfo("Length dependence of ub1")]
		public GivenParameter<double> Lub1 { get; } = new GivenParameter<double>();
		[ParameterName("luc"), ParameterInfo("Length dependence of uc")]
		public GivenParameter<double> Luc { get; } = new GivenParameter<double>();
		[ParameterName("luc1"), ParameterInfo("Length dependence of uc1")]
		public GivenParameter<double> Luc1 { get; } = new GivenParameter<double>();
		[ParameterName("lu0"), ParameterInfo("Length dependence of u0")]
		public GivenParameter<double> Lu0 { get; } = new GivenParameter<double>();
		[ParameterName("lute"), ParameterInfo("Length dependence of ute")]
		public GivenParameter<double> Lute { get; } = new GivenParameter<double>();
		[ParameterName("lvoff"), ParameterInfo("Length dependence of voff")]
		public GivenParameter<double> Lvoff { get; } = new GivenParameter<double>();
		[ParameterName("ldelta"), ParameterInfo("Length dependence of delta")]
		public GivenParameter<double> Ldelta { get; } = new GivenParameter<double>();
		[ParameterName("lrdsw"), ParameterInfo("Length dependence of rdsw ")]
		public GivenParameter<double> Lrdsw { get; } = new GivenParameter<double>();
		[ParameterName("lprwb"), ParameterInfo("Length dependence of prwb ")]
		public GivenParameter<double> Lprwb { get; } = new GivenParameter<double>();
		[ParameterName("lprwg"), ParameterInfo("Length dependence of prwg ")]
		public GivenParameter<double> Lprwg { get; } = new GivenParameter<double>();
		[ParameterName("lprt"), ParameterInfo("Length dependence of prt ")]
		public GivenParameter<double> Lprt { get; } = new GivenParameter<double>();
		[ParameterName("leta0"), ParameterInfo("Length dependence of eta0")]
		public GivenParameter<double> Leta0 { get; } = new GivenParameter<double>();
		[ParameterName("letab"), ParameterInfo("Length dependence of etab")]
		public GivenParameter<double> Letab { get; } = new GivenParameter<double>();
		[ParameterName("lpclm"), ParameterInfo("Length dependence of pclm")]
		public GivenParameter<double> Lpclm { get; } = new GivenParameter<double>();
		[ParameterName("lpdiblc1"), ParameterInfo("Length dependence of pdiblc1")]
		public GivenParameter<double> Lpdibl1 { get; } = new GivenParameter<double>();
		[ParameterName("lpdiblc2"), ParameterInfo("Length dependence of pdiblc2")]
		public GivenParameter<double> Lpdibl2 { get; } = new GivenParameter<double>();
		[ParameterName("lpdiblcb"), ParameterInfo("Length dependence of pdiblcb")]
		public GivenParameter<double> Lpdiblb { get; } = new GivenParameter<double>();
		[ParameterName("lpscbe1"), ParameterInfo("Length dependence of pscbe1")]
		public GivenParameter<double> Lpscbe1 { get; } = new GivenParameter<double>();
		[ParameterName("lpscbe2"), ParameterInfo("Length dependence of pscbe2")]
		public GivenParameter<double> Lpscbe2 { get; } = new GivenParameter<double>();
		[ParameterName("lpvag"), ParameterInfo("Length dependence of pvag")]
		public GivenParameter<double> Lpvag { get; } = new GivenParameter<double>();
		[ParameterName("lwr"), ParameterInfo("Length dependence of wr")]
		public GivenParameter<double> Lwr { get; } = new GivenParameter<double>();
		[ParameterName("ldwg"), ParameterInfo("Length dependence of dwg")]
		public GivenParameter<double> Ldwg { get; } = new GivenParameter<double>();
		[ParameterName("ldwb"), ParameterInfo("Length dependence of dwb")]
		public GivenParameter<double> Ldwb { get; } = new GivenParameter<double>();
		[ParameterName("lb0"), ParameterInfo("Length dependence of b0")]
		public GivenParameter<double> Lb0 { get; } = new GivenParameter<double>();
		[ParameterName("lb1"), ParameterInfo("Length dependence of b1")]
		public GivenParameter<double> Lb1 { get; } = new GivenParameter<double>();
		[ParameterName("lalpha0"), ParameterInfo("Length dependence of alpha0")]
		public GivenParameter<double> Lalpha0 { get; } = new GivenParameter<double>();
		[ParameterName("lalpha1"), ParameterInfo("Length dependence of alpha1")]
		public GivenParameter<double> Lalpha1 { get; } = new GivenParameter<double>();
		[ParameterName("lbeta0"), ParameterInfo("Length dependence of beta0")]
		public GivenParameter<double> Lbeta0 { get; } = new GivenParameter<double>();
		[ParameterName("lvfb"), ParameterInfo("Length dependence of vfb")]
		public GivenParameter<double> Lvfb { get; } = new GivenParameter<double>();
		[ParameterName("lelm"), ParameterInfo("Length dependence of elm")]
		public GivenParameter<double> Lelm { get; } = new GivenParameter<double>();
		[ParameterName("lcgsl"), ParameterInfo("Length dependence of cgsl")]
		public GivenParameter<double> Lcgsl { get; } = new GivenParameter<double>();
		[ParameterName("lcgdl"), ParameterInfo("Length dependence of cgdl")]
		public GivenParameter<double> Lcgdl { get; } = new GivenParameter<double>();
		[ParameterName("lckappa"), ParameterInfo("Length dependence of ckappa")]
		public GivenParameter<double> Lckappa { get; } = new GivenParameter<double>();
		[ParameterName("lcf"), ParameterInfo("Length dependence of cf")]
		public GivenParameter<double> Lcf { get; } = new GivenParameter<double>();
		[ParameterName("lclc"), ParameterInfo("Length dependence of clc")]
		public GivenParameter<double> Lclc { get; } = new GivenParameter<double>();
		[ParameterName("lcle"), ParameterInfo("Length dependence of cle")]
		public GivenParameter<double> Lcle { get; } = new GivenParameter<double>();
		[ParameterName("lvfbcv"), ParameterInfo("Length dependence of vfbcv")]
		public GivenParameter<double> Lvfbcv { get; } = new GivenParameter<double>();
		[ParameterName("lacde"), ParameterInfo("Length dependence of acde")]
		public GivenParameter<double> Lacde { get; } = new GivenParameter<double>();
		[ParameterName("lmoin"), ParameterInfo("Length dependence of moin")]
		public GivenParameter<double> Lmoin { get; } = new GivenParameter<double>();
		[ParameterName("lnoff"), ParameterInfo("Length dependence of noff")]
		public GivenParameter<double> Lnoff { get; } = new GivenParameter<double>();
		[ParameterName("lvoffcv"), ParameterInfo("Length dependence of voffcv")]
		public GivenParameter<double> Lvoffcv { get; } = new GivenParameter<double>();
		[ParameterName("wcdsc"), ParameterInfo("Width dependence of cdsc")]
		public GivenParameter<double> Wcdsc { get; } = new GivenParameter<double>();
		[ParameterName("wcdscb"), ParameterInfo("Width dependence of cdscb")]
		public GivenParameter<double> Wcdscb { get; } = new GivenParameter<double>();
		[ParameterName("wcdscd"), ParameterInfo("Width dependence of cdscd")]
		public GivenParameter<double> Wcdscd { get; } = new GivenParameter<double>();
		[ParameterName("wcit"), ParameterInfo("Width dependence of cit")]
		public GivenParameter<double> Wcit { get; } = new GivenParameter<double>();
		[ParameterName("wnfactor"), ParameterInfo("Width dependence of nfactor")]
		public GivenParameter<double> Wnfactor { get; } = new GivenParameter<double>();
		[ParameterName("wxj"), ParameterInfo("Width dependence of xj")]
		public GivenParameter<double> Wxj { get; } = new GivenParameter<double>();
		[ParameterName("wvsat"), ParameterInfo("Width dependence of vsat")]
		public GivenParameter<double> Wvsat { get; } = new GivenParameter<double>();
		[ParameterName("wa0"), ParameterInfo("Width dependence of a0")]
		public GivenParameter<double> Wa0 { get; } = new GivenParameter<double>();
		[ParameterName("wags"), ParameterInfo("Width dependence of ags")]
		public GivenParameter<double> Wags { get; } = new GivenParameter<double>();
		[ParameterName("wa1"), ParameterInfo("Width dependence of a1")]
		public GivenParameter<double> Wa1 { get; } = new GivenParameter<double>();
		[ParameterName("wa2"), ParameterInfo("Width dependence of a2")]
		public GivenParameter<double> Wa2 { get; } = new GivenParameter<double>();
		[ParameterName("wat"), ParameterInfo("Width dependence of at")]
		public GivenParameter<double> Wat { get; } = new GivenParameter<double>();
		[ParameterName("wketa"), ParameterInfo("Width dependence of keta")]
		public GivenParameter<double> Wketa { get; } = new GivenParameter<double>();
		[ParameterName("wnsub"), ParameterInfo("Width dependence of nsub")]
		public GivenParameter<double> Wnsub { get; } = new GivenParameter<double>();
		[ParameterName("wnch"), ParameterInfo("Width dependence of nch")]
		public GivenParameter<double> Wnpeak { get; } = new GivenParameter<double>();
		[ParameterName("wngate"), ParameterInfo("Width dependence of ngate")]
		public GivenParameter<double> Wngate { get; } = new GivenParameter<double>();
		[ParameterName("wgamma1"), ParameterInfo("Width dependence of gamma1")]
		public GivenParameter<double> Wgamma1 { get; } = new GivenParameter<double>();
		[ParameterName("wgamma2"), ParameterInfo("Width dependence of gamma2")]
		public GivenParameter<double> Wgamma2 { get; } = new GivenParameter<double>();
		[ParameterName("wvbx"), ParameterInfo("Width dependence of vbx")]
		public GivenParameter<double> Wvbx { get; } = new GivenParameter<double>();
		[ParameterName("wvbm"), ParameterInfo("Width dependence of vbm")]
		public GivenParameter<double> Wvbm { get; } = new GivenParameter<double>();
		[ParameterName("wxt"), ParameterInfo("Width dependence of xt")]
		public GivenParameter<double> Wxt { get; } = new GivenParameter<double>();
		[ParameterName("wk1"), ParameterInfo("Width dependence of k1")]
		public GivenParameter<double> Wk1 { get; } = new GivenParameter<double>();
		[ParameterName("wkt1"), ParameterInfo("Width dependence of kt1")]
		public GivenParameter<double> Wkt1 { get; } = new GivenParameter<double>();
		[ParameterName("wkt1l"), ParameterInfo("Width dependence of kt1l")]
		public GivenParameter<double> Wkt1l { get; } = new GivenParameter<double>();
		[ParameterName("wkt2"), ParameterInfo("Width dependence of kt2")]
		public GivenParameter<double> Wkt2 { get; } = new GivenParameter<double>();
		[ParameterName("wk2"), ParameterInfo("Width dependence of k2")]
		public GivenParameter<double> Wk2 { get; } = new GivenParameter<double>();
		[ParameterName("wk3"), ParameterInfo("Width dependence of k3")]
		public GivenParameter<double> Wk3 { get; } = new GivenParameter<double>();
		[ParameterName("wk3b"), ParameterInfo("Width dependence of k3b")]
		public GivenParameter<double> Wk3b { get; } = new GivenParameter<double>();
		[ParameterName("wnlx"), ParameterInfo("Width dependence of nlx")]
		public GivenParameter<double> Wnlx { get; } = new GivenParameter<double>();
		[ParameterName("ww0"), ParameterInfo("Width dependence of w0")]
		public GivenParameter<double> Ww0 { get; } = new GivenParameter<double>();
		[ParameterName("wdvt0"), ParameterInfo("Width dependence of dvt0")]
		public GivenParameter<double> Wdvt0 { get; } = new GivenParameter<double>();
		[ParameterName("wdvt1"), ParameterInfo("Width dependence of dvt1")]
		public GivenParameter<double> Wdvt1 { get; } = new GivenParameter<double>();
		[ParameterName("wdvt2"), ParameterInfo("Width dependence of dvt2")]
		public GivenParameter<double> Wdvt2 { get; } = new GivenParameter<double>();
		[ParameterName("wdvt0w"), ParameterInfo("Width dependence of dvt0w")]
		public GivenParameter<double> Wdvt0w { get; } = new GivenParameter<double>();
		[ParameterName("wdvt1w"), ParameterInfo("Width dependence of dvt1w")]
		public GivenParameter<double> Wdvt1w { get; } = new GivenParameter<double>();
		[ParameterName("wdvt2w"), ParameterInfo("Width dependence of dvt2w")]
		public GivenParameter<double> Wdvt2w { get; } = new GivenParameter<double>();
		[ParameterName("wdrout"), ParameterInfo("Width dependence of drout")]
		public GivenParameter<double> Wdrout { get; } = new GivenParameter<double>();
		[ParameterName("wdsub"), ParameterInfo("Width dependence of dsub")]
		public GivenParameter<double> Wdsub { get; } = new GivenParameter<double>();
		[ParameterName("wvth0"), ParameterName("wvtho"), ParameterInfo("Width dependence of vto")]
		public GivenParameter<double> Wvth0 { get; } = new GivenParameter<double>();
		[ParameterName("wua"), ParameterInfo("Width dependence of ua")]
		public GivenParameter<double> Wua { get; } = new GivenParameter<double>();
		[ParameterName("wua1"), ParameterInfo("Width dependence of ua1")]
		public GivenParameter<double> Wua1 { get; } = new GivenParameter<double>();
		[ParameterName("wub"), ParameterInfo("Width dependence of ub")]
		public GivenParameter<double> Wub { get; } = new GivenParameter<double>();
		[ParameterName("wub1"), ParameterInfo("Width dependence of ub1")]
		public GivenParameter<double> Wub1 { get; } = new GivenParameter<double>();
		[ParameterName("wuc"), ParameterInfo("Width dependence of uc")]
		public GivenParameter<double> Wuc { get; } = new GivenParameter<double>();
		[ParameterName("wuc1"), ParameterInfo("Width dependence of uc1")]
		public GivenParameter<double> Wuc1 { get; } = new GivenParameter<double>();
		[ParameterName("wu0"), ParameterInfo("Width dependence of u0")]
		public GivenParameter<double> Wu0 { get; } = new GivenParameter<double>();
		[ParameterName("wute"), ParameterInfo("Width dependence of ute")]
		public GivenParameter<double> Wute { get; } = new GivenParameter<double>();
		[ParameterName("wvoff"), ParameterInfo("Width dependence of voff")]
		public GivenParameter<double> Wvoff { get; } = new GivenParameter<double>();
		[ParameterName("wdelta"), ParameterInfo("Width dependence of delta")]
		public GivenParameter<double> Wdelta { get; } = new GivenParameter<double>();
		[ParameterName("wrdsw"), ParameterInfo("Width dependence of rdsw ")]
		public GivenParameter<double> Wrdsw { get; } = new GivenParameter<double>();
		[ParameterName("wprwb"), ParameterInfo("Width dependence of prwb ")]
		public GivenParameter<double> Wprwb { get; } = new GivenParameter<double>();
		[ParameterName("wprwg"), ParameterInfo("Width dependence of prwg ")]
		public GivenParameter<double> Wprwg { get; } = new GivenParameter<double>();
		[ParameterName("wprt"), ParameterInfo("Width dependence of prt")]
		public GivenParameter<double> Wprt { get; } = new GivenParameter<double>();
		[ParameterName("weta0"), ParameterInfo("Width dependence of eta0")]
		public GivenParameter<double> Weta0 { get; } = new GivenParameter<double>();
		[ParameterName("wetab"), ParameterInfo("Width dependence of etab")]
		public GivenParameter<double> Wetab { get; } = new GivenParameter<double>();
		[ParameterName("wpclm"), ParameterInfo("Width dependence of pclm")]
		public GivenParameter<double> Wpclm { get; } = new GivenParameter<double>();
		[ParameterName("wpdiblc1"), ParameterInfo("Width dependence of pdiblc1")]
		public GivenParameter<double> Wpdibl1 { get; } = new GivenParameter<double>();
		[ParameterName("wpdiblc2"), ParameterInfo("Width dependence of pdiblc2")]
		public GivenParameter<double> Wpdibl2 { get; } = new GivenParameter<double>();
		[ParameterName("wpdiblcb"), ParameterInfo("Width dependence of pdiblcb")]
		public GivenParameter<double> Wpdiblb { get; } = new GivenParameter<double>();
		[ParameterName("wpscbe1"), ParameterInfo("Width dependence of pscbe1")]
		public GivenParameter<double> Wpscbe1 { get; } = new GivenParameter<double>();
		[ParameterName("wpscbe2"), ParameterInfo("Width dependence of pscbe2")]
		public GivenParameter<double> Wpscbe2 { get; } = new GivenParameter<double>();
		[ParameterName("wpvag"), ParameterInfo("Width dependence of pvag")]
		public GivenParameter<double> Wpvag { get; } = new GivenParameter<double>();
		[ParameterName("wwr"), ParameterInfo("Width dependence of wr")]
		public GivenParameter<double> Wwr { get; } = new GivenParameter<double>();
		[ParameterName("wdwg"), ParameterInfo("Width dependence of dwg")]
		public GivenParameter<double> Wdwg { get; } = new GivenParameter<double>();
		[ParameterName("wdwb"), ParameterInfo("Width dependence of dwb")]
		public GivenParameter<double> Wdwb { get; } = new GivenParameter<double>();
		[ParameterName("wb0"), ParameterInfo("Width dependence of b0")]
		public GivenParameter<double> Wb0 { get; } = new GivenParameter<double>();
		[ParameterName("wb1"), ParameterInfo("Width dependence of b1")]
		public GivenParameter<double> Wb1 { get; } = new GivenParameter<double>();
		[ParameterName("walpha0"), ParameterInfo("Width dependence of alpha0")]
		public GivenParameter<double> Walpha0 { get; } = new GivenParameter<double>();
		[ParameterName("walpha1"), ParameterInfo("Width dependence of alpha1")]
		public GivenParameter<double> Walpha1 { get; } = new GivenParameter<double>();
		[ParameterName("wbeta0"), ParameterInfo("Width dependence of beta0")]
		public GivenParameter<double> Wbeta0 { get; } = new GivenParameter<double>();
		[ParameterName("wvfb"), ParameterInfo("Width dependence of vfb")]
		public GivenParameter<double> Wvfb { get; } = new GivenParameter<double>();
		[ParameterName("welm"), ParameterInfo("Width dependence of elm")]
		public GivenParameter<double> Welm { get; } = new GivenParameter<double>();
		[ParameterName("wcgsl"), ParameterInfo("Width dependence of cgsl")]
		public GivenParameter<double> Wcgsl { get; } = new GivenParameter<double>();
		[ParameterName("wcgdl"), ParameterInfo("Width dependence of cgdl")]
		public GivenParameter<double> Wcgdl { get; } = new GivenParameter<double>();
		[ParameterName("wckappa"), ParameterInfo("Width dependence of ckappa")]
		public GivenParameter<double> Wckappa { get; } = new GivenParameter<double>();
		[ParameterName("wcf"), ParameterInfo("Width dependence of cf")]
		public GivenParameter<double> Wcf { get; } = new GivenParameter<double>();
		[ParameterName("wclc"), ParameterInfo("Width dependence of clc")]
		public GivenParameter<double> Wclc { get; } = new GivenParameter<double>();
		[ParameterName("wcle"), ParameterInfo("Width dependence of cle")]
		public GivenParameter<double> Wcle { get; } = new GivenParameter<double>();
		[ParameterName("wvfbcv"), ParameterInfo("Width dependence of vfbcv")]
		public GivenParameter<double> Wvfbcv { get; } = new GivenParameter<double>();
		[ParameterName("wacde"), ParameterInfo("Width dependence of acde")]
		public GivenParameter<double> Wacde { get; } = new GivenParameter<double>();
		[ParameterName("wmoin"), ParameterInfo("Width dependence of moin")]
		public GivenParameter<double> Wmoin { get; } = new GivenParameter<double>();
		[ParameterName("wnoff"), ParameterInfo("Width dependence of noff")]
		public GivenParameter<double> Wnoff { get; } = new GivenParameter<double>();
		[ParameterName("wvoffcv"), ParameterInfo("Width dependence of voffcv")]
		public GivenParameter<double> Wvoffcv { get; } = new GivenParameter<double>();
		[ParameterName("pcdsc"), ParameterInfo("Cross-term dependence of cdsc")]
		public GivenParameter<double> Pcdsc { get; } = new GivenParameter<double>();
		[ParameterName("pcdscb"), ParameterInfo("Cross-term dependence of cdscb")]
		public GivenParameter<double> Pcdscb { get; } = new GivenParameter<double>();
		[ParameterName("pcdscd"), ParameterInfo("Cross-term dependence of cdscd")]
		public GivenParameter<double> Pcdscd { get; } = new GivenParameter<double>();
		[ParameterName("pcit"), ParameterInfo("Cross-term dependence of cit")]
		public GivenParameter<double> Pcit { get; } = new GivenParameter<double>();
		[ParameterName("pnfactor"), ParameterInfo("Cross-term dependence of nfactor")]
		public GivenParameter<double> Pnfactor { get; } = new GivenParameter<double>();
		[ParameterName("pxj"), ParameterInfo("Cross-term dependence of xj")]
		public GivenParameter<double> Pxj { get; } = new GivenParameter<double>();
		[ParameterName("pvsat"), ParameterInfo("Cross-term dependence of vsat")]
		public GivenParameter<double> Pvsat { get; } = new GivenParameter<double>();
		[ParameterName("pa0"), ParameterInfo("Cross-term dependence of a0")]
		public GivenParameter<double> Pa0 { get; } = new GivenParameter<double>();
		[ParameterName("pags"), ParameterInfo("Cross-term dependence of ags")]
		public GivenParameter<double> Pags { get; } = new GivenParameter<double>();
		[ParameterName("pa1"), ParameterInfo("Cross-term dependence of a1")]
		public GivenParameter<double> Pa1 { get; } = new GivenParameter<double>();
		[ParameterName("pa2"), ParameterInfo("Cross-term dependence of a2")]
		public GivenParameter<double> Pa2 { get; } = new GivenParameter<double>();
		[ParameterName("pat"), ParameterInfo("Cross-term dependence of at")]
		public GivenParameter<double> Pat { get; } = new GivenParameter<double>();
		[ParameterName("pketa"), ParameterInfo("Cross-term dependence of keta")]
		public GivenParameter<double> Pketa { get; } = new GivenParameter<double>();
		[ParameterName("pnsub"), ParameterInfo("Cross-term dependence of nsub")]
		public GivenParameter<double> Pnsub { get; } = new GivenParameter<double>();
		[ParameterName("pnch"), ParameterInfo("Cross-term dependence of nch")]
		public GivenParameter<double> Pnpeak { get; } = new GivenParameter<double>();
		[ParameterName("pngate"), ParameterInfo("Cross-term dependence of ngate")]
		public GivenParameter<double> Pngate { get; } = new GivenParameter<double>();
		[ParameterName("pgamma1"), ParameterInfo("Cross-term dependence of gamma1")]
		public GivenParameter<double> Pgamma1 { get; } = new GivenParameter<double>();
		[ParameterName("pgamma2"), ParameterInfo("Cross-term dependence of gamma2")]
		public GivenParameter<double> Pgamma2 { get; } = new GivenParameter<double>();
		[ParameterName("pvbx"), ParameterInfo("Cross-term dependence of vbx")]
		public GivenParameter<double> Pvbx { get; } = new GivenParameter<double>();
		[ParameterName("pvbm"), ParameterInfo("Cross-term dependence of vbm")]
		public GivenParameter<double> Pvbm { get; } = new GivenParameter<double>();
		[ParameterName("pxt"), ParameterInfo("Cross-term dependence of xt")]
		public GivenParameter<double> Pxt { get; } = new GivenParameter<double>();
		[ParameterName("pk1"), ParameterInfo("Cross-term dependence of k1")]
		public GivenParameter<double> Pk1 { get; } = new GivenParameter<double>();
		[ParameterName("pkt1"), ParameterInfo("Cross-term dependence of kt1")]
		public GivenParameter<double> Pkt1 { get; } = new GivenParameter<double>();
		[ParameterName("pkt1l"), ParameterInfo("Cross-term dependence of kt1l")]
		public GivenParameter<double> Pkt1l { get; } = new GivenParameter<double>();
		[ParameterName("pkt2"), ParameterInfo("Cross-term dependence of kt2")]
		public GivenParameter<double> Pkt2 { get; } = new GivenParameter<double>();
		[ParameterName("pk2"), ParameterInfo("Cross-term dependence of k2")]
		public GivenParameter<double> Pk2 { get; } = new GivenParameter<double>();
		[ParameterName("pk3"), ParameterInfo("Cross-term dependence of k3")]
		public GivenParameter<double> Pk3 { get; } = new GivenParameter<double>();
		[ParameterName("pk3b"), ParameterInfo("Cross-term dependence of k3b")]
		public GivenParameter<double> Pk3b { get; } = new GivenParameter<double>();
		[ParameterName("pnlx"), ParameterInfo("Cross-term dependence of nlx")]
		public GivenParameter<double> Pnlx { get; } = new GivenParameter<double>();
		[ParameterName("pw0"), ParameterInfo("Cross-term dependence of w0")]
		public GivenParameter<double> Pw0 { get; } = new GivenParameter<double>();
		[ParameterName("pdvt0"), ParameterInfo("Cross-term dependence of dvt0")]
		public GivenParameter<double> Pdvt0 { get; } = new GivenParameter<double>();
		[ParameterName("pdvt1"), ParameterInfo("Cross-term dependence of dvt1")]
		public GivenParameter<double> Pdvt1 { get; } = new GivenParameter<double>();
		[ParameterName("pdvt2"), ParameterInfo("Cross-term dependence of dvt2")]
		public GivenParameter<double> Pdvt2 { get; } = new GivenParameter<double>();
		[ParameterName("pdvt0w"), ParameterInfo("Cross-term dependence of dvt0w")]
		public GivenParameter<double> Pdvt0w { get; } = new GivenParameter<double>();
		[ParameterName("pdvt1w"), ParameterInfo("Cross-term dependence of dvt1w")]
		public GivenParameter<double> Pdvt1w { get; } = new GivenParameter<double>();
		[ParameterName("pdvt2w"), ParameterInfo("Cross-term dependence of dvt2w")]
		public GivenParameter<double> Pdvt2w { get; } = new GivenParameter<double>();
		[ParameterName("pdrout"), ParameterInfo("Cross-term dependence of drout")]
		public GivenParameter<double> Pdrout { get; } = new GivenParameter<double>();
		[ParameterName("pdsub"), ParameterInfo("Cross-term dependence of dsub")]
		public GivenParameter<double> Pdsub { get; } = new GivenParameter<double>();
		[ParameterName("pvth0"), ParameterName("pvtho"), ParameterInfo("Cross-term dependence of vto")]
		public GivenParameter<double> Pvth0 { get; } = new GivenParameter<double>();
		[ParameterName("pua"), ParameterInfo("Cross-term dependence of ua")]
		public GivenParameter<double> Pua { get; } = new GivenParameter<double>();
		[ParameterName("pua1"), ParameterInfo("Cross-term dependence of ua1")]
		public GivenParameter<double> Pua1 { get; } = new GivenParameter<double>();
		[ParameterName("pub"), ParameterInfo("Cross-term dependence of ub")]
		public GivenParameter<double> Pub { get; } = new GivenParameter<double>();
		[ParameterName("pub1"), ParameterInfo("Cross-term dependence of ub1")]
		public GivenParameter<double> Pub1 { get; } = new GivenParameter<double>();
		[ParameterName("puc"), ParameterInfo("Cross-term dependence of uc")]
		public GivenParameter<double> Puc { get; } = new GivenParameter<double>();
		[ParameterName("puc1"), ParameterInfo("Cross-term dependence of uc1")]
		public GivenParameter<double> Puc1 { get; } = new GivenParameter<double>();
		[ParameterName("pu0"), ParameterInfo("Cross-term dependence of u0")]
		public GivenParameter<double> Pu0 { get; } = new GivenParameter<double>();
		[ParameterName("pute"), ParameterInfo("Cross-term dependence of ute")]
		public GivenParameter<double> Pute { get; } = new GivenParameter<double>();
		[ParameterName("pvoff"), ParameterInfo("Cross-term dependence of voff")]
		public GivenParameter<double> Pvoff { get; } = new GivenParameter<double>();
		[ParameterName("pdelta"), ParameterInfo("Cross-term dependence of delta")]
		public GivenParameter<double> Pdelta { get; } = new GivenParameter<double>();
		[ParameterName("prdsw"), ParameterInfo("Cross-term dependence of rdsw ")]
		public GivenParameter<double> Prdsw { get; } = new GivenParameter<double>();
		[ParameterName("pprwb"), ParameterInfo("Cross-term dependence of prwb ")]
		public GivenParameter<double> Pprwb { get; } = new GivenParameter<double>();
		[ParameterName("pprwg"), ParameterInfo("Cross-term dependence of prwg ")]
		public GivenParameter<double> Pprwg { get; } = new GivenParameter<double>();
		[ParameterName("pprt"), ParameterInfo("Cross-term dependence of prt ")]
		public GivenParameter<double> Pprt { get; } = new GivenParameter<double>();
		[ParameterName("peta0"), ParameterInfo("Cross-term dependence of eta0")]
		public GivenParameter<double> Peta0 { get; } = new GivenParameter<double>();
		[ParameterName("petab"), ParameterInfo("Cross-term dependence of etab")]
		public GivenParameter<double> Petab { get; } = new GivenParameter<double>();
		[ParameterName("ppclm"), ParameterInfo("Cross-term dependence of pclm")]
		public GivenParameter<double> Ppclm { get; } = new GivenParameter<double>();
		[ParameterName("ppdiblc1"), ParameterInfo("Cross-term dependence of pdiblc1")]
		public GivenParameter<double> Ppdibl1 { get; } = new GivenParameter<double>();
		[ParameterName("ppdiblc2"), ParameterInfo("Cross-term dependence of pdiblc2")]
		public GivenParameter<double> Ppdibl2 { get; } = new GivenParameter<double>();
		[ParameterName("ppdiblcb"), ParameterInfo("Cross-term dependence of pdiblcb")]
		public GivenParameter<double> Ppdiblb { get; } = new GivenParameter<double>();
		[ParameterName("ppscbe1"), ParameterInfo("Cross-term dependence of pscbe1")]
		public GivenParameter<double> Ppscbe1 { get; } = new GivenParameter<double>();
		[ParameterName("ppscbe2"), ParameterInfo("Cross-term dependence of pscbe2")]
		public GivenParameter<double> Ppscbe2 { get; } = new GivenParameter<double>();
		[ParameterName("ppvag"), ParameterInfo("Cross-term dependence of pvag")]
		public GivenParameter<double> Ppvag { get; } = new GivenParameter<double>();
		[ParameterName("pwr"), ParameterInfo("Cross-term dependence of wr")]
		public GivenParameter<double> Pwr { get; } = new GivenParameter<double>();
		[ParameterName("pdwg"), ParameterInfo("Cross-term dependence of dwg")]
		public GivenParameter<double> Pdwg { get; } = new GivenParameter<double>();
		[ParameterName("pdwb"), ParameterInfo("Cross-term dependence of dwb")]
		public GivenParameter<double> Pdwb { get; } = new GivenParameter<double>();
		[ParameterName("pb0"), ParameterInfo("Cross-term dependence of b0")]
		public GivenParameter<double> Pb0 { get; } = new GivenParameter<double>();
		[ParameterName("pb1"), ParameterInfo("Cross-term dependence of b1")]
		public GivenParameter<double> Pb1 { get; } = new GivenParameter<double>();
		[ParameterName("palpha0"), ParameterInfo("Cross-term dependence of alpha0")]
		public GivenParameter<double> Palpha0 { get; } = new GivenParameter<double>();
		[ParameterName("palpha1"), ParameterInfo("Cross-term dependence of alpha1")]
		public GivenParameter<double> Palpha1 { get; } = new GivenParameter<double>();
		[ParameterName("pbeta0"), ParameterInfo("Cross-term dependence of beta0")]
		public GivenParameter<double> Pbeta0 { get; } = new GivenParameter<double>();
		[ParameterName("pvfb"), ParameterInfo("Cross-term dependence of vfb")]
		public GivenParameter<double> Pvfb { get; } = new GivenParameter<double>();
		[ParameterName("pelm"), ParameterInfo("Cross-term dependence of elm")]
		public GivenParameter<double> Pelm { get; } = new GivenParameter<double>();
		[ParameterName("pcgsl"), ParameterInfo("Cross-term dependence of cgsl")]
		public GivenParameter<double> Pcgsl { get; } = new GivenParameter<double>();
		[ParameterName("pcgdl"), ParameterInfo("Cross-term dependence of cgdl")]
		public GivenParameter<double> Pcgdl { get; } = new GivenParameter<double>();
		[ParameterName("pckappa"), ParameterInfo("Cross-term dependence of ckappa")]
		public GivenParameter<double> Pckappa { get; } = new GivenParameter<double>();
		[ParameterName("pcf"), ParameterInfo("Cross-term dependence of cf")]
		public GivenParameter<double> Pcf { get; } = new GivenParameter<double>();
		[ParameterName("pclc"), ParameterInfo("Cross-term dependence of clc")]
		public GivenParameter<double> Pclc { get; } = new GivenParameter<double>();
		[ParameterName("pcle"), ParameterInfo("Cross-term dependence of cle")]
		public GivenParameter<double> Pcle { get; } = new GivenParameter<double>();
		[ParameterName("pvfbcv"), ParameterInfo("Cross-term dependence of vfbcv")]
		public GivenParameter<double> Pvfbcv { get; } = new GivenParameter<double>();
		[ParameterName("pacde"), ParameterInfo("Cross-term dependence of acde")]
		public GivenParameter<double> Pacde { get; } = new GivenParameter<double>();
		[ParameterName("pmoin"), ParameterInfo("Cross-term dependence of moin")]
		public GivenParameter<double> Pmoin { get; } = new GivenParameter<double>();
		[ParameterName("pnoff"), ParameterInfo("Cross-term dependence of noff")]
		public GivenParameter<double> Pnoff { get; } = new GivenParameter<double>();
		[ParameterName("pvoffcv"), ParameterInfo("Cross-term dependence of voffcv")]
		public GivenParameter<double> Pvoffcv { get; } = new GivenParameter<double>();
		[ParameterName("cgso"), ParameterInfo("Gate-source overlap capacitance per width")]
		public GivenParameter<double> Cgso { get; } = new GivenParameter<double>();
		[ParameterName("cgdo"), ParameterInfo("Gate-drain overlap capacitance per width")]
		public GivenParameter<double> Cgdo { get; } = new GivenParameter<double>();
		[ParameterName("cgbo"), ParameterInfo("Gate-bulk overlap capacitance per length")]
		public GivenParameter<double> Cgbo { get; } = new GivenParameter<double>();
		[ParameterName("xpart"), ParameterInfo("Channel charge partitioning")]
		public GivenParameter<double> Xpart { get; } = new GivenParameter<double>();
		[ParameterName("rsh"), ParameterInfo("Source-drain sheet resistance")]
		public GivenParameter<double> SheetResistance { get; } = new GivenParameter<double>();
		[ParameterName("js"), ParameterInfo("Source/drain junction reverse saturation current density")]
		public GivenParameter<double> JctSatCurDensity { get; } = new GivenParameter<double>(0.0001);
		[ParameterName("jsw"), ParameterInfo("Sidewall junction reverse saturation current density")]
		public GivenParameter<double> JctSidewallSatCurDensity { get; } = new GivenParameter<double>();
		[ParameterName("pb"), ParameterInfo("Source/drain junction built-in potential")]
		public GivenParameter<double> BulkJctPotential { get; } = new GivenParameter<double>(1);
		[ParameterName("mj"), ParameterInfo("Source/drain bottom junction capacitance grading coefficient")]
		public GivenParameter<double> BulkJctBotGradingCoeff { get; } = new GivenParameter<double>(0.5);
		[ParameterName("pbsw"), ParameterInfo("Source/drain sidewall junction capacitance built in potential")]
		public GivenParameter<double> SidewallJctPotential { get; } = new GivenParameter<double>(1);
		[ParameterName("mjsw"), ParameterInfo("Source/drain sidewall junction capacitance grading coefficient")]
		public GivenParameter<double> BulkJctSideGradingCoeff { get; } = new GivenParameter<double>(0.33);
		[ParameterName("cj"), ParameterInfo("Source/drain bottom junction capacitance per unit area")]
		public GivenParameter<double> UnitAreaJctCap { get; } = new GivenParameter<double>(0.0005);
		[ParameterName("cjsw"), ParameterInfo("Source/drain sidewall junction capacitance per unit periphery")]
		public GivenParameter<double> UnitLengthSidewallJctCap { get; } = new GivenParameter<double>(5e-10);
		[ParameterName("nj"), ParameterInfo("Source/drain junction emission coefficient")]
		public GivenParameter<double> JctEmissionCoeff { get; } = new GivenParameter<double>(1);
		[ParameterName("pbswg"), ParameterInfo("Source/drain (gate side) sidewall junction capacitance built in potential")]
		public GivenParameter<double> GatesidewallJctPotential { get; } = new GivenParameter<double>();
		[ParameterName("mjswg"), ParameterInfo("Source/drain (gate side) sidewall junction capacitance grading coefficient")]
		public GivenParameter<double> BulkJctGateSideGradingCoeff { get; } = new GivenParameter<double>();
		[ParameterName("cjswg"), ParameterInfo("Source/drain (gate side) sidewall junction capacitance per unit width")]
		public GivenParameter<double> UnitLengthGateSidewallJctCap { get; } = new GivenParameter<double>();
		[ParameterName("xti"), ParameterInfo("Junction current temperature exponent")]
		public GivenParameter<double> JctTempExponent { get; } = new GivenParameter<double>(3);
		[ParameterName("lint"), ParameterInfo("Length reduction parameter")]
		public GivenParameter<double> Lint { get; } = new GivenParameter<double>();
		[ParameterName("ll"), ParameterInfo("Length reduction parameter")]
		public GivenParameter<double> Ll { get; } = new GivenParameter<double>();
		[ParameterName("llc"), ParameterInfo("Length reduction parameter for CV")]
		public GivenParameter<double> Llc { get; } = new GivenParameter<double>();
		[ParameterName("lln"), ParameterInfo("Length reduction parameter")]
		public GivenParameter<double> Lln { get; } = new GivenParameter<double>(1);
		[ParameterName("lw"), ParameterInfo("Length reduction parameter")]
		public GivenParameter<double> Lw { get; } = new GivenParameter<double>();
		[ParameterName("lwc"), ParameterInfo("Length reduction parameter for CV")]
		public GivenParameter<double> Lwc { get; } = new GivenParameter<double>();
		[ParameterName("lwn"), ParameterInfo("Length reduction parameter")]
		public GivenParameter<double> Lwn { get; } = new GivenParameter<double>(1);
		[ParameterName("lwl"), ParameterInfo("Length reduction parameter")]
		public GivenParameter<double> Lwl { get; } = new GivenParameter<double>();
		[ParameterName("lwlc"), ParameterInfo("Length reduction parameter for CV")]
		public GivenParameter<double> Lwlc { get; } = new GivenParameter<double>();
		[ParameterName("lmin"), ParameterInfo("Minimum length for the model")]
		public GivenParameter<double> Lmin { get; } = new GivenParameter<double>();
		[ParameterName("lmax"), ParameterInfo("Maximum length for the model")]
		public GivenParameter<double> Lmax { get; } = new GivenParameter<double>(1);
		[ParameterName("wint"), ParameterInfo("Width reduction parameter")]
		public GivenParameter<double> Wint { get; } = new GivenParameter<double>();
		[ParameterName("wl"), ParameterInfo("Width reduction parameter")]
		public GivenParameter<double> Wl { get; } = new GivenParameter<double>();
		[ParameterName("wlc"), ParameterInfo("Width reduction parameter for CV")]
		public GivenParameter<double> Wlc { get; } = new GivenParameter<double>();
		[ParameterName("wln"), ParameterInfo("Width reduction parameter")]
		public GivenParameter<double> Wln { get; } = new GivenParameter<double>(1);
		[ParameterName("ww"), ParameterInfo("Width reduction parameter")]
		public GivenParameter<double> Ww { get; } = new GivenParameter<double>();
		[ParameterName("wwc"), ParameterInfo("Width reduction parameter for CV")]
		public GivenParameter<double> Wwc { get; } = new GivenParameter<double>();
		[ParameterName("wwn"), ParameterInfo("Width reduction parameter")]
		public GivenParameter<double> Wwn { get; } = new GivenParameter<double>(1);
		[ParameterName("wwl"), ParameterInfo("Width reduction parameter")]
		public GivenParameter<double> Wwl { get; } = new GivenParameter<double>();
		[ParameterName("wwlc"), ParameterInfo("Width reduction parameter for CV")]
		public GivenParameter<double> Wwlc { get; } = new GivenParameter<double>();
		[ParameterName("wmin"), ParameterInfo("Minimum width for the model")]
		public GivenParameter<double> Wmin { get; } = new GivenParameter<double>();
		[ParameterName("wmax"), ParameterInfo("Maximum width for the model")]
		public GivenParameter<double> Wmax { get; } = new GivenParameter<double>(1);
		[ParameterName("noia"), ParameterInfo("Flicker noise parameter")]
		public GivenParameter<double> OxideTrapDensityA { get; } = new GivenParameter<double>();
		[ParameterName("noib"), ParameterInfo("Flicker noise parameter")]
		public GivenParameter<double> OxideTrapDensityB { get; } = new GivenParameter<double>();
		[ParameterName("noic"), ParameterInfo("Flicker noise parameter")]
		public GivenParameter<double> OxideTrapDensityC { get; } = new GivenParameter<double>();
		[ParameterName("em"), ParameterInfo("Flicker noise parameter")]
		public GivenParameter<double> Em { get; } = new GivenParameter<double>(41000000);
		[ParameterName("ef"), ParameterInfo("Flicker noise frequency exponent")]
		public GivenParameter<double> Ef { get; } = new GivenParameter<double>(1);
		[ParameterName("af"), ParameterInfo("Flicker noise exponent")]
		public GivenParameter<double> Af { get; } = new GivenParameter<double>(1);
		[ParameterName("kf"), ParameterInfo("Flicker noise coefficient")]
		public GivenParameter<double> Kf { get; } = new GivenParameter<double>();

	    public double B3Type { get; private set; } = 1.0;
		public double Cox { get; private set; }

	    [ParameterName("tnom"), ParameterInfo("Parameter measurement temperature")]
	    public double TnomCelsius
	    {
	        get => Tnom.Value - Constants.CelsiusKelvin;
	        set => Tnom.Value = value + Constants.CelsiusKelvin;
	    }
	    public GivenParameter<double> Tnom { get; } = new GivenParameter<double>(300.15);

	    [ParameterName("nmos"), ParameterInfo("Flag to indicate NMOS")]
	    public void SetNmos(bool flag = true)
	    {
	        if (flag)
	            B3Type = 1.0;
	    }
	    [ParameterName("pmos"), ParameterInfo("Flag to indicate PMOS")]
        public void SetPmos(bool flag = true)
	    {
	        if (flag)
	            B3Type = -1.0;
	    }

        /// <summary>
        /// Clone the parameter set.
        /// </summary>
        /// <returns></returns>
	    public override ParameterSet Clone()
	    {
	        var clone = (ModelBaseParameters) base.Clone();
	        clone.B3Type = B3Type;
	        clone.Cox = Cox;
	        return clone;
	    }

        /// <summary>
        /// Copy from another parameter set.
        /// </summary>
        /// <param name="source">The source.</param>
        public override void CopyFrom(ParameterSet source)
        {
            if (source is ModelBaseParameters mbp)
            {
                B3Type = mbp.B3Type;
                Cox = mbp.Cox;
                base.CopyFrom(mbp);
            }
        }

        /// <summary>
        /// Calculate default parameters
        /// </summary>
        public override void CalculateDefaults()
		{
			// Parameter set for Npeak (BSIM3npeak)
			if (Npeak > 1.0e20)
				Npeak.RawValue *= 1.0e-6;
			// Parameter set for Ngate.RawValue (BSIM3ngate.RawValue)
			if (Ngate.RawValue > 1.0e23)
				Ngate.RawValue *= 1.0e-6;
			// Parameter set for Lnpeak.RawValue (BSIM3lnpeak.RawValue)
			if (Lnpeak.RawValue > 1.0e20)
				Lnpeak.RawValue *= 1.0e-6;
			// Parameter set for Lngate.RawValue (BSIM3lngate.RawValue)
			if (Lngate.RawValue > 1.0e23)
				Lngate.RawValue *= 1.0e-6;
			// Parameter set for Wnpeak.RawValue (BSIM3wnpeak.RawValue)
			if (Wnpeak.RawValue > 1.0e20)
				Wnpeak.RawValue *= 1.0e-6;
			// Parameter set for Wngate.RawValue (BSIM3wngate.RawValue)
			if (Wngate.RawValue > 1.0e23)
				Wngate.RawValue *= 1.0e-6;
			// Parameter set for Pnpeak.RawValue (BSIM3pnpeak.RawValue)
			if (Pnpeak.RawValue > 1.0e20)
				Pnpeak.RawValue *= 1.0e-6;
			// Parameter set for Pngate.RawValue (BSIM3pngate.RawValue)
			if (Pngate.RawValue > 1.0e23)
				Pngate.RawValue *= 1.0e-6;
			Cox = 3.453133e-11 / Tox;
			if (!Toxm.Given)
				Toxm.RawValue = Tox;
			if (!Dsub.Given)
				Dsub.RawValue = Drout;
			if (!Vth0.Given)
				Vth0.RawValue = B3Type > 1 ? 0.7 : -0.7;
			if (!Uc.Given)
				Uc.RawValue = MobMod.RawValue == 3 ? -0.0465 : -0.0465e-9;
			if (!Uc1.Given)
				Uc1.RawValue = MobMod.RawValue == 3 ? -0.056 : -0.056e-9;
			if (!U0.Given)
				U0.RawValue = B3Type > 1 ? 0.067 : 0.025;
			if (!Llc.Given)
				Llc.RawValue = Ll;
			if (!Lwc.Given)
				Lwc.RawValue = Lw;
			if (!Lwlc.Given)
				Lwlc.RawValue = Lwl;
			if (!Wlc.Given)
				Wlc.RawValue = Wl;
			if (!Wwc.Given)
				Wwc.RawValue = Ww;
			if (!Wwlc.Given)
				Wwlc.RawValue = Wwl;
			if (!Dwc.Given)
				Dwc.RawValue = Wint;
			if (!Dlc.Given)
				Dlc.RawValue = Lint;
			if (!Cf.Given)
				Cf.RawValue = 2.0 * 3.453133e-11 / 3.141592654 * Math.Log(1.0 + 0.4e-6 / Tox);
			if (!Cgdo.Given)
			{
				if (Dlc.Given && Dlc > 0.0)
				{
					Cgdo.RawValue = Dlc * Cox - Cgdl;
				}
				else
					Cgdo.RawValue = 0.6 * Xj * Cox;
			}
			if (!Cgso.Given)
			{
				if (Dlc.Given && Dlc > 0.0)
				{
					Cgso.RawValue = Dlc * Cox - Cgsl;
				}
				else
					Cgso.RawValue = 0.6 * Xj * Cox;
			}
			if (!Cgbo.Given)
				Cgbo.RawValue = 2.0 * Dwc * Cox;
			if (!UnitLengthGateSidewallJctCap.Given)
				UnitLengthGateSidewallJctCap.RawValue = UnitLengthSidewallJctCap;
			if (!GatesidewallJctPotential.Given)
				GatesidewallJctPotential.RawValue = SidewallJctPotential;
			if (!BulkJctGateSideGradingCoeff.Given)
				BulkJctGateSideGradingCoeff.RawValue = BulkJctSideGradingCoeff;
			if (!OxideTrapDensityA.Given)
			{
				if (B3Type > 0)
					OxideTrapDensityA.RawValue = 1e20;
				else
					OxideTrapDensityA.RawValue = 9.9e18;
			}
			if (!OxideTrapDensityB.Given)
			{
				if (B3Type > 0)
					OxideTrapDensityB.RawValue = 5e4;
				else
					OxideTrapDensityB.RawValue = 2.4e3;
			}
			if (!OxideTrapDensityC.Given)
			{
				if (B3Type > 0)
					OxideTrapDensityC.RawValue = -1.4e-12;
				else
					OxideTrapDensityC.RawValue = 1.4e-12;
			}

		    if (BulkJctPotential < 0.1)
		    {
		        BulkJctPotential.RawValue = 0.1;
		        CircuitWarning.Warning(this, "Given pb is less than 0.1. Pb is set to 0.1.");
		    }

		    if (SidewallJctPotential < 0.1)
		    {
		        SidewallJctPotential.RawValue = 0.1;
		        CircuitWarning.Warning(this, "Given pbsw is less than 0.1. Pbsw is set to 0.1.");
		    }

		    if (GatesidewallJctPotential < 0.1)
		    {
		        GatesidewallJctPotential.RawValue = 0.1;
		        CircuitWarning.Warning(this, "Given pbswg is less than 0.1. Pbswg is set to 0.1.");
		    }
		}
	}
}