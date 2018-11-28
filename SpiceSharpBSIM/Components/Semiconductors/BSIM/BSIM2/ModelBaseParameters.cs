using System;
using SpiceSharp.Attributes;

namespace SpiceSharp.Components.BSIM2Behaviors
{
	/// <summary>
	/// Base parameters for a <see cref="BSIM2Model" />
	/// </summary>
	public class ModelBaseParameters : ParameterSet
	{
		/// <summary>
		/// Properties
		/// </summary>
		public double Type { get; private set; } = 1;
		[ParameterName("vfb"), ParameterInfo("Flat band voltage")]
		public GivenParameter<double> Vfb0 { get; } = new GivenParameter<double>(-1);
		[ParameterName("lvfb"), ParameterInfo("Length dependence of vfb")]
		public GivenParameter<double> VfbL { get; } = new GivenParameter<double>();
		[ParameterName("wvfb"), ParameterInfo("Width dependence of vfb")]
		public GivenParameter<double> VfbW { get; } = new GivenParameter<double>();
		[ParameterName("phi"), ParameterInfo("Strong inversion surface potential ")]
		public GivenParameter<double> Phi0 { get; } = new GivenParameter<double>(0.75);
		[ParameterName("lphi"), ParameterInfo("Length dependence of phi")]
		public GivenParameter<double> PhiL { get; } = new GivenParameter<double>();
		[ParameterName("wphi"), ParameterInfo("Width dependence of phi")]
		public GivenParameter<double> PhiW { get; } = new GivenParameter<double>();
		[ParameterName("k1"), ParameterInfo("Bulk effect coefficient 1")]
		public GivenParameter<double> K10 { get; } = new GivenParameter<double>(0.8);
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
		[ParameterName("eta0"), ParameterInfo("VDS dependence of threshold voltage at VDD=0")]
		public GivenParameter<double> Eta00 { get; } = new GivenParameter<double>();
		[ParameterName("leta0"), ParameterInfo("Length dependence of eta0")]
		public GivenParameter<double> Eta0L { get; } = new GivenParameter<double>();
		[ParameterName("weta0"), ParameterInfo("Width dependence of eta0")]
		public GivenParameter<double> Eta0W { get; } = new GivenParameter<double>();
		[ParameterName("etab"), ParameterInfo("VBS dependence of eta")]
		public GivenParameter<double> EtaB0 { get; } = new GivenParameter<double>();
		[ParameterName("letab"), ParameterInfo("Length dependence of etab")]
		public GivenParameter<double> EtaBL { get; } = new GivenParameter<double>();
		[ParameterName("wetab"), ParameterInfo("Width dependence of etab")]
		public GivenParameter<double> EtaBW { get; } = new GivenParameter<double>();
		[ParameterName("dl"), ParameterInfo("Channel length reduction in um")]
		public GivenParameter<double> DeltaL { get; } = new GivenParameter<double>();
		[ParameterName("dw"), ParameterInfo("Channel width reduction in um")]
		public GivenParameter<double> DeltaW { get; } = new GivenParameter<double>();
		[ParameterName("mu0"), ParameterInfo("Low-field mobility, at VDS=0 VGS=VTH")]
		public GivenParameter<double> Mob00 { get; } = new GivenParameter<double>(400);
		[ParameterName("mu0b"), ParameterInfo("VBS dependence of low-field mobility")]
		public GivenParameter<double> Mob0B0 { get; } = new GivenParameter<double>();
		[ParameterName("lmu0b"), ParameterInfo("Length dependence of mu0b")]
		public GivenParameter<double> Mob0BL { get; } = new GivenParameter<double>();
		[ParameterName("wmu0b"), ParameterInfo("Width dependence of mu0b")]
		public GivenParameter<double> Mob0BW { get; } = new GivenParameter<double>();
		[ParameterName("mus0"), ParameterInfo("Mobility at VDS=VDD VGS=VTH")]
		public GivenParameter<double> Mobs00 { get; } = new GivenParameter<double>(500);
		[ParameterName("lmus0"), ParameterInfo("Length dependence of mus0")]
		public GivenParameter<double> Mobs0L { get; } = new GivenParameter<double>();
		[ParameterName("wmus0"), ParameterInfo("Width dependence of mus")]
		public GivenParameter<double> Mobs0W { get; } = new GivenParameter<double>();
		[ParameterName("musb"), ParameterInfo("VBS dependence of mus")]
		public GivenParameter<double> MobsB0 { get; } = new GivenParameter<double>();
		[ParameterName("lmusb"), ParameterInfo("Length dependence of musb")]
		public GivenParameter<double> MobsBL { get; } = new GivenParameter<double>();
		[ParameterName("wmusb"), ParameterInfo("Width dependence of musb")]
		public GivenParameter<double> MobsBW { get; } = new GivenParameter<double>();
		[ParameterName("mu20"), ParameterInfo("VDS dependence of mu in tanh term")]
		public GivenParameter<double> Mob200 { get; } = new GivenParameter<double>(1.5);
		[ParameterName("lmu20"), ParameterInfo("Length dependence of mu20")]
		public GivenParameter<double> Mob20L { get; } = new GivenParameter<double>();
		[ParameterName("wmu20"), ParameterInfo("Width dependence of mu20")]
		public GivenParameter<double> Mob20W { get; } = new GivenParameter<double>();
		[ParameterName("mu2b"), ParameterInfo("VBS dependence of mu2")]
		public GivenParameter<double> Mob2B0 { get; } = new GivenParameter<double>();
		[ParameterName("lmu2b"), ParameterInfo("Length dependence of mu2b")]
		public GivenParameter<double> Mob2BL { get; } = new GivenParameter<double>();
		[ParameterName("wmu2b"), ParameterInfo("Width dependence of mu2b")]
		public GivenParameter<double> Mob2BW { get; } = new GivenParameter<double>();
		[ParameterName("mu2g"), ParameterInfo("VGS dependence of mu2")]
		public GivenParameter<double> Mob2G0 { get; } = new GivenParameter<double>();
		[ParameterName("lmu2g"), ParameterInfo("Length dependence of mu2g")]
		public GivenParameter<double> Mob2GL { get; } = new GivenParameter<double>();
		[ParameterName("wmu2g"), ParameterInfo("Width dependence of mu2g")]
		public GivenParameter<double> Mob2GW { get; } = new GivenParameter<double>();
		[ParameterName("mu30"), ParameterInfo("VDS dependence of mu in linear term")]
		public GivenParameter<double> Mob300 { get; } = new GivenParameter<double>(10);
		[ParameterName("lmu30"), ParameterInfo("Length dependence of mu30")]
		public GivenParameter<double> Mob30L { get; } = new GivenParameter<double>();
		[ParameterName("wmu30"), ParameterInfo("Width dependence of mu30")]
		public GivenParameter<double> Mob30W { get; } = new GivenParameter<double>();
		[ParameterName("mu3b"), ParameterInfo("VBS dependence of mu3")]
		public GivenParameter<double> Mob3B0 { get; } = new GivenParameter<double>();
		[ParameterName("lmu3b"), ParameterInfo("Length dependence of mu3b")]
		public GivenParameter<double> Mob3BL { get; } = new GivenParameter<double>();
		[ParameterName("wmu3b"), ParameterInfo("Width dependence of mu3b")]
		public GivenParameter<double> Mob3BW { get; } = new GivenParameter<double>();
		[ParameterName("mu3g"), ParameterInfo("VGS dependence of mu3")]
		public GivenParameter<double> Mob3G0 { get; } = new GivenParameter<double>();
		[ParameterName("lmu3g"), ParameterInfo("Length dependence of mu3g")]
		public GivenParameter<double> Mob3GL { get; } = new GivenParameter<double>();
		[ParameterName("wmu3g"), ParameterInfo("Width dependence of mu3g")]
		public GivenParameter<double> Mob3GW { get; } = new GivenParameter<double>();
		[ParameterName("mu40"), ParameterInfo("VDS dependence of mu in linear term")]
		public GivenParameter<double> Mob400 { get; } = new GivenParameter<double>();
		[ParameterName("lmu40"), ParameterInfo("Length dependence of mu40")]
		public GivenParameter<double> Mob40L { get; } = new GivenParameter<double>();
		[ParameterName("wmu40"), ParameterInfo("Width dependence of mu40")]
		public GivenParameter<double> Mob40W { get; } = new GivenParameter<double>();
		[ParameterName("mu4b"), ParameterInfo("VBS dependence of mu4")]
		public GivenParameter<double> Mob4B0 { get; } = new GivenParameter<double>();
		[ParameterName("lmu4b"), ParameterInfo("Length dependence of mu4b")]
		public GivenParameter<double> Mob4BL { get; } = new GivenParameter<double>();
		[ParameterName("wmu4b"), ParameterInfo("Width dependence of mu4b")]
		public GivenParameter<double> Mob4BW { get; } = new GivenParameter<double>();
		[ParameterName("mu4g"), ParameterInfo("VGS dependence of mu4")]
		public GivenParameter<double> Mob4G0 { get; } = new GivenParameter<double>();
		[ParameterName("lmu4g"), ParameterInfo("Length dependence of mu4g")]
		public GivenParameter<double> Mob4GL { get; } = new GivenParameter<double>();
		[ParameterName("wmu4g"), ParameterInfo("Width dependence of mu4g")]
		public GivenParameter<double> Mob4GW { get; } = new GivenParameter<double>();
		[ParameterName("ua0"), ParameterInfo("Linear VGS dependence of mobility")]
		public GivenParameter<double> Ua00 { get; } = new GivenParameter<double>(0.2);
		[ParameterName("lua0"), ParameterInfo("Length dependence of ua0")]
		public GivenParameter<double> Ua0L { get; } = new GivenParameter<double>();
		[ParameterName("wua0"), ParameterInfo("Width dependence of ua0")]
		public GivenParameter<double> Ua0W { get; } = new GivenParameter<double>();
		[ParameterName("uab"), ParameterInfo("VBS dependence of ua")]
		public GivenParameter<double> UaB0 { get; } = new GivenParameter<double>();
		[ParameterName("luab"), ParameterInfo("Length dependence of uab")]
		public GivenParameter<double> UaBL { get; } = new GivenParameter<double>();
		[ParameterName("wuab"), ParameterInfo("Width dependence of uab")]
		public GivenParameter<double> UaBW { get; } = new GivenParameter<double>();
		[ParameterName("ub0"), ParameterInfo("Quadratic VGS dependence of mobility")]
		public GivenParameter<double> Ub00 { get; } = new GivenParameter<double>();
		[ParameterName("lub0"), ParameterInfo("Length dependence of ub0")]
		public GivenParameter<double> Ub0L { get; } = new GivenParameter<double>();
		[ParameterName("wub0"), ParameterInfo("Width dependence of ub0")]
		public GivenParameter<double> Ub0W { get; } = new GivenParameter<double>();
		[ParameterName("ubb"), ParameterInfo("VBS dependence of ub")]
		public GivenParameter<double> UbB0 { get; } = new GivenParameter<double>();
		[ParameterName("lubb"), ParameterInfo("Length dependence of ubb")]
		public GivenParameter<double> UbBL { get; } = new GivenParameter<double>();
		[ParameterName("wubb"), ParameterInfo("Width dependence of ubb")]
		public GivenParameter<double> UbBW { get; } = new GivenParameter<double>();
		[ParameterName("u10"), ParameterInfo("VDS depence of mobility")]
		public GivenParameter<double> U100 { get; } = new GivenParameter<double>(0.1);
		[ParameterName("lu10"), ParameterInfo("Length dependence of u10")]
		public GivenParameter<double> U10L { get; } = new GivenParameter<double>();
		[ParameterName("wu10"), ParameterInfo("Width dependence of u10")]
		public GivenParameter<double> U10W { get; } = new GivenParameter<double>();
		[ParameterName("u1b"), ParameterInfo("VBS depence of u1")]
		public GivenParameter<double> U1B0 { get; } = new GivenParameter<double>();
		[ParameterName("lu1b"), ParameterInfo("Length depence of u1b")]
		public GivenParameter<double> U1BL { get; } = new GivenParameter<double>();
		[ParameterName("wu1b"), ParameterInfo("Width depence of u1b")]
		public GivenParameter<double> U1BW { get; } = new GivenParameter<double>();
		[ParameterName("u1d"), ParameterInfo("VDS depence of u1")]
		public GivenParameter<double> U1D0 { get; } = new GivenParameter<double>();
		[ParameterName("lu1d"), ParameterInfo("Length depence of u1d")]
		public GivenParameter<double> U1DL { get; } = new GivenParameter<double>();
		[ParameterName("wu1d"), ParameterInfo("Width depence of u1d")]
		public GivenParameter<double> U1DW { get; } = new GivenParameter<double>();
		[ParameterName("n0"), ParameterInfo("Subthreshold slope at VDS=0 VBS=0")]
		public GivenParameter<double> N00 { get; } = new GivenParameter<double>(1.4);
		[ParameterName("ln0"), ParameterInfo("Length dependence of n0")]
		public GivenParameter<double> N0L { get; } = new GivenParameter<double>();
		[ParameterName("wn0"), ParameterInfo("Width dependence of n0")]
		public GivenParameter<double> N0W { get; } = new GivenParameter<double>();
		[ParameterName("nb"), ParameterInfo("VBS dependence of n")]
		public GivenParameter<double> NB0 { get; } = new GivenParameter<double>(0.5);
		[ParameterName("lnb"), ParameterInfo("Length dependence of nb")]
		public GivenParameter<double> NBL { get; } = new GivenParameter<double>();
		[ParameterName("wnb"), ParameterInfo("Width dependence of nb")]
		public GivenParameter<double> NBW { get; } = new GivenParameter<double>();
		[ParameterName("nd"), ParameterInfo("VDS dependence of n")]
		public GivenParameter<double> ND0 { get; } = new GivenParameter<double>();
		[ParameterName("lnd"), ParameterInfo("Length dependence of nd")]
		public GivenParameter<double> NDL { get; } = new GivenParameter<double>();
		[ParameterName("wnd"), ParameterInfo("Width dependence of nd")]
		public GivenParameter<double> NDW { get; } = new GivenParameter<double>();
		[ParameterName("vof0"), ParameterInfo("Threshold voltage offset AT VDS=0 VBS=0")]
		public GivenParameter<double> Vof00 { get; } = new GivenParameter<double>(1.8);
		[ParameterName("lvof0"), ParameterInfo("Length dependence of vof0")]
		public GivenParameter<double> Vof0L { get; } = new GivenParameter<double>();
		[ParameterName("wvof0"), ParameterInfo("Width dependence of vof0")]
		public GivenParameter<double> Vof0W { get; } = new GivenParameter<double>();
		[ParameterName("vofb"), ParameterInfo("VBS dependence of vof")]
		public GivenParameter<double> VofB0 { get; } = new GivenParameter<double>();
		[ParameterName("lvofb"), ParameterInfo("Length dependence of vofb")]
		public GivenParameter<double> VofBL { get; } = new GivenParameter<double>();
		[ParameterName("wvofb"), ParameterInfo("Width dependence of vofb")]
		public GivenParameter<double> VofBW { get; } = new GivenParameter<double>();
		[ParameterName("vofd"), ParameterInfo("VDS dependence of vof")]
		public GivenParameter<double> VofD0 { get; } = new GivenParameter<double>();
		[ParameterName("lvofd"), ParameterInfo("Length dependence of vofd")]
		public GivenParameter<double> VofDL { get; } = new GivenParameter<double>();
		[ParameterName("wvofd"), ParameterInfo("Width dependence of vofd")]
		public GivenParameter<double> VofDW { get; } = new GivenParameter<double>();
		[ParameterName("ai0"), ParameterInfo("Pre-factor of hot-electron effect.")]
		public GivenParameter<double> Ai00 { get; } = new GivenParameter<double>();
		[ParameterName("lai0"), ParameterInfo("Length dependence of ai0")]
		public GivenParameter<double> Ai0L { get; } = new GivenParameter<double>();
		[ParameterName("wai0"), ParameterInfo("Width dependence of ai0")]
		public GivenParameter<double> Ai0W { get; } = new GivenParameter<double>();
		[ParameterName("aib"), ParameterInfo("VBS dependence of ai")]
		public GivenParameter<double> AiB0 { get; } = new GivenParameter<double>();
		[ParameterName("laib"), ParameterInfo("Length dependence of aib")]
		public GivenParameter<double> AiBL { get; } = new GivenParameter<double>();
		[ParameterName("waib"), ParameterInfo("Width dependence of aib")]
		public GivenParameter<double> AiBW { get; } = new GivenParameter<double>();
		[ParameterName("bi0"), ParameterInfo("Exponential factor of hot-electron effect.")]
		public GivenParameter<double> Bi00 { get; } = new GivenParameter<double>();
		[ParameterName("lbi0"), ParameterInfo("Length dependence of bi0")]
		public GivenParameter<double> Bi0L { get; } = new GivenParameter<double>();
		[ParameterName("wbi0"), ParameterInfo("Width dependence of bi0")]
		public GivenParameter<double> Bi0W { get; } = new GivenParameter<double>();
		[ParameterName("bib"), ParameterInfo("VBS dependence of bi")]
		public GivenParameter<double> BiB0 { get; } = new GivenParameter<double>();
		[ParameterName("lbib"), ParameterInfo("Length dependence of bib")]
		public GivenParameter<double> BiBL { get; } = new GivenParameter<double>();
		[ParameterName("wbib"), ParameterInfo("Width dependence of bib")]
		public GivenParameter<double> BiBW { get; } = new GivenParameter<double>();
		[ParameterName("vghigh"), ParameterInfo("Upper bound of the cubic spline function.")]
		public GivenParameter<double> Vghigh0 { get; } = new GivenParameter<double>(0.2);
		[ParameterName("lvghigh"), ParameterInfo("Length dependence of vghigh")]
		public GivenParameter<double> VghighL { get; } = new GivenParameter<double>();
		[ParameterName("wvghigh"), ParameterInfo("Width dependence of vghigh")]
		public GivenParameter<double> VghighW { get; } = new GivenParameter<double>();
		[ParameterName("vglow"), ParameterInfo("Lower bound of the cubic spline function.")]
		public GivenParameter<double> Vglow0 { get; } = new GivenParameter<double>(-0.15);
		[ParameterName("lvglow"), ParameterInfo("Length dependence of vglow")]
		public GivenParameter<double> VglowL { get; } = new GivenParameter<double>();
		[ParameterName("wvglow"), ParameterInfo("Width dependence of vglow")]
		public GivenParameter<double> VglowW { get; } = new GivenParameter<double>();
		[ParameterName("tox"), ParameterInfo("Gate oxide thickness in um")]
		public GivenParameter<double> Tox { get; } = new GivenParameter<double>(0.03);
		[ParameterName("temp"), ParameterInfo("Temperature in degree Celcius")]
		public GivenParameter<double> Temp { get; } = new GivenParameter<double>(27);
		[ParameterName("vdd"), ParameterInfo("Maximum Vds ")]
		public GivenParameter<double> Vdd { get; } = new GivenParameter<double>(5);
		[ParameterName("vgg"), ParameterInfo("Maximum Vgs ")]
		public GivenParameter<double> Vgg { get; } = new GivenParameter<double>(5);
		[ParameterName("vbb"), ParameterInfo("Maximum Vbs ")]
		public GivenParameter<double> Vbb { get; } = new GivenParameter<double>(5);
		[ParameterName("cgso"), ParameterInfo("Gate source overlap capacitance per unit channel width(m)")]
		public GivenParameter<double> GateSourceOverlapCap { get; } = new GivenParameter<double>();
		[ParameterName("cgdo"), ParameterInfo("Gate drain overlap capacitance per unit channel width(m)")]
		public GivenParameter<double> GateDrainOverlapCap { get; } = new GivenParameter<double>();
		[ParameterName("cgbo"), ParameterInfo("Gate bulk overlap capacitance per unit channel length(m)")]
		public GivenParameter<double> GateBulkOverlapCap { get; } = new GivenParameter<double>();
		[ParameterName("xpart"), ParameterInfo("Flag for channel charge partitioning")]
		public GivenParameter<double> ChannelChargePartitionFlag { get; } = new GivenParameter<double>();
		[ParameterName("rsh"), ParameterInfo("Source drain diffusion sheet resistance in ohm per square")]
		public GivenParameter<double> SheetResistance { get; } = new GivenParameter<double>();
		[ParameterName("js"), ParameterInfo("Source drain junction saturation current per unit area")]
		public GivenParameter<double> JctSatCurDensity { get; } = new GivenParameter<double>();
		[ParameterName("pb"), ParameterInfo("Source drain junction built in potential")]
		public GivenParameter<double> BulkJctPotential { get; } = new GivenParameter<double>();
		[ParameterName("mj"), ParameterInfo("Source drain bottom junction capacitance grading coefficient")]
		public GivenParameter<double> BulkJctBotGradingCoeff { get; } = new GivenParameter<double>();
		[ParameterName("pbsw"), ParameterInfo("Source drain side junction capacitance built in potential")]
		public GivenParameter<double> SidewallJctPotential { get; } = new GivenParameter<double>();
		[ParameterName("mjsw"), ParameterInfo("Source drain side junction capacitance grading coefficient")]
		public GivenParameter<double> BulkJctSideGradingCoeff { get; } = new GivenParameter<double>();
		[ParameterName("cj"), ParameterInfo("Source drain bottom junction capacitance per unit area")]
		public GivenParameter<double> UnitAreaJctCap { get; } = new GivenParameter<double>();
		[ParameterName("cjsw"), ParameterInfo("Source drain side junction capacitance per unit area")]
		public GivenParameter<double> UnitLengthSidewallJctCap { get; } = new GivenParameter<double>();
		[ParameterName("wdf"), ParameterInfo("Default width of source drain diffusion in um")]
		public GivenParameter<double> DefaultWidth { get; } = new GivenParameter<double>(10);
		[ParameterName("dell"), ParameterInfo("Length reduction of source drain diffusion")]
		public GivenParameter<double> DeltaLength { get; } = new GivenParameter<double>();
        
	    public double Cox { get; private set; }
	    public double Vdd2 { get; private set; }
	    public double Vgg2 { get; private set; }
	    public double Vbb2 { get; private set; }

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
	        Cox = 3.453e-13 / (Tox * 1.0e-4);
	        Vdd2 = 2.0 * Vdd;
	        Vgg2 = 2.0 * Vgg;
	        Vbb2 = 2.0 * Vbb;
	    }
	}
}