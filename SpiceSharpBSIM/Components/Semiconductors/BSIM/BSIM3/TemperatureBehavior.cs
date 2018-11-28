using System;
using System.IO;
using SpiceSharp.Behaviors;
using SpiceSharp.Simulations;
using SpiceSharp.Simulations.Behaviors;

namespace SpiceSharp.Components.BSIM3Behaviors
{
	/// <summary>
	/// Temperature behavior for a <see cref="BSIM3" />
	/// </summary>
	public class TemperatureBehavior : ExportingBehavior, ITemperatureBehavior
	{
        /// <summary>
        /// Gets the model temperature behavior.
        /// </summary>
        /// <value>
        /// The model temperature behavior.
        /// </value>
        protected ModelTemperatureBehavior ModelTemperature { get; private set; }

        /// <summary>
        /// Gets the base parameters.
        /// </summary>
        /// <value>
        /// The base parameters.
        /// </value>
        protected BaseParameters BaseParameters { get; private set; }

        /// <summary>
        /// Gets the model parameters.
        /// </summary>
        /// <value>
        /// The model parameters.
        /// </value>
        protected ModelBaseParameters ModelParameters { get; private set; }

        /// <summary>
        /// Size dependent parameters
        /// </summary>
	    public BSIM3SizeDependParams Param { get; private set; }

		/// <summary>
		/// Properties
		/// </summary>
		public double Vth0 { get; private set; }
		public double Vfb { get; private set; }
		public double Vfbzb { get; private set; }
		public double U0temp { get; private set; }
		public double Tconst { get; private set; }
		public double DrainConductance { get; private set; }
		public double SourceConductance { get; private set; }
		public double Cgso { get; set; }
		public double Cgdo { get; set; }
		public double Vjsm { get; private set; }
		public double IsEvjsm { get; private set; }
		public double Vjdm { get; private set; }
		public double IsEvjdm { get; private set; }
		
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

            // Get behaviors
			ModelTemperature = provider.GetBehavior<ModelTemperatureBehavior>("model");

            // Get parameters
			BaseParameters = provider.GetParameterSet<BaseParameters>();
			ModelParameters = provider.GetParameterSet<ModelBaseParameters>("model");
		}
		
		/// <summary>
		/// Temperature behavior
		/// </summary>
		public void Temperature(BaseSimulation simulation)
		{
			double tmp, tmp1, tmp2, tmp3, t0 = 0.0, t1, t2, t3, t4, t5, ldrn, wdrn, inv_L, inv_W, inv_LW, nvtm, sourceSatCurrent, drainSatCurrent;

		    var size = new Tuple<double, double>(BaseParameters.Width, BaseParameters.Length);
		    if (!ModelTemperature.SizeDependParams.TryGetValue(size, out var pParam))
            {
                pParam = new BSIM3SizeDependParams();
                Param = pParam;
                ModelTemperature.SizeDependParams.Add(size, pParam);
				ldrn = BaseParameters.Length;
				wdrn = BaseParameters.Width;
				pParam.Length = ldrn;
				pParam.Width = wdrn;
				t0 = Math.Pow(ldrn, ModelParameters.Lln);
				t1 = Math.Pow(wdrn, ModelParameters.Lwn);
				tmp1 = ModelParameters.Ll / t0 + ModelParameters.Lw / t1 + ModelParameters.Lwl / (t0 * t1);
				pParam.BSIM3dl = ModelParameters.Lint + tmp1;
				tmp2 = ModelParameters.Llc / t0 + ModelParameters.Lwc / t1 + ModelParameters.Lwlc / (t0 * t1);
				pParam.BSIM3dlc = ModelParameters.Dlc + tmp2;
				t2 = Math.Pow(ldrn, ModelParameters.Wln);
				t3 = Math.Pow(wdrn, ModelParameters.Wwn);
				tmp1 = ModelParameters.Wl / t2 + ModelParameters.Ww / t3 + ModelParameters.Wwl / (t2 * t3);
				pParam.BSIM3dw = ModelParameters.Wint + tmp1;
				tmp2 = ModelParameters.Wlc / t2 + ModelParameters.Wwc / t3 + ModelParameters.Wwlc / (t2 * t3);
				pParam.BSIM3dwc = ModelParameters.Dwc + tmp2;
				pParam.BSIM3leff = BaseParameters.Length + ModelParameters.Xl - 2.0 * pParam.BSIM3dl;
				if (pParam.BSIM3leff <= 0.0)
                    throw new CircuitException("BSIM3: mosfet {0}, model {1}: Effective channel length <= 0".FormatString(ModelTemperature.Name, Name));
				pParam.BSIM3weff = BaseParameters.Width + ModelParameters.Xw - 2.0 * pParam.BSIM3dw;
				if (pParam.BSIM3weff <= 0.0)
                    throw new CircuitException("BSIM3: mosfet {0}, model {1}: Effective channel width <= 0".FormatString(ModelTemperature.Name, Name));
				pParam.BSIM3leffCV = BaseParameters.Length + ModelParameters.Xl - 2.0 * pParam.BSIM3dlc;
				if (pParam.BSIM3leffCV <= 0.0)
                    throw new CircuitException("BSIM3: mosfet {0}, model {1}: Effective channel length for C-V <= 0".FormatString(ModelTemperature.Name, Name));
				pParam.BSIM3weffCV = BaseParameters.Width + ModelParameters.Xw - 2.0 * pParam.BSIM3dwc;
				if (pParam.BSIM3weffCV <= 0.0)
                    throw new CircuitException("BSIM3: mosfet {0}, model {1}: Effective channel width for C-V <= 0".FormatString(ModelTemperature.Name, Name));
				if (ModelParameters.BinUnit == 1)
				{
					inv_L = 1.0e-6 / pParam.BSIM3leff;
					inv_W = 1.0e-6 / pParam.BSIM3weff;
					inv_LW = 1.0e-12 / (pParam.BSIM3leff * pParam.BSIM3weff);
				}
				else
				{
					inv_L = 1.0 / pParam.BSIM3leff;
					inv_W = 1.0 / pParam.BSIM3weff;
					inv_LW = 1.0 / (pParam.BSIM3leff * pParam.BSIM3weff);
				}
				pParam.BSIM3cdsc = ModelParameters.Cdsc + ModelParameters.Lcdsc * inv_L + ModelParameters.Wcdsc * inv_W + ModelParameters.Pcdsc * inv_LW;
				pParam.BSIM3cdscb = ModelParameters.Cdscb + ModelParameters.Lcdscb * inv_L + ModelParameters.Wcdscb * inv_W + ModelParameters.Pcdscb * inv_LW;
				pParam.BSIM3cdscd = ModelParameters.Cdscd + ModelParameters.Lcdscd * inv_L + ModelParameters.Wcdscd * inv_W + ModelParameters.Pcdscd * inv_LW;
				pParam.BSIM3cit = ModelParameters.Cit + ModelParameters.Lcit * inv_L + ModelParameters.Wcit * inv_W + ModelParameters.Pcit * inv_LW;
				pParam.BSIM3nfactor = ModelParameters.Nfactor + ModelParameters.Lnfactor * inv_L + ModelParameters.Wnfactor * inv_W + ModelParameters.Pnfactor * inv_LW;
				pParam.BSIM3xj = ModelParameters.Xj + ModelParameters.Lxj * inv_L + ModelParameters.Wxj * inv_W + ModelParameters.Pxj * inv_LW;
				pParam.BSIM3vsat = ModelParameters.Vsat + ModelParameters.Lvsat * inv_L + ModelParameters.Wvsat * inv_W + ModelParameters.Pvsat * inv_LW;
				pParam.BSIM3at = ModelParameters.At + ModelParameters.Lat * inv_L + ModelParameters.Wat * inv_W + ModelParameters.Pat * inv_LW;
				pParam.BSIM3a0 = ModelParameters.A0 + ModelParameters.La0 * inv_L + ModelParameters.Wa0 * inv_W + ModelParameters.Pa0 * inv_LW;
				pParam.BSIM3ags = ModelParameters.Ags + ModelParameters.Lags * inv_L + ModelParameters.Wags * inv_W + ModelParameters.Pags * inv_LW;
				pParam.BSIM3a1 = ModelParameters.A1 + ModelParameters.La1 * inv_L + ModelParameters.Wa1 * inv_W + ModelParameters.Pa1 * inv_LW;
				pParam.BSIM3a2 = ModelParameters.A2 + ModelParameters.La2 * inv_L + ModelParameters.Wa2 * inv_W + ModelParameters.Pa2 * inv_LW;
				pParam.BSIM3keta = ModelParameters.Keta + ModelParameters.Lketa * inv_L + ModelParameters.Wketa * inv_W + ModelParameters.Pketa * inv_LW;
				pParam.BSIM3nsub = ModelParameters.Nsub + ModelParameters.Lnsub * inv_L + ModelParameters.Wnsub * inv_W + ModelParameters.Pnsub * inv_LW;
				pParam.BSIM3npeak = ModelParameters.Npeak + ModelParameters.Lnpeak * inv_L + ModelParameters.Wnpeak * inv_W + ModelParameters.Pnpeak * inv_LW;
				pParam.BSIM3ngate = ModelParameters.Ngate + ModelParameters.Lngate * inv_L + ModelParameters.Wngate * inv_W + ModelParameters.Pngate * inv_LW;
				pParam.BSIM3gamma1 = ModelParameters.Gamma1 + ModelParameters.Lgamma1 * inv_L + ModelParameters.Wgamma1 * inv_W + ModelParameters.Pgamma1 * inv_LW;
				pParam.BSIM3gamma2 = ModelParameters.Gamma2 + ModelParameters.Lgamma2 * inv_L + ModelParameters.Wgamma2 * inv_W + ModelParameters.Pgamma2 * inv_LW;
				pParam.BSIM3vbx = ModelParameters.Vbx + ModelParameters.Lvbx * inv_L + ModelParameters.Wvbx * inv_W + ModelParameters.Pvbx * inv_LW;
				pParam.BSIM3vbm = ModelParameters.Vbm + ModelParameters.Lvbm * inv_L + ModelParameters.Wvbm * inv_W + ModelParameters.Pvbm * inv_LW;
				pParam.BSIM3xt = ModelParameters.Xt + ModelParameters.Lxt * inv_L + ModelParameters.Wxt * inv_W + ModelParameters.Pxt * inv_LW;
				pParam.BSIM3vfb = ModelParameters.Vfb + ModelParameters.Lvfb * inv_L + ModelParameters.Wvfb * inv_W + ModelParameters.Pvfb * inv_LW;
				pParam.BSIM3k1 = ModelParameters.K1 + ModelParameters.Lk1 * inv_L + ModelParameters.Wk1 * inv_W + ModelParameters.Pk1 * inv_LW;
				pParam.BSIM3kt1 = ModelParameters.Kt1 + ModelParameters.Lkt1 * inv_L + ModelParameters.Wkt1 * inv_W + ModelParameters.Pkt1 * inv_LW;
				pParam.BSIM3kt1l = ModelParameters.Kt1l + ModelParameters.Lkt1l * inv_L + ModelParameters.Wkt1l * inv_W + ModelParameters.Pkt1l * inv_LW;
				pParam.BSIM3k2 = ModelParameters.K2 + ModelParameters.Lk2 * inv_L + ModelParameters.Wk2 * inv_W + ModelParameters.Pk2 * inv_LW;
				pParam.BSIM3kt2 = ModelParameters.Kt2 + ModelParameters.Lkt2 * inv_L + ModelParameters.Wkt2 * inv_W + ModelParameters.Pkt2 * inv_LW;
				pParam.BSIM3k3 = ModelParameters.K3 + ModelParameters.Lk3 * inv_L + ModelParameters.Wk3 * inv_W + ModelParameters.Pk3 * inv_LW;
				pParam.BSIM3k3b = ModelParameters.K3b + ModelParameters.Lk3b * inv_L + ModelParameters.Wk3b * inv_W + ModelParameters.Pk3b * inv_LW;
				pParam.BSIM3w0 = ModelParameters.W0 + ModelParameters.Lw0 * inv_L + ModelParameters.Ww0 * inv_W + ModelParameters.Pw0 * inv_LW;
				pParam.BSIM3nlx = ModelParameters.Nlx + ModelParameters.Lnlx * inv_L + ModelParameters.Wnlx * inv_W + ModelParameters.Pnlx * inv_LW;
				pParam.BSIM3dvt0 = ModelParameters.Dvt0 + ModelParameters.Ldvt0 * inv_L + ModelParameters.Wdvt0 * inv_W + ModelParameters.Pdvt0 * inv_LW;
				pParam.BSIM3dvt1 = ModelParameters.Dvt1 + ModelParameters.Ldvt1 * inv_L + ModelParameters.Wdvt1 * inv_W + ModelParameters.Pdvt1 * inv_LW;
				pParam.BSIM3dvt2 = ModelParameters.Dvt2 + ModelParameters.Ldvt2 * inv_L + ModelParameters.Wdvt2 * inv_W + ModelParameters.Pdvt2 * inv_LW;
				pParam.BSIM3dvt0w = ModelParameters.Dvt0w + ModelParameters.Ldvt0w * inv_L + ModelParameters.Wdvt0w * inv_W + ModelParameters.Pdvt0w * inv_LW;
				pParam.BSIM3dvt1w = ModelParameters.Dvt1w + ModelParameters.Ldvt1w * inv_L + ModelParameters.Wdvt1w * inv_W + ModelParameters.Pdvt1w * inv_LW;
				pParam.BSIM3dvt2w = ModelParameters.Dvt2w + ModelParameters.Ldvt2w * inv_L + ModelParameters.Wdvt2w * inv_W + ModelParameters.Pdvt2w * inv_LW;
				pParam.BSIM3drout = ModelParameters.Drout + ModelParameters.Ldrout * inv_L + ModelParameters.Wdrout * inv_W + ModelParameters.Pdrout * inv_LW;
				pParam.BSIM3dsub = ModelParameters.Dsub + ModelParameters.Ldsub * inv_L + ModelParameters.Wdsub * inv_W + ModelParameters.Pdsub * inv_LW;
				pParam.BSIM3vth0 = ModelParameters.Vth0 + ModelParameters.Lvth0 * inv_L + ModelParameters.Wvth0 * inv_W + ModelParameters.Pvth0 * inv_LW;
				pParam.BSIM3ua = ModelParameters.Ua + ModelParameters.Lua * inv_L + ModelParameters.Wua * inv_W + ModelParameters.Pua * inv_LW;
				pParam.BSIM3ua1 = ModelParameters.Ua1 + ModelParameters.Lua1 * inv_L + ModelParameters.Wua1 * inv_W + ModelParameters.Pua1 * inv_LW;
				pParam.BSIM3ub = ModelParameters.Ub + ModelParameters.Lub * inv_L + ModelParameters.Wub * inv_W + ModelParameters.Pub * inv_LW;
				pParam.BSIM3ub1 = ModelParameters.Ub1 + ModelParameters.Lub1 * inv_L + ModelParameters.Wub1 * inv_W + ModelParameters.Pub1 * inv_LW;
				pParam.BSIM3uc = ModelParameters.Uc + ModelParameters.Luc * inv_L + ModelParameters.Wuc * inv_W + ModelParameters.Puc * inv_LW;
				pParam.BSIM3uc1 = ModelParameters.Uc1 + ModelParameters.Luc1 * inv_L + ModelParameters.Wuc1 * inv_W + ModelParameters.Puc1 * inv_LW;
				pParam.BSIM3u0 = ModelParameters.U0 + ModelParameters.Lu0 * inv_L + ModelParameters.Wu0 * inv_W + ModelParameters.Pu0 * inv_LW;
				pParam.BSIM3ute = ModelParameters.Ute + ModelParameters.Lute * inv_L + ModelParameters.Wute * inv_W + ModelParameters.Pute * inv_LW;
				pParam.BSIM3voff = ModelParameters.Voff + ModelParameters.Lvoff * inv_L + ModelParameters.Wvoff * inv_W + ModelParameters.Pvoff * inv_LW;
				pParam.BSIM3delta = ModelParameters.Delta + ModelParameters.Ldelta * inv_L + ModelParameters.Wdelta * inv_W + ModelParameters.Pdelta * inv_LW;
				pParam.BSIM3rdsw = ModelParameters.Rdsw + ModelParameters.Lrdsw * inv_L + ModelParameters.Wrdsw * inv_W + ModelParameters.Prdsw * inv_LW;
				pParam.BSIM3prwg = ModelParameters.Prwg + ModelParameters.Lprwg * inv_L + ModelParameters.Wprwg * inv_W + ModelParameters.Pprwg * inv_LW;
				pParam.BSIM3prwb = ModelParameters.Prwb + ModelParameters.Lprwb * inv_L + ModelParameters.Wprwb * inv_W + ModelParameters.Pprwb * inv_LW;
				pParam.BSIM3prt = ModelParameters.Prt + ModelParameters.Lprt * inv_L + ModelParameters.Wprt * inv_W + ModelParameters.Pprt * inv_LW;
				pParam.BSIM3eta0 = ModelParameters.Eta0 + ModelParameters.Leta0 * inv_L + ModelParameters.Weta0 * inv_W + ModelParameters.Peta0 * inv_LW;
				pParam.BSIM3etab = ModelParameters.Etab + ModelParameters.Letab * inv_L + ModelParameters.Wetab * inv_W + ModelParameters.Petab * inv_LW;
				pParam.BSIM3pclm = ModelParameters.Pclm + ModelParameters.Lpclm * inv_L + ModelParameters.Wpclm * inv_W + ModelParameters.Ppclm * inv_LW;
				pParam.BSIM3pdibl1 = ModelParameters.Pdibl1 + ModelParameters.Lpdibl1 * inv_L + ModelParameters.Wpdibl1 * inv_W + ModelParameters.Ppdibl1 * inv_LW;
				pParam.BSIM3pdibl2 = ModelParameters.Pdibl2 + ModelParameters.Lpdibl2 * inv_L + ModelParameters.Wpdibl2 * inv_W + ModelParameters.Ppdibl2 * inv_LW;
				pParam.BSIM3pdiblb = ModelParameters.Pdiblb + ModelParameters.Lpdiblb * inv_L + ModelParameters.Wpdiblb * inv_W + ModelParameters.Ppdiblb * inv_LW;
				pParam.BSIM3pscbe1 = ModelParameters.Pscbe1 + ModelParameters.Lpscbe1 * inv_L + ModelParameters.Wpscbe1 * inv_W + ModelParameters.Ppscbe1 * inv_LW;
				pParam.BSIM3pscbe2 = ModelParameters.Pscbe2 + ModelParameters.Lpscbe2 * inv_L + ModelParameters.Wpscbe2 * inv_W + ModelParameters.Ppscbe2 * inv_LW;
				pParam.BSIM3pvag = ModelParameters.Pvag + ModelParameters.Lpvag * inv_L + ModelParameters.Wpvag * inv_W + ModelParameters.Ppvag * inv_LW;
				pParam.BSIM3wr = ModelParameters.Wr + ModelParameters.Lwr * inv_L + ModelParameters.Wwr * inv_W + ModelParameters.Pwr * inv_LW;
				pParam.BSIM3dwg = ModelParameters.Dwg + ModelParameters.Ldwg * inv_L + ModelParameters.Wdwg * inv_W + ModelParameters.Pdwg * inv_LW;
				pParam.BSIM3dwb = ModelParameters.Dwb + ModelParameters.Ldwb * inv_L + ModelParameters.Wdwb * inv_W + ModelParameters.Pdwb * inv_LW;
				pParam.BSIM3b0 = ModelParameters.B0 + ModelParameters.Lb0 * inv_L + ModelParameters.Wb0 * inv_W + ModelParameters.Pb0 * inv_LW;
				pParam.BSIM3b1 = ModelParameters.B1 + ModelParameters.Lb1 * inv_L + ModelParameters.Wb1 * inv_W + ModelParameters.Pb1 * inv_LW;
				pParam.BSIM3alpha0 = ModelParameters.Alpha0 + ModelParameters.Lalpha0 * inv_L + ModelParameters.Walpha0 * inv_W + ModelParameters.Palpha0 * inv_LW;
				pParam.BSIM3alpha1 = ModelParameters.Alpha1 + ModelParameters.Lalpha1 * inv_L + ModelParameters.Walpha1 * inv_W + ModelParameters.Palpha1 * inv_LW;
				pParam.BSIM3beta0 = ModelParameters.Beta0 + ModelParameters.Lbeta0 * inv_L + ModelParameters.Wbeta0 * inv_W + ModelParameters.Pbeta0 * inv_LW;
				pParam.BSIM3elm = ModelParameters.Elm + ModelParameters.Lelm * inv_L + ModelParameters.Welm * inv_W + ModelParameters.Pelm * inv_LW;
				pParam.BSIM3cgsl = ModelParameters.Cgsl + ModelParameters.Lcgsl * inv_L + ModelParameters.Wcgsl * inv_W + ModelParameters.Pcgsl * inv_LW;
				pParam.BSIM3cgdl = ModelParameters.Cgdl + ModelParameters.Lcgdl * inv_L + ModelParameters.Wcgdl * inv_W + ModelParameters.Pcgdl * inv_LW;
				pParam.BSIM3ckappa = ModelParameters.Ckappa + ModelParameters.Lckappa * inv_L + ModelParameters.Wckappa * inv_W + ModelParameters.Pckappa * inv_LW;
				pParam.BSIM3cf = ModelParameters.Cf + ModelParameters.Lcf * inv_L + ModelParameters.Wcf * inv_W + ModelParameters.Pcf * inv_LW;
				pParam.BSIM3clc = ModelParameters.Clc + ModelParameters.Lclc * inv_L + ModelParameters.Wclc * inv_W + ModelParameters.Pclc * inv_LW;
				pParam.BSIM3cle = ModelParameters.Cle + ModelParameters.Lcle * inv_L + ModelParameters.Wcle * inv_W + ModelParameters.Pcle * inv_LW;
				pParam.BSIM3vfbcv = ModelParameters.Vfbcv + ModelParameters.Lvfbcv * inv_L + ModelParameters.Wvfbcv * inv_W + ModelParameters.Pvfbcv * inv_LW;
				pParam.BSIM3acde = ModelParameters.Acde + ModelParameters.Lacde * inv_L + ModelParameters.Wacde * inv_W + ModelParameters.Pacde * inv_LW;
				pParam.BSIM3moin = ModelParameters.Moin + ModelParameters.Lmoin * inv_L + ModelParameters.Wmoin * inv_W + ModelParameters.Pmoin * inv_LW;
				pParam.BSIM3noff = ModelParameters.Noff + ModelParameters.Lnoff * inv_L + ModelParameters.Wnoff * inv_W + ModelParameters.Pnoff * inv_LW;
				pParam.BSIM3voffcv = ModelParameters.Voffcv + ModelParameters.Lvoffcv * inv_L + ModelParameters.Wvoffcv * inv_W + ModelParameters.Pvoffcv * inv_LW;
				pParam.BSIM3abulkCVfactor = 1.0 + Math.Pow(pParam.BSIM3clc / pParam.BSIM3leffCV, pParam.BSIM3cle);
				t0 = ModelTemperature.TRatio - 1.0;
				pParam.BSIM3ua = pParam.BSIM3ua + pParam.BSIM3ua1 * t0;
				pParam.BSIM3ub = pParam.BSIM3ub + pParam.BSIM3ub1 * t0;
				pParam.BSIM3uc = pParam.BSIM3uc + pParam.BSIM3uc1 * t0;
				if (pParam.BSIM3u0 > 1.0)
					pParam.BSIM3u0 = pParam.BSIM3u0 / 1.0e4;
				pParam.BSIM3u0temp = pParam.BSIM3u0 * Math.Pow(ModelTemperature.TRatio, pParam.BSIM3ute);
				pParam.BSIM3vsattemp = pParam.BSIM3vsat - pParam.BSIM3at * t0;
				pParam.BSIM3rds0 = (pParam.BSIM3rdsw + pParam.BSIM3prt * t0) / Math.Pow(pParam.BSIM3weff * 1E6, pParam.BSIM3wr);
				if (Check())
                    throw new CircuitException("Fatal error(s) detected during BSIM3V3.3 parameter checking for {0} in model {1}".FormatString(Name, ModelTemperature.Name));
				pParam.BSIM3cgdo = (ModelParameters.Cgdo + pParam.BSIM3cf) * pParam.BSIM3weffCV;
				pParam.BSIM3cgso = (ModelParameters.Cgso + pParam.BSIM3cf) * pParam.BSIM3weffCV;
				pParam.BSIM3cgbo = ModelParameters.Cgbo * pParam.BSIM3leffCV;
				t0 = pParam.BSIM3leffCV * pParam.BSIM3leffCV;
				pParam.BSIM3tconst = pParam.BSIM3u0temp * pParam.BSIM3elm / (ModelParameters.Cox * pParam.BSIM3weffCV * pParam.BSIM3leffCV * t0);
				if (!ModelParameters.Npeak.Given && ModelParameters.Gamma1.Given)
				{
					t0 = pParam.BSIM3gamma1 * ModelParameters.Cox;
					pParam.BSIM3npeak = 3.021E22 * t0 * t0;
				}
				pParam.BSIM3phi = 2.0 * ModelTemperature.Vtm0 * Math.Log(pParam.BSIM3npeak / ModelTemperature.Ni);
				pParam.BSIM3sqrtPhi = Math.Sqrt(pParam.BSIM3phi);
				pParam.BSIM3phis3 = pParam.BSIM3sqrtPhi * pParam.BSIM3phi;
				pParam.BSIM3Xdep0 = Math.Sqrt(2.0 * 1.03594e-10 / (1.60219e-19 * pParam.BSIM3npeak * 1.0e6)) * pParam.BSIM3sqrtPhi;
				pParam.BSIM3sqrtXdep0 = Math.Sqrt(pParam.BSIM3Xdep0);
				pParam.BSIM3litl = Math.Sqrt(3.0 * pParam.BSIM3xj * ModelParameters.Tox);
				pParam.BSIM3vbi = ModelTemperature.Vtm0 * Math.Log(1.0e20 * pParam.BSIM3npeak / (ModelTemperature.Ni * ModelTemperature.Ni));
				pParam.BSIM3cdep0 = Math.Sqrt(1.60219e-19 * 1.03594e-10 * pParam.BSIM3npeak * 1.0e6 / 2.0 / pParam.BSIM3phi);
				pParam.BSIM3ldeb = Math.Sqrt(1.03594e-10 * ModelTemperature.Vtm0 / (1.60219e-19 * pParam.BSIM3npeak * 1.0e6)) / 3.0;
				pParam.BSIM3acde *= Math.Pow(pParam.BSIM3npeak / 2.0e16, -0.25);
                if (ModelParameters.K1.Given || ModelParameters.K2.Given)
                {
                    if (!ModelParameters.K1.Given)
                    {
                        CircuitWarning.Warning(this, "Warning: k1 should be specified with k2.");
                        pParam.BSIM3k1 = 0.53;
                    }

                    if (!ModelParameters.K2.Given)
                    {
                        CircuitWarning.Warning(this, "Warning: k2 should be specified with k1.");
                        pParam.BSIM3k2 = -0.0186;
                    }

                    if (ModelParameters.Nsub.Given)
                        CircuitWarning.Warning(this, "Warning: nsub is ignored because k1 or k2 is given.");
                    if (ModelParameters.Xt.Given)
                        CircuitWarning.Warning(this, "Warning: xt is ignored because k1 or k2 is given.");
                    if (ModelParameters.Vbx.Given)
                        CircuitWarning.Warning(this, "Warning: vbx is ignored because k1 or k2 is given.");
                    if (ModelParameters.Gamma1.Given)
                        CircuitWarning.Warning(this, "Warning: gamma1 is ignored because k1 or k2 is given.");
                    if (ModelParameters.Gamma2.Given)
                        CircuitWarning.Warning(this, "Warning: gamma2 is ignored because k1 or k2 is given.");
                }
                else
				{
					if (!ModelParameters.Vbx.Given)
						pParam.BSIM3vbx = pParam.BSIM3phi - 7.7348e-4 * pParam.BSIM3npeak * pParam.BSIM3xt * pParam.BSIM3xt;
					if (pParam.BSIM3vbx > 0.0)
						pParam.BSIM3vbx = -pParam.BSIM3vbx;
					if (pParam.BSIM3vbm > 0.0)
						pParam.BSIM3vbm = -pParam.BSIM3vbm;
					if (!ModelParameters.Gamma1.Given)
						pParam.BSIM3gamma1 = 5.753e-12 * Math.Sqrt(pParam.BSIM3npeak) / ModelParameters.Cox;
					if (!ModelParameters.Gamma2.Given)
						pParam.BSIM3gamma2 = 5.753e-12 * Math.Sqrt(pParam.BSIM3nsub) / ModelParameters.Cox;
					t0 = pParam.BSIM3gamma1 - pParam.BSIM3gamma2;
					t1 = Math.Sqrt(pParam.BSIM3phi - pParam.BSIM3vbx) - pParam.BSIM3sqrtPhi;
					t2 = Math.Sqrt(pParam.BSIM3phi * (pParam.BSIM3phi - pParam.BSIM3vbm)) - pParam.BSIM3phi;
					pParam.BSIM3k2 = t0 * t1 / (2.0 * t2 + pParam.BSIM3vbm);
					pParam.BSIM3k1 = pParam.BSIM3gamma2 - 2.0 * pParam.BSIM3k2 * Math.Sqrt(pParam.BSIM3phi - pParam.BSIM3vbm);
				}
				if (pParam.BSIM3k2 < 0.0)
				{
					t0 = 0.5 * pParam.BSIM3k1 / pParam.BSIM3k2;
					pParam.BSIM3vbsc = 0.9 * (pParam.BSIM3phi - t0 * t0);
					if (pParam.BSIM3vbsc > -3.0)
						pParam.BSIM3vbsc = -3.0;
					else if (pParam.BSIM3vbsc < -30.0)
						pParam.BSIM3vbsc = -30.0;
				}
				else
				{
					pParam.BSIM3vbsc = -30.0;
				}
				if (pParam.BSIM3vbsc > pParam.BSIM3vbm)
					pParam.BSIM3vbsc = pParam.BSIM3vbm;
				if (!ModelParameters.Vfb.Given)
				{
					if (ModelParameters.Vth0.Given)
					{
						pParam.BSIM3vfb = ModelParameters.B3Type * pParam.BSIM3vth0 - pParam.BSIM3phi - pParam.BSIM3k1 * pParam.BSIM3sqrtPhi;
					}
					else
					{
						pParam.BSIM3vfb = -1.0;
					}
				}
				if (!ModelParameters.Vth0.Given)
				{
					pParam.BSIM3vth0 = ModelParameters.B3Type * (pParam.BSIM3vfb + pParam.BSIM3phi + pParam.BSIM3k1 * pParam.BSIM3sqrtPhi);
				}
				pParam.BSIM3k1ox = pParam.BSIM3k1 * ModelParameters.Tox / ModelParameters.Toxm;
				pParam.BSIM3k2ox = pParam.BSIM3k2 * ModelParameters.Tox / ModelParameters.Toxm;
				t1 = Math.Sqrt(1.03594e-10 / 3.453133e-11 * ModelParameters.Tox * pParam.BSIM3Xdep0);
				t0 = Math.Exp(-0.5 * pParam.BSIM3dsub * pParam.BSIM3leff / t1);
				pParam.BSIM3theta0vb0 = t0 + 2.0 * t0 * t0;
				t0 = Math.Exp(-0.5 * pParam.BSIM3drout * pParam.BSIM3leff / t1);
				t2 = t0 + 2.0 * t0 * t0;
				pParam.BSIM3thetaRout = pParam.BSIM3pdibl1 * t2 + pParam.BSIM3pdibl2;
				tmp = Math.Sqrt(pParam.BSIM3Xdep0);
				tmp1 = pParam.BSIM3vbi - pParam.BSIM3phi;
				tmp2 = ModelTemperature.Factor1 * tmp;
				t0 = -0.5 * pParam.BSIM3dvt1w * pParam.BSIM3weff * pParam.BSIM3leff / tmp2;
				if (t0 > -34.0)
				{
					t1 = Math.Exp(t0);
					t2 = t1 * (1.0 + 2.0 * t1);
				}
				else
				{
					t1 = 1.713908431e-15;
					t2 = t1 * (1.0 + 2.0 * t1);
				}
				t0 = pParam.BSIM3dvt0w * t2;
				t2 = t0 * tmp1;
				t0 = -0.5 * pParam.BSIM3dvt1 * pParam.BSIM3leff / tmp2;
				if (t0 > -34.0)
				{
					t1 = Math.Exp(t0);
					t3 = t1 * (1.0 + 2.0 * t1);
				}
				else
				{
					t1 = 1.713908431e-15;
					t3 = t1 * (1.0 + 2.0 * t1);
				}
				t3 = pParam.BSIM3dvt0 * t3 * tmp1;
				t4 = ModelParameters.Tox * pParam.BSIM3phi / (pParam.BSIM3weff + pParam.BSIM3w0);
				t0 = Math.Sqrt(1.0 + pParam.BSIM3nlx / pParam.BSIM3leff);
				t5 = pParam.BSIM3k1ox * (t0 - 1.0) * pParam.BSIM3sqrtPhi + (pParam.BSIM3kt1 + pParam.BSIM3kt1l / pParam.BSIM3leff) * (ModelTemperature.TRatio - 1.0);
				tmp3 = ModelParameters.B3Type * pParam.BSIM3vth0 - t2 - t3 + pParam.BSIM3k3 * t4 + t5;
				pParam.BSIM3vfbzb = tmp3 - pParam.BSIM3phi - pParam.BSIM3k1 * pParam.BSIM3sqrtPhi;
			}

		    Param = pParam;
			Vth0 = Param.BSIM3vth0 + BaseParameters.Delvto;
			Vfb = Param.BSIM3vfb + ModelParameters.B3Type * BaseParameters.Delvto;
			Vfbzb = Param.BSIM3vfbzb + ModelParameters.B3Type * BaseParameters.Delvto;
			U0temp = Param.BSIM3u0temp * BaseParameters.Mulu0;
			Tconst = U0temp * Param.BSIM3elm / (ModelParameters.Cox * Param.BSIM3weffCV * Param.BSIM3leffCV * t0);

		    var drainResistance = ModelParameters.SheetResistance * BaseParameters.DrainSquares;
		    var sourceResistance = ModelParameters.SheetResistance * BaseParameters.SourceSquares;

			if (drainResistance > 0.0)
				DrainConductance = 1.0 / drainResistance;
			else
				DrainConductance = 0.0;
			if (sourceResistance > 0.0)
				SourceConductance = 1.0 / sourceResistance;
			else
				SourceConductance = 0.0;
			Cgso = Param.BSIM3cgso;
			Cgdo = Param.BSIM3cgdo;
			nvtm = ModelTemperature.Vtm * ModelParameters.JctEmissionCoeff;
			if (BaseParameters.SourceArea <= 0.0 && BaseParameters.SourcePerimeter <= 0.0)
			{
				sourceSatCurrent = 1.0e-14;
			}
			else
			{
				sourceSatCurrent = BaseParameters.SourceArea * ModelTemperature.JctTempSatCurDensity + BaseParameters.SourcePerimeter * ModelTemperature.JctSidewallTempSatCurDensity;
			}
			if (sourceSatCurrent > 0.0 && ModelParameters.Ijth > 0.0)
			{
				Vjsm = nvtm * Math.Log(ModelParameters.Ijth / sourceSatCurrent + 1.0);
				IsEvjsm = sourceSatCurrent * Math.Exp(Vjsm / nvtm);
			}
			if (BaseParameters.DrainArea <= 0.0 && BaseParameters.DrainPerimeter <= 0.0)
			{
				drainSatCurrent = 1.0e-14;
			}
			else
			{
				drainSatCurrent = BaseParameters.DrainArea * ModelTemperature.JctTempSatCurDensity + BaseParameters.DrainPerimeter * ModelTemperature.JctSidewallTempSatCurDensity;
			}
			if (drainSatCurrent > 0.0 && ModelParameters.Ijth > 0.0)
			{
				Vjdm = nvtm * Math.Log(ModelParameters.Ijth / drainSatCurrent + 1.0);
				IsEvjdm = drainSatCurrent * Math.Exp(Vjdm / nvtm);
			}
		}

        /// <summary>
        /// Check the model and instance parameters
        /// </summary>
        /// <returns></returns>
	    private bool Check()
        {
            var fatal = false;
            using (var sw = new StreamWriter(ModelParameters.CheckPath))
            {
                sw.WriteLine("BSIM3v3.3.0 Parameter Checking.");
                if (ModelParameters.Version != "3.3.0" && ModelParameters.Version != "3.30" && ModelParameters.Version != "3.3")
                {
                    sw.WriteLine("Warning: This model is BSIM3v3.3.0; you specified a wrong version number.");
                    CircuitWarning.Warning(this, "Warning: This model is BSIM3v3.3.0; you specified a wrong version number.");
                }
                sw.WriteLine("Model = {0}".FormatString(ModelTemperature.Name));

                if (Param.BSIM3nlx < -Param.BSIM3leff)
                {
                    sw.WriteLine("Fatal: Nlx = {0} is less than -Leff.".FormatString(Param.BSIM3nlx));
                    CircuitWarning.Warning(this, "Fatal: Nlx = {0} is less than -Leff.".FormatString(Param.BSIM3nlx));
                    fatal = true;
                }

                if (ModelParameters.Tox <= 0.0)
                {
                    sw.WriteLine("Fatal: Tox = {0} is not positive.".FormatString(ModelParameters.Tox));
                    CircuitWarning.Warning(this, "Fatal: Tox = {0} is not positive.".FormatString(ModelParameters.Tox));
                    fatal = true;
                }

                if (ModelParameters.Toxm <= 0.0)
                {
                    sw.WriteLine("Fatal: Toxm = {0} is not positive.".FormatString(ModelParameters.Toxm));
                    CircuitWarning.Warning(this, "Fatal: Toxm = {0} is not positive.".FormatString(ModelParameters.Toxm));
                    fatal = true;
                }

                if (ModelParameters.Lintnoi > Param.BSIM3leff / 2)
                {
                    sw.WriteLine("Fatal: Lintnoi = {0} is too large - Leff for noise is negative.".FormatString(ModelParameters.Lintnoi));
                    CircuitWarning.Warning(this,
                        "Fatal: Lintnoi = {0} is too large - Leff for noise is negative.".FormatString(ModelParameters.Lintnoi));
                    fatal = true;
                }

                if (Param.BSIM3npeak <= 0.0)
                {
                    sw.WriteLine("Fatal: Nch = {0} is not positive.".FormatString(Param.BSIM3npeak));
                    CircuitWarning.Warning(this, "Fatal: Nch = {0} is not positive.".FormatString(Param.BSIM3npeak));
                    fatal = true;
                }
                if (Param.BSIM3nsub <= 0.0)
                {
                    sw.WriteLine("Fatal: Nsub = {0} is not positive.".FormatString(Param.BSIM3nsub));
                    CircuitWarning.Warning(this, "Fatal: Nsub = {0} is not positive.".FormatString(Param.BSIM3nsub));
                    fatal = true;
                }
                if (Param.BSIM3ngate < 0.0)
                {
                    sw.WriteLine("Fatal: Ngate = {0} is not positive.".FormatString(Param.BSIM3ngate));
                    CircuitWarning.Warning(this, "Fatal: Ngate = {0} Ngate is not positive.".FormatString(Param.BSIM3ngate));
                    fatal = true;
                }
                if (Param.BSIM3ngate > 1.0e25)
                {
                    sw.WriteLine("Fatal: Ngate = {0} is too high.".FormatString(Param.BSIM3ngate));
                    CircuitWarning.Warning(this, "Fatal: Ngate = {0} Ngate is too high".FormatString(Param.BSIM3ngate));
                    fatal = true;
                }
                if (Param.BSIM3xj <= 0.0)
                {
                    sw.WriteLine("Fatal: Xj = {0} is not positive.".FormatString(Param.BSIM3xj));
                    CircuitWarning.Warning(this, "Fatal: Xj = {0} is not positive.".FormatString(Param.BSIM3xj));
                    fatal = true;
                }

                if (Param.BSIM3dvt1 < 0.0)
                {
                    sw.WriteLine("Fatal: Dvt1 = {0} is negative.".FormatString(Param.BSIM3dvt1));
                    CircuitWarning.Warning(this, "Fatal: Dvt1 = {0} is negative.".FormatString(Param.BSIM3dvt1));
                    fatal = true;
                }

                if (Param.BSIM3dvt1w < 0.0)
                {
                    sw.WriteLine("Fatal: Dvt1w = {0} is negative.".FormatString(Param.BSIM3dvt1w));
                    CircuitWarning.Warning(this, "Fatal: Dvt1w = {0} is negative.".FormatString(Param.BSIM3dvt1w));
                    fatal = true;
                }

                if (Param.BSIM3w0 == -Param.BSIM3weff)
                {
                    sw.WriteLine("Fatal: (W0 + Weff) = 0 causing divided-by-zero.");
                    CircuitWarning.Warning(this, "Fatal: (W0 + Weff) = 0 causing divided-by-zero.");
                    fatal = true;
                }

                if (Param.BSIM3dsub < 0.0)
                {
                    sw.WriteLine("Fatal: Dsub = {0} is negative.".FormatString(Param.BSIM3dsub));
                    CircuitWarning.Warning(this, "Fatal: Dsub = {0} is negative.".FormatString(Param.BSIM3dsub));
                    fatal = true;
                }
                if (Param.BSIM3b1 == -Param.BSIM3weff)
                {
                    sw.WriteLine("Fatal: (B1 + Weff) = 0 causing divided-by-zero.");
                    CircuitWarning.Warning(this, "Fatal: (B1 + Weff) = 0 causing divided-by-zero.");
                    fatal = true;
                }
                if (Param.BSIM3u0temp <= 0.0)
                {
                    sw.WriteLine("Fatal: u0 at current temperature = {0} is not positive.".FormatString(Param.BSIM3u0temp));
                    CircuitWarning.Warning(this, "Fatal: u0 at current temperature = {0} is not positive.".FormatString(Param.BSIM3u0temp));
                    fatal = true;
                }

                /* Check delta parameter */
                if (Param.BSIM3delta < 0.0)
                {
                    sw.WriteLine("Fatal: Delta = {0} is less than zero.".FormatString(Param.BSIM3delta));
                    CircuitWarning.Warning(this, "Fatal: Delta = {0} is less than zero.".FormatString(Param.BSIM3delta));
                    fatal = true;
                }

                if (Param.BSIM3vsattemp <= 0.0)
                {
                    sw.WriteLine("Fatal: Vsat at current temperature = {0} is not positive.".FormatString(Param.BSIM3vsattemp));
                    CircuitWarning.Warning(this, "Fatal: Vsat at current temperature = {0} is not positive.".FormatString(Param.BSIM3vsattemp));
                    fatal = true;
                }
                /* Check Rout parameters */
                if (Param.BSIM3pclm <= 0.0)
                {
                    sw.WriteLine("Fatal: Pclm = {0} is not positive.".FormatString(Param.BSIM3pclm));
                    CircuitWarning.Warning(this, "Fatal: Pclm = {0} is not positive.".FormatString(Param.BSIM3pclm));
                    fatal = true;
                }

                if (Param.BSIM3drout < 0.0)
                {
                    sw.WriteLine("Fatal: Drout = {0} is negative.".FormatString(Param.BSIM3drout));
                    CircuitWarning.Warning(this, "Fatal: Drout = {0} is negative.".FormatString(Param.BSIM3drout));
                    fatal = true;
                }

                if (Param.BSIM3pscbe2 <= 0.0)
                {
                    sw.WriteLine("Warning: Pscbe2 = {0} is not positive.".FormatString(Param.BSIM3pscbe2));
                    CircuitWarning.Warning(this, "Warning: Pscbe2 = {0} is not positive.".FormatString(Param.BSIM3pscbe2));
                }

                if (Param.BSIM3noff < 0.1)
                {
                    sw.WriteLine("Warning: Noff = {0} is too small.".FormatString(Param.BSIM3noff));
                    CircuitWarning.Warning(this, "Warning: Noff = {0} is too small.".FormatString(Param.BSIM3noff));
                }
                if (Param.BSIM3noff > 4.0)
                {
                    sw.WriteLine("Warning: Noff = {0} is too large.".FormatString(Param.BSIM3noff));
                    CircuitWarning.Warning(this, "Warning: Noff = {0} is too large.".FormatString(Param.BSIM3noff));
                }

                if (Param.BSIM3voffcv < -0.5)
                {
                    sw.WriteLine("Warning: Voffcv = {0} is too small.".FormatString(Param.BSIM3voffcv));
                    CircuitWarning.Warning(this, "Warning: Voffcv = {0} is too small.".FormatString(Param.BSIM3voffcv));
                }
                if (Param.BSIM3voffcv > 0.5)
                {
                    sw.WriteLine("Warning: Voffcv = {0} is too large.".FormatString(Param.BSIM3voffcv));
                    CircuitWarning.Warning(this, "Warning: Voffcv = {0} is too large.".FormatString(Param.BSIM3voffcv));
                }

                if (ModelParameters.Ijth < 0.0)
                {
                    sw.WriteLine("Fatal: Ijth = {0} cannot be negative.".FormatString(ModelParameters.Ijth));
                    CircuitWarning.Warning(this, "Fatal: Ijth = {0} cannot be negative.".FormatString(ModelParameters.Ijth));
                    fatal = true;
                }

                /* Check capacitance parameters */
                if (Param.BSIM3clc < 0.0)
                {
                    sw.WriteLine("Fatal: Clc = {0} is negative.".FormatString(Param.BSIM3clc));
                    CircuitWarning.Warning(this, "Fatal: Clc = {0} is negative.".FormatString(Param.BSIM3clc));
                    fatal = true;
                }

                if (Param.BSIM3moin < 5.0)
                {
                    sw.WriteLine("Warning: Moin = {0} is too small.".FormatString(Param.BSIM3moin));
                    CircuitWarning.Warning(this, "Warning: Moin = {0} is too small.".FormatString(Param.BSIM3moin));
                }
                if (Param.BSIM3moin > 25.0)
                {
                    sw.WriteLine("Warning: Moin = {0} is too large.".FormatString(Param.BSIM3moin));
                    CircuitWarning.Warning(this, "Warning: Moin = {0} is too large.".FormatString(Param.BSIM3moin));
                }

                if (ModelParameters.CapMod == 3)
                {
                    if (Param.BSIM3acde < 0.4)
                    {
                        sw.WriteLine("Warning:  Acde = {0} is too small.".FormatString(Param.BSIM3acde));
                        CircuitWarning.Warning(this, "Warning: Acde = {0} is too small.".FormatString(Param.BSIM3acde));
                    }
                    if (Param.BSIM3acde > 1.6)
                    {
                        sw.WriteLine("Warning:  Acde = {0} is too large.".FormatString(Param.BSIM3acde));
                        CircuitWarning.Warning(this, "Warning: Acde = {0} is too large.".FormatString(Param.BSIM3acde));
                    }
                }

                if (ModelParameters.ParamChk == 1)
                {
                    /* Check L and W parameters */
                    if (Param.BSIM3leff <= 5.0e-8)
                    {
                        sw.WriteLine("Warning: Leff = {0} may be too small.".FormatString(Param.BSIM3leff));
                        CircuitWarning.Warning(this, "Warning: Leff = {0} may be too small.".FormatString(Param.BSIM3leff));
                    }

                    if (Param.BSIM3leffCV <= 5.0e-8)
                    {
                        sw.WriteLine("Warning: Leff for CV = {0} may be too small.".FormatString(Param.BSIM3leffCV));
                        CircuitWarning.Warning(this, "Warning: Leff for CV = {0} may be too small.".FormatString(Param.BSIM3leffCV));
                    }

                    if (Param.BSIM3weff <= 1.0e-7)
                    {
                        sw.WriteLine("Warning: Weff = {0} may be too small.".FormatString(Param.BSIM3weff));
                        CircuitWarning.Warning(this, "Warning: Weff = {0} may be too small.".FormatString(Param.BSIM3weff));
                    }

                    if (Param.BSIM3weffCV <= 1.0e-7)
                    {
                        sw.WriteLine("Warning: Weff for CV = {0} may be too small.".FormatString(Param.BSIM3weffCV));
                        CircuitWarning.Warning(this, "Warning: Weff for CV = {0} may be too small.".FormatString(Param.BSIM3weffCV));
                    }

                    /* Check threshold voltage parameters */
                    if (Param.BSIM3nlx < 0.0)
                    {
                        sw.WriteLine("Warning: Nlx = {0} is negative.".FormatString(Param.BSIM3nlx));
                        CircuitWarning.Warning(this, "Warning: Nlx = {0} is negative.".FormatString(Param.BSIM3nlx));
                    }
                    if (ModelParameters.Tox < 1.0e-9)
                    {
                        sw.WriteLine("Warning: Tox = {0} is less than 10A.".FormatString(ModelParameters.Tox));
                        CircuitWarning.Warning(this, "Warning: Tox = {0} is less than 10A.".FormatString(ModelParameters.Tox));
                    }

                    if (Param.BSIM3npeak <= 1.0e15)
                    {
                        sw.WriteLine("Warning: Nch = {0} may be too small.".FormatString(Param.BSIM3npeak));
                        CircuitWarning.Warning(this, "Warning: Nch = {0} may be too small.".FormatString(Param.BSIM3npeak));
                    }
                    else if (Param.BSIM3npeak >= 1.0e21)
                    {
                        sw.WriteLine("Warning: Nch = {0} may be too large.".FormatString(Param.BSIM3npeak));
                        CircuitWarning.Warning(this, "Warning: Nch = {0} may be too large.".FormatString(Param.BSIM3npeak));
                    }

                    if (Param.BSIM3nsub <= 1.0e14)
                    {
                        sw.WriteLine("Warning: Nsub = {0} may be too small.".FormatString(Param.BSIM3nsub));
                        CircuitWarning.Warning(this, "Warning: Nsub = {0} may be too small.".FormatString(Param.BSIM3nsub));
                    }
                    else if (Param.BSIM3nsub >= 1.0e21)
                    {
                        sw.WriteLine("Warning: Nsub = {0} may be too large.".FormatString(Param.BSIM3nsub));
                        CircuitWarning.Warning(this, "Warning: Nsub = {0} may be too large.".FormatString(Param.BSIM3nsub));
                    }

                    if ((Param.BSIM3ngate > 0.0) && (Param.BSIM3ngate <= 1.0e18))
                    {
                        sw.WriteLine("Warning: Ngate = {0} is less than 1.E18cm^-3.".FormatString(Param.BSIM3ngate));
                        CircuitWarning.Warning(this, "Warning: Ngate = {0} is less than 1.E18cm^-3.".FormatString(Param.BSIM3ngate));
                    }

                    if (Param.BSIM3dvt0 < 0.0)
                    {
                        sw.WriteLine("Warning: Dvt0 = {0} is negative.".FormatString(Param.BSIM3dvt0));
                        CircuitWarning.Warning(this, "Warning: Dvt0 = {0} is negative.".FormatString(Param.BSIM3dvt0));
                    }

                    if (Math.Abs(1.0e-6 / (Param.BSIM3w0 + Param.BSIM3weff)) > 10.0)
                    {
                        sw.WriteLine("Warning: (W0 + Weff) may be too small.");
                        CircuitWarning.Warning(this, "Warning: (W0 + Weff) may be too small.");
                    }

                    /* Check subthreshold parameters */
                    if (Param.BSIM3nfactor < 0.0)
                    {
                        sw.WriteLine("Warning: Nfactor = {0} is negative.".FormatString(Param.BSIM3nfactor));
                        CircuitWarning.Warning(this, "Warning: Nfactor = {0} is negative.".FormatString(Param.BSIM3nfactor));
                    }
                    if (Param.BSIM3cdsc < 0.0)
                    {
                        sw.WriteLine("Warning: Cdsc = {0} is negative.".FormatString(Param.BSIM3cdsc));
                        CircuitWarning.Warning(this, "Warning: Cdsc = {0} is negative.".FormatString(Param.BSIM3cdsc));
                    }
                    if (Param.BSIM3cdscd < 0.0)
                    {
                        sw.WriteLine("Warning: Cdscd = {0} is negative.".FormatString(Param.BSIM3cdscd));
                        CircuitWarning.Warning(this, "Warning: Cdscd = {0} is negative.".FormatString(Param.BSIM3cdscd));
                    }
                    /* Check DIBL parameters */
                    if (Param.BSIM3eta0 < 0.0)
                    {
                        sw.WriteLine("Warning: Eta0 = {0} is negative.".FormatString(Param.BSIM3eta0));
                        CircuitWarning.Warning(this, "Warning: Eta0 = {0} is negative.".FormatString(Param.BSIM3eta0));
                    }

                    /* Check Abulk parameters */
                    if (Math.Abs(1.0e-6 / (Param.BSIM3b1 + Param.BSIM3weff)) > 10.0)
                    {
                        sw.WriteLine("Warning: (B1 + Weff) may be too small.");
                        CircuitWarning.Warning(this, "Warning: (B1 + Weff) may be too small.");
                    }


                    /* Check Saturation parameters */
                    if (Param.BSIM3a2 < 0.01)
                    {
                        sw.WriteLine("Warning: A2 = {0} is too small. Set to 0.01.".FormatString(Param.BSIM3a2));
                        CircuitWarning.Warning(this, "Warning: A2 = {0} is too small. Set to 0.01.".FormatString(Param.BSIM3a2));
                        Param.BSIM3a2 = 0.01;
                    }
                    else if (Param.BSIM3a2 > 1.0)
                    {
                        sw.WriteLine("Warning: A2 = {0} is larger than 1. A2 is set to 1 and A1 is set to 0.".FormatString(Param.BSIM3a2));
                        CircuitWarning.Warning(this, "Warning: A2 = {0} is larger than 1. A2 is set to 1 and A1 is set to 0.".FormatString(Param.BSIM3a2));
                        Param.BSIM3a2 = 1.0;
                        Param.BSIM3a1 = 0.0;
                    }

                    if (Param.BSIM3rdsw < 0.0)
                    {
                        sw.WriteLine("Warning: Rdsw = {0} is negative. Set to zero.".FormatString(Param.BSIM3rdsw));
                        CircuitWarning.Warning(this, "Warning: Rdsw = {0} is negative. Set to zero.".FormatString(Param.BSIM3rdsw));
                        Param.BSIM3rdsw = 0.0;
                        Param.BSIM3rds0 = 0.0;
                    }
                    if (Param.BSIM3rds0 < 0.0)
                    {
                        sw.WriteLine("Warning: Rds at current temperature = {0} is negative. Set to zero.".FormatString(Param.BSIM3rds0));
                        CircuitWarning.Warning(this, "Warning: Rds at current temperature = {0} is negative. Set to zero.".FormatString(Param.BSIM3rds0));
                        Param.BSIM3rds0 = 0.0;
                    }

                    if (Param.BSIM3vsattemp < 1.0e3)
                    {
                        sw.WriteLine("Warning: Vsat at current temperature = {0} may be too small.".FormatString(Param.BSIM3vsattemp));
                        CircuitWarning.Warning(this, "Warning: Vsat at current temperature = {0} may be too small.".FormatString(Param.BSIM3vsattemp));
                    }

                    if (Param.BSIM3pdibl1 < 0.0)
                    {
                        sw.WriteLine("Warning: Pdibl1 = {0} is negative.".FormatString(Param.BSIM3pdibl1));
                        CircuitWarning.Warning(this, "Warning: Pdibl1 = {0} is negative.".FormatString(Param.BSIM3pdibl1));
                    }
                    if (Param.BSIM3pdibl2 < 0.0)
                    {
                        sw.WriteLine("Warning: Pdibl2 = {0} is negative.".FormatString(Param.BSIM3pdibl2));
                        CircuitWarning.Warning(this, "Warning: Pdibl2 = {0} is negative.".FormatString(Param.BSIM3pdibl2));
                    }
                    /* Check overlap capacitance parameters */
                    if (ModelParameters.Cgdo < 0.0)
                    {
                        sw.WriteLine("Warning: cgdo = {0} is negative. Set to zero.".FormatString(ModelParameters.Cgdo));
                        CircuitWarning.Warning(this, "Warning: cgdo = {0} is negative. Set to zero.".FormatString(ModelParameters.Cgdo));
                        ModelParameters.Cgdo.RawValue = 0.0;
                    }
                    if (ModelParameters.Cgso < 0.0)
                    {
                        sw.WriteLine("Warning: cgso = {0} is negative. Set to zero.".FormatString(ModelParameters.Cgso));
                        CircuitWarning.Warning(this, "Warning: cgso = {0} is negative. Set to zero.".FormatString(ModelParameters.Cgso));
                        ModelParameters.Cgso.RawValue = 0.0;
                    }
                    if (ModelParameters.Cgbo < 0.0)
                    {
                        sw.WriteLine("Warning: cgbo = {0} is negative. Set to zero.".FormatString(ModelParameters.Cgbo));
                        CircuitWarning.Warning(this, "Warning: cgbo = {0} is negative. Set to zero.".FormatString(ModelParameters.Cgbo));
                        ModelParameters.Cgbo.RawValue = 0.0;
                    }

                }
            }

            return fatal;
        }
	}
}