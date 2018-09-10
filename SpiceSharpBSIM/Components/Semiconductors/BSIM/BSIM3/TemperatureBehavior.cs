using System;
using System.IO;
using SpiceSharp.Behaviors;
using SpiceSharp.Simulations;
namespace SpiceSharp.Components.BSIM3Behaviors
{
	
	/// <summary>
	/// Temperature behavior for a <see cref="BSIM3" />
	/// </summary>
	public class TemperatureBehavior : BaseTemperatureBehavior
	{
		
		/// <summary>
		/// Necessary behaviors and parameters
		/// </summary>
		private ModelTemperatureBehavior _modelTemp;
		private BaseParameters _bp;
		private ModelBaseParameters _mbp;

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
			_modelTemp = provider.GetBehavior<ModelTemperatureBehavior>("model");
			_bp = provider.GetParameterSet<BaseParameters>("entity");
			_mbp = provider.GetParameterSet<ModelBaseParameters>("model");
		}
		
		/// <summary>
		/// Temperature behavior
		/// </summary>
		public override void Temperature(BaseSimulation simulation)
		{
			double tmp, tmp1, tmp2, tmp3, t0 = 0.0, t1, t2, t3, t4, t5, ldrn, wdrn, inv_L, inv_W, inv_LW, nvtm, sourceSatCurrent, drainSatCurrent;

		    var size = new Tuple<double, double>(_bp.Width, _bp.Length);
		    if (!_modelTemp.SizeDependParams.TryGetValue(size, out var pParam))
            {
                pParam = new BSIM3SizeDependParams();
                Param = pParam;
                _modelTemp.SizeDependParams.Add(size, pParam);
				ldrn = _bp.Length;
				wdrn = _bp.Width;
				pParam.Length = ldrn;
				pParam.Width = wdrn;
				t0 = Math.Pow(ldrn, _mbp.Lln);
				t1 = Math.Pow(wdrn, _mbp.Lwn);
				tmp1 = _mbp.Ll / t0 + _mbp.Lw / t1 + _mbp.Lwl / (t0 * t1);
				pParam.BSIM3dl = _mbp.Lint + tmp1;
				tmp2 = _mbp.Llc / t0 + _mbp.Lwc / t1 + _mbp.Lwlc / (t0 * t1);
				pParam.BSIM3dlc = _mbp.Dlc + tmp2;
				t2 = Math.Pow(ldrn, _mbp.Wln);
				t3 = Math.Pow(wdrn, _mbp.Wwn);
				tmp1 = _mbp.Wl / t2 + _mbp.Ww / t3 + _mbp.Wwl / (t2 * t3);
				pParam.BSIM3dw = _mbp.Wint + tmp1;
				tmp2 = _mbp.Wlc / t2 + _mbp.Wwc / t3 + _mbp.Wwlc / (t2 * t3);
				pParam.BSIM3dwc = _mbp.Dwc + tmp2;
				pParam.BSIM3leff = _bp.Length + _mbp.Xl - 2.0 * pParam.BSIM3dl;
				if (pParam.BSIM3leff <= 0.0)
                    throw new CircuitException("BSIM3: mosfet {0}, model {1}: Effective channel length <= 0".FormatString(_modelTemp.Name, Name));
				pParam.BSIM3weff = _bp.Width + _mbp.Xw - 2.0 * pParam.BSIM3dw;
				if (pParam.BSIM3weff <= 0.0)
                    throw new CircuitException("BSIM3: mosfet {0}, model {1}: Effective channel width <= 0".FormatString(_modelTemp.Name, Name));
				pParam.BSIM3leffCV = _bp.Length + _mbp.Xl - 2.0 * pParam.BSIM3dlc;
				if (pParam.BSIM3leffCV <= 0.0)
                    throw new CircuitException("BSIM3: mosfet {0}, model {1}: Effective channel length for C-V <= 0".FormatString(_modelTemp.Name, Name));
				pParam.BSIM3weffCV = _bp.Width + _mbp.Xw - 2.0 * pParam.BSIM3dwc;
				if (pParam.BSIM3weffCV <= 0.0)
                    throw new CircuitException("BSIM3: mosfet {0}, model {1}: Effective channel width for C-V <= 0".FormatString(_modelTemp.Name, Name));
				if (_mbp.BinUnit == 1)
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
				pParam.BSIM3cdsc = _mbp.Cdsc + _mbp.Lcdsc * inv_L + _mbp.Wcdsc * inv_W + _mbp.Pcdsc * inv_LW;
				pParam.BSIM3cdscb = _mbp.Cdscb + _mbp.Lcdscb * inv_L + _mbp.Wcdscb * inv_W + _mbp.Pcdscb * inv_LW;
				pParam.BSIM3cdscd = _mbp.Cdscd + _mbp.Lcdscd * inv_L + _mbp.Wcdscd * inv_W + _mbp.Pcdscd * inv_LW;
				pParam.BSIM3cit = _mbp.Cit + _mbp.Lcit * inv_L + _mbp.Wcit * inv_W + _mbp.Pcit * inv_LW;
				pParam.BSIM3nfactor = _mbp.Nfactor + _mbp.Lnfactor * inv_L + _mbp.Wnfactor * inv_W + _mbp.Pnfactor * inv_LW;
				pParam.BSIM3xj = _mbp.Xj + _mbp.Lxj * inv_L + _mbp.Wxj * inv_W + _mbp.Pxj * inv_LW;
				pParam.BSIM3vsat = _mbp.Vsat + _mbp.Lvsat * inv_L + _mbp.Wvsat * inv_W + _mbp.Pvsat * inv_LW;
				pParam.BSIM3at = _mbp.At + _mbp.Lat * inv_L + _mbp.Wat * inv_W + _mbp.Pat * inv_LW;
				pParam.BSIM3a0 = _mbp.A0 + _mbp.La0 * inv_L + _mbp.Wa0 * inv_W + _mbp.Pa0 * inv_LW;
				pParam.BSIM3ags = _mbp.Ags + _mbp.Lags * inv_L + _mbp.Wags * inv_W + _mbp.Pags * inv_LW;
				pParam.BSIM3a1 = _mbp.A1 + _mbp.La1 * inv_L + _mbp.Wa1 * inv_W + _mbp.Pa1 * inv_LW;
				pParam.BSIM3a2 = _mbp.A2 + _mbp.La2 * inv_L + _mbp.Wa2 * inv_W + _mbp.Pa2 * inv_LW;
				pParam.BSIM3keta = _mbp.Keta + _mbp.Lketa * inv_L + _mbp.Wketa * inv_W + _mbp.Pketa * inv_LW;
				pParam.BSIM3nsub = _mbp.Nsub + _mbp.Lnsub * inv_L + _mbp.Wnsub * inv_W + _mbp.Pnsub * inv_LW;
				pParam.BSIM3npeak = _mbp.Npeak + _mbp.Lnpeak * inv_L + _mbp.Wnpeak * inv_W + _mbp.Pnpeak * inv_LW;
				pParam.BSIM3ngate = _mbp.Ngate + _mbp.Lngate * inv_L + _mbp.Wngate * inv_W + _mbp.Pngate * inv_LW;
				pParam.BSIM3gamma1 = _mbp.Gamma1 + _mbp.Lgamma1 * inv_L + _mbp.Wgamma1 * inv_W + _mbp.Pgamma1 * inv_LW;
				pParam.BSIM3gamma2 = _mbp.Gamma2 + _mbp.Lgamma2 * inv_L + _mbp.Wgamma2 * inv_W + _mbp.Pgamma2 * inv_LW;
				pParam.BSIM3vbx = _mbp.Vbx + _mbp.Lvbx * inv_L + _mbp.Wvbx * inv_W + _mbp.Pvbx * inv_LW;
				pParam.BSIM3vbm = _mbp.Vbm + _mbp.Lvbm * inv_L + _mbp.Wvbm * inv_W + _mbp.Pvbm * inv_LW;
				pParam.BSIM3xt = _mbp.Xt + _mbp.Lxt * inv_L + _mbp.Wxt * inv_W + _mbp.Pxt * inv_LW;
				pParam.BSIM3vfb = _mbp.Vfb + _mbp.Lvfb * inv_L + _mbp.Wvfb * inv_W + _mbp.Pvfb * inv_LW;
				pParam.BSIM3k1 = _mbp.K1 + _mbp.Lk1 * inv_L + _mbp.Wk1 * inv_W + _mbp.Pk1 * inv_LW;
				pParam.BSIM3kt1 = _mbp.Kt1 + _mbp.Lkt1 * inv_L + _mbp.Wkt1 * inv_W + _mbp.Pkt1 * inv_LW;
				pParam.BSIM3kt1l = _mbp.Kt1l + _mbp.Lkt1l * inv_L + _mbp.Wkt1l * inv_W + _mbp.Pkt1l * inv_LW;
				pParam.BSIM3k2 = _mbp.K2 + _mbp.Lk2 * inv_L + _mbp.Wk2 * inv_W + _mbp.Pk2 * inv_LW;
				pParam.BSIM3kt2 = _mbp.Kt2 + _mbp.Lkt2 * inv_L + _mbp.Wkt2 * inv_W + _mbp.Pkt2 * inv_LW;
				pParam.BSIM3k3 = _mbp.K3 + _mbp.Lk3 * inv_L + _mbp.Wk3 * inv_W + _mbp.Pk3 * inv_LW;
				pParam.BSIM3k3b = _mbp.K3b + _mbp.Lk3b * inv_L + _mbp.Wk3b * inv_W + _mbp.Pk3b * inv_LW;
				pParam.BSIM3w0 = _mbp.W0 + _mbp.Lw0 * inv_L + _mbp.Ww0 * inv_W + _mbp.Pw0 * inv_LW;
				pParam.BSIM3nlx = _mbp.Nlx + _mbp.Lnlx * inv_L + _mbp.Wnlx * inv_W + _mbp.Pnlx * inv_LW;
				pParam.BSIM3dvt0 = _mbp.Dvt0 + _mbp.Ldvt0 * inv_L + _mbp.Wdvt0 * inv_W + _mbp.Pdvt0 * inv_LW;
				pParam.BSIM3dvt1 = _mbp.Dvt1 + _mbp.Ldvt1 * inv_L + _mbp.Wdvt1 * inv_W + _mbp.Pdvt1 * inv_LW;
				pParam.BSIM3dvt2 = _mbp.Dvt2 + _mbp.Ldvt2 * inv_L + _mbp.Wdvt2 * inv_W + _mbp.Pdvt2 * inv_LW;
				pParam.BSIM3dvt0w = _mbp.Dvt0w + _mbp.Ldvt0w * inv_L + _mbp.Wdvt0w * inv_W + _mbp.Pdvt0w * inv_LW;
				pParam.BSIM3dvt1w = _mbp.Dvt1w + _mbp.Ldvt1w * inv_L + _mbp.Wdvt1w * inv_W + _mbp.Pdvt1w * inv_LW;
				pParam.BSIM3dvt2w = _mbp.Dvt2w + _mbp.Ldvt2w * inv_L + _mbp.Wdvt2w * inv_W + _mbp.Pdvt2w * inv_LW;
				pParam.BSIM3drout = _mbp.Drout + _mbp.Ldrout * inv_L + _mbp.Wdrout * inv_W + _mbp.Pdrout * inv_LW;
				pParam.BSIM3dsub = _mbp.Dsub + _mbp.Ldsub * inv_L + _mbp.Wdsub * inv_W + _mbp.Pdsub * inv_LW;
				pParam.BSIM3vth0 = _mbp.Vth0 + _mbp.Lvth0 * inv_L + _mbp.Wvth0 * inv_W + _mbp.Pvth0 * inv_LW;
				pParam.BSIM3ua = _mbp.Ua + _mbp.Lua * inv_L + _mbp.Wua * inv_W + _mbp.Pua * inv_LW;
				pParam.BSIM3ua1 = _mbp.Ua1 + _mbp.Lua1 * inv_L + _mbp.Wua1 * inv_W + _mbp.Pua1 * inv_LW;
				pParam.BSIM3ub = _mbp.Ub + _mbp.Lub * inv_L + _mbp.Wub * inv_W + _mbp.Pub * inv_LW;
				pParam.BSIM3ub1 = _mbp.Ub1 + _mbp.Lub1 * inv_L + _mbp.Wub1 * inv_W + _mbp.Pub1 * inv_LW;
				pParam.BSIM3uc = _mbp.Uc + _mbp.Luc * inv_L + _mbp.Wuc * inv_W + _mbp.Puc * inv_LW;
				pParam.BSIM3uc1 = _mbp.Uc1 + _mbp.Luc1 * inv_L + _mbp.Wuc1 * inv_W + _mbp.Puc1 * inv_LW;
				pParam.BSIM3u0 = _mbp.U0 + _mbp.Lu0 * inv_L + _mbp.Wu0 * inv_W + _mbp.Pu0 * inv_LW;
				pParam.BSIM3ute = _mbp.Ute + _mbp.Lute * inv_L + _mbp.Wute * inv_W + _mbp.Pute * inv_LW;
				pParam.BSIM3voff = _mbp.Voff + _mbp.Lvoff * inv_L + _mbp.Wvoff * inv_W + _mbp.Pvoff * inv_LW;
				pParam.BSIM3delta = _mbp.Delta + _mbp.Ldelta * inv_L + _mbp.Wdelta * inv_W + _mbp.Pdelta * inv_LW;
				pParam.BSIM3rdsw = _mbp.Rdsw + _mbp.Lrdsw * inv_L + _mbp.Wrdsw * inv_W + _mbp.Prdsw * inv_LW;
				pParam.BSIM3prwg = _mbp.Prwg + _mbp.Lprwg * inv_L + _mbp.Wprwg * inv_W + _mbp.Pprwg * inv_LW;
				pParam.BSIM3prwb = _mbp.Prwb + _mbp.Lprwb * inv_L + _mbp.Wprwb * inv_W + _mbp.Pprwb * inv_LW;
				pParam.BSIM3prt = _mbp.Prt + _mbp.Lprt * inv_L + _mbp.Wprt * inv_W + _mbp.Pprt * inv_LW;
				pParam.BSIM3eta0 = _mbp.Eta0 + _mbp.Leta0 * inv_L + _mbp.Weta0 * inv_W + _mbp.Peta0 * inv_LW;
				pParam.BSIM3etab = _mbp.Etab + _mbp.Letab * inv_L + _mbp.Wetab * inv_W + _mbp.Petab * inv_LW;
				pParam.BSIM3pclm = _mbp.Pclm + _mbp.Lpclm * inv_L + _mbp.Wpclm * inv_W + _mbp.Ppclm * inv_LW;
				pParam.BSIM3pdibl1 = _mbp.Pdibl1 + _mbp.Lpdibl1 * inv_L + _mbp.Wpdibl1 * inv_W + _mbp.Ppdibl1 * inv_LW;
				pParam.BSIM3pdibl2 = _mbp.Pdibl2 + _mbp.Lpdibl2 * inv_L + _mbp.Wpdibl2 * inv_W + _mbp.Ppdibl2 * inv_LW;
				pParam.BSIM3pdiblb = _mbp.Pdiblb + _mbp.Lpdiblb * inv_L + _mbp.Wpdiblb * inv_W + _mbp.Ppdiblb * inv_LW;
				pParam.BSIM3pscbe1 = _mbp.Pscbe1 + _mbp.Lpscbe1 * inv_L + _mbp.Wpscbe1 * inv_W + _mbp.Ppscbe1 * inv_LW;
				pParam.BSIM3pscbe2 = _mbp.Pscbe2 + _mbp.Lpscbe2 * inv_L + _mbp.Wpscbe2 * inv_W + _mbp.Ppscbe2 * inv_LW;
				pParam.BSIM3pvag = _mbp.Pvag + _mbp.Lpvag * inv_L + _mbp.Wpvag * inv_W + _mbp.Ppvag * inv_LW;
				pParam.BSIM3wr = _mbp.Wr + _mbp.Lwr * inv_L + _mbp.Wwr * inv_W + _mbp.Pwr * inv_LW;
				pParam.BSIM3dwg = _mbp.Dwg + _mbp.Ldwg * inv_L + _mbp.Wdwg * inv_W + _mbp.Pdwg * inv_LW;
				pParam.BSIM3dwb = _mbp.Dwb + _mbp.Ldwb * inv_L + _mbp.Wdwb * inv_W + _mbp.Pdwb * inv_LW;
				pParam.BSIM3b0 = _mbp.B0 + _mbp.Lb0 * inv_L + _mbp.Wb0 * inv_W + _mbp.Pb0 * inv_LW;
				pParam.BSIM3b1 = _mbp.B1 + _mbp.Lb1 * inv_L + _mbp.Wb1 * inv_W + _mbp.Pb1 * inv_LW;
				pParam.BSIM3alpha0 = _mbp.Alpha0 + _mbp.Lalpha0 * inv_L + _mbp.Walpha0 * inv_W + _mbp.Palpha0 * inv_LW;
				pParam.BSIM3alpha1 = _mbp.Alpha1 + _mbp.Lalpha1 * inv_L + _mbp.Walpha1 * inv_W + _mbp.Palpha1 * inv_LW;
				pParam.BSIM3beta0 = _mbp.Beta0 + _mbp.Lbeta0 * inv_L + _mbp.Wbeta0 * inv_W + _mbp.Pbeta0 * inv_LW;
				pParam.BSIM3elm = _mbp.Elm + _mbp.Lelm * inv_L + _mbp.Welm * inv_W + _mbp.Pelm * inv_LW;
				pParam.BSIM3cgsl = _mbp.Cgsl + _mbp.Lcgsl * inv_L + _mbp.Wcgsl * inv_W + _mbp.Pcgsl * inv_LW;
				pParam.BSIM3cgdl = _mbp.Cgdl + _mbp.Lcgdl * inv_L + _mbp.Wcgdl * inv_W + _mbp.Pcgdl * inv_LW;
				pParam.BSIM3ckappa = _mbp.Ckappa + _mbp.Lckappa * inv_L + _mbp.Wckappa * inv_W + _mbp.Pckappa * inv_LW;
				pParam.BSIM3cf = _mbp.Cf + _mbp.Lcf * inv_L + _mbp.Wcf * inv_W + _mbp.Pcf * inv_LW;
				pParam.BSIM3clc = _mbp.Clc + _mbp.Lclc * inv_L + _mbp.Wclc * inv_W + _mbp.Pclc * inv_LW;
				pParam.BSIM3cle = _mbp.Cle + _mbp.Lcle * inv_L + _mbp.Wcle * inv_W + _mbp.Pcle * inv_LW;
				pParam.BSIM3vfbcv = _mbp.Vfbcv + _mbp.Lvfbcv * inv_L + _mbp.Wvfbcv * inv_W + _mbp.Pvfbcv * inv_LW;
				pParam.BSIM3acde = _mbp.Acde + _mbp.Lacde * inv_L + _mbp.Wacde * inv_W + _mbp.Pacde * inv_LW;
				pParam.BSIM3moin = _mbp.Moin + _mbp.Lmoin * inv_L + _mbp.Wmoin * inv_W + _mbp.Pmoin * inv_LW;
				pParam.BSIM3noff = _mbp.Noff + _mbp.Lnoff * inv_L + _mbp.Wnoff * inv_W + _mbp.Pnoff * inv_LW;
				pParam.BSIM3voffcv = _mbp.Voffcv + _mbp.Lvoffcv * inv_L + _mbp.Wvoffcv * inv_W + _mbp.Pvoffcv * inv_LW;
				pParam.BSIM3abulkCVfactor = 1.0 + Math.Pow(pParam.BSIM3clc / pParam.BSIM3leffCV, pParam.BSIM3cle);
				t0 = _modelTemp.TRatio - 1.0;
				pParam.BSIM3ua = pParam.BSIM3ua + pParam.BSIM3ua1 * t0;
				pParam.BSIM3ub = pParam.BSIM3ub + pParam.BSIM3ub1 * t0;
				pParam.BSIM3uc = pParam.BSIM3uc + pParam.BSIM3uc1 * t0;
				if (pParam.BSIM3u0 > 1.0)
					pParam.BSIM3u0 = pParam.BSIM3u0 / 1.0e4;
				pParam.BSIM3u0temp = pParam.BSIM3u0 * Math.Pow(_modelTemp.TRatio, pParam.BSIM3ute);
				pParam.BSIM3vsattemp = pParam.BSIM3vsat - pParam.BSIM3at * t0;
				pParam.BSIM3rds0 = (pParam.BSIM3rdsw + pParam.BSIM3prt * t0) / Math.Pow(pParam.BSIM3weff * 1E6, pParam.BSIM3wr);
				if (Check())
                    throw new CircuitException("Fatal error(s) detected during BSIM3V3.3 parameter checking for {0} in model {1}".FormatString(Name, _modelTemp.Name));
				pParam.BSIM3cgdo = (_mbp.Cgdo + pParam.BSIM3cf) * pParam.BSIM3weffCV;
				pParam.BSIM3cgso = (_mbp.Cgso + pParam.BSIM3cf) * pParam.BSIM3weffCV;
				pParam.BSIM3cgbo = _mbp.Cgbo * pParam.BSIM3leffCV;
				t0 = pParam.BSIM3leffCV * pParam.BSIM3leffCV;
				pParam.BSIM3tconst = pParam.BSIM3u0temp * pParam.BSIM3elm / (_mbp.Cox * pParam.BSIM3weffCV * pParam.BSIM3leffCV * t0);
				if (!_mbp.Npeak.Given && _mbp.Gamma1.Given)
				{
					t0 = pParam.BSIM3gamma1 * _mbp.Cox;
					pParam.BSIM3npeak = 3.021E22 * t0 * t0;
				}
				pParam.BSIM3phi = 2.0 * _modelTemp.Vtm0 * Math.Log(pParam.BSIM3npeak / _modelTemp.Ni);
				pParam.BSIM3sqrtPhi = Math.Sqrt(pParam.BSIM3phi);
				pParam.BSIM3phis3 = pParam.BSIM3sqrtPhi * pParam.BSIM3phi;
				pParam.BSIM3Xdep0 = Math.Sqrt(2.0 * 1.03594e-10 / (1.60219e-19 * pParam.BSIM3npeak * 1.0e6)) * pParam.BSIM3sqrtPhi;
				pParam.BSIM3sqrtXdep0 = Math.Sqrt(pParam.BSIM3Xdep0);
				pParam.BSIM3litl = Math.Sqrt(3.0 * pParam.BSIM3xj * _mbp.Tox);
				pParam.BSIM3vbi = _modelTemp.Vtm0 * Math.Log(1.0e20 * pParam.BSIM3npeak / (_modelTemp.Ni * _modelTemp.Ni));
				pParam.BSIM3cdep0 = Math.Sqrt(1.60219e-19 * 1.03594e-10 * pParam.BSIM3npeak * 1.0e6 / 2.0 / pParam.BSIM3phi);
				pParam.BSIM3ldeb = Math.Sqrt(1.03594e-10 * _modelTemp.Vtm0 / (1.60219e-19 * pParam.BSIM3npeak * 1.0e6)) / 3.0;
				pParam.BSIM3acde *= Math.Pow(pParam.BSIM3npeak / 2.0e16, -0.25);
                if (_mbp.K1.Given || _mbp.K2.Given)
                {
                    if (!_mbp.K1.Given)
                    {
                        CircuitWarning.Warning(this, "Warning: k1 should be specified with k2.");
                        pParam.BSIM3k1 = 0.53;
                    }

                    if (!_mbp.K2.Given)
                    {
                        CircuitWarning.Warning(this, "Warning: k2 should be specified with k1.");
                        pParam.BSIM3k2 = -0.0186;
                    }

                    if (_mbp.Nsub.Given)
                        CircuitWarning.Warning(this, "Warning: nsub is ignored because k1 or k2 is given.");
                    if (_mbp.Xt.Given)
                        CircuitWarning.Warning(this, "Warning: xt is ignored because k1 or k2 is given.");
                    if (_mbp.Vbx.Given)
                        CircuitWarning.Warning(this, "Warning: vbx is ignored because k1 or k2 is given.");
                    if (_mbp.Gamma1.Given)
                        CircuitWarning.Warning(this, "Warning: gamma1 is ignored because k1 or k2 is given.");
                    if (_mbp.Gamma2.Given)
                        CircuitWarning.Warning(this, "Warning: gamma2 is ignored because k1 or k2 is given.");
                }
                else
				{
					if (!_mbp.Vbx.Given)
						pParam.BSIM3vbx = pParam.BSIM3phi - 7.7348e-4 * pParam.BSIM3npeak * pParam.BSIM3xt * pParam.BSIM3xt;
					if (pParam.BSIM3vbx > 0.0)
						pParam.BSIM3vbx = -pParam.BSIM3vbx;
					if (pParam.BSIM3vbm > 0.0)
						pParam.BSIM3vbm = -pParam.BSIM3vbm;
					if (!_mbp.Gamma1.Given)
						pParam.BSIM3gamma1 = 5.753e-12 * Math.Sqrt(pParam.BSIM3npeak) / _mbp.Cox;
					if (!_mbp.Gamma2.Given)
						pParam.BSIM3gamma2 = 5.753e-12 * Math.Sqrt(pParam.BSIM3nsub) / _mbp.Cox;
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
				if (!_mbp.Vfb.Given)
				{
					if (_mbp.Vth0.Given)
					{
						pParam.BSIM3vfb = _mbp.B3Type * pParam.BSIM3vth0 - pParam.BSIM3phi - pParam.BSIM3k1 * pParam.BSIM3sqrtPhi;
					}
					else
					{
						pParam.BSIM3vfb = -1.0;
					}
				}
				if (!_mbp.Vth0.Given)
				{
					pParam.BSIM3vth0 = _mbp.B3Type * (pParam.BSIM3vfb + pParam.BSIM3phi + pParam.BSIM3k1 * pParam.BSIM3sqrtPhi);
				}
				pParam.BSIM3k1ox = pParam.BSIM3k1 * _mbp.Tox / _mbp.Toxm;
				pParam.BSIM3k2ox = pParam.BSIM3k2 * _mbp.Tox / _mbp.Toxm;
				t1 = Math.Sqrt(1.03594e-10 / 3.453133e-11 * _mbp.Tox * pParam.BSIM3Xdep0);
				t0 = Math.Exp(-0.5 * pParam.BSIM3dsub * pParam.BSIM3leff / t1);
				pParam.BSIM3theta0vb0 = t0 + 2.0 * t0 * t0;
				t0 = Math.Exp(-0.5 * pParam.BSIM3drout * pParam.BSIM3leff / t1);
				t2 = t0 + 2.0 * t0 * t0;
				pParam.BSIM3thetaRout = pParam.BSIM3pdibl1 * t2 + pParam.BSIM3pdibl2;
				tmp = Math.Sqrt(pParam.BSIM3Xdep0);
				tmp1 = pParam.BSIM3vbi - pParam.BSIM3phi;
				tmp2 = _modelTemp.Factor1 * tmp;
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
				t4 = _mbp.Tox * pParam.BSIM3phi / (pParam.BSIM3weff + pParam.BSIM3w0);
				t0 = Math.Sqrt(1.0 + pParam.BSIM3nlx / pParam.BSIM3leff);
				t5 = pParam.BSIM3k1ox * (t0 - 1.0) * pParam.BSIM3sqrtPhi + (pParam.BSIM3kt1 + pParam.BSIM3kt1l / pParam.BSIM3leff) * (_modelTemp.TRatio - 1.0);
				tmp3 = _mbp.B3Type * pParam.BSIM3vth0 - t2 - t3 + pParam.BSIM3k3 * t4 + t5;
				pParam.BSIM3vfbzb = tmp3 - pParam.BSIM3phi - pParam.BSIM3k1 * pParam.BSIM3sqrtPhi;
			}

		    Param = pParam;
			Vth0 = Param.BSIM3vth0 + _bp.Delvto;
			Vfb = Param.BSIM3vfb + _mbp.B3Type * _bp.Delvto;
			Vfbzb = Param.BSIM3vfbzb + _mbp.B3Type * _bp.Delvto;
			U0temp = Param.BSIM3u0temp * _bp.Mulu0;
			Tconst = U0temp * Param.BSIM3elm / (_mbp.Cox * Param.BSIM3weffCV * Param.BSIM3leffCV * t0);

		    var drainResistance = _mbp.SheetResistance * _bp.DrainSquares;
		    var sourceResistance = _mbp.SheetResistance * _bp.SourceSquares;

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
			nvtm = _modelTemp.Vtm * _mbp.JctEmissionCoeff;
			if (_bp.SourceArea <= 0.0 && _bp.SourcePerimeter <= 0.0)
			{
				sourceSatCurrent = 1.0e-14;
			}
			else
			{
				sourceSatCurrent = _bp.SourceArea * _modelTemp.JctTempSatCurDensity + _bp.SourcePerimeter * _modelTemp.JctSidewallTempSatCurDensity;
			}
			if (sourceSatCurrent > 0.0 && _mbp.Ijth > 0.0)
			{
				Vjsm = nvtm * Math.Log(_mbp.Ijth / sourceSatCurrent + 1.0);
				IsEvjsm = sourceSatCurrent * Math.Exp(Vjsm / nvtm);
			}
			if (_bp.DrainArea <= 0.0 && _bp.DrainPerimeter <= 0.0)
			{
				drainSatCurrent = 1.0e-14;
			}
			else
			{
				drainSatCurrent = _bp.DrainArea * _modelTemp.JctTempSatCurDensity + _bp.DrainPerimeter * _modelTemp.JctSidewallTempSatCurDensity;
			}
			if (drainSatCurrent > 0.0 && _mbp.Ijth > 0.0)
			{
				Vjdm = nvtm * Math.Log(_mbp.Ijth / drainSatCurrent + 1.0);
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
            using (var sw = new StreamWriter(_mbp.CheckPath))
            {
                sw.WriteLine("BSIM3v3.3.0 Parameter Checking.");
                if (_mbp.Version != "3.3.0" && _mbp.Version != "3.30" && _mbp.Version != "3.3")
                {
                    sw.WriteLine("Warning: This model is BSIM3v3.3.0; you specified a wrong version number.");
                    CircuitWarning.Warning(this, "Warning: This model is BSIM3v3.3.0; you specified a wrong version number.");
                }
                sw.WriteLine("Model = {0}".FormatString(_modelTemp.Name));

                if (Param.BSIM3nlx < -Param.BSIM3leff)
                {
                    sw.WriteLine("Fatal: Nlx = {0} is less than -Leff.".FormatString(Param.BSIM3nlx));
                    CircuitWarning.Warning(this, "Fatal: Nlx = {0} is less than -Leff.".FormatString(Param.BSIM3nlx));
                    fatal = true;
                }

                if (_mbp.Tox <= 0.0)
                {
                    sw.WriteLine("Fatal: Tox = {0} is not positive.".FormatString(_mbp.Tox));
                    CircuitWarning.Warning(this, "Fatal: Tox = {0} is not positive.".FormatString(_mbp.Tox));
                    fatal = true;
                }

                if (_mbp.Toxm <= 0.0)
                {
                    sw.WriteLine("Fatal: Toxm = {0} is not positive.".FormatString(_mbp.Toxm));
                    CircuitWarning.Warning(this, "Fatal: Toxm = {0} is not positive.".FormatString(_mbp.Toxm));
                    fatal = true;
                }

                if (_mbp.Lintnoi > Param.BSIM3leff / 2)
                {
                    sw.WriteLine("Fatal: Lintnoi = {0} is too large - Leff for noise is negative.".FormatString(_mbp.Lintnoi));
                    CircuitWarning.Warning(this,
                        "Fatal: Lintnoi = {0} is too large - Leff for noise is negative.".FormatString(_mbp.Lintnoi));
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

                if (_mbp.Ijth < 0.0)
                {
                    sw.WriteLine("Fatal: Ijth = {0} cannot be negative.".FormatString(_mbp.Ijth));
                    CircuitWarning.Warning(this, "Fatal: Ijth = {0} cannot be negative.".FormatString(_mbp.Ijth));
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

                if (_mbp.CapMod == 3)
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

                if (_mbp.ParamChk == 1)
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
                    if (_mbp.Tox < 1.0e-9)
                    {
                        sw.WriteLine("Warning: Tox = {0} is less than 10A.".FormatString(_mbp.Tox));
                        CircuitWarning.Warning(this, "Warning: Tox = {0} is less than 10A.".FormatString(_mbp.Tox));
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
                    if (_mbp.Cgdo < 0.0)
                    {
                        sw.WriteLine("Warning: cgdo = {0} is negative. Set to zero.".FormatString(_mbp.Cgdo));
                        CircuitWarning.Warning(this, "Warning: cgdo = {0} is negative. Set to zero.".FormatString(_mbp.Cgdo));
                        _mbp.Cgdo.RawValue = 0.0;
                    }
                    if (_mbp.Cgso < 0.0)
                    {
                        sw.WriteLine("Warning: cgso = {0} is negative. Set to zero.".FormatString(_mbp.Cgso));
                        CircuitWarning.Warning(this, "Warning: cgso = {0} is negative. Set to zero.".FormatString(_mbp.Cgso));
                        _mbp.Cgso.RawValue = 0.0;
                    }
                    if (_mbp.Cgbo < 0.0)
                    {
                        sw.WriteLine("Warning: cgbo = {0} is negative. Set to zero.".FormatString(_mbp.Cgbo));
                        CircuitWarning.Warning(this, "Warning: cgbo = {0} is negative. Set to zero.".FormatString(_mbp.Cgbo));
                        _mbp.Cgbo.RawValue = 0.0;
                    }

                }
            }

            return fatal;
        }
	}
}