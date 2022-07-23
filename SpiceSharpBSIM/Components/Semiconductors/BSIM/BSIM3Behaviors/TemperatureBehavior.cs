using System;
using System.Collections.Generic;
using System.IO;
using SpiceSharp;
using SpiceSharp.Attributes;
using SpiceSharp.Behaviors;
using SpiceSharp.Components;
using SpiceSharp.ParameterSets;
using SpiceSharp.Simulations;

namespace SpiceSharpBSIM.Components.Semiconductors.BSIM.BSIM3Behaviors
{
    /// <summary>
    /// Temperature behavior for a <see cref="BSIM3" />
    /// </summary>
    [BehaviorFor(typeof(BSIM3)), AddBehaviorIfNo(typeof(ITemperatureBehavior))]
    public class TemperatureBehavior : Behavior, ITemperatureBehavior, IParameterized<BaseParameters>
    {
        public const double Kb = 1.3806226e-23;
        public const double KboQ = 8.617087e-5;
        public const double EPSSI = 1.03594e-10;
        public const double EPSOX = 3.453133e-11;
        public const double Charge_q = 1.60219e-19;
        public const double MAX_EXP = 5.834617425e14;
        public const double MIN_EXP = 1.713908431e-15;
        public const double EXP_THRESHOLD = 34.0;

        /// <inheritdoc />
        public BaseParameters Parameters { get; }

        protected ModelTemperatureBehavior ModelTemperature { get; }
        protected ModelParameters ModelParameters { get; }

        /// <summary>
        /// Size dependent parameters
        /// </summary>
	    protected SizeDependParams Param { get; private set; }

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
        public TemperatureBehavior(ComponentBindingContext context)
            : base(context)
        {
            Parameters = context.GetParameterSet<BaseParameters>();
            ModelParameters = context.ModelBehaviors.GetParameterSet<ModelParameters>();
            ModelTemperature = context.ModelBehaviors.GetValue<ModelTemperatureBehavior>();
        }

        /// <summary>
        /// Temperature behavior
        /// </summary>
        void ITemperatureBehavior.Temperature()
        {
            double tmp, tmp1, tmp2, tmp3, ni, T0, T1, T2, T3, T4, T5, Ldrn, Wdrn;
            double TRatio, Inv_L, Inv_W, Inv_LW, Vtm0;
            double Nvtm, SourceSatCurrent, DrainSatCurrent;

            TRatio = ModelTemperature.TRatio;
            Vtm0 = ModelTemperature.Vtm0;
            ni = ModelTemperature.Ni;
            T0 = ModelTemperature.T0; // Why is temporary variable even needed?

            // Parameter defaulting in setup
            if (!Parameters.DrainSquares.Given)
            {
                if (ModelParameters.AcmMod.Value == 0)
                    Parameters.DrainSquares = new GivenParameter<double>(1.0, false);
                else
                    Parameters.DrainSquares = new GivenParameter<double>(0.0, false);
            }
            if (!Parameters.SourceSquares.Given)
            {
                if (ModelParameters.AcmMod.Value == 0)
                    Parameters.SourceSquares = new GivenParameter<double>(1.0, false);
                else
                    Parameters.SourceSquares = new GivenParameter<double>(0.0, false);
            }
            if (!Parameters.NqsMod.Given)
                Parameters.NqsMod = new GivenParameter<int>(ModelParameters.NqsMod, false);
            else if ((Parameters.NqsMod.Value != 0) && (Parameters.NqsMod.Value != 1))
            {
                Parameters.NqsMod = new GivenParameter<int>(ModelParameters.NqsMod, false);
                SpiceSharpWarning.Warning(this, string.Format("Warning: nqsMod has been set to its global value {0}.",
                    ModelParameters.NqsMod));
            }
            if (!Parameters.AcnqsMod.Given)
                Parameters.AcnqsMod = new GivenParameter<int>(ModelParameters.AcnqsMod, false);
            else if ((Parameters.AcnqsMod.Value != 0) && (Parameters.AcnqsMod.Value != 1))
            {
                Parameters.AcnqsMod = new GivenParameter<int>(Parameters.AcnqsMod, false);
                SpiceSharpWarning.Warning(this, string.Format("Warning: acnqsMod has been set to its global value %d.\n",
                    ModelParameters.AcnqsMod));
            }

            var key = Tuple.Create(Parameters.Width.Value, Parameters.Length.Value);
            if (ModelTemperature.SizeDependParams.TryGetValue(key, out var param))
                Param = param;
            else
            {
                Param = new SizeDependParams();
                ModelTemperature.SizeDependParams.Add(key, Param);

                Ldrn = Parameters.Length;
                Wdrn = Parameters.Width;
                Param.Length = Ldrn;
                Param.Width = Wdrn;

                T0 = Math.Pow(Ldrn, ModelParameters.Lln);
                T1 = Math.Pow(Wdrn, ModelParameters.Lwn);
                tmp1 = ModelParameters.Ll / T0 + ModelParameters.Lw / T1
                               + ModelParameters.Lwl / (T0 * T1);
                Param.BSIM3dl = ModelParameters.Lint + tmp1;
                tmp2 = ModelParameters.Llc / T0 + ModelParameters.Lwc / T1
                     + ModelParameters.Lwlc / (T0 * T1);
                Param.BSIM3dlc = ModelParameters.Dlc + tmp2;

                T2 = Math.Pow(Ldrn, ModelParameters.Wln);
                T3 = Math.Pow(Wdrn, ModelParameters.Wwn);
                tmp1 = ModelParameters.Wl / T2 + ModelParameters.Ww / T3
                               + ModelParameters.Wwl / (T2 * T3);
                Param.BSIM3dw = ModelParameters.Wint + tmp1;
                tmp2 = ModelParameters.Wlc / T2 + ModelParameters.Wwc / T3
                     + ModelParameters.Wwlc / (T2 * T3);
                Param.BSIM3dwc = ModelParameters.Dwc + tmp2;

                Param.BSIM3leff = Parameters.Length + ModelParameters.Xl - 2.0 * Param.BSIM3dl;
                if (Param.BSIM3leff <= 0.0)
                {
                    throw new SpiceSharpException("BSIM3: mosfet {0}, model {1}: Effective channel length <= 0".FormatString(Name, ModelTemperature.Name));
                }

                Param.BSIM3weff = Parameters.Width + ModelParameters.Xw - 2.0 * Param.BSIM3dw;
                if (Param.BSIM3weff <= 0.0)
                {
                    throw new SpiceSharpException("BSIM3: mosfet {0}, model {1}: Effective channel width <= 0".FormatString(Name, ModelTemperature.Name));
                }

                Param.BSIM3leffCV = Parameters.Length + ModelParameters.Xl - 2.0 * Param.BSIM3dlc;
                if (Param.BSIM3leffCV <= 0.0)
                {
                    throw new SpiceSharpException("BSIM3: mosfet {0}, model {1}: Effective channel length for C-V <= 0".FormatString(Name, ModelTemperature.Name));
                }

                Param.BSIM3weffCV = Parameters.Width + ModelParameters.Xw - 2.0 * Param.BSIM3dwc;
                if (Param.BSIM3weffCV <= 0.0)
                {
                    throw new SpiceSharpException("BSIM3: mosfet {0}, model {1}: Effective channel width for C-V <= 0".FormatString(Name, ModelTemperature.Name));
                }


                if (ModelParameters.BinUnit == 1)
                {
                    Inv_L = 1.0e-6 / Param.BSIM3leff;
                    Inv_W = 1.0e-6 / Param.BSIM3weff;
                    Inv_LW = 1.0e-12 / (Param.BSIM3leff
                           * Param.BSIM3weff);
                }
                else
                {
                    Inv_L = 1.0 / Param.BSIM3leff;
                    Inv_W = 1.0 / Param.BSIM3weff;
                    Inv_LW = 1.0 / (Param.BSIM3leff
                           * Param.BSIM3weff);
                }
                Param.BSIM3cdsc = ModelParameters.Cdsc
                                  + ModelParameters.Lcdsc * Inv_L
                                  + ModelParameters.Wcdsc * Inv_W
                                  + ModelParameters.Pcdsc * Inv_LW;
                Param.BSIM3cdscb = ModelParameters.Cdscb
                                   + ModelParameters.Lcdscb * Inv_L
                                   + ModelParameters.Wcdscb * Inv_W
                                   + ModelParameters.Pcdscb * Inv_LW;

                Param.BSIM3cdscd = ModelParameters.Cdscd
                                   + ModelParameters.Lcdscd * Inv_L
                                   + ModelParameters.Wcdscd * Inv_W
                                   + ModelParameters.Pcdscd * Inv_LW;

                Param.BSIM3cit = ModelParameters.Cit
                                 + ModelParameters.Lcit * Inv_L
                                 + ModelParameters.Wcit * Inv_W
                                 + ModelParameters.Pcit * Inv_LW;
                Param.BSIM3nfactor = ModelParameters.Nfactor
                                     + ModelParameters.Lnfactor * Inv_L
                                     + ModelParameters.Wnfactor * Inv_W
                                     + ModelParameters.Pnfactor * Inv_LW;
                Param.BSIM3xj = ModelParameters.Xj
                                + ModelParameters.Lxj * Inv_L
                                + ModelParameters.Wxj * Inv_W
                                + ModelParameters.Pxj * Inv_LW;
                Param.BSIM3vsat = ModelParameters.Vsat
                                  + ModelParameters.Lvsat * Inv_L
                                  + ModelParameters.Wvsat * Inv_W
                                  + ModelParameters.Pvsat * Inv_LW;
                Param.BSIM3at = ModelParameters.At
                                + ModelParameters.Lat * Inv_L
                                + ModelParameters.Wat * Inv_W
                                + ModelParameters.Pat * Inv_LW;
                Param.BSIM3a0 = ModelParameters.A0
                                + ModelParameters.La0 * Inv_L
                                + ModelParameters.Wa0 * Inv_W
                                + ModelParameters.Pa0 * Inv_LW;

                Param.BSIM3ags = ModelParameters.Ags
                                + ModelParameters.Lags * Inv_L
                                + ModelParameters.Wags * Inv_W
                                + ModelParameters.Pags * Inv_LW;

                Param.BSIM3a1 = ModelParameters.A1
                                + ModelParameters.La1 * Inv_L
                                + ModelParameters.Wa1 * Inv_W
                                + ModelParameters.Pa1 * Inv_LW;
                Param.BSIM3a2 = ModelParameters.A2
                                + ModelParameters.La2 * Inv_L
                                + ModelParameters.Wa2 * Inv_W
                                + ModelParameters.Pa2 * Inv_LW;
                Param.BSIM3keta = ModelParameters.Keta
                                  + ModelParameters.Lketa * Inv_L
                                  + ModelParameters.Wketa * Inv_W
                                  + ModelParameters.Pketa * Inv_LW;
                Param.BSIM3nsub = ModelParameters.Nsub
                                  + ModelParameters.Lnsub * Inv_L
                                  + ModelParameters.Wnsub * Inv_W
                                  + ModelParameters.Pnsub * Inv_LW;
                Param.BSIM3npeak = ModelParameters.Npeak
                                   + ModelParameters.Lnpeak * Inv_L
                                   + ModelParameters.Wnpeak * Inv_W
                                   + ModelParameters.Pnpeak * Inv_LW;
                Param.BSIM3ngate = ModelParameters.Ngate
                                   + ModelParameters.Lngate * Inv_L
                                   + ModelParameters.Wngate * Inv_W
                                   + ModelParameters.Pngate * Inv_LW;
                Param.BSIM3gamma1 = ModelParameters.Gamma1
                                    + ModelParameters.Lgamma1 * Inv_L
                                    + ModelParameters.Wgamma1 * Inv_W
                                    + ModelParameters.Pgamma1 * Inv_LW;
                Param.BSIM3gamma2 = ModelParameters.Gamma2
                                    + ModelParameters.Lgamma2 * Inv_L
                                    + ModelParameters.Wgamma2 * Inv_W
                                    + ModelParameters.Pgamma2 * Inv_LW;
                Param.BSIM3vbx = ModelParameters.Vbx
                                 + ModelParameters.Lvbx * Inv_L
                                 + ModelParameters.Wvbx * Inv_W
                                 + ModelParameters.Pvbx * Inv_LW;
                Param.BSIM3vbm = ModelParameters.Vbm
                                 + ModelParameters.Lvbm * Inv_L
                                 + ModelParameters.Wvbm * Inv_W
                                 + ModelParameters.Pvbm * Inv_LW;
                Param.BSIM3xt = ModelParameters.Xt
                                 + ModelParameters.Lxt * Inv_L
                                 + ModelParameters.Wxt * Inv_W
                                 + ModelParameters.Pxt * Inv_LW;
                Param.BSIM3vfb = ModelParameters.Vfb
                                 + ModelParameters.Lvfb * Inv_L
                                 + ModelParameters.Wvfb * Inv_W
                                 + ModelParameters.Pvfb * Inv_LW;
                Param.BSIM3k1 = ModelParameters.K1
                                + ModelParameters.Lk1 * Inv_L
                                + ModelParameters.Wk1 * Inv_W
                                + ModelParameters.Pk1 * Inv_LW;
                Param.BSIM3kt1 = ModelParameters.Kt1
                                 + ModelParameters.Lkt1 * Inv_L
                                 + ModelParameters.Wkt1 * Inv_W
                                 + ModelParameters.Pkt1 * Inv_LW;
                Param.BSIM3kt1l = ModelParameters.Kt1l
                                  + ModelParameters.Lkt1l * Inv_L
                                  + ModelParameters.Wkt1l * Inv_W
                                  + ModelParameters.Pkt1l * Inv_LW;
                Param.BSIM3k2 = ModelParameters.K2
                                + ModelParameters.Lk2 * Inv_L
                                + ModelParameters.Wk2 * Inv_W
                                + ModelParameters.Pk2 * Inv_LW;
                Param.BSIM3kt2 = ModelParameters.Kt2
                                 + ModelParameters.Lkt2 * Inv_L
                                 + ModelParameters.Wkt2 * Inv_W
                                 + ModelParameters.Pkt2 * Inv_LW;
                Param.BSIM3k3 = ModelParameters.K3
                                + ModelParameters.Lk3 * Inv_L
                                + ModelParameters.Wk3 * Inv_W
                                + ModelParameters.Pk3 * Inv_LW;
                Param.BSIM3k3b = ModelParameters.K3b
                                 + ModelParameters.Lk3b * Inv_L
                                 + ModelParameters.Wk3b * Inv_W
                                 + ModelParameters.Pk3b * Inv_LW;
                Param.BSIM3w0 = ModelParameters.W0
                                + ModelParameters.Lw0 * Inv_L
                                + ModelParameters.Ww0 * Inv_W
                                + ModelParameters.Pw0 * Inv_LW;
                Param.BSIM3nlx = ModelParameters.Nlx
                                 + ModelParameters.Lnlx * Inv_L
                                 + ModelParameters.Wnlx * Inv_W
                                 + ModelParameters.Pnlx * Inv_LW;
                Param.BSIM3dvt0 = ModelParameters.Dvt0
                                  + ModelParameters.Ldvt0 * Inv_L
                                  + ModelParameters.Wdvt0 * Inv_W
                                  + ModelParameters.Pdvt0 * Inv_LW;
                Param.BSIM3dvt1 = ModelParameters.Dvt1
                                  + ModelParameters.Ldvt1 * Inv_L
                                  + ModelParameters.Wdvt1 * Inv_W
                                  + ModelParameters.Pdvt1 * Inv_LW;
                Param.BSIM3dvt2 = ModelParameters.Dvt2
                                  + ModelParameters.Ldvt2 * Inv_L
                                  + ModelParameters.Wdvt2 * Inv_W
                                  + ModelParameters.Pdvt2 * Inv_LW;
                Param.BSIM3dvt0w = ModelParameters.Dvt0w
                                  + ModelParameters.Ldvt0w * Inv_L
                                  + ModelParameters.Wdvt0w * Inv_W
                                  + ModelParameters.Pdvt0w * Inv_LW;
                Param.BSIM3dvt1w = ModelParameters.Dvt1w
                                  + ModelParameters.Ldvt1w * Inv_L
                                  + ModelParameters.Wdvt1w * Inv_W
                                  + ModelParameters.Pdvt1w * Inv_LW;
                Param.BSIM3dvt2w = ModelParameters.Dvt2w
                                  + ModelParameters.Ldvt2w * Inv_L
                                  + ModelParameters.Wdvt2w * Inv_W
                                  + ModelParameters.Pdvt2w * Inv_LW;
                Param.BSIM3drout = ModelParameters.Drout
                                   + ModelParameters.Ldrout * Inv_L
                                   + ModelParameters.Wdrout * Inv_W
                                   + ModelParameters.Pdrout * Inv_LW;
                Param.BSIM3dsub = ModelParameters.Dsub
                                  + ModelParameters.Ldsub * Inv_L
                                  + ModelParameters.Wdsub * Inv_W
                                  + ModelParameters.Pdsub * Inv_LW;
                Param.BSIM3vth0 = ModelParameters.Vth0
                                  + ModelParameters.Lvth0 * Inv_L
                                  + ModelParameters.Wvth0 * Inv_W
                                  + ModelParameters.Pvth0 * Inv_LW;
                Param.BSIM3ua = ModelParameters.Ua
                                + ModelParameters.Lua * Inv_L
                                + ModelParameters.Wua * Inv_W
                                + ModelParameters.Pua * Inv_LW;
                Param.BSIM3ua1 = ModelParameters.Ua1
                                 + ModelParameters.Lua1 * Inv_L
                                 + ModelParameters.Wua1 * Inv_W
                                 + ModelParameters.Pua1 * Inv_LW;
                Param.BSIM3ub = ModelParameters.Ub
                                + ModelParameters.Lub * Inv_L
                                + ModelParameters.Wub * Inv_W
                                + ModelParameters.Pub * Inv_LW;
                Param.BSIM3ub1 = ModelParameters.Ub1
                                 + ModelParameters.Lub1 * Inv_L
                                 + ModelParameters.Wub1 * Inv_W
                                 + ModelParameters.Pub1 * Inv_LW;
                Param.BSIM3uc = ModelParameters.Uc
                                + ModelParameters.Luc * Inv_L
                                + ModelParameters.Wuc * Inv_W
                                + ModelParameters.Puc * Inv_LW;
                Param.BSIM3uc1 = ModelParameters.Uc1
                                 + ModelParameters.Luc1 * Inv_L
                                 + ModelParameters.Wuc1 * Inv_W
                                 + ModelParameters.Puc1 * Inv_LW;
                Param.BSIM3u0 = ModelParameters.U0
                                + ModelParameters.Lu0 * Inv_L
                                + ModelParameters.Wu0 * Inv_W
                                + ModelParameters.Pu0 * Inv_LW;
                Param.BSIM3ute = ModelParameters.Ute
                                 + ModelParameters.Lute * Inv_L
                                 + ModelParameters.Wute * Inv_W
                                 + ModelParameters.Pute * Inv_LW;
                Param.BSIM3voff = ModelParameters.Voff
                                  + ModelParameters.Lvoff * Inv_L
                                  + ModelParameters.Wvoff * Inv_W
                                  + ModelParameters.Pvoff * Inv_LW;
                Param.BSIM3delta = ModelParameters.Delta
                                   + ModelParameters.Ldelta * Inv_L
                                   + ModelParameters.Wdelta * Inv_W
                                   + ModelParameters.Pdelta * Inv_LW;
                Param.BSIM3rdsw = ModelParameters.Rdsw
                                  + ModelParameters.Lrdsw * Inv_L
                                  + ModelParameters.Wrdsw * Inv_W
                                  + ModelParameters.Prdsw * Inv_LW;
                Param.BSIM3prwg = ModelParameters.Prwg
                                  + ModelParameters.Lprwg * Inv_L
                                  + ModelParameters.Wprwg * Inv_W
                                  + ModelParameters.Pprwg * Inv_LW;
                Param.BSIM3prwb = ModelParameters.Prwb
                                  + ModelParameters.Lprwb * Inv_L
                                  + ModelParameters.Wprwb * Inv_W
                                  + ModelParameters.Pprwb * Inv_LW;
                Param.BSIM3prt = ModelParameters.Prt
                                  + ModelParameters.Lprt * Inv_L
                                  + ModelParameters.Wprt * Inv_W
                                  + ModelParameters.Pprt * Inv_LW;
                Param.BSIM3eta0 = ModelParameters.Eta0
                                  + ModelParameters.Leta0 * Inv_L
                                  + ModelParameters.Weta0 * Inv_W
                                  + ModelParameters.Peta0 * Inv_LW;
                Param.BSIM3etab = ModelParameters.Etab
                                  + ModelParameters.Letab * Inv_L
                                  + ModelParameters.Wetab * Inv_W
                                  + ModelParameters.Petab * Inv_LW;
                Param.BSIM3pclm = ModelParameters.Pclm
                                  + ModelParameters.Lpclm * Inv_L
                                  + ModelParameters.Wpclm * Inv_W
                                  + ModelParameters.Ppclm * Inv_LW;
                Param.BSIM3pdibl1 = ModelParameters.Pdibl1
                                    + ModelParameters.Lpdibl1 * Inv_L
                                    + ModelParameters.Wpdibl1 * Inv_W
                                    + ModelParameters.Ppdibl1 * Inv_LW;
                Param.BSIM3pdibl2 = ModelParameters.Pdibl2
                                    + ModelParameters.Lpdibl2 * Inv_L
                                    + ModelParameters.Wpdibl2 * Inv_W
                                    + ModelParameters.Ppdibl2 * Inv_LW;
                Param.BSIM3pdiblb = ModelParameters.Pdiblb
                                    + ModelParameters.Lpdiblb * Inv_L
                                    + ModelParameters.Wpdiblb * Inv_W
                                    + ModelParameters.Ppdiblb * Inv_LW;
                Param.BSIM3pscbe1 = ModelParameters.Pscbe1
                                    + ModelParameters.Lpscbe1 * Inv_L
                                    + ModelParameters.Wpscbe1 * Inv_W
                                    + ModelParameters.Ppscbe1 * Inv_LW;
                Param.BSIM3pscbe2 = ModelParameters.Pscbe2
                                    + ModelParameters.Lpscbe2 * Inv_L
                                    + ModelParameters.Wpscbe2 * Inv_W
                                    + ModelParameters.Ppscbe2 * Inv_LW;
                Param.BSIM3pvag = ModelParameters.Pvag
                                  + ModelParameters.Lpvag * Inv_L
                                  + ModelParameters.Wpvag * Inv_W
                                  + ModelParameters.Ppvag * Inv_LW;
                Param.BSIM3wr = ModelParameters.Wr
                                + ModelParameters.Lwr * Inv_L
                                + ModelParameters.Wwr * Inv_W
                                + ModelParameters.Pwr * Inv_LW;
                Param.BSIM3dwg = ModelParameters.Dwg
                                 + ModelParameters.Ldwg * Inv_L
                                 + ModelParameters.Wdwg * Inv_W
                                 + ModelParameters.Pdwg * Inv_LW;
                Param.BSIM3dwb = ModelParameters.Dwb
                                 + ModelParameters.Ldwb * Inv_L
                                 + ModelParameters.Wdwb * Inv_W
                                 + ModelParameters.Pdwb * Inv_LW;
                Param.BSIM3b0 = ModelParameters.B0
                                + ModelParameters.Lb0 * Inv_L
                                + ModelParameters.Wb0 * Inv_W
                                + ModelParameters.Pb0 * Inv_LW;
                Param.BSIM3b1 = ModelParameters.B1
                                + ModelParameters.Lb1 * Inv_L
                                + ModelParameters.Wb1 * Inv_W
                                + ModelParameters.Pb1 * Inv_LW;
                Param.BSIM3alpha0 = ModelParameters.Alpha0
                                    + ModelParameters.Lalpha0 * Inv_L
                                    + ModelParameters.Walpha0 * Inv_W
                                    + ModelParameters.Palpha0 * Inv_LW;
                Param.BSIM3alpha1 = ModelParameters.Alpha1
                                    + ModelParameters.Lalpha1 * Inv_L
                                    + ModelParameters.Walpha1 * Inv_W
                                    + ModelParameters.Palpha1 * Inv_LW;
                Param.BSIM3beta0 = ModelParameters.Beta0
                                   + ModelParameters.Lbeta0 * Inv_L
                                   + ModelParameters.Wbeta0 * Inv_W
                                   + ModelParameters.Pbeta0 * Inv_LW;
                /* CV model */
                Param.BSIM3elm = ModelParameters.Elm
                                + ModelParameters.Lelm * Inv_L
                                + ModelParameters.Welm * Inv_W
                                + ModelParameters.Pelm * Inv_LW;
                Param.BSIM3cgsl = ModelParameters.Cgsl
                                  + ModelParameters.Lcgsl * Inv_L
                                  + ModelParameters.Wcgsl * Inv_W
                                  + ModelParameters.Pcgsl * Inv_LW;
                Param.BSIM3cgdl = ModelParameters.Cgdl
                                  + ModelParameters.Lcgdl * Inv_L
                                  + ModelParameters.Wcgdl * Inv_W
                                  + ModelParameters.Pcgdl * Inv_LW;
                Param.BSIM3ckappa = ModelParameters.Ckappa
                                    + ModelParameters.Lckappa * Inv_L
                                    + ModelParameters.Wckappa * Inv_W
                                    + ModelParameters.Pckappa * Inv_LW;
                Param.BSIM3cf = ModelParameters.Cf
                                + ModelParameters.Lcf * Inv_L
                                + ModelParameters.Wcf * Inv_W
                                + ModelParameters.Pcf * Inv_LW;
                Param.BSIM3clc = ModelParameters.Clc
                                 + ModelParameters.Lclc * Inv_L
                                 + ModelParameters.Wclc * Inv_W
                                 + ModelParameters.Pclc * Inv_LW;
                Param.BSIM3cle = ModelParameters.Cle
                                 + ModelParameters.Lcle * Inv_L
                                 + ModelParameters.Wcle * Inv_W
                                 + ModelParameters.Pcle * Inv_LW;
                Param.BSIM3vfbcv = ModelParameters.Vfbcv
                                   + ModelParameters.Lvfbcv * Inv_L
                                   + ModelParameters.Wvfbcv * Inv_W
                                   + ModelParameters.Pvfbcv * Inv_LW;
                Param.BSIM3acde = ModelParameters.Acde
                                  + ModelParameters.Lacde * Inv_L
                                  + ModelParameters.Wacde * Inv_W
                                  + ModelParameters.Pacde * Inv_LW;
                Param.BSIM3moin = ModelParameters.Moin
                                  + ModelParameters.Lmoin * Inv_L
                                  + ModelParameters.Wmoin * Inv_W
                                  + ModelParameters.Pmoin * Inv_LW;
                Param.BSIM3noff = ModelParameters.Noff
                                  + ModelParameters.Lnoff * Inv_L
                                  + ModelParameters.Wnoff * Inv_W
                                  + ModelParameters.Pnoff * Inv_LW;
                Param.BSIM3voffcv = ModelParameters.Voffcv
                                    + ModelParameters.Lvoffcv * Inv_L
                                    + ModelParameters.Wvoffcv * Inv_W
                                    + ModelParameters.Pvoffcv * Inv_LW;

                Param.BSIM3abulkCVfactor = 1.0 + Math.Pow((Param.BSIM3clc
                                           / Param.BSIM3leffCV),
                                           Param.BSIM3cle);

                T0 = (TRatio - 1.0);
                Param.BSIM3ua = Param.BSIM3ua + Param.BSIM3ua1 * T0;
                Param.BSIM3ub = Param.BSIM3ub + Param.BSIM3ub1 * T0;
                Param.BSIM3uc = Param.BSIM3uc + Param.BSIM3uc1 * T0;
                if (Param.BSIM3u0 > 1.0)
                    Param.BSIM3u0 = Param.BSIM3u0 / 1.0e4;

                Param.BSIM3u0temp = Param.BSIM3u0
                                    * Math.Pow(TRatio, Param.BSIM3ute);
                Param.BSIM3vsattemp = Param.BSIM3vsat - Param.BSIM3at
                                      * T0;
                Param.BSIM3rds0 = (Param.BSIM3rdsw + Param.BSIM3prt * T0)
                                  / Math.Pow(Param.BSIM3weff * 1E6, Param.BSIM3wr);

                if (BSIM3checkModel())
                {
                    throw new SpiceSharpException("Fatal error(s) detected during BSIM3V3.3 parameter checking for {0} in model {1}".FormatString(Name, ModelTemperature.Name));
                }

                Param.BSIM3cgdo = (ModelParameters.Cgdo + Param.BSIM3cf)
                                  * Param.BSIM3weffCV;
                Param.BSIM3cgso = (ModelParameters.Cgso + Param.BSIM3cf)
                                  * Param.BSIM3weffCV;
                Param.BSIM3cgbo = ModelParameters.Cgbo * Param.BSIM3leffCV;

                T0 = Param.BSIM3leffCV * Param.BSIM3leffCV;
                Param.BSIM3tconst = Param.BSIM3u0temp * Param.BSIM3elm / (ModelTemperature.Cox
                                    * Param.BSIM3weffCV * Param.BSIM3leffCV * T0);

                if (!ModelParameters.Npeak.Given && ModelParameters.Gamma1.Given)
                {
                    T0 = Param.BSIM3gamma1 * ModelTemperature.Cox;
                    Param.BSIM3npeak = 3.021E22 * T0 * T0;
                }

                Param.BSIM3phi = 2.0 * Vtm0
                                 * Math.Log(Param.BSIM3npeak / ni);

                Param.BSIM3sqrtPhi = Math.Sqrt(Param.BSIM3phi);
                Param.BSIM3phis3 = Param.BSIM3sqrtPhi * Param.BSIM3phi;

                Param.BSIM3Xdep0 = Math.Sqrt(2.0 * EPSSI / (Charge_q
                                   * Param.BSIM3npeak * 1.0e6))
                                   * Param.BSIM3sqrtPhi;
                Param.BSIM3sqrtXdep0 = Math.Sqrt(Param.BSIM3Xdep0);
                Param.BSIM3litl = Math.Sqrt(3.0 * Param.BSIM3xj
                                  * ModelParameters.Tox);
                Param.BSIM3vbi = Vtm0 * Math.Log(1.0e20
                                 * Param.BSIM3npeak / (ni * ni));
                Param.BSIM3cdep0 = Math.Sqrt(Charge_q * EPSSI
                                   * Param.BSIM3npeak * 1.0e6 / 2.0
                                   / Param.BSIM3phi);

                Param.BSIM3ldeb = Math.Sqrt(EPSSI * Vtm0 / (Charge_q
                                  * Param.BSIM3npeak * 1.0e6)) / 3.0;
                Param.BSIM3acde *= Math.Pow((Param.BSIM3npeak / 2.0e16), -0.25);


                if (ModelParameters.K1.Given || ModelParameters.K2.Given)
                {
                    if (!ModelParameters.K1.Given)
                    {
                        SpiceSharpWarning.Warning(this, "Warning: k1 should be specified with k2.");
                        Param.BSIM3k1 = 0.53;
                    }
                    if (!ModelParameters.K2.Given)
                    {
                        SpiceSharpWarning.Warning(this, "Warning: k2 should be specified with k1.");
                        Param.BSIM3k2 = -0.0186;
                    }
                    if (ModelParameters.Nsub.Given)
                        SpiceSharpWarning.Warning(this, "Warning: nsub is ignored because k1 or k2 is given.");
                    if (ModelParameters.Xt.Given)
                        SpiceSharpWarning.Warning(this, "Warning: xt is ignored because k1 or k2 is given.");
                    if (ModelParameters.Vbx.Given)
                        SpiceSharpWarning.Warning(this, "Warning: vbx is ignored because k1 or k2 is given.");
                    if (ModelParameters.Gamma1.Given)
                        SpiceSharpWarning.Warning(this, "Warning: gamma1 is ignored because k1 or k2 is given.");
                    if (ModelParameters.Gamma2.Given)
                        SpiceSharpWarning.Warning(this, "Warning: gamma2 is ignored because k1 or k2 is given.");
                }
                else
                {
                    if (!ModelParameters.Vbx.Given)
                        Param.BSIM3vbx = Param.BSIM3phi - 7.7348e-4
                                         * Param.BSIM3npeak
                                         * Param.BSIM3xt * Param.BSIM3xt;
                    if (Param.BSIM3vbx > 0.0)
                        Param.BSIM3vbx = -Param.BSIM3vbx;
                    if (Param.BSIM3vbm > 0.0)
                        Param.BSIM3vbm = -Param.BSIM3vbm;

                    if (!ModelParameters.Gamma1.Given)
                        Param.BSIM3gamma1 = 5.753e-12
                                            * Math.Sqrt(Param.BSIM3npeak)
                                            / ModelTemperature.Cox;
                    if (!ModelParameters.Gamma2.Given)
                        Param.BSIM3gamma2 = 5.753e-12
                                            * Math.Sqrt(Param.BSIM3nsub)
                                            / ModelTemperature.Cox;

                    T0 = Param.BSIM3gamma1 - Param.BSIM3gamma2;
                    T1 = Math.Sqrt(Param.BSIM3phi - Param.BSIM3vbx)
                       - Param.BSIM3sqrtPhi;
                    T2 = Math.Sqrt(Param.BSIM3phi * (Param.BSIM3phi
                       - Param.BSIM3vbm)) - Param.BSIM3phi;
                    Param.BSIM3k2 = T0 * T1 / (2.0 * T2 + Param.BSIM3vbm);
                    Param.BSIM3k1 = Param.BSIM3gamma2 - 2.0
                                    * Param.BSIM3k2 * Math.Sqrt(Param.BSIM3phi
                                    - Param.BSIM3vbm);
                }

                if (Param.BSIM3k2 < 0.0)
                {
                    T0 = 0.5 * Param.BSIM3k1 / Param.BSIM3k2;
                    Param.BSIM3vbsc = 0.9 * (Param.BSIM3phi - T0 * T0);
                    if (Param.BSIM3vbsc > -3.0)
                        Param.BSIM3vbsc = -3.0;
                    else if (Param.BSIM3vbsc < -30.0)
                        Param.BSIM3vbsc = -30.0;
                }
                else
                {
                    Param.BSIM3vbsc = -30.0;
                }
                if (Param.BSIM3vbsc > Param.BSIM3vbm)
                    Param.BSIM3vbsc = Param.BSIM3vbm;

                if (!ModelParameters.Vfb.Given)
                {
                    if (ModelParameters.Vth0.Given)
                    {
                        Param.BSIM3vfb = ModelParameters.B3Type * Param.BSIM3vth0
                                         - Param.BSIM3phi - Param.BSIM3k1
                                         * Param.BSIM3sqrtPhi;
                    }
                    else
                    {
                        Param.BSIM3vfb = -1.0;
                    }
                }
                if (!ModelParameters.Vth0.Given)
                {
                    Param.BSIM3vth0 = ModelParameters.B3Type * (Param.BSIM3vfb
                                      + Param.BSIM3phi + Param.BSIM3k1
                                      * Param.BSIM3sqrtPhi);
                }

                Param.BSIM3k1ox = Param.BSIM3k1 * ModelParameters.Tox
                                  / ModelParameters.Toxm;
                Param.BSIM3k2ox = Param.BSIM3k2 * ModelParameters.Tox
                                  / ModelParameters.Toxm;

                T1 = Math.Sqrt(EPSSI / EPSOX * ModelParameters.Tox
                   * Param.BSIM3Xdep0);
                T0 = Math.Exp(-0.5 * Param.BSIM3dsub * Param.BSIM3leff / T1);
                Param.BSIM3theta0vb0 = (T0 + 2.0 * T0 * T0);

                T0 = Math.Exp(-0.5 * Param.BSIM3drout * Param.BSIM3leff / T1);
                T2 = (T0 + 2.0 * T0 * T0);
                Param.BSIM3thetaRout = Param.BSIM3pdibl1 * T2
                                       + Param.BSIM3pdibl2;

                tmp = Math.Sqrt(Param.BSIM3Xdep0);
                tmp1 = Param.BSIM3vbi - Param.BSIM3phi;
                tmp2 = ModelTemperature.Factor1 * tmp;

                T0 = -0.5 * Param.BSIM3dvt1w * Param.BSIM3weff
                   * Param.BSIM3leff / tmp2;
                if (T0 > -EXP_THRESHOLD)
                {
                    T1 = Math.Exp(T0);
                    T2 = T1 * (1.0 + 2.0 * T1);
                }
                else
                {
                    T1 = MIN_EXP;
                    T2 = T1 * (1.0 + 2.0 * T1);
                }
                T0 = Param.BSIM3dvt0w * T2;
                T2 = T0 * tmp1;

                T0 = -0.5 * Param.BSIM3dvt1 * Param.BSIM3leff / tmp2;
                if (T0 > -EXP_THRESHOLD)
                {
                    T1 = Math.Exp(T0);
                    T3 = T1 * (1.0 + 2.0 * T1);
                }
                else
                {
                    T1 = MIN_EXP;
                    T3 = T1 * (1.0 + 2.0 * T1);
                }
                T3 = Param.BSIM3dvt0 * T3 * tmp1;

                T4 = ModelParameters.Tox * Param.BSIM3phi
                   / (Param.BSIM3weff + Param.BSIM3w0);

                T0 = Math.Sqrt(1.0 + Param.BSIM3nlx / Param.BSIM3leff);
                T5 = Param.BSIM3k1ox * (T0 - 1.0) * Param.BSIM3sqrtPhi
                   + (Param.BSIM3kt1 + Param.BSIM3kt1l / Param.BSIM3leff)
                   * (TRatio - 1.0);

                tmp3 = ModelParameters.B3Type * Param.BSIM3vth0
                     - T2 - T3 + Param.BSIM3k3 * T4 + T5;
                Param.BSIM3vfbzb = tmp3 - Param.BSIM3phi - Param.BSIM3k1
                                   * Param.BSIM3sqrtPhi;
                /* End of vfbzb */
            }

            /* adding delvto  */
            Vth0 = Param.BSIM3vth0 + Parameters.Delvto;
            Vfb = Param.BSIM3vfb + ModelParameters.B3Type * Parameters.Delvto;
            Vfbzb = Param.BSIM3vfbzb + ModelParameters.B3Type * Parameters.Delvto;

            /* low field mobility multiplier */
            U0temp = Param.BSIM3u0temp * Parameters.Mulu0;
            Tconst = U0temp * Param.BSIM3elm / (ModelTemperature.Cox
                                  * Param.BSIM3weffCV * Param.BSIM3leffCV * T0);

            /* process source/drain series resistance */
            /* ACM model */
            if (ModelParameters.AcmMod.Value == 0)
            {
                DrainConductance = ModelParameters.SheetResistance
                                                * Parameters.DrainSquares;
                SourceConductance = ModelParameters.SheetResistance
                                             * Parameters.SourceSquares;
            }
            else /* ACM > 0 */
            {
                ACM.SourceDrainResistances(
                    ModelParameters.AcmMod,
                    ModelParameters.Ld,
                    ModelParameters.Ldif,
                    ModelParameters.Hdif,
                    ModelParameters.Wmlt,
                    Parameters.Width,
                    ModelParameters.Xw,
                    ModelParameters.SheetResistance,
                    ModelParameters.Rd,
                    ModelParameters.Rdc,
                    Parameters.DrainSquares,
                    ModelParameters.Rs,
                    ModelParameters.Rsc,
                    Parameters.SourceSquares,
                    out var dc,
                    out var sc
                );
                DrainConductance = dc;
                SourceConductance = sc;
            }
            if (DrainConductance > 0.0)
                DrainConductance = 1.0 / DrainConductance;
            else
                DrainConductance = 0.0;

            if (SourceConductance > 0.0)
                SourceConductance = 1.0 / SourceConductance;
            else
                SourceConductance = 0.0;
            Cgso = Param.BSIM3cgso;
            Cgdo = Param.BSIM3cgdo;

            Nvtm = ModelTemperature.Vtm * ModelParameters.JctEmissionCoeff;
            if ((Parameters.SourceArea <= 0.0) &&
                (Parameters.SourcePerimeter <= 0.0))
            {
                SourceSatCurrent = 1.0e-14;
            }
            else
            {
                SourceSatCurrent = Parameters.SourceArea
                                 * ModelTemperature.JctTempSatCurDensity
                                 + Parameters.SourcePerimeter
                                 * ModelTemperature.JctSidewallTempSatCurDensity;
            }
            if ((SourceSatCurrent > 0.0) && (ModelParameters.Ijth > 0.0))
            {
                Vjsm = Nvtm * Math.Log(ModelParameters.Ijth
                                / SourceSatCurrent + 1.0);
                IsEvjsm = SourceSatCurrent * Math.Exp(Vjsm
                                   / Nvtm);
            }

            if ((Parameters.DrainArea <= 0.0) &&
                (Parameters.DrainPerimeter <= 0.0))
            {
                DrainSatCurrent = 1.0e-14;
            }
            else
            {
                DrainSatCurrent = Parameters.DrainArea
                                * ModelTemperature.JctTempSatCurDensity
                                + Parameters.DrainPerimeter
                                * ModelTemperature.JctSidewallTempSatCurDensity;
            }
            if ((DrainSatCurrent > 0.0) && (ModelParameters.Ijth > 0.0))
            {
                Vjdm = Nvtm * Math.Log(ModelParameters.Ijth
                                / DrainSatCurrent + 1.0);
                IsEvjdm = DrainSatCurrent * Math.Exp(Vjdm
                                   / Nvtm);
            }
        }

        private bool BSIM3checkModel()
        {
            int Fatal_Flag = 0;
            List<string> words = new List<string>();

            if (ModelParameters.Version != "3.3.0" && ModelParameters.Version != "3.30" && ModelParameters.Version != "3.3")
            {
                SpiceSharpWarning.Warning(this, "Warning: This model is BSIM3v3.3.0; you specified a wrong version number.");
                words.Add("Warning: This model is BSIM3v3.3.0; you specified a wrong version number.");
            }

            if (Param.BSIM3nlx < -Param.BSIM3leff)
            {
                words.Add(string.Format("Fatal: Nlx = {0:g} is less than -Leff.",
                    Param.BSIM3nlx));
                Fatal_Flag = 1;
            }

            if (ModelParameters.Tox <= 0.0)
            {
                words.Add(string.Format("Fatal: Tox = {0:g} is not positive.",
                    ModelParameters.Tox));
                Fatal_Flag = 1;
            }

            if (ModelParameters.Toxm <= 0.0)
            {
                words.Add(string.Format("Fatal: Toxm = {0:g} is not positive.",
                        ModelParameters.Toxm));
                Fatal_Flag = 1;
            }

            if (ModelParameters.Lintnoi > Param.BSIM3leff / 2)
            {
                words.Add(string.Format("Fatal: Lintnoi = {0:g} is too large - Leff for noise is negative.",
                    ModelParameters.Lintnoi));
                Fatal_Flag = 1;
            }

            if (Param.BSIM3npeak <= 0.0)
            {
                words.Add(string.Format("Fatal: Nch = {0:g} is not positive.",
                    Param.BSIM3npeak));
                Fatal_Flag = 1;
            }
            if (Param.BSIM3nsub <= 0.0)
            {
                words.Add(string.Format("Fatal: Nsub = {0:g} is not positive.",
                    Param.BSIM3nsub));
                Fatal_Flag = 1;
            }
            if (Param.BSIM3ngate < 0.0)
            {
                words.Add(string.Format("Fatal: Ngate = {0:g} is not positive.",
                        Param.BSIM3ngate));
                Fatal_Flag = 1;
            }
            if (Param.BSIM3ngate > 1.0e25)
            {
                words.Add(string.Format("Fatal: Ngate = {0:g} is too high.",
                    Param.BSIM3ngate));
                Fatal_Flag = 1;
            }
            if (Param.BSIM3xj <= 0.0)
            {
                words.Add(string.Format("Fatal: Xj = {0:g} is not positive.",
                    Param.BSIM3xj));
                Fatal_Flag = 1;
            }

            if (Param.BSIM3dvt1 < 0.0)
            {
                words.Add(string.Format("Fatal: Dvt1 = {0:g} is negative.",
                    Param.BSIM3dvt1));
                Fatal_Flag = 1;
            }

            if (Param.BSIM3dvt1w < 0.0)
            {
                words.Add(string.Format("Fatal: Dvt1w = {0:g} is negative.",
                    Param.BSIM3dvt1w));
                Fatal_Flag = 1;
            }

            if (Param.BSIM3w0 == -Param.BSIM3weff)
            {
                words.Add(string.Format("Fatal: (W0 + Weff) = 0 causing divided-by-zero."));
                Fatal_Flag = 1;
            }

            if (Param.BSIM3dsub < 0.0)
            {
                words.Add(string.Format("Fatal: Dsub = {0:g} is negative.", Param.BSIM3dsub));
                Fatal_Flag = 1;
            }
            if (Param.BSIM3b1 == -Param.BSIM3weff)
            {
                words.Add(string.Format("Fatal: (B1 + Weff) = 0 causing divided-by-zero."));
                Fatal_Flag = 1;
            }
            if (Param.BSIM3u0temp <= 0.0)
            {
                words.Add(string.Format("Fatal: u0 at current temperature = {0:g} is not positive.", Param.BSIM3u0temp));
                Fatal_Flag = 1;
            }

            /* Check delta parameter */
            if (Param.BSIM3delta < 0.0)
            {
                words.Add(string.Format("Fatal: Delta = {0:g} is less than zero.",
                    Param.BSIM3delta));
                Fatal_Flag = 1;
            }

            if (Param.BSIM3vsattemp <= 0.0)
            {
                words.Add(string.Format("Fatal: Vsat at current temperature = {0:g} is not positive.", Param.BSIM3vsattemp));
                Fatal_Flag = 1;
            }
            /* Check Rout parameters */
            if (Param.BSIM3pclm <= 0.0)
            {
                words.Add(string.Format("Fatal: Pclm = {0:g} is not positive.", Param.BSIM3pclm));
                Fatal_Flag = 1;
            }

            if (Param.BSIM3drout < 0.0)
            {
                words.Add(string.Format("Fatal: Drout = {0:g} is negative.", Param.BSIM3drout));
                Fatal_Flag = 1;
            }

            if (Param.BSIM3pscbe2 <= 0.0)
            {
                words.Add(string.Format("Warning: Pscbe2 = {0:g} is not positive.",
                        Param.BSIM3pscbe2));
            }

            /* ACM model */
            if (ModelParameters.AcmMod.Value == 0)
            {
                if (ModelParameters.UnitLengthSidewallJctCap > 0.0 ||
                      ModelParameters.UnitLengthGateSidewallJctCap > 0.0)
                {
                    if (Parameters.DrainPerimeter < Param.BSIM3weff)
                    {
                        words.Add(string.Format("Warning: Pd = {0:g} is less than W.",
                             Parameters.DrainPerimeter));
                    }
                    if (Parameters.SourcePerimeter < Param.BSIM3weff)
                    {
                        words.Add(string.Format("Warning: Ps = {0:g} is less than W.",
                             Parameters.SourcePerimeter));
                    }
                }
            }

            if ((ModelParameters.Calcacm.Value > 0) && (ModelParameters.AcmMod.Value != 12))
            {
                words.Add(string.Format("Warning: CALCACM = {0} is wrong. Set back to 0.",
                    ModelParameters.Calcacm));
                ModelParameters.Calcacm = 0;
            }

            if (Param.BSIM3noff < 0.1)
            {
                words.Add(string.Format("Warning: Noff = {0:g} is too small.",
                        Param.BSIM3noff));
            }
            if (Param.BSIM3noff > 4.0)
            {
                words.Add(string.Format("Warning: Noff = {0:g} is too large.",
                        Param.BSIM3noff));
            }

            if (Param.BSIM3voffcv < -0.5)
            {
                words.Add(string.Format("Warning: Voffcv = {0:g} is too small.",
                        Param.BSIM3voffcv));
            }
            if (Param.BSIM3voffcv > 0.5)
            {
                words.Add(string.Format("Warning: Voffcv = {0:g} is too large.",
                        Param.BSIM3voffcv));
            }

            if (ModelParameters.Ijth < 0.0)
            {
                words.Add(string.Format("Fatal: Ijth = {0:g} cannot be negative.",
                        ModelParameters.Ijth));
                Fatal_Flag = 1;
            }

            /* Check capacitance parameters */
            if (Param.BSIM3clc < 0.0)
            {
                words.Add(string.Format("Fatal: Clc = {0:g} is negative.", Param.BSIM3clc));
                Fatal_Flag = 1;
            }

            if (Param.BSIM3moin < 5.0)
            {
                words.Add(string.Format("Warning: Moin = {0:g} is too small.",
                        Param.BSIM3moin));
            }
            if (Param.BSIM3moin > 25.0)
            {
                words.Add(string.Format("Warning: Moin = {0:g} is too large.",
                        Param.BSIM3moin));
            }

            if (ModelParameters.CapMod.Value == 3)
            {
                if (Param.BSIM3acde < 0.4)
                {
                    words.Add(string.Format("Warning:  Acde = {0:g} is too small.",
                            Param.BSIM3acde));
                }
                if (Param.BSIM3acde > 1.6)
                {
                    words.Add(string.Format("Warning:  Acde = {0:g} is too large.",
                            Param.BSIM3acde));
                }
            }

            if (ModelParameters.ParamChk.Value == 1)
            {
                /* Check L and W parameters */
                if (Param.BSIM3leff <= 5.0e-8)
                {
                    words.Add(string.Format("Warning: Leff = {0:g} may be too small.",
                            Param.BSIM3leff));
                }

                if (Param.BSIM3leffCV <= 5.0e-8)
                {
                    words.Add(string.Format("Warning: Leff for CV = {0:g} may be too small.",
                            Param.BSIM3leffCV));
                }

                if (Param.BSIM3weff <= 1.0e-7)
                {
                    words.Add(string.Format("Warning: Weff = {0:g} may be too small.",
                            Param.BSIM3weff));
                }

                if (Param.BSIM3weffCV <= 1.0e-7)
                {
                    words.Add(string.Format("Warning: Weff for CV = {0:g} may be too small.",
                            Param.BSIM3weffCV));
                }

                /* Check threshold voltage parameters */
                if (Param.BSIM3nlx < 0.0)
                {
                    words.Add(string.Format("Warning: Nlx = {0:g} is negative.", Param.BSIM3nlx));
                }
                if (ModelParameters.Tox < 1.0e-9)
                {
                    words.Add(string.Format("Warning: Tox = {0:g} is less than 10A.",
                            ModelParameters.Tox));
                }

                if (Param.BSIM3npeak <= 1.0e15)
                {
                    words.Add(string.Format("Warning: Nch = {0:g} may be too small.",
                            Param.BSIM3npeak));
                }
                else if (Param.BSIM3npeak >= 1.0e21)
                {
                    words.Add(string.Format("Warning: Nch = {0:g} may be too large.",
                            Param.BSIM3npeak));
                }

                if (Param.BSIM3nsub <= 1.0e14)
                {
                    words.Add(string.Format("Warning: Nsub = {0:g} may be too small.",
                            Param.BSIM3nsub));
                }
                else if (Param.BSIM3nsub >= 1.0e21)
                {
                    words.Add(string.Format("Warning: Nsub = {0:g} may be too large.",
                            Param.BSIM3nsub));
                }

                if ((Param.BSIM3ngate > 0.0) &&
                    (Param.BSIM3ngate <= 1.0e18))
                {
                    words.Add(string.Format("Warning: Ngate = {0:g} is less than 1.E18cm^-3.",
                            Param.BSIM3ngate));
                }

                if (Param.BSIM3dvt0 < 0.0)
                {
                    words.Add(string.Format("Warning: Dvt0 = {0:g} is negative.",
                            Param.BSIM3dvt0));
                }

                if (Math.Abs(1.0e-6 / (Param.BSIM3w0 + Param.BSIM3weff)) > 10.0)
                {
                    words.Add(string.Format("Warning: (W0 + Weff) may be too small."));
                }

                /* Check subthreshold parameters */
                if (Param.BSIM3nfactor < 0.0)
                {
                    words.Add(string.Format("Warning: Nfactor = {0:g} is negative.",
                            Param.BSIM3nfactor));
                }
                if (Param.BSIM3cdsc < 0.0)
                {
                    words.Add(string.Format("Warning: Cdsc = {0:g} is negative.",
                            Param.BSIM3cdsc));
                }
                if (Param.BSIM3cdscd < 0.0)
                {
                    words.Add(string.Format("Warning: Cdscd = {0:g} is negative.",
                            Param.BSIM3cdscd));
                }
                /* Check DIBL parameters */
                if (Param.BSIM3eta0 < 0.0)
                {
                    words.Add(string.Format("Warning: Eta0 = {0:g} is negative.",
                            Param.BSIM3eta0));
                }

                /* Check Abulk parameters */
                if (Math.Abs(1.0e-6 / (Param.BSIM3b1 + Param.BSIM3weff)) > 10.0)
                {
                    words.Add(string.Format("Warning: (B1 + Weff) may be too small."));
                }

                /* Check Saturation parameters */
                if (Param.BSIM3a2 < 0.01)
                {
                    words.Add(string.Format("Warning: A2 = {0:g} is too small. Set to 0.01.", Param.BSIM3a2));
                    Param.BSIM3a2 = 0.01;
                }
                else if (Param.BSIM3a2 > 1.0)
                {
                    words.Add(string.Format("Warning: A2 = {0:g} is larger than 1. A2 is set to 1 and A1 is set to 0.",
                            Param.BSIM3a2));
                    Param.BSIM3a2 = 1.0;
                    Param.BSIM3a1 = 0.0;

                }

                if (Param.BSIM3rdsw < 0.0)
                {
                    words.Add(string.Format("Warning: Rdsw = {0:g} is negative. Set to zero.",
                            Param.BSIM3rdsw));
                    Param.BSIM3rdsw = 0.0;
                    Param.BSIM3rds0 = 0.0;
                }
                else if ((Param.BSIM3rds0 > 0.0) && (Param.BSIM3rds0 < 0.001))
                {
                    words.Add(string.Format("Warning: Rds at current temperature = {0:g} is less than 0.001 ohm. Set to zero.",
                            Param.BSIM3rds0));
                    Param.BSIM3rds0 = 0.0;
                }
                if (Param.BSIM3vsattemp < 1.0e3)
                {
                    words.Add(string.Format("Warning: Vsat at current temperature = {0:g} may be too small.", Param.BSIM3vsattemp));
                }

                if (Param.BSIM3pdibl1 < 0.0)
                {
                    words.Add(string.Format("Warning: Pdibl1 = {0:g} is negative.",
                            Param.BSIM3pdibl1));
                }
                if (Param.BSIM3pdibl2 < 0.0)
                {
                    words.Add(string.Format("Warning: Pdibl2 = {0:g} is negative.",
                            Param.BSIM3pdibl2));
                }
                /* Check overlap capacitance parameters */
                if (ModelParameters.Cgdo < 0.0)
                {
                    words.Add(string.Format("Warning: cgdo = {0:g} is negative. Set to zero.", ModelParameters.Cgdo));
                    ModelParameters.Cgdo = 0.0;
                }
                if (ModelParameters.Cgso < 0.0)
                {
                    words.Add(string.Format("Warning: cgso = {0:g} is negative. Set to zero.", ModelParameters.Cgso));
                    ModelParameters.Cgso = 0.0;
                }
                if (ModelParameters.Cgbo < 0.0)
                {
                    words.Add(string.Format("Warning: cgbo = {0:g} is negative. Set to zero.", ModelParameters.Cgbo));
                    ModelParameters.Cgbo = 0.0;
                }
            }

            if (words.Count > 0)
            {
                string path = ModelParameters.CheckPath;
                if (string.IsNullOrWhiteSpace(path))
                    path = "b3v33check.log";
                using (var writer = new StreamWriter(path))
                {
                    foreach (string line in words)
                        writer.WriteLine(line);
                }
            }
            return Fatal_Flag != 0;
        }
    }
}