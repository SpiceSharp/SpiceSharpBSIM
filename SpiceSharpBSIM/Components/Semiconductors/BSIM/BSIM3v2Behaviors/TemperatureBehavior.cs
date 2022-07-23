using SpiceSharp;
using SpiceSharp.Attributes;
using SpiceSharp.Behaviors;
using SpiceSharp.Components;
using SpiceSharp.ParameterSets;
using System;
using System.Collections.Generic;
using System.IO;
using System.Text;

namespace SpiceSharpBSIM.Components.Semiconductors.BSIM.BSIM3v2Behaviors
{
    /// <summary>
    /// A temperature behavior for a <see cref="BSIM3v2"/>.
    /// </summary>
    [BehaviorFor(typeof(BSIM3v2)), AddBehaviorIfNo(typeof(ITemperatureBehavior))]
    public class TemperatureBehavior : Behavior, ITemperatureBehavior, IParameterized<BaseParameters>
    {
        /// <inheritdoc />
        public BaseParameters Parameters { get; }

        public const double Kb = 1.3806226e-23;
        public const double KboQ = 8.617087e-5; /* Kb / q  where q = 1.60219e-19 */
        public const double EPSOX = 3.453133e-11;
        public const double EPSSI = 1.03594e-10;
        public const double MAX_EXP = 5.834617425e14;
        public const double MIN_EXP = 1.713908431e-15;
        public const double EXP_THRESHOLD = 34.0;
        public const double Charge_q = 1.60219e-19;

        protected ModelTemperatureBehavior ModelTemperature { get; }
        protected ModelParameters ModelParameters { get; }
        protected SizeDependParams Param { get; private set; }
        public double Vth0 { get; private set; }
        public double Vfb { get; private set; }
        public double Vfbzb { get; private set; }
        public double U0temp { get; private set; }
        public double Tconst { get; private set; }
        public double DrainConductance { get; private set; }
        public double SourceConductance { get; private set; }
        public double Cgso { get; set; }
        public double Cgdo { get; set; }
        public double IsEvjsm { get; private set; }
        public double Vjsm { get; private set; }
        public double Vjdm { get; private set; }
        public double IsEvjdm { get; private set; }

        /// <summary>
        /// Creates a new <see cref="TemperatureBehavior"/>.
        /// </summary>
        /// <param name="context"></param>
        public TemperatureBehavior(ComponentBindingContext context)
            : base(context)
        {
            Parameters = context.GetParameterSet<BaseParameters>();
            ModelTemperature = context.ModelBehaviors.GetValue<ModelTemperatureBehavior>();
            ModelParameters = context.ModelBehaviors.GetParameterSet<ModelParameters>();
        }

        /// <inheritdoc />
        void ITemperatureBehavior.Temperature()
        {
            double tmp, tmp1, tmp2, tmp3, ni, T0, T1, T2, T3, T4, T5, Ldrn, Wdrn;
            double TRatio, Inv_L, Inv_W, Inv_LW, Vtm0;
            double Nvtm, SourceSatCurrent, DrainSatCurrent;

            TRatio = ModelTemperature.TRatio;
            Vtm0 = ModelTemperature.Vtm0;
            ni = ModelTemperature.Ni;
            T0 = ModelTemperature.T0;

            /* perform the parameter defaulting */
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
                SpiceSharpWarning.Warning(this, "Warning: nqsMod has been set to its global value {0}.".FormatString(ModelParameters.NqsMod.Value));
            }

            /* Channel length scaling with lmlt model parameter */

            if (ModelParameters.Lmlt.Given)
                Parameters.L *= ModelParameters.Lmlt;

            var key = Tuple.Create(Parameters.W.Value, Parameters.L.Value);
            if (ModelTemperature.SizeDependParams.TryGetValue(key, out var param))
                Param = param;
            else
            {
                Param = new SizeDependParams();
                ModelTemperature.SizeDependParams.Add(key, Param);

                Ldrn = Parameters.L;
                Wdrn = Parameters.W;
                Param.Length = Ldrn;
                Param.Width = Wdrn;

                T0 = Math.Pow(Ldrn, ModelParameters.Lln);
                T1 = Math.Pow(Wdrn, ModelParameters.Lwn);
                tmp1 = ModelParameters.Ll / T0 + ModelParameters.Lw / T1
                               + ModelParameters.Lwl / (T0 * T1);
                Param.BSIM3v32dl = ModelParameters.Lint + tmp1;
                tmp2 = ModelParameters.Llc / T0 + ModelParameters.Lwc / T1
                     + ModelParameters.Lwlc / (T0 * T1);
                Param.BSIM3v32dlc = ModelParameters.Dlc + tmp2;

                T2 = Math.Pow(Ldrn, ModelParameters.Wln);
                T3 = Math.Pow(Wdrn, ModelParameters.Wwn);
                tmp1 = ModelParameters.Wl / T2 + ModelParameters.Ww / T3
                               + ModelParameters.Wwl / (T2 * T3);
                Param.BSIM3v32dw = ModelParameters.Wint + tmp1;
                tmp2 = ModelParameters.Wlc / T2 + ModelParameters.Wwc / T3
                     + ModelParameters.Wwlc / (T2 * T3);
                Param.BSIM3v32dwc = ModelParameters.Dwc + tmp2;

                Param.BSIM3v32leff = Parameters.L + ModelParameters.Xl - 2.0 * Param.BSIM3v32dl;
                if (Param.BSIM3v32leff <= 0.0)
                    throw new SpiceSharpException("BSIM3v32: mosfet {0}, model {1}: Effective channel length <= 0".FormatString(Name, ModelTemperature.Name));

                Param.BSIM3v32weff = Parameters.W + ModelParameters.Xw - 2.0 * Param.BSIM3v32dw;
                if (Param.BSIM3v32weff <= 0.0)
                    throw new SpiceSharpException("BSIM3v32: mosfet {0}, model {1}:Effective channel width <= 0".FormatString(Name, ModelTemperature.Name));
                Param.BSIM3v32leffCV = Parameters.L + ModelParameters.Xl - 2.0 * Param.BSIM3v32dlc;
                if (Param.BSIM3v32leffCV <= 0.0)
                    throw new SpiceSharpException("BSIM3v32: mosfet {0}, model {1}: Effective channel length for C-V <= 0".FormatString(Name, ModelTemperature.Name));

                Param.BSIM3v32weffCV = Parameters.W + ModelParameters.Xw - 2.0 * Param.BSIM3v32dwc;
                if (Param.BSIM3v32weffCV <= 0.0)
                    throw new SpiceSharpException("BSIM3v32: mosfet {0}, model {1}: Effective channel width for C-V <= 0".FormatString(Name, ModelTemperature.Name));

                if (ModelParameters.BinUnit.Value == 1)
                {
                    Inv_L = 1.0e-6 / Param.BSIM3v32leff;
                    Inv_W = 1.0e-6 / Param.BSIM3v32weff;
                    Inv_LW = 1.0e-12 / (Param.BSIM3v32leff
                           * Param.BSIM3v32weff);
                }
                else
                {
                    Inv_L = 1.0 / Param.BSIM3v32leff;
                    Inv_W = 1.0 / Param.BSIM3v32weff;
                    Inv_LW = 1.0 / (Param.BSIM3v32leff
                           * Param.BSIM3v32weff);
                }
                Param.BSIM3v32cdsc = ModelParameters.Cdsc
                                  + ModelParameters.Lcdsc * Inv_L
                                  + ModelParameters.Wcdsc * Inv_W
                                  + ModelParameters.Pcdsc * Inv_LW;
                Param.BSIM3v32cdscb = ModelParameters.Cdscb
                                   + ModelParameters.Lcdscb * Inv_L
                                   + ModelParameters.Wcdscb * Inv_W
                                   + ModelParameters.Pcdscb * Inv_LW;

                Param.BSIM3v32cdscd = ModelParameters.Cdscd
                                   + ModelParameters.Lcdscd * Inv_L
                                   + ModelParameters.Wcdscd * Inv_W
                                   + ModelParameters.Pcdscd * Inv_LW;

                Param.BSIM3v32cit = ModelParameters.Cit
                                 + ModelParameters.Lcit * Inv_L
                                 + ModelParameters.Wcit * Inv_W
                                 + ModelParameters.Pcit * Inv_LW;
                Param.BSIM3v32nfactor = ModelParameters.Nfactor
                                     + ModelParameters.Lnfactor * Inv_L
                                     + ModelParameters.Wnfactor * Inv_W
                                     + ModelParameters.Pnfactor * Inv_LW;
                Param.BSIM3v32xj = ModelParameters.Xj
                                + ModelParameters.Lxj * Inv_L
                                + ModelParameters.Wxj * Inv_W
                                + ModelParameters.Pxj * Inv_LW;
                Param.BSIM3v32vsat = ModelParameters.Vsat
                                  + ModelParameters.Lvsat * Inv_L
                                  + ModelParameters.Wvsat * Inv_W
                                  + ModelParameters.Pvsat * Inv_LW;
                Param.BSIM3v32at = ModelParameters.At
                                + ModelParameters.Lat * Inv_L
                                + ModelParameters.Wat * Inv_W
                                + ModelParameters.Pat * Inv_LW;
                Param.BSIM3v32a0 = ModelParameters.A0
                                + ModelParameters.La0 * Inv_L
                                + ModelParameters.Wa0 * Inv_W
                                + ModelParameters.Pa0 * Inv_LW;

                Param.BSIM3v32ags = ModelParameters.Ags
                                + ModelParameters.Lags * Inv_L
                                + ModelParameters.Wags * Inv_W
                                + ModelParameters.Pags * Inv_LW;

                Param.BSIM3v32a1 = ModelParameters.A1
                                + ModelParameters.La1 * Inv_L
                                + ModelParameters.Wa1 * Inv_W
                                + ModelParameters.Pa1 * Inv_LW;
                Param.BSIM3v32a2 = ModelParameters.A2
                                + ModelParameters.La2 * Inv_L
                                + ModelParameters.Wa2 * Inv_W
                                + ModelParameters.Pa2 * Inv_LW;
                Param.BSIM3v32keta = ModelParameters.Keta
                                  + ModelParameters.Lketa * Inv_L
                                  + ModelParameters.Wketa * Inv_W
                                  + ModelParameters.Pketa * Inv_LW;
                Param.BSIM3v32nsub = ModelParameters.Nsub
                                  + ModelParameters.Lnsub * Inv_L
                                  + ModelParameters.Wnsub * Inv_W
                                  + ModelParameters.Pnsub * Inv_LW;
                Param.BSIM3v32npeak = ModelParameters.Npeak
                                   + ModelParameters.Lnpeak * Inv_L
                                   + ModelParameters.Wnpeak * Inv_W
                                   + ModelParameters.Pnpeak * Inv_LW;
                Param.BSIM3v32ngate = ModelParameters.Ngate
                                   + ModelParameters.Lngate * Inv_L
                                   + ModelParameters.Wngate * Inv_W
                                   + ModelParameters.Pngate * Inv_LW;
                Param.BSIM3v32gamma1 = ModelParameters.Gamma1
                                    + ModelParameters.Lgamma1 * Inv_L
                                    + ModelParameters.Wgamma1 * Inv_W
                                    + ModelParameters.Pgamma1 * Inv_LW;
                Param.BSIM3v32gamma2 = ModelParameters.Gamma2
                                    + ModelParameters.Lgamma2 * Inv_L
                                    + ModelParameters.Wgamma2 * Inv_W
                                    + ModelParameters.Pgamma2 * Inv_LW;
                Param.BSIM3v32vbx = ModelParameters.Vbx
                                 + ModelParameters.Lvbx * Inv_L
                                 + ModelParameters.Wvbx * Inv_W
                                 + ModelParameters.Pvbx * Inv_LW;
                Param.BSIM3v32vbm = ModelParameters.Vbm
                                 + ModelParameters.Lvbm * Inv_L
                                 + ModelParameters.Wvbm * Inv_W
                                 + ModelParameters.Pvbm * Inv_LW;
                Param.BSIM3v32xt = ModelParameters.Xt
                                 + ModelParameters.Lxt * Inv_L
                                 + ModelParameters.Wxt * Inv_W
                                 + ModelParameters.Pxt * Inv_LW;
                Param.BSIM3v32vfb = ModelParameters.Vfb
                           + ModelParameters.Lvfb * Inv_L
                           + ModelParameters.Wvfb * Inv_W
                           + ModelParameters.Pvfb * Inv_LW;
                Param.BSIM3v32k1 = ModelParameters.K1
                                + ModelParameters.Lk1 * Inv_L
                                + ModelParameters.Wk1 * Inv_W
                                + ModelParameters.Pk1 * Inv_LW;
                Param.BSIM3v32kt1 = ModelParameters.Kt1
                                 + ModelParameters.Lkt1 * Inv_L
                                 + ModelParameters.Wkt1 * Inv_W
                                 + ModelParameters.Pkt1 * Inv_LW;
                Param.BSIM3v32kt1l = ModelParameters.Kt1l
                                  + ModelParameters.Lkt1l * Inv_L
                                  + ModelParameters.Wkt1l * Inv_W
                                  + ModelParameters.Pkt1l * Inv_LW;
                Param.BSIM3v32k2 = ModelParameters.K2
                                + ModelParameters.Lk2 * Inv_L
                                + ModelParameters.Wk2 * Inv_W
                                + ModelParameters.Pk2 * Inv_LW;
                Param.BSIM3v32kt2 = ModelParameters.Kt2
                                 + ModelParameters.Lkt2 * Inv_L
                                 + ModelParameters.Wkt2 * Inv_W
                                 + ModelParameters.Pkt2 * Inv_LW;
                Param.BSIM3v32k3 = ModelParameters.K3
                                + ModelParameters.Lk3 * Inv_L
                                + ModelParameters.Wk3 * Inv_W
                                + ModelParameters.Pk3 * Inv_LW;
                Param.BSIM3v32k3b = ModelParameters.K3b
                                 + ModelParameters.Lk3b * Inv_L
                                 + ModelParameters.Wk3b * Inv_W
                                 + ModelParameters.Pk3b * Inv_LW;
                Param.BSIM3v32w0 = ModelParameters.W0
                                + ModelParameters.Lw0 * Inv_L
                                + ModelParameters.Ww0 * Inv_W
                                + ModelParameters.Pw0 * Inv_LW;
                Param.BSIM3v32nlx = ModelParameters.Nlx
                                 + ModelParameters.Lnlx * Inv_L
                                 + ModelParameters.Wnlx * Inv_W
                                 + ModelParameters.Pnlx * Inv_LW;
                Param.BSIM3v32dvt0 = ModelParameters.Dvt0
                                  + ModelParameters.Ldvt0 * Inv_L
                                  + ModelParameters.Wdvt0 * Inv_W
                                  + ModelParameters.Pdvt0 * Inv_LW;
                Param.BSIM3v32dvt1 = ModelParameters.Dvt1
                                  + ModelParameters.Ldvt1 * Inv_L
                                  + ModelParameters.Wdvt1 * Inv_W
                                  + ModelParameters.Pdvt1 * Inv_LW;
                Param.BSIM3v32dvt2 = ModelParameters.Dvt2
                                  + ModelParameters.Ldvt2 * Inv_L
                                  + ModelParameters.Wdvt2 * Inv_W
                                  + ModelParameters.Pdvt2 * Inv_LW;
                Param.BSIM3v32dvt0w = ModelParameters.Dvt0w
                                  + ModelParameters.Ldvt0w * Inv_L
                                  + ModelParameters.Wdvt0w * Inv_W
                                  + ModelParameters.Pdvt0w * Inv_LW;
                Param.BSIM3v32dvt1w = ModelParameters.Dvt1w
                                  + ModelParameters.Ldvt1w * Inv_L
                                  + ModelParameters.Wdvt1w * Inv_W
                                  + ModelParameters.Pdvt1w * Inv_LW;
                Param.BSIM3v32dvt2w = ModelParameters.Dvt2w
                                  + ModelParameters.Ldvt2w * Inv_L
                                  + ModelParameters.Wdvt2w * Inv_W
                                  + ModelParameters.Pdvt2w * Inv_LW;
                Param.BSIM3v32drout = ModelParameters.Drout
                                   + ModelParameters.Ldrout * Inv_L
                                   + ModelParameters.Wdrout * Inv_W
                                   + ModelParameters.Pdrout * Inv_LW;
                Param.BSIM3v32dsub = ModelParameters.Dsub
                                  + ModelParameters.Ldsub * Inv_L
                                  + ModelParameters.Wdsub * Inv_W
                                  + ModelParameters.Pdsub * Inv_LW;
                Param.BSIM3v32vth0 = ModelParameters.Vth0
                                  + ModelParameters.Lvth0 * Inv_L
                                  + ModelParameters.Wvth0 * Inv_W
                                  + ModelParameters.Pvth0 * Inv_LW;
                Param.BSIM3v32ua = ModelParameters.Ua
                                + ModelParameters.Lua * Inv_L
                                + ModelParameters.Wua * Inv_W
                                + ModelParameters.Pua * Inv_LW;
                Param.BSIM3v32ua1 = ModelParameters.Ua1
                                 + ModelParameters.Lua1 * Inv_L
                                 + ModelParameters.Wua1 * Inv_W
                                 + ModelParameters.Pua1 * Inv_LW;
                Param.BSIM3v32ub = ModelParameters.Ub
                                + ModelParameters.Lub * Inv_L
                                + ModelParameters.Wub * Inv_W
                                + ModelParameters.Pub * Inv_LW;
                Param.BSIM3v32ub1 = ModelParameters.Ub1
                                 + ModelParameters.Lub1 * Inv_L
                                 + ModelParameters.Wub1 * Inv_W
                                 + ModelParameters.Pub1 * Inv_LW;
                Param.BSIM3v32uc = ModelParameters.Uc
                                + ModelParameters.Luc * Inv_L
                                + ModelParameters.Wuc * Inv_W
                                + ModelParameters.Puc * Inv_LW;
                Param.BSIM3v32uc1 = ModelParameters.Uc1
                                 + ModelParameters.Luc1 * Inv_L
                                 + ModelParameters.Wuc1 * Inv_W
                                 + ModelParameters.Puc1 * Inv_LW;
                Param.BSIM3v32u0 = ModelParameters.U0
                                + ModelParameters.Lu0 * Inv_L
                                + ModelParameters.Wu0 * Inv_W
                                + ModelParameters.Pu0 * Inv_LW;
                Param.BSIM3v32ute = ModelParameters.Ute
                                 + ModelParameters.Lute * Inv_L
                                 + ModelParameters.Wute * Inv_W
                                 + ModelParameters.Pute * Inv_LW;
                Param.BSIM3v32voff = ModelParameters.Voff
                                  + ModelParameters.Lvoff * Inv_L
                                  + ModelParameters.Wvoff * Inv_W
                                  + ModelParameters.Pvoff * Inv_LW;
                Param.BSIM3v32delta = ModelParameters.Delta
                                   + ModelParameters.Ldelta * Inv_L
                                   + ModelParameters.Wdelta * Inv_W
                                   + ModelParameters.Pdelta * Inv_LW;
                Param.BSIM3v32rdsw = ModelParameters.Rdsw
                                  + ModelParameters.Lrdsw * Inv_L
                                  + ModelParameters.Wrdsw * Inv_W
                                  + ModelParameters.Prdsw * Inv_LW;
                Param.BSIM3v32prwg = ModelParameters.Prwg
                                  + ModelParameters.Lprwg * Inv_L
                                  + ModelParameters.Wprwg * Inv_W
                                  + ModelParameters.Pprwg * Inv_LW;
                Param.BSIM3v32prwb = ModelParameters.Prwb
                                  + ModelParameters.Lprwb * Inv_L
                                  + ModelParameters.Wprwb * Inv_W
                                  + ModelParameters.Pprwb * Inv_LW;
                Param.BSIM3v32prt = ModelParameters.Prt
                                  + ModelParameters.Lprt * Inv_L
                                  + ModelParameters.Wprt * Inv_W
                                  + ModelParameters.Pprt * Inv_LW;
                Param.BSIM3v32eta0 = ModelParameters.Eta0
                                  + ModelParameters.Leta0 * Inv_L
                                  + ModelParameters.Weta0 * Inv_W
                                  + ModelParameters.Peta0 * Inv_LW;
                Param.BSIM3v32etab = ModelParameters.Etab
                                  + ModelParameters.Letab * Inv_L
                                  + ModelParameters.Wetab * Inv_W
                                  + ModelParameters.Petab * Inv_LW;
                Param.BSIM3v32pclm = ModelParameters.Pclm
                                  + ModelParameters.Lpclm * Inv_L
                                  + ModelParameters.Wpclm * Inv_W
                                  + ModelParameters.Ppclm * Inv_LW;
                Param.BSIM3v32pdibl1 = ModelParameters.Pdibl1
                                    + ModelParameters.Lpdibl1 * Inv_L
                                    + ModelParameters.Wpdibl1 * Inv_W
                                    + ModelParameters.Ppdibl1 * Inv_LW;
                Param.BSIM3v32pdibl2 = ModelParameters.Pdibl2
                                    + ModelParameters.Lpdibl2 * Inv_L
                                    + ModelParameters.Wpdibl2 * Inv_W
                                    + ModelParameters.Ppdibl2 * Inv_LW;
                Param.BSIM3v32pdiblb = ModelParameters.Pdiblb
                                    + ModelParameters.Lpdiblb * Inv_L
                                    + ModelParameters.Wpdiblb * Inv_W
                                    + ModelParameters.Ppdiblb * Inv_LW;
                Param.BSIM3v32pscbe1 = ModelParameters.Pscbe1
                                    + ModelParameters.Lpscbe1 * Inv_L
                                    + ModelParameters.Wpscbe1 * Inv_W
                                    + ModelParameters.Ppscbe1 * Inv_LW;
                Param.BSIM3v32pscbe2 = ModelParameters.Pscbe2
                                    + ModelParameters.Lpscbe2 * Inv_L
                                    + ModelParameters.Wpscbe2 * Inv_W
                                    + ModelParameters.Ppscbe2 * Inv_LW;
                Param.BSIM3v32pvag = ModelParameters.Pvag
                                  + ModelParameters.Lpvag * Inv_L
                                  + ModelParameters.Wpvag * Inv_W
                                  + ModelParameters.Ppvag * Inv_LW;
                Param.BSIM3v32wr = ModelParameters.Wr
                                + ModelParameters.Lwr * Inv_L
                                + ModelParameters.Wwr * Inv_W
                                + ModelParameters.Pwr * Inv_LW;
                Param.BSIM3v32dwg = ModelParameters.Dwg
                                 + ModelParameters.Ldwg * Inv_L
                                 + ModelParameters.Wdwg * Inv_W
                                 + ModelParameters.Pdwg * Inv_LW;
                Param.BSIM3v32dwb = ModelParameters.Dwb
                                 + ModelParameters.Ldwb * Inv_L
                                 + ModelParameters.Wdwb * Inv_W
                                 + ModelParameters.Pdwb * Inv_LW;
                Param.BSIM3v32b0 = ModelParameters.B0
                                + ModelParameters.Lb0 * Inv_L
                                + ModelParameters.Wb0 * Inv_W
                                + ModelParameters.Pb0 * Inv_LW;
                Param.BSIM3v32b1 = ModelParameters.B1
                                + ModelParameters.Lb1 * Inv_L
                                + ModelParameters.Wb1 * Inv_W
                                + ModelParameters.Pb1 * Inv_LW;
                Param.BSIM3v32alpha0 = ModelParameters.Alpha0
                                    + ModelParameters.Lalpha0 * Inv_L
                                    + ModelParameters.Walpha0 * Inv_W
                                    + ModelParameters.Palpha0 * Inv_LW;
                Param.BSIM3v32alpha1 = ModelParameters.Alpha1
                              + ModelParameters.Lalpha1 * Inv_L
                              + ModelParameters.Walpha1 * Inv_W
                              + ModelParameters.Palpha1 * Inv_LW;
                Param.BSIM3v32beta0 = ModelParameters.Beta0
                                   + ModelParameters.Lbeta0 * Inv_L
                                   + ModelParameters.Wbeta0 * Inv_W
                                   + ModelParameters.Pbeta0 * Inv_LW;
                /* CV model */
                Param.BSIM3v32elm = ModelParameters.Elm
                                + ModelParameters.Lelm * Inv_L
                                + ModelParameters.Welm * Inv_W
                                + ModelParameters.Pelm * Inv_LW;
                Param.BSIM3v32cgsl = ModelParameters.Cgsl
                                  + ModelParameters.Lcgsl * Inv_L
                                  + ModelParameters.Wcgsl * Inv_W
                                  + ModelParameters.Pcgsl * Inv_LW;
                Param.BSIM3v32cgdl = ModelParameters.Cgdl
                                  + ModelParameters.Lcgdl * Inv_L
                                  + ModelParameters.Wcgdl * Inv_W
                                  + ModelParameters.Pcgdl * Inv_LW;
                Param.BSIM3v32ckappa = ModelParameters.Ckappa
                                    + ModelParameters.Lckappa * Inv_L
                                    + ModelParameters.Wckappa * Inv_W
                                    + ModelParameters.Pckappa * Inv_LW;
                Param.BSIM3v32cf = ModelParameters.Cf
                                + ModelParameters.Lcf * Inv_L
                                + ModelParameters.Wcf * Inv_W
                                + ModelParameters.Pcf * Inv_LW;
                Param.BSIM3v32clc = ModelParameters.Clc
                                 + ModelParameters.Lclc * Inv_L
                                 + ModelParameters.Wclc * Inv_W
                                 + ModelParameters.Pclc * Inv_LW;
                Param.BSIM3v32cle = ModelParameters.Cle
                                 + ModelParameters.Lcle * Inv_L
                                 + ModelParameters.Wcle * Inv_W
                                 + ModelParameters.Pcle * Inv_LW;
                Param.BSIM3v32vfbcv = ModelParameters.Vfbcv
                                   + ModelParameters.Lvfbcv * Inv_L
                                   + ModelParameters.Wvfbcv * Inv_W
                                   + ModelParameters.Pvfbcv * Inv_LW;
                Param.BSIM3v32acde = ModelParameters.Acde
                       + ModelParameters.Lacde * Inv_L
                       + ModelParameters.Wacde * Inv_W
                       + ModelParameters.Pacde * Inv_LW;
                Param.BSIM3v32moin = ModelParameters.Moin
                       + ModelParameters.Lmoin * Inv_L
                       + ModelParameters.Wmoin * Inv_W
                       + ModelParameters.Pmoin * Inv_LW;
                Param.BSIM3v32noff = ModelParameters.Noff
                       + ModelParameters.Lnoff * Inv_L
                       + ModelParameters.Wnoff * Inv_W
                       + ModelParameters.Pnoff * Inv_LW;
                Param.BSIM3v32voffcv = ModelParameters.Voffcv
                       + ModelParameters.Lvoffcv * Inv_L
                       + ModelParameters.Wvoffcv * Inv_W
                       + ModelParameters.Pvoffcv * Inv_LW;

                Param.BSIM3v32abulkCVfactor = 1.0 + Math.Pow((Param.BSIM3v32clc
                             / Param.BSIM3v32leffCV),
                             Param.BSIM3v32cle);

                T0 = (TRatio - 1.0);
                Param.BSIM3v32ua = Param.BSIM3v32ua + Param.BSIM3v32ua1 * T0;
                Param.BSIM3v32ub = Param.BSIM3v32ub + Param.BSIM3v32ub1 * T0;
                Param.BSIM3v32uc = Param.BSIM3v32uc + Param.BSIM3v32uc1 * T0;
                if (Param.BSIM3v32u0 > 1.0)
                    Param.BSIM3v32u0 = Param.BSIM3v32u0 / 1.0e4;

                Param.BSIM3v32u0temp = Param.BSIM3v32u0
                                    * Math.Pow(TRatio, Param.BSIM3v32ute);
                Param.BSIM3v32vsattemp = Param.BSIM3v32vsat - Param.BSIM3v32at
                                      * T0;
                Param.BSIM3v32rds0 = (Param.BSIM3v32rdsw + Param.BSIM3v32prt * T0)
                                  / Math.Pow(Param.BSIM3v32weff * 1E6, Param.BSIM3v32wr);

                if (BSIM3v2checkModel())
                {
                    throw new SpiceSharpException("Fatal error(s) detected during BSIM3V3.3 parameter checking for {0} in model {1}".FormatString(Name, ModelTemperature.Name));
                }

                Param.BSIM3v32cgdo = (ModelParameters.Cgdo + Param.BSIM3v32cf)
                  * Param.BSIM3v32weffCV;
                Param.BSIM3v32cgso = (ModelParameters.Cgso + Param.BSIM3v32cf)
                                  * Param.BSIM3v32weffCV;
                Param.BSIM3v32cgbo = ModelParameters.Cgbo * Param.BSIM3v32leffCV;

                T0 = Param.BSIM3v32leffCV * Param.BSIM3v32leffCV;
                Param.BSIM3v32tconst = Param.BSIM3v32u0temp * Param.BSIM3v32elm / (ModelTemperature.Cox
                                    * Param.BSIM3v32weffCV * Param.BSIM3v32leffCV * T0);

                if (!ModelParameters.Npeak.Given && ModelParameters.Gamma1.Given)
                {
                    T0 = Param.BSIM3v32gamma1 * ModelTemperature.Cox;
                    Param.BSIM3v32npeak = 3.021E22 * T0 * T0;
                }

                Param.BSIM3v32phi = 2.0 * Vtm0
                                 * Math.Log(Param.BSIM3v32npeak / ni);

                Param.BSIM3v32sqrtPhi = Math.Sqrt(Param.BSIM3v32phi);
                Param.BSIM3v32phis3 = Param.BSIM3v32sqrtPhi * Param.BSIM3v32phi;

                Param.BSIM3v32Xdep0 = Math.Sqrt(2.0 * EPSSI / (Charge_q
                                   * Param.BSIM3v32npeak * 1.0e6))
                                   * Param.BSIM3v32sqrtPhi;
                Param.BSIM3v32sqrtXdep0 = Math.Sqrt(Param.BSIM3v32Xdep0);
                Param.BSIM3v32litl = Math.Sqrt(3.0 * Param.BSIM3v32xj
                                  * ModelParameters.Tox);
                Param.BSIM3v32vbi = Vtm0 * Math.Log(1.0e20
                                 * Param.BSIM3v32npeak / (ni * ni));
                Param.BSIM3v32cdep0 = Math.Sqrt(Charge_q * EPSSI
                                   * Param.BSIM3v32npeak * 1.0e6 / 2.0
                                   / Param.BSIM3v32phi);

                Param.BSIM3v32ldeb = Math.Sqrt(EPSSI * Vtm0 / (Charge_q
                                  * Param.BSIM3v32npeak * 1.0e6)) / 3.0;
                Param.BSIM3v32acde *= Math.Pow((Param.BSIM3v32npeak / 2.0e16), -0.25);


                if (ModelParameters.K1.Given || ModelParameters.K2.Given)
                {
                    if (!ModelParameters.K1.Given)
                    {
                        SpiceSharpWarning.Warning(this, "Warning: k1 should be specified with k2.");
                        Param.BSIM3v32k1 = 0.53;
                    }
                    if (!ModelParameters.K2.Given)
                    {
                        SpiceSharpWarning.Warning(this, "Warning: k2 should be specified with k1.");
                        Param.BSIM3v32k2 = -0.0186;
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
                        Param.BSIM3v32vbx = Param.BSIM3v32phi - 7.7348e-4
                                         * Param.BSIM3v32npeak
                                         * Param.BSIM3v32xt * Param.BSIM3v32xt;
                    if (Param.BSIM3v32vbx > 0.0)
                        Param.BSIM3v32vbx = -Param.BSIM3v32vbx;
                    if (Param.BSIM3v32vbm > 0.0)
                        Param.BSIM3v32vbm = -Param.BSIM3v32vbm;

                    if (!ModelParameters.Gamma1.Given)
                        Param.BSIM3v32gamma1 = 5.753e-12
                                            * Math.Sqrt(Param.BSIM3v32npeak)
                                            / ModelTemperature.Cox;
                    if (!ModelParameters.Gamma2.Given)
                        Param.BSIM3v32gamma2 = 5.753e-12
                                            * Math.Sqrt(Param.BSIM3v32nsub)
                                            / ModelTemperature.Cox;

                    T0 = Param.BSIM3v32gamma1 - Param.BSIM3v32gamma2;
                    T1 = Math.Sqrt(Param.BSIM3v32phi - Param.BSIM3v32vbx)
                       - Param.BSIM3v32sqrtPhi;
                    T2 = Math.Sqrt(Param.BSIM3v32phi * (Param.BSIM3v32phi
                       - Param.BSIM3v32vbm)) - Param.BSIM3v32phi;
                    Param.BSIM3v32k2 = T0 * T1 / (2.0 * T2 + Param.BSIM3v32vbm);
                    Param.BSIM3v32k1 = Param.BSIM3v32gamma2 - 2.0
                                    * Param.BSIM3v32k2 * Math.Sqrt(Param.BSIM3v32phi
                                    - Param.BSIM3v32vbm);
                }

                if (Param.BSIM3v32k2 < 0.0)
                {
                    T0 = 0.5 * Param.BSIM3v32k1 / Param.BSIM3v32k2;
                    Param.BSIM3v32vbsc = 0.9 * (Param.BSIM3v32phi - T0 * T0);
                    if (Param.BSIM3v32vbsc > -3.0)
                        Param.BSIM3v32vbsc = -3.0;
                    else if (Param.BSIM3v32vbsc < -30.0)
                        Param.BSIM3v32vbsc = -30.0;
                }
                else
                {
                    Param.BSIM3v32vbsc = -30.0;
                }
                if (Param.BSIM3v32vbsc > Param.BSIM3v32vbm)
                    Param.BSIM3v32vbsc = Param.BSIM3v32vbm;

                if (!ModelParameters.Vfb.Given)
                {
                    if (ModelParameters.Vth0.Given)
                    {
                        Param.BSIM3v32vfb = ModelParameters.Type * Param.BSIM3v32vth0
                                         - Param.BSIM3v32phi - Param.BSIM3v32k1
                                         * Param.BSIM3v32sqrtPhi;
                    }
                    else
                    {
                        Param.BSIM3v32vfb = -1.0;
                    }
                }
                if (!ModelParameters.Vth0.Given)
                {
                    Param.BSIM3v32vth0 = ModelParameters.Type * (Param.BSIM3v32vfb
                                      + Param.BSIM3v32phi + Param.BSIM3v32k1
                                      * Param.BSIM3v32sqrtPhi);
                }

                Param.BSIM3v32k1ox = Param.BSIM3v32k1 * ModelParameters.Tox
                                  / ModelParameters.Toxm;
                Param.BSIM3v32k2ox = Param.BSIM3v32k2 * ModelParameters.Tox
                                  / ModelParameters.Toxm;

                T1 = Math.Sqrt(EPSSI / EPSOX * ModelParameters.Tox
                   * Param.BSIM3v32Xdep0);
                T0 = Math.Exp(-0.5 * Param.BSIM3v32dsub * Param.BSIM3v32leff / T1);
                Param.BSIM3v32theta0vb0 = (T0 + 2.0 * T0 * T0);

                T0 = Math.Exp(-0.5 * Param.BSIM3v32drout * Param.BSIM3v32leff / T1);
                T2 = (T0 + 2.0 * T0 * T0);
                Param.BSIM3v32thetaRout = Param.BSIM3v32pdibl1 * T2
                                       + Param.BSIM3v32pdibl2;

                tmp = Math.Sqrt(Param.BSIM3v32Xdep0);
                tmp1 = Param.BSIM3v32vbi - Param.BSIM3v32phi;
                tmp2 = ModelTemperature.Factor1 * tmp;

                T0 = -0.5 * Param.BSIM3v32dvt1w * Param.BSIM3v32weff
                   * Param.BSIM3v32leff / tmp2;
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
                T0 = Param.BSIM3v32dvt0w * T2;
                T2 = T0 * tmp1;

                T0 = -0.5 * Param.BSIM3v32dvt1 * Param.BSIM3v32leff / tmp2;
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
                T3 = Param.BSIM3v32dvt0 * T3 * tmp1;

                T4 = ModelParameters.Tox * Param.BSIM3v32phi
                   / (Param.BSIM3v32weff + Param.BSIM3v32w0);

                T0 = Math.Sqrt(1.0 + Param.BSIM3v32nlx / Param.BSIM3v32leff);
                T5 = Param.BSIM3v32k1ox * (T0 - 1.0) * Param.BSIM3v32sqrtPhi
                   + (Param.BSIM3v32kt1 + Param.BSIM3v32kt1l / Param.BSIM3v32leff)
                   * (TRatio - 1.0);

                tmp3 = ModelParameters.Type * Param.BSIM3v32vth0
                     - T2 - T3 + Param.BSIM3v32k3 * T4 + T5;
                Param.BSIM3v32vfbzb = tmp3 - Param.BSIM3v32phi - Param.BSIM3v32k1
                                   * Param.BSIM3v32sqrtPhi;
                /* End of vfbzb */
            }

            /* adding delvto  */
            this.Vth0 = Param.BSIM3v32vth0 + Parameters.Delvto;
            this.Vfb = Param.BSIM3v32vfb + ModelParameters.Type * Parameters.Delvto;
            this.Vfbzb = Param.BSIM3v32vfbzb + ModelParameters.Type * Parameters.Delvto;

            /* low field mobility multiplier */
            this.U0temp = Param.BSIM3v32u0temp * Parameters.Mulu0;

            this.Tconst = this.U0temp * Param.BSIM3v32elm / (ModelTemperature.Cox
                                  * Param.BSIM3v32weffCV * Param.BSIM3v32leffCV * T0);

            /* process source/drain series resistance */
            /* ACM model */

            double DrainResistance, SourceResistance;

            if (ModelParameters.AcmMod.Value == 0)
            {
                DrainResistance = ModelParameters.SheetResistance
                                                * Parameters.DrainSquares;
                SourceResistance = ModelParameters.SheetResistance
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
                    Parameters.W,
                    ModelParameters.Xw,
                    ModelParameters.SheetResistance,
                    ModelParameters.Rd,
                    ModelParameters.Rdc,
                    Parameters.DrainSquares,
                    ModelParameters.Rs,
                    ModelParameters.Rsc,
                    Parameters.SourceSquares,
                    out DrainResistance,
                    out SourceResistance
                );
            }
            if (DrainResistance > 0.0)
                this.DrainConductance = 1.0 / DrainResistance;
            else
                this.DrainConductance = 0.0;

            if (SourceResistance > 0.0)
                this.SourceConductance = 1.0 / SourceResistance;
            else
                this.SourceConductance = 0.0;

            this.Cgso = Param.BSIM3v32cgso;
            this.Cgdo = Param.BSIM3v32cgdo;

            Nvtm = ModelTemperature.Vtm * ModelParameters.JctEmissionCoeff;
            if (ModelParameters.AcmMod == 0)
            {
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
            }
            else /* ACM > 0 */
            {
                ACM.SaturationCurrents(
                    ModelParameters.AcmMod,
                    ModelParameters.Calcacm,
                    Parameters.Geo,
                    ModelParameters.Hdif,
                    ModelParameters.Wmlt,
                    Parameters.W,
                    ModelParameters.Xw,
                    ModelTemperature.JctTempSatCurDensity,
                    ModelTemperature.JctSidewallTempSatCurDensity,
                    Parameters.DrainArea,
                    Parameters.DrainPerimeter,
                    Parameters.SourceArea,
                    Parameters.SourcePerimeter,
                    out DrainSatCurrent,
                    out SourceSatCurrent
                );
            }

            if ((SourceSatCurrent > 0.0) && (ModelParameters.Ijth > 0.0))
            {
                this.Vjsm = Nvtm * Math.Log(ModelParameters.Ijth
                                / SourceSatCurrent + 1.0);
                /* Added revision dependent code */
                switch (ModelTemperature.IntVersion)
                {
                    case Version.BSIM3v32V324:
                    case Version.BSIM3v32V323:
                    case Version.BSIM3v32V322:
                        this.IsEvjsm =
                                SourceSatCurrent * Math.Exp(this.Vjsm / Nvtm);
                        break;
                    case Version.BSIM3v32V32:
                    default:
                        /* Do nothing */
                        break;
                }
            }

            if ((DrainSatCurrent > 0.0) && (ModelParameters.Ijth > 0.0))
            {
                this.Vjdm = Nvtm * Math.Log(ModelParameters.Ijth
                                / DrainSatCurrent + 1.0);
                /* Added revision dependent code */
                switch (ModelTemperature.IntVersion)
                {
                    case Version.BSIM3v32V324:
                    case Version.BSIM3v32V323:
                    case Version.BSIM3v32V322:
                        this.IsEvjdm =
                                DrainSatCurrent * Math.Exp(this.Vjdm / Nvtm);
                        break;
                    case Version.BSIM3v32V32:
                    default:
                        /* Do nothing */
                        break;
                }
            }
        }

        private bool BSIM3v2checkModel()
        {
            int Fatal_Flag = 0;
            List<string> words = new List<string>();

            if (ModelParameters.Version != "3.2.4" && ModelParameters.Version != "3.24" &&
                        ModelParameters.Version != "3.2.3" && ModelParameters.Version != "3.23" &&
                        ModelParameters.Version != "3.2.2" && ModelParameters.Version != "3.22" &&
                        ModelParameters.Version != "3.2" && ModelParameters.Version != "3.20")
            {
                SpiceSharpWarning.Warning(this, "Warning: This model supports BSIM3v3.2, BSIM3v3.2.2, BSIM3v3.2.3, BSIM3v3.2.4");
                SpiceSharpWarning.Warning(this, "You specified a wrong version number. Working now with BSIM3v3.2.4.");
                words.Add("Warning: This model supports BSIM3v3.2, BSIM3v3.2.2, BSIM3v3.2.3, BSIM3v3.2.4");
                words.Add("You specified a wrong version number. Working now with BSIM3v3.2.4.");
            }

            if (Param.BSIM3v32nlx < -Param.BSIM3v32leff)
            {
                words.Add("Fatal: Nlx = {0:g} is less than -Leff.".FormatString(
                            Param.BSIM3v32nlx));
                Fatal_Flag = 1;
            }

            if (ModelParameters.Tox <= 0.0)
            {
                words.Add("Fatal: Tox = {0:g} is not positive.".FormatString(ModelParameters.Tox));
                Fatal_Flag = 1;
            }

            if (ModelParameters.Toxm <= 0.0)
            {
                words.Add("Fatal: Toxm = {0:g} is not positive.".FormatString(ModelParameters.Toxm));
                Fatal_Flag = 1;
            }

            if (Param.BSIM3v32npeak <= 0.0)
            {
                words.Add("Fatal: Nch = {0:g} is not positive.".FormatString(Param.BSIM3v32npeak));
                Fatal_Flag = 1;
            }
            if (Param.BSIM3v32nsub <= 0.0)
            {
                words.Add("Fatal: Nsub = {0:g} is not positive.".FormatString(Param.BSIM3v32nsub));
                Fatal_Flag = 1;
            }
            if (Param.BSIM3v32ngate < 0.0)
            {
                words.Add("Fatal: Ngate = {0:g} is not positive.".FormatString(Param.BSIM3v32ngate));
                Fatal_Flag = 1;
            }
            if (Param.BSIM3v32ngate > 1.0e25)
            {
                words.Add("Fatal: Ngate = {0:g} is too high.".FormatString(Param.BSIM3v32ngate));
                Fatal_Flag = 1;
            }
            if (Param.BSIM3v32xj <= 0.0)
            {
                words.Add("Fatal: Xj = {0:g} is not positive.".FormatString(Param.BSIM3v32xj));
                Fatal_Flag = 1;
            }

            if (Param.BSIM3v32dvt1 < 0.0)
            {
                words.Add("Fatal: Dvt1 = {0:g} is negative.".FormatString(Param.BSIM3v32dvt1));
                Fatal_Flag = 1;
            }

            if (Param.BSIM3v32dvt1w < 0.0)
            {
                words.Add("Fatal: Dvt1w = {0:g} is negative.".FormatString(Param.BSIM3v32dvt1w));
                Fatal_Flag = 1;
            }

            if (Param.BSIM3v32w0 == -Param.BSIM3v32weff)
            {
                words.Add("Fatal: (W0 + Weff) = 0 causing divided-by-zero.");
                Fatal_Flag = 1;
            }

            if (Param.BSIM3v32dsub < 0.0)
            {
                words.Add("Fatal: Dsub = {0:g} is negative.".FormatString(Param.BSIM3v32dsub));
                Fatal_Flag = 1;
            }
            if (Param.BSIM3v32b1 == -Param.BSIM3v32weff)
            {
                words.Add("Fatal: (B1 + Weff) = 0 causing divided-by-zero.");
                Fatal_Flag = 1;
            }
            if (Param.BSIM3v32u0temp <= 0.0)
            {
                words.Add("Fatal: u0 at current temperature = {0:g} is not positive.".FormatString(Param.BSIM3v32u0temp));
                Fatal_Flag = 1;
            }

            /* Check delta parameter */
            if (Param.BSIM3v32delta < 0.0)
            {
                words.Add("Fatal: Delta = {0:g} is less than zero.".FormatString(Param.BSIM3v32delta));
                Fatal_Flag = 1;
            }

            if (Param.BSIM3v32vsattemp <= 0.0)
            {
                words.Add("Fatal: Vsat at current temperature = {0:g} is not positive.".FormatString(Param.BSIM3v32vsattemp));
                Fatal_Flag = 1;
            }
            /* Check Rout parameters */
            if (Param.BSIM3v32pclm <= 0.0)
            {
                words.Add("Fatal: Pclm = {0:g} is not positive.".FormatString(Param.BSIM3v32pclm));
                Fatal_Flag = 1;
            }

            if (Param.BSIM3v32drout < 0.0)
            {
                words.Add("Fatal: Drout = {0:g} is negative.".FormatString(Param.BSIM3v32drout));
                Fatal_Flag = 1;
            }

            if (Param.BSIM3v32pscbe2 <= 0.0)
            {
                words.Add("Warning: Pscbe2 = {0:g} is not positive.".FormatString(Param.BSIM3v32pscbe2));
            }

            /* ACM model */
            if (ModelParameters.AcmMod.Value == 0)
            {
                if (ModelParameters.UnitLengthSidewallJctCap > 0.0 ||
                      ModelParameters.UnitLengthGateSidewallJctCap > 0.0)
                {
                    if (Parameters.DrainPerimeter < Param.BSIM3v32weff)
                    {
                        words.Add("Warning: Pd = {0:g} is less than W.".FormatString(Parameters.DrainPerimeter));
                    }
                    if (Parameters.SourcePerimeter < Param.BSIM3v32weff)
                    {
                        words.Add("Warning: Ps = {0:g} is less than W.".FormatString(Parameters.SourcePerimeter));
                    }
                }
            }

            if ((ModelParameters.Calcacm.Value > 0) && (ModelParameters.AcmMod.Value != 12))
            {
                words.Add("Warning: CALCACM = %d is wrong. Set back to 0.".FormatString(ModelParameters.Calcacm.Value));
                ModelParameters.Calcacm = 0;
            }

            if (Param.BSIM3v32noff < 0.1)
            {
                words.Add("Warning: Noff = {0:g} is too small.".FormatString(Param.BSIM3v32noff));
            }
            if (Param.BSIM3v32noff > 4.0)
            {
                words.Add("Warning: Noff = {0:g} is too large.".FormatString(Param.BSIM3v32noff));
            }

            if (Param.BSIM3v32voffcv < -0.5)
            {
                words.Add("Warning: Voffcv = {0:g} is too small.".FormatString(Param.BSIM3v32voffcv));
            }
            if (Param.BSIM3v32voffcv > 0.5)
            {
                words.Add("Warning: Voffcv = {0:g} is too large.".FormatString(Param.BSIM3v32voffcv));
            }

            if (ModelParameters.Ijth < 0.0)
            {
                words.Add("Fatal: Ijth = {0:g} cannot be negative.".FormatString(ModelParameters.Ijth));
                Fatal_Flag = 1;
            }

            /* Check capacitance parameters */
            if (Param.BSIM3v32clc < 0.0)
            {
                words.Add("Fatal: Clc = {0:g} is negative.".FormatString(Param.BSIM3v32clc));
                Fatal_Flag = 1;
            }

            if (Param.BSIM3v32moin < 5.0)
            {
                words.Add("Warning: Moin = {0:g} is too small.".FormatString(Param.BSIM3v32moin));
            }
            if (Param.BSIM3v32moin > 25.0)
            {
                words.Add("Warning: Moin = {0:g} is too large.".FormatString(Param.BSIM3v32moin));
            }

            if (ModelParameters.CapMod.Value == 3)
            {
                if (Param.BSIM3v32acde < 0.4)
                {
                    words.Add("Warning:  Acde = {0:g} is too small.".FormatString(Param.BSIM3v32acde));
                }
                if (Param.BSIM3v32acde > 1.6)
                {
                    words.Add("Warning:  Acde = {0:g} is too large.".FormatString(Param.BSIM3v32acde));
                }
            }

            if (ModelParameters.ParamChk.Value == 1)
            {
                /* Check L and W parameters */
                if (Param.BSIM3v32leff <= 5.0e-8)
                {
                    words.Add("Warning: Leff = {0:g} may be too small.".FormatString(Param.BSIM3v32leff));
                }

                if (Param.BSIM3v32leffCV <= 5.0e-8)
                {
                    words.Add("Warning: Leff for CV = {0:g} may be too small.".FormatString(Param.BSIM3v32leffCV));
                }

                if (Param.BSIM3v32weff <= 1.0e-7)
                {
                    words.Add("Warning: Weff = {0:g} may be too small.".FormatString(Param.BSIM3v32weff));
                }

                if (Param.BSIM3v32weffCV <= 1.0e-7)
                {
                    words.Add("Warning: Weff for CV = {0:g} may be too small.".FormatString(Param.BSIM3v32weffCV));
                }

                /* Check threshold voltage parameters */
                if (Param.BSIM3v32nlx < 0.0)
                {
                    words.Add("Warning: Nlx = {0:g} is negative.".FormatString(Param.BSIM3v32nlx));
                }
                if (ModelParameters.Tox < 1.0e-9)
                {
                    words.Add("Warning: Tox = {0:g} is less than 10A.".FormatString(ModelParameters.Tox));
                }

                if (Param.BSIM3v32npeak <= 1.0e15)
                {
                    words.Add("Warning: Nch = {0:g} may be too small.".FormatString(Param.BSIM3v32npeak));
                }
                else if (Param.BSIM3v32npeak >= 1.0e21)
                {
                    words.Add("Warning: Nch = {0:g} may be too large.".FormatString(Param.BSIM3v32npeak));
                }

                if (Param.BSIM3v32nsub <= 1.0e14)
                {
                    words.Add("Warning: Nsub = {0:g} may be too small.".FormatString(Param.BSIM3v32nsub));
                }
                else if (Param.BSIM3v32nsub >= 1.0e21)
                {
                    words.Add("Warning: Nsub = {0:g} may be too large.".FormatString(Param.BSIM3v32nsub));
                }

                if ((Param.BSIM3v32ngate > 0.0) &&
                    (Param.BSIM3v32ngate <= 1.0e18))
                {
                    words.Add("Warning: Ngate = {0:g} is less than 1.E18cm^-3.".FormatString(Param.BSIM3v32ngate));
                }

                if (Param.BSIM3v32dvt0 < 0.0)
                {
                    words.Add("Warning: Dvt0 = {0:g} is negative.".FormatString(Param.BSIM3v32dvt0));
                }

                if (Math.Abs(1.0e-6 / (Param.BSIM3v32w0 + Param.BSIM3v32weff)) > 10.0)
                {
                    words.Add("Warning: (W0 + Weff) may be too small.");
                }

                /* Check subthreshold parameters */
                if (Param.BSIM3v32nfactor < 0.0)
                {
                    words.Add("Warning: Nfactor = {0:g} is negative.".FormatString(Param.BSIM3v32nfactor));
                }
                if (Param.BSIM3v32cdsc < 0.0)
                {
                    words.Add("Warning: Cdsc = {0:g} is negative.".FormatString(Param.BSIM3v32cdsc));
                }
                if (Param.BSIM3v32cdscd < 0.0)
                {
                    words.Add("Warning: Cdscd = {0:g} is negative.".FormatString(Param.BSIM3v32cdscd));
                }
                /* Check DIBL parameters */
                if (Param.BSIM3v32eta0 < 0.0)
                {
                    words.Add("Warning: Eta0 = {0:g} is negative.".FormatString(Param.BSIM3v32eta0));
                }

                /* Check Abulk parameters */
                if (Math.Abs(1.0e-6 / (Param.BSIM3v32b1 + Param.BSIM3v32weff)) > 10.0)
                {
                    words.Add("Warning: (B1 + Weff) may be too small.");
                }

                /* Check Saturation parameters */
                if (Param.BSIM3v32a2 < 0.01)
                {
                    words.Add("Warning: A2 = {0:g} is too small. Set to 0.01.".FormatString(Param.BSIM3v32a2));
                    Param.BSIM3v32a2 = 0.01;
                }
                else if (Param.BSIM3v32a2 > 1.0)
                {
                    words.Add("Warning: A2 = {0:g} is larger than 1. A2 is set to 1 and A1 is set to 0.".FormatString(Param.BSIM3v32a2));
                    Param.BSIM3v32a2 = 1.0;
                    Param.BSIM3v32a1 = 0.0;

                }

                if (Param.BSIM3v32rdsw < 0.0)
                {
                    words.Add("Warning: Rdsw = {0:g} is negative. Set to zero.".FormatString(Param.BSIM3v32rdsw));
                    Param.BSIM3v32rdsw = 0.0;
                    Param.BSIM3v32rds0 = 0.0;
                }
                else if ((Param.BSIM3v32rds0 > 0.0) && (Param.BSIM3v32rds0 < 0.001))
                {
                    words.Add("Warning: Rds at current temperature = {0:g} is less than 0.001 ohm. Set to zero.".FormatString(Param.BSIM3v32rds0));
                    Param.BSIM3v32rds0 = 0.0;
                }
                if (Param.BSIM3v32vsattemp < 1.0e3)
                {
                    words.Add("Warning: Vsat at current temperature = {0:g} may be too small.".FormatString(Param.BSIM3v32vsattemp));
                }

                if (Param.BSIM3v32pdibl1 < 0.0)
                {
                    words.Add("Warning: Pdibl1 = {0:g} is negative.".FormatString(Param.BSIM3v32pdibl1));
                }
                if (Param.BSIM3v32pdibl2 < 0.0)
                {
                    words.Add("Warning: Pdibl2 = {0:g} is negative.".FormatString(Param.BSIM3v32pdibl2));
                }
                /* Check overlap capacitance parameters */
                if (ModelParameters.Cgdo < 0.0)
                {
                    words.Add("Warning: cgdo = {0:g} is negative. Set to zero.".FormatString(ModelParameters.Cgdo));
                    ModelParameters.Cgdo = 0.0;
                }
                if (ModelParameters.Cgso < 0.0)
                {
                    words.Add("Warning: cgso = {0:g} is negative. Set to zero.".FormatString(ModelParameters.Cgso));
                    ModelParameters.Cgso = 0.0;
                }
                if (ModelParameters.Cgbo < 0.0)
                {
                    words.Add("Warning: cgbo = {0:g} is negative. Set to zero.".FormatString(ModelParameters.Cgbo));
                    ModelParameters.Cgbo = 0.0;
                }

            }

            if (words.Count > 0)
            {
                string path = ModelParameters.CheckPath;
                if (string.IsNullOrWhiteSpace(path))
                {
                    path = "b3v32check.log";
                }
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
