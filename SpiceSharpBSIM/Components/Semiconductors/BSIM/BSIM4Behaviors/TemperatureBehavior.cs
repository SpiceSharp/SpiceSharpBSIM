using SpiceSharp;
using SpiceSharp.Attributes;
using SpiceSharp.Behaviors;
using SpiceSharp.Components;
using SpiceSharp.ParameterSets;
using SpiceSharp.Simulations;
using System;
using System.Collections.Generic;
using System.IO;
using System.Text;

namespace SpiceSharpBSIM.Components.Semiconductors.BSIM.BSIM4Behaviors
{
    /// <summary>
    /// A temperature behavior for a <see cref="BSIM4"/>.
    /// </summary>
    [BehaviorFor(typeof(BSIM4)), AddBehaviorIfNo(typeof(ITemperatureBehavior))]
    public partial class TemperatureBehavior : Behavior, IParameterized<BaseParameters>, ITemperatureBehavior
    {
        private readonly ITemperatureSimulationState _temperature;

        public const double EPS0 = 8.85418e-12;
        public const double MAX_EXP = 5.834617425e14;
        public const double MIN_EXP = 1.713908431e-15;
        public const double EXP_THRESHOLD = 34.0;
        public const double DELTA = 1.0E-9;
        public const double KboQ = 8.617087e-5;

        /// <inheritdoc />
        public BaseParameters Parameters { get; }

        protected ModelTemperatureBehavior ModelTemperature { get; }
        protected ModelParameters ModelParameters { get; }
        protected SizeDependentProperties Param { get; private set; }

        protected double _u0temp, _vsattemp, _vth0, _eta0, _k2, _vfb, _vtfbphi1, _vtfbphi2, _vbsc, _k2ox, _vfbzb,
            _cgso, _cgdo, _grbdb, _grbpb, _grbps, _grbsb, _grbpd, _grgeltd, _pseff, _pdeff, _aseff, _adeff,
            _sourceConductance, _drainConductance, _xExpBVS, _vjsmFwd, _iVjsmFwd, _sslpFwd, _vjsmRev, _iVjsmRev,
            _xExpBVD, _vjdmFwd, _iVjdmFwd, _dslpFwd, _vjdmRev, _iVjdmRev, _dslpRev, _sjctTempRevSatCur, _djctTempRevSatCur,
            _sswTempRevSatCur, _dswTempRevSatCur, _sswgTempRevSatCur, _dswgTempRevSatCur, _toxp, _coxp,
            _sslpRev;

        /// <summary>
        /// Creates a new <see cref="TemperatureBehavior"/>.
        /// </summary>
        /// <param name="context"></param>
        public TemperatureBehavior(ComponentBindingContext context)
            : base(context)
        {
            _temperature = context.GetState<ITemperatureSimulationState>();
            Parameters = context.GetParameterSet<BaseParameters>();
            ModelTemperature = context.ModelBehaviors.GetValue<ModelTemperatureBehavior>();
            ModelParameters = context.ModelBehaviors.GetParameterSet<ModelParameters>();
            Setup();
        }

        private void Setup()
        {
            // Setup
            if (!Parameters.Rbdb.Given)
                Parameters.Rbdb = new GivenParameter<double>(ModelParameters.Rbdb, false);
            if (!Parameters.Rbsb.Given)
                Parameters.Rbsb = new GivenParameter<double>(ModelParameters.Rbsb, false);
            if (!Parameters.Rbpb.Given)
                Parameters.Rbpb = new GivenParameter<double>(ModelParameters.Rbpb, false);
            if (!Parameters.Rbps.Given)
                Parameters.Rbps = new GivenParameter<double>(ModelParameters.Rbps, false);
            if (!Parameters.Rbpd.Given)
                Parameters.Rbpd = new GivenParameter<double>(ModelParameters.Rbpd, false);
            if (!Parameters.Xgw.Given)
                Parameters.Xgw = new GivenParameter<double>(ModelParameters.Xgw, false);
            if (!Parameters.Ngcon.Given)
                Parameters.Ngcon = new GivenParameter<double>(ModelParameters.Ngcon, false);
            if (!Parameters.RbodyMod.Given)
                Parameters.RbodyMod = new GivenParameter<int>(ModelParameters.RbodyMod, false);
            else if ((Parameters.RbodyMod.Value != 0) && (Parameters.RbodyMod.Value != 1) && (Parameters.RbodyMod.Value != 2))
            {
                Parameters.RbodyMod = new GivenParameter<int>(ModelParameters.RbodyMod, false);
                SpiceSharpWarning.Warning(this, "Warning: rbodyMod has been set to its global value {0}.".FormatString(ModelParameters.RbodyMod));
            }
            if (!Parameters.RgateMod.Given)
                Parameters.RgateMod = new GivenParameter<int>(ModelParameters.RgateMod, false);
            else if ((Parameters.RgateMod.Value != 0) && (Parameters.RgateMod.Value != 1)
                && (Parameters.RgateMod.Value != 2) && (Parameters.RgateMod.Value != 3))
            {
                Parameters.RgateMod = new GivenParameter<int>(ModelParameters.RgateMod, false);
                SpiceSharpWarning.Warning(this, "Warning: rgateMod has been set to its global value {0}.".FormatString(ModelParameters.RgateMod));
            }
            if (!Parameters.GeoMod.Given)
                Parameters.GeoMod = new GivenParameter<int>(ModelParameters.GeoMod, false);
            if (!Parameters.RgeoMod.Given)
                Parameters.RgeoMod = new GivenParameter<int>(ModelParameters.RgeoMod, false);
            else if ((Parameters.RgeoMod.Value != 0) && (Parameters.RgeoMod.Value != 1))
            {
                Parameters.RgeoMod = new GivenParameter<int>(ModelParameters.RgeoMod, false);
                SpiceSharpWarning.Warning(this, "Warning: rgeoMod has been set to its global value {0}.".FormatString(ModelParameters.RgeoMod));
            }
            if (!Parameters.TrnqsMod.Given)
                Parameters.TrnqsMod = new GivenParameter<int>(ModelParameters.TrnqsMod, false);
            else if ((Parameters.TrnqsMod.Value != 0) && (Parameters.TrnqsMod.Value != 1))
            {
                Parameters.TrnqsMod = new GivenParameter<int>(ModelParameters.TrnqsMod, false);
                SpiceSharpWarning.Warning(this, "Warning: trnqsMod has been set to its global value {0}.".FormatString(ModelParameters.TrnqsMod));
            }
            if (!Parameters.AcnqsMod.Given)
                Parameters.AcnqsMod = new GivenParameter<int>(ModelParameters.AcnqsMod, false);
            else if ((Parameters.AcnqsMod.Value != 0) && (Parameters.AcnqsMod.Value != 1))
            {
                Parameters.AcnqsMod = new GivenParameter<int>(ModelParameters.AcnqsMod, false);
                SpiceSharpWarning.Warning(this, "Warning: acnqsMod has been set to its global value {0}.".FormatString(ModelParameters.AcnqsMod));
            }
            if (!Parameters.Sd.Given)
                Parameters.Sd = new GivenParameter<double>(2 * ModelParameters.Dmcg, false);
        }

        /// <inheritdoc />
        void ITemperatureBehavior.Temperature()
        {
            double tmp, tmp1, tmp2, tmp3, Eg0, ni, epssub;
            double T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, Lnew = 0.0, Wnew;
            double delTemp, TRatio, Inv_L, Inv_W, Inv_LW, Vtm0, Tnom;
            double dumPs, dumPd, dumAs, dumAd, PowWeffWr;
            double DMCGeff, DMCIeff, DMDGeff;
            double Nvtms, Nvtmd, SourceSatCurrent, DrainSatCurrent;
            double T10, T11;
            double Inv_saref, Inv_sbref, Inv_sa, Inv_sb, rho, Ldrn, dvth0_lod;
            double W_tmp, Inv_ODeff, OD_offset, dk2_lod, deta0_lod;
            double lnl, lnw, lnnf, rbpbx, rbpby, rbsbx, rbsby, rbdbx, rbdby, bodymode;
            double kvsat, wlod, sceff, Wdrn;
            double V0, lt1, ltw, Theta0, Delt_vth, Vth_NarrowW, Lpe_Vb, Vth;
            double n, n0, Vgsteff, Vgs_eff, niter, toxpf, toxpi, Tcen, toxe, epsrox, vddeot;
            double vtfbphi2eot, phieot, TempRatioeot, Vtm0eot, Vtmeot, vbieot;

            TRatio = ModelTemperature.TRatio;
            delTemp = ModelTemperature.DelTemp;
            Tnom = ModelParameters.Tnom;
            Vtm0 = ModelTemperature.Vtm0;
            ni = ModelTemperature.Ni;
            epssub = ModelTemperature.Epssub;
            epsrox = ModelTemperature.Epsrox;
            toxe = ModelTemperature.Toxe;
            Eg0 = ModelTemperature.Eg0;

            /* Temperature */

            /* stress effect */
            Ldrn = Parameters.L;
            Wdrn = Parameters.W / Parameters.Nf;

            var key = Tuple.Create(Parameters.W.Value, Parameters.L.Value, Parameters.Nf.Value);
            if (ModelTemperature.SizeDependentProperties.TryGetValue(key, out var param))
                Param = param;
            else
            {
                Param = new SizeDependentProperties();
                ModelTemperature.SizeDependentProperties.Add(key, Param);

                Param.Length = Parameters.L;
                Param.Width = Parameters.W;
                Param.NFinger = Parameters.Nf;
                Lnew = Parameters.L + ModelParameters.Xl;
                Wnew = Parameters.W / Parameters.Nf + ModelParameters.Xw;

                T0 = Math.Pow(Lnew, ModelParameters.Lln);
                T1 = Math.Pow(Wnew, ModelParameters.Lwn);
                tmp1 = ModelParameters.Ll / T0 + ModelParameters.Lw / T1
                               + ModelParameters.Lwl / (T0 * T1);
                Param.BSIM4dl = ModelParameters.Lint + tmp1;
                tmp2 = ModelParameters.Llc / T0 + ModelParameters.Lwc / T1
                     + ModelParameters.Lwlc / (T0 * T1);
                Param.BSIM4dlc = ModelParameters.Dlc + tmp2;

                T2 = Math.Pow(Lnew, ModelParameters.Wln);
                T3 = Math.Pow(Wnew, ModelParameters.Wwn);
                tmp1 = ModelParameters.Wl / T2 + ModelParameters.Ww / T3
                               + ModelParameters.Wwl / (T2 * T3);
                Param.BSIM4dw = ModelParameters.Wint + tmp1;
                tmp2 = ModelParameters.Wlc / T2 + ModelParameters.Wwc / T3
                     + ModelParameters.Wwlc / (T2 * T3);
                Param.BSIM4dwc = ModelParameters.Dwc + tmp2;
                Param.BSIM4dwj = ModelParameters.Dwj + tmp2;

                Param.BSIM4leff = Lnew - 2.0 * Param.BSIM4dl;
                if (Param.BSIM4leff <= 0.0)
                {
                    throw new SpiceSharpException("BSIM4: mosfet {0}, model {1}: Effective channel length <= 0".FormatString(Name, ModelTemperature.Name));
                }

                Param.BSIM4weff = Wnew - 2.0 * Param.BSIM4dw;
                if (Param.BSIM4weff <= 0.0)
                {
                    throw new SpiceSharpException("BSIM4: mosfet {0}, model {1}: Effective channel width <= 0".FormatString(Name, ModelTemperature.Name));
                }

                Param.BSIM4leffCV = Lnew - 2.0 * Param.BSIM4dlc;
                if (Param.BSIM4leffCV <= 0.0)
                {
                    throw new SpiceSharpException("BSIM4: mosfet {0}, model {1}: Effective channel length for C-V <= 0".FormatString(Name, ModelTemperature.Name));
                }

                Param.BSIM4weffCV = Wnew - 2.0 * Param.BSIM4dwc;
                if (Param.BSIM4weffCV <= 0.0)
                {
                    throw new SpiceSharpException("BSIM4: mosfet {0}, model {1}: Effective channel width for C-V <= 0".FormatString(Name, ModelTemperature.Name));
                }

                Param.BSIM4weffCJ = Wnew - 2.0 * Param.BSIM4dwj;
                if (Param.BSIM4weffCJ <= 0.0)
                {
                    throw new SpiceSharpException("BSIM4: mosfet {0}, model {1}: Effective channel width for S/D junctions <= 0".FormatString(Name, ModelTemperature.Name));
                }


                if (ModelParameters.BinUnit == 1)
                {
                    Inv_L = 1.0e-6 / Param.BSIM4leff;
                    Inv_W = 1.0e-6 / Param.BSIM4weff;
                    Inv_LW = 1.0e-12 / (Param.BSIM4leff
                           * Param.BSIM4weff);
                }
                else
                {
                    Inv_L = 1.0 / Param.BSIM4leff;
                    Inv_W = 1.0 / Param.BSIM4weff;
                    Inv_LW = 1.0 / (Param.BSIM4leff
                           * Param.BSIM4weff);
                }
                Param.BSIM4cdsc = ModelParameters.Cdsc
                                  + ModelParameters.Lcdsc * Inv_L
                                  + ModelParameters.Wcdsc * Inv_W
                                  + ModelParameters.Pcdsc * Inv_LW;
                Param.BSIM4cdscb = ModelParameters.Cdscb
                                   + ModelParameters.Lcdscb * Inv_L
                                   + ModelParameters.Wcdscb * Inv_W
                                   + ModelParameters.Pcdscb * Inv_LW;

                Param.BSIM4cdscd = ModelParameters.Cdscd
                               + ModelParameters.Lcdscd * Inv_L
                               + ModelParameters.Wcdscd * Inv_W
                               + ModelParameters.Pcdscd * Inv_LW;

                Param.BSIM4cit = ModelParameters.Cit
                                 + ModelParameters.Lcit * Inv_L
                                 + ModelParameters.Wcit * Inv_W
                                 + ModelParameters.Pcit * Inv_LW;
                Param.BSIM4nfactor = ModelParameters.Nfactor
                                     + ModelParameters.Lnfactor * Inv_L
                                     + ModelParameters.Wnfactor * Inv_W
                                     + ModelParameters.Pnfactor * Inv_LW;
                Param.BSIM4tnfactor = ModelParameters.Tnfactor                        /* v4.7 */
                                     + ModelParameters.Ltnfactor * Inv_L
                                     + ModelParameters.Wtnfactor * Inv_W
                                     + ModelParameters.Ptnfactor * Inv_LW;
                Param.BSIM4xj = ModelParameters.Xj
                                + ModelParameters.Lxj * Inv_L
                                + ModelParameters.Wxj * Inv_W
                                + ModelParameters.Pxj * Inv_LW;
                Param.BSIM4vsat = ModelParameters.Vsat
                                  + ModelParameters.Lvsat * Inv_L
                                  + ModelParameters.Wvsat * Inv_W
                                  + ModelParameters.Pvsat * Inv_LW;
                Param.BSIM4at = ModelParameters.At
                                + ModelParameters.Lat * Inv_L
                                + ModelParameters.Wat * Inv_W
                                + ModelParameters.Pat * Inv_LW;
                Param.BSIM4a0 = ModelParameters.A0
                                + ModelParameters.La0 * Inv_L
                                + ModelParameters.Wa0 * Inv_W
                                + ModelParameters.Pa0 * Inv_LW;

                Param.BSIM4ags = ModelParameters.Ags
                                + ModelParameters.Lags * Inv_L
                                + ModelParameters.Wags * Inv_W
                                + ModelParameters.Pags * Inv_LW;

                Param.BSIM4a1 = ModelParameters.A1
                                + ModelParameters.La1 * Inv_L
                                + ModelParameters.Wa1 * Inv_W
                                + ModelParameters.Pa1 * Inv_LW;
                Param.BSIM4a2 = ModelParameters.A2
                                + ModelParameters.La2 * Inv_L
                                + ModelParameters.Wa2 * Inv_W
                                + ModelParameters.Pa2 * Inv_LW;
                Param.BSIM4keta = ModelParameters.Keta
                                  + ModelParameters.Lketa * Inv_L
                                  + ModelParameters.Wketa * Inv_W
                                  + ModelParameters.Pketa * Inv_LW;
                Param.BSIM4nsub = ModelParameters.Nsub
                                  + ModelParameters.Lnsub * Inv_L
                                  + ModelParameters.Wnsub * Inv_W
                                  + ModelParameters.Pnsub * Inv_LW;
                Param.BSIM4ndep = ModelParameters.Ndep
                                  + ModelParameters.Lndep * Inv_L
                                  + ModelParameters.Wndep * Inv_W
                                  + ModelParameters.Pndep * Inv_LW;
                Param.BSIM4nsd = ModelParameters.Nsd
                                 + ModelParameters.Lnsd * Inv_L
                                 + ModelParameters.Wnsd * Inv_W
                                 + ModelParameters.Pnsd * Inv_LW;
                Param.BSIM4phin = ModelParameters.Phin
                                  + ModelParameters.Lphin * Inv_L
                                  + ModelParameters.Wphin * Inv_W
                                  + ModelParameters.Pphin * Inv_LW;
                Param.BSIM4ngate = ModelParameters.Ngate
                                   + ModelParameters.Lngate * Inv_L
                                   + ModelParameters.Wngate * Inv_W
                                   + ModelParameters.Pngate * Inv_LW;
                Param.BSIM4gamma1 = ModelParameters.Gamma1
                                    + ModelParameters.Lgamma1 * Inv_L
                                    + ModelParameters.Wgamma1 * Inv_W
                                    + ModelParameters.Pgamma1 * Inv_LW;
                Param.BSIM4gamma2 = ModelParameters.Gamma2
                                    + ModelParameters.Lgamma2 * Inv_L
                                    + ModelParameters.Wgamma2 * Inv_W
                                    + ModelParameters.Pgamma2 * Inv_LW;
                Param.BSIM4vbx = ModelParameters.Vbx
                                 + ModelParameters.Lvbx * Inv_L
                                 + ModelParameters.Wvbx * Inv_W
                                 + ModelParameters.Pvbx * Inv_LW;
                Param.BSIM4vbm = ModelParameters.Vbm
                                 + ModelParameters.Lvbm * Inv_L
                                 + ModelParameters.Wvbm * Inv_W
                                 + ModelParameters.Pvbm * Inv_LW;
                Param.BSIM4xt = ModelParameters.Xt
                                 + ModelParameters.Lxt * Inv_L
                                 + ModelParameters.Wxt * Inv_W
                                 + ModelParameters.Pxt * Inv_LW;
                Param.BSIM4vfb = ModelParameters.Vfb
                                 + ModelParameters.Lvfb * Inv_L
                                 + ModelParameters.Wvfb * Inv_W
                                 + ModelParameters.Pvfb * Inv_LW;
                Param.BSIM4k1 = ModelParameters.K1
                                + ModelParameters.Lk1 * Inv_L
                                + ModelParameters.Wk1 * Inv_W
                                + ModelParameters.Pk1 * Inv_LW;
                Param.BSIM4kt1 = ModelParameters.Kt1
                                 + ModelParameters.Lkt1 * Inv_L
                                 + ModelParameters.Wkt1 * Inv_W
                                 + ModelParameters.Pkt1 * Inv_LW;
                Param.BSIM4kt1l = ModelParameters.Kt1l
                                  + ModelParameters.Lkt1l * Inv_L
                                  + ModelParameters.Wkt1l * Inv_W
                                  + ModelParameters.Pkt1l * Inv_LW;
                Param.BSIM4k2 = ModelParameters.K2
                                + ModelParameters.Lk2 * Inv_L
                                + ModelParameters.Wk2 * Inv_W
                                + ModelParameters.Pk2 * Inv_LW;
                Param.BSIM4kt2 = ModelParameters.Kt2
                                 + ModelParameters.Lkt2 * Inv_L
                                 + ModelParameters.Wkt2 * Inv_W
                                 + ModelParameters.Pkt2 * Inv_LW;
                Param.BSIM4k3 = ModelParameters.K3
                                + ModelParameters.Lk3 * Inv_L
                                + ModelParameters.Wk3 * Inv_W
                                + ModelParameters.Pk3 * Inv_LW;
                Param.BSIM4k3b = ModelParameters.K3b
                                 + ModelParameters.Lk3b * Inv_L
                                 + ModelParameters.Wk3b * Inv_W
                                 + ModelParameters.Pk3b * Inv_LW;
                Param.BSIM4w0 = ModelParameters.W0
                                + ModelParameters.Lw0 * Inv_L
                                + ModelParameters.Ww0 * Inv_W
                                + ModelParameters.Pw0 * Inv_LW;
                Param.BSIM4lpe0 = ModelParameters.Lpe0
                                  + ModelParameters.Llpe0 * Inv_L
                                   + ModelParameters.Wlpe0 * Inv_W
                                  + ModelParameters.Plpe0 * Inv_LW;
                Param.BSIM4lpeb = ModelParameters.Lpeb
                                  + ModelParameters.Llpeb * Inv_L
                                  + ModelParameters.Wlpeb * Inv_W
                                  + ModelParameters.Plpeb * Inv_LW;
                Param.BSIM4dvtp0 = ModelParameters.Dvtp0
                                   + ModelParameters.Ldvtp0 * Inv_L
                                   + ModelParameters.Wdvtp0 * Inv_W
                                   + ModelParameters.Pdvtp0 * Inv_LW;
                Param.BSIM4dvtp1 = ModelParameters.Dvtp1
                                   + ModelParameters.Ldvtp1 * Inv_L
                                   + ModelParameters.Wdvtp1 * Inv_W
                                   + ModelParameters.Pdvtp1 * Inv_LW;
                Param.BSIM4dvtp2 = ModelParameters.Dvtp2                 /* v4.7  */
                                   + ModelParameters.Ldvtp2 * Inv_L
                                   + ModelParameters.Wdvtp2 * Inv_W
                                   + ModelParameters.Pdvtp2 * Inv_LW;
                Param.BSIM4dvtp3 = ModelParameters.Dvtp3                 /* v4.7  */
                                   + ModelParameters.Ldvtp3 * Inv_L
                                   + ModelParameters.Wdvtp3 * Inv_W
                                   + ModelParameters.Pdvtp3 * Inv_LW;
                Param.BSIM4dvtp4 = ModelParameters.Dvtp4                 /* v4.7  */
                                   + ModelParameters.Ldvtp4 * Inv_L
                                   + ModelParameters.Wdvtp4 * Inv_W
                                   + ModelParameters.Pdvtp4 * Inv_LW;
                Param.BSIM4dvtp5 = ModelParameters.Dvtp5                 /* v4.7  */
                                   + ModelParameters.Ldvtp5 * Inv_L
                                   + ModelParameters.Wdvtp5 * Inv_W
                                   + ModelParameters.Pdvtp5 * Inv_LW;
                Param.BSIM4dvt0 = ModelParameters.Dvt0
                                  + ModelParameters.Ldvt0 * Inv_L
                                  + ModelParameters.Wdvt0 * Inv_W
                                  + ModelParameters.Pdvt0 * Inv_LW;
                Param.BSIM4dvt1 = ModelParameters.Dvt1
                                  + ModelParameters.Ldvt1 * Inv_L
                                  + ModelParameters.Wdvt1 * Inv_W
                                  + ModelParameters.Pdvt1 * Inv_LW;
                Param.BSIM4dvt2 = ModelParameters.Dvt2
                                  + ModelParameters.Ldvt2 * Inv_L
                                  + ModelParameters.Wdvt2 * Inv_W
                                  + ModelParameters.Pdvt2 * Inv_LW;
                Param.BSIM4dvt0w = ModelParameters.Dvt0w
                                  + ModelParameters.Ldvt0w * Inv_L
                                  + ModelParameters.Wdvt0w * Inv_W
                                  + ModelParameters.Pdvt0w * Inv_LW;
                Param.BSIM4dvt1w = ModelParameters.Dvt1w
                                  + ModelParameters.Ldvt1w * Inv_L
                                  + ModelParameters.Wdvt1w * Inv_W
                                  + ModelParameters.Pdvt1w * Inv_LW;
                Param.BSIM4dvt2w = ModelParameters.Dvt2w
                                  + ModelParameters.Ldvt2w * Inv_L
                                  + ModelParameters.Wdvt2w * Inv_W
                                  + ModelParameters.Pdvt2w * Inv_LW;
                Param.BSIM4drout = ModelParameters.Drout
                                   + ModelParameters.Ldrout * Inv_L
                                   + ModelParameters.Wdrout * Inv_W
                                   + ModelParameters.Pdrout * Inv_LW;
                Param.BSIM4dsub = ModelParameters.Dsub
                                  + ModelParameters.Ldsub * Inv_L
                                  + ModelParameters.Wdsub * Inv_W
                                  + ModelParameters.Pdsub * Inv_LW;
                Param.BSIM4vth0 = ModelParameters.Vth0
                                  + ModelParameters.Lvth0 * Inv_L
                                  + ModelParameters.Wvth0 * Inv_W
                                  + ModelParameters.Pvth0 * Inv_LW;
                Param.BSIM4ua = ModelParameters.Ua
                                + ModelParameters.Lua * Inv_L
                                + ModelParameters.Wua * Inv_W
                                + ModelParameters.Pua * Inv_LW;
                Param.BSIM4ua1 = ModelParameters.Ua1
                                 + ModelParameters.Lua1 * Inv_L
                                 + ModelParameters.Wua1 * Inv_W
                                 + ModelParameters.Pua1 * Inv_LW;
                Param.BSIM4ub = ModelParameters.Ub
                                + ModelParameters.Lub * Inv_L
                                + ModelParameters.Wub * Inv_W
                                + ModelParameters.Pub * Inv_LW;
                Param.BSIM4ub1 = ModelParameters.Ub1
                                 + ModelParameters.Lub1 * Inv_L
                                 + ModelParameters.Wub1 * Inv_W
                                 + ModelParameters.Pub1 * Inv_LW;
                Param.BSIM4uc = ModelParameters.Uc
                                + ModelParameters.Luc * Inv_L
                                + ModelParameters.Wuc * Inv_W
                                + ModelParameters.Puc * Inv_LW;
                Param.BSIM4uc1 = ModelParameters.Uc1
                                 + ModelParameters.Luc1 * Inv_L
                                 + ModelParameters.Wuc1 * Inv_W
                                 + ModelParameters.Puc1 * Inv_LW;
                Param.BSIM4ud = ModelParameters.Ud
                                + ModelParameters.Lud * Inv_L
                                + ModelParameters.Wud * Inv_W
                                + ModelParameters.Pud * Inv_LW;
                Param.BSIM4ud1 = ModelParameters.Ud1
                                + ModelParameters.Lud1 * Inv_L
                                + ModelParameters.Wud1 * Inv_W
                                + ModelParameters.Pud1 * Inv_LW;
                Param.BSIM4up = ModelParameters.Up
                                + ModelParameters.Lup * Inv_L
                                + ModelParameters.Wup * Inv_W
                                + ModelParameters.Pup * Inv_LW;
                Param.BSIM4lp = ModelParameters.Lp
                                + ModelParameters.Llp * Inv_L
                                + ModelParameters.Wlp * Inv_W
                                + ModelParameters.Plp * Inv_LW;
                Param.BSIM4eu = ModelParameters.Eu
                                + ModelParameters.Leu * Inv_L
                                + ModelParameters.Weu * Inv_W
                                + ModelParameters.Peu * Inv_LW;
                Param.BSIM4u0 = ModelParameters.U0
                                + ModelParameters.Lu0 * Inv_L
                                + ModelParameters.Wu0 * Inv_W
                                + ModelParameters.Pu0 * Inv_LW;
                Param.BSIM4ute = ModelParameters.Ute
                                 + ModelParameters.Lute * Inv_L
                                 + ModelParameters.Wute * Inv_W
                                 + ModelParameters.Pute * Inv_LW;
                /*high k mobility*/
                Param.BSIM4ucs = ModelParameters.Ucs
                                 + ModelParameters.Lucs * Inv_L
                                 + ModelParameters.Wucs * Inv_W
                                 + ModelParameters.Pucs * Inv_LW;
                Param.BSIM4ucste = ModelParameters.Ucste
                         + ModelParameters.Lucste * Inv_L
                                 + ModelParameters.Wucste * Inv_W
                                 + ModelParameters.Pucste * Inv_LW;

                Param.BSIM4voff = ModelParameters.Voff
                                  + ModelParameters.Lvoff * Inv_L
                                  + ModelParameters.Wvoff * Inv_W
                                  + ModelParameters.Pvoff * Inv_LW;
                Param.BSIM4tvoff = ModelParameters.Tvoff
                                  + ModelParameters.Ltvoff * Inv_L
                                  + ModelParameters.Wtvoff * Inv_W
                                  + ModelParameters.Ptvoff * Inv_LW;
                Param.BSIM4minv = ModelParameters.Minv
                                  + ModelParameters.Lminv * Inv_L
                                  + ModelParameters.Wminv * Inv_W
                                  + ModelParameters.Pminv * Inv_LW;
                Param.BSIM4minvcv = ModelParameters.Minvcv
                                  + ModelParameters.Lminvcv * Inv_L
                                  + ModelParameters.Wminvcv * Inv_W
                                  + ModelParameters.Pminvcv * Inv_LW;
                Param.BSIM4fprout = ModelParameters.Fprout
                                   + ModelParameters.Lfprout * Inv_L
                                   + ModelParameters.Wfprout * Inv_W
                                   + ModelParameters.Pfprout * Inv_LW;
                Param.BSIM4pdits = ModelParameters.Pdits
                                   + ModelParameters.Lpdits * Inv_L
                                   + ModelParameters.Wpdits * Inv_W
                                   + ModelParameters.Ppdits * Inv_LW;
                Param.BSIM4pditsd = ModelParameters.Pditsd
                                    + ModelParameters.Lpditsd * Inv_L
                                    + ModelParameters.Wpditsd * Inv_W
                                    + ModelParameters.Ppditsd * Inv_LW;
                Param.BSIM4delta = ModelParameters.Delta
                                   + ModelParameters.Ldelta * Inv_L
                                   + ModelParameters.Wdelta * Inv_W
                                   + ModelParameters.Pdelta * Inv_LW;
                Param.BSIM4rdsw = ModelParameters.Rdsw
                                  + ModelParameters.Lrdsw * Inv_L
                                  + ModelParameters.Wrdsw * Inv_W
                                  + ModelParameters.Prdsw * Inv_LW;
                Param.BSIM4rdw = ModelParameters.Rdw
                                  + ModelParameters.Lrdw * Inv_L
                                  + ModelParameters.Wrdw * Inv_W
                                  + ModelParameters.Prdw * Inv_LW;
                Param.BSIM4rsw = ModelParameters.Rsw
                                  + ModelParameters.Lrsw * Inv_L
                                  + ModelParameters.Wrsw * Inv_W
                                  + ModelParameters.Prsw * Inv_LW;
                Param.BSIM4prwg = ModelParameters.Prwg
                                  + ModelParameters.Lprwg * Inv_L
                                  + ModelParameters.Wprwg * Inv_W
                                  + ModelParameters.Pprwg * Inv_LW;
                Param.BSIM4prwb = ModelParameters.Prwb
                                  + ModelParameters.Lprwb * Inv_L
                                  + ModelParameters.Wprwb * Inv_W
                                  + ModelParameters.Pprwb * Inv_LW;
                Param.BSIM4prt = ModelParameters.Prt
                                  + ModelParameters.Lprt * Inv_L
                                  + ModelParameters.Wprt * Inv_W
                                  + ModelParameters.Pprt * Inv_LW;
                Param.BSIM4eta0 = ModelParameters.Eta0
                                  + ModelParameters.Leta0 * Inv_L
                                  + ModelParameters.Weta0 * Inv_W
                                  + ModelParameters.Peta0 * Inv_LW;
                Param.BSIM4teta0 = ModelParameters.Teta0                 /* v4.7  */
                                  + ModelParameters.Lteta0 * Inv_L
                                  + ModelParameters.Wteta0 * Inv_W
                                  + ModelParameters.Pteta0 * Inv_LW;
                Param.BSIM4tvoffcv = ModelParameters.Tvoffcv /* v4.8.0  */
                                  + ModelParameters.Ltvoffcv * Inv_L
                                  + ModelParameters.Wtvoffcv * Inv_W
                                  + ModelParameters.Ptvoffcv * Inv_LW;
                Param.BSIM4etab = ModelParameters.Etab
                                  + ModelParameters.Letab * Inv_L
                                  + ModelParameters.Wetab * Inv_W
                                  + ModelParameters.Petab * Inv_LW;
                Param.BSIM4pclm = ModelParameters.Pclm
                                  + ModelParameters.Lpclm * Inv_L
                                  + ModelParameters.Wpclm * Inv_W
                                  + ModelParameters.Ppclm * Inv_LW;
                Param.BSIM4pdibl1 = ModelParameters.Pdibl1
                                    + ModelParameters.Lpdibl1 * Inv_L
                                    + ModelParameters.Wpdibl1 * Inv_W
                                    + ModelParameters.Ppdibl1 * Inv_LW;
                Param.BSIM4pdibl2 = ModelParameters.Pdibl2
                                    + ModelParameters.Lpdibl2 * Inv_L
                                    + ModelParameters.Wpdibl2 * Inv_W
                                    + ModelParameters.Ppdibl2 * Inv_LW;
                Param.BSIM4pdiblb = ModelParameters.Pdiblb
                                    + ModelParameters.Lpdiblb * Inv_L
                                    + ModelParameters.Wpdiblb * Inv_W
                                    + ModelParameters.Ppdiblb * Inv_LW;
                Param.BSIM4pscbe1 = ModelParameters.Pscbe1
                                    + ModelParameters.Lpscbe1 * Inv_L
                                    + ModelParameters.Wpscbe1 * Inv_W
                                    + ModelParameters.Ppscbe1 * Inv_LW;
                Param.BSIM4pscbe2 = ModelParameters.Pscbe2
                                    + ModelParameters.Lpscbe2 * Inv_L
                                    + ModelParameters.Wpscbe2 * Inv_W
                                    + ModelParameters.Ppscbe2 * Inv_LW;
                Param.BSIM4pvag = ModelParameters.Pvag
                                  + ModelParameters.Lpvag * Inv_L
                                  + ModelParameters.Wpvag * Inv_W
                                  + ModelParameters.Ppvag * Inv_LW;
                Param.BSIM4wr = ModelParameters.Wr
                                + ModelParameters.Lwr * Inv_L
                                + ModelParameters.Wwr * Inv_W
                                + ModelParameters.Pwr * Inv_LW;
                Param.BSIM4dwg = ModelParameters.Dwg
                                 + ModelParameters.Ldwg * Inv_L
                                 + ModelParameters.Wdwg * Inv_W
                                 + ModelParameters.Pdwg * Inv_LW;
                Param.BSIM4dwb = ModelParameters.Dwb
                                 + ModelParameters.Ldwb * Inv_L
                                 + ModelParameters.Wdwb * Inv_W
                                 + ModelParameters.Pdwb * Inv_LW;
                Param.BSIM4b0 = ModelParameters.B0
                                + ModelParameters.Lb0 * Inv_L
                                + ModelParameters.Wb0 * Inv_W
                                + ModelParameters.Pb0 * Inv_LW;
                Param.BSIM4b1 = ModelParameters.B1
                                + ModelParameters.Lb1 * Inv_L
                                + ModelParameters.Wb1 * Inv_W
                                + ModelParameters.Pb1 * Inv_LW;
                Param.BSIM4alpha0 = ModelParameters.Alpha0
                                    + ModelParameters.Lalpha0 * Inv_L
                                    + ModelParameters.Walpha0 * Inv_W
                                    + ModelParameters.Palpha0 * Inv_LW;
                Param.BSIM4alpha1 = ModelParameters.Alpha1
                                    + ModelParameters.Lalpha1 * Inv_L
                                    + ModelParameters.Walpha1 * Inv_W
                                    + ModelParameters.Palpha1 * Inv_LW;
                Param.BSIM4beta0 = ModelParameters.Beta0
                                   + ModelParameters.Lbeta0 * Inv_L
                                   + ModelParameters.Wbeta0 * Inv_W
                                   + ModelParameters.Pbeta0 * Inv_LW;
                Param.BSIM4agidl = ModelParameters.Agidl
                                   + ModelParameters.Lagidl * Inv_L
                                   + ModelParameters.Wagidl * Inv_W
                                   + ModelParameters.Pagidl * Inv_LW;
                Param.BSIM4bgidl = ModelParameters.Bgidl
                                   + ModelParameters.Lbgidl * Inv_L
                                   + ModelParameters.Wbgidl * Inv_W
                                   + ModelParameters.Pbgidl * Inv_LW;
                Param.BSIM4cgidl = ModelParameters.Cgidl
                                   + ModelParameters.Lcgidl * Inv_L
                                   + ModelParameters.Wcgidl * Inv_W
                                   + ModelParameters.Pcgidl * Inv_LW;
                Param.BSIM4egidl = ModelParameters.Egidl
                                   + ModelParameters.Legidl * Inv_L
                                   + ModelParameters.Wegidl * Inv_W
                                   + ModelParameters.Pegidl * Inv_LW;
                Param.BSIM4rgidl = ModelParameters.Rgidl                /* v4.7 New GIDL/GISL */
                                   + ModelParameters.Lrgidl * Inv_L
                                   + ModelParameters.Wrgidl * Inv_W
                                   + ModelParameters.Prgidl * Inv_LW;
                Param.BSIM4kgidl = ModelParameters.Kgidl                /* v4.7 New GIDL/GISL */
                                   + ModelParameters.Lkgidl * Inv_L
                                   + ModelParameters.Wkgidl * Inv_W
                                   + ModelParameters.Pkgidl * Inv_LW;
                Param.BSIM4fgidl = ModelParameters.Fgidl                /* v4.7 New GIDL/GISL */
                                   + ModelParameters.Lfgidl * Inv_L
                                   + ModelParameters.Wfgidl * Inv_W
                                   + ModelParameters.Pfgidl * Inv_LW;
                Param.BSIM4agisl = ModelParameters.Agisl
                                   + ModelParameters.Lagisl * Inv_L
                                   + ModelParameters.Wagisl * Inv_W
                                   + ModelParameters.Pagisl * Inv_LW;
                Param.BSIM4bgisl = ModelParameters.Bgisl
                                   + ModelParameters.Lbgisl * Inv_L
                                   + ModelParameters.Wbgisl * Inv_W
                                   + ModelParameters.Pbgisl * Inv_LW;
                Param.BSIM4cgisl = ModelParameters.Cgisl
                                   + ModelParameters.Lcgisl * Inv_L
                                   + ModelParameters.Wcgisl * Inv_W
                                   + ModelParameters.Pcgisl * Inv_LW;
                Param.BSIM4egisl = ModelParameters.Egisl
                                   + ModelParameters.Legisl * Inv_L
                                   + ModelParameters.Wegisl * Inv_W
                                   + ModelParameters.Pegisl * Inv_LW;
                Param.BSIM4rgisl = ModelParameters.Rgisl                /* v4.7 New GIDL/GISL */
                                  + ModelParameters.Lrgisl * Inv_L
                                  + ModelParameters.Wrgisl * Inv_W
                                  + ModelParameters.Prgisl * Inv_LW;
                Param.BSIM4kgisl = ModelParameters.Kgisl                /* v4.7 New GIDL/GISL */
                                   + ModelParameters.Lkgisl * Inv_L
                                   + ModelParameters.Wkgisl * Inv_W
                                   + ModelParameters.Pkgisl * Inv_LW;
                Param.BSIM4fgisl = ModelParameters.Fgisl                /* v4.7 New GIDL/GISL */
                                   + ModelParameters.Lfgisl * Inv_L
                                   + ModelParameters.Wfgisl * Inv_W
                                   + ModelParameters.Pfgisl * Inv_LW;
                Param.BSIM4aigc = ModelParameters.Aigc
                                   + ModelParameters.Laigc * Inv_L
                                   + ModelParameters.Waigc * Inv_W
                                   + ModelParameters.Paigc * Inv_LW;
                Param.BSIM4bigc = ModelParameters.Bigc
                                   + ModelParameters.Lbigc * Inv_L
                                   + ModelParameters.Wbigc * Inv_W
                                   + ModelParameters.Pbigc * Inv_LW;
                Param.BSIM4cigc = ModelParameters.Cigc
                                   + ModelParameters.Lcigc * Inv_L
                                   + ModelParameters.Wcigc * Inv_W
                                   + ModelParameters.Pcigc * Inv_LW;
                Param.BSIM4aigsd = ModelParameters.Aigsd
                                   + ModelParameters.Laigsd * Inv_L
                                   + ModelParameters.Waigsd * Inv_W
                                   + ModelParameters.Paigsd * Inv_LW;
                Param.BSIM4bigsd = ModelParameters.Bigsd
                                   + ModelParameters.Lbigsd * Inv_L
                                   + ModelParameters.Wbigsd * Inv_W
                                   + ModelParameters.Pbigsd * Inv_LW;
                Param.BSIM4cigsd = ModelParameters.Cigsd
                                   + ModelParameters.Lcigsd * Inv_L
                                   + ModelParameters.Wcigsd * Inv_W
                                   + ModelParameters.Pcigsd * Inv_LW;
                Param.BSIM4aigs = ModelParameters.Aigs
                                   + ModelParameters.Laigs * Inv_L
                                   + ModelParameters.Waigs * Inv_W
                                   + ModelParameters.Paigs * Inv_LW;
                Param.BSIM4bigs = ModelParameters.Bigs
                                   + ModelParameters.Lbigs * Inv_L
                                   + ModelParameters.Wbigs * Inv_W
                                   + ModelParameters.Pbigs * Inv_LW;
                Param.BSIM4cigs = ModelParameters.Cigs
                                   + ModelParameters.Lcigs * Inv_L
                                   + ModelParameters.Wcigs * Inv_W
                                   + ModelParameters.Pcigs * Inv_LW;
                Param.BSIM4aigd = ModelParameters.Aigd
                                   + ModelParameters.Laigd * Inv_L
                                   + ModelParameters.Waigd * Inv_W
                                   + ModelParameters.Paigd * Inv_LW;
                Param.BSIM4bigd = ModelParameters.Bigd
                                   + ModelParameters.Lbigd * Inv_L
                                   + ModelParameters.Wbigd * Inv_W
                                   + ModelParameters.Pbigd * Inv_LW;
                Param.BSIM4cigd = ModelParameters.Cigd
                                   + ModelParameters.Lcigd * Inv_L
                                   + ModelParameters.Wcigd * Inv_W
                                   + ModelParameters.Pcigd * Inv_LW;
                Param.BSIM4aigbacc = ModelParameters.Aigbacc
                                     + ModelParameters.Laigbacc * Inv_L
                                     + ModelParameters.Waigbacc * Inv_W
                                     + ModelParameters.Paigbacc * Inv_LW;
                Param.BSIM4bigbacc = ModelParameters.Bigbacc
                                     + ModelParameters.Lbigbacc * Inv_L
                                     + ModelParameters.Wbigbacc * Inv_W
                                     + ModelParameters.Pbigbacc * Inv_LW;
                Param.BSIM4cigbacc = ModelParameters.Cigbacc
                                     + ModelParameters.Lcigbacc * Inv_L
                                     + ModelParameters.Wcigbacc * Inv_W
                                     + ModelParameters.Pcigbacc * Inv_LW;
                Param.BSIM4aigbinv = ModelParameters.Aigbinv
                                     + ModelParameters.Laigbinv * Inv_L
                                     + ModelParameters.Waigbinv * Inv_W
                                     + ModelParameters.Paigbinv * Inv_LW;
                Param.BSIM4bigbinv = ModelParameters.Bigbinv
                                     + ModelParameters.Lbigbinv * Inv_L
                                     + ModelParameters.Wbigbinv * Inv_W
                                     + ModelParameters.Pbigbinv * Inv_LW;
                Param.BSIM4cigbinv = ModelParameters.Cigbinv
                                     + ModelParameters.Lcigbinv * Inv_L
                                     + ModelParameters.Wcigbinv * Inv_W
                                     + ModelParameters.Pcigbinv * Inv_LW;
                Param.BSIM4nigc = ModelParameters.Nigc
                                     + ModelParameters.Lnigc * Inv_L
                                     + ModelParameters.Wnigc * Inv_W
                                     + ModelParameters.Pnigc * Inv_LW;
                Param.BSIM4nigbacc = ModelParameters.Nigbacc
                                     + ModelParameters.Lnigbacc * Inv_L
                                     + ModelParameters.Wnigbacc * Inv_W
                                     + ModelParameters.Pnigbacc * Inv_LW;
                Param.BSIM4nigbinv = ModelParameters.Nigbinv
                                     + ModelParameters.Lnigbinv * Inv_L
                                     + ModelParameters.Wnigbinv * Inv_W
                                     + ModelParameters.Pnigbinv * Inv_LW;
                Param.BSIM4ntox = ModelParameters.Ntox
                                  + ModelParameters.Lntox * Inv_L
                                  + ModelParameters.Wntox * Inv_W
                                  + ModelParameters.Pntox * Inv_LW;
                Param.BSIM4eigbinv = ModelParameters.Eigbinv
                                     + ModelParameters.Leigbinv * Inv_L
                                     + ModelParameters.Weigbinv * Inv_W
                                     + ModelParameters.Peigbinv * Inv_LW;
                Param.BSIM4pigcd = ModelParameters.Pigcd
                                   + ModelParameters.Lpigcd * Inv_L
                                   + ModelParameters.Wpigcd * Inv_W
                                   + ModelParameters.Ppigcd * Inv_LW;
                Param.BSIM4poxedge = ModelParameters.Poxedge
                                     + ModelParameters.Lpoxedge * Inv_L
                                     + ModelParameters.Wpoxedge * Inv_W
                                     + ModelParameters.Ppoxedge * Inv_LW;
                Param.BSIM4xrcrg1 = ModelParameters.Xrcrg1
                                    + ModelParameters.Lxrcrg1 * Inv_L
                                    + ModelParameters.Wxrcrg1 * Inv_W
                                    + ModelParameters.Pxrcrg1 * Inv_LW;
                Param.BSIM4xrcrg2 = ModelParameters.Xrcrg2
                                    + ModelParameters.Lxrcrg2 * Inv_L
                                    + ModelParameters.Wxrcrg2 * Inv_W
                                    + ModelParameters.Pxrcrg2 * Inv_LW;
                Param.BSIM4lambda = ModelParameters.Lambda
                                    + ModelParameters.Llambda * Inv_L
                                    + ModelParameters.Wlambda * Inv_W
                                    + ModelParameters.Plambda * Inv_LW;
                Param.BSIM4vtl = ModelParameters.Vtl
                                    + ModelParameters.Lvtl * Inv_L
                                    + ModelParameters.Wvtl * Inv_W
                                    + ModelParameters.Pvtl * Inv_LW;
                Param.BSIM4xn = ModelParameters.Xn
                                    + ModelParameters.Lxn * Inv_L
                                    + ModelParameters.Wxn * Inv_W
                                    + ModelParameters.Pxn * Inv_LW;
                Param.BSIM4vfbsdoff = ModelParameters.Vfbsdoff
                                    + ModelParameters.Lvfbsdoff * Inv_L
                                    + ModelParameters.Wvfbsdoff * Inv_W
                                    + ModelParameters.Pvfbsdoff * Inv_LW;
                Param.BSIM4tvfbsdoff = ModelParameters.Tvfbsdoff
                                    + ModelParameters.Ltvfbsdoff * Inv_L
                                    + ModelParameters.Wtvfbsdoff * Inv_W
                                    + ModelParameters.Ptvfbsdoff * Inv_LW;

                Param.BSIM4cgsl = ModelParameters.Cgsl
                                  + ModelParameters.Lcgsl * Inv_L
                                  + ModelParameters.Wcgsl * Inv_W
                                  + ModelParameters.Pcgsl * Inv_LW;
                Param.BSIM4cgdl = ModelParameters.Cgdl
                                  + ModelParameters.Lcgdl * Inv_L
                                  + ModelParameters.Wcgdl * Inv_W
                                  + ModelParameters.Pcgdl * Inv_LW;
                Param.BSIM4ckappas = ModelParameters.Ckappas
                                     + ModelParameters.Lckappas * Inv_L
                                     + ModelParameters.Wckappas * Inv_W
                                      + ModelParameters.Pckappas * Inv_LW;
                Param.BSIM4ckappad = ModelParameters.Ckappad
                                     + ModelParameters.Lckappad * Inv_L
                                     + ModelParameters.Wckappad * Inv_W
                                     + ModelParameters.Pckappad * Inv_LW;
                Param.BSIM4cf = ModelParameters.Cf
                                + ModelParameters.Lcf * Inv_L
                                + ModelParameters.Wcf * Inv_W
                                + ModelParameters.Pcf * Inv_LW;
                Param.BSIM4clc = ModelParameters.Clc
                                 + ModelParameters.Lclc * Inv_L
                                 + ModelParameters.Wclc * Inv_W
                                 + ModelParameters.Pclc * Inv_LW;
                Param.BSIM4cle = ModelParameters.Cle
                                 + ModelParameters.Lcle * Inv_L
                                 + ModelParameters.Wcle * Inv_W
                                 + ModelParameters.Pcle * Inv_LW;
                Param.BSIM4vfbcv = ModelParameters.Vfbcv
                                   + ModelParameters.Lvfbcv * Inv_L
                                   + ModelParameters.Wvfbcv * Inv_W
                                   + ModelParameters.Pvfbcv * Inv_LW;
                Param.BSIM4acde = ModelParameters.Acde
                                  + ModelParameters.Lacde * Inv_L
                                  + ModelParameters.Wacde * Inv_W
                                  + ModelParameters.Pacde * Inv_LW;
                Param.BSIM4moin = ModelParameters.Moin
                                  + ModelParameters.Lmoin * Inv_L
                                  + ModelParameters.Wmoin * Inv_W
                                  + ModelParameters.Pmoin * Inv_LW;
                Param.BSIM4noff = ModelParameters.Noff
                                  + ModelParameters.Lnoff * Inv_L
                                  + ModelParameters.Wnoff * Inv_W
                                  + ModelParameters.Pnoff * Inv_LW;
                Param.BSIM4voffcv = ModelParameters.Voffcv
                                    + ModelParameters.Lvoffcv * Inv_L
                                    + ModelParameters.Wvoffcv * Inv_W
                                    + ModelParameters.Pvoffcv * Inv_LW;
                Param.BSIM4kvth0we = ModelParameters.Kvth0we
                                    + ModelParameters.Lkvth0we * Inv_L
                                    + ModelParameters.Wkvth0we * Inv_W
                                    + ModelParameters.Pkvth0we * Inv_LW;
                Param.BSIM4k2we = ModelParameters.K2we
                                    + ModelParameters.Lk2we * Inv_L
                                    + ModelParameters.Wk2we * Inv_W
                                    + ModelParameters.Pk2we * Inv_LW;
                Param.BSIM4ku0we = ModelParameters.Ku0we
                                    + ModelParameters.Lku0we * Inv_L
                                    + ModelParameters.Wku0we * Inv_W
                                    + ModelParameters.Pku0we * Inv_LW;

                Param.BSIM4abulkCVfactor = 1.0 + Math.Pow((Param.BSIM4clc
                                           / Param.BSIM4leffCV),
                                           Param.BSIM4cle);

                T0 = (TRatio - 1.0);

                PowWeffWr = Math.Pow(Param.BSIM4weffCJ * 1.0e6, Param.BSIM4wr) * Parameters.Nf;

                T1 = T2 = T3 = T4 = 0.0;
                Param.BSIM4ucs = Param.BSIM4ucs * Math.Pow(TRatio, Param.BSIM4ucste);
                if (ModelParameters.TempMod.Value == 0)
                {
                    Param.BSIM4ua = Param.BSIM4ua + Param.BSIM4ua1 * T0;
                    Param.BSIM4ub = Param.BSIM4ub + Param.BSIM4ub1 * T0;
                    Param.BSIM4uc = Param.BSIM4uc + Param.BSIM4uc1 * T0;
                    Param.BSIM4ud = Param.BSIM4ud + Param.BSIM4ud1 * T0;
                    Param.BSIM4vsattemp = Param.BSIM4vsat - Param.BSIM4at * T0;
                    T10 = Param.BSIM4prt * T0;
                    if (ModelParameters.RdsMod.Value != 0)
                    {
                        /* External Rd(V) */
                        T1 = Param.BSIM4rdw + T10;
                        T2 = ModelParameters.Rdwmin + T10;
                        /* External Rs(V) */
                        T3 = Param.BSIM4rsw + T10;
                        T4 = ModelParameters.Rswmin + T10;
                    }
                    /* Internal Rds(V) in IV */
                    Param.BSIM4rds0 = (Param.BSIM4rdsw + T10)
                                          * Parameters.Nf / PowWeffWr;
                    Param.BSIM4rdswmin = (ModelParameters.Rdswmin + T10)
                                             * Parameters.Nf / PowWeffWr;
                }
                else
                {
                    if (ModelParameters.TempMod.Value == 3)
                    {
                        Param.BSIM4ua = Param.BSIM4ua * Math.Pow(TRatio, Param.BSIM4ua1);
                        Param.BSIM4ub = Param.BSIM4ub * Math.Pow(TRatio, Param.BSIM4ub1);
                        Param.BSIM4uc = Param.BSIM4uc * Math.Pow(TRatio, Param.BSIM4uc1);
                        Param.BSIM4ud = Param.BSIM4ud * Math.Pow(TRatio, Param.BSIM4ud1);
                    }
                    else
                    {
                        /* tempMod = 1, 2 */
                        Param.BSIM4ua = Param.BSIM4ua * (1.0 + Param.BSIM4ua1 * delTemp);
                        Param.BSIM4ub = Param.BSIM4ub * (1.0 + Param.BSIM4ub1 * delTemp);
                        Param.BSIM4uc = Param.BSIM4uc * (1.0 + Param.BSIM4uc1 * delTemp);
                        Param.BSIM4ud = Param.BSIM4ud * (1.0 + Param.BSIM4ud1 * delTemp);
                    }
                    Param.BSIM4vsattemp = Param.BSIM4vsat * (1.0 - Param.BSIM4at * delTemp);
                    T10 = 1.0 + Param.BSIM4prt * delTemp;
                    if (ModelParameters.RdsMod.Value != 0)
                    {
                        /* External Rd(V) */
                        T1 = Param.BSIM4rdw * T10;
                        T2 = ModelParameters.Rdwmin * T10;

                        /* External Rs(V) */
                        T3 = Param.BSIM4rsw * T10;
                        T4 = ModelParameters.Rswmin * T10;
                    }
                    /* Internal Rds(V) in IV */
                    Param.BSIM4rds0 = Param.BSIM4rdsw * T10 * Parameters.Nf / PowWeffWr;
                    Param.BSIM4rdswmin = ModelParameters.Rdswmin * T10 * Parameters.Nf / PowWeffWr;
                }

                if (T1 < 0.0)
                {
                    T1 = 0.0;
                    SpiceSharpWarning.Warning(this, "Warning: Rdw at current temperature is negative; set to 0.");
                }
                if (T2 < 0.0)
                {
                    T2 = 0.0;
                    SpiceSharpWarning.Warning(this, "Warning: Rdwmin at current temperature is negative; set to 0.");
                }
                Param.BSIM4rd0 = T1 / PowWeffWr;
                Param.BSIM4rdwmin = T2 / PowWeffWr;
                if (T3 < 0.0)
                {
                    T3 = 0.0;
                    SpiceSharpWarning.Warning(this, "Warning: Rsw at current temperature is negative; set to 0.");
                }
                if (T4 < 0.0)
                {
                    T4 = 0.0;
                    SpiceSharpWarning.Warning(this, "Warning: Rswmin at current temperature is negative; set to 0.");
                }
                Param.BSIM4rs0 = T3 / PowWeffWr;
                Param.BSIM4rswmin = T4 / PowWeffWr;

                if (Param.BSIM4u0 > 1.0)
                    Param.BSIM4u0 = Param.BSIM4u0 / 1.0e4;

                /* mobility channel length dependence */
                T5 = 1.0 - Param.BSIM4up * Math.Exp(-Param.BSIM4leff / Param.BSIM4lp);
                Param.BSIM4u0temp = Param.BSIM4u0 * T5
                                    * Math.Pow(TRatio, Param.BSIM4ute);
                if (Param.BSIM4eu < 0.0)
                {
                    Param.BSIM4eu = 0.0;
                    SpiceSharpWarning.Warning(this, "Warning: eu has been negative; reset to 0.0.");
                }
                if (Param.BSIM4ucs < 0.0)
                {
                    Param.BSIM4ucs = 0.0;
                    SpiceSharpWarning.Warning(this, "Warning: ucs has been negative; reset to 0.0.");
                }

                Param.BSIM4vfbsdoff = Param.BSIM4vfbsdoff * (1.0 + Param.BSIM4tvfbsdoff * delTemp);
                Param.BSIM4voff = Param.BSIM4voff * (1.0 + Param.BSIM4tvoff * delTemp);

                Param.BSIM4nfactor = Param.BSIM4nfactor + Param.BSIM4tnfactor * delTemp / Tnom;  /* v4.7 temp dep of leakage currents */
                Param.BSIM4voffcv = Param.BSIM4voffcv * (1.0 + Param.BSIM4tvoffcv * delTemp);   /*         v4.7 temp dep of leakage currents */
                Param.BSIM4eta0 = Param.BSIM4eta0 + Param.BSIM4teta0 * delTemp / Tnom;   /*         v4.7 temp dep of leakage currents */

                /* Source End Velocity Limit  */
                if ((ModelParameters.Vtl.Given) && (ModelParameters.Vtl > 0.0))
                {
                    if (ModelParameters.Lc < 0.0) Param.BSIM4lc = 0.0;
                    else Param.BSIM4lc = ModelParameters.Lc;
                    T0 = Param.BSIM4leff / (Param.BSIM4xn * Param.BSIM4leff + Param.BSIM4lc);
                    Param.BSIM4tfactor = (1.0 - T0) / (1.0 + T0);
                }

                Param.BSIM4cgdo = (ModelParameters.Cgdo + Param.BSIM4cf)
                                  * Param.BSIM4weffCV;
                Param.BSIM4cgso = (ModelParameters.Cgso + Param.BSIM4cf)
                                  * Param.BSIM4weffCV;
                Param.BSIM4cgbo = ModelParameters.Cgbo * Param.BSIM4leffCV * Parameters.Nf;

                if (!ModelParameters.Ndep.Given && ModelParameters.Gamma1.Given)
                {
                    T0 = Param.BSIM4gamma1 * ModelTemperature.Coxe;
                    Param.BSIM4ndep = 3.01248e22 * T0 * T0;
                }

                Param.BSIM4phi = ModelTemperature.Vtm0 * Math.Log(Param.BSIM4ndep / ni)
                                 + Param.BSIM4phin + 0.4;

                Param.BSIM4sqrtPhi = Math.Sqrt(Param.BSIM4phi);
                Param.BSIM4phis3 = Param.BSIM4sqrtPhi * Param.BSIM4phi;

                Param.BSIM4Xdep0 = Math.Sqrt(2.0 * epssub / (Constants.Charge
                                   * Param.BSIM4ndep * 1.0e6))
                                   * Param.BSIM4sqrtPhi;
                Param.BSIM4sqrtXdep0 = Math.Sqrt(Param.BSIM4Xdep0);

                if (ModelParameters.MtrlMod.Value == 0)
                    Param.BSIM4litl = Math.Sqrt(3.0 * 3.9 / epsrox * Param.BSIM4xj * toxe);
                else
                    Param.BSIM4litl = Math.Sqrt(ModelParameters.Epsrsub / epsrox * Param.BSIM4xj * toxe);

                Param.BSIM4vbi = Vtm0 * Math.Log(Param.BSIM4nsd
                                 * Param.BSIM4ndep / (ni * ni));

                if (ModelParameters.MtrlMod.Value == 0)
                {
                    if (Param.BSIM4ngate > 0.0)
                    {
                        Param.BSIM4vfbsd = Vtm0 * Math.Log(Param.BSIM4ngate
                                         / Param.BSIM4nsd);
                    }
                    else
                        Param.BSIM4vfbsd = 0.0;
                }
                else
                {
                    T0 = Vtm0 * Math.Log(Param.BSIM4nsd / ni);
                    T1 = 0.5 * Eg0;
                    if (T0 > T1)
                        T0 = T1;
                    T2 = ModelParameters.Easub + T1 - ModelParameters.Type * T0;
                    Param.BSIM4vfbsd = ModelParameters.Phig - T2;
                }

                Param.BSIM4cdep0 = Math.Sqrt(Constants.Charge * epssub
                                   * Param.BSIM4ndep * 1.0e6 / 2.0
                                   / Param.BSIM4phi);

                Param.BSIM4ToxRatio = Math.Exp(Param.BSIM4ntox
                                      * Math.Log(ModelParameters.Toxref / toxe))
                                      / toxe / toxe;
                Param.BSIM4ToxRatioEdge = Math.Exp(Param.BSIM4ntox
                                          * Math.Log(ModelParameters.Toxref
                                          / (toxe * Param.BSIM4poxedge)))
                                          / toxe / toxe
                                          / Param.BSIM4poxedge / Param.BSIM4poxedge;
                Param.BSIM4Aechvb = (ModelParameters.Type > 0) ? 4.97232e-7 : 3.42537e-7;
                Param.BSIM4Bechvb = (ModelParameters.Type > 0) ? 7.45669e11 : 1.16645e12;

                if (ModelParameters.Version.Value != "4.8.1" && ModelParameters.Version.Value != "4.81")
                {
                    Param.BSIM4AechvbEdgeS = Param.BSIM4Aechvb * Param.BSIM4weff
                                            * ModelParameters.Dlcig * Param.BSIM4ToxRatioEdge;
                    Param.BSIM4AechvbEdgeD = Param.BSIM4Aechvb * Param.BSIM4weff
                                            * ModelParameters.Dlcigd * Param.BSIM4ToxRatioEdge;
                }
                else
                {
                    if (ModelParameters.Dlcig < 0.0)
                    {
                        SpiceSharpWarning.Warning(this, "Warning: dlcig = {0:g} is negative. Set to zero.".FormatString(ModelParameters.Dlcig));
                        ModelParameters.Dlcig = 0.0;
                    }
                    Param.BSIM4AechvbEdgeS = Param.BSIM4Aechvb * Param.BSIM4weff
                        * ModelParameters.Dlcig * Param.BSIM4ToxRatioEdge;
                    if (ModelParameters.Dlcigd < 0.0)
                    {
                        SpiceSharpWarning.Warning(this, "Warning: dlcigd = {0:g} is negative. Set to zero.".FormatString(ModelParameters.Dlcigd));
                        ModelParameters.Dlcigd = 0.0;
                    }
                    Param.BSIM4AechvbEdgeD = Param.BSIM4Aechvb * Param.BSIM4weff
                        * ModelParameters.Dlcigd * Param.BSIM4ToxRatioEdge;
                }

                Param.BSIM4BechvbEdge = -Param.BSIM4Bechvb
                                        * toxe * Param.BSIM4poxedge;
                Param.BSIM4Aechvb *= Param.BSIM4weff * Param.BSIM4leff
                                     * Param.BSIM4ToxRatio;
                Param.BSIM4Bechvb *= -toxe;


                Param.BSIM4mstar = 0.5 + Math.Atan(Param.BSIM4minv) / Math.PI;
                Param.BSIM4mstarcv = 0.5 + Math.Atan(Param.BSIM4minvcv) / Math.PI;
                Param.BSIM4voffcbn = Param.BSIM4voff + ModelParameters.Voffl / Param.BSIM4leff;
                Param.BSIM4voffcbncv = Param.BSIM4voffcv + ModelParameters.Voffcvl / Param.BSIM4leff;

                Param.BSIM4ldeb = Math.Sqrt(epssub * Vtm0 / (Constants.Charge
                                  * Param.BSIM4ndep * 1.0e6)) / 3.0;
                Param.BSIM4acde *= Math.Pow((Param.BSIM4ndep / 2.0e16), -0.25);


                if (ModelParameters.K1.Given || ModelParameters.K2.Given)
                {
                    if (!ModelParameters.K1.Given)
                    {
                        SpiceSharpWarning.Warning(this, "Warning: k1 should be specified with k2.");
                        Param.BSIM4k1 = 0.53;
                    }
                    if (!ModelParameters.K2.Given)
                    {
                        SpiceSharpWarning.Warning(this, "Warning: k2 should be specified with k1.");
                        Param.BSIM4k2 = -0.0186;
                    }
                    /* don't print in sensitivity */
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
                        Param.BSIM4vbx = Param.BSIM4phi - 7.7348e-4
                                         * Param.BSIM4ndep
                                         * Param.BSIM4xt * Param.BSIM4xt;
                    if (Param.BSIM4vbx > 0.0)
                        Param.BSIM4vbx = -Param.BSIM4vbx;
                    if (Param.BSIM4vbm > 0.0)
                        Param.BSIM4vbm = -Param.BSIM4vbm;

                    if (!ModelParameters.Gamma1.Given)
                        Param.BSIM4gamma1 = 5.753e-12
                                            * Math.Sqrt(Param.BSIM4ndep)
                                            / ModelTemperature.Coxe;
                    if (!ModelParameters.Gamma2.Given)
                        Param.BSIM4gamma2 = 5.753e-12
                                            * Math.Sqrt(Param.BSIM4nsub)
                                            / ModelTemperature.Coxe;

                    T0 = Param.BSIM4gamma1 - Param.BSIM4gamma2;
                    T1 = Math.Sqrt(Param.BSIM4phi - Param.BSIM4vbx)
                       - Param.BSIM4sqrtPhi;
                    T2 = Math.Sqrt(Param.BSIM4phi * (Param.BSIM4phi
                       - Param.BSIM4vbm)) - Param.BSIM4phi;
                    Param.BSIM4k2 = T0 * T1 / (2.0 * T2 + Param.BSIM4vbm);
                    Param.BSIM4k1 = Param.BSIM4gamma2 - 2.0
                                    * Param.BSIM4k2 * Math.Sqrt(Param.BSIM4phi
                                    - Param.BSIM4vbm);
                }

                if (!ModelParameters.Vfb.Given)
                {
                    if (ModelParameters.Vth0.Given)
                    {
                        Param.BSIM4vfb = ModelParameters.Type * Param.BSIM4vth0
                                         - Param.BSIM4phi - Param.BSIM4k1
                                         * Param.BSIM4sqrtPhi;
                    }
                    else
                    {
                        if ((ModelParameters.MtrlMod.Value != 0) && (ModelParameters.Phig.Given) &&
                            (ModelParameters.Nsub.Given))
                        {
                            T0 = Vtm0 * Math.Log(Param.BSIM4nsub / ni);
                            T1 = 0.5 * Eg0;
                            if (T0 > T1)
                                T0 = T1;
                            T2 = ModelParameters.Easub + T1 + ModelParameters.Type * T0;
                            Param.BSIM4vfb = ModelParameters.Phig - T2;
                        }
                        else
                        {
                            Param.BSIM4vfb = -1.0;
                        }
                    }
                }
                if (!ModelParameters.Vth0.Given)
                {
                    Param.BSIM4vth0 = ModelParameters.Type * (Param.BSIM4vfb
                                      + Param.BSIM4phi + Param.BSIM4k1
                                      * Param.BSIM4sqrtPhi);
                }

                Param.BSIM4k1ox = Param.BSIM4k1 * toxe
                                  / ModelParameters.Toxm;

                tmp = Math.Sqrt(epssub / (epsrox * EPS0) * toxe * Param.BSIM4Xdep0);
                T0 = Param.BSIM4dsub * Param.BSIM4leff / tmp;
                if (T0 < EXP_THRESHOLD)
                {
                    T1 = Math.Exp(T0);
                    T2 = T1 - 1.0;
                    T3 = T2 * T2;
                    T4 = T3 + 2.0 * T1 * MIN_EXP;
                    Param.BSIM4theta0vb0 = T1 / T4;
                }
                else
                    Param.BSIM4theta0vb0 = 1.0 / (MAX_EXP - 2.0);

                T0 = Param.BSIM4drout * Param.BSIM4leff / tmp;
                if (T0 < EXP_THRESHOLD)
                {
                    T1 = Math.Exp(T0);
                    T2 = T1 - 1.0;
                    T3 = T2 * T2;
                    T4 = T3 + 2.0 * T1 * MIN_EXP;
                    T5 = T1 / T4;
                }
                else
                    T5 = 1.0 / (MAX_EXP - 2.0); /* 3.0 * MIN_EXP omitted */
                Param.BSIM4thetaRout = Param.BSIM4pdibl1 * T5
                                       + Param.BSIM4pdibl2;

                tmp = Math.Sqrt(Param.BSIM4Xdep0);
                tmp1 = Param.BSIM4vbi - Param.BSIM4phi;
                tmp2 = ModelTemperature.Factor1 * tmp;

                T0 = Param.BSIM4dvt1w * Param.BSIM4weff
                   * Param.BSIM4leff / tmp2;
                if (T0 < EXP_THRESHOLD)
                {
                    T1 = Math.Exp(T0);
                    T2 = T1 - 1.0;
                    T3 = T2 * T2;
                    T4 = T3 + 2.0 * T1 * MIN_EXP;
                    T8 = T1 / T4;
                }
                else
                    T8 = 1.0 / (MAX_EXP - 2.0);
                T0 = Param.BSIM4dvt0w * T8;
                T8 = T0 * tmp1;

                T0 = Param.BSIM4dvt1 * Param.BSIM4leff / tmp2;
                if (T0 < EXP_THRESHOLD)
                {
                    T1 = Math.Exp(T0);
                    T2 = T1 - 1.0;
                    T3 = T2 * T2;
                    T4 = T3 + 2.0 * T1 * MIN_EXP;
                    T9 = T1 / T4;
                }
                else
                    T9 = 1.0 / (MAX_EXP - 2.0);
                T9 = Param.BSIM4dvt0 * T9 * tmp1;

                T4 = toxe * Param.BSIM4phi
                   / (Param.BSIM4weff + Param.BSIM4w0);

                T0 = Math.Sqrt(1.0 + Param.BSIM4lpe0 / Param.BSIM4leff);
                if ((ModelParameters.TempMod.Value == 1) || (ModelParameters.TempMod.Value == 0))
                    T3 = (Param.BSIM4kt1 + Param.BSIM4kt1l / Param.BSIM4leff)
                               * (TRatio - 1.0);
                if ((ModelParameters.TempMod.Value == 2) || (ModelParameters.TempMod.Value == 3))
                    T3 = -Param.BSIM4kt1 * (TRatio - 1.0);

                T5 = Param.BSIM4k1ox * (T0 - 1.0) * Param.BSIM4sqrtPhi
                   + T3;
                Param.BSIM4vfbzbfactor = -T8 - T9 + Param.BSIM4k3 * T4 + T5
                                           - Param.BSIM4phi - Param.BSIM4k1 * Param.BSIM4sqrtPhi;

                /* stress effect */

                wlod = ModelParameters.Wlod;
                if (ModelParameters.Wlod < 0.0)
                {
                    SpiceSharpWarning.Warning(this, "Warning: WLOD = {0:g} is less than 0. 0.0 is used".FormatString(ModelParameters.Wlod));
                    wlod = 0.0;
                }
                T0 = Math.Pow(Lnew, ModelParameters.Llodku0);
                W_tmp = Wnew + wlod;
                T1 = Math.Pow(W_tmp, ModelParameters.Wlodku0);
                tmp1 = ModelParameters.Lku0 / T0 + ModelParameters.Wku0 / T1
                       + ModelParameters.Pku0 / (T0 * T1);
                Param.BSIM4ku0 = 1.0 + tmp1;

                T0 = Math.Pow(Lnew, ModelParameters.Llodvth);
                T1 = Math.Pow(W_tmp, ModelParameters.Wlodvth);
                tmp1 = ModelParameters.Lkvth0 / T0 + ModelParameters.Wkvth0 / T1
                     + ModelParameters.Pkvth0 / (T0 * T1);
                Param.BSIM4kvth0 = 1.0 + tmp1;
                Param.BSIM4kvth0 = Math.Sqrt(Param.BSIM4kvth0 * Param.BSIM4kvth0 + DELTA);

                T0 = (TRatio - 1.0);
                Param.BSIM4ku0temp = Param.BSIM4ku0 * (1.0 + ModelParameters.Tku0 * T0) + DELTA;

                Inv_saref = 1.0 / (ModelParameters.Saref + 0.5 * Ldrn);
                Inv_sbref = 1.0 / (ModelParameters.Sbref + 0.5 * Ldrn);
                Param.BSIM4inv_od_ref = Inv_saref + Inv_sbref;
                Param.BSIM4rho_ref = ModelParameters.Ku0 / Param.BSIM4ku0temp * Param.BSIM4inv_od_ref;

                /*high k*/
                /*Calculate VgsteffVth for mobMod=3*/
                if (ModelParameters.MobMod.Value == 3)
                {        /*Calculate n @ Vbs=Vds=0*/
                    lt1 = ModelTemperature.Factor1 * Param.BSIM4sqrtXdep0;
                    T0 = Param.BSIM4dvt1 * Param.BSIM4leff / lt1;
                    if (T0 < EXP_THRESHOLD)
                    {
                        T1 = Math.Exp(T0);
                        T2 = T1 - 1.0;
                        T3 = T2 * T2;
                        T4 = T3 + 2.0 * T1 * MIN_EXP;
                        Theta0 = T1 / T4;
                    }
                    else
                        Theta0 = 1.0 / (MAX_EXP - 2.0);

                    tmp1 = epssub / Param.BSIM4Xdep0;
                    tmp2 = Param.BSIM4nfactor * tmp1;
                    tmp3 = (tmp2 + Param.BSIM4cdsc * Theta0 + Param.BSIM4cit) / ModelTemperature.Coxe;
                    if (tmp3 >= -0.5)
                        n0 = 1.0 + tmp3;
                    else
                    {
                        T0 = 1.0 / (3.0 + 8.0 * tmp3);
                        n0 = (1.0 + 3.0 * tmp3) * T0;
                    }

                    T0 = n0 * ModelTemperature.Vtm;
                    T1 = Param.BSIM4voffcbn;
                    T2 = T1 / T0;
                    if (T2 < -EXP_THRESHOLD)
                    {
                        T3 = ModelTemperature.Coxe * MIN_EXP / Param.BSIM4cdep0;
                        T4 = Param.BSIM4mstar + T3 * n0;
                    }
                    else if (T2 > EXP_THRESHOLD)
                    {
                        T3 = ModelTemperature.Coxe * MAX_EXP / Param.BSIM4cdep0;
                        T4 = Param.BSIM4mstar + T3 * n0;
                    }
                    else
                    {
                        T3 = Math.Exp(T2) * ModelTemperature.Coxe / Param.BSIM4cdep0;
                        T4 = Param.BSIM4mstar + T3 * n0;
                    }
                    Param.BSIM4VgsteffVth = T0 * Math.Log(2.0) / T4;
                }

                /* New DITS term added in 4.7 */
                T0 = -Param.BSIM4dvtp3 * Math.Log(Param.BSIM4leff);
                T1 = DEXP(T0);
                Param.BSIM4dvtp2factor = Param.BSIM4dvtp5 + Param.BSIM4dvtp2 * T1;

            } /* End of SizeNotFound */

            /*  stress effect */
            if ((Parameters.Sa > 0.0) && (Parameters.Sb > 0.0) &&
                      ((Parameters.Nf == 1.0) || ((Parameters.Nf > 1.0) && (Parameters.Sd > 0.0))))
            {
                Inv_sa = 0;
                Inv_sb = 0;

                kvsat = ModelParameters.Kvsat;
                if (ModelParameters.Kvsat < -1.0)
                {
                    SpiceSharpWarning.Warning(this, "Warning: KVSAT = {0:g} is too small; -1.0 is used.".FormatString(ModelParameters.Kvsat));
                    kvsat = -1.0;
                }
                if (ModelParameters.Kvsat > 1.0)
                {
                    SpiceSharpWarning.Warning(this, "Warning: KVSAT = {0:g} is too big; 1.0 is used.".FormatString(ModelParameters.Kvsat));
                    kvsat = 1.0;
                }

                for (int i = 0; i < Parameters.Nf; i++)
                {
                    T0 = 1.0 / Parameters.Nf / (Parameters.Sa + 0.5 * Ldrn + i * (Parameters.Sd + Ldrn));
                    T1 = 1.0 / Parameters.Nf / (Parameters.Sb + 0.5 * Ldrn + i * (Parameters.Sd + Ldrn));
                    Inv_sa += T0;
                    Inv_sb += T1;
                }
                Inv_ODeff = Inv_sa + Inv_sb;
                rho = ModelParameters.Ku0 / Param.BSIM4ku0temp * Inv_ODeff;
                T0 = (1.0 + rho) / (1.0 + Param.BSIM4rho_ref);
                this._u0temp = Param.BSIM4u0temp * T0;

                T1 = (1.0 + kvsat * rho) / (1.0 + kvsat * Param.BSIM4rho_ref);
                this._vsattemp = Param.BSIM4vsattemp * T1;

                OD_offset = Inv_ODeff - Param.BSIM4inv_od_ref;
                dvth0_lod = ModelParameters.Kvth0 / Param.BSIM4kvth0 * OD_offset;
                dk2_lod = ModelParameters.Stk2 / Math.Pow(Param.BSIM4kvth0, ModelParameters.Lodk2) *
                                 OD_offset;
                deta0_lod = ModelParameters.Steta0 / Math.Pow(Param.BSIM4kvth0, ModelParameters.Lodeta0) *
                                   OD_offset;
                this._vth0 = Param.BSIM4vth0 + dvth0_lod;

                this._eta0 = Param.BSIM4eta0 + deta0_lod;
                this._k2 = Param.BSIM4k2 + dk2_lod;
            }
            else
            {
                this._u0temp = Param.BSIM4u0temp;
                this._vth0 = Param.BSIM4vth0;
                this._vsattemp = Param.BSIM4vsattemp;
                this._eta0 = Param.BSIM4eta0;
                this._k2 = Param.BSIM4k2;
            }

            /*  Well Proximity Effect  */
            if (ModelParameters.Wpemod.Value != 0)
            {
                if ((!Parameters.Sca.Given) && (!Parameters.Scb.Given) && (!Parameters.Scc.Given))
                {
                    if ((Parameters.Sc.Given) && (Parameters.Sc > 0.0))
                    {
                        T1 = Parameters.Sc + Wdrn;
                        T2 = 1.0 / ModelParameters.Scref;
                        Parameters.Sca = ModelParameters.Scref * ModelParameters.Scref
                                        / (Parameters.Sc * T1);
                        Parameters.Scb = ((0.1 * Parameters.Sc + 0.01 * ModelParameters.Scref)
                                        * Math.Exp(-10.0 * Parameters.Sc * T2)
                                        - (0.1 * T1 + 0.01 * ModelParameters.Scref)
                                        * Math.Exp(-10.0 * T1 * T2)) / Wdrn;
                        Parameters.Scc = ((0.05 * Parameters.Sc + 0.0025 * ModelParameters.Scref)
                                        * Math.Exp(-20.0 * Parameters.Sc * T2)
                                        - (0.05 * T1 + 0.0025 * ModelParameters.Scref)
                                        * Math.Exp(-20.0 * T1 * T2)) / Wdrn;
                    }
                    else
                    {
                        SpiceSharpWarning.Warning(this, "Warning: No WPE as none of SCA, SCB, SCC, SC is given and/or SC not positive.");
                    }
                }

                if (Parameters.Sca < 0.0)
                {
                    SpiceSharpWarning.Warning(this, "Warning: SCA = {0:g} is negative. Set to 0.0.".FormatString(Parameters.Sca));
                    Parameters.Sca = 0.0;
                }
                if (Parameters.Scb < 0.0)
                {
                    SpiceSharpWarning.Warning(this, "Warning: SCB = {0:g} is negative. Set to 0.0.".FormatString(Parameters.Scb));
                    Parameters.Scb = 0.0;
                }
                if (Parameters.Scc < 0.0)
                {
                    SpiceSharpWarning.Warning(this, "Warning: SCC = {0:g} is negative. Set to 0.0.".FormatString(Parameters.Scc));
                    Parameters.Scc = 0.0;
                }
                if (Parameters.Sc < 0.0)
                {
                    SpiceSharpWarning.Warning(this, "Warning: SC = {0:g} is negative. Set to 0.0.".FormatString(Parameters.Sc));
                    Parameters.Sc = 0.0;
                }
                /*4.6.2*/
                sceff = Parameters.Sca + ModelParameters.Web * Parameters.Scb
                      + ModelParameters.Wec * Parameters.Scc;
                this._vth0 += Param.BSIM4kvth0we * sceff;
                this._k2 += Param.BSIM4k2we * sceff;
                T3 = 1.0 + Param.BSIM4ku0we * sceff;
                if (T3 <= 0.0)
                {
                    T3 = 0.0;
                    SpiceSharpWarning.Warning(this, "Warning: ku0we = {0:g} is negatively too high. Negative mobility! ".FormatString(Param.BSIM4ku0we));
                }
                this._u0temp *= T3;
            }

            /* adding delvto  */
            this._vth0 += Parameters.Delvto;
            this._vfb = Param.BSIM4vfb + ModelParameters.Type * Parameters.Delvto;

            /* low field mobility multiplier */
            this._u0temp = Param.BSIM4u0temp * Parameters.Mulu0;

            /* Instance variables calculation  */
            T3 = ModelParameters.Type * this._vth0
               - this._vfb - Param.BSIM4phi;
            T4 = T3 + T3;
            T5 = 2.5 * T3;
            this._vtfbphi1 = (ModelParameters.Type > 0) ? T4 : T5;
            if (this._vtfbphi1 < 0.0)
                this._vtfbphi1 = 0.0;

            this._vtfbphi2 = 4.0 * T3;
            if (this._vtfbphi2 < 0.0)
                this._vtfbphi2 = 0.0;

            if (this._k2 < 0.0)
            {
                T0 = 0.5 * Param.BSIM4k1 / this._k2;
                this._vbsc = 0.9 * (Param.BSIM4phi - T0 * T0);
                if (this._vbsc > -3.0)
                    this._vbsc = -3.0;
                else if (this._vbsc < -30.0)
                    this._vbsc = -30.0;
            }
            else
                this._vbsc = -30.0;
            if (this._vbsc > Param.BSIM4vbm)
                this._vbsc = Param.BSIM4vbm;
            this._k2ox = this._k2 * toxe
                              / ModelParameters.Toxm;

            this._vfbzb = Param.BSIM4vfbzbfactor
                                + ModelParameters.Type * this._vth0;

            this._cgso = Param.BSIM4cgso;
            this._cgdo = Param.BSIM4cgdo;

            lnl = Math.Log(Param.BSIM4leff * 1.0e6);
            lnw = Math.Log(Param.BSIM4weff * 1.0e6);
            lnnf = Math.Log(Parameters.Nf);

            bodymode = 5;
            if ((!ModelParameters.Rbps0.Given) ||
                (!ModelParameters.Rbpd0.Given))
                bodymode = 1;
            else
              if ((!ModelParameters.Rbsbx0.Given && !ModelParameters.Rbsby0.Given) ||
                    (!ModelParameters.Rbdbx0.Given && !ModelParameters.Rbdby0.Given))
                bodymode = 3;

            if (Parameters.RbodyMod.Value == 2)
            {
                if (bodymode == 5)
                {
                    /*rbsbx =  Math.Exp( Math.Log(ModelParameters.Rbsbx0) + ModelParameters.Rbsdbxl * lnl +
                          ModelParameters.Rbsdbxw * lnw + ModelParameters.Rbsdbxnf * lnnf );
                    rbsby =  Math.Exp( Math.Log(ModelParameters.Rbsby0) + ModelParameters.Rbsdbyl * lnl +
                          ModelParameters.Rbsdbyw * lnw + ModelParameters.Rbsdbynf * lnnf );
                           */
                    rbsbx = ModelParameters.Rbsbx0 * Math.Exp(ModelParameters.Rbsdbxl * lnl +
                          ModelParameters.Rbsdbxw * lnw + ModelParameters.Rbsdbxnf * lnnf);
                    rbsby = ModelParameters.Rbsby0 * Math.Exp(ModelParameters.Rbsdbyl * lnl +
                          ModelParameters.Rbsdbyw * lnw + ModelParameters.Rbsdbynf * lnnf);
                    Parameters.Rbsb = rbsbx * rbsby / (rbsbx + rbsby);


                    /*rbdbx =  Math.Exp( Math.Log(ModelParameters.Rbdbx0) + ModelParameters.Rbsdbxl * lnl +
                          ModelParameters.Rbsdbxw * lnw + ModelParameters.Rbsdbxnf * lnnf );
                    rbdby =  Math.Exp( Math.Log(ModelParameters.Rbdby0) + ModelParameters.Rbsdbyl * lnl +
                          ModelParameters.Rbsdbyw * lnw + ModelParameters.Rbsdbynf * lnnf );
                            */

                    rbdbx = ModelParameters.Rbdbx0 * Math.Exp(ModelParameters.Rbsdbxl * lnl +
                          ModelParameters.Rbsdbxw * lnw + ModelParameters.Rbsdbxnf * lnnf);
                    rbdby = ModelParameters.Rbdby0 * Math.Exp(ModelParameters.Rbsdbyl * lnl +
                          ModelParameters.Rbsdbyw * lnw + ModelParameters.Rbsdbynf * lnnf);

                    Parameters.Rbdb = rbdbx * rbdby / (rbdbx + rbdby);
                }

                if ((bodymode == 3) || (bodymode == 5))
                {
                    /*Parameters.Rbps = Math.Exp( Math.Log(ModelParameters.Rbps0) + ModelParameters.Rbpsl * lnl +
                               ModelParameters.Rbpsw * lnw + ModelParameters.Rbpsnf * lnnf );
                    Parameters.Rbpd = Math.Exp( Math.Log(ModelParameters.Rbpd0) + ModelParameters.Rbpdl * lnl +
                                             ModelParameters.Rbpdw * lnw + ModelParameters.Rbpdnf * lnnf );
                           */
                    Parameters.Rbps = ModelParameters.Rbps0 * Math.Exp(ModelParameters.Rbpsl * lnl +
                               ModelParameters.Rbpsw * lnw + ModelParameters.Rbpsnf * lnnf);
                    Parameters.Rbpd = ModelParameters.Rbpd0 * Math.Exp(ModelParameters.Rbpdl * lnl +
                                             ModelParameters.Rbpdw * lnw + ModelParameters.Rbpdnf * lnnf);

                }

                /*rbpbx =  Math.Exp( Math.Log(ModelParameters.Rbpbx0) + ModelParameters.Rbpbxl * lnl +
                          ModelParameters.Rbpbxw * lnw + ModelParameters.Rbpbxnf * lnnf );
                rbpby =  Math.Exp( Math.Log(ModelParameters.Rbpby0) + ModelParameters.Rbpbyl * lnl +
                          ModelParameters.Rbpbyw * lnw + ModelParameters.Rbpbynf * lnnf );
                        */
                rbpbx = ModelParameters.Rbpbx0 * Math.Exp(ModelParameters.Rbpbxl * lnl +
                  ModelParameters.Rbpbxw * lnw + ModelParameters.Rbpbxnf * lnnf);
                rbpby = ModelParameters.Rbpby0 * Math.Exp(ModelParameters.Rbpbyl * lnl +
                          ModelParameters.Rbpbyw * lnw + ModelParameters.Rbpbynf * lnnf);

                Parameters.Rbpb = rbpbx * rbpby / (rbpbx + rbpby);
            }

            if ((Parameters.RbodyMod.Value == 1) || ((Parameters.RbodyMod.Value == 2) && (bodymode == 5)))
            {
                if (Parameters.Rbdb < 1.0e-3)
                    this._grbdb = 1.0e3; /* in mho */
                else
                    this._grbdb = ModelParameters.Gbmin + 1.0 / Parameters.Rbdb;
                if (Parameters.Rbpb < 1.0e-3)
                    this._grbpb = 1.0e3;
                else
                    this._grbpb = ModelParameters.Gbmin + 1.0 / Parameters.Rbpb;
                if (Parameters.Rbps < 1.0e-3)
                    this._grbps = 1.0e3;
                else
                    this._grbps = ModelParameters.Gbmin + 1.0 / Parameters.Rbps;
                if (Parameters.Rbsb < 1.0e-3)
                    this._grbsb = 1.0e3;
                else
                    this._grbsb = ModelParameters.Gbmin + 1.0 / Parameters.Rbsb;
                if (Parameters.Rbpd < 1.0e-3)
                    this._grbpd = 1.0e3;
                else
                    this._grbpd = ModelParameters.Gbmin + 1.0 / Parameters.Rbpd;

            }

            if ((Parameters.RbodyMod.Value == 2) && (bodymode == 3))
            {
                this._grbdb = this._grbsb = ModelParameters.Gbmin;
                if (Parameters.Rbpb < 1.0e-3)
                    this._grbpb = 1.0e3;
                else
                    this._grbpb = ModelParameters.Gbmin + 1.0 / Parameters.Rbpb;
                if (Parameters.Rbps < 1.0e-3)
                    this._grbps = 1.0e3;
                else
                    this._grbps = ModelParameters.Gbmin + 1.0 / Parameters.Rbps;
                if (Parameters.Rbpd < 1.0e-3)
                    this._grbpd = 1.0e3;
                else
                    this._grbpd = ModelParameters.Gbmin + 1.0 / Parameters.Rbpd;
            }

            if ((Parameters.RbodyMod.Value == 2) && (bodymode == 1))
            {
                this._grbdb = this._grbsb = ModelParameters.Gbmin;
                this._grbps = this._grbpd = 1.0e3;
                if (Parameters.Rbpb < 1.0e-3)
                    this._grbpb = 1.0e3;
                else
                    this._grbpb = ModelParameters.Gbmin + 1.0 / Parameters.Rbpb;
            }


            /*
             * Process geomertry dependent parasitics
             */

            this._grgeltd = ModelParameters.Rshg * (Parameters.Xgw
                    + Param.BSIM4weffCJ / 3.0 / Parameters.Ngcon) /
                    (Parameters.Ngcon * Parameters.Nf *
                    (Lnew - ModelParameters.Xgl));
            if (this._grgeltd > 0.0)
                this._grgeltd = 1.0 / this._grgeltd;
            else
            {
                this._grgeltd = 1.0e3; /* mho */
                if (Parameters.RgateMod.Value != 0)
                    SpiceSharpWarning.Warning(this, "Warning: The gate conductance reset to 1.0e3 mho.");
            }

            DMCGeff = ModelParameters.Dmcg - ModelParameters.Dmcgt;
            DMCIeff = ModelParameters.Dmci;
            DMDGeff = ModelParameters.Dmdg - ModelParameters.Dmcgt;

            /*              if (Parameters.SourcePerimeter.Given)
                          {   if (ModelParameters.PerMod.Value == 0)
                                  this._pseff = Parameters.SourcePerimeter;
                              else
                                  this._pseff = Parameters.SourcePerimeter
                                                   - Param.BSIM4weffCJ * Parameters.Nf;
                          }
                          else
                              BSIM4PAeffGeo(Parameters.Nf, Parameters.GeoMod, Parameters.Min,
                                                Param.BSIM4weffCJ, DMCGeff, DMCIeff, DMDGeff,
                                                out (this._pseff), &dumPd, &dumAs, &dumAd);
                          if (this._pseff < 0.0) /4.6.2/
                                  this._pseff = 0.0; */

            /* New Diode Model v4.7*/
            if (Parameters.SourcePerimeter.Given)
            {   /* given */
                if (Parameters.SourcePerimeter == 0.0)
                    this._pseff = 0.0;
                else if (Parameters.SourcePerimeter < 0.0)
                {
                    SpiceSharpWarning.Warning(this, "Warning: Source Perimeter is specified as negative, it is set to zero.");
                    this._pseff = 0.0;
                }
                else
                {
                    if (ModelParameters.PerMod.Value == 0)
                        this._pseff = Parameters.SourcePerimeter;
                    else
                        this._pseff = Parameters.SourcePerimeter
                                - Param.BSIM4weffCJ * Parameters.Nf;
                }
            }
            else /* not given */
                BSIM4PAeffGeo(Parameters.Nf, Parameters.GeoMod, Parameters.Min,
                                  Param.BSIM4weffCJ, DMCGeff, DMCIeff, DMDGeff,
                                  out this._pseff, out dumPd, out dumAs, out dumAd);

            if (this._pseff < 0.0)
            { /* v4.7 final check */
                this._pseff = 0.0;
                SpiceSharpWarning.Warning(this, "Warning: Pseff is negative, it is set to zero.");
            }
            /*  if (Parameters.DrainPerimeter.Given)
            {   if (ModelParameters.PerMod.Value == 0)
                    this._pdeff = Parameters.DrainPerimeter;
                else
                    this._pdeff = Parameters.DrainPerimeter
                                     - Param.BSIM4weffCJ * Parameters.Nf;
            }
            else
                BSIM4PAeffGeo(Parameters.Nf, Parameters.GeoMod, Parameters.Min,
                                  Param.BSIM4weffCJ, DMCGeff, DMCIeff, DMDGeff,
                                  &dumPs, out (this._pdeff), &dumAs, &dumAd);
             if (this._pdeff < 0.0) /4.6.2/
                    this._pdeff = 0.0; */

            if (Parameters.DrainPerimeter.Given)
            {   /* given */
                if (Parameters.DrainPerimeter == 0.0)
                    this._pdeff = 0.0;
                else if (Parameters.DrainPerimeter < 0.0)
                {
                    SpiceSharpWarning.Warning(this, "Warning: Drain Perimeter is specified as negative, it is set to zero.");
                    this._pdeff = 0.0;
                }
                else
                {
                    if (ModelParameters.PerMod.Value == 0)
                        this._pdeff = Parameters.DrainPerimeter;
                    else
                        this._pdeff = Parameters.DrainPerimeter
                                - Param.BSIM4weffCJ * Parameters.Nf;
                }
            }
            else /* not given */
                BSIM4PAeffGeo(Parameters.Nf, Parameters.GeoMod, Parameters.Min,
                  Param.BSIM4weffCJ, DMCGeff, DMCIeff, DMDGeff,
                                  out dumPs, out (this._pdeff), out dumAs, out dumAd);

            if (this._pdeff < 0.0)
            {
                this._pdeff = 0.0; /*New Diode v4.7*/
                SpiceSharpWarning.Warning(this, "Warning: Pdeff is negative, it is set to zero.");
            }
            if (Parameters.SourceArea.Given)
                this._aseff = Parameters.SourceArea;
            else
                BSIM4PAeffGeo(Parameters.Nf, Parameters.GeoMod, Parameters.Min,
                                  Param.BSIM4weffCJ, DMCGeff, DMCIeff, DMDGeff,
                                  out dumPs, out dumPd, out (this._aseff), out dumAd);
            if (this._aseff < 0.0)
            {
                this._aseff = 0.0; /* v4.7 */
                SpiceSharpWarning.Warning(this, "Warning: Aseff is negative, it is set to zero.");
            }
            if (Parameters.DrainArea.Given)
                this._adeff = Parameters.DrainArea;
            else
                BSIM4PAeffGeo(Parameters.Nf, Parameters.GeoMod, Parameters.Min,
                                  Param.BSIM4weffCJ, DMCGeff, DMCIeff, DMDGeff,
                                  out dumPs, out dumPd, out dumAs, out (this._adeff));
            if (this._adeff < 0.0)
            {
                this._adeff = 0.0; /* v4.7 */
                SpiceSharpWarning.Warning(this, "Warning: Adeff is negative, it is set to zero.");
            }

            /* Processing S/D resistance and conductance below */
            // 20220724 - merged setup and temperature-dependent code - Sven Boulanger
            bool createNode = false;
            if ((ModelParameters.RdsMod != 0) || (ModelParameters.TnoiMod == 1 && (this is INoiseBehavior)))
            {
                createNode = true;
            }
            else if (ModelParameters.SheetResistance > 0)
            {
                if (Parameters.DrainSquares.Given && Parameters.DrainSquares > 0)
                {
                    createNode = true;
                }
                else if (!Parameters.DrainSquares.Given && (Parameters.RgeoMod != 0))
                {
                    BSIM4RdseffGeo(Parameters.Nf * Parameters.M, Parameters.GeoMod,
                            Parameters.RgeoMod, Parameters.Min,
                            Parameters.W, ModelParameters.SheetResistance,
                            DMCGeff, DMCIeff, DMDGeff, 0, out var Rtot);
                    if (Rtot > 0)
                        createNode = true;
                }
            }
            if (createNode)
            {
                this._drainConductance = 0.0;
                if (Parameters.DrainSquares.Given)
                {
                    this._drainConductance = ModelParameters.SheetResistance
                                              * Parameters.DrainSquares;
                }
                else if (Parameters.RgeoMod > 0)
                {
                    BSIM4RdseffGeo(Parameters.Nf, Parameters.GeoMod,
                      Parameters.RgeoMod, Parameters.Min,
                      Param.BSIM4weffCJ, ModelParameters.SheetResistance,
                  DMCGeff, DMCIeff, DMDGeff, 0, out this._drainConductance);
                }
                else
                {
                    this._drainConductance = 0.0;
                }

                if (this._drainConductance > 0.0)
                    this._drainConductance = 1.0
                                          / this._drainConductance;
                else
                {
                    this._drainConductance = 1.0e3; /* mho */
                    SpiceSharpWarning.Warning(this, "Warning: Drain conductance reset to 1.0e3 mho.");
                }
            }
            else
            {
                this._drainConductance = 0.0;
            }

            createNode = false;
            if ((ModelParameters.RdsMod != 0) || (ModelParameters.TnoiMod == 1 && (this is INoiseBehavior)))
            {
                createNode = true;
            } else if (ModelParameters.SheetResistance > 0)
            {
                if (Parameters.SourceSquares.Given && Parameters.SourceSquares > 0)
                {
                    createNode = true;
                }
                else if (!Parameters.SourceSquares.Given && (Parameters.RgeoMod != 0))
                {
                    BSIM4RdseffGeo(Parameters.Nf * Parameters.M, Parameters.GeoMod,
                            Parameters.RgeoMod, Parameters.Min,
                            Parameters.W, ModelParameters.SheetResistance,
                            DMCGeff, DMCIeff, DMDGeff, 1, out var Rtot);
                    if (Rtot > 0)
                        createNode = true;
                }
            }
            if (createNode)
            {
                this._sourceConductance = 0.0;
                if (Parameters.SourceSquares.Given)
                {
                    this._sourceConductance = ModelParameters.SheetResistance
                                               * Parameters.SourceSquares;
                }
                else if (Parameters.RgeoMod > 0)
                {
                    BSIM4RdseffGeo(Parameters.Nf, Parameters.GeoMod,
                      Parameters.RgeoMod, Parameters.Min,
                      Param.BSIM4weffCJ, ModelParameters.SheetResistance,
                  DMCGeff, DMCIeff, DMDGeff, 1, out this._sourceConductance);
                }
                else
                {
                    this._sourceConductance = 0.0;
                }

                if (this._sourceConductance > 0.0)
                    this._sourceConductance = 1.0
                                           / this._sourceConductance;
                else
                {
                    this._sourceConductance = 1.0e3; /* mho */
                    SpiceSharpWarning.Warning(this, "Warning: Source conductance reset to 1.0e3 mho.");
                }
            }
            else
            {
                this._sourceConductance = 0.0;
            }

            /* End of Rsd processing */


            Nvtms = ModelTemperature.Vtm * ModelParameters.SjctEmissionCoeff;
            if ((this._aseff <= 0.0) && (this._pseff <= 0.0))
            {
                SourceSatCurrent = 0.0; /* v4.7 */
                /* SourceSatCurrent = 1.0e-14; */
            }
            else
            {
                SourceSatCurrent = this._aseff * ModelTemperature.SjctTempSatCurDensity
                                 + this._pseff * ModelTemperature.SjctSidewallTempSatCurDensity
                                 + Param.BSIM4weffCJ * Parameters.Nf
                                 * ModelTemperature.SjctGateSidewallTempSatCurDensity;
            }
            if (SourceSatCurrent > 0.0)
            {
                switch (ModelParameters.DioMod)
                {
                    case 0:
                        if ((ModelParameters.Bvs / Nvtms) > EXP_THRESHOLD)
                            this._xExpBVS = ModelParameters.Xjbvs * MIN_EXP;
                        else
                            this._xExpBVS = ModelParameters.Xjbvs * Math.Exp(-ModelParameters.Bvs / Nvtms);
                        break;
                    case 1:
                        BSIM4DioIjthVjmEval(Nvtms, ModelParameters.Ijthsfwd, SourceSatCurrent,
                                            0.0, out (this._vjsmFwd));
                        this._iVjsmFwd = SourceSatCurrent * Math.Exp(this._vjsmFwd / Nvtms);
                        break;
                    case 2:
                        if ((ModelParameters.Bvs / Nvtms) > EXP_THRESHOLD)
                        {
                            this._xExpBVS = ModelParameters.Xjbvs * MIN_EXP;
                            tmp = MIN_EXP;
                        }
                        else
                        {
                            this._xExpBVS = Math.Exp(-ModelParameters.Bvs / Nvtms);
                            tmp = this._xExpBVS;
                            this._xExpBVS *= ModelParameters.Xjbvs;
                        }

                        BSIM4DioIjthVjmEval(Nvtms, ModelParameters.Ijthsfwd, SourceSatCurrent,
                                                   this._xExpBVS, out (this._vjsmFwd));
                        T0 = Math.Exp(this._vjsmFwd / Nvtms);
                        this._iVjsmFwd = SourceSatCurrent * (T0 - this._xExpBVS / T0
                                              + this._xExpBVS - 1.0);
                        this._sslpFwd = SourceSatCurrent
                                             * (T0 + this._xExpBVS / T0) / Nvtms;

                        T2 = ModelParameters.Ijthsrev / SourceSatCurrent;
                        if (T2 < 1.0)
                        {
                            T2 = 10.0;
                            SpiceSharpWarning.Warning(this, "Warning: ijthsrev too small and set to 10 times IsbSat.");
                        }
                        this._vjsmRev = -ModelParameters.Bvs
                                           - Nvtms * Math.Log((T2 - 1.0) / ModelParameters.Xjbvs);
                        T1 = ModelParameters.Xjbvs * Math.Exp(-(ModelParameters.Bvs
                           + this._vjsmRev) / Nvtms);
                        this._iVjsmRev = SourceSatCurrent * (1.0 + T1);
                        this._sslpRev = -SourceSatCurrent * T1 / Nvtms;
                        break;
                    default:
                        SpiceSharpWarning.Warning(this, "Specified dioMod = {0} not matched".FormatString(ModelParameters.DioMod));
                        break;
                }
            }

            Nvtmd = ModelTemperature.Vtm * ModelParameters.DjctEmissionCoeff;
            if ((this._adeff <= 0.0) && (this._pdeff <= 0.0))
            {  /* DrainSatCurrent = 1.0e-14;         v4.7 */
                DrainSatCurrent = 0.0;
            }
            else
            {
                DrainSatCurrent = this._adeff * ModelTemperature.DjctTempSatCurDensity
                                + this._pdeff * ModelTemperature.DjctSidewallTempSatCurDensity
                                + Param.BSIM4weffCJ * Parameters.Nf
                                * ModelTemperature.DjctGateSidewallTempSatCurDensity;
            }
            if (DrainSatCurrent > 0.0)
            {
                switch (ModelParameters.DioMod)
                {
                    case 0:
                        if ((ModelParameters.Bvd / Nvtmd) > EXP_THRESHOLD)
                            this._xExpBVD = ModelParameters.Xjbvd * MIN_EXP;
                        else
                            this._xExpBVD = ModelParameters.Xjbvd * Math.Exp(-ModelParameters.Bvd / Nvtmd);
                        break;
                    case 1:
                        BSIM4DioIjthVjmEval(Nvtmd, ModelParameters.Ijthdfwd, DrainSatCurrent,
                                            0.0, out this._vjdmFwd);
                        this._iVjdmFwd = DrainSatCurrent * Math.Exp(this._vjdmFwd / Nvtmd);
                        break;
                    case 2:
                        if ((ModelParameters.Bvd / Nvtmd) > EXP_THRESHOLD)
                        {
                            this._xExpBVD = ModelParameters.Xjbvd * MIN_EXP;
                            tmp = MIN_EXP;
                        }
                        else
                        {
                            this._xExpBVD = Math.Exp(-ModelParameters.Bvd / Nvtmd);
                            tmp = this._xExpBVD;
                            this._xExpBVD *= ModelParameters.Xjbvd;
                        }

                        BSIM4DioIjthVjmEval(Nvtmd, ModelParameters.Ijthdfwd, DrainSatCurrent,
                                            this._xExpBVD, out (this._vjdmFwd));
                        T0 = Math.Exp(this._vjdmFwd / Nvtmd);
                        this._iVjdmFwd = DrainSatCurrent * (T0 - this._xExpBVD / T0
                                            + this._xExpBVD - 1.0);
                        this._dslpFwd = DrainSatCurrent
                                             * (T0 + this._xExpBVD / T0) / Nvtmd;

                        T2 = ModelParameters.Ijthdrev / DrainSatCurrent;
                        if (T2 < 1.0)
                        {
                            T2 = 10.0;
                            SpiceSharpWarning.Warning(this, "Warning: ijthdrev too small and set to 10 times IdbSat.");
                        }
                        this._vjdmRev = -ModelParameters.Bvd
                                           - Nvtmd * Math.Log((T2 - 1.0) / ModelParameters.Xjbvd); /* bugfix */
                        T1 = ModelParameters.Xjbvd * Math.Exp(-(ModelParameters.Bvd
                           + this._vjdmRev) / Nvtmd);
                        this._iVjdmRev = DrainSatCurrent * (1.0 + T1);
                        this._dslpRev = -DrainSatCurrent * T1 / Nvtmd;
                        break;
                    default:
                        SpiceSharpWarning.Warning(this, "Specified dioMod = {0} not matched".FormatString(ModelParameters.DioMod));
                        break;
                }
            }

            /* GEDL current reverse bias */
            T0 = (TRatio - 1.0);
            /*
            // Moved to model temperature - 20220724 - Sven Boulanger
            ModelTemperature.Njtsstemp = ModelParameters.Njts * (1.0 + ModelParameters.Tnjts * T0);
            ModelTemperature.Njtsswstemp = ModelParameters.Njtssw * (1.0 + ModelParameters.Tnjtssw * T0);
            ModelTemperature.Njtsswgstemp = ModelParameters.Njtsswg * (1.0 + ModelParameters.Tnjtsswg * T0);
            ModelTemperature.Njtsdtemp = ModelParameters.Njtsd * (1.0 + ModelParameters.Tnjtsd * T0);
            ModelTemperature.Njtsswdtemp = ModelParameters.Njtsswd * (1.0 + ModelParameters.Tnjtsswd * T0);
            ModelTemperature.Njtsswgdtemp = ModelParameters.Njtsswgd * (1.0 + ModelParameters.Tnjtsswgd * T0);
            */
            T7 = Eg0 / ModelTemperature.Vtm * T0;
            T9 = ModelParameters.Xtss * T7;
            T1 = DEXP(T9);
            T9 = ModelParameters.Xtsd * T7;
            T2 = DEXP(T9);
            T9 = ModelParameters.Xtssws * T7;
            T3 = DEXP(T9);
            T9 = ModelParameters.Xtsswd * T7;
            T4 = DEXP(T9);
            T9 = ModelParameters.Xtsswgs * T7;
            T5 = DEXP(T9);
            T9 = ModelParameters.Xtsswgd * T7;
            T6 = DEXP(T9);
            /*IBM TAT*/
            if (ModelParameters.Jtweff < 0.0)
            {
                ModelParameters.Jtweff = 0.0;
                SpiceSharpWarning.Warning(this, "TAT width dependence effect is negative. Jtweff is clamped to zero.");
            }
            T11 = Math.Sqrt(ModelParameters.Jtweff / Param.BSIM4weffCJ) + 1.0;

            T10 = Param.BSIM4weffCJ * Parameters.Nf;
            this._sjctTempRevSatCur = T1 * this._aseff * ModelParameters.Jtss;
            this._djctTempRevSatCur = T2 * this._adeff * ModelParameters.Jtsd;
            this._sswTempRevSatCur = T3 * this._pseff * ModelParameters.Jtssws;
            this._dswTempRevSatCur = T4 * this._pdeff * ModelParameters.Jtsswd;
            this._sswgTempRevSatCur = T5 * T10 * T11 * ModelParameters.Jtsswgs;
            this._dswgTempRevSatCur = T6 * T10 * T11 * ModelParameters.Jtsswgd;

            if (ModelParameters.MtrlMod.Value != 0 && ModelParameters.MtrlCompatMod.Value == 0)
            {
                /* Calculate TOXP from EOT */
                /* Calculate Vgs_eff @ Vgs = VDD with Poly Depletion Effect */
                Vtm0eot = KboQ * ModelParameters.Tempeot;
                Vtmeot = Vtm0eot;
                vbieot = Vtm0eot * Math.Log(Param.BSIM4nsd
                               * Param.BSIM4ndep / (ni * ni));
                phieot = Vtm0eot * Math.Log(Param.BSIM4ndep / ni)
                               + Param.BSIM4phin + 0.4;
                tmp2 = this._vfb + phieot;
                vddeot = ModelParameters.Type * ModelParameters.Vddeot;
                T0 = ModelParameters.Epsrgate * EPS0;
                if ((Param.BSIM4ngate > 1.0e18) && (Param.BSIM4ngate < 1.0e25)
                    && (vddeot > tmp2) && (T0 != 0))
                {
                    T1 = 1.0e6 * Constants.Charge * T0 * Param.BSIM4ngate /
                      (ModelTemperature.Coxe * ModelTemperature.Coxe);
                    T8 = vddeot - tmp2;
                    T4 = Math.Sqrt(1.0 + 2.0 * T8 / T1);
                    T2 = 2.0 * T8 / (T4 + 1.0);
                    T3 = 0.5 * T2 * T2 / T1;
                    T7 = 1.12 - T3 - 0.05;
                    T6 = Math.Sqrt(T7 * T7 + 0.224);
                    T5 = 1.12 - 0.5 * (T7 + T6);
                    Vgs_eff = vddeot - T5;
                }
                else
                    Vgs_eff = vddeot;

                /* Calculate Vth @ Vds=Vbs=0 */

                V0 = vbieot - phieot;
                lt1 = ModelTemperature.Factor1 * Param.BSIM4sqrtXdep0;
                ltw = lt1;
                T0 = Param.BSIM4dvt1 * ModelParameters.Leffeot / lt1;
                if (T0 < EXP_THRESHOLD)
                {
                    T1 = Math.Exp(T0);
                    T2 = T1 - 1.0;
                    T3 = T2 * T2;
                    T4 = T3 + 2.0 * T1 * MIN_EXP;
                    Theta0 = T1 / T4;
                }
                else
                    Theta0 = 1.0 / (MAX_EXP - 2.0);
                Delt_vth = Param.BSIM4dvt0 * Theta0 * V0;
                T0 = Param.BSIM4dvt1w * ModelParameters.Weffeot * ModelParameters.Leffeot / ltw;
                if (T0 < EXP_THRESHOLD)
                {
                    T1 = Math.Exp(T0);
                    T2 = T1 - 1.0;
                    T3 = T2 * T2;
                    T4 = T3 + 2.0 * T1 * MIN_EXP;
                    T5 = T1 / T4;
                }
                else
                    T5 = 1.0 / (MAX_EXP - 2.0); /* 3.0 * MIN_EXP omitted */
                T2 = Param.BSIM4dvt0w * T5 * V0;
                TempRatioeot = ModelParameters.Tempeot / ModelParameters.Tnom - 1.0;
                T0 = Math.Sqrt(1.0 + Param.BSIM4lpe0 / ModelParameters.Leffeot);
                T1 = Param.BSIM4k1ox * (T0 - 1.0) * Math.Sqrt(phieot)
                  + (Param.BSIM4kt1 + Param.BSIM4kt1l / ModelParameters.Leffeot) * TempRatioeot;
                Vth_NarrowW = toxe * phieot
                  / (ModelParameters.Weffeot + Param.BSIM4w0);
                Lpe_Vb = Math.Sqrt(1.0 + Param.BSIM4lpeb / ModelParameters.Leffeot);
                Vth = ModelParameters.Type * this._vth0 +
                  (Param.BSIM4k1ox - Param.BSIM4k1) * Math.Sqrt(phieot) * Lpe_Vb
                  - Delt_vth - T2 + Param.BSIM4k3 * Vth_NarrowW + T1;

                /* Calculate n */
                tmp1 = epssub / Param.BSIM4Xdep0;
                tmp2 = Param.BSIM4nfactor * tmp1;
                tmp3 = (tmp2 + Param.BSIM4cdsc * Theta0 + Param.BSIM4cit) / ModelTemperature.Coxe;
                if (tmp3 >= -0.5)
                    n = 1.0 + tmp3;
                else
                {
                    T0 = 1.0 / (3.0 + 8.0 * tmp3);
                    n = (1.0 + 3.0 * tmp3) * T0;
                }


                /* Vth correction for Pocket implant */
                if (Param.BSIM4dvtp0 > 0.0)
                {
                    T3 = ModelParameters.Leffeot + Param.BSIM4dvtp0 * 2.0;
                    if (ModelParameters.TempMod < 2)
                        T4 = Vtmeot * Math.Log(ModelParameters.Leffeot / T3);
                    else
                        T4 = Vtm0eot * Math.Log(ModelParameters.Leffeot / T3);
                    Vth -= n * T4;
                }
                Vgsteff = Vgs_eff - Vth;
                /* calculating Toxp */
                T3 = ModelParameters.Type * this._vth0
            - this._vfb - phieot;
                T4 = T3 + T3;
                T5 = 2.5 * T3;

                vtfbphi2eot = 4.0 * T3;
                if (vtfbphi2eot < 0.0)
                    vtfbphi2eot = 0.0;


                niter = 0;
                toxpf = toxe;
                do
                {
                    toxpi = toxpf;
                    tmp2 = 2.0e8 * toxpf;
                    T0 = (Vgsteff + vtfbphi2eot) / tmp2;
                    T1 = 1.0 + Math.Exp(ModelParameters.Bdos * 0.7 * Math.Log(T0));
                    Tcen = ModelParameters.Ados * 1.9e-9 / T1;
                    toxpf = toxe - epsrox / ModelParameters.Epsrsub * Tcen;
                    niter++;
                } while ((niter <= 4) && (Math.Abs(toxpf - toxpi) > 1e-12));
                this._toxp = toxpf;
                this._coxp = epsrox * EPS0 / this._toxp;
            }
            else
            {
                this._toxp = ModelParameters.Toxp;
                this._coxp = ModelTemperature.Coxp;
            }

            if (BSIM4checkModel())
            {
                throw new SpiceSharpException("Detected errors during BSIM4.8.1 parameter checking for model {0} of device instance {1}".FormatString(ModelTemperature.Name, Name));
            }

            /* End of temperature */
        }

        private static double DEXP(double a)
        {
            if (a > EXP_THRESHOLD)
                return MAX_EXP * (1.0 + (a) - EXP_THRESHOLD);
            else if (a < -EXP_THRESHOLD)
                return MIN_EXP;
            else
                return Math.Exp(a);
        }

        private static void BSIM4DioIjthVjmEval(double Nvtm, double Ijth, double Isb, double XExpBV, out double Vjm)
        {
            double Tb, Tc, EVjmovNv;
            Tc = XExpBV;
            Tb = 1.0 + Ijth / Isb - Tc;
            EVjmovNv = 0.5 * (Tb + Math.Sqrt(Tb * Tb + 4.0 * Tc));
            Vjm = Nvtm * Math.Log(EVjmovNv);
        }

        private bool BSIM4checkModel()
        {
            List<string> words = new List<string>();
            int Fatal_Flag = 0;

            if (ModelParameters.Version != "4.8.1" && ModelParameters.Version != "4.81" && ModelParameters.Version != "4.8")
            {
                SpiceSharpWarning.Warning(this, "Warning: This model supports BSIM4 version 4.8");
                SpiceSharpWarning.Warning(this, "You specified a wrong version number. Working now with BSIM4.8.1");
                words.Add("Warning: This model supports BSIM4 version 4.8");
                words.Add("You specified a wrong version number. Working now with BSIM4.8.1.");
            }

            if ((Parameters.RgateMod.Value == 2) || (Parameters.RgateMod.Value == 3))
            {
                if ((Parameters.TrnqsMod.Value == 1) || (Parameters.AcnqsMod.Value == 1))
                {
                    words.Add("Warning: You've selected both Rg and charge deficit NQS; select one only.");
                }
            }

            if (ModelParameters.Toxe <= 0.0)
            {
                words.Add("Fatal: Toxe = {0:g} is not positive.".FormatString(ModelParameters.Toxe));
                Fatal_Flag = 1;
            }

            if (ModelParameters.Toxp <= 0.0)
            {
                words.Add("Fatal: Toxp = {0:g} is not positive.".FormatString(ModelParameters.Toxp));
                Fatal_Flag = 1;
            }

            if (ModelParameters.Eot <= 0.0)
            {
                words.Add("Fatal: EOT = {0:g} is not positive.".FormatString(ModelParameters.Eot));
                Fatal_Flag = 1;
            }
            if (ModelParameters.Epsrgate < 0.0)
            {
                words.Add("Fatal: Epsrgate = {0:g} is not positive.".FormatString(ModelParameters.Epsrgate));
                Fatal_Flag = 1;
            }
            if (ModelParameters.Epsrsub < 0.0)
            {
                words.Add("Fatal: Epsrsub = {0:g} is not positive.".FormatString(ModelParameters.Epsrsub));
                Fatal_Flag = 1;
            }
            if (ModelParameters.Easub < 0.0)
            {
                words.Add("Fatal: Easub = {0:g} is not positive.".FormatString(ModelParameters.Easub));
                Fatal_Flag = 1;
            }
            if (ModelParameters.Ni0sub <= 0.0)
            {
                words.Add("Fatal: Easub = {0:g} is not positive.".FormatString(ModelParameters.Ni0sub));
                Fatal_Flag = 1;
            }

            if (ModelParameters.Toxm <= 0.0)
            {
                words.Add("Fatal: Toxm = {0:g} is not positive.".FormatString(ModelParameters.Toxm));
                Fatal_Flag = 1;
            }
            if (ModelParameters.Toxref <= 0.0)
            {
                words.Add("Fatal: Toxref = {0:g} is not positive.".FormatString(ModelParameters.Toxref));
                Fatal_Flag = 1;
            }

            if (Param.BSIM4lpe0 < -Param.BSIM4leff)
            {
                words.Add("Fatal: Lpe0 = {0:g} is less than -Leff.".FormatString(Param.BSIM4lpe0));
                Fatal_Flag = 1;
            }
            if (ModelParameters.Lintnoi > Param.BSIM4leff / 2)
            {
                words.Add("Fatal: Lintnoi = {0:g} is too large - Leff for noise is negative.".FormatString(ModelParameters.Lintnoi));
                Fatal_Flag = 1;
            }
            if (Param.BSIM4lpeb < -Param.BSIM4leff)
            {
                words.Add("Fatal: Lpeb = {0:g} is less than -Leff.".FormatString(Param.BSIM4lpeb));
                Fatal_Flag = 1;
            }
            if (Param.BSIM4ndep <= 0.0)
            {
                words.Add("Fatal: Ndep = {0:g} is not positive.".FormatString(Param.BSIM4ndep));
                Fatal_Flag = 1;
            }
            if (Param.BSIM4phi <= 0.0)
            {
                words.Add("Fatal: Phi = {0:g} is not positive. Please check Phin and Ndep".FormatString(Param.BSIM4phi));
                words.Add("	   Phin = {0:g}  Ndep = %g ".FormatString(Param.BSIM4phin, Param.BSIM4ndep));
                Fatal_Flag = 1;
            }
            if (Param.BSIM4nsub <= 0.0)
            {
                words.Add("Fatal: Nsub = {0:g} is not positive.".FormatString(Param.BSIM4nsub));
                Fatal_Flag = 1;
            }
            if (Param.BSIM4ngate < 0.0)
            {
                words.Add("Fatal: Ngate = {0:g} Ngate is not positive.".FormatString(Param.BSIM4ngate));
                Fatal_Flag = 1;
            }
            if (Param.BSIM4ngate > 1.0e25)
            {
                words.Add("Fatal: Ngate = {0:g} Ngate is too high".FormatString(Param.BSIM4ngate));
                Fatal_Flag = 1;
            }
            if (Param.BSIM4xj <= 0.0)
            {
                words.Add("Fatal: Xj = {0:g} is not positive.".FormatString(Param.BSIM4xj));
                Fatal_Flag = 1;
            }

            if (Param.BSIM4dvt1 < 0.0)
            {
                words.Add("Fatal: Dvt1 = {0:g} is negative.".FormatString(Param.BSIM4dvt1));
                Fatal_Flag = 1;
            }

            if (Param.BSIM4dvt1w < 0.0)
            {
                words.Add("Fatal: Dvt1w = {0:g} is negative.".FormatString(Param.BSIM4dvt1w));
                Fatal_Flag = 1;
            }

            if (Param.BSIM4w0 == -Param.BSIM4weff)
            {
                words.Add("Fatal: (W0 + Weff) = 0 causing divided-by-zero.");
                Fatal_Flag = 1;
            }

            if (Param.BSIM4dsub < 0.0)
            {
                words.Add("Fatal: Dsub = {0:g} is negative.".FormatString(Param.BSIM4dsub));
                Fatal_Flag = 1;
            }
            if (Param.BSIM4b1 == -Param.BSIM4weff)
            {
                words.Add("Fatal: (B1 + Weff) = 0 causing divided-by-zero.");
                Fatal_Flag = 1;
            }
            if (this._u0temp <= 0.0)
            {
                words.Add("Fatal: u0 at current temperature = {0:g} is not positive.".FormatString(this._u0temp));
                Fatal_Flag = 1;
            }

            if (Param.BSIM4delta < 0.0)
            {
                words.Add("Fatal: Delta = {0:g} is less than zero.".FormatString(Param.BSIM4delta));
                Fatal_Flag = 1;
            }

            if (this._vsattemp <= 0.0)
            {
                words.Add("Fatal: Vsat at current temperature = {0:g} is not positive.".FormatString(this._vsattemp));
                Fatal_Flag = 1;
            }

            if (Param.BSIM4pclm <= 0.0)
            {
                words.Add("Fatal: Pclm = {0:g} is not positive.".FormatString(Param.BSIM4pclm));
                Fatal_Flag = 1;
            }

            if (Param.BSIM4drout < 0.0)
            {
                words.Add("Fatal: Drout = {0:g} is negative.".FormatString(Param.BSIM4drout));
                Fatal_Flag = 1;
            }

            if (Parameters.M <= 0.0)
            {
                words.Add("Fatal: multiplier = {0:g} is not positive.".FormatString(Parameters.M));
                Fatal_Flag = 1;
            }
            if (Parameters.Nf < 1.0)
            {
                words.Add("Fatal: Number of finger = {0:g} is smaller than one.".FormatString(Parameters.Nf));
                Fatal_Flag = 1;
            }

            if ((Parameters.Sa > 0.0) && (Parameters.Sb > 0.0) &&
                ((Parameters.Nf == 1.0) || ((Parameters.Nf > 1.0) && (Parameters.Sd > 0.0))))
            {
                if (ModelParameters.Saref <= 0.0)
                {
                    words.Add("Fatal: SAref = {0:g} is not positive.".FormatString(ModelParameters.Saref));
                    Fatal_Flag = 1;
                }
                if (ModelParameters.Sbref <= 0.0)
                {
                    words.Add("Fatal: SBref = {0:g} is not positive.".FormatString(ModelParameters.Sbref));
                    Fatal_Flag = 1;
                }
            }

            if ((Parameters.L + ModelParameters.Xl) <= ModelParameters.Xgl)
            {
                words.Add("Fatal: The parameter xgl must be smaller than Ldrawn+XL.");
                Fatal_Flag = 1;
            }
            if (Parameters.Ngcon < 1.0)
            {
                words.Add("Fatal: The parameter ngcon cannot be smaller than one.");
                Fatal_Flag = 1;
            }
            if ((Parameters.Ngcon != 1.0) && (Parameters.Ngcon != 2.0))
            {
                Parameters.Ngcon = 1.0;
                words.Add("Warning: Ngcon must be equal to one or two; reset to 1.0.");
            }

            if (ModelParameters.Gbmin < 1.0e-20)
            {
                words.Add("Warning: Gbmin = {0:g} is too small.".FormatString(ModelParameters.Gbmin));
            }

            /* Check saturation parameters */
            if (Param.BSIM4fprout < 0.0)
            {
                words.Add("Fatal: fprout = {0:g} is negative.".FormatString(Param.BSIM4fprout));
                Fatal_Flag = 1;
            }
            if (Param.BSIM4pdits < 0.0)
            {
                words.Add("Fatal: pdits = {0:g} is negative.".FormatString(Param.BSIM4pdits));
                Fatal_Flag = 1;
            }
            if (ModelParameters.Pditsl < 0.0)
            {
                words.Add("Fatal: pditsl = {0:g} is negative.".FormatString(ModelParameters.Pditsl));
                Fatal_Flag = 1;
            }

            /* Check gate current parameters */
            if (ModelParameters.IgbMod.Value != 0)
            {
                if (Param.BSIM4nigbinv <= 0.0)
                {
                    words.Add("Fatal: nigbinv = {0:g} is non-positive.".FormatString(Param.BSIM4nigbinv));
                    Fatal_Flag = 1;
                }
                if (Param.BSIM4nigbacc <= 0.0)
                {
                    words.Add("Fatal: nigbacc = {0:g} is non-positive.".FormatString(Param.BSIM4nigbacc));
                    Fatal_Flag = 1;
                }
            }
            if (ModelParameters.IgcMod.Value != 0)
            {
                if (Param.BSIM4nigc <= 0.0)
                {
                    words.Add("Fatal: nigc = {0:g} is non-positive.".FormatString(Param.BSIM4nigc));
                    Fatal_Flag = 1;
                }
                if (Param.BSIM4poxedge <= 0.0)
                {
                    words.Add("Fatal: poxedge = {0:g} is non-positive.".FormatString(Param.BSIM4poxedge));
                    Fatal_Flag = 1;
                }
                if (Param.BSIM4pigcd <= 0.0)
                {
                    words.Add("Fatal: pigcd = {0:g} is non-positive.".FormatString(Param.BSIM4pigcd));
                    Fatal_Flag = 1;
                }
            }

            /* Check capacitance parameters */
            if (Param.BSIM4clc < 0.0)
            {
                words.Add("Fatal: Clc = {0:g} is negative.".FormatString(Param.BSIM4clc));
                Fatal_Flag = 1;
            }

            /* Check overlap capacitance parameters */
            if (Param.BSIM4ckappas < 0.02)
            {
                words.Add("Warning: ckappas = {0:g} is too small.".FormatString(Param.BSIM4ckappas));
                Param.BSIM4ckappas = 0.02;
            }
            if (Param.BSIM4ckappad < 0.02)
            {
                words.Add("Warning: ckappad = {0:g} is too small.".FormatString(Param.BSIM4ckappad));
                Param.BSIM4ckappad = 0.02;
            }

            if (ModelParameters.Vtss < 0.0)
            {
                words.Add("Fatal: Vtss = {0:g} is negative.".FormatString(ModelParameters.Vtss));
                Fatal_Flag = 1;
            }
            if (ModelParameters.Vtsd < 0.0)
            {
                words.Add("Fatal: Vtsd = {0:g} is negative.".FormatString(ModelParameters.Vtsd));
                Fatal_Flag = 1;
            }
            if (ModelParameters.Vtssws < 0.0)
            {
                words.Add("Fatal: Vtssws = {0:g} is negative.".FormatString(ModelParameters.Vtssws));
                Fatal_Flag = 1;
            }
            if (ModelParameters.Vtsswd < 0.0)
            {
                words.Add("Fatal: Vtsswd = {0:g} is negative.".FormatString(ModelParameters.Vtsswd));
                Fatal_Flag = 1;
            }
            if (ModelParameters.Vtsswgs < 0.0)
            {
                words.Add("Fatal: Vtsswgs = {0:g} is negative.".FormatString(ModelParameters.Vtsswgs));
                Fatal_Flag = 1;
            }
            if (ModelParameters.Vtsswgd < 0.0)
            {
                words.Add("Fatal: Vtsswgd = {0:g} is negative.".FormatString(ModelParameters.Vtsswgd));
                Fatal_Flag = 1;
            }


            if (ModelParameters.ParamChk == 1)
            {
                /* Check L and W parameters */
                if (Param.BSIM4leff <= 1.0e-9)
                {
                    words.Add("Warning: Leff = {0:g} <= 1.0e-9. Recommended Leff >= 1e-8 ".FormatString(Param.BSIM4leff));
                }

                if (Param.BSIM4leffCV <= 1.0e-9)
                {
                    words.Add("Warning: Leff for CV = {0:g} <= 1.0e-9. Recommended LeffCV >=1e-8 ".FormatString(Param.BSIM4leffCV));
                }

                if (Param.BSIM4weff <= 1.0e-9)
                {
                    words.Add("Warning: Weff = {0:g} <= 1.0e-9. Recommended Weff >=1e-7 ".FormatString(Param.BSIM4weff));
                }

                if (Param.BSIM4weffCV <= 1.0e-9)
                {
                    words.Add("Warning: Weff for CV = {0:g} <= 1.0e-9. Recommended WeffCV >= 1e-7 ".FormatString(Param.BSIM4weffCV));
                }

                /* Check threshold voltage parameters */
                if (ModelParameters.Toxe < 1.0e-10)
                {
                    words.Add("Warning: Toxe = {0:g} is less than 1A. Recommended Toxe >= 5A".FormatString(ModelParameters.Toxe));
                }
                if (ModelParameters.Toxp < 1.0e-10)
                {
                    words.Add("Warning: Toxp = {0:g} is less than 1A. Recommended Toxp >= 5A".FormatString(ModelParameters.Toxp));
                }
                if (ModelParameters.Toxm < 1.0e-10)
                {
                    words.Add("Warning: Toxm = {0:g} is less than 1A. Recommended Toxm >= 5A".FormatString(ModelParameters.Toxm));
                }

                if (Param.BSIM4ndep <= 1.0e12)
                {
                    words.Add("Warning: Ndep = {0:g} may be too small.".FormatString(Param.BSIM4ndep));
                }
                else if (Param.BSIM4ndep >= 1.0e21)
                {
                    words.Add("Warning: Ndep = {0:g} may be too large.".FormatString(Param.BSIM4ndep));
                }

                if (Param.BSIM4nsub <= 1.0e14)
                {
                    words.Add("Warning: Nsub = {0:g} may be too small.".FormatString(Param.BSIM4nsub));
                }
                else if (Param.BSIM4nsub >= 1.0e21)
                {
                    words.Add("Warning: Nsub = {0:g} may be too large.".FormatString(Param.BSIM4nsub));
                }

                if ((Param.BSIM4ngate > 0.0) &&
                    (Param.BSIM4ngate <= 1.0e18))
                {
                    words.Add("Warning: Ngate = {0:g} is less than 1.E18cm^-3.".FormatString(Param.BSIM4ngate));
                }

                if (Param.BSIM4dvt0 < 0.0)
                {
                    words.Add("Warning: Dvt0 = {0:g} is negative.".FormatString(Param.BSIM4dvt0));
                }

                if (Math.Abs(1.0e-8 / (Param.BSIM4w0 + Param.BSIM4weff)) > 10.0)
                {
                    words.Add("Warning: (W0 + Weff) may be too small.");
                }

                /* Check subthreshold parameters */
                if (Param.BSIM4nfactor < 0.0)
                {
                    words.Add("Warning: Nfactor = {0:g} is negative.".FormatString(Param.BSIM4nfactor));
                }
                if (Param.BSIM4cdsc < 0.0)
                {
                    words.Add("Warning: Cdsc = {0:g} is negative.".FormatString(Param.BSIM4cdsc));
                }
                if (Param.BSIM4cdscd < 0.0)
                {
                    words.Add("Warning: Cdscd = {0:g} is negative.".FormatString(Param.BSIM4cdscd));
                }
                /* Check DIBL parameters */
                if (this._eta0 < 0.0)
                {
                    words.Add("Warning: Eta0 = {0:g} is negative.".FormatString(this._eta0));
                }

                /* Check Abulk parameters */
                if (Math.Abs(1.0e-8 / (Param.BSIM4b1 + Param.BSIM4weff)) > 10.0)
                {
                    words.Add("Warning: (B1 + Weff) may be too small.");
                }


                /* Check Saturation parameters */
                if (Param.BSIM4a2 < 0.01)
                {
                    words.Add("Warning: A2 = {0:g} is too small. Set to 0.01.".FormatString(Param.BSIM4a2));
                    Param.BSIM4a2 = 0.01;
                }
                else if (Param.BSIM4a2 > 1.0)
                {
                    words.Add("Warning: A2 = {0:g} is larger than 1. A2 is set to 1 and A1 is set to 0.".FormatString(Param.BSIM4a2));
                    Param.BSIM4a2 = 1.0;
                    Param.BSIM4a1 = 0.0;
                }

                if (Param.BSIM4prwg < 0.0)
                {
                    words.Add("Warning: Prwg = {0:g} is negative. Set to zero.".FormatString(Param.BSIM4prwg));
                    Param.BSIM4prwg = 0.0;
                }

                if (Param.BSIM4rdsw < 0.0)
                {
                    words.Add("Warning: Rdsw = {0:g} is negative. Set to zero.".FormatString(Param.BSIM4rdsw));
                    Param.BSIM4rdsw = 0.0;
                    Param.BSIM4rds0 = 0.0;
                }

                if (Param.BSIM4rds0 < 0.0)
                {
                    words.Add("Warning: Rds at current temperature = {0:g} is negative. Set to zero.".FormatString(Param.BSIM4rds0));
                    Param.BSIM4rds0 = 0.0;
                }

                if (Param.BSIM4rdswmin < 0.0)
                {
                    words.Add("Warning: Rdswmin at current temperature = {0:g} is negative. Set to zero.".FormatString(Param.BSIM4rdswmin));
                    Param.BSIM4rdswmin = 0.0;
                }

                if (Param.BSIM4pscbe2 <= 0.0)
                {
                    words.Add("Warning: Pscbe2 = {0:g} is not positive.".FormatString(Param.BSIM4pscbe2));
                }

                if (Param.BSIM4vsattemp < 1.0e3)
                {
                    words.Add("Warning: Vsat at current temperature = {0:g} may be too small.".FormatString(Param.BSIM4vsattemp));
                }

                if ((ModelParameters.Lambda.Given) && (Param.BSIM4lambda > 0.0))
                {
                    if (Param.BSIM4lambda > 1.0e-9)
                    {
                        words.Add("Warning: Lambda = {0:g} may be too large.".FormatString(Param.BSIM4lambda));
                    }
                }

                if ((ModelParameters.Vtl.Given) && (Param.BSIM4vtl > 0.0))
                {
                    if (Param.BSIM4vtl < 6.0e4)
                    {
                        words.Add("Warning: Thermal velocity vtl = {0:g} may be too small.".FormatString(Param.BSIM4vtl));
                    }

                    if (Param.BSIM4xn < 3.0)
                    {
                        words.Add("Warning: back scattering coeff xn = {0:g} is too small. Reset to 3.0 ".FormatString(Param.BSIM4xn));
                        Param.BSIM4xn = 3.0;
                    }

                    if (ModelParameters.Lc < 0.0)
                    {
                        words.Add("Warning: back scattering coeff lc = {0:g} is too small. Reset to 0.0".FormatString(ModelParameters.Lc));
                        Param.BSIM4lc = 0.0;
                    }
                }

                if (Param.BSIM4pdibl1 < 0.0)
                {
                    words.Add("Warning: Pdibl1 = {0:g} is negative.".FormatString(Param.BSIM4pdibl1));
                }
            }

            if (Param.BSIM4pdibl2 < 0.0)
            {
                words.Add("Warning: Pdibl2 = {0:g} is negative.".FormatString(Param.BSIM4pdibl2));
            }

            /* Check stress effect parameters */
            if ((Parameters.Sa > 0.0) && (Parameters.Sb > 0.0) &&
                ((Parameters.Nf == 1.0) || ((Parameters.Nf > 1.0) && (Parameters.Sd > 0.0))))
            {
                if (ModelParameters.Lodk2 <= 0.0)
                {
                    words.Add("Warning: LODK2 = {0:g} is not positive.".FormatString(ModelParameters.Lodk2));
                }
                if (ModelParameters.Lodeta0 <= 0.0)
                {
                    words.Add("Warning: LODETA0 = {0:g} is not positive.".FormatString(ModelParameters.Lodeta0));
                }
            }

            /* Check gate resistance parameters */
            if (Parameters.RgateMod.Value == 1)
            {
                if (ModelParameters.Rshg <= 0.0)
                    words.Add("Warning: rshg should be positive for rgateMod = 1.");
            }
            else if (Parameters.RgateMod.Value == 2)
            {
                if (ModelParameters.Rshg <= 0.0)
                    words.Add("Warning: rshg <= 0.0 for rgateMod = 2.");
                else if (Param.BSIM4xrcrg1 <= 0.0)
                    words.Add("Warning: xrcrg1 <= 0.0 for rgateMod = 2.");
            }
            if (Parameters.RgateMod.Value == 3)
            {
                if (ModelParameters.Rshg <= 0.0)
                    words.Add("Warning: rshg should be positive for rgateMod = 3.");
                else if (Param.BSIM4xrcrg1 <= 0.0)
                    words.Add("Warning: xrcrg1 should be positive for rgateMod = 3.");
            }

            /* Check body resistance parameters */
            if (ModelParameters.Rbps0 <= 0.0)
            {
                words.Add("Fatal: RBPS0 = {0:g} is not positive.".FormatString(ModelParameters.Rbps0));
                Fatal_Flag = 1;
            }
            if (ModelParameters.Rbpd0 <= 0.0)
            {
                words.Add("Fatal: RBPD0 = {0:g} is not positive.".FormatString(ModelParameters.Rbpd0));
                Fatal_Flag = 1;
            }
            if (ModelParameters.Rbpbx0 <= 0.0)
            {
                words.Add("Fatal: RBPBX0 = {0:g} is not positive.".FormatString(ModelParameters.Rbpbx0));
                Fatal_Flag = 1;
            }
            if (ModelParameters.Rbpby0 <= 0.0)
            {
                words.Add("Fatal: RBPBY0 = {0:g} is not positive.".FormatString(ModelParameters.Rbpby0));
                Fatal_Flag = 1;
            }
            if (ModelParameters.Rbdbx0 <= 0.0)
            {
                words.Add("Fatal: RBDBX0 = {0:g} is not positive.".FormatString(ModelParameters.Rbdbx0));
                Fatal_Flag = 1;
            }
            if (ModelParameters.Rbdby0 <= 0.0)
            {
                words.Add("Fatal: RBDBY0 = {0:g} is not positive.".FormatString(ModelParameters.Rbdby0));
                Fatal_Flag = 1;
            }
            if (ModelParameters.Rbsbx0 <= 0.0)
            {
                words.Add("Fatal: RBSBX0 = {0:g} is not positive.".FormatString(ModelParameters.Rbsbx0));
                Fatal_Flag = 1;
            }
            if (ModelParameters.Rbsby0 <= 0.0)
            {
                words.Add("Fatal: RBSBY0 = {0:g} is not positive.".FormatString(ModelParameters.Rbsby0));
                Fatal_Flag = 1;
            }

            /* Check capacitance parameters */
            if (Param.BSIM4noff < 0.1)
            {
                words.Add("Warning: Noff = {0:g} is too small.".FormatString(Param.BSIM4noff));
            }

            if (Param.BSIM4voffcv < -0.5)
            {
                words.Add("Warning: Voffcv = {0:g} is too small.".FormatString(Param.BSIM4voffcv));
            }

            if (Param.BSIM4moin < 5.0)
            {
                words.Add("Warning: Moin = {0:g} is too small.".FormatString(Param.BSIM4moin));
            }
            if (Param.BSIM4moin > 25.0)
            {
                words.Add("Warning: Moin = {0:g} is too large.".FormatString(Param.BSIM4moin));
            }
            if (ModelParameters.CapMod.Value == 2)
            {
                if (Param.BSIM4acde < 0.1)
                {
                    words.Add("Warning: Acde = {0:g} is too small.".FormatString(Param.BSIM4acde));
                }
                if (Param.BSIM4acde > 1.6)
                {
                    words.Add("Warning: Acde = {0:g} is too large.".FormatString(Param.BSIM4acde));
                }
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
            if (ModelParameters.TnoiMod.Value == 1)
            {
                words.Add("Warning: TNOIMOD=1 is not supported and may be removed from future version.");
            }

            if (ModelParameters.Version != "4.8.1" && ModelParameters.Version != "4.81")
            {
                /* v4.7 */
                if (ModelParameters.TnoiMod.Value == 1 || ModelParameters.TnoiMod.Value == 2)
                {
                    if (ModelParameters.Tnoia < 0.0)
                    {
                        words.Add("Warning: tnoia = {0:g} is negative. Set to zero.".FormatString(ModelParameters.Tnoia));
                        ModelParameters.Tnoia = 0.0;
                    }
                    if (ModelParameters.Tnoib < 0.0)
                    {
                        words.Add("Warning: tnoib = {0:g} is negative. Set to zero.".FormatString(ModelParameters.Tnoib));
                        ModelParameters.Tnoib = 0.0;
                    }

                    if (ModelParameters.Rnoia < 0.0)
                    {
                        words.Add("Warning: rnoia = {0:g} is negative. Set to zero.".FormatString(ModelParameters.Rnoia));
                        ModelParameters.Rnoia = 0.0;
                    }
                    if (ModelParameters.Rnoib < 0.0)
                    {
                        words.Add("Warning: rnoib = {0:g} is negative. Set to zero.".FormatString(ModelParameters.Rnoib));
                        ModelParameters.Rnoib = 0.0;
                    }
                }

                /* v4.7 */
                if (ModelParameters.TnoiMod.Value == 2)
                {
                    if (ModelParameters.Tnoic < 0.0)
                    {
                        words.Add("Warning: tnoic = {0:g} is negative. Set to zero.".FormatString(ModelParameters.Tnoic));
                        ModelParameters.Tnoic = 0.0;
                    }
                    if (ModelParameters.Rnoic < 0.0)
                    {
                        words.Add("Warning: rnoic = {0:g} is negative. Set to zero.".FormatString(ModelParameters.Rnoic));
                        ModelParameters.Rnoic = 0.0;
                    }
                }
            }
            else
            {
                if (ModelParameters.TnoiMod.Value == 1)
                {
                    if (ModelParameters.Tnoia < 0.0)
                    {
                        words.Add("Warning: tnoia = {0:g} is negative. Set to zero.".FormatString(ModelParameters.Tnoia));
                        ModelParameters.Tnoia = 0.0;
                    }
                    if (ModelParameters.Tnoib < 0.0)
                    {
                        words.Add("Warning: tnoib = {0:g} is negative. Set to zero.".FormatString(ModelParameters.Tnoib));
                        ModelParameters.Tnoib = 0.0;
                    }
                    if (ModelParameters.Rnoia < 0.0)
                    {
                        words.Add("Warning: rnoia = {0:g} is negative. Set to zero.".FormatString(ModelParameters.Rnoia));
                        ModelParameters.Rnoia = 0.0;
                    }
                    if (ModelParameters.Rnoib < 0.0)
                    {
                        words.Add("Warning: rnoib = {0:g} is negative. Set to zero.".FormatString(ModelParameters.Rnoib));
                        ModelParameters.Rnoib = 0.0;
                    }
                }
            }

            /* Limits of Njs and Njd modified in BSIM4.7 */
            if (ModelParameters.SjctEmissionCoeff < 0.1)
            {
                words.Add("Warning: Njs = {0:g} is less than 0.1. Setting Njs to 0.1.".FormatString(ModelParameters.SjctEmissionCoeff));
                ModelParameters.SjctEmissionCoeff = 0.1;
            }
            else if (ModelParameters.SjctEmissionCoeff < 0.7)
            {
                words.Add("Warning: Njs = {0:g} is less than 0.7.".FormatString(ModelParameters.SjctEmissionCoeff));
            }
            if (ModelParameters.DjctEmissionCoeff < 0.1)
            {
                words.Add("Warning: Njd = {0:g} is less than 0.1. Setting Njd to 0.1.".FormatString(ModelParameters.DjctEmissionCoeff));
                ModelParameters.DjctEmissionCoeff = 0.1;
            }
            else if (ModelParameters.DjctEmissionCoeff < 0.7)
            {
                words.Add("Warning: Njd = {0:g} is less than 0.7.".FormatString(ModelParameters.DjctEmissionCoeff));
            }

            if (ModelTemperature.Njtsstemp < 0.0)
            {
                words.Add("Warning: Njts = {0:g} is negative at temperature = %g.".FormatString(ModelTemperature.Njtsstemp, _temperature.Temperature));
            }
            if (ModelTemperature.Njtsswstemp < 0.0)
            {
                words.Add("Warning: Njtssw = {0:g} is negative at temperature = %g.".FormatString(ModelTemperature.Njtsswstemp, _temperature.Temperature));
            }
            if (ModelTemperature.Njtsswgstemp < 0.0)
            {
                words.Add("Warning: Njtsswg = {0:g} is negative at temperature = %g.".FormatString(ModelTemperature.Njtsswgstemp, _temperature.Temperature));
            }

            if (ModelParameters.Njtsd.Given && ModelTemperature.Njtsdtemp < 0.0)
            {
                words.Add("Warning: Njtsd = {0:g} is negative at temperature = %g.".FormatString(ModelTemperature.Njtsdtemp, _temperature.Temperature));
            }
            if (ModelParameters.Njtsswd.Given && ModelTemperature.Njtsswdtemp < 0.0)
            {
                words.Add("Warning: Njtsswd = {0:g} is negative at temperature = %g.".FormatString(ModelTemperature.Njtsswdtemp, _temperature.Temperature));
            }
            if (ModelParameters.Njtsswgd.Given && ModelTemperature.Njtsswgdtemp < 0.0)
            {
                words.Add("Warning: Njtsswgd = {0:g} is negative at temperature = %g.".FormatString(ModelTemperature.Njtsswgdtemp, _temperature.Temperature));
            }

            if (ModelParameters.Ntnoi < 0.0)
            {
                words.Add("Warning: ntnoi = {0:g} is negative. Set to zero.".FormatString(ModelParameters.Ntnoi));
                ModelParameters.Ntnoi = 0.0;
            }

            /* diode model */
            if (ModelParameters.SbulkJctBotGradingCoeff >= 0.99)
            {
                words.Add("Warning: MJS = {0:g} is too big. Set to 0.99.".FormatString(ModelParameters.SbulkJctBotGradingCoeff));
                ModelParameters.SbulkJctBotGradingCoeff = 0.99;
            }
            if (ModelParameters.SbulkJctSideGradingCoeff >= 0.99)
            {
                words.Add("Warning: MJSWS = {0:g} is too big. Set to 0.99.".FormatString(ModelParameters.SbulkJctSideGradingCoeff));
                ModelParameters.SbulkJctSideGradingCoeff = 0.99;
            }
            if (ModelParameters.SbulkJctGateSideGradingCoeff >= 0.99)
            {
                words.Add("Warning: MJSWGS = {0:g} is too big. Set to 0.99.".FormatString(ModelParameters.SbulkJctGateSideGradingCoeff));
                ModelParameters.SbulkJctGateSideGradingCoeff = 0.99;
            }

            if (ModelParameters.DbulkJctBotGradingCoeff >= 0.99)
            {
                words.Add("Warning: MJD = {0:g} is too big. Set to 0.99.".FormatString(ModelParameters.DbulkJctBotGradingCoeff));
                ModelParameters.DbulkJctBotGradingCoeff = 0.99;
            }
            if (ModelParameters.DbulkJctSideGradingCoeff >= 0.99)
            {
                words.Add("Warning: MJSWD = {0:g} is too big. Set to 0.99.".FormatString(ModelParameters.DbulkJctSideGradingCoeff));
                ModelParameters.DbulkJctSideGradingCoeff = 0.99;
            }
            if (ModelParameters.DbulkJctGateSideGradingCoeff >= 0.99)
            {
                words.Add("Warning: MJSWGD = {0:g} is too big. Set to 0.99.".FormatString(ModelParameters.DbulkJctGateSideGradingCoeff));
                ModelParameters.DbulkJctGateSideGradingCoeff = 0.99;
            }
            if (ModelParameters.Wpemod == 1)
            {
                if (ModelParameters.Scref <= 0.0)
                {
                    words.Add("Warning: SCREF = {0:g} is not positive. Set to 1e-6.".FormatString(ModelParameters.Scref));
                    ModelParameters.Scref = 1e-6;
                }
                if (Parameters.Sca < 0.0)
                {
                    words.Add("Warning: SCA = {0:g} is negative. Set to 0.0.".FormatString(Parameters.Sca));
                    Parameters.Sca = 0.0;
                }
                if (Parameters.Scb < 0.0)
                {
                    words.Add("Warning: SCB = {0:g} is negative. Set to 0.0.".FormatString(Parameters.Scb));
                    Parameters.Scb = 0.0;
                }
                if (Parameters.Scc < 0.0)
                {
                    words.Add("Warning: SCC = {0:g} is negative. Set to 0.0.".FormatString(Parameters.Scc));
                    Parameters.Scc = 0.0;
                }
                if (Parameters.Sc < 0.0)
                {
                    words.Add("Warning: SC = {0:g} is negative. Set to 0.0.".FormatString(Parameters.Sc));
                    Parameters.Sc = 0.0;
                }
            }

            if (words.Count > 0)
            {
                string path = ModelParameters.CheckPath;
                if (string.IsNullOrWhiteSpace(path))
                    path = "b4v8check.log";
                using (var writer = new StreamWriter(path))
                {
                    foreach (string line in words)
                        writer.WriteLine(line);
                }
            }
            return Fatal_Flag > 0;
        }

        private int BSIM4PAeffGeo(double nf, int geo, int minSD, double Weffcj, double DMCG, double DMCI, double DMDG,
            out double Ps, out double Pd, out double As, out double Ad)
        {
            double T0, T1, T2;
            double ADiso, ADsha, ADmer, ASiso, ASsha, ASmer;
            double PDiso, PDsha, PDmer, PSiso, PSsha, PSmer;
            double nuIntD = 0.0, nuEndD = 0.0, nuIntS = 0.0, nuEndS = 0.0;

            if (geo < 9) /* For geo = 9 and 10, the numbers of S/D diffusions already known */
                BSIM4NumFingerDiff(nf, minSD, out nuIntD, out nuEndD, out nuIntS, out nuEndS);

            T0 = DMCG + DMCI;
            T1 = DMCG + DMCG;
            T2 = DMDG + DMDG;

            PSiso = PDiso = T0 + T0 + Weffcj;
            PSsha = PDsha = T1;
            PSmer = PDmer = T2;

            ASiso = ADiso = T0 * Weffcj;
            ASsha = ADsha = DMCG * Weffcj;
            ASmer = ADmer = DMDG * Weffcj;

            switch (geo)
            {
                case 0:
                    Ps = nuEndS * PSiso + nuIntS * PSsha;
                    Pd = nuEndD * PDiso + nuIntD * PDsha;
                    As = nuEndS * ASiso + nuIntS * ASsha;
                    Ad = nuEndD * ADiso + nuIntD * ADsha;
                    break;
                case 1:
                    Ps = nuEndS * PSiso + nuIntS * PSsha;
                    Pd = (nuEndD + nuIntD) * PDsha;
                    As = nuEndS * ASiso + nuIntS * ASsha;
                    Ad = (nuEndD + nuIntD) * ADsha;
                    break;
                case 2:
                    Ps = (nuEndS + nuIntS) * PSsha;
                    Pd = nuEndD * PDiso + nuIntD * PDsha;
                    As = (nuEndS + nuIntS) * ASsha;
                    Ad = nuEndD * ADiso + nuIntD * ADsha;
                    break;
                case 3:
                    Ps = (nuEndS + nuIntS) * PSsha;
                    Pd = (nuEndD + nuIntD) * PDsha;
                    As = (nuEndS + nuIntS) * ASsha;
                    Ad = (nuEndD + nuIntD) * ADsha;
                    break;
                case 4:
                    Ps = nuEndS * PSiso + nuIntS * PSsha;
                    Pd = nuEndD * PDmer + nuIntD * PDsha;
                    As = nuEndS * ASiso + nuIntS * ASsha;
                    Ad = nuEndD * ADmer + nuIntD * ADsha;
                    break;
                case 5:
                    Ps = (nuEndS + nuIntS) * PSsha;
                    Pd = nuEndD * PDmer + nuIntD * PDsha;
                    As = (nuEndS + nuIntS) * ASsha;
                    Ad = nuEndD * ADmer + nuIntD * ADsha;
                    break;
                case 6:
                    Ps = nuEndS * PSmer + nuIntS * PSsha;
                    Pd = nuEndD * PDiso + nuIntD * PDsha;
                    As = nuEndS * ASmer + nuIntS * ASsha;
                    Ad = nuEndD * ADiso + nuIntD * ADsha;
                    break;
                case 7:
                    Ps = nuEndS * PSmer + nuIntS * PSsha;
                    Pd = (nuEndD + nuIntD) * PDsha;
                    As = nuEndS * ASmer + nuIntS * ASsha;
                    Ad = (nuEndD + nuIntD) * ADsha;
                    break;
                case 8:
                    Ps = nuEndS * PSmer + nuIntS * PSsha;
                    Pd = nuEndD * PDmer + nuIntD * PDsha;
                    As = nuEndS * ASmer + nuIntS * ASsha;
                    Ad = nuEndD * ADmer + nuIntD * ADsha;
                    break;
                case 9: /* geo = 9 and 10 happen only when nf = even */
                    Ps = PSiso + (nf - 1.0) * PSsha;
                    Pd = nf * PDsha;
                    As = ASiso + (nf - 1.0) * ASsha;
                    Ad = nf * ADsha;
                    break;
                case 10:
                    Ps = nf * PSsha;
                    Pd = PDiso + (nf - 1.0) * PDsha;
                    As = nf * ASsha;
                    Ad = ADiso + (nf - 1.0) * ADsha;
                    break;
                default:
                    Ps = 0;
                    Pd = 0;
                    As = 0;
                    Ad = 0;
                    SpiceSharpWarning.Warning(this, "Warning: Specified GEO = {0} not matched\n".FormatString(geo));
                    break;
            }
            return 0;
        }

        private int BSIM4RdsEndIso(double Weffcj, double Rsh, double DMCG, double DMCI, double DMDG, double nuEnd, int rgeo, int Type,
            out double Rend)
        {
            Rend = 0;
            if (Type == 1)
            {
                switch (rgeo)
                {
                    case 1:
                    case 2:
                    case 5:
                        if (nuEnd == 0.0)
                            Rend = 0.0;
                        else
                            Rend = Rsh * DMCG / (Weffcj * nuEnd);
                        break;
                    case 3:
                    case 4:
                    case 6:
                        if ((DMCG + DMCI) == 0.0)
                            SpiceSharpWarning.Warning(this, "(DMCG + DMCI) can not be equal to zero\n");
                        if ((nuEnd == 0.0) || ((DMCG + DMCI) == 0.0))
                            Rend = 0.0;
                        else
                            Rend = Rsh * Weffcj / (3.0 * nuEnd * (DMCG + DMCI));
                        break;
                    default:
                        SpiceSharpWarning.Warning(this, "Warning: Specified RGEO = %d not matched\n".FormatString(rgeo));
                        break;
                }
            }
            else
            {
                switch (rgeo)
                {
                    case 1:
                    case 3:
                    case 7:
                        if (nuEnd == 0.0)
                            Rend = 0.0;
                        else
                            Rend = Rsh * DMCG / (Weffcj * nuEnd);
                        break;
                    case 2:
                    case 4:
                    case 8:
                        if ((DMCG + DMCI) == 0.0)
                            SpiceSharpWarning.Warning(this, "(DMCG + DMCI) can not be equal to zero\n");
                        if ((nuEnd == 0.0) || ((DMCG + DMCI) == 0.0))
                            Rend = 0.0;
                        else
                            Rend = Rsh * Weffcj / (3.0 * nuEnd * (DMCG + DMCI));
                        break;
                    default:
                        SpiceSharpWarning.Warning(this, "Warning: Specified RGEO = %d not matched\n".FormatString(rgeo));
                        break;
                }
            }
            return 0;
        }

        private static int BSIM4NumFingerDiff(double nf, int minSD, 
            out double nuIntD, out double nuEndD, out double nuIntS, out double nuEndS)
        {
            int NF;
            NF = (int)nf;
            if ((NF % 2) != 0)
            {
                nuEndD = nuEndS = 1.0;
                nuIntD = nuIntS = 2.0 * Math.Max((nf - 1.0) / 2.0, 0.0);
            }
            else
            {
                if (minSD == 1) /* minimize # of source */
                {
                    nuEndD = 2.0;
                    nuIntD = 2.0 * Math.Max((nf / 2.0 - 1.0), 0.0);
                    nuEndS = 0.0;
                    nuIntS = nf;
                }
                else
                {
                    nuEndD = 0.0;
                    nuIntD = nf;
                    nuEndS = 2.0;
                    nuIntS = 2.0 * Math.Max((nf / 2.0 - 1.0), 0.0);
                }
            }
            return 0;
        }

        private int BSIM4RdseffGeo(double nf, int geo, int rgeo, int minSD, double Weffcj, double Rsh, double DMCG, double DMCI, double DMDG, int Type,
            out double Rtot)
        {
            double Rint = 0.0, Rend = 0.0;
            double nuIntD = 0.0, nuEndD = 0.0, nuIntS = 0.0, nuEndS = 0.0;

            if (geo < 9) /* since geo = 9 and 10 only happen when nf = even */
            {
                BSIM4NumFingerDiff(nf, minSD, out nuIntD, out nuEndD, out nuIntS, out nuEndS);

                /* Internal S/D resistance -- assume shared S or D and all wide contacts */
                if (Type == 1)
                {
                    if (nuIntS == 0.0)
                        Rint = 0.0;
                    else
                        Rint = Rsh * DMCG / (Weffcj * nuIntS);
                }
                else
                {
                    if (nuIntD == 0.0)
                        Rint = 0.0;
                    else
                        Rint = Rsh * DMCG / (Weffcj * nuIntD);
                }
            }

            /* End S/D resistance  -- geo dependent */
            switch (geo)
            {
                case 0:
                    if (Type == 1) BSIM4RdsEndIso(Weffcj, Rsh, DMCG, DMCI, DMDG,
                                      nuEndS, rgeo, 1, out Rend);
                    else BSIM4RdsEndIso(Weffcj, Rsh, DMCG, DMCI, DMDG,
                                           nuEndD, rgeo, 0, out Rend);
                    break;
                case 1:
                    if (Type == 1) BSIM4RdsEndIso(Weffcj, Rsh, DMCG, DMCI, DMDG,
                                                  nuEndS, rgeo, 1, out Rend);
                    else BSIM4RdsEndSha(Weffcj, Rsh, DMCG, DMCI, DMDG,
                              nuEndD, rgeo, 0, out Rend);
                    break;
                case 2:
                    if (Type == 1) BSIM4RdsEndSha(Weffcj, Rsh, DMCG, DMCI, DMDG,
                              nuEndS, rgeo, 1, out Rend);
                    else BSIM4RdsEndIso(Weffcj, Rsh, DMCG, DMCI, DMDG,
                              nuEndD, rgeo, 0, out Rend);
                    break;
                case 3:
                    if (Type == 1) BSIM4RdsEndSha(Weffcj, Rsh, DMCG, DMCI, DMDG,
                                                  nuEndS, rgeo, 1, out Rend);
                    else BSIM4RdsEndSha(Weffcj, Rsh, DMCG, DMCI, DMDG,
                                                  nuEndD, rgeo, 0, out Rend);
                    break;
                case 4:
                    if (Type == 1) BSIM4RdsEndIso(Weffcj, Rsh, DMCG, DMCI, DMDG,
                                                  nuEndS, rgeo, 1, out Rend);
                    else Rend = Rsh * DMDG / Weffcj;
                    break;
                case 5:
                    if (Type == 1) BSIM4RdsEndSha(Weffcj, Rsh, DMCG, DMCI, DMDG,
                                                  nuEndS, rgeo, 1, out Rend);
                    else Rend = Rsh * DMDG / (Weffcj * nuEndD);
                    break;
                case 6:
                    if (Type == 1) Rend = Rsh * DMDG / Weffcj;
                    else BSIM4RdsEndIso(Weffcj, Rsh, DMCG, DMCI, DMDG,
                                                  nuEndD, rgeo, 0, out Rend);
                    break;
                case 7:
                    if (Type == 1) Rend = Rsh * DMDG / (Weffcj * nuEndS);
                    else BSIM4RdsEndSha(Weffcj, Rsh, DMCG, DMCI, DMDG,
                                                  nuEndD, rgeo, 0, out Rend);
                    break;
                case 8:
                    Rend = Rsh * DMDG / Weffcj;
                    break;
                case 9: /* all wide contacts assumed for geo = 9 and 10 */
                    if (Type == 1)
                    {
                        Rend = 0.5 * Rsh * DMCG / Weffcj;
                        if (nf == 2.0)
                            Rint = 0.0;
                        else
                            Rint = Rsh * DMCG / (Weffcj * (nf - 2.0));
                    }
                    else
                    {
                        Rend = 0.0;
                        Rint = Rsh * DMCG / (Weffcj * nf);
                    }
                    break;
                case 10:
                    if (Type == 1)
                    {
                        Rend = 0.0;
                        Rint = Rsh * DMCG / (Weffcj * nf);
                    }
                    else
                    {
                        Rend = 0.5 * Rsh * DMCG / Weffcj; ;
                        if (nf == 2.0)
                            Rint = 0.0;
                        else
                            Rint = Rsh * DMCG / (Weffcj * (nf - 2.0));
                    }
                    break;
                default:
                    SpiceSharpWarning.Warning(this, "Warning: Specified GEO = %d not matched\n".FormatString(geo));
                break;
            }

            if (Rint <= 0.0)
                Rtot = Rend;
            else if (Rend <= 0.0)
                Rtot = Rint;
            else
                Rtot = Rint * Rend / (Rint + Rend);
            if (Rtot == 0.0)
                SpiceSharpWarning.Warning(this, "Warning: Zero resistance returned from RdseffGeo\n");
            return 0;
        }

        private int BSIM4RdsEndSha(double Weffcj, double Rsh, double DMCG, double DMCI, double DMDG, double nuEnd, int rgeo, int Type,
            out double Rend)
        {
            if (Type == 1)
            {
                switch (rgeo)
                {
                    case 1:
                    case 2:
                    case 5:
                        if (nuEnd == 0.0)
                            Rend = 0.0;
                        else
                            Rend = Rsh * DMCG / (Weffcj * nuEnd);
                        break;
                    case 3:
                    case 4:
                    case 6:
                        if (DMCG == 0.0)
                            SpiceSharpWarning.Warning(this, "DMCG can not be equal to zero\n");
                        if (nuEnd == 0.0)
                            Rend = 0.0;
                        else
                            Rend = Rsh * Weffcj / (6.0 * nuEnd * DMCG);
                        break;
                    default:
                        Rend = 0;
                        SpiceSharpWarning.Warning(this, "Warning: Specified RGEO = %d not matched\n".FormatString(rgeo));
                        break;
                }
            }
            else
            {
                switch (rgeo)
                {
                    case 1:
                    case 3:
                    case 7:
                        if (nuEnd == 0.0)
                            Rend = 0.0;
                        else
                            Rend = Rsh * DMCG / (Weffcj * nuEnd);
                        break;
                    case 2:
                    case 4:
                    case 8:
                        if (DMCG == 0.0)
                            SpiceSharpWarning.Warning(this, "DMCG can not be equal to zero\n");
                        if (nuEnd == 0.0)
                            Rend = 0.0;
                        else
                            Rend = Rsh * Weffcj / (6.0 * nuEnd * DMCG);
                        break;
                    default:
                        Rend = 0;
                        SpiceSharpWarning.Warning(this, "Warning: Specified RGEO = %d not matched\n".FormatString(rgeo));
                        break;
                }
            }
            return 0;
        }

    }
}
