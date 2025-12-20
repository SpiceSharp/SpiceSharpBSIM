using SpiceSharp;
using SpiceSharp.Attributes;
using SpiceSharp.Behaviors;
using SpiceSharp.Components;
using SpiceSharp.ParameterSets;
using SpiceSharp.Simulations;
using System;
using System.Collections.Generic;
using System.IO;

namespace SpiceSharpBSIM.Components.Semiconductors.BSIM.BSIM3v1Behaviors
{
    /// <summary>
    /// Temperature behavior for a <see cref="BSIM3v1"/>.
    /// </summary>
    [BehaviorFor(typeof(BSIM3v1)), AddBehaviorIfNo(typeof(ITemperatureBehavior))]
    public class TemperatureBehavior : Behavior, ITemperatureBehavior, IParameterized<BaseParameters>
    {
        private readonly ITemperatureSimulationState _temperature;

        public const double Kb = 1.3806226e-23;
        public const double KboQ = 8.617087e-5;  /* Kb / q  where q = 1.60219e-19 */
        public const double EPSOX = 3.453133e-11;
        public const double EPSSI = 1.03594e-10;
        public const double PI = 3.141592654;
        public const double MAX_EXP = 5.834617425e14;
        public const double MIN_EXP = 1.713908431e-15;
        public const double EXP_THRESHOLD = 34.0;
        public const double Charge_q = 1.60219e-19;

        /// <inheritdoc />
        public BaseParameters Parameters { get; }

        protected ModelTemperatureBehavior ModelTemperature { get; }
        protected ModelParameters ModelParameters { get; }

        protected SizeDependentProperties Param { get; private set; }

        protected double _drainConductance, _sourceConductance, _cgso, _cgdo;

        /// <summary>
        /// Gets the name of the model
        /// </summary>
        [ParameterName("model"), ParameterInfo("The name of the model.")]
        public string ModelName => ModelTemperature.Name;

        /// <summary>
        /// Creates a new <see cref="TemperatureBehavior"/>.
        /// </summary>
        /// <param name="context">The context.</param>
        public TemperatureBehavior(ComponentBindingContext context)
            : base(context)
        {
            _temperature = context.GetState<ITemperatureSimulationState>();
            Parameters = context.GetParameterSet<BaseParameters>();
            if (context.ModelBehaviors.TryGetValue<AggregateModelTemperatureBehavior>(out var aggregateBehavior))
            {
                ModelTemperature = aggregateBehavior.GetModel(Parameters.W, Parameters.L);
                ModelParameters = ModelTemperature.Parameters;
            }
            else
            {
                ModelParameters = context.ModelBehaviors.GetParameterSet<ModelParameters>();
                ModelTemperature = context.ModelBehaviors.GetValue<ModelTemperatureBehavior>();
            }
            Setup();
        }

        /// <summary>
        /// Set up the device.
        /// </summary>
        private void Setup()
        {
            if (!Parameters.DrainArea.Given)
            {
                if (ModelParameters.Hdif.Given)
                    Parameters.DrainArea = new GivenParameter<double>(Parameters.W * 2 * ModelParameters.Hdif, false);
                else
                    Parameters.DrainArea = new GivenParameter<double>(0.0, false);
            }
            if (!Parameters.DrainPerimeter.Given)
            {
                if (ModelParameters.Hdif.Given)
                    Parameters.DrainPerimeter = new GivenParameter<double>(2 * Parameters.W + 4 * ModelParameters.Hdif, false);
                else
                    Parameters.DrainPerimeter = new GivenParameter<double>(0.0, false);
            }
            if (!Parameters.SourceArea.Given)
            {
                if (ModelParameters.Hdif.Given)
                    Parameters.SourceArea = new GivenParameter<double>(Parameters.W * 2 * ModelParameters.Hdif, false);
                else
                    Parameters.SourceArea = new GivenParameter<double>(0.0, false);
            }
            if (!Parameters.SourcePerimeter.Given)
            {
                if (ModelParameters.Hdif.Given)
                    Parameters.SourcePerimeter = new GivenParameter<double>(2 * Parameters.W + 4 * ModelParameters.Hdif, false);
                else
                    Parameters.SourcePerimeter = new GivenParameter<double>(0.0, false);
            }
            if (!Parameters.NqsMod.Given)
                Parameters.NqsMod = new GivenParameter<int>(ModelParameters.NqsMod, false);
        }

        /// <inheritdoc />
        void ITemperatureBehavior.Temperature()
        {
            double tmp1, tmp2, ni, T0, T1, T2, T3, Ldrn, Wdrn;
            double TRatio, Inv_L, Inv_W, Inv_LW, Vtm0, Tnom;

            Tnom = ModelParameters.Tnom;
            TRatio = _temperature.Temperature / Tnom;
            Vtm0 = ModelTemperature.Vtm0;
            ni = ModelTemperature.Ni;

            var key = Tuple.Create(Parameters.W.Value, Parameters.L.Value);
            if (ModelTemperature.SizeDependentProperties.TryGetValue(key, out var param))
                Param = param;
            else
            {
                Param = new SizeDependentProperties();
                ModelTemperature.SizeDependentProperties.Add(key, Param);
                Ldrn = Parameters.L;
                Wdrn = Parameters.W;
                Param.Length = Ldrn;
                Param.Width = Wdrn;

                T0 = Math.Pow(Ldrn, ModelParameters.Lln);
                T1 = Math.Pow(Wdrn, ModelParameters.Lwn);
                tmp1 = ModelParameters.Ll / T0 + ModelParameters.Lw / T1
                     + ModelParameters.Lwl / (T0 * T1);
                Param.BSIM3v1dl = ModelParameters.Lint + tmp1;
                Param.BSIM3v1dlc = ModelParameters.Dlc + tmp1;

                T2 = Math.Pow(Ldrn, ModelParameters.Wln);
                T3 = Math.Pow(Wdrn, ModelParameters.Wwn);
                tmp2 = ModelParameters.Wl / T2 + ModelParameters.Ww / T3
                     + ModelParameters.Wwl / (T2 * T3);
                Param.BSIM3v1dw = ModelParameters.Wint + tmp2;
                Param.BSIM3v1dwc = ModelParameters.Dwc + tmp2;

                Param.BSIM3v1leff = Parameters.L - 2.0 * Param.BSIM3v1dl;
                if (Param.BSIM3v1leff <= 0.0)
                    throw new SpiceSharpException("BSIM3v1: mosfet {0}, model {1}: Effective channel length <= 0".FormatString(Name, ModelTemperature.Name));

                Param.BSIM3v1weff = Parameters.W - 2.0 * Param.BSIM3v1dw;
                if (Param.BSIM3v1weff <= 0.0)
                    throw new SpiceSharpException("BSIM3v1: mosfet {0}, model {1}: Effective channel width <= 0".FormatString(Name, ModelTemperature.Name));

                Param.BSIM3v1leffCV = Parameters.L - 2.0 * Param.BSIM3v1dlc;
                if (Param.BSIM3v1leffCV <= 0.0)
                    throw new SpiceSharpException("BSIM3v1: mosfet {0}, model {1}: Effective channel length for C-V <= 0".FormatString(Name, ModelTemperature.Name));

                Param.BSIM3v1weffCV = Parameters.W - 2.0 * Param.BSIM3v1dwc;
                if (Param.BSIM3v1weffCV <= 0.0)
                    throw new SpiceSharpException("BSIM3v1: mosfet {0}, model {1}: Effective channel width for C-V <= 0".FormatString(Name, ModelTemperature.Name));

                if (ModelParameters.BinUnit == 1)
                {
                    Inv_L = 1.0e-6 / Param.BSIM3v1leff;
                    Inv_W = 1.0e-6 / Param.BSIM3v1weff;
                    Inv_LW = 1.0e-12 / (Param.BSIM3v1leff
                           * Param.BSIM3v1weff);
                }
                else
                {
                    Inv_L = 1.0 / Param.BSIM3v1leff;
                    Inv_W = 1.0 / Param.BSIM3v1weff;
                    Inv_LW = 1.0 / (Param.BSIM3v1leff
                           * Param.BSIM3v1weff);
                }
                Param.BSIM3v1cdsc = ModelParameters.Cdsc
                                  + ModelParameters.Lcdsc * Inv_L
                                  + ModelParameters.Wcdsc * Inv_W
                                  + ModelParameters.Pcdsc * Inv_LW;
                Param.BSIM3v1cdscb = ModelParameters.Cdscb
                                   + ModelParameters.Lcdscb * Inv_L
                                   + ModelParameters.Wcdscb * Inv_W
                                   + ModelParameters.Pcdscb * Inv_LW;

                Param.BSIM3v1cdscd = ModelParameters.Cdscd
                               + ModelParameters.Lcdscd * Inv_L
                               + ModelParameters.Wcdscd * Inv_W
                               + ModelParameters.Pcdscd * Inv_LW;

                Param.BSIM3v1cit = ModelParameters.Cit
                                 + ModelParameters.Lcit * Inv_L
                                 + ModelParameters.Wcit * Inv_W
                                 + ModelParameters.Pcit * Inv_LW;
                Param.BSIM3v1nfactor = ModelParameters.Nfactor
                                     + ModelParameters.Lnfactor * Inv_L
                                     + ModelParameters.Wnfactor * Inv_W
                                     + ModelParameters.Pnfactor * Inv_LW;
                Param.BSIM3v1xj = ModelParameters.Xj
                                + ModelParameters.Lxj * Inv_L
                                + ModelParameters.Wxj * Inv_W
                                + ModelParameters.Pxj * Inv_LW;
                Param.BSIM3v1vsat = ModelParameters.Vsat
                                  + ModelParameters.Lvsat * Inv_L
                                  + ModelParameters.Wvsat * Inv_W
                                  + ModelParameters.Pvsat * Inv_LW;
                Param.BSIM3v1at = ModelParameters.At
                                + ModelParameters.Lat * Inv_L
                                + ModelParameters.Wat * Inv_W
                                + ModelParameters.Pat * Inv_LW;
                Param.BSIM3v1a0 = ModelParameters.A0
                                + ModelParameters.La0 * Inv_L
                                + ModelParameters.Wa0 * Inv_W
                                + ModelParameters.Pa0 * Inv_LW;

                Param.BSIM3v1ags = ModelParameters.Ags
                                + ModelParameters.Lags * Inv_L
                                + ModelParameters.Wags * Inv_W
                                + ModelParameters.Pags * Inv_LW;

                Param.BSIM3v1a1 = ModelParameters.A1
                                + ModelParameters.La1 * Inv_L
                                + ModelParameters.Wa1 * Inv_W
                                + ModelParameters.Pa1 * Inv_LW;
                Param.BSIM3v1a2 = ModelParameters.A2
                                + ModelParameters.La2 * Inv_L
                                + ModelParameters.Wa2 * Inv_W
                                + ModelParameters.Pa2 * Inv_LW;
                Param.BSIM3v1keta = ModelParameters.Keta
                                  + ModelParameters.Lketa * Inv_L
                                  + ModelParameters.Wketa * Inv_W
                                  + ModelParameters.Pketa * Inv_LW;
                Param.BSIM3v1nsub = ModelParameters.Nsub
                                  + ModelParameters.Lnsub * Inv_L
                                  + ModelParameters.Wnsub * Inv_W
                                  + ModelParameters.Pnsub * Inv_LW;
                Param.BSIM3v1npeak = ModelParameters.Npeak
                                   + ModelParameters.Lnpeak * Inv_L
                                   + ModelParameters.Wnpeak * Inv_W
                                   + ModelParameters.Pnpeak * Inv_LW;
                Param.BSIM3v1ngate = ModelParameters.Ngate
                                   + ModelParameters.Lngate * Inv_L
                                   + ModelParameters.Wngate * Inv_W
                                   + ModelParameters.Pngate * Inv_LW;
                Param.BSIM3v1gamma1 = ModelParameters.Gamma1
                                    + ModelParameters.Lgamma1 * Inv_L
                                    + ModelParameters.Wgamma1 * Inv_W
                                    + ModelParameters.Pgamma1 * Inv_LW;
                Param.BSIM3v1gamma2 = ModelParameters.Gamma2
                                    + ModelParameters.Lgamma2 * Inv_L
                                    + ModelParameters.Wgamma2 * Inv_W
                                    + ModelParameters.Pgamma2 * Inv_LW;
                Param.BSIM3v1vbx = ModelParameters.Vbx
                                 + ModelParameters.Lvbx * Inv_L
                                 + ModelParameters.Wvbx * Inv_W
                                 + ModelParameters.Pvbx * Inv_LW;
                Param.BSIM3v1vbm = ModelParameters.Vbm
                                 + ModelParameters.Lvbm * Inv_L
                                 + ModelParameters.Wvbm * Inv_W
                                 + ModelParameters.Pvbm * Inv_LW;
                Param.BSIM3v1xt = ModelParameters.Xt
                                 + ModelParameters.Lxt * Inv_L
                                 + ModelParameters.Wxt * Inv_W
                                 + ModelParameters.Pxt * Inv_LW;
                Param.BSIM3v1k1 = ModelParameters.K1
                                + ModelParameters.Lk1 * Inv_L
                                + ModelParameters.Wk1 * Inv_W
                                + ModelParameters.Pk1 * Inv_LW;
                Param.BSIM3v1kt1 = ModelParameters.Kt1
                                 + ModelParameters.Lkt1 * Inv_L
                                 + ModelParameters.Wkt1 * Inv_W
                                 + ModelParameters.Pkt1 * Inv_LW;
                Param.BSIM3v1kt1l = ModelParameters.Kt1l
                                  + ModelParameters.Lkt1l * Inv_L
                                  + ModelParameters.Wkt1l * Inv_W
                                  + ModelParameters.Pkt1l * Inv_LW;
                Param.BSIM3v1k2 = ModelParameters.K2
                                + ModelParameters.Lk2 * Inv_L
                                + ModelParameters.Wk2 * Inv_W
                                + ModelParameters.Pk2 * Inv_LW;
                Param.BSIM3v1kt2 = ModelParameters.Kt2
                                 + ModelParameters.Lkt2 * Inv_L
                                 + ModelParameters.Wkt2 * Inv_W
                                 + ModelParameters.Pkt2 * Inv_LW;
                Param.BSIM3v1k3 = ModelParameters.K3
                                + ModelParameters.Lk3 * Inv_L
                                + ModelParameters.Wk3 * Inv_W
                                + ModelParameters.Pk3 * Inv_LW;
                Param.BSIM3v1k3b = ModelParameters.K3b
                                 + ModelParameters.Lk3b * Inv_L
                                 + ModelParameters.Wk3b * Inv_W
                                 + ModelParameters.Pk3b * Inv_LW;
                Param.BSIM3v1w0 = ModelParameters.W0
                                + ModelParameters.Lw0 * Inv_L
                                + ModelParameters.Ww0 * Inv_W
                                + ModelParameters.Pw0 * Inv_LW;
                Param.BSIM3v1nlx = ModelParameters.Nlx
                                 + ModelParameters.Lnlx * Inv_L
                                 + ModelParameters.Wnlx * Inv_W
                                 + ModelParameters.Pnlx * Inv_LW;
                Param.BSIM3v1dvt0 = ModelParameters.Dvt0
                                  + ModelParameters.Ldvt0 * Inv_L
                                  + ModelParameters.Wdvt0 * Inv_W
                                  + ModelParameters.Pdvt0 * Inv_LW;
                Param.BSIM3v1dvt1 = ModelParameters.Dvt1
                                  + ModelParameters.Ldvt1 * Inv_L
                                  + ModelParameters.Wdvt1 * Inv_W
                                  + ModelParameters.Pdvt1 * Inv_LW;
                Param.BSIM3v1dvt2 = ModelParameters.Dvt2
                                  + ModelParameters.Ldvt2 * Inv_L
                                  + ModelParameters.Wdvt2 * Inv_W
                                  + ModelParameters.Pdvt2 * Inv_LW;
                Param.BSIM3v1dvt0w = ModelParameters.Dvt0w
                                  + ModelParameters.Ldvt0w * Inv_L
                                  + ModelParameters.Wdvt0w * Inv_W
                                  + ModelParameters.Pdvt0w * Inv_LW;
                Param.BSIM3v1dvt1w = ModelParameters.Dvt1w
                                  + ModelParameters.Ldvt1w * Inv_L
                                  + ModelParameters.Wdvt1w * Inv_W
                                  + ModelParameters.Pdvt1w * Inv_LW;
                Param.BSIM3v1dvt2w = ModelParameters.Dvt2w
                                  + ModelParameters.Ldvt2w * Inv_L
                                  + ModelParameters.Wdvt2w * Inv_W
                                  + ModelParameters.Pdvt2w * Inv_LW;
                Param.BSIM3v1drout = ModelParameters.Drout
                                   + ModelParameters.Ldrout * Inv_L
                                   + ModelParameters.Wdrout * Inv_W
                                   + ModelParameters.Pdrout * Inv_LW;
                Param.BSIM3v1dsub = ModelParameters.Dsub
                                  + ModelParameters.Ldsub * Inv_L
                                  + ModelParameters.Wdsub * Inv_W
                                  + ModelParameters.Pdsub * Inv_LW;
                Param.BSIM3v1vth0 = ModelParameters.Vth0
                                  + ModelParameters.Lvth0 * Inv_L
                                  + ModelParameters.Wvth0 * Inv_W
                                  + ModelParameters.Pvth0 * Inv_LW;
                Param.BSIM3v1ua = ModelParameters.Ua
                                + ModelParameters.Lua * Inv_L
                                + ModelParameters.Wua * Inv_W
                                + ModelParameters.Pua * Inv_LW;
                Param.BSIM3v1ua1 = ModelParameters.Ua1
                                 + ModelParameters.Lua1 * Inv_L
                                 + ModelParameters.Wua1 * Inv_W
                                 + ModelParameters.Pua1 * Inv_LW;
                Param.BSIM3v1ub = ModelParameters.Ub
                                + ModelParameters.Lub * Inv_L
                                + ModelParameters.Wub * Inv_W
                                + ModelParameters.Pub * Inv_LW;
                Param.BSIM3v1ub1 = ModelParameters.Ub1
                                 + ModelParameters.Lub1 * Inv_L
                                 + ModelParameters.Wub1 * Inv_W
                                 + ModelParameters.Pub1 * Inv_LW;
                Param.BSIM3v1uc = ModelParameters.Uc
                                + ModelParameters.Luc * Inv_L
                                + ModelParameters.Wuc * Inv_W
                                + ModelParameters.Puc * Inv_LW;
                Param.BSIM3v1uc1 = ModelParameters.Uc1
                                 + ModelParameters.Luc1 * Inv_L
                                 + ModelParameters.Wuc1 * Inv_W
                                 + ModelParameters.Puc1 * Inv_LW;
                Param.BSIM3v1u0 = ModelParameters.U0
                                + ModelParameters.Lu0 * Inv_L
                                + ModelParameters.Wu0 * Inv_W
                                + ModelParameters.Pu0 * Inv_LW;
                Param.BSIM3v1ute = ModelParameters.Ute
                                 + ModelParameters.Lute * Inv_L
                                 + ModelParameters.Wute * Inv_W
                                 + ModelParameters.Pute * Inv_LW;
                Param.BSIM3v1voff = ModelParameters.Voff
                                  + ModelParameters.Lvoff * Inv_L
                                  + ModelParameters.Wvoff * Inv_W
                                  + ModelParameters.Pvoff * Inv_LW;
                Param.BSIM3v1delta = ModelParameters.Delta
                                   + ModelParameters.Ldelta * Inv_L
                                   + ModelParameters.Wdelta * Inv_W
                                   + ModelParameters.Pdelta * Inv_LW;
                Param.BSIM3v1rdsw = ModelParameters.Rdsw
                                  + ModelParameters.Lrdsw * Inv_L
                                  + ModelParameters.Wrdsw * Inv_W
                                  + ModelParameters.Prdsw * Inv_LW;
                Param.BSIM3v1prwg = ModelParameters.Prwg
                                  + ModelParameters.Lprwg * Inv_L
                                  + ModelParameters.Wprwg * Inv_W
                                  + ModelParameters.Pprwg * Inv_LW;
                Param.BSIM3v1prwb = ModelParameters.Prwb
                                  + ModelParameters.Lprwb * Inv_L
                                  + ModelParameters.Wprwb * Inv_W
                                  + ModelParameters.Pprwb * Inv_LW;
                Param.BSIM3v1prt = ModelParameters.Prt
                                  + ModelParameters.Lprt * Inv_L
                                  + ModelParameters.Wprt * Inv_W
                                  + ModelParameters.Pprt * Inv_LW;
                Param.BSIM3v1eta0 = ModelParameters.Eta0
                                  + ModelParameters.Leta0 * Inv_L
                                  + ModelParameters.Weta0 * Inv_W
                                  + ModelParameters.Peta0 * Inv_LW;
                Param.BSIM3v1etab = ModelParameters.Etab
                                  + ModelParameters.Letab * Inv_L
                                  + ModelParameters.Wetab * Inv_W
                                  + ModelParameters.Petab * Inv_LW;
                Param.BSIM3v1pclm = ModelParameters.Pclm
                                  + ModelParameters.Lpclm * Inv_L
                                  + ModelParameters.Wpclm * Inv_W
                                  + ModelParameters.Ppclm * Inv_LW;
                Param.BSIM3v1pdibl1 = ModelParameters.Pdibl1
                                    + ModelParameters.Lpdibl1 * Inv_L
                                    + ModelParameters.Wpdibl1 * Inv_W
                                    + ModelParameters.Ppdibl1 * Inv_LW;
                Param.BSIM3v1pdibl2 = ModelParameters.Pdibl2
                                    + ModelParameters.Lpdibl2 * Inv_L
                                    + ModelParameters.Wpdibl2 * Inv_W
                                    + ModelParameters.Ppdibl2 * Inv_LW;
                Param.BSIM3v1pdiblb = ModelParameters.Pdiblb
                                    + ModelParameters.Lpdiblb * Inv_L
                                    + ModelParameters.Wpdiblb * Inv_W
                                    + ModelParameters.Ppdiblb * Inv_LW;
                Param.BSIM3v1pscbe1 = ModelParameters.Pscbe1
                                    + ModelParameters.Lpscbe1 * Inv_L
                                    + ModelParameters.Wpscbe1 * Inv_W
                                    + ModelParameters.Ppscbe1 * Inv_LW;
                Param.BSIM3v1pscbe2 = ModelParameters.Pscbe2
                                    + ModelParameters.Lpscbe2 * Inv_L
                                    + ModelParameters.Wpscbe2 * Inv_W
                                    + ModelParameters.Ppscbe2 * Inv_LW;
                Param.BSIM3v1pvag = ModelParameters.Pvag
                                  + ModelParameters.Lpvag * Inv_L
                                  + ModelParameters.Wpvag * Inv_W
                                  + ModelParameters.Ppvag * Inv_LW;
                Param.BSIM3v1wr = ModelParameters.Wr
                                + ModelParameters.Lwr * Inv_L
                                + ModelParameters.Wwr * Inv_W
                                + ModelParameters.Pwr * Inv_LW;
                Param.BSIM3v1dwg = ModelParameters.Dwg
                                 + ModelParameters.Ldwg * Inv_L
                                 + ModelParameters.Wdwg * Inv_W
                                 + ModelParameters.Pdwg * Inv_LW;
                Param.BSIM3v1dwb = ModelParameters.Dwb
                                 + ModelParameters.Ldwb * Inv_L
                                 + ModelParameters.Wdwb * Inv_W
                                 + ModelParameters.Pdwb * Inv_LW;
                Param.BSIM3v1b0 = ModelParameters.B0
                                + ModelParameters.Lb0 * Inv_L
                                + ModelParameters.Wb0 * Inv_W
                                + ModelParameters.Pb0 * Inv_LW;
                Param.BSIM3v1b1 = ModelParameters.B1
                                + ModelParameters.Lb1 * Inv_L
                                + ModelParameters.Wb1 * Inv_W
                                + ModelParameters.Pb1 * Inv_LW;
                Param.BSIM3v1alpha0 = ModelParameters.Alpha0
                                    + ModelParameters.Lalpha0 * Inv_L
                                    + ModelParameters.Walpha0 * Inv_W
                                    + ModelParameters.Palpha0 * Inv_LW;
                Param.BSIM3v1beta0 = ModelParameters.Beta0
                                   + ModelParameters.Lbeta0 * Inv_L
                                   + ModelParameters.Wbeta0 * Inv_W
                                   + ModelParameters.Pbeta0 * Inv_LW;
                /* CV model */
                Param.BSIM3v1elm = ModelParameters.Elm
                                + ModelParameters.Lelm * Inv_L
                                + ModelParameters.Welm * Inv_W
                                + ModelParameters.Pelm * Inv_LW;
                Param.BSIM3v1cgsl = ModelParameters.Cgsl
                                  + ModelParameters.Lcgsl * Inv_L
                                  + ModelParameters.Wcgsl * Inv_W
                                  + ModelParameters.Pcgsl * Inv_LW;
                Param.BSIM3v1cgdl = ModelParameters.Cgdl
                                  + ModelParameters.Lcgdl * Inv_L
                                  + ModelParameters.Wcgdl * Inv_W
                                  + ModelParameters.Pcgdl * Inv_LW;
                Param.BSIM3v1ckappa = ModelParameters.Ckappa
                                    + ModelParameters.Lckappa * Inv_L
                                    + ModelParameters.Wckappa * Inv_W
                                    + ModelParameters.Pckappa * Inv_LW;
                Param.BSIM3v1cf = ModelParameters.Cf
                                + ModelParameters.Lcf * Inv_L
                                + ModelParameters.Wcf * Inv_W
                                + ModelParameters.Pcf * Inv_LW;
                Param.BSIM3v1clc = ModelParameters.Clc
                                 + ModelParameters.Lclc * Inv_L
                                 + ModelParameters.Wclc * Inv_W
                                 + ModelParameters.Pclc * Inv_LW;
                Param.BSIM3v1cle = ModelParameters.Cle
                                 + ModelParameters.Lcle * Inv_L
                                 + ModelParameters.Wcle * Inv_W
                                 + ModelParameters.Pcle * Inv_LW;
                Param.BSIM3v1vfbcv = ModelParameters.Vfbcv
                                + ModelParameters.Lvfbcv * Inv_L
                                + ModelParameters.Wvfbcv * Inv_W
                                + ModelParameters.Pvfbcv * Inv_LW;
                Param.BSIM3v1abulkCVfactor = 1.0 + Math.Pow((Param.BSIM3v1clc
                                           / Param.BSIM3v1leff),
                                           Param.BSIM3v1cle);

                T0 = (TRatio - 1.0);
                Param.BSIM3v1ua = Param.BSIM3v1ua + Param.BSIM3v1ua1 * T0;
                Param.BSIM3v1ub = Param.BSIM3v1ub + Param.BSIM3v1ub1 * T0;
                Param.BSIM3v1uc = Param.BSIM3v1uc + Param.BSIM3v1uc1 * T0;
                if (Param.BSIM3v1u0 > 1.0)
                    Param.BSIM3v1u0 = Param.BSIM3v1u0 / 1.0e4;

                Param.BSIM3v1u0temp = Param.BSIM3v1u0
                                    * Math.Pow(TRatio, Param.BSIM3v1ute);
                Param.BSIM3v1vsattemp = Param.BSIM3v1vsat - Param.BSIM3v1at
                                      * T0;
                Param.BSIM3v1rds0 = (Param.BSIM3v1rdsw + Param.BSIM3v1prt * T0)
                                  / Math.Pow(Param.BSIM3v1weff * 1E6, Param.BSIM3v1wr);

                if (BSIM3v1checkModel())
                    throw new SpiceSharpException("Fatal error(s) detected during BSIM3V3.1 parameter checking for {0} in model {1}".FormatString(Name, ModelTemperature.Name));

                Param.BSIM3v1cgdo = (ModelParameters.Cgdo + Param.BSIM3v1cf)
                                  * Param.BSIM3v1weffCV;
                Param.BSIM3v1cgso = (ModelParameters.Cgso + Param.BSIM3v1cf)
                                  * Param.BSIM3v1weffCV;
                Param.BSIM3v1cgbo = ModelParameters.Cgbo * Param.BSIM3v1leffCV;

                if (!ModelParameters.Npeak.Given && ModelParameters.Gamma1.Given)
                {
                    T0 = Param.BSIM3v1gamma1 * ModelTemperature.Cox;
                    Param.BSIM3v1npeak = 3.021E22 * T0 * T0;
                }

                Param.BSIM3v1phi = 2.0 * Vtm0
                                 * Math.Log(Param.BSIM3v1npeak / ni);

                Param.BSIM3v1sqrtPhi = Math.Sqrt(Param.BSIM3v1phi);
                Param.BSIM3v1phis3 = Param.BSIM3v1sqrtPhi * Param.BSIM3v1phi;

                Param.BSIM3v1Xdep0 = Math.Sqrt(2.0 * EPSSI / (Charge_q
                                   * Param.BSIM3v1npeak * 1.0e6))
                                   * Param.BSIM3v1sqrtPhi;
                Param.BSIM3v1sqrtXdep0 = Math.Sqrt(Param.BSIM3v1Xdep0);
                Param.BSIM3v1litl = Math.Sqrt(3.0 * Param.BSIM3v1xj
                                  * ModelParameters.Tox);
                Param.BSIM3v1vbi = Vtm0 * Math.Log(1.0e20
                                 * Param.BSIM3v1npeak / (ni * ni));
                Param.BSIM3v1cdep0 = Math.Sqrt(Charge_q * EPSSI
                                   * Param.BSIM3v1npeak * 1.0e6 / 2.0
                                   / Param.BSIM3v1phi);

                if (ModelParameters.K1.Given || ModelParameters.K2.Given)
                {
                    if (!ModelParameters.K1.Given)
                    {
                        SpiceSharpWarning.Warning(this, "Warning: k1 should be specified with k2.");
                        Param.BSIM3v1k1 = 0.53;
                    }
                    if (!ModelParameters.K2.Given)
                    {
                        SpiceSharpWarning.Warning(this, "Warning: k2 should be specified with k1.");
                        Param.BSIM3v1k2 = -0.0186;
                    }
                    if (ModelParameters.Nsub.Given)
                        SpiceSharpWarning.Warning(this, "Warning: nsub is ignored because k1 or k2 is given.");
                    if (ModelParameters.Xt.Given)
                        SpiceSharpWarning.Warning(this, "Warning: xt is ignored because k1 or k2 is given.");
                    if (ModelParameters.Vbx.Given)
                        SpiceSharpWarning.Warning(this, "Warning: vbx is ignored because k1 or k2 is given.");
                    if (ModelParameters.Vbm.Given)
                        SpiceSharpWarning.Warning(this, "Warning: vbm is ignored because k1 or k2 is given.");
                    if (ModelParameters.Gamma1.Given)
                        SpiceSharpWarning.Warning(this, "Warning: gamma1 is ignored because k1 or k2 is given.");
                    if (ModelParameters.Gamma2.Given)
                        SpiceSharpWarning.Warning(this, "Warning: gamma2 is ignored because k1 or k2 is given.");
                }
                else
                {
                    if (!ModelParameters.Vbx.Given)
                        Param.BSIM3v1vbx = Param.BSIM3v1phi - 7.7348e-4
                                         * Param.BSIM3v1npeak
                                         * Param.BSIM3v1xt * Param.BSIM3v1xt;
                    if (Param.BSIM3v1vbx > 0.0)
                        Param.BSIM3v1vbx = -Param.BSIM3v1vbx;
                    if (Param.BSIM3v1vbm > 0.0)
                        Param.BSIM3v1vbm = -Param.BSIM3v1vbm;

                    if (!ModelParameters.Gamma1.Given)
                        Param.BSIM3v1gamma1 = 5.753e-12
                                            * Math.Sqrt(Param.BSIM3v1npeak)
                                            / ModelTemperature.Cox;
                    if (!ModelParameters.Gamma2.Given)
                        Param.BSIM3v1gamma2 = 5.753e-12
                                            * Math.Sqrt(Param.BSIM3v1nsub)
                                            / ModelTemperature.Cox;

                    T0 = Param.BSIM3v1gamma1 - Param.BSIM3v1gamma2;
                    T1 = Math.Sqrt(Param.BSIM3v1phi - Param.BSIM3v1vbx)
                       - Param.BSIM3v1sqrtPhi;
                    T2 = Math.Sqrt(Param.BSIM3v1phi * (Param.BSIM3v1phi
                       - Param.BSIM3v1vbm)) - Param.BSIM3v1phi;
                    Param.BSIM3v1k2 = T0 * T1 / (2.0 * T2 + Param.BSIM3v1vbm);
                    Param.BSIM3v1k1 = Param.BSIM3v1gamma2 - 2.0
                                    * Param.BSIM3v1k2 * Math.Sqrt(Param.BSIM3v1phi
                                    - Param.BSIM3v1vbm);
                }

                if (Param.BSIM3v1k2 < 0.0)
                {
                    T0 = 0.5 * Param.BSIM3v1k1 / Param.BSIM3v1k2;
                    Param.BSIM3v1vbsc = 0.9 * (Param.BSIM3v1phi - T0 * T0);
                    if (Param.BSIM3v1vbsc > -3.0)
                        Param.BSIM3v1vbsc = -3.0;
                    else if (Param.BSIM3v1vbsc < -30.0)
                        Param.BSIM3v1vbsc = -30.0;
                }
                else
                {
                    Param.BSIM3v1vbsc = -30.0;
                }
                if (Param.BSIM3v1vbsc > Param.BSIM3v1vbm)
                    Param.BSIM3v1vbsc = Param.BSIM3v1vbm;

                if (ModelParameters.Vth0.Given)
                {
                    Param.BSIM3v1vfb = ModelParameters.Type * Param.BSIM3v1vth0
                                     - Param.BSIM3v1phi - Param.BSIM3v1k1
                                     * Param.BSIM3v1sqrtPhi;
                }
                else
                {
                    Param.BSIM3v1vfb = -1.0;
                    Param.BSIM3v1vth0 = ModelParameters.Type * (Param.BSIM3v1vfb
                                      + Param.BSIM3v1phi + Param.BSIM3v1k1
                                      * Param.BSIM3v1sqrtPhi);
                }
                T1 = Math.Sqrt(EPSSI / EPSOX * ModelParameters.Tox
                   * Param.BSIM3v1Xdep0);
                T0 = Math.Exp(-0.5 * Param.BSIM3v1dsub * Param.BSIM3v1leff / T1);
                Param.BSIM3v1theta0vb0 = (T0 + 2.0 * T0 * T0);

                T0 = Math.Exp(-0.5 * Param.BSIM3v1drout * Param.BSIM3v1leff / T1);
                T2 = (T0 + 2.0 * T0 * T0);
                Param.BSIM3v1thetaRout = Param.BSIM3v1pdibl1 * T2
                                       + Param.BSIM3v1pdibl2;
            }

            /* process source/drain series resistance */
            this._drainConductance = ModelParameters.SheetResistance
                                            * Parameters.DrainSquares;
            if (this._drainConductance > 0.0)
                this._drainConductance = 1.0
                                            / this._drainConductance;
            else
                this._drainConductance = 0.0;

            this._sourceConductance = ModelParameters.SheetResistance
                                         * Parameters.SourceSquares;
            if (this._sourceConductance > 0.0)
                this._sourceConductance = 1.0
                                             / this._sourceConductance;
            else
                this._sourceConductance = 0.0;
            this._cgso = Param.BSIM3v1cgso;
            this._cgdo = Param.BSIM3v1cgdo;
        }

        private bool BSIM3v1checkModel()
        {
            int Fatal_Flag = 0;
            List<string> words = new List<string>();
            words.Add("BSIM3V3.1 Parameter Check");
            words.Add("Model = {0}".FormatString(ModelTemperature.Name));
            words.Add("W = {0}, L = {1}".FormatString(Parameters.W, Parameters.L));

            if (Param.BSIM3v1nlx < -Param.BSIM3v1leff)
            {
                SpiceSharpWarning.Warning(this, "Fatal: Nlx = {0:g} is less than -Leff.".FormatString(Param.BSIM3v1nlx));
                words.Add("Fatal: Nlx = {0:g} is less than -Leff.".FormatString(Param.BSIM3v1nlx));
                Fatal_Flag = 1;
            }

            if (ModelParameters.Tox <= 0.0)
            {
                SpiceSharpWarning.Warning(this, "Fatal: Tox = {0:g} is not positive.".FormatString(ModelParameters.Tox));
                words.Add("Fatal: Tox = {0:g} is not positive.".FormatString(ModelParameters.Tox));
                Fatal_Flag = 1;
            }

            if (Param.BSIM3v1npeak <= 0.0)
            {
                SpiceSharpWarning.Warning(this, "Fatal: Nch = {0:g} is not positive.".FormatString(Param.BSIM3v1npeak));
                words.Add("Fatal: Nch = {0:g} is not positive.".FormatString(Param.BSIM3v1npeak));
                Fatal_Flag = 1;
            }
            if (Param.BSIM3v1nsub <= 0.0)
            {
                SpiceSharpWarning.Warning(this, "Fatal: Nsub = {0:g} is not positive.".FormatString(Param.BSIM3v1nsub));
                words.Add("Fatal: Nsub = {0:g} is not positive.".FormatString(Param.BSIM3v1nsub));
                Fatal_Flag = 1;
            }
            if (Param.BSIM3v1ngate < 0.0)
            {
                SpiceSharpWarning.Warning(this, "Fatal: Ngate = {0:g} is not positive.".FormatString(Param.BSIM3v1ngate));
                words.Add("Fatal: Ngate = {0:g} Ngate is not positive.".FormatString(Param.BSIM3v1ngate));
                Fatal_Flag = 1;
            }
            if (Param.BSIM3v1ngate > 1.0e25)
            {
                SpiceSharpWarning.Warning(this, "Fatal: Ngate = {0:g} is too high.".FormatString(Param.BSIM3v1ngate));
                words.Add("Fatal: Ngate = {0:g} Ngate is too high".FormatString(Param.BSIM3v1ngate));
                Fatal_Flag = 1;
            }
            if (Param.BSIM3v1xj <= 0.0)
            {
                SpiceSharpWarning.Warning(this, "Fatal: Xj = {0:g} is not positive.".FormatString(Param.BSIM3v1xj));
                words.Add("Fatal: Xj = {0:g} is not positive.".FormatString(Param.BSIM3v1xj));
                Fatal_Flag = 1;
            }

            if (Param.BSIM3v1dvt1 < 0.0)
            {
                SpiceSharpWarning.Warning(this, "Fatal: Dvt1 = {0:g} is negative.".FormatString(Param.BSIM3v1dvt1));
                words.Add("Fatal: Dvt1 = {0:g} is negative.".FormatString(Param.BSIM3v1dvt1));
                Fatal_Flag = 1;
            }

            if (Param.BSIM3v1dvt1w < 0.0)
            {
                SpiceSharpWarning.Warning(this, "Fatal: Dvt1w = {0:g} is negative.".FormatString(Param.BSIM3v1dvt1w));
                words.Add("Fatal: Dvt1w = {0:g} is negative.".FormatString(Param.BSIM3v1dvt1w));
                Fatal_Flag = 1;
            }

            if (Param.BSIM3v1w0 == -Param.BSIM3v1weff)
            {
                SpiceSharpWarning.Warning(this, "Fatal: (W0 + Weff) = 0 cauing divided-by-zero.");
                words.Add("Fatal: (W0 + Weff) = 0 cauing divided-by-zero.");
                Fatal_Flag = 1;
            }

            if (Param.BSIM3v1dsub < 0.0)
            {
                SpiceSharpWarning.Warning(this, "Fatal: Dsub = {0:g} is negative.".FormatString(Param.BSIM3v1dsub));
                words.Add("Fatal: Dsub = {0:g} is negative.".FormatString(Param.BSIM3v1dsub));
                Fatal_Flag = 1;
            }
            if (Param.BSIM3v1b1 == -Param.BSIM3v1weff)
            {
                SpiceSharpWarning.Warning(this, "Fatal: (B1 + Weff) = 0 causing divided-by-zero.");
                words.Add("Fatal: (B1 + Weff) = 0 causing divided-by-zero.");
                Fatal_Flag = 1;
            }
            if (Param.BSIM3v1u0temp <= 0.0)
            {
                SpiceSharpWarning.Warning(this, "Fatal: u0 at current temperature = {0:g} is not positive.".FormatString(Param.BSIM3v1u0temp));
                words.Add("Fatal: u0 at current temperature = {0:g} is not positive.".FormatString(Param.BSIM3v1u0temp));
                Fatal_Flag = 1;
            }

            /* Check delta parameter */
            if (Param.BSIM3v1delta < 0.0)
            {
                SpiceSharpWarning.Warning(this, "Fatal: Delta = {0:g} is less than zero.".FormatString(Param.BSIM3v1delta));
                words.Add("Fatal: Delta = {0:g} is less than zero.".FormatString(Param.BSIM3v1delta));
                Fatal_Flag = 1;
            }

            if (Param.BSIM3v1vsattemp <= 0.0)
            {
                SpiceSharpWarning.Warning(this, "Fatal: Vsat at current temperature = {0:g} is not positive.".FormatString(Param.BSIM3v1vsattemp));
                words.Add("Fatal: Vsat at current temperature = {0:g} is not positive.".FormatString(Param.BSIM3v1vsattemp));
                Fatal_Flag = 1;
            }
            /* Check Rout parameters */
            if (Param.BSIM3v1pclm <= 0.0)
            {
                SpiceSharpWarning.Warning(this, "Fatal: Pclm = {0:g} is not positive.".FormatString(Param.BSIM3v1pclm));
                words.Add("Fatal: Pclm = {0:g} is not positive.".FormatString(Param.BSIM3v1pclm));
                Fatal_Flag = 1;
            }

            if (Param.BSIM3v1drout < 0.0)
            {
                SpiceSharpWarning.Warning(this, "Fatal: Drout = {0:g} is negative.".FormatString(Param.BSIM3v1drout));
                words.Add("Fatal: Drout = {0:g} is negative.".FormatString(Param.BSIM3v1drout));
                Fatal_Flag = 1;
            }
            if (ModelParameters.UnitLengthSidewallJctCap > 0.0 ||
                  ModelParameters.UnitLengthGateSidewallJctCap > 0.0)
            {
                if (Parameters.DrainPerimeter < Param.BSIM3v1weff)
                {
                    SpiceSharpWarning.Warning(this, "Warning: Pd = {0:g} is less than W.".FormatString(Parameters.DrainPerimeter));
                    words.Add("Warning: Pd = {0:g} is less than W.".FormatString(Parameters.DrainPerimeter));
                    Parameters.DrainPerimeter = Param.BSIM3v1weff;
                }
                if (Parameters.SourcePerimeter < Param.BSIM3v1weff)
                {
                    SpiceSharpWarning.Warning(this, "Warning: Ps = {0:g} is less than W.".FormatString(Parameters.SourcePerimeter));
                    words.Add("Warning: Ps = {0:g} is less than W.".FormatString(Parameters.SourcePerimeter));
                    Parameters.SourcePerimeter = Param.BSIM3v1weff;
                }
            }
            /* Check capacitance parameters */
            if (Param.BSIM3v1clc < 0.0)
            {
                SpiceSharpWarning.Warning(this, "Fatal: Clc = {0:g} is negative.".FormatString(Param.BSIM3v1clc));
                words.Add("Fatal: Clc = {0:g} is negative.".FormatString(Param.BSIM3v1clc));
                Fatal_Flag = 1;
            }
            if (ModelParameters.ParamChk == 1)
            {
                /* Check L and W parameters */
                if (Param.BSIM3v1leff <= 5.0e-8)
                {
                    SpiceSharpWarning.Warning(this, "Warning: Leff = {0:g} may be too small.".FormatString(Param.BSIM3v1leff));
                    words.Add("Warning: Leff = {0:g} may be too small.".FormatString(Param.BSIM3v1leff));
                }

                if (Param.BSIM3v1leffCV <= 5.0e-8)
                {
                    SpiceSharpWarning.Warning(this, "Warning: Leff for CV = {0:g} may be too small.".FormatString(Param.BSIM3v1leffCV));
                    words.Add("Warning: Leff for CV = {0:g} may be too small.".FormatString(Param.BSIM3v1leffCV));
                }

                if (Param.BSIM3v1weff <= 1.0e-7)
                {
                    SpiceSharpWarning.Warning(this, "Warning: Weff = {0:g} may be too small.".FormatString(Param.BSIM3v1weff));
                    words.Add("Warning: Weff = {0:g} may be too small.".FormatString(Param.BSIM3v1weff));
                }

                if (Param.BSIM3v1weffCV <= 1.0e-7)
                {
                    SpiceSharpWarning.Warning(this, "Warning: Weff for CV = {0:g} may be too small.".FormatString(Param.BSIM3v1weffCV));
                    words.Add("Warning: Weff for CV = {0:g} may be too small.".FormatString(Param.BSIM3v1weffCV));
                }

                /* Check threshold voltage parameters */
                if (Param.BSIM3v1nlx < 0.0)
                {
                    SpiceSharpWarning.Warning(this, "Warning: Nlx = {0:g} is negative.".FormatString(Param.BSIM3v1nlx));
                    words.Add("Warning: Nlx = {0:g} is negative.".FormatString(Param.BSIM3v1nlx));
                }
                if (ModelParameters.Tox < 1.0e-9)
                {
                    SpiceSharpWarning.Warning(this, "Warning: Tox = {0:g} is less than 10A.".FormatString(ModelParameters.Tox));
                    words.Add("Warning: Tox = {0:g} is less than 10A.".FormatString(ModelParameters.Tox));
                }

                if (Param.BSIM3v1npeak <= 1.0e15)
                {
                    SpiceSharpWarning.Warning(this, "Warning: Nch = {0:g} may be too small.".FormatString(Param.BSIM3v1npeak));
                    words.Add("Warning: Nch = {0:g} may be too small.".FormatString(Param.BSIM3v1npeak));
                }
                else if (Param.BSIM3v1npeak >= 1.0e21)
                {
                    SpiceSharpWarning.Warning(this, "Warning: Nch = {0:g} may be too large.".FormatString(Param.BSIM3v1npeak));
                    words.Add("Warning: Nch = {0:g} may be too large.".FormatString(Param.BSIM3v1npeak));
                }

                if (Param.BSIM3v1nsub <= 1.0e14)
                {
                    SpiceSharpWarning.Warning(this, "Warning: Nsub = {0:g} may be too small.".FormatString(Param.BSIM3v1nsub));
                    words.Add("Warning: Nsub = {0:g} may be too small.".FormatString(Param.BSIM3v1nsub));
                }
                else if (Param.BSIM3v1nsub >= 1.0e21)
                {
                    SpiceSharpWarning.Warning(this, "Warning: Nsub = {0:g} may be too large.".FormatString(Param.BSIM3v1nsub));
                    words.Add("Warning: Nsub = {0:g} may be too large.".FormatString(Param.BSIM3v1nsub));
                }

                if ((Param.BSIM3v1ngate > 0.0) &&
                    (Param.BSIM3v1ngate <= 1.0e18))
                {
                    SpiceSharpWarning.Warning(this, "Warning: Ngate = {0:g} is less than 1.0E18cm^-3.".FormatString(Param.BSIM3v1ngate));
                    words.Add("Warning: Ngate = {0:g} is less than 1.0E18cm^-3.".FormatString(Param.BSIM3v1ngate));
                }

                if (Param.BSIM3v1dvt0 < 0.0)
                {
                    SpiceSharpWarning.Warning(this, "Warning: Dvt0 = {0:g} is negative.".FormatString(Param.BSIM3v1dvt0));
                    words.Add("Warning: Dvt0 = {0:g} is negative.".FormatString(Param.BSIM3v1dvt0));
                }

                if (Math.Abs(1.0e-6 / (Param.BSIM3v1w0 + Param.BSIM3v1weff)) > 10.0)
                {
                    SpiceSharpWarning.Warning(this, "Warning: (W0 + Weff) may be too small.");
                    words.Add("Warning: (W0 + Weff) may be too small.");
                }

                /* Check subthreshold parameters */
                if (Param.BSIM3v1nfactor < 0.0)
                {
                    SpiceSharpWarning.Warning(this, "Warning: Nfactor = {0:g} is negative.".FormatString(Param.BSIM3v1nfactor));
                    words.Add("Warning: Nfactor = {0:g} is negative.".FormatString(Param.BSIM3v1nfactor));
                }
                if (Param.BSIM3v1cdsc < 0.0)
                {
                    SpiceSharpWarning.Warning(this, "Warning: Cdsc = {0:g} is negative.".FormatString(Param.BSIM3v1cdsc));
                    words.Add("Warning: Cdsc = {0:g} is negative.".FormatString(Param.BSIM3v1cdsc));
                }
                if (Param.BSIM3v1cdscd < 0.0)
                {
                    SpiceSharpWarning.Warning(this, "Warning: Cdscd = {0:g} is negative.".FormatString(Param.BSIM3v1cdscd));
                    words.Add("Warning: Cdscd = {0:g} is negative.".FormatString(Param.BSIM3v1cdscd));
                }
                /* Check DIBL parameters */
                if (Param.BSIM3v1eta0 < 0.0)
                {
                    SpiceSharpWarning.Warning(this, "Warning: Eta0 = {0:g} is negative.".FormatString(Param.BSIM3v1eta0));
                    words.Add("Warning: Eta0 = {0:g} is negative.".FormatString(Param.BSIM3v1eta0));
                }

                /* Check Abulk parameters */
                if (Math.Abs(1.0e-6 / (Param.BSIM3v1b1 + Param.BSIM3v1weff)) > 10.0)
                {
                    SpiceSharpWarning.Warning(this, "Warning: (B1 + Weff) may be too small.");
                    words.Add("Warning: (B1 + Weff) may be too small.");
                }


                /* Check Saturation parameters */
                if (Param.BSIM3v1a2 < 0.01)
                {
                    SpiceSharpWarning.Warning(this, "Warning: A2 = {0:g} is too small. Set to 0.01.".FormatString(Param.BSIM3v1a2));
                    words.Add("Warning: A2 = {0:g} is too small. Set to 0.01.".FormatString(Param.BSIM3v1a2));
                    Param.BSIM3v1a2 = 0.01;
                }
                else if (Param.BSIM3v1a2 > 1.0)
                {
                    SpiceSharpWarning.Warning(this, "Warning: A2 = {0:g} is larger than 1. A2 is set to 1 and A1 is set to 0.".FormatString(Param.BSIM3v1a2));
                    words.Add("Warning: A2 = {0:g} is larger than 1. A2 is set to 1 and A1 is set to 0.".FormatString(Param.BSIM3v1a2));
                    Param.BSIM3v1a2 = 1.0;
                    Param.BSIM3v1a1 = 0.0;

                }

                if (Param.BSIM3v1rdsw < 0.0)
                {
                    SpiceSharpWarning.Warning(this, "Warning: Rdsw = {0:g} is negative. Set to zero.".FormatString(Param.BSIM3v1rdsw));
                    words.Add("Warning: Rdsw = {0:g} is negative. Set to zero.".FormatString(Param.BSIM3v1rdsw));
                    Param.BSIM3v1rdsw = 0.0;
                    Param.BSIM3v1rds0 = 0.0;
                }
                else if ((Param.BSIM3v1rds0 > 0.0) && (Param.BSIM3v1rds0 < 0.001))
                {
                    SpiceSharpWarning.Warning(this, "Warning: Rds at current temperature = {0:g} is less than 0.001 ohm. Set to zero.".FormatString(Param.BSIM3v1rds0));
                    words.Add("Warning: Rds at current temperature = {0:g} is less than 0.001 ohm. Set to zero.".FormatString(Param.BSIM3v1rds0));
                    Param.BSIM3v1rds0 = 0.0;
                }
                if (Param.BSIM3v1vsattemp < 1.0e3)
                {
                    SpiceSharpWarning.Warning(this, "Warning: Vsat at current temperature = {0:g} may be too small.".FormatString(Param.BSIM3v1vsattemp));
                    words.Add("Warning: Vsat at current temperature = {0:g} may be too small.".FormatString(Param.BSIM3v1vsattemp));
                }

                if (Param.BSIM3v1pdibl1 < 0.0)
                {
                    SpiceSharpWarning.Warning(this, "Warning: Pdibl1 = {0:g} is negative.".FormatString(Param.BSIM3v1pdibl1));
                    words.Add("Warning: Pdibl1 = {0:g} is negative.".FormatString(Param.BSIM3v1pdibl1));
                }
                if (Param.BSIM3v1pdibl2 < 0.0)
                {
                    SpiceSharpWarning.Warning(this, "Warning: Pdibl2 = {0:g} is negative.".FormatString(Param.BSIM3v1pdibl2));
                    words.Add("Warning: Pdibl2 = {0:g} is negative.".FormatString(Param.BSIM3v1pdibl2));
                }
                /* Check overlap capacitance parameters */
                if (ModelParameters.Cgdo < 0.0)
                {
                    SpiceSharpWarning.Warning(this, "Warning: cgdo = {0:g} is negative. Set to zero.".FormatString(ModelParameters.Cgdo));
                    words.Add("Warning: cgdo = {0:g} is negative. Set to zero.".FormatString(ModelParameters.Cgdo));
                    ModelParameters.Cgdo = 0.0;
                }
                if (ModelParameters.Cgso < 0.0)
                {
                    SpiceSharpWarning.Warning(this, "Warning: cgso = {0:g} is negative. Set to zero.".FormatString(ModelParameters.Cgso));
                    words.Add("Warning: cgso = {0:g} is negative. Set to zero.".FormatString(ModelParameters.Cgso));
                    ModelParameters.Cgso = 0.0;
                }
                if (ModelParameters.Cgbo < 0.0)
                {
                    SpiceSharpWarning.Warning(this, "Warning: cgbo = {0:g} is negative. Set to zero.".FormatString(ModelParameters.Cgbo));
                    words.Add("Warning: cgbo = {0:g} is negative. Set to zero.".FormatString(ModelParameters.Cgbo));
                    ModelParameters.Cgbo = 0.0;
                }
            }

            if (words.Count > 0)
            {
                string path = ModelParameters.CheckPath;
                if (string.IsNullOrWhiteSpace(path))
                    path = "b3v1check.log";
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
