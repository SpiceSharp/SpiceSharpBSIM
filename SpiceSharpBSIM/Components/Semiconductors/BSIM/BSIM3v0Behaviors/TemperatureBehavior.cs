using SpiceSharp;
using SpiceSharp.Attributes;
using SpiceSharp.Behaviors;
using SpiceSharp.Components;
using SpiceSharp.ParameterSets;
using SpiceSharp.Simulations;
using System;

namespace SpiceSharpBSIM.Components.Semiconductors.BSIM.BSIM3v0Behaviors
{
    /// <summary>
    /// Temperature behavior for a <see cref="BSIM3v0"/>.
    /// </summary>
    [BehaviorFor(typeof(BSIM3v0)), AddBehaviorIfNo(typeof(ITemperatureBehavior))]
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
        }

        /// <summary>
        /// Set up the device.
        /// </summary>
        private void Setup()
        {
            if (!Parameters.NqsMod.Given)
                Parameters.NqsMod = new GivenParameter<int>(ModelParameters.NqsMod, false);
        }

        /// <inheritdoc />
        void ITemperatureBehavior.Temperature()
        {
            double tmp1, tmp2, Eg, ni, T0, T1, T2, T3, Ldrn, Wdrn;
            double Temp, TRatio, Inv_L, Inv_W, Inv_LW, Vtm0, Tnom;
            var key = Tuple.Create(Parameters.W.Value, Parameters.L.Value);

            Temp = _temperature.Temperature;
            Tnom = ModelParameters.Tnom;
            TRatio = Temp / Tnom;

            if (ModelTemperature.SizeDependentProperties.TryGetValue(key, out var param))
                Param = param;
            else
            {
                Param = new SizeDependentProperties();
                ModelTemperature.SizeDependentProperties.Add(key, Param);

                Ldrn = Parameters.L;
                Wdrn = Parameters.W;

                T0 = Math.Pow(Ldrn, ModelParameters.Lln);
                T1 = Math.Pow(Wdrn, ModelParameters.Lwn);
                tmp1 = ModelParameters.Ll / T0 + ModelParameters.Lw / T1
                     + ModelParameters.Lwl / (T0 * T1);
                Param.BSIM3v0dl = ModelParameters.Lint + tmp1;
                Param.BSIM3v0dlc = ModelParameters.Dlc + tmp1;

                T2 = Math.Pow(Ldrn, ModelParameters.Wln);
                T3 = Math.Pow(Wdrn, ModelParameters.Wwn);
                tmp2 = ModelParameters.Wl / T2 + ModelParameters.Ww / T3
                     + ModelParameters.Wwl / (T2 * T3);
                Param.BSIM3v0dw = ModelParameters.Wint + tmp2;
                Param.BSIM3v0dwc = ModelParameters.Dwc + tmp2;

                Param.BSIM3v0leff = Parameters.L - 2.0 * Param.BSIM3v0dl;
                if (Param.BSIM3v0leff <= 0.0)
                    throw new SpiceSharpException("BSIM3v0: mosfet {0}, model {1}: Effective channel length <= 0".FormatString(Name, ModelTemperature.Name));

                Param.BSIM3v0weff = Parameters.W - 2.0 * Param.BSIM3v0dw;
                if (Param.BSIM3v0weff <= 0.0)
                    throw new SpiceSharpException("BSIM3v0: mosfet {0}, model {1}: Effective channel width <= 0".FormatString(Name, ModelTemperature.Name));

                Param.BSIM3v0leffCV = Parameters.L - 2.0 * Param.BSIM3v0dlc;
                if (Param.BSIM3v0leffCV <= 0.0)
                    throw new SpiceSharpException("BSIM3v0: mosfet {0}, model {1}: Effective channel length for C-V <= 0".FormatString(Name, ModelTemperature.Name));

                Param.BSIM3v0weffCV = Parameters.W - 2.0 * Param.BSIM3v0dwc;
                if (Param.BSIM3v0weffCV <= 0.0)
                    throw new SpiceSharpException("BSIM3v0: mosfet {0}, model {1}: Effective channel width for C-V <= 0".FormatString(Name, ModelTemperature));

                if (ModelParameters.BinUnit == 1)
                {
                    Inv_L = 1.0e-6 / Param.BSIM3v0leff;
                    Inv_W = 1.0e-6 / Param.BSIM3v0weff;
                    Inv_LW = 1.0e-12 / (Param.BSIM3v0leff
                           * Param.BSIM3v0weff);
                }
                else
                {
                    Inv_L = 1.0 / Param.BSIM3v0leff;
                    Inv_W = 1.0 / Param.BSIM3v0weff;
                    Inv_LW = 1.0 / (Param.BSIM3v0leff
                           * Param.BSIM3v0weff);
                }
                Param.BSIM3v0cdsc = ModelParameters.Cdsc
                                  + ModelParameters.Lcdsc * Inv_L
                                  + ModelParameters.Wcdsc * Inv_W
                                  + ModelParameters.Pcdsc * Inv_LW;
                Param.BSIM3v0cdscb = ModelParameters.Cdscb
                                   + ModelParameters.Lcdscb * Inv_L
                                   + ModelParameters.Wcdscb * Inv_W
                                   + ModelParameters.Pcdscb * Inv_LW;

                Param.BSIM3v0cdscd = ModelParameters.Cdscd
                               + ModelParameters.Lcdscd * Inv_L
                               + ModelParameters.Wcdscd * Inv_W
                               + ModelParameters.Pcdscd * Inv_LW;

                Param.BSIM3v0cit = ModelParameters.Cit
                                 + ModelParameters.Lcit * Inv_L
                                 + ModelParameters.Wcit * Inv_W
                                 + ModelParameters.Pcit * Inv_LW;
                Param.BSIM3v0nfactor = ModelParameters.Nfactor
                                     + ModelParameters.Lnfactor * Inv_L
                                     + ModelParameters.Wnfactor * Inv_W
                                     + ModelParameters.Pnfactor * Inv_LW;
                Param.BSIM3v0xj = ModelParameters.Xj
                                + ModelParameters.Lxj * Inv_L
                                + ModelParameters.Wxj * Inv_W
                                + ModelParameters.Pxj * Inv_LW;
                Param.BSIM3v0vsat = ModelParameters.Vsat
                                  + ModelParameters.Lvsat * Inv_L
                                  + ModelParameters.Wvsat * Inv_W
                                  + ModelParameters.Pvsat * Inv_LW;
                Param.BSIM3v0at = ModelParameters.At
                                + ModelParameters.Lat * Inv_L
                                + ModelParameters.Wat * Inv_W
                                + ModelParameters.Pat * Inv_LW;
                Param.BSIM3v0a0 = ModelParameters.A0
                                + ModelParameters.La0 * Inv_L
                                + ModelParameters.Wa0 * Inv_W
                                + ModelParameters.Pa0 * Inv_LW;

                Param.BSIM3v0ags = ModelParameters.Ags
                                + ModelParameters.Lags * Inv_L
                                + ModelParameters.Wags * Inv_W
                                + ModelParameters.Pags * Inv_LW;

                Param.BSIM3v0a1 = ModelParameters.A1
                                + ModelParameters.La1 * Inv_L
                                + ModelParameters.Wa1 * Inv_W
                                + ModelParameters.Pa1 * Inv_LW;
                Param.BSIM3v0a2 = ModelParameters.A2
                                + ModelParameters.La2 * Inv_L
                                + ModelParameters.Wa2 * Inv_W
                                + ModelParameters.Pa2 * Inv_LW;
                Param.BSIM3v0keta = ModelParameters.Keta
                                  + ModelParameters.Lketa * Inv_L
                                  + ModelParameters.Wketa * Inv_W
                                  + ModelParameters.Pketa * Inv_LW;
                Param.BSIM3v0nsub = ModelParameters.Nsub
                                  + ModelParameters.Lnsub * Inv_L
                                  + ModelParameters.Wnsub * Inv_W
                                  + ModelParameters.Pnsub * Inv_LW;
                Param.BSIM3v0npeak = ModelParameters.Npeak
                                   + ModelParameters.Lnpeak * Inv_L
                                   + ModelParameters.Wnpeak * Inv_W
                                   + ModelParameters.Pnpeak * Inv_LW;
                Param.BSIM3v0ngate = ModelParameters.Ngate
                                   + ModelParameters.Lngate * Inv_L
                                   + ModelParameters.Wngate * Inv_W
                                   + ModelParameters.Pngate * Inv_LW;
                Param.BSIM3v0gamma1 = ModelParameters.Gamma1
                                    + ModelParameters.Lgamma1 * Inv_L
                                    + ModelParameters.Wgamma1 * Inv_W
                                    + ModelParameters.Pgamma1 * Inv_LW;
                Param.BSIM3v0gamma2 = ModelParameters.Gamma2
                                    + ModelParameters.Lgamma2 * Inv_L
                                    + ModelParameters.Wgamma2 * Inv_W
                                    + ModelParameters.Pgamma2 * Inv_LW;
                Param.BSIM3v0vbx = ModelParameters.Vbx
                                 + ModelParameters.Lvbx * Inv_L
                                 + ModelParameters.Wvbx * Inv_W
                                 + ModelParameters.Pvbx * Inv_LW;
                Param.BSIM3v0vbm = ModelParameters.Vbm
                                 + ModelParameters.Lvbm * Inv_L
                                 + ModelParameters.Wvbm * Inv_W
                                 + ModelParameters.Pvbm * Inv_LW;
                Param.BSIM3v0xt = ModelParameters.Xt
                                 + ModelParameters.Lxt * Inv_L
                                 + ModelParameters.Wxt * Inv_W
                                 + ModelParameters.Pxt * Inv_LW;
                Param.BSIM3v0k1 = ModelParameters.K1
                                + ModelParameters.Lk1 * Inv_L
                                + ModelParameters.Wk1 * Inv_W
                                + ModelParameters.Pk1 * Inv_LW;
                Param.BSIM3v0kt1 = ModelParameters.Kt1
                                 + ModelParameters.Lkt1 * Inv_L
                                 + ModelParameters.Wkt1 * Inv_W
                                 + ModelParameters.Pkt1 * Inv_LW;
                Param.BSIM3v0kt1l = ModelParameters.Kt1l
                                  + ModelParameters.Lkt1l * Inv_L
                                  + ModelParameters.Wkt1l * Inv_W
                                  + ModelParameters.Pkt1l * Inv_LW;
                Param.BSIM3v0k2 = ModelParameters.K2
                                + ModelParameters.Lk2 * Inv_L
                                + ModelParameters.Wk2 * Inv_W
                                + ModelParameters.Pk2 * Inv_LW;
                Param.BSIM3v0kt2 = ModelParameters.Kt2
                                 + ModelParameters.Lkt2 * Inv_L
                                 + ModelParameters.Wkt2 * Inv_W
                                 + ModelParameters.Pkt2 * Inv_LW;
                Param.BSIM3v0k3 = ModelParameters.K3
                                + ModelParameters.Lk3 * Inv_L
                                + ModelParameters.Wk3 * Inv_W
                                + ModelParameters.Pk3 * Inv_LW;
                Param.BSIM3v0k3b = ModelParameters.K3b
                                 + ModelParameters.Lk3b * Inv_L
                                 + ModelParameters.Wk3b * Inv_W
                                 + ModelParameters.Pk3b * Inv_LW;
                Param.BSIM3v0w0 = ModelParameters.W0
                                + ModelParameters.Lw0 * Inv_L
                                + ModelParameters.Ww0 * Inv_W
                                + ModelParameters.Pw0 * Inv_LW;
                Param.BSIM3v0nlx = ModelParameters.Nlx
                                 + ModelParameters.Lnlx * Inv_L
                                 + ModelParameters.Wnlx * Inv_W
                                 + ModelParameters.Pnlx * Inv_LW;
                Param.BSIM3v0dvt0 = ModelParameters.Dvt0
                                  + ModelParameters.Ldvt0 * Inv_L
                                  + ModelParameters.Wdvt0 * Inv_W
                                  + ModelParameters.Pdvt0 * Inv_LW;
                Param.BSIM3v0dvt1 = ModelParameters.Dvt1
                                  + ModelParameters.Ldvt1 * Inv_L
                                  + ModelParameters.Wdvt1 * Inv_W
                                  + ModelParameters.Pdvt1 * Inv_LW;
                Param.BSIM3v0dvt2 = ModelParameters.Dvt2
                                  + ModelParameters.Ldvt2 * Inv_L
                                  + ModelParameters.Wdvt2 * Inv_W
                                  + ModelParameters.Pdvt2 * Inv_LW;
                Param.BSIM3v0dvt0w = ModelParameters.Dvt0w
                                  + ModelParameters.Ldvt0w * Inv_L
                                  + ModelParameters.Wdvt0w * Inv_W
                                  + ModelParameters.Pdvt0w * Inv_LW;
                Param.BSIM3v0dvt1w = ModelParameters.Dvt1w
                                  + ModelParameters.Ldvt1w * Inv_L
                                  + ModelParameters.Wdvt1w * Inv_W
                                  + ModelParameters.Pdvt1w * Inv_LW;
                Param.BSIM3v0dvt2w = ModelParameters.Dvt2w
                                  + ModelParameters.Ldvt2w * Inv_L
                                  + ModelParameters.Wdvt2w * Inv_W
                                  + ModelParameters.Pdvt2w * Inv_LW;
                Param.BSIM3v0drout = ModelParameters.Drout
                                   + ModelParameters.Ldrout * Inv_L
                                   + ModelParameters.Wdrout * Inv_W
                                   + ModelParameters.Pdrout * Inv_LW;
                Param.BSIM3v0dsub = ModelParameters.Dsub
                                  + ModelParameters.Ldsub * Inv_L
                                  + ModelParameters.Wdsub * Inv_W
                                  + ModelParameters.Pdsub * Inv_LW;
                Param.BSIM3v0vth0 = ModelParameters.Vth0
                                  + ModelParameters.Lvth0 * Inv_L
                                  + ModelParameters.Wvth0 * Inv_W
                                  + ModelParameters.Pvth0 * Inv_LW;
                Param.BSIM3v0ua = ModelParameters.Ua
                                + ModelParameters.Lua * Inv_L
                                + ModelParameters.Wua * Inv_W
                                + ModelParameters.Pua * Inv_LW;
                Param.BSIM3v0ua1 = ModelParameters.Ua1
                                 + ModelParameters.Lua1 * Inv_L
                                 + ModelParameters.Wua1 * Inv_W
                                 + ModelParameters.Pua1 * Inv_LW;
                Param.BSIM3v0ub = ModelParameters.Ub
                                + ModelParameters.Lub * Inv_L
                                + ModelParameters.Wub * Inv_W
                                + ModelParameters.Pub * Inv_LW;
                Param.BSIM3v0ub1 = ModelParameters.Ub1
                                 + ModelParameters.Lub1 * Inv_L
                                 + ModelParameters.Wub1 * Inv_W
                                 + ModelParameters.Pub1 * Inv_LW;
                Param.BSIM3v0uc = ModelParameters.Uc
                                + ModelParameters.Luc * Inv_L
                                + ModelParameters.Wuc * Inv_W
                                + ModelParameters.Puc * Inv_LW;
                Param.BSIM3v0uc1 = ModelParameters.Uc1
                                 + ModelParameters.Luc1 * Inv_L
                                 + ModelParameters.Wuc1 * Inv_W
                                 + ModelParameters.Puc1 * Inv_LW;
                Param.BSIM3v0u0 = ModelParameters.U0
                                + ModelParameters.Lu0 * Inv_L
                                + ModelParameters.Wu0 * Inv_W
                                + ModelParameters.Pu0 * Inv_LW;
                Param.BSIM3v0ute = ModelParameters.Ute
                                 + ModelParameters.Lute * Inv_L
                                 + ModelParameters.Wute * Inv_W
                                 + ModelParameters.Pute * Inv_LW;
                Param.BSIM3v0voff = ModelParameters.Voff
                                  + ModelParameters.Lvoff * Inv_L
                                  + ModelParameters.Wvoff * Inv_W
                                  + ModelParameters.Pvoff * Inv_LW;
                Param.BSIM3v0delta = ModelParameters.Delta
                                   + ModelParameters.Ldelta * Inv_L
                                   + ModelParameters.Wdelta * Inv_W
                                   + ModelParameters.Pdelta * Inv_LW;
                Param.BSIM3v0rdsw = ModelParameters.Rdsw
                                  + ModelParameters.Lrdsw * Inv_L
                                  + ModelParameters.Wrdsw * Inv_W
                                  + ModelParameters.Prdsw * Inv_LW;
                Param.BSIM3v0prwg = ModelParameters.Prwg
                                  + ModelParameters.Lprwg * Inv_L
                                  + ModelParameters.Wprwg * Inv_W
                                  + ModelParameters.Pprwg * Inv_LW;
                Param.BSIM3v0prwb = ModelParameters.Prwb
                                  + ModelParameters.Lprwb * Inv_L
                                  + ModelParameters.Wprwb * Inv_W
                                  + ModelParameters.Pprwb * Inv_LW;
                Param.BSIM3v0prt = ModelParameters.Prt
                                  + ModelParameters.Lprt * Inv_L
                                  + ModelParameters.Wprt * Inv_W
                                  + ModelParameters.Pprt * Inv_LW;
                Param.BSIM3v0eta0 = ModelParameters.Eta0
                                  + ModelParameters.Leta0 * Inv_L
                                  + ModelParameters.Weta0 * Inv_W
                                  + ModelParameters.Peta0 * Inv_LW;
                Param.BSIM3v0etab = ModelParameters.Etab
                                  + ModelParameters.Letab * Inv_L
                                  + ModelParameters.Wetab * Inv_W
                                  + ModelParameters.Petab * Inv_LW;
                Param.BSIM3v0pclm = ModelParameters.Pclm
                                  + ModelParameters.Lpclm * Inv_L
                                  + ModelParameters.Wpclm * Inv_W
                                  + ModelParameters.Ppclm * Inv_LW;
                Param.BSIM3v0pdibl1 = ModelParameters.Pdibl1
                                    + ModelParameters.Lpdibl1 * Inv_L
                                    + ModelParameters.Wpdibl1 * Inv_W
                                    + ModelParameters.Ppdibl1 * Inv_LW;
                Param.BSIM3v0pdibl2 = ModelParameters.Pdibl2
                                    + ModelParameters.Lpdibl2 * Inv_L
                                    + ModelParameters.Wpdibl2 * Inv_W
                                    + ModelParameters.Ppdibl2 * Inv_LW;
                Param.BSIM3v0pdiblb = ModelParameters.Pdiblb
                                    + ModelParameters.Lpdiblb * Inv_L
                                    + ModelParameters.Wpdiblb * Inv_W
                                    + ModelParameters.Ppdiblb * Inv_LW;
                Param.BSIM3v0pscbe1 = ModelParameters.Pscbe1
                                    + ModelParameters.Lpscbe1 * Inv_L
                                    + ModelParameters.Wpscbe1 * Inv_W
                                    + ModelParameters.Ppscbe1 * Inv_LW;
                Param.BSIM3v0pscbe2 = ModelParameters.Pscbe2
                                    + ModelParameters.Lpscbe2 * Inv_L
                                    + ModelParameters.Wpscbe2 * Inv_W
                                    + ModelParameters.Ppscbe2 * Inv_LW;
                Param.BSIM3v0pvag = ModelParameters.Pvag
                                  + ModelParameters.Lpvag * Inv_L
                                  + ModelParameters.Wpvag * Inv_W
                                  + ModelParameters.Ppvag * Inv_LW;
                Param.BSIM3v0wr = ModelParameters.Wr
                                + ModelParameters.Lwr * Inv_L
                                + ModelParameters.Wwr * Inv_W
                                + ModelParameters.Pwr * Inv_LW;
                Param.BSIM3v0dwg = ModelParameters.Dwg
                                 + ModelParameters.Ldwg * Inv_L
                                 + ModelParameters.Wdwg * Inv_W
                                 + ModelParameters.Pdwg * Inv_LW;
                Param.BSIM3v0dwb = ModelParameters.Dwb
                                 + ModelParameters.Ldwb * Inv_L
                                 + ModelParameters.Wdwb * Inv_W
                                 + ModelParameters.Pdwb * Inv_LW;
                Param.BSIM3v0b0 = ModelParameters.B0
                                + ModelParameters.Lb0 * Inv_L
                                + ModelParameters.Wb0 * Inv_W
                                + ModelParameters.Pb0 * Inv_LW;
                Param.BSIM3v0b1 = ModelParameters.B1
                                + ModelParameters.Lb1 * Inv_L
                                + ModelParameters.Wb1 * Inv_W
                                + ModelParameters.Pb1 * Inv_LW;
                Param.BSIM3v0alpha0 = ModelParameters.Alpha0
                                    + ModelParameters.Lalpha0 * Inv_L
                                    + ModelParameters.Walpha0 * Inv_W
                                    + ModelParameters.Palpha0 * Inv_LW;
                Param.BSIM3v0beta0 = ModelParameters.Beta0
                                   + ModelParameters.Lbeta0 * Inv_L
                                   + ModelParameters.Wbeta0 * Inv_W
                                   + ModelParameters.Pbeta0 * Inv_LW;
                /* CV model */
                Param.BSIM3v0elm = ModelParameters.Elm
                                + ModelParameters.Lelm * Inv_L
                                + ModelParameters.Welm * Inv_W
                                + ModelParameters.Pelm * Inv_LW;
                Param.BSIM3v0cgsl = ModelParameters.Cgsl
                                  + ModelParameters.Lcgsl * Inv_L
                                  + ModelParameters.Wcgsl * Inv_W
                                  + ModelParameters.Pcgsl * Inv_LW;
                Param.BSIM3v0cgdl = ModelParameters.Cgdl
                                  + ModelParameters.Lcgdl * Inv_L
                                  + ModelParameters.Wcgdl * Inv_W
                                  + ModelParameters.Pcgdl * Inv_LW;
                Param.BSIM3v0ckappa = ModelParameters.Ckappa
                                    + ModelParameters.Lckappa * Inv_L
                                    + ModelParameters.Wckappa * Inv_W
                                    + ModelParameters.Pckappa * Inv_LW;
                Param.BSIM3v0cf = ModelParameters.Cf
                                + ModelParameters.Lcf * Inv_L
                                + ModelParameters.Wcf * Inv_W
                                + ModelParameters.Pcf * Inv_LW;
                Param.BSIM3v0clc = ModelParameters.Clc
                                 + ModelParameters.Lclc * Inv_L
                                 + ModelParameters.Wclc * Inv_W
                                 + ModelParameters.Pclc * Inv_LW;
                Param.BSIM3v0cle = ModelParameters.Cle
                                 + ModelParameters.Lcle * Inv_L
                                 + ModelParameters.Wcle * Inv_W
                                 + ModelParameters.Pcle * Inv_LW;
                Param.BSIM3v0abulkCVfactor = 1.0 + Math.Pow((Param.BSIM3v0clc
                                           / Param.BSIM3v0leff),
                                           Param.BSIM3v0cle);

                Param.BSIM3v0cgdo = (ModelParameters.Cgdo + Param.BSIM3v0cf)
                                  * Param.BSIM3v0weffCV;
                Param.BSIM3v0cgso = (ModelParameters.Cgso + Param.BSIM3v0cf)
                                  * Param.BSIM3v0weffCV;
                Param.BSIM3v0cgbo = ModelParameters.Cgbo * Param.BSIM3v0leffCV;

                T0 = (TRatio - 1.0);
                Param.BSIM3v0ua = Param.BSIM3v0ua + Param.BSIM3v0ua1 * T0;
                Param.BSIM3v0ub = Param.BSIM3v0ub + Param.BSIM3v0ub1 * T0;
                Param.BSIM3v0uc = Param.BSIM3v0uc + Param.BSIM3v0uc1 * T0;

                Param.BSIM3v0u0temp = Param.BSIM3v0u0
                                    * Math.Pow(TRatio, Param.BSIM3v0ute);
                Param.BSIM3v0vsattemp = Param.BSIM3v0vsat - Param.BSIM3v0at
                                      * T0;
                Param.BSIM3v0rds0 = (Param.BSIM3v0rdsw + Param.BSIM3v0prt * T0)
                                  / Math.Pow(Param.BSIM3v0weff * 1E6, Param.BSIM3v0wr);

                if (!ModelParameters.Npeak.Given && ModelParameters.Gamma1.Given)
                {
                    T0 = Param.BSIM3v0gamma1 * ModelTemperature.Cox;
                    Param.BSIM3v0npeak = 3.021E22 * T0 * T0;
                }

                Vtm0 = KboQ * Tnom;
                Eg = 1.16 - 7.02e-4 * Tnom * Tnom / (Tnom + 1108.0);
                ni = 1.45e10 * (Tnom / 300.15) * Math.Sqrt(Tnom / 300.15)
                   * Math.Exp(21.5565981 - Eg / (2.0 * Vtm0));

                Param.BSIM3v0phi = 2.0 * Vtm0
                                 * Math.Log(Param.BSIM3v0npeak / ni);

                Param.BSIM3v0sqrtPhi = Math.Sqrt(Param.BSIM3v0phi);
                Param.BSIM3v0phis3 = Param.BSIM3v0sqrtPhi * Param.BSIM3v0phi;

                Param.BSIM3v0Xdep0 = Math.Sqrt(2.0 * EPSSI / (Charge_q
                                   * Param.BSIM3v0npeak * 1.0e6))
                                   * Param.BSIM3v0sqrtPhi;
                Param.BSIM3v0sqrtXdep0 = Math.Sqrt(Param.BSIM3v0Xdep0);
                Param.BSIM3v0litl = Math.Sqrt(3.0 * Param.BSIM3v0xj
                                  * ModelParameters.Tox);
                Param.BSIM3v0vbi = Vtm0 * Math.Log(1.0e20
                                 * Param.BSIM3v0npeak / (ni * ni));
                Param.BSIM3v0cdep0 = Math.Sqrt(Charge_q * EPSSI
                                   * Param.BSIM3v0npeak * 1.0e6 / 2.0
                                   / Param.BSIM3v0phi);

                if (ModelParameters.K1.Given || ModelParameters.K2.Given)
                {
                    if (!ModelParameters.K1.Given)
                    {
                        SpiceSharpWarning.Warning(this, "Warning: k1 should be specified with k2.");
                        Param.BSIM3v0k1 = 0.53;
                    }
                    if (!ModelParameters.K2.Given)
                    {
                        SpiceSharpWarning.Warning(this, "Warning: k2 should be specified with k1.");
                        Param.BSIM3v0k2 = -0.0186;
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
                        Param.BSIM3v0vbx = Param.BSIM3v0phi - 7.7348e-4
                                         * Param.BSIM3v0npeak
                                         * Param.BSIM3v0xt * Param.BSIM3v0xt;
                    if (Param.BSIM3v0vbx > 0.0)
                        Param.BSIM3v0vbx = -Param.BSIM3v0vbx;
                    if (Param.BSIM3v0vbm > 0.0)
                        Param.BSIM3v0vbm = -Param.BSIM3v0vbm;

                    if (!ModelParameters.Gamma1.Given)
                        Param.BSIM3v0gamma1 = 5.753e-12
                                            * Math.Sqrt(Param.BSIM3v0npeak)
                                            / ModelTemperature.Cox;
                    if (!ModelParameters.Gamma2.Given)
                        Param.BSIM3v0gamma2 = 5.753e-12
                                            * Math.Sqrt(Param.BSIM3v0nsub)
                                            / ModelTemperature.Cox;

                    T0 = Param.BSIM3v0gamma1 - Param.BSIM3v0gamma2;
                    T1 = Math.Sqrt(Param.BSIM3v0phi - Param.BSIM3v0vbx)
                       - Param.BSIM3v0sqrtPhi;
                    T2 = Math.Sqrt(Param.BSIM3v0phi * (Param.BSIM3v0phi
                       - Param.BSIM3v0vbm)) - Param.BSIM3v0phi;
                    Param.BSIM3v0k2 = T0 * T1 / (2.0 * T2 + Param.BSIM3v0vbm);
                    Param.BSIM3v0k1 = Param.BSIM3v0gamma2 - 2.0
                                    * Param.BSIM3v0k2 * Math.Sqrt(Param.BSIM3v0phi
                                    - Param.BSIM3v0vbm);
                }

                if (Param.BSIM3v0k2 > 0.0)
                {
                    T0 = 0.5 * Param.BSIM3v0k1 / Param.BSIM3v0k2;
                    Param.BSIM3v0vbsc = 0.9 * (Param.BSIM3v0phi - T0 * T0);
                    if (Param.BSIM3v0vbsc > -3.0)
                        Param.BSIM3v0vbsc = -3.0;
                    else if (Param.BSIM3v0vbsc < -30.0)
                        Param.BSIM3v0vbsc = -30.0;
                }
                else
                {
                    Param.BSIM3v0vbsc = -10.0;
                }

                if (ModelParameters.Vth0.Given)
                    Param.BSIM3v0vfb = ModelParameters.Type * Param.BSIM3v0vth0
                                     - Param.BSIM3v0phi - Param.BSIM3v0k1
                                     * Param.BSIM3v0sqrtPhi;
                else
                    Param.BSIM3v0vth0 = ModelParameters.Type * (-1.0
                                      + Param.BSIM3v0phi + Param.BSIM3v0k1
                                      * Param.BSIM3v0sqrtPhi);

                T1 = Math.Sqrt(EPSSI / EPSOX * ModelParameters.Tox
                   * Param.BSIM3v0Xdep0);
                T0 = Math.Exp(-0.5 * Param.BSIM3v0dsub * Param.BSIM3v0leff / T1);
                Param.BSIM3v0theta0vb0 = (T0 + 2.0 * T0 * T0);

                T0 = Math.Exp(-0.5 * Param.BSIM3v0drout * Param.BSIM3v0leff / T1);
                T2 = (T0 + 2.0 * T0 * T0);
                Param.BSIM3v0thetaRout = Param.BSIM3v0pdibl1 * T2
                                       + Param.BSIM3v0pdibl2;

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
            }
            this._cgso = Param.BSIM3v0cgso;
            this._cgdo = Param.BSIM3v0cgdo;
        }
    }
}
