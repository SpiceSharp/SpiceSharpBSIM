using System;
using SpiceSharp;
using SpiceSharp.Attributes;
using SpiceSharp.Behaviors;
using SpiceSharp.Components;
using SpiceSharp.ParameterSets;

namespace SpiceSharpBSIM.Components.Semiconductors.BSIM.BSIM2Behaviors
{
    /// <summary>
    /// Temperature behavior for a <see cref="BSIM2" />
    /// </summary>
    [BehaviorFor(typeof(BSIM2)), AddBehaviorIfNo(typeof(ITemperatureBehavior))]
    [GeneratedParameters]
    public class TemperatureBehavior : Behavior, ITemperatureBehavior, IParameterized<BaseParameters>
    {
        /// <inheritdoc />
        public BaseParameters Parameters { get; }

        protected ModelTemperatureBehavior ModelTemperature { get; }
        protected ModelParameters ModelParameters { get; }

        /// <summary>
        /// Properties
        /// </summary>
        public BSIM2SizeDependParams Param { get; private set; }
        public double DrainConductance { get; private set; }
        public double SourceConductance { get; private set; }
        public double Von { get; protected set; }

        /// <summary>
        /// Constructor
        /// </summary>
        public TemperatureBehavior(ComponentBindingContext context)
            : base(context)
        {
            Parameters = context.GetParameterSet<BaseParameters>();
            ModelTemperature = context.ModelBehaviors.GetValue<ModelTemperatureBehavior>();
            ModelParameters = context.ModelBehaviors.GetParameterSet<ModelParameters>();
        }

        /// <summary>
        /// Temperature behavior
        /// </summary>
        void ITemperatureBehavior.Temperature()
        {
            double EffectiveLength;
            double EffectiveWidth;
            double CoxWoverL, Inv_L, Inv_W, tmp;

            var size = Tuple.Create(Parameters.Width.Value, Parameters.Length.Value);
            if (ModelTemperature.SizedParams.TryGetValue(size, out var sdp))
                Param = sdp;
            else
            {
                // There are no sized parameters yet, let's create them
                Param = new BSIM2SizeDependParams();
                ModelTemperature.SizedParams.Add(size, Param);

                EffectiveLength = Parameters.Length - ModelParameters.DeltaL * 1.0e-6;
                EffectiveWidth = Parameters.Width - ModelParameters.DeltaW * 1.0e-6;

                if (EffectiveLength <= 0)
                    throw new SpiceSharpException("B2: mosfet {0}, model {1}: Effective channel length <=0".FormatString(Name, ModelTemperature.Name));
                if (EffectiveWidth <= 0)
                    throw new SpiceSharpException("B2: mosfet {0}, model {1}: Effective channel width <=0".FormatString(Name, ModelTemperature.Name));

                Inv_L = 1.0e-6 / EffectiveLength;
                Inv_W = 1.0e-6 / EffectiveWidth;
                Param.Width = Parameters.Width;
                Param.Length = Parameters.Length;
                Param.B2vfb = ModelParameters.Vfb0 + ModelParameters.VfbW * Inv_W
                   + ModelParameters.VfbL * Inv_L;
                Param.B2phi = ModelParameters.Phi0 + ModelParameters.PhiW * Inv_W
                   + ModelParameters.PhiL * Inv_L;
                Param.B2k1 = ModelParameters.K10 + ModelParameters.K1W * Inv_W
                   + ModelParameters.K1L * Inv_L;
                Param.B2k2 = ModelParameters.K20 + ModelParameters.K2W * Inv_W
                   + ModelParameters.K2L * Inv_L;
                Param.B2eta0 = ModelParameters.Eta00
                   + ModelParameters.Eta0W * Inv_W
                   + ModelParameters.Eta0L * Inv_L;
                Param.B2etaB = ModelParameters.EtaB0 + ModelParameters.EtaBW
                   * Inv_W + ModelParameters.EtaBL * Inv_L;
                Param.B2beta0 = ModelParameters.Mob00;
                Param.B2beta0B = ModelParameters.Mob0B0
                   + ModelParameters.Mob0BW * Inv_W
                   + ModelParameters.Mob0BL * Inv_L;
                Param.B2betas0 = ModelParameters.Mobs00
                   + ModelParameters.Mobs0W * Inv_W
                   + ModelParameters.Mobs0L * Inv_L;
                if (Param.B2betas0 < 1.01 * Param.B2beta0)
                    Param.B2betas0 = 1.01 * Param.B2beta0;
                Param.B2betasB = ModelParameters.MobsB0
                   + ModelParameters.MobsBW * Inv_W
                   + ModelParameters.MobsBL * Inv_L;
                tmp = (Param.B2betas0 - Param.B2beta0
                       - Param.B2beta0B * ModelParameters.Vbb);
                if ((-Param.B2betasB * ModelParameters.Vbb) > tmp)
                    Param.B2betasB = -tmp / ModelParameters.Vbb;
                Param.B2beta20 = ModelParameters.Mob200
                  + ModelParameters.Mob20W * Inv_W
                  + ModelParameters.Mob20L * Inv_L;
                Param.B2beta2B = ModelParameters.Mob2B0
                  + ModelParameters.Mob2BW * Inv_W
                  + ModelParameters.Mob2BL * Inv_L;
                Param.B2beta2G = ModelParameters.Mob2G0
                  + ModelParameters.Mob2GW * Inv_W
                  + ModelParameters.Mob2GL * Inv_L;
                Param.B2beta30 = ModelParameters.Mob300
                  + ModelParameters.Mob30W * Inv_W
                  + ModelParameters.Mob30L * Inv_L;
                Param.B2beta3B = ModelParameters.Mob3B0
                  + ModelParameters.Mob3BW * Inv_W
                  + ModelParameters.Mob3BL * Inv_L;
                Param.B2beta3G = ModelParameters.Mob3G0
                  + ModelParameters.Mob3GW * Inv_W
                  + ModelParameters.Mob3GL * Inv_L;
                Param.B2beta40 = ModelParameters.Mob400
                  + ModelParameters.Mob40W * Inv_W
                  + ModelParameters.Mob40L * Inv_L;
                Param.B2beta4B = ModelParameters.Mob4B0
                  + ModelParameters.Mob4BW * Inv_W
                  + ModelParameters.Mob4BL * Inv_L;
                Param.B2beta4G = ModelParameters.Mob4G0
                  + ModelParameters.Mob4GW * Inv_W
                  + ModelParameters.Mob4GL * Inv_L;

                CoxWoverL = ModelTemperature.Cox * EffectiveWidth / EffectiveLength;

                Param.B2beta0 *= CoxWoverL;
                Param.B2beta0B *= CoxWoverL;
                Param.B2betas0 *= CoxWoverL;
                Param.B2betasB *= CoxWoverL;
                Param.B2beta30 *= CoxWoverL;
                Param.B2beta3B *= CoxWoverL;
                Param.B2beta3G *= CoxWoverL;
                Param.B2beta40 *= CoxWoverL;
                Param.B2beta4B *= CoxWoverL;
                Param.B2beta4G *= CoxWoverL;

                Param.B2ua0 = ModelParameters.Ua00 + ModelParameters.Ua0W * Inv_W
                   + ModelParameters.Ua0L * Inv_L;
                Param.B2uaB = ModelParameters.UaB0 + ModelParameters.UaBW * Inv_W
                   + ModelParameters.UaBL * Inv_L;
                Param.B2ub0 = ModelParameters.Ub00 + ModelParameters.Ub0W * Inv_W
                   + ModelParameters.Ub0L * Inv_L;
                Param.B2ubB = ModelParameters.UbB0 + ModelParameters.UbBW * Inv_W
                   + ModelParameters.UbBL * Inv_L;
                Param.B2u10 = ModelParameters.U100 + ModelParameters.U10W * Inv_W
                   + ModelParameters.U10L * Inv_L;
                Param.B2u1B = ModelParameters.U1B0 + ModelParameters.U1BW * Inv_W
                   + ModelParameters.U1BL * Inv_L;
                Param.B2u1D = ModelParameters.U1D0 + ModelParameters.U1DW * Inv_W
                   + ModelParameters.U1DL * Inv_L;
                Param.B2n0 = ModelParameters.N00 + ModelParameters.N0W * Inv_W
                   + ModelParameters.N0L * Inv_L;
                Param.B2nB = ModelParameters.NB0 + ModelParameters.NBW * Inv_W
                   + ModelParameters.NBL * Inv_L;
                Param.B2nD = ModelParameters.ND0 + ModelParameters.NDW * Inv_W
                   + ModelParameters.NDL * Inv_L;
                if (Param.B2n0 < 0.0)
                    Param.B2n0 = 0.0;

                Param.B2vof0 = ModelParameters.Vof00
                   + ModelParameters.Vof0W * Inv_W
                   + ModelParameters.Vof0L * Inv_L;
                Param.B2vofB = ModelParameters.VofB0
                   + ModelParameters.VofBW * Inv_W
                   + ModelParameters.VofBL * Inv_L;
                Param.B2vofD = ModelParameters.VofD0
                   + ModelParameters.VofDW * Inv_W
                   + ModelParameters.VofDL * Inv_L;
                Param.B2ai0 = ModelParameters.Ai00 + ModelParameters.Ai0W * Inv_W
                   + ModelParameters.Ai0L * Inv_L;
                Param.B2aiB = ModelParameters.AiB0 + ModelParameters.AiBW * Inv_W
                   + ModelParameters.AiBL * Inv_L;
                Param.B2bi0 = ModelParameters.Bi00 + ModelParameters.Bi0W * Inv_W
                   + ModelParameters.Bi0L * Inv_L;
                Param.B2biB = ModelParameters.BiB0 + ModelParameters.BiBW * Inv_W
                   + ModelParameters.BiBL * Inv_L;
                Param.B2vghigh = ModelParameters.Vghigh0
                   + ModelParameters.VghighW * Inv_W
                   + ModelParameters.VghighL * Inv_L;
                Param.B2vglow = ModelParameters.Vglow0
                   + ModelParameters.VglowW * Inv_W
                   + ModelParameters.VglowL * Inv_L;

                Param.CoxWL = ModelTemperature.Cox * EffectiveLength
                           * EffectiveWidth * 1.0e4;
                Param.One_Third_CoxWL = Param.CoxWL / 3.0;
                Param.Two_Third_CoxWL = 2.0
                           * Param.One_Third_CoxWL;
                Param.B2GSoverlapCap = ModelParameters.GateSourceOverlapCap
                   * EffectiveWidth;
                Param.B2GDoverlapCap = ModelParameters.GateDrainOverlapCap
                   * EffectiveWidth;
                Param.B2GBoverlapCap = ModelParameters.GateBulkOverlapCap
                   * EffectiveLength;
                Param.SqrtPhi = Math.Sqrt(Param.B2phi);
                Param.Phis3 = Param.SqrtPhi
                               * Param.B2phi;
                Param.Arg = Param.B2betasB
                   - Param.B2beta0B - ModelParameters.Vdd
                           * (Param.B2beta3B - ModelParameters.Vdd
                           * Param.B2beta4B);
            }
            /* process drain series resistance */
            if ((DrainConductance = ModelParameters.SheetResistance *
                    Parameters.DrainSquares) != 0.0)
            {
                DrainConductance = 1.0 / DrainConductance;
            }

            /* process source series resistance */
            if ((SourceConductance = ModelParameters.SheetResistance *
                    Parameters.SourceSquares) != 0.0)
            {
                SourceConductance = 1.0 / SourceConductance;
            }


            Param.B2vt0 = Param.B2vfb
              + Param.B2phi
              + Param.B2k1 * Param.SqrtPhi
              - Param.B2k2 * Param.B2phi;
            Von = Param.B2vt0; /* added for initialization*/
        }
    }
}