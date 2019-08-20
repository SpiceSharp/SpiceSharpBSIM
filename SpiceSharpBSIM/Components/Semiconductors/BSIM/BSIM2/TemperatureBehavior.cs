using System;
using SpiceSharp.Behaviors;
using SpiceSharp.Simulations;

namespace SpiceSharp.Components.BSIM2Behaviors
{
	/// <summary>
	/// Temperature behavior for a <see cref="BSIM2" />
	/// </summary>
	public class TemperatureBehavior : Behavior, ITemperatureBehavior
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
		/// Properties
		/// </summary>
		public BSIM2SizeDependParams Param { get; private set; }
		public double DrainConductance { get; private set; }
		public double SourceConductance { get; private set; }
		public double Von { get; protected set; }
		
		/// <summary>
		/// Constructor
		/// </summary>
		public TemperatureBehavior(string name) : base(name)
		{
		}
		
		/// <summary>
		/// Setup the behavior
		/// </summary>
		public override void Bind(Simulation simulation, BindingContext context)
		{
            base.Bind(simulation, context);

            // Get behaviors
			ModelTemperature = context.GetBehavior<ModelTemperatureBehavior>("model");

            // Get parameter sets
			BaseParameters = context.GetParameterSet<BaseParameters>();
			ModelParameters = context.GetParameterSet<ModelBaseParameters>("model");
        }
		
		/// <summary>
		/// Temperature behavior
		/// </summary>
		void ITemperatureBehavior.Temperature()
		{
		    double effectiveLength, effectiveWidth, coxWoverL, inv_L, inv_W, tmp;
		    var size = new Tuple<double, double>(BaseParameters.Width, BaseParameters.Length);
		    if (ModelTemperature.Params.TryGetValue(size, out var sdp))
		        Param = sdp;
            else
            {
				Param = new BSIM2SizeDependParams();
				effectiveLength = BaseParameters.Length - ModelParameters.DeltaL * 1.0e-6;
				effectiveWidth = BaseParameters.Width - ModelParameters.DeltaW * 1.0e-6;
				if (effectiveLength <= 0)
				{
					throw new CircuitException("B2: mosfet {0}, model {1}: Effective channel length <=0".FormatString(ModelTemperature.Name, Name));
				}
				if (effectiveWidth <= 0)
				{
					throw new CircuitException("B2: mosfet {0}, model {1}: Effective channel width <=0".FormatString(ModelTemperature.Name, Name));
				}
				inv_L = 1.0e-6 / effectiveLength;
				inv_W = 1.0e-6 / effectiveWidth;
				Param.Width = BaseParameters.Width;
				Param.Length = BaseParameters.Length;
				Param.B2vfb = ModelParameters.Vfb0 + ModelParameters.VfbW * inv_W + ModelParameters.VfbL * inv_L;
				Param.B2phi = ModelParameters.Phi0 + ModelParameters.PhiW * inv_W + ModelParameters.PhiL * inv_L;
				Param.B2k1 = ModelParameters.K10 + ModelParameters.K1W * inv_W + ModelParameters.K1L * inv_L;
				Param.B2k2 = ModelParameters.K20 + ModelParameters.K2W * inv_W + ModelParameters.K2L * inv_L;
				Param.B2eta0 = ModelParameters.Eta00 + ModelParameters.Eta0W * inv_W + ModelParameters.Eta0L * inv_L;
				Param.B2etaB = ModelParameters.EtaB0 + ModelParameters.EtaBW * inv_W + ModelParameters.EtaBL * inv_L;
				Param.B2beta0 = ModelParameters.Mob00;
				Param.B2beta0B = ModelParameters.Mob0B0 + ModelParameters.Mob0BW * inv_W + ModelParameters.Mob0BL * inv_L;
				Param.B2betas0 = ModelParameters.Mobs00 + ModelParameters.Mobs0W * inv_W + ModelParameters.Mobs0L * inv_L;
				if (Param.B2betas0 < 1.01 * Param.B2beta0)
					Param.B2betas0 = 1.01 * Param.B2beta0;
				Param.B2betasB = ModelParameters.MobsB0 + ModelParameters.MobsBW * inv_W + ModelParameters.MobsBL * inv_L;
				tmp = (Param.B2betas0 - Param.B2beta0 - Param.B2beta0B * ModelParameters.Vbb);
				if ((-Param.B2betasB * ModelParameters.Vbb) > tmp)
					Param.B2betasB = -tmp / ModelParameters.Vbb;
				Param.B2beta20 = ModelParameters.Mob200 + ModelParameters.Mob20W * inv_W + ModelParameters.Mob20L * inv_L;
				Param.B2beta2B = ModelParameters.Mob2B0 + ModelParameters.Mob2BW * inv_W + ModelParameters.Mob2BL * inv_L;
				Param.B2beta2G = ModelParameters.Mob2G0 + ModelParameters.Mob2GW * inv_W + ModelParameters.Mob2GL * inv_L;
				Param.B2beta30 = ModelParameters.Mob300 + ModelParameters.Mob30W * inv_W + ModelParameters.Mob30L * inv_L;
				Param.B2beta3B = ModelParameters.Mob3B0 + ModelParameters.Mob3BW * inv_W + ModelParameters.Mob3BL * inv_L;
				Param.B2beta3G = ModelParameters.Mob3G0 + ModelParameters.Mob3GW * inv_W + ModelParameters.Mob3GL * inv_L;
				Param.B2beta40 = ModelParameters.Mob400 + ModelParameters.Mob40W * inv_W + ModelParameters.Mob40L * inv_L;
				Param.B2beta4B = ModelParameters.Mob4B0 + ModelParameters.Mob4BW * inv_W + ModelParameters.Mob4BL * inv_L;
				Param.B2beta4G = ModelParameters.Mob4G0 + ModelParameters.Mob4GW * inv_W + ModelParameters.Mob4GL * inv_L;
				coxWoverL = ModelParameters.Cox * effectiveWidth / effectiveLength;
				Param.B2beta0 *= coxWoverL;
				Param.B2beta0B *= coxWoverL;
				Param.B2betas0 *= coxWoverL;
				Param.B2betasB *= coxWoverL;
				Param.B2beta30 *= coxWoverL;
				Param.B2beta3B *= coxWoverL;
				Param.B2beta3G *= coxWoverL;
				Param.B2beta40 *= coxWoverL;
				Param.B2beta4B *= coxWoverL;
				Param.B2beta4G *= coxWoverL;
				Param.B2ua0 = ModelParameters.Ua00 + ModelParameters.Ua0W * inv_W + ModelParameters.Ua0L * inv_L;
				Param.B2uaB = ModelParameters.UaB0 + ModelParameters.UaBW * inv_W + ModelParameters.UaBL * inv_L;
				Param.B2ub0 = ModelParameters.Ub00 + ModelParameters.Ub0W * inv_W + ModelParameters.Ub0L * inv_L;
				Param.B2ubB = ModelParameters.UbB0 + ModelParameters.UbBW * inv_W + ModelParameters.UbBL * inv_L;
				Param.B2u10 = ModelParameters.U100 + ModelParameters.U10W * inv_W + ModelParameters.U10L * inv_L;
				Param.B2u1B = ModelParameters.U1B0 + ModelParameters.U1BW * inv_W + ModelParameters.U1BL * inv_L;
				Param.B2u1D = ModelParameters.U1D0 + ModelParameters.U1DW * inv_W + ModelParameters.U1DL * inv_L;
				Param.B2n0 = ModelParameters.N00 + ModelParameters.N0W * inv_W + ModelParameters.N0L * inv_L;
				Param.B2nB = ModelParameters.NB0 + ModelParameters.NBW * inv_W + ModelParameters.NBL * inv_L;
				Param.B2nD = ModelParameters.ND0 + ModelParameters.NDW * inv_W + ModelParameters.NDL * inv_L;
				if (Param.B2n0 < 0.0)
					Param.B2n0 = 0.0;
				Param.B2vof0 = ModelParameters.Vof00 + ModelParameters.Vof0W * inv_W + ModelParameters.Vof0L * inv_L;
				Param.B2vofB = ModelParameters.VofB0 + ModelParameters.VofBW * inv_W + ModelParameters.VofBL * inv_L;
				Param.B2vofD = ModelParameters.VofD0 + ModelParameters.VofDW * inv_W + ModelParameters.VofDL * inv_L;
				Param.B2ai0 = ModelParameters.Ai00 + ModelParameters.Ai0W * inv_W + ModelParameters.Ai0L * inv_L;
				Param.B2aiB = ModelParameters.AiB0 + ModelParameters.AiBW * inv_W + ModelParameters.AiBL * inv_L;
				Param.B2bi0 = ModelParameters.Bi00 + ModelParameters.Bi0W * inv_W + ModelParameters.Bi0L * inv_L;
				Param.B2biB = ModelParameters.BiB0 + ModelParameters.BiBW * inv_W + ModelParameters.BiBL * inv_L;
				Param.B2vghigh = ModelParameters.Vghigh0 + ModelParameters.VghighW * inv_W + ModelParameters.VghighL * inv_L;
				Param.B2vglow = ModelParameters.Vglow0 + ModelParameters.VglowW * inv_W + ModelParameters.VglowL * inv_L;
				Param.CoxWL = ModelParameters.Cox * effectiveLength * effectiveWidth * 1.0e4;
				Param.One_Third_CoxWL = Param.CoxWL / 3.0;
				Param.Two_Third_CoxWL = 2.0 * Param.One_Third_CoxWL;
				Param.B2GSoverlapCap = ModelParameters.GateSourceOverlapCap * effectiveWidth;
				Param.B2GDoverlapCap = ModelParameters.GateDrainOverlapCap * effectiveWidth;
				Param.B2GBoverlapCap = ModelParameters.GateBulkOverlapCap * effectiveLength;
				Param.SqrtPhi = Math.Sqrt(Param.B2phi);
				Param.Phis3 = Param.SqrtPhi * Param.B2phi;
				Param.Arg = Param.B2betasB - Param.B2beta0B - ModelParameters.Vdd * (Param.B2beta3B - ModelParameters.Vdd * Param.B2beta4B);
			}
			if ((DrainConductance = ModelParameters.SheetResistance * BaseParameters.DrainSquares) > 0.0)
				DrainConductance = 1.0 / DrainConductance;
			if ((SourceConductance = ModelParameters.SheetResistance * BaseParameters.SourceSquares) > 0.0)
				SourceConductance = 1.0 / SourceConductance;
			Param.B2vt0 = Param.B2vfb + Param.B2phi + Param.B2k1 * Param.SqrtPhi - Param.B2k2 * Param.B2phi;
			Von = Param.B2vt0;
		}
	}
}