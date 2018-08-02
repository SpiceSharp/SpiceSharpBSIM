using System;
using SpiceSharp.Behaviors;
using SpiceSharp.Simulations;
namespace SpiceSharp.Components.BSIM2Behaviors
{
	
	/// <summary>
	/// Temperature behavior for a <see cref="BSIM2" />
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
		/// Properties
		/// </summary>
		public BSIM2SizeDependParams Param { get; private set; }
		public double DrainConductance { get; private set; }
		public double SourceConductance { get; private set; }
		public double Von { get; set; }
		
		/// <summary>
		/// Constructor
		/// </summary>
		public TemperatureBehavior(Identifier name) : base(name)
		{
			
		}
		
		/// <summary>
		/// Setup the behavior
		/// </summary>
		public override void Setup(SetupDataProvider provider)
		{
			if (provider == null)
				throw new ArgumentNullException(nameof(provider));
			_modelTemp = provider.GetBehavior<ModelTemperatureBehavior>("model");
			_bp = provider.GetParameterSet<BaseParameters>("instance");
			_mbp = provider.GetParameterSet<ModelBaseParameters>("model");
		}
		
		/// <summary>
		/// Temperature behavior
		/// </summary>
		public override void Temperature(BaseSimulation simulation)
		{
		    double effectiveLength, effectiveWidth, coxWoverL, inv_L, inv_W, tmp;
		    var size = new Tuple<double, double>(_bp.Width, _bp.Length);
		    if (_modelTemp.Params.TryGetValue(size, out var sdp))
		        Param = sdp;
            else
            {
				Param = new BSIM2SizeDependParams();
				effectiveLength = _bp.Length - _mbp.DeltaL * 1.0e-6;
				effectiveWidth = _bp.Width - _mbp.DeltaW * 1.0e-6;
				if (effectiveLength <= 0)
				{
					throw new CircuitException("B2: mosfet {0}, model {1}: Effective channel length <=0".FormatString(_modelTemp.Name, Name));
				}
				if (effectiveWidth <= 0)
				{
					throw new CircuitException("B2: mosfet {0}, model {1}: Effective channel width <=0".FormatString(_modelTemp.Name, Name));
				}
				inv_L = 1.0e-6 / effectiveLength;
				inv_W = 1.0e-6 / effectiveWidth;
				Param.Width = _bp.Width;
				Param.Length = _bp.Length;
				Param.B2vfb = _mbp.Vfb0 + _mbp.VfbW * inv_W + _mbp.VfbL * inv_L;
				Param.B2phi = _mbp.Phi0 + _mbp.PhiW * inv_W + _mbp.PhiL * inv_L;
				Param.B2k1 = _mbp.K10 + _mbp.K1W * inv_W + _mbp.K1L * inv_L;
				Param.B2k2 = _mbp.K20 + _mbp.K2W * inv_W + _mbp.K2L * inv_L;
				Param.B2eta0 = _mbp.Eta00 + _mbp.Eta0W * inv_W + _mbp.Eta0L * inv_L;
				Param.B2etaB = _mbp.EtaB0 + _mbp.EtaBW * inv_W + _mbp.EtaBL * inv_L;
				Param.B2beta0 = _mbp.Mob00;
				Param.B2beta0B = _mbp.Mob0B0 + _mbp.Mob0BW * inv_W + _mbp.Mob0BL * inv_L;
				Param.B2betas0 = _mbp.Mobs00 + _mbp.Mobs0W * inv_W + _mbp.Mobs0L * inv_L;
				if (Param.B2betas0 < 1.01 * Param.B2beta0)
					Param.B2betas0 = 1.01 * Param.B2beta0;
				Param.B2betasB = _mbp.MobsB0 + _mbp.MobsBW * inv_W + _mbp.MobsBL * inv_L;
				tmp = (Param.B2betas0 - Param.B2beta0 - Param.B2beta0B * _mbp.Vbb);
				if ((-Param.B2betasB * _mbp.Vbb) > tmp)
					Param.B2betasB = -tmp / _mbp.Vbb;
				Param.B2beta20 = _mbp.Mob200 + _mbp.Mob20W * inv_W + _mbp.Mob20L * inv_L;
				Param.B2beta2B = _mbp.Mob2B0 + _mbp.Mob2BW * inv_W + _mbp.Mob2BL * inv_L;
				Param.B2beta2G = _mbp.Mob2G0 + _mbp.Mob2GW * inv_W + _mbp.Mob2GL * inv_L;
				Param.B2beta30 = _mbp.Mob300 + _mbp.Mob30W * inv_W + _mbp.Mob30L * inv_L;
				Param.B2beta3B = _mbp.Mob3B0 + _mbp.Mob3BW * inv_W + _mbp.Mob3BL * inv_L;
				Param.B2beta3G = _mbp.Mob3G0 + _mbp.Mob3GW * inv_W + _mbp.Mob3GL * inv_L;
				Param.B2beta40 = _mbp.Mob400 + _mbp.Mob40W * inv_W + _mbp.Mob40L * inv_L;
				Param.B2beta4B = _mbp.Mob4B0 + _mbp.Mob4BW * inv_W + _mbp.Mob4BL * inv_L;
				Param.B2beta4G = _mbp.Mob4G0 + _mbp.Mob4GW * inv_W + _mbp.Mob4GL * inv_L;
				coxWoverL = _modelTemp.Cox * effectiveWidth / effectiveLength;
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
				Param.B2ua0 = _mbp.Ua00 + _mbp.Ua0W * inv_W + _mbp.Ua0L * inv_L;
				Param.B2uaB = _mbp.UaB0 + _mbp.UaBW * inv_W + _mbp.UaBL * inv_L;
				Param.B2ub0 = _mbp.Ub00 + _mbp.Ub0W * inv_W + _mbp.Ub0L * inv_L;
				Param.B2ubB = _mbp.UbB0 + _mbp.UbBW * inv_W + _mbp.UbBL * inv_L;
				Param.B2u10 = _mbp.U100 + _mbp.U10W * inv_W + _mbp.U10L * inv_L;
				Param.B2u1B = _mbp.U1B0 + _mbp.U1BW * inv_W + _mbp.U1BL * inv_L;
				Param.B2u1D = _mbp.U1D0 + _mbp.U1DW * inv_W + _mbp.U1DL * inv_L;
				Param.B2n0 = _mbp.N00 + _mbp.N0W * inv_W + _mbp.N0L * inv_L;
				Param.B2nB = _mbp.NB0 + _mbp.NBW * inv_W + _mbp.NBL * inv_L;
				Param.B2nD = _mbp.ND0 + _mbp.NDW * inv_W + _mbp.NDL * inv_L;
				if (Param.B2n0 < 0.0)
					Param.B2n0 = 0.0;
				Param.B2vof0 = _mbp.Vof00 + _mbp.Vof0W * inv_W + _mbp.Vof0L * inv_L;
				Param.B2vofB = _mbp.VofB0 + _mbp.VofBW * inv_W + _mbp.VofBL * inv_L;
				Param.B2vofD = _mbp.VofD0 + _mbp.VofDW * inv_W + _mbp.VofDL * inv_L;
				Param.B2ai0 = _mbp.Ai00 + _mbp.Ai0W * inv_W + _mbp.Ai0L * inv_L;
				Param.B2aiB = _mbp.AiB0 + _mbp.AiBW * inv_W + _mbp.AiBL * inv_L;
				Param.B2bi0 = _mbp.Bi00 + _mbp.Bi0W * inv_W + _mbp.Bi0L * inv_L;
				Param.B2biB = _mbp.BiB0 + _mbp.BiBW * inv_W + _mbp.BiBL * inv_L;
				Param.B2vghigh = _mbp.Vghigh0 + _mbp.VghighW * inv_W + _mbp.VghighL * inv_L;
				Param.B2vglow = _mbp.Vglow0 + _mbp.VglowW * inv_W + _mbp.VglowL * inv_L;
				Param.CoxWL = _modelTemp.Cox * effectiveLength * effectiveWidth * 1.0e4;
				Param.One_Third_CoxWL = Param.CoxWL / 3.0;
				Param.Two_Third_CoxWL = 2.0 * Param.One_Third_CoxWL;
				Param.B2GSoverlapCap = _mbp.GateSourceOverlapCap * effectiveWidth;
				Param.B2GDoverlapCap = _mbp.GateDrainOverlapCap * effectiveWidth;
				Param.B2GBoverlapCap = _mbp.GateBulkOverlapCap * effectiveLength;
				Param.SqrtPhi = Math.Sqrt(Param.B2phi);
				Param.Phis3 = Param.SqrtPhi * Param.B2phi;
				Param.Arg = Param.B2betasB - Param.B2beta0B - _mbp.Vdd * (Param.B2beta3B - _mbp.Vdd * Param.B2beta4B);
			}
			if ((DrainConductance = _mbp.SheetResistance * _bp.DrainSquares) != 0.0)
			{
				DrainConductance = 1.0 / DrainConductance;
			}
			if ((SourceConductance = _mbp.SheetResistance * _bp.SourceSquares) != 0.0)
			{
				SourceConductance = 1.0 / SourceConductance;
			}
			Param.B2vt0 = Param.B2vfb + Param.B2phi + Param.B2k1 * Param.SqrtPhi - Param.B2k2 * Param.B2phi;
			Von = Param.B2vt0;
		}
	}
}