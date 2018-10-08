using System;
using SpiceSharp.Behaviors;
using SpiceSharp.Simulations;
using System.Collections.Generic;

namespace SpiceSharp.Components.BSIM3v24Behaviors
{
	
	/// <summary>
	/// Temperature behavior for a <see cref="BSIM3v24Model" />
	/// </summary>
	public class ModelTemperatureBehavior : BaseTemperatureBehavior
	{
		
		/// <summary>
		/// Necessary behaviors and parameters
		/// </summary>
		private ModelBaseParameters _mbp;

        /// <summary>
        /// Size-dependent parameters cache
        /// </summary>
        public Dictionary<Tuple<double, double>, BSIM3SizeDependParams> SizeDependParams { get; } = new Dictionary<Tuple<double, double>, BSIM3SizeDependParams>();
		
		/// <summary>
		/// Properties
		/// </summary>
		public double PSizeDependParamKnot { get; private set; }
		public double Vcrit { get; private set; }
		public double Factor1 { get; private set; }
		public double Vtm { get; private set; }
        public double Vtm0 { get; private set; }
		public double JctTempSatCurDensity { get; private set; }
		public double JctSidewallTempSatCurDensity { get; private set; }
		public double UnitAreaTempJctCap { get; private set; }
		public double UnitLengthSidewallTempJctCap { get; private set; }
		public double UnitLengthGateSidewallTempJctCap { get; private set; }
		public double PhiB { get; private set; }
		public double PhiBSW { get; private set; }
		public double PhiBSWG { get; private set; }
        public double Ni { get; private set; }
        public double TRatio { get; private set; }
		
		/// <summary>
		/// Constructor
		/// </summary>
		public ModelTemperatureBehavior(string name) : base(name)
		{
			
		}
		
		/// <summary>
		/// Setup the behavior
		/// </summary>
		public override void Setup(Simulation simulation, SetupDataProvider provider)
		{
			if (provider == null)
				throw new ArgumentNullException(nameof(provider));
			_mbp = provider.GetParameterSet<ModelBaseParameters>("model");
		}
		
		/// <summary>
		/// Temperature behavior
		/// </summary>
		public override void Temperature(BaseSimulation simulation)
		{
			double eg, eg0, t0, t1, delTemp, temp, tnom;
		    var state = simulation.RealState;

		    // Update nominal temperature
		    if (!_mbp.Tnom.Given)
		        _mbp.Tnom.RawValue = state.NominalTemperature;

		    temp = state.Temperature;
		    if (_mbp.BulkJctPotential < 0.1)
		    {
		        _mbp.BulkJctPotential.RawValue = 0.1;
		        CircuitWarning.Warning(this, "Given pb is less than 0.1. Pb is set to 0.1.");
		    }
		    if (_mbp.SidewallJctPotential < 0.1)
		    {
		        _mbp.SidewallJctPotential.RawValue = 0.1;
		        CircuitWarning.Warning(this, "Given pbsw is less than 0.1. Pbsw is set to 0.1.");
		    }
		    if (_mbp.GatesidewallJctPotential < 0.1)
		    {
		        _mbp.GatesidewallJctPotential.RawValue = 0.1;
		        CircuitWarning.Warning(this, "Given pbswg is less than 0.1. Pbswg is set to 0.1.");
		    }
			
            SizeDependParams.Clear();

			tnom = _mbp.Tnom;
			TRatio = temp / tnom;
			Vcrit = Circuit.Vt0 * Math.Log(Circuit.Vt0 / (Circuit.Root2 * 1.0e-14));
			Factor1 = Math.Sqrt(1.03594e-10 / 3.453133e-11 * _mbp.Tox);
			Vtm0 = 8.617087e-5 * tnom;
			eg0 = 1.16 - 7.02e-4 * tnom * tnom / (tnom + 1108.0);
			Ni = 1.45e10 * (tnom / 300.15) * Math.Sqrt(tnom / 300.15) * Math.Exp(21.5565981 - eg0 / (2.0 * Vtm0));
			Vtm = 8.617087e-5 * temp;
			eg = 1.16 - 7.02e-4 * temp * temp / (temp + 1108.0);
			if (temp != tnom)
			{
				t0 = eg0 / Vtm0 - eg / Vtm + _mbp.JctTempExponent * Math.Log(temp / tnom);
				t1 = Math.Exp(t0 / _mbp.JctEmissionCoeff);
				JctTempSatCurDensity = _mbp.JctSatCurDensity * t1;
				JctSidewallTempSatCurDensity = _mbp.JctSidewallSatCurDensity * t1;
			}
			else
			{
				JctTempSatCurDensity = _mbp.JctSatCurDensity;
				JctSidewallTempSatCurDensity = _mbp.JctSidewallSatCurDensity;
			}
			if (JctTempSatCurDensity < 0.0)
				JctTempSatCurDensity = 0.0;
			if (JctSidewallTempSatCurDensity < 0.0)
				JctSidewallTempSatCurDensity = 0.0;
			delTemp = state.Temperature - _mbp.Tnom;
			t0 = _mbp.Tcj * delTemp;
			if (t0 >= -1.0)
			{
				UnitAreaTempJctCap = _mbp.UnitAreaJctCap * (1.0 + t0);
			}
			else if (_mbp.UnitAreaJctCap > 0.0)
			{
				UnitAreaTempJctCap = 0.0;
			    CircuitWarning.Warning(this, "Temperature effect has caused cj to be negative. Cj is clamped to zero.");
			}
			t0 = _mbp.Tcjsw * delTemp;
			if (t0 >= -1.0)
			{
				UnitLengthSidewallTempJctCap = _mbp.UnitLengthSidewallJctCap * (1.0 + t0);
			}
			else if (_mbp.UnitLengthSidewallJctCap > 0.0)
			{
				UnitLengthSidewallTempJctCap = 0.0;
			    CircuitWarning.Warning(this, "Temperature effect has caused cjsw to be negative. Cjsw is clamped to zero.");
			}
			t0 = _mbp.Tcjswg * delTemp;
			if (t0 >= -1.0)
			{
				UnitLengthGateSidewallTempJctCap = _mbp.UnitLengthGateSidewallJctCap * (1.0 + t0);
			}
			else if (_mbp.UnitLengthGateSidewallJctCap > 0.0)
			{
				UnitLengthGateSidewallTempJctCap = 0.0;
			    CircuitWarning.Warning(this, "Temperature effect has caused cjswg to be negative. Cjswg is clamped to zero.");
			}
			PhiB = _mbp.BulkJctPotential - _mbp.Tpb * delTemp;
			if (PhiB < 0.01)
			{
				PhiB = 0.01;
			    CircuitWarning.Warning(this, "Temperature effect has caused pb to be less than 0.01. Pb is clamped to 0.01.");
			}
			PhiBSW = _mbp.SidewallJctPotential - _mbp.Tpbsw * delTemp;
			if (PhiBSW <= 0.01)
			{
				PhiBSW = 0.01;
			    CircuitWarning.Warning(this, "Temperature effect has caused pbsw to be less than 0.01. Pbsw is clamped to 0.01.");
			}
			PhiBSWG = _mbp.GatesidewallJctPotential - _mbp.Tpbswg * delTemp;
			if (PhiBSWG <= 0.01)
			{
				PhiBSWG = 0.01;
			    CircuitWarning.Warning(this, "Temperature effect has caused pbswg to be less than 0.01. Pbswg is clamped to 0.01.");
			}
		}
	}
}