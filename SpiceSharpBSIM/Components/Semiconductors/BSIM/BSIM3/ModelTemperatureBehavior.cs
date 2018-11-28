using System;
using System.Collections.Generic;
using SpiceSharp.Behaviors;
using SpiceSharp.Simulations;
using SpiceSharp.Simulations.Behaviors;

namespace SpiceSharp.Components.BSIM3Behaviors
{
	
	/// <summary>
	/// Temperature behavior for a <see cref="BSIM3Model" />
	/// </summary>
	public class ModelTemperatureBehavior : ExportingBehavior, ITemperatureBehavior
	{
        /// <summary>
        /// Gets the model parameters.
        /// </summary>
        /// <value>
        /// The model parameters.
        /// </value>
        protected ModelBaseParameters ModelParameters { get; private set; }
		
        /// <summary>
        /// Size-dependent parameters cache
        /// </summary>
        public Dictionary<Tuple<double, double>, BSIM3SizeDependParams> SizeDependParams { get; } = new Dictionary<Tuple<double, double>, BSIM3SizeDependParams>();

		/// <summary>
		/// Properties
		/// </summary>
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

            // Get parameters
			ModelParameters = provider.GetParameterSet<ModelBaseParameters>();
		}
		
		/// <summary>
		/// Temperature behavior
		/// </summary>
		public void Temperature(BaseSimulation simulation)
		{
		    var state = simulation.RealState;
			double eg, eg0, t0, t1, delTemp, temp, tnom;

            // Update nominal temperature
		    if (!ModelParameters.Tnom.Given)
		        ModelParameters.Tnom.RawValue = state.NominalTemperature;

		    temp = state.Temperature;
            SizeDependParams.Clear();

			tnom = ModelParameters.Tnom;
			TRatio = temp / tnom;
			Vcrit = Circuit.Vt0 * Math.Log(Circuit.Vt0 / (Circuit.Root2 * 1.0e-14));
			Factor1 = Math.Sqrt(1.03594e-10 / 3.453133e-11 * ModelParameters.Tox);
			Vtm0 = 8.617087e-5 * tnom;
			eg0 = 1.16 - 7.02e-4 * tnom * tnom / (tnom + 1108.0);
			Ni = 1.45e10 * (tnom / 300.15) * Math.Sqrt(tnom / 300.15) * Math.Exp(21.5565981 - eg0 / (2.0 * Vtm0));
			Vtm = 8.617087e-5 * temp;
			eg = 1.16 - 7.02e-4 * temp * temp / (temp + 1108.0);

            // It is common that the nominal temperature matches the given temperature, this avoids some computations
			if (!temp.Equals(tnom))
			{
				t0 = eg0 / Vtm0 - eg / Vtm + ModelParameters.JctTempExponent * Math.Log(temp / tnom);
				t1 = Math.Exp(t0 / ModelParameters.JctEmissionCoeff);
				JctTempSatCurDensity = ModelParameters.JctSatCurDensity * t1;
				JctSidewallTempSatCurDensity = ModelParameters.JctSidewallSatCurDensity * t1;
			}
			else
			{
				JctTempSatCurDensity = ModelParameters.JctSatCurDensity;
				JctSidewallTempSatCurDensity = ModelParameters.JctSidewallSatCurDensity;
			}
			if (JctTempSatCurDensity < 0.0)
				JctTempSatCurDensity = 0.0;
			if (JctSidewallTempSatCurDensity < 0.0)
				JctSidewallTempSatCurDensity = 0.0;
			delTemp = state.Temperature - ModelParameters.Tnom;
			t0 = ModelParameters.Tcj * delTemp;
			if (t0 >= -1.0)
			{
				UnitAreaTempJctCap = ModelParameters.UnitAreaJctCap * (1.0 + t0);
			}
			else if (ModelParameters.UnitAreaJctCap > 0.0)
			{
				UnitAreaTempJctCap = 0.0;
				CircuitWarning.Warning(this, "Temperature effect has caused cj to be negative. Cj is clamped to zero.");
			}
			t0 = ModelParameters.Tcjsw * delTemp;
			if (t0 >= -1.0)
			{
				UnitLengthSidewallTempJctCap = ModelParameters.UnitLengthSidewallJctCap * (1.0 + t0);
			}
			else if (ModelParameters.UnitLengthSidewallJctCap > 0.0)
			{
				UnitLengthSidewallTempJctCap = 0.0;
				CircuitWarning.Warning(this, "Temperature effect has caused cjsw to be negative. Cjsw is clamped to zero.");
			}
			t0 = ModelParameters.Tcjswg * delTemp;
			if (t0 >= -1.0)
			{
				UnitLengthGateSidewallTempJctCap = ModelParameters.UnitLengthGateSidewallJctCap * (1.0 + t0);
			}
			else if (ModelParameters.UnitLengthGateSidewallJctCap > 0.0)
			{
				UnitLengthGateSidewallTempJctCap = 0.0;
				CircuitWarning.Warning(this, "Temperature effect has caused cjswg to be negative. Cjswg is clamped to zero.");
			}
			PhiB = ModelParameters.BulkJctPotential - ModelParameters.Tpb * delTemp;
			if (PhiB < 0.01)
			{
				PhiB = 0.01;
				CircuitWarning.Warning(this, "Temperature effect has caused pb to be less than 0.01. Pb is clamped to 0.01.");
			}
			PhiBSW = ModelParameters.SidewallJctPotential - ModelParameters.Tpbsw * delTemp;
			if (PhiBSW <= 0.01)
			{
				PhiBSW = 0.01;
				CircuitWarning.Warning(this, "Temperature effect has caused pbsw to be less than 0.01. Pbsw is clamped to 0.01.");
			}
			PhiBSWG = ModelParameters.GatesidewallJctPotential - ModelParameters.Tpbswg * delTemp;
			if (PhiBSWG <= 0.01)
			{
				PhiBSWG = 0.01;
				CircuitWarning.Warning(this, "Temperature effect has caused pbswg to be less than 0.01. Pbswg is clamped to 0.01.");
			}
		}
	}
}