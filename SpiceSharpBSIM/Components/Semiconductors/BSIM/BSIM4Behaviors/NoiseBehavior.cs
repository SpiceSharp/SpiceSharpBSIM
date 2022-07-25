using SpiceSharp;
using SpiceSharp.Attributes;
using SpiceSharp.Behaviors;
using SpiceSharp.Components;
using SpiceSharp.Components.NoiseSources;
using SpiceSharp.Simulations;
using System;

namespace SpiceSharpBSIM.Components.Semiconductors.BSIM.BSIM4Behaviors
{
    /// <summary>
    /// Noise behavior for a <see cref="BSIM4"/>.
    /// </summary>
    [BehaviorFor(typeof(BSIM4)), AddBehaviorIfNo(typeof(INoiseBehavior))]
    public class NoiseBehavior : FrequencyBehavior, INoiseBehavior
    {
        private readonly INoiseSimulationState _state;
        private readonly ITemperatureSimulationState _temperature;
        private readonly NoiseThermal _rd, _rs, _id, _rg, _rbps, _rbpd, _rbpb, _rbsb, _rbdb;
        private readonly NoiseShot _nigs, _nigd, _nigb;
        private readonly NoiseGain _fl;
        private CorrelatedNoiseThermal _corl;

        /// <inheritdoc/>
        [ParameterName("noise"), ParameterInfo("The total output noise density")]
        public double OutputNoiseDensity =>
            _rd.OutputNoiseDensity +
            _rs.OutputNoiseDensity +
            _rg.OutputNoiseDensity +
            (_rbps?.OutputNoiseDensity ?? 0.0) +
            (_rbpd?.OutputNoiseDensity ?? 0.0) +
            (_rbpb?.OutputNoiseDensity ?? 0.0) +
            (_rbsb?.OutputNoiseDensity ?? 0.0) +
            (_rbdb?.OutputNoiseDensity ?? 0.0) +
            _id.OutputNoiseDensity + 
            _fl.OutputNoiseDensity +
            _nigs.OutputNoiseDensity + 
            _nigd.OutputNoiseDensity +
            _nigb.OutputNoiseDensity +
            (_corl?.OutputNoiseDensity ?? 0.0);

        /// <inheritdoc/>
        [ParameterName("onoise"), ParameterInfo("The total integrated output noise")]
        public double TotalOutputNoise =>
            _rd.TotalOutputNoise +
            _rs.TotalOutputNoise +
            _rg.TotalOutputNoise +
            (_rbps?.TotalOutputNoise ?? 0.0) +
            (_rbpd?.TotalOutputNoise ?? 0.0) +
            (_rbpb?.TotalOutputNoise ?? 0.0) +
            (_rbsb?.TotalOutputNoise ?? 0.0) +
            (_rbdb?.TotalOutputNoise ?? 0.0) +
            _id.TotalOutputNoise +
            _fl.TotalOutputNoise +
            _nigs.TotalOutputNoise +
            _nigd.TotalOutputNoise +
            _nigb.TotalOutputNoise +
            (_corl?.TotalOutputNoise ?? 0.0);

        /// <inheritdoc/>
        [ParameterName("inoise"), ParameterInfo("The total integrated input noise")]
        public double TotalInputNoise =>
            _rd.TotalInputNoise +
            _rs.TotalInputNoise +
            _rg.TotalInputNoise +
            (_rbps?.TotalInputNoise ?? 0.0) +
            (_rbpd?.TotalInputNoise ?? 0.0) +
            (_rbpb?.TotalInputNoise ?? 0.0) +
            (_rbsb?.TotalInputNoise ?? 0.0) +
            (_rbdb?.TotalInputNoise ?? 0.0) +
            _id.TotalInputNoise +
            _fl.TotalInputNoise +
            _nigs.TotalInputNoise +
            _nigd.TotalInputNoise +
            _nigb.TotalInputNoise +
            (_corl?.TotalInputNoise ?? 0.0);

        [ParameterName("rd"), ParameterInfo("The noise due to rd.")]
        public INoiseSource Rd => _rd;
        [ParameterName("rs"), ParameterInfo("The noise due to rs.")]
        public INoiseSource Rs => _rs;
        [ParameterName("rg"), ParameterInfo("The noise due to rgeltd.")]
        public INoiseSource Rg => _rg;
        [ParameterName("rbps"), ParameterInfo("The noise due to rbps.")]
        public INoiseSource Rbps => _rbps;
        [ParameterName("rbpd"), ParameterInfo("The noise due to rbpd.")]
        public INoiseSource Rbpd => _rbpd;
        [ParameterName("rbpb"), ParameterInfo("The noise due to rbpb.")]
        public INoiseSource Rbpb => _rbpb;
        [ParameterName("rbsb"), ParameterInfo("The noise due to rbsb.")]
        public INoiseSource Rbsb => _rbsb;
        [ParameterName("rbdb"), ParameterInfo("The noise due to rbdb.")]
        public INoiseSource Rbdb => _rbdb;
        [ParameterName("id"), ParameterInfo("Noise due to id (for tnoiMod2: uncorrelated portion only).")]
        public INoiseSource Id => _id;
        [ParameterName("1overf"), ParameterInfo("Flicker (1/f) noise.")]
        public INoiseSource Flicker => _fl;
        [ParameterName("igs"), ParameterInfo("The shot noise due to igs.")]
        public INoiseSource NoiseIgs => _nigs;
        [ParameterName("igd"), ParameterInfo("The shot noise due to igd.")]
        public INoiseSource NoiseIgd => _nigd;
        [ParameterName("igb"), ParameterInfo("The shot noise due to igb.")]
        public INoiseSource NoiseIgb => _nigb;
        [ParameterName("corl"), ParameterInfo("The contribution of correlated drain and induced gate noise.")]
        public INoiseSource Corl => _corl;

        /// <summary>
        /// Creates a new <see cref="NoiseBehavior"/>.
        /// </summary>
        /// <param name="context">The context.</param>
        public NoiseBehavior(ComponentBindingContext context)
            : base(context)
        {
            _state = context.GetState<INoiseSimulationState>();
            _temperature = context.GetState<ITemperatureSimulationState>();
            _rd = new NoiseThermal("rd", _drainPrime, _drain);
            _rs = new NoiseThermal("rs", _sourcePrime, _source);
            if (Parameters.RgateMod.Value == 1 || Parameters.RgateMod.Value == 2)
                _rg = new NoiseThermal("rg", _gatePrime, _gate);
            else if (Parameters.RgateMod.Value == 3)
                _rg = new NoiseThermal("rg", _gateMid, _gate);
            if (Parameters.RbodyMod.Value != 0)
            {
                _rbps = new NoiseThermal("rbps", _bulkPrime, _sourceBulk);
                _rbpd = new NoiseThermal("rbpd", _bulkPrime, _drainBulk);
                _rbpb = new NoiseThermal("rbpb", _bulkPrime, _bulk);
                _rbsb = new NoiseThermal("rbsb", _bulk, _sourceBulk);
                _rbdb = new NoiseThermal("rbdb", _bulk, _drainBulk);
            }
            _id = new NoiseThermal("id", _drainPrime, _sourcePrime);
            _fl = new NoiseGain("1overf", _drainPrime, _sourcePrime);

            _nigs = new NoiseShot("igs", _gatePrime, _sourcePrime);
            _nigd = new NoiseShot("igd", _gatePrime, _drainPrime);
            _nigb = new NoiseShot("igb", _gatePrime, _bulkPrime);
        }

        /// <inheritdoc />
        void INoiseSource.Initialize()
        {
            _rd.Initialize();
            _rs.Initialize();
            _rg?.Initialize();
            _rbps?.Initialize();
            _rbpd?.Initialize();
            _rbpb?.Initialize();
            _rbsb?.Initialize();
            _rbdb?.Initialize();
            _id.Initialize();
            _fl.Initialize();
            _nigs.Initialize();
            _nigd.Initialize();
            _nigb.Initialize();
            if (this._mode > 0)
                _corl = new CorrelatedNoiseThermal("corl", _drainPrime, _sourcePrime, _gatePrime, _sourcePrime);
            else
                _corl = new CorrelatedNoiseThermal("corl", _sourcePrime, _drainPrime, _gatePrime, _drainPrime);
        }

        /// <inheritdoc />
        void INoiseBehavior.Compute()
        {
            double T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11;
            double Vds, Ssi, Swi;
            double tmp = 0.0, gdpr, gspr, npart_theta = 0.0, npart_beta = 0.0, igsquare, bodymode;

            /* tnoiMod=2 (v4.7) */
            double eta, Leff, Lvsat, gamma, delta, epsilon, GammaGd0 = 0.0;
            double npart_c, sigrat = 0.0, C0, omega, ctnoi = 0.0;

            double m = Parameters.M;

            if (ModelParameters.TnoiMod == 0)
            {
                if (ModelParameters.RdsMod == 0)
                {
                    gspr = this._sourceConductance;
                    gdpr = this._drainConductance;
                    if (this._grdsw > 0.0)
                        tmp = 1.0 / this._grdsw; /* tmp used below */
                    else
                        tmp = 0.0;
                }
                else
                {
                    gspr = this._gstot;
                    gdpr = this._gdtot;
                    tmp = 0.0;
                }
            }
            else if (ModelParameters.TnoiMod == 1)
            {
                T5 = this._vgsteff / this._esatL;
                T5 *= T5;
                npart_beta = ModelParameters.Rnoia * (1.0 + T5
                           * ModelParameters.Tnoia * Param.BSIM4leff);
                npart_theta = ModelParameters.Rnoib * (1.0 + T5
                            * ModelParameters.Tnoib * Param.BSIM4leff);
                if (npart_theta > 0.9)
                    npart_theta = 0.9;
                if (npart_theta > 0.9 * npart_beta)
                    npart_theta = 0.9 * npart_beta; //4.6.2

                if (ModelParameters.RdsMod == 0)
                {
                    gspr = this._sourceConductance;
                    gdpr = this._drainConductance;
                }
                else
                {
                    gspr = this._gstot;
                    gdpr = this._gdtot;
                }

                if (this._vds >= 0.0)
                    gspr = gspr * (1.0 + npart_theta * npart_theta * gspr
                         / this._idovVds);
                else
                    gdpr = gdpr * (1.0 + npart_theta * npart_theta * gdpr
                         / this._idovVds);
            }
            else
            {   /* tnoiMod=2 (v4.7) */

                if (ModelParameters.RdsMod == 0)
                {
                    gspr = this._sourceConductance;
                    gdpr = this._drainConductance;
                }
                else
                {
                    gspr = this._gstot;
                    gdpr = this._gdtot;
                }

            }

            _rd.Compute(gdpr * m, _temperature.Temperature);
            _rs.Compute(gspr * m, _temperature.Temperature);

            if (Parameters.RgateMod == 1)
                _rg.Compute(this._grgeltd * m, _temperature.Temperature);
            else if (Parameters.RgateMod == 2)
            {
                T0 = 1.0 + this._grgeltd / this._gcrg;
                T1 = T0 * T0;
                _rg.Compute(this._grgeltd * m / T1, _temperature.Temperature);
            }
            else if (Parameters.RgateMod == 3)
                _rg.Compute(this._grgeltd * m, _temperature.Temperature);
            else
                _rg?.Compute(0.0, _temperature.Temperature);

            bodymode = 5;
            if (Parameters.RbodyMod == 2)
            {
                if ((!ModelParameters.Rbps0.Given) ||
              (!ModelParameters.Rbpd0.Given))
                    bodymode = 1;
                else
             if ((!ModelParameters.Rbsbx0.Given && !ModelParameters.Rbsby0.Given) ||
                  (!ModelParameters.Rbdbx0.Given && !ModelParameters.Rbdby0.Given))
                    bodymode = 3;
            }

            if (Parameters.RbodyMod.Value != 0)
            {
                if (bodymode == 5)
                {
                    _rbps.Compute(this._grbps * m, _temperature.Temperature);
                    _rbpd.Compute(this._grbpd * m, _temperature.Temperature);
                    _rbpb.Compute(this._grbpb * m, _temperature.Temperature);
                    _rbsb.Compute(this._grbsb * m, _temperature.Temperature);
                    _rbdb.Compute(this._grbdb * m, _temperature.Temperature);
                }
                if (bodymode == 3)
                {
                    _rbps.Compute(this._grbps * m, _temperature.Temperature);
                    _rbpd.Compute(this._grbpd * m, _temperature.Temperature);
                    _rbpb.Compute(this._grbpb * m, _temperature.Temperature);
                    _rbsb.Compute(0.0, _temperature.Temperature);
                    _rbdb.Compute(0.0, _temperature.Temperature);
                }
                if (bodymode == 1)
                {
                    _rbpb.Compute(this._grbpb * m, _temperature.Temperature);
                    _rbps.Compute(0.0, _temperature.Temperature);
                    _rbsb.Compute(0.0, _temperature.Temperature);
                    _rbpd.Compute(0.0, _temperature.Temperature);
                    _rbdb.Compute(0.0, _temperature.Temperature);
                }
            }

            if (ModelParameters.TnoiMod == 2)
            {
                eta = 1.0 - this._vdseff * this._abovVgst2Vtm;
                T0 = 1.0 - eta;
                T1 = 1.0 + eta;
                T2 = T1 + 2.0 * this._abulk * ModelTemperature.Vtm / this._vgsteff;
                Leff = Param.BSIM4leff;
                Lvsat = Leff * (1.0 + this._vdseff / this._esatL);
                T6 = Leff / Lvsat;
                /*Unwanted code for T5 commented*/
                /*T5 = this._vgsteff / this._esatL;
                T5 = T5 * T5;
                */
                gamma = T6 * (0.5 * T1 + T0 * T0 / (6.0 * T2));
                T3 = T2 * T2;
                T4 = T0 * T0;
                T5 = T3 * T3;
                delta = (T1 / T3 - (5.0 * T1 + T2) * T4 / (15.0 * T5) + T4 * T4 / (9.0 * T5 * T2)) / (6.0 * T6 * T6 * T6);
                T7 = T0 / T2;
                epsilon = (T7 - T7 * T7 * T7 / 3.0) / (6.0 * T6);

                T8 = this._vgsteff / this._esatL;
                T8 *= T8;
                if (ModelParameters.Version.Value != "4.8.1" && ModelParameters.Version.Value != "4.81")
                {
                    npart_c = ModelParameters.Rnoic * (1.0 + T8
                            * ModelParameters.Tnoic * Leff);
                    ctnoi = epsilon / Math.Sqrt(gamma * delta)
                        * (2.5316 * npart_c);

                    npart_beta = ModelParameters.Rnoia * (1.0 + T8
                        * ModelParameters.Tnoia * Leff);
                    npart_theta = ModelParameters.Rnoib * (1.0 + T8
                        * ModelParameters.Tnoib * Leff);
                    gamma = gamma * (3.0 * npart_beta * npart_beta);
                    delta = delta * (3.75 * npart_theta * npart_theta);

                    GammaGd0 = gamma * this._noiGd0;
                    C0 = this._coxeff * Param.BSIM4weffCV * Parameters.Nf * Param.BSIM4leffCV;
                    T0 = C0 / this._noiGd0;
                    sigrat = T0 * Math.Sqrt(delta / gamma);
                }
                else
                {
                    npart_c = ModelParameters.Rnoic * (1.0 + T8
                           * ModelParameters.Tnoic * Leff);
                    /* Limits added for rnoia, rnoib, rnoic, tnoia, tnoib and tnoic in BSIM4.8.1 */
                    T9 = gamma * delta;
                    if (T9 > 0)
                        ctnoi = epsilon / Math.Sqrt(gamma * delta) * (2.5316 * npart_c);
                    else
                        ctnoi = 1.0;
                    if (ctnoi > 1)
                        ctnoi = 1;
                    if (ctnoi < 0)
                        ctnoi = 0;

                    npart_beta = ModelParameters.Rnoia * (1.0 + T8
                        * ModelParameters.Tnoia * Leff);
                    npart_theta = ModelParameters.Rnoib * (1.0 + T8
                        * ModelParameters.Tnoib * Leff);
                    gamma = gamma * (3.0 * npart_beta * npart_beta);
                    delta = delta * (3.75 * npart_theta * npart_theta);

                    GammaGd0 = gamma * this._noiGd0;
                    C0 = this._coxeff * Param.BSIM4weffCV * Parameters.Nf * Param.BSIM4leffCV;
                    T0 = C0 / this._noiGd0;

                    if (gamma > 0 && delta > 0)
                        sigrat = T0 * Math.Sqrt(delta / gamma);
                    else
                        sigrat = 0.0;
                }
            }

            switch (ModelParameters.TnoiMod)
            {
                case 0:
                    if (ModelParameters.Version.Value != "4.8.1" && ModelParameters.Version.Value != "4.81")
                    {
                        T0 = this._ueff * Math.Abs(this._qinv);
                        T1 = T0 * tmp + Param.BSIM4leff
                           * Param.BSIM4leff;
                        _id.Compute((T0 / T1) * ModelParameters.Ntnoi * m, _temperature.Temperature);
                    }
                    else
                    {
                        T0 = this._ueff * Math.Abs(this._qinv);
                        T1 = T0 * tmp + Param.BSIM4leff
                                * Param.BSIM4leff;
                        _id.Compute((T0 / T1) * ModelParameters.Ntnoi * m, _temperature.Temperature);
                        _corl?.Compute(0.0, 0.0, 0.0, _temperature.Temperature);
                    }
                    break;
                case 1:
                    if (ModelParameters.Version.Value != "4.8.1" && ModelParameters.Version.Value != "4.81")
                    {
                        T0 = this._gm + this._gmbs + this._gds;
                        T0 *= T0;
                        igsquare = npart_theta * npart_theta * T0 / this._idovVds;
                        T1 = npart_beta * (this._gm
                           + this._gmbs) + this._gds;
                        T2 = T1 * T1 / this._idovVds;
                        _id.Compute((T2 - igsquare) * m, _temperature.Temperature);
                    }
                    else
                    {
                        T0 = this._gm + this._gmbs + this._gds;
                        T0 *= T0;
                        igsquare = npart_theta * npart_theta * T0 / this._idovVds;
                        T1 = npart_beta * (this._gm
                        + this._gmbs) + this._gds;
                        T2 = T1 * T1 / this._idovVds;
                        _id.Compute((T2 - igsquare) * m, _temperature.Temperature);
                        _corl?.Compute(0.0, 0.0, 0.0, _temperature.Temperature);
                    }
                    break;
                case 2:
                    T2 = GammaGd0;
                    T3 = ctnoi * ctnoi;
                    T4 = 1.0 - T3;
                    _id.Compute(T2 * T4 * m, _temperature.Temperature);

                    /* Evaluate output noise due to two correlated noise sources */
                    omega = 2.0 * Math.PI * _state.Point.Value.Frequency;
                    T5 = omega * sigrat;
                    T6 = T5 * T5;
                    T7 = T6 / (1.0 + T6);

                    // We have already reconnected the source at this point
                    // We don't need to discrimate based on Mode - Sven Boulanger 20220725
                    _corl?.Compute(T2 * T3 * m, T2 * T7 * m, 0.5 * Math.PI, _temperature.Temperature);
                    break;
            }

            switch (ModelParameters.FnoiMod)
            {
                case 0:
                    _fl.Compute(m * ModelParameters.Kf
                          * Math.Exp(ModelParameters.Af
                          * Math.Log(Math.Max(Math.Abs(this._cd),
                          1e-38)))
                          / (Math.Pow(_state.Point.Value.Frequency, ModelParameters.Ef)
                          * Param.BSIM4leff
                          * Param.BSIM4leff
                          * ModelTemperature.Coxe));
                    break;
                case 1:
                    Vds = this._vds;
                    if (Vds < 0.0)
                        Vds = -Vds;

                    Ssi = Eval1ovFNoise(Vds, _state.Point.Value.Frequency, _temperature.Temperature);
                    T10 = ModelParameters.OxideTrapDensityA
                        * Constants.Boltzmann * _temperature.Temperature;
                    T11 = Param.BSIM4weff * Parameters.Nf * Param.BSIM4leff
                        * Math.Pow(_state.Point.Value.Frequency, ModelParameters.Ef) * 1.0e10
                        * this._nstar * this._nstar;
                    Swi = T10 / T11 * this._cd
                        * this._cd;
                    T1 = Swi + Ssi;
                    if (T1 > 0.0)
                        _fl.Compute(m * (Ssi * Swi) / T1);
                    else
                        _fl.Compute(0.0);
                    break;
            }

            if (this._mode >= 0)
            {
                /* bugfix  */
                _nigs.Compute(m * (this._igs + this._igcs));
                _nigd.Compute(m * (this._igd + this._igcd));
            }
            else
            {
                _nigs.Compute(m * (this._igs + this._igcd));
                _nigd.Compute(m * (this._igd + this._igcs));
            }
            _nigb.Compute(m * this._igb);
        }

        private double Eval1ovFNoise(double Vds, double freq, double temp)
        {
            double cd, esat, DelClm, EffFreq, N0, Nl, Leff, Leffsq;
            double T0 = 0.0, T1, T2, T3, T4, T5, T6, T7, T8, T9, Ssi;

            cd = Math.Abs(this._cd);
            Leff = Param.BSIM4leff - 2.0 * ModelParameters.Lintnoi;
            Leffsq = Leff * Leff;
            esat = 2.0 * this._vsattemp / this._ueff;
            if (ModelParameters.Em.Value <= 0.0)
                DelClm = 0.0; /* flicker noise modified -JX  */
            else
            {
                T0 = ((((Vds - this._vdseff) / Param.BSIM4litl)
                           + ModelParameters.Em) / esat);
                DelClm = Param.BSIM4litl * Math.Log(Math.Max(T0, 1e-38));
                if (DelClm < 0.0) DelClm = 0.0;  /* bugfix */
            }
            EffFreq = Math.Pow(freq, ModelParameters.Ef);
            T1 = Constants.Charge * Constants.Charge * Constants.Boltzmann * cd * temp * this._ueff;
            T2 = 1.0e10 * EffFreq * this._abulk * ModelTemperature.Coxe * Leffsq;
            N0 = ModelTemperature.Coxe * this._vgsteff / Constants.Charge;
            Nl = ModelTemperature.Coxe * this._vgsteff
               * (1.0 - this._abovVgst2Vtm * this._vdseff) / Constants.Charge;

            T3 = ModelParameters.OxideTrapDensityA
               * Math.Log(Math.Max(((N0 + this._nstar) / (Nl + this._nstar)), 1e-38));
            T4 = ModelParameters.OxideTrapDensityB * (N0 - Nl);
            T5 = ModelParameters.OxideTrapDensityC * 0.5 * (N0 * N0 - Nl * Nl);

            T6 = Constants.Boltzmann * temp * cd * cd;
            T7 = 1.0e10 * EffFreq * Leffsq * Param.BSIM4weff * Parameters.Nf;
            T8 = ModelParameters.OxideTrapDensityA + ModelParameters.OxideTrapDensityB * Nl
               + ModelParameters.OxideTrapDensityC * Nl * Nl;
            T9 = (Nl + this._nstar) * (Nl + this._nstar);
            Ssi = T1 / T2 * (T3 + T4 + T5) + T6 / T7 * DelClm * T8 / T9;
            return Ssi;
        }
    }
}
