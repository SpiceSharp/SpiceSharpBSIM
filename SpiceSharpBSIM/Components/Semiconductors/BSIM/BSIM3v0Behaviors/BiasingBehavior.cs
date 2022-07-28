using SpiceSharp;
using SpiceSharp.Algebra;
using SpiceSharp.Attributes;
using SpiceSharp.Behaviors;
using SpiceSharp.Components;
using SpiceSharp.Components.Mosfets;
using SpiceSharp.Components.Semiconductors;
using SpiceSharp.Simulations;
using System;

namespace SpiceSharpBSIM.Components.Semiconductors.BSIM.BSIM3v0Behaviors
{
    /// <summary>
    /// Biasing behavior for a <see cref="BSIM3v0"/>.
    /// </summary>
    [BehaviorFor(typeof(BSIM3v0)), AddBehaviorIfNo(typeof(IBiasingBehavior))]
    public partial class BiasingBehavior : TemperatureBehavior, IBiasingBehavior, ITimeBehavior
    {
        private readonly ITemperatureSimulationState _temperature;
        private readonly IIntegrationMethod _method;
        private readonly IIterationSimulationState _iteration;
        private readonly ITimeSimulationState _time;
        private readonly IDerivative _qcheq, _qcdump, _qb, _qg, _qd;

        public const double DELTA_1 = 0.02;
        public const double DELTA_2 = 0.02;
        public const double DELTA_3 = 0.02;
        public const double DELTA_4 = 0.02;

        [ParameterName("gmbs"), ParameterInfo("Gmb")]
        public double Gmbs => _gmbs * Parameters.M;
        [ParameterName("gm"), ParameterInfo("Gm")]
        public double Gm => _gm * Parameters.M;
        [ParameterName("gds"), ParameterInfo("Gds")]
        public double Gds => _gds * Parameters.M;
        [ParameterName("vdsat"), ParameterInfo("Vdsat")]
        public double Vdsat => _vdsat;
        [ParameterName("vth"), ParameterInfo("Vth")]
        public double Von => _von;
        [ParameterName("id"), ParameterInfo("Ids")]
        public double Id => _cd * Parameters.M;
        [ParameterName("vbs"), ParameterInfo("Vbs")]
        public double Vbs => _vbs;
        [ParameterName("vgs"), ParameterInfo("Vgs")]
        public double Vgs => _vgs;
        [ParameterName("vds"), ParameterInfo("Vds")]
        public double Vds => _vds;
        [ParameterName("ibd"), ParameterInfo("Ibd")]
        public double Ibd => _cbd * Parameters.M;
        [ParameterName("ibs"), ParameterInfo("Ibs")]
        public double Ibs => _cbs * Parameters.M;
        [ParameterName("gbd"), ParameterInfo("gbd")]
        public double Gbd => _gbd * Parameters.M;
        [ParameterName("gbs"), ParameterInfo("gbs")]
        public double Gbs => _gbs * Parameters.M;
        [ParameterName("qb"), ParameterInfo("Qbulk")]
        public double Qb => _qb.Value * Parameters.M;
        [ParameterName("cqb"), ParameterInfo("CQbulk")]
        public double CQb => _qb.Derivative * Parameters.M;
        [ParameterName("qg"), ParameterInfo("Qgate")]
        public double Qg => _qg.Value * Parameters.M;
        [ParameterName("cqg"), ParameterInfo("CQgate")]
        public double CQg => _qg.Derivative * Parameters.M;
        [ParameterName("qd"), ParameterInfo("Qdrain")]
        public double Qd => _qd.Value * Parameters.M;
        [ParameterName("cqd"), ParameterInfo("CQdrain")]
        public double CQd => _qd.Derivative * Parameters.M;
        [ParameterName("cgg"), ParameterInfo("Cggb")]
        public double Cgg => _cggb * Parameters.M;
        [ParameterName("cgd"), ParameterInfo("Cgdb")]
        public double Cgd => _cgdb * Parameters.M;
        [ParameterName("cgs"), ParameterInfo("Cgsb")]
        public double Cgs => _cgsb * Parameters.M;
        [ParameterName("cdg"), ParameterInfo("Cdgb")]
        public double Cdg => _cdgb * Parameters.M;
        [ParameterName("cdd"), ParameterInfo("Cddb")]
        public double Cdd => _cddb * Parameters.M;
        [ParameterName("cds"), ParameterInfo("Cdsb")]
        public double Cds => _cdsb * Parameters.M;
        [ParameterName("cbg"), ParameterInfo("Cbgb")]
        public double Cbg => _cbgb * Parameters.M;
        [ParameterName("cbd"), ParameterInfo("Cbdb")]
        public double Cbdb => _cbdb * Parameters.M;
        [ParameterName("cbs"), ParameterInfo("Cbsb")]
        public double Cbs => _cbsb * Parameters.M;
        [ParameterName("capbd"), ParameterInfo("Capbd")]
        public double Capbd => _capbd * Parameters.M;
        [ParameterName("capbs"), ParameterInfo("Capbs")]
        public double Capbs => _capbs * Parameters.M;

        protected bool InitializeSmallSignal { get; set; }
        protected bool InitializeTransient { get; set; }

        private readonly IVariable<double> _drain, _gate, _source, _bulk, _drainPrime, _sourcePrime, _q;
        protected double _vbs, _vgs, _vds, _vbd, _mode, _cd, _gbd, _gmbs, _gm, _gds, _cbs, _cbd,
            _gbs, _von, _thetavth, _ueff, _vdsat, _gbbs, _gbgs, _gbds, _csub, _qinv, _cggb, _cgsb, _cgdb, _cdgb,
            _cdsb, _cddb, _cbgb, _cbsb, _cbdb, _cqdb, _cqsb, _cqgb, _cqbb, _gtau, _tconst, _gtb, _gtd, _gts,
            _capbd, _capbs, _gtg, _qbs, _qbd;
        private readonly Element<double> _dpPtr, _gPtr, _spPtr, _bPtr, _qPtr;
        private readonly Element<double> _ddPtr, _ggPtr, _ssPtr, _bbPtr, _dpdpPtr, _spspPtr, _ddpPtr, _gbPtr, _gdpPtr,
            _gspPtr, _sspPtr, _bgPtr, _bdpPtr, _bspPtr, _dpdPtr, _dpgPtr, _dpbPtr, _dpspPtr, _spgPtr, _spsPtr,
            _spbPtr, _spdpPtr, _qqPtr, _dpqPtr, _spqPtr, _gqPtr, _qgPtr, _qdpPtr, _qspPtr, _qbPtr;

        /// <summary>
        /// Creates a new <see cref="BiasingBehavior"/>.
        /// </summary>
        /// <param name="context">The context.</param>
        public BiasingBehavior(ComponentBindingContext context)
            : base(context)
        {
            _temperature = context.GetState<ITemperatureSimulationState>();
            _iteration = context.GetState<IIterationSimulationState>();
            context.TryGetState(out _time);
            context.TryGetState(out _method);
            var state = context.GetState<IBiasingSimulationState>();

            _drain = state.GetSharedVariable(context.Nodes[0]);
            _gate = state.GetSharedVariable(context.Nodes[1]);
            _source = state.GetSharedVariable(context.Nodes[2]);
            _bulk = state.GetSharedVariable(context.Nodes[3]);
            if (!_drainConductance.Equals(0.0))
                _drainPrime = state.CreatePrivateVariable(Name.Combine("drain"), Units.Volt);
            else
                _drainPrime = _drain;
            if (!_sourceConductance.Equals(0.0))
                _sourcePrime = state.CreatePrivateVariable(Name.Combine("source"), Units.Volt);
            else
                _sourcePrime = _source;

            if (Parameters.NqsMod.Value != 0)
                _q = state.CreatePrivateVariable(Name.Combine("charge"), Units.Coulomb);
            else
                _q = state.GetSharedVariable(Constants.Ground);

            int drain = state.Map[_drain];
            int drainPrime = state.Map[_drainPrime];
            int gate = state.Map[_gate];
            int source = state.Map[_source];
            int sourcePrime = state.Map[_sourcePrime];
            int bulk = state.Map[_bulk];
            int q = state.Map[_q];

            _dpPtr = state.Solver.GetElement(drainPrime);
            _gPtr = state.Solver.GetElement(gate);
            _spPtr = state.Solver.GetElement(sourcePrime);
            _bPtr = state.Solver.GetElement(bulk);
            _qPtr = state.Solver.GetElement(q);
            _ddPtr = state.Solver.GetElement(new MatrixLocation(drain, drain));
            _ggPtr = state.Solver.GetElement(new MatrixLocation(gate, gate));
            _ssPtr = state.Solver.GetElement(new MatrixLocation(source, source));
            _bbPtr = state.Solver.GetElement(new MatrixLocation(bulk, bulk));
            _dpdpPtr = state.Solver.GetElement(new MatrixLocation(drainPrime, drainPrime));
            _spspPtr = state.Solver.GetElement(new MatrixLocation(sourcePrime, sourcePrime));
            _ddpPtr = state.Solver.GetElement(new MatrixLocation(drain, drainPrime));
            _gbPtr = state.Solver.GetElement(new MatrixLocation(gate, bulk));
            _gdpPtr = state.Solver.GetElement(new MatrixLocation(gate, drainPrime));
            _gspPtr = state.Solver.GetElement(new MatrixLocation(gate, sourcePrime));
            _sspPtr = state.Solver.GetElement(new MatrixLocation(source, sourcePrime));
            _bdpPtr = state.Solver.GetElement(new MatrixLocation(bulk, drainPrime));
            _bspPtr = state.Solver.GetElement(new MatrixLocation(bulk, sourcePrime));
            _dpspPtr = state.Solver.GetElement(new MatrixLocation(drainPrime, sourcePrime));
            _dpdPtr = state.Solver.GetElement(new MatrixLocation(drainPrime, drain));
            _bgPtr = state.Solver.GetElement(new MatrixLocation(bulk, gate));
            _dpgPtr = state.Solver.GetElement(new MatrixLocation(drainPrime, gate));
            _spgPtr = state.Solver.GetElement(new MatrixLocation(sourcePrime, gate));
            _spsPtr = state.Solver.GetElement(new MatrixLocation(sourcePrime, source));
            _dpbPtr = state.Solver.GetElement(new MatrixLocation(drainPrime, bulk));
            _spbPtr = state.Solver.GetElement(new MatrixLocation(sourcePrime, bulk));
            _spdpPtr = state.Solver.GetElement(new MatrixLocation(sourcePrime, drainPrime));

            _qqPtr = state.Solver.GetElement(new MatrixLocation(q, q));

            _qdpPtr = state.Solver.GetElement(new MatrixLocation(q, drainPrime));
            _qspPtr = state.Solver.GetElement(new MatrixLocation(q, sourcePrime));
            _qgPtr = state.Solver.GetElement(new MatrixLocation(q, gate));
            _qbPtr = state.Solver.GetElement(new MatrixLocation(q, bulk));
            _dpqPtr = state.Solver.GetElement(new MatrixLocation(drainPrime, q));
            _spqPtr = state.Solver.GetElement(new MatrixLocation(sourcePrime, q));
            _gqPtr = state.Solver.GetElement(new MatrixLocation(gate, q));
            // Element is never used...
            // _bqPtr = state.Solver.GetElement(new MatrixLocation(bulk, q));

            _qg = _method?.CreateDerivative() ?? new SimpleDerivative();
            _qd = _method?.CreateDerivative() ?? new SimpleDerivative();
            _qb = _method?.CreateDerivative() ?? new SimpleDerivative();
            if (Parameters.NqsMod.Value != 0)
            {
                _qcdump = _method?.CreateDerivative() ?? new SimpleDerivative();
                _qcheq = _method?.CreateDerivative() ?? new SimpleDerivative();
            }
            else
            {
                _qcdump = new SimpleDerivative();
                _qcheq = new SimpleDerivative();
            }
        }

        /// <inheritdoc />
        void ITimeBehavior.InitializeStates()
        {
            InitializeTransient = true;
            ((IBiasingBehavior)this).Load();
            InitializeTransient = false;
        }

        /// <inheritdoc />
        void IBiasingBehavior.Load()
        {
            double SourceSatCurrent, DrainSatCurrent;
            double ag0, qgd, qgs, qgb, von, VgstNVt, ExpVgst = 0.0;
            double cdrain, cdreq, ceqbd, ceqbs, ceqqb, ceqqd, ceqqg;
            double czbd, czbdsw, czbs, czbssw, evbd, evbs, arg, sarg;
            double Vfbeff, dVfbeff_dVg, dVfbeff_dVd, dVfbeff_dVb, V3, V4;
            double gcbdb, gcbgb, gcbsb, gcddb, gcdgb, gcdsb, gcgdb, gcggb, gcgsb, gcsdb;

            double gcsgb, gcssb, PhiB, PhiBSW, MJ, MJSW;
            double vbd, vbs, vds, vgb, vgd, vgs, vgdo;

            double qgate, qbulk, qdrn, qsrc, cqgate, cqbulk, cqdrn;
            double Vds, Vgs, Vbs, Gmbs, FwdSum, RevSum;
            double Vgs_eff, Vfb, dVfb_dVb, dVfb_dVd;
            double Phis, dPhis_dVb, sqrtPhis, dsqrtPhis_dVb, Vth, dVth_dVb, dVth_dVd;
            double Vgst, dVgs_eff_dVg;
            double n, dn_dVb, Vtm;
            double ExpArg;
            double Denomi, dDenomi_dVg, dDenomi_dVd, dDenomi_dVb;
            double ueff, dueff_dVg, dueff_dVd, dueff_dVb;
            double Esat, Vdsat;
            double EsatL, dEsatL_dVg, dEsatL_dVd, dEsatL_dVb;
            double dVdsat_dVg, dVdsat_dVb, dVdsat_dVd, Vasat;
            double dVasat_dVg, dVasat_dVb, dVasat_dVd, Va, dVa_dVd, dVa_dVg, dVa_dVb;
            double Vbseff, dVbseff_dVb, VbseffCV, dVbseffCV_dVb;
            double Arg1, One_Third_CoxWL, Two_Third_CoxWL, CoxWL;
            double T0, dT0_dVg, dT0_dVd, dT0_dVb;
            double T1, dT1_dVg, dT1_dVd, dT1_dVb;
            double T2, dT2_dVg, dT2_dVd, dT2_dVb;
            double T3, dT3_dVg, dT3_dVd, dT3_dVb;
            double T4, dT4_dVg, dT4_dVd, dT4_dVb;
            double T5;
            double T6;
            double T7;
            double T8;
            double T10;
            double Abulk, dAbulk_dVb, Abulk0, dAbulk0_dVb;
            double VACLM, dVACLM_dVg, dVACLM_dVd, dVACLM_dVb;
            double VADIBL, dVADIBL_dVg, dVADIBL_dVd, dVADIBL_dVb;
            double Xdep, dXdep_dVb, lt1, dlt1_dVb, ltw, dltw_dVb, Delt_vth, dDelt_vth_dVb;
            double Theta0, dTheta0_dVb;
            double TempRatio, tmp1, tmp2, tmp3, tmp4;
            double DIBL_Sft, dDIBL_Sft_dVd, Pmos_factor;

            double a1;

            double Vgsteff, dVgsteff_dVg, dVgsteff_dVd, dVgsteff_dVb;
            double Vdseff, dVdseff_dVg, dVdseff_dVd, dVdseff_dVb;
            double VdseffCV, dVdseffCV_dVg, dVdseffCV_dVd, dVdseffCV_dVb;
            double diffVds;
            double dAbulk_dVg, dn_dVd;
            double beta, dbeta_dVg, dbeta_dVd, dbeta_dVb;
            double gche, dgche_dVg, dgche_dVd, dgche_dVb;
            double fgche1, dfgche1_dVg, dfgche1_dVd, dfgche1_dVb;
            double fgche2, dfgche2_dVg, dfgche2_dVd, dfgche2_dVb;
            double Idl, dIdl_dVg, dIdl_dVd, dIdl_dVb;
            double Idsa, dIdsa_dVg, dIdsa_dVd, dIdsa_dVb;
            double Ids, Gm, Gds, Gmb;
            double Isub, Gbd, Gbg, Gbb;
            double VASCBE, dVASCBE_dVg, dVASCBE_dVd, dVASCBE_dVb;
            double CoxWovL;
            double Rds, dRds_dVg, dRds_dVb, WVCox, WVCoxRds;
            double Vgst2Vtm, VdsatCV, dVdsatCV_dVg, dVdsatCV_dVb;
            double Leff, Weff, dWeff_dVg, dWeff_dVb;
            double AbulkCV, dAbulkCV_dVb;
            double qgdo, qgso, cgdo, cgso;

            double qcheq, qdef, gqdef = 0, cqdef = 0, cqcheq, gtau_diff, gtau_drift;
            double gcqdb, gcqsb, gcqgb, gcqbb;
            double dxpart, sxpart;

            double gbspsp, gbbdp, gbbsp, gbspg, gbspb, gbspdp;
            double gbdpdp, gbdpg, gbdpb, gbdpsp;
            double Cgg, Cgd, Cgb;
            double Csg, Csd, Csb, Cbg, Cbd, Cbb;
            double Cgg1, Cgb1, Cgd1, Cbg1, Cbb1, Cbd1, Qac0, Qsub0;
            double dQac0_dVg, dQac0_dVd, dQac0_dVb, dQsub0_dVg, dQsub0_dVd, dQsub0_dVb;


            bool Check, ChargeComputationNeeded;
            double m;

            Check = true;
            if (InitializeSmallSignal || InitializeTransient)
            {
                vbs = this._vbs;
                vgs = this._vgs;
                vds = this._vds;
                qdef = this._qcdump.Value;
            }
            else if (_iteration.Mode == IterationModes.Junction && !Parameters.Off)
            {
                vds = ModelParameters.Type * Parameters.IcVDS;
                vgs = ModelParameters.Type * Parameters.IcVGS;
                vbs = ModelParameters.Type * Parameters.IcVBS;
                qdef = 0.0;

                if (vds.Equals(0.0) && vgs.Equals(0.0) && vbs.Equals(0.0))
                {
                    vbs = 0.0;
                    vgs = Param.BSIM3v0vth0 + 0.1;
                    vds = 0.1;
                }
            }
            else if ((_iteration.Mode == IterationModes.Junction || _iteration.Mode == IterationModes.Fix) && Parameters.Off)
            {
                qdef = vbs = vgs = vds = 0.0;
            }
            else
            {
                /* PREDICTOR */
                vbs = ModelParameters.Type
                    * (this._bulk.Value
                    - this._sourcePrime.Value);
                vgs = ModelParameters.Type
                    * (this._gate.Value
                    - this._sourcePrime.Value);
                vds = ModelParameters.Type
                    * (this._drainPrime.Value
                    - this._sourcePrime.Value);
                qdef = this._q.Value;
                /* PREDICTOR */

                vbd = vbs - vds;
                vgd = vgs - vds;
                vgdo = this._vgs
                     - this._vds;
                
                /*NOBYPASS*/
                Check = false;
                von = this._von;
                if (this._vds >= 0.0)
                {
                    vgs = Transistor.LimitFet(vgs, this._vgs, von);
                    vds = vgs - vgd;
                    vds = Transistor.LimitVds(vds, this._vds);
                    vgd = vgs - vds;

                }
                else
                {
                    vgd = Transistor.LimitFet(vgd, vgdo, von);
                    vds = vgs - vgd;
                    vds = -Transistor.LimitVds(-vds, -(this._vds));
                    vgs = vgd + vds;
                }

                if (vds >= 0.0)
                {
                    vbs = Semiconductor.LimitJunction(vbs, this._vbs,
                                        Constants.Vt0, ModelTemperature.Vcrit, ref Check);
                    vbd = vbs - vds;

                }
                else
                {
                    vbd = Semiconductor.LimitJunction(vbd, this._vbd,
                                        Constants.Vt0, ModelTemperature.Vcrit, ref Check);
                    vbs = vbd + vds;
                }
            }

            /* determine DC current and derivatives */
            vbd = vbs - vds;
            vgd = vgs - vds;
            vgb = vgs - vbs;

            m = Parameters.M;

            SourceSatCurrent = Parameters.SourceArea * ModelParameters.JctSatCurDensity;
            if (Parameters.SourceArea <= 0.0)
                SourceSatCurrent = 1.0e-14;
            if (SourceSatCurrent <= 0.0)
            {
                this._gbs = _iteration.Gmin;
                this._cbs = this._gbs * vbs;
            }
            else if (vbs <= -0.5)
            {
                this._gbs = _iteration.Gmin;
                this._cbs = -SourceSatCurrent + this._gbs * vbs;
            }
            else if (vbs < 0.5)
            {
                evbs = Math.Exp(vbs / Constants.Vt0);
                this._gbs = SourceSatCurrent * evbs / Constants.Vt0
                               + _iteration.Gmin;
                this._cbs = SourceSatCurrent * (evbs - 1.0)
                               + _iteration.Gmin * vbs;
            }
            else
            {
                evbs = Math.Exp(0.5 / Constants.Vt0);
                this._gbs = SourceSatCurrent * evbs / Constants.Vt0
                               + _iteration.Gmin;
                this._cbs = SourceSatCurrent * (evbs - 1.0)
                               + _iteration.Gmin * 0.5;
            }

            DrainSatCurrent = Parameters.DrainArea * ModelParameters.JctSatCurDensity;
            if (Parameters.DrainArea <= 0.0)
                DrainSatCurrent = 1.0e-14;
            if (DrainSatCurrent <= 0.0)
            {
                this._gbd = _iteration.Gmin;
                this._cbd = this._gbd * vbd;
            }
            else if (vbd <= -0.5)
            {
                this._gbd = _iteration.Gmin;
                this._cbd = -DrainSatCurrent + this._gbd * vbd;
            }
            else if (vbd < 0.5)
            {
                evbd = Math.Exp(vbd / Constants.Vt0);
                this._gbd = DrainSatCurrent * evbd / Constants.Vt0
                               + _iteration.Gmin;
                this._cbd = DrainSatCurrent * (evbd - 1.0)
                               + _iteration.Gmin * vbd;
            }
            else
            {
                evbd = Math.Exp(0.5 / Constants.Vt0);
                this._gbd = DrainSatCurrent * evbd / Constants.Vt0
                               + _iteration.Gmin;
                this._cbd = DrainSatCurrent * (evbd - 1.0)
                               + _iteration.Gmin * 0.5;
            }

            if (vds >= 0.0)
            {   /* normal mode */
                this._mode = 1;
                Vds = vds;
                Vgs = vgs;
                Vbs = vbs;
            }
            else
            {   /* inverse mode */
                this._mode = -1;
                Vds = -vds;
                Vgs = vgd;
                Vbs = vbd;
            }

            ChargeComputationNeeded = _time != null || InitializeSmallSignal || InitializeTransient;

            T0 = Vbs - Param.BSIM3v0vbsc - 0.001;
            T1 = Math.Sqrt(T0 * T0 - 0.004 * Param.BSIM3v0vbsc);
            Vbseff = Param.BSIM3v0vbsc + 0.5 * (T0 + T1);
            dVbseff_dVb = 0.5 * (1.0 + T0 / T1);

            if (Vbseff > 0.0)
            {
                T0 = Param.BSIM3v0phi / (Param.BSIM3v0phi + Vbseff);
                Phis = Param.BSIM3v0phi * T0;
                dPhis_dVb = -T0 * T0;
                sqrtPhis = Param.BSIM3v0phis3 / (Param.BSIM3v0phi + 0.5 * Vbseff);
                dsqrtPhis_dVb = -0.5 * sqrtPhis * sqrtPhis / Param.BSIM3v0phis3;
            }
            else
            {
                Phis = Param.BSIM3v0phi - Vbseff;
                dPhis_dVb = -1.0;
                sqrtPhis = Math.Sqrt(Phis);
                dsqrtPhis_dVb = -0.5 / sqrtPhis;
            }
            Xdep = Param.BSIM3v0Xdep0 * sqrtPhis / Param.BSIM3v0sqrtPhi;
            dXdep_dVb = (Param.BSIM3v0Xdep0 / Param.BSIM3v0sqrtPhi)
              * dsqrtPhis_dVb;

            Leff = Param.BSIM3v0leff;
            /* Vth Calculation */
            if ((T1 = 1.0 + Param.BSIM3v0dvt2 * Vbseff) < 1.0e-10)
                T1 = 1.0E-10;
            if ((T2 = 1.0 + Param.BSIM3v0dvt2w * Vbseff) < 1.0E-10)
                T2 = 1.0E-10;

            T3 = Math.Sqrt(Xdep);
            lt1 = ModelTemperature.Factor1 * T3 * T1;
            dlt1_dVb = ModelTemperature.Factor1 * (0.5 / T3 * T1 * dXdep_dVb
                     + T3 * Param.BSIM3v0dvt2);

            ltw = ModelTemperature.Factor1 * T3 * T2;
            dltw_dVb = ModelTemperature.Factor1 * (0.5 / T3 * T2 * dXdep_dVb
                     + T3 * Param.BSIM3v0dvt2w);

            T0 = -0.5 * Param.BSIM3v0dvt1 * Leff / lt1;
            if (T0 > -EXP_THRESHOLD)
            {
                T1 = Math.Exp(T0);
                dT1_dVb = -T0 / lt1 * T1 * dlt1_dVb;
            }
            else
            {
                T1 = MIN_EXP;
                dT1_dVb = 0.0;
            }

            Theta0 = T1 * (1.0 + 2.0 * T1);
            dTheta0_dVb = (1.0 + 4.0 * T1) * dT1_dVb;

            this._thetavth = Param.BSIM3v0dvt0 * Theta0;
            T0 = Param.BSIM3v0vbi - Param.BSIM3v0phi;
            Delt_vth = this._thetavth * T0;
            dDelt_vth_dVb = Param.BSIM3v0dvt0 * dTheta0_dVb * T0;

            T0 = -0.5 * Param.BSIM3v0dvt1w * Param.BSIM3v0weff * Leff / ltw;
            if (T0 > -EXP_THRESHOLD)
            {
                T1 = Math.Exp(T0);
                dT1_dVb = -T0 / ltw * T1 * dltw_dVb;
            }
            else
            {
                T1 = MIN_EXP;
                dT1_dVb = 0.0;
            }

            T2 = T1 * (1.0 + 2.0 * T1);
            dT2_dVb = (1.0 + 4.0 * T1) * dT1_dVb;

            T0 = Param.BSIM3v0dvt0w * T2;
            T1 = Param.BSIM3v0vbi - Param.BSIM3v0phi;
            T2 = T0 * T1;
            dT2_dVb = Param.BSIM3v0dvt0w * dT2_dVb * T1;

            TempRatio = _temperature.Temperature / ModelParameters.Tnom - 1.0;
            T0 = Math.Sqrt(1.0 + Param.BSIM3v0nlx / Leff);
            T1 = Param.BSIM3v0k1 * (T0 - 1.0) * Param.BSIM3v0sqrtPhi
               + (Param.BSIM3v0kt1 + Param.BSIM3v0kt1l / Leff
               + Param.BSIM3v0kt2 * Vbseff) * TempRatio;
            tmp2 = ModelParameters.Tox / (Param.BSIM3v0weff
             + Param.BSIM3v0w0) * Param.BSIM3v0phi;

            dDIBL_Sft_dVd = (Param.BSIM3v0eta0 + Param.BSIM3v0etab
                          * Vbseff) * Param.BSIM3v0theta0vb0;
            DIBL_Sft = dDIBL_Sft_dVd * Vds;

            Vth = ModelParameters.Type * Param.BSIM3v0vth0 + Param.BSIM3v0k1
                * (sqrtPhis - Param.BSIM3v0sqrtPhi) - Param.BSIM3v0k2
                * Vbseff - Delt_vth - T2 + (Param.BSIM3v0k3 + Param.BSIM3v0k3b
                * Vbseff) * tmp2 + T1 - DIBL_Sft;

            this._von = Vth;

            dVth_dVb = Param.BSIM3v0k1 * dsqrtPhis_dVb - Param.BSIM3v0k2
                     - dDelt_vth_dVb - dT2_dVb + Param.BSIM3v0k3b * tmp2
                     - Param.BSIM3v0etab * Vds * Param.BSIM3v0theta0vb0
                     + Param.BSIM3v0kt2 * TempRatio;
            dVth_dVd = -dDIBL_Sft_dVd;

            /* Calculate n */
            tmp2 = Param.BSIM3v0nfactor * EPSSI / Xdep;
            tmp3 = Param.BSIM3v0cdsc + Param.BSIM3v0cdscb * Vbseff
                 + Param.BSIM3v0cdscd * Vds;
            n = 1.0 + (tmp2 + tmp3 * Theta0 + Param.BSIM3v0cit) / ModelTemperature.Cox;
            if (n > 1.0)
            {
                dn_dVb = (-tmp2 / Xdep * dXdep_dVb + tmp3 * dTheta0_dVb
                           + Param.BSIM3v0cdscb * Theta0) / ModelTemperature.Cox;
                dn_dVd = Param.BSIM3v0cdscd * Theta0 / ModelTemperature.Cox;
            }
            else
            {
                n = 1.0;
                dn_dVb = dn_dVd = 0.0;
            }

            /* Poly Gate Si Depletion Effect */
            T0 = Param.BSIM3v0vfb + Param.BSIM3v0phi;
            if (ModelParameters.Ngate.Given && (Vgs > T0))
            {
                T1 = 1.0e6 * Charge_q * EPSSI * Param.BSIM3v0ngate
                       / (ModelTemperature.Cox * ModelTemperature.Cox);
                T4 = Math.Sqrt(1.0 + 2.0 * (Vgs - T0) / T1);
                T2 = T1 * (T4 - 1.0);
                T3 = 0.5 * T2 * T2 / T1;

                if (T3 < 1.12)
                {
                    Vgs_eff = T0 + T2;
                    dVgs_eff_dVg = 1.0 / T4;
                }
                else
                {
                    Vgs_eff = Vgs - 1.12;
                    dVgs_eff_dVg = 1.0;
                }
            }
            else
            {
                Vgs_eff = Vgs;
                dVgs_eff_dVg = 1.0;
            }
            Vgst = Vgs_eff - Vth;

            /* Effective Vgst (Vgsteff) Calculation */
            Vtm = ModelTemperature.Vtm;

            T10 = 2.0 * n * Vtm;
            VgstNVt = Vgst / T10;
            if (VgstNVt < -EXP_THRESHOLD)
            {
                T1 = T10 * MIN_EXP;
                dT1_dVg = dT1_dVd = dT1_dVb = 0.0;
            }
            else if (VgstNVt > EXP_THRESHOLD)
            {
                T1 = Vgst;
                dT1_dVg = dVgs_eff_dVg;
                dT1_dVd = -dVth_dVd;
                dT1_dVb = -dVth_dVb;
            }
            else
            {
                ExpVgst = Math.Exp(VgstNVt);
                T1 = T10 * Math.Log(1.0 + ExpVgst);
                dT1_dVg = ExpVgst / (1.0 + ExpVgst);
                dT1_dVb = -dT1_dVg * (dVth_dVb + Vgst / n * dn_dVb)
                + T1 / n * dn_dVb;
                dT1_dVd = -dT1_dVg * (dVth_dVd + Vgst / n * dn_dVd)
                    + T1 / n * dn_dVd;
                dT1_dVg *= dVgs_eff_dVg;
            }

            T2 = ModelParameters.Tox / (Param.BSIM3v0weff + Param.BSIM3v0w0);
            ExpArg = (2.0 * Param.BSIM3v0voff - Vgst) / T10;
            if (ExpArg < -EXP_THRESHOLD)
            {
                T2 = 1.0;
                dT2_dVg = dT2_dVd = dT2_dVb = 0.0;

            }
            else if (ExpArg > EXP_THRESHOLD)
            {
                T2 = 1.0 + 2.0 * n * ModelTemperature.Cox / Param.BSIM3v0cdep0
               * MAX_EXP;
                dT2_dVg = dT2_dVd = dT2_dVb = 0.0;
            }
            else
            {
                dT2_dVg = -ModelTemperature.Cox / Vtm / Param.BSIM3v0cdep0
                    * Math.Exp(ExpArg);
                T2 = 1.0 - T10 * dT2_dVg;
                dT2_dVd = -dT2_dVg * (dVth_dVd - 2.0 * Vtm
                * ExpArg * dn_dVd) + (T2 - 1.0) / n * dn_dVd;
                dT2_dVb = -dT2_dVg * (dVth_dVb - 2.0 * Vtm
                * ExpArg * dn_dVb) + (T2 - 1.0) / n * dn_dVb;
                dT2_dVg *= dVgs_eff_dVg;
            }

            Vgsteff = T1 / T2;
            dVgsteff_dVg = (T2 * dT1_dVg - T1 * dT2_dVg) / (T2 * T2);
            dVgsteff_dVd = (T2 * dT1_dVd - T1 * dT2_dVd) / (T2 * T2);
            dVgsteff_dVb = (T2 * dT1_dVb - T1 * dT2_dVb) / (T2 * T2);

            /* Calculate Effective Channel Geometry */
            Weff = Param.BSIM3v0weff - 2.0 * (Param.BSIM3v0dwg * Vgsteff
                 + Param.BSIM3v0dwb * (sqrtPhis - Param.BSIM3v0sqrtPhi));
            dWeff_dVg = -2.0 * Param.BSIM3v0dwg;
            dWeff_dVb = -2.0 * Param.BSIM3v0dwb * dsqrtPhis_dVb;

            if (Weff < 1.0e-8)
            {
                Weff = 1.0e-8;
                dWeff_dVg = dWeff_dVb = 0;
            }

            Rds = Param.BSIM3v0rds0 * (1.0 + Param.BSIM3v0prwg * Vgsteff
                + Param.BSIM3v0prwb * (sqrtPhis - Param.BSIM3v0sqrtPhi));
            if (Rds > 0.0)
            {
                dRds_dVg = Param.BSIM3v0rds0 * Param.BSIM3v0prwg;
                dRds_dVb = Param.BSIM3v0rds0 * Param.BSIM3v0prwb * dsqrtPhis_dVb;
            }
            else
            {
                Rds = dRds_dVg = dRds_dVb = 0.0;
            }

            WVCox = Weff * Param.BSIM3v0vsattemp * ModelTemperature.Cox;
            WVCoxRds = WVCox * Rds;

            /* Calculate Abulk */
            T0 = 1.0 / (1.0 + Param.BSIM3v0keta * Vbseff);
            dT0_dVb = -Param.BSIM3v0keta * T0 * T0;

            T1 = 0.5 * Param.BSIM3v0k1 / sqrtPhis;
            dT1_dVb = -T1 / sqrtPhis * dsqrtPhis_dVb;

            tmp1 = Leff + 2.0 * Math.Sqrt(Param.BSIM3v0xj * Xdep);
            T5 = Leff / tmp1;
            tmp2 = Param.BSIM3v0a0 * T5;
            tmp3 = Param.BSIM3v0weff + Param.BSIM3v0b1;
            tmp4 = Param.BSIM3v0b0 / tmp3;
            T2 = tmp2 + tmp4;
            dT2_dVb = -tmp2 / tmp1 * Math.Sqrt(Param.BSIM3v0xj / Xdep) * dXdep_dVb;
            T6 = T5 * T5;
            T7 = T5 * T6;
            Abulk0 = T0 * (1.0 + T1 * T2);
            if (Abulk0 < 0.01)
                Abulk0 = 0.01;
            T8 = Param.BSIM3v0ags * Param.BSIM3v0a0 * T7;
            dAbulk_dVg = -T1 * T0 * T8;
            Abulk = Abulk0 + dAbulk_dVg * Vgsteff;
            if (Abulk < 0.01)
                Abulk = 0.01;

            dAbulk0_dVb = T0 * T1 * dT2_dVb + T0 * T2 * dT1_dVb
                           + (1.0 + T1 * T2) * dT0_dVb;
            dAbulk_dVb = dAbulk0_dVb - T8 * Vgsteff * (T1 * (3.0 * T0 * dT2_dVb
               / tmp2 + dT0_dVb) + T0 * dT1_dVb);
            /* Mobility calculation */

            if (ModelParameters.MobMod.Value == 1)
            {
                T0 = Vgsteff + Vth + Vth;
                T2 = Param.BSIM3v0ua + Param.BSIM3v0uc * Vbseff;
                T3 = T0 / ModelParameters.Tox;
                Denomi = 1.0 + T3 * (T2
                       + Param.BSIM3v0ub * T3);
                T1 = T2 / ModelParameters.Tox + 2.0 * Param.BSIM3v0ub * T3
                   / ModelParameters.Tox;
                dDenomi_dVg = T1;
                dDenomi_dVd = T1 * 2.0 * dVth_dVd;
                dDenomi_dVb = T1 * 2.0 * dVth_dVb + Param.BSIM3v0uc * T3;
            }
            else if (ModelParameters.MobMod.Value == 2)
            {
                Denomi = 1.0 + Vgsteff / ModelParameters.Tox * (Param.BSIM3v0ua
                   + Param.BSIM3v0uc * Vbseff + Param.BSIM3v0ub * Vgsteff
                           / ModelParameters.Tox);
                T1 = (Param.BSIM3v0ua + Param.BSIM3v0uc * Vbseff) / ModelParameters.Tox
                   + 2.0 * Param.BSIM3v0ub / (ModelParameters.Tox * ModelParameters.Tox)
           * Vgsteff;
                dDenomi_dVg = T1;
                dDenomi_dVd = 0.0;
                dDenomi_dVb = Vgsteff * Param.BSIM3v0uc / ModelParameters.Tox;
            }
            else
            {
                T0 = Vgsteff + Vth + Vth;
                T2 = 1.0 + Param.BSIM3v0uc * Vbseff;
                T3 = T0 / ModelParameters.Tox;
                T4 = T3 * (Param.BSIM3v0ua + Param.BSIM3v0ub * T3);
                Denomi = 1.0 + T4 * T2;

                T1 = (Param.BSIM3v0ua / ModelParameters.Tox + 2.0 * Param.BSIM3v0ub
           * T3 / ModelParameters.Tox) * T2;
                dDenomi_dVg = T1;
                dDenomi_dVd = T1 * 2.0 * dVth_dVd;
                dDenomi_dVb = T1 * 2.0 * dVth_dVb + Param.BSIM3v0uc * T4;
            }

            this._ueff = ueff = Param.BSIM3v0u0temp / Denomi;
            dueff_dVg = -ueff / Denomi * dDenomi_dVg;
            dueff_dVd = -ueff / Denomi * dDenomi_dVd;
            dueff_dVb = -ueff / Denomi * dDenomi_dVb;


            /* Saturation Drain Voltage  Vdsat */
            Esat = 2.0 * Param.BSIM3v0vsattemp / ueff;
            EsatL = Esat * Leff;
            T0 = -EsatL / ueff;
            dEsatL_dVg = T0 * dueff_dVg;
            dEsatL_dVd = T0 * dueff_dVd;
            dEsatL_dVb = T0 * dueff_dVb;

            a1 = Param.BSIM3v0a1;
            if ((Pmos_factor = a1 * Vgsteff + Param.BSIM3v0a2) > 1.0)
            {
                Pmos_factor = 1.0;
                a1 = 0.0;
            }

            Vgst2Vtm = Vgsteff + 2.0 * Vtm;
            if ((Rds == 0.0) && (Pmos_factor == 1.0))
            {
                T0 = 1.0 / (Abulk * EsatL + Vgst2Vtm);
                this._vdsat = Vdsat = EsatL * Vgst2Vtm * T0;

                dT0_dVg = -(Abulk * dEsatL_dVg + EsatL * dAbulk_dVg + 1.0) * T0 * T0;
                dT0_dVd = -(Abulk * dEsatL_dVd) * T0 * T0;
                dT0_dVb = -(Abulk * dEsatL_dVb + dAbulk_dVb * EsatL) * T0 * T0;

                dVdsat_dVg = EsatL * Vgst2Vtm * dT0_dVg + EsatL * T0
                           + Vgst2Vtm * T0 * dEsatL_dVg;
                dVdsat_dVd = EsatL * Vgst2Vtm * dT0_dVd
               + Vgst2Vtm * T0 * dEsatL_dVd;
                dVdsat_dVb = EsatL * Vgst2Vtm * dT0_dVb
               + Vgst2Vtm * T0 * dEsatL_dVb;
            }
            else
            {
                tmp1 = a1 / (Pmos_factor * Pmos_factor);
                if (Rds > 0)
                {
                    tmp2 = dRds_dVg / Rds + dWeff_dVg / Weff;
                    tmp3 = dRds_dVb / Rds + dWeff_dVb / Weff;
                }
                else
                {
                    tmp2 = dWeff_dVg / Weff;
                    tmp3 = dWeff_dVb / Weff;
                }

                T0 = 2.0 * Abulk * (Abulk * WVCoxRds - 1.0 + 1.0 / Pmos_factor);
                dT0_dVg = 2.0 * (Abulk * Abulk * WVCoxRds * tmp2 - Abulk * tmp1
                + (2.0 * WVCoxRds * Abulk + 1.0 / Pmos_factor - 1.0)
                * dAbulk_dVg);

                dT0_dVb = 2.0 * (Abulk * Abulk * WVCoxRds * (2.0 / Abulk
                * dAbulk_dVb + tmp3) + (1.0 / Pmos_factor - 1.0)
                * dAbulk_dVb);
                dT0_dVd = 0.0;
                T1 = Vgst2Vtm * (2.0 / Pmos_factor - 1.0) + Abulk
                   * EsatL + 3.0 * Abulk * Vgst2Vtm * WVCoxRds;

                dT1_dVg = (2.0 / Pmos_factor - 1.0) - 2.0 * Vgst2Vtm * tmp1
                + Abulk * dEsatL_dVg + EsatL * dAbulk_dVg + 3.0 * (Abulk
                * WVCoxRds + Abulk * Vgst2Vtm * WVCoxRds * tmp2 + Vgst2Vtm
                * WVCoxRds * dAbulk_dVg);
                dT1_dVb = Abulk * dEsatL_dVb + EsatL * dAbulk_dVb
                    + 3.0 * (Vgst2Vtm * WVCoxRds * dAbulk_dVb
                + Abulk * Vgst2Vtm * WVCoxRds * tmp3);
                dT1_dVd = Abulk * dEsatL_dVd;

                T2 = Vgst2Vtm * (EsatL + 2.0 * Vgst2Vtm * WVCoxRds);
                dT2_dVg = EsatL + Vgst2Vtm * dEsatL_dVg + Vgst2Vtm * WVCoxRds
                * (4.0 + 2.0 * Vgst2Vtm * tmp2);
                dT2_dVb = Vgst2Vtm * dEsatL_dVb + 2.0 * Vgst2Vtm * WVCoxRds
                * Vgst2Vtm * tmp3;
                dT2_dVd = Vgst2Vtm * dEsatL_dVd;

                T3 = Math.Sqrt(T1 * T1 - 2.0 * T0 * T2);
                this._vdsat = Vdsat = (T1 - T3) / T0;

                dT3_dVg = (T1 * dT1_dVg - 2.0 * (T0 * dT2_dVg + T2 * dT0_dVg))
                    / T3;
                dT3_dVd = (T1 * dT1_dVd - 2.0 * (T0 * dT2_dVd + T2 * dT0_dVd))
                / T3;
                dT3_dVb = (T1 * dT1_dVb - 2.0 * (T0 * dT2_dVb + T2 * dT0_dVb))
                / T3;

                T4 = T1 - T3;
                dT4_dVg = -dT1_dVg - dT3_dVg;
                dT4_dVd = -dT1_dVd - dT3_dVd;
                dT4_dVb = -dT1_dVb - dT3_dVb;

                dVdsat_dVg = (dT1_dVg - (T1 * dT1_dVg - dT0_dVg * T2
               - T0 * dT2_dVg) / T3 - Vdsat * dT0_dVg) / T0;
                dVdsat_dVb = (dT1_dVb - (T1 * dT1_dVb - dT0_dVb * T2
               - T0 * dT2_dVb) / T3 - Vdsat * dT0_dVb) / T0;
                dVdsat_dVd = (dT1_dVd - (T1 * dT1_dVd - T0 * dT2_dVd) / T3) / T0;
            }

            /* Effective Vds (Vdseff) Calculation */
            T1 = Vdsat - Vds - Param.BSIM3v0delta;
            dT1_dVg = dVdsat_dVg;
            dT1_dVd = dVdsat_dVd - 1.0;
            dT1_dVb = dVdsat_dVb;

            T2 = Math.Sqrt(T1 * T1 + 4.0 * Param.BSIM3v0delta * Vdsat);
            T0 = T1 / T2;
            T3 = 2.0 * Param.BSIM3v0delta / T2;
            dT2_dVg = T0 * dT1_dVg + T3 * dVdsat_dVg;
            dT2_dVd = T0 * dT1_dVd + T3 * dVdsat_dVd;
            dT2_dVb = T0 * dT1_dVb + T3 * dVdsat_dVb;

            Vdseff = Vdsat - 0.5 * (T1 + T2);
            dVdseff_dVg = dVdsat_dVg - 0.5 * (dT1_dVg + dT2_dVg);
            dVdseff_dVd = dVdsat_dVd - 0.5 * (dT1_dVd + dT2_dVd);
            dVdseff_dVb = dVdsat_dVb - 0.5 * (dT1_dVb + dT2_dVb);

            /* Calculate VAsat */
            tmp1 = a1 / (Pmos_factor * Pmos_factor);
            if (Rds > 0)
            {
                tmp2 = dRds_dVg / Rds + dWeff_dVg / Weff;
                tmp3 = dRds_dVb / Rds + dWeff_dVb / Weff;
            }
            else
            {
                tmp2 = dWeff_dVg / Weff;
                tmp3 = dWeff_dVb / Weff;
            }
            tmp4 = 1.0 - 0.5 * Abulk * Vdsat / Vgst2Vtm;
            T0 = EsatL + Vdsat + 2.0 * WVCoxRds * Vgsteff * tmp4;

            dT0_dVg = dEsatL_dVg + dVdsat_dVg + 2.0 * WVCoxRds * tmp4
            * (1.0 + tmp2 * Vgsteff) - WVCoxRds * Vgsteff / Vgst2Vtm
            * (Abulk * dVdsat_dVg - Abulk * Vdsat / Vgst2Vtm
            + Vdsat * dAbulk_dVg);

            dT0_dVb = dEsatL_dVb + dVdsat_dVb + 2.0 * WVCoxRds * tmp4 * tmp3
            * Vgsteff - WVCoxRds * Vgsteff / Vgst2Vtm * (dAbulk_dVb
            * Vdsat + Abulk * dVdsat_dVb);
            dT0_dVd = dEsatL_dVd + dVdsat_dVd - WVCoxRds * Vgsteff / Vgst2Vtm
                    * Abulk * dVdsat_dVd;

            T1 = 2.0 / Pmos_factor - 1.0 + WVCoxRds * Abulk;
            dT1_dVg = -2.0 * tmp1 + WVCoxRds * (Abulk * tmp2 + dAbulk_dVg);
            dT1_dVb = dAbulk_dVb * WVCoxRds + Abulk * WVCoxRds * tmp3;

            Vasat = T0 / T1;
            dVasat_dVg = (dT0_dVg - T0 / T1 * dT1_dVg) / T1;
            dVasat_dVb = (dT0_dVb - T0 / T1 * dT1_dVb) / T1;
            dVasat_dVd = dT0_dVd / T1;

            diffVds = Vds - Vdseff;
            /* Calculate VACLM */
            if (Param.BSIM3v0pclm > 0.0)
            {
                T0 = 1.0 / (Param.BSIM3v0pclm * Abulk * Param.BSIM3v0litl);
                dT0_dVb = -T0 / Abulk * dAbulk_dVb;
                dT0_dVg = -T0 / Abulk * dAbulk_dVg;

                T2 = Vgsteff / EsatL;
                T1 = Leff * (Abulk + T2);
                dT1_dVg = Leff * ((1.0 - T2 * dEsatL_dVg) / EsatL + dAbulk_dVg);
                dT1_dVb = Leff * (dAbulk_dVb - T2 * dEsatL_dVb / EsatL);
                dT1_dVd = -T2 * dEsatL_dVd / Esat;

                VACLM = T0 * T1 * diffVds;
                dVACLM_dVg = T0 * dT1_dVg * diffVds - T0 * T1 * dVdseff_dVg
                           + T1 * diffVds * dT0_dVg;
                dVACLM_dVb = (dT0_dVb * T1 + T0 * dT1_dVb) * diffVds
               - T0 * T1 * dVdseff_dVb;
                dVACLM_dVd = T0 * dT1_dVd * diffVds
               + T0 * T1 * (1.0 - dVdseff_dVd);
            }
            else
            {
                VACLM = MAX_EXP;
                dVACLM_dVd = dVACLM_dVg = dVACLM_dVb = 0.0;
            }

            /* Calculate VADIBL */
            if (Param.BSIM3v0thetaRout > 0.0)
            {
                T0 = Vgst2Vtm * Abulk * Vdsat;
                dT0_dVg = Vgst2Vtm * Abulk * dVdsat_dVg + Abulk * Vdsat
                + Vgst2Vtm * Vdsat * dAbulk_dVg;
                dT0_dVb = Vgst2Vtm * (dAbulk_dVb * Vdsat + Abulk * dVdsat_dVb);
                dT0_dVd = Vgst2Vtm * Abulk * dVdsat_dVd;

                T1 = Vgst2Vtm + Abulk * Vdsat;
                dT1_dVg = 1.0 + Abulk * dVdsat_dVg + Vdsat * dAbulk_dVg;
                dT1_dVb = Abulk * dVdsat_dVb + dAbulk_dVb * Vdsat;
                dT1_dVd = Abulk * dVdsat_dVd;

                T2 = Param.BSIM3v0thetaRout * (1.0 + Param.BSIM3v0pdiblb * Vbseff);
                VADIBL = (Vgst2Vtm - T0 / T1) / T2;
                dVADIBL_dVg = (1.0 - dT0_dVg / T1 + T0 * dT1_dVg / (T1 * T1)) / T2;
                dVADIBL_dVb = ((-dT0_dVb / T1 + T0 * dT1_dVb / (T1 * T1)) - VADIBL
                * Param.BSIM3v0thetaRout * Param.BSIM3v0pdiblb) / T2;
                dVADIBL_dVd = (-dT0_dVd / T1 + T0 * dT1_dVd / (T1 * T1)) / T2;
            }
            else
            {
                VADIBL = MAX_EXP;
                dVADIBL_dVd = dVADIBL_dVg = dVADIBL_dVb = 0.0;
            }

            /* Calculate VA */
            T0 = 1.0 + Param.BSIM3v0pvag * Vgsteff / EsatL;
            dT0_dVg = Param.BSIM3v0pvag * (1.0 - Vgsteff * dEsatL_dVg
                    / EsatL) / EsatL;
            dT0_dVb = -Param.BSIM3v0pvag * Vgsteff * dEsatL_dVb
                    / EsatL / EsatL;
            dT0_dVd = -Param.BSIM3v0pvag * Vgsteff * dEsatL_dVd
                    / EsatL / EsatL;

            tmp1 = VACLM * VACLM;
            tmp2 = VADIBL * VADIBL;
            tmp3 = VACLM + VADIBL;

            T1 = VACLM * VADIBL / tmp3;
            tmp3 *= tmp3;
            dT1_dVg = (tmp1 * dVADIBL_dVg + tmp2 * dVACLM_dVg) / tmp3;
            dT1_dVd = (tmp1 * dVADIBL_dVd + tmp2 * dVACLM_dVd) / tmp3;
            dT1_dVb = (tmp1 * dVADIBL_dVb + tmp2 * dVACLM_dVb) / tmp3;

            Va = Vasat + T0 * T1;
            dVa_dVg = dVasat_dVg + T1 * dT0_dVg + T0 * dT1_dVg;
            dVa_dVd = dVasat_dVd + T1 * dT0_dVd + T0 * dT1_dVd;
            dVa_dVb = dVasat_dVb + T1 * dT0_dVb + T0 * dT1_dVb;

            /* Calculate VASCBE */
            if ((!diffVds.Equals(0.0)) && ((T0 = Param.BSIM3v0pscbe1 * Param.BSIM3v0litl
            / diffVds) > 0.0) && (T0 < EXP_THRESHOLD))
            {
                VASCBE = Leff * Math.Exp(T0) / Param.BSIM3v0pscbe2;
                T1 = T0 * VASCBE / diffVds;
                dVASCBE_dVg = T1 * dVdseff_dVg;
                dVASCBE_dVd = -T1 * (1.0 - dVdseff_dVd);
                dVASCBE_dVb = T1 * dVdseff_dVb;
            }
            else
            {
                VASCBE = MAX_EXP * Leff / Param.BSIM3v0pscbe2;
                dVASCBE_dVg = dVASCBE_dVd = dVASCBE_dVb = 0.0;
            }

            /* Calculate Ids */
            CoxWovL = ModelTemperature.Cox * Weff / Leff;
            beta = ueff * CoxWovL;
            dbeta_dVg = CoxWovL * dueff_dVg + beta * dWeff_dVg / Weff;
            dbeta_dVd = CoxWovL * dueff_dVd;
            dbeta_dVb = CoxWovL * dueff_dVb + beta * dWeff_dVb / Weff;

            T0 = 1.0 - 0.5 * Abulk * Vdseff / Vgst2Vtm;
            dT0_dVg = -0.5 * (Abulk * dVdseff_dVg
            - Abulk * Vdseff / Vgst2Vtm + Vdseff * dAbulk_dVg) / Vgst2Vtm;
            dT0_dVd = -0.5 * Abulk * dVdseff_dVd / Vgst2Vtm;
            dT0_dVb = -0.5 * (Abulk * dVdseff_dVb + dAbulk_dVb * Vdseff)
                    / Vgst2Vtm;

            fgche1 = Vgsteff * T0;
            dfgche1_dVg = Vgsteff * dT0_dVg + T0;
            dfgche1_dVd = Vgsteff * dT0_dVd;
            dfgche1_dVb = Vgsteff * dT0_dVb;

            fgche2 = 1.0 + Vdseff / EsatL;
            dfgche2_dVg = (dVdseff_dVg - Vdseff / EsatL * dEsatL_dVg) / EsatL;
            dfgche2_dVd = (dVdseff_dVd - Vdseff / EsatL * dEsatL_dVd) / EsatL;
            dfgche2_dVb = (dVdseff_dVb - Vdseff / EsatL * dEsatL_dVb) / EsatL;

            gche = beta * fgche1 / fgche2;
            dgche_dVg = (beta * dfgche1_dVg + fgche1 * dbeta_dVg
              - gche * dfgche2_dVg) / fgche2;
            dgche_dVd = (beta * dfgche1_dVd + fgche1 * dbeta_dVd
              - gche * dfgche2_dVd) / fgche2;
            dgche_dVb = (beta * dfgche1_dVb + fgche1 * dbeta_dVb
              - gche * dfgche2_dVb) / fgche2;

            T0 = 1.0 + gche * Rds;
            Idl = gche * Vdseff / T0;

            dIdl_dVg = (gche * dVdseff_dVg + Vdseff * dgche_dVg / T0) / T0
                     - Idl * gche / T0 * dRds_dVg;

            dIdl_dVd = (gche * dVdseff_dVd + Vdseff * dgche_dVd / T0) / T0;
            dIdl_dVb = (gche * dVdseff_dVb + Vdseff * dgche_dVb / T0
                     - Idl * dRds_dVb * gche) / T0;

            T0 = 1.0 + diffVds / Va;
            Idsa = Idl * T0;
            dIdsa_dVg = T0 * dIdl_dVg - Idl * (dVdseff_dVg + diffVds / Va
              * dVa_dVg) / Va;
            dIdsa_dVd = T0 * dIdl_dVd + Idl * (1.0 - dVdseff_dVd - diffVds / Va
              * dVa_dVd) / Va;
            dIdsa_dVb = T0 * dIdl_dVb - Idl * (dVdseff_dVb + diffVds / Va
              * dVa_dVb) / Va;

            T0 = 1.0 + diffVds / VASCBE;
            Ids = Idsa * T0;

            Gm = T0 * dIdsa_dVg - Idsa * (dVdseff_dVg + diffVds / VASCBE
           * dVASCBE_dVg) / VASCBE;
            Gds = T0 * dIdsa_dVd + Idsa * (1.0 - dVdseff_dVd - diffVds / VASCBE
            * dVASCBE_dVd) / VASCBE;
            Gmb = T0 * dIdsa_dVb - Idsa * (dVdseff_dVb + diffVds / VASCBE
            * dVASCBE_dVb) / VASCBE;

            Gds += Gm * dVgsteff_dVd;
            Gmb += Gm * dVgsteff_dVb;
            Gm *= dVgsteff_dVg;
            Gmb *= dVbseff_dVb;

            /* calculate substrate current Isub */
            if ((Param.BSIM3v0alpha0 <= 0.0) || (Param.BSIM3v0beta0 <= 0.0))
            {
                Isub = Gbd = Gbb = Gbg = 0.0;
            }
            else
            {
                T2 = Param.BSIM3v0alpha0 / Leff;
                if (diffVds < 0.0)
                {
                    diffVds = 0.0;
                    Vdseff = Vds;
                } /* added to avoid the hardwrae problem
                                           when Vds=0 */

                if ((diffVds != 0.0) && ((T0 = -Param.BSIM3v0beta0 / diffVds)
                > -EXP_THRESHOLD))
                {
                    T1 = T2 * diffVds * Math.Exp(T0);
                    dT1_dVg = T1 / diffVds * (T0 - 1.0) * dVdseff_dVg;
                    dT1_dVd = T1 / diffVds * (1.0 - T0) * (1.0 - dVdseff_dVd);
                    dT1_dVb = T1 / diffVds * (T0 - 1.0) * dVdseff_dVb;
                }
                else
                {
                    T1 = T2 * diffVds * MIN_EXP;
                    dT1_dVg = -T2 * MIN_EXP * dVdseff_dVg;
                    dT1_dVd = T2 * MIN_EXP * (1.0 - dVdseff_dVd);
                    dT1_dVb = -T2 * MIN_EXP * dVdseff_dVb;
                }
                Isub = T1 * Idsa;
                Gbg = T1 * dIdsa_dVg + Idsa * dT1_dVg;
                Gbd = T1 * dIdsa_dVd + Idsa * dT1_dVd;
                Gbb = T1 * dIdsa_dVb + Idsa * dT1_dVb;

                Gbd += Gbg * dVgsteff_dVd;
                Gbb += Gbg * dVgsteff_dVb;
                Gbg *= dVgsteff_dVg;
                Gbg *= dVbseff_dVb;
            }

            cdrain = Ids;
            this._gds = Gds;
            this._gm = Gm;
            this._gmbs = Gmb;

            this._gbbs = Gbb;
            this._gbgs = Gbg;
            this._gbds = Gbd;

            this._csub = Isub - (Gbb * Vbseff + Gbd * Vds + Gbg * Vgs);

            /* Calculate Qinv for Noise analysis */

            T1 = Vgsteff * (1.0 - 0.5 * Abulk * Vdseff / Vgst2Vtm);
            this._qinv = -ModelTemperature.Cox * Weff * Leff * T1;

            if ((ModelParameters.Xpart < 0) || (!ChargeComputationNeeded))
            {
                qgate = qdrn = qsrc = qbulk = 0.0;
                this._cggb = this._cgsb = this._cgdb = 0.0;
                this._cdgb = this._cdsb = this._cddb = 0.0;
                this._cbgb = this._cbsb = this._cbdb = 0.0;
                this._cqdb = this._cqsb = this._cqgb
                            = this._cqbb = 0.0;
                this._gtau = 0.0;
                goto finished;
            }
            else
            {
                if (Vbseff < 0.0)
                {
                    VbseffCV = Vbseff;
                    dVbseffCV_dVb = 1.0;
                }
                else
                {
                    VbseffCV = Param.BSIM3v0phi - Phis;
                    dVbseffCV_dVb = -dPhis_dVb;
                }

                CoxWL = ModelTemperature.Cox * Param.BSIM3v0weffCV
              * Param.BSIM3v0leffCV;
                Vfb = Vth - Param.BSIM3v0phi - Param.BSIM3v0k1 * sqrtPhis;

                dVfb_dVb = dVth_dVb - Param.BSIM3v0k1 * dsqrtPhis_dVb;
                dVfb_dVd = dVth_dVd;

                if ((VgstNVt > -EXP_THRESHOLD) && (VgstNVt < EXP_THRESHOLD))
                {
                    ExpVgst *= ExpVgst;
                    Vgsteff = n * Vtm * Math.Log(1.0 + ExpVgst);
                    dVgsteff_dVg = ExpVgst / (1.0 + ExpVgst);
                    dVgsteff_dVd = -dVgsteff_dVg * (dVth_dVd + (Vgs_eff - Vth)
                     / n * dn_dVd) + Vgsteff / n * dn_dVd;
                    dVgsteff_dVb = -dVgsteff_dVg * (dVth_dVb + (Vgs_eff - Vth)
                     / n * dn_dVb) + Vgsteff / n * dn_dVb;
                    dVgsteff_dVg *= dVgs_eff_dVg;
                }

                if (ModelParameters.CapMod.Value == 1)
                {
                    Arg1 = Vgs_eff - VbseffCV - Vfb;

                    if (Arg1 <= 0.0)
                    {
                        qgate = CoxWL * (Arg1 - Vgsteff);
                        Cgg = CoxWL * (dVgs_eff_dVg - dVgsteff_dVg);
                        Cgd = -CoxWL * (dVfb_dVd + dVgsteff_dVd);
                        Cgb = -CoxWL * (dVfb_dVb + dVbseffCV_dVb + dVgsteff_dVb);
                    }
                    else
                    {
                        T0 = 0.5 * Param.BSIM3v0k1;
                        T1 = Math.Sqrt(T0 * T0 + Arg1 - Vgsteff);
                        qgate = CoxWL * Param.BSIM3v0k1 * (T1 - T0);

                        T2 = CoxWL * T0 / T1;
                        Cgg = T2 * (dVgs_eff_dVg - dVgsteff_dVg);
                        Cgd = -T2 * (dVfb_dVd + dVgsteff_dVd);
                        Cgb = -T2 * (dVfb_dVb + dVbseffCV_dVb + dVgsteff_dVb);
                    }
                    qbulk = -qgate;
                    Cbg = -Cgg;
                    Cbd = -Cgd;
                    Cbb = -Cgb;

                    One_Third_CoxWL = CoxWL / 3.0;
                    Two_Third_CoxWL = 2.0 * One_Third_CoxWL;
                    AbulkCV = Abulk0 * Param.BSIM3v0abulkCVfactor;
                    dAbulkCV_dVb = Param.BSIM3v0abulkCVfactor * dAbulk0_dVb;
                    VdsatCV = Vgsteff / AbulkCV;
                    if (VdsatCV <= Vds)
                    {
                        dVdsatCV_dVg = 1.0 / AbulkCV;
                        dVdsatCV_dVb = -VdsatCV * dAbulkCV_dVb / AbulkCV;
                        T0 = Vgsteff - VdsatCV / 3.0;
                        dT0_dVg = 1.0 - dVdsatCV_dVg / 3.0;
                        dT0_dVb = -dVdsatCV_dVb / 3.0;
                        qgate += CoxWL * T0;
                        Cgg1 = CoxWL * dT0_dVg;
                        Cgb1 = CoxWL * dT0_dVb + Cgg1 * dVgsteff_dVb;
                        Cgd1 = Cgg1 * dVgsteff_dVd;
                        Cgg1 *= dVgsteff_dVg;
                        Cgg += Cgg1;
                        Cgb += Cgb1;
                        Cgd += Cgd1;

                        T0 = VdsatCV - Vgsteff;
                        dT0_dVg = dVdsatCV_dVg - 1.0;
                        dT0_dVb = dVdsatCV_dVb;
                        qbulk += One_Third_CoxWL * T0;
                        Cbg1 = One_Third_CoxWL * dT0_dVg;
                        Cbb1 = One_Third_CoxWL * dT0_dVb + Cbg1 * dVgsteff_dVb;
                        Cbd1 = Cbg1 * dVgsteff_dVd;
                        Cbg1 *= dVgsteff_dVg;
                        Cbg += Cbg1;
                        Cbb += Cbb1;
                        Cbd += Cbd1;

                        if (ModelParameters.Xpart > 0.5)
                            T0 = -Two_Third_CoxWL;
                        else if (ModelParameters.Xpart < 0.5)
                            T0 = -0.4 * CoxWL;
                        else
                            T0 = -One_Third_CoxWL;

                        qsrc = T0 * Vgsteff;
                        Csg = T0 * dVgsteff_dVg;
                        Csb = T0 * dVgsteff_dVb;
                        Csd = T0 * dVgsteff_dVd;
                        Cgb *= dVbseff_dVb;
                        Cbb *= dVbseff_dVb;
                        Csb *= dVbseff_dVb;
                    }
                    else
                    {
                        T0 = AbulkCV * Vds;
                        T1 = 12.0 * (Vgsteff - 0.5 * T0);
                        T2 = Vds / T1;
                        T3 = T0 * T2;
                        dT3_dVg = -12.0 * T2 * T2 * AbulkCV;
                        dT3_dVd = 6.0 * T0 * (4.0 * Vgsteff - T0) / T1 / T1 - 0.5;
                        dT3_dVb = 12.0 * T2 * T2 * dAbulkCV_dVb * Vgsteff;

                        qgate += CoxWL * (Vgsteff - 0.5 * Vds + T3);
                        Cgg1 = CoxWL * (1.0 + dT3_dVg);
                        Cgb1 = CoxWL * dT3_dVb + Cgg1 * dVgsteff_dVb;
                        Cgd1 = CoxWL * dT3_dVd + Cgg1 * dVgsteff_dVd;
                        Cgg1 *= dVgsteff_dVg;
                        Cgg += Cgg1;
                        Cgb += Cgb1;
                        Cgd += Cgd1;

                        qbulk += CoxWL * (1.0 - AbulkCV) * (0.5 * Vds - T3);
                        Cbg1 = -CoxWL * ((1.0 - AbulkCV) * dT3_dVg);
                        Cbb1 = -CoxWL * ((1.0 - AbulkCV) * dT3_dVb
                             + (0.5 * Vds - T3) * dAbulkCV_dVb)
                             + Cbg1 * dVgsteff_dVb;
                        Cbd1 = -CoxWL * (1.0 - AbulkCV) * dT3_dVd
                             + Cbg1 * dVgsteff_dVd;
                        Cbg1 *= dVgsteff_dVg;
                        Cbg += Cbg1;
                        Cbb += Cbb1;
                        Cbd += Cbd1;

                        if (ModelParameters.Xpart > 0.5)
                        {   /* 0/100 Charge petition model */
                            T1 = T1 + T1;
                            qsrc = -CoxWL * (0.5 * Vgsteff + 0.25 * T0
                     - T0 * T0 / T1);
                            Csg = -CoxWL * (0.5 + 24.0 * T0 * Vds / T1 / T1
                            * AbulkCV);
                            Csb = -CoxWL * (0.25 * Vds * dAbulkCV_dVb
                    - 12.0 * T0 * Vds / T1 / T1 * (4.0 * Vgsteff - T0)
                    * dAbulkCV_dVb) + Csg * dVgsteff_dVb;
                            Csd = -CoxWL * (0.25 * AbulkCV - 12.0 * AbulkCV * T0
                            / T1 / T1 * (4.0 * Vgsteff - T0))
                            + Csg * dVgsteff_dVd;
                            Csg *= dVgsteff_dVg;
                        }
                        else if (ModelParameters.Xpart < 0.5)
                        {   /* 40/60 Charge petition model */
                            T1 = T1 / 12.0;
                            T2 = 0.5 * CoxWL / (T1 * T1);
                            T3 = Vgsteff * (2.0 * T0 * T0 / 3.0 + Vgsteff
                           * (Vgsteff - 4.0 * T0 / 3.0))
                           - 2.0 * T0 * T0 * T0 / 15.0;
                            qsrc = -T2 * T3;
                            T4 = 4.0 / 3.0 * Vgsteff * (Vgsteff - T0)
                           + 0.4 * T0 * T0;
                            Csg = -2.0 * qsrc / T1 - T2 * (Vgsteff * (3.0
                            * Vgsteff - 8.0 * T0 / 3.0)
                            + 2.0 * T0 * T0 / 3.0);
                            Csb = (qsrc / T1 * Vds + T2 * T4 * Vds) * dAbulkCV_dVb
                                + Csg * dVgsteff_dVb;
                            Csd = (qsrc / T1 + T2 * T4) * AbulkCV
                            + Csg * dVgsteff_dVd;
                            Csg *= dVgsteff_dVg;
                        }
                        else
                        {   /* 50/50 Charge petition model */
                            qsrc = -0.5 * (qgate + qbulk);
                            Csg = -0.5 * (Cgg1 + Cbg1);
                            Csb = -0.5 * (Cgb1 + Cbb1);
                            Csd = -0.5 * (Cgd1 + Cbd1);
                        }
                        Cgb *= dVbseff_dVb;
                        Cbb *= dVbseff_dVb;
                        Csb *= dVbseff_dVb;
                    }
                    qdrn = -(qgate + qbulk + qsrc);
                    this._cggb = Cgg;
                    this._cgsb = -(Cgg + Cgd + Cgb);
                    this._cgdb = Cgd;
                    this._cdgb = -(Cgg + Cbg + Csg);
                    this._cdsb = (Cgg + Cgd + Cgb + Cbg + Cbd + Cbb
                            + Csg + Csd + Csb);
                    this._cddb = -(Cgd + Cbd + Csd);
                    this._cbgb = Cbg;
                    this._cbsb = -(Cbg + Cbd + Cbb);
                    this._cbdb = Cbd;
                }
                else
                {
                    V3 = Vfb - Vgs_eff + VbseffCV - DELTA_3;
                    if (Vfb <= 0.0)
                    {
                        T0 = Math.Sqrt(V3 * V3 - 4.0 * DELTA_3 * Vfb);
                        T2 = -DELTA_3 / T0;
                    }
                    else
                    {
                        T0 = Math.Sqrt(V3 * V3 + 4.0 * DELTA_3 * Vfb);
                        T2 = DELTA_3 / T0;
                    }

                    T1 = 0.5 * (1.0 + V3 / T0);
                    Vfbeff = Vfb - 0.5 * (V3 + T0);
                    dVfbeff_dVd = (1.0 - T1 - T2) * dVfb_dVd;
                    dVfbeff_dVg = T1 * dVgs_eff_dVg;
                    dVfbeff_dVb = (1.0 - T1 - T2) * dVfb_dVb
                                - T1 * dVbseffCV_dVb;
                    Qac0 = CoxWL * (Vfbeff - Vfb);
                    dQac0_dVg = CoxWL * dVfbeff_dVg;
                    dQac0_dVd = CoxWL * (dVfbeff_dVd - dVfb_dVd);
                    dQac0_dVb = CoxWL * (dVfbeff_dVb - dVfb_dVb);

                    T0 = 0.5 * Param.BSIM3v0k1;
                    T1 = Math.Sqrt(T0 * T0 + Vgs_eff - Vfbeff - VbseffCV - Vgsteff);

                    Qsub0 = CoxWL * Param.BSIM3v0k1 * (T1 - T0);

                    T2 = CoxWL * T0 / T1;
                    dQsub0_dVg = T2 * (dVgs_eff_dVg - dVfbeff_dVg - dVgsteff_dVg);
                    dQsub0_dVd = -T2 * (dVfbeff_dVd + dVgsteff_dVd);
                    dQsub0_dVb = -T2 * (dVfbeff_dVb + dVbseffCV_dVb
                                 + dVgsteff_dVb);

                    One_Third_CoxWL = CoxWL / 3.0;
                    Two_Third_CoxWL = 2.0 * One_Third_CoxWL;
                    AbulkCV = Abulk0 * Param.BSIM3v0abulkCVfactor;
                    dAbulkCV_dVb = Param.BSIM3v0abulkCVfactor * dAbulk0_dVb;
                    VdsatCV = Vgsteff / AbulkCV;

                    V4 = VdsatCV - Vds - DELTA_4;
                    T0 = Math.Sqrt(V4 * V4 + 4.0 * DELTA_4 * VdsatCV);
                    VdseffCV = VdsatCV - 0.5 * (V4 + T0);
                    T1 = 0.5 * (1.0 + V4 / T0);
                    T2 = DELTA_4 / T0;
                    T3 = (1.0 - T1 - T2) / AbulkCV;
                    dVdseffCV_dVg = T3;
                    dVdseffCV_dVd = T1;
                    dVdseffCV_dVb = -T3 * VdsatCV * dAbulkCV_dVb;

                    T0 = AbulkCV * VdseffCV;
                    T1 = 12.0 * (Vgsteff - 0.5 * T0 + 1e-20);
                    T2 = VdseffCV / T1;
                    T3 = T0 * T2;

                    T4 = (1.0 - 12.0 * T2 * T2 * AbulkCV);
                    T5 = (6.0 * T0 * (4.0 * Vgsteff - T0) / (T1 * T1) - 0.5);
                    T6 = 12.0 * T2 * T2 * Vgsteff;

                    qgate = CoxWL * (Vgsteff - 0.5 * VdseffCV + T3);
                    Cgg1 = CoxWL * (T4 + T5 * dVdseffCV_dVg);
                    Cgd1 = CoxWL * T5 * dVdseffCV_dVd + Cgg1 * dVgsteff_dVd;
                    Cgb1 = CoxWL * (T5 * dVdseffCV_dVb + T6 * dAbulkCV_dVb)
                 + Cgg1 * dVgsteff_dVb;
                    Cgg1 *= dVgsteff_dVg;

                    T7 = 1.0 - AbulkCV;
                    qbulk = CoxWL * T7 * (0.5 * VdseffCV - T3);
                    T4 = -T7 * (T4 - 1.0);
                    T5 = -T7 * T5;
                    T6 = -(T7 * T6 + (0.5 * VdseffCV - T3));
                    Cbg1 = CoxWL * (T4 + T5 * dVdseffCV_dVg);
                    Cbd1 = CoxWL * T5 * dVdseffCV_dVd + Cbg1 * dVgsteff_dVd;
                    Cbb1 = CoxWL * (T5 * dVdseffCV_dVb + T6 * dAbulkCV_dVb)
                 + Cbg1 * dVgsteff_dVb;
                    Cbg1 *= dVgsteff_dVg;

                    if (ModelParameters.Xpart > 0.5)
                    {   /* 0/100 Charge petition model */
                        T1 = T1 + T1;
                        qsrc = -CoxWL * (0.5 * Vgsteff + 0.25 * T0
                 - T0 * T0 / T1);
                        T7 = (4.0 * Vgsteff - T0) / (T1 * T1);
                        T4 = -(0.5 + 24.0 * T0 * T0 / (T1 * T1));
                        T5 = -(0.25 * AbulkCV - 12.0 * AbulkCV * T0 * T7);
                        T6 = -(0.25 * VdseffCV - 12.0 * T0 * VdseffCV * T7);
                        Csg = CoxWL * (T4 + T5 * dVdseffCV_dVg);
                        Csd = CoxWL * T5 * dVdseffCV_dVd + Csg * dVgsteff_dVd;
                        Csb = CoxWL * (T5 * dVdseffCV_dVb + T6 * dAbulkCV_dVb)
                + Csg * dVgsteff_dVb;
                        Csg *= dVgsteff_dVg;
                    }
                    else if (ModelParameters.Xpart < 0.5)
                    {   /* 40/60 Charge petition model */
                        T1 = T1 / 12.0;
                        T2 = 0.5 * CoxWL / (T1 * T1);
                        T3 = Vgsteff * (2.0 * T0 * T0 / 3.0 + Vgsteff
                       * (Vgsteff - 4.0 * T0 / 3.0))
                       - 2.0 * T0 * T0 * T0 / 15.0;
                        qsrc = -T2 * T3;
                        T7 = 4.0 / 3.0 * Vgsteff * (Vgsteff - T0)
                       + 0.4 * T0 * T0;
                        T4 = -2.0 * qsrc / T1 - T2 * (Vgsteff * (3.0
                       * Vgsteff - 8.0 * T0 / 3.0)
                       + 2.0 * T0 * T0 / 3.0);
                        T5 = (qsrc / T1 + T2 * T7) * AbulkCV;
                        T6 = (qsrc / T1 * VdseffCV + T2 * T7 * VdseffCV);
                        Csg = (T4 + T5 * dVdseffCV_dVg);
                        Csd = T5 * dVdseffCV_dVd + Csg * dVgsteff_dVd;
                        Csb = (T5 * dVdseffCV_dVb + T6 * dAbulkCV_dVb)
                + Csg * dVgsteff_dVb;
                        Csg *= dVgsteff_dVg;
                    }
                    else
                    {   /* 50/50 Charge petition model */
                        qsrc = -0.5 * (qgate + qbulk);
                        Csg = -0.5 * (Cgg1 + Cbg1);
                        Csb = -0.5 * (Cgb1 + Cbb1);
                        Csd = -0.5 * (Cgd1 + Cbd1);
                    }

                    qgate += Qac0 + Qsub0;
                    qbulk -= (Qac0 + Qsub0);
                    qdrn = -(qgate + qbulk + qsrc);

                    Cgg = dQac0_dVg + dQsub0_dVg + Cgg1;
                    Cgd = dQac0_dVd + dQsub0_dVd + Cgd1;
                    Cgb = dQac0_dVb + dQsub0_dVb + Cgb1;

                    Cbg = Cbg1 - dQac0_dVg - dQsub0_dVg;
                    Cbd = Cbd1 - dQac0_dVd - dQsub0_dVd;
                    Cbb = Cbb1 - dQac0_dVb - dQsub0_dVb;

                    Cgb *= dVbseff_dVb;
                    Cbb *= dVbseff_dVb;
                    Csb *= dVbseff_dVb;

                    this._cggb = Cgg;
                    this._cgsb = -(Cgg + Cgd + Cgb);
                    this._cgdb = Cgd;
                    this._cdgb = -(Cgg + Cbg + Csg);
                    this._cdsb = (Cgg + Cgd + Cgb + Cbg + Cbd + Cbb
                            + Csg + Csd + Csb);
                    this._cddb = -(Cgd + Cbd + Csd);
                    this._cbgb = Cbg;
                    this._cbsb = -(Cbg + Cbd + Cbb);
                    this._cbdb = Cbd;

                }

                /* Non-quasi-static Model */

                if (Parameters.NqsMod.Value != 0)
                {
                    qcheq = -qbulk - qgate;
                    qbulk = qgate = qdrn = qsrc = 0.0;

                    this._cqgb = -(this._cggb + this._cbgb);
                    this._cqdb = -(this._cgdb + this._cbdb);
                    this._cqsb = -(this._cgsb + this._cbsb);
                    this._cqbb = this._cggb + this._cgdb
                                    + this._cgsb + this._cbgb
                                    + this._cbdb + this._cbsb;

                    this._cggb = this._cgsb = this._cgdb = 0.0;
                    this._cdgb = this._cdsb = this._cddb = 0.0;
                    this._cbgb = this._cbsb = this._cbdb = 0.0;

                    T0 = Param.BSIM3v0leffCV * Param.BSIM3v0leffCV;
                    this._tconst = Param.BSIM3v0u0temp * Param.BSIM3v0elm
                      / CoxWL / T0;

                    if (qcheq == 0.0)
                        this._tconst = 0.0;
                    else if (qcheq < 0.0)
                        this._tconst = -this._tconst;

                    gtau_drift = Math.Abs(this._tconst * qcheq);
                    gtau_diff = 16.0 * Param.BSIM3v0u0temp * ModelTemperature.Vtm / T0;

                    this._gtau = gtau_drift + gtau_diff;

                    this._qcheq.Value = qcheq;
                    this._qcheq.Derive();
                }
                else
                {
                    this._cqgb = this._cqdb = this._cqsb
                                    = this._cqbb = 0.0;
                    this._gtau = 0.0;
                }
            }

        finished: /* returning Values to Calling Routine */
            /*
             *  COMPUTE EQUIVALENT DRAIN CURRENT SOURCE
             */
            this._cd = this._mode * cdrain - this._cbd;
            if (ChargeComputationNeeded)
            {   /*  charge storage elements
               *  bulk-drain and bulk-source depletion capacitances
               *  czbd : zero bias drain junction capacitance
               *  czbs : zero bias source junction capacitance
               *  czbdsw:zero bias drain junction sidewall capacitance
               *  czbssw:zero bias source junction sidewall capacitance
               */

                czbd = ModelParameters.UnitAreaJctCap * Parameters.DrainArea;
                czbs = ModelParameters.UnitAreaJctCap * Parameters.SourceArea;
                czbdsw = ModelParameters.UnitLengthSidewallJctCap
               * Parameters.DrainPerimeter;
                czbssw = ModelParameters.UnitLengthSidewallJctCap
               * Parameters.SourcePerimeter;
                PhiB = ModelParameters.BulkJctPotential;
                PhiBSW = ModelParameters.SidewallJctPotential;
                MJ = ModelParameters.BulkJctBotGradingCoeff;
                MJSW = ModelParameters.BulkJctSideGradingCoeff;

                /* Source Bulk Junction */
                if (vbs == 0.0)
                {
                    this._qbs = 0.0;
                    this._capbs = czbs + czbssw;
                }
                else if (vbs < 0.0)
                {
                    if (czbs > 0.0)
                    {
                        arg = 1.0 - vbs / PhiB;
                        if (MJ == 0.5)
                            sarg = 1.0 / Math.Sqrt(arg);
                        else
                            sarg = Math.Exp(-MJ * Math.Log(arg));
                        this._qbs = PhiB * czbs
                      * (1.0 - arg * sarg) / (1.0 - MJ);
                        this._capbs = czbs * sarg;
                    }
                    else
                    {
                        this._qbs = 0.0;
                        this._capbs = 0.0;
                    }
                    if (czbssw > 0.0)
                    {
                        arg = 1.0 - vbs / PhiBSW;
                        if (MJSW == 0.5)
                            sarg = 1.0 / Math.Sqrt(arg);
                        else
                            sarg = Math.Exp(-MJSW * Math.Log(arg));
                        this._qbs += PhiBSW * czbssw
                      * (1.0 - arg * sarg) / (1.0 - MJSW);
                        this._capbs += czbssw * sarg;
                    }
                }
                else
                {
                    this._qbs = vbs * (czbs + czbssw)
                             + vbs * vbs * (czbs * MJ * 0.5 / PhiB
                                         + czbssw * MJSW * 0.5 / PhiBSW);
                    this._capbs = czbs + czbssw + vbs * (czbs * MJ / PhiB
                             + czbssw * MJSW / PhiBSW);
                }

                /* Drain Bulk Junction */
                if (vbd == 0.0)
                {
                    this._qbd = 0.0;
                    this._capbd = czbd + czbdsw;
                }
                else if (vbd < 0.0)
                {
                    if (czbd > 0.0)
                    {
                        arg = 1.0 - vbd / PhiB;
                        if (MJ == 0.5)
                            sarg = 1.0 / Math.Sqrt(arg);
                        else
                            sarg = Math.Exp(-MJ * Math.Log(arg));
                        this._qbd = PhiB * czbd
                         * (1.0 - arg * sarg) / (1.0 - MJ);
                        this._capbd = czbd * sarg;
                    }
                    else
                    {
                        this._qbd = 0.0;
                        this._capbd = 0.0;
                    }
                    if (czbdsw > 0.0)
                    {
                        arg = 1.0 - vbd / PhiBSW;
                        if (MJSW == 0.5)
                            sarg = 1.0 / Math.Sqrt(arg);
                        else
                            sarg = Math.Exp(-MJSW * Math.Log(arg));
                        this._qbd += PhiBSW * czbdsw
                             * (1.0 - arg * sarg) / (1.0 - MJSW);
                        this._capbd += czbdsw * sarg;
                    }
                }
                else
                {
                    this._qbd = vbd * (czbd + czbdsw)
                        + vbd * vbd * (czbd * MJ * 0.5 / PhiB
                                    + czbdsw * MJSW * 0.5 / PhiBSW);
                    this._capbd = czbd + czbdsw + vbd * (czbd * MJ / PhiB
                    + czbdsw * MJSW / PhiBSW);
                }
            }

            /*
             *  check convergence
             */
            if (!Parameters.Off || _iteration.Mode != IterationModes.Fix)
            {
                if (Check)
                    _iteration.IsConvergent = false;
            }
            this._vbs = vbs;
            this._vbd = vbd;
            this._vgs = vgs;
            this._vds = vds;

            /* bulk and channel charge plus overlaps */

            if (!ChargeComputationNeeded)
                goto line850;

            if (_method != null)
                ag0 = _method.Slope;
            else
                ag0 = 0.0;

            if (ModelParameters.CapMod.Value == 1)
            {
                if (vgd < 0.0)
                {
                    T1 = Math.Sqrt(1.0 - 4.0 * vgd / Param.BSIM3v0ckappa);
                    cgdo = Param.BSIM3v0cgdo + Param.BSIM3v0weffCV
                     * Param.BSIM3v0cgdl / T1;
                    qgdo = Param.BSIM3v0cgdo * vgd - Param.BSIM3v0weffCV * 0.5
                     * Param.BSIM3v0cgdl * Param.BSIM3v0ckappa * (T1 - 1.0);
                }
                else
                {
                    cgdo = Param.BSIM3v0cgdo + Param.BSIM3v0weffCV
                     * Param.BSIM3v0cgdl;
                    qgdo = (Param.BSIM3v0weffCV * Param.BSIM3v0cgdl
                     + Param.BSIM3v0cgdo) * vgd;
                }

                if (vgs < 0.0)
                {
                    T1 = Math.Sqrt(1.0 - 4.0 * vgs / Param.BSIM3v0ckappa);
                    cgso = Param.BSIM3v0cgso + Param.BSIM3v0weffCV
                     * Param.BSIM3v0cgsl / T1;
                    qgso = Param.BSIM3v0cgso * vgs - Param.BSIM3v0weffCV * 0.5
                     * Param.BSIM3v0cgsl * Param.BSIM3v0ckappa * (T1 - 1.0);
                }
                else
                {
                    cgso = Param.BSIM3v0cgso + Param.BSIM3v0weffCV
                     * Param.BSIM3v0cgsl;
                    qgso = (Param.BSIM3v0weffCV * Param.BSIM3v0cgsl
                     + Param.BSIM3v0cgso) * vgs;
                }
            }
            else
            {
                T0 = vgd + DELTA_1;
                T1 = Math.Sqrt(T0 * T0 + 4.0 * DELTA_1);
                T2 = 0.5 * (T0 - T1);

                T3 = Param.BSIM3v0weffCV * Param.BSIM3v0cgdl;
                T4 = Math.Sqrt(1.0 - 4.0 * T2 / Param.BSIM3v0ckappa);
                cgdo = Param.BSIM3v0cgdo + T3 - T3 * (1.0 - 1.0 / T4)
                 * (0.5 - 0.5 * T0 / T1);
                qgdo = (Param.BSIM3v0cgdo + T3) * vgd - T3 * (T2
                 + 0.5 * Param.BSIM3v0ckappa * (T4 - 1.0));

                T0 = vgs + DELTA_1;
                T1 = Math.Sqrt(T0 * T0 + 4.0 * DELTA_1);
                T2 = 0.5 * (T0 - T1);
                T3 = Param.BSIM3v0weffCV * Param.BSIM3v0cgsl;
                T4 = Math.Sqrt(1.0 - 4.0 * T2 / Param.BSIM3v0ckappa);
                cgso = Param.BSIM3v0cgso + T3 - T3 * (1.0 - 1.0 / T4)
                 * (0.5 - 0.5 * T0 / T1);
                qgso = (Param.BSIM3v0cgso + T3) * vgs - T3 * (T2
                 + 0.5 * Param.BSIM3v0ckappa * (T4 - 1.0));
            }

            if (this._mode > 0)
            {
                gcdgb = (this._cdgb - cgdo) * ag0;
                gcddb = (this._cddb + this._capbd + cgdo) * ag0;
                gcdsb = this._cdsb * ag0;
                gcsgb = -(this._cggb + this._cbgb + this._cdgb
                  + cgso) * ag0;
                gcsdb = -(this._cgdb + this._cbdb + this._cddb)
              * ag0;
                gcssb = (this._capbs + cgso - (this._cgsb
              + this._cbsb + this._cdsb)) * ag0;
                gcggb = (this._cggb + cgdo + cgso + Param.BSIM3v0cgbo) * ag0;
                gcgdb = (this._cgdb - cgdo) * ag0;
                gcgsb = (this._cgsb - cgso) * ag0;
                gcbgb = (this._cbgb - Param.BSIM3v0cgbo) * ag0;
                gcbdb = (this._cbdb - this._capbd) * ag0;
                gcbsb = (this._cbsb - this._capbs) * ag0;

                gcqgb = this._cqgb * ag0;
                gcqdb = this._cqdb * ag0;
                gcqsb = this._cqsb * ag0;
                gcqbb = this._cqbb * ag0;

                T0 = this._tconst * qdef;
                this._gtg = T0 * this._cqgb;
                this._gtb = T0 * this._cqbb;
                this._gtd = T0 * this._cqdb;
                this._gts = T0 * this._cqsb;

                sxpart = 0.6;
                dxpart = 0.4;

                /* compute total terminal charge */
                qgd = qgdo;
                qgs = qgso;
                qgb = Param.BSIM3v0cgbo * vgb;
                qgate += qgd + qgs + qgb;
                qbulk -= qgb;
                qdrn -= qgd;
                qsrc = -(qgate + qbulk + qdrn);
            }
            else
            {
                gcsgb = (this._cdgb - cgso) * ag0;
                gcsdb = this._cdsb * ag0;
                gcssb = (this._cddb + this._capbs + cgso) * ag0;

                gcdgb = -(this._cggb + this._cbgb + this._cdgb
              + cgdo) * ag0;
                gcdsb = -(this._cgdb + this._cbdb + this._cddb)
              * ag0;
                gcddb = (this._capbd + cgdo - (this._cgsb
              + this._cbsb + this._cdsb)) * ag0;
                gcggb = (this._cggb + cgdo + cgso + Param.BSIM3v0cgbo) * ag0;
                gcgdb = (this._cgsb - cgdo) * ag0;
                gcgsb = (this._cgdb - cgso) * ag0;
                gcbgb = (this._cbgb - Param.BSIM3v0cgbo) * ag0;
                gcbdb = (this._cbsb - this._capbd) * ag0;
                gcbsb = (this._cbdb - this._capbs) * ag0;

                gcqgb = this._cqgb * ag0;
                gcqdb = this._cqsb * ag0;
                gcqsb = this._cqdb * ag0;
                gcqbb = this._cqbb * ag0;

                T0 = this._tconst * qdef;
                this._gtg = T0 * this._cqgb;
                this._gtb = T0 * this._cqbb;
                this._gtd = T0 * this._cqdb;
                this._gts = T0 * this._cqsb;

                dxpart = 0.6;
                sxpart = 0.4;

                /* compute total terminal charge */
                qgd = qgdo;
                qgs = qgso;
                qgb = Param.BSIM3v0cgbo * vgb;
                qgate += qgd + qgs + qgb;
                qbulk -= qgb;
                qsrc = qdrn - qgs;
                qdrn = -(qgate + qbulk + qsrc);
            }

            this._cgdo = cgdo;
            this._cgso = cgso;

            /* added by Mansun 11/1/93 */


            if (Parameters.NqsMod.Value != 0)
            {
                this._qcdump.Value = qdef;
                _qcdump.Derive();
                var r = this._qcdump.GetContributions(1.0);
                gqdef = r.Jacobian;
                cqdef = r.Rhs;
            }
            else
            {
                gqdef = cqdef = 0.0;
            }

            /* End added by Mansun 11/1/93 */

            this._qg.Value = qgate;
            this._qd.Value = qdrn - this._qbd;
            this._qb.Value = qbulk + this._qbd + this._qbs;

            /* store small signal parameters */
            if (InitializeSmallSignal)
                return;
            if (!ChargeComputationNeeded)
                goto line850;

            this._qb.Derive();
            this._qg.Derive();
            this._qd.Derive();

            goto line860;

        line850:
            /* initialize to zero charge conductance and current */
            ceqqg = ceqqb = ceqqd = 0.0;

            cqcheq = cqdef = 0.0;

            gcdgb = gcddb = gcdsb = 0.0;
            gcsgb = gcsdb = gcssb = 0.0;
            gcggb = gcgdb = gcgsb = 0.0;
            gcbgb = gcbdb = gcbsb = 0.0;

            gcqgb = gcqdb = gcqsb = gcqbb = 0.0;
            this._gtg = this._gtd = this._gts = this._gtb = 0.0;
            gqdef = 0.0;
            sxpart = (1.0 - (dxpart = (this._mode > 0) ? 0.4 : 0.6));
            if (Parameters.NqsMod.Value != 0)
                this._gtau = 16.0 * Param.BSIM3v0u0temp * ModelTemperature.Vtm
                                / Leff / Leff;
            else
                this._gtau = 0.0;

            goto line900;

        line860:
            /* evaluate equivalent charge current */

            cqgate = this._qg.Derivative;
            cqbulk = this._qb.Derivative;
            cqdrn = this._qd.Derivative;

            ceqqg = cqgate - gcggb * vgb + gcgdb * vbd + gcgsb * vbs
                    + (this._gtg * vgb - this._gtd * vbd - this._gts * vbs);
            ceqqb = cqbulk - gcbgb * vgb + gcbdb * vbd + gcbsb * vbs;
            ceqqd = cqdrn - gcdgb * vgb + gcddb * vbd + gcdsb * vbs
                    - dxpart * (this._gtg * vgb - this._gtd * vbd
                    - this._gts * vbs);

            cqcheq = this._qcheq.Derivative
                    - (gcqgb * vgb - gcqdb * vbd - gcqsb * vbs)
                    + (this._gtg * vgb - this._gtd * vbd - this._gts * vbs);

        /*
         *  load current vector
         */
        line900:

            if (this._mode >= 0)
            {
                Gm = this._gm;
                Gmbs = this._gmbs;
                FwdSum = Gm + Gmbs;
                RevSum = 0.0;
                cdreq = ModelParameters.Type * (cdrain - this._gds * vds
              - Gm * vgs - Gmbs * vbs);
                ceqbs = -this._csub;
                ceqbd = 0.0;

                gbspsp = -this._gbds - this._gbgs - this._gbbs;
                gbbdp = -this._gbds;
                gbbsp = this._gbds + this._gbgs + this._gbbs;
                gbspg = this._gbgs;
                gbspb = this._gbbs;
                gbspdp = this._gbds;
                gbdpdp = 0.0;
                gbdpg = 0.0;
                gbdpb = 0.0;
                gbdpsp = 0.0;
            }
            else
            {
                Gm = -this._gm;
                Gmbs = -this._gmbs;
                FwdSum = 0.0;
                RevSum = -(Gm + Gmbs);
                cdreq = -ModelParameters.Type * (cdrain + this._gds * vds
                      + Gm * vgd + Gmbs * vbd);
                ceqbs = 0.0;
                ceqbd = -this._csub;

                gbspsp = 0.0;
                gbbdp = this._gbds + this._gbgs + this._gbbs;
                gbbsp = -this._gbds;
                gbspg = 0.0;
                gbspb = 0.0;
                gbspdp = 0.0;
                gbdpdp = -this._gbds - this._gbgs - this._gbbs;
                gbdpg = this._gbgs;
                gbdpb = this._gbbs;
                gbdpsp = this._gbds;
            }

            if (ModelParameters.Type > 0)
            {
                ceqbs += (this._cbs - (this._gbs - _iteration.Gmin) * vbs);
                ceqbd += (this._cbd - (this._gbd - _iteration.Gmin) * vbd);
                /*
                ceqqg = ceqqg;
                ceqqb = ceqqb;
                ceqqd = ceqqd;
                cqcheq = cqcheq;
                */
            }
            else
            {
                ceqbs = -ceqbs - (this._cbs - (this._gbs
                  - _iteration.Gmin) * vbs);
                ceqbd = -ceqbd - (this._cbd - (this._gbd
              - _iteration.Gmin) * vbd);
                ceqqg = -ceqqg;
                ceqqb = -ceqqb;
                ceqqd = -ceqqd;
                cqcheq = -cqcheq;
            }

           this._gPtr.Value -= m * ceqqg;
            this._bPtr.Value -= m * (ceqbs + ceqbd + ceqqb);
            this._dpPtr.Value += m * (ceqbd - cdreq - ceqqd);
            this._spPtr.Value += m * (cdreq + ceqbs + ceqqg)
                           + ceqqb + ceqqd;

            this._qPtr.Value += m * (cqcheq - cqdef);

            /*
             *  load y matrix
             */

            _ddPtr.Value += m * this._drainConductance;
            _ggPtr.Value += m * (gcggb - this._gtg);
            _ssPtr.Value += m * this._sourceConductance;
            _bbPtr.Value += m * ((this._gbd + this._gbs
                                - gcbgb - gcbdb - gcbsb) - this._gbbs);
            _dpdpPtr.Value += m * ((this._drainConductance
                                  + this._gds + this._gbd
                                  + RevSum + gcddb) + dxpart * this._gtd + gbdpdp);
            _spspPtr.Value += m * ((this._sourceConductance
                                  + this._gds + this._gbs
                                  + FwdSum + gcssb) + sxpart * this._gts + gbspsp);
            _ddpPtr.Value -= m * this._drainConductance;
            _gbPtr.Value -= m * (gcggb + gcgdb + gcgsb + this._gtb);
            _gdpPtr.Value += m * (gcgdb - this._gtd);
            _gspPtr.Value += m * (gcgsb - this._gts);
            _sspPtr.Value -= m * this._sourceConductance;
            _bgPtr.Value += m * (gcbgb - this._gbgs);
            _bdpPtr.Value += m * (gcbdb - this._gbd + gbbdp);
            _bspPtr.Value += m * (gcbsb - this._gbs + gbbsp);
            _dpdPtr.Value -= m * this._drainConductance;
            _dpgPtr.Value += m * ((Gm + gcdgb) + dxpart * this._gtg + gbdpg);
            _dpbPtr.Value -= m * ((this._gbd - Gmbs + gcdgb + gcddb
                                 + gcdsb - dxpart * this._gtb) - gbdpb);
            _dpspPtr.Value -= m * ((this._gds + FwdSum - gcdsb
                  - 0.5 * this._gts) - gbdpsp);
            _spgPtr.Value += m * (gcsgb - Gm + sxpart * this._gtg + gbspg);
            _spsPtr.Value -= m * this._sourceConductance;
            _spbPtr.Value -= m * ((this._gbs + Gmbs + gcsgb + gcsdb
                                 + gcssb - sxpart * this._gtg) - gbspb);
            _spdpPtr.Value -= m * ((this._gds + RevSum - gcsdb
                                  - sxpart * this._gtd - this._gbd) - gbspdp);

            _qqPtr.Value += m * (gqdef + this._gtau);

            _dpqPtr.Value += m * (dxpart * this._gtau);
            _spqPtr.Value += m * (sxpart * this._gtau);
            _gqPtr.Value -= m * this._gtau;

            _qgPtr.Value += m * (-gcqgb + this._gtg);
            _qdpPtr.Value += m * (-gcqdb + this._gtd);
            _qspPtr.Value += m * (-gcqsb + this._gts);
            _qbPtr.Value += m * (-gcqbb + this._gtb);
        }
    }
}
