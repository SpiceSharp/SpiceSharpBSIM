using System;
using SpiceSharp.Simulations;
using SpiceSharp.Behaviors;
using SpiceSharp.Components.Semiconductors;
using SpiceSharp;
using SpiceSharp.Components;
using SpiceSharp.Algebra;
using SpiceSharp.Components.Mosfets;
using SpiceSharp.Attributes;

namespace SpiceSharpBSIM.Components.Semiconductors.BSIM.BSIM3Behaviors
{
    /// <summary>
    /// Load behavior for a <see cref="BSIM3" />
    /// </summary>
    [BehaviorFor(typeof(BSIM3)), AddBehaviorIfNo(typeof(IBiasingBehavior))]
    [GeneratedParameters]
    public partial class BiasingBehavior : TemperatureBehavior, IBiasingBehavior, ITimeBehavior
    {
        public const double ScalingFactor = 1.0e-9;
        public const double DELTA_1 = 0.02;
        public const double DELTA_3 = 0.02;
        public const double DELTA_4 = 0.02;

        private readonly ITemperatureSimulationState _temperature;
        private readonly IIntegrationMethod _method;
        private readonly IIterationSimulationState _iteration;
        private readonly ITimeSimulationState _time;
        private readonly IDerivative _qb, _qg, _qd, _qcdump, _qcheq;

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

        protected double _gbs, _cbs, _gbd, _cbd, _mode,
            _thetavth, _von, _vgsteff, _rds, _abulk, _ueff,
            _abovVgst2Vtm, _vdsat, _vdseff, _gds, _gm, _gmbs,
            _gbbs, _gbgs, _gbds, _csub, _cggb, _cgsb, _cgdb, _cdgb,
            _cdsb, _cddb, _cbgb, _cbsb, _cbdb, _cqdb, _cqsb, _cqgb,
            _cqbb, _gtau, _qinv, _qgate, _qbulk, _qdrn, _cd, _capbs,
            _capbd, _taunet, _gtg, _gtd, _gts, _gtb, _vbd, _vbs, _vgs,
            _vds, _qdef, _qbs, _qbd;
        private readonly IVariable<double> _drain, _gate, _source, _bulk, _drainPrime, _sourcePrime, _q;
        private readonly Element<double> _dpPtr, _spPtr, _gPtr, _bPtr, _qPtr;
        private readonly Element<double> _ddPtr, _ggPtr, _ssPtr, _bbPtr, _dpdpPtr, _spspPtr, _ddpPtr,
            _gbPtr, _gdpPtr, _gspPtr, _sspPtr, _bdpPtr, _bspPtr, _dpspPtr, _dpdPtr, _bgPtr, _dpgPtr,
            _spgPtr, _spsPtr, _dpbPtr, _spbPtr, _spdpPtr, _qqPtr, _qdpPtr, _qspPtr, _qgPtr, _qbPtr,
            _dpqPtr, _spqPtr, _gqPtr; //, _bqPtr;

        /// <summary>
        /// Constructor
        /// </summary>
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
            if (!DrainConductance.Equals(0.0))
                _drainPrime = state.CreatePrivateVariable(Name.Combine("drain"), Units.Volt);
            else
                _drainPrime = _drain;
            if (!SourceConductance.Equals(0.0))
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
            double ag0, qgd, qgs, qgb, von, VgstNVt, ExpVgst;
            double cdrain, cdreq, ceqbd, ceqbs, ceqqb, ceqqd, ceqqg;
            double czbd, czbdsw, czbdswg, czbs, czbssw, czbsswg, evbd, evbs, arg, sarg;
            double Vfbeff, dVfbeff_dVg, dVfbeff_dVb, V3, V4;
            double gcbdb, gcbgb, gcbsb, gcddb, gcdgb, gcdsb, gcgdb, gcggb, gcgsb, gcsdb;
            double gcsgb, gcssb, MJ, MJSW, MJSWG;
            double vbd, vbs, vds, vgb, vgd, vgs, vgdo;
            double qgate = 0.0, qbulk = 0.0, qdrn = 0.0, qsrc, qinoi, cqgate, cqbulk, cqdrn;
            double Vds, Vgs, Vbs, Gmbs, FwdSum, RevSum;
            double Vgs_eff, Vfb;
            double Phis, dPhis_dVb, sqrtPhis, dsqrtPhis_dVb, Vth, dVth_dVb, dVth_dVd;
            double Vgst, dVgst_dVg, dVgst_dVb, dVgs_eff_dVg, Nvtm;
            double Vtm;
            double n, dn_dVb, dn_dVd, voffcv, noff, dnoff_dVd, dnoff_dVb;
            double ExpArg, V0, CoxWLcen, QovCox, LINK;
            double DeltaPhi, dDeltaPhi_dVg, VgDP, dVgDP_dVg;
            double Cox, Tox, Tcen, dTcen_dVg, dTcen_dVd, dTcen_dVb;
            double Ccen, Coxeff, dCoxeff_dVg, dCoxeff_dVd, dCoxeff_dVb;
            double Denomi, dDenomi_dVg, dDenomi_dVd, dDenomi_dVb;
            double ueff, dueff_dVg, dueff_dVd, dueff_dVb;
            double Esat, Vdsat;
            double EsatL, dEsatL_dVg, dEsatL_dVd, dEsatL_dVb;

            double dVdsat_dVg, dVdsat_dVb, dVdsat_dVd, Vasat, dAlphaz_dVg, dAlphaz_dVb;
            double dVasat_dVg, dVasat_dVb, dVasat_dVd, Va, dVa_dVd, dVa_dVg, dVa_dVb;
            double Vbseff, dVbseff_dVb, VbseffCV, dVbseffCV_dVb;
            double Arg1, One_Third_CoxWL, Two_Third_CoxWL, Alphaz, CoxWL;

            double T0, dT0_dVg, dT0_dVd, dT0_dVb;
            double T1, dT1_dVg, dT1_dVd, dT1_dVb;
            double T2, dT2_dVg, dT2_dVd, dT2_dVb;
            double T3, dT3_dVg, dT3_dVd, dT3_dVb;
            double T4;
            double T5;
            double T6;
            double T7;
            double T8;
            double T9;
            double T10;
            double T11, T12;
            double tmp, Abulk, dAbulk_dVb, Abulk0, dAbulk0_dVb;

            double VACLM, dVACLM_dVg, dVACLM_dVd, dVACLM_dVb;
            double VADIBL, dVADIBL_dVg, dVADIBL_dVd, dVADIBL_dVb;

            double Xdep, dXdep_dVb, lt1, dlt1_dVb, ltw, dltw_dVb, Delt_vth, dDelt_vth_dVb;
            double Theta0, dTheta0_dVb;
            double TempRatio, tmp1, tmp2, tmp3, tmp4;
            double DIBL_Sft, dDIBL_Sft_dVd, Lambda, dLambda_dVg;
            double a1;

            double Vgsteff, dVgsteff_dVg, dVgsteff_dVd, dVgsteff_dVb;
            double Vdseff, dVdseff_dVg, dVdseff_dVd, dVdseff_dVb;
            double VdseffCV, dVdseffCV_dVg, dVdseffCV_dVd, dVdseffCV_dVb;
            double diffVds, dAbulk_dVg;
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

            double qcheq = 0.0, qdef, gqdef = 0.0, cqdef, cqcheq, gtau_diff, gtau_drift;
            double gcqdb = 0.0, gcqsb = 0.0, gcqgb = 0.0, gcqbb = 0.0;
            double dxpart, sxpart, ggtg, ggtd, ggts, ggtb;
            double ddxpart_dVd, ddxpart_dVg, ddxpart_dVb, ddxpart_dVs;
            double dsxpart_dVd, dsxpart_dVg, dsxpart_dVb, dsxpart_dVs;

            double gbspsp, gbbdp, gbbsp, gbspg, gbspb, gbspdp;
            double gbdpdp, gbdpg, gbdpb, gbdpsp;
            double Cgg, Cgd, Cgb, Cdg, Cdd, Cds;
            double Csg, Csd, Css, Csb, Cbg, Cbd, Cbb;
            double Cgg1, Cgb1, Cgd1, Cbg1, Cbb1, Cbd1, Qac0, Qsub0;
            double dQac0_dVg, dQac0_dVb, dQsub0_dVg, dQsub0_dVd, dQsub0_dVb;

            double m;

            bool Check, ChargeComputationNeeded;
            /* double junk[50]; */


            ChargeComputationNeeded = _time != null || InitializeSmallSignal || InitializeTransient;

            Check = true;
            if (InitializeSmallSignal || InitializeTransient)
            {
                vbs = this._vbs;
                vgs = this._vgs;
                vds = this._vds;
                qdef = this._qdef;
            }
            else if (_iteration.Mode == IterationModes.Junction && !Parameters.Off)
            {
                vds = ModelParameters.B3Type * Parameters.IcVDS;
                vgs = ModelParameters.B3Type * Parameters.IcVGS;
                vbs = ModelParameters.B3Type * Parameters.IcVBS;
                qdef = 0.0;

                if (vds.Equals(0.0) && vgs.Equals(0.0) && vbs.Equals(0.0))
                {
                    vbs = 0.0;
                    vgs = ModelParameters.B3Type * this.Vth0 + 0.1;
                    vds = 0.1;
                }
            }
            else if ((_iteration.Mode == IterationModes.Junction || _iteration.Mode == IterationModes.Fix) && Parameters.Off)
            {
                qdef = vbs = vgs = vds = 0.0;
            }
            else
            {
                vbs = ModelParameters.B3Type * (_bulk.Value - _sourcePrime.Value);
                vgs = ModelParameters.B3Type * (_gate.Value - _sourcePrime.Value);
                vds = ModelParameters.B3Type * (_drainPrime.Value - _sourcePrime.Value);
                qdef = ModelParameters.B3Type * _q.Value;

                vbd = vbs - vds;
                vgd = vgs - vds;
                vgdo = this._vgs - this._vds;
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

            /* Source/drain junction diode DC model begins */
            Nvtm = ModelTemperature.Vtm * ModelParameters.JctEmissionCoeff;
            /* acm model */
            if (ModelParameters.AcmMod.Value == 0)
            {
                if ((Parameters.SourceArea <= 0.0) && (Parameters.SourcePerimeter <= 0.0))
                {
                    SourceSatCurrent = 1.0e-14;
                }
                else
                {
                    SourceSatCurrent = Parameters.SourceArea
                                     * ModelTemperature.JctTempSatCurDensity
                                     + Parameters.SourcePerimeter
                                     * ModelTemperature.JctSidewallTempSatCurDensity;
                }
                if ((Parameters.DrainArea <= 0.0) && (Parameters.DrainPerimeter <= 0.0))
                {
                    DrainSatCurrent = 1.0e-14;
                }
                else
                {
                    DrainSatCurrent = Parameters.DrainArea
                                    * ModelTemperature.JctTempSatCurDensity
                                    + Parameters.DrainPerimeter
                                    * ModelTemperature.JctSidewallTempSatCurDensity;
                }
            }
            else
            {
                ACM.SaturationCurrents(
                    ModelParameters.AcmMod,
                    ModelParameters.Calcacm,
                    Parameters.Geo,
                    ModelParameters.Hdif,
                    ModelParameters.Wmlt,
                    Parameters.Width,
                    ModelParameters.Xw,
                    ModelTemperature.JctTempSatCurDensity,
                    ModelTemperature.JctSidewallTempSatCurDensity,
                    Parameters.DrainArea,
                    Parameters.DrainPerimeter,
                    Parameters.SourceArea,
                    Parameters.SourcePerimeter,
                    out DrainSatCurrent,
                    out SourceSatCurrent
                );
            }

            if (SourceSatCurrent <= 0.0)
            {
                this._gbs = _iteration.Gmin;
                this._cbs = this._gbs * vbs;
            }
            else
            {
                if (ModelParameters.Ijth == 0.0)
                {
                    evbs = Math.Exp(vbs / Nvtm);
                    this._gbs = SourceSatCurrent * evbs / Nvtm + _iteration.Gmin;
                    this._cbs = SourceSatCurrent * (evbs - 1.0)
                                   + _iteration.Gmin * vbs;
                }
                else
                {
                    if (vbs < this.Vjsm)
                    {
                        evbs = Math.Exp(vbs / Nvtm);
                        this._gbs = SourceSatCurrent * evbs / Nvtm + _iteration.Gmin;
                        this._cbs = SourceSatCurrent * (evbs - 1.0)
                                       + _iteration.Gmin * vbs;
                    }
                    else
                    {
                        T0 = this.IsEvjsm / Nvtm;
                        this._gbs = T0 + _iteration.Gmin;
                        this._cbs = this.IsEvjsm - SourceSatCurrent
                                       + T0 * (vbs - this.Vjsm)
                                       + _iteration.Gmin * vbs;
                    }
                }
            }

            if (DrainSatCurrent <= 0.0)
            {
                this._gbd = _iteration.Gmin;
                this._cbd = this._gbd * vbd;
            }
            else
            {
                if (ModelParameters.Ijth == 0.0)
                {
                    evbd = Math.Exp(vbd / Nvtm);
                    this._gbd = DrainSatCurrent * evbd / Nvtm + _iteration.Gmin;
                    this._cbd = DrainSatCurrent * (evbd - 1.0)
                                   + _iteration.Gmin * vbd;
                }
                else
                {
                    if (vbd < this.Vjdm)
                    {
                        evbd = Math.Exp(vbd / Nvtm);
                        this._gbd = DrainSatCurrent * evbd / Nvtm + _iteration.Gmin;
                        this._cbd = DrainSatCurrent * (evbd - 1.0)
                                       + _iteration.Gmin * vbd;
                    }
                    else
                    {
                        T0 = this.IsEvjdm / Nvtm;
                        this._gbd = T0 + _iteration.Gmin;
                        this._cbd = this.IsEvjdm - DrainSatCurrent
                                       + T0 * (vbd - this.Vjdm)
                                       + _iteration.Gmin * vbd;
                    }
                }
            }
            /* End of diode DC model */

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

            T0 = Vbs - Param.BSIM3vbsc - 0.001;
            T1 = Math.Sqrt(T0 * T0 - 0.004 * Param.BSIM3vbsc);
            Vbseff = Param.BSIM3vbsc + 0.5 * (T0 + T1);
            dVbseff_dVb = 0.5 * (1.0 + T0 / T1);
            if (Vbseff < Vbs)
            {
                Vbseff = Vbs;
            }

            if (Vbseff > 0.0)
            {
                T0 = Param.BSIM3phi / (Param.BSIM3phi + Vbseff);
                Phis = Param.BSIM3phi * T0;
                dPhis_dVb = -T0 * T0;
                sqrtPhis = Param.BSIM3phis3 / (Param.BSIM3phi + 0.5 * Vbseff);
                dsqrtPhis_dVb = -0.5 * sqrtPhis * sqrtPhis / Param.BSIM3phis3;
            }
            else
            {
                Phis = Param.BSIM3phi - Vbseff;
                dPhis_dVb = -1.0;
                sqrtPhis = Math.Sqrt(Phis);
                dsqrtPhis_dVb = -0.5 / sqrtPhis;
            }
            Xdep = Param.BSIM3Xdep0 * sqrtPhis / Param.BSIM3sqrtPhi;
            dXdep_dVb = (Param.BSIM3Xdep0 / Param.BSIM3sqrtPhi)
                      * dsqrtPhis_dVb;

            Leff = Param.BSIM3leff;
            Vtm = ModelTemperature.Vtm;
            /* Vth Calculation */
            T3 = Math.Sqrt(Xdep);
            V0 = Param.BSIM3vbi - Param.BSIM3phi;

            T0 = Param.BSIM3dvt2 * Vbseff;
            if (T0 >= -0.5)
            {
                T1 = 1.0 + T0;
                T2 = Param.BSIM3dvt2;
            }
            else /* Added to avoid any discontinuity problems caused by dvt2 */
            {
                T4 = 1.0 / (3.0 + 8.0 * T0);
                T1 = (1.0 + 3.0 * T0) * T4;
                T2 = Param.BSIM3dvt2 * T4 * T4;
            }
            lt1 = ModelTemperature.Factor1 * T3 * T1;
            dlt1_dVb = ModelTemperature.Factor1 * (0.5 / T3 * T1 * dXdep_dVb + T3 * T2);

            T0 = Param.BSIM3dvt2w * Vbseff;
            if (T0 >= -0.5)
            {
                T1 = 1.0 + T0;
                T2 = Param.BSIM3dvt2w;
            }
            else /* Added to avoid any discontinuity problems caused by dvt2w */
            {
                T4 = 1.0 / (3.0 + 8.0 * T0);
                T1 = (1.0 + 3.0 * T0) * T4;
                T2 = Param.BSIM3dvt2w * T4 * T4;
            }
            ltw = ModelTemperature.Factor1 * T3 * T1;
            dltw_dVb = ModelTemperature.Factor1 * (0.5 / T3 * T1 * dXdep_dVb + T3 * T2);

            T0 = -0.5 * Param.BSIM3dvt1 * Leff / lt1;
            if (T0 > -EXP_THRESHOLD)
            {
                T1 = Math.Exp(T0);
                Theta0 = T1 * (1.0 + 2.0 * T1);
                dT1_dVb = -T0 / lt1 * T1 * dlt1_dVb;
                dTheta0_dVb = (1.0 + 4.0 * T1) * dT1_dVb;
            }
            else
            {
                T1 = MIN_EXP;
                Theta0 = T1 * (1.0 + 2.0 * T1);
                dTheta0_dVb = 0.0;
            }

            this._thetavth = Param.BSIM3dvt0 * Theta0;
            Delt_vth = this._thetavth * V0;
            dDelt_vth_dVb = Param.BSIM3dvt0 * dTheta0_dVb * V0;

            T0 = -0.5 * Param.BSIM3dvt1w * Param.BSIM3weff * Leff / ltw;
            if (T0 > -EXP_THRESHOLD)
            {
                T1 = Math.Exp(T0);
                T2 = T1 * (1.0 + 2.0 * T1);
                dT1_dVb = -T0 / ltw * T1 * dltw_dVb;
                dT2_dVb = (1.0 + 4.0 * T1) * dT1_dVb;
            }
            else
            {
                T1 = MIN_EXP;
                T2 = T1 * (1.0 + 2.0 * T1);
                dT2_dVb = 0.0;
            }

            T0 = Param.BSIM3dvt0w * T2;
            T2 = T0 * V0;
            dT2_dVb = Param.BSIM3dvt0w * dT2_dVb * V0;

            TempRatio = _temperature.Temperature / ModelParameters.Tnom - 1.0;
            T0 = Math.Sqrt(1.0 + Param.BSIM3nlx / Leff);
            T1 = Param.BSIM3k1ox * (T0 - 1.0) * Param.BSIM3sqrtPhi
               + (Param.BSIM3kt1 + Param.BSIM3kt1l / Leff
               + Param.BSIM3kt2 * Vbseff) * TempRatio;
            tmp2 = ModelParameters.Tox * Param.BSIM3phi
                 / (Param.BSIM3weff + Param.BSIM3w0);

            T3 = Param.BSIM3eta0 + Param.BSIM3etab * Vbseff;
            if (T3 < 1.0e-4) /* avoid  discontinuity problems caused by etab */
            {
                T9 = 1.0 / (3.0 - 2.0e4 * T3);
                T3 = (2.0e-4 - T3) * T9;
                T4 = T9 * T9;
            }
            else
            {
                T4 = 1.0;
            }
            dDIBL_Sft_dVd = T3 * Param.BSIM3theta0vb0;
            DIBL_Sft = dDIBL_Sft_dVd * Vds;

            Vth = ModelParameters.B3Type * this.Vth0 - Param.BSIM3k1
                * Param.BSIM3sqrtPhi + Param.BSIM3k1ox * sqrtPhis
                - Param.BSIM3k2ox * Vbseff - Delt_vth - T2 + (Param.BSIM3k3
                + Param.BSIM3k3b * Vbseff) * tmp2 + T1 - DIBL_Sft;

            this._von = Vth;

            dVth_dVb = Param.BSIM3k1ox * dsqrtPhis_dVb - Param.BSIM3k2ox
                     - dDelt_vth_dVb - dT2_dVb + Param.BSIM3k3b * tmp2
                     - Param.BSIM3etab * Vds * Param.BSIM3theta0vb0 * T4
                     + Param.BSIM3kt2 * TempRatio;
            dVth_dVd = -dDIBL_Sft_dVd;

            /* Calculate n */
            tmp2 = Param.BSIM3nfactor * EPSSI / Xdep;
            tmp3 = Param.BSIM3cdsc + Param.BSIM3cdscb * Vbseff
                 + Param.BSIM3cdscd * Vds;
            tmp4 = (tmp2 + tmp3 * Theta0 + Param.BSIM3cit) / ModelTemperature.Cox;
            if (tmp4 >= -0.5)
            {
                n = 1.0 + tmp4;
                dn_dVb = (-tmp2 / Xdep * dXdep_dVb + tmp3 * dTheta0_dVb
                       + Param.BSIM3cdscb * Theta0) / ModelTemperature.Cox;
                dn_dVd = Param.BSIM3cdscd * Theta0 / ModelTemperature.Cox;
            }
            else
            /* avoid  discontinuity problems caused by tmp4 */
            {
                T0 = 1.0 / (3.0 + 8.0 * tmp4);
                n = (1.0 + 3.0 * tmp4) * T0;
                T0 *= T0;
                dn_dVb = (-tmp2 / Xdep * dXdep_dVb + tmp3 * dTheta0_dVb
                       + Param.BSIM3cdscb * Theta0) / ModelTemperature.Cox * T0;
                dn_dVd = Param.BSIM3cdscd * Theta0 / ModelTemperature.Cox * T0;
            }

            /* Poly Gate Si Depletion Effect */
            T0 = this.Vfb + Param.BSIM3phi;
            if ((Param.BSIM3ngate > 1.0e18) && (Param.BSIM3ngate < 1.0e25)
                 && (Vgs > T0))
            /* added to avoid the problem caused by ngate */
            {
                T1 = 1.0e6 * Charge_q * EPSSI * Param.BSIM3ngate
                   / (ModelTemperature.Cox * ModelTemperature.Cox);
                T4 = Math.Sqrt(1.0 + 2.0 * (Vgs - T0) / T1);
                T2 = T1 * (T4 - 1.0);
                T3 = 0.5 * T2 * T2 / T1; /* T3 = Vpoly */
                T7 = 1.12 - T3 - 0.05;
                T6 = Math.Sqrt(T7 * T7 + 0.224);
                T5 = 1.12 - 0.5 * (T7 + T6);
                Vgs_eff = Vgs - T5;
                dVgs_eff_dVg = 1.0 - (0.5 - 0.5 / T4) * (1.0 + T7 / T6);
            }
            else
            {
                Vgs_eff = Vgs;
                dVgs_eff_dVg = 1.0;
            }
            Vgst = Vgs_eff - Vth;

            /* Effective Vgst (Vgsteff) Calculation */

            T10 = 2.0 * n * Vtm;
            VgstNVt = Vgst / T10;
            ExpArg = (2.0 * Param.BSIM3voff - Vgst) / T10;

            /* MCJ: Very small Vgst */
            if (VgstNVt > EXP_THRESHOLD)
            {
                Vgsteff = Vgst;
                dVgsteff_dVg = dVgs_eff_dVg;
                dVgsteff_dVd = -dVth_dVd;
                dVgsteff_dVb = -dVth_dVb;
            }
            else if (ExpArg > EXP_THRESHOLD)
            {
                T0 = (Vgst - Param.BSIM3voff) / (n * Vtm);
                ExpVgst = Math.Exp(T0);
                Vgsteff = Vtm * Param.BSIM3cdep0 / ModelTemperature.Cox * ExpVgst;
                dVgsteff_dVg = Vgsteff / (n * Vtm);
                dVgsteff_dVd = -dVgsteff_dVg * (dVth_dVd + T0 * Vtm * dn_dVd);
                dVgsteff_dVb = -dVgsteff_dVg * (dVth_dVb + T0 * Vtm * dn_dVb);
                dVgsteff_dVg *= dVgs_eff_dVg;
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

                dT2_dVg = -ModelTemperature.Cox / (Vtm * Param.BSIM3cdep0)
                        * Math.Exp(ExpArg);
                T2 = 1.0 - T10 * dT2_dVg;
                dT2_dVd = -dT2_dVg * (dVth_dVd - 2.0 * Vtm * ExpArg * dn_dVd)
                        + (T2 - 1.0) / n * dn_dVd;
                dT2_dVb = -dT2_dVg * (dVth_dVb - 2.0 * Vtm * ExpArg * dn_dVb)
                        + (T2 - 1.0) / n * dn_dVb;

                Vgsteff = T1 / T2;
                T3 = T2 * T2;
                dVgsteff_dVg = (T2 * dT1_dVg - T1 * dT2_dVg) / T3 * dVgs_eff_dVg;
                dVgsteff_dVd = (T2 * dT1_dVd - T1 * dT2_dVd) / T3;
                dVgsteff_dVb = (T2 * dT1_dVb - T1 * dT2_dVb) / T3;
            }
            this._vgsteff = Vgsteff;

            /* Calculate Effective Channel Geometry */
            T9 = sqrtPhis - Param.BSIM3sqrtPhi;
            Weff = Param.BSIM3weff - 2.0 * (Param.BSIM3dwg * Vgsteff
                 + Param.BSIM3dwb * T9);
            dWeff_dVg = -2.0 * Param.BSIM3dwg;
            dWeff_dVb = -2.0 * Param.BSIM3dwb * dsqrtPhis_dVb;

            if (Weff < 2.0e-8) /* to avoid the discontinuity problem due to Weff*/
            {
                T0 = 1.0 / (6.0e-8 - 2.0 * Weff);
                Weff = 2.0e-8 * (4.0e-8 - Weff) * T0;
                T0 *= T0 * 4.0e-16;
                dWeff_dVg *= T0;
                dWeff_dVb *= T0;
            }

            T0 = Param.BSIM3prwg * Vgsteff + Param.BSIM3prwb * T9;
            if (T0 >= -0.9)
            {
                Rds = Param.BSIM3rds0 * (1.0 + T0);
                dRds_dVg = Param.BSIM3rds0 * Param.BSIM3prwg;
                dRds_dVb = Param.BSIM3rds0 * Param.BSIM3prwb * dsqrtPhis_dVb;
            }
            else
            /* to avoid the discontinuity problem due to prwg and prwb*/
            {
                T1 = 1.0 / (17.0 + 20.0 * T0);
                Rds = Param.BSIM3rds0 * (0.8 + T0) * T1;
                T1 *= T1;
                dRds_dVg = Param.BSIM3rds0 * Param.BSIM3prwg * T1;
                dRds_dVb = Param.BSIM3rds0 * Param.BSIM3prwb * dsqrtPhis_dVb
                         * T1;
            }
            this._rds = Rds; /* Noise Bugfix */

            /* Calculate Abulk */
            T1 = 0.5 * Param.BSIM3k1ox / sqrtPhis;
            dT1_dVb = -T1 / sqrtPhis * dsqrtPhis_dVb;

            T9 = Math.Sqrt(Param.BSIM3xj * Xdep);
            tmp1 = Leff + 2.0 * T9;
            T5 = Leff / tmp1;
            tmp2 = Param.BSIM3a0 * T5;
            tmp3 = Param.BSIM3weff + Param.BSIM3b1;
            tmp4 = Param.BSIM3b0 / tmp3;
            T2 = tmp2 + tmp4;
            dT2_dVb = -T9 / tmp1 / Xdep * dXdep_dVb;
            T6 = T5 * T5;
            T7 = T5 * T6;

            Abulk0 = 1.0 + T1 * T2;
            dAbulk0_dVb = T1 * tmp2 * dT2_dVb + T2 * dT1_dVb;

            T8 = Param.BSIM3ags * Param.BSIM3a0 * T7;
            dAbulk_dVg = -T1 * T8;
            Abulk = Abulk0 + dAbulk_dVg * Vgsteff;
            dAbulk_dVb = dAbulk0_dVb - T8 * Vgsteff * (dT1_dVb
                       + 3.0 * T1 * dT2_dVb);

            if (Abulk0 < 0.1) /* added to avoid the problems caused by Abulk0 */
            {
                T9 = 1.0 / (3.0 - 20.0 * Abulk0);
                Abulk0 = (0.2 - Abulk0) * T9;
                dAbulk0_dVb *= T9 * T9;
            }

            if (Abulk < 0.1)
            /* added to avoid the problems caused by Abulk */
            {
                T9 = 1.0 / (3.0 - 20.0 * Abulk);
                Abulk = (0.2 - Abulk) * T9;
                T10 = T9 * T9;
                dAbulk_dVb *= T10;
                dAbulk_dVg *= T10;
            }
            this._abulk = Abulk;

            T2 = Param.BSIM3keta * Vbseff;
            if (T2 >= -0.9)
            {
                T0 = 1.0 / (1.0 + T2);
                dT0_dVb = -Param.BSIM3keta * T0 * T0;
            }
            else
            /* added to avoid the problems caused by Keta */
            {
                T1 = 1.0 / (0.8 + T2);
                T0 = (17.0 + 20.0 * T2) * T1;
                dT0_dVb = -Param.BSIM3keta * T1 * T1;
            }
            dAbulk_dVg *= T0;
            dAbulk_dVb = dAbulk_dVb * T0 + Abulk * dT0_dVb;
            dAbulk0_dVb = dAbulk0_dVb * T0 + Abulk0 * dT0_dVb;
            Abulk *= T0;
            Abulk0 *= T0;


            /* Mobility calculation */
            if (ModelParameters.MobMod.Value == 1)
            {
                T0 = Vgsteff + Vth + Vth;
                T2 = Param.BSIM3ua + Param.BSIM3uc * Vbseff;
                T3 = T0 / ModelParameters.Tox;
                T5 = T3 * (T2 + Param.BSIM3ub * T3);
                dDenomi_dVg = (T2 + 2.0 * Param.BSIM3ub * T3) / ModelParameters.Tox;
                dDenomi_dVd = dDenomi_dVg * 2.0 * dVth_dVd;
                dDenomi_dVb = dDenomi_dVg * 2.0 * dVth_dVb + Param.BSIM3uc * T3;
            }
            else if (ModelParameters.MobMod.Value == 2)
            {
                T5 = Vgsteff / ModelParameters.Tox * (Param.BSIM3ua
                   + Param.BSIM3uc * Vbseff + Param.BSIM3ub * Vgsteff
                   / ModelParameters.Tox);
                dDenomi_dVg = (Param.BSIM3ua + Param.BSIM3uc * Vbseff
                            + 2.0 * Param.BSIM3ub * Vgsteff / ModelParameters.Tox)
                            / ModelParameters.Tox;
                dDenomi_dVd = 0.0;
                dDenomi_dVb = Vgsteff * Param.BSIM3uc / ModelParameters.Tox;
            }
            else
            {
                T0 = Vgsteff + Vth + Vth;
                T2 = 1.0 + Param.BSIM3uc * Vbseff;
                T3 = T0 / ModelParameters.Tox;
                T4 = T3 * (Param.BSIM3ua + Param.BSIM3ub * T3);
                T5 = T4 * T2;
                dDenomi_dVg = (Param.BSIM3ua + 2.0 * Param.BSIM3ub * T3) * T2
                            / ModelParameters.Tox;
                dDenomi_dVd = dDenomi_dVg * 2.0 * dVth_dVd;
                dDenomi_dVb = dDenomi_dVg * 2.0 * dVth_dVb + Param.BSIM3uc * T4;
            }

            if (T5 >= -0.8)
            {
                Denomi = 1.0 + T5;
            }
            else /* Added to avoid the discontinuity problem caused by ua and ub*/
            {
                T9 = 1.0 / (7.0 + 10.0 * T5);
                Denomi = (0.6 + T5) * T9;
                T9 *= T9;
                dDenomi_dVg *= T9;
                dDenomi_dVd *= T9;
                dDenomi_dVb *= T9;
            }

            this._ueff = ueff = this.U0temp / Denomi;
            T9 = -ueff / Denomi;
            dueff_dVg = T9 * dDenomi_dVg;
            dueff_dVd = T9 * dDenomi_dVd;
            dueff_dVb = T9 * dDenomi_dVb;

            /* Saturation Drain Voltage  Vdsat */
            WVCox = Weff * Param.BSIM3vsattemp * ModelTemperature.Cox;
            WVCoxRds = WVCox * Rds;

            Esat = 2.0 * Param.BSIM3vsattemp / ueff;
            EsatL = Esat * Leff;
            T0 = -EsatL / ueff;
            dEsatL_dVg = T0 * dueff_dVg;
            dEsatL_dVd = T0 * dueff_dVd;
            dEsatL_dVb = T0 * dueff_dVb;

            /* Sqrt() */
            a1 = Param.BSIM3a1;
            if (a1 == 0.0)
            {
                Lambda = Param.BSIM3a2;
                dLambda_dVg = 0.0;
            }
            else if (a1 > 0.0)
            /* Added to avoid the discontinuity problem
               caused by a1 and a2 (Lambda) */
            {
                T0 = 1.0 - Param.BSIM3a2;
                T1 = T0 - Param.BSIM3a1 * Vgsteff - 0.0001;
                T2 = Math.Sqrt(T1 * T1 + 0.0004 * T0);
                Lambda = Param.BSIM3a2 + T0 - 0.5 * (T1 + T2);
                dLambda_dVg = 0.5 * Param.BSIM3a1 * (1.0 + T1 / T2);
            }
            else
            {
                T1 = Param.BSIM3a2 + Param.BSIM3a1 * Vgsteff - 0.0001;
                T2 = Math.Sqrt(T1 * T1 + 0.0004 * Param.BSIM3a2);
                Lambda = 0.5 * (T1 + T2);
                dLambda_dVg = 0.5 * Param.BSIM3a1 * (1.0 + T1 / T2);
            }

            Vgst2Vtm = Vgsteff + 2.0 * Vtm;
            this._abovVgst2Vtm = Abulk / Vgst2Vtm;

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
            if ((Rds == 0.0) && (Lambda == 1.0))
            {
                T0 = 1.0 / (Abulk * EsatL + Vgst2Vtm);
                tmp1 = 0.0;
                T1 = T0 * T0;
                T2 = Vgst2Vtm * T0;
                T3 = EsatL * Vgst2Vtm;
                Vdsat = T3 * T0;

                dT0_dVg = -(Abulk * dEsatL_dVg + EsatL * dAbulk_dVg + 1.0) * T1;
                dT0_dVd = -(Abulk * dEsatL_dVd) * T1;
                dT0_dVb = -(Abulk * dEsatL_dVb + dAbulk_dVb * EsatL) * T1;

                dVdsat_dVg = T3 * dT0_dVg + T2 * dEsatL_dVg + EsatL * T0;
                dVdsat_dVd = T3 * dT0_dVd + T2 * dEsatL_dVd;
                dVdsat_dVb = T3 * dT0_dVb + T2 * dEsatL_dVb;
            }
            else
            {
                tmp1 = dLambda_dVg / (Lambda * Lambda);
                T9 = Abulk * WVCoxRds;
                T8 = Abulk * T9;
                T7 = Vgst2Vtm * T9;
                T6 = Vgst2Vtm * WVCoxRds;
                T0 = 2.0 * Abulk * (T9 - 1.0 + 1.0 / Lambda);
                dT0_dVg = 2.0 * (T8 * tmp2 - Abulk * tmp1
                        + (2.0 * T9 + 1.0 / Lambda - 1.0) * dAbulk_dVg);

                dT0_dVb = 2.0 * (T8 * (2.0 / Abulk * dAbulk_dVb + tmp3)
                        + (1.0 / Lambda - 1.0) * dAbulk_dVb);
                dT0_dVd = 0.0;
                T1 = Vgst2Vtm * (2.0 / Lambda - 1.0) + Abulk * EsatL + 3.0 * T7;

                dT1_dVg = (2.0 / Lambda - 1.0) - 2.0 * Vgst2Vtm * tmp1
                        + Abulk * dEsatL_dVg + EsatL * dAbulk_dVg + 3.0 * (T9
                        + T7 * tmp2 + T6 * dAbulk_dVg);
                dT1_dVb = Abulk * dEsatL_dVb + EsatL * dAbulk_dVb
                        + 3.0 * (T6 * dAbulk_dVb + T7 * tmp3);
                dT1_dVd = Abulk * dEsatL_dVd;

                T2 = Vgst2Vtm * (EsatL + 2.0 * T6);
                dT2_dVg = EsatL + Vgst2Vtm * dEsatL_dVg
                        + T6 * (4.0 + 2.0 * Vgst2Vtm * tmp2);
                dT2_dVb = Vgst2Vtm * (dEsatL_dVb + 2.0 * T6 * tmp3);
                dT2_dVd = Vgst2Vtm * dEsatL_dVd;

                T3 = Math.Sqrt(T1 * T1 - 2.0 * T0 * T2);
                Vdsat = (T1 - T3) / T0;

                dT3_dVg = (T1 * dT1_dVg - 2.0 * (T0 * dT2_dVg + T2 * dT0_dVg))
                        / T3;
                dT3_dVd = (T1 * dT1_dVd - 2.0 * (T0 * dT2_dVd + T2 * dT0_dVd))
                        / T3;
                dT3_dVb = (T1 * dT1_dVb - 2.0 * (T0 * dT2_dVb + T2 * dT0_dVb))
                        / T3;

                dVdsat_dVg = (dT1_dVg - (T1 * dT1_dVg - dT0_dVg * T2
                           - T0 * dT2_dVg) / T3 - Vdsat * dT0_dVg) / T0;
                dVdsat_dVb = (dT1_dVb - (T1 * dT1_dVb - dT0_dVb * T2
                           - T0 * dT2_dVb) / T3 - Vdsat * dT0_dVb) / T0;
                dVdsat_dVd = (dT1_dVd - (T1 * dT1_dVd - T0 * dT2_dVd) / T3) / T0;
            }
            this._vdsat = Vdsat;

            /* Effective Vds (Vdseff) Calculation */
            T1 = Vdsat - Vds - Param.BSIM3delta;
            dT1_dVg = dVdsat_dVg;
            dT1_dVd = dVdsat_dVd - 1.0;
            dT1_dVb = dVdsat_dVb;

            T2 = Math.Sqrt(T1 * T1 + 4.0 * Param.BSIM3delta * Vdsat);
            T0 = T1 / T2;
            T3 = 2.0 * Param.BSIM3delta / T2;
            dT2_dVg = T0 * dT1_dVg + T3 * dVdsat_dVg;
            dT2_dVd = T0 * dT1_dVd + T3 * dVdsat_dVd;
            dT2_dVb = T0 * dT1_dVb + T3 * dVdsat_dVb;

            Vdseff = Vdsat - 0.5 * (T1 + T2);
            dVdseff_dVg = dVdsat_dVg - 0.5 * (dT1_dVg + dT2_dVg);
            dVdseff_dVd = dVdsat_dVd - 0.5 * (dT1_dVd + dT2_dVd);
            dVdseff_dVb = dVdsat_dVb - 0.5 * (dT1_dVb + dT2_dVb);
            /* Added to eliminate non-zero Vdseff at Vds=0.0 */
            if (Vds == 0.0)
            {
                Vdseff = 0.0;
                dVdseff_dVg = 0.0;
                dVdseff_dVb = 0.0;
            }

            /* Calculate VAsat */
            tmp4 = 1.0 - 0.5 * Abulk * Vdsat / Vgst2Vtm;
            T9 = WVCoxRds * Vgsteff;
            T8 = T9 / Vgst2Vtm;
            T0 = EsatL + Vdsat + 2.0 * T9 * tmp4;

            T7 = 2.0 * WVCoxRds * tmp4;
            dT0_dVg = dEsatL_dVg + dVdsat_dVg + T7 * (1.0 + tmp2 * Vgsteff)
                    - T8 * (Abulk * dVdsat_dVg - Abulk * Vdsat / Vgst2Vtm
                    + Vdsat * dAbulk_dVg);

            dT0_dVb = dEsatL_dVb + dVdsat_dVb + T7 * tmp3 * Vgsteff
                    - T8 * (dAbulk_dVb * Vdsat + Abulk * dVdsat_dVb);
            dT0_dVd = dEsatL_dVd + dVdsat_dVd - T8 * Abulk * dVdsat_dVd;

            T9 = WVCoxRds * Abulk;
            T1 = 2.0 / Lambda - 1.0 + T9;
            dT1_dVg = -2.0 * tmp1 + WVCoxRds * (Abulk * tmp2 + dAbulk_dVg);
            dT1_dVb = dAbulk_dVb * WVCoxRds + T9 * tmp3;

            Vasat = T0 / T1;
            dVasat_dVg = (dT0_dVg - Vasat * dT1_dVg) / T1;
            dVasat_dVb = (dT0_dVb - Vasat * dT1_dVb) / T1;
            dVasat_dVd = dT0_dVd / T1;

            if (Vdseff > Vds)
                Vdseff = Vds;
            diffVds = Vds - Vdseff;
            this._vdseff = Vdseff;

            /* Calculate VACLM */
            if ((Param.BSIM3pclm > 0.0) && (diffVds > 1.0e-10))
            {
                T0 = 1.0 / (Param.BSIM3pclm * Abulk * Param.BSIM3litl);
                dT0_dVb = -T0 / Abulk * dAbulk_dVb;
                dT0_dVg = -T0 / Abulk * dAbulk_dVg;

                T2 = Vgsteff / EsatL;
                T1 = Leff * (Abulk + T2);
                dT1_dVg = Leff * ((1.0 - T2 * dEsatL_dVg) / EsatL + dAbulk_dVg);
                dT1_dVb = Leff * (dAbulk_dVb - T2 * dEsatL_dVb / EsatL);
                dT1_dVd = -T2 * dEsatL_dVd / Esat;

                T9 = T0 * T1;
                VACLM = T9 * diffVds;
                dVACLM_dVg = T0 * dT1_dVg * diffVds - T9 * dVdseff_dVg
                           + T1 * diffVds * dT0_dVg;
                dVACLM_dVb = (dT0_dVb * T1 + T0 * dT1_dVb) * diffVds
                           - T9 * dVdseff_dVb;
                dVACLM_dVd = T0 * dT1_dVd * diffVds + T9 * (1.0 - dVdseff_dVd);
            }
            else
            {
                VACLM = MAX_EXP;
                dVACLM_dVd = dVACLM_dVg = dVACLM_dVb = 0.0;
            }

            /* Calculate VADIBL */
            if (Param.BSIM3thetaRout > 0.0)
            {
                T8 = Abulk * Vdsat;
                T0 = Vgst2Vtm * T8;
                dT0_dVg = Vgst2Vtm * Abulk * dVdsat_dVg + T8
                        + Vgst2Vtm * Vdsat * dAbulk_dVg;
                dT0_dVb = Vgst2Vtm * (dAbulk_dVb * Vdsat + Abulk * dVdsat_dVb);
                dT0_dVd = Vgst2Vtm * Abulk * dVdsat_dVd;

                T1 = Vgst2Vtm + T8;
                dT1_dVg = 1.0 + Abulk * dVdsat_dVg + Vdsat * dAbulk_dVg;
                dT1_dVb = Abulk * dVdsat_dVb + dAbulk_dVb * Vdsat;
                dT1_dVd = Abulk * dVdsat_dVd;

                T9 = T1 * T1;
                T2 = Param.BSIM3thetaRout;
                VADIBL = (Vgst2Vtm - T0 / T1) / T2;
                dVADIBL_dVg = (1.0 - dT0_dVg / T1 + T0 * dT1_dVg / T9) / T2;
                dVADIBL_dVb = (-dT0_dVb / T1 + T0 * dT1_dVb / T9) / T2;
                dVADIBL_dVd = (-dT0_dVd / T1 + T0 * dT1_dVd / T9) / T2;

                T7 = Param.BSIM3pdiblb * Vbseff;
                if (T7 >= -0.9)
                {
                    T3 = 1.0 / (1.0 + T7);
                    VADIBL *= T3;
                    dVADIBL_dVg *= T3;
                    dVADIBL_dVb = (dVADIBL_dVb - VADIBL * Param.BSIM3pdiblb)
                                * T3;
                    dVADIBL_dVd *= T3;
                }
                else
                /* Added to avoid the discontinuity problem caused by pdiblcb */
                {
                    T4 = 1.0 / (0.8 + T7);
                    T3 = (17.0 + 20.0 * T7) * T4;
                    dVADIBL_dVg *= T3;
                    dVADIBL_dVb = dVADIBL_dVb * T3
                                - VADIBL * Param.BSIM3pdiblb * T4 * T4;
                    dVADIBL_dVd *= T3;
                    VADIBL *= T3;
                }
            }
            else
            {
                VADIBL = MAX_EXP;
                dVADIBL_dVd = dVADIBL_dVg = dVADIBL_dVb = 0.0;
            }

            /* Calculate VA */

            T8 = Param.BSIM3pvag / EsatL;
            T9 = T8 * Vgsteff;
            if (T9 > -0.9)
            {
                T0 = 1.0 + T9;
                dT0_dVg = T8 * (1.0 - Vgsteff * dEsatL_dVg / EsatL);
                dT0_dVb = -T9 * dEsatL_dVb / EsatL;
                dT0_dVd = -T9 * dEsatL_dVd / EsatL;
            }
            else /* Added to avoid the discontinuity problems caused by pvag */
            {
                T1 = 1.0 / (17.0 + 20.0 * T9);
                T0 = (0.8 + T9) * T1;
                T1 *= T1;
                dT0_dVg = T8 * (1.0 - Vgsteff * dEsatL_dVg / EsatL) * T1;

                T9 *= T1 / EsatL;
                dT0_dVb = -T9 * dEsatL_dVb;
                dT0_dVd = -T9 * dEsatL_dVd;
            }

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
            if (Param.BSIM3pscbe2 > 0.0)
            {
                if (diffVds > Param.BSIM3pscbe1 * Param.BSIM3litl
                    / EXP_THRESHOLD)
                {
                    T0 = Param.BSIM3pscbe1 * Param.BSIM3litl / diffVds;
                    VASCBE = Leff * Math.Exp(T0) / Param.BSIM3pscbe2;
                    T1 = T0 * VASCBE / diffVds;
                    dVASCBE_dVg = T1 * dVdseff_dVg;
                    dVASCBE_dVd = -T1 * (1.0 - dVdseff_dVd);
                    dVASCBE_dVb = T1 * dVdseff_dVb;
                }
                else
                {
                    VASCBE = MAX_EXP * Leff / Param.BSIM3pscbe2;
                    dVASCBE_dVg = dVASCBE_dVd = dVASCBE_dVb = 0.0;
                }
            }
            else
            {
                VASCBE = MAX_EXP;
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

            T9 = Vdseff / EsatL;
            fgche2 = 1.0 + T9;
            dfgche2_dVg = (dVdseff_dVg - T9 * dEsatL_dVg) / EsatL;
            dfgche2_dVd = (dVdseff_dVd - T9 * dEsatL_dVd) / EsatL;
            dfgche2_dVb = (dVdseff_dVb - T9 * dEsatL_dVb) / EsatL;

            gche = beta * fgche1 / fgche2;
            dgche_dVg = (beta * dfgche1_dVg + fgche1 * dbeta_dVg
                      - gche * dfgche2_dVg) / fgche2;
            dgche_dVd = (beta * dfgche1_dVd + fgche1 * dbeta_dVd
                      - gche * dfgche2_dVd) / fgche2;
            dgche_dVb = (beta * dfgche1_dVb + fgche1 * dbeta_dVb
                      - gche * dfgche2_dVb) / fgche2;

            T0 = 1.0 + gche * Rds;
            T9 = Vdseff / T0;
            Idl = gche * T9;

            dIdl_dVg = (gche * dVdseff_dVg + T9 * dgche_dVg) / T0
                     - Idl * gche / T0 * dRds_dVg;

            dIdl_dVd = (gche * dVdseff_dVd + T9 * dgche_dVd) / T0;
            dIdl_dVb = (gche * dVdseff_dVb + T9 * dgche_dVb
                     - Idl * dRds_dVb * gche) / T0;

            T9 = diffVds / Va;
            T0 = 1.0 + T9;
            Idsa = Idl * T0;
            dIdsa_dVg = T0 * dIdl_dVg - Idl * (dVdseff_dVg + T9 * dVa_dVg) / Va;
            dIdsa_dVd = T0 * dIdl_dVd + Idl * (1.0 - dVdseff_dVd
                      - T9 * dVa_dVd) / Va;
            dIdsa_dVb = T0 * dIdl_dVb - Idl * (dVdseff_dVb + T9 * dVa_dVb) / Va;

            T9 = diffVds / VASCBE;
            T0 = 1.0 + T9;
            Ids = Idsa * T0;

            Gm = T0 * dIdsa_dVg - Idsa * (dVdseff_dVg + T9 * dVASCBE_dVg) / VASCBE;
            Gds = T0 * dIdsa_dVd + Idsa * (1.0 - dVdseff_dVd
                - T9 * dVASCBE_dVd) / VASCBE;
            Gmb = T0 * dIdsa_dVb - Idsa * (dVdseff_dVb
                + T9 * dVASCBE_dVb) / VASCBE;

            Gds += Gm * dVgsteff_dVd;
            Gmb += Gm * dVgsteff_dVb;
            Gm *= dVgsteff_dVg;
            Gmb *= dVbseff_dVb;

            /* Substrate current begins */
            tmp = Param.BSIM3alpha0 + Param.BSIM3alpha1 * Leff;
            if ((tmp <= 0.0) || (Param.BSIM3beta0 <= 0.0))
            {
                Isub = Gbd = Gbb = Gbg = 0.0;
            }
            else
            {
                T2 = tmp / Leff;
                if (diffVds > Param.BSIM3beta0 / EXP_THRESHOLD)
                {
                    T0 = -Param.BSIM3beta0 / diffVds;
                    T1 = T2 * diffVds * Math.Exp(T0);
                    T3 = T1 / diffVds * (T0 - 1.0);
                    dT1_dVg = T3 * dVdseff_dVg;
                    dT1_dVd = T3 * (dVdseff_dVd - 1.0);
                    dT1_dVb = T3 * dVdseff_dVb;
                }
                else
                {
                    T3 = T2 * MIN_EXP;
                    T1 = T3 * diffVds;
                    dT1_dVg = -T3 * dVdseff_dVg;
                    dT1_dVd = T3 * (1.0 - dVdseff_dVd);
                    dT1_dVb = -T3 * dVdseff_dVb;
                }
                Isub = T1 * Idsa;
                Gbg = T1 * dIdsa_dVg + Idsa * dT1_dVg;
                Gbd = T1 * dIdsa_dVd + Idsa * dT1_dVd;
                Gbb = T1 * dIdsa_dVb + Idsa * dT1_dVb;

                Gbd += Gbg * dVgsteff_dVd;
                Gbb += Gbg * dVgsteff_dVb;
                Gbg *= dVgsteff_dVg;
                Gbb *= dVbseff_dVb; /* bug fixing */
            }

            cdrain = Ids;
            this._gds = Gds;
            this._gm = Gm;
            this._gmbs = Gmb;

            this._gbbs = Gbb;
            this._gbgs = Gbg;
            this._gbds = Gbd;

            this._csub = Isub;

            /* BSIM3 thermal noise Qinv calculated from all capMod
             * 0, 1, 2 & 3 stored in this.Qinv 1/1998 */

            if ((ModelParameters.Xpart.Value < 0) || (!ChargeComputationNeeded))
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
            else if (ModelParameters.CapMod.Value == 0)
            {
                if (Vbseff < 0.0)
                {
                    Vbseff = Vbs;
                    dVbseff_dVb = 1.0;
                }
                else
                {
                    Vbseff = Param.BSIM3phi - Phis;
                    dVbseff_dVb = -dPhis_dVb;
                }

                Vfb = Param.BSIM3vfbcv;
                Vth = Vfb + Param.BSIM3phi + Param.BSIM3k1ox * sqrtPhis;
                Vgst = Vgs_eff - Vth;
                dVth_dVb = Param.BSIM3k1ox * dsqrtPhis_dVb;
                dVgst_dVb = -dVth_dVb;
                dVgst_dVg = dVgs_eff_dVg;

                CoxWL = ModelTemperature.Cox * Param.BSIM3weffCV
                      * Param.BSIM3leffCV;
                Arg1 = Vgs_eff - Vbseff - Vfb;

                if (Arg1 <= 0.0)
                {
                    qgate = CoxWL * Arg1;
                    qbulk = -qgate;
                    qdrn = 0.0;

                    this._cggb = CoxWL * dVgs_eff_dVg;
                    this._cgdb = 0.0;
                    this._cgsb = CoxWL * (dVbseff_dVb - dVgs_eff_dVg);

                    this._cdgb = 0.0;
                    this._cddb = 0.0;
                    this._cdsb = 0.0;

                    this._cbgb = -CoxWL * dVgs_eff_dVg;
                    this._cbdb = 0.0;
                    this._cbsb = -this._cgsb;
                    this._qinv = 0.0;
                }
                else if (Vgst <= 0.0)
                {
                    T1 = 0.5 * Param.BSIM3k1ox;
                    T2 = Math.Sqrt(T1 * T1 + Arg1);
                    qgate = CoxWL * Param.BSIM3k1ox * (T2 - T1);
                    qbulk = -qgate;
                    qdrn = 0.0;

                    T0 = CoxWL * T1 / T2;
                    this._cggb = T0 * dVgs_eff_dVg;
                    this._cgdb = 0.0;
                    this._cgsb = T0 * (dVbseff_dVb - dVgs_eff_dVg);

                    this._cdgb = 0.0;
                    this._cddb = 0.0;
                    this._cdsb = 0.0;

                    this._cbgb = -this._cggb;
                    this._cbdb = 0.0;
                    this._cbsb = -this._cgsb;
                    this._qinv = 0.0;
                }
                else
                {
                    One_Third_CoxWL = CoxWL / 3.0;
                    Two_Third_CoxWL = 2.0 * One_Third_CoxWL;

                    AbulkCV = Abulk0 * Param.BSIM3abulkCVfactor;
                    dAbulkCV_dVb = Param.BSIM3abulkCVfactor * dAbulk0_dVb;
                    Vdsat = Vgst / AbulkCV;
                    dVdsat_dVg = dVgs_eff_dVg / AbulkCV;
                    dVdsat_dVb = -(Vdsat * dAbulkCV_dVb + dVth_dVb) / AbulkCV;

                    if (ModelParameters.Xpart.Value > 0.5)
                    {   /* 0/100 Charge partition model */
                        if (Vdsat <= Vds)
                        {   /* saturation region */
                            T1 = Vdsat / 3.0;
                            qgate = CoxWL * (Vgs_eff - Vfb
                                  - Param.BSIM3phi - T1);
                            T2 = -Two_Third_CoxWL * Vgst;
                            qbulk = -(qgate + T2);
                            qdrn = 0.0;

                            this._cggb = One_Third_CoxWL * (3.0
                                            - dVdsat_dVg) * dVgs_eff_dVg;
                            T2 = -One_Third_CoxWL * dVdsat_dVb;
                            this._cgsb = -(this._cggb + T2);
                            this._cgdb = 0.0;

                            this._cdgb = 0.0;
                            this._cddb = 0.0;
                            this._cdsb = 0.0;

                            this._cbgb = -(this._cggb
                                            - Two_Third_CoxWL * dVgs_eff_dVg);
                            T3 = -(T2 + Two_Third_CoxWL * dVth_dVb);
                            this._cbsb = -(this._cbgb + T3);
                            this._cbdb = 0.0;
                            this._qinv = -(qgate + qbulk);
                        }
                        else
                        {   /* linear region */
                            Alphaz = Vgst / Vdsat;
                            T1 = 2.0 * Vdsat - Vds;
                            T2 = Vds / (3.0 * T1);
                            T3 = T2 * Vds;
                            T9 = 0.25 * CoxWL;
                            T4 = T9 * Alphaz;
                            T7 = 2.0 * Vds - T1 - 3.0 * T3;
                            T8 = T3 - T1 - 2.0 * Vds;
                            qgate = CoxWL * (Vgs_eff - Vfb
                                  - Param.BSIM3phi - 0.5 * (Vds - T3));
                            T10 = T4 * T8;
                            qdrn = T4 * T7;
                            qbulk = -(qgate + qdrn + T10);

                            T5 = T3 / T1;
                            this._cggb = CoxWL * (1.0 - T5 * dVdsat_dVg)
                                            * dVgs_eff_dVg;
                            T11 = -CoxWL * T5 * dVdsat_dVb;
                            this._cgdb = CoxWL * (T2 - 0.5 + 0.5 * T5);
                            this._cgsb = -(this._cggb + T11
                                            + this._cgdb);
                            T6 = 1.0 / Vdsat;
                            dAlphaz_dVg = T6 * (1.0 - Alphaz * dVdsat_dVg);
                            dAlphaz_dVb = -T6 * (dVth_dVb + Alphaz * dVdsat_dVb);
                            T7 = T9 * T7;
                            T8 = T9 * T8;
                            T9 = 2.0 * T4 * (1.0 - 3.0 * T5);
                            this._cdgb = (T7 * dAlphaz_dVg - T9
                                            * dVdsat_dVg) * dVgs_eff_dVg;
                            T12 = T7 * dAlphaz_dVb - T9 * dVdsat_dVb;
                            this._cddb = T4 * (3.0 - 6.0 * T2 - 3.0 * T5);
                            this._cdsb = -(this._cdgb + T12
                                            + this._cddb);

                            T9 = 2.0 * T4 * (1.0 + T5);
                            T10 = (T8 * dAlphaz_dVg - T9 * dVdsat_dVg)
                                * dVgs_eff_dVg;
                            T11 = T8 * dAlphaz_dVb - T9 * dVdsat_dVb;
                            T12 = T4 * (2.0 * T2 + T5 - 1.0);
                            T0 = -(T10 + T11 + T12);

                            this._cbgb = -(this._cggb
                                            + this._cdgb + T10);
                            this._cbdb = -(this._cgdb
                                            + this._cddb + T12);
                            this._cbsb = -(this._cgsb
                                            + this._cdsb + T0);
                            this._qinv = -(qgate + qbulk);
                        }
                    }
                    else if (ModelParameters.Xpart.Value < 0.5)
                    {   /* 40/60 Charge partition model */
                        if (Vds >= Vdsat)
                        {   /* saturation region */
                            T1 = Vdsat / 3.0;
                            qgate = CoxWL * (Vgs_eff - Vfb
                                  - Param.BSIM3phi - T1);
                            T2 = -Two_Third_CoxWL * Vgst;
                            qbulk = -(qgate + T2);
                            qdrn = 0.4 * T2;

                            this._cggb = One_Third_CoxWL * (3.0
                                            - dVdsat_dVg) * dVgs_eff_dVg;
                            T2 = -One_Third_CoxWL * dVdsat_dVb;
                            this._cgsb = -(this._cggb + T2);
                            this._cgdb = 0.0;

                            T3 = 0.4 * Two_Third_CoxWL;
                            this._cdgb = -T3 * dVgs_eff_dVg;
                            this._cddb = 0.0;
                            T4 = T3 * dVth_dVb;
                            this._cdsb = -(T4 + this._cdgb);

                            this._cbgb = -(this._cggb
                                            - Two_Third_CoxWL * dVgs_eff_dVg);
                            T3 = -(T2 + Two_Third_CoxWL * dVth_dVb);
                            this._cbsb = -(this._cbgb + T3);
                            this._cbdb = 0.0;
                            this._qinv = -(qgate + qbulk);
                        }
                        else
                        {   /* linear region  */
                            Alphaz = Vgst / Vdsat;
                            T1 = 2.0 * Vdsat - Vds;
                            T2 = Vds / (3.0 * T1);
                            T3 = T2 * Vds;
                            T9 = 0.25 * CoxWL;
                            T4 = T9 * Alphaz;
                            qgate = CoxWL * (Vgs_eff - Vfb - Param.BSIM3phi
                                  - 0.5 * (Vds - T3));

                            T5 = T3 / T1;
                            this._cggb = CoxWL * (1.0 - T5 * dVdsat_dVg)
                                            * dVgs_eff_dVg;
                            tmp = -CoxWL * T5 * dVdsat_dVb;
                            this._cgdb = CoxWL * (T2 - 0.5 + 0.5 * T5);
                            this._cgsb = -(this._cggb
                                            + this._cgdb + tmp);

                            T6 = 1.0 / Vdsat;
                            dAlphaz_dVg = T6 * (1.0 - Alphaz * dVdsat_dVg);
                            dAlphaz_dVb = -T6 * (dVth_dVb + Alphaz * dVdsat_dVb);

                            T6 = 8.0 * Vdsat * Vdsat - 6.0 * Vdsat * Vds
                               + 1.2 * Vds * Vds;
                            T8 = T2 / T1;
                            T7 = Vds - T1 - T8 * T6;
                            qdrn = T4 * T7;
                            T7 *= T9;
                            tmp = T8 / T1;
                            tmp1 = T4 * (2.0 - 4.0 * tmp * T6
                                 + T8 * (16.0 * Vdsat - 6.0 * Vds));

                            this._cdgb = (T7 * dAlphaz_dVg - tmp1
                                            * dVdsat_dVg) * dVgs_eff_dVg;
                            T10 = T7 * dAlphaz_dVb - tmp1 * dVdsat_dVb;
                            this._cddb = T4 * (2.0 - (1.0 / (3.0 * T1
                                            * T1) + 2.0 * tmp) * T6 + T8
                                            * (6.0 * Vdsat - 2.4 * Vds));
                            this._cdsb = -(this._cdgb
                                            + T10 + this._cddb);

                            T7 = 2.0 * (T1 + T3);
                            qbulk = -(qgate - T4 * T7);
                            T7 *= T9;
                            T0 = 4.0 * T4 * (1.0 - T5);
                            T12 = (-T7 * dAlphaz_dVg - this._cdgb
                                - T0 * dVdsat_dVg) * dVgs_eff_dVg;
                            T11 = -T7 * dAlphaz_dVb - T10 - T0 * dVdsat_dVb;
                            T10 = -4.0 * T4 * (T2 - 0.5 + 0.5 * T5)
                                - this._cddb;
                            tmp = -(T10 + T11 + T12);

                            this._cbgb = -(this._cggb
                                            + this._cdgb + T12);
                            this._cbdb = -(this._cgdb
                                            + this._cddb + T10); /* bug fix */
                            this._cbsb = -(this._cgsb
                                            + this._cdsb + tmp);
                            this._qinv = -(qgate + qbulk);
                        }
                    }
                    else
                    {   /* 50/50 partitioning */
                        if (Vds >= Vdsat)
                        {   /* saturation region */
                            T1 = Vdsat / 3.0;
                            qgate = CoxWL * (Vgs_eff - Vfb
                                  - Param.BSIM3phi - T1);
                            T2 = -Two_Third_CoxWL * Vgst;
                            qbulk = -(qgate + T2);
                            qdrn = 0.5 * T2;

                            this._cggb = One_Third_CoxWL * (3.0
                                            - dVdsat_dVg) * dVgs_eff_dVg;
                            T2 = -One_Third_CoxWL * dVdsat_dVb;
                            this._cgsb = -(this._cggb + T2);
                            this._cgdb = 0.0;

                            this._cdgb = -One_Third_CoxWL * dVgs_eff_dVg;
                            this._cddb = 0.0;
                            T4 = One_Third_CoxWL * dVth_dVb;
                            this._cdsb = -(T4 + this._cdgb);

                            this._cbgb = -(this._cggb
                                            - Two_Third_CoxWL * dVgs_eff_dVg);
                            T3 = -(T2 + Two_Third_CoxWL * dVth_dVb);
                            this._cbsb = -(this._cbgb + T3);
                            this._cbdb = 0.0;
                            this._qinv = -(qgate + qbulk);
                        }
                        else
                        {   /* linear region */
                            Alphaz = Vgst / Vdsat;
                            T1 = 2.0 * Vdsat - Vds;
                            T2 = Vds / (3.0 * T1);
                            T3 = T2 * Vds;
                            T9 = 0.25 * CoxWL;
                            T4 = T9 * Alphaz;
                            qgate = CoxWL * (Vgs_eff - Vfb - Param.BSIM3phi
                                  - 0.5 * (Vds - T3));

                            T5 = T3 / T1;
                            this._cggb = CoxWL * (1.0 - T5 * dVdsat_dVg)
                                            * dVgs_eff_dVg;
                            tmp = -CoxWL * T5 * dVdsat_dVb;
                            this._cgdb = CoxWL * (T2 - 0.5 + 0.5 * T5);
                            this._cgsb = -(this._cggb
                                            + this._cgdb + tmp);

                            T6 = 1.0 / Vdsat;
                            dAlphaz_dVg = T6 * (1.0 - Alphaz * dVdsat_dVg);
                            dAlphaz_dVb = -T6 * (dVth_dVb + Alphaz * dVdsat_dVb);

                            T7 = T1 + T3;
                            qdrn = -T4 * T7;
                            qbulk = -(qgate + qdrn + qdrn);
                            T7 *= T9;
                            T0 = T4 * (2.0 * T5 - 2.0);

                            this._cdgb = (T0 * dVdsat_dVg - T7
                                            * dAlphaz_dVg) * dVgs_eff_dVg;
                            T12 = T0 * dVdsat_dVb - T7 * dAlphaz_dVb;
                            this._cddb = T4 * (1.0 - 2.0 * T2 - T5);
                            this._cdsb = -(this._cdgb + T12
                                            + this._cddb);

                            this._cbgb = -(this._cggb
                                            + 2.0 * this._cdgb);
                            this._cbdb = -(this._cgdb
                                            + 2.0 * this._cddb);
                            this._cbsb = -(this._cgsb
                                            + 2.0 * this._cdsb);
                            this._qinv = -(qgate + qbulk);
                        }
                    }
                }
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
                    VbseffCV = Param.BSIM3phi - Phis;
                    dVbseffCV_dVb = -dPhis_dVb;
                }

                CoxWL = ModelTemperature.Cox * Param.BSIM3weffCV
                      * Param.BSIM3leffCV;

                /* Seperate VgsteffCV with noff and voffcv */
                noff = n * Param.BSIM3noff;
                dnoff_dVd = Param.BSIM3noff * dn_dVd;
                dnoff_dVb = Param.BSIM3noff * dn_dVb;
                T0 = Vtm * noff;
                voffcv = Param.BSIM3voffcv;
                VgstNVt = (Vgst - voffcv) / T0;

                if (VgstNVt > EXP_THRESHOLD)
                {
                    Vgsteff = Vgst - voffcv;
                    dVgsteff_dVg = dVgs_eff_dVg;
                    dVgsteff_dVd = -dVth_dVd;
                    dVgsteff_dVb = -dVth_dVb;
                }
                else if (VgstNVt < -EXP_THRESHOLD)
                {
                    Vgsteff = T0 * Math.Log(1.0 + MIN_EXP);
                    dVgsteff_dVg = 0.0;
                    dVgsteff_dVd = Vgsteff / noff;
                    dVgsteff_dVb = dVgsteff_dVd * dnoff_dVb;
                    dVgsteff_dVd *= dnoff_dVd;
                }
                else
                {
                    ExpVgst = Math.Exp(VgstNVt);
                    Vgsteff = T0 * Math.Log(1.0 + ExpVgst);
                    dVgsteff_dVg = ExpVgst / (1.0 + ExpVgst);
                    dVgsteff_dVd = -dVgsteff_dVg * (dVth_dVd + (Vgst - voffcv)
                                 / noff * dnoff_dVd) + Vgsteff / noff * dnoff_dVd;
                    dVgsteff_dVb = -dVgsteff_dVg * (dVth_dVb + (Vgst - voffcv)
                                 / noff * dnoff_dVb) + Vgsteff / noff * dnoff_dVb;
                    dVgsteff_dVg *= dVgs_eff_dVg;
                } /* End of VgsteffCV */

                if (ModelParameters.CapMod.Value == 1)
                {
                    Vfb = this.Vfbzb;
                    Arg1 = Vgs_eff - VbseffCV - Vfb - Vgsteff;

                    if (Arg1 <= 0.0)
                    {
                        qgate = CoxWL * Arg1;
                        Cgg = CoxWL * (dVgs_eff_dVg - dVgsteff_dVg);
                        Cgd = -CoxWL * dVgsteff_dVd;
                        Cgb = -CoxWL * (dVbseffCV_dVb + dVgsteff_dVb);
                    }
                    else
                    {
                        T0 = 0.5 * Param.BSIM3k1ox;
                        T1 = Math.Sqrt(T0 * T0 + Arg1);
                        T2 = CoxWL * T0 / T1;

                        qgate = CoxWL * Param.BSIM3k1ox * (T1 - T0);

                        Cgg = T2 * (dVgs_eff_dVg - dVgsteff_dVg);
                        Cgd = -T2 * dVgsteff_dVd;
                        Cgb = -T2 * (dVbseffCV_dVb + dVgsteff_dVb);
                    }
                    qbulk = -qgate;
                    Cbg = -Cgg;
                    Cbd = -Cgd;
                    Cbb = -Cgb;

                    One_Third_CoxWL = CoxWL / 3.0;
                    Two_Third_CoxWL = 2.0 * One_Third_CoxWL;
                    AbulkCV = Abulk0 * Param.BSIM3abulkCVfactor;
                    dAbulkCV_dVb = Param.BSIM3abulkCVfactor * dAbulk0_dVb;
                    VdsatCV = Vgsteff / AbulkCV;
                    if (VdsatCV < Vds)
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

                        if (ModelParameters.Xpart.Value > 0.5)
                            T0 = -Two_Third_CoxWL;
                        else if (ModelParameters.Xpart.Value < 0.5)
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
                        T1 = 12.0 * (Vgsteff - 0.5 * T0 + 1.0e-20);
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

                        if (ModelParameters.Xpart.Value > 0.5)
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
                        else if (ModelParameters.Xpart.Value < 0.5)
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
                    this._qinv = -(qgate + qbulk);
                }

                else if (ModelParameters.CapMod.Value == 2)
                {
                    Vfb = this.Vfbzb;
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
                    dVfbeff_dVg = T1 * dVgs_eff_dVg;
                    dVfbeff_dVb = -T1 * dVbseffCV_dVb;
                    Qac0 = CoxWL * (Vfbeff - Vfb);
                    dQac0_dVg = CoxWL * dVfbeff_dVg;
                    dQac0_dVb = CoxWL * dVfbeff_dVb;

                    T0 = 0.5 * Param.BSIM3k1ox;
                    T3 = Vgs_eff - Vfbeff - VbseffCV - Vgsteff;
                    if (Param.BSIM3k1ox == 0.0)
                    {
                        T1 = 0.0;
                        T2 = 0.0;
                    }
                    else if (T3 < 0.0)
                    {
                        T1 = T0 + T3 / Param.BSIM3k1ox;
                        T2 = CoxWL;
                    }
                    else
                    {
                        T1 = Math.Sqrt(T0 * T0 + T3);
                        T2 = CoxWL * T0 / T1;
                    }

                    Qsub0 = CoxWL * Param.BSIM3k1ox * (T1 - T0);

                    dQsub0_dVg = T2 * (dVgs_eff_dVg - dVfbeff_dVg - dVgsteff_dVg);
                    dQsub0_dVd = -T2 * dVgsteff_dVd;
                    dQsub0_dVb = -T2 * (dVfbeff_dVb + dVbseffCV_dVb
                               + dVgsteff_dVb);

                    AbulkCV = Abulk0 * Param.BSIM3abulkCVfactor;
                    dAbulkCV_dVb = Param.BSIM3abulkCVfactor * dAbulk0_dVb;
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
                    /* Added to eliminate non-zero VdseffCV at Vds=0.0 */
                    if (Vds == 0.0)
                    {
                        VdseffCV = 0.0;
                        dVdseffCV_dVg = 0.0;
                        dVdseffCV_dVb = 0.0;
                    }

                    T0 = AbulkCV * VdseffCV;
                    T1 = 12.0 * (Vgsteff - 0.5 * T0 + 1e-20);
                    T2 = VdseffCV / T1;
                    T3 = T0 * T2;

                    T4 = (1.0 - 12.0 * T2 * T2 * AbulkCV);
                    T5 = (6.0 * T0 * (4.0 * Vgsteff - T0) / (T1 * T1) - 0.5);
                    T6 = 12.0 * T2 * T2 * Vgsteff;

                    qinoi = -CoxWL * (Vgsteff - 0.5 * T0 + AbulkCV * T3);
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

                    if (ModelParameters.Xpart.Value > 0.5)
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
                    else if (ModelParameters.Xpart.Value < 0.5)
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
                    Cgd = dQsub0_dVd + Cgd1;
                    Cgb = dQac0_dVb + dQsub0_dVb + Cgb1;

                    Cbg = Cbg1 - dQac0_dVg - dQsub0_dVg;
                    Cbd = Cbd1 - dQsub0_dVd;
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
                    this._qinv = qinoi;
                }

                /* New Charge-Thickness capMod (CTM) begins */
                else if (ModelParameters.CapMod.Value == 3)
                {
                    V3 = this.Vfbzb - Vgs_eff + VbseffCV - DELTA_3;
                    if (this.Vfbzb <= 0.0)
                    {
                        T0 = Math.Sqrt(V3 * V3 - 4.0 * DELTA_3 * this.Vfbzb);
                        T2 = -DELTA_3 / T0;
                    }
                    else
                    {
                        T0 = Math.Sqrt(V3 * V3 + 4.0 * DELTA_3 * this.Vfbzb);
                        T2 = DELTA_3 / T0;
                    }

                    T1 = 0.5 * (1.0 + V3 / T0);
                    Vfbeff = this.Vfbzb - 0.5 * (V3 + T0);
                    dVfbeff_dVg = T1 * dVgs_eff_dVg;
                    dVfbeff_dVb = -T1 * dVbseffCV_dVb;

                    Cox = ModelTemperature.Cox;
                    Tox = 1.0e8 * ModelParameters.Tox;
                    T0 = (Vgs_eff - VbseffCV - this.Vfbzb) / Tox;
                    dT0_dVg = dVgs_eff_dVg / Tox;
                    dT0_dVb = -dVbseffCV_dVb / Tox;

                    tmp = T0 * Param.BSIM3acde;
                    if ((-EXP_THRESHOLD < tmp) && (tmp < EXP_THRESHOLD))
                    {
                        Tcen = Param.BSIM3ldeb * Math.Exp(tmp);
                        dTcen_dVg = Param.BSIM3acde * Tcen;
                        dTcen_dVb = dTcen_dVg * dT0_dVb;
                        dTcen_dVg *= dT0_dVg;
                    }
                    else if (tmp <= -EXP_THRESHOLD)
                    {
                        Tcen = Param.BSIM3ldeb * MIN_EXP;
                        dTcen_dVg = dTcen_dVb = 0.0;
                    }
                    else
                    {
                        Tcen = Param.BSIM3ldeb * MAX_EXP;
                        dTcen_dVg = dTcen_dVb = 0.0;
                    }

                    LINK = 1.0e-3 * ModelParameters.Tox;
                    V3 = Param.BSIM3ldeb - Tcen - LINK;
                    V4 = Math.Sqrt(V3 * V3 + 4.0 * LINK * Param.BSIM3ldeb);
                    Tcen = Param.BSIM3ldeb - 0.5 * (V3 + V4);
                    T1 = 0.5 * (1.0 + V3 / V4);
                    dTcen_dVg *= T1;
                    dTcen_dVb *= T1;

                    Ccen = EPSSI / Tcen;
                    T2 = Cox / (Cox + Ccen);
                    Coxeff = T2 * Ccen;
                    T3 = -Ccen / Tcen;
                    dCoxeff_dVg = T2 * T2 * T3;
                    dCoxeff_dVb = dCoxeff_dVg * dTcen_dVb;
                    dCoxeff_dVg *= dTcen_dVg;
                    CoxWLcen = CoxWL * Coxeff / Cox;

                    Qac0 = CoxWLcen * (Vfbeff - this.Vfbzb);
                    QovCox = Qac0 / Coxeff;
                    dQac0_dVg = CoxWLcen * dVfbeff_dVg
                              + QovCox * dCoxeff_dVg;
                    dQac0_dVb = CoxWLcen * dVfbeff_dVb
                              + QovCox * dCoxeff_dVb;

                    T0 = 0.5 * Param.BSIM3k1ox;
                    T3 = Vgs_eff - Vfbeff - VbseffCV - Vgsteff;
                    if (Param.BSIM3k1ox == 0.0)
                    {
                        T1 = 0.0;
                        T2 = 0.0;
                    }
                    else if (T3 < 0.0)
                    {
                        T1 = T0 + T3 / Param.BSIM3k1ox;
                        T2 = CoxWLcen;
                    }
                    else
                    {
                        T1 = Math.Sqrt(T0 * T0 + T3);
                        T2 = CoxWLcen * T0 / T1;
                    }

                    Qsub0 = CoxWLcen * Param.BSIM3k1ox * (T1 - T0);
                    QovCox = Qsub0 / Coxeff;
                    dQsub0_dVg = T2 * (dVgs_eff_dVg - dVfbeff_dVg - dVgsteff_dVg)
                               + QovCox * dCoxeff_dVg;
                    dQsub0_dVd = -T2 * dVgsteff_dVd;
                    dQsub0_dVb = -T2 * (dVfbeff_dVb + dVbseffCV_dVb + dVgsteff_dVb)
                               + QovCox * dCoxeff_dVb;

                    /* Gate-bias dependent delta Phis begins */
                    if (Param.BSIM3k1ox <= 0.0)
                    {
                        Denomi = 0.25 * Param.BSIM3moin * Vtm;
                        T0 = 0.5 * Param.BSIM3sqrtPhi;
                    }
                    else
                    {
                        Denomi = Param.BSIM3moin * Vtm
                               * Param.BSIM3k1ox * Param.BSIM3k1ox;
                        T0 = Param.BSIM3k1ox * Param.BSIM3sqrtPhi;
                    }
                    T1 = 2.0 * T0 + Vgsteff;

                    DeltaPhi = Vtm * Math.Log(1.0 + T1 * Vgsteff / Denomi);
                    dDeltaPhi_dVg = 2.0 * Vtm * (T1 - T0) / (Denomi + T1 * Vgsteff);
                    /* End of delta Phis */

                    /* VgDP = Vgsteff - DeltaPhi */
                    T0 = Vgsteff - DeltaPhi - 0.001;
                    dT0_dVg = 1.0 - dDeltaPhi_dVg;
                    T1 = Math.Sqrt(T0 * T0 + Vgsteff * 0.004);
                    VgDP = 0.5 * (T0 + T1);
                    dVgDP_dVg = 0.5 * (dT0_dVg + (T0 * dT0_dVg + 0.002) / T1);

                    T3 = 4.0 * (Vth - this.Vfbzb - Param.BSIM3phi);
                    Tox += Tox;
                    if (T3 >= 0.0)
                    {
                        T0 = (Vgsteff + T3) / Tox;
                        dT0_dVd = (dVgsteff_dVd + 4.0 * dVth_dVd) / Tox;
                        dT0_dVb = (dVgsteff_dVb + 4.0 * dVth_dVb) / Tox;
                    }
                    else
                    {
                        T0 = (Vgsteff + 1.0e-20) / Tox;
                        dT0_dVd = dVgsteff_dVd / Tox;
                        dT0_dVb = dVgsteff_dVb / Tox;
                    }
                    tmp = Math.Exp(0.7 * Math.Log(T0));
                    T1 = 1.0 + tmp;
                    T2 = 0.7 * tmp / (T0 * Tox);
                    Tcen = 1.9e-9 / T1;
                    dTcen_dVg = -1.9e-9 * T2 / T1 / T1;
                    dTcen_dVd = Tox * dTcen_dVg;
                    dTcen_dVb = dTcen_dVd * dT0_dVb;
                    dTcen_dVd *= dT0_dVd;
                    dTcen_dVg *= dVgsteff_dVg;

                    Ccen = EPSSI / Tcen;
                    T0 = Cox / (Cox + Ccen);
                    Coxeff = T0 * Ccen;
                    T1 = -Ccen / Tcen;
                    dCoxeff_dVg = T0 * T0 * T1;
                    dCoxeff_dVd = dCoxeff_dVg * dTcen_dVd;
                    dCoxeff_dVb = dCoxeff_dVg * dTcen_dVb;
                    dCoxeff_dVg *= dTcen_dVg;
                    CoxWLcen = CoxWL * Coxeff / Cox;

                    AbulkCV = Abulk0 * Param.BSIM3abulkCVfactor;
                    dAbulkCV_dVb = Param.BSIM3abulkCVfactor * dAbulk0_dVb;
                    VdsatCV = VgDP / AbulkCV;
                    T0 = VdsatCV - Vds - DELTA_4;
                    dT0_dVg = dVgDP_dVg / AbulkCV;
                    dT0_dVb = -VdsatCV * dAbulkCV_dVb / AbulkCV;
                    T1 = Math.Sqrt(T0 * T0 + 4.0 * DELTA_4 * VdsatCV);
                    dT1_dVg = (T0 + DELTA_4 + DELTA_4) / T1;
                    dT1_dVd = -T0 / T1;
                    dT1_dVb = dT1_dVg * dT0_dVb;
                    dT1_dVg *= dT0_dVg;
                    if (T0 >= 0.0)
                    {
                        VdseffCV = VdsatCV - 0.5 * (T0 + T1);
                        dVdseffCV_dVg = 0.5 * (dT0_dVg - dT1_dVg);
                        dVdseffCV_dVd = 0.5 * (1.0 - dT1_dVd);
                        dVdseffCV_dVb = 0.5 * (dT0_dVb - dT1_dVb);
                    }
                    else
                    {
                        T3 = (DELTA_4 + DELTA_4) / (T1 - T0);
                        T4 = 1.0 - T3;
                        T5 = VdsatCV * T3 / (T1 - T0);
                        VdseffCV = VdsatCV * T4;
                        dVdseffCV_dVg = dT0_dVg * T4 + T5 * (dT1_dVg - dT0_dVg);
                        dVdseffCV_dVd = T5 * (dT1_dVd + 1.0);
                        dVdseffCV_dVb = dT0_dVb * (1.0 - T5) + T5 * dT1_dVb;
                    }

                    /* Added to eliminate non-zero VdseffCV at Vds=0.0 */
                    if (Vds == 0.0)
                    {
                        VdseffCV = 0.0;
                        dVdseffCV_dVg = 0.0;
                        dVdseffCV_dVb = 0.0;
                    }

                    T0 = AbulkCV * VdseffCV;
                    T1 = VgDP;
                    T2 = 12.0 * (T1 - 0.5 * T0 + 1.0e-20);
                    T3 = T0 / T2;
                    T4 = 1.0 - 12.0 * T3 * T3;
                    T5 = AbulkCV * (6.0 * T0 * (4.0 * T1 - T0) / (T2 * T2) - 0.5);
                    T6 = T5 * VdseffCV / AbulkCV;

                    qgate = qinoi = CoxWLcen * (T1 - T0 * (0.5 - T3));
                    QovCox = qgate / Coxeff;
                    Cgg1 = CoxWLcen * (T4 * dVgDP_dVg
                         + T5 * dVdseffCV_dVg);
                    Cgd1 = CoxWLcen * T5 * dVdseffCV_dVd + Cgg1
                         * dVgsteff_dVd + QovCox * dCoxeff_dVd;
                    Cgb1 = CoxWLcen * (T5 * dVdseffCV_dVb + T6 * dAbulkCV_dVb)
                         + Cgg1 * dVgsteff_dVb + QovCox * dCoxeff_dVb;
                    Cgg1 = Cgg1 * dVgsteff_dVg + QovCox * dCoxeff_dVg;


                    T7 = 1.0 - AbulkCV;
                    T8 = T2 * T2;
                    T9 = 12.0 * T7 * T0 * T0 / (T8 * AbulkCV);
                    T10 = T9 * dVgDP_dVg;
                    T11 = -T7 * T5 / AbulkCV;
                    T12 = -(T9 * T1 / AbulkCV + VdseffCV * (0.5 - T0 / T2));

                    qbulk = CoxWLcen * T7 * (0.5 * VdseffCV - T0 * VdseffCV / T2);
                    QovCox = qbulk / Coxeff;
                    Cbg1 = CoxWLcen * (T10 + T11 * dVdseffCV_dVg);
                    Cbd1 = CoxWLcen * T11 * dVdseffCV_dVd + Cbg1
                         * dVgsteff_dVd + QovCox * dCoxeff_dVd;
                    Cbb1 = CoxWLcen * (T11 * dVdseffCV_dVb + T12 * dAbulkCV_dVb)
                         + Cbg1 * dVgsteff_dVb + QovCox * dCoxeff_dVb;
                    Cbg1 = Cbg1 * dVgsteff_dVg + QovCox * dCoxeff_dVg;

                    if (ModelParameters.Xpart.Value > 0.5)
                    {   /* 0/100 partition */
                        qsrc = -CoxWLcen * (T1 / 2.0 + T0 / 4.0
                             - 0.5 * T0 * T0 / T2);
                        QovCox = qsrc / Coxeff;
                        T2 += T2;
                        T3 = T2 * T2;
                        T7 = -(0.25 - 12.0 * T0 * (4.0 * T1 - T0) / T3);
                        T4 = -(0.5 + 24.0 * T0 * T0 / T3) * dVgDP_dVg;
                        T5 = T7 * AbulkCV;
                        T6 = T7 * VdseffCV;

                        Csg = CoxWLcen * (T4 + T5 * dVdseffCV_dVg);
                        Csd = CoxWLcen * T5 * dVdseffCV_dVd + Csg * dVgsteff_dVd
                            + QovCox * dCoxeff_dVd;
                        Csb = CoxWLcen * (T5 * dVdseffCV_dVb + T6 * dAbulkCV_dVb)
                            + Csg * dVgsteff_dVb + QovCox * dCoxeff_dVb;
                        Csg = Csg * dVgsteff_dVg + QovCox * dCoxeff_dVg;
                    }
                    else if (ModelParameters.Xpart.Value < 0.5)
                    {   /* 40/60 partition */
                        T2 = T2 / 12.0;
                        T3 = 0.5 * CoxWLcen / (T2 * T2);
                        T4 = T1 * (2.0 * T0 * T0 / 3.0 + T1 * (T1 - 4.0
                           * T0 / 3.0)) - 2.0 * T0 * T0 * T0 / 15.0;
                        qsrc = -T3 * T4;
                        QovCox = qsrc / Coxeff;
                        T8 = 4.0 / 3.0 * T1 * (T1 - T0) + 0.4 * T0 * T0;
                        T5 = -2.0 * qsrc / T2 - T3 * (T1 * (3.0 * T1 - 8.0
                           * T0 / 3.0) + 2.0 * T0 * T0 / 3.0);
                        T6 = AbulkCV * (qsrc / T2 + T3 * T8);
                        T7 = T6 * VdseffCV / AbulkCV;

                        Csg = T5 * dVgDP_dVg + T6 * dVdseffCV_dVg;
                        Csd = Csg * dVgsteff_dVd + T6 * dVdseffCV_dVd
                            + QovCox * dCoxeff_dVd;
                        Csb = Csg * dVgsteff_dVb + T6 * dVdseffCV_dVb
                            + T7 * dAbulkCV_dVb + QovCox * dCoxeff_dVb;
                        Csg = Csg * dVgsteff_dVg + QovCox * dCoxeff_dVg;
                    }
                    else
                    {   /* 50/50 partition */
                        qsrc = -0.5 * qgate;
                        Csg = -0.5 * Cgg1;
                        Csd = -0.5 * Cgd1;
                        Csb = -0.5 * Cgb1;
                    }

                    qgate += Qac0 + Qsub0 - qbulk;
                    qbulk -= (Qac0 + Qsub0);
                    qdrn = -(qgate + qbulk + qsrc);

                    Cbg = Cbg1 - dQac0_dVg - dQsub0_dVg;
                    Cbd = Cbd1 - dQsub0_dVd;
                    Cbb = Cbb1 - dQac0_dVb - dQsub0_dVb;

                    Cgg = Cgg1 - Cbg;
                    Cgd = Cgd1 - Cbd;
                    Cgb = Cgb1 - Cbb;

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
                    this._qinv = -qinoi;
                }  /* End of CTM */
            }

        finished:
            /* Returning Values to Calling Routine */
            /*
             *  COMPUTE EQUIVALENT DRAIN CURRENT SOURCE
             */

            this._qgate = qgate;
            this._qbulk = qbulk;
            this._qdrn = qdrn;
            this._cd = cdrain;

            if (ChargeComputationNeeded)
            {   /*  charge storage elements
               *  bulk-drain and bulk-source depletion capacitances
               *  czbd : zero bias drain junction capacitance
               *  czbs : zero bias source junction capacitance
               *  czbdsw: zero bias drain junction sidewall capacitance
                          along field oxide
               *  czbssw: zero bias source junction sidewall capacitance
                          along field oxide
               *  czbdswg: zero bias drain junction sidewall capacitance
                           along gate side
               *  czbsswg: zero bias source junction sidewall capacitance
                           along gate side
               */

                if (ModelParameters.AcmMod.Value == 0)
                {
                    czbd = ModelTemperature.UnitAreaTempJctCap * Parameters.DrainArea; /*bug fix */
                    czbs = ModelTemperature.UnitAreaTempJctCap * Parameters.SourceArea;
                    if (Parameters.DrainPerimeter < Param.BSIM3weff)
                    {
                        czbdswg = ModelTemperature.UnitLengthGateSidewallTempJctCap
                                * Parameters.DrainPerimeter;
                        czbdsw = 0.0;
                    }
                    else
                    {
                        czbdsw = ModelTemperature.UnitLengthSidewallTempJctCap
                               * (Parameters.DrainPerimeter - Param.BSIM3weff);
                        czbdswg = ModelTemperature.UnitLengthGateSidewallTempJctCap
                                * Param.BSIM3weff;
                    }
                    if (Parameters.SourcePerimeter < Param.BSIM3weff)
                    {
                        czbssw = 0.0;
                        czbsswg = ModelTemperature.UnitLengthGateSidewallTempJctCap
                                * Parameters.SourcePerimeter;
                    }
                    else
                    {
                        czbssw = ModelTemperature.UnitLengthSidewallTempJctCap
                               * (Parameters.SourcePerimeter - Param.BSIM3weff);
                        czbsswg = ModelTemperature.UnitLengthGateSidewallTempJctCap
                                * Param.BSIM3weff;
                    }
                }
                else
                {
                    ACM.JunctionCapacitances(
                        ModelParameters.AcmMod,
                        ModelParameters.Calcacm,
                        Parameters.Geo,
                        ModelParameters.Hdif,
                        ModelParameters.Wmlt,
                        Parameters.Width,
                        ModelParameters.Xw,
                        Parameters.DrainArea,
                        Parameters.DrainPerimeter,
                        Parameters.SourceArea,
                        Parameters.SourcePerimeter,
                        ModelTemperature.UnitAreaTempJctCap,
                        ModelTemperature.UnitLengthSidewallTempJctCap,
                        ModelTemperature.UnitLengthGateSidewallTempJctCap, // Changed to temperature-dependent version - Sven Boulanger 20220721
                        out czbd,
                        out czbdsw,
                        out czbdswg,
                        out czbs,
                        out czbssw,
                        out czbsswg
                    );
                }

                MJ = ModelParameters.BulkJctBotGradingCoeff;
                MJSW = ModelParameters.BulkJctSideGradingCoeff;
                MJSWG = ModelParameters.BulkJctGateSideGradingCoeff;

                /* Source Bulk Junction */
                if (vbs == 0.0)
                {
                    this._qbs = 0.0;
                    this._capbs = czbs + czbssw + czbsswg;
                }
                else if (vbs < 0.0)
                {
                    if (czbs > 0.0)
                    {
                        arg = 1.0 - vbs / ModelTemperature.PhiB;
                        if (MJ == 0.5)
                            sarg = 1.0 / Math.Sqrt(arg);
                        else
                            sarg = Math.Exp(-MJ * Math.Log(arg));
                        this._qbs = ModelTemperature.PhiB * czbs
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
                        arg = 1.0 - vbs / ModelTemperature.PhiBSW;
                        if (MJSW == 0.5)
                            sarg = 1.0 / Math.Sqrt(arg);
                        else
                            sarg = Math.Exp(-MJSW * Math.Log(arg));
                        this._qbs += ModelTemperature.PhiBSW * czbssw
                                         * (1.0 - arg * sarg) / (1.0 - MJSW);
                        this._capbs += czbssw * sarg;
                    }
                    if (czbsswg > 0.0)
                    {
                        arg = 1.0 - vbs / ModelTemperature.PhiBSWG;
                        if (MJSWG == 0.5)
                            sarg = 1.0 / Math.Sqrt(arg);
                        else
                            sarg = Math.Exp(-MJSWG * Math.Log(arg));
                        this._qbs += ModelTemperature.PhiBSWG * czbsswg
                                         * (1.0 - arg * sarg) / (1.0 - MJSWG);
                        this._capbs += czbsswg * sarg;
                    }

                }
                else
                {
                    T0 = czbs + czbssw + czbsswg;
                    T1 = vbs * (czbs * MJ / ModelTemperature.PhiB + czbssw * MJSW
                       / ModelTemperature.PhiBSW + czbsswg * MJSWG / ModelTemperature.PhiBSWG);
                    this._qbs = vbs * (T0 + 0.5 * T1);
                    this._capbs = T0 + T1;
                }

                /* Drain Bulk Junction */
                if (vbd == 0.0)
                {
                    this._qbd = 0.0;
                    this._capbd = czbd + czbdsw + czbdswg;
                }
                else if (vbd < 0.0)
                {
                    if (czbd > 0.0)
                    {
                        arg = 1.0 - vbd / ModelTemperature.PhiB;
                        if (MJ == 0.5)
                            sarg = 1.0 / Math.Sqrt(arg);
                        else
                            sarg = Math.Exp(-MJ * Math.Log(arg));
                        this._qbd = ModelTemperature.PhiB * czbd
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
                        arg = 1.0 - vbd / ModelTemperature.PhiBSW;
                        if (MJSW == 0.5)
                            sarg = 1.0 / Math.Sqrt(arg);
                        else
                            sarg = Math.Exp(-MJSW * Math.Log(arg));
                        this._qbd += ModelTemperature.PhiBSW * czbdsw
                                         * (1.0 - arg * sarg) / (1.0 - MJSW);
                        this._capbd += czbdsw * sarg;
                    }
                    if (czbdswg > 0.0)
                    {
                        arg = 1.0 - vbd / ModelTemperature.PhiBSWG;
                        if (MJSWG == 0.5)
                            sarg = 1.0 / Math.Sqrt(arg);
                        else
                            sarg = Math.Exp(-MJSWG * Math.Log(arg));
                        this._qbd += ModelTemperature.PhiBSWG * czbdswg
                                         * (1.0 - arg * sarg) / (1.0 - MJSWG);
                        this._capbd += czbdswg * sarg;
                    }
                }
                else
                {
                    T0 = czbd + czbdsw + czbdswg;
                    T1 = vbd * (czbd * MJ / ModelTemperature.PhiB + czbdsw * MJSW
                       / ModelTemperature.PhiBSW + czbdswg * MJSWG / ModelTemperature.PhiBSWG);
                    this._qbd = vbd * (T0 + 0.5 * T1);
                    this._capbd = T0 + T1;
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
            this._qdef = qdef;

            /* bulk and channel charge plus overlaps */

            if (!ChargeComputationNeeded)
                goto line850;

            /* NQS begins */
            if ((Parameters.NqsMod.Value != 0) || (Parameters.AcnqsMod.Value != 0))
            {
                qcheq = -(qbulk + qgate);

                this._cqgb = -(this._cggb + this._cbgb);
                this._cqdb = -(this._cgdb + this._cbdb);
                this._cqsb = -(this._cgsb + this._cbsb);
                this._cqbb = -(this._cqgb + this._cqdb
                                + this._cqsb);

                gtau_drift = Math.Abs(this.Tconst * qcheq) * ScalingFactor;
                T0 = Param.BSIM3leffCV * Param.BSIM3leffCV;
                gtau_diff = 16.0 * this.U0temp * ModelTemperature.Vtm / T0
                          * ScalingFactor;
                this._gtau = gtau_drift + gtau_diff;
                if (Parameters.AcnqsMod.Value != 0)
                    this._taunet = ScalingFactor / this._gtau;

            }

            if (ModelParameters.CapMod.Value == 0) /* code merge -JX */
            {
                cgdo = Param.BSIM3cgdo;
                qgdo = Param.BSIM3cgdo * vgd;
                cgso = Param.BSIM3cgso;
                qgso = Param.BSIM3cgso * vgs;
            }
            else if (ModelParameters.CapMod.Value == 1)
            {
                if (vgd < 0.0)
                {
                    T1 = Math.Sqrt(1.0 - 4.0 * vgd / Param.BSIM3ckappa);
                    cgdo = Param.BSIM3cgdo + Param.BSIM3weffCV
                         * Param.BSIM3cgdl / T1;
                    qgdo = Param.BSIM3cgdo * vgd - Param.BSIM3weffCV * 0.5
                         * Param.BSIM3cgdl * Param.BSIM3ckappa * (T1 - 1.0);
                }
                else
                {
                    cgdo = Param.BSIM3cgdo + Param.BSIM3weffCV
                         * Param.BSIM3cgdl;
                    qgdo = (Param.BSIM3weffCV * Param.BSIM3cgdl
                         + Param.BSIM3cgdo) * vgd;
                }

                if (vgs < 0.0)
                {
                    T1 = Math.Sqrt(1.0 - 4.0 * vgs / Param.BSIM3ckappa);
                    cgso = Param.BSIM3cgso + Param.BSIM3weffCV
                         * Param.BSIM3cgsl / T1;
                    qgso = Param.BSIM3cgso * vgs - Param.BSIM3weffCV * 0.5
                         * Param.BSIM3cgsl * Param.BSIM3ckappa * (T1 - 1.0);
                }
                else
                {
                    cgso = Param.BSIM3cgso + Param.BSIM3weffCV
                         * Param.BSIM3cgsl;
                    qgso = (Param.BSIM3weffCV * Param.BSIM3cgsl
                         + Param.BSIM3cgso) * vgs;
                }
            }
            else
            {
                T0 = vgd + DELTA_1;
                T1 = Math.Sqrt(T0 * T0 + 4.0 * DELTA_1);
                T2 = 0.5 * (T0 - T1);

                T3 = Param.BSIM3weffCV * Param.BSIM3cgdl;
                T4 = Math.Sqrt(1.0 - 4.0 * T2 / Param.BSIM3ckappa);
                cgdo = Param.BSIM3cgdo + T3 - T3 * (1.0 - 1.0 / T4)
                     * (0.5 - 0.5 * T0 / T1);
                qgdo = (Param.BSIM3cgdo + T3) * vgd - T3 * (T2
                     + 0.5 * Param.BSIM3ckappa * (T4 - 1.0));

                T0 = vgs + DELTA_1;
                T1 = Math.Sqrt(T0 * T0 + 4.0 * DELTA_1);
                T2 = 0.5 * (T0 - T1);
                T3 = Param.BSIM3weffCV * Param.BSIM3cgsl;
                T4 = Math.Sqrt(1.0 - 4.0 * T2 / Param.BSIM3ckappa);
                cgso = Param.BSIM3cgso + T3 - T3 * (1.0 - 1.0 / T4)
                     * (0.5 - 0.5 * T0 / T1);
                qgso = (Param.BSIM3cgso + T3) * vgs - T3 * (T2
                     + 0.5 * Param.BSIM3ckappa * (T4 - 1.0));
            }

            this.Cgdo = cgdo;
            this.Cgso = cgso;

            if (_method != null)
                ag0 = _method.Slope;
            else
                ag0 = 0;
            if (this._mode > 0)
            {
                if (Parameters.NqsMod.Value == 0)
                {
                    gcggb = (this._cggb + cgdo + cgso
                          + Param.BSIM3cgbo) * ag0;
                    gcgdb = (this._cgdb - cgdo) * ag0;
                    gcgsb = (this._cgsb - cgso) * ag0;

                    gcdgb = (this._cdgb - cgdo) * ag0;
                    gcddb = (this._cddb + this._capbd + cgdo) * ag0;
                    gcdsb = this._cdsb * ag0;

                    gcsgb = -(this._cggb + this._cbgb
                          + this._cdgb + cgso) * ag0;
                    gcsdb = -(this._cgdb + this._cbdb
                          + this._cddb) * ag0;
                    gcssb = (this._capbs + cgso - (this._cgsb
                          + this._cbsb + this._cdsb)) * ag0;

                    gcbgb = (this._cbgb - Param.BSIM3cgbo) * ag0;
                    gcbdb = (this._cbdb - this._capbd) * ag0;
                    gcbsb = (this._cbsb - this._capbs) * ag0;

                    qgd = qgdo;
                    qgs = qgso;
                    qgb = Param.BSIM3cgbo * vgb;
                    qgate += qgd + qgs + qgb;
                    qbulk -= qgb;
                    qdrn -= qgd;
                    qsrc = -(qgate + qbulk + qdrn);

                    ggtg = ggtd = ggtb = ggts = 0.0;
                    sxpart = 0.6;
                    dxpart = 0.4;
                    ddxpart_dVd = ddxpart_dVg = ddxpart_dVb = ddxpart_dVs = 0.0;
                    dsxpart_dVd = dsxpart_dVg = dsxpart_dVb = dsxpart_dVs = 0.0;
                }
                else
                {
                    if (qcheq > 0.0)
                        T0 = this.Tconst * qdef * ScalingFactor;
                    else
                        T0 = -this.Tconst * qdef * ScalingFactor;
                    ggtg = this._gtg = T0 * this._cqgb;
                    ggtd = this._gtd = T0 * this._cqdb;
                    ggts = this._gts = T0 * this._cqsb;
                    ggtb = this._gtb = T0 * this._cqbb;
                    gqdef = ScalingFactor * ag0;

                    gcqgb = this._cqgb * ag0;
                    gcqdb = this._cqdb * ag0;
                    gcqsb = this._cqsb * ag0;
                    gcqbb = this._cqbb * ag0;

                    gcggb = (cgdo + cgso + Param.BSIM3cgbo) * ag0;
                    gcgdb = -cgdo * ag0;
                    gcgsb = -cgso * ag0;

                    gcdgb = -cgdo * ag0;
                    gcddb = (this._capbd + cgdo) * ag0;
                    gcdsb = 0.0;

                    gcsgb = -cgso * ag0;
                    gcsdb = 0.0;
                    gcssb = (this._capbs + cgso) * ag0;

                    gcbgb = -Param.BSIM3cgbo * ag0;
                    gcbdb = -this._capbd * ag0;
                    gcbsb = -this._capbs * ag0;

                    CoxWL = ModelTemperature.Cox * Param.BSIM3weffCV
                          * Param.BSIM3leffCV;
                    if (Math.Abs(qcheq) <= 1.0e-5 * CoxWL)
                    {
                        if (ModelParameters.Xpart.Value < 0.5)
                        {
                            dxpart = 0.4;
                        }
                        else if (ModelParameters.Xpart.Value > 0.5)
                        {
                            dxpart = 0.0;
                        }
                        else
                        {
                            dxpart = 0.5;
                        }
                        ddxpart_dVd = ddxpart_dVg = ddxpart_dVb
                                    = ddxpart_dVs = 0.0;
                    }
                    else
                    {
                        dxpart = qdrn / qcheq;
                        Cdd = this._cddb;
                        Csd = -(this._cgdb + this._cddb
                            + this._cbdb);
                        ddxpart_dVd = (Cdd - dxpart * (Cdd + Csd)) / qcheq;
                        Cdg = this._cdgb;
                        Csg = -(this._cggb + this._cdgb
                            + this._cbgb);
                        ddxpart_dVg = (Cdg - dxpart * (Cdg + Csg)) / qcheq;

                        Cds = this._cdsb;
                        Css = -(this._cgsb + this._cdsb
                            + this._cbsb);
                        ddxpart_dVs = (Cds - dxpart * (Cds + Css)) / qcheq;

                        ddxpart_dVb = -(ddxpart_dVd + ddxpart_dVg + ddxpart_dVs);
                    }
                    sxpart = 1.0 - dxpart;
                    dsxpart_dVd = -ddxpart_dVd;
                    dsxpart_dVg = -ddxpart_dVg;
                    dsxpart_dVs = -ddxpart_dVs;
                    dsxpart_dVb = -(dsxpart_dVd + dsxpart_dVg + dsxpart_dVs);

                    qgd = qgdo;
                    qgs = qgso;
                    qgb = Param.BSIM3cgbo * vgb;
                    qgate = qgd + qgs + qgb;
                    qbulk = -qgb;
                    qdrn = -qgd;
                    qsrc = -(qgate + qbulk + qdrn);
                }
            }
            else
            {
                if (Parameters.NqsMod.Value == 0)
                {
                    gcggb = (this._cggb + cgdo + cgso
                          + Param.BSIM3cgbo) * ag0;
                    gcgdb = (this._cgsb - cgdo) * ag0;
                    gcgsb = (this._cgdb - cgso) * ag0;

                    gcdgb = -(this._cggb + this._cbgb
                          + this._cdgb + cgdo) * ag0;
                    gcddb = (this._capbd + cgdo - (this._cgsb
                          + this._cbsb + this._cdsb)) * ag0;
                    gcdsb = -(this._cgdb + this._cbdb
                          + this._cddb) * ag0;

                    gcsgb = (this._cdgb - cgso) * ag0;
                    gcsdb = this._cdsb * ag0;
                    gcssb = (this._cddb + this._capbs + cgso) * ag0;

                    gcbgb = (this._cbgb - Param.BSIM3cgbo) * ag0;
                    gcbdb = (this._cbsb - this._capbd) * ag0;
                    gcbsb = (this._cbdb - this._capbs) * ag0;

                    qgd = qgdo;
                    qgs = qgso;
                    qgb = Param.BSIM3cgbo * vgb;
                    qgate += qgd + qgs + qgb;
                    qbulk -= qgb;
                    qsrc = qdrn - qgs;
                    qdrn = -(qgate + qbulk + qsrc);

                    ggtg = ggtd = ggtb = ggts = 0.0;
                    sxpart = 0.4;
                    dxpart = 0.6;
                    ddxpart_dVd = ddxpart_dVg = ddxpart_dVb = ddxpart_dVs = 0.0;
                    dsxpart_dVd = dsxpart_dVg = dsxpart_dVb = dsxpart_dVs = 0.0;
                }
                else
                {
                    if (qcheq > 0.0)
                        T0 = this.Tconst * qdef * ScalingFactor;
                    else
                        T0 = -this.Tconst * qdef * ScalingFactor;
                    ggtg = this._gtg = T0 * this._cqgb;
                    ggts = this._gtd = T0 * this._cqdb;
                    ggtd = this._gts = T0 * this._cqsb;
                    ggtb = this._gtb = T0 * this._cqbb;
                    gqdef = ScalingFactor * ag0;

                    gcqgb = this._cqgb * ag0;
                    gcqdb = this._cqsb * ag0;
                    gcqsb = this._cqdb * ag0;
                    gcqbb = this._cqbb * ag0;

                    gcggb = (cgdo + cgso + Param.BSIM3cgbo) * ag0;
                    gcgdb = -cgdo * ag0;
                    gcgsb = -cgso * ag0;

                    gcdgb = -cgdo * ag0;
                    gcddb = (this._capbd + cgdo) * ag0;
                    gcdsb = 0.0;

                    gcsgb = -cgso * ag0;
                    gcsdb = 0.0;
                    gcssb = (this._capbs + cgso) * ag0;

                    gcbgb = -Param.BSIM3cgbo * ag0;
                    gcbdb = -this._capbd * ag0;
                    gcbsb = -this._capbs * ag0;

                    CoxWL = ModelTemperature.Cox * Param.BSIM3weffCV
                          * Param.BSIM3leffCV;
                    if (Math.Abs(qcheq) <= 1.0e-5 * CoxWL)
                    {
                        if (ModelParameters.Xpart.Value < 0.5)
                        {
                            sxpart = 0.4;
                        }
                        else if (ModelParameters.Xpart.Value > 0.5)
                        {
                            sxpart = 0.0;
                        }
                        else
                        {
                            sxpart = 0.5;
                        }
                        dsxpart_dVd = dsxpart_dVg = dsxpart_dVb
                                    = dsxpart_dVs = 0.0;
                    }
                    else
                    {
                        sxpart = qdrn / qcheq;
                        Css = this._cddb;
                        Cds = -(this._cgdb + this._cddb
                            + this._cbdb);
                        dsxpart_dVs = (Css - sxpart * (Css + Cds)) / qcheq;
                        Csg = this._cdgb;
                        Cdg = -(this._cggb + this._cdgb
                            + this._cbgb);
                        dsxpart_dVg = (Csg - sxpart * (Csg + Cdg)) / qcheq;

                        Csd = this._cdsb;
                        Cdd = -(this._cgsb + this._cdsb
                            + this._cbsb);
                        dsxpart_dVd = (Csd - sxpart * (Csd + Cdd)) / qcheq;

                        dsxpart_dVb = -(dsxpart_dVd + dsxpart_dVg + dsxpart_dVs);
                    }
                    dxpart = 1.0 - sxpart;
                    ddxpart_dVd = -dsxpart_dVd;
                    ddxpart_dVg = -dsxpart_dVg;
                    ddxpart_dVs = -dsxpart_dVs;
                    ddxpart_dVb = -(ddxpart_dVd + ddxpart_dVg + ddxpart_dVs);

                    qgd = qgdo;
                    qgs = qgso;
                    qgb = Param.BSIM3cgbo * vgb;
                    qgate = qgd + qgs + qgb;
                    qbulk = -qgb;
                    qsrc = -qgs;
                    qdrn = -(qgate + qbulk + qsrc);
                }
            }

            cqdef = cqcheq = 0.0;

            _qg.Value = qgate;
            _qd.Value = qdrn - this._qbd;
            _qb.Value = qbulk + this._qbd + this._qbs;

            if (Parameters.NqsMod.Value != 0)
            {
                _qcdump.Value = qdef * ScalingFactor;
                _qcheq.Value = qcheq;
            }

            /* store small signal parameters */
            if (InitializeSmallSignal)
                return;
            if (!ChargeComputationNeeded)
                goto line850;

            _qb.Derive();
            _qg.Derive();
            _qd.Derive();
            if (Parameters.NqsMod.Value != 0)
            {
                _qcdump.Derive();
                _qcheq.Derive();
            }

            goto line860;

        line850:
            /* initialize to zero charge conductance and current */
            ceqqg = ceqqb = ceqqd = 0.0;
            cqcheq = cqdef = 0.0;

            gcdgb = gcddb = gcdsb = 0.0;
            gcsgb = gcsdb = gcssb = 0.0;
            gcggb = gcgdb = gcgsb = 0.0;
            gcbgb = gcbdb = gcbsb = 0.0;

            gqdef = gcqgb = gcqdb = gcqsb = gcqbb = 0.0;
            ggtg = ggtd = ggtb = ggts = 0.0;
            sxpart = (1.0 - (dxpart = (this._mode > 0) ? 0.4 : 0.6));
            ddxpart_dVd = ddxpart_dVg = ddxpart_dVb = ddxpart_dVs = 0.0;
            dsxpart_dVd = dsxpart_dVg = dsxpart_dVb = dsxpart_dVs = 0.0;

            if (Parameters.NqsMod.Value != 0)
                this._gtau = 16.0 * this.U0temp * ModelTemperature.Vtm
                                / Param.BSIM3leffCV / Param.BSIM3leffCV
                                * ScalingFactor;
            else
                this._gtau = 0.0;

            goto line900;

        line860:
            /* evaluate equivalent charge current */

            cqgate = _qg.Derivative;
            cqbulk = _qb.Derivative;
            cqdrn = _qd.Derivative;

            ceqqg = cqgate - gcggb * vgb + gcgdb * vbd + gcgsb * vbs;
            ceqqb = cqbulk - gcbgb * vgb + gcbdb * vbd + gcbsb * vbs;
            ceqqd = cqdrn - gcdgb * vgb + gcddb * vbd + gcdsb * vbs;

            if (Parameters.NqsMod.Value != 0)
            {
                T0 = ggtg * vgb - ggtd * vbd - ggts * vbs;
                ceqqg += T0;
                T1 = qdef * this._gtau;
                ceqqd -= dxpart * T0 + T1 * (ddxpart_dVg * vgb - ddxpart_dVd
                      * vbd - ddxpart_dVs * vbs);
                cqdef = _qcdump.Derivative - gqdef * qdef;
                cqcheq = _qcheq.Derivative - (gcqgb * vgb - gcqdb * vbd - gcqsb * vbs) + T0;
            }

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
                cdreq = ModelParameters.B3Type * (cdrain - this._gds * vds
                      - Gm * vgs - Gmbs * vbs);

                ceqbd = -ModelParameters.B3Type * (this._csub
                      - this._gbds * vds - this._gbgs * vgs
                      - this._gbbs * vbs);
                ceqbs = 0.0;

                gbbdp = -this._gbds;
                gbbsp = (this._gbds + this._gbgs + this._gbbs);

                gbdpg = this._gbgs;
                gbdpdp = this._gbds;
                gbdpb = this._gbbs;
                gbdpsp = -(gbdpg + gbdpdp + gbdpb);

                gbspg = 0.0;
                gbspdp = 0.0;
                gbspb = 0.0;
                gbspsp = 0.0;
            }
            else
            {
                Gm = -this._gm;
                Gmbs = -this._gmbs;
                FwdSum = 0.0;
                RevSum = -(Gm + Gmbs);
                cdreq = -ModelParameters.B3Type * (cdrain + this._gds * vds
                      + Gm * vgd + Gmbs * vbd);

                ceqbs = -ModelParameters.B3Type * (this._csub
                      + this._gbds * vds - this._gbgs * vgd
                      - this._gbbs * vbd);
                ceqbd = 0.0;

                gbbsp = -this._gbds;
                gbbdp = (this._gbds + this._gbgs + this._gbbs);

                gbdpg = 0.0;
                gbdpsp = 0.0;
                gbdpb = 0.0;
                gbdpdp = 0.0;

                gbspg = this._gbgs;
                gbspsp = this._gbds;
                gbspb = this._gbbs;
                gbspdp = -(gbspg + gbspsp + gbspb);
            }

            if (ModelParameters.B3Type > 0)
            {
                ceqbs += (this._cbs - this._gbs * vbs);
                ceqbd += (this._cbd - this._gbd * vbd);
                /*
                ceqqg = ceqqg;
                ceqqb = ceqqb;
                ceqqd = ceqqd;
                cqdef = cqdef;
                cqcheq = cqcheq;
                */
            }
            else
            {
                ceqbs -= (this._cbs - this._gbs * vbs);
                ceqbd -= (this._cbd - this._gbd * vbd);
                ceqqg = -ceqqg;
                ceqqb = -ceqqb;
                ceqqd = -ceqqd;
                cqdef = -cqdef;
                cqcheq = -cqcheq;
            }

            m = Parameters.M;

            _gPtr.Value -= m * ceqqg;
            _bPtr.Value -= m * (ceqbs + ceqbd + ceqqb);
            _dpPtr.Value += m * (ceqbd - cdreq - ceqqd);
            _spPtr.Value += m * (cdreq + ceqbs + ceqqg
                                                     + ceqqb + ceqqd);
            if (Parameters.NqsMod.Value != 0)
                _qPtr.Value += m * (cqcheq - cqdef);

            /*
             *  load y matrix
             */

            T1 = qdef * this._gtau;

            _ddPtr.Value += m * this.DrainConductance;
            _ggPtr.Value += m * (gcggb - ggtg);
            _ssPtr.Value += m * this.SourceConductance;
            _bbPtr.Value += m * (this._gbd + this._gbs
                                 - gcbgb - gcbdb - gcbsb - this._gbbs);
            _dpdpPtr.Value += m * (this.DrainConductance
                                   + this._gds + this._gbd
                                   + RevSum + gcddb + dxpart * ggtd
                                   + T1 * ddxpart_dVd + gbdpdp);
            _spspPtr.Value += m * (this.SourceConductance
                                   + this._gds + this._gbs
                                   + FwdSum + gcssb + sxpart * ggts
                                   + T1 * dsxpart_dVs + gbspsp);
            _ddpPtr.Value -= m * this.DrainConductance;
            _gbPtr.Value -= m * (gcggb + gcgdb + gcgsb + ggtb);
            _gdpPtr.Value += m * (gcgdb - ggtd);
            _gspPtr.Value += m * (gcgsb - ggts);
            _sspPtr.Value -= m * this.SourceConductance;
            _bgPtr.Value += m * (gcbgb - this._gbgs);
            _bdpPtr.Value += m * (gcbdb - this._gbd + gbbdp);
            _bspPtr.Value += m * (gcbsb - this._gbs + gbbsp);
            _dpdPtr.Value -= m * this.DrainConductance;
            _dpgPtr.Value += m * (Gm + gcdgb + dxpart * ggtg
                                 + T1 * ddxpart_dVg + gbdpg);
            _dpbPtr.Value -= m * (this._gbd - Gmbs + gcdgb + gcddb
                                  + gcdsb - dxpart * ggtb
                                  - T1 * ddxpart_dVb - gbdpb);
            _dpspPtr.Value -= m * (this._gds + FwdSum - gcdsb
                                   - dxpart * ggts - T1 * ddxpart_dVs - gbdpsp);
            _spgPtr.Value += m * (gcsgb - Gm + sxpart * ggtg
                                  + T1 * dsxpart_dVg + gbspg);
            _spsPtr.Value -= m * this.SourceConductance;
            _spbPtr.Value -= m * (this._gbs + Gmbs + gcsgb + gcsdb
                                  + gcssb - sxpart * ggtb
                                  - T1 * dsxpart_dVb - gbspb);
            _spdpPtr.Value -= m * (this._gds + RevSum - gcsdb
                                   - sxpart * ggtd - T1 * dsxpart_dVd - gbspdp);

            if (Parameters.NqsMod.Value != 0)
            {
                _qqPtr.Value += m * (gqdef + this._gtau);

                _dpqPtr.Value += m * (dxpart * this._gtau);
                _spqPtr.Value += m * (sxpart * this._gtau);
                _gqPtr.Value -= m * this._gtau;

                _qgPtr.Value += m * (ggtg - gcqgb);
                _qdpPtr.Value += m * (ggtd - gcqdb);
                _qspPtr.Value += m * (ggts - gcqsb);
                _qbPtr.Value += m * (ggtb - gcqbb);
            }
        }
    }
}