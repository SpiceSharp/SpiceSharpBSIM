using SpiceSharp.Behaviors;
using System;
using SpiceSharp.Components;
using SpiceSharp.Attributes;
using SpiceSharp.Simulations;
using SpiceSharp.Algebra;
using SpiceSharp;
using SpiceSharp.Components.Mosfets;
using SpiceSharp.Components.Semiconductors;

namespace SpiceSharpBSIM.Components.Semiconductors.BSIM.BSIM4Behaviors
{
    /// <summary>
    /// Biasing behavior for a <see cref="BSIM4"/>
    /// </summary>
    [BehaviorFor(typeof(BSIM4)), AddBehaviorIfNo(typeof(IBiasingBehavior))]
    [GeneratedParameters]
    public partial class BiasingBehavior : TemperatureBehavior, IBiasingBehavior, ITimeBehavior
    {
        private readonly ITemperatureSimulationState _temperature;
        private readonly IIntegrationMethod _method;
        private readonly IIterationSimulationState _iteration;
        private readonly ITimeSimulationState _time;
        private readonly IDerivative _qcheq, _qcdump, _qb, _qgmid, _qs, _qg, _qd, _qbs, _qbd;

        public const double ScalingFactor = 1.0e-9;
        public const double EPSSI = 1.03594e-10;
        public const double EXPL_THRESHOLD = 100.0;
        public const double MAX_EXPL = 2.688117142e+43;
        public const double MIN_EXPL = 3.720075976e-44;
        public const double MM = 3.0; // Smooth factor
        public const double DELTA_1 = 0.02;
        public const double DELTA_2 = 0.02;
        public const double DELTA_3 = 0.02;
        public const double DELTA_4 = 0.02;

        protected bool InitializeSmallSignal { get; set; }
        protected bool InitializeTransient { get; set; }

        protected double _vds, _vgs, _vbs, _vges, _vgms, _vdbs, _vsbs, _vses, _vdes, _qdef, _vbd, _vdbd, _von, _gbs, _cbs,
            _gbd, _cbd, _mode, _thetavth, _nstar, _vgs_eff, _vgd_eff, _dvgs_eff_dvg, _dvgd_eff_dvg, _vgsteff, _grdsw, _abulk,
            _ueff, _esatL, _vdsat, _vdseff, _coxeff, _abovVgst2Vtm, _csub, _gbbs, _gbgs, _gbds, _gds, _gm, _gmbs, _idovVds,
            _gcrg, _gcrgd, _gcrgb, _gcrgg, _gcrgs, _gstot, _gdtot, _gstotd, _gstotg, _gstots, _gstotb, _gdtotd, _gdtotg,
            _gdtots, _gdtotb, _igidl, _ggidld, _ggidlg, _ggidlb, _igisl, _ggisls, _ggislg, _ggislb, _igcs, _gIgcsg, _gIgcsd,
            _gIgcsb, _igcd, _gIgcdg, _gIgcdd, _gIgcdb, _igs, _gIgsg, _gIgss, _igd, _gIgdg, _gIgdd, _igb, _gIgbg, _gIgbd, _gIgbb,
            _gIgbs, _ggidls, _ggisld, _gIgcss, _gIgcds, _cd, _qinv, _noiGd0, _cggb, _cgsb, _cgdb, _cdgb, _cdsb, _cddb, _cbgb, _cbsb,
            _cbdb, _csgb, _cssb, _csdb, _cgbb, _csbb, _cdbb, _cbbb, _cqdb, _cqsb, _cqgb, _cqbb, _gtau, _qgate, _qbulk, _qdrn, _qsrc,
            _qchqs, _taunet, _capbs, _capbd, _qgdo, _qgso, _gtg, _gtd, _gts, _gtb;
        private readonly IVariable<double> _drain, _gate, _source, _bulk, _drainPrime, _gatePrime, _gateMid, _sourcePrime,
            _bulkPrime, _drainBulk, _sourceBulk, _q;
        private readonly Element<double> _dnodeprimePtr, _gnodeprimePtr, _gnodeextPtr, _gnodemidPtr, _snodeprimePtr, _qPtr, _bnodeprimePtr, _dbnodePtr, _sbnodePtr, _snodePtr, _dnodePtr;
        private readonly Element<double> _dpbpPtr, _gpbpPtr, _spbpPtr, _bpdpPtr, _bpgpPtr, _bpspPtr, _bpbpPtr, _ddPtr,
            _gpgpPtr, _ssPtr, _dpdpPtr, _spspPtr, _ddpPtr, _gpdpPtr, _gpspPtr, _sspPtr, _dpspPtr, _dpdPtr,
            _dpgpPtr, _spgpPtr, _spsPtr, _spdpPtr, _qqPtr, _qbpPtr, _qdpPtr, _qspPtr, _qgpPtr, _dpqPtr, _spqPtr,
            _gpqPtr, _gegePtr, _gegpPtr, _gpgePtr, _gedpPtr, _gespPtr, _gebpPtr, _gmdpPtr, _gmgpPtr, _gmgmPtr,
            _gmgePtr, _gmspPtr, _gmbpPtr, _dpgmPtr, _gpgmPtr, _gegmPtr, _spgmPtr, _bpgmPtr, _dpdbPtr, _spsbPtr,
            _dbdpPtr, _dbdbPtr, _dbbpPtr, _dbbPtr, _bpdbPtr, _bpbPtr, _bpsbPtr, _sbspPtr, _sbbpPtr, _sbbPtr,
            _sbsbPtr, _bdbPtr, _bbpPtr, _bsbPtr, _bbPtr, _dgpPtr, _dspPtr, _dbpPtr, _sdpPtr, _sgpPtr, _sbpPtr;

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

            if (this._drainConductance > 0)
                _drainPrime = state.CreatePrivateVariable(Name.Combine("drain"), Units.Volt);
            else
                _drainPrime = _drain;

            if (this._sourceConductance > 0)
                _sourcePrime = state.CreatePrivateVariable(Name.Combine("source"), Units.Volt);
            else
                _sourcePrime = _source;

            if (Parameters.RgateMod.Value > 0)
                _gatePrime = state.CreatePrivateVariable(Name.Combine("gate"), Units.Volt);
            else
                _gatePrime = _gate;

            if (Parameters.RgateMod.Value == 3)
                _gateMid = state.CreatePrivateVariable(Name.Combine("gatemid"), Units.Volt);
            else
                _gateMid = _gate;

            if (Parameters.RbodyMod.Value == 1 || Parameters.RbodyMod.Value == 2)
            {
                _drainBulk = state.CreatePrivateVariable(Name.Combine("dbody"), Units.Volt);
                _bulkPrime = state.CreatePrivateVariable(Name.Combine("body"), Units.Volt);
                _sourceBulk = state.CreatePrivateVariable(Name.Combine("sbody"), Units.Volt);
            }
            else
            {
                _drainBulk = _bulk;
                _bulkPrime = _bulk;
                _sourceBulk = _bulk;
            }

            if (Parameters.TrnqsMod.Value != 0)
                _q = state.CreatePrivateVariable(Name.Combine("charge"), Units.Coulomb);
            else
                _q = state.GetSharedVariable(Constants.Ground);

            int drain = state.Map[_drain];
            int gate = state.Map[_gate];
            int source = state.Map[_source];
            int bulk = state.Map[_bulk];
            int drainPrime = state.Map[_drainPrime];
            int gatePrime = state.Map[_gatePrime];
            int gateMid = state.Map[_gateMid];
            int sourcePrime = state.Map[_sourcePrime];
            int bulkPrime = state.Map[_bulkPrime];
            int q = state.Map[_q];
            int dbulk = state.Map[_drainBulk];
            int sbulk = state.Map[_sourceBulk];

            _dnodeprimePtr = state.Solver.GetElement(drainPrime);
            _snodeprimePtr = state.Solver.GetElement(sourcePrime);
            _gnodeprimePtr = state.Solver.GetElement(gatePrime);
            _qPtr = state.Solver.GetElement(q);
            _bnodeprimePtr = state.Solver.GetElement(bulkPrime);

            _dpbpPtr = state.Solver.GetElement(new MatrixLocation(drainPrime, bulkPrime));
            _gpbpPtr = state.Solver.GetElement(new MatrixLocation(gatePrime, bulkPrime));
            _spbpPtr = state.Solver.GetElement(new MatrixLocation(sourcePrime, bulkPrime));
            _bpdpPtr = state.Solver.GetElement(new MatrixLocation(bulkPrime, drainPrime));
            _bpgpPtr = state.Solver.GetElement(new MatrixLocation(bulkPrime, gatePrime));
            _bpspPtr = state.Solver.GetElement(new MatrixLocation(bulkPrime, sourcePrime));
            _bpbpPtr = state.Solver.GetElement(new MatrixLocation(bulkPrime, bulkPrime));
            _ddPtr = state.Solver.GetElement(new MatrixLocation(drain, drain));
            _gpgpPtr = state.Solver.GetElement(new MatrixLocation(gatePrime, gatePrime));
            _ssPtr = state.Solver.GetElement(new MatrixLocation(source, source));
            _dpdpPtr = state.Solver.GetElement(new MatrixLocation(drainPrime, drainPrime));
            _spspPtr = state.Solver.GetElement(new MatrixLocation(sourcePrime, sourcePrime));
            _ddpPtr = state.Solver.GetElement(new MatrixLocation(drain, drainPrime));
            _gpdpPtr = state.Solver.GetElement(new MatrixLocation(gatePrime, drainPrime));
            _gpspPtr = state.Solver.GetElement(new MatrixLocation(gatePrime, sourcePrime));
            _sspPtr = state.Solver.GetElement(new MatrixLocation(source, sourcePrime));
            _dpspPtr = state.Solver.GetElement(new MatrixLocation(drainPrime, sourcePrime));
            _dpdPtr = state.Solver.GetElement(new MatrixLocation(drainPrime, drain));
            _dpgpPtr = state.Solver.GetElement(new MatrixLocation(drainPrime, gatePrime));
            _spgpPtr = state.Solver.GetElement(new MatrixLocation(sourcePrime, gatePrime));
            _spsPtr = state.Solver.GetElement(new MatrixLocation(sourcePrime, source));
            _spdpPtr = state.Solver.GetElement(new MatrixLocation(sourcePrime, drainPrime));
            _qqPtr = state.Solver.GetElement(new MatrixLocation(q, q));
            _qbpPtr = state.Solver.GetElement(new MatrixLocation(q, bulkPrime));
            _qdpPtr = state.Solver.GetElement(new MatrixLocation(q, drainPrime));
            _qspPtr = state.Solver.GetElement(new MatrixLocation(q, sourcePrime));
            _qgpPtr = state.Solver.GetElement(new MatrixLocation(q, gatePrime));
            _dpqPtr = state.Solver.GetElement(new MatrixLocation(drainPrime, q));
            _spqPtr = state.Solver.GetElement(new MatrixLocation(sourcePrime, q));
            _gpqPtr = state.Solver.GetElement(new MatrixLocation(gatePrime, q));
            if (Parameters.RgateMod.Value != 0)
            {
                _gnodeextPtr = state.Solver.GetElement(gate);
                _gnodemidPtr = state.Solver.GetElement(gateMid);
                _gegePtr = state.Solver.GetElement(new MatrixLocation(gate, gate));
                _gegpPtr = state.Solver.GetElement(new MatrixLocation(gate, gatePrime));
                _gpgePtr = state.Solver.GetElement(new MatrixLocation(gatePrime, gate));
                _gedpPtr = state.Solver.GetElement(new MatrixLocation(gate, drainPrime));
                _gespPtr = state.Solver.GetElement(new MatrixLocation(gate, sourcePrime));
                _gebpPtr = state.Solver.GetElement(new MatrixLocation(gate, bulkPrime));
                _gmdpPtr = state.Solver.GetElement(new MatrixLocation(gateMid, drainPrime));
                _gmgpPtr = state.Solver.GetElement(new MatrixLocation(gateMid, gatePrime));
                _gmgmPtr = state.Solver.GetElement(new MatrixLocation(gateMid, gateMid));
                _gmgePtr = state.Solver.GetElement(new MatrixLocation(gateMid, gate));
                _gmspPtr = state.Solver.GetElement(new MatrixLocation(gateMid, sourcePrime));
                _gmbpPtr = state.Solver.GetElement(new MatrixLocation(gateMid, bulkPrime));
                _dpgmPtr = state.Solver.GetElement(new MatrixLocation(drainPrime, gateMid));
                _gpgmPtr = state.Solver.GetElement(new MatrixLocation(gatePrime, gateMid));
                _gegmPtr = state.Solver.GetElement(new MatrixLocation(gate, gateMid));
                _spgmPtr = state.Solver.GetElement(new MatrixLocation(sourcePrime, gateMid));
                _bpgmPtr = state.Solver.GetElement(new MatrixLocation(bulkPrime, gateMid));
            }
            if ((Parameters.RbodyMod.Value == 1) || (Parameters.RbodyMod.Value == 2))
            {
                _dbnodePtr = state.Solver.GetElement(dbulk);
                _sbnodePtr = state.Solver.GetElement(sbulk);
                _dpdbPtr = state.Solver.GetElement(new MatrixLocation(drainPrime, dbulk));
                _spsbPtr = state.Solver.GetElement(new MatrixLocation(sourcePrime, sbulk));
                _dbdpPtr = state.Solver.GetElement(new MatrixLocation(dbulk, drainPrime));
                _dbdbPtr = state.Solver.GetElement(new MatrixLocation(dbulk, dbulk));
                _dbbpPtr = state.Solver.GetElement(new MatrixLocation(dbulk, bulkPrime));
                _dbbPtr = state.Solver.GetElement(new MatrixLocation(dbulk, bulk));
                _bpdbPtr = state.Solver.GetElement(new MatrixLocation(bulkPrime, dbulk));
                _bpbPtr = state.Solver.GetElement(new MatrixLocation(bulkPrime, bulk));
                _bpsbPtr = state.Solver.GetElement(new MatrixLocation(bulkPrime, sbulk));
                _sbspPtr = state.Solver.GetElement(new MatrixLocation(sbulk, sourcePrime));
                _sbbpPtr = state.Solver.GetElement(new MatrixLocation(sbulk, bulkPrime));
                _sbbPtr = state.Solver.GetElement(new MatrixLocation(sbulk, bulk));
                _sbsbPtr = state.Solver.GetElement(new MatrixLocation(sbulk, sbulk));
                _bdbPtr = state.Solver.GetElement(new MatrixLocation(bulk, dbulk));
                _bbpPtr = state.Solver.GetElement(new MatrixLocation(bulk, bulkPrime));
                _bsbPtr = state.Solver.GetElement(new MatrixLocation(bulk, sbulk));
                _bbPtr = state.Solver.GetElement(new MatrixLocation(bulk, bulk));
            }
            if (ModelParameters.RdsMod.Value != 0)
            {
                _snodePtr = state.Solver.GetElement(source);
                _dnodePtr = state.Solver.GetElement(drain);
                _dgpPtr = state.Solver.GetElement(new MatrixLocation(drain, gatePrime));
                _dspPtr = state.Solver.GetElement(new MatrixLocation(drain, sourcePrime));
                _dbpPtr = state.Solver.GetElement(new MatrixLocation(drain, bulkPrime));
                _sdpPtr = state.Solver.GetElement(new MatrixLocation(source, drainPrime));
                _sgpPtr = state.Solver.GetElement(new MatrixLocation(source, gatePrime));
                _sbpPtr = state.Solver.GetElement(new MatrixLocation(source, bulkPrime));
            }

            if (_method != null)
            {
                _qcheq = _method.CreateDerivative();
                _qcdump = _method.CreateDerivative();
                _qb = _method.CreateDerivative();
                _qgmid = _method.CreateDerivative();
                _qg = _method.CreateDerivative();
                _qs = _method.CreateDerivative();
                _qd = _method.CreateDerivative();
                _qbs = _method.CreateDerivative();
                _qbd = _method.CreateDerivative();
            }
            else
            {
                _qcheq = new SimpleDerivative();
                _qcdump = new SimpleDerivative();
                _qb = new SimpleDerivative();
                _qgmid = new SimpleDerivative();
                _qg = new SimpleDerivative();
                _qs = new SimpleDerivative();
                _qd = new SimpleDerivative();
                _qbs = new SimpleDerivative();
                _qbd = new SimpleDerivative();
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
            double ceqgstot, dgstot_dvd, dgstot_dvg, dgstot_dvs, dgstot_dvb;
            double ceqgdtot, dgdtot_dvd, dgdtot_dvg, dgdtot_dvs, dgdtot_dvb;
            double gstot, gstotd, gstotg, gstots, gstotb, gspr, Rs, Rd;
            double gdtot, gdtotd, gdtotg, gdtots, gdtotb, gdpr;
            double vgs_eff, vgd_eff, dvgs_eff_dvg, dvgd_eff_dvg;
            double dRs_dvg, dRd_dvg, dRs_dvb, dRd_dvb;
            double dT0_dvg, dT1_dvb, dT3_dvg, dT3_dvb;
            double vses, vdes, vdedo, delvses, delvded, delvdes;

            double geltd, gcrg, gcrgg, gcrgd, gcrgs, gcrgb, ceqgcrg;
            double vges, vgms, vgedo, vgmdo, vged, vgmd, delvged, delvgmd;
            double delvges, delvgms, vgmb;
            double gcgmgmb = 0.0, gcgmdb = 0.0, gcgmsb = 0.0, gcdgmb, gcsgmb;
            double gcgmbb = 0.0, gcbgmb, qgmb, qgmid = 0.0, ceqqgmid;

            double vbd, vbs, vds, vgb, vgd, vgs, vgdo;

            double vdbs, vdbd, vsbs, vsbdo, vsbd;
            double delvdbs, delvdbd, delvsbs;
            double delvbd_jct, delvbs_jct, vbs_jct, vbd_jct;

            double SourceSatCurrent, DrainSatCurrent;
            double ag0, qgb, von, VgstNVt, ExpVgst;
            double ceqqb, ceqqd, ceqqg, ceqqjd = 0.0, ceqqjs = 0.0;
            double cdrain, ceqdrn, ceqbd, ceqbs, ceqjd, ceqjs, gjbd, gjbs;
            double czbd, czbdsw, czbdswg, czbs, czbssw, czbsswg, evbd, evbs, arg, sarg;
            double delvbd, delvbs, delvds, delvgd, delvgs;
            double Vfbeff, dVfbeff_dVg, dVfbeff_dVb, V3, V4;
            double gcbdb, gcbgb, gcbsb, gcddb, gcdgb, gcdsb, gcgdb, gcggb, gcgsb, gcsdb;
            double gcgbb, gcdbb, gcsbb, gcbbb;
            double gcdbdb, gcsbsb;
            double gcsgb, gcssb, MJD, MJSWD, MJSWGD, MJS, MJSWS, MJSWGS;
            double qgate = 0.0, qbulk = 0.0, qdrn = 0.0, qsrc, cqgate, cqbody, cqdrn;
            double Vdb, Vds, Vgs, Vbs, Gmbs, FwdSum, RevSum;
            double Igidl, Ggidld, Ggidlg, Ggidlb;
            double Voxacc = 0.0, dVoxacc_dVg = 0.0, dVoxacc_dVb = 0.0;
            double Voxdepinv = 0.0, dVoxdepinv_dVg = 0.0, dVoxdepinv_dVd = 0.0, dVoxdepinv_dVb = 0.0;
            double VxNVt = 0.0, ExpVxNVt, Vaux = 0.0, dVaux_dVg = 0.0, dVaux_dVd = 0.0, dVaux_dVb = 0.0;
            double Igc, dIgc_dVg, dIgc_dVd, dIgc_dVb;
            double Igcs, dIgcs_dVg, dIgcs_dVd, dIgcs_dVb;
            double Igcd, dIgcd_dVg, dIgcd_dVd, dIgcd_dVb;
            double Igs, dIgs_dVg, dIgs_dVs, Igd, dIgd_dVg, dIgd_dVd;
            double Igbacc, dIgbacc_dVg, dIgbacc_dVb;
            double Igbinv, dIgbinv_dVg, dIgbinv_dVd, dIgbinv_dVb;
            double Pigcd, dPigcd_dVg, dPigcd_dVd, dPigcd_dVb;
            double Istoteq, gIstotg, gIstotd, gIstots, gIstotb;
            double Idtoteq, gIdtotg, gIdtotd, gIdtots, gIdtotb;
            double Ibtoteq, gIbtotg, gIbtotd, gIbtots, gIbtotb;
            double Igtoteq, gIgtotg, gIgtotd, gIgtots, gIgtotb;
            double Vgs_eff, Vfb = 0.0, Vth_NarrowW;
            /* double Vgd_eff, dVgd_eff_dVg;          v4.7.0 */
            double Phis, dPhis_dVb, sqrtPhis, dsqrtPhis_dVb, Vth, dVth_dVb, dVth_dVd;
            double Vgst, dVgst_dVg, dVgst_dVb, dVgs_eff_dVg, Nvtms, Nvtmd;
            double Vtm, Vtm0;
            double n, dn_dVb, dn_dVd, voffcv, noff, dnoff_dVd, dnoff_dVb;
            double V0, CoxWLcen, QovCox, LINK;
            double DeltaPhi, dDeltaPhi_dVg, VgDP, dVgDP_dVg;
            double Cox, Tox, Tcen, dTcen_dVg, dTcen_dVd, dTcen_dVb;
            double Ccen, Coxeff, dCoxeff_dVd, dCoxeff_dVg, dCoxeff_dVb;
            double Denomi, dDenomi_dVg, dDenomi_dVd, dDenomi_dVb;
            double ueff, dueff_dVg, dueff_dVd, dueff_dVb;
            double Esat, Vdsat;
            double EsatL, dEsatL_dVg, dEsatL_dVd, dEsatL_dVb;
            double dVdsat_dVg, dVdsat_dVb, dVdsat_dVd, Vasat, dAlphaz_dVg, dAlphaz_dVb;
            double dVasat_dVg, dVasat_dVb, dVasat_dVd, Va, dVa_dVd, dVa_dVg, dVa_dVb;
            double Vbseff, dVbseff_dVb, VbseffCV, dVbseffCV_dVb;
            double VgsteffVth, dT11_dVg;
            double Arg1, One_Third_CoxWL, Two_Third_CoxWL, Alphaz, CoxWL;
            double T0 = 0.0, dT0_dVg, dT0_dVd, dT0_dVb;
            double T1, dT1_dVg, dT1_dVd, dT1_dVb;
            double T2, dT2_dVg, dT2_dVd, dT2_dVb;
            double T3, dT3_dVg, dT3_dVd, dT3_dVb;
            double T4, dT4_dVd, dT4_dVb;
            double T5, dT5_dVg, dT5_dVd, dT5_dVb;
            double T6, dT6_dVg, dT6_dVd, dT6_dVb;
            double T7, dT7_dVg, dT7_dVd, dT7_dVb;
            double T8, dT8_dVg, dT8_dVd, dT8_dVb;
            double T9, dT9_dVg, dT9_dVd, dT9_dVb;
            double T10, dT10_dVg, dT10_dVb, dT10_dVd;
            double T11, T12, T13, T14;
            double tmp, Abulk, dAbulk_dVb, Abulk0, dAbulk0_dVb;
            double Cclm, dCclm_dVg, dCclm_dVd, dCclm_dVb;
            double FP, dFP_dVg, PvagTerm, dPvagTerm_dVg, dPvagTerm_dVd, dPvagTerm_dVb;
            double VADITS, dVADITS_dVg, dVADITS_dVd;
            double Lpe_Vb, dDITS_Sft_dVb, dDITS_Sft_dVd;
            double DITS_Sft2, dDITS_Sft2_dVd;        /* v4.7 New DITS */
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
            double Ids, Gm, Gds, Gmb, devbs_dvb, devbd_dvb;
            double Isub, Gbd, Gbg, Gbb;
            double VASCBE, dVASCBE_dVg, dVASCBE_dVd, dVASCBE_dVb;
            double CoxeffWovL;
            double Rds, dRds_dVg, dRds_dVb, WVCox, WVCoxRds;
            double Vgst2Vtm, VdsatCV;
            double Leff, Weff, dWeff_dVg, dWeff_dVb;
            double AbulkCV, dAbulkCV_dVb;
            double qcheq, qdef, gqdef = 0.0, cqdef = 0.0, cqcheq = 0.0;
            double gcqdb = 0.0, gcqsb = 0.0, gcqgb = 0.0, gcqbb = 0.0;
            double dxpart, sxpart, ggtg, ggtd, ggts, ggtb;
            double ddxpart_dVd, ddxpart_dVg, ddxpart_dVb, ddxpart_dVs;
            double dsxpart_dVd, dsxpart_dVg, dsxpart_dVb, dsxpart_dVs;
            double gbspsp, gbbdp, gbbsp, gbspg, gbspb, gbspdp;
            double gbdpdp, gbdpg, gbdpb, gbdpsp;
            double qgdo, qgso, cgdo, cgso;
            double Cgg, Cgd, Cgb, Cdg, Cdd, Cds;
            double Csg, Csd, Css, Csb, Cbg, Cbd, Cbb;
            double Cgg1, Cgd1, Cgb1, Cbg1, Cbb1, Cbd1, Qac0, Qsub0;
            double dQac0_dVg, dQac0_dVb, dQsub0_dVg, dQsub0_dVd, dQsub0_dVb;
            double ggidld, ggidlg, ggidlb, ggislg, ggislb, ggisls;
            double Igisl, Ggislg, Ggislb, Ggisls;
            double Nvtmrss, Nvtmrssws, Nvtmrsswgs;
            double Nvtmrsd, Nvtmrsswd, Nvtmrsswgd;

            double vs, Fsevl, dvs_dVg, dvs_dVd, dvs_dVb, dFsevl_dVg, dFsevl_dVd, dFsevl_dVb;
            double vgdx, vgsx, epssub, toxe, epsrox;
            bool ChargeComputationNeeded, Check, Check1, Check2;

            double m;


            ChargeComputationNeeded = _time != null || InitializeSmallSignal || InitializeTransient;

            Check = Check1 = Check2 = true;

            if (InitializeSmallSignal || InitializeTransient)
            {
                vds = this._vds;
                vgs = this._vgs;
                vbs = this._vbs;
                vges = this._vges;
                vgms = this._vgms;
                vdbs = this._vdbs;
                vsbs = this._vsbs;
                vses = this._vses;
                vdes = this._vdes;

                qdef = this._qdef;
            }
            else if (_iteration.Mode == IterationModes.Junction && !Parameters.Off)
            {
                vds = ModelParameters.Type * Parameters.IcVDS;
                vgs = vges = vgms = ModelParameters.Type * Parameters.IcVGS;
                vbs = vdbs = vsbs = ModelParameters.Type * Parameters.IcVBS;
                if (vds > 0.0)
                {
                    vdes = vds + 0.01;
                    vses = -0.01;
                }
                else if (vds < 0.0)
                {
                    vdes = vds - 0.01;
                    vses = 0.01;
                }
                else
                    vdes = vses = 0.0;

                qdef = 0.0;

                if (vds.Equals(0.0) && vgs.Equals(0.0) && vbs.Equals(0.0))
                {
                    vds = 0.1;
                    vdes = 0.11;
                    vses = -0.01;
                    vgs = vges = vgms = ModelParameters.Type
                                      * this._vth0 + 0.1;
                    vbs = vdbs = vsbs = 0.0;
                }
            }
            else if ((_iteration.Mode == IterationModes.Junction || _iteration.Mode == IterationModes.Fix) && (Parameters.Off))
            {
                vds = vgs = vbs = vges = vgms = 0.0;
                vdbs = vsbs = vdes = vses = qdef = 0.0;
            }
            else
            {
                vds = ModelParameters.Type * (_drainPrime.Value - _sourcePrime.Value);
                vgs = ModelParameters.Type * (_gatePrime.Value - _sourcePrime.Value);
                vbs = ModelParameters.Type * (_bulkPrime.Value - _sourcePrime.Value);
                vges = ModelParameters.Type * (_gate.Value - _sourcePrime.Value);
                vgms = ModelParameters.Type * (_gateMid.Value - _sourcePrime.Value);
                vdbs = ModelParameters.Type * (_drainBulk.Value - _sourcePrime.Value);
                vsbs = ModelParameters.Type * (_sourceBulk.Value - _sourcePrime.Value);
                vses = ModelParameters.Type * (_source.Value - _sourcePrime.Value);
                vdes = ModelParameters.Type * (_drain.Value - _sourcePrime.Value);
                qdef = ModelParameters.Type * (_q.Value);

                vgdo = this._vgs - this._vds;
                vgedo = this._vges - this._vds;
                vgmdo = this._vgms - this._vds;

                vbd = vbs - vds;
                vdbd = vdbs - vds;
                vgd = vgs - vds;
                vged = vges - vds;
                vgmd = vgms - vds;

                delvbd = vbd - this._vbd;
                delvdbd = vdbd - this._vdbd;
                delvgd = vgd - vgdo;
                delvged = vged - vgedo;
                delvgmd = vgmd - vgmdo;

                delvds = vds - this._vds;
                delvgs = vgs - this._vgs;
                delvges = vges - this._vges;
                delvgms = vgms - this._vgms;
                delvbs = vbs - this._vbs;
                delvdbs = vdbs - this._vdbs;
                delvsbs = vsbs - this._vsbs;

                delvses = vses - (this._vses);
                vdedo = this._vdes
                      - this._vds;
                delvdes = vdes - this._vdes;
                delvded = vdes - vds - vdedo;

                delvbd_jct = (Parameters.RbodyMod.Value == 0) ? delvbd : delvdbd;
                delvbs_jct = (Parameters.RbodyMod.Value == 0) ? delvbs : delvsbs;
                
                von = this._von;
                if (this._vds >= 0.0)
                {
                    vgs = Transistor.LimitFet(vgs, this._vgs, von);
                    vds = vgs - vgd;
                    vds = Transistor.LimitVds(vds, this._vds);
                    vgd = vgs - vds;
                    if (Parameters.RgateMod == 3)
                    {
                        vges = Transistor.LimitFet(vges, this._vges, von);
                        vgms = Transistor.LimitFet(vgms, this._vgms, von);
                        vged = vges - vds;
                        vgmd = vgms - vds;
                    }
                    else if ((Parameters.RgateMod == 1) || (Parameters.RgateMod == 2))
                    {
                        vges = Transistor.LimitFet(vges, this._vges, von);
                        vged = vges - vds;
                    }

                    if (ModelParameters.RdsMod.Value != 0)
                    {
                        vdes = Transistor.LimitVds(vdes, this._vdes);
                        vses = -Transistor.LimitVds(-vses, -(this._vses));
                    }

                }
                else
                {
                    vgd = Transistor.LimitFet(vgd, vgdo, von);
                    vds = vgs - vgd;
                    vds = -Transistor.LimitVds(-vds, -(this._vds));
                    vgs = vgd + vds;

                    if (Parameters.RgateMod == 3)
                    {
                        vged = Transistor.LimitFet(vged, vgedo, von);
                        vges = vged + vds;
                        vgmd = Transistor.LimitFet(vgmd, vgmdo, von);
                        vgms = vgmd + vds;
                    }
                    if ((Parameters.RgateMod == 1) || (Parameters.RgateMod == 2))
                    {
                        vged = Transistor.LimitFet(vged, vgedo, von);
                        vges = vged + vds;
                    }

                    if (ModelParameters.RdsMod.Value != 0)
                    {
                        vdes = -Transistor.LimitVds(-vdes, -(this._vdes));
                        vses = Transistor.LimitVds(vses, this._vses);
                    }
                }

                Check = false;
                Check1 = false;
                Check2 = false;
                if (vds >= 0.0)
                {
                    vbs = Semiconductor.LimitJunction(vbs, this._vbs,
                                    Constants.Vt0, ModelTemperature.Vcrit, ref Check);
                    vbd = vbs - vds;
                    if (Parameters.RbodyMod.Value != 0)
                    {
                        vdbs = Semiconductor.LimitJunction(vdbs, this._vdbs,
                                         Constants.Vt0, ModelTemperature.Vcrit, ref Check1);
                        vdbd = vdbs - vds;
                        vsbs = Semiconductor.LimitJunction(vsbs, this._vsbs,
                                         Constants.Vt0, ModelTemperature.Vcrit, ref Check2);
                        // 20220724 - this would override the Check parameter... - Sven Boulanger
                        if (Check1 || Check2)
                            Check = true;
                    }
                }
                else
                {
                    vbd = Semiconductor.LimitJunction(vbd, this._vbd,
                                    Constants.Vt0, ModelTemperature.Vcrit, ref Check);
                    vbs = vbd + vds;
                    if (Parameters.RbodyMod.Value != 0)
                    {
                        vdbd = Semiconductor.LimitJunction(vdbd, this._vdbd,
                                         Constants.Vt0, ModelTemperature.Vcrit, ref Check1);
                        vdbs = vdbd + vds;
                        vsbdo = this._vsbs
                              - this._vds;
                        vsbd = vsbs - vds;
                        vsbd = Semiconductor.LimitJunction(vsbd, vsbdo, Constants.Vt0, ModelTemperature.Vcrit, ref Check2);
                        vsbs = vsbd + vds;
                        // 20220724 - this would override the Check parameter... - Sven Boulanger
                        if (Check1 || Check2)
                            Check = true;
                    }
                }
            }

            /* Calculate DC currents and their derivatives */
            vbd = vbs - vds;
            vgd = vgs - vds;
            vgb = vgs - vbs;
            vged = vges - vds;
            vgmd = vgms - vds;
            vgmb = vgms - vbs;
            vdbd = vdbs - vds;

            vbs_jct = (Parameters.RbodyMod.Value == 0) ? vbs : vsbs;
            vbd_jct = (Parameters.RbodyMod.Value == 0) ? vbd : vdbd;

            /* Source/drain junction diode DC model begins */
            Nvtms = ModelTemperature.Vtm * ModelParameters.SjctEmissionCoeff;
            /*          if ((this._aseff <= 0.0) && (this._pseff <= 0.0))
                      {   SourceSatCurrent = 1.0e-14;
                      } v4.7 */
            if ((this._aseff <= 0.0) && (this._pseff <= 0.0))
            {
                SourceSatCurrent = 0.0;
            }
            else
            {
                SourceSatCurrent = this._aseff * ModelTemperature.SjctTempSatCurDensity
                                 + this._pseff * ModelTemperature.SjctSidewallTempSatCurDensity
                                 + Param.BSIM4weffCJ * Parameters.Nf
                                 * ModelTemperature.SjctGateSidewallTempSatCurDensity;
            }

            if (SourceSatCurrent <= 0.0)
            {
                this._gbs = _iteration.Gmin;
                this._cbs = this._gbs * vbs_jct;
            }
            else
            {
                switch (ModelParameters.DioMod)
                {
                    case 0:
                        evbs = Math.Exp(vbs_jct / Nvtms);
                        T1 = ModelParameters.Xjbvs * Math.Exp(-(ModelParameters.Bvs + vbs_jct) / Nvtms);
                        /* WDLiu: Magic T1 in this form; different from BSIM4 beta. */
                        this._gbs = SourceSatCurrent * (evbs + T1) / Nvtms + _iteration.Gmin;
                        this._cbs = SourceSatCurrent * (evbs + this._xExpBVS
                                       - T1 - 1.0) + _iteration.Gmin * vbs_jct;
                        break;
                    case 1:
                        T2 = vbs_jct / Nvtms;
                        if (T2 < -EXP_THRESHOLD)
                        {
                            this._gbs = _iteration.Gmin;
                            this._cbs = SourceSatCurrent * (MIN_EXP - 1.0)
                                           + _iteration.Gmin * vbs_jct;
                        }
                        else if (vbs_jct <= this._vjsmFwd)
                        {
                            evbs = Math.Exp(T2);
                            this._gbs = SourceSatCurrent * evbs / Nvtms + _iteration.Gmin;
                            this._cbs = SourceSatCurrent * (evbs - 1.0)
                                           + _iteration.Gmin * vbs_jct;
                        }
                        else
                        {
                            T0 = this._iVjsmFwd / Nvtms;
                            this._gbs = T0 + _iteration.Gmin;
                            this._cbs = this._iVjsmFwd - SourceSatCurrent + T0
                                           * (vbs_jct - this._vjsmFwd) + _iteration.Gmin * vbs_jct;
                        }
                        break;
                    case 2:
                        if (vbs_jct < this._vjsmRev)
                        {
                            T0 = vbs_jct / Nvtms;
                            if (T0 < -EXP_THRESHOLD)
                            {
                                evbs = MIN_EXP;
                                devbs_dvb = 0.0;
                            }
                            else
                            {
                                evbs = Math.Exp(T0);
                                devbs_dvb = evbs / Nvtms;
                            }

                            T1 = evbs - 1.0;
                            T2 = this._iVjsmRev + this._sslpRev
                               * (vbs_jct - this._vjsmRev);
                            this._gbs = devbs_dvb * T2 + T1 * this._sslpRev + _iteration.Gmin;
                            this._cbs = T1 * T2 + _iteration.Gmin * vbs_jct;
                        }
                        else if (vbs_jct <= this._vjsmFwd)
                        {
                            T0 = vbs_jct / Nvtms;
                            if (T0 < -EXP_THRESHOLD)
                            {
                                evbs = MIN_EXP;
                                devbs_dvb = 0.0;
                            }
                            else
                            {
                                evbs = Math.Exp(T0);
                                devbs_dvb = evbs / Nvtms;
                            }

                            T1 = (ModelParameters.Bvs + vbs_jct) / Nvtms;
                            if (T1 > EXP_THRESHOLD)
                            {
                                T2 = MIN_EXP;
                                T3 = 0.0;
                            }
                            else
                            {
                                T2 = Math.Exp(-T1);
                                T3 = -T2 / Nvtms;
                            }
                            this._gbs = SourceSatCurrent * (devbs_dvb - ModelParameters.Xjbvs * T3)
                                           + _iteration.Gmin;
                            this._cbs = SourceSatCurrent * (evbs + this._xExpBVS - 1.0
                                           - ModelParameters.Xjbvs * T2) + _iteration.Gmin * vbs_jct;
                        }
                        else
                        {
                            this._gbs = this._sslpFwd + _iteration.Gmin;
                            this._cbs = this._iVjsmFwd + this._sslpFwd * (vbs_jct
                                           - this._vjsmFwd) + _iteration.Gmin * vbs_jct;
                        }
                        break;
                    default: break;
                }
            }

            Nvtmd = ModelTemperature.Vtm * ModelParameters.DjctEmissionCoeff;
            /*          if ((this._adeff <= 0.0) && (this._pdeff <= 0.0))
                      {   DrainSatCurrent = 1.0e-14;
                      } v4.7 */
            if ((this._adeff <= 0.0) && (this._pdeff <= 0.0))
            {
                DrainSatCurrent = 0.0;
            }
            else
            {
                DrainSatCurrent = this._adeff * ModelTemperature.DjctTempSatCurDensity
                                + this._pdeff * ModelTemperature.DjctSidewallTempSatCurDensity
                                + Param.BSIM4weffCJ * Parameters.Nf
                                * ModelTemperature.DjctGateSidewallTempSatCurDensity;
            }

            if (DrainSatCurrent <= 0.0)
            {
                this._gbd = _iteration.Gmin;
                this._cbd = this._gbd * vbd_jct;
            }
            else
            {
                switch (ModelParameters.DioMod)
                {
                    case 0:
                        evbd = Math.Exp(vbd_jct / Nvtmd);
                        T1 = ModelParameters.Xjbvd * Math.Exp(-(ModelParameters.Bvd + vbd_jct) / Nvtmd);
                        /* WDLiu: Magic T1 in this form; different from BSIM4 beta. */
                        this._gbd = DrainSatCurrent * (evbd + T1) / Nvtmd + _iteration.Gmin;
                        this._cbd = DrainSatCurrent * (evbd + this._xExpBVD
                                       - T1 - 1.0) + _iteration.Gmin * vbd_jct;
                        break;
                    case 1:
                        T2 = vbd_jct / Nvtmd;
                        if (T2 < -EXP_THRESHOLD)
                        {
                            this._gbd = _iteration.Gmin;
                            this._cbd = DrainSatCurrent * (MIN_EXP - 1.0)
                                           + _iteration.Gmin * vbd_jct;
                        }
                        else if (vbd_jct <= this._vjdmFwd)
                        {
                            evbd = Math.Exp(T2);
                            this._gbd = DrainSatCurrent * evbd / Nvtmd + _iteration.Gmin;
                            this._cbd = DrainSatCurrent * (evbd - 1.0)
                                           + _iteration.Gmin * vbd_jct;
                        }
                        else
                        {
                            T0 = this._iVjdmFwd / Nvtmd;
                            this._gbd = T0 + _iteration.Gmin;
                            this._cbd = this._iVjdmFwd - DrainSatCurrent + T0
                                           * (vbd_jct - this._vjdmFwd) + _iteration.Gmin * vbd_jct;
                        }
                        break;
                    case 2:
                        if (vbd_jct < this._vjdmRev)
                        {
                            T0 = vbd_jct / Nvtmd;
                            if (T0 < -EXP_THRESHOLD)
                            {
                                evbd = MIN_EXP;
                                devbd_dvb = 0.0;
                            }
                            else
                            {
                                evbd = Math.Exp(T0);
                                devbd_dvb = evbd / Nvtmd;
                            }

                            T1 = evbd - 1.0;
                            T2 = this._iVjdmRev + this._dslpRev
                               * (vbd_jct - this._vjdmRev);
                            this._gbd = devbd_dvb * T2 + T1 * this._dslpRev + _iteration.Gmin;
                            this._cbd = T1 * T2 + _iteration.Gmin * vbd_jct;
                        }
                        else if (vbd_jct <= this._vjdmFwd)
                        {
                            T0 = vbd_jct / Nvtmd;
                            if (T0 < -EXP_THRESHOLD)
                            {
                                evbd = MIN_EXP;
                                devbd_dvb = 0.0;
                            }
                            else
                            {
                                evbd = Math.Exp(T0);
                                devbd_dvb = evbd / Nvtmd;
                            }

                            T1 = (ModelParameters.Bvd + vbd_jct) / Nvtmd;
                            if (T1 > EXP_THRESHOLD)
                            {
                                T2 = MIN_EXP;
                                T3 = 0.0;
                            }
                            else
                            {
                                T2 = Math.Exp(-T1);
                                T3 = -T2 / Nvtmd;
                            }
                            this._gbd = DrainSatCurrent * (devbd_dvb - ModelParameters.Xjbvd * T3)
                                           + _iteration.Gmin;
                            this._cbd = DrainSatCurrent * (evbd + this._xExpBVD - 1.0
                                           - ModelParameters.Xjbvd * T2) + _iteration.Gmin * vbd_jct;
                        }
                        else
                        {
                            this._gbd = this._dslpFwd + _iteration.Gmin;
                            this._cbd = this._iVjdmFwd + this._dslpFwd * (vbd_jct
                                           - this._vjdmFwd) + _iteration.Gmin * vbd_jct;
                        }
                        break;
                    default: break;
                }
            }

            /* trap-assisted tunneling and recombination current for reverse bias  */
            Nvtmrssws = ModelTemperature.Vtm0 * ModelTemperature.Njtsswstemp;
            Nvtmrsswgs = ModelTemperature.Vtm0 * ModelTemperature.Njtsswgstemp;
            Nvtmrss = ModelTemperature.Vtm0 * ModelTemperature.Njtsstemp;
            Nvtmrsswd = ModelTemperature.Vtm0 * ModelTemperature.Njtsswdtemp;
            Nvtmrsswgd = ModelTemperature.Vtm0 * ModelTemperature.Njtsswgdtemp;
            Nvtmrsd = ModelTemperature.Vtm0 * ModelTemperature.Njtsdtemp;

            if ((ModelParameters.Vtss - vbs_jct) < (ModelParameters.Vtss * 1e-3))
            {
                T9 = 1.0e3;
                T0 = -vbs_jct / Nvtmrss * T9;
                DEXP(T0, out T1, out T10);
                dT1_dVb = T10 / Nvtmrss * T9;
            }
            else
            {
                T9 = 1.0 / (ModelParameters.Vtss - vbs_jct);
                T0 = -vbs_jct / Nvtmrss * ModelParameters.Vtss * T9;
                dT0_dVb = ModelParameters.Vtss / Nvtmrss * (T9 + vbs_jct * T9 * T9);
                DEXP(T0, out T1, out T10);
                dT1_dVb = T10 * dT0_dVb;
            }

            if ((ModelParameters.Vtsd - vbd_jct) < (ModelParameters.Vtsd * 1e-3))
            {
                T9 = 1.0e3;
                T0 = -vbd_jct / Nvtmrsd * T9;
                DEXP(T0, out T2, out T10);
                dT2_dVb = T10 / Nvtmrsd * T9;
            }
            else
            {
                T9 = 1.0 / (ModelParameters.Vtsd - vbd_jct);
                T0 = -vbd_jct / Nvtmrsd * ModelParameters.Vtsd * T9;
                dT0_dVb = ModelParameters.Vtsd / Nvtmrsd * (T9 + vbd_jct * T9 * T9);
                DEXP(T0, out T2, out T10);
                dT2_dVb = T10 * dT0_dVb;
            }

            if ((ModelParameters.Vtssws - vbs_jct) < (ModelParameters.Vtssws * 1e-3))
            {
                T9 = 1.0e3;
                T0 = -vbs_jct / Nvtmrssws * T9;
                DEXP(T0, out T3, out T10);
                dT3_dVb = T10 / Nvtmrssws * T9;
            }
            else
            {
                T9 = 1.0 / (ModelParameters.Vtssws - vbs_jct);
                T0 = -vbs_jct / Nvtmrssws * ModelParameters.Vtssws * T9;
                dT0_dVb = ModelParameters.Vtssws / Nvtmrssws * (T9 + vbs_jct * T9 * T9);
                DEXP(T0, out T3, out T10);
                dT3_dVb = T10 * dT0_dVb;
            }

            if ((ModelParameters.Vtsswd - vbd_jct) < (ModelParameters.Vtsswd * 1e-3))
            {
                T9 = 1.0e3;
                T0 = -vbd_jct / Nvtmrsswd * T9;
                DEXP(T0, out T4, out T10);
                dT4_dVb = T10 / Nvtmrsswd * T9;
            }
            else
            {
                T9 = 1.0 / (ModelParameters.Vtsswd - vbd_jct);
                T0 = -vbd_jct / Nvtmrsswd * ModelParameters.Vtsswd * T9;
                dT0_dVb = ModelParameters.Vtsswd / Nvtmrsswd * (T9 + vbd_jct * T9 * T9);
                DEXP(T0, out T4, out T10);
                dT4_dVb = T10 * dT0_dVb;
            }

            if ((ModelParameters.Vtsswgs - vbs_jct) < (ModelParameters.Vtsswgs * 1e-3))
            {
                T9 = 1.0e3;
                T0 = -vbs_jct / Nvtmrsswgs * T9;
                DEXP(T0, out T5, out T10);
                dT5_dVb = T10 / Nvtmrsswgs * T9;
            }
            else
            {
                T9 = 1.0 / (ModelParameters.Vtsswgs - vbs_jct);
                T0 = -vbs_jct / Nvtmrsswgs * ModelParameters.Vtsswgs * T9;
                dT0_dVb = ModelParameters.Vtsswgs / Nvtmrsswgs * (T9 + vbs_jct * T9 * T9);
                DEXP(T0, out T5, out T10);
                dT5_dVb = T10 * dT0_dVb;
            }

            if ((ModelParameters.Vtsswgd - vbd_jct) < (ModelParameters.Vtsswgd * 1e-3))
            {
                T9 = 1.0e3;
                T0 = -vbd_jct / Nvtmrsswgd * T9;
                DEXP(T0, out T6, out T10);
                dT6_dVb = T10 / Nvtmrsswgd * T9;
            }
            else
            {
                T9 = 1.0 / (ModelParameters.Vtsswgd - vbd_jct);
                T0 = -vbd_jct / Nvtmrsswgd * ModelParameters.Vtsswgd * T9;
                dT0_dVb = ModelParameters.Vtsswgd / Nvtmrsswgd * (T9 + vbd_jct * T9 * T9);
                DEXP(T0, out T6, out T10);
                dT6_dVb = T10 * dT0_dVb;
            }

            this._gbs += this._sjctTempRevSatCur * dT1_dVb
                                    + this._sswTempRevSatCur * dT3_dVb
                                    + this._sswgTempRevSatCur * dT5_dVb;
            this._cbs -= this._sjctTempRevSatCur * (T1 - 1.0)
                                    + this._sswTempRevSatCur * (T3 - 1.0)
                                    + this._sswgTempRevSatCur * (T5 - 1.0);
            this._gbd += this._djctTempRevSatCur * dT2_dVb
                                    + this._dswTempRevSatCur * dT4_dVb
                                    + this._dswgTempRevSatCur * dT6_dVb;
            this._cbd -= this._djctTempRevSatCur * (T2 - 1.0)
                                    + this._dswTempRevSatCur * (T4 - 1.0)
                                    + this._dswgTempRevSatCur * (T6 - 1.0);

            /* End of diode DC model */

            if (vds >= 0.0)
            {
                this._mode = 1;
                Vds = vds;
                Vgs = vgs;
                Vbs = vbs;
                Vdb = vds - vbs;  /* WDLiu: for GIDL */

            }
            else
            {
                this._mode = -1;
                Vds = -vds;
                Vgs = vgd;
                Vbs = vbd;
                Vdb = -vbs;
            }


            /* dunga */
            if (ModelParameters.MtrlMod.Value != 0)
            {
                epsrox = 3.9;
                toxe = ModelParameters.Eot;
                epssub = EPS0 * ModelParameters.Epsrsub;
            }
            else
            {
                epsrox = ModelParameters.Epsrox;
                toxe = ModelParameters.Toxe;
                epssub = EPSSI;
            }


            T0 = Vbs - this._vbsc - 0.001;
            T1 = Math.Sqrt(T0 * T0 - 0.004 * this._vbsc);
            if (T0 >= 0.0)
            {
                Vbseff = this._vbsc + 0.5 * (T0 + T1);
                dVbseff_dVb = 0.5 * (1.0 + T0 / T1);
            }
            else
            {
                T2 = -0.002 / (T1 - T0);
                Vbseff = this._vbsc * (1.0 + T2);
                dVbseff_dVb = T2 * this._vbsc / T1;
            }

            /* JX: Correction to forward body bias  */
            T9 = 0.95 * Param.BSIM4phi;
            T0 = T9 - Vbseff - 0.001;
            T1 = Math.Sqrt(T0 * T0 + 0.004 * T9);
            Vbseff = T9 - 0.5 * (T0 + T1);
            dVbseff_dVb *= 0.5 * (1.0 + T0 / T1);
            Phis = Param.BSIM4phi - Vbseff;
            dPhis_dVb = -1.0;
            sqrtPhis = Math.Sqrt(Phis);
            dsqrtPhis_dVb = -0.5 / sqrtPhis;

            Xdep = Param.BSIM4Xdep0 * sqrtPhis / Param.BSIM4sqrtPhi;
            dXdep_dVb = (Param.BSIM4Xdep0 / Param.BSIM4sqrtPhi)
                      * dsqrtPhis_dVb;

            Leff = Param.BSIM4leff;
            Vtm = ModelTemperature.Vtm;
            Vtm0 = ModelTemperature.Vtm0;

            /* Vth Calculation */
            T3 = Math.Sqrt(Xdep);
            V0 = Param.BSIM4vbi - Param.BSIM4phi;

            T0 = Param.BSIM4dvt2 * Vbseff;
            if (T0 >= -0.5)
            {
                T1 = 1.0 + T0;
                T2 = Param.BSIM4dvt2;
            }
            else
            {
                T4 = 1.0 / (3.0 + 8.0 * T0);
                T1 = (1.0 + 3.0 * T0) * T4;
                T2 = Param.BSIM4dvt2 * T4 * T4;
            }
            lt1 = ModelTemperature.Factor1 * T3 * T1;
            dlt1_dVb = ModelTemperature.Factor1 * (0.5 / T3 * T1 * dXdep_dVb + T3 * T2);

            T0 = Param.BSIM4dvt2w * Vbseff;
            if (T0 >= -0.5)
            {
                T1 = 1.0 + T0;
                T2 = Param.BSIM4dvt2w;
            }
            else
            {
                T4 = 1.0 / (3.0 + 8.0 * T0);
                T1 = (1.0 + 3.0 * T0) * T4;
                T2 = Param.BSIM4dvt2w * T4 * T4;
            }
            ltw = ModelTemperature.Factor1 * T3 * T1;
            dltw_dVb = ModelTemperature.Factor1 * (0.5 / T3 * T1 * dXdep_dVb + T3 * T2);

            T0 = Param.BSIM4dvt1 * Leff / lt1;
            if (T0 < EXP_THRESHOLD)
            {
                T1 = Math.Exp(T0);
                T2 = T1 - 1.0;
                T3 = T2 * T2;
                T4 = T3 + 2.0 * T1 * MIN_EXP;
                Theta0 = T1 / T4;
                dT1_dVb = -T0 * T1 * dlt1_dVb / lt1;
                dTheta0_dVb = dT1_dVb * (T4 - 2.0 * T1 * (T2 + MIN_EXP)) / T4 / T4;
            }
            else
            {
                Theta0 = 1.0 / (MAX_EXP - 2.0); /* 3.0 * MIN_EXP omitted */
                dTheta0_dVb = 0.0;
            }
            this._thetavth = Param.BSIM4dvt0 * Theta0;
            Delt_vth = this._thetavth * V0;
            dDelt_vth_dVb = Param.BSIM4dvt0 * dTheta0_dVb * V0;

            T0 = Param.BSIM4dvt1w * Param.BSIM4weff * Leff / ltw;
            if (T0 < EXP_THRESHOLD)
            {
                T1 = Math.Exp(T0);
                T2 = T1 - 1.0;
                T3 = T2 * T2;
                T4 = T3 + 2.0 * T1 * MIN_EXP;
                T5 = T1 / T4;
                dT1_dVb = -T0 * T1 * dltw_dVb / ltw;
                dT5_dVb = dT1_dVb * (T4 - 2.0 * T1 * (T2 + MIN_EXP)) / T4 / T4;
            }
            else
            {
                T5 = 1.0 / (MAX_EXP - 2.0); /* 3.0 * MIN_EXP omitted */
                dT5_dVb = 0.0;
            }
            T0 = Param.BSIM4dvt0w * T5;
            T2 = T0 * V0;
            dT2_dVb = Param.BSIM4dvt0w * dT5_dVb * V0;

            TempRatio = _temperature.Temperature / ModelParameters.Tnom - 1.0;
            T0 = Math.Sqrt(1.0 + Param.BSIM4lpe0 / Leff);
            T1 = Param.BSIM4k1ox * (T0 - 1.0) * Param.BSIM4sqrtPhi
               + (Param.BSIM4kt1 + Param.BSIM4kt1l / Leff
               + Param.BSIM4kt2 * Vbseff) * TempRatio;
            Vth_NarrowW = toxe * Param.BSIM4phi
                        / (Param.BSIM4weff + Param.BSIM4w0);

            T3 = this._eta0 + Param.BSIM4etab * Vbseff;
            if (T3 < 1.0e-4)
            {
                T9 = 1.0 / (3.0 - 2.0e4 * T3);
                T3 = (2.0e-4 - T3) * T9;
                T4 = T9 * T9;
            }
            else
            {
                T4 = 1.0;
            }
            dDIBL_Sft_dVd = T3 * Param.BSIM4theta0vb0;
            DIBL_Sft = dDIBL_Sft_dVd * Vds;

            Lpe_Vb = Math.Sqrt(1.0 + Param.BSIM4lpeb / Leff);

            Vth = ModelParameters.Type * this._vth0 + (Param.BSIM4k1ox * sqrtPhis
                - Param.BSIM4k1 * Param.BSIM4sqrtPhi) * Lpe_Vb
                - this._k2ox * Vbseff - Delt_vth - T2 + (Param.BSIM4k3
                + Param.BSIM4k3b * Vbseff) * Vth_NarrowW + T1 - DIBL_Sft;

            dVth_dVb = Lpe_Vb * Param.BSIM4k1ox * dsqrtPhis_dVb - this._k2ox
                     - dDelt_vth_dVb - dT2_dVb + Param.BSIM4k3b * Vth_NarrowW
                     - Param.BSIM4etab * Vds * Param.BSIM4theta0vb0 * T4
                     + Param.BSIM4kt2 * TempRatio;
            dVth_dVd = -dDIBL_Sft_dVd;


            /* Calculate n */
            tmp1 = epssub / Xdep;
            this._nstar = ModelTemperature.Vtm / Constants.Charge * (ModelTemperature.Coxe
                             + tmp1 + Param.BSIM4cit);
            tmp2 = Param.BSIM4nfactor * tmp1;
            tmp3 = Param.BSIM4cdsc + Param.BSIM4cdscb * Vbseff
                 + Param.BSIM4cdscd * Vds;
            tmp4 = (tmp2 + tmp3 * Theta0 + Param.BSIM4cit) / ModelTemperature.Coxe;
            if (tmp4 >= -0.5)
            {
                n = 1.0 + tmp4;
                dn_dVb = (-tmp2 / Xdep * dXdep_dVb + tmp3 * dTheta0_dVb
                       + Param.BSIM4cdscb * Theta0) / ModelTemperature.Coxe;
                dn_dVd = Param.BSIM4cdscd * Theta0 / ModelTemperature.Coxe;
            }
            else
            {
                T0 = 1.0 / (3.0 + 8.0 * tmp4);
                n = (1.0 + 3.0 * tmp4) * T0;
                T0 *= T0;
                dn_dVb = (-tmp2 / Xdep * dXdep_dVb + tmp3 * dTheta0_dVb
                       + Param.BSIM4cdscb * Theta0) / ModelTemperature.Coxe * T0;
                dn_dVd = Param.BSIM4cdscd * Theta0 / ModelTemperature.Coxe * T0;
            }


            /* Vth correction for Pocket implant */
            if (Param.BSIM4dvtp0 > 0.0)
            {
                T0 = -Param.BSIM4dvtp1 * Vds;
                if (T0 < -EXP_THRESHOLD)
                {
                    T2 = MIN_EXP;
                    dT2_dVd = 0.0;
                }
                else
                {
                    T2 = Math.Exp(T0);
                    dT2_dVd = -Param.BSIM4dvtp1 * T2;
                }

                T3 = Leff + Param.BSIM4dvtp0 * (1.0 + T2);
                dT3_dVd = Param.BSIM4dvtp0 * dT2_dVd;
                if (ModelParameters.TempMod < 2)
                {
                    T4 = Vtm * Math.Log(Leff / T3);
                    dT4_dVd = -Vtm * dT3_dVd / T3;
                }
                else
                {
                    T4 = ModelTemperature.Vtm0 * Math.Log(Leff / T3);
                    dT4_dVd = -ModelTemperature.Vtm0 * dT3_dVd / T3;
                }
                dDITS_Sft_dVd = dn_dVd * T4 + n * dT4_dVd;
                dDITS_Sft_dVb = T4 * dn_dVb;

                Vth -= n * T4;
                dVth_dVd -= dDITS_Sft_dVd;
                dVth_dVb -= dDITS_Sft_dVb;
            }

            /* v4.7 DITS_SFT2  */
            if ((Param.BSIM4dvtp4 == 0.0) || (Param.BSIM4dvtp2factor == 0.0))
            {
                T0 = 0.0;
                DITS_Sft2 = 0.0;
            }
            else
            {
                //T0 = Math.Exp(2.0 * Param.BSIM4dvtp4 * Vds);   /* beta code */
                T1 = 2.0 * Param.BSIM4dvtp4 * Vds;
                DEXP(T1, out T0, out T10);
                DITS_Sft2 = Param.BSIM4dvtp2factor * (T0 - 1) / (T0 + 1);
                //dDITS_Sft2_dVd = Param.BSIM4dvtp2factor * Param.BSIM4dvtp4 * 4.0 * T0 / ((T0+1) * (T0+1));   /* beta code */
                dDITS_Sft2_dVd = Param.BSIM4dvtp2factor * Param.BSIM4dvtp4 * 4.0 * T10 / ((T0 + 1) * (T0 + 1));
                Vth -= DITS_Sft2;
                dVth_dVd -= dDITS_Sft2_dVd;
            }



            this._von = Vth;


            /* Poly Gate Si Depletion Effect */
            T0 = this._vfb + Param.BSIM4phi;
            if (ModelParameters.MtrlMod == 0)
                T1 = EPSSI;
            else
                T1 = ModelParameters.Epsrgate * EPS0;


            BSIM4polyDepletion(T0, Param.BSIM4ngate, T1, ModelTemperature.Coxe, vgs, out vgs_eff, out dvgs_eff_dvg);

            BSIM4polyDepletion(T0, Param.BSIM4ngate, T1, ModelTemperature.Coxe, vgd, out vgd_eff, out dvgd_eff_dvg);

            if (this._mode > 0)
            {
                Vgs_eff = vgs_eff;
                dVgs_eff_dVg = dvgs_eff_dvg;
            }
            else
            {
                Vgs_eff = vgd_eff;
                dVgs_eff_dVg = dvgd_eff_dvg;
            }
            this._vgs_eff = vgs_eff;
            this._vgd_eff = vgd_eff;
            this._dvgs_eff_dvg = dvgs_eff_dvg;
            this._dvgd_eff_dvg = dvgd_eff_dvg;


            Vgst = Vgs_eff - Vth;

            /* Calculate Vgsteff */
            T0 = n * Vtm;
            T1 = Param.BSIM4mstar * Vgst;
            T2 = T1 / T0;
            if (T2 > EXP_THRESHOLD)
            {
                T10 = T1;
                dT10_dVg = Param.BSIM4mstar * dVgs_eff_dVg;
                dT10_dVd = -dVth_dVd * Param.BSIM4mstar;
                dT10_dVb = -dVth_dVb * Param.BSIM4mstar;
            }
            else if (T2 < -EXP_THRESHOLD)
            {
                T10 = Vtm * Math.Log(1.0 + MIN_EXP);
                dT10_dVg = 0.0;
                dT10_dVd = T10 * dn_dVd;
                dT10_dVb = T10 * dn_dVb;
                T10 *= n;
            }
            else
            {
                ExpVgst = Math.Exp(T2);
                T3 = Vtm * Math.Log(1.0 + ExpVgst);
                T10 = n * T3;
                dT10_dVg = Param.BSIM4mstar * ExpVgst / (1.0 + ExpVgst);
                dT10_dVb = T3 * dn_dVb - dT10_dVg * (dVth_dVb + Vgst * dn_dVb / n);
                dT10_dVd = T3 * dn_dVd - dT10_dVg * (dVth_dVd + Vgst * dn_dVd / n);
                dT10_dVg *= dVgs_eff_dVg;
            }

            T1 = Param.BSIM4voffcbn - (1.0 - Param.BSIM4mstar) * Vgst;
            T2 = T1 / T0;
            if (T2 < -EXP_THRESHOLD)
            {
                T3 = ModelTemperature.Coxe * MIN_EXP / Param.BSIM4cdep0;
                T9 = Param.BSIM4mstar + T3 * n;
                dT9_dVg = 0.0;
                dT9_dVd = dn_dVd * T3;
                dT9_dVb = dn_dVb * T3;
            }
            else if (T2 > EXP_THRESHOLD)
            {
                T3 = ModelTemperature.Coxe * MAX_EXP / Param.BSIM4cdep0;
                T9 = Param.BSIM4mstar + T3 * n;
                dT9_dVg = 0.0;
                dT9_dVd = dn_dVd * T3;
                dT9_dVb = dn_dVb * T3;
            }
            else
            {
                ExpVgst = Math.Exp(T2);
                T3 = ModelTemperature.Coxe / Param.BSIM4cdep0;
                T4 = T3 * ExpVgst;
                T5 = T1 * T4 / T0;
                T9 = Param.BSIM4mstar + n * T4;
                dT9_dVg = T3 * (Param.BSIM4mstar - 1.0) * ExpVgst / Vtm;
                dT9_dVb = T4 * dn_dVb - dT9_dVg * dVth_dVb - T5 * dn_dVb;
                dT9_dVd = T4 * dn_dVd - dT9_dVg * dVth_dVd - T5 * dn_dVd;
                dT9_dVg *= dVgs_eff_dVg;
            }
            this._vgsteff = Vgsteff = T10 / T9;
            T11 = T9 * T9;
            dVgsteff_dVg = (T9 * dT10_dVg - T10 * dT9_dVg) / T11;
            dVgsteff_dVd = (T9 * dT10_dVd - T10 * dT9_dVd) / T11;
            dVgsteff_dVb = (T9 * dT10_dVb - T10 * dT9_dVb) / T11;

            /* Calculate Effective Channel Geometry */
            T9 = sqrtPhis - Param.BSIM4sqrtPhi;
            Weff = Param.BSIM4weff - 2.0 * (Param.BSIM4dwg * Vgsteff
                 + Param.BSIM4dwb * T9);
            dWeff_dVg = -2.0 * Param.BSIM4dwg;
            dWeff_dVb = -2.0 * Param.BSIM4dwb * dsqrtPhis_dVb;

            if (Weff < 2.0e-8) /* to avoid the discontinuity problem due to Weff*/
            {
                T0 = 1.0 / (6.0e-8 - 2.0 * Weff);
                Weff = 2.0e-8 * (4.0e-8 - Weff) * T0;
                T0 *= T0 * 4.0e-16;
                dWeff_dVg *= T0;
                dWeff_dVb *= T0;
            }

            if (ModelParameters.RdsMod == 1)
                Rds = dRds_dVg = dRds_dVb = 0.0;
            else
            {
                T0 = 1.0 + Param.BSIM4prwg * Vgsteff;
                dT0_dVg = -Param.BSIM4prwg / T0 / T0;
                T1 = Param.BSIM4prwb * T9;
                dT1_dVb = Param.BSIM4prwb * dsqrtPhis_dVb;

                T2 = 1.0 / T0 + T1;
                T3 = T2 + Math.Sqrt(T2 * T2 + 0.01); /* 0.01 = 4.0 * 0.05 * 0.05 */
                dT3_dVg = 1.0 + T2 / (T3 - T2);
                dT3_dVb = dT3_dVg * dT1_dVb;
                dT3_dVg *= dT0_dVg;

                T4 = Param.BSIM4rds0 * 0.5;
                Rds = Param.BSIM4rdswmin + T3 * T4;
                dRds_dVg = T4 * dT3_dVg;
                dRds_dVb = T4 * dT3_dVb;

                if (Rds > 0.0)
                    this._grdsw = 1.0 / Rds * Parameters.Nf; /*4.6.2*/
                else
                    this._grdsw = 0.0;
            }

            /* Calculate Abulk */
            T9 = 0.5 * Param.BSIM4k1ox * Lpe_Vb / sqrtPhis;
            T1 = T9 + this._k2ox - Param.BSIM4k3b * Vth_NarrowW;
            dT1_dVb = -T9 / sqrtPhis * dsqrtPhis_dVb;

            T9 = Math.Sqrt(Param.BSIM4xj * Xdep);
            tmp1 = Leff + 2.0 * T9;
            T5 = Leff / tmp1;
            tmp2 = Param.BSIM4a0 * T5;
            tmp3 = Param.BSIM4weff + Param.BSIM4b1;
            tmp4 = Param.BSIM4b0 / tmp3;
            T2 = tmp2 + tmp4;
            dT2_dVb = -T9 / tmp1 / Xdep * dXdep_dVb;
            T6 = T5 * T5;
            T7 = T5 * T6;

            Abulk0 = 1.0 + T1 * T2;
            dAbulk0_dVb = T1 * tmp2 * dT2_dVb + T2 * dT1_dVb;

            T8 = Param.BSIM4ags * Param.BSIM4a0 * T7;
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
            {
                T9 = 1.0 / (3.0 - 20.0 * Abulk);
                Abulk = (0.2 - Abulk) * T9;
                T10 = T9 * T9;
                dAbulk_dVb *= T10;
                dAbulk_dVg *= T10;
            }
            this._abulk = Abulk;

            T2 = Param.BSIM4keta * Vbseff;
            if (T2 >= -0.9)
            {
                T0 = 1.0 / (1.0 + T2);
                dT0_dVb = -Param.BSIM4keta * T0 * T0;
            }
            else
            {
                T1 = 1.0 / (0.8 + T2);
                T0 = (17.0 + 20.0 * T2) * T1;
                dT0_dVb = -Param.BSIM4keta * T1 * T1;
            }
            dAbulk_dVg *= T0;
            dAbulk_dVb = dAbulk_dVb * T0 + Abulk * dT0_dVb;
            dAbulk0_dVb = dAbulk0_dVb * T0 + Abulk0 * dT0_dVb;
            Abulk *= T0;
            Abulk0 *= T0;

            /* Mobility calculation */
            if (ModelParameters.MtrlMod.Value != 0 && ModelParameters.MtrlCompatMod == 0)
                T14 = 2.0 * ModelParameters.Type * (ModelParameters.Phig - ModelParameters.Easub - 0.5 * ModelTemperature.Eg0 + 0.45);
            else
                T14 = 0.0;

            if (ModelParameters.MobMod == 0)
            {
                T0 = Vgsteff + Vth + Vth - T14;
                T2 = Param.BSIM4ua + Param.BSIM4uc * Vbseff;
                T3 = T0 / toxe;
                T12 = Math.Sqrt(Vth * Vth + 0.0001);
                T9 = 1.0 / (Vgsteff + 2 * T12);
                T10 = T9 * toxe;
                T8 = Param.BSIM4ud * T10 * T10 * Vth;
                T6 = T8 * Vth;
                T5 = T3 * (T2 + Param.BSIM4ub * T3) + T6;
                T7 = -2.0 * T6 * T9;
                T11 = T7 * Vth / T12;
                dDenomi_dVg = (T2 + 2.0 * Param.BSIM4ub * T3) / toxe;
                T13 = 2.0 * (dDenomi_dVg + T11 + T8);
                dDenomi_dVd = T13 * dVth_dVd;
                dDenomi_dVb = T13 * dVth_dVb + Param.BSIM4uc * T3;
                dDenomi_dVg += T7;
            }
            else if (ModelParameters.MobMod == 1)
            {
                T0 = Vgsteff + Vth + Vth - T14;
                T2 = 1.0 + Param.BSIM4uc * Vbseff;
                T3 = T0 / toxe;
                T4 = T3 * (Param.BSIM4ua + Param.BSIM4ub * T3);
                T12 = Math.Sqrt(Vth * Vth + 0.0001);
                T9 = 1.0 / (Vgsteff + 2 * T12);
                T10 = T9 * toxe;
                T8 = Param.BSIM4ud * T10 * T10 * Vth;
                T6 = T8 * Vth;
                T5 = T4 * T2 + T6;
                T7 = -2.0 * T6 * T9;
                T11 = T7 * Vth / T12;
                dDenomi_dVg = (Param.BSIM4ua + 2.0 * Param.BSIM4ub * T3) * T2 / toxe;
                T13 = 2.0 * (dDenomi_dVg + T11 + T8);
                dDenomi_dVd = T13 * dVth_dVd;
                dDenomi_dVb = T13 * dVth_dVb + Param.BSIM4uc * T4;
                dDenomi_dVg += T7;
            }
            else if (ModelParameters.MobMod == 2)
            {
                T0 = (Vgsteff + this._vtfbphi1) / toxe;
                T1 = Math.Exp(Param.BSIM4eu * Math.Log(T0));
                dT1_dVg = T1 * Param.BSIM4eu / T0 / toxe;
                T2 = Param.BSIM4ua + Param.BSIM4uc * Vbseff;

                T12 = Math.Sqrt(Vth * Vth + 0.0001);
                T9 = 1.0 / (Vgsteff + 2 * T12);
                T10 = T9 * toxe;
                T8 = Param.BSIM4ud * T10 * T10 * Vth;
                T6 = T8 * Vth;
                T5 = T1 * T2 + T6;
                T7 = -2.0 * T6 * T9;
                T11 = T7 * Vth / T12;
                dDenomi_dVg = T2 * dT1_dVg + T7;
                T13 = 2.0 * (T11 + T8);
                dDenomi_dVd = T13 * dVth_dVd;
                dDenomi_dVb = T13 * dVth_dVb + T1 * Param.BSIM4uc;
            }
            else if (ModelParameters.MobMod == 4) /* Synopsys 08/30/2013 add */
            {
                T0 = Vgsteff + this._vtfbphi1 - T14;
                T2 = Param.BSIM4ua + Param.BSIM4uc * Vbseff;
                T3 = T0 / toxe;
                T12 = Math.Sqrt(this._vtfbphi1 * this._vtfbphi1 + 0.0001);
                T9 = 1.0 / (Vgsteff + 2 * T12);
                T10 = T9 * toxe;
                T8 = Param.BSIM4ud * T10 * T10 * this._vtfbphi1;
                T6 = T8 * this._vtfbphi1;
                T5 = T3 * (T2 + Param.BSIM4ub * T3) + T6;
                T7 = -2.0 * T6 * T9;
                dDenomi_dVg = (T2 + 2.0 * Param.BSIM4ub * T3) / toxe;
                dDenomi_dVd = 0.0;
                dDenomi_dVb = Param.BSIM4uc * T3;
                dDenomi_dVg += T7;
            }
            else if (ModelParameters.MobMod == 5) /* Synopsys 08/30/2013 add */
            {
                T0 = Vgsteff + this._vtfbphi1 - T14;
                T2 = 1.0 + Param.BSIM4uc * Vbseff;
                T3 = T0 / toxe;
                T4 = T3 * (Param.BSIM4ua + Param.BSIM4ub * T3);
                T12 = Math.Sqrt(this._vtfbphi1 * this._vtfbphi1 + 0.0001);
                T9 = 1.0 / (Vgsteff + 2 * T12);
                T10 = T9 * toxe;
                T8 = Param.BSIM4ud * T10 * T10 * this._vtfbphi1;
                T6 = T8 * this._vtfbphi1;
                T5 = T4 * T2 + T6;
                T7 = -2.0 * T6 * T9;
                dDenomi_dVg = (Param.BSIM4ua + 2.0 * Param.BSIM4ub * T3) * T2
                            / toxe;
                dDenomi_dVd = 0.0;
                dDenomi_dVb = Param.BSIM4uc * T4;
                dDenomi_dVg += T7;
            }
            else if (ModelParameters.MobMod == 6) /* Synopsys 08/30/2013 modify */
            {
                T0 = (Vgsteff + this._vtfbphi1) / toxe;
                T1 = Math.Exp(Param.BSIM4eu * Math.Log(T0));
                dT1_dVg = T1 * Param.BSIM4eu / T0 / toxe;
                T2 = Param.BSIM4ua + Param.BSIM4uc * Vbseff;

                T12 = Math.Sqrt(this._vtfbphi1 * this._vtfbphi1 + 0.0001);
                T9 = 1.0 / (Vgsteff + 2 * T12);
                T10 = T9 * toxe;
                T8 = Param.BSIM4ud * T10 * T10 * this._vtfbphi1;
                T6 = T8 * this._vtfbphi1;
                T5 = T1 * T2 + T6;
                T7 = -2.0 * T6 * T9;
                dDenomi_dVg = T2 * dT1_dVg + T7;
                dDenomi_dVd = 0;
                dDenomi_dVb = T1 * Param.BSIM4uc;
            }

            /*high K mobility*/
            else
            {


                /*univsersal mobility*/
                T0 = (Vgsteff + this._vtfbphi1) * 1.0e-8 / toxe / 6.0;
                T1 = Math.Exp(Param.BSIM4eu * Math.Log(T0));
                dT1_dVg = T1 * Param.BSIM4eu * 1.0e-8 / T0 / toxe / 6.0;
                T2 = Param.BSIM4ua + Param.BSIM4uc * Vbseff;

                /*Coulombic*/
                VgsteffVth = Param.BSIM4VgsteffVth;

                T10 = Math.Exp(Param.BSIM4ucs * Math.Log(0.5 + 0.5 * Vgsteff / VgsteffVth));
                T11 = Param.BSIM4ud / T10;
                dT11_dVg = -0.5 * Param.BSIM4ucs * T11 / (0.5 + 0.5 * Vgsteff / VgsteffVth) / VgsteffVth;

                dDenomi_dVg = T2 * dT1_dVg + dT11_dVg;
                dDenomi_dVd = 0.0;
                dDenomi_dVb = T1 * Param.BSIM4uc;

                T5 = T1 * T2 + T11;
            }





            if (T5 >= -0.8)
            {
                Denomi = 1.0 + T5;
            }
            else
            {
                T9 = 1.0 / (7.0 + 10.0 * T5);
                Denomi = (0.6 + T5) * T9;
                T9 *= T9;
                dDenomi_dVg *= T9;
                dDenomi_dVd *= T9;
                dDenomi_dVb *= T9;
            }


            this._ueff = ueff = this._u0temp / Denomi;
            T9 = -ueff / Denomi;
            dueff_dVg = T9 * dDenomi_dVg;
            dueff_dVd = T9 * dDenomi_dVd;
            dueff_dVb = T9 * dDenomi_dVb;

            /* Saturation Drain Voltage  Vdsat */
            WVCox = Weff * this._vsattemp * ModelTemperature.Coxe;
            WVCoxRds = WVCox * Rds;

            Esat = 2.0 * this._vsattemp / ueff;
            this._esatL = EsatL = Esat * Leff;
            T0 = -EsatL / ueff;
            dEsatL_dVg = T0 * dueff_dVg;
            dEsatL_dVd = T0 * dueff_dVd;
            dEsatL_dVb = T0 * dueff_dVb;

            /* Sqrt() */
            a1 = Param.BSIM4a1;
            if (a1 == 0.0)
            {
                Lambda = Param.BSIM4a2;
                dLambda_dVg = 0.0;
            }
            else if (a1 > 0.0)
            {
                T0 = 1.0 - Param.BSIM4a2;
                T1 = T0 - Param.BSIM4a1 * Vgsteff - 0.0001;
                T2 = Math.Sqrt(T1 * T1 + 0.0004 * T0);
                Lambda = Param.BSIM4a2 + T0 - 0.5 * (T1 + T2);
                dLambda_dVg = 0.5 * Param.BSIM4a1 * (1.0 + T1 / T2);
            }
            else
            {
                T1 = Param.BSIM4a2 + Param.BSIM4a1 * Vgsteff - 0.0001;
                T2 = Math.Sqrt(T1 * T1 + 0.0004 * Param.BSIM4a2);
                Lambda = 0.5 * (T1 + T2);
                dLambda_dVg = 0.5 * Param.BSIM4a1 * (1.0 + T1 / T2);
            }

            Vgst2Vtm = Vgsteff + 2.0 * Vtm;
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

            /* Calculate Vdseff */
            T1 = Vdsat - Vds - Param.BSIM4delta;
            dT1_dVg = dVdsat_dVg;
            dT1_dVd = dVdsat_dVd - 1.0;
            dT1_dVb = dVdsat_dVb;

            T2 = Math.Sqrt(T1 * T1 + 4.0 * Param.BSIM4delta * Vdsat);
            T0 = T1 / T2;
            T9 = 2.0 * Param.BSIM4delta;
            T3 = T9 / T2;
            dT2_dVg = T0 * dT1_dVg + T3 * dVdsat_dVg;
            dT2_dVd = T0 * dT1_dVd + T3 * dVdsat_dVd;
            dT2_dVb = T0 * dT1_dVb + T3 * dVdsat_dVb;

            if (T1 >= 0.0)
            {
                Vdseff = Vdsat - 0.5 * (T1 + T2);
                dVdseff_dVg = dVdsat_dVg - 0.5 * (dT1_dVg + dT2_dVg);
                dVdseff_dVd = dVdsat_dVd - 0.5 * (dT1_dVd + dT2_dVd);
                dVdseff_dVb = dVdsat_dVb - 0.5 * (dT1_dVb + dT2_dVb);
            }
            else
            {
                T4 = T9 / (T2 - T1);
                T5 = 1.0 - T4;
                T6 = Vdsat * T4 / (T2 - T1);
                Vdseff = Vdsat * T5;
                dVdseff_dVg = dVdsat_dVg * T5 + T6 * (dT2_dVg - dT1_dVg);
                dVdseff_dVd = dVdsat_dVd * T5 + T6 * (dT2_dVd - dT1_dVd);
                dVdseff_dVb = dVdsat_dVb * T5 + T6 * (dT2_dVb - dT1_dVb);
            }

            if (Vds == 0.0)
            {
                Vdseff = 0.0;
                dVdseff_dVg = 0.0;
                dVdseff_dVb = 0.0;
            }

            if (Vdseff > Vds)
                Vdseff = Vds;
            diffVds = Vds - Vdseff;
            this._vdseff = Vdseff;

            /* Velocity Overshoot */
            if ((ModelParameters.Lambda.Given) && (ModelParameters.Lambda > 0.0))
            {
                T1 = Leff * ueff;
                T2 = Param.BSIM4lambda / T1;
                T3 = -T2 / T1 * Leff;
                dT2_dVd = T3 * dueff_dVd;
                dT2_dVg = T3 * dueff_dVg;
                dT2_dVb = T3 * dueff_dVb;
                T5 = 1.0 / (Esat * Param.BSIM4litl);
                T4 = -T5 / EsatL;
                dT5_dVg = dEsatL_dVg * T4;
                dT5_dVd = dEsatL_dVd * T4;
                dT5_dVb = dEsatL_dVb * T4;
                T6 = 1.0 + diffVds * T5;
                dT6_dVg = dT5_dVg * diffVds - dVdseff_dVg * T5;
                dT6_dVd = dT5_dVd * diffVds + (1.0 - dVdseff_dVd) * T5;
                dT6_dVb = dT5_dVb * diffVds - dVdseff_dVb * T5;
                T7 = 2.0 / (T6 * T6 + 1.0);
                T8 = 1.0 - T7;
                T9 = T6 * T7 * T7;
                dT8_dVg = T9 * dT6_dVg;
                dT8_dVd = T9 * dT6_dVd;
                dT8_dVb = T9 * dT6_dVb;
                T10 = 1.0 + T2 * T8;
                dT10_dVg = dT2_dVg * T8 + T2 * dT8_dVg;
                dT10_dVd = dT2_dVd * T8 + T2 * dT8_dVd;
                dT10_dVb = dT2_dVb * T8 + T2 * dT8_dVb;
                if (T10 == 1.0)
                    dT10_dVg = dT10_dVd = dT10_dVb = 0.0;

                dEsatL_dVg *= T10;
                dEsatL_dVg += EsatL * dT10_dVg;
                dEsatL_dVd *= T10;
                dEsatL_dVd += EsatL * dT10_dVd;
                dEsatL_dVb *= T10;
                dEsatL_dVb += EsatL * dT10_dVb;
                EsatL *= T10;
                Esat = EsatL / Leff;  /* bugfix by Wenwei Yang (4.6.4) */
                this._esatL = EsatL;
            }

            /* Calculate Vasat */
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

            /* Calculate Idl first */

            tmp1 = this._vtfbphi2;
            tmp2 = 2.0e8 * this._toxp;
            dT0_dVg = 1.0 / tmp2;
            T0 = (Vgsteff + tmp1) * dT0_dVg;

            tmp3 = Math.Exp(ModelParameters.Bdos * 0.7 * Math.Log(T0));
            T1 = 1.0 + tmp3;
            T2 = ModelParameters.Bdos * 0.7 * tmp3 / T0;
            Tcen = ModelParameters.Ados * 1.9e-9 / T1;
            dTcen_dVg = -Tcen * T2 * dT0_dVg / T1;

            Coxeff = epssub * this._coxp
                   / (epssub + this._coxp * Tcen);
            this._coxeff = Coxeff;
            dCoxeff_dVg = -Coxeff * Coxeff * dTcen_dVg / epssub;

            CoxeffWovL = Coxeff * Weff / Leff;
            beta = ueff * CoxeffWovL;
            T3 = ueff / Leff;
            dbeta_dVg = CoxeffWovL * dueff_dVg + T3
                      * (Weff * dCoxeff_dVg + Coxeff * dWeff_dVg);
            dbeta_dVd = CoxeffWovL * dueff_dVd;
            dbeta_dVb = CoxeffWovL * dueff_dVb + T3 * Coxeff * dWeff_dVb;

            this._abovVgst2Vtm = Abulk / Vgst2Vtm;
            T0 = 1.0 - 0.5 * Vdseff * this._abovVgst2Vtm;
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
            Idl = gche / T0;
            T1 = (1.0 - Idl * Rds) / T0;
            T2 = Idl * Idl;
            dIdl_dVg = T1 * dgche_dVg - T2 * dRds_dVg;
            dIdl_dVd = T1 * dgche_dVd;
            dIdl_dVb = T1 * dgche_dVb - T2 * dRds_dVb;

            /* Calculate degradation factor due to pocket implant */

            if (Param.BSIM4fprout <= 0.0)
            {
                FP = 1.0;
                dFP_dVg = 0.0;
            }
            else
            {
                T9 = Param.BSIM4fprout * Math.Sqrt(Leff) / Vgst2Vtm;
                FP = 1.0 / (1.0 + T9);
                dFP_dVg = FP * FP * T9 / Vgst2Vtm;
            }

            /* Calculate VACLM */
            T8 = Param.BSIM4pvag / EsatL;
            T9 = T8 * Vgsteff;
            if (T9 > -0.9)
            {
                PvagTerm = 1.0 + T9;
                dPvagTerm_dVg = T8 * (1.0 - Vgsteff * dEsatL_dVg / EsatL);
                dPvagTerm_dVb = -T9 * dEsatL_dVb / EsatL;
                dPvagTerm_dVd = -T9 * dEsatL_dVd / EsatL;
            }
            else
            {
                T4 = 1.0 / (17.0 + 20.0 * T9);
                PvagTerm = (0.8 + T9) * T4;
                T4 *= T4;
                dPvagTerm_dVg = T8 * (1.0 - Vgsteff * dEsatL_dVg / EsatL) * T4;
                T9 *= T4 / EsatL;
                dPvagTerm_dVb = -T9 * dEsatL_dVb;
                dPvagTerm_dVd = -T9 * dEsatL_dVd;
            }

            if ((Param.BSIM4pclm > MIN_EXP) && (diffVds > 1.0e-10))
            {
                T0 = 1.0 + Rds * Idl;
                dT0_dVg = dRds_dVg * Idl + Rds * dIdl_dVg;
                dT0_dVd = Rds * dIdl_dVd;
                dT0_dVb = dRds_dVb * Idl + Rds * dIdl_dVb;

                T2 = Vdsat / Esat;
                T1 = Leff + T2;
                dT1_dVg = (dVdsat_dVg - T2 * dEsatL_dVg / Leff) / Esat;
                dT1_dVd = (dVdsat_dVd - T2 * dEsatL_dVd / Leff) / Esat;
                dT1_dVb = (dVdsat_dVb - T2 * dEsatL_dVb / Leff) / Esat;

                Cclm = FP * PvagTerm * T0 * T1 / (Param.BSIM4pclm * Param.BSIM4litl);
                dCclm_dVg = Cclm * (dFP_dVg / FP + dPvagTerm_dVg / PvagTerm
                          + dT0_dVg / T0 + dT1_dVg / T1);
                dCclm_dVb = Cclm * (dPvagTerm_dVb / PvagTerm + dT0_dVb / T0
                          + dT1_dVb / T1);
                dCclm_dVd = Cclm * (dPvagTerm_dVd / PvagTerm + dT0_dVd / T0
                          + dT1_dVd / T1);
                VACLM = Cclm * diffVds;

                dVACLM_dVg = dCclm_dVg * diffVds - dVdseff_dVg * Cclm;
                dVACLM_dVb = dCclm_dVb * diffVds - dVdseff_dVb * Cclm;
                dVACLM_dVd = dCclm_dVd * diffVds + (1.0 - dVdseff_dVd) * Cclm;
            }
            else
            {
                VACLM = Cclm = MAX_EXP;
                dVACLM_dVd = dVACLM_dVg = dVACLM_dVb = 0.0;
                dCclm_dVd = dCclm_dVg = dCclm_dVb = 0.0;
            }

            /* Calculate VADIBL */
            if (Param.BSIM4thetaRout > MIN_EXP)
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
                T2 = Param.BSIM4thetaRout;
                VADIBL = (Vgst2Vtm - T0 / T1) / T2;
                dVADIBL_dVg = (1.0 - dT0_dVg / T1 + T0 * dT1_dVg / T9) / T2;
                dVADIBL_dVb = (-dT0_dVb / T1 + T0 * dT1_dVb / T9) / T2;
                dVADIBL_dVd = (-dT0_dVd / T1 + T0 * dT1_dVd / T9) / T2;

                T7 = Param.BSIM4pdiblb * Vbseff;
                if (T7 >= -0.9)
                {
                    T3 = 1.0 / (1.0 + T7);
                    VADIBL *= T3;
                    dVADIBL_dVg *= T3;
                    dVADIBL_dVb = (dVADIBL_dVb - VADIBL * Param.BSIM4pdiblb)
                                * T3;
                    dVADIBL_dVd *= T3;
                }
                else
                {
                    T4 = 1.0 / (0.8 + T7);
                    T3 = (17.0 + 20.0 * T7) * T4;
                    dVADIBL_dVg *= T3;
                    dVADIBL_dVb = dVADIBL_dVb * T3
                                - VADIBL * Param.BSIM4pdiblb * T4 * T4;
                    dVADIBL_dVd *= T3;
                    VADIBL *= T3;
                }

                dVADIBL_dVg = dVADIBL_dVg * PvagTerm + VADIBL * dPvagTerm_dVg;
                dVADIBL_dVb = dVADIBL_dVb * PvagTerm + VADIBL * dPvagTerm_dVb;
                dVADIBL_dVd = dVADIBL_dVd * PvagTerm + VADIBL * dPvagTerm_dVd;
                VADIBL *= PvagTerm;
            }
            else
            {
                VADIBL = MAX_EXP;
                dVADIBL_dVd = dVADIBL_dVg = dVADIBL_dVb = 0.0;
            }

            /* Calculate Va */
            Va = Vasat + VACLM;
            dVa_dVg = dVasat_dVg + dVACLM_dVg;
            dVa_dVb = dVasat_dVb + dVACLM_dVb;
            dVa_dVd = dVasat_dVd + dVACLM_dVd;

            /* Calculate VADITS */
            T0 = Param.BSIM4pditsd * Vds;
            if (T0 > EXP_THRESHOLD)
            {
                T1 = MAX_EXP;
                dT1_dVd = 0;
            }
            else
            {
                T1 = Math.Exp(T0);
                dT1_dVd = T1 * Param.BSIM4pditsd;
            }

            if (Param.BSIM4pdits > MIN_EXP)
            {
                T2 = 1.0 + ModelParameters.Pditsl * Leff;
                VADITS = (1.0 + T2 * T1) / Param.BSIM4pdits;
                dVADITS_dVg = VADITS * dFP_dVg;
                dVADITS_dVd = FP * T2 * dT1_dVd / Param.BSIM4pdits;
                VADITS *= FP;
            }
            else
            {
                VADITS = MAX_EXP;
                dVADITS_dVg = dVADITS_dVd = 0;
            }

            /* Calculate VASCBE */
            if ((Param.BSIM4pscbe2 > 0.0) && (Param.BSIM4pscbe1 >= 0.0)) /*4.6.2*/
            {
                if (diffVds > Param.BSIM4pscbe1 * Param.BSIM4litl
                    / EXP_THRESHOLD)
                {
                    T0 = Param.BSIM4pscbe1 * Param.BSIM4litl / diffVds;
                    VASCBE = Leff * Math.Exp(T0) / Param.BSIM4pscbe2;
                    T1 = T0 * VASCBE / diffVds;
                    dVASCBE_dVg = T1 * dVdseff_dVg;
                    dVASCBE_dVd = -T1 * (1.0 - dVdseff_dVd);
                    dVASCBE_dVb = T1 * dVdseff_dVb;
                }
                else
                {
                    VASCBE = MAX_EXP * Leff / Param.BSIM4pscbe2;
                    dVASCBE_dVg = dVASCBE_dVd = dVASCBE_dVb = 0.0;
                }
            }
            else
            {
                VASCBE = MAX_EXP;
                dVASCBE_dVg = dVASCBE_dVd = dVASCBE_dVb = 0.0;
            }

            /* Add DIBL to Ids */
            T9 = diffVds / VADIBL;
            T0 = 1.0 + T9;
            Idsa = Idl * T0;
            dIdsa_dVg = T0 * dIdl_dVg - Idl * (dVdseff_dVg + T9 * dVADIBL_dVg) / VADIBL;
            dIdsa_dVd = T0 * dIdl_dVd + Idl
                      * (1.0 - dVdseff_dVd - T9 * dVADIBL_dVd) / VADIBL;
            dIdsa_dVb = T0 * dIdl_dVb - Idl * (dVdseff_dVb + T9 * dVADIBL_dVb) / VADIBL;

            /* Add DITS to Ids */
            T9 = diffVds / VADITS;
            T0 = 1.0 + T9;
            dIdsa_dVg = T0 * dIdsa_dVg - Idsa * (dVdseff_dVg + T9 * dVADITS_dVg) / VADITS;
            dIdsa_dVd = T0 * dIdsa_dVd + Idsa
                      * (1.0 - dVdseff_dVd - T9 * dVADITS_dVd) / VADITS;
            dIdsa_dVb = T0 * dIdsa_dVb - Idsa * dVdseff_dVb / VADITS;
            Idsa *= T0;

            /* Add CLM to Ids */
            T0 = Math.Log(Va / Vasat);
            dT0_dVg = dVa_dVg / Va - dVasat_dVg / Vasat;
            dT0_dVb = dVa_dVb / Va - dVasat_dVb / Vasat;
            dT0_dVd = dVa_dVd / Va - dVasat_dVd / Vasat;
            T1 = T0 / Cclm;
            T9 = 1.0 + T1;
            dT9_dVg = (dT0_dVg - T1 * dCclm_dVg) / Cclm;
            dT9_dVb = (dT0_dVb - T1 * dCclm_dVb) / Cclm;
            dT9_dVd = (dT0_dVd - T1 * dCclm_dVd) / Cclm;

            dIdsa_dVg = dIdsa_dVg * T9 + Idsa * dT9_dVg;
            dIdsa_dVb = dIdsa_dVb * T9 + Idsa * dT9_dVb;
            dIdsa_dVd = dIdsa_dVd * T9 + Idsa * dT9_dVd;
            Idsa *= T9;

            /* Substrate current begins */
            tmp = Param.BSIM4alpha0 + Param.BSIM4alpha1 * Leff;
            if ((tmp <= 0.0) || (Param.BSIM4beta0 <= 0.0))
            {
                Isub = Gbd = Gbb = Gbg = 0.0;
            }
            else
            {
                T2 = tmp / Leff;
                if (diffVds > Param.BSIM4beta0 / EXP_THRESHOLD)
                {
                    T0 = -Param.BSIM4beta0 / diffVds;
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
                T4 = Idsa * Vdseff;
                Isub = T1 * T4;
                Gbg = T1 * (dIdsa_dVg * Vdseff + Idsa * dVdseff_dVg)
                    + T4 * dT1_dVg;
                Gbd = T1 * (dIdsa_dVd * Vdseff + Idsa * dVdseff_dVd)
                    + T4 * dT1_dVd;
                Gbb = T1 * (dIdsa_dVb * Vdseff + Idsa * dVdseff_dVb)
                    + T4 * dT1_dVb;

                Gbd += Gbg * dVgsteff_dVd;
                Gbb += Gbg * dVgsteff_dVb;
                Gbg *= dVgsteff_dVg;
                Gbb *= dVbseff_dVb;
            }
            this._csub = Isub;
            this._gbbs = Gbb;
            this._gbgs = Gbg;
            this._gbds = Gbd;

            /* Add SCBE to Ids */
            T9 = diffVds / VASCBE;
            T0 = 1.0 + T9;
            Ids = Idsa * T0;

            Gm = T0 * dIdsa_dVg - Idsa
               * (dVdseff_dVg + T9 * dVASCBE_dVg) / VASCBE;
            Gds = T0 * dIdsa_dVd + Idsa
                * (1.0 - dVdseff_dVd - T9 * dVASCBE_dVd) / VASCBE;
            Gmb = T0 * dIdsa_dVb - Idsa
                * (dVdseff_dVb + T9 * dVASCBE_dVb) / VASCBE;


            tmp1 = Gds + Gm * dVgsteff_dVd;
            tmp2 = Gmb + Gm * dVgsteff_dVb;
            tmp3 = Gm;

            Gm = (Ids * dVdseff_dVg + Vdseff * tmp3) * dVgsteff_dVg;
            Gds = Ids * (dVdseff_dVd + dVdseff_dVg * dVgsteff_dVd)
                + Vdseff * tmp1;
            Gmb = (Ids * (dVdseff_dVb + dVdseff_dVg * dVgsteff_dVb)
                + Vdseff * tmp2) * dVbseff_dVb;

            cdrain = Ids * Vdseff;

            /* Source End Velocity Limit  */
            if ((ModelParameters.Vtl.Given) && (ModelParameters.Vtl > 0.0))
            {
                T12 = 1.0 / Leff / CoxeffWovL;
                T11 = T12 / Vgsteff;
                T10 = -T11 / Vgsteff;
                vs = cdrain * T11; /* vs */
                dvs_dVg = Gm * T11 + cdrain * T10 * dVgsteff_dVg;
                dvs_dVd = Gds * T11 + cdrain * T10 * dVgsteff_dVd;
                dvs_dVb = Gmb * T11 + cdrain * T10 * dVgsteff_dVb;
                T0 = 2 * MM;
                T1 = vs / (Param.BSIM4vtl * Param.BSIM4tfactor);
                if (T1 > 0.0)
                {
                    T2 = 1.0 + Math.Exp(T0 * Math.Log(T1));
                    T3 = (T2 - 1.0) * T0 / vs;
                    Fsevl = 1.0 / Math.Exp(Math.Log(T2) / T0);
                    dT2_dVg = T3 * dvs_dVg;
                    dT2_dVd = T3 * dvs_dVd;
                    dT2_dVb = T3 * dvs_dVb;
                    T4 = -1.0 / T0 * Fsevl / T2;
                    dFsevl_dVg = T4 * dT2_dVg;
                    dFsevl_dVd = T4 * dT2_dVd;
                    dFsevl_dVb = T4 * dT2_dVb;
                }
                else
                {
                    Fsevl = 1.0;
                    dFsevl_dVg = 0.0;
                    dFsevl_dVd = 0.0;
                    dFsevl_dVb = 0.0;
                }
                Gm *= Fsevl;
                Gm += cdrain * dFsevl_dVg;
                Gmb *= Fsevl;
                Gmb += cdrain * dFsevl_dVb;
                Gds *= Fsevl;
                Gds += cdrain * dFsevl_dVd;

                cdrain *= Fsevl;
            }

            this._gds = Gds;
            this._gm = Gm;
            this._gmbs = Gmb;
            this._idovVds = Ids;
            if (this._idovVds <= 1.0e-9) this._idovVds = 1.0e-9;

            /* Calculate Rg */
            if ((Parameters.RgateMod > 1) ||
                (Parameters.TrnqsMod != 0) || (Parameters.AcnqsMod != 0))
            {
                T9 = Param.BSIM4xrcrg2 * ModelTemperature.Vtm;
                T0 = T9 * beta;
                dT0_dVd = (dbeta_dVd + dbeta_dVg * dVgsteff_dVd) * T9;
                dT0_dVb = (dbeta_dVb + dbeta_dVg * dVgsteff_dVb) * T9;
                dT0_dVg = dbeta_dVg * T9;

                this._gcrg = Param.BSIM4xrcrg1 * (T0 + Ids);
                this._gcrgd = Param.BSIM4xrcrg1 * (dT0_dVd + tmp1);
                this._gcrgb = Param.BSIM4xrcrg1 * (dT0_dVb + tmp2)
                                  * dVbseff_dVb;
                this._gcrgg = Param.BSIM4xrcrg1 * (dT0_dVg + tmp3)
                                 * dVgsteff_dVg;

                if (Parameters.Nf != 1.0)
                {
                    this._gcrg *= Parameters.Nf;
                    this._gcrgg *= Parameters.Nf;
                    this._gcrgd *= Parameters.Nf;
                    this._gcrgb *= Parameters.Nf;
                }

                if (Parameters.RgateMod == 2)
                {
                    T10 = this._grgeltd * this._grgeltd;
                    T11 = this._grgeltd + this._gcrg;
                    this._gcrg = this._grgeltd * this._gcrg / T11;
                    T12 = T10 / T11 / T11;
                    this._gcrgg *= T12;
                    this._gcrgd *= T12;
                    this._gcrgb *= T12;
                }
                this._gcrgs = -(this._gcrgg + this._gcrgd
                                 + this._gcrgb);
            }


            /* Calculate bias-dependent external S/D resistance */
            if (ModelParameters.RdsMod.Value != 0)
            {   /* Rs(V) */
                T0 = vgs - Param.BSIM4vfbsd;
                T1 = Math.Sqrt(T0 * T0 + 1.0e-4);
                vgs_eff = 0.5 * (T0 + T1);
                dvgs_eff_dvg = vgs_eff / T1;

                T0 = 1.0 + Param.BSIM4prwg * vgs_eff;
                dT0_dvg = -Param.BSIM4prwg / T0 / T0 * dvgs_eff_dvg;
                T1 = -Param.BSIM4prwb * vbs;
                dT1_dvb = -Param.BSIM4prwb;

                T2 = 1.0 / T0 + T1;
                T3 = T2 + Math.Sqrt(T2 * T2 + 0.01);
                dT3_dvg = T3 / (T3 - T2);
                dT3_dvb = dT3_dvg * dT1_dvb;
                dT3_dvg *= dT0_dvg;

                T4 = Param.BSIM4rs0 * 0.5;
                Rs = Param.BSIM4rswmin + T3 * T4;
                dRs_dvg = T4 * dT3_dvg;
                dRs_dvb = T4 * dT3_dvb;

                T0 = 1.0 + this._sourceConductance * Rs;
                this._gstot = this._sourceConductance / T0;
                T0 = -this._gstot * this._gstot;
                dgstot_dvd = 0.0; /* place holder */
                dgstot_dvg = T0 * dRs_dvg;
                dgstot_dvb = T0 * dRs_dvb;
                dgstot_dvs = -(dgstot_dvg + dgstot_dvb + dgstot_dvd);

                /* Rd(V) */
                T0 = vgd - Param.BSIM4vfbsd;
                T1 = Math.Sqrt(T0 * T0 + 1.0e-4);
                vgd_eff = 0.5 * (T0 + T1);
                dvgd_eff_dvg = vgd_eff / T1;

                T0 = 1.0 + Param.BSIM4prwg * vgd_eff;
                dT0_dvg = -Param.BSIM4prwg / T0 / T0 * dvgd_eff_dvg;
                T1 = -Param.BSIM4prwb * vbd;
                dT1_dvb = -Param.BSIM4prwb;

                T2 = 1.0 / T0 + T1;
                T3 = T2 + Math.Sqrt(T2 * T2 + 0.01);
                dT3_dvg = T3 / (T3 - T2);
                dT3_dvb = dT3_dvg * dT1_dvb;
                dT3_dvg *= dT0_dvg;

                T4 = Param.BSIM4rd0 * 0.5;
                Rd = Param.BSIM4rdwmin + T3 * T4;
                dRd_dvg = T4 * dT3_dvg;
                dRd_dvb = T4 * dT3_dvb;

                T0 = 1.0 + this._drainConductance * Rd;
                this._gdtot = this._drainConductance / T0;
                T0 = -this._gdtot * this._gdtot;
                dgdtot_dvs = 0.0;
                dgdtot_dvg = T0 * dRd_dvg;
                dgdtot_dvb = T0 * dRd_dvb;
                dgdtot_dvd = -(dgdtot_dvg + dgdtot_dvb + dgdtot_dvs);

                this._gstotd = vses * dgstot_dvd;
                this._gstotg = vses * dgstot_dvg;
                this._gstots = vses * dgstot_dvs;
                this._gstotb = vses * dgstot_dvb;

                T2 = vdes - vds;
                this._gdtotd = T2 * dgdtot_dvd;
                this._gdtotg = T2 * dgdtot_dvg;
                this._gdtots = T2 * dgdtot_dvs;
                this._gdtotb = T2 * dgdtot_dvb;
            }
            else /* WDLiu: for bypass */
            {
                this._gstot = this._gstotd = this._gstotg = 0.0;
                this._gstots = this._gstotb = 0.0;
                this._gdtot = this._gdtotd = this._gdtotg = 0.0;
                this._gdtots = this._gdtotb = 0.0;
            }

            /* GIDL/GISL Models */

            if (ModelParameters.MtrlMod == 0)
                T0 = 3.0 * toxe;
            else
                T0 = ModelParameters.Epsrsub * toxe / epsrox;

            /* Calculate GIDL current */

            vgs_eff = this._vgs_eff;
            dvgs_eff_dvg = this._dvgs_eff_dvg;
            vgd_eff = this._vgd_eff;
            dvgd_eff_dvg = this._dvgd_eff_dvg;

            if (ModelParameters.GidlMod == 0)
            {

                if (ModelParameters.MtrlMod == 0)
                    T1 = (vds - vgs_eff - Param.BSIM4egidl) / T0;
                else
                    T1 = (vds - vgs_eff - Param.BSIM4egidl + Param.BSIM4vfbsd) / T0;

                if ((Param.BSIM4agidl <= 0.0) || (Param.BSIM4bgidl <= 0.0)
                    || (T1 <= 0.0) || (Param.BSIM4cgidl <= 0.0) || (vbd > 0.0))
                    Igidl = Ggidld = Ggidlg = Ggidlb = 0.0;
                else
                {
                    dT1_dVd = 1.0 / T0;
                    dT1_dVg = -dvgs_eff_dvg * dT1_dVd;
                    T2 = Param.BSIM4bgidl / T1;
                    if (T2 < 100.0)
                    {
                        Igidl = Param.BSIM4agidl * Param.BSIM4weffCJ * T1 * Math.Exp(-T2);
                        T3 = Igidl * (1.0 + T2) / T1;
                        Ggidld = T3 * dT1_dVd;
                        Ggidlg = T3 * dT1_dVg;
                    }
                    else
                    {
                        Igidl = Param.BSIM4agidl * Param.BSIM4weffCJ * 3.720075976e-44;
                        Ggidld = Igidl * dT1_dVd;
                        Ggidlg = Igidl * dT1_dVg;
                        Igidl *= T1;
                    }

                    T4 = vbd * vbd;
                    T5 = -vbd * T4;
                    T6 = Param.BSIM4cgidl + T5;
                    T7 = T5 / T6;
                    T8 = 3.0 * Param.BSIM4cgidl * T4 / T6 / T6;
                    Ggidld = Ggidld * T7 + Igidl * T8;
                    Ggidlg = Ggidlg * T7;
                    Ggidlb = -Igidl * T8;
                    Igidl *= T7;
                }
                this._igidl = Igidl;
                this._ggidld = Ggidld;
                this._ggidlg = Ggidlg;
                this._ggidlb = Ggidlb;
                /* Calculate GISL current  */

                if (ModelParameters.MtrlMod == 0)
                    T1 = (-vds - vgd_eff - Param.BSIM4egisl) / T0;
                else
                    T1 = (-vds - vgd_eff - Param.BSIM4egisl + Param.BSIM4vfbsd) / T0;

                if ((Param.BSIM4agisl <= 0.0) || (Param.BSIM4bgisl <= 0.0)
                    || (T1 <= 0.0) || (Param.BSIM4cgisl <= 0.0) || (vbs > 0.0))
                    Igisl = Ggisls = Ggislg = Ggislb = 0.0;
                else
                {
                    dT1_dVd = 1.0 / T0;
                    dT1_dVg = -dvgd_eff_dvg * dT1_dVd;
                    T2 = Param.BSIM4bgisl / T1;
                    if (T2 < 100.0)
                    {
                        Igisl = Param.BSIM4agisl * Param.BSIM4weffCJ * T1 * Math.Exp(-T2);
                        T3 = Igisl * (1.0 + T2) / T1;
                        Ggisls = T3 * dT1_dVd;
                        Ggislg = T3 * dT1_dVg;
                    }
                    else
                    {
                        Igisl = Param.BSIM4agisl * Param.BSIM4weffCJ * 3.720075976e-44;
                        Ggisls = Igisl * dT1_dVd;
                        Ggislg = Igisl * dT1_dVg;
                        Igisl *= T1;
                    }

                    T4 = vbs * vbs;
                    T5 = -vbs * T4;
                    T6 = Param.BSIM4cgisl + T5;
                    T7 = T5 / T6;
                    T8 = 3.0 * Param.BSIM4cgisl * T4 / T6 / T6;
                    Ggisls = Ggisls * T7 + Igisl * T8;
                    Ggislg = Ggislg * T7;
                    Ggislb = -Igisl * T8;
                    Igisl *= T7;
                }
                this._igisl = Igisl;
                this._ggisls = Ggisls;
                this._ggislg = Ggislg;
                this._ggislb = Ggislb;
            }
            else
            {
                /* v4.7 New Gidl/GISL model */

                /* GISL */
                if (ModelParameters.MtrlMod == 0)
                    T1 = (-vds - Param.BSIM4rgisl * vgd_eff - Param.BSIM4egisl) / T0;
                else
                    T1 = (-vds - Param.BSIM4rgisl * vgd_eff - Param.BSIM4egisl + Param.BSIM4vfbsd) / T0;

                if ((Param.BSIM4agisl <= 0.0) ||
                        (Param.BSIM4bgisl <= 0.0) || (T1 <= 0.0) ||
                        (Param.BSIM4cgisl < 0.0))
                    Igisl = Ggisls = Ggislg = Ggislb = 0.0;
                else
                {
                    dT1_dVd = 1 / T0;
                    dT1_dVg = -Param.BSIM4rgisl * dT1_dVd * dvgd_eff_dvg;
                    T2 = Param.BSIM4bgisl / T1;
                    if (T2 < EXPL_THRESHOLD)
                    {
                        Igisl = Param.BSIM4weffCJ * Param.BSIM4agisl * T1 * Math.Exp(-T2);
                        T3 = Igisl / T1 * (T2 + 1);
                        Ggisls = T3 * dT1_dVd;
                        Ggislg = T3 * dT1_dVg;
                    }
                    else
                    {
                        T3 = Param.BSIM4weffCJ * Param.BSIM4agisl * MIN_EXPL;
                        Igisl = T3 * T1;
                        Ggisls = T3 * dT1_dVd;
                        Ggislg = T3 * dT1_dVg;

                    }
                    T4 = vbs - Param.BSIM4fgisl;

                    if (T4 == 0)
                        T5 = EXPL_THRESHOLD;
                    else
                        T5 = Param.BSIM4kgisl / T4;
                    if (T5 < EXPL_THRESHOLD)
                    {
                        T6 = Math.Exp(T5);
                        Ggislb = -Igisl * T6 * T5 / T4;
                    }
                    else
                    {
                        T6 = MAX_EXPL;
                        Ggislb = 0.0;
                    }
                    Ggisls *= T6;
                    Ggislg *= T6;
                    Igisl *= T6;

                }
                this._igisl = Igisl;
                this._ggisls = Ggisls;
                this._ggislg = Ggislg;
                this._ggislb = Ggislb;
                /* End of GISL */

                /* GIDL */
                if (ModelParameters.MtrlMod == 0)
                    T1 = (vds - Param.BSIM4rgidl * vgs_eff - Param.BSIM4egidl) / T0;
                else
                    T1 = (vds - Param.BSIM4rgidl * vgs_eff - Param.BSIM4egidl + Param.BSIM4vfbsd) / T0;



                if ((Param.BSIM4agidl <= 0.0) ||
                        (Param.BSIM4bgidl <= 0.0) || (T1 <= 0.0) ||
                        (Param.BSIM4cgidl < 0.0))
                    Igidl = Ggidld = Ggidlg = Ggidlb = 0.0;
                else
                {
                    dT1_dVd = 1 / T0;
                    dT1_dVg = -Param.BSIM4rgidl * dT1_dVd * dvgs_eff_dvg;
                    T2 = Param.BSIM4bgidl / T1;
                    if (T2 < EXPL_THRESHOLD)
                    {
                        Igidl = Param.BSIM4weffCJ * Param.BSIM4agidl * T1 * Math.Exp(-T2);
                        T3 = Igidl / T1 * (T2 + 1);
                        Ggidld = T3 * dT1_dVd;
                        Ggidlg = T3 * dT1_dVg;

                    }
                    else
                    {
                        T3 = Param.BSIM4weffCJ * Param.BSIM4agidl * MIN_EXPL;
                        Igidl = T3 * T1;
                        Ggidld = T3 * dT1_dVd;
                        Ggidlg = T3 * dT1_dVg;
                    }
                    T4 = vbd - Param.BSIM4fgidl;
                    if (T4 == 0)
                        T5 = EXPL_THRESHOLD;
                    else
                        T5 = Param.BSIM4kgidl / T4;
                    if (T5 < EXPL_THRESHOLD)
                    {
                        T6 = Math.Exp(T5);
                        Ggidlb = -Igidl * T6 * T5 / T4;
                    }
                    else
                    {
                        T6 = MAX_EXPL;
                        Ggidlb = 0.0;
                    }
                    Ggidld *= T6;
                    Ggidlg *= T6;
                    Igidl *= T6;
                }
                this._igidl = Igidl;
                this._ggidld = Ggidld;
                this._ggidlg = Ggidlg;
                this._ggidlb = Ggidlb;
                /* End of New GIDL */
            }
            /*End of Gidl*/



            /* Calculate gate tunneling current */
            if ((ModelParameters.IgcMod != 0) || (ModelParameters.IgbMod != 0))
            {
                Vfb = this._vfbzb;
                V3 = Vfb - Vgs_eff + Vbseff - DELTA_3;
                if (Vfb <= 0.0)
                    T0 = Math.Sqrt(V3 * V3 - 4.0 * DELTA_3 * Vfb);
                else
                    T0 = Math.Sqrt(V3 * V3 + 4.0 * DELTA_3 * Vfb);
                T1 = 0.5 * (1.0 + V3 / T0);
                Vfbeff = Vfb - 0.5 * (V3 + T0);
                dVfbeff_dVg = T1 * dVgs_eff_dVg;
                dVfbeff_dVb = -T1; /* WDLiu: -No surprise? No. -Good! */

                Voxacc = Vfb - Vfbeff;
                dVoxacc_dVg = -dVfbeff_dVg;
                dVoxacc_dVb = -dVfbeff_dVb;
                if (Voxacc < 0.0) /* WDLiu: Avoiding numerical instability. */
                    Voxacc = dVoxacc_dVg = dVoxacc_dVb = 0.0;

                T0 = 0.5 * Param.BSIM4k1ox;
                T3 = Vgs_eff - Vfbeff - Vbseff - Vgsteff;
                if (Param.BSIM4k1ox == 0.0)
                    Voxdepinv = dVoxdepinv_dVg = dVoxdepinv_dVd
                              = dVoxdepinv_dVb = 0.0;
                else if (T3 < 0.0)
                {
                    Voxdepinv = -T3;
                    dVoxdepinv_dVg = -dVgs_eff_dVg + dVfbeff_dVg
                                   + dVgsteff_dVg;
                    dVoxdepinv_dVd = dVgsteff_dVd;
                    dVoxdepinv_dVb = dVfbeff_dVb + 1.0 + dVgsteff_dVb;
                }
                else
                {
                    T1 = Math.Sqrt(T0 * T0 + T3);
                    T2 = T0 / T1;
                    Voxdepinv = Param.BSIM4k1ox * (T1 - T0);
                    dVoxdepinv_dVg = T2 * (dVgs_eff_dVg - dVfbeff_dVg
                                   - dVgsteff_dVg);
                    dVoxdepinv_dVd = -T2 * dVgsteff_dVd;
                    dVoxdepinv_dVb = -T2 * (dVfbeff_dVb + 1.0 + dVgsteff_dVb);
                }

                Voxdepinv += Vgsteff;
                dVoxdepinv_dVg += dVgsteff_dVg;
                dVoxdepinv_dVd += dVgsteff_dVd;
                dVoxdepinv_dVb += dVgsteff_dVb;
            }

            if (ModelParameters.TempMod < 2)
                tmp = Vtm;
            else /* ModelParameters.TempMod = 2 , 3*/
                tmp = Vtm0;
            if (ModelParameters.IgcMod.Value != 0)
            {
                T0 = tmp * Param.BSIM4nigc;
                if (ModelParameters.IgcMod == 1)
                {
                    VxNVt = (Vgs_eff - ModelParameters.Type * this._vth0) / T0;
                    if (VxNVt > EXP_THRESHOLD)
                    {
                        Vaux = Vgs_eff - ModelParameters.Type * this._vth0;
                        dVaux_dVg = dVgs_eff_dVg;
                        dVaux_dVd = 0.0;
                        dVaux_dVb = 0.0;
                    }
                }
                else if (ModelParameters.IgcMod == 2)
                {
                    VxNVt = (Vgs_eff - this._von) / T0;
                    if (VxNVt > EXP_THRESHOLD)
                    {
                        Vaux = Vgs_eff - this._von;
                        dVaux_dVg = dVgs_eff_dVg;
                        dVaux_dVd = -dVth_dVd;
                        dVaux_dVb = -dVth_dVb;
                    }
                }
                if (VxNVt < -EXP_THRESHOLD)
                {
                    Vaux = T0 * Math.Log(1.0 + MIN_EXP);
                    dVaux_dVg = dVaux_dVd = dVaux_dVb = 0.0;
                }
                else if ((VxNVt >= -EXP_THRESHOLD) && (VxNVt <= EXP_THRESHOLD))
                {
                    ExpVxNVt = Math.Exp(VxNVt);
                    Vaux = T0 * Math.Log(1.0 + ExpVxNVt);
                    dVaux_dVg = ExpVxNVt / (1.0 + ExpVxNVt);
                    if (ModelParameters.IgcMod == 1)
                    {
                        dVaux_dVd = 0.0;
                        dVaux_dVb = 0.0;
                    }
                    else if (ModelParameters.IgcMod == 2)
                    {
                        dVaux_dVd = -dVaux_dVg * dVth_dVd;  /* Synopsys 08/30/2013 modify */
                        dVaux_dVb = -dVaux_dVg * dVth_dVb;  /* Synopsys 08/30/2013 modify */
                    }
                    dVaux_dVg *= dVgs_eff_dVg;
                }

                T2 = Vgs_eff * Vaux;
                dT2_dVg = dVgs_eff_dVg * Vaux + Vgs_eff * dVaux_dVg;
                dT2_dVd = Vgs_eff * dVaux_dVd;
                dT2_dVb = Vgs_eff * dVaux_dVb;

                T11 = Param.BSIM4Aechvb;
                T12 = Param.BSIM4Bechvb;
                T3 = Param.BSIM4aigc * Param.BSIM4cigc
                   - Param.BSIM4bigc;
                T4 = Param.BSIM4bigc * Param.BSIM4cigc;
                T5 = T12 * (Param.BSIM4aigc + T3 * Voxdepinv
                   - T4 * Voxdepinv * Voxdepinv);

                if (T5 > EXP_THRESHOLD)
                {
                    T6 = MAX_EXP;
                    dT6_dVg = dT6_dVd = dT6_dVb = 0.0;
                }
                else if (T5 < -EXP_THRESHOLD)
                {
                    T6 = MIN_EXP;
                    dT6_dVg = dT6_dVd = dT6_dVb = 0.0;
                }
                else
                {
                    T6 = Math.Exp(T5);
                    dT6_dVg = T6 * T12 * (T3 - 2.0 * T4 * Voxdepinv);
                    dT6_dVd = dT6_dVg * dVoxdepinv_dVd;
                    dT6_dVb = dT6_dVg * dVoxdepinv_dVb;
                    dT6_dVg *= dVoxdepinv_dVg;
                }

                Igc = T11 * T2 * T6;
                dIgc_dVg = T11 * (T2 * dT6_dVg + T6 * dT2_dVg);
                dIgc_dVd = T11 * (T2 * dT6_dVd + T6 * dT2_dVd);
                dIgc_dVb = T11 * (T2 * dT6_dVb + T6 * dT2_dVb);

                if (ModelParameters.Pigcd.Given)
                {
                    Pigcd = Param.BSIM4pigcd;
                    dPigcd_dVg = dPigcd_dVd = dPigcd_dVb = 0.0;
                }
                else
                {  /* T11 = Param.BSIM4Bechvb * toxe; v4.7 */
                    T11 = -Param.BSIM4Bechvb;
                    T12 = Vgsteff + 1.0e-20;
                    T13 = T11 / T12 / T12;
                    T14 = -T13 / T12;
                    Pigcd = T13 * (1.0 - 0.5 * Vdseff / T12);
                    dPigcd_dVg = T14 * (2.0 + 0.5 * (dVdseff_dVg
                                - 3.0 * Vdseff / T12));
                    dPigcd_dVd = 0.5 * T14 * dVdseff_dVd;
                    dPigcd_dVb = 0.5 * T14 * dVdseff_dVb;
                }

                T7 = -Pigcd * Vdseff; /* bugfix */
                dT7_dVg = -Vdseff * dPigcd_dVg - Pigcd * dVdseff_dVg;
                dT7_dVd = -Vdseff * dPigcd_dVd - Pigcd * dVdseff_dVd + dT7_dVg * dVgsteff_dVd;
                dT7_dVb = -Vdseff * dPigcd_dVb - Pigcd * dVdseff_dVb + dT7_dVg * dVgsteff_dVb;
                dT7_dVg *= dVgsteff_dVg;
                /*dT7_dVb *= dVbseff_dVb;*/ /* Synopsys, 2013/08/30 */
                T8 = T7 * T7 + 2.0e-4;
                dT8_dVg = 2.0 * T7;
                dT8_dVd = dT8_dVg * dT7_dVd;
                dT8_dVb = dT8_dVg * dT7_dVb;
                dT8_dVg *= dT7_dVg;

                if (T7 > EXP_THRESHOLD)
                {
                    T9 = MAX_EXP;
                    dT9_dVg = dT9_dVd = dT9_dVb = 0.0;
                }
                else if (T7 < -EXP_THRESHOLD)
                {
                    T9 = MIN_EXP;
                    dT9_dVg = dT9_dVd = dT9_dVb = 0.0;
                }
                else
                {
                    T9 = Math.Exp(T7);
                    dT9_dVg = T9 * dT7_dVg;
                    dT9_dVd = T9 * dT7_dVd;
                    dT9_dVb = T9 * dT7_dVb;
                }

                T0 = T8 * T8;
                T1 = T9 - 1.0 + 1.0e-4;
                T10 = (T1 - T7) / T8;
                dT10_dVg = (dT9_dVg - dT7_dVg - T10 * dT8_dVg) / T8;
                dT10_dVd = (dT9_dVd - dT7_dVd - T10 * dT8_dVd) / T8;
                dT10_dVb = (dT9_dVb - dT7_dVb - T10 * dT8_dVb) / T8;

                Igcs = Igc * T10;
                dIgcs_dVg = dIgc_dVg * T10 + Igc * dT10_dVg;
                dIgcs_dVd = dIgc_dVd * T10 + Igc * dT10_dVd;
                dIgcs_dVb = dIgc_dVb * T10 + Igc * dT10_dVb;

                T1 = T9 - 1.0 - 1.0e-4;
                T10 = (T7 * T9 - T1) / T8;
                dT10_dVg = (dT7_dVg * T9 + (T7 - 1.0) * dT9_dVg
                         - T10 * dT8_dVg) / T8;
                dT10_dVd = (dT7_dVd * T9 + (T7 - 1.0) * dT9_dVd
                         - T10 * dT8_dVd) / T8;
                dT10_dVb = (dT7_dVb * T9 + (T7 - 1.0) * dT9_dVb
                         - T10 * dT8_dVb) / T8;
                Igcd = Igc * T10;
                dIgcd_dVg = dIgc_dVg * T10 + Igc * dT10_dVg;
                dIgcd_dVd = dIgc_dVd * T10 + Igc * dT10_dVd;
                dIgcd_dVb = dIgc_dVb * T10 + Igc * dT10_dVb;

                this._igcs = Igcs;
                this._gIgcsg = dIgcs_dVg;
                this._gIgcsd = dIgcs_dVd;
                this._gIgcsb = dIgcs_dVb * dVbseff_dVb;
                this._igcd = Igcd;
                this._gIgcdg = dIgcd_dVg;
                this._gIgcdd = dIgcd_dVd;
                this._gIgcdb = dIgcd_dVb * dVbseff_dVb;

                T0 = vgs - (Param.BSIM4vfbsd + Param.BSIM4vfbsdoff);
                vgs_eff = Math.Sqrt(T0 * T0 + 1.0e-4);
                dvgs_eff_dvg = T0 / vgs_eff;

                T2 = vgs * vgs_eff;
                dT2_dVg = vgs * dvgs_eff_dvg + vgs_eff;
                T11 = Param.BSIM4AechvbEdgeS;
                T12 = Param.BSIM4BechvbEdge;
                T3 = Param.BSIM4aigs * Param.BSIM4cigs
                   - Param.BSIM4bigs;
                T4 = Param.BSIM4bigs * Param.BSIM4cigs;
                T5 = T12 * (Param.BSIM4aigs + T3 * vgs_eff
                   - T4 * vgs_eff * vgs_eff);
                if (T5 > EXP_THRESHOLD)
                {
                    T6 = MAX_EXP;
                    dT6_dVg = 0.0;
                }
                else if (T5 < -EXP_THRESHOLD)
                {
                    T6 = MIN_EXP;
                    dT6_dVg = 0.0;
                }
                else
                {
                    T6 = Math.Exp(T5);
                    dT6_dVg = T6 * T12 * (T3 - 2.0 * T4 * vgs_eff)
                            * dvgs_eff_dvg;
                }
                Igs = T11 * T2 * T6;
                dIgs_dVg = T11 * (T2 * dT6_dVg + T6 * dT2_dVg);
                dIgs_dVs = -dIgs_dVg;


                T0 = vgd - (Param.BSIM4vfbsd + Param.BSIM4vfbsdoff);
                vgd_eff = Math.Sqrt(T0 * T0 + 1.0e-4);
                dvgd_eff_dvg = T0 / vgd_eff;

                T2 = vgd * vgd_eff;
                dT2_dVg = vgd * dvgd_eff_dvg + vgd_eff;
                T11 = Param.BSIM4AechvbEdgeD;
                T3 = Param.BSIM4aigd * Param.BSIM4cigd
                   - Param.BSIM4bigd;
                T4 = Param.BSIM4bigd * Param.BSIM4cigd;
                T5 = T12 * (Param.BSIM4aigd + T3 * vgd_eff
                   - T4 * vgd_eff * vgd_eff);
                if (T5 > EXP_THRESHOLD)
                {
                    T6 = MAX_EXP;
                    dT6_dVg = 0.0;
                }
                else if (T5 < -EXP_THRESHOLD)
                {
                    T6 = MIN_EXP;
                    dT6_dVg = 0.0;
                }
                else
                {
                    T6 = Math.Exp(T5);
                    dT6_dVg = T6 * T12 * (T3 - 2.0 * T4 * vgd_eff)
                            * dvgd_eff_dvg;
                }
                Igd = T11 * T2 * T6;
                dIgd_dVg = T11 * (T2 * dT6_dVg + T6 * dT2_dVg);
                dIgd_dVd = -dIgd_dVg;

                this._igs = Igs;
                this._gIgsg = dIgs_dVg;
                this._gIgss = dIgs_dVs;
                this._igd = Igd;
                this._gIgdg = dIgd_dVg;
                this._gIgdd = dIgd_dVd;
            }
            else
            {
                this._igcs = this._gIgcsg = this._gIgcsd
                                 = this._gIgcsb = 0.0;
                this._igcd = this._gIgcdg = this._gIgcdd
                                      = this._gIgcdb = 0.0;
                this._igs = this._gIgsg = this._gIgss = 0.0;
                this._igd = this._gIgdg = this._gIgdd = 0.0;
            }

            if (ModelParameters.IgbMod.Value != 0)
            {
                T0 = tmp * Param.BSIM4nigbacc;
                T1 = -Vgs_eff + Vbseff + Vfb;
                VxNVt = T1 / T0;
                if (VxNVt > EXP_THRESHOLD)
                {
                    Vaux = T1;
                    dVaux_dVg = -dVgs_eff_dVg;
                    dVaux_dVb = 1.0;
                }
                else if (VxNVt < -EXP_THRESHOLD)
                {
                    Vaux = T0 * Math.Log(1.0 + MIN_EXP);
                    dVaux_dVg = dVaux_dVb = 0.0;
                }
                else
                {
                    ExpVxNVt = Math.Exp(VxNVt);
                    Vaux = T0 * Math.Log(1.0 + ExpVxNVt);
                    dVaux_dVb = ExpVxNVt / (1.0 + ExpVxNVt);
                    dVaux_dVg = -dVaux_dVb * dVgs_eff_dVg;
                }

                T2 = (Vgs_eff - Vbseff) * Vaux;
                dT2_dVg = dVgs_eff_dVg * Vaux + (Vgs_eff - Vbseff) * dVaux_dVg;
                dT2_dVb = -Vaux + (Vgs_eff - Vbseff) * dVaux_dVb;

                T11 = 4.97232e-7 * Param.BSIM4weff
                    * Param.BSIM4leff * Param.BSIM4ToxRatio;
                T12 = -7.45669e11 * toxe;
                T3 = Param.BSIM4aigbacc * Param.BSIM4cigbacc
                   - Param.BSIM4bigbacc;
                T4 = Param.BSIM4bigbacc * Param.BSIM4cigbacc;
                T5 = T12 * (Param.BSIM4aigbacc + T3 * Voxacc
                   - T4 * Voxacc * Voxacc);

                if (T5 > EXP_THRESHOLD)
                {
                    T6 = MAX_EXP;
                    dT6_dVg = dT6_dVb = 0.0;
                }
                else if (T5 < -EXP_THRESHOLD)
                {
                    T6 = MIN_EXP;
                    dT6_dVg = dT6_dVb = 0.0;
                }
                else
                {
                    T6 = Math.Exp(T5);
                    dT6_dVg = T6 * T12 * (T3 - 2.0 * T4 * Voxacc);
                    dT6_dVb = dT6_dVg * dVoxacc_dVb;
                    dT6_dVg *= dVoxacc_dVg;
                }

                Igbacc = T11 * T2 * T6;
                dIgbacc_dVg = T11 * (T2 * dT6_dVg + T6 * dT2_dVg);
                dIgbacc_dVb = T11 * (T2 * dT6_dVb + T6 * dT2_dVb);


                T0 = tmp * Param.BSIM4nigbinv;
                T1 = Voxdepinv - Param.BSIM4eigbinv;
                VxNVt = T1 / T0;
                if (VxNVt > EXP_THRESHOLD)
                {
                    Vaux = T1;
                    dVaux_dVg = dVoxdepinv_dVg;
                    dVaux_dVd = dVoxdepinv_dVd;
                    dVaux_dVb = dVoxdepinv_dVb;
                }
                else if (VxNVt < -EXP_THRESHOLD)
                {
                    Vaux = T0 * Math.Log(1.0 + MIN_EXP);
                    dVaux_dVg = dVaux_dVd = dVaux_dVb = 0.0;
                }
                else
                {
                    ExpVxNVt = Math.Exp(VxNVt);
                    Vaux = T0 * Math.Log(1.0 + ExpVxNVt);
                    dVaux_dVg = ExpVxNVt / (1.0 + ExpVxNVt);
                    dVaux_dVd = dVaux_dVg * dVoxdepinv_dVd;
                    dVaux_dVb = dVaux_dVg * dVoxdepinv_dVb;
                    dVaux_dVg *= dVoxdepinv_dVg;
                }

                T2 = (Vgs_eff - Vbseff) * Vaux;
                dT2_dVg = dVgs_eff_dVg * Vaux + (Vgs_eff - Vbseff) * dVaux_dVg;
                dT2_dVd = (Vgs_eff - Vbseff) * dVaux_dVd;
                dT2_dVb = -Vaux + (Vgs_eff - Vbseff) * dVaux_dVb;

                T11 *= 0.75610;
                T12 *= 1.31724;
                T3 = Param.BSIM4aigbinv * Param.BSIM4cigbinv
                   - Param.BSIM4bigbinv;
                T4 = Param.BSIM4bigbinv * Param.BSIM4cigbinv;
                T5 = T12 * (Param.BSIM4aigbinv + T3 * Voxdepinv
                   - T4 * Voxdepinv * Voxdepinv);

                if (T5 > EXP_THRESHOLD)
                {
                    T6 = MAX_EXP;
                    dT6_dVg = dT6_dVd = dT6_dVb = 0.0;
                }
                else if (T5 < -EXP_THRESHOLD)
                {
                    T6 = MIN_EXP;
                    dT6_dVg = dT6_dVd = dT6_dVb = 0.0;
                }
                else
                {
                    T6 = Math.Exp(T5);
                    dT6_dVg = T6 * T12 * (T3 - 2.0 * T4 * Voxdepinv);
                    dT6_dVd = dT6_dVg * dVoxdepinv_dVd;
                    dT6_dVb = dT6_dVg * dVoxdepinv_dVb;
                    dT6_dVg *= dVoxdepinv_dVg;
                }

                Igbinv = T11 * T2 * T6;
                dIgbinv_dVg = T11 * (T2 * dT6_dVg + T6 * dT2_dVg);
                dIgbinv_dVd = T11 * (T2 * dT6_dVd + T6 * dT2_dVd);
                dIgbinv_dVb = T11 * (T2 * dT6_dVb + T6 * dT2_dVb);

                this._igb = Igbinv + Igbacc;
                this._gIgbg = dIgbinv_dVg + dIgbacc_dVg;
                this._gIgbd = dIgbinv_dVd;
                this._gIgbb = (dIgbinv_dVb + dIgbacc_dVb) * dVbseff_dVb;
            }
            else
            {
                this._igb = this._gIgbg = this._gIgbd
                              = this._gIgbs = this._gIgbb = 0.0;
            } /* End of Gate current */

            if (Parameters.Nf != 1.0)
            {
                cdrain *= Parameters.Nf;
                this._gds *= Parameters.Nf;
                this._gm *= Parameters.Nf;
                this._gmbs *= Parameters.Nf;
                this._idovVds *= Parameters.Nf;

                this._gbbs *= Parameters.Nf;
                this._gbgs *= Parameters.Nf;
                this._gbds *= Parameters.Nf;
                this._csub *= Parameters.Nf;

                this._igidl *= Parameters.Nf;
                this._ggidld *= Parameters.Nf;
                this._ggidlg *= Parameters.Nf;
                this._ggidlb *= Parameters.Nf;

                this._igisl *= Parameters.Nf;
                this._ggisls *= Parameters.Nf;
                this._ggislg *= Parameters.Nf;
                this._ggislb *= Parameters.Nf;

                this._igcs *= Parameters.Nf;
                this._gIgcsg *= Parameters.Nf;
                this._gIgcsd *= Parameters.Nf;
                this._gIgcsb *= Parameters.Nf;
                this._igcd *= Parameters.Nf;
                this._gIgcdg *= Parameters.Nf;
                this._gIgcdd *= Parameters.Nf;
                this._gIgcdb *= Parameters.Nf;

                this._igs *= Parameters.Nf;
                this._gIgsg *= Parameters.Nf;
                this._gIgss *= Parameters.Nf;
                this._igd *= Parameters.Nf;
                this._gIgdg *= Parameters.Nf;
                this._gIgdd *= Parameters.Nf;

                this._igb *= Parameters.Nf;
                this._gIgbg *= Parameters.Nf;
                this._gIgbd *= Parameters.Nf;
                this._gIgbb *= Parameters.Nf;
            }

            this._ggidls = -(this._ggidld + this._ggidlg
                              + this._ggidlb);
            this._ggisld = -(this._ggisls + this._ggislg
                              + this._ggislb);
            this._gIgbs = -(this._gIgbg + this._gIgbd
                             + this._gIgbb);
            this._gIgcss = -(this._gIgcsg + this._gIgcsd
                              + this._gIgcsb);
            this._gIgcds = -(this._gIgcdg + this._gIgcdd
                              + this._gIgcdb);
            this._cd = cdrain;


            /* Calculations for noise analysis */

            if (ModelParameters.TnoiMod == 0)
            {
                Abulk = Abulk0 * Param.BSIM4abulkCVfactor;
                Vdsat = Vgsteff / Abulk;
                T0 = Vdsat - Vds - DELTA_4;
                T1 = Math.Sqrt(T0 * T0 + 4.0 * DELTA_4 * Vdsat);
                if (T0 >= 0.0)
                    Vdseff = Vdsat - 0.5 * (T0 + T1);
                else
                {
                    T3 = (DELTA_4 + DELTA_4) / (T1 - T0);
                    T4 = 1.0 - T3;
                    T5 = Vdsat * T3 / (T1 - T0);
                    Vdseff = Vdsat * T4;
                }
                if (Vds == 0.0)
                    Vdseff = 0.0;

                T0 = Abulk * Vdseff;
                T1 = 12.0 * (Vgsteff - 0.5 * T0 + 1.0e-20);
                T2 = Vdseff / T1;
                T3 = T0 * T2;
                this._qinv = Coxeff * Param.BSIM4weffCV * Parameters.Nf
                                * Param.BSIM4leffCV
                                * (Vgsteff - 0.5 * T0 + Abulk * T3);
            }
            else if (ModelParameters.TnoiMod == 2)
            {
                this._noiGd0 = Parameters.Nf * beta * Vgsteff / (1.0 + gche * Rds);
            }

            /*
             *  BSIM4 C-V begins
             */

            if ((ModelParameters.Xpart < 0) || (!ChargeComputationNeeded))
            {
                qgate = qdrn = qsrc = qbulk = 0.0;
                this._cggb = this._cgsb = this._cgdb = 0.0;
                this._cdgb = this._cdsb = this._cddb = 0.0;
                this._cbgb = this._cbsb = this._cbdb = 0.0;
                this._csgb = this._cssb = this._csdb = 0.0;
                this._cgbb = this._csbb = this._cdbb = this._cbbb = 0.0;
                this._cqdb = this._cqsb = this._cqgb
                                = this._cqbb = 0.0;
                this._gtau = 0.0;
                goto finished;
            }
            else if (ModelParameters.CapMod == 0)
            {
                if (Vbseff < 0.0)
                {
                    VbseffCV = Vbs; /*4.6.2*/
                    dVbseffCV_dVb = 1.0;
                }
                else
                {
                    VbseffCV = Param.BSIM4phi - Phis;
                    dVbseffCV_dVb = -dPhis_dVb * dVbseff_dVb; /*4.6.2*/
                }

                Vfb = Param.BSIM4vfbcv;
                Vth = Vfb + Param.BSIM4phi + Param.BSIM4k1ox * sqrtPhis;
                Vgst = Vgs_eff - Vth;
                dVth_dVb = Param.BSIM4k1ox * dsqrtPhis_dVb * dVbseff_dVb; /*4.6.2*/
                dVgst_dVb = -dVth_dVb;
                dVgst_dVg = dVgs_eff_dVg;

                CoxWL = ModelTemperature.Coxe * Param.BSIM4weffCV
                      * Param.BSIM4leffCV * Parameters.Nf;
                Arg1 = Vgs_eff - VbseffCV - Vfb;

                if (Arg1 <= 0.0)
                {
                    qgate = CoxWL * Arg1;
                    qbulk = -qgate;
                    qdrn = 0.0;

                    this._cggb = CoxWL * dVgs_eff_dVg;
                    this._cgdb = 0.0;
                    this._cgsb = CoxWL * (dVbseffCV_dVb - dVgs_eff_dVg);

                    this._cdgb = 0.0;
                    this._cddb = 0.0;
                    this._cdsb = 0.0;

                    this._cbgb = -CoxWL * dVgs_eff_dVg;
                    this._cbdb = 0.0;
                    this._cbsb = -this._cgsb;
                } /* Arg1 <= 0.0, end of accumulation */
                else if (Vgst <= 0.0)
                {
                    T1 = 0.5 * Param.BSIM4k1ox;
                    T2 = Math.Sqrt(T1 * T1 + Arg1);
                    qgate = CoxWL * Param.BSIM4k1ox * (T2 - T1);
                    qbulk = -qgate;
                    qdrn = 0.0;

                    T0 = CoxWL * T1 / T2;
                    this._cggb = T0 * dVgs_eff_dVg;
                    this._cgdb = 0.0;
                    this._cgsb = T0 * (dVbseffCV_dVb - dVgs_eff_dVg);

                    this._cdgb = 0.0;
                    this._cddb = 0.0;
                    this._cdsb = 0.0;

                    this._cbgb = -this._cggb;
                    this._cbdb = 0.0;
                    this._cbsb = -this._cgsb;
                } /* Vgst <= 0.0, end of depletion */
                else
                {
                    One_Third_CoxWL = CoxWL / 3.0;
                    Two_Third_CoxWL = 2.0 * One_Third_CoxWL;

                    AbulkCV = Abulk0 * Param.BSIM4abulkCVfactor;
                    dAbulkCV_dVb = Param.BSIM4abulkCVfactor * dAbulk0_dVb * dVbseff_dVb;

                    dVdsat_dVg = 1.0 / AbulkCV;  /*4.6.2*/
                    Vdsat = Vgst * dVdsat_dVg;
                    dVdsat_dVb = -(Vdsat * dAbulkCV_dVb + dVth_dVb) * dVdsat_dVg;

                    if (ModelParameters.Xpart > 0.5)
                    {   /* 0/100 Charge partition model */
                        if (Vdsat <= Vds)
                        {   /* saturation region */
                            T1 = Vdsat / 3.0;
                            qgate = CoxWL * (Vgs_eff - Vfb
                                  - Param.BSIM4phi - T1);
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
                                  - Param.BSIM4phi - 0.5 * (Vds - T3));
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
                        }
                    }
                    else if (ModelParameters.Xpart < 0.5)
                    {   /* 40/60 Charge partition model */
                        if (Vds >= Vdsat)
                        {   /* saturation region */
                            T1 = Vdsat / 3.0;
                            qgate = CoxWL * (Vgs_eff - Vfb
                                  - Param.BSIM4phi - T1);
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
                        }
                        else
                        {   /* linear region  */
                            Alphaz = Vgst / Vdsat;
                            T1 = 2.0 * Vdsat - Vds;
                            T2 = Vds / (3.0 * T1);
                            T3 = T2 * Vds;
                            T9 = 0.25 * CoxWL;
                            T4 = T9 * Alphaz;
                            qgate = CoxWL * (Vgs_eff - Vfb - Param.BSIM4phi
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
                            T12 = (-T7 * dAlphaz_dVg - T0 * dVdsat_dVg) * dVgs_eff_dVg
                                    - this._cdgb;  /*4.6.2*/
                            T11 = -T7 * dAlphaz_dVb - T10 - T0 * dVdsat_dVb;
                            T10 = -4.0 * T4 * (T2 - 0.5 + 0.5 * T5)
                                - this._cddb;
                            tmp = -(T10 + T11 + T12);

                            this._cbgb = -(this._cggb
                                            + this._cdgb + T12);
                            this._cbdb = -(this._cgdb
                                            + this._cddb + T10);
                            this._cbsb = -(this._cgsb
                                            + this._cdsb + tmp);
                        }
                    }
                    else
                    {   /* 50/50 partitioning */
                        if (Vds >= Vdsat)
                        {   /* saturation region */
                            T1 = Vdsat / 3.0;
                            qgate = CoxWL * (Vgs_eff - Vfb
                                  - Param.BSIM4phi - T1);
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
                        }
                        else
                        {   /* linear region */
                            Alphaz = Vgst / Vdsat;
                            T1 = 2.0 * Vdsat - Vds;
                            T2 = Vds / (3.0 * T1);
                            T3 = T2 * Vds;
                            T9 = 0.25 * CoxWL;
                            T4 = T9 * Alphaz;
                            qgate = CoxWL * (Vgs_eff - Vfb - Param.BSIM4phi
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
                        } /* end of linear region */
                    } /* end of 50/50 partition */
                } /* end of inversion */
            } /* end of capMod=0 */
            else
            {
                if (Vbseff < 0.0)
                {
                    VbseffCV = Vbseff;
                    dVbseffCV_dVb = 1.0;
                }
                else
                {
                    VbseffCV = Param.BSIM4phi - Phis;
                    dVbseffCV_dVb = -dPhis_dVb;
                }

                CoxWL = ModelTemperature.Coxe * Param.BSIM4weffCV
                      * Param.BSIM4leffCV * Parameters.Nf;

                if (ModelParameters.CvchargeMod == 0)
                {
                    /* Seperate VgsteffCV with noff and voffcv */
                    noff = n * Param.BSIM4noff;
                    dnoff_dVd = Param.BSIM4noff * dn_dVd;
                    dnoff_dVb = Param.BSIM4noff * dn_dVb;
                    T0 = Vtm * noff;
                    voffcv = Param.BSIM4voffcv;
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
                    }
                    /* End of VgsteffCV for cvchargeMod = 0 */
                }
                else
                {
                    T0 = n * Vtm;
                    T1 = Param.BSIM4mstarcv * Vgst;
                    T2 = T1 / T0;
                    if (T2 > EXP_THRESHOLD)
                    {
                        T10 = T1;
                        dT10_dVg = Param.BSIM4mstarcv * dVgs_eff_dVg;
                        dT10_dVd = -dVth_dVd * Param.BSIM4mstarcv;
                        dT10_dVb = -dVth_dVb * Param.BSIM4mstarcv;
                    }
                    else if (T2 < -EXP_THRESHOLD)
                    {
                        T10 = Vtm * Math.Log(1.0 + MIN_EXP);
                        dT10_dVg = 0.0;
                        dT10_dVd = T10 * dn_dVd;
                        dT10_dVb = T10 * dn_dVb;
                        T10 *= n;
                    }
                    else
                    {
                        ExpVgst = Math.Exp(T2);
                        T3 = Vtm * Math.Log(1.0 + ExpVgst);
                        T10 = n * T3;
                        dT10_dVg = Param.BSIM4mstarcv * ExpVgst / (1.0 + ExpVgst);
                        dT10_dVb = T3 * dn_dVb - dT10_dVg * (dVth_dVb + Vgst * dn_dVb / n);
                        dT10_dVd = T3 * dn_dVd - dT10_dVg * (dVth_dVd + Vgst * dn_dVd / n);
                        dT10_dVg *= dVgs_eff_dVg;
                    }

                    T1 = Param.BSIM4voffcbncv - (1.0 - Param.BSIM4mstarcv) * Vgst;
                    T2 = T1 / T0;
                    if (T2 < -EXP_THRESHOLD)
                    {
                        T3 = ModelTemperature.Coxe * MIN_EXP / Param.BSIM4cdep0;
                        T9 = Param.BSIM4mstarcv + T3 * n;
                        dT9_dVg = 0.0;
                        dT9_dVd = dn_dVd * T3;
                        dT9_dVb = dn_dVb * T3;
                    }
                    else if (T2 > EXP_THRESHOLD)
                    {
                        T3 = ModelTemperature.Coxe * MAX_EXP / Param.BSIM4cdep0;
                        T9 = Param.BSIM4mstarcv + T3 * n;
                        dT9_dVg = 0.0;
                        dT9_dVd = dn_dVd * T3;
                        dT9_dVb = dn_dVb * T3;
                    }
                    else
                    {
                        ExpVgst = Math.Exp(T2);
                        T3 = ModelTemperature.Coxe / Param.BSIM4cdep0;
                        T4 = T3 * ExpVgst;
                        T5 = T1 * T4 / T0;
                        T9 = Param.BSIM4mstarcv + n * T4;
                        dT9_dVg = T3 * (Param.BSIM4mstarcv - 1.0) * ExpVgst / Vtm;
                        dT9_dVb = T4 * dn_dVb - dT9_dVg * dVth_dVb - T5 * dn_dVb;
                        dT9_dVd = T4 * dn_dVd - dT9_dVg * dVth_dVd - T5 * dn_dVd;
                        dT9_dVg *= dVgs_eff_dVg;
                    }

                    Vgsteff = T10 / T9;
                    T11 = T9 * T9;
                    dVgsteff_dVg = (T9 * dT10_dVg - T10 * dT9_dVg) / T11;
                    dVgsteff_dVd = (T9 * dT10_dVd - T10 * dT9_dVd) / T11;
                    dVgsteff_dVb = (T9 * dT10_dVb - T10 * dT9_dVb) / T11;
                    /* End of VgsteffCV for cvchargeMod = 1 */
                }


                if (ModelParameters.CapMod == 1)
                {
                    Vfb = this._vfbzb;
                    V3 = Vfb - Vgs_eff + VbseffCV - DELTA_3;
                    if (Vfb <= 0.0)
                        T0 = Math.Sqrt(V3 * V3 - 4.0 * DELTA_3 * Vfb);
                    else
                        T0 = Math.Sqrt(V3 * V3 + 4.0 * DELTA_3 * Vfb);

                    T1 = 0.5 * (1.0 + V3 / T0);
                    Vfbeff = Vfb - 0.5 * (V3 + T0);
                    dVfbeff_dVg = T1 * dVgs_eff_dVg;
                    dVfbeff_dVb = -T1 * dVbseffCV_dVb;
                    Qac0 = CoxWL * (Vfbeff - Vfb);
                    dQac0_dVg = CoxWL * dVfbeff_dVg;
                    dQac0_dVb = CoxWL * dVfbeff_dVb;

                    T0 = 0.5 * Param.BSIM4k1ox;
                    T3 = Vgs_eff - Vfbeff - VbseffCV - Vgsteff;
                    if (Param.BSIM4k1ox == 0.0)
                    {
                        T1 = 0.0;
                        T2 = 0.0;
                    }
                    else if (T3 < 0.0)
                    {
                        T1 = T0 + T3 / Param.BSIM4k1ox;
                        T2 = CoxWL;
                    }
                    else
                    {
                        T1 = Math.Sqrt(T0 * T0 + T3);
                        T2 = CoxWL * T0 / T1;
                    }

                    Qsub0 = CoxWL * Param.BSIM4k1ox * (T1 - T0);

                    dQsub0_dVg = T2 * (dVgs_eff_dVg - dVfbeff_dVg - dVgsteff_dVg);
                    dQsub0_dVd = -T2 * dVgsteff_dVd;
                    dQsub0_dVb = -T2 * (dVfbeff_dVb + dVbseffCV_dVb
                               + dVgsteff_dVb);

                    AbulkCV = Abulk0 * Param.BSIM4abulkCVfactor;
                    dAbulkCV_dVb = Param.BSIM4abulkCVfactor * dAbulk0_dVb;
                    VdsatCV = Vgsteff / AbulkCV;

                    T0 = VdsatCV - Vds - DELTA_4;
                    dT0_dVg = 1.0 / AbulkCV;
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
                        dVdseffCV_dVb = dT0_dVb * (T4 - T5) + T5 * dT1_dVb;
                    }

                    if (Vds == 0.0)
                    {
                        VdseffCV = 0.0;
                        dVdseffCV_dVg = 0.0;
                        dVdseffCV_dVb = 0.0;
                    }

                    T0 = AbulkCV * VdseffCV;
                    T1 = 12.0 * (Vgsteff - 0.5 * T0 + 1.0e-20);
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
                }

                /* Charge-Thickness capMod (CTM) begins */
                else if (ModelParameters.CapMod == 2)
                {
                    V3 = this._vfbzb - Vgs_eff + VbseffCV - DELTA_3;
                    if (this._vfbzb <= 0.0)
                        T0 = Math.Sqrt(V3 * V3 - 4.0 * DELTA_3 * this._vfbzb);
                    else
                        T0 = Math.Sqrt(V3 * V3 + 4.0 * DELTA_3 * this._vfbzb);

                    T1 = 0.5 * (1.0 + V3 / T0);
                    Vfbeff = this._vfbzb - 0.5 * (V3 + T0);
                    dVfbeff_dVg = T1 * dVgs_eff_dVg;
                    dVfbeff_dVb = -T1 * dVbseffCV_dVb;

                    Cox = this._coxp;
                    Tox = 1.0e8 * this._toxp;
                    T0 = (Vgs_eff - VbseffCV - this._vfbzb) / Tox;
                    dT0_dVg = dVgs_eff_dVg / Tox;
                    dT0_dVb = -dVbseffCV_dVb / Tox;

                    tmp = T0 * Param.BSIM4acde;
                    if ((-EXP_THRESHOLD < tmp) && (tmp < EXP_THRESHOLD))
                    {
                        Tcen = Param.BSIM4ldeb * Math.Exp(tmp);
                        dTcen_dVg = Param.BSIM4acde * Tcen;
                        dTcen_dVb = dTcen_dVg * dT0_dVb;
                        dTcen_dVg *= dT0_dVg;
                    }
                    else if (tmp <= -EXP_THRESHOLD)
                    {
                        Tcen = Param.BSIM4ldeb * MIN_EXP;
                        dTcen_dVg = dTcen_dVb = 0.0;
                    }
                    else
                    {
                        Tcen = Param.BSIM4ldeb * MAX_EXP;
                        dTcen_dVg = dTcen_dVb = 0.0;
                    }

                    LINK = 1.0e-3 * this._toxp;
                    V3 = Param.BSIM4ldeb - Tcen - LINK;
                    V4 = Math.Sqrt(V3 * V3 + 4.0 * LINK * Param.BSIM4ldeb);
                    Tcen = Param.BSIM4ldeb - 0.5 * (V3 + V4);
                    T1 = 0.5 * (1.0 + V3 / V4);
                    dTcen_dVg *= T1;
                    dTcen_dVb *= T1;

                    Ccen = epssub / Tcen;
                    T2 = Cox / (Cox + Ccen);
                    Coxeff = T2 * Ccen;
                    T3 = -Ccen / Tcen;
                    dCoxeff_dVg = T2 * T2 * T3;
                    dCoxeff_dVb = dCoxeff_dVg * dTcen_dVb;
                    dCoxeff_dVg *= dTcen_dVg;
                    CoxWLcen = CoxWL * Coxeff / ModelTemperature.Coxe;

                    Qac0 = CoxWLcen * (Vfbeff - this._vfbzb);
                    QovCox = Qac0 / Coxeff;
                    dQac0_dVg = CoxWLcen * dVfbeff_dVg
                              + QovCox * dCoxeff_dVg;
                    dQac0_dVb = CoxWLcen * dVfbeff_dVb
                              + QovCox * dCoxeff_dVb;

                    T0 = 0.5 * Param.BSIM4k1ox;
                    T3 = Vgs_eff - Vfbeff - VbseffCV - Vgsteff;
                    if (Param.BSIM4k1ox == 0.0)
                    {
                        T1 = 0.0;
                        T2 = 0.0;
                    }
                    else if (T3 < 0.0)
                    {
                        T1 = T0 + T3 / Param.BSIM4k1ox;
                        T2 = CoxWLcen;
                    }
                    else
                    {
                        T1 = Math.Sqrt(T0 * T0 + T3);
                        T2 = CoxWLcen * T0 / T1;
                    }

                    Qsub0 = CoxWLcen * Param.BSIM4k1ox * (T1 - T0);
                    QovCox = Qsub0 / Coxeff;
                    dQsub0_dVg = T2 * (dVgs_eff_dVg - dVfbeff_dVg - dVgsteff_dVg)
                               + QovCox * dCoxeff_dVg;
                    dQsub0_dVd = -T2 * dVgsteff_dVd;
                    dQsub0_dVb = -T2 * (dVfbeff_dVb + dVbseffCV_dVb + dVgsteff_dVb)
                               + QovCox * dCoxeff_dVb;

                    /* Gate-bias dependent delta Phis begins */
                    if (Param.BSIM4k1ox <= 0.0)
                    {
                        Denomi = 0.25 * Param.BSIM4moin * Vtm;
                        T0 = 0.5 * Param.BSIM4sqrtPhi;
                    }
                    else
                    {
                        Denomi = Param.BSIM4moin * Vtm
                               * Param.BSIM4k1ox * Param.BSIM4k1ox;
                        T0 = Param.BSIM4k1ox * Param.BSIM4sqrtPhi;
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

                    Tox += Tox; /* WDLiu: Tcen reevaluated below due to different Vgsteff */
                    T0 = (Vgsteff + this._vtfbphi2) / Tox;
                    tmp = Math.Exp(ModelParameters.Bdos * 0.7 * Math.Log(T0));
                    T1 = 1.0 + tmp;
                    T2 = ModelParameters.Bdos * 0.7 * tmp / (T0 * Tox);
                    Tcen = ModelParameters.Ados * 1.9e-9 / T1;
                    dTcen_dVg = -Tcen * T2 / T1;
                    dTcen_dVd = dTcen_dVg * dVgsteff_dVd;
                    dTcen_dVb = dTcen_dVg * dVgsteff_dVb;
                    dTcen_dVg *= dVgsteff_dVg;

                    Ccen = epssub / Tcen;
                    T0 = Cox / (Cox + Ccen);
                    Coxeff = T0 * Ccen;
                    T1 = -Ccen / Tcen;
                    dCoxeff_dVg = T0 * T0 * T1;
                    dCoxeff_dVd = dCoxeff_dVg * dTcen_dVd;
                    dCoxeff_dVb = dCoxeff_dVg * dTcen_dVb;
                    dCoxeff_dVg *= dTcen_dVg;
                    CoxWLcen = CoxWL * Coxeff / ModelTemperature.Coxe;

                    AbulkCV = Abulk0 * Param.BSIM4abulkCVfactor;
                    dAbulkCV_dVb = Param.BSIM4abulkCVfactor * dAbulk0_dVb;
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
                        dVdseffCV_dVb = dT0_dVb * (T4 - T5) + T5 * dT1_dVb;
                    }

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

                    qgate = CoxWLcen * (T1 - T0 * (0.5 - T3));
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

                    if (ModelParameters.Xpart > 0.5)
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
                    else if (ModelParameters.Xpart < 0.5)
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
                }  /* End of CTM */
            }

            this._csgb = -this._cggb - this._cdgb - this._cbgb;
            this._csdb = -this._cgdb - this._cddb - this._cbdb;
            this._cssb = -this._cgsb - this._cdsb - this._cbsb;
            this._cgbb = -this._cgdb - this._cggb - this._cgsb;
            this._cdbb = -this._cddb - this._cdgb - this._cdsb;
            this._cbbb = -this._cbgb - this._cbdb - this._cbsb;
            this._csbb = -this._cgbb - this._cdbb - this._cbbb;
            this._qgate = qgate;
            this._qbulk = qbulk;
            this._qdrn = qdrn;
            this._qsrc = -(qgate + qbulk + qdrn);

            /* NQS begins */
            if ((Parameters.TrnqsMod.Value != 0) || (Parameters.AcnqsMod.Value != 0))
            {
                this._qchqs = qcheq = -(qbulk + qgate);
                this._cqgb = -(this._cggb + this._cbgb);
                this._cqdb = -(this._cgdb + this._cbdb);
                this._cqsb = -(this._cgsb + this._cbsb);
                this._cqbb = -(this._cqgb + this._cqdb
                                + this._cqsb);

                CoxWL = ModelTemperature.Coxe * Param.BSIM4weffCV * Parameters.Nf
                      * Param.BSIM4leffCV;
                T1 = this._gcrg / CoxWL; /* 1 / tau */
                this._gtau = T1 * ScalingFactor;

                if (Parameters.AcnqsMod.Value != 0)
                    this._taunet = 1.0 / T1;

                this._qcheq.Value = qcheq;
                if (Parameters.TrnqsMod.Value != 0)
                    _qcheq.Derive();
            }


        finished:

            /* Calculate junction C-V */
            if (ChargeComputationNeeded)
            {
                czbd = ModelTemperature.DunitAreaTempJctCap * this._adeff; /* bug fix */
                czbs = ModelTemperature.SunitAreaTempJctCap * this._aseff;
                czbdsw = ModelTemperature.DunitLengthSidewallTempJctCap * this._pdeff;
                czbdswg = ModelTemperature.DunitLengthGateSidewallTempJctCap
                        * Param.BSIM4weffCJ * Parameters.Nf;
                czbssw = ModelTemperature.SunitLengthSidewallTempJctCap * this._pseff;
                czbsswg = ModelTemperature.SunitLengthGateSidewallTempJctCap
                        * Param.BSIM4weffCJ * Parameters.Nf;

                MJS = ModelParameters.SbulkJctBotGradingCoeff;
                MJSWS = ModelParameters.SbulkJctSideGradingCoeff;
                MJSWGS = ModelParameters.SbulkJctGateSideGradingCoeff;

                MJD = ModelParameters.DbulkJctBotGradingCoeff;
                MJSWD = ModelParameters.DbulkJctSideGradingCoeff;
                MJSWGD = ModelParameters.DbulkJctGateSideGradingCoeff;

                /* Source Bulk Junction */
                if (vbs_jct == 0.0)
                {
                    this._qbs.Value = 0.0;
                    this._capbs = czbs + czbssw + czbsswg;
                }
                else if (vbs_jct < 0.0)
                {
                    if (czbs > 0.0)
                    {
                        arg = 1.0 - vbs_jct / ModelTemperature.PhiBS;
                        if (MJS == 0.5)
                            sarg = 1.0 / Math.Sqrt(arg);
                        else
                            sarg = Math.Exp(-MJS * Math.Log(arg));
                        this._qbs.Value = ModelTemperature.PhiBS * czbs
                                         * (1.0 - arg * sarg) / (1.0 - MJS);
                        this._capbs = czbs * sarg;
                    }
                    else
                    {
                        this._qbs.Value = 0.0;
                        this._capbs = 0.0;
                    }
                    if (czbssw > 0.0)
                    {
                        arg = 1.0 - vbs_jct / ModelTemperature.PhiBSWS;
                        if (MJSWS == 0.5)
                            sarg = 1.0 / Math.Sqrt(arg);
                        else
                            sarg = Math.Exp(-MJSWS * Math.Log(arg));
                        this._qbs.Value += ModelTemperature.PhiBSWS * czbssw
                                         * (1.0 - arg * sarg) / (1.0 - MJSWS);
                        this._capbs += czbssw * sarg;
                    }
                    if (czbsswg > 0.0)
                    {
                        arg = 1.0 - vbs_jct / ModelTemperature.PhiBSWGS;
                        if (MJSWGS == 0.5)
                            sarg = 1.0 / Math.Sqrt(arg);
                        else
                            sarg = Math.Exp(-MJSWGS * Math.Log(arg));
                        this._qbs.Value += ModelTemperature.PhiBSWGS * czbsswg
                                         * (1.0 - arg * sarg) / (1.0 - MJSWGS);
                        this._capbs += czbsswg * sarg;
                    }

                }
                else
                {
                    T0 = czbs + czbssw + czbsswg;
                    T1 = vbs_jct * (czbs * MJS / ModelTemperature.PhiBS + czbssw * MJSWS
                       / ModelTemperature.PhiBSWS + czbsswg * MJSWGS / ModelTemperature.PhiBSWGS);
                    this._qbs.Value = vbs_jct * (T0 + 0.5 * T1);
                    this._capbs = T0 + T1;
                }

                /* Drain Bulk Junction */
                if (vbd_jct == 0.0)
                {
                    this._qbd.Value = 0.0;
                    this._capbd = czbd + czbdsw + czbdswg;
                }
                else if (vbd_jct < 0.0)
                {
                    if (czbd > 0.0)
                    {
                        arg = 1.0 - vbd_jct / ModelTemperature.PhiBD;
                        if (MJD == 0.5)
                            sarg = 1.0 / Math.Sqrt(arg);
                        else
                            sarg = Math.Exp(-MJD * Math.Log(arg));
                        this._qbd.Value = ModelTemperature.PhiBD * czbd
                                         * (1.0 - arg * sarg) / (1.0 - MJD);
                        this._capbd = czbd * sarg;
                    }
                    else
                    {
                        this._qbd.Value = 0.0;
                        this._capbd = 0.0;
                    }
                    if (czbdsw > 0.0)
                    {
                        arg = 1.0 - vbd_jct / ModelTemperature.PhiBSWD;
                        if (MJSWD == 0.5)
                            sarg = 1.0 / Math.Sqrt(arg);
                        else
                            sarg = Math.Exp(-MJSWD * Math.Log(arg));
                        this._qbd.Value += ModelTemperature.PhiBSWD * czbdsw
                                         * (1.0 - arg * sarg) / (1.0 - MJSWD);
                        this._capbd += czbdsw * sarg;
                    }
                    if (czbdswg > 0.0)
                    {
                        arg = 1.0 - vbd_jct / ModelTemperature.PhiBSWGD;
                        if (MJSWGD == 0.5)
                            sarg = 1.0 / Math.Sqrt(arg);
                        else
                            sarg = Math.Exp(-MJSWGD * Math.Log(arg));
                        this._qbd.Value += ModelTemperature.PhiBSWGD * czbdswg
                                         * (1.0 - arg * sarg) / (1.0 - MJSWGD);
                        this._capbd += czbdswg * sarg;
                    }
                }
                else
                {
                    T0 = czbd + czbdsw + czbdswg;
                    T1 = vbd_jct * (czbd * MJD / ModelTemperature.PhiBD + czbdsw * MJSWD
                       / ModelTemperature.PhiBSWD + czbdswg * MJSWGD / ModelTemperature.PhiBSWGD);
                    this._qbd.Value = vbd_jct * (T0 + 0.5 * T1);
                    this._capbd = T0 + T1;
                }
            }


            /*
             *  check convergence
             */

            if (Parameters.Off || _iteration.Mode != IterationModes.Fix)
            {
                if (Check)
                    _iteration.IsConvergent = false;
            }
            this._vds = vds;
            this._vgs = vgs;
            this._vbs = vbs;
            this._vbd = vbd;
            this._vges = vges;
            this._vgms = vgms;
            this._vdbs = vdbs;
            this._vdbd = vdbd;
            this._vsbs = vsbs;
            this._vses = vses;
            this._vdes = vdes;
            this._qdef = qdef;


            if (!ChargeComputationNeeded)
                goto line850;

            if (Parameters.RgateMod == 3)
            {
                vgdx = vgmd;
                vgsx = vgms;
            }
            else  /* For rgateMod == 0, 1 and 2 */
            {
                vgdx = vgd;
                vgsx = vgs;
            }
            if (ModelParameters.CapMod == 0)
            {
                cgdo = Param.BSIM4cgdo;
                qgdo = Param.BSIM4cgdo * vgdx;
                cgso = Param.BSIM4cgso;
                qgso = Param.BSIM4cgso * vgsx;
            }
            else /* For both capMod == 1 and 2 */
            {
                T0 = vgdx + DELTA_1;
                T1 = Math.Sqrt(T0 * T0 + 4.0 * DELTA_1);
                T2 = 0.5 * (T0 - T1);

                T3 = Param.BSIM4weffCV * Param.BSIM4cgdl;
                T4 = Math.Sqrt(1.0 - 4.0 * T2 / Param.BSIM4ckappad);
                cgdo = Param.BSIM4cgdo + T3 - T3 * (1.0 - 1.0 / T4)
                     * (0.5 - 0.5 * T0 / T1);
                qgdo = (Param.BSIM4cgdo + T3) * vgdx - T3 * (T2
                     + 0.5 * Param.BSIM4ckappad * (T4 - 1.0));

                T0 = vgsx + DELTA_1;
                T1 = Math.Sqrt(T0 * T0 + 4.0 * DELTA_1);
                T2 = 0.5 * (T0 - T1);
                T3 = Param.BSIM4weffCV * Param.BSIM4cgsl;
                T4 = Math.Sqrt(1.0 - 4.0 * T2 / Param.BSIM4ckappas);
                cgso = Param.BSIM4cgso + T3 - T3 * (1.0 - 1.0 / T4)
                     * (0.5 - 0.5 * T0 / T1);
                qgso = (Param.BSIM4cgso + T3) * vgsx - T3 * (T2
                     + 0.5 * Param.BSIM4ckappas * (T4 - 1.0));
            }

            if (Parameters.Nf != 1.0)
            {
                cgdo *= Parameters.Nf;
                cgso *= Parameters.Nf;
                qgdo *= Parameters.Nf;
                qgso *= Parameters.Nf;
            }
            this._cgdo = cgdo;
            this._qgdo = qgdo;
            this._cgso = cgso;
            this._qgso = qgso;

            if (_method != null)
                ag0 = _method.Slope;
            else
                ag0 = 0.0;
            if (this._mode > 0)
            {
                if (Parameters.TrnqsMod == 0)
                {
                    qdrn -= qgdo;
                    if (Parameters.RgateMod == 3)
                    {
                        gcgmgmb = (cgdo + cgso + Param.BSIM4cgbo) * ag0;
                        gcgmdb = -cgdo * ag0;
                        gcgmsb = -cgso * ag0;
                        gcgmbb = -Param.BSIM4cgbo * ag0;

                        gcdgmb = gcgmdb;
                        gcsgmb = gcgmsb;
                        gcbgmb = gcgmbb;

                        gcggb = this._cggb * ag0;
                        gcgdb = this._cgdb * ag0;
                        gcgsb = this._cgsb * ag0;
                        gcgbb = -(gcggb + gcgdb + gcgsb);

                        gcdgb = this._cdgb * ag0;
                        gcsgb = -(this._cggb + this._cbgb
                              + this._cdgb) * ag0;
                        gcbgb = this._cbgb * ag0;

                        qgmb = Param.BSIM4cgbo * vgmb;
                        qgmid = qgdo + qgso + qgmb;
                        qbulk -= qgmb;
                        qsrc = -(qgate + qgmid + qbulk + qdrn);
                    }
                    else
                    {
                        gcggb = (this._cggb + cgdo + cgso
                              + Param.BSIM4cgbo) * ag0;
                        gcgdb = (this._cgdb - cgdo) * ag0;
                        gcgsb = (this._cgsb - cgso) * ag0;
                        gcgbb = -(gcggb + gcgdb + gcgsb);

                        gcdgb = (this._cdgb - cgdo) * ag0;
                        gcsgb = -(this._cggb + this._cbgb
                              + this._cdgb + cgso) * ag0;
                        gcbgb = (this._cbgb - Param.BSIM4cgbo) * ag0;

                        gcdgmb = gcsgmb = gcbgmb = 0.0;

                        qgb = Param.BSIM4cgbo * vgb;
                        qgate += qgdo + qgso + qgb;
                        qbulk -= qgb;
                        qsrc = -(qgate + qbulk + qdrn);
                    }
                    gcddb = (this._cddb + this._capbd + cgdo) * ag0;
                    gcdsb = this._cdsb * ag0;

                    gcsdb = -(this._cgdb + this._cbdb
                          + this._cddb) * ag0;
                    gcssb = (this._capbs + cgso - (this._cgsb
                          + this._cbsb + this._cdsb)) * ag0;

                    if (Parameters.RbodyMod.Value == 0)
                    {
                        gcdbb = -(gcdgb + gcddb + gcdsb + gcdgmb);
                        gcsbb = -(gcsgb + gcsdb + gcssb + gcsgmb);
                        gcbdb = (this._cbdb - this._capbd) * ag0;
                        gcbsb = (this._cbsb - this._capbs) * ag0;
                        gcdbdb = 0.0; gcsbsb = 0.0;
                    }
                    else
                    {
                        gcdbb = -(this._cddb + this._cdgb
                               + this._cdsb) * ag0;
                        gcsbb = -(gcsgb + gcsdb + gcssb + gcsgmb)
                              + this._capbs * ag0;
                        gcbdb = this._cbdb * ag0;
                        gcbsb = this._cbsb * ag0;

                        gcdbdb = -this._capbd * ag0;
                        gcsbsb = -this._capbs * ag0;
                    }
                    gcbbb = -(gcbdb + gcbgb + gcbsb + gcbgmb);

                    ggtg = ggtd = ggtb = ggts = 0.0;
                    sxpart = 0.6;
                    dxpart = 0.4;
                    ddxpart_dVd = ddxpart_dVg = ddxpart_dVb = ddxpart_dVs = 0.0;
                    dsxpart_dVd = dsxpart_dVg = dsxpart_dVb = dsxpart_dVs = 0.0;
                }
                else
                {
                    qcheq = this._qchqs;
                    CoxWL = ModelTemperature.Coxe * Param.BSIM4weffCV * Parameters.Nf
                          * Param.BSIM4leffCV;
                    T0 = qdef * ScalingFactor / CoxWL;

                    ggtg = this._gtg = T0 * this._gcrgg;
                    ggtd = this._gtd = T0 * this._gcrgd;
                    ggts = this._gts = T0 * this._gcrgs;
                    ggtb = this._gtb = T0 * this._gcrgb;
                    gqdef = ScalingFactor * ag0;

                    gcqgb = this._cqgb * ag0;
                    gcqdb = this._cqdb * ag0;
                    gcqsb = this._cqsb * ag0;
                    gcqbb = this._cqbb * ag0;

                    if (Math.Abs(qcheq) <= 1.0e-5 * CoxWL)
                    {
                        if (ModelParameters.Xpart < 0.5)
                        {
                            dxpart = 0.4;
                        }
                        else if (ModelParameters.Xpart > 0.5)
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

                    if (Parameters.RgateMod == 3)
                    {
                        gcgmgmb = (cgdo + cgso + Param.BSIM4cgbo) * ag0;
                        gcgmdb = -cgdo * ag0;
                        gcgmsb = -cgso * ag0;
                        gcgmbb = -Param.BSIM4cgbo * ag0;

                        gcdgmb = gcgmdb;
                        gcsgmb = gcgmsb;
                        gcbgmb = gcgmbb;

                        gcdgb = gcsgb = gcbgb = 0.0;
                        gcggb = gcgdb = gcgsb = gcgbb = 0.0;

                        qgmb = Param.BSIM4cgbo * vgmb;
                        qgmid = qgdo + qgso + qgmb;
                        qgate = 0.0;
                        qbulk = -qgmb;
                        qdrn = -qgdo;
                        qsrc = -(qgmid + qbulk + qdrn);
                    }
                    else
                    {
                        gcggb = (cgdo + cgso + Param.BSIM4cgbo) * ag0;
                        gcgdb = -cgdo * ag0;
                        gcgsb = -cgso * ag0;
                        gcgbb = -Param.BSIM4cgbo * ag0;

                        gcdgb = gcgdb;
                        gcsgb = gcgsb;
                        gcbgb = gcgbb;
                        gcdgmb = gcsgmb = gcbgmb = 0.0;

                        qgb = Param.BSIM4cgbo * vgb;
                        qgate = qgdo + qgso + qgb;
                        qbulk = -qgb;
                        qdrn = -qgdo;
                        qsrc = -(qgate + qbulk + qdrn);
                    }

                    gcddb = (this._capbd + cgdo) * ag0;
                    gcdsb = gcsdb = 0.0;
                    gcssb = (this._capbs + cgso) * ag0;

                    if (Parameters.RbodyMod.Value == 0)
                    {
                        gcdbb = -(gcdgb + gcddb + gcdgmb);
                        gcsbb = -(gcsgb + gcssb + gcsgmb);
                        gcbdb = -this._capbd * ag0;
                        gcbsb = -this._capbs * ag0;
                        gcdbdb = 0.0; gcsbsb = 0.0;
                    }
                    else
                    {
                        gcdbb = gcsbb = gcbdb = gcbsb = 0.0;
                        gcdbdb = -this._capbd * ag0;
                        gcsbsb = -this._capbs * ag0;
                    }
                    gcbbb = -(gcbdb + gcbgb + gcbsb + gcbgmb);
                }
            }
            else
            {
                if (Parameters.TrnqsMod == 0)
                {
                    qsrc = qdrn - qgso;
                    if (Parameters.RgateMod == 3)
                    {
                        gcgmgmb = (cgdo + cgso + Param.BSIM4cgbo) * ag0;
                        gcgmdb = -cgdo * ag0;
                        gcgmsb = -cgso * ag0;
                        gcgmbb = -Param.BSIM4cgbo * ag0;

                        gcdgmb = gcgmdb;
                        gcsgmb = gcgmsb;
                        gcbgmb = gcgmbb;

                        gcggb = this._cggb * ag0;
                        gcgdb = this._cgsb * ag0;
                        gcgsb = this._cgdb * ag0;
                        gcgbb = -(gcggb + gcgdb + gcgsb);

                        gcdgb = -(this._cggb + this._cbgb
                              + this._cdgb) * ag0;
                        gcsgb = this._cdgb * ag0;
                        gcbgb = this._cbgb * ag0;

                        qgmb = Param.BSIM4cgbo * vgmb;
                        qgmid = qgdo + qgso + qgmb;
                        qbulk -= qgmb;
                        qdrn = -(qgate + qgmid + qbulk + qsrc);
                    }
                    else
                    {
                        gcggb = (this._cggb + cgdo + cgso
                              + Param.BSIM4cgbo) * ag0;
                        gcgdb = (this._cgsb - cgdo) * ag0;
                        gcgsb = (this._cgdb - cgso) * ag0;
                        gcgbb = -(gcggb + gcgdb + gcgsb);

                        gcdgb = -(this._cggb + this._cbgb
                              + this._cdgb + cgdo) * ag0;
                        gcsgb = (this._cdgb - cgso) * ag0;
                        gcbgb = (this._cbgb - Param.BSIM4cgbo) * ag0;

                        gcdgmb = gcsgmb = gcbgmb = 0.0;

                        qgb = Param.BSIM4cgbo * vgb;
                        qgate += qgdo + qgso + qgb;
                        qbulk -= qgb;
                        qdrn = -(qgate + qbulk + qsrc);
                    }
                    gcddb = (this._capbd + cgdo - (this._cgsb
                          + this._cbsb + this._cdsb)) * ag0;
                    gcdsb = -(this._cgdb + this._cbdb
                          + this._cddb) * ag0;

                    gcsdb = this._cdsb * ag0;
                    gcssb = (this._cddb + this._capbs + cgso) * ag0;

                    if (Parameters.RbodyMod.Value == 0)
                    {
                        gcdbb = -(gcdgb + gcddb + gcdsb + gcdgmb);
                        gcsbb = -(gcsgb + gcsdb + gcssb + gcsgmb);
                        gcbdb = (this._cbsb - this._capbd) * ag0;
                        gcbsb = (this._cbdb - this._capbs) * ag0;
                        gcdbdb = 0.0; gcsbsb = 0.0;
                    }
                    else
                    {
                        gcdbb = -(gcdgb + gcddb + gcdsb + gcdgmb)
                              + this._capbd * ag0;
                        gcsbb = -(this._cddb + this._cdgb
                              + this._cdsb) * ag0;
                        gcbdb = this._cbsb * ag0;
                        gcbsb = this._cbdb * ag0;
                        gcdbdb = -this._capbd * ag0;
                        gcsbsb = -this._capbs * ag0;
                    }
                    gcbbb = -(gcbgb + gcbdb + gcbsb + gcbgmb);

                    ggtg = ggtd = ggtb = ggts = 0.0;
                    sxpart = 0.4;
                    dxpart = 0.6;
                    ddxpart_dVd = ddxpart_dVg = ddxpart_dVb = ddxpart_dVs = 0.0;
                    dsxpart_dVd = dsxpart_dVg = dsxpart_dVb = dsxpart_dVs = 0.0;
                }
                else
                {
                    qcheq = this._qchqs;
                    CoxWL = ModelTemperature.Coxe * Param.BSIM4weffCV * Parameters.Nf
                          * Param.BSIM4leffCV;
                    T0 = qdef * ScalingFactor / CoxWL;
                    ggtg = this._gtg = T0 * this._gcrgg;
                    ggts = this._gts = T0 * this._gcrgd;
                    ggtd = this._gtd = T0 * this._gcrgs;
                    ggtb = this._gtb = T0 * this._gcrgb;
                    gqdef = ScalingFactor * ag0;

                    gcqgb = this._cqgb * ag0;
                    gcqdb = this._cqsb * ag0;
                    gcqsb = this._cqdb * ag0;
                    gcqbb = this._cqbb * ag0;

                    if (Math.Abs(qcheq) <= 1.0e-5 * CoxWL)
                    {
                        if (ModelParameters.Xpart < 0.5)
                        {
                            sxpart = 0.4;
                        }
                        else if (ModelParameters.Xpart > 0.5)
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

                    if (Parameters.RgateMod == 3)
                    {
                        gcgmgmb = (cgdo + cgso + Param.BSIM4cgbo) * ag0;
                        gcgmdb = -cgdo * ag0;
                        gcgmsb = -cgso * ag0;
                        gcgmbb = -Param.BSIM4cgbo * ag0;

                        gcdgmb = gcgmdb;
                        gcsgmb = gcgmsb;
                        gcbgmb = gcgmbb;

                        gcdgb = gcsgb = gcbgb = 0.0;
                        gcggb = gcgdb = gcgsb = gcgbb = 0.0;

                        qgmb = Param.BSIM4cgbo * vgmb;
                        qgmid = qgdo + qgso + qgmb;
                        qgate = 0.0;
                        qbulk = -qgmb;
                        qdrn = -qgdo;
                        qsrc = -qgso;
                    }
                    else
                    {
                        gcggb = (cgdo + cgso + Param.BSIM4cgbo) * ag0;
                        gcgdb = -cgdo * ag0;
                        gcgsb = -cgso * ag0;
                        gcgbb = -Param.BSIM4cgbo * ag0;

                        gcdgb = gcgdb;
                        gcsgb = gcgsb;
                        gcbgb = gcgbb;
                        gcdgmb = gcsgmb = gcbgmb = 0.0;

                        qgb = Param.BSIM4cgbo * vgb;
                        qgate = qgdo + qgso + qgb;
                        qbulk = -qgb;
                        qdrn = -qgdo;
                        qsrc = -qgso;
                    }

                    gcddb = (this._capbd + cgdo) * ag0;
                    gcdsb = gcsdb = 0.0;
                    gcssb = (this._capbs + cgso) * ag0;
                    if (Parameters.RbodyMod.Value == 0)
                    {
                        gcdbb = -(gcdgb + gcddb + gcdgmb);
                        gcsbb = -(gcsgb + gcssb + gcsgmb);
                        gcbdb = -this._capbd * ag0;
                        gcbsb = -this._capbs * ag0;
                        gcdbdb = 0.0; gcsbsb = 0.0;
                    }
                    else
                    {
                        gcdbb = gcsbb = gcbdb = gcbsb = 0.0;
                        gcdbdb = -this._capbd * ag0;
                        gcsbsb = -this._capbs * ag0;
                    }
                    gcbbb = -(gcbdb + gcbgb + gcbsb + gcbgmb);
                }
            }


            if (Parameters.TrnqsMod.Value != 0)
            {
                this._qcdump.Value = qdef * ScalingFactor;
                _qcdump.Derive();
            }

            this._qg.Value = qgate;
            this._qd.Value = qdrn - this._qbd.Value;
            this._qs.Value = qsrc - this._qbs.Value;
            if (Parameters.RgateMod.Value == 3)
                this._qgmid.Value = qgmid;

            if (Parameters.RbodyMod.Value == 0)
            {
                this._qb.Value = qbulk + this._qbd.Value + this._qbs.Value;
            }
            else
                this._qb.Value = qbulk;


            /* Store small signal parameters */
            if (InitializeSmallSignal)
                goto line1000;

            if (!ChargeComputationNeeded)
                goto line850;

            _qb.Derive();
            _qg.Derive();
            _qd.Derive();
            if (Parameters.RgateMod == 3)
                _qgmid.Derive();

            if (Parameters.RbodyMod.Value != 0)
            {
                _qbs.Derive();
                _qbd.Derive();
            }

            goto line860;


        line850:
            /* Zero gcap and ceqcap if (!ChargeComputationNeeded) */
            ceqqg = ceqqb = ceqqd = 0.0;
            ceqqjd = ceqqjs = 0.0;
            cqcheq = cqdef = 0.0;

            gcdgb = gcddb = gcdsb = gcdbb = 0.0;
            gcsgb = gcsdb = gcssb = gcsbb = 0.0;
            gcggb = gcgdb = gcgsb = gcgbb = 0.0;
            gcbdb = gcbgb = gcbsb = gcbbb = 0.0;

            gcgmgmb = gcgmdb = gcgmsb = gcgmbb = 0.0;
            gcdgmb = gcsgmb = gcbgmb = ceqqgmid = 0.0;
            gcdbdb = gcsbsb = 0.0;

            gqdef = gcqgb = gcqdb = gcqsb = gcqbb = 0.0;
            ggtg = ggtd = ggtb = ggts = 0.0;
            sxpart = (1.0 - (dxpart = (this._mode > 0) ? 0.4 : 0.6));
            ddxpart_dVd = ddxpart_dVg = ddxpart_dVb = ddxpart_dVs = 0.0;
            dsxpart_dVd = dsxpart_dVg = dsxpart_dVb = dsxpart_dVs = 0.0;

            if (Parameters.TrnqsMod.Value != 0)
            {
                CoxWL = ModelTemperature.Coxe * Param.BSIM4weffCV * Parameters.Nf
                      * Param.BSIM4leffCV;
                T1 = this._gcrg / CoxWL;
                this._gtau = T1 * ScalingFactor;
            }
            else
                this._gtau = 0.0;

            goto line900;


        line860:
            /* Calculate equivalent charge current */

            cqgate = this._qg.Derivative;
            cqbody = this._qb.Derivative;
            cqdrn = this._qd.Derivative;

            ceqqg = cqgate - gcggb * vgb + gcgdb * vbd + gcgsb * vbs;
            ceqqd = cqdrn - gcdgb * vgb - gcdgmb * vgmb + (gcddb + gcdbdb)
                  * vbd - gcdbdb * vbd_jct + gcdsb * vbs;
            ceqqb = cqbody - gcbgb * vgb - gcbgmb * vgmb
                  + gcbdb * vbd + gcbsb * vbs;


            if (Parameters.RgateMod == 3)
                ceqqgmid = this._qgmid.Derivative
                         + gcgmdb * vbd + gcgmsb * vbs - gcgmgmb * vgmb;
            else
                ceqqgmid = 0.0;

            if (Parameters.RbodyMod.Value != 0)
            {
                ceqqjs = this._qbs.Derivative + gcsbsb * vbs_jct;
                ceqqjd = this._qbd.Derivative + gcdbdb * vbd_jct;
            }

            if (Parameters.TrnqsMod.Value != 0)
            {
                T0 = ggtg * vgb - ggtd * vbd - ggts * vbs;
                ceqqg += T0;
                T1 = qdef * this._gtau;
                ceqqd -= dxpart * T0 + T1 * (ddxpart_dVg * vgb - ddxpart_dVd
                       * vbd - ddxpart_dVs * vbs);
                cqdef = this._qcdump.Derivative - gqdef * qdef;
                cqcheq = this._qcheq.Derivative
                       - (gcqgb * vgb - gcqdb * vbd - gcqsb * vbs) + T0;
            }

        /*
         *  Load current vector
         */

        line900:
            if (this._mode >= 0)
            {
                Gm = this._gm;
                Gmbs = this._gmbs;
                FwdSum = Gm + Gmbs;
                RevSum = 0.0;

                ceqdrn = ModelParameters.Type * (cdrain - this._gds * vds
                       - Gm * vgs - Gmbs * vbs);
                ceqbd = ModelParameters.Type * (this._csub + this._igidl
                      - (this._gbds + this._ggidld) * vds
                      - (this._gbgs + this._ggidlg) * vgs
                      - (this._gbbs + this._ggidlb) * vbs);
                ceqbs = ModelParameters.Type * (this._igisl + this._ggisls * vds
                            - this._ggislg * vgd - this._ggislb * vbd);

                gbbdp = -(this._gbds);
                gbbsp = this._gbds + this._gbgs + this._gbbs;

                gbdpg = this._gbgs;
                gbdpdp = this._gbds;
                gbdpb = this._gbbs;
                gbdpsp = -(gbdpg + gbdpdp + gbdpb);

                gbspg = 0.0;
                gbspdp = 0.0;
                gbspb = 0.0;
                gbspsp = 0.0;

                if (ModelParameters.IgcMod.Value != 0)
                {
                    gIstotg = this._gIgsg + this._gIgcsg;
                    gIstotd = this._gIgcsd;
                    gIstots = this._gIgss + this._gIgcss;
                    gIstotb = this._gIgcsb;
                    Istoteq = ModelParameters.Type * (this._igs + this._igcs
                             - gIstotg * vgs - this._gIgcsd * vds
                           - this._gIgcsb * vbs);

                    gIdtotg = this._gIgdg + this._gIgcdg;
                    gIdtotd = this._gIgdd + this._gIgcdd;
                    gIdtots = this._gIgcds;
                    gIdtotb = this._gIgcdb;
                    Idtoteq = ModelParameters.Type * (this._igd + this._igcd
                            - this._gIgdg * vgd - this._gIgcdg * vgs
                            - this._gIgcdd * vds - this._gIgcdb * vbs);
                }
                else
                {
                    gIstotg = gIstotd = gIstots = gIstotb = Istoteq = 0.0;
                    gIdtotg = gIdtotd = gIdtots = gIdtotb = Idtoteq = 0.0;
                }

                if (ModelParameters.IgbMod.Value != 0)
                {
                    gIbtotg = this._gIgbg;
                    gIbtotd = this._gIgbd;
                    gIbtots = this._gIgbs;
                    gIbtotb = this._gIgbb;
                    Ibtoteq = ModelParameters.Type * (this._igb
                            - this._gIgbg * vgs - this._gIgbd * vds
                            - this._gIgbb * vbs);
                }
                else
                    gIbtotg = gIbtotd = gIbtots = gIbtotb = Ibtoteq = 0.0;

                if ((ModelParameters.IgcMod != 0) || (ModelParameters.IgbMod != 0))
                {
                    gIgtotg = gIstotg + gIdtotg + gIbtotg;
                    gIgtotd = gIstotd + gIdtotd + gIbtotd;
                    gIgtots = gIstots + gIdtots + gIbtots;
                    gIgtotb = gIstotb + gIdtotb + gIbtotb;
                    Igtoteq = Istoteq + Idtoteq + Ibtoteq;
                }
                else
                    gIgtotg = gIgtotd = gIgtots = gIgtotb = Igtoteq = 0.0;


                if (Parameters.RgateMod == 2)
                    T0 = vges - vgs;
                else if (Parameters.RgateMod == 3)
                    T0 = vgms - vgs;
                if (Parameters.RgateMod > 1)
                {
                    gcrgd = this._gcrgd * T0;
                    gcrgg = this._gcrgg * T0;
                    gcrgs = this._gcrgs * T0;
                    gcrgb = this._gcrgb * T0;
                    ceqgcrg = -(gcrgd * vds + gcrgg * vgs
                            + gcrgb * vbs);
                    gcrgg -= this._gcrg;
                    gcrg = this._gcrg;
                }
                else
                    ceqgcrg = gcrg = gcrgd = gcrgg = gcrgs = gcrgb = 0.0;
            }
            else
            {
                Gm = -this._gm;
                Gmbs = -this._gmbs;
                FwdSum = 0.0;
                RevSum = -(Gm + Gmbs);

                ceqdrn = -ModelParameters.Type * (cdrain + this._gds * vds
                       + Gm * vgd + Gmbs * vbd);

                ceqbs = ModelParameters.Type * (this._csub + this._igisl
                      + (this._gbds + this._ggisls) * vds
                      - (this._gbgs + this._ggislg) * vgd
                      - (this._gbbs + this._ggislb) * vbd);
                ceqbd = ModelParameters.Type * (this._igidl - this._ggidld * vds
                                - this._ggidlg * vgs - this._ggidlb * vbs);

                gbbsp = -(this._gbds);
                gbbdp = this._gbds + this._gbgs + this._gbbs;

                gbdpg = 0.0;
                gbdpsp = 0.0;
                gbdpb = 0.0;
                gbdpdp = 0.0;

                gbspg = this._gbgs;
                gbspsp = this._gbds;
                gbspb = this._gbbs;
                gbspdp = -(gbspg + gbspsp + gbspb);

                if (ModelParameters.IgcMod.Value != 0)
                {
                    gIstotg = this._gIgsg + this._gIgcdg;
                    gIstotd = this._gIgcds;
                    gIstots = this._gIgss + this._gIgcdd;
                    gIstotb = this._gIgcdb;
                    Istoteq = ModelParameters.Type * (this._igs + this._igcd
                            - this._gIgsg * vgs - this._gIgcdg * vgd
                            + this._gIgcdd * vds - this._gIgcdb * vbd);

                    gIdtotg = this._gIgdg + this._gIgcsg;
                    gIdtotd = this._gIgdd + this._gIgcss;
                    gIdtots = this._gIgcsd;
                    gIdtotb = this._gIgcsb;
                    Idtoteq = ModelParameters.Type * (this._igd + this._igcs
                            - (this._gIgdg + this._gIgcsg) * vgd
                            + this._gIgcsd * vds - this._gIgcsb * vbd);
                }
                else
                {
                    gIstotg = gIstotd = gIstots = gIstotb = Istoteq = 0.0;
                    gIdtotg = gIdtotd = gIdtots = gIdtotb = Idtoteq = 0.0;
                }

                if (ModelParameters.IgbMod.Value != 0)
                {
                    gIbtotg = this._gIgbg;
                    gIbtotd = this._gIgbs;
                    gIbtots = this._gIgbd;
                    gIbtotb = this._gIgbb;
                    Ibtoteq = ModelParameters.Type * (this._igb
                            - this._gIgbg * vgd + this._gIgbd * vds
                            - this._gIgbb * vbd);
                }
                else
                    gIbtotg = gIbtotd = gIbtots = gIbtotb = Ibtoteq = 0.0;

                if ((ModelParameters.IgcMod != 0) || (ModelParameters.IgbMod != 0))
                {
                    gIgtotg = gIstotg + gIdtotg + gIbtotg;
                    gIgtotd = gIstotd + gIdtotd + gIbtotd;
                    gIgtots = gIstots + gIdtots + gIbtots;
                    gIgtotb = gIstotb + gIdtotb + gIbtotb;
                    Igtoteq = Istoteq + Idtoteq + Ibtoteq;
                }
                else
                    gIgtotg = gIgtotd = gIgtots = gIgtotb = Igtoteq = 0.0;


                if (Parameters.RgateMod == 2)
                    T0 = vges - vgs;
                else if (Parameters.RgateMod == 3)
                    T0 = vgms - vgs;
                if (Parameters.RgateMod > 1)
                {
                    gcrgd = this._gcrgs * T0;
                    gcrgg = this._gcrgg * T0;
                    gcrgs = this._gcrgd * T0;
                    gcrgb = this._gcrgb * T0;
                    ceqgcrg = -(gcrgg * vgd - gcrgs * vds
                            + gcrgb * vbd);
                    gcrgg -= this._gcrg;
                    gcrg = this._gcrg;
                }
                else
                    ceqgcrg = gcrg = gcrgd = gcrgg = gcrgs = gcrgb = 0.0;
            }

            if (ModelParameters.RdsMod == 1)
            {
                ceqgstot = ModelParameters.Type * (this._gstotd * vds
                         + this._gstotg * vgs + this._gstotb * vbs);
                /* WDLiu: ceqgstot flowing away from sNodePrime */
                gstot = this._gstot;
                gstotd = this._gstotd;
                gstotg = this._gstotg;
                gstots = this._gstots - gstot;
                gstotb = this._gstotb;

                ceqgdtot = -ModelParameters.Type * (this._gdtotd * vds
                         + this._gdtotg * vgs + this._gdtotb * vbs);
                /* WDLiu: ceqgdtot defined as flowing into dNodePrime */
                gdtot = this._gdtot;
                gdtotd = this._gdtotd - gdtot;
                gdtotg = this._gdtotg;
                gdtots = this._gdtots;
                gdtotb = this._gdtotb;
            }
            else
            {
                gstot = gstotd = gstotg = gstots = gstotb = ceqgstot = 0.0;
                gdtot = gdtotd = gdtotg = gdtots = gdtotb = ceqgdtot = 0.0;
            }

            if (ModelParameters.Type > 0)
            {
                ceqjs = (this._cbs - this._gbs * vbs_jct);
                ceqjd = (this._cbd - this._gbd * vbd_jct);
            }
            else
            {
                ceqjs = -(this._cbs - this._gbs * vbs_jct);
                ceqjd = -(this._cbd - this._gbd * vbd_jct);
                ceqqg = -ceqqg;
                ceqqd = -ceqqd;
                ceqqb = -ceqqb;
                ceqgcrg = -ceqgcrg;

                if (Parameters.TrnqsMod.Value != 0)
                {
                    cqdef = -cqdef;
                    cqcheq = -cqcheq;
                }

                if (Parameters.RbodyMod.Value != 0)
                {
                    ceqqjs = -ceqqjs;
                    ceqqjd = -ceqqjd;
                }

                if (Parameters.RgateMod == 3)
                    ceqqgmid = -ceqqgmid;
            }


            /*
             *  Loading RHS
             */

            m = Parameters.M;


            this._dnodeprimePtr.Value += m * (ceqjd - ceqbd + ceqgdtot
                                                        - ceqdrn - ceqqd + Idtoteq);
            this._gnodeprimePtr.Value -= m * (ceqqg - ceqgcrg + Igtoteq);

            if (Parameters.RgateMod == 2)
                this._gnodeextPtr.Value -= m * ceqgcrg;
            else if (Parameters.RgateMod == 3)
                this._gnodemidPtr.Value -= m * (ceqqgmid + ceqgcrg);

            if (Parameters.RbodyMod.Value == 0)
            {
                this._bnodeprimePtr.Value += m * (ceqbd + ceqbs - ceqjd
                                                            - ceqjs - ceqqb + Ibtoteq);
                this._snodeprimePtr.Value += m * (ceqdrn - ceqbs + ceqjs
                                  + ceqqg + ceqqb + ceqqd + ceqqgmid - ceqgstot + Istoteq);
            }
            else
            {
                this._dbnodePtr.Value -= m * (ceqjd + ceqqjd);
                this._bnodeprimePtr.Value += m * (ceqbd + ceqbs - ceqqb + Ibtoteq);
                this._sbnodePtr.Value -= m * (ceqjs + ceqqjs);
                this._snodeprimePtr.Value += m * (ceqdrn - ceqbs + ceqjs + ceqqd
                    + ceqqg + ceqqb + ceqqjd + ceqqjs + ceqqgmid - ceqgstot + Istoteq);
            }

            if (ModelParameters.RdsMod.Value != 0)
            {
                this._dnodePtr.Value -= m * ceqgdtot;
                this._snodePtr.Value += m * ceqgstot;
            }

            if (Parameters.TrnqsMod.Value != 0)
                _qPtr.Value += m * (cqcheq - cqdef);

            /*
             *  Loading matrix
             */

            if (Parameters.RbodyMod.Value == 0)
            {
                gjbd = this._gbd;
                gjbs = this._gbs;
            }
            else
                gjbd = gjbs = 0.0;

            if (ModelParameters.RdsMod.Value == 0)
            {
                gdpr = this._drainConductance;
                gspr = this._sourceConductance;
            }
            else
                gdpr = gspr = 0.0;

            geltd = this._grgeltd;

            T1 = qdef * this._gtau;

            if (Parameters.RgateMod == 1)
            {
                this._gegePtr.Value += m * geltd;
                this._gpgePtr.Value -= m * geltd;
                this._gegpPtr.Value -= m * geltd;
                this._gpgpPtr.Value += m * (gcggb + geltd - ggtg + gIgtotg);
                this._gpdpPtr.Value += m * (gcgdb - ggtd + gIgtotd);
                this._gpspPtr.Value += m * (gcgsb - ggts + gIgtots);
                this._gpbpPtr.Value += m * (gcgbb - ggtb + gIgtotb);
            } /* WDLiu: gcrg already subtracted from all gcrgg below */
            else if (Parameters.RgateMod == 2)
            {
                this._gegePtr.Value += m * gcrg;
                this._gegpPtr.Value += m * gcrgg;
                this._gedpPtr.Value += m * gcrgd;
                this._gespPtr.Value += m * gcrgs;
                this._gebpPtr.Value += m * gcrgb;

                this._gpgePtr.Value -= m * gcrg;
                this._gpgpPtr.Value += m * (gcggb - gcrgg - ggtg + gIgtotg);
                this._gpdpPtr.Value += m * (gcgdb - gcrgd - ggtd + gIgtotd);
                this._gpspPtr.Value += m * (gcgsb - gcrgs - ggts + gIgtots);
                this._gpbpPtr.Value += m * (gcgbb - gcrgb - ggtb + gIgtotb);
            }
            else if (Parameters.RgateMod == 3)
            {
                this._gegePtr.Value += m * geltd;
                this._gegmPtr.Value -= m * geltd;
                this._gmgePtr.Value -= m * geltd;
                this._gmgmPtr.Value += m * (geltd + gcrg + gcgmgmb);

                this._gmdpPtr.Value += m * (gcrgd + gcgmdb);
                this._gmgpPtr.Value += m * gcrgg;
                this._gmspPtr.Value += m * (gcrgs + gcgmsb);
                this._gmbpPtr.Value += m * (gcrgb + gcgmbb);

                this._dpgmPtr.Value += m * gcdgmb;
                this._gpgmPtr.Value -= m * gcrg;
                this._spgmPtr.Value += m * gcsgmb;
                this._bpgmPtr.Value += m * gcbgmb;

                this._gpgpPtr.Value += m * (gcggb - gcrgg - ggtg + gIgtotg);
                this._gpdpPtr.Value += m * (gcgdb - gcrgd - ggtd + gIgtotd);
                this._gpspPtr.Value += m * (gcgsb - gcrgs - ggts + gIgtots);
                this._gpbpPtr.Value += m * (gcgbb - gcrgb - ggtb + gIgtotb);
            }
            else
            {
                this._gpgpPtr.Value += m * (gcggb - ggtg + gIgtotg);
                this._gpdpPtr.Value += m * (gcgdb - ggtd + gIgtotd);
                this._gpspPtr.Value += m * (gcgsb - ggts + gIgtots);
                this._gpbpPtr.Value += m * (gcgbb - ggtb + gIgtotb);
            }

            if (ModelParameters.RdsMod.Value != 0)
            {
                this._dgpPtr.Value += m * gdtotg;
                this._dspPtr.Value += m * gdtots;
                this._dbpPtr.Value += m * gdtotb;
                this._sdpPtr.Value += m * gstotd;
                this._sgpPtr.Value += m * gstotg;
                this._sbpPtr.Value += m * gstotb;
            }

            this._dpdpPtr.Value += m * (gdpr + this._gds + this._gbd + T1 * ddxpart_dVd
                                    - gdtotd + RevSum + gcddb + gbdpdp + dxpart * ggtd - gIdtotd);
            this._dpdPtr.Value -= m * (gdpr + gdtot);
            this._dpgpPtr.Value += m * (Gm + gcdgb - gdtotg + gbdpg - gIdtotg
                                    + dxpart * ggtg + T1 * ddxpart_dVg);
            this._dpspPtr.Value -= m * (this._gds + gdtots - dxpart * ggts + gIdtots
                                    - T1 * ddxpart_dVs + FwdSum - gcdsb - gbdpsp);
            this._dpbpPtr.Value -= m * (gjbd + gdtotb - Gmbs - gcdbb - gbdpb + gIdtotb
                                    - T1 * ddxpart_dVb - dxpart * ggtb);

            this._ddpPtr.Value -= m * (gdpr - gdtotd);
            this._ddPtr.Value += m * (gdpr + gdtot);

            this._spdpPtr.Value -= m * (this._gds + gstotd + RevSum - gcsdb - gbspdp
                                    - T1 * dsxpart_dVd - sxpart * ggtd + gIstotd);
            this._spgpPtr.Value += m * (gcsgb - Gm - gstotg + gbspg + sxpart * ggtg
                                    + T1 * dsxpart_dVg - gIstotg);
            this._spspPtr.Value += m * (gspr + this._gds + this._gbs + T1 * dsxpart_dVs
                                    - gstots + FwdSum + gcssb + gbspsp + sxpart * ggts - gIstots);
            this._spsPtr.Value -= m * (gspr + gstot);
            this._spbpPtr.Value -= m * (gjbs + gstotb + Gmbs - gcsbb - gbspb - sxpart * ggtb
                                    - T1 * dsxpart_dVb + gIstotb);

            this._sspPtr.Value -= m * (gspr - gstots);
            this._ssPtr.Value += m * (gspr + gstot);

            this._bpdpPtr.Value += m * (gcbdb - gjbd + gbbdp - gIbtotd);
            this._bpgpPtr.Value += m * (gcbgb - this._gbgs - gIbtotg);
            this._bpspPtr.Value += m * (gcbsb - gjbs + gbbsp - gIbtots);
            this._bpbpPtr.Value += m * (gjbd + gjbs + gcbbb - this._gbbs
                                    - gIbtotb);

            ggidld = this._ggidld;
            ggidlg = this._ggidlg;
            ggidlb = this._ggidlb;
            ggislg = this._ggislg;
            ggisls = this._ggisls;
            ggislb = this._ggislb;

            /* stamp gidl */
            this._dpdpPtr.Value += m * ggidld;
            this._dpgpPtr.Value += m * ggidlg;
            this._dpspPtr.Value -= m * (ggidlg + ggidld + ggidlb);
            this._dpbpPtr.Value += m * ggidlb;
            this._bpdpPtr.Value -= m * ggidld;
            this._bpgpPtr.Value -= m * ggidlg;
            this._bpspPtr.Value += m * (ggidlg + ggidld + ggidlb);
            this._bpbpPtr.Value -= m * ggidlb;
            /* stamp gisl */
            this._spdpPtr.Value -= m * (ggisls + ggislg + ggislb);
            this._spgpPtr.Value += m * ggislg;
            this._spspPtr.Value += m * ggisls;
            this._spbpPtr.Value += m * ggislb;
            this._bpdpPtr.Value += m * (ggislg + ggisls + ggislb);
            this._bpgpPtr.Value -= m * ggislg;
            this._bpspPtr.Value -= m * ggisls;
            this._bpbpPtr.Value -= m * ggislb;


            if (Parameters.RbodyMod.Value != 0)
            {
                this._dpdbPtr.Value += m * (gcdbdb - this._gbd);
                this._spsbPtr.Value -= m * (this._gbs - gcsbsb);

                this._dbdpPtr.Value += m * (gcdbdb - this._gbd);
                this._dbdbPtr.Value += m * (this._gbd - gcdbdb
                                        + this._grbpd + this._grbdb);
                this._dbbpPtr.Value -= m * this._grbpd;
                this._dbbPtr.Value -= m * this._grbdb;

                this._bpdbPtr.Value -= m * this._grbpd;
                this._bpbPtr.Value -= m * this._grbpb;
                this._bpsbPtr.Value -= m * this._grbps;
                this._bpbpPtr.Value += m * (this._grbpd + this._grbps
                                        + this._grbpb);
                /* WDLiu: (gcbbb - this._gbbs) already added to BPbpPtr */

                this._sbspPtr.Value += m * (gcsbsb - this._gbs);
                this._sbbpPtr.Value -= m * this._grbps;
                this._sbbPtr.Value -= m * this._grbsb;
                this._sbsbPtr.Value += m * (this._gbs - gcsbsb
                                        + this._grbps + this._grbsb);

                this._bdbPtr.Value -= m * this._grbdb;
                this._bbpPtr.Value -= m * this._grbpb;
                this._bsbPtr.Value -= m * this._grbsb;
                this._bbPtr.Value += m * (this._grbsb + this._grbdb
                                      + this._grbpb);
            }

            if (Parameters.TrnqsMod.Value != 0)
            {
                this._qqPtr.Value += m * (gqdef + this._gtau);
                this._qgpPtr.Value += m * (ggtg - gcqgb);
                this._qdpPtr.Value += m * (ggtd - gcqdb);
                this._qspPtr.Value += m * (ggts - gcqsb);
                this._qbpPtr.Value += m * (ggtb - gcqbb);

                this._dpqPtr.Value += m * dxpart * this._gtau;
                this._spqPtr.Value += m * sxpart * this._gtau;
                this._gpqPtr.Value -= m * this._gtau;
            }

        line1000:;
        }

        private static void DEXP(double A, out double B, out double C)
        {
            if (A > EXP_THRESHOLD)
            {
                B = MAX_EXP * (1.0 + (A) - EXP_THRESHOLD);
                C = MAX_EXP;
            }
            else if (A < -EXP_THRESHOLD)
            {
                B = MIN_EXP;
                C = 0;
            }
            else
            {
                B = Math.Exp(A);
                C = B;
            }
        }

        /* function to compute poly depletion effect */
        private int BSIM4polyDepletion(double phi, double ngate, double epsgate, double coxe, double Vgs,
            out double Vgs_eff, out double dVgs_eff_dVg)
        {
            double T1, T2, T3, T4, T5, T6, T7, T8;

            /* Poly Gate Si Depletion Effect */
            if ((ngate > 1.0e18) &&
                (ngate < 1.0e25) && (Vgs > phi) && (epsgate != 0)
               )
            {
                T1 = 1.0e6 * Constants.Charge * epsgate * ngate / (coxe * coxe);
                T8 = Vgs - phi;
                T4 = Math.Sqrt(1.0 + 2.0 * T8 / T1);
                T2 = 2.0 * T8 / (T4 + 1.0);
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
            return (0);
        }

    }
}
