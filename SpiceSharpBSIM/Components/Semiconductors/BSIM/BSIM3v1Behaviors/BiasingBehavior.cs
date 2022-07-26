using SpiceSharp;
using SpiceSharp.Algebra;
using SpiceSharp.Attributes;
using SpiceSharp.Behaviors;
using SpiceSharp.Components;
using SpiceSharp.Components.Mosfets;
using SpiceSharp.Components.Semiconductors;
using SpiceSharp.Simulations;
using System;
using System.Collections.Generic;
using System.Text;

namespace SpiceSharpBSIM.Components.Semiconductors.BSIM.BSIM3v1Behaviors
{
    /// <summary>
    /// Biasing behavior for a <see cref="BSIM3v1"/>.
    /// </summary>
    [BehaviorFor(typeof(BSIM3v1)), AddBehaviorIfNo(typeof(IBiasingBehavior))]
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
            double ag0, qgd, qgs, qgb, von, cbhat, VgstNVt, ExpVgst = 0.0;
            double cdrain, cdhat, cdreq, ceqbd, ceqbs, ceqqb, ceqqd, ceqqg;
            double czbd, czbdsw, czbdswg, czbs, czbssw, czbsswg, evbd, evbs, arg, sarg;
            double delvbd, delvbs, delvds, delvgd, delvgs;
            double Vfbeff, dVfbeff_dVg, dVfbeff_dVd, dVfbeff_dVb, V3, V4;
            double gcbdb, gcbgb, gcbsb, gcddb, gcdgb, gcdsb, gcgdb, gcggb, gcgsb, gcsdb;
            double gcsgb, gcssb, PhiB, PhiBSW, MJ, MJSW, PhiBSWG, MJSWG;
            double vbd, vbs, vds, vgb, vgd, vgs, vgdo;

            double qgate = 0.0, qbulk = 0.0, qdrn = 0.0, qsrc, cqgate, cqbulk, cqdrn;
            double Vds, Vgs, Vbs, Gmbs, FwdSum, RevSum;
            double Vgs_eff, Vfb, dVfb_dVb, dVfb_dVd;
            double Phis, dPhis_dVb, sqrtPhis, dsqrtPhis_dVb, Vth, dVth_dVb, dVth_dVd;
            double Vgst, dVgst_dVg, dVgst_dVb, dVgs_eff_dVg, Nvtm;
            double n, dn_dVb, Vtm;
            double ExpArg, V0;
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
            double Leff = 0.0, Weff, dWeff_dVg, dWeff_dVb;
            double AbulkCV, dAbulkCV_dVb;
            double qgdo, qgso, cgdo, cgso;

            double qcheq, qdef, gqdef, cqdef, cqcheq, gtau_diff, gtau_drift;
            double gcqdb, gcqsb, gcqgb, gcqbb;
            double dxpart, sxpart;

            double gbspsp, gbbdp, gbbsp, gbspg, gbspb, gbspdp;
            double gbdpdp, gbdpg, gbdpb, gbdpsp;
            double Cgg, Cgd, Cgb;
            double Csg, Csd, Csb, Cbg, Cbd, Cbb;
            double Cgg1, Cgb1, Cgd1, Cbg1, Cbb1, Cbd1, Qac0, Qsub0;
            double dQac0_dVg, dQac0_dVd, dQac0_dVb, dQsub0_dVg, dQsub0_dVd, dQsub0_dVb;

            double m;

            bool Check, ChargeComputationNeeded;

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
                    vgs = Param.BSIM3v1vth0 + 0.1;
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
                delvbs = vbs - this._vbs;
                delvbd = vbd - this._vbd;
                delvgs = vgs - this._vgs;
                delvds = vds - this._vds;
                delvgd = vgd - vgdo;

                if (this._mode >= 0)
                {
                    cdhat = this._cd - this._gbd * delvbd
                              + this._gmbs * delvbs + this._gm * delvgs
                              + this._gds * delvds;
                }
                else
                {
                    cdhat = this._cd - (this._gbd - this._gmbs)
                              * delvbd - this._gm * delvgd
                              + this._gds * delvds;

                }
                cbhat = this._cbs + this._cbd + this._gbd
                      * delvbd + this._gbs * delvbs;

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

            /* the following code has been changed for the calculation of S/B and D/B diodes*/

            if ((Parameters.SourceArea <= 0.0) &&
                (Parameters.SourcePerimeter <= 0.0))
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

            Nvtm = ModelTemperature.Vtm * ModelParameters.JctEmissionCoeff;
            if (SourceSatCurrent <= 0.0)
            {
                this._gbs = _iteration.Gmin;
                this._cbs = this._gbs * vbs;
            }
            else if (vbs < 0.5)
            {
                evbs = Math.Exp(vbs / Nvtm);
                this._gbs = SourceSatCurrent * evbs / Nvtm + _iteration.Gmin;
                this._cbs = SourceSatCurrent * (evbs - 1.0)
                               + _iteration.Gmin * vbs;
            }
            else
            {
                evbs = Math.Exp(0.5 / Nvtm);
                T0 = SourceSatCurrent * evbs / Nvtm;
                this._gbs = T0 + _iteration.Gmin;
                this._cbs = SourceSatCurrent * (evbs - 1.0)
                   + T0 * (vbs - 0.5) + _iteration.Gmin * vbs;
            }

            if ((Parameters.DrainArea <= 0.0) &&
                (Parameters.DrainPerimeter <= 0.0))
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

            if (DrainSatCurrent <= 0.0)
            {
                this._gbd = _iteration.Gmin;
                this._cbd = this._gbd * vbd;
            }
            else if (vbd < 0.5)
            {
                evbd = Math.Exp(vbd / Nvtm);
                this._gbd = DrainSatCurrent * evbd / Nvtm + _iteration.Gmin;
                this._cbd = DrainSatCurrent * (evbd - 1.0)
                               + _iteration.Gmin * vbd;
            }
            else
            {
                evbd = Math.Exp(0.5 / Nvtm);
                T0 = DrainSatCurrent * evbd / Nvtm;
                this._gbd = T0 + _iteration.Gmin;
                this._cbd = DrainSatCurrent * (evbd - 1.0)
                   + T0 * (vbd - 0.5) + _iteration.Gmin * vbd;
            }
            /* S/B and D/B diodes code change ends */

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
            T0 = Vbs - Param.BSIM3v1vbsc - 0.001;
            T1 = Math.Sqrt(T0 * T0 - 0.004 * Param.BSIM3v1vbsc);
            Vbseff = Param.BSIM3v1vbsc + 0.5 * (T0 + T1);
            dVbseff_dVb = 0.5 * (1.0 + T0 / T1);
            if (Vbseff < Vbs)
            {
                Vbseff = Vbs;
            } /* Added to avoid the possible numerical problems due to computer accuracy. See comments for diffVds */

            if (Vbseff > 0.0)
            {
                T0 = Param.BSIM3v1phi / (Param.BSIM3v1phi + Vbseff);
                Phis = Param.BSIM3v1phi * T0;
                dPhis_dVb = -T0 * T0;
                sqrtPhis = Param.BSIM3v1phis3 / (Param.BSIM3v1phi + 0.5 * Vbseff);
                dsqrtPhis_dVb = -0.5 * sqrtPhis * sqrtPhis / Param.BSIM3v1phis3;
            }
            else
            {
                Phis = Param.BSIM3v1phi - Vbseff;
                dPhis_dVb = -1.0;
                sqrtPhis = Math.Sqrt(Phis);
                dsqrtPhis_dVb = -0.5 / sqrtPhis;
            }
            Xdep = Param.BSIM3v1Xdep0 * sqrtPhis / Param.BSIM3v1sqrtPhi;
            dXdep_dVb = (Param.BSIM3v1Xdep0 / Param.BSIM3v1sqrtPhi)
              * dsqrtPhis_dVb;

            Leff = Param.BSIM3v1leff;
            Vtm = ModelTemperature.Vtm;
            /* Vth Calculation */
            T3 = Math.Sqrt(Xdep);
            V0 = Param.BSIM3v1vbi - Param.BSIM3v1phi;

            T0 = Param.BSIM3v1dvt2 * Vbseff;
            if (T0 >= -0.5)
            {
                T1 = 1.0 + T0;
                T2 = Param.BSIM3v1dvt2;
            }
            else /* Added to avoid any discontinuity problems caused by dvt2 */
            {
                T4 = 1.0 / (3.0 + 8.0 * T0);
                T1 = (1.0 + 3.0 * T0) * T4;
                T2 = Param.BSIM3v1dvt2 * T4 * T4;
            }
            lt1 = ModelTemperature.Factor1 * T3 * T1;
            dlt1_dVb = ModelTemperature.Factor1 * (0.5 / T3 * T1 * dXdep_dVb + T3 * T2);

            T0 = Param.BSIM3v1dvt2w * Vbseff;
            if (T0 >= -0.5)
            {
                T1 = 1.0 + T0;
                T2 = Param.BSIM3v1dvt2w;
            }
            else /* Added to avoid any discontinuity problems caused by dvt2w */
            {
                T4 = 1.0 / (3.0 + 8.0 * T0);
                T1 = (1.0 + 3.0 * T0) * T4;
                T2 = Param.BSIM3v1dvt2w * T4 * T4;
            }
            ltw = ModelTemperature.Factor1 * T3 * T1;
            dltw_dVb = ModelTemperature.Factor1 * (0.5 / T3 * T1 * dXdep_dVb + T3 * T2);

            T0 = -0.5 * Param.BSIM3v1dvt1 * Leff / lt1;
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

            this._thetavth = Param.BSIM3v1dvt0 * Theta0;
            Delt_vth = this._thetavth * V0;
            dDelt_vth_dVb = Param.BSIM3v1dvt0 * dTheta0_dVb * V0;

            T0 = -0.5 * Param.BSIM3v1dvt1w * Param.BSIM3v1weff * Leff / ltw;
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

            T0 = Param.BSIM3v1dvt0w * T2;
            T2 = T0 * V0;
            dT2_dVb = Param.BSIM3v1dvt0w * dT2_dVb * V0;

            TempRatio = _temperature.Temperature / ModelParameters.Tnom - 1.0;
            T0 = Math.Sqrt(1.0 + Param.BSIM3v1nlx / Leff);
            T1 = Param.BSIM3v1k1 * (T0 - 1.0) * Param.BSIM3v1sqrtPhi
               + (Param.BSIM3v1kt1 + Param.BSIM3v1kt1l / Leff
               + Param.BSIM3v1kt2 * Vbseff) * TempRatio;
            tmp2 = ModelParameters.Tox * Param.BSIM3v1phi
             / (Param.BSIM3v1weff + Param.BSIM3v1w0);

            T3 = Param.BSIM3v1eta0 + Param.BSIM3v1etab * Vbseff;
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
            dDIBL_Sft_dVd = T3 * Param.BSIM3v1theta0vb0;
            DIBL_Sft = dDIBL_Sft_dVd * Vds;

            Vth = ModelParameters.Type * Param.BSIM3v1vth0 + Param.BSIM3v1k1
                * (sqrtPhis - Param.BSIM3v1sqrtPhi) - Param.BSIM3v1k2
                * Vbseff - Delt_vth - T2 + (Param.BSIM3v1k3 + Param.BSIM3v1k3b
                * Vbseff) * tmp2 + T1 - DIBL_Sft;

            this._von = Vth;

            dVth_dVb = Param.BSIM3v1k1 * dsqrtPhis_dVb - Param.BSIM3v1k2
                     - dDelt_vth_dVb - dT2_dVb + Param.BSIM3v1k3b * tmp2
                     - Param.BSIM3v1etab * Vds * Param.BSIM3v1theta0vb0 * T4
                     + Param.BSIM3v1kt2 * TempRatio;
            dVth_dVd = -dDIBL_Sft_dVd;

            /* Calculate n */
            tmp2 = Param.BSIM3v1nfactor * EPSSI / Xdep;
            tmp3 = Param.BSIM3v1cdsc + Param.BSIM3v1cdscb * Vbseff
                 + Param.BSIM3v1cdscd * Vds;
            tmp4 = (tmp2 + tmp3 * Theta0 + Param.BSIM3v1cit) / ModelTemperature.Cox;
            if (tmp4 >= -0.5)
            {
                n = 1.0 + tmp4;
                dn_dVb = (-tmp2 / Xdep * dXdep_dVb + tmp3 * dTheta0_dVb
                           + Param.BSIM3v1cdscb * Theta0) / ModelTemperature.Cox;
                dn_dVd = Param.BSIM3v1cdscd * Theta0 / ModelTemperature.Cox;
            }
            else
            /* avoid  discontinuity problems caused by tmp4 */
            {
                T0 = 1.0 / (3.0 + 8.0 * tmp4);
                n = (1.0 + 3.0 * tmp4) * T0;
                T0 *= T0;
                dn_dVb = (-tmp2 / Xdep * dXdep_dVb + tmp3 * dTheta0_dVb
                           + Param.BSIM3v1cdscb * Theta0) / ModelTemperature.Cox * T0;
                dn_dVd = Param.BSIM3v1cdscd * Theta0 / ModelTemperature.Cox * T0;
            }

            /* Poly Gate Si Depletion Effect */
            T0 = Param.BSIM3v1vfb + Param.BSIM3v1phi;
            if ((Param.BSIM3v1ngate > 1.0e18) && (Param.BSIM3v1ngate < 1.0e25)
                 && (Vgs > T0))
            /* added to avoid the problem caused by ngate */
            {
                T1 = 1.0e6 * Charge_q * EPSSI * Param.BSIM3v1ngate
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
            ExpArg = (2.0 * Param.BSIM3v1voff - Vgst) / T10;

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
                T0 = (Vgst - Param.BSIM3v1voff) / (n * Vtm);
                ExpVgst = Math.Exp(T0);
                Vgsteff = Vtm * Param.BSIM3v1cdep0 / ModelTemperature.Cox * ExpVgst;
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

                dT2_dVg = -ModelTemperature.Cox / (Vtm * Param.BSIM3v1cdep0)
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

            /* Calculate Effective Channel Geometry */
            T9 = sqrtPhis - Param.BSIM3v1sqrtPhi;
            Weff = Param.BSIM3v1weff - 2.0 * (Param.BSIM3v1dwg * Vgsteff
                 + Param.BSIM3v1dwb * T9);
            dWeff_dVg = -2.0 * Param.BSIM3v1dwg;
            dWeff_dVb = -2.0 * Param.BSIM3v1dwb * dsqrtPhis_dVb;

            if (Weff < 2.0e-8) /* to avoid the discontinuity problem due to Weff*/
            {
                T0 = 1.0 / (6.0e-8 - 2.0 * Weff);
                Weff = 2.0e-8 * (4.0e-8 - Weff) * T0;
                T0 *= T0 * 4.0e-16;
                dWeff_dVg *= T0;
                dWeff_dVb *= T0;
            }

            T0 = Param.BSIM3v1prwg * Vgsteff + Param.BSIM3v1prwb * T9;
            if (T0 >= -0.9)
            {
                Rds = Param.BSIM3v1rds0 * (1.0 + T0);
                dRds_dVg = Param.BSIM3v1rds0 * Param.BSIM3v1prwg;
                dRds_dVb = Param.BSIM3v1rds0 * Param.BSIM3v1prwb * dsqrtPhis_dVb;
            }
            else
            /* to avoid the discontinuity problem due to prwg and prwb*/
            {
                T1 = 1.0 / (17.0 + 20.0 * T0);
                Rds = Param.BSIM3v1rds0 * (0.8 + T0) * T1;
                T1 *= T1;
                dRds_dVg = Param.BSIM3v1rds0 * Param.BSIM3v1prwg * T1;
                dRds_dVb = Param.BSIM3v1rds0 * Param.BSIM3v1prwb * dsqrtPhis_dVb
                 * T1;
            }

            /* Calculate Abulk */
            T1 = 0.5 * Param.BSIM3v1k1 / sqrtPhis;
            dT1_dVb = -T1 / sqrtPhis * dsqrtPhis_dVb;

            T9 = Math.Sqrt(Param.BSIM3v1xj * Xdep);
            tmp1 = Leff + 2.0 * T9;
            T5 = Leff / tmp1;
            tmp2 = Param.BSIM3v1a0 * T5;
            tmp3 = Param.BSIM3v1weff + Param.BSIM3v1b1;
            tmp4 = Param.BSIM3v1b0 / tmp3;
            T2 = tmp2 + tmp4;
            dT2_dVb = -T9 / tmp1 / Xdep * dXdep_dVb;
            T6 = T5 * T5;
            T7 = T5 * T6;

            Abulk0 = 1.0 + T1 * T2;
            dAbulk0_dVb = T1 * tmp2 * dT2_dVb + T2 * dT1_dVb;

            T8 = Param.BSIM3v1ags * Param.BSIM3v1a0 * T7;
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
                dAbulk_dVb *= T9 * T9;
            }

            T2 = Param.BSIM3v1keta * Vbseff;
            if (T2 >= -0.9)
            {
                T0 = 1.0 / (1.0 + T2);
                dT0_dVb = -Param.BSIM3v1keta * T0 * T0;
            }
            else
            /* added to avoid the problems caused by Keta */
            {
                T1 = 1.0 / (0.8 + T2);
                T0 = (17.0 + 20.0 * T2) * T1;
                dT0_dVb = -Param.BSIM3v1keta * T1 * T1;
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
                T2 = Param.BSIM3v1ua + Param.BSIM3v1uc * Vbseff;
                T3 = T0 / ModelParameters.Tox;
                T5 = T3 * (T2 + Param.BSIM3v1ub * T3);
                dDenomi_dVg = (T2 + 2.0 * Param.BSIM3v1ub * T3) / ModelParameters.Tox;
                dDenomi_dVd = dDenomi_dVg * 2.0 * dVth_dVd;
                dDenomi_dVb = dDenomi_dVg * 2.0 * dVth_dVb + Param.BSIM3v1uc * T3;
            }
            else if (ModelParameters.MobMod.Value == 2)
            {
                T5 = Vgsteff / ModelParameters.Tox * (Param.BSIM3v1ua
               + Param.BSIM3v1uc * Vbseff + Param.BSIM3v1ub * Vgsteff
                       / ModelParameters.Tox);
                dDenomi_dVg = (Param.BSIM3v1ua + Param.BSIM3v1uc * Vbseff
                    + 2.0 * Param.BSIM3v1ub * Vgsteff / ModelParameters.Tox)
                    / ModelParameters.Tox;
                dDenomi_dVd = 0.0;
                dDenomi_dVb = Vgsteff * Param.BSIM3v1uc / ModelParameters.Tox;
            }
            else
            {
                T0 = Vgsteff + Vth + Vth;
                T2 = 1.0 + Param.BSIM3v1uc * Vbseff;
                T3 = T0 / ModelParameters.Tox;
                T4 = T3 * (Param.BSIM3v1ua + Param.BSIM3v1ub * T3);
                T5 = T4 * T2;
                dDenomi_dVg = (Param.BSIM3v1ua + 2.0 * Param.BSIM3v1ub * T3) * T2
                    / ModelParameters.Tox;
                dDenomi_dVd = dDenomi_dVg * 2.0 * dVth_dVd;
                dDenomi_dVb = dDenomi_dVg * 2.0 * dVth_dVb + Param.BSIM3v1uc * T4;
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

            this._ueff = ueff = Param.BSIM3v1u0temp / Denomi;
            T9 = -ueff / Denomi;
            dueff_dVg = T9 * dDenomi_dVg;
            dueff_dVd = T9 * dDenomi_dVd;
            dueff_dVb = T9 * dDenomi_dVb;

            /* Saturation Drain Voltage  Vdsat */
            WVCox = Weff * Param.BSIM3v1vsattemp * ModelTemperature.Cox;
            WVCoxRds = WVCox * Rds;

            Esat = 2.0 * Param.BSIM3v1vsattemp / ueff;
            EsatL = Esat * Leff;
            T0 = -EsatL / ueff;
            dEsatL_dVg = T0 * dueff_dVg;
            dEsatL_dVd = T0 * dueff_dVd;
            dEsatL_dVb = T0 * dueff_dVb;

            /* Sqrt() */
            a1 = Param.BSIM3v1a1;
            if (a1 == 0.0)
            {
                Lambda = Param.BSIM3v1a2;
                dLambda_dVg = 0.0;
            }
            else if (a1 > 0.0)
            /* Added to avoid the discontinuity problem
               caused by a1 and a2 (Lambda) */
            {
                T0 = 1.0 - Param.BSIM3v1a2;
                T1 = T0 - Param.BSIM3v1a1 * Vgsteff - 0.0001;
                T2 = Math.Sqrt(T1 * T1 + 0.0004 * T0);
                Lambda = Param.BSIM3v1a2 + T0 - 0.5 * (T1 + T2);
                dLambda_dVg = 0.5 * Param.BSIM3v1a1 * (1.0 + T1 / T2);
            }
            else
            {
                T1 = Param.BSIM3v1a2 + Param.BSIM3v1a1 * Vgsteff - 0.0001;
                T2 = Math.Sqrt(T1 * T1 + 0.0004 * Param.BSIM3v1a2);
                Lambda = 0.5 * (T1 + T2);
                dLambda_dVg = 0.5 * Param.BSIM3v1a1 * (1.0 + T1 / T2);
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

            /* Effective Vds (Vdseff) Calculation */
            T1 = Vdsat - Vds - Param.BSIM3v1delta;
            dT1_dVg = dVdsat_dVg;
            dT1_dVd = dVdsat_dVd - 1.0;
            dT1_dVb = dVdsat_dVb;

            T2 = Math.Sqrt(T1 * T1 + 4.0 * Param.BSIM3v1delta * Vdsat);
            T0 = T1 / T2;
            T3 = 2.0 * Param.BSIM3v1delta / T2;
            dT2_dVg = T0 * dT1_dVg + T3 * dVdsat_dVg;
            dT2_dVd = T0 * dT1_dVd + T3 * dVdsat_dVd;
            dT2_dVb = T0 * dT1_dVb + T3 * dVdsat_dVb;

            Vdseff = Vdsat - 0.5 * (T1 + T2);
            dVdseff_dVg = dVdsat_dVg - 0.5 * (dT1_dVg + dT2_dVg);
            dVdseff_dVd = dVdsat_dVd - 0.5 * (dT1_dVd + dT2_dVd);
            dVdseff_dVb = dVdsat_dVb - 0.5 * (dT1_dVb + dT2_dVb);

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
                Vdseff = Vds; /* This code is added to fixed the problem
			       caused by computer precision when
			       Vds is very close to Vdseff. */
            diffVds = Vds - Vdseff;
            /* Calculate VACLM */
            if ((Param.BSIM3v1pclm > 0.0) && (diffVds > 1.0e-10))
            {
                T0 = 1.0 / (Param.BSIM3v1pclm * Abulk * Param.BSIM3v1litl);
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
            if (Param.BSIM3v1thetaRout > 0.0)
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
                T2 = Param.BSIM3v1thetaRout;
                VADIBL = (Vgst2Vtm - T0 / T1) / T2;
                dVADIBL_dVg = (1.0 - dT0_dVg / T1 + T0 * dT1_dVg / T9) / T2;
                dVADIBL_dVb = (-dT0_dVb / T1 + T0 * dT1_dVb / T9) / T2;
                dVADIBL_dVd = (-dT0_dVd / T1 + T0 * dT1_dVd / T9) / T2;

                T7 = Param.BSIM3v1pdiblb * Vbseff;
                if (T7 >= -0.9)
                {
                    T3 = 1.0 / (1.0 + T7);
                    VADIBL *= T3;
                    dVADIBL_dVg *= T3;
                    dVADIBL_dVb = (dVADIBL_dVb - VADIBL * Param.BSIM3v1pdiblb)
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
                    - VADIBL * Param.BSIM3v1pdiblb * T4 * T4;
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

            T8 = Param.BSIM3v1pvag / EsatL;
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
            if (Param.BSIM3v1pscbe2 > 0.0)
            {
                if (diffVds > Param.BSIM3v1pscbe1 * Param.BSIM3v1litl
                / EXP_THRESHOLD)
                {
                    T0 = Param.BSIM3v1pscbe1 * Param.BSIM3v1litl / diffVds;
                    VASCBE = Leff * Math.Exp(T0) / Param.BSIM3v1pscbe2;
                    T1 = T0 * VASCBE / diffVds;
                    dVASCBE_dVg = T1 * dVdseff_dVg;
                    dVASCBE_dVd = -T1 * (1.0 - dVdseff_dVd);
                    dVASCBE_dVb = T1 * dVdseff_dVb;
                }
                else
                {
                    VASCBE = MAX_EXP * Leff / Param.BSIM3v1pscbe2;
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

            /* calculate substrate current Isub */
            if ((Param.BSIM3v1alpha0 <= 0.0) || (Param.BSIM3v1beta0 <= 0.0))
            {
                Isub = Gbd = Gbb = Gbg = 0.0;
            }
            else
            {
                T2 = Param.BSIM3v1alpha0 / Leff;
                if (diffVds > Param.BSIM3v1beta0 / EXP_THRESHOLD)
                {
                    T0 = -Param.BSIM3v1beta0 / diffVds;
                    T1 = T2 * diffVds * Math.Exp(T0);
                    T3 = T1 / diffVds * (T0 - 1.0);
                    dT1_dVg = T3 * dVdseff_dVg;
                    dT1_dVd = -T3 * (1.0 - dVdseff_dVd);
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
            else if (ModelParameters.CapMod.Value == 0)
            {
                if (Vbseff < 0.0)
                {
                    Vbseff = Vbs;
                    dVbseff_dVb = 1.0;
                }
                else
                {
                    Vbseff = Param.BSIM3v1phi - Phis;
                    dVbseff_dVb = -dPhis_dVb;
                }

                Vfb = Param.BSIM3v1vfbcv;
                Vth = Vfb + Param.BSIM3v1phi + Param.BSIM3v1k1 * sqrtPhis;
                Vgst = Vgs_eff - Vth;
                dVth_dVb = Param.BSIM3v1k1 * dsqrtPhis_dVb;
                dVgst_dVb = -dVth_dVb;
                dVgst_dVg = dVgs_eff_dVg;

                CoxWL = ModelTemperature.Cox * Param.BSIM3v1weffCV
                        * Param.BSIM3v1leffCV;
                Arg1 = Vgs_eff - Vbseff - Vfb;
                if (Arg1 <= 0.0)
                {
                    qgate = CoxWL * Arg1;
                    qbulk = -qgate;
                    qdrn = 0.0;

                    this._cggb = CoxWL * dVgs_eff_dVg;
                    this._cgdb = 0.0;
                    this._cgsb = CoxWL * (dVbseff_dVb
                    - dVgs_eff_dVg);

                    this._cdgb = 0.0;
                    this._cddb = 0.0;
                    this._cdsb = 0.0;

                    this._cbgb = -CoxWL * dVgs_eff_dVg;
                    this._cbdb = 0.0;
                    this._cbsb = -this._cgsb;

                }
                else if (Vgst <= 0.0)
                {
                    T1 = 0.5 * Param.BSIM3v1k1;
                    T2 = Math.Sqrt(T1 * T1 + Arg1);
                    qgate = CoxWL * Param.BSIM3v1k1 * (T2 - T1);
                    qbulk = -qgate;
                    qdrn = 0.0;

                    T0 = CoxWL * T1 / T2;
                    this._cggb = T0 * dVgs_eff_dVg;
                    this._cgdb = 0.0;
                    this._cgsb = T0 * (dVbseff_dVb
                    - dVgs_eff_dVg);

                    this._cdgb = 0.0;
                    this._cddb = 0.0;
                    this._cdsb = 0.0;

                    this._cbgb = -this._cggb;
                    this._cbdb = 0.0;
                    this._cbsb = -this._cgsb;
                }
                else
                {
                    One_Third_CoxWL = CoxWL / 3.0;
                    Two_Third_CoxWL = 2.0 * One_Third_CoxWL;

                    AbulkCV = Abulk0 * Param.BSIM3v1abulkCVfactor;
                    dAbulkCV_dVb = Param.BSIM3v1abulkCVfactor * dAbulk0_dVb;
                    Vdsat = Vgst / AbulkCV;
                    dVdsat_dVg = dVgs_eff_dVg / AbulkCV;
                    dVdsat_dVb = -(Vdsat * dAbulkCV_dVb + dVth_dVb) / AbulkCV;

                    if (ModelParameters.Xpart > 0.5)
                    {   /* 0/100 Charge petition model */
                        if (Vdsat <= Vds)
                        {   /* saturation region */
                            T1 = Vdsat / 3.0;
                            qgate = CoxWL * (Vgs_eff - Vfb
                          - Param.BSIM3v1phi - T1);
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
                          - Param.BSIM3v1phi - 0.5 * (Vds - T3));
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
                    {   /* 40/60 Charge petition model */
                        if (Vds >= Vdsat)
                        {   /* saturation region */
                            T1 = Vdsat / 3.0;
                            qgate = CoxWL * (Vgs_eff - Vfb
                          - Param.BSIM3v1phi - T1);
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
                            qgate = CoxWL * (Vgs_eff - Vfb - Param.BSIM3v1phi
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
                        + this._cddb + T11);
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
                          - Param.BSIM3v1phi - T1);
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
                            qgate = CoxWL * (Vgs_eff - Vfb - Param.BSIM3v1phi
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
                    VbseffCV = Param.BSIM3v1phi - Phis;
                    dVbseffCV_dVb = -dPhis_dVb;
                }

                CoxWL = ModelTemperature.Cox * Param.BSIM3v1weffCV
              * Param.BSIM3v1leffCV;
                Vfb = Vth - Param.BSIM3v1phi - Param.BSIM3v1k1 * sqrtPhis;

                dVfb_dVb = dVth_dVb - Param.BSIM3v1k1 * dsqrtPhis_dVb;
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
                    Arg1 = Vgs_eff - VbseffCV - Vfb - Vgsteff;

                    if (Arg1 <= 0.0)
                    {
                        qgate = CoxWL * Arg1;
                        Cgg = CoxWL * (dVgs_eff_dVg - dVgsteff_dVg);
                        Cgd = -CoxWL * (dVfb_dVd + dVgsteff_dVd);
                        Cgb = -CoxWL * (dVfb_dVb + dVbseffCV_dVb + dVgsteff_dVb);
                    }
                    else
                    {
                        T0 = 0.5 * Param.BSIM3v1k1;
                        T1 = Math.Sqrt(T0 * T0 + Arg1);
                        T2 = CoxWL * T0 / T1;

                        qgate = CoxWL * Param.BSIM3v1k1 * (T1 - T0);

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
                    AbulkCV = Abulk0 * Param.BSIM3v1abulkCVfactor;
                    dAbulkCV_dVb = Param.BSIM3v1abulkCVfactor * dAbulk0_dVb;
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
                else if (ModelParameters.CapMod.Value == 2)

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

                    T0 = 0.5 * Param.BSIM3v1k1;
                    T3 = Vgs_eff - Vfbeff - VbseffCV - Vgsteff;
                    if (Param.BSIM3v1k1 == 0.0)
                    {
                        T1 = 0.0;
                        T2 = 0.0;
                    }
                    else if (T3 < 0.0)
                    {
                        T1 = T0 + T3 / Param.BSIM3v1k1;
                        T2 = CoxWL;
                    }
                    else
                    {
                        T1 = Math.Sqrt(T0 * T0 + T3);
                        T2 = CoxWL * T0 / T1;
                    }

                    Qsub0 = CoxWL * Param.BSIM3v1k1 * (T1 - T0);

                    dQsub0_dVg = T2 * (dVgs_eff_dVg - dVfbeff_dVg - dVgsteff_dVg);
                    dQsub0_dVd = -T2 * (dVfbeff_dVd + dVgsteff_dVd);
                    dQsub0_dVb = -T2 * (dVfbeff_dVb + dVbseffCV_dVb
                               + dVgsteff_dVb);

                    One_Third_CoxWL = CoxWL / 3.0;
                    Two_Third_CoxWL = 2.0 * One_Third_CoxWL;
                    AbulkCV = Abulk0 * Param.BSIM3v1abulkCVfactor;
                    dAbulkCV_dVb = Param.BSIM3v1abulkCVfactor * dAbulk0_dVb;
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

                    T0 = Param.BSIM3v1leffCV * Param.BSIM3v1leffCV;
                    this._tconst = Param.BSIM3v1u0temp * Param.BSIM3v1elm
                      / CoxWL / T0;

                    if (qcheq == 0.0)
                        this._tconst = 0.0;
                    else if (qcheq < 0.0)
                        this._tconst = -this._tconst;

                    gtau_drift = Math.Abs(this._tconst * qcheq);
                    gtau_diff = 16.0 * Param.BSIM3v1u0temp * ModelTemperature.Vtm / T0;

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
               *  czbdsw: zero bias drain junction sidewall capacitance
			  along field oxide
               *  czbssw: zero bias source junction sidewall capacitance
			  along field oxide
	       *  czbdswg: zero bias drain junction sidewall capacitance
			   along gate side
	       *  czbsswg: zero bias source junction sidewall capacitance
			   along gate side
               */

                czbd = ModelParameters.UnitAreaJctCap * Parameters.DrainArea;
                czbs = ModelParameters.UnitAreaJctCap * Parameters.SourceArea;
                if (Parameters.DrainPerimeter < Param.BSIM3v1weff)
                {
                    czbdswg = ModelParameters.UnitLengthGateSidewallJctCap
                   * Parameters.DrainPerimeter;
                    czbdsw = 0.0;
                }
                else
                {
                    czbdsw = ModelParameters.UnitLengthSidewallJctCap
                   * (Parameters.DrainPerimeter - Param.BSIM3v1weff);
                    czbdswg = ModelParameters.UnitLengthGateSidewallJctCap
                        * Param.BSIM3v1weff;
                }
                if (Parameters.SourcePerimeter < Param.BSIM3v1weff)
                {
                    czbssw = 0.0;
                    czbsswg = ModelParameters.UnitLengthGateSidewallJctCap
                          * Parameters.SourcePerimeter;
                }
                else
                {
                    czbssw = ModelParameters.UnitLengthSidewallJctCap
                   * (Parameters.SourcePerimeter - Param.BSIM3v1weff);
                    czbsswg = ModelParameters.UnitLengthGateSidewallJctCap
                        * Param.BSIM3v1weff;
                }
                PhiB = ModelParameters.BulkJctPotential;
                PhiBSW = ModelParameters.SidewallJctPotential;
                PhiBSWG = ModelParameters.GatesidewallJctPotential;
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
                    if (czbsswg > 0.0)
                    {
                        arg = 1.0 - vbs / PhiBSWG;
                        if (MJSWG == 0.5)
                            sarg = 1.0 / Math.Sqrt(arg);
                        else
                            sarg = Math.Exp(-MJSWG * Math.Log(arg));
                        this._qbs += PhiBSWG * czbsswg
                      * (1.0 - arg * sarg) / (1.0 - MJSWG);
                        this._capbs += czbsswg * sarg;
                    }

                }
                else
                {
                    this._qbs = vbs * (czbs + czbssw
                        + czbsswg) + vbs * vbs * (czbs * MJ * 0.5 / PhiB
                                    + czbssw * MJSW * 0.5 / PhiBSW
                        + czbsswg * MJSWG * 0.5 / PhiBSWG);
                    this._capbs = czbs + czbssw + czbsswg + vbs * (czbs * MJ /
                                         PhiB + czbssw * MJSW / PhiBSW + czbsswg * MJSWG / PhiBSWG);
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
                    if (czbdswg > 0.0)
                    {
                        arg = 1.0 - vbd / PhiBSWG;
                        if (MJSWG == 0.5)
                            sarg = 1.0 / Math.Sqrt(arg);
                        else
                            sarg = Math.Exp(-MJSWG * Math.Log(arg));
                        this._qbd += PhiBSWG * czbdswg
                      * (1.0 - arg * sarg) / (1.0 - MJSWG);
                        this._capbd += czbdswg * sarg;
                    }
                }
                else
                {
                    this._qbd = vbd * (czbd + czbdsw
                        + czbdswg) + vbd * vbd * (czbd * MJ * 0.5 / PhiB
                                    + czbdsw * MJSW * 0.5 / PhiBSW
                        + czbdswg * MJSWG * 0.5 / PhiBSWG);
                    this._capbd = czbd + czbdsw + czbdswg + vbd * (czbd * MJ /
                                         PhiB + czbdsw * MJSW / PhiBSW + czbdswg * MJSWG / PhiBSWG);
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
            if (ModelParameters.CapMod.Value == 0)
            {
                cgdo = Param.BSIM3v1cgdo;
                qgdo = Param.BSIM3v1cgdo * vgd;
                cgso = Param.BSIM3v1cgso;
                qgso = Param.BSIM3v1cgso * vgs;
            }
            else if (ModelParameters.CapMod.Value == 1)
            {
                if (vgd < 0.0)
                {
                    T1 = Math.Sqrt(1.0 - 4.0 * vgd / Param.BSIM3v1ckappa);
                    cgdo = Param.BSIM3v1cgdo + Param.BSIM3v1weffCV
                     * Param.BSIM3v1cgdl / T1;
                    qgdo = Param.BSIM3v1cgdo * vgd - Param.BSIM3v1weffCV * 0.5
                     * Param.BSIM3v1cgdl * Param.BSIM3v1ckappa * (T1 - 1.0);
                }
                else
                {
                    cgdo = Param.BSIM3v1cgdo + Param.BSIM3v1weffCV
                     * Param.BSIM3v1cgdl;
                    qgdo = (Param.BSIM3v1weffCV * Param.BSIM3v1cgdl
                     + Param.BSIM3v1cgdo) * vgd;
                }

                if (vgs < 0.0)
                {
                    T1 = Math.Sqrt(1.0 - 4.0 * vgs / Param.BSIM3v1ckappa);
                    cgso = Param.BSIM3v1cgso + Param.BSIM3v1weffCV
                     * Param.BSIM3v1cgsl / T1;
                    qgso = Param.BSIM3v1cgso * vgs - Param.BSIM3v1weffCV * 0.5
                     * Param.BSIM3v1cgsl * Param.BSIM3v1ckappa * (T1 - 1.0);
                }
                else
                {
                    cgso = Param.BSIM3v1cgso + Param.BSIM3v1weffCV
                     * Param.BSIM3v1cgsl;
                    qgso = (Param.BSIM3v1weffCV * Param.BSIM3v1cgsl
                     + Param.BSIM3v1cgso) * vgs;
                }
            }
            else
            {
                T0 = vgd + DELTA_1;
                T1 = Math.Sqrt(T0 * T0 + 4.0 * DELTA_1);
                T2 = 0.5 * (T0 - T1);

                T3 = Param.BSIM3v1weffCV * Param.BSIM3v1cgdl;
                T4 = Math.Sqrt(1.0 - 4.0 * T2 / Param.BSIM3v1ckappa);
                cgdo = Param.BSIM3v1cgdo + T3 - T3 * (1.0 - 1.0 / T4)
                 * (0.5 - 0.5 * T0 / T1);
                qgdo = (Param.BSIM3v1cgdo + T3) * vgd - T3 * (T2
                 + 0.5 * Param.BSIM3v1ckappa * (T4 - 1.0));

                T0 = vgs + DELTA_1;
                T1 = Math.Sqrt(T0 * T0 + 4.0 * DELTA_1);
                T2 = 0.5 * (T0 - T1);
                T3 = Param.BSIM3v1weffCV * Param.BSIM3v1cgsl;
                T4 = Math.Sqrt(1.0 - 4.0 * T2 / Param.BSIM3v1ckappa);
                cgso = Param.BSIM3v1cgso + T3 - T3 * (1.0 - 1.0 / T4)
                 * (0.5 - 0.5 * T0 / T1);
                qgso = (Param.BSIM3v1cgso + T3) * vgs - T3 * (T2
                 + 0.5 * Param.BSIM3v1ckappa * (T4 - 1.0));
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
                gcggb = (this._cggb + cgdo + cgso + Param.BSIM3v1cgbo) * ag0;
                gcgdb = (this._cgdb - cgdo) * ag0;
                gcgsb = (this._cgsb - cgso) * ag0;
                gcbgb = (this._cbgb - Param.BSIM3v1cgbo) * ag0;
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
                qgb = Param.BSIM3v1cgbo * vgb;
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
                gcggb = (this._cggb + cgdo + cgso + Param.BSIM3v1cgbo) * ag0;
                gcgdb = (this._cgsb - cgdo) * ag0;
                gcgsb = (this._cgdb - cgso) * ag0;
                gcbgb = (this._cbgb - Param.BSIM3v1cgbo) * ag0;
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
                qgb = Param.BSIM3v1cgbo * vgb;
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
                this._qcdump.Derive();
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
                goto line1000;
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
                this._gtau = 16.0 * Param.BSIM3v1u0temp * ModelTemperature.Vtm
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
            this._spPtr.Value += m * (cdreq + ceqbs + ceqqg
                           + ceqqb + ceqqd);

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
            _dpgPtr.Value += m * ((Gm + gcdgb) + dxpart *
                                       this._gtg + gbdpg);
            _dpbPtr.Value -= m * ((this._gbd - Gmbs + gcdgb + gcddb
                                           + gcdsb - dxpart * this._gtb) - gbdpb);
            _dpspPtr.Value -= m * ((this._gds + FwdSum - gcdsb
                           - dxpart * this._gts) - gbdpsp);
            _spgPtr.Value += m * (gcsgb - Gm + sxpart *
                                        this._gtg + gbspg);
            _spsPtr.Value -= m * this._sourceConductance;
            _spbPtr.Value -= m * ((this._gbs + Gmbs + gcsgb + gcsdb
                                             + gcssb - sxpart * this._gtb) - gbspb);
            _spdpPtr.Value -= m * ((this._gds + RevSum - gcsdb
                                              - sxpart * this._gtd - this._gbd) - gbspdp);

            _qqPtr.Value += m * ((gqdef + this._gtau));

            _dpqPtr.Value += m * ((dxpart * this._gtau));
            _spqPtr.Value += m * ((sxpart * this._gtau));
            _gqPtr.Value -= m * this._gtau;

            _qgPtr.Value += m * (-gcqgb + this._gtg);
            _qdpPtr.Value += m * (-gcqdb + this._gtd);
            _qspPtr.Value += m * (-gcqsb + this._gts);
            _qbPtr.Value += m * (-gcqbb + this._gtb);

        line1000:;
        }
    }
}
