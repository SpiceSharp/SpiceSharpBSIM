using SpiceSharp;
using SpiceSharp.Algebra;
using SpiceSharp.Attributes;
using SpiceSharp.Behaviors;
using SpiceSharp.Components;
using SpiceSharp.Simulations;
using System;
using System.Collections.Generic;
using System.Numerics;
using System.Text;

namespace SpiceSharpBSIM.Components.Semiconductors.BSIM.BSIM3v2Behaviors
{
    /// <summary>
    /// Frequency behavior for a <see cref="BSIM3v2"/>.
    /// </summary>
    [BehaviorFor(typeof(BSIM3v2)), AddBehaviorIfNo(typeof(IFrequencyBehavior))]
    public partial class FrequencyBehavior : BiasingBehavior, IFrequencyBehavior
    {
        private readonly IComplexSimulationState _state;
        protected readonly IVariable<Complex> _drain, _gate, _source, _bulk, _drainPrime, _sourcePrime, _q;
        private readonly Element<Complex> _ddPtr, _ggPtr, _ssPtr, _bbPtr, _dpdpPtr, _spspPtr, _ddpPtr,
            _gbPtr, _gdpPtr, _gspPtr, _sspPtr, _bdpPtr, _bspPtr, _dpspPtr, _dpdPtr, _bgPtr, _dpgPtr,
            _spgPtr, _spsPtr, _dpbPtr, _spbPtr, _spdpPtr, _qqPtr, _qdpPtr, _qspPtr, _qgPtr, _qbPtr,
            _dpqPtr, _spqPtr, _gqPtr; //, _bqPtr;

        [ParameterName("vbs"), ParameterInfo("Vbs")]
        public Complex ComplexVbs => _bulk.Value - _sourcePrime.Value;
        [ParameterName("vgs"), ParameterInfo("Vgs")]
        public Complex ComplexVgs => _gate.Value - _sourcePrime.Value;
        [ParameterName("vds"), ParameterInfo("Vds")]
        public Complex ComplexVds => _drainPrime.Value - _sourcePrime.Value;

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="name">Name</param>
        public FrequencyBehavior(ComponentBindingContext context)
            : base(context)
        {
            _state = context.GetState<IComplexSimulationState>();
            _drain = _state.GetSharedVariable(context.Nodes[0]);
            _gate = _state.GetSharedVariable(context.Nodes[1]);
            _source = _state.GetSharedVariable(context.Nodes[2]);
            _bulk = _state.GetSharedVariable(context.Nodes[3]);
            if (!DrainConductance.Equals(0.0))
                _drainPrime = _state.CreatePrivateVariable(Name.Combine("drain"), Units.Volt);
            else
                _drainPrime = _drain;
            if (!SourceConductance.Equals(0.0))
                _sourcePrime = _state.CreatePrivateVariable(Name.Combine("source"), Units.Volt);
            else
                _sourcePrime = _source;

            if (Parameters.NqsMod.Value != 0)
                _q = _state.CreatePrivateVariable(Name.Combine("charge"), Units.Coulomb);
            else
                _q = _state.GetSharedVariable(Constants.Ground);

            int drain = _state.Map[_drain];
            int drainPrime = _state.Map[_drainPrime];
            int gate = _state.Map[_gate];
            int source = _state.Map[_source];
            int sourcePrime = _state.Map[_sourcePrime];
            int bulk = _state.Map[_bulk];
            int q = _state.Map[_q];

            _ddPtr = _state.Solver.GetElement(new MatrixLocation(drain, drain));
            _ggPtr = _state.Solver.GetElement(new MatrixLocation(gate, gate));
            _ssPtr = _state.Solver.GetElement(new MatrixLocation(source, source));
            _bbPtr = _state.Solver.GetElement(new MatrixLocation(bulk, bulk));
            _dpdpPtr = _state.Solver.GetElement(new MatrixLocation(drainPrime, drainPrime));
            _spspPtr = _state.Solver.GetElement(new MatrixLocation(sourcePrime, sourcePrime));
            _ddpPtr = _state.Solver.GetElement(new MatrixLocation(drain, drainPrime));
            _gbPtr = _state.Solver.GetElement(new MatrixLocation(gate, bulk));
            _gdpPtr = _state.Solver.GetElement(new MatrixLocation(gate, drainPrime));
            _gspPtr = _state.Solver.GetElement(new MatrixLocation(gate, sourcePrime));
            _sspPtr = _state.Solver.GetElement(new MatrixLocation(source, sourcePrime));
            _bdpPtr = _state.Solver.GetElement(new MatrixLocation(bulk, drainPrime));
            _bspPtr = _state.Solver.GetElement(new MatrixLocation(bulk, sourcePrime));
            _dpspPtr = _state.Solver.GetElement(new MatrixLocation(drainPrime, sourcePrime));
            _dpdPtr = _state.Solver.GetElement(new MatrixLocation(drainPrime, drain));
            _bgPtr = _state.Solver.GetElement(new MatrixLocation(bulk, gate));
            _dpgPtr = _state.Solver.GetElement(new MatrixLocation(drainPrime, gate));
            _spgPtr = _state.Solver.GetElement(new MatrixLocation(sourcePrime, gate));
            _spsPtr = _state.Solver.GetElement(new MatrixLocation(sourcePrime, source));
            _dpbPtr = _state.Solver.GetElement(new MatrixLocation(drainPrime, bulk));
            _spbPtr = _state.Solver.GetElement(new MatrixLocation(sourcePrime, bulk));
            _spdpPtr = _state.Solver.GetElement(new MatrixLocation(sourcePrime, drainPrime));

            _qqPtr = _state.Solver.GetElement(new MatrixLocation(q, q));

            _qdpPtr = _state.Solver.GetElement(new MatrixLocation(q, drainPrime));
            _qspPtr = _state.Solver.GetElement(new MatrixLocation(q, sourcePrime));
            _qgPtr = _state.Solver.GetElement(new MatrixLocation(q, gate));
            _qbPtr = _state.Solver.GetElement(new MatrixLocation(q, bulk));
            _dpqPtr = _state.Solver.GetElement(new MatrixLocation(drainPrime, q));
            _spqPtr = _state.Solver.GetElement(new MatrixLocation(sourcePrime, q));
            _gqPtr = _state.Solver.GetElement(new MatrixLocation(gate, q));
            // Element is never used...
            // _bqPtr = _state.Solver.GetElement(new MatrixLocation(bulk, q));
        }

        /// <inheritdoc />
        void IFrequencyBehavior.InitializeParameters()
        {
            InitializeSmallSignal = true;
            ((IBiasingBehavior)this).Load();
            InitializeSmallSignal = false;
        }

        /// <inheritdoc />
        void IFrequencyBehavior.Load()
        {
            double xcggb, xcgdb, xcgsb, xcbgb, xcbdb, xcbsb, xcddb, xcssb, xcdgb;
            double gdpr, gspr, gds, gbd, gbs, capbd, capbs, xcsgb, xcdsb, xcsdb;
            double cggb, cgdb, cgsb, cbgb, cbdb, cbsb, cddb, cdgb, cdsb, omega;
            double GSoverlapCap, GDoverlapCap, GBoverlapCap, FwdSum, RevSum, Gm, Gmbs;
            double dxpart, sxpart, xgtg, xgtd, xgts, xgtb;
            double xcqgb = 0.0, xcqdb = 0.0, xcqsb = 0.0, xcqbb = 0.0;
            double gbspsp, gbbdp, gbbsp, gbspg, gbspb;
            double gbspdp, gbdpdp, gbdpg, gbdpb, gbdpsp;
            double ddxpart_dVd, ddxpart_dVg, ddxpart_dVb, ddxpart_dVs;
            double dsxpart_dVd, dsxpart_dVg, dsxpart_dVb, dsxpart_dVs;
            double T1, CoxWL, qcheq, Cdg, Cdd, Cds, Csg, Csd, Css;
            double ScalingFactor = 1.0e-9;
            double m;
            omega = _state.Laplace.Imaginary;

            if (this._mode >= 0)
            {
                Gm = this._gm;
                Gmbs = this._gmbs;
                FwdSum = Gm + Gmbs;
                RevSum = 0.0;

                gbbdp = -this._gbds;
                gbbsp = this._gbds + this._gbgs + this._gbbs;

                gbdpg = this._gbgs;
                gbdpb = this._gbbs;
                gbdpdp = this._gbds;
                gbdpsp = -(gbdpg + gbdpb + gbdpdp);

                gbspdp = 0.0;
                gbspg = 0.0;
                gbspb = 0.0;
                gbspsp = 0.0;

                if (Parameters.NqsMod.Value == 0)
                {
                    cggb = this._cggb;
                    cgsb = this._cgsb;
                    cgdb = this._cgdb;

                    cbgb = this._cbgb;
                    cbsb = this._cbsb;
                    cbdb = this._cbdb;

                    cdgb = this._cdgb;
                    cdsb = this._cdsb;
                    cddb = this._cddb;

                    xgtg = xgtd = xgts = xgtb = 0.0;
                    sxpart = 0.6;
                    dxpart = 0.4;
                    ddxpart_dVd = ddxpart_dVg = ddxpart_dVb
                                = ddxpart_dVs = 0.0;
                    dsxpart_dVd = dsxpart_dVg = dsxpart_dVb
                                = dsxpart_dVs = 0.0;
                }
                else
                {
                    cggb = cgdb = cgsb = 0.0;
                    cbgb = cbdb = cbsb = 0.0;
                    cdgb = cddb = cdsb = 0.0;

                    xgtg = this._gtg;
                    xgtd = this._gtd;
                    xgts = this._gts;
                    xgtb = this._gtb;

                    xcqgb = this._cqgb * omega;
                    xcqdb = this._cqdb * omega;
                    xcqsb = this._cqsb * omega;
                    xcqbb = this._cqbb * omega;

                    CoxWL = ModelTemperature.Cox * Param.BSIM3v32weffCV
                          * Param.BSIM3v32leffCV;
                    qcheq = -(this._qgate + this._qbulk);
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
                        dxpart = this._qdrn / qcheq;
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

                        ddxpart_dVb = -(ddxpart_dVd + ddxpart_dVg
                                    + ddxpart_dVs);
                    }
                    sxpart = 1.0 - dxpart;
                    dsxpart_dVd = -ddxpart_dVd;
                    dsxpart_dVg = -ddxpart_dVg;
                    dsxpart_dVs = -ddxpart_dVs;
                    dsxpart_dVb = -(dsxpart_dVd + dsxpart_dVg + dsxpart_dVs);
                }
            }
            else
            {
                Gm = -this._gm;
                Gmbs = -this._gmbs;
                FwdSum = 0.0;
                RevSum = -(Gm + Gmbs);

                gbbsp = -this._gbds;
                gbbdp = this._gbds + this._gbgs + this._gbbs;

                gbdpg = 0.0;
                gbdpsp = 0.0;
                gbdpb = 0.0;
                gbdpdp = 0.0;

                gbspg = this._gbgs;
                gbspsp = this._gbds;
                gbspb = this._gbbs;
                gbspdp = -(gbspg + gbspsp + gbspb);

                if (Parameters.NqsMod.Value == 0)
                {
                    cggb = this._cggb;
                    cgsb = this._cgdb;
                    cgdb = this._cgsb;

                    cbgb = this._cbgb;
                    cbsb = this._cbdb;
                    cbdb = this._cbsb;

                    cdgb = -(this._cdgb + cggb + cbgb);
                    cdsb = -(this._cddb + cgsb + cbsb);
                    cddb = -(this._cdsb + cgdb + cbdb);

                    xgtg = xgtd = xgts = xgtb = 0.0;
                    sxpart = 0.4;
                    dxpart = 0.6;
                    ddxpart_dVd = ddxpart_dVg = ddxpart_dVb
                                = ddxpart_dVs = 0.0;
                    dsxpart_dVd = dsxpart_dVg = dsxpart_dVb
                                = dsxpart_dVs = 0.0;
                }
                else
                {
                    cggb = cgdb = cgsb = 0.0;
                    cbgb = cbdb = cbsb = 0.0;
                    cdgb = cddb = cdsb = 0.0;

                    xgtg = this._gtg;
                    xgtd = this._gts;
                    xgts = this._gtd;
                    xgtb = this._gtb;

                    xcqgb = this._cqgb * omega;
                    xcqdb = this._cqsb * omega;
                    xcqsb = this._cqdb * omega;
                    xcqbb = this._cqbb * omega;

                    CoxWL = ModelTemperature.Cox * Param.BSIM3v32weffCV
                          * Param.BSIM3v32leffCV;
                    qcheq = -(this._qgate + this._qbulk);
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
                        sxpart = this._qdrn / qcheq;
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

                        dsxpart_dVb = -(dsxpart_dVd + dsxpart_dVg
                                    + dsxpart_dVs);
                    }
                    dxpart = 1.0 - sxpart;
                    ddxpart_dVd = -dsxpart_dVd;
                    ddxpart_dVg = -dsxpart_dVg;
                    ddxpart_dVs = -dsxpart_dVs;
                    ddxpart_dVb = -(ddxpart_dVd + ddxpart_dVg + ddxpart_dVs);
                }
            }

            T1 = this._qdef * this._gtau;
            gdpr = this.DrainConductance;
            gspr = this.SourceConductance;
            gds = this._gds;
            gbd = this._gbd;
            gbs = this._gbs;
            capbd = this._capbd;
            capbs = this._capbs;

            GSoverlapCap = this.Cgso;
            GDoverlapCap = this.Cgdo;
            GBoverlapCap = Param.BSIM3v32cgbo;

            xcdgb = (cdgb - GDoverlapCap) * omega;
            xcddb = (cddb + capbd + GDoverlapCap) * omega;
            xcdsb = cdsb * omega;
            xcsgb = -(cggb + cbgb + cdgb + GSoverlapCap) * omega;
            xcsdb = -(cgdb + cbdb + cddb) * omega;
            xcssb = (capbs + GSoverlapCap - (cgsb + cbsb + cdsb)) * omega;
            xcggb = (cggb + GDoverlapCap + GSoverlapCap + GBoverlapCap)
                  * omega;
            xcgdb = (cgdb - GDoverlapCap) * omega;
            xcgsb = (cgsb - GSoverlapCap) * omega;
            xcbgb = (cbgb - GBoverlapCap) * omega;
            xcbdb = (cbdb - capbd) * omega;
            xcbsb = (cbsb - capbs) * omega;

            m = Parameters.M;

            this._ggPtr.Value += new Complex(0, m * xcggb);
            this._bbPtr.Value -= new Complex(0,
                    m * (xcbgb + xcbdb + xcbsb));
            this._dpdpPtr.Value += new Complex(0, m * xcddb);
            this._spspPtr.Value += new Complex(0, m * xcssb);
            this._gbPtr.Value -= new Complex(0,
                    m * (xcggb + xcgdb + xcgsb));
            this._gdpPtr.Value += new Complex(0, m * xcgdb);
            this._gspPtr.Value += new Complex(0, m * xcgsb);
            this._bgPtr.Value += new Complex(0, m * xcbgb);
            this._bdpPtr.Value += new Complex(0, m * xcbdb);
            this._bspPtr.Value += new Complex(0, m * xcbsb);
            this._dpgPtr.Value += new Complex(0, m * xcdgb);
            this._dpbPtr.Value -= new Complex(0,
                    m * (xcdgb + xcddb + xcdsb));
            this._dpspPtr.Value += new Complex(0, m * xcdsb);
            this._spgPtr.Value += new Complex(0, m * xcsgb);
            this._spbPtr.Value -= new Complex(0,
                    m * (xcsgb + xcsdb + xcssb));
            this._spdpPtr.Value += new Complex(0, m * xcsdb);

            this._ddPtr.Value += m * gdpr;
            this._ssPtr.Value += m * gspr;
            this._bbPtr.Value +=
                    m * (gbd + gbs - this._gbbs);
            this._dpdpPtr.Value +=
                    m * (gdpr + gds + gbd + RevSum +
                         dxpart * xgtd + T1 * ddxpart_dVd +
                         gbdpdp);
            this._spspPtr.Value +=
                    m * (gspr + gds + gbs + FwdSum +
                         sxpart * xgts + T1 * dsxpart_dVs +
                         gbspsp);

            this._ddpPtr.Value -= m * gdpr;
            this._sspPtr.Value -= m * gspr;

            this._bgPtr.Value -= m * this._gbgs;
            this._bdpPtr.Value -= m * (gbd - gbbdp);
            this._bspPtr.Value -= m * (gbs - gbbsp);

            this._dpdPtr.Value -= m * gdpr;
            this._dpgPtr.Value +=
                    m * (Gm + dxpart * xgtg + T1 * ddxpart_dVg +
                         gbdpg);
            this._dpbPtr.Value -=
                    m * (gbd - Gmbs - dxpart * xgtb -
                         T1 * ddxpart_dVb - gbdpb);
            this._dpspPtr.Value -=
                    m * (gds + FwdSum - dxpart * xgts -
                         T1 * ddxpart_dVs - gbdpsp);

            this._spgPtr.Value -=
                    m * (Gm - sxpart * xgtg - T1 * dsxpart_dVg -
                         gbspg);
            this._spsPtr.Value -= m * gspr;
            this._spbPtr.Value -=
                    m * (gbs + Gmbs - sxpart * xgtb -
                         T1 * dsxpart_dVb - gbspb);
            this._spdpPtr.Value -=
                    m * (gds + RevSum - sxpart * xgtd -
                         T1 * dsxpart_dVd - gbspdp);

            this._ggPtr.Value -= m * xgtg;
            this._gbPtr.Value -= m * xgtb;
            this._gdpPtr.Value -= m * xgtd;
            this._gspPtr.Value -= m * xgts;

            if (Parameters.NqsMod.Value != 0)
            {
                this._qqPtr.Value += new Complex(0,
                        m * omega * ScalingFactor);
                this._qgPtr.Value -= new Complex(0, m * xcqgb);
                this._qdpPtr.Value -= new Complex(0, m * xcqdb);
                this._qspPtr.Value -= new Complex(0, m * xcqsb);
                this._qbPtr.Value -= new Complex(0, m * xcqbb);

                this._qqPtr.Value += m * this._gtau;

                this._dpqPtr.Value +=
                        m * (dxpart * this._gtau);
                this._spqPtr.Value +=
                        m * (sxpart * this._gtau);
                this._gqPtr.Value -= m * this._gtau;

                this._qgPtr.Value += m * xgtg;
                this._qdpPtr.Value += m * xgtd;
                this._qspPtr.Value += m * xgts;
                this._qbPtr.Value += m * xgtb;
            }
        }
    }
}
