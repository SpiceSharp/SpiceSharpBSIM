using System;
using System.Numerics;
using SpiceSharp;
using SpiceSharp.Algebra;
using SpiceSharp.Attributes;
using SpiceSharp.Behaviors;
using SpiceSharp.Components;
using SpiceSharp.Simulations;

namespace SpiceSharpBSIM.Components.Semiconductors.BSIM.BSIM3Behaviors
{
    /// <summary>
    /// Frequency behavior for a <see cref="BSIM3"/>
    /// </summary>
    [BehaviorFor(typeof(BSIM3)), AddBehaviorIfNo(typeof(IFrequencyBehavior))]
    public class FrequencyBehavior : BiasingBehavior, IFrequencyBehavior
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

        /// <summary>
        /// Initializes the parameters.
        /// </summary>
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
            double dxpart, sxpart, xgtg, xgtd, xgts, xgtb, xcqgb = 0.0, xcqdb = 0.0, xcqsb = 0.0, xcqbb = 0.0;
            double gbspsp, gbbdp, gbbsp, gbspg, gbspb;
            double gbspdp, gbdpdp, gbdpg, gbdpb, gbdpsp;
            double ddxpart_dVd, ddxpart_dVg, ddxpart_dVb, ddxpart_dVs;
            double dsxpart_dVd, dsxpart_dVg, dsxpart_dVb, dsxpart_dVs;
            double T1, CoxWL, qcheq, Cdg, Cdd, Cds, Csg, Csd, Css;
            double ScalingFactor = 1.0e-9;
            /* For ACNQSMOD */
            double T0, T2, T3, gmr, gmbsr, gmi, gmbsi, gdsi;
            double Cddr, Cdgr, Cdsr, Csdr, Csgr, Cssr, Cgdr, Cggr, Cgsr;
            double Cddi, Cdgi, Cdsi, Cdbi, Csdi, Csgi, Cssi, Csbi;
            double Cgdi, Cggi, Cgsi, Cgbi, Gmi, Gmbsi, FwdSumi, RevSumi;
            double xcdgbi, xcsgbi, xcddbi, xcdsbi, xcsdbi, xcssbi, xcdbbi;
            double xcsbbi, xcggbi, xcgdbi, xcgsbi, xcgbbi;
            double m;

            omega = _state.Laplace.Imaginary;

            Csd = -(this._cddb + this._cgdb + this._cbdb);
            Csg = -(this._cdgb + this._cggb + this._cbgb);
            Css = -(this._cdsb + this._cgsb + this._cbsb);

            if (Parameters.AcnqsMod.Value != 0)
            {
                T0 = omega * this._taunet;
                T1 = T0 * T0;
                T2 = 1.0 / (1.0 + T1);
                T3 = T0 * T2;

                gmr = this._gm * T2;
                gmbsr = this._gmbs * T2;
                gds = this._gds * T2;

                gmi = -this._gm * T3;
                gmbsi = -this._gmbs * T3;
                gdsi = -this._gds * T3;

                Cddr = this._cddb * T2;
                Cdgr = this._cdgb * T2;
                Cdsr = this._cdsb * T2;

                Cddi = this._cddb * T3 * omega;
                Cdgi = this._cdgb * T3 * omega;
                Cdsi = this._cdsb * T3 * omega;
                Cdbi = -(Cddi + Cdgi + Cdsi);

                Csdr = Csd * T2;
                Csgr = Csg * T2;
                Cssr = Css * T2;

                Csdi = Csd * T3 * omega;
                Csgi = Csg * T3 * omega;
                Cssi = Css * T3 * omega;
                Csbi = -(Csdi + Csgi + Cssi);

                Cgdr = -(Cddr + Csdr + this._cbdb);
                Cggr = -(Cdgr + Csgr + this._cbgb);
                Cgsr = -(Cdsr + Cssr + this._cbsb);

                Cgdi = -(Cddi + Csdi);
                Cggi = -(Cdgi + Csgi);
                Cgsi = -(Cdsi + Cssi);
                Cgbi = -(Cgdi + Cggi + Cgsi);
            }
            else /* QS */
            {
                gmr = this._gm;
                gmbsr = this._gmbs;
                gds = this._gds;
                gmi = gmbsi = gdsi = 0.0;

                Cddr = this._cddb;
                Cdgr = this._cdgb;
                Cdsr = this._cdsb;
                Cddi = Cdgi = Cdsi = Cdbi = 0.0;

                Csdr = Csd;
                Csgr = Csg;
                Cssr = Css;
                Csdi = Csgi = Cssi = Csbi = 0.0;

                Cgdr = this._cgdb;
                Cggr = this._cggb;
                Cgsr = this._cgsb;
                Cgdi = Cggi = Cgsi = Cgbi = 0.0;
            }

            if (this._mode >= 0)
            {
                Gm = gmr;
                Gmbs = gmbsr;
                FwdSum = Gm + Gmbs;
                RevSum = 0.0;
                Gmi = gmi;
                Gmbsi = gmbsi;
                FwdSumi = Gmi + Gmbsi;
                RevSumi = 0.0;

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

                if (Parameters.NqsMod.Value == 0 || Parameters.AcnqsMod.Value == 1)
                {
                    cggb = Cggr;
                    cgsb = Cgsr;
                    cgdb = Cgdr;

                    cbgb = this._cbgb;
                    cbsb = this._cbsb;
                    cbdb = this._cbdb;

                    cdgb = Cdgr;
                    cdsb = Cdsr;
                    cddb = Cddr;

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

                    CoxWL = ModelTemperature.Cox * Param.BSIM3weffCV
                                  * Param.BSIM3leffCV;
                    qcheq = -(this._qgate + this._qbulk);
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
                xcdgbi = Cdgi;
                xcsgbi = Csgi;
                xcddbi = Cddi;
                xcdsbi = Cdsi;
                xcsdbi = Csdi;
                xcssbi = Cssi;
                xcdbbi = Cdbi;
                xcsbbi = Csbi;
                xcggbi = Cggi;
                xcgdbi = Cgdi;
                xcgsbi = Cgsi;
                xcgbbi = Cgbi;
            }
            else
            {
                Gm = -gmr;
                Gmbs = -gmbsr;
                FwdSum = 0.0;
                RevSum = -(Gm + Gmbs);
                Gmi = -gmi;
                Gmbsi = -gmbsi;
                FwdSumi = 0.0;
                RevSumi = -(Gmi + Gmbsi);

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

                if (Parameters.NqsMod.Value == 0 || Parameters.AcnqsMod.Value == 1)
                {
                    cggb = Cggr;
                    cgsb = Cgdr;
                    cgdb = Cgsr;

                    cbgb = this._cbgb;
                    cbsb = this._cbdb;
                    cbdb = this._cbsb;

                    cdgb = -(Cdgr + cggb + cbgb);
                    cdsb = -(Cddr + cgsb + cbsb);
                    cddb = -(Cdsr + cgdb + cbdb);

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

                    CoxWL = ModelTemperature.Cox * Param.BSIM3weffCV
                                  * Param.BSIM3leffCV;
                    qcheq = -(this._qgate + this._qbulk);
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
                xcdgbi = Csgi;
                xcsgbi = Cdgi;
                xcddbi = Cssi;
                xcdsbi = Csdi;
                xcsdbi = Cdsi;
                xcssbi = Cddi;
                xcdbbi = Csbi;
                xcsbbi = Cdbi;
                xcggbi = Cggi;
                xcgdbi = Cgsi;
                xcgsbi = Cgdi;
                xcgbbi = Cgbi;
            }

            T1 = this._qdef * this._gtau;
            gdpr = this.DrainConductance;
            gspr = this.SourceConductance;
            gbd = this._gbd;
            gbs = this._gbs;
            capbd = this._capbd;
            capbs = this._capbs;

            GSoverlapCap = this.Cgso;
            GDoverlapCap = this.Cgdo;
            GBoverlapCap = Param.BSIM3cgbo;

            xcdgb = (cdgb - GDoverlapCap) * omega;
            xcddb = (cddb + capbd + GDoverlapCap) * omega;
            xcdsb = cdsb * omega;
            xcsgb = -(cggb + cbgb + cdgb + GSoverlapCap) * omega;
            xcsdb = -(cgdb + cbdb + cddb) * omega;
            xcssb = (capbs + GSoverlapCap - (cgsb + cbsb + cdsb)) * omega;
            xcggb = (cggb + GDoverlapCap + GSoverlapCap + GBoverlapCap) * omega;
            xcgdb = (cgdb - GDoverlapCap) * omega;
            xcgsb = (cgsb - GSoverlapCap) * omega;
            xcbgb = (cbgb - GBoverlapCap) * omega;
            xcbdb = (cbdb - capbd) * omega;
            xcbsb = (cbsb - capbs) * omega;

            m = Parameters.M;

            _ggPtr.Value += new Complex(0, m * xcggb);
            _bbPtr.Value -= new Complex(0, m * (xcbgb + xcbdb + xcbsb));
            _dpdpPtr.Value += new Complex(0, m * (xcddb + gdsi + RevSumi));
            _spspPtr.Value += new Complex(0, m * (xcssb + gdsi + FwdSumi));
            _gbPtr.Value -= new Complex(0, m * (xcggb + xcgdb + xcgsb));
            _gdpPtr.Value += new Complex(0, m * xcgdb);
            _gspPtr.Value += new Complex(0, m * xcgsb);
            _bgPtr.Value += new Complex(0, m * xcbgb);
            _bdpPtr.Value += new Complex(0, m * xcbdb);
            _bspPtr.Value += new Complex(0, m * xcbsb);
            _dpgPtr.Value += new Complex(0, m * (xcdgb + Gmi));
            _dpbPtr.Value -= new Complex(0, m * (xcdgb + xcddb + xcdsb + Gmbsi));
            _dpspPtr.Value += new Complex(0, m * (xcdsb - gdsi - FwdSumi));
            _spgPtr.Value += new Complex(0, m * (xcsgb - Gmi));
            _spbPtr.Value -= new Complex(0, m * (xcsgb + xcsdb + xcssb - Gmbsi));
            _spdpPtr.Value += new Complex(0, m * (xcsdb - gdsi - RevSumi));

            _ddPtr.Value += m * gdpr;
            _ssPtr.Value += m * gspr;
            _bbPtr.Value += m * (gbd + gbs - this._gbbs);
            _dpdpPtr.Value += m * (gdpr + gds + gbd + RevSum + xcddbi
                                  + dxpart * xgtd + T1 * ddxpart_dVd + gbdpdp);
            _spspPtr.Value += m * (gspr + gds + gbs + FwdSum + xcssbi
                                  + sxpart * xgts + T1 * dsxpart_dVs + gbspsp);

            _ddpPtr.Value -= m * gdpr;
            _sspPtr.Value -= m * gspr;

            _bgPtr.Value -= m * this._gbgs;
            _bdpPtr.Value -= m * (gbd - gbbdp);
            _bspPtr.Value -= m * (gbs - gbbsp);

            _dpdPtr.Value -= m * gdpr;
            _dpgPtr.Value += m * (Gm + dxpart * xgtg + T1 * ddxpart_dVg
                 + gbdpg + xcdgbi);
            _dpbPtr.Value -= m * (gbd - Gmbs - dxpart * xgtb
                 - T1 * ddxpart_dVb - gbdpb - xcdbbi);
            _dpspPtr.Value -= m * (gds + FwdSum - dxpart * xgts
                  - T1 * ddxpart_dVs - gbdpsp - xcdsbi);

            _spgPtr.Value -= m * (Gm - sxpart * xgtg - T1 * dsxpart_dVg
                 - gbspg - xcsgbi);
            _spsPtr.Value -= m * gspr;
            _spbPtr.Value -= m * (gbs + Gmbs - sxpart * xgtb
                 - T1 * dsxpart_dVb - gbspb - xcsbbi);
            _spdpPtr.Value -= m * (gds + RevSum - sxpart * xgtd
                  - T1 * dsxpart_dVd - gbspdp - xcsdbi);

            _ggPtr.Value -= m * (xgtg - xcggbi);
            _gbPtr.Value -= m * (xgtb - xcgbbi);
            _gdpPtr.Value -= m * (xgtd - xcgdbi);
            _gspPtr.Value -= m * (xgts - xcgsbi);

            if (Parameters.NqsMod.Value != 0)
            {
                if (Parameters.AcnqsMod.Value != 0)
                {
                    _qqPtr.Value += m * 1.0;
                    _qgPtr.Value += 0.0;
                    _qdpPtr.Value += 0.0;
                    _qspPtr.Value += 0.0;
                    _qbPtr.Value += 0.0;

                    _dpqPtr.Value += 0.0;
                    _spqPtr.Value += 0.0;
                    _gqPtr.Value += 0.0;

                }
                else
                {
                    _qqPtr.Value += new Complex(0, m * omega * ScalingFactor);
                    _qgPtr.Value -= new Complex(0, m * xcqgb);
                    _qdpPtr.Value -= new Complex(0, m * xcqdb);
                    _qspPtr.Value -= new Complex(0, m * xcqsb);
                    _qbPtr.Value -= new Complex(0, m * xcqbb);

                    _qqPtr.Value += m * this._gtau;

                    _dpqPtr.Value += m * dxpart * this._gtau;
                    _spqPtr.Value += m * sxpart * this._gtau;
                    _gqPtr.Value -= m * this._gtau;

                    _qgPtr.Value += m * xgtg;
                    _qdpPtr.Value += m * xgtd;
                    _qspPtr.Value += m * xgts;
                    _qbPtr.Value += m * xgtb;
                }
            }
        }
    }
}