using System;
using System.Numerics;
using SpiceSharp.Algebra;
using SpiceSharp.Behaviors;
using SpiceSharp.Simulations;

namespace SpiceSharp.Components.BSIM3Behaviors
{
    /// <summary>
    /// Frequency behavior for a <see cref="BSIM3"/>
    /// </summary>
    public class FrequencyBehavior : BaseFrequencyBehavior, IConnectedBehavior
    {
        private const double ScalingFactor = 1.0e-9;

        /// <summary>
        /// Necessary behaviors and parameter sets
        /// </summary>
        private BaseParameters _bp;
        private ModelBaseParameters _mbp;
        private TemperatureBehavior _temp;
        private LoadBehavior _load;

        /// <summary>
        /// Nodes
        /// </summary>
        private int _drainNode, _gateNode, _sourceNode, _bulkNode, _drainNodePrime, _sourceNodePrime, _qNode;
        protected MatrixElement<Complex> DdPtr { get; private set; }
        protected MatrixElement<Complex> GgPtr { get; private set; }
        protected MatrixElement<Complex> SsPtr { get; private set; }
        protected MatrixElement<Complex> BbPtr { get; private set; }
        protected MatrixElement<Complex> DPdpPtr { get; private set; }
        protected MatrixElement<Complex> SPspPtr { get; private set; }
        protected MatrixElement<Complex> DdpPtr { get; private set; }
        protected MatrixElement<Complex> GbPtr { get; private set; }
        protected MatrixElement<Complex> GdpPtr { get; private set; }
        protected MatrixElement<Complex> GspPtr { get; private set; }
        protected MatrixElement<Complex> SspPtr { get; private set; }
        protected MatrixElement<Complex> BdpPtr { get; private set; }
        protected MatrixElement<Complex> BspPtr { get; private set; }
        protected MatrixElement<Complex> DPspPtr { get; private set; }
        protected MatrixElement<Complex> DPdPtr { get; private set; }
        protected MatrixElement<Complex> BgPtr { get; private set; }
        protected MatrixElement<Complex> DPgPtr { get; private set; }
        protected MatrixElement<Complex> SPgPtr { get; private set; }
        protected MatrixElement<Complex> SPsPtr { get; private set; }
        protected MatrixElement<Complex> DPbPtr { get; private set; }
        protected MatrixElement<Complex> SPbPtr { get; private set; }
        protected MatrixElement<Complex> SPdpPtr { get; private set; }
        protected MatrixElement<Complex> QqPtr { get; private set; }
        protected MatrixElement<Complex> QdpPtr { get; private set; }
        protected MatrixElement<Complex> QspPtr { get; private set; }
        protected MatrixElement<Complex> QgPtr { get; private set; }
        protected MatrixElement<Complex> QbPtr { get; private set; }
        protected MatrixElement<Complex> DPqPtr { get; private set; }
        protected MatrixElement<Complex> SPqPtr { get; private set; }
        protected MatrixElement<Complex> GqPtr { get; private set; }
        protected MatrixElement<Complex> BqPtr { get; private set; }
        protected VectorElement<Complex> GateNodePtr { get; private set; }
        protected VectorElement<Complex> BulkNodePtr { get; private set; }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="name">Name</param>
        public FrequencyBehavior(Identifier name) : base(name)
        {
        }

        public override void Setup(SetupDataProvider provider)
        {
            if (provider == null)
                throw new ArgumentNullException(nameof(provider));

            // Get parameters
            _mbp = provider.GetParameterSet<ModelBaseParameters>("model");
            _bp = provider.GetParameterSet<BaseParameters>("entity");

            // Get behaviors
            _temp = provider.GetBehavior<TemperatureBehavior>("entity");
            _load = provider.GetBehavior<LoadBehavior>("entity");
        }

        /// <summary>
        /// Connect the behavior
        /// </summary>
        /// <param name="pins">Pins</param>
        public void Connect(params int[] pins)
        {
            _drainNode = pins[0];
            _gateNode = pins[1];
            _sourceNode = pins[2];
            _bulkNode = pins[3];
        }

        /// <summary>
        /// Get equation pointers
        /// </summary>
        /// <param name="solver">Solver</param>
        public override void GetEquationPointers(Solver<Complex> solver)
        {
            _sourceNodePrime = _load.SourceNodePrime;
            _drainNodePrime = _load.DrainNodePrime;
            _qNode = _load.QNode;

            DdPtr = solver.GetMatrixElement(_drainNode, _drainNode);
            GgPtr = solver.GetMatrixElement(_gateNode, _gateNode);
            SsPtr = solver.GetMatrixElement(_sourceNode, _sourceNode);
            BbPtr = solver.GetMatrixElement(_bulkNode, _bulkNode);
            DPdpPtr = solver.GetMatrixElement(_drainNodePrime, _drainNodePrime);
            SPspPtr = solver.GetMatrixElement(_sourceNodePrime, _sourceNodePrime);
            DdpPtr = solver.GetMatrixElement(_drainNode, _drainNodePrime);
            GbPtr = solver.GetMatrixElement(_gateNode, _bulkNode);
            GdpPtr = solver.GetMatrixElement(_gateNode, _drainNodePrime);
            GspPtr = solver.GetMatrixElement(_gateNode, _sourceNodePrime);
            SspPtr = solver.GetMatrixElement(_sourceNode, _sourceNodePrime);
            BdpPtr = solver.GetMatrixElement(_bulkNode, _drainNodePrime);
            BspPtr = solver.GetMatrixElement(_bulkNode, _sourceNodePrime);
            DPspPtr = solver.GetMatrixElement(_drainNodePrime, _sourceNodePrime);
            DPdPtr = solver.GetMatrixElement(_drainNodePrime, _drainNode);
            BgPtr = solver.GetMatrixElement(_bulkNode, _gateNode);
            DPgPtr = solver.GetMatrixElement(_drainNodePrime, _gateNode);
            SPgPtr = solver.GetMatrixElement(_sourceNodePrime, _gateNode);
            SPsPtr = solver.GetMatrixElement(_sourceNodePrime, _sourceNode);
            DPbPtr = solver.GetMatrixElement(_drainNodePrime, _bulkNode);
            SPbPtr = solver.GetMatrixElement(_sourceNodePrime, _bulkNode);
            SPdpPtr = solver.GetMatrixElement(_sourceNodePrime, _drainNodePrime);
            QqPtr = solver.GetMatrixElement(_qNode, _qNode);
            QdpPtr = solver.GetMatrixElement(_qNode, _drainNodePrime);
            QspPtr = solver.GetMatrixElement(_qNode, _sourceNodePrime);
            QgPtr = solver.GetMatrixElement(_qNode, _gateNode);
            QbPtr = solver.GetMatrixElement(_qNode, _bulkNode);
            DPqPtr = solver.GetMatrixElement(_drainNodePrime, _qNode);
            SPqPtr = solver.GetMatrixElement(_sourceNodePrime, _qNode);
            GqPtr = solver.GetMatrixElement(_gateNode, _qNode);
            BqPtr = solver.GetMatrixElement(_bulkNode, _qNode);
        }

        /// <summary>
        /// Load frequency behavior
        /// </summary>
        /// <param name="simulation">Simulation</param>
        public override void Load(FrequencySimulation simulation)
        {
            var state = simulation.ComplexState;
            double xcggb, xcgdb, xcgsb, xcbgb, xcbdb, xcbsb, xcddb, xcssb, xcdgb;
            double gdpr, gspr, gds, gbd, gbs, capbd, capbs, xcsgb, xcdsb, xcsdb;
            double cggb, cgdb, cgsb, cbgb, cbdb, cbsb, cddb, cdgb, cdsb, omega;
            double GSoverlapCap, GDoverlapCap, GBoverlapCap, FwdSum, RevSum, Gm, Gmbs;
            double dxpart, sxpart, xgtg, xgtd, xgts, xgtb, xcqgb = 0, xcqdb = 0, xcqsb = 0, xcqbb = 0;
            double gbspsp, gbbdp, gbbsp, gbspg, gbspb;
            double gbspdp, gbdpdp, gbdpg, gbdpb, gbdpsp;
            double ddxpart_dVd, ddxpart_dVg, ddxpart_dVb, ddxpart_dVs;
            double dsxpart_dVd, dsxpart_dVg, dsxpart_dVb, dsxpart_dVs;
            double T1, CoxWL, qcheq, Cdg, Cdd, Cds, Csg, Csd, Css;

            // For ACNQSMOD
            double T0, T2, T3, gmr, gmbsr, gmi, gmbsi, gdsi;
            double Cddr, Cdgr, Cdsr, Csdr, Csgr, Cssr, Cgdr, Cggr, Cgsr;
            double Cddi, Cdgi, Cdsi, Cdbi, Csdi, Csgi, Cssi, Csbi;
            double Cgdi, Cggi, Cgsi, Cgbi, Gmi, Gmbsi, FwdSumi, RevSumi;
            double xcdgbi, xcsgbi, xcddbi, xcdsbi, xcsdbi, xcssbi, xcdbbi;
            double xcsbbi, xcggbi, xcgdbi, xcgsbi, xcgbbi;

            var pParam = _temp.Param;
            omega = state.Laplace.Imaginary;
            Csd = -(_load.Cddb + _load.Cgdb + _load.Cbdb);
            Csg = -(_load.Cdgb + _load.Cggb + _load.Cbgb);
            Css = -(_load.Cdsb + _load.Cgsb + _load.Cbsb);

            if (_bp.AcnqsMod > 0)
            {
                T0 = omega * _load.Taunet;
                T1 = T0 * T0;
                T2 = 1.0 / (1.0 + T1);
                T3 = T0 * T2;

                gmr = _load.Gm * T2;
                gmbsr = _load.Gmbs * T2;
                gds = _load.Gds * T2;

                gmi = -_load.Gm * T3;
                gmbsi = -_load.Gmbs * T3;
                gdsi = -_load.Gds * T3;

                Cddr = _load.Cddb * T2;
                Cdgr = _load.Cdgb * T2;
                Cdsr = _load.Cdsb * T2;

                Cddi = _load.Cddb * T3 * omega;
                Cdgi = _load.Cdgb * T3 * omega;
                Cdsi = _load.Cdsb * T3 * omega;
                Cdbi = -(Cddi + Cdgi + Cdsi);

                Csdr = Csd * T2;
                Csgr = Csg * T2;
                Cssr = Css * T2;

                Csdi = Csd * T3 * omega;
                Csgi = Csg * T3 * omega;
                Cssi = Css * T3 * omega;
                Csbi = -(Csdi + Csgi + Cssi);

                Cgdr = -(Cddr + Csdr + _load.Cbdb);
                Cggr = -(Cdgr + Csgr + _load.Cbgb);
                Cgsr = -(Cdsr + Cssr + _load.Cbsb);

                Cgdi = -(Cddi + Csdi);
                Cggi = -(Cdgi + Csgi);
                Cgsi = -(Cdsi + Cssi);
                Cgbi = -(Cgdi + Cggi + Cgsi);
            }
            else /* QS */
            {
                gmr = _load.Gm;
                gmbsr = _load.Gmbs;
                gds = _load.Gds;
                gmi = gmbsi = gdsi = 0.0;

                Cddr = _load.Cddb;
                Cdgr = _load.Cdgb;
                Cdsr = _load.Cdsb;
                Cddi = Cdgi = Cdsi = Cdbi = 0.0;

                Csdr = Csd;
                Csgr = Csg;
                Cssr = Css;
                Csdi = Csgi = Cssi = Csbi = 0.0;

                Cgdr = _load.Cgdb;
                Cggr = _load.Cggb;
                Cgsr = _load.Cgsb;
                Cgdi = Cggi = Cgsi = Cgbi = 0.0;
            }

            if (_load.Mode >= 0)
            {
                Gm = gmr;
                Gmbs = gmbsr;
                FwdSum = Gm + Gmbs;
                RevSum = 0.0;
                Gmi = gmi;
                Gmbsi = gmbsi;
                FwdSumi = Gmi + Gmbsi;
                RevSumi = 0.0;

                gbbdp = -_load.Gbds;
                gbbsp = _load.Gbds + _load.Gbgs + _load.Gbbs;

                gbdpg = _load.Gbgs;
                gbdpb = _load.Gbbs;
                gbdpdp = _load.Gbds;
                gbdpsp = -(gbdpg + gbdpb + gbdpdp);

                gbspdp = 0.0;
                gbspg = 0.0;
                gbspb = 0.0;
                gbspsp = 0.0;

                if (_bp.NqsMod == 0 || _bp.AcnqsMod == 1)
                {
                    cggb = Cggr;
                    cgsb = Cgsr;
                    cgdb = Cgdr;

                    cbgb = _load.Cbgb;
                    cbsb = _load.Cbsb;
                    cbdb = _load.Cbdb;

                    cdgb = Cdgr;
                    cdsb = Cdsr;
                    cddb = Cddr;

                    xgtg = xgtd = xgts = xgtb = 0.0;
                    sxpart = 0.6;
                    dxpart = 0.4;
                    ddxpart_dVd = ddxpart_dVg = ddxpart_dVb = ddxpart_dVs = 0.0;
                    dsxpart_dVd = dsxpart_dVg = dsxpart_dVb = dsxpart_dVs = 0.0;
                }
                else
                {
                    cggb = cgdb = cgsb = 0.0;
                    cbgb = cbdb = cbsb = 0.0;
                    cdgb = cddb = cdsb = 0.0;

                    xgtg = _load.Gtg;
                    xgtd = _load.Gtd;
                    xgts = _load.Gts;
                    xgtb = _load.Gtb;

                    xcqgb = _load.Cqgb * omega;
                    xcqdb = _load.Cqdb * omega;
                    xcqsb = _load.Cqsb * omega;
                    xcqbb = _load.Cqbb * omega;

                    CoxWL = _mbp.Cox * pParam.BSIM3weffCV * pParam.BSIM3leffCV;
                    qcheq = -(_load.Qgate + _load.Qbulk);
                    if (Math.Abs(qcheq) <= 1.0e-5 * CoxWL)
                    {
                        if (_mbp.Xpart < 0.5)
                        {
                            dxpart = 0.4;
                        }
                        else if (_mbp.Xpart > 0.5)
                        {
                            dxpart = 0.0;
                        }
                        else
                        {
                            dxpart = 0.5;
                        }

                        ddxpart_dVd = ddxpart_dVg = ddxpart_dVb = ddxpart_dVs = 0.0;
                    }
                    else
                    {
                        dxpart = _load.Qdrn / qcheq;
                        Cdd = _load.Cddb;
                        Csd = -(_load.Cgdb + _load.Cddb + _load.Cbdb);
                        ddxpart_dVd = (Cdd - dxpart * (Cdd + Csd)) / qcheq;
                        Cdg = _load.Cdgb;
                        Csg = -(_load.Cggb + _load.Cdgb + _load.Cbgb);
                        ddxpart_dVg = (Cdg - dxpart * (Cdg + Csg)) / qcheq;

                        Cds = _load.Cdsb;
                        Css = -(_load.Cgsb + _load.Cdsb + _load.Cbsb);
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

                gbbsp = -_load.Gbds;
                gbbdp = _load.Gbds + _load.Gbgs + _load.Gbbs;

                gbdpg = 0.0;
                gbdpsp = 0.0;
                gbdpb = 0.0;
                gbdpdp = 0.0;

                gbspg = _load.Gbgs;
                gbspsp = _load.Gbds;
                gbspb = _load.Gbbs;
                gbspdp = -(gbspg + gbspsp + gbspb);

                if (_bp.NqsMod == 0 || _bp.AcnqsMod == 1)
                {
                    cggb = Cggr;
                    cgsb = Cgdr;
                    cgdb = Cgsr;

                    cbgb = _load.Cbgb;
                    cbsb = _load.Cbdb;
                    cbdb = _load.Cbsb;

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

                    xgtg = _load.Gtg;
                    xgtd = _load.Gts;
                    xgts = _load.Gtd;
                    xgtb = _load.Gtb;

                    xcqgb = _load.Cqgb * omega;
                    xcqdb = _load.Cqsb * omega;
                    xcqsb = _load.Cqdb * omega;
                    xcqbb = _load.Cqbb * omega;

                    CoxWL = _mbp.Cox * pParam.BSIM3weffCV * pParam.BSIM3leffCV;
                    qcheq = -(_load.Qgate + _load.Qbulk);
                    if (Math.Abs(qcheq) <= 1.0e-5 * CoxWL)
                    {
                        if (_mbp.Xpart < 0.5)
                        {
                            sxpart = 0.4;
                        }
                        else if (_mbp.Xpart > 0.5)
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
                        sxpart = _load.Qdrn / qcheq;
                        Css = _load.Cddb;
                        Cds = -(_load.Cgdb + _load.Cddb + _load.Cbdb);
                        dsxpart_dVs = (Css - sxpart * (Css + Cds)) / qcheq;
                        Csg = _load.Cdgb;
                        Cdg = -(_load.Cggb + _load.Cdgb + _load.Cbgb);
                        dsxpart_dVg = (Csg - sxpart * (Csg + Cdg)) / qcheq;

                        Csd = _load.Cdsb;
                        Cdd = -(_load.Cgsb + _load.Cdsb + _load.Cbsb);
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

            T1 = _load.Qdef * _load.Gtau;
            gdpr = _temp.DrainConductance;
            gspr = _temp.SourceConductance;
            gbd = _load.Gbd;
            gbs = _load.Gbs;
            capbd = _load.Capbd;
            capbs = _load.Capbs;

            GSoverlapCap = _temp.Cgso;
            GDoverlapCap = _temp.Cgdo;
            GBoverlapCap = pParam.BSIM3cgbo;

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

            var m = _bp.Multiplier;
            GgPtr.Value += new Complex(0.0, m * (xcggb));
            BbPtr.Value -= new Complex(0.0, m * (xcbgb + xcbdb + xcbsb));
            DPdpPtr.Value += new Complex(0.0, m * (xcddb + gdsi + RevSumi));
            SPspPtr.Value += new Complex(0.0, m * (xcssb + gdsi + FwdSumi));
            GbPtr.Value -= new Complex(0.0, m * (xcggb + xcgdb + xcgsb));
            GdpPtr.Value += new Complex(0.0, m * (xcgdb));
            GspPtr.Value += new Complex(0.0, m * (xcgsb));
            BgPtr.Value += new Complex(0.0, m * (xcbgb));
            BdpPtr.Value += new Complex(0.0, m * (xcbdb));
            BspPtr.Value += new Complex(0.0, m * (xcbsb));
            DPgPtr.Value += new Complex(0.0, m * (xcdgb + Gmi));
            DPbPtr.Value -= new Complex(0.0, m * (xcdgb + xcddb + xcdsb + Gmbsi));
            DPspPtr.Value += new Complex(0.0, m * (xcdsb - gdsi - FwdSumi));
            SPgPtr.Value += new Complex(0.0, m * (xcsgb - Gmi));
            SPbPtr.Value -= new Complex(0.0, m * (xcsgb + xcsdb + xcssb - Gmbsi));
            SPdpPtr.Value += new Complex(0.0, m * (xcsdb - gdsi - RevSumi));

            DdPtr.Value += m * (gdpr);
            SsPtr.Value += m * (gspr);
            BbPtr.Value += m * (gbd + gbs - _load.Gbbs);
            DPdpPtr.Value += m * (gdpr + gds + gbd + RevSum + xcddbi + dxpart * xgtd + T1 * ddxpart_dVd + gbdpdp);
            SPspPtr.Value += m * (gspr + gds + gbs + FwdSum + xcssbi + sxpart * xgts + T1 * dsxpart_dVs + gbspsp);

            DdpPtr.Value -= m * (gdpr);
            SspPtr.Value -= m * (gspr);

            BgPtr.Value -= m * (_load.Gbgs);
            BdpPtr.Value -= m * (gbd - gbbdp);
            BspPtr.Value -= m * (gbs - gbbsp);

            DPdPtr.Value -= m * (gdpr);
            DPgPtr.Value += m * (Gm + dxpart * xgtg + T1 * ddxpart_dVg + gbdpg + xcdgbi);
            DPbPtr.Value -= m * (gbd - Gmbs - dxpart * xgtb - T1 * ddxpart_dVb - gbdpb - xcdbbi);
            DPspPtr.Value -= m * (gds + FwdSum - dxpart * xgts - T1 * ddxpart_dVs - gbdpsp - xcdsbi);

            SPgPtr.Value -= m * (Gm - sxpart * xgtg - T1 * dsxpart_dVg - gbspg - xcsgbi);
            SPsPtr.Value -= m * (gspr);
            SPbPtr.Value -= m * (gbs + Gmbs - sxpart * xgtb - T1 * dsxpart_dVb - gbspb - xcsbbi);
            SPdpPtr.Value -= m * (gds + RevSum - sxpart * xgtd - T1 * dsxpart_dVd - gbspdp - xcsdbi);

            GgPtr.Value -= m * (xgtg - xcggbi);
            GbPtr.Value -= m * (xgtb - xcgbbi);
            GdpPtr.Value -= m * (xgtd - xcgdbi);
            GspPtr.Value -= m * (xgts - xcgsbi);

            if (_bp.NqsMod > 0)
            {
                if (_bp.AcnqsMod > 0)
                {
                    QqPtr.Value += m;
                    QgPtr.Value += 0.0;
                    QdpPtr.Value += 0.0;
                    QspPtr.Value += 0.0;
                    QbPtr.Value += 0.0;

                    DPqPtr.Value += 0.0;
                    SPqPtr.Value += 0.0;
                    GqPtr.Value += 0.0;

                }
                else
                {
                    QqPtr.Value += new Complex(0.0, m * (omega * ScalingFactor));
                    QgPtr.Value -= new Complex(0.0, m * (xcqgb));
                    QdpPtr.Value -= new Complex(0.0, m * (xcqdb));
                    QspPtr.Value -= new Complex(0.0, m * (xcqsb));
                    QbPtr.Value -= new Complex(0.0, m * (xcqbb));

                    QqPtr.Value += m * (_load.Gtau);

                    DPqPtr.Value += m * (dxpart * _load.Gtau);
                    SPqPtr.Value += m * (sxpart * _load.Gtau);
                    GqPtr.Value -= m * (_load.Gtau);

                    QgPtr.Value += m * (xgtg);
                    QdpPtr.Value += m * (xgtd);
                    QspPtr.Value += m * (xgts);
                    QbPtr.Value += m * (xgtb);
                }
            }
        }
    }
}
