using System;
using System.Numerics;
using SpiceSharp.Algebra;
using SpiceSharp.Behaviors;
using SpiceSharp.Simulations;

namespace SpiceSharp.Components.BSIM3v24Behaviors
{
    /// <summary>
    /// Frequency behavior for a <see cref="BSIM3v24"/>
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
        private int _drainNode, _gateNode, _sourceNode, _bulkNode, _sourceNodePrime, _drainNodePrime, _qNode;
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

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="name">Name</param>
        public FrequencyBehavior(Identifier name) : base(name)
        {
        }

        /// <summary>
        /// Setup behavior
        /// </summary>
        /// <param name="provider">Provider</param>
        public override void Setup(SetupDataProvider provider)
        {
            if (provider == null)
                throw new ArgumentNullException(nameof(provider));

            // Get parameter sets
            _mbp = provider.GetParameterSet<ModelBaseParameters>("model");
            _bp = provider.GetParameterSet<BaseParameters>("entity");

            // Get behaviors
            _load = provider.GetBehavior<LoadBehavior>("entity");
            _temp = provider.GetBehavior<TemperatureBehavior>("entity");
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
            _drainNodePrime = _load.DrainNodePrime;
            _sourceNodePrime = _load.SourceNodePrime;
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
        /// <param name="simulation"></param>
        public override void Load(FrequencySimulation simulation)
        {
            double xcggb, xcgdb, xcgsb, xcbgb, xcbdb, xcbsb, xcddb, xcssb, xcdgb;
            double gdpr, gspr, gds, gbd, gbs, capbd, capbs, xcsgb, xcdsb, xcsdb;
            double cggb, cgdb, cgsb, cbgb, cbdb, cbsb, cddb, cdgb, cdsb, omega;
            double gSoverlapCap, gDoverlapCap, gBoverlapCap, fwdSum, revSum, gm, gmbs;
            double dxpart, sxpart, xgtg, xgtd, xgts, xgtb, xcqgb = 0.0, xcqdb = 0.0, xcqsb = 0.0, xcqbb = 0.0;
            double gbspsp, gbbdp, gbbsp, gbspg, gbspb;
            double gbspdp, gbdpdp, gbdpg, gbdpb, gbdpsp;
            double ddxpartDVd, ddxpartDVg, ddxpartDVb, ddxpartDVs;
            double dsxpartDVd, dsxpartDVg, dsxpartDVb, dsxpartDVs;
            double t1, coxWl, qcheq, cdg, cdd, cds, csg, csd, css;
            var pParam = _temp.Param;
            omega = simulation.ComplexState.Laplace.Imaginary;

            if (_load.Mode >= 0)
            {
                gm = _load.Gm;
                gmbs = _load.Gmbs;
                fwdSum = gm + gmbs;
                revSum = 0.0;

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

                if (_bp.NqsMod == 0)
                {
                    cggb = _load.Cggb;
                    cgsb = _load.Cgsb;
                    cgdb = _load.Cgdb;

                    cbgb = _load.Cbgb;
                    cbsb = _load.Cbsb;
                    cbdb = _load.Cbdb;

                    cdgb = _load.Cdgb;
                    cdsb = _load.Cdsb;
                    cddb = _load.Cddb;

                    xgtg = xgtd = xgts = xgtb = 0.0;
                    sxpart = 0.6;
                    dxpart = 0.4;
                    ddxpartDVd = ddxpartDVg = ddxpartDVb = ddxpartDVs = 0.0;
                    dsxpartDVd = dsxpartDVg = dsxpartDVb = dsxpartDVs = 0.0;
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

                    coxWl = _mbp.Cox * pParam.BSIM3weffCV * pParam.BSIM3leffCV;
                    qcheq = -(_load.Qgate + _load.Qbulk);
                    if (Math.Abs(qcheq) <= 1.0e-5 * coxWl)
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

                        ddxpartDVd = ddxpartDVg = ddxpartDVb = ddxpartDVs = 0.0;
                    }
                    else
                    {
                        dxpart = _load.Qdrn / qcheq;
                        cdd = _load.Cddb;
                        csd = -(_load.Cgdb + _load.Cddb + _load.Cbdb);
                        ddxpartDVd = (cdd - dxpart * (cdd + csd)) / qcheq;
                        cdg = _load.Cdgb;
                        csg = -(_load.Cggb + _load.Cdgb + _load.Cbgb);
                        ddxpartDVg = (cdg - dxpart * (cdg + csg)) / qcheq;

                        cds = _load.Cdsb;
                        css = -(_load.Cgsb + _load.Cdsb + _load.Cbsb);
                        ddxpartDVs = (cds - dxpart * (cds + css)) / qcheq;

                        ddxpartDVb = -(ddxpartDVd + ddxpartDVg
                                                  + ddxpartDVs);
                    }

                    sxpart = 1.0 - dxpart;
                    dsxpartDVd = -ddxpartDVd;
                    dsxpartDVg = -ddxpartDVg;
                    dsxpartDVs = -ddxpartDVs;
                    dsxpartDVb = -(dsxpartDVd + dsxpartDVg + dsxpartDVs);
                }
            }
            else
            {
                gm = -_load.Gm;
                gmbs = -_load.Gmbs;
                fwdSum = 0.0;
                revSum = -(gm + gmbs);

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

                if (_bp.NqsMod == 0)
                {
                    cggb = _load.Cggb;
                    cgsb = _load.Cgdb;
                    cgdb = _load.Cgsb;

                    cbgb = _load.Cbgb;
                    cbsb = _load.Cbdb;
                    cbdb = _load.Cbsb;

                    cdgb = -(_load.Cdgb + cggb + cbgb);
                    cdsb = -(_load.Cddb + cgsb + cbsb);
                    cddb = -(_load.Cdsb + cgdb + cbdb);

                    xgtg = xgtd = xgts = xgtb = 0.0;
                    sxpart = 0.4;
                    dxpart = 0.6;
                    ddxpartDVd = ddxpartDVg = ddxpartDVb = ddxpartDVs = 0.0;
                    dsxpartDVd = dsxpartDVg = dsxpartDVb = dsxpartDVs = 0.0;
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

                    coxWl = _mbp.Cox * pParam.BSIM3weffCV * pParam.BSIM3leffCV;
                    qcheq = -(_load.Qgate + _load.Qbulk);
                    if (Math.Abs(qcheq) <= 1.0e-5 * coxWl)
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

                        dsxpartDVd = dsxpartDVg = dsxpartDVb = dsxpartDVs = 0.0;
                    }
                    else
                    {
                        sxpart = _load.Qdrn / qcheq;
                        css = _load.Cddb;
                        cds = -(_load.Cgdb + _load.Cddb + _load.Cbdb);
                        dsxpartDVs = (css - sxpart * (css + cds)) / qcheq;
                        csg = _load.Cdgb;
                        cdg = -(_load.Cggb + _load.Cdgb + _load.Cbgb);
                        dsxpartDVg = (csg - sxpart * (csg + cdg)) / qcheq;

                        csd = _load.Cdsb;
                        cdd = -(_load.Cgsb + _load.Cdsb + _load.Cbsb);
                        dsxpartDVd = (csd - sxpart * (csd + cdd)) / qcheq;

                        dsxpartDVb = -(dsxpartDVd + dsxpartDVg + dsxpartDVs);
                    }

                    dxpart = 1.0 - sxpart;
                    ddxpartDVd = -dsxpartDVd;
                    ddxpartDVg = -dsxpartDVg;
                    ddxpartDVs = -dsxpartDVs;
                    ddxpartDVb = -(ddxpartDVd + ddxpartDVg + ddxpartDVs);
                }
            }

            t1 = _load.Qdef * _load.Gtau;
            gdpr = _temp.DrainConductance;
            gspr = _temp.SourceConductance;
            gds = _load.Gds;
            gbd = _load.Gbd;
            gbs = _load.Gbs;
            capbd = _load.Capbd;
            capbs = _load.Capbs;

            gSoverlapCap = _temp.Cgso;
            gDoverlapCap = _temp.Cgdo;
            gBoverlapCap = pParam.BSIM3cgbo;

            xcdgb = (cdgb - gDoverlapCap) * omega;
            xcddb = (cddb + capbd + gDoverlapCap) * omega;
            xcdsb = cdsb * omega;
            xcsgb = -(cggb + cbgb + cdgb + gSoverlapCap) * omega;
            xcsdb = -(cgdb + cbdb + cddb) * omega;
            xcssb = (capbs + gSoverlapCap - (cgsb + cbsb + cdsb)) * omega;
            xcggb = (cggb + gDoverlapCap + gSoverlapCap + gBoverlapCap) * omega;
            xcgdb = (cgdb - gDoverlapCap) * omega;
            xcgsb = (cgsb - gSoverlapCap) * omega;
            xcbgb = (cbgb - gBoverlapCap) * omega;
            xcbdb = (cbdb - capbd) * omega;
            xcbsb = (cbsb - capbs) * omega;

            GgPtr.Value += new Complex(0.0, xcggb);
            BbPtr.Value -= new Complex(0.0, xcbgb + xcbdb + xcbsb);
            DPdpPtr.Value += new Complex(0.0, xcddb);
            SPspPtr.Value += new Complex(0.0, xcssb);
            GbPtr.Value -= new Complex(0.0, xcggb + xcgdb + xcgsb);
            GdpPtr.Value += new Complex(0.0, xcgdb);
            GspPtr.Value += new Complex(0.0, xcgsb);
            BgPtr.Value += new Complex(0.0, xcbgb);
            BdpPtr.Value += new Complex(0.0, xcbdb);
            BspPtr.Value += new Complex(0.0, xcbsb);
            DPgPtr.Value += new Complex(0.0, xcdgb);
            DPbPtr.Value -= new Complex(0.0, xcdgb + xcddb + xcdsb);
            DPspPtr.Value += new Complex(0.0, xcdsb);
            SPgPtr.Value += new Complex(0.0, xcsgb);
            SPbPtr.Value -= new Complex(0.0, xcsgb + xcsdb + xcssb);
            SPdpPtr.Value += new Complex(0.0, xcsdb);

            DdPtr.Value += gdpr;
            SsPtr.Value += gspr;
            BbPtr.Value += gbd + gbs - _load.Gbbs;
            DPdpPtr.Value += gdpr + gds + gbd + revSum + dxpart * xgtd + t1 * ddxpartDVd + gbdpdp;
            SPspPtr.Value += gspr + gds + gbs + fwdSum + sxpart * xgts + t1 * dsxpartDVs + gbspsp;

            DdpPtr.Value -= gdpr;
            SspPtr.Value -= gspr;

            BgPtr.Value -= _load.Gbgs;
            BdpPtr.Value -= gbd - gbbdp;
            BspPtr.Value -= gbs - gbbsp;

            DPdPtr.Value -= gdpr;
            DPgPtr.Value += gm + dxpart * xgtg + t1 * ddxpartDVg + gbdpg;
            DPbPtr.Value -= gbd - gmbs - dxpart * xgtb - t1 * ddxpartDVb - gbdpb;
            DPspPtr.Value -= gds + fwdSum - dxpart * xgts - t1 * ddxpartDVs - gbdpsp;

            SPgPtr.Value -= gm - sxpart * xgtg - t1 * dsxpartDVg - gbspg;
            SPsPtr.Value -= gspr;
            SPbPtr.Value -= gbs + gmbs - sxpart * xgtb - t1 * dsxpartDVb - gbspb;
            SPdpPtr.Value -= gds + revSum - sxpart * xgtd - t1 * dsxpartDVd - gbspdp;

            GgPtr.Value -= xgtg;
            GbPtr.Value -= xgtb;
            GdpPtr.Value -= xgtd;
            GspPtr.Value -= xgts;

            if (_bp.NqsMod > 0)
            {
                QqPtr.Value += new Complex(0.0, omega * ScalingFactor);
                QgPtr.Value -= new Complex(0.0, xcqgb);
                QdpPtr.Value -= new Complex(0.0, xcqdb);
                QspPtr.Value -= new Complex(0.0, xcqsb);
                QbPtr.Value -= new Complex(0.0, xcqbb);

                QqPtr.Value += _load.Gtau;

                DPqPtr.Value += dxpart * _load.Gtau;
                SPqPtr.Value += sxpart * _load.Gtau;
                GqPtr.Value -= _load.Gtau;

                QgPtr.Value += xgtg;
                QdpPtr.Value += xgtd;
                QspPtr.Value += xgts;
                QbPtr.Value += xgtb;
            }
        }
    }
}
