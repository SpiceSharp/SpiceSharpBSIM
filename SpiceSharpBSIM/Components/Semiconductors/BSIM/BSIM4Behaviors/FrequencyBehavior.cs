using SpiceSharp;
using SpiceSharp.Algebra;
using SpiceSharp.Attributes;
using SpiceSharp.Behaviors;
using SpiceSharp.Components;
using SpiceSharp.Simulations;
using System.Numerics;

namespace SpiceSharpBSIM.Components.Semiconductors.BSIM.BSIM4Behaviors
{
    /// <summary>
    /// A frequency behavior for a <see cref="BSIM4"/>.
    /// </summary>
    [BehaviorFor(typeof(BSIM4)), AddBehaviorIfNo(typeof(IFrequencyBehavior))]
    public class FrequencyBehavior : BiasingBehavior, IFrequencyBehavior
    {
        private readonly IComplexSimulationState _state;
        protected readonly IVariable<Complex> _drain, _gate, _source, _bulk, _drainPrime, _gatePrime, _gateMid, _sourcePrime,
            _bulkPrime, _drainBulk, _sourceBulk, _q;
        private readonly Element<Complex> _dpbpPtr, _gpbpPtr, _spbpPtr, _bpdpPtr, _bpgpPtr, _bpspPtr, _bpbpPtr, _ddPtr,
            _gpgpPtr, _ssPtr, _dpdpPtr, _spspPtr, _ddpPtr, _gpdpPtr, _gpspPtr, _sspPtr, _dpspPtr, _dpdPtr,
            _dpgpPtr, _spgpPtr, _spsPtr, _spdpPtr, _qqPtr, _qbpPtr, _qdpPtr, _qspPtr, _qgpPtr, _dpqPtr, _spqPtr,
            _gpqPtr, _gegePtr, _gegpPtr, _gpgePtr, _gedpPtr, _gespPtr, _gebpPtr, _gmdpPtr, _gmgpPtr, _gmgmPtr,
            _gmgePtr, _gmspPtr, _gmbpPtr, _dpgmPtr, _gpgmPtr, _gegmPtr, _spgmPtr, _bpgmPtr, _dpdbPtr, _spsbPtr,
            _dbdpPtr, _dbdbPtr, _dbbpPtr, _dbbPtr, _bpdbPtr, _bpbPtr, _bpsbPtr, _sbspPtr, _sbbpPtr, _sbbPtr,
            _sbsbPtr, _bdbPtr, _bbpPtr, _bsbPtr, _bbPtr, _dgpPtr, _dspPtr, _dbpPtr, _sdpPtr, _sgpPtr, _sbpPtr;

        /// <summary>
        /// Creates a new <see cref="FrequencyBehavior"/>.
        /// </summary>
        /// <param name="context">The binding context</param>
        public FrequencyBehavior(ComponentBindingContext context)
            : base(context)
        {
            _state = context.GetState<IComplexSimulationState>();
            _drain = _state.GetSharedVariable(context.Nodes[0]);
            _gate = _state.GetSharedVariable(context.Nodes[1]);
            _source = _state.GetSharedVariable(context.Nodes[2]);
            _bulk = _state.GetSharedVariable(context.Nodes[3]);

            if (this._drainConductance > 0)
                _drainPrime = _state.CreatePrivateVariable(Name.Combine("drain"), Units.Volt);
            else
                _drainPrime = _drain;

            if (this._sourceConductance > 0)
                _sourcePrime = _state.CreatePrivateVariable(Name.Combine("source"), Units.Volt);
            else
                _sourcePrime = _source;

            if (Parameters.RgateMod.Value > 0)
                _gatePrime = _state.CreatePrivateVariable(Name.Combine("gate"), Units.Volt);
            else
                _gatePrime = _gate;

            if (Parameters.RgateMod.Value == 3)
                _gateMid = _state.CreatePrivateVariable(Name.Combine("gatemid"), Units.Volt);
            else
                _gateMid = _gate;

            if (Parameters.RbodyMod.Value == 1 || Parameters.RbodyMod.Value == 2)
            {
                _drainBulk = _state.CreatePrivateVariable(Name.Combine("dbody"), Units.Volt);
                _bulkPrime = _state.CreatePrivateVariable(Name.Combine("body"), Units.Volt);
                _sourceBulk = _state.CreatePrivateVariable(Name.Combine("sbody"), Units.Volt);
            }
            else
            {
                _drainBulk = _bulk;
                _bulkPrime = _bulk;
                _sourceBulk = _bulk;
            }

            if (Parameters.TrnqsMod.Value != 0)
                _q = _state.CreatePrivateVariable(Name.Combine("charge"), Units.Coulomb);
            else
                _q = _state.GetSharedVariable(Constants.Ground);

            int drain = _state.Map[_drain];
            int gate = _state.Map[_gate];
            int source = _state.Map[_source];
            int bulk = _state.Map[_bulk];
            int drainPrime = _state.Map[_drainPrime];
            int gatePrime = _state.Map[_gatePrime];
            int gateMid = _state.Map[_gateMid];
            int sourcePrime = _state.Map[_sourcePrime];
            int bulkPrime = _state.Map[_bulkPrime];
            int q = _state.Map[_q];
            int dbulk = _state.Map[_drainBulk];
            int sbulk = _state.Map[_sourceBulk];

            _dpbpPtr = _state.Solver.GetElement(new MatrixLocation(drainPrime, bulkPrime));
            _gpbpPtr = _state.Solver.GetElement(new MatrixLocation(gatePrime, bulkPrime));
            _spbpPtr = _state.Solver.GetElement(new MatrixLocation(sourcePrime, bulkPrime));
            _bpdpPtr = _state.Solver.GetElement(new MatrixLocation(bulkPrime, drainPrime));
            _bpgpPtr = _state.Solver.GetElement(new MatrixLocation(bulkPrime, gatePrime));
            _bpspPtr = _state.Solver.GetElement(new MatrixLocation(bulkPrime, sourcePrime));
            _bpbpPtr = _state.Solver.GetElement(new MatrixLocation(bulkPrime, bulkPrime));
            _ddPtr = _state.Solver.GetElement(new MatrixLocation(drain, drain));
            _gpgpPtr = _state.Solver.GetElement(new MatrixLocation(gatePrime, gatePrime));
            _ssPtr = _state.Solver.GetElement(new MatrixLocation(source, source));
            _dpdpPtr = _state.Solver.GetElement(new MatrixLocation(drainPrime, drainPrime));
            _spspPtr = _state.Solver.GetElement(new MatrixLocation(sourcePrime, sourcePrime));
            _ddpPtr = _state.Solver.GetElement(new MatrixLocation(drain, drainPrime));
            _gpdpPtr = _state.Solver.GetElement(new MatrixLocation(gatePrime, drainPrime));
            _gpspPtr = _state.Solver.GetElement(new MatrixLocation(gatePrime, sourcePrime));
            _sspPtr = _state.Solver.GetElement(new MatrixLocation(source, sourcePrime));
            _dpspPtr = _state.Solver.GetElement(new MatrixLocation(drainPrime, sourcePrime));
            _dpdPtr = _state.Solver.GetElement(new MatrixLocation(drainPrime, drain));
            _dpgpPtr = _state.Solver.GetElement(new MatrixLocation(drainPrime, gatePrime));
            _spgpPtr = _state.Solver.GetElement(new MatrixLocation(sourcePrime, gatePrime));
            _spsPtr = _state.Solver.GetElement(new MatrixLocation(sourcePrime, source));
            _spdpPtr = _state.Solver.GetElement(new MatrixLocation(sourcePrime, drainPrime));
            _qqPtr = _state.Solver.GetElement(new MatrixLocation(q, q));
            _qbpPtr = _state.Solver.GetElement(new MatrixLocation(q, bulkPrime));
            _qdpPtr = _state.Solver.GetElement(new MatrixLocation(q, drainPrime));
            _qspPtr = _state.Solver.GetElement(new MatrixLocation(q, sourcePrime));
            _qgpPtr = _state.Solver.GetElement(new MatrixLocation(q, gatePrime));
            _dpqPtr = _state.Solver.GetElement(new MatrixLocation(drainPrime, q));
            _spqPtr = _state.Solver.GetElement(new MatrixLocation(sourcePrime, q));
            _gpqPtr = _state.Solver.GetElement(new MatrixLocation(gatePrime, q));
            if (Parameters.RgateMod.Value != 0)
            {
                _gegePtr = _state.Solver.GetElement(new MatrixLocation(gate, gate));
                _gegpPtr = _state.Solver.GetElement(new MatrixLocation(gate, gatePrime));
                _gpgePtr = _state.Solver.GetElement(new MatrixLocation(gatePrime, gate));
                _gedpPtr = _state.Solver.GetElement(new MatrixLocation(gate, drainPrime));
                _gespPtr = _state.Solver.GetElement(new MatrixLocation(gate, sourcePrime));
                _gebpPtr = _state.Solver.GetElement(new MatrixLocation(gate, bulkPrime));
                _gmdpPtr = _state.Solver.GetElement(new MatrixLocation(gateMid, drainPrime));
                _gmgpPtr = _state.Solver.GetElement(new MatrixLocation(gateMid, gatePrime));
                _gmgmPtr = _state.Solver.GetElement(new MatrixLocation(gateMid, gateMid));
                _gmgePtr = _state.Solver.GetElement(new MatrixLocation(gateMid, gate));
                _gmspPtr = _state.Solver.GetElement(new MatrixLocation(gateMid, sourcePrime));
                _gmbpPtr = _state.Solver.GetElement(new MatrixLocation(gateMid, bulkPrime));
                _dpgmPtr = _state.Solver.GetElement(new MatrixLocation(drainPrime, gateMid));
                _gpgmPtr = _state.Solver.GetElement(new MatrixLocation(gatePrime, gateMid));
                _gegmPtr = _state.Solver.GetElement(new MatrixLocation(gate, gateMid));
                _spgmPtr = _state.Solver.GetElement(new MatrixLocation(sourcePrime, gateMid));
                _bpgmPtr = _state.Solver.GetElement(new MatrixLocation(bulkPrime, gateMid));
            }
            if ((Parameters.RbodyMod.Value == 1) || (Parameters.RbodyMod.Value == 2))
            {
                _dpdbPtr = _state.Solver.GetElement(new MatrixLocation(drainPrime, dbulk));
                _spsbPtr = _state.Solver.GetElement(new MatrixLocation(sourcePrime, sbulk));
                _dbdpPtr = _state.Solver.GetElement(new MatrixLocation(dbulk, drainPrime));
                _dbdbPtr = _state.Solver.GetElement(new MatrixLocation(dbulk, dbulk));
                _dbbpPtr = _state.Solver.GetElement(new MatrixLocation(dbulk, bulkPrime));
                _dbbPtr = _state.Solver.GetElement(new MatrixLocation(dbulk, bulk));
                _bpdbPtr = _state.Solver.GetElement(new MatrixLocation(bulkPrime, dbulk));
                _bpbPtr = _state.Solver.GetElement(new MatrixLocation(bulkPrime, bulk));
                _bpsbPtr = _state.Solver.GetElement(new MatrixLocation(bulkPrime, sbulk));
                _sbspPtr = _state.Solver.GetElement(new MatrixLocation(sbulk, sourcePrime));
                _sbbpPtr = _state.Solver.GetElement(new MatrixLocation(sbulk, bulkPrime));
                _sbbPtr = _state.Solver.GetElement(new MatrixLocation(sbulk, bulk));
                _sbsbPtr = _state.Solver.GetElement(new MatrixLocation(sbulk, sbulk));
                _bdbPtr = _state.Solver.GetElement(new MatrixLocation(bulk, dbulk));
                _bbpPtr = _state.Solver.GetElement(new MatrixLocation(bulk, bulkPrime));
                _bsbPtr = _state.Solver.GetElement(new MatrixLocation(bulk, sbulk));
                _bbPtr = _state.Solver.GetElement(new MatrixLocation(bulk, bulk));
            }
            if (ModelParameters.RdsMod.Value != 0)
            {
                _dgpPtr = _state.Solver.GetElement(new MatrixLocation(drain, gatePrime));
                _dspPtr = _state.Solver.GetElement(new MatrixLocation(drain, sourcePrime));
                _dbpPtr = _state.Solver.GetElement(new MatrixLocation(drain, bulkPrime));
                _sdpPtr = _state.Solver.GetElement(new MatrixLocation(source, drainPrime));
                _sgpPtr = _state.Solver.GetElement(new MatrixLocation(source, gatePrime));
                _sbpPtr = _state.Solver.GetElement(new MatrixLocation(source, bulkPrime));
            }
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
            double gjbd, gjbs, geltd, gcrg, gcrgg, gcrgd, gcrgs, gcrgb;
            double xcbgb, xcbdb, xcbsb, xcbbb;
            double xcggbr, xcgdbr, xcgsbr, xcgbbr, xcggbi, xcgdbi, xcgsbi, xcgbbi;
            double Cggr, Cgdr, Cgsr, Cgbr, Cggi, Cgdi, Cgsi, Cgbi;
            double xcddbr, xcdgbr, xcdsbr, xcdbbr, xcsdbr, xcsgbr, xcssbr, xcsbbr;
            double xcddbi, xcdgbi, xcdsbi, xcdbbi, xcsdbi, xcsgbi, xcssbi, xcsbbi;
            double xcdbdb, xcsbsb = 0.0, xcgmgmb = 0.0, xcgmdb = 0.0, xcgmsb = 0.0, xcdgmb, xcsgmb;
            double xcgmbb = 0.0, xcbgmb;
            double capbd, capbs, omega;
            double gstot, gstotd, gstotg, gstots, gstotb, gspr;
            double gdtot, gdtotd, gdtotg, gdtots, gdtotb, gdpr;
            double gIstotg, gIstotd, gIstots, gIstotb;
            double gIdtotg, gIdtotd, gIdtots, gIdtotb;
            double gIbtotg, gIbtotd, gIbtots, gIbtotb;
            double gIgtotg, gIgtotd, gIgtots, gIgtotb;
            double cgso, cgdo, cgbo;
            double gbspsp, gbbdp, gbbsp, gbspg, gbspb;
            double gbspdp, gbdpdp, gbdpg, gbdpb, gbdpsp;
            double T0 = 0.0, T1, T2, T3;
            double Csg, Csd, Css;
            double Cdgr, Cddr, Cdsr, Cdbr, Csgr, Csdr, Cssr, Csbr;
            double Cdgi, Cddi, Cdsi, Cdbi, Csgi, Csdi, Cssi, Csbi;
            double gmr, gmi, gmbsr, gmbsi, gdsr, gdsi;
            double FwdSumr, RevSumr, Gmr, Gmbsr;
            double FwdSumi, RevSumi, Gmi, Gmbsi;
            double ggidld, ggidlg, ggidlb, ggislg, ggislb, ggisls;

            double m;

            omega = _state.Laplace.Imaginary;
            capbd = this._capbd;
            capbs = this._capbs;
            cgso = this._cgso;
            cgdo = this._cgdo;
            cgbo = Param.BSIM4cgbo;

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
                gdsr = this._gds * T2;

                gmi = -this._gm * T3;
                gmbsi = -this._gmbs * T3;
                gdsi = -this._gds * T3;

                Cddr = this._cddb * T2;
                Cdgr = this._cdgb * T2;
                Cdsr = this._cdsb * T2;
                Cdbr = -(Cddr + Cdgr + Cdsr);

                /* WDLiu: Cxyi mulitplied by jomega below, and actually to be of conductance */
                Cddi = this._cddb * T3 * omega;
                Cdgi = this._cdgb * T3 * omega;
                Cdsi = this._cdsb * T3 * omega;
                Cdbi = -(Cddi + Cdgi + Cdsi);

                Csdr = Csd * T2;
                Csgr = Csg * T2;
                Cssr = Css * T2;
                Csbr = -(Csdr + Csgr + Cssr);

                Csdi = Csd * T3 * omega;
                Csgi = Csg * T3 * omega;
                Cssi = Css * T3 * omega;
                Csbi = -(Csdi + Csgi + Cssi);

                Cgdr = -(Cddr + Csdr + this._cbdb);
                Cggr = -(Cdgr + Csgr + this._cbgb);
                Cgsr = -(Cdsr + Cssr + this._cbsb);
                Cgbr = -(Cgdr + Cggr + Cgsr);

                Cgdi = -(Cddi + Csdi);
                Cggi = -(Cdgi + Csgi);
                Cgsi = -(Cdsi + Cssi);
                Cgbi = -(Cgdi + Cggi + Cgsi);
            }
            else /* QS */
            {
                gmr = this._gm;
                gmbsr = this._gmbs;
                gdsr = this._gds;
                gmi = gmbsi = gdsi = 0.0;

                Cddr = this._cddb;
                Cdgr = this._cdgb;
                Cdsr = this._cdsb;
                Cdbr = -(Cddr + Cdgr + Cdsr);
                Cddi = Cdgi = Cdsi = Cdbi = 0.0;

                Csdr = Csd;
                Csgr = Csg;
                Cssr = Css;
                Csbr = -(Csdr + Csgr + Cssr);
                Csdi = Csgi = Cssi = Csbi = 0.0;

                Cgdr = this._cgdb;
                Cggr = this._cggb;
                Cgsr = this._cgsb;
                Cgbr = -(Cgdr + Cggr + Cgsr);
                Cgdi = Cggi = Cgsi = Cgbi = 0.0;
            }


            if (this._mode >= 0)
            {
                Gmr = gmr;
                Gmbsr = gmbsr;
                FwdSumr = Gmr + Gmbsr;
                RevSumr = 0.0;
                Gmi = gmi;
                Gmbsi = gmbsi;
                FwdSumi = Gmi + Gmbsi;
                RevSumi = 0.0;

                gbbdp = -(this._gbds);
                gbbsp = this._gbds + this._gbgs + this._gbbs;
                gbdpg = this._gbgs;
                gbdpdp = this._gbds;
                gbdpb = this._gbbs;
                gbdpsp = -(gbdpg + gbdpdp + gbdpb);

                gbspdp = 0.0;
                gbspg = 0.0;
                gbspb = 0.0;
                gbspsp = 0.0;

                if (ModelParameters.IgcMod.Value != 0)
                {
                    gIstotg = this._gIgsg + this._gIgcsg;
                    gIstotd = this._gIgcsd;
                    gIstots = this._gIgss + this._gIgcss;
                    gIstotb = this._gIgcsb;

                    gIdtotg = this._gIgdg + this._gIgcdg;
                    gIdtotd = this._gIgdd + this._gIgcdd;
                    gIdtots = this._gIgcds;
                    gIdtotb = this._gIgcdb;
                }
                else
                {
                    gIstotg = gIstotd = gIstots = gIstotb = 0.0;
                    gIdtotg = gIdtotd = gIdtots = gIdtotb = 0.0;
                }

                if (ModelParameters.IgbMod.Value != 0)
                {
                    gIbtotg = this._gIgbg;
                    gIbtotd = this._gIgbd;
                    gIbtots = this._gIgbs;
                    gIbtotb = this._gIgbb;
                }
                else
                    gIbtotg = gIbtotd = gIbtots = gIbtotb = 0.0;

                if ((ModelParameters.IgcMod != 0) || (ModelParameters.IgbMod != 0))
                {
                    gIgtotg = gIstotg + gIdtotg + gIbtotg;
                    gIgtotd = gIstotd + gIdtotd + gIbtotd;
                    gIgtots = gIstots + gIdtots + gIbtots;
                    gIgtotb = gIstotb + gIdtotb + gIbtotb;
                }
                else
                    gIgtotg = gIgtotd = gIgtots = gIgtotb = 0.0;

                if (Parameters.RgateMod == 2)
                    T0 = this._vges - this._vgs;
                else if (Parameters.RgateMod == 3)
                    T0 = this._vgms - this._vgs;
                if (Parameters.RgateMod > 1)
                {
                    gcrgd = this._gcrgd * T0;
                    gcrgg = this._gcrgg * T0;
                    gcrgs = this._gcrgs * T0;
                    gcrgb = this._gcrgb * T0;
                    gcrgg -= this._gcrg;
                    gcrg = this._gcrg;
                }
                else
                    gcrg = gcrgd = gcrgg = gcrgs = gcrgb = 0.0;

                if (Parameters.RgateMod == 3)
                {
                    xcgmgmb = (cgdo + cgso + Param.BSIM4cgbo) * omega;
                    xcgmdb = -cgdo * omega;
                    xcgmsb = -cgso * omega;
                    xcgmbb = -Param.BSIM4cgbo * omega;

                    xcdgmb = xcgmdb;
                    xcsgmb = xcgmsb;
                    xcbgmb = xcgmbb;

                    xcggbr = Cggr * omega;
                    xcgdbr = Cgdr * omega;
                    xcgsbr = Cgsr * omega;
                    xcgbbr = -(xcggbr + xcgdbr + xcgsbr);

                    xcdgbr = Cdgr * omega;
                    xcsgbr = Csgr * omega;
                    xcbgb = this._cbgb * omega;
                }
                else
                {
                    xcggbr = (Cggr + cgdo + cgso + Param.BSIM4cgbo) * omega;
                    xcgdbr = (Cgdr - cgdo) * omega;
                    xcgsbr = (Cgsr - cgso) * omega;
                    xcgbbr = -(xcggbr + xcgdbr + xcgsbr);

                    xcdgbr = (Cdgr - cgdo) * omega;
                    xcsgbr = (Csgr - cgso) * omega;
                    xcbgb = (this._cbgb - Param.BSIM4cgbo) * omega;

                    xcdgmb = xcsgmb = xcbgmb = 0.0;
                }
                xcddbr = (Cddr + this._capbd + cgdo) * omega;
                xcdsbr = Cdsr * omega;
                xcsdbr = Csdr * omega;
                xcssbr = (this._capbs + cgso + Cssr) * omega;

                if (Parameters.RbodyMod.Value == 0)
                {
                    xcdbbr = -(xcdgbr + xcddbr + xcdsbr + xcdgmb);
                    xcsbbr = -(xcsgbr + xcsdbr + xcssbr + xcsgmb);

                    xcbdb = (this._cbdb - this._capbd) * omega;
                    xcbsb = (this._cbsb - this._capbs) * omega;
                    xcdbdb = 0.0;
                }
                else
                {
                    xcdbbr = Cdbr * omega;
                    xcsbbr = -(xcsgbr + xcsdbr + xcssbr + xcsgmb)
                           + this._capbs * omega;

                    xcbdb = this._cbdb * omega;
                    xcbsb = this._cbsb * omega;

                    xcdbdb = -this._capbd * omega;
                    xcsbsb = -this._capbs * omega;
                }
                xcbbb = -(xcbdb + xcbgb + xcbsb + xcbgmb);

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
            else /* Reverse mode */
            {
                Gmr = -gmr;
                Gmbsr = -gmbsr;
                FwdSumr = 0.0;
                RevSumr = -(Gmr + Gmbsr);
                Gmi = -gmi;
                Gmbsi = -gmbsi;
                FwdSumi = 0.0;
                RevSumi = -(Gmi + Gmbsi);

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

                    gIdtotg = this._gIgdg + this._gIgcsg;
                    gIdtotd = this._gIgdd + this._gIgcss;
                    gIdtots = this._gIgcsd;
                    gIdtotb = this._gIgcsb;
                }
                else
                {
                    gIstotg = gIstotd = gIstots = gIstotb = 0.0;
                    gIdtotg = gIdtotd = gIdtots = gIdtotb = 0.0;
                }

                if (ModelParameters.IgbMod.Value != 0)
                {
                    gIbtotg = this._gIgbg;
                    gIbtotd = this._gIgbs;
                    gIbtots = this._gIgbd;
                    gIbtotb = this._gIgbb;
                }
                else
                    gIbtotg = gIbtotd = gIbtots = gIbtotb = 0.0;

                if ((ModelParameters.IgcMod != 0) || (ModelParameters.IgbMod != 0))
                {
                    gIgtotg = gIstotg + gIdtotg + gIbtotg;
                    gIgtotd = gIstotd + gIdtotd + gIbtotd;
                    gIgtots = gIstots + gIdtots + gIbtots;
                    gIgtotb = gIstotb + gIdtotb + gIbtotb;
                }
                else
                    gIgtotg = gIgtotd = gIgtots = gIgtotb = 0.0;

                if (Parameters.RgateMod == 2)
                    T0 = this._vges - this._vgs;
                else if (Parameters.RgateMod == 3)
                    T0 = this._vgms - this._vgs;
                if (Parameters.RgateMod > 1)
                {
                    gcrgd = this._gcrgs * T0;
                    gcrgg = this._gcrgg * T0;
                    gcrgs = this._gcrgd * T0;
                    gcrgb = this._gcrgb * T0;
                    gcrgg -= this._gcrg;
                    gcrg = this._gcrg;
                }
                else
                    gcrg = gcrgd = gcrgg = gcrgs = gcrgb = 0.0;

                if (Parameters.RgateMod == 3)
                {
                    xcgmgmb = (cgdo + cgso + Param.BSIM4cgbo) * omega;
                    xcgmdb = -cgdo * omega;
                    xcgmsb = -cgso * omega;
                    xcgmbb = -Param.BSIM4cgbo * omega;

                    xcdgmb = xcgmdb;
                    xcsgmb = xcgmsb;
                    xcbgmb = xcgmbb;

                    xcggbr = Cggr * omega;
                    xcgdbr = Cgsr * omega;
                    xcgsbr = Cgdr * omega;
                    xcgbbr = -(xcggbr + xcgdbr + xcgsbr);

                    xcdgbr = Csgr * omega;
                    xcsgbr = Cdgr * omega;
                    xcbgb = this._cbgb * omega;
                }
                else
                {
                    xcggbr = (Cggr + cgdo + cgso + Param.BSIM4cgbo) * omega;
                    xcgdbr = (Cgsr - cgdo) * omega;
                    xcgsbr = (Cgdr - cgso) * omega;
                    xcgbbr = -(xcggbr + xcgdbr + xcgsbr);

                    xcdgbr = (Csgr - cgdo) * omega;
                    xcsgbr = (Cdgr - cgso) * omega;
                    xcbgb = (this._cbgb - Param.BSIM4cgbo) * omega;

                    xcdgmb = xcsgmb = xcbgmb = 0.0;
                }
                xcddbr = (this._capbd + cgdo + Cssr) * omega;
                xcdsbr = Csdr * omega;
                xcsdbr = Cdsr * omega;
                xcssbr = (Cddr + this._capbs + cgso) * omega;

                if (Parameters.RbodyMod.Value == 0)
                {
                    xcdbbr = -(xcdgbr + xcddbr + xcdsbr + xcdgmb);
                    xcsbbr = -(xcsgbr + xcsdbr + xcssbr + xcsgmb);

                    xcbdb = (this._cbsb - this._capbd) * omega;
                    xcbsb = (this._cbdb - this._capbs) * omega;
                    xcdbdb = 0.0;
                }
                else
                {
                    xcdbbr = -(xcdgbr + xcddbr + xcdsbr + xcdgmb)
                           + this._capbd * omega;
                    xcsbbr = Cdbr * omega;

                    xcbdb = this._cbsb * omega;
                    xcbsb = this._cbdb * omega;
                    xcdbdb = -this._capbd * omega;
                    xcsbsb = -this._capbs * omega;
                }
                xcbbb = -(xcbgb + xcbdb + xcbsb + xcbgmb);

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

            if (ModelParameters.RdsMod == 1)
            {
                gstot = this._gstot;
                gstotd = this._gstotd;
                gstotg = this._gstotg;
                gstots = this._gstots - gstot;
                gstotb = this._gstotb;

                gdtot = this._gdtot;
                gdtotd = this._gdtotd - gdtot;
                gdtotg = this._gdtotg;
                gdtots = this._gdtots;
                gdtotb = this._gdtotb;
            }
            else
            {
                gstot = gstotd = gstotg = gstots = gstotb = 0.0;
                gdtot = gdtotd = gdtotg = gdtots = gdtotb = 0.0;
            }


            /*
             * Loading AC matrix
             */
            m = Parameters.M;

            if (ModelParameters.RdsMod.Value == 0)
            {
                gdpr = this._drainConductance;
                gspr = this._sourceConductance;
            }
            else
                gdpr = gspr = 0.0;

            if (Parameters.RbodyMod.Value == 0)
            {
                gjbd = this._gbd;
                gjbs = this._gbs;
            }
            else
                gjbd = gjbs = 0.0;

            geltd = this._grgeltd;

            if (Parameters.RgateMod == 1)
            {
                this._gegePtr.Value += m * geltd;
                this._gpgePtr.Value -= m * geltd;
                this._gegpPtr.Value -= m * geltd;

                this._gpgpPtr.Value += new Complex(0, m * xcggbr);
                this._gpgpPtr.Value += m * (geltd + xcggbi + gIgtotg);
                this._gpdpPtr.Value += new Complex(0, m * xcgdbr);
                this._gpdpPtr.Value += m * (xcgdbi + gIgtotd);
                this._gpspPtr.Value += new Complex(0, m * xcgsbr);
                this._gpspPtr.Value += m * (xcgsbi + gIgtots);
                this._gpbpPtr.Value += new Complex(0, m * xcgbbr);
                this._gpbpPtr.Value += m * (xcgbbi + gIgtotb);
            } /* WDLiu: gcrg already subtracted from all gcrgg below */
            else if (Parameters.RgateMod == 2)
            {
                this._gegePtr.Value += m * gcrg;
                this._gegpPtr.Value += m * gcrgg;
                this._gedpPtr.Value += m * gcrgd;
                this._gespPtr.Value += m * gcrgs;
                this._gebpPtr.Value += m * gcrgb;

                this._gpgePtr.Value -= m * gcrg;
                this._gpgpPtr.Value += new Complex(0, m * xcggbr);
                this._gpgpPtr.Value -= m * (gcrgg - xcggbi - gIgtotg);
                this._gpdpPtr.Value += new Complex(0, m * xcgdbr);
                this._gpdpPtr.Value -= m * (gcrgd - xcgdbi - gIgtotd);
                this._gpspPtr.Value += new Complex(0, m * xcgsbr);
                this._gpspPtr.Value -= m * (gcrgs - xcgsbi - gIgtots);
                this._gpbpPtr.Value += new Complex(0, m * xcgbbr);
                this._gpbpPtr.Value -= m * (gcrgb - xcgbbi - gIgtotb);
            }
            else if (Parameters.RgateMod == 3)
            {
                this._gegePtr.Value += m * geltd;
                this._gegmPtr.Value -= m * geltd;
                this._gmgePtr.Value -= m * geltd;
                this._gmgmPtr.Value += m * (geltd + gcrg);
                this._gmgmPtr.Value += new Complex(0, m * xcgmgmb);

                this._gmdpPtr.Value += m * gcrgd;
                this._gmdpPtr.Value += new Complex(0, m * xcgmdb);
                this._gmgpPtr.Value += m * gcrgg;
                this._gmspPtr.Value += m * gcrgs;
                this._gmspPtr.Value += new Complex(0, m * xcgmsb);
                this._gmbpPtr.Value += m * gcrgb;
                this._gmbpPtr.Value += new Complex(0, m * xcgmbb);

                this._dpgmPtr.Value += new Complex(0, m * xcdgmb);
                this._gpgmPtr.Value -= m * gcrg;
                this._spgmPtr.Value += new Complex(0, m * xcsgmb);
                this._bpgmPtr.Value += new Complex(0, m * xcbgmb);

                this._gpgpPtr.Value -= m * (gcrgg - xcggbi - gIgtotg);
                this._gpgpPtr.Value += new Complex(0, m * xcggbr);
                this._gpdpPtr.Value -= m * (gcrgd - xcgdbi - gIgtotd);
                this._gpdpPtr.Value += new Complex(0, m * xcgdbr);
                this._gpspPtr.Value -= m * (gcrgs - xcgsbi - gIgtots);
                this._gpspPtr.Value += new Complex(0, m * xcgsbr);
                this._gpbpPtr.Value -= m * (gcrgb - xcgbbi - gIgtotb);
                this._gpbpPtr.Value += new Complex(0, m * xcgbbr);
            }
            else
            {
                this._gpgpPtr.Value += new Complex(0, m * xcggbr);
                this._gpgpPtr.Value += m * (xcggbi + gIgtotg);
                this._gpdpPtr.Value += new Complex(0, m * xcgdbr);
                this._gpdpPtr.Value += m * (xcgdbi + gIgtotd);
                this._gpspPtr.Value += new Complex(0, m * xcgsbr);
                this._gpspPtr.Value += m * (xcgsbi + gIgtots);
                this._gpbpPtr.Value += new Complex(0, m * xcgbbr);
                this._gpbpPtr.Value += m * (xcgbbi + gIgtotb);
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

            this._dpdpPtr.Value += new Complex(0, m * (xcddbr + gdsi + RevSumi));
            this._dpdpPtr.Value += m * (gdpr + xcddbi + gdsr + this._gbd
                                   - gdtotd + RevSumr + gbdpdp - gIdtotd);
            this._dpdPtr.Value -= m * (gdpr + gdtot);
            this._dpgpPtr.Value += new Complex(0, m * (xcdgbr + Gmi));
            this._dpgpPtr.Value += m * (Gmr + xcdgbi - gdtotg + gbdpg - gIdtotg);
            this._dpspPtr.Value += new Complex(0, m * (xcdsbr - gdsi - FwdSumi));
            this._dpspPtr.Value -= m * (gdsr - xcdsbi + FwdSumr + gdtots - gbdpsp + gIdtots);
            this._dpbpPtr.Value += new Complex(0, m * (xcdbbr + Gmbsi));
            this._dpbpPtr.Value -= m * (gjbd + gdtotb - xcdbbi - Gmbsr - gbdpb + gIdtotb);

            this._ddpPtr.Value -= m * (gdpr - gdtotd);
            this._ddPtr.Value += m * (gdpr + gdtot);

            this._spdpPtr.Value += new Complex(0, m * (xcsdbr - gdsi - RevSumi));
            this._spdpPtr.Value -= m * (gdsr - xcsdbi + gstotd + RevSumr - gbspdp + gIstotd);
            this._spgpPtr.Value += new Complex(0, m * (xcsgbr - Gmi));
            this._spgpPtr.Value -= m * (Gmr - xcsgbi + gstotg - gbspg + gIstotg);
            this._spspPtr.Value += new Complex(0, m * (xcssbr + gdsi + FwdSumi));
            this._spspPtr.Value += m * (gspr + xcssbi + gdsr + this._gbs
                                   - gstots + FwdSumr + gbspsp - gIstots);
            this._spsPtr.Value -= m * (gspr + gstot);
            this._spbpPtr.Value += new Complex(0, m * (xcsbbr - Gmbsi));
            this._spbpPtr.Value -= m * (gjbs + gstotb - xcsbbi + Gmbsr - gbspb + gIstotb);

            this._sspPtr.Value -= m * (gspr - gstots);
            this._ssPtr.Value += m * (gspr + gstot);

            this._bpdpPtr.Value += new Complex(0, m * xcbdb);
            this._bpdpPtr.Value -= m * (gjbd - gbbdp + gIbtotd);
            this._bpgpPtr.Value += new Complex(0, m * xcbgb);
            this._bpgpPtr.Value -= m * (this._gbgs + gIbtotg);
            this._bpspPtr.Value += new Complex(0, m * xcbsb);
            this._bpspPtr.Value -= m * (gjbs - gbbsp + gIbtots);
            this._bpbpPtr.Value += new Complex(0, m * xcbbb);
            this._bpbpPtr.Value += m * (gjbd + gjbs - this._gbbs
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
            this._dpspPtr.Value -= m * ((ggidlg + ggidld) + ggidlb);
            this._dpbpPtr.Value += m * ggidlb;
            this._bpdpPtr.Value -= m * ggidld;
            this._bpgpPtr.Value -= m * ggidlg;
            this._bpspPtr.Value += m * ((ggidlg + ggidld) + ggidlb);
            this._bpbpPtr.Value -= m * ggidlb;
            /* stamp gisl */
            this._spdpPtr.Value -= m * ((ggisls + ggislg) + ggislb);
            this._spgpPtr.Value += m * ggislg;
            this._spspPtr.Value += m * ggisls;
            this._spbpPtr.Value += m * ggislb;
            this._bpdpPtr.Value += m * ((ggislg + ggisls) + ggislb);
            this._bpgpPtr.Value -= m * ggislg;
            this._bpspPtr.Value -= m * ggisls;
            this._bpbpPtr.Value -= m * ggislb;

            if (Parameters.RbodyMod.Value != 0)
            {
                this._dpdbPtr.Value += new Complex(0, m * xcdbdb);
                this._dpdbPtr.Value -= m * this._gbd;
                this._spsbPtr.Value += new Complex(0, m * xcsbsb);
                this._spsbPtr.Value -= m * this._gbs;

                this._dbdpPtr.Value += new Complex(0, m * xcdbdb);
                this._dbdpPtr.Value -= m * this._gbd;
                this._dbdbPtr.Value -= new Complex(0, m * xcdbdb);
                this._dbdbPtr.Value += m * (this._gbd + this._grbpd
                                        + this._grbdb);
                this._dbbpPtr.Value -= m * this._grbpd;
                this._dbbPtr.Value -= m * this._grbdb;

                this._bpdbPtr.Value -= m * this._grbpd;
                this._bpbPtr.Value -= m * this._grbpb;
                this._bpsbPtr.Value -= m * this._grbps;
                this._bpbpPtr.Value += m * (this._grbpd + this._grbps
                                        + this._grbpb);
                /* WDLiu: (-this._gbbs) already added to BPbpPtr */

                this._sbspPtr.Value += new Complex(0, m * xcsbsb);
                this._sbspPtr.Value -= m * this._gbs;
                this._sbbpPtr.Value -= m * this._grbps;
                this._sbbPtr.Value -= m * this._grbsb;
                this._sbsbPtr.Value -= new Complex(0, m * xcsbsb);
                this._sbsbPtr.Value += m * (this._gbs
                                        + this._grbps + this._grbsb);

                this._bdbPtr.Value -= m * this._grbdb;
                this._bbpPtr.Value -= m * this._grbpb;
                this._bsbPtr.Value -= m * this._grbsb;
                this._bbPtr.Value += m * (this._grbsb + this._grbdb
                                      + this._grbpb);
            }


            /*
             * WDLiu: The internal charge node generated for transient NQS is not needed for
             *        AC NQS. The following is not doing a real job, but we have to keep it;
             *        otherwise a singular AC NQS matrix may occur if the transient NQS is on.
             *        The charge node is isolated from the instance.
             */
            if (Parameters.TrnqsMod.Value != 0)
            {
                this._qqPtr.Value += m * 1.0;
                this._qgpPtr.Value += 0.0;
                this._qdpPtr.Value += 0.0;
                this._qspPtr.Value += 0.0;
                this._qbpPtr.Value += 0.0;

                this._dpqPtr.Value += 0.0;
                this._spqPtr.Value += 0.0;
                this._gpqPtr.Value += 0.0;
            }
        }
    }
}
