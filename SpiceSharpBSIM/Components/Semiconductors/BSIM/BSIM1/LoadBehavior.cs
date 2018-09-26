using System;
using SpiceSharp.Simulations;
using SpiceSharp.Behaviors;
using SpiceSharp.Algebra;
using SpiceSharp.Components.MosfetBehaviors;
using SpiceSharp.Components.Semiconductors;
namespace SpiceSharp.Components.BSIM1Behaviors
{

    /// <summary>
    /// Load behavior for a <see cref="BSIM1" />
    /// </summary>
    public class LoadBehavior : BaseLoadBehavior, IConnectedBehavior
    {

        /// <summary>
        /// Necessary behaviors and parameters
        /// </summary>
        private ModelBaseParameters _mbp;
        private BaseParameters _bp;
        private TemperatureBehavior _temp;
        private ModelTemperatureBehavior _modelTemp;

        /// <summary>
        /// Properties
        /// </summary>
        public double Mode { get; private set; }
        public double Vdsat { get; private set; }
        public TransientBehavior TranBehavior { get; set; }
        public double Vbd { get; set; }
        public double Vbs { get; set; }
        public double Vgs { get; set; }
        public double Vds { get; set; }
        public double Gm { get; set; }
        public double Gds { get; set; }
        public double Gmbs { get; set; }
        public double Gbd { get; set; }
        public double Gbs { get; set; }
        public double Vono { get; set; }
        public double Vdsato { get; set; }
        public double Cd { get; set; }
        public double Cbs { get; set; }
        public double Cbd { get; set; }

        public double Cggb { get; private set; }
        public double Cgdb { get; private set; }
        public double Cgsb { get; private set; }
        public double Cbgb { get; private set; }
        public double Cbdb { get; private set; }
        public double Cbsb { get; private set; }
        public double Capbd { get; private set; }
        public double Iqbd { get; private set; }
        public double Capbs { get; private set; }
        public double Iqbs { get; private set; }
        public double Cdgb { get; private set; }
        public double Cddb { get; private set; }
        public double Cdsb { get; private set; }
        public double Qbs { get; private set; }
        public double Qbd { get; private set; }

        /// <summary>
        /// Nodes
        /// </summary>
        private int _drainNode, _gateNode, _sourceNode, _bulkNode;
        protected VectorElement<double> DrainNodePrimePtr { get; private set; }
        public int DrainNodePrime { get; private set; }
        protected VectorElement<double> SourceNodePrimePtr { get; private set; }
        public int SourceNodePrime { get; private set; }
        protected MatrixElement<double> DdPtr { get; private set; }
        protected MatrixElement<double> GgPtr { get; private set; }
        protected MatrixElement<double> SsPtr { get; private set; }
        protected MatrixElement<double> BbPtr { get; private set; }
        protected MatrixElement<double> DPdpPtr { get; private set; }
        protected MatrixElement<double> SPspPtr { get; private set; }
        protected MatrixElement<double> DdpPtr { get; private set; }
        protected MatrixElement<double> GbPtr { get; private set; }
        protected MatrixElement<double> GdpPtr { get; private set; }
        protected MatrixElement<double> GspPtr { get; private set; }
        protected MatrixElement<double> SspPtr { get; private set; }
        protected MatrixElement<double> BdpPtr { get; private set; }
        protected MatrixElement<double> BspPtr { get; private set; }
        protected MatrixElement<double> DPspPtr { get; private set; }
        protected MatrixElement<double> DPdPtr { get; private set; }
        protected MatrixElement<double> BgPtr { get; private set; }
        protected MatrixElement<double> DPgPtr { get; private set; }
        protected MatrixElement<double> SPgPtr { get; private set; }
        protected MatrixElement<double> SPsPtr { get; private set; }
        protected MatrixElement<double> DPbPtr { get; private set; }
        protected MatrixElement<double> SPbPtr { get; private set; }
        protected MatrixElement<double> SPdpPtr { get; private set; }
        protected VectorElement<double> GateNodePtr { get; private set; }
        protected VectorElement<double> BulkNodePtr { get; private set; }

        /// <summary>
        /// Constructor
        /// </summary>
        public LoadBehavior(Identifier name) : base(name)
        {

        }

        /// <summary>
        /// Setup the behavior
        /// </summary>
        public override void Setup(Simulation simulation, SetupDataProvider provider)
        {
            if (provider == null)
                throw new ArgumentNullException(nameof(provider));

            // Get parameters
            _bp = provider.GetParameterSet<BaseParameters>();
            _mbp = provider.GetParameterSet<ModelBaseParameters>("model");

            // Get behaviors
            _temp = provider.GetBehavior<TemperatureBehavior>();
            _modelTemp = provider.GetBehavior<ModelTemperatureBehavior>("model");
        }

        /// <summary>
        /// Connect
        /// </summary>
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
        public override void GetEquationPointers(VariableSet variables, Solver<double> solver)
        {
            if (_mbp.SheetResistance > 0 && _bp.DrainSquares > 0.0)
                DrainNodePrime = variables.Create(new SubIdentifier(Name, "drain")).Index;
            else
                DrainNodePrime = _drainNode;
            DrainNodePrimePtr = solver.GetRhsElement(DrainNodePrime);

            if (_mbp.SheetResistance > 0 && _bp.SourceSquares > 0.0)
                SourceNodePrime = variables.Create(new SubIdentifier(Name, "source")).Index;
            else
                SourceNodePrime = _sourceNode;
            SourceNodePrimePtr = solver.GetRhsElement(SourceNodePrime);

            DdPtr = solver.GetMatrixElement(_drainNode, _drainNode);
            GgPtr = solver.GetMatrixElement(_gateNode, _gateNode);
            SsPtr = solver.GetMatrixElement(_sourceNode, _sourceNode);
            BbPtr = solver.GetMatrixElement(_bulkNode, _bulkNode);
            DPdpPtr = solver.GetMatrixElement(DrainNodePrime, DrainNodePrime);
            SPspPtr = solver.GetMatrixElement(SourceNodePrime, SourceNodePrime);
            DdpPtr = solver.GetMatrixElement(_drainNode, DrainNodePrime);
            GbPtr = solver.GetMatrixElement(_gateNode, _bulkNode);
            GdpPtr = solver.GetMatrixElement(_gateNode, DrainNodePrime);
            GspPtr = solver.GetMatrixElement(_gateNode, SourceNodePrime);
            SspPtr = solver.GetMatrixElement(_sourceNode, SourceNodePrime);
            BdpPtr = solver.GetMatrixElement(_bulkNode, DrainNodePrime);
            BspPtr = solver.GetMatrixElement(_bulkNode, SourceNodePrime);
            DPspPtr = solver.GetMatrixElement(DrainNodePrime, SourceNodePrime);
            DPdPtr = solver.GetMatrixElement(DrainNodePrime, _drainNode);
            BgPtr = solver.GetMatrixElement(_bulkNode, _gateNode);
            DPgPtr = solver.GetMatrixElement(DrainNodePrime, _gateNode);
            SPgPtr = solver.GetMatrixElement(SourceNodePrime, _gateNode);
            SPsPtr = solver.GetMatrixElement(SourceNodePrime, _sourceNode);
            DPbPtr = solver.GetMatrixElement(DrainNodePrime, _bulkNode);
            SPbPtr = solver.GetMatrixElement(SourceNodePrime, _bulkNode);
            SPdpPtr = solver.GetMatrixElement(SourceNodePrime, DrainNodePrime);
            GateNodePtr = solver.GetRhsElement(_gateNode);
            BulkNodePtr = solver.GetRhsElement(_bulkNode);
        }

        /// <summary>
        /// Load the behavior
        /// </summary>
        /// <param name="simulation">Simulation</param>
        public override void Load(BaseSimulation simulation)
        {
            var state = simulation.RealState;
            double DrainSatCurrent = 0.0;
            double EffectiveLength = 0.0;
            double GateBulkOverlapCap = 0.0;
            double GateDrainOverlapCap = 0.0;
            double GateSourceOverlapCap = 0.0;
            double SourceSatCurrent = 0.0;
            double DrainArea = 0.0;
            double SourceArea = 0.0;
            double DrainPerimeter = 0.0;
            double SourcePerimeter = 0.0;
            double arg = 0.0;
            double capbd = 0.0;
            double capbs = 0.0;
            double cbd = 0.0;
            double cbhat = 0.0;
            double cbs = 0.0;
            double cd = 0.0;
            double cdrain = 0.0;
            double cdhat = 0.0;
            double cdreq = 0.0;
            double ceq = 0.0;
            double ceqbd = 0.0;
            double ceqbs = 0.0;
            double ceqqb = 0.0;
            double ceqqd = 0.0;
            double ceqqg = 0.0;
            double czbd = 0.0;
            double czbdsw = 0.0;
            double czbs = 0.0;
            double czbssw = 0.0;
            double delvbd = 0.0;
            double delvbs = 0.0;
            double delvds = 0.0;
            double delvgd = 0.0;
            double delvgs = 0.0;
            double evbd = 0.0;
            double evbs = 0.0;
            double gbd = 0.0;
            double gbs = 0.0;
            double gcbdb = 0.0;
            double gcbgb = 0.0;
            double gcbsb = 0.0;
            double gcddb = 0.0;
            double gcdgb = 0.0;
            double gcdsb = 0.0;
            double gcgdb = 0.0;
            double gcggb = 0.0;
            double gcgsb = 0.0;
            double gcsdb = 0.0;
            double gcsgb = 0.0;
            double gcssb = 0.0;
            double gds = 0.0;
            double geq = 0.0;
            double gm = 0.0;
            double gmbs = 0.0;
            double sarg = 0.0;
            double sargsw = 0.0;
            double tol = 0.0;
            double vbd = 0.0;
            double vbs = 0.0;
            double vcrit = 0.0;
            double vds = 0.0;
            double vdsat = 0.0;
            double vgb = 0.0;
            double vgd = 0.0;
            double vgdo = 0.0;
            double vgs = 0.0;
            double von = 0.0;
            double xfact = 0.0;
            double xnrm = 0.0;
            double xrev = 0.0;
            bool Check = false;
            double cgdb = 0.0;
            double cgsb = 0.0;
            double cbdb = 0.0;
            double cdgb = 0.0;
            double cddb = 0.0;
            double cdsb = 0.0;
            double cggb = 0.0;
            double cbgb = 0.0;
            double cbsb = 0.0;
            double csgb = 0.0;
            double cssb = 0.0;
            double csdb = 0.0;
            double PhiB = 0.0;
            double PhiBSW = 0.0;
            double MJ = 0.0;
            double MJSW = 0.0;
            double argsw = 0.0;
            double qgate = 0.0;
            double qbulk = 0.0;
            double qdrn = 0.0;
            double qsrc = 0.0;
            double cqgate = 0.0;
            double cqbulk = 0.0;
            double cqdrn = 0.0;
            double vt0 = 0.0;
            double[] args = new double[8];
            int ByPass = 0;
            int error = 0;

            EffectiveLength = _bp.Length - _mbp.DeltaL * 1.0e-6; /* m */
            DrainArea = _bp.DrainArea;
            SourceArea = _bp.SourceArea;
            DrainPerimeter = _bp.DrainPerimeter;
            SourcePerimeter = _bp.SourcePerimeter;
            if ((DrainSatCurrent = DrainArea * _mbp.JctSatCurDensity)
                < 1e-15)
            {
                DrainSatCurrent = 1.0e-15;
            }

            if ((SourceSatCurrent = SourceArea * _mbp.JctSatCurDensity)
                < 1.0e-15)
            {
                SourceSatCurrent = 1.0e-15;
            }

            GateSourceOverlapCap = _mbp.GateSourceOverlapCap * _bp.Width;
            GateDrainOverlapCap = _mbp.GateDrainOverlapCap * _bp.Width;
            GateBulkOverlapCap = _mbp.GateBulkOverlapCap * EffectiveLength;
            von = _mbp.Type * _temp.Von;
            vdsat = _mbp.Type * this.Vdsat;
            vt0 = _mbp.Type * _temp.Vt0;

            Check = true;
            if (state.Domain == RealState.DomainType.Frequency)
            {
                vbs = this.Vbs;
                vgs = this.Vgs;
                vds = this.Vds;
            }
            else if (state.Init == RealState.InitializationStates.InitJunction && !_bp.Off)
            {
                vds = _mbp.Type * _bp.IcVDS;
                vgs = _mbp.Type * _bp.IcVGS;
                vbs = _mbp.Type * _bp.IcVBS;
                if ((vds == 0) && (vgs == 0) && (vbs == 0) &&
                    (state.Domain == RealState.DomainType.None || state.UseDc || TranBehavior != null || !state.UseIc))
                {
                    vbs = -1;
                    vgs = vt0;
                    vds = 0;
                }
            }
            else if ((state.Init == RealState.InitializationStates.InitJunction ||
                      state.Init == RealState.InitializationStates.InitFix) && (_bp.Off))
            {
                vbs = vgs = vds = 0;
            }
            else
            {
                vbs = _mbp.Type * (state.Solution[_bulkNode] - state.Solution[SourceNodePrime]);
                vgs = _mbp.Type * (state.Solution[_gateNode] - state.Solution[SourceNodePrime]);
                vds = _mbp.Type * (state.Solution[DrainNodePrime] - state.Solution[SourceNodePrime]);
                vbd = vbs - vds;
                vgd = vgs - vds;
                vgdo = this.Vgs - this.Vds;

                von = _mbp.Type * _temp.Von;
                if (this.Vds >= 0)
                {
                    vgs = Transistor.LimitFet(vgs, this.Vgs, von);
                    vds = vgs - vgd;
                    vds = Transistor.LimitVoltageDs(vds, this.Vds);
                    vgd = vgs - vds;
                }
                else
                {
                    vgd = Transistor.LimitFet(vgd, vgdo, von);
                    vds = vgs - vgd;
                    vds = -Transistor.LimitVoltageDs(-vds, -this.Vds);
                    vgs = vgd + vds;
                }

                Check = false;
                if (vds >= 0)
                {
                    vcrit = Circuit.Vt0 * Math.Log(Circuit.Vt0 / (Circuit.Root2 * SourceSatCurrent));
                    vbs = Semiconductor.LimitJunction(vbs, this.Vbs, Circuit.Vt0, vcrit, ref Check); /* B1 test */
                    vbd = vbs - vds;
                }
                else
                {
                    vcrit = Circuit.Vt0 * Math.Log(Circuit.Vt0 / (Circuit.Root2 * DrainSatCurrent));
                    vbd = Semiconductor.LimitJunction(vbd, this.Vbd, Circuit.Vt0, vcrit, ref Check); /* B1 test*/
                    vbs = vbd + vds;
                }
            }

            /* determine DC current and derivatives */
            vbd = vbs - vds;
            vgd = vgs - vds;
            vgb = vgs - vbs;


            if (vbs <= 0.0)
            {
                gbs = SourceSatCurrent / Circuit.Vt0 + state.Gmin;
                cbs = gbs * vbs;
            }
            else
            {
                evbs = Math.Exp(vbs / Circuit.Vt0);
                gbs = SourceSatCurrent * evbs / Circuit.Vt0 + state.Gmin;
                cbs = SourceSatCurrent * (evbs - 1) + state.Gmin * vbs;
            }

            if (vbd <= 0.0)
            {
                gbd = DrainSatCurrent / Circuit.Vt0 + state.Gmin;
                cbd = gbd * vbd;
            }
            else
            {
                evbd = Math.Exp(vbd / Circuit.Vt0);
                gbd = DrainSatCurrent * evbd / Circuit.Vt0 + state.Gmin;
                cbd = DrainSatCurrent * (evbd - 1) + state.Gmin * vbd;
            }

            /* line 400 */
            if (vds >= 0)
            {
                /* normal mode */
                this.Mode = 1;
            }
            else
            {
                /* inverse mode */
                this.Mode = -1;
            }

            /* call B1evaluate to calculate drain current and its 
             * derivatives and charge and capacitances related to gate
             * drain, and bulk
             */
            if (vds >= 0)
            {
                Evaluate(TranBehavior != null || state.Domain == RealState.DomainType.Frequency,
                    vds, vbs, vgs, out gm, out gds, out gmbs, out qgate,
                    out qbulk, out qdrn, out cggb, out cgdb, out cgsb, out cbgb, out cbdb, out cbsb, out cdgb,
                    out cddb, out cdsb, out cdrain, out von, out vdsat);
            }
            else
            {
                Evaluate(TranBehavior != null || state.Domain == RealState.DomainType.Frequency,
                    -vds, vbd, vgd, out gm, out gds, out gmbs, out qgate,
                    out qbulk, out qsrc, out cggb, out cgsb, out cgdb, out cbgb, out cbsb, out cbdb, out csgb,
                    out cssb, out csdb, out cdrain, out von, out vdsat);
            }

            _temp.Von = _mbp.Type * von;
            this.Vdsat = _mbp.Type * vdsat;

            /*
             *  COMPUTE EQUIVALENT DRAIN CURRENT SOURCE
             */
            cd = this.Mode * cdrain - cbd;
            if (TranBehavior != null || state.Domain == RealState.DomainType.Frequency)
            {
                /*
                 *  charge storage elements
                 *
                 *   bulk-drain and bulk-source depletion capacitances
                 *  czbd : zero bias drain junction capacitance
                 *  czbs : zero bias source junction capacitance
                 * czbdsw:zero bias drain junction sidewall capacitance
                 * czbssw:zero bias source junction sidewall capacitance
                 */

                czbd = _mbp.UnitAreaJctCap * DrainArea;
                czbs = _mbp.UnitAreaJctCap * SourceArea;
                czbdsw = _mbp.UnitLengthSidewallJctCap * DrainPerimeter;
                czbssw = _mbp.UnitLengthSidewallJctCap * SourcePerimeter;
                PhiB = _mbp.BulkJctPotential;
                PhiBSW = _mbp.SidewallJctPotential;
                MJ = _mbp.BulkJctBotGradingCoeff;
                MJSW = _mbp.BulkJctSideGradingCoeff;

                /* Source Bulk Junction */
                if (vbs < 0)
                {
                    arg = 1 - vbs / PhiB;
                    argsw = 1 - vbs / PhiBSW;
                    sarg = Math.Exp(-MJ * Math.Log(arg));
                    sargsw = Math.Exp(-MJSW * Math.Log(argsw));
                    this.Qbs = PhiB * czbs * (1 - arg * sarg) / (1 - MJ) +
                               PhiBSW * czbssw * (1 - argsw * sargsw) / (1 - MJSW);
                    capbs = czbs * sarg + czbssw * sargsw;
                }
                else
                {
                    this.Qbs = vbs * (czbs + czbssw) +
                               vbs * vbs * (czbs * MJ * 0.5 / PhiB + czbssw * MJSW * 0.5 / PhiBSW);
                    capbs = czbs + czbssw + vbs * (czbs * MJ / PhiB + czbssw * MJSW / PhiBSW);
                }

                /* Drain Bulk Junction */
                if (vbd < 0)
                {
                    arg = 1 - vbd / PhiB;
                    argsw = 1 - vbd / PhiBSW;
                    sarg = Math.Exp(-MJ * Math.Log(arg));
                    sargsw = Math.Exp(-MJSW * Math.Log(argsw));

                    this.Qbd = PhiB * czbd * (1 - arg * sarg) / (1 - MJ) +
                               PhiBSW * czbdsw * (1 - argsw * sargsw) / (1 - MJSW);
                    capbd = czbd * sarg + czbdsw * sargsw;
                }
                else
                {
                    this.Qbd = vbd * (czbd + czbdsw) +
                               vbd * vbd * (czbd * MJ * 0.5 / PhiB + czbdsw * MJSW * 0.5 / PhiBSW);
                    capbd = czbd + czbdsw + vbd * (czbd * MJ / PhiB + czbdsw * MJSW / PhiBSW);
                }
            }

            /*
             *  check convergence
             */
            if ((!_bp.Off) || state.Init != RealState.InitializationStates.InitFix)
            {
                if (Check)
                {
                    state.IsConvergent = false;
                }
            }


            this.Vbs = vbs;
            this.Vbd = vbd;
            this.Vgs = vgs;
            this.Vds = vds;
            this.Cd = cd;
            this.Cbs = cbs;
            this.Cbd = cbd;
            this.Gm = gm;
            this.Gds = gds;
            this.Gmbs = gmbs;
            this.Gbd = gbd;
            this.Gbs = gbs;

            this.Cggb = cggb;
            this.Cgdb = cgdb;
            this.Cgsb = cgsb;

            this.Cbgb = cbgb;
            this.Cbdb = cbdb;
            this.Cbsb = cbsb;

            this.Cdgb = cdgb;
            this.Cddb = cddb;
            this.Cdsb = cdsb;

            this.Capbs = capbs;
            this.Capbd = capbd;

            /* bulk and channel charge plus overlaps */
            // if (!(ckt->CKTmode & (MODETRAN | MODEAC)) && (!(ckt->CKTmode & MODETRANOP) || !(ckt->CKTmode & MODEUIC)) && !(ckt->CKTmode & MODEINITSMSIG))
            if (TranBehavior == null && state.Domain != RealState.DomainType.Frequency)
                goto line850;

            line755:
            if (this.Mode > 0)
            {

                args[0] = GateDrainOverlapCap;
                args[1] = GateSourceOverlapCap;
                args[2] = GateBulkOverlapCap;
                args[3] = capbd;
                args[4] = capbs;
                args[5] = cggb;
                args[6] = cgdb;
                args[7] = cgsb;

                MosCap(vgd, vgs, vgb,
                    args,
                    /*
                    GateDrainOverlapCap,
                    GateSourceOverlapCap,GateBulkOverlapCap,
                    capbd,capbs,
                                cggb,cgdb,cgsb,
                    */
                    cbgb, cbdb, cbsb,
                    cdgb, cddb, cdsb,
                    out gcggb, out gcgdb, out gcgsb,
                    out gcbgb, out gcbdb, out gcbsb,
                    out gcdgb, out gcddb, out gcdsb, out gcsgb, out gcsdb, out gcssb,
                    ref qgate, ref qbulk,
                    ref qdrn, out qsrc);
            }
            else
            {

                args[0] = GateSourceOverlapCap;
                args[1] = GateDrainOverlapCap;
                args[2] = GateBulkOverlapCap;
                args[3] = capbs;
                args[4] = capbd;
                args[5] = cggb;
                args[6] = cgsb;
                args[7] = cgdb;

                MosCap(vgs, vgd, vgb,
                    args,
                    /*
                    GateSourceOverlapCap,
                    GateDrainOverlapCap,GateBulkOverlapCap,
                    capbs,capbd,
                    cggb,cgsb,cgdb,
                    */
                    cbgb, cbsb, cbdb,
                    csgb, cssb, csdb,
                    out gcggb, out gcgsb, out gcgdb,
                    out gcbgb, out gcbsb, out gcbdb,
                    out gcsgb, out gcssb, out gcsdb, out gcdgb, out gcdsb, out gcddb,
                    ref qgate, ref qbulk,
                    ref qsrc, out qdrn);
            }



            // store small signal parameters
            // if ((!(ckt->CKTmode & (MODEAC | MODETRAN))) && (ckt->CKTmode & MODETRANOP) && (ckt->CKTmode & MODEUIC))
                // goto line850;

            if (state.Domain == RealState.DomainType.Frequency)
            {
                this.Cggb = cggb;
                this.Cgdb = cgdb;
                this.Cgsb = cgsb;
                this.Cbgb = cbgb;
                this.Cbdb = cbdb;
                this.Cbsb = cbsb;
                this.Cdgb = cdgb;
                this.Cddb = cddb;
                this.Cdsb = cdsb;
                this.Capbd = capbd;
                this.Capbs = capbs;
                goto line1000;
            }

            if (TranBehavior != null)
            {
                TranBehavior.Qg.Current = qgate;
                TranBehavior.Qd.Current = qdrn - this.Qbd;
                TranBehavior.Qb.Current = qbulk + this.Qbd + this.Qbs;

                if (!state.UseDc)
                {
                    TranBehavior.Qb.Integrate();
                    TranBehavior.Qg.Integrate();
                    TranBehavior.Qd.Integrate();
                }
            }

            goto line860;

            line850:
            /* initialize to zero charge conductance and current */
            ceqqg = ceqqb = ceqqd = 0.0;
            gcdgb = gcddb = gcdsb = 0.0;
            gcsgb = gcsdb = gcssb = 0.0;
            gcggb = gcgdb = gcgsb = 0.0;
            gcbgb = gcbdb = gcbsb = 0.0;
            goto line900;

            line860:
            /* evaluate equivalent charge current */
            cqgate = TranBehavior.Qg.Derivative;
            cqbulk = TranBehavior.Qb.Derivative;
            cqdrn = TranBehavior.Qd.Derivative;
            ceqqg = cqgate - gcggb * vgb + gcgdb * vbd + gcgsb * vbs;
            ceqqb = cqbulk - gcbgb * vgb + gcbdb * vbd + gcbsb * vbs;
            ceqqd = cqdrn - gcdgb * vgb + gcddb * vbd + gcdsb * vbs;

            /*
             *  load current vector
             */
            line900:
            ceqbs = _mbp.Type * (cbs - (gbs - state.Gmin) * vbs);
            ceqbd = _mbp.Type * (cbd - (gbd - state.Gmin) * vbd);

            ceqqg = _mbp.Type * ceqqg;
            ceqqb = _mbp.Type * ceqqb;
            ceqqd = _mbp.Type * ceqqd;
            if (this.Mode >= 0)
            {
                xnrm = 1;
                xrev = 0;
                cdreq = _mbp.Type * (cdrain - gds * vds - gm * vgs - gmbs * vbs);
            }
            else
            {
                xnrm = 0;
                xrev = 1;
                cdreq = -(_mbp.Type) * (cdrain + gds * vds - gm * vgd - gmbs * vbd);
            }

            var m = _bp.Multiplier;
            GateNodePtr.Value -= m * ceqqg;
            BulkNodePtr.Value -= m * (ceqbs + ceqbd + ceqqb);
            DrainNodePrimePtr.Value += m * (ceqbd - cdreq - ceqqd);
            SourceNodePrimePtr.Value += m * (cdreq + ceqbs + ceqqg + ceqqb + ceqqd);

            /*
             *  load y matrix
             */
            DdPtr.Value += m * (_temp.DrainConductance);
            GgPtr.Value += m * (gcggb);
            SsPtr.Value += m * (_temp.SourceConductance);
            BbPtr.Value += m * (gbd + gbs - gcbgb - gcbdb - gcbsb);
            DPdpPtr.Value += m * (_temp.DrainConductance + gds + gbd + xrev * (gm + gmbs) + gcddb);
            SPspPtr.Value += m * (_temp.SourceConductance + gds + gbs + xnrm * (gm + gmbs) + gcssb);
            DdpPtr.Value += m * (-_temp.DrainConductance);
            GbPtr.Value += m * (-gcggb - gcgdb - gcgsb);
            GdpPtr.Value += m * (gcgdb);
            GspPtr.Value += m * (gcgsb);
            SspPtr.Value += m * (-_temp.SourceConductance);
            BgPtr.Value += m * (gcbgb);
            BdpPtr.Value += m * (-gbd + gcbdb);
            BspPtr.Value += m * (-gbs + gcbsb);
            DPdPtr.Value += m * (-_temp.DrainConductance);
            DPgPtr.Value += m * ((xnrm - xrev) * gm + gcdgb);
            DPbPtr.Value += m * (-gbd + (xnrm - xrev) * gmbs - gcdgb - gcddb - gcdsb);
            DPspPtr.Value += m * (-gds - xnrm * (gm + gmbs) + gcdsb);
            SPgPtr.Value += m * (-(xnrm - xrev) * gm + gcsgb);
            SPsPtr.Value += m * (-_temp.SourceConductance);
            SPbPtr.Value += m * (-gbs - (xnrm - xrev) * gmbs - gcsgb - gcsdb - gcssb);
            SPdpPtr.Value += m * (-gds - xrev * (gm + gmbs) + gcsdb);

            line1000: ;
        }

        /// <summary>
        /// Helper method Evaluate
        /// </summary>
        public void Evaluate(bool chargeComputationNeeded, double vds, double vbs, double vgs, out double gmPointer, out double gdsPointer, out double gmbsPointer, out double qgPointer, out double qbPointer, out double qdPointer, out double cggbPointer, out double cgdbPointer, out double cgsbPointer, out double cbgbPointer, out double cbdbPointer, out double cbsbPointer, out double cdgbPointer, out double cddbPointer, out double cdsbPointer, out double cdrainPointer, out double vonPointer, out double vdsatPointer)
        {
            double gm = 0, gds = 0, gmbs = 0, qg = 0, qb = 0, qd = 0, cggb = 0, cgdb = 0, cgsb = 0, cbgb = 0, cbdb = 0, cbsb = 0, cdgb = 0, cddb = 0, cdsb = 0, vfb, phi, k1, k2, vdd, ugs, uds, dUgsdVbs, leff, dUdsdVbs, dUdsdVds, eta, dEtadVds, dEtadVbs, vpb, sqrtVpb, von, vth, dVthdVbs, dVthdVds, vgs_Vth, drainCurrent = 0, g, a, arg, dGdVbs, dAdVbs, beta, beta_Vds_0, betaVdd, dBetaVdd_dVds, beta0, dBeta0dVds, dBeta0dVbs, vddSquare, c1, c2, dBetaVdd_dVbs, dBeta_Vds_0_dVbs, dC1dVbs, dC2dVbs, dBetadVgs, dBetadVds, dBetadVbs, vdsSat = 0, argl1, argl2, vc, term1, k, args1, dVcdVgs, dVcdVds, dVcdVbs, dKdVc, dKdVgs, dKdVds, dKdVbs, args2, args3, warg1, n, n0, nB, nD, warg2, wds, wgs, ilimit, iexp, temp1, vth0, arg1, arg2, arg3, arg5, ent, vcom, vgb, vgb_Vfb, vdsPinchoff, entSquare, vgs_VthSquare, argl3, argl4, argl5, argl6, argl7, argl8, argl9, dEntdVds, dEntdVbs, cgbb, cdbb, cbbb, wLCox, vtsquare, temp3, co4v15;
            vfb = _temp.Vfb;
            phi = _temp.Phi;
            k1 = _temp.K1;
            k2 = _temp.K2;
            vdd = _mbp.Vdd;
            if ((ugs = _temp.Ugs + _temp.UgsB * vbs) <= 0)
            {
                ugs = 0;
                dUgsdVbs = 0.0;
            }
            else
            {
                dUgsdVbs = _temp.UgsB;
            }
            if ((uds = _temp.Uds + _temp.UdsB * vbs + _temp.UdsD * (vds - vdd)) <= 0)
            {
                uds = 0.0;
                dUdsdVbs = dUdsdVds = 0.0;
            }
            else
            {
                leff = _bp.Length * 1.0e6 - _mbp.DeltaL;
                uds = uds / leff;
                dUdsdVbs = _temp.UdsB / leff;
                dUdsdVds = _temp.UdsD / leff;
            }
            eta = _temp.Eta + _temp.EtaB * vbs + _temp.EtaD * (vds - vdd);
            if (eta <= 0)
            {
                eta = 0;
                dEtadVds = dEtadVbs = 0.0;
            }
            else if (eta > 1)
            {
                eta = 1;
                dEtadVds = dEtadVbs = 0;
            }
            else
            {
                dEtadVds = _temp.EtaD;
                dEtadVbs = _temp.EtaB;
            }
            if (vbs < 0)
            {
                vpb = phi - vbs;
            }
            else
            {
                vpb = phi;
            }
            sqrtVpb = Math.Sqrt(vpb);
            von = vfb + phi + k1 * sqrtVpb - k2 * vpb - eta * vds;
            vth = von;
            dVthdVds = -eta - dEtadVds * vds;
            dVthdVbs = k2 - 0.5 * k1 / sqrtVpb - dEtadVbs * vds;
            vgs_Vth = vgs - vth;
            g = 1.0 - 1.0 / (1.744 + 0.8364 * vpb);
            a = 1.0 + 0.5 * g * k1 / sqrtVpb;
            a = Math.Max(a, 1.0);
            arg = Math.Max(1 + ugs * vgs_Vth, 1.0);
            dGdVbs = -0.8364 * (1 - g) * (1 - g);
            dAdVbs = 0.25 * k1 / sqrtVpb * (2 * dGdVbs + g / vpb);
            if (vgs_Vth < 0)
            {
                goto SubthresholdComputation;
            }
            beta_Vds_0 = _temp.BetaZero + _temp.BetaZeroB * vbs;
            betaVdd = _temp.BetaVdd + _temp.BetaVddB * vbs;
            dBetaVdd_dVds = Math.Max(_temp.BetaVddD, 0.0);
            if (vds > vdd)
            {
                beta0 = betaVdd + dBetaVdd_dVds * (vds - vdd);
                dBeta0dVds = dBetaVdd_dVds;
                dBeta0dVbs = _temp.BetaVddB;
            }
            else
            {
                vddSquare = vdd * vdd;
                c1 = (-betaVdd + beta_Vds_0 + dBetaVdd_dVds * vdd) / vddSquare;
                c2 = 2 * (betaVdd - beta_Vds_0) / vdd - dBetaVdd_dVds;
                dBeta_Vds_0_dVbs = _temp.BetaZeroB;
                dBetaVdd_dVbs = _temp.BetaVddB;
                dC1dVbs = (dBeta_Vds_0_dVbs - dBetaVdd_dVbs) / vddSquare;
                dC2dVbs = dC1dVbs * -2 * vdd;
                beta0 = (c1 * vds + c2) * vds + beta_Vds_0;
                dBeta0dVds = 2 * c1 * vds + c2;
                dBeta0dVbs = dC1dVbs * vds * vds + dC2dVbs * vds + dBeta_Vds_0_dVbs;
            }
            beta = beta0 / arg;
            dBetadVgs = -beta * ugs / arg;
            dBetadVds = dBeta0dVds / arg - dBetadVgs * dVthdVds;
            dBetadVbs = dBeta0dVbs / arg + beta * ugs * dVthdVbs / arg - beta * vgs_Vth * dUgsdVbs / arg;
            if ((vc = uds * vgs_Vth / a) < 0.0)
                vc = 0.0;
            term1 = Math.Sqrt(1 + 2 * vc);
            k = 0.5 * (1 + vc + term1);
            vdsSat = Math.Max(vgs_Vth / (a * Math.Sqrt(k)), 0.0);
            if (vds < vdsSat)
            {
                argl1 = Math.Max(1 + uds * vds, 1.0);
                argl2 = vgs_Vth - 0.5 * a * vds;
                drainCurrent = beta * argl2 * vds / argl1;
                gm = (dBetadVgs * argl2 * vds + beta * vds) / argl1;
                gds = (dBetadVds * argl2 * vds + beta * (vgs_Vth - vds * dVthdVds - a * vds) - drainCurrent * (vds * dUdsdVds + uds)) / argl1;
                gmbs = (dBetadVbs * argl2 * vds + beta * vds * (-dVthdVbs - 0.5 * vds * dAdVbs) - drainCurrent * vds * dUdsdVbs) / argl1;
            }
            else
            {
                args1 = 1.0 + 1.0 / term1;
                dVcdVgs = uds / a;
                dVcdVds = vgs_Vth * dUdsdVds / a - dVcdVgs * dVthdVds;
                dVcdVbs = (vgs_Vth * dUdsdVbs - uds * (dVthdVbs + vgs_Vth * dAdVbs / a)) / a;
                dKdVc = 0.5 * args1;
                dKdVgs = dKdVc * dVcdVgs;
                dKdVds = dKdVc * dVcdVds;
                dKdVbs = dKdVc * dVcdVbs;
                args2 = vgs_Vth / a / k;
                args3 = args2 * vgs_Vth;
                drainCurrent = 0.5 * beta * args3;
                gm = 0.5 * args3 * dBetadVgs + beta * args2 - drainCurrent * dKdVgs / k;
                gds = 0.5 * args3 * dBetadVds - beta * args2 * dVthdVds - drainCurrent * dKdVds / k;
                gmbs = 0.5 * dBetadVbs * args3 - beta * args2 * dVthdVbs - drainCurrent * (dAdVbs / a + dKdVbs / k);
            }
            SubthresholdComputation:
            n0 = _temp.SubthSlope;
            if (n0 >= 200)
            {
                goto ChargeComputation;
            }
            nB = _temp.SubthSlopeB;
            nD = _temp.SubthSlopeD;
            n = n0 + nB * vbs + nD * vds;
            if (n < 0.5)
                n = 0.5;
            warg1 = Math.Exp(-vds / Circuit.Vt0);
            wds = 1 - warg1;
            wgs = Math.Exp(vgs_Vth / (n * Circuit.Vt0));
            vtsquare = Circuit.Vt0 * Circuit.Vt0;
            warg2 = 6.04965 * vtsquare * _temp.BetaZero;
            ilimit = 4.5 * vtsquare * _temp.BetaZero;
            iexp = warg2 * wgs * wds;
            drainCurrent = drainCurrent + ilimit * iexp / (ilimit + iexp);
            temp1 = ilimit / (ilimit + iexp);
            temp1 = temp1 * temp1;
            temp3 = ilimit / (ilimit + wgs * warg2);
            temp3 = temp3 * temp3 * warg2 * wgs;
            gm = gm + temp1 * iexp / (n * Circuit.Vt0);
            gds = gds + temp3 * (-wds / n / Circuit.Vt0 * (dVthdVds + vgs_Vth * nD / n) + warg1 / Circuit.Vt0);
            gmbs = gmbs - temp1 * iexp * (dVthdVbs + vgs_Vth * nB / n) / (n * Circuit.Vt0);
            ChargeComputation:
            if (drainCurrent < 0.0)
                drainCurrent = 0.0;
            if (gm < 0.0)
                gm = 0.0;
            if (gds < 0.0)
                gds = 0.0;
            if (gmbs < 0.0)
                gmbs = 0.0;
            wLCox = _modelTemp.Cox * (_bp.Length - _mbp.DeltaL * 1.0e-6) * (_bp.Width - _mbp.DeltaW * 1.0e-6) * 1.0e4;
            if (!chargeComputationNeeded)
            {
                qg = 0;
                qd = 0;
                qb = 0;
                cggb = 0;
                cgsb = 0;
                cgdb = 0;
                cdgb = 0;
                cdsb = 0;
                cddb = 0;
                cbgb = 0;
                cbsb = 0;
                cbdb = 0;
                goto finished;
            }
            g = 1.0 - 1.0 / (1.744 + 0.8364 * vpb);
            a = 1.0 + 0.5 * g * k1 / sqrtVpb;
            a = Math.Max(a, 1.0);
            phi = Math.Max(0.1, phi);
            if (_mbp.ChannelChargePartitionFlag >= 1)
            {
                vth0 = vfb + phi + k1 * sqrtVpb;
                vgs_Vth = vgs - vth0;
                arg1 = a * vds;
                arg2 = vgs_Vth - 0.5 * arg1;
                arg3 = vds - arg1;
                dVthdVbs = -0.5 * k1 / sqrtVpb;
                dAdVbs = 0.5 * k1 * (0.5 * g / vpb - 0.8364 * (1 - g) * (1 - g)) / sqrtVpb;
                ent = Math.Max(arg2, 1.0e-8);
                dEntdVds = -0.5 * a;
                dEntdVbs = -dVthdVbs - 0.5 * vds * dAdVbs;
                vdsPinchoff = Math.Max(vgs_Vth / a, 0.0);
                vgb = vgs - vbs;
                vgb_Vfb = vgb - vfb;
                if (vgb_Vfb < 0)
                {
                    qg = wLCox * vgb_Vfb;
                    qb = -qg;
                    qd = 0.0;
                    cggb = wLCox;
                    cgdb = 0.0;
                    cgsb = 0.0;
                    cbgb = -wLCox;
                    cbdb = 0.0;
                    cbsb = 0.0;
                    cdgb = 0.0;
                    cddb = 0.0;
                    cdsb = 0.0;
                    goto finished;
                }
                else if (vgs < vth0)
                {
                    qg = 0.5 * wLCox * k1 * k1 * (-1 + Math.Sqrt(1 + 4 * vgb_Vfb / (k1 * k1)));
                    qb = -qg;
                    qd = 0.0;
                    cggb = wLCox / Math.Sqrt(1 + 4 * vgb_Vfb / (k1 * k1));
                    cgdb = cgsb = 0.0;
                    cbgb = -cggb;
                    cbdb = cbsb = cdgb = cddb = cdsb = 0.0;
                    goto finished;
                }
                else if (vds < vdsPinchoff)
                {
                    entSquare = ent * ent;
                    argl1 = 1.2e1 * entSquare;
                    argl2 = 1.0 - a;
                    argl3 = arg1 * vds;
                    if (ent > 1.0e-8)
                    {
                        argl5 = arg1 / ent;
                    }
                    else
                    {
                        argl5 = 2.0;
                    }
                    argl7 = argl5 / 1.2e1;
                    argl8 = 6.0 * ent;
                    argl9 = 0.125 * argl5 * argl5;
                    qg = wLCox * (vgs - vfb - phi - 0.5 * vds + vds * argl7);
                    qb = wLCox * (-vth0 + vfb + phi + 0.5 * arg3 - arg3 * argl7);
                    qd = -wLCox * (0.5 * vgs_Vth - 0.75 * arg1 + 0.125 * arg1 * argl5);
                    cggb = wLCox * (1.0 - argl3 / argl1);
                    cgdb = wLCox * (-0.5 + arg1 / argl8 - argl3 * dEntdVds / argl1);
                    cgbb = wLCox * (vds * vds * dAdVbs * ent - argl3 * dEntdVbs) / argl1;
                    cgsb = -(cggb + cgdb + cgbb);
                    cbgb = wLCox * argl3 * argl2 / argl1;
                    cbdb = wLCox * argl2 * (0.5 - arg1 / argl8 + argl3 * dEntdVds / argl1);
                    cbbb = -wLCox * (dVthdVbs + 0.5 * vds * dAdVbs + vds * vds * ((1.0 - 2.0 * a) * dAdVbs * ent - argl2 * a * dEntdVbs) / argl1);
                    cbsb = -(cbgb + cbdb + cbbb);
                    cdgb = -wLCox * (0.5 - argl9);
                    cddb = wLCox * (0.75 * a - 0.25 * a * arg1 / ent + argl9 * dEntdVds);
                    cdbb = wLCox * (0.5 * dVthdVbs + vds * dAdVbs * (0.75 - 0.25 * argl5) + argl9 * dEntdVbs);
                    cdsb = -(cdgb + cddb + cdbb);
                    goto finished;
                }
                else if (vds >= vdsPinchoff)
                {
                    args1 = 1.0 / (3.0 * a);
                    qg = wLCox * (vgs - vfb - phi - vgs_Vth * args1);
                    qb = wLCox * (vfb + phi - vth0 + (1.0 - a) * vgs_Vth * args1);
                    qd = 0.0;
                    cggb = wLCox * (1.0 - args1);
                    cgdb = 0.0;
                    cgbb = wLCox * args1 * (dVthdVbs + vgs_Vth * dAdVbs / a);
                    cgsb = -(cggb + cgdb + cgbb);
                    cbgb = wLCox * (args1 - 1.0 / 3.0);
                    cbdb = 0.0;
                    cbbb = -wLCox * ((2.0 / 3.0 + args1) * dVthdVbs + vgs_Vth * args1 * dAdVbs / a);
                    cbsb = -(cbgb + cbdb + cbbb);
                    cdgb = 0.0;
                    cddb = 0.0;
                    cdsb = 0.0;
                    goto finished;
                }
                goto finished;
            }
            else
            {
                co4v15 = 4.0 / 15.0;
                vth0 = vfb + phi + k1 * sqrtVpb;
                vgs_Vth = vgs - vth0;
                arg1 = a * vds;
                arg2 = vgs_Vth - 0.5 * arg1;
                arg3 = vds - arg1;
                arg5 = arg1 * arg1;
                dVthdVbs = -0.5 * k1 / sqrtVpb;
                dAdVbs = 0.5 * k1 * (0.5 * g / vpb - 0.8364 * (1 - g) * (1 - g)) / sqrtVpb;
                ent = Math.Max(arg2, 1.0e-8);
                dEntdVds = -0.5 * a;
                dEntdVbs = -dVthdVbs - 0.5 * vds * dAdVbs;
                vcom = vgs_Vth * vgs_Vth / 6.0 - 1.25e-1 * arg1 * vgs_Vth + 2.5e-2 * arg5;
                vdsPinchoff = Math.Max(vgs_Vth / a, 0.0);
                vgb = vgs - vbs;
                vgb_Vfb = vgb - vfb;
                if (vgb_Vfb < 0)
                {
                    qg = wLCox * vgb_Vfb;
                    qb = -qg;
                    qd = 0.0;
                    cggb = wLCox;
                    cgdb = 0.0;
                    cgsb = 0.0;
                    cbgb = -wLCox;
                    cbdb = 0.0;
                    cbsb = 0.0;
                    cdgb = 0.0;
                    cddb = 0.0;
                    cdsb = 0.0;
                    goto finished;
                }
                else if (vgs < vth0)
                {
                    qg = 0.5 * wLCox * k1 * k1 * (-1 + Math.Sqrt(1 + 4 * vgb_Vfb / (k1 * k1)));
                    qb = -qg;
                    qd = 0.0;
                    cggb = wLCox / Math.Sqrt(1 + 4 * vgb_Vfb / (k1 * k1));
                    cgdb = cgsb = 0.0;
                    cbgb = -cggb;
                    cbdb = cbsb = cdgb = cddb = cdsb = 0.0;
                    goto finished;
                }
                else if (vds < vdsPinchoff)
                {
                    vgs_VthSquare = vgs_Vth * vgs_Vth;
                    entSquare = ent * ent;
                    argl1 = 1.2e1 * entSquare;
                    argl2 = 1.0 - a;
                    argl3 = arg1 * vds;
                    argl4 = vcom / ent / entSquare;
                    if (ent > 1.0e-8)
                    {
                        argl5 = arg1 / ent;
                        argl6 = vcom / entSquare;
                    }
                    else
                    {
                        argl5 = 2.0;
                        argl6 = 4.0 / 1.5e1;
                    }
                    argl7 = argl5 / 1.2e1;
                    argl8 = 6.0 * ent;
                    qg = wLCox * (vgs - vfb - phi - 0.5 * vds + vds * argl7);
                    qb = wLCox * (-vth0 + vfb + phi + 0.5 * arg3 - arg3 * argl7);
                    qd = -wLCox * (0.5 * (vgs_Vth - arg1) + arg1 * argl6);
                    cggb = wLCox * (1.0 - argl3 / argl1);
                    cgdb = wLCox * (-0.5 + arg1 / argl8 - argl3 * dEntdVds / argl1);
                    cgbb = wLCox * (vds * vds * dAdVbs * ent - argl3 * dEntdVbs) / argl1;
                    cgsb = -(cggb + cgdb + cgbb);
                    cbgb = wLCox * argl3 * argl2 / argl1;
                    cbdb = wLCox * argl2 * (0.5 - arg1 / argl8 + argl3 * dEntdVds / argl1);
                    cbbb = -wLCox * (dVthdVbs + 0.5 * vds * dAdVbs + vds * vds * ((1.0 - 2.0 * a) * dAdVbs * ent - argl2 * a * dEntdVbs) / argl1);
                    cbsb = -(cbgb + cbdb + cbbb);
                    cdgb = -wLCox * (0.5 + arg1 * (4.0 * vgs_Vth - 1.5 * arg1) / argl1 - 2.0 * arg1 * argl4);
                    cddb = wLCox * (0.5 * a + 2.0 * arg1 * dEntdVds * argl4 - a * (2.0 * vgs_VthSquare - 3.0 * arg1 * vgs_Vth + 0.9 * arg5) / argl1);
                    cdbb = wLCox * (0.5 * dVthdVbs + 0.5 * vds * dAdVbs + 2.0 * arg1 * dEntdVbs * argl4 - vds * (2.0 * vgs_VthSquare * dAdVbs - 4.0 * a * vgs_Vth * dVthdVbs - 3.0 * arg1 * vgs_Vth * dAdVbs + 1.5 * a * arg1 * dVthdVbs + 0.9 * arg5 * dAdVbs) / argl1);
                    cdsb = -(cdgb + cddb + cdbb);
                    goto finished;
                }
                else if (vds >= vdsPinchoff)
                {
                    args1 = 1.0 / (3.0 * a);
                    qg = wLCox * (vgs - vfb - phi - vgs_Vth * args1);
                    qb = wLCox * (vfb + phi - vth0 + (1.0 - a) * vgs_Vth * args1);
                    qd = -co4v15 * wLCox * vgs_Vth;
                    cggb = wLCox * (1.0 - args1);
                    cgdb = 0.0;
                    cgbb = wLCox * args1 * (dVthdVbs + vgs_Vth * dAdVbs / a);
                    cgsb = -(cggb + cgdb + cgbb);
                    cbgb = wLCox * (args1 - 1.0 / 3.0);
                    cbdb = 0.0;
                    cbbb = -wLCox * ((2.0 / 3.0 + args1) * dVthdVbs + vgs_Vth * args1 * dAdVbs / a);
                    cbsb = -(cbgb + cbdb + cbbb);
                    cdgb = -co4v15 * wLCox;
                    cddb = 0.0;
                    cdbb = co4v15 * wLCox * dVthdVbs;
                    cdsb = -(cdgb + cddb + cdbb);
                    goto finished;
                }
            }
            finished:
            gmPointer = Math.Max(gm, 0.0);
            gdsPointer = Math.Max(gds, 0.0);
            gmbsPointer = Math.Max(gmbs, 0.0);
            qgPointer = qg;
            qbPointer = qb;
            qdPointer = qd;
            cggbPointer = cggb;
            cgdbPointer = cgdb;
            cgsbPointer = cgsb;
            cbgbPointer = cbgb;
            cbdbPointer = cbdb;
            cbsbPointer = cbsb;
            cdgbPointer = cdgb;
            cddbPointer = cddb;
            cdsbPointer = cdsb;
            cdrainPointer = Math.Max(drainCurrent, 0.0);
            vonPointer = von;
            vdsatPointer = vdsSat;
        }

        /// <summary>
        /// Helper method MosCap
        /// </summary>
        public void MosCap(double vgd, double vgs, double vgb, double[] args, double cbgb, double cbdb, double cbsb, double cdgb, double cddb, double cdsb, out double gcggbPointer, out double gcgdbPointer, out double gcgsbPointer, out double gcbgbPointer, out double gcbdbPointer, out double gcbsbPointer, out double gcdgbPointer, out double gcddbPointer, out double gcdsbPointer, out double gcsgbPointer, out double gcsdbPointer, out double gcssbPointer, ref double qGatePointer, ref double qBulkPointer, ref double qDrainPointer, out double qSourcePointer)
        {
            double ag0, qgb, qgs, qgd;
            ag0 = TranBehavior?.Qb.Jacobian(1.0) ?? 0.0; // ckt->CKTag[0];
            gcdgbPointer = (cdgb - args[0]) * ag0;
            gcddbPointer = (cddb + args[3] + args[0]) * ag0;
            gcdsbPointer = cdsb * ag0;
            gcsgbPointer = -(args[5] + cbgb + cdgb + args[1]) * ag0;
            gcsdbPointer = -(args[6] + cbdb + cddb) * ag0;
            gcssbPointer = (args[4] + args[1] - (args[7] + cbsb + cdsb)) * ag0;
            gcggbPointer = (args[5] + args[0] + args[1] + args[2]) * ag0;
            gcgdbPointer = (args[6] - args[0]) * ag0;
            gcgsbPointer = (args[7] - args[1]) * ag0;
            gcbgbPointer = (cbgb - args[2]) * ag0;
            gcbdbPointer = (cbdb - args[3]) * ag0;
            gcbsbPointer = (cbsb - args[4]) * ag0;
            qgd = args[0] * vgd;
            qgs = args[1] * vgs;
            qgb = args[2] * vgb;
            qGatePointer = qGatePointer + qgd + qgs + qgb;
            qBulkPointer = qBulkPointer - qgb;
            qDrainPointer = qDrainPointer - qgd;
            qSourcePointer = -(qGatePointer + qBulkPointer + qDrainPointer);
        }
    }
}