using System;
using SpiceSharp.Simulations;
using SpiceSharp.Behaviors;
using SpiceSharp.Components.Semiconductors;
using SpiceSharp;
using SpiceSharp.Components;
using SpiceSharp.Algebra;
using SpiceSharp.Components.Mosfets;
using SpiceSharp.Attributes;

namespace SpiceSharpBSIM.Components.Semiconductors.BSIM.BSIM1Behaviors
{
    /// <summary>
    /// Load behavior for a <see cref="BSIM1" />
    /// </summary>
    [BehaviorFor(typeof(BSIM1)), AddBehaviorIfNo(typeof(IBiasingBehavior))]
    public class BiasingBehavior : TemperatureBehavior, IBiasingBehavior, ITimeBehavior
    {
        private readonly IIterationSimulationState _iteration;
        private readonly IIntegrationMethod _method;
        private readonly ITimeSimulationState _time;
        private readonly IDerivative _qg, _qb, _qd;

        /// <summary>
        /// Gets or sets whether the small-signal parameters need to be computed.
        /// </summary>
        protected bool ComputeSmallSignal { get; set; } = false;

        /// <summary>
        /// Properties
        /// </summary>
        public double Mode { get; private set; }
        public double Vdsat { get; private set; }
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

        private readonly IVariable<double> _drain, _gate, _source, _bulk, _drainPrime, _sourcePrime;
        private readonly Element<double> _dpPtr, _spPtr, _gPtr, _bPtr;
        private readonly Element<double> _ddPtr, _ggPtr, _ssPtr, _bbPtr,
            _dpdpPtr, _spspPtr, _ddpPtr, _gbPtr, _gdpPtr, _gspPtr, _sspPtr, _bdpPtr, _bspPtr,
            _spsPtr, _dpspPtr, _dpdPtr, _bgPtr, _dpgPtr, _spgPtr, _dpbPtr, _spbPtr, _spdpPtr;

        /// <summary>
        /// Constructor
        /// </summary>
        public BiasingBehavior(ComponentBindingContext context)
            : base(context)
        {
            // Basically the setup at this point
            _iteration = context.GetState<IIterationSimulationState>();
            context.TryGetState(out _time);
            context.TryGetState(out _method);
            var state = context.GetState<IBiasingSimulationState>();
            _drain = state.GetSharedVariable(context.Nodes[0]);
            _gate = state.GetSharedVariable(context.Nodes[1]);
            _source = state.GetSharedVariable(context.Nodes[2]);
            _bulk = state.GetSharedVariable(context.Nodes[3]);

            if (!ModelParameters.SheetResistance.Value.Equals(0.0) && !Parameters.DrainSquares.Value.Equals(0.0))
                _drainPrime = state.CreatePrivateVariable(context.Behaviors.Name.Combine("drain"), Units.Volt);
            else
                _drainPrime = _drain;
            if (!ModelParameters.SheetResistance.Value.Equals(0.0) && !Parameters.SourceSquares.Value.Equals(0.0))
                _sourcePrime = state.CreatePrivateVariable(context.Behaviors.Name.Combine("source"), Units.Volt);
            else
                _sourcePrime = _source;

            var drain = state.Map[_drain];
            var gate = state.Map[_gate];
            var source = state.Map[_source];
            var bulk = state.Map[_bulk];
            var drainPrime = state.Map[_drainPrime];
            var sourcePrime = state.Map[_sourcePrime];
            _dpPtr = state.Solver.GetElement(drainPrime);
            _gPtr = state.Solver.GetElement(gate);
            _spPtr = state.Solver.GetElement(sourcePrime);
            _bPtr = state.Solver.GetElement(bulk);
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

            _qg = _method?.CreateDerivative() ?? new SimpleDerivative();
            _qd = _method?.CreateDerivative() ?? new SimpleDerivative();
            _qb = _method?.CreateDerivative() ?? new SimpleDerivative();
        }

        /// <inheritdoc />
        public void InitializeStates()
        {
            // Since we are trying to stick to the original code as closely as possible,
            // we leave this method empty and do it in the Load() method.
        }

        /// <summary>
        /// Load the behavior
        /// </summary>
        /// <param name="simulation">Simulation</param>
        void IBiasingBehavior.Load()
        {
            double capbd = 0.0;
            double capbs = 0.0;
            bool Check;
            double cdgb = 0.0;
            double cddb = 0.0;
            double cdsb = 0.0;
            double csgb = 0.0;
            double cssb = 0.0;
            double csdb = 0.0;
            double qdrn = 0.0;
            double qsrc = 0.0;
            double[] args = new double[8];

            double m; /* parallel multiplier */

            double EffectiveLength = Parameters.Length - ModelParameters.DeltaL * 1.0e-6;
            double DrainArea = Parameters.DrainArea;
            double SourceArea = Parameters.SourceArea;
            double DrainPerimeter = Parameters.DrainPerimeter;
            double SourcePerimeter = Parameters.SourcePerimeter;
            double DrainSatCurrent;
            if ((DrainSatCurrent = DrainArea * ModelParameters.JctSatCurDensity) < 1e-15)
                DrainSatCurrent = 1.0e-15;
            double SourceSatCurrent;
            if ((SourceSatCurrent = SourceArea * ModelParameters.JctSatCurDensity) < 1.0e-15)
                SourceSatCurrent = 1.0e-15;
            double GateSourceOverlapCap = ModelParameters.GateSourceOverlapCap * Parameters.Width;
            double GateDrainOverlapCap = ModelParameters.GateDrainOverlapCap * Parameters.Width;
            double GateBulkOverlapCap = ModelParameters.GateBulkOverlapCap * EffectiveLength;
            double von = ModelParameters.Type * Von;
            double vdsat = ModelParameters.Type * Vdsat;
            double vt0 = ModelParameters.Type * Vt0;

            Check = true;
            double vbd;
            double vbs;
            double vds;
            double vgd;
            double vgs;
            if (_iteration.Mode == IterationModes.Float || (_time != null && !_time.UseDc && _method != null && _method.BaseTime.Equals(0.0)) ||
                _iteration.Mode == IterationModes.Fix && !Parameters.Off)
            {
                vbs = ModelParameters.Type * (_bulk.Value - _sourcePrime.Value);
                vgs = ModelParameters.Type * (_gate.Value - _sourcePrime.Value);
                vds = ModelParameters.Type * (_drainPrime.Value - _sourcePrime.Value);

                vbd = vbs - vds;
                vgd = vgs - vds;
                double vgdo = Vgs - Vds;

                von = ModelParameters.Type * this.Von;
                if (Vds >= 0)
                {
                    vgs = Transistor.LimitFet(vgs, Vgs, von);
                    vds = vgs - vgd;
                    vds = Transistor.LimitVds(vds, Vds);
                    vgd = vgs - vds;
                }
                else
                {
                    vgd = Transistor.LimitFet(vgd, vgdo, von);
                    vds = vgs - vgd;
                    vds = -Transistor.LimitVds(-vds, -Vds);
                    vgs = vgd + vds;
                }

                double vcrit;
                Check = false;
                if (vds >= 0)
                {
                    vcrit = Constants.Vt0 * Math.Log(Constants.Vt0 / (Constants.Root2 * SourceSatCurrent));
                    vbs = Semiconductor.LimitJunction(vbs, Vbs, Constants.Vt0, vcrit, ref Check); /* B1 test */
                    vbd = vbs - vds;
                }
                else
                {
                    vcrit = Constants.Vt0 * Math.Log(Constants.Vt0 / (Constants.Root2 * DrainSatCurrent));
                    vbd = Semiconductor.LimitJunction(vbd, Vbd, Constants.Vt0, vcrit, ref Check); /* B1 test*/
                    vbs = vbd + vds;
                }
            }
            else
            {
                if (_iteration.Mode == IterationModes.Junction && !Parameters.Off)
                {
                    vds = ModelParameters.Type * Parameters.IcVDS;
                    vgs = ModelParameters.Type * Parameters.IcVGS;
                    vbs = ModelParameters.Type * Parameters.IcVBS;

                    // TODO: At some point, check what this is supposed to do
                    if (vds.Equals(0) && vgs.Equals(0) && vbs.Equals(0) && (_time == null || _time.UseDc || !_time.UseIc))
                    {
                        vbs = -1;
                        vgs = vt0;
                        vds = 0;
                    }
                }
                else
                {
                    vbs = vgs = vds = 0;
                }
            }


            /* determine DC current and derivatives */
            vbd = vbs - vds;
            vgd = vgs - vds;
            double vgb = vgs - vbs;
            double cbs;
            double gbs;
            if (vbs <= 0.0)
            {
                gbs = SourceSatCurrent / Constants.Vt0 + _iteration.Gmin;
                cbs = gbs * vbs;
            }
            else
            {
                double evbs = Math.Exp(vbs / Constants.Vt0);
                gbs = SourceSatCurrent * evbs / Constants.Vt0 + _iteration.Gmin;
                cbs = SourceSatCurrent * (evbs - 1) + _iteration.Gmin * vbs;
            }
            double cbd;
            double gbd;
            if (vbd <= 0.0)
            {
                gbd = DrainSatCurrent / Constants.Vt0 + _iteration.Gmin;
                cbd = gbd * vbd;
            }
            else
            {
                double evbd = Math.Exp(vbd / Constants.Vt0);
                gbd = DrainSatCurrent * evbd / Constants.Vt0 + _iteration.Gmin;
                cbd = DrainSatCurrent * (evbd - 1) + _iteration.Gmin * vbd;
            }
            /* line 400 */
            if (vds >= 0)
            {
                /* normal mode */
                Mode = 1;
            }
            else
            {
                /* inverse mode */
                Mode = -1;
            }

            double cdrain;
            double gds;
            double gm;
            double gmbs;
            double cgdb;
            double cgsb;
            double cbdb;
            double cggb;
            double cbgb;
            double cbsb;
            double qgate;
            double qbulk;
            /* call B1evaluate to calculate drain current and its 
             * derivatives and charge and capacitances related to gate
             * drain, and bulk
             */
            if (vds >= 0)
            {
                B1Evaluate(vds, vbs, vgs, out gm, out gds, out gmbs, out qgate,
                    out qbulk, out qdrn, out cggb, out cgdb, out cgsb, out cbgb, out cbdb, out cbsb, out cdgb,
                    out cddb, out cdsb, out cdrain, out von, out vdsat);
            }
            else
            {
                B1Evaluate(-vds, vbd, vgd, out gm, out gds, out gmbs, out qgate,
                    out qbulk, out qsrc, out cggb, out cgsb, out cgdb, out cbgb, out cbsb, out cbdb, out csgb,
                    out cssb, out csdb, out cdrain, out von, out vdsat);
            }

            this.Von = ModelParameters.Type * von;
            this.Vdsat = ModelParameters.Type * vdsat;

            /*
             *  COMPUTE EQUIVALENT DRAIN CURRENT SOURCE
             */
            double cd = Mode * cdrain - cbd;
            if (_method != null || ComputeSmallSignal)
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

                double czbd = ModelParameters.UnitAreaJctCap * DrainArea;
                double czbs = ModelParameters.UnitAreaJctCap * SourceArea;
                double czbdsw = ModelParameters.UnitLengthSidewallJctCap * DrainPerimeter;
                double czbssw = ModelParameters.UnitLengthSidewallJctCap * SourcePerimeter;
                double PhiB = ModelParameters.BulkJctPotential;
                double PhiBSW = ModelParameters.SidewallJctPotential;
                double MJ = ModelParameters.BulkJctBotGradingCoeff;
                double MJSW = ModelParameters.BulkJctSideGradingCoeff;
                double arg;
                double sarg;
                double sargsw;
                double argsw;
                /* Source Bulk Junction */
                if (vbs < 0)
                {
                    arg = 1 - vbs / PhiB;
                    argsw = 1 - vbs / PhiBSW;
                    sarg = Math.Exp(-MJ * Math.Log(arg));
                    sargsw = Math.Exp(-MJSW * Math.Log(argsw));
                    Qbs =
                        PhiB * czbs * (1 - arg * sarg) / (1 - MJ) + PhiBSW *
                    czbssw * (1 - argsw * sargsw) / (1 - MJSW);
                    capbs = czbs * sarg + czbssw * sargsw;
                }
                else
                {
                    Qbs =
                        vbs * (czbs + czbssw) + vbs * vbs * (czbs * MJ * 0.5 / PhiB
                        + czbssw * MJSW * 0.5 / PhiBSW);
                    capbs = czbs + czbssw + vbs * (czbs * MJ / PhiB +
                        czbssw * MJSW / PhiBSW);
                }

                /* Drain Bulk Junction */
                if (vbd < 0)
                {
                    arg = 1 - vbd / PhiB;
                    argsw = 1 - vbd / PhiBSW;
                    sarg = Math.Exp(-MJ * Math.Log(arg));
                    sargsw = Math.Exp(-MJSW * Math.Log(argsw));
                    Qbd =
                        PhiB * czbd * (1 - arg * sarg) / (1 - MJ) + PhiBSW *
                    czbdsw * (1 - argsw * sargsw) / (1 - MJSW);
                    capbd = czbd * sarg + czbdsw * sargsw;
                }
                else
                {
                    Qbd =
                        vbd * (czbd + czbdsw) + vbd * vbd * (czbd * MJ * 0.5 / PhiB
                        + czbdsw * MJSW * 0.5 / PhiBSW);
                    capbd = czbd + czbdsw + vbd * (czbd * MJ / PhiB +
                        czbdsw * MJSW / PhiBSW);
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
            Vbs = vbs;
            Vbs = vbs;
            Vbd = vbd;
            Vgs = vgs;
            Vds = vds;
            Cd = cd;
            Cbs = cbs;
            Cbd = cbd;
            Gm = gm;
            Gds = gds;
            Gmbs = gmbs;
            Gbd = gbd;
            Gbs = gbs;

            Cggb = cggb;
            Cgdb = cgdb;
            Cgsb = cgsb;

            Cbgb = cbgb;
            Cbdb = cbdb;
            Cbsb = cbsb;

            Cdgb = cdgb;
            Cddb = cddb;
            Cdsb = cdsb;

            Capbs = capbs;
            Capbd = capbd;

            /* bulk and channel charge plus overlaps */
            if (_method == null && (_time != null && !_time.UseIc))
                goto line850;
            double gcbdb;
            double gcbgb;
            double gcbsb;
            double gcddb;
            double gcdgb;
            double gcdsb;
            double gcgdb;
            double gcggb;
            double gcgsb;
            double gcsdb;
            double gcsgb;
            double gcssb;
            if (Mode > 0)
            {

                args[0] = GateDrainOverlapCap;
                args[1] = GateSourceOverlapCap;
                args[2] = GateBulkOverlapCap;
                args[3] = capbd;
                args[4] = capbs;
                args[5] = cggb;
                args[6] = cgdb;
                args[7] = cgsb;

                B1mosCap(vgd, vgs, vgb,
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
                    ref qgate, ref qbulk, ref qdrn, ref qsrc);
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

                B1mosCap(vgs, vgd, vgb,
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
                    ref qgate, ref qbulk, ref qsrc, ref qdrn);
            }

            _qg.Value = qgate;
            _qd.Value = qdrn - Qbd;
            _qb.Value = qbulk + Qbd + Qbs;

            /* store small signal parameters */
            if (_method == null)
                goto line850;

            _qb.Derive();
            _qg.Derive();
            _qd.Derive();
            goto line860;
        line850:
            double ceqqb;
            double ceqqd;
            double ceqqg;
            /* initialize to zero charge conductance and current */
            ceqqg = ceqqb = ceqqd = 0.0;
            gcdgb = gcddb = gcdsb = 0.0;
            gcsgb = gcsdb = gcssb = 0.0;
            gcggb = gcgdb = gcgsb = 0.0;
            gcbgb = gcbdb = gcbsb = 0.0;
            goto line900;

        line860:
            /* evaluate equivalent charge current */
            double cqgate = _qg.Derivative;
            double cqbulk = _qb.Derivative;
            double cqdrn = _qd.Derivative;
            ceqqg = cqgate - gcggb * vgb + gcgdb * vbd + gcgsb * vbs;
            ceqqb = cqbulk - gcbgb * vgb + gcbdb * vbd + gcbsb * vbs;
            ceqqd = cqdrn - gcdgb * vgb + gcddb * vbd + gcdsb * vbs;

        /*
         *  load current vector
         */
        line900:
            m = Parameters.Multiplier;
            double ceqbs = ModelParameters.Type * (cbs - (gbs - _iteration.Gmin) * vbs);
            double ceqbd = ModelParameters.Type * (cbd - (gbd - _iteration.Gmin) * vbd);

            ceqqg = ModelParameters.Type * ceqqg;
            ceqqb = ModelParameters.Type * ceqqb;
            ceqqd = ModelParameters.Type * ceqqd;
            double cdreq;
            double xnrm;
            double xrev;
            if (Mode >= 0)
            {
                xnrm = 1;
                xrev = 0;
                cdreq = ModelParameters.Type * (cdrain - gds * vds - gm * vgs - gmbs * vbs);
            }
            else
            {
                xnrm = 0;
                xrev = 1;
                cdreq = -(ModelParameters.Type) * (cdrain + gds * vds - gm * vgd - gmbs * vbd);
            }

            _gPtr.Value -= m * ceqqg;
            _bPtr.Value -= m * (ceqbs + ceqbd + ceqqb);
            _dpPtr.Value += m * (ceqbd - cdreq - ceqqd);
            _spPtr.Value += m * (cdreq + ceqbs + ceqqg + ceqqb + ceqqd);

            /*
             *  load y matrix
             */
            _ddPtr.Value += m * (this.DrainConductance);
            _ggPtr.Value += m * (gcggb);
            _ssPtr.Value += m * (this.SourceConductance);
            _bbPtr.Value += m * (gbd + gbs - gcbgb - gcbdb - gcbsb);
            _dpdpPtr.Value += m * (this.DrainConductance + gds + gbd + xrev * (gm + gmbs) + gcddb);
            _spspPtr.Value += m * (this.SourceConductance + gds + gbs + xnrm * (gm + gmbs) + gcssb);
            _ddpPtr.Value += m * (-this.DrainConductance);
            _gbPtr.Value += m * (-gcggb - gcgdb - gcgsb);
            _gdpPtr.Value += m * (gcgdb);
            _gspPtr.Value += m * (gcgsb);
            _sspPtr.Value += m * (-this.SourceConductance);
            _bgPtr.Value += m * (gcbgb);
            _bdpPtr.Value += m * (-gbd + gcbdb);
            _bspPtr.Value += m * (-gbs + gcbsb);
            _dpdPtr.Value += m * (-this.DrainConductance);
            _dpgPtr.Value += m * ((xnrm - xrev) * gm + gcdgb);
            _dpbPtr.Value += m * (-gbd + (xnrm - xrev) * gmbs - gcdgb - gcddb - gcdsb);
            _dpspPtr.Value += m * (-gds - xnrm * (gm + gmbs) + gcdsb);
            _spgPtr.Value += m * (-(xnrm - xrev) * gm + gcsgb);
            _spsPtr.Value += m * (-this.SourceConductance);
            _spbPtr.Value += m * (-gbs - (xnrm - xrev) * gmbs - gcsgb - gcsdb - gcssb);
            _spdpPtr.Value += m * (-gds - xrev * (gm + gmbs) + gcsdb);
        }

        protected void B1Evaluate(double vds, double vbs, double vgs,
           out double gmPointer, out double gdsPointer, out double gmbsPointer,
           out double qgPointer, out double qbPointer, out double qdPointer,
           out double cggbPointer, out double cgdbPointer, out double cgsbPointer,
           out double cbgbPointer, out double cbdbPointer, out double cbsbPointer,
           out double cdgbPointer, out double cddbPointer, out double cdsbPointer,
           out double cdrainPointer, out double vonPointer, out double vdsatPointer)
        {
            double gm;
            double gds;
            double gmbs;
            double qg = 0.0;
            double qb = 0.0;
            double qd = 0.0;
            double cggb = 0.0;
            double cgdb = 0.0;
            double cgsb = 0.0;
            double cbgb = 0.0;
            double cbdb = 0.0;
            double cbsb = 0.0;
            double cdgb = 0.0;
            double cddb = 0.0;
            double cdsb = 0.0;
            double Vfb;
            double Phi;
            double K1;
            double K2;
            double Vdd;
            double Ugs;
            double Uds;
            double dUgsdVbs;
            double Leff;
            double dUdsdVbs;
            double dUdsdVds;
            double Eta;
            double dEtadVds;
            double dEtadVbs;
            double Vpb;
            double SqrtVpb;
            double Von;
            double Vth;
            double dVthdVbs;
            double dVthdVds;
            double Vgs_Vth;
            double DrainCurrent;
            double G;
            double A;
            double Arg;
            double dGdVbs;
            double dAdVbs;
            double Beta;
            double Beta_Vds_0;
            double BetaVdd;
            double dBetaVdd_dVds;
            double Beta0;
            double dBeta0dVds;
            double dBeta0dVbs;
            double VddSquare;
            double C1;
            double C2;
            double dBetaVdd_dVbs;
            double dBeta_Vds_0_dVbs;
            double dC1dVbs;
            double dC2dVbs;
            double dBetadVgs;
            double dBetadVds;
            double dBetadVbs;
            double VdsSat = 0.0;
            double Argl1;
            double Argl2;
            double Vc;
            double Term1;
            double K;
            double Args1;
            double dVcdVgs;
            double dVcdVds;
            double dVcdVbs;
            double dKdVc;
            double dKdVgs;
            double dKdVds;
            double dKdVbs;
            double Args2;
            double Args3;
            double Warg1;
            double Vcut;
            double N;
            double N0;
            double NB;
            double ND;
            double Warg2;
            double Wds;
            double Wgs;
            double Ilimit;
            double Iexp;
            double Temp1;
            double Vth0;
            double Arg1;
            double Arg2;
            double Arg3;
            double Arg5;
            double Ent;
            double Vcom;
            double Vgb;
            double Vgb_Vfb;
            double VdsPinchoff;
            double EntSquare;
            double Vgs_VthSquare;
            double Argl3;
            double Argl4;
            double Argl5;
            double Argl6;
            double Argl7;
            double Argl8;
            double Argl9;
            double dEntdVds;
            double dEntdVbs;
            double cgbb;
            double cdbb;
            double cbbb;
            double WLCox;
            double Vtsquare;
            double Temp3;
            bool ChargeComputationNeeded;
            double co4v15;

            if (_method != null || ComputeSmallSignal)
                ChargeComputationNeeded = true;
            else
                ChargeComputationNeeded = false;

            Vfb = this.Vfb;
            Phi = this.Phi;
            K1 = this.K1;
            K2 = this.K2;
            Vdd = ModelParameters.Vdd;
            if ((Ugs = this.Ugs + this.UgsB * vbs) <= 0)
            {
                Ugs = 0;
                dUgsdVbs = 0.0;
            }
            else
            {
                dUgsdVbs = this.UgsB;
            }
            if ((Uds = this.Uds + this.UdsB * vbs +
                    this.UdsD * (vds - Vdd)) <= 0)
            {
                Uds = 0.0;
                dUdsdVbs = dUdsdVds = 0.0;
            }
            else
            {
                Leff = Parameters.Length * 1.0e6 - ModelParameters.DeltaL; /* Leff in um */
                Uds = Uds / Leff;
                dUdsdVbs = this.UdsB / Leff;
                dUdsdVds = this.UdsD / Leff;
            }
            Eta = this.Eta + this.EtaB * vbs + this.EtaD *
                (vds - Vdd);
            if (Eta <= 0)
            {
                Eta = 0;
                dEtadVds = dEtadVbs = 0.0;
            }
            else if (Eta > 1)
            {
                Eta = 1;
                dEtadVds = dEtadVbs = 0;
            }
            else
            {
                dEtadVds = this.EtaD;
                dEtadVbs = this.EtaB;
            }
            if (vbs < 0)
            {
                Vpb = Phi - vbs;
            }
            else
            {
                Vpb = Phi;
            }
            SqrtVpb = Math.Sqrt(Vpb);
            Von = Vfb + Phi + K1 * SqrtVpb - K2 * Vpb - Eta * vds;
            Vth = Von;
            dVthdVds = -Eta - dEtadVds * vds;
            dVthdVbs = K2 - 0.5 * K1 / SqrtVpb - dEtadVbs * vds;
            Vgs_Vth = vgs - Vth;

            G = 1.0 - 1.0 / (1.744 + 0.8364 * Vpb);
            A = 1.0 + 0.5 * G * K1 / SqrtVpb;
            A = Math.Max(A, 1.0);   /* Modified */
            Arg = Math.Max((1 + Ugs * Vgs_Vth), 1.0);
            dGdVbs = -0.8364 * (1 - G) * (1 - G);
            dAdVbs = 0.25 * K1 / SqrtVpb * (2 * dGdVbs + G / Vpb);

            if (Vgs_Vth < 0)
            {
                /* cutoff */
                DrainCurrent = 0;
                gm = 0;
                gds = 0;
                gmbs = 0;
                goto SubthresholdComputation;
            }

            /* Quadratic Interpolation for Beta0 (Beta at vgs  =  0, vds=Vds) */

            Beta_Vds_0 = (this.BetaZero + this.BetaZeroB * vbs);
            BetaVdd = (this.BetaVdd + this.BetaVddB * vbs);
            dBetaVdd_dVds = Math.Max(this.BetaVddD, 0.0); /* Modified */
            if (vds > Vdd)
            {
                Beta0 = BetaVdd + dBetaVdd_dVds * (vds - Vdd);
                dBeta0dVds = dBetaVdd_dVds;
                dBeta0dVbs = this.BetaVddB;
            }
            else
            {
                VddSquare = Vdd * Vdd;
                C1 = (-BetaVdd + Beta_Vds_0 + dBetaVdd_dVds * Vdd) / VddSquare;
                C2 = 2 * (BetaVdd - Beta_Vds_0) / Vdd - dBetaVdd_dVds;
                dBeta_Vds_0_dVbs = this.BetaZeroB;
                dBetaVdd_dVbs = this.BetaVddB;
                dC1dVbs = (dBeta_Vds_0_dVbs - dBetaVdd_dVbs) / VddSquare;
                dC2dVbs = dC1dVbs * (-2) * Vdd;
                Beta0 = (C1 * vds + C2) * vds + Beta_Vds_0;
                dBeta0dVds = 2 * C1 * vds + C2;
                dBeta0dVbs = dC1dVbs * vds * vds + dC2dVbs * vds + dBeta_Vds_0_dVbs;
            }

            /*Beta  =  Beta0 / ( 1 + Ugs * Vgs_Vth );*/

            Beta = Beta0 / Arg;
            dBetadVgs = -Beta * Ugs / Arg;
            dBetadVds = dBeta0dVds / Arg - dBetadVgs * dVthdVds;
            dBetadVbs = dBeta0dVbs / Arg + Beta * Ugs * dVthdVbs / Arg -
                Beta * Vgs_Vth * dUgsdVbs / Arg;

            /*VdsSat  = Math.Max( Vgs_Vth / ( A + Uds * Vgs_Vth ),  0.0);*/

            if ((Vc = Uds * Vgs_Vth / A) < 0.0) Vc = 0.0;
            Term1 = Math.Sqrt(1 + 2 * Vc);
            K = 0.5 * (1 + Vc + Term1);
            VdsSat = Math.Max(Vgs_Vth / (A * Math.Sqrt(K)), 0.0);

            if (vds < VdsSat)
            {
                /* Triode Region */
                /*Argl1  =  1 + Uds * vds;*/
                Argl1 = Math.Max((1 + Uds * vds), 1.0);
                Argl2 = Vgs_Vth - 0.5 * A * vds;
                DrainCurrent = Beta * Argl2 * vds / Argl1;
                gm = (dBetadVgs * Argl2 * vds + Beta * vds) / Argl1;
                gds = (dBetadVds * Argl2 * vds + Beta *
                    (Vgs_Vth - vds * dVthdVds - A * vds) -
                    DrainCurrent * (vds * dUdsdVds + Uds)) / Argl1;
                gmbs = (dBetadVbs * Argl2 * vds + Beta * vds *
                    (-dVthdVbs - 0.5 * vds * dAdVbs) -
                    DrainCurrent * vds * dUdsdVbs) / Argl1;
            }
            else
            {
                /* Pinchoff (Saturation) Region */
                Args1 = 1.0 + 1.0 / Term1;
                dVcdVgs = Uds / A;
                dVcdVds = Vgs_Vth * dUdsdVds / A - dVcdVgs * dVthdVds;
                dVcdVbs = (Vgs_Vth * dUdsdVbs - Uds *
                    (dVthdVbs + Vgs_Vth * dAdVbs / A)) / A;
                dKdVc = 0.5 * Args1;
                dKdVgs = dKdVc * dVcdVgs;
                dKdVds = dKdVc * dVcdVds;
                dKdVbs = dKdVc * dVcdVbs;
                Args2 = Vgs_Vth / A / K;
                Args3 = Args2 * Vgs_Vth;
                DrainCurrent = 0.5 * Beta * Args3;
                gm = 0.5 * Args3 * dBetadVgs + Beta * Args2 -
                    DrainCurrent * dKdVgs / K;
                gds = 0.5 * Args3 * dBetadVds - Beta * Args2 * dVthdVds -
                    DrainCurrent * dKdVds / K;
                gmbs = 0.5 * dBetadVbs * Args3 - Beta * Args2 * dVthdVbs -
                    DrainCurrent * (dAdVbs / A + dKdVbs / K);
            }

        SubthresholdComputation:

            N0 = this.SubthSlope;
            Vcut = -40.0 * N0 * Constants.Vt0;

            /* The following 'if' statement has been modified so that subthreshold  *
             * current computation is always executed unless N0 >= 200. This should *
             * get rid of the Ids kink seen on Ids-Vgs plots at low Vds.            *
             *                                                Peter M. Lee          *
             *                                                6/8/90                *
             *  Old 'if' statement:                                                 *
             *  if( (N0 >=  200) || (Vgs_Vth < Vcut ) || (Vgs_Vth > (-0.5*Vcut)))   */

            if (N0 >= 200)
            {
                goto ChargeComputation;
            }

            NB = this.SubthSlopeB;
            ND = this.SubthSlopeD;
            N = N0 + NB * vbs + ND * vds; /* subthreshold slope */
            if (N < 0.5) N = 0.5;
            Warg1 = Math.Exp(-vds / Constants.Vt0);
            Wds = 1 - Warg1;
            Wgs = Math.Exp(Vgs_Vth / (N * Constants.Vt0));
            Vtsquare = Constants.Vt0 * Constants.Vt0;
            Warg2 = 6.04965 * Vtsquare * this.BetaZero;
            Ilimit = 4.5 * Vtsquare * this.BetaZero;
            Iexp = Warg2 * Wgs * Wds;
            DrainCurrent = DrainCurrent + Ilimit * Iexp / (Ilimit + Iexp);
            Temp1 = Ilimit / (Ilimit + Iexp);
            Temp1 = Temp1 * Temp1;
            Temp3 = Ilimit / (Ilimit + Wgs * Warg2);
            Temp3 = Temp3 * Temp3 * Warg2 * Wgs;
            /*    if ( Temp3 > Ilimit ) Temp3=Ilimit;*/
            gm = gm + Temp1 * Iexp / (N * Constants.Vt0);
            /* gds term has been modified to prevent blow up at Vds=0 */
            gds = gds + Temp3 * (-Wds / N / Constants.Vt0 * (dVthdVds +
                Vgs_Vth * ND / N) + Warg1 / Constants.Vt0);
            gmbs = gmbs - Temp1 * Iexp * (dVthdVbs + Vgs_Vth * NB / N) /
                (N * Constants.Vt0);

        ChargeComputation:

            /* Some Limiting of DC Parameters */
            if (DrainCurrent < 0.0) DrainCurrent = 0.0;
            if (gm < 0.0) gm = 0.0;
            if (gds < 0.0) gds = 0.0;
            if (gmbs < 0.0) gmbs = 0.0;

            WLCox = ModelTemperature.Cox *
                (Parameters.Length - ModelParameters.DeltaL * 1.0e-6) *
                (Parameters.Width - ModelParameters.DeltaW * 1.0e-6) * 1.0e4;   /* F */

            if (!ChargeComputationNeeded)
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
            G = 1.0 - 1.0 / (1.744 + 0.8364 * Vpb);
            A = 1.0 + 0.5 * G * K1 / SqrtVpb;
            A = Math.Max(A, 1.0);   /* Modified */
            /*Arg  =  1 + Ugs * Vgs_Vth;*/
            dGdVbs = -0.8364 * (1 - G) * (1 - G);
            dAdVbs = 0.25 * K1 / SqrtVpb * (2 * dGdVbs + G / Vpb);
            Phi = Math.Max(0.1, Phi);

            if (!ModelParameters.ChannelChargePartitionFlag.Equals(0.0))
            {

                /*0/100 partitioning for drain/source chArges at the saturation region*/
                Vth0 = Vfb + Phi + K1 * SqrtVpb;
                Vgs_Vth = vgs - Vth0;
                Arg1 = A * vds;
                Arg2 = Vgs_Vth - 0.5 * Arg1;
                Arg3 = vds - Arg1;
                Arg5 = Arg1 * Arg1;
                dVthdVbs = -0.5 * K1 / SqrtVpb;
                dAdVbs = 0.5 * K1 * (0.5 * G / Vpb - 0.8364 * (1 - G) * (1 - G)) /
                    SqrtVpb;
                Ent = Math.Max(Arg2, 1.0e-8);
                dEntdVds = -0.5 * A;
                dEntdVbs = -dVthdVbs - 0.5 * vds * dAdVbs;
                Vcom = Vgs_Vth * Vgs_Vth / 6.0 - 1.25e-1 * Arg1 *
                    Vgs_Vth + 2.5e-2 * Arg5;
                VdsPinchoff = Math.Max(Vgs_Vth / A, 0.0);
                Vgb = vgs - vbs;
                Vgb_Vfb = Vgb - Vfb;

                if (Vgb_Vfb < 0)
                {
                    /* Accumulation Region */
                    qg = WLCox * Vgb_Vfb;
                    qb = -qg;
                    qd = 0.0;
                    cggb = WLCox;
                    cgdb = 0.0;
                    cgsb = 0.0;
                    cbgb = -WLCox;
                    cbdb = 0.0;
                    cbsb = 0.0;
                    cdgb = 0.0;
                    cddb = 0.0;
                    cdsb = 0.0;
                    goto finished;
                }
                else if (vgs < Vth0)
                {
                    /* Subthreshold Region */
                    qg = 0.5 * WLCox * K1 * K1 * (-1 +
                        Math.Sqrt(1 + 4 * Vgb_Vfb / (K1 * K1)));
                    qb = -qg;
                    qd = 0.0;
                    cggb = WLCox / Math.Sqrt(1 + 4 * Vgb_Vfb / (K1 * K1));
                    cgdb = cgsb = 0.0;
                    cbgb = -cggb;
                    cbdb = cbsb = cdgb = cddb = cdsb = 0.0;
                    goto finished;
                }
                else if (vds < VdsPinchoff)
                {    /* triode region  */
                    /*Vgs_Vth2 = Vgs_Vth*Vgs_Vth;*/
                    EntSquare = Ent * Ent;
                    Argl1 = 1.2e1 * EntSquare;
                    Argl2 = 1.0 - A;
                    Argl3 = Arg1 * vds;
                    /*Argl4 = Vcom/Ent/EntSquare;*/
                    if (Ent > 1.0e-8)
                    {
                        Argl5 = Arg1 / Ent;
                        /*Argl6 = Vcom/EntSquare;*/
                    }
                    else
                    {
                        Argl5 = 2.0;
                    }
                    Argl7 = Argl5 / 1.2e1;
                    Argl8 = 6.0 * Ent;
                    Argl9 = 0.125 * Argl5 * Argl5;
                    qg = WLCox * (vgs - Vfb - Phi - 0.5 * vds + vds * Argl7);
                    qb = WLCox * (-Vth0 + Vfb + Phi + 0.5 * Arg3 - Arg3 * Argl7);
                    qd = -WLCox * (0.5 * Vgs_Vth - 0.75 * Arg1 +
                        0.125 * Arg1 * Argl5);
                    cggb = WLCox * (1.0 - Argl3 / Argl1);
                    cgdb = WLCox * (-0.5 + Arg1 / Argl8 - Argl3 * dEntdVds /
                        Argl1);
                    cgbb = WLCox * (vds * vds * dAdVbs * Ent - Argl3 * dEntdVbs) /
                        Argl1;
                    cgsb = -(cggb + cgdb + cgbb);
                    cbgb = WLCox * Argl3 * Argl2 / Argl1;
                    cbdb = WLCox * Argl2 * (0.5 - Arg1 / Argl8 + Argl3 * dEntdVds /
                        Argl1);
                    cbbb = -WLCox * (dVthdVbs + 0.5 * vds * dAdVbs + vds *
                        vds * ((1.0 - 2.0 * A) * dAdVbs * Ent - Argl2 *
                        A * dEntdVbs) / Argl1);
                    cbsb = -(cbgb + cbdb + cbbb);
                    cdgb = -WLCox * (0.5 - Argl9);
                    cddb = WLCox * (0.75 * A - 0.25 * A * Arg1 / Ent +
                        Argl9 * dEntdVds);
                    cdbb = WLCox * (0.5 * dVthdVbs + vds * dAdVbs *
                        (0.75 - 0.25 * Argl5) + Argl9 * dEntdVbs);
                    cdsb = -(cdgb + cddb + cdbb);
                    goto finished;
                }
                else if (vds >= VdsPinchoff)
                {    /* saturation region   */
                    Args1 = 1.0 / (3.0 * A);
                    qg = WLCox * (vgs - Vfb - Phi - Vgs_Vth * Args1);
                    qb = WLCox * (Vfb + Phi - Vth0 + (1.0 - A) * Vgs_Vth * Args1);
                    qd = 0.0;
                    cggb = WLCox * (1.0 - Args1);
                    cgdb = 0.0;
                    cgbb = WLCox * Args1 * (dVthdVbs + Vgs_Vth * dAdVbs / A);
                    cgsb = -(cggb + cgdb + cgbb);
                    cbgb = WLCox * (Args1 - 1.0 / 3.0);
                    cbdb = 0.0;
                    cbbb = -WLCox * ((2.0 / 3.0 + Args1) * dVthdVbs +
                        Vgs_Vth * Args1 * dAdVbs / A);      /* Modified */
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
                /* ChannelChargePartionFlag  < = 0 */

                /*40/60 partitioning for drain/source chArges at the saturation region*/
                co4v15 = 4.0/ 15.0;
                Vth0 = Vfb + Phi + K1 * SqrtVpb;
                Vgs_Vth = vgs - Vth0;
                Arg1 = A * vds;
                Arg2 = Vgs_Vth - 0.5 * Arg1;
                Arg3 = vds - Arg1;
                Arg5 = Arg1 * Arg1;
                dVthdVbs = -0.5 * K1 / SqrtVpb;
                dAdVbs = 0.5 * K1 * (0.5 * G / Vpb - 0.8364 * (1 - G) * (1 - G)) / SqrtVpb;
                Ent = Math.Max(Arg2, 1.0e-8);
                dEntdVds = -0.5 * A;
                dEntdVbs = -dVthdVbs - 0.5 * vds * dAdVbs;
                Vcom = Vgs_Vth * Vgs_Vth / 6.0 - 1.25e-1 * Arg1 * Vgs_Vth + 2.5e-2 * Arg5;
                VdsPinchoff = Math.Max(Vgs_Vth / A, 0.0);
                Vgb = vgs - vbs;
                Vgb_Vfb = Vgb - Vfb;

                if (Vgb_Vfb < 0)
                {           /* Accumulation Region */
                    qg = WLCox * Vgb_Vfb;
                    qb = -qg;
                    qd = 0.0;
                    cggb = WLCox;
                    cgdb = 0.0;
                    cgsb = 0.0;
                    cbgb = -WLCox;
                    cbdb = 0.0;
                    cbsb = 0.0;
                    cdgb = 0.0;
                    cddb = 0.0;
                    cdsb = 0.0;
                    goto finished;
                }
                else if (vgs < Vth0)
                {    /* Subthreshold Region */
                    qg = 0.5 * WLCox * K1 * K1 * (-1 + Math.Sqrt(1 + 4 * Vgb_Vfb / (K1 * K1)));
                    qb = -qg;
                    qd = 0.0;
                    cggb = WLCox / Math.Sqrt(1 + 4 * Vgb_Vfb / (K1 * K1));
                    cgdb = cgsb = 0.0;
                    cbgb = -cggb;
                    cbdb = cbsb = cdgb = cddb = cdsb = 0.0;
                    goto finished;
                }
                else if (vds < VdsPinchoff)
                {      /* triode region */

                    Vgs_VthSquare = Vgs_Vth * Vgs_Vth;
                    EntSquare = Ent * Ent;
                    Argl1 = 1.2e1 * EntSquare;
                    Argl2 = 1.0 - A;
                    Argl3 = Arg1 * vds;
                    Argl4 = Vcom / Ent / EntSquare;
                    if (Ent > 1.0e-8)
                    {
                        Argl5 = Arg1 / Ent;
                        Argl6 = Vcom / EntSquare;
                    }
                    else
                    {
                        Argl5 = 2.0;
                        Argl6 = 4.0 / 1.5e1;
                    }
                    Argl7 = Argl5 / 1.2e1;
                    Argl8 = 6.0 * Ent;
                    qg = WLCox * (vgs - Vfb - Phi - 0.5 * vds + vds * Argl7);
                    qb = WLCox * (-Vth0 + Vfb + Phi + 0.5 * Arg3 - Arg3 * Argl7);
                    qd = -WLCox * (0.5 * (Vgs_Vth - Arg1) + Arg1 * Argl6);
                    cggb = WLCox * (1.0 - Argl3 / Argl1);
                    cgdb = WLCox * (-0.5 + Arg1 / Argl8 - Argl3 * dEntdVds / Argl1);
                    cgbb = WLCox * (vds * vds * dAdVbs * Ent - Argl3 * dEntdVbs) / Argl1;
                    cgsb = -(cggb + cgdb + cgbb);
                    cbgb = WLCox * Argl3 * Argl2 / Argl1;
                    cbdb = WLCox * Argl2 * (0.5 - Arg1 / Argl8 + Argl3 * dEntdVds / Argl1);
                    cbbb = -WLCox * (dVthdVbs + 0.5 * vds * dAdVbs + vds * vds * ((1.0 - 2.0 * A)
                        * dAdVbs * Ent - Argl2 * A * dEntdVbs) / Argl1);
                    cbsb = -(cbgb + cbdb + cbbb);
                    cdgb = -WLCox * (0.5 + Arg1 * (4.0 * Vgs_Vth - 1.5 * Arg1) / Argl1 -
                        2.0 * Arg1 * Argl4);
                    cddb = WLCox * (0.5 * A + 2.0 * Arg1 * dEntdVds * Argl4 - A * (2.0 * Vgs_VthSquare
                        - 3.0 * Arg1 * Vgs_Vth + 0.9 * Arg5) / Argl1);
                    cdbb = WLCox * (0.5 * dVthdVbs + 0.5 * vds * dAdVbs + 2.0 * Arg1 * dEntdVbs
                        * Argl4 - vds * (2.0 * Vgs_VthSquare * dAdVbs - 4.0 * A * Vgs_Vth * dVthdVbs - 3.0
                        * Arg1 * Vgs_Vth * dAdVbs + 1.5 * A * Arg1 * dVthdVbs + 0.9 * Arg5 * dAdVbs)
                        / Argl1);
                    cdsb = -(cdgb + cddb + cdbb);
                    goto finished;
                }
                else if (vds >= VdsPinchoff)
                {      /* saturation region */

                    Args1 = 1.0 / (3.0 * A);
                    qg = WLCox * (vgs - Vfb - Phi - Vgs_Vth * Args1);
                    qb = WLCox * (Vfb + Phi - Vth0 + (1.0 - A) * Vgs_Vth * Args1);
                    qd = -co4v15 * WLCox * Vgs_Vth;
                    cggb = WLCox * (1.0 - Args1);
                    cgdb = 0.0;
                    cgbb = WLCox * Args1 * (dVthdVbs + Vgs_Vth * dAdVbs / A);
                    cgsb = -(cggb + cgdb + cgbb);
                    cbgb = WLCox * (Args1 - 1.0 / 3.0);
                    cbdb = 0.0;
                    cbbb = -WLCox * ((2.0 / 3.0 + Args1) * dVthdVbs + Vgs_Vth * Args1 * dAdVbs / A);
                    cbsb = -(cbgb + cbdb + cbbb);
                    cdgb = -co4v15 * WLCox;
                    cddb = 0.0;
                    cdbb = co4v15 * WLCox * dVthdVbs;
                    cdsb = -(cdgb + cddb + cdbb);
                    goto finished;
                }
            }

        finished:       /* returning Values to Calling Routine */

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
            cdrainPointer = Math.Max(DrainCurrent, 0.0);
            vonPointer = Von;
            vdsatPointer = VdsSat;
        }

        protected void B1mosCap(double vgd,
            double vgs,
            double vgb,
            double[] args,
            /*
           double GateDrainOverlapCap,
           double GateSourceOverlapCap,
           double GateBulkOverlapCap,
           double capbd,
           double capbs,
           double cggb,
           double cgdb,
           double cgsb,
            */
            double cbgb,
            double cbdb,
            double cbsb,
            double cdgb,
            double cddb,
            double cdsb,
            out double gcggbPointer,
            out double gcgdbPointer,
            out double gcgsbPointer,
            out double gcbgbPointer,
            out double gcbdbPointer,
            out double gcbsbPointer,
            out double gcdgbPointer,
            out double gcddbPointer,
            out double gcdsbPointer,
            out double gcsgbPointer,
            out double gcsdbPointer,
            out double gcssbPointer,
            ref double qGatePointer,
            ref double qBulkPointer,
            ref double qDrainPointer,
            ref double qSourcePointer)
        {
            double qgd;
            double qgs;
            double qgb;
            double ag0;

            ag0 = _method?.Slope ?? 0.0;

            /* compute equivalent conductance */
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

            /* compute total terminal charge */
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