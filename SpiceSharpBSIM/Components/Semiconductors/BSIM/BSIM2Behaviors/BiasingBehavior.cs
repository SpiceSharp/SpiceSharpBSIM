using System;
using SpiceSharp.Simulations;
using SpiceSharp.Behaviors;
using SpiceSharp.Algebra;
using SpiceSharp.Components.Semiconductors;
using SpiceSharp;
using SpiceSharp.Components;
using SpiceSharp.Components.Mosfets;
using SpiceSharp.Attributes;

namespace SpiceSharpBSIM.Components.Semiconductors.BSIM.BSIM2Behaviors
{
    /// <summary>
    /// Load behavior for a <see cref="BSIM2" />
    /// </summary>
    [BehaviorFor(typeof(BSIM2)), AddBehaviorIfNo(typeof(IBiasingBehavior))]
    public class BiasingBehavior : TemperatureBehavior, IBiasingBehavior, ITimeBehavior
    {
        private readonly IIntegrationMethod _method;
        private readonly IIterationSimulationState _iteration;
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
        public double Cggb { get; set; }
        public double Cgdb { get; set; }
        public double Cgsb { get; set; }
        public double Cbgb { get; set; }
        public double Cbdb { get; set; }
        public double Cbsb { get; set; }
        public double Capbd { get; set; }
        public double Iqbd { get; set; }
        public double Capbs { get; set; }
        public double Iqbs { get; set; }
        public double Cdgb { get; set; }
        public double Cddb { get; set; }
        public double Cdsb { get; set; }
        public double Qbs { get; set; }
        public double Qbd { get; set; }

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
            ComputeSmallSignal = true;

            // This doesn't do much except calculate the small-signal stuff again
            Load();
        }

        /// <summary>
        /// Temperature behavior
        /// </summary>
        public void Load()
        {
            double DrainSatCurrent;
            double EffectiveLength;
            double GateBulkOverlapCap;
            double GateDrainOverlapCap;
            double GateSourceOverlapCap;
            double SourceSatCurrent;
            double DrainArea;
            double SourceArea;
            double DrainPerimeter;
            double SourcePerimeter;
            double arg;
            double capbd = 0.0;
            double capbs = 0.0;
            double cbd;
            double cbs;
            double cd;
            double cdrain;
            double cdreq;
            double ceqbd;
            double ceqbs;
            double ceqqb;
            double ceqqd;
            double ceqqg;
            double czbd;
            double czbdsw;
            double czbs;
            double czbssw;
            double evbd;
            double evbs;
            double gbd;
            double gbs;
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
            double gds;
            double gm;
            double gmbs;
            double sarg;
            double sargsw;
            double vbd;
            double vbs;
            double vcrit;
            double vds;
            double vdsat;
            double vgb;
            double vgd;
            double vgdo;
            double vgs;
            double von;
            double xnrm;
            double xrev;
            bool Check;
            double cgdb;
            double cgsb;
            double cbdb;
            double cdgb = 0;
            double cddb = 0;
            double cdsb = 0;
            double cggb = 0;
            double cbgb = 0;
            double cbsb = 0;
            double csgb = 0;
            double cssb = 0;
            double csdb = 0;
            double PhiB;
            double PhiBSW;
            double MJ;
            double MJSW;
            double argsw;
            double qgate;
            double qbulk;
            double qdrn = 0;
            double qsrc = 0;
            double cqgate;
            double cqbulk;
            double cqdrn;
            double vt0;
            double[] args = new double[8];

            double m;

            EffectiveLength = Parameters.Length - ModelParameters.DeltaL * 1.0e-6; /* m */
            DrainArea = Parameters.DrainArea;
            SourceArea = Parameters.SourceArea;
            DrainPerimeter = Parameters.DrainPerimeter;
            SourcePerimeter = Parameters.SourcePerimeter;
            if ((DrainSatCurrent = DrainArea * ModelParameters.JctSatCurDensity) < 1e-15)
                DrainSatCurrent = 1.0e-15;
            if ((SourceSatCurrent = SourceArea * ModelParameters.JctSatCurDensity) < 1.0e-15)
                SourceSatCurrent = 1.0e-15;

            GateSourceOverlapCap = ModelParameters.GateSourceOverlapCap * Parameters.Width;
            GateDrainOverlapCap = ModelParameters.GateDrainOverlapCap * Parameters.Width;
            GateBulkOverlapCap = ModelParameters.GateBulkOverlapCap * EffectiveLength;
            von = ModelParameters.Type * this.Von;
            vdsat = ModelParameters.Type * this.Vdsat;
            vt0 = ModelParameters.Type * Param.B2vt0;

            Check = true;

            if (_iteration.Mode == IterationModes.Float || (_time != null && !_time.UseDc && _method != null && _method.BaseTime.Equals(0.0)) ||
                _iteration.Mode == IterationModes.Fix && !Parameters.Off)
            {
                vbs = ModelParameters.Type * (_bulk.Value - _sourcePrime.Value);
                vgs = ModelParameters.Type * (_gate.Value - _sourcePrime.Value);
                vds = ModelParameters.Type * (_drainPrime.Value - _sourcePrime.Value);

                vbd = vbs - vds;
                vgd = vgs - vds;
                vgdo = Vgs - Vds;

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
            vgb = vgs - vbs;


            if (vbs <= 0.0)
            {
                gbs = SourceSatCurrent / Constants.Vt0 + _iteration.Gmin;
                cbs = gbs * vbs;
            }
            else
            {
                evbs = Math.Exp(vbs / Constants.Vt0);
                gbs = SourceSatCurrent * evbs / Constants.Vt0 + _iteration.Gmin;
                cbs = SourceSatCurrent * (evbs - 1) + _iteration.Gmin * vbs;
            }
            if (vbd <= 0.0)
            {
                gbd = DrainSatCurrent / Constants.Vt0 + _iteration.Gmin;
                cbd = gbd * vbd;
            }
            else
            {
                evbd = Math.Exp(vbd / Constants.Vt0);
                gbd = DrainSatCurrent * evbd / Constants.Vt0 + _iteration.Gmin;
                cbd = DrainSatCurrent * (evbd - 1) + _iteration.Gmin * vbd;
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
            /* call B2evaluate to calculate drain current and its 
             * derivatives and charge and capacitances related to gate
             * drain, and bulk
             */
            if (vds >= 0)
            {
                B2Evaluate(vds, vbs, vgs, out gm, out gds, out gmbs, out qgate,
                    out qbulk, out qdrn, out cggb, out cgdb, out cgsb, out cbgb, out cbdb, out cbsb, out cdgb,
                    out cddb, out cdsb, out cdrain, out von, out vdsat);
            }
            else
            {
                B2Evaluate(-vds, vbd, vgd, out gm, out gds, out gmbs, out qgate,
                    out qbulk, out qsrc, out cggb, out cgsb, out cgdb, out cbgb, out cbsb, out cbdb, out csgb,
                    out cssb, out csdb, out cdrain, out von, out vdsat);
            }

            this.Von = ModelParameters.Type * von;
            this.Vdsat = ModelParameters.Type * vdsat;



            /*
             *  COMPUTE EQUIVALENT DRAIN CURRENT SOURCE
             */
            cd = this.Mode * cdrain - cbd;
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

                czbd = ModelParameters.UnitAreaJctCap * DrainArea;
                czbs = ModelParameters.UnitAreaJctCap * SourceArea;
                czbdsw = ModelParameters.UnitLengthSidewallJctCap * DrainPerimeter;
                czbssw = ModelParameters.UnitLengthSidewallJctCap * SourcePerimeter;
                PhiB = ModelParameters.BulkJctPotential;
                PhiBSW = ModelParameters.SidewallJctPotential;
                MJ = ModelParameters.BulkJctBotGradingCoeff;
                MJSW = ModelParameters.BulkJctSideGradingCoeff;

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

            if (_method == null && (_time != null && !_time.UseIc))
                goto line850;

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

                B2mosCap(vgd, vgs, vgb,
                    args,
                    /*
		            GateDrainOverlapCap,
		            GateSourceOverlapCap,GateBulkOverlapCap,
		            capbd,capbs,cggb,cgdb,cgsb,
		            */
                    cbgb, cbdb, cbsb, cdgb, cddb, cdsb,
                    out gcggb, out gcgdb, out gcgsb, out gcbgb, out gcbdb, out gcbsb, out gcdgb,
                    out gcddb, out gcdsb, out gcsgb, out gcsdb, out gcssb, ref qgate, ref qbulk,
                    ref qdrn, ref qsrc);
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

                B2mosCap(vgs, vgd, vgb, args,
                    /*
		            GateSourceOverlapCap,
                            GateDrainOverlapCap,GateBulkOverlapCap,
		            capbs,capbd,cggb,cgsb,cgdb,
		            */
                    cbgb, cbsb, cbdb, csgb, cssb, csdb,
                    out gcggb, out gcgsb, out gcgdb, out gcbgb, out gcbsb, out gcbdb, out gcsgb,
                    out gcssb, out gcsdb, out gcdgb, out gcdsb, out gcddb, ref qgate, ref qbulk,
                    ref qsrc, ref qdrn);
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
            /* initialize to zero charge conductance and current */
            ceqqg = ceqqb = ceqqd = 0.0;
            gcdgb = gcddb = gcdsb = 0.0;
            gcsgb = gcsdb = gcssb = 0.0;
            gcggb = gcgdb = gcgsb = 0.0;
            gcbgb = gcbdb = gcbsb = 0.0;
            goto line900;

        line860:
            /* evaluate equivalent charge current */
            cqgate = _qg.Derivative;
            cqbulk = _qb.Derivative;
            cqdrn = _qd.Derivative;
            ceqqg = cqgate - gcggb * vgb + gcgdb * vbd + gcgsb * vbs;
            ceqqb = cqbulk - gcbgb * vgb + gcbdb * vbd + gcbsb * vbs;
            ceqqd = cqdrn - gcdgb * vgb + gcddb * vbd + gcdsb * vbs;

        /*
         *  load current vector
         */
        line900:

            m = Parameters.Multiplier;

            ceqbs = ModelParameters.Type * (cbs - (gbs - _iteration.Gmin) * vbs);
            ceqbd = ModelParameters.Type * (cbd - (gbd - _iteration.Gmin) * vbd);

            ceqqg = ModelParameters.Type * ceqqg;
            ceqqb = ModelParameters.Type * ceqqb;
            ceqqd = ModelParameters.Type * ceqqd;
            if (this.Mode >= 0)
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

            _gPtr.Value -= m * (ceqqg);
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

        protected void B2Evaluate(double Vds, double Vbs, double Vgs,
           out double gm, out double gds, out double gmb, out double qg, out double qb, out double qd,
           out double cgg, out double cgd, out double cgs, out double cbg, out double cbd, out double cbs,
           out double cdg, out double cdd, out double cds, out double Ids, out double von,
           out double vdsat)
        {
            double Vth, Vdsat = 0.0;
            double Phisb, T1s, Eta, Gg, Aa, Inv_Aa, U1, U1s, Vc, Kk, SqrtKk;
            double dPhisb_dVb, dT1s_dVb, dVth_dVb, dVth_dVd, dAa_dVb, dVc_dVd;
            double dVc_dVg, dVc_dVb, dKk_dVc;
            double dVdsat_dVd = 0.0, dVdsat_dVg = 0.0, dVdsat_dVb = 0.0;
            double dUvert_dVg, dUvert_dVd, dUvert_dVb, Inv_Kk;
            double dUtot_dVd, dUtot_dVb, dUtot_dVg, Ai, Bi, Vghigh, Vglow, Vgeff, Vof;
            double Vbseff, Vgst, Vgdt, Qbulk, Utot;
            double T0, T1, T2, T3, T4, T5, Arg1, Arg2, Exp0 = 0.0, Exp1 = 0.0;
            double tmp, tmp1, tmp2, tmp3, Uvert, Beta1, Beta2, Beta0, dGg_dVb;
            double T6, T7, T8, T9, n = 0.0, ExpArg, ExpArg1;
            double Beta, dQbulk_dVb, dVgdt_dVg, dVgdt_dVd;
            double dVbseff_dVb, Ua, Ub, dVgdt_dVb, dQbulk_dVd;
            double Con1, Con3, Con4, SqrVghigh, SqrVglow, CubVghigh, CubVglow;
            double delta, Coeffa, Coeffb, Coeffc, Coeffd, Inv_Uvert, Inv_Utot;
            double Inv_Vdsat, tanh, Sqrsech, dBeta1_dVb, dU1_dVd, dU1_dVg, dU1_dVb;
            double Betaeff, FR, dFR_dVd, dFR_dVg, dFR_dVb, Betas, Beta3, Beta4;
            double dBeta_dVd, dBeta_dVg, dBeta_dVb, dVgeff_dVg, dVgeff_dVd, dVgeff_dVb;
            double dCon3_dVd, dCon3_dVb, dCon4_dVd, dCon4_dVb, dCoeffa_dVd, dCoeffa_dVb;
            double dCoeffb_dVd, dCoeffb_dVb, dCoeffc_dVd, dCoeffc_dVb;
            double dCoeffd_dVd, dCoeffd_dVb;
            bool ChargeComputationNeeded;
            int valuetypeflag;          /* added  3/19/90 JSD   */

            if (_method != null || ComputeSmallSignal)
                ChargeComputationNeeded = true;
            else
                ChargeComputationNeeded = false;

            if (Vbs < ModelTemperature.Vbb2) Vbs = ModelTemperature.Vbb2;
            if (Vgs > ModelTemperature.Vgg2) Vgs = ModelTemperature.Vgg2;
            if (Vds > ModelTemperature.Vdd2) Vds = ModelTemperature.Vdd2;

            /* Threshold Voltage. */
            if (Vbs <= 0.0)
            {
                Phisb = Param.B2phi - Vbs;
                dPhisb_dVb = -1.0;
                T1s = Math.Sqrt(Phisb);
                dT1s_dVb = -0.5 / T1s;
            }
            else
            {
                tmp = Param.B2phi / (Param.B2phi + Vbs);
                Phisb = Param.B2phi * tmp;
                dPhisb_dVb = -tmp * tmp;
                T1s = Param.Phis3 / (Param.B2phi + 0.5 * Vbs);
                dT1s_dVb = -0.5 * T1s * T1s / Param.Phis3;
            }

            Eta = Param.B2eta0 + Param.B2etaB * Vbs;
            Ua = Param.B2ua0 + Param.B2uaB * Vbs;
            Ub = Param.B2ub0 + Param.B2ubB * Vbs;
            U1s = Param.B2u10 + Param.B2u1B * Vbs;

            Vth = Param.B2vfb + Param.B2phi + Param.B2k1
                * T1s - Param.B2k2 * Phisb - Eta * Vds;
            dVth_dVd = -Eta;
            dVth_dVb = Param.B2k1 * dT1s_dVb + Param.B2k2
                 - Param.B2etaB * Vds;

            Vgst = Vgs - Vth;

            tmp = 1.0 / (1.744 + 0.8364 * Phisb);
            Gg = 1.0 - tmp;
            dGg_dVb = 0.8364 * tmp * tmp * dPhisb_dVb;
            T0 = Gg / T1s;
            tmp1 = 0.5 * T0 * Param.B2k1;
            Aa = 1.0 + tmp1;
            dAa_dVb = (Aa - 1.0) * (dGg_dVb / Gg - dT1s_dVb / T1s);
            Inv_Aa = 1.0 / Aa;

            Vghigh = Param.B2vghigh;
            Vglow = Param.B2vglow;

            if ((Vgst >= Vghigh) || (Param.B2n0 == 0.0))
            {
                Vgeff = Vgst;
                dVgeff_dVg = 1.0;
                dVgeff_dVd = -dVth_dVd;
                dVgeff_dVb = -dVth_dVb;
            }
            else
            {
                Vof = Param.B2vof0 + Param.B2vofB * Vbs
                + Param.B2vofD * Vds;
                n = Param.B2n0 + Param.B2nB / T1s
                  + Param.B2nD * Vds;
                tmp = 0.5 / (n * ModelTemperature.Vtm);

                ExpArg1 = -Vds / ModelTemperature.Vtm;
                ExpArg1 = Math.Max(ExpArg1, -30.0);
                Exp1 = Math.Exp(ExpArg1);
                tmp1 = 1.0 - Exp1;
                tmp1 = Math.Max(tmp1, 1.0e-18);
                tmp2 = 2.0 * Aa * tmp1;

                if (Vgst <= Vglow)
                {
                    ExpArg = Vgst * tmp;
                    ExpArg = Math.Max(ExpArg, -30.0);
                    Exp0 = Math.Exp(0.5 * Vof + ExpArg);
                    Vgeff = Math.Sqrt(tmp2) * ModelTemperature.Vtm * Exp0;
                    T0 = n * ModelTemperature.Vtm;
                    dVgeff_dVg = Vgeff * tmp;
                    dVgeff_dVd = dVgeff_dVg * (n / tmp1 * Exp1 - dVth_dVd - Vgst
                           * Param.B2nD / n + T0 * Param.B2vofD);
                    dVgeff_dVb = dVgeff_dVg * (Param.B2vofB * T0
                           - dVth_dVb + Param.B2nB * Vgst
                           / (n * T1s * T1s) * dT1s_dVb + T0 * Inv_Aa * dAa_dVb);
                }
                else
                {
                    ExpArg = Vglow * tmp;
                    ExpArg = Math.Max(ExpArg, -30.0);
                    Exp0 = Math.Exp(0.5 * Vof + ExpArg);
                    Vgeff = Math.Sqrt(2.0 * Aa * (1.0 - Exp1)) * ModelTemperature.Vtm * Exp0;
                    Con1 = Vghigh;
                    Con3 = Vgeff;
                    Con4 = Con3 * tmp;
                    SqrVghigh = Vghigh * Vghigh;
                    SqrVglow = Vglow * Vglow;
                    CubVghigh = Vghigh * SqrVghigh;
                    CubVglow = Vglow * SqrVglow;
                    T0 = 2.0 * Vghigh;
                    T1 = 2.0 * Vglow;
                    T2 = 3.0 * SqrVghigh;
                    T3 = 3.0 * SqrVglow;
                    T4 = Vghigh - Vglow;
                    T5 = SqrVghigh - SqrVglow;
                    T6 = CubVghigh - CubVglow;
                    T7 = Con1 - Con3;
                    delta = (T1 - T0) * T6 + (T2 - T3) * T5 + (T0 * T3 - T1 * T2) * T4;
                    delta = 1.0 / delta;
                    Coeffb = (T1 - Con4 * T0) * T6 + (Con4 * T2 - T3) * T5
                   + (T0 * T3 - T1 * T2) * T7;
                    Coeffc = (Con4 - 1.0) * T6 + (T2 - T3) * T7 + (T3 - Con4 * T2) * T4;
                    Coeffd = (T1 - T0) * T7 + (1.0 - Con4) * T5 + (Con4 * T0 - T1) * T4;
                    Coeffa = SqrVghigh * (Coeffc + Coeffd * T0);
                    Vgeff = (Coeffa + Vgst * (Coeffb + Vgst * (Coeffc + Vgst * Coeffd)))
                  * delta;
                    dVgeff_dVg = (Coeffb + Vgst * (2.0 * Coeffc + 3.0 * Vgst * Coeffd))
                           * delta;
                    T7 = Con3 * tmp;
                    T8 = dT1s_dVb * Param.B2nB / (T1s * T1s * n);
                    T9 = n * ModelTemperature.Vtm;
                    dCon3_dVd = T7 * (n * Exp1 / tmp1 - Vglow * Param.B2nD
                          / n + T9 * Param.B2vofD);
                    dCon3_dVb = T7 * (T9 * Inv_Aa * dAa_dVb + Vglow * T8
                          + T9 * Param.B2vofB);
                    dCon4_dVd = tmp * dCon3_dVd - T7 * Param.B2nD / n;
                    dCon4_dVb = tmp * dCon3_dVb + T7 * T8;

                    dCoeffb_dVd = dCon4_dVd * (T2 * T5 - T0 * T6) + dCon3_dVd
                    * (T1 * T2 - T0 * T3);
                    dCoeffc_dVd = dCon4_dVd * (T6 - T2 * T4) + dCon3_dVd * (T3 - T2);
                    dCoeffd_dVd = dCon4_dVd * (T0 * T4 - T5) + dCon3_dVd * (T0 - T1);
                    dCoeffa_dVd = SqrVghigh * (dCoeffc_dVd + dCoeffd_dVd * T0);

                    dVgeff_dVd = -dVgeff_dVg * dVth_dVd + (dCoeffa_dVd + Vgst
                           * (dCoeffb_dVd + Vgst * (dCoeffc_dVd + Vgst
                           * dCoeffd_dVd))) * delta;

                    dCoeffb_dVb = dCon4_dVb * (T2 * T5 - T0 * T6) + dCon3_dVb
                    * (T1 * T2 - T0 * T3);
                    dCoeffc_dVb = dCon4_dVb * (T6 - T2 * T4) + dCon3_dVb * (T3 - T2);
                    dCoeffd_dVb = dCon4_dVb * (T0 * T4 - T5) + dCon3_dVb * (T0 - T1);
                    dCoeffa_dVb = SqrVghigh * (dCoeffc_dVb + dCoeffd_dVb * T0);

                    dVgeff_dVb = -dVgeff_dVg * dVth_dVb + (dCoeffa_dVb + Vgst
                           * (dCoeffb_dVb + Vgst * (dCoeffc_dVb + Vgst
                           * dCoeffd_dVb))) * delta;
                }
            }

            if (Vgeff > 0.0)
            {
                Uvert = 1.0 + Vgeff * (Ua + Vgeff * Ub);
                Uvert = Math.Max(Uvert, 0.2);
                Inv_Uvert = 1.0 / Uvert;
                T8 = Ua + 2.0 * Ub * Vgeff;
                dUvert_dVg = T8 * dVgeff_dVg;
                dUvert_dVd = T8 * dVgeff_dVd;
                dUvert_dVb = T8 * dVgeff_dVb + Vgeff * (Param.B2uaB
                       + Vgeff * Param.B2ubB);

                T8 = U1s * Inv_Aa * Inv_Uvert;
                Vc = T8 * Vgeff;
                T9 = Vc * Inv_Uvert;
                dVc_dVg = T8 * dVgeff_dVg - T9 * dUvert_dVg;
                dVc_dVd = T8 * dVgeff_dVd - T9 * dUvert_dVd;
                dVc_dVb = T8 * dVgeff_dVb + Param.B2u1B * Vgeff * Inv_Aa
                    * Inv_Uvert - Vc * Inv_Aa * dAa_dVb - T9 * dUvert_dVb;


                tmp2 = Math.Sqrt(1.0 + 2.0 * Vc);
                Kk = 0.5 * (1.0 + Vc + tmp2);
                Inv_Kk = 1.0 / Kk;
                dKk_dVc = 0.5 + 0.5 / tmp2;
                SqrtKk = Math.Sqrt(Kk);

                T8 = Inv_Aa / SqrtKk;
                Vdsat = Vgeff * T8;
                Vdsat = Math.Max(Vdsat, 1.0e-18);
                Inv_Vdsat = 1.0 / Vdsat;
                T9 = 0.5 * Vdsat * Inv_Kk * dKk_dVc;
                dVdsat_dVd = T8 * dVgeff_dVd - T9 * dVc_dVd;
                dVdsat_dVg = T8 * dVgeff_dVg - T9 * dVc_dVg;
                dVdsat_dVb = T8 * dVgeff_dVb - T9 * dVc_dVb - Vdsat * Inv_Aa * dAa_dVb;

                Beta0 = Param.B2beta0 + Param.B2beta0B * Vbs;
                Betas = Param.B2betas0 + Param.B2betasB * Vbs;
                Beta2 = Param.B2beta20 + Param.B2beta2B * Vbs
                  + Param.B2beta2G * Vgs;
                Beta3 = Param.B2beta30 + Param.B2beta3B * Vbs
                  + Param.B2beta3G * Vgs;
                Beta4 = Param.B2beta40 + Param.B2beta4B * Vbs
                  + Param.B2beta4G * Vgs;
                Beta1 = Betas - (Beta0 + ModelParameters.Vdd * (Beta3 - ModelParameters.Vdd
                  * Beta4));

                T0 = Vds * Beta2 * Inv_Vdsat;
                T0 = Math.Min(T0, 30.0);
                T1 = Math.Exp(T0);
                T2 = T1 * T1;
                T3 = T2 + 1.0;
                tanh = (T2 - 1.0) / T3;
                Sqrsech = 4.0 * T2 / (T3 * T3);

                Beta = Beta0 + Beta1 * tanh + Vds * (Beta3 - Beta4 * Vds);
                T4 = Beta1 * Sqrsech * Inv_Vdsat;
                T5 = ModelParameters.Vdd * tanh;
                dBeta_dVd = Beta3 - 2.0 * Beta4 * Vds + T4 * (Beta2 - T0 * dVdsat_dVd);
                dBeta_dVg = T4 * (Param.B2beta2G * Vds - T0 * dVdsat_dVg)
                      + Param.B2beta3G * (Vds - T5)
                  - Param.B2beta4G * (Vds * Vds - ModelParameters.Vdd * T5);
                dBeta1_dVb = Param.Arg;
                dBeta_dVb = Param.B2beta0B + dBeta1_dVb * tanh + Vds
                      * (Param.B2beta3B - Vds * Param.B2beta4B)
                      + T4 * (Param.B2beta2B * Vds - T0 * dVdsat_dVb);


                if (Vgst > Vglow)
                {
                    if (Vds <= Vdsat) /* triode region */
                    {
                        T3 = Vds * Inv_Vdsat;
                        T4 = T3 - 1.0;
                        T2 = 1.0 - Param.B2u1D * T4 * T4;
                        U1 = U1s * T2;
                        Utot = Uvert + U1 * Vds;
                        Utot = Math.Max(Utot, 0.5);
                        Inv_Utot = 1.0 / Utot;
                        T5 = 2.0 * U1s * Param.B2u1D * Inv_Vdsat * T4;
                        dU1_dVd = T5 * (T3 * dVdsat_dVd - 1.0);
                        dU1_dVg = T5 * T3 * dVdsat_dVg;
                        dU1_dVb = T5 * T3 * dVdsat_dVb + Param.B2u1B * T2;
                        dUtot_dVd = dUvert_dVd + U1 + Vds * dU1_dVd;
                        dUtot_dVg = dUvert_dVg + Vds * dU1_dVg;
                        dUtot_dVb = dUvert_dVb + Vds * dU1_dVb;

                        tmp1 = (Vgeff - 0.5 * Aa * Vds);
                        tmp3 = tmp1 * Vds;
                        Betaeff = Beta * Inv_Utot;
                        Ids = Betaeff * tmp3;
                        T6 = Ids / Betaeff * Inv_Utot;
                        gds = T6 * (dBeta_dVd - Betaeff * dUtot_dVd) + Betaeff * (tmp1
                            + (dVgeff_dVd - 0.5 * Aa) * Vds);
                        gm = T6 * (dBeta_dVg - Betaeff * dUtot_dVg) + Betaeff * Vds
                       * dVgeff_dVg;
                        gmb = T6 * (dBeta_dVb - Betaeff * dUtot_dVb) + Betaeff * Vds
                            * (dVgeff_dVb - 0.5 * Vds * dAa_dVb);
                    }
                    else  /* Saturation */
                    {
                        tmp1 = Vgeff * Inv_Aa * Inv_Kk;
                        tmp3 = 0.5 * Vgeff * tmp1;
                        Betaeff = Beta * Inv_Uvert;
                        Ids = Betaeff * tmp3;
                        T0 = Ids / Betaeff * Inv_Uvert;
                        T1 = Betaeff * Vgeff * Inv_Aa * Inv_Kk;
                        T2 = Ids * Inv_Kk * dKk_dVc;

                        if (Param.B2ai0 != 0.0)
                        {
                            Ai = Param.B2ai0 + Param.B2aiB * Vbs;
                            Bi = Param.B2bi0 + Param.B2biB * Vbs;
                            T5 = Bi / (Vds - Vdsat);
                            T5 = Math.Min(T5, 30.0);
                            T6 = Math.Exp(-T5);
                            FR = 1.0 + Ai * T6;
                            T7 = T5 / (Vds - Vdsat);
                            T8 = (1.0 - FR) * T7;
                            dFR_dVd = T8 * (dVdsat_dVd - 1.0);
                            dFR_dVg = T8 * dVdsat_dVg;
                            dFR_dVb = T8 * dVdsat_dVb + T6 * (Param.B2aiB - Ai
                                * Param.B2biB / (Vds - Vdsat));

                            gds = (T0 * (dBeta_dVd - Betaeff * dUvert_dVd) + T1
                             * dVgeff_dVd - T2 * dVc_dVd) * FR + Ids * dFR_dVd;
                            gm = (T0 * (dBeta_dVg - Betaeff * dUvert_dVg)
                            + T1 * dVgeff_dVg - T2 * dVc_dVg) * FR + Ids * dFR_dVg;
                            gmb = (T0 * (dBeta_dVb - Betaeff * dUvert_dVb) + T1
                             * dVgeff_dVb - T2 * dVc_dVb - Ids * Inv_Aa * dAa_dVb)
                             * FR + Ids * dFR_dVb;
                            Ids *= FR;
                        }
                        else
                        {
                            gds = T0 * (dBeta_dVd - Betaeff * dUvert_dVd) + T1
                            * dVgeff_dVd - T2 * dVc_dVd;
                            gm = T0 * (dBeta_dVg - Betaeff * dUvert_dVg) + T1 * dVgeff_dVg
                            - T2 * dVc_dVg;
                            gmb = T0 * (dBeta_dVb - Betaeff * dUvert_dVb) + T1
                             * dVgeff_dVb - T2 * dVc_dVb - Ids * Inv_Aa * dAa_dVb;
                        }
                    } /* end of Saturation */
                }
                else
                {
                    T0 = Exp0 * Exp0;
                    T1 = Exp1;
                    Ids = Beta * ModelTemperature.Vtm * ModelTemperature.Vtm * T0 * (1.0 - T1);
                    T2 = Ids / Beta;
                    T4 = n * ModelTemperature.Vtm;
                    T3 = Ids / T4;
                    if ((Vds > Vdsat) && Param.B2ai0 != 0.0)
                    {
                        Ai = Param.B2ai0 + Param.B2aiB * Vbs;
                        Bi = Param.B2bi0 + Param.B2biB * Vbs;
                        T5 = Bi / (Vds - Vdsat);
                        T5 = Math.Min(T5, 30.0);
                        T6 = Math.Exp(-T5);
                        FR = 1.0 + Ai * T6;
                        T7 = T5 / (Vds - Vdsat);
                        T8 = (1.0 - FR) * T7;
                        dFR_dVd = T8 * (dVdsat_dVd - 1.0);
                        dFR_dVg = T8 * dVdsat_dVg;
                        dFR_dVb = T8 * dVdsat_dVb + T6 * (Param.B2aiB - Ai
                            * Param.B2biB / (Vds - Vdsat));
                    }
                    else
                    {
                        FR = 1.0;
                        dFR_dVd = 0.0;
                        dFR_dVg = 0.0;
                        dFR_dVb = 0.0;
                    }

                    gds = (T2 * dBeta_dVd + T3 * (Param.B2vofD * T4 - dVth_dVd
                     - Param.B2nD * Vgst / n) + Beta * ModelTemperature.Vtm
                 * T0 * T1) * FR + Ids * dFR_dVd;
                    gm = (T2 * dBeta_dVg + T3) * FR + Ids * dFR_dVg;
                    gmb = (T2 * dBeta_dVb + T3 * (Param.B2vofB * T4 - dVth_dVb
                     + Param.B2nB * Vgst / (n * T1s * T1s) * dT1s_dVb)) * FR
                     + Ids * dFR_dVb;
                    Ids *= FR;
                }
            }
            else
            {
                Ids = 0.0;
                gm = 0.0;
                gds = 0.0;
                gmb = 0.0;
            }

            /* Some Limiting of DC Parameters */
            gds = Math.Max(gds, 1.0e-20);


            if ((ModelParameters.ChannelChargePartitionFlag > 1) || ((!ChargeComputationNeeded) && (ModelParameters.ChannelChargePartitionFlag > -5)))
            {
                qg = 0.0;
                qd = 0.0;
                qb = 0.0;
                cgg = 0.0;
                cgs = 0.0;
                cgd = 0.0;
                cdg = 0.0;
                cds = 0.0;
                cdd = 0.0;
                cbg = 0.0;
                cbs = 0.0;
                cbd = 0.0;
                goto finished;
            }
            else
            {
                if (Vbs < 0.0)
                {
                    Vbseff = Vbs;
                    dVbseff_dVb = 1.0;
                }
                else
                {
                    Vbseff = Param.B2phi - Phisb;
                    dVbseff_dVb = -dPhisb_dVb;
                }
                Arg1 = Vgs - Vbseff - Param.B2vfb;
                Arg2 = Arg1 - Vgst;
                Qbulk = Param.One_Third_CoxWL * Arg2;
                dQbulk_dVb = Param.One_Third_CoxWL * (dVth_dVb - dVbseff_dVb);
                dQbulk_dVd = Param.One_Third_CoxWL * dVth_dVd;
                if (Arg1 <= 0.0)
                {
                    qg = Param.CoxWL * Arg1;
                    qb = -(qg);
                    qd = 0.0;

                    cgg = Param.CoxWL;
                    cgd = 0.0;
                    cgs = -cgg * (1.0 - dVbseff_dVb);

                    cdg = 0.0;
                    cdd = 0.0;
                    cds = 0.0;

                    cbg = -Param.CoxWL;
                    cbd = 0.0;
                    cbs = -cgs;
                }
                else if (Vgst <= 0.0)
                {
                    T2 = Arg1 / Arg2;
                    T3 = T2 * T2 * (Param.CoxWL - Param.Two_Third_CoxWL
                   * T2);

                    qg = Param.CoxWL * Arg1 * (1.0 - T2 * (1.0 - T2 / 3.0));
                    qb = -(qg);
                    qd = 0.0;

                    cgg = Param.CoxWL * (1.0 - T2 * (2.0 - T2));
                    tmp = T3 * dVth_dVb - (cgg + T3) * dVbseff_dVb;
                    cgd = T3 * dVth_dVd;
                    cgs = -(cgg + cgd + tmp);

                    cdg = 0.0;
                    cdd = 0.0;
                    cds = 0.0;

                    cbg = -cgg;
                    cbd = -cgd;
                    cbs = -cgs;
                }
                else
                {
                    if (Vgst < Param.B2vghigh)
                    {
                        Uvert = 1.0 + Vgst * (Ua + Vgst * Ub);
                        Uvert = Math.Max(Uvert, 0.2);
                        Inv_Uvert = 1.0 / Uvert;
                        dUvert_dVg = Ua + 2.0 * Ub * Vgst;
                        dUvert_dVd = -dUvert_dVg * dVth_dVd;
                        dUvert_dVb = -dUvert_dVg * dVth_dVb + Vgst
                     * (Param.B2uaB + Vgst * Param.B2ubB);

                        T8 = U1s * Inv_Aa * Inv_Uvert;
                        Vc = T8 * Vgst;
                        T9 = Vc * Inv_Uvert;
                        dVc_dVg = T8 - T9 * dUvert_dVg;
                        dVc_dVd = -T8 * dVth_dVd - T9 * dUvert_dVd;
                        dVc_dVb = -T8 * dVth_dVb + Param.B2u1B * Vgst * Inv_Aa
                            * Inv_Uvert - Vc * Inv_Aa * dAa_dVb - T9 * dUvert_dVb;

                        tmp2 = Math.Sqrt(1.0 + 2.0 * Vc);
                        Kk = 0.5 * (1.0 + Vc + tmp2);
                        Inv_Kk = 1.0 / Kk;
                        dKk_dVc = 0.5 + 0.5 / tmp2;
                        SqrtKk = Math.Sqrt(Kk);

                        T8 = Inv_Aa / SqrtKk;
                        Vdsat = Vgst * T8;
                        T9 = 0.5 * Vdsat * Inv_Kk * dKk_dVc;
                        dVdsat_dVd = -T8 * dVth_dVd - T9 * dVc_dVd;
                        dVdsat_dVg = T8 - T9 * dVc_dVg;
                        dVdsat_dVb = -T8 * dVth_dVb - T9 * dVc_dVb
                       - Vdsat * Inv_Aa * dAa_dVb;
                    }
                    if (Vds >= Vdsat)
                    {       /* saturation region */
                        cgg = Param.Two_Third_CoxWL;
                        cgd = -cgg * dVth_dVd + dQbulk_dVd;
                        tmp = -cgg * dVth_dVb + dQbulk_dVb;
                        cgs = -(cgg + cgd + tmp);

                        cbg = 0.0;
                        cbd = -dQbulk_dVd;
                        cbs = dQbulk_dVd + dQbulk_dVb;

                        cdg = -0.4 * cgg;
                        tmp = -cdg * dVth_dVb;
                        cdd = -cdg * dVth_dVd;
                        cds = -(cdg + cdd + tmp);

                        qb = -Qbulk;
                        qg = Param.Two_Third_CoxWL * Vgst + Qbulk;
                        qd = cdg * Vgst;
                    }
                    else
                    {       /* linear region  */
                        T7 = Vds / Vdsat;
                        T8 = Vgst / Vdsat;
                        T6 = T7 * T8;
                        T9 = 1.0 - T7;
                        Vgdt = Vgst * T9;
                        T0 = Vgst / (Vgst + Vgdt);
                        T1 = Vgdt / (Vgst + Vgdt);
                        T5 = T0 * T1;
                        T2 = 1.0 - T1 + T5;
                        T3 = 1.0 - T0 + T5;

                        dVgdt_dVg = T9 + T6 * dVdsat_dVg;
                        dVgdt_dVd = T6 * dVdsat_dVd - T8 - T9 * dVth_dVd;
                        dVgdt_dVb = T6 * dVdsat_dVb - T9 * dVth_dVb;

                        qg = Param.Two_Third_CoxWL * (Vgst + Vgdt
                           - Vgdt * T0) + Qbulk;
                        qb = -Qbulk;
                        qd = -Param.One_Third_CoxWL * (0.2 * Vgdt
                       + 0.8 * Vgst + Vgdt * T1
                       + 0.2 * T5 * (Vgdt - Vgst));

                        cgg = Param.Two_Third_CoxWL * (T2 + T3 * dVgdt_dVg);
                        tmp = dQbulk_dVb + Param.Two_Third_CoxWL * (T3 * dVgdt_dVb
                            - T2 * dVth_dVb);
                        cgd = Param.Two_Third_CoxWL * (T3 * dVgdt_dVd
                         - T2 * dVth_dVd) + dQbulk_dVd;
                        cgs = -(cgg + cgd + tmp);

                        T2 = 0.8 - 0.4 * T1 * (2.0 * T1 + T0 + T0 * (T1 - T0));
                        T3 = 0.2 + T1 + T0 * (1.0 - 0.4 * T0 * (T1 + 3.0 * T0));
                        cdg = -Param.One_Third_CoxWL * (T2 + T3 * dVgdt_dVg);
                        tmp = Param.One_Third_CoxWL * (T2 * dVth_dVb
                        - T3 * dVgdt_dVb);
                        cdd = Param.One_Third_CoxWL * (T2 * dVth_dVd
                         - T3 * dVgdt_dVd);
                        cds = -(cdg + tmp + cdd);

                        cbg = 0.0;
                        cbd = -dQbulk_dVd;
                        cbs = dQbulk_dVd + dQbulk_dVb;
                    }
                }
            }

        finished:       /* returning Values to Calling Routine */
            valuetypeflag = (int)ModelParameters.ChannelChargePartitionFlag;
            switch (valuetypeflag)
            {
                case 0:
                    Ids = Math.Max(Ids, 1e-50);
                    break;
                case -1:
                    Ids = Math.Max(Ids, 1e-50);
                    break;
                case -2:
                    Ids = gm;
                    break;
                case -3:
                    Ids = gds;
                    break;
                case -4:
                    Ids = 1.0 / gds;
                    break;
                case -5:
                    Ids = gmb;
                    break;
                case -6:
                    Ids = qg / 1.0e-12;
                    break;
                case -7:
                    Ids = qb / 1.0e-12;
                    break;
                case -8:
                    Ids = qd / 1.0e-12;
                    break;
                case -9:
                    Ids = -(qb + qg + qd) / 1.0e-12;
                    break;
                case -10:
                    Ids = cgg / 1.0e-12;
                    break;
                case -11:
                    Ids = cgd / 1.0e-12;
                    break;
                case -12:
                    Ids = cgs / 1.0e-12;
                    break;
                case -13:
                    Ids = -(cgg + cgd + cgs) / 1.0e-12;
                    break;
                case -14:
                    Ids = cbg / 1.0e-12;
                    break;
                case -15:
                    Ids = cbd / 1.0e-12;
                    break;
                case -16:
                    Ids = cbs / 1.0e-12;
                    break;
                case -17:
                    Ids = -(cbg + cbd + cbs) / 1.0e-12;
                    break;
                case -18:
                    Ids = cdg / 1.0e-12;
                    break;
                case -19:
                    Ids = cdd / 1.0e-12;
                    break;
                case -20:
                    Ids = cds / 1.0e-12;
                    break;
                case -21:
                    Ids = -(cdg + cdd + cds) / 1.0e-12;
                    break;
                case -22:
                    Ids = -(cgg + cdg + cbg) / 1.0e-12;
                    break;
                case -23:
                    Ids = -(cgd + cdd + cbd) / 1.0e-12;
                    break;
                case -24:
                    Ids = -(cgs + cds + cbs) / 1.0e-12;
                    break;
                default:
                    Ids = Math.Max(Ids, 1.0e-50);
                    break;
            }
            von = Vth;
            vdsat = Vdsat;
        }

        protected void B2mosCap(double vgd,
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
            qBulkPointer -= qgb;
            qDrainPointer -= qgd;
            qSourcePointer = -(qGatePointer + qBulkPointer + qDrainPointer);
        }
    }
}