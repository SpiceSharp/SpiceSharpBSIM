using System;
using System.Numerics;
using SpiceSharp;
using SpiceSharp.Algebra;
using SpiceSharp.Attributes;
using SpiceSharp.Behaviors;
using SpiceSharp.Components;
using SpiceSharp.Simulations;

namespace SpiceSharpBSIM.Components.Semiconductors.BSIM.BSIM2Behaviors
{
    /// <summary>
    /// Frequency behavior for a <see cref="BSIM2"/>
    /// </summary>
    [BehaviorFor(typeof(BSIM2)), AddBehaviorIfNo(typeof(IFrequencyBehavior))]
    public class FrequencyBehavior : BiasingBehavior, IFrequencyBehavior
    {
        private readonly IComplexSimulationState _state;
        private readonly IVariable<Complex> _drain, _gate, _source, _bulk, _drainPrime, _sourcePrime;
        private readonly Element<Complex> _ddPtr, _ggPtr, _ssPtr, _bbPtr,
            _dpdpPtr, _spspPtr, _ddpPtr, _gbPtr, _gdpPtr, _gspPtr, _sspPtr, _bdpPtr, _bspPtr,
            _spsPtr, _dpspPtr, _dpdPtr, _bgPtr, _dpgPtr, _spgPtr, _dpbPtr, _spbPtr, _spdpPtr;


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

            if (!ModelParameters.SheetResistance.Equals(0.0) && !Parameters.DrainSquares.Equals(0.0))
                _drainPrime = _state.CreatePrivateVariable(context.Behaviors.Name.Combine("drain"), Units.Volt);
            else
                _drainPrime = _drain;
            if (!ModelParameters.SheetResistance.Equals(0.0) && !Parameters.SourceSquares.Equals(0.0))
                _sourcePrime = _state.CreatePrivateVariable(context.Behaviors.Name.Combine("source"), Units.Volt);
            else
                _sourcePrime = _source;

            var drain = _state.Map[_drain];
            var gate = _state.Map[_gate];
            var source = _state.Map[_source];
            var bulk = _state.Map[_bulk];
            var drainPrime = _state.Map[_drainPrime];
            var sourcePrime = _state.Map[_sourcePrime];
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
        }

        /// <summary>
        /// Initializes the parameters.
        /// </summary>
        void IFrequencyBehavior.InitializeParameters()
        {
            ComputeSmallSignal = true;
            ((IBiasingBehavior)this).Load();
        }

        /// <summary>
        /// Load frequency behavior
        /// </summary>
        /// <param name="simulation">Simulation</param>
        void IFrequencyBehavior.Load()
        {
            int xnrm;
            int xrev;
            double gdpr;
            double gspr;
            double gm;
            double gds;
            double gmbs;
            double gbd;
            double gbs;
            double capbd;
            double capbs;
            double xcggb;
            double xcgdb;
            double xcgsb;
            double xcbgb;
            double xcbdb;
            double xcbsb;
            double xcddb;
            double xcssb;
            double xcdgb;
            double xcsgb;
            double xcdsb;
            double xcsdb;
            double cggb;
            double cgdb;
            double cgsb;
            double cbgb;
            double cbdb;
            double cbsb;
            double cddb;
            double cdgb;
            double cdsb;
            double omega; /* angular fequency of the signal */

            double m;     /* parallel multiplier */


            omega = _state.Laplace.Imaginary;
            if (Mode >= 0)
            {
                xnrm = 1;
                xrev = 0;
            }
            else
            {
                xnrm = 0;
                xrev = 1;
            }
            gdpr = DrainConductance;
            gspr = SourceConductance;
            gm = Gm;
            gds = Gds;
            gmbs = Gmbs;
            gbd = Gbd;
            gbs = Gbs;
            capbd = Capbd;
            capbs = Capbs;

            /*
             *    charge oriented model parameters
             */

            cggb = Cggb;
            cgsb = Cgsb;
            cgdb = Cgdb;

            cbgb = Cbgb;
            cbsb = Cbsb;
            cbdb = Cbdb;

            cdgb = Cdgb;
            cdsb = Cdsb;
            cddb = Cddb;

            xcdgb = (cdgb - Param.B2GDoverlapCap) * omega;
            xcddb = (cddb + capbd + Param.B2GDoverlapCap) * omega;
            xcdsb = cdsb * omega;
            xcsgb = -(cggb + cbgb + cdgb + Param.B2GSoverlapCap)
          * omega;
            xcsdb = -(cgdb + cbdb + cddb) * omega;
            xcssb = (capbs + Param.B2GSoverlapCap
          - (cgsb + cbsb + cdsb)) * omega;
            xcggb = (cggb + Param.B2GDoverlapCap
          + Param.B2GSoverlapCap
          + Param.B2GBoverlapCap) * omega;
            xcgdb = (cgdb - Param.B2GDoverlapCap) * omega;
            xcgsb = (cgsb - Param.B2GSoverlapCap) * omega;
            xcbgb = (cbgb - Param.B2GBoverlapCap) * omega;
            xcbdb = (cbdb - capbd) * omega;
            xcbsb = (cbsb - capbs) * omega;

            m = Parameters.Multiplier;

            _ggPtr.Value += new Complex(0, m * (xcggb));
            _bbPtr.Value += new Complex(m * (gbd + gbs), m * (-xcbgb - xcbdb - xcbsb));
            _dpdpPtr.Value += new Complex(m * (gdpr + gds + gbd + xrev * (gm + gmbs)), m * (xcddb));
            _spspPtr.Value += new Complex(m * (gspr + gds + gbs + xnrm * (gm + gmbs)), m * (xcssb));
            _gbPtr.Value += new Complex(0, m * (-xcggb - xcgdb - xcgsb));
            _gdpPtr.Value += new Complex(0, m * (xcgdb));
            _gspPtr.Value += new Complex(0, m * (xcgsb));
            _bgPtr.Value += new Complex(0, m * (xcbgb));
            _bdpPtr.Value += new Complex(-m * (gbd), m * (xcbdb));
            _bspPtr.Value += new Complex(-m * (gbs), m * (xcbsb));
            _dpgPtr.Value += new Complex(m * ((xnrm - xrev) * gm), m * (xcdgb));
            _dpbPtr.Value += new Complex(m * (-gbd + (xnrm - xrev) * gmbs), m * (-xcdgb - xcddb - xcdsb));
            _dpspPtr.Value += new Complex(m * (-gds - xnrm * (gm + gmbs)), m * (xcdsb));
            _spgPtr.Value += new Complex(m * (-(xnrm - xrev) * gm), m * (xcsgb));
            _spbPtr.Value += new Complex(m * (-gbs - (xnrm - xrev) * gmbs), m * (-xcsgb - xcsdb - xcssb));
            _spdpPtr.Value += new Complex(m * (-gds - xrev * (gm + gmbs)), m * (xcsdb));
            _ddPtr.Value += m * (gdpr);
            _ssPtr.Value += m * (gspr);
            _ddpPtr.Value -= m * (gdpr);
            _sspPtr.Value -= m * (gspr);
            _dpdPtr.Value -= m * (gdpr);
            _spsPtr.Value -= m * (gspr);
        }
    }
}
