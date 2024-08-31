using System.Numerics;
using SpiceSharp;
using SpiceSharp.Algebra;
using SpiceSharp.Attributes;
using SpiceSharp.Behaviors;
using SpiceSharp.Components;
using SpiceSharp.Simulations;

namespace SpiceSharpBSIM.Components.Semiconductors.BSIM.BSIM1Behaviors
{
    /// <summary>
    /// Frequency behavior for a <see cref="BSIM1"/>
    /// </summary>
    [BehaviorFor(typeof(BSIM1)), AddBehaviorIfNo(typeof(IFrequencyBehavior))]
    public class FrequencyBehavior : BiasingBehavior, IFrequencyBehavior
    {
        private readonly IComplexSimulationState _state;
        protected readonly IVariable<Complex> _drain, _gate, _source, _bulk, _drainPrime, _sourcePrime;
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
            InitializeSmallSignal = true;
            ((IBiasingBehavior)this).Load();
            InitializeSmallSignal = false;
        }

        /// <summary>
        /// Load frequency behavior
        /// </summary>
        /// <param name="simulation">Simulation</param>
        void IFrequencyBehavior.Load()
        {
            int xnrm;
            int xrev;

            var omega = _state.Laplace.Imaginary;

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

            var gdpr = DrainConductance;
            var gspr = SourceConductance;
            var gm = Gm;
            var gds = Gds;
            var gmbs = Gmbs;
            var gbd = Gbd;
            var gbs = Gbs;
            var capbd = Capbd;
            var capbs = Capbs;

            /*
             *    charge oriented model parameters
             */
            var cggb = Cggb;
            var cgsb = Cgsb;
            var cgdb = Cgdb;

            var cbgb = Cbgb;
            var cbsb = Cbsb;
            var cbdb = Cbdb;

            var cdgb = Cdgb;
            var cdsb = Cdsb;
            var cddb = Cddb;

            var xcdgb = (cdgb - GDoverlapCap) * omega;
            var xcddb = (cddb + capbd + GDoverlapCap) * omega;
            var xcdsb = cdsb * omega;
            var xcsgb = -(cggb + cbgb + cdgb + GSoverlapCap) * omega;
            var xcsdb = -(cgdb + cbdb + cddb) * omega;
            var xcssb = (capbs + GSoverlapCap - (cgsb + cbsb + cdsb)) * omega;
            var xcggb = (cggb + GDoverlapCap + GSoverlapCap +
                            GBoverlapCap) * omega;
            var xcgdb = (cgdb - GDoverlapCap) * omega;
            var xcgsb = (cgsb - GSoverlapCap) * omega;
            var xcbgb = (cbgb - GBoverlapCap) * omega;
            var xcbdb = (cbdb - capbd) * omega;
            var xcbsb = (cbsb - capbs) * omega;

            var m = Parameters.Multiplier;
            _ggPtr.Value += new Complex(0.0, m * xcggb);
            _gbPtr.Value += new Complex(0.0, m * (-xcggb - xcgdb - xcgsb));
            _gdpPtr.Value += new Complex(0.0, m * xcgdb);
            _gspPtr.Value += new Complex(0.0, m * xcgsb);
            _bgPtr.Value += new Complex(0.0, m * xcbgb);
            _ddPtr.Value += m * gdpr;
            _ssPtr.Value += m * gspr;
            _bbPtr.Value += new Complex(m * (gbd + gbs), m * (-xcbgb - xcbdb - xcbsb));
            _dpdpPtr.Value += new Complex(m * (gdpr + gds + gbd + xrev * (gm + gmbs)), m * xcddb);
            _spspPtr.Value += new Complex(m * (gspr + gds + gbs + xnrm * (gm + gmbs)), m * xcssb);
            _ddpPtr.Value -= m * gdpr;
            _sspPtr.Value -= m * gspr;
            _bdpPtr.Value += new Complex(m * -gbd, m * xcbdb);
            _bspPtr.Value += new Complex(m * -gbs, m * xcbsb);
            _dpdPtr.Value -= m * gdpr;
            _dpgPtr.Value += new Complex(m * (xnrm - xrev) * gm, m * xcdgb);
            _dpbPtr.Value += new Complex(m * (-gbd + (xnrm - xrev) * gmbs), m * (-xcdgb - xcddb - xcdsb));
            _dpspPtr.Value += new Complex(m * (-gds - xnrm * (gm + gmbs)), m * xcdsb);
            _spgPtr.Value += new Complex(m * (-(xnrm - xrev) * gm), m * xcsgb);
            _spsPtr.Value -= m * gspr;
            _spbPtr.Value += new Complex(m * (-gbs - (xnrm - xrev) * gmbs), m * (-xcsgb - xcsdb - xcssb));
            _spdpPtr.Value += new Complex(m * (-gds - xrev * (gm + gmbs)), m * xcsdb);
        }
    }
}
