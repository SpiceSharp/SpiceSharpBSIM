using SpiceSharp.Attributes;
using SpiceSharp.Behaviors;
using SpiceSharp.Components;
using SpiceSharp.Simulations;
using System.Numerics;
using SpiceSharp.Algebra;
using SpiceSharp;

namespace SpiceSharpBSIM.Components.Semiconductors.BSIM.BSIM3v1Behaviors
{
    /// <summary>
    /// Frequency behavior for a <see cref="BSIM3v1"/>.
    /// </summary>
    [BehaviorFor(typeof(BSIM3v1)), AddBehaviorIfNo(typeof(IFrequencyBehavior))]
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
            if (!_drainConductance.Equals(0.0))
                _drainPrime = _state.CreatePrivateVariable(Name.Combine("drain"), Units.Volt);
            else
                _drainPrime = _drain;
            if (!_sourceConductance.Equals(0.0))
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

			double dxpart, sxpart, cqgb, cqdb, cqsb, cqbb, xcqgb, xcqdb, xcqsb, xcqbb;

			double m;
			omega = _state.Laplace.Imaginary;

			if (this._mode >= 0)
			{
				Gm = this._gm;
				Gmbs = this._gmbs;
				FwdSum = Gm + Gmbs;
				RevSum = 0.0;
				cggb = this._cggb;
				cgsb = this._cgsb;
				cgdb = this._cgdb;

				cbgb = this._cbgb;
				cbsb = this._cbsb;
				cbdb = this._cbdb;

				cdgb = this._cdgb;
				cdsb = this._cdsb;
				cddb = this._cddb;

				cqgb = this._cqgb;
				cqdb = this._cqdb;
				cqsb = this._cqsb;
				cqbb = this._cqbb;
				sxpart = 0.6;
				dxpart = 0.4;

			}
			else
			{
				Gm = -this._gm;
				Gmbs = -this._gmbs;
				FwdSum = 0.0;
				RevSum = -Gm - Gmbs;
				cggb = this._cggb;
				cgsb = this._cgdb;
				cgdb = this._cgsb;

				cbgb = this._cbgb;
				cbsb = this._cbdb;
				cbdb = this._cbsb;

				cdgb = -(this._cdgb + cggb + cbgb);
				cdsb = -(this._cddb + cgsb + cbsb);
				cddb = -(this._cdsb + cgdb + cbdb);

				cqgb = this._cqgb;
				cqdb = this._cqsb;
				cqsb = this._cqdb;
				cqbb = this._cqbb;
				sxpart = 0.4;
				dxpart = 0.6;
			}

			gdpr = this._drainConductance;
			gspr = this._sourceConductance;
			gds = this._gds;
			gbd = this._gbd;
			gbs = this._gbs;
			capbd = this._capbd;
			capbs = this._capbs;

			GSoverlapCap = this._cgso;
			GDoverlapCap = this._cgdo;
			GBoverlapCap = Param.BSIM3v1cgbo;

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
			xcqgb = cqgb * omega;
			xcqdb = cqdb * omega;
			xcqsb = cqsb * omega;
			xcqbb = cqbb * omega;

			m = Parameters.M;

			this._ggPtr.Value += new Complex(0, m * xcggb);
			this._bbPtr.Value -= new Complex(0, m * (xcbgb + xcbdb + xcbsb));
			this._dpdpPtr.Value += new Complex(0, m * xcddb);
			this._spspPtr.Value += new Complex(0, m * xcssb);
			this._gbPtr.Value -= new Complex(0, m * (xcggb + xcgdb + xcgsb));
			this._gdpPtr.Value += new Complex(0, m * xcgdb);
			this._gspPtr.Value += new Complex(0, m * xcgsb);
			this._bgPtr.Value += new Complex(0, m * xcbgb);
			this._bdpPtr.Value += new Complex(0, m * xcbdb);
			this._bspPtr.Value += new Complex(0, m * xcbsb);
			this._dpgPtr.Value += new Complex(0, m * xcdgb);
			this._dpbPtr.Value -= new Complex(0, m * (xcdgb + xcddb + xcdsb));
			this._dpspPtr.Value += new Complex(0, m * xcdsb);
			this._spgPtr.Value += new Complex(0, m * xcsgb);
			this._spbPtr.Value -= new Complex(0, m * (xcsgb + xcsdb + xcssb));
			this._spdpPtr.Value += new Complex(0, m * xcsdb);

			this._qqPtr.Value += new Complex(0, m * omega);

			this._qgPtr.Value -= new Complex(0, m * xcqgb);
			this._qdpPtr.Value -= new Complex(0, m * xcqdb);
			this._qspPtr.Value -= new Complex(0, m * xcqsb);
			this._qbPtr.Value -= new Complex(0, m * xcqbb);


			_ddPtr.Value += m * gdpr;
			_ssPtr.Value += m * gspr;
			_bbPtr.Value += m * (gbd + gbs);
			_dpdpPtr.Value += m * (gdpr + gds + gbd + RevSum + dxpart * this._gtd);
			_spspPtr.Value += m * (gspr + gds + gbs + FwdSum + sxpart * this._gts);
			_ddpPtr.Value -= m * gdpr;
			_sspPtr.Value -= m * gspr;
			_bdpPtr.Value -= m * gbd;
			_bspPtr.Value -= m * gbs;
			_dpdPtr.Value -= m * gdpr;
			_dpgPtr.Value += m * (Gm + dxpart * this._gtg);
			_dpbPtr.Value -= m * (gbd - Gmbs - dxpart * this._gtb);
			_dpspPtr.Value -= m * (gds + FwdSum - dxpart * this._gts);
			_spgPtr.Value -= m * (Gm - sxpart * this._gtg);
			_spsPtr.Value -= m * gspr;
			_spbPtr.Value -= m * (gbs + Gmbs - sxpart * this._gtg);
			_spdpPtr.Value -= m * (gds + RevSum - sxpart * this._gtd);
			_ggPtr.Value -= m * this._gtg;
			_gbPtr.Value -= m * this._gtb;
			_gdpPtr.Value -= m * this._gtd;
			_gspPtr.Value -= m * this._gts;

			_qqPtr.Value += m * this._gtau;

			_dpqPtr.Value += m * dxpart * this._gtau;
			_spqPtr.Value += m * sxpart * this._gtau;
			_gqPtr.Value -= m * this._gtau;

			_qgPtr.Value += m * this._gtg;
			_qdpPtr.Value += m * this._gtd;
			_qspPtr.Value += m * this._gts;
			_qbPtr.Value += m * this._gtb;
		}
    }
}
