﻿using System;
using System.Numerics;
using SpiceSharp.Algebra;
using SpiceSharp.Behaviors;
using SpiceSharp.Simulations;

namespace SpiceSharp.Components.BSIM1Behaviors
{
    /// <summary>
    /// Frequency behavior for a <see cref="BSIM1"/>
    /// </summary>
    public class FrequencyBehavior : Behavior, IFrequencyBehavior
    {
        /// <summary>
        /// Necessary behaviors and parameters
        /// </summary>
        private ModelBaseParameters _mbp;
        private BaseParameters _bp;
        private TemperatureBehavior _temp;
        private BiasingBehavior _load;

        /// <summary>
        /// Nodes
        /// </summary>
        private int _drainNode, _gateNode, _sourceNode, _bulkNode, _drainNodePrime, _sourceNodePrime;

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

        private ComplexSimulationState _state;

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="name">Name</param>
        public FrequencyBehavior(string name) : base(name)
        {
        }

        /// <summary>
        /// Set up the behavior.
        /// </summary>
        /// <param name="simulation">The simulation.</param>
        /// <param name="context">The provider.</param>
        /// <exception cref="ArgumentNullException">provider</exception>
        public override void Bind(Simulation simulation, BindingContext context)
        {
            base.Bind(simulation, context);

            // Get parameters
            _mbp = context.GetParameterSet<ModelBaseParameters>("model");
            _bp = context.GetParameterSet<BaseParameters>();

            // Get behaviors
            _temp = context.GetBehavior<TemperatureBehavior>();
            _load = context.GetBehavior<BiasingBehavior>();

            if (context is ComponentBindingContext cc)
            {
                _drainNode = cc.Pins[0];
                _gateNode = cc.Pins[1];
                _sourceNode = cc.Pins[2];
                _bulkNode = cc.Pins[3];
            }

            _state = ((FrequencySimulation)simulation).ComplexState;
            var solver = _state.Solver;
            _drainNodePrime = _load.DrainNodePrime;
            _sourceNodePrime = _load.SourceNodePrime;
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
        }

        /// <summary>
        /// Initializes the parameters.
        /// </summary>
        void IFrequencyBehavior.InitializeParameters()
        {
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

            if (_load.Mode >= 0)
            {
                xnrm = 1;
                xrev = 0;
            }
            else
            {
                xnrm = 0;
                xrev = 1;
            }

            var gdpr = _temp.DrainConductance;
            var gspr = _temp.SourceConductance;
            var gm = _load.Gm;
            var gds = _load.Gds;
            var gmbs = _load.Gmbs;
            var gbd = _load.Gbd;
            var gbs = _load.Gbs;
            var capbd = _load.Capbd;
            var capbs = _load.Capbs;

            /*
             *    charge oriented model parameters
             */
            var cggb = _load.Cggb;
            var cgsb = _load.Cgsb;
            var cgdb = _load.Cgdb;

            var cbgb = _load.Cbgb;
            var cbsb = _load.Cbsb;
            var cbdb = _load.Cbdb;

            var cdgb = _load.Cdgb;
            var cdsb = _load.Cdsb;
            var cddb = _load.Cddb;

            var xcdgb = (cdgb - _temp.GDoverlapCap) * omega;
            var xcddb = (cddb + capbd + _temp.GDoverlapCap) * omega;
            var xcdsb = cdsb * omega;
            var xcsgb = -(cggb + cbgb + cdgb + _temp.GSoverlapCap) * omega;
            var xcsdb = -(cgdb + cbdb + cddb) * omega;
            var xcssb = (capbs + _temp.GSoverlapCap - (cgsb + cbsb + cdsb)) * omega;
            var xcggb = (cggb + _temp.GDoverlapCap + _temp.GSoverlapCap +
                            _temp.GBoverlapCap) * omega;
            var xcgdb = (cgdb - _temp.GDoverlapCap) * omega;
            var xcgsb = (cgsb - _temp.GSoverlapCap) * omega;
            var xcbgb = (cbgb - _temp.GBoverlapCap) * omega;
            var xcbdb = (cbdb - capbd) * omega;
            var xcbsb = (cbsb - capbs) * omega;

            var m = _bp.Multiplier;
            GgPtr.Value += new Complex(0.0, m * xcggb);
            GbPtr.Value += new Complex(0.0, m * (-xcggb - xcgdb - xcgsb));
            GdpPtr.Value += new Complex(0.0, m * xcgdb);
            GspPtr.Value += new Complex(0.0, m *xcgsb);
            BgPtr.Value += new Complex(0.0, m * xcbgb);
            DdPtr.Value += m * gdpr;
            SsPtr.Value += m * gspr;
            BbPtr.Value += new Complex(m * (gbd + gbs), m * (-xcbgb - xcbdb - xcbsb));
            DPdpPtr.Value += new Complex(m * (gdpr + gds + gbd + xrev * (gm + gmbs)), m * xcddb);
            SPspPtr.Value += new Complex(m * (gspr + gds + gbs + xnrm * (gm + gmbs)), m * xcssb);
            DdpPtr.Value -= m * gdpr;
            SspPtr.Value -= m * gspr;
            BdpPtr.Value += new Complex(m * -gbd, m * xcbdb);
            BspPtr.Value += new Complex(m * -gbs, m * xcbsb);
            DPdPtr.Value -= m * gdpr;
            DPgPtr.Value += new Complex(m * (xnrm - xrev) * gm, m * xcdgb);
            DPbPtr.Value += new Complex(m * (-gbd + (xnrm - xrev) * gmbs), m * (-xcdgb - xcddb - xcdsb));
            DPspPtr.Value += new Complex(m * (-gds - xnrm * (gm + gmbs)), m * xcdsb);
            SPgPtr.Value += new Complex(m * (-(xnrm - xrev) * gm), m * xcsgb);
            SPsPtr.Value -= m * gspr;
            SPbPtr.Value += new Complex(m * (-gbs - (xnrm - xrev) * gmbs), m * (-xcsgb - xcsdb - xcssb));
            SPdpPtr.Value += new Complex(m * (-gds - xrev * (gm + gmbs)), m * xcsdb);
        }
    }
}
