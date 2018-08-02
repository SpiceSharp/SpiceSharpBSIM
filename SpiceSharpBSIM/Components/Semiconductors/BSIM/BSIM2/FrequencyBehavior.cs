using System;
using System.Numerics;
using SpiceSharp.Algebra;
using SpiceSharp.Behaviors;
using SpiceSharp.Simulations;

namespace SpiceSharp.Components.BSIM2Behaviors
{
    /// <summary>
    /// Frequency behavior for a <see cref="BSIM2"/>
    /// </summary>
    public class FrequencyBehavior : BaseFrequencyBehavior, IConnectedBehavior
    {
        /// <summary>
        /// Necessary behaviors and parameters
        /// </summary>
        private LoadBehavior _load;
        private TemperatureBehavior _temp;

        /// <summary>
        /// Nodes
        /// </summary>
        private int _drainNode, _gateNode, _sourceNode, _bulkNode, _sourceNodePrime, _drainNodePrime;
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

            // Get behaviors
            _temp = provider.GetBehavior<TemperatureBehavior>("entity");
            _load = provider.GetBehavior<LoadBehavior>("entity");
        }

        /// <summary>
        /// Connect
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
        /// Load frequency behavior
        /// </summary>
        /// <param name="simulation">Simulation</param>
        public override void Load(FrequencySimulation simulation)
        {
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
            double omega; // angular fequency of the signal
            int xnrm, xrev;
            var state = simulation.ComplexState;
            var pParam = _temp.Param;

            omega = state.Laplace.Imaginary;

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

            gdpr = _temp.DrainConductance;
            gspr = _temp.SourceConductance;
            gm = _load.Gm;
            gds = _load.Gds;
            gmbs = _load.Gmbs;
            gbd = _load.Gbd;
            gbs = _load.Gbs;
            capbd = _load.Capbd;
            capbs = _load.Capbs;

            // Charge oriented model parameters
            cggb = _load.Cggb;
            cgsb = _load.Cgsb;
            cgdb = _load.Cgdb;

            cbgb = _load.Cbgb;
            cbsb = _load.Cbsb;
            cbdb = _load.Cbdb;

            cdgb = _load.Cdgb;
            cdsb = _load.Cdsb;
            cddb = _load.Cddb;

            xcdgb = (cdgb - pParam.B2GDoverlapCap) * omega;
            xcddb = (cddb + capbd + pParam.B2GDoverlapCap) * omega;
            xcdsb = cdsb * omega;
            xcsgb = -(cggb + cbgb + cdgb + pParam.B2GSoverlapCap) * omega;
            xcsdb = -(cgdb + cbdb + cddb) * omega;
            xcssb = (capbs + pParam.B2GSoverlapCap - (cgsb + cbsb + cdsb)) * omega;
            xcggb = (cggb + pParam.B2GDoverlapCap + pParam.B2GSoverlapCap + pParam.B2GBoverlapCap) * omega;
            xcgdb = (cgdb - pParam.B2GDoverlapCap) * omega;
            xcgsb = (cgsb - pParam.B2GSoverlapCap) * omega;
            xcbgb = (cbgb - pParam.B2GBoverlapCap) * omega;
            xcbdb = (cbdb - capbd) * omega;
            xcbsb = (cbsb - capbs) * omega;

            GgPtr.Value += new Complex(0.0, xcggb);
            BbPtr.Value += new Complex(0.0, -xcbgb - xcbdb - xcbsb);
            DPdpPtr.Value += new Complex(0.0, xcddb);
            SPspPtr.Value += new Complex(0.0, xcssb);
            GbPtr.Value += new Complex(0.0, -xcggb - xcgdb - xcgsb);
            GdpPtr.Value += new Complex(0.0, xcgdb);
            GspPtr.Value += new Complex(0.0, xcgsb);
            BgPtr.Value += new Complex(0.0, xcbgb);
            BdpPtr.Value += new Complex(0.0, xcbdb);
            BspPtr.Value += new Complex(0.0, xcbsb);
            DPgPtr.Value += new Complex(0.0, xcdgb);
            DPbPtr.Value += new Complex(0.0, -xcdgb - xcddb - xcdsb);
            DPspPtr.Value += new Complex(0.0, xcdsb);
            SPgPtr.Value += new Complex(0.0, xcsgb);
            SPbPtr.Value += new Complex(0.0, -xcsgb - xcsdb - xcssb);
            SPdpPtr.Value += new Complex(0.0, xcsdb);
            DdPtr.Value += gdpr;
            SsPtr.Value += gspr;
            BbPtr.Value += gbd + gbs;
            DPdpPtr.Value += gdpr + gds + gbd + xrev * (gm + gmbs);
            SPspPtr.Value += gspr + gds + gbs + xnrm * (gm + gmbs);
            DdpPtr.Value -= gdpr;
            SspPtr.Value -= gspr;
            BdpPtr.Value -= gbd;
            BspPtr.Value -= gbs;
            DPdPtr.Value -= gdpr;
            DPgPtr.Value += (xnrm - xrev) * gm;
            DPbPtr.Value += -gbd + (xnrm - xrev) * gmbs;
            DPspPtr.Value += -gds - xnrm * (gm + gmbs);
            SPgPtr.Value += -(xnrm - xrev) * gm;
            SPsPtr.Value -= gspr;
            SPbPtr.Value += -gbs - (xnrm - xrev) * gmbs;
            SPdpPtr.Value += -gds - xrev * (gm + gmbs);
        }
    }
}
