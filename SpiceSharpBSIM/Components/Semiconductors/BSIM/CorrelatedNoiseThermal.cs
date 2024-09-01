using SpiceSharp;
using SpiceSharp.Simulations;
using System;
using System.Numerics;

namespace SpiceSharpBSIM.Components.Semiconductors.BSIM
{
    /// <summary>
    /// A class that describes correlated thermal noise between multiple ports.
    /// </summary>
    public class CorrelatedNoiseThermal : NoiseSource
    {
        private IVariable<Complex> _n1, _n2, _n3, _n4;

        /// <summary>
        /// Creates a new <see cref="CorrelatedNoiseThermal"/>.
        /// </summary>
        /// <param name="name">The name of the noise source.</param>
        /// <param name="n1">Node 1.</param>
        /// <param name="n2">Node 2.</param>
        /// <param name="n3">Node 3.</param>
        /// <param name="n4">Node 4.</param>
        /// <exception cref="ArgumentNullException">Thrown if any node is <c>null</c>.</exception>
        public CorrelatedNoiseThermal(string name, IVariable<Complex> n1, IVariable<Complex> n2,
            IVariable<Complex> n3, IVariable<Complex> n4)
            : base(name)
        {
            _n1 = n1 ?? throw new ArgumentNullException(nameof(n1));
            _n2 = n2 ?? throw new ArgumentNullException(nameof(n2));
            _n3 = n3 ?? throw new ArgumentNullException(nameof(n3));
            _n4 = n4 ?? throw new ArgumentNullException(nameof(n4));
        }

        /// <summary>
        /// Connects the correlated noise.
        /// </summary>
        /// <param name="n1">Node 1.</param>
        /// <param name="n2">Node 2.</param>
        /// <param name="n3">Node 3.</param>
        /// <param name="n4">Node 4.</param>
        /// <exception cref="ArgumentNullException">Thrown if any node is <c>null</c>.</exception>
        public void Connect(IVariable<Complex> n1, IVariable<Complex> n2,
            IVariable<Complex> n3, IVariable<Complex> n4)
        {
            _n1 = n1 ?? throw new ArgumentNullException(nameof(n1));
            _n2 = n2 ?? throw new ArgumentNullException(nameof(n2));
            _n3 = n3 ?? throw new ArgumentNullException(nameof(n3));
            _n4 = n4 ?? throw new ArgumentNullException(nameof(n4));
        }

        /// <summary>
        /// Computes the noise density.
        /// </summary>
        /// <param name="param1">Parameter 1.</param>
        /// <param name="param2">Parameter 2.</param>
        /// <param name="phi21">Phi 21.</param>
        /// <param name="temperature">The temperature.</param>
        public void Compute(double param1, double param2, double phi21, double temperature)
        {
            var val1 = _n1.Value - _n2.Value;
            var val2 = _n3.Value - _n4.Value;
            double T0 = Math.Sqrt(param1);
            double T1 = Math.Sqrt(param2);
            double T2 = T1 * Math.Cos(phi21);
            double T3 = T1 * Math.Sin(phi21);
            var @out = T0 * val1 + T2 * val2 + T3 * new Complex(-val2.Imaginary, val2.Real);
            double param_gain = @out.Real * @out.Real + @out.Imaginary * @out.Imaginary;

            OutputNoiseDensity = 4.0 * Constants.Boltzmann * temperature * param_gain;
        }
    }
}
