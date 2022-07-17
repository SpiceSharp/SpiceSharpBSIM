using SpiceSharp.Simulations;

namespace SpiceSharpBSIM.Components
{
    /// <summary>
    /// A simple state that can be used for DC analysis to avoid needing an
    /// integration method.
    /// </summary>
    public class SimpleDerivative : IDerivative
    {
        /// <inheritdoc />
        public double Value { get; set; }

        /// <inheritdoc />
        public double Derivative => 0.0;

        /// <inheritdoc />
        public void Accept()
        {
        }

        /// <inheritdoc />
        public double GetPreviousValue(int index) => Value;

        /// <inheritdoc />
        public double GetPreviousDerivative(int index) => 0;

        /// <inheritdoc />
        public JacobianInfo GetContributions(double coefficient, double currentValue) => default;

        /// <inheritdoc />
        public JacobianInfo GetContributions(double coefficient) => default;

        /// <inheritdoc />
        public void Derive()
        {
        }
    }
}
