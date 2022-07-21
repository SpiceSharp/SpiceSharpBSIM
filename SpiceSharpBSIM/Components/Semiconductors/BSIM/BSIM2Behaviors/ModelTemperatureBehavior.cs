using System;
using System.Collections.Generic;
using SpiceSharp;
using SpiceSharp.Attributes;
using SpiceSharp.Behaviors;
using SpiceSharp.Components;
using SpiceSharp.Entities;
using SpiceSharp.ParameterSets;

namespace SpiceSharpBSIM.Components.Semiconductors.BSIM.BSIM2Behaviors
{
    /// <summary>
    /// Temperature behavior for a <see cref="BSIM2Model" />
    /// </summary>
    [BehaviorFor(typeof(BSIM2Model)), AddBehaviorIfNo(typeof(ITemperatureBehavior))]
    public class ModelTemperatureBehavior : Behavior, ITemperatureBehavior, IParameterized<ModelParameters>
    {
        /// <inheritdoc />
        public ModelParameters Parameters { get; }

        public double Vtm { get; private set; }
        public double Cox { get; private set; }
        public double Vdd2 { get; private set; }
        public double Vgg2 { get; private set; }
        public double Vbb2 { get; private set; }

        /// <summary>
        /// Size-dependent parameters for this model
        /// </summary>
        public Dictionary<Tuple<double, double>, BSIM2SizeDependParams> SizedParams { get; } = new Dictionary<Tuple<double, double>, BSIM2SizeDependParams>();

        /// <summary>
        /// Constructor
        /// </summary>
        public ModelTemperatureBehavior(BindingContext context)
            : base(context)
        {
            Parameters = context.GetParameterSet<ModelParameters>();
        }

        /// <summary>
        /// Temperature behavior
        /// </summary>
        void ITemperatureBehavior.Temperature()
        {
            SizedParams.Clear();
            Cox = 3.453e-13 / (Parameters.Tox * 1.0e-4);/*in F/cm**2 */
            Vdd2 = 2.0 * Parameters.Vdd;
            Vgg2 = 2.0 * Parameters.Vgg;
            Vbb2 = 2.0 * Parameters.Vbb;
            Vtm = 8.625e-5 * (Parameters.Temp + 273.0);
        }
    }
}