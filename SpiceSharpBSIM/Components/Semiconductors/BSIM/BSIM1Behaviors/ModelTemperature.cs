using SpiceSharp.Attributes;
using SpiceSharp.Behaviors;
using SpiceSharp.Entities;
using SpiceSharp.ParameterSets;
using SpiceSharp.Components;

namespace SpiceSharpBSIM.Components.Semiconductors.BSIM.BSIM1Behaviors
{
    /// <summary>
    /// Temperature behavior for a <see cref="BSIM1Model"/>.
    /// </summary>
    [BehaviorFor(typeof(BSIM1Model)), AddBehaviorIfNo(typeof(ITemperatureBehavior))]
    [GeneratedParameters]
    public partial class ModelTemperature : Behavior, ITemperatureBehavior, IParameterized<ModelParameters>
    {
        /// <inheritdoc />
        public ModelParameters Parameters { get; }

        [ParameterName("cox")]
        public double Cox { get; private set; }

        /// <summary>
        /// Creates a new <see cref="ModelTemperature"/>.
        /// </summary>
        /// <param name="context">The binding context.</param>
        public ModelTemperature(IBindingContext context)
            : base(context)
        {
            Parameters = context.GetParameterSet<ModelParameters>();
        }

        /// <inheritdoc />
        public void Temperature()
        {
            Cox = 3.453e-13 / (Parameters.OxideThickness * 1.0e-4);
        }
    }
}
