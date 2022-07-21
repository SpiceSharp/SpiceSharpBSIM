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
    public class ModelTemperatureBehavior : Behavior, ITemperatureBehavior, IParameterized<ModelParameters>
    {
        /// <inheritdoc />
        public ModelParameters Parameters { get; }

        public double Cox { get; private set; }

        /// <summary>
        /// Creates a new <see cref="ModelTemperatureBehavior"/>.
        /// </summary>
        /// <param name="context">The binding context.</param>
        public ModelTemperatureBehavior(IBindingContext context)
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
