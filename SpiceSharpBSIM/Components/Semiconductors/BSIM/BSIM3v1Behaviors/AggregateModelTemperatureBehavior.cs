using SpiceSharp;
using SpiceSharp.Attributes;
using SpiceSharp.Behaviors;
using SpiceSharp.Components;
using System.Collections.Generic;

namespace SpiceSharpBSIM.Components.Semiconductors.BSIM.BSIM3v1Behaviors
{
    /// <summary>
    /// A temperature behavior for a <see cref="BSIM3v1AggregateModel"/>.
    /// </summary>
    [BehaviorFor(typeof(BSIM3v1AggregateModel)), AddBehaviorIfNo(typeof(ITemperatureBehavior))]
    public class AggregateModelTemperatureBehavior : Behavior, ITemperatureBehavior
    {
        private readonly List<ModelTemperatureBehavior> _modelTemperatureBehaviors;

        /// <summary>
        /// Creates a new <see cref="AggregateModelTemperatureBehavior"/>.
        /// </summary>
        /// <param name="context">The binding context.</param>
        public AggregateModelTemperatureBehavior(AggregateModelBindingContext context)
            : base(context)
        {
            _modelTemperatureBehaviors = [];
            foreach (var behaviors in context.Models)
            {
                // We require the model temperature behavior to exist.
                var behavior = behaviors.GetValue<ModelTemperatureBehavior>();
                _modelTemperatureBehaviors.Add(behavior);
            }
        }

        /// <summary>
        /// Gets a model for a given size.
        /// </summary>
        /// <param name="width">The width.</param>
        /// <param name="length">The length.</param>
        /// <returns>Returns the model behavior.</returns>
        /// <exception cref="SpiceSharpException">Thrown if no valid model could be found.</exception>
        public ModelTemperatureBehavior GetModel(double width, double length)
        {
            foreach (var behavior in _modelTemperatureBehaviors)
            {
                if (behavior.Parameters.Lmin.Given && length < behavior.Parameters.Lmin)
                    continue;
                if (behavior.Parameters.Lmax.Given && length > behavior.Parameters.Lmax)
                    continue;
                if (behavior.Parameters.Wmin.Given && width < behavior.Parameters.Wmin)
                    continue;
                if (behavior.Parameters.Wmax.Given && width > behavior.Parameters.Wmax)
                    continue;
                return behavior;
            }
            throw new SpiceSharpException($"There is no model that can be used for sizes W={width} and L={length}.");
        }

        /// <inheritdoc />
        public void Temperature()
        {
            // The models have been created already, so no need to do anything here.
        }
    }
}
