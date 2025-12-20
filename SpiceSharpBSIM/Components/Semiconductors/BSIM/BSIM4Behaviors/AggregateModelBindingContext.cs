using SpiceSharp.Attributes;
using SpiceSharp.Behaviors;
using SpiceSharp.Components;
using SpiceSharp.Entities;
using SpiceSharp.Simulations;
using System.Collections.Generic;

namespace SpiceSharpBSIM.Components.Semiconductors.BSIM.BSIM4Behaviors
{
    /// <summary>
    /// A binding context for an aggregate model.
    /// </summary>
    [BindingContextFor(typeof(BSIM4AggregateModel))]
    public class AggregateModelBindingContext : BindingContext
    {
        private readonly HashSet<IBehaviorContainer> _set;

        /// <summary>
        /// Gets the models
        /// </summary>
        public IEnumerable<IBehaviorContainer> Models => _set;

        /// <summary>
        /// Creates a binding context for aggregate models.
        /// </summary>
        /// <param name="aggregateModel">The model</param>
        /// <param name="simulation">The simulation.</param>
        /// <param name="behaviors">The behaviors.</param>
        public AggregateModelBindingContext(BSIM4AggregateModel aggregateModel,
            ISimulation simulation,
            IBehaviorContainer behaviors)
            : base(aggregateModel, simulation, behaviors)
        {
            // First make the behaviors for all sub-models
            _set = new HashSet<IBehaviorContainer>();
            foreach (var model in aggregateModel.Models)
            {
                // Create the behaviors for the different models
                model.CreateBehaviors(simulation);
                _set.Add(simulation.EntityBehaviors[model.Name]);
            }
        }
    }
}
