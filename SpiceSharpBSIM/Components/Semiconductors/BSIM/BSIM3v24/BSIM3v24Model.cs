using SpiceSharp.Components.BSIM3v24Behaviors;

namespace SpiceSharp.Components
{
    /// <summary>
    /// Model for a <see cref="BSIM3v24"/>
    /// </summary>
    public class BSIM3v24Model : Model
    {
        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="name"></param>
        public BSIM3v24Model(Identifier name) : base(name)
        {
            // Add parameter sets
            ParameterSets.Add(new ModelBaseParameters());

            // Add behaviors
            Behaviors.Add(typeof(ModelTemperatureBehavior), () => new ModelTemperatureBehavior(Name));
        }
    }
}
