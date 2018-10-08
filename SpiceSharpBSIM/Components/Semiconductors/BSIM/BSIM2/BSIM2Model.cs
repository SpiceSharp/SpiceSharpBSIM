using SpiceSharp.Components.BSIM2Behaviors;

namespace SpiceSharp.Components
{
    /// <summary>
    /// Model for the <see cref="BSIM2"/> component
    /// </summary>
    public class BSIM2Model : Model
    {
        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="name">Name</param>
        public BSIM2Model(string name) : base(name)
        {
            // Add parameters
            ParameterSets.Add(new ModelBaseParameters());

            // Add behaviors
            Behaviors.Add(typeof(ModelTemperatureBehavior), () => new ModelTemperatureBehavior(Name));
        }
    }
}