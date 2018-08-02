using SpiceSharp.Attributes;
using SpiceSharp.Behaviors;
using SpiceSharp.Components.BSIM2Behaviors;

namespace SpiceSharp.Components
{
    /// <summary>
    /// BSIM2 model transistor
    /// </summary>
    public class BSIM2 : Component
    {
        /// <summary>
        /// Set the model
        /// </summary>
        /// <param name="model">Model</param>
        public void SetModel(BSIM2Model model) => Model = model;

        /// <summary>
        /// Number of pins
        /// </summary>
        [ParameterName("pincount"), ParameterInfo("Number of pins")]
        public const int BSIM2PinCount = 4;

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="name">Name</param>
        public BSIM2(Identifier name)
            : base(name, BSIM2PinCount)
        {
            // Add parameters
            ParameterSets.Add(new BaseParameters());

            // Add behaviors
            Behaviors.Add(typeof(BaseTemperatureBehavior), () => new TemperatureBehavior(Name));
            Behaviors.Add(typeof(BaseLoadBehavior), () => new LoadBehavior(Name));
            Behaviors.Add(typeof(TransientBehavior), () => new TransientBehavior(Name));
            Behaviors.Add(typeof(BaseFrequencyBehavior), () => new FrequencyBehavior(Name));
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="name">Name</param>
        /// <param name="drain">Drain</param>
        /// <param name="gate">Gate</param>
        /// <param name="source">Source</param>
        /// <param name="bulk">Bulk</param>
        public BSIM2(Identifier name, Identifier drain, Identifier gate, Identifier source, Identifier bulk)
            : base(name, BSIM2PinCount)
        {
            // Add parameters
            ParameterSets.Add(new BaseParameters());

            // Add behaviors
            Behaviors.Add(typeof(BaseTemperatureBehavior), () => new TemperatureBehavior(Name));
            Behaviors.Add(typeof(BaseLoadBehavior), () => new LoadBehavior(Name));
            Behaviors.Add(typeof(TransientBehavior), () => new TransientBehavior(Name));
            Behaviors.Add(typeof(BaseFrequencyBehavior), () => new FrequencyBehavior(Name));

            Connect(drain, gate, source, bulk);
        }
    }
}
