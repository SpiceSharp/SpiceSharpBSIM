using SpiceSharp.Attributes;
using SpiceSharp.Components.BSIM3Behaviors;

namespace SpiceSharp.Components
{
    /// <summary>
    /// BSIM3 model transistor
    /// </summary>
    public class BSIM3 : Component
    {
        /// <summary>
        /// Set the model
        /// </summary>
        /// <param name="model">Model</param>
        public void SetModel(BSIM3Model model) => Model = model;

        /// <summary>
        /// Number of pins
        /// </summary>
        [ParameterName("pincount"), ParameterInfo("Number of pins")]
        public const int BSIM3PinCount = 4;

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="name">Name</param>
        public BSIM3(string name)
            : base(name, BSIM3PinCount)
        {
            // Add parameters
            ParameterSets.Add(new BaseParameters());

            // Add behaviors
            Behaviors.Add(typeof(TemperatureBehavior), () => new TemperatureBehavior(Name));
            Behaviors.Add(typeof(LoadBehavior), () => new LoadBehavior(Name));
            Behaviors.Add(typeof(TransientBehavior), () => new TransientBehavior(Name));
            Behaviors.Add(typeof(FrequencyBehavior), () => new FrequencyBehavior(Name));
            Behaviors.Add(typeof(NoiseBehavior), () => new NoiseBehavior(Name));
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="name">Name</param>
        /// <param name="drain">Drain</param>
        /// <param name="gate">Gate</param>
        /// <param name="source">Source</param>
        /// <param name="bulk">Bulk</param>
        public BSIM3(string name, string drain, string gate, string source, string bulk)
            : base(name, BSIM3PinCount)
        {
            // Add parameters
            ParameterSets.Add(new BaseParameters());

            // Add behaviors
            Behaviors.Add(typeof(TemperatureBehavior), () => new TemperatureBehavior(Name));
            Behaviors.Add(typeof(LoadBehavior), () => new LoadBehavior(Name));
            Behaviors.Add(typeof(TransientBehavior), () => new TransientBehavior(Name));
            Behaviors.Add(typeof(FrequencyBehavior), () => new FrequencyBehavior(Name));
            Behaviors.Add(typeof(NoiseBehavior), () => new NoiseBehavior(Name));

            Connect(drain, gate, source, bulk);
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="name">Name</param>
        /// <param name="drain">Drain</param>
        /// <param name="gate">Gate</param>
        /// <param name="source">Source</param>
        /// <param name="bulk">Bulk</param>
        /// <param name="width">Transistor width</param>
        /// <param name="length">Transistor length</param>
        public BSIM3(string name, string drain, string gate, string source, string bulk, double width, double length)
            : base(name, BSIM3PinCount)
        {
            // Add parameters
            var bp = new BaseParameters();
            ParameterSets.Add(bp);
            bp.Width.Value = width;
            bp.Length.Value = length;

            // Add behaviors
            Behaviors.Add(typeof(TemperatureBehavior), () => new TemperatureBehavior(Name));
            Behaviors.Add(typeof(LoadBehavior), () => new LoadBehavior(Name));
            Behaviors.Add(typeof(TransientBehavior), () => new TransientBehavior(Name));
            Behaviors.Add(typeof(FrequencyBehavior), () => new FrequencyBehavior(Name));

            Connect(drain, gate, source, bulk);
        }
    }
}
