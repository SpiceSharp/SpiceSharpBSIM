using SpiceSharp.Attributes;
using SpiceSharp.Components.BSIM3v24Behaviors;

namespace SpiceSharp.Components
{
    /// <summary>
    /// BSIM 3.2.4
    /// </summary>
    public class BSIM3v24 : Component
    {
        static BSIM3v24()
        {
            RegisterBehaviorFactory(typeof(BSIM3v24), new Behaviors.BehaviorFactoryDictionary()
            {
                { typeof(TemperatureBehavior), e => new TemperatureBehavior(e.Name) },
                { typeof(BiasingBehavior), e => new BiasingBehavior(e.Name) },
                { typeof(TransientBehavior), e => new TransientBehavior(e.Name) },
                { typeof(FrequencyBehavior), e => new FrequencyBehavior(e.Name) },
                { typeof(NoiseBehavior), e => new NoiseBehavior(e.Name) }
            });
        }

        /// <summary>
        /// Number of pins
        /// </summary>
        [ParameterName("pincount"), ParameterInfo("Number of pins")]
        public const int BSIM3PinCount = 4;

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="name">Name</param>
        public BSIM3v24(string name)
            : base(name, BSIM3PinCount)
        {
            ParameterSets.Add(new BaseParameters());
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="name">Name</param>
        /// <param name="drain">Drain</param>
        /// <param name="gate">Gate</param>
        /// <param name="source">Source</param>
        /// <param name="bulk">Bulk</param>
        public BSIM3v24(string name, string drain, string gate, string source, string bulk)
            : this(name)
        {
            ParameterSets.Add(new BaseParameters());
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
        public BSIM3v24(string name, string drain, string gate, string source, string bulk, double width, double length)
            : base(name, BSIM3PinCount)
        {
            var bp = new BaseParameters();
            ParameterSets.Add(bp);
            bp.Width.Value = width;
            bp.Length.Value = length;
            Connect(drain, gate, source, bulk);
        }
    }
}
