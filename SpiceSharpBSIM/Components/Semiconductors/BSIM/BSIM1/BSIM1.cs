using SpiceSharp.Attributes;
using SpiceSharp.Components.BSIM1Behaviors;

namespace SpiceSharp.Components
{
    /// <summary>
    /// BSIM1 model transistor
    /// </summary>
    public class BSIM1 : Component
    {
        static BSIM1()
        {
            RegisterBehaviorFactory(typeof(BSIM1), new Behaviors.BehaviorFactoryDictionary()
            {
                { typeof(TemperatureBehavior), e => new TemperatureBehavior(e.Name) },
                { typeof(BiasingBehavior), e => new BiasingBehavior(e.Name) },
                { typeof(TransientBehavior), e => new TransientBehavior(e.Name) },
                { typeof(FrequencyBehavior), e => new FrequencyBehavior(e.Name) }
            });
        }

        /// <summary>
        /// Number of pins
        /// </summary>
        [ParameterName("pincount"), ParameterInfo("Number of pins")]
        public const int BSIM1PinCount = 4;

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="name">Name</param>
        public BSIM1(string name) 
            : base(name, BSIM1PinCount)
        {
            // Add parameters
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
        public BSIM1(string name, string drain, string gate, string source, string bulk)
            : base(name, BSIM1PinCount)
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
        public BSIM1(string name, string drain, string gate, string source, string bulk, double width, double length)
            : base(name, BSIM1PinCount)
        {
            var bp = new BaseParameters();
            ParameterSets.Add(bp);
            bp.Width.Value = width;
            bp.Length.Value = length;
            Connect(drain, gate, source, bulk);
        }
    }
}
