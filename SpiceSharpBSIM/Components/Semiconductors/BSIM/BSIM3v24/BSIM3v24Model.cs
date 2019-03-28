using SpiceSharp.Components.BSIM3v24Behaviors;

namespace SpiceSharp.Components
{
    /// <summary>
    /// Model for a <see cref="BSIM3v24"/>
    /// </summary>
    public class BSIM3v24Model : Model
    {
        static BSIM3v24Model()
        {
            RegisterBehaviorFactory(typeof(BSIM3v24Model), new Behaviors.BehaviorFactoryDictionary()
            {
                { typeof(ModelTemperatureBehavior), e => new ModelTemperatureBehavior(e.Name) }
            });
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="name"></param>
        public BSIM3v24Model(string name) : base(name)
        {
            ParameterSets.Add(new ModelBaseParameters());
        }
    }
}
