using SpiceSharp.Components.BSIM3Behaviors;

namespace SpiceSharp.Components
{
    /// <summary>
    /// Model for the <see cref="BSIM3"/> component
    /// </summary>
    public class BSIM3Model : Model
    {
        static BSIM3Model()
        {
            RegisterBehaviorFactory(typeof(BSIM3Model), new Behaviors.BehaviorFactoryDictionary()
            {
                { typeof(ModelTemperatureBehavior), e => new ModelTemperatureBehavior(e.Name) }
            });
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="name">Name</param>
        public BSIM3Model(string name) : base(name)
        {
            ParameterSets.Add(new ModelBaseParameters());
        }
    }
}
