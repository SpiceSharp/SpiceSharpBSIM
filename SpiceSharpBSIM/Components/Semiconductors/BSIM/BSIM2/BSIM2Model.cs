using SpiceSharp.Components.BSIM2Behaviors;

namespace SpiceSharp.Components
{
    /// <summary>
    /// Model for the <see cref="BSIM2"/> component
    /// </summary>
    public class BSIM2Model : Model
    {
        static BSIM2Model()
        {
            RegisterBehaviorFactory(typeof(BSIM2Model), new Behaviors.BehaviorFactoryDictionary()
            {
                { typeof(ModelTemperatureBehavior), e => new ModelTemperatureBehavior(e.Name) }
            });
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="name">Name</param>
        public BSIM2Model(string name) : base(name)
        {
            ParameterSets.Add(new ModelBaseParameters());
        }
    }
}