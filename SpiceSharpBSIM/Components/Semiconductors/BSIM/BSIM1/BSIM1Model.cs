﻿using SpiceSharp.Components.BSIM1Behaviors;

namespace SpiceSharp.Components
{
    /// <summary>
    /// Model for the <see cref="BSIM1"/> component
    /// </summary>
    public class BSIM1Model : Model
    {
        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="name">Name</param>
        public BSIM1Model(Identifier name) : base(name)
        {
            // Add parameters
            ParameterSets.Add(new ModelBaseParameters());

            // Add behaviors
            Behaviors.Add(typeof(ModelTemperatureBehavior), () => new ModelTemperatureBehavior(Name));
        }
    }
}
