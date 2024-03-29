﻿using SpiceSharp.Entities;
using SpiceSharpBSIM.Components.Semiconductors.BSIM.BSIM3v1Behaviors;
using SpiceSharp.Attributes;

namespace SpiceSharp.Components
{
    /// <summary>
    /// A model for a BSIM3.0 transistor.
    /// </summary>
    [AutoGeneratedBehaviors]
    public partial class BSIM3v1Model : Entity<ModelParameters>
    {
        /// <summary>
        /// Creates a new <see cref="BSIM3v0Model"/>.
        /// </summary>
        /// <param name="name">The name of the model.</param>
        public BSIM3v1Model(string name)
            : base(name)
        {
        }
    }
}
