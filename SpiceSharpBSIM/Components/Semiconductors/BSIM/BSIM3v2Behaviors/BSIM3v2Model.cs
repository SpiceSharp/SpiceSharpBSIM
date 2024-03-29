﻿using SpiceSharp.Entities;
using SpiceSharpBSIM.Components.Semiconductors.BSIM.BSIM3v2Behaviors;
using SpiceSharp.Attributes;

namespace SpiceSharp.Components
{
    /// <summary>
    /// A model for a BSIM3.2 transistor.
    /// </summary>
    [AutoGeneratedBehaviors]
    public partial class BSIM3v2Model : Entity<ModelParameters>
    {
        /// <summary>
        /// Creates a new <see cref="BSIM3v2Model"/>.
        /// </summary>
        /// <param name="name">The name of the model.</param>
        public BSIM3v2Model(string name)
            : base(name)
        {
        }
    }
}
