﻿using SpiceSharp.Attributes;
using SpiceSharp.Entities;
using SpiceSharpBSIM.Components.Semiconductors.BSIM.BSIM1Behaviors;

namespace SpiceSharp.Components
{
    /// <summary>
    /// Model for the <see cref="BSIM1"/> component
    /// </summary>
    [AutoGeneratedBehaviors]
    public partial class BSIM1Model : Entity<ModelParameters>
    {
        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="name">Name</param>
        public BSIM1Model(string name)
            : base(name)
        {
        }
    }
}
