﻿using SpiceSharp.Attributes;
using SpiceSharpBSIM.Components.Semiconductors.BSIM.BSIM1Behaviors;

namespace SpiceSharp.Components
{
    /// <summary>
    /// BSIM1 model transistor
    /// </summary>
    [AutoGeneratedBehaviors]
    public partial class BSIM1 : Component<BaseParameters>
    {
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
            : this(name)
        {
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
            : this(name)
        {
            Parameters.Width = width;
            Parameters.Length = length;
            Connect(drain, gate, source, bulk);
        }
    }
}
