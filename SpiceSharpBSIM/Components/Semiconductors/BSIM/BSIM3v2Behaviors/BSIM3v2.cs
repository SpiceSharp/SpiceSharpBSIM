﻿using SpiceSharp.Attributes;
using SpiceSharpBSIM.Components.Semiconductors.BSIM.BSIM3v2Behaviors;

namespace SpiceSharp.Components
{
    /// <summary>
    /// BSIM3v2 model transistor
    /// </summary>
    [Pin(0, "drain"), Pin(1, "gate"), Pin(2, "source"), Pin(3, "bulk"), Connected(0, 2)]
    [AutoGeneratedBehaviors]
    public partial class BSIM3v2 : Component<BaseParameters>
    {
        /// <summary>
        /// Number of pins
        /// </summary>
        [ParameterName("pincount"), ParameterInfo("Number of pins")]
        public const int BSIM3v2PinCount = 4;

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="name">Name</param>
        public BSIM3v2(string name)
            : base(name, BSIM3v2PinCount)
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
        public BSIM3v2(string name, string drain, string gate, string source, string bulk)
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
        public BSIM3v2(string name, string drain, string gate, string source, string bulk, double width, double length)
            : this(name)
        {
            Parameters.W = width;
            Parameters.L = length;
            Connect(drain, gate, source, bulk);
        }
    }
}
