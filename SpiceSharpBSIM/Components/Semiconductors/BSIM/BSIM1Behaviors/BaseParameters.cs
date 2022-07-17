using SpiceSharp;
using SpiceSharp.Attributes;
using SpiceSharp.ParameterSets;

namespace SpiceSharpBSIM.Components.Semiconductors.BSIM.BSIM1Behaviors
{
    /// <summary>
    /// Base parameters for a <see cref="BSIM1" />
    /// </summary>
    [GeneratedParameters]
    public partial class BaseParameters : ParameterSet<BaseParameters>
    {
        [ParameterName("w"), ParameterInfo("Width")]
        [Finite, GreaterThan(0)]
        private GivenParameter<double> _width = new GivenParameter<double>(5e-6);

        [ParameterName("l"), ParameterInfo("Length")]
        [Finite, GreaterThan(0)]
        private GivenParameter<double> _length = new GivenParameter<double>(5e-6);

        [ParameterName("as"), ParameterInfo("Source area")]
        [Finite, GreaterThanOrEquals(0)]
        private GivenParameter<double> _sourceArea = new GivenParameter<double>();

        [ParameterName("ad"), ParameterInfo("Drain area")]
        [Finite, GreaterThanOrEquals(0)]
        private GivenParameter<double> _drainArea = new GivenParameter<double>();

        [ParameterName("ps"), ParameterInfo("Source perimeter")]
        [Finite, GreaterThanOrEquals(0)]
        private GivenParameter<double> _sourcePerimeter = new GivenParameter<double>();

        [ParameterName("pd"), ParameterInfo("Drain perimeter")]
        [Finite, GreaterThanOrEquals(0)]
        private GivenParameter<double> _drainPerimeter = new GivenParameter<double>();

        [ParameterName("nrs"), ParameterInfo("Number of squares in source")]
        [Finite, GreaterThan(0)]
        private GivenParameter<double> _sourceSquares = new GivenParameter<double>(1);

        [ParameterName("nrd"), ParameterInfo("Number of squares in drain")]
        [Finite, GreaterThan(0)]
        private GivenParameter<double> _drainSquares = new GivenParameter<double>(1);

        [ParameterName("off"), ParameterInfo("Device is initially off")]
        public bool Off { get; set; }

        [ParameterName("m"), ParameterInfo("Multiplier")]
        public double Multiplier { get; set; } = 1.0;

        [ParameterName("vbs"), ParameterInfo("Initial B-S voltage")]
        [Finite]
        private GivenParameter<double> _icVBS = new GivenParameter<double>();

        [ParameterName("vds"), ParameterInfo("Initial D-S voltage")]
        [Finite]
        private GivenParameter<double> _icVDS = new GivenParameter<double>();

        [ParameterName("vgs"), ParameterInfo("Initial G-S voltage")]
        [Finite]
        private GivenParameter<double> _icVGS = new GivenParameter<double>();

        [ParameterName("ic"), ParameterInfo("Vector of DS,GS,BS initial voltages")]
        public void SetIc(double[] value)
        {
            switch (value.Length)
            {
                case 3:
                    IcVBS = value[2];
                    goto case 2;
                case 2:
                    IcVGS = value[1];
                    goto case 1;
                case 1:
                    IcVDS = value[0];
                    break;
                default:
                    throw new SpiceSharpException("Invalid parameter");
            }
        }
    }
}