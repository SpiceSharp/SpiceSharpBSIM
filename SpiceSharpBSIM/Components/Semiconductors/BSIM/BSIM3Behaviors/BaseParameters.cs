using SpiceSharp;
using SpiceSharp.Attributes;
using SpiceSharp.ParameterSets;

namespace SpiceSharpBSIM.Components.Semiconductors.BSIM.BSIM3Behaviors
{
    /// <summary>
    /// Base parameters for a <see cref="BSIM3" />
    /// </summary>
    [GeneratedParameters]
    public partial class BaseParameters : ParameterSet<BaseParameters>
    {
        /// <summary>
        /// Properties
        /// </summary>
        [ParameterName("w"), ParameterInfo("Width")]
        [Finite]
        private GivenParameter<double> _width = new GivenParameter<double>(5e-06);

        [ParameterName("l"), ParameterInfo("Length")]
        [Finite]
        private GivenParameter<double> _length = new GivenParameter<double>(5e-06);

        [ParameterName("m"), ParameterInfo("Parallel multiplier")]
        [Finite]
        private GivenParameter<double> _m = new GivenParameter<double>(1);

        [ParameterName("as"), ParameterInfo("Source area")]
        [Finite]
        private GivenParameter<double> _sourceArea = new GivenParameter<double>();

        [ParameterName("ad"), ParameterInfo("Drain area")]
        [Finite]
        private GivenParameter<double> _drainArea = new GivenParameter<double>();

        [ParameterName("ps"), ParameterInfo("Source perimeter")]
        [Finite]
        private GivenParameter<double> _sourcePerimeter = new GivenParameter<double>();

        [ParameterName("pd"), ParameterInfo("Drain perimeter")]
        [Finite]
        private GivenParameter<double> _drainPerimeter = new GivenParameter<double>();

        [ParameterName("nrs"), ParameterInfo("Number of squares in source")]
        [Finite]
        private GivenParameter<double> _sourceSquares = new GivenParameter<double>();

        [ParameterName("nrd"), ParameterInfo("Number of squares in drain")]
        [Finite]
        private GivenParameter<double> _drainSquares = new GivenParameter<double>();

        [ParameterName("off"), ParameterInfo("Device is initially off")]
        public bool Off { get; set; }

        [ParameterName("icvbs"), ParameterInfo("The initial base-source voltage")]
        [Finite]
        private GivenParameter<double> _icVBS = new GivenParameter<double>();

        [ParameterName("icvds"), ParameterInfo("The initial drain-source voltage")]
        [Finite]
        private GivenParameter<double> _icVDS = new GivenParameter<double>();

        [ParameterName("icvgs"), ParameterInfo("The initial gate-source voltage")]
        [Finite]
        private GivenParameter<double> _icVGS = new GivenParameter<double>();

        [ParameterName("nqsmod"), ParameterInfo("Non-quasi-static model selector")]
        [Finite]
        private GivenParameter<int> _nqsMod = new GivenParameter<int>();

        [ParameterName("acnqsmod"), ParameterInfo("AC NQS model selector")]
        [Finite]
        private GivenParameter<int> _acnqsMod = new GivenParameter<int>();

        [ParameterName("geo"), ParameterInfo("ACM model drain/source connection")]
        [Finite]
        private GivenParameter<int> _geo = new GivenParameter<int>();

        [ParameterName("delvto"), ParameterInfo("Zero bias threshold voltage variation")]
        [Finite]
        private GivenParameter<double> _delvto = new GivenParameter<double>();

        [ParameterName("mulu0"), ParameterInfo("Low field mobility multiplier")]
        [Finite]
        private GivenParameter<double> _mulu0 = new GivenParameter<double>(1);

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