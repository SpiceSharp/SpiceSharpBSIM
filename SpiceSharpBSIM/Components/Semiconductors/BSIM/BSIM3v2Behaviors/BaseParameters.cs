using SpiceSharp;
using SpiceSharp.Attributes;
using SpiceSharp.ParameterSets;

namespace SpiceSharpBSIM.Components.Semiconductors.BSIM.BSIM3v2Behaviors
{
    /// <summary>
    /// Parameters for a <see cref="BSIM3v2"/>.
    /// </summary>
    [GeneratedParameters]
    public partial class BaseParameters : ParameterSet<BaseParameters>
    {
        [ParameterName("l"), ParameterInfo("Length")]
        [Finite]
        private GivenParameter<double> _l = new GivenParameter<double>(5.0e-6);

        [ParameterName("w"), ParameterInfo("Width")]
        [Finite]
        private GivenParameter<double> _w = new GivenParameter<double>(5.0e-6);

        [ParameterName("m"), ParameterInfo("Parallel multiplier")]
        [Finite]
        private GivenParameter<double> _m = new GivenParameter<double>(1);

        [ParameterName("ad"), ParameterInfo("Drain area")]
        [Finite]
        private GivenParameter<double> _drainArea = new GivenParameter<double>(0.0);

        [ParameterName("as"), ParameterInfo("Source area")]
        [Finite]
        private GivenParameter<double> _sourceArea = new GivenParameter<double>(0.0);

        [ParameterName("pd"), ParameterInfo("Drain perimeter")]
        [Finite]
        private GivenParameter<double> _drainPerimeter = new GivenParameter<double>(0.0);

        [ParameterName("ps"), ParameterInfo("Source perimeter")]
        [Finite]
        private GivenParameter<double> _sourcePerimeter = new GivenParameter<double>(0.0);

        [ParameterName("nrd"), ParameterInfo("Number of squares in drain")]
        [Finite]
        private GivenParameter<double> _drainSquares = new GivenParameter<double>();

        [ParameterName("nrs"), ParameterInfo("Number of squares in source")]
        [Finite]
        private GivenParameter<double> _sourceSquares = new GivenParameter<double>();

        [ParameterName("off"), ParameterInfo("Device is initially off")]
        private bool _off;

        [ParameterName("nqsmod"), ParameterInfo("Non-quasi-static model selector")]
        private GivenParameter<int> _nqsMod = new GivenParameter<int>();

        [ParameterName("geo"), ParameterInfo("ACM model drain/source connection")]
        private GivenParameter<int> _geo = new GivenParameter<int>(0);

        [ParameterName("delvto"), ParameterInfo("Zero bias threshold voltage variation")]
        [Finite]
        private GivenParameter<double> _delvto = new GivenParameter<double>(0.0);

        [ParameterName("mulu0"), ParameterInfo("Low field mobility multiplier")]
        [Finite]
        private GivenParameter<double> _mulu0 = new GivenParameter<double>(1.0);

        [ParameterName("icvbs"), ParameterInfo("The initial base-source voltage")]
        [Finite]
        private GivenParameter<double> _icVBS = new GivenParameter<double>();

        [ParameterName("icvds"), ParameterInfo("The initial drain-source voltage")]
        [Finite]
        private GivenParameter<double> _icVDS = new GivenParameter<double>();

        [ParameterName("icvgs"), ParameterInfo("The initial gate-source voltage")]
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
