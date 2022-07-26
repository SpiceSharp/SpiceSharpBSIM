using SpiceSharp.Attributes;
using SpiceSharp.ParameterSets;
using SpiceSharp.Components;
using SpiceSharp;

namespace SpiceSharpBSIM.Components.Semiconductors.BSIM.BSIM3v1Behaviors
{
    /// <summary>
    /// Parameters for a <see cref="BSIM3v1"/>.
    /// </summary>
    [GeneratedParameters]
    public partial class BaseParameters : ParameterSet<BaseParameters>
    {
        [ParameterName("l"), ParameterInfo("Length")]
        [Finite]
        private GivenParameter<double> _l = new GivenParameter<double>(5e-6);

        [ParameterName("w"), ParameterInfo("Width")]
        [Finite]
        private GivenParameter<double> _w = new GivenParameter<double>(5e-6);

        [ParameterName("m"), ParameterInfo("Parallel multiplier")]
        [Finite]
        private GivenParameter<double> _m = new GivenParameter<double>(1);

        [ParameterName("ad"), ParameterInfo("Drain area")]
        [Finite]
        private GivenParameter<double> _drainArea = new GivenParameter<double>();

        [ParameterName("as"), ParameterInfo("Source area")]
        [Finite]
        private GivenParameter<double> _sourceArea = new GivenParameter<double>();

        [ParameterName("pd"), ParameterInfo("Drain perimeter")]
        [Finite]
        private GivenParameter<double> _drainPerimeter = new GivenParameter<double>();

        [ParameterName("ps"), ParameterInfo("Source perimeter")]
        [Finite]
        private GivenParameter<double> _sourcePerimeter = new GivenParameter<double>();

        [ParameterName("nrd"), ParameterInfo("Number of squares in drain")]
        [Finite]
        private GivenParameter<double> _drainSquares = new GivenParameter<double>(1.0);

        [ParameterName("nrs"), ParameterInfo("Number of squares in source")]
        [Finite]
        private GivenParameter<double> _sourceSquares = new GivenParameter<double>(1);

        [ParameterName("off"), ParameterInfo("Device is initially off")]
        private bool _off;

        [ParameterName("nqsmod"), ParameterInfo("Non-quasi-static model selector")]
        private GivenParameter<int> _nqsMod = new GivenParameter<int>();



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
