using SpiceSharp.ParameterSets;
using System;
using System.Collections.Generic;
using System.Text;
using SpiceSharp.Components;
using SpiceSharp.Attributes;
using SpiceSharp;

namespace SpiceSharpBSIM.Components.Semiconductors.BSIM.BSIM4Behaviors
{
    /// <summary>
    /// Base parameters for a <see cref="BSIM4"/>.
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

        [ParameterName("m"), ParameterInfo("Separate Parallel multiplier")]
        [Finite]
        private GivenParameter<double> _m = new GivenParameter<double>(1.0);

        [ParameterName("nf"), ParameterInfo("Number of fingers")]
        [Finite]
        private GivenParameter<double> _nf = new GivenParameter<double>(1.0);

        [ParameterName("sa"), ParameterInfo("distance between  OD edge to poly of one side ")]
        [Finite]
        private GivenParameter<double> _sa = new GivenParameter<double>(0.0);

        [ParameterName("sb"), ParameterInfo("distance between  OD edge to poly of the other side")]
        [Finite]
        private GivenParameter<double> _sb = new GivenParameter<double>(0.0);

        [ParameterName("sd"), ParameterInfo("distance between neighbour fingers")]
        [Finite]
        private GivenParameter<double> _sd = new GivenParameter<double>();

        [ParameterName("sca"), ParameterInfo("Integral of the first distribution function for scattered well dopant")]
        [Finite]
        private GivenParameter<double> _sca = new GivenParameter<double>(0.0);

        [ParameterName("scb"), ParameterInfo("Integral of the second distribution function for scattered well dopant")]
        [Finite]
        private GivenParameter<double> _scb = new GivenParameter<double>(0.0);

        [ParameterName("scc"), ParameterInfo("Integral of the third distribution function for scattered well dopant")]
        [Finite]
        private GivenParameter<double> _scc = new GivenParameter<double>(0.0);

        [ParameterName("sc"), ParameterInfo("Distance to a single well edge ")]
        [Finite]
        private GivenParameter<double> _sc = new GivenParameter<double>(0.0);

        [ParameterName("min"), ParameterInfo("Minimize either D or S")]
        private GivenParameter<int> _min = new GivenParameter<int>(0);

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
        private GivenParameter<double> _drainSquares = new GivenParameter<double>(1.0);

        [ParameterName("nrs"), ParameterInfo("Number of squares in source")]
        [Finite]
        private GivenParameter<double> _sourceSquares = new GivenParameter<double>(1.0);

        [ParameterName("off"), ParameterInfo("Device is initially off")]
        private bool _off;

        [ParameterName("rbdb"), ParameterInfo("Body resistance")]
        [Finite]
        private GivenParameter<double> _rbdb = new GivenParameter<double>();

        [ParameterName("rbsb"), ParameterInfo("Body resistance")]
        [Finite]
        private GivenParameter<double> _rbsb = new GivenParameter<double>();

        [ParameterName("rbpb"), ParameterInfo("Body resistance")]
        [Finite]
        private GivenParameter<double> _rbpb = new GivenParameter<double>();

        [ParameterName("rbps"), ParameterInfo("Body resistance")]
        [Finite]
        private GivenParameter<double> _rbps = new GivenParameter<double>();

        [ParameterName("rbpd"), ParameterInfo("Body resistance")]
        [Finite]
        private GivenParameter<double> _rbpd = new GivenParameter<double>();

        [ParameterName("delvto"), ParameterName("delvt0"), ParameterInfo("Zero bias threshold voltage variation")]
        [Finite]
        private GivenParameter<double> _delvto = new GivenParameter<double>(0.0);

        [ParameterName("mulu0"), ParameterInfo("Low field mobility multiplier")]
        [Finite]
        private GivenParameter<double> _mulu0 = new GivenParameter<double>(1.0);

        [ParameterName("xgw"), ParameterInfo("Distance from gate contact center to device edge")]
        [Finite]
        private GivenParameter<double> _xgw = new GivenParameter<double>();

        [ParameterName("ngcon"), ParameterInfo("Number of gate contacts")]
        [Finite]
        private GivenParameter<double> _ngcon = new GivenParameter<double>();

        [ParameterName("wnflag"), ParameterInfo("W/NF device flag for bin selection")]
        private GivenParameter<int> _wnflag = new GivenParameter<int>();

        [ParameterName("trnqsmod"), ParameterInfo("Transient NQS model selector")]
        private GivenParameter<int> _trnqsMod = new GivenParameter<int>();

        [ParameterName("acnqsmod"), ParameterInfo("AC NQS model selector")]
        private GivenParameter<int> _acnqsMod = new GivenParameter<int>();

        [ParameterName("rbodymod"), ParameterInfo("Distributed body R model selector")]
        private GivenParameter<int> _rbodyMod = new GivenParameter<int>();

        [ParameterName("rgatemod"), ParameterInfo("Gate resistance model selector")]
        private GivenParameter<int> _rgateMod = new GivenParameter<int>();

        [ParameterName("geomod"), ParameterInfo("Geometry dependent parasitics model selector")]
        private GivenParameter<int> _geoMod = new GivenParameter<int>();

        [ParameterName("rgeomod"), ParameterInfo("S/D resistance and contact model selector")]
        private GivenParameter<int> _rgeoMod = new GivenParameter<int>();

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
