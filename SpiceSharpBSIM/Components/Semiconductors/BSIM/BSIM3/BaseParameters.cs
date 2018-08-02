using System;
using SpiceSharp.Attributes;
namespace SpiceSharp.Components.BSIM3Behaviors
{
	
	/// <summary>
	/// Base parameters for a <see cref="BSIM3" />
	/// </summary>
	public class BaseParameters : ParameterSet
	{
		
		/// <summary>
		/// Properties
		/// </summary>
		[ParameterName("w"), ParameterInfo("Width")]
		public GivenParameter<double> Width { get; } = new GivenParameter<double>(5e-06);
		[ParameterName("l"), ParameterInfo("Length")]
		public GivenParameter<double> Length { get; } = new GivenParameter<double>(5e-06);
		[ParameterName("m"), ParameterInfo("Parallel multiplier")]
		public GivenParameter<double> M { get; } = new GivenParameter<double>(1);
		[ParameterName("as"), ParameterInfo("Source area")]
		public GivenParameter<double> SourceArea { get; } = new GivenParameter<double>();
		[ParameterName("ad"), ParameterInfo("Drain area")]
		public GivenParameter<double> DrainArea { get; } = new GivenParameter<double>();
		[ParameterName("ps"), ParameterInfo("Source perimeter")]
		public GivenParameter<double> SourcePerimeter { get; } = new GivenParameter<double>();
		[ParameterName("pd"), ParameterInfo("Drain perimeter")]
		public GivenParameter<double> DrainPerimeter { get; } = new GivenParameter<double>();
		[ParameterName("nrs"), ParameterInfo("Number of squares in source")]
		public GivenParameter<double> SourceSquares { get; } = new GivenParameter<double>();
		[ParameterName("nrd"), ParameterInfo("Number of squares in drain")]
		public GivenParameter<double> DrainSquares { get; } = new GivenParameter<double>();
		[ParameterName("off"), ParameterInfo("Device is initially off")]
		public bool Off { get; set; }
		[ParameterInfo("")]
		public GivenParameter<double> IcVBS { get; } = new GivenParameter<double>();
		[ParameterInfo("")]
		public GivenParameter<double> IcVDS { get; } = new GivenParameter<double>();
		[ParameterInfo("")]
		public GivenParameter<double> IcVGS { get; } = new GivenParameter<double>();
		[ParameterName("nqsmod"), ParameterInfo("Non-quasi-static model selector")]
		public GivenParameter<int> NqsMod { get; } = new GivenParameter<int>();
		[ParameterName("acnqsmod"), ParameterInfo("AC NQS model selector")]
		public GivenParameter<int> AcnqsMod { get; } = new GivenParameter<int>();
		[ParameterName("geo"), ParameterInfo("ACM model drain/source connection")]
		public GivenParameter<double> Geo { get; } = new GivenParameter<double>();
		[ParameterName("delvto"), ParameterInfo("Zero bias threshold voltage variation")]
		public GivenParameter<double> Delvto { get; } = new GivenParameter<double>();
		[ParameterName("mulu0"), ParameterInfo("Low field mobility multiplier")]
		public GivenParameter<double> Mulu0 { get; } = new GivenParameter<double>(1);
		[ParameterName("ic"), ParameterInfo("Vector of DS,GS,BS initial voltages")]
		public void SetIc(double[] value)
		{
		    switch (value.Length)
		    {
		        case 3:
		            IcVBS.Value = value[2];
		            goto case 2;
		        case 2:
		            IcVGS.Value = value[1];
		            goto case 1;
		        case 1:
		            IcVDS.Value = value[0];
		            break;
		        default:
		            throw new BadParameterException(nameof(value));
		    }
		}
    }
}