using System;
using SpiceSharp.Attributes;
namespace SpiceSharp.Components.BSIM3v24Behaviors
{
	
	/// <summary>
	/// Base parameters for a <see cref="BSIM3v24" />
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
		[ParameterName("as"), ParameterInfo("Source area")]
		public GivenParameter<double> SourceArea { get; } = new GivenParameter<double>();
		[ParameterName("ad"), ParameterInfo("Drain area")]
		public GivenParameter<double> DrainArea { get; } = new GivenParameter<double>();
		[ParameterName("ps"), ParameterInfo("Source perimeter")]
		public GivenParameter<double> SourcePerimeter { get; } = new GivenParameter<double>();
		[ParameterName("pd"), ParameterInfo("Drain perimeter")]
		public GivenParameter<double> DrainPerimeter { get; } = new GivenParameter<double>();
		[ParameterName("nrs"), ParameterInfo("Number of squares in source")]
		public GivenParameter<double> SourceSquares { get; } = new GivenParameter<double>(1);
		[ParameterName("nrd"), ParameterInfo("Number of squares in drain")]
		public GivenParameter<double> DrainSquares { get; } = new GivenParameter<double>(1);
        [ParameterName("m"), ParameterInfo("Multiplier")]
	    public GivenParameter<double> Multiplier { get; } = new GivenParameter<double>(1.0);
	    [ParameterName("off"), ParameterInfo("Device is initially off")]
		public bool Off { get; set; }
		public GivenParameter<double> IcVBS { get; } = new GivenParameter<double>();
		public GivenParameter<double> IcVDS { get; } = new GivenParameter<double>();
		public GivenParameter<double> IcVGS { get; } = new GivenParameter<double>();
		[ParameterName("nqsmod"), ParameterInfo("Non-quasi-static model selector")]
		public GivenParameter<int> NqsMod { get; } = new GivenParameter<int>();
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
					throw new ArgumentException(nameof(value));
			}
		}
	}
}