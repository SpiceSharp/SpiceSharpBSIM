using NUnit.Framework;
using SpiceSharpParser;
using SpiceSharpBSIM.Parser;

namespace SpiceSharpBSIMTests.BSIM1Tests;

[TestFixture]
public class BSIM1ParserTests
{
    [Test]
    public void When_BSIM1Netlist_Expect_Reference()
    {
        string netlistContent = """
            V1 g 0 0
            V2 d 0 0
            M1 d g 0 0 mod w=100u l=100u
            .MODEL mod NMOS(LEVEL=4 temp=25 muz=600 vdd=5 vfb=-0.3 phi=0.6 k1=0.5 u0=670 x2e=-0.07 mus=1082 n0=0.5 tox=1e-7 mj=0.5 mjsw=0.33 pb=0.8 pbsw=1.0 xpart=1.0)
            .DC V2 0 3.3 0.3 V1 0 3.3 0.3

            .END
            """;

        var parser = new SpiceNetlistParser();
        parser.Settings.Lexing.HasTitle = false;
        var parseResult = parser.ParseNetlist(netlistContent);

        // Convert to Spice#
        var spiceSharpReader = new SpiceSharpReader();
        spiceSharpReader.Settings.CaseSensitivity.IsModelTypeCaseSensitive = false;
        spiceSharpReader.UseBSIM1();
        var spiceSharpModel = spiceSharpReader.Read(parseResult.FinalModel);

        Assert.That(spiceSharpModel.ValidationResult.HasError || spiceSharpModel.ValidationResult.HasWarning, Is.False);
        Assert.That(spiceSharpModel.Simulations.Count, Is.EqualTo(1));

        var simulation = spiceSharpModel.Simulations[0];
        var codes = simulation.Run(spiceSharpModel.Circuit, -1);
        foreach (int _ in simulation.InvokeEvents(codes)) ;
    }
}
