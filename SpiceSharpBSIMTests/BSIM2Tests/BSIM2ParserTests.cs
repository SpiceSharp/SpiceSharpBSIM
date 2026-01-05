using NUnit.Framework;
using SpiceSharpParser;
using SpiceSharpBSIM.Parser;

namespace SpiceSharpBSIMTests.BSIM2Tests;

[TestFixture]
public class BSIM2ParserTests
{
    [Test]
    public void When_BSIM2Netlist_Expect_Reference()
    {
        string netlistContent = """
            V1 g 0 0
            V2 d 0 0
            M1 d g 0 0 mod w=100u l=100u
            .MODEL mod NMOS(LEVEL=5 vfb=-0.3 phi=0.8 k1=0.6 mu0=250 n0=1.3 tox=1e-7 mj=0.5 mjsw=0.33 pb=0.8 pbsw=1.0 xpart=1.0)
            .DC V2 0 3.3 0.3 V1 0 3.3 0.3

            .END
            """;

        var parser = new SpiceNetlistParser();
        parser.Settings.Lexing.HasTitle = false;
        var parseResult = parser.ParseNetlist(netlistContent);

        // Convert to Spice#
        var spiceSharpReader = new SpiceSharpReader();
        spiceSharpReader.Settings.CaseSensitivity.IsModelTypeCaseSensitive = false;
        spiceSharpReader.UseBSIM2();
        var spiceSharpModel = spiceSharpReader.Read(parseResult.FinalModel);

        Assert.That(spiceSharpModel.ValidationResult.HasError || spiceSharpModel.ValidationResult.HasWarning, Is.False);
        Assert.That(spiceSharpModel.Simulations.Count, Is.EqualTo(1));

        var simulation = spiceSharpModel.Simulations[0];
        var codes = simulation.Run(spiceSharpModel.Circuit, -1);
        foreach (int _ in simulation.InvokeEvents(codes)) ;
    }
}
