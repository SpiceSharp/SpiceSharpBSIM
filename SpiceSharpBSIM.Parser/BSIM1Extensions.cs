using System;
using SpiceSharp.Components;
using SpiceSharpParser;
using SpiceSharpParser.ModelReaders.Netlist.Spice.Readers.EntityGenerators.Components.Semiconductors;
using SpiceSharpParser.ModelReaders.Netlist.Spice.Readers.EntityGenerators.Models;

namespace SpiceSharpBSIM.Parser;

/// <summary>
/// A helper class that can extend Spice#.Parser with BSIM1 models.
/// </summary>
public static class BSIM1Extensions
{
    /// <summary>
    /// Extends a <see cref="SpiceSharpReader"/> with the BSIM1 models.
    /// </summary>
    /// <param name="reader">The reader.</param>
    /// <param name="level">The level.</param>
    /// <exception cref="NotImplementedException"></exception>
    public static void UseBSIM1(this SpiceSharpReader reader, int level = 4)
    {
        // Register the model level
        var nmosModelGenerator = reader.Settings.Mappings.Models.GetValue("NMOS", false);
        var pmosModelGenerator = reader.Settings.Mappings.Models.GetValue("PMOS", false);
        if (ReferenceEquals(nmosModelGenerator, pmosModelGenerator) &&
            nmosModelGenerator is MosfetModelGenerator modelGenerator)
            modelGenerator.AddGenericLevel<BSIM1Model, Components.Semiconductors.BSIM.BSIM1Behaviors.ModelParameters>(level);
        else
            throw new NotImplementedException();

        // Register the mosfet component
        if (reader.Settings.Mappings.Components.TryGetValue("M", false, out var generator) &&
            generator is MosfetGenerator mosGenerator)
            mosGenerator.AddMosfet<BSIM1Model, BSIM1>();
        else
            throw new NotImplementedException();
    }
}