using System;
using SpiceSharp.Components;
using SpiceSharpParser;
using SpiceSharpParser.ModelReaders.Netlist.Spice.Readers.EntityGenerators.Components.Semiconductors;
using SpiceSharpParser.ModelReaders.Netlist.Spice.Readers.EntityGenerators.Models;

namespace SpiceSharpBSIM.Parser;

/// <summary>
/// A helper class that can extend Spice#.Parser with BSIM1 models.
/// </summary>
public static class BSIM2Extensions
{
    /// <summary>
    /// Extends a <see cref="SpiceSharpReader"/> with the BSIM1 models.
    /// </summary>
    /// <param name="reader">The reader.</param>
    /// <param name="level">The level.</param>
    /// <exception cref="NotImplementedException"></exception>
    public static void UseBSIM2(this SpiceSharpReader reader, int level = 5)
    {
        // Register the model level
        var nmosModelGenerator = reader.Settings.Mappings.Models.GetValue("NMOS", false);
        var pmosModelGenerator = reader.Settings.Mappings.Models.GetValue("PMOS", false);
        if (ReferenceEquals(nmosModelGenerator, pmosModelGenerator) &&
            nmosModelGenerator is MosfetModelGenerator modelGenerator)
            modelGenerator.AddGenericLevel<BSIM2Model, Components.Semiconductors.BSIM.BSIM2Behaviors.ModelParameters>(level);
        else
            throw new NotImplementedException();

        // Register the mosfet component
        if (reader.Settings.Mappings.Components.TryGetValue("M", false, out var generator) &&
            generator is MosfetGenerator mosGenerator)
            mosGenerator.AddMosfet<BSIM2Model, BSIM2>();
        else
            throw new NotImplementedException();
    }
}