using System;
using SpiceSharp;
using SpiceSharp.Components;

namespace SpiceSharpBSIM.Parser
{
    /// <summary>
    /// A modified mosfet model generator that adds the BSIM components
    /// </summary>
    public class MosfetGenerator : SpiceSharpParser.ModelsReaders.Netlist.Spice.Readers.EntityGenerators.Components.Semiconductors.MosfetGenerator
    {
        /// <summary>
        /// Constructor
        /// </summary>
        public MosfetGenerator()
        {
            AddGenerator("BSIM1Model", "BSIM1");
            AddGenerator("BSIM2Model", "BSIM2");
            AddGenerator("BSIM3Model", "BSIM3");
            AddGenerator("BSIM3v24Model", "BSIM3v24");
        }

        /// <summary>
        /// Add a generator by class names
        /// </summary>
        /// <param name="model">Model class name</param>
        /// <param name="component">Component class name</param>
        private void AddGenerator(string model, string component)
        {
            var modelType = Type.GetType("SpiceSharp.Components." + model, false);
            var bsim1Type = Type.GetType("SpiceSharp.Components." + component, false);
            if (modelType != null && bsim1Type != null)
            {
                var method = bsim1Type.GetMethod("SetModel") ?? throw new CircuitException("The found {0} class doesn't work with this method".FormatString(component));
                Mosfets.Add(modelType, (id, entity) =>
                {
                    var c = (Component) Activator.CreateInstance(bsim1Type, id);
                    method.Invoke(c, new object[] {entity});
                    return c;
                });
            }
        }
    }
}
