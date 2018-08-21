using System;
using SpiceSharp;
using SpiceSharp.Components;

namespace SpiceSharpBSIM.Parser
{
    /// <summary>
    /// A modified mosfet model generator that adds the BSIM components
    /// </summary>
    public class MosfetModelGenerator : SpiceSharpParser.ModelsReaders.Netlist.Spice.Readers.EntityGenerators.Models.MosfetModelGenerator
    {
        /// <summary>
        /// Constructor
        /// </summary>
        public MosfetModelGenerator()
        {
            // Level 4 = BSIM1
            var bsim1Class = Type.GetType("SpiceSharp.Components.BSIM1Model", false);
            if (bsim1Class != null && bsim1Class.BaseType == typeof(Model))
            {
                Levels.Add(4, (name, type, version) =>
                {
                    var model = (Model) Activator.CreateInstance(bsim1Class, name);
                    switch (type)
                    {
                        case "nmos":
                            model.SetParameter("nmos", true);
                            break;
                        case "pmos":
                            model.SetParameter("pmos", true);
                            break;
                    }

                    return model;
                });
            }

            // Level 5 = BSIM2
            var bsim2Class = Type.GetType("SpiceSharp.Components.BSIM2Model", false);
            if (bsim2Class != null && bsim2Class.BaseType == typeof(Model))
            {
                Levels.Add(5, (name, type, version) =>
                {
                    var model = (Model) Activator.CreateInstance(bsim2Class, name);
                    switch (type)
                    {
                        case "nmos":
                            model.SetParameter("nmos", true);
                            break;
                        case "pmos":
                            model.SetParameter("pmos", true);
                            break;
                    }

                    return model;
                });
            }

            // Level x = BSIM3
            var bsim3V30Class = Type.GetType("SpiceSharp.Components.BSIM3Model", false);
            var bsim3V24Class = Type.GetType("SpiceSharp.Components.BSIM3v24Model", false);
            if (bsim3V30Class != null && bsim3V30Class.BaseType != typeof(Model))
                bsim3V30Class = null;
            if (bsim3V24Class != null && bsim3V24Class.BaseType != typeof(Model))
                bsim3V24Class = null;
            if (bsim3V24Class != null || bsim3V30Class != null)
            {
                Levels.Add(49, (name, type, version) =>
                {
                    Model model = null;
                    switch (version)
                    {
                        case null:
                        case "":
                        case "3.3.0":
                        case "3.30":
                        case "3.3":
                            if (bsim3V30Class != null)
                                model = (Model) Activator.CreateInstance(bsim3V30Class, name);
                            break;
                        case "3.2.4":
                        case "3.24":
                            if (bsim3V24Class != null)
                                model = (Model) Activator.CreateInstance(bsim3V24Class, name);
                            break;
                    }

                    if (model == null)
                        throw new CircuitException("Invalid BSIM3 model");

                    switch (type)
                    {
                        case "nmos":
                            model.SetParameter("nmos", true);
                            break;
                        case "pmos":
                            model.SetParameter("pmos", true);
                            break;
                    }

                    return model;
                });
            }
        }
    }
}
