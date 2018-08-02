using System.Collections.Generic;
using System.Diagnostics;
using System.Numerics;
using NUnit;
using NUnit.Framework;
using SpiceSharp;
using SpiceSharp.Behaviors;
using SpiceSharp.Components;
using SpiceSharp.IntegrationMethods;
using SpiceSharp.Simulations;

namespace SpiceSharpTest.Models
{
    [TestFixture]
    public class BSIM2Tests : Framework
    {
        /// <summary>
        /// Generate a BSIM1 transistor
        /// </summary>
        private BSIM2 Create(Identifier name, Identifier drain, Identifier gate, Identifier source, Identifier bulk, double w, double l, Identifier model, string parameters)
        {
            // Create the model
            var m = new BSIM2Model(model);
            ApplyParameters(m, parameters);

            // Create the device
            var e = new BSIM2(name, drain, gate, source, bulk);
            e.SetParameter("w", w);
            e.SetParameter("l", l);
            e.SetModel(m);
            return e;
        }

        [Test]
        public void When_BSIM2DC_Expect_Reference()
        {
            var ckt = new Circuit();
            ckt.Objects.Add(
                new VoltageSource("V1", "g", "0", 0.0),
                new VoltageSource("V2", "d", "0", 0.0),
                Create("M1", "d", "g", "0", "0", 100e-6, 100e-6, "mod", "temp=25 muz=600 vdd=5 vfb=-0.3 phi=0.6 k1=0.5 u0=670 x2e=-0.07 mus=1082 n0=0.5 tox=1e-7 mj=0.5 mjsw=0.33 pb=0.8 pbsw=1.0 xpart=1.0"));

            // Create simulation
            var dc = new DC("dc", new[]
            {
                new SweepConfiguration("V1", 0, 3.3, 0.3),
                new SweepConfiguration("V2", 0, 3.3, 0.3),
            });

            // Create exports
            Export<double>[] exports = {new RealPropertyExport(dc, "V2", "i")};

            // Create references
            double[][] references =
            {
                new[]
                {
                    0.0
                }
            };

            // Run simulation
            AnalyzeDC(dc, ckt, exports, references);
        }
    }
}
