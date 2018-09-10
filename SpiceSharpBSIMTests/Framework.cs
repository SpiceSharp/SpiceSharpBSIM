using System;
using System.Numerics;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;
using SpiceSharp;
using SpiceSharp.Simulations;
using SpiceSharp.Circuits;
using NUnit.Framework;

namespace SpiceSharpTest.Models
{
    /// <summary>
    /// Framework for testing models
    /// </summary>
    public class Framework
    {
        /// <summary>
        /// Absolute tolerance used
        /// </summary>
        public double AbsTol = 1e-12;

        /// <summary>
        /// Relative tolerance used
        /// </summary>
        public double RelTol = 1e-3;

        /// <summary>
        /// Apply a parameter definition to an entity
        /// Parameters are a series of assignments [name]=[value] delimited by spaces.
        /// </summary>
        /// <param name="entity">Entity</param>
        /// <param name="definition">Definition string</param>
        protected void ApplyParameters(Entity entity, string definition)
        {
            // Get all assignments
            definition = Regex.Replace(definition, @"\s*\=\s*", "=");
            string[] assignments = definition.Split(new[] { ',', ';', ' ' }, StringSplitOptions.RemoveEmptyEntries);
            foreach (var assignment in assignments)
            {
                // Get the name and value
                string[] parts = assignment.Split('=');
                if (parts.Length != 2)
                    throw new Exception("Invalid assignment");
                string name = parts[0].ToLower();
                double value = double.Parse(parts[1], System.Globalization.CultureInfo.InvariantCulture);

                // Set the entity parameter
                entity.SetParameter(name, value);
            }
        }

        /// <summary>
        /// Perform a test for OP analysis
        /// </summary>
        /// <param name="sim">Simulation</param>
        /// <param name="ckt">Circuit</param>
        /// <param name="exports">Exports</param>
        /// <param name="references">References</param>
        protected void AnalyzeOp(OP sim, Circuit ckt, IEnumerable<Export<double>> exports, IEnumerable<double> references)
        {
            if (exports == null)
                throw new ArgumentNullException(nameof(exports));
            if (references == null)
                throw new ArgumentNullException(nameof(references));

            sim.ExportSimulationData += (sender, data) =>
            {
                using (var exportIt = exports.GetEnumerator())
                using (var referencesIt = references.GetEnumerator())
                {
                    while (exportIt.MoveNext() && referencesIt.MoveNext())
                    {
                        double actual = exportIt.Current?.Value ?? throw new ArgumentNullException();
                        double expected = referencesIt.Current;
                        double tol = Math.Max(Math.Abs(actual), Math.Abs(expected)) * RelTol + AbsTol;
                        Assert.AreEqual(expected, actual, tol);
                    }
                }
            };
            sim.Run(ckt);
        }

        /// <summary>
        /// Perform a test for DC analysis
        /// </summary>
        /// <param name="sim">Simulation</param>
        /// <param name="ckt">Circuit</param>
        /// <param name="exports">Exports</param>
        /// <param name="references">References</param>
        protected void AnalyzeDC(DC sim, Circuit ckt, IEnumerable<Export<double>> exports, IEnumerable<double[]> references)
        {
            if (exports == null)
                throw new ArgumentNullException(nameof(exports));
            if (references == null)
                throw new ArgumentNullException(nameof(references));

            int index = 0;
            sim.ExportSimulationData += (sender, data) =>
            {
                using (var exportIt = exports.GetEnumerator())
                using (var referencesIt = references.GetEnumerator())
                {
                    while (exportIt.MoveNext() && referencesIt.MoveNext())
                    {
                        double actual = exportIt.Current?.Value ?? double.NaN;
                        double expected = referencesIt.Current?[index] ?? double.NaN;
                        double tol = Math.Max(Math.Abs(actual), Math.Abs(expected)) * RelTol + AbsTol;

                        try
                        {
                            Assert.AreEqual(expected, actual, tol);
                        }
                        catch (Exception ex)
                        {
                            string[] sweeps = new string[sim.Sweeps.Count];
                            for (int k = 0; k < sim.Sweeps.Count; k++)
                                sweeps[k] += $"{sim.Sweeps[k].Parameter}={sim.Sweeps[k].CurrentValue}";
                            string msg = ex.Message + " at " + string.Join(" ", sweeps);
                            throw new Exception(msg, ex);
                        }
                    }

                    index++;
                }
            };
            sim.Run(ckt);
        }

        /// <summary>
        /// Perform a test for DC analysis
        /// </summary>
        /// <param name="sim">Simulation</param>
        /// <param name="ckt">Circuit</param>
        /// <param name="exports">Exports</param>
        /// <param name="references">References</param>
        protected void AnalyzeDC(DC sim, Circuit ckt, IEnumerable<Export<double>> exports, IEnumerable<Func<double, double>> references)
        {
            sim.ExportSimulationData += (sender, data) =>
            {
                using (var exportIt = exports.GetEnumerator())
                using (var referencesIt = references.GetEnumerator())
                {
                    while (exportIt.MoveNext() && referencesIt.MoveNext())
                    {
                        double actual = exportIt.Current?.Value ?? throw new ArgumentNullException();
                        double expected = referencesIt.Current?.Invoke(sim.Sweeps[0].CurrentValue) ?? throw new ArgumentNullException();
                        double tol = Math.Max(Math.Abs(actual), Math.Abs(expected)) * RelTol + AbsTol;

                        try
                        {
                            Assert.AreEqual(expected, actual, tol);
                        }
                        catch (Exception ex)
                        {
                            string[] sweeps = new string[sim.Sweeps.Count];
                            for (int k = 0; k < sim.Sweeps.Count; k++)
                                sweeps[k] += $"{sim.Sweeps[k].Parameter}={sim.Sweeps[k].CurrentValue}";
                            string msg = ex.Message + " at " + string.Join(" ", sweeps);
                            throw new Exception(msg, ex);
                        }
                    }
                }
            };
            sim.Run(ckt);
        }

        /// <summary>
        /// Perform a test for AC analysis
        /// </summary>
        /// <param name="sim">Simulation</param>
        /// <param name="ckt">Circuit</param>
        /// <param name="exports">Exports</param>
        /// <param name="references">References</param>
        protected void AnalyzeAC(AC sim, Circuit ckt, IEnumerable<Export<double>> exports, IEnumerable<double[]> references)
        {
            int index = 0;
            sim.ExportSimulationData += (sender, data) =>
            {
                using (var exportsIt = exports.GetEnumerator())
                using (var referencesIt = references.GetEnumerator())
                {
                    while (exportsIt.MoveNext() && referencesIt.MoveNext())
                    {
                        // Test export
                        double actual = exportsIt.Current?.Value ?? throw new ArgumentNullException();
                        double expected = referencesIt.Current?[index] ?? throw new ArgumentNullException();
                        double tol = Math.Max(Math.Abs(actual), Math.Abs(expected)) * RelTol + AbsTol;

                        try
                        {
                            Assert.AreEqual(expected, actual, tol);
                        }
                        catch (Exception ex)
                        {
                            string msg = $"{ex.Message} at {data.Frequency} Hz";
                            throw new Exception(msg, ex);
                        }
                    }

                    index++;
                }
            };
            sim.Run(ckt);
        }

        /// <summary>
        /// Perform a test for AC analysis
        /// </summary>
        /// <param name="sim">Simulation</param>
        /// <param name="ckt">Circuit</param>
        /// <param name="exports">Exports</param>
        /// <param name="references">References</param>
        protected void AnalyzeAC(AC sim, Circuit ckt, IEnumerable<Export<Complex>> exports, IEnumerable<Complex[]> references)
        {
            int index = 0;
            sim.ExportSimulationData += (sender, data) =>
            {
                using (var exportsIt = exports.GetEnumerator())
                using (var referencesIt = references.GetEnumerator())
                {
                    while (exportsIt.MoveNext() && referencesIt.MoveNext())
                    {
                        // Test export
                        Complex actual = exportsIt.Current?.Value ?? throw new ArgumentNullException();
                        Complex expected = referencesIt.Current?[index] ?? throw new ArgumentNullException();

                        // Test real part
                        double rtol = Math.Max(Math.Abs(actual.Real), Math.Abs(expected.Real)) * RelTol + AbsTol;
                        double itol = Math.Max(Math.Abs(actual.Imaginary), Math.Abs(expected.Imaginary)) * RelTol +
                                      AbsTol;

                        try
                        {
                            Assert.AreEqual(expected.Real, actual.Real, rtol);
                            Assert.AreEqual(expected.Imaginary, actual.Imaginary, itol);
                        }
                        catch (Exception ex)
                        {
                            string msg = $"{ex.Message} at {data.Frequency} Hz";
                            throw new Exception(msg, ex);
                        }
                    }

                    index++;
                }
            };
            sim.Run(ckt);
        }

        /// <summary>
        /// Perform a test for AC analysis
        /// </summary>
        /// <param name="sim">Simulation</param>
        /// <param name="ckt">Circuit</param>
        /// <param name="exports">Exports</param>
        /// <param name="references">References</param>
        protected void AnalyzeAC(AC sim, Circuit ckt, IEnumerable<Export<Complex>> exports, IEnumerable<Func<double, Complex>> references)
        {
            sim.ExportSimulationData += (sender, data) =>
            {
                using (var exportsIt = exports.GetEnumerator())
                using (var referencesIt = references.GetEnumerator())
                {
                    while (exportsIt.MoveNext() && referencesIt.MoveNext())
                    {
                        // Test export
                        Complex actual = exportsIt.Current?.Value ?? throw new ArgumentNullException();
                        Complex expected = referencesIt.Current?.Invoke(data.Frequency) ?? throw new ArgumentNullException();

                        // Test real part
                        double rtol = Math.Max(Math.Abs(actual.Real), Math.Abs(expected.Real)) * RelTol + AbsTol;
                        double itol = Math.Max(Math.Abs(actual.Imaginary), Math.Abs(expected.Imaginary)) * RelTol +
                                      AbsTol;

                        try
                        {
                            Assert.AreEqual(expected.Real, actual.Real, rtol);
                            Assert.AreEqual(expected.Imaginary, actual.Imaginary, itol);
                        }
                        catch (Exception ex)
                        {
                            string msg = $"{ex.Message} at {data.Frequency} Hz";
                            throw new Exception(msg, ex);
                        }
                    }
                }
            };
            sim.Run(ckt);
        }

        /// <summary>
        /// Perform a test for transient analysis
        /// </summary>
        /// <param name="sim">Simulation</param>
        /// <param name="ckt">Circuit</param>
        /// <param name="exports">Exports</param>
        /// <param name="references">References</param>
        protected void AnalyzeTransient(Transient sim, Circuit ckt, IEnumerable<Export<double>> exports, IEnumerable<double[]> references)
        {
            int index = 0;
            sim.ExportSimulationData += (sender, data) =>
            {
                using (var exportsIt = exports.GetEnumerator())
                using (var referencesIt = references.GetEnumerator())
                {
                    while (exportsIt.MoveNext() && referencesIt.MoveNext())
                    {
                        double actual = exportsIt.Current?.Value ?? throw new ArgumentNullException();
                        double expected = referencesIt.Current?[index] ?? throw new ArgumentNullException();
                        double tol = Math.Max(Math.Abs(actual), Math.Abs(expected)) * RelTol + AbsTol;

                        try
                        {
                            Assert.AreEqual(expected, actual, tol);
                        }
                        catch (Exception ex)
                        {
                            string msg = $"{ex.Message} at t={data.Time} s";
                            throw new Exception(msg, ex);
                        }
                    }

                    index++;
                }
            };
            sim.Run(ckt);
        }

        /// <summary>
        /// Perform a test for transient analysis where the reference is a function in time
        /// </summary>
        /// <param name="sim">Simulation</param>
        /// <param name="ckt">Circuit</param>
        /// <param name="exports">Exports</param>
        /// <param name="references">References</param>
        protected void AnalyzeTransient(Transient sim, Circuit ckt, IEnumerable<Export<double>> exports, IEnumerable<Func<double, double>> references)
        {
            sim.ExportSimulationData += (sender, data) =>
            {
                using (var exportsIt = exports.GetEnumerator())
                using (var referencesIt = references.GetEnumerator())
                {
                    while (exportsIt.MoveNext() && referencesIt.MoveNext())
                    {
                        double t = data.Time;
                        double actual = exportsIt.Current?.Value ?? throw new ArgumentNullException();
                        double expected = referencesIt.Current?.Invoke(t) ?? throw new ArgumentNullException();
                        double tol = Math.Max(Math.Abs(actual), Math.Abs(expected)) * RelTol + AbsTol;

                        try
                        {
                            Assert.AreEqual(expected, actual, tol);
                        }
                        catch (Exception ex)
                        {
                            string msg = $"{ex.Message} at t={data.Time} s";
                            throw new Exception(msg, ex);
                        }
                    }
                }
            };
            sim.Run(ckt);
        }

        /// <summary>
        /// Perform a test for transient analysis where the reference is the solve matrix
        /// </summary>
        /// <param name="sim"></param>
        /// <param name="ckt"></param>
        /// <param name="path"></param>
        protected void AnalyzeTransient(Transient sim, Circuit ckt, string path)
        {
            // Get all lines
            var filename = Path.Combine(TestContext.CurrentContext.TestDirectory, path);
            var lines = System.IO.File.ReadAllText(filename).Split('\n');
            var regReference = new Regex(@"(?<row>\d+)\s+(?<col>\d+)\s+(?<value>[\+\-]?\d+(\.\d+)?([eE][\+\-]?\d+)?)");

            int index = 0;
            sim.ExportSimulationData += (sender, data) =>
            {
                // Read until the first empty line
                while (!string.IsNullOrWhiteSpace(lines[index]))
                {
                    var m = regReference.Match(lines[index]);
                    if (!m.Success)
                        throw new Exception("Could not read reference");
                    var row = int.Parse(m.Groups["row"].Value);
                    var col = int.Parse(m.Groups["col"].Value);
                    var value = double.Parse(m.Groups["value"].Value);

                    // Check the solver
                    var eltIndex = sim.RealState.Solver.InternalToExternal(new Tuple<int, int>(row, col));
                    var elt = sim.RealState.Solver.FindMatrixElement(eltIndex.Item1, eltIndex.Item2);
                    Assert.NotNull(elt);

                    double tol = Math.Max(Math.Abs(value), Math.Abs(elt.Value)) * RelTol + AbsTol;
                    Assert.AreEqual(value, elt.Value, tol);
                    index++;
                }
            };
            sim.Run(ckt);
        }

        /// <summary>
        /// Perform a test for noise analysis
        /// </summary>
        /// <param name="sim">Simulation</param>
        /// <param name="ckt">Circuit</param>
        /// <param name="exports">Exports</param>
        /// <param name="references">References</param>
        protected void AnalyzeNoise(Noise sim, Circuit ckt, IEnumerable<Export<double>> exports, IEnumerable<double[]> references)
        {
            int index = 0;
            sim.ExportSimulationData += (sender, data) =>
            {
                using (var exportsIt = exports.GetEnumerator())
                using (var referencesIt = references.GetEnumerator())
                {
                    while (exportsIt.MoveNext() && referencesIt.MoveNext())
                    {
                        // Test export
                        double actual = exportsIt.Current?.Value ?? throw new ArgumentNullException();
                        double expected = referencesIt.Current?[index] ?? throw new ArgumentNullException();
                        double tol = Math.Max(Math.Abs(actual), Math.Abs(expected)) * 1e-6 + 1e-12;

                        try
                        {
                            Assert.AreEqual(expected, actual, tol);
                        }
                        catch (Exception ex)
                        {
                            string msg = $"{ex.Message} at {data.Frequency} Hz";
                            throw new Exception(msg, ex);
                        }
                    }

                    index++;
                }
            };
            sim.Run(ckt);
        }

        /// <summary>
        /// Output for Matlab plotting
        /// </summary>
        /// <param name="sim">Simulation</param>
        /// <param name="ckt">Circuit</param>
        /// <param name="exports">Exports</param>
        /// <param name="references">References</param>
        protected void MatlabTransient(Transient sim, Circuit ckt, IEnumerable<Export<double>> exports,
            IEnumerable<double[]> references)
        {
            var refQuantity = new List<List<double>>();
            var simQuantity = new List<List<double>>();

            int index = 0;
            sim.ExportSimulationData += (sender, data) =>
            {
                using (var exportsIt = exports.GetEnumerator())
                using (var referencesIt = references.GetEnumerator())
                {
                    int cIndex = 0;
                    while (exportsIt.MoveNext() && referencesIt.MoveNext())
                    {
                        double actual = exportsIt.Current?.Value ?? throw new ArgumentNullException();
                        double expected = referencesIt.Current?[index] ?? throw new ArgumentNullException();
                        double tol = Math.Max(Math.Abs(actual), Math.Abs(expected)) * RelTol + AbsTol;

                        // Add new series
                        if (index == 0)
                        {
                            refQuantity.Add(new List<double>());
                            simQuantity.Add(new List<double>());
                        }
                        refQuantity[cIndex].Add(expected);
                        simQuantity[cIndex].Add(actual);

                        Console.Write(expected.ToString() + ", ");
                    }

                    index++;
                }
            };
            sim.Run(ckt);

            // Show the vectors
            foreach (var vec in refQuantity)
                Console.WriteLine("v = [" + string.Join(",", vec) + "];");
            foreach (var vec in simQuantity)
                Console.WriteLine("actual = [" + string.Join(",", vec) + "];");
        }
    }
}