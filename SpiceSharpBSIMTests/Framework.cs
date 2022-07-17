using System;
using System.Numerics;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;
using SpiceSharp;
using SpiceSharp.Simulations;
using NUnit.Framework;
using SpiceSharp.Entities;
using SpiceSharp.Algebra;

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
        protected void ApplyParameters(IEntity entity, string definition)
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
        protected void AnalyzeOp(OP sim, Circuit ckt, IEnumerable<IExport<double>> exports, IEnumerable<double> references)
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
        protected void AnalyzeDC(DC sim, Circuit ckt, IEnumerable<IExport<double>> exports, IEnumerable<double[]> references)
        {
            if (exports == null)
                throw new ArgumentNullException(nameof(exports));
            if (references == null)
                throw new ArgumentNullException(nameof(references));

            int index = 0;
            void LogWarnings(object sender, WarningEventArgs e)
            {
                Console.WriteLine(e.Message);
            }
            try
            {
                SpiceSharpWarning.WarningGenerated += LogWarnings;
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
                                throw new Exception($"{ex.Message} at {BuildMessage(sim)}", ex);
                            }
                        }

                        index++;
                    }
                };
                sim.Run(ckt);
            }
            finally
            {
                SpiceSharpWarning.WarningGenerated -= LogWarnings;
            }
        }

        private string BuildMessage(DC sim)
        {
            var values = sim.GetCurrentSweepValue();
            var names = sim.DCParameters.Sweeps.Select(s => s.Name).ToList();
            return string.Join(" ", names.Zip(values, (name, value) => $"{name} = {value}"));
        }

        /// <summary>
        /// Perform a test for AC analysis
        /// </summary>
        /// <param name="sim">Simulation</param>
        /// <param name="ckt">Circuit</param>
        /// <param name="exports">Exports</param>
        /// <param name="references">References</param>
        protected void AnalyzeAC(AC sim, Circuit ckt, IEnumerable<IExport<double>> exports, IEnumerable<double[]> references)
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
        protected void AnalyzeAC(AC sim, Circuit ckt, IEnumerable<IExport<Complex>> exports, IEnumerable<Complex[]> references)
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
        protected void AnalyzeAC(AC sim, Circuit ckt, IEnumerable<IExport<Complex>> exports, IEnumerable<Func<double, Complex>> references)
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
        protected void AnalyzeTransient(Transient sim, Circuit ckt, IEnumerable<IExport<double>> exports, IEnumerable<double[]> references)
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
        protected void AnalyzeTransient(Transient sim, Circuit ckt, IEnumerable<IExport<double>> exports, IEnumerable<Func<double, double>> references)
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
        /// Perform a test for noise analysis
        /// </summary>
        /// <param name="sim">Simulation</param>
        /// <param name="ckt">Circuit</param>
        /// <param name="exports">Exports</param>
        /// <param name="references">References</param>
        protected void AnalyzeNoise(Noise sim, Circuit ckt, IEnumerable<IExport<double>> exports, IEnumerable<double[]> references)
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
        protected void MatlabTransient(Transient sim, Circuit ckt, IEnumerable<IExport<double>> exports,
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

        /// <summary>
        /// Prints the current solver to the output stream. If <c>null</c>, the console output is used instead.
        /// </summary>
        /// <param name="state">The biasing simulation state.</param>
        public string PrintSolver(IBiasingSimulationState state)
        {
            var writer = new StringWriter();

            // Make a list of all our variables
            var variables = new string[state.Solver.Size + 1];
            var columnWidths = new int[state.Solver.Size + 1];
            var leadWidth = 0;
            foreach (var p in state.Map)
            {
                if (p.Key.Unit == Units.Volt)
                    variables[p.Value] = @"V({0})".FormatString(p.Key.Name);
                else if (p.Key.Unit == Units.Ampere)
                    variables[p.Value] = @"I({0})".FormatString(p.Key.Name);
                else
                    variables[p.Value] = @"?({0})".FormatString(p.Key.Name);
                columnWidths[p.Value] = Math.Max(variables[p.Value].Length, 6);
                leadWidth = Math.Max(leadWidth, columnWidths[p.Value]);
            }

            // Write our columns
            writer.Write(new string(' ', leadWidth + 1));
            for (var i = 1; i < variables.Length; i++)
                writer.Write($"{{0,{columnWidths[i] + 1}}}".FormatString(variables[i]));
            writer.WriteLine();

            // Write each row in the solver
            for (var row = 1; row < variables.Length; row++)
            {
                writer.Write($"{{0,{leadWidth + 1}}}".FormatString(variables[row]));
                for (var col = 1; col < variables.Length; col++)
                {
                    var elt = state.Solver.FindElement(new MatrixLocation(row, col));
                    var value = elt?.Value.ToString("g6") ?? ".";
                    writer.Write($"{{0,{columnWidths[col] + 1}:g}}".FormatString(value));
                }
                var rhsElt = state.Solver.FindElement(row);
                var rhsValue = rhsElt?.Value.ToString("g6") ?? ".";
                writer.WriteLine("{0,7}".FormatString(rhsValue));
            }

            return writer.ToString();
        }
    }
}