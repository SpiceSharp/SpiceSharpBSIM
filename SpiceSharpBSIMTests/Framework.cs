using System;
using System.Numerics;
using System.Collections.Generic;
using System.Linq;
using System.Text.RegularExpressions;
using SpiceSharp;
using SpiceSharp.Simulations;
using NUnit.Framework;
using SpiceSharp.Entities;

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
                if (!entity.TrySetParameter(name, value) && !entity.TrySetParameter(name, (int)value))
                    throw new ArgumentException("Could not find parameter " + name);
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

            foreach (int code in sim.Run(ckt))
            {
                if (code == OP.ExportOperatingPoint)
                {
                    using (var exportIt = exports.GetEnumerator())
                    using (var referencesIt = references.GetEnumerator())
                    {
                        while (exportIt.MoveNext() && referencesIt.MoveNext())
                        {
                            double actual = exportIt.Current?.Value ?? throw new ArgumentNullException();
                            double expected = referencesIt.Current;
                            double tol = Math.Max(Math.Abs(actual), Math.Abs(expected)) * RelTol + AbsTol;
                            Assert.That(actual, Is.EqualTo(actual).Within(tol));
                        }
                    }
                }
            }
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
                foreach (int code in sim.Run(ckt))
                {
                    if (code == DC.ExportSweep)
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
                                    Assert.That(actual, Is.EqualTo(expected).Within(tol));
                                }
                                catch (Exception ex)
                                {
                                    throw new Exception($"{ex.Message} at {BuildMessage(sim)}", ex);
                                }
                            }

                            index++;
                        }
                    }
                }
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
            foreach (int code in sim.Run(ckt))
            {
                if (code == AC.ExportSmallSignal)
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
                                Assert.That(actual, Is.EqualTo(expected).Within(tol));
                            }
                            catch (Exception ex)
                            {
                                string msg = $"{ex.Message} at {sim.Frequency} Hz";
                                throw new Exception(msg, ex);
                            }
                        }

                        index++;
                    }
                }
            }
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
            foreach (int code in sim.Run(ckt))
            {
                if (code == AC.ExportSmallSignal)
                {
                    using (var exportsIt = exports.GetEnumerator())
                    using (var referencesIt = references.GetEnumerator())
                    {
                        while (exportsIt.MoveNext() && referencesIt.MoveNext())
                        {
                            // Test export
                            var actual = exportsIt.Current?.Value ?? throw new ArgumentNullException();
                            var expected = referencesIt.Current?[index] ?? throw new ArgumentNullException();

                            // Test real part
                            double rtol = Math.Max(Math.Abs(actual.Real), Math.Abs(expected.Real)) * RelTol + AbsTol;
                            double itol = Math.Max(Math.Abs(actual.Imaginary), Math.Abs(expected.Imaginary)) * RelTol +
                                          AbsTol;

                            try
                            {
                                Assert.That(actual.Real, Is.EqualTo(expected.Real).Within(rtol));
                                Assert.That(actual.Imaginary, Is.EqualTo(expected.Imaginary).Within(itol));
                            }
                            catch (Exception ex)
                            {
                                string msg = $"{ex.Message} at {sim.Frequency} Hz";
                                throw new Exception(msg, ex);
                            }
                        }

                        index++;
                    }
                }
            }
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
            foreach (int code in sim.Run(ckt))
            {
                if (code == AC.ExportSmallSignal)
                {
                    using (var exportsIt = exports.GetEnumerator())
                    using (var referencesIt = references.GetEnumerator())
                    {
                        while (exportsIt.MoveNext() && referencesIt.MoveNext())
                        {
                            // Test export
                            var actual = exportsIt.Current?.Value ?? throw new ArgumentNullException();
                            var expected = referencesIt.Current?.Invoke(sim.Frequency) ?? throw new ArgumentNullException();

                            // Test real part
                            double rtol = Math.Max(Math.Abs(actual.Real), Math.Abs(expected.Real)) * RelTol + AbsTol;
                            double itol = Math.Max(Math.Abs(actual.Imaginary), Math.Abs(expected.Imaginary)) * RelTol +
                                          AbsTol;

                            try
                            {
                                Assert.That(actual.Real, Is.EqualTo(expected.Real).Within(rtol));
                                Assert.That(actual.Imaginary, Is.EqualTo(expected.Imaginary).Within(itol));
                            }
                            catch (Exception ex)
                            {
                                string msg = $"{ex.Message} at {sim.Frequency} Hz";
                                throw new Exception(msg, ex);
                            }
                        }
                    }
                }
            }
        }

        /// <summary>
        /// Perform a test for transient analysis
        /// </summary>
        /// <param name="sim">Simulation</param>
        /// <param name="ckt">Circuit</param>
        /// <param name="exports">Exports</param>
        /// <param name="references">References</param>
        protected void AnalyzeTransient(Transient sim, Circuit ckt, IEnumerable<IExport<double>> exports, IEnumerable<double[]> references,
            Action<int> additional = null)
        {
            int index = 0;
            foreach (int code in sim.Run(ckt))
            {
                additional?.Invoke(code);
                if (code == Transient.ExportTransient)
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
                                Assert.That(actual, Is.EqualTo(expected).Within(tol));
                            }
                            catch (Exception ex)
                            {
                                string msg = $"{ex.Message} at t={sim.Time} s";
                                throw new Exception(msg, ex);
                            }
                        }

                        index++;
                    }
                }
            }
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
            foreach (int code in sim.Run(ckt))
            {
                if (code == Transient.ExportTransient)
                {
                    using (var exportsIt = exports.GetEnumerator())
                    using (var referencesIt = references.GetEnumerator())
                    {
                        while (exportsIt.MoveNext() && referencesIt.MoveNext())
                        {
                            double t = sim.Time;
                            double actual = exportsIt.Current?.Value ?? throw new ArgumentNullException();
                            double expected = referencesIt.Current?.Invoke(t) ?? throw new ArgumentNullException();
                            double tol = Math.Max(Math.Abs(actual), Math.Abs(expected)) * RelTol + AbsTol;

                            try
                            {
                                Assert.That(actual, Is.EqualTo(expected).Within(tol));
                            }
                            catch (Exception ex)
                            {
                                string msg = $"{ex.Message} at t={sim.Time} s";
                                throw new Exception(msg, ex);
                            }
                        }
                    }
                }
            }
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
            foreach (int code in sim.Run(ckt))
            {
                if (code == Noise.ExportNoise)
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
                                Assert.That(actual, Is.EqualTo(expected).Within(tol));
                            }
                            catch (Exception ex)
                            {
                                string msg = $"{ex.Message} at {sim.Frequency} Hz";
                                throw new Exception(msg, ex);
                            }
                        }

                        index++;
                    }
                }
            }
        }
    }
}