/* ******************************************************************************
   *  BSIM4 4.8.1 released by Chetan Kumar Dabhi 2/15/2017                      *
   *  BSIM4 Model Equations                                                     *
   ******************************************************************************

   ******************************************************************************
   *  Copyright 2017 Regents of the University of California.                   *
   *  All rights reserved.                                                      *
   *                                                                            *
   *  Project Director: Prof. Chenming Hu.                                      *
   *  Authors: Gary W. Ng, Weidong Liu, Xuemei Xi, Mohan Dunga, Wenwei Yang     *
   *           Ali Niknejad, Shivendra Singh Parihar, Chetan Kumar Dabhi        *
   *           Yogesh Singh Chauhan, Sayeef Salahuddin, Chenming Hu             *
   ******************************************************************************

   ******************************************************************************
   *                          CMC In-Code Statement                             *
   *                                                                            *
   *  The Developer agrees that the following statement will appear in the      *
   *  model code that has been adopted as a CMC Standard.                       *
   *                                                                            *
   *  Software is distributed as is, completely without warranty or service     *
   *  support. The University of California and its employees are not liable    *
   *  for the condition or performance of the software.                         *
   *                                                                            *
   *  The University of California owns the copyright and grants users a        *
   *  perpetual, irrevocable, worldwide, non-exclusive, royalty-free license    *
   *  with respect to the software as set forth below.                          *
   *                                                                            *
   *  The University of California hereby disclaims all implied warranties.     *
   *                                                                            *
   *  The University of California grants the users the right to modify,        *
   *  copy, and redistribute the software and documentation, both within        *
   *  the user's organization and externally, subject to the following          *
   *  restrictions:                                                             *
   *                                                                            *
   *  1. The users agree not to charge for the University of California code    *
   *     itself but may charge for additions, extensions, or support.           *
   *                                                                            *
   *  2. In any product based on the software, the users agree to               *
   *     acknowledge the University of California that developed the            *
   *     software. This acknowledgment shall appear in the product              *
   *     documentation.                                                         *
   *                                                                            *
   *  3. Redistributions to others of source code and documentation must        *
   *     retain the copyright notice, disclaimer, and list of conditions.       *
   *                                                                            *
   *  4. Redistributions to others in binary form must reproduce the            *
   *     copyright notice, disclaimer, and list of conditions in the            *
   *     documentation and/or other materials provided with the                 *
   *     distribution.                                                          *
   *                                                                            *
   *  Agreed to on ______Feb. 15, 2017______________                            *
   *                                                                            *
   *  By: ____University of California, Berkeley___                             *
   *      ____Chenming Hu__________________________                             *
   *      ____Professor in Graduate School ________                             *
   *                                                                            *
   ****************************************************************************** */

using SpiceSharp;
using SpiceSharp.Attributes;
using SpiceSharp.Behaviors;
using SpiceSharp.Components;
using SpiceSharp.Entities;
using SpiceSharp.ParameterSets;
using SpiceSharp.Simulations;
using System;
using System.Collections.Generic;

namespace SpiceSharpBSIM.Components.Semiconductors.BSIM.BSIM4Behaviors
{
    /// <summary>
    /// Temperature behavior for a <see cref="BSIM4Model"/>.
    /// </summary>
    [BehaviorFor(typeof(BSIM4Model)), AddBehaviorIfNo(typeof(ITemperatureBehavior))]
    public partial class ModelTemperatureBehavior: Behavior, IParameterized<ModelParameters>, ITemperatureBehavior
    {
        private readonly ITemperatureSimulationState _temperature;

        public const double EPS0 = 8.85418e-12;
        public const double EPSSI = 1.03594e-10;
        public const double KboQ = 8.617087e-5;

        /// <summary>
        /// Gets the size-dependent properties.
        /// </summary>
        public Dictionary<Tuple<double, double, double>, SizeDependentProperties> SizeDependentProperties = new Dictionary<Tuple<double, double, double>, SizeDependentProperties>();

        /// <inheritdoc />
        public ModelParameters Parameters { get; }
        public double DMCGeff { get; private set; }
        public double DMCIeff { get; private set; }
        public double DMDGeff { get; private set; }
        public double Coxe { get; private set; }
        public double Coxp { get; private set; }
        public double Vcrit { get; private set; }
        public double Factor1 { get; private set; }
        public double Vtm0 { get; private set; }
        public double Eg0 { get; private set; }
        public double Vtm { get; private set; }
        public double SjctTempSatCurDensity { get; private set; }
        public double SjctSidewallTempSatCurDensity { get; private set; }
        public double SjctGateSidewallTempSatCurDensity { get; private set; }
        public double DjctTempSatCurDensity { get; private set; }
        public double DjctSidewallTempSatCurDensity { get; private set; }
        public double DjctGateSidewallTempSatCurDensity { get; private set; }
        public double SunitAreaTempJctCap { get; private set; }
        public double DunitAreaTempJctCap { get; private set; }
        public double SunitLengthSidewallTempJctCap { get; private set; }
        public double DunitLengthSidewallTempJctCap { get; private set; }
        public double SunitLengthGateSidewallTempJctCap { get; private set; }
        public double DunitLengthGateSidewallTempJctCap { get; private set; }
        public double PhiBS { get; private set; }
        public double PhiBD { get; private set; }
        public double PhiBSWS { get; private set; }
        public double PhiBSWD { get; private set; }
        public double PhiBSWGS { get; private set; }
        public double PhiBSWGD { get; private set; }
        public double TRatio { get; private set; }
        public double DelTemp { get; private set; }
        public double Ni { get; private set; }
        public double Epssub { get; private set; }
        public double Epsrox { get; private set; }
        public double Toxe { get; private set; }
        public double Njtsstemp { get; private set; }
        public double Njtsswstemp { get; private set; }
        public double Njtsswgstemp { get; private set; }
        public double Njtsdtemp { get; private set; }
        public double Njtsswdtemp { get; private set; }
        public double Njtsswgdtemp { get; private set; }

        /// <summary>
        /// Creates a new <see cref="ModelTemperatureBehavior"/>.
        /// </summary>
        /// <param name="context">The context.</param>
        public ModelTemperatureBehavior(BindingContext context)
            : base(context)
        {
            Parameters = context.GetParameterSet<ModelParameters>();
            _temperature = context.GetState<ITemperatureSimulationState>();
            Setup();
        }

        private void Setup()
        {
            if (!Parameters.MobMod.Given)
                Parameters.MobMod = new GivenParameter<int>(0, false);
            else if ((Parameters.MobMod != 0) && (Parameters.MobMod != 1) && (Parameters.MobMod != 2) && (Parameters.MobMod != 3)
                     && (Parameters.MobMod != 4) && (Parameters.MobMod != 5) && (Parameters.MobMod != 6))
            {
                Parameters.MobMod = new GivenParameter<int>(0, false);
                SpiceSharpWarning.Warning(this, "Warning: mobMod has been set to its default value: 0.");
            }
            if (!Parameters.DioMod.Given)
                Parameters.DioMod = new GivenParameter<int>(1, false);
            else if ((Parameters.DioMod != 0) && (Parameters.DioMod != 1)
                && (Parameters.DioMod != 2))
            {
                Parameters.DioMod = new GivenParameter<int>(1, false);
                SpiceSharpWarning.Warning(this, "Warning: dioMod has been set to its default value: 1.");
            }
            if (!Parameters.CapMod.Given)
                Parameters.CapMod = new GivenParameter<int>(2, false);
            else if ((Parameters.CapMod != 0) && (Parameters.CapMod != 1)
                && (Parameters.CapMod != 2))
            {
                Parameters.CapMod = new GivenParameter<int>(2, false);
                SpiceSharpWarning.Warning(this, "Warning: capMod has been set to its default value: 2.");
            }
            if (!Parameters.RdsMod.Given)
                Parameters.RdsMod = new GivenParameter<int>(0, false);
            else if ((Parameters.RdsMod != 0) && (Parameters.RdsMod != 1))
            {
                Parameters.RdsMod = new GivenParameter<int>(0, false);
                SpiceSharpWarning.Warning(this, "Warning: rdsMod has been set to its default value: 0.");
            }
            if (!Parameters.RbodyMod.Given)
                Parameters.RbodyMod = new GivenParameter<int>(0, false);
            else if ((Parameters.RbodyMod != 0) && (Parameters.RbodyMod != 1) && (Parameters.RbodyMod != 2))
            {
                Parameters.RbodyMod = new GivenParameter<int>(0, false);
                SpiceSharpWarning.Warning(this, "Warning: rbodyMod has been set to its default value: 0.");
            }
            if (!Parameters.RgateMod.Given)
                Parameters.RgateMod = new GivenParameter<int>(0, false);
            else if ((Parameters.RgateMod != 0) && (Parameters.RgateMod != 1)
                && (Parameters.RgateMod != 2) && (Parameters.RgateMod != 3))
            {
                Parameters.RgateMod = new GivenParameter<int>(0, false);
                SpiceSharpWarning.Warning(this, "Warning: rgateMod has been set to its default value: 0.");
            }
            if (!Parameters.PerMod.Given)
                Parameters.PerMod = new GivenParameter<int>(1, false);
            else if ((Parameters.PerMod != 0) && (Parameters.PerMod != 1))
            {
                Parameters.PerMod = new GivenParameter<int>(1, false);
                SpiceSharpWarning.Warning(this, "Warning: perMod has been set to its default value: 1.");
            }
            if (!Parameters.RgeoMod.Given)
                Parameters.RgeoMod = new GivenParameter<int>(0, false);
            else if ((Parameters.RgeoMod != 0) && (Parameters.RgeoMod != 1))
            {
                Parameters.RgeoMod = new GivenParameter<int>(1, false);
                SpiceSharpWarning.Warning(this, "Warning: rgeoMod has been set to its default value: 1.");
            }
            if (!Parameters.FnoiMod.Given)
                Parameters.FnoiMod = new GivenParameter<int>(1, false);
            else if ((Parameters.FnoiMod != 0) && (Parameters.FnoiMod != 1))
            {
                Parameters.FnoiMod = new GivenParameter<int>(1, false);
                SpiceSharpWarning.Warning(this, "Warning: fnoiMod has been set to its default value: 1.");
            }
            if (!Parameters.TnoiMod.Given)
                Parameters.TnoiMod = new GivenParameter<int>(0, false);
            else if ((Parameters.TnoiMod != 0) && (Parameters.TnoiMod != 1) && (Parameters.TnoiMod != 2))
            {
                Parameters.TnoiMod = new GivenParameter<int>(0, false);
                SpiceSharpWarning.Warning(this, "Warning: tnoiMod has been set to its default value: 0.");
            }
            if (!Parameters.TrnqsMod.Given)
                Parameters.TrnqsMod = new GivenParameter<int>(0, false);
            else if ((Parameters.TrnqsMod != 0) && (Parameters.TrnqsMod != 1))
            {
                Parameters.TrnqsMod = new GivenParameter<int>(0, false);
                SpiceSharpWarning.Warning(this, "Warning: trnqsMod has been set to its default value: 0.");
            }
            if (!Parameters.AcnqsMod.Given)
                Parameters.AcnqsMod = new GivenParameter<int>(0, false);
            else if ((Parameters.AcnqsMod != 0) && (Parameters.AcnqsMod != 1))
            {
                Parameters.AcnqsMod = new GivenParameter<int>(0, false);
                SpiceSharpWarning.Warning(this, "Warning: acnqsMod has been set to its default value: 0.");
            }
            if (!Parameters.MtrlMod.Given)
                Parameters.MtrlMod = new GivenParameter<int>(0, false);
            else if ((Parameters.MtrlMod != 0) && (Parameters.MtrlMod != 1))
            {
                Parameters.MtrlMod = new GivenParameter<int>(0, false);
                SpiceSharpWarning.Warning(this, "Warning: mtrlMod has been set to its default value: 0.");
            }
            if (!Parameters.MtrlCompatMod.Given)
                Parameters.MtrlCompatMod = new GivenParameter<int>(0, false);
            else if ((Parameters.MtrlCompatMod != 0) && (Parameters.MtrlCompatMod != 1))
            {
                Parameters.MtrlCompatMod = new GivenParameter<int>(0, false);
                SpiceSharpWarning.Warning(this, "Warning: mtrlCompatMod has been set to its default value: 0.");
            }
            if (!Parameters.IgcMod.Given)
                Parameters.IgcMod = new GivenParameter<int>(0, false);
            else if ((Parameters.IgcMod != 0) && (Parameters.IgcMod != 1)
                      && (Parameters.IgcMod != 2))
            {
                Parameters.IgcMod = new GivenParameter<int>(0, false);
                SpiceSharpWarning.Warning(this, "Warning: igcMod has been set to its default value: 0.");
            }
            if (!Parameters.IgbMod.Given)
                Parameters.IgbMod = new GivenParameter<int>(0, false);
            else if ((Parameters.IgbMod != 0) && (Parameters.IgbMod != 1))
            {
                Parameters.IgbMod = new GivenParameter<int>(0, false);
                SpiceSharpWarning.Warning(this, "Warning: igbMod has been set to its default value: 0.");
            }
            if (!Parameters.TempMod.Given)
                Parameters.TempMod = new GivenParameter<int>(0, false);
            else if ((Parameters.TempMod != 0) && (Parameters.TempMod != 1)
                      && (Parameters.TempMod != 2) && (Parameters.TempMod != 3))
            {
                Parameters.TempMod = new GivenParameter<int>(0, false);
                SpiceSharpWarning.Warning(this, "Warning: tempMod has been set to its default value: 0.");
            }
            if (!Parameters.Version.Given)
                Parameters.Version = new GivenParameter<string>("4.8.1", false);
            if (!Parameters.Vddeot.Given)
                Parameters.Vddeot = new GivenParameter<double>((Parameters.Type > 0) ? 1.5 : -1.5, false);
            if (!Parameters.Toxp.Given)
                Parameters.Toxp = new GivenParameter<double>(Parameters.Toxe, false);
            if (!Parameters.Toxm.Given)
                Parameters.Toxm = new GivenParameter<double>(Parameters.Toxe, false);
            if (!Parameters.Dsub.Given)
                Parameters.Dsub = new GivenParameter<double>(Parameters.Drout, false);
            if (!Parameters.Vth0.Given)
                Parameters.Vth0 = new GivenParameter<double>((Parameters.Type > 0) ? 0.7 : -0.7, false);
            if (!Parameters.Eu.Given)
                Parameters.Eu = new GivenParameter<double>((Parameters.Type > 0) ? 1.67 : 1.0, false);
            if (!Parameters.Ucs.Given)
                Parameters.Ucs = new GivenParameter<double>((Parameters.Type > 0) ? 1.67 : 1.0, false);
            if (Parameters.Version.Value != "4.8.1" && Parameters.Version.Value != "4.81")
            {
                if (!Parameters.Ua.Given)
                    Parameters.Ua = new GivenParameter<double>(((Parameters.MobMod == 2)) ? 1.0e-15 : 1.0e-9, false);
                if (!Parameters.Uc.Given)
                    Parameters.Uc = new GivenParameter<double>((Parameters.MobMod == 1) ? -0.0465 : -0.0465e-9, false);
                if (!Parameters.Uc1.Given)
                    Parameters.Uc1 = new GivenParameter<double>((Parameters.MobMod == 1) ? -0.056 : -0.056e-9, false);
            }
            else
            {
                if (!Parameters.Ua.Given)
                    Parameters.Ua = new GivenParameter<double>(((Parameters.MobMod == 2 || Parameters.MobMod == 6)) ? 1.0e-15 : 1.0e-9, false);
                if (!Parameters.Uc.Given)
                    Parameters.Uc = new GivenParameter<double>((Parameters.MobMod == 1 || Parameters.MobMod == 5) ? -0.0465 : -0.0465e-9, false);
                if (!Parameters.Uc1.Given)
                    Parameters.Uc1 = new GivenParameter<double>((Parameters.MobMod == 1 || Parameters.MobMod == 5) ? -0.056 : -0.056e-9, false);
            }
            if (!Parameters.U0.Given)
                Parameters.U0 = new GivenParameter<double>((Parameters.Type > 0) ? 0.067 : 0.025, false);
            if (!Parameters.Agisl.Given)
                Parameters.Agisl = new GivenParameter<double>(Parameters.Agidl, false);
            if (!Parameters.Bgisl.Given)
                Parameters.Bgisl = new GivenParameter<double>(Parameters.Bgidl, false);
            if (!Parameters.Cgisl.Given)
                Parameters.Cgisl = new GivenParameter<double>(Parameters.Cgidl, false);
            if (!Parameters.Egisl.Given)
                Parameters.Egisl = new GivenParameter<double>(Parameters.Egidl, false);
            if (!Parameters.Rgisl.Given)
                Parameters.Rgisl = new GivenParameter<double>(Parameters.Rgidl, false);
            if (!Parameters.Kgisl.Given)
                Parameters.Kgisl = new GivenParameter<double>(Parameters.Kgidl, false);
            if (!Parameters.Fgisl.Given)
                Parameters.Fgisl = new GivenParameter<double>(Parameters.Fgidl, false);
            if (!Parameters.Aigc.Given)
                Parameters.Aigc = new GivenParameter<double>((Parameters.Type > 0) ? 1.36e-2 : 9.80e-3, false);
            if (!Parameters.Bigc.Given)
                Parameters.Bigc = new GivenParameter<double>((Parameters.Type > 0) ? 1.71e-3 : 7.59e-4, false);
            if (!Parameters.Cigc.Given)
                Parameters.Cigc = new GivenParameter<double>((Parameters.Type > 0) ? 0.075 : 0.03, false);
            if (Parameters.Aigsd.Given)
            {
                Parameters.Aigs = new GivenParameter<double>(Parameters.Aigd = Parameters.Aigsd, false);
            }
            else
            {
                Parameters.Aigsd = new GivenParameter<double>((Parameters.Type > 0) ? 1.36e-2 : 9.80e-3, false);
                if (!Parameters.Aigs.Given)
                    Parameters.Aigs = new GivenParameter<double>((Parameters.Type > 0) ? 1.36e-2 : 9.80e-3, false);
                if (!Parameters.Aigd.Given)
                    Parameters.Aigd = new GivenParameter<double>((Parameters.Type > 0) ? 1.36e-2 : 9.80e-3, false);
            }
            if (Parameters.Bigsd.Given)
            {
                Parameters.Bigs = new GivenParameter<double>(Parameters.Bigd = Parameters.Bigsd, false);
            }
            else
            {
                Parameters.Bigsd = new GivenParameter<double>((Parameters.Type > 0) ? 1.71e-3 : 7.59e-4, false);
                if (!Parameters.Bigs.Given)
                    Parameters.Bigs = new GivenParameter<double>((Parameters.Type > 0) ? 1.71e-3 : 7.59e-4, false);
                if (!Parameters.Bigd.Given)
                    Parameters.Bigd = new GivenParameter<double>((Parameters.Type > 0) ? 1.71e-3 : 7.59e-4, false);
            }
            if (Parameters.Cigsd.Given)
            {
                Parameters.Cigs = new GivenParameter<double>(Parameters.Cigd = Parameters.Cigsd, false);
            }
            else
            {
                Parameters.Cigsd = new GivenParameter<double>((Parameters.Type > 0) ? 0.075 : 0.03, false);
                if (!Parameters.Cigs.Given)
                    Parameters.Cigs = new GivenParameter<double>((Parameters.Type > 0) ? 0.075 : 0.03, false);
                if (!Parameters.Cigd.Given)
                    Parameters.Cigd = new GivenParameter<double>((Parameters.Type > 0) ? 0.075 : 0.03, false);
            }
            if (!Parameters.Ijthdfwd.Given)
                Parameters.Ijthdfwd = new GivenParameter<double>(Parameters.Ijthsfwd, false);
            if (!Parameters.Ijthdrev.Given)
                Parameters.Ijthdrev = new GivenParameter<double>(Parameters.Ijthsrev, false);
            if (!Parameters.Xjbvd.Given)
                Parameters.Xjbvd = new GivenParameter<double>(Parameters.Xjbvs, false);
            if (!Parameters.Bvd.Given)
                Parameters.Bvd = new GivenParameter<double>(Parameters.Bvs, false);
            if (!Parameters.Ckappad.Given)
                Parameters.Ckappad = new GivenParameter<double>(Parameters.Ckappas, false);
            if (!Parameters.Dmci.Given)
                Parameters.Dmci = new GivenParameter<double>(Parameters.Dmcg, false);
            if (!Parameters.Lagisl.Given)
                Parameters.Lagisl = new GivenParameter<double>(Parameters.Lagidl, false);
            if (!Parameters.Lbgisl.Given)
                Parameters.Lbgisl = new GivenParameter<double>(Parameters.Lbgidl, false);
            if (!Parameters.Lcgisl.Given)
                Parameters.Lcgisl = new GivenParameter<double>(Parameters.Lcgidl, false);
            if (!Parameters.Legisl.Given)
                Parameters.Legisl = new GivenParameter<double>(Parameters.Legidl, false);
            if (!Parameters.Lrgisl.Given)
                Parameters.Lrgisl = new GivenParameter<double>(Parameters.Lrgidl, false);
            if (!Parameters.Lkgisl.Given)
                Parameters.Lkgisl = new GivenParameter<double>(Parameters.Lkgidl, false);
            if (!Parameters.Lfgisl.Given)
                Parameters.Lfgisl = new GivenParameter<double>(Parameters.Lfgidl, false);
            if (!Parameters.Aigsd.Given && (Parameters.Aigs.Given || Parameters.Aigd.Given))
            {
                if (!Parameters.Laigs.Given)
                    Parameters.Laigs = new GivenParameter<double>(0.0, false);
            }
            else
            {
                if (!Parameters.Laigsd.Given)
                    Parameters.Laigsd = new GivenParameter<double>(0.0, false);
                Parameters.Laigs = new GivenParameter<double>(Parameters.Laigd = Parameters.Laigsd, false);
            }
            if (!Parameters.Bigsd.Given && (Parameters.Bigs.Given || Parameters.Bigd.Given))
            {
                if (!Parameters.Lbigs.Given)
                    Parameters.Lbigs = new GivenParameter<double>(0.0, false);
            }
            else
            {
                if (!Parameters.Lbigsd.Given)
                    Parameters.Lbigsd = new GivenParameter<double>(0.0, false);
                Parameters.Lbigs = new GivenParameter<double>(Parameters.Lbigd = Parameters.Lbigsd, false);
            }
            if (!Parameters.Cigsd.Given && (Parameters.Cigs.Given || Parameters.Cigd.Given))
            {
                if (!Parameters.Lcigs.Given)
                    Parameters.Lcigs = new GivenParameter<double>(0.0, false);
            }
            else
            {
                if (!Parameters.Lcigsd.Given)
                    Parameters.Lcigsd = new GivenParameter<double>(0.0, false);
                Parameters.Lcigs = new GivenParameter<double>(Parameters.Lcigd = Parameters.Lcigsd, false);
            }
            if (!Parameters.Wagisl.Given)
                Parameters.Wagisl = new GivenParameter<double>(Parameters.Wagidl, false);
            if (!Parameters.Wbgisl.Given)
                Parameters.Wbgisl = new GivenParameter<double>(Parameters.Wbgidl, false);
            if (!Parameters.Wcgisl.Given)
                Parameters.Wcgisl = new GivenParameter<double>(Parameters.Wcgidl, false);
            if (!Parameters.Wegisl.Given)
                Parameters.Wegisl = new GivenParameter<double>(Parameters.Wegidl, false);
            if (!Parameters.Wrgisl.Given)
                Parameters.Wrgisl = new GivenParameter<double>(Parameters.Wrgidl, false);
            if (!Parameters.Wkgisl.Given)
                Parameters.Wkgisl = new GivenParameter<double>(Parameters.Wkgidl, false);
            if (!Parameters.Wfgisl.Given)
                Parameters.Wfgisl = new GivenParameter<double>(Parameters.Wfgidl, false);
            if (!Parameters.Aigsd.Given && (Parameters.Aigs.Given || Parameters.Aigd.Given))
            {
                if (!Parameters.Waigs.Given)
                    Parameters.Waigs = new GivenParameter<double>(0.0, false);
            }
            else
            {
                if (!Parameters.Waigsd.Given)
                    Parameters.Waigsd = new GivenParameter<double>(0.0, false);
                Parameters.Waigs = new GivenParameter<double>(Parameters.Waigd = Parameters.Waigsd, false);
            }
            if (!Parameters.Bigsd.Given && (Parameters.Bigs.Given || Parameters.Bigd.Given))
            {
                if (!Parameters.Wbigs.Given)
                    Parameters.Wbigs = new GivenParameter<double>(0.0, false);
            }
            else
            {
                if (!Parameters.Wbigsd.Given)
                    Parameters.Wbigsd = new GivenParameter<double>(0.0, false);
                Parameters.Wbigs = new GivenParameter<double>(Parameters.Wbigd = Parameters.Wbigsd, false);
            }
            if (!Parameters.Cigsd.Given && (Parameters.Cigs.Given || Parameters.Cigd.Given))
            {
                if (!Parameters.Wcigs.Given)
                    Parameters.Wcigs = new GivenParameter<double>(0.0, false);
            }
            else
            {
                if (!Parameters.Wcigsd.Given)
                    Parameters.Wcigsd = new GivenParameter<double>(0.0, false);
                Parameters.Wcigs = new GivenParameter<double>(Parameters.Wcigd = Parameters.Wcigsd, false);
            }
            if (!Parameters.Pagisl.Given)
                Parameters.Pagisl = new GivenParameter<double>(Parameters.Pagidl, false);
            if (!Parameters.Pbgisl.Given)
                Parameters.Pbgisl = new GivenParameter<double>(Parameters.Pbgidl, false);
            if (!Parameters.Pcgisl.Given)
                Parameters.Pcgisl = new GivenParameter<double>(Parameters.Pcgidl, false);
            if (!Parameters.Pegisl.Given)
                Parameters.Pegisl = new GivenParameter<double>(Parameters.Pegidl, false);
            if (!Parameters.Prgisl.Given)
                Parameters.Prgisl = new GivenParameter<double>(Parameters.Prgidl, false);
            if (!Parameters.Pkgisl.Given)
                Parameters.Pkgisl = new GivenParameter<double>(Parameters.Pkgidl, false);
            if (!Parameters.Pfgisl.Given)
                Parameters.Pfgisl = new GivenParameter<double>(Parameters.Pfgidl, false);
            if (!Parameters.Aigsd.Given && (Parameters.Aigs.Given || Parameters.Aigd.Given))
            {
                if (!Parameters.Paigs.Given)
                    Parameters.Paigs = new GivenParameter<double>(0.0, false);
            }
            else
            {
                if (!Parameters.Paigsd.Given)
                    Parameters.Paigsd = new GivenParameter<double>(0.0, false);
                Parameters.Paigs = new GivenParameter<double>(Parameters.Paigd = Parameters.Paigsd, false);
            }
            if (!Parameters.Bigsd.Given && (Parameters.Bigs.Given || Parameters.Bigd.Given))
            {
                if (!Parameters.Pbigs.Given)
                    Parameters.Pbigs = new GivenParameter<double>(0.0, false);
            }
            else
            {
                if (!Parameters.Pbigsd.Given)
                    Parameters.Pbigsd = new GivenParameter<double>(0.0, false);
                Parameters.Pbigs = new GivenParameter<double>(Parameters.Pbigd = Parameters.Pbigsd, false);
            }
            if (!Parameters.Cigsd.Given && (Parameters.Cigs.Given || Parameters.Cigd.Given))
            {
                if (!Parameters.Pcigs.Given)
                    Parameters.Pcigs = new GivenParameter<double>(0.0, false);
            }
            else
            {
                if (!Parameters.Pcigsd.Given)
                    Parameters.Pcigsd = new GivenParameter<double>(0.0, false);
                Parameters.Pcigs = new GivenParameter<double>(Parameters.Pcigd = Parameters.Pcigsd, false);
            }
            if (!Parameters.Tnom.Given)
                Parameters.Tnom = new GivenParameter<double>(_temperature.NominalTemperature, false);
            if (!Parameters.Llc.Given)
                Parameters.Llc = new GivenParameter<double>(Parameters.Ll, false);
            if (!Parameters.Lwc.Given)
                Parameters.Lwc = new GivenParameter<double>(Parameters.Lw, false);
            if (!Parameters.Lwlc.Given)
                Parameters.Lwlc = new GivenParameter<double>(Parameters.Lwl, false);
            if (!Parameters.Wlc.Given)
                Parameters.Wlc = new GivenParameter<double>(Parameters.Wl, false);
            if (!Parameters.Wwc.Given)
                Parameters.Wwc = new GivenParameter<double>(Parameters.Ww, false);
            if (!Parameters.Wwlc.Given)
                Parameters.Wwlc = new GivenParameter<double>(Parameters.Wwl, false);
            if (!Parameters.Dwc.Given)
                Parameters.Dwc = new GivenParameter<double>(Parameters.Wint, false);
            if (!Parameters.Dlc.Given)
                Parameters.Dlc = new GivenParameter<double>(Parameters.Lint, false);
            if (!Parameters.Dlcig.Given)
                Parameters.Dlcig = new GivenParameter<double>(Parameters.Lint, false);
            if (!Parameters.Dlcigd.Given)
            {
                if (Parameters.Dlcig.Given)
                    Parameters.Dlcigd = new GivenParameter<double>(Parameters.Dlcig, false);
                else
                    Parameters.Dlcigd = new GivenParameter<double>(Parameters.Lint, false);
            }
            if (!Parameters.Dwj.Given)
                Parameters.Dwj = new GivenParameter<double>(Parameters.Dwc, false);
            if (!Parameters.Cf.Given)
                Parameters.Cf = new GivenParameter<double>(2.0 * Parameters.Epsrox * EPS0 / Math.PI
                               * Math.Log(1.0 + 0.4e-6 / Parameters.Toxe), false);
            if (!Parameters.DunitAreaJctCap.Given)
                Parameters.DunitAreaJctCap = new GivenParameter<double>(Parameters.SunitAreaJctCap, false);
            if (!Parameters.DunitLengthSidewallJctCap.Given)
                Parameters.DunitLengthSidewallJctCap = new GivenParameter<double>(Parameters.SunitLengthSidewallJctCap, false);
            if (!Parameters.SunitLengthGateSidewallJctCap.Given)
                Parameters.SunitLengthGateSidewallJctCap = new GivenParameter<double>(Parameters.SunitLengthSidewallJctCap, false);
            if (!Parameters.DunitLengthGateSidewallJctCap.Given)
                Parameters.DunitLengthGateSidewallJctCap = new GivenParameter<double>(Parameters.SunitLengthGateSidewallJctCap, false);
            if (!Parameters.DjctSatCurDensity.Given)
                Parameters.DjctSatCurDensity = new GivenParameter<double>(Parameters.SjctSatCurDensity, false);
            if (!Parameters.DjctSidewallSatCurDensity.Given)
                Parameters.DjctSidewallSatCurDensity = new GivenParameter<double>(Parameters.SjctSidewallSatCurDensity, false);
            if (!Parameters.DjctGateSidewallSatCurDensity.Given)
                Parameters.DjctGateSidewallSatCurDensity = new GivenParameter<double>(Parameters.SjctGateSidewallSatCurDensity, false);
            if (!Parameters.DbulkJctPotential.Given)
                Parameters.DbulkJctPotential = new GivenParameter<double>(Parameters.SbulkJctPotential, false);
            if (!Parameters.DsidewallJctPotential.Given)
                Parameters.DsidewallJctPotential = new GivenParameter<double>(Parameters.SsidewallJctPotential, false);
            if (!Parameters.SGatesidewallJctPotential.Given)
                Parameters.SGatesidewallJctPotential = new GivenParameter<double>(Parameters.SsidewallJctPotential, false);
            if (!Parameters.DGatesidewallJctPotential.Given)
                Parameters.DGatesidewallJctPotential = new GivenParameter<double>(Parameters.SGatesidewallJctPotential, false);
            if (!Parameters.DbulkJctBotGradingCoeff.Given)
                Parameters.DbulkJctBotGradingCoeff = new GivenParameter<double>(Parameters.SbulkJctBotGradingCoeff, false);
            if (!Parameters.DbulkJctSideGradingCoeff.Given)
                Parameters.DbulkJctSideGradingCoeff = new GivenParameter<double>(Parameters.SbulkJctSideGradingCoeff, false);
            if (!Parameters.SbulkJctGateSideGradingCoeff.Given)
                Parameters.SbulkJctGateSideGradingCoeff = new GivenParameter<double>(Parameters.SbulkJctSideGradingCoeff, false);
            if (!Parameters.DbulkJctGateSideGradingCoeff.Given)
                Parameters.DbulkJctGateSideGradingCoeff = new GivenParameter<double>(Parameters.SbulkJctGateSideGradingCoeff, false);
            if (!Parameters.DjctEmissionCoeff.Given)
                Parameters.DjctEmissionCoeff = new GivenParameter<double>(Parameters.SjctEmissionCoeff, false);
            if (!Parameters.DjctTempExponent.Given)
                Parameters.DjctTempExponent = new GivenParameter<double>(Parameters.SjctTempExponent, false);
            if (!Parameters.Jtsd.Given)
                Parameters.Jtsd = new GivenParameter<double>(Parameters.Jtss, false);
            if (!Parameters.Jtsswd.Given)
                Parameters.Jtsswd = new GivenParameter<double>(Parameters.Jtssws, false);
            if (!Parameters.Jtsswgd.Given)
                Parameters.Jtsswgd = new GivenParameter<double>(Parameters.Jtsswgs, false);
            if (!Parameters.Njtsd.Given)
            {
                if (Parameters.Njts.Given)
                    Parameters.Njtsd = new GivenParameter<double>(Parameters.Njts, false);
                else
                    Parameters.Njtsd = new GivenParameter<double>(20.0, false);
            }
            if (!Parameters.Njtsswd.Given)
            {
                if (Parameters.Njtssw.Given)
                    Parameters.Njtsswd = new GivenParameter<double>(Parameters.Njtssw, false);
                else
                    Parameters.Njtsswd = new GivenParameter<double>(20.0, false);
            }
            if (!Parameters.Njtsswgd.Given)
            {
                if (Parameters.Njtsswg.Given)
                    Parameters.Njtsswgd = new GivenParameter<double>(Parameters.Njtsswg, false);
                else
                    Parameters.Njtsswgd = new GivenParameter<double>(20.0, false);
            }
            if (!Parameters.Xtsd.Given)
                Parameters.Xtsd = new GivenParameter<double>(Parameters.Xtss, false);
            if (!Parameters.Xtsswd.Given)
                Parameters.Xtsswd = new GivenParameter<double>(Parameters.Xtssws, false);
            if (!Parameters.Xtsswgd.Given)
                Parameters.Xtsswgd = new GivenParameter<double>(Parameters.Xtsswgs, false);
            if (!Parameters.Tnjtsd.Given)
            {
                if (Parameters.Tnjts.Given)
                    Parameters.Tnjtsd = new GivenParameter<double>(Parameters.Tnjts, false);
                else
                    Parameters.Tnjtsd = new GivenParameter<double>(0.0, false);
            }
            if (!Parameters.Tnjtsswd.Given)
            {
                if (Parameters.Tnjtssw.Given)
                    Parameters.Tnjtsswd = new GivenParameter<double>(Parameters.Tnjtssw, false);
                else
                    Parameters.Tnjtsswd = new GivenParameter<double>(0.0, false);
            }
            if (!Parameters.Tnjtsswgd.Given)
            {
                if (Parameters.Tnjtsswg.Given)
                    Parameters.Tnjtsswgd = new GivenParameter<double>(Parameters.Tnjtsswg, false);
                else
                    Parameters.Tnjtsswgd = new GivenParameter<double>(0.0, false);
            }
            if (!Parameters.Vtsd.Given)
                Parameters.Vtsd = new GivenParameter<double>(Parameters.Vtss, false);
            if (!Parameters.Vtsswd.Given)
                Parameters.Vtsswd = new GivenParameter<double>(Parameters.Vtssws, false);
            if (!Parameters.Vtsswgd.Given)
                Parameters.Vtsswgd = new GivenParameter<double>(Parameters.Vtsswgs, false);
            if (!Parameters.OxideTrapDensityA.Given)
            {
                if (Parameters.Type > 0)
                    Parameters.OxideTrapDensityA = new GivenParameter<double>(6.25e41, false);
                else
                    Parameters.OxideTrapDensityA = new GivenParameter<double>(6.188e40, false);
            }
            if (!Parameters.OxideTrapDensityB.Given)
            {
                if (Parameters.Type > 0)
                    Parameters.OxideTrapDensityB = new GivenParameter<double>(3.125e26, false);
                else
                    Parameters.OxideTrapDensityB = new GivenParameter<double>(1.5e25, false);
            }
            if (!Parameters.Wpemod.Given)
                Parameters.Wpemod = new GivenParameter<double>(0, false);
            else if ((Parameters.Wpemod != 0) && (Parameters.Wpemod != 1))
            {
                Parameters.Wpemod = new GivenParameter<double>(0, false);
                SpiceSharpWarning.Warning(this, "Warning: wpemod has been set to its default value: 0.");
            }
            this.DMCGeff = Parameters.Dmcg - Parameters.Dmcgt;
            this.DMCIeff = Parameters.Dmci;
            this.DMDGeff = Parameters.Dmdg - Parameters.Dmcgt;
        }

        /// <inheritdoc />
        void ITemperatureBehavior.Temperature()
        {
            double Eg, Eg0, ni, epssub;
            double T0, T1, T2, T3;
            double delTemp, Temp, TRatio, Vtm0, Tnom;
            double toxe, epsrox;

            Temp = _temperature.Temperature;
            if (Parameters.SbulkJctPotential < 0.1)
            {
                Parameters.SbulkJctPotential = 0.1;
                SpiceSharpWarning.Warning(this, "Given pbs is less than 0.1. Pbs is set to 0.1.");
            }
            if (Parameters.SsidewallJctPotential < 0.1)
            {
                Parameters.SsidewallJctPotential = 0.1;
                SpiceSharpWarning.Warning(this, "Given pbsws is less than 0.1. Pbsws is set to 0.1.");
            }
            if (Parameters.SGatesidewallJctPotential < 0.1)
            {
                Parameters.SGatesidewallJctPotential = 0.1;
                SpiceSharpWarning.Warning(this, "Given pbswgs is less than 0.1. Pbswgs is set to 0.1.");
            }

            if (Parameters.DbulkJctPotential < 0.1)
            {
                Parameters.DbulkJctPotential = 0.1;
                SpiceSharpWarning.Warning(this, "Given pbd is less than 0.1. Pbd is set to 0.1.");
            }
            if (Parameters.DsidewallJctPotential < 0.1)
            {
                Parameters.DsidewallJctPotential = 0.1;
                SpiceSharpWarning.Warning(this, "Given pbswd is less than 0.1. Pbswd is set to 0.1.");
            }
            if (Parameters.DGatesidewallJctPotential < 0.1)
            {
                Parameters.DGatesidewallJctPotential = 0.1;
                SpiceSharpWarning.Warning(this, "Given pbswgd is less than 0.1. Pbswgd is set to 0.1.");
            }

            if (Parameters.MtrlMod == 0)
            {
                if ((Parameters.Toxe.Given) && (Parameters.Toxp.Given) && (Parameters.Dtox.Given)
                    && (Parameters.Toxe != (Parameters.Toxp + Parameters.Dtox)))
                    SpiceSharpWarning.Warning(this, "Warning: toxe, toxp and dtox all given and toxe != toxp + dtox; dtox ignored.");
                else if ((Parameters.Toxe.Given) && (!Parameters.Toxp.Given))
                    Parameters.Toxp = Parameters.Toxe - Parameters.Dtox;
                else if ((!Parameters.Toxe.Given) && (Parameters.Toxp.Given))
                {
                    Parameters.Toxe = Parameters.Toxp + Parameters.Dtox;
                    if (!Parameters.Toxm.Given)                        /* v4.7 */
                        Parameters.Toxm = Parameters.Toxe;
                }
            }
            else if (Parameters.MtrlCompatMod != 0) /* v4.7 */
            {
                T0 = Parameters.Epsrox / 3.9;
                if ((Parameters.Eot.Given) && (Parameters.Toxp.Given) && (Parameters.Dtox.Given)
                    && (Math.Abs(Parameters.Eot * T0 - (Parameters.Toxp + Parameters.Dtox)) > 1.0e-20))
                {
                    SpiceSharpWarning.Warning(this, "Warning: eot, toxp and dtox all given and eot * EPSROX / 3.9 != toxp + dtox; dtox ignored.");
                }
                else if ((Parameters.Eot.Given) && (!Parameters.Toxp.Given))
                    Parameters.Toxp = T0 * Parameters.Eot - Parameters.Dtox;
                else if ((!Parameters.Eot.Given) && (Parameters.Toxp.Given))
                {
                    Parameters.Eot = (Parameters.Toxp + Parameters.Dtox) / T0;
                    if (!Parameters.Toxm.Given)
                        Parameters.Toxm = Parameters.Eot;
                }
            }

            if (Parameters.MtrlMod.Value != 0)
            {
                epsrox = 3.9;
                toxe = Parameters.Eot;
                epssub = EPS0 * Parameters.Epsrsub;
            }
            else
            {
                epsrox = Parameters.Epsrox;
                toxe = Parameters.Toxe;
                epssub = EPSSI;
            }


            this.Coxe = epsrox * EPS0 / toxe;
            if (Parameters.MtrlMod == 0 || Parameters.MtrlCompatMod != 0)
                this.Coxp = Parameters.Epsrox * EPS0 / Parameters.Toxp;

            if (!Parameters.Cgdo.Given)
            {
                if (Parameters.Dlc.Given && (Parameters.Dlc > 0.0))
                    Parameters.Cgdo = Parameters.Dlc * this.Coxe
                                     - Parameters.Cgdl;
                else
                    Parameters.Cgdo = 0.6 * Parameters.Xj * this.Coxe;
            }
            if (!Parameters.Cgso.Given)
            {
                if (Parameters.Dlc.Given && (Parameters.Dlc > 0.0))
                    Parameters.Cgso = Parameters.Dlc * this.Coxe
                                     - Parameters.Cgsl;
                else
                    Parameters.Cgso = 0.6 * Parameters.Xj * this.Coxe;
            }
            if (!Parameters.Cgbo.Given)
                Parameters.Cgbo = 2.0 * Parameters.Dwc * this.Coxe;

            SizeDependentProperties.Clear();

            Tnom = Parameters.Tnom;
            TRatio = Temp / Tnom;

            this.Vcrit = Constants.Vt0 * Math.Log(Constants.Vt0 / (Constants.Root2 * 1.0e-14));
            this.Factor1 = Math.Sqrt(epssub / (epsrox * EPS0) * toxe);

            Vtm0 = this.Vtm0 = KboQ * Tnom;

            if (Parameters.MtrlMod == 0)
            {
                Eg0 = 1.16 - 7.02e-4 * Tnom * Tnom / (Tnom + 1108.0);
                ni = 1.45e10 * (Tnom / 300.15) * Math.Sqrt(Tnom / 300.15)
                    * Math.Exp(21.5565981 - Eg0 / (2.0 * Vtm0));
            }
            else
            {
                Eg0 = Parameters.Bg0sub - Parameters.Tbgasub * Tnom * Tnom
                                           / (Tnom + Parameters.Tbgbsub);
                T0 = Parameters.Bg0sub - Parameters.Tbgasub * 90090.0225
                                           / (300.15 + Parameters.Tbgbsub);
                ni = Parameters.Ni0sub * (Tnom / 300.15) * Math.Sqrt(Tnom / 300.15)
                      * Math.Exp((T0 - Eg0) / (2.0 * Vtm0));
            }

            this.Eg0 = Eg0;
            this.Vtm = KboQ * Temp;
            if (Parameters.MtrlMod == 0)
                Eg = 1.16 - 7.02e-4 * Temp * Temp / (Temp + 1108.0);
            else
                Eg = Parameters.Bg0sub - Parameters.Tbgasub * Temp * Temp
                                           / (Temp + Parameters.Tbgbsub);
            if (Temp != Tnom)
            {
                T0 = Eg0 / Vtm0 - Eg / this.Vtm;
                T1 = Math.Log(Temp / Tnom);
                T2 = T0 + Parameters.SjctTempExponent * T1;
                T3 = Math.Exp(T2 / Parameters.SjctEmissionCoeff);
                this.SjctTempSatCurDensity = Parameters.SjctSatCurDensity
                                                  * T3;
                this.SjctSidewallTempSatCurDensity
                            = Parameters.SjctSidewallSatCurDensity * T3;
                this.SjctGateSidewallTempSatCurDensity
                            = Parameters.SjctGateSidewallSatCurDensity * T3;

                T2 = T0 + Parameters.DjctTempExponent * T1;
                T3 = Math.Exp(T2 / Parameters.DjctEmissionCoeff);
                this.DjctTempSatCurDensity = Parameters.DjctSatCurDensity
                                                  * T3;
                this.DjctSidewallTempSatCurDensity
                            = Parameters.DjctSidewallSatCurDensity * T3;
                this.DjctGateSidewallTempSatCurDensity
                            = Parameters.DjctGateSidewallSatCurDensity * T3;
            }
            else
            {
                this.SjctTempSatCurDensity = Parameters.SjctSatCurDensity;
                this.SjctSidewallTempSatCurDensity
                           = Parameters.SjctSidewallSatCurDensity;
                this.SjctGateSidewallTempSatCurDensity
                           = Parameters.SjctGateSidewallSatCurDensity;
                this.DjctTempSatCurDensity = Parameters.DjctSatCurDensity;
                this.DjctSidewallTempSatCurDensity
                           = Parameters.DjctSidewallSatCurDensity;
                this.DjctGateSidewallTempSatCurDensity
                           = Parameters.DjctGateSidewallSatCurDensity;
            }

            if (this.SjctTempSatCurDensity < 0.0)
                this.SjctTempSatCurDensity = 0.0;
            if (this.SjctSidewallTempSatCurDensity < 0.0)
                this.SjctSidewallTempSatCurDensity = 0.0;
            if (this.SjctGateSidewallTempSatCurDensity < 0.0)
                this.SjctGateSidewallTempSatCurDensity = 0.0;
            if (this.DjctTempSatCurDensity < 0.0)
                this.DjctTempSatCurDensity = 0.0;
            if (this.DjctSidewallTempSatCurDensity < 0.0)
                this.DjctSidewallTempSatCurDensity = 0.0;
            if (this.DjctGateSidewallTempSatCurDensity < 0.0)
                this.DjctGateSidewallTempSatCurDensity = 0.0;

            /* Temperature dependence of D/B and S/B diode capacitance begins */
            delTemp = _temperature.Temperature - Parameters.Tnom;
            T0 = Parameters.Tcj * delTemp;
            if (T0 >= -1.0)
            {
                this.SunitAreaTempJctCap = Parameters.SunitAreaJctCap * (1.0 + T0); /*bug_fix -JX */
                this.DunitAreaTempJctCap = Parameters.DunitAreaJctCap * (1.0 + T0);
            }
            else
            {
                if (Parameters.SunitAreaJctCap > 0.0)
                {
                    this.SunitAreaTempJctCap = 0.0;
                    SpiceSharpWarning.Warning(this, "Temperature effect has caused cjs to be negative. Cjs is clamped to zero.");
                }
                if (Parameters.DunitAreaJctCap > 0.0)
                {
                    this.DunitAreaTempJctCap = 0.0;
                    SpiceSharpWarning.Warning(this, "Temperature effect has caused cjd to be negative. Cjd is clamped to zero.");
                }
            }
            T0 = Parameters.Tcjsw * delTemp;
            if (Parameters.SunitLengthSidewallJctCap < 0.0)/*4.6.2*/
            {
                Parameters.SunitLengthSidewallJctCap = 0.0;
                SpiceSharpWarning.Warning(this, "CJSWS is negative. Cjsws is clamped to zero.");
            }
            if (Parameters.DunitLengthSidewallJctCap < 0.0)
            {
                Parameters.DunitLengthSidewallJctCap = 0.0;
                SpiceSharpWarning.Warning(this, "CJSWD is negative. Cjswd is clamped to zero.");
            }
            if (T0 >= -1.0)
            {
                this.SunitLengthSidewallTempJctCap = Parameters.SunitLengthSidewallJctCap * (1.0 + T0);
                this.DunitLengthSidewallTempJctCap = Parameters.DunitLengthSidewallJctCap * (1.0 + T0);
            }
            else
            {
                if (Parameters.SunitLengthSidewallJctCap > 0.0)
                {
                    this.SunitLengthSidewallTempJctCap = 0.0;
                    SpiceSharpWarning.Warning(this, "Temperature effect has caused cjsws to be negative. Cjsws is clamped to zero.");
                }
                if (Parameters.DunitLengthSidewallJctCap > 0.0)
                {
                    this.DunitLengthSidewallTempJctCap = 0.0;
                    SpiceSharpWarning.Warning(this, "Temperature effect has caused cjswd to be negative. Cjswd is clamped to zero.");
                }
            }
            T0 = Parameters.Tcjswg * delTemp;
            if (T0 >= -1.0)
            {
                this.SunitLengthGateSidewallTempJctCap = Parameters.SunitLengthGateSidewallJctCap * (1.0 + T0);
                this.DunitLengthGateSidewallTempJctCap = Parameters.DunitLengthGateSidewallJctCap * (1.0 + T0);
            }
            else
            {
                if (Parameters.SunitLengthGateSidewallJctCap > 0.0)
                {
                    this.SunitLengthGateSidewallTempJctCap = 0.0;
                    SpiceSharpWarning.Warning(this, "Temperature effect has caused cjswgs to be negative. Cjswgs is clamped to zero.");
                }
                if (Parameters.DunitLengthGateSidewallJctCap > 0.0)
                {
                    this.DunitLengthGateSidewallTempJctCap = 0.0;
                    SpiceSharpWarning.Warning(this, "Temperature effect has caused cjswgd to be negative. Cjswgd is clamped to zero.");
                }
            }

            this.PhiBS = Parameters.SbulkJctPotential
                              - Parameters.Tpb * delTemp;
            if (this.PhiBS < 0.01)
            {
                this.PhiBS = 0.01;
                SpiceSharpWarning.Warning(this, "Temperature effect has caused pbs to be less than 0.01. Pbs is clamped to 0.01.");
            }
            this.PhiBD = Parameters.DbulkJctPotential
                              - Parameters.Tpb * delTemp;
            if (this.PhiBD < 0.01)
            {
                this.PhiBD = 0.01;
                SpiceSharpWarning.Warning(this, "Temperature effect has caused pbd to be less than 0.01. Pbd is clamped to 0.01.");
            }

            this.PhiBSWS = Parameters.SsidewallJctPotential
                                - Parameters.Tpbsw * delTemp;
            if (this.PhiBSWS <= 0.01)
            {
                this.PhiBSWS = 0.01;
                SpiceSharpWarning.Warning(this, "Temperature effect has caused pbsws to be less than 0.01. Pbsws is clamped to 0.01.");
            }
            this.PhiBSWD = Parameters.DsidewallJctPotential
                                - Parameters.Tpbsw * delTemp;
            if (this.PhiBSWD <= 0.01)
            {
                this.PhiBSWD = 0.01;
                SpiceSharpWarning.Warning(this, "Temperature effect has caused pbswd to be less than 0.01. Pbswd is clamped to 0.01.");
            }

            this.PhiBSWGS = Parameters.SGatesidewallJctPotential
                                 - Parameters.Tpbswg * delTemp;
            if (this.PhiBSWGS <= 0.01)
            {
                this.PhiBSWGS = 0.01;
                SpiceSharpWarning.Warning(this, "Temperature effect has caused pbswgs to be less than 0.01. Pbswgs is clamped to 0.01.");
            }
            this.PhiBSWGD = Parameters.DGatesidewallJctPotential
                                 - Parameters.Tpbswg * delTemp;
            if (this.PhiBSWGD <= 0.01)
            {
                this.PhiBSWGD = 0.01;
                SpiceSharpWarning.Warning(this, "Temperature effect has caused pbswgd to be less than 0.01. Pbswgd is clamped to 0.01.");
            } /* End of junction capacitance */


            if (Parameters.Ijthdfwd <= 0.0)
            {
                Parameters.Ijthdfwd = 0.0;
                SpiceSharpWarning.Warning(this, "Ijthdfwd reset to {0:g}.".FormatString(Parameters.Ijthdfwd));
            }
            if (Parameters.Ijthsfwd <= 0.0)
            {
                Parameters.Ijthsfwd = 0.0;
                SpiceSharpWarning.Warning(this, "Ijthsfwd reset to {0:g}.".FormatString(Parameters.Ijthsfwd));
            }
            if (Parameters.Ijthdrev <= 0.0)
            {
                Parameters.Ijthdrev = 0.0;
                SpiceSharpWarning.Warning(this, "Ijthdrev reset to {0:g}.".FormatString(Parameters.Ijthdrev));
            }
            if (Parameters.Ijthsrev <= 0.0)
            {
                Parameters.Ijthsrev = 0.0;
                SpiceSharpWarning.Warning(this, "Ijthsrev reset to {0:g}.".FormatString(Parameters.Ijthsrev));
            }

            if ((Parameters.Xjbvd <= 0.0) && (Parameters.DioMod == 2))
            {
                Parameters.Xjbvd = 0.0;
                SpiceSharpWarning.Warning(this, "Xjbvd reset to {0:g}.".FormatString(Parameters.Xjbvd));
            }
            else if ((Parameters.Xjbvd < 0.0) && (Parameters.DioMod == 0))
            {
                Parameters.Xjbvd = 0.0;
                SpiceSharpWarning.Warning(this, "Xjbvd reset to {0:g}.".FormatString(Parameters.Xjbvd));
            }

            if (Parameters.Bvd <= 0.0)   /*4.6.2*/
            {
                Parameters.Bvd = 0.0;
                SpiceSharpWarning.Warning(this, "BVD reset to {0:g}.".FormatString(Parameters.Bvd));
            }

            if ((Parameters.Xjbvs <= 0.0) && (Parameters.DioMod == 2))
            {
                Parameters.Xjbvs = 0.0;
                SpiceSharpWarning.Warning(this, "Xjbvs reset to {0:g}.".FormatString(Parameters.Xjbvs));
            }
            else if ((Parameters.Xjbvs < 0.0) && (Parameters.DioMod == 0))
            {
                Parameters.Xjbvs = 0.0;
                SpiceSharpWarning.Warning(this, "Xjbvs reset to {0:g}.".FormatString(Parameters.Xjbvs));
            }

            if (Parameters.Bvs <= 0.0)
            {
                Parameters.Bvs = 0.0;
                SpiceSharpWarning.Warning(this, "BVS reset to {0:g}.".FormatString(Parameters.Bvs));
            }

            // Moved to model temperature - 20220724 - Sven Boulanger
            T0 = (TRatio - 1.0);
            this.Njtsstemp = Parameters.Njts * (1.0 + Parameters.Tnjts * T0);
            this.Njtsswstemp = Parameters.Njtssw * (1.0 + Parameters.Tnjtssw * T0);
            this.Njtsswgstemp = Parameters.Njtsswg * (1.0 + Parameters.Tnjtsswg * T0);
            this.Njtsdtemp = Parameters.Njtsd * (1.0 + Parameters.Tnjtsd * T0);
            this.Njtsswdtemp = Parameters.Njtsswd * (1.0 + Parameters.Tnjtsswd * T0);
            this.Njtsswgdtemp = Parameters.Njtsswgd * (1.0 + Parameters.Tnjtsswgd * T0);

            this.TRatio = TRatio;
            this.DelTemp = delTemp;
            this.Ni = ni;
            this.Epssub = epssub;
            this.Epsrox = epsrox;
            this.Toxe = toxe;
            this.Eg0 = Eg0;
        }
    }
}
