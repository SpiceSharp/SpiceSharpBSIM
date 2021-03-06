using System;
using SpiceSharp.Simulations;
using SpiceSharp.Behaviors;
using SpiceSharp.Algebra;
using SpiceSharp.Components.MosfetBehaviors;
using SpiceSharp.Components.Semiconductors;

namespace SpiceSharp.Components.BSIM3v24Behaviors
{
    /// <summary>
    /// Load behavior for a <see cref="BSIM3v24" />
    /// </summary>
    public class BiasingBehavior : TemperatureBehavior, IBiasingBehavior
    {
        private const double ScalingFactor = 1.0e-9;
        private const double EPSSI = 1.03594e-10;
        private const double EXP_THRESHOLD = 34.0;
        private const double MIN_EXP = 1.713908431e-15;
        private const double MAX_EXP = 5.834617425e14;

        private const double DELTA_1 = 0.02;
        private const double DELTA_3 = 0.02;
        private const double DELTA_4 = 0.02;

        /// <summary>
        /// Necessary behaviors and parameters
        /// </summary>
        protected BaseConfiguration BaseConfiguration { get; private set; }

        /// <summary>
        /// Properties
        /// </summary>
        public double Gbs { get; private set; }
        public double Cbs { get; private set; }
        public double Gbd { get; private set; }
        public double Cbd { get; private set; }
        public double Mode { get; private set; }
        public double Thetavth { get; private set; }
        public double Von { get; private set; }
        public double Vgsteff { get; private set; }
        public double Rds { get; private set; }
        public double Abulk { get; private set; }
        public double Ueff { get; private set; }
        public double AbovVgst2Vtm { get; private set; }
        public double Vdsat { get; private set; }
        public double Vdseff { get; private set; }
        public double Gds { get; private set; }
        public double Gm { get; private set; }
        public double Gmbs { get; private set; }
        public double Gbbs { get; private set; }
        public double Gbgs { get; private set; }
        public double Gbds { get; private set; }
        public double Csub { get; private set; }
        public double Cggb { get; private set; }
        public double Cgsb { get; private set; }
        public double Cgdb { get; private set; }
        public double Cdgb { get; private set; }
        public double Cdsb { get; private set; }
        public double Cddb { get; private set; }
        public double Cbgb { get; private set; }
        public double Cbsb { get; private set; }
        public double Cbdb { get; private set; }
        public double Cqdb { get; private set; }
        public double Cqsb { get; private set; }
        public double Cqgb { get; private set; }
        public double Cqbb { get; private set; }
        public double Gtau { get; private set; }
        public double Qinv { get; private set; }
        public double Qgate { get; private set; }
        public double Qbulk { get; private set; }
        public double Qdrn { get; private set; }
        public double Cd { get; private set; }
        public double Capbs { get; private set; }
        public double Capbd { get; private set; }
        public double Gtg { get; private set; }
        public double Gtd { get; private set; }
        public double Gts { get; private set; }
        public double Gtb { get; private set; }
        public double Vbd { get; set; }
        public double Vbs { get; set; }
        public double Vgs { get; set; }
        public double Vds { get; set; }
        public double Qbs { get; set; }
        public double Qbd { get; set; }
        public double Qdef { get; set; }
        public TransientBehavior TranBehavior { get; set; }

        /// <summary>
        /// Nodes
        /// </summary>
        protected int DrainNode { get; private set; }
        protected int GateNode { get; private set; }
        protected int SourceNode { get; private set; }
        protected int BulkNode { get; private set; }
        protected VectorElement<double> DrainNodePrimePtr { get; private set; }
        public int DrainNodePrime { get; private set; }
        protected VectorElement<double> SourceNodePrimePtr { get; private set; }
        public int SourceNodePrime { get; private set; }
        protected VectorElement<double> QNodePtr { get; private set; }
        public int QNode { get; private set; }
        protected MatrixElement<double> DdPtr { get; private set; }
        protected MatrixElement<double> GgPtr { get; private set; }
        protected MatrixElement<double> SsPtr { get; private set; }
        protected MatrixElement<double> BbPtr { get; private set; }
        protected MatrixElement<double> DPdpPtr { get; private set; }
        protected MatrixElement<double> SPspPtr { get; private set; }
        protected MatrixElement<double> DdpPtr { get; private set; }
        protected MatrixElement<double> GbPtr { get; private set; }
        protected MatrixElement<double> GdpPtr { get; private set; }
        protected MatrixElement<double> GspPtr { get; private set; }
        protected MatrixElement<double> SspPtr { get; private set; }
        protected MatrixElement<double> BdpPtr { get; private set; }
        protected MatrixElement<double> BspPtr { get; private set; }
        protected MatrixElement<double> DPspPtr { get; private set; }
        protected MatrixElement<double> DPdPtr { get; private set; }
        protected MatrixElement<double> BgPtr { get; private set; }
        protected MatrixElement<double> DPgPtr { get; private set; }
        protected MatrixElement<double> SPgPtr { get; private set; }
        protected MatrixElement<double> SPsPtr { get; private set; }
        protected MatrixElement<double> DPbPtr { get; private set; }
        protected MatrixElement<double> SPbPtr { get; private set; }
        protected MatrixElement<double> SPdpPtr { get; private set; }
        protected MatrixElement<double> QqPtr { get; private set; }
        protected MatrixElement<double> QdpPtr { get; private set; }
        protected MatrixElement<double> QspPtr { get; private set; }
        protected MatrixElement<double> QgPtr { get; private set; }
        protected MatrixElement<double> QbPtr { get; private set; }
        protected MatrixElement<double> DPqPtr { get; private set; }
        protected MatrixElement<double> SPqPtr { get; private set; }
        protected MatrixElement<double> GqPtr { get; private set; }
        protected MatrixElement<double> BqPtr { get; private set; }
        protected VectorElement<double> GateNodePtr { get; private set; }
        protected VectorElement<double> BulkNodePtr { get; private set; }

        private BaseSimulationState _state;

        /// <summary>
        /// Constructor
        /// </summary>
        public BiasingBehavior(string name) : base(name)
        {
        }

        /// <summary>
        /// Setup the behavior
        /// </summary>
        public override void Bind(Simulation simulation, BindingContext context)
        {
            base.Bind(simulation, context);

            // Get configuration
            BaseConfiguration = simulation.Configurations.Get<BaseConfiguration>();

            if (context is ComponentBindingContext cc)
            {
                DrainNode = cc.Pins[0];
                GateNode = cc.Pins[1];
                SourceNode = cc.Pins[2];
                BulkNode = cc.Pins[3];
            }

            _state = ((BaseSimulation)simulation).RealState;
            var solver = _state.Solver;
            var variables = simulation.Variables;
            if (ModelParameters.SheetResistance > 0.0 && BaseParameters.DrainSquares > 0.0)
                DrainNodePrime = variables.Create(Name.Combine("drain")).Index;
            else
                DrainNodePrime = DrainNode;

            DrainNodePrimePtr = solver.GetRhsElement(DrainNodePrime);
            if (ModelParameters.SheetResistance > 0.0 && BaseParameters.SourceSquares > 0.0)
                SourceNodePrime = variables.Create(Name.Combine("source")).Index;
            else
                SourceNodePrime = SourceNode;

            SourceNodePrimePtr = solver.GetRhsElement(SourceNodePrime);
            if (BaseParameters.NqsMod > 0)
                QNode = variables.Create(Name.Combine("charge")).Index;
            else
                QNode = 0;

            QNodePtr = solver.GetRhsElement(QNode);
            DdPtr = solver.GetMatrixElement(DrainNode, DrainNode);
            GgPtr = solver.GetMatrixElement(GateNode, GateNode);
            SsPtr = solver.GetMatrixElement(SourceNode, SourceNode);
            BbPtr = solver.GetMatrixElement(BulkNode, BulkNode);
            DPdpPtr = solver.GetMatrixElement(DrainNodePrime, DrainNodePrime);
            SPspPtr = solver.GetMatrixElement(SourceNodePrime, SourceNodePrime);
            DdpPtr = solver.GetMatrixElement(DrainNode, DrainNodePrime);
            GbPtr = solver.GetMatrixElement(GateNode, BulkNode);
            GdpPtr = solver.GetMatrixElement(GateNode, DrainNodePrime);
            GspPtr = solver.GetMatrixElement(GateNode, SourceNodePrime);
            SspPtr = solver.GetMatrixElement(SourceNode, SourceNodePrime);
            BdpPtr = solver.GetMatrixElement(BulkNode, DrainNodePrime);
            BspPtr = solver.GetMatrixElement(BulkNode, SourceNodePrime);
            DPspPtr = solver.GetMatrixElement(DrainNodePrime, SourceNodePrime);
            DPdPtr = solver.GetMatrixElement(DrainNodePrime, DrainNode);
            BgPtr = solver.GetMatrixElement(BulkNode, GateNode);
            DPgPtr = solver.GetMatrixElement(DrainNodePrime, GateNode);
            SPgPtr = solver.GetMatrixElement(SourceNodePrime, GateNode);
            SPsPtr = solver.GetMatrixElement(SourceNodePrime, SourceNode);
            DPbPtr = solver.GetMatrixElement(DrainNodePrime, BulkNode);
            SPbPtr = solver.GetMatrixElement(SourceNodePrime, BulkNode);
            SPdpPtr = solver.GetMatrixElement(SourceNodePrime, DrainNodePrime);
            QqPtr = solver.GetMatrixElement(QNode, QNode);
            QdpPtr = solver.GetMatrixElement(QNode, DrainNodePrime);
            QspPtr = solver.GetMatrixElement(QNode, SourceNodePrime);
            QgPtr = solver.GetMatrixElement(QNode, GateNode);
            QbPtr = solver.GetMatrixElement(QNode, BulkNode);
            DPqPtr = solver.GetMatrixElement(DrainNodePrime, QNode);
            SPqPtr = solver.GetMatrixElement(SourceNodePrime, QNode);
            GqPtr = solver.GetMatrixElement(GateNode, QNode);
            BqPtr = solver.GetMatrixElement(BulkNode, QNode);
            GateNodePtr = solver.GetRhsElement(GateNode);
            BulkNodePtr = solver.GetRhsElement(BulkNode);
        }

        /// <summary>
        /// Load the behavior
        /// </summary>
        void IBiasingBehavior.Load()
        {
            double SourceSatCurrent, DrainSatCurrent;
            double ag0, qgd, qgs, qgb, von, VgstNVt, ExpVgst;
            double cdrain, cdreq, ceqbd, ceqbs, ceqqb, ceqqd, ceqqg;
            double czbd, czbdsw, czbdswg, czbs, czbssw, czbsswg, evbd, evbs, arg, sarg;
            double Vfbeff, dVfbeff_dVg, dVfbeff_dVb, V3, V4;
            double gcbdb, gcbgb, gcbsb, gcddb, gcdgb, gcdsb, gcgdb, gcggb, gcgsb, gcsdb;
            double gcsgb, gcssb, MJ, MJSW, MJSWG;
            double vbd, vbs, vds, vgb, vgd, vgs, vgdo;
            double qgate = 0.0, qbulk = 0.0, qdrn = 0.0, qsrc = 0.0, qinoi, cqgate, cqbulk, cqdrn;
            double Vds, Vgs, Vbs, Gmbs, FwdSum, RevSum;
            double Vgs_eff, Vfb;
            double Phis, dPhis_dVb, sqrtPhis, dsqrtPhis_dVb, Vth, dVth_dVb, dVth_dVd;
            double Vgst, dVgst_dVg, dVgst_dVb, dVgs_eff_dVg, Nvtm;
            double Vtm;
            double n, dn_dVb, dn_dVd, voffcv, noff, dnoff_dVd, dnoff_dVb;
            double ExpArg, V0, CoxWLcen, QovCox, LINK;
            double DeltaPhi, dDeltaPhi_dVg, dDeltaPhi_dVd, dDeltaPhi_dVb;
            double Cox, Tox, Tcen, dTcen_dVg, dTcen_dVd, dTcen_dVb;
            double Ccen, Coxeff, dCoxeff_dVg, dCoxeff_dVd, dCoxeff_dVb;
            double Denomi, dDenomi_dVg, dDenomi_dVd, dDenomi_dVb;
            double ueff, dueff_dVg, dueff_dVd, dueff_dVb;
            double Esat, Vdsat;
            double EsatL, dEsatL_dVg, dEsatL_dVd, dEsatL_dVb;
            double dVdsat_dVg, dVdsat_dVb, dVdsat_dVd, Vasat, dAlphaz_dVg, dAlphaz_dVb;
            double dVasat_dVg, dVasat_dVb, dVasat_dVd, Va, dVa_dVd, dVa_dVg, dVa_dVb;
            double Vbseff, dVbseff_dVb, VbseffCV, dVbseffCV_dVb;
            double Arg1, One_Third_CoxWL, Two_Third_CoxWL, Alphaz, CoxWL;
            double T0, dT0_dVg, dT0_dVd, dT0_dVb;
            double T1, dT1_dVg, dT1_dVd, dT1_dVb;
            double T2, dT2_dVg, dT2_dVd, dT2_dVb;
            double T3, dT3_dVg, dT3_dVd, dT3_dVb;
            double T4;
            double T5;
            double T6;
            double T7;
            double T8;
            double T9;
            double T10;
            double T11, T12;
            double tmp, Abulk, dAbulk_dVb, Abulk0, dAbulk0_dVb;
            double VACLM, dVACLM_dVg, dVACLM_dVd, dVACLM_dVb;
            double VADIBL, dVADIBL_dVg, dVADIBL_dVd, dVADIBL_dVb;
            double Xdep, dXdep_dVb, lt1, dlt1_dVb, ltw, dltw_dVb, Delt_vth, dDelt_vth_dVb;
            double Theta0, dTheta0_dVb;
            double TempRatio, tmp1, tmp2, tmp3, tmp4;
            double DIBL_Sft, dDIBL_Sft_dVd, Lambda, dLambda_dVg;
            double a1;

            double Vgsteff, dVgsteff_dVg, dVgsteff_dVd, dVgsteff_dVb;
            double Vdseff, dVdseff_dVg, dVdseff_dVd, dVdseff_dVb;
            double VdseffCV, dVdseffCV_dVg, dVdseffCV_dVd, dVdseffCV_dVb;
            double diffVds, dAbulk_dVg;
            double beta, dbeta_dVg, dbeta_dVd, dbeta_dVb;
            double gche, dgche_dVg, dgche_dVd, dgche_dVb;
            double fgche1, dfgche1_dVg, dfgche1_dVd, dfgche1_dVb;
            double fgche2, dfgche2_dVg, dfgche2_dVd, dfgche2_dVb;
            double Idl, dIdl_dVg, dIdl_dVd, dIdl_dVb;
            double Idsa, dIdsa_dVg, dIdsa_dVd, dIdsa_dVb;
            double Ids, Gm, Gds, Gmb;
            double Isub, Gbd, Gbg, Gbb;
            double VASCBE, dVASCBE_dVg, dVASCBE_dVd, dVASCBE_dVb;
            double CoxWovL;
            double Rds, dRds_dVg, dRds_dVb, WVCox, WVCoxRds;
            double Vgst2Vtm, VdsatCV, dVdsatCV_dVg, dVdsatCV_dVb;
            double Leff, Weff, dWeff_dVg, dWeff_dVb;
            double AbulkCV, dAbulkCV_dVb;
            double qgdo, qgso, cgdo, cgso;

            double qcheq = 0.0, qdef, gqdef = 0.0, cqdef, cqcheq, gtau_diff, gtau_drift;
            double gcqdb = 0.0, gcqsb = 0.0, gcqgb = 0.0, gcqbb = 0.0;
            double dxpart, sxpart, ggtg, ggtd, ggts, ggtb;
            double ddxpart_dVd, ddxpart_dVg, ddxpart_dVb, ddxpart_dVs;
            double dsxpart_dVd, dsxpart_dVg, dsxpart_dVb, dsxpart_dVs;

            double gbspsp, gbbdp, gbbsp, gbspg, gbspb, gbspdp;
            double gbdpdp, gbdpg, gbdpb, gbdpsp;
            double Cgg, Cgd, Cgb, Cdg, Cdd, Cds;
            double Csg, Csd, Css, Csb, Cbg, Cbd, Cbb;
            double Cgg1, Cgb1, Cgd1, Cbg1, Cbb1, Cbd1, Qac0, Qsub0;
            double dQac0_dVg, dQac0_dVb, dQsub0_dVg, dQsub0_dVd, dQsub0_dVb;

            bool check, chargeComputationNeeded;

            var state = _state;
            chargeComputationNeeded = TranBehavior != null;
            check = true;
            var pParam = base.Param;
            if (Simulation is FrequencySimulation && !state.UseDc)
            {
                vbs = this.Vbs;
                vgs = this.Vgs;
                vds = this.Vds;
                qdef = Qdef;
                chargeComputationNeeded = true;
            }
            else if (state.Init == InitializationModes.Junction && !BaseParameters.Off)
            {
                vds = ModelParameters.B3Type * BaseParameters.IcVDS;
                vgs = ModelParameters.B3Type * BaseParameters.IcVGS;
                vbs = ModelParameters.B3Type * BaseParameters.IcVBS;
                qdef = 0.0;

                if (vds == 0.0 && vgs == 0.0 && vbs == 0.0 &&
                    (state.UseDc || TranBehavior != null || !state.UseIc))
                {
                    vbs = 0.0;
                    vgs = ModelParameters.B3Type * pParam.BSIM3vth0 + 0.1;
                    vds = 0.1;
                }
            }
            else if ((state.Init == InitializationModes.Junction ||
                      state.Init == InitializationModes.Fix) && BaseParameters.Off)
            {
                qdef = vbs = vgs = vds = 0.0;
            }
            else
            {
                vbs = ModelParameters.B3Type * (state.Solution[BulkNode] - state.Solution[SourceNodePrime]);
                vgs = ModelParameters.B3Type * (state.Solution[GateNode] - state.Solution[SourceNodePrime]);
                vds = ModelParameters.B3Type * (state.Solution[DrainNodePrime] - state.Solution[SourceNodePrime]);
                qdef = ModelParameters.B3Type * state.Solution[QNode];

                vbd = vbs - vds;
                vgd = vgs - vds;
                vgdo = this.Vgs - this.Vds;

                von = Von;
                if (this.Vds >= 0.0)
                {
                    vgs = Transistor.LimitFet(vgs, this.Vgs, von);
                    vds = vgs - vgd;
                    vds = Transistor.LimitVds(vds, this.Vds);
                    vgd = vgs - vds;
                }
                else
                {
                    vgd = Transistor.LimitFet(vgd, vgdo, von);
                    vds = vgs - vgd;
                    vds = -Transistor.LimitVds(-vds, -this.Vds);
                    vgs = vgd + vds;
                }

                if (vds >= 0.0)
                {
                    vbs = Semiconductor.LimitJunction(vbs, this.Vbs, Constants.Vt0, ModelTemperature.Vcrit, ref check);
                    vbd = vbs - vds;
                }
                else
                {
                    vbd = Semiconductor.LimitJunction(vbd, Vbd, Constants.Vt0, ModelTemperature.Vcrit, ref check);
                    vbs = vbd + vds;
                }
            }

            /* determine DC current and derivatives */
            vbd = vbs - vds;
            vgd = vgs - vds;
            vgb = vgs - vbs;

            /* Source/drain junction diode DC model begins */
            Nvtm = ModelTemperature.Vtm * ModelParameters.JctEmissionCoeff;
            if (BaseParameters.SourceArea <= 0.0 && BaseParameters.SourcePerimeter <= 0.0)
            {
                SourceSatCurrent = 1.0e-14;
            }
            else
            {
                SourceSatCurrent = BaseParameters.SourceArea * ModelTemperature.JctTempSatCurDensity +
                                   BaseParameters.SourcePerimeter * ModelTemperature.JctSidewallTempSatCurDensity;
            }

            if (SourceSatCurrent <= 0.0)
            {
                Gbs = BaseConfiguration.Gmin;
                Cbs = Gbs * vbs;
            }
            else
            {
                if (ModelParameters.Ijth == 0.0)
                {
                    evbs = Math.Exp(vbs / Nvtm);
                    Gbs = SourceSatCurrent * evbs / Nvtm + BaseConfiguration.Gmin;
                    Cbs = SourceSatCurrent * (evbs - 1.0)
                          + BaseConfiguration.Gmin * vbs;
                }
                else
                {
                    if (vbs < base.Vjsm)
                    {
                        evbs = Math.Exp(vbs / Nvtm);
                        Gbs = SourceSatCurrent * evbs / Nvtm + BaseConfiguration.Gmin;
                        Cbs = SourceSatCurrent * (evbs - 1.0)
                              + BaseConfiguration.Gmin * vbs;
                    }
                    else
                    {
                        T0 = base.IsEvjsm / Nvtm;
                        Gbs = T0 + BaseConfiguration.Gmin;
                        Cbs = base.IsEvjsm - SourceSatCurrent
                              + T0 * (vbs - base.Vjsm)
                              + BaseConfiguration.Gmin * vbs;
                    }
                }
            }

            if (BaseParameters.DrainArea <= 0.0 && BaseParameters.DrainPerimeter <= 0.0)
            {
                DrainSatCurrent = 1.0e-14;
            }
            else
            {
                DrainSatCurrent = BaseParameters.DrainArea
                                  * ModelTemperature.JctTempSatCurDensity
                                  + BaseParameters.DrainPerimeter
                                  * ModelTemperature.JctSidewallTempSatCurDensity;
            }

            if (DrainSatCurrent <= 0.0)
            {
                this.Gbd = BaseConfiguration.Gmin;
                this.Cbd = this.Gbd * vbd;
            }
            else
            {
                if (ModelParameters.Ijth == 0.0)
                {
                    evbd = Math.Exp(vbd / Nvtm);
                    this.Gbd = DrainSatCurrent * evbd / Nvtm + BaseConfiguration.Gmin;
                    this.Cbd = DrainSatCurrent * (evbd - 1.0)
                               + BaseConfiguration.Gmin * vbd;
                }
                else
                {
                    if (vbd < base.Vjdm)
                    {
                        evbd = Math.Exp(vbd / Nvtm);
                        this.Gbd = DrainSatCurrent * evbd / Nvtm + BaseConfiguration.Gmin;
                        this.Cbd = DrainSatCurrent * (evbd - 1.0)
                                   + BaseConfiguration.Gmin * vbd;
                    }
                    else
                    {
                        T0 = base.IsEvjdm / Nvtm;
                        this.Gbd = T0 + BaseConfiguration.Gmin;
                        this.Cbd = base.IsEvjdm - DrainSatCurrent
                                   + T0 * (vbd - base.Vjdm)
                                   + BaseConfiguration.Gmin * vbd;
                    }
                }
            }
            /* End of diode DC model */

            if (vds >= 0.0)
            {
                /* normal mode */
                Mode = 1;
                Vds = vds;
                Vgs = vgs;
                Vbs = vbs;
            }
            else
            {
                /* inverse mode */
                Mode = -1;
                Vds = -vds;
                Vgs = vgd;
                Vbs = vbd;
            }

            T0 = Vbs - pParam.BSIM3vbsc - 0.001;
            T1 = Math.Sqrt(T0 * T0 - 0.004 * pParam.BSIM3vbsc);
            Vbseff = pParam.BSIM3vbsc + 0.5 * (T0 + T1);
            dVbseff_dVb = 0.5 * (1.0 + T0 / T1);
            if (Vbseff < Vbs)
            {
                Vbseff = Vbs;
            }

            if (Vbseff > 0.0)
            {
                T0 = pParam.BSIM3phi / (pParam.BSIM3phi + Vbseff);
                Phis = pParam.BSIM3phi * T0;
                dPhis_dVb = -T0 * T0;
                sqrtPhis = pParam.BSIM3phis3 / (pParam.BSIM3phi + 0.5 * Vbseff);
                dsqrtPhis_dVb = -0.5 * sqrtPhis * sqrtPhis / pParam.BSIM3phis3;
            }
            else
            {
                Phis = pParam.BSIM3phi - Vbseff;
                dPhis_dVb = -1.0;
                sqrtPhis = Math.Sqrt(Phis);
                dsqrtPhis_dVb = -0.5 / sqrtPhis;
            }

            Xdep = pParam.BSIM3Xdep0 * sqrtPhis / pParam.BSIM3sqrtPhi;
            dXdep_dVb = pParam.BSIM3Xdep0 / pParam.BSIM3sqrtPhi
                        * dsqrtPhis_dVb;

            Leff = pParam.BSIM3leff;
            Vtm = ModelTemperature.Vtm;
            /* Vth Calculation */
            T3 = Math.Sqrt(Xdep);
            V0 = pParam.BSIM3vbi - pParam.BSIM3phi;

            T0 = pParam.BSIM3dvt2 * Vbseff;
            if (T0 >= -0.5)
            {
                T1 = 1.0 + T0;
                T2 = pParam.BSIM3dvt2;
            }
            else /* Added to avoid any discontinuity problems caused by dvt2 */
            {
                T4 = 1.0 / (3.0 + 8.0 * T0);
                T1 = (1.0 + 3.0 * T0) * T4;
                T2 = pParam.BSIM3dvt2 * T4 * T4;
            }

            lt1 = ModelTemperature.Factor1 * T3 * T1;
            dlt1_dVb = ModelTemperature.Factor1 * (0.5 / T3 * T1 * dXdep_dVb + T3 * T2);

            T0 = pParam.BSIM3dvt2w * Vbseff;
            if (T0 >= -0.5)
            {
                T1 = 1.0 + T0;
                T2 = pParam.BSIM3dvt2w;
            }
            else /* Added to avoid any discontinuity problems caused by dvt2w */
            {
                T4 = 1.0 / (3.0 + 8.0 * T0);
                T1 = (1.0 + 3.0 * T0) * T4;
                T2 = pParam.BSIM3dvt2w * T4 * T4;
            }

            ltw = ModelTemperature.Factor1 * T3 * T1;
            dltw_dVb = ModelTemperature.Factor1 * (0.5 / T3 * T1 * dXdep_dVb + T3 * T2);

            T0 = -0.5 * pParam.BSIM3dvt1 * Leff / lt1;
            if (T0 > -EXP_THRESHOLD)
            {
                T1 = Math.Exp(T0);
                Theta0 = T1 * (1.0 + 2.0 * T1);
                dT1_dVb = -T0 / lt1 * T1 * dlt1_dVb;
                dTheta0_dVb = (1.0 + 4.0 * T1) * dT1_dVb;
            }
            else
            {
                T1 = MIN_EXP;
                Theta0 = T1 * (1.0 + 2.0 * T1);
                dTheta0_dVb = 0.0;
            }

            Thetavth = pParam.BSIM3dvt0 * Theta0;
            Delt_vth = Thetavth * V0;
            dDelt_vth_dVb = pParam.BSIM3dvt0 * dTheta0_dVb * V0;

            T0 = -0.5 * pParam.BSIM3dvt1w * pParam.BSIM3weff * Leff / ltw;
            if (T0 > -EXP_THRESHOLD)
            {
                T1 = Math.Exp(T0);
                T2 = T1 * (1.0 + 2.0 * T1);
                dT1_dVb = -T0 / ltw * T1 * dltw_dVb;
                dT2_dVb = (1.0 + 4.0 * T1) * dT1_dVb;
            }
            else
            {
                T1 = MIN_EXP;
                T2 = T1 * (1.0 + 2.0 * T1);
                dT2_dVb = 0.0;
            }

            T0 = pParam.BSIM3dvt0w * T2;
            T2 = T0 * V0;
            dT2_dVb = pParam.BSIM3dvt0w * dT2_dVb * V0;

            TempRatio = state.Temperature / ModelParameters.Tnom - 1.0;
            T0 = Math.Sqrt(1.0 + pParam.BSIM3nlx / Leff);
            T1 = pParam.BSIM3k1ox * (T0 - 1.0) * pParam.BSIM3sqrtPhi
                 + (pParam.BSIM3kt1 + pParam.BSIM3kt1l / Leff
                                    + pParam.BSIM3kt2 * Vbseff) * TempRatio;
            tmp2 = ModelParameters.Tox * pParam.BSIM3phi
                   / (pParam.BSIM3weff + pParam.BSIM3w0);

            T3 = pParam.BSIM3eta0 + pParam.BSIM3etab * Vbseff;
            if (T3 < 1.0e-4) /* avoid  discontinuity problems caused by etab */
            {
                T9 = 1.0 / (3.0 - 2.0e4 * T3);
                T3 = (2.0e-4 - T3) * T9;
                T4 = T9 * T9;
            }
            else
            {
                T4 = 1.0;
            }

            dDIBL_Sft_dVd = T3 * pParam.BSIM3theta0vb0;
            DIBL_Sft = dDIBL_Sft_dVd * Vds;

            Vth = ModelParameters.B3Type * pParam.BSIM3vth0 - pParam.BSIM3k1
                  * pParam.BSIM3sqrtPhi + pParam.BSIM3k1ox * sqrtPhis
                  - pParam.BSIM3k2ox * Vbseff - Delt_vth - T2 + (pParam.BSIM3k3
                                                                 + pParam.BSIM3k3b * Vbseff) * tmp2 + T1 - DIBL_Sft;

            Von = Vth;

            dVth_dVb = pParam.BSIM3k1ox * dsqrtPhis_dVb - pParam.BSIM3k2ox
                                                        - dDelt_vth_dVb - dT2_dVb + pParam.BSIM3k3b * tmp2
                       - pParam.BSIM3etab * Vds * pParam.BSIM3theta0vb0 * T4
                       + pParam.BSIM3kt2 * TempRatio;
            dVth_dVd = -dDIBL_Sft_dVd;

            /* Calculate n */
            tmp2 = pParam.BSIM3nfactor * EPSSI / Xdep;
            tmp3 = pParam.BSIM3cdsc + pParam.BSIM3cdscb * Vbseff
                                    + pParam.BSIM3cdscd * Vds;
            tmp4 = (tmp2 + tmp3 * Theta0 + pParam.BSIM3cit) / ModelParameters.Cox;
            if (tmp4 >= -0.5)
            {
                n = 1.0 + tmp4;
                dn_dVb = (-tmp2 / Xdep * dXdep_dVb + tmp3 * dTheta0_dVb
                                                   + pParam.BSIM3cdscb * Theta0) / ModelParameters.Cox;
                dn_dVd = pParam.BSIM3cdscd * Theta0 / ModelParameters.Cox;
            }
            else
                /* avoid  discontinuity problems caused by tmp4 */
            {
                T0 = 1.0 / (3.0 + 8.0 * tmp4);
                n = (1.0 + 3.0 * tmp4) * T0;
                T0 *= T0;
                dn_dVb = (-tmp2 / Xdep * dXdep_dVb + tmp3 * dTheta0_dVb
                                                   + pParam.BSIM3cdscb * Theta0) / ModelParameters.Cox * T0;
                dn_dVd = pParam.BSIM3cdscd * Theta0 / ModelParameters.Cox * T0;
            }

            /* Poly Gate Si Depletion Effect */
            T0 = pParam.BSIM3vfb + pParam.BSIM3phi;
            if (pParam.BSIM3ngate > 1.0e18 && pParam.BSIM3ngate < 1.0e25
                                           && Vgs > T0)
                /* added to avoid the problem caused by ngate */
            {
                T1 = 1.0e6 * Constants.Charge * EPSSI * pParam.BSIM3ngate
                     / (ModelParameters.Cox * ModelParameters.Cox);
                T4 = Math.Sqrt(1.0 + 2.0 * (Vgs - T0) / T1);
                T2 = T1 * (T4 - 1.0);
                T3 = 0.5 * T2 * T2 / T1; /* T3 = Vpoly */
                T7 = 1.12 - T3 - 0.05;
                T6 = Math.Sqrt(T7 * T7 + 0.224);
                T5 = 1.12 - 0.5 * (T7 + T6);
                Vgs_eff = Vgs - T5;
                dVgs_eff_dVg = 1.0 - (0.5 - 0.5 / T4) * (1.0 + T7 / T6);
            }
            else
            {
                Vgs_eff = Vgs;
                dVgs_eff_dVg = 1.0;
            }

            Vgst = Vgs_eff - Vth;

/* Effective Vgst (Vgsteff) Calculation */

            T10 = 2.0 * n * Vtm;
            VgstNVt = Vgst / T10;
            ExpArg = (2.0 * pParam.BSIM3voff - Vgst) / T10;

            /* MCJ: Very small Vgst */
            if (VgstNVt > EXP_THRESHOLD)
            {
                Vgsteff = Vgst;
                dVgsteff_dVg = dVgs_eff_dVg;
                dVgsteff_dVd = -dVth_dVd;
                dVgsteff_dVb = -dVth_dVb;
            }
            else if (ExpArg > EXP_THRESHOLD)
            {
                T0 = (Vgst - pParam.BSIM3voff) / (n * Vtm);
                ExpVgst = Math.Exp(T0);
                Vgsteff = Vtm * pParam.BSIM3cdep0 / ModelParameters.Cox * ExpVgst;
                dVgsteff_dVg = Vgsteff / (n * Vtm);
                dVgsteff_dVd = -dVgsteff_dVg * (dVth_dVd + T0 * Vtm * dn_dVd);
                dVgsteff_dVb = -dVgsteff_dVg * (dVth_dVb + T0 * Vtm * dn_dVb);
                dVgsteff_dVg *= dVgs_eff_dVg;
            }
            else
            {
                ExpVgst = Math.Exp(VgstNVt);
                T1 = T10 * Math.Log(1.0 + ExpVgst);
                dT1_dVg = ExpVgst / (1.0 + ExpVgst);
                dT1_dVb = -dT1_dVg * (dVth_dVb + Vgst / n * dn_dVb)
                          + T1 / n * dn_dVb;
                dT1_dVd = -dT1_dVg * (dVth_dVd + Vgst / n * dn_dVd)
                          + T1 / n * dn_dVd;

                dT2_dVg = -ModelParameters.Cox / (Vtm * pParam.BSIM3cdep0)
                          * Math.Exp(ExpArg);
                T2 = 1.0 - T10 * dT2_dVg;
                dT2_dVd = -dT2_dVg * (dVth_dVd - 2.0 * Vtm * ExpArg * dn_dVd)
                          + (T2 - 1.0) / n * dn_dVd;
                dT2_dVb = -dT2_dVg * (dVth_dVb - 2.0 * Vtm * ExpArg * dn_dVb)
                          + (T2 - 1.0) / n * dn_dVb;

                Vgsteff = T1 / T2;
                T3 = T2 * T2;
                dVgsteff_dVg = (T2 * dT1_dVg - T1 * dT2_dVg) / T3 * dVgs_eff_dVg;
                dVgsteff_dVd = (T2 * dT1_dVd - T1 * dT2_dVd) / T3;
                dVgsteff_dVb = (T2 * dT1_dVb - T1 * dT2_dVb) / T3;
            }

            this.Vgsteff = Vgsteff;

/* Calculate Effective Channel Geometry */
            T9 = sqrtPhis - pParam.BSIM3sqrtPhi;
            Weff = pParam.BSIM3weff - 2.0 * (pParam.BSIM3dwg * Vgsteff
                                             + pParam.BSIM3dwb * T9);
            dWeff_dVg = -2.0 * pParam.BSIM3dwg;
            dWeff_dVb = -2.0 * pParam.BSIM3dwb * dsqrtPhis_dVb;

            if (Weff < 2.0e-8) /* to avoid the discontinuity problem due to Weff*/
            {
                T0 = 1.0 / (6.0e-8 - 2.0 * Weff);
                Weff = 2.0e-8 * (4.0e-8 - Weff) * T0;
                T0 *= T0 * 4.0e-16;
                dWeff_dVg *= T0;
                dWeff_dVb *= T0;
            }

            T0 = pParam.BSIM3prwg * Vgsteff + pParam.BSIM3prwb * T9;
            if (T0 >= -0.9)
            {
                Rds = pParam.BSIM3rds0 * (1.0 + T0);
                dRds_dVg = pParam.BSIM3rds0 * pParam.BSIM3prwg;
                dRds_dVb = pParam.BSIM3rds0 * pParam.BSIM3prwb * dsqrtPhis_dVb;
            }
            else
                /* to avoid the discontinuity problem due to prwg and prwb*/
            {
                T1 = 1.0 / (17.0 + 20.0 * T0);
                Rds = pParam.BSIM3rds0 * (0.8 + T0) * T1;
                T1 *= T1;
                dRds_dVg = pParam.BSIM3rds0 * pParam.BSIM3prwg * T1;
                dRds_dVb = pParam.BSIM3rds0 * pParam.BSIM3prwb * dsqrtPhis_dVb
                           * T1;
            }

            this.Rds = Rds; /* Noise Bugfix */

/* Calculate Abulk */
            T1 = 0.5 * pParam.BSIM3k1ox / sqrtPhis;
            dT1_dVb = -T1 / sqrtPhis * dsqrtPhis_dVb;

            T9 = Math.Sqrt(pParam.BSIM3xj * Xdep);
            tmp1 = Leff + 2.0 * T9;
            T5 = Leff / tmp1;
            tmp2 = pParam.BSIM3a0 * T5;
            tmp3 = pParam.BSIM3weff + pParam.BSIM3b1;
            tmp4 = pParam.BSIM3b0 / tmp3;
            T2 = tmp2 + tmp4;
            dT2_dVb = -T9 / tmp1 / Xdep * dXdep_dVb;
            T6 = T5 * T5;
            T7 = T5 * T6;

            Abulk0 = 1.0 + T1 * T2;
            dAbulk0_dVb = T1 * tmp2 * dT2_dVb + T2 * dT1_dVb;

            T8 = pParam.BSIM3ags * pParam.BSIM3a0 * T7;
            dAbulk_dVg = -T1 * T8;
            Abulk = Abulk0 + dAbulk_dVg * Vgsteff;
            dAbulk_dVb = dAbulk0_dVb - T8 * Vgsteff * (dT1_dVb
                                                       + 3.0 * T1 * dT2_dVb);

            if (Abulk0 < 0.1) /* added to avoid the problems caused by Abulk0 */
            {
                T9 = 1.0 / (3.0 - 20.0 * Abulk0);
                Abulk0 = (0.2 - Abulk0) * T9;
                dAbulk0_dVb *= T9 * T9;
            }

            if (Abulk < 0.1)
                /* added to avoid the problems caused by Abulk */
            {
                T9 = 1.0 / (3.0 - 20.0 * Abulk);
                Abulk = (0.2 - Abulk) * T9;
                T10 = T9 * T9;
                dAbulk_dVb *= T10;
                dAbulk_dVg *= T10;
            }

            this.Abulk = Abulk;

            T2 = pParam.BSIM3keta * Vbseff;
            if (T2 >= -0.9)
            {
                T0 = 1.0 / (1.0 + T2);
                dT0_dVb = -pParam.BSIM3keta * T0 * T0;
            }
            else
                /* added to avoid the problems caused by Keta */
            {
                T1 = 1.0 / (0.8 + T2);
                T0 = (17.0 + 20.0 * T2) * T1;
                dT0_dVb = -pParam.BSIM3keta * T1 * T1;
            }

            dAbulk_dVg *= T0;
            dAbulk_dVb = dAbulk_dVb * T0 + Abulk * dT0_dVb;
            dAbulk0_dVb = dAbulk0_dVb * T0 + Abulk0 * dT0_dVb;
            Abulk *= T0;
            Abulk0 *= T0;


/* Mobility calculation */
            if (ModelParameters.MobMod == 1)
            {
                T0 = Vgsteff + Vth + Vth;
                T2 = pParam.BSIM3ua + pParam.BSIM3uc * Vbseff;
                T3 = T0 / ModelParameters.Tox;
                T5 = T3 * (T2 + pParam.BSIM3ub * T3);
                dDenomi_dVg = (T2 + 2.0 * pParam.BSIM3ub * T3) / ModelParameters.Tox;
                dDenomi_dVd = dDenomi_dVg * 2.0 * dVth_dVd;
                dDenomi_dVb = dDenomi_dVg * 2.0 * dVth_dVb + pParam.BSIM3uc * T3;
            }
            else if (ModelParameters.MobMod == 2)
            {
                T5 = Vgsteff / ModelParameters.Tox * (pParam.BSIM3ua
                                           + pParam.BSIM3uc * Vbseff + pParam.BSIM3ub * Vgsteff
                                           / ModelParameters.Tox);
                dDenomi_dVg = (pParam.BSIM3ua + pParam.BSIM3uc * Vbseff
                                              + 2.0 * pParam.BSIM3ub * Vgsteff / ModelParameters.Tox)
                              / ModelParameters.Tox;
                dDenomi_dVd = 0.0;
                dDenomi_dVb = Vgsteff * pParam.BSIM3uc / ModelParameters.Tox;
            }
            else
            {
                T0 = Vgsteff + Vth + Vth;
                T2 = 1.0 + pParam.BSIM3uc * Vbseff;
                T3 = T0 / ModelParameters.Tox;
                T4 = T3 * (pParam.BSIM3ua + pParam.BSIM3ub * T3);
                T5 = T4 * T2;
                dDenomi_dVg = (pParam.BSIM3ua + 2.0 * pParam.BSIM3ub * T3) * T2
                              / ModelParameters.Tox;
                dDenomi_dVd = dDenomi_dVg * 2.0 * dVth_dVd;
                dDenomi_dVb = dDenomi_dVg * 2.0 * dVth_dVb + pParam.BSIM3uc * T4;
            }

            if (T5 >= -0.8)
            {
                Denomi = 1.0 + T5;
            }
            else /* Added to avoid the discontinuity problem caused by ua and ub*/
            {
                T9 = 1.0 / (7.0 + 10.0 * T5);
                Denomi = (0.6 + T5) * T9;
                T9 *= T9;
                dDenomi_dVg *= T9;
                dDenomi_dVd *= T9;
                dDenomi_dVb *= T9;
            }

            Ueff = ueff = pParam.BSIM3u0temp / Denomi;
            T9 = -ueff / Denomi;
            dueff_dVg = T9 * dDenomi_dVg;
            dueff_dVd = T9 * dDenomi_dVd;
            dueff_dVb = T9 * dDenomi_dVb;

/* Saturation Drain Voltage  Vdsat */
            WVCox = Weff * pParam.BSIM3vsattemp * ModelParameters.Cox;
            WVCoxRds = WVCox * Rds;

            Esat = 2.0 * pParam.BSIM3vsattemp / ueff;
            EsatL = Esat * Leff;
            T0 = -EsatL / ueff;
            dEsatL_dVg = T0 * dueff_dVg;
            dEsatL_dVd = T0 * dueff_dVd;
            dEsatL_dVb = T0 * dueff_dVb;

            /* Sqrt() */
            a1 = pParam.BSIM3a1;
            if (a1 == 0.0)
            {
                Lambda = pParam.BSIM3a2;
                dLambda_dVg = 0.0;
            }
            else if (a1 > 0.0)
/* Added to avoid the discontinuity problem
   caused by a1 and a2 (Lambda) */
            {
                T0 = 1.0 - pParam.BSIM3a2;
                T1 = T0 - pParam.BSIM3a1 * Vgsteff - 0.0001;
                T2 = Math.Sqrt(T1 * T1 + 0.0004 * T0);
                Lambda = pParam.BSIM3a2 + T0 - 0.5 * (T1 + T2);
                dLambda_dVg = 0.5 * pParam.BSIM3a1 * (1.0 + T1 / T2);
            }
            else
            {
                T1 = pParam.BSIM3a2 + pParam.BSIM3a1 * Vgsteff - 0.0001;
                T2 = Math.Sqrt(T1 * T1 + 0.0004 * pParam.BSIM3a2);
                Lambda = 0.5 * (T1 + T2);
                dLambda_dVg = 0.5 * pParam.BSIM3a1 * (1.0 + T1 / T2);
            }

            Vgst2Vtm = Vgsteff + 2.0 * Vtm;
            AbovVgst2Vtm = Abulk / Vgst2Vtm;

            if (Rds > 0)
            {
                tmp2 = dRds_dVg / Rds + dWeff_dVg / Weff;
                tmp3 = dRds_dVb / Rds + dWeff_dVb / Weff;
            }
            else
            {
                tmp2 = dWeff_dVg / Weff;
                tmp3 = dWeff_dVb / Weff;
            }

            if (Rds == 0.0 && Lambda == 1.0)
            {
                T0 = 1.0 / (Abulk * EsatL + Vgst2Vtm);
                tmp1 = 0.0;
                T1 = T0 * T0;
                T2 = Vgst2Vtm * T0;
                T3 = EsatL * Vgst2Vtm;
                Vdsat = T3 * T0;

                dT0_dVg = -(Abulk * dEsatL_dVg + EsatL * dAbulk_dVg + 1.0) * T1;
                dT0_dVd = -(Abulk * dEsatL_dVd) * T1;
                dT0_dVb = -(Abulk * dEsatL_dVb + dAbulk_dVb * EsatL) * T1;

                dVdsat_dVg = T3 * dT0_dVg + T2 * dEsatL_dVg + EsatL * T0;
                dVdsat_dVd = T3 * dT0_dVd + T2 * dEsatL_dVd;
                dVdsat_dVb = T3 * dT0_dVb + T2 * dEsatL_dVb;
            }
            else
            {
                tmp1 = dLambda_dVg / (Lambda * Lambda);
                T9 = Abulk * WVCoxRds;
                T8 = Abulk * T9;
                T7 = Vgst2Vtm * T9;
                T6 = Vgst2Vtm * WVCoxRds;
                T0 = 2.0 * Abulk * (T9 - 1.0 + 1.0 / Lambda);
                dT0_dVg = 2.0 * (T8 * tmp2 - Abulk * tmp1
                                 + (2.0 * T9 + 1.0 / Lambda - 1.0) * dAbulk_dVg);

                dT0_dVb = 2.0 * (T8 * (2.0 / Abulk * dAbulk_dVb + tmp3)
                                 + (1.0 / Lambda - 1.0) * dAbulk_dVb);
                dT0_dVd = 0.0;
                T1 = Vgst2Vtm * (2.0 / Lambda - 1.0) + Abulk * EsatL + 3.0 * T7;

                dT1_dVg = 2.0 / Lambda - 1.0 - 2.0 * Vgst2Vtm * tmp1
                          + Abulk * dEsatL_dVg + EsatL * dAbulk_dVg + 3.0 * (T9
                                                                             + T7 * tmp2 + T6 * dAbulk_dVg);
                dT1_dVb = Abulk * dEsatL_dVb + EsatL * dAbulk_dVb
                                             + 3.0 * (T6 * dAbulk_dVb + T7 * tmp3);
                dT1_dVd = Abulk * dEsatL_dVd;

                T2 = Vgst2Vtm * (EsatL + 2.0 * T6);
                dT2_dVg = EsatL + Vgst2Vtm * dEsatL_dVg
                                + T6 * (4.0 + 2.0 * Vgst2Vtm * tmp2);
                dT2_dVb = Vgst2Vtm * (dEsatL_dVb + 2.0 * T6 * tmp3);
                dT2_dVd = Vgst2Vtm * dEsatL_dVd;

                T3 = Math.Sqrt(T1 * T1 - 2.0 * T0 * T2);
                Vdsat = (T1 - T3) / T0;

                dT3_dVg = (T1 * dT1_dVg - 2.0 * (T0 * dT2_dVg + T2 * dT0_dVg))
                          / T3;
                dT3_dVd = (T1 * dT1_dVd - 2.0 * (T0 * dT2_dVd + T2 * dT0_dVd))
                          / T3;
                dT3_dVb = (T1 * dT1_dVb - 2.0 * (T0 * dT2_dVb + T2 * dT0_dVb))
                          / T3;

                dVdsat_dVg = (dT1_dVg - (T1 * dT1_dVg - dT0_dVg * T2
                                                      - T0 * dT2_dVg) / T3 - Vdsat * dT0_dVg) / T0;
                dVdsat_dVb = (dT1_dVb - (T1 * dT1_dVb - dT0_dVb * T2
                                                      - T0 * dT2_dVb) / T3 - Vdsat * dT0_dVb) / T0;
                dVdsat_dVd = (dT1_dVd - (T1 * dT1_dVd - T0 * dT2_dVd) / T3) / T0;
            }

            this.Vdsat = Vdsat;

/* Effective Vds (Vdseff) Calculation */
            T1 = Vdsat - Vds - pParam.BSIM3delta;
            dT1_dVg = dVdsat_dVg;
            dT1_dVd = dVdsat_dVd - 1.0;
            dT1_dVb = dVdsat_dVb;

            T2 = Math.Sqrt(T1 * T1 + 4.0 * pParam.BSIM3delta * Vdsat);
            T0 = T1 / T2;
            T3 = 2.0 * pParam.BSIM3delta / T2;
            dT2_dVg = T0 * dT1_dVg + T3 * dVdsat_dVg;
            dT2_dVd = T0 * dT1_dVd + T3 * dVdsat_dVd;
            dT2_dVb = T0 * dT1_dVb + T3 * dVdsat_dVb;

            Vdseff = Vdsat - 0.5 * (T1 + T2);
            dVdseff_dVg = dVdsat_dVg - 0.5 * (dT1_dVg + dT2_dVg);
            dVdseff_dVd = dVdsat_dVd - 0.5 * (dT1_dVd + dT2_dVd);
            dVdseff_dVb = dVdsat_dVb - 0.5 * (dT1_dVb + dT2_dVb);
            /* Added to eliminate non-zero Vdseff at Vds=0.0 */
            if (Vds == 0.0)
            {
                Vdseff = 0.0;
                dVdseff_dVg = 0.0;
                dVdseff_dVb = 0.0;
            }

/* Calculate VAsat */
            tmp4 = 1.0 - 0.5 * Abulk * Vdsat / Vgst2Vtm;
            T9 = WVCoxRds * Vgsteff;
            T8 = T9 / Vgst2Vtm;
            T0 = EsatL + Vdsat + 2.0 * T9 * tmp4;

            T7 = 2.0 * WVCoxRds * tmp4;
            dT0_dVg = dEsatL_dVg + dVdsat_dVg + T7 * (1.0 + tmp2 * Vgsteff)
                      - T8 * (Abulk * dVdsat_dVg - Abulk * Vdsat / Vgst2Vtm
                              + Vdsat * dAbulk_dVg);

            dT0_dVb = dEsatL_dVb + dVdsat_dVb + T7 * tmp3 * Vgsteff
                      - T8 * (dAbulk_dVb * Vdsat + Abulk * dVdsat_dVb);
            dT0_dVd = dEsatL_dVd + dVdsat_dVd - T8 * Abulk * dVdsat_dVd;

            T9 = WVCoxRds * Abulk;
            T1 = 2.0 / Lambda - 1.0 + T9;
            dT1_dVg = -2.0 * tmp1 + WVCoxRds * (Abulk * tmp2 + dAbulk_dVg);
            dT1_dVb = dAbulk_dVb * WVCoxRds + T9 * tmp3;

            Vasat = T0 / T1;
            dVasat_dVg = (dT0_dVg - Vasat * dT1_dVg) / T1;
            dVasat_dVb = (dT0_dVb - Vasat * dT1_dVb) / T1;
            dVasat_dVd = dT0_dVd / T1;

            if (Vdseff > Vds)
                Vdseff = Vds;
            diffVds = Vds - Vdseff;
            this.Vdseff = Vdseff;

/* Calculate VACLM */
            if (pParam.BSIM3pclm > 0.0 && diffVds > 1.0e-10)
            {
                T0 = 1.0 / (pParam.BSIM3pclm * Abulk * pParam.BSIM3litl);
                dT0_dVb = -T0 / Abulk * dAbulk_dVb;
                dT0_dVg = -T0 / Abulk * dAbulk_dVg;

                T2 = Vgsteff / EsatL;
                T1 = Leff * (Abulk + T2);
                dT1_dVg = Leff * ((1.0 - T2 * dEsatL_dVg) / EsatL + dAbulk_dVg);
                dT1_dVb = Leff * (dAbulk_dVb - T2 * dEsatL_dVb / EsatL);
                dT1_dVd = -T2 * dEsatL_dVd / Esat;

                T9 = T0 * T1;
                VACLM = T9 * diffVds;
                dVACLM_dVg = T0 * dT1_dVg * diffVds - T9 * dVdseff_dVg
                             + T1 * diffVds * dT0_dVg;
                dVACLM_dVb = (dT0_dVb * T1 + T0 * dT1_dVb) * diffVds
                             - T9 * dVdseff_dVb;
                dVACLM_dVd = T0 * dT1_dVd * diffVds + T9 * (1.0 - dVdseff_dVd);
            }
            else
            {
                VACLM = MAX_EXP;
                dVACLM_dVd = dVACLM_dVg = dVACLM_dVb = 0.0;
            }

/* Calculate VADIBL */
            if (pParam.BSIM3thetaRout > 0.0)
            {
                T8 = Abulk * Vdsat;
                T0 = Vgst2Vtm * T8;
                dT0_dVg = Vgst2Vtm * Abulk * dVdsat_dVg + T8
                                                        + Vgst2Vtm * Vdsat * dAbulk_dVg;
                dT0_dVb = Vgst2Vtm * (dAbulk_dVb * Vdsat + Abulk * dVdsat_dVb);
                dT0_dVd = Vgst2Vtm * Abulk * dVdsat_dVd;

                T1 = Vgst2Vtm + T8;
                dT1_dVg = 1.0 + Abulk * dVdsat_dVg + Vdsat * dAbulk_dVg;
                dT1_dVb = Abulk * dVdsat_dVb + dAbulk_dVb * Vdsat;
                dT1_dVd = Abulk * dVdsat_dVd;

                T9 = T1 * T1;
                T2 = pParam.BSIM3thetaRout;
                VADIBL = (Vgst2Vtm - T0 / T1) / T2;
                dVADIBL_dVg = (1.0 - dT0_dVg / T1 + T0 * dT1_dVg / T9) / T2;
                dVADIBL_dVb = (-dT0_dVb / T1 + T0 * dT1_dVb / T9) / T2;
                dVADIBL_dVd = (-dT0_dVd / T1 + T0 * dT1_dVd / T9) / T2;

                T7 = pParam.BSIM3pdiblb * Vbseff;
                if (T7 >= -0.9)
                {
                    T3 = 1.0 / (1.0 + T7);
                    VADIBL *= T3;
                    dVADIBL_dVg *= T3;
                    dVADIBL_dVb = (dVADIBL_dVb - VADIBL * pParam.BSIM3pdiblb)
                                  * T3;
                    dVADIBL_dVd *= T3;
                }
                else
/* Added to avoid the discontinuity problem caused by pdiblcb */
                {
                    T4 = 1.0 / (0.8 + T7);
                    T3 = (17.0 + 20.0 * T7) * T4;
                    dVADIBL_dVg *= T3;
                    dVADIBL_dVb = dVADIBL_dVb * T3
                                  - VADIBL * pParam.BSIM3pdiblb * T4 * T4;
                    dVADIBL_dVd *= T3;
                    VADIBL *= T3;
                }
            }
            else
            {
                VADIBL = MAX_EXP;
                dVADIBL_dVd = dVADIBL_dVg = dVADIBL_dVb = 0.0;
            }

/* Calculate VA */

            T8 = pParam.BSIM3pvag / EsatL;
            T9 = T8 * Vgsteff;
            if (T9 > -0.9)
            {
                T0 = 1.0 + T9;
                dT0_dVg = T8 * (1.0 - Vgsteff * dEsatL_dVg / EsatL);
                dT0_dVb = -T9 * dEsatL_dVb / EsatL;
                dT0_dVd = -T9 * dEsatL_dVd / EsatL;
            }
            else /* Added to avoid the discontinuity problems caused by pvag */
            {
                T1 = 1.0 / (17.0 + 20.0 * T9);
                T0 = (0.8 + T9) * T1;
                T1 *= T1;
                dT0_dVg = T8 * (1.0 - Vgsteff * dEsatL_dVg / EsatL) * T1;

                T9 *= T1 / EsatL;
                dT0_dVb = -T9 * dEsatL_dVb;
                dT0_dVd = -T9 * dEsatL_dVd;
            }

            tmp1 = VACLM * VACLM;
            tmp2 = VADIBL * VADIBL;
            tmp3 = VACLM + VADIBL;

            T1 = VACLM * VADIBL / tmp3;
            tmp3 *= tmp3;
            dT1_dVg = (tmp1 * dVADIBL_dVg + tmp2 * dVACLM_dVg) / tmp3;
            dT1_dVd = (tmp1 * dVADIBL_dVd + tmp2 * dVACLM_dVd) / tmp3;
            dT1_dVb = (tmp1 * dVADIBL_dVb + tmp2 * dVACLM_dVb) / tmp3;

            Va = Vasat + T0 * T1;
            dVa_dVg = dVasat_dVg + T1 * dT0_dVg + T0 * dT1_dVg;
            dVa_dVd = dVasat_dVd + T1 * dT0_dVd + T0 * dT1_dVd;
            dVa_dVb = dVasat_dVb + T1 * dT0_dVb + T0 * dT1_dVb;

/* Calculate VASCBE */
            if (pParam.BSIM3pscbe2 > 0.0)
            {
                if (diffVds > pParam.BSIM3pscbe1 * pParam.BSIM3litl
                    / EXP_THRESHOLD)
                {
                    T0 = pParam.BSIM3pscbe1 * pParam.BSIM3litl / diffVds;
                    VASCBE = Leff * Math.Exp(T0) / pParam.BSIM3pscbe2;
                    T1 = T0 * VASCBE / diffVds;
                    dVASCBE_dVg = T1 * dVdseff_dVg;
                    dVASCBE_dVd = -T1 * (1.0 - dVdseff_dVd);
                    dVASCBE_dVb = T1 * dVdseff_dVb;
                }
                else
                {
                    VASCBE = MAX_EXP * Leff / pParam.BSIM3pscbe2;
                    dVASCBE_dVg = dVASCBE_dVd = dVASCBE_dVb = 0.0;
                }
            }
            else
            {
                VASCBE = MAX_EXP;
                dVASCBE_dVg = dVASCBE_dVd = dVASCBE_dVb = 0.0;
            }

/* Calculate Ids */
            CoxWovL = ModelParameters.Cox * Weff / Leff;
            beta = ueff * CoxWovL;
            dbeta_dVg = CoxWovL * dueff_dVg + beta * dWeff_dVg / Weff;
            dbeta_dVd = CoxWovL * dueff_dVd;
            dbeta_dVb = CoxWovL * dueff_dVb + beta * dWeff_dVb / Weff;

            T0 = 1.0 - 0.5 * Abulk * Vdseff / Vgst2Vtm;
            dT0_dVg = -0.5 * (Abulk * dVdseff_dVg
                              - Abulk * Vdseff / Vgst2Vtm + Vdseff * dAbulk_dVg) / Vgst2Vtm;
            dT0_dVd = -0.5 * Abulk * dVdseff_dVd / Vgst2Vtm;
            dT0_dVb = -0.5 * (Abulk * dVdseff_dVb + dAbulk_dVb * Vdseff)
                      / Vgst2Vtm;

            fgche1 = Vgsteff * T0;
            dfgche1_dVg = Vgsteff * dT0_dVg + T0;
            dfgche1_dVd = Vgsteff * dT0_dVd;
            dfgche1_dVb = Vgsteff * dT0_dVb;

            T9 = Vdseff / EsatL;
            fgche2 = 1.0 + T9;
            dfgche2_dVg = (dVdseff_dVg - T9 * dEsatL_dVg) / EsatL;
            dfgche2_dVd = (dVdseff_dVd - T9 * dEsatL_dVd) / EsatL;
            dfgche2_dVb = (dVdseff_dVb - T9 * dEsatL_dVb) / EsatL;

            gche = beta * fgche1 / fgche2;
            dgche_dVg = (beta * dfgche1_dVg + fgche1 * dbeta_dVg
                         - gche * dfgche2_dVg) / fgche2;
            dgche_dVd = (beta * dfgche1_dVd + fgche1 * dbeta_dVd
                         - gche * dfgche2_dVd) / fgche2;
            dgche_dVb = (beta * dfgche1_dVb + fgche1 * dbeta_dVb
                         - gche * dfgche2_dVb) / fgche2;

            T0 = 1.0 + gche * Rds;
            T9 = Vdseff / T0;
            Idl = gche * T9;

            dIdl_dVg = (gche * dVdseff_dVg + T9 * dgche_dVg) / T0
                       - Idl * gche / T0 * dRds_dVg;

            dIdl_dVd = (gche * dVdseff_dVd + T9 * dgche_dVd) / T0;
            dIdl_dVb = (gche * dVdseff_dVb + T9 * dgche_dVb
                        - Idl * dRds_dVb * gche) / T0;

            T9 = diffVds / Va;
            T0 = 1.0 + T9;
            Idsa = Idl * T0;
            dIdsa_dVg = T0 * dIdl_dVg - Idl * (dVdseff_dVg + T9 * dVa_dVg) / Va;
            dIdsa_dVd = T0 * dIdl_dVd + Idl * (1.0 - dVdseff_dVd
                                                   - T9 * dVa_dVd) / Va;
            dIdsa_dVb = T0 * dIdl_dVb - Idl * (dVdseff_dVb + T9 * dVa_dVb) / Va;

            T9 = diffVds / VASCBE;
            T0 = 1.0 + T9;
            Ids = Idsa * T0;

            Gm = T0 * dIdsa_dVg - Idsa * (dVdseff_dVg + T9 * dVASCBE_dVg) / VASCBE;
            Gds = T0 * dIdsa_dVd + Idsa * (1.0 - dVdseff_dVd
                                               - T9 * dVASCBE_dVd) / VASCBE;
            Gmb = T0 * dIdsa_dVb - Idsa * (dVdseff_dVb
                                           + T9 * dVASCBE_dVb) / VASCBE;

            Gds += Gm * dVgsteff_dVd;
            Gmb += Gm * dVgsteff_dVb;
            Gm *= dVgsteff_dVg;
            Gmb *= dVbseff_dVb;

            /* Substrate current begins */
            tmp = pParam.BSIM3alpha0 + pParam.BSIM3alpha1 * Leff;
            if (tmp <= 0.0 || pParam.BSIM3beta0 <= 0.0)
            {
                Isub = Gbd = Gbb = Gbg = 0.0;
            }
            else
            {
                T2 = tmp / Leff;
                if (diffVds > pParam.BSIM3beta0 / EXP_THRESHOLD)
                {
                    T0 = -pParam.BSIM3beta0 / diffVds;
                    T1 = T2 * diffVds * Math.Exp(T0);
                    T3 = T1 / diffVds * (T0 - 1.0);
                    dT1_dVg = T3 * dVdseff_dVg;
                    dT1_dVd = T3 * (dVdseff_dVd - 1.0);
                    dT1_dVb = T3 * dVdseff_dVb;
                }
                else
                {
                    T3 = T2 * MIN_EXP;
                    T1 = T3 * diffVds;
                    dT1_dVg = -T3 * dVdseff_dVg;
                    dT1_dVd = T3 * (1.0 - dVdseff_dVd);
                    dT1_dVb = -T3 * dVdseff_dVb;
                }

                Isub = T1 * Idsa;
                Gbg = T1 * dIdsa_dVg + Idsa * dT1_dVg;
                Gbd = T1 * dIdsa_dVd + Idsa * dT1_dVd;
                Gbb = T1 * dIdsa_dVb + Idsa * dT1_dVb;

                Gbd += Gbg * dVgsteff_dVd;
                Gbb += Gbg * dVgsteff_dVb;
                Gbg *= dVgsteff_dVg;
                Gbb *= dVbseff_dVb; /* bug fixing */
            }

            cdrain = Ids;
            this.Gds = Gds;
            this.Gm = Gm;
            this.Gmbs = Gmb;

            Gbbs = Gbb;
            Gbgs = Gbg;
            Gbds = Gbd;

            Csub = Isub;

            /* BSIM3 thermal noise Qinv calculated from all capMod 
             * 0, 1, 2 & 3 stored in BaseParameters.qinv 1/1998 */

            if (ModelParameters.Xpart < 0 || !chargeComputationNeeded)
            {
                qgate = qdrn = qsrc = qbulk = 0.0;
                Cggb = Cgsb = Cgdb = 0.0;
                Cdgb = Cdsb = Cddb = 0.0;
                Cbgb = Cbsb = Cbdb = 0.0;
                Cqdb = Cqsb = Cqgb
                    = Cqbb = 0.0;
                Gtau = 0.0;
                goto finished;
            }
            else if (ModelParameters.CapMod == 0)
            {
                if (Vbseff < 0.0)
                {
                    Vbseff = Vbs;
                    dVbseff_dVb = 1.0;
                }
                else
                {
                    Vbseff = pParam.BSIM3phi - Phis;
                    dVbseff_dVb = -dPhis_dVb;
                }

                Vfb = pParam.BSIM3vfbcv;
                Vth = Vfb + pParam.BSIM3phi + pParam.BSIM3k1ox * sqrtPhis;
                Vgst = Vgs_eff - Vth;
                dVth_dVb = pParam.BSIM3k1ox * dsqrtPhis_dVb;
                dVgst_dVb = -dVth_dVb;
                dVgst_dVg = dVgs_eff_dVg;

                CoxWL = ModelParameters.Cox * pParam.BSIM3weffCV
                                 * pParam.BSIM3leffCV;
                Arg1 = Vgs_eff - Vbseff - Vfb;

                if (Arg1 <= 0.0)
                {
                    qgate = CoxWL * Arg1;
                    qbulk = -qgate;
                    qdrn = 0.0;

                    Cggb = CoxWL * dVgs_eff_dVg;
                    Cgdb = 0.0;
                    Cgsb = CoxWL * (dVbseff_dVb - dVgs_eff_dVg);

                    Cdgb = 0.0;
                    Cddb = 0.0;
                    Cdsb = 0.0;

                    Cbgb = -CoxWL * dVgs_eff_dVg;
                    Cbdb = 0.0;
                    Cbsb = -Cgsb;
                    Qinv = 0.0;
                }
                else if (Vgst <= 0.0)
                {
                    T1 = 0.5 * pParam.BSIM3k1ox;
                    T2 = Math.Sqrt(T1 * T1 + Arg1);
                    qgate = CoxWL * pParam.BSIM3k1ox * (T2 - T1);
                    qbulk = -qgate;
                    qdrn = 0.0;

                    T0 = CoxWL * T1 / T2;
                    Cggb = T0 * dVgs_eff_dVg;
                    Cgdb = 0.0;
                    Cgsb = T0 * (dVbseff_dVb - dVgs_eff_dVg);

                    Cdgb = 0.0;
                    Cddb = 0.0;
                    Cdsb = 0.0;

                    Cbgb = -Cggb;
                    Cbdb = 0.0;
                    Cbsb = -Cgsb;
                    Qinv = 0.0;
                }
                else
                {
                    One_Third_CoxWL = CoxWL / 3.0;
                    Two_Third_CoxWL = 2.0 * One_Third_CoxWL;

                    AbulkCV = Abulk0 * pParam.BSIM3abulkCVfactor;
                    dAbulkCV_dVb = pParam.BSIM3abulkCVfactor * dAbulk0_dVb;
                    Vdsat = Vgst / AbulkCV;
                    dVdsat_dVg = dVgs_eff_dVg / AbulkCV;
                    dVdsat_dVb = -(Vdsat * dAbulkCV_dVb + dVth_dVb) / AbulkCV;

                    if (ModelParameters.Xpart > 0.5)
                    {
                        /* 0/100 Charge partition model */
                        if (Vdsat <= Vds)
                        {
                            /* saturation region */
                            T1 = Vdsat / 3.0;
                            qgate = CoxWL * (Vgs_eff - Vfb
                                                     - pParam.BSIM3phi - T1);
                            T2 = -Two_Third_CoxWL * Vgst;
                            qbulk = -(qgate + T2);
                            qdrn = 0.0;

                            Cggb = One_Third_CoxWL * (3.0
                                                      - dVdsat_dVg) * dVgs_eff_dVg;
                            T2 = -One_Third_CoxWL * dVdsat_dVb;
                            Cgsb = -(Cggb + T2);
                            Cgdb = 0.0;

                            Cdgb = 0.0;
                            Cddb = 0.0;
                            Cdsb = 0.0;

                            Cbgb = -(Cggb
                                     - Two_Third_CoxWL * dVgs_eff_dVg);
                            T3 = -(T2 + Two_Third_CoxWL * dVth_dVb);
                            Cbsb = -(Cbgb + T3);
                            Cbdb = 0.0;
                            Qinv = -(qgate + qbulk);
                        }
                        else
                        {
                            /* linear region */
                            Alphaz = Vgst / Vdsat;
                            T1 = 2.0 * Vdsat - Vds;
                            T2 = Vds / (3.0 * T1);
                            T3 = T2 * Vds;
                            T9 = 0.25 * CoxWL;
                            T4 = T9 * Alphaz;
                            T7 = 2.0 * Vds - T1 - 3.0 * T3;
                            T8 = T3 - T1 - 2.0 * Vds;
                            qgate = CoxWL * (Vgs_eff - Vfb
                                                     - pParam.BSIM3phi - 0.5 * (Vds - T3));
                            T10 = T4 * T8;
                            qdrn = T4 * T7;
                            qbulk = -(qgate + qdrn + T10);

                            T5 = T3 / T1;
                            Cggb = CoxWL * (1.0 - T5 * dVdsat_dVg)
                                         * dVgs_eff_dVg;
                            T11 = -CoxWL * T5 * dVdsat_dVb;
                            Cgdb = CoxWL * (T2 - 0.5 + 0.5 * T5);
                            Cgsb = -(Cggb + T11
                                          + Cgdb);
                            T6 = 1.0 / Vdsat;
                            dAlphaz_dVg = T6 * (1.0 - Alphaz * dVdsat_dVg);
                            dAlphaz_dVb = -T6 * (dVth_dVb + Alphaz * dVdsat_dVb);
                            T7 = T9 * T7;
                            T8 = T9 * T8;
                            T9 = 2.0 * T4 * (1.0 - 3.0 * T5);
                            Cdgb = (T7 * dAlphaz_dVg - T9
                                    * dVdsat_dVg) * dVgs_eff_dVg;
                            T12 = T7 * dAlphaz_dVb - T9 * dVdsat_dVb;
                            Cddb = T4 * (3.0 - 6.0 * T2 - 3.0 * T5);
                            Cdsb = -(Cdgb + T12
                                          + Cddb);

                            T9 = 2.0 * T4 * (1.0 + T5);
                            T10 = (T8 * dAlphaz_dVg - T9 * dVdsat_dVg)
                                  * dVgs_eff_dVg;
                            T11 = T8 * dAlphaz_dVb - T9 * dVdsat_dVb;
                            T12 = T4 * (2.0 * T2 + T5 - 1.0);
                            T0 = -(T10 + T11 + T12);

                            Cbgb = -(Cggb
                                     + Cdgb + T10);
                            Cbdb = -(Cgdb
                                     + Cddb + T12);
                            Cbsb = -(Cgsb
                                     + Cdsb + T0);
                            Qinv = -(qgate + qbulk);
                        }
                    }
                    else if (ModelParameters.Xpart < 0.5)
                    {
                        /* 40/60 Charge partition model */
                        if (Vds >= Vdsat)
                        {
                            /* saturation region */
                            T1 = Vdsat / 3.0;
                            qgate = CoxWL * (Vgs_eff - Vfb
                                                     - pParam.BSIM3phi - T1);
                            T2 = -Two_Third_CoxWL * Vgst;
                            qbulk = -(qgate + T2);
                            qdrn = 0.4 * T2;

                            Cggb = One_Third_CoxWL * (3.0
                                                      - dVdsat_dVg) * dVgs_eff_dVg;
                            T2 = -One_Third_CoxWL * dVdsat_dVb;
                            Cgsb = -(Cggb + T2);
                            Cgdb = 0.0;

                            T3 = 0.4 * Two_Third_CoxWL;
                            Cdgb = -T3 * dVgs_eff_dVg;
                            Cddb = 0.0;
                            T4 = T3 * dVth_dVb;
                            Cdsb = -(T4 + Cdgb);

                            Cbgb = -(Cggb
                                     - Two_Third_CoxWL * dVgs_eff_dVg);
                            T3 = -(T2 + Two_Third_CoxWL * dVth_dVb);
                            Cbsb = -(Cbgb + T3);
                            Cbdb = 0.0;
                            Qinv = -(qgate + qbulk);
                        }
                        else
                        {
                            /* linear region  */
                            Alphaz = Vgst / Vdsat;
                            T1 = 2.0 * Vdsat - Vds;
                            T2 = Vds / (3.0 * T1);
                            T3 = T2 * Vds;
                            T9 = 0.25 * CoxWL;
                            T4 = T9 * Alphaz;
                            qgate = CoxWL * (Vgs_eff - Vfb - pParam.BSIM3phi
                                             - 0.5 * (Vds - T3));

                            T5 = T3 / T1;
                            Cggb = CoxWL * (1.0 - T5 * dVdsat_dVg)
                                         * dVgs_eff_dVg;
                            tmp = -CoxWL * T5 * dVdsat_dVb;
                            Cgdb = CoxWL * (T2 - 0.5 + 0.5 * T5);
                            Cgsb = -(Cggb
                                     + Cgdb + tmp);

                            T6 = 1.0 / Vdsat;
                            dAlphaz_dVg = T6 * (1.0 - Alphaz * dVdsat_dVg);
                            dAlphaz_dVb = -T6 * (dVth_dVb + Alphaz * dVdsat_dVb);

                            T6 = 8.0 * Vdsat * Vdsat - 6.0 * Vdsat * Vds
                                 + 1.2 * Vds * Vds;
                            T8 = T2 / T1;
                            T7 = Vds - T1 - T8 * T6;
                            qdrn = T4 * T7;
                            T7 *= T9;
                            tmp = T8 / T1;
                            tmp1 = T4 * (2.0 - 4.0 * tmp * T6
                                         + T8 * (16.0 * Vdsat - 6.0 * Vds));

                            Cdgb = (T7 * dAlphaz_dVg - tmp1
                                    * dVdsat_dVg) * dVgs_eff_dVg;
                            T10 = T7 * dAlphaz_dVb - tmp1 * dVdsat_dVb;
                            Cddb = T4 * (2.0 - (1.0 / (3.0 * T1
                                                           * T1) + 2.0 * tmp) * T6 + T8
                                         * (6.0 * Vdsat - 2.4 * Vds));
                            Cdsb = -(Cdgb
                                     + T10 + Cddb);

                            T7 = 2.0 * (T1 + T3);
                            qbulk = -(qgate - T4 * T7);
                            T7 *= T9;
                            T0 = 4.0 * T4 * (1.0 - T5);
                            T12 = (-T7 * dAlphaz_dVg - Cdgb
                                                     - T0 * dVdsat_dVg) * dVgs_eff_dVg;
                            T11 = -T7 * dAlphaz_dVb - T10 - T0 * dVdsat_dVb;
                            T10 = -4.0 * T4 * (T2 - 0.5 + 0.5 * T5)
                                  - Cddb;
                            tmp = -(T10 + T11 + T12);

                            Cbgb = -(Cggb
                                     + Cdgb + T12);
                            Cbdb = -(Cgdb
                                     + Cddb + T10); /* bug fix */
                            Cbsb = -(Cgsb
                                     + Cdsb + tmp);
                            Qinv = -(qgate + qbulk);
                        }
                    }
                    else
                    {
                        /* 50/50 partitioning */
                        if (Vds >= Vdsat)
                        {
                            /* saturation region */
                            T1 = Vdsat / 3.0;
                            qgate = CoxWL * (Vgs_eff - Vfb
                                                     - pParam.BSIM3phi - T1);
                            T2 = -Two_Third_CoxWL * Vgst;
                            qbulk = -(qgate + T2);
                            qdrn = 0.5 * T2;

                            Cggb = One_Third_CoxWL * (3.0
                                                      - dVdsat_dVg) * dVgs_eff_dVg;
                            T2 = -One_Third_CoxWL * dVdsat_dVb;
                            Cgsb = -(Cggb + T2);
                            Cgdb = 0.0;

                            Cdgb = -One_Third_CoxWL * dVgs_eff_dVg;
                            Cddb = 0.0;
                            T4 = One_Third_CoxWL * dVth_dVb;
                            Cdsb = -(T4 + Cdgb);

                            Cbgb = -(Cggb
                                     - Two_Third_CoxWL * dVgs_eff_dVg);
                            T3 = -(T2 + Two_Third_CoxWL * dVth_dVb);
                            Cbsb = -(Cbgb + T3);
                            Cbdb = 0.0;
                            Qinv = -(qgate + qbulk);
                        }
                        else
                        {
                            /* linear region */
                            Alphaz = Vgst / Vdsat;
                            T1 = 2.0 * Vdsat - Vds;
                            T2 = Vds / (3.0 * T1);
                            T3 = T2 * Vds;
                            T9 = 0.25 * CoxWL;
                            T4 = T9 * Alphaz;
                            qgate = CoxWL * (Vgs_eff - Vfb - pParam.BSIM3phi
                                             - 0.5 * (Vds - T3));

                            T5 = T3 / T1;
                            Cggb = CoxWL * (1.0 - T5 * dVdsat_dVg)
                                         * dVgs_eff_dVg;
                            tmp = -CoxWL * T5 * dVdsat_dVb;
                            Cgdb = CoxWL * (T2 - 0.5 + 0.5 * T5);
                            Cgsb = -(Cggb
                                     + Cgdb + tmp);

                            T6 = 1.0 / Vdsat;
                            dAlphaz_dVg = T6 * (1.0 - Alphaz * dVdsat_dVg);
                            dAlphaz_dVb = -T6 * (dVth_dVb + Alphaz * dVdsat_dVb);

                            T7 = T1 + T3;
                            qdrn = -T4 * T7;
                            qbulk = -(qgate + qdrn + qdrn);
                            T7 *= T9;
                            T0 = T4 * (2.0 * T5 - 2.0);

                            Cdgb = (T0 * dVdsat_dVg - T7
                                    * dAlphaz_dVg) * dVgs_eff_dVg;
                            T12 = T0 * dVdsat_dVb - T7 * dAlphaz_dVb;
                            Cddb = T4 * (1.0 - 2.0 * T2 - T5);
                            Cdsb = -(Cdgb + T12
                                          + Cddb);

                            Cbgb = -(Cggb
                                     + 2.0 * Cdgb);
                            Cbdb = -(Cgdb
                                     + 2.0 * Cddb);
                            Cbsb = -(Cgsb
                                     + 2.0 * Cdsb);
                            Qinv = -(qgate + qbulk);
                        }
                    }
                }
            }
            else
            {
                if (Vbseff < 0.0)
                {
                    VbseffCV = Vbseff;
                    dVbseffCV_dVb = 1.0;
                }
                else
                {
                    VbseffCV = pParam.BSIM3phi - Phis;
                    dVbseffCV_dVb = -dPhis_dVb;
                }

                CoxWL = ModelParameters.Cox * pParam.BSIM3weffCV
                                 * pParam.BSIM3leffCV;

                /* Seperate VgsteffCV with noff and voffcv */
                noff = n * pParam.BSIM3noff;
                dnoff_dVd = pParam.BSIM3noff * dn_dVd;
                dnoff_dVb = pParam.BSIM3noff * dn_dVb;
                T0 = Vtm * noff;
                voffcv = pParam.BSIM3voffcv;
                VgstNVt = (Vgst - voffcv) / T0;

                if (VgstNVt > EXP_THRESHOLD)
                {
                    Vgsteff = Vgst - voffcv;
                    dVgsteff_dVg = dVgs_eff_dVg;
                    dVgsteff_dVd = -dVth_dVd;
                    dVgsteff_dVb = -dVth_dVb;
                }
                else if (VgstNVt < -EXP_THRESHOLD)
                {
                    Vgsteff = T0 * Math.Log(1.0 + MIN_EXP);
                    dVgsteff_dVg = 0.0;
                    dVgsteff_dVd = Vgsteff / noff;
                    dVgsteff_dVb = dVgsteff_dVd * dnoff_dVb;
                    dVgsteff_dVd *= dnoff_dVd;
                }
                else
                {
                    ExpVgst = Math.Exp(VgstNVt);
                    Vgsteff = T0 * Math.Log(1.0 + ExpVgst);
                    dVgsteff_dVg = ExpVgst / (1.0 + ExpVgst);
                    dVgsteff_dVd = -dVgsteff_dVg * (dVth_dVd + (Vgst - voffcv)
                                                    / noff * dnoff_dVd) + Vgsteff / noff * dnoff_dVd;
                    dVgsteff_dVb = -dVgsteff_dVg * (dVth_dVb + (Vgst - voffcv)
                                                    / noff * dnoff_dVb) + Vgsteff / noff * dnoff_dVb;
                    dVgsteff_dVg *= dVgs_eff_dVg;
                } /* End of VgsteffCV */

                if (ModelParameters.CapMod == 1)
                {
                    Vfb = pParam.BSIM3vfbzb;
                    Arg1 = Vgs_eff - VbseffCV - Vfb - Vgsteff;

                    if (Arg1 <= 0.0)
                    {
                        qgate = CoxWL * Arg1;
                        Cgg = CoxWL * (dVgs_eff_dVg - dVgsteff_dVg);
                        Cgd = -CoxWL * dVgsteff_dVd;
                        Cgb = -CoxWL * (dVbseffCV_dVb + dVgsteff_dVb);
                    }
                    else
                    {
                        T0 = 0.5 * pParam.BSIM3k1ox;
                        T1 = Math.Sqrt(T0 * T0 + Arg1);
                        T2 = CoxWL * T0 / T1;

                        qgate = CoxWL * pParam.BSIM3k1ox * (T1 - T0);

                        Cgg = T2 * (dVgs_eff_dVg - dVgsteff_dVg);
                        Cgd = -T2 * dVgsteff_dVd;
                        Cgb = -T2 * (dVbseffCV_dVb + dVgsteff_dVb);
                    }

                    qbulk = -qgate;
                    Cbg = -Cgg;
                    Cbd = -Cgd;
                    Cbb = -Cgb;

                    One_Third_CoxWL = CoxWL / 3.0;
                    Two_Third_CoxWL = 2.0 * One_Third_CoxWL;
                    AbulkCV = Abulk0 * pParam.BSIM3abulkCVfactor;
                    dAbulkCV_dVb = pParam.BSIM3abulkCVfactor * dAbulk0_dVb;
                    VdsatCV = Vgsteff / AbulkCV;
                    if (VdsatCV < Vds)
                    {
                        dVdsatCV_dVg = 1.0 / AbulkCV;
                        dVdsatCV_dVb = -VdsatCV * dAbulkCV_dVb / AbulkCV;
                        T0 = Vgsteff - VdsatCV / 3.0;
                        dT0_dVg = 1.0 - dVdsatCV_dVg / 3.0;
                        dT0_dVb = -dVdsatCV_dVb / 3.0;
                        qgate += CoxWL * T0;
                        Cgg1 = CoxWL * dT0_dVg;
                        Cgb1 = CoxWL * dT0_dVb + Cgg1 * dVgsteff_dVb;
                        Cgd1 = Cgg1 * dVgsteff_dVd;
                        Cgg1 *= dVgsteff_dVg;
                        Cgg += Cgg1;
                        Cgb += Cgb1;
                        Cgd += Cgd1;

                        T0 = VdsatCV - Vgsteff;
                        dT0_dVg = dVdsatCV_dVg - 1.0;
                        dT0_dVb = dVdsatCV_dVb;
                        qbulk += One_Third_CoxWL * T0;
                        Cbg1 = One_Third_CoxWL * dT0_dVg;
                        Cbb1 = One_Third_CoxWL * dT0_dVb + Cbg1 * dVgsteff_dVb;
                        Cbd1 = Cbg1 * dVgsteff_dVd;
                        Cbg1 *= dVgsteff_dVg;
                        Cbg += Cbg1;
                        Cbb += Cbb1;
                        Cbd += Cbd1;

                        if (ModelParameters.Xpart > 0.5)
                            T0 = -Two_Third_CoxWL;
                        else if (ModelParameters.Xpart < 0.5)
                            T0 = -0.4 * CoxWL;
                        else
                            T0 = -One_Third_CoxWL;

                        qsrc = T0 * Vgsteff;
                        Csg = T0 * dVgsteff_dVg;
                        Csb = T0 * dVgsteff_dVb;
                        Csd = T0 * dVgsteff_dVd;
                        Cgb *= dVbseff_dVb;
                        Cbb *= dVbseff_dVb;
                        Csb *= dVbseff_dVb;
                    }
                    else
                    {
                        T0 = AbulkCV * Vds;
                        T1 = 12.0 * (Vgsteff - 0.5 * T0 + 1.0e-20);
                        T2 = Vds / T1;
                        T3 = T0 * T2;
                        dT3_dVg = -12.0 * T2 * T2 * AbulkCV;
                        dT3_dVd = 6.0 * T0 * (4.0 * Vgsteff - T0) / T1 / T1 - 0.5;
                        dT3_dVb = 12.0 * T2 * T2 * dAbulkCV_dVb * Vgsteff;

                        qgate += CoxWL * (Vgsteff - 0.5 * Vds + T3);
                        Cgg1 = CoxWL * (1.0 + dT3_dVg);
                        Cgb1 = CoxWL * dT3_dVb + Cgg1 * dVgsteff_dVb;
                        Cgd1 = CoxWL * dT3_dVd + Cgg1 * dVgsteff_dVd;
                        Cgg1 *= dVgsteff_dVg;
                        Cgg += Cgg1;
                        Cgb += Cgb1;
                        Cgd += Cgd1;

                        qbulk += CoxWL * (1.0 - AbulkCV) * (0.5 * Vds - T3);
                        Cbg1 = -CoxWL * ((1.0 - AbulkCV) * dT3_dVg);
                        Cbb1 = -CoxWL * ((1.0 - AbulkCV) * dT3_dVb
                                         + (0.5 * Vds - T3) * dAbulkCV_dVb)
                               + Cbg1 * dVgsteff_dVb;
                        Cbd1 = -CoxWL * (1.0 - AbulkCV) * dT3_dVd
                               + Cbg1 * dVgsteff_dVd;
                        Cbg1 *= dVgsteff_dVg;
                        Cbg += Cbg1;
                        Cbb += Cbb1;
                        Cbd += Cbd1;

                        if (ModelParameters.Xpart > 0.5)
                        {
                            /* 0/100 Charge petition model */
                            T1 = T1 + T1;
                            qsrc = -CoxWL * (0.5 * Vgsteff + 0.25 * T0
                                             - T0 * T0 / T1);
                            Csg = -CoxWL * (0.5 + 24.0 * T0 * Vds / T1 / T1
                                            * AbulkCV);
                            Csb = -CoxWL * (0.25 * Vds * dAbulkCV_dVb
                                            - 12.0 * T0 * Vds / T1 / T1 * (4.0 * Vgsteff - T0)
                                                                        * dAbulkCV_dVb) + Csg * dVgsteff_dVb;
                            Csd = -CoxWL * (0.25 * AbulkCV - 12.0 * AbulkCV * T0
                                            / T1 / T1 * (4.0 * Vgsteff - T0))
                                  + Csg * dVgsteff_dVd;
                            Csg *= dVgsteff_dVg;
                        }
                        else if (ModelParameters.Xpart < 0.5)
                        {
                            /* 40/60 Charge petition model */
                            T1 = T1 / 12.0;
                            T2 = 0.5 * CoxWL / (T1 * T1);
                            T3 = Vgsteff * (2.0 * T0 * T0 / 3.0 + Vgsteff
                                            * (Vgsteff - 4.0 * T0 / 3.0))
                                 - 2.0 * T0 * T0 * T0 / 15.0;
                            qsrc = -T2 * T3;
                            T4 = 4.0 / 3.0 * Vgsteff * (Vgsteff - T0)
                                 + 0.4 * T0 * T0;
                            Csg = -2.0 * qsrc / T1 - T2 * (Vgsteff * (3.0
                                                                      * Vgsteff - 8.0 * T0 / 3.0)
                                                           + 2.0 * T0 * T0 / 3.0);
                            Csb = (qsrc / T1 * Vds + T2 * T4 * Vds) * dAbulkCV_dVb
                                  + Csg * dVgsteff_dVb;
                            Csd = (qsrc / T1 + T2 * T4) * AbulkCV
                                  + Csg * dVgsteff_dVd;
                            Csg *= dVgsteff_dVg;
                        }
                        else
                        {
                            /* 50/50 Charge petition model */
                            qsrc = -0.5 * (qgate + qbulk);
                            Csg = -0.5 * (Cgg1 + Cbg1);
                            Csb = -0.5 * (Cgb1 + Cbb1);
                            Csd = -0.5 * (Cgd1 + Cbd1);
                        }

                        Cgb *= dVbseff_dVb;
                        Cbb *= dVbseff_dVb;
                        Csb *= dVbseff_dVb;
                    }

                    qdrn = -(qgate + qbulk + qsrc);
                    Cggb = Cgg;
                    Cgsb = -(Cgg + Cgd + Cgb);
                    Cgdb = Cgd;
                    Cdgb = -(Cgg + Cbg + Csg);
                    Cdsb = Cgg + Cgd + Cgb + Cbg + Cbd + Cbb
                           + Csg + Csd + Csb;
                    Cddb = -(Cgd + Cbd + Csd);
                    Cbgb = Cbg;
                    Cbsb = -(Cbg + Cbd + Cbb);
                    Cbdb = Cbd;
                    Qinv = -(qgate + qbulk);
                }

                else if (ModelParameters.CapMod == 2)
                {
                    Vfb = pParam.BSIM3vfbzb;
                    V3 = Vfb - Vgs_eff + VbseffCV - DELTA_3;
                    if (Vfb <= 0.0)
                    {
                        T0 = Math.Sqrt(V3 * V3 - 4.0 * DELTA_3 * Vfb);
                        T2 = -DELTA_3 / T0;
                    }
                    else
                    {
                        T0 = Math.Sqrt(V3 * V3 + 4.0 * DELTA_3 * Vfb);
                        T2 = DELTA_3 / T0;
                    }

                    T1 = 0.5 * (1.0 + V3 / T0);
                    Vfbeff = Vfb - 0.5 * (V3 + T0);
                    dVfbeff_dVg = T1 * dVgs_eff_dVg;
                    dVfbeff_dVb = -T1 * dVbseffCV_dVb;
                    Qac0 = CoxWL * (Vfbeff - Vfb);
                    dQac0_dVg = CoxWL * dVfbeff_dVg;
                    dQac0_dVb = CoxWL * dVfbeff_dVb;

                    T0 = 0.5 * pParam.BSIM3k1ox;
                    T3 = Vgs_eff - Vfbeff - VbseffCV - Vgsteff;
                    if (pParam.BSIM3k1ox == 0.0)
                    {
                        T1 = 0.0;
                        T2 = 0.0;
                    }
                    else if (T3 < 0.0)
                    {
                        T1 = T0 + T3 / pParam.BSIM3k1ox;
                        T2 = CoxWL;
                    }
                    else
                    {
                        T1 = Math.Sqrt(T0 * T0 + T3);
                        T2 = CoxWL * T0 / T1;
                    }

                    Qsub0 = CoxWL * pParam.BSIM3k1ox * (T1 - T0);

                    dQsub0_dVg = T2 * (dVgs_eff_dVg - dVfbeff_dVg - dVgsteff_dVg);
                    dQsub0_dVd = -T2 * dVgsteff_dVd;
                    dQsub0_dVb = -T2 * (dVfbeff_dVb + dVbseffCV_dVb
                                                    + dVgsteff_dVb);

                    AbulkCV = Abulk0 * pParam.BSIM3abulkCVfactor;
                    dAbulkCV_dVb = pParam.BSIM3abulkCVfactor * dAbulk0_dVb;
                    VdsatCV = Vgsteff / AbulkCV;

                    V4 = VdsatCV - Vds - DELTA_4;
                    T0 = Math.Sqrt(V4 * V4 + 4.0 * DELTA_4 * VdsatCV);
                    VdseffCV = VdsatCV - 0.5 * (V4 + T0);
                    T1 = 0.5 * (1.0 + V4 / T0);
                    T2 = DELTA_4 / T0;
                    T3 = (1.0 - T1 - T2) / AbulkCV;
                    dVdseffCV_dVg = T3;
                    dVdseffCV_dVd = T1;
                    dVdseffCV_dVb = -T3 * VdsatCV * dAbulkCV_dVb;
                    /* Added to eliminate non-zero VdseffCV at Vds=0.0 */
                    if (Vds == 0.0)
                    {
                        VdseffCV = 0.0;
                        dVdseffCV_dVg = 0.0;
                        dVdseffCV_dVb = 0.0;
                    }

                    T0 = AbulkCV * VdseffCV;
                    T1 = 12.0 * (Vgsteff - 0.5 * T0 + 1e-20);
                    T2 = VdseffCV / T1;
                    T3 = T0 * T2;

                    T4 = 1.0 - 12.0 * T2 * T2 * AbulkCV;
                    T5 = 6.0 * T0 * (4.0 * Vgsteff - T0) / (T1 * T1) - 0.5;
                    T6 = 12.0 * T2 * T2 * Vgsteff;

                    qinoi = -CoxWL * (Vgsteff - 0.5 * T0 + AbulkCV * T3);
                    qgate = CoxWL * (Vgsteff - 0.5 * VdseffCV + T3);
                    Cgg1 = CoxWL * (T4 + T5 * dVdseffCV_dVg);
                    Cgd1 = CoxWL * T5 * dVdseffCV_dVd + Cgg1 * dVgsteff_dVd;
                    Cgb1 = CoxWL * (T5 * dVdseffCV_dVb + T6 * dAbulkCV_dVb)
                           + Cgg1 * dVgsteff_dVb;
                    Cgg1 *= dVgsteff_dVg;

                    T7 = 1.0 - AbulkCV;
                    qbulk = CoxWL * T7 * (0.5 * VdseffCV - T3);
                    T4 = -T7 * (T4 - 1.0);
                    T5 = -T7 * T5;
                    T6 = -(T7 * T6 + (0.5 * VdseffCV - T3));
                    Cbg1 = CoxWL * (T4 + T5 * dVdseffCV_dVg);
                    Cbd1 = CoxWL * T5 * dVdseffCV_dVd + Cbg1 * dVgsteff_dVd;
                    Cbb1 = CoxWL * (T5 * dVdseffCV_dVb + T6 * dAbulkCV_dVb)
                           + Cbg1 * dVgsteff_dVb;
                    Cbg1 *= dVgsteff_dVg;

                    if (ModelParameters.Xpart > 0.5)
                    {
                        /* 0/100 Charge petition model */
                        T1 = T1 + T1;
                        qsrc = -CoxWL * (0.5 * Vgsteff + 0.25 * T0
                                         - T0 * T0 / T1);
                        T7 = (4.0 * Vgsteff - T0) / (T1 * T1);
                        T4 = -(0.5 + 24.0 * T0 * T0 / (T1 * T1));
                        T5 = -(0.25 * AbulkCV - 12.0 * AbulkCV * T0 * T7);
                        T6 = -(0.25 * VdseffCV - 12.0 * T0 * VdseffCV * T7);
                        Csg = CoxWL * (T4 + T5 * dVdseffCV_dVg);
                        Csd = CoxWL * T5 * dVdseffCV_dVd + Csg * dVgsteff_dVd;
                        Csb = CoxWL * (T5 * dVdseffCV_dVb + T6 * dAbulkCV_dVb)
                              + Csg * dVgsteff_dVb;
                        Csg *= dVgsteff_dVg;
                    }
                    else if (ModelParameters.Xpart < 0.5)
                    {
                        /* 40/60 Charge petition model */
                        T1 = T1 / 12.0;
                        T2 = 0.5 * CoxWL / (T1 * T1);
                        T3 = Vgsteff * (2.0 * T0 * T0 / 3.0 + Vgsteff
                                        * (Vgsteff - 4.0 * T0 / 3.0))
                             - 2.0 * T0 * T0 * T0 / 15.0;
                        qsrc = -T2 * T3;
                        T7 = 4.0 / 3.0 * Vgsteff * (Vgsteff - T0)
                             + 0.4 * T0 * T0;
                        T4 = -2.0 * qsrc / T1 - T2 * (Vgsteff * (3.0
                                                                 * Vgsteff - 8.0 * T0 / 3.0)
                                                      + 2.0 * T0 * T0 / 3.0);
                        T5 = (qsrc / T1 + T2 * T7) * AbulkCV;
                        T6 = qsrc / T1 * VdseffCV + T2 * T7 * VdseffCV;
                        Csg = T4 + T5 * dVdseffCV_dVg;
                        Csd = T5 * dVdseffCV_dVd + Csg * dVgsteff_dVd;
                        Csb = T5 * dVdseffCV_dVb + T6 * dAbulkCV_dVb
                                                 + Csg * dVgsteff_dVb;
                        Csg *= dVgsteff_dVg;
                    }
                    else
                    {
                        /* 50/50 Charge petition model */
                        qsrc = -0.5 * (qgate + qbulk);
                        Csg = -0.5 * (Cgg1 + Cbg1);
                        Csb = -0.5 * (Cgb1 + Cbb1);
                        Csd = -0.5 * (Cgd1 + Cbd1);
                    }

                    qgate += Qac0 + Qsub0;
                    qbulk -= Qac0 + Qsub0;
                    qdrn = -(qgate + qbulk + qsrc);

                    Cgg = dQac0_dVg + dQsub0_dVg + Cgg1;
                    Cgd = dQsub0_dVd + Cgd1;
                    Cgb = dQac0_dVb + dQsub0_dVb + Cgb1;

                    Cbg = Cbg1 - dQac0_dVg - dQsub0_dVg;
                    Cbd = Cbd1 - dQsub0_dVd;
                    Cbb = Cbb1 - dQac0_dVb - dQsub0_dVb;

                    Cgb *= dVbseff_dVb;
                    Cbb *= dVbseff_dVb;
                    Csb *= dVbseff_dVb;

                    Cggb = Cgg;
                    Cgsb = -(Cgg + Cgd + Cgb);
                    Cgdb = Cgd;
                    Cdgb = -(Cgg + Cbg + Csg);
                    Cdsb = Cgg + Cgd + Cgb + Cbg + Cbd + Cbb
                           + Csg + Csd + Csb;
                    Cddb = -(Cgd + Cbd + Csd);
                    Cbgb = Cbg;
                    Cbsb = -(Cbg + Cbd + Cbb);
                    Cbdb = Cbd;
                    Qinv = qinoi;
                }
                /* New Charge-Thickness capMod (CTM) begins */
                else if (ModelParameters.CapMod == 3)
                {
                    V3 = pParam.BSIM3vfbzb - Vgs_eff + VbseffCV - DELTA_3;
                    if (pParam.BSIM3vfbzb <= 0.0)
                    {
                        T0 = Math.Sqrt(V3 * V3 - 4.0 * DELTA_3 * pParam.BSIM3vfbzb);
                        T2 = -DELTA_3 / T0;
                    }
                    else
                    {
                        T0 = Math.Sqrt(V3 * V3 + 4.0 * DELTA_3 * pParam.BSIM3vfbzb);
                        T2 = DELTA_3 / T0;
                    }

                    T1 = 0.5 * (1.0 + V3 / T0);
                    Vfbeff = pParam.BSIM3vfbzb - 0.5 * (V3 + T0);
                    dVfbeff_dVg = T1 * dVgs_eff_dVg;
                    dVfbeff_dVb = -T1 * dVbseffCV_dVb;

                    Cox = ModelParameters.Cox;
                    Tox = 1.0e8 * ModelParameters.Tox;
                    T0 = (Vgs_eff - VbseffCV - pParam.BSIM3vfbzb) / Tox;
                    dT0_dVg = dVgs_eff_dVg / Tox;
                    dT0_dVb = -dVbseffCV_dVb / Tox;

                    tmp = T0 * pParam.BSIM3acde;
                    if (-EXP_THRESHOLD < tmp && tmp < EXP_THRESHOLD)
                    {
                        Tcen = pParam.BSIM3ldeb * Math.Exp(tmp);
                        dTcen_dVg = pParam.BSIM3acde * Tcen;
                        dTcen_dVb = dTcen_dVg * dT0_dVb;
                        dTcen_dVg *= dT0_dVg;
                    }
                    else if (tmp <= -EXP_THRESHOLD)
                    {
                        Tcen = pParam.BSIM3ldeb * MIN_EXP;
                        dTcen_dVg = dTcen_dVb = 0.0;
                    }
                    else
                    {
                        Tcen = pParam.BSIM3ldeb * MAX_EXP;
                        dTcen_dVg = dTcen_dVb = 0.0;
                    }

                    LINK = 1.0e-3 * ModelParameters.Tox;
                    V3 = pParam.BSIM3ldeb - Tcen - LINK;
                    V4 = Math.Sqrt(V3 * V3 + 4.0 * LINK * pParam.BSIM3ldeb);
                    Tcen = pParam.BSIM3ldeb - 0.5 * (V3 + V4);
                    T1 = 0.5 * (1.0 + V3 / V4);
                    dTcen_dVg *= T1;
                    dTcen_dVb *= T1;

                    Ccen = EPSSI / Tcen;
                    T2 = Cox / (Cox + Ccen);
                    Coxeff = T2 * Ccen;
                    T3 = -Ccen / Tcen;
                    dCoxeff_dVg = T2 * T2 * T3;
                    dCoxeff_dVb = dCoxeff_dVg * dTcen_dVb;
                    dCoxeff_dVg *= dTcen_dVg;
                    CoxWLcen = CoxWL * Coxeff / Cox;

                    Qac0 = CoxWLcen * (Vfbeff - pParam.BSIM3vfbzb);
                    QovCox = Qac0 / Coxeff;
                    dQac0_dVg = CoxWLcen * dVfbeff_dVg
                                + QovCox * dCoxeff_dVg;
                    dQac0_dVb = CoxWLcen * dVfbeff_dVb
                                + QovCox * dCoxeff_dVb;

                    T0 = 0.5 * pParam.BSIM3k1ox;
                    T3 = Vgs_eff - Vfbeff - VbseffCV - Vgsteff;
                    if (pParam.BSIM3k1ox == 0.0)
                    {
                        T1 = 0.0;
                        T2 = 0.0;
                    }
                    else if (T3 < 0.0)
                    {
                        T1 = T0 + T3 / pParam.BSIM3k1ox;
                        T2 = CoxWLcen;
                    }
                    else
                    {
                        T1 = Math.Sqrt(T0 * T0 + T3);
                        T2 = CoxWLcen * T0 / T1;
                    }

                    Qsub0 = CoxWLcen * pParam.BSIM3k1ox * (T1 - T0);
                    QovCox = Qsub0 / Coxeff;
                    dQsub0_dVg = T2 * (dVgs_eff_dVg - dVfbeff_dVg - dVgsteff_dVg)
                                 + QovCox * dCoxeff_dVg;
                    dQsub0_dVd = -T2 * dVgsteff_dVd;
                    dQsub0_dVb = -T2 * (dVfbeff_dVb + dVbseffCV_dVb + dVgsteff_dVb)
                                 + QovCox * dCoxeff_dVb;

                    /* Gate-bias dependent delta Phis begins */
                    if (pParam.BSIM3k1ox <= 0.0)
                    {
                        Denomi = 0.25 * pParam.BSIM3moin * Vtm;
                        T0 = 0.5 * pParam.BSIM3sqrtPhi;
                    }
                    else
                    {
                        Denomi = pParam.BSIM3moin * Vtm
                                                  * pParam.BSIM3k1ox * pParam.BSIM3k1ox;
                        T0 = pParam.BSIM3k1ox * pParam.BSIM3sqrtPhi;
                    }

                    T1 = 2.0 * T0 + Vgsteff;

                    DeltaPhi = Vtm * Math.Log(1.0 + T1 * Vgsteff / Denomi);
                    dDeltaPhi_dVg = 2.0 * Vtm * (T1 - T0) / (Denomi + T1 * Vgsteff);
                    dDeltaPhi_dVd = dDeltaPhi_dVg * dVgsteff_dVd;
                    dDeltaPhi_dVb = dDeltaPhi_dVg * dVgsteff_dVb;
                    /* End of delta Phis */

                    T3 = 4.0 * (Vth - pParam.BSIM3vfbzb - pParam.BSIM3phi);
                    Tox += Tox;
                    if (T3 >= 0.0)
                    {
                        T0 = (Vgsteff + T3) / Tox;
                        dT0_dVd = (dVgsteff_dVd + 4.0 * dVth_dVd) / Tox;
                        dT0_dVb = (dVgsteff_dVb + 4.0 * dVth_dVb) / Tox;
                    }
                    else
                    {
                        T0 = (Vgsteff + 1.0e-20) / Tox;
                        dT0_dVd = dVgsteff_dVd / Tox;
                        dT0_dVb = dVgsteff_dVb / Tox;
                    }

                    tmp = Math.Exp(0.7 * Math.Log(T0));
                    T1 = 1.0 + tmp;
                    T2 = 0.7 * tmp / (T0 * Tox);
                    Tcen = 1.9e-9 / T1;
                    dTcen_dVg = -1.9e-9 * T2 / T1 / T1;
                    dTcen_dVd = Tox * dTcen_dVg;
                    dTcen_dVb = dTcen_dVd * dT0_dVb;
                    dTcen_dVd *= dT0_dVd;
                    dTcen_dVg *= dVgsteff_dVg;

                    Ccen = EPSSI / Tcen;
                    T0 = Cox / (Cox + Ccen);
                    Coxeff = T0 * Ccen;
                    T1 = -Ccen / Tcen;
                    dCoxeff_dVg = T0 * T0 * T1;
                    dCoxeff_dVd = dCoxeff_dVg * dTcen_dVd;
                    dCoxeff_dVb = dCoxeff_dVg * dTcen_dVb;
                    dCoxeff_dVg *= dTcen_dVg;
                    CoxWLcen = CoxWL * Coxeff / Cox;

                    AbulkCV = Abulk0 * pParam.BSIM3abulkCVfactor;
                    dAbulkCV_dVb = pParam.BSIM3abulkCVfactor * dAbulk0_dVb;
                    VdsatCV = (Vgsteff - DeltaPhi) / AbulkCV;
                    V4 = VdsatCV - Vds - DELTA_4;
                    T0 = Math.Sqrt(V4 * V4 + 4.0 * DELTA_4 * VdsatCV);
                    VdseffCV = VdsatCV - 0.5 * (V4 + T0);
                    T1 = 0.5 * (1.0 + V4 / T0);
                    T2 = DELTA_4 / T0;
                    T3 = (1.0 - T1 - T2) / AbulkCV;
                    T4 = T3 * (1.0 - dDeltaPhi_dVg);
                    dVdseffCV_dVg = T4;
                    dVdseffCV_dVd = T1;
                    dVdseffCV_dVb = -T3 * VdsatCV * dAbulkCV_dVb;
                    /* Added to eliminate non-zero VdseffCV at Vds=0.0 */
                    if (Vds == 0.0)
                    {
                        VdseffCV = 0.0;
                        dVdseffCV_dVg = 0.0;
                        dVdseffCV_dVb = 0.0;
                    }

                    T0 = AbulkCV * VdseffCV;
                    T1 = Vgsteff - DeltaPhi;
                    T2 = 12.0 * (T1 - 0.5 * T0 + 1.0e-20);
                    T3 = T0 / T2;
                    T4 = 1.0 - 12.0 * T3 * T3;
                    T5 = AbulkCV * (6.0 * T0 * (4.0 * T1 - T0) / (T2 * T2) - 0.5);
                    T6 = T5 * VdseffCV / AbulkCV;

                    qgate = qinoi = CoxWLcen * (T1 - T0 * (0.5 - T3));
                    QovCox = qgate / Coxeff;
                    Cgg1 = CoxWLcen * (T4 * (1.0 - dDeltaPhi_dVg)
                                       + T5 * dVdseffCV_dVg);
                    Cgd1 = CoxWLcen * T5 * dVdseffCV_dVd + Cgg1
                           * dVgsteff_dVd + QovCox * dCoxeff_dVd;
                    Cgb1 = CoxWLcen * (T5 * dVdseffCV_dVb + T6 * dAbulkCV_dVb)
                           + Cgg1 * dVgsteff_dVb + QovCox * dCoxeff_dVb;
                    Cgg1 = Cgg1 * dVgsteff_dVg + QovCox * dCoxeff_dVg;


                    T7 = 1.0 - AbulkCV;
                    T8 = T2 * T2;
                    T9 = 12.0 * T7 * T0 * T0 / (T8 * AbulkCV);
                    T10 = T9 * (1.0 - dDeltaPhi_dVg);
                    T11 = -T7 * T5 / AbulkCV;
                    T12 = -(T9 * T1 / AbulkCV + VdseffCV * (0.5 - T0 / T2));

                    qbulk = CoxWLcen * T7 * (0.5 * VdseffCV - T0 * VdseffCV / T2);
                    QovCox = qbulk / Coxeff;
                    Cbg1 = CoxWLcen * (T10 + T11 * dVdseffCV_dVg);
                    Cbd1 = CoxWLcen * T11 * dVdseffCV_dVd + Cbg1
                           * dVgsteff_dVd + QovCox * dCoxeff_dVd;
                    Cbb1 = CoxWLcen * (T11 * dVdseffCV_dVb + T12 * dAbulkCV_dVb)
                           + Cbg1 * dVgsteff_dVb + QovCox * dCoxeff_dVb;
                    Cbg1 = Cbg1 * dVgsteff_dVg + QovCox * dCoxeff_dVg;

                    if (ModelParameters.Xpart > 0.5)
                    {
                        /* 0/100 partition */
                        qsrc = -CoxWLcen * (T1 / 2.0 + T0 / 4.0
                                            - 0.5 * T0 * T0 / T2);
                        QovCox = qsrc / Coxeff;
                        T2 += T2;
                        T3 = T2 * T2;
                        T7 = -(0.25 - 12.0 * T0 * (4.0 * T1 - T0) / T3);
                        T4 = -(0.5 + 24.0 * T0 * T0 / T3) * (1.0 - dDeltaPhi_dVg);
                        T5 = T7 * AbulkCV;
                        T6 = T7 * VdseffCV;

                        Csg = CoxWLcen * (T4 + T5 * dVdseffCV_dVg);
                        Csd = CoxWLcen * T5 * dVdseffCV_dVd + Csg * dVgsteff_dVd
                                                            + QovCox * dCoxeff_dVd;
                        Csb = CoxWLcen * (T5 * dVdseffCV_dVb + T6 * dAbulkCV_dVb)
                              + Csg * dVgsteff_dVb + QovCox * dCoxeff_dVb;
                        Csg = Csg * dVgsteff_dVg + QovCox * dCoxeff_dVg;
                    }
                    else if (ModelParameters.Xpart < 0.5)
                    {
                        /* 40/60 partition */
                        T2 = T2 / 12.0;
                        T3 = 0.5 * CoxWLcen / (T2 * T2);
                        T4 = T1 * (2.0 * T0 * T0 / 3.0 + T1 * (T1 - 4.0
                                                               * T0 / 3.0)) - 2.0 * T0 * T0 * T0 / 15.0;
                        qsrc = -T3 * T4;
                        QovCox = qsrc / Coxeff;
                        T8 = 4.0 / 3.0 * T1 * (T1 - T0) + 0.4 * T0 * T0;
                        T5 = -2.0 * qsrc / T2 - T3 * (T1 * (3.0 * T1 - 8.0
                                                            * T0 / 3.0) + 2.0 * T0 * T0 / 3.0);
                        T6 = AbulkCV * (qsrc / T2 + T3 * T8);
                        T7 = T6 * VdseffCV / AbulkCV;

                        Csg = T5 * (1.0 - dDeltaPhi_dVg) + T6 * dVdseffCV_dVg;
                        Csd = Csg * dVgsteff_dVd + T6 * dVdseffCV_dVd
                                                 + QovCox * dCoxeff_dVd;
                        Csb = Csg * dVgsteff_dVb + T6 * dVdseffCV_dVb
                                                 + T7 * dAbulkCV_dVb + QovCox * dCoxeff_dVb;
                        Csg = Csg * dVgsteff_dVg + QovCox * dCoxeff_dVg;
                    }
                    else
                    {
                        /* 50/50 partition */
                        qsrc = -0.5 * qgate;
                        Csg = -0.5 * Cgg1;
                        Csd = -0.5 * Cgd1;
                        Csb = -0.5 * Cgb1;
                    }

                    qgate += Qac0 + Qsub0 - qbulk;
                    qbulk -= Qac0 + Qsub0;
                    qdrn = -(qgate + qbulk + qsrc);

                    Cbg = Cbg1 - dQac0_dVg - dQsub0_dVg;
                    Cbd = Cbd1 - dQsub0_dVd;
                    Cbb = Cbb1 - dQac0_dVb - dQsub0_dVb;

                    Cgg = Cgg1 - Cbg;
                    Cgd = Cgd1 - Cbd;
                    Cgb = Cgb1 - Cbb;

                    Cgb *= dVbseff_dVb;
                    Cbb *= dVbseff_dVb;
                    Csb *= dVbseff_dVb;

                    Cggb = Cgg;
                    Cgsb = -(Cgg + Cgd + Cgb);
                    Cgdb = Cgd;
                    Cdgb = -(Cgg + Cbg + Csg);
                    Cdsb = Cgg + Cgd + Cgb + Cbg + Cbd + Cbb
                           + Csg + Csd + Csb;
                    Cddb = -(Cgd + Cbd + Csd);
                    Cbgb = Cbg;
                    Cbsb = -(Cbg + Cbd + Cbb);
                    Cbdb = Cbd;
                    Qinv = -qinoi;
                } /* End of CTM */
            }

            finished:
            /* Returning Values to Calling Routine */
            /*
             *  COMPUTE EQUIVALENT DRAIN CURRENT SOURCE
             */

            Qgate = qgate;
            Qbulk = qbulk;
            Qdrn = qdrn;
            Cd = cdrain;

            if (chargeComputationNeeded)
            {
                /*  charge storage elements
                            *  bulk-drain and bulk-source depletion capacitances
                            *  czbd : zero bias drain junction capacitance
                            *  czbs : zero bias source junction capacitance
                            *  czbdsw: zero bias drain junction sidewall capacitance
                           along field oxide
                            *  czbssw: zero bias source junction sidewall capacitance
                           along field oxide
                        *  czbdswg: zero bias drain junction sidewall capacitance
                            along gate side
                        *  czbsswg: zero bias source junction sidewall capacitance
                            along gate side
                            */

                czbd = ModelTemperature.UnitAreaTempJctCap * BaseParameters.DrainArea; /*bug fix */
                czbs = ModelTemperature.UnitAreaTempJctCap * BaseParameters.SourceArea;
                if (BaseParameters.DrainPerimeter < pParam.BSIM3weff)
                {
                    czbdswg = ModelTemperature.UnitLengthGateSidewallTempJctCap * BaseParameters.DrainPerimeter;
                    czbdsw = 0.0;
                }
                else
                {
                    czbdsw = ModelTemperature.UnitLengthSidewallTempJctCap * (BaseParameters.DrainPerimeter - pParam.BSIM3weff);
                    czbdswg = ModelTemperature.UnitLengthGateSidewallTempJctCap * pParam.BSIM3weff;
                }

                if (BaseParameters.SourcePerimeter < pParam.BSIM3weff)
                {
                    czbssw = 0.0;
                    czbsswg = ModelTemperature.UnitLengthGateSidewallTempJctCap
                              * BaseParameters.SourcePerimeter;
                }
                else
                {
                    czbssw = ModelTemperature.UnitLengthSidewallTempJctCap * (BaseParameters.SourcePerimeter - pParam.BSIM3weff);
                    czbsswg = ModelTemperature.UnitLengthGateSidewallTempJctCap * pParam.BSIM3weff;
                }

                MJ = ModelParameters.BulkJctBotGradingCoeff;
                MJSW = ModelParameters.BulkJctSideGradingCoeff;
                MJSWG = ModelParameters.BulkJctGateSideGradingCoeff;

                /* Source Bulk Junction */
                if (vbs == 0.0)
                {
                    Qbs = 0.0;
                    Capbs = czbs + czbssw + czbsswg;
                }
                else if (vbs < 0.0)
                {
                    if (czbs > 0.0)
                    {
                        arg = 1.0 - vbs / ModelTemperature.PhiB;
                        if (MJ == 0.5)
                            sarg = 1.0 / Math.Sqrt(arg);
                        else
                            sarg = Math.Exp(-MJ * Math.Log(arg));
                        Qbs = ModelTemperature.PhiB * czbs * (1.0 - arg * sarg) / (1.0 - MJ);
                        Capbs = czbs * sarg;
                    }
                    else
                    {
                        Qbs = 0.0;
                        Capbs = 0.0;
                    }

                    if (czbssw > 0.0)
                    {
                        arg = 1.0 - vbs / ModelTemperature.PhiBSW;
                        if (MJSW == 0.5)
                            sarg = 1.0 / Math.Sqrt(arg);
                        else
                            sarg = Math.Exp(-MJSW * Math.Log(arg));
                        Qbs += ModelTemperature.PhiBSW * czbssw * (1.0 - arg * sarg) / (1.0 - MJSW);
                        Capbs += czbssw * sarg;
                    }

                    if (czbsswg > 0.0)
                    {
                        arg = 1.0 - vbs / ModelTemperature.PhiBSWG;
                        if (MJSWG == 0.5)
                            sarg = 1.0 / Math.Sqrt(arg);
                        else
                            sarg = Math.Exp(-MJSWG * Math.Log(arg));
                        Qbs += ModelTemperature.PhiBSWG * czbsswg * (1.0 - arg * sarg) / (1.0 - MJSWG);
                        Capbs += czbsswg * sarg;
                    }
                }
                else
                {
                    T0 = czbs + czbssw + czbsswg;
                    T1 = vbs * (czbs * MJ / ModelTemperature.PhiB + czbssw * MJSW
                                / ModelTemperature.PhiBSW + czbsswg * MJSWG / ModelTemperature.PhiBSWG);
                    Qbs = vbs * (T0 + 0.5 * T1);
                    Capbs = T0 + T1;
                }

                /* Drain Bulk Junction */
                if (vbd == 0.0)
                {
                    Qbd = 0.0;
                    Capbd = czbd + czbdsw + czbdswg;
                }
                else if (vbd < 0.0)
                {
                    if (czbd > 0.0)
                    {
                        arg = 1.0 - vbd / ModelTemperature.PhiB;
                        if (MJ == 0.5)
                            sarg = 1.0 / Math.Sqrt(arg);
                        else
                            sarg = Math.Exp(-MJ * Math.Log(arg));
                        Qbd = ModelTemperature.PhiB * czbd * (1.0 - arg * sarg) / (1.0 - MJ);
                        Capbd = czbd * sarg;
                    }
                    else
                    {
                        Qbd = 0.0;
                        Capbd = 0.0;
                    }

                    if (czbdsw > 0.0)
                    {
                        arg = 1.0 - vbd / ModelTemperature.PhiBSW;
                        if (MJSW == 0.5)
                            sarg = 1.0 / Math.Sqrt(arg);
                        else
                            sarg = Math.Exp(-MJSW * Math.Log(arg));
                        Qbd += ModelTemperature.PhiBSW * czbdsw * (1.0 - arg * sarg) / (1.0 - MJSW);
                        Capbd += czbdsw * sarg;
                    }

                    if (czbdswg > 0.0)
                    {
                        arg = 1.0 - vbd / ModelTemperature.PhiBSWG;
                        if (MJSWG == 0.5)
                            sarg = 1.0 / Math.Sqrt(arg);
                        else
                            sarg = Math.Exp(-MJSWG * Math.Log(arg));
                        Qbd += ModelTemperature.PhiBSWG * czbdswg * (1.0 - arg * sarg) / (1.0 - MJSWG);
                        Capbd += czbdswg * sarg;
                    }
                }
                else
                {
                    T0 = czbd + czbdsw + czbdswg;
                    T1 = vbd * (czbd * MJ / ModelTemperature.PhiB + czbdsw * MJSW / ModelTemperature.PhiBSW +
                                czbdswg * MJSWG / ModelTemperature.PhiBSWG);
                    Qbd = vbd * (T0 + 0.5 * T1);
                    Capbd = T0 + T1;
                }
            }

            /*
             *  check convergence
             */
            if (!BaseParameters.Off || state.Init == InitializationModes.Fix)
            {
                if (check)
                {
                    state.IsConvergent = false;
                }
            }

            this.Vbs = vbs;
            Vbd = vbd;
            this.Vgs = vgs;
            this.Vds = vds;
            Qdef = qdef;

            /* bulk and channel charge plus overlaps */

            if (!chargeComputationNeeded)
                goto line850;

            // line755:
            /* NQS begins */
            if (BaseParameters.NqsMod > 0)
            {
                qcheq = -(qbulk + qgate);

                Cqgb = -(Cggb + Cbgb);
                Cqdb = -(Cgdb + Cbdb);
                Cqsb = -(Cgsb + Cbsb);
                Cqbb = -(Cqgb + Cqdb + Cqsb);

                gtau_drift = Math.Abs(pParam.BSIM3tconst * qcheq) * ScalingFactor;
                T0 = pParam.BSIM3leffCV * pParam.BSIM3leffCV;
                gtau_diff = 16.0 * pParam.BSIM3u0temp * ModelTemperature.Vtm / T0
                            * ScalingFactor;
                Gtau = gtau_drift + gtau_diff;
            }

            if (ModelParameters.CapMod == 0) /* code merge -JX */
            {
                cgdo = pParam.BSIM3cgdo;
                qgdo = pParam.BSIM3cgdo * vgd;
                cgso = pParam.BSIM3cgso;
                qgso = pParam.BSIM3cgso * vgs;
            }
            else if (ModelParameters.CapMod == 1)
            {
                if (vgd < 0.0)
                {
                    T1 = Math.Sqrt(1.0 - 4.0 * vgd / pParam.BSIM3ckappa);
                    cgdo = pParam.BSIM3cgdo + pParam.BSIM3weffCV
                           * pParam.BSIM3cgdl / T1;
                    qgdo = pParam.BSIM3cgdo * vgd - pParam.BSIM3weffCV * 0.5
                                                                       * pParam.BSIM3cgdl * pParam.BSIM3ckappa *
                                                                       (T1 - 1.0);
                }
                else
                {
                    cgdo = pParam.BSIM3cgdo + pParam.BSIM3weffCV
                           * pParam.BSIM3cgdl;
                    qgdo = (pParam.BSIM3weffCV * pParam.BSIM3cgdl
                            + pParam.BSIM3cgdo) * vgd;
                }

                if (vgs < 0.0)
                {
                    T1 = Math.Sqrt(1.0 - 4.0 * vgs / pParam.BSIM3ckappa);
                    cgso = pParam.BSIM3cgso + pParam.BSIM3weffCV
                           * pParam.BSIM3cgsl / T1;
                    qgso = pParam.BSIM3cgso * vgs - pParam.BSIM3weffCV * 0.5
                                                                       * pParam.BSIM3cgsl * pParam.BSIM3ckappa *
                                                                       (T1 - 1.0);
                }
                else
                {
                    cgso = pParam.BSIM3cgso + pParam.BSIM3weffCV
                           * pParam.BSIM3cgsl;
                    qgso = (pParam.BSIM3weffCV * pParam.BSIM3cgsl
                            + pParam.BSIM3cgso) * vgs;
                }
            }
            else
            {
                T0 = vgd + DELTA_1;
                T1 = Math.Sqrt(T0 * T0 + 4.0 * DELTA_1);
                T2 = 0.5 * (T0 - T1);

                T3 = pParam.BSIM3weffCV * pParam.BSIM3cgdl;
                T4 = Math.Sqrt(1.0 - 4.0 * T2 / pParam.BSIM3ckappa);
                cgdo = pParam.BSIM3cgdo + T3 - T3 * (1.0 - 1.0 / T4)
                                                  * (0.5 - 0.5 * T0 / T1);
                qgdo = (pParam.BSIM3cgdo + T3) * vgd - T3 * (T2
                                                             + 0.5 * pParam.BSIM3ckappa * (T4 - 1.0));

                T0 = vgs + DELTA_1;
                T1 = Math.Sqrt(T0 * T0 + 4.0 * DELTA_1);
                T2 = 0.5 * (T0 - T1);
                T3 = pParam.BSIM3weffCV * pParam.BSIM3cgsl;
                T4 = Math.Sqrt(1.0 - 4.0 * T2 / pParam.BSIM3ckappa);
                cgso = pParam.BSIM3cgso + T3 - T3 * (1.0 - 1.0 / T4)
                                                  * (0.5 - 0.5 * T0 / T1);
                qgso = (pParam.BSIM3cgso + T3) * vgs - T3 * (T2
                                                             + 0.5 * pParam.BSIM3ckappa * (T4 - 1.0));
            }

            base.Cgdo = cgdo;
            base.Cgso = cgso;

            ag0 = TranBehavior?.Qg.Jacobian(1.0) ?? 0.0;
            if (Mode > 0)
            {
                if (BaseParameters.NqsMod == 0)
                {
                    gcggb = (Cggb + cgdo + cgso
                             + pParam.BSIM3cgbo) * ag0;
                    gcgdb = (Cgdb - cgdo) * ag0;
                    gcgsb = (Cgsb - cgso) * ag0;

                    gcdgb = (Cdgb - cgdo) * ag0;
                    gcddb = (Cddb + Capbd + cgdo) * ag0;
                    gcdsb = Cdsb * ag0;

                    gcsgb = -(Cggb + Cbgb + Cdgb + cgso) * ag0;
                    gcsdb = -(Cgdb + Cbdb + Cddb) * ag0;
                    gcssb = (Capbs + cgso - (Cgsb + Cbsb + Cdsb)) * ag0;

                    gcbgb = (Cbgb - pParam.BSIM3cgbo) * ag0;
                    gcbdb = (Cbdb - Capbd) * ag0;
                    gcbsb = (Cbsb - Capbs) * ag0;

                    qgd = qgdo;
                    qgs = qgso;
                    qgb = pParam.BSIM3cgbo * vgb;
                    qgate += qgd + qgs + qgb;
                    qbulk -= qgb;
                    qdrn -= qgd;
                    qsrc = -(qgate + qbulk + qdrn);

                    ggtg = ggtd = ggtb = ggts = 0.0;
                    sxpart = 0.6;
                    dxpart = 0.4;
                    ddxpart_dVd = ddxpart_dVg = ddxpart_dVb = ddxpart_dVs = 0.0;
                    dsxpart_dVd = dsxpart_dVg = dsxpart_dVb = dsxpart_dVs = 0.0;
                }
                else
                {
                    if (qcheq > 0.0)
                        T0 = pParam.BSIM3tconst * qdef * ScalingFactor;
                    else
                        T0 = -pParam.BSIM3tconst * qdef * ScalingFactor;
                    ggtg = Gtg = T0 * Cqgb;
                    ggtd = Gtd = T0 * Cqdb;
                    ggts = Gts = T0 * Cqsb;
                    ggtb = Gtb = T0 * Cqbb;
                    gqdef = ScalingFactor * ag0;

                    gcqgb = Cqgb * ag0;
                    gcqdb = Cqdb * ag0;
                    gcqsb = Cqsb * ag0;
                    gcqbb = Cqbb * ag0;

                    gcggb = (cgdo + cgso + pParam.BSIM3cgbo) * ag0;
                    gcgdb = -cgdo * ag0;
                    gcgsb = -cgso * ag0;

                    gcdgb = -cgdo * ag0;
                    gcddb = (Capbd + cgdo) * ag0;
                    gcdsb = 0.0;

                    gcsgb = -cgso * ag0;
                    gcsdb = 0.0;
                    gcssb = (Capbs + cgso) * ag0;

                    gcbgb = -pParam.BSIM3cgbo * ag0;
                    gcbdb = -Capbd * ag0;
                    gcbsb = -Capbs * ag0;

                    CoxWL = ModelParameters.Cox * pParam.BSIM3weffCV
                                     * pParam.BSIM3leffCV;
                    if (Math.Abs(qcheq) <= 1.0e-5 * CoxWL)
                    {
                        if (ModelParameters.Xpart < 0.5)
                        {
                            dxpart = 0.4;
                        }
                        else if (ModelParameters.Xpart > 0.5)
                        {
                            dxpart = 0.0;
                        }
                        else
                        {
                            dxpart = 0.5;
                        }

                        ddxpart_dVd = ddxpart_dVg = ddxpart_dVb
                            = ddxpart_dVs = 0.0;
                    }
                    else
                    {
                        dxpart = qdrn / qcheq;
                        Cdd = Cddb;
                        Csd = -(Cgdb + Cddb + Cbdb);
                        ddxpart_dVd = (Cdd - dxpart * (Cdd + Csd)) / qcheq;
                        Cdg = Cdgb;
                        Csg = -(Cggb + Cdgb + Cbgb);
                        ddxpart_dVg = (Cdg - dxpart * (Cdg + Csg)) / qcheq;

                        Cds = Cdsb;
                        Css = -(Cgsb + Cdsb + Cbsb);
                        ddxpart_dVs = (Cds - dxpart * (Cds + Css)) / qcheq;

                        ddxpart_dVb = -(ddxpart_dVd + ddxpart_dVg + ddxpart_dVs);
                    }

                    sxpart = 1.0 - dxpart;
                    dsxpart_dVd = -ddxpart_dVd;
                    dsxpart_dVg = -ddxpart_dVg;
                    dsxpart_dVs = -ddxpart_dVs;
                    dsxpart_dVb = -(dsxpart_dVd + dsxpart_dVg + dsxpart_dVs);

                    qgd = qgdo;
                    qgs = qgso;
                    qgb = pParam.BSIM3cgbo * vgb;
                    qgate = qgd + qgs + qgb;
                    qbulk = -qgb;
                    qdrn = -qgd;
                    qsrc = -(qgate + qbulk + qdrn);
                }
            }
            else
            {
                if (BaseParameters.NqsMod == 0)
                {
                    gcggb = (Cggb + cgdo + cgso + pParam.BSIM3cgbo) * ag0;
                    gcgdb = (Cgsb - cgdo) * ag0;
                    gcgsb = (Cgdb - cgso) * ag0;

                    gcdgb = -(Cggb + Cbgb + Cdgb + cgdo) * ag0;
                    gcddb = (Capbd + cgdo - (Cgsb + Cbsb + Cdsb)) * ag0;
                    gcdsb = -(Cgdb + Cbdb + Cddb) * ag0;

                    gcsgb = (Cdgb - cgso) * ag0;
                    gcsdb = Cdsb * ag0;
                    gcssb = (Cddb + Capbs + cgso) * ag0;

                    gcbgb = (Cbgb - pParam.BSIM3cgbo) * ag0;
                    gcbdb = (Cbsb - Capbd) * ag0;
                    gcbsb = (Cbdb - Capbs) * ag0;

                    qgd = qgdo;
                    qgs = qgso;
                    qgb = pParam.BSIM3cgbo * vgb;
                    qgate += qgd + qgs + qgb;
                    qbulk -= qgb;
                    qsrc = qdrn - qgs;
                    qdrn = -(qgate + qbulk + qsrc);

                    ggtg = ggtd = ggtb = ggts = 0.0;
                    sxpart = 0.4;
                    dxpart = 0.6;
                    ddxpart_dVd = ddxpart_dVg = ddxpart_dVb = ddxpart_dVs = 0.0;
                    dsxpart_dVd = dsxpart_dVg = dsxpart_dVb = dsxpart_dVs = 0.0;
                }
                else
                {
                    if (qcheq > 0.0)
                        T0 = pParam.BSIM3tconst * qdef * ScalingFactor;
                    else
                        T0 = -pParam.BSIM3tconst * qdef * ScalingFactor;
                    ggtg = Gtg = T0 * Cqgb;
                    ggts = Gtd = T0 * Cqdb;
                    ggtd = Gts = T0 * Cqsb;
                    ggtb = Gtb = T0 * Cqbb;
                    gqdef = ScalingFactor * ag0;

                    gcqgb = Cqgb * ag0;
                    gcqdb = Cqsb * ag0;
                    gcqsb = Cqdb * ag0;
                    gcqbb = Cqbb * ag0;

                    gcggb = (cgdo + cgso + pParam.BSIM3cgbo) * ag0;
                    gcgdb = -cgdo * ag0;
                    gcgsb = -cgso * ag0;

                    gcdgb = -cgdo * ag0;
                    gcddb = (Capbd + cgdo) * ag0;
                    gcdsb = 0.0;

                    gcsgb = -cgso * ag0;
                    gcsdb = 0.0;
                    gcssb = (Capbs + cgso) * ag0;

                    gcbgb = -pParam.BSIM3cgbo * ag0;
                    gcbdb = -Capbd * ag0;
                    gcbsb = -Capbs * ag0;

                    CoxWL = ModelParameters.Cox * pParam.BSIM3weffCV
                                     * pParam.BSIM3leffCV;
                    if (Math.Abs(qcheq) <= 1.0e-5 * CoxWL)
                    {
                        if (ModelParameters.Xpart < 0.5)
                        {
                            sxpart = 0.4;
                        }
                        else if (ModelParameters.Xpart > 0.5)
                        {
                            sxpart = 0.0;
                        }
                        else
                        {
                            sxpart = 0.5;
                        }

                        dsxpart_dVd = dsxpart_dVg = dsxpart_dVb
                            = dsxpart_dVs = 0.0;
                    }
                    else
                    {
                        sxpart = qdrn / qcheq;
                        Css = Cddb;
                        Cds = -(Cgdb + Cddb + Cbdb);
                        dsxpart_dVs = (Css - sxpart * (Css + Cds)) / qcheq;
                        Csg = Cdgb;
                        Cdg = -(Cggb + Cdgb + Cbgb);
                        dsxpart_dVg = (Csg - sxpart * (Csg + Cdg)) / qcheq;

                        Csd = Cdsb;
                        Cdd = -(Cgsb + Cdsb + Cbsb);
                        dsxpart_dVd = (Csd - sxpart * (Csd + Cdd)) / qcheq;

                        dsxpart_dVb = -(dsxpart_dVd + dsxpart_dVg + dsxpart_dVs);
                    }

                    dxpart = 1.0 - sxpart;
                    ddxpart_dVd = -dsxpart_dVd;
                    ddxpart_dVg = -dsxpart_dVg;
                    ddxpart_dVs = -dsxpart_dVs;
                    ddxpart_dVb = -(ddxpart_dVd + ddxpart_dVg + ddxpart_dVs);

                    qgd = qgdo;
                    qgs = qgso;
                    qgb = pParam.BSIM3cgbo * vgb;
                    qgate = qgd + qgs + qgb;
                    qbulk = -qgb;
                    qsrc = -qgs;
                    qdrn = -(qgate + qbulk + qsrc);
                }
            }

            cqdef = cqcheq = 0.0;

            if (TranBehavior != null)
            {
                TranBehavior.Qg.Current = qgate;
                TranBehavior.Qd.Current = qdrn - Qbd;
                TranBehavior.Qb.Current = qbulk + Qbd + Qbs;

                if (BaseParameters.NqsMod > 0)
                {
                    TranBehavior.Qcdump.Current = qdef * ScalingFactor;
                    TranBehavior.Qcheq.Current = qcheq;
                }
            }

            /* store small signal parameters */
            if (Simulation is FrequencySimulation && !state.UseDc)
                goto line1000;
            if (!chargeComputationNeeded)
                goto line850;

            if (TranBehavior != null)
            {
                TranBehavior.Qb.Integrate();
                TranBehavior.Qg.Integrate();
                TranBehavior.Qd.Integrate();
                if (BaseParameters.NqsMod > 0)
                {
                    TranBehavior.Qcdump.Integrate();
                    TranBehavior.Qcheq.Integrate();
                }
            }

            goto line860;

            line850:
            /* initialize to zero charge conductance and current */
            ceqqg = ceqqb = ceqqd = 0.0;
            cqcheq = cqdef = 0.0;

            gcdgb = gcddb = gcdsb = 0.0;
            gcsgb = gcsdb = gcssb = 0.0;
            gcggb = gcgdb = gcgsb = 0.0;
            gcbgb = gcbdb = gcbsb = 0.0;

            gqdef = gcqgb = gcqdb = gcqsb = gcqbb = 0.0;
            ggtg = ggtd = ggtb = ggts = 0.0;
            sxpart = 1.0 - (dxpart = Mode > 0 ? 0.4 : 0.6);
            ddxpart_dVd = ddxpart_dVg = ddxpart_dVb = ddxpart_dVs = 0.0;
            dsxpart_dVd = dsxpart_dVg = dsxpart_dVb = dsxpart_dVs = 0.0;

            if (BaseParameters.NqsMod > 0)
                Gtau = 16.0 * pParam.BSIM3u0temp * ModelTemperature.Vtm / pParam.BSIM3leffCV / pParam.BSIM3leffCV *
                       ScalingFactor;
            else
                Gtau = 0.0;

            goto line900;

            line860:
            /* evaluate equivalent charge current */

            cqgate = TranBehavior.Qg.Derivative;
            cqbulk = TranBehavior.Qb.Derivative;
            cqdrn = TranBehavior.Qd.Derivative;

            ceqqg = cqgate - gcggb * vgb + gcgdb * vbd + gcgsb * vbs;
            ceqqb = cqbulk - gcbgb * vgb + gcbdb * vbd + gcbsb * vbs;
            ceqqd = cqdrn - gcdgb * vgb + gcddb * vbd + gcdsb * vbs;

            if (BaseParameters.NqsMod > 0)
            {
                T0 = ggtg * vgb - ggtd * vbd - ggts * vbs;
                ceqqg += T0;
                T1 = qdef * Gtau;
                ceqqd -= dxpart * T0 + T1 * (ddxpart_dVg * vgb - ddxpart_dVd
                                             * vbd - ddxpart_dVs * vbs);
                cqdef = TranBehavior.Qcdump.Derivative - gqdef * qdef;
                cqcheq = TranBehavior.Qcheq.Derivative - (gcqgb * vgb - gcqdb * vbd - gcqsb * vbs) + T0;
            }

            /*
             *  load current vector
             */
            line900:

            if (Mode >= 0)
            {
                Gm = this.Gm;
                Gmbs = this.Gmbs;
                FwdSum = Gm + Gmbs;
                RevSum = 0.0;
                cdreq = ModelParameters.B3Type * (cdrain - this.Gds * vds - Gm * vgs - Gmbs * vbs);

                ceqbd = -ModelParameters.B3Type * (Csub - Gbds * vds - Gbgs * vgs - Gbbs * vbs);
                ceqbs = 0.0;

                gbbdp = -Gbds;
                gbbsp = Gbds + Gbgs + Gbbs;

                gbdpg = Gbgs;
                gbdpdp = Gbds;
                gbdpb = Gbbs;
                gbdpsp = -(gbdpg + gbdpdp + gbdpb);

                gbspg = 0.0;
                gbspdp = 0.0;
                gbspb = 0.0;
                gbspsp = 0.0;
            }
            else
            {
                Gm = -this.Gm;
                Gmbs = -this.Gmbs;
                FwdSum = 0.0;
                RevSum = -(Gm + Gmbs);
                cdreq = -ModelParameters.B3Type * (cdrain + this.Gds * vds + Gm * vgd + Gmbs * vbd);

                ceqbs = -ModelParameters.B3Type * (Csub + Gbds * vds - Gbgs * vgd - Gbbs * vbd);
                ceqbd = 0.0;

                gbbsp = -Gbds;
                gbbdp = Gbds + Gbgs + Gbbs;

                gbdpg = 0.0;
                gbdpsp = 0.0;
                gbdpb = 0.0;
                gbdpdp = 0.0;

                gbspg = Gbgs;
                gbspsp = Gbds;
                gbspb = Gbbs;
                gbspdp = -(gbspg + gbspsp + gbspb);
            }

            if (ModelParameters.B3Type > 0)
            {
                ceqbs += Cbs - Gbs * vbs;
                ceqbd += this.Cbd - this.Gbd * vbd;
                /*
                    ceqqg = ceqqg;
                    ceqqb = ceqqb;
                    ceqqd = ceqqd;
                    cqdef = cqdef;
                    cqcheq = cqcheq;
                */
            }
            else
            {
                ceqbs -= Cbs - Gbs * vbs;
                ceqbd -= this.Cbd - this.Gbd * vbd;
                ceqqg = -ceqqg;
                ceqqb = -ceqqb;
                ceqqd = -ceqqd;
                cqdef = -cqdef;
                cqcheq = -cqcheq;
            }

            GateNodePtr.Value -= ceqqg;
            BulkNodePtr.Value -= ceqbs + ceqbd + ceqqb;
            DrainNodePrimePtr.Value += ceqbd - cdreq - ceqqd;
            SourceNodePrimePtr.Value += cdreq + ceqbs + ceqqg + ceqqb + ceqqd;
            if (BaseParameters.NqsMod > 0)
                QNodePtr.Value += cqcheq - cqdef;

            /*
             *  load y matrix
             */

            T1 = qdef * Gtau;
            DdPtr.Value += base.DrainConductance;
            GgPtr.Value += gcggb - ggtg;
            SsPtr.Value += base.SourceConductance;
            BbPtr.Value += this.Gbd + Gbs - gcbgb - gcbdb - gcbsb - Gbbs;
            DPdpPtr.Value += base.DrainConductance + this.Gds + this.Gbd + RevSum + gcddb + dxpart * ggtd +
                             T1 * ddxpart_dVd + gbdpdp;
            SPspPtr.Value += base.SourceConductance + this.Gds + Gbs + FwdSum + gcssb + sxpart * ggts +
                             T1 * dsxpart_dVs + gbspsp;
            DdpPtr.Value -= base.DrainConductance;
            GbPtr.Value -= gcggb + gcgdb + gcgsb + ggtb;
            GdpPtr.Value += gcgdb - ggtd;
            GspPtr.Value += gcgsb - ggts;
            SspPtr.Value -= base.SourceConductance;
            BgPtr.Value += gcbgb - Gbgs;
            BdpPtr.Value += gcbdb - this.Gbd + gbbdp;
            BspPtr.Value += gcbsb - Gbs + gbbsp;
            DPdPtr.Value -= base.DrainConductance;
            DPgPtr.Value += Gm + gcdgb + dxpart * ggtg + T1 * ddxpart_dVg + gbdpg;
            DPbPtr.Value -= this.Gbd - Gmbs + gcdgb + gcddb + gcdsb - dxpart * ggtb - T1 * ddxpart_dVb - gbdpb;
            DPspPtr.Value -= this.Gds + FwdSum - gcdsb - dxpart * ggts - T1 * ddxpart_dVs - gbdpsp;
            SPgPtr.Value += gcsgb - Gm + sxpart * ggtg + T1 * dsxpart_dVg + gbspg;
            SPsPtr.Value -= base.SourceConductance;
            SPbPtr.Value -= Gbs + Gmbs + gcsgb + gcsdb + gcssb - sxpart * ggtb - T1 * dsxpart_dVb - gbspb;
            SPdpPtr.Value -= this.Gds + RevSum - gcsdb - sxpart * ggtd - T1 * dsxpart_dVd - gbspdp;

            if (BaseParameters.NqsMod > 0)
            {
                QqPtr.Value += gqdef + Gtau;

                DPqPtr.Value += dxpart * Gtau;
                SPqPtr.Value += sxpart * Gtau;
                GqPtr.Value -= Gtau;

                QgPtr.Value += ggtg - gcqgb;
                QdpPtr.Value += ggtd - gcqdb;
                QspPtr.Value += ggts - gcqsb;
                QbPtr.Value += ggtb - gcqbb;
            }

            line1000: ;
        }

        /// <summary>
        /// Determines whether the specified simulation is convergent.
        /// </summary>
        /// <returns>
        ///   <c>true</c> if the specified simulation is convergent; otherwise, <c>false</c>.
        /// </returns>
        bool IBiasingBehavior.IsConvergent() => true;
    }
}