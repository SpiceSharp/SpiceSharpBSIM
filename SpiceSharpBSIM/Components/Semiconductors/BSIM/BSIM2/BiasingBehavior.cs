using System;
using SpiceSharp.Simulations;
using SpiceSharp.Behaviors;
using SpiceSharp.Algebra;
using SpiceSharp.Components.MosfetBehaviors;
using SpiceSharp.Components.Semiconductors;

namespace SpiceSharp.Components.BSIM2Behaviors
{	
	/// <summary>
	/// Load behavior for a <see cref="BSIM2" />
	/// </summary>
	public class BiasingBehavior : TemperatureBehavior, IBiasingBehavior, IConnectedBehavior
	{
        /// <summary>
        /// Gets the base configuration.
        /// </summary>
        /// <value>
        /// The base configuration.
        /// </value>
        protected BaseConfiguration BaseConfiguration { get; private set; }
		
		/// <summary>
		/// Properties
		/// </summary>
		public double Mode { get; private set; }
		public double Vdsat { get; private set; }
		public TransientBehavior TranBehavior { get; set; }
		public double Vbd { get; set; }
		public double Vbs { get; set; }
		public double Vgs { get; set; }
		public double Vds { get; set; }
		public double Gm { get; set; }
		public double Gds { get; set; }
		public double Gmbs { get; set; }
		public double Gbd { get; set; }
		public double Gbs { get; set; }
		public double Vono { get; set; }
		public double Vdsato { get; set; }
		public double Cd { get; set; }
		public double Cbs { get; set; }
		public double Cbd { get; set; }
	    public double Cggb { get; set; }
	    public double Cgdb { get; set; }
	    public double Cgsb { get; set; }
	    public double Cbgb { get; set; }
	    public double Cbdb { get; set; }
	    public double Cbsb { get; set; }
	    public double Capbd { get; set; }
	    public double Iqbd { get; set; }
	    public double Capbs { get; set; }
	    public double Iqbs { get; set; }
	    public double Cdgb { get; set; }
	    public double Cddb { get; set; }
	    public double Cdsb { get; set; }
	    public double Qbs { get; set; }
	    public double Qbd { get; set; }
		
		/// <summary>
		/// Nodes
		/// </summary>
		protected int DrainNode { get; private set; }
	    protected int GateNode { get; private set; }
	    protected int SourceNode { get;private set; }
	    protected int BulkNode { get; private set; }
		protected VectorElement<double> DrainNodePrimePtr { get; private set; }
		public int DrainNodePrime { get; private set; }
		protected VectorElement<double> SourceNodePrimePtr { get; private set; }
		public int SourceNodePrime { get; private set; }
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
		protected VectorElement<double> GateNodePtr { get; private set; }
		protected VectorElement<double> BulkNodePtr { get; private set; }

	    /// <summary>
	    /// Constructor
	    /// </summary>
	    public BiasingBehavior(string name) : base(name)
	    {
	    }
		
		/// <summary>
		/// Setup the behavior
		/// </summary>
		public override void Setup(Simulation simulation, SetupDataProvider provider)
		{
			if (provider == null)
				throw new ArgumentNullException(nameof(provider));
            base.Setup(simulation, provider);

            // Get configuration
		    BaseConfiguration = simulation.Configurations.Get<BaseConfiguration>();
		}
		
		/// <summary>
		/// Connect
		/// </summary>
		public void Connect(params int[] pins)
		{
			DrainNode = pins[0];
			GateNode = pins[1];
			SourceNode = pins[2];
			BulkNode = pins[3];
		}
		
		/// <summary>
		/// Get equation pointers
		/// </summary>
		public void GetEquationPointers(VariableSet variables, Solver<double> solver)
		{
			if ((ModelParameters.SheetResistance > 0) && (BaseParameters.DrainSquares > 0.0))
				DrainNodePrime = variables.Create(Name.Combine("drain")).Index;
			else
				DrainNodePrime = DrainNode;
			DrainNodePrimePtr = solver.GetRhsElement(DrainNodePrime);
			if ((ModelParameters.SheetResistance > 0) && (BaseParameters.SourceSquares > 0.0))
					SourceNodePrime = variables.Create(Name.Combine("source")).Index;
			else
				SourceNodePrime = SourceNode;

			SourceNodePrimePtr = solver.GetRhsElement(SourceNodePrime);
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
			GateNodePtr = solver.GetRhsElement(GateNode);
			BulkNodePtr = solver.GetRhsElement(BulkNode);
		}

	    /// <summary>
	    /// Temperature behavior
	    /// </summary>
	    public void Load(BaseSimulation simulation)
	    {
	        double DrainSatCurrent;
	        double EffectiveLength;
	        double GateBulkOverlapCap;
	        double GateDrainOverlapCap;
	        double GateSourceOverlapCap;
	        double SourceSatCurrent;
	        double DrainArea;
	        double SourceArea;
	        double DrainPerimeter;
	        double SourcePerimeter;
	        double arg;
	        double capbd = 0.0;
	        double capbs = 0.0;
	        double cbd;
	        double cbs;
	        double cd;
	        double cdrain;
	        double cdreq;
	        double ceqbd;
	        double ceqbs;
	        double ceqqb;
	        double ceqqd;
	        double ceqqg;
	        double czbd;
	        double czbdsw;
	        double czbs;
	        double czbssw;
	        double evbd;
	        double evbs;
	        double gbd;
	        double gbs;
	        double gcbdb;
	        double gcbgb;
	        double gcbsb;
	        double gcddb;
	        double gcdgb;
	        double gcdsb;
	        double gcgdb;
	        double gcggb;
	        double gcgsb;
	        double gcsdb;
	        double gcsgb;
	        double gcssb;
	        double gds;
	        double gm;
	        double gmbs;
	        double sarg;
	        double sargsw;
	        double vbd;
	        double vbs;
	        double vcrit;
	        double vds;
	        double vdsat;
	        double vgb;
	        double vgd;
	        double vgdo;
	        double vgs;
	        double von;
	        double xnrm;
	        double xrev;
	        bool Check;
	        double cgdb;
	        double cgsb;
	        double cbdb;
	        double cdgb = 0.0;
	        double cddb = 0.0;
	        double cdsb = 0.0;
	        double cggb;
	        double cbgb;
	        double cbsb;
	        double csgb = 0.0;
	        double cssb = 0.0;
	        double csdb = 0.0;
	        double PhiB;
	        double PhiBSW;
	        double MJ;
	        double MJSW;
	        double argsw;
	        double qgate;
	        double qbulk;
	        double qdrn = 0.0;
	        double qsrc = 0.0;
	        double cqgate;
	        double cqbulk;
	        double cqdrn;
	        double vt0;
	        double[] args = new double[8];
	        var pParam = base.Param;
	        var state = simulation.RealState;
	        bool chargeComputationNeeded = TranBehavior != null;

	        EffectiveLength = BaseParameters.Length - ModelParameters.DeltaL * 1.0e-6; /* m */
	        DrainArea = BaseParameters.DrainArea;
	        SourceArea = BaseParameters.SourceArea;
	        DrainPerimeter = BaseParameters.DrainPerimeter;
	        SourcePerimeter = BaseParameters.SourcePerimeter;
	        if ((DrainSatCurrent = DrainArea * ModelParameters.JctSatCurDensity)
	            < 1e-15)
	        {
	            DrainSatCurrent = 1.0e-15;
	        }

	        if ((SourceSatCurrent = SourceArea * ModelParameters.JctSatCurDensity)
	            < 1.0e-15)
	        {
	            SourceSatCurrent = 1.0e-15;
	        }

	        GateSourceOverlapCap = ModelParameters.GateSourceOverlapCap * BaseParameters.Width;
	        GateDrainOverlapCap = ModelParameters.GateDrainOverlapCap * BaseParameters.Width;
	        GateBulkOverlapCap = ModelParameters.GateBulkOverlapCap * EffectiveLength;
	        von = ModelParameters.Type * base.Von;
	        vdsat = ModelParameters.Type * this.Vdsat;
	        vt0 = ModelParameters.Type * pParam.B2vt0;

	        Check = true;
	        if (simulation is FrequencySimulation && !state.UseDc)
	        {
	            vbs = this.Vbs;
	            vgs = this.Vgs;
	            vds = this.Vds;
	            chargeComputationNeeded = true;
	        }
	        else if (state.Init == InitializationModes.Junction && !BaseParameters.Off)
	        {
	            vds = ModelParameters.Type * BaseParameters.IcVDS;
	            vgs = ModelParameters.Type * BaseParameters.IcVGS;
	            vbs = ModelParameters.Type * BaseParameters.IcVBS;
	            if ((vds == 0) && (vgs == 0) && (vbs == 0) &&
	                (state.UseDc || TranBehavior != null || !state.UseIc))
	            {
	                vbs = -1;
	                vgs = vt0;
	                vds = 0;
	            }
	        }
	        else if ((state.Init == InitializationModes.Junction ||
	                  state.Init == InitializationModes.Fix) && (BaseParameters.Off))
	        {
	            vbs = vgs = vds = 0;
	        }
	        else
	        {
	            vbs = ModelParameters.Type * (state.Solution[BulkNode] - state.Solution[SourceNodePrime]);
	            vgs = ModelParameters.Type * (state.Solution[GateNode] - state.Solution[SourceNodePrime]);
	            vds = ModelParameters.Type * (state.Solution[DrainNodePrime] - state.Solution[SourceNodePrime]);
	            vbd = vbs - vds;
	            vgd = vgs - vds;
	            vgdo = this.Vgs - this.Vds;

	            von = ModelParameters.Type * base.Von;
	            if (this.Vds >= 0)
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

	            Check = false;
	            if (vds >= 0)
	            {
	                vcrit = Circuit.Vt0 * Math.Log(Circuit.Vt0 / (Circuit.Root2 * SourceSatCurrent));
	                vbs = Semiconductor.LimitJunction(vbs, this.Vbs, Circuit.Vt0, vcrit, ref Check); /* B2 test */
	                vbd = vbs - vds;
	            }
	            else
	            {
	                vcrit = Circuit.Vt0 * Math.Log(Circuit.Vt0 / (Circuit.Root2 * DrainSatCurrent));
	                vbd = Semiconductor.LimitJunction(vbd, this.Vbd, Circuit.Vt0, vcrit, ref Check); /* B2 test*/
	                vbs = vbd + vds;
	            }
	        }

	        /* determine DC current and derivatives */
	        vbd = vbs - vds;
	        vgd = vgs - vds;
	        vgb = vgs - vbs;

	        if (vbs <= 0.0)
	        {
	            gbs = SourceSatCurrent / Circuit.Vt0 + BaseConfiguration.Gmin;
	            cbs = gbs * vbs;
	        }
	        else
	        {
	            evbs = Math.Exp(vbs / Circuit.Vt0);
	            gbs = SourceSatCurrent * evbs / Circuit.Vt0 + BaseConfiguration.Gmin;
	            cbs = SourceSatCurrent * (evbs - 1) + BaseConfiguration.Gmin * vbs;
	        }

	        if (vbd <= 0.0)
	        {
	            gbd = DrainSatCurrent / Circuit.Vt0 + BaseConfiguration.Gmin;
	            cbd = gbd * vbd;
	        }
	        else
	        {
	            evbd = Math.Exp(vbd / Circuit.Vt0);
	            gbd = DrainSatCurrent * evbd / Circuit.Vt0 + BaseConfiguration.Gmin;
	            cbd = DrainSatCurrent * (evbd - 1) + BaseConfiguration.Gmin * vbd;
	        }

	        /* line 400 */
	        if (vds >= 0)
	        {
	            /* normal mode */
	            this.Mode = 1;
	        }
	        else
	        {
	            /* inverse mode */
	            this.Mode = -1;
	        }

	        /* call B2evaluate to calculate drain current and its 
             * derivatives and charge and capacitances related to gate
             * drain, and bulk
             */
	        if (vds >= 0)
	        {
	            Evaluate(chargeComputationNeeded, vds, vbs, vgs, out gm, out gds, out gmbs, out qgate,
	                out qbulk, out qdrn, out cggb, out cgdb, out cgsb, out cbgb, out cbdb, out cbsb, out cdgb,
	                out cddb, out cdsb, out cdrain, out von, out vdsat);
	        }
	        else
	        {
	            Evaluate(chargeComputationNeeded, -vds, vbd, vgd, out gm, out gds, out gmbs, out qgate,
	                out qbulk, out qsrc, out cggb, out cgsb, out cgdb, out cbgb, out cbsb, out cbdb, out csgb,
	                out cssb, out csdb, out cdrain, out von, out vdsat);
	        }

	        base.Von = ModelParameters.Type * von;
	        this.Vdsat = ModelParameters.Type * vdsat;

	        /*
             *  COMPUTE EQUIVALENT DRAIN CURRENT SOURCE
             */
	        cd = this.Mode * cdrain - cbd;
	        if (chargeComputationNeeded)
	        {
	            /*
                 *  charge storage elements
                 *
                 *   bulk-drain and bulk-source depletion capacitances
                 *  czbd : zero bias drain junction capacitance
                 *  czbs : zero bias source junction capacitance
                 * czbdsw:zero bias drain junction sidewall capacitance
                 * czbssw:zero bias source junction sidewall capacitance
                 */

	            czbd = ModelParameters.UnitAreaJctCap * DrainArea;
	            czbs = ModelParameters.UnitAreaJctCap * SourceArea;
	            czbdsw = ModelParameters.UnitLengthSidewallJctCap * DrainPerimeter;
	            czbssw = ModelParameters.UnitLengthSidewallJctCap * SourcePerimeter;
	            PhiB = ModelParameters.BulkJctPotential;
	            PhiBSW = ModelParameters.SidewallJctPotential;
	            MJ = ModelParameters.BulkJctBotGradingCoeff;
	            MJSW = ModelParameters.BulkJctSideGradingCoeff;

	            /* Source Bulk Junction */
	            if (vbs < 0)
	            {
	                arg = 1 - vbs / PhiB;
	                argsw = 1 - vbs / PhiBSW;
	                sarg = Math.Exp(-MJ * Math.Log(arg));
	                sargsw = Math.Exp(-MJSW * Math.Log(argsw));

	                this.Qbs = PhiB * czbs * (1 - arg * sarg) / (1 - MJ) +
	                           PhiBSW * czbssw * (1 - argsw * sargsw) / (1 - MJSW);
	                capbs = czbs * sarg + czbssw * sargsw;
	            }
	            else
	            {
	                this.Qbs = vbs * (czbs + czbssw) +
	                           vbs * vbs * (czbs * MJ * 0.5 / PhiB + czbssw * MJSW * 0.5 / PhiBSW);
	                capbs = czbs + czbssw + vbs * (czbs * MJ / PhiB + czbssw * MJSW / PhiBSW);
	            }

	            /* Drain Bulk Junction */
	            if (vbd < 0)
	            {
	                arg = 1 - vbd / PhiB;
	                argsw = 1 - vbd / PhiBSW;
	                sarg = Math.Exp(-MJ * Math.Log(arg));
	                sargsw = Math.Exp(-MJSW * Math.Log(argsw));
	                this.Qbd = PhiB * czbd * (1 - arg * sarg) / (1 - MJ) +
	                           PhiBSW * czbdsw * (1 - argsw * sargsw) / (1 - MJSW);
	                capbd = czbd * sarg + czbdsw * sargsw;
	            }
	            else
	            {
	                this.Qbd = vbd * (czbd + czbdsw) +
	                           vbd * vbd * (czbd * MJ * 0.5 / PhiB + czbdsw * MJSW * 0.5 / PhiBSW);
	                capbd = czbd + czbdsw + vbd * (czbd * MJ / PhiB + czbdsw * MJSW / PhiBSW);
	            }
	        }

	        /*
             *  check convergence
             */
	        if (!BaseParameters.Off || state.Init != InitializationModes.Fix)
	        {
	            if (Check)
	            {
	                state.IsConvergent = false;
	            }
	        }

	        this.Vbs = vbs;
	        this.Vbd = vbd;
	        this.Vgs = vgs;
	        this.Vds = vds;
	        this.Cd = cd;
	        this.Cbs = cbs;
	        this.Cbd = cbd;
	        this.Gm = gm;
	        this.Gds = gds;
	        this.Gmbs = gmbs;
	        this.Gbd = gbd;
	        this.Gbs = gbs;

	        this.Cggb = cggb;
	        this.Cgdb = cgdb;
	        this.Cgsb = cgsb;

	        this.Cbgb = cbgb;
	        this.Cbdb = cbdb;
	        this.Cbsb = cbsb;

	        this.Cdgb = cdgb;
	        this.Cddb = cddb;
	        this.Cdsb = cdsb;

	        this.Capbs = capbs;
	        this.Capbd = capbd;

	        /* bulk and channel charge plus overlaps */

	        // if((!(ckt->CKTmode & (MODETRAN | MODEAC))) && ((!(ckt->CKTmode & MODETRANOP)) ||  (!(ckt->CKTmode & MODEUIC)))  && (!(ckt->CKTmode  &  MODEINITSMSIG)))
	        if (!chargeComputationNeeded)
	            goto line850;

	        line755:
	        if (this.Mode > 0)
	        {

	            args[0] = GateDrainOverlapCap;
	            args[1] = GateSourceOverlapCap;
	            args[2] = GateBulkOverlapCap;
	            args[3] = capbd;
	            args[4] = capbs;
	            args[5] = cggb;
	            args[6] = cgdb;
	            args[7] = cgsb;

	            MosCap(vgd, vgs, vgb,
	                args,
	                /*
                    GateDrainOverlapCap,
                    GateSourceOverlapCap,GateBulkOverlapCap,
                    capbd,capbs,cggb,cgdb,cgsb,
                    */
	                cbgb, cbdb, cbsb, cdgb, cddb, cdsb, out gcggb, out gcgdb, out gcgsb, out gcbgb, out gcbdb, out gcbsb,
	                out gcdgb, out gcddb, out gcdsb, out gcsgb, out gcsdb, out gcssb, ref qgate, ref qbulk, ref qdrn,
	                out qsrc);
	        }
	        else
	        {

	            args[0] = GateSourceOverlapCap;
	            args[1] = GateDrainOverlapCap;
	            args[2] = GateBulkOverlapCap;
	            args[3] = capbs;
	            args[4] = capbd;
	            args[5] = cggb;
	            args[6] = cgsb;
	            args[7] = cgdb;

	            MosCap(vgs, vgd, vgb, args,
	                /*
                    GateSourceOverlapCap,
                            GateDrainOverlapCap,GateBulkOverlapCap,
                    capbs,capbd,cggb,cgsb,cgdb,
                    */
	                cbgb, cbsb, cbdb, csgb, cssb, csdb, out gcggb, out gcgsb, out gcgdb, out gcbgb, out gcbsb, out gcbdb,
	                out gcsgb, out gcssb, out gcsdb, out gcdgb, out gcdsb, out gcddb, ref qgate, ref qbulk, ref qsrc,
	                out qdrn);
	        }



	        /* store small signal parameters */
	        // if ((!(ckt->CKTmode & (MODEAC | MODETRAN))) && (ckt->CKTmode & MODETRANOP) && (ckt->CKTmode & MODEUIC))
	        // goto line850;

	        if (simulation is FrequencySimulation && !state.UseDc)
	        {
	            this.Cggb = cggb;
	            this.Cgdb = cgdb;
	            this.Cgsb = cgsb;
	            this.Cbgb = cbgb;
	            this.Cbdb = cbdb;
	            this.Cbsb = cbsb;
	            this.Cdgb = cdgb;
	            this.Cddb = cddb;
	            this.Cdsb = cdsb;
	            this.Capbd = capbd;
	            this.Capbs = capbs;

	            goto line1000;
	        }

	        if (TranBehavior != null)
	        {
	            TranBehavior.Qg.Current = qgate;
	            TranBehavior.Qd.Current = qdrn - this.Qbd;
	            TranBehavior.Qb.Current = qbulk + this.Qbd + this.Qbs;
	            TranBehavior.Qb.Integrate();
	            TranBehavior.Qg.Integrate();
	            TranBehavior.Qd.Integrate();
	        }

	        goto line860;

	        line850:
	        /* initialize to zero charge conductance and current */
	        ceqqg = ceqqb = ceqqd = 0.0;
	        gcdgb = gcddb = gcdsb = 0.0;
	        gcsgb = gcsdb = gcssb = 0.0;
	        gcggb = gcgdb = gcgsb = 0.0;
	        gcbgb = gcbdb = gcbsb = 0.0;
	        goto line900;

	        line860:
	        /* evaluate equivalent charge current */
	        cqgate = TranBehavior.Qg.Derivative;
	        cqbulk = TranBehavior.Qb.Derivative;
	        cqdrn = TranBehavior.Qd.Derivative;
	        ceqqg = cqgate - gcggb * vgb + gcgdb * vbd + gcgsb * vbs;
	        ceqqb = cqbulk - gcbgb * vgb + gcbdb * vbd + gcbsb * vbs;
	        ceqqd = cqdrn - gcdgb * vgb + gcddb * vbd + gcdsb * vbs;

	        /*
             *  load current vector
             */
	        line900:

	        ceqbs = ModelParameters.Type * (cbs - (gbs - BaseConfiguration.Gmin) * vbs);
	        ceqbd = ModelParameters.Type * (cbd - (gbd - BaseConfiguration.Gmin) * vbd);

	        ceqqg = ModelParameters.Type * ceqqg;
	        ceqqb = ModelParameters.Type * ceqqb;
	        ceqqd = ModelParameters.Type * ceqqd;
	        if (this.Mode >= 0)
	        {
	            xnrm = 1;
	            xrev = 0;
	            cdreq = ModelParameters.Type * (cdrain - gds * vds - gm * vgs - gmbs * vbs);
	        }
	        else
	        {
	            xnrm = 0;
	            xrev = 1;
	            cdreq = -ModelParameters.Type * (cdrain + gds * vds - gm * vgd - gmbs * vbd);
	        }

	        var m = BaseParameters.Multiplier;
	        GateNodePtr.Value -= m * ceqqg;
	        BulkNodePtr.Value -= m * (ceqbs + ceqbd + ceqqb);
	        DrainNodePrimePtr.Value += m * (ceqbd - cdreq - ceqqd);
	        SourceNodePrimePtr.Value += m * (cdreq + ceqbs + ceqqg + ceqqb + ceqqd);

	        /*
             *  load y matrix
             */
	        DdPtr.Value += m * (base.DrainConductance);
	        GgPtr.Value += m * (gcggb);
	        SsPtr.Value += m * (base.SourceConductance);
	        BbPtr.Value += m * (gbd + gbs - gcbgb - gcbdb - gcbsb);
	        DPdpPtr.Value += m * (base.DrainConductance + gds + gbd + xrev * (gm + gmbs) + gcddb);
	        SPspPtr.Value += m * (base.SourceConductance + gds + gbs + xnrm * (gm + gmbs) + gcssb);
	        DdpPtr.Value += m * (-base.DrainConductance);
	        GbPtr.Value += m * (-gcggb - gcgdb - gcgsb);
	        GdpPtr.Value += m * (gcgdb);
	        GspPtr.Value += m * (gcgsb);
	        SspPtr.Value += m * (-base.SourceConductance);
	        BgPtr.Value += m * (gcbgb);
	        BdpPtr.Value += m * (-gbd + gcbdb);
	        BspPtr.Value += m * (-gbs + gcbsb);
	        DPdPtr.Value += m * (-base.DrainConductance);
	        DPgPtr.Value += m * ((xnrm - xrev) * gm + gcdgb);
	        DPbPtr.Value += m * (-gbd + (xnrm - xrev) * gmbs - gcdgb - gcddb - gcdsb);
	        DPspPtr.Value += m * (-gds - xnrm * (gm + gmbs) + gcdsb);
	        SPgPtr.Value += m * (-(xnrm - xrev) * gm + gcsgb);
	        SPsPtr.Value += m * (-base.SourceConductance);
	        SPbPtr.Value += m * (-gbs - (xnrm - xrev) * gmbs - gcsgb - gcsdb - gcssb);
	        SPdpPtr.Value += m * (-gds - xrev * (gm + gmbs) + gcsdb);

	        line1000: ;
	    }

        /// <summary>
        /// Determines whether the specified simulation is convergent.
        /// </summary>
        /// <param name="simulation">The simulation.</param>
        /// <returns>
        ///   <c>true</c> if the specified simulation is convergent; otherwise, <c>false</c>.
        /// </returns>
        public bool IsConvergent(BaseSimulation simulation) => true;

	    /// <summary>
        /// Helper method Evaluate
        /// </summary>
        public void Evaluate(bool chargeComputationNeeded, double vds, double vbs, double vgs, out double gm, out double gds, out double gmb, out double qg, out double qb, out double qd, out double cgg, out double cgd, out double cgs, out double cbg, out double cbd, out double cbs, out double cdg, out double cdd, out double cds, out double Ids, out double von, out double Vdsat)
        {
            double vth, vdsat = 0, phisb, t1s, eta, gg, aa, inv_Aa, u1, u1s, vc, kk, sqrtKk, dPhisb_dVb, dT1s_dVb, dVth_dVb, dVth_dVd, dAa_dVb, dVc_dVd, dVc_dVg, dVc_dVb, dKk_dVc, dVdsat_dVd = 0, dVdsat_dVg = 0, dVdsat_dVb = 0, dUvert_dVg, dUvert_dVd, dUvert_dVb, inv_Kk, dUtot_dVd, dUtot_dVb, dUtot_dVg, ai, bi, vghigh, vglow, vgeff, vof, vbseff, vgst, vgdt, qbulk, utot, t0, t1, t2, t3, t4, t5, arg1, arg2, exp0 = 0, tmp, tmp1, tmp2, tmp3, uvert, beta1, beta2, beta0, dGg_dVb, exp1 = 0, t6, t7, t8, t9, n = 0, expArg, expArg1, beta, dQbulk_dVb, dVgdt_dVg, dVgdt_dVd, dVbseff_dVb, ua, ub, dVgdt_dVb, dQbulk_dVd, con1, con3, con4, sqrVghigh, sqrVglow, cubVghigh, cubVglow, delta, coeffa, coeffb, coeffc, coeffd, inv_Uvert, inv_Utot, inv_Vdsat, tanh, sqrsech, dBeta1_dVb, dU1_dVd, dU1_dVg, dU1_dVb, betaeff, fR, dFR_dVd, dFR_dVg, dFR_dVb, betas, beta3, beta4, dBeta_dVd, dBeta_dVg, dBeta_dVb, dVgeff_dVg, dVgeff_dVd, dVgeff_dVb, dCon3_dVd, dCon3_dVb, dCon4_dVd, dCon4_dVb, dCoeffa_dVd, dCoeffa_dVb, dCoeffb_dVd, dCoeffb_dVb, dCoeffc_dVd, dCoeffc_dVb, dCoeffd_dVd, dCoeffd_dVb;

            if (vbs < ModelParameters.Vbb2)
                vbs = ModelParameters.Vbb2;
            if (vgs > ModelParameters.Vgg2)
                vgs = ModelParameters.Vgg2;
            if (vds > ModelParameters.Vdd2)
                vds = ModelParameters.Vdd2;
            if (vbs <= 0.0)
            {
                phisb = base.Param.B2phi - vbs;
                dPhisb_dVb = -1.0;
                t1s = Math.Sqrt(phisb);
                dT1s_dVb = -0.5 / t1s;
            }
            else
            {
                tmp = base.Param.B2phi / (base.Param.B2phi + vbs);
                phisb = base.Param.B2phi * tmp;
                dPhisb_dVb = -tmp * tmp;
                t1s = base.Param.Phis3 / (base.Param.B2phi + 0.5 * vbs);
                dT1s_dVb = -0.5 * t1s * t1s / base.Param.Phis3;
            }
            eta = base.Param.B2eta0 + base.Param.B2etaB * vbs;
            ua = base.Param.B2ua0 + base.Param.B2uaB * vbs;
            ub = base.Param.B2ub0 + base.Param.B2ubB * vbs;
            u1s = base.Param.B2u10 + base.Param.B2u1B * vbs;
            vth = base.Param.B2vfb + base.Param.B2phi + base.Param.B2k1 * t1s - base.Param.B2k2 * phisb - eta * vds;
            dVth_dVd = -eta;
            dVth_dVb = base.Param.B2k1 * dT1s_dVb + base.Param.B2k2 - base.Param.B2etaB * vds;
            vgst = vgs - vth;
            tmp = 1.0 / (1.744 + 0.8364 * phisb);
            gg = 1.0 - tmp;
            dGg_dVb = 0.8364 * tmp * tmp * dPhisb_dVb;
            t0 = gg / t1s;
            tmp1 = 0.5 * t0 * base.Param.B2k1;
            aa = 1.0 + tmp1;
            dAa_dVb = (aa - 1.0) * (dGg_dVb / gg - dT1s_dVb / t1s);
            inv_Aa = 1.0 / aa;
            vghigh = base.Param.B2vghigh;
            vglow = base.Param.B2vglow;
            if ((vgst >= vghigh) || (base.Param.B2n0 == 0.0))
            {
                vgeff = vgst;
                dVgeff_dVg = 1.0;
                dVgeff_dVd = -dVth_dVd;
                dVgeff_dVb = -dVth_dVb;
            }
            else
            {
                vof = base.Param.B2vof0 + base.Param.B2vofB * vbs + base.Param.B2vofD * vds;
                n = base.Param.B2n0 + base.Param.B2nB / t1s + base.Param.B2nD * vds;
                tmp = 0.5 / (n * ModelTemperature.Vtm);
                expArg1 = -vds / ModelTemperature.Vtm;
                expArg1 = Math.Max(expArg1, -30.0);
                exp1 = Math.Exp(expArg1);
                tmp1 = 1.0 - exp1;
                tmp1 = Math.Max(tmp1, 1.0e-18);
                tmp2 = 2.0 * aa * tmp1;
                if (vgst <= vglow)
                {
                    expArg = vgst * tmp;
                    expArg = Math.Max(expArg, -30.0);
                    exp0 = Math.Exp(0.5 * vof + expArg);
                    vgeff = Math.Sqrt(tmp2) * ModelTemperature.Vtm * exp0;
                    t0 = n * ModelTemperature.Vtm;
                    dVgeff_dVg = vgeff * tmp;
                    dVgeff_dVd = dVgeff_dVg * (n / tmp1 * exp1 - dVth_dVd - vgst * base.Param.B2nD / n + t0 * base.Param.B2vofD);
                    dVgeff_dVb = dVgeff_dVg * (base.Param.B2vofB * t0 - dVth_dVb + base.Param.B2nB * vgst / (n * t1s * t1s) * dT1s_dVb + t0 * inv_Aa * dAa_dVb);
                }
                else
                {
                    expArg = vglow * tmp;
                    expArg = Math.Max(expArg, -30.0);
                    exp0 = Math.Exp(0.5 * vof + expArg);
                    vgeff = Math.Sqrt(2.0 * aa * (1.0 - exp1)) * ModelTemperature.Vtm * exp0;
                    con1 = vghigh;
                    con3 = vgeff;
                    con4 = con3 * tmp;
                    sqrVghigh = vghigh * vghigh;
                    sqrVglow = vglow * vglow;
                    cubVghigh = vghigh * sqrVghigh;
                    cubVglow = vglow * sqrVglow;
                    t0 = 2.0 * vghigh;
                    t1 = 2.0 * vglow;
                    t2 = 3.0 * sqrVghigh;
                    t3 = 3.0 * sqrVglow;
                    t4 = vghigh - vglow;
                    t5 = sqrVghigh - sqrVglow;
                    t6 = cubVghigh - cubVglow;
                    t7 = con1 - con3;
                    delta = (t1 - t0) * t6 + (t2 - t3) * t5 + (t0 * t3 - t1 * t2) * t4;
                    delta = 1.0 / delta;
                    coeffb = (t1 - con4 * t0) * t6 + (con4 * t2 - t3) * t5 + (t0 * t3 - t1 * t2) * t7;
                    coeffc = (con4 - 1.0) * t6 + (t2 - t3) * t7 + (t3 - con4 * t2) * t4;
                    coeffd = (t1 - t0) * t7 + (1.0 - con4) * t5 + (con4 * t0 - t1) * t4;
                    coeffa = sqrVghigh * (coeffc + coeffd * t0);
                    vgeff = (coeffa + vgst * (coeffb + vgst * (coeffc + vgst * coeffd))) * delta;
                    dVgeff_dVg = (coeffb + vgst * (2.0 * coeffc + 3.0 * vgst * coeffd)) * delta;
                    t7 = con3 * tmp;
                    t8 = dT1s_dVb * base.Param.B2nB / (t1s * t1s * n);
                    t9 = n * ModelTemperature.Vtm;
                    dCon3_dVd = t7 * (n * exp1 / tmp1 - vglow * base.Param.B2nD / n + t9 * base.Param.B2vofD);
                    dCon3_dVb = t7 * (t9 * inv_Aa * dAa_dVb + vglow * t8 + t9 * base.Param.B2vofB);
                    dCon4_dVd = tmp * dCon3_dVd - t7 * base.Param.B2nD / n;
                    dCon4_dVb = tmp * dCon3_dVb + t7 * t8;
                    dCoeffb_dVd = dCon4_dVd * (t2 * t5 - t0 * t6) + dCon3_dVd * (t1 * t2 - t0 * t3);
                    dCoeffc_dVd = dCon4_dVd * (t6 - t2 * t4) + dCon3_dVd * (t3 - t2);
                    dCoeffd_dVd = dCon4_dVd * (t0 * t4 - t5) + dCon3_dVd * (t0 - t1);
                    dCoeffa_dVd = sqrVghigh * (dCoeffc_dVd + dCoeffd_dVd * t0);
                    dVgeff_dVd = -dVgeff_dVg * dVth_dVd + (dCoeffa_dVd + vgst * (dCoeffb_dVd + vgst * (dCoeffc_dVd + vgst * dCoeffd_dVd))) * delta;
                    dCoeffb_dVb = dCon4_dVb * (t2 * t5 - t0 * t6) + dCon3_dVb * (t1 * t2 - t0 * t3);
                    dCoeffc_dVb = dCon4_dVb * (t6 - t2 * t4) + dCon3_dVb * (t3 - t2);
                    dCoeffd_dVb = dCon4_dVb * (t0 * t4 - t5) + dCon3_dVb * (t0 - t1);
                    dCoeffa_dVb = sqrVghigh * (dCoeffc_dVb + dCoeffd_dVb * t0);
                    dVgeff_dVb = -dVgeff_dVg * dVth_dVb + (dCoeffa_dVb + vgst * (dCoeffb_dVb + vgst * (dCoeffc_dVb + vgst * dCoeffd_dVb))) * delta;
                }
            }
            if (vgeff > 0.0)
            {
                uvert = 1.0 + vgeff * (ua + vgeff * ub);
                uvert = Math.Max(uvert, 0.2);
                inv_Uvert = 1.0 / uvert;
                t8 = ua + 2.0 * ub * vgeff;
                dUvert_dVg = t8 * dVgeff_dVg;
                dUvert_dVd = t8 * dVgeff_dVd;
                dUvert_dVb = t8 * dVgeff_dVb + vgeff * (base.Param.B2uaB + vgeff * base.Param.B2ubB);
                t8 = u1s * inv_Aa * inv_Uvert;
                vc = t8 * vgeff;
                t9 = vc * inv_Uvert;
                dVc_dVg = t8 * dVgeff_dVg - t9 * dUvert_dVg;
                dVc_dVd = t8 * dVgeff_dVd - t9 * dUvert_dVd;
                dVc_dVb = t8 * dVgeff_dVb + base.Param.B2u1B * vgeff * inv_Aa * inv_Uvert - vc * inv_Aa * dAa_dVb - t9 * dUvert_dVb;
                tmp2 = Math.Sqrt(1.0 + 2.0 * vc);
                kk = 0.5 * (1.0 + vc + tmp2);
                inv_Kk = 1.0 / kk;
                dKk_dVc = 0.5 + 0.5 / tmp2;
                sqrtKk = Math.Sqrt(kk);
                t8 = inv_Aa / sqrtKk;
                vdsat = vgeff * t8;
                vdsat = Math.Max(vdsat, 1.0e-18);
                inv_Vdsat = 1.0 / vdsat;
                t9 = 0.5 * vdsat * inv_Kk * dKk_dVc;
                dVdsat_dVd = t8 * dVgeff_dVd - t9 * dVc_dVd;
                dVdsat_dVg = t8 * dVgeff_dVg - t9 * dVc_dVg;
                dVdsat_dVb = t8 * dVgeff_dVb - t9 * dVc_dVb - vdsat * inv_Aa * dAa_dVb;
                beta0 = base.Param.B2beta0 + base.Param.B2beta0B * vbs;
                betas = base.Param.B2betas0 + base.Param.B2betasB * vbs;
                beta2 = base.Param.B2beta20 + base.Param.B2beta2B * vbs + base.Param.B2beta2G * vgs;
                beta3 = base.Param.B2beta30 + base.Param.B2beta3B * vbs + base.Param.B2beta3G * vgs;
                beta4 = base.Param.B2beta40 + base.Param.B2beta4B * vbs + base.Param.B2beta4G * vgs;
                beta1 = betas - (beta0 + ModelParameters.Vdd * (beta3 - ModelParameters.Vdd * beta4));
                t0 = vds * beta2 * inv_Vdsat;
                t0 = Math.Min(t0, 30.0);
                t1 = Math.Exp(t0);
                t2 = t1 * t1;
                t3 = t2 + 1.0;
                tanh = (t2 - 1.0) / t3;
                sqrsech = 4.0 * t2 / (t3 * t3);
                beta = beta0 + beta1 * tanh + vds * (beta3 - beta4 * vds);
                t4 = beta1 * sqrsech * inv_Vdsat;
                t5 = ModelParameters.Vdd * tanh;
                dBeta_dVd = beta3 - 2.0 * beta4 * vds + t4 * (beta2 - t0 * dVdsat_dVd);
                dBeta_dVg = t4 * (base.Param.B2beta2G * vds - t0 * dVdsat_dVg) + base.Param.B2beta3G * (vds - t5) - base.Param.B2beta4G * (vds * vds - ModelParameters.Vdd * t5);
                dBeta1_dVb = base.Param.Arg;
                dBeta_dVb = base.Param.B2beta0B + dBeta1_dVb * tanh + vds * (base.Param.B2beta3B - vds * base.Param.B2beta4B) + t4 * (base.Param.B2beta2B * vds - t0 * dVdsat_dVb);
                if (vgst > vglow)
                {
                    if (vds <= vdsat)
                    {
                        t3 = vds * inv_Vdsat;
                        t4 = t3 - 1.0;
                        t2 = 1.0 - base.Param.B2u1D * t4 * t4;
                        u1 = u1s * t2;
                        utot = uvert + u1 * vds;
                        utot = Math.Max(utot, 0.5);
                        inv_Utot = 1.0 / utot;
                        t5 = 2.0 * u1s * base.Param.B2u1D * inv_Vdsat * t4;
                        dU1_dVd = t5 * (t3 * dVdsat_dVd - 1.0);
                        dU1_dVg = t5 * t3 * dVdsat_dVg;
                        dU1_dVb = t5 * t3 * dVdsat_dVb + base.Param.B2u1B * t2;
                        dUtot_dVd = dUvert_dVd + u1 + vds * dU1_dVd;
                        dUtot_dVg = dUvert_dVg + vds * dU1_dVg;
                        dUtot_dVb = dUvert_dVb + vds * dU1_dVb;
                        tmp1 = (vgeff - 0.5 * aa * vds);
                        tmp3 = tmp1 * vds;
                        betaeff = beta * inv_Utot;
                        Ids = betaeff * tmp3;
                        t6 = Ids / betaeff * inv_Utot;
                        gds = t6 * (dBeta_dVd - betaeff * dUtot_dVd) + betaeff * (tmp1 + (dVgeff_dVd - 0.5 * aa) * vds);
                        gm = t6 * (dBeta_dVg - betaeff * dUtot_dVg) + betaeff * vds * dVgeff_dVg;
                        gmb = t6 * (dBeta_dVb - betaeff * dUtot_dVb) + betaeff * vds * (dVgeff_dVb - 0.5 * vds * dAa_dVb);
                    }
                    else
                    {
                        tmp1 = vgeff * inv_Aa * inv_Kk;
                        tmp3 = 0.5 * vgeff * tmp1;
                        betaeff = beta * inv_Uvert;
                        Ids = betaeff * tmp3;
                        t0 = Ids / betaeff * inv_Uvert;
                        t1 = betaeff * vgeff * inv_Aa * inv_Kk;
                        t2 = Ids * inv_Kk * dKk_dVc;
                        if (base.Param.B2ai0 != 0.0)
                        {
                            ai = base.Param.B2ai0 + base.Param.B2aiB * vbs;
                            bi = base.Param.B2bi0 + base.Param.B2biB * vbs;
                            t5 = bi / (vds - vdsat);
                            t5 = Math.Min(t5, 30.0);
                            t6 = Math.Exp(-t5);
                            fR = 1.0 + ai * t6;
                            t7 = t5 / (vds - vdsat);
                            t8 = (1.0 - fR) * t7;
                            dFR_dVd = t8 * (dVdsat_dVd - 1.0);
                            dFR_dVg = t8 * dVdsat_dVg;
                            dFR_dVb = t8 * dVdsat_dVb + t6 * (base.Param.B2aiB - ai * base.Param.B2biB / (vds - vdsat));
                            gds = (t0 * (dBeta_dVd - betaeff * dUvert_dVd) + t1 * dVgeff_dVd - t2 * dVc_dVd) * fR + Ids * dFR_dVd;
                            gm = (t0 * (dBeta_dVg - betaeff * dUvert_dVg) + t1 * dVgeff_dVg - t2 * dVc_dVg) * fR + Ids * dFR_dVg;
                            gmb = (t0 * (dBeta_dVb - betaeff * dUvert_dVb) + t1 * dVgeff_dVb - t2 * dVc_dVb - Ids * inv_Aa * dAa_dVb) * fR + Ids * dFR_dVb;
                            Ids *= fR;
                        }
                        else
                        {
                            gds = t0 * (dBeta_dVd - betaeff * dUvert_dVd) + t1 * dVgeff_dVd - t2 * dVc_dVd;
                            gm = t0 * (dBeta_dVg - betaeff * dUvert_dVg) + t1 * dVgeff_dVg - t2 * dVc_dVg;
                            gmb = t0 * (dBeta_dVb - betaeff * dUvert_dVb) + t1 * dVgeff_dVb - t2 * dVc_dVb - Ids * inv_Aa * dAa_dVb;
                        }
                    }
                }
                else
                {
                    t0 = exp0 * exp0;
                    t1 = exp1;
                    Ids = beta * ModelTemperature.Vtm * ModelTemperature.Vtm * t0 * (1.0 - t1);
                    t2 = Ids / beta;
                    t4 = n * ModelTemperature.Vtm;
                    t3 = Ids / t4;
                    if ((vds > vdsat) && base.Param.B2ai0 != 0.0)
                    {
                        ai = base.Param.B2ai0 + base.Param.B2aiB * vbs;
                        bi = base.Param.B2bi0 + base.Param.B2biB * vbs;
                        t5 = bi / (vds - vdsat);
                        t5 = Math.Min(t5, 30.0);
                        t6 = Math.Exp(-t5);
                        fR = 1.0 + ai * t6;
                        t7 = t5 / (vds - vdsat);
                        t8 = (1.0 - fR) * t7;
                        dFR_dVd = t8 * (dVdsat_dVd - 1.0);
                        dFR_dVg = t8 * dVdsat_dVg;
                        dFR_dVb = t8 * dVdsat_dVb + t6 * (base.Param.B2aiB - ai * base.Param.B2biB / (vds - vdsat));
                    }
                    else
                    {
                        fR = 1.0;
                        dFR_dVd = 0.0;
                        dFR_dVg = 0.0;
                        dFR_dVb = 0.0;
                    }
                    gds = (t2 * dBeta_dVd + t3 * (base.Param.B2vofD * t4 - dVth_dVd - base.Param.B2nD * vgst / n) + beta * ModelTemperature.Vtm * t0 * t1) * fR + Ids * dFR_dVd;
                    gm = (t2 * dBeta_dVg + t3) * fR + Ids * dFR_dVg;
                    gmb = (t2 * dBeta_dVb + t3 * (base.Param.B2vofB * t4 - dVth_dVb + base.Param.B2nB * vgst / (n * t1s * t1s) * dT1s_dVb)) * fR + Ids * dFR_dVb;
                    Ids *= fR;
                }
            }
            else
            {
                Ids = 0.0;
                gm = 0.0;
                gds = 0.0;
                gmb = 0.0;
            }

            ChargeComputation:
            gds = Math.Max(gds, 1.0e-20);
            if ((ModelParameters.ChannelChargePartitionFlag > 1) || ((!chargeComputationNeeded) && (ModelParameters.ChannelChargePartitionFlag > -5)))
            {
                qg = 0.0;
                qd = 0.0;
                qb = 0.0;
                cgg = 0.0;
                cgs = 0.0;
                cgd = 0.0;
                cdg = 0.0;
                cds = 0.0;
                cdd = 0.0;
                cbg = 0.0;
                cbs = 0.0;
                cbd = 0.0;
                goto finished;
            }
            else
            {
                if (vbs < 0.0)
                {
                    vbseff = vbs;
                    dVbseff_dVb = 1.0;
                }
                else
                {
                    vbseff = base.Param.B2phi - phisb;
                    dVbseff_dVb = -dPhisb_dVb;
                }
                arg1 = vgs - vbseff - base.Param.B2vfb;
                arg2 = arg1 - vgst;
                qbulk = base.Param.One_Third_CoxWL * arg2;
                dQbulk_dVb = base.Param.One_Third_CoxWL * (dVth_dVb - dVbseff_dVb);
                dQbulk_dVd = base.Param.One_Third_CoxWL * dVth_dVd;
                if (arg1 <= 0.0)
                {
                    qg = base.Param.CoxWL * arg1;
                    qb = -(qg);
                    qd = 0.0;
                    cgg = base.Param.CoxWL;
                    cgd = 0.0;
                    cgs = -cgg * (1.0 - dVbseff_dVb);
                    cdg = 0.0;
                    cdd = 0.0;
                    cds = 0.0;
                    cbg = -base.Param.CoxWL;
                    cbd = 0.0;
                    cbs = -cgs;
                }
                else if (vgst <= 0.0)
                {
                    t2 = arg1 / arg2;
                    t3 = t2 * t2 * (base.Param.CoxWL - base.Param.Two_Third_CoxWL * t2);
                    qg = base.Param.CoxWL * arg1 * (1.0 - t2 * (1.0 - t2 / 3.0));
                    qb = -(qg);
                    qd = 0.0;
                    cgg = base.Param.CoxWL * (1.0 - t2 * (2.0 - t2));
                    tmp = t3 * dVth_dVb - (cgg + t3) * dVbseff_dVb;
                    cgd = t3 * dVth_dVd;
                    cgs = -(cgg + cgd + tmp);
                    cdg = 0.0;
                    cdd = 0.0;
                    cds = 0.0;
                    cbg = -cgg;
                    cbd = -cgd;
                    cbs = -cgs;
                }
                else
                {
                    if (vgst < base.Param.B2vghigh)
                    {
                        uvert = 1.0 + vgst * (ua + vgst * ub);
                        uvert = Math.Max(uvert, 0.2);
                        inv_Uvert = 1.0 / uvert;
                        dUvert_dVg = ua + 2.0 * ub * vgst;
                        dUvert_dVd = -dUvert_dVg * dVth_dVd;
                        dUvert_dVb = -dUvert_dVg * dVth_dVb + vgst * (base.Param.B2uaB + vgst * base.Param.B2ubB);
                        t8 = u1s * inv_Aa * inv_Uvert;
                        vc = t8 * vgst;
                        t9 = vc * inv_Uvert;
                        dVc_dVg = t8 - t9 * dUvert_dVg;
                        dVc_dVd = -t8 * dVth_dVd - t9 * dUvert_dVd;
                        dVc_dVb = -t8 * dVth_dVb + base.Param.B2u1B * vgst * inv_Aa * inv_Uvert - vc * inv_Aa * dAa_dVb - t9 * dUvert_dVb;
                        tmp2 = Math.Sqrt(1.0 + 2.0 * vc);
                        kk = 0.5 * (1.0 + vc + tmp2);
                        inv_Kk = 1.0 / kk;
                        dKk_dVc = 0.5 + 0.5 / tmp2;
                        sqrtKk = Math.Sqrt(kk);
                        t8 = inv_Aa / sqrtKk;
                        vdsat = vgst * t8;
                        t9 = 0.5 * vdsat * inv_Kk * dKk_dVc;
                        dVdsat_dVd = -t8 * dVth_dVd - t9 * dVc_dVd;
                        dVdsat_dVg = t8 - t9 * dVc_dVg;
                        dVdsat_dVb = -t8 * dVth_dVb - t9 * dVc_dVb - vdsat * inv_Aa * dAa_dVb;
                    }
                    if (vds >= vdsat)
                    {
                        cgg = base.Param.Two_Third_CoxWL;
                        cgd = -cgg * dVth_dVd + dQbulk_dVd;
                        tmp = -cgg * dVth_dVb + dQbulk_dVb;
                        cgs = -(cgg + cgd + tmp);
                        cbg = 0.0;
                        cbd = -dQbulk_dVd;
                        cbs = dQbulk_dVd + dQbulk_dVb;
                        cdg = -0.4 * cgg;
                        tmp = -cdg * dVth_dVb;
                        cdd = -cdg * dVth_dVd;
                        cds = -(cdg + cdd + tmp);
                        qb = -qbulk;
                        qg = base.Param.Two_Third_CoxWL * vgst + qbulk;
                        qd = cdg * vgst;
                    }
                    else
                    {
                        t7 = vds / vdsat;
                        t8 = vgst / vdsat;
                        t6 = t7 * t8;
                        t9 = 1.0 - t7;
                        vgdt = vgst * t9;
                        t0 = vgst / (vgst + vgdt);
                        t1 = vgdt / (vgst + vgdt);
                        t5 = t0 * t1;
                        t2 = 1.0 - t1 + t5;
                        t3 = 1.0 - t0 + t5;
                        dVgdt_dVg = t9 + t6 * dVdsat_dVg;
                        dVgdt_dVd = t6 * dVdsat_dVd - t8 - t9 * dVth_dVd;
                        dVgdt_dVb = t6 * dVdsat_dVb - t9 * dVth_dVb;
                        qg = base.Param.Two_Third_CoxWL * (vgst + vgdt - vgdt * t0) + qbulk;
                        qb = -qbulk;
                        qd = -base.Param.One_Third_CoxWL * (0.2 * vgdt + 0.8 * vgst + vgdt * t1 + 0.2 * t5 * (vgdt - vgst));
                        cgg = base.Param.Two_Third_CoxWL * (t2 + t3 * dVgdt_dVg);
                        tmp = dQbulk_dVb + base.Param.Two_Third_CoxWL * (t3 * dVgdt_dVb - t2 * dVth_dVb);
                        cgd = base.Param.Two_Third_CoxWL * (t3 * dVgdt_dVd - t2 * dVth_dVd) + dQbulk_dVd;
                        cgs = -(cgg + cgd + tmp);
                        t2 = 0.8 - 0.4 * t1 * (2.0 * t1 + t0 + t0 * (t1 - t0));
                        t3 = 0.2 + t1 + t0 * (1.0 - 0.4 * t0 * (t1 + 3.0 * t0));
                        cdg = -base.Param.One_Third_CoxWL * (t2 + t3 * dVgdt_dVg);
                        tmp = base.Param.One_Third_CoxWL * (t2 * dVth_dVb - t3 * dVgdt_dVb);
                        cdd = base.Param.One_Third_CoxWL * (t2 * dVth_dVd - t3 * dVgdt_dVd);
                        cds = -(cdg + tmp + cdd);
                        cbg = 0.0;
                        cbd = -dQbulk_dVd;
                        cbs = dQbulk_dVd + dQbulk_dVb;
                    }
                }
            }
            finished:
            switch (-ModelParameters.ChannelChargePartitionFlag)
            {
                case 0:
                    Ids = Math.Max(Ids, 1e-50);
                    break;
                case 1:
                    Ids = Math.Max(Ids, 1e-50);
                    break;
                case 2:
                    Ids = gm;
                    break;
                case 3:
                    Ids = gds;
                    break;
                case 4:
                    Ids = 1.0 / gds;
                    break;
                case 5:
                    Ids = gmb;
                    break;
                case 6:
                    Ids = qg / 1.0e-12;
                    break;
                case 7:
                    Ids = qb / 1.0e-12;
                    break;
                case 8:
                    Ids = qd / 1.0e-12;
                    break;
                case 9:
                    Ids = -(qb + qg + qd) / 1.0e-12;
                    break;
                case 10:
                    Ids = cgg / 1.0e-12;
                    break;
                case 11:
                    Ids = cgd / 1.0e-12;
                    break;
                case 12:
                    Ids = cgs / 1.0e-12;
                    break;
                case 13:
                    Ids = -(cgg + cgd + cgs) / 1.0e-12;
                    break;
                case 14:
                    Ids = cbg / 1.0e-12;
                    break;
                case 15:
                    Ids = cbd / 1.0e-12;
                    break;
                case 16:
                    Ids = cbs / 1.0e-12;
                    break;
                case 17:
                    Ids = -(cbg + cbd + cbs) / 1.0e-12;
                    break;
                case 18:
                    Ids = cdg / 1.0e-12;
                    break;
                case 19:
                    Ids = cdd / 1.0e-12;
                    break;
                case 20:
                    Ids = cds / 1.0e-12;
                    break;
                case 21:
                    Ids = -(cdg + cdd + cds) / 1.0e-12;
                    break;
                case 22:
                    Ids = -(cgg + cdg + cbg) / 1.0e-12;
                    break;
                case 23:
                    Ids = -(cgd + cdd + cbd) / 1.0e-12;
                    break;
                case 24:
                    Ids = -(cgs + cds + cbs) / 1.0e-12;
                    break;
                default:
                    Ids = Math.Max(Ids, 1.0e-50);
                    break;
            }
            von = vth;
            Vdsat = vdsat;
        }

	    /// <summary>
	    /// Helper method MosCap
	    /// </summary>
	    public void MosCap(double vgd, double vgs, double vgb, double[] args, double cbgb, double cbdb, double cbsb, double cdgb, double cddb, double cdsb, out double gcggbPointer, out double gcgdbPointer, out double gcgsbPointer, out double gcbgbPointer, out double gcbdbPointer, out double gcbsbPointer, out double gcdgbPointer, out double gcddbPointer, out double gcdsbPointer, out double gcsgbPointer, out double gcsdbPointer, out double gcssbPointer, ref double qGatePointer, ref double qBulkPointer, ref double qDrainPointer, out double qSourcePointer)
	    {
	        double qgd, qgs, qgb, ag0;
	        ag0 = TranBehavior?.Qg.Jacobian(1.0) ?? 0.0; // ckt.CKTag[0];
	        gcdgbPointer = (cdgb - args[0]) * ag0;
	        gcddbPointer = (cddb + args[3] + args[0]) * ag0;
	        gcdsbPointer = cdsb * ag0;
	        gcsgbPointer = -(args[5] + cbgb + cdgb + args[1]) * ag0;
	        gcsdbPointer = -(args[6] + cbdb + cddb) * ag0;
	        gcssbPointer = (args[4] + args[1] - (args[7] + cbsb + cdsb)) * ag0;
	        gcggbPointer = (args[5] + args[0] + args[1] + args[2]) * ag0;
	        gcgdbPointer = (args[6] - args[0]) * ag0;
	        gcgsbPointer = (args[7] - args[1]) * ag0;
	        gcbgbPointer = (cbgb - args[2]) * ag0;
	        gcbdbPointer = (cbdb - args[3]) * ag0;
	        gcbsbPointer = (cbsb - args[4]) * ag0;
	        qgd = args[0] * vgd;
	        qgs = args[1] * vgs;
	        qgb = args[2] * vgb;
	        qGatePointer = qGatePointer + qgd + qgs + qgb;
	        qBulkPointer = qBulkPointer - qgb;
	        qDrainPointer = qDrainPointer - qgd;
	        qSourcePointer = -(qGatePointer + qBulkPointer + qDrainPointer);
	    }
    }
}