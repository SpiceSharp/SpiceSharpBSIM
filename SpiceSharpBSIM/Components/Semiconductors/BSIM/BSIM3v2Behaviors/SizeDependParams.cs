using System;
using System.Collections.Generic;
using System.Text;
using SpiceSharp.Components;

namespace SpiceSharpBSIM.Components.Semiconductors.BSIM.BSIM3v2Behaviors
{
    /// <summary>
    /// Size-dependent parameters for <see cref="BSIM3v2Model"/>
    /// </summary>
    public class SizeDependParams
    {
        public double Width { get; set; }
        public double Length { get; set; }
        public double BSIM3v32cdsc { get; set; }
        public double BSIM3v32cdscb { get; set; }
        public double BSIM3v32cdscd { get; set; }
        public double BSIM3v32cit { get; set; }
        public double BSIM3v32nfactor { get; set; }
        public double BSIM3v32xj { get; set; }
        public double BSIM3v32vsat { get; set; }
        public double BSIM3v32at { get; set; }
        public double BSIM3v32a0 { get; set; }
        public double BSIM3v32ags { get; set; }
        public double BSIM3v32a1 { get; set; }
        public double BSIM3v32a2 { get; set; }
        public double BSIM3v32keta { get; set; }
        public double BSIM3v32nsub { get; set; }
        public double BSIM3v32npeak { get; set; }
        public double BSIM3v32ngate { get; set; }
        public double BSIM3v32gamma1 { get; set; }
        public double BSIM3v32gamma2 { get; set; }
        public double BSIM3v32vbx { get; set; }
        public double BSIM3v32vbi { get; set; }
        public double BSIM3v32vbm { get; set; }
        public double BSIM3v32vbsc { get; set; }
        public double BSIM3v32xt { get; set; }
        public double BSIM3v32phi { get; set; }
        public double BSIM3v32litl { get; set; }
        public double BSIM3v32k1 { get; set; }
        public double BSIM3v32kt1 { get; set; }
        public double BSIM3v32kt1l { get; set; }
        public double BSIM3v32kt2 { get; set; }
        public double BSIM3v32k2 { get; set; }
        public double BSIM3v32k3 { get; set; }
        public double BSIM3v32k3b { get; set; }
        public double BSIM3v32w0 { get; set; }
        public double BSIM3v32nlx { get; set; }
        public double BSIM3v32dvt0 { get; set; }
        public double BSIM3v32dvt1 { get; set; }
        public double BSIM3v32dvt2 { get; set; }
        public double BSIM3v32dvt0w { get; set; }
        public double BSIM3v32dvt1w { get; set; }
        public double BSIM3v32dvt2w { get; set; }
        public double BSIM3v32drout { get; set; }
        public double BSIM3v32dsub { get; set; }
        public double BSIM3v32vth0 { get; set; }
        public double BSIM3v32ua { get; set; }
        public double BSIM3v32ua1 { get; set; }
        public double BSIM3v32ub { get; set; }
        public double BSIM3v32ub1 { get; set; }
        public double BSIM3v32uc { get; set; }
        public double BSIM3v32uc1 { get; set; }
        public double BSIM3v32u0 { get; set; }
        public double BSIM3v32ute { get; set; }
        public double BSIM3v32voff { get; set; }
        public double BSIM3v32vfb { get; set; }
        public double BSIM3v32delta { get; set; }
        public double BSIM3v32rdsw { get; set; }
        public double BSIM3v32rds0 { get; set; }
        public double BSIM3v32prwg { get; set; }
        public double BSIM3v32prwb { get; set; }
        public double BSIM3v32prt { get; set; }
        public double BSIM3v32eta0 { get; set; }
        public double BSIM3v32etab { get; set; }
        public double BSIM3v32pclm { get; set; }
        public double BSIM3v32pdibl1 { get; set; }
        public double BSIM3v32pdibl2 { get; set; }
        public double BSIM3v32pdiblb { get; set; }
        public double BSIM3v32pscbe1 { get; set; }
        public double BSIM3v32pscbe2 { get; set; }
        public double BSIM3v32pvag { get; set; }
        public double BSIM3v32wr { get; set; }
        public double BSIM3v32dwg { get; set; }
        public double BSIM3v32dwb { get; set; }
        public double BSIM3v32b0 { get; set; }
        public double BSIM3v32b1 { get; set; }
        public double BSIM3v32alpha0 { get; set; }
        public double BSIM3v32alpha1 { get; set; }
        public double BSIM3v32beta0 { get; set; }
        public double BSIM3v32elm { get; set; }
        public double BSIM3v32cgsl { get; set; }
        public double BSIM3v32cgdl { get; set; }
        public double BSIM3v32ckappa { get; set; }
        public double BSIM3v32cf { get; set; }
        public double BSIM3v32clc { get; set; }
        public double BSIM3v32cle { get; set; }
        public double BSIM3v32vfbcv { get; set; }
        public double BSIM3v32noff { get; set; }
        public double BSIM3v32voffcv { get; set; }
        public double BSIM3v32acde { get; set; }
        public double BSIM3v32moin { get; set; }
        public double BSIM3v32dw { get; set; }
        public double BSIM3v32dl { get; set; }
        public double BSIM3v32leff { get; set; }
        public double BSIM3v32weff { get; set; }
        public double BSIM3v32dwc { get; set; }
        public double BSIM3v32dlc { get; set; }
        public double BSIM3v32leffCV { get; set; }
        public double BSIM3v32weffCV { get; set; }
        public double BSIM3v32abulkCVfactor { get; set; }
        public double BSIM3v32cgso { get; set; }
        public double BSIM3v32cgdo { get; set; }
        public double BSIM3v32cgbo { get; set; }
        public double BSIM3v32tconst { get; set; }
        public double BSIM3v32u0temp { get; set; }
        public double BSIM3v32vsattemp { get; set; }
        public double BSIM3v32sqrtPhi { get; set; }
        public double BSIM3v32phis3 { get; set; }
        public double BSIM3v32Xdep0 { get; set; }
        public double BSIM3v32sqrtXdep0 { get; set; }
        public double BSIM3v32theta0vb0 { get; set; }
        public double BSIM3v32thetaRout { get; set; }
        public double BSIM3v32cof1 { get; set; }
        public double BSIM3v32cof2 { get; set; }
        public double BSIM3v32cof3 { get; set; }
        public double BSIM3v32cof4 { get; set; }
        public double BSIM3v32cdep0 { get; set; }
        public double BSIM3v32vfbzb { get; set; }
        public double BSIM3v32ldeb { get; set; }
        public double BSIM3v32k1ox { get; set; }
        public double BSIM3v32k2ox { get; set; }
    }
}
