﻿using System;
using System.Collections.Generic;
using System.Text;

namespace SpiceSharpBSIM.Components.Semiconductors.BSIM.BSIM4Behaviors
{
    public class SizeDependentProperties
    {
        public double Width;
        public double Length;
        public double NFinger;
        public double BSIM4cdsc;
        public double BSIM4cdscb;
        public double BSIM4cdscd;
        public double BSIM4cit;
        public double BSIM4nfactor;
        public double BSIM4xj;
        public double BSIM4vsat;
        public double BSIM4at;
        public double BSIM4a0;
        public double BSIM4ags;
        public double BSIM4a1;
        public double BSIM4a2;
        public double BSIM4keta;
        public double BSIM4nsub;
        public double BSIM4ndep;
        public double BSIM4nsd;
        public double BSIM4phin;
        public double BSIM4ngate;
        public double BSIM4gamma1;
        public double BSIM4gamma2;
        public double BSIM4vbx;
        public double BSIM4vbi;
        public double BSIM4vbm;
        public double BSIM4xt;
        public double BSIM4phi;
        public double BSIM4litl;
        public double BSIM4k1;
        public double BSIM4kt1;
        public double BSIM4kt1l;
        public double BSIM4kt2;
        public double BSIM4k2;
        public double BSIM4k3;
        public double BSIM4k3b;
        public double BSIM4w0;
        public double BSIM4dvtp0;
        public double BSIM4dvtp1;
        public double BSIM4dvtp2;  /* New DIBL/Rout */
        public double BSIM4dvtp3;
        public double BSIM4dvtp4;
        public double BSIM4dvtp5;
        public double BSIM4lpe0;
        public double BSIM4lpeb;
        public double BSIM4dvt0;
        public double BSIM4dvt1;
        public double BSIM4dvt2;
        public double BSIM4dvt0w;
        public double BSIM4dvt1w;
        public double BSIM4dvt2w;
        public double BSIM4drout;
        public double BSIM4dsub;
        public double BSIM4vth0;
        public double BSIM4ua;
        public double BSIM4ua1;
        public double BSIM4ub;
        public double BSIM4ub1;
        public double BSIM4uc;
        public double BSIM4uc1;
        public double BSIM4ud;
        public double BSIM4ud1;
        public double BSIM4up;
        public double BSIM4lp;
        public double BSIM4u0;
        public double BSIM4eu;
        public double BSIM4ucs;
        public double BSIM4ute;
        public double BSIM4ucste;
        public double BSIM4voff;
        public double BSIM4tvoff;
        public double BSIM4tnfactor;   /* v4.7 Temp dep of leakage current */
        public double BSIM4teta0;      /* v4.7 temp dep of leakage current */
        public double BSIM4tvoffcv;    /* v4.7 temp dep of leakage current */
        public double BSIM4minv;
        public double BSIM4minvcv;
        public double BSIM4vfb;
        public double BSIM4delta;
        public double BSIM4rdsw;
        public double BSIM4rds0;
        public double BSIM4rs0;
        public double BSIM4rd0;
        public double BSIM4rsw;
        public double BSIM4rdw;
        public double BSIM4prwg;
        public double BSIM4prwb;
        public double BSIM4prt;
        public double BSIM4eta0;
        public double BSIM4etab;
        public double BSIM4pclm;
        public double BSIM4pdibl1;
        public double BSIM4pdibl2;
        public double BSIM4pdiblb;
        public double BSIM4fprout;
        public double BSIM4pdits;
        public double BSIM4pditsd;
        public double BSIM4pscbe1;
        public double BSIM4pscbe2;
        public double BSIM4pvag;
        public double BSIM4wr;
        public double BSIM4dwg;
        public double BSIM4dwb;
        public double BSIM4b0;
        public double BSIM4b1;
        public double BSIM4alpha0;
        public double BSIM4alpha1;
        public double BSIM4beta0;
        public double BSIM4agidl;
        public double BSIM4bgidl;
        public double BSIM4cgidl;
        public double BSIM4egidl;
        public double BSIM4fgidl; /* v4.7 New GIDL/GISL */
        public double BSIM4kgidl; /* v4.7 New GIDL/GISL */
        public double BSIM4rgidl; /* v4.7 New GIDL/GISL */
        public double BSIM4agisl;
        public double BSIM4bgisl;
        public double BSIM4cgisl;
        public double BSIM4egisl;
        public double BSIM4fgisl; /* v4.7 New GIDL/GISL */
        public double BSIM4kgisl; /* v4.7 New GIDL/GISL */
        public double BSIM4rgisl; /* v4.7 New GIDL/GISL */
        public double BSIM4aigc;
        public double BSIM4bigc;
        public double BSIM4cigc;
        public double BSIM4aigsd;
        public double BSIM4bigsd;
        public double BSIM4cigsd;
        public double BSIM4aigs;
        public double BSIM4bigs;
        public double BSIM4cigs;
        public double BSIM4aigd;
        public double BSIM4bigd;
        public double BSIM4cigd;
        public double BSIM4aigbacc;
        public double BSIM4bigbacc;
        public double BSIM4cigbacc;
        public double BSIM4aigbinv;
        public double BSIM4bigbinv;
        public double BSIM4cigbinv;
        public double BSIM4nigc;
        public double BSIM4nigbacc;
        public double BSIM4nigbinv;
        public double BSIM4ntox;
        public double BSIM4eigbinv;
        public double BSIM4pigcd;
        public double BSIM4poxedge;
        public double BSIM4xrcrg1;
        public double BSIM4xrcrg2;
        public double BSIM4lambda; /* overshoot */
        public double BSIM4vtl; /* thermal velocity limit */
        public double BSIM4xn; /* back scattering parameter */
        public double BSIM4lc; /* back scattering parameter */
        public double BSIM4tfactor;  /* ballistic transportation factor  */
        public double BSIM4vfbsdoff;  /* S/D flatband offset voltage  */
        public double BSIM4tvfbsdoff;

        /* added for stress effect */
        public double BSIM4ku0;
        public double BSIM4kvth0;
        public double BSIM4ku0temp;
        public double BSIM4rho_ref;
        public double BSIM4inv_od_ref;
        /* added for well proximity effect */
        public double BSIM4kvth0we;
        public double BSIM4k2we;
        public double BSIM4ku0we;

        /* CV model */
        public double BSIM4cgsl;
        public double BSIM4cgdl;
        public double BSIM4ckappas;
        public double BSIM4ckappad;
        public double BSIM4cf;
        public double BSIM4clc;
        public double BSIM4cle;
        public double BSIM4vfbcv;
        public double BSIM4noff;
        public double BSIM4voffcv;
        public double BSIM4acde;
        public double BSIM4moin;

        /* Pre-calculated constants */

        public double BSIM4dw;
        public double BSIM4dl;
        public double BSIM4leff;
        public double BSIM4weff;

        public double BSIM4dwc;
        public double BSIM4dlc;
        public double BSIM4dwj;
        public double BSIM4leffCV;
        public double BSIM4weffCV;
        public double BSIM4weffCJ;
        public double BSIM4abulkCVfactor;
        public double BSIM4cgso;
        public double BSIM4cgdo;
        public double BSIM4cgbo;

        public double BSIM4u0temp;
        public double BSIM4vsattemp;
        public double BSIM4sqrtPhi;
        public double BSIM4phis3;
        public double BSIM4Xdep0;
        public double BSIM4sqrtXdep0;
        public double BSIM4theta0vb0;
        public double BSIM4thetaRout;
        public double BSIM4mstar;
        public double BSIM4VgsteffVth;
        public double BSIM4mstarcv;
        public double BSIM4voffcbn;
        public double BSIM4voffcbncv;
        public double BSIM4rdswmin;
        public double BSIM4rdwmin;
        public double BSIM4rswmin;
        public double BSIM4vfbsd;

        public double BSIM4cof1;
        public double BSIM4cof2;
        public double BSIM4cof3;
        public double BSIM4cof4;
        public double BSIM4cdep0;
        public double BSIM4ToxRatio;
        public double BSIM4Aechvb;
        public double BSIM4Bechvb;
        public double BSIM4ToxRatioEdge;
        public double BSIM4AechvbEdgeS;
        public double BSIM4AechvbEdgeD;
        public double BSIM4BechvbEdge;
        public double BSIM4ldeb;
        public double BSIM4k1ox;
        public double BSIM4k2ox;
        public double BSIM4vfbzbfactor;
        public double BSIM4dvtp2factor; /* v4.7 */
    }
}
