namespace SpiceSharp.Components.BSIM2Behaviors
{
    /// <summary>
    /// Size-dependent parameters
    /// </summary>
    public class BSIM2SizeDependParams
    {
        public double Width;
        public double Length;

        // flat band voltage at given L and W 
        public double B2vfb;

        // surface potential at strong inversion 
        public double B2phi;

        // bulk effect coefficient 1             
        public double B2k1;

        // bulk effect coefficient 2             
        public double B2k2;

        // drain induced barrier lowering        
        public double B2eta0;

        // Vbs dependence of Eta                 
        public double B2etaB;

        // Beta at Vds = 0 and Vgs = Vth         
        public double B2beta0;

        // Vbs dependence of Beta0               
        public double B2beta0B;

        // Beta at Vds=Vdd and Vgs=Vth           
        public double B2betas0;

        // Vbs dependence of Betas               
        public double B2betasB;

        // Vds dependence of Beta in tanh term   
        public double B2beta20;

        // Vbs dependence of Beta2               
        public double B2beta2B;

        // Vgs dependence of Beta2               
        public double B2beta2G;

        // Vds dependence of Beta in linear term 
        public double B2beta30;

        // Vbs dependence of Beta3               
        public double B2beta3B;

        // Vgs dependence of Beta3               
        public double B2beta3G;

        // Vds dependence of Beta in quadra term 
        public double B2beta40;

        // Vbs dependence of Beta4               
        public double B2beta4B;

        // Vgs dependence of Beta4               
        public double B2beta4G;

        // Linear Vgs dependence of Mobility     
        public double B2ua0;

        // Vbs dependence of Ua                  
        public double B2uaB;

        // Quadratic Vgs dependence of Mobility  
        public double B2ub0;

        // Vbs dependence of Ub                  
        public double B2ubB;

        // Drift Velocity Saturation due to Vds  
        public double B2u10;

        // Vbs dependence of U1                  
        public double B2u1B;

        // Vds dependence of U1                  
        public double B2u1D;

        // Subthreshold slope at Vds=0, Vbs=0    
        public double B2n0;

        // Vbs dependence of n                   
        public double B2nB;

        // Vds dependence of n                   
        public double B2nD;

        // Vth offset at Vds=0, Vbs=0            
        public double B2vof0;

        // Vbs dependence of Vof                 
        public double B2vofB;

        // Vds dependence of Vof                 
        public double B2vofD;

        // Pre-factor in hot-electron effects    
        public double B2ai0;

        // Vbs dependence of Ai                  
        public double B2aiB;

        // Exp-factor in hot-electron effects    
        public double B2bi0;

        // Vbs dependence of Bi                  
        public double B2biB;

        // Upper bound of cubic spline function  
        public double B2vghigh;

        // Lower bound of cubic spline function  
        public double B2vglow;

        // Gate Drain Overlap Capacitance
        public double B2GDoverlapCap;

        // Gate Source Overlap Capacitance
        public double B2GSoverlapCap;

        // Gate Bulk Overlap Capacitance
        public double B2GBoverlapCap;

        public double SqrtPhi;
        public double Phis3;
        public double CoxWL;
        public double One_Third_CoxWL;
        public double Two_Third_CoxWL;
        public double Arg;
        public double B2vt0;
    }
}
