using NUnit.Framework;
using SpiceSharp.Components;
using System.IO;
using SpiceSharp;
using SpiceSharp.Simulations;
using System.Numerics;

namespace SpiceSharpBSIMTests.BSIM4Tests;

public class BSIM4Tests : Framework
{
    private BSIM4 CreateMosfet(string name, string drain, string gate, string source, string bulk, double w, double l, string model)
    {
        var e = new BSIM4(name, drain, gate, source, bulk);
        e.SetParameter("w", w);
        e.SetParameter("l", l);
        e.Model = model;
        return e;
    }

    private BSIM4Model CreateModel(string name, string parameters)
    {
        var m = new BSIM4Model(name);
        m.Parameters.CheckPath = Path.Combine(TestContext.CurrentContext.WorkDirectory, m.Parameters.CheckPath);
        ApplyParameters(m, parameters);
        return m;
    }

    [Test]
    public void When_BSIM4DC_Expect_Reference()
    {
        var ckt = new Circuit(
            new VoltageSource("V1", "g", "0", 0.0),
            new VoltageSource("V2", "d", "0", 0.0),
            CreateMosfet("M1", "d", "g", "0", "0", 2e-6, 100e-9, "mod"),
            CreateModel("mod", "binunit=1 paramchk=1 mobmod=0 capmod=2 igcmod=1 igbmod=1 geomod=1 diomod=1 rdsmod=0 rbodymod=1 rgatemod=1 permod=1 acnqsmod=0 trnqsmod=0 tnom=27 toxe=1.8e-009 toxp=1.5e-009 toxm=1.8e-009 dtox=3e-010 epsrox=3.9 wint=5e-009 lint=0 ll=0 wl=0 lln=1 wln=1 lw=0 ww=0 lwn=1 wwn=1 lwl=0 wwl=0 xpart=0 toxref=1.8e-009 vth0=0.62261 k1=0.4 k2=0 k3=0 k3b=0 w0=2.5e-006 dvt0=1 dvt1=2 dvt2=0 dvt0w=0 dvt1w=0 dvt2w=0 dsub=0.1 minv=0.05 voffl=0 dvtp0=1e-010 dvtp1=0.1 lpe0=0 lpeb=0 xj=1.4e-008 ngate=1e023 ndep=3.24e018 nsd=2e020 phin=0 cdsc=0 cdscb=0 cdscd=0 cit=0 voff=-0.13 nfactor=1.6 eta0=0.0125 etab=0 vfb=-0.55 u0=0.049 ua=6e-010 ub=1.2e-018 uc=0 vsat=130000 a0=1 ags=0 a1=0 a2=1 b0=0 b1=0 keta=0.04 dwg=0 dwb=0 pclm=0.02 pdiblc1=0.001 pdiblc2=0.001 pdiblcb=-0.005 drout=0.5 pvag=1e-020 delta=0.01 pscbe1=8.14e008 pscbe2=1e-007 fprout=0.2 pdits=0.08 pditsd=0.23 pditsl=2300000 rsh=5 rdsw=210 rsw=80 rdw=80 rdswmin=0 rdwmin=0 rswmin=0 prwg=0 prwb=0 wr=1 alpha0=0.074 alpha1=0.005 beta0=30 agidl=0.0002 bgidl=2.1e009 cgidl=0.0002 egidl=0.8 aigbacc=0.012 bigbacc=0.0028 cigbacc=0.002 nigbacc=1 aigbinv=0.014 bigbinv=0.004 cigbinv=0.004 eigbinv=1.1 nigbinv=3 aigc=0.015211 bigc=0.0027432 cigc=0.002 aigsd=0.015211 bigsd=0.0027432 cigsd=0.002 nigc=1 poxedge=1 pigcd=1 ntox=1 xrcrg1=12 xrcrg2=5 cgso=1.1e-010 cgdo=1.1e-010 cgbo=2.56e-011 cgdl=2.653e-010 cgsl=2.653e-010 ckappas=0.03 ckappad=0.03 acde=1 moin=15 noff=0.9 voffcv=0.02 kt1=-0.11 kt1l=0 kt2=0.022 ute=-1.5 ua1=4.31e-009 ub1=7.61e-018 uc1=-5.6e-011 prt=0 at=33000 fnoimod=1 tnoimod=0 jss=0.0001 jsws=1e-011 jswgs=1e-010 njs=1 ijthsfwd=0.01 ijthsrev=0.001 bvs=10 xjbvs=1 jsd=0.0001 jswd=1e-011 jswgd=1e-010 njd=1 ijthdfwd=0.01 ijthdrev=0.001 bvd=10 xjbvd=1 pbs=1 cjs=0.0005 mjs=0.5 pbsws=1 cjsws=5e-010 mjsws=0.33 pbswgs=1 cjswgs=3e-010 mjswgs=0.33 pbd=1 cjd=0.0005 mjd=0.5 pbswd=1 cjswd=5e-010 mjswd=0.33 pbswgd=1 cjswgd=5e-010 mjswgd=0.33 tpb=0.005 tcj=0.001 tpbsw=0.005 tcjsw=0.001 tpbswg=0.005 tcjswg=0.001 xtis=3 xtid=3 dmcg=0 dmci=0 dmdg=0 dmcgt=0 dwj=0 xgw=0 xgl=0 rshg=0.4 gbmin=1e-010 rbpb=5 rbpd=15 rbps=15 rbdb=15 rbsb=15 ngcon=1"));
        ckt["M1"].SetParameter("m", 4.0);

        // Create simulation
        var dc = new DC("dc", new[]
        {
            new ParameterSweep("V1", new LinearSweep(0, 1, 0.2)),
            new ParameterSweep("V2", new LinearSweep(0, 1, 0.2))
        });

        // Create exports
        var exports = new IExport<double>[] { new RealPropertyExport(dc, "V2", "i") };

        // Create references
        // 20220724 - Sven Boulanger - using ngSpice
        double[][] references =
        {
            new[]
            {
                -1.867172903871107e-57, -1.518642094248029e-11, -1.725534315605076e-11, -1.936879123959608e-11, -2.157333104510600e-11, -2.388963914948112e-11, 1.031224154996069e-19, -2.866382476494586e-09, -3.118541441210817e-09, -3.379289861110353e-09, -3.657866462301296e-09, -3.958343148212254e-09, 5.440260251875172e-16, -4.882015011054686e-07, -5.285692882330689e-07, -5.699228343063414e-07, -6.137916853345052e-07, -6.607994011262966e-07, 1.929548160065753e-12, -3.775429363181440e-05, -4.028266507593196e-05, -4.271941532569848e-05, -4.522917997873378e-05, -4.785123418633101e-05, 9.011203663572709e-11, -4.958907285937110e-04, -5.318196157311574e-04, -5.500745218157370e-04, -5.669286129158061e-04, -5.835255272207452e-04, 4.954382075632889e-10, -1.236292175877507e-03, -1.546919838983124e-03, -1.599734815280453e-03, -1.634488112664665e-03, -1.665751515632273e-03,
            }
        };

        // Run simulation
        AnalyzeDC(dc, ckt, exports, references);
    }

    [Test]
    public void When_BSIM4SmallSignal_Expect_Reference()
    {
        var ckt = new Circuit(
            new VoltageSource("Vsupply", "vdd", "0", 0.9),
            new VoltageSource("V1", "in", "0", 0.0),
            new Resistor("R1", "out", "g", 100.0e3),
            new CurrentSource("I1", "vdd", "out", 20e-6),
            new Capacitor("C1", "in", "g", 1e-12),
            CreateMosfet("M1", "out", "g", "0", "0", 1e-6, 100e-9, "mod"),
            CreateModel("mod", "binunit=1 paramchk=1 mobmod=0 capmod=2 igcmod=1 igbmod=1 geomod=1 diomod=1 rdsmod=0 rbodymod=1 rgatemod=1 permod=1 acnqsmod=0 trnqsmod=0 tnom=27 toxe=1.8e-009 toxp=1.5e-009 toxm=1.8e-009 dtox=3e-010 epsrox=3.9 wint=5e-009 lint=0 ll=0 wl=0 lln=1 wln=1 lw=0 ww=0 lwn=1 wwn=1 lwl=0 wwl=0 xpart=0 toxref=1.8e-009 vth0=0.62261 k1=0.4 k2=0 k3=0 k3b=0 w0=2.5e-006 dvt0=1 dvt1=2 dvt2=0 dvt0w=0 dvt1w=0 dvt2w=0 dsub=0.1 minv=0.05 voffl=0 dvtp0=1e-010 dvtp1=0.1 lpe0=0 lpeb=0 xj=1.4e-008 ngate=1e023 ndep=3.24e018 nsd=2e020 phin=0 cdsc=0 cdscb=0 cdscd=0 cit=0 voff=-0.13 nfactor=1.6 eta0=0.0125 etab=0 vfb=-0.55 u0=0.049 ua=6e-010 ub=1.2e-018 uc=0 vsat=130000 a0=1 ags=0 a1=0 a2=1 b0=0 b1=0 keta=0.04 dwg=0 dwb=0 pclm=0.02 pdiblc1=0.001 pdiblc2=0.001 pdiblcb=-0.005 drout=0.5 pvag=1e-020 delta=0.01 pscbe1=8.14e008 pscbe2=1e-007 fprout=0.2 pdits=0.08 pditsd=0.23 pditsl=2300000 rsh=5 rdsw=210 rsw=80 rdw=80 rdswmin=0 rdwmin=0 rswmin=0 prwg=0 prwb=0 wr=1 alpha0=0.074 alpha1=0.005 beta0=30 agidl=0.0002 bgidl=2.1e009 cgidl=0.0002 egidl=0.8 aigbacc=0.012 bigbacc=0.0028 cigbacc=0.002 nigbacc=1 aigbinv=0.014 bigbinv=0.004 cigbinv=0.004 eigbinv=1.1 nigbinv=3 aigc=0.015211 bigc=0.0027432 cigc=0.002 aigsd=0.015211 bigsd=0.0027432 cigsd=0.002 nigc=1 poxedge=1 pigcd=1 ntox=1 xrcrg1=12 xrcrg2=5 cgso=1.1e-010 cgdo=1.1e-010 cgbo=2.56e-011 cgdl=2.653e-010 cgsl=2.653e-010 ckappas=0.03 ckappad=0.03 acde=1 moin=15 noff=0.9 voffcv=0.02 kt1=-0.11 kt1l=0 kt2=0.022 ute=-1.5 ua1=4.31e-009 ub1=7.61e-018 uc1=-5.6e-011 prt=0 at=33000 fnoimod=1 tnoimod=0 jss=0.0001 jsws=1e-011 jswgs=1e-010 njs=1 ijthsfwd=0.01 ijthsrev=0.001 bvs=10 xjbvs=1 jsd=0.0001 jswd=1e-011 jswgd=1e-010 njd=1 ijthdfwd=0.01 ijthdrev=0.001 bvd=10 xjbvd=1 pbs=1 cjs=0.0005 mjs=0.5 pbsws=1 cjsws=5e-010 mjsws=0.33 pbswgs=1 cjswgs=3e-010 mjswgs=0.33 pbd=1 cjd=0.0005 mjd=0.5 pbswd=1 cjswd=5e-010 mjswd=0.33 pbswgd=1 cjswgd=5e-010 mjswgd=0.33 tpb=0.005 tcj=0.001 tpbsw=0.005 tcjsw=0.001 tpbswg=0.005 tcjswg=0.001 xtis=3 xtid=3 dmcg=0 dmci=0 dmdg=0 dmcgt=0 dwj=0 xgw=0 xgl=0 rshg=0.4 gbmin=1e-010 rbpb=5 rbpd=15 rbps=15 rbdb=15 rbsb=15 ngcon=1"));
        ckt["V1"].SetParameter("acmag", 1.0);

        // Create simulation
        var ac = new AC("ac 1", new DecadeSweep(0.1, 1.0e9, 20));

        // Create exports
        var exports = new IExport<Complex>[] { new ComplexVoltageExport(ac, "out") };

        // Create reference
        // 20220724 - Sven Boulanger - using ngSpice
        double[] riref =
        {
            -1.973331280539948e-16,-5.954974936854565e-08, -2.484276894960090e-16,-6.681591774056490e-08, -3.127519312998352e-16,-7.496869274603577e-08, -3.937313539010662e-16,-8.411625675594586e-08, -4.956784068461717e-16,-9.437999238697346e-08, -6.240221424562981e-16,-1.058960931750630e-07, -7.855973326604757e-16,-1.188173707809009e-07, -9.890084455239884e-16,-1.333152827078097e-07, -1.245087864549192e-15,-1.495822074386459e-07, -1.567472752597510e-15,-1.678339971813851e-07, -1.973331280539962e-15,-1.883128420967778e-07, -2.484276894960084e-15,-2.112904840146363e-07, -3.127519312998361e-15,-2.370718222828159e-07, -3.937313539010662e-15,-2.659989595963149e-07, -4.956784068461764e-15,-2.984557414921875e-07, -6.240221424562981e-15,-3.348728497466105e-07, -7.855973326604662e-15,-3.757335172603856e-07, -9.890084455239846e-15,-4.215799402659380e-07, -1.245087864549200e-14,-4.730204729419015e-07, -1.567472752597515e-14,-5.307376999034562e-07, -1.973331280539959e-14,-5.954974936854550e-07, -2.484276894960065e-14,-6.681591774056472e-07, -3.127519312998340e-14,-7.496869274603554e-07, -3.937313539010667e-14,-8.411625675594552e-07, -4.956784068461701e-14,-9.437999238697310e-07, -6.240221424562865e-14,-1.058960931750626e-06, -7.855973326604662e-14,-1.188173707809002e-06, -9.890084455239750e-14,-1.333152827078088e-06, -1.245087864549178e-13,-1.495822074386447e-06, -1.567472752597494e-13,-1.678339971813834e-06, -1.973331280539925e-13,-1.883128420967755e-06, -2.484276894960020e-13,-2.112904840146332e-06, -3.127519312998291e-13,-2.370718222828116e-06, -3.937313539010582e-13,-2.659989595963087e-06, -4.956784068461584e-13,-2.984557414921791e-06, -6.240221424562732e-13,-3.348728497465988e-06, -7.855973326604392e-13,-3.757335172603690e-06, -9.890084455239245e-13,-4.215799402659147e-06, -1.245087864549101e-12,-4.730204729418689e-06, -1.567472752597371e-12,-5.307376999034098e-06, -1.973331280539750e-12,-5.954974936853898e-06, -2.484276894959741e-12,-6.681591774055556e-06, -3.127519312997818e-12,-7.496869274602265e-06, -3.937313539009791e-12,-8.411625675592739e-06, -4.956784068460382e-12,-9.437999238694744e-06, -6.240221424560816e-12,-1.058960931750264e-05, -7.855973326601364e-12,-1.188173707808492e-05, -9.890084455234508e-12,-1.333152827077368e-05, -1.245087864548346e-11,-1.495822074385429e-05, -1.567472752596161e-11,-1.678339971812398e-05, -1.973331280537824e-11,-1.883128420965728e-05, -2.484276894956698e-11,-2.112904840143469e-05, -3.127519312993005e-11,-2.370718222824073e-05, -3.937313539002201e-11,-2.659989595957379e-05, -4.956784068448344e-11,-2.984557414913729e-05, -6.240221424541767e-11,-3.348728497454601e-05, -7.855973326571094e-11,-3.757335172587608e-05, -9.890084455186527e-11,-4.215799402636431e-05, -1.245087864540750e-10,-4.730204729386603e-05, -1.567472752584136e-10,-5.307376998988779e-05, -1.973331280518783e-10,-5.954974936789887e-05, -2.484276894926488e-10,-6.681591773965142e-05, -3.127519312945103e-10,-7.496869274474551e-05, -3.937313538926305e-10,-8.411625675412339e-05, -4.956784068328113e-10,-9.437999238439931e-05, -6.240221424351148e-10,-1.058960931714270e-04, -7.855973326268921e-10,-1.188173707757651e-04, -9.890084454707759e-10,-1.333152827005553e-04, -1.245087864464862e-09,-1.495822074283989e-04, -1.567472752463874e-09,-1.678339971669110e-04, -1.973331280328165e-09,-1.883128420763329e-04, -2.484276894624426e-09,-2.112904839857574e-04, -3.127519312466360e-09,-2.370718222420237e-04, -3.937313538167506e-09,-2.659989595386946e-04, -4.956784067125449e-09,-2.984557414107971e-04, -6.240221422445152e-09,-3.348728496316439e-04, -7.855973323248177e-09,-3.757335170979913e-04, -9.890084449920057e-09,-4.215799400365503e-04, -1.245087863706086e-08,-4.730204726178835e-04, -1.567472751261267e-08,-5.307376994457691e-04, -1.973331278422167e-08,-5.954974930389553e-04, -2.484276891600741e-08,-6.681591764920554e-04, -3.127519307675085e-08,-7.496869261699864e-04, -3.937313530575068e-08,-8.411625657368881e-04, -4.956784055093706e-08,-9.437999212954286e-04, -6.240221403377926e-08,-1.058960928114487e-03, -7.855973293022007e-08,-1.188173702672312e-03, -9.890084402020757e-08,-1.333152819822721e-03, -1.245087856113848e-07,-1.495822064137553e-03, -1.567472739226179e-07,-1.678339957335699e-03, -1.973331259348499e-07,-1.883128400517169e-03, -2.484276861378342e-07,-2.112904811261017e-03, -3.127519259776892e-07,-2.370718182027289e-03, -3.937313454659509e-07,-2.659989538330103e-03, -4.956783934770355e-07,-2.984557333511882e-03, -6.240221212677179e-07,-3.348728382471648e-03, -7.855972990787328e-07,-3.757335010169621e-03, -9.890083923004249e-07,-4.215799173214738e-03, -1.245087780195855e-06,-4.730204405320444e-03, -1.567472618906065e-06,-5.307376541232462e-03, -1.973331068655241e-06,-5.954974290194816e-03, -2.484276559143037e-06,-6.681590860622139e-03, -3.127518780766154e-06,-7.496867984345656e-03, -3.937312695472380e-06,-8.411623853049211e-03, -4.956782731554374e-06,-9.437996664293902e-03, -6.240219305706348e-06,-1.058960568106391e-02, -7.855969968442045e-06,-1.188173194147804e-02, -9.890079132899291e-06,-1.333152101511621e-02, -1.245087021017869e-05,-1.495821049498351e-02, -1.567471415689679e-05,-1.678338524120609e-02, -1.973329161685688e-05,-1.883126376047849e-02, -2.484273536802588e-05,-2.112901951620664e-02, -3.127513990681334e-05,-2.370714142679389e-02, -3.937305103703364e-05,-2.659983832599980e-02, -4.956770699414913e-05,-2.984549273961620e-02, -6.240200236062149e-05,-3.348716998060926e-02, -7.855939745136063e-05,-3.757318929279611e-02, -9.890031232253032e-05,-4.215776458378438e-02, -1.245079429286295e-04,-4.730172319806483e-02, -1.567459383627423e-04,-5.307331219315574e-02, -1.973310092200182e-04,-5.954910271430239e-02, -2.484243313796717e-04,-6.681500431969037e-02, -3.127466090632929e-04,-7.496740250933828e-02, -3.937229187605425e-04,-8.411443425606746e-02, -4.956650381245011e-04,-9.437741805184148e-02, -6.240009546079426e-04,-1.058924568553995e-01, -7.855637524791011e-04,-1.188122343880079e-01, -9.889552251134762e-04,-1.333080274403883e-01, -1.245003517055577e-03,-1.495719592436457e-01, -1.567339073176795e-03,-1.678195214752370e-01, -1.973119417613613e-03,-1.882923950698027e-01, -2.483941124182369e-03,-2.112616026240936e-01, -3.126987170803466e-03,-2.370310276619751e-01, -3.936470187611243e-03,-2.659413382007573e-01, -4.955447520688189e-03,-2.983743536209033e-01, -6.238103286959729e-03,-3.347578943476106e-01, -7.852616599717296e-03,-3.755711527115154e-01, -9.884764990095718e-03,-4.213506196125434e-01, -1.244244903520592e-02,-4.726965940088476e-01, -1.566136983530297e-02,-5.302802889245333e-01, -1.971214696209267e-02,-5.948515260207911e-01, -2.480923266125773e-02,-6.672469771341284e-01, -3.122206026444080e-02,-7.483988604984393e-01, -3.928896249493941e-02,-8.393439245101301e-01, -4.943450943093997e-02,-9.412324433073040e-01, -6.219104550538906e-02,-1.055336792883550e+00, -7.822534632199030e-02,-1.183058956606662e+00, -9.837146032803767e-02,-1.325936001728316e+00, -1.236709300813643e-01,-1.485642144731053e+00, -1.554216718267490e-01,-1.663985447063651e+00, -1.952367780522909e-01,-1.862896412866533e+00, -2.451143115429754e-01,-2.084404753941481e+00, -3.075186514848643e-01,-2.330599310535641e+00, -3.854729403764414e-01,-2.603564563637852e+00, -4.826603942649869e-01,-2.905285404086179e+00, -6.035293204336464e-01,-3.237510018289058e+00, -7.533922079119794e-01,-3.601559310062671e+00, -9.385027894076231e-01,-3.998070948085211e+00, -1.166086267521213e+00,-4.426668150217980e+00, -1.444287224057762e+00,-4.885549540890601e+00, -1.781988455794924e+00,-5.371009257998529e+00, -2.188445352414248e+00,-5.876918449350357e+00, -2.672680358066598e+00,-6.394231771787498e+00, -3.242599226236356e+00,-6.910623478116531e+00, -3.903834986201087e+00,-7.410398594605336e+00, -4.658402849793489e+00,-7.874847817629554e+00, -5.503353667906635e+00,-8.283194196188862e+00, -6.429717946945069e+00,-8.614190238451057e+00, -7.422087151257924e+00,-8.848259607127540e+00, -8.459127348457898e+00,-8.969871591850088e+00, -9.515133562842177e+00,-8.969667688127512e+00, -1.056244785177438e+01,-8.845823348408784e+00, -1.157428648689068e+01,-8.604278615024704e+00, -1.252737979512019e+01,-8.257772182658043e+00, -1.340389371282710e+01,-7.823942349860427e+00, -1.419234144983875e+01,-7.322976376962988e+00, -1.488749137045521e+01,-6.775322948847078e+00, -1.548951017352869e+01,-6.199853907561818e+00, -1.600268330089853e+01,-5.612658437940168e+00, -1.643403438134978e+01,-5.026465947312281e+00, -1.679207324869834e+01,-4.450574702833648e+00, -1.708579306590721e+01,-3.891118824743717e+00, -1.732394761754509e+01,-3.351516017690559e+00, -1.751458310042834e+01,-2.832975254248971e+00, -1.766477218513856e+01,-2.334986195045595e+00, -1.778049259597332e+01,-1.855748585800912e+00, -1.786659851452369e+01,-1.392526109838460e+00, -1.792684374276608e+01,-9.419256383387000e-01, -1.796392660147298e+01,-5.001118203552337e-01, -1.797953604782318e+01,-6.297106883928294e-02, -1.797438593864741e+01,3.737596395210603e-01, -1.794822996449392e+01,8.143836167062132e-01, -1.789985410078243e+01,1.263139627179675e+00, -1.782704716133626e+01,1.724086042025338e+00, -1.772655393037542e+01,2.200958395490180e+00, -1.759402012474940e+01,2.696984435679622e+00, -1.742394478186552e+01,3.214641018113864e+00, -1.720966408705294e+01,3.755339163495082e+00, -1.694340119918808e+01,4.319029053093189e+00, -1.661642840005847e+01,4.903728069997145e+00, -1.621939826732769e+01,5.504994943745728e+00,
        };
        var references = new Complex[1][];
        references[0] = new Complex[riref.Length / 2];
        for (var i = 0; i < riref.Length; i += 2)
            references[0][i / 2] = new Complex(riref[i], riref[i + 1]);

        // run simulation
        AnalyzeAC(ac, ckt, exports, references);
    }

    [Test]
    public void When_BSIM4Noise_Expect_NoException()
    {
        var ckt = new Circuit(
            new VoltageSource("Vsupply", "vdd", "0", 0.9),
            new VoltageSource("V1", "in", "0", 0.0),
            new Resistor("R1", "out", "g", 100.0e3),
            new CurrentSource("I1", "vdd", "out", 20e-6),
            new Capacitor("C1", "in", "g", 1e-12),
            CreateMosfet("M1", "out", "g", "0", "0", 1e-6, 100e-9, "mod"),
            CreateModel("mod", "binunit=1 paramchk=1 mobmod=0 capmod=2 igcmod=1 igbmod=1 geomod=1 diomod=1 rdsmod=0 rbodymod=1 rgatemod=1 permod=1 acnqsmod=0 trnqsmod=0 tnom=27 toxe=1.8e-009 toxp=1.5e-009 toxm=1.8e-009 dtox=3e-010 epsrox=3.9 wint=5e-009 lint=0 ll=0 wl=0 lln=1 wln=1 lw=0 ww=0 lwn=1 wwn=1 lwl=0 wwl=0 xpart=0 toxref=1.8e-009 vth0=0.62261 k1=0.4 k2=0 k3=0 k3b=0 w0=2.5e-006 dvt0=1 dvt1=2 dvt2=0 dvt0w=0 dvt1w=0 dvt2w=0 dsub=0.1 minv=0.05 voffl=0 dvtp0=1e-010 dvtp1=0.1 lpe0=0 lpeb=0 xj=1.4e-008 ngate=1e023 ndep=3.24e018 nsd=2e020 phin=0 cdsc=0 cdscb=0 cdscd=0 cit=0 voff=-0.13 nfactor=1.6 eta0=0.0125 etab=0 vfb=-0.55 u0=0.049 ua=6e-010 ub=1.2e-018 uc=0 vsat=130000 a0=1 ags=0 a1=0 a2=1 b0=0 b1=0 keta=0.04 dwg=0 dwb=0 pclm=0.02 pdiblc1=0.001 pdiblc2=0.001 pdiblcb=-0.005 drout=0.5 pvag=1e-020 delta=0.01 pscbe1=8.14e008 pscbe2=1e-007 fprout=0.2 pdits=0.08 pditsd=0.23 pditsl=2300000 rsh=5 rdsw=210 rsw=80 rdw=80 rdswmin=0 rdwmin=0 rswmin=0 prwg=0 prwb=0 wr=1 alpha0=0.074 alpha1=0.005 beta0=30 agidl=0.0002 bgidl=2.1e009 cgidl=0.0002 egidl=0.8 aigbacc=0.012 bigbacc=0.0028 cigbacc=0.002 nigbacc=1 aigbinv=0.014 bigbinv=0.004 cigbinv=0.004 eigbinv=1.1 nigbinv=3 aigc=0.015211 bigc=0.0027432 cigc=0.002 aigsd=0.015211 bigsd=0.0027432 cigsd=0.002 nigc=1 poxedge=1 pigcd=1 ntox=1 xrcrg1=12 xrcrg2=5 cgso=1.1e-010 cgdo=1.1e-010 cgbo=2.56e-011 cgdl=2.653e-010 cgsl=2.653e-010 ckappas=0.03 ckappad=0.03 acde=1 moin=15 noff=0.9 voffcv=0.02 kt1=-0.11 kt1l=0 kt2=0.022 ute=-1.5 ua1=4.31e-009 ub1=7.61e-018 uc1=-5.6e-011 prt=0 at=33000 fnoimod=1 tnoimod=0 jss=0.0001 jsws=1e-011 jswgs=1e-010 njs=1 ijthsfwd=0.01 ijthsrev=0.001 bvs=10 xjbvs=1 jsd=0.0001 jswd=1e-011 jswgd=1e-010 njd=1 ijthdfwd=0.01 ijthdrev=0.001 bvd=10 xjbvd=1 pbs=1 cjs=0.0005 mjs=0.5 pbsws=1 cjsws=5e-010 mjsws=0.33 pbswgs=1 cjswgs=3e-010 mjswgs=0.33 pbd=1 cjd=0.0005 mjd=0.5 pbswd=1 cjswd=5e-010 mjswd=0.33 pbswgd=1 cjswgd=5e-010 mjswgd=0.33 tpb=0.005 tcj=0.001 tpbsw=0.005 tcjsw=0.001 tpbswg=0.005 tcjswg=0.001 xtis=3 xtid=3 dmcg=0 dmci=0 dmdg=0 dmcgt=0 dwj=0 xgw=0 xgl=0 rshg=0.4 gbmin=1e-010 rbpb=5 rbpd=15 rbps=15 rbdb=15 rbsb=15 ngcon=1"));
        ckt["V1"].SetParameter("acmag", 1.0);

        // Create simulation
        var noise = new Noise("ac 1", "V1", "out", new DecadeSweep(0.1, 1.0e9, 20));
        noise.Run(ckt);
    }

    [Test]
    public void When_BSIM4Transient_Expect_NoErrors()
    {
        var ckt = new Circuit(
            new VoltageSource("V1", "g", "0", new Pulse(0, 5, 0.1e-7, 1e-9, 1e-9, 0.6e-6, 1e-6)),
            new VoltageSource("V2", "vdd", "0", 3.3),
            new Resistor("R1", "vdd", "d", 1e6),
            CreateMosfet("M1", "d", "g", "0", "0", 2e-6, 100e-9, "mod"),
            CreateModel("mod", "binunit=1 paramchk=1 mobmod=0 capmod=2 igcmod=1 igbmod=1 geomod=1 diomod=1 rdsmod=0 rbodymod=1 rgatemod=1 permod=1 acnqsmod=0 trnqsmod=0 tnom=27 toxe=1.8e-009 toxp=1.5e-009 toxm=1.8e-009 dtox=3e-010 epsrox=3.9 wint=5e-009 lint=0 ll=0 wl=0 lln=1 wln=1 lw=0 ww=0 lwn=1 wwn=1 lwl=0 wwl=0 xpart=0 toxref=1.8e-009 vth0=0.62261 k1=0.4 k2=0 k3=0 k3b=0 w0=2.5e-006 dvt0=1 dvt1=2 dvt2=0 dvt0w=0 dvt1w=0 dvt2w=0 dsub=0.1 minv=0.05 voffl=0 dvtp0=1e-010 dvtp1=0.1 lpe0=0 lpeb=0 xj=1.4e-008 ngate=1e023 ndep=3.24e018 nsd=2e020 phin=0 cdsc=0 cdscb=0 cdscd=0 cit=0 voff=-0.13 nfactor=1.6 eta0=0.0125 etab=0 vfb=-0.55 u0=0.049 ua=6e-010 ub=1.2e-018 uc=0 vsat=130000 a0=1 ags=0 a1=0 a2=1 b0=0 b1=0 keta=0.04 dwg=0 dwb=0 pclm=0.02 pdiblc1=0.001 pdiblc2=0.001 pdiblcb=-0.005 drout=0.5 pvag=1e-020 delta=0.01 pscbe1=8.14e008 pscbe2=1e-007 fprout=0.2 pdits=0.08 pditsd=0.23 pditsl=2300000 rsh=5 rdsw=210 rsw=80 rdw=80 rdswmin=0 rdwmin=0 rswmin=0 prwg=0 prwb=0 wr=1 alpha0=0.074 alpha1=0.005 beta0=30 agidl=0.0002 bgidl=2.1e009 cgidl=0.0002 egidl=0.8 aigbacc=0.012 bigbacc=0.0028 cigbacc=0.002 nigbacc=1 aigbinv=0.014 bigbinv=0.004 cigbinv=0.004 eigbinv=1.1 nigbinv=3 aigc=0.015211 bigc=0.0027432 cigc=0.002 aigsd=0.015211 bigsd=0.0027432 cigsd=0.002 nigc=1 poxedge=1 pigcd=1 ntox=1 xrcrg1=12 xrcrg2=5 cgso=1.1e-010 cgdo=1.1e-010 cgbo=2.56e-011 cgdl=2.653e-010 cgsl=2.653e-010 ckappas=0.03 ckappad=0.03 acde=1 moin=15 noff=0.9 voffcv=0.02 kt1=-0.11 kt1l=0 kt2=0.022 ute=-1.5 ua1=4.31e-009 ub1=7.61e-018 uc1=-5.6e-011 prt=0 at=33000 fnoimod=1 tnoimod=0 jss=0.0001 jsws=1e-011 jswgs=1e-010 njs=1 ijthsfwd=0.01 ijthsrev=0.001 bvs=10 xjbvs=1 jsd=0.0001 jswd=1e-011 jswgd=1e-010 njd=1 ijthdfwd=0.01 ijthdrev=0.001 bvd=10 xjbvd=1 pbs=1 cjs=0.0005 mjs=0.5 pbsws=1 cjsws=5e-010 mjsws=0.33 pbswgs=1 cjswgs=3e-010 mjswgs=0.33 pbd=1 cjd=0.0005 mjd=0.5 pbswd=1 cjswd=5e-010 mjswd=0.33 pbswgd=1 cjswgd=5e-010 mjswgd=0.33 tpb=0.005 tcj=0.001 tpbsw=0.005 tcjsw=0.001 tpbswg=0.005 tcjswg=0.001 xtis=3 xtid=3 dmcg=0 dmci=0 dmdg=0 dmcgt=0 dwj=0 xgw=0 xgl=0 rshg=0.4 gbmin=1e-010 rbpb=5 rbpd=15 rbps=15 rbdb=15 rbsb=15 ngcon=1"));
        ckt["M1"].SetParameter("m", 4.0);

        var tran = new Transient("Transient 1", 1e-8, 10e-6);
        var export = new RealVoltageExport(tran, "d");

        foreach (int _ in tran.Run(ckt, Transient.ExportTransient))
        {
            Assert.That(export.Value, Is.LessThan(4.0));
            Assert.That(export.Value, Is.GreaterThan(-1.0));
        }
    }

    [Test]
    public void When_BSIM3v2AggregateModel_Expect_NoErrors()
    {
        var ckt = new Circuit(
            new VoltageSource("V1", "g", "0", 0.0),
            new VoltageSource("V2", "d", "0", 0.0),
            CreateMosfet("M1", "d", "g", "0", "0", 10e-6, 1e-6, "mod"),
            CreateMosfet("M2", "d", "g", "0", "0", 10e-6, 1.5e-6, "mod"),
            new BSIM4AggregateModel(
                "mod",
                [
                    CreateModel("mod.1", "lmin=0.5e-6 lmax=1e-6 wmin=0.5e-6 wmax=100e-6 binunit=1 paramchk=1 mobmod=0 capmod=2 igcmod=1 igbmod=1 geomod=1 diomod=1 rdsmod=0 rbodymod=1 rgatemod=1 permod=1 acnqsmod=0 trnqsmod=0 tnom=27 toxe=1.8e-009 toxp=1.5e-009 toxm=1.8e-009 dtox=3e-010 epsrox=3.9 wint=5e-009 lint=0 ll=0 wl=0 lln=1 wln=1 lw=0 ww=0 lwn=1 wwn=1 lwl=0 wwl=0 xpart=0 toxref=1.8e-009 vth0=0.62261 k1=0.4 k2=0 k3=0 k3b=0 w0=2.5e-006 dvt0=1 dvt1=2 dvt2=0 dvt0w=0 dvt1w=0 dvt2w=0 dsub=0.1 minv=0.05 voffl=0 dvtp0=1e-010 dvtp1=0.1 lpe0=0 lpeb=0 xj=1.4e-008 ngate=1e023 ndep=3.24e018 nsd=2e020 phin=0 cdsc=0 cdscb=0 cdscd=0 cit=0 voff=-0.13 nfactor=1.6 eta0=0.0125 etab=0 vfb=-0.55 u0=0.049 ua=6e-010 ub=1.2e-018 uc=0 vsat=130000 a0=1 ags=0 a1=0 a2=1 b0=0 b1=0 keta=0.04 dwg=0 dwb=0 pclm=0.02 pdiblc1=0.001 pdiblc2=0.001 pdiblcb=-0.005 drout=0.5 pvag=1e-020 delta=0.01 pscbe1=8.14e008 pscbe2=1e-007 fprout=0.2 pdits=0.08 pditsd=0.23 pditsl=2300000 rsh=5 rdsw=210 rsw=80 rdw=80 rdswmin=0 rdwmin=0 rswmin=0 prwg=0 prwb=0 wr=1 alpha0=0.074 alpha1=0.005 beta0=30 agidl=0.0002 bgidl=2.1e009 cgidl=0.0002 egidl=0.8 aigbacc=0.012 bigbacc=0.0028 cigbacc=0.002 nigbacc=1 aigbinv=0.014 bigbinv=0.004 cigbinv=0.004 eigbinv=1.1 nigbinv=3 aigc=0.015211 bigc=0.0027432 cigc=0.002 aigsd=0.015211 bigsd=0.0027432 cigsd=0.002 nigc=1 poxedge=1 pigcd=1 ntox=1 xrcrg1=12 xrcrg2=5 cgso=1.1e-010 cgdo=1.1e-010 cgbo=2.56e-011 cgdl=2.653e-010 cgsl=2.653e-010 ckappas=0.03 ckappad=0.03 acde=1 moin=15 noff=0.9 voffcv=0.02 kt1=-0.11 kt1l=0 kt2=0.022 ute=-1.5 ua1=4.31e-009 ub1=7.61e-018 uc1=-5.6e-011 prt=0 at=33000 fnoimod=1 tnoimod=0 jss=0.0001 jsws=1e-011 jswgs=1e-010 njs=1 ijthsfwd=0.01 ijthsrev=0.001 bvs=10 xjbvs=1 jsd=0.0001 jswd=1e-011 jswgd=1e-010 njd=1 ijthdfwd=0.01 ijthdrev=0.001 bvd=10 xjbvd=1 pbs=1 cjs=0.0005 mjs=0.5 pbsws=1 cjsws=5e-010 mjsws=0.33 pbswgs=1 cjswgs=3e-010 mjswgs=0.33 pbd=1 cjd=0.0005 mjd=0.5 pbswd=1 cjswd=5e-010 mjswd=0.33 pbswgd=1 cjswgd=5e-010 mjswgd=0.33 tpb=0.005 tcj=0.001 tpbsw=0.005 tcjsw=0.001 tpbswg=0.005 tcjswg=0.001 xtis=3 xtid=3 dmcg=0 dmci=0 dmdg=0 dmcgt=0 dwj=0 xgw=0 xgl=0 rshg=0.4 gbmin=1e-010 rbpb=5 rbpd=15 rbps=15 rbdb=15 rbsb=15 ngcon=1"),
                    CreateModel("mod.2", "lmin=1e-6 lmax=2e-6 wmin=0.5e-6 wmax=100e-6 binunit=1 paramchk=1 mobmod=0 capmod=2 igcmod=1 igbmod=1 geomod=1 diomod=1 rdsmod=0 rbodymod=1 rgatemod=1 permod=1 acnqsmod=0 trnqsmod=0 tnom=27 toxe=1.8e-009 toxp=1.5e-009 toxm=1.8e-009 dtox=3e-010 epsrox=3.9 wint=5e-009 lint=0 ll=0 wl=0 lln=1 wln=1 lw=0 ww=0 lwn=1 wwn=1 lwl=0 wwl=0 xpart=0 toxref=1.8e-009 vth0=0.62261 k1=0.4 k2=0 k3=0 k3b=0 w0=2.5e-006 dvt0=1 dvt1=2 dvt2=0 dvt0w=0 dvt1w=0 dvt2w=0 dsub=0.1 minv=0.05 voffl=0 dvtp0=1e-010 dvtp1=0.1 lpe0=0 lpeb=0 xj=1.4e-008 ngate=1e023 ndep=3.24e018 nsd=2e020 phin=0 cdsc=0 cdscb=0 cdscd=0 cit=0 voff=-0.13 nfactor=1.6 eta0=0.0125 etab=0 vfb=-0.55 u0=0.049 ua=6e-010 ub=1.2e-018 uc=0 vsat=130000 a0=1 ags=0 a1=0 a2=1 b0=0 b1=0 keta=0.04 dwg=0 dwb=0 pclm=0.02 pdiblc1=0.001 pdiblc2=0.001 pdiblcb=-0.005 drout=0.5 pvag=1e-020 delta=0.01 pscbe1=8.14e008 pscbe2=1e-007 fprout=0.2 pdits=0.08 pditsd=0.23 pditsl=2300000 rsh=5 rdsw=210 rsw=80 rdw=80 rdswmin=0 rdwmin=0 rswmin=0 prwg=0 prwb=0 wr=1 alpha0=0.074 alpha1=0.005 beta0=30 agidl=0.0002 bgidl=2.1e009 cgidl=0.0002 egidl=0.8 aigbacc=0.012 bigbacc=0.0028 cigbacc=0.002 nigbacc=1 aigbinv=0.014 bigbinv=0.004 cigbinv=0.004 eigbinv=1.1 nigbinv=3 aigc=0.015211 bigc=0.0027432 cigc=0.002 aigsd=0.015211 bigsd=0.0027432 cigsd=0.002 nigc=1 poxedge=1 pigcd=1 ntox=1 xrcrg1=12 xrcrg2=5 cgso=1.1e-010 cgdo=1.1e-010 cgbo=2.56e-011 cgdl=2.653e-010 cgsl=2.653e-010 ckappas=0.03 ckappad=0.03 acde=1 moin=15 noff=0.9 voffcv=0.02 kt1=-0.11 kt1l=0 kt2=0.022 ute=-1.5 ua1=4.31e-009 ub1=7.61e-018 uc1=-5.6e-011 prt=0 at=33000 fnoimod=1 tnoimod=0 jss=0.0001 jsws=1e-011 jswgs=1e-010 njs=1 ijthsfwd=0.01 ijthsrev=0.001 bvs=10 xjbvs=1 jsd=0.0001 jswd=1e-011 jswgd=1e-010 njd=1 ijthdfwd=0.01 ijthdrev=0.001 bvd=10 xjbvd=1 pbs=1 cjs=0.0005 mjs=0.5 pbsws=1 cjsws=5e-010 mjsws=0.33 pbswgs=1 cjswgs=3e-010 mjswgs=0.33 pbd=1 cjd=0.0005 mjd=0.5 pbswd=1 cjswd=5e-010 mjswd=0.33 pbswgd=1 cjswgd=5e-010 mjswgd=0.33 tpb=0.005 tcj=0.001 tpbsw=0.005 tcjsw=0.001 tpbswg=0.005 tcjswg=0.001 xtis=3 xtid=3 dmcg=0 dmci=0 dmdg=0 dmcgt=0 dwj=0 xgw=0 xgl=0 rshg=0.4 gbmin=1e-010 rbpb=5 rbpd=15 rbps=15 rbdb=15 rbsb=15 ngcon=1")
                ]));
        ckt["M1"].SetParameter("m", 4.0);

        // Create simulation
        var dc = new DC("dc",
        [
            new ParameterSweep("V1", new LinearSweep(0, 3.3, 0.3)),
            new ParameterSweep("V2", new LinearSweep(0, 3.3, 0.3))
        ]);

        // Make sure the right models are being used
        bool hit = false;
        foreach (var state in dc.Run(ckt, -1))
        {
            if (state == BiasingSimulation.AfterTemperature)
            {
                // Make sure that the right models are being used
                var m1Behavior = dc.EntityBehaviors["M1"].GetValue<SpiceSharpBSIM.Components.Semiconductors.BSIM.BSIM4Behaviors.BiasingBehavior>();
                Assert.That(m1Behavior.ModelName, Is.EqualTo("mod.1"));
                var m2Behavior = dc.EntityBehaviors["M2"].GetValue<SpiceSharpBSIM.Components.Semiconductors.BSIM.BSIM4Behaviors.BiasingBehavior>();
                Assert.That(m2Behavior.ModelName, Is.EqualTo("mod.2"));
                hit = true;
            }
        }
        Assert.That(hit, Is.True);
    }
}
