﻿using System.Numerics;
using NUnit.Framework;
using SpiceSharp;
using SpiceSharp.Components;
using SpiceSharp.Simulations;

namespace SpiceSharpTest.Models
{
    [TestFixture]
    public class BSIM2Tests : Framework
    {
        private BSIM2 CreateMosfet(string name, string drain, string gate, string source, string bulk, double w, double l, string model)
        {
            var e = new BSIM2(name, drain, gate, source, bulk);
            e.SetParameter("w", w);
            e.SetParameter("l", l);
            e.Model = model;
            return e;
        }

        private BSIM2Model CreateModel(string name, string parameters)
        {
            var m = new BSIM2Model(name);
            ApplyParameters(m, parameters);
            return m;
        }

        [Test]
        public void When_BSIM2DC_Expect_Reference()
        {
            var ckt = new Circuit(
                new VoltageSource("V1", "g", "0", 0.0),
                new VoltageSource("V2", "d", "0", 0.0),
                CreateMosfet("M1", "d", "g", "0", "0", 10e-6, 1e-6, "mod"),
                CreateModel("mod", "vfb=-0.3 phi=0.8 k1=0.6 mu0=250 n0=1.3 tox=1e-7 mj=0.5 mjsw=0.33 pb=0.8 pbsw=1.0 xpart=1.0")
            );
            ckt["M1"].SetParameter("m", 3.0);

            // Create simulation
            var dc = new DC("dc", new[]
            {
                new ParameterSweep("V1", new LinearSweep(0, 3.3, 0.3)),
                new ParameterSweep("V2", new LinearSweep(0, 3.3, 0.3))
            });

            // Create exports
            var exports = new IExport<double>[] { new RealPropertyExport(dc, "V2", "i") };

            // Create references
            // 20220723 - Sven Boulanger - using ngSpice
            double[][] references =
            {
                new[]
                {
                    -3.000000000000000e-50, -3.013116867112740e-02, -3.033099246666570e-02, -3.053053847413867e-02, -3.073008447903643e-02, -3.092963048393419e-02, -3.112917648883196e-02, -3.132872249372974e-02, -3.152826849862750e-02, -3.172781450352526e-02, -3.192736050842302e-02, -3.212690651332079e-02, -3.000000000000000e-50, -2.628379371797995e-01, -2.645810250224200e-01, -2.663216896851476e-01, -2.680623543254114e-01, -2.698030189656752e-01, -2.715436836059391e-01, -2.732843482462029e-01, -2.750250128864669e-01, -2.767656775267306e-01, -2.785063421669944e-01, -2.802470068072583e-01, -3.000000000000000e-50, -1.239885715521496e+00, -1.248102176622370e+00, -1.256313375203478e+00, -1.264524573736018e+00, -1.272735772268559e+00, -1.280946970801100e+00, -1.289158169333640e+00, -1.297369367866180e+00, -1.305580566398721e+00, -1.313791764931261e+00, -1.322002963463802e+00, -3.000000000000000e-50, -1.108156644201510e+01, -1.118547369945625e+01, -1.125915752477491e+01, -1.133274708557860e+01, -1.140633635328996e+01, -1.147992562009001e+01, -1.155351488688724e+01, -1.162710415368447e+01, -1.170069342048166e+01, -1.177428268727888e+01, -1.184787195407610e+01, -3.000000000000000e-50, -3.677115458160296e+01, -3.851437512904238e+01, -3.884417435489301e+01, -3.910178024421657e+01, -3.935586883115890e+01, -3.960978645494171e+01, -3.986369576978210e+01, -4.011760468080372e+01, -4.037151357219948e+01, -4.062542246264145e+01, -4.087933135303705e+01, -3.000000000000000e-50, -6.406055062103707e+01, -7.822718864321652e+01, -7.972696704474419e+01, -8.037627931833049e+01, -8.091466121887007e+01, -8.143879998964275e+01, -8.196111369619905e+01, -8.248319360639155e+01, -8.300524356754701e+01, -8.352728969228272e+01, -8.404933532558002e+01, -3.000000000000000e-50, -8.547784869906626e+01, -1.263072740185677e+02, -1.309832090819341e+02, -1.326818776478385e+02, -1.337221945242090e+02, -1.346225750713827e+02, -1.354934618216065e+02, -1.363581433714179e+02, -1.372215198992645e+02, -1.380846219883043e+02, -1.389476663653595e+02, -3.000000000000000e-50, -1.030341528520895e+02, -1.752136820184235e+02, -1.897769951431887e+02, -1.937351643651651e+02, -1.957819147955629e+02, -1.972705484603432e+02, -1.985993138169629e+02, -1.998825371855746e+02, -2.011528072281528e+02, -2.024193946490373e+02, -2.036849352357177e+02, -3.000000000000000e-50, -1.179781804054846e+02, -2.165866179831287e+02, -2.537277909154678e+02, -2.614081256257277e+02, -2.652790583627922e+02, -2.677568830157111e+02, -2.697412462251713e+02, -2.715527746953744e+02, -2.733040058393493e+02, -2.750342297687113e+02, -2.767571384161292e+02, -3.000000000000000e-50, -1.310384863240151e+02, -2.517318776413583e+02, -3.186053470978439e+02, -3.338577711467806e+02, -3.405554871309926e+02, -3.446170207291166e+02, -3.475836657504483e+02, -3.501037847839783e+02, -3.524431409832210e+02, -3.547095392334439e+02, -3.569465263202291e+02, -3.000000000000000e-50, -1.426602305723909e+02, -2.819553869386261e+02, -3.758486553269672e+02, -4.096113560838789e+02, -4.201556798989479e+02, -4.265254367488924e+02, -4.309286323796449e+02, -4.344287286608684e+02, -4.375188748653064e+02, -4.404239016937602e+02, -4.432455320506682e+02, -3.000000000000000e-50, -1.531338419512231e+02, -3.082847772420649e+02, -4.259028812269407e+02, -4.864848222262390e+02, -5.028645256583067e+02, -5.123095968828391e+02, -5.186962438287508e+02, -5.235360562221041e+02, -5.276060135979194e+02, -5.312958445564543e+02, -5.347987139058062e+02,
                }
            };

            // Run simulation
            AnalyzeDC(dc, ckt, exports, references);
        }

        [Test]
        public void When_BSIM2SmallSignal_Expect_Reference()
        {
            // Build the circuit
            var ckt = new Circuit(
                new VoltageSource("Vsupply", "vdd", "0", 5.0),
                new VoltageSource("V1", "in", "0", 0.0),
                new Resistor("R1", "vdd", "out", 10),
                new Resistor("R2", "out", "g", 1.0e6),
                new Capacitor("C1", "in", "g", 1e-9),
                CreateMosfet("M1", "out", "g", "0", "0", 10e-6, 1e-6, "Nch4"),
                CreateModel("Nch4", "vfb=-0.3 phi=0.8 k1=0.6 mu0=250 n0=1.3 tox=1e-7 mj=0.5 mjsw=0.33 pb=0.8 pbsw=1.0 xpart=1.0")
            );
            ckt["V1"].SetParameter("acmag", 1.0);

            // AC simulation
            var ac = new AC("ac 1", new DecadeSweep(0.1, 1.0e9, 20));

            // Create exports
            var exports = new IExport<Complex>[] { new ComplexVoltageExport(ac, "out") };

            // Reference
            // 20220723 - Sven Boulanger - using ngSpice
            var riref = new[]
            {
                -6.554049667309590e-08,-5.949569069278144e-04, -8.251059650411451e-08,-6.675526269898970e-04, -1.038746862704300e-07,-7.490063637375844e-04, -1.307704815367668e-07,-8.403989583180630e-04, -1.646302812870750e-07,-9.429431342973807e-04, -2.072572430273195e-07,-1.057999589688593e-03, -2.609214074329220e-07,-1.187095052479249e-03, -3.284805861894588e-07,-1.331942539231032e-03, -4.135325507227128e-07,-1.494464085552403e-03, -5.206066263876775e-07,-1.676816250037754e-03, -6.554048951498711e-07,-1.881418730058017e-03, -8.251058515927704e-07,-2.110986469064865e-03, -1.038746682900752e-06,-2.368565681379483e-03, -1.307704530398266e-06,-2.657574272389166e-03, -1.646302361224716e-06,-2.981847190356611e-03, -2.072571714462533e-06,-3.345687311425817e-03, -2.609212939845900e-06,-3.753922532745187e-03, -3.284804063859944e-06,-4.211969830885801e-03, -4.135322657534759e-06,-4.725907134987520e-03, -5.206061747419727e-06,-5.302553967521622e-03, -6.554041793398648e-06,-5.949561921559532e-03, -8.251047171107563e-06,-6.675516173481279e-03, -1.038744884868714e-05,-7.490049375812489e-03, -1.307701680711090e-05,-8.403969438197084e-03, -1.646297844778025e-05,-9.429402887446231e-03, -2.072564556383134e-05,-1.057995570241692e-02, -2.609201595066990e-05,-1.187089374865292e-02, -3.284786083621821e-05,-1.331934519398285e-02, -4.135294160827125e-05,-1.494452757255566e-02, -5.206016583280283e-05,-1.676800248425198e-02, -6.553970213258022e-05,-1.881396127236483e-02, -8.250933724622074e-05,-2.110954541832155e-02, -1.038726904890693e-04,-2.368520583144761e-02, -1.307673184522420e-04,-2.657510569759824e-02, -1.646252681674041e-04,-2.981757208570729e-02, -2.072492978308493e-04,-3.345560209788090e-02, -2.609088152703628e-04,-3.753742998711576e-02, -3.284606292066192e-04,-4.211716235526571e-02, -4.135009215350278e-04,-4.725548927718522e-02, -5.205564984981538e-04,-5.302047996436003e-02, -6.553254497836525e-04,-5.948847236418829e-02, -8.249799431323313e-04,-6.674506685919969e-02, -1.038547139338917e-03,-7.488623493690040e-02, -1.307388290926650e-03,-8.401955427444396e-02, -1.645801186883031e-03,-9.426558201725745e-02, -2.071777469390437e-03,-1.057593779722655e-01, -2.607954271379040e-03,-1.186521887601098e-01, -3.282809458408465e-03,-1.331133023544509e-01, -4.132161918876167e-03,-1.493320794204404e-01, -5.201053307410592e-03,-1.675201627975785e-01, -6.546105931679620e-03,-1.879138584377164e-01, -8.238473627796358e-03,-2.107766688247573e-01, -1.036752900099877e-02,-2.364019415881902e-01, -1.304546167581720e-02,-2.651155692259869e-01, -1.641299822489913e-02,-2.972786372309494e-01, -2.064649459229738e-02,-3.332898629696233e-01, -2.596669428557188e-02,-3.735875904411628e-01, -3.264948658305430e-02,-4.186509988464130e-01, -4.103903145202754e-02,-4.690000360905471e-01, -5.156362966952391e-02,-5.251933904541369e-01, -6.475468789396806e-02,-5.878235521255866e-01, -8.126902788969644e-02,-6.575076794086726e-01, -1.019145695497472e-01,-7.348725656546508e-01, -1.276790106184019e-01,-8.205315089614079e-01, -1.597604389703666e-01,-9.150503512790030e-01, -1.995977284842867e-01,-1.018899460320625e+00, -2.488969570842776e-01,-1.132388142815550e+00, -3.096478382714614e-01,-1.255578192214069e+00, -3.841212666229644e-01,-1.388174445705050e+00, -4.748357685980908e-01,-1.529392993360592e+00, -5.844775871837812e-01,-1.677812798564466e+00, -7.157576599174186e-01,-1.831224578225676e+00, -8.711910641417192e-01,-1.986501833623950e+00, -1.052793439256749e+00,-2.139531425706404e+00, -1.261706996127091e+00,-2.285251159890267e+00, -1.497796232345845e+00,-2.417842841541054e+00, -1.759285592270818e+00,-2.531113177438989e+00, -2.042536351566931e+00,-2.619056499869221e+00, -2.342061294085567e+00,-2.676536350136931e+00, -2.650839494221056e+00,-2.699964725507677e+00, -2.960920466001533e+00,-2.687825738841905e+00, -3.264219509080916e+00,-2.640908510242945e+00, -3.553341587351890e+00,-2.562185076221188e+00, -3.822261057345647e+00,-2.456367276571195e+00, -4.066734152418690e+00,-2.329260215740997e+00, -4.284405514425669e+00,-2.187065927353993e+00, -4.474651146404934e+00,-2.035772897182377e+00, -4.638249227555456e+00,-1.880714042130676e+00, -4.776979882114945e+00,-1.726315862522970e+00, -4.893235643563924e+00,-1.576016159835154e+00, -4.989692896542604e+00,-1.432304776475938e+00, -5.069064764055015e+00,-1.296838207730490e+00, -5.133934647577032e+00,-1.170586949195543e+00, -5.186657984457888e+00,-1.053986727958434e+00, -5.229315682246922e+00,-9.470765337446351e-01, -5.263703229303631e+00,-8.496154511879788e-01, -5.291342220551417e+00,-7.611763088646935e-01, -5.313504367732735e+00,-6.812175953866666e-01, -5.331241128564367e+00,-6.091366965501382e-01, -5.345414548649014e+00,-5.443079689575597e-01, -5.356726716068980e+00,-4.861090011875207e-01, -5.365746467744401e+00,-4.339379637207876e-01, -5.372932786854940e+00,-3.872244146138868e-01, -5.378654815445882e+00,-3.454354151841080e-01, -5.383208677083229e+00,-3.080783655209749e-01, -5.386831436912339e+00,-2.747016069150974e-01, -5.389712574070368e+00,-2.448935538457110e-01, -5.392003339920421e+00,-2.182809013414752e-01, -5.393824348089624e+00,-1.945262917373779e-01, -5.395271703076938e+00,-1.733257061464961e-01, -5.396421931674954e+00,-1.544057600959194e-01, -5.397335940293758e+00,-1.375210214367302e-01, -5.398062183822407e+00,-1.224514253476046e-01, -5.398639198856820e+00,-1.089998310806112e-01, -5.399097626107978e+00,-9.698974436529532e-02, -5.399461823303692e+00,-8.626321538285681e-02, -5.399751150436109e+00,-7.667891297213034e-02, -5.399980993244221e+00,-6.811036980794687e-02, -5.400163577820808e+00,-6.044438967896443e-02, -5.400308618704545e+00,-5.357960595911005e-02, -5.400423834326141e+00,-4.742517939997318e-02, -5.400515356851048e+00,-4.189962310921026e-02, -5.400588057987425e+00,-3.692974276972596e-02, -5.400645807947679e+00,-3.244968062111211e-02, -5.400691681251838e+00,-2.840005234841261e-02, -5.400728120267916e+00,-2.472716672220105e-02, -5.400757065157653e+00,-2.138231855369220e-02, -5.400780057121924e+00,-1.832114623342011e-02, -5.400798320427825e+00,-1.550304578913677e-02, -5.400812827575387e+00,-1.289063401408483e-02, -5.400824351067841e+00,-1.044925377269410e-02, -5.400833504538303e+00,-8.146515082979019e-03, -5.400840775420448e+00,-5.951866001891301e-03, -5.400846550901379e+00,-3.836187701534530e-03, -5.400851138537755e+00,-1.771408421201963e-03, -5.400854782632413e+00,2.698687862826878e-04, -5.400857677243193e+00,2.314729424402266e-03, -5.400859976516478e+00,4.390306840056429e-03, -5.400861802895561e+00,6.524142180881721e-03, -5.400863253640913e+00,8.744549788870655e-03, -5.400864406009464e+00,1.108099287153084e-02, -5.400865321368692e+00,1.356447443033695e-02, -5.400866048464592e+00,1.622794863019467e-02, -5.400866626017535e+00,1.910675806594559e-02, -5.400867084784231e+00,2.223910272633751e-02, -5.400867449195627e+00,2.566654687702737e-02, -5.400867738657925e+00,2.943457058772546e-02, -5.400867968586022e+00,3.359317322125989e-02, -5.400868151224413e+00,3.819753689203313e-02, -5.400868296299255e+00,4.330875869723461e-02, -5.400868411536304e+00,4.899466143684154e-02, -5.400868503072349e+00,5.533069358000627e-02, -5.400868575782016e+00,6.240093041972123e-02, -5.400868633537358e+00,7.029918970039392e-02, -5.400868679414057e+00,7.913027652197359e-02, -5.400868715855216e+00,8.901137403971020e-02, -5.400868744801459e+00,1.000735984132572e-01, -5.400868767794277e+00,1.124637386383277e-01, -5.400868786058122e+00,1.263462043473980e-01, -5.400868800565608e+00,1.419052074255804e-01, -5.400868812089316e+00,1.593472063903755e-01, -5.400868821242920e+00,1.789036459707314e-01, -5.400868828513887e+00,2.008340282379587e-01, -5.400868834289422e+00,2.254293560405370e-01, -5.400868838877093e+00,2.530159944350999e-01, -5.400868842521208e+00,2.839600013524421e-01, -5.400868845415833e+00,3.186719849638514e-01, -5.400868847715115e+00,3.576125522020669e-01, -5.400868849541500e+00,4.012984207354255e-01, -5.400868850992248e+00,4.503092754973801e-01, -5.400868852144618e+00,5.052954607533686e-01, -5.400868853059978e+00,5.669866097740897e-01, -5.400868853787076e+00,6.362013266257051e-01, -5.400868854364630e+00,7.138580485484424e-01, -5.400868854823397e+00,8.009872330607631e-01, -5.400868855187808e+00,8.987450315045631e-01, -5.400868855477269e+00,1.008428630471043e+00, -5.400868855707198e+00,1.131493464678650e+00, -5.400868855889837e+00,1.269572529707519e+00, -5.400868856034912e+00,1.424498050858674e+00, -5.400868856150151e+00,1.598325795670576e+00, -5.400868856241686e+00,1.793362352705347e+00, -5.400868856314396e+00,2.012195738577593e+00, -5.400868856372150e+00,2.257729739362485e+00, -5.400868856418026e+00,2.533222442072750e+00, -5.400868856454468e+00,2.842329467493824e+00, -5.400868856483415e+00,3.189152478050652e+00, -5.400868856506408e+00,3.578293604376049e+00, -5.400868856524670e+00,4.014916512788097e+00, -5.400868856539178e+00,4.504814924004941e+00, -5.400868856550703e+00,5.054489492298646e+00, -5.400868856559858e+00,5.671234065227827e+00, -5.400868856567127e+00,6.363232468563275e+00, -5.400868856572901e+00,7.139667100683735e+00, -5.400868856577491e+00,8.010840777423491e+00, -5.400868856581135e+00,8.988313444178857e+00, -5.400868856584029e+00,1.008505556936019e+01, -5.400868856586328e+00,1.131562025462730e+01, -5.400868856588156e+00,1.269633634570651e+01, -5.400868856589605e+00,1.424552510625265e+01, -5.400868856590759e+00,1.598374332988642e+01,
            };
            var references = new Complex[1][];
            references[0] = new Complex[riref.Length / 2];
            for (var i = 0; i < riref.Length; i += 2)
                references[0][i / 2] = new Complex(riref[i], riref[i + 1]);

            // Run test
            AnalyzeAC(ac, ckt, exports, references);
        }
        /*
        [Test]
        public void When_BSIM2Netlist_Expect_Parameters()
        {
            // Create the parser
            var parser = new ParserFacade();
            var settings = new ParserSettings();
            settings.SpiceNetlistModelReaderSettings.EvaluatorMode = SpiceEvaluatorMode.Spice3f5;
            settings.SpiceNetlistParserSettings.HasTitle = false;
            settings.SpiceNetlistModelReaderSettings.Context.Models.Add(new MosfetModelGenerator(), true);
            settings.SpiceNetlistModelReaderSettings.Context.Components.Add(new MosfetGenerator(), true);

            var netlist = "M1 d g 0 0 nmodel w=10u l=40u\r\n"
                          + ".model nmodel nmos(level=5 vfb=-0.3 phi=0.8 k1=0.6 mu0=250 n0=1.3 tox=1e-7 mj=0.5 mjsw=0.33 pb=0.8 pbsw=1.0 xpart=1.0)\r\n";
            var result = parser.ParseNetlist(netlist, settings);

            // Find back the model
            var entity = (BSIM2) result.ReaderResult.Circuit.Objects["M1"];
            var model = (BSIM2Model) entity.Model;

            // Check a few component parameters
            Assert.AreEqual(entity.ParameterSets.GetParameter<double>("w"), 10e-6, 1e-12);
            Assert.AreEqual(entity.ParameterSets.GetParameter<double>("l"), 40e-6, 1e-12);

            // Check a few model parameters
            Assert.AreEqual(model.ParameterSets.GetParameter<double>("vfb"), -0.3, 1e-12);
            Assert.AreEqual(model.ParameterSets.GetParameter<double>("mu0"), 250.0, 1e-12);
        }
        */
    }
}
