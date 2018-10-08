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
        /// <summary>
        /// Generate a BSIM1 transistor
        /// </summary>
        private BSIM2 Create(string name, string drain, string gate, string source, string bulk, double w, double l, string model, string parameters)
        {
            // Create the model
            var m = new BSIM2Model(model);
            ApplyParameters(m, parameters);

            // Create the device
            var e = new BSIM2(name, drain, gate, source, bulk);
            e.SetParameter("w", w);
            e.SetParameter("l", l);
            e.SetModel(m);
            return e;
        }

        [Test]
        public void When_BSIM2DC_Expect_Reference()
        {
            var ckt = new Circuit();
            ckt.Entities.Add(
                new VoltageSource("V1", "g", "0", 0.0),
                new VoltageSource("V2", "d", "0", 0.0),
                Create("M1", "d", "g", "0", "0", 10e-6, 1e-6, "mod", "vfb=-0.3 phi=0.8 k1=0.6 mu0=250 n0=1.3 tox=1e-7 mj=0.5 mjsw=0.33 pb=0.8 pbsw=1.0 xpart=1.0"));
            ckt.Entities["M1"].SetParameter("m", 3.0);

            // Create simulation
            var dc = new DC("dc", new[]
            {
                new SweepConfiguration("V1", 0, 3.3, 0.3),
                new SweepConfiguration("V2", 0, 3.3, 0.3),
            });

            // Create exports
            Export<double>[] exports = {new RealPropertyExport(dc, "V2", "i")};

            // Create references
            double[][] references =
            {
                new[]
                {
                    -1.000000000000000e-50, -1.004372289037580e-02, -1.011033082222190e-02, -1.017684615804623e-02,
                    -1.024336149301214e-02, -1.030987682797807e-02, -1.037639216294399e-02, -1.044290749790992e-02,
                    -1.050942283287584e-02, -1.057593816784175e-02, -1.064245350280768e-02, -1.070896883777360e-02,
                    -1.000000000000000e-50, -8.761264572659984e-02, -8.819367500747333e-02, -8.877389656171586e-02,
                    -8.935411810847045e-02, -8.993433965522507e-02, -9.051456120197969e-02, -9.109478274873432e-02,
                    -9.167500429548893e-02, -9.225522584224355e-02, -9.283544738899814e-02, -9.341566893575277e-02,
                    -1.000000000000000e-50, -4.132952385071651e-01, -4.160340588741233e-01, -4.187711250678259e-01,
                    -4.215081912453394e-01, -4.242452574228530e-01, -4.269823236003666e-01, -4.297193897778800e-01,
                    -4.324564559553935e-01, -4.351935221329070e-01, -4.379305883104205e-01, -4.406676544879341e-01,
                    -1.000000000000000e-50, -3.693855480671697e+00, -3.728491233152088e+00, -3.753052508258300e+00,
                    -3.777582361859537e+00, -3.802112117763318e+00, -3.826641873363339e+00, -3.851171628962412e+00,
                    -3.875701384561484e+00, -3.900231140160553e+00, -3.924760895759626e+00, -3.949290651358698e+00,
                    -1.000000000000000e-50, -1.225705152720098e+01, -1.283812504301412e+01, -1.294805811829767e+01,
                    -1.303392674807220e+01, -1.311862294371964e+01, -1.320326215164724e+01, -1.328789858992737e+01,
                    -1.337253489360123e+01, -1.345717119073316e+01, -1.354180748754715e+01, -1.362644378434568e+01,
                    -1.000000000000000e-50, -2.135351687367902e+01, -2.607572954773885e+01, -2.657565568158140e+01,
                    -2.679209310611017e+01, -2.697155373962337e+01, -2.714626666321426e+01, -2.732037123206636e+01,
                    -2.749439786879719e+01, -2.766841452251566e+01, -2.784242989742756e+01, -2.801644510852665e+01,
                    -1.000000000000000e-50, -2.849261623302208e+01, -4.210242467285588e+01, -4.366106969397802e+01,
                    -4.422729254927948e+01, -4.457406484140301e+01, -4.487419169046089e+01, -4.516448727386884e+01,
                    -4.545271445713931e+01, -4.574050663308816e+01, -4.602820732943476e+01, -4.631588878845319e+01,
                    -1.000000000000000e-50, -3.434471761736317e+01, -5.840456067280783e+01, -6.325899838106288e+01,
                    -6.457838812172167e+01, -6.526063826518765e+01, -6.575684948678106e+01, -6.619977127232094e+01,
                    -6.662751239519153e+01, -6.705093574271756e+01, -6.747313154967910e+01, -6.789497841190590e+01,
                    -1.000000000000000e-50, -3.932606013516154e+01, -7.219553932770953e+01, -8.457593030515592e+01,
                    -8.713604187524258e+01, -8.842635278759738e+01, -8.925229433857037e+01, -8.991374874172371e+01,
                    -9.051759156512476e+01, -9.110133527978309e+01, -9.167807658957040e+01, -9.225237947204305e+01,
                    -1.000000000000000e-50, -4.367949544133837e+01, -8.391062588045278e+01, -1.062017823659480e+02,
                    -1.112859237155935e+02, -1.135184957103308e+02, -1.148723402430389e+02, -1.158612219168161e+02,
                    -1.167012615946594e+02, -1.174810469944070e+02, -1.182365130778147e+02, -1.189821754400764e+02,
                    -1.000000000000000e-50, -4.755341019079696e+01, -9.398512897954207e+01, -1.252828851089890e+02,
                    -1.365371186946263e+02, -1.400518932996493e+02, -1.421751455829641e+02, -1.436428774598817e+02,
                    -1.448095762202895e+02, -1.458396249551021e+02, -1.468079672312534e+02, -1.477485106835561e+02,
                    -1.000000000000000e-50, -5.104461398374103e+01, -1.027615924140216e+02, -1.419676270756469e+02,
                    -1.621616074087463e+02, -1.676215085527689e+02, -1.707698656276131e+02, -1.728987479429170e+02,
                    -1.745120187407013e+02, -1.758686711993065e+02, -1.770986148521515e+02, -1.782662379686020e+02
                }
            };
            for (var i = 0; i < references[0].Length; i++)
                references[0][i] *= 3.0;

            // Run simulation
            AnalyzeDC(dc, ckt, exports, references);
        }

                [Test]
        public void When_BSIM2SmallSignal_Expect_Reference()
        {
            // Build the circuit
            var ckt = new Circuit();
            ckt.Entities.Add(
                new VoltageSource("Vsupply", "vdd", "0", 5.0),
                new VoltageSource("V1", "in", "0", 0.0),
                new Resistor("R1", "vdd", "out", 10),
                new Resistor("R2", "out", "g", 1.0e6),
                new Capacitor("C1", "in", "g", 1e-9),
                Create("M1", "out", "g", "0", "0", 10e-6, 1e-6, "Nch4", "vfb=-0.3 phi=0.8 k1=0.6 mu0=250 n0=1.3 tox=1e-7 mj=0.5 mjsw=0.33 pb=0.8 pbsw=1.0 xpart=1.0")
            );
            ckt.Entities["V1"].SetParameter("acmag", 1.0);

            // AC simulation
            var ac = new AC("ac 1", new DecadeSweep(0.1, 1.0e9, 20));

            // Create exports
            var exports = new Export<Complex>[] { new ComplexVoltageExport(ac, "out") };

            // Reference
            var riref = new[]
            {
                -6.554049663587416e-08, -5.949569069478866e-04, -8.251059645725513e-08, -6.675526270124183e-04,
                -1.038746862114376e-07, -7.490063637628539e-04, -1.307704814624998e-07, -8.403989583464159e-04,
                -1.646302811935783e-07, -9.429431343291931e-04, -2.072572429096141e-07, -1.057999589724287e-03,
                -2.609214072847397e-07, -1.187095052519298e-03, -3.284805860029082e-07, -1.331942539275968e-03,
                -4.135325504878597e-07, -1.494464085602822e-03, -5.206066260920149e-07, -1.676816250094326e-03,
                -6.554048947776540e-07, -1.881418730121491e-03, -8.251058511241771e-07, -2.110986469136084e-03,
                -1.038746682310829e-06, -2.368565681459393e-03, -1.307704529655596e-06, -2.657574272478826e-03,
                -1.646302360289749e-06, -2.981847190457211e-03, -2.072571713285479e-06, -3.345687311538692e-03,
                -2.609212938364079e-06, -3.753922532871836e-03, -3.284804061994442e-06, -4.211969831027904e-03,
                -4.135322655186233e-06, -4.725907135146964e-03, -5.206061744463112e-06, -5.302553967700523e-03,
                -6.554041789676491e-06, -5.949561921760263e-03, -8.251047166421652e-06, -6.675516173706504e-03,
                -1.038744884278793e-05, -7.490049376065201e-03, -1.307701679968424e-05, -8.403969438480634e-03,
                -1.646297843843067e-05, -9.429402887764390e-03, -2.072564555206094e-05, -1.057995570277391e-02,
                -2.609201593585190e-05, -1.187089374905348e-02, -3.284786081756352e-05, -1.331934519443231e-02,
                -4.135294158478650e-05, -1.494452757305999e-02, -5.206016580323747e-05, -1.676800248481787e-02,
                -6.553970209535990e-05, -1.881396127299983e-02, -8.250933719936364e-05, -2.110954541903411e-02,
                -1.038726904300803e-04, -2.368520583224723e-02, -1.307673183779805e-04, -2.657510569849558e-02,
                -1.646252680739163e-04, -2.981757208671434e-02, -2.072492977131580e-04, -3.345560209901115e-02,
                -2.609088151222028e-04, -3.753742998838435e-02, -3.284606290201041e-04, -4.211716235668971e-02,
                -4.135009213002307e-04, -4.725548927878383e-02, -5.205564982025801e-04, -5.302047996615495e-02,
                -6.553254494115759e-04, -5.948847236620394e-02, -8.249799426639609e-04, -6.674506686146375e-02,
                -1.038547138749346e-03, -7.488623493944417e-02, -1.307388290184540e-03, -8.401955427730302e-02,
                -1.645801185948952e-03, -9.426558202047225e-02, -2.071777468214791e-03, -1.057593779758824e-01,
                -2.607954269899447e-03, -1.186521887641817e-01, -3.282809456546495e-03, -1.331133023590392e-01,
                -4.132161916533236e-03, -1.493320794256159e-01, -5.201053304462840e-03, -1.675201628034243e-01,
                -6.546105927971507e-03, -1.879138584443301e-01, -8.238473623132691e-03, -2.107766688322551e-01,
                -1.036752899513480e-02, -2.364019415967118e-01, -1.304546166844635e-02, -2.651155692357017e-01,
                -1.641299821563789e-02, -2.972786372420657e-01, -2.064649458066681e-02, -3.332898629824007e-01,
                -2.596669427097513e-02, -3.735875904559280e-01, -3.264948656474952e-02, -4.186509988635825e-01,
                -4.103903142909593e-02, -4.690000361106582e-01, -5.156362964083230e-02, -5.251933904778892e-01,
                -6.475468785812684e-02, -5.878235521538989e-01, -8.126902784501369e-02, -6.575076794427603e-01,
                -1.019145694941821e-01, -7.348725656961298e-01, -1.276790105495226e-01, -8.205315090124319e-01,
                -1.597604388853222e-01, -9.150503513424452e-01, -1.995977283798075e-01, -1.018899460400311e+00,
                -2.488969569567266e-01, -1.132388142916550e+00, -3.096478381169689e-01, -1.255578192343050e+00,
                -3.841212664376883e-01, -1.388174445870686e+00, -4.748357683786569e-01, -1.529392993573990e+00,
                -5.844775869279566e-01, -1.677812798839559e+00, -7.157576596250684e-01, -1.831224578579483e+00,
                -8.711910638160524e-01, -1.986501834076553e+00, -1.052793438905793e+00, -2.139531426280441e+00,
                -1.261706995765219e+00, -2.285251160609769e+00, -1.497796231995045e+00, -2.417842842429484e+00,
                -1.759285591961271e+00, -2.531113178516503e+00, -2.042536351336446e+00, -2.619056501149458e+00,
                -2.342061293977578e+00, -2.676536351623883e+00, -2.650839494281197e+00, -2.699964727193392e+00,
                -2.960920466273266e+00, -2.687825740705734e+00, -3.264219509601063e+00, -2.640908512252703e+00,
                -3.553341588147002e+00, -2.562185078336087e+00, -3.822261058429944e+00, -2.456367278745926e+00,
                -4.066734153793851e+00, -2.329260217930071e+00, -4.284405516082215e+00, -2.187065929515503e+00,
                -4.474651148324704e+00, -2.035772899280658e+00, -4.638249229714478e+00, -1.880714044137697e+00,
                -4.776979884486195e+00, -1.726315864418592e+00, -4.893235646119626e+00, -1.576016161606568e+00,
                -4.989692899255917e+00, -1.432304778116617e+00, -5.069064766901131e+00, -1.296838209238940e+00,
                -5.133934650533771e+00, -1.170586950574075e+00, -5.186657987505918e+00, -1.053986729212051e+00,
                -5.229315685369721e+00, -9.470765348801146e-01, -5.263703232487290e+00, -8.496154522131409e-01,
                -5.291342223784380e+00, -7.611763097878562e-01, -5.313504371005474e+00, -6.812175962162471e-01,
                -5.331241131869097e+00, -6.091366972943786e-01, -5.345414551979412e+00, -5.443079696243484e-01,
                -5.356726719419926e+00, -4.861090017842805e-01, -5.365746471111772e+00, -4.339379642544191e-01,
                -5.372932790235423e+00, -3.872244150907445e-01, -5.378654818836823e+00, -3.454354156100020e-01,
                -5.383208680482504e+00, -3.080783659011887e-01, -5.386831440318249e+00, -2.747016072544146e-01,
                -5.389712577481559e+00, -2.448935541484482e-01, -5.392003343335814e+00, -2.182809016115179e-01,
                -5.393824351508360e+00, -1.945262919782159e-01, -5.395271706498331e+00, -1.733257063612588e-01,
                -5.396421935098457e+00, -1.544057602874095e-01, -5.397335943718941e+00, -1.375210216074555e-01,
                -5.398062187248925e+00, -1.224514254998068e-01, -5.398639202284399e+00, -1.089998312162933e-01,
                -5.399097629536397e+00, -9.698974448624577e-02, -5.399461826732780e+00, -8.626321549067227e-02,
                -5.399751153865731e+00, -7.667891306823554e-02, -5.399980996674267e+00, -6.811036989361283e-02,
                -5.400163581251189e+00, -6.044438975532491e-02, -5.400308622135189e+00, -5.357960602717618e-02,
                -5.400423837757001e+00, -4.742517946064671e-02, -5.400515360282073e+00, -4.189962316329521e-02,
                -5.400588061418585e+00, -3.692974281793922e-02, -5.400645811378946e+00, -3.244968066409298e-02,
                -5.400691684683188e+00, -2.840005238673101e-02, -5.400728123699333e+00, -2.472716675636516e-02,
                -5.400757068589124e+00, -2.138231858415520e-02, -5.400780060553437e+00, -1.832114626058606e-02,
                -5.400798323859373e+00, -1.550304581336610e-02, -5.400812831006961e+00, -1.289063403569893e-02,
                -5.400824354499436e+00, -1.044925379197976e-02, -5.400833507969915e+00, -8.146515100192108e-03,
                -5.400840778852073e+00, -5.951866017260194e-03, -5.400846554333015e+00, -3.836187715263154e-03,
                -5.400851141969399e+00, -1.771408433472475e-03, -5.400854786064063e+00, 2.698687753074774e-04,
                -5.400857680674849e+00, 2.314729414576728e-03, -5.400859979948137e+00, 4.390306831250178e-03,
                -5.400861806327225e+00, 6.524142172977919e-03, -5.400863257072579e+00, 8.744549781764425e-03,
                -5.400864409441133e+00, 1.108099286512788e-02, -5.400865324800362e+00, 1.356447442455230e-02,
                -5.400866051896264e+00, 1.622794862495158e-02, -5.400866629449205e+00, 1.910675806117447e-02,
                -5.400867088215905e+00, 2.223910272197508e-02, -5.400867452627302e+00, 2.566654687301572e-02,
                -5.400867742089599e+00, 2.943457058401135e-02, -5.400867972017696e+00, 3.359317321779406e-02,
                -5.400868154656090e+00, 3.819753688876958e-02, -5.400868299730932e+00, 4.330875869413003e-02,
                -5.400868414967980e+00, 4.899466143385476e-02, -5.400868506504025e+00, 5.533069357709764e-02,
                -5.400868579213691e+00, 6.240093041685213e-02, -5.400868636969034e+00, 7.029918969752630e-02,
                -5.400868682845734e+00, 7.913027651906941e-02, -5.400868719286893e+00, 8.901137403673090e-02,
                -5.400868748233136e+00, 1.000735984101633e-01, -5.400868771225953e+00, 1.124637386350781e-01,
                -5.400868789489799e+00, 1.263462043439496e-01, -5.400868803997285e+00, 1.419052074218874e-01,
                -5.400868815520992e+00, 1.593472063863888e-01, -5.400868824674597e+00, 1.789036459663983e-01,
                -5.400868831945564e+00, 2.008340282332217e-01, -5.400868837721100e+00, 2.254293560353331e-01,
                -5.400868842308770e+00, 2.530159944293601e-01, -5.400868845952885e+00, 2.839600013460902e-01,
                -5.400868848847511e+00, 3.186719849568032e-01, -5.400868851146791e+00, 3.576125521942287e-01,
                -5.400868852973176e+00, 4.012984207266936e-01, -5.400868854423925e+00, 4.503092754876385e-01,
                -5.400868855576296e+00, 5.052954607424880e-01, -5.400868856491656e+00, 5.669866097619258e-01,
                -5.400868857218754e+00, 6.362013266120964e-01, -5.400868857796307e+00, 7.138580485332082e-01,
                -5.400868858255073e+00, 8.009872330437015e-01, -5.400868858619485e+00, 8.987450314854474e-01,
                -5.400868858908948e+00, 1.008428630449620e+00, -5.400868859138876e+00, 1.131493464654635e+00,
                -5.400868859321514e+00, 1.269572529680594e+00, -5.400868859466590e+00, 1.424498050828481e+00,
                -5.400868859581826e+00, 1.598325795636714e+00, -5.400868859673362e+00, 1.793362352667367e+00,
                -5.400868859746073e+00, 2.012195738534992e+00, -5.400868859803826e+00, 2.257729739314697e+00,
                -5.400868859849704e+00, 2.533222442019141e+00, -5.400868859886144e+00, 2.842329467433681e+00,
                -5.400868859915090e+00, 3.189152477983179e+00, -5.400868859938083e+00, 3.578293604300350e+00,
                -5.400868859956347e+00, 4.014916512703168e+00, -5.400868859970855e+00, 4.504814923909656e+00,
                -5.400868859982380e+00, 5.054489492191738e+00, -5.400868859991534e+00, 5.671234065107879e+00,
                -5.400868859998805e+00, 6.363232468428696e+00, -5.400868860004578e+00, 7.139667100532737e+00,
                -5.400868860009167e+00, 8.010840777254073e+00, -5.400868860012811e+00, 8.988313443988769e+00,
                -5.400868860015705e+00, 1.008505556914691e+01, -5.400868860018004e+00, 1.131562025438800e+01,
                -5.400868860019831e+00, 1.269633634543801e+01, -5.400868860021282e+00, 1.424552510595140e+01,
                -5.400868860022435e+00, 1.598374332954840e+01
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
