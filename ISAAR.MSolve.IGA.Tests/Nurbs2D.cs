﻿using System.Collections.Generic;
using System.Linq;
using ISAAR.MSolve.IGA.Elements;
using ISAAR.MSolve.IGA.Entities;
using ISAAR.MSolve.IGA.SupportiveClasses;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Materials;
using Xunit;

namespace ISAAR.MSolve.IGA.Tests
{
    public class Nurbs2DTests
	{
		private List<ControlPoint> ElementControlPoints()
		{
			return new List<ControlPoint>
			{
				new ControlPoint {ID = 0, X = 0.0, Y =  0.0, Z = 0.0, WeightFactor = 1.0},
				new ControlPoint {ID = 1, X = 0.0, Y =  0.037037037, Z = 0.0, WeightFactor = 1.0},
				new ControlPoint {ID = 2, X = 0.0, Y =  0.111111111, Z = 0.0, WeightFactor = 1.0},
				new ControlPoint {ID = 3, X = 0.0, Y =  0.222222222, Z = 0.0, WeightFactor = 1.0},
				new ControlPoint {ID = 12, X = 0.037037037, Y =  0.0, Z = 0.0, WeightFactor = 1.0},
				new ControlPoint {ID = 13, X = 0.037037037, Y =  0.037037037, Z = 0.0, WeightFactor = 1.0},
				new ControlPoint {ID = 14, X = 0.037037037, Y =  0.111111111, Z = 0.0, WeightFactor = 1.0},
				new ControlPoint {ID = 15, X = 0.037037037, Y =  0.222222222, Z = 0.0, WeightFactor = 1.0},
				new ControlPoint {ID = 24, X = 0.111111111, Y =  0.0, Z = 0.0, WeightFactor = 1.0},
				new ControlPoint {ID = 25, X = 0.111111111, Y =  0.037037037, Z = 0.0, WeightFactor = 1.0},
				new ControlPoint {ID = 26, X = 0.111111111, Y =  0.111111111, Z = 0.0, WeightFactor = 1.0},
				new ControlPoint {ID = 27, X = 0.111111111, Y =  0.222222222, Z = 0.0, WeightFactor = 1.0},
				new ControlPoint {ID = 36, X = 0.222222222, Y =  0.0, Z = 0.0, WeightFactor = 1.0},
				new ControlPoint {ID = 37, X = 0.222222222, Y =  0.037037037, Z = 0.0, WeightFactor = 1.0},
				new ControlPoint {ID = 38, X = 0.222222222, Y =  0.111111111, Z = 0.0, WeightFactor = 1.0},
				new ControlPoint {ID = 39, X = 0.222222222, Y =  0.222222222, Z = 0.0, WeightFactor = 1.0},
			};
		}

		private List<Knot> ElementKnot()
		{
			return new List<Knot>()
			{
				new Knot(){ID=0,Ksi=0.0,Heta=0.0,Zeta =0.0 },
				new Knot(){ID=1,Ksi=0.0,Heta=0.111111111,Zeta =0.0 },
				new Knot(){ID=2,Ksi=0.111111111,Heta=0.0,Zeta =0.0 },
				new Knot(){ID=3,Ksi=0.111111111,Heta=0.111111111,Zeta =0.0 }
			};
		}

		private double[] KnotValueVector()
		{
			return new double[16]
			{
				0.0, 0.0, 0.0, 0.0, 0.111111111, 0.22222222, 0.33333333, 0.44444444, 0.55555555, 0.66666666, 0.77777777,
				0.88888888, 1.0, 1.0, 1.0, 1.0
			};
		}

		private NURBSElement2D Element
		{
			get
			{
                var degreeKsi = 3;
                var degreeHeta = 3;
                var numberOfControlPointsHeta = 12;
                var knotValueVectorKsi = KnotValueVector();
                var knotValueVectorHeta = KnotValueVector();
				var gauss = new GaussQuadrature();
                var parametricPointKsi = gauss.CalculateElementGaussPoints(degreeKsi, new List<Knot>()
                {
					new Knot(){Ksi = ElementKnot()[0].Ksi},
					new Knot(){Ksi = ElementKnot()[2].Ksi},
                }).Select(x => x.Ksi).ToArray();
                var parametricPointHeta = gauss.CalculateElementGaussPoints(degreeHeta, new List<Knot>()
                {
					new Knot(){Ksi = ElementKnot()[0].Heta},
                    new Knot(){Ksi = ElementKnot()[1].Heta},
				}).Select(x => x.Ksi).ToArray();
                var gaussPoints = gauss.CalculateElementGaussPoints(degreeKsi, degreeHeta, ElementKnot()).ToArray();

                var nurbs = new Nurbs2D(degreeKsi, knotValueVectorKsi, degreeHeta, knotValueVectorHeta,
                    ElementControlPoints().ToArray(), parametricPointKsi, parametricPointHeta);


				var element = new NURBSElement2D(null,nurbs, gaussPoints, 1);
				var patch = new Patch();
				foreach (var controlPoint in ElementControlPoints())
					element.ControlPointsDictionary.Add(controlPoint.ID, controlPoint);
				foreach (var knot in ElementKnot())
					element.KnotsDictionary.Add(knot.ID, knot);
				element.Patch = patch;
				return element;
			}
		}

		private readonly Matrix _nurbsDerivativeKsiExpectedValues= Matrix.CreateFromArray(new double[16,16]
		{
			{ -18.841032499144923, -7.0317936031190476, -0.8403097033832355, -0.00782593507823648, -9.766650411253282, -3.645079954554914, -0.4355924289446699, -0.0040567400992365275, -2.3695242185857377, -0.8843467173809627, -0.10568073662468015, -9.842211514577343E-4, -0.1048879812532743, -0.03914597757072591, -0.004678002037277168, -4.3566961195606604E-5},
			{ -4.376607599350048, -13.1598003124014, -12.071332733959101, -7.133581537359361, -2.2687077479446156, -6.821662726761882, -6.2574323788172155, -3.697843897818317, -0.5504198192098536, -1.6550295507107928, -1.5181394790875182, -0.8971479784240177, -0.024364563664671737, -0.07326057574215399, -0.06720107942911306, -0.03971263074851178 },
			{ -0.16189710210327568, -3.0491959878628467, -9.297233464177204, -13.099261910983246, -0.08392281043565063, -1.5806156722145708, -4.8194189485393135, -6.790281357235874, -0.02036083236786359, -0.3834791825120721, -1.169257569096251, -1.6474159972433124, -9.012807663878547E-4, -0.016974866511537113, -0.051757675665714695, -0.07292355861141347},
			{ -0.001304322536517215, -0.1400516197514684, -1.1719656216152226, -3.1401721397139175, -6.76123362042973E-4, -0.07259873946422228, -0.6075133366943926, -1.6277750978421652, -1.640368615289521E-4, -0.017613456421155572, -0.147391122216534, -0.39492071020619557, -7.261160329968631E-6, -7.796670202468678E-4, -0.006524329712558918, -0.017481330523542993},
			{ 17.409037442795718, 6.497348705930349, 0.7764427501734436, 0.007231132200854201, 4.363623339647517, 1.6285784065999005, 0.1946175207957304, 0.0018125032671913687, -4.882784401186915, -1.8223381398675393, -0.2177720946042568, -0.002028145417532123, -6.011143204568842, -2.2434608260050535, -0.2680968765090044, -0.0024968279454468927},
			{ 4.043967631443085, 12.15960199573501, 11.153862377737367, 6.591400343413321, 1.0136316611029301, 3.0478378396142336, 2.79574642530581, 1.6521527094182684, -1.1342282498155225, -3.4104536304496285, -3.1283696993563113, -1.848717189824605, -1.3963361631636362, -4.19857267508905, -3.8513021904251175, -2.275935789850242},
			{ 0.1495922642521826, 2.817444698189227, 8.59060592882157, 12.103664758905582, 0.03749571438779763, 0.7062002987576046, 2.153259113600955, 3.033817016637299, -0.04195675819693619, -0.7902203133705189, -2.4094426107028735, -3.3947646833534026, -0.05165251241892718, -0.9728317034999107, -2.9662388068177896, -4.179258181559211},
			{ 0.00120518872183576, 0.12940712735823554, 1.0828914704804382, 2.9015062926930653, 3.020838833027371E-4, 0.03243625404980762, 0.2714297393190522, 0.7272705696987896, -3.380242423397893E-4, -0.036295349754025756, -0.30372302877827106, -0.8137974148461722, -4.161379983982982E-4, -0.044682813555791394, -0.37391014439789144, -1.0018572187949033},
			{ 1.4145137258256606, 0.5279205674761392, 0.06308728618954511, 5.875417170523119E-4, 5.008106361399569, 1.8691104257384794, 0.22336145172778393, 0.0020801999705162065, 5.624533521930485, 2.0991715205541093, 0.2508540917647146, 0.0023362432069468824, 2.975859046108198, 1.1106411819605566, 0.1327232587948132, 0.0012360723701252802},
			{ 0.32857920722309125, 0.9879882204895831, 0.906270176115278, 0.5355624220354847, 1.1633394486554502, 3.4979866237307604, 3.2086627024423997, 1.8961665226379345, 1.3065301042284831, 3.928539372901785, 3.6036037632019955, 2.1295578408520384, 0.6912661138068102, 2.078533159284901, 1.9066144446470485, 1.1267181429714486},
			{ 0.012154624387324581, 0.22892214520410584, 0.6980012559227395, 0.9834432254236596, 0.0430336239254599, 0.8105021752211832, 2.471283569914544, 3.4818950027824416, 0.048330455240422884, 0.9102635457716486, 2.7754636739176477, 3.910467100686559, 0.025570942349077776, 0.481607229553385, 1.468457543895603, 2.06897138237563},
			{ 9.792362127134351E-5, 0.01051454788751938, 0.08798676282978495, 0.2357522918811513, 3.467000013266334E-4, 0.03722690929138168, 0.311518409897078, 0.834684408590914, 3.893738748332089E-4, 0.041809016046681385, 0.3498619263898662, 0.9374222705286817, 2.0601206539576331E-4, 0.02212054353063865, 0.18510686699201637, 0.4959765166122767},
			{ 0.017481330523542962, 0.006524329712558904, 7.796670202468663E-4, 7.261160329968639E-6, 0.3949207102061956, 0.14739112221653394, 0.017613456421155572, 1.6403686152895253E-4, 1.627775097842166, 0.6075133366943926, 0.0725987394642223, 6.761233620429752E-4, 3.140172139713919, 1.1719656216152226, 0.14005161975146838, 0.001304322536517219},
			{ 0.004060760683872307, 0.01221009617680737, 0.011200180106455395, 0.006618771910556511, 0.09173663818623505, 0.27583826341688744, 0.25302325106900486, 0.1495246657621136, 0.37811796479689275, 1.1369438082586358, 1.0429054152418333, 0.6163073273965847, 0.7294346130214976, 2.193300091546302, 2.011888825207182, 1.1889302776273052},
			{ 1.502134637684845E-4, 0.002829144469514113, 0.008626279432892129, 0.012153926654006235, 0.003393472122393095, 0.0639131982357829, 0.1948762650238146, 0.27456933781613346, 0.013987135324376879, 0.2634359501109423, 0.8032365058814762, 1.1317135799101568, 0.026982850836237263, 0.5081993404580626, 1.5495389385879013, 2.183210357794994},
			{ 1.210193410111584E-6, 1.299445057134788E-4, 0.0010873883049994737, 0.002913555139701152, 2.73394774136026E-5, 0.0029355761230329645, 0.024565187478262373, 0.0658201195524614, 1.1268722903553233E-4, 0.012099790128499937, 0.10125222460493882, 0.2712958545236862, 2.1738709333250355E-4, 0.0233419370453996, 0.195327607118434, 0.5233620327061695}
		});

		private readonly Matrix _nurbsDerivativeHetaExpectedValues = Matrix.CreateFromArray(new double[16, 16]
		{
			{-18.841032499144923, -9.766650411253282, -2.3695242185857373, -0.1048879812532743, -7.0317936031190476, -3.645079954554914, -0.8843467173809627, -0.039145977570725896, -0.8403097033832356, -0.4355924289446699, -0.10568073662468015, -0.004678002037277168, -0.00782593507823648, -0.004056740099236528, -9.842211514577343E-4, -4.3566961195606604E-5},
			{17.409037442795718, 4.363623339647517, -4.882784401186914, -6.011143204568842, 6.497348705930349, 1.6285784065999003, -1.8223381398675391, -2.243460826005053, 0.7764427501734437, 0.19461752079573036, -0.21777209460425678, -0.26809687650900443, 0.007231132200854201, 0.0018125032671913689, -0.0020281454175321233, -0.0024968279454468927},
			{1.4145137258256606, 5.008106361399569, 5.624533521930485, 2.9758590461081975, 0.5279205674761392, 1.8691104257384794, 2.0991715205541097, 1.1106411819605566, 0.06308728618954512, 0.2233614517277839, 0.2508540917647147, 0.1327232587948132, 5.875417170523119E-4, 0.0020801999705162065, 0.0023362432069468824, 0.0012360723701252804},
			{0.017481330523542962, 0.3949207102061956, 1.6277750978421657, 3.140172139713919, 0.006524329712558904, 0.14739112221653394, 0.6075133366943926, 1.1719656216152223, 7.796670202468664E-4, 0.017613456421155572, 0.0725987394642223, 0.14005161975146838, 7.261160329968639E-6, 1.640368615289526E-4, 6.761233620429752E-4, 0.001304322536517219},
			{-4.376607599350048, -2.2687077479446156, -0.5504198192098535, -0.024364563664671737, -13.1598003124014, -6.821662726761882, -1.655029550710793, -0.07326057574215397, -12.071332733959103, -6.257432378817215, -1.5181394790875182, -0.06720107942911306, -7.133581537359363, -3.697843897818317, -0.8971479784240176, -0.03971263074851178},
			{4.043967631443085, 1.0136316611029301, -1.1342282498155223, -1.3963361631636362, 12.15960199573501, 3.0478378396142327, -3.410453630449628, -4.198572675089049, 11.15386237773737, 2.79574642530581, -3.128369699356311, -3.851302190425118, 6.591400343413321, 1.6521527094182684, -1.848717189824605, -2.275935789850242},
			{0.32857920722309125, 1.1633394486554502, 1.3065301042284831, 0.6912661138068101, 0.9879882204895831, 3.4979866237307604, 3.9285393729017857, 2.078533159284901, 0.9062701761152782, 3.2086627024423993, 3.603603763201996, 1.9066144446470488, 0.5355624220354847, 1.8961665226379347, 2.129557840852038, 1.1267181429714488},
			{0.004060760683872307, 0.09173663818623505, 0.3781179647968927, 0.7294346130214976, 0.01221009617680737, 0.27583826341688744, 1.1369438082586358, 2.1933000915463015, 0.011200180106455399, 0.2530232510690048, 1.0429054152418333, 2.011888825207182, 0.006618771910556511, 0.1495246657621136, 0.6163073273965847, 1.1889302776273052},
			{-0.16189710210327568, -0.08392281043565063, -0.020360832367863582, -9.012807663878547E-4, -3.0491959878628467, -1.5806156722145708, -0.38347918251207214, -0.01697486651153711, -9.297233464177205, -4.8194189485393135, -1.169257569096251, -0.0517576756657147, -13.099261910983248, -6.790281357235874, -1.6474159972433122, -0.07292355861141347},
			{0.1495922642521826, 0.03749571438779763, -0.04195675819693618, -0.05165251241892718, 2.817444698189227, 0.7062002987576045, -0.7902203133705189, -0.9728317034999104, 8.590605928821573, 2.153259113600955, -2.409442610702873, -2.9662388068177896, 12.103664758905582, 3.0338170166372986, -3.3947646833534026, -4.179258181559211},
			{0.012154624387324581, 0.0430336239254599, 0.04833045524042289, 0.025570942349077773, 0.22892214520410584, 0.8105021752211832, 0.9102635457716488, 0.481607229553385, 0.6980012559227398, 2.4712835699145437, 2.775463673917648, 1.468457543895603, 0.9834432254236596, 3.481895002782441, 3.9104671006865583, 2.0689713823756306},
			{1.502134637684845E-4, 0.003393472122393095, 0.013987135324376876, 0.026982850836237263, 0.002829144469514113, 0.0639131982357829, 0.2634359501109423, 0.5081993404580626, 0.00862627943289213, 0.19487626502381458, 0.8032365058814762, 1.5495389385879015, 0.012153926654006235, 0.27456933781613346, 1.1317135799101568, 2.183210357794994},
			{-0.001304322536517215, -6.76123362042973E-4, -1.6403686152895205E-4, -7.261160329968631E-6, -0.1400516197514684, -0.07259873946422228, -0.017613456421155572, -7.796670202468677E-4, -1.1719656216152228, -0.6075133366943926, -0.147391122216534, -0.006524329712558918, -3.1401721397139184, -1.627775097842165, -0.3949207102061955, -0.017481330523542993},
			{0.00120518872183576, 3.020838833027371E-4, -3.380242423397892E-4, -4.161379983982982E-4, 0.12940712735823554, 0.03243625404980762, -0.036295349754025756, -0.04468281355579138, 1.0828914704804384, 0.2714297393190521, -0.30372302877827106, -0.37391014439789144, 2.9015062926930653, 0.7272705696987897, -0.8137974148461722, -1.0018572187949033},
			{9.792362127134351E-5, 3.467000013266334E-4, 3.8937387483320895E-4, 2.060120653957633E-4, 0.01051454788751938, 0.03722690929138168, 0.04180901604668139, 0.02212054353063865, 0.08798676282978497, 0.311518409897078, 0.34986192638986624, 0.18510686699201637, 0.2357522918811513, 0.834684408590914, 0.9374222705286815, 0.49597651661227676},
			{1.210193410111584E-6, 2.73394774136026E-5, 1.1268722903553231E-4, 2.1738709333250355E-4, 1.299445057134788E-4, 0.0029355761230329645, 0.012099790128499939, 0.023341937045399597, 0.001087388304999474, 0.02456518747826237, 0.10125222460493882, 0.19532760711843403, 0.002913555139701152, 0.0658201195524614, 0.2712958545236862, 0.5233620327061695}
		});

		private readonly Matrix _nurbsExpectedValues = Matrix.CreateFromArray(new double[16, 16]
		{
			{0.6493653647595652, 0.24235419254282428, 0.028961683340508082, 2.697246654021931E-4, 0.24235419254282428, 0.09045070438093325, 0.01080899253576399, 1.0066583011651002E-4, 0.028961683340508085, 0.010808992535763989, 0.001291690545008403, 1.2029715122541487E-5, 2.697246654021931E-4, 1.0066583011651004E-4, 1.2029715122541485E-5, 1.1203460959649736E-7},
			{0.15084191327043306, 0.45355892944897813, 0.41604436403775547, 0.24586236328925679, 0.05629676615999721, 0.16927590241572968, 0.1552748289026726, 0.09176001333810073, 0.006727546561144669, 0.020228719922296495, 0.018555570994036534, 0.010965456886615581, 6.265468839804624E-5, 1.8839321763195488E-4, 1.7281092120473892E-4, 1.0212300697272503E-4},
			{0.005579862503054512, 0.1050920253426569, 0.3204336811156322, 0.45147244394874647, 0.002082499536957078, 0.03922213029338267, 0.11959129677182241, 0.16849719055934947, 2.4886176514337906E-4, 0.00468710239991568, 0.014291336291920061, 0.020135662707311455, 2.317688358977807E-6, 4.365171428950837E-5, 1.3309743961245632E-4, 1.8752656130262488E-4},
			{4.495411170953938E-5, 0.004826947310303064, 0.040392366150866645, 0.10822756274135698, 1.677763865474911E-5, 0.0018014987883876426, 0.015075117668174333, 0.04039236615086663, 2.004952556512721E-6, 2.1528176137647218E-4, 0.0018014987883876433, 0.0048269473103030635, 1.867243526885334E-8, 2.0049525565127274E-6, 1.6777638654749164E-5, 4.4954111709539504E-5},
			{0.15084191327043306, 0.05629676615999721, 0.006727546561144668, 6.265468839804624E-5, 0.45355892944897813, 0.16927590241572968, 0.0202287199222965, 1.8839321763195485E-4, 0.4160443640377555, 0.15527482890267258, 0.018555570994036534, 1.7281092120473895E-4, 0.24586236328925679, 0.09176001333810073, 0.01096545688661558, 1.0212300697272503E-4},
			{0.0350392614603175, 0.1053577853267627, 0.09664347882192954, 0.05711168364126321, 0.1053577853267627, 0.3167950027009396, 0.29059239466230824, 0.17172623662575537, 0.09664347882192956, 0.2905923946623082, 0.2665570451415591, 0.15752249215412392, 0.05711168364126321, 0.17172623662575537, 0.1575224921541239, 0.09308827504922461},
			{0.0012961534159715011, 0.024411961327843054, 0.07443395068301234, 0.10487311293438388, 0.003897338233145458, 0.07340309338131852, 0.22381168638400276, 0.3153376657668297, 0.0035749833183064566, 0.06733180921239043, 0.20529987324684518, 0.28925559633562303, 0.0021126445238391234, 0.03978988581130588, 0.12132242708345792, 0.1709362526138337},
			{1.0442448254305073E-5, 0.0011212577803734179, 0.009382794528889838, 0.025140319330931256, 3.1398870170505374E-5, 0.003371453381068909, 0.02821265091044934, 0.07559315626883982, 2.8801820719222513E-5, 0.003092595221339469, 0.025879138616269927, 0.06934072858775901, 1.7020501468489817E-5, 0.0018275761806656348, 0.015293335831630166, 0.04097706128579655},
			{0.005579862503054512, 0.002082499536957078, 2.4886176514337906E-4, 2.317688358977807E-6, 0.1050920253426569, 0.03922213029338267, 0.004687102399915681, 4.365171428950837E-5, 0.3204336811156322, 0.1195912967718224, 0.014291336291920061, 1.3309743961245634E-4, 0.45147244394874647, 0.16849719055934947, 0.02013566270731145, 1.8752656130262488E-4},
			{0.0012961534159715011, 0.003897338233145458, 0.003574983318306456, 0.0021126445238391234, 0.024411961327843054, 0.07340309338131852, 0.06733180921239044, 0.03978988581130588, 0.07443395068301234, 0.22381168638400273, 0.20529987324684518, 0.12132242708345793, 0.10487311293438388, 0.3153376657668297, 0.289255596335623, 0.1709362526138337},
			{4.7946606398574696E-5, 9.030340751183527E-4, 0.002753420460967857, 0.0038794094940445972, 9.030340751183527E-4, 0.01700788860938237, 0.051858362585927895, 0.07306542063348324, 0.002753420460967857, 0.05185836258592789, 0.15812014247377087, 0.2227820961629703, 0.0038794094940445972, 0.07306542063348324, 0.22278209616297026, 0.3138870329502765},
			{3.862814001160318E-7, 4.1476961603815144E-5, 3.4708326240701174E-4, 9.299771005820063E-4, 7.275285845873183E-6, 7.811837473805872E-4, 0.006537021833230174, 0.017515337872385562, 2.2182907001368314E-5, 0.002381889424860686, 0.019931883154077634, 0.05340561448057308, 3.125442744635042E-5, 0.0033559438449588985, 0.028082865598716488, 0.0752454086792084},
			{4.495411170953938E-5, 1.677763865474911E-5, 2.004952556512721E-6, 1.867243526885334E-8, 0.004826947310303064, 0.0018014987883876426, 2.152817613764722E-4, 2.0049525565127274E-6, 0.040392366150866645, 0.01507511766817433, 0.0018014987883876433, 1.6777638654749168E-5, 0.10822756274135698, 0.04039236615086663, 0.004826947310303063, 4.4954111709539504E-5},
			{1.0442448254305073E-5, 3.1398870170505374E-5, 2.880182071922251E-5, 1.7020501468489817E-5, 0.0011212577803734179, 0.003371453381068909, 0.0030925952213394694, 0.0018275761806656346, 0.009382794528889838, 0.028212650910449336, 0.025879138616269927, 0.015293335831630167, 0.025140319330931256, 0.07559315626883982, 0.069340728587759, 0.04097706128579655},
			{3.862814001160318E-7, 7.275285845873183E-6, 2.2182907001368314E-5, 3.125442744635042E-5, 4.1476961603815144E-5, 7.811837473805872E-4, 0.0023818894248606864, 0.003355943844958898, 3.4708326240701174E-4, 0.006537021833230173, 0.019931883154077634, 0.02808286559871649, 9.299771005820063E-4, 0.017515337872385562, 0.05340561448057307, 0.0752454086792084},
			{3.1120726008261867E-9, 3.341587654336446E-7, 2.796273159457769E-6, 7.492352086452233E-6, 3.341587654336446E-7, 3.588029420859065E-5, 3.002498034691037E-4, 8.044912328646861E-4, 2.796273159457769E-6, 3.0024980346910363E-4, 0.002512519656587743, 0.006732061146321539, 7.492352086452233E-6, 8.044912328646863E-4, 0.006732061146321538, 0.01803792744824218}
		});
		
		private const double Tolerance = 1e-14;

		#region ShellElement
		private List<ControlPoint> ShellElementControlPoints()
		{
			return new List<ControlPoint>
			{
				new ControlPoint {ID = 0, X = 0.0, Y =  0.0, Z = 0.0, WeightFactor = 1.0},
				new ControlPoint {ID = 1, X = 0.0, Y =  0.5, Z = 0.0, WeightFactor = 1.0},
				new ControlPoint {ID = 2, X = 0.0, Y =  1.0, Z = 0.0, WeightFactor = 1.0},
				new ControlPoint {ID = 3, X = 16.66666667, Y =  0.0, Z = 0.0, WeightFactor = 1.0},
				new ControlPoint {ID = 4, X = 16.66666667, Y =  0.5, Z = 0.0, WeightFactor = 1.0},
				new ControlPoint {ID = 5, X = 16.66666667, Y =  1.0, Z = 0.0, WeightFactor = 1.0},
				new ControlPoint {ID = 6, X = 33.33333333, Y =  0.0, Z = 0.0, WeightFactor = 1.0},
				new ControlPoint {ID = 7, X = 33.33333333, Y =  0.5, Z = 0.0, WeightFactor = 1.0},
				new ControlPoint {ID = 8, X = 33.33333333, Y =  1.0, Z = 0.0, WeightFactor = 1.0},
				new ControlPoint {ID = 9, X = 50.0, Y =  0.0, Z = 0.0, WeightFactor = 1.0},
				new ControlPoint {ID = 10, X = 50.0, Y =  0.5, Z = 0.0, WeightFactor = 1.0},
				new ControlPoint {ID = 11, X = 50.0, Y =  1.0, Z = 0.0, WeightFactor = 1.0},
			};
		}

		private List<Knot> ShellElementKnot()
		{
			return new List<Knot>()
			{
				new Knot(){ID=0,Ksi=0.0,Heta=0.0,Zeta =0.0 },
				new Knot(){ID=1,Ksi=0.0,Heta=1.0,Zeta =0.0 },
				new Knot(){ID=2,Ksi=1.0,Heta=0.0,Zeta =0.0 },
				new Knot(){ID=3,Ksi=1.0,Heta=1.0,Zeta =0.0 }
			};
		}

		private double[] ShellKnotValueVectorKsi()
		{
			return new double[8]
			{
				0, 0, 0, 0, 1, 1, 1, 1
			};
		}

		private double[] ShellKnotValueVectorHeta()
		{
			return new double[6]
			{
				0, 0, 0, 1, 1, 1
			};
		}

		private NurbsKirchhoffLoveShellElement ShellElement
		{
			get
			{
                var thickness = 1;
                var degreeKsi = 3;
                var degreeHeta = 2;
                var numberOfControlPointsHeta = 3;
                var knotValueVectorKsi = ShellKnotValueVectorKsi();
                var knotValueVectorHeta = ShellKnotValueVectorHeta();
				var gauss= new GaussQuadrature();
				var gaussPoints = gauss.CalculateElementGaussPoints(degreeKsi, degreeHeta, ShellElementKnot()).ToArray();

                var parametricGaussPointKsi = new double[degreeKsi + 1];
                for (int i = 0; i < degreeKsi + 1; i++)
                {
                    parametricGaussPointKsi[i] = gaussPoints[i * (degreeHeta + 1)].Ksi;
                }

                var parametricGaussPointHeta = new double[degreeHeta + 1];
                for (int i = 0; i < degreeHeta + 1; i++)
                {
                    parametricGaussPointHeta[i] = gaussPoints[i].Heta;
                }

				var nurbs = new Nurbs2D(degreeKsi, knotValueVectorKsi, degreeHeta, knotValueVectorHeta,
                    ShellElementControlPoints().ToArray(), parametricGaussPointKsi, parametricGaussPointHeta);
				var material = new ShellElasticMaterial2Dtransformationb()
                {
                    YoungModulus = 100,
                    PoissonRatio = 0.0
                };
				var element = new NurbsKirchhoffLoveShellElement(material,nurbs, gaussPoints, thickness);
				var patch = new Patch();
				
				foreach (var controlPoint in ShellElementControlPoints())
					element.ControlPointsDictionary.Add(controlPoint.ID, controlPoint);
				foreach (var knot in ShellElementKnot())
					element.KnotsDictionary.Add(knot.ID, knot);
				
				element.Patch = patch;
				return element;
			}
		}


		private readonly Matrix _nurbsSecondDerivativesKsiExpectedValues = Matrix.CreateFromArray(new double[12, 12]
		{
			{ 4.395808555860562, 1.3958522336955388, 0.07091859196516387, 3.164894532113163, 1.0049857826886424, 0.05105997249049246, 1.5588954756112883, 0.4950142173113578, 0.025150019785057237, 0.3279814518638883, 0.10414776630446065, 0.005291400310385801},
			{ 1.11668178695643, 2.7917044673910776, 1.11668178695643, 0.8039886261509129, 2.0099715653772847, 0.8039886261509129, 0.39601137384908586, 0.9900284346227156, 0.39601137384908586, 0.08331821304356843, 0.2082955326089213, 0.08331821304356844},
			{ 0.07091859196516387, 1.3958522336955388, 4.395808555860562, 0.05105997249049246, 1.0049857826886424, 3.164894532113163, 0.02515001978505724, 0.4950142173113578, 1.558895475611288, 0.0052914003103858, 0.10414776630446065, 0.3279814518638884},
			{ -8.463635659857237, -2.687556701086618, -0.136545783619942, -4.7708935886150385, -1.514957348065927, -0.07696992519592769, 0.04710358089058723, 0.014957348065926615, 7.599329203779763E-4, 3.7398456521327854, 1.187556701086618, 0.06033579134439229},
			{ -2.150045360869292, -5.375113402173236, -2.150045360869292, -1.2119658784527403, -3.029914696131854, -1.2119658784527403, 0.011965878452741288, 0.02991469613185323, 0.011965878452741163, 0.9500453608692933, 2.375113402173236, 0.9500453608692935},
			{ -0.136545783619942, -2.687556701086618, -8.463635659857237, -0.07696992519592769, -1.514957348065927, -4.7708935886150385, 7.599329203779842E-4, 0.014957348065926615, 0.047103580890586735, 0.06033579134439228, 1.187556701086618, 3.7398456521327863},
			{ 3.7398456521327876, 1.1875567010866184, 0.06033579134439231, 0.04710358089058695, 0.014957348065926619, 7.599329203779809E-4, -4.7708935886150385, -1.5149573480659266, -0.07696992519592769, -8.463635659857237, -2.6875567010866175, -0.13654578361994202},
			{ 0.9500453608692938, 2.3751134021732367, 0.9500453608692938, 0.011965878452741217, 0.029914696131853238, 0.011965878452741236, -1.2119658784527403, -3.0299146961318533, -1.2119658784527403, -2.150045360869292, -5.375113402173235, -2.1500453608692927},
			{ 0.06033579134439231, 1.1875567010866184, 3.7398456521327876, 7.599329203779797E-4, 0.014957348065926619, 0.04710358089058703, -0.07696992519592769, -1.5149573480659266, -4.7708935886150385, -0.136545783619942, -2.6875567010866175, -8.463635659857239},
			{ 0.32798145186388816, 0.1041477663044606, 0.005291400310385798, 1.5588954756112883, 0.4950142173113579, 0.02515001978505724, 3.164894532113163, 1.0049857826886421, 0.05105997249049247, 4.395808555860564, 1.3958522336955388, 0.07091859196516392},
			{ 0.08331821304356839, 0.2082955326089212, 0.08331821304356839, 0.39601137384908586, 0.9900284346227158, 0.39601137384908586, 0.8039886261509129, 2.0099715653772843, 0.8039886261509132, 1.1166817869564307, 2.7917044673910776, 1.116681786956431},
			{ 0.005291400310385798, 0.1041477663044606, 0.32798145186388816, 0.02515001978505724, 0.4950142173113579, 1.5588954756112883, 0.05105997249049246, 1.0049857826886421, 3.164894532113164, 0.07091859196516392, 1.3958522336955388, 4.395808555860565},
			
		});
		
		private readonly Matrix _nurbsSecondDerivativesHetaExpectedValues = Matrix.CreateFromArray(new double[12, 12]
		{
			{ 1.6116641892895254, 1.6116641892895252, 1.6116641892895252, 0.6015004717568666, 0.6015004717568667, 0.6015004717568666, 0.07188019323870527, 0.07188019323870527, 0.07188019323870527, 6.694314291889707E-4, 6.694314291889707E-4, 6.694314291889707E-4},
			{ -3.2233283785790507, -3.2233283785790503, -3.2233283785790507, -1.2030009435137332, -1.2030009435137334, -1.2030009435137332, -0.14376038647741055, -0.14376038647741055, -0.14376038647741055, -0.0013388628583779411, -0.0013388628583779413, -0.0013388628583779416},
			{ 1.6116641892895252, 1.6116641892895252, 1.6116641892895258, 0.6015004717568666, 0.6015004717568667, 0.6015004717568666, 0.07188019323870527, 0.07188019323870527, 0.07188019323870527, 6.694314291889706E-4, 6.694314291889707E-4, 6.69431429188971E-4},
			{ 0.36074998763229715, 0.36074998763229704, 0.36074998763229704, 0.8888223804795408, 0.888822380479541, 0.8888223804795408, 0.4377969545248873, 0.4377969545248873, 0.4377969545248873, 0.02691639164898922, 0.02691639164898922, 0.02691639164898922},
			{ -0.7214999752645943, -0.7214999752645941, -0.7214999752645943, -1.7776447609590815, -1.777644760959082, -1.7776447609590815, -0.8755939090497746, -0.8755939090497746, -0.8755939090497746, -0.053832783297978436, -0.05383278329797844, -0.053832783297978457},
			{ 0.36074998763229704, 0.36074998763229704, 0.3607499876322972, 0.8888223804795408, 0.888822380479541, 0.8888223804795408, 0.4377969545248873, 0.4377969545248873, 0.4377969545248873, 0.026916391648989214, 0.02691639164898922, 0.02691639164898923},
			{ 0.026916391648989187, 0.026916391648989183, 0.026916391648989183, 0.4377969545248873, 0.4377969545248874, 0.4377969545248873, 0.8888223804795408, 0.8888223804795408, 0.8888223804795408, 0.36074998763229726, 0.36074998763229726, 0.36074998763229726},
			{ -0.05383278329797837, -0.053832783297978366, -0.05383278329797837, -0.8755939090497746, -0.8755939090497749, -0.8755939090497746, -1.7776447609590815, -1.7776447609590815, -1.7776447609590815, -0.7214999752645944, -0.7214999752645945, -0.7214999752645946},
			{ 0.026916391648989183, 0.026916391648989183, 0.026916391648989194, 0.4377969545248873, 0.4377969545248874, 0.4377969545248873, 0.8888223804795408, 0.8888223804795408, 0.8888223804795408, 0.3607499876322972, 0.36074998763229726, 0.36074998763229743},
			{ 6.694314291889694E-4, 6.694314291889691E-4, 6.694314291889691E-4, 0.07188019323870527, 0.07188019323870529, 0.07188019323870527, 0.6015004717568666, 0.6015004717568666, 0.6015004717568666, 1.6116641892895247, 1.6116641892895247, 1.6116641892895247},
			{ -0.0013388628583779387, -0.0013388628583779383, -0.0013388628583779387, -0.14376038647741055, -0.14376038647741057, -0.14376038647741055, -1.2030009435137332, -1.2030009435137332, -1.2030009435137332, -3.223328378579049, -3.2233283785790494, -3.2233283785790503},
			{ 6.694314291889691E-4, 6.694314291889691E-4, 6.694314291889695E-4, 0.07188019323870527, 0.07188019323870529, 0.07188019323870527, 0.6015004717568666, 0.6015004717568666, 0.6015004717568666, 1.6116641892895243, 1.6116641892895247, 1.6116641892895254},
		});

		private readonly Matrix _nurbsSecondDerivativesKsiHetaExpectedValues = Matrix.CreateFromArray(new double[12, 12]
		{
			{ 4.61017371661404, 2.5978712777504356, 0.585568838886831, 2.389781718563515, 1.3466618978750706, 0.3035420771866259, 0.579794035910726, 0.3267187671205015, 0.07364349833027709, 0.02566482566093424, 0.014462342968278067, 0.0032598602756218845},
			{ -4.024604877727209, 0.0, 4.024604877727209, -2.086239641376889, 0.0, 2.086239641376889, -0.5061505375804488, 0.0, 0.5061505375804488, -0.02240496538531236, 0.0, 0.022404965385312364},
			{ -0.585568838886831, -2.5978712777504356, -4.610173716614041, -0.3035420771866259, -1.3466618978750706, -2.389781718563515, -0.07364349833027709, -0.3267187671205015, -0.579794035910726, -0.0032598602756218837, -0.014462342968278067, -0.02566482566093425},
			{ -3.922222251164564, -2.2102048984691494, -0.49818754577373475, -0.035567465313304646, -0.02004256287064224, -0.004517660427979713, 1.774420217339484, 0.9999005678839266, 0.22538091842836905, 0.6622866397885423, 0.3732040363130084, 0.08412143283747439},
			{ 3.424034705390829, 0.0, -3.424034705390829, 0.031049804885324873, 0.0, -0.031049804885324862, -1.549039298911115, 0.0, 1.549039298911115, -0.578165206951068, 0.0, 0.5781652069510681},
			{ 0.49818754577373475, 2.2102048984691494, 3.922222251164564, 0.004517660427979713, 0.02004256287064224, 0.03556746531330472, -0.22538091842836905, -0.9999005678839266, -1.774420217339484, -0.08412143283747438, -0.3732040363130084, -0.6622866397885425},
			{ -0.6622866397885421, -0.3732040363130082, -0.08412143283747434, -1.774420217339484, -0.9999005678839268, -0.22538091842836905, 0.03556746531330479, 0.02004256287064221, 0.0045176604279797035, 3.9222222511645626, 2.210204898469149, 0.4981875457737347},
			{ 0.5781652069510678, 0.0, -0.5781652069510678, 1.549039298911115, 0.0, -1.549039298911115, -0.031049804885324942, 0.0, 0.031049804885324817, -3.4240347053908278, 0.0, 3.4240347053908287},
			{ 0.08412143283747434, 0.3732040363130082, 0.6622866397885422, 0.22538091842836905, 0.9999005678839268, 1.774420217339484, -0.004517660427979724, -0.02004256287064221, -0.03556746531330452, -0.4981875457737346, -2.210204898469149, -3.9222222511645635},
			{ -0.025664825660934216, -0.014462342968278048, -0.0032598602756218798, -0.579794035910726, -0.3267187671205016, -0.07364349833027709, -2.389781718563515, -1.3466618978750704, -0.3035420771866259, -4.61017371661404, -2.597871277750435, -0.5855688388868311},
			{ 0.022404965385312333, 0.0, -0.022404965385312333, 0.5061505375804488, 0.0, -0.5061505375804488, 2.086239641376889, 0.0, -2.086239641376889, 4.024604877727209, 0.0, -4.02460487772721},
			{ 0.0032598602756218798, 0.014462342968278048, 0.025664825660934216, 0.07364349833027709, 0.3267187671205016, 0.579794035910726, 0.3035420771866259, 1.3466618978750704, 2.389781718563515, 0.585568838886831, 2.597871277750435, 4.6101737166140415},
		});

		#endregion

		[Fact]
		public void TestShapeNurbs2DPartitionOfUnity()
		{
			var element = Element;

            var nurbs2D = element._nurbs;
			
			for (var p = 0; p < 16; p++)
			{
				var sum = 0.0;
				for (var f = 0; f < nurbs2D.NurbsValues.GetLength(0); f++)
					sum += nurbs2D.NurbsValues[f, p];
				Assert.True(Utilities.AreValuesEqual(1.0, sum, Tolerance));
			}
		}

		[Fact]
		public void TestShapeNurbs2DValues()
		{
			var element = Element;
            var nurbs2D = element._nurbs;

			for (var i = 0; i < 16; i++)
			{
				for (var j = 0; j < 16; j++)
				{
					Assert.True(Utilities.AreValuesEqual(_nurbsExpectedValues[i, j], nurbs2D.NurbsValues[i, j],
						Tolerance));
				}
			}

		}


		[Fact]
		public void TestShapeNurbs2DDerivativeValuesKsi()
		{
			var element = Element;
            var nurbs2D = element._nurbs;

			for (var i = 0; i < 16; i++)
			{
				for (var j = 0; j < 16; j++)
				{
					Assert.True(Utilities.AreValuesEqual(_nurbsDerivativeKsiExpectedValues[i, j], nurbs2D.NurbsDerivativeValuesKsi[i, j],
						Tolerance));
				}
			}

		}

		[Fact]
		public void TestShapeNurbs2DDerivativeValuesHeta()
		{
			var element = Element;
            var nurbs2D = element._nurbs;

			for (var i = 0; i < 16; i++)
			{
				for (var j = 0; j < 16; j++)
				{
					Assert.True(Utilities.AreValuesEqual(_nurbsDerivativeHetaExpectedValues[i, j], nurbs2D.NurbsDerivativeValuesHeta[i, j],
						Tolerance));
				}
			}

		}

		[Fact]
		public void TestShapeNurbs2DSecondDerivativeValuesKsi()
		{
			var element = ShellElement;
            var nurbs2D = element._nurbs;

			for (var i = 0; i < 12; i++)
			{
				for (var j = 0; j < 12; j++)
				{
					Assert.True(Utilities.AreValuesEqual(_nurbsSecondDerivativesKsiExpectedValues[i, j], nurbs2D.NurbsSecondDerivativeValueKsi[i, j],
						Tolerance));
				}
			}

		}

		[Fact]
		public void TestShapeNurbs2DSecondDerivativeValuesHeta()
		{
			var element = ShellElement;
            var nurbs2D = element._nurbs;

			for (var i = 0; i < 12; i++)
			{
				for (var j = 0; j < 12; j++)
				{
					Assert.True(Utilities.AreValuesEqual(_nurbsSecondDerivativesHetaExpectedValues[i, j], nurbs2D.NurbsSecondDerivativeValueHeta[i, j],
						Tolerance));
				}
			}

		}

		[Fact]
		public void TestShapeNurbs2DSecondDerivativeValuesKsiHeta()
		{
			var element = ShellElement;
            var nurbs2D = element._nurbs;

			for (var i = 0; i < 12; i++)
			{
				for (var j = 0; j < 12; j++)
				{
					Assert.True(Utilities.AreValuesEqual(_nurbsSecondDerivativesKsiHetaExpectedValues[i, j], nurbs2D.NurbsSecondDerivativeValueKsiHeta[i, j],
						Tolerance));
				}
			}

		}
	}
}
