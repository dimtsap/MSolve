﻿using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.IGA.Elements;
using ISAAR.MSolve.IGA.Entities;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using Xunit;

namespace ISAAR.MSolve.IGA.Tests
{
	public class NurbsElement2D
	{
		private List<ControlPoint> ElementControlPoints()
		{
			return new List<ControlPoint>
			{
				new ControlPoint {ID = 0, X = 0.0, Y = 0.0, Z = 0.0, WeightFactor = 1.0},
				new ControlPoint {ID = 1, X = 0.0, Y = 0.037037037, Z = 0.0, WeightFactor = 1.0},
				new ControlPoint {ID = 2, X = 0.0, Y = 0.111111111, Z = 0.0, WeightFactor = 1.0},
				new ControlPoint {ID = 3, X = 0.0, Y = 0.222222222, Z = 0.0, WeightFactor = 1.0},
				new ControlPoint {ID = 12, X = 0.037037037, Y = 0.0, Z = 0.0, WeightFactor = 1.0},
				new ControlPoint {ID = 13, X = 0.037037037, Y = 0.037037037, Z = 0.0, WeightFactor = 1.0},
				new ControlPoint {ID = 14, X = 0.037037037, Y = 0.111111111, Z = 0.0, WeightFactor = 1.0},
				new ControlPoint {ID = 15, X = 0.037037037, Y = 0.222222222, Z = 0.0, WeightFactor = 1.0},
				new ControlPoint {ID = 24, X = 0.111111111, Y = 0.0, Z = 0.0, WeightFactor = 1.0},
				new ControlPoint {ID = 25, X = 0.111111111, Y = 0.037037037, Z = 0.0, WeightFactor = 1.0},
				new ControlPoint {ID = 26, X = 0.111111111, Y = 0.111111111, Z = 0.0, WeightFactor = 1.0},
				new ControlPoint {ID = 27, X = 0.111111111, Y = 0.222222222, Z = 0.0, WeightFactor = 1.0},
				new ControlPoint {ID = 36, X = 0.222222222, Y = 0.0, Z = 0.0, WeightFactor = 1.0},
				new ControlPoint {ID = 37, X = 0.222222222, Y = 0.037037037, Z = 0.0, WeightFactor = 1.0},
				new ControlPoint {ID = 38, X = 0.222222222, Y = 0.111111111, Z = 0.0, WeightFactor = 1.0},
				new ControlPoint {ID = 39, X = 0.222222222, Y = 0.222222222, Z = 0.0, WeightFactor = 1.0},
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

		private Vector KnotValueVector()
		{
			return new Vector(new double[16]
			{
				0.0, 0.0, 0.0, 0.0, 0.111111111, 0.22222222, 0.33333333, 0.44444444, 0.55555555, 0.66666666, 0.77777777,
				0.88888888, 1.0, 1.0, 1.0, 1.0
			});
		}

		private NURBSElement2D Element
		{
			get
			{
				var element = new NURBSElement2D();
				var patch = new Patch();
				patch.Material= new ElasticMaterial2D(StressState2D.PlaneStrain)
				{
					YoungModulus = 1,
					PoissonRatio = 0.3
				};
				foreach (var controlPoint in ElementControlPoints())
					element.ControlPointsDictionary.Add(controlPoint.ID, controlPoint);
				foreach (var knot in ElementKnot())
					element.KnotsDictionary.Add(knot.ID, knot);
				patch.Thickness = 1;
				patch.DegreeKsi = 3;
				patch.DegreeHeta = 3;
				patch.NumberOfControlPointsHeta = 12;
				patch.KnotValueVectorKsi = KnotValueVector();
				patch.KnotValueVectorHeta = KnotValueVector();
				element.Patch = patch;
				return element;
			}
		}

		private readonly Matrix2D _expectedStiffnessMatrix = new Matrix2D(new double[32, 32]
		{
			{ 0.44505494480126356, 0.24038461538461528, 0.1419642858312804, -0.0372596152980769, 0.018612637733054365, -0.010016025713141015, 1.3736266324764615E-4, -8.012820657051282E-4, -0.18461538395973207, 0.03725961529807692, -0.1930889418144461, -0.14438100894471154, -0.047656250328230586, -0.03881209954827724, -0.0037259616019179272, -0.0031049679973958345, -0.07857142878676508, 0.010016025713141032, -0.06499828326825709, -0.03881209954827725, -0.01516998649003986, -0.010433360192975432, -0.0011160714573417883, -8.346688244524579E-4, -0.008791208903667585, 8.012820657051294E-4, -0.0064732143783847625, -0.0031049679973958345, -0.0014594780583278557, -8.346688244524578E-4, -1.0302198173518456E-4, -6.67735066773505E-5},
			{ 0.24038461538461528, 0.4450549448012636, 0.03725961529807693, -0.18461538395973207, 0.010016025713141028, -0.07857142878676508, 8.012820657051287E-4, -0.008791208903667583, -0.037259615298076885, 0.1419642858312804, -0.14438100894471156, -0.1930889418144461, -0.038812099548277246, -0.06499828326825709, -0.003104967997395834, -0.006473214378384763, -0.01001602571314101, 0.018612637733054365, -0.038812099548277246, -0.04765625032823058, -0.01043336019297543, -0.015169986490039861, -8.346688244524576E-4, -0.0014594780583278562, -8.012820657051277E-4, 1.373626632476452E-4, -0.0031049679973958345, -0.0037259616019179272, -8.346688244524578E-4, -0.0011160714573417883, -6.677350667735051E-5, -1.0302198173518459E-4},
			{ 0.14196428583128043, 0.037259615298076934, 0.5911401068520181, -0.01502403792067307, 0.29210164935090555, -0.055488781575320494, 0.0348901103581034, -0.014823718013621794, -0.1930889418144461, 0.14438100894471148, -0.31802884333354975, -0.0023287258722956626, -0.2078125002385139, -0.08461037657371795, -0.03197115424443938, -0.02018229187129407, -0.06499828326825707, 0.03881209954827724, -0.12603021950902876, -6.260015845352545E-4, -0.0773523357650367, -0.022744725101996537, -0.011366758458338528, -0.005425347328892898, -0.006473214378384764, 0.0031049679973958336, -0.013667582525064515, -5.008012730368613E-5, -0.00813873640624927, -0.0018195780278111656, -0.0011675824509985785, -4.34027790998932E-4},
			{ -0.03725961529807691, -0.1846153839597321, -0.01502403792067306, 0.4139423064802436, -0.014623396830128194, 0.08125000032082072, -0.005208333268429487, -0.007692307872101654, 0.14438100894471148, -0.19308894181444616, 0.0023287258722956773, 0.059224760241795564, -0.07374298894611378, -0.060727163508334896, -0.017077323929787665, -0.019951923384999997, 0.038812099548277246, -0.04765625032823058, 6.260015845352559E-4, -0.004356970769804685, -0.01982338431406918, -0.022385817477778597, -0.004590678519464478, -0.00552884628151786, 0.003104967997395834, -0.0037259616019179272, 5.0080127303685465E-5, -0.0018629807569222182, -0.001585870762252939, -0.0023437500408370546, -3.6725428552350465E-4, -4.8076924623626416E-4},
			{ 0.018612637733054365, 0.01001602571314103, 0.29210164935090555, -0.014623396830128194, 0.2919642900153031, -0.08179754315571579, 0.053571429996890706, -0.025774573195245738, -0.04765625032823058, 0.038812099548277246, -0.2078125002385139, 0.07374298894611378, -0.1786859000130428, -0.012678619159688842, -0.030689103427123676, -0.012937366748965017, -0.015169986490039858, 0.01043336019297543, -0.07735233576503668, 0.019823384314069183, -0.06867750445742977, -0.0034082309893607565, -0.011977259255141962, -0.0034777867935919735, -0.0014594780583278557, 8.346688244524578E-4, -0.008138736406249266, 0.0015858707622529388, -0.0073412700618291915, -2.72658482093572E-4, -0.0012896825951882482, -2.782229464921657E-4},
			{ -0.010016025713141016, -0.07857142878676508, -0.055488781575320494, 0.08125000032082072, -0.08179754315571579, 0.15185439694562838, -0.02096688075774573, 0.032967033547788474, 0.03881209954827724, -0.06499828326825707, 0.08461037657371792, -0.06072716350833491, 0.012678619159688822, -0.009136046904170188, -0.005692441403078262, 0.0020489925359304727, 0.010433360192975426, -0.015169986490039858, 0.022744725101996523, -0.022385817477778597, 0.003408230989360756, -0.010782490368232103, -0.0015302261981948913, -0.0011408730832775168, 8.346688244524577E-4, -0.0011160714573417883, 0.0018195780278111653, -0.0023437500408370546, 2.726584820935722E-4, -0.0015272054827213213, -1.2241809717770675E-4, -2.2130648241256784E-4},
			{ 1.373626632476463E-4, 8.012820657051285E-4, 0.03489011035810341, -0.005208333268429487, 0.053571429996890706, -0.02096688075774573, 0.012362637849065947, -0.0066773506677350455, -0.003725961601917927, 0.0031049679973958345, -0.03197115424443939, 0.01707732392978766, -0.03068910342712369, 0.0056924414030782595, -0.0051282053506776564, -0.0010349893510950866, -0.0011160714573417883, 8.346688244524578E-4, -0.01136675845833853, 0.004590678519464478, -0.011977259255141964, 0.0015302261981948904, -0.0021825397863492082, -2.782229464921657E-4, -1.0302198173518456E-4, 6.677350667735049E-5, -0.0011675824509985782, 3.672542855235046E-4, -0.0012896825951882482, 1.2241809717770672E-4, -2.4420025805555595E-4, -2.225783595975789E-5},
			{ -8.012820657051284E-4, -0.008791208903667585, -0.014823718013621796, -0.0076923078721016554, -0.025774573195245738, 0.032967033547788474, -0.0066773506677350455, 0.012362637761497262, 0.0031049679973958336, -0.006473214378384763, 0.020182291871294073, -0.019951923384999997, 0.01293736674896501, 0.002048992535930472, 0.0010349893510950851, 0.003943452477694681, 8.346688244524577E-4, -0.0014594780583278562, 0.005425347328892896, -0.00552884628151786, 0.0034777867935919713, -0.0011408730832775172, 2.7822294649216547E-4, 5.170177210742364E-4, 6.67735066773505E-5, -1.0302198173518456E-4, 4.3402779099893196E-4, -4.807692462362641E-4, 2.782229464921656E-4, -2.2130648241256795E-4, 2.2257835959757862E-5, 3.815628676140337E-6},

			{ -0.1846153839597321, -0.037259615298076885, -0.1930889418144461, 0.14438100894471148, -0.04765625032823058, 0.03881209954827724, -0.003725961601917927, 0.0031049679973958336, 0.41394230648024355, -0.015024037920673042, 0.05922476024179558, 0.002328725872295669, -0.004356970769804689, 6.260015845352511E-4, -0.0018629807569222185, 5.008012730368508E-5, 0.08125000032082072, -0.014623396830128192, -0.06072716350833491, -0.0737429889461138, -0.022385817477778593, -0.019823384314069183, -0.002343750040837055, -0.001585870762252939, -0.007692307872101654, -0.005208333268429488, -0.019951923385000007, -0.017077323929787662, -0.005528846281517861, -0.004590678519464479, -4.807692462362642E-4, -3.6725428552350465E-4},
			{ 0.03725961529807692, 0.14196428583128046, 0.14438100894471148, -0.19308894181444614, 0.038812099548277246, -0.06499828326825709, 0.0031049679973958345, -0.006473214378384763, -0.015024037920673042, 0.5911401068520181, -0.002328725872295669, -0.3180288433335498, -6.260015845352576E-4, -0.12603021950902876, -5.00801273036865E-5, -0.013667582525064515, -0.05548878157532047, 0.2921016493509056, -0.08461037657371795, -0.20781250023851386, -0.022744725101996533, -0.07735233576503671, -0.0018195780278111658, -0.008138736406249267, -0.014823718013621796, 0.03489011035810342, -0.020182291871294077, -0.03197115424443939, -0.005425347328892899, -0.011366758458338528, -4.3402779099893207E-4, -0.0011675824509985782},
			{ -0.1930889418144461, -0.14438100894471156, -0.31802884333354964, 0.002328725872295678, -0.20781250023851386, 0.08461037657371794, -0.03197115424443938, 0.02018229187129407, 0.059224760241795585, -0.002328725872295662, 0.5103064870541438, 9.390023362379752E-4, 0.22466947174899374, 0.003468048723607768, 0.022956731075979994, 9.264823424979948E-4, -0.060727163508334896, 0.07374298894611378, 0.060907452162476804, 9.139622689803607E-4, -0.0026141826607984243, -0.03830712516649974, -0.004927884715819457, -0.008905916329443781, -0.019951923385, 0.017077323929787662, -0.019711538705667074, 3.2552081755809125E-4, -0.016346154198960343, -0.00860543555210003, -0.002884615477860578, -0.001986511817040599},
			{ -0.14438100894471154, -0.1930889418144461, -0.0023287258722956574, 0.059224760241795564, 0.07374298894611378, -0.06072716350833491, 0.01707732392978766, -0.019951923384999997, 0.002328725872295668, -0.31802884333354975, 9.390023362379734E-4, 0.5103064870541439, 9.139622689803669E-4, 0.06090745216247679, 3.255208175580923E-4, -0.019711538705667074, 0.08461037657371794, -0.2078125002385139, 0.0034680487236077653, 0.2246694717489938, -0.038307125166499735, -0.0026141826607984195, -0.008605435552100029, -0.016346154198960343, 0.020182291871294073, -0.03197115424443939, 9.264823424979963E-4, 0.02295673107598, -0.00890591632944378, -0.004927884715819455, -0.0019865118170405995, -0.0028846154778605785},
			{ -0.047656250328230586, -0.038812099548277246, -0.2078125002385139, -0.07374298894611378, -0.17868590001304277, 0.012678619159688822, -0.030689103427123682, 0.01293736674896501, -0.0043569707698046896, -6.260015845352598E-4, 0.22466947174899377, 9.139622689803662E-4, 0.24031450620931513, 0.005112346263187768, 0.0452323728595726, 0.001610910766710072, -0.022385817477778586, 0.01982338431406918, -0.002614182660798432, 0.04257228250948184, 0.016129140057261813, 0.004976017002509567, 0.004313568435682189, -0.0033358931640680694, -0.00552884628151786, 0.004590678519464478, -0.016346154198960343, 0.010124532731954463, -0.012553419275649524, 0.0017722800796218395, -0.002029914639405146, -5.953971216390682E-4},
			{ -0.03881209954827724, -0.06499828326825709, -0.08461037657371795, -0.060727163508334896, -0.012678619159688842, -0.009136046904170195, 0.005692441403078258, 0.00204899253593047, 6.26001584535252E-4, -0.12603021950902873, 0.003468048723607768, 0.060907452162476776, 0.005112346263187761, 0.1704598791670373, 0.0013104300001836257, 0.039194139621822305, 0.022744725101996526, -0.0773523357650367, 0.05449135287626869, -0.0026141826607984225, 0.01888159938315193, 0.0624937999551716, -6.399136842837744E-5, 0.016170635143339784, 0.005425347328892897, -0.011366758458338526, 0.013229500683376737, -0.0049278847158194545, 0.0050441818500211455, 0.004397512204557419, 1.3911145181178726E-4, 0.0014804639994483145},
			{ -0.0037259616019179272, -0.0031049679973958336, -0.03197115424443938, -0.017077323929787662, -0.030689103427123676, -0.005692441403078262, -0.0051282053506776564, 0.0010349893510950851, -0.0018629807569222185, -5.0080127303686366E-5, 0.022956731075979994, 3.255208175580923E-4, 0.0452323728595726, 0.0013104300001836257, 0.011498397841682119, 4.1733440170940174E-4, -0.0023437500408370546, 0.0015858707622529384, -0.004927884715819456, 0.010124532731954463, 0.0043135684356821885, 0.006179331495086584, 0.0022569445202131047, 4.0620548212695823E-4, -4.8076924623626405E-4, 3.6725428552350454E-4, -0.0028846154778605785, 0.002420539589409723, -0.002029914639405146, 0.0016081286113336908, -2.1367523189064468E-4, 1.4467592933137458E-4},
			{ -0.0031049679973958345, -0.006473214378384763, -0.020182291871294073, -0.019951923384999997, -0.012937366748965017, 0.002048992535930474, -0.0010349893510950866, 0.00394345247769468, 5.008012730368509E-5, -0.013667582525064514, 9.264823424979961E-4, -0.019711538705667074, 0.0016109107667100703, 0.039194139621822305, 4.1733440170940185E-4, 0.01642055904328428, 0.0018195780278111651, -0.008138736406249266, 0.013229500683376739, -0.016346154198960343, 0.010853476896695826, 0.016170635143339784, 0.001541355099247685, 0.00811393494736662, 4.340277909989319E-4, -0.001167582450998578, 0.003221821651509083, -0.0028846154778605785, 0.002743278221242879, 0.0014804639994483145, 4.1176995964654586E-4, 9.691697592986609E-4},

			{ -0.07857142878676508, -0.01001602571314101, -0.06499828326825709, 0.03881209954827725, -0.015169986490039858, 0.010433360192975426, -0.0011160714573417883, 8.346688244524577E-4, 0.08125000032082072, -0.05548878157532047, -0.0607271635083349, 0.08461037657371794, -0.02238581747777859, 0.022744725101996526, -0.0023437500408370546, 0.0018195780278111647, 0.15185439694562847, -0.08179754315571582, -0.009136046904170214, 0.012678619159688825, -0.01078249036823211, 0.003408230989360752, -0.0015272054827213217, 2.7265848209357183E-4, 0.032967033547788474, -0.02096688075774573, 0.0020489925359304653, -0.005692441403078263, -0.0011408730832775176, -0.0015302261981948917, -2.213064824125679E-4, -1.224180971777068E-4},
			{ 0.010016025713141032, 0.01861263773305436, 0.03881209954827724, -0.04765625032823058, 0.01043336019297543, -0.01516998649003986, 8.346688244524578E-4, -0.0014594780583278562, -0.014623396830128185, 0.2921016493509056, 0.07374298894611377, -0.2078125002385139, 0.01982338431406918, -0.07735233576503668, 0.0015858707622529384, -0.008138736406249266, -0.08179754315571582, 0.2919642900153032, -0.01267861915968884, -0.17868590001304288, -0.0034082309893607578, -0.06867750445742979, -2.726584820935722E-4, -0.0073412700618291915, -0.02577457319524573, 0.05357142999689071, -0.012937366748965019, -0.03068910342712369, -0.003477786793591972, -0.011977259255141962, -2.7822294649216563E-4, -0.0012896825951882484},
			{ -0.06499828326825709, -0.038812099548277246, -0.12603021950902876, 6.260015845352548E-4, -0.07735233576503668, 0.022744725101996526, -0.01136675845833853, 0.005425347328892896, -0.06072716350833491, -0.08461037657371794, 0.06090745216247679, 0.0034680487236077636, -0.002614182660798425, 0.05449135287626869, -0.004927884715819456, 0.013229500683376739, -0.009136046904170221, -0.01267861915968884, 0.1704598791670373, 0.00511234626318777, 0.06249379995517158, 0.01888159938315193, 0.004397512204557419, 0.005044181850021146, 0.0020489925359304683, 0.0056924414030782595, 0.03919413962182229, 0.0013104300001836279, 0.01617063514333978, -6.3991368428377E-5, 0.0014804639994483145, 1.3911145181178723E-4},
			{ -0.03881209954827725, -0.04765625032823058, -6.260015845352558E-4, -0.004356970769804682, 0.019823384314069183, -0.022385817477778597, 0.004590678519464478, -0.00552884628151786, -0.0737429889461138, -0.20781250023851386, 9.13962268980364E-4, 0.22466947174899382, 0.04257228250948184, -0.0026141826607984225, 0.010124532731954463, -0.016346154198960343, 0.012678619159688825, -0.17868590001304288, 0.005112346263187769, 0.24031450620931516, 0.004976017002509566, 0.01612914005726181, 0.0017722800796218388, -0.012553419275649524, 0.012937366748965015, -0.030689103427123682, 0.0016109107667100695, 0.04523237285957261, -0.0033358931640680685, 0.0043135684356821885, -5.953971216390681E-4, -0.002029914639405146},
			{ -0.015169986490039861, -0.01043336019297543, -0.0773523357650367, -0.01982338431406918, -0.06867750445742977, 0.0034082309893607556, -0.011977259255141962, 0.0034777867935919713, -0.022385817477778593, -0.022744725101996537, -0.0026141826607984295, -0.038307125166499735, 0.01612914005726181, 0.018881599383151933, 0.00431356843568219, 0.010853476896695824, -0.010782490368232114, -0.003408230989360759, 0.06249379995517158, 0.004976017002509567, 0.07491844228056799, 0.027833886355854556, 0.014646291639303444, 0.0087705145351544, -0.0011408730832775168, 0.0015302261981948906, 0.01617063514333978, 0.006179331495086584, 0.01799450597868881, 0.007134563627869408, 0.003434066067719293, 0.0016711924874317438},
			{ -0.010433360192975432, -0.015169986490039863, -0.022744725101996537, -0.022385817477778593, -0.0034082309893607573, -0.010782490368232103, 0.0015302261981948906, -0.0011408730832775172, -0.019823384314069183, -0.0773523357650367, -0.038307125166499735, -0.0026141826607984247, 0.004976017002509565, 0.0624937999551716, 0.006179331495086584, 0.016170635143339784, 0.003408230989360753, -0.06867750445742977, 0.018881599383151933, 0.016129140057261807, 0.027833886355854556, 0.07491844228056799, 0.007134563627869408, 0.017994505978688808, 0.003477786793591972, -0.01197725925514196, 0.010853476896695828, 0.004313568435682191, 0.0087705145351544, 0.014646291639303444, 0.0016711924874317438, 0.003434066067719292},
			{ -0.0011160714573417883, -8.346688244524576E-4, -0.011366758458338528, -0.004590678519464478, -0.011977259255141962, -0.0015302261981948915, -0.0021825397863492082, 2.7822294649216536E-4, -0.002343750040837055, -0.0018195780278111653, -0.004927884715819457, -0.008605435552100029, 0.004313568435682189, -6.399136842837722E-5, 0.0022569445202131047, 0.001541355099247685, -0.0015272054827213217, -2.726584820935722E-4, 0.004397512204557419, 0.0017722800796218388, 0.014646291639303442, 0.007134563627869408, 0.004218177845871684, 0.002272154058345206, -2.2130648241256787E-4, 1.224180971777067E-4, 0.0014804639994483145, 0.0016081286113336906, 0.003434066067719292, 0.0024057010771085614, 9.157509661664389E-4, 5.824133753487069E-4},
			{ -8.346688244524579E-4, -0.0014594780583278562, -0.005425347328892898, -0.00552884628151786, -0.003477786793591973, -0.0011408730832775163, -2.782229464921657E-4, 5.170177210742363E-4, -0.0015858707622529392, -0.008138736406249267, -0.008905916329443781, -0.016346154198960343, -0.0033358931640680694, 0.016170635143339784, 4.0620548212695845E-4, 0.008113934947366622, 2.726584820935719E-4, -0.007341270061829191, 0.005044181850021146, -0.012553419275649524, 0.0087705145351544, 0.017994505978688808, 0.002272154058345206, 0.008110119417935973, 2.782229464921656E-4, -0.0012896825951882484, 0.002743278221242879, -0.002029914639405146, 0.003340530181408776, 0.0034340660677192924, 7.159603923091772E-4, 0.0014880953242802386},

			{ -0.008791208903667585, -8.012820657051279E-4, -0.006473214378384764, 0.003104967997395834, -0.0014594780583278557, 8.346688244524578E-4, -1.0302198173518456E-4, 6.67735066773505E-5, -0.007692307872101654, -0.014823718013621798, -0.019951923385, 0.020182291871294073, -0.005528846281517859, 0.005425347328892897, -4.8076924623626405E-4, 4.340277909989319E-4, 0.032967033547788474, -0.02577457319524573, 0.00204899253593047, 0.012937366748965012, -0.0011408730832775168, 0.0034777867935919713, -2.213064824125679E-4, 2.782229464921656E-4, 0.01236263776149726, -0.006677350667735046, 0.003943452477694682, 0.0010349893510950858, 5.170177210742364E-4, 2.782229464921653E-4, 3.815628676140336E-6, 2.2257835959757865E-5},
			{ 8.012820657051292E-4, 1.3736266324764529E-4, 0.0031049679973958336, -0.0037259616019179272, 8.346688244524578E-4, -0.0011160714573417883, 6.677350667735049E-5, -1.0302198173518458E-4, -0.005208333268429488, 0.03489011035810342, 0.017077323929787662, -0.03197115424443938, 0.004590678519464478, -0.011366758458338526, 3.6725428552350454E-4, -0.001167582450998578, -0.02096688075774573, 0.053571429996890706, 0.0056924414030782595, -0.030689103427123682, 0.0015302261981948906, -0.011977259255141962, 1.224180971777067E-4, -0.0012896825951882484, -0.006677350667735046, 0.01236263784906594, -0.0010349893510950869, -0.0051282053506776564, -2.782229464921658E-4, -0.0021825397863492082, -2.225783595975788E-5, -2.4420025805555595E-4},
			{ -0.006473214378384764, -0.0031049679973958345, -0.013667582525064514, 5.0080127303685505E-5, -0.008138736406249266, 0.0018195780278111651, -0.0011675824509985785, 4.3402779099893196E-4, -0.019951923385000007, -0.020182291871294077, -0.019711538705667078, 9.264823424979969E-4, -0.016346154198960343, 0.013229500683376737, -0.002884615477860578, 0.003221821651509082, 0.0020489925359304653, -0.012937366748965017, 0.03919413962182229, 0.0016109107667100695, 0.016170635143339784, 0.010853476896695828, 0.0014804639994483145, 0.002743278221242879, 0.003943452477694681, -0.0010349893510950864, 0.01642055904328428, 4.173344017094025E-4, 0.008113934947366622, 0.0015413550992476857, 9.691697592986609E-4, 4.117699596465459E-4},
			{ -0.0031049679973958345, -0.0037259616019179272, -5.008012730368607E-5, -0.0018629807569222182, 0.001585870762252939, -0.0023437500408370546, 3.672542855235046E-4, -4.807692462362641E-4, -0.017077323929787662, -0.03197115424443939, 3.2552081755809185E-4, 0.02295673107598, 0.010124532731954463, -0.0049278847158194545, 0.002420539589409723, -0.0028846154778605785, -0.005692441403078264, -0.03068910342712369, 0.0013104300001836279, 0.045232372859572614, 0.006179331495086584, 0.004313568435682191, 0.0016081286113336906, -0.002029914639405146, 0.0010349893510950862, -0.0051282053506776564, 4.173344017094026E-4, 0.011498397841682124, 4.062054821269584E-4, 0.002256944520213105, 1.4467592933137464E-4, -2.1367523189064463E-4},
			{ -0.0014594780583278557, -8.346688244524579E-4, -0.00813873640624927, -0.0015858707622529392, -0.0073412700618291915, 2.7265848209357226E-4, -0.0012896825951882482, 2.782229464921656E-4, -0.005528846281517861, -0.005425347328892899, -0.016346154198960343, -0.008905916329443781, -0.012553419275649527, 0.005044181850021146, -0.002029914639405146, 0.002743278221242879, -0.0011408730832775174, -0.0034777867935919726, 0.01617063514333978, -0.0033358931640680685, 0.01799450597868881, 0.008770514535154398, 0.003434066067719293, 0.003340530181408776, 5.170177210742365E-4, -2.782229464921657E-4, 0.008113934947366624, 4.062054821269585E-4, 0.008110119417935973, 0.002272154058345206, 0.0014880953242802386, 7.159603923091774E-4},
			{ -8.346688244524578E-4, -0.0011160714573417883, -0.0018195780278111656, -0.0023437500408370546, -2.7265848209357194E-4, -0.0015272054827213213, 1.224180971777067E-4, -2.2130648241256795E-4, -0.004590678519464479, -0.011366758458338528, -0.00860543555210003, -0.004927884715819455, 0.0017722800796218388, 0.004397512204557419, 0.0016081286113336908, 0.0014804639994483145, -0.0015302261981948915, -0.011977259255141962, -6.399136842837635E-5, 0.004313568435682189, 0.007134563627869408, 0.014646291639303444, 0.0024057010771085614, 0.0034340660677192924, 2.7822294649216536E-4, -0.0021825397863492082, 0.001541355099247686, 0.002256944520213105, 0.002272154058345206, 0.004218177845871684, 5.82413375348707E-4, 9.157509661664391E-4},
			{ -1.0302198173518456E-4, -6.677350667735051E-5, -0.0011675824509985782, -3.672542855235047E-4, -0.0012896825951882484, -1.2241809717770678E-4, -2.4420025805555595E-4, 2.2257835959757862E-5, -4.8076924623626427E-4, -4.34027790998932E-4, -0.002884615477860578, -0.0019865118170405995, -0.002029914639405146, 1.3911145181178726E-4, -2.1367523189064455E-4, 4.1176995964654586E-4, -2.2130648241256793E-4, -2.782229464921657E-4, 0.0014804639994483145, -5.953971216390679E-4, 0.003434066067719293, 0.0016711924874317438, 9.157509661664389E-4, 7.159603923091773E-4, 3.815628676140336E-6, -2.225783595975788E-5, 9.69169759298661E-4, 1.4467592933137467E-4, 0.0014880953242802383, 5.824133753487069E-4, 3.434066181936821E-4, 1.8548196966999093E-4},
			{ -6.67735066773505E-5, -1.0302198173518459E-4, -4.34027790998932E-4, -4.8076924623626416E-4, -2.7822294649216563E-4, -2.2130648241256784E-4, -2.225783595975789E-5, 3.815628676140342E-6, -3.672542855235047E-4, -0.0011675824509985782, -0.001986511817040599, -0.002884615477860579, -5.953971216390681E-4, 0.0014804639994483147, 1.446759293313746E-4, 9.691697592986611E-4, -1.224180971777068E-4, -0.0012896825951882484, 1.3911145181178723E-4, -0.002029914639405146, 0.0016711924874317436, 0.003434066067719292, 5.824133753487069E-4, 0.0014880953242802386, 2.225783595975786E-5, -2.4420025805555595E-4, 4.11769959646546E-4, -2.1367523189064452E-4, 7.159603923091773E-4, 9.157509661664391E-4, 1.8548196966999093E-4, 3.434066181936821E-4},
		});

		private const double Tolerance = 1e-14;

		[Fact]
		private void TestNurbsElement2DStiffnessMatrix()
		{
			var element = Element;

			var stiffnessMatrix = element.StiffnessMatrix(element);

			for (var i = 0; i < 32; i++)
			{
				for (var j = 0; j < 32; j++)
				{
					Assert.True(Utilities.AreValuesEqual(_expectedStiffnessMatrix[i, j], stiffnessMatrix[i, j],
						Tolerance));
				}
			}
		}
	}
}
