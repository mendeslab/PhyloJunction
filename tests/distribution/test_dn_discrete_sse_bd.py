import unittest
import math
import statistics

# pj imports
import phylojunction.utility.helper_functions as pjh
import phylojunction.calculation.math_utils as pjmath
import phylojunction.calculation.discrete_sse as sseobj
import phylojunction.distribution.dn_discrete_sse as distsse

__author__ = "Fabio K. Mendes"
__email__ = "fmendes@lsu.edu"


class TestBDTrees(unittest.TestCase):

    @classmethod
    def setUpClass(cls) -> None:
        # not state-dependent (just state 0, and no transition)
        total_n_states = 1

        # not state-dependent (just state 0, and no transition)
        rates_t0_s0 = [
            sseobj.DiscreteStateDependentRate(
                name="lambda",
                val=1.0,
                event=sseobj.MacroevolEvent.W_SPECIATION,
                states=[0,0,0]),
            sseobj.DiscreteStateDependentRate(
                name="mu",
                val=0.8,
                event=sseobj.MacroevolEvent.EXTINCTION,
                states=[0])
        ]

        matrix_atomic_rate_params = [rates_t0_s0]

        state_dep_par_manager = \
            sseobj.DiscreteStateDependentParameterManager(
                matrix_atomic_rate_params,
                total_n_states)

        event_handler = sseobj.MacroevolEventHandler(state_dep_par_manager)

        cls.sse_stash = sseobj.SSEStash(event_handler)


    def test_tree_size_total_count_max_taxa_bd(self):
        """
        Test if birth-death trees simulated with PJ have similar total number
        of taxa and root ages as trees simulated with geiger.
        """

        # note that the method used by both PhyloJunction and geiger, when
        # simulating for max taxa, is what Stadler (2011) calls the "simple
        # sampling approach" (SSA).
        #
        # trees simulated with TreeSim (which simulates using the GSA approach,
        # see paper above) for the parameters below will have twice the number
        # of tips, and larger root ages.
        
        stop_condition = "size"
        stop_condition_value = [ 10 ] # 10 taxa

        start_at_origin = False

        # simulation initialization
        n_batches = i = 100
        n_sim = 100
        start_states_list = [0 for i in range(n_sim)]

        # simulations
        print(("\n\nRunning TestBDTrees.test_tree_size_total_count_max_"
               "taxa_bd"))
        sim_batches = list()
        for i in range(n_batches):
            sse_sim = distsse.DnSSE(
                self.sse_stash,
                n=100,
                origin=start_at_origin,
                start_states_list=start_states_list,
                stop=stop_condition,
                stop_value=stop_condition_value,
                condition_on_speciation=True,
                condition_on_survival=True,
                epsilon=1e-12,
                runtime_limit=3600,
                max_n_failed_attempts=400,
                debug=False)

            trs = sse_sim.generate()

            sim_batches.append(trs)

            # printing progress
            pjh.print_progress(i , n_batches)

        n_total_mean_maxtaxa_geiger = [
             23.8, 24.98, 23.92, 22.52, 23.65, 22.65, 23.13, 23.67,  23.2, 24.85,
            23.25,  23.2, 24.35, 23.57, 24.16, 24.36, 23.39, 24.96, 21.43, 22.88,
            23.41, 23.24, 22.43, 24.17, 23.04, 23.08, 22.75, 23.12, 22.22,  22.8,
            22.69, 22.84, 22.89, 22.66, 21.98, 22.11, 22.69, 22.72,  23.7, 22.21,
            23.35, 23.02, 24.37, 21.75, 22.57, 22.34, 22.02, 25.04, 25.14, 22.14,
            24.42, 20.03,  24.1, 22.82, 22.25, 22.88, 22.43, 22.81, 23.47, 24.14,
            22.76, 23.29,  23.9, 21.73, 24.07, 22.59, 22.53, 22.73, 23.66, 23.37,
             24.6, 22.87, 24.63, 22.15, 23.51, 23.59, 21.92,  22.5, 24.58, 26.73,
            23.89, 24.31, 24.09, 24.68, 23.86, 24.47, 23.56, 22.92, 22.33, 21.51,
            24.14, 23.73, 23.04, 26.01, 23.55, 22.88, 24.84, 21.85, 24.47, 23.43 ]

        n_total_ci_width_maxtaxa_geiger = [
            1.9625720497615,  2.41564402889653, 2.22999494894568, 1.95652831930499, 2.15507739672119, 1.67876158058387, 2.21742649305929, 2.21907084423124, 2.07033676077838, 3.00034203437413,
            2.19910221350772, 2.07501718624987, 2.79296960819899, 2.09256919116142, 2.27803315625785, 2.45792932610383, 2.37774414064757, 2.45056140638521, 1.75976143975497, 2.36175995187975,
            2.13948553015377, 2.64845454428263, 2.09812493527921, 1.97342752068424, 1.90781469315926, 2.00727756169817,  1.8607491811503, 2.23895583744303, 1.97588881051418, 1.80348403560072,
            2.60182939183998, 2.17307027372459, 1.97339605919963, 1.79468829920773, 1.94349421694872, 1.84746364117255, 2.09783639924085, 1.88657591215494, 2.38891402570754, 1.98191168156791,
            2.29054519556574, 1.79268934973309,  2.0493612627974, 1.83449650153425, 2.00661534254308, 1.96241782183987, 1.71299329062824, 2.20905727285362, 2.72569937765956, 2.05831014498041,
            2.21156074912508, 1.65127980534964, 2.41877761931873, 1.84882524579662, 1.86095770934205, 2.08979042243402, 1.77665937846024, 2.20074615041941, 2.38744435983608,  2.4132397115609,
            2.16867305665623, 1.86052394445188, 2.27680296854509, 1.77460513564983, 2.14205220244041, 1.81337478517684, 2.04742901249643, 1.86105779457476, 1.90604433327351, 2.00823908230997,
            2.23231155278265, 2.39720892622224, 2.45905549117682, 1.94975250436634, 2.02130490865118, 1.90153469926932, 1.90063353010452, 1.58142986061362, 2.40359839197379, 2.74553462137554,
            2.27840361586082, 2.27479012738238, 2.60650822396864, 2.48217896326122,  2.4299050000967, 2.50095748646212, 2.19785697791133,  1.8975686031436, 1.84728719953768, 1.98057985550707,
            1.96376196545551,  2.0721830467876, 1.98867938881941, 2.12864837809718, 2.17157961406024,  1.9156010449027, 2.70608226083839, 1.70695588221793,  2.2142392782656, 2.14784129391312 ]

        root_ages_mean_maxtaxa_geiger = [
            4.75919719451304, 5.20987497009789, 4.42214969938176, 4.43818599404481,  4.6341960624793, 4.56861106538485, 4.55582644597555, 4.67982237467629, 4.40979759110583, 4.99341363381576,
            4.55996289562049, 4.50965287737105, 4.69654251483261, 4.64455748507909, 4.57677392768319, 4.97252209486094, 4.53562940240093, 4.75938806627352, 4.11714581057345, 4.56897685166536,
            4.64975074226505, 4.11709826436588, 4.35904065979678, 4.64585378689555, 4.27256247840422, 4.56245530719198, 4.41722586824256, 4.29226959962754, 4.28826416290394, 4.45937865210144,
            4.23259606950589, 4.34458727965486, 4.46082107400556, 4.30648950647023, 4.30799585907166, 4.33045073363683, 4.59935978781897, 4.72726692528105, 4.82620413495524, 4.13918894336496,
            4.31545169858186, 4.49386490867832, 5.00361665513474, 4.16135893393947, 4.34627889452647, 4.30483483844008, 4.09957268764572, 5.02150644071729, 5.12630050983125, 4.49644068753224,
            4.98329431323828, 3.91000977984065, 4.83814008400218, 4.33097594846902, 4.16691981144242, 4.46563459118089, 4.24987491114575, 4.40662976466601, 4.53295515954265,   4.744221501693,
            4.35788952690176, 4.51863200248317, 4.60701762531806, 4.37647196687126, 4.63421851612157, 4.46773796194996, 4.34109041336949, 4.79488683669726, 4.68675602453657, 4.49092584503989,
            5.01849250872929, 4.48414026925606, 5.28540153335028, 4.33217417686163, 4.74537177593351, 4.65408516441171, 3.94970210217709, 4.34211315377471, 4.92819389472028, 5.52793215140959,
            4.61179397728341, 4.69854098392395, 4.57292240865579, 4.98772454725259, 4.52493291357441, 4.90914823109389,  4.6153793937955, 4.56524278088623, 4.24927755153827, 4.19642432218669,
            5.23912055257292, 4.90263672601394, 4.54769510896059, 5.24802393411313, 4.68467516502108, 4.67428216360859, 4.88085386904166, 4.39089461469432, 4.99809318643099, 4.63679151115732 ]

        root_ages_ci_width_maxtaxa_geiger = [
            0.619115457862184, 0.659348117762013, 0.552695697073964,  0.57988064011635, 0.591304223615948, 0.570029632634927,  0.62872276043324, 0.706141294603482, 0.560972530473922, 0.728209581245204,
            0.690854829292233,  0.62382486111516, 0.709760133569344, 0.555341236321892, 0.588085560062609, 0.785273696063883, 0.637854418338117, 0.671656670316297, 0.504890557351337,  0.72943422695993,
            0.627312584198747, 0.643008224537711,  0.58330204553214, 0.544502154431687, 0.490865579702361, 0.555597492069546, 0.516232329047186, 0.605253640135854,  0.57341136020713, 0.541551335301634,
            0.610126477168408, 0.605760537233207, 0.596884224999702, 0.494793683397535, 0.562977643563764, 0.575670256549305,  0.61929039408224, 0.559856930114848, 0.721510760059612, 0.487405223740379,
            0.616191874773647, 0.512522423552798, 0.605159284373032, 0.583661127548458, 0.566322294130923,  0.55107634250668, 0.487998184073022, 0.678234270200986, 0.824078338786173, 0.609469603220884,
             0.61207594510121,   0.5382495910301,  0.69314756385037, 0.493458637154855, 0.484680235488658, 0.605011788711805, 0.481504638024688, 0.641078962009783, 0.582764844970514,  0.65224599068827,
            0.614509853443098, 0.525553150801019, 0.591429787878082, 0.512378027859671, 0.570389362705274, 0.557995721388504, 0.607549948111608, 0.560879395086076, 0.535574527023085, 0.555391763035823,
            0.664220703683588, 0.691811036106428, 0.659830461052841, 0.578856569421069, 0.573030721386176, 0.516653132158309, 0.465381466128567, 0.492384156090673, 0.665019850961892, 0.711665078455223,
            0.583835501040612, 0.594321905934378, 0.776095583545907, 0.820489880439986, 0.758580099035892, 0.826335190188717, 0.688525921107229, 0.523448232801541, 0.509288749598238, 0.593465097365611,
            0.633914615289665, 0.655894415698242, 0.561403510626239, 0.612560449185574,  0.58621918300181, 0.572559192666974, 0.688377916004347, 0.586286714350306, 0.673372163799723, 0.602760466055246 ]

        # parsing simulations
        n_total_ci_overlap_count = 0
        root_age_ci_overlap_count = 0
        global_mean_total = 0.0
        global_mean_root_age = 0.0
        for i, batch in enumerate(sim_batches):
            n_total = [len(ann_tr.tree.leaf_nodes()) for ann_tr in batch]
            root_ages = [ann_tr.root_age for ann_tr in batch]

            mean_total = statistics.mean(n_total)
            mean_root_ages = statistics.mean(root_ages)
            
            global_mean_total += mean_total
            global_mean_root_age += mean_root_ages
            
            stdevs_n_total = statistics.stdev(n_total)
            stdevs_root_ages = statistics.stdev(root_ages)

            sterr_n_total = stdevs_n_total / math.sqrt(n_sim)
            sterr_root_ages = stdevs_root_ages / math.sqrt(n_sim)

            n_total_ci_width_maxtaxa = 1.96 * sterr_n_total
            root_ages_ci_width_maxtaxa = 1.96 * sterr_root_ages

            if abs(mean_total - n_total_mean_maxtaxa_geiger[i]) \
                <= (n_total_ci_width_maxtaxa
                    + n_total_ci_width_maxtaxa_geiger[i]):
                n_total_ci_overlap_count += 1

            if abs(mean_root_ages - root_ages_mean_maxtaxa_geiger[i]) \
                <= (root_ages_ci_width_maxtaxa 
                    + root_ages_ci_width_maxtaxa_geiger[i]):
                root_age_ci_overlap_count += 1

        # [==== * ====][.... + ....] if we take '+' to be the "truth" of the '*' interval, + cannot be more than '====' away from '*' 95% of the time
        # then abs('+' - '*') can be at most ('====' + '....'). '....' can be added because we still are guaranteed to see '+' falling within that range
        # 95% of the time

        print("\n\nPJ global mean total taxon count = " \
              + str(global_mean_total / 100.0))
        print("geiger global mean total taxon count = " \
              + str(statistics.mean(n_total_mean_maxtaxa_geiger)))
        print("\nPJ global mean root age = " \
              + str(global_mean_root_age / 100.0))
        print("geiger global mean total root age = " \
              + str(statistics.mean(root_ages_mean_maxtaxa_geiger)))

        print("\n95% CIs of simulations here and from geiger overlapped " \
              + str(n_total_ci_overlap_count) \
              + " times for the total taxon count.")
        print("\n95% CIs of simulations here and from geiger overlapped " \
              + str(root_age_ci_overlap_count) + " times for root age.")
        exp_count = int(0.95 * n_batches)
        a_delta = math.ceil(0.07 * exp_count)
        self.assertAlmostEqual(
            n_total_ci_overlap_count,
            exp_count,
            msg="Mean absolute difference must be 1.96 * (stderr_python + " \
                + "stderr_geiger) apart " + str(exp_count) + " (+/- " \
                + str(a_delta) + ") out of 100 times.",
            delta=a_delta)
        self.assertAlmostEqual(
            root_age_ci_overlap_count,
            exp_count,
            msg="Mean absolute difference must be 1.96 * (stderr_python + " \
                + "stderr_geiger) apart " + str(exp_count) + " (+/- " \
                + str(a_delta) + ") out of 100 times.",
            delta=a_delta)


    def test_tree_size_total_count_max_t_from_root_bd(self):
        """
        Test if birth-death trees simulated with PJ have similar total
        number of taxa and number of extant taxa as trees simulated
        with phytools.
        """

        stop_condition = "age"
        stop_condition_value = [3.0]  # 3.0 time units

        start_at_origin = False

        # simulation initialization
        n_batches = i = 100
        n_sim = 100
        start_states_list = [0 for i in range(n_sim)]

        # simulations
        print(("\n\nRunning TestBDTrees.test_tree_size_total_count_max_t_from_"
              "root_bd"))
        sim_batches = list()
        for i in range(n_batches):
            sse_sim = distsse.DnSSE(
                self.sse_stash,
                n=100,
                origin=start_at_origin,
                stop=stop_condition,
                stop_value=stop_condition_value,
                condition_on_speciation=True,
                condition_on_survival=True,
                start_states_list=start_states_list, 
                epsilon=1e-12,
                runtime_limit=3600,
                debug=False)

            trs = sse_sim.generate()

            sim_batches.append(trs)

            # printing progress
            pjh.print_progress(i , n_batches)

        # phytools from root
        n_total_mean_maxt_phytools = [
            15.56, 15.12, 15.54, 15.09, 14.12, 15.41, 16.01, 14.33, 15.28, 14.87,
            16.61, 15.49, 15.02, 15.18, 15.12, 14.51, 14.87, 14.51, 15.05, 16.22,
            17.31, 15.43, 13.03, 14.64, 15.57, 13.88, 13.67, 15.28, 13.53, 15.15,
            13.54, 14.49, 14.23, 13.65, 13.74, 14.56, 15.88, 16.32, 14.09, 16.82,
            15.94, 15.28, 15.78, 15.41, 15.31, 18.09, 16.55, 14.44, 13.84, 15.74,
            13.78, 13.42, 13.71, 14.32,  15.0, 15.88, 15.71, 14.16, 14.29, 15.94,
            13.91, 16.13, 15.82, 16.38, 15.08, 15.42, 14.44, 15.78,  14.0, 15.43,
             16.0, 15.48, 14.08, 14.66, 14.13, 15.29, 15.95, 15.94, 17.01, 14.77,
            15.45, 14.43, 16.11, 15.62, 14.77, 14.56, 14.66, 14.28, 14.33, 14.94,
             15.2, 12.27, 13.76, 14.63, 15.01, 14.07, 14.86, 17.61, 17.02, 14.42 ]

        n_total_ci_width_maxt_phytools = [
             2.2328885903622, 1.93394692708142, 2.04058802537233, 1.98890864708881, 1.83299512296353, 2.04014585297575, 2.46108081668657, 1.81464970535123, 2.22291893498507, 1.67280592002199,
             2.2664503921776, 2.22822195691833, 1.86691279589453, 1.95967924648166, 1.88702836494341, 1.80208386625321, 1.87782945020052, 1.69099552245702, 1.67181279460411,  2.0691331190241,
            2.63163607413082, 2.18456375267351, 1.67646703601004, 1.81282368239008, 2.16475721959715, 1.61704469482155, 1.98010174659503, 1.89355628904925, 1.59923087119838, 1.83639923257829,
             1.8234736454793,  2.0016276605086, 2.07237029997103, 2.11142113087825, 1.80297187841171, 1.89647220135309, 2.14332442846541, 2.49838443993253, 1.97579356014823, 2.65566210458661,
            2.21360216111133, 2.29438841691566, 2.21723836455455, 2.08417801794977, 1.74368801722503, 2.37875574605589, 2.45174160320554, 1.85111996195151,  1.9739514773032, 1.91931255538344,
            1.75010428133742, 1.26139024536862, 1.84712754734523, 1.92824013593177, 2.11887410545022, 1.99401085313602, 2.34279158667301, 1.69081767034279, 1.72007237761741, 1.84261060476924,
            1.49842434091767, 2.24213820423021, 1.97486733005362,  2.0597725699698, 2.10650936735337, 1.92689135147886, 1.67388770920582, 2.12410166052438, 2.00482480657301, 2.08737048871821,
            2.15140557296711, 1.98175210550582, 1.44880993922598, 1.92407801097649,  1.9380307952229, 1.85760173787754, 2.08571905097402, 2.30143690345968, 2.32018978638569,  2.6700123166969,
            1.91946519277184, 1.80246280495114, 2.24027522804091, 2.07160713720268, 1.84766526812361, 1.73132100624542, 1.69780295057757, 1.78837666518149, 2.01968398808258, 2.05355385041777,
            2.01678940773486, 1.37436099884692, 1.51260680678066, 1.89551744743628, 1.78673012053762,  1.5608701282207, 1.88450558486674, 2.34520856180934, 2.05150449685251, 1.51697702372351 ]

        n_extant_mean_maxt_phytools = [
            6.66, 6.55, 6.63, 6.03, 5.65, 5.91, 6.69, 5.54, 6.23, 6.01,
            6.98, 6.16, 5.98, 6.38, 6.25, 5.87, 6.21,  5.7, 6.25, 6.92,
            7.05, 6.15, 5.62, 6.13, 6.59, 5.69, 5.35, 6.17, 5.11, 6.61,
            5.43, 6.09,  5.6, 5.44, 5.29, 5.99, 6.56, 6.84, 5.47, 7.11,
            6.66, 6.61,  6.5, 6.76, 6.04, 7.59,  7.1, 5.79, 5.46, 6.65,
            5.74, 5.41, 5.65, 5.71, 6.28, 6.73,  6.5, 5.63, 5.48, 6.68,
            5.61, 6.34, 5.94,  6.9, 6.03, 6.23, 6.09, 6.44, 5.73, 6.55,
            6.27, 6.54, 5.32, 5.56, 5.75, 6.14, 6.79, 6.73, 7.07, 6.17,
            6.52,  5.8, 6.62, 6.36, 5.96, 6.46, 6.22, 5.98, 5.66,  6.1,
            6.46, 4.65, 5.74, 5.74, 5.85, 6.15, 6.61, 7.19, 6.68, 5.67 ]

        n_extant_ci_width_maxt_phytools = [
             1.29419979972806, 0.985872305128916,  1.15826420530591, 0.987083852015163, 0.991758747953056,  1.01179262377987,  1.39449697038181, 0.955102561735829,  1.19947542136746, 0.966040005027292,
             1.19562758703063,  1.15156623708054, 0.933969649756707,  1.02045563679499,  1.00651700673584,  1.03775354190548,  1.06965978467431, 0.894721942021852, 0.879356307279089,  1.11804196520506,
             1.39607662375456,  1.17744271772998, 0.998934022553112,  1.05591684015993,  1.22293967342976, 0.901718349720257,  1.03315155781431,  1.09962565604589, 0.877819317365136,  1.04677820598504,
             1.00685621414361,  1.13540475157496, 0.984937058954028,  1.24813562070665, 0.941548313506358,  1.03244520895539,  1.27123899731043,  1.21677642941993,  0.96923211592758,  1.45421412231327,
             1.19602996077469,  1.22509542254611,  1.09341322694052,  1.25210876605038, 0.951235058099509,  1.28453163668964,  1.40470119760234,  1.06347481456409,  1.08411821190022,  1.09052594113246,
            0.998871867901408, 0.707563317936203, 0.980346403424356,  1.07183419417276,  1.14289373367109,  1.11183820565564,  1.22686948257454, 0.908560527383394,  0.89644640394812, 0.979358335386445,
            0.830574710606451,  1.24996856138242, 0.942999953083324,  1.09021454658051,  1.03954682511956, 0.997019048681791, 0.965879319061082,  1.07834011163623,  1.10588915439938,  1.27430454790371,
             1.10131819712249,  1.03091249880968, 0.745372879094007,  1.03874956995943,  1.06277401457723, 0.983367789021406,    1.096883857281,  1.33250154115172,  1.26632572827887,  1.31497175819425,
             1.07731765365127,  1.01481980288287,  1.19296331219465,   1.0298203463996,  1.08787000411749,   1.0022846499852, 0.985165537851227, 0.931889962724281,  1.08440451954988,  1.14676454832408,
             1.03803955394559, 0.844235282442763, 0.879464413691364,    1.087406198625, 0.939931375100858, 0.924110339294999,  1.17629528279392,  1.28923558214262,   1.0254249681242, 0.807914998029127 ]

        # parsing simulations
        n_total_ci_overlap_count = 0
        n_extant_ci_overlap_count = 0
        # root_age_ci_overlap_count = 0
        global_mean_total = 0.0
        global_mean_extant = 0.0
        for i, batch in enumerate(sim_batches):
            n_total = [len(ann_tr.tree.leaf_nodes()) for ann_tr in batch]
            n_extant = [ann_tr.n_extant_terminal_nodes for ann_tr in batch]

            mean_total = statistics.mean(n_total)
            mean_extant = statistics.mean(n_extant)

            global_mean_total += mean_total
            global_mean_extant += mean_extant

            stdevs_n_total = statistics.stdev(n_total)
            stdevs_n_extant = statistics.stdev(n_extant)

            sterr_n_total = stdevs_n_total / math.sqrt(n_sim)
            sterr_n_extant = stdevs_n_extant / math.sqrt(n_sim)

            n_total_ci_width_maxt = 1.96 * sterr_n_total
            n_extant_ci_width_maxt = 1.96 * sterr_n_extant

            if abs(mean_total - n_total_mean_maxt_phytools[i]) \
                <= (n_total_ci_width_maxt \
                    + n_total_ci_width_maxt_phytools[i]):
                n_total_ci_overlap_count += 1

            if abs(mean_extant - n_extant_mean_maxt_phytools[i]) \
                <= (n_extant_ci_width_maxt \
                    + n_extant_ci_width_maxt_phytools[i]):
                n_extant_ci_overlap_count += 1

        # [==== * ====][.... + ....] if we take '+' to be the "truth" of the '*' interval, + cannot be more than '====' away from '*' 95% of the time
        # then abs('+' - '*') can be at most ('====' + '....'). '....' can be added because we still are guaranteed to see '+' falling within that range
        # 95% of the time

        print("\n\nPJ global mean total taxon count = " \
              + str(global_mean_total / 100.0))
        print("phytools global mean total taxon count = " \
              + str(statistics.mean(n_total_mean_maxt_phytools)))
        print("\nPJ global mean extant taxon count = " \
              + str(global_mean_extant / 100.0))
        print("phytools global mean extant taxon count = " \
              + str(statistics.mean(n_extant_mean_maxt_phytools)))

        print("\n95% CIs of simulations here and from phytools overlapped " \
              + str(n_total_ci_overlap_count) \
              + " times for the sampled ancestor count.")
        print("\n95% CIs of simulations here and from phytools overlapped " \
              + str(n_extant_ci_overlap_count) \
              + " times for the sampled ancestor count.")

        exp_count = int(0.95 * n_batches)
        a_delta = math.ceil(0.07 * exp_count)
        self.assertAlmostEqual(
            n_total_ci_overlap_count,
            exp_count,
            msg="Mean absolute difference must be 1.96 * (stderr_python + " \
                + "stderr_treesim) apart " + str(exp_count) + " (+/- " \
                + str(a_delta) + ") out of 100 times.",
            delta=a_delta)
        self.assertAlmostEqual(
            n_extant_ci_overlap_count,
            exp_count,
            msg="Mean absolute difference must be 1.96 * (stderr_python + " \
                + "stderr_treesim) apart " + str(exp_count) + " (+/- " \
                + str(a_delta) + ") out of 100 times.",
            delta=a_delta)


    def test_tree_size_total_count_max_t_from_origin_bd(self):
        """
        Test if birth-death trees simulated with PJ have similar total
        number of taxa and number of extant taxa as trees simulated
        with TreeSim, starting from the origin.
        """

        stop_condition = "age"
        stop_condition_value = [ 3.0 ] # 3.0 time units

        start_at_origin = True

        # simulation initialization
        n_batches = i = 100
        n_sim = 100
        start_states_list = [0 for i in range(n_sim)]

        # simulations
        print(("\n\nRunning TestBDTrees.test_tree_size_total_count_max_t_"
               "from_origin_bd"))
        sim_batches = list()
        for i in range(n_batches):
            sse_sim = distsse.DnSSE(
                self.sse_stash,
                n=100,
                origin=start_at_origin,
                start_states_list=start_states_list,
                stop=stop_condition,
                stop_value=stop_condition_value,
                condition_on_speciation=True,
                condition_on_survival=True,
                epsilon=1e-12,
                runtime_limit=3600,
                max_n_failed_attempts=500,
                debug=False)

            trs = sse_sim.generate()

            sim_batches.append(trs)

            # printing progress
            pjh.print_progress(i , n_batches)

        # treesim from origin
        n_total_mean_maxt_treesim = [
            12.86, 11.36, 11.82, 12.37, 10.83, 11.74, 12.82, 11.53, 11.21, 10.57,
            12.51, 10.12,  9.89, 12.38, 10.78, 12.45, 10.51, 11.57, 10.66, 12.68,
            12.31, 11.16, 11.01, 12.15, 11.18, 10.53, 12.43, 12.59, 11.27, 12.36,
            10.71,    12, 12.07, 12.22,  11.3,  9.85, 11.14, 11.18, 12.53, 11.29,
            11.5,  11.02,  11.6,  10.7, 11.61, 10.94, 10.93, 12.57, 11.72, 11.35,
            12.33, 12.31, 11.45, 13.19, 12.18,  12.7, 10.52, 13.22, 11.35, 12.27,
            10.53, 11.73, 10.48, 11.39, 11.83, 11.68, 11.78, 12.35,  11.3, 10.47,
            11.54,  10.5,  10.9,  11.3, 13.14, 11.88, 12.34,    12, 11.34, 12.46,
            10.53, 13.04, 10.54, 11.33,  10.7, 13.26, 12.52,  9.68, 13.07, 10.53,
            10.63,  9.87, 11.21, 11.07, 11.19, 10.61, 10.45, 12.46,  11.6,   12.7
        ]

        n_total_ci_width_maxt_treesim = [
            2.00598965887543, 1.53334677411764, 1.84672520462957, 1.81520559806109, 1.38622464458418, 1.74856040019535,  1.8082929395293, 1.58998374478699, 1.40740440743256, 1.46520900070863,
            1.60626173647404, 1.56311217974419,  1.4053626464656, 1.72940587964933,  1.6328771241665, 1.79428933626164, 1.43692696734905, 1.55613946091544, 1.56483407659494, 1.68602813449332,
             1.7642628535923, 1.35498704395759, 1.59875522235585, 1.82325125293064, 1.52849538534806, 1.45276179531931, 1.82641540010409, 1.99630876298525, 1.49135343371128,  1.8181671330316,
            1.51136850636111, 1.91918720169543, 1.65451954074501, 1.75254154999088, 1.52165603567596, 1.26421573085507, 1.48644461121907, 1.71376331168713, 1.84253373681817, 1.26756310161098,
            1.78586227935081,  1.6572671656018, 1.34038318055214, 1.51501121774691, 1.41279805701275, 1.34264515886598, 1.38488034767403, 1.81148214376287, 1.78468420815232, 1.55940764459464,
            2.02658883861727, 1.73565808623796, 1.55666801635137,   1.693105371914, 1.57451484963374, 1.86119227563105, 1.57108544913536, 2.25110677730495, 1.34049173850871, 1.66191427704392,
            1.30293585568914,  1.3732311647786, 1.32743208983786, 1.65724726317865, 1.62202238697201, 1.45794108373293, 1.65201399754385, 1.80011899718891, 1.61493639365628, 1.40662116313016,
            1.24505394568941, 1.17451421292534, 1.59608435006589, 1.53359476032259, 2.04032178250432, 1.92469504067423,  2.0369336606227, 1.66567970777364, 1.66129891690602,  2.0020008109061,
            1.55718642468143, 1.71998777726453, 1.43393305790244, 1.42378571251137, 1.77299614019075, 1.92677857449981, 1.72754253952213, 1.42618886136543,  1.8242895618998, 1.47581701583984,
            1.44262988302403, 1.08777547519902, 1.79556916341478, 1.47681582173036, 1.32544718674434, 1.56276831761636, 1.26175780864505, 1.63734818076379, 1.96810043268409, 1.74919941283394
        ]

        root_ages_mean_maxt_treesim = [
            2.41134926469317, 2.38920913840397, 2.43975399368484, 2.36167155882785,  2.3587007364241, 2.33778283109888,  2.4577397815274, 2.42626338160149, 2.39713681308128,  2.4084376193809,
             2.4974528389355, 2.40511012325899, 2.30647267157316, 2.42565744285912, 2.26765674607618,  2.4834094327036, 2.24709301128232, 2.45818899585454, 2.40231025551501, 2.42716561100311,
            2.51587695830744, 2.43829360327129, 2.36654427609215, 2.29366342006442, 2.44129237874595, 2.44364909306366, 2.45982058094011,  2.3500470152416,   2.418030790723, 2.36690045733178,
            2.33344159000591,  2.3488046990303, 2.43983650212711, 2.45543345103501, 2.42322554561376, 2.37134138005089, 2.39025548763228, 2.38353641595131, 2.41053305434273, 2.45320698395875,
            2.56030978698553, 2.38549730708735, 2.40636439184482, 2.33156546281527, 2.45372703782439, 2.43273393188503, 2.43442793638736, 2.45165518183406, 2.35961928551049, 2.40281359050419,
             2.3125562566843, 2.39925423527952, 2.39144320514972,  2.4277004420588, 2.44469274408377, 2.41753460371012, 2.31053533778282, 2.34655967858245, 2.36056707080135,  2.4257683669853,
            2.36392051223435,  2.4433124600778,  2.3708687917585, 2.37442558977751, 2.40860605171841,  2.4345556612665, 2.35992545428336, 2.44112488756847, 2.38450054327198, 2.44956570399421,
            2.46092076374808, 2.43811520578555, 2.29448856565299, 2.37649818075683, 2.38469577201613,  2.4643347062554, 2.44166991306082, 2.45571811529939, 2.48262532651773,  2.4040188589081,
            2.49845301372491, 2.52427137935501, 2.30222139001482, 2.31649414607209, 2.40718261057246, 2.47965263780969, 2.52176339091588, 2.34595828935434, 2.41761244731798, 2.32563882473638,
            2.35136820012941, 2.43721756305005, 2.37140250618802, 2.34120221673433, 2.42689883102106, 2.30395083398858, 2.32984965307698, 2.48754373987955, 2.35823136681131, 2.43976798550804
        ]

        root_ages_ci_width_maxt_treesim = [
             0.118614537742372, 0.0973025403085302, 0.0924813356752944,  0.110194862975739,  0.114440741277904,  0.110178798503782, 0.0957698644204067,  0.116719073615701,  0.110946947854556,  0.106754991316681,
            0.0871260862329165,  0.117324236823729,   0.12459877289986,  0.119319780757935,  0.135223244008495, 0.0845561395512595,  0.125747814556681,  0.109656611487338,  0.122301468085595, 0.0910437169307571,
             0.110916473580753,  0.117261807568291,  0.105129241091335,  0.130308337844548,  0.117198030247258,  0.111154148813296, 0.0899255607893612,   0.12119151921046,  0.116914638975389,  0.121049494696862,
             0.120881228390186,  0.127553739150464,  0.112069672234155, 0.0981225418422778,   0.10700017134518,  0.124603154132274,   0.13080762953373,  0.111977296852409,  0.105382701132416,  0.101228706099052,
            0.0880208207778674,  0.122068533115286, 0.0958017554961159,  0.122574672405576, 0.0934483208090163,  0.104521531144681,  0.105018649124541,  0.110618124060909,  0.129357799956822, 0.0991629983683035,
             0.113002372302917,  0.120914396167643, 0.0998846417472681,  0.112128717510314, 0.0962641661378887,  0.107585124592922,    0.1271274597584,  0.120671574088347,  0.117222788974941,  0.113849859692386,
            0.0988910493234009,  0.109599399979178,  0.109263277034012,   0.12104537332385,  0.118114588466149, 0.0992280983891192,  0.111306108533451,  0.101338565822798, 0.0992505816563806,  0.107205384355046,
            0.0846739513590895,  0.101923577275351,  0.120585421036524,   0.13533080458462,  0.119862388629098,  0.107469642266994,   0.10187092799781,  0.102821551426194, 0.0848564790648721,  0.109156482672756,
             0.110050492222876, 0.0881731057795677,  0.113992856348529,  0.121231280141388,  0.104053085645244,   0.10507398152035, 0.0917024098969028,  0.119101188258717,  0.104633851513134,  0.123410247752883,
             0.123102428889098,  0.091508461525494,     0.115638713124,  0.117966708499831,  0.103408689657657,  0.111697405092111,  0.118632053538949, 0.0972386360773996,  0.113507182558513, 0.0976551271056787
        ]

        # parsing simulations
        n_total_ci_overlap_count = 0
        root_age_ci_overlap_count = 0
        global_mean_total = 0.0
        global_mean_root_age = 0.0
        for i, batch in enumerate(sim_batches):
            n_total = [len(ann_tr.tree.leaf_nodes()) for ann_tr in batch]
            root_ages = [ann_tr.root_age for ann_tr in batch]

            mean_total = statistics.mean(n_total)
            mean_root_ages = statistics.mean(root_ages)

            global_mean_total += mean_total
            global_mean_root_age += mean_root_ages

            stdevs_n_total = statistics.stdev(n_total)
            stdevs_root_ages = statistics.stdev(root_ages)

            sterr_n_total = stdevs_n_total / math.sqrt(n_sim)
            sterr_root_ages = stdevs_root_ages / math.sqrt(n_sim)

            n_total_ci_width_maxt = 1.96 * sterr_n_total
            root_ages_ci_width_maxt = 1.96 * sterr_root_ages

            if abs(mean_total - n_total_mean_maxt_treesim[i]) \
                <= (n_total_ci_width_maxt + n_total_ci_width_maxt_treesim[i]):
                n_total_ci_overlap_count += 1

            if abs(mean_root_ages - root_ages_mean_maxt_treesim[i]) \
                <= (root_ages_ci_width_maxt \
                    + root_ages_ci_width_maxt_treesim[i]):
                root_age_ci_overlap_count += 1

        print("\n\nPJ global mean total taxon count = " \
              + str(global_mean_total / 100.0))
        print("TreeSim global mean total taxon count = " \
              + str(statistics.mean(n_total_mean_maxt_treesim)))
        print("\nPJ global mean root age = " \
              + str(global_mean_root_age / 100.0))
        print("TreeSim global mean root age = " \
              + str(statistics.mean(root_ages_mean_maxt_treesim)))

        # [==== * ====][.... + ....] if we take '+' to be the "truth" of the '*' interval, + cannot be more than '====' away from '*' 95% of the time
        # then abs('+' - '*') can be at most ('====' + '....'). '....' can be added because we still are guaranteed to see '+' falling within that range
        # 95% of the time

        print("\n95% CIs of simulations here and from TreeSim overlapped " \
              + str(n_total_ci_overlap_count) \
              + " times for the sampled ancestor count.")
        print("\n95% CIs of simulations here and from TreeSim overlapped " \
              + str(root_age_ci_overlap_count) \
              + " times for the sampled ancestor count.")

        exp_count = int(0.95 * n_batches)
        a_delta = math.ceil(0.07 * exp_count)
        self.assertAlmostEqual(
            n_total_ci_overlap_count,
            exp_count,
            msg="Mean absolute difference must be 1.96 * (stderr_python + " \
                + "stderr_treesim) apart " + str(exp_count) + " (+/- " \
                + str(a_delta) + ") out of 100 times.",
            delta=a_delta)
        self.assertAlmostEqual(
            root_age_ci_overlap_count,
            exp_count,
            msg="Mean absolute difference must be 1.96 * (stderr_python + " \
                + "stderr_treesim) apart " + str(exp_count) + " (+/- " \
                + str(a_delta) + ") out of 100 times.",
            delta=a_delta)


    def test_expected_size_bd(self):
        """
        Test if birth-death trees simulated with PJ have the expected
        number of extant observable nodes (against the theoretical
        expectation).
        """
        
        # setting up stopping conditions
        stop_condition = "age"
        stop_condition_value = [5.0] ## 5 time units
        start_at_origin = True

        # simulation initialization
        n_batches = i = 100
        n_sim = 100
        start_states_list = [0 for i in range(n_sim)]
        # seeds_list = [i+1 for i in range(n_sim)]

        # simulations
        print("\n\nRunning TestBDTrees.test_expected_size_bd")
        sim_batches = list()
        for i in range(n_batches):
            # print("Doing batch " + str(n_batches - i))
            sse_sim = distsse.DnSSE(
                self.sse_stash,
                n=n_sim,
                origin=start_at_origin,
                start_states_list=start_states_list,
                stop=stop_condition,
                stop_value=stop_condition_value,
                condition_on_speciation=False,
                condition_on_survival=False,
                epsilon=1e-12,
                runtime_limit=3000,
                debug=False)

            trs = sse_sim.generate()

            sim_batches.append(trs)

            # printing progress
            pjh.print_progress(i , n_batches)

        # parsing results
        exp_size = pjmath.exp_extant_count_bd(1.0, 0.8, 5.0)

        # origin_heights_when_gt_0 = set()
        batch_cis = list()
        in_ci_count = 0
        global_avg_size = list()
        for i, batch in enumerate(sim_batches):
            n_extant_obs_nodes = \
                [ann_tr.n_extant_terminal_nodes for ann_tr in batch]
            stdevs = statistics.stdev(n_extant_obs_nodes)
            sterr = stdevs / math.sqrt(n_sim)
            diff = 1.96 * sterr
            avg_size = statistics.mean(n_extant_obs_nodes)
            global_avg_size.append(avg_size)
            batch_cis = (avg_size - diff, avg_size + diff)

            if exp_size >= batch_cis[0] and exp_size <= batch_cis[1]:
                in_ci_count += 1

        print("\n\nPJ global mean extant taxon count = " \
              + str(statistics.mean(global_avg_size)))
        print("Expected BD tree size =", exp_size)

        print("\n\n95% CI includes expectation " + str(in_ci_count) + " times.")
        exp_count = int(0.95 * n_batches)
        a_delta = math.ceil(0.07 * exp_count)
        self.assertAlmostEqual(
            in_ci_count,
            exp_count,
            msg="Truth should be contained within 95%-CI " \
                + str(exp_count) + " (+/- " + str(a_delta) \
                    + ") out of 100 times.",
            delta=a_delta)


if __name__ == '__main__':
    # Assuming you opened the PhyloJunction/ (repo root) folder
    # on VSCode and that you want to run this as a standalone script,
    # i.e., "Run Without Debugging", you will need to configure your
    # launch.json to have:
    # 
    # "env": {"PYTHONPATH": "${workspaceRoot}/src/phylojunction/"}
    # 
    # and your settings.json to have:
    #   
    # "python.analysis.extraPaths": [ "${workspaceFolder}/src/phylojunction/" ]
    # 
    # If you want to run this as a standalone from PhyloJunction/
    # on the terminal, remember to add "src/" to
    # PYTHONPATH (system variable), or to set it if it does not
    # exist -- don't forget to export it!
    # 
    # Then you can do:
    # $ python3 tests/distribution/test_dn_discrete_sse_bd.py
    # 
    # or
    #
    # $ python3 -m tests.distribution.test_dn_discrete_sse_bd
    #
    # or 
    #
    # $ python3 -m unittest tests.distribution.test_dn_discrete_sse_bd.TestBDTrees.test_expected_size_bd

    #############
    # Debugging #
    #############
    # total_n_states = 1

    # # not state-dependent (just state 0, and no transition)
    # rates_t0_s0 = [ sseobj.DiscreteStateDependentRate(name="lambda", val=1.0, event=sseobj.MacroevolEvent.W_SPECIATION, states=[0,0,0]),
    #                 sseobj.DiscreteStateDependentRate(name="mu", val=0.8, event=sseobj.MacroevolEvent.EXTINCTION, states=[0]) ]

    # matrix_atomic_rate_params = [ rates_t0_s0 ] # 1D: time slices (i) , 2D: all rates from all states in i-th time slice

    # state_dep_par_manager = sseobj.DiscreteStateDependentParameterManager(matrix_atomic_rate_params, total_n_states)

    # event_handler = sseobj.MacroevolEventHandler(state_dep_par_manager)

    # # stop_condition = "size"
    # stop_condition = "age"
        
    # # stop_condition_value = [ 10 ] # 10 living taxa
    # stop_condition_value = [ 3.0 ] # origin age 

    # start_at_origin = False

    # # simulation initialization
    # n_sim = 100
    # start_states_list = [0 for i in range(n_sim)]

    # sse_sim = distsse.DnSSE(
    #     event_handler,
    #     stop_condition_value,
    #     n=n_sim,
    #     stop=stop_condition,
    #     origin=start_at_origin,
    #     start_states_list=start_states_list,
    #     epsilon=1e-12,
    #     runtime_limit=3600,
    #     condition_on_speciation=True,
    #     condition_on_survival=True,
    #     debug=False)

    # trs = sse_sim.generate()

    # print(trs[0].tree.as_string(schema="newick"))

    # n_leaves = float()
    # root_age = float()
    # origin_age = float()
    # for tr in trs:
        # print(tr.tree.as_string(schema="newick"))
        # print("root age = " + str(tr.root_age) + "\n")
        # print("total taxa = " + str(tr.n_extant_terminal_nodes + tr.n_extinct_terminal_nodes) + "\n")
        # n_leaves += tr.n_extant_terminal_nodes + tr.n_extinct_terminal_nodes
        # n_leaves += len(tr.tree.leaf_nodes())
        # root_age += tr.root_age
        # origin_age += tr.origin_age

    # print("mean extant and fossil leaves = " + str(n_leaves / n_sim))
    # print("mean root age = " + str(root_age / n_sim))

    unittest.main()