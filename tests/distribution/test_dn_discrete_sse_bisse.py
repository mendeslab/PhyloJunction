import unittest
import math
import statistics

# pj imports
import phylojunction.utility.helper_functions as pjh
import phylojunction.calculation.discrete_sse as sseobj
import phylojunction.distribution.dn_discrete_sse as distsse

__author__ = "Fabio K. Mendes"
__email__ = "fmendes@lsu.edu"


class TestBiSSETrees(unittest.TestCase):

    @classmethod
    def setUpClass(cls) -> None:
        total_n_states = 2

        rates_t0_s0 = [
            sseobj.DiscreteStateDependentRate(
                name="lambda0",
                val=0.5,
                event=sseobj.MacroevolEvent.W_SPECIATION,
                states=[0,0,0]),
            sseobj.DiscreteStateDependentRate(
                name="mu0",
                val=0.25,
                event=sseobj.MacroevolEvent.EXTINCTION,
                states=[0]),
            sseobj.DiscreteStateDependentRate(
                name="q01",
                val=0.8,
                event=sseobj.MacroevolEvent.ANAGENETIC_TRANSITION,
                states=[0,1]) ]

        rates_t0_s1 = [
            sseobj.DiscreteStateDependentRate(
                name="lambda1",
                val=1.5,
                event=sseobj.MacroevolEvent.W_SPECIATION,
                states=[1,1,1]),
            sseobj.DiscreteStateDependentRate(
                name="mu1",
                val=0.25,
                event=sseobj.MacroevolEvent.EXTINCTION,
                states=[1]),
            sseobj.DiscreteStateDependentRate(
                name="q10",
                val=0.4,
                event=sseobj.MacroevolEvent.ANAGENETIC_TRANSITION,
                states=[1,0]) ]

        rates_t0 = rates_t0_s0 + rates_t0_s1

        matrix_atomic_rates = [rates_t0] 

        state_dep_par_manager = \
            sseobj.DiscreteStateDependentParameterManager(
                matrix_atomic_rates,
                total_n_states)

        event_handler = sseobj.MacroevolEventHandler(state_dep_par_manager)

        cls.sse_stash = sseobj.SSEStash(event_handler)


    def test_tree_size_state_count_max_taxa_bisse(self):
        """
        Test if BiSSE trees simulated with PJ have similar root ages
        and number of tips for both states as BiSSE trees simulated
        with diversitree.
        """
        stop_condition = "size"
        stop_condition_value = [ 50 ] # 50 living taxa

        start_at_origin = True

        # simulation initialization
        n_batches = i = 100
        n_sim = 100
        start_states_list = [0 for i in range(n_sim)]

        # simulations
        print(("\n\nRunning TestBiSSETrees.test_tree_size_state_count_max"
               "_taxa_bisse"))
        sim_batches = list()
        for i in range(n_batches):
            sse_sim = distsse.DnSSE(
                self.sse_stash,
                origin=start_at_origin,
                start_states_list=start_states_list,
                n=n_sim,
                stop=stop_condition,
                stop_value=stop_condition_value,
                condition_on_speciation=True,
                condition_on_survival=True,
                epsilon=1e-12,
                runtime_limit=3600,
                debug=False)

            trs = sse_sim.generate()

            sim_batches.append(trs)

            # printing progress
            pjh.print_progress(i , n_batches)


        # "expectations" from diversitree
        n0_mean_maxtaxa_divtree = [ 
            12.06, 13.43, 13.42, 12.35, 12.45, 12.13,  12.5, 12.72, 11.99, 12.29,
            12.36, 12.81, 11.59, 12.04,  12.9, 12.56, 12.06, 11.92, 12.21, 12.15,
            12.17,  12.6, 12.51, 12.11, 13.06, 12.37, 12.34, 12.79, 11.77, 12.51,
            12.51, 12.53, 12.35, 12.59, 12.26, 12.71, 12.19, 12.24, 12.23, 12.32,
            12.06, 11.93,  12.4, 12.23, 12.51, 12.81, 11.51, 12.16,  12.4, 12.16,
            12.72, 12.83, 12.14, 12.59, 12.46, 12.86, 13.22, 12.74, 11.87,  12.8,
            11.85, 12.67, 12.37, 12.77, 12.17, 12.55, 12.81, 13.14, 12.69,  13.0,
            12.45, 12.19, 12.72,  12.9, 12.34,  12.7, 12.48, 12.31, 12.28, 12.47,
            13.35, 12.75, 12.21, 12.97, 12.24, 12.18, 12.26, 12.49, 12.23, 12.36,
            13.26,  13.1, 12.63, 12.75, 12.02, 12.06, 12.64, 11.94, 12.81, 12.16
        ]

        n0_ci_width_maxtaxa_divtree = [
            0.807321611405531, 0.895915987222918, 0.946843718825228, 0.836384982160631, 0.892496464639463, 0.723592977837512, 0.803312636474239, 0.812238126003288, 0.838294274889869, 0.827016405978488,
            0.768765931197329, 0.830948382599939, 0.735528416112587, 0.810248281190227, 0.889502323683521, 0.983880639882376, 0.826793503466094, 0.762257292122509, 0.858781527232131,  0.89119116746616,
            0.780555006894376,  0.85411706433952, 0.850700577836475,  0.86085752783624, 0.967775718236452, 0.812035060197464, 0.845541933921791, 0.901546199536657, 0.763363710921285, 0.685502041491409,
            0.886003305281474, 0.920778670035269, 0.758529438513041, 0.840882479207638, 0.796579682331297, 0.900684955066325, 0.846218577680803, 0.852352489326221, 0.855866666666667, 0.873825145812677,
            0.810200388231769, 0.929085442352811, 0.824057465878871, 0.879568095112251, 0.845209147746487, 0.832814227650908, 0.861127941549438, 0.813498387623889, 0.897968796560303,  0.82862187548076,
             0.91215169392337, 0.816989705806591, 0.826417952810364, 0.805529194525141, 0.775249959299541, 0.940197618599516, 0.857553877095318, 0.793162420092168, 0.692911446234822, 0.739162288600338,
            0.767176826310781, 0.830182173995685, 0.721444713296476, 0.876474453867783, 0.910011487106985, 0.792553092253386, 0.846677011628521, 0.944315821164196,  0.81634349969541, 0.734422331915473,
            0.794020557960968, 0.766832803079516, 0.817950941701459, 0.952286944278049, 0.764057265821447,  0.79895335573891, 0.701717714690965, 0.774386053593426, 0.863640980402333, 0.879391608793353,
            0.819985981649083, 0.805664065140419, 0.811863011781538, 0.849312773643931, 0.791448255932238,  0.88200219977279, 0.887808081909052,  0.78523392795788, 0.791808538228858,  0.71597269176522,
            0.942506028917298, 0.939156982822275, 0.890006043334379, 0.885512646134795, 0.883584614687881, 0.884392311137993, 0.818633800544504, 0.768917343550086, 0.742145989068801, 0.796141141742943
        ]

        n1_mean_maxtaxa_divtree = [ 
            49.04,  49.1, 48.66, 49.38, 49.87, 50.29, 48.78, 48.86, 49.32, 49.12,
            48.94, 49.04, 49.76, 49.42, 48.72,  48.9, 49.65, 49.25, 49.59, 49.37,
            50.06, 48.94, 48.55, 48.88, 48.94, 50.46, 49.56, 49.06, 49.11, 48.56,
            48.25, 49.58, 48.83, 49.05,  49.0, 48.72, 49.53, 49.05, 49.69, 50.59,
            49.28, 48.92, 49.88, 49.89, 50.17, 49.17, 49.41, 48.89, 50.11, 49.94,
            48.83, 48.63, 48.99, 50.27, 48.81, 48.62, 48.48, 48.15, 50.48, 49.46,
            49.55,  48.4,  49.0, 48.69, 49.17, 50.29, 49.42, 48.85, 49.06, 49.06,
            50.06,  49.6,  49.6, 49.66, 50.39, 49.16, 49.71, 50.06, 48.92, 49.59,
            49.75, 49.44,  50.0, 48.79, 49.06, 49.52,  49.4,  49.7, 49.21, 49.37,
            49.16,  48.7, 49.04, 48.72, 49.78, 49.79,  49.5, 49.61, 48.92, 49.15
        ]

        n1_ci_width_maxtaxa_divtree = [ 
            0.895770879796228, 0.975634722173781,   1.0099212442542, 0.895276905191738, 0.973307206789243, 0.894199181481182,  1.01312830142658, 0.929822317337125, 0.894886732632921, 0.774679138271276,
            0.883075035076724, 0.908246559324753, 0.973492576429649, 0.999477710879765, 0.999827068683744,  1.08450470937035,  1.04472993948898,  1.03089556045492,  1.07284741042732, 0.963288585308608,
            0.871129918988882, 0.936807183679726, 0.971199882681836, 0.970362469816451, 0.919667541159644,  0.99294943486666,  1.00417992190422,  1.01589752445905, 0.916740981230067, 0.898539028792087,
            0.937451068083292,  1.03869726955344,  1.02851389191327, 0.979554444168832,  1.07208939774987,  1.10596810080581,  1.07080552403461, 0.729850987973891, 0.960076713264901,  1.10105038481397,
            0.871690998907157,  1.01124980058931,  1.15175492342614,  1.06040539532022,  1.13914432604282,  1.06920985533593, 0.991257802343604, 0.983713006517612, 0.875162995727601, 0.872465230356957,
            0.938973099947087, 0.925067231845526,  1.04514585440009, 0.959737144990621, 0.996178020681489, 0.910320583779432,  1.05106398886408, 0.902010929629856, 0.976024420956274,  1.00769025544018,
            0.906730724118581, 0.992394349769559, 0.858196182040381, 0.991101204914936,   1.1178250245251, 0.895933311912187,  1.02174770927766, 0.977571739081335,  0.97893825312649, 0.912891610944768,
              1.1051748718523,  0.95655594333383,  0.90399839098811,  1.01833918581862,  1.06478757135929, 0.981669527392035, 0.953018094945174, 0.991150144269196, 0.915379116469612,  1.06413139539564,
             1.03690064591307,  1.01264559117854, 0.951267692199321, 0.971167918362507, 0.949561058828996,  1.01002113869268, 0.959795769521787, 0.962420097693541, 0.988199676383291, 0.850773557452632,
            0.900437194821594,  1.00347638162781, 0.960975585096476,  1.00176572676107, 0.988704133824889,  1.03349704196847, 0.918265646738349, 0.982133878837094,  1.06834756649149,  1.04584362544737
        ]

        root_ages_mean_maxtaxa_divtree = [
            3.58470486649072, 3.85648231959483, 3.92664981403215, 3.68434412057153, 3.61223060143573,  3.7452459151552, 3.90701900083615,  3.7665154330315, 3.55863820275877, 3.65311801813261,
            3.67536474959641, 3.64862531645963, 3.56622872783685, 3.67970361037552, 3.84644392824895, 3.85599604802811, 3.75087838718462,  3.5238275460679, 3.91273915002564, 3.83723690356374,
            3.97276063269651, 3.68559273347413, 3.65435672063552, 3.75721161327911, 4.16839306091379, 4.03346928004677, 3.94400950964582, 3.75794574247996, 3.52046092543478, 3.62512969791142,
            3.47864518191692, 3.80655942399656, 3.66304010962505, 3.80669984304012,  3.6642584264747, 3.81121892251389, 3.72417088963375, 3.75738192869735, 3.88342446010284, 3.72772876708416,
            3.61770546813194, 3.97052013870234, 3.89747946825895, 3.94995345057505, 3.86518800145563, 3.57420897607499, 3.60359823130246, 3.52785867559025, 3.75879564346931, 4.04197168826976,
            3.68516356799044, 3.73052188612906, 3.64989131588324, 3.70445177562523, 3.77360283581457, 3.72594592122096, 3.96390199660019, 3.69527498262034, 3.68089623718042, 3.86416701513073,
            3.67010142438172, 3.69627874115369, 3.59360409029375, 3.89275023918086,  3.5932595515561, 3.82250111786492, 3.88768221962166, 3.72389032277866, 3.76772369292507, 3.77911840953097,
            3.88060116318094, 3.61123239688932, 3.84259880893419, 3.95622110628335, 4.07911697105035, 3.63466187512128, 3.77548614081508, 3.93399335737022, 3.92033212217965, 3.75252935325512,
            3.86221369432535, 3.76495651099909, 3.90129707767582, 3.56930476797875, 3.48979251174607, 3.68836590187977, 3.67405199553446, 3.83940936528713, 4.00296821261155, 3.75910889969945,
            3.81057146446284, 3.69346930467036, 3.66916284438519, 3.74477895803878, 3.72956960160533, 4.04022901958621, 3.74396968836535, 3.88997600068596, 3.65947470213476,  3.5366200517605
        ]

        root_ages_ci_width_maxtaxa_divtree = [
            0.235981773497821, 0.297611003879549, 0.313392315848257, 0.222901161852675, 0.221304112196267, 0.200804659966008, 0.269289633469476, 0.240916281204074,  0.23503869994916, 0.259639200403118,
            0.225924219945494, 0.219262678273446, 0.213858901666017, 0.297585618566676, 0.247833980528564, 0.261814541099581, 0.250140579078444, 0.227262975001599, 0.264444527079894,  0.24548251496905,
            0.315155377165016,  0.23773591902188, 0.214339558808837, 0.293664442150587, 0.364639679540608, 0.340743266803247, 0.304090755283025, 0.312687981021168, 0.250059656251136,  0.21885499662489,
            0.226297841766879, 0.256158970342579, 0.241717778763277, 0.216916813277661,  0.22385189294512, 0.308934587059059,  0.22691575626854, 0.288971585364274, 0.293269963702371, 0.212153223774821,
            0.247718960787939, 0.342280857063048, 0.282928481390018, 0.249671228228243, 0.268597326550702, 0.192888020680325, 0.27491469680958,  0.186448989305866,  0.23712893282878, 0.302271877649638,
            0.285730805047942, 0.222022469366532, 0.208187059518006, 0.251901509793011, 0.270003160749116, 0.243540725198734, 0.351429262566948, 0.254749202625518, 0.210088114830728, 0.281348236108657,
            0.262923768248279, 0.231129559314536, 0.229906761532691, 0.255410339793626, 0.255706891275041, 0.261304527137823, 0.267359412156601, 0.226901126196345, 0.249303280057289, 0.259029462170184,
            0.266474421352723, 0.229068007756948, 0.243765041668662,  0.31955680888726, 0.278484940625921, 0.210729085099957, 0.241768871570644, 0.239317077031962, 0.251404810421938, 0.222450009673356,
            0.238258795470196,  0.23401031124416, 0.270562077243049, 0.240563530484513, 0.214953435719659, 0.237058897207485, 0.275607578383277, 0.276542232915678, 0.281455567529177, 0.234305348235221,
            0.226721794503856, 0.258408045027351, 0.210593647376075, 0.220252290480401, 0.246032098286528, 0.339821674447148, 0.245674275468596, 0.293126844889882, 0.224168068761252, 0.200622710415004
        ]

        # parsing simulations
        n0_ci_overlap_count = 0
        n1_ci_overlap_count = 0
        root_age_ci_overlap_count = 0
        global_mean_n0 = 0.0
        global_mean_n1 = 0.0
        global_mean_root_age = 0.0
        for i, batch in enumerate(sim_batches):
            n0s = [ann_tr.state_count_dict[0] for ann_tr in batch]
            n1s = [ann_tr.state_count_dict[1] for ann_tr in batch]
            root_ages = [ann_tr.root_age for ann_tr in batch]

            mean_n0 = statistics.mean(n0s)
            mean_n1 = statistics.mean(n1s)
            mean_root_ages = statistics.mean(root_ages)

            global_mean_n0 += mean_n0
            global_mean_n1 += mean_n1
            global_mean_root_age += mean_root_ages

            stdevs_n0 = statistics.stdev(n0s)
            stdevs_n1 = statistics.stdev(n1s)
            stdevs_root_ages = statistics.stdev(root_ages)

            sterr_n0 = stdevs_n0 / math.sqrt(n_sim)
            sterr_n1 = stdevs_n1 / math.sqrt(n_sim)
            sterr_root_ages = stdevs_root_ages / math.sqrt(n_sim)

            n0_ci_width_maxtaxa = 1.96 * sterr_n0
            n1_ci_width_maxtaxa = 1.96 * sterr_n1
            root_ages_ci_width_maxtaxa = 1.96 * sterr_root_ages

            if abs(mean_n0 - n0_mean_maxtaxa_divtree[i]) \
                    <= (n0_ci_width_maxtaxa + n0_ci_width_maxtaxa_divtree[i]):
                n0_ci_overlap_count += 1

            if abs(mean_n1 - n1_mean_maxtaxa_divtree[i]) \
                    <= (n1_ci_width_maxtaxa + n1_ci_width_maxtaxa_divtree[i]):
                n1_ci_overlap_count += 1

            if abs(mean_root_ages - root_ages_mean_maxtaxa_divtree[i]) \
                <= (root_ages_ci_width_maxtaxa \
                    + root_ages_ci_width_maxtaxa_divtree[i]):
                root_age_ci_overlap_count += 1

        print("\n\nPJ global mean n0 taxon count = " \
              + str(global_mean_n0 / 100.0))
        print("diversitree global mean n0 taxon count = " \
              + str(statistics.mean(n0_mean_maxtaxa_divtree)))
        print("\nPJ global mean n1 taxon count = " \
              + str(global_mean_n1 / 100.0))
        print("diversitree global mean n1 taxon count = " \
              + str(statistics.mean(n1_mean_maxtaxa_divtree)))
        print("\nPJ global mean root age = " \
              + str(global_mean_root_age / 100.0))
        print("diversitree global mean root age = " \
              + str(statistics.mean(root_ages_mean_maxtaxa_divtree)))

        # [==== * ====][.... + ....] if we take '+' to be the "truth" of the '*' interval, + cannot be more than '====' away from '*' 95% of the time
        # then abs('+' - '*') can be at most ('====' + '....'). '....' can be added because we still are guaranteed to see '+' falling within that range
        # 95% of the time

        print("\n95% CIs of simulations here and from diversitree " \
              + "overlapped " + str(n0_ci_overlap_count) \
                + " times for state 0 count.")
        print("\n95% CIs of simulations here and from diversitree " \
              + "overlapped " + str(n1_ci_overlap_count) \
                + " times for state 1 count.")
        print("\n95% CIs of simulations here and from diversitree " \
              + "overlapped " + str(root_age_ci_overlap_count) \
                + " times for root age.")
        exp_count = int(0.95 * n_batches)
        a_delta = math.ceil(0.07 * exp_count)
        self.assertAlmostEqual(
            n0_ci_overlap_count,
            exp_count,
            msg="Mean absolute difference must be 1.96 * (stderr_python " \
                + "+ stderr_divtree) apart " + str(exp_count) + " (+/- " \
                + str(a_delta) + ") out of 100 times.",
            delta=a_delta)
        self.assertAlmostEqual(
            n1_ci_overlap_count,
            exp_count,
            msg="Mean absolute difference must be 1.96 * (stderr_python " \
                + "+ stderr_divtree) apart " + str(exp_count) + " (+/- " \
                + str(a_delta) + ") out of 100 times.",
            delta=a_delta)
        self.assertAlmostEqual(
            root_age_ci_overlap_count,
            exp_count,
            msg="Mean absolute difference must be 1.96 * (stderr_python " \
                + "+ stderr_divtree) apart " + str(exp_count) + " (+/- " \
                + str(a_delta) + ") out of 100 times.", delta=a_delta)


    def test_tree_size_state_count_max_t_bisse(self):
        """
        Test if BiSSE trees simulated with PJ have similar root ages
        and number of tips (for both states) as BiSSE trees simulated
        with diversitree.
        """

        stop_condition = "age"
        stop_condition_value = [ 3.0 ] # 3.0 time units

        start_at_origin = True

        # simulation initialization
        n_batches = i = 100
        n_sim = 100
        start_states_list = [0 for i in range(n_sim)]
        # seeds_list = [i+1 for i in range(n_sim)]

        # simulations
        print(("\n\nRunning TestBiSSETrees.test_tree_size_state_count_max_t_"
              "bisse"))
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
                debug=False)

            trs = sse_sim.generate()

            sim_batches.append(trs)

            # printing progress
            pjh.print_progress(i , n_batches)


        # "expectations" from diversitree
        n0_mean_maxt_divtree = [
            5.22, 4.36, 4.35, 5.46, 4.46, 3.87, 5.07, 4.19, 5.59,  4.7,
             5.0, 3.62, 4.37, 4.28, 4.45, 4.48, 5.56, 4.57, 4.89, 4.15,
            5.69, 5.66, 5.03, 3.93, 4.88, 4.21, 5.02, 3.92, 4.03, 5.03,
             4.4, 4.07, 5.19, 4.06, 4.35, 4.41, 4.25, 5.05,  4.5,  4.0,
            4.28, 4.39, 4.34, 4.32, 4.55, 5.49, 4.03, 4.26, 4.49, 4.87,
            4.54, 4.26, 4.57,  5.0,  4.4, 4.73, 5.03, 4.36, 4.23, 4.57,
             4.4, 4.43, 4.52, 3.57, 4.96, 4.62, 4.51, 4.32, 4.76, 5.47,
            4.42, 4.45, 3.61, 3.83, 4.78, 3.74, 4.11, 5.24, 5.63, 4.58,
             5.6,  4.0, 5.25, 4.93, 3.98, 3.79, 4.61, 4.17, 4.05,  4.3,
            4.75, 4.73, 3.92, 4.22, 4.48, 4.87,  5.5,  4.1, 5.21, 5.32
        ]

        n0_ci_width_maxt_divtree = [
             0.98161418576596, 0.801871854544796,  1.14739883249378,   1.3215550583011, 0.982341284067731,  1.29132873765237,  1.07254354605387, 0.979681200671675,  1.20246092374402,  1.09199274796867,
            0.865848751576198, 0.630224274333825, 0.857587813685466,  0.97188087227756, 0.918635329544248,  0.99804017732532,  1.19476397330987,  1.02670133825753, 0.879144468425774, 0.800348484705067,
             1.29823372722521,  1.34420491703828,   1.1508027525834, 0.974741406577647,  0.82654002549022, 0.985447124277845,  1.05566875429174, 0.811569011361338,  1.02223181541392,  1.15147693722875,
              1.0303778653749, 0.941943875314461,   1.0223229156379,  0.69802510909658,  0.86601679597401, 0.903180301473464, 0.898994522516768,  1.03277590047699, 0.932522678252141, 0.890156449760661,
            0.947212489680146, 0.905668772596563,  1.03346888182324, 0.835226625993016,   1.1203631986859,  1.25871565069039, 0.738350759490639, 0.761207890742925, 0.949019439934576,  1.00778074495348,
            0.905465232503561,   1.0012000266347,  1.05393051809391,  1.01022089793582, 0.811243813691374,  1.06440849680976,   1.1977272018707, 0.936268548197276, 0.988812058211593, 0.954629126568496,
            0.923112521458098,  1.13135103002177, 0.880284706928687, 0.756931659378637,  1.05159548546506, 0.887878011458964,  1.01041101667899, 0.956101491590836,  0.93941311917859,  1.34917108319041,
            0.954192057409955, 0.938692040807179, 0.688552020730841, 0.836700408350165,  1.03360404352338, 0.674732927636086, 0.770971082794323,  1.31469581626053,  1.35811499252799,  1.05757843160207,
             1.11119060725717, 0.825938878067276,  1.27795346179157, 0.951779493179152, 0.772019778885165, 0.861938673672157, 0.797200535206377, 0.915114132971203, 0.794997363568464, 0.946974853804141,
             1.12623579711669,  1.36387432164864, 0.746309369057788, 0.895190214875074,  1.08019534035902,  1.08456016761253,  1.11345816756136, 0.743089156797147,  1.24001815496648,  1.16402060171293
        ]

        n1_mean_maxt_divtree = [
            18.78, 16.13, 16.33, 22.09, 16.78, 14.22, 19.43, 16.52, 19.33, 16.52,
            17.54, 14.32, 16.62, 16.44, 15.34, 15.74, 21.75, 16.59,  18.2, 14.68,
            21.18,  21.4, 18.59, 15.36, 16.66, 14.99, 18.96, 14.44, 13.92, 18.72,
            17.36, 15.63, 18.77, 14.61, 14.21, 17.18, 15.11, 18.28, 16.11, 15.75,
            14.98, 17.23, 16.32, 14.91, 15.82, 20.33, 17.79, 16.24, 15.07, 17.83,
             16.3, 16.38, 15.68, 18.61, 16.46, 17.88, 16.26, 15.32, 15.33,  15.9,
            17.57,  15.4,  17.8, 15.61, 18.37, 16.53, 18.69, 17.34, 16.26, 21.91,
            16.05, 17.75, 14.45, 13.67, 19.23, 14.97, 16.29,  21.1, 22.29, 16.98,
             20.3, 12.44, 19.26, 20.15, 14.31, 14.24,  18.2, 16.39, 16.45,  16.6,
            18.19, 16.12, 14.67, 16.43,  15.8, 17.94, 21.41, 15.57, 20.64, 17.34
        ]

        n1_ci_width_maxt_divtree = [
            4.08810154571542, 3.19892697438919, 4.06911866684279, 5.48729960177098, 4.06334754088052, 3.70405273344228, 4.79066379087196, 3.73855744793565, 4.73246052746776, 4.43554066643925,
             3.9676629223071, 2.74058696831702, 3.92591315628696, 4.25285922388478,  3.5725088700613, 3.52850811155839, 4.83758531885193, 3.01037694045719, 3.60750593495105, 3.09275004322997,
            5.15093870052495, 5.16618030946453, 4.81521172430652, 4.07164164722712, 3.48174567709935, 3.95846632866209, 4.12733711379493, 3.22702659249763, 3.80105880751427, 4.81327605371777,
            4.23779582578138, 3.90657923803298, 3.84371406477084, 3.15977817377656, 2.91305621287917, 3.80347751055969, 3.72612455080085, 4.32297211939315, 3.87688992900445,  4.3483270067436,
            3.27378381051342, 3.83461739147335, 4.67187366232946, 3.55801965552341, 4.21142132988113,  5.4776709363072, 3.49476932294623, 3.59890340739475, 3.52307511476018, 4.21958393487253,
            3.97011500736656, 4.32168832627886,  4.4048499580725, 4.55044289623759, 3.67437446258722, 5.07159614362583, 4.04890554090756, 3.63204191138547, 3.50461089896315, 3.69239754034687,
            3.71082196802593, 4.11727253056062, 3.79240418934395, 3.14451333075131, 4.08801374426913, 3.64734893423714, 4.32966456143492, 3.91635745915497,  3.8022775362972, 5.95521942622986,
            3.79244255930878, 4.38228372270243, 3.29038269856719, 3.39733961537881, 4.50587810404197, 2.92663303383726, 3.21340692962369, 5.81873606324753, 4.74901434264697, 4.35845321390259,
            4.37147867177432, 2.65346065295261, 4.74740643265011,  4.0489778982111, 3.94119533505901, 3.74432810096508, 3.80680403084962, 3.71684449179161,  3.9726432337458, 3.85330741314549,
            4.74684036777852, 5.10735443679311,  3.3016471214755, 3.85121775855876, 4.63536117821743, 3.79074828337082, 4.63905103827901, 3.26717206999409, 4.74895101722901,  3.9566744620178
        ]

        root_ages_mean_maxt_divtree = [
            2.16019074739235, 2.23985325286739, 2.10912011582493, 2.08439540715148, 2.13755192815195, 2.04205603439848, 2.08936733990997, 2.08124115321389, 2.08807639506063, 2.11766684110982,
             2.0757680302778, 2.00991995256452, 2.09897928747077, 2.04872196530228, 1.94931979163084, 1.97173449766771, 2.13918824968117, 2.03011397192554, 2.09541540835165, 2.06150834429441,
            2.17898464720682, 2.02624688600771, 2.09667511152114,  1.9754571847295, 2.07751262302644,   2.158776386014, 2.11504982259296, 2.10446610835592, 1.98656140946664, 2.13145930566523,
            2.05392749413261, 1.99715972917505,   2.094930002603,  2.0021278866453, 1.99800273456582, 2.06648552997892, 1.98789166119901, 2.09592859576434, 2.00193105859188, 2.00967553866122,
            2.12256784098586, 2.09696623504926,  1.9873263882158, 2.00317034938936, 2.05783703333234, 2.16276211673741, 2.01938809176279,  2.0778645042394, 2.07649292566626, 2.19972519294765,
            2.10734467222416, 1.95657366059353, 2.05713384191161, 2.11220638537172, 2.12465050375679, 2.22008365896472, 2.10599092168607, 1.92049420120249, 2.01865311783623, 1.98159012747726,
            2.03027863148591, 1.95820335112912, 1.96147032624551,  2.0424904697804, 2.02598968913367, 2.05281521660041, 1.87214418547339, 2.04883756427977,  2.1634103222118, 2.05967692385619,
            2.05446496840829, 2.00338268809876, 2.04269943590095, 1.96483549376368, 1.94699243666359, 2.00933205297941, 2.09417603053312, 1.92748936307839, 2.04613235181204, 2.02311834630362,
            2.17709563733529, 1.92738468601729, 2.16764177144626, 2.10799911987284,  2.0678850656439, 1.92466950356008, 2.07922029495231, 1.98672184902005,  2.0120457850184, 2.05198182016563,
            1.99223230590896, 1.98443208217947,  2.1411956180265, 2.07942520081614, 1.98297995887198, 2.09256649157022, 2.13358578217028, 2.05623500706743, 2.15917821269505, 2.04823393358692
        ]

        root_ages_ci_width_maxt_divtree = [
            0.131509336256971, 0.113368190652701, 0.132905496094447, 0.150075373915352, 0.124348500246804, 0.134352098826473, 0.142229926075057, 0.147769194139888, 0.142024033059407, 0.142960742226984,
            0.135376437999344, 0.127962850153195, 0.141706314053254, 0.138853892488107, 0.167679239140398, 0.141929760304129, 0.130405591631781, 0.140592563673068,  0.13479959633956, 0.140810658626993,
            0.133925886541875, 0.140655099303131, 0.125554463378381, 0.163574574741462, 0.137954565306727, 0.133925102730878, 0.133117618975698, 0.123658101466697, 0.145921842031752, 0.136109701567449,
             0.12919792470185, 0.146469249860248, 0.138961220099624, 0.141154937526522, 0.144905559356234, 0.136989581099801, 0.153565856318093,  0.14619780010837, 0.142149733645523, 0.125054194902086,
            0.126578945088363, 0.130274514759025, 0.153630694718841, 0.137166225385444, 0.151670094693458, 0.131659161345938, 0.129825151242959, 0.153244384030126, 0.125697814332391, 0.122163352340789,
             0.14169065050112,  0.14423459316781, 0.133861971841433, 0.129031343436785, 0.133393151436868, 0.124243690494388, 0.144132490964085, 0.149497756055504, 0.141231842519813, 0.154111281742736,
            0.150148408014678, 0.155317703312488, 0.160034800733246, 0.143822107792082, 0.136440771932919, 0.144500713812307, 0.151579473914853, 0.139653180814706, 0.140066921781386, 0.141832052213146,
            0.135947837835367, 0.137563663037994, 0.131445063292238, 0.157989987160941, 0.157077862816222,  0.13224158335642, 0.140424798436539, 0.146143851988532, 0.146891598902084, 0.150635164827912,
            0.127420233505672, 0.157284650622664, 0.138457705182025, 0.140304209487322, 0.138073021199882, 0.149866940133464, 0.136411943841081,   0.1514282969283, 0.137968831355123, 0.141772791438584,
            0.163207912189686, 0.148880499430412, 0.141555578522862, 0.143785731101053, 0.143693129538231, 0.126868118721876, 0.119099619176478, 0.138228325030658, 0.138096034395849, 0.148795113183592
        ]

        # parsing simulations
        n0_ci_overlap_count = 0
        n1_ci_overlap_count = 0
        root_age_ci_overlap_count = 0
        global_mean_n0 = 0.0
        global_mean_n1 = 0.0
        global_mean_root_age = 0.0
        for i, batch in enumerate(sim_batches):
            n0s = [ann_tr.state_count_dict[0] for ann_tr in batch]
            n1s = [ann_tr.state_count_dict[1] for ann_tr in batch]
            root_ages = [ann_tr.root_age for ann_tr in batch]

            mean_n0 = statistics.mean(n0s)
            mean_n1 = statistics.mean(n1s)
            mean_root_ages = statistics.mean(root_ages)

            global_mean_n0 += mean_n0
            global_mean_n1 += mean_n1
            global_mean_root_age += mean_root_ages

            stdevs_n0 = statistics.stdev(n0s)
            stdevs_n1 = statistics.stdev(n1s)
            stdevs_root_ages = statistics.stdev(root_ages)

            sterr_n0 = stdevs_n0 / math.sqrt(n_sim)
            sterr_n1 = stdevs_n1 / math.sqrt(n_sim)
            sterr_root_ages = stdevs_root_ages / math.sqrt(n_sim)

            n0_ci_width_maxt = 1.96 * sterr_n0
            n1_ci_width_maxt = 1.96 * sterr_n1
            root_ages_ci_width_maxt = 1.96 * sterr_root_ages

            if abs(mean_n0 - n0_mean_maxt_divtree[i]) \
                    <= (n0_ci_width_maxt + n0_ci_width_maxt_divtree[i]):
                n0_ci_overlap_count += 1

            if abs(mean_n1 - n1_mean_maxt_divtree[i]) \
                    <= (n1_ci_width_maxt + n1_ci_width_maxt_divtree[i]):
                n1_ci_overlap_count += 1

            if abs(mean_root_ages - root_ages_mean_maxt_divtree[i]) \
                <= (root_ages_ci_width_maxt \
                    + root_ages_ci_width_maxt_divtree[i]):
                root_age_ci_overlap_count += 1

        print("\n\nPJ global mean n0 taxon count = " \
              + str(global_mean_n0 / 100.0))
        print("diversitree global mean n0 taxon count = " \
              + str(statistics.mean(n0_mean_maxt_divtree)))
        print("\nPJ global mean n1 taxon count = " \
              + str(global_mean_n1 / 100.0))
        print("diversitree global mean n1 taxon count = " \
              + str(statistics.mean(n1_mean_maxt_divtree)))
        print("\nPJ global mean root age = " \
              + str(global_mean_root_age / 100.0))
        print("diversitree global mean root age = " \
              + str(statistics.mean(root_ages_mean_maxt_divtree)))

        # [==== * ====][.... + ....] if we take '+' to be the "truth" of the '*' interval, + cannot be more than '====' away from '*' 95% of the time
        # then abs('+' - '*') can be at most ('====' + '....'). '....' can be added because we still are guaranteed to see '+' falling within that range
        # 95% of the time

        print("\n95% CIs of simulations here and from diversitree " \
              + "overlapped " + str(n0_ci_overlap_count) \
              + " times for state 0 count.")
        print("\n95% CIs of simulations here and from diversitree " \
              + "overlapped " + str(n1_ci_overlap_count) \
              + " times for state 1 count.")
        print("\n95% CIs of simulations here and from diversitree " \
              + "overlapped " + str(root_age_ci_overlap_count) \
              + " times for root age.")
        exp_count = int(0.95 * n_batches)
        a_delta = math.ceil(0.07 * exp_count)
        self.assertAlmostEqual(
            n0_ci_overlap_count,
            exp_count,
            msg="Mean absolute difference must be 1.96 * (stderr_python " \
                + "stderr_divtree) apart " + str(exp_count) + " (+/- " \
                + str(a_delta) + ") out of 100 times.",
        delta=a_delta)
        self.assertAlmostEqual(
            n1_ci_overlap_count,
            exp_count,
            msg="Mean absolute difference must be 1.96 * (stderr_python " \
                + "stderr_divtree) apart " + str(exp_count) + " (+/- " \
                + str(a_delta) + ") out of 100 times.",
            delta=a_delta)
        self.assertAlmostEqual(
            root_age_ci_overlap_count,
            exp_count,
            msg="Mean absolute difference must be 1.96 * (stderr_python " \
                + "+ stderr_divtree) apart " + str(exp_count) + " (+/- " \
                + str(a_delta) + ") out of 100 times.",
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
    # $ python3 tests/distribution/test_dn_discrete_sse_bisse.py
    # 
    # or
    #
    # $ python3 -m tests.distribution.test_dn_discrete_sse_bisse
    #
    # or 
    #
    # $ python3 -m unittest tests.distribution.test_dn_discrete_sse_bisse.TestBiSSETrees.test_tree_size_state_count_max_taxa_bisse

    # total_n_states = 2

    # rates_t0_s0 = [ sseobj.DiscreteStateDependentRate(name="lambda0", val=0.5, event=sseobj.MacroevolEvent.W_SPECIATION, states=[0,0,0]),
    #                 sseobj.DiscreteStateDependentRate(name="mu0", val=0.25, event=sseobj.MacroevolEvent.EXTINCTION, states=[0]),
    #                 sseobj.DiscreteStateDependentRate(name="q01", val=0.8, event=sseobj.MacroevolEvent.ANAGENETIC_TRANSITION, states=[0,1]) ]

    # rates_t0_s1 = [ sseobj.DiscreteStateDependentRate(name="lambda1", val=1.5, event=sseobj.MacroevolEvent.W_SPECIATION, states=[1,1,1]),
    #                 sseobj.DiscreteStateDependentRate(name="mu1", val=0.25, event=sseobj.MacroevolEvent.EXTINCTION, states=[1]),
    #                 sseobj.DiscreteStateDependentRate(name="q10", val=0.4, event=sseobj.MacroevolEvent.ANAGENETIC_TRANSITION, states=[1,0]) ]

    # rates_t0 = rates_t0_s0 + rates_t0_s1

    # matrix_atomic_rate_params = [ rates_t0 ] # 1D: time slices (i) , 2D: all rates from all states in i-th time slice

    # state_dep_par_manager = sseobj.DiscreteStateDependentParameterManager(matrix_atomic_rate_params, total_n_states)

    # event_handler = sseobj.MacroevolEventHandler(state_dep_par_manager)

    # # stop_condition = "size"
    # stop_condition = "age"
        
    # # stop_condition_value = [ 10 ] # 10 living taxa
    # stop_condition_value = [ 3.0 ] # origin age 

    # start_at_origin = True

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

    unittest.main()