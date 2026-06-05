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


class TestFBDTrees(unittest.TestCase):

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
                val=0.5,
                event=sseobj.MacroevolEvent.EXTINCTION,
                states=[0]),
            sseobj.DiscreteStateDependentRate(
                name="psi",
                val=0.75,
                event=sseobj.MacroevolEvent.ANCESTOR_SAMPLING,
                states=[0]) ]

        matrix_atomic_rate_params = [rates_t0_s0]

        state_dep_par_manager = sseobj.DiscreteStateDependentParameterManager(matrix_atomic_rate_params, total_n_states)

        event_handler = sseobj.MacroevolEventHandler(state_dep_par_manager)

        cls.sse_stash = sseobj.SSEStash(event_handler)

    
    # NOTE: activate this test after I eventually implement the GSM approach as a block in
    # dn_discrete_sse.is_tr_ok()
    # def test_tree_size_sa_count_max_taxa_fbd(self):
    #     """
    #     Test if FBD trees simulated here have similar root ages and number of sampled ancestors
    #     as trees simulated with FossilSim
    
    #     Note: FossilSim builds on top of TreeSim, which uses the GSA approach (see Stadler, 2011).
    #     PhyloJunction uses what that paper refers to as SSA (simple sampling approach).
    #     """

    #     stop_condition = "size"
    #     stop_condition_value = [ 10 ] # 10 taxa

    #     start_at_origin = True

    #     # simulation initialization
    #     n_batches = i = 100
    #     n_sim = 100
    #     start_states_list = [0 for i in range(n_sim)]

    #     # simulations
    #     sim_batches = list()
    #     for i in range(n_batches):
    #         # print("Doing batch " + str(n_batches - i))
    #         sse_sim = distsse.DnSSE(self.event_handler, stop_condition_value, n=100, stop=stop_condition, origin=start_at_origin,
    #                                 start_states_list=start_states_list, epsilon=1e-12, runtime_limit=3600,
    #                                 condition_on_speciation=True, condition_on_survival=True,
    #                                 debug=False)

    #         trs = sse_sim.generate()

    #         sim_batches.append(trs)

    #         # printing progress
    #         pjh.print_progress(i , n_batches)

    #     n_sa_mean_maxtaxa_fossilsim = [
    #         11.79, 11.43, 12.31, 11.86, 13.68, 11.66, 11.75, 11.7, 11.82, 12.38, 11.98, 12.36, 11.95, 12.93, 11.88, 11.39, 11.44, 11.86, 13.44, 10.6, 12.51, 10.84, 11.73, 12.21, 11.81, 11.73, 12.26, 12.11, 12, 12.6, 11.77, 12.97, 12.35, 11.2, 11.58, 11.06, 11.66, 12.77, 12.48, 11.33, 11.79, 10.97, 13.07, 11.86, 12.33, 12.28, 13.29, 11.29, 12.66, 11.6, 12.35, 11.26, 11.77, 11.69, 11.92, 11.3, 12.24, 11.4, 12.84, 11.8, 11.34, 12.29, 12.21, 12.32, 11.44, 11.65, 13.44, 12.05, 12.22, 11.61, 11.89, 10.98, 11.78, 11.9, 12.22, 13.22, 9.65, 12.55, 11.21, 12.72, 11.86, 10.57, 12.34, 11.67, 12.22, 11.75, 12.85, 11.86, 11.55, 12.26, 11.82, 12.35, 11.59, 12.7, 11.93, 12.52, 11.47, 12.49, 11.33, 12.63
    #     ]

    #     n_sa_ci_width_maxtaxa_fossilsim = [
    #         1.33379389560741, 1.24593720418313, 1.56990316319771, 1.11572807519475, 1.12088433773778, 1.20861690315247, 1.23313962442756, 1.28887585503685, 1.20931019723283, 1.27989681149724, 1.2932279844629, 1.32948845877661, 1.17248886622183, 1.41481811716081, 1.28519152744726, 1.18057597926956, 1.32504457801893, 1.31732008581242, 1.39633926241354, 1.14048736702492, 1.39471956502033, 1.11144201101019, 1.30216129336553, 1.20928292233977, 1.30062271891956, 1.14213804186126, 1.51504195308115, 1.22572874485328, 1.4264935607243, 1.37016038285577, 1.37718152878708, 1.19350802242631, 1.3215418451509, 1.47464518588114, 1.19504325503897, 1.25170581886068, 1.34910349222262, 1.68878193914425, 1.42961298485452, 1.17568814046769, 1.19637868333542, 1.18698769304249, 1.53782837886036, 1.21750969813407, 1.30905657199636, 1.22356935449223, 1.30705118584707, 1.29960793404096, 1.25677968722534, 1.44218465033316, 1.30350754473916, 1.20327547850512, 1.25295608184417, 1.08421663868249, 1.4198625262804, 1.2987731273643, 1.26536470504265, 1.40318103254748, 1.44640815676613, 1.40179763377622, 1.23025538549333, 1.34162600518094, 1.2137669947629, 1.32900970255862, 1.17643547492698, 1.21347283580323, 1.59701279932633, 1.10009136443343, 1.27223370955071, 1.25760534285754, 1.44538162672309, 1.10838639362623, 1.2651990962419, 1.30800234763327, 1.27680062023731, 1.49751380007469, 0.910574177831402, 1.51625933924193, 1.3502750690212, 1.33721810796999, 1.00406398765872, 1.2727444940249, 1.42246736498397, 1.20083338805123, 1.23477643075149, 1.25311711545669, 1.33410806280163, 1.18717238389836, 1.43094101975595, 1.33650406265561, 1.14265604181617, 1.40466666666667, 1.25920880658684, 1.27403807246919, 1.31680007057206, 1.37707586317141, 1.31980242843428, 1.29847282418675, 1.15739282838279, 1.27894901992048
    #     ]

    #     root_ages_mean_maxtaxa_fossilsim = [
            
    #     ]

    #     root_ages_ci_width_maxtaxa_fossilsim = [
            
    #     ]

    #     # parsing simulations
    #     n_sa_ci_overlap_count = 0
    #     root_age_ci_overlap_count = 0
    #     for i, batch in enumerate(sim_batches):
    #         n_sas = [ann_tr.n_sa for ann_tr in batch]
    #         root_ages = [ann_tr.root_age for ann_tr in batch]

    #         mean_sa = statistics.mean(n_sas)
    #         mean_root_ages = statistics.mean(root_ages)

    #         # debugging
    #         # print("PJ's vs. FossilSim mean sa count: " + str(mean_sa) + " <-> " + str(n_sa_mean_maxtaxa_fossilsim[i]))
    #         # print("PJ's vs. FossilSim mean root age: " + str(mean_root_ages) + " <-> " + str(root_ages_mean_maxtaxa_fossilsim[i]))

    #         stdevs_n_sa = statistics.stdev(n_sas)
    #         stdevs_root_ages = statistics.stdev(root_ages)

    #         sterr_n_sa = stdevs_n_sa / math.sqrt(n_sim)
    #         sterr_root_ages = stdevs_root_ages / math.sqrt(n_sim)

    #         n_sa_ci_width_maxtaxa = 1.96 * sterr_n_sa
    #         root_ages_ci_width_maxtaxa = 1.96 * sterr_root_ages

    #         if abs(mean_sa - n_sa_mean_maxtaxa_fossilsim[i]) <= (n_sa_ci_width_maxtaxa + n_sa_ci_width_maxtaxa_fossilsim[i]):
    #             n_sa_ci_overlap_count += 1

    #         if abs(mean_root_ages - root_ages_mean_maxtaxa_fossilsim[i]) <= (root_ages_ci_width_maxtaxa + root_ages_ci_width_maxtaxa_fossilsim[i]):
    #             root_age_ci_overlap_count += 1

    #     # [==== * ====][.... + ....] if we take '+' to be the "truth" of the '*' interval, + cannot be more than '====' away from '*' 95% of the time
    #     # then abs('+' - '*') can be at most ('====' + '....'). '....' can be added because we still are guaranteed to see '+' falling within that range
    #     # 95% of the time

    #     print("\n95% CIs of simulations here and from FossilSim overlapped " + str(n_sa_ci_overlap_count) + " times for the sampled ancestor count.")
    #     print("\n95% CIs of simulations here and from FossilSim overlapped " + str(root_age_ci_overlap_count) + " times for root age.")
    #     exp_count = int(0.95 * n_batches)
    #     a_delta = math.ceil(0.07 * exp_count)
    #     self.assertAlmostEqual(n_sa_ci_overlap_count, exp_count,
    #                             msg="Mean absolute difference must be 1.96 * (stderr_python + stderr_divtree) apart " + str(exp_count) + " (+/- " + str(a_delta) + ") out of 100 times.", delta=a_delta)
    #     self.assertAlmostEqual(root_age_ci_overlap_count, exp_count,
    #                             msg="Mean absolute difference must be 1.96 * (stderr_python + stderr_divtree) apart " + str(exp_count) + " (+/- " + str(a_delta) + ") out of 100 times.", delta=a_delta)

    def test_tree_size_sa_count_max_t_fbd(self):
        """
        Test if FBD trees simulated here have similar root ages and
        number of sampled ancestors as trees simulated with FossilSim,
        starting from the origin.

        Note: condition_on_speciation=True to match FossilSim!
        """

        stop_condition = "age"
        stop_condition_value = [ 3.0 ] # 3.0 time units

        start_at_origin = True

        # simulation initialization
        n_batches = i = 100
        n_sim = 100
        start_states_list = [0 for i in range(n_sim)]

        # simulations
        print("\n\nRunning TestFBDTrees.test_tree_size_sa_count_max_t_fbd")
        sim_batches = list()
        for i in range(n_batches):
            # print("Doing batch " + str(n_batches - i))
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


        # "expectations" from FossilSim
        n_sa_mean_maxt_fossilsim = [
            9.56, 8.73,  8.35, 9.66,  9.86,  9.6, 9.49,  8.24,  8.15,  8.41,
            9.31, 9.07,  8.43,  8.8,   6.8, 9.44, 8.42,  9.07,  9.23, 10.48,
            9.77, 7.74, 10.14, 9.47,  9.18, 9.19, 8.49,  9.78,  10.0,  8.84,
            7.85,  9.2,   9.0, 8.37,  8.69, 9.59, 8.72,  9.04,  9.41,   8.7,
            7.57, 7.82,  9.13, 9.79, 10.82, 9.63, 8.17,  8.86, 10.01,  8.02,
            7.87, 9.29,  9.34, 9.22,  8.65,  9.6, 9.82,  8.42,  8.25, 10.01,
            9.47,  9.4,  9.95, 7.54,  8.73, 9.75, 9.93,  9.48,   8.5,  10.7,
            9.47, 7.48,   8.4, 9.49,  8.81, 8.12, 7.91,  7.54,  8.96,  8.33,
            8.65, 9.75,  9.68, 8.55,  9.57, 8.53, 8.82, 10.95,  8.79,  9.13,
             9.8, 8.11,  8.95, 9.14,  9.83, 8.81,  8.2,   9.5,  7.66,  9.24
        ]

        n_sa_ci_width_maxt_fossilsim = [
            1.57622429004417, 1.22349799634928, 1.28763334702222, 1.67618115459824,  1.52891167683276, 1.51385819225667, 1.59681232915626, 1.12131353277657, 1.21379257043385, 1.24245662095282,
            1.32456860993466, 1.64605467759007, 1.29899868259122, 1.30339590642632, 0.848647741939487, 1.24595744791714, 1.19990399346587, 1.33610768996811, 1.24331831723811, 1.33733417686308,
            1.50146668819444, 1.21099042794138,  1.3965337777237, 1.44499497940479,  1.50186075885546, 1.15492259673298, 1.25377342305683,  1.4591914827857, 1.67983029445888, 1.29494915870715,
            1.44201647533983, 1.24429948650685, 1.38312660010542, 1.41025790663705,  1.54699596233305, 1.52027707150421, 1.25642765414528, 1.50921671248781, 1.33140613388002, 1.33006623971147,
            1.16780643961405, 1.14469179902167, 1.27712728837304, 1.47918921771969,  1.57155465766234, 1.50296489794155, 1.17568814046769, 1.27023128666246, 1.34226216349697, 1.06736643410146,
            1.10230431259964, 1.64818438656639, 1.29119801484333, 1.80403906696931,  1.29514541908067,  1.3808803516071, 1.40770352480108, 1.14156882022642, 1.33759384453177, 1.75230129737088,
            1.44392041583763, 1.78966072062542, 1.46603505683827, 1.00344544546288,  1.31402711778054, 1.39746568595953, 1.62994555328289, 1.60480797605196, 1.25840887190487, 1.52751004890146,
            1.56836003442314, 1.27585813213434, 1.18028176742041, 1.45277247949282,  1.23701036733282, 1.29301792842649, 1.13813557907819, 1.29275981285673,  1.6672910268099, 1.32438110478331,
            1.46471102400384,  1.5246812947318, 1.64643771000585, 1.22080554627598,  1.46732616081845, 1.28191894731092, 1.33846532108889, 1.66832754819934, 1.51239514851527, 1.42067534431319,
            1.68029223047599, 1.24239415596979, 1.56784532065346, 1.43922729212027,  1.48148547400809,  1.3460727997483, 1.20146173598259, 1.44420122272389, 1.06890677202117, 1.49142238338325 
        ]

        root_ages_mean_maxt_fossilsim = [
            2.38366417139404, 2.26615642051174, 2.29884814593732, 2.34319135802387, 2.30531653051718, 2.20545251810196, 2.27447065726496, 2.32239234552556, 2.25607107827111, 2.32211751895091,
            2.37952066965928,  2.3173093272037, 2.37465571232764, 2.37183161929928, 2.30471754871234, 2.43816409639795, 2.36059092928566, 2.43390279468663, 2.38232248137719, 2.41715008765125,
            2.34640735126263, 2.25027910543559, 2.40834882562394, 2.41269512761795, 2.29759679613712, 2.37162834997766, 2.30683456339556,  2.3974272556093, 2.37777730683781, 2.31736769418913,
            2.33926033265122, 2.34683059987156, 2.29919055639335, 2.29952418743682, 2.39068376550115, 2.30738423115837,  2.4570738071157, 2.42941813035316, 2.34424512582289, 2.26253949740992,
            2.36782879841323,  2.3521675671179, 2.37370624566779,   2.380202603813, 2.44488733800071, 2.35086178444249, 2.44422217802661, 2.40403700842255, 2.44012025078346, 2.29854223695015,
            2.39087134015631,  2.3091933640457, 2.36665381711774,  2.2726937198092, 2.27661941118283, 2.31968496265258,  2.3155811217305, 2.32310267116237, 2.27482730441775, 2.33210792621038,
            2.39265551613513, 2.41343698564524, 2.38306926453957, 2.32035198852619, 2.32813429086543, 2.48211288264269, 2.35236607042751, 2.28842650494241, 2.25189889378323,  2.4066149647916,
            2.34775538181992, 2.27872815738236, 2.23724297390451, 2.32609428065826, 2.36578487901597,  2.3523509671141, 2.22649453744094, 2.29937247130154, 2.22490395826341, 2.27926207982628,
             2.3027658651316, 2.33315857575794, 2.30739468710185,  2.3470409763579, 2.35613153211857, 2.24839377618693,  2.3741445695071, 2.45614703318293, 2.37054526030235, 2.39808966833028,
            2.36546295861751, 2.27631512497324, 2.26860948389909, 2.39696502331146, 2.31670146868658, 2.31528162539359, 2.26796780515744, 2.47311972707381, 2.30022715489888, 2.33487333561942
        ]

        root_ages_ci_width_maxt_fossilsim = [
            0.120686995522846, 0.129156792405392,  0.12613483330819, 0.122613587372116, 0.122050275350315, 0.130636705421805, 0.116992493091365, 0.120010906523596, 0.147955699426646, 0.126754757452484,
            0.104255053599428, 0.116532093465973, 0.105508158848639,  0.10487191897401, 0.122595304390188, 0.106357697597291, 0.108591964947324, 0.117301788446831, 0.112587871953908, 0.112741571771597,
            0.114828590253191, 0.140856411965727, 0.106706364060663, 0.113488572167187, 0.135832198096519, 0.106245030796297, 0.125570490788327, 0.111768594202013,  0.11463782500689, 0.130322825881461,
            0.118413208578857, 0.104421097137126, 0.117824808904298, 0.105881655723213, 0.109074761228308, 0.128488631876412, 0.100781365411774, 0.102610763773376, 0.112463724421251, 0.110425321436047,
            0.110664541708584, 0.112015661649281, 0.131090960708247, 0.111976426672395, 0.103891220383613, 0.113850620461398, 0.109230124881383, 0.099802434333246, 0.115571383361228, 0.115674480948108,
            0.118394064964809, 0.110938265051153, 0.127783994062688, 0.122497744797376, 0.124482379026199,  0.12346457224317, 0.116363348780944, 0.103766063669536, 0.115449383873525,  0.12577007207312,
            0.117651917786888, 0.105079038445749, 0.115918517564665, 0.113523273421784, 0.126768929573245, 0.106076104074851, 0.114653778429165,  0.13538295957738, 0.142059337381011, 0.119141184048604,
             0.10483380103667, 0.127815794407722, 0.143563176589759, 0.121832791035764, 0.108104713897136, 0.111948597672012, 0.140092974729279,  0.11515903583842, 0.123081105582667, 0.133272715793624,
            0.121880013528989,   0.1108907552003, 0.115950895173066,  0.11334004114695, 0.112492213408179, 0.137926538786301, 0.114848133809852, 0.109357832522376, 0.113671108835167, 0.112710666320603,
             0.12396737086016, 0.127439129639021, 0.140587177544125,  0.11585982199213, 0.124690999666968, 0.113161658009068, 0.132245953801939, 0.101205188373861, 0.125043524955854, 0.113374220136956
        ]

        # parsing simulations
        n_sa_ci_overlap_count = 0
        root_age_ci_overlap_count = 0
        global_mean_sa = 0.0
        global_mean_root_age = 0.0
        for i, batch in enumerate(sim_batches):
            n_sas = [ann_tr.n_sa_nodes for ann_tr in batch]
            root_ages = [ann_tr.root_age for ann_tr in batch]

            mean_sa = statistics.mean(n_sas)
            mean_root_ages = statistics.mean(root_ages)

            global_mean_sa += mean_sa
            global_mean_root_age += mean_root_ages

            stdevs_n_sa = statistics.stdev(n_sas)
            stdevs_root_ages = statistics.stdev(root_ages)

            sterr_n_sa = stdevs_n_sa / math.sqrt(n_sim)
            sterr_root_ages = stdevs_root_ages / math.sqrt(n_sim)

            n_sa_ci_width_maxt = 1.96 * sterr_n_sa
            root_ages_ci_width_maxt = 1.96 * sterr_root_ages

            if abs(mean_sa - n_sa_mean_maxt_fossilsim[i]) \
                    <= (n_sa_ci_width_maxt + n_sa_ci_width_maxt_fossilsim[i]):
                n_sa_ci_overlap_count += 1

            if abs(mean_root_ages - root_ages_mean_maxt_fossilsim[i]) \
                    <= (root_ages_ci_width_maxt
                        + root_ages_ci_width_maxt_fossilsim[i]):
                root_age_ci_overlap_count += 1

        print("\n\nPJ global mean sa count = " + str(global_mean_sa / 100.0))
        print("FossilSim global mean sa count = " \
              + str(statistics.mean(n_sa_mean_maxt_fossilsim)))
        print("\nPJ global mean root age = " \
              + str(global_mean_root_age / 100.0))
        print("FossilSim global mean total root age = " \
              + str(statistics.mean(root_ages_mean_maxt_fossilsim)))

        # [==== * ====][.... + ....] if we take '+' to be the "truth" of the '*' interval, + cannot be more than '====' away from '*' 95% of the time
        # then abs('+' - '*') can be at most ('====' + '....'). '....' can be added because we still are guaranteed to see '+' falling within that range
        # 95% of the time
        print(("\n95% CIs of simulations here and from FossilSim "
               "overlapped ") + str(n_sa_ci_overlap_count) \
                + " times for the sampled ancestor count.")
        print(("\n95% CIs of simulations here and from FossilSim "
               "overlapped ") + str(root_age_ci_overlap_count) \
                + " times for root age.")

        exp_count = int(0.95 * n_batches)
        a_delta = math.ceil(0.07 * exp_count)
        self.assertAlmostEqual(
            n_sa_ci_overlap_count,
            exp_count,
            msg=("Mean absolute difference must be 1.96 * (stderr_python "
                 "+ stderr_fossilsim) apart ") + str(exp_count) \
                    + " (+/- " + str(a_delta) \
                    + ") out of 100 times.",
            delta=a_delta)
        self.assertAlmostEqual(
            root_age_ci_overlap_count,
            exp_count,
            msg=("Mean absolute difference must be 1.96 * (stderr_python "
                 "+ stderr_fossilsim) apart ") + str(exp_count) \
                    + " (+/- " + str(a_delta) \
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
    # $ python3 tests/distribution/test_dn_discrete_sse_fbd.py
    # 
    # or
    #
    # $ python3 -m tests.distribution.test_dn_discrete_sse_fbd
    #
    # or 
    #
    # $ python3 -m unittest tests.distribution.test_dn_discrete_sse_fbd.TestFBDTrees.test_tree_size_sa_count_max_t_fbd

    #############
    # Debugging #
    #############
    # total_n_states = 1

    # # not state-dependent (just state 0, and no transition)
    # rates_t0_s0 = [ sseobj.DiscreteStateDependentRate(name="lambda", val=1.0, event=sseobj.MacroevolEvent.W_SPECIATION, states=[0,0,0]),
    #                 sseobj.DiscreteStateDependentRate(name="mu", val=1.0, event=sseobj.MacroevolEvent.EXTINCTION, states=[0]),
    #                 sseobj.DiscreteStateDependentRate(name="psi", val=2.0, event=sseobj.MacroevolEvent.ANCESTOR_SAMPLING, states=[0]) ]

    # matrix_atomic_rate_params = [ rates_t0_s0 ] # 1D: time slices (i) , 2D: all rates from all states in i-th time slice

    # state_dep_par_manager = sseobj.DiscreteStateDependentParameterManager(matrix_atomic_rate_params, total_n_states)

    # event_handler = sseobj.MacroevolEventHandler(state_dep_par_manager)

    # # stop_condition = "size" # don't use this against FossilSim, PJ does not have GSA implemented
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

    # n_sa = float()
    # n_extant = float()
    # root_age = float()
    # origin_age = float()
    # for tr in trs:
    #     # print(tr.tree.as_string(schema="newick"))
    #     # print(tr.n_sa)
    #     # print(str(tr.root_age) + "\n")
    #     n_sa += tr.n_sa
    #     n_extant += tr.n_extant_terminal_nodes
    #     root_age += tr.root_age
    #     origin_age += tr.origin_age
    
    # print("mean sa = " + str(n_sa / n_sim))
    # print("mean extant count = " + str(n_extant / n_sim))
    # print("mean root age = " + str(root_age / n_sim))
    # print("mean origin age = " + str(origin_age / n_sim))

    unittest.main()