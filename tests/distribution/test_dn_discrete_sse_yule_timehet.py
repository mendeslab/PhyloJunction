import unittest
import math
import statistics

# pj imports
import phylojunction.utility.helper_functions as pjh
import phylojunction.calculation.discrete_sse as sseobj
import phylojunction.distribution.dn_discrete_sse as distsse

__author__ = "Fabio K. Mendes"
__email__ = "fmendes@lsu.edu"


class TestYuleTimeHetTrees(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.seed_age_for_time_slicing = 3.0

        # (time units = t.u., with the 3.0 coming from the seed age)
        #
        # 3.0 -> 2.2 t.u. ago: lambda = 1.0
        # 2.2 -> 1.2 t.u. ago: lambda = 0.25
        # 1.2 -> 0.7 t.u. ago: lambda = 3.0
        # 0.7 -> 0.0 (present) t.u. ago: lambda = 0.4
        #
        # in MASTER, we do it forward in time: 1.0:0.0,0.25:0.8,3.0:1.8,0.4:2.3
        # (there is an implicit 3.0 boundary for the last interval,
        # and it comes, again, from the seed age)
        #
        # which we read as
        # lambda = 1.0 from 0.0 -> 0.8 fwd-t.u.
        # lambda = 0.25 from 0.8 -> 1.8 fwd-t.u.
        # lambda = 3.0 from 1.8 -> 2.3 fwd-t.u.
        # lambda = 0.4 from 2.3 -> 3.0 fwd-t.u
        cls.time_slice_age_ends = [2.2, 1.2, 0.7] # if 2 slices, 1 age end; 0.0 is assumed to be the last age end

        # not state-dependent (just state 0, and no transition)
        total_n_states = 1

        rates_t0_s0 = [
            sseobj.DiscreteStateDependentRate(
                name="lambda_t0",
                val=1.0,
                event=sseobj.MacroevolEvent.W_SPECIATION,
                states=[0,0,0],
                epoch_idx=1)]
        rates_t1_s0 = [
            sseobj.DiscreteStateDependentRate(
                name="lambda_t1",
                val=0.25,
                event=sseobj.MacroevolEvent.W_SPECIATION,
                states=[0,0,0],
                epoch_idx=2)]
        rates_t2_s0 = [
            sseobj.DiscreteStateDependentRate(
                name="lambda_t2",
                val=3.0,
                event=sseobj.MacroevolEvent.W_SPECIATION,
                states=[0,0,0],
                epoch_idx=3)]
        rates_t3_s0 = [
            sseobj.DiscreteStateDependentRate(
                name="lambda_t3",
                val=0.4,
                event=sseobj.MacroevolEvent.W_SPECIATION,
                states=[0,0,0],
                epoch_idx=4)]

        matrix_atomic_rate_params = [rates_t0_s0, rates_t1_s0, rates_t2_s0, rates_t3_s0]

        state_dep_par_manager = sseobj.DiscreteStateDependentParameterManager(
            matrix_atomic_rate_params,
            total_n_states,
            seed_age_for_time_slicing=cls.seed_age_for_time_slicing,
            list_time_slice_age_ends=cls.time_slice_age_ends
        )

        event_handler = sseobj.MacroevolEventHandler(state_dep_par_manager)

        cls.sse_stash = sseobj.SSEStash(event_handler)


    def test_tree_size_state_count_max_taxa_timehet_yule(self):
        """
        Test if time-het Yule trees simulated here have similar root ages and number of tips for both states
        as time-het Yule trees simulated with MASTER
        """

        # setting up stopping conditions
        runtime_limit = 3600 # 1h

        stop_condition = "age"
        stop_condition_value = [3.0] # 3.0 time units

        start_at_origin = True

        # simulation initialization
        n_batches = i = 100
        n_sim = 100
        start_states_list = [0 for i in range(n_sim)]
        # seeds_list = [i+1 for i in range(n_sim)]

        # simulations
        print(("Running TestYuleTimeHetTrees.test_tree_size_state_count_max_"
               "taxa_timehet_yule"))
        sim_batches = list()
        for i in range(n_batches):
            # print("Doing batch " + str(n_batches - i))
            sse_sim = distsse.DnSSE(self.sse_stash,
                                    n=n_sim,
                                    origin=start_at_origin,
                                    start_states_list=start_states_list,
                                    stop=stop_condition,
                                    stop_value=stop_condition_value,
                                    condition_on_speciation=True,
                                    condition_on_survival=True,
                                    epsilon=1e-12,
                                    runtime_limit=runtime_limit,
                                    debug=False)

            trs = sse_sim.generate()

            sim_batches.append(trs)

            # printing progress
            pjh.print_progress(i , n_batches)


        # "expectations" from MASTER
        n_mean_master = [
            18.83, 16.68, 16.81, 14.21, 20.18, 14.85,  16.6, 16.36, 17.25, 16.33,
            13.77, 15.77, 15.34, 15.65, 16.86,  17.5, 18.71, 21.19, 16.23, 16.33,
            19.26, 17.25, 17.35, 15.41, 17.56, 18.18, 16.39,  16.6, 15.25, 17.41,
            19.05, 17.32, 16.57, 16.27, 17.03, 18.29, 16.93, 18.88, 17.44, 16.87,
            17.11, 14.17, 18.02, 16.99, 15.58, 16.79, 15.58, 19.93, 16.26, 15.77,
            14.86, 18.29,  14.2, 15.47, 16.53, 17.01, 16.46, 16.27, 16.31, 17.15,
            18.89, 19.31, 16.58, 18.39, 15.77, 14.97, 14.75, 19.76, 18.46, 18.51,
            14.68, 15.05,  18.6, 17.78, 16.41, 16.12, 18.31, 19.62, 19.77,  17.5,
            13.45, 17.57, 15.42, 17.43, 16.49,  15.4, 16.07,  16.4, 16.63,  17.1,
            16.33, 17.16, 20.43, 19.53, 17.17, 17.93, 18.09, 16.19,  17.6, 16.16
        ]

        n_ci_width_maxt_master = [
            3.96801841210293, 3.21866580703461, 2.88957033302444, 2.17138661986513, 3.59007923178654,  2.9409072674186, 3.02271390533993, 3.32745734309262, 3.6383926909746,  3.17486696096041,
            2.26866819686678,  2.8674838499805, 3.91715003469914,  2.9619433628698, 3.17962579692719, 3.77080422791379, 3.89629896577599, 3.86383615082266, 3.03486647129422, 4.00665374091597,
            3.57689436401372, 3.11830201788529, 3.32755704972761, 2.83881962068531, 3.50995086765233, 3.52925145032144, 3.23562182779837, 3.16764088363513, 3.17135968036144, 3.28881619278107,
            4.17319354495288, 2.98107500805269, 3.66695704657338, 2.63902190972337, 3.02905602950804, 3.48342538002206,  3.3217072039978, 3.62417009791012, 3.61742268070917, 3.21863868105631,
            3.10677252492296,  2.4299105893681,  3.3848582611508, 3.52277110826294, 3.28422148541363, 3.53561345732193, 2.70958744811639, 3.87692596146832, 3.78365798612396,  2.7913686215572,
            2.95665177926627, 3.90495386157252,  3.3343865487621, 3.02469728773381, 3.71825989300361,  2.8599824818281, 3.30251261161846, 2.94004685484549, 2.84722732084491, 3.35635277639225,
            3.25916193701957, 3.51981779578937, 3.36372142062263, 3.55974239880621,  2.4999332422602, 2.79325857743318, 2.50296130875635, 3.30951984042118, 3.06167110479545, 3.43860077227804,
            2.63357951315718, 2.27217529377831, 3.22927202604268, 3.18859039588456, 3.07554005283067, 2.82376902824018, 3.43467140480569, 3.97463585646725, 3.72017963528482, 3.03149482745541,
            2.64614664763843, 3.79355563091521,  2.8098415889782, 3.71144933459051, 3.30357873876684, 2.57142149109174, 3.61140358869014, 3.71058615533471, 3.61799492012084,  4.0561712354028,
            2.92326332413319, 3.25127350191312, 4.15370320858996, 3.72618286901399, 3.40395786535188,  2.6550592988489, 3.23492617285749, 3.26663043616272,  3.5389817858825, 3.01277676383214
        ]

        root_ages_mean_maxt_master = [
            1.97629921227542, 1.85051399695374, 2.06264636332191, 1.82110224062833, 1.92549484066431, 1.85971869734166, 1.85068873761823,  1.8774340764665, 1.86874749572255, 1.96802999405385,
             1.9554688468196, 1.93509987348782, 1.76558501337048, 1.85155996729765, 1.93061553574989, 1.85518944251469, 1.90602975221464, 2.10411816085201, 1.71192190185486, 1.75396687591692,
            2.02581535383311, 1.95287658616189, 1.94509642496361, 1.90167514805919, 1.95631995034473, 1.87520939540817, 1.78875152768178, 1.80106018727001, 1.82006619612074, 1.88814676673972,
            1.89305364769844, 1.92906429250561, 1.74007228918927, 1.83654701027872, 2.02824927930184,  1.9082733984022, 1.85190867163633, 1.88862272783124, 1.98059593009719, 1.92702280584655,
            1.91069542019431, 1.80267147674487, 1.97899874481849, 1.77758229128923, 1.76096220498027, 1.88003060321508, 1.96061908532571, 1.99309670005335, 1.74188040078966, 1.85788239277638,
            1.75480942928007, 2.00326210348625,  1.7596798342906, 1.78320697908013, 1.85631853752761, 1.93054577432492, 1.89082595160985, 1.91863994092888, 1.94759215630421, 1.95304511308185,
            1.97315113732319, 2.01784206767499, 1.84558889540385, 1.88070998861862, 1.83327537511832,  1.9679431728126, 1.81992943624496, 2.06192169759929,   1.973801931514, 1.97067716168212,
            2.01371723590042, 1.98784300900212, 1.97701923875934, 1.98791768774367, 1.89352334403068, 1.92302779522217, 2.09894328854021, 1.97861459551644, 1.93049298276508, 2.01800223047691,
            1.66398364207401, 1.93049953503647, 1.84508304271114, 2.00407254461266, 1.79776414252716, 1.96066358520783, 1.88929660975402, 1.83221379188934,  1.8942664262065, 1.96388242283886,
            1.69958409917955, 2.02914786937395, 1.94903494451456, 1.97041156511769, 1.99116796872001,   2.083744480241, 1.99747021932381, 1.89962843477574, 1.85273467363469, 1.88303096021827
        ]

        root_ages_ci_width_maxt_master = [
            0.193115987862331, 0.174152327359508, 0.173308445598097, 0.179322183620376,  0.19239906750395, 0.186469681259913,    0.178453407341, 0.182823343971814, 0.182943342722912,  0.15780269765017,
            0.173818929441031, 0.171837524324625, 0.197054328205507, 0.193388219944991,   0.1820117488964, 0.180742384015258, 0.172125054807962, 0.155005053486548, 0.202612512638381, 0.170862901208229,
            0.172240759094807, 0.175162482043405, 0.180249084477777,  0.17802368898265, 0.185132447930469, 0.187342784125557, 0.173590530711008, 0.183188068305647, 0.180641014285398,  0.18599117964152,
            0.196358346189945, 0.171864021003998,  0.17988864330182, 0.184099531521293, 0.177232820746839, 0.181066027676122, 0.185453557520067, 0.179267880046955, 0.157094417169529, 0.181888968156683,
            0.190035074367999, 0.167875959663171, 0.181521941527234, 0.198378137206061, 0.190860020924085, 0.166736545494386, 0.178762489022622, 0.174341485150705, 0.192988701379276, 0.180676911707256,
            0.194453872884801, 0.165671506197077, 0.194674760353365, 0.184157017945081, 0.175946556729693,  0.18380004359209, 0.179967637979974, 0.190630531227159, 0.172613393501352, 0.187062497758941,
             0.17590233529447, 0.175738411675054, 0.179340657872063, 0.180011911879459, 0.172194969854192, 0.177233295000952,  0.18962206776576, 0.168989733375307, 0.172189347689699, 0.182955457750696,
            0.177192454693026,  0.16184412542311,  0.17946791200871, 0.168318200249559, 0.179270650624793, 0.179716083537765, 0.166950374654242, 0.184646858384092, 0.180184310389181, 0.166432947126615,
            0.196508293275151, 0.175292949028916, 0.183120250138278, 0.166917676788866, 0.177969834005972, 0.180554322401007, 0.180269181019248, 0.185841222113395, 0.183256946983878, 0.186319365595564,
            0.184660256878766, 0.174707888037051, 0.188163529695133, 0.172267423984443, 0.164119688371032, 0.166581290895668, 0.166569686809021, 0.172029278852548, 0.190093925860537, 0.169940952610519
        ]

        # parsing simulations
        n_ci_overlap_count = 0
        root_age_ci_overlap_count = 0
        for i, batch in enumerate(sim_batches):
            ns = [ann_tr.n_extant_terminal_nodes for ann_tr in batch]
            root_ages = [ann_tr.root_age for ann_tr in batch]

            mean_n = statistics.mean(ns)
            mean_root_ages = statistics.mean(root_ages)

            # print("mean_n = " + str(mean_n))
            # print("mean_n MASTER = " + str(n_mean_master[i]))
            # print("mean_root_ages = " + str(mean_root_ages))
            # print("mean_root_ages MASTER = " + str(root_ages_mean_maxt_master[i]))

            stdevs_n = statistics.stdev(ns)
            stdevs_root_ages = statistics.stdev(root_ages)

            sterr_ns = stdevs_n / math.sqrt(n_sim)
            sterr_root_ages = stdevs_root_ages / math.sqrt(n_sim)

            n_ci_width_maxt = 1.96 * sterr_ns
            root_ages_ci_width_maxt = 1.96 * sterr_root_ages

            if abs(mean_n - n_mean_master[i]) <= (n_ci_width_maxt + n_ci_width_maxt_master[i]):
                n_ci_overlap_count += 1

            if abs(mean_root_ages - root_ages_mean_maxt_master[i]) <= (root_ages_ci_width_maxt + root_ages_ci_width_maxt_master[i]):
                root_age_ci_overlap_count += 1

        # [==== * ====][.... + ....] if we take '+' to be the "truth" of the '*' interval, + cannot be more than '====' away from '*' 95% of the time
        # then abs('+' - '*') can be at most ('====' + '....'). '....' can be added because we still are guaranteed to see '+' falling within that range
        # 95% of the time

        print("\n95% CIs of simulations here and from MASTER overlapped " + str(n_ci_overlap_count) + " times for the number of extant taxa.")
        print("\n95% CIs of simulations here and from MASTER overlapped " + str(root_age_ci_overlap_count) + " times for root age.")
        exp_count = int(0.95 * n_batches)
        a_delta = math.ceil(0.07 * exp_count)
        self.assertAlmostEqual(n_ci_overlap_count, exp_count,
                                msg="Mean absolute difference must be 1.96 * (stderr_python + stderr_divtree) apart " + str(exp_count) + " (+/- " + str(a_delta) + ") out of 100 times.", delta=a_delta)
        self.assertAlmostEqual(root_age_ci_overlap_count, exp_count,
                                msg="Mean absolute difference must be 1.96 * (stderr_python + stderr_divtree) apart " + str(exp_count) + " (+/- " + str(a_delta) + ") out of 100 times.", delta=a_delta)

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
    # on the terminal, remember to add "src/phylojunction" to
    # PYTHONPATH (system variable), or to set it if it does not
    # exist -- don't forget to export it!
    # 
    # Then you can do:
    # $ python3 tests/distribution/test_dn_discrete_sse_yule_timehet.py
    # 
    # or
    #
    # $ python3 -m tests.distribution.test_dn_discrete_sse_yule_timehet
    #
    # or 
    #
    # $ python3 -m unittest tests.distribution.test_dn_discrete_sse_yule_timehet.TestYuleTimeHetTrees.test_tree_size_state_count_max_taxa_timehet_yule

    unittest.main()