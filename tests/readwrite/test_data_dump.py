import unittest
import io

# pj imports
import phylojunction.calculation.discrete_sse as sseobj
import phylojunction.pgm.pgm as pgm
import phylojunction.readwrite.pj_write as pjwrite
import phylojunction.distribution.dn_discrete_sse as dnsse

__author__ = "Fabio K. Mendes"
__email__ = "fmendes@lsu.edu"


class TestDataDump(unittest.TestCase):

    @classmethod
    def setUpClass(cls) -> None:
        #########################
        # BiSSE pgm ingredients #
        #########################
        total_n_states = 2

        l0 = [ 1.0, 1.1, 0.9, 0.95, 1.05 ]
        l0rate = sseobj.DiscreteStateDependentRate(
            name="lambda0",
            val=l0,
            event=sseobj.MacroevolEvent.W_SPECIATION,
            states=[0,0,0]
        )

        mu0 = [ 0.23, 0.24, 0.25, 0.26, 0.27 ]
        mu0rate = sseobj.DiscreteStateDependentRate(
            name="mu0",
            val=mu0,
            event=sseobj.MacroevolEvent.EXTINCTION,
            states=[0]
        )
        
        q01rate = sseobj.DiscreteStateDependentRate(
            name="q01",
            val=0.75,
            event=sseobj.MacroevolEvent.ANAGENETIC_TRANSITION,
            states=[0,1]
        )
        
        rates_t0_s0 = [ l0rate, mu0rate, q01rate ]

        l1rate = sseobj.DiscreteStateDependentRate(
            name="lambda1",
            val=1.5,
            event=sseobj.MacroevolEvent.W_SPECIATION,
            states=[1,1,1]
        )

        mu1rate = sseobj.DiscreteStateDependentRate(
            name="mu1",
            val=0.25,
            event=sseobj.MacroevolEvent.EXTINCTION,
            states=[1]
        )

        q10rate = sseobj.DiscreteStateDependentRate(
            name="q10",
            val=0.75,
            event=sseobj.MacroevolEvent.ANAGENETIC_TRANSITION,
            states=[1,0]
        )

        rates_t0_s1 = [ l1rate, mu1rate, q10rate ]
        
        rates_t0 = rates_t0_s0 + rates_t0_s1

        # 1D: time slices (i)
        # 2D: all rates from all states in i-th time slice
        matrix_atomic_rate_params = [ rates_t0 ]

        state_dep_param_manager = sseobj.DiscreteStateDependentParameterManager(
            matrix_atomic_rate_params, total_n_states
        )

        meh = sseobj.MacroevolEventHandler(state_dep_param_manager)

        sse_stash = sseobj.SSEStash(meh)

        ########
        # Tree #
        ########
        n_sim = 5
        n_repl = 2
        stop_condition = "size"
        stop_condition_value = [50] # 50 living taxa
        start_at_origin = False
        start_states_list = [0 for i in range(n_sim)]

        sse_sim = dnsse.DnSSE(
            sse_stash,
            n=n_sim,
            n_replicates=n_repl,
            origin=start_at_origin,
            start_states_list=start_states_list,
            stop=stop_condition,
            stop_value=stop_condition_value,
            condition_on_speciation=True,
            condition_on_survival=True,
            epsilon=1e-12,
            runtime_limit=3600,
            debug=False
        )

        trs = sse_sim.generate()

        ####################
        # Initializing PGM #
        ####################
        cls.bisse_pgm = pgm.DirectedAcyclicGraph()

        # rv
        cls.bisse_pgm.add_node(
            pgm.StochasticNodeDAG(
                "l0", n_sim, value=l0, sampled_from="Log-normal"
            )
        )
        cls.bisse_pgm.add_node(
            pgm.StochasticNodeDAG(
                "mu0", n_sim, value=mu0, sampled_from="Log-normal"
            )
        )

        # deterministic
        cls.bisse_pgm.add_node(pgm.DeterministicNodeDAG("l0r", value=l0rate))
        cls.bisse_pgm.add_node(pgm.DeterministicNodeDAG("mu0r", value=mu0rate))
        cls.bisse_pgm.add_node(pgm.DeterministicNodeDAG("q01r", value=q01rate))
        cls.bisse_pgm.add_node(pgm.DeterministicNodeDAG("l1r", value=l1rate))
        cls.bisse_pgm.add_node(pgm.DeterministicNodeDAG("mu1r", value=mu1rate))
        cls.bisse_pgm.add_node(pgm.DeterministicNodeDAG("q10r", value=q10rate))
        cls.bisse_pgm.add_node(pgm.DeterministicNodeDAG("meh", value=meh))

        # more rv
        cls.bisse_pgm.add_node(
            pgm.StochasticNodeDAG(
                "trs", n_sim, value=trs, sampled_from="DnSSE",
                replicate_size=n_repl
            )
        )

        # sorted_node_dag_list = bisse_pgm.get_sorted_node_dag_list()

        ###################
        # Output handling #
        ###################
        # populating dataframe to be dumped
        # data_df_names_list, data_df_list = \
        #     pjwrite.prep_data_df(sorted_node_dag_list)
        cls.scalar_output_stash, cls.tree_output_stash = \
            pjwrite.prep_data_df(cls.bisse_pgm, write_nex_states=True)

        cls.output_fp_list, cls.output_df_str_list = \
            pjwrite.prep_data_filepaths_dfs(
                cls.scalar_output_stash,
                cls.tree_output_stash
            )
        
        cls.scalar_constant_df, \
        cls.scalar_value_df_dict, \
        cls.scalar_repl_summary_df = cls.scalar_output_stash

        cls.tree_value_df_dict, cls.tree_ann_value_df_dict, \
        cls.tree_rec_value_df_dict, cls.tree_rec_ann_value_df_dict, \
        cls.tree_summary_df_dict, cls.tree_repl_summary_df_dict, \
        cls.tree_living_nd_states_str_dict, \
        cls.tree_living_nd_states_str_nexus_dict, \
        cls.tree_internal_nd_states_str_dict = cls.tree_output_stash 

    def test_bisse_data_dump_file_names(self):
        """Test output file names."""

        self.assertEqual(
            [
                "scalar_rvs_repl1.csv",
                "trs_complete.tsv",
                "trs_annotated_complete.tsv",
                "trs_reconstructed.tsv",
                "trs_annotated_reconstructed.tsv",
                "trs_stats.csv",
                "trs_stats_summary.csv",
                "trs_sample1_repl1.tsv", "trs_sample1_repl2.tsv",
                "trs_sample2_repl1.tsv", "trs_sample2_repl2.tsv",
                "trs_sample3_repl1.tsv", "trs_sample3_repl2.tsv",
                "trs_sample4_repl1.tsv", "trs_sample4_repl2.tsv",
                "trs_sample5_repl1.tsv", "trs_sample5_repl2.tsv",
                "trs_sample1_repl1.nex", "trs_sample1_repl2.nex",
                "trs_sample2_repl1.nex", "trs_sample2_repl2.nex",
                "trs_sample3_repl1.nex", "trs_sample3_repl2.nex",
                "trs_sample4_repl1.nex", "trs_sample4_repl2.nex",
                "trs_sample5_repl1.nex", "trs_sample5_repl2.nex",
                "trs_sample1_repl1_anc_states.tsv",
                "trs_sample1_repl2_anc_states.tsv",
                "trs_sample2_repl1_anc_states.tsv",
                "trs_sample2_repl2_anc_states.tsv",
                "trs_sample3_repl1_anc_states.tsv",
                "trs_sample3_repl2_anc_states.tsv",
                "trs_sample4_repl1_anc_states.tsv",
                "trs_sample4_repl2_anc_states.tsv",
                "trs_sample5_repl1_anc_states.tsv",
                "trs_sample5_repl2_anc_states.tsv"
            ],
            self.output_fp_list
        )

    def test_bisse_data_dump_scalar_val_table(self):
        """Test scalar variables values."""

        # writing to file handle
        scalar_constant_df_outfile = io.StringIO()
        scalar_value_df_dict_outfile = io.StringIO()

        pjwrite.write_data_df(scalar_constant_df_outfile, self.scalar_constant_df)
        scalar_constant_df_outfile.seek(0)  # move to the start of the file handle
        csvstring_dump_scalar_constant_df = scalar_constant_df_outfile.read()

        # this dataframe is empty
        self.assertEqual(csvstring_dump_scalar_constant_df, "\n")

        pjwrite.write_data_df(
            scalar_value_df_dict_outfile,
            self.scalar_value_df_dict[1]
        )  # [1] accesses replicate "1"
        
        scalar_value_df_dict_outfile.seek(0)  # move to the start of the file handle
        csvstring_dump_scalar_value_df_dict = scalar_value_df_dict_outfile.read()

        # debugging
        # print("csvstring_dump_scalar_value_df_dict = ")
        # print(csvstring_dump_scalar_value_df_dict)
        
        self.assertEqual(
            csvstring_dump_scalar_value_df_dict,
            "sample,replicate,l0,mu0\n1,1,1.0,0.23\n2,1,1.1,0.24\n" \
            + "3,1,0.9,0.25\n4,1,0.95,0.26\n5,1,1.05,0.27\n"
        )

    def test_bisse_data_dump_scalar_summary_table(self):
        """Test scalar variables summary."""

        # writing to file handle
        scalar_repl_summary_df_outfile = io.StringIO()
        scalar_value_df_dict_outfile = io.StringIO()

        pjwrite.write_data_df(
            scalar_repl_summary_df_outfile,
            self.scalar_repl_summary_df
        )

        # move to the start of the file handle
        scalar_value_df_dict_outfile.seek(0)
        scalar_repl_summary_df_outfile.seek(0)

        csvstring_dump_scalar_repl_summary_df_outfile = \
            scalar_value_df_dict_outfile.read()
        
        # this dataframe is empty
        self.assertEqual(csvstring_dump_scalar_repl_summary_df_outfile, "")

    def test_bisse_data_dump_tree_val_table(self):
        """Test tree values."""

        # writing to file handle
        tree_value_df_dict_outfile = io.StringIO()
       
        # tree_value_df_dict:
        #     key: tree node name
        #     value: dataframe
        
        pjwrite.write_data_df(
            tree_value_df_dict_outfile,
            self.tree_value_df_dict["trs"]
        )

        tree_value_df_dict_outfile.seek(0)  # move to the start of the file handle
        csvstring_dump_tree_value_df_dict = \
            tree_value_df_dict_outfile.read()
        csvstring_dump_tree_value_df_dict_list = \
            csvstring_dump_tree_value_df_dict.split("\n")

        # debugging
        # print(csvstring_dump_tree_value_df_dict_list)

        self.assertEqual(
            csvstring_dump_tree_value_df_dict_list[0],
            "sample,replicate,trs"
        )

        self.assertEqual(
            csvstring_dump_tree_value_df_dict_list[1][:4],
            "1,1,"
        )

        self.assertEqual(
            csvstring_dump_tree_value_df_dict_list[9][:4],
            "5,1,"
        )

        self.assertEqual(
            csvstring_dump_tree_value_df_dict_list[1][-2],
            ";"
        )

        self.assertEqual(
            csvstring_dump_tree_value_df_dict_list[9][-2],
            ";"
        )
       
    def test_bisse_data_dump_rec_tree_val_table(self):
        """Test rec. tree values."""

        # writing to file handle
        tree_rec_value_df_dict_outfile = io.StringIO()

        pjwrite.write_data_df(
            tree_rec_value_df_dict_outfile,
            self.tree_rec_value_df_dict["trs"]
        )

        tree_rec_value_df_dict_outfile.seek(0)  # move to the start of the file handle
        
        csvstring_dump_tree_rec_value_df_dict = \
            tree_rec_value_df_dict_outfile.read()
        csvstring_dump_tree_rec_value_df_dict_list = \
            csvstring_dump_tree_rec_value_df_dict.split("\n")
        
        self.assertEqual(
            csvstring_dump_tree_rec_value_df_dict_list[0],
            "sample,replicate,trs"
        )

        self.assertEqual(
            csvstring_dump_tree_rec_value_df_dict_list[1][:4],
            "1,1,"
        )

        self.assertEqual(
            csvstring_dump_tree_rec_value_df_dict_list[9][:4],
            "5,1,"
        )

        self.assertEqual(
            csvstring_dump_tree_rec_value_df_dict_list[1][-2],
            ";"
        )

        self.assertEqual(
            csvstring_dump_tree_rec_value_df_dict_list[9][-2],
            ";"
        )

    def test_bisse_data_dump_ann_tree_val_table(self):
        """Test annotated tree values."""

        # writing to file handle
        tree_ann_value_df_dict_outfile = io.StringIO()
        
        pjwrite.write_data_df(
            tree_ann_value_df_dict_outfile,
            self.tree_ann_value_df_dict["trs"]
        )

        tree_ann_value_df_dict_outfile.seek(0)  # move to the start of the file handle
        
        csvstring_dump_tree_ann_value_df_dict = \
            tree_ann_value_df_dict_outfile.read()
        csvstring_dump_tree_ann_value_df_dict_list = \
            csvstring_dump_tree_ann_value_df_dict.split("\n")
        
        self.assertEqual(
            csvstring_dump_tree_ann_value_df_dict_list[0],
            "sample,replicate,trs"
        )

        self.assertIn(
            "&state=",
            csvstring_dump_tree_ann_value_df_dict_list[1]
        )    

    def test_bisse_data_dump_rec_ann_tree_val_table(self):
        """Test rec. ann. tree values."""

        # writing to file handle
        tree_rec_ann_value_df_dict_outfile = io.StringIO()

        pjwrite.write_data_df(
            tree_rec_ann_value_df_dict_outfile,
            self.tree_rec_ann_value_df_dict["trs"]
        )
        
        tree_rec_ann_value_df_dict_outfile.seek(0)  # move to the start of the file handle

        csvstring_dump_tree_rec_ann_value_df_dict = \
            tree_rec_ann_value_df_dict_outfile.read()
        csvstring_dump_tree_rec_ann_value_df_dict_list = \
            csvstring_dump_tree_rec_ann_value_df_dict.split("\n")

        self.assertEqual(
            csvstring_dump_tree_rec_ann_value_df_dict_list[0],
            "sample,replicate,trs"
        )

        self.assertIn(
            "&state=",
            csvstring_dump_tree_rec_ann_value_df_dict_list[1]
        )

        # debugging
        # print(csvstring_dump_tree_rec_ann_value_df_dict_list[0])

    def test_bisse_data_dump_tree_summary_table(self):
        """Test tree summary."""

        # writing to file handle
        tree_summary_df_dict_outfile = io.StringIO()

        pjwrite.write_data_df(
            tree_summary_df_dict_outfile,
            self.tree_summary_df_dict["trs"]
        )

        tree_summary_df_dict_outfile.seek(0)  # move to the start of the file handle
        
        csvstring_dump_tree_summary_df_dict = \
            tree_summary_df_dict_outfile.read()
        csvstring_dump_tree_summary_df_dict_list = \
            csvstring_dump_tree_summary_df_dict.split("\n")

        self.assertEqual(
            csvstring_dump_tree_summary_df_dict_list[0],
            "sample,replicate,origin_age,root_age,n_total,n_extant,n_extinct,n_sa,n_0,n_1"
        )

        self.assertEqual(
            csvstring_dump_tree_summary_df_dict_list[1][:4],
            "1,1,"
        )

    def test_bisse_data_dump_tree_repl_summary_table(self):
        """Test tree replicate summary."""
        
        # writing to file handle
        tree_repl_summary_df_dict_outfile = io.StringIO()

        pjwrite.write_data_df(
            tree_repl_summary_df_dict_outfile,
            self.tree_repl_summary_df_dict["trs"]
        )

        tree_repl_summary_df_dict_outfile.seek(0)  # move to the start of the file handle
        
        csvstring_dump_tree_repl_summary_df_dict = \
            tree_repl_summary_df_dict_outfile.read()
        csvstring_dump_tree_repl_summary_df_dict_list = \
            csvstring_dump_tree_repl_summary_df_dict.split("\n")
        
        self.assertEqual(
            csvstring_dump_tree_repl_summary_df_dict_list[0],
            "sample,summary,origin_age,root_age,n_total,n_extant,n_extinct,n_sa"
        )

        self.assertEqual(
            csvstring_dump_tree_repl_summary_df_dict_list[1][:10],
            "1,average,"
        )

    def test_bisse_data_dump_tree_nd_states_table(self):
        """Test tree states summary."""

        tree_living_nd_states_str_dict_outfile = \
            io.StringIO(self
                        .tree_living_nd_states_str_dict[
                            "trs_sample1_repl1.tsv"
                        ])
        
        tree_living_nd_states_str_nexus_dict_outfile = \
            io.StringIO(self
                        .tree_living_nd_states_str_nexus_dict[
                            "trs_sample1_repl1.nex"])

        tree_internal_nd_states_str_dict_outfile = \
            io.StringIO(self
                        .tree_internal_nd_states_str_dict[
                            "trs_sample5_repl2_anc_states.tsv"])
        
        # move to the start of the file handle
        tree_living_nd_states_str_dict_outfile.seek(0)
        tree_living_nd_states_str_nexus_dict_outfile.seek(0)
        tree_internal_nd_states_str_dict_outfile.seek(0)
        
        csvstring_dump_tree_living_nd_states_str_dict = \
            tree_living_nd_states_str_dict_outfile.read()
        csvstring_dump_tree_living_nd_states_str_dict_list = \
            csvstring_dump_tree_living_nd_states_str_dict.split("\n")

        # debugging
        # print(csvstring_dump_tree_living_nd_states_str_dict_list)

        self.assertIn(
            "nd",
            csvstring_dump_tree_living_nd_states_str_dict_list[0]
        )

        csvstring_dump_tree_living_nd_states_str_nexus_dict = \
            tree_living_nd_states_str_nexus_dict_outfile.read()
        csvstring_dump_tree_living_nd_states_str_nexus_dict_list = \
            csvstring_dump_tree_living_nd_states_str_nexus_dict.split("\n")

        # debugging
        # print(csvstring_dump_tree_living_nd_states_str_nexus_dict_list[3])

        self.assertEqual(
            csvstring_dump_tree_living_nd_states_str_nexus_dict_list[3],
            "Dimensions ntax=50 nchar=1;"
        )

        csvstring_dump_tree_internal_nd_states_str_dict = \
            tree_internal_nd_states_str_dict_outfile.read()
        csvstring_dump_tree_internal_nd_states_str_dict_list = \
            csvstring_dump_tree_internal_nd_states_str_dict.split("\n")

        trs_node = self.bisse_pgm.get_node_dag_by_name("trs")
        last_rec_tr = trs_node.value[9].tree_reconstructed
        last_rec_tr_seed_node_label = last_rec_tr.seed_node.label

        # debugging
        # print(last_rec_tr)
        # print("seed node of last_rec_tr = " + last_rec_tr.seed_node.label)

        self.assertIn(
            last_rec_tr_seed_node_label + "\t",
            csvstring_dump_tree_internal_nd_states_str_dict_list[0]
        )


if __name__ == "__main__":
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
    # $ python3.11 tests/readwrite/test_data_dump.py
    # 
    # or
    #
    # $ python3.11 -m tests.readwrite.test_data_dump
    #
    # or 
    #
    # $ python3.11 -m unittest tests.readwrite.test_data_dump.TestDataDump.test_bisse_data_dump_file_names

    unittest.main()