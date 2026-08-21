import re
import typing as ty
import numpy as np
import matplotlib.pyplot as plt  # type: ignore
import random

# pj imports
import phylojunction.pgm.pgm as pgm
import phylojunction.data.tree as pjtr
# import user_interface.phylojunction_inference as pjinf
# import user_interface.phylojunction_io as pjio
import phylojunction.interface.cmdbox.cmd_parse_utils as cmdu
import phylojunction.interface.grammar.ct_fn_grammar as ctgrammar
import phylojunction.interface.grammar.dn_grammar as dngrammar
import phylojunction.interface.grammar.det_fn_grammar as detgrammar
import phylojunction.utility.exception_classes as ec

__author__ = "Fabio K. Mendes"
__email__ = "fmendes@lsu.edu"


def script2dag(script_file_path_or_model_spec: str,
               in_pj_file: bool = True,
               random_seed: ty.Optional[int] = None) -> pgm.DirectedAcyclicGraph:
    """Go through .pj lines and populate and return DAG.

    Args:
        script_file_path_or_model_spec (str): Path to script .pj
            file or full string specifying model directly.
        random_seed (int): Random seed (integer) for simulation.
            Defaults to None.

    Returns:
        (DirectedAcyclicGraph): DAG object (model) built from .pj
            script commands.
    """

    def _execute_spec_lines(
            all_lines_list: ty.List[str],
            dag_obj: pgm.DirectedAcyclicGraph) -> None:

        for line in all_lines_list:
            # clear padding whitespaces
            line = line.lstrip().rstrip()
            _ = cmdline2dag(dag_obj, line)

    # initialize model DAG
    dag = pgm.DirectedAcyclicGraph()
    if random_seed is not None:
        dag.random_seed = random_seed

    # get script strings
    all_lines_list: ty.List[str] = list()
    if in_pj_file:
        with open(script_file_path_or_model_spec, "r") as infile:
            all_lines_list = infile.readlines()

    else:
        all_lines_list = script_file_path_or_model_spec.split("\n")

    # side-effect: populates DAG
    _execute_spec_lines(all_lines_list, dag)

    return dag


def cmdline2dag(dag_obj: pgm.DirectedAcyclicGraph,
                cmd_line: str) -> str:
    """Update PGM object with user-entered .pj command.

    Args:
        dag_obj (DirectedAcyclicGraph): DAG object that will hold
            all stochastic nodes (random variables) created by user via
            commands.
        cmd_line (str): Command line string provided by user through
            static script or via GUI.

    Returns:
        (str): Command string if there was nothing wrong with it.
    """

    # skip lines without useful commands
    if cmd_line.startswith("#") or \
            cmd_line.startswith("\n") or not \
            cmd_line:
        return ""

    ##################################
    # Checking command line is valid #
    ##################################
    if re.match(cmdu.cannot_start_with_this_regex, cmd_line):
        raise ec.ScriptSyntaxError(
            cmd_line,
            ("A special character or digit was detected as the first "
             "character in this line, but this is not allowed. Exiting..."))

    if re.search(cmdu.cannot_end_with_this_regex, cmd_line):
        raise ec.ScriptSyntaxError(
            cmd_line,
            ("Cannot end command with special characters, apart from \']\' "
             "and \')\'. Exiting..."))

    assignment_call_count = len(re.findall(cmdu.assign_regex, cmd_line))
    sampled_as_call_count = len(re.findall(cmdu.sampled_as_regex, cmd_line))
    deterministic_call_count = \
        len(re.findall(cmdu.deterministic_regex, cmd_line))

    if (assignment_call_count +
            + sampled_as_call_count +
            + deterministic_call_count) > 1:
        print("\n\nCommand line with problem:\n  " + cmd_line)
        raise ec.ScriptSyntaxError(
            cmd_line,
            ("Only _one_ of the following can be done:\n    (1) [variable] <- "
             "[value]\n    (2) [variable] ~ [sampling distribution]"
             "\n    (3) [variable] := [deterministic function]"))

    elif (assignment_call_count +
          + sampled_as_call_count +
          + deterministic_call_count) == 0:
        print("\n\nCommand line with problem:\n  " + cmd_line)
        raise ec.ScriptSyntaxError(
            cmd_line,
            ("One of the following must be done:\n    (1) [variable]"
             " <- [value]\n    (2) [variable] ~ [sampling distribution]"
             "\n    (3) [variable] := [deterministic function]"))

    ##########################
    # Executing command line #
    ##########################

    # (1) variable assignment
    elif assignment_call_count == 1:
        try:
            rv_name, operator_str, rv_spec = \
                re.split(cmdu.assign_regex, cmd_line)

        except ValueError:
            raise ec.ScriptSyntaxError(
                cmd_line,
                ("Something went wrong during constant value assignment. "
                 "Could not tokenize constant node name, \'<-\' operator, "
                 "and the value being assigned"))

        try:
            parse_variable_assignment(
                dag_obj,
                rv_name,
                rv_spec,
                cmd_line)

        except ec.ScriptSyntaxError:
            # re-raises the ScriptSyntaxError exception
            # if raised within parse_variable_assignment
            raise

    # (2) sampling distribution assignment
    elif sampled_as_call_count == 1:
        rv_name: str = ""
        operator_str: str = ""
        rv_spec: str = ""

        try:
            rv_name, operator_str, rv_dn_spec = \
                re.split(cmdu.sampled_as_regex, cmd_line)

        except ValueError:
            raise ec.ScriptSyntaxError(
                cmd_line,
                ("Something went wrong during sampling distribution "
                 "assignment. Could not tokenize stochastic node "
                 "name, \'~\' operator, and sampling distribution "
                 "specification"))

        try:
            parse_samp_dn_assignment(
                dag_obj,
                rv_name,
                rv_dn_spec,
                cmd_line)

        except ec.ParseDnInitFailError as e:
            raise ec.ScriptSyntaxError(
                cmd_line,
                e.message)

        except ec.ScriptSyntaxError as e:
            # re-raises the ScriptSyntaxError exception
            # if raised within parse_samp_dn_assignment
            raise

    # (3) deterministic function assignment
    elif deterministic_call_count == 1:
        try:
            det_nd_name, operator_str, det_fn_spec = \
                re.split(cmdu.deterministic_regex, cmd_line)

        except ValueError:
            raise ec.ScriptSyntaxError(
                cmd_line,
                ("Something went wrong during deterministic function "
                 "assignment. Could not tokenize deterministic node "
                 "name, \':=\' operator, and deterministic function "
                 "specification"))

        try:
            parse_deterministic_function_assignment(
                dag_obj,
                det_nd_name,
                det_fn_spec,
                cmd_line)

        except ec.ParseDnInitFailError as e:
            raise ec.ScriptSyntaxError(
                cmd_line,
                e.message)

        except ec.ScriptSyntaxError as e:
            # re-raises the ScriptSyntaxError exception
            # if raised within parse_deterministic_function_assignment
            raise

    # cmd_line is valid, returning (used by GUI, when displaying command in History window)!
    return cmd_line


#####################
# Parsing functions #
#####################

# function called by (1. variable assignment)
def parse_variable_assignment(
        dag_obj: pgm.DirectedAcyclicGraph,
        stoch_node_name: str,
        stoch_node_spec: str,
        cmd_line: str) -> None:
    """Create and add constant node to DAG.

    Called when '<-' operator is used.
    This node will not be sampled and is in practice a constant node.

    Args:
        dag_obj (DirectedAcyclicGraph): Object that will hold all nodes
            created by user via commands.
        stoch_node_name (str): Name of stochastic node being created.
        stoch_node_spec (str): Stochastic node specification string
            (whatever is right of '<-' operator in PJ command).
        cmd_line (str): Command line string provided by user through
            static script or via GUI.
    """

    def create_add_stoch_node_dag(a_stoch_node_name: str,
                                  sample_size: int,
                                  a_val_obj_list: ty.List[ty.Any],
                                  a_ct_fn_obj: pgm.ConstantFn):

        replicate_size: int = 1
        if a_ct_fn_obj is not None:
            replicate_size = a_ct_fn_obj.n_repl

        stoch_node = pgm.StochasticNodeDAG(
            a_stoch_node_name,
            sample_size=sample_size,
            replicate_size=replicate_size,
            value=a_val_obj_list,
            returned_from=a_ct_fn_obj,
            clamped=True)

        dag_obj.add_node(stoch_node)

    #############
    # IMPORTANT #
    #############
    # arguments of dn/det parameters will come out as lists
    values_list: ty.List[ty.Any] = list()
    ct_fn_obj: pgm.ConstantFn = None

    # if argument is a vector of values
    if re.match(cmdu.vector_value_regex, stoch_node_spec):
        values_list = cmdu.parse_val_vector(stoch_node_spec)

    # if scalar variable
    elif re.match(cmdu.int_or_float_regex, stoch_node_spec):
        if len(re.findall(cmdu.int_or_float_regex, stoch_node_spec)) > 1:
            raise ec.ScriptSyntaxError(
                stoch_node_spec,
                ("Something went wrong during variable assignment. If a value, "
                 "it must be a single integer or float."))

        if len(re.findall(r"\.", stoch_node_spec)) > 1:
            raise ec.ScriptSyntaxError(
                stoch_node_spec,
                ("Something went wrong during variable assignment. It looks "
                 "like there were two dots \".\" when only one is allowed."))

        values_list.append(stoch_node_spec)
    
    # if the value is the return of a function (e.g., read_tree())
    elif re.search(cmdu.sampling_dn_spec_regex, stoch_node_spec) is not None:
        constant_fn_name, constant_fn_spec = \
            re.search(cmdu.sampling_dn_spec_regex, stoch_node_spec).groups()

        # parses, e.g., "par1=arg1, par2=arg2" into { par1:arg1, par2:arg2 }
        spec_dict, parent_pgm_nodes = \
            cmdu.parse_spec(dag_obj, constant_fn_spec, cmd_line)

        ################################################
        # Create the constant function object          #
        #                                              #
        # Inside PJCtFnGrammar, we check that constant #
        # function is ok                               #
        ################################################

        try:
            # exists outside try!
            ct_fn_obj = \
                ctgrammar.PJCtFnGrammar.create_ct_fn_obj(
                    constant_fn_name,
                    spec_dict)

        except (ec.ParseNotAParameterError,
                ec.ParseRequireSingleValueError,
                ec.ParseRequireIntegerError,
                ec.ParseRequireNumericError,
                ec.ParseMissingParameterError) as e:
            raise ec.ParseCtFnInitFailError(constant_fn_name, e.message)
        
    elif re.match(cmdu.character_value_regex, stoch_node_spec):
        if len(re.findall(cmdu.character_value_regex, stoch_node_spec)) > 1:
            raise ec.ScriptSyntaxError(
                stoch_node_spec,
                ("Something went wrong during variable assignment. If being "
                 "assigned the value of another variable, it must be a single "
                 "one."))
        
        values_list.append(stoch_node_spec)

    else:
        raise ec.ScriptSyntaxError(
            stoch_node_spec,
            ("Something went wrong during variable assignment "
             "Could not find both the name of a function "
             "(e.g., 'read_tree') and its specification (e.g., "
             "'(file=\"examples/geosse_dummy_tree1\", node_name_attr=\"index\")')."))

    val_obj_list = cmdu.val_or_obj(dag_obj, values_list)
    n_samples = len(val_obj_list)

    create_add_stoch_node_dag(stoch_node_name,
                              n_samples,
                              val_obj_list,
                              ct_fn_obj)


# function called by (2. sampling distribution assignment)
def parse_samp_dn_assignment(
        dag_obj: pgm.DirectedAcyclicGraph,
        stoch_node_name: str,
        stoch_node_dn_spec: str,
        cmd_line: str) -> None:
    """Create and add stochastic node to graphical model.

    Called when '~' operator is used. This node will be sampled from 
    distribution and is a random variable.

    Args:
        dag_obj (DirectedAcyclicGraph): Object that will hold all nodes
            created by user via commands.
        stoch_node_name (str): Name of stochastic node being created
        stoch_node_dn_spec (str): Specification string for distribution
            from which stochastic node will be sampled (whatever is
            right of '~' operator in PJ command)
        cmd_line (str): Command line string provided by user through
            static script or via GUI
    """

    def create_add_rv_pgm(
            a_stoch_node_name: str,
            sample_size: int,
            replicate_size: int,
            a_dn_obj: pgm.DistrForSampling,
            parent_pgm_nodes: ty.List[pgm.NodeDAG],
            clamped: bool):

        # set dn inside rv, then call .sample
        stoch_node_dag = pgm.StochasticNodeDAG(
            a_stoch_node_name,
            sample_size=sample_size,
            replicate_size=replicate_size,
            sampled_from=a_dn_obj,
            parent_nodes=parent_pgm_nodes,
            clamped=clamped)

        dag_obj.add_node(stoch_node_dag)


    if re.search(cmdu.sampling_dn_spec_regex, stoch_node_dn_spec) is None:
        raise ec.ScriptSyntaxError(
            cmd_line,
            ("Something went wrong during sampling distribution "
             "specification. Could not find either the name of a "
             "distribution (e.g., \'normal\') or its specification "
             "(e.g., \'(mean=0.0, sd=1.0)\'), or both."))

    # ok!
    else:
        dn_name, dn_spec = \
            re.search(
                cmdu.sampling_dn_spec_regex,
                stoch_node_dn_spec).groups()

        if dn_name not in dngrammar.PJDnGrammar.dn_grammar_dict:
            raise ec.ScriptSyntaxError(
                cmd_line,
                ("Something went wrong during sampling distribution "
                 "specification. Distribution name not recognized."))

        # parses, e.g., "par1=arg1, par2=arg2" into { par1:arg1, par2:arg2 }
        spec_dict, parent_pgm_nodes = \
            cmdu.parse_spec(dag_obj, dn_spec, cmd_line)

        ####################################################
        # Create the sampling distribution object          #
        #                                                  #
        # Inside create_dn_obj, we check that distribution #
        # are ok                                           #
        ####################################################
        try:
            # exists outside try!
            dn_obj = \
                dngrammar.PJDnGrammar.create_dn_obj(dn_name, spec_dict)

        except (ec.ParseNotAParameterError,
                ec.ParseRequireSingleValueError,
                ec.ParseRequireIntegerError,
                ec.ParseRequireNumericError,
                ec.ParseMissingParameterError) as e:
            raise ec.ParseDnInitFailError(dn_name, e.message)

        ####################
        # Clamping, if any #
        ####################

        # if not provided by user, no clamping
        clamped = False

        if "clamp" in spec_dict:
            if isinstance(spec_dict["clamp"], list):
                clamp = spec_dict["clamp"][0]

            if clamp in ("\"true\"", "\"T\"", "\"True\""):
                clamped = True

            if clamp in ("\"false\"", "\"F\"", "\"False\""):
                clamped = False

        ###################################################
        # Get other quantities we need to initialize node #
        # and then add it to the PGM object               #
        ###################################################

        # if not provided by user, we sample once with one replicate
        sample_size = 1
        repl_size = 1

        if "n" in spec_dict:
            # if user passes 'n' as a simple string,
            # e.g., 'normal(n=2...)'
            if isinstance(spec_dict["n"][0], str):
                # this check was already done inside
                # .create_dn_obj() above
                # try:
                sample_size = int(spec_dict["n"][0])
                # except:
                    # pass

            # if user passes 'r' as a node,
            # e.g., 'n_sim <- 2', then 'normal(n=n_sim...)
            elif isinstance(spec_dict["n"][0], pgm.StochasticNodeDAG):
                # this check and the try-except below were already done inside
                # .create_dn_obj() above
                #
                # if len(spec_dict["n"][0].value) > 1:
                #     raise ec.RequireScalarError(dn_name, "n")
                # try:
                # [0] because it is vectorized
                sample_size = int(spec_dict["n"][0].value[0])
                # except:
                #     raise ec.FunctionArgError(
                #         "n",
                #         ("Was expecting either a scalar string, or a node "
                #          "storing a scalar constant. Distribution "
                #          "discrete_sse() could not be initialized."))

        if "nr" in spec_dict:
            # if user passes 'nr' as a simple string,
            # e.g., 'normal(...nr=2...)'
            if isinstance(spec_dict["nr"][0], str):
                # this check was already done inside
                # .create_dn_obj() above
                # try:
                repl_size = int(spec_dict["nr"][0])
                # except ValueError:
                    # pass  # user did not provide replicate size

            # if user passes 'nr' as a node,
            # e.g., 'n_rep <- 2', then 'normal(...nr=n_rep...)
            elif isinstance(spec_dict["nr"][0], pgm.StochasticNodeDAG):
                # this check and the try-except below were already done inside
                # .create_dn_obj() above
                #
                # if len(spec_dict["nr"][0].value) > 1:
                #     raise ec.RequireScalarError(dn_name, "nr")
                # try:
                # [0] because it is vectorized
                repl_size = int(spec_dict["nr"][0].value[0])
                # except:
                #     raise ec.FunctionArgError(
                #         "nr",
                #         ("Was expecting either a scalar string, or a node "
                #          "storing a scalar constant. Distribution "
                #          "discrete_sse() could not be initialized."))

        try:
            create_add_rv_pgm(
                stoch_node_name,
                sample_size,
                repl_size,
                dn_obj,
                parent_pgm_nodes,
                clamped)

        except ec.DAGCannotAddNodeError as e:
            raise ec.ParseDnInitFailError(dn_name, e.message)


# function called by (3. deterministic function assignment)
def parse_deterministic_function_assignment(
        dag_obj: pgm.DirectedAcyclicGraph,
        det_nd_name: str,
        det_node_fn_spec: str,
        cmd_line: str) -> None:
    """
    Create DeterministicNodeDAG instance from command string with
    ':=' operator, then add it to ProbabilisticGraphiclModel
    instance. This node is not sampled (not a random variable) and is
    deterministically initialized via a deterministic function call.

    Args:
        dag_obj (DirectedAcyclicGraph): DAG object that will hold all
            nodes created by user via commands.
        det_nd_name (str): Name of deterministic node being created
        det_node_fn_spec (str): Specification string for deterministic
            function that will specify the deterministic node (whatever
            is right of ':=' operator in PJ command).
        cmd_line (str): Command line string provided by user through static
            script or via GUI.
    """

    def create_add_det_nd_pgm(det_nd_name: str,
                              det_obj: ty.Any,
                              parent_pgm_nodes: ty.List[pgm.NodeDAG]):

        det_nd_pgm = pgm.DeterministicNodeDAG(det_nd_name,
                                              value=det_obj,
                                              parent_nodes=parent_pgm_nodes)
        dag_obj.add_node(det_nd_pgm)

        # deterministic node is of class DeterministicNodeDAG, which
        # derives NodeDAG -- we do not need to initialize NodeDAG
        # as when a new StochasticNodeDAG is created (see above in
        # create_add_rv_pgm())
        #
        # this check is just to make sure we are adding a class
        # deriving from NodeDAG
        # if isinstance(det_obj, NodeDAG):
        #     det_obj.node_name = det_nd_name
        #     det_obj.parent_nodes = parent_pgm_nodes
        #     det_nd_dag = det_obj
        #     dag_obj.add_node(det_nd_dag)

    if re.search(cmdu.sampling_dn_spec_regex, det_node_fn_spec) is None:
        raise ec.ScriptSyntaxError(
            det_node_fn_spec,
            ("Something went wrong during deterministic function "
             "specification. Could not find both the name of a function "
             "(e.g., \'sse_prob\') and its specification (e.g., "
             "\'(name=prob1, value=1.0, state=0)\')."))

    # ok!
    else:
        det_fn_matches = re.search(cmdu.sampling_dn_spec_regex, det_node_fn_spec)

        if det_fn_matches is not None:
            det_fn_name, det_fn_spec = det_fn_matches.groups()

            if det_fn_name not in detgrammar.PJDetFnGrammar.det_fn_grammar_dict:
                raise ec.ScriptSyntaxError(
                    det_node_fn_spec,
                    ("Something went wrong during sampling distribution "
                     "specification. Distribution name not recognized."))

            # parses, e.g., "par1=arg1, par2=arg2" into { par1:arg1, par2:arg2 }
            spec_dict, parent_pgm_nodes = cmdu.parse_spec(dag_obj, det_fn_spec, cmd_line)

            try:
                # exists outside try!
                det_fn_obj = \
                    detgrammar.PJDetFnGrammar.create_det_fn_obj(
                        det_fn_name, spec_dict)

            except (ec.ParseMissingSpecificationError,
                    ec.ParseNotAParameterError,
                    ec.ParseMissingParameterError) as e:
                raise ec.ParseDetFnInitFailError(det_fn_name, e.message)

            # except ec.ParseDetFnInitFailError as e:
            #     # re-raises the DetFnInitFailError exception from
            #     # within .create_det_fn_obj()
            #     raise

            create_add_det_nd_pgm(det_nd_name, det_fn_obj, parent_pgm_nodes)


if __name__ == "__main__":

    script_str1 = "a <- 1.0"
    script_str2 = "a <- [1, 2, 3]"
    script_str3 = "a <- 1\nb <- 2\nc <- [a, b, 3]"
    script_str4 = "m <- [1,2,3,4,5,6,7,8,9,10]\nrv ~ lognormal(n=10, mean=m, sd=1.0)"
    script_str5 = ("l0 <- 1.0\na := sse_rate(name=\"lambda\", value=l0,"
                   " event=\"b_speciation\", states=[2, 0, 1])")
    script_str6 = ("l0rate := sse_rate(name=\"lambda\", value=1.0, "
                   "event=\"w_speciation\")")
    script_str6_2 = ("l0rate := sse_rate(name=\"lambda\", "
                     "value=[1.0, 1.1], event=\"w_speciation\")")
    script_str7 = script_str6 + \
        "\nstash := sse_stash(flat_rate_mat=[l0rate], n_states=1, n_epochs=1)"
    script_str8 = script_str7 + \
        ("\ntr ~ discrete_sse(n=2, stash=stash, start_state=[0,0], "
         "stop=\"age\", stop_value=2.0, origin=\"true\", cond_spn=\"false\","
         " cond_surv=\"true\")")
    script_str9 = "n <- 1000\nrv ~ unif(n=n, min=0.0, max=1.0)"

    # bisse
    script_str10 = \
        "l0rate := sse_rate(name=\"lambda0\", value=1.0, states=[0,0,0], event=\"w_speciation\")\n" + \
        "l1rate := sse_rate(name=\"lambda1\", value=1.0, states=[1,1,1], event=\"w_speciation\")\n" + \
        "m0rate := sse_rate(name=\"mu0\", value=0.8, states=[0], event=\"extinction\")\n" + \
        "m1rate := sse_rate(name=\"mu1\", value=0.8, states=[1], event=\"extinction\")\n" + \
        "q01rate := sse_rate(name=\"q01\", value=0.6, states=[0,1], event=\"transition\")\n" + \
        "q10rate := sse_rate(name=\"q10\", value=0.6, states=[1,0], event=\"transition\")\n" + \
        "stash := sse_stash(flat_rate_mat=[l0rate, l1rate, m0rate, m1rate, q01rate, q10rate], n_states=2, n_epochs=1)\n" + \
        "tr ~ discrete_sse(n=1, stash=stash, start_state=[0], stop=\"age\", stop_value=1.0, origin=\"true\", cond_spn=\"false\", cond_surv=\"true\")"

    script_str11 = "m <- [1.0, 1.1]\ns <- [1.0, 0.9]\na ~ lognormal(n=2, mean=m, sd=s, log_space=\"false\")"
    script_str12 = "a ~ lognormal(n=3, mean=[0.0, 0.01, 0.011], sd=[1.0], log_space=\"true\")"
    script_str13 = "a ~ lognormal(mean=0.0, sd=1.0)"
    script_str14 = "a ~ lognormal(n=3, mean=[0.0, 0.01, 0.011], sd=1.0)"
    script_str15 = "a ~ normal(n=3, mean=[0.0, 0.01, 0.011], sd=[1.0])"
    script_str16 = "a ~ exponential(n=3, rate=1.0)"
    script_str17 = "a ~ exponential(rate=1.0, rate_parameterization=\"true\")"
    script_str18 = "l0 ~ unif(n=1, min=0.8, max=0.8)\nl0rate := sse_rate(name=\"lambda0\", value=l0, states=[0,0,0], event=\"w_speciation\")"
    script_str19 = "l0 ~ unif(n=3, min=0.8, max=0.8)\nl0rate := sse_rate(name=\"lambda0\", value=l0, states=[0,0,0], event=\"w_speciation\")"
    script_str20 = script_str18 + "\nstash :=stash(flat_rate_mat=[l0rate], n_states=1, n_epochs=1)" + \
        "\ntr ~ discrete_sse(n=3, stash=stash, start_state=[0], stop=\"age\", stop_value=2.0, origin=\"true\", cond_spn=\"false\", cond_surv=\"true\")"
    script_str21 = script_str19 + ("\nstash :=stash(flat_rate_mat=[l0rate], n_states=1, n_epochs=1)"
                                   "\ntr ~ discrete_sse(n=3, stash=stash, start_state=[0,0,0], stop=\"age\", stop_value=2.0, origin=\"true\", cond_spn=\"false\", cond_surv=\"true\")")
    script_str22 = "rv1 ~ lognormal(n=2, mean=0.0, sd=0.25)\nrv2 ~ lognormal(n=2, mean=-1.5, sd=0.25)\nrv3 ~ normal(n=2, nr=10, mean=rv1, sd=rv2, clamp=\"true\")"

    # time-het yule
    script_str23 = \
        ("l0ratet1 := sse_rate(name=\"lambdat1\", value=1.0, event=\"w_speciation\", epoch=1)\n"
         "l0ratet2 := sse_rate(name=\"lambdat2\", value=0.25, event=\"w_speciation\", epoch=2)\n"
         "l0ratet3 := sse_rate(name=\"lambdat3\", value=3.0, event=\"w_speciation\", epoch=3)\n"
         "l0ratet4 := sse_rate(name=\"lambdat4\", value=0.4, event=\"w_speciation\", epoch=4)\n"
         "stash := sse_stash(flat_rate_mat=[l0ratet1, l0ratet2, l0ratet3, l0ratet4], n_states=1, n_epochs=4, seed_age=3.0, epoch_age_ends=[2.2, 1.2, 0.7])")
    script_str24 = script_str23 \
        + ("\ntr ~ discrete_sse(n=2, stash=stash, start_state=[0], stop=\"age\", stop_value=3.0, origin=\"true\", cond_spn=\"false\", cond_surv=\"true\")")

    script_str25 = "a <- [1,2,3]\nb <- a"
    script_str26 = "a <- [1,2,3]\nb <- [a, 4, 5, 6]"
    script_str27 = "rv9 ~ exponential(n=5, rate=[0.8, 0.9, 1.0, 1.1, 1.2], clamp=\"true\")"

    script_str28 = \
        "birth_rate := sse_rate(name=\"lambda\", value=1.0, states=[0,0,0], event=\"w_speciation\")\n" + \
        "death_rate := sse_rate(name=\"mu\", value=0.5, states=[0], event=\"extinction\")\n" + \
        "fossil_rate := sse_rate(name=\"psi\", value=0.8, states=[0], event=\"anc_sampling\")\n" + \
        "stash := sse_stash(flat_rate_mat=[birth_rate,death_rate,fossil_rate], n_states=1, n_epochs=1)\n" + \
        "trs ~ discrete_sse(n=5, stash=stash, start_state=[0], stop=\"age\", stop_value=3.0, origin=\"true\", cond_spn=\"true\", cond_surv=\"true\")"

    # time-het bd
    script_str29 = \
        "birth_rate_t0 <- 0.5\n" + \
        "death_rate_t0 <- 0.25\n" + \
        "birth_rate_t1 <- 1.0\n" + \
        "death_rate_t1 <- 0.4\n" + \
        "det_birth_rate_t0 := sse_rate(name=\"lambda_t0\", value=birth_rate_t0, states=[0,0,0], event=\"w_speciation\", epoch=1)\n" + \
        "det_death_rate_t0 := sse_rate(name=\"mu_t0\", value=death_rate_t0, states=[0], event=\"extinction\", epoch=1)\n" + \
        "det_birth_rate_t1 := sse_rate(name=\"lambda_t1\", value=birth_rate_t1, states=[0,0,0], event=\"w_speciation\", epoch=2)\n" + \
        "det_death_rate_t1 := sse_rate(name=\"mu_t1\", value=death_rate_t1, states=[0], event=\"extinction\", epoch=2)\n" + \
        "stash := sse_stash(flat_rate_mat=[det_birth_rate_t0, det_death_rate_t0, det_birth_rate_t1, det_death_rate_t1], n_states=1, n_epochs=2, epoch_age_ends=[1.5], seed_age=3.0)"

    script_str30 = "a ~ lognormal(n=10, mean=[0.0, 0.01], sd=[1.0, 0.9])"
    script_str31 = "a ~ lognormal(n=2, mean=[0.0, 0.01], sd=[1.0, 0.9, 0.8])"
    script_str32 = "l0 ~ unif(n=3, nr=2, min=0.9, max=1.1)\nl0rate := sse_rate(name=\"lambda\", value=l0, event=\"w_speciation\")"
    script_str33 = script_str23 + "\ntr ~ discrete_sse(n=1, stash=stash, start_state=[0], stop=\"size\", stop_value=3.0, origin=\"true\", cond_spn=\"false\", cond_surv=\"true\")"
    script_str34 = script_str23 + "\ntr ~ discrete_sse(n=1, stash=stash, start_state=[0], stop=\"age\", stop_value=2.0, origin=\"true\", cond_spn=\"false\", cond_surv=\"true\")"

    # time-het geosse
    script_str35 = \
        "w_birth_rate0_t0 <- 1.5\n" + \
        "death_rate0_t0 <- 0.25\n" + \
        "trans_rate02_t0 <- 0.0\n" + \
        "w_birth_rate1_t0 <- 1.0\n" + \
        "death_rate1_t0 <- 0.5\n" + \
        "trans_rate12_t0 <- 0.4\n" + \
        "b_birth_rate201_t0 <- 2.0\n" + \
        "b_birth_rate202_t0 <- 1.0\n" + \
        "b_birth_rate212_t0 <- 1.0\n" + \
        "trans_rate20_t0 <- 0.6\n" + \
        "trans_rate21_t0 <- 0.6\n" + \
        "w_birth_rate0_t1 <- 0.85\n" + \
        "death_rate0_t1 <- 0.5\n" + \
        "trans_rate02_t1 <- 0.8\n" + \
        "w_birth_rate1_t1 <- 1.5\n" + \
        "death_rate1_t1 <- 0.75\n" + \
        "trans_rate12_t1 <- 0.6\n" + \
        "b_birth_rate201_t1 <- 2.5\n" + \
        "b_birth_rate202_t1 <- 1.25\n" + \
        "b_birth_rate212_t1 <- 1.25\n" + \
        "trans_rate20_t1 <- 0.75\n" + \
        "trans_rate21_t1 <- 0.75\n" + \
        "det_w_birth_rate0_t0 := sse_rate(name=\"lambda0_t0\", value=w_birth_rate0_t0, states=[0,0,0], event=\"w_speciation\", epoch=1)\n" + \
        "det_death_rate0_t0 := sse_rate(name=\"mu0_t0\", value=death_rate0_t0, states=[0], event=\"extinction\", epoch=1)\n" + \
        "det_trans_rate02_t0 := sse_rate(name=\"q02_t0\", value=trans_rate02_t0, states=[0,2], event=\"transition\", epoch=1)\n" + \
        "det_w_birth_rate1_t0 := sse_rate(name=\"lambda1_t0\", value=w_birth_rate1_t0, states=[1,1,1], event=\"w_speciation\", epoch=1)\n" + \
        "det_death_rate1_t0 := sse_rate(name=\"mu1_t0\", value=death_rate1_t0, states=[1], event=\"extinction\", epoch=1)\n" + \
        "det_trans_rate12_t0 := sse_rate(name=\"q12_t0\", value=trans_rate12_t0, states=[1,2], event=\"transition\", epoch=1)\n" + \
        "det_b_birth_rate201_t0 := sse_rate(name=\"lambda201_t0\", value=b_birth_rate201_t0, states=[2,0,1], event=\"bw_speciation\", epoch=1)\n" + \
        "det_b_birth_rate202_t0 := sse_rate(name=\"lambda202_t0\", value=b_birth_rate202_t0, states=[2,0,2], event=\"bw_speciation\", epoch=1)\n" + \
        "det_b_birth_rate212_t0 := sse_rate(name=\"lambda212_t0\", value=b_birth_rate212_t0, states=[2,1,2], event=\"bw_speciation\", epoch=1)\n" + \
        "det_trans_rate20_t0 := sse_rate(name=\"q20_t0\", value=trans_rate20_t0, states=[2,0], event=\"transition\", epoch=1)\n" + \
        "det_trans_rate21_t0 := sse_rate(name=\"q21_t0\", value=trans_rate21_t0, states=[2,1], event=\"transition\", epoch=1)\n" + \
        "det_w_birth_rate0_t1 := sse_rate(name=\"lambda0_t1\", value=w_birth_rate0_t1, states=[0,0,0], event=\"w_speciation\", epoch=2)\n" + \
        "det_death_rate0_t1 := sse_rate(name=\"mu0_t1\", value=death_rate0_t1, states=[0], event=\"extinction\", epoch=2)\n" + \
        "det_trans_rate02_t1 := sse_rate(name=\"q02_t1\", value=trans_rate02_t1, states=[0,2], event=\"transition\", epoch=2)\n" + \
        "det_w_birth_rate1_t1 := sse_rate(name=\"lambda1_t1\", value=w_birth_rate1_t1, states=[1,1,1], event=\"w_speciation\", epoch=2)\n" + \
        "det_death_rate1_t1 := sse_rate(name=\"mu1_t1\", value=death_rate1_t1, states=[1], event=\"extinction\", epoch=2)\n" + \
        "det_trans_rate12_t1 := sse_rate(name=\"q12_t1\", value=trans_rate12_t1, states=[1,2], event=\"transition\", epoch=2)\n" + \
        "det_b_birth_rate201_t1 := sse_rate(name=\"lambda201_t0\", value=b_birth_rate201_t1, states=[2,0,1], event=\"bw_speciation\", epoch=2)\n" + \
        "det_b_birth_rate202_t1 := sse_rate(name=\"lambda202_t0\", value=b_birth_rate202_t1, states=[2,0,2], event=\"bw_speciation\", epoch=2)\n" + \
        "det_b_birth_rate212_t1 := sse_rate(name=\"lambda212_t0\", value=b_birth_rate212_t1, states=[2,1,2], event=\"bw_speciation\", epoch=2)\n" + \
        "det_trans_rate20_t1 := sse_rate(name=\"q20_t1\", value=trans_rate20_t1, states=[2,0], event=\"transition\", epoch=2)\n" + \
        "det_trans_rate21_t1 := sse_rate(name=\"q21_t1\", value=trans_rate21_t1, states=[2,1], event=\"transition\", epoch=2)\n" + \
        "stash := sse_stash(flat_rate_mat=[det_w_birth_rate0_t0, det_death_rate0_t0, det_trans_rate02_t0, det_w_birth_rate1_t0, det_death_rate1_t0, det_trans_rate12_t0, det_b_birth_rate201_t0, det_b_birth_rate202_t0, det_b_birth_rate212_t0, det_trans_rate20_t0, det_trans_rate21_t0, det_w_birth_rate0_t1, det_death_rate0_t1, det_trans_rate02_t1, det_w_birth_rate1_t1, det_death_rate1_t1, det_trans_rate12_t1, det_b_birth_rate201_t1, det_b_birth_rate202_t1, det_b_birth_rate212_t1, det_trans_rate20_t1, det_trans_rate21_t1], n_states=3, n_epochs=2, epoch_age_ends=[1.0], seed_age=3.0)\n" + \
        "trs ~ discrete_sse(n=2, stash=stash, start_state=[0], stop=\"age\", stop_value=3.0, origin=\"true\")\n"

    script_str36 = "a ~ lognormal(n=2, nr=10, mean=10.0, sd=1.0)"

    # yule with n_sim as a model node
    script_str37 = \
        "n_sim <- 2\n" \
        + "n_rep <- 2\n" \
        + "birth_rate <- [0.8, 0.9]\n" \
        + "det_birth_rate := sse_rate(name=\"lambda\", value=birth_rate, states=[0,0,0], event=\"w_speciation\")\n" \
        + "stash := sse_stash(flat_rate_mat=[det_birth_rate], n_states=1, n_epochs=1)\n" \
        + "trs ~ discrete_sse(n=n_sim, nr=n_rep, stash=stash, start_state=[0,0], stop=\"age\", stop_value=1.0, origin=\"true\")"

    # yule with incomplete sampling
    script_str38 = \
        "n_sim <- 2\n" \
        + "n_rep <- 2\n" \
        + "birth_rate <- [0.8, 0.9]\n" \
        + "det_birth_rate := sse_rate(name=\"lambda\", value=birth_rate, states=[0,0,0], event=\"w_speciation\")\n" \
        + "sampling_rate <- [0.5, 0.5]\n" \
        + "det_sampling_rate := sse_prob(name=\"rho\", value=sampling_rate, state=[0])\n" \
        + "stash := sse_stash(flat_rate_mat=[det_birth_rate], flat_prob_mat=[det_sampling_rate], n_states=1, n_epochs=1)\n" \
        + "trs ~ discrete_sse(n=n_sim, nr=n_rep, stash=stash, start_state=[0,0], stop=\"age\", stop_value=1.0, origin=\"true\")"

    # fbd with incomplete sampling
    script_str39 = \
        "birth_rate <- 1.0\n" \
        + "death_rate <- 0.2\n" \
        + "fossil_rate <- 0.5\n" \
        + "det_birth_rate := sse_rate(name=\"lambda\", value=birth_rate, states=[0, 0, 0], event=\"w_speciation\")\n" \
        + "det_death_rate := sse_rate(name=\"mu\", value=death_rate, states=[0], event=\"extinction\")\n" \
        + "det_fossil_rate := sse_rate(name=\"psi\", value=fossil_rate, states=[0], event=\"anc_sampling\")\n" \
        + "sampling_prob <- 0.5\n" \
        + "det_sampling_prob := sse_prob(name=\"rho\", value=sampling_prob, state=0)\n" \
        + "stash := sse_stash(flat_rate_mat=[det_birth_rate, det_death_rate, det_fossil_rate]," \
        + "flat_prob_mat=[det_sampling_prob], n_states=1, n_epochs=1)\n" \
        + "ori <- 3\n" \
        + "trs ~ discrete_sse(n=1, stash=stash, start_state=[0], stop=\"age\", stop_value=ori, origin=\"true\"," \
        + "cond_spn=\"true\", cond_surv=\"true\")"

    script_str40 = \
        'tr <- read_tree(string="((sp1[&index=1]:1.0,sp2[&index=2]:1.0)[&index=4]:1.0,sp3[&index=3]:2.0)[&index=5];", node_name_attr="index")'
    
    script_str41 = \
        'tr <- read_tree(string="((sp1:1.0,sp2:1.0):1.0,sp3:2.0);")'

    script_str42 = \
        'tr <- read_tree(file_path="examples/trees_maps_files/tree_to_read.tre", node_name_attr="index")'
    
    script_str43 = \
        'tr <- read_tree(file_path="examples/trees_maps_files/trees_to_read.tre", node_name_attr="index")'
        # tr <- read_tree(file_path="examples/trees_maps_files/trees_to_read.tre", node_name_attr="index")
    
    # for copying and pasting in GUI:
    #
    # vectorized atomic rate parameter
    # l0rate := sse_rate(name="lambda", value=[1.0, 1.1], event="w_speciation")
    # stash := stash(flat_rate_mat=[l0rate], n_states=1, n_epochs=1)

    # yule
    # l0 <- [0.9, 1.0, 1.1]
    # l0 ~ unif(n=3, min=0.9, max=1.1)
    # l0rate := sse_rate(name="lambda", value=l0, event="w_speciation")
    # stash := sse_stash(flat_rate_mat=[l0rate], n_states=1, n_epochs=1)
    # trs ~ discrete_sse(n=3, stash=stash, start_state=[0], stop="age", stop_value=2.0, origin="true", cond_spn="false", cond_surv="true")
    # trs ~ discrete_sse(n=3, nr=2, stash=stash, start_state=[0], stop="age", stop_value=2.0, origin="true", cond_spn="false", cond_surv="true")

    # bisse
    # l0rate := sse_rate(name="lambda0", value=1.0, states=[0,0,0], event="w_speciation")
    # l1rate := sse_rate(name="lambda1", value=1.0, states=[1,1,1], event="w_speciation")
    # m0rate := sse_rate(name="mu0", value=0.8, states=[0,0,0], event="extinction")
    # m1rate := sse_rate(name="mu1", value=0.8, states=[1,1,1], event="extinction")
    # q01rate := sse_rate(name="q01", value=0.6, states=[0,1], event="transition")
    # q10rate := sse_rate(name="q10", value=0.6, states=[1,0], event="transition")
    # stash := sse_stash(flat_rate_mat=[l0rate, l1rate, m0rate, m1rate, q01rate, q10rate], n_states=2, n_epochs=1)
    # trs ~ discrete_sse(n=5, stash=stash, start_state=[0], stop="age", stop_value=2.5, origin="true", cond_spn="false", cond_surv="true")

    # time-het yule
    # l0ratet1 := sse_rate(name="lambdat1", value=1.0, event="w_speciation")
    # l0ratet2 := sse_rate(name="lambdat2", value=0.25, event="w_speciation")
    # l0ratet3 := sse_rate(name="lambdat3", value=3.0, event="w_speciation")
    # l0ratet4 := sse_rate(name="lambdat4", value=0.4, event="w_speciation")
    # stash := sse_stash(flat_rate_mat=[l0ratet1, l0ratet2, l0ratet3, l0ratet4], n_states=1, n_epochs=4, seed_age=3.0, epoch_age_ends=[2.2, 1.2, 0.7])
    # tr ~ discrete_sse(n=10, stash=stash, start_state=[0], stop="age", stop_value=3.0, origin="true", cond_spn="false", cond_surv="true")

    # simple hierarchical
    # rv1 ~ lognormal(n=2, mean=0.0, sd=0.25)
    # rv2 ~ lognormal(n=2, mean=-1.5, sd=0.25)
    # rv3 ~ normal(n=2, nr=10, mean=rv1, sd=rv2, clamp="true")

    # should throw exceptions (in tests)
    script_str_error1 = "u ~ unif(n=1, nr=1, max=1.0)"
    script_str_error2 = "u ~ unif(n=1, rate=1.0)"
    script_str_error3 = "u ~ unif(n=1, min=\"a\", max=2)"
    script_str_error4 = "u ~ unif()"
    # script_str36 = "1.0 + 1.0"
    # script_str37 = "~ lognormal(n=10, mean=0.0, sd=1.0)"
    # script_str38 = "a <-- 1.0"
    # script_str39 = "a <- 1.0 <- b"
    # script_str40 = "a ~ lognormal(n=10, mean=0.0, sd=1.0) ~"
    # script_str41 = "a ~ lognormal(n=10, mean=0.0, sd=1.0) :="
    # script_str42 = "a ~ lognormal(n=10, mean=0.0, sd=1.0) ("
    # script_str43 = "a ~ lognormal(n=10, mean=0.0, sd=1.0) ["
    # script_str44 = "a ~ lognormal(n=10, mean=0.0, sd=1.0) {"
    # script_str45 = "a <- 1.0.0"

    # TODO: this is not caught!
    # script_str47 = "a ~ normal(mean=0.0, sd=0.5, sd=1.0)"

    # file_handle_exception = io.StringIO(script_str36)

    dag = script2dag(script_str39 , in_pj_file=False, random_seed=456)

    for node_name, node_dag in dag.name_node_dict.items():
        if isinstance(node_dag, pgm.StochasticNodeDAG):
            if isinstance(node_dag.value[0], pjtr.AnnotatedTree):
                ann_tr = node_dag.value[0]
                print(ann_tr.tree.as_string(schema="newick"))

                rec_tr = ann_tr.extract_reconstructed_tree(plotting_overhead=True)
                print(rec_tr.as_string(schema="newick"))


                # print(node_dag.value[0].tree.as_string(schema="newick"))
                # print(node_dag.get_node_stats_str(0, len(node_dag.value), 0))

                ######################
                # Preparing plotting #
                ######################
                # fig = Figure(figsize=(11,4.5))

                # note that pjgui uses matplotlib.figure.Figure
                # (which is part of Matplotlib's OOP class library)
                #
                # here, we instead use pyplot's figure, which is the
                # Matlab-like state-machine API
                fig = plt.figure()

                ax = fig.add_axes([0.25, 0.2, 0.5, 0.6])
                ax.patch.set_alpha(0.0)
                ax.xaxis.set_ticks([])
                ax.yaxis.set_ticks([])
                ax.spines['left'].set_visible(False)
                ax.spines['bottom'].set_visible(False)
                ax.spines['right'].set_visible(False)
                ax.spines['top'].set_visible(False)

                pjtr.plot_ann_tree(ann_tr,
                                   ax,
                                   use_age=False,
                                   start_at_origin=True,
                                   sa_along_branches=True,
                                   attr_of_interest="state",
                                   draw_reconstructed=True)
                # sa_along_branches=False

                plt.show()