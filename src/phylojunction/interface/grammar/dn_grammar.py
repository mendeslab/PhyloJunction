import typing as ty

# pj imports
import phylojunction.distribution.dn_parametric as dnpar
import phylojunction.pgm.pgm as pgm
import phylojunction.interface.grammar.dn_discrete_sse_makers as make_dnsse
import phylojunction.utility.exception_classes as ec
from phylojunction.interface.grammar.dn_quasse_makers import make_quasse
# from user_interface.dn_discrete_sse import make_discrete_SSE_dn # https://stackoverflow.com/questions/16981921/relative-imports-in-python-3

__author__ = "Fabio K. Mendes"
__email__ = "fmendes@lsu.edu"


class PJDnGrammar():

    ################################
    #  All available distributions #
    ################################

    dn_grammar_dict: ty.Dict[str, ty.Tuple[str, ...]] = \
        {
            "quasse": ("n", "nr", "birth_rate", "death_rate", "start_trait", "drift", "diffusion",
                       "stop", "stop_value", "method", "k", "runtime_limit", "max_steps", "max_alive"),
            "lognormal":
            tuple(["n", "nr", "meanlog", "sdlog", "log_space"]),
            "normal":
            tuple(["n", "nr", "mean", "sd"]),
            "exponential":
            tuple(["n", "nr", "rate", "rate_parameterization"]),
            "gamma":
            tuple(["n", "nr", "shape", "scale", "rate_parameterization"]),
            "unif":
            tuple(["n", "nr", "min", "max"]),
            "categorical":
            tuple(["n", "nr", "cats", "probs"]),
            "discrete_sse": tuple(["n", "nr", "stop", "stop_value", "origin",
                                   "stash", "start_state", "eps",
                                   "runtime_limit", "max_n_attempts",
                                   "cond_spn", "cond_surv", "cond_obs_both_sides",
                                   "min_rec_taxa", "max_rec_taxa",
                                   "abort_at_alive_count"])
        }

    def __init__(self) -> None:
        pass

    @classmethod
    def grammar_check(cls, dn_id: str, dn_param: str) -> bool:
        "Check if specified distribution is supported."

        # we ignore the clamp parameter-arg pair;
        # they are handled by cmd_parse()
        if dn_param in cls.dn_grammar_dict[dn_id] or dn_param == "clamp":
            return True

        return False

    @classmethod
    def init_return_parametric_dn(
        cls,
        dn_id: str,
        dn_param_dict: ty.Dict[str, ty.List[ty.Union[str, pgm.NodeDAG]]]) \
            -> pgm.DistrForSampling:
        """Create and return parametric distributions for sampling.

        Args:
            dn_id (str): Name of parametric distribution.
            dn_param_dict (dict): Dictionary with distribution
                parameter names as keys, and their values as values.

        Returns:
            DistrForSampling: A DistrForSampling instance.
        """

        n_samples: int = 1
        n_repl: int = 1
        extracted_val_list: ty.List[str] = list()

        # key: parameter name
        # value: paramenter's parent (DAG node) name
        # e.g., { lambda: node_dag1_name }
        parent_node_tracker: ty.Dict[str, str] = dict()

        if not dn_param_dict:
            raise ec.ParseMissingSpecificationError(dn_id)

        # Default values inside each if block! #

        if dn_id == "lognormal":
            ln_mean: ty.List[float] = list()
            ln_sd: ty.List[float] = list()
            ln_log_space: bool = True

            # { mean: node_dag1_name, sd: node_dag2_name, ... }
            if dn_param_dict:
                # val is list
                for arg, val in dn_param_dict.items():
                    # val = val[0] # TODO: deal with vectorization later

                    # needed for building inference specifications
                    if isinstance(val[0], pgm.StochasticNodeDAG):
                        parent_node_tracker[arg] = val[0].node_name

                    # if element in val is string, it remains unchanged,
                    # if NodeDAG, we get its string-fied value
                    extracted_val_list = \
                        pgm.extract_vals_as_str_from_node_dag(val)

                    if not cls.grammar_check("lognormal", arg):
                        raise ec.ParseNotAParameterError(arg)

                    elif arg == "n":
                        if len(extracted_val_list) > 1:
                            raise ec.ParseRequireSingleValueError(
                                dnpar.DnLogNormal.DN_NAME, arg)

                        # only one element always
                        try:
                            n_samples = int(extracted_val_list[0])

                        except ValueError:
                            raise ec.ParseRequireIntegerError(
                                dnpar.DnNormal.DN_NAME, arg)

                    elif arg == "nr":
                        if len(extracted_val_list) > 1:
                            raise ec.ParseRequireSingleValueError(
                                dnpar.DnLogNormal.DN_NAME, arg)

                        # only one element always
                        try:
                            n_repl = int(extracted_val_list[0])

                        except ValueError:
                            raise ec.ParseRequireIntegerError(
                                dnpar.DnNormal.DN_NAME, arg)

                    elif arg == "meanlog":
                        try:
                            ln_mean = [float(v) for v in extracted_val_list]

                        except ValueError:
                            raise ec.ParseRequireNumericError(
                                dnpar.DnNormal.DN_NAME, arg)

                    elif arg == "sdlog":
                        try:
                            ln_sd = [float(v) for v in extracted_val_list]

                        except ValueError:
                            raise ec.ParseRequireNumericError(
                                dnpar.DnNormal.DN_NAME, arg)

                    elif arg == "log_space":
                        if extracted_val_list[0] in ("\"true\"", "\"T\"", "\"True\""):
                            # pars[4] = True
                            ln_log_space = True

                        elif extracted_val_list[0] in ("\"false\"", "\"F\"", "\"False\""):
                            # pars[4] = False
                            ln_log_space = False

            # making sure essential parameters of distribution have been specified
            for par_obj, par_name in ((ln_mean, "meanlog"), (ln_sd, "sdlog")):
                if not par_obj:
                    raise ec.ParseMissingParameterError(par_name)

            # return dnpar.DnLogNormal(pars, parent_node_tracker)
            return dnpar.DnLogNormal(
                n_samples,
                n_repl,
                ln_mean,
                ln_sd,
                ln_log_space,
                parent_node_tracker)

        elif dn_id == "normal":
            norm_mean: ty.List[float] = list()
            norm_sd: ty.List[float] = list()

            # { mean: node_dag1_name, sd: node_dag2_name, ... }
            if dn_param_dict:
                # val is list
                for arg, val in dn_param_dict.items():

                    # needed for building inference specifications
                    if isinstance(val[0], pgm.StochasticNodeDAG):
                        parent_node_tracker[arg] = val[0].node_name

                    extracted_val_list = pgm.extract_vals_as_str_from_node_dag(val)

                    if not cls.grammar_check("normal", arg):
                        raise ec.ParseNotAParameterError(arg)

                    elif arg == "n":
                        if len(extracted_val_list) > 1:
                            raise ec.ParseRequireSingleValueError(
                                dnpar.DnNormal.DN_NAME, arg)

                        try:
                            n_samples = int(extracted_val_list[0])

                        except ValueError:
                            raise ec.ParseRequireIntegerError(
                                dnpar.DnNormal.DN_NAME, arg)

                    elif arg == "nr":
                        if len(extracted_val_list) > 1:
                            raise ec.ParseRequireSingleValueError(
                                dnpar.DnNormal.DN_NAME, arg)

                        try:
                            n_repl = int(extracted_val_list[0])

                        except ValueError:
                            raise ec.ParseRequireIntegerError(
                                dnpar.DnNormal.DN_NAME, arg)

                    elif arg == "mean":
                        try:
                            norm_mean = [float(v) for v in extracted_val_list]

                        except ValueError:
                            raise ec.ParseRequireNumericError(
                                dnpar.DnNormal.DN_NAME, arg)

                    elif arg == "sd":
                        try:
                            norm_sd = [float(v) for v in extracted_val_list]

                        except ValueError:
                            raise ec.ParseRequireNumericError(
                                dnpar.DnNormal.DN_NAME, arg)

            # making sure essential parameters of distribution have been specified
            for par_obj, par_name in ((norm_mean, "mean"), (norm_sd, "sd")):
                if not par_obj:
                    raise ec.ParseMissingParameterError(par_name)

            # return dnpar.DnNormal(pars, parent_node_tracker)
            return dnpar.DnNormal(
                n_samples,
                n_repl,
                norm_mean,
                norm_sd,
                parent_node_tracker)

        elif dn_id == "exponential":
            exp_scale_or_rate: ty.List[float] = list()
            exp_rate_parameterization: bool = True

            if dn_param_dict:
                # val is list
                for arg, val in dn_param_dict.items():

                    # needed for building inference specifications
                    if isinstance(val[0], pgm.StochasticNodeDAG):
                        # needed for building inference specifications
                        parent_node_tracker[arg] = val[0].node_name

                    extracted_val_list = pgm.extract_vals_as_str_from_node_dag(val)

                    if not cls.grammar_check("exponential", arg):
                        raise ec.ParseNotAParameterError(arg)

                    elif arg == "n":
                        if len(extracted_val_list) > 1:
                            raise ec.ParseRequireSingleValueError(
                                dnpar.DnExponential.DN_NAME, arg)
                        n_samples = int(extracted_val_list[0])

                    elif arg == "nr":
                        if len(extracted_val_list) > 1:
                            raise ec.ParseRequireSingleValueError(
                                dnpar.DnExponential.DN_NAME, arg)
                        n_repl = int(extracted_val_list[0])

                    elif arg == "rate":
                        # pars[2] = [float(v) for v in _extracted_val]
                        exp_scale_or_rate = [float(v) for v in extracted_val_list]

                    elif arg == "rate_parameterization":
                        if extracted_val_list[0] in ("\"true\"", "\"T\"", "\"True\""):
                            # pars[3] = True
                            exp_rate_parameterization = True

                        elif extracted_val_list[0] in ("\"false\"", "\"F\"", "\"False\""):
                            # pars[3] = False
                            exp_rate_parameterization = False

            # making sure essential parameters of distribution have been specified
            if not exp_scale_or_rate:
                raise ec.ParseMissingParameterError("rate")

            # return dnpar.DnExponential(pars, parent_node_tracker)
            return dnpar.DnExponential(
                n_samples,
                n_repl,
                exp_scale_or_rate,
                exp_rate_parameterization,
                parent_node_tracker)

        elif dn_id == "gamma":
            gamma_shape: ty.List[float] = []
            gamma_scale_or_rate: ty.List[float] = []
            gamma_rate_parameterization: bool = False

            # { mean: node_dag1_name, sd: node_dag2_name, ... }
            if dn_param_dict:
                # val is list
                for arg, val in dn_param_dict.items():

                    # needed for building inference specifications
                    if isinstance(val[0], pgm.StochasticNodeDAG):
                        parent_node_tracker[arg] = val[0].node_name

                    extracted_val_list = pgm.extract_vals_as_str_from_node_dag(val)

                    if not cls.grammar_check("gamma", arg):
                        raise ec.ParseNotAParameterError(arg)

                    elif arg == "n":
                        if len(extracted_val_list) > 1:
                            raise ec.ParseRequireSingleValueError(
                                dnpar.DnGamma.DN_NAME, arg)

                        try:
                            n_samples = int(extracted_val_list[0])

                        except ValueError:
                            raise ec.ParseRequireIntegerError(
                                dnpar.DnGamma.DN_NAME, arg)

                    elif arg == "nr":
                        if len(extracted_val_list) > 1:
                            raise ec.ParseRequireSingleValueError(
                                dnpar.DnGamma.DN_NAME, arg)
                        
                        try:
                            n_repl = int(extracted_val_list[0])

                        except ValueError:
                            raise ec.ParseRequireIntegerError(
                                dnpar.DN_NAME, arg)

                    elif arg == "shape":
                        try:
                            gamma_shape = [float(v) for v in extracted_val_list]

                        except ValueError:
                            raise ec.ParseRequireNumericError(
                                dnpar.DN_NAME, arg)

                    elif arg == "scale":
                        try:
                            gamma_scale_or_rate = \
                                [float(v) for v in extracted_val_list]

                        except ValueError:
                            raise ec.ParseRequireNumericError(
                                dnpar.DN_NAME, arg)

                    elif arg == "rate_parameterization":
                        if extracted_val_list[0] in \
                            ("\"true\"", "\"T\"", "\"True\""):
                            gamma_rate_parameterization = True

                        elif extracted_val_list[0] in \
                            ("\"false\"", "\"F\"", "\"False\""):
                            gamma_rate_parameterization = False

            # making sure essential parameters of distribution have been specified
            for par_obj, par_name in \
                    ((gamma_shape, "shape"), (gamma_scale_or_rate, "scale")):
                if not par_obj:
                    raise ec.ParseMissingParameterError(par_name)

            # return dnpar.DnGamma(pars, parent_node_tracker)
            return dnpar.DnGamma(
                n_samples,
                n_repl,
                gamma_shape,
                gamma_scale_or_rate,
                gamma_rate_parameterization,
                parent_node_tracker)

        elif dn_id == "unif":
            unif_min: ty.List[float] = list()
            unif_max: ty.List[float] = list()

            # { mean: node_dag1_name, sd: node_dag2_name, ... }
            if dn_param_dict:
                # val is list
                for arg, val in dn_param_dict.items():
                    if isinstance(val[0], pgm.StochasticNodeDAG):
                        # needed for building inference specifications
                        parent_node_tracker[arg] = val[0].node_name

                    extracted_val_list = pgm.extract_vals_as_str_from_node_dag(val)

                    if not cls.grammar_check("unif", arg):
                        raise ec.ParseNotAParameterError(arg)

                    elif arg == "n":
                        if len(extracted_val_list) > 1:
                            raise ec.ParseRequireSingleValueError(
                                dnpar.DnUnif.DN_NAME, arg)

                        try:
                            n_samples = int(extracted_val_list[0])

                        except ValueError:
                            raise ec.ParseRequireIntegerError(
                                dnpar.DnUnif.DN_NAME, arg)

                    elif arg == "nr":
                        if len(extracted_val_list) > 1:
                            raise ec.ParseRequireSingleValueError(
                                dnpar.DnUnif.DN_NAME, arg)
                        
                        try:
                            n_repl = int(extracted_val_list[0])

                        except ValueError:
                            raise ec.ParseRequireIntegerError(
                                dnpar.DnUnif.DN_NAME, arg)

                    elif arg == "min":
                        try:
                            unif_min = [float(v) for v in extracted_val_list]

                        except ValueError:
                            raise ec.ParseRequireNumericError(
                                dnpar.DnUnif.DN_NAME, arg)

                    elif arg == "max":
                        try:
                            unif_max = [float(v) for v in extracted_val_list]

                        except ValueError:
                            raise ec.ParseRequireNumericError(
                                dnpar.DnUnif.DN_NAME, arg)

            # making sure essential parameters of distribution have been specified
            for par_obj, par_name in \
                ((unif_min, "min"), (unif_max, "max")):
                if not par_obj:
                    raise ec.ParseMissingParameterError(par_name)

            return dnpar.DnUnif(
                n_samples,
                n_repl,
                unif_min,
                unif_max,
                parent_node_tracker)

        elif dn_id == "categorical":
            n_cats: int = 0
            cats: ty.List[int] = list()

            n_probs: int = 0
            probs: ty.List[float] = list()

            if dn_param_dict:
                # val is list
                for arg, val in dn_param_dict.items():
                    if isinstance(val[0], pgm.StochasticNodeDAG):
                        # needed for building inference specifications
                        parent_node_tracker[arg] = val[0].node_name

                    # extracting values and removing empty values
                    extracted_val_list = \
                        [v for v in pgm.extract_vals_as_str_from_node_dag(val) \
                         if v]

                    if not cls.grammar_check("categorical", arg):
                        raise ec.ParseNotAParameterError(arg)

                    elif arg == "n":
                        if len(extracted_val_list) > 1:
                            raise ec.ParseRequireSingleValueError(
                                dnpar.DnCategorical.DN_NAME, arg)

                        try:
                            n_samples = int(extracted_val_list[0])

                        except ValueError:
                            raise ec.ParseRequireIntegerError(
                                dnpar.DnCategorical.DN_NAME, arg)

                    elif arg == "nr":
                        if len(extracted_val_list) > 1:
                            raise ec.ParseRequireSingleValueError(
                                dnpar.DnCategorical.DN_NAME, arg)

                        try:
                            n_repl = int(extracted_val_list[0])

                        except ValueError:
                            raise ec.ParseRequireIntegerError(
                                dnpar.DnCategorical.DN_NAME, arg)

                    elif arg == "cats":
                        n_cats = len(extracted_val_list)

                        # check 1:
                        # it only makes sense to use this distribution
                        # if there are 2+ categories
                        if n_cats in (0, 1):
                            raise ec.IncorrectDimensionError("\'cats\'",
                                                             n_cats,
                                                             "> 1")

                        # check 2: categories must get integer IDs
                        try:
                            cats = [int(v) for v in extracted_val_list]

                        except ValueError:
                            raise ec.ParseRequireIntegerError(
                                dnpar.DnCategorical.DN_NAME, arg)

                    elif arg == "probs":
                        n_probs = len(extracted_val_list)

                        # check 1:
                        # it only makes sense to use this distribution
                        # if there are 2+ categories (meaning 2+ probs)
                        if n_probs in (0, 1):
                            raise ec.IncorrectDimensionError("\'probs\'",
                                                             n_probs,
                                                             "> 1")

                        # check 2: there must be 1 prob per cat
                        elif n_probs != n_cats:
                            raise ec.IncorrectDimensionError("\'probs\'",
                                                             n_probs,
                                                             "# of categories (" +
                                                             str(n_cats) +
                                                             ")")

                        # check 3: values must be floats
                        try:
                            probs = [float(v) for v in extracted_val_list]

                        except ValueError:
                            raise ec.ParseRequireNumericError(
                                dnpar.DnCategorical.DN_NAME, arg)

                        # check 4: floats must add up to 1.0
                        prob_sum = sum(probs)
                        if prob_sum > 1.0 or (1.0 - prob_sum) > 1e-3:
                            arg_str = "[" + ", ".join(str(v) for v in val) + "]"
                            raise ec.ParseInvalidArgumentError(
                                "probs",
                                arg_str,
                                ("Sum of probabilities must be within "
                                "1e-3 of unity.")
                            )

            # making sure essential parameters of distribution have been specified
            for par_obj, par_name in \
                    ((cats, "cats"), (probs, "probs")):
                if not par_obj:
                    raise ec.ParseMissingParameterError(par_name)

            # this distribution requires a vector parameter, which could in
            # principle be vectorized itself (i.e., be a 2D parameter)
            cats = [cats]
            probs = [probs]

            return dnpar.DnCategorical(
                n_samples,
                n_repl,
                cats,
                probs,
                parent_node_tracker)


    @classmethod
    def init_return_discrete_SSE_dn(
        cls,
        dn_param_dict: ty.Dict[str, ty.List[ty.Union[str, pgm.NodeDAG]]]) \
            -> pgm.DistrForSampling:
        """Create and return SSE distribution for sampling.

        Args:
            dn_param_dict (dict): Dictionary with distribution
                parameter names as keys, and their values as values.

        Returns:
            DistributionDAG: a DistributionDAG instance.
        """

        if not dn_param_dict:
            raise ec.ParseMissingSpecificationError("discrete_sse")
        
        for arg, val in dn_param_dict.items():
            if not cls.grammar_check("discrete_sse", arg):
                raise ec.ParseNotAParameterError(arg)

        return make_dnsse.make_discrete_SSE_dn("discrete_sse", dn_param_dict)

    @classmethod
    def create_dn_obj(
        cls,
        dn_id: str,
        dn_param_dict: ty.Dict[str, ty.List[ty.Union[str, pgm.NodeDAG]]]) \
            -> pgm.DistrForSampling:
        """Create and return prob. distribution (for sampling) object.

        Args:
            dn_id (str): Name of the sampling distribution to create.
            dn_param_dict (dict): Dictionary containing distribution
                parameter names (str) as keys and lists (of either
                strings or StochasticNodeDAGs) as values.

        Returns:
            DistributionDAG: a DistributionDAG instance.
        """

        #############################
        #  Parametric distributions #
        #############################

        if dn_id == "quasse":
            for name in dn_param_dict:
                if not cls.grammar_check(dn_id, name):
                    raise ec.ParseNotAParameterError(name)
            return make_quasse(dn_param_dict)

        if dn_id in ("lognormal",
                     "normal",
                     "exponential",
                     "gamma",
                     "unif",
                     "categorical"):
            return cls.init_return_parametric_dn(dn_id, dn_param_dict)

        #################################
        #  Non-parametric distributions #
        #################################

        # discrete_sse
        elif dn_id == "discrete_sse":
            # dn_param_dict is validated inside
            return cls.init_return_discrete_SSE_dn(dn_param_dict)
