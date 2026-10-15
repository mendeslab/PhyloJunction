import typing as ty

# pj imports
import phylojunction.calculation.discrete_sse as sseobj
import phylojunction.utility.exception_classes as ec
import phylojunction.interface.grammar.det_fn_discrete_sse_makers as detsse
import phylojunction.interface.grammar.det_fn_map_attribute as detatt
import phylojunction.interface.grammar.det_fn_clock_makers as detclock
import phylojunction.interface.grammar.det_fn_branch_count_makers as detbrcount
import phylojunction.pgm.pgm as pgm
import phylojunction.data.tree as pjtr

__author__ = "Fabio K. Mendes"
__email__ = "fmendes@lsu.edu"


class PJDetFnGrammar():

    det_fn_grammar_dict: ty.Dict[str, ty.Set[str]]

    #########################################
    # All available deterministic functions #
    #########################################
    det_fn_grammar_dict = {
        "sse_prob": set(["name", "value", "state", "epoch"]),
        "sse_rate": set(["name", "value", "event", "states", "epoch"]),
        "sse_stash": set(["flat_rate_mat", "flat_prob_mat", "n_states",
                          "seed_age", "epoch_age_ends", "n_epochs"]),
        "map_attr": set(["tree", "fun", "tip_attr_file_path",
                         "maps_file_path", "attr_name", "n_regions",
                         "geo"]),
        "gamma_spike_clock": set(["tree", "attr_name", "clock_rate",
                                  "br_specific_rel_rate",
                                  "spike_mean", "spike_shape",
                                  "strict_clock_when_no_stubs"]),
        "n_rec_branches": set(["tree"])
    }

    @classmethod
    def grammar_check(cls, det_fn_id: str, fn_param: str) -> bool:
        if fn_param in cls.det_fn_grammar_dict[det_fn_id]:
            return True
        return False

    @classmethod
    def init_return_state_dep_rate(
        cls,
        det_fn_param_dict:
            ty.Dict[str, ty.List[ty.Union[str, pgm.NodeDAG]]]) \
            -> sseobj.DiscreteStateDependentRate:

        try:
            return detsse.make_DiscreteStateDependentRate(
                "sse_rate", det_fn_param_dict)

        except (ec.ParseMissingSpecificationError,
                ec.ParseMissingArgumentError,
                ec.ParseFunctionArgError,
                ec.ParseRequireSingleValueError,
                ec.ParseRequireIntegerError) as e:
            raise ec.ParseDetFnInitFailError("sse_rate", e.message)

    @classmethod
    def init_return_state_dep_prob(
        cls,
        det_fn_param_dict:
            ty.Dict[str, ty.List[ty.Union[str, pgm.NodeDAG]]]) \
            -> sseobj.DiscreteStateDependentProbability:

        try:
            return detsse.make_DiscreteStateDependentProbability(
                "sse_prob",
                det_fn_param_dict)

        except (ec.ParseMissingArgumentError,
                ec.ParseRequireSingleValueError,
                ec.ParseRequireIntegerError) as e:
            raise ec.ParseDetFnInitFailError("sse_prob", e.message)

    @classmethod
    def init_return_sse_stash(
        cls,
        det_fn_param_dict:
            ty.Dict[str, ty.List[ty.Union[str, pgm.NodeDAG]]]) \
            -> sseobj.SSEStash:

        try:
            return detsse.make_SSEStash("sse_stash", det_fn_param_dict)

        except (ec.ParseRequireSingleValueError,
                ec.ParseRequireIntegerError,
                ec.ParseRequireNumericError,
                ec.IncorrectDimensionError,
                ec.ParseMissingParameterError,
                ec.ParseInvalidArgumentError) as e:
            raise ec.ParseDetFnInitFailError("sse_stash", e.message)

    @classmethod
    def init_return_tree_mapped_with_attr(
        cls,
        det_fn_param_dict: \
            ty.Dict[str, ty.List[ty.Union[str, pgm.NodeDAG]]]) \
                -> ty.List[pjtr.AnnotatedTree]:

        try:
            return detatt.make_mapped_ann_tree("map_attr", det_fn_param_dict)

        except (ec.ParseInvalidArgumentError,
                ec.ParseRequireIntegerError,
                ec.PJIOFileDoesNotExistError) as e:
            raise ec.ParseDetFnInitFailError("map_attr", e.message)

    @classmethod
    def init_return_gamma_spike_clock_tree(
        cls,
        det_fn_param_dict:
            ty.Dict[str, ty.List[ty.Union[str, pgm.NodeDAG]]]) \
            -> ty.List[pjtr.AnnotatedTree]:

        try:
            return detclock.make_gamma_spike_clock_tree(
                "gamma_spike_clock", det_fn_param_dict)

        except (ec.ParseMissingParameterError,
                ec.ParseMissingArgumentError,
                ec.ParseInvalidArgumentError,
                ec.ParseRequireSingleValueError,
                ec.ParseRequireNumericError,
                ValueError) as e:
            raise ec.ParseDetFnInitFailError("gamma_spike_clock", str(e))

    @classmethod
    def init_return_n_rec_branches(
        cls,
        det_fn_param_dict:
            ty.Dict[str, ty.List[ty.Union[str, pgm.NodeDAG]]]) \
            -> ty.List[int]:

        try:
            return detbrcount.make_n_rec_branches(
                "n_rec_branches", det_fn_param_dict)

        except (ec.ParseMissingParameterError,
                ec.ParseMissingArgumentError,
                ec.ParseInvalidArgumentError) as e:
            raise ec.ParseDetFnInitFailError("n_rec_branches", str(e))

    @classmethod
    def create_det_fn_obj(
        cls,
        det_fn_id: str,
        det_fn_param_dict: ty.Dict[str, ty.List[ty.Union[str, pgm.NodeDAG]]]) \
            -> ty.Optional[
                ty.Union[
                    sseobj.DiscreteStateDependentRate,
                    sseobj.DiscreteStateDependentProbability,
                    sseobj.SSEStash,
                    ty.List[pjtr.AnnotatedTree],
                    ty.List[int]
                ]]:
        """
        Build and return deterministic function object.

        Args:
            det_fn_id (str): Name of deterministic function being called
            det_fn_param_dict (dict): Dictionary containing deterministic
                function parameter names (str) as keys and lists (of either
                strings or NodeDAGs) as values

        Returns:
            Object: one of a variety of objects containing information for
                later building a sampling distribution
        """

        if not det_fn_param_dict:
            raise ec.ParseMissingSpecificationError(det_fn_id)

        for arg in det_fn_param_dict:
            if not cls.grammar_check(det_fn_id, arg):
                raise ec.ParseNotAParameterError(arg)
            
        # validate input
        if det_fn_id == "sse_rate":
            # health checks inside
            return cls.init_return_state_dep_rate(det_fn_param_dict)

        if det_fn_id == "sse_prob":
            # health checks inside
            return cls.init_return_state_dep_prob(det_fn_param_dict)

        if det_fn_id == "sse_stash":
            return cls.init_return_sse_stash(det_fn_param_dict)
        
        if det_fn_id == "map_attr":
            return cls.init_return_tree_mapped_with_attr(det_fn_param_dict)

        if det_fn_id == "gamma_spike_clock":
            return cls.init_return_gamma_spike_clock_tree(det_fn_param_dict)

        if det_fn_id == "n_rec_branches":
            return cls.init_return_n_rec_branches(det_fn_param_dict)

        return None


if __name__ == "__main__":
    # can be called from interface/grammar/
    # $ python3 det_fn_grammar.py
    #
    # can also be called from phylojunction/
    # $ python3 interface/grammar/det_fn_grammar.py
    # or
    # $ python3 -m interface.grammar.det_fn_grammar
    #
    # can also be called from VS Code, if open folder is phylojuction/

    print(PJDetFnGrammar.grammar_check("sse_rate", "name"))  # True
    print(PJDetFnGrammar.grammar_check("sse_rate", "cloud"))  # False

    # TODO: add code to run init_return_state_dep_rate and init_return_sse_stash
