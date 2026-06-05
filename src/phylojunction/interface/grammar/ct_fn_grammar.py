import typing as ty

# pj imports
import phylojunction.utility.exception_classes as ec
import phylojunction.pgm.pgm as pgm
import phylojunction.interface.grammar.ct_fn_treereader_makers as make_tree

__author__ = "Fabio K. Mendes"
__email__ = "fmendes@lsu.edu"


class PJCtFnGrammar():

    ct_fn_grammar_dict: ty.Dict[str, ty.Set[str]]

    #########################################
    # All available deterministic functions #
    #########################################
    ct_fn_grammar_dict = {
        "read_tree": set(["n", "nr", "file_path", "string", "node_name_attr", "epoch_age_ends"])
    }

    @classmethod
    def grammar_check(cls, ct_fn_id: str, fn_param: str) -> bool:
        if fn_param in cls.ct_fn_grammar_dict[ct_fn_id]:
            return True
        return False

    @classmethod
    def init_return_ann_tr(
        cls,
        ct_fn_param_dict:
            ty.Dict[str, ty.List[ty.Union[str, pgm.NodeDAG]]]) \
            -> pgm.ConstantFn:
        
        if not ct_fn_param_dict:
            raise ec.ParseMissingSpecificationError("read_tree")
        
        for arg in ct_fn_param_dict:
            if not cls.grammar_check("read_tree", arg):
                raise ec.ParseNotAParameterError(arg)

        # debugging    
        # print(ct_fn_param_dict)

        return make_tree.make_tree_reader("read_tree", ct_fn_param_dict)
        
    @classmethod
    def create_ct_fn_obj(
        cls,
        ct_fn_id: str,
        ct_fn_param_dict: ty.Dict[str, ty.List[ty.Union[str, pgm.NodeDAG]]]) \
            -> ty.Optional[
                # ty.Union[
                    pgm.ConstantFn
                # ]
                ]:
        """
        Build and return constant function object.

        Args:
            ct_fn_id (str): Name of constant function to being called
            ct_fn_param_dict (dict): Dictionary containing constant
                function parameter names (str) as keys and lists (of either
                strings or NodeDAGs) as values

        Returns:
            Object: one of a variety of objects to be stored within a clamped
            stochastic node
        """

        # validate input
        if ct_fn_id == "read_tree":
            # health checks inside
            return cls.init_return_ann_tr(ct_fn_param_dict)

        return None