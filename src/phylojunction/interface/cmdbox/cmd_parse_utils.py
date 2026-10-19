import sys
import typing as ty
import ast
import io
import tokenize
from phylojunction.interface.cmdbox.expressions import parse_expression
import re

# pj imports
import phylojunction.pgm.pgm as pgm
import phylojunction.utility.exception_classes as ec

__author__ = "Fabio K. Mendes"
__email__ = "fmendes@lsu.edu"


# note that in string literals, backslash "\" has a special
# meaning in that it escapes characters; importantly in string
# literals, some characters cannot be escaped, like "s" or "\n",
# i.e., "\s" does not escape "s"
# 
# if one wants to match "\s" or other non-escapable characters,
# e.g., "[", either one must use double backslashes, or make the
# string raw with "r"
cannot_start_with_this_regex = re.compile(r"[\(\[\{|^&+\-%*/=!><_~0-9]")
cannot_end_with_this_regex = re.compile(r"[\(\[\{|^&+\-%*/=!><_~]$")
whitespace_regex = re.compile(r"\s+")
int_or_float_regex = re.compile(r"(\d+(?:\.\d+)?)")
vector_value_regex = re.compile(r"[\[\]]")

# assignment
# Would two '-' make into this regex?
assign_regex = re.compile(r"\s*(<-)\s*")

character_value_regex = re.compile(r"\s*([a-zA-Z]+)\s*")
quoted_character_value_regex = re.compile(r"\s*(\"[a-zA-Z]+\")\s*")

# sampling dn
sampled_as_regex = re.compile(r"\s*(~)\s*")
sampling_dn_spec_regex = re.compile(r"([a-zA-Z][a-zA-Z_]*)\((.+)\)")

# det functions
deterministic_regex = re.compile(r"\s*(:=)\s*")


def val_or_obj(dag_obj: pgm.DirectedAcyclicGraph,
               val: ty.List[str]) -> ty.List[ty.Union[pgm.NodeDAG, str]]:
    """Return list of strings with values or node names.

    Checks if provided values are directly accessible as values
    (e.g., 1.0, \"a_string_in_quotes\")) or if they are names of nodes
    potentially in the graphical model. If the latter, check if they
    are indeed nodes in the graphical model, and if so, append to
    return

    Args:
        dag_obj (DirectedAcyclicGraph): DAG object being interrogated.
        val (str): List of strings, one for each argument coming
            from a parsed command.

    Raises:
        ec.InexistentVariableError: _description_

    Returns:
        (str): List of strings (representing values or quoted-enclosed
            strings) and/or stochastic node objects.
    """

    val_or_obj_list: ty.List[ty.Union[pgm.NodeDAG, str]] = list()

    for v in val:
        if isinstance(v, str):
            # checking if string could potentially be a node object that
            # has a name
            if re.match(character_value_regex, v):

                # if it does find a node with that name, we add the node object
                try:
                    # appending StochasticNodeDAG
                    val_or_obj_list.append(dag_obj.name_node_dict[v])

                except KeyError:
                    raise ec.InexistentVariableError(v)

            # if string is instead number, or a word enclosed in quotes,
            # we add the string
            else:
                val_or_obj_list.append(v)

    return val_or_obj_list


def parse_spec(
        dag_obj: pgm.DirectedAcyclicGraph,
        fn_spec_str: str,
        cmd_line: str) \
        -> ty.Tuple[
            ty.Dict[str, ty.List[ty.Union[str, pgm.NodeDAG]]],
            ty.List[pgm.NodeDAG]]:

    spec_dict: ty.Dict[str, str] = tokenize_fn_spec(fn_spec_str, cmd_line)
    spec_dict_return: \
        ty.Dict[str, ty.List[ty.Union[str, pgm.NodeDAG]]] = dict()

    parent_pgm_nodes: ty.List[pgm.NodeDAG] = []
    for param_name, an_arg in spec_dict.items():

        #############
        # IMPORTANT #
        #############
        # arguments of dn/det parameters will come out as lists
        arg_list = list()

        # if argument is a list
        if re.match(vector_value_regex, an_arg):
            arg_list = parse_val_vector(an_arg)

        # if scalar variable
        else:
            arg_list.append(an_arg)

        val_obj_list = val_or_obj(dag_obj, arg_list)
        spec_dict_return[param_name] = val_obj_list
        # { param_name_str:
        #       [ number_or_quoted_str_str1,
        #         a_NodeDAG1,
        #         number_or_quoted_str_str2,
        #         aNodeDAG2
        #       ]
        # }

        for vo in val_obj_list:
            if isinstance(vo, pgm.NodeDAG):
                parent_pgm_nodes.append(vo)

    # values in spec_dict will be lists of strings or lists of NodeDAGs
    return spec_dict_return, parent_pgm_nodes


def parse_val_vector(vec_str: str) -> ty.List[str]:
    """Parse string specifying vector when a variable is assigned

    Args:
        vec_str (str): A string specifying a vector, e.g., "[1, 2, 3]"

    Raises:
        ec.ScriptSyntaxError: If there are more than two squared
            brackets, this is thrown

    Returns:
        (str): Vector of strings, each being one of the values inside
            the vector being specified by \'vec_str\'.
    """

    if len(re.findall(vector_value_regex, vec_str)) > 2:
        raise ec.ScriptSyntaxError(
            vec_str,
            ("Something went wrong during variable assignment. If a vector "
             "of values is specified, there can only be one left and one "
             "right squared brackets. Exiting."))

    # ok!
    else:
        scalar_values_str = \
            re.sub(whitespace_regex, "",
                   re.sub(vector_value_regex, "", vec_str))
        return scalar_values_str.split(",")


def tokenize_fn_spec(fn_spec_str: str, cmd_line: str) -> ty.Dict[str, str]:
    """Tokenize string containing function specification into dictionary

    Args:
        fn_spec (str): String containing the specification of function
            (i.e., what appears between \'(\' and \')\'.

    Returns:
        (dict): Dictionary containing function parameter strings as
            keys, argument strings as values (e.g., {"mean": "0.0"}).
    """

    # Transitional string arguments for the existing factories. Structural parsing handles
    # nested delimiters and quoted text; remove this conversion when factories accept values.
    try:
        source = "constructor(" + fn_spec_str + ")"
        call = parse_expression(source)
        if call.args:
            raise ValueError("Use named arguments in this legacy entry point.")
        result = {}
        for keyword in call.keywords:
            text = ast.get_source_segment(source, keyword.value)
            # Strip formatting whitespace without changing whitespace inside quoted strings.
            result[keyword.arg] = ''.join(token.string for token in
                tokenize.generate_tokens(io.StringIO(text).readline))
        return result
    except (SyntaxError, ValueError) as error:
        raise ec.ScriptSyntaxError(cmd_line, str(error)) from error


if __name__ == "__main__":
    # can be called from interface/
    # $ python3 cmd_parse_utils.py
    #
    # can also be called from phylojunction/
    # $ python3 interface/cmd_parse_utils.py
    # or
    # $ python3 -m interface.cmd_parse_utils
    #
    # can also be called from VS Code, if open folder is phylojuction/

    fn_spec1 = \
        "name=\"lambda0\", value=1.0, event=\"w_speciation\", " + \
        "states=[0, 0, 0]"
    cmd_line1 = "a <- sse_rate(" + fn_spec1 + ")"
    token_dict1 = tokenize_fn_spec(fn_spec1, cmd_line1)

    print(token_dict1)
