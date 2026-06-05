import typing as ty
import enum

__author__ = "Fabio K. Mendes"
__email__ = "fmendes@lsu.edu"


class ScriptSyntaxError(Exception):
    cmd_line: str
    message: str

    def __init__(self, cmd_line: str, message: str) -> None:
        self.cmd_line = cmd_line
        self.message = "ERROR: " \
            + self.cmd_line \
            + ("\n\nThe line above had a syntax problem and could not "
               "be tokenized. ") \
            + message
        # super().__init__(self.message)

    def __str__(self) -> str:
        return self.message


class InexistentVariableError(Exception):
    rv_name: str
    message: str

    def __init__(self, rv_name: str) -> None:
        self.rv_name = rv_name
        self.message = ""
        super().__init__(self.message)

    def __str__(self) -> str:
        return "ERROR: There was a problem when assigning \'" \
            + self.rv_name + "\' as an argument. It seems this " \
            + "variable has not been not initialized."


class VariableAssignmentError(Exception):
    rv_name: str
    message: str

    def __init__(self, rv_name: str, message: str = "") -> None:
        self.message = message
        self.rv_name = rv_name
        super().__init__(self.message)

    def __str__(self) -> str:
        return "ERROR: Cannot assign values of different types to a " \
            + "variable (\'" + self.rv_name + "\'). Exiting..."


class VariableMisspec(Exception):
    rv_name: str
    message: str

    def __init__(self, rv_name: str) -> None:
        self.message = ""
        self.rv_name = rv_name
        super().__init__(self.message)

    def __str__(self) -> str:
        return "ERROR: Could not parse value out of variable \'" \
            + self.rv_name + "\'."


class FunctionArgsMismatchError(Exception):
    fn_name: str
    message: str

    def __init__(self, fn_name: str, message: str) -> None:
        self.fn_name = fn_name
        self.message = message
        super().__init__(self.message)

    def __str__(self) -> str:
        return "ERROR: Could not call function " + self.fn_name \
            + " because argument values are incongruent in some way. " \
            + self.message


class NotBetweenZeroAndOneError(Exception):
    par_name: str
    message: str

    def __init__(self, par_name: str, message: str = "") -> None:
        self.par_name = par_name
        self.message = message
        super().__init__(self.message)

    def __str__(self) -> str:
        return "ERROR: \'" + self.par_name + "\' must be contained in " \
            + "[0,1], but one or more of its values was " + self.message \
            + "."


class IncorrectDimensionError(Exception):
    container_name: str
    obs_len: int
    exp_len: int
    message: str

    def __init__(self,
                 container_name: str,
                 obs_len: int,
                 exp_len: int = 0) -> None:
        self.obj_name = container_name
        self.obs_len = obs_len
        self.exp_len = exp_len
        self.message = \
            "Incorrect dimension of container " + self.obj_name \
            + ", which was of size " + str(obs_len) + ". "

        if exp_len != 0:
            self.message += "The expected dimension was " + str(exp_len) \
                + "."

        super().__init__(self.message)

    def __str__(self) -> str:
        return self.message


class DimensionalityError(Exception):
    dn_name: str
    par_name: str
    message: str

    def __init__(self, dn_name: str, par_name: ty.Optional[str] = "") -> None:
        self.dn_name = dn_name
        self.par_name = par_name
        self.message = \
            "ERROR: When specifying distribution " + self.dn_name
        if self.par_name:
            self.message += " (parameter " + self.par_name + ")"
        self.message += \
            ":\n (i)   the number of values provided for the parameter(s) " \
            + "differed, or\n (ii)  the number(s) was larger than the " \
            + "specified number of samples (i.e., draws), or\n (iii) the " \
            + "number(s) was larger than 1, but smaller than the specified " \
            + "number of samples."

        super().__init__(self.message)

    def __str__(self) -> str:
        return self.message


class NodeInferenceDimensionalityError(Exception):
    node_name: str
    message: str

    def __init__(self, node_name: str) -> None:
        self.node_name = node_name
        self.message = ""
        super().__init__(self.message)

    def __str__(self) -> str:
        return "ERROR: When converting values of node " + self.node_name \
            + " to inference script, the number of values as larger than " \
            + "1, but smaller than the number of samples specified for " \
            + "other nodes. Unsure how to parse values."


class NoPlatingAllowedError(Exception):
    det_name: str
    message: str
    node_dag_name: str

    def __init__(self,
                 det_name: str,
                 problematic_node_dag_name: str,
                 message: str = "") -> None:
        self.det_name = det_name
        self.message = message
        self.node_dag_name = problematic_node_dag_name
        super().__init__(self.message)

    def __str__(self) -> str:
        return "ERROR: When executing " + self.det_name + "(), replicates " \
            + "were detected for argument " + self.node_dag_name \
            + ". Plating is not supported for this deterministic function."


class ReplicateNumberError(Exception):
    node_name: str
    message: str

    def __init__(self, node_name, message="") -> None:
        self.node_name = node_name
        self.message = message
        super().__init__(self.message)

    def __str__(self) -> str:
        return "ERROR: When going through values in " \
            + self.node_name + ", the number of replicates " \
            + "differed among simulations. One or more of the " \
            + "simulations failed to generate the specified " \
            + "number of replicates. Exiting..."


class DimensionalityWarning(Exception):
    rv_name: str
    dn_name: str
    message: str

    def __init__(self, rv_name: str, dn_name: str, message: str = "") -> None:
        self.rv_name = rv_name
        self.dn_name = dn_name
        self.message = message
        super().__init__(self.message)

    def __str__(self) -> str:
        return "Warning: Distribution " + self.dn_name \
            + " will be called to sample " + self.rv_name \
            + " using vectorization."


class StateDependentParameterMisspec(Exception):
    message: str

    def __init__(self, message: str = "") -> None:
        self.message = message
        super().__init__(self.message)

    def __str__(self) -> str:
        return "ERROR: Misspecified state-dependent parameter input. " \
            + self.message


class SSEStashMisspec(Exception):
    message: str

    def __init__(self, message: str) -> None:
        self.message = message
        super().__init__(self.message)

    def __str__(self) -> str:
        return "ERROR: Misspecified SSE stash parameter input. " + self.message


class InvalidMCMCChainLength(Exception):
    message: str

    def __init__(self, message: str) -> None:
        self.message = message
        super().__init__(self.message)

    def __str__(self) -> str:
        return "ERROR: " + self.message


# Initialization errors
class ObjInitInvalidArgError(Exception):
    obj_name: str
    message: str

    def __init__(self, obj_name: str, par_name: str, message: ty.Optional[str] = "") -> None:
        self.obj_name = obj_name
        self.message = "Could not initialize the object of \'" + obj_name + \
            "\'. The argument provided to " + par_name + " was invalid."

        if message:
            self.message += " " + message

        super().__init__(self.message)

    def __str__(self) -> str:
        return self.message


class ObjInitMissingParameterError(Exception):
    obj_name: str
    message: str

    def __init__(self, obj_name: str, par_name: str, message: ty.Optional[str] = "") -> None:
        self.obj_name = obj_name
        self.message = "Could not initialize the object of " + obj_name + \
            " . Parameter " + par_name + " was not provided an argument."
        if message:
            self.message += " " + message

        super().__init__(self.message)

    def __str__(self) -> str:
        return self.message


class ObjInitIncorrectDimensionError(Exception):
    obj_name: str
    message: str
    obs_len: int
    exp_len: int
    at_least: bool

    def __init__(self,
                 obj_name: str,
                 container_name: str,
                 obs_len: int,
                 exp_len: int = 0,
                 at_least: bool = False) -> None:
        self.obj_name = obj_name
        self.obs_len = obs_len
        self.exp_len = exp_len
        self.message = "Could not initialize the object of " + obj_name \
            + ". Incorrect dimension of container " + container_name \
            + ", which was of size " + str(obs_len) + ". "

        if exp_len != 0:
            if at_least:
                self.message += "The expected dimension was at least " + str(exp_len) \
                    + "."

            else:
                self.message += "The expected dimension was " + str(exp_len) \
                    + "."

        super().__init__(self.message)

    def __str__(self) -> str:
        return self.message


class ObjInitRequireSameParameterTypeError(Exception):
    message: str
    obj_name: str
    n_diff_par: int

    def __init__(self, obj_name: str, n_diff_par: int) -> None:
        self.obj_name = obj_name
        self.n_diff_par = n_diff_par
        self.message = "When specifying object " + self.obj_name \
            + " only one type of parameter is allowed. Found " \
            + str(self.n_diff_par) + "."

        super().__init__(self.message)

    def __str__(self) -> str:
        return self.message


class ObjInitMissingStateDependentParameterError(Exception):
    symmetric_diff_str: str
    message: str

    def __init__(self,
                 epoch_missing_param: int,
                 symmetric_diff_set: ty.Set[ty.Any]) -> None:
        self.symmetric_diff_str = "( "

        for i in symmetric_diff_set:
            ith_str = ""

            if isinstance(i, int):
                ith_str = str(i)

            elif isinstance(i, tuple):
                ith_str = "(" + ", ".join(str(j) for j in i) + ")"

            if self.symmetric_diff_str == "( ":
                self.symmetric_diff_str += ith_str

            else:
                self.symmetric_diff_str += ", " + ith_str

        self.symmetric_diff_str += " )."

        self.message = \
            ("One or more state-dependent parameters were missing "
             "at least from time slice ") + str(epoch_missing_param) \
            + ". The symmetric difference between time-slice state " \
            + "sets is " + self.symmetric_diff_str

        super().__init__(self.message)

    def __str__(self) -> str:
        return self.message


class ObjInitRepeatedStateDependentParameterError(Exception):
    message: str

    def __init__(self,
                 epoch_w_repeated_param: int,
                 repeated_param: \
                    ty.Tuple[ty.Tuple[int], enum.Enum]) -> None:
        repeated_state_str = ""

        if isinstance(repeated_param[0], int):
            repeated_state_str = str(repeated_param[0])

        elif isinstance(repeated_param[0], tuple):
            repeated_state_str = "("

            for i in repeated_param[0]:
                ith_str = str(i)

                if repeated_state_str == "(":
                    repeated_state_str += ith_str

                else:
                    repeated_state_str += ", " + ith_str

            repeated_state_str += ")"

        self.message = "State-dependent parameter defined by states "

        if len(repeated_param) == 1:
            self.message += repeated_state_str + " were repeated in epoch " \
                + str(epoch_w_repeated_param) + "."
            
        else:
            self.message += "(event: " + str(repeated_param[1]) \
                .replace("MacroevolEvent.", "") + ") " \
                + repeated_state_str + " were repeated in epoch " \
                + str(epoch_w_repeated_param) + "."

        super().__init__(self.message)

    def __str__(self) -> str:
        return self.message


class ObjInitRequireNonZeroStateDependentParameterError(Exception):
    message: str
    obj_name: str

    def __init__(self, obj_name: str, dimension_idx: int) -> None:
        self.obj_name = obj_name
        self.message = "When specifying object " + self.obj_name \
            + ", one of its dimensions (" + str(dimension_idx) + (") "
            "had zero-valued state-dependent parameters. At least one "
            "non-zero parameter must be provided.")

        super().__init__(self.message)

    def __str__(self) -> str:
        return self.message


# Syntax errors: constant functions,
# sampling distribution,
# det. function exceptions
class ParseFunctionArgError(Exception):
    par_name: str
    message: str

    def __init__(self, par_name: str, message: str) -> None:
        self.par_name = par_name
        self.message = "Could not set value for function parameter " \
            + self.par_name + ". " + message
        super().__init__(self.message)

    def __str__(self) -> str:
        return self.message


class ParseNotAParameterError(Exception):
    par_name: str
    message: str

    def __init__(self, par_name: str) -> None:
        self.par_name = par_name
        self.message = "\'" + par_name + "\' is not a valid parameter."
        super().__init__(self.message)

    def __str__(self) -> str:
        return self.message


class ParseRequireSingleValueError(Exception):
    obj_name: str
    message: str

    def __init__(self, obj_name: str, arg: str) -> None:
        self.obj_name = obj_name
        self.arg = arg
        self.message = \
            "When specifying object " + obj_name + "'s parameter \'" \
            + self.arg + ("\', more than one value was "
                          "provided. A scalar is required.")

        super().__init__(self.message)

    def __str__(self) -> str:
        return self.message


class ParseRequireIntegerError(Exception):
    obj_name: str
    message: str

    def __init__(self, obj_name: str, arg: str) -> None:
        self.obj_name = obj_name
        self.message = \
            "When specifying object " + obj_name + "'s parameter \'" \
            + arg + ("\', something that could not be read as an integer "
                          "was provided. An integer is required.")
        super().__init__(self.message)

    def __str__(self) -> str:
        return self.message


class ParseRequirePositiveIntegerError(Exception):
    obj_name: str
    message: str

    def __init__(self, obj_name: str, arg: str) -> None:
        self.obj_name = obj_name
        self.message = \
            "When specifying object " + obj_name + "'s parameter \'" \
            + arg + ("\', something other than a positive integer was "
                          "provided. A positive integer is required.")
        super().__init__(self.message)

    def __str__(self) -> str:
        return self.message


class ParseRequireNumericError(Exception):
    obj_name: str
    message: str

    def __init__(self, obj_name: str, arg: str) -> None:
        self.obj_name = obj_name
        self.arg = arg
        self.message = \
            "When specifying object " + obj_name + "'s parameter \'" \
            + self.arg + ("\', something other than a number was "
                          "provided. A number is required.")
        super().__init__(self.message)

    def __str__(self) -> str:
        return self.message


class ParseMissingParameterError(Exception):
    par_name: str
    message: str

    def __init__(self, par_name: str) -> None:
        self.par_name = par_name
        self.message = "Parameter \'" + par_name + "\' is missing."
        super().__init__(self.message)

    def __str__(self) -> str:
        return self.message


class ParseMissingArgumentError(Exception):
    par_name: str
    message: str

    def __init__(self, par_name: str, arg_name: str = "") -> None:
        self.par_name = par_name
        self.message = "Argument "

        if arg_name:
            self.message += "\'" + arg_name + "\' "

        self.message += "for \'" + par_name + "\' is missing."

        super().__init__(self.message)

    def __str__(self) -> str:
        return self.message


class ParseInvalidArgumentError(Exception):
    par_name: str
    message: str

    def __init__(self, par_name: str, invalid_arg: str, message: str) -> None:
        self.par_name = par_name
        self.message = "Argument " + invalid_arg +  " is invalid " \
            + "for parameter \'" + self.par_name + "\'. " + message
        
        super().__init__(self.message)

    def __str__(self) -> str:
        return self.message


class ParseMissingSpecificationError(Exception):
    message: str

    def __init__(self, obj2spec_name: str) -> None:
        self.message = "Missing specification when calling \'" \
            + obj2spec_name + "\'."
        super().__init__(self.message)

    def __str__(self) -> str:
        return self.message


class ParseMutuallyExclusiveParametersError(Exception):
    message: str
    par_name: str
    mutually_exclusive_par_name: str

    def __init__(self,
                 par_name: str,
                 mutually_exclusive_par_name: str,
                 message: str = "") -> None:
        self.par_name = par_name
        self.mutually_exclusive_par_name = mutually_exclusive_par_name
        self.message = "Argument " + self.par_name \
            + " cannot be specified at the same time as \'" \
            + self.mutually_exclusive_par_name + "\'. " + message
        
        super().__init__(self.message)

    def __str__(self) -> str:
        return self.message


class ParsePathDoesNotExistError(Exception):
    message: str
    par_name: str
    path_str: str

    def __init__(self, par_name: str, path_str: str, message: ty.Optional[str] = "") -> None:
        super().__init__(par_name, path_str, message)
        self.par_name = par_name
        self.path_str = path_str
        self.message = "Path provided to argument \'" + self.par_name \
            + "\' does not seem to store a file." + message
        
    def __str__(self) -> str:
        return self.message


class ParseInvalidNewickStringError(Exception):
    message: str
    par_name: str

    def __init__(self,
                 par_name: str,
                 message: str = "") -> None:
        self.par_name = par_name
        self.message = "Newick string(s) provided in \'" + par_name \
            + "\' were invalid. " + message
        super().__init__(self.message)

    def __str__(self) -> str:
        return self.message


class ParseCtFnInitFailError(Exception):
    ct_fn_name: str
    message: str

    def __init__(self, ct_fn_name: str, message: str) -> None:
        super().__init__(ct_fn_name, message)
        self.ct_fn_name = ct_fn_name
        self.message = "Constant output from \'" + self.ct_fn_name \
            + "\' could not be instantiated. " + message

    def __str__(self) -> str:
        return self.message


class ParseDnInitFailError(Exception):
    dn_name: str
    message: str

    def __init__(self, dn_name: str, message: str) -> None:
        self.dn_name = dn_name
        self.message = "Parsing the specification of \'" + self.dn_name \
            + "\' failed. " + message
        super().__init__(self.message)

    def __str__(self) -> str:
        return self.message


class ParseDetFnInitFailError(Exception):
    det_fn_name: str
    message: str

    def __init__(self, det_fn_name: str, message: str) -> None:
        self.det_fn_name = det_fn_name
        self.message = "Deterministic output from \'" + det_fn_name \
            + "\' could not be instantiated. " + message
        super().__init__(self.message)

    def __str__(self) -> str:
        return self.message


# PGM exceptions #
class NodeDAGStatCantFloatError(Exception):
    message: str

    def __init__(self, node_name: str) -> None:
        self.message = "ERROR: When summarizing " + node_name \
            + " into string, a node statistic could not be float-ified"

    def __str__(self) -> str:
        return self.message


class DAGCannotInitialize(Exception):
    message: str

    def __init__(self, message: str) -> None:
        self.message = "ERROR: When initializing DAG object, " + message

    def __str__(self) -> str:
        return self.message
    

class DAGCannotAddNodeError(Exception):
    message: str

    def __init__(self, node_name: str, message: str) -> None:
        self.message = "Could not add \'" + node_name + "\' to DAG. " \
            + message

    def __str__(self) -> str:
        return self.message


# Tree exceptions #
class AnnotatedTreeMisspecError(Exception):
    message: str

    def __init__(self, message: str) -> None:
        self.message = ("Could not initialize AnnotatedTree because one of "
                        "its arguments was misspecified. ") + message
        super().__init__(self.message)

    def __str__(self) -> str:
        return self.message


class AnnotatedTreeIncorrectAnnotationError(Exception):
    message: str

    def __init__(self, message: str) -> None:
        self.message = ("Some aspect of the tree annotation was wrong. ") \
            + message
        super().__init__(self.message)

    def __str__(self) -> str:
        return self.message


class AnnotatedTreeNodeMissingAttrError(Exception):
    message: str

    def __init__(self, nd_name: str, attr_name: str, message: str) -> None:
        self.message = "Node " + nd_name + " did not have attribute " \
            + attr_name + ". " + message
        super().__init__(self.message)

    def __str__(self) -> str:
        return self.message


class AnnotatedTreeMissingNodeName(Exception):
    message: str

    def __init__(self) -> None:
        self.message = "Could not find the name of a node."
        super().__init__(self.message)
    
    def __str__(self) -> str:
        return self.message


# Generation exceptions #
class GenerateFailError(Exception):
    message: str

    def __init__(self, dn_name: str, message: str) -> None:
        self.message = "ERROR: Could not generate a sample from " \
            + dn_name + ". " + message

    def __str__(self) -> str:
        return self.message


class MissingColumnName(Exception):
    message: str

    def __init__(self, col_name: str, message: str) -> None:
        self.message = "ERROR: " + message + ". " + col_name \
            + " was not found in DataFrame."

    def __str__(self) -> str:
        return self.message


class RunTimeLimit(Exception):
    def __init__(self, runtime_limit: float) -> None:
        self.message = "ERROR: Hit runtime limit of " + str(runtime_limit) \
            + " seconds."

        super().__init__(self.message)

    def __str__(self) -> str:
        return self.message


class MaxNFailedAttemptsLimit(Exception):
    def __init__(self, max_n_attempts: int) -> None:
        self.message = \
            "ERROR: Hit maximum number of (failed) sampling attempts (" + \
            str(max_n_attempts) + ")."

        super().__init__(self.message)

    def __str__(self) -> str:
        return self.message


# CLI/GUI exceptions #
class PJCLIInvalidInputError(Exception):
    par_name: str
    message: str

    def __init__(self, par_name: str, message: str) -> None:
        self.par_name = par_name
        self.message = "ERROR: The argument to " + par_name \
            + " was invalid. " + message

        super().__init__(self.message)

    def __str__(self) -> str:
        return self.message


class PJCLISampleOutOfRangeError(Exception):
    def __init__(self, range_str: str) -> None:
        self.message = "ERROR: The range passed to -f, " \
            + range_str + ", was invalid."

        super().__init__(self.message)

    def __str__(self) -> str:
        return self.message


class PJIOFileDoesNotExistError(Exception):
    def __init__(self, fn_name: str, file_path: str) -> None:
        self.message = "ERROR: When calling " + fn_name + ", could not " \
            + "find file/directory " + file_path + "."

        super().__init__(self.message)

    def __str__(self) -> str:
        return self.message


