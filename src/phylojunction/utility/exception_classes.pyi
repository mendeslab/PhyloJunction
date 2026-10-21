import typing as ty
from _typeshed import Incomplete

class ScriptEvaluationError(Exception):
    source: str
    source_name: str
    line_number: int
    variable: str | None
    execution: int | None
    def __init__(self, message, *, source='', source_name='<script>', line_number=1,
                 variable=None, execution=None): ...
    def __str__(self): ...

class ScriptSyntaxError(Exception):
    cmd_line: str
    message: str
    def __init__(self, cmd_line: str, message: str) -> None: ...

class InexistentVariableError(Exception):
    rv_name: str
    message: str
    def __init__(self, rv_name: str) -> None: ...

class VariableAssignmentError(Exception):
    rv_name: str
    message: str
    def __init__(self, rv_name: str, message: str = ...) -> None: ...

class VariableMisspec(Exception):
    rv_name: str
    message: str
    def __init__(self, rv_name: str) -> None: ...

class FunctionArgsMismatchError(Exception):
    fn_name: str
    message: str
    def __init__(self, fn_name: str, message: str) -> None: ...

class NotBetweenZeroAndOneError(Exception):
    par_name: str
    message: str
    def __init__(self, par_name: str, message: str="") -> None: ...
    def __str__(self) -> str: ...

class IncorrectDimensionError(Exception):
    container_name: str
    obs_len: int
    exp_len: int
    message: str
    obj_name: Incomplete
    def __init__(self, container_name: str, obs_len: int, exp_len: int = ...) -> None: ...

class DimensionalityError(Exception):
    dn_name: str
    par_name: str
    message: str
    def __init__(self, dn_name: str, par_name: ty.Optional[str]="") -> None: ...
    def __str__(self) -> str: ...

class NodeInferenceDimensionalityError(Exception):
    node_name: str
    message: str
    def __init__(self, node_name: str) -> None: ...

class NoPlatingAllowedError(Exception):
    det_name: str
    message: str
    node_dag_name: str
    def __init__(self, det_name: str, problematic_node_dag_name: str, message: str = ...) -> None: ...

class ObjInitRequireSameParameterTypeError(Exception):
    message: str
    obj_name: str
    n_diff_par: int
    def __init__(self, obj_name: str, n_diff_par: int) -> None: ...

class ReplicateNumberError(Exception):
    node_name: str
    message: str
    def __init__(self, node_name, message: str = ...) -> None: ...

class DimensionalityWarning(Exception):
    rv_name: str
    dn_name: str
    message: str
    def __init__(self, rv_name: str, dn_name: str, message: str = ...) -> None: ...

class StateDependentParameterMisspec(Exception):
    message: str
    def __init__(self, message: str = ...) -> None: ...

class SSEStashMisspec(Exception):
    message: str
    def __init__(self, message: str) -> None: ...

class InvalidMCMCChainLength(Exception):
    message: str
    def __init__(self, message: str) -> None: ...

# Initialization errors
class ObjInitInvalidArgError(Exception):
    dn_name: str
    message: str
    def __init__(self, dn_name: str, par_name: str, message: str = "") -> None: ...
    def __str__(self) -> str: ...

class ObjInitMissingParameterError(Exception):
    dn_name: str
    message: str
    def __init__(self, dn_name: str, par_name: str, message: str = "") -> None: ...
    def __str__(self) -> str: ...

class ObjInitIncorrectDimensionError(Exception):
    obj_name: str
    message: str
    obs_len: int
    exp_len: int
    at_least: bool
    def __init__(self, obj_name: str, container_name: str, obs_len: int, exp_len: int = 0, at_least: bool = False) -> None: ...
    def __str__(self) -> str: ...

class ObjInitMissingStateDependentParameterError(Exception):
    message: str
    def __init__(self, epoch_missing_param: int, symmetric_diff_set: ty.Set[ty.Any]) -> None: ...

class ObjInitRepeatedStateDependentParameterError(Exception):
    message: str
    def __init__(self, epoch_w_repeated_param: int, repeated_state: ty.Union[int, ty.Tuple[int]]) -> None: ...
    def __str__(self) -> str: ...

class ObjInitRequireNonZeroStateDependentParameterError(Exception):
    message: str
    obj_name: str
    def __init__(self, obj_name: str, dimension_idx: int) -> None: ...
    def __str__(self) -> str: ...

# Syntax errors: Sampling distribution, det. functions exceptions #
class ParseFunctionArgError(Exception):
    par_name: str
    message: str
    def __init__(self, par_name: str, message: str) -> None: ...

class ParseNotAParameterError(Exception):
    par_name: str
    message: str
    def __init__(self, par_name: str) -> None: ...

class ParseRequireSingleValueError(Exception):
    dn_name: str
    message: str
    arg: Incomplete
    def __init__(self, dn_name: str, arg: str) -> None: ...

class ParseRequireIntegerError(Exception):
    obj_name: str
    message: str
    def __init__(self, obj_name: str, arg: str) -> None: ...

class ParseRequirePositiveIntegerError(Exception):
    obj_name: str
    message: str
    def __init__(self, obj_name: str, arg: str) -> None: ...

class ParseRequireNumericError(Exception):
    dn_name: str
    message: str
    arg: Incomplete
    def __init__(self, dn_name: str, arg: str) -> None: ...

class ParseMissingParameterError(Exception):
    par_name: str
    message: str
    def __init__(self, par_name: str) -> None: ...

class ParseMissingArgumentError(Exception):
    par_name: str
    message: str
    def __init__(self, par_name: str) -> None: ...

class ParseInvalidArgumentError(Exception):
    par_name: str
    message: str
    def __init__(self, par_name: str, invalid_arg: str, message: str) -> None: ...
    def __str__(self) -> str: ...

class ParseMissingSpecificationError(Exception):
    message: str
    def __init__(self, obj2spec_name: str) -> None: ...
    def __str__(self) -> str: ...

class ParseMutuallyExclusiveParametersError(Exception):
    message: str
    par_name: str
    mutually_exclusive_par_name: str
    def __init__(self, par_name: str, mutually_exclusive_par_name: str, message: str) -> None: ...
    def __str__(self) -> str: ...

class ParsePathDoesNotExistError(Exception):
    message: str
    par_name: str
    path_str: str
    def __init__(self, par_name: str, path_str: str, message: ty.Optional[str] = "") -> None: ...
    def __str__(self) -> str: ...

class ParseInvalidNewickStringError(Exception):
    message: str
    par_name: str
    def __init__(self, par_name: str, message: str = "") -> None: ...
    def __str__(self) -> str: ...

class ParseCtFnInitFailError(Exception):
    ct_fn_name: str
    def __init__(self, ct_fn_name: str, message: str) -> None: ...
    def __str__(self) -> str: ...

class ParseDnInitFailError(Exception):
    dn_name: str
    def __init__(self, dn_name: str, message: str) -> None: ...
    def __str__(self) -> str: ...

class ParseDetFnInitFailError(Exception):
    det_fn_name: str
    def __init__(self, det_fn_name: str, message: str) -> None: ...
    def __str__(self) -> str: ...

# PGM exceptions #
class NodeDAGNodeStatCantFloatError(Exception):
    message: str
    def __init__(self, node_name: str) -> None: ...
    def __str__(self) -> str: ...

class DAGCannotInitialize(Exception):
    message: str
    def __init__(self, message: str) -> None: ...
    def __str__(self) -> str: ...

class DAGCannotAddNodeError(Exception):
    message: str
    def __init__(self, node_name: str, message: str) -> None: ...
    def __str__(self) -> str: ...

# Tree exceptions #
class AnnotatedTreeMisspecError(Exception):
    message: str
    def __init__(self, message: str) -> None: ...
    def __str__(self) -> str: ...

class AnnotatedTreeIncorrectAnnotationError(Exception):
    message: str
    def __init__(self, message: str) -> None: ...
    def __str__(self) -> str: ...

class AnnotatedTreeMissingNodeAnnotationError(Exception):
    message: str
    def __init__(self, nd_name: str, attr_name: str) -> None: ...
    def __str__(self) -> str: ...

class AnnotatedTreeNodeMissingAttrError(Exception):
    message: str
    def __init__(self, nd_name: str, attr_name: str, message: str) -> None: ...
    def __str__(self) -> str: ...

class AnnotatedTreeMissingNodeName(Exception):
    message: str
    def __init__(self) -> None: ...
    def __str__(self) -> str: ...
    
# Generation exceptions #
class GenerateFailError(Exception):
    message: str
    def __init__(self, dn_name: str, message: str) -> None: ...
    def __str__(self) -> str: ...

class MissingColumnName(Exception):
    message: str
    def __init__(self, col_name: str, message: str) -> None: ...
    def __str__(self) -> str: ...

class RunTimeLimit(Exception):
    def __init__(self, runtime_limit: float) -> None: ...
    def __str__(self) -> str: ...


class MaxNFailedAttemptsLimit(Exception):
    def __init__(self, max_n_attempts: int) -> None: ...
    def __str__(self) -> str: ...


# CLI/GUI exceptions #
class PJCLIInvalidInputError(Exception):
    par_name: str
    message: str
    def __init__(self, par_name: str, message: str) -> None: ...
    def __str__(self) -> str: ...

class PJCLISampleOutOfRangeError(Exception):
    def __init__(self, range_str: str) -> None: ...
    def __str__(self) -> str: ...

class PJIOFileDoesNotExistError(Exception):
    def __init__(self, fn_name: str, file_path: str) -> None: ...
    def __str__(self) -> str: ...