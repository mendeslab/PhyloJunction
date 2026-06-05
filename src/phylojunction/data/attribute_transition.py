import typing as ty

__author__ = "Fabio K. Mendes"
__email__ = "fmendes@lsu.edu"


class AttributeTransition():

    attr_label: str
    from_state: int
    to_state: int  # if happening at speciation, child 1
    to_state2: int  # if happening at speciation, child 2
    at_speciation: ty.Optional[bool]

    # (i) if object stored in an AnnotatedTree's at_dict, this is
    # the node subtending a branch (note that this transition may be
    # one that happen right at speciation, but it is nonetheless
    # 'framed' as an anagenetic transition happening along a branch
    # for the purposes of painting branches when tree plotting)
    #
    # (ii) if object is stored in an AnnotatedTree's clado_at_dict,
    # this node is the speciating node; the clado_at_dict is not used
    # for plotting, but instead assists when tabulating and parsing
    # stochastic maps
    #
    # (parent can be dummy node)
    subtending_or_speciating_node_label: str

    # attr transition (forward) occurrence time (i.e., origin or root = 0.0)
    # but note that this will be of either the complete OR the
    # reconstructed process!!!
    global_time: float

    # attribute transition occurrence time (i.e., present time = 0.0)
    age: float

    str_representation: str

    def __init__(self,
                 attr_label: str,
                 subtending_node_label: str,
                 global_time: float,
                 from_state: int,
                 to_state: int,
                 age: ty.Optional[float] = None,
                 to_state2: ty.Optional[int] = None,
                 at_speciation: bool = False) -> None:

        self.attr_label = attr_label
        self.subtending_or_speciating_node_label = subtending_node_label
        self.global_time = global_time
        self.age = age
        self.from_state = from_state
        self.to_state = to_state
        self.to_state2 = to_state2
        self.at_speciation = at_speciation

    def update_daughter_members(self,
                                daughter_node_label: str,
                                daughter_node_time: float) -> None:

        if not self.daughter_node_label:
            # TODO: do stuff
            pass

        # we only need to set daughter members once
        else:
            pass

    def prep_str_representation(self) -> None:
        self.str_representation = \
            "Attribute (\'" + self.attr_label + "\') transition:" + \
            "\n    Time: " + str(self.global_time) + \
            "\n    Age: " + str(self.age) + \
            "\n    Subtending node: " + self.subtending_or_speciating_node_label + \
            "\n    Departing state: " + str(self.from_state) + \
            "\n    Arriving state: " + str(self.to_state)

        # if second arriving state provided, this is a cladogenetic
        # attribute transition, we must adjust the string representation
        if self.to_state2 or self.at_speciation:
            self.str_representation = \
                self.str_representation.replace("Subtending", "Speciating")
            self.str_representation = \
                self.str_representation.replace("Departing", "Parent")
            self.str_representation = \
                self.str_representation.replace("Arriving", "Left child")
            self.str_representation += "\n    Right child state: " + \
                                       str(self.to_state2)

    def __str__(self) -> str:
        self.prep_str_representation()

        return self.str_representation
