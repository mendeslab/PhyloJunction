from __future__ import annotations

import copy
import typing as ty
from abc import abstractmethod

# pj imports
import phylojunction.pgm.pgm as pgm
import phylojunction.data.tree as pjtr

__author__ = "Fabio K. Mendes"
__email__ = "fmendes@lsu.edu"


class ClockModel(pgm.DeterministicFn):
    """Abstract base class for molecular clock models.

    A clock model is a deterministic function (PhyloJunction's ':='
    operator) that takes a chronogram DAG node -- one or more
    phylogenetic trees whose branch lengths are measured in time --
    and, for each tree, returns a new tree whose reconstructed
    topology has been annotated with molecular-rate-related
    attributes.

    Every concrete clock model must, at a minimum, populate a node
    attribute, named after 'self.attr_name' (float, in expected
    number of character changes per unit of branch length), on every
    branch of non-zero length in the reconstructed tree. The name of
    this attribute is chosen by the user (e.g., 'mutation_rate' for a
    molecular clock, 'morph_rate' for a morphological one) rather
    than hard-coded, so that PhyloJunction's tree-plotting machinery
    (which can color branches by any such attribute) is not tied to
    any one clock model's naming convention. Concrete subclasses are
    free to also populate additional, model-specific attributes under
    their own, fixed names (e.g., a punctuated, 'Gamma-spike' clock
    model also populates 'branch_raw_rate' and 'spike_effect').

    A clock model does not originate its own sample size the way a
    sampling distribution (i.e., a 'pgm.DistrForSampling' instance)
    does -- it inherits 'n_samples' and 'n_repl' from the tree DAG
    node it is applied to ('tree_dag_node' below).

    Parameters:
        attr_name (str): See '__init__'.
        tree_dag_node (pgm.NodeDAG): See '__init__'.
        n_samples (int): Number of tree samples this clock model will
            be applied to. Mirrors 'tree_dag_node.sample_size'.
        n_repl (int): Number of (plated) replicates per tree sample
            this clock model will be applied to. Mirrors
            'tree_dag_node.repl_size'.
    """

    attr_name: str
    tree_dag_node: pgm.NodeDAG
    n_samples: int
    n_repl: int

    @property
    @abstractmethod
    def CLOCK_MODEL_NAME(self) -> str:
        raise NotImplementedError

    @property
    def DET_FN_NAME(self) -> str:
        """See 'pgm.DeterministicFn.DET_FN_NAME'.

        Every clock model is itself a deterministic function, so this
        concretely satisfies 'pgm.DeterministicFn's contract on
        behalf of every clock model at once, by forwarding to the
        more specific 'CLOCK_MODEL_NAME' each concrete subclass must
        provide.
        """

        return self.CLOCK_MODEL_NAME

    def __init__(self, *, tree_dag_node: pgm.NodeDAG, attr_name: str) -> None:
        """Validate and store state shared by all clock models.

        Concrete subclasses must call this (via 'super().__init__()')
        before storing their own, model-specific hyperparameters.

        Args:
            tree_dag_node (pgm.NodeDAG): DAG node holding the
                chronogram(s) (as a list of 'AnnotatedTree'
                instances, under '.value') this clock model will be
                applied to.
            attr_name (str): Name of the node attribute this clock
                model will populate on every non-zero-length branch
                of a tree's reconstructed topology (e.g.,
                'mutation_rate', 'morph_rate'). Required -- there is
                no default -- so that scripts are always explicit
                about which attribute a clock model is writing to.
                Must be a valid Python identifier, since it is both
                set (via 'setattr') here, and later read (via
                'getattr') by, e.g., PhyloJunction's tree-plotting
                machinery.

        Raises:
            ValueError: Is raised if 'attr_name' is empty, or is not
                a valid Python identifier.
        """

        if not attr_name or not attr_name.isidentifier():
            raise ValueError(
                "'attr_name' must be a non-empty, valid Python "
                "identifier (e.g., 'mutation_rate', 'morph_rate'). "
                "Got: '" + str(attr_name) + "'.")

        self.attr_name = attr_name
        self.tree_dag_node = tree_dag_node
        self.n_samples = tree_dag_node.sample_size
        self.n_repl = tree_dag_node.repl_size

    @abstractmethod
    def populate_branch_rates(
            self, ann_tr: pjtr.AnnotatedTree, tree_idx: int) -> None:
        """Populate per-branch molecular-rate attributes (side-effect).

        This method mutates 'ann_tr' in place; 'generate()' below is
        responsible for deep-copying each input tree before calling
        this method on the copy, so that 'tree_dag_node's own trees
        (this clock model's DAG parent) are never mutated.

        Args:
            ann_tr (AnnotatedTree): Tree instance whose reconstructed
                tree's branches will be annotated with molecular-rate
                attributes, in place.
            tree_idx (int): 'ann_tr's position (0-based) within
                'self.tree_dag_node.value' -- i.e., which (sample,
                replicate) slot this tree is. A concrete clock model
                that consumes another, equally-shaped DAG node for
                per-tree data (e.g., pre-generated per-branch base
                rates) uses this to index into that node's own value
                list, so the two stay in lockstep without either
                model needing to track its own call count.
        """

        pass

    def generate(self) -> ty.List[pjtr.AnnotatedTree]:
        """Apply this clock model to every tree in 'tree_dag_node'.

        This method is the same for every clock model -- what varies
        between concrete clock models is only 'populate_branch_rates'
        -- so it is implemented here, concretely, rather than being
        abstract. Concrete subclasses should not need to override it.

        Returns:
            (list of AnnotatedTree): One new, deep-copied tree per
                tree held by 'self.tree_dag_node', each annotated by
                'populate_branch_rates'.
        """

        ann_tr_list: ty.List[pjtr.AnnotatedTree] = []

        for tree_idx, ann_tr in enumerate(self.tree_dag_node.value):
            ann_tr_copy = copy.deepcopy(ann_tr)
            self.populate_branch_rates(ann_tr_copy, tree_idx)

            # generic marker, set here (once, for every clock model)
            # rather than by each concrete 'populate_branch_rates' --
            # see 'tree.AnnotatedTree.clock_model_attr_name'
            ann_tr_copy.clock_model_attr_name = self.attr_name

            ann_tr_list.append(ann_tr_copy)

        return ann_tr_list


if __name__ == "__main__":
    pass
