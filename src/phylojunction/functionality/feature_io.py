import sys
import os
import enum
import math
import numpy as np
import typing as ty

# pj imports
import phylojunction.utility.helper_functions as pjh
import phylojunction.functionality.evol_event as pjev

__author__ = "Fabio K. Mendes"
__email__ = "fmendes@lsu.edu"


class GeoGraph():
    """
    Parameters:
        n_nodes (int):
        node_set (set):
        edge_set (set):
    """

    _n_nodes: int
    _node_set: ty.Set[int]
    _edge_set: ty.Set[ty.Tuple[int]]
    _edge_dict: ty.Dict[int, ty.Set[int]]
    comm_class_set_list: ty.List[ty.Set[int]]

    # one comm class per node idx (starting at 0)
    _node_comm_class_dict: ty.Dict[int, int]
    _n_comm_classes: int

    def __init__(self, n_nodes: int) -> None:
        self._n_nodes = n_nodes
        self._node_set = set(i for i in range(n_nodes))
        self._edge_set = set()
        self._edge_dict = dict()
        self.comm_class_set_list = list()
        self._node_comm_class_dict = dict()
        self._n_comm_classes = 0

    def add_node(self, node_idx: int) -> None:
        self._n_nodes += 1
        self._node_set.add(node_idx)

    def add_edge(self,
                 node1_idx: int,
                 node2_idx: int,
                 is_directed: bool = False) -> None:

        if not ((node1_idx in self._node_set) \
                and (node2_idx in self._node_set)):
            # TODO: later add exception
            exit("Either or both nodes could not be found in the connectivity graph. Exiting....")

        edge_tup = (node1_idx, node2_idx)
        self._edge_set.add(edge_tup)
        if node1_idx not in self._edge_dict:
            self._edge_dict[node1_idx] = set([node2_idx])

        else:
            self._edge_dict[node1_idx].add(node2_idx)

        if not is_directed:
            edge_tup2 = (node2_idx, node1_idx)
            self._edge_set.add(edge_tup2)

            if node2_idx not in self._edge_dict:
                self._edge_dict[node2_idx] = set([node1_idx])

            else:
                self._edge_dict[node2_idx].add(node1_idx)

    @property
    def edge_set(self) -> ty.Set[ty.Tuple[int]]:
        return self._edge_set

    @property
    def edge_dict(self) -> ty.Dict[int, ty.Set[int]]:
        return self._edge_dict

    @property
    def n_comm_classes(self) -> int:
        return self._n_comm_classes

    def populate_comm_class_members(self) -> None:
        """Populate communicating-class-related members.

        Side-effect is to populate:
            (i)   self.comm_class_set_list
            (ii)  self._node_comm_class_dict
            (iii) self._n_comm_classes

        Communicating classes are mutually exclusive sets of nodes
        that are connected to each other directly or indirectly.
        A connection is established by a single directed edge (even
        if it is missing in the other direction).
        """

        # Treat each directed edge as an undirected connection for component membership.
        # Mark nodes when queued: each is visited once, and traversal exhausts exactly
        # the nodes reachable from its seed. Sorted seeds order components by minimum node;
        # neighbor order cannot change membership. Isolated nodes form singleton components.
        neighbors = {node: set() for node in self._node_set}
        for source, target in self._edge_set:
            neighbors[source].add(target)
            neighbors[target].add(source)

        self.comm_class_set_list = []
        self._node_comm_class_dict = {}
        visited = set()
        for seed in sorted(self._node_set):
            if seed in visited:
                continue
            component = set()
            pending = [seed]
            visited.add(seed)
            while pending:
                node = pending.pop()
                component.add(node)
                for neighbor in neighbors[node]:
                    if neighbor not in visited:
                        visited.add(neighbor)
                        pending.append(neighbor)
            component_idx = len(self.comm_class_set_list)
            self.comm_class_set_list.append(component)
            for node in component:
                self._node_comm_class_dict[node] = component_idx
        self._n_comm_classes = len(self.comm_class_set_list)


    def over_barrier(self,
                     node1_idx: int,
                     node2_idx: int,
                     directed: bool = False) -> bool:
        """Return boolean for two regions depending on barrier existence

        If two nodes (i.e., each node being a region) have an edge
        between them, it means they are are connected directly, and
        thus there is no barrier between them. If they are not
        directly connected, there is a barrier between them.

        When edge direction matters, only edges from node1 to node 2
        are considered. Otherwise the edge on the opposite direction
        is additionally considered.

        Args:
            node1_idx (int): Index of node 1 (region 1).
            node2_idx (int): Index of node 2 (region 2).
            directed (bool): If direction of connectivity edge
                matters. Defaults to 'False'.

        Returns:
            (bool): True if two regions have a barrier between them.
        """

        if node1_idx == node2_idx:
            exit("Can only compare different nodes. Exiting...")

        barrier_exists = node2_idx not in self._edge_dict[node1_idx]

        if directed:
            return barrier_exists

        else:
            barrier_exists = barrier_exists or node1_idx not in self._edge_dict[node2_idx]
            return barrier_exists

    def are_connected(self,
                      node1_idx: int,
                      node2_idx: int,
                      range_bit_pattern: str = "") -> bool:
        """Return boolean for two nodes in the same comm class.

        If two nodes are in the same communicating class AND
        the necessary (in-between) regions in the range are occupied,
        it means they are connected in the connectivity graph.

        Args:
            node1_idx (int): Index of node 1 (region 1).
            node2_idx (int): Index of node 2 (region 2).
            range_bit_pattern (str): Bit pattern of range

        Returns:
            (bool): If two nodes are connected
        """

        if node1_idx == node2_idx:
            exit("Can only compare different nodes. Exiting...")

        # debugging
        # print('self._node_comm_class_dict', self._node_comm_class_dict)
        # print('node1_idx', node1_idx, 'node2_idx', node2_idx)

        cc1 = self._node_comm_class_dict[node1_idx]
        cc2 = self._node_comm_class_dict[node2_idx]

        if cc1 == cc2:
            return True

        return False

    def __str__(self) -> str:
        str_representation = ""
        for e_tup in self._edge_set:
            e_str = " -> ".join(str(i) for i in e_tup)
            str_representation += e_str + "\n"

        return str_representation


class BarrierAppearance(pjev.EvolRelevantEvent):
    """Barrier appearance event.

    Parameters:
        from_node_idx (int): Integer representation of source node
            in barrier.
        to_node_idx (int): Integer representation of target node
            in barrier.
    """

    from_node_idx: int
    to_node_idx: int
    destabilized_range: ty.Optional[str]
    range_bit_patt: ty.Optional[str]

    def __init__(self,
                 n_chars: int,
                 age: ty.Optional[float],
                 from_node_idx: int,
                 to_node_idx: int,
                 time: ty.Optional[int] = None) -> None:

        # initializing EvolRelevantEvent members
        super().__init__(n_chars, age, time=time)

        self.from_node_idx = from_node_idx
        self.to_node_idx = to_node_idx
        self.destabilized_range = None
        self.range_bit_patt = None

    def short_str(self) -> str:
        short_str = "b+(" + str(self.age) + ")" \
                    + "_reg(" + str(self.from_node_idx) + "|" \
                    + str(self.to_node_idx) + ")"

        if self.destabilized_range is not None and \
                self.range_bit_patt is not None:
            short_str += "_" + self.destabilized_range \
                + "(" + self.range_bit_patt + ")"

        return short_str

    def __str__(self) -> str:
        str_representation = \
            "Barrier (at age = " + str(self.age) \
            + ") between region " + str(self.from_node_idx) \
            + " and " + str(self.to_node_idx)

        if self.destabilized_range is not None and \
                self.range_bit_patt is not None:
            str1 = "\n  Destabilized range (" + self.range_bit_patt + "): "
            str2 = "True" if self.destabilized_range == "destab" \
                else "N/A"
            suffix = str1 + str2
            str_representation += suffix

        return str_representation

    def __lt__(self, other) -> bool:
        return super().__lt__(other)


class BarrierDisappearance(pjev.EvolRelevantEvent):
    """Barrier disappearance event.

    Parameters:
        from_node_idx (int): Integer representation of source node
            in barrier.
        to_node_idx (int): Integer representation of target node
            in barrier.
    """

    from_node_idx: int
    to_node_idx: int
    restabilized_range: ty.Optional[str]
    range_bit_patt: ty.Optional[str]  # range at time of event

    def __init__(self,
                 n_chars: int,
                 age: ty.Optional[float],
                 from_node_idx: int,
                 to_node_idx: int,
                 time: ty.Optional[int] = None) -> None:

        # initializing EvolRelevantEvent members
        super().__init__(n_chars, age, time=time)

        self.from_node_idx = from_node_idx
        self.to_node_idx = to_node_idx
        self.restabilized_range = None
        self.range_bit_patt = None

    def short_str(self) -> str:
        short_str = "b-(" + str(self.age) + ")" \
            + "_reg(" + str(self.from_node_idx) + "|" \
            + str(self.to_node_idx) + ")"

        if self.restabilized_range is not None and \
                self.range_bit_patt is not None:
            short_str += "_" + self.restabilized_range \
                + "(" + self.range_bit_patt + ")"

        return short_str

    def __str__(self) -> str:
        str_representation = \
            "Connectivity restablished (at age = " \
            + str(self.age) + ") between region " \
            + str(self.from_node_idx) + " and " \
            + str(self.to_node_idx)

        if self.restabilized_range is not None and \
                self.range_bit_patt is not None:
            str1 = "\n  Restabilized range (" + self.range_bit_patt + "): "
            str2 = "True" if self.restabilized_range == "stab" \
                else "N/A"
            suffix = str1 + str2
            str_representation += suffix

        return str_representation

    def __lt__(self, other) -> bool:
        return super().__lt__(other)


class GeoFeatureRelationship(enum.Enum):
    WITHIN = "within"
    BETWEEN = "between"


class GeoFeatureType(enum.Enum):
    CATEGORICAL = "categorical"
    QUANTITATIVE = "quantitative"


class GeoFeature():
    """Geographic feature class.

    Parameters:
        time_idx (int): Index of epoch when feature is measured.
        feat_idx (int): Index of feature in its container.
        n_regions (int): Number of atomic (i.e., unsplittable) regions
            in the system.
        feat_rel (GeoFeatureRelationship): Enum object that
            characterizes how the feature relates one or more regions
            ("within" or "between)".
        feat_type (GeoFeatureType): Enum object carrying the type of
            the feature ("categorical" or "quantitative").
        n_dim (int): Number of dimensions in the container containing
            geographic feature values.
        value (NumPy.array): Value(s) measured for geographic feature.
        name (str): Geographic feature name.
        str_representation (str): String representation of geographic
            feature for printing.
    """

    time_idx: int  # epoch index
    feat_idx: int  # should be unique to this feature
    n_regions: int  # automatically computed
    feat_rel: GeoFeatureRelationship
    feat_type: GeoFeatureType
    value: np.array
    name: str
    str_representation: str
    n_dim: int

    def __init__(self,
                 time_idx: int,
                 feat_idx: int,
                 feat_rel: GeoFeatureRelationship,
                 feat_type: GeoFeatureType,
                 feat_vals: ty.Optional[np.array] = None,
                 feat_name: ty.Optional[str] = ""):

        self.time_idx = time_idx
        self.feat_idx = feat_idx
        self.feat_rel = feat_rel
        self.feat_type = feat_type
        self.value = feat_vals
        self.name = feat_name
        self.n_regions = None

        if feat_type == GeoFeatureRelationship.BETWEEN:
            self.n_dim = 2
        
        if self.value is not None:
            # initialized here or when self.value is 
            # set by setter
            self._initialize_n_regions()

        self._initialize_str_representation()

    # init methods
    def _initialize_str_representation(self) -> None:
        self.str_representation = \
            "Feature | " + self.feat_rel.value + "-" \
            + self.feat_type.value + " " + str(self.feat_idx) + " | " \
            + "Epoch " + str(self.time_idx) 
        
        if self.name:
            str_replacement = "Feature (" + self.name + ")"
            self.str_representation = \
                self.str_representation.replace("Feature", str_replacement)

        if self.n_regions:
            self.str_representation += " | " + str(self.n_regions) \
                + " regions"

        if self.value is not None:
            if self.feat_rel == GeoFeatureRelationship.WITHIN:
                self.str_representation += "\n" + np.array2string(self.value)

            elif self.feat_rel == GeoFeatureRelationship.BETWEEN:
                self.str_representation += "\n" + np.array2string(
                    self.value.reshape(self.n_regions,
                                       self.n_regions))

    # setters
    def _set_attr_values(self,
                         attr_name: str,
                         val: ty.Tuple[int,
                                       GeoFeatureRelationship,
                                       GeoFeatureType,
                                       str,
                                       np.array]) -> None:
        self.__setattr__(attr_name, val)

        # automatically setting self.n_regions
        if attr_name == "value":
            self._initialize_n_regions()

        # force update
        self._initialize_str_representation()

    # getters
    def get_val_from_to(self, from_idx: int, to_idx: int) -> float:
        if not self.n_regions:
            exit("ERROR: Feature " + self.name + " is not aware of " \
                 + "how many regions the feature is scored for. Exiting.")
        
        return self.value[from_idx * self.n_regions + to_idx]

    # health
    def _initialize_n_regions(self) -> None:
        if self.feat_rel == GeoFeatureRelationship.WITHIN:
            self.n_regions = len(self.value)

        elif self.feat_rel == GeoFeatureRelationship.BETWEEN:
            self.n_regions = int(math.sqrt(len(self.value)))
    
    # other
    def __str__(self) -> None:
        return self.str_representation


class GeoFeatureCollection():
    """Collection of geographic features.

    This class organizes

    Parameters:
        geofeat_list (GeoFeature): List of GeoFeature objects.
        epoch_age_end_list_young2old (list[float]): List of age ends
            (float), one per epoch, going from the youngest epoch
            to the oldest.
        epoch_mid_age_list_young2old (list[float]): List of age mid points
            (float), one per epoch, going from the youngest epoch
            to the oldest.
        epoch_age_start_list_young2old (list[float]): List of age starts
        n_epochs (int): Number of epochs
        is_timehet (bool): Flag specifying if the geographic data
            is time-heterogeneous.
        feat_name_epochs_dict (dict): Nested dictionaries. The outer
            dictionary has feature names (str) as keys (i.e., features
            must have different names), and the inner dictionaries as
            values. The inner dictionary has epoch indices as keys
            (int), and GeoFeature objects as values.
        feat_type_rel_featid_epochs_dict (dict): Quadruple nested
            dictionaries. The outer keys are GeoFeatureType. The next
            nested keys are GeoFeatureRelationship. The third nested
            keys are the index (int) of the feature. The fourth and
            last nested keys are epoch indices (int), with values being
            GeoFeatures.
        region_name_idx_dict (dict): Dictionary for converting between
            region names and their indices. Names (str) for keys,
            indices (int) for values.
        region_idx_name_dict (dict): Dictionary for converting between
            region indices and their names. Indices (int) for keys,
            names (str) for values.
    """

    geofeat_list: ty.List[GeoFeature]
    epoch_age_end_list_young2old: ty.List[float]
    epoch_mid_age_list_young2old: ty.List[float]
    epoch_age_start_list_young2old: ty.List[float]
    n_epochs: int
    is_timehet: bool

    # ways in which to store (for later grabbing) features values
    #
    # NOTE: time_idx is user provided, class is unaware of what is
    # old and young
    #
    # (1) by name (assumes name cannot be repeated across feature
    # types and/or relationships); final dict is k: time_idx, v: feat
    #
    # TODO: add check here!!!!
    feat_name_epochs_dict: ty.Dict[str, ty.Dict[int, GeoFeature]]
    # (2) by feature type, relationship and feature number (index);
    # final dict is k: time_idx, v: feat
    feat_type_rel_featid_epochs_dict: \
        ty.Dict[GeoFeatureType, 
                ty.Dict[GeoFeatureRelationship,
                        ty.Dict[int,
                                ty.Dict[int, GeoFeature]
                                ]]]

    region_name_idx_dict: ty.Dict[str, int]
    region_idx_name_dict: ty.Dict[int, str]

    def __init__(self,
                 feat_summary_fp: str,
                 age_summary_fp: str = "") -> None:

        # starting class members off
        self.geofeat_list = list()
        self.epoch_age_end_list_young2old = [0.0]  # present epoch ends at 0.0
        self.epoch_mid_age_list_young2old = [math.inf]  # oldest epoch midpoint
        self.epoch_age_start_list_young2old = list()
        self.is_timehet = False
        self.feat_name_epochs_dict = pjh.autovivify(2)
        self.feat_type_rel_featid_epochs_dict = pjh.autovivify(4)
        self.region_name_idx_dict = dict()
        self.region_idx_name_dict = dict()

        self._check_filepaths(feat_summary_fp)


        # if only one epoch, age_summary not necessary
        if age_summary_fp is not None and age_summary_fp != "":
            self._check_filepaths(age_summary_fp)

            # initializes
            #     self.epoch_age_end_list_young2old
            #     self.epoch_mid_age_list_young2old
            #     self.epoch_age_end_list_old2young
            #     self.epoch_mid_age_list_old2young
            self._read_age_summary(age_summary_fp)
            self.n_epochs = \
                len(self.epoch_age_end_list_young2old)

        # initializes:
        #     self.is_timehet
        #     self.region_name_idx_dict
        #     self.region_idx_name_dict
        self._read_feat_summary_init_feats(feat_summary_fp)

        # initializes self.feat_name_epochs_dict
        self._init_feat_name_epochs_dict()

        # initializes self.feat_type_rel_featid_epochs_dict
        self._init_feat_type_rel_featid_epochs_dict()

        # debugging
        # self._debug_print_dicts()

    # init methods
    def _read_feat_summary_init_feats(
            self,
            feat_summary_fp):
        """Populate feature list with initial info
        
        Open, read and parse feature summary .csv file, and populate
        class member list of features with initial information
        """

        with open(feat_summary_fp, "r") as infile:
            header = infile.readline()
            if header != ("time_index,feature_index,feature_relationship,"
                          "feature_type,feature_path,name\n"):
                exit(("ERROR: Expecting header:\ntime_index,feature_index,"
                      "feature_relationship,feature_type,feature_path,name\n"))
            
            # linecount = 0
            # feat_unique_names = set([])
            for line in infile:
                line = line.rstrip()
                time_idx, feat_idx, feat_rel, feat_type, feat_fp, feat_name \
                    = line.split(",")

                if int(time_idx) > 1:
                    self.is_timehet = True

                geofeat = GeoFeature(
                    int(time_idx),
                    int(feat_idx),
                    (GeoFeatureRelationship.WITHIN \
                    if feat_rel == GeoFeatureRelationship.WITHIN.value \
                    else GeoFeatureRelationship.BETWEEN),
                    (GeoFeatureType.CATEGORICAL \
                    if feat_type == GeoFeatureType.CATEGORICAL.value \
                    else GeoFeatureType.QUANTITATIVE),
                    feat_name = feat_name)
                
                self._check_filepaths(feat_fp)
                
                with open(feat_fp, "r") as infile2:
                    regions_str = infile2.readline().rstrip()
                    regions_list = regions_str.split(",")

                    for idx, reg in enumerate(regions_list):
                        self.region_name_idx_dict[reg] = idx
                        self.region_idx_name_dict[idx] = reg

                    vals_list = list()
                    for line in infile2:
                        line = line.rstrip().split(",")
                        vals_list.extend(int(v) \
                                         if feat_type == GeoFeatureType.CATEGORICAL.value \
                                         else float(v) \
                                         for v in line)
                    
                    vals_array = np.array(vals_list)
                
                geofeat._set_attr_values("value", vals_array)
                
                # debugging
                # print(geofeat)

                self.geofeat_list.append(geofeat)

                # linecount += 1
            
            # if len(feat_unique_names) != linecount:
            #     exit("ERROR: Detected " + str(len(feat_unique_names)) \
            #          + " unique feature names, but summary file had " \
            #          + str(linecount) + " lines minus the header. Exiting.")

    def _read_age_summary(
            self,
            age_summary_fp: str) -> None:
        """Populate class members holding information about epochs
        
        Open, read and parse age summary .csv file. then populate class
        members.

        Side-effects are populating:
            (i)   a
            (ii)  a
            (iii) a

        Args:
            age_summary_fp (str): Path to .csv file containing
                information about the age ends of the different epochs
                should geographic data be time-heterogeneous.
        """

        ############################
        # Doing end age dictionary #
        ############################

        with open(age_summary_fp, "r") as infile:
            header = infile.readline()
            
            if header != ("index,age_end\n"):
                exit("ERROR: Expecting header:\nindex,age_end\n")

            content = infile.readlines()
            n_epochs_minus_one = len(content)
            epoch_age_ends_from_file = \
                [-1.0 for i in range(n_epochs_minus_one)]

            for line in content:
                idx, age_end = line.rstrip().split(",")
                # note the offset!
                epoch_age_ends_from_file[int(idx)-1] = float(age_end)

        # actual initialization
        self.epoch_age_end_list_young2old += epoch_age_ends_from_file
        self.epoch_age_end_list_old2young \
            = self.epoch_age_end_list_young2old[::-1]

        # debugging
        # print(self.epoch_age_end_list_young2old)

        ############################
        # Doing mid age dictionary #
        ############################

        n_epochs = len(self.epoch_age_end_list_young2old)
        mid_ages_minus_last_epoch = list()
        for idx in range(len(self.epoch_age_end_list_young2old)):
            if idx < n_epochs-1:
                younger_age = self.epoch_age_end_list_young2old[idx]
                older_age = self.epoch_age_end_list_young2old[idx+1]
                mid_age = (younger_age + older_age) / 2.0

                mid_ages_minus_last_epoch.append(mid_age)

                # debugging
                # print(younger_age, older_age, mid_age)

        # actual initialization
        self.epoch_mid_age_list_young2old \
            = mid_ages_minus_last_epoch \
            + self.epoch_mid_age_list_young2old
        
        self.epoch_mid_age_list_old2young \
            = self.epoch_mid_age_list_young2old[::-1]
        
        # debugging
        # print(self.epoch_mid_age_list_young2old)

        ##############################
        # Doing start age dictionary #
        ##############################

        for age_end in self.epoch_age_end_list_young2old[1:]:
            self.epoch_age_start_list_young2old.append(age_end)

        self.epoch_age_start_list_young2old.append(math.inf)
        self.epoch_age_start_list_old2young = \
            self.epoch_age_start_list_young2old[::-1]


    def _init_feat_name_epochs_dict(self):
        for geofeat in self.geofeat_list:
            self.feat_name_epochs_dict[geofeat.name][geofeat.time_idx] \
                = geofeat

    def _init_feat_type_rel_featid_epochs_dict(self):
        for geofeat in self.geofeat_list:
            fty = geofeat.feat_type
            fr = geofeat.feat_rel
            fid = geofeat.feat_idx
            fti = geofeat.time_idx
            self.feat_type_rel_featid_epochs_dict[fty][fr][fid][fti] \
                = geofeat

    # health check methods
    def _check_filepaths(self,
                         feat_summary_fp) -> None:
        """Check files exist in initialization file paths"""

        if not os.path.isfile(feat_summary_fp):
            exit("Could not find file with feature summary table " + feat_summary_fp + ". Exiting.")

    # getters
    def get_feat_by_name(self,
                         feat_name: str,
                         time_idx: ty.Optional[int] = None) \
            -> ty.Tuple[float, ty.List[float]]:

        feat_dict_all_epochs = self.feat_name_epochs_dict[feat_name]

        if time_idx:
            return feat_dict_all_epochs[time_idx]

        else:
            return [v for k, v in feat_dict_all_epochs.items()]
    
    def get_feat_by_type_rel_id(self,
                                feat_type: GeoFeatureType,
                                feat_rel: GeoFeatureRelationship,
                                feat_id: int,
                                time_idx: ty.Optional[int] = None) \
            -> ty.Tuple[float, ty.List[float]]:
        """

        Args:
            feat_type:
            feat_rel:
            feat_id:
            time_idx: Index of epoch for which we want to get features.
                Note that these indices are provided by the user in
                the feature summary .csv file.

        Returns:
            (GeoFeature): Geographic feature (GeoFeature object) or
                list of geographic features, depending on whether
                time_idx was provided.
        """

        # digging into the 4-nested dictionary
        feat_rel_dict_all_epochs = \
            self.feat_type_rel_featid_epochs_dict[feat_type]
        feat_id_dict_all_epochs = \
            feat_rel_dict_all_epochs[feat_rel]
        feat_dict_all_epochs = \
            feat_id_dict_all_epochs[feat_id]
        
        # print(self.feat_type_rel_featid_epochs_dict.keys())
        # print(feat_rel_dict_all_epochs.keys())
        # print(feat_id_dict_all_epochs.keys())
        # print(feat_dict_all_epochs.keys())
        # print(feat_dict_all_epochs[1])
        
        if time_idx:
            return feat_dict_all_epochs[time_idx]

        else:
            # returning GeoFeature list sorted by the
            # epoch indices provided by the user
            return [feat_dict_all_epochs[k] for k in sorted(feat_dict_all_epochs)]

    # debug
    def _debug_print_dicts(self) -> None:
        for feat_name, feat_dict in self.feat_name_epochs_dict.items():
            for t, feat in feat_dict.items():
                print(feat_name, ":\n", feat, sep="")
        
        for fty, feat_dict1 in self.feat_type_rel_featid_epochs_dict.items():
            for fr, feat_dict2 in feat_dict1.items():
                for fid, feat_dict3 in feat_dict2.items():
                    for fti, feat in feat_dict3.items():
                        print(fty.value + "-" + fr.value + " feature " \
                              + str(fid) + " epoch " + str(fti),
                              feat,
                              sep="\n")


# callable type for type hinting requirement functions
class MyCallableType(ty.Protocol):
    def __call__(self, feat_col: GeoFeatureCollection) -> None:
        ...


class GeoFeatureQuery():
    """

    Parameters:
        n_regions (int): Number of regions in the system.
        n_epochs (int): Number of time slices (i.e., epoch)  for which
            one has geographic features scored.
        feat_col (GeoFeatureCollection):
        geo_cond_bit_dict (dict): Dictionary with keys being the name of
            a paleogeographical condition (e.g., 'land_bridge'), and
            values being 1- or 2-D lists with bit patterns. Each bit
            in the pattern corresponds to a time slice (i.e., epoch),
            in whatever temporal order the user specifies in the feature
            and age summary files (conventionally from young to old).
            If the paleogeographical condition is a function of a
            within-region feature, the value of the dictionary is a
            list; if of a between-region feature, the value will be a
            2-D list.
        conn_graph_list (GeoGraph): List of connectivity graphs, one
            per epoch.
        geo_cond_name (str): Name of geographic condition verified for
            determining the bit value within or between regions (e.g.,
            a barrier).
    """

    n_regions: int
    n_time_slices: int

    feat_coll: GeoFeatureCollection

    # bit patterns are in whatever order the time slices appear
    # in the user-specified feature summary file (according
    # to each time slice index)
    #
    # the value is either a list or a nested list depending on
    # whether the feature was between or within
    geo_cond_bit_dict: ty.Dict[str,
                               ty.Union[ty.List[ty.List[str]],
                               ty.List[str]]]

    # the value is a 2-D or 3-D list of floats, with all ages
    # of geographic condition change as indicated by a bit in
    # the bit pattern flipping, '01' (0 in epoch 1, 1 in epoch 2);
    # if no flips, no changes
    geo_cond_change_times_dict: \
        ty.Dict[str,
                ty.Union[ty.List[ty.List[float]],
                ty.List[ty.List[ty.List[float]]]]]

    # same as above, but with the bit flipping the other way around, '10'
    geo_cond_change_back_times_dict: \
        ty.Dict[str,
                ty.Union[ty.List[ty.List[float]],
                ty.List[ty.List[ty.List[float]]]]]

    geo_cond_name: str

    # DEPRECATED
    # this member below helps us tell if geographic connectivity had
    # previously happened, e.g., maybe a barrier never "appears" during the
    # considered epochs, but that's because it appeared in the unobservable
    # past
    # (the geog. condition bit pattern could be '111111'; there is no '01' in there
    # but that is not because a barrier never appeared between the focal two
    # regions -- instead, it always existed for the analyzed time window)
    #
    # also, if a barrier between a specific pair of regions is a 'component' of a
    # 'full' barrier between complex ranges, we need to know if the 'component'
    # barrier was/wasn't always there, so we can say the 'full' barrier is
    # manifested
    # geo_oldest_cond_bit_dict: ty.Dict[str,
    #                               ty.Union[ty.List[ty.List[str]],
    #                                        ty.List[str]]]

    # connectivity graph list, one graph per epoch
    conn_graph_list: ty.List[GeoGraph]

    def __init__(self, feat_coll) -> None:
        self.feat_coll = feat_coll
        self.n_regions = len(self.feat_coll.region_idx_name_dict)
        self.n_time_slices = self.feat_coll.n_epochs

        self.geo_cond_bit_dict = dict()
        self.geo_oldest_cond_bit_dict = dict()
        self.geo_cond_change_times_dict = dict()
        self.geo_cond_change_back_times_dict = dict()
        self.conn_graph_list = list()

    # class methods take 'cls' as first argument, and can know about class state
    # static methods do not know about state at all

    # requirement function 1
    @staticmethod
    def cw_feature_equals_value(
            feat_col: GeoFeatureCollection,
            val2match: int,
            feat_name: str = "",
            feat_id: int = None) -> ty.List[str]:
        """Determine geographic condition from (feat.) value requirement.

        The feature must be 'categorical-within'.

        Args:
            feat_col (GeoFeatureCollection):
            val2match (int): Categorical (feature) value that, if
                matched, indicate the manifestation of the geographic
                condition.
            feat_name (str): Name of the feature.
            feat_id (int, optional): Feature index. Defaults to 'None'.

        Returns:
            (list): List of strings. These strings are bit patterns
                denoting the manifestation of the geographic condition
                for each region, in the different epochs.
        """
    
        ft_all_epochs_list: ty.List[GeoFeature] = list()
        
        if feat_name:
            # young -> old
            ft_all_epochs_list = feat_col.get_feat_by_name(feat_name)
    
        elif feat_id:
            # young -> old
            ft_all_epochs_list = feat_col.get_feat_by_type_rel_id(
                GeoFeatureType.CATEGORICAL,
                GeoFeatureRelationship.WITHIN,
                feat_id
            )

        else:
            exit(("ERROR: Function \'cw_feature_equals_value\' needs a "
                    "categorical-within feature name or feature id."
                    ". Exiting."))

        n_regions = ft_all_epochs_list[0].n_regions  # just grabbing 1st
        # mid_ages = feat_col.epoch_mid_age_list_young2old
        
        # mid_ages_match_list = list()
        bits_match_list = list()

        for region_idx in range(n_regions):
            # mid_ages_match_this_region = list()
            bits_match_this_region_str = str()
            
            # young to old
            # for epoch_idx, ft_this_epoch in enumerate(ft_all_epochs_list):
            for ft_this_epoch in ft_all_epochs_list:
                if ft_this_epoch.value[region_idx] == val2match:
                    # mid_ages_match_this_region.append(
                    #     (epoch_idx, mid_ages[ft_this_epoch.time_idx-1]))
                    bits_match_this_region_str += "1"

                else:
                    bits_match_this_region_str += "0"

            # mid_ages_match_list.append(mid_ages_match_this_region)
            bits_match_list.append(bits_match_this_region_str)

        # return mid_ages_match_list
        return bits_match_list
    
    # requirement function 2
    @staticmethod
    def qw_feature_threshold(
            feat_col: GeoFeatureCollection,
            threshold: float,
            above: bool,
            feat_name: str = "",
            feat_id: int = None) \
            -> ty.List[ty.List[float]]:
    
        ft_all_epochs_list: ty.List[GeoFeature] = list()

        if feat_name:
            # young -> old
            ft_all_epochs_list = feat_col.get_feat_by_name(feat_name)

        elif feat_id:
            # young -> old
            ft_all_epochs_list = feat_col.get_feat_by_type_rel_id(
                GeoFeatureType.QUANTITATIVE,
                GeoFeatureRelationship.WITHIN,
                feat_id
            )

        else:
            exit(("ERROR: Function \'qw_feature_threshold\' needs a "
                    "quantitative-within feature name or feature id."
                    ". Exiting."))
            
        n_regions = ft_all_epochs_list[0].n_regions  # just grabbing 1st
        # mid_ages = feat_col.epoch_mid_age_list_young2old
        
        # mid_ages_match_list = list()
        bits_match_list = list()

        for region_idx in range(n_regions):
            # mid_ages_match_this_region = list()
            bits_match_this_region_str = str()

            # for epoch_idx, ft_this_epoch in enumerate(ft_all_epochs_list):
            for ft_this_epoch in ft_all_epochs_list:
                if above:
                    if ft_this_epoch.value[region_idx] > threshold:
                        # mid_ages_match_this_region.append(
                        #     (epoch_idx, mid_ages[ft_this_epoch.time_idx-1]))
                        bits_match_this_region_str += "1"

                    else:
                        bits_match_this_region_str += "0"

                else:
                    if ft_this_epoch.value[region_idx] < threshold:
                        # mid_ages_match_this_region.append(
                        #     (epoch_idx, mid_ages[ft_this_epoch.time_idx-1]))
                        bits_match_this_region_str += "1"

                    else:
                        bits_match_this_region_str += "0"

            # mid_ages_match_list.append(mid_ages_match_this_region)
            bits_match_list.append(bits_match_this_region_str)

        return bits_match_list
        # return mid_ages_match_list
        
    # requirement function 3
    @staticmethod
    def cb_feature_equals_value_is_connected(
            feat_col: GeoFeatureCollection,
            val2match: int,
            feat_name: str = "",
            feat_id: ty.Optional[int] = None) \
            -> ty.List[ty.List[str]]:
        """Determine geographic connection from (feat.) value.

        This function builds bit patterns for all pairs of regions.
        The bit pattern will have as many bits as there are time
        slices (epochs), and will be sorted according to the time
        indices provided in the feature summary .csv file specified
        by the user. So if the youngest epoch has time index 1 and the
        oldest has time index 3, a bit pattern '100' indicates that
        for a given pair of regions, at the present moment (youngest
        epoch) there is connectivity.

        Args:
            feat_col (GeoFeatureCollection): A GeoFeatureCollection
                object containing all information available for
                features over time.
            val2match (int): Categorical (feature) value that, if
                matched, indicate the geographic connection between
                two regions.
            feat_name (str): Name of the feature.
            feat_id (int, optional): Feature index. Defaults to 'None'.

        Returns:
            (list): Two-dimensional list of strings. These strings are
                bit patterns denoting the manifestation of a geographic
                condition (e.g., a barrier), for each pair of regions,
                for each epoch.
        """
        
        ft_all_epochs_list: ty.List[GeoFeature] = list()

        if feat_name != "":
            # young -> old
            ft_all_epochs_list = feat_col.get_feat_by_name(feat_name)

        elif feat_id:
            # young -> old
            ft_all_epochs_list = feat_col.get_feat_by_type_rel_id(
                GeoFeatureType.CATEGORICAL,
                GeoFeatureRelationship.BETWEEN,
                feat_id
            )

        else:
            exit(("ERROR: Function \'cb_feature_threshold\' needs a "
                    "categorical-between feature name or feature id."
                    ". Exiting."))

        n_regions = ft_all_epochs_list[0].n_regions  # just grabbing 1st

        # 1D: "from" region
        # 2D: "to" region
        # values: bit pattern where 1 means requirement met, 0 otherwise
        bits_match_2d_list = list()

        # nested iterations to check all pairs of regions
        for region_idx1 in range(n_regions):
            bits_match_from_region = list()
            
            for region_idx2 in range(n_regions):
                bits_match_to_region_str = str()   

                # for epoch_idx, ft_this_epoch in enumerate(ft_all_epochs_list):
                for ft_this_epoch in ft_all_epochs_list:
                    # bit is 1 if requirement is met
                    if ft_this_epoch \
                        .get_val_from_to(region_idx1, region_idx2) \
                            == val2match:
                        bits_match_to_region_str += "1"

                    # bit is 0 if requirement is not met
                    else:
                        bits_match_to_region_str += "0"

                bits_match_from_region.append(bits_match_to_region_str)

            bits_match_2d_list.append(bits_match_from_region)

        return bits_match_2d_list
        
    # requirement function 4
    @classmethod
    def qb_feature_threshold(cls,
                             feat_col: GeoFeatureCollection,
                             threshold: float,
                             above: bool,
                             feat_name: str = "",
                             feat_id: int = None
                             ) \
            -> ty.List[ty.List[str]]:
            # -> ty.List[ty.List[float]]:
        
        ft_all_epochs_list: ty.List[GeoFeature] = list()

        if feat_name:
            # young -> old
            ft_all_epochs_list = feat_col.get_feat_by_name(feat_name)

        elif feat_id:
            # young -> old
            ft_all_epochs_list = feat_col.get_feat_by_type_rel_id(
                GeoFeatureType.QUANTITATIVE,
                GeoFeatureRelationship.BETWEEN,
                feat_id
            )

        else:
            exit(("ERROR: Function \'qw_feature_threshold\' needs a "
                    "quantitative-within feature name or feature id."
                    ". Exiting."))
            
        n_regions = ft_all_epochs_list[0].n_regions  # just grabbing 1st
        # mid_ages = feat_col.epoch_mid_age_list_young2old

        # mid_ages_match_2d_list = list()
        bits_match_2d_list = list()

        for region_idx1 in range(n_regions):
            # mid_ages_match_from_region = list()
            bits_match_from_region = list()
            
            for region_idx2 in range(n_regions):
                # mid_ages_match_to_region = list()
                bits_match_to_region_str = str()

                # for epoch_idx, ft_this_epoch in enumerate(ft_all_epochs_list):
                for ft_this_epoch in ft_all_epochs_list:
                    if above:
                        if ft_this_epoch .get_val_from_to(region_idx1, region_idx2) \
                                > threshold:
                            # second element is a True/False bit (1 = True, 0 = False)
                            # mid_ages_match_to_region.append(
                            #     (epoch_idx, 1, mid_ages[ft_this_epoch.time_idx-1]))
                            bits_match_to_region_str += "1"

                        else:
                            # mid_ages_match_to_region.append(
                            #     (epoch_idx, 0, None))
                            bits_match_to_region_str += "0"

                    else: 
                        if ft_this_epoch.get_val_from_to(region_idx1, region_idx2) \
                                < threshold:
                            # second element is a True/False bit (1 = True, 0 = False)
                            # mid_ages_match_to_region.append(
                            #     (epoch_idx, 1, mid_ages[ft_this_epoch.time_idx-1]))
                            bits_match_to_region_str += "1"

                        else:
                            # mid_ages_match_to_region.append(
                            #     (epoch_idx, 0, None))
                            bits_match_to_region_str += "0"

                # mid_ages_match_from_region.append(mid_ages_match_to_region)
                # reverse so it's old to young
                bits_match_from_region.append(bits_match_to_region_str[::-1])

            # mid_ages_match_2d_list.append(mid_ages_match_from_region)
            bits_match_2d_list.append(bits_match_from_region)

        # return mid_ages_match_2d_list
        return bits_match_2d_list
                          
    def populate_geo_cond_member_dicts(
            self,
            geo_cond_name: str,
            w_bw_bit_patterns: ty.Union[ty.List[ty.List[str]],
                               ty.List[str]],
            is_directed: bool = False,
            user_young2old: bool = True) -> None:
        """Populate geographic connectivity bit dictionaries

        e.g., for 2 regions, 3 epochs {"altitude": ["010", "010"]}
        the 0's represent no geographic connectivity, while 1's
        represent connectivities (depending on the requirement
        function provided by the user).

        Populates self.geo_cond_name.
        
        Args:
            geo_cond_name (str): Name of the geographic
                condition (e.g., \"mountain range\", \"ancient sea\",
                \"riverine barrier\").
            w_bw_bit_patterns (list): Within or between bit patterns
                returned from connectivity function, one bit per epoch.
                If within-region features, this will be a list of
                strings. If between-region features, this will be a 2-D
                list.
            is_directed (bool): Flag specifying if connectivity graph
                is directed or not. Defaults to False.
            user_young2old (bool): Flag specifying if epochs (their
                ages specifically) are passed from young to old.
                Defaults to True.
        """

        self.geo_cond_name = geo_cond_name

        self.geo_cond_bit_dict[geo_cond_name] = w_bw_bit_patterns

        # DEPRECATED
        # self._populate_oldest_geo_cond_bit_dict(geo_cond_name)

        self._populate_geo_cond_change_times_dict(geo_cond_name)

        self._populate_conn_graph_list(geo_cond_name,
                                       is_directed,
                                       user_young2old)

    # internal
    def _populate_geo_cond_change_times_dict(self,
                                             geo_cond_name: str) -> None:
        """Populate class member holding geog. condition change times.

        This method has the side-effect of populating:
            (i)  self.geo_cond_change_times_dict, and
            (ii) self.geo_cond_change_back_times_dict,

        whose keys are names of the geographic condition, and values
        are the epoch times at which geographic conditions either
        change (in (i)) or change back (in (ii)).

        Args:
            geo_cond_name (str): Name of the geographic condition
            (e.g., 'mountain range', 'ancient sea', 'riverine barrier').
        """
        
        age_starts_old2young = self.feat_coll.epoch_age_start_list_old2young

        # initializing dictionary values
        self.geo_cond_change_times_dict[geo_cond_name] = list()
        self.geo_cond_change_back_times_dict[geo_cond_name] = list()

        # either 1d or 2d list, depending on
        # if within or between feature condition
        geo_cond_bits_list = self.geo_cond_bit_dict[geo_cond_name]

        for region1_str_or_list in geo_cond_bits_list:
            region1_times_list = list()
            region1_back_times_list = list()

            # between
            if type(region1_str_or_list) == list:
                # region2_str is a bit that will be in the order
                # epochs are specified in the age summary file (which
                # should be time index 1 -- top of file -- for the youngest
                # epoch, and a larger time index for oldest epochs)
                for region2_str in region1_str_or_list:
                    # make it old -> young
                    region2_str = region2_str[::-1]

                    # scanning the geographic condition bit pattern
                    # (sliding window of size 2) for a '01' string
                    # which indicates when the geographic condition
                    # manifested -- note that each position in the
                    # bit pattern is a different epoch
                    #
                    # we get the position in the bit pattern where
                    # we had the departing '0', and use that position
                    # to grab the start of the epoch with the '1'
                    # [age_starts_young2old[idx + 1] for idx, (i, j)
                    region2_times_list = \
                        [age_starts_old2young[idx + 1] for idx, (i, j) \
                           in enumerate(zip(region2_str,
                                            region2_str[1:])) \
                         if (i, j) == ('0', '1')]

                    # now doing change back times ('10')
                    # [age_starts_young2old[idx + 1] for idx, (i, j)
                    region2_back_times_list = \
                        [age_starts_old2young[idx + 1] for idx, (i, j) \
                         in enumerate(zip(region2_str,
                                          region2_str[1:])) \
                         if (i, j) == ('1', '0')]
                    
                    region1_times_list.append(region2_times_list)
                    region1_back_times_list.append(region2_back_times_list)

                self.geo_cond_change_times_dict[geo_cond_name].\
                    append(region1_times_list)

                self.geo_cond_change_back_times_dict[geo_cond_name].\
                    append(region1_back_times_list)

            # within
            else:
                # old to young bit
                region1_times_list.append(
                    [age_starts_old2young[idx+1] for idx, (i, j) \
                     in enumerate(zip(region1_str_or_list,
                                      region1_str_or_list[1:])) \
                                        if (i, j) == ('0', '1')]
                )

                region1_back_times_list.append(
                    [age_starts_old2young[idx + 1] for idx, (i, j) \
                     in enumerate(zip(region1_str_or_list,
                                      region1_str_or_list[1:])) \
                     if (i, j) == ('1', '0')]
                )

                self.geo_cond_change_times_dict[geo_cond_name].\
                    append(region1_times_list)

                self.geo_cond_change_back_times_dict[geo_cond_name]. \
                    append(region1_back_times_list)

    def _populate_oldest_geo_cond_bit_dict(self, geo_cond_name: str) -> None:
        """Populate class member holding oldest geog. condition bit.

        This method has the side-effect of populating
        self.geo_oldest_cond_bit_dict, whose keys are names of the
        geographic condition, and values are whether the condition is
        met in each region or region pair, at the oldest epoch.

        Args:
            ggeo_cond_name (str): Name of the geographic condition
                (e.g., 'mountain range', 'ancient sea',
                'riverine barrier').
        """

        # adding initial values to dictionary
        self.geo_oldest_cond_bit_dict[geo_cond_name] = list()

        # either 1d or 2d list, depending on
        # if within or between feature condition
        geo_cond_bits_list = self.geo_cond_bit_dict[geo_cond_name]

        for region1_str_or_list in geo_cond_bits_list:
            region1_oldest_bit_list = list()

            # between
            if type(region1_str_or_list) == list:
                for region2_str in region1_str_or_list:                    
                    # region2_oldest_bit = region2_str[0]
                    region2_oldest_bit = region2_str[-1]

                    region1_oldest_bit_list.append(region2_oldest_bit)

                self.geo_oldest_cond_bit_dict[geo_cond_name].append(region1_oldest_bit_list)

            # within
            else:
                # time directionality is given by user
                # region1_oldest_bit_list = region1_str_or_list[0]
                region1_oldest_bit_list = region1_str_or_list[-1]

                self.geo_oldest_cond_bit_dict[geo_cond_name].append(region1_oldest_bit_list)

    def _populate_conn_graph_list(self,
                                  geo_cond_name: str,
                                  is_directed: bool = False,
                                  user_young2old: bool = True):
        for idx, ep_age in \
                enumerate(self.feat_coll.epoch_age_start_list_old2young):

            g = GeoGraph(self.n_regions)

            for from_region_idx in range(self.n_regions):
                for to_region_idx in range(self.n_regions):
                    if from_region_idx != to_region_idx:
                        cond_change_bit_patt = \
                            self.geo_cond_bit_dict[geo_cond_name]\
                                [from_region_idx][to_region_idx]

                        # geog condition is met in this time slice
                        if cond_change_bit_patt[idx] == '1':
                            # for now assume that condition met here means pair is connected
                            g.add_edge(from_region_idx, to_region_idx, is_directed=is_directed)

            g.populate_comm_class_members()

            # one graph per epoch, and we want it from
            # old to young (so if user provided young -> old, we want
            # to reverse it)
            if user_young2old:
                self.conn_graph_list.insert(0, g)

            else:
                self.conn_graph_list.append(g)



    def find_epoch_idx(self, an_age: float) -> int:
        """

        Args:
            an_age (float): The age of an event.

        Returns:
            (int): Index of epoch containing 'an_age' (with 0 being
                present).
        """

        time_slice_index = 0
        for time_slice_index in range(0, self.n_time_slices):
            # old -> young
            time_slice_age_end = \
                self.feat_coll.epoch_age_end_list_old2young[time_slice_index]

            if an_age < time_slice_age_end or \
                    (abs(an_age - time_slice_age_end) <= 1e-12):
                continue

            break

        return time_slice_index

    def get_comm_classes(self, an_age: float) -> ty.List[ty.Set[int]]:
        epoch_idx = self.find_epoch_idx(an_age)
        g = self.conn_graph_list[epoch_idx]

        if len(g.comm_class_set_list) == 0:
            g.populate_comm_class_members()

        return g.comm_class_set_list


    # getters
    def get_geo_condition_change_times(self,
                                       geo_cond_name: str) -> \
            ty.Dict[str,
                    ty.Union[ty.List[ty.List[float]],
                    ty.List[ty.List[ty.List[float]]]]]:

        if not self.geo_cond_change_times_dict[geo_cond_name]:
            self._populate_geo_cond_change_times_dict(geo_cond_name)

        return self.geo_cond_change_times_dict[geo_cond_name]

    def get_geo_condition_change_back_times(self,
                                            geo_cond_name: str) -> \
            ty.Dict[str,
                    ty.Union[ty.List[ty.List[float]],
                    ty.List[ty.List[ty.List[float]]]]]:

        if not self.geo_cond_change_back_times_dict[geo_cond_name]:
            self._populate_geo_cond_change_times_dict(geo_cond_name)

        return self.geo_cond_change_back_times_dict[geo_cond_name]
    
    def get_geo_oldest_condition_bit(
            self,
            geo_cond_name):

        if not self.geo_oldest_cond_bit_dict[geo_cond_name]:
            self._populate_oldest_geo_cond_bit_dict(geo_cond_name)

        return self.geo_oldest_cond_bit_dict[geo_cond_name]


if __name__ == "__main__":
    # sys.argv[1]: feature summary file path
    # sys.argv[2]: age summary file path

    fc = GeoFeatureCollection(sys.argv[1],
                              age_summary_fp=sys.argv[2])

    fq = GeoFeatureQuery(fc)

    g = GeoGraph(name_list=['A', 'B', 'C', 'D'])

    test_basic_stuff = False
    if test_basic_stuff:
        # (1) can do it like this
        # print(fc.get_feat_by_name("qb_1", time_idx=1))
        # feat_list = fc.get_feat_by_name("qb_1")
        # for feat in feat_list:
        #     print(feat)
        
        # (2) or like this
        print(fc.get_feat_by_type_rel_id(
            GeoFeatureType.QUANTITATIVE,
            GeoFeatureRelationship.BETWEEN,
            1,
            time_idx=1))
        feat_list = fc.get_feat_by_type_rel_id(
            GeoFeatureType.QUANTITATIVE,
            GeoFeatureRelationship.BETWEEN,
            1)
        for feat in feat_list:
            print(feat)
        
        # checking getter
        # print(feat_list[0].get_val_from_to(1, 2))
    
    # testing querying
    test_querying = False
    if test_querying:
        requirement_fn1 = \
            GeoFeatureQuery.cw_feature_equals_value(fc, 1, feat_name="cw_1")
        
        # for all regions, gives all times when it happened
        fq.populate_geo_cond_member_dicts("ancient_sea", requirement_fn1)
        print("\nancient_sea 1:")
        print(" ".join(fq.geo_cond_bit_dict["ancient_sea"]))
        # ancient_sea 1:
        # 101 101 010 010

        print("\nancient_sea times 1:")
        for k in fq.geo_cond_change_times_dict["ancient_sea"]:
            print(*k)
        # ancient_sea times 1:
        # [5.0]
        # [5.0]
        # [15.0]
        # [15.0]

        print("\nancient_sea oldest bit 1:")
        print(" ".join(fq.geo_oldest_cond_bit_dict["ancient_sea"]))
        # ancient_sea oldest bit 1:
        # 1 1 0 0

        requirement_fn1_1 = \
            GeoFeatureQuery.cw_feature_equals_value(fc, 1, feat_id=1)
        
        # for all regions, gives all times when it happened
        fq.populate_geo_cond_member_dicts("ancient_sea", requirement_fn1_1)
        print("\nancient_sea 2:")
        print(" ".join(fq.geo_cond_bit_dict["ancient_sea"]))
        # ancient_sea 2:
        # 101 101 010 010

        print("\nancient_sea times 2:")
        for k in fq.geo_cond_change_times_dict["ancient_sea"]:
            print(*k)

        print("\nancient_sea oldest bit 2:")
        print(" ".join(fq.geo_oldest_cond_bit_dict["ancient_sea"]))
        # ancient_sea oldest bit 2:
        # 1 1 0 0

        # for all regions, gives all times when it happened
        thresh = 25.0
        requirement_fn2 = \
            GeoFeatureQuery.qw_feature_threshold(fc, thresh, True, feat_name="qw_1")
        
        fq.populate_geo_cond_member_dicts("altitude", requirement_fn2)
        print("\naltitude:")
        print(" ".join(fq.geo_cond_bit_dict["altitude"]))
        # altitude:
        # 010 010 101 101

        print("\naltitude times:")
        for k in fq.geo_cond_change_times_dict["altitude"]:
            print(*k)

        print("\naltitude oldest bit:")
        print(" ".join(fq.geo_oldest_cond_bit_dict["altitude"]))
        # altitude oldest bit:
        # 0 0 1 1

        # for all pairs of regions, gives all times when it happened
        requirement_fn3 = \
            GeoFeatureQuery.cb_feature_equals_value_is_connected(fc, 0, feat_name="cb_1")
        
        fq.populate_geo_cond_member_dicts("land_bridge", requirement_fn3)
        print("\nland_bridge 1:")
        for k in fq.geo_cond_bit_dict["land_bridge"]:
            print(*k)
        # land_bridge 1:
        # 111 000 010 111
        # 000 111 010 111
        # 111 010 111 010
        # 111 111 010 111

        print("\nland_bridge times 1:")
        for k in fq.geo_cond_change_times_dict["land_bridge"]:
            print(*k)
        # land_bridge times 1:
        # [] [] [15.0] []
        # [] [] [15.0] []
        # [] [15.0] [] [15.0]
        # [] [] [15.0] []

        print("\nland_bridge oldest bit 1:")
        for k in fq.geo_oldest_cond_bit_dict["land_bridge"]:
            print(*k)
        # land_bridge oldest bit 1:
        # 1 0 0 1
        # 0 1 0 1
        # 1 0 1 0
        # 1 1 0 1

        # for all pairs of regions, gives all times when it happened
        requirement_fn3_1 = \
            GeoFeatureQuery.cb_feature_equals_value_is_connected(fc, 0, feat_id=1)
        
        fq.populate_geo_cond_member_dicts("land_bridge", requirement_fn3_1)
        print("\nland_bridge 2:")
        for k in fq.geo_cond_bit_dict["land_bridge"]:
            print(*k)
        # land_bridge 2:
        # 111 000 010 111
        # 000 111 010 111
        # 111 010 111 010
        # 111 111 010 111

        print("\nland_bridge times 2:")
        for k in fq.geo_cond_change_times_dict["land_bridge"]:
            print(*k)
        # land_bridge times 2
        # [] [] [15.0] []
        # [] [] [15.0] []
        # [] [15.0] [] [15.0]
        # [] [] [15.0] []

        print("\nland_bridge oldest bit 2:")
        for k in fq.geo_oldest_cond_bit_dict["land_bridge"]:
            print(*k)
        # land_bridge oldest bit 2:
        # 1 0 0 1
        # 0 1 0 1
        # 1 0 1 0
        # 1 1 0 1

        # for all pairs of regions, gives all times when it happened
        thresh = 15.0
        requirement_fn4 = \
            GeoFeatureQuery.qb_feature_threshold(fc, thresh, True, feat_name="qb_1")
        
        fq.populate_geo_cond_member_dicts("distance", requirement_fn4)
        print("\ndistance 1:")
        for k in fq.geo_cond_bit_dict["distance"]:
            print(*k)
        # distance 1:
        # 000 010 101 111
        # 010 000 101 000
        # 101 101 000 101
        # 000 000 101 000

        print("\ndistance barriers times 1:")
        for k in fq.geo_cond_change_times_dict["distance"]:
            print(*k)
        # distance barriers times 1:
        # [] [15.0] [5.0] []
        # [15.0] [] [5.0] []
        # [5.0] [5.0] [] [5.0]
        # [] [] [5.0] []

        print("\ndistance oldest bit 1:")
        for k in fq.geo_oldest_cond_bit_dict["distance"]:
            print(*k)
        # distance oldest bit 1:
        # 0 0 1 1
        # 0 0 1 0
        # 1 1 0 1
        # 0 0 1 0
         
        requirement_fn4_1 = \
            GeoFeatureQuery.qb_feature_threshold(fc, thresh, True, feat_id=2)
        
        fq.populate_geo_cond_member_dicts("distance", requirement_fn4_1)
        print("\ndistance 2:")
        for k in fq.geo_cond_bit_dict["distance"]:
            print(*k)
        # distance 2:
        # 000 101 010 110
        # 101 000 010 000
        # 010 010 000 010
        # 000 000 010 000
        
        print("\ndistance barriers times 2:")
        for k in fq.geo_cond_change_times_dict["distance"]:
            print(*k)
        # distance barriers times 2:
        # [] [5.0] [15.0] [15.0]
        # [5.0] [] [15.0] []
        # [15.0] [15.0] [] [15.0]
        # [] [] [15.0] []

        print("\ndistance oldest bit 2:")
        for k in fq.geo_oldest_cond_bit_dict["distance"]:
            print(*k)
        # distance oldest bit 2:
        # 0 1 0 0
        # 1 0 0 0
        # 0 0 0 0
        # 0 0 0 0

    test_graph = True
    if test_graph:
        requirement_fn = \
            GeoFeatureQuery.cb_feature_equals_value_is_connected(fc, 0, feat_id=1)

        # all members, including graph, are populated here!
        fq.populate_geo_cond_member_dicts("land_bridge", requirement_fn)


