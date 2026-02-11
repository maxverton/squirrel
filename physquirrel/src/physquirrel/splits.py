import itertools

from .utils import Partition, validate, id_generator

############################################################

class Split(Partition):
    """
    Class for 2-partitions of sets, child-class of the general Partition class. 
    Takes as input two sets of elements that form the split.
        self.elements: set containing elements of both sides
        self.set1, self.set2: sets making up the two sides of the split.
        self.weight: weight value associated to the split (default 1)
        self.penalty: penalty value associated to the split (default 0)
    """
    
    def __init__(self, set1:set, set2:set, weight = 1, penalty = 0):
        self.set1 = set1
        self.set2 = set2
        self.elements = self.set1 | self.set2
        self.parts = [frozenset(self.set1), frozenset(self.set2)]
        self.penalty = penalty
        self.weight = weight

        if validate():
            if not super()._is_valid():
                raise ValueError("Invalid split: sets overlap")
    
    def __repr__(self):
        return f"Split({set(self.set1)}, {set(self.set2)})"
    
    def is_trivial(self) -> bool:
        """Returns whether this is a trivial split."""
        return len(self.set1) == 1 or len(self.set2) == 1

    def is_subsplit(self, other) -> bool:
        """Returns whether the split is a subsplit of 'other'. Raises an error 
        if 'other' is not of type Split."""
        if validate():
            if not isinstance(other, Split):
                raise ValueError("Not a Split instance")
        return (self.set1.issubset(other.set1) and self.set2.issubset(other.set2)) or \
                (self.set1.issubset(other.set2) and self.set2.issubset(other.set1))
        
    def induced_quartetsplits(self, include_trivial=False):
        """Returns a QuartetSplitSet of all subsplits of size 4 of the split."""
        res = []
        for s1 in itertools.combinations(self.set1, 2):
            for s2 in itertools.combinations(self.set2, 2):
                split = Split(set(s1), set(s2))
                res.append(split)
                
        if include_trivial == True:
            for s1 in itertools.combinations(self.set1, 1):
                for s2 in itertools.combinations(self.set2, 3):
                    split = Split(set(s1), set(s2))
                    res.append(split)
            
            for s1 in itertools.combinations(self.set1, 3):
                for s2 in itertools.combinations(self.set2, 1):
                    split = Split(set(s1), set(s2))
                    res.append(split)
                    
        return QuartetSplitSet(set(res))

############################################################






############################################################

class SplitSet:
    """
    Class for a set of splits. Takes as input a set of splits.
        self.splits: set of the splits
        self.elements: set containing all elements appearing in the splits
    """
    
    def __init__(self, splits):
        self.splits = set(splits)
        self.elements = set()
        for split in self.splits:
            self.elements.update(split.elements)
    
    def __repr__(self):
        return f"SplitSet({list(self.splits)})"
    
    def __iter__(self):
        return iter(self.splits)
    
    def __contains__(self, split):
        return split in self.splits    
    
    def __len__(self):
        return len(self.splits)  

############################################################




############################################################

class SplitSystem(SplitSet):
    """
    Class for a split system: set of full splits (complete partitions of elements);
    child class of the SplitSet class. Raises an error if not all splits cover
    the complete set of elements. Does not check if the system 
    is compatible.
        self.splits: set of the splits
        self.elements: set containing all elements appearing in the splits
    """
    
    def __init__(self, splits):
        self.splits = set(splits)
        if len(self.splits) == 0:
            self.elements = set()
        else:
            self.elements = list(splits)[0].elements

        if validate():
            if not self._is_valid():
                raise ValueError("Not a set of full splits.")
    
    def add(self, split):
        """Adds a split to the system. Raises an error if the split is not a full split on self.elements."""
        if validate():
            if split.elements != self.elements:
                raise ValueError("Split does not cover all elements.")
        self.splits.add(split)
    
    def remove(self, split):
        """Removes a split from the system. Raises an error if the split is not in the system."""
        if validate():
            if split not in self.splits:
                raise ValueError("Split not in system.")
        self.splits.remove(split)

    def __repr__(self):
        return f"SplitSystem({list(self.splits)})"

    def _is_valid(self):
        for split in self.splits:
            if split.elements != self.elements:
                return False
        return True
    
    def displayed_tree(self):
        """Returns the tree that displays the splitsystem. Raises an error if
        no such tree exists (i.e. the system is incompatible)."""
        
        from .sdnetwork import SemiDirectedNetwork
        
        tree = SemiDirectedNetwork()
        center_node = id_generator()
        tree.add_node(center_node)
        tree.add_leaves_from(self.elements)
        tree.add_edges_from([(center_node,leaf) for leaf in self.elements])
        
        for split in self.splits:
            if not split.is_trivial():
                tree.create_split(split)     
        
        return tree
    
    def induced_quartetsplits(self):
        """Returns a QuartetSplitSet of all subsplits of size 4 of all the splits in the set."""
        res = set()
        for split in self.splits:
            res = res | split.induced_quartetsplits().splits
        return QuartetSplitSet(res, elements=self.elements)
    
    def list_sorted_penalties(self):
        """Returns a list of the splits sorted by their penalty values (lowest
        penalty first)."""
        return sorted(self.splits, key=lambda split: split.penalty)
    
    def remove_trivial_splits(self):
        """Removes all trivial splits from the system."""
        self.splits = {split for split in self.splits if not split.is_trivial()}

        
    def build_trees_bstar_greedy(self):
        """Returns the trees built from the B*-set of the splitsystem. Raises an error if no such tree exists."""
        self.remove_trivial_splits()
        split_penalty_list = self.list_sorted_penalties()
        
        tree_list = []
        #Set up all the trivial splits first
        tree_split_system = SplitSystem([Split({leaf}, self.elements - {leaf}) for leaf in self.elements])
        for split in split_penalty_list:
            tree_split_system.add(split)
            try:
                tree = tree_split_system.displayed_tree()
                tree_list.append({'tree' : tree, 'split_system' : tree_split_system})
            except ValueError:
                tree_split_system.remove(split)
        return tree_list

            
############################################################




############################################################

class QuartetSplitSet(SplitSet):
    """
    Class for a set of (non-trivial) quartet-splits (2|2 splits of some large 
    set of elements); child class of the SplitSet class. The optional
    argument is the total set of elements which should be a superset of the elements of the splits.
    Raises an error otherwise. If elements=None, it is assumed to be self.elements.
    Raises an error if not all splits are quartet-splits, or if more than one quartet-split exists on 
    the same element set.
        self.splits: set of the splits
        self.elements: set containing all elements appearing in the splits
    """
    
    def __init__(self, splits, elements=None):
        super().__init__(splits)
        if elements is None: elements = self.elements
        
        if validate():
            if not self.elements.issubset(elements):
                raise ValueError
            if not self._is_valid():
                raise ValueError("Not a set of quartet-splits.")        
            
        self.elements = elements

    def __repr__(self):
        return "Quartet" + super().__repr__()
    
    def _is_valid(self):
        if not hasattr(self, '_split_dict'):
            self._split_dict = {frozenset(split.elements):split for split in self.splits}

        covered_elements = []
        for split in self.splits:
            if len(split.set1) != 2 or len(split.set2) != 2:
                return False
            covered_elements.append(frozenset(split.elements))
        return len(covered_elements) == len(set(covered_elements))
    
    def obtain_split(self, four_elements):
        """Returns the split on the given four leaves, if it exists. Otherwise,
        returns None. Raises an error if the argument is not a size four subset
        of self.elements."""
        if validate():
            if len(four_elements) != 4:
                raise ValueError("Method needs size 4 subset.")
            if not four_elements.issubset(self.elements):
                raise ValueError("Argument is not a subset of self.elements.")
        
        if not hasattr(self, '_split_dict'):
            self._split_dict = {frozenset(split.elements):split for split in self.splits}

        if frozenset(four_elements) in self._split_dict.keys():
            return self._split_dict[frozenset(four_elements)]
        else:
            return None
     
    def bstar(self):
        """Returns the B*-set of the QuartetSplitSet as a SplitSystem.
        Uses the incremental O(n^5) algorithm from 
        ---  Berry, Vincent, and Olivier Gascuel. "Inferring evolutionary trees  ---
        ---  with strong combinatorial evidence." Theoretical computer science   ---
        ---  240.2 (2000): 271-298.                                              ---
        """
        order = list(self.elements)
        a, b, c, d = order[0:4]
        bstar = [Split({a},{b,c,d}),Split({b},{a,c,d}),Split({c},{a,b,d}),Split({d},{a,b,c})]
        abcd_split = self.obtain_split({a,b,c,d})
        if abcd_split is not None: bstar.append(abcd_split)
        
        for i, element in enumerate(order):
            if i < 4: continue
            new_bstar = [Split({element},set(order[0:i]))]
            for split in bstar:
                candidate_split1 = Split(split.set1 | {element}, split.set2)
                candidate_split2 = Split(split.set1, split.set2 | {element})

                add1 = not any(
                    self.obtain_split({element, x, y, z}) != Split({x, element}, {y, z})
                    for x in split.set1 for y, z in itertools.combinations(split.set2, 2))
                
                add2 = not any(
                    self.obtain_split({element, x, y, z}) != Split({x, y}, {z, element})
                    for x, y in itertools.combinations(split.set1, 2) for z in split.set2)
                
                if add1 == True:
                    new_bstar.append(candidate_split1)
                if add2 == True:
                    new_bstar.append(candidate_split2)
                
                bstar = new_bstar

        return SplitSystem(bstar)
    
    def bstar_penalty(self, threshold=3):
        taxa = list(self.elements)
        a, b, c, d = taxa[0:4]
        bstar = [Split({a},{b,c,d}),Split({b},{a,c,d}),Split({c},{a,b,d}),Split({d},{a,b,c})]
        abcd_splits = [Split({a,b},{c,d}), Split({a,c},{b,d}), Split({a,d},{b,c})]
        abcd_split = self.obtain_split({a,b,c,d})

        if abcd_split in abcd_splits:
            for split in abcd_splits:
                if split != abcd_split:
                    split.penalty = abcd_split.weight

        bstar.extend(abcd_splits)
        for i, element in enumerate(taxa):
            if i < 4: continue
            new_bstar = [Split({element},set(taxa[0:i]))]

            
            for split in bstar:
                candidate_split1 = Split(split.set1 | {element}, split.set2, penalty = split.penalty)
                candidate_split2 = Split(split.set1, split.set2 | {element}, penalty = split.penalty)

                ### For the split where the new element is added to set1 we compute the penalty for the
                # newly induced quartet trees by the split containing the new element.  
                for x in split.set1:
                    for y, z in itertools.combinations(split.set2, 2):
                        q_split = self.obtain_split({element, x, y, z})
                        if q_split is None: continue

                        if q_split != Split({x, element}, {y, z}):
                            candidate_split1.penalty += q_split.weight

                for x, y in itertools.combinations(split.set1, 2):
                    for z in split.set2:
                        q_split = self.obtain_split({element, x, y, z})
                        if q_split is None: continue

                        if q_split != Split({x, y}, {z, element}):
                            candidate_split2.penalty += q_split.weight

                if candidate_split1.penalty < threshold:
                    new_bstar.append(candidate_split1)
                if candidate_split2.penalty < threshold:
                    new_bstar.append(candidate_split2)
                
                bstar = new_bstar

        return SplitSystem(bstar)

    def get_blobtrees_tstar_penalty(self, threshold=3):
        bstar_splitsystem = self.bstar_penalty(threshold=threshold)
        for split in bstar_splitsystem.splits:
            print(f"Split: {split}, Penalty: {split.penalty}")
        tstar_list = bstar_splitsystem.build_trees_bstar_greedy()
        return tstar_list
