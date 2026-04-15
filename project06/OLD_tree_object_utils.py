from collections import defaultdict
import numpy as np

class Node:
    """A class to represent the nodes on our tree structure"""
    def __init__(self, name: str, branch_length: float) -> None:
        self.seq_id = name
        self.branch_length = branch_length
        self.parent_name = name + "_parent"
        self.balded_distances = {}

    def __repr__(self) -> str:
        """string representation of nodes in such a way that makes Newick representation esier to generate"""
        return f"{self.seq_id}: {round(self.branch_length, 4)}"


class Tree_Graph:
    """The constructor for our graph structure"""
    def __init__(self, distance_matrix: np.matrix) -> None:
        self.nodes = {}  # dict mapping node names to node objects
        self.edges = defaultdict(list)  # dict mapping parent node to list of children nodes
        self.edges_in = defaultdict(int)  # maps node object to number of edges going into it
        self.edges_out = defaultdict(int)  # maps node object to number of edges going out of it
        self.leaves = []  # list of terminal nodes
        self.matrix_labels = []  # list of seqIDs in order they appear on the distance matrix
        self.distance_matrix = distance_matrix  # is a numpy matrix so no named rows or columns
        self.last_nodes = []  # this is tracking the last two sequences in the distance matrix after trimming

    def add_edge(self, parent_node: Node, child_node: Node) -> None:
        # create the association between the parent and child
        self.edges[parent_node].append(child_node)
        # update who the child thinks its parent is
        child_node.parent = parent_node.seq_id
        # update the balances for all the nodes
        self.edges_out[parent_node] += 1
        self.edges_in[child_node] -= 1

    def get_matrix_idx(self, node_name: str):
        """Gets the row/col number on the dist matrix for a node"""
        return self.matrix_labels.index(node_name)

    def get_leaf_branch_length(self, seq_one: str, other1: str, other2: str):
        """Use the distances between a sequence and 2 others """
        # get the rows/columns of each sequence being looked at
        node_idx = self.get_matrix_idx(seq_one)
        other1_idx = self.get_matrix_idx(other1)
        other2_idx = self.get_matrix_idx(other2)
        # get the branch length for 
        return (self.distance_matrix[node_idx, other1_idx] + self.distance_matrix[node_idx, other2_idx] - self.distance_matrix[other1_idx, other2_idx]) / 2

    def get_internal_branch_length(self, leaves: list[Node]) -> float | None:
        """
        get the branch length for an internal node, NOT FOR LEAVES
        
        @param leaves: list of nodes, these are the leaves that are descendents of the node
                       for which we're calculating the branch length. Don't try putting in only one leaf I swear to god
        @return: float, the branch length for whatever node is the ancestor to all the leaves in node 
        """

        # select a node on the distance matrix NOT in the list of nodes (create copy of self.matrix_labels and remove all members of nodes)
        # shallow-copy the list of sequences on the distance matrix
        outgroups = self.leaves.copy()

        # remove the leaves that descend from the node for which we're calculating branch length, as we need an outgroup to calc branch length
        for leaf in leaves:
            outgroups.remove(leaf)

        # grab two of the leaves from the param list
        leaf_one = leaves[0]
        leaf_two = leaves[-1]

        # initialize a list of distances to outgroups
        dists_to_outgroups = []

        # pick ANY outgroup, doesn't matter which
        try:
            for outgroup in outgroups:
                # calculate how far the outgroup is from the shared ancestor of the two leaves we selected
                dists_to_outgroups.append(self.get_leaf_branch_length(outgroup.name, leaf_one.seq_id, leaf_two.seq_id))

        except IndexError:
            # this only runs if there are no outgroups left to make
            print("Attempted to calculate branch length for an internal node that is ancestral to ALL leaves on the graph. Cannot compute branch length")
            return None
        
        # now, problem: we have the distances from our internal node to all leaves that aren't descended from it
        # now, we do branch length calc for the internal node using these values

        # grab two outgroups
        outgroup1 = outgroups[0]  # this is the same outgroup as the one used for dists_tooutgroups[0]
        outgroup2 = outgroups[-1]  # this is the same outgroup as the one used for dists_tooutgroups[-1]
        between_outgroups = self.get_dist(outgroup1, outgroup2)

        # do the branch length formula
        return (dists_to_outgroups[0] + dists_to_outgroups[-1] - between_outgroups) / 2
    

    def make_parent_name(self, children: list[Node]):
        """Function to generate the name for a node based on what it's children are"""
        parent_name = []
        for child in children:
            # append the string representation of the child node, which is formatted as a newick string
            parent_name.append(str(child))  # not to fear, str(child) will include the child's branch length
        
        # return the newick representation for this new internal node's name
        return "(" + ",".join(parent_name) + ")"
        

    def add_parent(self, children):
        """Function to make a node that is the parent to whatever nodes are in the given list called children"""
        # dude forgive me for what I'm about to do, I'm 16 hours into my day and sleep deprived
        parent_name = self.make_parent_name(children)

        # get the leaves that descend from this node (does this work? Idk man, it's past midnight)
        terminal_descendents = [leaf for leaf in self.leaves if leaf.name in self.matrix_labels]

        # attempt to calculate the branch length for the parent we're creating
        branch_length = self.get_internal_branch_length(terminal_descendents)

        # check that a branch length could be calculated
        if branch_length is not None:

            # create the node object
            new_parent = Node(parent_name, branch_length)

            # add this node to the graph
            self.nodes[new_parent.seq_id] = new_parent

            # add the connection between this new ancestor and its children
            for child in children:
                self.add_edge(new_parent, child)

        else:
            # tell the user that this parent can't be made into a node because it doesn't have a branch length
            print(f"This parent {parent_name} can't be made into a node object because it's ancestral to EVERYTHING ELSE on the tree, so I can't calculate branch length")


    def get_path_to_leaf(self, node: Node) -> tuple[Node, float]:
        """Function to traverse the graph to find a leaf that the given node is ancestral to"""

        # initialize our current position on the graph
        curr_position = node
        branches_hit = []

        # now we keep going until we reach a terminal node (skipped completely if we passed in a leaf)
        while curr_position not in self.leaves:
            # add this position to the list of branches hit along the way to the leaf
            branches_hit.append(curr_position)

            # access any one of the nodes that this one is parent to
            curr_position = self.edges[curr_position][0]
        
        # get the distance from the original ancestral node to the leaf we identified
        dist_to_leaf = sum(descendent.branch_length for descendent in branches_hit)

        return curr_position, dist_to_leaf
    

    def get_all_leaves(self, node: Node) -> tuple[Node]:
        """Function to get a list of all leaves that descend from node by walking the graph in a truly deranged manner"""

        # make a shallow copy of what nodes we need to hit, e.g. ABCD -> ABC and D so it's [ABC, D]
        nodes_to_hit = self.edges[node].copy()
        # make a list of the leaves descending from node
        leaves = []

        while nodes_to_hit != []:
            # access one of the nodes we want to hit, [ABC, D] --> access ABC
            curr_position = nodes_to_hit[0]
            # check if the child is a leaf
            if curr_position.name in self.leaves:
                # if the child IS a leaf, add it to the leaves list and remove it from the list of nodes to hit
                # e.g. nodes_to_hit is [D, AB, C]
                leaves.append(curr_position)  # e.g. now leaves is [D]
                nodes_to_hit.remove(curr_position)  # e.g. now nodes_to_hit is [AB, C]
            else:
                # if the child isn't a leaf, add its children to the list of nodes we need to check and remove the current node from the list
                # e.g. ABC --> AB and C, so now nodes_to_hit is [D, AB, C]
                nodes_to_hit.append(child for child in self.edges[curr_position])  # e.g. nodes_to_hit is [ABC, D, AB, C]
                nodes_to_hit.remove(curr_position)  # e.g. nodes_to_hit is now [D, AB, C]

        # return as tuple to save memory because we don't need to edit it later
        return tuple(leaves)
    
    

    def get_dist(self, node_one: Node, node_two: Node) -> float:
        """Function to get the distance between any two nodes, regardless of either or both being internal nodes"""
        # minus any branch lengths they hit. E.g. Internal node ABC to D will be matrix[A, D] - (A.branch_length + AB.branch_length)
        # identify the leaves we can use to find distances between these two nodes
        node_one_leaf, node_one_descendent_dist = self.get_path_to_leaf(node_one)
        node_two_leaf, node_two_descendent_dist = self.get_path_to_leaf(node_two)
        
        # get the row and column number for these two leaves being measured on the distance matrix
        des1_idx = self.get_matrix_idx(node_one_leaf.seq_id)
        des2_idx = self.get_matrix_idx(node_two_leaf.seq_id)
        # get the distance between the leaves that descend from node 1 and node 2
        leaf_distance = self.distance_matrix[des1_idx, des2_idx]

        # now we can get the distance between node 1 and node 2, regardless of if they're internal
        return leaf_distance - (node_one_descendent_dist + node_two_descendent_dist)
        

    def shares_parent(self, node_one: Node, node_two: Node) -> tuple[bool, bool]:
        """Function to do a quick check to make sure two nodes share a parent"""
        # if these nodes share parents, their distance will be (less than or) equal to the sum of their branch lengths
        sum_branch_lengths = node_one.branch_length + node_two.branch_length

        # Regardless of if these nodes share parents, their distance should be equal to the distance between their leaves
        distance = self.get_dist(node_one, node_two)

        # initialize the booleans we spit out
        share_parent = False
        update_branch_lengths = False

        if sum_branch_lengths == distance:
            # this one makes sense
            share_parent = True
        elif sum_branch_lengths > distance:
            # this only happens if two nodes share a parent and one of them calculated their branch length using
            # two other sequences on the distance matrix that were more related to each other than to it, so the internal node
            # that the balded distances are relative to is a more distant ancestor to one of them idk man
            share_parent = True
            update_branch_lengths = True

        return share_parent, update_branch_lengths
   

    def __repr__(self) -> str:
        """This should give us the newick string representation of our tree structure"""
        # actually since each node's string representation is conducive to newick format, just find the most ancestral nodes on the graph <3
        # "what if there are multiple equally ancestral nodes on the graph because it's unrooted?" Idk man pick one and then express the other as being some distance from it or something
        parentless_nodes = []
        for node in self.nodes.keys():
            if self.edges_in[node] == 0:
                parentless_nodes.append(node)
        
        # now since the last two parentless nodes can't have a parent node made for them, we just put them next to each other
        parent_name = self.make_parent_name(parentless_nodes)

        return parent_name