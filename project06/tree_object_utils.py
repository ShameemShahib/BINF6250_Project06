from collections import defaultdict
import numpy as np

class Node:
    """A class to represent the nodes on our tree structure"""
    def __init__(self, name: str, branch_length: float) -> None:
        self.name = name
        self.branch_length = branch_length
        self.parent_name = "default_parent"

    def __repr__(self) -> str:
        """string representation of nodes in such a way that makes Newick representation esier to generate"""
        return f"{self.name}:{round(self.branch_length, 3)}"


class Tree:
    """
    The class representing the tree graph structure

    Attributes
    ----------
    self.nodes: dict[str, Node]
        A dictionary mapping a node's name to the node object representing it (enables Node lookup by name)
    self.top_layer: list[Node]
        List of nodes that currently have no parent
    self.edges: dict[Node, list[Node]]
        Dictionary mapping a parent node to a list of its children

    Methods
    -------
    add_node(new_node: Node) -> None
        Takes a Node object as input and adds it to the tree without any connections 
        (creates entry in self.nodes and adds it to the top layer)
    get_parentless_nodes() -> list[Node]
        Identifies nodes on the graph that have no parents (and therefore belong on self.top_layer)
    make_parent_name(children: list[Node]) -> str
        generate the name for an internal node based on what its children are
    create_parent(children: list[Node], dist_to_outgroup: float) -> str
        create an internal node and its connections, adding it to the graph 
        returns the name of the new internal node so the Neighbor Joining function can access it from self.nodes
    """
    def __init__(self, nodes_dict: dict[str, Node] = {},
                 top_layer: list[Node] = [],
                 node_mapping: dict[Node, list[Node]] = {}) -> None:
        # completely overhauled datastructure (make copies to not get owned by mutable defaults)
        self.nodes = nodes_dict.copy()
        self.top_layer = top_layer.copy()
        self.edges = node_mapping.copy()
    

    def add_node(self, new_node: Node) -> None:
        """Helper method to add a Node to the graph with no connections"""
        # Give the node an entry in the dict that lets us look up nodes by name
        self.nodes[new_node.name] = new_node
        # now add the node to the top layer of the graph
        self.top_layer.append(new_node)


    def get_parentless_nodes(self) -> list[Node]:
        """
        Function to identify nodes on the graph that have no parents (and therefore belong on self.top_layer)

        Returns
        -------
        list[Node]
            list of nodes that have no edges going into them (and no parents)
        """

        # make a list of just the node objects (we don't care about their names)
        parentless_nodes = [node for node in self.nodes.values()]

        # iterate through the edges in self.edges
        for edge_list in self.edges.values():
            # edge_list will be a list of Nodes that are the children of... something (it doesn't matter to us what their parents are)
            for child in edge_list:
                # check if the current child is still in the parentless nodes list
                if child in parentless_nodes:
                    # since we know this node appears as the child of something, it can't be parentless
                    parentless_nodes.remove(child)
        
        # return the list of whatever nodes remain after removing the ones that are listed as the child of literally anything at any level
        return parentless_nodes
    
    
    def make_parent_name(self, children: list[Node]) -> str:
        """
        Function to generate the name for an internal node based on what its children are

        NOTE: Yes, this is a shortcut to bypass recursively traversing the tree to make the newick string for the whole tree 
              This works entirely off the good-will of the Neighbor-Joining algorithm building the tree from bottom to top

        Parameters
        ----------
        children : list[Node]
            list of children that will descend from the node that will be given this name

        Returns
        -------
        str
            the name for an internal node, which will happen to also be the Newick string representation of the 
            subtree descending from this internal node
        """
        # initialize the parent name
        parent_name = "("

        # concatenate the names of the children to the name
        for child_node in children[:-1]:
            # add this child's name + a comma
            parent_name += str(child_node) + ","
        
        # now add the last child's name (no comma needed)
        return parent_name + str(children[-1]) + ")"

    
    def create_parent(self, children: list[Node], dist_to_outgroup: float) -> str:
        """
        Function to create an internal node and its connections, adding it to the graph

        Parameters
        ----------
        children : list[Node]
            the list of Nodes that this one will be parent to
        dist_to_outgroup : float
            Distance to ANY leaf that doesn't descend to this one, used as an approximation for this node's limb length
            This gets updated later if we make a node that's parental to this one, but otherwise this value will be our best guess

        Returns
        -------
        str
            The name of this new internal node
            This only gets returned so the Neighbor Joining algorithm can access this node by name later
        """

        # make the name for this parent
        parent_name = self.make_parent_name(children)

        # create a node object
        parent = Node(parent_name, dist_to_outgroup)
        self.add_node(parent)

        # connect this parent to its children
        self.edges[parent] = children

        # iterate through the nodes that this internal node will be parental to
        for child in children:
            # attempt to remove this child from the top layer
            try:
                self.top_layer.remove(child)
            except ValueError:
                print(f"Attempted to remove Node {child.name} from tree's top layer but it isn't in the top layer")
            
            # update this child's name for its parent
            child.parent_name = parent_name

        # now return the name for the parent so the Neighbor Joining func can access this node by name
        return parent_name
   

    def __repr__(self) -> str:
        """This should give us the newick string representation of our tree structure"""
        # Since the internal nodes' names are Newick strings representing everything that descends from them, we just find
        # the nodes that don't have any parents and use them as the basis for the whole-tree Newick-String
        if self.top_layer:
            return self.make_parent_name(self.top_layer) + ";"
        else:
            # this only runs if the tree is empty
            return "();"
