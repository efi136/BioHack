import matplotlib as plt
import networkx as nx
import numpy as np
from Bio import Phylo
from io import StringIO

class Tree:
    def __init__(self, left = None, right = None, transition_matrix=None):
        self.parent = None
        self.left  = left
        self.right = right
        left.parent = self
        right.parent = self
        self.transition_matrix = transition_matrix
        self.name = self.generate_parent_name()

    def generate_parent_name(self):
        new_name = ''
        dct = self.transition_matrix
        for i in range(len(self.left.name)):
            l_letter = self.left.name[i]
            r_letter = self.right.name[i]
            max_score = 0
            for key in dct.keys():
                cur_score = dct[key][l_letter] * dct[key][r_letter]
                if cur_score > max_score:
                    max_score = cur_score
                    max_letter = key
            new_name = new_name + max_letter
        return new_name

    def get_graph(self, g=None):
        if g is None:
            g = nx.Graph()
        self.left.get_graph(g)
        self.right.get_graph(g)
        g.add_node(self.name)
        g.add_edge(self.left.name, self.name)
        g.add_edge(self.right.name, self.name)
        return g

    def get_phylo_string(self):
        return f'({self.left.get_phylo_string()}, {self.right.get_phylo_string()})'

    def get_phylo_graph(self):
        phylo_str = self.get_phylo_string()
        phylo_f = StringIO(phylo_str)
        return Phylo.read(phylo_f, "newick")
    
    def draw_graph(self):
        g = self.get_graph()
        nx.draw(g)
        plt.pyplot.show()

    
    def draw(self):
        tree = self.get_phylo_graph()
        tree.rooted = True
        Phylo.draw(tree)

class Leaf(Tree):
    def __init__(self, name):
        self.name = name
    
    def __repr__(self):
        return self.name
    
    def __eq__(self, other):
        return self.name == other.name

    def get_graph(self, g=None):
        if g is None:
            g = nx.Graph()
        g.add_node(self.name)
        return g
    
    def get_phylo_string(self):
        return self.name