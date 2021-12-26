import matplotlib as plt
import networkx as nx

class Tree:
    def __init__(self, left = None, right = None):
        self.parent = None
        self.left  = left
        self.right = right
        left.parent = self
        right.parent = self
        self.name = f'{left.name}_{right.name}'

    def get_graph(self, g=None):
        if g is None:
            g = nx.Graph()
        self.left.get_graph(g)
        self.right.get_graph(g)
        g.add_node(self.name)
        g.add_edge(self.left.name, self.name)
        g.add_edge(self.right.name, self.name)
        return g
    
    def draw(self):
        g = self.get_graph()
        nx.draw(g)
        plt.pyplot.show()
    
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