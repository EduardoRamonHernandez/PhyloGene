"""
    Author: Eduardo Hernandez
    File: tree.py
    Description: This program creates tree objects to hold genome as leaves and
    their similar genomes.
"""


class Tree:
    """
    This class creates a leaf node for genomes or holds the leaves of genomes
    if it is the root of the genomes. Each node has a name if it has a genome,
    a reference to its left or right unless it is the most outer leaf, and all
    the genome under the node.
    """
    def __init__(self, name):
        self._name = name
        self._left = None
        self._right = None
        self._under = []

    def is_leaf(self):
        """
        This method returns if the node is a leaf
        return: True or False if the node is a leaf
        """
        return self.name() != None

    def name(self):
        """
        This method returns the name of the leaf node
        return: string representing the name of the leaf
        """
        return self._name

    def left(self):
        """
        This method returns the left side of the leaf node
        return: reference to the left of the leaf node
        """
        return self._left

    def right(self):
        """
        This method returns the right side of the leaf node
        return: reference to the right of the leaf node
        """
        return self._right

    def __str__(self):
        if self.is_leaf():
            return self.name()
        else:
            return "({}, {})".format(str(self.left()), str(self.right()))