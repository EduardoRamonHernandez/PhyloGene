"""
    Author: Eduardo Hernandez
    File: phylo.py
    Description: This program reads in a fasta file containing names of
    organisms and their genetic sequence then creates a tree containing all the
    organisms with the most similar being closest to each other. This is done
    by importing two ADT's: GenomeData that holds the information of each
    organism and Tree which creates tree nodes and connects the organisms. The
    program reads in an integer N to create n-grams of length N out of the
    genetic sequences of these organisms to determine the similarity between
    them.
"""
from genome import GenomeData
from tree import Tree


def seq_sim(set1, set2):
    """
    This function returns the Jaccard similarity of two sets of n-grams
    param set1: set of n-grams of one organism
    param set2: set of n-grams of another organism
    return: Jaccard similarity of these two sets
    """
    return len(set1.intersection(set2)) / len(set1.union(set2))


def read_file(fname):
    """
    This function reads in the file inputted by the user and returns the list
    of lines within the file
    param fname: string representing the name of the file
    return: list of lines within the file
    """
    file = open(fname, 'r')
    lines = file.readlines()
    file.close()
    return lines


def get_ngrams(s, n):
    """
    This function returns the set of n-grams of an organism. This is done by
    list comprehension using a set operator and indexing from the beginning of
    the sequence up to n and moving the index up one.
    param s: string representing the genetic sequence of an organism
    param n: integer representing the length of individual n-grams
    return: set of n-grams from a genetic sequence
    """
    return set(s[i:i+n] for i in range(len(s)-n+1))


def make_genomes(lines, n):
    """
    This function reads the lines from an input fasta file and creates
    GenomeData objects containing the name of the organism, the genetic
    sequence of that object, and the n-grams of that organism's genetic
    sequence.
    param lines: lines from an inputted fasta file
    param n: integer representing length of n-grams
    return: list of all Genomes in the file
    """
    # create list to add genomes to
    genomes = []
    # temporary list to get name of organisms and their data
    temp = []
    i = 0
    # while loop to go through the file
    while i < len(lines):
        # name of genome has > as first character in line
        if lines[i][0] == ">":
            info = lines[i].split()
            name = info[0][1:]
            temp.append(name)
            i += 1
        # if line is not a name, it is either an empty line or genetic sequence
        else:
            # string to add entire genetic sequence to
            sequence = ""
            # counter to advance through file correctly
            counter = 0
            # look at all likes through the rest of the file
            for lines[i] in lines[i:]:
                counter += 1
                # once empty line is reached, that is end of genetic sequence
                if lines[i].strip() == '':
                    break
                # add genetic sequence to the sequence variable without \n
                else:
                    sequence += lines[i].strip("\n")
            # if sequence has genetic sequence then create ngrams and genome
            if sequence != '':
                ngram = get_ngrams(sequence, n)
                genomes.append(GenomeData(temp[-1], sequence, ngram))
            i += counter
    return genomes


def find_similarities(genomes):
    """
    This function finds the similarities using the Jaccard Index between two
    sets of n-grams of different genomes. Both orders of the genomes are stored
    within the dictionary for ease later on.
    param genomes: list of all GenomeData objects that have been created
    return: dictionary containing the similarity between all pairs of genomes
    """
    # create dictionary for similarities
    sim = {}
    # loop through all genomes twice for every pair
    for genome in genomes:
        for g in genomes:
            # only look at pairs of different genomes
            if g != genome:
                # get similarities and add both orders of the pair to sim
                sim[(genome.name(), g.name())] = \
                    seq_sim(genome.ngrams(), g.ngrams())
                sim[(g.name(), genome.name())] = \
                    seq_sim(genome.ngrams(), g.ngrams())
    return sim


def make_tree(leaves, sim):
    """
    This function combines all the genomes into one Tree object by pairing the
    most similar genomes and adding the most similar genomes next to the tree.
    This is done by creating new trees that have the genomes as leaves and
    remove the old trees until only one tree is left.
    param leaves: list of all Genome Tree nodes
    param sim: dictionary of the similarities between all the genomes
    return: Tree node containing all genomes
    """
    # index through the trees until only one is left
    while len(leaves) > 1:
        # list to keep track of genomes with the highest similarity
        max_sim_leaves = []
        # integer to keep track of the highest similarity
        max_sim = 0
        # loop through all the pairs of tree nodes
        for curr in leaves:
            for leaf in leaves:
                # make sure the nodes are not the same
                if curr != leaf:
                    # if the trees have one genome, check their names in sim
                    if curr.is_leaf() and leaf.is_leaf():
                        # only keep track of them if they higher than max_sim
                        if sim[(curr.name(), leaf.name())] > max_sim:
                            max_sim = sim[(curr.name(), leaf.name())]
                            max_sim_leaves = [curr, leaf]
                    # if one tree node has multiple nodes, loop through nodes
                    elif not curr.is_leaf() and leaf.is_leaf():
                        for node in curr._under:
                            if sim[(node.name(), leaf.name())] > max_sim:
                                max_sim = sim[(node.name(), leaf.name())]
                                max_sim_leaves = [curr, leaf]
                    elif curr.is_leaf() and not leaf.is_leaf():
                        for node in leaf._under:
                            if sim[(curr.name(), node.name())] > max_sim:
                                max_sim = sim[(curr.name(), node.name())]
                                max_sim_leaves = [curr, leaf]
                    # if trees have multiple nodes, loop through their nodes
                    else:
                        for node in curr._under:
                            for tree in leaf._under:
                                if sim[(node.name(), tree.name())] > max_sim:
                                    max_sim = sim[(node.name(), tree.name())]
                                    max_sim_leaves = [curr, leaf]
        # create new tree once max similarity is found and set its name to None
        # this is done to say that the tree is not a leaf
        new = Tree(None)
        # check if the trees with the highest similarity are leaves
        # if they are, then add their trees to the list of leaves in new tree
        if max_sim_leaves[0].is_leaf():
            new._under.append(max_sim_leaves[0])
        if max_sim_leaves[1].is_leaf():
            new._under.append(max_sim_leaves[1])
        # if not, add all trees in either of the highest similarity trees
        else:
            for node in max_sim_leaves[0]._under:
                new._under.append(node)
            for leaf in max_sim_leaves[1]._under:
                new._under.append(leaf)
        # the smaller tree goes on the left of the new tree
        if str(max_sim_leaves[0]) < str(max_sim_leaves[1]):
            new._left = max_sim_leaves[0]
            new._right = max_sim_leaves[1]
        else:
            new._left = max_sim_leaves[1]
            new._right = max_sim_leaves[0]
        # add the new tree to the list
        leaves.append(new)
        # get rid of the old trees that were added to the new tree
        for leaf in max_sim_leaves:
            if leaf in leaves:
                leaves.remove(leaf)
    return leaves[0]


def main():
    # get name of fasta file and integer for n-grams
    fname = input('FASTA file: ')
    N = int(input('n-gram size: '))
    # get lines within the file
    lines = read_file(fname)
    # make all the genomes from organisms in the file
    genomes = make_genomes(lines, N)
    # get the similarities between every pair of genomes
    sim = find_similarities(genomes)
    # make every genome a tree node and add them to a list
    leaves = []
    for genome in genomes:
        leaves.append(Tree(genome.name()))
    # make the tree of genomes and print it
    tree = make_tree(leaves, sim)
    print(str(tree))


main()
