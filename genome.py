"""
    Author: Eduardo Hernandez
    File: genome.py
    Description: This program creates the Genome Data objects that hold the
    information of a genome. These objects have the name of the genome, their
    sequence, and the set of ngrams attributed to the genome.
"""


class GenomeData:
    """
    This class is used to hold the data of a genome. These objects simply store
    information and their information can be accessed with getter methods.
    """
    def __init__(self, name, sequence, ngrams):
        self._name = name
        self._sequence = sequence
        self._ngrams = ngrams

    def name(self):
        """
        This method just returns the name of the genome
        return: string representing the name of the genome
        """
        return self._name

    def sequence(self):
        """
        This method returns the sequence of the genome.
        return: string representing the sequence of the genome
        """
        return self._sequence

    def ngrams(self):
        """
        This method returns the set of ngrams of the genome
        return: set containing the ngrams of the genome
        """
        return self._ngrams

    def __str__(self):
        return str(self._name) + ": " + str(self._sequence) + \
               "; " + str(self._ngrams)

