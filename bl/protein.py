# -*- coding: utf8 -*-
from __future__ import print_function
import argparse
import sys


class Protein(object):

    """ A protein in the context of LocNuclei-Predictions """

    def __init__(self, verbose):
        self.verbose = verbose
        self.fasta_file = None
        self.blast_file = None
        self.location_prediction = None
        self.has_blast_hit = False
        self.has_prediction = False
        self.reliability = None

def error(*objs):
    print("ERROR: ", *objs, file=sys.stderr)
