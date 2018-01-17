# -*- coding: utf8 -*-
""" Taking dictionary of query proteins as input and writing prediction-results to given file
"""
from __future__ import print_function
import sys


class ResultWriter(object):

    HEADER1 = '# Sub-nuclear Localization Prediction using LocNuclei\n' \
             '#\n' \
             '# NOTATION Protein Id: Fasta sequences Id truncated by whitespace\n' \
             '# NOTATION Localization: Predicted sub-nuclear localization class\n' \
             '# NOTATION Source: s == svm, b == blast\n' \
             '#\n' \
             '# Protein Id\tLocalization\tSource\n'
    HEADER2 = '# Sub-nuclear Localization Prediction using LocNuclei\n' \
             '#\n' \
             '# NOTATION Protein Id: Fasta sequences Id truncated by whitespace\n' \
             '# NOTATION Localization: Predicted sub-nuclear localization class\n' \
             '# NOTATION Source: s == svm, b == blast\n' \
             '# NOTATION Reliability Index (RI) between 0 and 100\n' \
             '#\n' \
             '# Protein Id\tLocalization\tSource\tRI\n'


    def write_results_to_file(self, proteins, out_file,ri):
        with open(out_file, 'w') as target:
            if ri == True:
                target.write(self.HEADER2)
            else:
                target.write(self.HEADER)
            for protein_name in proteins:
                ac = protein_name
                protein = proteins[ac]
                if protein.has_prediction:
                    loc = protein.location_prediction
                    if protein.has_blast_hit:
                        source = 'b'
                    else:
                        source = 's'
                    if ri == True:
                        reliability = protein.reliability
                else:
                    loc = 'unknown'
                    source = 'NA'
                    reliability = 'NA'

                result_line = ""
                if ri == True:
                    result_line = '{ac} \t {loc} \t {src} \t {r}\n'.format(ac=ac, loc=loc, src=source, r=reliability)
                else:
                    result_line = '{ac} \t {loc} \t {src}\n'.format(ac=ac, loc=loc, src=source)
                target.write(result_line)

    def __init__(self, is_verbose):
        self.verbose = is_verbose


def error(*objs):
    print("ERROR: ", *objs, file=sys.stderr)
