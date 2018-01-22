# -*- coding: utf8 -*-
"""
DESCRIPTION:

Sub-nuclear localization prediction for proteins of all (eukaryotic) species
"""
from __future__ import print_function
import argparse
import sys

from bl.helper import Helper
from bl.locnuclei_predictor import LocNucleiPredictor


def main():

    usage_string = 'python locnuclei.py example/ example/ -t results.log'

    # parse command line options
    parser = argparse.ArgumentParser(description=__doc__, usage=usage_string)
    parser.add_argument('fasta_folder', help='Folder with protein sequences in fasta-format. '
                                             'Every file may only contain one sequence.')
    parser.add_argument('--fasta_suffix', help='Suffix of files in given Fasta-folder (default: "*.fasta")',
                        default='*.fasta')
    parser.add_argument('blast_folder', help='Folder with Blast-Profiles for the fasta files. '
                                             'Blast-Files need to have the same name in front of the suffix '
                                             'as the fasta files! (e.g. Q9XLZ3.fasta <-> Q9XLZ3.profile)')
    parser.add_argument('--blast_suffix', help='Suffix of files in given Blast-Folder (default: "*.profile")',
                        default='*.profile')
    parser.add_argument('--temp_folder', help='Folder to work in, will be automatically deleted afterwards. '
                                              'If not given, a temporary directory will be created.')
    parser.add_argument('output_file', help='Where should results be written to?')
    parser.add_argument('-t', '--traveller', help='Predict nuclear travelling proteins '
                                                  'instead of sub-nuclear localization', action='store_true')

    parser.add_argument('-d', '--debug',
                        help='Toggles clean up of temporary files off, '
                             'i.e. no files will be deleted, that were created during prediction',
                        action='store_true')
    parser.add_argument('-v', '--verbose', help='Toggles verbose mode on', action='store_true')
    parser.add_argument('-b', '--only_blast', help='Only run the blast search', action='store_true')
    args = parser.parse_args()
    print(args)
    helper = Helper(args.verbose)
    if helper.folder_existence_check(args.fasta_folder):  # check if fasta-folder exists and is reachable
        if helper.folder_existence_check(args.blast_folder):  # check if blast-profile-folder exists and is reachable
            if helper.file_not_there_check(args.output_file):  # check if output-file doesn't exist yet -
                                                               # no overwriting of existing files
                with LocNucleiPredictor(args.verbose, args.debug, args.traveller) as loc_nuclei:
                    loc_nuclei.predict_given_files(args.fasta_folder, args.fasta_suffix, args.blast_folder,
                                                   args.blast_suffix, args.temp_folder, args.output_file, 
                                                   args.only_blast)


def error(*objs):
    print("ERROR: ", *objs, file=sys.stderr)


if __name__ == "__main__":
    main()
