# -*- coding: utf8 -*-
from __future__ import print_function
import argparse
import os
import sys
from bl.helper import Helper


class ExternalFileManager(object):

    STRING_KERNEL = 'my-string-kernel'

    @property
    def amino_file(self):
        return self._amino_file

    @property
    def train_id_file(self):
        return self._train_id_file

    @property
    def train_kernel_input(self):
        return self._train_kernel_input

    @property
    def test_id_file(self):
        return self._test_id_file

    @test_id_file.setter
    def test_id_file(self, value):
        self._test_id_file = value

    @property
    def globals_file(self):
        return self._globals_file

    @property
    def my_string_kernel(self):
        return self._my_string_kernel

    @property
    def blast_db(self):
        return self._blast_db

    @property
    def lookup_fasta(self):
        return self._lookup_fasta

    @property
    def psiblast2hssp(self):
        return self._psiblast2hssp

    @property
    def blast_prediction_printer(self):
        return self._blast_prediction_printer

    @property
    def files_to_remove(self):
        return self._remove_list

    @property
    def train_fasta_file(self):
        return self._train_fasta_file

    @property
    def best_params(self):
        return self._best_params

    def matrix_file_for_params(self, k_mer, sub_score):
        matrix_path = os.path.join(self.matrix_folder, 'l{k}_y{sub}.matrix'.format(k=k_mer, sub=sub_score))
        return matrix_path

    def normalized_matrix_file_for_params(self, k_mer, sub_score):
        matrix_path = os.path.join(self.matrix_folder, 'l{k}_y{sub}.norm.matrix'.format(k=k_mer, sub=sub_score))
        return matrix_path

    @property
    def matrix_folder(self):
        return self._matrix_folder

    def add_file_to_delete(self, file_to_remove):
        self._remove_list.append(file_to_remove)

    def add_file_list_to_deletion(self, files_to_remove):
        self._remove_list.extend(files_to_remove)

    def is_predictor_setup_sane(self, predict_traveller=False):
        if self.verbose:
            print('Running selfcheck - preparing to start predictions')

        script_folder = os.path.dirname(os.path.abspath(__file__))
        data_folder = os.path.join(script_folder, 'data')

        if predict_traveller:
            target_class_abbreviation = 'tr'
        else:
            target_class_abbreviation = 'sn'

        target_folder = os.path.join(data_folder, target_class_abbreviation)
        self._matrix_folder = os.path.join(target_folder, 'matrices')
        blast_script_folder = os.path.join(script_folder, 'blast_scripts')
        helper = Helper(self.verbose)
        # Check for needed files
        # 1) check for matrices - parameters taken from best results from parameter-optimization # TODO use getter
        if predict_traveller:
            helper.file_check(os.path.join(self._matrix_folder, 'l3_y6.matrix'))
        else:
            helper.file_check(os.path.join(self._matrix_folder, 'l3_y5.matrix'))
            helper.file_check(os.path.join(self._matrix_folder, 'l3_y7.matrix'))
            helper.file_check(os.path.join(self._matrix_folder, 'l4_y6.matrix'))
            helper.file_check(os.path.join(self._matrix_folder, 'l4_y7.matrix'))
            helper.file_check(os.path.join(self._matrix_folder, 'l4_y8.matrix'))
            helper.file_check(os.path.join(self._matrix_folder, 'l4_y9.matrix'))
            helper.file_check(os.path.join(self._matrix_folder, 'l5_y8.matrix'))
        # 2) check train-files:
        self._train_fasta_file = \
            os.path.join(self._matrix_folder, '{abr}_train.fasta'.format(abr=target_class_abbreviation))
        helper.file_check(self.train_fasta_file)
        self._train_id_file = \
            os.path.join(self._matrix_folder, '{abr}_train.idList'.format(abr=target_class_abbreviation))
        helper.file_check(self.train_id_file)
        self._train_kernel_input = \
            os.path.join(self._matrix_folder, '{abr}_train.psiBlastMat'.format(abr=target_class_abbreviation))
        helper.file_check(self.train_kernel_input)
        self._globals_file = os.path.join(self._matrix_folder, '{abr}.globals'.format(abr=target_class_abbreviation))
        helper.file_check(self.globals_file)
        self._best_params = os.path.join(target_folder, '{abr}_best_params'.format(abr=target_class_abbreviation))
        helper.file_check(self.best_params)
        # 3) check for kernel-creation-script
        self._my_string_kernel = os.path.join(data_folder, self.STRING_KERNEL)
        helper.file_check(self.my_string_kernel)
        # 4) check for BLAST-DB
        self._blast_db = os.path.join(target_folder, '{abr}_blastdb'.format(abr=target_class_abbreviation))
        helper.file_check(os.path.join(target_folder, '{abr}_blastdb.phr'.format(abr=target_class_abbreviation)))
        helper.file_check(os.path.join(target_folder, '{abr}_blastdb.pin'.format(abr=target_class_abbreviation)))
        helper.file_check(os.path.join(target_folder, '{abr}_blastdb.psq'.format(abr=target_class_abbreviation)))
        # 4b) check for mandatory blast-scripts
        self._psiblast2hssp = os.path.join(blast_script_folder, 'psi-blast2hssp.pl')
        self._blast_prediction_printer = os.path.join(blast_script_folder, 'PrintBlastPredictions.jar')
        # 5) check for Blast-Background-Fasta
        self._lookup_fasta = os.path.join(target_folder, '{abr}_lookup.fa'.format(abr=target_class_abbreviation))
        helper.file_check(self.lookup_fasta)
        # 6) check for Amino-File
        helper.file_check(self.amino_file)
        return True  # if there's an error before, we'd have exited already

    def __init__(self, verbose):
        self.verbose = verbose
        self._lookup_fasta = None
        self._train_id_file = None
        self._train_kernel_input = None
        self._train_fasta_file = None
        self._best_params = None
        self._test_id_file = None  # dynamically created while running
        self._test_kernel_input = None  # dynamically created while running
        self._globals_file = None
        self._matrix_folder = None
        self._my_string_kernel = None
        self._blast_db = None
        self._psiblast2hssp = None
        self._blast_prediction_printer = None
        self._amino_file = '/usr/share/fastprofkernel/data/Amino.txt'

        self._remove_list = list()


def main():
    """
    When accessing this script from the commandline a self-check will be started.
    :return:
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--verbose', help='Toggles verbose mode on', action='store_true')
    args = parser.parse_args()

    efm = ExternalFileManager(args.verbose)
    if efm.is_predictor_setup_sane():
        exit(0)
    else:
        exit(404)


def error(*objs):
    print("ERROR: ", *objs, file=sys.stderr)


if __name__ == "__main__":
     main()