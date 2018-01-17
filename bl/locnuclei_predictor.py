# -*- coding: utf8 -*-
from __future__ import print_function
import glob
import os
import sys
import tempfile
import subprocess
from bl.blast_predictor import BlastPredictor
from bl.external_file_manager import ExternalFileManager
from bl.result_writer import ResultWriter
from bl.svm_predictor import SVMPredictor

from bl.helper import Helper
from bl.protein import Protein


class LocNucleiPredictor(object):

    def __init__(self, verbose, debug, predict_traveller):
        self.verbose = verbose
        self.debug = debug
        self.predict_traveller = predict_traveller

    def __enter__(self):
        # using encapsulated class in 'PackageResource' as in
        # http://stackoverflow.com/questions/865115/how-do-i-correctly-clean-up-a-python-object to force users to use
        # a with-statement. That way we can ensure a nice 'clean-up' of used resources (temp folders etc.)
        class LocNuclei(object):
            """
            The initial starting point for all locnuclei related methods and sub-procedures
            """
            def get_fasta_files(self, fasta_folder, fasta_suffix):
                """ Read all Fasta-Files in a big dictionary.
                :param fasta_folder: Folder where the fasta-files to predict are stored
                :param fasta_suffix: Suffix which defines the Fasta files
                :return:
                """

                fasta_pattern = os.path.join(fasta_folder, fasta_suffix)  # a pattern like '/mnt/home/foo/*.fasta'
                fasta_counter = 0
                for fasta_file in glob.iglob(fasta_pattern):
                    fasta_counter += 1
                    fasta_name = os.path.basename(fasta_file)
                    suffix_start = fasta_name.find('.')  # split file-name at the first 'dot'
                    fasta_name = fasta_name[:suffix_start]  # reduces the file-name to the prefix only
                    self.all_query_proteins[fasta_name] = Protein(self.verbose)
                    self.all_query_proteins[fasta_name].fasta_file = fasta_file

                if len(self.all_query_proteins) == fasta_counter:  # sanity check
                    if self.verbose:
                        print('Read {nr} FastaFiles from {dire}'.format(nr=fasta_counter, dire=fasta_folder))
                else:
                    error('Read {nr} of FastaFiles but only {nr2} unique prefixes'.format(
                        nr=fasta_counter, nr2=len(self.all_query_proteins)))
                    exit(500)

            def get_blast_files(self, blast_folder, blast_suffix):
                """ Read all Fasta-Files in a big dictionary.
                :param blast_folder: Folder where the fasta-files to predict are stored
                :param blast_suffix: Suffix which defines the Fasta files
                :return:
                """
                blast_pattern = os.path.join(blast_folder, blast_suffix)  # a pattern like '/mnt/home/foo/bar/*.blast'
                for blast_file in glob.iglob(blast_pattern):
                    blast_name = os.path.basename(blast_file)
                    suffix_start = blast_name.find('.')
                    blast_name = blast_name[:suffix_start]  # reduces the file-name to the prefix only
                    if blast_name in self.all_query_proteins:
                        self.all_query_proteins[blast_name].blast_file = blast_file
                    else:
                        error('Found Blast-File {bf} for which no fasta-file was provided'.format(bf=blast_file))
                        exit(404)

            def __create_temp_dir_from_input(self, temp_folder):
                helper = Helper(self.verbose)
                if temp_folder and helper.folder_existence_check(temp_folder):
                    # work in given directory
                    self.working_directory = temp_folder
                    self.remove_working_dir = False  # folder was given by user, we don't delete it
                else:
                    # create temporary directory
                    self.working_directory = tempfile.mkdtemp()
                    self.remove_working_dir = True  # we created it, we've to clean it up

            def prepare_temporary_directory(self, temp_folder):
                self.__create_temp_dir_from_input(temp_folder)
                work_dir = self.working_directory
                
                if self.verbose:
                    print('Setting working-directory to {dr}'.format(dr=work_dir))

                return work_dir

            def write_results_to_output_file(self, result_file, ri):
                writer = ResultWriter(self.verbose)
                writer.write_results_to_file(self.all_query_proteins, result_file, ri)

            def clean_up(self):
                if self.verbose:
                    print('Cleaning up temp-files and -folders')
                for file_to_delete in self.files_to_remove:
                    if self.verbose:
                        print('Deleting {fl}'.format(fl=file_to_delete))
                    os.remove(file_to_delete)

                for file_to_delete in self.file_manager.files_to_remove:
                    if self.verbose:
                        print('Deleting {fl}'.format(fl=file_to_delete))
                    os.remove(file_to_delete)

                if self.remove_working_dir:
                    if self.verbose:
                        print('Deleting {fl}'.format(fl=self.working_directory))
                    os.removedirs(self.working_directory)

            def predict_given_files(self, fasta_folder, fasta_suffix, blast_folder, blast_suffix, tmp_folder, out_file, ri, only_blast):
                # 0) Read files and determine the workload
                self.get_fasta_files(fasta_folder, fasta_suffix)
                self.get_blast_files(blast_folder, blast_suffix)
                self.prepare_temporary_directory(tmp_folder)
                # self.prepare_query_files(tmp_folder)
                # 1) Ask Blast for homologues proteins - if we've a hit we don't need to run the whole SVM-process:
                blaster = BlastPredictor(self.verbose, self.working_directory, self.file_manager, self.predict_traveller)
                self.all_query_proteins = blaster.predict_all_query_proteins(self.all_query_proteins)

                # 2) Run SVMs for proteins without an blast-hit
                if only_blast == False:
                    profkernel_predictor = SVMPredictor(self.verbose, self.working_directory, self.file_manager)
                    self.all_query_proteins = profkernel_predictor.predict_all_query_proteins_without_blast_hit(self.all_query_proteins, ri)

                # 3) Write out results
                self.write_results_to_output_file(out_file, ri)

            def __init__(self, is_verbose, predict_traveller):
                self.verbose = is_verbose
                self.all_query_proteins = dict()

                self.remove_working_dir = False
                self.working_directory = None
                self.files_to_remove = list()

                self.predict_traveller = predict_traveller

                self.file_manager = ExternalFileManager(is_verbose)
                self.file_manager.is_predictor_setup_sane(predict_traveller)

        self.package_obj = LocNuclei(self.verbose, self.predict_traveller)
        return self.package_obj

    def __exit__(self, type, value, traceback):
        if not self.debug:
            self.package_obj.clean_up()


def error(*objs):
    print("ERROR: ", *objs, file=sys.stderr)
