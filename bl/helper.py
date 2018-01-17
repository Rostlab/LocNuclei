# -*- coding: utf8 -*-
from __future__ import print_function
import os
import re
import sys
import datetime
import numpy


class Helper(object):
    def clean_fasta_input(self, fasta_file, proteinname):
        clean_fasta = list()
        with open(fasta_file, 'r') as fasta_src:
            for line in fasta_src:
                if line and not line.isspace():  # skip empty lines
                    if line.startswith('>'):  # header
                        line = '>{p_name}\n'.format(p_name=proteinname)
                    else:
                        line = line.replace('\n', '').strip()  # remove surplus linebreaks

                    clean_fasta.append(line)

        return clean_fasta

    def file_check(self, file_name):
        """
        Checks if file is available, if not exits the program!
        :param file_name:
        :return:
        """
        if os.path.isfile(file_name):
            file_found = True
        else:
            error('File {fl} not available - exit!'.format(fl=file_name))
            file_found = False
            exit(404)

        return file_found

    def file_not_there_check(self, file_name):
        """
        Checks if given file_name doesn't already exist, if it exists, quits the program!
        :param file_name:
        :return:
        """
        if not os.path.isfile(file_name):
            file_is_free = True
        else:
            error('File {fl} exists already - exit!'.format(fl=file_name))
            file_is_free = False
            exit(500)

        return file_is_free

    def folder_existence_check(self, folder_name=None):
        """
        Check if folder exists and is reachable, if not exit the program.
        :param folder_name: Folder to check, if it exists (default: None)
        :type folder_name: String'
        :return: True if folder exists and is reachable, False otherwise
        """
        if os.path.exists(folder_name):
            folder_is_there = True
        else:
            error('Folder {fl} does not exist or is not reachable - exit!'.format(fl=folder_name))
            folder_is_there = False
            exit(404)

        return folder_is_there

    def read_protein_numbers(self, proteinnumbers_file):
        """

        :param proteinnumbers_file:
        :return: tupel(train, test)
        """
        number_regex = re.compile(r"([0-9]+) Trainset([0-9]+) Testsettrain;.*")

        with open(proteinnumbers_file, 'r') as number_src:
            for line in number_src:
                if number_regex.match(line):
                    for train, test in number_regex.findall(line):
                        if self.verbose:
                            print(line)
                            print(train, test)
                        return int(train), int(test)

    def read_fasta_file(self, class_name, fasta_file):
        """

        :param class_name:
        :param fasta_file:
        :return: Array of y-values read from the fasta
        """
        ac_set = set()
        ac_list = list()
        y_list = list()
        with open(fasta_file, 'r') as fasta_src:
            for line in fasta_src:
                if line.startswith(">"):
                    cols = line.split('#')
                    ac_set.add(cols[1].strip())
                    ac_list.append(cols[1].strip())
                    locations = cols[4].strip().split('.')
                    locations = filter(None, locations)  # remove empty entries
                    locations = [location.strip() for location in locations]  # strip whitespace

                    if class_name.replace(' ', '_') in locations or class_name.replace('_', ' ') in locations:
                        # translate from class name to binary representation
                        y_list.append(1)
                    else:
                        y_list.append(0)

        # sanity check:
        if not len(ac_set) == len(y_list):
            error('AC set {anr} is not as long as y-list {ynr}'.format(anr=len(ac_set), ynr=len(y_list)))

        if self.verbose:
            print('For {cl} Read {pnr} positives and {nnr} negatives'.format(cl=class_name,
                                                                             pnr=y_list.count(1), nnr=y_list.count(0)))

        y_array = numpy.asarray(y_list)  # adapt to numpy-format
        return y_array, ac_list

    def read_fasta_file_for_prediction_only(self, class_name, fasta_file):
        """

        :param class_name:
        :param fasta_file:
        :return: Array of y-values read from the fasta
        """
        ac_list = list()
        y_list = list()
        prediction_set_hit = False
        with open(fasta_file, 'r') as fasta_src:
            for line in fasta_src:
                if line.startswith(">"):
                    cols = line.split('#')
                    if len(cols) < 3:
                        # no location given, we predict only
                        prediction_set_hit = True
                        ac_list.append(cols[1].replace('\n', '').strip())
                    else:
                        # location is given, we usually use it for training
                        if prediction_set_hit:
                            error('Found protein with locations after protein without location - is fasta mixed up?')
                            exit(500)
                        ac_list.append(cols[1].replace('\n', '').strip())
                        locations = cols[4].strip().split('.')
                        locations = filter(None, locations)  # remove empty entries
                        locations = [location.strip() for location in locations]  # strip whitespace

                        if class_name.replace(' ', '_') in locations or class_name.replace('_', ' ') in locations:
                            # translate from class name to binary representation
                            y_list.append(1)
                        else:
                            y_list.append(0)

        if self.verbose:
            print('For {cl} Read {pnr} positives and {nnr} negatives'.format(cl=class_name,
                                                                             pnr=y_list.count(1), nnr=y_list.count(0)))
        if self.verbose:
            print('Read {nr} ACs'.format(nr=len(ac_list)))

        y_array = numpy.asarray(y_list)  # adapt to numpy-format
        return y_array, ac_list

    def read_matrix_file(self, matrix_file, header_row_count=1):
        # gram = numpy.genfromtxt(matrix_file, skip_header=header_row_count)

        if not self.file_check(matrix_file):
            exit(404)

        gram_train = numpy.loadtxt(matrix_file, skiprows=header_row_count)

        if self.verbose:
            print("gram_train_shape: " + str(gram_train.shape))

        return gram_train

    def read_param_file(self, param_file):
        """ Reads best parameters for predictions from file """
        if self.verbose:
            print('Reading parameters from file {fl}'.format(fl=param_file))
        params_from_file = dict()
        with open(param_file, 'r') as param_src:
            for line in param_src:
                if line and not line.isspace():
                    if line.startswith('#'):
                        continue  # skip comments
                    elif line.replace('\n', '') == 'class;fold;f1;C;tol;l;y;class_weights':
                        continue  # skip header
                    else:
                        # looks like this: Chromatin;fold_5;0.68622294052;1.0;0.1;4;7;auto
                        splitted_line = line.split(';')
                        classname = splitted_line[0].strip()
                        leslie_l = int(splitted_line[5].strip())
                        leslie_y = int(splitted_line[6].strip())
                        c = float(splitted_line[3])
                        tol = float(splitted_line[4])
                        cw = splitted_line[7].replace('\n', '').strip()

                        params_from_file[classname] = dict()
                        params_from_file[classname]['l'] = leslie_l
                        params_from_file[classname]['y'] = leslie_y
                        params_from_file[classname]['C'] = c
                        params_from_file[classname]['tol'] = tol
                        params_from_file[classname]['cw'] = cw

        return params_from_file


    def print_result_dictionary_to_file(self, result_dict, class_name, out_file):
        true_positive = result_dict['TP']
        true_negative = result_dict['TN']
        false_negative = result_dict['FN']
        false_positive = result_dict['FP']
        condition_positive = result_dict['CP']
        condition_negative = result_dict['CN']

        # #2014-08-23 20:35
        # cl;tp;tn;fn;fp;cp;cn
        # Nuclear matrix;0;14;0;0;2;29
        with open(out_file, 'w') as out:
            now = datetime.datetime.now()
            out.write('#'+now.strftime("%Y-%m-%d %H:%M")+'\n')
            out.write('cl;tp;tn;fn;fp;cp;cn\n')
            result_string = '{cl};{tp};{tn};{fn};{fp};{cp};{cn}\n'.format(cl=class_name, tp=true_positive,
                                                                          tn=true_negative, fn=false_negative,
                                                                          fp=false_positive, cp=condition_positive,
                                                                          cn=condition_negative)
            out.write(result_string)

    def __init__(self, verbose):
        self.verbose = verbose


def error(*objs):
    print("ERROR: ", *objs, file=sys.stderr)
