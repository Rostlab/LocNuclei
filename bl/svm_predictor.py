# -*- coding: utf8 -*-
from __future__ import print_function
from math import sqrt
import os
from sklearn import svm
import sys
import subprocess
import numpy
from bl.helper import Helper
from collections import defaultdict

class SVMPredictor(object):

    QUERY_IDS_FILENAME = 'query_ids.lst'
    QUERY_KERNEL_INPUT_FILE = 'query_kernel_input.psiBlastMat'

    def __init__(self, is_verbose, working_directory, file_manager):
        self.verbose = is_verbose
        self.all_lookup_proteins = dict()
        self.working_directory = working_directory
        self.fm = file_manager

    def predict_all_query_proteins_without_blast_hit(self, all_query_proteins, ri):
        if self.verbose:
            print('Starting to create SVM-Predictions for all query-proteins:')

        # 1) Iterate over all left query proteins, i.e. the ones without predictions from blast
        # 1a) Create test-id-file
        test_id_file, test_id_list = self.create_test_id_file_from_protein_list(all_query_proteins)
        # 1b) Create test-kernel-input-file (sequences & profiles)
        test_kernel_input = self.create_test_kernel_input_file_from_protein_list(all_query_proteins, test_id_list)

        # 2) get best params:
        helper = Helper(self.verbose)
        all_params = helper.read_param_file(self.fm.best_params)
        # 2) Iterate over all parameter-combinations / classes
        for class_name in all_params:
            class_params = all_params[class_name]

            max_kmer_length = class_params['l']  # subprocess arguments need to be strings
            max_sub_score = class_params['y']
            # print(max_kmer_length)
            # print(max_sub_score)
            # 3) call my-string-kernel with the created files (for all parameters)
            normalized_query_matrix = self.call_string_kernel(max_kmer_length, max_sub_score, class_name)
            # print(numpy.matrix(normalized_query_matrix))
            # 4) predict
            normalized_train_matrix = self.fm.normalized_matrix_file_for_params(max_kmer_length, max_sub_score)
            query_acs = numpy.asarray(test_id_list)
            class_results = self.predict_query_matrix(normalized_train_matrix, normalized_query_matrix,
                                                      query_acs, class_name, class_params, ri)

            for result in class_results:
                if class_results[result][0]:
                    print(class_results)
                    if ri == True:
                        new_loc_prediction = class_results[result][0].replace('_', ' ')
                        r_index = str(class_results[result][1])                        
                        ri_prediction = all_query_proteins[result].reliability
                        print(r_index)
                    else:
                        new_loc_prediction = class_results[result].replace('_', ' ')
                    loc_prediction = all_query_proteins[result].location_prediction
                    if not loc_prediction or loc_prediction.isspace():
                        loc_prediction = new_loc_prediction
                        loc_prediction += '.'
                    else:
                        loc_prediction += ' '
                        loc_prediction += new_loc_prediction
                        loc_prediction += '.'
                        loc_prediction.strip()
                    if ri == True:
                        if not ri_prediction or ri_prediction.isspace():
                            ri_prediction = r_index
                            ri_prediction += '.'
                        else:
                            ri_prediction += ' '
                            ri_prediction += r_index
                            ri_prediction += '.'
                            ri_prediction.strip()
                    all_query_proteins[result].location_prediction = loc_prediction
                    all_query_proteins[result].has_prediction = True
                    all_query_proteins[result].reliability = ri_prediction

        return all_query_proteins

    def create_test_id_file_from_protein_list(self, all_query_proteins):
        if self.verbose:
            print('Creating file of all query-ids')
        test_id_list = list()
        for protein_id, protein in all_query_proteins.items():
            if protein.has_blast_hit:
                continue  # skip already predicted proteins
            else:
                test_id_list.append(protein_id)

        helper = Helper(self.verbose)

        query_id_file_path = None

        if self.working_directory and helper.folder_existence_check(self.working_directory):
            query_id_file_path = os.path.join(self.working_directory, self.QUERY_IDS_FILENAME)
            helper.file_not_there_check(query_id_file_path)  # quit if file already exists
            with open(query_id_file_path, 'w') as id_lst_file:
                for name in test_id_list:
                    id_lst_file.write('{protein_name}\n'.format(protein_name=name))

            self.fm.test_id_file = query_id_file_path
            self.fm.add_file_to_delete(query_id_file_path)
        else:
            error('Temporary working-directory {dr} was '
                  'not found or is not accessible'.format(dr=self.working_directory))
            exit(500)

        return query_id_file_path, test_id_list

    def create_test_kernel_input_file_from_protein_list(self, all_query_proteins, test_ids):
        if self.verbose:
            print('Creating file of all query-proteins used as input for string-kernel')
        helper = Helper(self.verbose)
        query_kernel_input_path = None
        if self.working_directory and helper.folder_existence_check(self.working_directory):
            query_kernel_input_path = os.path.join(self.working_directory, self.QUERY_KERNEL_INPUT_FILE)
            helper.file_not_there_check(query_kernel_input_path)  # quit if file already exists
            with open(query_kernel_input_path, 'w') as kernel_input:
                self.fm.add_file_to_delete(query_kernel_input_path)
                self.fm.test_kernel_input = query_kernel_input_path
                for protein in test_ids:  # loop over id-list to make sure we have always the same order of IDs
                    fasta_file = all_query_proteins[protein].fasta_file
                    cleaned_fasta = helper.clean_fasta_input(fasta_file, protein)
                    for line in cleaned_fasta:
                        kernel_input.write(line)
                    kernel_input.write('\n')  # add linebreak at end of sequence

                    blast_file = all_query_proteins[protein].blast_file
                    with open(blast_file, 'r') as blast_src:
                        kernel_input.write(blast_src.read())
        else:
            error('Temporary working-directory {dr} was '
                  'not found or is not accessible'.format(dr=self.working_directory))
            exit(500)

        return query_kernel_input_path

    def call_string_kernel(self, max_kmer_length, max_sub_score, class_name):
        max_kmer_length = str(max_kmer_length)  # subprocess arguments need to be strings
        max_sub_score = str(max_sub_score)

        kernel_call = list()
        kernel_call.append(self.fm.my_string_kernel)
        kernel_call.append('-o')
        kernel_call.append(self.fm.test_id_file)
        kernel_call.append('-O')
        kernel_call.append(self.fm.train_id_file)
        kernel_call.append('-p')
        kernel_call.append(self.fm.test_kernel_input)
        kernel_call.append('-P')
        kernel_call.append(self.fm.train_kernel_input)
        kernel_call.append('-K')
        kernel_call.append('-L')
        kernel_call.append(max_kmer_length)
        kernel_call.append('-Y')
        kernel_call.append(max_sub_score)
        kernel_call.append('-i')
        kernel_call.append(self.fm.amino_file)
        kernel_call.append('-g')
        kernel_call.append(self.fm.globals_file)

        if self.verbose:
            print('Calling String-Kernel for class {cl}'.format(cl=class_name))
            call = ''
            for c in kernel_call:
                call += ' {ca}'.format(ca=c)
            print(call)

        #debug:
        import datetime
        print('Start: {dt}'.format(dt=datetime.datetime.utcnow()))

        sp = subprocess.Popen(kernel_call, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = sp.communicate()
        if out:
            print('String-Kernel finished for parameters k-mer-length {k} and sub-scores {l}'.format(
                k=max_kmer_length, l=max_sub_score))
        if err:
            if self.verbose:
                print(err)
        if sp.returncode:
            code = int(sp.returncode)
            print(out)
            print(err)
            if code != 0:
                error('String-kernel returned with exit-code {nr}'.format(nr=code))
                exit(666)

        print('Finish: {dt}'.format(dt=datetime.datetime.utcnow()))

        if self.verbose:
            print('Normalizing query-matrix:')
            print('\t Reading diagonal values of train-matrix')
        train_diagonal_values = self.__get_diagonal_values_from_matrix_files(max_kmer_length, max_sub_score)

        if self.verbose:
            print('\t Calcualting normalized values for queries')


        normalized_rows = list()
        if out:
            # print(out)
            row_index = len(train_diagonal_values)
            start = True
            lines = out.split(b'\n')
            for line in lines:
                # print(line)
                if not line:
                    continue
                if start:                    
                        if line.strip().startswith(b'Read in all data files.'):
                            continue
                        else:
                            # print(line)
                            header_values = line.replace(b'\n', b'').strip().split()
                            try:
                                row_count = int(header_values[0])
                                col_count = int(header_values[1])
                                start = False
                                continue
                            except ValueError:
                                error('Was not able to parse row/column-counts from my-string-kernel output')
                                exit(712)
                else:
                    current_row_values = line.replace(b'\n', b'').strip().split()  # TODO fix to use correct diagonal values.
                    # last column in my-string-kernel output is diagonal value/self-hit
                    normalized_row = list()
                    column_index = 0
                    for cell_value in current_row_values:
                        if row_index < (row_count+len(train_diagonal_values)):
                            if column_index < col_count:  # the my-string-kernel output has +1 columsn for the diagonal values
                                self_hit_value = float(current_row_values[col_count])
                                cell_value = float(cell_value)
                                diag_product = (sqrt(float(train_diagonal_values[column_index])) * sqrt(self_hit_value))
                                normalized_value = cell_value / diag_product
                                normalized_row.append(float(normalized_value))
                                column_index += 1
                            
                            # $tmp_ar[$i]/(sqrt($diags[$i]*$diags[$ctr]));
                            # Spalte[$i] / sqrt(diagonalwert[$i] * diagonalwert[Zeilennr])
                            # Example:
                            # AC / (sqrt(CC * AA))

                    normalized_rows.append(normalized_row)
                    row_index += 1

            if self.verbose:
                print('Normalized query-matrix:')
                for r in normalized_rows:
                    print(r)

        # call to kernel:
        # my $kernelCommand = "$stringKernelExePath -o $combinedIDFile -O $trainIDsFilePath
        # -p $kernelInputFile -P $trainInputFilePath -K -L $maxKmerLength -Y $maxSubScore -i $aminoFilePath
        # -g $globalsFilePath  1> $kernelMatrixFilePath 2> $kernelMatrixErrorFilePath";

        return normalized_rows

    def __get_diagonal_values_from_matrix_files(self, max_kmer_length, max_sub_score):
        if self.verbose:
            print('Reading diagonal values for matrix with K={k} and L={l}'.format(k=max_kmer_length, l=max_sub_score))
        matrix_file = self.fm.matrix_file_for_params(max_kmer_length, max_sub_score)
        header = True
        row_count = 0
        row_counter = 0  # used for looping over the matrix

        diags = dict()

        with open(matrix_file, 'r') as matrix_src:
            for m_row in matrix_src:

                if header:
                    header_values = m_row.replace('\n', '').strip().split()
                    try:
                        row_count = int(header_values[0])
                        col_count = int(header_values[1])
                        if row_count == col_count:
                            header = False
                        else:
                            error('Precalculated Matrix-File {fl} was not in square-format, '
                                  'or the header is broken'.format(fl=matrix_file))
                            exit(701)
                    except ValueError:
                        error('Was not able to parse row/column-counts from matrix-file {fl}'.format(fl=matrix_file))
                        exit(702)
                else:
                    if m_row and row_counter < row_count:  # otherwise it runs over the max index (starts at 0)
                        diag = m_row.split()[row_counter]  # automatically gets the diagonal - 0/0, 1/1, 2/2...
                        diag_value = -1
                        try:
                            diag_value = int(diag)
                        except ValueError:
                            error('Was not able to parse diagonal value from matrix-file {fl}'
                                  ' \n Diagonal position {d}/{d}'.format(fl=matrix_file, d=row_counter))
                            exit(703)

                        diags[row_counter] = diag_value

                    row_counter += 1

        return diags

    def predict_query_matrix(self, train_matrix, query_matrix, query_acs, class_name, class_params, ri):
        helper = Helper(self.verbose)

        # 1) Read fasta file:
        train_fasta_file = self.fm.train_fasta_file
        y_values_train, ac_list_train = helper.read_fasta_file(class_name, train_fasta_file)

        # 2) Read matrix into array
        gram_train = helper.read_matrix_file(train_matrix)

        # 3) Train predictor
        c = class_params['C']
        tol = class_params['tol']
        class_weight = class_params['cw']

        if class_weight == 'auto':
            class_weight_auto = True
        else:
            class_weight_auto = False

        classifier = self.train_predictor(gram_train, y_values_train, c, tol, class_weight_auto)
        if self.verbose:
            print(classifier)

        # 4) Create numpy array from query-matrix
        gram_test = numpy.asarray(query_matrix)

        # 6) Predict
        results_for_class = None
        if ri == True:
            results_for_class = self.predict_with_reliability_index(classifier, gram_test, query_acs, class_name)
        else:
            results_for_class = self.predict_with_given_classifier(classifier, gram_test, query_acs, class_name)
        return results_for_class

    def train_predictor(self, gram_train, y_train, c, tol, cw_auto):
        if cw_auto:
            class_weights = 'balanced'
        else:
            class_weights = None

        c_param = float(c)
        tolerance = float(tol)

        classifier = svm.SVC(kernel='precomputed', probability=True, verbose=self.verbose,
                             C=c_param, tol=tolerance, class_weight=class_weights)
        cl = classifier.fit(gram_train, y_train)

        return cl

    def predict_with_given_classifier(self, classifier, gram_query, query_acs, class_name):
        if self.verbose:
            print('Predicting for {cl} for query proteins'.format(cl=class_name))
        # print(numpy.matrix(gram_query))
        predictions = classifier.decision_function(gram_query)

        if self.verbose:
            print('Got {nr} predictions'.format(nr=len(predictions)))
        results = dict()
        # print(predictions)

        for idx, prediction in enumerate(predictions):
            # print(prediction)
            prediction_class = None
            # if len(predictions) > 1 and prediction[0] > 0.0:
                # prediction_class = class_name
            if prediction > 0.0:
                prediction_class = class_name
            if self.verbose:
                print('Query AC {ac} is predicted with {p} for {c}'.format(
                    ac=query_acs[idx], p=prediction, c=class_name))

            results[query_acs[idx]] = prediction_class

        return results


    def predict_with_reliability_index(self, classifier, gram_query, query_acs, class_name):
        if self.verbose:
            print('Predicting for {cl} for query proteins'.format(cl=class_name))
        predictions = classifier.decision_function(gram_query)
        proba_predictions = classifier.predict_proba(gram_query)

        if self.verbose:
            print('Got {nr} predictions'.format(nr=len(predictions)))
        results = defaultdict(list)
       
        # predictions = numpy.array(predictions)
        # pos_predictions = predictions[:1] 
        # print(predictions)
        # print(pos_predictions)

        for idx, p in enumerate(predictions):
            prediction_class = None
            # if len(predictions) > 1 and prediction[0] > 0.0:
                # prediction_class = class_name
            # prediction = p[1]
            # print(prediction)
            if p > 0.0:
                prediction_class = class_name
            if self.verbose:
                print('Query AC {ac} is predicted with {p} for {c}'.format(
                    ac=query_acs[idx], p=p, c=class_name))

            # get reliability index ranging between 0 and 100
            prediction = proba_predictions[idx,1]
            # print(prediction)
            prediction = prediction * 100
            prediction = round(prediction,0)
            if prediction_class == None:
                prediction = None
            results[query_acs[idx]].append(prediction_class)
            results[query_acs[idx]].append(prediction)

        return results
        

def error(*objs):
    print("ERROR: ", *objs, file=sys.stderr)
