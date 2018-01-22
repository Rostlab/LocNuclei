# -*- coding: utf8 -*-
from __future__ import print_function
import os
import sys
import subprocess
from bl.helper import Helper


class BlastPredictor(object):
    SUBNUCLEAR_EVALUE = -20
    TRAVELLER_EVALUE = -5
    BLAST_MIN = 20.0

    def predict_all_query_proteins(self, all_query_proteins):
        if self.verbose:
            print('Starting to create BLAST-Predictions for all query-proteins:')

        # 1) initialize needed information
        self.get_locations_from_lookup_fasta(self.fm.lookup_fasta)

        # 2) query blast-db with all proteins:
        blast_prediction_file = self.blast_query_protein_against_db(all_query_proteins)

        # 3) extract predictions from blast-prediction-file
        all_query_proteins = self.extract_predicitons_from_file(blast_prediction_file, all_query_proteins)

        return all_query_proteins

    def get_locations_from_lookup_fasta(self, lookup_fasta):
        """
        :param lookup_fasta: Open LookUp-Fasta which was used for to build the blastDB.
                Read all locations and corresponding ACs.
        :return: None
        """
        with open(lookup_fasta, 'r') as lookup:
            for line in lookup:
                if line.startswith('>'):  # header found!
                    header = line.split('#')
                    ac = header[1].strip()
                    locations = header[4].replace('\n', '').strip()

                    if self.predict_traveller:
                        if 'Traveller' in locations:
                            locations = 'Traveller.'
                        else:
                            locations = 'NOT Traveller.'

                    if ac in self.all_lookup_proteins:
                        error('Internal error. AC {u_ac} is multiple times in lookup-fasta!'.format(u_ac=ac))
                        exit(-1)
                    else:
                        self.all_lookup_proteins[ac] = locations

    def blast_query_protein_against_db(self, query_proteins):
        if self.predict_traveller:
            e_value = self.TRAVELLER_EVALUE
        else:
            e_value = self.SUBNUCLEAR_EVALUE

        # 1) calling:
        # /usr/bin/blastpgp -F F -a 1 -j 3 -b 150 -e 1e${e_param} -h 1e-10 -d $DB -i $(pwd)/${DIR}/$(basename $i)
        # -o ${folder}${prefix}.blastPsiOutTmp
        # based on runPsiBlastProfileCreatorAli.sh
        blast_files = list()
        hssp_files = list()
        for cur_query_protein in query_proteins:
            fasta = query_proteins[cur_query_protein].fasta_file
            if self.verbose:
                print('############')
                print('Blast Prediciton for {ac}:'.format(ac=cur_query_protein))
            blast_file = self.__call_blast(e_value, cur_query_protein, fasta)
            hssp_file = '{blast_file}.psiBlast2hssp'.format(blast_file=blast_file)
            blast_files.append(blast_file)
            hssp_files.append(hssp_file)

        self.fm.add_file_list_to_deletion(blast_files)
        self.fm.add_file_list_to_deletion(hssp_files)

        # 2) calling perl psi-blast2hssp.pl {blast_dir}       
        self.__call_hssp_calculations()  # works on full directory
        # 3) Create Blast-Predictions
        blast_prediction_file = self.__create_blast_prediction_file()  # works on full directory

        return blast_prediction_file

    def extract_predicitons_from_file(self, blast_prediction_file, all_query_proteins):
        with open(blast_prediction_file, 'r') as pred_src:
            for line in pred_src:
                if line and not line.isspace():
                    # example-line:
                    # A4Q9E5	A4Q9E5	100.0
                    predictions = line.split()
                    if len(predictions) > 1:
                        query = predictions[0].strip()
                        hit = predictions[1].strip()
                        hit_locations = self.all_lookup_proteins[hit]
                        seq_id = (float(predictions[2].strip())-self.BLAST_MIN)*100/(100-self.BLAST_MIN)
                        seq_id = round(seq_id,2)
                        if query in all_query_proteins:
                            all_query_proteins[query].has_blast_hit = True
                            all_query_proteins[query].has_prediction = True
                            all_query_proteins[query].reliability = seq_id
                            all_query_proteins[query].location_prediction = hit_locations.replace('\n', '').strip()
                        else:
                            error('We built a blast-query for {q} but it is not in the query-proteins - '
                                  'This should not be possible'.format(q=query))
                        if self.verbose:
                            print('Found BLAST-HIT: -> {q} {h} {hl}'.format(q=query, h=hit, hl=hit_locations))
        
        return all_query_proteins

    def __create_blast_prediction_file(self):
        predictor_call, prediction_file = self.__build_blast_prediction_call()
        # 3) write blast-predictions to result-file
        sp = subprocess.Popen(predictor_call, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = sp.communicate()
        if out:
            if self.verbose:
                print('JAR-Out: \n{o}'.format(o=out))
        if err:
            if self.verbose:
                print('JAR-Error: \n{e}'.format(e=err))
        if sp.returncode:
            code = int(sp.returncode)
            print('JAR-Out: \n{o}'.format(o=out))
            print('JAR-Error: \n{e}'.format(e=err))
            if code != 0:
                error('PrintBlastPredictions.jar returned with exit-code {nr}'.format(nr=code))
                exit(9)

        self.fm.add_file_to_delete(prediction_file)  # clean up!

        return prediction_file

    def __build_blast_prediction_call(self):
        java_predictor = self.fm.blast_prediction_printer
        helper = Helper(self.verbose)
        prediction_file = os.path.join(self.working_directory, 'blast.predictions')
        helper.file_not_there_check(prediction_file)
        predictor_call = ['java', '-jar', java_predictor,  self.working_directory, 'maxSeqId', prediction_file]

        return predictor_call, prediction_file

    def __call_hssp_calculations(self):
        hssp_call = self.__build_hssp_call()
        # 2) build hssp-files:
        sp = subprocess.Popen(hssp_call, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = sp.communicate()
        if out:
            if self.verbose:
                print('HSSP-Out: \n{o}'.format(o=out))
        if err:
            if self.verbose:
                print('HSSP-Error: \n{e}'.format(e=err))
            if 'Use of uninitialized value $query in concatenation' in str(err):
                # we are using the wrong separator character! (# instead of |)
                error('HSSP calculation failed, because query-ac could not be determined in input fasta!')
                exit(8)
        if sp.returncode:
            code = int(sp.returncode)
            print('HSSP-Out: \n{o}'.format(o=out))
            print('HSSP-Error: \n{e}'.format(e=err))
            if code != 0:
                error('psi-blast2hssp.pl returned with exit-code {nr}'.format(nr=code))
                exit(8)

    def __build_hssp_call(self):
        path_to_perlscript = self.fm.psiblast2hssp
        hssp_call = ['perl', path_to_perlscript, self.working_directory]
        if self.verbose:
            call_string = ''
            for c in hssp_call:
                call_string += ' {sub_call}'.format(sub_call=c)
        return hssp_call

    def __call_blast(self, evalue, protein_name, query_fasta):
        fasta_in = query_fasta
        blast_out = protein_name+'.blastPsiOutTmp'
        blast_out = os.path.join(self.working_directory, blast_out)

        blast_call = self.__build_blast_call(evalue, fasta_in, blast_out)

        sp = subprocess.Popen(blast_call, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = sp.communicate()
        if out:
            if self.verbose:
                print(out)
        if err:
            if self.verbose:
                print(err)
        if sp.returncode:
            code = int(sp.returncode)
            print(out)
            print(err)
            if code != 0:
                error('blastpgp returned with exit-code {nr}'.format(nr=code))
                exit(7)

        return blast_out

    def __build_blast_call(self, e_value, fasta_in, blast_out):
        e_value = '1e{e_val}'.format(e_val=e_value)
        blast_call = ['/usr/bin/blastpgp', '-F', 'F', '-a', '1', '-j', '3', '-b', '150', '-e', e_value, '-h', '1e-10',
                      '-d', self.fm.blast_db, '-i', fasta_in, '-o', blast_out]

        return blast_call

    def __init__(self, is_verbose, working_directory, file_manager, predict_traveller):
        self.verbose = is_verbose
        self.all_lookup_proteins = dict()
        self.working_directory = working_directory
        self.fm = file_manager
        self.predict_traveller = predict_traveller


def error(*objs):
    print("ERROR: ", *objs, file=sys.stderr)
