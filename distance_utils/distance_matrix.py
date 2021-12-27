from pathlib import Path
from typing import Dict, List
from Bio import SeqIO
import pickle
from scipy.spatial import distance

import os
import Bio
import numpy as np


class Fasta2DistancesMatrix:
    def __init__(self):
        self.pickle_path = Path('./dist_matrix_as_pickle.pickle')
        self.distance_calculator = DistanceCalculator()

    def _distance_matrix_gen(self, filename='./sequences.fasta'):
        sequences_obj = SeqIO.parse(filename, 'fasta')
        seqs_dct = {sequence.id: sequence.seq for sequence in sequences_obj}
        distance_matrix = [[] for _ in range(len(seqs_dct))]
        for idx_row, (id_row, seq_row) in enumerate(seqs_dct.items()):
            for idx_col, (id_col, seq_col) in enumerate(seqs_dct.items()):
                current_distance = self.distance_calculator.find_distance(seq_row, seq_col)
                distance_matrix[idx_row].append(current_distance)
        np_array = np.array(distance_matrix)
        return np_array, seqs_dct.keys(), self.distance_calculator.dist_matrix

    def _matrix_saver(self, np_matrix, ids):
        with open(self.pickle_path, 'wb') as handle:
            pickle.dump({'dist_matrix': np_matrix, 'ids': list(ids)}, handle, protocol=pickle.HIGHEST_PROTOCOL)

    def _matrix_loader(self):
        with open(self.pickle_path, 'rb') as handle:
            content = pickle.load(handle)
        dist_matrix, ids = content['dist_matrix'], content['ids']
        return dist_matrix, ids

    def get_distance_matrix(self, load_if_exists=True):
        if load_if_exists and self.pickle_path.exists():
            dist_matrix, ids = self._matrix_loader()
            return dist_matrix, ids
        dist_matrix, ids = self._distance_matrix_gen()
        self._matrix_saver(dist_matrix, list(ids))
        return dist_matrix, ids


class MatrixLoader:
    def __init__(self, matrix_name='PAM10'):
        self.matrix_name = matrix_name
        dir_name = os.path.dirname(os.path.abspath(__file__))
        self.matrix_path = Path(os.path.join(dir_name, 'distance_matrices', f'distance_matrix_{self.matrix_name}'))
        self.matrix = self.__matrix_loader()

    def get_matrix(self):
        return self.matrix

    def __matrix_loader(self):
        with open(self.matrix_path) as file:
            content = file.read()
        lines = content.strip().split('\n')
        columns_names_0 = lines[0].split(' ')
        columns_names = self.empty_strings_cleaner(columns_names_0)
        entries_0 = [line.split(' ')[1:] for line in lines[1:]]
        entries = [self.empty_strings_cleaner(line) for line in entries_0]
        matrix = self.lst_2_dct_matrix_converter(columns_names, entries)
        return matrix

    def empty_strings_cleaner(self, lst):
        new_lst = list(filter(lambda name: name != '', lst))
        return new_lst

    def lst_2_dct_matrix_converter(self, columns_names, entries):
        matrix = {column_name: {} for column_name in columns_names}
        for row_idx, column_name_row in enumerate(columns_names):
            for col_idx, column_name_col in enumerate(columns_names):
                matrix[column_name_row][column_name_col] = float(entries[row_idx][col_idx])
        return matrix


class DistanceCalculator:
    def __init__(self, is_hamming=False, dist_matrix: Dict[str, Dict[str, int]] = MatrixLoader()):
        self.is_hamming = is_hamming
        self.dist_matrix = dist_matrix.get_matrix()

    def find_distance(self, str1, str2):
        if self.is_hamming:
            dist = self.hamming(str1, str2)
            return dist
        dist = self.matrix_distance(str1, str2)
        return dist

    def hamming(self, str1, str2):
        dist = sum([1 for char1, char2 in zip(str1, str2) if char1 != char2])
        return dist

    def matrix_distance(self, str1, str2):
        dist = sum([self.dist_matrix[char1][char2] for char1, char2 in zip(str1, str2)])
        return dist


#########################################
#                   TESTS               #


def basic_tests():
    dc = DistanceCalculator()
    d0 = dc.find_distance('NNN', 'NNN')
    d1 = dc.find_distance('ARN', 'NRN')
    d2 = dc.find_distance('ARN', 'NNN')
    print(f'd0 == {d0}')
    print(f'd1 == {d1}')
    print(f'd2 == {d2}')
    ############################################
    print('Now Hamming')
    dc_ham = DistanceCalculator(is_hamming=True)
    d0 = dc_ham.find_distance('NNN', 'NNN')
    d1 = dc_ham.find_distance('ARN', 'NRN')
    d2 = dc_ham.find_distance('ARN', 'NNN')
    print(f'd0 == {d0}')
    print(f'd1 == {d1}')
    print(f'd2 == {d2}')


def general_test():
    F2D = Fasta2DistancesMatrix()
    matrix = F2D.get_distance_matrix()
    return matrix


if __name__ == '__main__':
    basic_tests_ = False
    general_test_ = True
    if basic_tests_:
        basic_tests()
    if general_test_:
        dist_mat = general_test()
        print(dist_mat)
