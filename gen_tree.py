from numpy import hstack, argmin, shape, unravel_index, delete, vstack, identity
import numpy as np
import argparse
import pickle
from distancetree import *
from distance_utils.distance_matrix import Fasta2DistancesMatrix
from distance_utils.distance_matrix import DistanceCalculator
import matplotlib.pyplot as plt
USE_CACHE = False

map2    = lambda fn, mat: map(lambda arr: map(fn, arr), mat)
idxmin  = lambda mat: unravel_index(argmin(mat), shape(mat))
exists  = lambda x: x != None
wraparr = lambda x: [x]

# Remove rows and columns of matrix with the listed indices
withoutIndices = lambda m, ids: delete(delete(m, ids, axis=0), ids, axis=1)
# Append a vector as both a row and a column

def appendRowCol(m, v):
    m = vstack((m, np.atleast_2d(v)))
    v = np.append(v, 100000000)
    return hstack((m, np.atleast_2d(v).transpose()))


def get_saito_nei_matrix(D):
    n = D.shape[0]
    r = np.sum(D, 1) / (shape(D)[0] - 2)
    Q = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i == j:
                Q[i, j] = 1000000000
            else:
                Q[i, j] = r[i] + r[j] - D[i,j]
    return Q


def getNeighbors(D, use_saito_nei):
    # returns the indecies of the neighbors to join.
    if use_saito_nei:
        D = get_saito_nei_matrix(D)
    return idxmin(D)


def neighborJoin(D, forest, use_saito_nei, transition_matrix):
    while True:
        if len(D) == 2:
            return Tree(forest[0], forest[1], transition_matrix=transition_matrix)
        i, j = getNeighbors(D, use_saito_nei)
        u = np.array([(D[i,k] + D[j,k] - D[i,j]) / 2 for k in range(len(D))])
        forest = hstack((forest, [Tree(forest[i], forest[j], transition_matrix=transition_matrix)]))
        D = appendRowCol(D, u)
        D = withoutIndices(D, (i, j))
        forest = delete(forest, (i, j))


def get_argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--seq-file', help='The path to the sequences', required=True)
    parser.add_argument('--saito', action="store_true", help='if set, will use the saito-nei distance')
    return parser


def histogram(origin_path, others):
    f = open(origin_path, 'r')
    origin = ''.join([s.split('\n')[0] for s in f.readlines()])

    dc = DistanceCalculator()
    ls = []
    for seq in others:
        ls.append(dc.find_distance(origin, seq))
    plt.hist(ls, bins=30)
    plt.show()


def print_aligned(origin_path, s2):
    f = open(origin_path, 'r')
    origin = ''.join([s.split('\n')[0] for s in f.readlines()])
    i = 0
    while i + 50 < len(origin):
        print(origin[i:i + 50])
        print(s2[i:i + 50])
        print()
        i += 50
    print(origin[i:])
    print(s2[i:])
    print()

    dc = DistanceCalculator()
    score = dc.find_distance(origin, s2)
    print("===== Score:" + str(score) + " =====")
    print("score of origin vs origin:" + str(dc.find_distance(origin, origin)))



def main(args):
    calculator = Fasta2DistancesMatrix()
    if not USE_CACHE:
        D, seqs, names = calculator.distance_matrix_gen(args.seq_file)
    else:
        with open('./cache/distance.pickle', 'rb') as f:
            D = pickle.load(f)
        with open('./cache/names.pickle', 'rb') as f:
            names = pickle.load(f)
    transition_matrix = calculator.distance_calculator.dist_matrix
    D = D.max() - D
    print(names)
    # print(seqs)
    forest = [Leaf(seq, name) for name, seq in zip(names, seqs)]
    tree = neighborJoin(D, forest, args.saito, transition_matrix)
    tree.draw('saito' if args.saito else 'regular')


if __name__ == "__main__":
    parser = get_argparser()
    args = parser.parse_args()
    main(args)