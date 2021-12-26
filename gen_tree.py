from numpy import hstack, argmin, shape, unravel_index, delete, vstack, identity
import numpy as np
import argparse
from distancetree import *
from distance_utils.distance_matrix import Fasta2DistancesMatrix

map2    = lambda fn, mat: map(lambda arr: map(fn, arr), mat)
idxmin  = lambda mat: unravel_index(argmin(mat), shape(mat))
exists  = lambda x: x != None
wraparr = lambda x: [x]

# Remove rows and columns of matrix with the listed indices
withoutIndices = lambda m, ids: delete(delete(m, ids, axis=0), ids, axis=1)
# Append a vector as both a row and a column
appendRowCol   = lambda m, v: hstack((vstack((m, [v])), map(wraparr, v + [0])))


def get_saito_nei_matrix(D):
    n = D.shape()[0]
    r = np.sum(D, 1) / (shape(D)[0] - 2)
    Q = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            Q[i, j] = r[i] + r[j] - D[i,j]
    return Q


def getNeighbors(D, use_saito_nei):
    # returns the indecies of the neighbors to join.
    if use_saito_nei:
        D = get_saito_nei_matrix(D)
    return idxmin(D)


def neighborJoin(D, forest, use_saito_nei, transition_matrix):
    if len(D) == 2:
        return Tree(forest[0], forest[1])
    i, j = getNeighbors(D, use_saito_nei)
    u = [(D[i,k] + D[j,k] - D[i,j]) / 2 for k in range(len(D))]
    forest = hstack((forest, [Tree(forest[i], forest[j], transition_matrix=transition_matrix)]))
    D = appendRowCol(D, u)
    D = withoutIndices(D, (i, j))
    forest = delete(forest, (i, j))
    return neighborJoin(D, forest)


def get_argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--seq-file', help='The path to the sequences', required=True)
    parser.add_argument('--saito', action="store_true", help='if set, will use the saito-nei distance')
    return parser


def main(args):
    D, names, transition_matrix = Fasta2DistancesMatrix().distance_matrix_gen(args.seq_file)
    forest = [Leaf(seq) for seq in names]
    tree = neighborJoin(D, forest, True, transition_matrix)
    print(tree)


if __name__ == "__main__":
    parser = get_argparser()
    args = parser.parse_args()
    main(args)