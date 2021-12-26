from numpy import hstack, argmin, shape, unravel_index, delete, vstack, identity
import numpy as np
from distancetree import *

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


def neighborJoin(D, forest, use_saito_nei):
    if len(D) == 2:
        return Tree(forest[0], forest[1])
    i, j = getNeighbors(D, use_saito_nei)
    u = [(D[i,k] + D[j,k] - D[i,j]) / 2 for k in range(len(D))]
    forest = hstack((forest, [Tree(forest[i], forest[j])]))
    D = appendRowCol(D, u)
    D = withoutIndices(D, (i, j))
    forest = delete(forest, (i, j))
    return neighborJoin(D, forest)


if __name__ == "__main__":
    tree = neighborJoin(D, forest)
    print(tree)
