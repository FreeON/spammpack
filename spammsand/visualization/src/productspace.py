#!/usr/bin/env python

from mayavi import mlab
import itertools
import math
import numpy
import random

class SpAMM:

    @staticmethod
    def index(i, j, k):
        """Return the linear octree index.
        """
        return i*4+j*2+k

    def __init__(self):
        self.children = [None for i in range(8)]

    def multiply(self, A, B, C, tolerance, symbolic_only=True, tier=0):
        """Multiply two quadtree nodes under the SpAMM condition.
        """

        if tier == 0:
            self.depth = A.depth
            self.N_padded = A.N_padded
            self.N_b = A.N_b

        if tier < A.depth:
            for (i, j, k) in itertools.product(range(2), range(2), range(2)):
                if A.children[Quadtree.index(i, k)] is not None:
                    A_child = A.children[Quadtree.index(i, k)]
                    if B.children[Quadtree.index(k, j)] is not None:
                        B_child = B.children[Quadtree.index(k, j)]
                        if (A_child.norm2*B_child.norm2) > tolerance**2:
                            if self.children[SpAMM.index(i, j, k)] is None:
                                self.children[SpAMM.index(i, j, k)] = SpAMM()
                            if C.children[Quadtree.index(i, j)] is None:
                                C.children[Quadtree.index(i, j)] = Quadtree(C.N, N_b=C.N_b,
                                                                            i_lower=A_child.i_lower,
                                                                            i_upper=A_child.i_upper,
                                                                            j_lower=B_child.j_lower,
                                                                            j_upper=B_child.j_upper,
                                                                            tier=tier+1,
                                                                            depth=C.depth)
                            C_child = C.children[Quadtree.index(i, j)]
                            self.children[SpAMM.index(i, j, k)].multiply(A_child, B_child, C_child,
                                                                         tolerance=tolerance,
                                                                         symbolic_only=symbolic_only,
                                                                         tier=tier+1)
        else:
            if not symbolic_only:
                if C.submatrix is None:
                    C.submatrix = numpy.matrix(numpy.zeros((C.N_b, C.N_b)))
                C.submatrix = C.submatrix+numpy.dot(A.submatrix, B.submatrix)

    def get_nodes(self, depth, tier=0,
                  i_lower=0, i_upper=0,
                  j_lower=0, j_upper=0,
                  k_lower=0, k_upper=0):
        """Get a list of nodes on tier == depth.

        The function retusn a tuple of ([i], [j], [k], width).  The
        ([i], [j], [k]) part is a list of the centers of the nodes, and
        the width is the index width on this tier.
        """

        if tier == 0:
            i_upper = self.N_padded
            j_upper = self.N_padded
            k_upper = self.N_padded

        half = (i_upper-i_lower)/2
        if tier < depth:
            x = []
            y = []
            z = []
            width = 2*half

            for (i, j, k) in itertools.product(range(2), range(2), range(2)):
                if self.children[SpAMM.index(i, j, k)] is not None:
                    (x_temp, y_temp, z_temp, width) = self.children[SpAMM.index(i, j, k)].get_nodes(
                        depth,
                        tier=tier+1,
                        i_lower=i_lower+i*half, i_upper=i_lower+(i+1)*half,
                        j_lower=j_lower+j*half, j_upper=j_lower+(j+1)*half,
                        k_lower=k_lower+k*half, k_upper=k_lower+(k+1)*half)
                    x = x+x_temp
                    y = y+y_temp
                    z = z+z_temp
            return (x, y, z, width)
        else:
            return ([i_lower+half], [j_lower+half], [k_lower+half], 2*half)

class Quadtree:

    @staticmethod
    def index(i, j):
        """Return the linear Quadtree index
        """
        return i*2+j

    @staticmethod
    def random(N, N_b=16, gamma=-1):
        """Generate a random matrix.
        """

        A = Quadtree(N, N_b=N_b, gamma=gamma)
        A_dense = numpy.matrix(numpy.zeros((A.N_padded, A.N_padded)))
        for i in range(N):
            for j in range(i, N):
                A_dense[i, j] = random.uniform(0.6, 1.0)*math.exp(-abs(i-j)/A.gamma)
                A_dense[j, i] = A_dense[i, j]
        for i, j in itertools.product(range(A.N_padded/N_b), range(A.N_padded/N_b)):
            A.set_submatrix(i*N_b, j*N_b, A_dense[i*N_b:(i+1)*N_b, j*N_b:(j+1)*N_b])
        return A

    @staticmethod
    def zeros(N, N_b=16):
        """Generate a zero matrix.
        """

        return Quadtree(N, N_b=N_b)

    def __init__(self, N, N_b=16,
                 N_padded=-1,
                 i_lower=-1, i_upper=-1,
                 j_lower=-1, j_upper=-1,
                 gamma=-1, tier=0, depth=-1):
        """Generate a matrix of size N by N with exponential decay.

        The decay constant is gamma.
        """

        edge = 12 # The minimum size of matrix elements, 10^{-edge}.

        self.N = N
        self.N_b = N_b
        if gamma < 0:
            self.gamma = N/(edge*math.log(10))
        else:
            self.gamma = gamma
        self.tier = tier
        self.children = [ None for i in range(4) ]
        self.submatrix = None
        self.norm2 = 0
        self.nnonzero = 0

        if depth < 0:
            self.depth = 0
            while self.N_b*2**self.depth < self.N:
                self.depth += 1
            self.N_padded = self.N_b*2**self.depth
            self.i_lower = 0
            self.j_lower = 0
            self.i_upper = self.N_padded
            self.j_upper = self.N_padded
        else:
            self.depth = depth
            self.N_padded = N_padded
            self.i_lower = i_lower
            self.j_lower = j_lower
            self.i_upper = i_upper
            self.j_upper = j_upper

    def print_tree_info(self):
        """Print some information about the Quadtree.
        """

        print("N =        %d" % (self.N))
        print("N_b =      %d" % (self.N_b))
        print("N_padded = %d" % (self.N_padded))
        print("depth =    %d" % (self.depth))
        print("gamma =    %1.3f" % (self.gamma))
        print("norm =     %e" % (math.sqrt(self.norm2)))
        print("nnonzero = %d" % (self.nnonzero))
        print("fillin =   %1.2f%%" % (100*self.nnonzero/self.N**2))

    def get_dense(self, tier=0, A=None):
        """Return the matrix as a dense matrix.
        """

        if tier == 0:
            A = numpy.matrix(numpy.zeros((self.N, self.N)))

        if tier < self.depth:
            for i in range(4):
                if self.children[i] is not None:
                    A = self.children[i].get_dense(tier+1, A)
        else:
            A[self.i_lower:self.i_upper, self.j_lower:self.j_upper] = self.submatrix

        return A

    def print_node(self):
        """Print information about this node.
        """

        print("node %s" % (repr(self)))
        print("N = %d, N_b = %d, N_pad = %d, tier = %d, depth = %d" % (
            self.N, self.N_b, self.N_padded, self.tier, self.depth))

    def get_nodes(self, depth, tier=0):
        """Get a list of nodes on tier == depth.

        The function returns a tuple of ([i], [j], width, norm).  The
        ([i], [j]) part is a list of the centers of the nodes, the width
        is the index width on this tier, and the norm is the norm of
        that submatrix..
        """

        half = (self.i_upper-self.i_lower)/2
        if tier < depth:
            x = []
            y = []
            norm = []
            width = 2*half

            for (i, j) in itertools.product(range(2), range(2)):
                if self.children[Quadtree.index(i, j)] is not None:
                    (x_temp, y_temp, norm_temp, width) = self.children[Quadtree.index(i, j)].get_nodes(
                        depth, tier=tier+1)
                    x = x+x_temp
                    y = y+y_temp
                    norm = norm+norm_temp
            return (x, y, norm, width)
        else:
            return ([self.i_lower+half], [self.j_lower+half], [math.sqrt(self.norm2)], 2*half)

    def get_element(self, i, j, i_lower=0, i_upper=0, j_lower=0, j_upper=0):
        """Get a matrix element Aij <- A(i,j).
        """

        if self.tier == 0:
            i_lower = 0
            i_upper = self.N_padded
            j_lower = 0
            j_upper = self.N_padded

        if self.tier < self.depth:
            i_half = (i_upper-i_lower)/2
            j_half = (j_upper-j_lower)/2

            for (i_child, j_child) in itertools.product(range(2), range(2)):
                i_l = i_lower+i_child*i_half
                i_u = i_lower+(i_child+1)*i_half
                j_l = j_lower+j_child*j_half
                j_u = j_lower+(j_child+1)*j_half

                if i >= i_l and i < i_u and j >= j_l and j < j_u:
                    if self.children[Quadtree.index(i_child, j_child)] is not None:
                        return self.children[Quadtree.index(i_child, j_child)].get_element(
                                i, j, i_lower=i_l, i_upper=i_u, j_lower=j_l, j_upper=j_u)
        else:
            if self.submatrix is not None:
                return self.submatrix[i-i_lower, j-j_lower]

        return 0

    def set_submatrix(self, i, j, A_b):
        """Set a submatrix A_b.

        The indices i and j are the lower left corner of the
        submatrix. The matrix size has to be (N_b, N_b).
        """

        if self.tier == 0:
            if A_b.shape != (self.N_b, self.N_b):
                raise Exception("submatrix shape incorrect")

        if self.tier < self.depth:
            i_half = (self.i_upper-self.i_lower)/2
            j_half = (self.j_upper-self.j_lower)/2

            for (i_child, j_child) in itertools.product(range(2), range(2)):
                i_l = self.i_lower+i_child*i_half
                i_u = self.i_lower+(i_child+1)*i_half
                j_l = self.j_lower+j_child*j_half
                j_u = self.j_lower+(j_child+1)*j_half

                if i >= i_l and i < i_u and j >= j_l and j < j_u:
                    if self.children[Quadtree.index(i_child, j_child)] is None:
                        self.children[Quadtree.index(i_child, j_child)] = \
                          Quadtree(N=self.N, N_b=self.N_b, N_padded=self.N_padded,
                                   i_lower=i_l, i_upper=i_u,
                                   j_lower=j_l, j_upper=j_u,
                                   gamma=self.gamma, tier=self.tier+1, depth=self.depth)
                    self.children[Quadtree.index(i_child, j_child)].set_submatrix(i, j, A_b)

            self.norm2 = 0
            self.nnonzero = 0
            for k in range(4):
                if self.children[k] is not None:
                    self.norm2 += self.children[k].norm2
                    self.nnonzero += self.children[k].nnonzero
        else:
            self.submatrix = A_b
            self.norm2 = numpy.linalg.norm(self.submatrix, 'fro')**2
            self.nnonzero = numpy.count_nonzero(self.submatrix)

    def print_matrix(self):
        """Print a matrix.
        """

        for i in range(self.N):
            for j in range(self.N):
                print("A(%d,%d) = %e" % (i, j, self.get_element(i, j)))

@mlab.show
def plot_octree(A, B, C, C_multiplication):
    """Set-up the 3D multiplication space.

    This forms the product space octree between two matrices A and B.
    """

    (x, y, norm, width) = A.get_nodes(A.depth)
    print(x, y, norm, width)
    mlab.points3d(x, y, numpy.zeros((len(x))), mode='2dsquare', scale_factor=width)

    (x, y, z, width) = C_multiplication.get_nodes(C_multiplication.depth)
    op = 0.6
    mlab.points3d(x, y, z, mode='cube', scale_factor=width, opacity=op)

    mlab.axes(xlabel='i', ylabel='j', zlabel='k')

if __name__ == "__main__":
    N = 32
    N_b = 16
    debug = True
    print_matrices = False

    #A = Quadtree(16000, N_b=16)

    A = Quadtree.random(N, N_b=N_b)
    if debug:
        A.print_tree_info()
        A_dense = A.get_dense()
        C_ref = numpy.dot(A_dense, A_dense)
        if print_matrices:
            print("A =")
            print(A_dense)
            print("A**2 = ")
            print(C_ref)

    C_multiplication = SpAMM()
    C = Quadtree.zeros(N, N_b=N_b)
    C_multiplication.multiply(A, A, C, tolerance=1e-6, symbolic_only=False)

    if debug:
        C_dense = C.get_dense()
        if print_matrices:
            print("C = ")
            print(C.get_dense())
        print("max diff. = %e" % (numpy.amax(numpy.absolute(C_ref-C_dense))))

    plot_octree(A, A, C, C_multiplication)
