      subroutine ecn(N,IC,NE,IA,JA,AR,NN,IERR)

c*********************************************************************72
c
cc ECN generates sparse (square) matrices of the type E(N,C).
c
c   PURPOSE
c
c   The subroutine generates sparse (square) matrices of the type
c   E(N,C).  This type of matrix has the following characteristics:
c   Symmetric, positive-definite, N x N matrices with 4 in the diagonal
c   and -1 in the two sidediagonal and in the two bands at the distance
c   C from the diagonal. These matrices are similar to matrices obtained
c   from using the five-point formula in the discretization of the
c   elliptic PDE.
c
c
c   Note: If A is the sparse matrix of type E(N,C), then
c
c       min|A(i,j)| = 1,     max|A(i,j)| = 4
c
c
c   CONTRIBUTOR: Ernest E. Rothman
c                Cornell Theory Center/Cornell National Supercomputer
c                Facility.
c                e-mail address: BITNET:   eer@cornellf
c                                INTERNET: eer@cornellf.tn.cornell.edu
c
c   REFERENCE
c
c   1) Zlatev, Zahari; Schaumburg, Kjeld; Wasniewski, Jerzy;
c      "A Testing Scheme for Subroutines Solving Large Linear Problems",
c       Computers and Chemistry, Vol. 5, No. 2-3, pp. 91-100, 1981.
c   2) Osterby, Ole and Zletev, Zahari;
c      "Direct Methods for Sparse Matrices";
c       Springer-Verlag 1983.
c
c
c   INPUT PARAMETERS
c
c   N    - Integer. The size of the square matrix.
c          N > 2 must be specified.
c
c   NN   - Integer. The dimension of integer arrays IA and JA and
c          double precision array AR. Must be at least NE.
c
c   NN  - Integer. The dimension of integer array JA. Must be at least
c          NE.
c
c   IC   - Integer. The sparsity pattern can be changed by means of this
c          parameter.  1 < IC < N   must be specified.
c
c
c
c   OUTPUT PARAMETERS
c
c   NE   - Integer. The number of nonzero elements in the sparse matrix
c          of the type E(N,C). NE = 5*N - 2*IC - 2 .
c
c   AR(NN)  - Real array.
c             Stored entries of the sparse matrix A.
c             NE is the number of nonzeros including a mandatory
c             diagonal entry for each row.
c
c   IA(NN)  - Integer array.(Double precision)
c             Pointers to specify rows for the stored nonzero entries
c             in AR.
c
c   JA(NN) - Integer array.
c             Pointers to specify columns for the stored nonzero entries
c             in AR.
c
c   IERR    - Error parameter is returned as zero on successful
c             execution of the subroutine.
c             Error diagnostics are given by means of positive values
c             of this parameter as follows:
c             IERR = 1    -  N       is out of range.
c             IERR = 2    -  IC      is out of range.
c             IERR = 3    -  NN      is out of range.
c
      double precision ar(nn)
      integer ia(nn), ja(nn), n, ne, ierr
      ierr = 0
c
c  check the input parameters:
c
      if(n.le.2)then
         ierr = 1
         return
      end if
      if(ic.le.1.or.ic.ge.n)then
         ierr = 2
         return
      end if

      ne = 5*n-2*ic-2
      if(nn.lt.ne)then
         ierr = 3
         return
      end if
c
c Begin to generate the nonzero elements as well as the row and column
c pointers:
c
      do 20 i=1,n
      ar(i) = 4.0
      ia(i) = i
      ja(i) = i
20    continue
      ilast = n
      do 30 i=1,n-1
      it = ilast + i
      ar(it) = -1.0
      ia(it) = i+1
      ja(it) = i
30    continue
      ilast = ilast + n - 1
      do 40 i=1,n-1
      it = ilast + i
      ar(it) = -1.0
      ia(it) = i
      ja(it) = i+1
40    continue
      ilast = ilast + n-1
      do 50 i=1,n-ic
      it = ilast + i
      ar(it) = -1.0
      ia(it) = i+ic
      ja(it) = i
50    continue
      ilast = ilast + n-ic
      do 60 I=1,n-ic
      it = ilast + i
      ar(it) = -1.0
      ia(it) = i
      ja(it) = i+ic
60    continue
c      ilast = ilast + n-ic
c      if(ilast.ne.5*n-2*ic-2) then
c      write(*,*)' ilast equal to ', ilast
c      write(*,*)' ILAST, the no. of nonzeros, should = ', 5*n-2*ic-2
c      stop
c      end if
c
      return
      end
