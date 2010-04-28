      subroutine amuxj (n, x, y, jdiag, a, ja, ia)

c*********************************************************************72
c
cc AMUXJ multiplies a JAD matrix times a vector.
c
c multiplies a matrix by a vector when the original matrix is stored
c in the jagged diagonal storage format.
c
c on entry:
c
c n      = row dimension of A
c x      = double precision array of length equal to the column dimension of
c         the A matrix.
c jdiag  = integer. The number of jadded-diagonals in the data-structure.
c a      = double precision array containing the jadded diagonals of A stored
c          in succession (in decreasing lengths)
c j      = integer array containing the colum indices of the
c          corresponding elements in a.
c ia     = integer array containing the lengths of the  jagged diagonals
c
c on return:
c
c y      = double precision array of length n, containing the product y=A*x
c
c Note:
c
c Permutation related to the JAD format is not performed.
c this can be done by:
c     call permvec (n,y,y,iperm)
c after the call to amuxj, where iperm is the permutation produced
c by csrjad.
c
      integer n, jdiag, ja(*), ia(*)
      double precision x(n), y(n), a(*)
      integer i, ii, k1, len, j

      do 1 i=1, n
         y(i) = 0.0
 1    continue
      do 70 ii=1, jdiag
         k1 = ia(ii)-1
         len = ia(ii+1)-k1-1
         do 60 j=1,len
            y(j)= y(j)+a(k1+j)*x(ja(k1+j))
 60      continue
 70   continue
      return
      end
