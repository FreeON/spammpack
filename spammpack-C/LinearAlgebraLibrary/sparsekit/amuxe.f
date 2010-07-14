      subroutine amuxe(n,x,y,na,ncol,a,ja)

c*********************************************************************72
c
cc AMUXE multiplies an ELL matrix times a vector.
c
c multiplies a matrix by a vector when the original matrix is stored
c in the ellpack-itpack sparse format.
c
c on entry:
c
c n     = row dimension of A
c x     = double precision array of length equal to the column dimension of
c         the A matrix.
c na    = integer. The first dimension of arrays a and ja
c         as declared by the calling program.
c ncol  = integer. The number of active columns in array a.
c         (i.e., the number of generalized diagonals in matrix.)
c a, ja = the real and integer arrays of the itpack format
c         (a(i,k),k=1,ncol contains the elements of row i in matrix
c          ja(i,k),k=1,ncol contains their column numbers)
c
c on return:
c
c y     = double precision array of length n, containing the product y=y=A*x
c
      INTEGER N,NA
      double precision x(n), y(n), a(na,*)
      integer ncol, ja(na,*)
      integer i, j

      do 1 i=1, n
         y(i) = 0.0
 1    continue
      do 10 j=1,ncol
         do 25 i = 1,n
            y(i) = y(i)+a(i,j)*x(ja(i,j))
 25      continue
 10   continue
      return
      end
