      subroutine amuxd (n,x,y,diag,ndiag,idiag,ioff)

c*********************************************************************72
c
cc AMUXD multiplies a DIA matrix times a vector.
c
c multiplies a matrix by a vector when the original matrix is stored
c in the diagonal storage format.
c
c on entry:
c
c n     = row dimension of A
c x     = double precision array of length equal to the column dimension of
c         the A matrix.
c ndiag  = integer. The first dimension of array adiag as declared in
c         the calling program.
c idiag  = integer. The number of diagonals in the matrix.
c diag   = double precision array containing the diagonals stored of A.
c idiag  = number of diagonals in matrix.
c diag   = double precision array of size (ndiag x idiag) containing the diagonals
c
c ioff   = integer array of length idiag, containing the offsets of the
c            diagonals of the matrix:
c          diag(i,k) contains the element a(i,i+ioff(k)) of the matrix.
c
c on return:
c
c y     = double precision array of length n, containing the product y=A*x
c
      integer n, ndiag, idiag, ioff(idiag)
      double precision x(n), y(n), diag(ndiag,idiag)
      integer j, k, io, i1, i2

      do 1 j=1, n
         y(j) = 0.0
 1    continue
      do 10 j=1, idiag
         io = ioff(j)
         i1 = max0(1,1-io)
         i2 = min0(n,n-io)
         do 9 k=i1, i2
            y(k) = y(k)+diag(k,j)*x(k+io)
 9       continue
 10   continue

      return
      end
