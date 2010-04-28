      subroutine amux(n, x, y, a,ja,ia)

c*********************************************************************72
c
cc AMUX multiplies a CSR matrix A times a vector.
c
c multiplies a matrix by a vector using the dot product form
c Matrix A is stored in compressed sparse row storage.
c
c on entry:
c
c n     = row dimension of A
c x     = double precision array of length equal to the column dimension of
c         the A matrix.
c a, ja,
c    ia = input matrix in compressed sparse row format.
c
c on return:
c
c y     = double precision array of length n, containing the product y=Ax
c
      double precision  x(*), y(*), a(*)
      integer n, ja(*), ia(*)
      double precision t
      integer i, k

      do 100 i = 1,n
c
c     compute the inner product of row i with vector x
c
         t = 0.0
         do 99 k=ia(i), ia(i+1)-1
            t = t + a(k)*x(ja(k))
 99      continue
c
c     store result in y(i)
c
         y(i) = t
 100  continue

      return
      end
