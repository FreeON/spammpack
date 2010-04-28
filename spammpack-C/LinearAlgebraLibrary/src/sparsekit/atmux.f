      subroutine atmux (n, x, y, a, ja, ia)

c*********************************************************************72
c
cc ATMUX computes A' * x for a CSR matrix A.
c
c multiplies the transpose of a matrix by a vector when the original
c matrix is stored in compressed sparse row storage. Can also be
c viewed as the product of a matrix by a vector when the original
c matrix is stored in the compressed sparse column format.
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
c y     = double precision array of length n, containing the product y=transp(A)*x
c
c
      double precision x(*), y(*), a(*)
      integer n, ia(*), ja(*)
      integer i, k
c
c     zero out output vector
c
      do 1 i=1,n
         y(i) = 0.0
 1    continue
c
c loop over the rows
c
      do 100 i = 1,n
         do 99 k=ia(i), ia(i+1)-1
            y(ja(k)) = y(ja(k)) + x(i)*a(k)
 99      continue
 100  continue

      return
      end
