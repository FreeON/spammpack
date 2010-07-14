      subroutine ldsol (n,x,y,al,jal)

c*********************************************************************72
c
cc LDSOL solves L * x = y, for L a triangular matrix in MSR format.
c
c solves a (non-unit) lower triangular system by standard (sequential)
c forward elimination - matrix stored in MSR format
c with diagonal elements already inverted (otherwise do inversion,
c al(1:n) = 1.0/al(1:n),  before calling ldsol).
c
c On entry:
c
c n      = integer. dimension of problem.
c y      = double precision array containg the right hand side.
c
c al,
c jal,   = Lower triangular matrix stored in Modified Sparse Row
c          format.
c
c On return:
c
c      x = The solution of  L x = y .
c
      integer n, jal(*)
      double precision x(n), y(n), al(*)
      integer k, j
      double precision t

      x(1) = y(1)*al(1)
      do 150 k = 2, n
         t = y(k)
         do 100 j = jal(k), jal(k+1)-1
            t = t - al(j)*x(jal(j))
 100     continue
         x(k) = al(k)*t
 150  continue
      return
      end
