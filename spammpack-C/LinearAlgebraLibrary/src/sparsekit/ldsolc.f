      subroutine ldsolc (n,x,y,al,jal)

c*********************************************************************72
c
cc LDSOLC solves L*x = y;    L = nonunit Low. Triang. MSC format
c
c solves a (non-unit) lower triangular system by standard (sequential)
c forward elimination - matrix stored in Modified Sparse Column format
c with diagonal elements already inverted (otherwise do inversion,
c al(1:n) = 1.0/al(1:n),  before calling ldsol).
c
c
c On entry:
c
c n      = integer. dimension of problem.
c y      = double precision array containg the right hand side.
c
c al,
c jal,
c ial,    = Lower triangular matrix stored in Modified Sparse Column
c           format.
c
c On return:
c
c      x = The solution of  L x = y .
c
      integer n, jal(*)
      double precision x(n), y(n), al(*)
      integer k, j
      double precision t

      do 140 k=1,n
         x(k) = y(k)
 140  continue
      do 150 k = 1, n
         x(k) = x(k)*al(k)
         t = x(k)
         do 100 j = jal(k), jal(k+1)-1
            x(jal(j)) = x(jal(j)) - t*al(j)
 100     continue
 150  continue
c
      return
      end
