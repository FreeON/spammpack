      subroutine udsolc (n,x,y,au,jau)

c*********************************************************************72
c
cc UDSOLC solves U * x = y, for upper triangular U in MSC format.
c
c solves a (non-unit) upper triangular system by standard (sequential)
c forward elimination - matrix stored in Modified Sparse Column format
c with diagonal elements already inverted (otherwise do inversion,
c auuuul(1:n) = 1.0/au(1:n),  before calling ldsol).
c
c
c On entry:
c
c n      = integer. dimension of problem.
c y      = double precision array containg the right hand side.
c
c au,
c jau,   = Upper triangular matrix stored in Modified Sparse Column
c          format.
c
c On return:
c
c      x = The solution of  U x = y .
c
      integer n, jau(*)
      double precision x(n), y(n), au(*)
      integer k, j
      double precision t

      do k=1,n
         x(k) = y(k)
      end do

      do 150 k = n,1,-1
         x(k) = x(k)*au(k)
         t = x(k)
         do 100 j = jau(k), jau(k+1)-1
            x(jau(j)) = x(jau(j)) - t*au(j)
 100     continue
 150  continue
c
      return
      end
