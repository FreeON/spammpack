      subroutine udsol (n,x,y,au,jau)

c*********************************************************************72
c
cc UDSOL solves U*x = y;   U = upper triangular in MSR format
c
c solves a non-unit upper triangular matrix by standard (sequential )
c backward elimination - matrix stored in MSR format.
c with diagonal elements already inverted (otherwise do inversion,
c au(1:n) = 1.0/au(1:n),  before calling).
c
c On entry:
c
c n      = integer. dimension of problem.
c y      = double precision array containg the right side.
c
c au,
c jau,    = Lower triangular matrix stored in modified sparse row
c          format.
c
c On return:
c
c      x = The solution of  U x = y .
c
      integer n, jau(*)
      double precision  x(n), y(n),au(*)
      integer k, j
      double precision t

      x(n) = y(n)*au(n)
      do k = n-1,1,-1
         t = y(k)
         do j = jau(k), jau(k+1)-1
            t = t - au(j)*x(jau(j))
         end do
         x(k) = au(k)*t
      end do

      return
      end
