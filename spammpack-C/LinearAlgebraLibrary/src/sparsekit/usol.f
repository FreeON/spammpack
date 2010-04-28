      subroutine usol (n,x,y,au,jau,iau)

c*********************************************************************72
c
cc USOL solves   U x = y    U = unit upper triangular.
c
c solves a unit upper triangular system by standard (sequential )
c backward elimination - matrix stored in CSR format.
c
c
c On entry:
c
c n      = integer. dimension of problem.
c y      = double precision array containg the right side.
c
c au,
c jau,
c iau,    = Lower triangular matrix stored in compressed sparse row
c          format.
c
c On return:
c
c      x = The solution of  U x = y .
c
      integer n, jau(*),iau(n+1)
      double precision  x(n), y(n), au(*)
      integer k, j
      double precision  t

      x(n) = y(n)
      do 150 k = n-1,1,-1
         t = y(k)
         do 100 j = iau(k), iau(k+1)-1
            t = t - au(j)*x(jau(j))
 100     continue
         x(k) = t
 150  continue

      return
      end
