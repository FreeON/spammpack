      subroutine usolc (n,x,y,au,jau,iau)

c*********************************************************************72
c
cc USOLC solves U * x = y for unit upper triangular U in CSC format.
c
c solves a unit upper triangular system by standard (sequential )
c forward elimination - matrix stored in CSC format.
c
c
c On entry:
c
c n      = integer. dimension of problem.
c y      = double precision array containg the right side.
c
c au,
c jau,
c iau,    = Uower triangular matrix stored in compressed sparse column
c          format.
c
c On return:
c
c      x  = The solution of  U x  = y.
c
      double precision  x(*), y(*), au(*)
      integer n, jau(*),iau(*)
      integer k, j
      double precision t
c
      do 140 k=1,n
         x(k) = y(k)
 140  continue
      do 150 k = n,1,-1
         t = x(k)
         do 100 j = iau(k), iau(k+1)-1
            x(jau(j)) = x(jau(j)) - t*au(j)
 100     continue
 150  continue

      return
      end
