      subroutine lsol (n,x,y,al,jal,ial)

c*********************************************************************72
c
cc LSOL solves L*x = y ; L = lower unit triang. /  CSR format
c
c solves a unit lower triangular system by standard (sequential )
c forward elimination - matrix stored in CSR format.
c
c On entry:
c
c n      = integer. dimension of problem.
c y      = double precision array containg the right side.
c
c al,
c jal,
c ial,    = Lower triangular matrix stored in compressed sparse row
c          format.
c
c On return:
c
c      x  = The solution of  L x  = y.
c
      integer n, jal(*),ial(n+1)
      double precision  x(n), y(n), al(*)
      integer k, j
      double precision  t

      x(1) = y(1)
      do 150 k = 2, n
         t = y(k)
         do 100 j = ial(k), ial(k+1)-1
            t = t-al(j)*x(jal(j))
 100     continue
         x(k) = t
 150  continue

      return
      end
