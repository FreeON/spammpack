      subroutine lsolc (n,x,y,al,jal,ial)

c*********************************************************************72
c
cc LSOLC solves L*x = y where L = unit lower triang. CSC format
c
c solves a unit lower triangular system by standard (sequential )
c forward elimination - matrix stored in CSC format.
c
c On entry:
c
c n      = integer. dimension of problem.
c y      = double precision array containg the right side.
c
c al,
c jal,
c ial,    = Lower triangular matrix stored in compressed sparse column
c          format.
c
c On return:
c
c      x  = The solution of  L x  = y.
c
      integer n, jal(*),ial(*)
      double precision  x(n), y(n), al(*)
      integer k, j
      double precision t

      do 140 k=1,n
         x(k) = y(k)
 140  continue
      do 150 k = 1, n-1
         t = x(k)
         do 100 j = ial(k), ial(k+1)-1
            x(jal(j)) = x(jal(j)) - t*al(j)
 100     continue
 150  continue
c
      return
      end
