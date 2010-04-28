      subroutine ldsoll (n,x,y,al,jal,nlev,lev,ilev)

c*********************************************************************72
c
cc LDSOLL solves L*x = y; L = triangular.
c
c On entry:
c
c n      = integer. dimension of problem.
c y      = double precision array containg the right hand side.
c
c al,
c jal,   = Lower triangular matrix stored in Modified Sparse Row
c          format.
c nlev   = number of levels in matrix
c lev    = integer array of length n, containing the permutation
c          that defines the levels in the level scheduling ordering.
c ilev   = pointer to beginning of levels in lev.
c          the numbers lev(i) to lev(i+1)-1 contain the row numbers
c          that belong to level number i, in the level shcheduling
c          ordering.
c
c On return:
c
c      x = The solution of  L x = y .
      integer n, nlev, jal(*), ilev(nlev+1), lev(n)
      double precision x(n), y(n), al(*)
      integer ii, jrow, i
      double precision t
c
c     outer loop goes through the levels. (SEQUENTIAL loop)
c
      do 150 ii=1, nlev
c
c     next loop executes within the same level. PARALLEL loop
c
         do 100 i=ilev(ii), ilev(ii+1)-1
            jrow = lev(i)
c
c compute inner product of row jrow with x
c
            t = y(jrow)
            do 130 k=jal(jrow), jal(jrow+1)-1
               t = t - al(k)*x(jal(k))
 130        continue
            x(jrow) = t*al(jrow)
 100     continue
 150  continue
      return
      end
