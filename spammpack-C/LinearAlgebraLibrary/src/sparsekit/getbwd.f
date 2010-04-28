      subroutine getbwd(n,a,ja,ia,ml,mu)

c*********************************************************************72
c
cc GETBWD gets the bandwidth of lower part and upper part of A.
c
c  Discussion:
c
c    This routine does not assume that A is sorted.
c
c on entry:
c
c n      = integer = the row dimension of the matrix
c a, ja,
c    ia = matrix in compressed sparse row format.
c
c on return:
c
c ml      = integer. The bandwidth of the strict lower part of A
c mu      = integer. The bandwidth of the strict upper part of A
c
c Notes:
c ===== ml and mu are allowed to be negative or return. This may be
c       useful since it will tell us whether a band is confined
c       in the strict  upper/lower triangular part.
c       indeed the definitions of ml and mu are
c
c       ml = max ( (i-j)  s.t. a(i,j) .ne. 0  )
c       mu = max ( (j-i)  s.t. a(i,j) .ne. 0  )
c
c Y. Saad, Sep. 21 1989                                                c
c
      double precision a(*)
      integer ja(*),ia(n+1),ml,mu,ldist,i,k
      ml = - n
      mu = - n
      do 3 i=1,n
         do 31 k=ia(i),ia(i+1)-1
            ldist = i-ja(k)
            ml = max(ml,ldist)
            mu = max(mu,-ldist)
 31      continue
 3    continue
      return
      end
