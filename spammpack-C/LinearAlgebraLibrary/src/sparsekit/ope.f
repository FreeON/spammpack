      subroutine ope (n, x, y, a, ja, ia)

c*********************************************************************72
c
cc OPE sparse matrix * vector multiplication
c

      INTEGER N
      double precision  x(n), y(n), a(*)
      integer k1, k2, ja(*), ia(n+1)

      do 100 i=1,n
           k1 = ia(i)
           k2 = ia(i+1) -1
           y(i) = 0.0
           do 99 k=k1, k2
 99           y(i) = y(i) + a(k) * x(ja(k))
 100      continue
      return
      end
