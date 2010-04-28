      subroutine lusol0 (n, y, x, alu, jlu, ju)

c*********************************************************************72
c
cc LUSOL0 performs forward and backward solves for LU matrix produced by ILUT.
c
c performs a forward followed by a backward solve
c for LU matrix as produced by  ILUT
c
      INTEGER N
        double precision x(n), y(n), alu(*)
      integer jlu(*), ju(*)
        integer i,k
c
c forward solve
c
        do 40 i = 1, n
           x(i) = y(i)
           do 41 k=jlu(i),ju(i)-1
              x(i) = x(i) - alu(k)* x(jlu(k))
 41        continue
 40     continue
c
c     backward solve.
c
      do 90 i = n, 1, -1
         do 91 k=ju(i),jlu(i+1)-1
              x(i) = x(i) - alu(k)*x(jlu(k))
 91         continue
           x(i) = alu(i)*x(i)
 90     continue

        return
      end
