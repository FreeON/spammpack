      subroutine getl (n,a,ja,ia,ao,jao,iao)

c*********************************************************************72
c
cc GETL extracts the lower triangular part of a matrix.
c
c This routine extracts the lower triangular part of a matrix
c and writes the result ao, jao, iao. The routine is in place in
c that ao, jao, iao can be the same as a, ja, ia if desired.
c
c on input:
c
c n     = dimension of the matrix a.
c a, ja,
c    ia = matrix stored in compressed sparse row format.
c On return:
c ao, jao,
c    iao = lower triangular matrix (lower part of a)
c      stored in a, ja, ia, format
c note: the diagonal element is the last element in each row.
c i.e. in  a(ia(i+1)-1 )
c ao, jao, iao may be the same as a, ja, ia on entry -- in which case
c getl will overwrite the result on a, ja, ia.
c
      integer n, ia(*), ja(*), iao(*), jao(*)
      double precision a(*), ao(*)
      double precision t
      integer ko, kold, kdiag, k, i
c
c  inititialize ko (pointer for output matrix)
c
      ko = 0
      do  7 i=1, n
         kold = ko
         kdiag = 0
         do 71 k = ia(i), ia(i+1) -1
            if (ja(k)  .gt. i) goto 71
            ko = ko+1
            ao(ko) = a(k)
            jao(ko) = ja(k)
            if (ja(k)  .eq. i) kdiag = ko
 71      continue
         if (kdiag .eq. 0 .or. kdiag .eq. ko) goto 72
c
c  exchange
c
         t = ao(kdiag)
         ao(kdiag) = ao(ko)
         ao(ko) = t

         k = jao(kdiag)
         jao(kdiag) = jao(ko)
         jao(ko) = k
 72      iao(i) = kold+1
 7    continue
c     redefine iao(n+1)
      iao(n+1) = ko+1
      return
      end
