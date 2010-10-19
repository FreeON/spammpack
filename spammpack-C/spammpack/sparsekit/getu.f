      subroutine getu (n,a,ja,ia,ao,jao,iao)

c*********************************************************************72
c
cc GETU extracts the upper triangular part of a matrix.
c
c This routine extracts the upper triangular part of a matrix
c and writes the result ao, jao, iao. The routine is in place in
c that ao, jao, iao can be the same as a, ja, ia if desired.
c
c on input:
c
c n     = dimension of the matrix a.
c a, ja,
c    ia = matrix stored in a, ja, ia, format
c On return:
c ao, jao,
c    iao = upper triangular matrix (upper part of a)
c      stored in compressed sparse row format
c note: the diagonal element is the last element in each row.
c i.e. in  a(ia(i+1)-1 )
c ao, jao, iao may be the same as a, ja, ia on entry -- in which case
c getu will overwrite the result on a, ja, ia.
c
      integer n, ia(*), ja(*), iao(*), jao(*)
      double precision a(*), ao(*)
      double precision t
      integer ko, k, i, kdiag, kfirst

      ko = 0
      do  7 i=1, n
         kfirst = ko+1
         kdiag = 0
         do 71 k = ia(i), ia(i+1) -1
            if (ja(k)  .lt. i) goto 71
            ko = ko+1
            ao(ko) = a(k)
            jao(ko) = ja(k)
            if (ja(k)  .eq. i) kdiag = ko
 71      continue
         if (kdiag .eq. 0 .or. kdiag .eq. kfirst) goto 72
c     exchange
         t = ao(kdiag)
         ao(kdiag) = ao(kfirst)
         ao(kfirst) = t

         k = jao(kdiag)
         jao(kdiag) = jao(kfirst)
         jao(kfirst) = k
 72      iao(i) = kfirst
 7    continue
c     redefine iao(n+1)
      iao(n+1) = ko+1
      return
      end
