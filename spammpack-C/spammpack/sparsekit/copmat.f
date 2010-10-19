      subroutine copmat (nrow,a,ja,ia,ao,jao,iao,ipos)

c*********************************************************************72
c
cc COPMAT copies the matrix A, JA, IA, into the matrix AO, JAO, IAO.
c
c on entry:
c
c nrow      = row dimension of the matrix
c a,
c ja,
c ia    = input matrix in compressed sparse row format.
c ipos  = integer. indicates the position in the array ao, jao
c         where the first element should be copied. Thus
c         iao(1) = ipos on return.
c
c on return:
c
c ao,
c jao,
c iao   = output matrix containing the same data as a, ja, ia.
c
c           Y. Saad, March 1990.
c
      double precision a(*),ao(*)
      integer nrow, ia(*),ja(*),jao(*),iao(*), ipos
      integer kst, i, k

      kst    = ipos -ia(1)
      do 100 i = 1, nrow+1
         iao(i) = ia(i) + kst
 100  continue

      do 200 k=ia(1), ia(nrow+1)-1
         ao(kst+k) = a(k)
         jao(kst+k)= ja(k)
 200  continue

      return
      end
