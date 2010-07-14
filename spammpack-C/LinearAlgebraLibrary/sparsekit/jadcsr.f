      subroutine jadcsr (nrow, idiag, a, ja, ia, iperm, ao, jao, iao)

c*********************************************************************72
c
cc JADSCR converts Jagged Diagonal Storage to Compressed Sparse Row.
c
c this subroutine converts a matrix stored in the jagged diagonal format
c to the compressed sparse row format.
c
c on entry:
c
c nrow         = integer. the row dimension of the matrix A.
c
c idiag   = integer. The  number of jagged diagonals in the data
c           structure a, ja, ia.
c
c a,
c ja,
c ia      = input matrix in jagged diagonal format.
c
c iperm   = permutation of the rows used to obtain the JAD ordering.
c
c on return:
c
c
c ao, jao,
c iao     = matrix in CSR format.
c
c determine first the pointers for output matrix. Go through the
c structure once:
c
      integer ja(*), jao(*), ia(idiag+1), iperm(nrow), iao(nrow+1)
      double precision a(*), ao(*)

      do 137 j=1,nrow
         jao(j) = 0
 137  continue
c
c     compute the lengths of each row of output matrix
c
      do 140 i=1, idiag
         len = ia(i+1)-ia(i)
         do 138 k=1,len
            jao(iperm(k)) = jao(iperm(k))+1
 138     continue
 140  continue
c
c     remember to permute
c
      kpos = 1
      iao(1) = 1
      do 141 i=1, nrow
         kpos = kpos+jao(i)
         iao(i+1) = kpos
 141  continue
c
c     copy elemnts one at a time.
c
      do 200 jj = 1, idiag
         k1 = ia(jj)-1
         len = ia(jj+1)-k1-1
         do 160 k=1,len
            kpos = iao(iperm(k))
            ao(kpos) = a(k1+k)
            jao(kpos) = ja(k1+k)
            iao(iperm(k)) = kpos+1
 160     continue
 200  continue
c
c     rewind pointers
c
      do 5 j=nrow,1,-1
         iao(j+1) = iao(j)
 5    continue
      iao(1) = 1
      return
      end
