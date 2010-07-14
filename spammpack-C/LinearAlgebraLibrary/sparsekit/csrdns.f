      subroutine csrdns(nrow,ncol,a,ja,ia,dns,ndns,ierr)

c*********************************************************************72
c
cc CSRDNS converts Compressed Sparse Row to Dense format.
c
c converts a row-stored sparse matrix into a densely stored one
c
c On entry:
c
c
c nrow      = row-dimension of a
c ncol      = column dimension of a
c a,
c ja,
c ia    = input matrix in compressed sparse row format.
c         (a=value array, ja=column array, ia=pointer array)
c dns   = array where to store dense matrix
c ndns      = first dimension of array dns
c
c on return:
c
c dns   = the sparse matrix a, ja, ia has been stored in dns(ndns,*)
c
c ierr  = integer error indicator.
c         ierr .eq. 0  means normal return
c         ierr .eq. i  means that the code has stopped when processing
c         row number i, because it found a column number .gt. ncol.
c
      double precision dns(ndns,*),a(*)
      integer ja(*),ia(*)

      ierr = 0
      do 1 i=1, nrow
         do 2 j=1,ncol
          dns(i,j) = 0.0
 2       continue
 1    continue

      do 4 i=1,nrow
         do 3 k=ia(i),ia(i+1)-1
            j = ja(k)
          if (j .gt. ncol) then
               ierr = i
               return
          end if
          dns(i,j) = a(k)
 3       continue
 4    continue
      return
      end
