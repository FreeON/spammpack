      subroutine ssrcsr (nrow,a,ja,ia,nzmax,ao,jao,iao,indu,ierr)

c*********************************************************************72
c
cc SSRCSR converts Symmetric Sparse Row to (regular) Compressed Sparse Row.
c
c this subroutine converts  a symmetric  matrix in which only the lower
c part is  stored in compressed sparse row format, i.e.,
c a matrix stored in symmetric sparse format, into a fully stored matrix
c i.e., a matrix where both the lower and upper parts are stored in
c compressed sparse row format. the algorithm is in place (i.e. result
c may be overwritten onto the input matrix a, ja, ia ).
c the output matrix delivered by ssrcsr is such that each row starts with
c the elements of the lower part followed by those of the upper part.
c
c on entry:
c
c
c nrow  = row dimension of inout matrix
c a,
c ia,
c ja    = matrix in compressed sparse row format. This is assumed to be
c         a lower triangular matrix.
c
c nzmax      = size of arrays ao and jao. ssrcsr will abort if the storage
c         provided in a, ja is not sufficient to store A. See ierr.
c
c on return:
c
c ao, iao,
c   jao = output matrix in compressed sparse row format. The resulting
c         matrix is symmetric and is equal to A+A**T - D, if
c         A is the original lower triangular matrix. ao, jao, iao,
c         can be the same as a, ja, ia in the calling sequence.
c
c indu  = integer array of length nrow+1. If the input matrix is such
c         that the last element in each row is its diagonal element then
c         on return, indu will contain the pointers to the diagonal
c         element in each row of the output matrix. Otherwise used as
c         work array.
c ierr  = integer. Serving as error message. If the length of the arrays
c         ao, jao exceeds nzmax, ierr returns the minimum value
c         needed for nzmax. otherwise ierr=0 (normal return).
c
      integer ia(nrow+1),iao(nrow+1),ja(*),jao(nzmax),indu(nrow+1)
      double precision a(*),ao(nzmax)

      ierr = 0
      do 1 i=1,nrow+1
         indu(i) = 0
 1    continue
c
c     compute  number of elements in each row of strict upper part.
c     put result in indu(i+1)  for row i.
c
      do 3 i=1, nrow
         do 2 k=ia(i),ia(i+1)-1
            j = ja(k)
            if (j .lt. i) indu(j+1) = indu(j+1)+1
 2       continue
 3    continue
c
c     find addresses of first elements of ouput matrix. result in indu
c
      indu(1) = 1
      do 4 i=1,nrow
         lenrow = ia(i+1)-ia(i)
         indu(i+1) = indu(i) + indu(i+1) + lenrow
 4    continue
c  enough storage in a, ja ?
c
      nnz = indu(nrow+1)-1
      if (nnz .gt. nzmax) then
         ierr = nnz
         return
      end if
c
c     now copy lower part (backwards).
c
      kosav = indu(nrow+1)
      do 6 i=nrow,1,-1
         klast = ia(i+1)-1
         kfirst = ia(i)
         iao(i+1) = kosav
         ko = indu(i)
         kosav = ko
         do 5 k = kfirst, klast
            ao(ko) = a(k)
            jao(ko) = ja(k)
          ko = ko+1
 5       continue
         indu(i) = ko
 6    continue
      iao(1) = 1
c
c     now copy upper part. Go through the structure of ao, jao, iao
c     that has already been copied (lower part). indu(i) is the address
c     of the next free location in row i for ao, jao.
c
      do 8 i=1,nrow
c     i-th row is now in ao, jao, iao structure: lower half part
         do 9 k=iao(i), iao(i+1)-1
            j = jao(k)
            if (j .ge. i)  goto 8
            ipos = indu(j)
            ao(ipos) = ao(k)
            jao(ipos) = i
            indu(j) = indu(j) + 1
 9       continue
 8    continue
      return
      end
