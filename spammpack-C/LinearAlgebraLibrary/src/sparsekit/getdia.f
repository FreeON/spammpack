      subroutine getdia (nrow,ncol,job,a,ja,ia,len,diag,idiag,ioff)

c*********************************************************************72
c
cc GETDIA extracts a given diagonal from a matrix stored in CSR format.
c
c this subroutine extracts a given diagonal from a matrix stored in CSR
c format. The output matrix may be transformed with the diagonal removed
c from it if desired (as indicated by job.)
c
c Our definition of a diagonal of matrix is a vector of length nrow
c (always) which contains the elements in rows 1 to nrow of
c the matrix that are contained in the diagonal offset by ioff
c with respect to the main diagonal. If the diagonal element
c falls outside the matrix then it is defined as a zero entry.
c Thus the proper definition of diag(*) with offset ioff is
c
c     diag(k) = a(k,ioff+k) k=1,2,...,nrow
c     with elements falling outside the matrix being defined as zero.
c
c on entry:
c
c
c nrow      = integer. The row dimension of the matrix A.
c ncol      = integer. The column dimension of the matrix A.
c job   = integer. Job indicator.  If job = 0 then
c         the matrix a, ja, ia, is not altered on return.
c         if job.ne.1  then getdia will remove the entries
c         collected in diag from the original matrix.
c         This is done in place.
c
c a,ja,
c    ia = matrix stored in compressed sparse row a,ja,ia,format
c ioff  = integer,containing the offset of the wanted diagonal
c        the diagonal extracted is the one corresponding to the
c        entries a(i,j) with j-i = ioff.
c        thus ioff = 0 means the main diagonal
c
c on return:
c
c len   = number of nonzero elements found in diag.
c         (len .le. min(nrow,ncol-ioff)-max(1,1-ioff) + 1 )
c
c diag  = double precision array of length nrow containing the wanted diagonal.
c        diag contains the diagonal (a(i,j),j-i = ioff ) as defined
c         above.
c
c idiag = integer array of  length len, containing the poisitions
c         in the original arrays a and ja of the diagonal elements
c         collected in diag. A zero entry in idiag(i) means that
c         there was no entry found in row i belonging to the diagonal.
c
c a, ja,
c    ia = if job .ne. 0 the matrix is unchanged. otherwise the nonzero
c         diagonal entries collected in diag are removed from the
c         matrix. the structure is modified since the diagonal elements
c        are removed from a,ja,ia. Thus, the  returned matrix will
c         have len fewer elements if the diagonal is full.
c
c           Y. Saad, Sep. 21 1989 - modified and tested May 9, 1990.
c
      implicit double precision (a-h,o-z)
      double precision diag(*),a(*)
      integer nrow, ncol, job, len, ia(*), ja(*), idiag(*)

      integer istart, max, iend, i, kold, k, kdiag, ko

      istart = max(0,-ioff)
      iend = min0(nrow,ncol-ioff)
      len = 0
      do 1 i=1,nrow
         idiag(i) = 0
       diag(i) = 0.0
 1    continue
c
c     extract  diagonal elements
c
      do 6 i=istart+1, iend
         do 51 k= ia(i),ia(i+1) -1
            if (ja(k)-i .eq. ioff) then
               diag(i)= a(k)
               idiag(i) = k
               len = len+1
               goto 6
            end if
 51      continue
 6    continue
      if (job .eq. 0 .or. len .eq.0) return
c
c  rewind the structure
c
      ko = 0
      do  7 i=istart+1,iend
         kold = ko
         kdiag = idiag(i)
         if (kdiag .eq. 0) goto 7
         do 71 k= ia(i), ia(i+1)-1
            if (ja(k) .eq. kdiag) goto 71
            ko = ko+1
            a(ko) = a(k)
            ja(ko) = ja(k)
 71      continue
         ia(i) = kold+1
 7    continue
c
c redefine ia(nrow+1)
c
      ia(nrow+1) = ko+1
      return
      end
