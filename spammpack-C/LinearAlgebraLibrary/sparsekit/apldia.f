      subroutine apldia (nrow, job, a, ja, ia, diag, b, jb, ib, iw)

c*********************************************************************72
c
cc APLDIA adds a diagonal matrix to a general sparse matrix:  B = A + Diag.
c
c on entry:
c
c nrow      = integer. The row dimension of A
c
c job   = integer. job indicator. Job=0 means get array b only
c         (i.e. assume that a has already been copied into array b,
c         or that algorithm is used in place. ) For all practical
c         puposes enter job=0 for an in-place call and job=1 otherwise.
c
c         Note: in case there are missing diagonal elements in A,
c         then the option job =0 will be ignored, since the algorithm
c         must modify the data structure (i.e. jb, ib) in this
c         situation.
c
c a,
c ja,
c ia   = Matrix A in compressed sparse row format.
c
c diag = diagonal matrix stored as a vector dig(1:n)
c
c on return:
c
c
c b,
c jb,
c ib      = resulting matrix B in compressed sparse row sparse format.
c
c
c iw    = integer work array of length n. On return iw will
c         contain  the positions of the diagonal entries in the
c         output matrix. (i.e., a(iw(k)), ja(iw(k)), k=1,...n,
c         are the values/column indices of the diagonal elements
c         of the output matrix. ).
c
c Notes:
c
c 1)        The column dimension of A is not needed.
c 2)        algorithm in place (b, jb, ib, can be the same as
c           a, ja, ia, on entry). See comments for parameter job.
c
c coded by Y. Saad. Latest version July, 19, 1990
c
      double precision a(*), b(*), diag(nrow)
      integer ja(*),jb(*), ia(nrow+1),ib(nrow+1), iw(*)

      logical test
c
c     copy integer arrays into b's data structure if required
c
      if (job .ne. 0) then
         nnz = ia(nrow+1)-1
         do 2  k=1, nnz
            jb(k) = ja(k)
 2       continue
         do 3 k=1, nrow+1
            ib(k) = ia(k)
 3       continue
      end if
c
c     get positions of diagonal elements in data structure.
c
      call diapos (nrow,ja,ia,iw)
c
c     count number of holes in diagonal and add diag(*) elements to
c     valid diagonal entries.
c
      icount = 0
      do 1 j=1, nrow
         if (iw(j) .eq. 0) then
            icount = icount+1
         else
            b(iw(j)) = a(iw(j)) + diag(j)
         end if
 1    continue
c
c     if no diagonal elements to insert return
c
      if (icount .eq. 0) return
c
c     shift the nonzero elements if needed, to allow for created
c     diagonal elements.
c
      ko = ib(nrow+1)+icount
c
c     copy rows backward
c
      do 5 ii=nrow, 1, -1
c
c     go through  row ii
c
         k1 = ib(ii)
         k2 = ib(ii+1)-1
         ib(ii+1) = ko
         test = (iw(ii) .eq. 0)
         do 4 k = k2,k1,-1
            j = jb(k)
            if (test .and. (j .lt. ii)) then
               test = .false.
               ko = ko - 1
               b(ko) = diag(ii)
               jb(ko) = ii
               iw(ii) = ko
            end if
            ko = ko-1
            b(ko) = a(k)
            jb(ko) = j
 4       continue
c     diagonal element has not been added yet.
         if (test) then
            ko = ko-1
            b(ko) =  diag(ii)
            jb(ko) = ii
            iw(ii) = ko
         end if
 5    continue
      ib(1) = ko
      return
      end
