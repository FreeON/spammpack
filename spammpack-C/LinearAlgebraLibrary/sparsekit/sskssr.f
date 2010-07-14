      subroutine sskssr (n,imod,asky,isky,ao,jao,iao,nzmax,ierr)

c*********************************************************************72
c
cc SSKSSR converts Symmetric Skyline Format to Symmetric Sparse Row format.
c
c  tests for exact zeros in skyline matrix (and ignores them in
c  output matrix).  In place routine (a, isky :: ao, iao)
c
c this subroutine translates a  symmetric skyline format into a
c symmetric sparse row format. Each element is tested to see if it is
c a zero element. Only the actual nonzero elements are retained. Note
c that the test used is simple and does take into account the smallness
c of a value. the subroutine filter (see unary module) can be used
c for this purpose.
c
c Coded by Y. Saad, Oct 5, 1989. Revised Feb 18, 1991./
c
c on entry:
c
c n      = integer equal to the dimension of A.
c imod  = integer indicating the variant of skyline format used:
c         imod = 0 means the pointer iao points to the `zeroth'
c         element of the row, i.e., to the position of the diagonal
c         element of previous row (for i=1, iao(1)= 0)
c         imod = 1 means that itpr points to the beginning of the row.
c         imod = 2 means that iao points to the end of the row
c                  (diagonal element)
c asky  = double precision array containing the values of the matrix. asky contains
c         the sequence of active rows from i=1, to n, an active row
c         being the row of elemnts of the matrix contained between the
c         leftmost nonzero element and the diagonal element.
c isky       = integer array of size n+1 containing the pointer array to
c         each row. isky (k) contains the address of the beginning of the
c         k-th active row in the array asky.
c nzmax = integer. equal to the number of available locations in the
c         output array ao.
c
c on return:
c
c ao      = double precision array of size nna containing the nonzero elements
c jao      = integer array of size      nnz containing the column positions
c         of the corresponding elements in a.
c iao      = integer of size n+1. iao(k) contains the position in a, ja of
c        the beginning of the k-th row.
c ierr  = integer. Serving as error message. If the length of the
c         output arrays ao, jao exceeds nzmax then ierr returns
c         the row number where the algorithm stopped: rows
c         i, to ierr-1 have been processed succesfully.
c         ierr = 0 means normal return.
c         ierr = -1  : illegal value for imod
c Notes:
c
c This module is in place: ao and iao can be the same as asky, and isky.
c
      INTEGER NZMAX
      double precision asky(*),ao(nzmax)
      integer n, imod,ierr, isky(n+1),iao(n+1),jao(nzmax)
      integer next, kend, kstart, i, j
      ierr = 0
c
c check for validity of imod
c
      if (imod.ne.0 .and. imod.ne.1 .and. imod .ne. 2) then
         ierr =-1
         return
      end if
c
c next  = pointer to next available position in output matrix
c kend  = pointer to end of current row in skyline matrix.
c
      next = 1
c
c set kend = start position -1 in  skyline matrix.
c
      kend = 0
      if (imod .eq. 1) kend = isky(1)-1
      if (imod .eq. 0) kend = isky(1)
c
c loop through all rows
c
      do 50 i=1,n
c
c save value of pointer to ith row in output matrix
c
         iao(i) = next
c
c get beginnning and end of skyline  row
c
         kstart = kend+1
         if (imod .eq. 0) kend = isky(i+1)
         if (imod .eq. 1) kend = isky(i+1)-1
         if (imod .eq. 2) kend = isky(i)
c
c copy element into output matrix unless it is a zero element.
c
         do 40 k=kstart,kend
            if (asky(k) .eq. 0.0) goto 40
            j = i-(kend-k)
            jao(next) = j
            ao(next)  = asky(k)
            next=next+1
            if (next .gt. nzmax+1) then
               ierr = i
               return
            end if
 40      continue
 50    continue
      iao(n+1) = next
      return
      end
