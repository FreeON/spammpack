      subroutine csrssk(n,imod,a,ja,ia,asky,isky,nzmax,ierr)

c*********************************************************************72
c
cc CSRSSK converts Compressed Sparse Row to Symmetric Skyline Format.
c
c CSRSSK translates a compressed sparse row or a symmetric
c sparse row format into a symmetric skyline format.
c the input matrix can be in either compressed sparse row or the
c symmetric sparse row format. The output matrix is in a symmetric
c skyline format: a double precision array containing the (active portions) of the
c rows in  sequence and a pointer to the beginning of each row.
c
c This module is NOT  in place.
c
c Coded by Y. Saad, Oct 5, 1989. Revised Feb. 18, 1991.
c
c on entry:
c
c n      = integer equal to the dimension of A.
c imod  = integer indicating the variant of skyline format wanted:
c         imod = 0 means the pointer isky points to the `zeroth'
c         element of the row, i.e., to the position of the diagonal
c         element of previous row (for i=1, isky(1)= 0)
c         imod = 1 means that itpr points to the beginning of the row.
c         imod = 2 means that isky points to the end of the row (diagonal
c                  element)
c
c a      = double precision array of size nna containing the nonzero elements
c ja      = integer array of size      nnz containing the column positions
c         of the corresponding elements in a.
c ia      = integer of size n+1. ia(k) contains the position in a, ja of
c        the beginning of the k-th row.
c nzmax = integer. must be set to the number of available locations
c         in the output array asky.
c
c on return:
c
c
c asky    = double precision array containing the values of the matrix stored in skyline
c         format. asky contains the sequence of active rows from
c         i=1, to n, an active row being the row of elemnts of
c         the matrix contained between the leftmost nonzero element
c         and the diagonal element.
c isky      = integer array of size n+1 containing the pointer array to
c         each row. The meaning of isky depends on the input value of
c         imod (see above).
c ierr  =  integer.  Error message. If the length of the
c         output array asky exceeds nzmax. ierr returns the minimum value
c         needed for nzmax. otherwise ierr=0 (normal return).
c
c Notes:
c         1) This module is NOT  in place.
c         2) even when imod = 2, length of  isky is  n+1, not n.
c
c
c first determine individial bandwidths and pointers.
c
      INTEGER NZMAX
      double precision a(*),asky(nzmax)
      integer n, imod, ierr, ia(n+1), isky(n+1), ja(*)
c
      ierr = 0
      isky(1) = 0
      do 3 i=1,n
         ml = 0
         do 31 k=ia(i),ia(i+1)-1
            ml = max(ml,i-ja(k)+1)
 31      continue
         isky(i+1) = isky(i)+ml
 3    continue
c
c     test if there is enough space  asky to do the copying.
c
      nnz = isky(n+1)
      if (nnz .gt. nzmax) then
         ierr = nnz
         return
      end if
c
c   fill asky with zeros.
c
      do 1 k=1, nnz
         asky(k) = 0.0
 1    continue
c
c     copy nonzero elements.
c
      do 4 i=1,n
         kend = isky(i+1)
         do 41 k=ia(i),ia(i+1)-1
            j = ja(k)
            if (j .le. i) asky(kend+j-i) = a(k)
 41      continue
 4    continue
c
c modify pointer according to imod if necessary.
c
      if (imod .eq. 0) return
      if (imod .eq. 1) then
         do 50 k=1, n+1
            isky(k) = isky(k)+1
 50      continue
      end if
      if (imod .eq. 2) then
         do 60 k=1, n
            isky(k) = isky(k+1)
 60      continue
      end if

      return
      end
