      subroutine lnkcsr (n, a, jcol, istart, link, ao, jao, iao)

c*********************************************************************72
c
cc LNKCSR converts linked list storage to Compressed Sparse Row format.
c
c this subroutine translates a matrix stored in linked list storage
c format into the compressed sparse row format.
c
c Coded by Y. Saad, Feb 21, 1991.
c
c
c on entry:
c
c n      = integer equal to the dimension of A.
c
c a      = double precision array of size nna containing the nonzero elements
c jcol      = integer array of size      nnz containing the column positions
c         of the corresponding elements in a.
c istart= integer array of size n poiting to the beginning of the rows.
c         istart(i) contains the position of the first element of
c         row i in data structure. (a, jcol, link).
c         if a row is empty istart(i) must be zero.
c link      = integer array of size nnz containing the links in the linked
c         list data structure. link(k) points to the next element
c         of the row after element ao(k), jcol(k). if link(k) = 0,
c         then there is no next element, i.e., ao(k), jcol(k) is
c         the last element of the current row.
c
c on return:
c
c ao, jao, iao = matrix stored in csr format:
c
c ao    = double precision array containing the values of the nonzero elements of
c         the matrix stored row-wise.
c jao      = integer array of size nnz containing the column indices.
c iao      = integer array of size n+1 containing the pointers array to the
c         beginning of each row. iao(i) is the address in ao,jao of
c         first element of row i.
c
c  NZMAX is not provided in the calling sequence, so the following line
c  had to be commented out.
c
c     double precision a(*), ao(nzmax)
      double precision A(*), AO(*)
      integer n, jcol(*), istart(n), link(*), jao(*), iao(*)
      integer irow, ipos, next
c
c first determine individual bandwidths and pointers.
c
      ipos = 1
      iao(1) = ipos
c
c     loop through all rows
c
      do 100 irow =1, n
c
c     unroll i-th row.
c
         next = istart(irow)
 10      if (next .eq. 0) goto 99
         jao(ipos) = jcol(next)
         ao(ipos)  = a(next)
         ipos = ipos+1
         next = link(next)
         goto 10
 99      iao(irow+1) = ipos
 100  continue
c
      return
      end
