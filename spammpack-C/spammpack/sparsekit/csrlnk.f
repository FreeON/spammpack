      subroutine csrlnk (n,a,ja,ia,link)

c*********************************************************************72
c
cc CSRLNK converts Compressed Sparse Row to Linked storage format.
c
c CSRLNK translates a matrix stored in compressed sparse
c row into one with a linked list storage format. Only the link
c array needs to be obtained since the arrays a, ja, and ia may
c be unchanged and have carry the same meaning for the output matrix.
c in  other words a, ja, ia, link   ia the output linked list data
c structure with a, ja, ia being the same.
c
c Coded by Y. Saad, Feb 21, 1991.
c
c
c on entry:
c
c n      = integer equal to the dimension of A.
c
c a      = double precision array of size nna containing the nonzero elements
c ja      = integer array of size      nnz containing the column positions
c         of the corresponding elements in a.
c ia      = integer of size n+1 containing the pointers to the beginning
c         of each row. ia(k) contains the position in a, ja of the
c         beginning of the k-th row.
c
c on return:
c
c a, ja, ia are not changed.
c
c a     = nonzero elements.
c ja    = column positions.
c ia    = points to the first row of matrix in structure.
c link      = integer array of size containing the linked list information.
c         link(k) points to the next element of the row after element
c         ao(k), jcol(k). if link(k) = 0, then there is no next element,
c         i.e., ao(k), jcol(k) is the last element of the current row.
c
      double precision a(*)
      integer n, ja(*), ia(n+1), link(*)
      integer i, k
c
c loop through all rows
c
      do 100 i =1, n
         do 99  k=ia(i), ia(i+1)-2
            link(k) = k+1
 99      continue
         link(ia(i+1)-1) = 0
 100  continue

      return
      end
