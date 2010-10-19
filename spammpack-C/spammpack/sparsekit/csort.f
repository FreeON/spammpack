      subroutine csort (n,a,ja,ia,iwork,values)

c*********************************************************************72
c
cc CSORT sorts the elements of a CSR matrix.
c
c This routine sorts the elements of  a matrix (stored in Compressed
c Sparse Row Format) in increasing order of their column indices within
c each row. It uses a form of bucket sort with a cost of O(nnz) where
c nnz = number of nonzero elements.
c requires an integer work array of size length 2*nnz.
c
c on entry:
c
c n     = the row dimension of the matrix
c a     = the matrix A in compressed sparse row format.
c ja    = the array of column indices of the elements in array a.
c ia    = the array of pointers to the rows.
c iwork = integer work array of length max ( n+1, 2*nnz )
c         where nnz = 2* (ia(n+1)-ia(1))  ) .
c values= logical indicating whether or not the real values a(*) must
c         also be permuted. if (.not. values) then the array a is not
c         touched by csort and can be a dummy array.
c
c on return:
c
c the matrix stored in the structure a, ja, ia is permuted in such a
c way that the column indices are in increasing order within each row.
c iwork(1:nnz) contains the permutation used  to rearrange the elements.
c
c Y. Saad - Feb. 1, 1991.
c
      logical values
      integer n, ja(*), ia(n+1), iwork(*)
      double precision a(*)
      integer i, k, j, ifirst, nnz, next
c
c count the number of elements in each column
c
      do 1 i=1,n+1
         iwork(i) = 0
 1    continue
      do 3 i=1, n
         do 2 k=ia(i), ia(i+1)-1
            j = ja(k)+1
            iwork(j) = iwork(j)+1
 2       continue
 3    continue
c
c compute pointers from lengths.
c
      iwork(1) = 1
      do 4 i=1,n
         iwork(i+1) = iwork(i) + iwork(i+1)
 4    continue
c
c get the positions of the nonzero elements in order of columns.
c
      ifirst = ia(1)
      nnz = ia(n+1)-ifirst
      do 5 i=1,n
         do 51 k=ia(i),ia(i+1)-1
            j = ja(k)
            next = iwork(j)
            iwork(nnz+next) = k
            iwork(j) = next+1
 51      continue
 5    continue
c
c convert to coordinate format
c
      do 6 i=1, n
         do 61 k=ia(i), ia(i+1)-1
            iwork(k) = i
 61      continue
 6    continue
c
c loop to find permutation: for each element find the correct
c position in (sorted) arrays a, ja. Record this in iwork.
c
      do 7 k=1, nnz
         ko = iwork(nnz+k)
         irow = iwork(ko)
         next = ia(irow)
c
c the current element should go in next position in row. iwork
c records this position.
c
         iwork(ko) = next
         ia(irow)  = next+1
 7       continue
c
c perform an in-place permutation of the  arrays.
c
         call ivperm (nnz, ja(ifirst), iwork)
         if (values) call dvperm (nnz, a(ifirst), iwork)
c
c reshift the pointers of the original matrix back.
c
      do 8 i=n,1,-1
         ia(i+1) = ia(i)
 8    continue
      ia(1) = ifirst
c
      return
      end
