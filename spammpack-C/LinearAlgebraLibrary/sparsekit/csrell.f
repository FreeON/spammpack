      subroutine csrell (nrow,a,ja,ia,maxcol,coef,jcoef,ncoef,
     *                   ndiag,ierr)

c*********************************************************************72
c
cc CSRELL converts Compressed Sparse Row to Ellpack/Itpack format.
c
c on entry:
c
c nrow         = row dimension of the matrix A.
c
c a,
c ia,
c ja      = input matrix in compressed sparse row format.
c
c ncoef  = first dimension of arrays coef, and jcoef.
c
c maxcol = integer equal to the number of columns available in coef.
c
c on return:
c
c coef      = double precision array containing the values of the matrix A in
c         itpack-ellpack format.
c jcoef = integer array containing the column indices of coef(i,j)
c         in A.
c ndiag = number of active 'diagonals' found.
c
c ierr       = error message. 0 = correct return. If ierr .ne. 0 on
c        return this means that the number of diagonals found
c         (ndiag) exceeds maxcol.
c
      integer ia(nrow+1), ja(*), jcoef(ncoef,1)
      double precision a(*), coef(ncoef,1)
c
c first determine the length of each row of lower-part-of(A)
      ierr = 0
      ndiag = 0
      do 3 i=1, nrow
         k = ia(i+1)-ia(i)
         ndiag = max0(ndiag,k)
 3    continue
c  check whether sufficient columns are available.
      if (ndiag .gt. maxcol) then
         ierr = 1
         return
      end if
c
c fill coef with zero elements and jcoef with row numbers.
c
      do 4 j=1,ndiag
         do 41 i=1,nrow
            coef(i,j) = 0.0
            jcoef(i,j) = i
 41      continue
 4    continue
c
c  copy elements row by row.
c
      do 6 i=1, nrow
         k1 = ia(i)
         k2 = ia(i+1)-1
         do 5 k=k1,k2
            coef(i,k-k1+1) = a(k)
            jcoef(i,k-k1+1) = ja(k)
 5       continue
 6    continue
      return
      end
