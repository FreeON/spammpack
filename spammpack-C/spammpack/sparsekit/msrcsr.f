      subroutine msrcsr (n,a,ja,ao,jao,iao,wk)

c*********************************************************************72
c
cc MSRCSR converts Modified Sparse Row to Compressed Sparse Row.
c
c converts a compressed matrix using a separated diagonal
c (modified sparse row format) in the Compressed Sparse Row
c format.
c does not check for zero elements in the diagonal.
c
c
c on entry :
c
c n        = row dimension of matrix
c ao, jao  = sparse matrix in msr sparse storage format
c           see routine csrmsr for details
c
c on return :
c
c a, ja, ia = matrix in csr format. note that the
c           algorithm is in place: ao, jao can be the same
c            as a, ja, in which case it will be overwritten on it
c            upon return.
c
c             here nnz = number of nonzero elements+1
c work arrays:
c
c wk      = ouble precision work array of length n
c
c notes:
c  In place algorithm (see a, ja, ia).
c

      double precision a(*),ao(*),wk(n)
      integer ja(*),jao(*),iao(n+1)

      logical added
      do 1 i=1,n
         wk(i) = a(i)
 1    continue
      iao(1) = 1
      iptr = 1

      do 500 ii=1,n
         added = .false.
         idiag = iptr + (ja(ii+1)-ja(ii))
         do 100 k=ja(ii),ja(ii+1)-1
            j = ja(k)
            if (j .lt. ii) then
               ao(iptr) = a(k)
               jao(iptr) = j
               iptr = iptr+1
            elseif (added) then
               ao(iptr) = a(k)
               jao(iptr) = j
               iptr = iptr+1
            else
c add diag element - only reserve a position for it.
               idiag = iptr
               iptr = iptr+1
               added = .true.
c     then other element
               ao(iptr) = a(k)
               jao(iptr) = j
               iptr = iptr+1
            end if
 100     continue
         ao(idiag) = wk(ii)
         jao(idiag) = ii
         if (.not. added) iptr = iptr+1
         iao(ii+1) = iptr
 500  continue
      return
      end
