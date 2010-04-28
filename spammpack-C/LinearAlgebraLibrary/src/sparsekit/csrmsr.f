      subroutine csrmsr (n,a,ja,ia,ao,jao,wk,iwk)

c*********************************************************************72
c
cc CSRMSR converts Compressed Sparse Row to Modified Sparse Row.
c
c converts a general sparse matrix a, ja, ia into
c a compressed matrix using a separated diagonal (referred to as
c the bell-labs format as it is used by bell labs semi conductor
c group. We refer to it here as the modified sparse row format.
c Note: this has been coded in such a way that one can overwrite
c the output matrix onto the input matrix if desired by a call of
c the form
c
c     call csrmsr (n, a, ja, ia, a, ja, wk,iwk)
c
c In case ao, jao, are different from a, ja, then one can
c use ao, jao as the work arrays in the calling sequence:
c
c     call csrmsr (n, a, ja, ia, ao, jao, ao,jao)
c
c
c
c on entry :
c
c a, ja, ia = matrix in csr format. note that the
c           algorithm is in place: ao, jao can be the same
c            as a, ja, in which case it will be overwritten on it
c            upon return.
c
c on return :
c
c ao, jao  = sparse matrix in modified sparse row storage format:
c         +  ao(1:n) contains the diagonal of the matrix.
c         +  ao(n+2:nnz) contains the nondiagonal elements of the
c             matrix, stored rowwise.
c         +  jao(n+2:nnz) : their column indices
c         +  jao(1:n+1) contains the pointer array for the nondiagonal
c             elements in ao(n+1:nnz) and jao(n+2:nnz).
c             i.e., for i .le. n+1 jao(i) points to beginning of row i
c            in arrays ao, jao.
c             here nnz = number of nonzero elements+1
c work arrays:
c
c wk      = real work array of length n
c iwk   = integer work array of length n+1
c
c notes:
c
c        Algorithm is in place.  i.e. both:
c
c          call csrmsr (n, a, ja, ia, ao, jao, ao,jao)
c          (in which  ao, jao, are different from a, ja)
c           and
c          call csrmsr (n, a, ja, ia, a, ja, wk,iwk)
c          (in which  wk, jwk, are different from a, ja)
c        are OK.
c
c coded by Y. Saad Sep. 1989. Rechecked Feb 27, 1990.
c
      double precision a(*),ao(*),wk(n)
      integer ia(n+1),ja(*),jao(*),iwk(n+1)

      icount = 0
c
c store away diagonal elements and count nonzero diagonal elements.
c
      do 1 i=1,n
         wk(i) = 0.0
         iwk(i+1) = ia(i+1)-ia(i)
         do 2 k=ia(i),ia(i+1)-1
            if (ja(k) .eq. i) then
               wk(i) = a(k)
               icount = icount + 1
               iwk(i+1) = iwk(i+1)-1
            end if
 2       continue
 1    continue
c
c compute total length
c
      iptr = n + ia(n+1) - icount
c
c     copy backwards (to avoid collisions)
c
      do 500 ii=n,1,-1
         do 100 k=ia(ii+1)-1,ia(ii),-1
            j = ja(k)
            if (j .ne. ii) then
               ao(iptr) = a(k)
               jao(iptr) = j
               iptr = iptr-1
            end if
 100     continue
 500  continue
c
c compute pointer values and copy wk(*)
c
      jao(1) = n+2
      do 600 i=1,n
         ao(i) = wk(i)
         jao(i+1) = jao(i)+iwk(i+1)
 600  continue
      return
      end
