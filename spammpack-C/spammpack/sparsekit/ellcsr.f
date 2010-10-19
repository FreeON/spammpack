      subroutine ellcsr(nrow,coef,jcoef,ncoef,ndiag,a,ja,ia,nzmax,ierr)

c*********************************************************************72
c
cc ELLCSR converts Ellpack/Itpack to Compressed Sparse Row.
c
c this subroutine converts a matrix stored in ellpack-itpack format
c coef-jcoef into the compressed sparse row format. It actually checks
c whether an entry in the input matrix is a nonzero element before
c putting it in the output matrix. The test does not account for small
c values but only for exact zeros.
c
c on entry:
c
c
c nrow       = row dimension of the matrix A.
c coef      = array containing the values of the matrix A in ellpack format.
c jcoef = integer arraycontains the column indices of coef(i,j) in A.
c ncoef = first dimension of arrays coef, and jcoef.
c ndiag = number of active columns in coef, jcoef.
c
c ndiag = on entry the number of columns made available in coef.
c
c on return:
c
c a, ia,
c    ja = matrix in a, ia, ja format where.
c
c nzmax      = size of arrays a and ja. ellcsr will abort if the storage
c         provided in a, ja is not sufficient to store A. See ierr.
c
c ierr       = integer. serves are output error message.
c         ierr = 0 means normal return.
c         ierr = 1 means that there is not enough space in
c         a and ja to store output matrix.
c
      integer ia(nrow+1), ja(*), jcoef(ncoef,1)
      double precision a(*), coef(ncoef,1)
c
c first determine the length of each row of lower-part-of(A)
      ierr = 0
c  check whether sufficient columns are available.
c
c copy elements row by row.
      kpos = 1
      do 6 i=1, nrow
         do 5 k=1,ndiag
            if (coef(i,k) .ne. 0.0) then
               if (kpos .gt. nzmax) then
                  ierr = kpos
                  return
               end if
               a(kpos) = coef(i,k)
               ja(kpos) = jcoef(i,k)
               kpos = kpos+1
          end if
 5       continue
         ia(i+1) = kpos
 6    continue
      return
      end
