      subroutine csrcoo (nrow,job,nzmax,a,ja,ia,nnz,ao,ir,jc,ierr)

c*********************************************************************72
c
cc CSRCOO converts Compressed Sparse Row to Coordinate format.
c
c converts a matrix that is stored in coordinate format
c  a, ir, jc into a row general sparse ao, jao, iao format.
c
c on entry:
c
c nrow      = dimension of the matrix.
c job   = integer serving as a job indicator.
c         if job = 1 fill in only the array ir, ignore jc, and ao.
c         if job = 2 fill in ir, and jc but not ao
c         if job = 3 fill in everything.
c         The reason why these options are provided is that on return
c         ao and jc are the same as a, ja. So when job = 3, a and ja are
c         simply copied into ao, jc.  When job=2, only jc and ir are
c         returned. With job=1 only the array ir is returned. Moreover,
c         the algorithm is in place:
c           call csrcoo (nrow,1,nzmax,a,ja,ia,nnz,a,ia,ja,ierr)
c         will write the output matrix in coordinate format on a, ja,ia.
c         (Important: note the order in the output arrays a, ja, ia. )
c         i.e., ao can be the same as a, ir can be the same as ia
c         and jc can be the same as ja.
c
c a,
c ja,
c ia    = matrix in compressed sparse row format.
c nzmax = length of space available in ao, ir, jc.
c         the code will stop immediatly if the number of
c         nonzero elements found in input matrix exceeds nzmax.
c
c on return:
c
c ao, ir, jc = matrix in coordinate format.
c
c nnz        = number of nonzero elements in matrix.
c ierr       = integer error indicator.
c         ierr .eq. 0 means normal retur
c         ierr .eq. 1 means that the the code stopped
c         because there was no space in ao, ir, jc
c         (according to the value of  nzmax).
c
      double precision a(*),ao(*)
      integer ir(*),jc(*),ja(*),ia(*)

      ierr = 0
      nnz = ia(nrow+1)-1
      if (nnz .gt. nzmax) then
         ierr = 1
         return
      end if

      goto (3,2,1) job
 1    do 10 k=1,nnz
         ao(k) = a(k)
 10   continue
 2    do 11 k=1,nnz
         jc(k) = ja(k)
 11   continue
c copy backward to allow
 3    do 13 i=nrow,1,-1
         k1 = ia(i+1)-1
         k2 = ia(i)
         do 12 k=k1,k2,-1
            ir(k) = i
 12      continue
 13   continue
      return
      end
