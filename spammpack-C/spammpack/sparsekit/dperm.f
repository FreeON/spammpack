      subroutine dperm (nrow,a,ja,ia,ao,jao,iao,perm,qperm,job)

c*********************************************************************72
c
cc DPERM permutes the rows and columns of a matrix stored in CSR format.
c
c This routine permutes the rows and columns of a matrix stored in CSR
c format. i.e., it computes P A Q, where P, Q are permutation matrices.
c P maps row i into row perm(i) and Q maps column j into column qperm(j):
c      a(i,j)    becomes   a(perm(i),qperm(j)) in new matrix
c In the particular case where Q is the transpose of P (symmetric
c permutation of A) then qperm is not needed.
c note that qperm should be of length ncol (number of columns) but this
c is not checked.
c
c Y. Saad, Sep. 21 1989 / recoded Jan. 28 1991.
c
c on entry:
c
c n       = dimension of the matrix
c a, ja,
c    ia = input matrix in a, ja, ia format
c perm       = integer array of length n containing the permutation arrays
c        for the rows: perm(i) is the destination of row i in the
c         permuted matrix -- also the destination of column i in case
c         permutation is symmetric (job .le. 2)
c
c qperm      = same thing for the columns. This should be provided only
c         if job=3 or job=4, i.e., only in the case of a nonsymmetric
c        permutation of rows and columns. Otherwise qperm is a dummy
c
c job      = integer indicating the work to be done:
c * job = 1,2 permutation is symmetric  Ao :== P * A * transp(P)
c             job = 1      permute a, ja, ia into ao, jao, iao
c             job = 2 permute matrix ignoring real values.
c * job = 3,4 permutation is non-symmetric  Ao :== P * A * Q
c             job = 3      permute a, ja, ia into ao, jao, iao
c             job = 4 permute matrix ignoring real values.
c
c on return:
c
c ao, jao, iao = input matrix in a, ja, ia format
c
c in case job .eq. 2 or job .eq. 4, a and ao are never referred to
c and can be dummy arguments.
c Notes:
c
c  1) algorithm is in place
c  2) column indices may not be sorted on return even  though they may be
c     on entry.
c
      integer nrow,ja(*),ia(nrow+1),jao(*),iao(nrow+1),perm(nrow),
     +        qperm(*),job
      double precision a(*),ao(*)
      integer locjob, mod
c
c     locjob indicates whether or not real values must be copied.
c
      locjob = mod(job,2)
c
c permute rows first
c
      call rperm (nrow,a,ja,ia,ao,jao,iao,perm,locjob)
c
c then permute columns
c
      locjob = 0

      if (job .le. 2) then
         call cperm (nrow,ao,jao,iao,ao,jao,iao,perm,locjob)
      else
         call cperm (nrow,ao,jao,iao,ao,jao,iao,qperm,locjob)
      end if

      return
      end
