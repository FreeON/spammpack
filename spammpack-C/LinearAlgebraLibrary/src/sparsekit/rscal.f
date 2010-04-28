      subroutine rscal(nrow, job, nrm, a, ja, ia, diag, b, jb, ib)

c*********************************************************************72
c
cc RSCAL normalizes the rows of A.
c
c scales the rows of A such that their norms are one on return
c 3 choices of norms: 1-norm, 2-norm, max-norm.
c
c on entry:
c
c nrow      = integer. The row dimension of A
c
c job   = integer. job indicator. Job=0 means get array b only
c         job = 1 means get b, and the integer arrays ib, jb.
c
c nrm   = integer. norm indicator. nrm = 1, means 1-norm, nrm =2
c                  means the 2-nrm, nrm = 0 means max norm
c
c a,
c ja,
c ia   = Matrix A in compressed sparse row format.
c
c on return:
c
c
c diag = diagonal matrix stored as a vector containing the matrix
c        by which the rows have been scaled, i.e., on return
c        we have B = Diag*A.
c
c b,
c jb,
c ib      = resulting matrix B in compressed sparse row sparse format.
c
c Notes:
c
c 1)        The column dimension of A is not needed.
c 2)        algorithm in place (B can take the place of A).
c
      double precision a(*), b(*), diag(nrow)
      integer ja(*),jb(*), ia(nrow+1),ib(nrow+1)

      call rnrms (nrow,nrm,a,ja,ia,diag)

      do j=1, nrow
         diag(j) = 1.0/diag(j)
      end do

      call diamua(nrow,job,a,ja,ia,diag,b,jb,ib)
      return
      end
