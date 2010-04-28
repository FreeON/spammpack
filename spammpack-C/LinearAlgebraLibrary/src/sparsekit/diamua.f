      subroutine diamua (nrow,job, a, ja, ia, diag, b, jb, ib)

c*********************************************************************72
c
cc DIAMUA performs the matrix by matrix product B = Diag * A.
c
c on entry:
c
c nrow      = integer. The row dimension of A
c
c job   = integer. job indicator. Job=0 means get array b only
c         job = 1 means get b, and the integer arrays ib, jb.
c
c a,
c ja,
c ia   = Matrix A in compressed sparse row format.
c
c diag = diagonal matrix stored as a vector dig(1:n)
c
c on return:
c
c
c b,
c jb,
c ib      = resulting matrix B in compressed sparse row sparse format.
c
c Notes:
c
c 1)        The column dimension of A is not needed.
c 2)        algorithm in place (B can take the place of A).
c           in this case use job=0.
c
      double precision a(*), b(*), diag(nrow), scal
      integer ja(*),jb(*), ia(nrow+1),ib(nrow+1)

      do 1 ii=1,nrow
c
c     normalize each row
c
         k1 = ia(ii)
         k2 = ia(ii+1)-1
         scal = diag(ii)
         do 2 k=k1, k2
            b(k) = a(k)*scal
 2       continue
 1    continue

      if (job .eq. 0) return

      ib(1) = ia(1)
      do 3 ii=1, nrow
         ib(ii) = ia(ii)
         do 31 k=ia(ii),ia(ii+1)-1
            jb(k) = ja(k)
 31      continue
 3    continue
      return
      end
