      subroutine aplb (nrow,ncol,job,a,ja,ia,b,jb,ib,
     *     c,jc,ic,nzmax,iw,ierr)

c*********************************************************************72
c
cc APLB performs the CSR matrix sum C = A + B.
c
c on entry:
c
c nrow      = integer. The row dimension of A and B
c ncol  = integer. The column dimension of A and B.
c job   = integer. Job indicator. When job = 0, only the structure
c                  (i.e. the arrays jc, ic) is computed and the
c                  real values are ignored.
c
c a,
c ja,
c ia   = Matrix A in compressed sparse row format.
c
c b,
c jb,
c ib      =  Matrix B in compressed sparse row format.
c
c nzmax      = integer. The  length of the arrays c and jc.
c         amub will stop if the result matrix C  has a number
c         of elements that exceeds exceeds nzmax. See ierr.
c
c on return:
c
c c,
c jc,
c ic      = resulting matrix C in compressed sparse row sparse format.
c
c ierr      = integer. serving as error message.
c         ierr = 0 means normal return,
c         ierr .gt. 0 means that amub stopped while computing the
c         i-th row  of C with i=ierr, because the number
c         of elements in C exceeds nzmax.
c
c work arrays:
c
c iw      = integer work array of length equal to the number of
c         columns in A.
c
      double precision a(*), b(*), c(*)
      integer ja(*),jb(*),jc(*),ia(nrow+1),ib(nrow+1),ic(nrow+1),
     *     iw(ncol)
      logical values

      values = (job .ne. 0)
      ierr = 0
      len = 0
      ic(1) = 1
      do 1 j=1, ncol
         iw(j) = 0
 1    continue

      do 500 ii=1, nrow
c     row i
         do 200 ka=ia(ii), ia(ii+1)-1
            len = len+1
            jcol    = ja(ka)
            if (len .gt. nzmax) goto 999
            jc(len) = jcol
            if (values) c(len)  = a(ka)
            iw(jcol)= len
 200     continue

         do 300 kb=ib(ii),ib(ii+1)-1
            jcol = jb(kb)
            jpos = iw(jcol)
            if (jpos .eq. 0) then
               len = len+1
               if (len .gt. nzmax) goto 999
               jc(len) = jcol
               if (values) c(len)  = b(kb)
               iw(jcol)= len
            else
               if (values) c(jpos) = c(jpos) + b(kb)
            end if
 300     continue
         do 301 k=ic(ii), len
          iw(jc(k)) = 0
 301     continue
         ic(ii+1) = len+1
 500  continue
      return
 999  ierr = ii
      return
      end
