      subroutine amub (nrow,ncol,job,a,ja,ia,b,jb,ib,
     *                 c,jc,ic,nzmax,iw,ierr)

c*********************************************************************72
c
cc AMUB performs the matrix product C = A * B.
c
c on entry:
c
c nrow  = integer. The row dimension of A
c ncol  = integer. The column dimension of A
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
c ib    =  Matrix B in compressed sparse row format.
c
c nzmax = integer. The  length of the arrays c and jc.
c         amub will stop if the result matrix C  has a number
c         of elements that exceeds exceeds nzmax. See ierr.
c
c on return:
c
c c,
c jc,
c ic    = resulting matrix C in compressed sparse row sparse format.
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
c Notes:
c
c         The column dimension of B is not needed.
c
      double precision a(*), b(*), c(*)
      integer ja(*),jb(*),jc(*),ia(nrow+1),ib(ncol+1),
     *     ic(ncol+1),iw(ncol)

      double precision scal
      logical values
      values = (job .ne. 0)
      len = 0
      ic(1) = 1
      ierr = 0
c
c  Initialize array iw.
c
      do 1 j=1, ncol
         iw(j) = 0
 1    continue

      do 500 ii=1, nrow
c     row i
         do 200 ka=ia(ii), ia(ii+1)-1
          if (values) scal = a(ka)
          jj   = ja(ka)
          do 100 kb=ib(jj),ib(jj+1)-1
               jcol = jb(kb)
               jpos = iw(jcol)
               if (jpos .eq. 0) then
                  len = len+1
                  if (len .gt. nzmax) then
                     ierr = ii
                     return
                  end if
                  jc(len) = jcol
                  iw(jcol)= len
                  if (values) c(len)  = scal*b(kb)
               else
                  if (values) c(jpos) = c(jpos) + scal*b(kb)
               end if
 100          continue
 200     continue
         do 201 k=ic(ii), len
          iw(jc(k)) = 0
 201     continue
         ic(ii+1) = len+1
 500  continue
      return
      end
