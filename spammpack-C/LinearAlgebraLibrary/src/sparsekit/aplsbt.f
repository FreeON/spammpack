      subroutine aplsbt(nrow,ncol,a,ja,ia,s,b,jb,ib,
     *     c,jc,ic,nzmax,iw,ierr)

c*********************************************************************72
c
cc APLSBT performs the matrix sum C = A + B'.
c
c on entry:
c
c nrow      = integer. The row dimension of A and transp(B)
c ncol  = integer. The column dimension of A. Also the row
c                  dimension of B.
c
c a,
c ja,
c ia    = Matrix A in compressed sparse row format.
c
c s      = real. scalar factor for B.
c
c
c b,
c jb,
c ib      =  Matrix B in compressed sparse row format.
c
c nzmax      = integer. The  length of the arrays c, jc, and ic.
c         amub will stop if the result matrix C  has a number
c         of elements that exceeds exceeds nzmax. See ierr.
c
c on return:
c
c c,
c jc,
c ic      = resulting matrix C in compressed sparse row format.
c
c ierr      = integer. serving as error message.
c         ierr = 0 means normal return.
c         ierr = -1 means that nzmax was .lt. either the number of
c         nonzero elements of A or the number of nonzero elements in B.
c         ierr .gt. 0 means that amub stopped while computing the
c         i-th row  of C with i=ierr, because the number
c         of elements in C exceeds nzmax.
c
c work arrays:
c
c iw      = integer work array of length equal to the number of
c         columns in A.
c
c Notes:
c        It is important to note that here all of three arrays c, ic,
c        and jc are assumed to be of length nnz(c). This is because
c        the matrix is internally converted in coordinate format.
c
      double precision a(*), b(*), c(*), s
      integer ja(*),jb(*),jc(*),ia(nrow+1),ib(ncol+1),ic(*),
     *     iw(ncol)

      ierr = 0
      do 1 j=1, ncol
         iw(j) = 0
 1    continue

      nnza = ia(nrow+1)-1
      nnzb = ib(ncol+1)-1
      len = nnzb
      if (nzmax .lt. nnzb .or. nzmax .lt. nnza) then
         ierr = -1
         return
      end if
c
c     trasnpose matrix b into c
c
      ljob = 1
      ipos = 1
      call csrcsc (ncol,ljob,ipos,b,jb,ib,c,jc,ic)
      do 2 k=1,len
 2       c(k) = c(k)*s
c
c     main loop. add rows from ii = 1 to nrow.
c
         do 500 ii=1, nrow
c     iw is used as a system to recognize whether there
c     was a nonzero element in c.
            do 200 k = ic(ii),ic(ii+1)-1
               iw(jc(k)) = k
 200        continue
c
            do 300 ka = ia(ii), ia(ii+1)-1
               jcol = ja(ka)
               jpos = iw(jcol)
           if (jpos .eq. 0) then
c
c     if fill-in append in coordinate format to matrix.
c
              len = len+1
              if (len .gt. nzmax) goto 999
              jc(len) = jcol
              ic(len) = ii
              c(len)  = a(ka)
           else
c     else do addition.
              c(jpos) = c(jpos) + a(ka)
           end if
 300    continue
        do 301 k=ic(ii), ic(ii+1)-1
           iw(jc(k)) = 0
 301    continue
 500  continue
c
c     convert matrix without fill-ins into coordinate format
c
      ljob = 3
      call csrcoo (nrow,ljob,nnzb,c,jc,ic,nnzb,c,ic,jc,ierr)
      if (ierr .ne. 0) ierr = -ierr
c     convert the whole thing back to csr format.
      ljob = 1
      call coicsr (nrow,len,1,c,jc,ic,iw)
      return
 999  ierr = ii
      return
      end
