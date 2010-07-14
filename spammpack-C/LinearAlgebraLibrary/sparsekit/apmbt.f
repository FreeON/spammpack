      subroutine apmbt (nrow,ncol,job,a,ja,ia,b,jb,ib,
     *     c,jc,ic,nzmax,iw,ierr)

c*********************************************************************72
c
cc APMBT performs the matrix sum  C = A + B' or C = A - B'.
c
c on entry:
c
c nrow      = integer. The row dimension of A and transp(B)
c ncol  = integer. The column dimension of A. Also the row
c                  dimension of B.
c
c job      = integer. if job = -1, apmbt will compute C= A - transp(B)
c         (structure + values)
c         if (job .eq. 1)  it will compute C=A+transp(A)
c         (structure+ values)
c         if (job .eq. 0) it will compute the structure of
c         C= A+/-transp(B) only (ignoring all real values).
c         any other value of job will be treated as  job=1
c a,
c ja,
c ia    = Matrix A in compressed sparse row format.
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
c  It is important to note that here all of three arrays c, ic,
c        and jc are assumed to be of length nnz(c). This is because
c        the matrix is internally converted in coordinate format.
c
      double precision a(*), b(*), c(*)
      integer ja(*),jb(*),jc(*),ia(nrow+1),ib(ncol+1),ic(*),
     *     iw(ncol)

      logical values
      values = (job .ne. 0)

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
c trasnpose matrix b into c
c
      ljob = 0
      if (values) ljob = 1
      ipos = 1
      call csrcsc (ncol,ljob,ipos,b,jb,ib,c,jc,ic)

      if (job .eq. -1) then
         do 2 k=1,len
          c(k) = -c(k)
 2       continue
      end if
c
c  main loop
c
      do 500 ii=1, nrow
         do 200 k = ic(ii),ic(ii+1)-1
            iw(jc(k)) = k
 200     continue

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
               if (values) c(len)  = a(ka)
            else
c     else do addition.
               if (values) c(jpos) = c(jpos) + a(ka)
            end if
 300     continue
         do 301 k=ic(ii), ic(ii+1)-1
          iw(jc(k)) = 0
 301     continue
 500  continue
c
c     convert first part of matrix (without fill-ins) into coo format
c
      ljob = 2
      if (values) ljob = 3
      call csrcoo (nrow,ljob,nnzb,c,jc,ic,nnzb,c,ic,jc,ierr)
      if (ierr .ne. 0) ierr = -ierr
c     convert the whole thing back to csr format.
      ljob = 0
      if (values) ljob = 1
      call coicsr (nrow,len,ljob,c,jc,ic,iw)
      return
 999  ierr = ii
      return
      end
