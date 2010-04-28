      subroutine aplsb1 (nrow,ncol,a,ja,ia,s,b,jb,ib,c,jc,ic,
     *     nzmax,ierr)

c*********************************************************************72
c
cc APLSB1 performs the operation C = A + s * B for sorted CSR matrices.
c
c the difference with aplsb is that the resulting matrix is such that
c the elements of each row are sorted with increasing column indices in
c each row, provided the original matrices are sorted in the same way.
c
c on entry:
c
c nrow      = integer. The row dimension of A and B
c ncol  = integer. The column dimension of A and B.
c
c a,
c ja,
c ia   = Matrix A in compressed sparse row format with entries sorted
c
c s      = real. scalar factor for B.
c
c b,
c jb,
c ib      =  Matrix B in compressed sparse row format with entries sorted
c        ascendly in each row
c
c nzmax      = integer. The  length of the arrays c and jc.
c         amub will stop if the result matrix C  has a number
c         of elements that exceeds exceeds nzmax. See ierr.
c
c on return:
c
c c,
c jc,
c ic      = resulting matrix C in compressed sparse row sparse format
c         with entries sorted ascendly in each row.
c
c ierr      = integer. serving as error message.
c         ierr = 0 means normal return,
c         ierr .gt. 0 means that amub stopped while computing the
c         i-th row  of C with i=ierr, because the number
c         of elements in C exceeds nzmax.
c
c Notes:
c
c     this will not work if any of the two input matrices is not sorted
c
      double precision a(*), b(*), c(*), s
      integer ja(*),jb(*),jc(*),ia(nrow+1),ib(nrow+1),ic(nrow+1)

      ierr = 0
      kc = 1
      ic(1) = kc

      do 6 i=1, nrow
         ka = ia(i)
         kb = ib(i)
         kamax = ia(i+1)-1
         kbmax = ib(i+1)-1
 5       continue
         if (ka .le. kamax) then
            j1 = ja(ka)
         else
            j1 = ncol+1
         end if
         if (kb .le. kbmax) then
            j2 = jb(kb)
         else
            j2 = ncol+1
         end if
c
c     three cases
c
         if (j1 .eq. j2) then
            c(kc) = a(ka)+s*b(kb)
            jc(kc) = j1
            ka = ka+1
            kb = kb+1
            kc = kc+1
         else if (j1 .lt. j2) then
            jc(kc) = j1
            c(kc) = a(ka)
            ka = ka+1
            kc = kc+1
         else if (j1 .gt. j2) then
            jc(kc) = j2
            c(kc) = s*b(kb)
            kb = kb+1
            kc = kc+1
         end if
         if (kc .gt. nzmax) goto 999
         if (ka .le. kamax .or. kb .le. kbmax) goto 5
         ic(i+1) = kc
 6    continue
      return
 999  ierr = i
      return
      end
