      subroutine dnscsr_single(nrow,ncol,nzmax,dns,ndns,a,ja,ia,ierr)

c*********************************************************************72
c
cc DNSCSR converts Dense to Compressed Row Sparse format.
c
c converts a densely stored matrix into a row orientied
c compactly sparse matrix. ( reverse of csrdns )
c Note: this routine does not check whether an element
c is small. It considers that a(i,j) is zero if it is exactly
c equal to zero: see test below.
c
c on entry:
c
c
c nrow      = row-dimension of a
c ncol      = column dimension of a
c nzmax = maximum number of nonzero elements allowed. This
c         should be set to be the lengths of the arrays a and ja.
c dns   = input nrow x ncol (dense) matrix.
c ndns      = first dimension of dns.
c
c on return:
c
c
c a, ja, ia = value, column, pointer  arrays for output matrix
c
c ierr      = integer error indicator:
c         ierr .eq. 0 means normal retur
c         ierr .eq. i means that the the code stopped while
c         processing row number i, because there was no space left in
c         a, and ja (as defined by parameter nzmax).
c
      implicit none
      real dns(ndns,*), a(*)
      integer nrow, ncol, nzmax, ndns, ierr, ia(*),ja(*)
      integer i, j, next

      ierr = 0
      next = 1
      ia(1) = 1
      do 4 i=1,nrow
         do 3 j=1, ncol
            if (dns(i,j) .eq. 0.0) goto 3
            if (next .gt. nzmax) then
               ierr = i
               return
            end if
            ja(next) = j
            a(next) = dns(i,j)
            next = next+1
 3       continue
         ia(i+1) = next
 4    continue
      return
      end
