      subroutine aplsca (nrow, a, ja, ia, scal,iw)

c*********************************************************************72
c
cc APLSCA adds a scalar to the diagonal entries of a sparse matrix A :=A + s I
c
c on entry:
c
c nrow      = integer. The row dimension of A
c
c a,
c ja,
c ia    = Matrix A in compressed sparse row format.
c
c scal  = real. scalar to add to the diagonal entries.
c
c on return:
c
c
c a,
c ja,
c ia      = matrix A with diagonal elements shifted (or created).
c
c iw    = integer work array of length n. On return iw will
c         contain  the positions of the diagonal entries in the
c         output matrix. (i.e., a(iw(k)), ja(iw(k)), k=1,...n,
c         are the values/column indices of the diagonal elements
c         of the output matrix. ).
c
c Notes:
c
c     The column dimension of A is not needed.
c     important: the matrix a may be expanded slightly to allow for
c     additions of nonzero elements to previously nonexisting diagonals.
c     The is no checking as to whether there is enough space appended
c     to the arrays a and ja. if not sure allow for n additional
c     elemnts.
c coded by Y. Saad. Latest version July, 19, 1990
c
      double precision a(*), scal
      integer ja(*), ia(nrow+1),iw(*)
      logical test

      call diapos (nrow,ja,ia,iw)
      icount = 0
      do 1 j=1, nrow
         if (iw(j) .eq. 0) then
            icount = icount+1
         else
            a(iw(j)) = a(iw(j)) + scal
         end if
 1    continue
c
c     if no diagonal elements to insert in data structure return.
c
      if (icount .eq. 0) return
c
c shift the nonzero elements if needed, to allow for created
c diagonal elements.
c
      ko = ia(nrow+1)+icount
c
c     copy rows backward
c
      do 5 ii=nrow, 1, -1
c
c     go through  row ii
c
         k1 = ia(ii)
         k2 = ia(ii+1)-1
         ia(ii+1) = ko
         test = (iw(ii) .eq. 0)
         do 4 k = k2,k1,-1
            j = ja(k)
            if (test .and. (j .lt. ii)) then
               test = .false.
               ko = ko - 1
               a(ko) = scal
               ja(ko) = ii
               iw(ii) = ko
            end if
            ko = ko-1
            a(ko) = a(k)
            ja(ko) = j
 4       continue
c     diagonal element has not been added yet.
         if (test) then
            ko = ko-1
            a(ko) = scal
            ja(ko) = ii
            iw(ii) = ko
         end if
 5    continue
      ia(1) = ko
      return
      end
