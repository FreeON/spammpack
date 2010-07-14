      subroutine bndcsr (n,abd,nabd,lowd,ml,mu,a,ja,ia,len,ierr)

c*********************************************************************72
c
cc BNDCSR converts Banded Linpack format to Compressed Sparse Row format.
c
c Banded (Linpack ) format   to    Compressed Sparse Row  format.
c
c on entry:
c
c n      = integer,the actual row dimension of the matrix.
c
c nabd  = first dimension of array abd.
c
c abd   = double precision array containing the values of the matrix stored in
c         banded form. The j-th column of abd contains the elements
c         of the j-th column of  the original matrix,comprised in the
c         band ( i in (j-ml,j+mu) ) with the lowest diagonal located
c         in row lowd (see below).
c
c lowd  = integer. this should be set to the row number in abd where
c         the lowest diagonal (leftmost) of A is located.
c         lowd should be s.t.  ( 1  .le.  lowd  .le. nabd).
c         The routines dgbco, ... of linpack use lowd=2*ml+mu+1.
c
c ml      = integer. equal to the bandwidth of the strict lower part of A
c mu      = integer. equal to the bandwidth of the strict upper part of A
c         thus the total bandwidth of A is ml+mu+1.
c         if ml+mu+1 is found to be larger than nabd then an error
c         message is set. see ierr.
c
c len   = integer. length of arrays a and ja. bndcsr will stop if the
c         length of the arrays a and ja is insufficient to store the
c         matrix. see ierr.
c
c on return:
c
c a,
c ja,
c ia    = input matrix stored in compressed sparse row format.
c
c lowd  = if on entry lowd was zero then lowd is reset to the default
c         value ml+mu+l.
c
c ierr  = integer. used for error message output.
c         ierr .eq. 0 :means normal return
c         ierr .eq. -1 : means invalid value for lowd.
c        ierr .gt. 0 : means that there was not enough storage in a and ja
c         for storing the ourput matrix. The process ran out of space
c         (as indicated by len) while trying to fill row number ierr.
c         This should give an idea of much more storage might be required.
c         Moreover, the first irow-1 rows are correctly filled.
c
c notes:  the values in abd found to be equal to zero
c         (actual test: if (abd(...) .eq. 0.0) are removed.
c         The resulting may not be identical to a csr matrix
c         originally transformed to a bnd format.
c
      double precision a(*),abd(nabd,*), t
      integer ia(n+1),ja(*)

      ierr = 0

      if (lowd .gt. nabd .or. lowd .le. 0) then
         ierr = -1
         return
      end if

      ko = 1
      ia(1) = 1
      do 30 irow=1,n

         i = lowd
          do  20 j=irow-ml,irow+mu
             if (j .le. 0 ) goto 19
             if (j .gt. n) goto 21
             t = abd(i,j)
             if (t .eq. 0.0) goto 19
             if (ko .gt. len) then
               ierr = irow
               return
            end if
            a(ko) = t
            ja(ko) = j
            ko = ko+1
 19         i = i-1
 20      continue
c     end for row irow
 21      ia(irow+1) = ko
 30   continue
      return
      end
