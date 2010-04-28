      subroutine diacsr (n,job,idiag,diag,ndiag,ioff,a,ja,ia)

c*********************************************************************72
c
cc DIACSR converts diagonal format to compressed sparse row.
c
c DIACSR extracts the idiag most important diagonals from the
c input matrix a, ja, ia, i.e, those diagonals of the matrix which have
c the largest number of nonzero elements. If requested (see job),
c the rest of the matrix is put in a the output matrix ao, jao, iao
c
c on entry:
c
c n      = integer. dimension of the matrix a.
c job      = integer. job indicator with the following meaning.
c         if (job .eq. 0) then check for each entry in diag
c         whether this entry is zero. If it is then do not include
c         in the output matrix. Note that the test is a test for
c         an exact arithmetic zero. Be sure that the zeros are
c         actual zeros otherwise this would not
c         work.
c
c idiag = integer equal to the number of diagonals to be extracted.
c         Note: on return idiag may be modified.
c
c diag  = double precision array of size (ndiag x idiag) containing the diagonals
c         of A on return.
c
c ndiag = integer equal to the first dimension of array diag.
c
c ioff  = integer array of length idiag, containing the offsets of the
c           diagonals to be extracted.
c
c on return:
c
c a,
c ja,
c ia    = matrix stored in a, ja, ia, format
c
c Note:
c the arrays a and ja should be of length n*idiag.
c
      double precision diag(ndiag,idiag), a(*), t
      integer ia(*), ja(*), ioff(*)

      ia(1) = 1
      ko = 1
      do 80 i=1, n
         do 70 jj = 1, idiag
            j = i+ioff(jj)
            if (j .lt. 1 .or. j .gt. n) goto 70
            t = diag(i,jj)
            if (job .eq. 0 .and. t .eq. 0.0) goto 70
            a(ko) = t
            ja(ko) = j
            ko = ko+1
 70      continue
         ia(i+1) = ko
 80   continue
      return
      end
