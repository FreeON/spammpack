      subroutine csrdia (n,idiag,job,a,ja,ia,ndiag,
     *                   diag,ioff,ao,jao,iao,ind)

c*********************************************************************72
c
cc CSRDIA converts Compressed Sparse Row to diagonal format.
c
c CSRDIA extracts  idiag diagonals  from the  input matrix a,
c a, ia, and puts the rest of  the matrix  in the  output matrix ao,
c jao, iao.  The diagonals to be extracted depend  on the  value of job
c (see below for details.)  In  the first  case, the  diagonals to be
c extracted are simply identified by  their offsets  provided in ioff
c by the caller.  In the second case, the  code internally determines
c the idiag most significant diagonals, i.e., those  diagonals of the
c matrix which  have  the  largest  number  of  nonzero elements, and
c extracts them.
c
c on entry:
c
c n      = dimension of the matrix a.
c idiag = integer equal to the number of diagonals to be extracted.
c         Note: on return idiag may be modified.
c a, ja,
c    ia = matrix stored in a, ja, ia, format
c job      = integer. serves as a job indicator.  Job is better thought
c         of as a two-digit number job=xy. If the first (x) digit
c         is one on entry then the diagonals to be extracted are
c         internally determined. In this case csrdia exctracts the
c         idiag most important diagonals, i.e. those having the largest
c         number on nonzero elements. If the first digit is zero
c         then csrdia assumes that ioff(*) contains the offsets
c         of the diagonals to be extracted. there is no verification
c         that ioff(*) contains valid entries.
c         The second (y) digit of job determines whether or not
c         the remainder of the matrix is to be written on ao,jao,iao.
c         If it is zero  then ao, jao, iao is not filled, i.e.,
c         the diagonals are found  and put in array diag and the rest is
c         is discarded. if it is one, ao, jao, iao contains matrix
c         of the remaining elements.
c         Thus:
c         job= 0 means do not select diagonals internally (pick those
c                defined by ioff) and do not fill ao,jao,iao
c         job= 1 means do not select diagonals internally
c                      and fill ao,jao,iao
c         job=10 means  select diagonals internally
c                      and do not fill ao,jao,iao
c         job=11 means select diagonals internally
c                      and fill ao,jao,iao
c
c ndiag = integer equal to the first dimension of array diag.
c
c on return:
c
c
c idiag = number of diagonals found. This may be smaller than its value
c         on entry.
c diag  = double precision array of size (ndiag x idiag) containing the diagonals
c         of A on return
c
c ioff  = integer array of length idiag, containing the offsets of the
c           diagonals to be extracted.
c ao, jao
c  iao  = remainder of the matrix in a, ja, ia format.
c work arrays:
c
c ind   = integer array of length 2*n-1 used as integer work space.
c         needed only when job.ge.10 i.e., in case the diagonals are to
c         be selected internally.
c
c Notes:
c
c    1) The algorithm is in place: ao, jao, iao can be overwritten on
c       a, ja, ia if desired
c    2) When the code is required to select the diagonals (job .ge. 10)
c       the selection of the diagonals is done from left to right
c       as a result if several diagonals have the same weight (number
c       of nonzero elemnts) the leftmost one is selected first.
c
      double precision diag(ndiag,idiag), a(*), ao(*)
      integer ia(*), ind(*), ja(*), jao(*), iao(*), ioff(*)

      job1 = job/10
      job2 = job-job1*10
      if (job1 .eq. 0) goto 50
      n2 = n+n-1
      call infdia(n,ja,ia,ind,idum)
c  determine diagonals to  accept.
c
      ii = 0
 4    ii=ii+1
      jmax = 0
      do 41 k=1, n2
         j = ind(k)
         if (j .le. jmax) goto 41
         i = k
         jmax = j
 41   continue
      if (jmax .le. 0) then
         ii = ii-1
         goto 42
      end if
      ioff(ii) = i-n
      ind(i) = - jmax
      if (ii .lt.  idiag) goto 4
 42   idiag = ii
c  initialize diago to zero
 50   continue
      do 55 j=1,idiag
         do 54 i=1,n
            diag(i,j) = 0.0
 54      continue
 55   continue

      ko = 1
c
c extract diagonals and accumulate remaining matrix.
c
      do 6 i=1, n
         do 51 k=ia(i),ia(i+1)-1
            j = ja(k)
            do 52 l=1,idiag
               if (j-i .ne. ioff(l)) goto 52
               diag(i,l) = a(k)
               goto 51
 52         continue
c  append element not in any diagonal to ao,jao,iao
            if (job2 .eq. 0) goto 51
            ao(ko) = a(k)
            jao(ko) = j
            ko = ko+1
 51      continue
         if (job2 .ne. 0 ) ind(i+1) = ko
 6    continue
      if (job2 .eq. 0) return
c     finish with iao
      iao(1) = 1
      do 7 i=2,n+1
         iao(i) = ind(i)
 7    continue
      return
      end
