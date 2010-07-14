      subroutine bsort2 (w, ind, n, ncut)

c*********************************************************************72
c
cc BSORT2 returns the NCUT largest elements of an array, using bubble sort.
c
c simple bubble sort for getting the ncut largest
c elements in modulus, in array w. ind is sorted accordingly.
c (Ought to be replaced by a more efficient sort especially
c if ncut is not that small).
c
      integer n, ncut, ind(*)
      double precision w(*)
      logical test
      integer i, j, iswp
      double precision wswp

      i = 1
 1    test = .false.
      do 2 j = n-1,i,-1
         if (abs(w(j+1))  .gt. abs(w(j)) ) then
c  swap.
            wswp = w(j)
            w(j) = w(j+1)
            w(j+1) = wswp
c
c  reorder original ind array accordingly.
c
            iswp = ind(j)
            ind(j) = ind(j+1)
            ind(j+1) = iswp
c  set indicator that sequence is still unsorted--
            test = .true.
         end if
 2    continue
      i = i+ 1
      if (test .and. i .le. ncut) goto 1
      return
      end
