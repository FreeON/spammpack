      subroutine dcsort(ival, n, icnt, index, ilo, ihi)

c*********************************************************************72
c
cc DCSORT computes a sorting permutation for a vector.
c
c     Specifications for arguments:
c
c    This routine computes a permutation which, when applied to the
c    input vector ival, sorts the integers in ival in descending
c    order.  The permutation is represented by the vector index.  The
c    permuted ival can be interpreted as follows:
c      ival(index(i-1)) .ge. ival(index(i)) .ge. ival(index(i+1))
c
c    A specialized sort, the distribution counting sort, is used
c    which takes advantage of the knowledge that
c        1)  The values are in the (small) range [ ilo, ihi ]
c        2)  Values are likely to be repeated often
c
c    contributed to SPARSKIT by Mike Heroux. (Cray Research)
c
c Usage:
c
c     call dcsort( ival, n, icnt, index, ilo, ihi )
c
c Arguments:
c
c    ival  integer array (input)
c          On entry, ia is an n dimensional array that contains
c          the values to be sorted.  ival is unchanged on exit.
c
c    n     integer (input)
c          On entry, n is the number of elements in ival and index.
c
c    icnt  integer (work)
c          On entry, is an integer work vector of length
c          (ihi - ilo + 1).
c
c    index integer array (output)
c          On exit, index is an n-length integer vector containing
c          the permutation which sorts the vector ival.
c
c    ilo   integer (input)
c          On entry, ilo is .le. to the minimum value in ival.
c
c    ihi   integer (input)
c          On entry, ihi is .ge. to the maximum value in ival.
c
c Remarks:
c
c         The permutation is NOT applied to the vector ival.
c
c Author:
c
c    Michael Heroux
c    Sandra Carney
c       Mathematical Software Research Group
c       Cray Research, Inc.
c
c References:
c    Knuth, Donald E., "The Art of Computer Programming, Volume 3:
c    Sorting and Searching," Addison-Wesley, Reading, Massachusetts,
c    1973, pp. 78-79.
c
c Revision history:
c    05/09/90: Original implementation.  A variation of the
c              Distribution Counting Sort recommended by
c              Sandra Carney. (Mike Heroux)
c
c
      integer n, ilo, ihi, ival(n), icnt(ilo:ihi), index(n)
      integer i, j, ivalj

      do i = ilo, ihi
        icnt(i) = 0
      end do

      do i = 1, n
        icnt(ival(i)) = icnt(ival(i)) + 1
      end do

      do i = ihi-1,ilo,-1
        icnt(i) = icnt(i) + icnt(i+1)
      end do

      do j = n, 1, -1
        ivalj = ival(j)
        index(icnt(ivalj)) = j
        icnt(ivalj) = icnt(ivalj) - 1
      end do

      return
      end
