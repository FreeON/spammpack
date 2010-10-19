      subroutine csrncf ( nrow, a, ja, ia, maxnz, nonz, coef, jcoef,
     &  ierr )

c*********************************************************************72
c
cc CSRNCF converts CSR to NSPCG NCF format.
c
c  Discussion:
c
c    This routine converts a matrix stored in the general A, JA, IA
c    compressed sparse row format into the Nonsymmetric Coordinate Format
c    used as storage format 5 by NSPCG.
c
c  Modified:
c
c    25 October 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer NROW, the row dimension of the matrix A.
c
c    Input, real A(*), integer JA(*), IA(NROW+1), the matrix, stored in
c    compressed sparse row format.
c
c    Input, integer MAXNZ, the maximum number of nonzeros allowed for
c    in the storage of COEF and JCOEF.
c
c    Output, integer NONZ, the actual number of nonzeros encountered.
c
c    Output, real COEF(MAXNZ), the values of the matrix A in NCF format.
c
c    Output, integer JCOEF(MAXNZ,2), the row and column indices of each
c    entry in COEF.
c
c    Output, integer IERR, an error flag.
c    0 = correct return.
c    nonzero means that MAXNZ < NONZ.
c
      implicit none

      integer maxnz
      integer nrow

      double precision a(*)
      double precision coef(maxnz)
      integer i
      integer ia(nrow+1)
      integer ierr
      integer j
      integer ja(*)
      integer jcoef(maxnz,2)
      integer k
      integer k1
      integer k2
      integer nonz

      ierr = 0
c
c  Initialize COEF and JCOEF.
c
      do i = 1, maxnz
        coef(i) = 0.0D+00
      end do

      do j = 1, 2
        do i = 1, maxnz
          jcoef(i,j) = 0
        end do
      end do
c
c  The first N entries are reserved for the diagonals.
c
      do i = 1, nrow
        jcoef(i,1:2) = i
      end do

      nonz = nrow

      if ( maxnz .lt. nonz ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CSRNCF - Fatal error!'
        write ( *, '(a)' ) '  MAXNZ < NONZ.'
        ierr = 1
        return
      end if

      do i = 1, nrow

        k1 = ia(i)
        k2 = ia(i+1) - 1

        do k = k1, k2

          if ( ja(k) .eq. i ) then

            coef(i) = coef(i) + a(k)

          else if ( 0.0D+00 .ne. a(k) ) then

            nonz = nonz + 1

            if ( maxnz .lt. nonz ) then
              write ( *, '(a)' ) ' '
              write ( *, '(a)' ) 'CSRNCF - Fatal error!'
              write ( *, '(a)' ) '  MAXNZ < NONZ.'
              ierr = 1
              return
            end if

            coef(nonz) = a(k)
            jcoef(nonz,1) = i
            jcoef(nonz,2) = ja(k)

          end if

        end do

      end do

      return
      end
