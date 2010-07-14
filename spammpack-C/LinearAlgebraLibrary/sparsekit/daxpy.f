      subroutine daxpy ( n, da, dx, incx, dy, incy )

c*********************************************************************72
c
cc DAXPY computes constant times a vector plus a vector.
c
c  Discussion:
c
c    This routine uses double precision real arithmetic.
c
c    This routine uses unrolled loops for increments equal to one.
c
c  Modified:
c
c    08 July 2007
c
c  Author:
c
c    Jack Dongarra
c
c  Reference:
c
c    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
c    LINPACK User's Guide,
c    SIAM, 1979,
c    ISBN13: 978-0-898711-72-1,
c    LC: QA214.L56.
c
c    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
c    Basic Linear Algebra Subprograms for FORTRAN usage,
c    ACM Transactions on Mathematical Software,
c    Volume 5, Number 3, pages 308-323, 1979.
c
c  Parameters:
c
c    Input, integer N, the number of elements in DX and DY.
c
c    Input, double precision DA, the multiplier of DX.
c
c    Input, double precision DX(*), the first vector.
c
c    Input, integer INCX, the increment between successive entries of DX.
c
c    Input/output, double precision DY(*), the second vector.
c    On output, DY(*) has been replaced by DY(*) + DA * DX(*).
c
c    Input, integer INCY, the increment between successive entries of DY.
c
      implicit none

      double precision da
      double precision dx(*)
      double precision dy(*)
      integer i
      integer incx
      integer incy
      integer ix
      integer iy
      integer m
      integer n

      if ( n .le. 0 ) then
        return
      end if

      if ( da .eq. 0.0d0 ) then
        return
      end if

      if(incx.eq.1.and.incy.eq.1)go to 20
c
c  code for unequal increments or equal increments
c  not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1

      do i = 1,n
        dy(iy) = dy(iy) + da*dx(ix)
        ix = ix + incx
        iy = iy + incy
      end do

      return
c
c  code for both increments equal to 1
c
c
c  clean-up loop
c
   20 m = mod(n,4)
      if( m .eq. 0 ) go to 40
      do i = 1,m
        dy(i) = dy(i) + da*dx(i)
      end do

      if( n .lt. 4 ) return

   40 continue

      do i = m + 1, n, 4
        dy(i) = dy(i) + da*dx(i)
        dy(i + 1) = dy(i + 1) + da*dx(i + 1)
        dy(i + 2) = dy(i + 2) + da*dx(i + 2)
        dy(i + 3) = dy(i + 3) + da*dx(i + 3)
      end do

      return
      end
