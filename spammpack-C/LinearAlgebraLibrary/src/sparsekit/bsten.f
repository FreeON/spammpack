      subroutine bsten (nx,ny,nz,kx,ky,kz,nfree,stencil,h)

c*********************************************************************72
c
cc BSTEN calculates block stencil values.
c
c  This routine calculates the correct block-stencil values for
c     centered difference discretization of the elliptic operator
c     (block version of stencil)
c
c L u = delx( a delx u ) + dely ( b dely u) + delz ( c delz u ) +
c       d delx ( u ) + e dely (u) + f delz( u ) + g u
c
c   For 2-D problems the discretization formula that is used is:
c
c h**2 * Lu == a(i+1/2,j)*{u(i+1,j) - u(i,j)} +
c             a(i-1/2,j)*{u(i-1,j) - u(i,j)} +
c              b(i,j+1/2)*{u(i,j+1) - u(i,j)} +
c              b(i,j-1/2)*{u(i,j-1) - u(i,j)} +
c              (h/2)*d(i,j)*{u(i+1,j) - u(i-1,j)} +
c              (h/2)*e(i,j)*{u(i,j+1) - u(i,j-1)} +
c              (h/2)*e(i,j)*{u(i,j+1) - u(i,j-1)} +
c              (h**2)*g(i,j)*u(i,j)
c
c      implicit double precision (a-h,o-z)
      double precision stencil(7,*)
      double precision cntr(225), coeff(225),h, hhalf, x, y, z

      if (nfree .gt. 15) then
        WRITE(*,*)'BSTEN  - FATAL ERROR'
        WRITE(*,*)'         Input value of NFREE is greater than 15.'
        STOP
      end if

      nfree2 = nfree*nfree
      do 200 k=1, nfree2
         cntr(k) = 0.0
         do 199 i=1,7
            stencil(i,k) = 0.0
 199     continue
 200  continue

      hhalf = h*0.5
      h2 = h*h
      x = h*dble(kx)
      y = h*dble(ky)
      z = h*dble(kz)
c differentiation wrt x:
      call afunbl(nfree,x+hhalf,y,z,coeff)
      do 1 k=1, nfree2
      stencil(3,k) = stencil(3,k) + coeff(k)
      cntr(k) = cntr(k) + coeff(k)
 1    continue

      call afunbl(nfree,x-hhalf,y,z,coeff)
      do 2 k=1, nfree2
         stencil(2,k) = stencil(2,k) + coeff(k)
         cntr(k) = cntr(k) + coeff(k)
 2    continue

      call dfunbl(nfree,x,y,z,coeff)
      do 3 k=1, nfree2
         stencil(3,k) = stencil(3,k) + coeff(k)*hhalf
         stencil(2,k) = stencil(2,k) - coeff(k)*hhalf
 3    continue
      if (ny .le. 1) goto 99
c
c differentiation wrt y:
c
      call bfunbl(nfree,x,y+hhalf,z,coeff)
      do 4 k=1,nfree2
         stencil(5,k) = stencil(5,k) + coeff(k)
         cntr(k) = cntr(k) + coeff(k)
 4    continue
c
      call bfunbl(nfree,x,y-hhalf,z,coeff)
      do 5 k=1, nfree2
         stencil(4,k) = stencil(4,k) + coeff(k)
         cntr(k) = cntr(k) + coeff(k)
 5    continue
c
      call efunbl(nfree,x,y,z,coeff)
      do 6 k=1, nfree2
         stencil(5,k) = stencil(5,k) + coeff(k)*hhalf
         stencil(4,k) = stencil(4,k) - coeff(k)*hhalf
 6    continue
      if (nz .le. 1) goto 99
c
c differentiation wrt z:
c
      call cfunbl(nfree,x,y,z+hhalf,coeff)
      do 7 k=1, nfree2
         stencil(7,k) = stencil(7,k) + coeff(k)
         cntr(k) = cntr(k) + coeff(k)
 7    continue
c
      call cfunbl(nfree,x,y,z-hhalf,coeff)
      do 8 k=1, nfree2
         stencil(6,k) = stencil(6,k) + coeff(k)
         cntr(k) = cntr(k) + coeff(k)
 8    continue
c
      call ffunbl(nfree,x,y,z,coeff)
      do 9 k=1, nfree2
         stencil(7,k) = stencil(7,k) + coeff(k)*hhalf
         stencil(6,k) = stencil(6,k) - coeff(k)*hhalf
 9    continue
c
c discretization of  product by g:
c
 99   call gfunbl(nfree,x,y,z,coeff)
      do 10 k=1, nfree2
         stencil(1,k) = h2*coeff(k) - cntr(k)
 10   continue
c
      return
      end
