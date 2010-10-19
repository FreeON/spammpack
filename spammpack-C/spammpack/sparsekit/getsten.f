      subroutine getsten(nx,ny,nz,kx,ky,kz,stencil,h)

c*********************************************************************72
c
cc GETSTEN calculates the stencil for centered elliptic discretization.
c
c     This subroutine calculates the correct stencil values for
c     centered difference discretization of the elliptic operator
c
c L u = delx( a delx u ) + dely ( b dely u) + delz ( c delz u ) +
c      delx ( d u ) + dely (e u) + delz( f u ) + g u
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
c               (h**2)*g(i,j)*u(i,j)
c
      double precision stencil(*), h, hhalf,cntr, afun, bfun, cfun, dfun,
     *      efun, ffun, gfun, x, y, z, coeff

      do k=1,7
         stencil(k) = 0.0
      end do

      hhalf = h*0.5
      x = h*dble(kx)
      y = h*dble(ky)
      z = h*dble(kz)
      cntr = 0.0
c     differentiation wrt x:
      coeff = afun(x+hhalf,y,z)
      stencil(3) = stencil(3) + coeff
      cntr = cntr + coeff

      coeff = afun(x-hhalf,y,z)
      stencil(2) = stencil(2) + coeff
      cntr = cntr + coeff

      coeff = dfun(x,y,z)*hhalf
      stencil(3) = stencil(3) + coeff
      stencil(2) = stencil(2) - coeff
      if (ny .le. 1) goto 99
c
c     differentiation wrt y:
c
      coeff = bfun(x,y+hhalf,z)
      stencil(5) = stencil(5) + coeff
      cntr = cntr + coeff

      coeff = bfun(x,y-hhalf,z)
      stencil(4) = stencil(4) + coeff
      cntr = cntr + coeff

      coeff = efun(x,y,z)*hhalf
      stencil(5) = stencil(5) + coeff
      stencil(4) = stencil(4) - coeff
      if (nz .le. 1) goto 99
c
c differentiation wrt z:
c
      coeff = cfun(x,y,z+hhalf)
      stencil(7) = stencil(7) + coeff
      cntr = cntr + coeff

      coeff = cfun(x,y,z-hhalf)
      stencil(6) = stencil(6) + coeff
      cntr = cntr + coeff

      coeff = ffun(x,y,z)*hhalf
      stencil(7) = stencil(7) + coeff
      stencil(6) = stencil(6) - coeff
c
c discretization of  product by g:
c
 99   coeff = gfun(x,y,z)
      stencil(1) = h*h*coeff - cntr
      return
      end
