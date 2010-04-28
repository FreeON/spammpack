      subroutine gradi3(nel, xe, ye, dn, det,ierr)

c*********************************************************************72
c
cc GRADI3 constructs the first derivative of the shape functions.
c
c arguments:
c nel      = element nuumber
c xy, ye= coordinates of the three nodal points in an element.
c dn      = gradients (1-st derivatives) of the shape functions.
c area      = area of the triangle
c
c
      PARAMETER (TOL=1.0E-17)
c
      dimension xe(3), ye(3), dn(3,2)
c
c compute area
c
      ierr = 0
      if (det .le. TOL) goto 100
c
      dn(1,1) = (ye(2)-ye(3))/det
      dn(2,1) = (ye(3)-ye(1))/det
      dn(3,1) = (ye(1)-ye(2))/det
      dn(1,2) = (xe(3)-xe(2))/det
      dn(2,2) = (xe(1)-xe(3))/det
      dn(3,2) = (xe(2)-xe(1))/det
c
      return
c
 100      continue
      ierr = 3
c       write(iout,*) ' ** Error-negative area encountered at elmt: '
c      write(iout,*) nel,(xe(i),ye(i),i=1,3)
      return
      end
