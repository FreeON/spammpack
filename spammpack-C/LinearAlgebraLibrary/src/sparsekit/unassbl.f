      subroutine unassbl (a,na,f,nx,nelx,ijk,nodcode,
     *                     node,x,y,ierr,xyk)

c*********************************************************************72
c
cc UNASSBL ?
c
c a      = un-assembled matrix on output
c na       = 1-st dimension of a.  a(na,node,node)
c
c f      = right hand side (global load vector) in un-assembled form
c nx     = number of nodes at input
c nelx       = number of elements at input
c ijk       = connectivity matrix: for node k, ijk(*,k) point to the
c          nodes of element k.
c node       = total number of nodal points in each element
c         also second dimension of a.
c
c nodcode= boundary information list for each node with the
c         following meaning:
c      nodcode(i) = 0 -->  node i is internal
c      nodcode(i) = 1 -->  node i is a boundary but not a corner point
c      nodcode(i) = 2 -->  node i is a corner point (corner points
c
c x,y   = double precision arrays containing the $x$ and $y$ coordinates
c        resp. of the nodes.
c         K11, K22, and K12 at that element.
c ierr      = error message integer .
c        ierr = 0 --> normal return
c        ierr = 1 --> negative area encountered (due to bad
c                 numbering of nodes of an element-
c               message printed in unit iout).
c iout      = output unit (not used here).
c
c xyk      = subroutine defining the material properties at each
c         element. Form:
c       call xyk(nel,xyke,x,y,ijk,node)
c
      implicit double precision (a-h,o-z)
        dimension a(na,node,node),ijk(node,1),x(1),y(1),f(node,1),
     *          ske(3,3),fe(3),xe(3),ye(3),xyke(2,2)
            integer nodcode(1)
       external xyk
c max number of nonzeros allowed  = 200
c
c   initialize
c
      do i=1, node
        do j=1, nx
          f(i,j) = 0.0
        end do
      end do
c
c main loop
c
      do 102 nel=1, nelx
c
c get coordinetes of nodal points
c
      do 104 i=1, node
      j = ijk(i,nel)
      xe(i) = x(j)
      ye(i) = y(j)
 104      continue
c
c compute determinant
c
       det=xe(2)*(ye(3)-ye(1))+xe(3)*(ye(1)-ye(2))+xe(1)*(ye(2)-ye(3))
c        print *, 'nel', nel, ' det = ' , del
c
c set material properties
c
      call xyk(nel,xyke,x,y,ijk,node)
c
c construct element stiffness matrix
c
      ierr = 0
      call estif3(nel,ske,fe,det,xe,ye,xyke,ierr)
      if (ierr .ne. 0) return
c      write (8,'(9f8.4)') ((ske(i,j),j=1,3),i=1,3)
c assemble: add element stiffness matrix to global matrix
c
      do 120 ka=1, node
            f(ka,nel) = fe(ka)
        do 108 kb = 1,node
            a(nel,ka,kb) = ske(ka,kb)
 108      continue
 120      continue
 102      continue
        return
      end
