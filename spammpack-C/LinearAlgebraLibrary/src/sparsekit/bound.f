      subroutine bound (nx,nelx,ijk,nodcode,node,nint,iperm,
     *             x,y,wk,iwk)

c*********************************************************************72
c
cc BOUND counts the number of boundary points.
c
c  this routine counts the number of boundary points and
c reorders the points in such a way that the boundary nodes
c are last.
c
c nx, nelx, ijk, nodcode, node: see other routines
c nint = on return the number of points on the boundary
c iperm = permutation array from old orderin to new ordering,
c iwk   = reverse permutation array or return.
c wk      = double precision work array
c On return
c x, y, nodecode, are permuted
c ijk  is updated according to new oerdering.
c nint = number of interior points.
c
      implicit double precision  (a-h,o-z)
      dimension ijk(node,1),x(1),y(1),wk(1),iwk(1),iperm(1),
     *              nodcode(1)
c max number of nonzeros allowed  = 200
c put all boundary points at the end, backwards
      nint = 1
      nbound = nx
      do 1 j=1, nx
      if (nodcode(j) .eq. 0) then
        iperm(nint) = j
        nint = nint+1
      else
      iperm(nbound) = j
      nbound = nbound-1
       end if
 1      continue

      nint = nint-1
c
c permute x's
c
      do 2 k=1, nx
      wk(k) = x(k)
 2      continue
      do 3 k=1,nx
        x(k) = wk(iperm(k))
 3      continue
c
c permute the y's
c
      do 4 k=1, nx
      wk(k) = y(k)
 4      continue
      do 5 k=1, nx
      y(k) = wk(iperm(k))
 5      continue
c
c permute the boundary information
c
      do 6 k=1, nx
      iwk(k) = nodcode(k)
 6      continue
      do 7 k=1,nx
      nodcode(k) = iwk(iperm(k))
 7      continue
c
c get reverse permutation
c
      do 8 k=1, nx
        iwk(iperm(k)) = k
 8      continue
c
c update the elements connectivity matrix
c
      do 10 nel = 1, nelx
         do 9 j=1, node
            knod = ijk(j,nel)
              ijk(j,nel) = iwk(knod)
 9         continue
 10      continue
      return
      end
