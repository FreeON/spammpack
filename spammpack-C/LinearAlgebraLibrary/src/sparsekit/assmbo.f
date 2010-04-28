      subroutine assmbo (nx,nelx,node,ijk,nodcode,x,y,
     *             a,ja,ia,f,iwk,jwk,ierr,xyk)

c*********************************************************************72
c
cc ASSMBO assembles a finite element matrix.
c
c
c a,ja,ia= assembled matrix on output
c f      = right hand side (global load vector)
c nx     = number of nodes at input
c nelx       = number of elements at input
c ijk       = connectivity matrix: for node k, ijk(*,k) point to the
c          nodes of element k.
c node       = total number of nodal points in each element
c
c nodcode= boundary information list for each node with the
c         following meaning:
c      nodcode(i) = 0 -->  node i is internal
c      nodcode(i) = 1 -->  node i is a boundary but not a corner point
c      nodcode(i) = 2 -->  node i is a corner point (corner points
c
c mp      = group material number for each element.
c
c x,y   = double precision arrays containing the $x$ and $y$ coordinates
c        resp. of the nodes.
c iwk,jwk = two integer work arrays.
c ierr      = error message integer .
c        ierr = 0 --> normal return
c        ierr = 1 --> negative area encountered (due to bad
c                 numbering of nodes of an element- see
c               message printed in unit iout ).
c iout      = output unit (not used here).
c
c xyk      = routine defining the material properties at each
c         element. Form:
c       call xyk(nel,xyke,x,y,ijk,node) with on return
c         xyke =  material constant matrices.
c         for each element nel, xyke(1,nel),xyke(2,nel)
c         and xyke(3,nel) represent the constants
c         K11, K22, and K12 at that element.
c
      implicit double precision  (a-h,o-z)
      dimension a(1),ijk(node,1),x(1),y(1),f(1),ske(3,3),fe(3),
     *      xe(3),ye(3),xyke(2,2),iwk(1),jwk(1)
      integer ia(1), ja(1), nodcode(1)
      external xyk
c max number of nonzeros allowed  = 200
c
c   initialize
c
      do 100 i=1,nx
 100          f(i) = 0.0
c initialize  pointer arrays.
      do 5 k=1,nx+1
      ia(k) = 1
      jwk(k) = 0
 5      continue
      do 6 k=1,nelx
      do 59 j=1,node
      knod = ijk(j,k)
 59      ia(knod) = ia(knod) + 1
 6      continue

      do 7 k=1, nx
      if (nodcode(k) .ge.1 ) ia(k)=ia(k)+1
 7      continue

      ksav = ia(1)
      ia(1) = 1
      do 101 j=2, nx+1
            ksavn = ia(j)
            ia(j) = ia(j-1) +  ksav
            iwk(j-1) = ia(j-1)-1
 101            ksav = ksavn
c
c main loop
c
      do 102 nel=1, nelx
c
c get coordinates of nodal points
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
c
c assemble: add element stiffness matrix to global matrix
c
      do 120 ka=1, node
            ii = ijk(ka,nel)
          f(ii) = f(ii) + fe(ka)
c
c unpack row into jwk1
c
         irowst = ia(ii)
         ilast  = iwk(ii)
           do 109 k=irowst,ilast
         jwk(ja(k)) = k
 109         continue
c
          do 108 kb = 1,node
c
c column number = jj
c
            jj = ijk(kb,nel)
           k = jwk(jj)
          if (k .eq. 0) then
             ilast = ilast+1
             jwk(jj) = ilast
             ja(ilast) = jj
             a(ilast) = ske(ka,kb)
          else
              a(k) = a(k) + ske(ka,kb)
          end if
 108         continue
c refresh jwk
           do 119 k=irowst,ilast
         jwk(ja(k)) = 0
 119         continue
           iwk(ii) = ilast
 120    continue
c
 102      continue
        return
      end
