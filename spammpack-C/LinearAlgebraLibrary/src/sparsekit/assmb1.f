      subroutine assmb1 (u,nu,a,ja,ia,fu,f,nx,nelx,ijk,nodcode,
     *     node,iwk,jwk)

c*********************************************************************72
c
cc ASSMB1 assembles a finite element matrix in the CSR format.
c
c u       = unassembled matrix u(na,node,node)
c nu       = 1-st dimension of u
c a,ja,ia= assembled matrix on output
c fu       = unassembled right hand side
c f      = right hand side (global load vector) assembled
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
c x,y   = double precision arrays containing the $x$ and $y$ coordinates
c        resp. of the nodes.
c         K11, K22, and K12 at that element.
c iwk,jwk = two integer work arrays.
c ierr      = error message integer .
c        ierr = 0 --> normal return
c        ierr = 1 --> negative area encountered (due to bad
c                 numbering of nodes of an element- see
c               message printed in unit iout). not used..
c iout      = output unit (not used here).
c
      implicit double precision (a-h,o-z)
      double precision u(nu,node,node),a(*),fu(node,*),f(*)
      integer ja(*),ia(*),ijk(node,*),iwk(*),jwk(*),nodcode(*)
c     max number of nonzeros per row allowed  = 200
c
c     initialize
c
      do 100 i=1,nx
         f(i) = 0.0
 100  continue
c
c     initialize  pointer arrays.
c
      do 5 k=1,nx+1
         ia(k) = 1
         jwk(k) = 0
 5    continue
      do 6 k=1,nelx
         do 59 j=1,node
            knod = ijk(j,k)
            ia(knod) = ia(knod) + 1
 59      continue
 6    continue

      do 7 k=1, nx
         if (nodcode(k) .ge.1 ) ia(k)=ia(k)+1
 7    continue

      ksav = ia(1)
      ia(1) = 1
      do 101 j=2, nx+1
         ksavn = ia(j)
         ia(j) = ia(j-1) +  ksav
         iwk(j-1) = ia(j-1)-1
         ksav = ksavn
 101  continue
c
c     main loop
c
      do 102 nel=1, nelx
c
c     get nodal points
c
         do 120 ka=1, node
            ii = ijk(ka,nel)
            f(ii) = f(ii) + fu(ka,nel)
c
c     unpack row into jwk1
c
            irowst = ia(ii)
            ilast  = iwk(ii)
            do 109 k=irowst,ilast
               jwk(ja(k)) = k
 109        continue

            do 108 kb = 1,node
c
c     column number = jj
c
               jj = ijk(kb,nel)
               k = jwk(jj)
               if (k .eq. 0) then
                  ilast = ilast+1
                  jwk(jj) = ilast
                  ja(ilast) = jj
                  a(ilast) = u(nel,ka,kb)
               else
                  a(k) = a(k) + u(nel,ka,kb)
               end if
 108        continue
c     refresh jwk
            do 119 k=irowst,ilast
               jwk(ja(k)) = 0
 119        continue
            iwk(ii) = ilast
 120     continue

 102  continue
      return
      end
