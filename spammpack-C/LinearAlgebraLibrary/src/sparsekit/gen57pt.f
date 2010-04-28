      subroutine gen57pt(nx,ny,nz,a,ja,ia,iau,stencil)

c*********************************************************************72
c
cc GEN57PT computes the compressed sparse matrix for an elliptic operator.
c
c L u = delx( a delx u ) + dely ( b dely u) + delz ( c delz u ) +
c       d delx ( u ) + e dely (u) + f delz( u ) + g u
c
c with Dirichlet Boundary conditions, on a rectangular 1-D,
c 2-D or 3-D grid using centered difference schemes.
c
c The functions a, b, ..., g are known through the
c subroutines  afun, bfun, ..., gfun.
c note that to obtain the correct matrix, any function that is not
c needed should be set to zero. For example for two-dimensional
c problems, nz should be set to 1 and the functions cfun and ffun
c should be zero functions.
c
c uses natural ordering, first x direction, then y, then z
c mesh size h is uniform and determined by grid points
c in the x-direction.
c
c parameters:
c
c nx      = number of points in x direction
c ny        = number of points in y direction
c nz        = number of points in z direction
c
c a, ja, ia =  resulting matrix in row-sparse format
c
c iau     = integer*n containing the poisition of the diagonal element
c           in the a, ja, ia structure
c
c stencil =  work array of size 7, used to store local stencils.
c
c
c     stencil [1:7] has the following meaning:
c
c     center point = stencil(1)
c     west point = stencil(2)
c     east point = stencil(3)
c     south point = stencil(4)
c     north point = stencil(5)
c     front point = stencil(6)
c     back point = stencil(7)
c
c
c                           st(5)
c                            |
c                            |
c                            |
c                            |          .st(7)
c                            |     .
c                            | .
c         st(2) ----------- st(1) ---------- st(3)
c                       .    |
c                   .        |
c               .            |
c            st(6)           |
c                            |
c                            |
c                           st(4)
c
c
      integer ja(*),ia(*),iau(*)
      double precision a(*), stencil(*), h

      h = 1.0/dble(nx+1)
      kx = 1
      ky = nx
      kz = nx*ny
      iedge = 1
      node = 1
      do 100 iz = 1,nz
         do 90 iy = 1,ny
            do 80 ix = 1,nx
               ia(node) = iedge
               call getsten(nx,ny,nz,ix,iy,iz,stencil,h)
c     west
               if (ix.gt.1) then
                  ja(iedge)=node-kx
              a(iedge) = stencil(2)
                  iedge=iedge + 1
               end if
c     south
               if (iy.gt.1) then
                  ja(iedge)=node-ky
              a(iedge) = stencil(4)
                  iedge=iedge + 1
               end if
c     front plane
               if (iz.gt.1) then
                  ja(iedge)=node-kz
              a(iedge) = stencil(6)
                  iedge=iedge + 1
               end if
c     center node
               ja(iedge) = node
               iau(node) = iedge
               a(iedge) = stencil(1)
               iedge = iedge + 1
c     -- upper part
c     east
               if (ix.lt.nx) then
                  ja(iedge)=node+kx
              a(iedge) = stencil(3)
                  iedge=iedge + 1
               end if
c     north
               if (iy.lt.ny) then
                  ja(iedge)=node+ky
              a(iedge) = stencil(5)
                  iedge=iedge + 1
               end if
c     back plane
               if (iz.lt.nz) then
                  ja(iedge)=node+kz
                  a(iedge) = stencil(7)
                  iedge=iedge + 1
               end if
c  next node
               node=node+1
 80         continue
 90      continue
 100  continue
      ia(node)=iedge
      return
      end
