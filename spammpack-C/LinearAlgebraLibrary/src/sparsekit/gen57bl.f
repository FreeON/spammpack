      subroutine gen57bl(nx,ny,nz,nfree,na,n,a,ja,ia,iau,stencil)

c*********************************************************************72
c
cc GEN57BL computes the sparse matrix for an elliptic operator.
c
c This subroutine computes the sparse matrix in compressed
c format for the elliptic operator
c
c L u = delx( a . delx u ) + dely ( b . dely u) + delz ( c . delz u ) +
c      delx ( d . u ) + dely (e . u) + delz( f . u ) + g . u
c
c Here u is a vector of nfree componebts and each of the functions
c a, b, c, d, e, f, g   is an (nfree x nfree) matrix depending of
c the coordinate (x,y,z).
c with Dirichlet Boundary conditions, on a rectangular 1-D,
c 2-D or 3-D grid using centered difference schemes.
c
c The functions a, b, ..., g are known through the
c subroutines  afunbl, bfunbl, ..., gfunbl. (user supplied) .
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
c nfree   = number of degrees of freedom per point
c n        = dimension of matrix (output)
c na        = first dimension of array a as declared in calling
c           program. Must be .ge. nfree**2
c
c a, ja, ia = resulting matrix in  row-sparse block-reduced format
c           a(1:nfree**2, j ) contains a nonzero block.
c           ja(j) contains the column number of (1,1) entry of the block.
c
c iau     = integer*n containing the position of the diagonal element
c           in the a, ja, ia structure
c
c stencil =  work array of size (7,nfree**2), used to store
c            local stencils.
c
c
c     stencil (1:7,*) has the following meaning:
c
c     center point = stencil(1)
c     west point   = stencil(2)
c     east point   = stencil(3)
c     south point  = stencil(4)
c     north point  = stencil(5)
c     front point  = stencil(6)
c     back point   = stencil(7)
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
c     implicit double precision (a-h,o-z)
      integer ja(*),ia(*),iau(*)
      double precision a(na,1), stencil(7,1)

      h = 1.0/dble(nx+1)
      kx = 1
      ky = nx
      kz = nx*ny
      nfree2 = nfree*nfree
      iedge = 1
      node = 1
      do 100 iz = 1,nz
         do 90 iy = 1,ny
            do 80 ix = 1,nx
               ia(node) = iedge
               call bsten(nx,ny,nz,ix,iy,iz,nfree,stencil,h)
c     west
               if (ix.gt.1) then
                  ja(iedge)=node-kx
                do 4 k=1,nfree2
                 a(iedge,k) = stencil(2,k)
 4              continue
                  iedge=iedge + 1
               end if
c     south
               if (iy.gt.1) then
                  ja(iedge)=node-ky
                do 5 k=1,nfree2
                 a(iedge,k) = stencil(4,k)
 5              continue
                  iedge=iedge + 1
               end if
c     front plane
               if (iz.gt.1) then
                  ja(iedge)=node-kz
                do 6 k=1,nfree2
                 a(iedge,k) = stencil(6,k)
 6              continue
                  iedge=iedge + 1
               end if
c     center node
               ja(iedge) = node
               iau(node) = iedge
               do 7 k=1,nfree2
                  a(iedge,k) = stencil(1,k)
 7             continue
               iedge = iedge + 1
c     -- upper part
c     east
               if (ix.lt.nx) then
                  ja(iedge)=node+kx
                do 8 k=1,nfree2
                 a(iedge,k) = stencil(3,k)
 8              continue
                  iedge=iedge + 1
               end if
c     north
               if (iy.lt.ny) then
                  ja(iedge)=node+ky
                do 9 k=1,nfree2
                 a(iedge,k) = stencil(5,k)
 9              continue
                  iedge=iedge + 1
               end if
c     back plane
               if (iz.lt.nz) then
                  ja(iedge)=node+kz
                do 10 k=1,nfree2
                     a(iedge,k) = stencil(7,k)
 10              continue
                  iedge=iedge + 1
               end if
c  next node
               node=node+1
 80         continue
 90      continue
 100  continue
c     change numbering of nodes so that each ja(k) will contain the
c     actual column number in the original matrix of entry (1,1) of each
c     block (k).
      do 101 k=1,iedge-1
         ja(k) = (ja(k)-1)*nfree+1
 101  continue

      n = (node-1)*nfree
      ia(node)=iedge
      return
      end
