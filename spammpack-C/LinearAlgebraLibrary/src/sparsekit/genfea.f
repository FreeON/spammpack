      subroutine genfea (nx,nelx,node,job,x,y,ijk,nodcode,fs,nint,
     *     a,ja,ia,f,iwk,jwk,ierr,xyk)

c*********************************************************************72
c
cc GENFEA generates finite element matrices for heat conduction problems.
c
c                  - Div ( K(x,y) Grad u ) = f
c                    u = 0 on boundary
c
c (with Dirichlet boundary conditions). The matrix is returned
c assembled in compressed sparse row format. See genfeu for
c matrices in unassembled form. The user must provide the grid,
c (coordinates x, y and connectivity matrix ijk) as well as some
c information on the nodes (nodcode) and the material properties
c (the function K(x,y) above) in the form of a subroutine xyk.
c
c
c
c on entry:
c
c nx          = integer . the number of nodes in the grid .
c nelx          = integer . the number of elements in the grid.
c node      = integer = the number of nodes per element (should be
c             set to three in this version). also the first dimension
c             of ijk
c job          = integer. If job=0, it is assumed that there is no heat
c             source (i.e. fs = 0) and the right hand side
c             produced will therefore be a zero vector.
c             If job = 1 on entry then the contributions from the
c             heat source in each element are taken into account.
c
c x, y      = two double precision arrays containing the coordinates of the nodes.
c
c ijk       =  an integer array containing the connectivity matrix.
c              ijk(i,nel), i=1,2,..node, is the list of the nodes
c              constituting the element nel, ans listed in
c              counter clockwise order.
c
c nodcode   = an integer array containing the boundary information for
c             each node with the following meaning.
c      nodcode(i) = 0 -->  node i is internal
c      nodcode(i) = 1 -->  node i is a boundary but not a corner point
c      nodcode(i) = 2 -->  node i is a corner node. [This node and the
c             corresponmding element are discarded.]
c
c fs          = double precision array of length nelx on entry containing the heat
c             source for each element (job = 1 only)
c
c xyk          = subroutine defining the material properties at each
c            element. Form:
c             call xyk(nel,xyke,x,y,ijk,node) with on return
c             xyke =  material constant matrices.
c            for each element nel, xyke(1,nel),xyke(2,nel)
c             and xyke(3,nel) represent the constants
c             K11, K22, and K12 at that element.
c
c on return
c
c nint          = integer. The number of active (nonboundary) nodes. Also
c             equal to the dimension of the assembled matrix.
c
c a, ja, ia = assembled matrix in compressed sparse row format.
c
c f          = double precision array containing the right hand for the linears
c             system to solve.
c
c ierr          = integer. Error message. If (ierr .ne. 0) on return
c             it means that one of the elements has a negative or zero
c             area probably because of a bad ordering of the nodes
c             (see ijk above). Use the subroutine chkelmt to reorder
c             the nodes properly if necessary.
c iwk, jwk  = two integer work arrays of length nx each.
c
      double precision a(*),x(*),y(*),f(*),fs(*)
      integer ijk(node,*), nodcode(*),ia(*),ja(*),iwk(*),jwk(*)
      external xyk

      ierr = 0
c
c     take into boundary conditions to remove boundary nodes.
c
      call bound (nx,nelx,ijk,nodcode,node,nint,jwk,
     *     x,y,f,iwk)
c
c     assemble the matrix
c
      call assmbo (nx,nelx,node,ijk,nodcode,x,y,
     *     a,ja,ia,f,iwk,jwk,ierr,xyk)
c
c     if applicable (job .eq. 1) get heat source function
c
      indic = 1
      if (job .eq. 1)
     *     call hsourc (indic,nx,nelx,node,x,y,ijk,fs,f)
c
c     call diric for Dirichlet conditions
c
      call diric(nx,nint,a,ja,ia,f)
c     done
      return
      end
