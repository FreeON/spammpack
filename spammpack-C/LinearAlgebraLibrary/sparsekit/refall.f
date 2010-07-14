      subroutine refall(nx, nelx,ijk,node,ndeg,x,y,
     *              ichild,iparnts,nodcode,nxmax,nelmax,ierr)

c*********************************************************************72
c
cc REFALL refines a finite element grid using triangular elements.
c
c REFALL refines a finite element grid using triangular elements.
c uses mid points to refine all the elements of the grid.
c
c nx      = number of nodes at input
c nelx      = number of elements at input
c ijk      = connectivity matrix: for node k, ijk(*,k) point to the
c         nodes of element k.
c ndeg      = first dimension of array ichild which is at least as large
c         as the max degree of each node
c x,y   = double precision arrays containing the x(*) and y(*) coordinates
c        resp. of the nodes.
c ichild= list of the children of a node: ichild(1,k) stores
c         the position in ichild(*,k)  of the last child so far.
c         (local use)
c iparnts= list of the 2 parents of each node.
c         (local use)
c nodcode= boundary information list for each node with the
c         following meaning:
c      nodcode(i) = 0 -->  node i is internal
c      nodcode(i) = 1 -->  node i is a boundary but not a corner point
c      nodcode(i) = 2 -->  node i is a corner point.
c corner elements are used only to generate the grid by refinement
c since they do not  correspond to real elements.
c nxmax  = maximum number of nodes allowed. If during the algorithm
c          the number of nodes being created exceeds nxmax then
c         refall  quits without modifying the (x,y) xoordinates
c         and nx, nelx. ijk is modified. Also ierr is set to 1.
c nelmax = same as above for number of elements allowed. See ierr..
c ierr       = error message:
c         0 --> normal return
c         1 --> refall quit because nxmax  was exceeded.
c         2 --> refall quit because nelmax was exceeded.
c
       implicit double precision  (a-h,o-z)
       integer ichild(ndeg,1),iparnts(2,nx),ijk(node,1),nodcode(nx)
       integer midnode(10),inod(10)
       double precision  x(1),y(1)
c
c inilitialize lists of children and parents --
c data structure is as follows
c ichild(1,k) stores the position of last child of node k so far in list
c ichild(j,k) , j .ge. 2 = list of children of node k.
c iparnts(1,k) and iparnts(2,k) are the two parents of node k.
c
c  do a first check :
      if (nx .ge. nxmax) goto 800
      if (nelx .ge. nelmax) goto 900
c  initialize
        do 1 k=1,nx
      do 2 j=2,ndeg
       ichild(j,k) = 0
 2      continue
      ichild(1,k) = 1
      iparnts(1,k)= 0
      iparnts(2,k)= 0
 1      continue
c  initialize nelxnew and nxnew
      nelxnew = nelx
      nxnew   = nx
      ierr    = 0
c
c main loop: scan all elements
c
c      do 100 nel = nelx,1,-1
      do 100 nel = 1, nelx
c note : interesting question which order is best for parallelism?
c alternative order: do 100 nel = nelx, 1, -1
c
c  unpack nodes of element
      do 101 i=1,node
      inod(i) = ijk(i,nel)
c convention: node after last node = first node.
      inod(node+i) = inod(i)
      midnode(i) = 0
 101      continue
c
c for each new potential node determine if it has already been
c numbered. a potential node is the middle of any two nodes ..
c
      do 80 ii=1,node
              k1 = inod(ii)
              k2 = inod(ii+1)
c  test for current pair :
      last = ichild(1,k1)
      do 21 k=2,last
            jchild = ichild(k,k1)
            ipar1 = iparnts(1,jchild)
            ipar2 = iparnts(2,jchild)
              if( (ipar1 .eq. k1 .and. ipar2 .eq. k2) .or.
     *                (ipar2 .eq. k1 .and. ipar1 .eq. k2)) then
c node has already been created and numbered ....
          midnode(ii) = jchild
c therefore it must be an internal node
          nodcode(jchild) = 0
c  and no new node to create.
          goto 80
              end if
 21      continue
c
c else  create a new node
c
      nxnew = nxnew + 1
      if (nxnew .gt. nxmax) goto 800

      x(nxnew) = (x(k1) + x(k2))*0.5
      y(nxnew) = (y(k1) + y(k2))*0.5
      midnode(ii) = nxnew
c
c update nodcode information -- normally min0(nodcode(k1),nodcode(k2))
c
       nodcode(nxnew) = min0(1,nodcode(k1),nodcode(k2))
c
c update parents and children's lists
c
      iparnts(1,nxnew) = k1
      iparnts(2,nxnew) = k2
c
      last = last+1
      ichild(last,k1) = nxnew
      ichild(1,k1) = last

      last = ichild(1,k2)+1
      ichild(last,k2) = nxnew
      ichild(1,k2) = last
c
 80     continue
c
c  replace current element by new one
c
      do 81 i=1,node
      jnod = midnode(i)
        ijk(i,nel) = jnod
 81     continue
c  create new elements
      do 82 ii=1, node
      nelxnew = nelxnew+1
      if (nelxnew .gt. nelmax) goto 900
      ijk(1,nelxnew) = inod(ii)
      k = ii
      do 82 jj=2,node
      ijk(jj,nelxnew) = midnode(k)
       k = k+2
      if (k .gt. node) k =  k-node
 82      continue
c  done !
 100      continue
      nx = nxnew
      nelx = nelxnew
        return
 800       ierr = 1
      return
 900    ierr = 2
       return
      end
