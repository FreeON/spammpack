      subroutine checkref(nx,nelx,ijk,node,nodcode,
     *               nbound,  nxnew,nelxnew)

c*********************************************************************72
c
cc CHECKREF returns the expected number of new nodes and elements.
c
c returns the expected the new number of nodes and
c elemnts if refall is applied to current grid once.
c
c nx      = number of nodes at input
c nelx      = number of elements at input
c ijk      = connectivity matrix: for node k, ijk(*,k) point to the
c         nodes of element k.
c nbound  = number of boundary points on entry - enter zero if
c           unknown
c
c nodcode= boundary information list for each node with the
c         following meaning:
c      nodcode(i) = 0 -->  node i is internal
c      nodcode(i) = 1 -->  node i is a boundary but not a corner point
c      nodcode(i) = 2 -->  node i is a corner point.
c
c nxnew  = new number of nodes if refall were to be applied
c nelxnew = same for nelx.
c
       integer ijk(node,1),nodcode(nx)

      nelxnew = nelx*4
c
c count the number of boundary nodes
c
      if (nbound .ne. 0) goto 2
      do 1 j=1, nx
      if (nodcode(j) .ge. 1) nbound = nbound+1
 1      continue
c number of edges=[3*(number of elmts) + number of bound nodes ]/ 2
 2      continue
      nxnew = nx + (3*nelx+nbound)/2
      nbound = 2*nbound
      return
      end
