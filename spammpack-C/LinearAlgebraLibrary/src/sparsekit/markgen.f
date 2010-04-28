      subroutine markgen (m, n, a, ja, ia)

c*********************************************************************72
c
cc MARKGEN is a matrix generator for a Markov random walk on a triang. grid
c
c this subroutine generates a test matrix that models a random
c walk on a triangular grid. This test example was used by
c G. W. Stewart ["{SRRIT} - a FORTRAN subroutine to calculate the
c dominant invariant subspaces of a real matrix",
c Tech. report. TR-514, University of Maryland (1978).] and in a few
c papers on eigenvalue problems by Y. Saad [see e.g. LAA, vol. 34,
c pp. 269-295 (1980) ]. These matrices provide reasonably easy
c test problems for eigenvalue algorithms. The transpose of the
c matrix  is stochastic and so it is known that one is an exact
c eigenvalue. One seeks the eigenvector of the transpose associated
c with the eigenvalue unity. The problem is to calculate the
c steady state probability distribution of the system, which is
c the eigevector associated with the eigenvalue one and scaled in
c such a way that the sum all the components is equal to one.
c
c parameters
c
c on entry :
c
c m     = integer. number of points in each direction.
c
c on return:
c
c n     = integer. The dimension of the matrix. (In fact n is known
c         to be equal to (m(m+1))/2      )
c a,
c ja,
c ia    = the matrix stored in CSR format.
c
c
c Notes: 1) the code will actually compute the transpose of the
c stochastic matrix that contains the transition probibilities.
c        2) It should also be possible to have a matrix generator
c with an additional parameter (basically redefining `half' below
c to be another parameter and changing the rest accordingly, but
c this is not as simple as it sounds). This is not likely to provide
c any more interesting matrices.
c
      double precision a(*), cst, pd, pu, half
      integer ja(*), ia(*)
c
      cst = 0.5/dble(m-1)
c
c  ix counts the grid point (natural ordering used), i.e.,
c  the row number of the matrix.
c
      ix = 0
      jax = 1
      ia(1) = jax
c
c     sweep y coordinates
c
      do 20 i=1,m
         jmax = m-i+1
c
c     sweep x coordinates
c
         do 10 j=1,jmax
            ix = ix + 1
            if (j .eq. jmax) goto 2
            pd = cst*dble(i+j-1)
c
c     north
c
            a(jax) = pd
            if (i.eq. 1) a(jax) = a(jax)+pd
            ja(jax) =  ix + 1
            jax = jax+1
c     east
            a(jax) = pd
            if (j .eq. 1) a(jax) = a(jax)+pd
            ja(jax) = ix + jmax
            jax = jax+1
c     south
 2          pu = 0.5 - cst*dble(i+j-3)
            if ( j .gt. 1) then
               a(jax) = pu
               ja(jax) = ix-1
               jax = jax+1
            end if
c     west
            if ( i .gt. 1) then
               a(jax) = pu
               ja(jax) = ix - jmax - 1
               jax = jax+1
            end if
            ia(ix+1) = jax
 10      continue
 20   continue
      n = ix
      return
      end
