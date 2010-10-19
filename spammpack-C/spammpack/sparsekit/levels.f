      subroutine levels (n, jal, ial, nlev, lev, ilev, levnum)

c*********************************************************************72
c
cc LEVELS gets the level structure of a lower triangular matrix.
c
c levels gets the level structure of a lower triangular matrix
c for level scheduling in the parallel solution of triangular systems
c strict lower matrices (e.g. unit) as well matrices with their main
c diagonal are accepted.
c
c on entry:
c
c n        = integer. The row dimension of the matrix
c jal, ial =
c
c on return:
c
c nlev     = integer. number of levels found
c lev      = integer array of length n containing the level
c            scheduling permutation.
c ilev     = integer array. pointer to beginning of levels in lev.
c            the numbers lev(i) to lev(i+1)-1 contain the row numbers
c            that belong to level number i, in the level scheduling
c            ordering. The equations of the same level can be solved
c            in parallel, once those of all the previous levels have
c            been solved.
c work arrays:
c
c levnum   = integer array of length n (containing the level numbers
c            of each unknown on return)
      integer jal(*),ial(*), levnum(*), ilev(*), lev(*)

      do 10 i = 1, n
         levnum(i) = 0
 10   continue
c
c     compute level of each node --
c
      nlev = 0
      do 20 i = 1, n
         levi = 0
         do 15 j = ial(i), ial(i+1) - 1
            levi = max (levi, levnum(jal(j)))
 15      continue
         levi = levi+1
         levnum(i) = levi
         nlev = max(nlev,levi)
 20   continue
c  set data structure.
      do 21 j=1, nlev+1
         ilev(j) = 0
 21   continue
c  count  number   of elements in each level.
      do 22 j=1, n
         i = levnum(j)+1
         ilev(i) = ilev(i)+1
 22   continue
c  set up pointer for  each  level.
      ilev(1) = 1
      do 23 j=1, nlev
         ilev(j+1) = ilev(j)+ilev(j+1)
 23   continue

c  determine elements of each level.
      do 30 j=1,n
         i = levnum(j)
         lev(ilev(i)) = j
         ilev(i) = ilev(i)+1
 30   continue
c     reset pointers backwards
      do 35 j=nlev, 1, -1
         ilev(j+1) = ilev(j)
 35   continue
      return
      end
