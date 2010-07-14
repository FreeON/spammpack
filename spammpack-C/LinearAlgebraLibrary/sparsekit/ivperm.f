      subroutine ivperm (n, ix, perm)

c*********************************************************************72
c
cc IVPERM performs an in-place permutation of an integer vector.
c
c this subroutine performs an in-place permutation of an integer vector
c ix according to the permutation array perm(*), i.e., on return,
c the vector x satisfies,
c
c      ix(perm(j)) :== ix(j), j=1,2,.., n
c
c on entry:
c
c n       = length of vector x.
c perm       = integer array of length n containing the permutation  array.
c ix      = input vector
c
c on return:
c
c ix      = vector x permuted according to ix(perm(*)) :=  ix(*)
c
c
c           Y. Saad, Sep. 21 1989
c
      integer n, perm(n), ix(n)
      integer tmp, tmp1

      init      = 1
      tmp      = ix(init)
      ii        = perm(init)
      perm(init)= -perm(init)
      k         = 0
c
c loop
c
 6    k = k+1
c
c save the chased element --
c
      tmp1        = ix(ii)
      ix(ii)     = tmp
      next        = perm(ii)
      if (next .lt. 0 ) goto 65
c
c test for end
c
      if (k .gt. n) goto 101
      tmp       = tmp1
      perm(ii)  = - perm(ii)
      ii        = next
c
c end loop
c
      goto 6
c
c reinitilaize cycle --
c
 65   init      = init+1
      if (init .gt. n) goto 101
      if (perm(init) .lt. 0) goto 65
      tmp      = ix(init)
      ii      = perm(init)
      perm(init)=-perm(init)
      goto 6

 101  continue
      do 200 j=1, n
         perm(j) = -perm(j)
 200  continue

      return
      end
