      subroutine csrbnd (n,a,ja,ia,job,abd,nabd,lowd,ml,mu,ierr)

c*********************************************************************72
c
cc CSRBND converts Compressed Sparse Row to Banded Linpack format.
c
c CSRBND converts a general sparse matrix stored in
c compressed sparse row format into the banded format. for the
c banded format,the Linpack conventions are assumed (see below).
c
c on entry:
c
c n      = integer,the actual row dimension of the matrix.
c
c a,
c ja,
c ia    = input matrix stored in compressed sparse row format.
c
c job      = integer. if job=1 then the values of the lower bandwith ml
c         and the upper bandwidth mu are determined internally.
c         otherwise it is assumed that the values of ml and mu
c         are the correct bandwidths on input. See ml and mu below.
c
c nabd  = integer. first dimension of array abd.
c
c lowd  = integer. this should be set to the row number in abd where
c         the lowest diagonal (leftmost) of A is located.
c         lowd should be  ( 1  .le.  lowd  .le. nabd).
c         if it is not known in advance what lowd should be
c         enter lowd = 0 and the default value lowd = ml+mu+1
c         will be chosen. Alternative: call routine getbwd from unary
c         first to detrermione ml and mu then define lowd accordingly.
c         (Note: the banded solvers in linpack use lowd=2*ml+mu+1. )
c
c ml      = integer. equal to the bandwidth of the strict lower part of A
c mu      = integer. equal to the bandwidth of the strict upper part of A
c         thus the total bandwidth of A is ml+mu+1.
c         if ml+mu+1 is found to be larger than lowd then an error
c         flag is raised (unless lowd = 0). see ierr.
c
c note:   ml and mu are assumed to have       the correct bandwidth values
c         as defined above if job is set to zero on entry.
c
c on return:
c
c
c abd   = double precision array of dimension abd(nabd,n).
c         on return contains the values of the matrix stored in
c         banded form. The j-th column of abd contains the elements
c         of the j-th column of  the original matrix comprised in the
c         band ( i in (j-ml,j+mu) ) with the lowest diagonal at
c         the bottom row (row lowd). See details below for this format.
c
c ml      = integer. equal to the bandwidth of the strict lower part of A
c mu      = integer. equal to the bandwidth of the strict upper part of A
c         if job=1 on entry then these two values are internally computed.
c
c lowd  = integer. row number in abd where the lowest diagonal
c         (leftmost) of A is located on return. In case lowd = 0
c         on return, then it is defined to ml+mu+1 on return and the
c         lowd will contain this value on return. `
c
c ierr  = integer. used for error messages. On return:
c         ierr .eq. 0  :means normal return
c         ierr .eq. -1 : means invalid value for lowd. (either .lt. 0
c         or larger than nabd).
c         ierr .eq. -2 : means that lowd is not large enough and as
c         result the matrix cannot be stored in array abd.
c         lowd should be at least ml+mu+1, where ml and mu are as
c         provided on output.
c
c Additional details on banded format.  (this closely follows the
c format used in linpack. may be useful for converting a matrix into
c this storage format in order to use the linpack  banded solvers).
c
c band storage format  for matrix abd
c uses ml+mu+1 rows of abd(nabd,*) to store the diagonals of
c a in rows of abd starting from the lowest (sub)-diagonal  which  is
c stored in row number lowd of abd. the minimum number of rows needed
c in abd is ml+mu+1, i.e., the minimum value for lowd is ml+mu+1. the
c j-th  column  of  abd contains the elements of the j-th column of a,
c from bottom to top: the element a(j+ml,j) is stored in  position
c abd(lowd,j), then a(j+ml-1,j) in position abd(lowd-1,j) and so on.
c Generally, the element a(j+k,j) of original matrix a is stored in
c position abd(lowd+k-ml,j), for k=ml,ml-1,..,0,-1, -mu.
c The first dimension nabd of abd must be .ge. lowd
c
c     example [from linpack ]:   if the original matrix is
c
c              11 12 13  0  0  0
c              21 22 23 24  0  0
c               0 32 33 34 35  0     original banded matrix
c               0  0 43 44 45 46
c               0  0  0 54 55 56
c               0  0  0  0 65 66
c
c then  n = 6, ml = 1, mu = 2. lowd should be .ge. 4 (=ml+mu+1)  and
c if lowd = 5 for example, abd  should be:
c
c untouched --> x  x  x  x  x  x
c               *  * 13 24 35 46
c               * 12 23 34 45 56    resulting abd matrix in banded
c              11 22 33 44 55 66    format
c  row lowd--> 21 32 43 54 65  *
c
c * = not used
c
      double precision a(*),abd(nabd,n)
      integer ia(n+1),ja(*)
c
c first determine ml and mu.
c
      ierr = 0

      if (job .eq. 1) call getbwd(n,a,ja,ia,ml,mu)
      m = ml+mu+1
      if (lowd .eq. 0) lowd = m
      if (m .gt. lowd)  ierr = -2
      if (lowd .gt. nabd .or. lowd .lt. 0) ierr = -1
      if (ierr .lt. 0) return

      do 15  i=1,m
         ii = lowd -i+1
         do 10 j=1,n
          abd(ii,j) = 0.0
 10      continue
 15   continue

      mdiag = lowd-ml
      do 30 i=1,n
         do 20 k=ia(i),ia(i+1)-1
            j = ja(k)
            abd(i-j+mdiag,j) = a(k)
 20      continue
 30   continue
      return
      end
