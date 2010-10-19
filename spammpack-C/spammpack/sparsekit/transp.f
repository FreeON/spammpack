      subroutine transp(nrow,ncol,a,ja,ia,iwk,ierr)

c*********************************************************************72
c
cc TRANSP carries out in-place transposition routine.
c
c this subroutine transposes a matrix stored in compressed sparse row
c format. the transposition is done in place in that the arrays a,ja,ia
c of the transpose are overwritten onto the original arrays.
c
c on entry:
c
c nrow      = integer. The row dimension of A.
c ncol      = integer. The column dimension of A.
c a      = double precision array of size nnz (number of nonzero elements in A).
c         containing the nonzero elements
c ja      = integer array of length nnz containing the column positions
c         of the corresponding elements in a.
c ia      = integer of size n+1, where n = max(nrow,ncol). On entry
c         ia(k) contains the position in a,ja of  the beginning of
c         the k-th row.
c
c iwk      = integer work array of same length as ja.
c
c on return:
c
c
c ncol      = actual row dimension of the transpose of the input matrix.
c         Note that this may be .le. the input value for ncol, in
c         case some of the last columns of the input matrix are zero
c         columns. In the case where the actual number of rows found
c         in transp(A) exceeds the input value of ncol, transp will
c         return without completing the transposition. see ierr.
c a,
c ja,
c ia      = contains the transposed matrix in compressed sparse
c         row format. The row dimension of a, ja, ia is now ncol.
c
c ierr      = integer. error message. If the number of rows for the
c         transposed matrix exceeds the input value of ncol,
c         then ierr is  set to that number and transp quits.
c         Otherwise ierr is set to 0 (normal return).
c
c Note:
c      1) If you do not need the transposition to be done in place
c         it is preferrable to use the conversion routine csrcsc
c         (see conversion routines in formats).
c      2) the entries of the output matrix are not sorted (the column
c         indices in each are not in increasing order) use csrcsc
c         if you want them sorted.
c
c           Y. Saad, Sep. 21 1989                                      c
c  modified Oct. 11, 1989.                                             c
c
      integer nrow, ncol, ia(*), ja(*), iwk(*), ierr
      double precision a(*)
      double precision t, t1
      ierr = 0
      nnz = ia(nrow+1)-1
c
c     determine column dimension
c
      jcol = 0

      do k=1, nnz
         jcol = max(jcol,ja(k))
      end do

      if (jcol .gt. ncol) then
         ierr = jcol
         return
      end if
c
c     convert to coordinate format. use iwk for row indices.
c
      ncol = jcol
c
      do 3 i=1,nrow
         do 2 k=ia(i),ia(i+1)-1
            iwk(k) = i
 2       continue
 3    continue
c     find pointer array for transpose.
      do 35 i=1,ncol+1
         ia(i) = 0
 35   continue
      do 4 k=1,nnz
         i = ja(k)
         ia(i+1) = ia(i+1)+1
 4    continue
      ia(1) = 1

      do i=1,ncol
         ia(i+1) = ia(i) + ia(i+1)
      end do
c
c  loop for a cycle in chasing process.
c
      init = 1
      k = 0
 5    t = a(init)
      i = ja(init)
      j = iwk(init)
      iwk(init) = -1

 6    k = k+1
c
c     current row number is i.  determine  where to go.
c
      l = ia(i)
c
c     save the chased element.
c
      t1 = a(l)
      inext = ja(l)
c
c     then occupy its location.
c
      a(l)  = t
      ja(l) = j
c
c     update pointer information for next element to be put in row i.
c
      ia(i) = l+1
c
c     determine  next element to be chased
c
      if (iwk(l) .lt. 0) goto 65
      t = t1
      i = inext
      j = iwk(l)
      iwk(l) = -1
      if (k .lt. nnz) goto 6
      goto 70
 65   init = init+1
      if (init .gt. nnz) goto 70
      if (iwk(init) .lt. 0) goto 65
c
c     restart chasing
c
      goto 5
 70   continue
      do 80 i=ncol,1,-1
         ia(i+1) = ia(i)
 80   continue
      ia(1) = 1

      return
      end
