      subroutine coicsr (n,nnz,job,a,ja,ia,iwk)

c*********************************************************************72
c
cc COICSR converts COO to CSR in place.
c
c IN-PLACE coo-csr conversion routine.
c
c COICSR converts a matrix stored in coordinate format into
c the csr format. The conversion is done in place in that the arrays
c a,ja,ia of the result are overwritten onto the original arrays.
c
c on entry:
c
c n      = integer. row dimension of A.
c nnz      = integer. number of nonzero elements in A.
c job   = integer. Job indicator. when job=1, the real values in a are
c         filled. Otherwise a is not touched and the structure of the
c         array only (i.e. ja, ia)  is obtained.
c a      = double precision array of size nnz (number of nonzero elements in A)
c         containing the nonzero elements
c ja      = integer array of length nnz containing the column positions
c         of the corresponding elements in a.
c ia      = integer array of length nnz containing the row positions
c         of the corresponding elements in a.
c iwk      = integer work array of length n.
c on return:
c
c a
c ja
c ia      = contains the compressed sparse row data structure for the
c         resulting matrix.
c Note:
c
c         the entries of the output matrix are not sorted (the column
c         indices in each are not in increasing order) use coocsr
c         if you want them sorted.
c
c  Coded by Y. Saad, Sep. 26 1989
c
      integer ia(nnz),ja(nnz),iwk(n)
      double precision a(*)
      double precision t,tnext
      logical values

      values = (job .eq. 1)
c
c find pointer array for resulting matrix.
c
      do 35 i=1,n+1
         iwk(i) = 0
 35   continue
      do 4 k=1,nnz
         i = ia(k)
         iwk(i+1) = iwk(i+1)+1
 4    continue

      iwk(1) = 1
      do 44 i=2,n
         iwk(i) = iwk(i-1) + iwk(i)
 44   continue
c
c  loop for a cycle in chasing process.
c
      init = 1
      k = 0
 5    if (values) t = a(init)
      i = ia(init)
      j = ja(init)
      ia(init) = -1

 6    k = k+1
c
c  current row number is i.  determine  where to go.
c
      ipos = iwk(i)
c
c  save the chased element.
c
      if (values) tnext = a(ipos)
      inext = ia(ipos)
      jnext = ja(ipos)
c
c  then occupy its location.
c
      if (values) a(ipos)  = t
      ja(ipos) = j
c
c  update pointer information for next element to come in row i.
c
      iwk(i) = ipos+1
c
c  determine  next element to be chased,
c
      if (ia(ipos) .lt. 0) goto 65
      t = tnext
      i = inext
      j = jnext
      ia(ipos) = -1
      if (k .lt. nnz) goto 6
      goto 70
 65   init = init+1
      if (init .gt. nnz) goto 70
      if (ia(init) .lt. 0) goto 65
c
c  restart chasing.
c
      goto 5
 70   do 80 i=1,n
         ia(i+1) = iwk(i)
 80   continue
      ia(1) = 1
      return
      end
