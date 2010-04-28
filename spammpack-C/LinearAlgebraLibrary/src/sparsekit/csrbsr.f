      subroutine csrbsr (n,nblk,na,a,ja,ia,ao,jao,iao)

c*********************************************************************72
c
cc CSRBSR converts Compressed Sparse Row to Block Sparse Row.
c
c CSRBSR does the reverse of bsrcsr. It converts
c a matrix stored in a general compressed a, ja, ia format into a
c a block reduced matrix a(*,*),ja(*),ia(*) format. The code
c assumes that the original matrix is indeed a block format
c and that the elements are ordered in such a way that their
c column numbers are increasing. (This can be achieved
c by transposing a, ja, ia twice, putting the resulting matrix
c into a, ja, ia).
c
c See routine bsrcsr for more details on data structure for blocked
c matrices. The input matrix is a, ja, ia (in compressed format) and
c the output matrix is the matrix ao, jao, iao in block-reduced
c format.
c
c on entry:
c
c n      = integer, the actual row dimension of the matrix.
c nblk  = integer equal to the dimension of each block.
c         nblk must divide n.
c na      = first dimension of array a as declared in calling program
c
c a, ja,
c    ia = input matrix stored in compressed sparse row format.
c
c on return:
c
c ao    = double precision array containing the values of the matrix. For details
c         on the format see below. Each row of a contains the nblk x nblk
c         block matrix unpacked column-wise (this allows the user to
c         declare the array a as a(na,nblk,nblk) on entry if desired).
c         the block rows are stored in sequence just as for the compressed
c         sparse row format.
c jao   = integer array of length n/nblk. ja(k) contains the column index
c         of the leading element, i.e., the element (1,1) of the block
c         that is held in the row a(k,*) of the value array.
c iao   = integer array of length n/nblk+1. ia(i) points to the beginning
c        of block row number i in the arrays a and ja.
c
c Notes:
c  1) this code is not in place.
c        2) see routine bsrcsr for details on data sctructure for
c           block sparse row format.
c        3) The routine assumes that  the input matrix has been
c           sorted in such a way that the column indices are always
c           in increasing order for the same row.
c           for all k "in the SAME ROW."
c        4) THERE IS NO CHECKING AS TO WHETHER the input is correct.
c           it is recommended to use the routine blchk to check
c           if the matrix is a block-matrix before calling csrbsr.
c
      double precision a(*),ao(na,*)
      integer ia(n+1),ja(*),jao(*),iao(*)
c
c     nr is the dimension of the reduced matrix.
c
      nr = n/nblk
      iao(1) = 1
      ibrow = 1
      irow  = 1
c  main loop
      do 2 ii=1, nr
c     i1= starting position for group of nblk rows in original matrix
         i1 = ia(irow)
c     lena = length of each row in that group  in the original matrix
         lena = ia(irow+1)-i1
c     len = length of each block-row in that group in the output matrix
         len = lena/nblk
         k1 = iao(ibrow)
c  copy the real values of A
c     for each block
         do 7 k=0, len-1
c     store column positions of the (1,1) elements of each block
            jao(k1+k) = ja(i1+nblk*k)
c     for each column
            do 5 j=1, nblk
               j1 = (j-1)*nblk
               j2 = i1+k*nblk+j-1
c     for each row
               do 6 i = 1, nblk
                  ao(k1+k,j1+i) = a(j2+(i-1)*lena)
 6             continue
 5          continue
 7       continue
c     done with a whole block row. now update iao(*), ibrow and irow
         iao(ibrow+1) = iao(ibrow)+len
         ibrow = ibrow + 1
         irow  = irow + nblk
 2    continue
      return
      end
