      subroutine blkchk (nrow,ja,ia,nblk,imsg)

c*********************************************************************72
c
cc BLKCHK checks whether the input matrix is a block matrix.
c
c BLKCHK checks whether the input matrix is a block
c matrix with block size of nblk. A block matrix is one which is
c comprised of small square dense blocks. If there are zero
c elements within the square blocks and the data structure
c takes them into account then blkchk may fail to find the
c correct block size.
c
c on entry
c
c nrow      = integer equal to the row dimension of the matrix.
c ja    = integer array containing the column indices of the entries
c         nonzero entries of the matrix stored by row.
c ia    = integer array of length nrow + 1 containing the pointers
c         beginning of each row in array ja.
c
c nblk  = integer containing the value of nblk to be checked.
c
c on return
c
c
c imsg  = integer containing a message  with the following meaning.
c          imsg = 0 means that the output value of nblk is a correct
c                   block size. nblk .lt. 0 means nblk not correct
c                   block size.
c          imsg = -1 : nblk does not divide nrow
c          imsg = -2 : a starting element in a row is at wrong position
c             (j .ne. mult*nblk +1 )
c          imsg = -3 : nblk does divide a row length -
c          imsg = -4 : an element is isolated outside a block or
c             two rows in same group have different lengths
c
c           Y. Saad, Sep. 21 1989                                      c
c
      integer ia(nrow+1),ja(*)
c
c first part of code will find candidate block sizes.
c this is not guaranteed to work . so a check is done at the end
c the criterion used here is a simple one:
c scan rows and determine groups of rows that have the same length
c and such that the first column number and the last column number
c are identical.
c
      imsg = 0
      if (nblk .le. 1) return
      nr = nrow/nblk
      if (nr*nblk .ne. nrow) goto 101
c   main loop
      irow = 1
      do 20 ii=1, nr
c     i1= starting position for group of nblk rows in original matrix
         i1 = ia(irow)
         j2 = i1
c     lena = length of each row in that group  in the original matrix
         lena = ia(irow+1)-i1
c     len = length of each block-row in that group in the output matrix
         len = lena/nblk
         if (len* nblk .ne. lena) goto 103
c
c     for each row
c
         do 6 i = 1, nblk
            irow = irow + 1
            if (ia(irow)-ia(irow-1) .ne. lena ) goto 104
c
c     for each block
c
            do 7 k=0, len-1
               jstart = ja(i1+nblk*k)-1
               if ( (jstart/nblk)*nblk .ne. jstart) goto 102
c
c     for each column
c
               do 5 j=1, nblk
                  if (jstart+j .ne. ja(j2) )  goto 104
                  j2 = j2+1
 5             continue
 7          continue
 6       continue
 20   continue
c     went through all loops successfully:
      return
 101  imsg = -1
      return
 102  imsg = -2
      return
 103  imsg = -3
      return
 104  imsg = -4
      return
      end
