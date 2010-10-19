      subroutine blkfnd (nrow,ja,ia,nblk)

c*********************************************************************72
c
cc BLKFND determines the block structure of a matrix.
c
c This routine attemptps to determine whether or not  the input
c matrix has a block structure and finds the blocks size
c if it does. A block matrix is one which is
c comprised of small square dense blocks. If there are zero
c elements within the square blocks and the original data structure
c takes these zeros into account then blkchk may fail to find the
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
c nblk  = integer containing the assumed value of nblk if job = 0
c
c on return
c
c nblk  = integer containing the value found for nblk when job = 1.
c         if imsg .ne. 0 this value is meaningless however.
c
c
c           Y. Saad, Sep. 21 1989                                      c
c
      integer ia(nrow+1),ja(*)
c
c first part of code will find candidate block sizes.
c criterion used here is a simple one: scan rows and  determine groups
c of rows that have the same length and such that the first column
c number and the last column number are identical.
c
      minlen = ia(2)-ia(1)
      irow   = 1
      do 1 i=2,nrow
         len = ia(i+1)-ia(i)
         if (len .lt. minlen) then
            minlen = len
            irow = i
         end if
 1    continue
c
c  candidates are all dividers of minlen
c
      nblk = 1
      if (minlen .le. 1) return

      do 99 iblk = minlen, 1, -1
         if (mod(minlen,iblk) .ne. 0) goto 99
         len = ia(2) - ia(1)
         len0 = len
         jfirst = ja(1)
         jlast = ja(ia(2)-1)
         do 10 jrow = irow+1,irow+nblk-1
            i1 = ia(jrow)
            i2 = ia(jrow+1)-1
            len = i2+1-i1
            jf = ja(i1)
            jl = ja(i2)
            if (len .ne. len0 .or. jf .ne. jfirst .or.
     *           jl .ne. jlast) goto 99
 10      continue
c
c     check for this candidate
c
         call blkchk (nrow,ja,ia,iblk,imsg)
         if (imsg .eq. 0) then
c
c     block size found
c
            nblk = iblk
            return
         end if
 99   continue
      end
