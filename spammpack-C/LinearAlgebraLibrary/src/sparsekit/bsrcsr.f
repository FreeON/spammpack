      subroutine bsrcsr (n, nblk, na, a, ja, ia, ao, jao, iao)

c*********************************************************************72
c
cc BSRCSR converts Block Sparse Row to Compressed Sparse Row (CSR) format.
c
c this routine converts a matrix stored in block-reduced
c a, ja, ia format to the general sparse row a, ja, ia format.
c A matrix that has a block structure is a matrix whose entries
c are blocks of the same size nblk (e.g. 3 x 3). Then it is often
c preferred to work with the reduced graph of the matrix, i.e.,
c Instead of storing one element at a time one can store the whole
c block. In this storage scheme a row of the array a will
c hold the nblk**2 entries of a block.
c
c on entry:
c
c n      = integer, the actual row dimension of the matrix.
c nblk  = integer equal to the dimension of each block.
c         nblk must divide n.
c na      = first dimension of array a as declared in calling program
c a      = double precision array containing the values of the matrix. For details
c         on the format see below. Each row of a contains the nblk x nblk
c         block matrix unpacked column-wise (this allows the user to
c         declare the array a as a(na,nblk,nblk) on entry if desired).
c         the block rows are stored in sequence just as for the compressed
c         sparse row format.
c ja      = integer array of length n/nblk. ja(k) contains the column index
c         of the leading element, i.e., the element (1,1) of the block
c         that is held in the row a(k,*) of the value array.
c ia    = integer array of length n/nblk+1. ia(i) points to the beginning
c        of block row number i in the arrays a and ja.
c
c on return:
c
c ao, jao,
c     iao = matrix stored in compressed sparse row format.
c
c Notes: this code is not in place.
c
c
c general picture: (nblk = 2)
c     --- A ---                                --- JA --  -- IA --
c A=  x x x x   1st block in block row 1           x         x
c     x x x x  2-nd block in block row 1           x
c     . . . .                                      .
c     x x x x  last block in block row 1           x
c     -------                                     ---
c     x x x x   1st block in block row 2           x          x
c     x x x x  2-nd block in block row 2           x
c     . . . .                                      x
c     x x x x   last block in block row 2          x
c     -------                                     ---
c     .......                                     ...         .
c     -------                                     ---
c     x x x x   1st block in block row n/nblk      x          x
c     x x x x  2-nd block in block row n/nblk      x
c     . . . .                                      x
c     x x x x  last block in block row n/nblk      x
c     -------                                     ---
c                                               end + 1       x
c
c
c example  with nblk = 2:
c
c
c             1   2   0   0   3   4
c             5   6   0   0   7   8
c             0   0   9  10  11  12
c             0   0  13  14  15  16
c             17 18   0   0   0   0
c             22 23   0   0   0   0
c THEN:
c
c  ---- A ----                                     -- JA --   -- IA --
c
c  1   5   2  6  Block row 1 (2 block matrices)      | 1  <--- | 1
c  3   7   4  8                                      | 5       |
c  ------------                                      | --      |
c  9  13  10 14  block row 2 (2 block matrices)      | 3  <--- | 3
c 11  15  12 16                                      | 5       |
c  ------------                                      | --      |
c 17  22  18 23  Block row 3 (1 block matrix)        | 1  <--- | 5
c  ------------                                      | --      |
c                                                   end+1 <--- | 6
c
c JA  =  1  5 | 3  5 | 1       column numbers of (1,1) entries of blocks
c IA  =  1      3      5  6    pointers to beginnings of BLOCK-rows
c
c
c get ia, ja data structure for output matrix
c
      double precision a(na,*), ao(*)
      integer ia(*), ja(*), jao(*), iao(n+1)

      nr = n/nblk
      do 1 k=1,n+1
         iao(k) = 0
 1    continue

      irow = 0
      krow = 1
      do 2 ii=1, nr
c     nr is the dimension of the reduced matrix.
         i1 = ia(ii)
         i2 = ia(ii+1)-1
c     create nblk rows for each k
         do 23 i=1,nblk
            do 21 k=i1, i2
               jst = ja(k)-1
               do 22  j=1,nblk
                  ij = (j-1)*nblk + i
                  ao(krow) = a(k,ij)
                  jao(krow) = jst+j
                  krow = krow+1
 22            continue
 21          continue
          iao(irow+i) = krow
 23      continue
         irow = irow + nblk
 2    continue
      do 3 jj=1,n
         j = n-jj+1
         iao(j+1)=iao(j)
 3    continue
      iao(1) = 1
      return
      end
