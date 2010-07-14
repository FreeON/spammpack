      subroutine amubdg(nrow,ncol,ncolb,ja,ia,jb,ib,ndegr,nnz,iw)

c*********************************************************************72
c
cc AMUBDG gets the number of nonzero elements in each row of A * B.
c
c on entry:
c
c nrow  = integer.  row dimension of matrix A
c ncol  = integer.  column dimension of matrix A = row dimension of
c                   matrix B.
c ncolb = integer. the colum dimension of the matrix B.
c
c ja, ia= row structure of input matrix A: ja = column indices of
c         the nonzero elements of A stored by rows.
c         ia = pointer to beginning of each row  in ja.
c
c jb, ib= row structure of input matrix B: jb = column indices of
c         the nonzero elements of A stored by rows.
c         ib = pointer to beginning of each row  in jb.
c
c on return:
c
c ndegr      = integer array of length nrow containing the degrees (i.e.,
c         the number of nonzeros in  each row of the matrix A * B
c
c nnz   = total number of nonzero elements found in A * B
c
c work arrays:
c
c iw      = integer work array of length ncolb.
c
      integer ja(*),jb(*),ia(nrow+1),ib(ncol+1),ndegr(nrow),iw(ncolb)

      do 1 k=1, ncolb
         iw(k) = 0
 1    continue

      do 2 k=1, nrow
         ndegr(k) = 0
 2    continue
c
c     method used: Transp(A) * A = sum [over i=1, nrow]  a(i)^T a(i)
c     where a(i) = i-th row of  A. We must be careful not to add  the
c     elements already accounted for.
c
c
      do 7 ii=1,nrow
c
c     for each row of A
c
         ldg = 0
c
c    end-of-linked list
c
         last = -1
         do 6 j = ia(ii),ia(ii+1)-1
c
c     row number to be added:
c
            jr = ja(j)
            do 5 k=ib(jr),ib(jr+1)-1
               jc = jb(k)
               if (iw(jc) .eq. 0) then
c
c     add one element to the linked list
c
                  ldg = ldg + 1
                  iw(jc) = last
                  last = jc
               end if
 5          continue
 6       continue
         ndegr(ii) = ldg
c
c     reset iw to zero
c
         do 61 k=1,ldg
            j = iw(last)
            iw(last) = 0
            last = j
 61      continue

 7    continue

      nnz = 0
      do 8 ii=1, nrow
         nnz = nnz+ndegr(ii)
 8    continue

      return
      end
