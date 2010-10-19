      subroutine aplbdg (nrow,ncol,ja,ia,jb,ib,ndegr,nnz,iw)
      integer ja(*),jb(*),ia(nrow+1),ib(nrow+1),iw(ncol),ndegr(nrow)

c*********************************************************************72
c
cc APLBDG gets the number of nonzero elements in each row of A + B.
c
c on entry:
c
c nrow      = integer. The row dimension of A and B
c ncol  = integer. The column dimension of A and B.
c
c a,
c ja,
c ia   = Matrix A in compressed sparse row format.
c
c b,
c jb,
c ib      =  Matrix B in compressed sparse row format.
c
c on return:
c
c ndegr      = integer array of length nrow containing the degrees (i.e.,
c         the number of nonzeros in  each row of the matrix A + B.
c
c nnz   = total number of nonzero elements found in A * B
c
c work arrays:
c
c iw      = integer work array of length equal to ncol.
c
      do 1 k=1, ncol
         iw(k) = 0
 1    continue

      do 2 k=1, nrow
         ndegr(k) = 0
 2    continue

      do 7 ii=1,nrow
         ldg = 0
c
c    end-of-linked list
c
         last = -1
c
c     row of A
c
         do 5 j = ia(ii),ia(ii+1)-1
            jr = ja(j)
c
c     add element to the linked list
c
            ldg = ldg + 1
            iw(jr) = last
            last = jr
 5       continue
c
c     row of B
c
         do 6 j=ib(ii),ib(ii+1)-1
            jc = jb(j)
            if (iw(jc) .eq. 0) then
c
c     add one element to the linked list
c
               ldg = ldg + 1
               iw(jc) = last
               last = jc
            end if
 6       continue
c     done with row ii.
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
c
      nnz = 0
      do 8 ii=1, nrow
         nnz = nnz+ndegr(ii)
 8    continue
      return
      end
