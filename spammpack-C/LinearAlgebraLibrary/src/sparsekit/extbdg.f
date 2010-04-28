      subroutine extbdg (n,a,ja,ia,bdiag,nblk,ao,jao,iao)

c*********************************************************************72
c
cc EXTBDG extracts the main diagonal blocks of a matrix.
c
c this subroutine extracts the main diagonal blocks of a
c matrix stored in compressed sparse row format and puts the result
c into the array bdiag and the remainder in ao,jao,iao.
c
c on entry:
c
c n      = integer. The row dimension of the matrix a.
c a,
c ja,
c ia    = matrix stored in csr format
c nblk  = dimension of each diagonal block. The diagonal blocks are
c         stored in compressed format rowwise,i.e.,we store in
c        succession the i nonzeros of the i-th row after those of
c        row number i-1..
c
c on return:
c
c bdiag = double precision array of size (n x nblk) containing the diagonal
c        blocks of A on return
c ao,
c jao,
c iao   = remainder of the matrix stored in csr format.
c
c           Y. Saad, Sep. 21 1989                                      c
c
      implicit double precision (a-h,o-z)
      double precision bdiag(*),a(*),ao(*)
      integer ia(*),ja(*),jao(*),iao(*)

      m = 1 + (n-1)/nblk
c this version is sequential -- there is a more parallel version
c that goes through the structure twice ....
      ltr =  ((nblk-1)*nblk)/2
      l = m * ltr
      do 1 i=1,l
         bdiag(i) = 0.0
 1    continue
      ko = 0
      kb = 1
      iao(1) = 1

      do 11 jj = 1,m
         j1 = (jj-1)*nblk+1
         j2 =  min0 (n,j1+nblk-1)
         do 12 j=j1,j2
            do 13 i=ia(j),ia(j+1) -1
               k = ja(i)
               if (k .lt. j1) then
                  ko = ko+1
                  ao(ko) = a(i)
                  jao(ko) = k
               else if (k .lt. j) then
c     kb = (jj-1)*ltr+((j-j1)*(j-j1-1))/2+k-j1+1
c     bdiag(kb) = a(i)
                  bdiag(kb+k-j1) = a(i)
               end if
 13         continue
            kb = kb + j-j1
            iao(j+1) = ko+1
 12      continue
 11   continue
      return
      end
