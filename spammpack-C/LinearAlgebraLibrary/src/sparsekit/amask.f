c  SPARSKIT.F  03 December 1991
c
      subroutine amask (nrow,ncol,a,ja,ia,jmask,imask,
     *                  c,jc,ic,iw,nzmax,ierr)

c*********************************************************************72
c
cc AMASK extracts a sparse matrix from a masked input matrix.
c
c  Discussion:
c
c    AMASK extracts a sparse matrix from the input matrix
c    by looking at positions defined by the mask jmask, imask
c
c    the  algorithm is in place: c, jc, ic can be the same as
c    a, ja, ia in which cas the code will overwrite the matrix c
c    on a, ja, ia
c
c On entry:
c
c nrow  = integer. row dimension of input matrix
c ncol      = integer. Column dimension of input matrix.
c
c a,
c ja,
c ia      = matrix in Compressed Sparse Row format
c
c jmask,
c imask = matrix defining mask (pattern only) stored in compressed
c         sparse row format.
c
c nzmax = length of arrays c and jc. see ierr.
c
c On return:
c
c
c a, ja, ia and jmask, imask are unchanged.
c
c c
c jc,
c ic      = the output matrix in Compressed Sparse Row format.
c
c ierr  = integer. serving as error message.c
c         ierr = 1  means normal return
c         ierr .gt. 1 means that amask stopped when processing
c         row number ierr, because there was not enough space in
c         c, jc according to the value of nzmax.
c
c work arrays:
c
c iw      = logical work array of length ncol.
c
      double precision a(*),c(*)
      integer ia(nrow+1),ja(*),jc(*),ic(nrow+1),jmask(*),imask(nrow+1)
      logical iw(ncol)

      ierr = 0
      len = 0

      do 1 j=1, ncol
         iw(j) = .false.
 1    continue
c
c     unpack the mask for row ii in iw
c
      do 100 ii=1, nrow
c
c     save pointer in order to be able to do things in place
c
         do 2 k=imask(ii), imask(ii+1)-1
            iw(jmask(k)) = .true.
 2       continue
c
c     add umasked elemnts of row ii
c
         k1 = ia(ii)
         k2 = ia(ii+1)-1
         ic(ii) = len+1
         do 200 k=k1,k2
            j = ja(k)
            if (iw(j)) then
               len = len+1
               if (len .gt. nzmax) then
                  ierr = ii
                  return
               end if
               jc(len) = j
               c(len) = a(k)
            end if
 200     continue

         do k=imask(ii), imask(ii+1)-1
            iw(jmask(k)) = .false.
         end do

 100  continue
      ic(nrow+1)=len+1
      return
      end
