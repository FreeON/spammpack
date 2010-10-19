      subroutine rperm (nrow,a,ja,ia,ao,jao,iao,perm,job)

c*********************************************************************72
c
cc RPERM permutes the rows of a matrix in CSR format.
c
c rperm  computes B = P A  where P is a permutation matrix.
c the permutation P is defined through the array perm: for each j,
c perm(j) represents the destination row number of row number j.
c Youcef Saad -- recoded Jan 28, 1991.
c
c on entry:
c
c n       = dimension of the matrix
c a, ja, ia = input matrix in csr format
c perm       = integer array of length nrow containing the permutation arrays
c        for the rows: perm(i) is the destination of row i in the
c         permuted matrix.
c         ---> a(i,j) in the original matrix becomes a(perm(i),j)
c         in the output  matrix.
c
c job      = integer indicating the work to be done:
c             job = 1      permute a, ja, ia into ao, jao, iao
c                       (including the copying of real values ao and
c                       the array iao).
c             job .ne. 1 :  ignore real values.
c                     (in which case arrays a and ao are not needed nor
c                      used).
c
c
c on return:
c
c ao, jao, iao = input matrix in a, ja, ia format
c note :
c        if (job.ne.1)  then the arrays a and ao are not used.
c
c           Y. Saad, May  2, 1990                                      c
c
      integer nrow,ja(*),ia(nrow+1),jao(*),iao(nrow+1),perm(nrow),job
      double precision a(*),ao(*)
c
      logical values
      values = (job .eq. 1)
c
c     determine pointers for output matix.
c
      do j=1,nrow
         i = perm(j)
         iao(i+1) = ia(j+1) - ia(j)
      end do
c
c get pointers from lengths
c
      iao(1) = 1
      do 51 j=1,nrow
         iao(j+1)=iao(j+1)+iao(j)
 51   continue
c
c copying
c
      do 100 ii=1,nrow
c
c old row = ii  -- new row = iperm(ii) -- ko = new pointer
c
         ko = iao(perm(ii))
         do 60 k=ia(ii), ia(ii+1)-1
            jao(ko) = ja(k)
            if (values) ao(ko) = a(k)
            ko = ko+1
 60      continue
 100  continue
c
      return
      end
