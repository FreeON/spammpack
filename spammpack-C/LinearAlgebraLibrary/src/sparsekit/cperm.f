      subroutine cperm (nrow,a,ja,ia,ao,jao,iao,perm,job)

c*********************************************************************72
c
cc CPERM permutes the columns of a matrix.
c
c  This routine permutes the columns of a matrix a, ja, ia.
c the result is written in the output matrix  ao, jao, iao.
c cperm computes B = A P, where  P is a permutation matrix
c that maps column j into column perm(j), i.e., on return
c      a(i,j) becomes a(i,perm(j)) in new matrix
c Y. Saad, May 2, 1990 / modified Jan. 28, 1991.
c
c on entry:
c
c nrow       = row dimension of the matrix
c
c a, ja, ia = input matrix in csr format.
c
c perm      = integer array of length ncol (number of columns of A
c         containing the permutation array  the columns:
c         a(i,j) in the original matrix becomes a(i,perm(j))
c         in the output matrix.
c
c job      = integer indicating the work to be done:
c             job = 1      permute a, ja, ia into ao, jao, iao
c                       (including the copying of real values ao and
c                       the array iao).
c             job .ne. 1 :  ignore real values ao and ignore iao.
c
c
c on return:
c
c ao, jao, iao = input matrix in a, ja, ia format (array ao not needed)
c
c Notes:
c
c 1. if job=1 then ao, iao are not used.
c 2. This routine is in place: ja, jao can be the same.
c 3. If the matrix is initially sorted (by increasing column number)
c    then ao,jao,iao  may not be on return.
c
c
      integer nrow,ja(*),ia(nrow+1),jao(*),iao(nrow+1),perm(*), job
      double precision a(*), ao(*)
      integer k, i, nnz

      nnz = ia(nrow+1)-1
      do 100 k=1,nnz
         jao(k) = perm(ja(k))
 100  continue
c
c     done with ja array. return if no need to touch values.
c
      if (job .ne. 1) return
c
c else get new pointers -- and copy values too.
c
      do 1 i=1, nrow+1
         iao(i) = ia(i)
 1    continue

      do 2 k=1, nnz
         ao(k) = a(k)
 2    continue

      return
      end
