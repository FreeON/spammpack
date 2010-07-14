      subroutine coocsr(nrow,nnz,a,ir,jc,ao,jao,iao)

c*********************************************************************72
c
cc COOCSR converts COO to CSR.
c
c converts a matrix that is stored in coordinate format
c  a, ir, jc into a row general sparse ao, jao, iao format.
c
c on entry:
c
c nrow      = dimension of the matrix
c nnz      = number of nonzero elements in matrix
c a,
c ir,
c jc    = matrix in coordinate format. a(k), ir(k), jc(k) store the nnz
c         nonzero elements of the matrix with a(k) = actual real value of
c         the elements, ir(k) = its row number and jc(k) = its column
c        number. The order of the elements is arbitrary.
c
c on return:
c
c ir       is destroyed
c
c ao, jao, iao = matrix in general sparse matrix format with ao
c       continung the real values, jao containing the column indices,
c      and iao being the pointer to the beginning of the row,
c      in arrays ao, jao.
c
c Notes:
c This routine is NOT in place.  See coicsr
c
      double precision a(*),ao(*),x
      integer ir(*),jc(*),jao(*),iao(*)

      do 1 k=1,nrow+1
         iao(k) = 0
 1    continue
c determine row-lengths.
      do 2 k=1, nnz
         iao(ir(k)) = iao(ir(k))+1
 2    continue
c starting position of each row.
      k = 1
      do 3 j=1,nrow+1
         k0 = iao(j)
         iao(j) = k
         k = k+k0
 3    continue
c go through the structure  once more. Fill in output matrix.
      do 4 k=1, nnz
         i = ir(k)
         j = jc(k)
         x = a(k)
         iad = iao(i)
         ao(iad) =  x
         jao(iad) = j
         iao(i) = iad+1
 4    continue
c shift back iao
      do 5 j=nrow,1,-1
         iao(j+1) = iao(j)
 5    continue
      iao(1) = 1
      return
      end
