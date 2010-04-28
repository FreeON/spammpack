      subroutine diapos  (n,ja,ia,idiag)

c*********************************************************************72
c
cc DIAPOS returns the positions of the diagonal elements of a sparse matrix.
c
c on entry:
c
c
c n      = integer. row dimension of the matrix a.
c a,ja,
c    ia = matrix stored compressed sparse row format. a array skipped.
c
c on return:
c
c idiag  = integer array of length n. The i-th entry of idiag
c          points to the diagonal element a(i,i) in the arrays
c          a, ja. (i.e., a(idiag(i)) = element A(i,i) of matrix A)
c          if no diagonal element is found the entry is set to 0.
c
c           Y. Saad, March, 1990
c
      integer ia(n+1), ja(*), idiag(n)

      do 1 i=1, n
         idiag(i) = 0
 1    continue
c
c     sweep through data structure.
c
      do  6 i=1,n
         do 51 k= ia(i),ia(i+1) -1
            if (ja(k) .eq. i) idiag(i) = k
 51      continue
 6    continue
      return
      end
