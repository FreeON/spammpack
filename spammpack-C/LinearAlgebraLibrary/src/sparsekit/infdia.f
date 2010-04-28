      subroutine infdia (n,ja,ia,ind,idiag)

c*********************************************************************72
c
cc INFDIA obtains information on the diagonals of A.
c
c this subroutine finds the lengths of each of the 2*n-1 diagonals of A
c it also outputs the number of nonzero diagonals found.
c
c on entry:
c
c n      = dimension of the matrix a.
c
c a,    not needed here.
c ja,
c ia    = matrix stored in csr format
c
c on return:
c
c
c idiag = integer. number of nonzero diagonals found.
c
c ind   = integer array of length at least 2*n-1. The k-th entry in
c         ind contains the number of nonzero elements in the diagonal
c         number k, the numbering beeing from the lowermost diagonal
c         (bottom-left). In other words ind(k) = length of diagonal
c         whose offset wrt the main diagonal is = - n + k.
c
c           Y. Saad, Sep. 21 1989
c
      integer ia(*), ind(*), ja(*)

      n2= n+n-1
      do 1 i=1,n2
         ind(i) = 0
 1    continue
      do 3 i=1, n
         do 2 k=ia(i),ia(i+1)-1
            j = ja(k)
            ind(n+j-i) = ind(n+j-i) +1
 2       continue
 3    continue
c     count the nonzero ones.
      idiag = 0
      do 41 k=1, n2
         if (ind(k) .ne. 0) idiag = idiag+1
 41   continue
      return
      end
