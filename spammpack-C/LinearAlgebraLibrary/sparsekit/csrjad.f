      subroutine csrjad (nrow, a, ja, ia, idiag, iperm, ao, jao, iao)

c*********************************************************************72
c
cc CSRJAD converts Compressed Sparse Row to Jagged Diagonal storage.
c
c CSRJAG converts  matrix stored in the compressed sparse
c row format to the jagged diagonal format. The data structure
c for the JAD (Jagged Diagonal storage) is as follows. The rows of
c the matrix are (implicitly) permuted so that their lengths are in
c decreasing order. The real entries ao(*) and their column indices
c jao(*) are stored in succession. The number of such diagonals is idiag.
c the lengths of each of these diagonals is stored in iao(*).
c For more details see [E. Anderson and Y. Saad,
c ``Solving sparse triangular systems on parallel computers'' in
c Inter. J. of High Speed Computing, Vol 1, pp. 73-96 (1989).]
c or  [Y. Saad, ``Krylov Subspace Methods on Supercomputers''
c SIAM J. on  Stat. Scient. Comput., volume 10, pp. 1200-1232 (1989).]
c
c on entry:
c
c nrow         = row dimension of the matrix A.
c
c a,
c ia,
c ja      = input matrix in compressed sparse row format.
c
c on return:
c
c
c idiag = integer. The number of jagged diagonals in the matrix.
c
c iperm = integer array of length nrow containing the permutation
c         of the rows that leads to a decreasing order of the
c         number of nonzero elements.
c
c ao    = double precision array containing the values of the matrix A in
c         jagged diagonal storage. The j-diagonals are stored
c         in ao in sequence.
c
c jao   = integer array containing the column indices of the
c         entries in ao.
c
c iao   = integer array containing pointers to the beginning
c         of each j-diagonal in ao, jao. iao is also used as
c         a work array and it should be of length n at least.
c
      integer ja(*), jao(*), ia(nrow+1), iperm(nrow), iao(nrow)
      double precision a(*), ao(*)
c
c  define initial iperm and get lengths of each row
c  jao is used a work vector to store tehse lengths
c
      idiag = 0
      ilo = nrow
      do 10 j=1, nrow
         iperm(j) = j
         len = ia(j+1) - ia(j)
         ilo = min(ilo,len)
         idiag = max(idiag,len)
         jao(j) = len
 10   continue
c
c     call sorter to get permutation. use iao as work array.
c
      call dcsort (jao, nrow, iao, iperm, ilo, idiag)
c
c     define output data structure. first lengths of j-diagonals
c
      do 20 j=1, nrow
         iao(j) = 0
 20   continue
      do 40 k=1, nrow
         len = jao(iperm(k))
         do 30 i=1,len
            iao(i) = iao(i)+1
 30      continue
 40   continue
c
c     get the output matrix itself
c
      k1 = 1
      k0 = k1
      do 60 jj=1, idiag
         len = iao(jj)
         do 50 k=1,len
            i = ia(iperm(k))+jj-1
            ao(k1) = a(i)
            jao(k1) = ja(i)
            k1 = k1+1
 50      continue
         iao(jj) = k0
         k0 = k1
 60   continue
      iao(idiag+1) = k1
      return
      end
