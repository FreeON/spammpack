      subroutine cooell(n,nnz,a,ja,ia,ac,jac,nac,ner,ncmax,ierr)

c*********************************************************************72
c
cc COOELL converts coordinate format to Ellpack/Itpack format.
c
c   DATE WRITTEN: June 4, 1989.
c
c   PURPOSE
c
c  COOELL takes a sparse matrix in coordinate format and
c  converts it into the Ellpack-Itpack storage.
c
c  Example:
c
c       (   11   0   13    0     0     0  )
c       |   21  22    0   24     0     0  |
c       |    0  32   33    0    35     0  |
c   A = |    0   0   43   44     0    46  |
c       |   51   0    0   54    55     0  |
c       (   61  62    0    0    65    66  )
c
c   Coordinate storage scheme:
c
c    A  = (11,22,33,44,55,66,13,21,24,32,35,43,46,51,54,61,62,65)
c    IA = (1, 2, 3, 4, 5, 6, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 6 )
c    JA = ( 1, 2, 3, 4, 5, 6, 3, 1, 4, 2, 5, 3, 6, 1, 4, 1, 2, 5)
c
c   Ellpack-Itpack storage scheme:
c
c       (   11  13    0    0   )          (   1   3   *    *  )
c       |   22  21   24    0   |          |   2   1   4    *  |
c  AC = |   33  32   35    0   |    JAC = |   3   2   5    *  |
c       |   44  43   46    0   |          |   4   3   6    *  |
c       |   55  51   54    0   |          |   5   1   4    *  |
c       (   66  61   62   65   )          (   6   1   2    5  )
c
c   Note: * means that you can store values from 1 to 6 (1 to n, where
c         n is the order of the matrix) in that position in the array.
c
c   Contributed by:
c
c   Ernest E. Rothman
c   Cornell Thoery Center/Cornell National Supercomputer Facility
c   e-mail address: BITNET:   EER@CORNELLF.BITNET
c                   INTERNET: eer@cornellf.tn.cornell.edu
c
c   checked and modified  04/13/90 Y.Saad.
c
c   REFERENCES
c
c   Kincaid, D. R.; Oppe, T. C.; Respess, J. R.; Young, D. M. 1984.
c   ITPACKV 2C User's Guide, CNA-191. Center for Numerical Analysis,
c   University of Texas at Austin.
c
c   "Engineering and Scientific Subroutine Library; Guide and
c   Reference; Release 3 (SC23-0184-3). Pp. 79-86.
c
c   INPUT PARAMETERS
c
c  N       - Integer. The size of the square matrix.
c
c  NNZ     - Integer. Must be greater than or equal to the number of
c            nonzero elements in the sparse matrix. Dimension of A, IA
c            and JA.
c
c  NCA     - Integer. First dimension of output arrays ca and jac.
c
c  A(NNZ)  - Real array.
c            Stored entries of the sparse matrix A.
c            NNZ is the number of nonzeros.
c
c  IA(NNZ) - Integer array.
c            Pointers to specify rows for the stored nonzero entries
c            in A.
c
c  JA(NNZ) - Integer array.
c            Pointers to specify columns for the stored nonzero
c            entries in A.
c
c  NER     - Integer. Must be set greater than or equal to the maximum
c            number of nonzeros in any row of the sparse matrix.
c
c  OUTPUT PARAMETERS
c
c  AC(NAC,*)  - Real array.
c               Stored entries of the sparse matrix A in compressed
c               storage mode.
c
c  JAC(NAC,*) - Integer array.
c               Contains the column numbers of the sparse matrix
c               elements stored in the corresponding positions in
c               array AC.
c
c  NCMAX   -  Integer. Equals the maximum number of nonzeros in any
c             row of the sparse matrix.
c
c  IERR    - Error parameter is returned as zero on successful
c             execution of the subroutin<e.
c             Error diagnostics are given by means of positive values
c             of this parameter as follows:
c
c             IERR = -1   -  NER is too small and should be set equal
c                            to NCMAX. The array AC may not be large
c                            enough to accomodate all the non-zeros of
c                            of the sparse matrix.
c             IERR =  1   -  The array AC has a zero column. (Warning)
c             IERR =  2   -  The array AC has a zero row.    (Warning)
c
c
      double precision a(nnz), ac(nac,ner)
      integer ja(nnz), ia(nnz), jac(nac,ner), ierr, ncmax, icount
c
c   Initial error parameter to zero:
c
      ierr = 0
c
c   Initial output arrays to zero:
c
      do 4 in = 1,ner
         do 4 innz =1,n
            jac(innz,in) = n
            ac(innz,in) = 0.0
 4    continue
c
c   Assign nonzero elements of the sparse matrix (stored in the one
c   dimensional array A to the two dimensional array AC.
c   Also, assign the correct values with information about their
c   column indices to the two dimensional array KA. And at the same
c   time count the number of nonzeros in each row so that the
c   parameter NCMAX equals the maximum number of nonzeros in any row
c   of the sparse matrix.
c
      ncmax = 1
      do 10 is = 1,n
         k = 0
         do 30 ii = 1,nnz
            if(ia(ii).eq.is)then
               k = k + 1
               if (k .le. ner) then
                  ac(is,k) = a(ii)
                  jac(is,k) = ja(ii)
               end if
            end if
 30      continue
         if (k.ge.ncmax) ncmax = k
 10   continue
c
c     Perform some simple error checks:
c
c  check maximum number of nonzeros in each row:
c
      if (ncmax.eq.ner) ierr = 0
      if (ncmax.gt.ner) then
         ierr = -1
         return
      end if
c
c  check if there are any zero columns in AC:
c
      do 45 in = 1,ncmax
         icount = 0
         do 44 inn =1,n
            if (ac(inn,in).ne.0.0) icount = 1
 44      continue
         if (icount.eq.0) then
            ierr = 1
            return
         end if
 45   continue
c
c  check if there are any zero rows in AC:
c
      do 55 inn = 1,n
         icount = 0
         do 54 in =1,ncmax
            if (ac(inn,in).ne.0.0) icount = 1
 54      continue
         if (icount.eq.0) then
            ierr = 2
            return
         end if
 55   continue
      return
      end
