      function getelm (i,j,a,ja,ia,iadd,sorted)

c*********************************************************************72
c
cc GETELM returns the element A(I,J) of a CSR matrix A.
c
c     Purpose:
c
c     This function returns the element a(i,j) of a matrix A,
c     for any pair (i,j).  The matrix is assumed to be stored
c     in Compressed Sparse Row (CSR) format. GETELM Performs a
c     binary search in the case where it is known that the elements
c     are sorted so that the column indices are in increasing order.
c     Also returns (in iadd) the address of the element a(i,j) in
c     arrays A and JA when the search is successsful (zero if not).
c
c     First contributed by Noel Nachtigal (MIT).
c     Recoded Jan. 20, 1991, by Y. Saad [In particular
c     added handling of the non-sorted case + the IADD output]
c
c     Parameters:
c
c On entry:
c
c     I      = the row index of the element sought (input).
c     J      = the column index of the element sought (input).
c     A      = the matrix A in compressed sparse row format (input).
c     JA     = the array of column indices (input).
c     IA     = the array of pointers to the rows' data (input).
c     SORTED = logical indicating whether the matrix is knonw to
c              have its column indices sorted in increasing order
c              (sorted=.true.) or not (sorted=.false.).
c              (input).
c On return:
c
c     GETELM = value of a(i,j).
c     IADD   = address of element a(i,j) in arrays a, ja if found,
c              zero if not found. (output)
c
c     Note: the inputs I and J are not checked for validity.
c
c     Noel M. Nachtigal October 28, 1990 -- Youcef Saad Jan 20, 1991.
c
      INTEGER I, IA(*), IADD, J, JA(*)
      double precision A(*)
      double precision getelm
      LOGICAL SORTED
      INTEGER IBEG, IEND, IMID, K
c
c     Initialization
c
      IADD = 0
      GETELM = 0.0
      IBEG = IA(I)
      IEND = IA(I+1)-1
c
c  case where matrix is not necessarily sorted
c
      IF (.NOT. SORTED) THEN
c
c scan the row - exit as soon as a(i,j) is found
c
         DO 5  K=IBEG, IEND
            IF (JA(K) .EQ.  J) THEN
               IADD = K
               GOTO 20
            end if
 5       CONTINUE
c
c     end unsorted case. begin sorted case
c
      ELSE
c
c     begin binary search.   Compute the middle index.
c
 10      IMID = ( IBEG + IEND ) / 2
c
c     test if  found
c
         IF (JA(IMID).EQ.J) THEN
            IADD = IMID
            GOTO 20
         end if
         IF (IBEG .GE. IEND) GOTO 20
c
c     else     Update the interval bounds.
c
         IF (JA(IMID).GT.J) THEN
            IEND = IMID -1
         ELSE
            IBEG = IMID +1
         end if
         GOTO 10
c
c     end both cases
c
      end if

 20   IF (IADD .NE. 0) GETELM = A(IADD)

      RETURN
      END
