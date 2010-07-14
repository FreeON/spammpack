      subroutine milu0 (n, a, ja, ia, alu, jlu, ju, iw, ierr)

c*********************************************************************72
c
cc MILU0 is a simple milu(0) preconditioner.
c
c Note that this has been coded in such a way that it can be used
c with pgmres. Normally, since the data structure of a, ja, ia is
c the same as that of a, ja, ia, savings can be made. In fact with
c some definitions (not correct for general sparse matrices) all we
c need in addition to a, ja, ia is an additional diagonal.
c Ilu0 is not recommended for serious problems. It is only provided
c here for comparison purposes.
c
c on entry:
c
c n       = dimension of matrix
c a, ja,
c ia      = original matrix in compressed sparse row storage.
c
c on return:
c
c alu,jlu = matrix stored in Modified Sparse Row (MSR) format containing
c           the L and U factors together. The diagonal (stored in
c           alu(1:n) ) is inverted. Each i-th row of the alu,jlu matrix
c           contains the i-th row of L (excluding the diagonal entry=1)
c           followed by the i-th row of U.
c
c ju        = pointer to the diagonal elements in alu, jlu.
c
c ierr        = integer indicating error code on return
c           ierr = 0 --> normal return
c           ierr = k --> code encountered a zero pivot at step k.
c work arrays:
c
c iw          = integer work array of length n.
c
c Note (IMPORTANT):
c
c it is assumed that the the elements in the input matrix are ordered
c    in such a way that in each row the lower part comes first and
c    then the upper part. To get the correct ILU factorization, it is
c    also necessary to have the elements of L ordered by increasing
c    column number. It may therefore be necessary to sort the
c    elements of a, ja, ia prior to calling milu0. This can be
c    achieved by transposing the matrix twice using csrcsc.
c
      implicit double precision (a-h,o-z)
      double precision a(*), alu(*)
      integer ja(*), ia(*), ju(*), jlu(*), iw(*)

          ju0 = n+2
          jlu(1) = ju0
c initialize work vector to zero's
      do 31 i=1, n
 31           iw(i) = 0
c
c  MAIN LOOP
c
      do 500 ii = 1, n
           js = ju0
c
c generating row number ii or L and U.
c
           do 100 j=ia(ii),ia(ii+1)-1
c
c copy row ii of a, ja, ia into row ii of alu, jlu (L/U) matrix.
c
              jcol = ja(j)
              if (jcol .eq. ii) then
                 alu(ii) = a(j)
                 iw(jcol) = ii
                 ju(ii)  = ju0
              else
                 alu(ju0) = a(j)
                 jlu(ju0) = ja(j)
                 iw(jcol) = ju0
                 ju0 = ju0+1
              end if
 100       continue
           jlu(ii+1) = ju0
           jf = ju0-1
           jm = ju(ii)-1
c s accumulates fill-in values
           s = 0.0
           do 150 j=js, jm
              jrow = jlu(j)
              tl = alu(j)*alu(jrow)
              alu(j) = tl
c  perform linear combination
              do 140 jj = ju(jrow), jlu(jrow+1)-1
                 jw = iw(jlu(jj))
                 if (jw .ne. 0) then
                       alu(jw) = alu(jw) - tl*alu(jj)
                    else
                       s = s + tl*alu(jj)
                    end if
 140          continue
 150       continue
c  invert and store diagonal element.
           alu(ii) = alu(ii)-s
           if (alu(ii) .eq. 0.0) goto 600
           alu(ii) = 1.0/alu(ii)
c  reset pointer iw to zero
           iw(ii) = 0
           do 201 i = js, jf
 201          iw(jlu(i)) = 0
 500       continue
           ierr = 0
           return
c     zero pivot :
 600       ierr = ii
           return
           end
