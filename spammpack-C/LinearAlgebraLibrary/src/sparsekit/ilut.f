       subroutine ilut (n,a,ja,ia,lfil,tol,alu,jlu,ju,iwk,
     *                  wu,wl,jr,jwl,jwu,ierr)

c*********************************************************************72
c
cc ILUT is an ILUT preconditioner.
c
c      incomplete LU factorization with dual truncation mechanism
c      VERSION 2 : sorting  done for both L and U.
c
c
c  coded by Youcef Saad May, 5, 1990.
c  Dual drop-off strategy works as follows.
c
c     1) Theresholding in L and U as set by tol. Any element whose size
c        is less than some tolerance (relative to the norm of current
c        row in u) is dropped.
c
c     2) Keeping only the largest lenl0+lfil elements in L and the
c        largest lenu0+lfil elements in U, where lenl0=initial number
c        of nonzero elements in a given row of lower part of A
c        and lenlu0 is similarly defined.
c
c Flexibility: one can use tol=0 to get a strategy based on keeping the
c largest elements in each row of L and U. Taking tol .ne. 0 but lfil=n
c will give the usual threshold strategy (however, fill-in is then
c impredictible).
c
c*
c PARAMETERS
c
c on entry:
c
c n       = integer. The dimension of the matrix A.
c
c a,ja,ia = matrix stored in Compressed Sparse Row format.
c
c lfil    = integer. The fill-in parameter. Each row of L and
c           each row of U will have a maximum of lfil elements
c           in addition to the original number of nonzero elements.
c           Thus storage can be determined beforehand.
c           lfil must be .ge. 0.
c
c iwk     = integer. The minimum length of arrays alu and jlu
c
c On return:
c
c
c alu,jlu = matrix stored in Modified Sparse Row (MSR) format containing
c           the L and U factors together. The diagonal (stored in
c           alu(1:n) ) is inverted. Each i-th row of the alu,jlu matrix
c           contains the i-th row of L (excluding the diagonal entry=1)
c           followed by the i-th row of U.
c
c ju      = integer array of length n containing the pointers to
c           the beginning of each row of U in the matrix alu,jlu.
c
c ierr    = integer. Error message with the following meaning.
c           ierr  = 0    --> successful return.
c           ierr .gt. 0  --> zero pivot encountered at step number ierr.
c           ierr  = -1   --> Error. input matrix may be wrong.
c                            (The elimination process has generated a
c                            row in L or U whose length is .gt.  n.)
c           ierr  = -2   --> The matrix L overflows the array al.
c           ierr  = -3   --> The matrix U overflows the array alu.
c           ierr  = -4   --> Illegal value for lfil.
c           ierr  = -5   --> zero pivot encountered.
c
c work arrays:
c
c jr,jwu,jwl         = integer work arrays of length n.
c wu, wl          = double precision work arrays of length n+1, and n resp.
c
c Notes:
c
c A must have all nonzero diagonal elements.
c
       implicit double precision (a-h,o-z)
      INTEGER N
       double precision a(*), alu(*), wu(n), wl(n), tol
       integer ja(*),ia(n+1),jlu(*),ju(n),jr(n), jwu(n),
     *      jwl(n), lfil, iwk, ierr

        if (lfil .lt. 0) goto 998
c
c initialize ju0 (points to next element to be added to alu,jlu)
c and pointer.
c
      ju0 = n+2
      jlu(1) = ju0
c
c  integer double pointer array.
c
      do 1 j=1, n
            jr(j)  = 0
 1           continue
c
c  beginning of main loop.
c
      do 500 ii = 1, n
           j1 = ia(ii)
           j2 = ia(ii+1) - 1
           lenu = 0
           lenl = 0
           tnorm = 0.0
           do 501 k=j1,j2
              tnorm = tnorm+abs(a(k))
 501          continue
              tnorm = tnorm/dble(j2-j1+1)
c
c  unpack L-part and U-part of row of A in arrays wl, wu --
c
      do 170  j = j1, j2
           k = ja(j)
           t = a(j)
           if (abs(t) .lt. tol*tnorm) goto 170
           if (k .lt. ii) then
              lenl = lenl+1
              jwl(lenl) = k
              wl(lenl) = t
              jr(k) = lenl
           else
              lenu = lenu+1
              jwu(lenu) = k
              wu(lenu) = t
              jr(k) = lenu
           end if
 170      continue
c      tnorm = 0.0
c      do 171 k=j1,j2
c           tnorm = tnorm + abs(a(k))
c 171      continue
c
c        tnorm = tnorm/dble(j2-j1+1)
        lenl0 = lenl
        lenu0 = lenu
        jj = 0
        nl = 0
c
c  eliminate previous rows
c
 150    jj = jj+1
        if (jj .gt. lenl) goto 160
c
c in order to do the elimination in the correct order we need to
c exchange the current row number with the one that has
c smallest column number, among jj,jj+1,...,lenl.
c
        jrow = jwl(jj)
        k = jj
c
c determine smallest column index
c
        do 151 j=jj+1,lenl
           if (jwl(j) .lt. jrow) then
              jrow = jwl(j)
              k = j
           end if
 151    continue
c
c exchange in jwl
c
        j = jwl(jj)
        jwl(jj) = jrow
        jwl(k) = j
c
c exchange in jr
c
        jr(jrow) = jj
        jr(j) = k
c
c exchange in wl
c
        s = wl(k)
        wl(k) = wl(jj)
        wl(jj) = s
c
        if (jrow .ge. ii) goto 160
c  get the multiplier for row to be eliminated: jrow
        fact = wl(jj)*alu(jrow)
        jr(jrow) = 0
        if (abs(fact)*wu(n+2-jrow) .le. tol*tnorm) goto 150
c
c  combine current row and row jrow
c
        do 203 k = ju(jrow), jlu(jrow+1)-1
           s = fact*alu(k)
           j = jlu(k)
           jpos = jr(j)
c
c if fill-in element and small disregard:
c
           if (abs(s) .lt. tol*tnorm .and. jpos .eq. 0) goto 203
           if (j .ge. ii) then
c
c     dealing with upper part.
c
              if (jpos .eq. 0) then
c     this is a fill-in element
                 lenu = lenu+1
                 if (lenu .gt. n) goto 995
                 jwu(lenu) = j
                 jr(j) = lenu
                 wu(lenu) = - s
              else
c     no fill-in element --
                 wu(jpos) = wu(jpos) - s
              end if
           else
c
c     dealing with lower part.
c
              if (jpos .eq. 0) then
c     this is a fill-in element
                 lenl = lenl+1
                 if (lenl .gt. n) goto 995
                 jwl(lenl) = j
                 jr(j) = lenl
                 wl(lenl) = - s
              else
c     no fill-in element --
                 wl(jpos) = wl(jpos) - s
              end if
           end if
 203      continue
        nl = nl+1
        wl(nl) = fact
        jwl(nl)  = jrow
      goto 150
c
c  update l-matrix
c
 160    len = min0(nl,lenl0+lfil)
        call bsort2 (wl,jwl,nl,len)
c
        do 204 k=1, len
           if (ju0 .gt. iwk) goto 996
           alu(ju0) =  wl(k)
           jlu(ju0) =  jwl(k)
           ju0 = ju0+1
 204    continue
c
c  save pointer to beginning of row ii of U
c
        ju(ii) = ju0
c
c  reset double-pointer jr to zero (L-part - except first
c  jj-1 elements which have already been reset)
c
      do 306 k= jj, lenl
              jr(jwl(k)) = 0
 306      continue
c
c
c be sure that the diagonal element is first in w, jw
c
        idiag = 0
        idiag = jr(ii)
        if (idiag .eq. 0) goto 900

        if (idiag .ne. 1) then
           s = wu(1)
           wu(j) = wu(idiag)
           wu(idiag) = s

           j = jwu(1)
           jwu(1) = jwu(idiag)
           jwu(idiag) = j

        end if

        len = min0(lenu,lenu0+lfil)
      call bsort2 (wu(2), jwu(2), lenu-1,len)
c
c  update u-matrix
c
        t = 0.0
        do 302 k=2, len
           if (ju0 .gt. iwk) goto 997
           jlu(ju0) = jwu(k)
           alu(ju0)  = wu(k)
           t = t+ abs(wu(k) )
           ju0 = ju0+1
 302      continue
c
c     save norm in wu (backwards). Norm is in fact average abs value
c
        wu(n+2-ii) = t / dble(len+1)
c
c     store inverse of diagonal element of u
c
        if (wu(1) .eq. 0.0) goto 999

        alu(ii) = 1.0/ wu(1)
c
c     update pointer to beginning of next row of U.
c
      jlu(ii+1) = ju0
c
c     reset double-pointer jr to zero (U-part)
c
      do 308 k=1, lenu
           jr(jwu(k)) = 0
 308      continue
c
c     end main loop
c
 500  continue
        ierr = 0
        return
c
c     zero pivot :
c
 900    ierr = ii
        return
c
c     incomprehensible error. Matrix must be wrong.
c
 995    ierr = -1
        return
c
c     insufficient storage in L.
c
 996    ierr = -2
        return
c
c     insufficient storage in U.
c
 997    ierr = -3
        return
c
c     illegal lfil entered.
c
 998    ierr = -4
        return
c
c     zero pivot encountered
c
 999    ierr = -5
        return
        end
