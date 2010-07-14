      subroutine submat (n,job,i1,i2,j1,j2,a,ja,ia,nr,nc,ao,jao,iao)

c*********************************************************************72
c
cc SUBMAT extracts the submatrix A(i1:i2,j1:j2).
c
c extracts the submatrix A(i1:i2,j1:j2) and puts the result in
c matrix ao,iao,jao
c  In place: ao,jao,iao may be the same as a,ja,ia.
c
c on input
c
c n      = row dimension of the matrix
c i1,i2 = two integers with i2 .ge. i1 indicating the range of rows to be
c          extracted.
c j1,j2 = two integers with j2 .ge. j1 indicating the range of columns
c         to be extracted.
c         * There is no checking whether the input values for i1, i2, j1,
c           j2 are between 1 and n.
c a,
c ja,
c ia    = matrix in compressed sparse row format.
c
c job      = job indicator: if job .ne. 1 then the real values in a are NOT
c         extracted, only the column indices (i.e. data structure) are.
c         otherwise values as well as column indices are extracted...
c
c on output
c
c nr      = number of rows of submatrix
c nc      = number of columns of submatrix
c        * if either of nr or nc is nonpositive the code will quit.
c
c ao,
c jao,iao = extracted matrix in general sparse format with jao containing
c      the column indices,and iao being the pointer to the beginning
c      of the row,in arrays a,ja.
c
c           Y. Saad, Sep. 21 1989                                      c
c
      integer n,job,i1,i2,j1,j2,nr,nc,ia(*),ja(*),jao(*),iao(*)
      double precision a(*),ao(*)

      nr = i2-i1+1
      nc = j2-j1+1
c
      if ( nr .le. 0 .or. nc .le. 0) return
c
      klen = 0
c
c     simple procedure that proceeds row-wise...
c
      do 100 i = 1,nr
         ii = i1+i-1
         k1 = ia(ii)
         k2 = ia(ii+1)-1
         iao(i) = klen+1

         do k=k1,k2
            j = ja(k)
            if (j .ge. j1 .and. j .le. j2) then
               klen = klen+1
               if (job .eq. 1) ao(klen) = a(k)
               jao(klen) =j
            end if
         end do

 100  continue
      iao(nr+1) = klen+1
      return
      end
