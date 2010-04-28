      subroutine dscaldg (n,a,ja,ia,diag,job)

c*********************************************************************72
c
cc DSCALDG scales rows by a diagonal factor.
c
c scales rows by diag where diag is either given (job=0)
c or to be computed:
c  job = 1 ,scale row i by  by  +/- max |a(i,j) | and put inverse of
c       scaling factor in diag(i),where +/- is the sign of a(i,i).
c  job = 2 scale by 2-norm of each row..
c if diag(i) = 0,then diag(i) is replaced by one
c (no scaling)..
c
c           Y. Saad, Sep. 21 1989
c
      double precision a(*), diag(*),t
      integer ia(*),ja(*)

      goto (12,11,10) job+1
 10   do 110 j=1,n
         k1= ia(j)
         k2 = ia(j+1)-1
         t = 0.0
         do 111 k = k1,k2
 111        t = t+a(k)*a(k)
 110        diag(j) = sqrt(t)
            goto 12
 11   continue
      call retmx (n,a,ja,ia,diag)

 12   do 1 j=1,n
         if (diag(j) .ne. 0.0) then
            diag(j) = 1.0/diag(j)
         else
            diag(j) = 1.0
         end if
 1    continue
      do 2 i=1,n
         t = diag(i)
         do 21 k=ia(i),ia(i+1) -1
            a(k) = a(k)*t
 21      continue
 2    continue
      return
      end
