      subroutine mgsr (n, i0, i1, ss, r)

c*********************************************************************72
c
cc MGSR is a modified Gram - Schmidt with partial reorthogonalization.
c
c modified gram - schmidt  with  partial  reortho. the vector ss(*,i1) is
c orthogonalized against the first i vectors  of ss  (which  are  already
c orthogonal).  the coefficients of the orthogonalization are returned in
c the array r
c
      implicit double precision (a-h,o-z)
      double precision ss(n,1), r(1), hinorm, tet, ddot, t, sqrt
      do 53 j=1, i1
         r(j) = 0.0
 53   continue
      i = i1-1
      it = 0
 54   hinorm = 0.0
      it = it +1
      if (i .eq. 0) goto 56
c
      do 55 j=i0, i
         t = ddot(n, ss(1,j),1,ss(1,i1),1)
         hinorm = hinorm + t**2
         r(j) = r(j) + t
         call daxpy(n,-t,ss(1,j),1,ss(1,i1),1)
 55   continue
      t = ddot(n, ss(1,i1), 1, ss(1,i1), 1)
 56   continue
c
c     test for reorthogonalization see daniel et. al.
c     two reorthogonalization allowed.
c
      if (t*10.0 .le. hinorm .and. it .lt. 2) goto 54
      t =sqrt(t)
      r(i1)= t
      if (t .eq. 0.0) return
      t = 1.0/t
      do 57  k=1,n
         ss(k,i1) = ss(k,i1)*t
 57   continue
      return
      end
