      subroutine project(n,m,u,v,w)

c*********************************************************************72
c
cc PROJECT computes the matrix-vector product w = U * v.
c
      implicit double precision (a-h,o-z)
      double precision u(n,*),v(*),w(*)

      do k=1,n
         w(k) = 0.0
      end do

      do j=1,m
         do k=1,n
            w(k) = w(k) + v(j) * u(k,j)
         end do
      end do

      return
      end
