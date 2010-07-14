      subroutine diric (nx,nint,a,ja,ia, f)

c*********************************************************************72
c
cc DIRIC accounts for Dirichlet boundary conditions.
c
c
      implicit double precision  (a-h,o-z)
      dimension a(*),ia(*),ja(*),f(*)
c call extract from UNARY
      call submat (nx,1,1,nint,1,nint,a,ja,ia,nr,nc,a,ja,ia)
      return
      END
