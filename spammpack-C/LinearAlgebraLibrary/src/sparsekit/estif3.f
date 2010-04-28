      subroutine estif3(nel,ske,fe,det,xe,ye,xyke,ierr)

c*********************************************************************72
c
cc ESTIF3 constructs an element stiffness matrix using 3 node triangles.
c
c constructs the element stiffness matrix using 3-node triangular elements
c arguments:
c nel      = element number
c ske      = element stiffness matrix
c fe      = element load vector
c det      = 2*area of the triangle
c xy, ye= coordinates of the three nodal points in an element.
c xyke  = material constants (kxx, kxy, kyx, kyy)
c
      implicit double precision (a-h,o-z)
      dimension ske(3,3), fe(3), xe(3), ye(3), dn(3,2),xyke(2,2)
c
c initialize
c
      area = 0.5*det

      do 200 i=1,3
        fe(i) = 0.0
      do 200 j=1,3
      ske(i,j) = 0.0
 200      continue
c
c get first gradient of shape function
c
      call gradi3(nel,xe,ye,dn,det,ierr)
      if (ierr .ne. 0) return

      do 100 i=1,3
      do 100 j=1,3
      t = 0.0
      do 102 k=1,2
      do 102 l=1,2
 102      t = t+xyke(k,l)*dn(i,k)*dn(j,l)
 100      ske(i,j) = t*area

      return
      end
