      subroutine hsourc (indic,nx,nelx,node,x,y,ijk,fs,f)

c*********************************************************************72
c
cc HSOURC assembles the load vector F from element contributions in FS.
c
c generates the load vector f in assembled form from the
c the element contributions fs.
c indic = indicates if f is to be assembled (1) or not (zero)
c note: f(*) not initilazed. because might use values from boundary
c conditions.
c
      implicit double precision (a-h,o-z)
        double precision x(*),y(*),fs(*),f(*),xe(3),ye(3),det,areao3
      integer ijk(node,*)

      jnod = 0
      do 130 nel = 1,nelx
c
c get coordinates of nodal points
c
      do 104 i=1, node
      j = ijk(i,nel)
      xe(i) = x(j)
      ye(i) = y(j)
 104      continue
c
c compute determinant
c
      det=xe(2)*(ye(3)-ye(1))+xe(3)*(ye(1)-ye(2))+xe(1)*(ye(2)-ye(3))
c area3 = area/3
      areao3 = det/6.0
c
c contributions to nodes in the element
c
      if (indic .eq. 0) then
         do 115 ka=1,node
         jnod = jnod+1
         f(jnod) = fs(nel)*areao3
 115      continue
      else
      do 120 ka=1, node
            ii = ijk(ka,nel)
          f(ii) = f(ii) + fs(nel)*areao3
 120      continue
      end if

 130      continue
      return
      end
