      subroutine chkelmt (nx, x, y, nelx, ijk, node)

c*********************************************************************72
c
cc CHKELMT checks the labeling within each element and reorders if necessary.
c
c  CHKELMT checks the labeling within each element and reorders
c the nodes if they are not correctly ordered.
c
      implicit double precision (a-h,o-z)
      dimension ijk(node,1),x(1),y(1)

      do 1 nel =1, nelx
       det = x(ijk(2,nel))*(y(ijk(3,nel))-y(ijk(1,nel)))+
     *        x(ijk(3,nel))*(y(ijk(1,nel))-y(ijk(2,nel)))+
     *        x(ijk(1,nel))*(y(ijk(2,nel))-y(ijk(3,nel)))
c
c if determinant negative exchange last two nodes of elements.
c
      if (det .lt. 0.0) then
          j = ijk(2,nel)
          ijk(2,nel) = ijk(3,nel)
          ijk(3,nel) = j
      end if
 1      continue

      return
      end
