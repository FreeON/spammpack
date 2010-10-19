      subroutine rnrms(nrow, nrm, a, ja, ia, diag)

c*********************************************************************72
c
cc RNRMS gets the norms of each row of A. (choice of three norms)
c
c on entry:
c
c nrow      = integer. The row dimension of A
c
c nrm   = integer. norm indicator. nrm = 1, means 1-norm, nrm =2
c                  means the 2-nrm, nrm = 0 means max norm
c
c a,
c ja,
c ia   = Matrix A in compressed sparse row format.
c
c on return:
c
c
c diag = double precision vector of length nrow containing the norms
c
c
      double precision a(*), diag(nrow), scal
      integer ja(*), ia(nrow+1)
c
      do 1 ii=1,nrow
c
c     compute the norm if each element.
c
         scal = 0.0
         k1 = ia(ii)
         k2 = ia(ii+1)-1
         if (nrm .eq. 0) then
            do 2 k=k1, k2
               scal = max(scal,abs(a(k) ) )
 2          continue
         elseif (nrm .eq. 1) then
            do 3 k=k1, k2
               scal = scal + abs(a(k) )
 3          continue
         else
            do 4 k=k1, k2
               scal = scal+a(k)**2
 4          continue
         end if
         if (nrm .eq. 2) scal = sqrt(scal)
         diag(ii) = scal
 1    continue
      return
      end
