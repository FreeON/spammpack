      subroutine cnrms   (nrow, nrm, a, ja, ia, diag)

c*********************************************************************72
c
cc CNRMS gets the norms of each column of A.
c
c gets the norms of each column of A. (choice of three norms)
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
      double precision a(*), diag(nrow)
      integer ja(*), ia(nrow+1)

      do 10 k=1, nrow
         diag(k) = 0.0
 10   continue
      do 1 ii=1,nrow
         k1 = ia(ii)
         k2 = ia(ii+1)-1
         do 2 k=k1, k2
            j = ja(k)
c     update the norm of each column
            if (nrm .eq. 0) then
               diag(j) = max(diag(j),abs(a(k) ) )
            elseif (nrm .eq. 1) then
               diag(j) = diag(j) + abs(a(k) )
            else
               diag(j) = diag(j)+a(k)**2
            end if
 2       continue
 1    continue
      if (nrm .ne. 2) return
      do 3 k=1, nrow
         diag(k) = sqrt(diag(k))
 3    continue
      return
      end
