      subroutine hes (ndg,m,hh,ih,dt,y,root,coef,coef0,w2)

c*********************************************************************72
c
cc HES computes exp ( H dt) * y  where H = Hessenberg matrix (hh)
c
c ndg      = number of poles as determined by getrat
c m     = dimension of hessenberg matrix
c hh      = hessenberg matrix (real)
c ih      = first dimenbsion of hh
c dt      = scaling factor used for hh (see (1))
c y      = double precision vector. on return exp(H dt ) y is computed
c         and overwritten on y.
c root  = poles of the rational approximation to exp as
c         computed by getrat
c coef,
c coef0 = coefficients of partial fraction expansion
c
c  exp(t) ~ coef0 +  sum     Real [   coef(i) / (t - root(i)  ]
c                  i=1,ndg
c
c valid for double precision real t.
c coef0 is double precision, coef(*) is a double complex array.
c
      parameter (mmax=70)
      implicit double precision (a-h,o-z)
      double precision  hh(ih,*), y(*)
      double complex hloc(mmax+1,mmax), coef(*),root(*),w2(*),t,zpiv
      double precision yloc(mmax)
c
c  loop associated with the poles.
c
      do 10 j=1,m
         yloc(j) = y(j)
         y(j) = y(j)*coef0
 10   continue

      do 8 ii = 1, ndg
c
c  copy Hessenberg matrix into temporary
c
         do 2 j=1, m
            do 1 i=1, j+1
               hloc(i,j) = CMPLX( dt*hh(i,j) )
 1          continue
            hloc(j,j) = hloc(j,j) - root(ii)
            w2(j)     = CMPLX(yloc(j))
 2       continue
c
c  forward solve
c
         do 4 i=2,m
            zpiv  = hloc(i,i-1) / hloc(i-1,i-1)
            do 3 j=i,m
               hloc(i,j) = hloc(i,j) - zpiv*hloc(i-1,j)
 3          continue
            w2(i)     = w2(i) - zpiv*w2(i-1)
 4       continue
c
c     backward solve
c
         do 6 i=m,1,-1
            t=w2(i)
            do 5 j=i+1,m
               t = t-hloc(i,j)*w2(j)
 5          continue
            w2(i) = t/hloc(i,i)
 6       continue
c
c     accumulate result in y.
c
         do 7 i=1,m
            y(i) = y(i) + coef(ii) * w2(i)
 7       continue
 8    continue
      return
      end
