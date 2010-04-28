      subroutine retmx(n,a,ja,ia,dd)

c*********************************************************************72
c
cc RETMX returns in dd(*) the max absolute value of elements in row *.
c
c RETMX returns in dd(*) the max absolute value of elements in row *.
c used for scaling purposes. superseded by rnrms  .
c
c on entry:
c n      = dimension of A
c a,ja,ia
c      = matrix stored in compressed sparse row format
c dd      = double precision array of length n. On output,entry dd(i) contains
c        the element of row i that has the largest absolute value.
c        Moreover the sign of dd is modified such that it is the
c        same as that of the diagonal element in row i.
c
c           Y. Saad, Sep. 21 1989                                      c
c
      double precision a(*),dd(*)
      integer n,ia(*),ja(*)
      integer k2, i, k1, k
      double precision t, t1, t2
c
c initialize
c
      k2 = 1
      do 11 i=1,n
         k1 = k2
         k2 = ia(i+1) - 1
         t = 0.0
         do 101  k=k1,k2
            t1 = abs(a(k))
            if (t1 .gt. t) t = t1
            if (ja(k) .eq. i)then
              if(a(k).lt.0.0)then
                t2=-1.0
              elseif(a(k).eq.0.0)then
                t2=0.0
              else
                t2=1.0
                end if
              end if
 101     continue
         dd(i) =  t2*t
c     we do not invert diag here..
 11   continue
      return
      end
