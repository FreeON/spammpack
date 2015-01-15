module test_schulz
contains

  real*8 function f3(x,a)
    real*8 :: x,a
    f3=(1d0/8d0)*(15d0-a*10d0*x+a*3*x**2)
  end function f3

end module test_schulz

program test_random_number

  use test_schulz

  integer, parameter :: MAX_RAND = 1000
  integer, parameter :: seed = 86456

  real*8 :: s(1:MAX_RAND), x(1:MAX_RAND), z(1:MAX_RAND), t(1:MAX_RAND)
  real*8 :: low1, low2, high1, high2, a

  call srand(seed)

  low1=1d10
  high1=-1d10
  do i=1,MAX_RAND
     s(i)=RAND()
     low1=min(low1,s(i))
     high1=max(high1,s(i))
  enddo
  low2=-5
  high2=3

  do i=1,MAX_RAND
     s(i)=low2 + (s(i) - low1) * (high2 - low2) / (high1 - low1)
     s(i)=10d0**s(i)
     ! WRITE(*,*)s(i)
  enddo

  write(*,*)' xo=Lmin/Lmax = '

  do i=1,MAX_RAND
     s(i)=s(i)/10d0**high2
     write(*,*)s(i)
  enddo

  z=1
  a=1
  XMin=1
  XMax=1
  do j=1,20
     do i=1,MAX_RAND
        x(i)=z(i)*s(i)*z(i)
        a=XMin/XMax
        t(i)=f3(x(i),1d0)
        z(i)=z(i)*t(i)
     enddo
     XMin= 1d10
     XMax=-1d10
     do i=1,MAX_RAND
        XMin=min(XMin,x(i))
        XMax=max(XMax,x(i))
     enddo

     write(*,*)j,XMin,sum(abs(x(:)-1d0)),sum(abs(z(:)-1d0/sqrt(s(:))))
  enddo

end program test_random_number
