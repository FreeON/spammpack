      subroutine exphes (n,m,dt,eps,u,w,job,z,wkc,beta,errst,hh,ih,
     *                   x, y, indic,ierr)

c*********************************************************************72
c
cc EXPHES computes the Arnoldi basis.
c
c this subroutine computes the Arnoldi basis and the corresponding
c coeffcient vector in the approximation
c
c              w  ::= beta  Vm  ym
c               where ym = exp(- Hm *dt) * e1
c
c to the vector exp(-A dt) w where A is an arbitary matrix and
c w is a given input vector. In case job = 0 the arnoldi basis
c is recomputed. Otherwise the
c code assumes assumes that  u(*) contains an already computed
c arnoldi basis and computes only the y-vector (which is stored in v(*))
c
c en entry:
c
c n      = dimension of matrix
c
c m      = dimension of Krylov subspace (= degree of polynomial
c         approximation to the exponential used. )
c
c dt      = scalar by which to multiply matrix. Can be viewed
c         as a time step. dt must be positive [to be fixed].
c
c eps   = scalar indicating the relative error tolerated for the result.
c         the code will try to compute an answer such that
c         norm2(exactanswer-approximation) / norm2(w) .le. eps
c
c u      = work array of size n*(m+1) to contain the Arnoldi basis
c
c w      = double precision array of length n = input vector to  which exp(-A) is
c         to be applied.
c
c y     = duble precision work array of  size (m+1)
c wkc   = double complex work array of size (m+1)
c
c job      = integer. job indicator. If job .lt.  0 then the Arnoldi
c         basis is recomputed. If job .gt. 0 then it is assumed
c         that the user wants to use a previously computed Krylov
c         subspace but a different dt. Thus the Arnoldi basis and
c         the Hessenberg matrix Hm are not recomputed.
c        In that case the user should not modify the values of beta
c         and the matrices hh and u(n,*) when recalling phipro.
c         job = -1 : recompute basis and get an initial estimate for
c                    time step dt to be used.
c         job = 0  : recompute basis and do not alter dt.
c         job = 1  : do not recompute arnoldi basis.
c
c hh    = work array of size size at least (m+1)*m
c
c ih      = first dimension of hh as declared in the calling program.
c         ih must be .ge. m.
c
c  entries specific to the matrix
c
c diagonal storage is used :
c         a(n,ndiag) is a rectangular array with a(*,k) containing the
c         the diagonal offset by ioff(k) (negative or positive or zero)
c         i.e.,
c        a(i,jdiag) contains the element A(i,i+ioff(jdiag)) in the
c         usual dense storage scheme.
c
c a      = matrix in diagonal storage form
c ioff      = offsets  of diagonals
c ndiag = number of diagonals.
c
c on return:
c
c w2      = resulting vector w2 = exp(-A *dt) * w
c beta  = real equal to the 2-norm of w. Needed if exppro will
c         be recalled with the same Krylov subspace and a different
c         dt.
c errst = rough estimates of the 2-norm of the error.
c hh      = work array of dimension at least (m+1) x m
c
      parameter (ndmx=20)
      implicit double precision (a-h,o-z)
      double precision hh(ih,ih), u(n,*), w(*),z(m+1), errst,x(*), y(*)
      double precision alp0
      double complex alp(ndmx+1), rd(ndmx+1),wkc(ih)
      save
c  use degree 14 chebyshev all the time
      if (indic .ge. 3) goto 60
c
c  input fraction expansion of rational function
c
      ldg= 7
      alp0 =  0.183216998528140087E-11
      alp(1)=( 0.557503973136501826E+02,-0.204295038779771857E+03)
      rd(1)=(-0.562314417475317895E+01, 0.119406921611247440E+01)
      alp(2)=(-0.938666838877006739E+02, 0.912874896775456363E+02)
      rd(2)=(-0.508934679728216110E+01, 0.358882439228376881E+01)
      alp(3)=( 0.469965415550370835E+02,-0.116167609985818103E+02)
      rd(3)=(-0.399337136365302569E+01, 0.600483209099604664E+01)
      alp(4)=(-0.961424200626061065E+01,-0.264195613880262669E+01)
      rd(4)=(-0.226978543095856366E+01, 0.846173881758693369E+01)
      alp(5)=( 0.752722063978321642E+00, 0.670367365566377770E+00)
      rd(5)=( 0.208756929753827868E+00, 0.109912615662209418E+02)
      alp(6)=(-0.188781253158648576E-01,-0.343696176445802414E-01)
      rd(6)=( 0.370327340957595652E+01, 0.136563731924991884E+02)
      alp(7)=( 0.143086431411801849E-03, 0.287221133228814096E-03)
      rd(7)=( 0.889777151877331107E+01, 0.166309842834712071E+02)
c
c     if job .gt. 0 skip arnoldi process:
c
      if (job .gt. 0) goto 2
c
c  normalize vector w and put in first column of u --
c
      beta = sqrt( ddot (n,w,1,w,1))
      if (beta .eq. 0.0) then
         ierr = -1
         indic = 1
         return
      end if

      t = 1.0/beta
      do 25 j=1, n
         u(j,1) = w(j)*t
 25   continue
c  Arnoldi loop
      i1 = 1
 58   i = i1
      i1 = i + 1
      do 59 k=1, n
         x(k) = u(k,i)
 59   continue
      indic = 3
      return
 60   continue
      do 61 k=1, n
         u(k,i1) = y(k)
 61   continue
      i0 =1
c
c switch  for Lanczos version
c
c     i0 = max0(1, i-1)
      call mgsr (n, i0, i1, u, hh(1,i))
      fnorm = fnorm + ddot(i1, hh(1,i),1, hh(1,i),1)
      if (hh(i1,i) .eq. 0.0) m = i
      if  (i .lt. m) goto 58
c  done with arnoldi loop
      rm = dble(m)
      fnorm = sqrt( fnorm / rm )
c  get  : beta*e1 into z
      m1 = m+1
      do 4 i=1,m1
         hh(i,m1) = 0.0
 4    continue
c     compute initial dt when  job .lt.
      if (job .ge. 0) goto 2
c
c     t = eps / beta
c
      t = eps
      do 41 k=1, m-1
         t = t*(1.0 - dble(m-k)/rm )
 41   continue
      t = 2.0*rm* (t**(1.0/rm) )  / fnorm
c      dt = min(dt,t)
      t = min(abs(dt),t)
      dt = sign(t, dt)

 2    continue
      z(1) = beta
      do 3 k=2, m1
         z(k) = 0.0
 3    continue
c
c  get  : exp(H) * beta*e1
c
      call hes(ldg,m1,hh,ih,dt,z,rd,alp,alp0,wkc)
c  error estimate
      errst = abs(z(m1))

      indic = 2
      return
      end
