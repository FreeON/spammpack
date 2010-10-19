      subroutine exppro(n, m, eps, tn, u, w, x, y, indic, ierr)

c*********************************************************************72
c
cc EXPPRO computes an approximation to the vector
c
c              w :=  exp( - A * tn ) * w
c
c where A is an arbitary matrix and w is a given input vector
c uses a dynamic estimation of internal time advancement (dt)
c
c THIS IS A REVERSE COMMUNICATION IMPLEMENTATION.
c
c USAGE: (see also comments on indic below).
c
c
c      indic = 0
c 1    continue
c      call exppro (n, m, eps, tn, u, w, x, y, indic)
c      if (indic .eq. 1) goto 2 <-- indic .eq.1 means job is finished
c      call matvec(n, x, y)     <--- user's matrix-vec. product
c                                    with x = input vector, and
c                                     y = result = A * x.
c      goto 1
c 2    continue
c      .....
c
c en entry:
c
c n      = dimension of matrix
c
c m      = dimension of Krylov subspace (= degree of polynomial
c         approximation to the exponential used. )
c
c eps   = scalar indicating the relative error tolerated for the result.
c         the code will try to compute an answer such that
c         norm2(exactanswer-approximation) / norm2(w) .le. eps
c
c tn      = scalar by which to multiply matrix. (may be .lt. 0)
c         the code will compute an approximation to exp(- tn * A) w
c         and overwrite the result onto w.
c
c u      = work array of size n*(m+1) (used to hold the Arnoldi basis )
c
c w      = double precision array of length n = input vector to  which exp(-A) is
c         to be applied. this is also an output argument
c
c x, y  = two double precision work vectors of length at least  n each.
c         see indic for usage.
c
c indic = integer used as indicator for the reverse communication.
c         in the first call enter indic = 0. See below for more.
c
c on return:
c
c w     = contains the resulting vector exp(-A * tn ) * w when
c         exppro has finished (see indic)
c
c indic = indicator for the reverse communication protocole.
c       * INDIC .eq. 1  means that exppro has finished and w contains the
c         result.
c       * INDIC .gt. 1 ,  means that exppro has not finished and that
c         it is requesting another matrix vector product before
c         continuing. The user must compute Ax where A is the matrix
c         and x is the vector provided by exppro, and return the
c         result in y. Then exppro must be called again without
c         changing any other argument. typically this must be
c         implemented in a loop with exppro being called as long
c         indic is returned with a value .ne. 1.
c
c ierr  = error indicator.
c         ierr = 1 means phipro was called with indic=1 (not allowed)
c         ierr = -1 means that the input is zero the solution has been
c         unchanged.
c
c NOTES:  im should not exceed 60 in this version  (see ih0 below)
c
c written by Y. Saad -- version feb, 1991.
c
c For reference see fololowing papers :
c (1) E. Gallopoulos and Y. Saad: Efficient solution of parabolic
c     equations by Krylov approximation methods. RIACS technical
c     report 90-14.
c (2) Y.Saad: Analysis of some Krylov subspace approximations to the
c     matrix exponential operator. RIACS Tech report. 90-14
c
      integer n, m, indic, ierr
      double precision eps, tn, u(*), w(n), x(*), y(*)
      parameter (ih0=60)
      double precision hh(ih0,ih0)
      double precision z(ih0),errst, tcur, told, dtl, beta, red
      double complex wkc(ih0)
      save
c
c indic = 3  means  passing through only with result of y= Ax to exphes
c indic = 2  means exphes has finished its job
c indic = 1  means exppro has finished its job (real end)/
c
      ierr = 0
      if (indic .eq. 3) goto 101
      if (indic .eq. 1) then
         ierr = 1
         return
      end if
      ih = ih0
      m  = min(m,ih0)
      tcur = 0.0
      dtl = tn-tcur
      job = -1
c outer loop
 100  continue
c
c  call exponential propagator
c
      told = tcur
 101  continue
c     if (told + dtl .gt. tn) dtl = tn-told
      call  exphes (n,m,dtl,eps,u,w,job,z,wkc,beta,errst,hh,ih,
     *              x,y,indic,ierr)

      if (ierr .ne. 0) return
      if (indic .ge. 3) return
      tcur = told + dtl
c
c     relative error
c
      errst = errst / beta

      if ((errst .le. eps) .and. ( (errst .gt. eps/100.0) .or.
     *     (tcur .eq. tn))) goto 102
c
c     use approximation :  [ new err ] = fact**m  * [cur. error]
c
      red =  (0.5*eps / errst)**(1.0 /dble(m) )
      dtl = dtl*red
      if (abs(told+dtl) .gt. abs(tn) )  dtl = tn-told
      job = 1
      goto 101

 102  continue
      call project(n,m,u,z,w)
      job = 0
      dtl = min(dtl, tn-tcur)
      if (abs(tcur+dtl) .gt. abs(tn)) dtl = tn-tcur
      if (abs(tcur) .lt. abs(tn)) goto 100
      indic = 1

      return
      end
