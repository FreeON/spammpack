      subroutine pgmres(n, im, rhs, sol, vv, eps, maxits, iout,
     *                    aa, ja, ia, alu, jlu, ju, ierr)

c*********************************************************************72
c
cc PGMRES is an ILUT - Preconditioned GMRES solver.
c
c                 *** ILUT - Preconditioned GMRES ***
c
c*
c This is a simple version of the ILUT preconditioned GMRES algorithm.
c The ILUT preconditioner uses a dual strategy for dropping elements
c instead  of the usual level of-fill-in approach. See details in ILUT
c subroutine documentation. PGMRES uses the L and U matrices generated
c from the subroutine ILUT to precondition the GMRES algorithm.
c The preconditioning is applied to the right. The stopping criterion
c utilized is based simply on reducing the residual norm by epsilon.
c This preconditioning is more reliable than ilu0 but requires more
c storage. It seems to be much less prone to difficulties related to
c strong nonsymmetries in the matrix. We recommend using a nonzero tol
c (tol=.005 or .001 usually give good results) in ILUT. Use a large
c lfil whenever possible (e.g. lfil = 5 to 10). The higher lfil the
c more reliable the code is. Efficiency may also be much improved.
c Note that lfil=n and tol=0.0 in ILUT  will yield the same factors as
c Gaussian elimination without pivoting.
c
c ILU(0) and MILU(0) are also provided for comparison purposes
c USAGE: first call ILUT or ILU0 or MILU0 to set up preconditioner and
c then call pgmres.
c
c Coded by Y. Saad - This version dated May, 7, 1990.
c*
c parameters
c
c on entry:
c
c
c n     == integer. The dimension of the matrix.
c im    == size of krylov subspace:  should not exceed 50 in this
c          version (can be reset by changing parameter command for
c          kmax below)
c rhs   == double precision vector of length n containing the right hand side.
c          Destroyed on return.
c sol   == double precision vector of length n containing an initial guess to the
c          solution on input. approximate solution on output
c eps   == tolerance for stopping criterion. process is stopped
c          as soon as ( ||.|| is the euclidean norm):
c          || current residual||/||initial residual|| <= eps
c maxits== maximum number of iterations allowed
c iout  == output unit number number for printing intermediate results
c          if (iout .le. 0) nothing is printed out.
c
c aa, ja,
c ia    == the input matrix in compressed sparse row format:
c          aa(1:nnz)  = nonzero elements of A stored row-wise in order
c          ja(1:nnz) = corresponding column indices.
c          ia(1:n+1) = pointer to beginning of each row in aa and ja.
c          here nnz = number of nonzero elements in A = ia(n+1)-ia(1)
c
c alu,jlu== A matrix stored in Modified Sparse Row format containing
c           the L and U factors, as computed by subroutine ilut.
c
c ju     == integer array of length n containing the pointers to
c           the beginning of each row of U in alu, jlu as computed
c           by subroutine ILUT.
c
c on return:
c
c sol   == contains an approximate solution (upon successful return).
c ierr  == integer. Error message with the following meaning.
c          ierr = 0 --> successful return.
c          ierr = 1 --> convergence not achieved in itmax iterations.
c          ierr =-1 --> the initial guess seems to be the exact
c                       solution (initial residual computed was zero)
c
c
c
c work arrays:
c
c vv    == work array of length  n x (im+1) (used to store the Arnoli
c          basis)
c
c subroutines called :
c ope    : matrix by vector multiplication delivers y=ax, given x
c lusol0 : combined forward and backward solves (Preconditioning ope.)
c BLAS2  routines.
c
       implicit double precision (a-h,o-z)
       integer n, im, maxits, iout, ierr, ja(*), ia(n+1), jlu(*), ju(n)
       double precision vv(n,*), rhs(n), sol(n), aa(*), alu(*), eps

       parameter (kmax=50)
      PARAMETER (EPSMAC=1.0E-16)
c
       double precision hh(kmax+1,kmax), c(kmax), s(kmax), rs(kmax+1),t
c
c arnoldi size should not exceed kmax=50 in this version..
c to reset modify paramter kmax accordingly.
c
       n1 = n + 1
       its = 0
c
c outer loop starts here..
c  compute initial residual vector
       call ope (n, sol, vv, aa, ja, ia)
       do 21 j=1,n
          vv(j,1) = rhs(j) - vv(j,1)
 21    continue

 20    ro = sqrt( ddot(n, vv, 1, vv, 1) )
       if (iout .gt. 0 .and. its .eq. 0)
     *      write(iout, 199) its, ro
       if (ro .eq. 0.0) goto 999
       t = 1.0/ ro
       do 210 j=1, n
          vv(j,1) = vv(j,1)*t
 210   continue
       if (its .eq. 0) eps1=eps*ro
c     ** initialize 1-st term  of rhs of hessenberg system..
       rs(1) = ro
       i = 0
 4     i=i+1
       its = its + 1
       i1 = i + 1
       call lusol0 (n, vv(1,i), rhs, alu, jlu, ju)
       call ope (n, rhs, vv(1,i1), aa, ja, ia)
c
c     modified gram - schmidt...
c
       do 55 j=1, i
          t = ddot(n, vv(1,j),1,vv(1,i1),1)
          hh(j,i) = t
          call daxpy(n, -t, vv(1,j), 1, vv(1,i1), 1)
 55    continue
       t = sqrt(ddot(n, vv(1,i1), 1, vv(1,i1), 1))
       hh(i1,i) = t
       if ( t .eq. 0.0) goto 58
       t = 1.0/t
       do 57  k=1,n
          vv(k,i1) = vv(k,i1)*t
 57    continue
c
c     done with modified gram schimd and arnoldi step..
c     now  update factorization of hh
c
 58    if (i .eq. 1) goto 121
c
c  perform previous transformations  on i-th column of h
c
       do 66 k=2,i
          k1 = k-1
          t = hh(k1,i)
          hh(k1,i) = c(k1)*t + s(k1)*hh(k,i)
          hh(k,i) = -s(k1)*t + c(k1)*hh(k,i)
 66    continue
 121   gam = sqrt(hh(i,i)**2 + hh(i1,i)**2)
c
c     if gamma is zero then any small value will do...
c     will affect only residual estimate
c
       if (gam .eq. 0.0) gam = epsmac
c
c     get  next plane rotation
c
       c(i) = hh(i,i)/gam
       s(i) = hh(i1,i)/gam
       rs(i1) = -s(i)*rs(i)
       rs(i) =  c(i)*rs(i)
c
c     detrermine residual norm and test for convergence-
c
       hh(i,i) = c(i)*hh(i,i) + s(i)*hh(i1,i)
       ro = abs(rs(i1))
 131   format(1h ,2e14.4)
       if (iout .gt. 0)
     *      write(iout, 199) its, ro
       if (i .lt. im .and. (ro .gt. eps1))  goto 4
c
c     now compute solution. first solve upper triangular system.
c
       rs(i) = rs(i)/hh(i,i)
       do 30 ii=2,i
          k=i-ii+1
          k1 = k+1
          t=rs(k)
          do 40 j=k1,i
             t = t-hh(k,j)*rs(j)
 40       continue
          rs(k) = t/hh(k,k)
 30    continue
c
c     form linear combination of v(*,i)'s to get solution
c
       t = rs(1)
       do 15 k=1, n
          rhs(k) = vv(k,1)*t
 15    continue
       do 16 j=2, i
          t = rs(j)
          do 161 k=1, n
             rhs(k) = rhs(k)+t*vv(k,j)
 161      continue
 16    continue
c
c     call preconditioner.
c
       call lusol0 (n, rhs, rhs, alu, jlu, ju)
       do 17 k=1, n
          sol(k) = sol(k) + rhs(k)
 17    continue
c
c     restart outer loop  when necessary
c
       if (ro .le. eps1) goto 990
       if (its .gt. maxits) goto 991
c
c     else compute residual vector and continue..
c
       do 24 j=1,i
          jj = i1-j+1
          rs(jj-1) = -s(jj-1)*rs(jj)
          rs(jj) = c(jj-1)*rs(jj)
 24    continue
       do 25  j=1,i1
          t = rs(j)
          if (j .eq. 1)  t = t-1.0
          call daxpy (n, t, vv(1,j), 1,  vv, 1)
 25    continue
 199   format(' its =', i4, ' res. norm =', G14.6)
c     restart outer loop.
       goto 20
 990   ierr = 0
       return
 991   ierr = 1
       return
 999   continue
       ierr = -1
       return
       end
