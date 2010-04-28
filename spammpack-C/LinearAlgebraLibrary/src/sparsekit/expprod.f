      subroutine expprod(n, m, eps, tn, u, w, x, y, a, ioff, ndiag)

c*********************************************************************72
c
cc EXPPROD computes an approximation to the vector
c
c              w :=  exp( - A * tn ) * w
c
c for matrices stored in diagonal (DIA) format.
c
c this routine constitutes an interface for the routine exppro for
c matrices stored in diagonal (DIA) format.
c
c ARGUMENTS
c
c see exppro for meaning of parameters n, m, eps, tn, u, w, x, y.
c
c a, ioff, and ndiag are the arguments of the matrix:
c
c a(n,ndiag) = a rectangular array with a(*,k) containing the diagonal
c              offset by ioff(k) (negative or positive or zero), i.e.,
c              a(i,jdiag) contains the element A(i,i+ioff(jdiag)) in
c              the usual dense storage scheme.
c
c ioff           = integer array containing the offsets  of the ndiag diagonals
c ndiag      = integer. the number of diagonals.
c
c
      double precision a(*), u(*), w(n), x(n), y(n), tn
      integer ioff(ndiag)

      indic = 0
 101  continue
      call exppro (n, m, eps, tn, u, w, x, y, indic,ierr)
      if (indic .eq. 1) goto 102
c
c     matrix vector-product for diagonal storage --
c
      call oped(n, x, y, a, ioff, ndiag)
      goto 101
 102  continue
      return
      end
