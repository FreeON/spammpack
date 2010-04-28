      subroutine prtmt (nrow,ncol,a,ja,ia,rhs,guesol,title,key,type,
     1     ifmt,job,iounit)

c*********************************************************************72
c
cc PRTMT writes a matrix in Harwell-Boeing format into a file.
c
c writes a matrix in Harwell-Boeing format into a file.
c assumes that the matrix is stored in COMPRESSED SPARSE COLUMN FORMAT.
c some limited functionality for right hand sides.
c Author: Youcef Saad - Date: Sept., 1989 - updated Oct. 31, 1989 to
c cope with new format.
c
c on entry:
c
c nrow   = number of rows in matrix
c ncol       = number of columns in matrix
c a       = double precision array containing the values of the matrix stored
c          columnwise
c ja        = integer array of the same length as a containing the row indices
c          of the corresponding matrix elements of array a.
c ia     = integer array of containing the pointers to the beginning of
c         the columns in arrays a, ja.
c rhs    = double precision array  containing the right hand side (s) and optionally
c          the associated initial guesses and/or exact solutions
c          in this order. See also guesol for details. the vector rhs will
c          be used only if job .gt. 2 (see below). Only full storage for
c          the right hand sides is supported.
c
c guesol = a 2-character string indicating whether an initial guess
c          (1-st character) and / or the exact solution (2-nd)
c          character) is provided with the right hand side.
c         if the first character of guesol is 'G' it means that an
c          an intial guess is provided for each right hand sides.
c          These are assumed to be appended to the right hand sides in
c          the array rhs.
c         if the second character of guesol is 'X' it means that an
c          exact solution is provided for each right hand side.
c          These are assumed to be appended to the right hand sides
c          and the initial guesses (if any) in the array rhs.
c
c title  = character*71 = title of matrix test ( character a*71 ).
c key    = character*8  = key of matrix
c type   = charatcer*3  = type of matrix.
c
c ifmt       = integer specifying the format chosen for the real values
c         to be output (i.e., for a, and for rhs-guess-sol if
c          applicable). the meaning of ifmt is as follows.
c        * if (ifmt .lt. 100) then the E descriptor is used,
c           format Ed.m, in which the length (m) of the mantissa is
c           precisely the integer ifmt (and d = ifmt+6)
c        * if (ifmt .gt. 100) then prtmt will use the
c           F- descriptor (format Fd.m) in which the length of the
c           mantissa (m) is the integer mod(ifmt,100) and the length
c           of the integer part is k=ifmt/100 (and d = k+m+2)
c          Thus  ifmt= 4   means  E10.4  +.xxxxD+ee    while
c                ifmt=104  means  F7.4   +x.xxxx
c                ifmt=205  means  F9.5   +xx.xxxxx
c          Note: formats for ja, and ia are internally computed.
c
c job       = integer to indicate whether matrix values and
c         a right hand side is available to be written
c          job = 1   write srtucture only, i.e., the arrays ja and ia.
c          job = 2   write matrix including values, i.e., a, ja, ia
c          job = 3   write matrix and one right hand side: a,ja,ia,rhs.
c         job = nrhs+2 write matrix and nrhs successive right hand sides
c         Note that there cannot be any right hand side if the matrix
c         has no values. Also the initial guess and exact solutions when
c          provided are for each right hand side. For example if nrhs=2
c          and guesol='GX' there are 6 vectors to write.
c
c
c iounit = logical unit number where to write the matrix into.
c
c on return:
c
c the matrix a, ja, ia will be written in output unit iounit
c in the Harwell-Boeing format. Noe of the inputs is modofied.
c
c Notes: 1) This code attempts to pack as many elements as possible per
c        80-character line.
c        2) this code attempts to avoid as much as possible to put
c        blanks in the formats that are written in the 4-line header
c       (This is done for purely esthetical reasons since blanks
c        are ignored in format descriptors.)
c        3) sparse formats for righr hand sides and guesses not
c        suported.
c
      character title*72,key*8,type*3,ptrfmt*16,indfmt*16,valfmt*20,
     *              guesol*2, rhstyp*3
      integer totcrd, ptrcrd, indcrd, valcrd, rhscrd, nrow, ncol,
     1     nnz, nrhs, len, nperli
      integer ja(*), ia(*)
      double precision a(*),rhs(*)
c
c     compute pointer format
c
      nnz    = ia(ncol+1) -1
      len    = int ( dlog10(0.1+dble(nnz+1))) + 1
      nperli = 80/len
      ptrcrd = ncol/nperli + 1
      if (len .gt. 9) then
         assign 101 to ix
      else
         assign 100 to ix
      end if
      write (ptrfmt,ix) nperli,len
 100  format(1h(,i2,1HI,i1,1h) )
 101  format(1h(,i2,1HI,i2,1h) )
c
c compute ROW index format
c
      len    = int ( dlog10(0.1+dble(nrow) )) + 1
      nperli = min0(80/len,nnz)
      indcrd = (nnz-1)/nperli+1
      write (indfmt,100) nperli,len
c
c compute values and rhs format (using the same for both)
c
      valcrd      = 0
      rhscrd  = 0
c quit this part if no values provided.
      if (job .le. 1) goto 20

      if (ifmt .ge. 100) then
         ihead = ifmt/100
         ifmt = ifmt-100*ihead
         len = ihead+ifmt+2
         nperli = 80/len

         if (len .le. 9 ) then
            assign 102 to ix
         elseif (ifmt .le. 9) then
            assign 103 to ix
         else
            assign 104 to ix
         end if

         write(valfmt,ix) nperli,len,ifmt
 102     format(1h(,i2,1hF,i1,1h.,i1,1h) )
 103     format(1h(,i2,1hF,i2,1h.,i1,1h) )
 104     format(1h(,i2,1hF,i2,1h.,i2,1h) )

      else
         len = ifmt + 6
         nperli = 80/len
c     try to minimize the blanks in the format strings.
         if (nperli .le. 9) then
          if (len .le. 9 ) then
             assign 105 to ix
          elseif (ifmt .le. 9) then
             assign 106 to ix
          else
             assign 107 to ix
          end if
       else
          if (len .le. 9 ) then
             assign 108 to ix
          elseif (ifmt .le. 9) then
             assign 109 to ix
          else
               assign 110 to ix
            end if
         end if

         write(valfmt,ix) nperli,len,ifmt
 105     format(1h(,i1,1hE,i1,1h.,i1,1h) )
 106     format(1h(,i1,1hE,i2,1h.,i1,1h) )
 107     format(1h(,i1,1hE,i2,1h.,i2,1h) )
 108     format(1h(,i2,1hE,i1,1h.,i1,1h) )
 109     format(1h(,i2,1hE,i2,1h.,i1,1h) )
 110     format(1h(,i2,1hE,i2,1h.,i2,1h) )

      end if
      valcrd = (nnz-1)/nperli+1
      nrhs   = job -2
      if (nrhs .ge. 1) then
         i = (nrhs*nrow-1)/nperli+1
         rhscrd = i
         if (guesol(1:1) .eq. 'G') rhscrd = rhscrd+i
         if (guesol(2:2) .eq. 'X') rhscrd = rhscrd+i
         rhstyp = 'F'//guesol
      end if
 20   continue
c
      totcrd = ptrcrd+indcrd+valcrd+rhscrd
c     write 4-line or five line header
      write(iounit,10) title,key,totcrd,ptrcrd,indcrd,valcrd,
     1     rhscrd,type,nrow,ncol,nnz,nrhs,ptrfmt,indfmt,valfmt,valfmt
c
      if (nrhs .ge. 1) write (iounit,11) rhstyp, nrhs
 10   format (a72, a8 / 5i14 / a3, 11x, 4i14 / 2a16, 2a20)
 11   format(A3,11x,i4)
c
      write(iounit,ptrfmt) (ia (i), i = 1, ncol+1)
      write(iounit,indfmt) (ja (i), i = 1, nnz)
      if (job .le. 1) return
      write(iounit,valfmt) (a(i), i = 1, nnz)
      if (job .le. 2) return
      len = nrow*nrhs
      next = 1
      iend = len
      write(iounit,valfmt) (rhs(i), i = next, iend)
c
c     write initial guesses if available
c
      if (guesol(1:1) .eq. 'G') then
         next = next+len
         iend = iend+ len
         write(iounit,valfmt) (rhs(i), i = next, iend)
      end if
c
c     write exact solutions if available
c
      if (guesol(2:2) .eq. 'X')then
         next = next+len
         iend = iend+ len
         write(iounit,valfmt) (rhs(i), i = next, iend)
      end if

      return
      end
