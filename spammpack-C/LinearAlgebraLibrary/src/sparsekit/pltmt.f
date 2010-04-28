      subroutine pltmt (nrow,ncol,mode,ja,ia,title,key,type,
     1     job, iounit)

c*********************************************************************72
c
cc PLTMT creates a 'pic' plot of a matrix.
c
c this subroutine creates a 'pic' file for plotting the pattern of
c a sparse matrix stored in general sparse format. It is not intended
c to be a means of plotting large matrices (It is very inefficient).
c It is however useful for small matrices and can be used for example
c for inserting matrix plots in a text. The size of the plot can be
c 7in x 7in or 5 in x 5in .. There is also an option for writing a
c 3-line header in troff (see description of parameter job).
c Author: Youcef Saad - Date: Sept., 1989
c See SPARSKIT/UNSUPP/ for a version of this to produce a post-script
c file.
c
c nrow   = number of rows in matrix
c
c ncol       = number of columns in matrix
c
c mode   = integer indicating whether the matrix is stored
c          row-wise (mode = 0) or column-wise (mode=1)
c
c ja     = column indices of nonzero elements when matrix is
c         stored rowise. Row indices if stores column-wise.
c ia     = integer array of containing the pointers to the
c         beginning of the columns in arrays a, ja.
c
c title  = character*71 = title of matrix test ( character a*71 ).
c key    = character*8  = key of matrix
c type   = character*3  = type of matrix.
c
c job    = this integer parameter allows to set a few minor
c          options. First it tells pltmt whether or not to
c          reduce the plot. The standard size of 7in is then
c          replaced by a 5in plot. It also tells pltmt whether or
c          not to append to the pic file a few 'troff' lines that
c          produce a centered caption includingg the title, key and
c          types as well as the size and number of nonzero elements.
c          job =  0 : do not reduce and do not make caption.
c          job =  1 : reduce and do not make caption.
c          job = 10 : do not reduce and make caption
c          job = 11 : reduce and make caption.
c          (i.e. trailing digit for reduction, leading digit for caption)
c
c iounit = logical unit number where to write the matrix into.
c
c example of usage .
c
c In the fortran code:
c  a) read a Harwell/Boeing matrix
c          call readmt (.....)
c         iout = 13
c  b) generate pic file:
c          call  pltmt (nrow,ncol,mode,ja,ia,title,key,type,iout)
c         stop
c
c Then in a unix environment plot the matrix by the command
c
c      pic FOR013.DAT | troff -me | lpr -Ppsx
c
c notes: 1) Plots square as well as rectangular matrices.
c            (however not as much tested with rectangular matrices.)
c        2) the dot-size is adapted according to the size of the
c            matrix.
c        3) This is not meant at all as a way of plotting large
c            matrices. The pic file generaled will have one line for
c            each nonzero element. It is  only meant for use in
c           such things as document poreparations etc..
c         4) The caption written will print the 71 character long
c            title. This may not be centered correctly if the
c            title has trailing blanks (a problem with Troff).
c            if you want the title centered then you can center
c            the string in title before calling pltmt.
c
      integer ja(*), ia(*)
      character key*8,title*72,type*3
      double precision x, y

      n = ncol
      if (mode .eq. 0) n = nrow
      nnz = ia(n+1) - ia(1)
      maxdim = max0 (nrow, ncol)
      xnrow = dble(nrow)
      xncol = dble(ncol)
      ptsize = 0.08
      hscale = (7.0 -2.0*ptsize)/dble(maxdim-1)
      vscale = hscale
      xwid  = ptsize + dble(ncol-1)*hscale + ptsize
      xht   = ptsize + dble(nrow-1)*vscale + ptsize
      xshift = (7.0-xwid)/2.0
      yshift = (7.0-xht)/2.0

      if (mod(job,10) .eq. 1) then
         write (iounit,88)
      else
         write (iounit,89)
      end if
 88   format('.PS 5in',/,'.po 1.8i')
 89   format('.PS',/,'.po 0.7i')
      write(iounit,90)
 90   format('box invisible wid 7.0 ht 7.0 with .sw at (0.0,0.0) ')
      write(iounit,91) xwid, xht, xshift, yshift
 91   format('box wid ',f5.2,' ht ',f5.2,
     *     ' with .sw at (',f5.2,',',f5.2,')' )
c
c     shift points slightly to account for size of dot , etc..
c
      tiny = 0.03
      if (mod(job,10) .eq. 1) tiny = 0.05
      xshift = xshift + ptsize - tiny
      yshift = yshift + ptsize + tiny

      ips = 8
      if (maxdim .le. 500) ips = 10
      if (maxdim .le. 300) ips = 12
      if (maxdim .le. 100) ips = 16
      if (maxdim .lt. 50) ips = 24
      write(iounit,92) ips
 92   format('.ps ',i2)
c
c  plottingloop
c
      do 1 ii=1, n
         istart = ia(ii)
         ilast  = ia(ii+1)-1
         if (mode .ne. 0) then
            x = dble(ii-1)
            do 2 k=istart, ilast
               y = xnrow-dble(ja(k))
               write(iounit,128) xshift+x*hscale, yshift+y*vscale
 2          continue
         else
            y = xnrow - dble(ii)
            do 3 k=istart, ilast
               x = dble(ja(k)-1)
               write(iounit,128) xshift+x*hscale, yshift+y*vscale
 3          continue
         end if
 1    continue

 128  format(7h"." at ,f6.3,',',f6.3,8h ljust  )
      write (iounit, 129)
 129  format('.PE')
c     quit if caption not desired.
      if ( (job/10) .ne. 1)  return

      write(iounit,127) key, type, title
      write(iounit,130) nrow,ncol,nnz
 127  format('.sp 4'/'.ll 7i'/'.ps 12'/'.po 0.7i'/'.ce 3'/,
     *     'Matrix:  ',a8,',  Type:  ',a3,/,a72)
 130  format('Dimension: ',i4,' x ',i4,',  Nonzero elements: ',i5)
      return
      end
