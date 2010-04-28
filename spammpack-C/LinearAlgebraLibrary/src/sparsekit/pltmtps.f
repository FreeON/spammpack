      subroutine pltmtps (nrow,ncol,mode,ja,ia,title,key,type,
     1                        job, iounit)

c*********************************************************************72
c
cc PLTMTPS creates a PostScript plot of a sparse matrix.
c
c this subroutine creates a 'PS' file for plotting the pattern of
c a sparse matrix stored in general sparse format. It can be used
c for inserting matrix plots in a text. The size of the plot can be
c 7in x 7in or 5 in x 5in ..
c
c Adapted from pltmt in module INOUT by Paul Frederickson. March, 1990
c         + slight modifications by Y. Saad.
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
c title  = character*72 = title of matrix test ( character a*72 ).
c key    = character*8  = key of matrix
c type   = character*3  = type of matrix.
c
c job    =  integer. tells pltmt whether or not to reduce the plot.
c           if enabled then the standard size of 7in will be
c           replaced by a 5in plot.
c          job =  0 : do not reduce
c          job =  1 : reduce plot to 5 inches.
c
c iounit = logical unit number where to write the matrix into.
c
c notes: 1) Plots square as well as rectangular matrices.
c        2) Does not writer a caption yet.
c        3) No bounding box put in yet
c
      integer ja(*), ia(*)
      character key*8,title*72,type*3
      double precision x, y, delta

      n = ncol
      if (mode .eq. 0) n = nrow
      nnz = ia(n+1) - ia(1)
      maxdim = max0 (nrow, ncol)
      m = 1 + maxdim
c keep this test as in old pltmt (for future changes).
       if (mod(job,10) .eq. 1) then
        delta = 72*5.0/(2.0+maxdim)
       else
        delta = 72*7.0/(2.0+maxdim)
       end if

      write(iounit,*)'%!PS'
       write(iounit,*)' gsave 50 50 translate'
      write(iounit,*) delta, delta, ' scale'
      write(iounit,*) ' 0.25 setlinewidth'

       if (mod(job,10) .eq. 1) then
          write (iounit,*) ' 23 55 translate'
       else
          write (iounit,*) ' 2 35 translate'
       end if

      write(iounit,*) ' newpath'
      write(iounit,*) 0,0,' moveto'
      write(iounit,*) m,0,' lineto'
      write(iounit,*) m,m,' lineto'
      write(iounit,*) 0,m,' lineto'
      write(iounit,*) ' closepath stroke'
      write(iounit,*) ' 1 1 translate'
      write(iounit,*) ' 0.5 setlinewidth'
      write(iounit,*) ' /p {moveto 0 -.25 rmoveto '
      write(iounit,*) '            0  .50 rlineto stroke} def'
c
c  plotting loop
c
           do 1 ii=1, n
           istart = ia(ii)
           ilast  = ia(ii+1)-1
             if (mode .ne. 0) then
                do 2 k=istart, ilast
            write(iounit,*) ii-1, nrow-ja(k), ' p'
 2              continue
          else
c             y = xnrow - dble(ii)
              do 3 k=istart, ilast
c               x = dble(ja(k)-1)
            write(iounit,*) ja(k)-1, nrow-ii, ' p'
 3              continue
            end if
 1      continue

      write(iounit,*)' showpage grestore'
c
c quit if caption not desired.
c      if ( (job/10) .ne. 1)  return
c
c      write(iounit,127) key, type, title
c      write(iounit,130) nrow,ncol,nnz
c 127         format('.sp 4'/'.ll 7i'/'.ps 12'/'.po 0.7i'/'.ce 3'/,
c     *        'Matrix:  ',a8,',  Type:  ',a3,/,a71)
 130    format('Dimension: ',i4,' x ',i4',  Nonzero elements: ',i5)
      return
      end
