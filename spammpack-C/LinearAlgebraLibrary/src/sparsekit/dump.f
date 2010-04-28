      subroutine dump(n,a,ja,ia,iout)

c*********************************************************************72
c
cc DUMP writes the matrix to a file.
c
c writes the matrix in a file, one row at a time in a nice readable
c format. This is a simple routine which is useful for debugging.
c
c on entry:
c
c n     = integer = size of matrix
c a,
c ja,
c ia    =  matrix in CSR format
c iout  = output unit number.
c
c on return:
c
c the output file iout will have written in it the matrix in
c one of two possible formats (depending on the max number of
c elements per row. the values are output with only two digits
c of accuracy (D9.2).
c
      integer ia(*),ja(*)
      double precision a(*)
c
c select mode horizontal or vertical
c
        maxr = 0
        do 1 i=1, n
           maxr = max0(maxr,ia(i+1)-ia(i))
 1      continue
        if (maxr .le. 8) then
c
c able to one row across line
c
      do 2 i=1, n
           write(iout,100) i
         k1=ia(i)
         k2 = ia(i+1)-1
         write (iout,101) (ja(k),k=k1,k2)
         write (iout,102) (a(k),k=k1,k2)
 2      continue
      else
c
c unable to one row across line. do three items at a time across line.
c
         do 3 i=1, n
            write(iout,200) i
            k1=ia(i)
            k2 = ia(i+1)-1
            write (iout,201) (ja(k),a(k),k=k1,k2)
 3       continue
      end if

 100  format(1X,35(1h-),' row',i3,1x,35(1h-) )
 101  format(' col:',8(i5,6h     :))
 102  format(' val:',8(E9.2,2h :) )
 200  format(1h ,31(1h-),' row',i3,1x,31(1h-),/
     *       3('  columns :   values   *') )
 201  format(3(1h ,i5,6h    : ,D9.2,3h  *) )
      return
      end
