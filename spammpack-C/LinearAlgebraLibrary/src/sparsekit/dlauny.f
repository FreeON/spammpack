      subroutine dlauny(x,y,nodes,elmnts,nemax,nelmnt)

c*********************************************************************72
c
cc DLAUNY is a simple, nonoptimal Delaunay triangulation code.
c
c code written by P.K. Sweby
c
c     Performs a Delaunay triangularisation of a region given a set
c     of mesh points.
c       X,Y    :- 1D arrays holding coordinates of mesh points.
c                 dimensioned AT LEAST NODES+3.
c       NODES  :- number of mesh points.
c       ELMNTS :- INTEGER array, dimensioned NEMAX x 3, which on exit
c                 contains the index of global nodes associated with
c                 each element.
c       NELMNT :- on exit contains the number of elements in the
c                 triangularisation.
c
c                                       P.K.Sweby
c
      IMPLICIT double precision (A-H,O-Z)

      INTEGER ELMNTS
      DIMENSION X(NODES),Y(NODES),ELMNTS(NEMAX,3)

      PI=4.0*ATAN(1.0)
c
c     Calculate artificial nodes NODES+i i=1,2,3,4 and construct first
c     two (artificial) elements.
c
      XMIN=X(1)
      XMAX=X(1)
      YMIN=Y(1)
      YMAX=Y(1)
      DO 10 I=2,NODES
      XMIN=MIN(XMIN,X(I))
      XMAX=MAX(XMAX,X(I))
      YMIN=MIN(YMIN,Y(I))
      YMAX=MAX(YMAX,Y(I))
 10   CONTINUE
      DX=XMAX-XMIN
      DY=YMAX-YMIN
      XL=XMIN-4.0*DX
      XR=XMAX+4.0*DX
      YL=YMIN-4.0*DY
      YR=YMAX+4.0*DY
      X(NODES+1)=XL
      Y(NODES+1)=YL
      X(NODES+2)=XL
      Y(NODES+2)=YR
      X(NODES+3)=XR
      Y(NODES+3)=YR
      X(NODES+4)=XR
      Y(NODES+4)=YL
      ELMNTS(1,1)=NODES+1
      ELMNTS(1,2)=NODES+2
      ELMNTS(1,3)=NODES+3
      ELMNTS(2,1)=NODES+3
      ELMNTS(2,2)=NODES+4
      ELMNTS(2,3)=NODES+1
      NELMNT=2
      DO 90 IN=1,NODES
c
c     Add one mesh point at a time and remesh locally if necessary
c
      NDEL=0
      NEWEL=0
      DO 40 IE=1,NELMNT
c
c     Is point IN insided circumcircle of element IE ?
c
      I1=ELMNTS(IE,1)
      I2=ELMNTS(IE,2)
      I3=ELMNTS(IE,3)
      X2=X(I2)-X(I1)
      X3=X(I3)-X(I1)
      Y2=Y(I2)-Y(I1)
      Y3=Y(I3)-Y(I1)
      Z=(X2*(X2-X3)+Y2*(Y2-Y3))/(Y2*X3-Y3*X2)
      CX=0.5*(X3-Z*Y3)
      CY=0.5*(Y3+Z*X3)
      R2=CX**2+CY**2
      RN2=((X(IN)-X(I1)-CX)**2+(Y(IN)-Y(I1)-CY)**2)
      IF(RN2.GT.R2)GOTO 40
c
c     Yes it is inside,create new elements and mark old for deletion.
c
      DO 30 J=1,3
      DO 20 K=1,3
      ELMNTS(NELMNT+NEWEL+J,K)=ELMNTS(IE,K)
 20   CONTINUE
      ELMNTS(NELMNT+NEWEL+J,J)=IN
 30   CONTINUE
      NEWEL=NEWEL+3
      ELMNTS(IE,1)=0
      NDEL=NDEL+1

 40   CONTINUE
c
c     If IN was inside circumcircle of more than 1 element then will
c     have created 2 identical new elements: delete them both.
c
      IF(NDEL.GT.1)THEN
          DO 60 IE=NELMNT+1,NELMNT+NEWEL-1
          DO 60 JE=IE+1,NELMNT+NEWEL
          MATCH=0
          DO 50 K=1,3
          DO 50 L=1,3
          IF(ELMNTS(IE,K).EQ.ELMNTS(JE,L))MATCH=MATCH+1
 50       CONTINUE
          IF(MATCH.EQ.3)THEN
              ELMNTS(IE,1)=0
              ELMNTS(JE,1)=0
              NDEL=NDEL+2
          end if
 60       CONTINUE
      end if
c
c     Delete any elements
c
      NN=NELMNT+NEWEL
      IE=1
 70   CONTINUE
      IF(ELMNTS(IE,1).EQ.0)THEN
          DO 80 J=IE,NN-1
          DO 80 K=1,3
          ELMNTS(J,K)=ELMNTS(J+1,K)
 80       CONTINUE
          NN=NN-1
          IE=IE-1
      end if
      IE=IE+1
      IF(IE.LE.NN)GOTO 70
      NELMNT=NN
 90   CONTINUE
c
c     Finally remove elements containing artificial nodes
c
      IE=1
 100  CONTINUE
      NART=0
      DO 110 L=1,3
      IF(ELMNTS(IE,L).GT.NODES)NART=NART+1
 110  CONTINUE
      IF(NART.GT.0)THEN
          DO 120 J=IE,NN-1
          DO 120 K=1,3
          ELMNTS(J,K)=ELMNTS(J+1,K)
 120      CONTINUE
          NELMNT=NELMNT-1
          IE=IE-1
      end if
      IE=IE+1
      IF(IE.LE.NELMNT)GOTO 100
      RETURN
      END
