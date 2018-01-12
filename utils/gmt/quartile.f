      program boxwhiskermonthly

      parameter (ik=1000)
      integer idata(ik)
      real rdata(ik),rquar(5)


      open (100)
        open (22)
          do i=1,ik
            read (100,*,end=201) ic,jk,(rdata(k),k=1,jk)
            do k=1,jk
              idata(k)=INT(rdata(k)*100.)
c              print*,k,idata(k)
            enddo
            call quartile(jk,idata,rquar,rmittel)
            write (22,'(I5,5(1x,F11.2))') ic,rquar(3)/100.,rquar(1)/100.,rquar(2)/100.,rquar(4)/100.,rquar(5)/100.
          enddo
 201      continue
        close (22)
      close (100)

      stop
      end

************************************************************************************


************************************************************************************

      subroutine quartile(ik,idata,rquar,rmittel)

      integer iord(ik),idata(ik)
      real x(ik),x2(ik),rquar(5)

      isum=0
      do k=1,ik
        isum=isum+idata(k)
      enddo
      rmittel=FLOAT(isum)/FLOAT(ik)
      
c Daten sortieren
      call qsorti(iord,ik,idata)
      do k=1,ik
        x(k)=FLOAT(idata(iord(k)))
      enddo

c Quartile berechnen
c 1: minimaler Wert (lowest value)
c 2: 1. Quartil = 25. Perzentil (first quartile = 25th percentile)
c 3: 2. Quartil = 50. Perzentil = Median (second quartile = 50th percentile = median)
c 4: 3. Quartil = 75. Perzentil (third quartil = 75th percentile)
c 5: maximaler Wert (highest value)
c There is no universal agreement on choosing the quartile values.
c One possible rule is as follows:
c 1. Use the median to divide the ordered data set into two halves.
c    Do not include the median into the halves.
c 2. The lower quartile value is the median of the lower half of the data.
c    The upper quartile value is the median of the upper half of the data.

      do j=1,5
        rquar(j)=-9.0
      enddo
      if (ik.eq.1) then
        rquar(3)=x(1)
      elseif (ik.eq.2) then
        rquar(1)=x(1)
        rquar(5)=x(2)
      elseif (ik.eq.3) then
        rquar(1)=x(1)
        rquar(3)=x(2)
        rquar(5)=x(3)
      else
        rquar(1)=x(1)
        rquar(5)=x(ik)
        call median(x,ik,xmed)
        rquar(3)=xmed
c i1: 1. Haelfte der Daten ohne den Wert des Medians 
c Bsp. 1: 7 Werte: INT(7./2.)=3, 4. Wert=Median)
c Bsp. 2: 8 Werte: INT(8./2.)=4, Median: Mittel zw. 4. und 5. Wert)
        i1=INT(FLOAT(ik)/2.)
        call median(x,i1,xmed)
        rquar(2)=xmed
c i1: 2. Haelfte der Daten ohne den Wert des Medians
c Bsp. 1: 7 Werte (ungerade): INT(7./2.)+2=5, 4. Wert=Median)
c Bsp. 2: 8 Werte: INT(8./2.)+1=5, Median: Mittel zw. 4. und 5. Wert)
        if (MOD(ik,2).eq.1) then
          i3=INT(FLOAT(ik)/2.)+2
        else
          i3=INT(FLOAT(ik)/2.)+1
        endif
        il=0
        do j=i3,ik
          il=il+1
          x2(il)=x(j)
        enddo
        call median(x2,il,xmed)
        rquar(4)=xmed
      endif

      return
      end

*************************************************************************************




*************************************************************************************

      SUBROUTINE QSORTI (ORD,N,A)

c ORD geordnetes integer Feld A
c n Felder ORD und A haben diese Dimension
c A ungeordnetes integer Feld


C
C==============SORTS THE ARRAY A(I),I=1,2,...,N BY PUTTING THE
C   ASCENDING ORDER VECTOR IN ORD.  THAT IS ASCENDING ORDERED A
C   IS A(ORD(I)),I=1,2,...,N; DESCENDING ORDER A IS A(ORD(N-I+1)),
C   I=1,2,...,N .  THIS SORT RUNS IN TIME PROPORTIONAL TO N LOG N .
C
C
C     ACM QUICKSORT - ALGORITHM #402 - IMPLEMENTED IN FORTRAN 66 BY
C                                 WILLIAM H. VERITY, WHV@PSUVM.PSU.EDU
C                                 CENTER FOR ACADEMIC COMPUTING
C                                 THE PENNSYLVANIA STATE UNIVERSITY
C                                 UNIVERSITY PARK, PA.  16802
C
      IMPLICIT INTEGER (A-Z)
C
      DIMENSION ORD(N),POPLST(2,20)
      INTEGER X,XX,Z,ZZ,Y
C
C     TO SORT DIFFERENT INPUT TYPES, CHANGE THE FOLLOWING
C     SPECIFICATION STATEMENTS; FOR EXAMPLE, FOR FORTRAN CHARACTER
C     USE THE FOLLOWING:  CHARACTER *(*) A(N)
C
      INTEGER A(N)
C
      NDEEP=0
      U1=N
      L1=1
      DO 1  I=1,N
 1        ORD(I)=I
 2         IF (U1.LE.L1) RETURN
C
 3          L=L1
      U=U1
C
C PART
C
 4     P=L
      Q=U
C     FOR CHARACTER SORTS, THE FOLLOWING 3 STATEMENTS WOULD BECOME
C     X = ORD(P)
C     Z = ORD(Q)
C     IF (A(X) .LE. A(Z)) GO TO 2
C
C     WHERE "CLE" IS A LOGICAL FUNCTION WHICH RETURNS "TRUE" IF THE
C     FIRST ARGUMENT IS LESS THAN OR EQUAL TO THE SECOND, BASED ON "LEN"
C     CHARACTERS.
C
      X=A(ORD(P))
      Z=A(ORD(Q))
      IF (X.LE.Z) GO TO 5
      Y=X
      X=Z
      Z=Y
      YP=ORD(P)
      ORD(P)=ORD(Q)
      ORD(Q)=YP
 5     IF (U-L.LE.1) GO TO 15
      XX=X
      IX=P
      ZZ=Z
      IZ=Q
C
C LEFT
C
 6     P=P+1
      IF (P.GE.Q) GO TO 7
      X=A(ORD(P))
      IF (X.GE.XX) GO TO 8
      GO TO 6
 7     P=Q-1
      GO TO 13
C
C RIGHT
C
 8     Q=Q-1
      IF (Q.LE.P) GO TO 9
      Z=A(ORD(Q))
      IF (Z.LE.ZZ) GO TO 10
      GO TO 8
 9     Q=P
      P=P-1
      Z=X
      X=A(ORD(P))
C
C DIST
C
 10    IF (X.LE.Z) GO TO 11
      Y=X
      X=Z
      Z=Y
      IP=ORD(P)
      ORD(P)=ORD(Q)
      ORD(Q)=IP
 11    IF (X.LE.XX) GO TO 12
      XX=X
      IX=P
 12    IF (Z.GE.ZZ) GO TO 6
      ZZ=Z
      IZ=Q
      GO TO 6
C
C OUT
C
 13    CONTINUE
      IF (.NOT.(P.NE.IX.AND.X.NE.XX)) GO TO 14
      IP=ORD(P)
      ORD(P)=ORD(IX)
      ORD(IX)=IP
 14    CONTINUE
      IF (.NOT.(Q.NE.IZ.AND.Z.NE.ZZ)) GO TO 15
      IQ=ORD(Q)
      ORD(Q)=ORD(IZ)
      ORD(IZ)=IQ
 15    CONTINUE
      IF (U-Q.LE.P-L) GO TO 16
      L1=L
      U1=P-1
      L=Q+1
      GO TO 17
 16    U1=U
      L1=Q+1
      U=P-1
 17    CONTINUE
      IF (U1.LE.L1) GO TO 18
C
C START RECURSIVE CALL
C
      NDEEP=NDEEP+1
      POPLST(1,NDEEP)=U
      POPLST(2,NDEEP)=L
      GO TO 3
 18    IF (U.GT.L) GO TO 4
C
C POP BACK UP IN THE RECURSION LIST
C
      IF (NDEEP.EQ.0) GO TO 2
      U=POPLST(1,NDEEP)
      L=POPLST(2,NDEEP)
      NDEEP=NDEEP-1
      GO TO 18
C
C END SORT
C END QSORT
C
      return

      END

*************************************************************************************




************************************************************************************************

      subroutine median(x,n,xmed)

c x real Feld mit n Eintraegen
c n Anzahl der Eintraege in Feld x
c xmed zu findender Median des Feldes x


c Find the median of X(1),... ,X(N),using as much of the quicksort
c algorithm as is needed to isolate it.
c N.B. On exit, the array X is partially ordered.

c Local variables

      real temp,xhi,xlo,xmax,xmin
      logical odd
      integer hi,lo,nby2,nby2p1,mid,i,j,k
      real x(n)

      nby2=INT(FLOAT(n)/2.)
      nby2p1=nby2+1
      odd=.true.

c     HI & LO are position limits encompassing the median.

      if (n.eq.2*nby2) odd=.false.
      lo=1
      hi=n
      if (n.lt.3) then
        if (n.lt.1) then
          xmed=0.0
          return
        endif
        xmed=x(1)
        if (n.eq.1) return
        xmed=0.5*(xmed+x(2))
        return
      endif

c     Find median of 1st, middle & last values.

 10      mid=(lo+hi)/2
      xmed=x(mid)
      xlo=x(lo)
      xhi=x(hi)

c Swap xhi & xlo

      if (xhi.lt.xlo) then
        temp=xhi
        xhi=xlo
        xlo=temp
      endif
      if (xmed.gt.xhi) then
        xmed=xhi
      else if (xmed.lt.xlo) then
        xmed=xlo
      endif

! The basic quicksort algorithm to move all values.le.the sort key (XMED)
! to the left-hand end, and all higher values to the other end.

      i=lo
      j=hi
 50   do k=1,n
        if (x(i).ge.xmed) goto 51
        i=i+1
      enddo
 51   continue

      do k=1,n
        if (x(j).le.xmed) goto 52
        j=j-1
      enddo
 52   continue

      if (i.lt.j) then
        temp=x(i)
        x(i)=x(j)
        x(j)=temp
        i=i+1
        j=j-1

c     Decide which half the median is in.

        if (i.le.j) goto 50
      endif

      if (odd.eqv..false.) then
        if (j.eq.nby2 .AND. i.eq.nby2p1) goto 130
        if (j.lt.nby2) lo=i
        if (i.gt.nby2p1) hi=j
        if (i.ne.j) goto 100
        if (i.eq.nby2) lo=nby2
        if (j.eq.nby2p1) hi=nby2p1
      else
        if (j.lt.nby2p1) lo=i
        if (i.gt.nby2p1) hi=j
        if (i.ne.j) goto 100

c Test whether median has been isolated.

        if (i.eq.nby2p1) return
      endif
 100    if (lo.lt.hi-1) goto 10

      if (odd.eqv..false.) then
        xmed=0.5*(x(nby2)+x(nby2p1))
        return
      endif
      temp=x(lo)
      if (temp.gt.x(hi)) then
        x(lo)=x(hi)
        x(hi)=temp
      endif
      xmed=x(nby2p1)
      return

c Special case, N even, J=N/2 & I=J+1, so the median is
c between the two halves of the series.   Find max. of the first
c half & min. of the second half, then average.

 130    xmax=x(1)
      do k=lo,j
        xmax=MAX(xmax,x(k))
      enddo
      xmin=x(n)
      do k=i,hi
        xmin=MIN(xmin,x(k))
      enddo
      xmed=0.5*(xmin+xmax)

      return
      end

************************************************************************************
