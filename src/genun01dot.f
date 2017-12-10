c generate N uniformly distributed  dots on sphera RI at origin 0,0,0
c
c Y.V. 2003
c
c INput:
c RI - radius sphere
c N  - required approximate number of dots
c
c RESULT:
c  U(3,k) - dot's xyz
c  AR(k)  - dot's area
c  TET(k) - dot's TETa angle
c  N      - corrected number of dots
c 
        SUBROUTINE genun01dot(RI,U,AR,TET,N,nMX)
c 
        implicit none
c
        include "charStringSiz.h"
        include "output.h"
        include "kanalUse.h"
        include "statusMessg_mDyn.h"
c
        integer N,nMX
        real*4 RI
        real*4 U(*)
        real*4 AR(*),TET(*)
c local 
        real*8 dri
        real*8 FI,FJ,Z,XY,X,Y
        real*8 PI,twoPI
        real*8 dNHOR,dNVERT
        real*8 aat,aaf,RI2
        real*8 area
        real*8 dtet,dtet2,sdtet,cdtet
        real*8 OPT_nhorkk
c 
        integer NEQUAT,NVERT,I,NHOR,J,NU
        integer Nlarg
        integer NHORmin
        integer n3
        character*18 char18
        logical CONTROL
        integer kanalp
c 
        data NHORmin/6/
c 
        kanalp = kanalRunOut
        CONTROL = .false.
c
        char18 = '                '
        OPT_nhorkk=1.3d0
        Nlarg=N*OPT_nhorkk
        OPT_nhorkk=0.5d0*OPT_nhorkk
        PI = dacos(-1.0d0)     ! 3.1415927d0
cx      PI = 3.1415927d0    ! 3.1415927d0
        twoPI=2.0d0*PI
        dri = dble(RI)
        RI2 = dRI**2  
        area = 0.0d0
        NEQUAT = DSQRT(dfloat(N ) * PI)
        NVERT = 0.5d0 * NEQUAT
C!:NVERT is even = 2*n  
        I=nvert/2
        nvert = I + I
        IF (NVERT .LT. 2) NVERT = 2
c
        if(CONTROL)then
        write(kanalp,*)
     &  'GENUN01: RI:',RI,' Ndot:',N,' Nvert:',nvert
        endif
c 
        NU = 0
        dNVERT=dfloat(NVERT)
        dtet = PI/dNVERT
        dtet2= dtet*0.5d0
        sdtet =2.0d0*DSIN(dtet2)
        cdtet = DCOS(dtet2)
c 
        DO 100 I = 0,NVERT
 
           FI = dtet*I
           Z = dRI*DCOS(FI)
           XY = DSIN(FI)

           if(I.eq.0.or.I.eq.NVERT)then
           aat=1.0d0 - cdtet
           else
           aat=XY*sdtet
           end if
 
           NHOR = NEQUAT * XY
C!: even NHOR
cx           NHOR=int(OPT_nhorkk*NHOR + 0.5d0)*2  !horosho dlja bol'shoi dotden ~ 4.
           NHOR=2*int(0.5*NHOR + 0.5)   ! luche dlja malih dotden ~1.
c
           IF (NHOR .LT. NHORmin) NHOR = NHORmin
           if(I.eq.0.or.I.eq.NVERT) NHOR = 1
 
        if(CONTROL)then
        write(kanalp,*)
     &  'GENUN01: iVert:',i,' Nhor:',NHOR,' XY:',XY,' FI:',FI
        endif
 
           dNHOR=dfloat(NHOR)
           aaf=twoPI/dNHOR
           XY = XY*dRI
 
           DO 50 J = 0,NHOR-1
 
              FJ = aaf*J
              X = DCOS(FJ) * XY
              Y = DSIN(FJ) * XY
 
         IF (NU .GE. Nlarg)then
         write(kanalp,*)
     &   'GENUN01: probeSphe is not completed:Nlarg,NU:', Nlarg,NU
c
          write(kanalPStat,*)mError,
     &    ' GENUN01: dot density is too large/ or sphera too large',
     &    ' GENUN01: probeSphe is not completed:Nlarg,NU:', Nlarg,NU
           stop
c
         GO TO 150
         end if
 
              NU = NU + 1
c
         if(nu .gt. nMX)then
         write(kanalp,*)
     &   'genun01:sasDataSiz.h:ERROR! param ndotAtMAX: is low ',nu,nMX
         write(kanalp,*)'genun01:NdotMX,Rat:',N,RI
c
          write(kanalPStat,*)mError,
     &    ' GENUN01: sasDataSiz.h: param ndotAtMAX: is low ',nu,nMX               
c
           stop
         end if
c
              n3 = 3*nu
              U(n3-2) = X 
              U(n3-1) = Y
              U(n3)   = Z 
              AR(NU) = aaf*aat*RI2
              TET(NU) = FI
              area=area + AR(NU)
 
50      CONTINUE
100     CONTINUE
150     CONTINUE
 
        N = NU
        aaf = 4.0d0*PI*RI2/area
c
        area = 0.0
        do i=1,N
        AR(i) = aaf*AR(i)
        area = area + AR(i)
        end do
c 
        if(CONTROL)then
        write(kanalp,'(a16,i5,a12,f16.14)')
     &  'GENUN01(fin):N:',NU,' area norm=',aaf
        write(kanalp,*)'dotN   x   y   z area TETA'
        do i=1,N
        n3=3*i 
        write(kanalp,7001)
     &  'dot   ',i,char18,U(n3-2),U(n3-1),U(n3),AR(i),TET(i)
        enddo
        write(kanalp,*)' area:',area,'areaEx: ',(4.0d0*PI*RI2/area)
        write(kanalp,*)'genun01 : finish !:'
        endif
c
7001    format(a6,i6,a18,3f8.3,f8.3,1x,f8.3) 
c 
        RETURN
        END
