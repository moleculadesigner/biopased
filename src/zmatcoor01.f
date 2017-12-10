c
c atomXYZ from amberstyle ZMatrix
c 
       subroutine zmcoord01(NAT,
     &    molZmatrx1,molZmatrx2,molZmatrx3,
     &    atXYZ,MNAT,atName)
c
       implicit none
       integer NAT
       character*4 molZmatrx1(*)
       integer molZmatrx2(*)
       real molZmatrx3(*)
       real atXYZ(*)
       integer  MNAT(*)  
       character*4 atName(*)
c
       include "charStringSiz.h"
       include "output.h"
       include "kanalUse.h"
       include "statusMessg_mDyn.h"
c local
       integer natzmcMAX
       parameter (natzmcMAX=1000 ) ! variable
       real DC(3,natzmcMAX),ETA(natzmcMAX)
       real C(3,natzmcMAX),R(natzmcMAX),ZETA(natzmcMAX)
       integer NI(natzmcMAX+1),NF(natzmcMAX),NBACK(natzmcMAX+1)
       real DEL(3),C0(3) 
c*
       integer M,NIT,NBOND,NFT
       integer NBACK1,NBACK2,NNIT
       integer i,i1,i2,i4,i3,i5
       integer k,j,j3,j5,ii,jj 
       integer ires
       real FAC,TEST,RS,ALOC,BLOC,CLOC
       real ZETA1,SINZ,COSZ,TANZ
c
       integer kanalp     
       logical CONTROL
c
       kanalp =  kanalRunOut  
       CONTROL = .true.

       write(kanalp,*)'In zmatcoor: start'   

       if(nat .gt. natzmcMAX-3) then
       write(kanalp, *)'ERROR:zmatcoor !! natzmcMAX < NAT: ',
     &       natzmcMAX, NAT
       stop 
       end if
c z-matrix
       if(CONTROL)then
       write(kanalp,*)'In Zmatrix: '
       do i = 1,NAT
       i3=3*i-3
       i4=4*i-4
       write(kanalp,'(3a4,4i6,3f8.3,1x,f9.4)')
     & (molZmatrx1(i3+k),k=1,3),(molZmatrx2(i4+k),k=1,4),
     & (molZmatrx3(i3+k),k=1,3) 
       end do!i
       end if !C
c init
       do i=1,nat*3
       atXYZ(i) = 0.0
       end do
c*
c copy to internal variables
       M=NAT+1
c
c       write(kanalp,*)'nbond, zeta(i),eta(i),r(i)'
c
       do i = 3,M
       J=2*I-5
       i4 = 4*i-8 
       i3 = 3*i-6 
       ni(i) = molZmatrx2(i4+2)
       nf(i) = molZmatrx2(i4+1)
c
c       write(kanalp,*)i, ni(i),nf(i)     
c  
       r(i) = molZmatrx3(i3+1)
       zeta(i) = molZmatrx3(i3+2)
       eta(i) = molZmatrx3(i3+3)
c
c       write(kanalp,*) i, zeta(i),eta(i),r(i)
c
       end do !i
c
        DC(1,1)=1.
        DC(2,1)=0.
        DC(3,1)=0.
        DC(1,2) =0.
        DC(2,2)=0.
        DC(3,2)=1.
c standart orienation in XY plane
        eta(1)= 0.0
        zeta(1)=0.0
        r(1) = 0.0 
        r(2) = 0.0
        eta(2)= 0.0
        zeta(2)=90.0
        eta(3) = 180.0
        zeta(3) = 90.0
        eta(4) = 270.0
c
        NI(2)=natzmcMAX
        NBACK(natzmcMAX)=1
c
        FAC=1.74532952E-2
        do i = 1,M
        eta(i) = eta(i)*FAC
        zeta(i) = zeta(i)*FAC
        end do !i
c
        NIT=NI(3)
        NBACK(NIT)=2
c
        do  I=1,3
        C(I,NIT)=0.0
        end do !i
c
        DO 10 NBOND=3,M
        
c        write(kanalp,*)'DO10:W1: nbond=',nbond
c
        NFT=NF(NBOND)
        NIT=NI(NBOND)
        NBACK(NFT)=NBOND
c
c        write(kanalp,*)'W2:NIT,NFT,NBACK: ',NIT,NFT,NBACK(NFT)
c
        TEST=ABS(ZETA(NBOND)-180.)
        IF(TEST.GT.0.00001)GOTO 35
        NBACK1=NBACK(NIT)
        DEL(1)=R(NBOND)*DC(1,NBACK1)
        DEL(2)=R(NBOND)*DC(2,NBACK1)
        DEL(3)=R(NBOND)*DC(3,NBACK1)
        GOTO 36
c
35      RS=R(NBOND)*SIN(ZETA(NBOND))
        ALOC=RS*COS(ETA(NBOND))
        BLOC=RS*SIN(ETA(NBOND))
        CLOC=-R(NBOND)*COS(ZETA(NBOND))
        NBACK1=NBACK(NIT)
c
37      NNIT=NI(NBACK1)
        NBACK2=NBACK(NNIT)
        ZETA1=ZETA(NBACK1)
c
        TEST=ABS(ZETA1-180.)
        IF(TEST.GT.0.00001)GOTO 38
        NBACK1=NBACK2
        GOTO 37
c
c38      ZETA1=ZETA1*FAC
38     SINZ=SIN(ZETA1)
       COSZ=COS(ZETA1)
       TANZ=-COSZ/SINZ
c
c       write(kanalp,*)'nback1, nback2,zeta(nback1):',
c     &  nback1, nback2,zeta(nback1)
c       write(kanalp,*)'zeta1, sinz,cosz,tanz:',
c     & zeta1, sinz,cosz,tanz
c
       DEL(1)=(TANZ*DC(1,NBACK1)-(DC(1,NBACK2)/SINZ))*ALOC+
     * (DC(3,NBACK1)*DC(2,NBACK2)-DC(2,NBACK1)*DC(3,NBACK2))*
     * BLOC/SINZ+DC(1,NBACK1)*CLOC
       DEL(2)=(TANZ*DC(2,NBACK1)-(DC(2,NBACK2)/SINZ))*ALOC+
     * (DC(1,NBACK1)*DC(3,NBACK2)-DC(3,NBACK1)*DC(1,NBACK2))*
     * BLOC/SINZ+DC(2,NBACK1)*CLOC
       DEL(3)=(TANZ*DC(3,NBACK1)-(DC(3,NBACK2)/SINZ))*ALOC+
     *  (DC(2,NBACK1)*DC(1,NBACK2)-DC(1,NBACK1)*DC(2,NBACK2))*
     *  BLOC/SINZ+DC(3,NBACK1)*CLOC
c
       j3=3*NFT-3
36     DO 6 I=1,3
       DC(I,NBOND)=DEL(I)/R(NBOND)
       C(I,NFT)=DEL(I)+C(I,NIT)
6      atXYZ(j3+I)=C(I,NFT)
c
c      PRINT 1001,NBOND,NIT,NFT,C(1,NFT),C(2,NFT),C(3,NFT)
c end DO10
   10 CONTINUE
c
c       write(kanalp,*)'End of DO10'
c
       DO 61 I=1,3
61     atXYZ(I)=0.0
c
       I2=NAT*5
       DO 55 I=1,I2
55      MNAT(I)=0
c
       DO 44 II=3,M
       J5=5*NI(II)-5
       I5=5*NF(II)-5
c
       MNAT(J5+1)=NI(II)
       MNAT(I5+1)=NF(II)
       DO 46 JJ=2,5
        I1=0
        IF(MNAT(J5+JJ).EQ.0)I1=1
        IF(I1.EQ.1)then
        MNAT(J5+JJ)=NF(II)
        GOTO49
        end if
46     CONTINUE
49     CONTINUE
c
       DO 48 JJ=2,5
        I1=0
        IF(MNAT(I5+JJ).EQ.0)I1=1
        IF(I1.EQ.1)MNAT(I5+JJ)=NI(II)
        IF(I1.EQ.1)GOTO52
        IF(MNAT(I5+JJ).EQ.0)GOTO52
48     CONTINUE
52     CONTINUE
44     CONTINUE
c
c assign atomName
       do i = 1,nat
       i3 = 3*i-3
       atName(i) = molZmatrx1(i3+1)
cx       atFFname(i) = molZmatrx1(i3+2)(1:2)   
       end do !i
c
       if(CONTROL)then
       write(kanalp,*)'ZMatcoorRes:'
       write(kanalp,'(a20)' )'nat  x   y   z  '
       ires = 1
       do i = 1,nat
       i3=3*i-3
c       write(kanalp,'(i6,1x,3f8.3)')i,(atXYZ(i*3-3+j),j=1,3)
       write(kanalp,7071)
     &    'ATOM  ',i, atName(i), molZmatrx1(i3+3), ' ',
     &  ires, (atXYZ(i3+j),j=1,3)
       end do !i
       write(kanalp,'(a20)' )'Connectivity matrix:'
       do i = 1,nat
       i5 = i*5 - 5
       write(kanalp,'(i6,2x,4i5)')(MNAT(i5+k),k=1,5)
       end do !i
       end if
c
1001   FORMAT(2X,'IN RATCOR,NB,NI,NFT,XYZ=',3I5,3F8.3)
7071    format(a6,1x,i4,2x,a4,a4,a1,i4,4x,3f8.3) !orig PDB 
c
       RETURN
       END
