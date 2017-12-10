c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c*                                                                           *
c*  allAtVDWEForce                                                           *
c* Defines  the VDW COUL1(short) energy and Forces                           *
c*                                                                           *
c*     Yury Vorobjev 2002                                                    *
c*                   2004 -  hBHxY128force                                   *
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c allAtNonBondEForce : VDW and COULOMBICR1 - SoftCore
c
        subroutine allAtVDWEForceR1SC(atomXYZ,atomQ,
     &           natom,nmoveatom,moveAtomList,
     &           resNumb,atomBlockName,
     &           nbpairListV,startnbPairLV,nnbPairLV,
     &           pair14List,startPairL14,nPairL14,
     &           nVDWtype,atomVDWtype,atomVDW12ab,aSoftCore,
     &           rcutV,rcutC,engVDW,vdwForce,engCOULR1,coulForceR1,
     &           hBHxYeng128, hBHxY128force)
c
c SHORTrange rij<R1 energy/forces:
c
cx	implicit none
        include 'xyzPDBsize.h'
c
        include "output.h"
        include "kanalUse.h"
        include "statusMessg_mDyn.h"
c
c Used for CoulEng Models
        include 'coulEnPar.h'
        include 'molGenPar.h'
c Used for hB128ang E/force 
        include 'hbond128.h'
c SAS for dielModLaZ
        include 'dataSASdr.h'
        include 'solvGBorn.h'
c PWat:
	include 'enForce_PW.h'
	include 'optionPar.h'
c
        real atomXYZ(*) 
        real atomQ(*)
        integer natom
        integer nmoveatom,moveAtomList(*)
        integer resNumb(*)
        character*4 atomBlockName(*)
        integer nbpairListV(*)
        integer startnbPairLV(*),nnbPairLV(*)
        integer pair14List(*)
        integer startPairL14(*),nPairL14(*)
        integer nVDWtype,atomVDWtype(*)
        real atomVDW12ab(*)
        real rcutV,rcutC,aSoftCore
        real engVDW,vdwForce(*)
        real engCOULR1,coulForceR1(*)
        real hBHxYeng128, hBHxY128force(*)
c
c local
        real cReps
        integer cVar
        real DCij,dDCij,rc1,rc2
        real epsSOL
        real sckal
        real sc14VDW,sc14COUL
        real evdw,fi(3),fci(3)
        real fs(3),qi,qj,xi(3)
        integer i,i1,i2,i3,ia,iam
        integer ih,ih3,ih4
        integer j,ja,j3,k
        integer t1,t2,t12,p4
        integer rcutR12
        real rij(3),dij2,dij1
        real A12,B12,ecoul
        real rm12,rSC12
        real evmax
        logical OPT_SHORTVDW
	integer kanalp
        logical OPT_removePairInNABase 
        logical remove14        
        logical OPT_do14    
c hB128
        real hbAngPar(6)
        real eHbAng,f1(3),f2(3),f3(3)
c
        logical CONTROL,CONTROL0
c
        epsSOL = epsSOL_OPT
c scaling coeff for 14 pair - amber ff94 standart
        sc14VDW = 0.50   ! amber=0.5
        sc14COUL = 0.25  ! amber=0.8        
c take 1/2 symmetric 14 lis:
        sc14VDW = sc14VDW*0.5
        sc14COUL = sc14COUL*0.5
c
        OPT_removePairInNABase = .true.  !remove i-j pair 
c                                        !inside Block-Unit for VDW/COUL     
        sckal = 332.0
        evmax = 0.00    ! bad VDW contact > evmax (kcal)kcal will be printed
        OPT_SHORTVDW = .false.  ! .true. print SHORT VDW CONTACTs
        OPT_do14 = .true. 
c
cx        kanalp = 6
        kanalp = kanalRunOut
c
        CONTROL0 = .false.  !print
        CONTROL = .false.   !print
c
        if(CONTROL0)then
        write(kanalp,*)'* * * * * * * * * * * * * * * *' 
        write(kanalp,*)'Update VDW+COUL r<R1 en/force:'
        write(kanalp,*)'* * * * * * * * * * * * * * * *' 
        end if
c initialize
        cVar = coulVar_OPT
        cReps = cReps_OPT
        rc1 = 1.0/rcutC
        rc2 = rc1**2
        rcutR12 = rcutV**2
        engVDW = 0.0
        engCOULR1 = 0.0        
        hBHxYeng128 = 0.0
c init PWat:
        do k=1,3
	engVDWR1_PWat(k)=0.0
	engCOULR1_PWat(k)=0.0
	hBHxYeng128_PWat(k)=0.0
	end do!k
c
	do i = 1,3*natom
        vdwForce(i) = 0.0
        coulForceR1(i) = 0.0
        hBHxY128force(i) = 0.0
        end do
c
	if(CONTROL0)then
        write(kanalp,*)'allAtNBondEForce: nmoveatom:',nmoveatom
        write(kanalp,*)'moveAtomList:'
        write(kanalp,*)(moveAtomList(i),i=1,nmoveatom)
        end if !C0
c
c loop over movingAtomList
        do iam = 1,nmoveatom
        ia = moveAtomList(iam)
        t1 = atomVDWtype(ia)
        qi = atomQ(ia)
        i3 = 3*ia-3
        xi(1) = atomXYZ(i3+1)
        xi(2) = atomXYZ(i3+2)
        xi(3) = atomXYZ(i3+3)
c
c loop over shortR1 neighbours        
	do j = startnbPairLV(ia),startnbPairLV(ia)+nnbPairLV(ia)-1
        ja = nbpairListV(j)
c atoms ia,ja interacts
        j3 = 3*ja-3
        rij(1) = atomXYZ(j3+1) - xi(1)
        rij(2) = atomXYZ(j3+2) - xi(2)
        rij(3) = atomXYZ(j3+3) - xi(3)
        dij2 = rij(1)**2 + rij(2)**2 + rij(3)**2 
c
        if( dij2 .lt. rcutR12)then  
         dij1 = sqrt(dij2)
c
c VDW e/force
c calculate energy and force
        t2 = atomVDWtype(ja) 
        qj = atomQ(ja)
c get pointer to the AB table
        call vdw12TablePos(nVDWtype,t1,t2,t12)
        p4 = 4*t12
        A12 = atomVDW12ab(p4-3)
        B12 = atomVDW12ab(p4-2)
c
        rm12 = atomVDW12ab(p4-1)
        rSC12 = aSoftCore*rm12
c
	call vdwenforceijSC(dij2,dij1,rij,A12,B12,rSC12,evdw,fi)
c collect
        if(CONTROL)then
        if( dij1 .eq. 0.0 )then
        write(kanalp,*)
     &  'allAtNBondEForceBS:ERROR!: ia,ja,dij,evdw:',
     &   ia,ja,dij1,evdw
         else
         write(kanalp,*)
     &  'allAtNBondEForceBS: ia,ja,dij,evdw:',
     &   ia,ja,dij1,evdw
        end if
        end if !C0
c distribute PWat:
        if(OPT_SolvateExWat)
     &  call engPWatDistribute(atomPWatTag(ia),atomPWatTag(ja),
     &  	evdw,engVDWR1_PWat)
c
        engVDW = engVDW + evdw
        vdwForce(i3+1) =  vdwForce(i3+1) + fi(1)
        vdwForce(i3+2) =  vdwForce(i3+2) + fi(2)
        vdwForce(i3+3) =  vdwForce(i3+3) + fi(3)
c 
        vdwForce(j3+1) =  vdwForce(j3+1) - fi(1)
        vdwForce(j3+2) =  vdwForce(j3+2) - fi(2)
        vdwForce(j3+3) =  vdwForce(j3+3) - fi(3)
c
	if(OPT_SHORTVDW)then      ! print SHORT VDW contacts
        if( evdw .gt. evmax )then
         write(kanalp, '(a20,2i6,f7.3,f9.2,3i5,2f11.2)')
     & 'SHORTVDW!!:vdwList',ia,ja,dij1, evdw,t1,t2,t12,A12,B12
         end if
        end if! OPT_SHORTVDW
c
c COUL shortRange
        if( qi .ne. 0.0 .and. qj .ne. 0.0 )then
c
        if(cVar .eq. 0)then
        call coulenforceij0SC(dij2,dij1,rij,qi,qj,rSC12,ecoul,fci)
        end if
c
        if(cVar .eq. 1)then 
        call coulenforceij01SC
     &               (rc1,rc2,dij2,dij1,rij,qi,qj,rSC12,ecoul,fci) 
        end if
c
        if(cVar .eq. 2)then
        call  coulenforceij02SC
     &         (cReps,rc1,rc2,dij2,dij1,rij,qi,qj,rSC12,ecoul,fci) 
        end if
c
        if(cVar .eq. 3)then
c
cx       call getDijConst04(atomDconst,epsMol_OPT,epsSol_OPT,
cx     &                        molRGyr_c,AdcPar_OPT,ndcPar_OPT,
cx     &                        ia,ja,dij1,DCij,dDCij)
c
        call getDijConst05(parDrLZmodel,
     &               atDistToSAS,atSASexp,ia,ja,dij1,DCij,dDCij)
c
        call coulenforceij03SC
     &   (DCij,dDCij,rc1,rc2,dij2,dij1,rij,qi,qj,rSC12,ecoul,fci) 
c
        end if
c
        if(cVar .eq. 4)then  ! calculation pairElect via GBorn approximation
c data atBornRad(ia),atBornRadDr() : solvGBorn.h
        call electPforceijGBSC(rc1,rc2,dij2,dij1,
     &              rij,qi,qj,rSC12,
     &              atBornRad(ia),atBornRad(ja),ecoul,fci)   
c
        end if ! cVar = 4
c
c collect
c
c distribute PWat:
        if(OPT_SolvateExWat)
     &  call engPWatDistribute(atomPWatTag(ia),atomPWatTag(ja),
     &          (ecoul*sckal),engCOULR1_PWat)
c
        engCOULR1 = engCOULR1  + ecoul
        coulForceR1(i3+1) = coulForceR1(i3+1) + fci(1)
        coulForceR1(i3+2) = coulForceR1(i3+2) + fci(2)
        coulForceR1(i3+3) = coulForceR1(i3+3) + fci(3)
c
        coulForceR1(j3+1) = coulForceR1(j3+1) - fci(1)
        coulForceR1(j3+2) = coulForceR1(j3+2) - fci(2)
        coulForceR1(j3+3) = coulForceR1(j3+3) - fci(3)
c
cx	if(CONTROL0)then
cx        if( evdw .gt. evmax )then
cx         write(kanalp, '(a20,2i6,f7.3,f9.5,3f9.5)')
cx     & 'SHORTCOL!!:vdwList',ia,ja,dij1, ecoul,qi,qj,evdw            
cx         end if
cx        end if! CONTROL
c
        end if ! qi*qj ne 0                
	end if !dij < rcut2
c
	end do !shortR1 j neigbours
c
c VDW + COUL for 14V list
c
        if(OPT_do14)then
c
        if(nPairL14(ia) .ge. 1 )then 
c
        do j = startPairL14(ia),startPairL14(ia)+nPairL14(ia)-1
        ja = pair14List(j)
c
        remove14 = .false.
c
        if(OPT_removePairInNABase)then
        if(atomBlockName(ia)(1:1) .eq. 'B' )then
        if(resNumb(ia) .eq. resNumb(ja) .and.
     &    atomBlockName(ja)(1:1) .eq. 'B' )then
          remove14 = .true.    ! skip 14 pair ia-ja
         end if ! atomBlockName
        end if ! atomBlockName(ia)(1:1) .eq. 'B'
        end if !OPT_removePairInNABase
c   
        if((.not. remove14) .and. (ja .gt. ia)) then   !include
c
c        if(CONTROL)then
c         write(kanalp,*)'AllAtNBondEForce: j,pair14List(j):',j,ja    
c         end if    
c 14atoms ia,ja interacts
        j3 = 3*ja-3
        rij(1) = atomXYZ(j3+1) - xi(1)
        rij(2) = atomXYZ(j3+2) - xi(2)
        rij(3) = atomXYZ(j3+3) - xi(3)
        dij2 = rij(1)**2 + rij(2)**2 + rij(3)**2 

        if( dij2 .lt. rcutR12)then
c calculate energy and force
        t2 = atomVDWtype(ja) 
        qj = atomQ(ja)
c get pointer to the AB table
        call vdw12TablePos(nVDWtype,t1,t2,t12)
        p4 = 4*t12
        A12 = atomVDW12ab(p4-3)
        B12 = atomVDW12ab(p4-2)
        dij1 = sqrt(dij2)
c
        rm12 = atomVDW12ab(p4-1)
        rSC12 = aSoftCore*rm12

	call vdwenforceijSC(dij2,dij1,rij,A12,B12,rSC12,evdw,fi)
c
        if(CONTROL)then
        if(dij1 .eq. 0.0)then
        write(kanalp,*)
     &  'allAtNBondEForceBS:14ERR: ia,ja,dij,evdw:',ia,ja,dij1,evdw
         else
         write(kanalp,*)
     &  'allAtNBondEForceBS: ia,ja,dij,evdw:',ia,ja,dij1,evdw
        end if !
        end if !C0
c collect
        engVDW = engVDW + evdw*sc14VDW 
c
        vdwForce(i3+1) =  vdwForce(i3+1) + fi(1)*sc14VDW
        vdwForce(i3+2) =  vdwForce(i3+2) + fi(2)*sc14VDW
        vdwForce(i3+3) =  vdwForce(i3+3) + fi(3)*sc14VDW
c 
        vdwForce(j3+1) =  vdwForce(j3+1) - fi(1)*sc14VDW
        vdwForce(j3+2) =  vdwForce(j3+2) - fi(2)*sc14VDW
        vdwForce(j3+3) =  vdwForce(j3+3) - fi(3)*sc14VDW
c
	if(OPT_SHORTVDW)then
        if( evdw .gt. evmax )then
         write(kanalp, '(a20,2i6,2f8.3,3i5,2f11.2)')
     & 'SHORTVDW!!: 14List',ia,ja,dij1, evdw,t1,t2,t12,A12,B12
         end if
        end if! OPT_SHORTVDW
c
c COUL shortRange
        if( qi .ne. 0.0 .and. qj .ne. 0.0 )then
c
        if(cVar .eq. 0)then
        call coulenforceij0SC(dij2,dij1,rij,qi,qj,rSC12,ecoul,fci)
        end if
c
        if(cVar .eq. 1)then 
        call coulenforceij01SC
     &               (rc1,rc2,dij2,dij1,rij,qi,qj,rSC12,ecoul,fci) 
        end if
c
        if(cVar .eq. 2)then
        call  coulenforceij02SC
     &          (cReps,rc1,rc2,dij2,dij1,rij,qi,qj,rSC12,ecoul,fci) 
        end if
c
        if(cVar .eq. 3)then
c
c        call getDijConst04(atomDconst,epsMol_OPT,epsSol_OPT,
c     &                        molRGyr_c,AdcPar_OPT,ndcPar_OPT,
c     &                        ia,ja,dij1,DCij,dDCij)
c

        call getDijConst05(parDrLZmodel,
     &               atDistToSAS,atSASexp,ia,ja,dij1,DCij,dDCij)
c
        call coulenforceij03SC
     &      (DCij,dDCij,rc1,rc2,dij2,dij1,rij,qi,qj,rSC12,ecoul,fci) 
c
        end if ! cVAr=3
c
         if(cVar .eq. 4)then  ! calculation pairElect via GBorn approximation
c data atBornRad(ia),atBornRadDr() : solvGBorn.h
c
        call electPforceijGBSC(rc1,rc2,dij2,dij1,
     &              rij,qi,qj,rSC12,
     &              atBornRad(ia),atBornRad(ja),ecoul,fci)
 
        end if ! cVar = 4
c
c collect
        engCOULR1 = engCOULR1  + ecoul*sc14COUL 
c
        coulForceR1(i3+1) = coulForceR1(i3+1) + fci(1)*sc14COUL
        coulForceR1(i3+2) = coulForceR1(i3+2) + fci(2)*sc14COUL
        coulForceR1(i3+3) = coulForceR1(i3+3) + fci(3)*sc14COUL
c
        coulForceR1(j3+1) = coulForceR1(j3+1) - fci(1)*sc14COUL
        coulForceR1(j3+2) = coulForceR1(j3+2) - fci(2)*sc14COUL
        coulForceR1(j3+3) = coulForceR1(j3+3) - fci(3)*sc14COUL
c
cx	if(CONTROL0)then
cx        if( evdw .gt. evmax )then
cx         write(kanalp, '(a20,2i6,f7.3,f9.5,3f9.5)')
cx     & 'SHORTCOL!!: 14List',ia,ja,dij1, ecoul,qi,qj,evdw       
cx         end if
cx        end if! CONTROL
c
        end if ! qi*qj ne 0                
	end if !dij < rcut2
c
        end if ! ja>ia
	end do !j 14neigbours
        end if ! n14 > 0
        end if ! OPT_do14
c
        end do !iam
c
c scale COUL to kcal/mol
        do i = 1,3*natom
        coulForceR1(i) = coulForceR1(i)*sckal
        end do
        engCOULR1 = engCOULR1*sckal
c
        if(CONTROL)then
        write(kanalp,*)'vwdEngy : ', engVDW
        write(kanalp,*)'vwdForce :ia: f1,f2,f2 fs1,fs2,fs2 '
        fs(1)= 0.0
        fs(2) = 0.0
        fs(3) = 0.0
        do ia = 1,natom
        i3=3*ia-3
        fs(1) = fs(1) + vdwForce(i3+1)
        fs(2) = fs(2) + vdwForce(i3+2)
        fs(3) = fs(3) + vdwForce(i3+3)
        write(kanalp,'(i6,2x,3f8.3,2x,3f8.3)')ia,
     &  vdwForce(i3+1),
     &  vdwForce(i3+2), vdwForce(i3+3), fs(1),fs(2),fs(3)
        end do!i
c
        write(kanalp,*)'coulEngyR1: ',  engCOULR1
        write(kanalp,*)'coulForce :ia: f1,f2,f2 fs1,fs2,fs2 '
        fs(1)= 0.0
        fs(2) = 0.0
        fs(3) = 0.0
        do ia = 1,natom
        i3=3*ia-3
        fs(1) = fs(1) + coulForceR1(i3+1)
        fs(2) = fs(2) + coulForceR1(i3+2)
        fs(3) = fs(3) + coulForceR1(i3+3)
        write(kanalp,'(i6,2x,3f8.3,2x,3f8.3)')ia,
     &  coulForceR1(i3+1),
     &  coulForceR1(3+2), coulForceR1(i3+3), fs(1),fs(2),fs(3)
        end do!i
        end if!control
c 
c hB128ang
c
c initialize
        do i = 1,nHBbondHxYMAX
        hBHxYeng128List(i) = 0.0
	end do !i
c
	do i = 1,3*natom
        hBHxY128force(i) = 0.0
        end do
        hBHxYeng128 = 0.0 
c
        if(nHBbondHxY .gt. 0)then
c
cd        CONTROL0 = .true.
        if(CONTROL0)then
        write(kanalp,*)'allAtNBondEForce: iHb128, eHb, X-H..Y (xyz):'
        end if!C
c
        do ih=1,nHBbondHxY
        ih3=ih*3-2
c get atoms Y,H,X
        i1=3*hB128List(ih3)-3
        i2=3*hB128List(ih3+1)-3
        i3=3*hB128List(ih3+2)-3
c
        ih4=ih*4-4
        do k=1,4
        hbAngPar(k)=hB128PotParList(ih4+k)
        end do !k
c
        hbAngPar(5) = hB128AgLpar(1)
        hbAngPar(6) = hB128AgLpar(2)
c
        call getHb128angEforce
     &       (atomXYZ(i1+1),atomXYZ(i2+1),atomXYZ(i3+1),
     &                    hbAngPar,eHbAng,f1,f2,f3)
c
        eHbAng = eHbAng*hB128WfList(ih)
	hBHxYeng128List(ih) = eHbAng
c
         if(CONTROL0)then
        write(kanalp,'(i4,1x,f7.3,3(i3,1x,3f7.2,1x))')
     &  ih, eHbAng,
     &  hB128List(ih3),(atomXYZ(i1+k),k=1,3),
     &  hB128List(ih3+1),(atomXYZ(i2+k),k=1,3), 
     &  hB128List(ih3+2),(atomXYZ(i3+k),k=1,3) 
         end if!
c collect
c PWat:define hBHxYeng128_PWat
         if(OPT_SolvateExWat)then
         call engPWatDistribute(atomPWatTag(hB128List(ih3)),
     &   atomPWatTag(hB128List(ih3+2)),eHbAng,hBHxYeng128_PWat)
         end if!OPT_SolvateExWat
c	
	hBHxYeng128 = hBHxYeng128 + eHbAng
        do k=1,3
        hBHxY128force(i1+k) = hBHxY128force(i1+k)+f1(k)*hB128WfList(ih)
        hBHxY128force(i2+k) = hBHxY128force(i2+k)+f2(k)*hB128WfList(ih)
        hBHxY128force(i3+k) = hBHxY128force(i3+k)+f3(k)*hB128WfList(ih)
        end do !k
c
        end do !ih
c
        end if ! nHBbondHxY .gt. 0
c
cd        stop
c
	return
	end
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c*  allAtVDWEForce                                                           *
c* Defines  the  COUL2(longRange) energy and Forces                          *
c                R1 < rij < R2                                               *
c*                                                                           *
c*     Yury Vorobjev 2002                                                    *
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
        subroutine allAtVDWEForceR2(atomXYZ,atomQ,
     &           natom,nmoveatom,moveAtomList,
     &           nbpairListC,startnbPairLC,nnbPairLC,
     &           rcutR1,rcutR2,engCOULR2,coulForceR2)
c
c LongRamge -   RCUT1 < rij < RCUT2

cx	implicit none
        include 'xyzPDBsize.h'
        include 'coulEnPar.h'
        include 'molGenPar.h'
c SAS for dielModLaZ
        include 'dataSASdr.h'
        include 'solvGBorn.h'
c PWat:
        include 'enForce_PW.h'
        include 'optionPar.h'
c
        real atomXYZ(*) 
        real atomQ(*)
        integer natom
        integer nmoveatom,moveAtomList(*)
        integer nbpairListC(*)
        integer startnbPairLC(*),nnbPairLC(*)
        real rcutR1,rcutR2
        real engCOULR2,coulForceR2(*)
c
c local
        integer cVar
        real cReps,epsSOL
        real DCij,dDCij,rc1,rc2,rSC12
        real fci(3)
        real fs(3),qi,qj,xi(3)
        integer i,i3,ia,iam
        integer j,ja,j3
        integer t1,t2,t12,p4
        integer rcutR12,rcutR22
        real rij(3),dij2,dij1
        real ecoul,sckal
	integer kanalp
        logical CONTROL,CONTROL0
c
        sckal = 332.0
c
        kanalp = 6
        CONTROL = .false.
        CONTROL0 = .false.
c
        if(CONTROL)then
        write(kanalp,*)'* * * * * * * * * * * * * * * * *' 
        write(kanalp,*)'Update COULR2 (R1<r<R2) en/force:'
        write(kanalp,*)'* * * * * * * * * * * * * * * * *' 
        end if
c initialize
        epsSOL = epsSOL_OPT 
        rSC12 = 1.0         !
        cReps = cReps_OPT
        cVar = coulVar_OPT 
        rc1 = 1.0/rcutR2
        rc2 = rc1**2
        rcutR22 = rcutR2**2
        rcutR12 = rcutR1**2
        engCOULR2 = 0.0
	do i = 1,3*natom
        coulForceR2(i) = 0.0
        end do
c 
        engCOULR2_PWat(1)=0.0
	engCOULR2_PWat(2)=0.0
	engCOULR2_PWat(3)=0.0
c
c loop over movingAtomList
        do iam = 1,nmoveatom
        ia = moveAtomList(iam)
        qi = atomQ(ia)
        i3 = 3*ia-3
        xi(1) = atomXYZ(i3+1)
        xi(2) = atomXYZ(i3+2)
        xi(3) = atomXYZ(i3+3)

c loop over shortR1 neighbours        
	do j = startnbPairLC(ia),startnbPairLC(ia)+nnbPairLC(ia)-1
        ja = nbpairListC(j)
c atoms ia,ja interacts
        j3 = 3*ja-3
        rij(1) = atomXYZ(j3+1) - xi(1)
        rij(2) = atomXYZ(j3+2) - xi(2)
        rij(3) = atomXYZ(j3+3) - xi(3)
        dij2 = rij(1)**2 + rij(2)**2 + rij(3)**2 

        if( dij2 .ge. rcutR12 .and. 
     &                    dij2 .lt. rcutR22)then
c calculate energy and force
        qj = atomQ(ja)
        dij1 = sqrt(dij2)

        if( qi .ne. 0.0 .and. qj .ne. 0.0 )then
c
        if(cVar .eq. 0)then
        call coulenforceij0(dij2,dij1,rij,qi,qj,ecoul,fci)
        end if
c
        if(cVar .eq. 1)then 
        call coulenforceij01(rc1,rc2,dij2,dij1,rij,qi,qj,ecoul,fci) 
        end if
c
        if(cVar .eq. 2)then
        call 
     &  coulenforceij02(cReps,rc1,rc2,dij2,dij1,rij,qi,qj,ecoul,fci) 
        end if
c
        if(cVar .eq. 3)then
c
cx        call getDijConst04(atomDconst,epsMol_OPT,epsSol_OPT,
cx     &                        molRGyr_c,AdcPar_OPT,ndcPar_OPT,
cx     &                        ia,ja,dij1,DCij,dDCij)
c
c
        call getDijConst05(parDrLZmodel,
     &               atDistToSAS,atSASexp,ia,ja,dij1,DCij,dDCij)
c
        call coulenforceij03
     &         (DCij,dDCij,rc1,rc2,dij2,dij1,rij,qi,qj,ecoul,fci) 
c
        end if !cVar=3
c
        if(cVar .eq. 4)then  ! calculation pairElect via GBorn approximation
c data atBornRad(ia),atBornRadDr() : solvGBorn.h
c
        call electPforceijGBSC(rc1,rc2,dij2,dij1,
     &              rij,qi,qj,rSC12,
     &              atBornRad(ia),atBornRad(ja),ecoul,fci)
 
        end if ! cVar = 4
c
c c distribute PWat:
        if(OPT_SolvateExWat)
     &  call engPWatDistribute(atomPWatTag(ia),atomPWatTag(ja),
     &          (ecoul*sckal),engCOULR2_PWat)
c
c collect
        engCOULR2 = engCOULR2  + ecoul
        coulForceR2(i3+1) = coulForceR2(i3+1) + fci(1)
        coulForceR2(i3+2) = coulForceR2(i3+2) + fci(2)
        coulForceR2(i3+3) = coulForceR2(i3+3) + fci(3)
c
        coulForceR2(j3+1) = coulForceR2(j3+1) - fci(1)
        coulForceR2(j3+2) = coulForceR2(j3+2) - fci(2)
        coulForceR2(j3+3) = coulForceR2(j3+3) - fci(3)
c
        end if ! qi*qj ne 0                
	end if !dij < rcut2
c
	end do !j neigbours
c
        end do !iam
c
c scale COUL to kcal/mol
        do i = 1,3*natom
        coulForceR2(i) = coulForceR2(i)*sckal
        end do
        engCOULR2 = engCOULR2*sckal
c
        if(CONTROL)then
        write(kanalp,*)'coulEngyR2: ',  engCOULR2
        write(kanalp,*)'coulForceR2 :ia: f1,f2,f2 fs1,fs2,fs2 '
        fs(1)= 0.0
        fs(2) = 0.0
        fs(3) = 0.0
        do ia = 1,natom
        i3=3*ia-3
        fs(1) = fs(1) + coulForceR2(i3+1)
        fs(2) = fs(2) + coulForceR2(i3+2)
        fs(3) = fs(3) + coulForceR2(i3+3)
        write(kanalp,'(i6,2x,3f8.3,2x,3f8.3)')ia,
     &  coulForceR2(i3+1),
     &  coulForceR2(3+2), coulForceR2(i3+3), fs(1),fs(2),fs(3)
        end do!i
        end if!control
c
c        write(kanalp,*)'allAtNBondEForceBS: STOP'
c	stop
c
	return
	end
