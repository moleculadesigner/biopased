c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c                                                                           *
c  initialize: all FAST forces on atoms                                     * 
c                                                                           *
c Yury Vorobjev    2002-2005                                                * 
c                                                                           *
c SoftCOre                                                                  *
c USED for MD                                                               *
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
	subroutine initAllForce01(fcall,fEngWF,atomXYZ,
     &              eVbondDef,vbdefForce,
     &              eVangDef,vAngdefForce,
     &              eImpDef,impDefForce,
     &              eTorsDef,torsAngForce,
     &              engVDWR1,vdwForceR1,
     &              engCOULR1,coulForceR1,
     &              hBHxYeng128,hBHxY128force,
     &              restr1Eng,restr1AtForce, 
     &              restr1MHWEng,restr1MHWAtForce,
     &              restr1RHWEng,restr1RHWAtForce, 
     &              restrDistA2Eng,restrDistA2Force,
     &              engFastF1,atomForceF1)
c
        include 'xyzPDBsize.h'
        include 'xyzPDBinfo.h'
c
        include "output.h"
        include 'pair1234array.h'
        include 'nbondPairVCS.h'
        include 'vdw12Par.h'
        include 'restrainInfo.h'
        include 'loopInfo.h'
        include 'movingAtom.h'
        include 'optionPar.h'
        include 'shake.h'
        include 'restrain1MHW.h'
        include 'restrain1RHW.h'
        include 'restrDistA2.h'
        include 'hbond128.h'
c
        real fEngWF(*)
        real atomXYZ(*)
        integer fcall    ! sequemtial number of call
        real eVbondDef,vbdefForce(*)
        real eVangDef,vAngdefForce(*)
        real eImpDef,impDefForce(*)
        real eTorsDef,torsAngForce(*)
        real engVDWR1,vdwForceR1(*)
        real engCOULR1,coulForceR1(*)
        real hBHxYeng128,hBHxY128force(*)
        real restr1Eng,restr1AtForce(*)
        real restr1MHWEng,restr1MHWAtForce(*)
        real restr1RHWEng,restr1RHWAtForce(*)
        real restrDistA2Eng,restrDistA2Force(*)
        real engFastF1,atomForceF1(*)
c
       integer i,i3,k
       integer kanalp
       logical CONTROL
       logical OPT_SoftCore  
c
       kanalp = kanalRunOut
       CONTROL = .false.
c       
       OPT_SoftCore = .true. 
c
c       write(kanalp,*)'initAllForce start!: fcall=',fcall
c
c init to zero All forces on Atoms
       if(fcall .eq. 0)then
       eVbondDef = 0.0
       eVangDef = 0.0      
       eImpDef = 0.0
       restr1Eng = 0.0            
       restr1MHWEng = 0.0
       restr1RHWEng = 0.0  
       restrDistA2Eng = 0.0
        engVDWR1 = 0.0
       engCOULR1 = 0.0
       eTorsDef = 0.0
       engFastF1 = 0.0        
       hBHxYeng128 = 0.0
c
       do i = 1,3*natom
       vbdefForce(i) = 0.0
       vAngdefForce(i) = 0.0
       impDefForce(i) = 0.0
       torsAngForce(i) = 0.0
       restr1AtForce(i) = 0.0
       restr1MHWAtForce(i) = 0.0
       restr1RHWAtForce(i) = 0.0
       restrDistA2Force(i) = 0.0
       vdwForceR1(i) = 0.0 
       coulForceR1(i) = 0.0
       atomForceF1(i) = 0.0
       atomForceF1(i) = 0.0
       hBHxY128force(i) = 0.0 
       end do !i
       end if !fcall
c
c all GeoDef forces are calculated at each step
c
 	call allAtVBondEForce(atomXYZ,
     &           natom,bond12List,nbond12,nbond12noH,iShake,
     &           bond12ParL,eVbondDef,vbdefForce )
c
	call allAtVangEForce(atomXYZ,
     &           natom,trip123List,nTrip123,ang123ParL,
     &           eVangDef,vAngdefForce )
c
	call allAtImpTEForce(atomXYZ,
     &           natom,quarImp1234L,nImp1234,impAng1234ParL,
     &           eImpDef,impDefForce )
c
c torsionEnForces
c
        call allAtTorsEForce(atomXYZ,
     &           natom,quar1234List,nQuar1234,
     &           quar1234ParL,quar1234nPar,
     &           eTorsDef,torsAngForce )
c
c restrain1atomForce
c
          call refAtPosRestr1EF(atomXYZ,
     &           natom,nRestr1atom,restr1atomList,
     &           refAtomPos,restr1AtConstList,restr1Eng,restr1AtForce)
c
c refMolPosRestr1HWEF - molecule CM restrain
c
         call refMolPosRestr1HWEF(atomXYZ,natom,
     &           nRestr1MHWatom,restr1MHWatList,atMs1MHWatList,
     &           ref1MHWPos,restr1MHWConst,sizeRestr1MHW,
     &           restr1MHWEng,restr1MHWAtForce)
c 
c USED for exWatShell:
         call refResPosRestr1RHWEF(atomXYZ,natom,
     &           atomMass,startAtInRes,stopAtInRes,
     &           nRestr1RHWAtom,restr1RHWatomList,
     &           nRestr1ResidHW,restr1ResidHWList,
     &           restr1RefResidPosHW,restr1ResidMass,
     &           restr1RHWConst,sizeRestr1RHW,
     &           restr1RHWEng,restr1RHWAtForce)
c
c       write(kanalp,*)'initAllForce : finish refResPosRestr1RHWEF:'
c
c
c restrDistA2EngForces - a1-a2 distant restrain
c
        if(CONTROL)then
        write(kanalp,*)'initAllForce01:nRestrDistA2:',nRestrDistA2
c
         write(kanalp,*)'initAllForce01:restrDistA2List:',
     &  (restrDistA2List(i),i=1,2*nRestrDistA2)
         end if !CL
c
          call allAtDistA2EnForce(atomXYZ,
     &           natom,nRestrDistA2,restrDistA2List,
     &           restrDistA2DistHK,restrDistA2Eng,restrDistA2Force)
c
        call  allAtVDWEForceR1SC(atomXYZ,atomQ,
     &        natom,nmoveatom,moveAtomList,
     &        resNumb,atomBlockName, 
     &        nbpairListV,startnbPairLV,nnbPairLV,
     &        pair14VList,startPairL14V,nPairL14V, 
     &        nVDWtype,atomVDWtype,atomVDW12ab,aSoftCore,
     &        rcutV,rcutC,engVDWR1,vdwForceR1,engCOULR1,coulForceR1,
     &        hBHxYeng128, hBHxY128force)
c
       fcall = fcall + 1
c
        do i = 1,3*natom
        atomForceF1(i) = vbdefForce(i)*fEngWF(1)  
     &                   + vAngdefForce(i)*fEngWF(2) 
     &                   + impDefForce(i)*fEngWF(3)
     &                   + torsAngForce(i)*fEngWF(4)
     &                   + vdwForceR1(i)*fEngWF(5)
     &                   +  coulForceR1(i)*fEngWF(6)
     &                   + (restr1AtForce(i)
     &                   + restr1MHWAtForce(i)
     &                   + restr1RHWAtForce(i)
     &                   + restrDistA2Force(i))*fEngWF(8)
     &                   + hBHxY128force(i)
        end do !i
c
        engFastF1 = eVbondDef*fEngWF(1) + eVangDef*fEngWF(2)
     &              + eImpDef*fEngWF(3) + eTorsDef*fEngWF(4) 
     &              + engVDWR1*fEngWF(5) + engCOULR1*fEngWF(6)
     &              + (restr1Eng+restr1MHWEng+restr1RHWEng
     &              + restrDistA2Eng)*fEngWF(8)
     &              + hBHxYeng128                
c                     !hb128Eng/F are scaled in getHbTripletAll
c
c	call zeroAllAtomForce(nmoveatom,moveAtomList,atomForceF1)
c
         if(CONTROL)then
         write(kanalp,*)
     &   'initAllForce01:Res:  atomXYZ  atomForceF1: fcall:',fcall
         do i = 1,natom
         i3=3*i-3
         write(kanalp,'(i5,2x,3f8.3,2x,3f8.3)') 
     &   i,(atomXYZ(i3+k),k=1,3),(atomForceF1(i3+k),k=1,3)
         end do! i 
         end if !C
c
c         STOP
c
         return
         end
