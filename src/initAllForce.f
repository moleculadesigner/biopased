c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c                                                                           *
c  initialize: all forces on atoms                                          * 
c                                                                           *
c Yury Vorobjev    2002     
c                  2005                                                     *
c                  2009 P-exWATint:PWat
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c USed for energy optimization 
c

	subroutine initAllForce(fcall,atomXYZ,
     &              doFastF1,doMedF2,doSlowF3,
     &              eVbondDef,vbdefForce,
     &              eVangDef,vAngdefForce,
     &              eImpDef,impDefForce,
     &              eTorsDef,torsAngForce,
     &              engVDWR1,vdwForceR1,
     &              engCOULR1,coulForceR1,
     &              hBHxYeng128,hBHxY128force,
     &              engCOULR2,coulForceR2,
     &              restr1Eng,restr1AtForce,
     &              restr1MHWEng,restr1MHWAtForce,
     &              restr1RHWEng,restr1RHWAtForce,
     &              restrDistA2Eng,restrDistA2Force,
     &              molSolEn,atomSolEn,atomSolFr,
     &              watBrgEnergySoL,watBrgForceAt,
     &              hpSoLEnergy,hpSoLForceAt,
     &              bornPolzEng,bornPolzEngAt,bornPolzForceAt)
c
        include 'xyzPDBsize.h'
        include 'xyzPDBinfo.h'
c
        include "output.h"
        include "kanalUse.h"
        include "statusMessg_mDyn.h"
c
        include 'pair1234array.h'
        include 'nbondPairVCS.h'
        include 'vdw12Par.h'
c        include 'enForce.h'
        include 'restrainInfo.h'
        include 'loopInfo.h'
        include 'movingAtom.h'
        include 'solvGSarray.h'
        include 'optionPar.h'
        include 'shake.h'
        include 'restrain1MHW.h'
        include 'restrain1RHW.h'
        include 'restrDistA2.h'
        include 'hbond128.h'
        include 'solvate01.h'
c
        real atomXYZ(*)
        integer fcall    ! sequemtial number of call
        integer doFastF1,doMedF2,doSlowF3 
c doFastF1,doMedF2,doSlowF3 - flags to update en/forces
        real eVbondDef,vbdefForce(*)
        real eVangDef,vAngdefForce(*)
        real eImpDef,impDefForce(*)
        real eTorsDef,torsAngForce(*)
        real engVDWR1,vdwForceR1(*)
        real engCOULR1,coulForceR1(*)
        real hBHxYeng128,hBHxY128force(*)
        real engCOULR2,coulForceR2(*)
        real restr1Eng,restr1AtForce(*)
        real restr1MHWEng,restr1MHWAtForce(*)
        real restr1RHWEng,restr1RHWAtForce(*) 
        real restrDistA2Eng,restrDistA2Force(*)
        real molSolEn, atomSolEn(*), atomSolFr(*)
        real watBrgEnergySoL,watBrgForceAt(*) 
        real hpSoLEnergy,hpSoLForceAt(*)
        real bornPolzEng,bornPolzEngAt(*),bornPolzForceAt(*)
c
       logical OPT_SoftCore
       integer i,kanalp
c
       kanalp =  kanalRunOut
       OPT_SoftCore = .true.      ! use subroutine allAtVDWEForceR1SC()
c
c       write(kanalp,*)'initAllForce start!: fcall=',fcall
c
c init to zero All forces on Atoms
       if(fcall .eq. 0)then
       eVbondDef = 0.0
       eVangDef = 0.0      
       eImpDef = 0.0
       eTorsDef = 0.0
       restr1Eng = 0.0            
       restr1MHWEng = 0.0
       restr1RHWEng = 0.0 
       restrDistA2Eng = 0.0 
       engVDWR1 = 0.0
       engCOULR1 = 0.0
       engCOULR2 = 0.0
       molSolEn = 0.0
       hBHxYeng128 = 0.0
       watBrgEnergySoL = 0.0
       hpSoLEnergy = 0.0 
       bornPolzEng = 0.0 
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
       coulForceR2(i) = 0.0
       atomSolFr(i) = 0.0
       hBHxY128force(i) = 0.0
       watBrgForceAt(i) = 0.0
       hpSoLForceAt(i) = 0.0
       bornPolzForceAt(i) = 0.0
       end do !i
       end if !fcall
c
c all GeoDef forces are calculated at each step
       if(doFastF1 .eq. 1) then
 	call allAtVBondEForce(atomXYZ,
     &           natom,bond12List,nbond12,nbond12noH,iShake,
     &           bond12ParL,eVbondDef,vbdefForce )
c
	call allAtVangEForce(atomXYZ,
     &           natom,trip123List,nTrip123,ang123ParL,
     &           eVangDef,vAngdefForce )
c
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
c restr1MolHWAtForce 
c
c        write(kanalp,*)'initAllForce : call refMolPosRestr1HWEF:'
c
        call refMolPosRestr1HWEF(atomXYZ,natom,
     &           nRestr1MHWatom,restr1MHWatList,atMs1MHWatList,
     &           ref1MHWPos,restr1MHWConst,sizeRestr1MHW,
     &           restr1MHWEng,restr1MHWAtForce)
c
c       write(kanalp,*)'initAllForce : finis refMolPosRestr1HWEF:'
c
c restr1ResHWAtForce
c
c        write(kanalp,*)'initAllForce : call refResPosRestr1HWEF:'
c
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
c restrDistA2EngForces - a1-a2 distant restrain 
        call allAtDistA2EnForce(atomXYZ,
     &           natom,nRestrDistA2,restrDistA2List,
     &           restrDistA2DistHK,restrDistA2Eng,restrDistA2Force)
c
c        write(kanalp,*)'initAllForce: call allAtDistA2EnForce: '        
c
         if (doSlowF3 .eq. 1) then 
c init BornRad before ElectrostaticEngForce calculation
         if (OPT_SolvGBorn)then
c init SAS + bornRadAt + bornRadAtDr
         call initBornRadSAS(atomXYZ)
         end if !OPT_SolvBorn 
         end if !doSlowF3 .eq. 1        
c
cx        if(OPT_SoftCore)then
        call  allAtVDWEForceR1SC(atomXYZ,atomQ,
     &        natom,nmoveatom,moveAtomList,
     &        resNumb,atomBlockName, 
     &        nbpairListV,startnbPairLV,nnbPairLV,
     &        pair14VList,startPairL14V,nPairL14V, 
     &        nVDWtype,atomVDWtype,atomVDW12ab,aSoftCore,
     &        rcutV,rcutC,engVDWR1,vdwForceR1,engCOULR1,coulForceR1,
     &        hBHxYeng128, hBHxY128force)
c
cx          else
c
cx        call  allAtVDWEForceR1(atomXYZ,atomQ,
cx     &        natom,nmoveatom,moveAtomList,
cx     &        resNumb,atomBlockName, 
cx     &        nbpairListV,startnbPairLV,nnbPairLV,
cx     &        pair14VList,startPairL14V,nPairL14V, 
cx     &        nVDWtype,atomVDWtype,atomVDW12ab,
cx     &        rcutV,rcutC,engVDWR1,vdwForceR1,engCOULR1,coulForceR1,
cx     &        hBHxYeng128, hBHxY128force)
cx
cx        end if ! OPT_SoftCore
c
        end if !doFastF1 (shortR1)
c
        if(doMedF2 .eq. 1) then

        call allAtVDWEForceR2(atomXYZ,atomQ,
     &           natom,nmoveatom,moveAtomList,
     &           nbpairListC,startnbPairLC,nnbPairLC,
     &           rcutV,rcutC,engCOULR2,coulForceR2)

         end if ! doMedF2:coulomb R2
c
c solvent Forces
         if (doSlowF3 .eq. 1) then
         if(OPT_SolvGS)then
         call SolventEnForces(natom, atomXYZ,
     &         atomName,startPairL12,nPairL12,pair12List,
     &         nbpairListS,startnbPairLS,nnbPairLS,
     &         atomSolPar, molSolEn, atomSolEn, atomSolFr)
         end if !GaussSolvModel
c
         if(OPT_SolvGBorn)then
cx         call initBornRadSAS(atomXYZ) ! it is alredy done
           call getBornPolzEforce(atomXYZ,
     &               bornPolzEng,bornPolzEngAt,bornPolzForceAt)
         end if !GBorn
c
         if(OPT_SolvEFS)then
c calculates dataSASdr.h
         call solvateMoL03(atomXYZ)      !take care donot init SAS twice!!
c                                         if OPT_SolvBorn=.true.
c
         call getWatBrgForceSAS03(atomXYZ,
     &                   watBrgEnergySoL,watBrgForceAt)
c
          end if !OPT_SolvEFS
c
         if(OPT_SolvSASHP)then
c use SAS data from dataSASdr.h         
         call getHPolvForceSAS03(atomXYZ,
     &                   hpSoLEnergy,hpSoLForceAt)
c
         end if !OPT_SolvSASHP
c
         end if !doSlowF3
c
       fcall = fcall + 1
c
         return
         end
