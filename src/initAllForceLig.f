c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c                                                                           *
c  calculate: all forces on LIGand atoms                                    * 
c             and energy of the LIGand with the  MOLECULE                   *
c                                                                           *
c Yury Vorobjev    2003                                                     *
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c USed for calculation Lig-Protein interactions
c

	subroutine initAllForceLig(fcall,atomXYZ,
     &              eVbondDef,vbdefForce,
     &              eVangDef,vAngdefForce,
     &              eImpDef,impDefForce,
     &              eTorsDef,torsAngForce,
     &              engVDWR1,vdwForceR1,
     &              hBHxYeng128,hBHxY128force,
     &              engCOULR1,coulForceR1,
     &              engCOULR2,coulForceR2,
     &              restr1Eng,restr1AtForce,
     &              restr1MHWEng,restr1MHWAtForce,
     &              molSolEn,atomSolEn,atomSolFr)
c
        include 'xyzPDBsize.h'
        include 'xyzPDBinfo.h'
c
        include "output.h"
cX        include 'xyzPDBcrd.h'
        include 'pair1234array.h'
        include 'nbondPairVCS.h'
        include 'vdw12Par.h'
cX        include 'enForce.h'
        include 'restrainInfo.h'
        include 'loopInfo.h'
cX        include 'movingAtom.h'
        include 'ligInfo.h'
        include 'solvGSarray.h'
        include 'optionPar.h'
        include 'shake.h'
        include 'restrain1MHW.h'
        include 'hbond128.h'
c
        real atomXYZ(*)
        integer fcall    ! sequemtial number of call
        integer k
        real eVbondDef,vbdefForce(*)
        real eVangDef,vAngdefForce(*)
        real eImpDef,impDefForce(*)
        real eTorsDef,torsAngForce(*)
        real engVDWR1,vdwForceR1(*)
        real engCOULR1,coulForceR1(*)
        real engCOULR2,coulForceR2(*)
        real restr1Eng,restr1AtForce(*)
        real restr1MHWEng,restr1MHWAtForce(*)
        real molSolEn, atomSolEn(*), atomSolFr(*)
        real hBHxYeng128,hBHxY128force(*)
c
       logical CONTROL
       logical OPT_SoftCore
       integer i,kanalp
c
       kanalp = kanalRunOut
       CONTROL = .false.
       OPT_SoftCore = .true.      ! use subroutine allAtVDWEForceR1SC()
c       OPT_SoftCore = .false.
c
       if(CONTROL)then
       write(kanalp,*)'initAllForceLig start!: fcall=',fcall
       end if
c
c init to zero All forces on Atoms
       if(fcall .eq. 0)then
       eVbondDef = 0.0
       eVangDef = 0.0      
       eImpDef = 0.0
       eTorsDef = 0.0
       restr1Eng = 0.0            
       restr1MHWEng = 0.0
       engVDWR1 = 0.0
       engCOULR1 = 0.0
       engCOULR2 = 0.0
       molSolEn = 0.0
       hBHxYeng128 = 0.0
       
c
       do i = 1,3*natom
       vbdefForce(i) = 0.0
       vAngdefForce(i) = 0.0
       impDefForce(i) = 0.0
       torsAngForce(i) = 0.0
       restr1AtForce(i) = 0.0
       restr1MHWAtForce(i) = 0.0
       vdwForceR1(i) = 0.0 
       coulForceR1(i) = 0.0
       coulForceR2(i) = 0.0
       atomSolFr(i) = 0.0
       hBHxY128force(i) = 0.0
       end do !i
       end if !fcall
c
c all GeoDef forces are calculated at each step
 	call allAtVBondEForceLig(atomXYZ,
     &           vLigFlag,bond12LigFlag,
     &           natom,bond12List,nbond12,nbond12noH,iShake,
     &           bond12ParL,eVbondDef,vbdefForce )
c
	call allAtVangEForceLig(atomXYZ,
     &           vLigFlag,trip123LigFlag,     
     &           natom,trip123List,nTrip123,ang123ParL,
     &           eVangDef,vAngdefForce )
c
c
	call allAtImpTEForceLig(atomXYZ,
     &           vLigFlag,qImp1234LigFlag,
     &           natom,quarImp1234L,nImp1234,impAng1234ParL,
     &           eImpDef,impDefForce )
c
c torsionEnForces
c
        call allAtTorsEForceLig(atomXYZ,
     &           vLigFlag,quar1234LigFlag,  
     &           natom,quar1234List,nQuar1234,
     &           quar1234ParL,quar1234nPar,
     &           eTorsDef,torsAngForce )
c
c restrain1atomForce
c
          restr1Eng = 0.0
cx        call refAtPosRestr1EF(atomXYZ,
cx     &           natom,nRestr1atom,restr1atomList,
cx     &           refAtomPos,restr1AtConst,restr1Eng,restr1AtForce)
c
c restr1MolHWAtForce 
c
        write(kanalp,*)'initAllForce : call refMolPosRestr1HWEF:'
c
        restr1MHWEng = 0.0
cx        call refMolPosRestr1HWEF(atomXYZ,natom,
cx     &           nRestr1MHWatom,restr1MHWatList,atMs1MHWatList,
cx     &           ref1MHWPos,restr1MHWConst,sizeRestr1MHW,
cx     &           restr1MHWEng,restr1MHWAtForce)
c
         if(CONTROL)then
         write(kanalp,*)'initAllForce : finis refMolPosRestr1HWEF:'
         write(kanalp,*)'nAtomInLig:',nAtomInLig
         write(kanalp,*)'atomInLigList:'
         write(kanalp,'(10i6)')(atomInLigList(k),k=1,nAtomInLig)
         end if !CL
c
cx        if(OPT_SoftCore)then
c
c USE: ligandAt-moleculeAt nonBonded pairList's from ligInfo.h
c
        call  allAtVDWEForceR1SC(atomXYZ,atomQ,
     &        natom,nAtomInLig,atomInLigList,
     &        resNumb,atomBlockName, 
     &        nbpairListVLig,startnbPairLVLig,nnbPairLVLig,
     &        pair14VList,startPairL14V,nPairL14V, 
     &        nVDWtype,atomVDWtype,atomVDW12ab,aSoftCore,
     &        rcutV,rcutC,engVDWR1,vdwForceR1,engCOULR1,coulForceR1,
     &        hBHxYeng128, hBHxY128force)
c
cx          else
c
cx        call  allAtVDWEForceR1(atomXYZ,atomQ,
cx     &        natom,nAtomInLig,atomInLigList, 
cx     &        resNumb,atomBlockName, 
cx     &        nbpairListVLig,startnbPairLVLig,nnbPairLVLig,
cx     &        pair14VList,startPairL14V,nPairL14V, 
cx     &        nVDWtype,atomVDWtype,atomVDW12ab,
cx     &        rcutV,rcutC,engVDWR1,vdwForceR1,engCOULR1,coulForceR1,
cx     &        hBHxYeng128, hBHxY128force)
c
cx        end if ! OPT_SoftCore
c
        call allAtVDWEForceR2(atomXYZ,atomQ,
     &           natom,nAtomInLig,atomInLigList, 
     &           nbpairListCLig,startnbPairLCLig,nnbPairLCLig,
     &           rcutV,rcutC,engCOULR2,coulForceR2)
c
c solvent Forces
         if (OPT_SolvGS ) then
         call SolventEnForces(natom, atomXYZ,
     &         atomName,startPairL12,nPairL12,pair12List,
     &         nbpairListSLig,startnbPairLSLig,nnbPairLSLig,
     &         atomSolPar, molSolEn, atomSolEn, atomSolFr)
c
         end if !OPT_SolvGS
c
         return
         end
