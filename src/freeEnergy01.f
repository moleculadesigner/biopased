c free Energy   
c
c Yury Vorobjev, 2005
c
c calculation of free energy over md trajectory
c
	subroutine mdFreeEnergy01 
c
cx        implicit none
        include 'xyzPDBsize.h'
        include 'xyzPDBinfo.h' 
        include 'xyzPDBcrd.h'
        include 'xyzPDBopt.h'
        include 'engXYZtra.h'
        include 'modeXYZtra.h'
        include 'filedat.h'
        include 'optionPar.h'
        include 'dataSASdr.h'
        include 'nbondPairVCS.h'
        include 'movingAtom.h'
        include 'enForce.h'
c
        include "output.h"
        include "kanalUse.h"
        include "statusMessg_mDyn.h"
c
	integer i,j,ia,ia3,ja,ja3
        integer i3,j3,is,iSH,k
        integer iwe,iwx
        integer na3,iam,im3
        real invSt,an1,an2
c
        logical CONTROL,CONTROL1
        integer kanalp
c
        kanalp = kanalRunOut
        CONTROL = .true.
        CONTROL1 = .false.
c
        iStatus = iStatus + 1
        write(kanalPStat,*)mStatus,iStatus,messageStatus
     &  ,' freeEnergy01 for mdTra start ...'
c
        engTablFile2 = 'engTable.dat'
        if(nMolNameLett .ge. 1)then
        engTablFile2=molNameArgLine(1:nMolNameLett)//'.engTable.dat'
        end if
c
        if(OPT_resFileDir)then 
        engTablFile2=resultFileDir(1:nResFileDirLett)//'engTable.dat'
        if(nMolNameLett .ge. 1)then
        engTablFile2=resultFileDir(1:nResFileDirLett)//
     &  molNameArgLine(1:nMolNameLett)//'.engTable.dat'
        end if
        end if
c
         open(unit=engTablkanaL2, file=engTablFile2, form='formatted',
     &        status= 'unknown')
c
        if(CONTROL)then
        write(kanalp,*)'mdFreeEnergy01 start: '
        write(kanalp,*)'nSnap:',is
        end if
c
c init OPT_ for GeneralizedBorn model
c
        OPT_SolvGBorn = .true.
        OPT_SolvEFS = .true.
        OPT_SolvSASHP = .true.
        OPT_SolvGS = .false.      ! Gaussian Shell solvation model is off
c
c init:
          engVDWR1TraSd = 0.0
          engVDWR1TraAv = 0.0
          engCOULTraAv = 0.0
          engCOULTraSd = 0.0
          eGeoDefTraAv = 0.0
          eGeoDefTraSd = 0.0
          engMolSolAv = 0.0
          engMolSolSd = 0.0
          engPOTENTTraAv = 0.0
          engPOTENTTraSd = 0.0
          hBHxYeng128TraAv = 0.0
          hBHxYeng128TraSd = 0.0
c
         essModeEntropyMdTra = essModeEntropy    ! T*EssModeEntropy (kcal/m) 
c
c nStepTra : accumulated md trajectory
c
        na3 = 3*natom
	do is = 1,nStepTra
        iSH = na3*(is-1) 
c take traj Snap XYZ
        do ia3=1,na3
        atomXYZ(ia3) = atomXYZtra(ia3+iSH)
        end do!ia3
c
c calculate energy terms for the snap XYZ
c
        nSAScall = 0            ! update SAS for each snapshot
        ncallNBPL = 0
c make atomXYZopt() = atomXYZ()
        call initSwap01(natom,atomXYZ,atomXYZopt)
c
        do iam=1,nmoveatom
        im3=iam*3-3
        i3 = 3*moveAtomList(iam)-3
        moveAtomXYZ(im3+1)=atomXYZopt(i3+1)
        moveAtomXYZ(im3+2)=atomXYZopt(i3+2)
        moveAtomXYZ(im3+3)=atomXYZopt(i3+3)
        end do ! iam
c
        iwe = 1
        iwx = 1
        call engOptFunct(moveAtomXYZ,engPOTENT,atForceTot,iwe,iwx)
c
        iStatus = iStatus + 1
        write(kanalPStat,*)mStatus,iStatus,messageStatus
     &  ,' snapXYZ :',is,' energy calculation done ...'
c
c energy terms for snapshot XYZ:
c
          eVbondDefTra(is) = eVbondDef
          eVangDefTra(is) = eVangDef
          eImpDefTra(is) = eImpDef
          eTorsDefTra(is) = eTorsDef
          engVDWR1Tra(is) = engVDWR1
          hBHxYeng128Tra(is) = hBHxYeng128
          solvMolEngF3 = molSolEn+hpSoLEnergy+watBrgEnergySoL
     &                  + bornPolzEng
          eMolSolvTra(is) = solvMolEngF3
          engCOULR1Tra(is) = engCOULR1
          engCOULTra(is) = engCOUL
          engPOTENTTra(is) =  engPOTENT
          eGeoDefTra(is) = eGeoDef
c
          invSt = 1.0/is
          an1 = 1.0 - invSt
          an2 = (is - 1.0)/(1.0+is)**2
c
          engVDWR1TraAv = engVDWR1TraAv
     &                    + (engVDWR1-engVDWR1TraAv)*invSt
          engVDWR1TraSd = an1*engVDWR1TraSd
     &                    + an2*(engVDWR1-engVDWR1TraAv)**2
          engCOULTraAv = engCOULTraAv
     &                   +  (engCOUL-engCOULTraAv)*invSt
          engCOULTraSd = an1*engCOULTraSd
     &                   +  an2*(engCOUL - engCOULTraAv)**2
          engPOTENTTraAv = engPOTENTTraAv
     &                   + (engPOTENT-engPOTENTTraAv)*invSt
          engPOTENTTraSD = an1*engPOTENTTraSD
     &                   + an2*(engPOTENT-engPOTENTTraAv)**2
          hBHxYeng128TraAv = hBHxYeng128TraAv
     &                   +  (hBHxYeng128 - hBHxYeng128TraAv)*invSt
          hBHxYeng128TraSd = an1*hBHxYeng128TraSd
     &                   + an2*(hBHxYeng128 - hBHxYeng128TraAv)**2
          engMolSolAv = engMolSolAv + (solvMolEngF3 - engMolSolAv)*invSt
          engMolSolSd = an1*engMolSolSd +
     &                              an2*(solvMolEngF3 - engMolSolAv)**2
c
        end do!is
c
       freeEngMdTra = engPOTENTTraAv - essModeEntropyMdTra
c
c atomicRMSD from the FIRST snap for mdTra
c
        do is=1,nStepTra
        iSH=na3*(is-1)
        atomXYZRmsdTra(is)=0.0
        do i=1,na3
        atomXYZRmsdTra(is)=atomXYZRmsdTra(is)+
     &     (atomXYZtra(i+iSH)-atomXYZtra(ia3))**2
        end do!i
        atomXYZRmsdTra(is)=sqrt(atomXYZRmsdTra(is)/natom)
        end do !is
c
        if(CONTROL)then
        write(kanalp,*)'freeEnergy01: iSnap atomXYZRmsdTra:'
        do is=1,is 
        write(kanalp,*)'atomXYZRmsdTra:',is,atomXYZRmsdTra(is)
        end do!is
        end if
c
c write engTablFile2
c
        write(engTablkanaL2,'(a27)')'#ENERGY TERMS vs snapShot  '
        write(engTablkanaL2,'(a50)')
     &  '# is     ePot     eVDW     eCOL     eHbond   eSOLV '
c
        do is = 1,nStepTra
        write(engTablkanaL2,'(i5,2x,5f9.1)')
     &  is,engPOTENTTra(is),engVDWR1Tra(is),engCOULTra(is),
     &  hBHxYeng128Tra(is),eMolSolvTra(is)
c
        end do!is
c 
        write(engTablkanaL2,'(a7,5f9.1)')
     & '#aver: ',engPOTENTTraAv,engVDWR1TraAv,engCOULTraAv,
     &  hBHxYeng128TraAv,engMolSolAv
c
         write(engTablkanaL2,'(a7,5f9.1)')
     & '#sdev: ',sqrt(engPOTENTTraSd),sqrt(engVDWR1TraSd),
     &  sqrt(engCOULTraSd),sqrt(hBHxYeng128TraSd),
     &  sqrt(engMolSolSd)
c 
        write(engTablkanaL2,*)            
c
        write(engTablkanaL2,'(a20,f9.1)')
     &  'essEntropy(-T*S) = ', -essModeEntropyMdTra
c
        write(engTablkanaL2,'(a20,f9.1,a9)')
     &  'mdTraFreeEng     = ',freeEngMdTra, ' kcal/mol'
c
	return
7071    format(a6,1x,i4,2x,a4,a4,a1,i4,4x,3f8.3,7x,f8.5) ! PDB rq  
        end
