c engergyOptimizationFunctional
c
c Yuri Vorobjev  2003                           
c                2004   hb128
c                2005   WBridge Solvation
c                2005   GBorn solvation
c
	subroutine engOptFunct(xyz,eVal,eGrad,iwe,iwx)
c 
c engOptFunct calculates energy and gradEn mdatomXYZvel.h 
c for a given set of optimized xyz
c
        include 'xyzPDBsize.h'
c
        include "output.h"
        include "kanalUse.h"
        include "statusMessg_mDyn.h"
c
        include 'xyzPDBinfo.h'
        include 'xyzPDBopt.h'
        include 'movingAtom.h'
        include 'enForce.h'
        include 'filedat.h'  
        include 'optionPar.h'
        include 'pair1234array.h' 
        include 'nbondPairVCS.h'
        include 'hbond128.h'
        include 'compGeoEF.h'
c
        real xyz(*)
        real eVal
        real eGrad(*)
        integer iwx,iwe    !write flag
clocal
cx       integer makeVdW,makeCL,makeSL
        integer doVDWef,doCLef,doSLef
        integer ia,iam,im3,k,k3,i3,i,j
        integer ich
        real atForceMX
c
        real eVbondDefLg,eVangDefLg
        real eImpDefLg,eTorsDefLg
        real engVDWR1Lg,engCOULR1Lg
        real engCOULR2Lg
        real restr1MHWEngLg,restr1EngLg
        real molSolEnLg,eGeoDefLg
        real eGeoRstLg,engCOULLg
        real engPOTENTLg
        real hBHxYeng128Lg
c
        integer kanalp
        logical CONTROL,CONTROL1,CONTROLx
c
cx        kanalp = 6
        kanalp = kanalRunOut
c
        CONTROL = .true.     
        CONTROL1 = .false.    
        CONTROLx = .false.
cx        kanalPDBfin = 24
        atForceMX = 0.0 
c
        if(CONTROL)then
        write(kanalp,*)'engOptFunct: starts: * * * * * * * * * * '
        write(kanalp,*) 'EnergyTerms:fEngWF(*):',fEngWF
        write(kanalp,*)'engOptFunctBS: nmoveatom :', nmoveatom
        end if
c
        makeVdW = 1
        makeCL = 1
        makeSL = 1
c
        doVDWef = 1
        doCLef = 1
        doSLef = 1
c
c update atomXYZopt()
        do iam=1,nmoveatom
        im3=iam*3-3
        i3 =3*moveAtomList(iam)-3
        atomXYZopt(i3+1)=xyz(im3+1)
        atomXYZopt(i3+2)=xyz(im3+2)
        atomXYZopt(i3+3)=xyz(im3+3)
        end do ! iam
c
        if(CONTROLx)then
        write(kanalp,*)'engOptFunctBS:Start: atomXYZopt:'
        i3 = 3*nmoveatom 
        call print9r(kanalp,i3,atomXYZopt)
        end if!Contr
c
       call initNonBondList(atomXYZopt,makeVdW,makeCL,makeSL) 
c
c update hB128 tripletList
        if(OPT_HBond128)then
        call getHbTripletAll(atomXYZopt,
     &           natom,nmoveatom,moveAtomList,
     &           nbpairListV,startnbPairLV,nnbPairLV,
     &           rcutHb128,nHBbondHxYmx,
     &           nHBbondHxY,hB128List,hB128TypeList,hB128PotParList)
c        
      do i = 1,nHBbondHxY
      if(hB128TypeList(i).eq.1)hB128WfList(i)=fEngWF(10)*hB128WfScaleALL
      if(hB128TypeList(i).eq.2)hB128WfList(i)=fEngWF(11)*hB128WfScaleALL
        end do !i
        end if ! hb128
c
        call initSASdielModLaZ(atomXYZopt)
c
c artificial Compactization forces
c
         if(OPT_CompactForce)then
         call compactGeoCntEF(atomXYZopt,
     &           natom,nmoveatom,moveAtomList,
     &           compactGeoFpar,compactGeoEn,compactGeoForce)
         end if ! OPT_CompactForce
c
c update forces/energy
	call initAllForce(fcall,atomXYZopt,
     &              doVdWef,doCLef,doSLef,
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
     &              molSolEn, atomSolEn, atomSolFr,
     &              watBrgEnergySoL,watBrgForceAt,
     &              hpSoLEnergy,hpSoLForceAt,
     &              bornPolzEng,bornPolzEngAt,bornPolzForceAt)
c
          eGeoDef   = fEngWF(1)*eVbondDef + fEngWF(2)*eVangDef 
     &                + fEngWF(3)*eImpDef + fEngWF(4)*eTorsDef
c
          eGeoRst   = fEngWF(8)*(restr1Eng+restr1MHWEng+restrDistA2Eng
     &                + restr1RHWEng) + compactGeoEn*fEngWF(9)
c
          engCOUL   = fEngWF(6)*engCOULR1 + fEngWF(7)*engCOULR2
c
          solvMolEngF3 = (molSolEn + watBrgEnergySoL+ hpSoLEnergy +
     &                   bornPolzEng)*fEngWF(9)
c
c eGeoRst is not accounted in the total potentialEnergy: engPOTENT
c
          engPOTENT = eGeoDef + fEngWF(5)*engVDWR1 + engCOUL  
     &        + solvMolEngF3 
     &        + hBHxYeng128       ! hBHxYeng128 is scaled in allAtNBondEForceBS.f
c
          eVal = engPOTENT
c
        if(CONTROL)then
       write(kanalp,'(a45,6f9.1)') 
     &  'optRun: ePOT,eDef,eVDWR1,eCOUL,eSolv,eRst: ',
     &  engPOTENT,eGeoDef,engVDWR1,engCOUL,solvMolEngF3,eGeoRst
        end if! C
c
c assign eGrad()
        do iam=1,nmoveatom
        im3=iam*3-3
        i3 =3*moveAtomList(iam)-3
        do k=1,3
        eGrad(im3+k) =   
     &  - fEngWF(1)*vbdefForce(i3+k)
     &  - fEngWF(2)*vAngdefForce(i3+k)
     &  - fEngWF(3)*impDefForce(i3+k)
     &  - fEngWF(4)*torsAngForce(i3+k)
     &  - fEngWF(5)*vdwForceR1(i3+k)
     &  - fEngWF(6)*coulForceR1(i3+k)
     &  - fEngWF(7)*coulForceR2(i3+k)
     &  - fEngWF(8)*(restr1AtForce(i3+k)
     &  -            restr1MHWAtForce(i3+k) 
     &  -            restr1RHWAtForce(i3+k)
     &  -            restrDistA2Force(i3+k))
     &  - fEngWF(9)*(atomSolFr(i3+k) 
     &  - compactGeoForce(i3+k)
     &  - watBrgForceAt(i3+k)+hpSoLForceAt(i3+k)
     &  - bornPolzForceAt(i3+k))
     &  - hBHxY128force(i3+k)  
c
        if(atForceMX .lt. abs(eGrad(im3+k)))
     &   atForceMX = abs(eGrad(im3+k))
        end do!k
        end do ! iam
c
c write energyTraj
          if(CONTROL1 .or. iwe .eq. 1)then
          write(kanalp,'(a8)')'$ENERGY: '
          write(kanalp,'(a10,f14.5)')'eVbondDef:',eVbondDef
          write(kanalp,'(a10,f14.5)')'eVangDef :',eVangDef
          write(kanalp,'(a10,f14.5)')'eImpDef  :',eImpDef
          write(kanalp,'(a10,f14.5)')'eTorsDef :',eTorsDef
          write(kanalp,'(a10,f14.5)')'engVDWR1 :',engVDWR1
          write(kanalp,'(a12,f14.5)')'engHb128 :',hBHxYeng128
          write(kanalp,'(a10,f14.5)')'engCOULR1:',engCOULR1
          write(kanalp,'(a10,f14.5)')'engCOULR2:',engCOULR2
          write(kanalp,'(a10,f14.5)')'eRstH1At :',restr1Eng
          write(kanalp,'(a10,f14.5)')'eRstHW1RS:',restr1RHWEng
          write(kanalp,'(a10,f14.5)')'eRstHW1ML:',restr1MHWEng
          write(kanalp,'(a10,f14.5)')'eRstH2At :',restrDistA2Eng
          write(kanalp,'(a10,f14.5)')'eCompcEng:',compactGeoEn
          write(kanalp,'(a10,f14.5)')'eGeoRst  :',eGeoRst 
          write(kanalp,'(a10,f14.5)')'eGeoDef  :',eGeoDef
          write(kanalp,'(a10,f14.5)')'eCOUL    :',engCOUL
          write(kanalp,'(a10,f14.5)')'eSolvGS  :',molSolEn
          write(kanalp,'(a10,f14.5)')'eSolvCav :',hpSoLEnergy
          write(kanalp,'(a10,f14.5)')'eSolvGBp :',bornPolzEng 
          write(kanalp,'(a10,f14.5)')'eSolvWbr :',watBrgEnergySoL 
          write(kanalp,'(a10,f14.5)')'eSolvTot :',solvMolEngF3
          write(kanalp,'(a10,f14.5)')'ePOTENT  :',engPOTENT
          write(kanalp,'(a10,f14.5)')'atForceMX:',atForceMX
          write(kanalp,'(a4)')'$END '
          end if !control1
c
c !write result PDB file    
        if(iwx .eq. 1) then   
        open(unit=kanalPDBfin, file=pdbOptFile, form='formatted',
     &       status='unknown')
c
          write(kanalPDBfin,'(a10,f14.5)')'eVbondDef:',eVbondDef
          write(kanalPDBfin,'(a10,f14.5)')'eVangDef :',eVangDef
          write(kanalPDBfin,'(a10,f14.5)')'eImpDef  :',eImpDef
          write(kanalPDBfin,'(a10,f14.5)')'eTorsDef :',eTorsDef
          write(kanalPDBfin,'(a10,f14.5)')'engVDWR1 :',engVDWR1
          write(kanalPDBfin,'(a10,f14.5)')'engHb128 :',hBHxYeng128
          write(kanalPDBfin,'(a10,f14.5)')'engCOULR1:',engCOULR1
          write(kanalPDBfin,'(a10,f14.5)')'engCOULR2:',engCOULR2
          write(kanalPDBfin,'(a10,f14.5)')'eRstH1At :',restr1Eng
          write(kanalPDBfin,'(a10,f14.5)')'eRstHW1RS:',restr1RHWEng
          write(kanalPDBfin,'(a10,f14.5)')'eRstHW1ML:',restr1MHWEng
          write(kanalPDBfin,'(a10,f14.5)')'eRstH2At :',restrDistA2Eng
          write(kanalPDBfin,'(a10,f14.5)')'eCompcEng:',compactGeoEn
          write(kanalPDBfin,'(a10,f14.5)')'eGeoRst  :',eGeoRst
          write(kanalPDBfin,'(a10,f14.5)')'eGeoDef  :',eGeoDef
          write(kanalPDBfin,'(a10,f14.5)')'eCOUL    :',engCOUL
          write(kanalPDBfin,'(a10,f14.5)')'eSolvGS  :',molSolEn
          write(kanalPDBfin,'(a10,f14.5)')'eSolvCav :',hpSoLEnergy
          write(kanalPDBfin,'(a10,f14.5)')'eSolvGBp :',bornPolzEng
          write(kanalPDBfin,'(a10,f14.5)')'eSolvWbr :',watBrgEnergySoL
          write(kanalPDBfin,'(a10,f14.5)')'eSolvTot :',solvMolEngF3
          write(kanalPDBfin,'(a10,f14.5)')'ePOTENT  :',engPOTENT
          write(kanalPDBfin,'(a10,f14.5)')'atForceMX:',atForceMX
          write(kanalPDBfin,'(a4)')'$END '
c
c write energy:Lig-Prot
c
          vLigFlag = 0
          if(OPT_LigRes)vLigFlag = 1
          if(vLigFlag .eq. 1) then
c
       call initNonBondListLig(atomXYZopt,makeVdW,makeCL,makeSL)
c
       call initAllForceLig(fcall,atomXYZopt,
     &              eVbondDefLg,vbdefForce,
     &              eVangDefLg,vAngdefForce,
     &              eImpDefLg,impDefForce,
     &              eTorsDefLg,torsAngForce,
     &              engVDWR1Lg,vdwForceR1,
     &              hBHxYeng128Lg,hBHxY128force, 
     &              engCOULR1Lg,coulForceR1,
     &              engCOULR2Lg,coulForceR2,
     &              restr1EngLg,restr1AtForce,
     &              restr1MHWEngLg,restr1MHWAtForce,
     &              molSolEnLg,atomSolEn,atomSolFr)
c
c sum up Ligand-Prot Energies:
        eGeoDefLg   = fEngWF(1)*eVbondDefLg + fEngWF(2)*eVangDefLg
     &                + fEngWF(3)*eImpDefLg + fEngWF(4)*eTorsDefLg
c
        eGeoRstLg   = fEngWF(8)*restr1EngLg+fEngWF(8)*restr1MHWEngLg
c
        engCOULLg   = fEngWF(6)*engCOULR1Lg + fEngWF(7)*engCOULR2Lg
        engPOTENTLg = eGeoDefLg + fEngWF(5)*engVDWR1Lg + engCOULLg +
     &                molSolEnLg*fEngWF(9) + hBHxYeng128Lg
c 
          write(kanalPDBfin,'(a8)')'$ENGLIG: '
          write(kanalPDBfin,'(a12,f14.5)')'eVbondDefLG:',eVbondDefLg
          write(kanalPDBfin,'(a12,f14.5)')'eVangDefLG :',eVangDefLg
          write(kanalPDBfin,'(a12,f14.5)')'eImpDefLG  :',eImpDefLg
          write(kanalPDBfin,'(a12,f14.5)')'eTorsDefLG :',eTorsDefLg
          write(kanalPDBfin,'(a12,f14.5)')'engVDWR1LG :',engVDWR1Lg
          write(kanalPDBfin,'(a12,f14.5)')'engCOULR1LG:',engCOULR1Lg
          write(kanalPDBfin,'(a14,f14.5)')'hBHxYeng128LG:',hBHxYeng128Lg
          write(kanalPDBfin,'(a12,f14.5)')'engCOULR2LG:',engCOULR2Lg
          write(kanalPDBfin,'(a12,f14.5)')'restr1EngLG:',restr1EngLg
          write(kanalPDBfin,'(a12,f14.5)')'eRstHW1MLLG:',restr1MHWEngLg
          write(kanalPDBfin,'(a12,f14.5)')'eGeoDefLG  :',eGeoDefLg
          write(kanalPDBfin,'(a12,f14.5)')'engCOULLG  :',engCOULLg
          write(kanalPDBfin,'(a12,f14.5)')'engSolvLG  :',molSolEnLg
          write(kanalPDBfin,'(a12,f14.5)')'engPOTENTLG:',engPOTENTLg
          write(kanalPDBfin,'(a7)')'$ENDLIG '
          end if ! vLigFlag .eq. 1
c
          write(kanalPDBfin,'(a12)')'REMARK: PDB: '
c
           do ich=1,nChain
           do k=startAtInCha(ich),stopAtInCha(ich)
           k3 = 3*k-3
           write(kanalPDBfin,7071)"ATOM  ",
     &     atomNumb(k),atomName(k),resName(k),chName(k),resNumb(k),
     &     (atomXYZopt(k3+j),j=1,3), atomQ(k)
           end do !k
c TERline
           if( (resName(stopAtInCha(ich))(1:3) .ne. 'HOH' .and.
     &        resName(stopAtInCha(ich))(1:3) .ne. 'ION') .or.
     &        ich .eq. nChain)
     &     write(kanalPDBfin,7071)"TER   ", stopAtInCha(ich)      
           end do !ich 
c TERminal line
           write(kanalPDBfin,7071)"END   ", natom
           close(kanalPDBfin)
           end if !iw
c
        if(CONTROLx)then
        write(kanalp,*)'engOptFunctBS:END: atomXYZopt:'
        i3 = 3*nmoveatom
        call print9r(kanalp,i3,atomXYZopt)
        end if!Contr
c
        if(CONTROL)then
        write(kanalp,*)'engOptFunc Finish:'
        end if
c
         return
c
7071    format(a6,1x,i4,2x,a4,a4,a1,i4,4x,3f8.3,7x,f8.5) ! PDB rq
c
	 end
