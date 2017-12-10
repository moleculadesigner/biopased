c mdRUN     
c
c Yuri Vorobjev  2002
c                2004 : hB128
c                2004 : RigBody
c                2005 : GBorn solvation
c                2009 : ExWatSH
c PWat:distribute
c MTS with Verlet Velocity integrator in the inner(smallest dt) loop
c
c Yuri Vorobjev  2002                           
c
	subroutine mdRun(ntimeMX,ntime0,ntime,ntimeR1,ntimeR2,
     &            ntimeF1,ntimeF2,ntimeF3,deltat,
     &            tempTg1,tempTg2,tauTRF,atype,optra,wtra,nwtra,cltra)	
c 
c MD RUN makes MDtraj for files in mdAtomXYZvel.h [ atomXYZ0(*),atomVel0(*) ]
c        call initMDStart(T) inits the MD start from the atomXYZ(*)-->atom0XYZ(*)
c 
c ntimeMX max number of time steps to be done in the call mdRun
c ntime0 - executed number of timesteps im the previous call
c ntime   executed number of timesteps in the call 
c ntimeR1, ntimeR2 - update frequency  for R1, R2 pairLists
c ntimeF1,ntimeF2 - update freq for R1=(vdw+coulR1), R2-coulR2 en/forces
c ntimeF3 -  SOLVation forces
c GeoEn/force ntimeF1=1 - standart
c deltat- Fastes timestep, 
c tempTg1 - initial(temp) of MD run
c tempTg1 -->tempTg2 
c run Md with slow change of the temperature from tempTg1 --> tempTg2
c           over the ntimeMX time steps for NTV ansemble[K]
c tauTRF - tau Relaxation Factor [ps]
c atype - ansamble type = 0/1 - NEV, NTV
c wtra - 0/1 - flag to write traj
c nwtra - write each nwtra snap
c optra - 0/1 - flag to open traj files
c cltra - 0/1 - flag to close traj files
c xyzTraFile, engTraFile - names of traj files
c
        include 'xyzPDBsize.h'
c
        include "output.h"
        include "kanalUse.h"
        include "statusMessg_mDyn.h"
c
        include 'xyzPDBinfo.h'
        include 'xyzPDBcrd.h'    !edyv
        include 'movingAtom.h'
        include 'enForce.h'
        include 'filedat.h'  
        include 'mdAtomXYZvel.h'
        include 'optionPar.h'
        include 'pair1234array.h'
        include 'nbondPairVCS.h'
        include 'solvGSarray.h'
        include 'engXYZtra.h'
        include 'ligInfo.h'
        include 'compGeoEF.h'
        include 'hbond128.h'
        include 'coulEnPar.h'
	include 'enForce_PW.h'
	include "solvate01.h"
	include "restrain1RHW.h"
	include "exWatShRestrSASPot.h"
c
        integer ntimeMX,ntime,ntime0
        integer ntimeR1,ntimeR2
        integer ntimeF1,ntimeF2,ntimeF3
        real deltat
        real tempTg1,tempTg2
        real tauTRF
        integer atype,optra,wtra,nwtra,cltra
c  local
        real atomXYZ0wr(3*natomMAX)         !kabsh rotated XYZ to write traj
        real sccal,sct,sc
        real lambTR(3),lambTRpw
	real tempTg
        real deltat2,deltat3
        real deltat05
        real invSt,an1,an2
cPWEx:
        integer Nfd(3),kPW
        integer ntimew
	integer ctype
c
        character*4 nRecPdbC4,nRecPdbC
        character*(charLenMAX) filePdbN,filePdb 
        character*(charLenMAX) fileFinPdbN,fileFinPdb
        integer nLfilePdb0,nLfilePdb
        integer nLfileFinPdb
        logical OPT_wrtieMdFinPDB
        integer i21,i1,i2,i,it,k,k1,k2,k3
        integer ia,iam,ia3,i3,na3,j,ich,sH
        integer kf12,kf23
        integer kf1,kf2,kf3,kf3MX
        integer kDrLaZnew_OPT
        integer kEFSsolvUpdate_OPT
        integer kSolvGBornNew_OPT
	integer kHPsolvUpdate_OPT
	integer kSolvateExWatRefPos_OPT
	real scaleCoulPW_OPT
	real scaleEngSOLVSheLL
	real*8 uRotMx(3,3),CMtmp(3),CMin(3)
        integer doOmega
        integer kanalp
        logical doWRpdb
        logical CONTROL
        logical CONTROL1,CONTROL2
        logical CONTROLFHb128
c yv:040109
	logical OPT_writeXYZForce  !yv:040109
	real atomForceWrite(3)     !yv:040109
c
        kanalp =  kanalRunOut
        CONTROL = .true.     
        CONTROL1 = .false.
        CONTROL2 = .false.
        CONTROLFHb128 = .false.
c
c        OPT_writeXYZForce= .true.
        OPT_writeXYZForce= .false.
        OPT_wrtieMdFinPDB = .true. 
c        OPT_WRxyzvtra = .true.
c
        kSolvateExWatRefPos_OPT = 10000! 1000  !updateFreq for ExWatRefPos
        kDrLaZnew_OPT = 5    ! update frequency for DrSASLaZ model
	kHPsolvUpdate_OPT = 5  !update frequency for HPsolv force
        kSolvGBornNew_OPT = 10  ! update frequency for bornRad(*)
        kEFSsolvUpdate_OPT = 10 ! update frequency for WatBrg EFS sas03 solvation
cx	scaleCoulPW_OPT = 1.0
	scaleCoulPW_OPT = 0.5  ! approximate Epolz ~ 0.5Ecoul(P-W)
	scaleEngSOLVSheLL=1.0 ! scaling for thin shell
c
        sccal = 2.390*0.001*0.5      ! scaling (MV^2) -> kcal 
        sct = 1.2028
c        tauTRF = 0.5     ! ps gromos
cx PWatEx
        lambTR(1) = 1.0     ! default scaling
	lambTR(2) = 1.0
	lambTR(3) = 1.0
        Nfd(1) = nAtFreeDg(1) ! number of freedom degrees
        Nfd(2) = nAtFreeDg(2)
	Nfd(3) = nAtFreeDg(3)
c
        if(CONTROL)then
        write(kanalp,*)'mdRun start: ntimeMX:',ntimeMX
        write(kanalp,*) 'EnergyTerms:fEngWF(*):',fEngWF
        write(kanalp,*) 'N atom freedom degrees:',Nfd
        end if
c
        makeVdW = 0
        makeCL = 0
        makeSL = 0
c
          engVDWR1TraSd = 0.0
          engCOULTraSd = 0.0 
          engPOTENTTraSd = 0.0                   
          tempT0traSd = 0.0
          engVDWR1TraAv = 0.0
          engCOULTraAv = 0.0
          engMolSolAv = 0.0
          engMolSolSd = 0.0 
          engPOTENTTraAv = 0.0
          tempT0traAv = 0.0
          hBHxYeng128TraAv = 0.0
          hBHxYeng128TraSd = 0.0
c        
        nLfilePdb0 = charLenMAX   ! n letters in file name
        nRecPdbC = '0000'
        tempTg = tempTg1        
c
c open Traj files in currentDir ./
c
        if( optra .eq. 1) then

	if(CONTROL)then
        write(kanalp,*)'XYZ tra file:',xyzTraFile
        write(kanalp,*)'Eng tra file:',engTraFile  
        write(kanalp,*)'PDB tra file:',pdbResFile
        write(kanalp,*)'PDB fin file:',pdbFinFile
        end if
c
cx        kanalXYZtra = 31
cx        kanalEngtra = 32
c
        if(OPT_WRxyzvtra)then
c
        open(unit=kanalXYZtra, file=xyzTraFile, form='formatted',
     &       status='unknown')
c
        end if ! OPT_WRxyzvtra
c
        open(unit=kanalEngtra, file=engTraFile, form='formatted',
     &       status='unknown')
c
        end if ! optra
c 
        call clean_spacelf(pdbResFile,nLfilePdb0,filePdb,nLfilePdb)
        call clean_spacelf
     &       (pdbFinFile,nLfilePdb0,fileFinPdb,nLfileFinPdb)
c
          if(wtra .eq. 1 )then
          write(kanalEngtra,'(a8)')'$ENERGY:'
          if(atype .eq. 0)then
          write(kanalEngtra,'(a30,f8.2)')
     &          '#NEV ansamble TempInitial(K): ',tempTg
          end if
          if(atype .eq. 1)then
          write(kanalEngtra,'(a30,f8.2)')
     &          '#NTV ansamble TempTarget(K):  ',tempTg
          else 
          write(kanalEngtra,'(a30,f8.2)')
     &          '#ENV ansamble TempRInit(K):  ',tempTg
          end if
          end if
c 
c adjust frequencies to update F1,F2,F3 forces (fast,medium,slow)
        kf12 = int(ntimeF2/ntimeF1 + 0.5)
        if(kf12 .lt. 1) kf12=1
        ntimeF2 = kf12*ntimeF1
c
        kf23 = int(ntimeF3/ntimeF2 + 0.5)
        if(kf23 .lt. 1) kf23=1 
        ntimeF3 = kf23*ntimeF2
c MAX Numb of largest timeSteps
        kf3MX = int(ntimeMX/(kf12*kf23) + 0.5)
c
c adjust frequencies to update pairList R1 and R2
c ntimeR1 and ntimeR2 should have kf12 and kf23 as a multiplyers
         i21 = int(ntimeR1/(kf12*kf23) + 0.5)
         if(i21 .lt. 1) i21 = 1
         ntimeR1 = i21*kf12*kf23
         i21 = int(ntimeR2/ntimeR1 +0.5) 
         if(i21 .lt. 1) i21 = 1 
         ntimeR2 = i21*ntimeR1
c
         write(kanalp,*)
     &   'mdRun: Adjusted Update freQ:PLR1,PLR2:',ntimeR1,ntimeR2
c
c  adjust frequencies to write traj file
        k = int(nwtra/(kf12*kf23) + 0.5)
        if(k .lt. 1) k=1
        nwtra = kf12*kf23*k
c
        write(kanalp,*)'mdRun: traj WR freQ:',nwtra
c
c long loop over time
c        fcall = 0        !initMDStart
        deltat05 = deltat*0.5
        deltat2 = kf12*deltat
        deltat3 = kf23*deltat2
c
        write(kanalp,'(a29,3(f6.4,1x))')
     &  'mdRun:MTS: dt1,dt2,dt3(ps): ',deltat,deltat2,deltat3
c
        it = 0
        ntime = it 
        nStepTra = 0   ! number of TRAinfo written in this mdRun call
c
c largest timeStep loop: dt3 = kf12*kf23*dt
        do kf3 = 0,kf3MX
c define current bathT(time): tempTg (t)
        tempTg = tempTg1 + (tempTg2-tempTg1)*kf3/kf3MX
c
c update SLOW force F3
c atomXYZ dependent Dconst   
        if(coulVar_OPT .eq. 3)then      
        if( (kf3/kDrLaZnew_OPT)*kDrLaZnew_OPT .eq. kf3)then !update
        call initSASdielModLaZ(atomXYZ0)
        end if !(kf3/kDrLaZnew_OPT
        end if !coulVar_OPT=3
c update exWatSheLLRefPos
cx        if(OPT_SolvateExWat .and. 
cx     &  (((kf3+1)/kSolvateExWatRefPos_OPT)*kSolvateExWatRefPos_OPT
cx     &                        .eq.(kf3+1)))then !update
cx        call initRestr1RHWAtList(atomXYZ0,nRestr1RHWSeg,resEndRestr1RHW,
cx     &        nRestr1ResidHW,restr1ResidHWList,
cx     &        nRestr1RHWAtom,restr1RHWatomList,
cx     &        restr1RefResidPosHW,restr1ResidMass)
cx        end if!OPT_SolvateExWat update
c
c update solvent forces/energy F3
         if (OPT_SolvGS) then
c
         call SolventEnForces(natom,atomXYZ0,
     &         atomName,startPairL12,nPairL12,pair12List,
     &         nbpairListS,startnbPairLS,nnbPairLS,
     &         atomSolPar,molSolEn,atomSolEn,atomSolFr)
c
         if(CONTROL2)then
         write(kanalp,*)'zeroAllAtomForce: atomSolF'
         call zeroAllAtomForce(nmoveatom,moveAtomList,atomSolFr) 
         end if !C2
         end if !SolvGS
c
         if(OPT_SolvGBorn)then
         if( (kf3/kSolvGBornNew_OPT)*kSolvGBornNew_OPT .eq. kf3)then !update 
         call initBornRadSAS(atomXYZ0) 
         end if ! update borRad
c
         call getBornPolzEforce(atomXYZ0,
     &               bornPolzEng,bornPolzEngAt,bornPolzForceAt)
c
         end if ! OPT_SolvGBorn
c
         if(OPT_SolvEFS)then
         if( (kf3/kEFSsolvUpdate_OPT)*kEFSsolvUpdate_OPT .eq. kf3 )then !update EFS WaterBRg
         call solvateMoL03(atomXYZ0)
         end if ! upDateEFSwatBrg
c
         call getWatBrgForceSAS03(atomXYZ0,
     &                   watBrgEnergySoL,watBrgForceAt)
c
        end if !OPT_SolvEFS
c
         if(OPT_SolvSASHP)then
	 if( (kf3/kHPsolvUpdate_OPT)*kHPsolvUpdate_OPT .eq. kf3 )then !update 
c update SAS04 :
          call getHPolvForceSAS03(atomXYZ0,
     &                   hpSoLEnergy,hpSoLForceAt)
c update exWatSASatNeigh List:
       if(OPT_SolvateExWat)then 
       call initExWatShRestrSASPot
     &     (ncall_exWSASPot,natomWSolv,atomXYZ0(3*natomSolutMoL+1),
     &                                  eng_exWatSASP,ff_exWatSASP)
       end if!OPT_SolvateExWat
         end if!update
         end if !OPT_SolSASHP
c
c artificial Compactization forces
c
         if(OPT_CompactForce)then
         call compactGeoCntEF(atomXYZ0,
     &           natom,nmoveatom,moveAtomList,
     &           compactGeoFpar,compactGeoEn,compactGeoForce)
         end if ! OPT_CompactForce
c
         if(OPT_SolvateExWat)then
	 call get_expLSASPot02(natomWSolv,atomXYZ0(3*natomSolutMoL+1),
     &           eng_exWatSASP,ff_exWatSASP)
         end if!OPT_SolvateExWat
c
c COLLECT F3 forces
c
         if(.not.OPT_SolvateExWat)then
         solvMolEngF3 = molSolEn+hpSoLEnergy+watBrgEnergySoL
     &                  + bornPolzEng
         end if!.not.OPT_SolvateExWat
	 if(OPT_SolvateExWat)then
c collect true solvationEng of Protein= allEng(P-W:2)
         solvMolEngF3 = engVDWR1_PWat(2)+engCOULR1_PWat(2)
     &                 + engCOULR2_PWat(2)+ molSolEn_PWat(1)
c
	 end if!OPT_SolvateExWat
c
         do i = 1,3*natom
         atomForceF3(i) = atomSolFr(i) + compactGeoForce(i)
     &    + watBrgForceAt(i) + hpSoLForceAt(i) 
     &    + bornPolzForceAt(i) + ff_exWatSASP(i)
         end do !i
c correct to ZEROsum
cx         call zeroAllAtomForce(nmoveatom,moveAtomList,atomForceF3)
         ctype = 1
	 call corr_force_to_zeroSum(ctype,natom,atomForceF3)
c
	do kf2 = 1,kf23
c
c update forces/energy F2
        call allAtVDWEForceR2(atomXYZ0,atomQ,
     &           natom,nmoveatom,moveAtomList,
     &           nbpairListC,startnbPairLC,nnbPairLC,
     &           rcutV,rcutC,engCOULR2,coulForceR2)
c COLLECT F2 forces
        do i = 1,3*natom
        atomForceF2(i) = coulForceR2(i)
        end do !i
c
        if(CONTROL2)then
c correct to ZERO summ of allForces
        write(kanalp,*)'zeroAllAtomForce: coulForceR2'
        call zeroAllAtomForce(nmoveatom,moveAtomList,coulForceR2)
        end if !C2
c
        do kf1 = 1,kf12
c fast force F1 loop
c current timeStep(fast) *dt
c it = number of time propogation done      
        ntimew = it + ntime0
c
        i1 = ntime/ntimeR1
        i2 = ntime/ntimeR2
c
c define flags to update nonBondPairList R1,R2,SOLV
        makeVdW = 0
        makeCL = 0
        makeSL = 0
        if(i1*ntimeR1 .eq. ntime )then
        makeVdW = 1
        if(i2*ntimeR2 .eq. ntime )then
        makeCL = 1
        makeSL = 1
        end if
        end if 
c
        if(makeVdW .eq. 1 .or. makeCL .eq. 1 ) then
        write(kanalp,*)
     &  '* * * * * * * * * * * * * * * * * * * * * * * * * * * * *'
        write(kanalp,*)
     &  'mdRun: tstep:', ntime,
     &  '  Update NonBondPairList:makeVdW,makeCL,makeSL:',
     &  makeVdW,makeCL,makeSL
        write(kanalp,*)
     &  '* * * * * * * * * * * * * * * * * * * * * * * * * * * * *'
c
        call initNonBondList(atomXYZ0,makeVdW,makeCL,makeSL) 
c
c update hB128 tripletList
        if(OPT_HBond128)then
        call getHbTripletAll(atomXYZ0,
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
        if(makeVdW .eq. 1)makeVdW = 0
        if(makeCL .eq. 1)then
        makeCL = 0
        makeSL = 0
        end if
c
        end if !makeVdW .eq. 1 .or. makeCL
c
c update F1 forces/energy
	call initAllForce01(fcall,fEngWF,atomXYZ0,
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
        if(CONTROL2)then
        write(kanalp,*)'zeroAllAtomForce: vbdefForce:'
        call zeroAllAtomForce(nmoveatom,moveAtomList,vbdefForce)
c 
        write(kanalp,*)'zeroAllAtomForce: vAngdefForce:'
        call zeroAllAtomForce(nmoveatom,moveAtomList,vAngdefForce)
c
         write(kanalp,*)'zeroAllAtomForce: impDefForce:'
         call zeroAllAtomForce(nmoveatom,moveAtomList,impDefForce)
c
         write(kanalp,*)'zeroAllAtomForce: torsAngForce:'
         call zeroAllAtomForce(nmoveatom,moveAtomList,torsAngForce)
c
         write(kanalp,*)'zeroAllAtomForce: vdwForceR1:'
         call zeroAllAtomForce(nmoveatom,moveAtomList,vdwForceR1)
c
         write(kanalp,*)'zeroAllAtomForce: coulForceR1:'
         call zeroAllAtomForce(nmoveatom,moveAtomList,coulForceR1)
c
         write(kanalp,*)'zeroAllAtomForce: hBHxY128force:'
         call zeroAllAtomForce(nmoveatom,moveAtomList,hBHxY128force)
c
        write(kanalp,*)'zeroAllAtomForce: atomForceF1:'
        call zeroAllAtomForce(nmoveatom,moveAtomList,atomForceF1)
        end if !CONTROL2
c
c COLLECT TOTAL Forces F1+F2+F3
        do i=1,3*natom
        atForceTot(i) = atomForceF1(i) + atomForceF2(i) + atomForceF3(i)
        end do! i 
c
c update RigidBody parameters
         if(OPT_RigidBodyMd)then
c
         doOmega = 0
         call getRigBodyParam01(doOmega,atomXYZ0,atomVel0)
c
         call getRigBRotTraForce(atForceTot)   
c
         call getRigBdOmegaDt(atomXYZ0,atForceTot) 
c
         call getRigBodyRotXV(atomXYZ0,deltat)  !testRig
c
         end if !OPT_RigidBodyMd
c
c make timeStep for transl with totForces=F1 + F2 + F3 and make Shake
c
	call mdTimeStepProp01(nmoveatom,moveAtomList,
     &                                moveFlag,deltat)
c
        it = it + 1
        ntime = it
        ntimew = ntime + ntime0 
c 
        end do ! kf1 - small timeStep
c
        end do ! kf2 - medium timeStep
c
c end of medium timeStep
c correct CM velocity
        if(OPT_fullProtMD)then
        call initZeroCMVelocity(natom,atomMass,atomVelm)
        end if ! correc CM Vel OPT_fullProtMD 
c sum up final Energies:
c
          eGeoDef   = fEngWF(1)*eVbondDef + fEngWF(2)*eVangDef 
     &                + fEngWF(3)*eImpDef + fEngWF(4)*eTorsDef
          if(OPT_SolvateExWat)then
	  eGeoDef_PWat(1) = fEngWF(1)*eVbondDef_PWat(1)+
     &         fEngWF(2)*eVangDef_PWat(1)+fEngWF(3)*eImpDef +
     &         fEngWF(4)*eTorsDef
          eGeoDef_PWat(2) = fEngWF(1)*eVbondDef_PWat(2)+
     &         fEngWF(2)*eVangDef_PWat(2)
          eGeoDef_PWat(3) = fEngWF(1)*eVbondDef_PWat(3)+
     &         fEngWF(2)*eVangDef_PWat(3)
c
          eVbondDef = eVbondDef_PWat(1)
	  eVangDef = eVangDef_PWat(1)
          eGeoDef   = eGeoDef_PWat(1)
	  end if!OPT_SolvateExWat
c
          eGeoRst = fEngWF(8)*(restr1Eng + restr1MHWEng + eng_exWatSASP
     &                + restrDistA2Eng) + compactGeoEn*fEngWF(9)
c
          engCOUL   = fEngWF(6)*engCOULR1 + fEngWF(7)*engCOULR2
c
c true ePOtential of Internal interactions ! no Artificial RESTRAInts ENergy = eGeoRst
c
          if(.not.OPT_SolvateExWat)
     &    solvMolEngF3 = (molSolEn + watBrgEnergySoL+ hpSoLEnergy +
     &                   bornPolzEng)*fEngWF(9)
          if(OPT_SolvateExWat)then
c collect true solvationEng of Protein= allEng(P-W:2)
c scaleCoulPW-OPT=0.5 to approximate engPOLz by 0.5*ELECT(P-W): linearResponce:
c
          solvMolEngF3 = (engVDWR1_PWat(2)+molSolEn_PWat(1)+
     &	  scaleCoulPW_OPT*(engCOULR1_PWat(2)+ engCOULR2_PWat(2)+
     &                  hBHxYeng128_PWat(2)))*fEngWF(9) 
          end if!OPT_SolvateExWat
c
          if(.not.OPT_SolvateExWat)
     &    engPOTENT = eGeoDef + fEngWF(5)*engVDWR1 + engCOUL + 
     &        + hBHxYeng128
     &        + solvMolEngF3
c
          if(OPT_SolvateExWat)then
c
        scaleEngSOLVSheLL=1.00
	call scalePWatShEng(dSOLVSheLL,scaleEngSOLVSheLL)
c scale (P-Wat) and add hpSoLEnergy[soluteMol]
	solvMolEngF3 = solvMolEngF3*scaleEngSOLVSheLL
     &                 + hpSoLEnergy
c
	  do k=1,3
	  engCOUL_PWat(k)=engCOULR1_PWat(k)+engCOULR2_PWat(k)
	  engPOTENT_PWat(k)=eGeoDef_PWat(k)+ engVDWR1_PWat(k)+
     &        engCOUL_PWat(k) + hBHxYeng128_PWat(k)
          end do!k
	  engPOTENT_PWat(1)=engPOTENT_PWat(1)+solvMolEngF3
c PW: redefine
	  engPOTENT = engPOTENT_PWat(1)
	  engVDWR1=engVDWR1_PWat(1)
	  engCOULR1=engCOULR1_PWat(1)
	  engCOULR2=engCOULR2_PWat(1)
	  engCOUL = engCOUL_PWat(1)         
	  hBHxYeng128 = hBHxYeng128_PWat(1)
	  eGeoDef = eGeoDef_PWat(1)
	  molSolEn = molSolEn_PWat(1)
	  end if!OPT_SolvateExWat
c
c calculate kinetic energy
        kinEng(1) = 0.0
	kinEng(2) = 0.0
	kinEng(3) = 0.0
        do iam = 1,nmoveatom
        ia = moveAtomList(iam)
        i3 = ia*3-3
        do k = 1,3
        atomVel0(i3+k)=(atomVelm(i3+k)+atomVelp(i3+k))*0.5
	if(atomPWatTag(ia).eq.1)kPW=2
	if(atomPWatTag(ia).eq.2)kPW=3
        kinEng(kPW) = kinEng(kPW) + atomMass(ia)*atomVel0(i3+k)**2 
        end do !k
        end do!iam
	kinEng(1) = kinEng(2)+kinEng(3)
c
c PWatEx: separate Temp0 of Protein and solvent
        do k=1,3
	tempT0(k) = 0.0
	lambTR(k) = 1.0
	if(Nfd(k) .ne. 0)
     &        tempT0(k) = kinEng(k)*sct/Nfd(k)      ! kineticTemperature =1(total), =2(Prot),=3(watSolv)
        kinEng(k) = kinEng(k)*sccal                 ! *0.0005! energy in Kcal/mol
	if(tempT0(k) .gt. 0.0)
     &	lambTR(k) = sqrt(1.0 + (deltat3/tauTRF)*
     &	            (tempTg/tempT0(k) -1.0))        !separate scaling for ProtAtom and Water
	end do!k
c
        if(CONTROL1)then
        write(kanalp,*)'mdRun: ia, mass, vel:'
        do iam = 1,nmoveatom
        ia = moveAtomList(iam)
        i3 = 3*ia-3
        write(kanalp,*)ia,atomMass(ia),
     &                 atomVel0(i3+1),atomVel0(i3+2),atomVel0(i3+3)
        end do !iam
        end if !Control1
c
        engTotal = engPOTENT + kinEng(1)
c
	if(atype .eq. 1)then
c scale velocity for NTV ansemble: Berendsen termostat
c PWatEx:
cx        lambTR = sqrt(1.0 + (deltat3/tauTRF)*(tempTg/tempT0 -1.0))
c
        lambTRpw=lambTR(1)
        do iam = 1,nmoveatom
        ia = moveAtomList(iam)
        i3 = ia*3-3
	if(atomPWatTag(ia).eq.1)lambTRpw=lambTR(2)
	if(atomPWatTag(ia).eq.2)lambTRpw=lambTR(3)
        do k=1,3
        atomVelp(i3+k) = lambTRpw*atomVelp(i3+k)
        atomVelm(i3+k) = atomVelp(i3+k) 
        end do!k
        end do!iam 
        end if !atype
c
c write Traj
c          write(kanalp,*)'mdStep Done: dt3 loop: ntime=',ntime,
c     &    ' nwtra=',nwtra
c
c Control print
          if(CONTROL)then
          write(kanalp,'(a57,i7,3f8.1,1x,5f8.1,2f7.1)')
     &  'mdRunE:nts,Tmd(K),ePot,eDef,eVDW,eCL,eSol,eRstA1,eRstM1: ',
     &   ntime,tempT0,engPOTENT,eGeoDef,engVDWR1,
     &   engCOUL,solvMolEngF3,restr1Eng,restr1MHWEng
          end if !C
c
          if (wtra .ge. 1) then
          k = ntime/nwtra
          if (ntime .eq. k*nwtra)then  !writeTraj
c
c Control print
        if(CONTROL1)then
        write(kanalp,'(a57,i7,3f8.1,1x,5f8.1,2f7.1)')
     &  'mdRunE:nts,Tmd(K),ePot,eDef,eVDW,eCL,eSol,eRstA1,eRstM1: ',
     &   ntime,tempT0,engPOTENT,eGeoDef,engVDWR1,
     &   engCOUL,solvMolEngF3,restr1Eng,restr1MHWEng
        end if !C
c
c do write Traj
c
          nRecPdb = nRecPdb + 1                    ! global
c
c WRITE XYZ tra
          if(OPT_WRxyzvtra)then
          write(kanalXYZtra,'(a7)')'$XYZTRA'
c
          write(kanalXYZtra,'(a7,i6,2x,f8.3,a5)')
     &         '$nstep: ', ntimew,ntimew*deltat,' (ps)'
c
          write(kanalXYZtra,'(a45,a25)')
     & '#   x      y        z       x       y       z',
     & '        x       y       z'
c
          na3 = 3*natom
          call print9r(kanalXYZtra,na3,atomXYZ0)
          write(kanalXYZtra,'(a4)')'END '
c
          write(kanalXYZTra,'(a7)')'$VELTRA'
c
          write(kanalXYZtra,'(a7,i4,f8.3,a5)')
     &         '#nstep: ', ntime,ntime*deltat,' (ps)'
c
          write(kanalXYZtra,'(a45,a25)')
     & '#   x      y        z       x       y       z',
     & '        x       y       z'
c
          call print9r(kanalXYZtra,na3,atomVelm)
          write(kanalXYZtra,'(a4)')'END '
c
          end if ! OPT_WRxyzvtra
c
c write energyTraj
          write(kanalEngtra,'(a8)')'$ENERGY: '
          write(kanalEngtra,'(a7,i5,f8.3,a5,a9,f8.1)')
     &    '#nstep: ', ntimew,ntimew*deltat,' (ps)',' tgT(K):',tempTg
          write(kanalEngtra,'(a10,f14.5)')'eVbondDef:',eVbondDef
          write(kanalEngtra,'(a10,f14.5)')'eVangDef :',eVangDef
          write(kanalEngtra,'(a10,f14.5)')'eImpDef  :',eImpDef
          write(kanalEngtra,'(a10,f14.5)')'eTorsDef :',eTorsDef
          write(kanalEngtra,'(a10,f14.5)')'engVDWR1 :',engVDWR1
          write(kanalEngtra,'(a10,f14.5)')'engHb128 :',hBHxYeng128 
          write(kanalEngtra,'(a10,f14.5)')'engCOULR1:',engCOULR1
          write(kanalEngtra,'(a10,f14.5)')'engCOULR2:',engCOULR2
          write(kanalEngtra,'(a10,f14.5)')'eRstH1At :',restr1Eng
          write(kanalEngtra,'(a10,f14.5)')'eRstHW1RS:',
     &	  restr1RHWEng+eng_exWatSASP
          write(kanalEngtra,'(a10,f14.5)')'eRstHW1ML:',restr1MHWEng
          write(kanalEngtra,'(a10,f14.5)')'eRstH2At :',restrDistA2Eng 
          write(kanalEngtra,'(a10,f14.5)')'eCompcEng:',compactGeoEn  
          write(kanalEngtra,'(a10,f14.5)')'eGeoRst  :',eGeoRst
          write(kanalEngtra,'(a10,f14.5)')'eGeoDef  :',eGeoDef
          write(kanalEngtra,'(a10,f14.5)')'eCOUL    :',engCOUL
          write(kanalEngtra,'(a10,f14.5)')'eSolvGS  :',molSolEn
          write(kanalEngtra,'(a10,f14.5)')'eSolvCav :',hpSoLEnergy
          write(kanalEngtra,'(a10,f14.5)')'eSolvWbr :',watBrgEnergySoL
          write(kanalEngtra,'(a10,f14.5)')'eSolvGBp :',bornPolzEng
          write(kanalEngtra,'(a10,f14.5)')'eSolvTot :',solvMolEngF3 
          write(kanalEngtra,'(a10,f14.5)')'ePOTENT  :',engPOTENT
          write(kanalEngtra,'(a10,3f12.3)')'kinEng   :',kinEng
          write(kanalEngtra,'(a10,3f12.3)')'Temp(K)  :',tempT0
          write(kanalEngTra,'(a10,3f12.3)')'lambTR   :',lambTR
          write(kanalEngtra,'(a10,f14.5)')'engTotal :',engTotal
          write(kanalEngtra,'(a4)')'$END '
c
c - - - - - Ligand - Prot - - - - - - - - - - - - - - - - - - - - - - 
c write energy:Lig-Prot
          vLigFlag = 0
          if(OPT_LigRes)vLigFlag = 1
          if(vLigFlag .eq. 1) then
          makeVdWLig=1
          makeCLLig =1
          makeSLLig =1
c
          call 
     &    initNonBondListLig(atomXYZ0,makeVdWLig,makeCLLig,makeSLLig)
c
          call initAllForceLig(fcall,atomXYZ0,
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
        eGeoRstLg   = fEngWF(8)*(restr1EngLg+restr1MHWEngLg)
c
        engCOULLg   = fEngWF(6)*engCOULR1Lg + fEngWF(7)*engCOULR2Lg
c
        engPOTENTLg = eGeoDefLg + fEngWF(5)*engVDWR1Lg + engCOULLg +
     &                molSolEnLg*fEngWF(9)+hBHxYeng128Lg
c
          write(kanalEngtra,'(a8)')'########'
          write(kanalEngtra,'(a8)')'$ENGLIG: '
          write(kanalEngtra,'(a12,f14.5)')'eVbondDefLG:',eVbondDefLg
          write(kanalEngtra,'(a12,f14.5)')'eVangDefLG :',eVangDefLg
          write(kanalEngtra,'(a12,f14.5)')'eImpDefLG  :',eImpDefLg
          write(kanalEngtra,'(a12,f14.5)')'eTorsDefLG :',eTorsDefLg
          write(kanalEngtra,'(a12,f14.5)')'engVDWR1LG :',engVDWR1Lg
          write(kanalEngtra,'(a14,f14.5)')'hBHxYeng128Lg:',hBHxYeng128Lg
          write(kanalEngtra,'(a12,f14.5)')'engCOULR1LG:',engCOULR1Lg
          write(kanalEngtra,'(a12,f14.5)')'engCOULR2LG:',engCOULR2Lg
          write(kanalEngtra,'(a12,f14.5)')'restr1EngLG:',restr1EngLg
          write(kanalEngtra,'(a12,f14.5)')'eRstHW1MLLG:',restr1MHWEngLg
          write(kanalEngtra,'(a12,f14.5)')'eGeoDefLG  :',eGeoDefLg
          write(kanalEngtra,'(a12,f14.5)')'engCOULLG  :',engCOULLg
          write(kanalEngtra,'(a12,f14.5)')'engSolvLG  :',molSolEnLg
          write(kanalEngtra,'(a12,f14.5)')'engPOTENTLG:',engPOTENTLg
          write(kanalEngtra,'(a7)')'$ENDLIG '
c          
          end if ! vLigFlag .eq. 1
c
c - - - END - - Lig-Prot - - - - - - - - - - - - - - - - - - - - - - - -
c
	  nStepTra = nStepTra + 1
          if(nStepTra .gt. nstepTraMAX) then
	  write(kanalp,*)
     &    'ERROR: mdRun: nstepTraMAX is SMALL (engXYZtra.h)!'
          stop
          end if !nStepTra .gt. nstepTraMAX
c
c calculate RMSD for current atomXYZtmp()
           if(.not.OPT_SolvateExWat)then
           call initSwap01(natom,atomXYZ0,atomXYZtmp)
	   call initSwap01(natom,atomXYZ0,atomXYZ0wr)
c
           call  kabsch8(atomXYZtmp,atomXYZin,
     &                    natom,atomMass,uRotMx,CMtmp,CMin)
           call xyz_rmsd(atomXYZin,atomXYZtmp,natom,rmsdTraj(nStepTra))
	   end if!not.OPT_SolvateExWat
c
           if(OPT_SolvateExWat)then
c RMSD for SoluteMolec
           call initSwap01(natom,atomXYZ0,atomXYZtmp)
	   call initSwap01(natom,atomXYZ0,atomXYZ0wr)
c
           call  kabsch8(atomXYZtmp,atomXYZin,
     &	                  natomSolutMol,atomMass,uRotMx,CMtmp,CMin)
c rotate WatSolvMolec by the uRotMx:
           call rotateMOLEC8(atomXYZtmp(3*natomSolutMol+1), 
     &	      (natom-natomSolutMol), uRotMx,CMtmp,CMin)
c
	   call xyz_rmsd(atomXYZin,atomXYZtmp,
     &	               natomSolutMol,rmsdTraj(nStepTra))
c
	   end if!OPT_SolvateExWat 
c REturn rotated XYZ to write traj:
           if(.not.OPT_SolvWbr)
     &       call initSwap01(natom,atomXYZtmp,atomXYZ0wr)
c
c NOTE: watBrg XYZ is not rotated! = atomXYZ0(*)
c geoControl
cx           call protBBondGeoAnalys(atomXYZ0)
cx
c write Md result PDB
c define fileName:       
         call conv_numb_to_char(nRecPdb,nRecPdbC)
c
         if(nRecPdb .le. 9) then
         nRecPdbC4 = '000'//nRecPdbC(1:1)
          else 
           if(nRecPdb .le. 99 )then
            nRecPdbC4 = '00'//nRecPdbC(1:2)
             else
              if(nRecPdb .le. 999 )then
               nRecPdbC4 = '0'//nRecPdbC(1:3)
               else 
                if(nRecPdb .le. 9999 )then
                 nRecPdbC4 = nRecPdbC(1:4)
                else 
                 write(kanalp,*)
     &           'mdRun: ERROR!:number of nRecPdb > 9999 '
                 stop
                end if
               end if
              end if
             end if
c
        filePdbN = filePdb(1:nLfilePdb-4)//nRecPdbC4//'.pdb'
c
c        write(kanalp,*)'mdRun: filePdbN:',filePdbN
c
cx        kanalPDBres = 33
        open(unit=kanalPDBres, file=filePdbN, form='formatted',
     &       status='unknown')
c
           write(kanalPDBres,'(a31,f10.4 )')
     &     "REMARK: Md result : MdTime(ps):", ntimew*deltat
          write(kanalPDBres,'(a9,i7)')'$nstep:  ', ntimew
          write(kanalPDBres,'(2(a9,i7))')
     &	  '$nRecPDB:', nRecPdb,' nStepTra:',nStepTra
c
          write(kanalPdbres,'(a18,f8.2)')
     &            'REMARK: RMSD(x0):  ',rmsdTraj(nStepTra)
c
          write(kanalPDBres,'(a8,2x,a6)')'$ENERGY: ',' :Kcal'
	  if(.not.OPT_SolvateExWat)
     &       write(kanalPDBres,'(a10,f14.5)')'eVbondDef:',eVbondDef
          if(OPT_SolvateExWat)
     &     call writePWatSh(kanalPDBres,"eVbondDef",eVbondDef_PWat)  
c
          if(.not.OPT_SolvateExWat)
     &    write(kanalPDBres,'(a10,f14.5)')'eVangDef :',eVangDef
	  if(OPT_SolvateExWat)
     &     call writePWatSh(kanalPDBres,"eVangDef ",eVangDef_PWat)
c
          write(kanalPDBres,'(a10,f14.5)')'eImpDef  :',eImpDef
          write(kanalPDBres,'(a10,f14.5)')'eTorsDef :',eTorsDef
c
          if(.not.OPT_SolvateExWat)
     &      write(kanalPDBres,'(a10,f14.5)')'engVDWR1 :',engVDWR1
          if(OPT_SolvateExWat)
     &      call writePWatSh(kanalPDBres,"engVDWR1 ",engVDWR1_PWat)
c
          if(.not.OPT_SolvateExWat)
     &    write(kanalPDBres,'(a10,f14.5)')'engHbHxY :',hBHxYeng128
          if(OPT_SolvateExWat)
     &      call writePWatSh(kanalPDBres,"engHbHxY ",hBHxYeng128_PWat)
c
          if(.not.OPT_SolvateExWat)
     &     write(kanalPDBres,'(a10,f14.5)')'eCOULR1  :',engCOULR1
          if(OPT_SolvateExWat)
     &    call writePWatSh(kanalPDBres,"eCOULR1  ",engCOULR1_PWat)
c
          if(.not.OPT_SolvateExWat)
     &    write(kanalPDBres,'(a10,f14.5)')'eCOULR2  :',engCOULR2
          if(OPT_SolvateExWat)
     &    call writePWatSh(kanalPDBres,"eCOULR2  ",engCOULR2_PWat)
c
          write(kanalPDBres,'(a10,f14.5)')'eRstHA1  :',restr1Eng
          write(kanalPDBres,'(a10,f14.5)')'eRstHW1RS:',restr1RHWEng
     &	  +eng_exWatSASP
          write(kanalPDBres,'(a10,f14.5)')'eRstHW1ML:',restr1MHWEng 
          write(kanalPDBres,'(a10,f14.5)')'eRstHA2  :',restrDistA2Eng 
          write(kanalPDBres,'(a10,f14.5)')'eCompcEng:',compactGeoEn  
c
          if(.not.OPT_SolvateExWat)
     &    write(kanalPDBres,'(a10,f14.5)')'eGeoDef  :',eGeoDef
          if(OPT_SolvateExWat)
     &    call writePWatSh(kanalPDBres,"eGeoDef  ",eGeoDef_PWat)
c
          if(.not.OPT_SolvateExWat)
     &    write(kanalPDBres,'(a10,f14.5)')'eCOULTot :',engCOUL
          if(OPT_SolvateExWat)
     &    call writePWatSh(kanalPDBres,"eCOULTot ",engCOUL_PWat)
c 
          if(.not.OPT_SolvateExWat)
     &    write(kanalPDBres,'(a10,f14.5)')'eSolvGS  :',molSolEn
          if(OPT_SolvateExWat)
     &    call writePWatSh(kanalPDBres,"eSolvGS  ",molSolEn_PWat)
c
          write(kanalPDBres,'(a10,f14.5)')'eSolvCav :',hpSoLEnergy
          write(kanalPDBres,'(a10,f14.5)')'eGBnPolz :',bornPolzEng
          write(kanalPDBres,'(a10,f14.5)')'eWatBridg:',watBrgEnergySoL
          write(kanalPDBres,'(a10,f14.5)')'eSolvTot :',solvMolEngF3
          write(kanalPDBres,'(a10,f14.5)')'engPOTENT:',engPOTENT
          if(OPT_SolvateExWat)
     &    call writePWatSh(kanalPDBres,"engPOTENT",engPOTENT_PWat)
c
          write(kanalPDBres,'(a10,3f12.5)')'kinEng   :',kinEng
          write(kanalPDBres,'(a10,3f12.5)')'Temp(K)  :',tempT0
          write(kanalPDBres,'(a10,f14.5)')'engTotal :',engTotal
c
          if(vLigFlag .eq. 1) then
          write(kanalPDBres,'(a8)')'########'
          write(kanalPDBres,'(a8)')'$ENELIG: '
          write(kanalPDBres,'(a12,f14.5)')'eVbondDefLG:',eVbondDefLg
          write(kanalPDBres,'(a12,f14.5)')'eVangDefLG :',eVangDefLg
          write(kanalPDBres,'(a12,f14.5)')'eImpDefLG  :',eImpDefLg
          write(kanalPDBres,'(a12,f14.5)')'eTorsDefLG :',eTorsDefLg
          write(kanalPDBres,'(a12,f14.5)')'engVDWR1LG :',engVDWR1Lg
          write(kanalPDBres,'(a14,f14.5)')'hBHxYeng128Lg:',hBHxYeng128Lg
          write(kanalPDBres,'(a12,f14.5)')'engCOULR1LG:',engCOULR1Lg
          write(kanalPDBres,'(a12,f14.5)')'engCOULR2LG:',engCOULR2Lg
          write(kanalPDBres,'(a12,f14.5)')'restr1EngLG:',restr1EngLg
          write(kanalPDBres,'(a12,f14.5)')'eRstHW1MLLG:',restr1MHWEngLg
          write(kanalPDBres,'(a12,f14.5)')'eGeoDefLG  :',eGeoDefLg
          write(kanalPDBres,'(a12,f14.5)')'engCOULLG  :',engCOULLg
          write(kanalPDBres,'(a12,f14.5)')'engSolvLG  :',molSolEnLg
          write(kanalPDBres,'(a12,f14.5)')'engPOTENTLG:',engPOTENTLg
          write(kanalPDBres,'(a7)')'$ENDLIG '
          end if !Lig
c
          write(kanalPDBres,'(a12)')'REMARK: PDB: '

           do ich=1,nChain
           do k = startAtInCha(ich),stopAtInCha(ich) 
           k3 = 3*k-3
c
           doWRpdb=.true.
           if(OPT_doLigDock .eq. 1)then
           doWRpdb=.false.
           if(resName(k) .eq. resName(atomInLigList(1))) doWRpdb=.true.
           end if !OPT_doLigDock .eq. 1
c
           if(doWRpdb)then
	   if(OPT_writeXYZForce)then
c write PDB+atForce: 
           do j=1,3
           atomForceWrite(j)= vdwForceR1(k3+j)*fEngWF(5) 
     &	   + coulForceR1(k3+j)*fEngWF(6) + hBHxY128force(k3+j) 
     &     + atomForceF2(k3+j)+atomForceF3(k3+j)
           end do!j
c
           write(kanalPDBres,7073)"ATOM  ",
     &     atomNumb(k),atomName(k),resName(k),chName(k),resNumb(k),
     &     (atomXYZ0wr(k3+j),j=1,3), 
     &     (atomForceWrite(j), j=1,3)
	   else
           write(kanalPDBres,7071)"ATOM  ",
     &     atomNumb(k),atomName(k),resName(k),chName(k),resNumb(k),
     &     (atomXYZ0wr(k3+j),j=1,3), atomQ0(k)
           end if !OPT_writeXYZForce
           end if !doWRpdb
           end do !k
c TERminal chainLine
cx           if( (resName(stopAtInCha(ich))(1:3) .ne. 'HOH' .and.
cx     &        resName(stopAtInCha(ich))(1:3) .ne. 'ION') .or.
cx     &        ich .eq. nChain) 
cx     &     write(kanalPDBres,7071)"TER   ", stopAtInCha(ich)

             write(kanalPDBres,7071)"TER   ", stopAtInCha(ich)
           end do !ich
c END      line
           write(kanalPDBres,7071)"END   " 
c
           if(OPT_SolvEFS) call writewBrgWmoLXYZ 
c
           iStatus = iStatus + 1
           write(kanalPStat,*)mStatus,iStatus,messageStatus
     &     ,' mdSnap : ', nRecPdb , ' is wrote ...'
c
           close(kanalPDBres)
c
c write arrayTra
cx          nStepTra = nStepTra + 1
cx        if(nStepTra .gt. nstepTraMAX) then
cx          write(kanalp,*)
cx     &    'ERROR: mdRun: nstepTraMAX is SMALL (engXYZtra.h)!'
cx          stop
cx          end if !nStepTra .gt. nstepTraMAX
c
          invSt = 1.0/nStepTra
          an1 = 1.0 - invSt
          an2 = (nStepTra - 1.0)/(1.0+nStepTra)**2
c
          eVbondDefTra(nStepTra) = eVbondDef
          eVangDefTra(nStepTra) = eVangDef
          eImpDefTra(nStepTra) = eImpDef
          eTorsDefTra(nStepTra) = eTorsDef
          engVDWR1Tra(nStepTra) = engVDWR1
          hBHxYeng128Tra(nStepTra) = hBHxYeng128 
cxPW:          solvMolEngF3 = molSolEn+hpSoLEnergy+watBrgEnergySoL
cxPW:     &                  + bornPolzEng
          eMolSolvTra(nStepTra) = solvMolEngF3
          
c
          engVDWR1TraAv = engVDWR1TraAv 
     &                    + (engVDWR1-engVDWR1TraAv)*invSt
          engVDWR1TraSd = an1*engVDWR1TraSd 
     &                    + an2*(engVDWR1-engVDWR1TraAv)**2
c
          engCOULR1Tra(nStepTra) = engCOULR1
          engCOULTra(nStepTra) = engCOUL
          engCOULTraAv = engCOULTraAv 
     &                   +  (engCOUL-engCOULTraAv)*invSt
          engCOULTraSd = an1*engCOULTraSd  
     &                   +  an2*(engCOUL - engCOULTraAv)**2
          engPOTENTTra(nStepTra) =  engPOTENT
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
cx PWex:
          kPW=1
          kinEngTra(nStepTra) = kinEng(kPW)
          tempT0tra(nStepTra) = tempT0(kPW)
          tempT0traAv = tempT0traAv + (tempT0(kPW) - tempT0traAv)*invSt
          tempT0traSd = tempT0traSd*an1 
     &                  + an2*(tempT0(kPW) - tempT0traAv)**2
c
          restr1EngTra(nStepTra) = restr1Eng
          restr1MHWEngTra(nStepTra) = restr1MHWEng
          eGeoDefTra(nStepTra) = eGeoDef
c
c write atom coords atomXYZtra()  for nStepTra
           sH = 3*natom*(nStepTra-1)
           do k=1,3*natom
           atomXYZtra(k+sH) = atomXYZ0(k)
	   atomXYZtmp(k) = atomXYZ0(k)
           end do !k
c
          engVDWR1TraSd = sqrt(engVDWR1TraSd)
          engCOULTraSd = sqrt(engCOULTraSd)
          engPOTENTTraSd = sqrt(engPOTENTTraSd)
          tempT0traSd = sqrt(tempT0traSd)
          hBHxYeng128TraSd = sqrt(hBHxYeng128TraSd)  
          engMolSolSd = sqrt(engMolSolSd)
c
          end if !nwtra !writeTraj
          end if !wtra
          end do ! kf3 largest timestep
          ntime0 = ntimew
c
c write Md Final PDB
        if(OPT_wrtieMdFinPDB)then
c open FinPDB file
        fileFinPdbN = pdbFinFile(1:nLfileFinPdb)
cx        kanalPDBfin = 34
        open(unit=kanalPDBfin, file=fileFinPdbN,form='formatted',
     &       status='unknown')
c
        rewind kanalPDBfin
c
           write(kanalPDBfin,'(a35,i5,a12,f10.4 )')
     &     "REMARK: Md result MdynSB03:ntimeMX:",ntimeMX,
     &     " MdTime(ps):", ntime*deltat
c
          write(kanalPDBfin,'(a9,i7)')'$nstep:  ', ntimew
          write(kanalPDBfin,'(a9,i7)')'$nRecPDB:', nRecPdb
	  write(kanalPdbfin,'(a18,f8.2)')
     &            'REMARK: RMSD(x0):  ',rmsdTraj(nStepTra)
c
          write(kanalPDBfin,'(a8,2x,a6)')'$ENERGY: ',' :Kcal'
          write(kanalPDBfin,'(a10,f14.5)')'eVbondDef:',eVbondDef
          write(kanalPDBfin,'(a10,f14.5)')'eVangDef :',eVangDef
          write(kanalPDBfin,'(a10,f14.5)')'eImpDef  :',eImpDef
          write(kanalPDBfin,'(a10,f14.5)')'eTorsDef :',eTorsDef
          write(kanalPDBfin,'(a10,f14.5)')'engVDWR1 :',engVDWR1
          write(kanalPDBfin,'(a10,f14.5)')'engHbHxY :',hBHxYeng128
          write(kanalPDBfin,'(a10,f14.5)')'eCOULR1  :',engCOULR1
          write(kanalPDBfin,'(a10,f14.5)')'eCOULR2  :',engCOULR2
          write(kanalPDBfin,'(a10,f14.5)')'eRstHA1  :',restr1Eng
          write(kanalPDBfin,'(a10,f14.5)')'eRstHW1RS:',restr1RHWEng
          write(kanalPDBfin,'(a10,f14.5)')'eRstHW1ML:',restr1MHWEng
          write(kanalPDBfin,'(a10,f14.5)')'eRstHA2  :',restrDistA2Eng
          write(kanalPDBfin,'(a10,f14.5)')'eCompcEng:',compactGeoEn
          write(kanalPDBfin,'(a10,f14.5)')'eGeoDef  :',eGeoDef
          write(kanalPDBfin,'(a10,f14.5)')'eCOUL    :',engCOUL
          write(kanalPDBfin,'(a10,f14.5)')'eSolvGS  :',molSolEn
          write(kanalPDBfin,'(a10,f14.5)')'eSolvCav :',hpSoLEnergy
          write(kanalPDBfin,'(a10,f14.5)')'eGBnPolz :',bornPolzEng
          write(kanalPDBfin,'(a10,f14.5)')'eWatBridg:',watBrgEnergySoL
          write(kanalPDBfin,'(a10,f14.5)')'eSolvTot :',solvMolEngF3
          write(kanalPDBfin,'(a10,f14.5)')'engPOTENT:',engPOTENT
          write(kanalPDBfin,'(a10,3f12.3)')'kinEng   :',kinEng
          write(kanalPDBfin,'(a10,3f12.5)')'Temp(K)  :',tempT0
          write(kanalPDBfin,'(a10,f14.5)')'engTotal :',engTotal
c
          write(kanalPDBfin,'(a30)')'- - - - - - - - - - - - - - - -'    
          write(kanalPDBfin,'(a19,i8)')'AVERAGE over nSnap:',nStepTra
          write(kanalPDBfin,'(a14,f14.5)')
     &                              'eVDWR1TraAv    :',engVDWR1TraAv
          write(kanalPDBfin,'(a14,f14.5)')
     &                              'eVDWR1TraSd    :',engVDWR1TraSd
          write(kanalPDBfin,'(a14,f14.5)')
     &                              'eCOULTraAv     :',engCOULTraAv
          write(kanalPDBfin,'(a14,f14.5)')
     &                              'eCOULTraSd     :',engCOULTraSd
          write(kanalPDBfin,'(a17,f14.5)')
     &                             'eHBHxY128TraAv :',hBHxYeng128TraAv
          write(kanalPDBfin,'(a17,f14.5)')
     &                             'eHBHxY128TraSd :',hBHxYeng128TraSd
          write(kanalPDBfin,'(a14,f14.5)')
     &                              'eSolvTotAv     :',engMolSolAv   
          write(kanalPDBfin,'(a14,f14.5)')
     &                              'eSolvTotSd     :',engMolSolSd
          write(kanalPDBfin,'(a14,f14.5)')
     &                              'engPOTTraAv    :',engPOTENTTraAv
          write(kanalPDBfin,'(a14,f14.5)')
     &                              'engPOTTraSd    :',engPOTENTTraSd
          write(kanalPDBfin,'(a14,f14.5)')
     &                              'tempT0traAv    :',tempT0traAv
          write(kanalPDBfin,'(a14,f14.5)')'tempT0traSd   :',tempT0traSd 
c
c control hb128
          if(CONTROLFHb128)then
          call printHB128EngForces
          end if ! CHb128
c
          if(vLigFlag .eq. 1) then
          write(kanalPDBfin,'(a8)')'########'
          write(kanalPDBfin,'(a8)')'$ENELIG: '
          write(kanalPDBfin,'(a12,f14.5)')'eVbondDefLG:',eVbondDefLg
          write(kanalPDBfin,'(a12,f14.5)')'eVangDefLG :',eVangDefLg
          write(kanalPDBfin,'(a12,f14.5)')'eImpDefLG  :',eImpDefLg
          write(kanalPDBfin,'(a12,f14.5)')'eTorsDefLG :',eTorsDefLg
          write(kanalPDBfin,'(a12,f14.5)')'engVDWR1LG :',engVDWR1Lg
          write(kanalPDBfin,'(a14,f14.5)')'hBHxYeng128Lg:',hBHxYeng128Lg 
          write(kanalPDBfin,'(a12,f14.5)')'engCOULR1LG:',engCOULR1Lg
          write(kanalPDBfin,'(a12,f14.5)')'engCOULR2LG:',engCOULR2Lg
          write(kanalPDBfin,'(a12,f14.5)')'restr1EngLG:',restr1EngLg
          write(kanalPDBfin,'(a12,f14.5)')'eRstHW1MLLG:',restr1MHWEngLg
          write(kanalPDBfin,'(a12,f14.5)')'eGeoDefLG  :',eGeoDefLg
          write(kanalPDBfin,'(a12,f14.5)')'engCOULLG  :',engCOULLg
          write(kanalPDBfin,'(a12,f14.5)')'engSolvLG  :',molSolEnLg
          write(kanalPDBfin,'(a12,f14.5)')'engPOTENTLG:',engPOTENTLg
          write(kanalPDBfin,'(a7)')'$ENDLIG '
          end if !ELig
c
          write(kanalPDBfin,'(a12)')'REMARK: PDB: '
c
           do ich=1,nChain
           do k = startAtInCha(ich),stopAtInCha(ich)
           k3 = 3*k-3
           doWRpdb=.true.
c
           if(OPT_doLigDock .eq. 1)then
           doWRpdb=.false.
           if(resName(k) .eq. resName(atomInLigList(1))) doWRpdb=.true.
           end if !OPT_doLigDock .eq. 1
c
           if(doWRpdb)then
           if(OPT_writeXYZForce)then
           write(kanalPDBfin,7072)"ATOM  ",
     &     atomNumb(k),atomName(k),resName(k),chName(k),resNumb(k),
     &     (atomXYZ0(k3+j),j=1,3),
     &     (atomVel0(k3+j),j=1,3)

           else
           write(kanalPDBfin,7071)"ATOM  ",
     &     atomNumb(k),atomName(k),resName(k),chName(k),resNumb(k),
     &     (atomXYZ0(k3+j),j=1,3), atomQ0(k)
           end if !OPT_writeXYZForce 
           end if !doWRpdb
c
           end do !k
c TERminal chainLine
           if( (resName(stopAtInCha(ich))(1:3) .ne. 'HOH' .and.
     &        resName(stopAtInCha(ich))(1:3) .ne. 'ION') .or.
     &        ich .eq. nChain)
     &     write(kanalPDBfin,7071)"TER   ", stopAtInCha(ich)           
           end do !ich
c END line
           write(kanalPDBfin,7071)"END   " 
c
           close(kanalPDBfin)
c
           end if !OPT_writeMdFinPDB
c
        if(cltra .eq. 1)then
        close(kanalXYZtra)
        close(kanalEngTra)
        end if
c
	if(OPT_SolvEFS) call writewBrgWmoLXYZ
c
        if(CONTROL)then
        write(kanalp,*)'mdRun Finish: StepDone=ntime:',ntime
        end if
c
         return
7071    format(a6,1x,i4,2x,a4,a4,a1,i4,4x,3f8.3,7x,f8.5) ! PDB rq
7072    format(a6,1x,i4,2x,a4,a4,a1,i4,4x,3f8.3,1x,3f8.3) ! PDB vel
7073    format(a6,1x,i4,2x,a4,a4,a1,i4,4x,3f8.3,1x,3f8.2) ! PDB vel
	 end
