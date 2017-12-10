c  Yuri Vorobjev 2002
c               -2004
c  initMolDynamics                     
c  init:
c  1) start XYZ = atomXYZ() in xyzPDBcrd.h
c              md works with atomXYZ0()   = XYZ at current time ti
c  2) atomVel0(),atomVelm,atomVelp
c
	subroutine initMdStart(temp)
c
        include 'xyzPDBsize.h'
        include 'xyzPDBinfo.h'
        include 'xyzPDBcrd.h'
c
        include "output.h"
        include "kanalUse.h"
        include "statusMessg_mDyn.h"
c
        include 'pair1234array.h'
        include 'nbondPairVCS.h'
        include 'movingAtom.h'
        include 'enForce.h'
        include 'restrainInfo.h'
        include 'mdAtomXYZvel.h'
        include 'mdRunPar.h'
        include 'optionPar.h'
        include 'solvGSarray.h'
        include 'filedat.h' 
        include 'hbond128.h'
	include "solvate01.h"
	include "exWatShRestrSASPot.h"
c
        real temp
clocal       
        real dv,dt12,scf
        integer i,i3,k,ia,iam,ia3
        integer kanalp
        logical doFastF1,doMediF2,doSlowF3
        logical CONTROL
        kanalp = kanalRunOut 
        CONTROL = .true.
c
        doFastF1 = .true.
        doMediF2 = .true.
        doSlowF3 = .true.
c
        ncall_exWSASPot=0   !cunter exWatShRestrSASPot.h
c
	if(CONTROL)then
        write(kanalp,*)'* * * * * * * * * * * * * * * *'
        write(kanalp,*)'*    initMDStart: start       *'
        write(kanalp,*)'* * * * * * * * * * * * * * * *'
        end if
c
c positions at time = 0 
c
       if(OPT_mdRestart)then
       call readPDBxyzV01(mdRestXYZVfile,natom,atomXYZ,atomVel0,
     &       nRecPdb,ntime0) 
       end if!OPT_mdRestart
c
       do i = 1,3*natom
       atomXYZ0(i) = atomXYZ(i)
       end do !i
c
c generate all pair Lists and forces
        makeVdW = 1
        makeCL = 1
        makeSL = 1
c
        ncallNBPL = 0
        call initNonBondList(atomXYZ,makeVdW,makeCL,makeSL) 
c
c update hB128 tripletList
        if(OPT_HBond128)then
        call getHbTripletAll(atomXYZ,
     &           natom,nmoveatom,moveAtomList,
     &           nbpairListV,startnbPairLV,nnbPairLV,
     &           rcutHb128,nHBbondHxYmx,
     &           nHBbondHxY,hB128List,hB128TypeList,hB128PotParList) 
c
        do i = 1,nHBbondHxY
        if(hB128TypeList(i) .eq. 1)hB128WfList(i) = fEngWF(10)
        if(hB128TypeList(i) .eq. 2)hB128WfList(i) = fEngWF(11)
        end do !i
        end if !hb128
c
cx        call initAtomDconst(natom,atomXYZ)
c
        call initSASdielModLaZ(atomXYZ)
c
        fcall = 0
	call initAllForce00(fcall)
c
        fcall = 0
        if(doFastF1)then
        call initAllForce01(fcall,fEngWF,atomXYZ,
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
        end if ! doFastF1
c
        if(doMediF2) then
c medium Range COULR2 
        call allAtVDWEForceR2(atomXYZ,atomQ,
     &           natom,nmoveatom,moveAtomList,
     &           nbpairListC,startnbPairLC,nnbPairLC,
     &           rcutV,rcutC,engCOULR2,coulForceR2)
 
         end if ! coulomb R2
c
c solvent Forces
         if (OPT_SolvGS .and. doSlowF3 ) then
         call SolventEnForces(natom, atomXYZ,
     &         atomName,startPairL12,nPairL12,pair12List,
     &         nbpairListS,startnbPairLS,nnbPairLS,
     &         atomSolPar, molSolEn, atomSolEn, atomSolFr)
 
         end if !doSlowF3
c
c generate Random initial Velocity at time t=0
        if(.not. OPT_mdRestart)then 
	call initVelocity(temp,natom,
     &       nmoveatom,moveAtomList,atomMass,atomVel0)
        end if ! init V
c 
c estimate Velocity at time (t-dt/2) = atomVelm()
c                           (t+dt/2) = atomVelp()
c      scf = 418.40              ! scaling coef if force[Kcal/mol*A]
c                                  to get deltaV in [A/ps]
       scf = 418.40*deltat*0.5   ! 
       do iam = 1,nmoveatom
       ia = moveAtomList(iam)
       dt12 = scf/atomMass(ia)
       i3=3*ia-2
       do i=i3,i3+2
c
       atomForceF2(i)=coulForceR2(i)*fEngWF(7)
       atomForceF3(i)=atomSolFr(i)*fEngWF(9)
       atForceTot(i)=atomForceF1(i)+atomForceF2(i)+atomForceF3(i)
c
        dv = dt12*atForceTot(i) 
        atomVelm(i) =  atomVel0(i) - dv
        atomVelp(i) =  atomVel0(i) + dv
       end do !i
       end do !iam
c
	if(CONTROL)then
        write(kanalp,*)'initMDStart: DONE initMDStart:'
        end if
c
         return
c
7071   format(a6,1x,i4,2x,a4,a4,a1,i4,4x,3f8.3,2x,f6.3,f8.3) ! PDB
	 end
