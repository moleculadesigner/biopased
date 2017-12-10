c md time step propagator for MultiTimeStep method        
c dynamics of nmoveatom is propagated for one smallest time step
c Velocity Verlet method 
c 
c IN: atomVelm()[t-dt/2]
c     atomXYZ0()[t]      : mdAtomXYZvel.h
c OUT:
c     atomXYZp()[t+dt]
c atomXYZp()[t+dt] -->SHAKE-->atomXYZ0()[t]
c atomVelp()[t+dt/2], atomVel0[t]
c atomVelm()[t-dt/2] --> atomVelp()[t+dt/2]-> atomVel0[t]
c atomVelm()[t+dt/2]=atomVelp()[t+dt/2]
c
c RESULT: all arrays in mdAtomXYZvel.h are defined
c
        subroutine mdTimeStepProp01(nmoveatom,moveAtomList,
     &                                      moveFlag,deltat)
c
        include 'xyzPDBsize.h'
        include 'xyzPDBinfo.h'
c
        include "output.h"
        include "kanalUse.h"
        include "statusMessg_mDyn.h"
c
        include 'enForce.h'
        include 'restrainInfo.h'
        include 'mdAtomXYZvel.h'
        include 'shake.h'
        include 'pair1234array.h'
        include 'rigBody01.h'
        include 'optionPar.h'
c
        integer nmoveatom,moveAtomList(*)
        integer moveFlag(*)
        real deltat
c
        real*8 uRotMx(3,3),tRt(3),tFx(3)
        integer ia,iam,k,i3
        real  dtm
        real rmsd0p
        integer kanalp
        logical CONTROL
        logical CONTROLF
        real scf
c
        kanalp = kanalRunOut
c
        CONTROL = .false. 
c        scf = 418.40            ! scaling coef if force[Kcal/mol*A]
c                                 to get deltaV in [A/ps]
        scf = 418.40*deltat        
c
        CONTROLF = .false.  !control forces for NAN
c
        if(CONTROL)then
        write(kanalp,*)'mdTimeStep01: start:'
        end if
c
        if(CONTROLF)then
        call printForces               ! 1
        end if !CF
c IN:
c atomXYZ(), all forces(*) at [time=t]
c atomVel()[t-dt/2] = atomVelm(*)
c OUT:
c propagate global atomXYZ to [t+dt] and atomVel[t+dt/2]
c
	do iam = 1,nmoveatom
        ia = moveAtomList(iam)
c
        i3=3*ia-2
        dtm = scf/atomMass(ia)
c
        do k = i3,i3+2
c calculate Vel for t+dt/2 
        if(.not. OPT_RigidBodyMd .or. atRigBodyFlag(ia) .eq. 0)then
        atomVelp(k) = atomVelm(k) + atForceTot(k)*dtm 
c calculate XYZp(t+dt) (unconstrained move)
        atomXYZp(k) = atomXYZ0(k) + atomVelp(k)*deltat
c
        else
c do RigBody translational move
        atVelRigBodyTraP(k) = atVelRigBodyTraM(k) + 
     &                        atForceRigBtra(k)*dtm
        atomXYZp(k) = atomXYZ0(k) + atVelRigBodyTraP(k)*deltat
c total Vel for RigBody atoms
        atomVelp(k) = atVelRigBodyTraP(k) + atVelRigBodyRotP(k)
c
        end if !.not. OPT_RigidBodyMd 
c
        end do !k 
        end do !iam
c
c update XYZ0
c
c eliminate WholeMolecule Translation and rotation
c kabsh rotation ! all atoms must be moving!!
         if(OPT_zeroROT)then
c rotate atomXYZp(*) to atomXYZ0(*)
c
cx         call xyz_rmsd(atomXYZp,atomXYZ0,natom,rmsd0p)
cx         write(kanalp,*)'mdTimeStep:xyz_rmsd beforeRot:',rmsd0p
c
         call  kabsch8(atomXYZp,atomXYZ0,natom,atomMass,uRotMx,
     &                                                 tRt,tFx)
c
cx         call xyz_rmsd(atomXYZp,atomXYZ0,natom,rmsd0p)
cx         write(kanalp,*)'mdTimeStep:xyz_rmsd afterRot:',rmsd0p
c
c update XYZ0
        do ia = 1,natom
        i3=3*ia-2
c redefine velocity after kabsch Rotation
        do k=i3,i3+2
        atomVelp(k) = (atomXYZp(k) - atomXYZ0(k))/deltat
        atomXYZp05(k) = 0.5*(atomXYZ0(k)+atomXYZp(k))
        atomXYZ0(k) = atomXYZp(k)
c velocity at time = t
        atomVel0(k) = (atomVelm(k) + atomVelp(k))*0.5
c update velocity
        atomVelm(k) = atomVelp(k)
        end do!k
        end do !ia   
        end if !OPT_zeroROT
c
c noSHAKE
c
        if(iShake .eq. 0)then
        do iam = 1,nmoveatom
        i3=3*moveAtomList(iam)-2
        do k=i3,i3+2
        atomXYZp05(k) = 0.5*(atomXYZ0(k)+atomXYZp(k))
        atomXYZ0(k) = atomXYZp(k)
c velocity at time = t
        atomVel0(k) = (atomVelm(k) + atomVelp(k))*0.5       
c update velocity
        atomVelm(k) = atomVelp(k)   
        end do!k
c
        end do !iam
        end if !iShake=0
c
        if(iShake .gt. 0) then 
c !doSHAKE
        call shakeBonds(iShake,nbond12,nbond12noH,nbondShake,
     &                     bond12List,bond12ParL,shitMX,shtol,
     &              moveFlag,atomMass,atomXYZ0,atomXYZp,
     &              shExitFlag)       
c 
c update XYZ0
        do iam = 1,nmoveatom
        i3=3*moveAtomList(iam)-2
c redefine velocity after shake
        do k=i3,i3+2
        atomVelp(k) = (atomXYZp(k) - atomXYZ0(k))/deltat
        atomXYZp05(k) = 0.5*(atomXYZ0(k)+atomXYZp(k)) 
        atomXYZ0(k) = atomXYZp(k)
c velocity at time = t
        atomVel0(k) = (atomVelm(k) + atomVelp(k))*0.5
c update velocity
        atomVelm(k) = atomVelp(k)
        end do!k
c
        end do !iam
c
        end if !SHAKE
         
        if(CONTROL)then
        write(kanalp,*)'mdTimeStep01: finish:'
        write(kanalp,*)'ia  atomVelp  NewatomXYZ  atForceF1:'
         do iam = 1,nmoveatom
         ia = moveAtomList(iam)
         i3=3*ia-3
         write(kanalp,'(i5,3f8.3,2x,3f8.3,2x,3f8.3)')
     &   ia,(atomVelp(i3+k),k=1,3),(atomXYZ0(i3+k),k=1,3),
     &      (atomForceF1(i3+k),k=1,3)
         end do
        end if !C
c
         return
	 end
c
c md time step propagator for MultiTimeStep method
c Velocity of nmoveatom is propagated for 1/2 of time step
c
c atomVel()[t] <--atomVel()[t-dt/2] + dt/2 * atForce
c
c RESULT: all arrays in mdAtomXYZvel.h are defined
c
        subroutine mdTimeStepProp02(nmoveatom,moveAtomList,
     &                 moveFlag,atomMass,deltat,atForce,atomVel)
c
        integer nmoveatom,moveAtomList(*)
        integer moveFlag(*)
        real atomMass(*)
        real atForce(*)
        real atomVel(*)
        real deltat
c
        integer ia,iam,k,i3
        real  dtm
        integer kanalp
        logical CONTROL
        logical CONTROLF
        real scf
c
        kanalp = 6
        CONTROL = .false.
c        scf = 418.40            ! scaling coef if force[Kcal/mol*A]
c                                 to get deltaV in [A/ps]
        scf = 418.40*deltat      
c
        CONTROLF = .false.  !control forces for NAN
c
        if(CONTROL)then
        write(kanalp,*)'mdTimeStep02: start:'
        write(kanalp,*)'atomVelp() ia x y z:'
        do iam = 1,nmoveatom
        ia = moveAtomList(iam)
        i3=3*ia-3
        write(kanalp,*) ia,atomVel(i3+1),atomVel(i3+2),atomVel(i3+3)
        end do
        end if
c
        if(CONTROLF)then
        call printForces
        end if !CF
c
c propagate atomVel(t-dt/2) to  atomVel[t+dt/2]
c
        do iam = 1,nmoveatom
        ia = moveAtomList(iam)
        dtm = scf/atomMass(ia)
        i3=3*ia-2
c
        do k = i3,i3+2
c update vel for +dt/2 (time=t)
        atomVel(k)=atomVel(k) + atForce(k)*dtm
        end do !k
        end do !iam
c
        if(CONTROL)then
        write(kanalp,*)'mdTimeStep02: finish:'
        write(kanalp,*)'atomVelp() ia x y z:'
         do iam = 1,nmoveatom
         ia = moveAtomList(iam)
         i3=3*ia-3
         write(kanalp,*) ia,atomVel(i3+1),atomVel(i3+2),atomVel(i3+3)
         end do
        end if !C
c
         return
         end
c
