c  initRigBodyMolDynamics                     
	subroutine initRigBodyMdStart
c
        include 'xyzPDBsize.h'
        include 'xyzPDBinfo.h'
        include 'xyzPDBcrd.h'
        include 'pair1234array.h'
        include 'nbondPairVCS.h'
        include 'movingAtom.h'
        include 'enForce.h'
        include 'restrainInfo.h'
        include 'mdAtomXYZvel.h'
        include 'mdRunPar.h'
        include 'shake.h'
        include 'optionPar.h'
        include 'solvGSarray.h'
        include 'rigBody01.h'
clocal       
        integer rgbFflag
        real dv,dt12,scf,dt05
        integer i,ia,ia3,ir
        integer k,k3
        integer doOmega
        integer kanalp
        logical doFastF1,doMediF2,doSlowF3
        logical CONTROL
c
        kanalp = 6
        CONTROL = .true.
c
        doFastF1 = .true.
        doMediF2 = .true.
        doSlowF3 = .true.
c
	if(CONTROL)then
        write(kanalp,*)'* * * * * * * * * * *  * * * * * * * *'
        write(kanalp,*)'*    initRigBodyMDStart: start       *'
        write(kanalp,*)'* * * * * * * * * * * * * *  * * * * *'
        end if
c
c positions at time = 0 
cx       atomXYZ0(i) = atomXYZ(i)
cx       atomVel0(i) 
c generate
cx       atomVelp(i) at t=t0+dt/2
cx       atomVelm(i) at t=t0-dt/2      

        if(nRigBody .eq. 0) return 
c
        doOmega = 1  
        call  getRigBodyParam01(doOmega,atomXYZ0,atomVel0)
c
        call getRigBRotTraForce(atForceTot)
c
        dt05 = deltat*0.5
        call initRigBodyVel(atomXYZ0,dt05)                
c
c correct number of freedom degres for RigibBody 
        nAtFreeDg = 3*nmoveatom - nbondShake - 3*nAtRigBody + 6*nRigBody
c
	if(CONTROL)then
        write(kanalp,*)'initRigBodyMDStart: nAtFreeDg:',nAtFreeDg
        write(kanalp,*)'initRigBodyMDStart: DONE :'
        end if
c
7071    format(a6,1x,i4,2x,a4,a4,a1,i4,4x,3f8.3,7x,f8.5) ! PDB rq 
c
        return
        end
