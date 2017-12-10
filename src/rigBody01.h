c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c YN Vorobjev 2004
c 
c rigidBodyMd data structures
c
        character*(charLenMAX) rigBodyFile  
c
        common /rigBodyD00/rigBodyFile
c
	integer nRigBodyMAX
        integer nRigBodyMX
        integer nRigBody
        parameter (nRigBodyMAX=nresMAX)
c
        integer rigBodyStEndRes(2*nRigBodyMAX)
        integer startAtInRGB(nRigBodyMAX)
	integer stopAtInRGB(nRigBodyMAX)
        integer nAtRigBody                !number of atoms in all rigBodySegm
        integer nAtRigBodySeg(nRigBodyMAX) !number of atoms in the RGBsegment
        integer atRigBodyFlag(natomMAX)
c
	common /rigBodyD010/nRigBody,natRigBody,nAtRigBodySeg,
     &         startAtInRGB,stopAtInRGB,atRigBodyFlag
c
        real cAcARigBodyHK 
        common/rigBody011/cAcARigBodyHK
c --------------------------------------------------------------------------------
c RigBodyMD
        real rgbCmXYZ(3*nRigBodyMAX)     !rigBodycenter of Mass XYZ
        real rgbMass(nRigBodyMAX)        !rigBodyMass
        real rgbTensorInrt(9*nRigBodyMAX)   !tensorInercii
        real rgbTensorInrt01(9*nRigBodyMAX) !invertedTensorInert
        real rgbCmVel(3*nRigBodyMAX)        ! rigBodyCmassVelocity 
        real rgbTransForce(3*nRigBodyMAX)
        real rgbMomForce(3*nRigBodyMAX)
c
        common/rigBodyD020/rgbCmXYZ,rgbMass,
     &         rgbTransForce,rgbMomForce,
     &         rgbTensorInrt,rgbTensorInrt01
c
        real atForceRigBrot(3*natomMAX)
        real atForceRigBtra(3*natomMAX)
c
        common/rigBodyD021/ atForceRigBrot,atForceRigBtra     
c
        real atomXYZrigBodyFrame(3*natomMAX) ! atXYZ(ia) in the respective 
c                                            ! rigBodyFrames= x(i) - xCMass
        real atVelRigBodyFrame(3*natomMAX)   ! = v(i) - vCMass
        common/rigBodyD03/atomXYZrigBodyFrame,
     &                    atVelRigBodyFrame
c
        real atVelRigBodyRot(3*natomMAX)     ! rotationalVelocity t=t
        real atVelRigBodyRotMod(natomMAX)    ! modul rotationalVelocity
        real atVelRigBodyRotM(3*natomMAX)    ! t=t-dt/2
        real atVelRigBodyRotP(3*natomMAX)    ! t=t+dt/2
c
        real atVelRigBodyTra(3*natomMAX)     ! t = t = rgbCMVel(ir)
        real atVelRigBodyTraM(3*natomMAX)    ! t-dt/2
        real atVelRigBodyTraP(3*natomMAX)    ! t+dt/2
c
        common/rigBodyD031/atVelRigBodyRot,atVelRigBodyRotM,
     &                     atVelRigBodyRotP,atVelRigBodyRotMod
        common/rigBodyD032/atVelRigBodyTra,atVelRigBodyTraM,
     &                     atVelRigBodyTraP
c
         real atomXYZrigBparOmg(3*natomMAX)  ! atXYZ(ia) parallel to vect Omega
         real atomXYZrigBppdOmg1(3*natomMAX)  ! atXYZ(ia) perpedicular to vect Omega,
c                                             ! normalized to 1.0
c
        real atomXYZrigBMod(natomMAX)       ! atomXYZrigBodyFrame(ia)=modul
        common/rigBody032/atomXYZrigBparOmg,
     &         atomXYZrigBppdOmg1,atomXYZrigBMod
c
        real rgbKmomForce(3*nRigBodyMAX)     ! momForce
        real rgbPmomVel(3*nRigBodyMAX)       ! momImpulse 
        real rgbOmegaRot(3*nRigBodyMAX)     ! rgbVectorOmega t=t
        real rgbOmegaRotP(3*nRigBodyMAX)    ! t = t+dt/2
        real rgbOmegaRotM(3*nRigBodyMAX)    ! t = t-dt/2
        real rgbOmegaRot1(3*nRigBodyMAX)    ! Unit rgbVectorOmega 
        real rgbOmegaRotDt(3*nRigBodyMAX)   ! dOmega/dt
        real rgbOmegaRotMod(nRigBodyMAX)    ! modul of rgbVectorOmega
c
        common/rigBodyD04/
     &       rgbOmegaRot,rgbOmegaRotDt,rgbOmegaRotMod,
     &       rgbOmegaRotP,rgbOmegaRotM
c
        real rgbAtForceRot0(3*natomMAX)    ! rgbAtomForceRot0
c                           ! due to rgbRotation=centrostremitel'naja sila
