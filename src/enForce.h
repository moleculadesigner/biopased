c global variables
c energy and force arrays
	real vbdefForce(3*natomMAX)
        real eVbondDef
c
        real vAngdefForce(3*natomMAX)
        real eVangDef
c
	real impDefForce(3*natomMAX)
        real eImpDef
c
        real eTorsDef
        real torsAngForce(3*natomMAX)
c
        real engVDWR1
        real vdwForceR1(3*natomMAX)
c
        real engCOULR1
c
        real coulForceR1(3*natomMAX)
        real engFastF1
        real atomForceF1(3*natomMAX) 
        real restr1Eng
        real restr1AtForce(3*natomMAX)
        real restr1MHWEng
        real restr1MHWAtForce(3*natomMAX)
        real restr1RHWEng
        real restr1RHWAtForce(3*natomMAX)
        real restrDistA2Eng
        real restrDistA2Force(3*natomMAX) 
c hb128 eng/force
        real hBHxYeng128,hBHxY128force(3*natomMAX)
c 
c in range R1 - R2 mediumForce F2
        real engCOULR2
        real coulForceR2(3*natomMAX)
        real atomForceF2(3*natomMAX) 
c
        real atForceTot(3*natomMAX)
        real eGeoDef,eGeoRst
        real engCOUL,engPOTENT,engPOTENT1
cx      real tempT0,kinEng,engTotal
c PWatEx:
        real tempT0(3),kinEng(3),engTotal
c
	common/enforceF1/ vbdefForce,eVbondDef,
     &                  vAngdefForce,eVangDef,
     &                  impDefForce,eImpDef,
     &                  torsAngForce,eTorsDef,
     &                  vdwForceR1,engVDWR1,
     &                  coulForceR1,engCOULR1,
     &                  restr1Eng,restr1AtForce,
     &                  restr1MHWEng,restr1MHWAtForce,
     &                  restrDistA2Eng,restrDistA2Force,
     &                  hBHxYeng128,hBHxY128force 
c
        common/enforceF01/engFastF1,atomForceF1
c
        common/enforceF2/coulForceR2,engCOULR2
        common/enforceF02/atomForceF2
c
c solvationEnForce
        real solvMolEngF3
        real atomForceF3(3*natomMAX)
        common/enforceF03/solvMolEngF3,atomForceF3
c GSsolvation
        real    atomSolEn(natomMAX)
        real    atomSolFr(3*natomMAX)
        real    molSolEn
        common/enforceSolGS/molSolEn,atomSolEn,atomSolFr
c
        real watBrgEnergySoL
        real watBrgForceAt(3*natomMAX)
        common/enforceFwBrg01/watBrgEnergySoL,watBrgForceAt
c
        real hpSoLEnergy
        real hpSoLForceAt(3*natomMAX)
        common/enforceHPSol01/hpSoLEnergy,hpSoLForceAt
c
        common/enforceFt/eGeoDef,engCOUL,engPOTENT,
     &                  tempT0,kinEng,engTotal,
     &                  eGeoRst,atForceTot 
c
        real bornPolzEng,bornPolzEngAt(natomMAX)
        real bornPolzForceAt(3*natomMAX)
        common/enforceGB01/bornPolzEng,bornPolzEngAt,
     &                     bornPolzForceAt
c
	integer nforceMAX
        parameter (nforceMAX=11)
        real fEngWF(nforceMAX)        !weighting factor for force/eng types
c                                     ! fEngWF(10)=hb128BB, fEngWF(11)=hb128SB,SS
        integer fcall                 !counter for number of forceUpdate
        common/enforceWF/fEngWF,fcall
c
        real eng_exWatSASP,ff_exWatSASP(3*natomMAX)
	common/exWatSASresr01/eng_exWatSASP,ff_exWatSASP
c
c
