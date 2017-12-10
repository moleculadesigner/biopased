c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c                                                                           *
c  initialize: all forces on atoms to ZERO                                  * 
c                                                                           *
c Yury Vorobjev    2002                                                     *
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

	subroutine initAllForce00(vcall)
c
        include 'xyzPDBsize.h'
c
        include "output.h"
c
        include 'enForce.h'
        include 'compGeoEF.h'
c
        integer vcall    
        integer i,kanalp
c
       kanalp = kanalRunOut
c
c       write(kanalp,*)'initAllForce00 start!: vcall=',vcall
c
       if(vcall .ne. 0 ) return
c init to zero All forces on Atoms
       if(vcall .eq. 0)then
       eVbondDef = 0.0
       eVangDef = 0.0      
       eImpDef = 0.0
       eTorsDef = 0.0
       restr1Eng = 0.0          
       restr1MHWEng = 0.0  
        engVDWR1 = 0.0
       engCOULR1 = 0.0
       engFastF1 = 0.0
       engCOULR2 = 0.0
       molSolEn = 0.0
       hBHxYeng128 = 0.0 
       compactGeoEn=0.0
       watBrgEnergySoL = 0.0
       hpSoLEnergy = 0.0
       bornPolzEng = 0.0 
       eng_exWatSASP=0.0
c
       do i = 1,3*natomMAX
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
       atForceTot(i) = 0.0
       hBHxY128force(i) = 0.0
       compactGeoForce(i) = 0.0
       watBrgForceAt(i) = 0.0
       hpSoLForceAt(i) = 0.0
       bornPolzForceAt(i) = 0.0
       ff_exWatSASP(i) = 0.0
       atomForceF1(i) = 0.0
       atomForceF2(i) = 0.0
       atomForceF3(i) = 0.0
       end do !i
       end if !vcall
c
       vcall = vcall + 1
c
         return
         end
