c global variables
c energy and XYZ short traj
        integer nStepTraMAX
cX        parameter (nStepTraMAX=1000)  ! max 1 ns in one mdRun call
cX        parameter (nStepTraMAX=2500)  ! max 2.5ms in one mdRun call : 2015
        parameter (nStepTraMAX=4500)  ! max 4.5ms
c
        integer nStepTra
        real eVbondDefTra(nStepTraMAX)
        real eVangDefTra(nStepTraMAX)
        real eImpDefTra(nStepTraMAX)
        real eTorsDefTra(nStepTraMAX)
        real engVDWR1Tra(nStepTraMAX)
        real engCOULR1Tra(nStepTraMAX)
        real restr1EngTra(nStepTraMAX)
        real restr1MHWEngTra(nStepTraMAX)
        real eGeoDefTra(nStepTraMAX)
        real engCOULTra(nStepTraMAX)
        real engPOTENTTra(nStepTraMAX)
        real kinEngTra(nStepTraMAX)
        real engTotalTra(nStepTraMAX)
        real tempT0tra(nStepTraMAX)
        real hBHxYeng128Tra(nStepTraMAX)
        real eMolSolvTra(nStepTraMAX)
c
        real eVbondDefTraAv
        real eVangDefTraAv
        real eImpDefTraAv             
        real eTorsDefTraAv           
        real engVDWR1TraAv  
        real engVDWR1TraSd
        real engCOULR1TraAv
        real restr1EngTraAv
        real restr1MHWEngTraAv
        real eGeoDefTraAv
        real engCOULTraAv
        real engCOULTraSd
        real engPOTENTTraAv
        real engPOTENTTraSd
        real engMolSolAv,engMolSolSd
        real kinEngTraAv
        real engTotalTraAv
        real tempT0traAv
        real tempT0traSd
        real hBHxYeng128TraAv
        real hBHxYeng128TraSd
        real eGeoDefTraSd
        real eRestrALL
        real eRestrALLAv
        real eRestrALLSd
c
        common/engTra01/ nStepTra
        common/engTra02/ eVbondDefTra,eVangDefTra,
     &                   eImpDefTra,eTorsDefTra,
     &                   engVDWR1Tra,engCOULR1Tra,
     &                   restr1EngTra,restr1MHWEngTra,
     &                   eGeoDefTra,engCOULTra,
     &                   engPOTENTTra,kinEngTra,
     &                   engTotalTra,tempT0tra,
     &                   eMolSolvTra,hBHxYeng128Tra
c
        real essModeEntropyMdTra
        real freeEngMdTra
       common/engTra03/  eVbondDefTraAv,eVangDefTraAv,eImpDefTraAv,
     &  eTorsDefTraAv,engVDWR1TraAv,engVDWR1TraSd,engCOULR1TraAv,
     &  restr1EngTraAv,restr1MHWEngTraAv,eGeoDefTraAv,
     & engCOULTraAv,engCOULTraSd,engPOTENTTraAv,engPOTENTTraSd,
     & engMolSolAv,engMolSolSd,kinEngTraAv,engTotalTraAv,
     & tempT0traAv,tempT0traSd,hBHxYeng128TraAv,hBHxYeng128TraSd,
     & eGeoDefTraSd,eRestrALL,eRestrALLAv,eRestrALLSd,
     & essModeEntropyMdTra,freeEngMdTra
c
        real atomXYZtra(3*natomMAX*nStepTraMAX) 
        common/xyzTra04/atomXYZtra
c
        real rmsdTraj(nstepTraMAX)
        common/engTra05/rmsdTraj
c
c Stucture of Minimal engPOT along mdTraj: is done for mdRunRepL
        real engPOTtraMin,engPOT1traMin
        real atomXYZVtraMin(6*natomMAX)
        common/engTra06/ engPOTtraMin,engPOT1traMin,
     &                   atomXYZVtraMin
c
