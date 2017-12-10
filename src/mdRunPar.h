c parameters of mdRUN
	integer ntimeMX           ! number of time step
        integer ntime             ! current time step
        integer ntime0            ! total timeStep done in all previous calls of mdRun
        integer ntimeR1           ! update frequency for R1 pairList
        integer ntimeR2           ! update frequency for R2 pairList
        integer ntimeF1           ! update frequency for R1 en/force
        integer ntimeF2           ! update frequency for R2 en/force
        integer ntimeF3           ! update frequency for Solv en/force
        real deltat               ! deltat time step [ps]
        integer atype             ! ensamble type 0-NEV, 1-NTV 
        integer optra,wtra,cltra  ! flags 0(1) to open,write,close TrajFiles
        integer nwtra             ! write each nwtra snap
c MD temperatureRelaxationPar
        real tempTg               ! target T of NTV ansemble
        real tauTRF               ! T realaxation factor
c
        common /mdRunGV/tauTRF,tempTg,
     &                  ntime,ntime0,ntimeMX,ntimeR1,ntimeR2,
     &                  ntimeF1,ntimeF2,ntimeF3,
     &                  deltat,atype,optra,wtra,cltra,
     &                  nwtra
c
