c data structure for mobile ions aroun molecule
c
        integer nIonTypeMAX
        integer nIonTypeMX
        integer nIonType
        parameter (nIonTypeMAX=1)
        real zarIon(nIonTypeMAX)
c
	integer nIonIASMAX
        parameter (nIonIASMAX = 100)
c
        real eIonBindMIN_OPT      ! 0.6
        real ionEpotMAX_OPT      ! -2.50 kcal
        common/ionSAS01/eIonBindMIN_OPT,ionEpotMAX_OPT
c
        integer nCounterIonMX
        integer nIonNow(nIonTypeMAX)
        character*4 ionResName(nIonTypeMAX)
        character*4 ionAtomName(nIonTypeMAX)
        character*2 ionffatomName(nIonTypeMAX)
        common/ionSAS02/ionResName,ionAtomName,ionffatomName
        real ePotdotMX(nIonIASMAX)
        integer iDotePotdotMX(nIonIASMAX)
        common/ionSAS03/nCounterIonMX,nIonNow,ePotdotMX,
     &  iDotePotdotMX
        real ionIASXYZ(3*nIonIASMAX)
        real ionAtomQ(nIonIASMAX)
        real ionIasEpot(nIonIASMAX)    !pot Eng of Ion in the ion build up
c                                   procedure
        common/ionSAS04/ionIASXYZ,ionIasEpot,ionAtomQ
c
        character*4 ionIAsAtomName(nIonIASMAX),ionIasResName(nIonIASMAX)
        common/ionSAS05/ionIAsResName,ionIasAtomName
c
