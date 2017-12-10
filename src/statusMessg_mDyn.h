c RUN-status variables
	integer iStatus
	character*7 mStatus
        character*25 messageStatus
c approximate CPU% for modules
        integer iMolTopoMX
        integer iFFbuildMX
        integer iSolvModMX
        integer iengOptimizMX
        integer imdStartRunMX
        integer imdSArunMX
        integer iCPUMX
        integer iStatRepFreq
c 
        character*6 mError
        character*(charLenMAX) messageError
c
        common/status00/mStatus,messageStatus,iStatus
c
        common/statu01/mError, messageError
c
        common/status02/iMolTopoMX,iFFbuildMX,iengOptimizMX,
     &         imdStartRunMX,imdSArunMX,iCPUMX,iStatRepFreq,
     &         iSolvModMX 
