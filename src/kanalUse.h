c
c mdynSB kanal in use list
c
c runControl -> output.h
cx	integer kanalPStat    !PrintStatus
cx        integer kanalRunOut
cx        integer kanalRunErr
cx        common/kanal00/kanalPStat,kanalRunOut,kanalRunErr
c inData
        integer kanalInPdb
        integer kanalInMdynPar
        integer kanalInMoveRes
        integer kanalInSAprot
        integer kanalInRestrALLType
        common/kanal01/kanalInPdb,kanalInMdynPar,
     &         kanalInMoveRes,kanalInSAprot,kanalInRestrALLType
c Lib:
        integer kanal_htg
        integer kanal_AmbZmConv
        integer kanal_ResTop01
        integer kanal_ResTop02
        integer kanal_ResTop03
        integer kanal_ResTop04
        integer kanal_solvAtType
        integer kanal_solvGSdat
        integer kanal_watBox0
c
cx Lib:
        common/kanal03/kanal_htg 
        common/kanal04/kanal_ResTop01,kanal_ResTop02,
     &                 kanal_ResTop03,kanal_ResTop04,
     &                 kanal_AmbZmConv,kanal_solvAtType,kanal_solvGSdat,
     &                 kanal_watBox0
c
cx Results:
        integer kanalXYZtra,kanalEngTra,kanalPDBres,kanalPDBfin
        integer kanalpHB128,kanalwWatBrg,kanalModFreq
        integer kanalModXYZ,engTablkanal1,engTablkanal2
        common/kanal05/kanalXYZtra,kanalEngTra,kanalPDBres,
     &                 kanalPDBfin,kanalpHB128,kanalwWatBrg,
     &                 kanalModXYZ,kanalModFreq,
     &                 engTablkanal1,engTablkanal2
c
