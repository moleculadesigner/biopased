c ouputFIle
c
cx        include "charStringSiz.h"
c
	character*(charLenMAX) outputFile
         
        common/outFile01/outputFile
c
        character*(charLenMAX) outErrFile
        common/outFile02/outErrFile
c
c runControl
        integer kanalSTD6
        integer kanalPStat     ! PrintStatus
        integer kanalRunOut
        integer kanalRunErr
        common/kanal00/kanalPStat,kanalRunOut,kanalRunErr,
     &                 kanalSTD6
c
        logical OPT_OUTshort,OPT_OUTfull
        logical OPT_WrPDBq,OPT_WrPDB,OPT_WrBabBondAng
        common/outFile03/OPT_OUTshort,OPT_OUTfull,
     &          OPT_WrPDBq,OPT_WrPDB,OPT_WrBabBondAng
c
