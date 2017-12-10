c controlCall.h
c parameters for controlPrint of subroutine
c
        integer nCaseMAX 
        parameter (nCaseMAX = 1)
 	integer ncallSubTot(nCaseMAX)
        integer nPrintNcall(nCaseMAX) 
c                        ! print each nPrintNcall of tot calls ncallSubTot
        common/controlCall01/ncallSubTot,nPrintNcall
c
