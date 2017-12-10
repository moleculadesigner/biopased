c global variables
c VDW parameters for atomPairs  type1,type2
	integer nVDWtypeMAX
        parameter (nVDWtypeMAX = 60)
        integer nVDWatomParMAX
        parameter (nVDWatomParMAX = 4)    
        real atomVDW12ab(nVDWatomParMAX*nVDWtypeMAX*(nVDWtypeMAX+1)/2)
        integer nVDWtypeMX
        integer nVDWatomParMX
        integer nVDWtype               !number of vdwTypes in ff
c  A12, B12, rm12, ev12
 	common/vdwPar12AB/nVDWtype,atomVDW12ab
c
        real aSoftCore,aSoftCoreMIN      ! rSoftCore12 = aSoftCore*rmvdw12
        common/vdwParSC/ aSoftCore,aSoftCoreMIN
c
c end
