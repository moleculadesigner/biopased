c global variables
c PWatEx:2009
c
c movingAtom List
c
        integer nmoveatom
cx        integer nAtFreeDg        ! number of freedom degrees
        integer nAtFreeDg(3)       !1=whole system, 2=Psolute, 3=watSolvent
	integer nMovingAtom(3)     !
cx	real kinetEnergy(3)
        integer moveAtomList(natomMAX)
        integer moveFlag(natomMAX)
c
	common/moveAtom/nmoveatom,moveAtomList,moveFlag,
     &	nAtFreeDg,nMovingAtom
c
        real moveAtomXYZ(3*natomMAX)
        common/moveAtom2/moveAtomXYZ
cend
