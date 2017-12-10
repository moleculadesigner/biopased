c data arays for restrain1RHW 
c harmonic well restraint for REsidue CMass
c
        character*(charLenMAX) restr1RHWFile
	integer nRestr1RHWAtomMAX
        parameter (nRestr1RHWAtomMAX = natomMAX)
        integer nRestr1RHWAtom
        integer restr1RHWatomList(nRestr1RHWAtomMAX)
        integer nRestr1RHWSegMAX
        integer nRestr1RHWSegMX
        integer nRestr1RHWSeg    ! number of restrained segments
        parameter (nRestr1RHWSegMAX = nresMAX)
        integer resEndRestr1RHW(2*nRestr1RHWSegMAX)      ! start/end Res numbers
c restrains of type1HW = oneAtom Harmonic Well restrains
        integer nRestr1ResidHWMAX
        parameter (nRestr1ResidHWMAX=nresMAX)
        integer nRestr1ResidHW
        integer restr1ResidHWList(nRestr1ResidHWMAX)
        real restr1ResidMass(nRestr1ResidHWMAX)
        real restr1RefResidPosHW(3*nRestr1ResidHWMAX)
        real restr1RHWConst
        real sizeRestr1RHW
c
        common/restr1RHW00/restr1RHWFile
        common/restr1RHW01/resEndRestr1RHW
        common/restr1RHW02/nRestr1RHWSeg,nRestr1RHWSegMX         
	common/restr1RHW03/nRestr1ResidHW,restr1RHWConst,
     &        sizeRestr1RHW,restr1ResidHWList,
     &        restr1RefResidPosHW,restr1RHWatomList,
     &        restr1ResidMass
c
cend
