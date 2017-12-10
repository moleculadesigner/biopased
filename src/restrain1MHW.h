c arays to keep restrain1HW Info
        character*(charLenMAX) restr1MHWFile
	integer nRestr1MHWAtomMAX
        parameter (nRestr1MHWAtomMAX = natomMAX)
        integer nRestr1MHWSegMAX
        integer nRestr1MHWSegMX
        integer nRestr1MHWSeg    ! number of restrained segments
        parameter (nRestr1MHWSegMAX = nresMAX)
        integer resEndRestr1MHW(2*nRestr1MHWSegMAX)      ! start/end Res numbers
c restrains of type1HW = oneAtom Harmonic Well restrains
        integer nRestr1MHWatom    
        integer restr1MHWatList(nrestr1MHWAtomMAX)
        real atMs1MHWatList(nrestr1MHWAtomMAX)
        real ref1MHWPos(3)
        real restr1MHWConst
        real sizeRestr1MHW
c
        common/restr1HW00/restr1MHWFile
	common/restr1HW01/nRestr1MHWatom,restr1MHWatList,
     &                   atMs1MHWatList,
     &                   ref1MHWPos,restr1MHWConst,sizeRestr1MHW
c
cend
