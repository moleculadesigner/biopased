c arays to keep restrainAt1 Info
c modified: 2009 - each res-res has own Kharm constant
c
        character*(charLenMAX) restr1File
	integer nrestr1AtomMAX
        parameter (nrestr1AtomMAX = natomMAX)
        integer nRestr1SegMAX
        integer nRestr1SegMX
        integer nRestr1Seg    ! number of restrained segments
        parameter (nRestr1SegMAX = nresMAX)
        integer resEndRestr1(2*nRestr1SegMAX)      ! start/end Res numbers
c
        real resResSegKhConst(nRestr1SegMAX)  !harmonicConst for res-res segment
c
c restrains of type1 = oneAtom restarins
        integer nRestr1atom    
        integer restr1atomList(nrestr1AtomMAX)
        real refAtomPos(3*nrestr1AtomMAX)
        real restr1AtConstList(nrestr1AtomMAX)
        real restr1AtConst
c
        common/restr1Inf00/restr1File,resEndRestr1,resResSegKhConst
c
	common/restr1Inf01/nRestr1atom,restr1atomList,
     &          refAtomPos,restr1AtConst,restr1AtConstList
c
        character*40 resRest1AtLine(nRestr1SegMAX)
        common/restr1Inf02/resRest1AtLine
c
c targAtPosition data:
        character*(charLenMAX) targAtPosFile
        real targAtPosHConst         ! harmonic constant
        common/targAtPos01/targAtPosFile,targAtPosHConst  
c
        integer natomTargAtPos
        character*4 atNameTargAtPos(natomMAX)
        character*4 resNameTargAtPos(natomMAX)
        integer resNumbTargAtPos(natomMAX)
        integer atNumbTargAtPos(natomMAX)
        real atXYZTargAtPos(3*natomMAX)
c
        common/targAtPos02/natomTargAtPos,atNameTargAtPos,
     &  resNameTargAtPos, resNumbTargAtPos,atNumbTargAtPos,
     &  atXYZTargAtPos 
c
cend
