c global variables
c        logical  OPT_LoopMD    !optionPar.h
c loopModeller Info: atoms in loops
c file name loopMD.inp
        integer nLoopMAX
        parameter (nLoopMAX = 20)
        integer resEndLoop(2*nLoopMAX)
        integer nAtomInLoopMAX
        parameter (nAtomInLoopMAX=natomMAX)  !5000 !
        integer atomInLoopList(nAtomInLoopMAX)  ! atomInLoopList(iaL) = global iaNumb=iaGL
c                                               ! iaL = 1,..,nAtomInLoop 
c                                               ! local numeration for LoopAtoms 
        integer nLoop,nAtomInLoop
        integer nLoopMX
        integer nAtomInLoopMX,nAnchorAtomMX
        integer startAtInLoop(nLoopMAX+1)       ! startAtInLoop(iLoop)=iL in local numeration
        integer stopAtInLoop(nLoopMAX)
        integer nMoveLoop                   ! number of moving Loops in movingLoops subSet
        integer moveLoopList(nLoopMAX)      ! list of movingLoops subSet
        integer loopType(nLoopMAX)          !1/2 : 1 - N/C end, 2-internalLoop
c
c anchor's for loop ends
        integer addAnchorAtF                ! 0/1 addAnchorAtFlag
        integer nAnchorAtomMAX
        parameter (nAnchorAtomMAX = natomMAX)
        integer nAnchorAtom 
        integer anchorAtomList(nAnchorAtomMAX)
        real refPosAnchorAtom(3*nAnchorAtomMAX) 
c
	common/loopInfo/ nLoop,nAtomInLoop,resEndLoop,
     &                   atomInLoopList,nAnchorAtom,
     &                   anchorAtomList,refPosAnchorAtom,
     &                   startAtInLoop,stopAtInLoop,
     &                   nMoveLoop,moveLoopList,loopType
c end  
