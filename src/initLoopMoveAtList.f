c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c                                                                   *
c  Yury Vorobjev 2003                                               *
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c 1) define: moveAtomList(ia) for a given list of moveLoopList_s(*)
c    of moving Loops
c
	subroutine initLoopMoveAtList(nMoveLoop_s,moveLoopList_s,
     &                OPT_LoopMD_s)
c
        include 'xyzPDBsize.h'
        include 'xyzPDBinfo.h'
c
        include "output.h" 
        include 'xyzPDBcrd.h'
        include 'filedat.h'
        include 'loopInfo.h'
        include 'movingAtom.h'
c
        integer nMoveLoop_s,moveLoopList_s(*)
        logical OPT_LoopMD_s
clocal
        integer kanalp
        integer i,ia,n,ilp,ilpm 
        integer i3,ia3
        logical CONTROL
c        
c initialize:
        kanalp = kanalRunOut
        CONTROL = .true.
c
         if(.not. OPT_LoopMD_s)then
c
c all realAtom(ia) [protCORE + LOOP] are moving 
         nmoveatom = 0    
         do i=1,natom    
         if(realAtomFlag(i) .eq. 1)then
         nmoveatom = nmoveatom + 1
         moveAtomList(nmoveatom) = i 
         moveFlag(i) = 1
         end if
         end do
c
	else      ! OPT_LoopMD_s
c
c define Loop atoms from a given loopList as a moving atoms 
         do i =1,natom
         moveFlag(i) = 0
         end do !i
         nmoveatom = 0 
c         
c define moving atom list
c
         if(nMoveLoop_s .ge. 1)then
         do ilp = 1,nMoveLoop_s
         ilpm = moveLoopList_s(ilp)
         do ia = startAtInLoop(ilpm),stopAtInLoop(ilpm)
         nmoveatom = nmoveatom + 1 
         moveAtomList(nmoveatom)=atomInLoopList(ia)
         moveFlag(moveAtomList(nmoveatom)) = 1
         end do !ia
         end do !ilp
         end if !nMoveLoop_s .ge. 1
c
        end if !OPT_fullProt
c
c define moveAtomXYZ()
        do i = 1,nmoveatom
        i3 = i*3-3
        ia3 =  3*moveAtomList(i)-3
        moveAtomXYZ(i3+1) = atomXYZ(ia3+1)
        moveAtomXYZ(i3+2) = atomXYZ(ia3+2)
        moveAtomXYZ(i3+3) = atomXYZ(ia3+3)
        end do !i
c
        if(CONTROL)then
        write(kanalp,*)'initMoveAtList:OPT_fullProtMD:'
     &  , (.not. OPT_LoopMD_s)      
        write(kanalp,*)
     &  'initMoveAtList: moveAtomList:nmoveatom:',nmoveatom
        write(kanalp,'(a21)')'i  iatmove  mvFlag:'
c
        do i = 1,nmoveatom
        ia = moveAtomList(i)
        write(kanalp,'(3i6)')i,ia,moveFlag(ia)
        end do !
        end if !CONTROL
c
        write(kanalp,*)'initMoveAtListinit: Finish:'
c
	return
        end
