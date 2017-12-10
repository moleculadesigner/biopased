c writeAtXYZinLigLocation  
c 
c write on All Prot Atoms if dist At - LigAtoms is lower then RLigLoc
c
c Yuri Vorobjev  2004                           
c
	subroutine wrAtXYZinLigLoc
c
        include 'xyzPDBsize.h'
c
        include "output.h"
        include "kanalUse.h"
        include "statusMessg_mDyn.h"
c
        include 'xyzPDBinfo.h'
c        include 'mdAtomXYZvel.h'
        include 'xyzPDBcrd.h'
        include 'ligInfo.h'
c local
        character*(charLenMAX) fileAtinLigLoc
        real rLigLoc
c
        integer LigLocResList(nresMAX)
        integer LigLocResIn(nresMAX)
        integer nLigLocRes
c  local
        integer kanalLigloc
        integer kanalp
        integer ich,k,k3,j,i
        integer iLm,ir
        real rKLig
        logical doWRpdb
        logical CONTROL
c
        kanalp = kanalRunOut
        CONTROL = .true.     

        fileAtinLigLoc = molAtXYZIntLigFile ! ligInfo.h 
        rLigLoc = rAtInLig   ! ligInfo.h
c
        if(CONTROL)then
        write(kanalp,*)'writeAtXYZinLigLoc: rLigLoc: ',rLigLoc
        end if
c
        kanalLigloc = 32
c
        open(unit=kanalLigloc, file=fileAtinLigLoc, form='formatted',
     &       status='unknown')
c
        write(kanalLigloc,'(a32,f8.2)')
     &       'REMARK: AtXYZinLigLoc: rLigLoc: ',rLigLoc
c init
         nLigLocRes = 0
         do i=1,nresMAX
         LigLocResList(i) = 0
         LigLocResIn(i) = 0
         end do
c
           do ich=1,nChain
           do k = startAtInCha(ich),stopAtInCha(ich) 
           k3 = 3*k-3
c
           doWRpdb=.false.
c calculate nearest dist from atom k to any LigAtom
           call distAtKtoLigAt(k,rKLig,iLm) 
           if(rKLig .le. rLigLoc)doWRpdb=.true.
c
	  if(doWRpdb)then
          ir = resNumb(k)
          if(LigLocResIn(ir) .eq. 0)then
          nLigLocRes = nLigLocRes + 1
          LigLocResList(nLigLocRes) = ir   
          LigLocResIn(ir) = LigLocResIn(ir) + 1
          end if  
          end if
c
          end do!k
          end do!ich
c
c write file
          if(nLigLocRes .gt. 0)then
          do i = 1,nLigLocRes
          ir = LigLocResList(i)
          do k = startAtInRes(ir),stopAtInRes(ir) 
          k3=3*k-3 
           write(kanalLigloc,7071)"ATOM  ",
     &     atomNumb(k),atomName(k),resName(k),chName(k),resNumb(k),
     &     (atomXYZ(k3+j),j=1,3), atomQ0(k)
           end do !k
           end do !i
c END      line
           write(kanalLigloc,7071)"END   " 
c
           end if !nLigLocRes .gt. 0
c
           close(kanalLigloc)
c
         return
7071    format(a6,1x,i4,2x,a4,a4,a1,i4,4x,3f8.3,7x,f8.5) ! PDB rq
	 end
