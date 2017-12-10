c LoopAtomList  for ALL LOOPs
c 
c Yurii Vorobjev, 2002
c
	subroutine inLoopAtomList(addAnchorAtF,
     &           nLoopMX,nLoop,resEndLoop,
     &           natom,atomName,resName,chName,resNumb,
     &           nres,resNameRes,chNameRes,
     &           atomNameEx,startAtInRes,
     &           nAtomInLoop,nAtomInLoopMX,atomInLoopList,
     &           startAtInLoop,stopAtInLoop,loopType,
     &           nAnchorAtom,nAnchorAtomMX,anchorAtomList,
     &           atomXYZ,refPosAnchorAtom)

c
c Input:
c addAnchorAtF : flag 0/1 - add anchor atoms
c nLoop - number of Loops
c resEndLoop(2*nLoop) - resNunber for N,C end of Loop
c  PDB data: natom, etc. 
c      atomXYZ() for all fixed + loop(initial xyz)
c 
c Result: 
c         nAtomInLoop - total number of atoms in Loops (including anchorAt) 
c         nAnchorAtom - achor atoms to adjast the loopEnds to fixed protein
c         atomInLoopList(iloop)=ia, il=1,nAtomInLoop
c         anchorAtomList(ianch)=ia, ianch=1,nAnchorAtom
c         refPosAnchorAtom(ianch) - fixedProt xyz for anchorAtoms
c
c         startAtInLoop(iLoop) - startAtom in the Local numeration 
c                                iaLoop=1,nLoopAtom
c         stopAtInLoop(iLoop) - stopAtom in the Local numeration 
c
 	implicit none
        integer natom
        integer addAnchorAtF
        character*4 atomName(*)
        character*4 resName(*)
        character*1 chName(*)
        integer  resNumb(*)
        integer nres
        character*4 resNameRes(*)
        character*1 chNameRes(*)
        character*8 atomNameEx(*)
        integer startAtInRes(*)
        integer nLoop, nLoopMX
        integer resEndLoop(*)
        integer natomInLoop
        integer natomInLoopMX,atomInLoopList(*)
        integer nAnchorAtom,nAnchorAtomMX
        integer anchorAtomList(*)
        real atomXYZ(*)
        real refPosAnchorAtom(*)
	integer startAtInLoop(*),stopAtInLoop(*)
        integer loopType(*)
c
        include "charStringSiz.h"
        include "output.h"
        include "kanalUse.h"
        include "statusMessg_mDyn.h"
c local
        integer nAnchorAtRes
        character*4 anchorAtomName(3)
	integer ia,i,ir1,ia1,k
        integer ilp,ilp2,ir,ir2
        integer ia3,k3
        integer kanalp
        logical CONTROL
        data anchorAtomName /'N   ','CA  ','C   '/
c add anchor atoms for both sides of loop

c initialize
        kanalp =  kanalRunOut
        CONTROL = .true.
c
        write(kanalp,*)'inLoopAtomList : START '
        write(kanalp,*)
c initialize:
        nAtomInLoop = 0
        nAnchorAtom = 0
        do i =1,nAtomInLoopMX
        atomInLoopList(i) = 0
        end do
        do i = 1,nAnchorAtomMX
        anchorAtomList(i) = 0
        end do
        do i = 1,nLoopMX
        startAtInLoop(i) = 0
        stopAtInLoop(i) = 0
        loopType(i) = 0
        end do
c anchor atom perResidue in Protein
        nAnchorAtRes = 3 
c
c collect moving atoms in Loop
        startAtInLoop(1) = 1
c
        if(nLoop .ge. 1)then
c
        do ilp=1,nLoop
c define loop
        ilp2=2*ilp
c first and last loop residues
        ir1=resEndLoop(ilp2-1)  !start Loop residue
        ir2=resEndLoop(ilp2)    !end Loop residue
c
c anchor res if Loop is Internal
        if(addAnchorAtF .eq. 1)then
        if(ir1 .gt. 1 )then  
c add anchor at N term of the (ir1-1) residue
        do ia1 = startAtInRes(ir1-1),startAtInRes(ir1)-1 
        do k=1,nanchorAtRes
        if(atomName(ia1) .eq. anchorAtomName(k))then
        nAtomInLoop = nAtomInLoop + 1
        nAnchorAtom = nAnchorAtom + 1
        stopAtInLoop(ilp) = nAtomInLoop
        startAtInLoop(ilp+1) = stopAtInLoop(ilp)+1
c
        if(nAnchorAtom .gt. nAnchorAtomMX-2*nanchorAtRes)then
        write(kanalp,*)'ERROR!: nAnchorAtomMX small !',nAnchorAtom
        write(kanalp,*)'increase param nAnchorAtomMAX in loopInfo.h'
        stop
        end if
c
        atomInLoopList(nAtomInLoop) = ia1
        anchorAtomList(nAnchorAtom) = ia1
        ia3=3*ia1
        k3 = 3*nAnchorAtom
        refPosAnchorAtom(k3-2)=atomXYZ(ia3-2)
        refPosAnchorAtom(k3-1)=atomXYZ(ia3-1)
        refPosAnchorAtom(k3)=atomXYZ(ia3)
        end if ! .eq. anchorAtomName(k)
        end do !k 
        end do !ia1
        end if !ir1 .gt. 1
        end if !addAnchorAtF
c
c add Loop residues 
        do ir = resEndLoop(ilp2-1),resEndLoop(ilp2)
c residue in the loop
        do ia = startAtInRes(ir),startAtInRes(ir+1)-1  
c atom in the res ir
        nAtomInLoop = nAtomInLoop + 1
        stopAtInLoop(ilp) = nAtomInLoop
        startAtInLoop(ilp+1) = stopAtInLoop(ilp)+1
c
        if(nAtomInLoop .gt. nAtomInLoopMX-2*nanchorAtRes)then
        write(kanalp,*)'ERROR!:nAtomInLoopMAX small !',nAtomInLoop
        write(kanalp,*)'increase param natomInLoopMAX in loopInfo.h'
        stop
        end if
        atomInLoopList(nAtomInLoop) = ia
        end do !ia 
        end do !ir
c
        if(addAnchorAtF .eq. 1)then
        if(ir2 .lt. nres )then
c add anchor at C term of the (ir2+1) residue
        do ia1 = startAtInRes(ir2+1),startAtInRes(ir2+2)-1 
        do k=1,nanchorAtRes
        if(atomName(ia1) .eq. anchorAtomName(k))then
        nAtomInLoop = nAtomInLoop + 1
        stopAtInLoop(ilp) = nAtomInLoop
        startAtInLoop(ilp+1) = stopAtInLoop(ilp)+1
c
        nAnchorAtom = nAnchorAtom + 1
        atomInLoopList(nAtomInLoop) = ia1
        anchorAtomList(nAnchorAtom) = ia1
        ia3=3*ia1
        k3 = 3*nAnchorAtom
        refPosAnchorAtom(k3-2)=atomXYZ(ia3-2)
        refPosAnchorAtom(k3-1)=atomXYZ(ia3-1)
        refPosAnchorAtom(k3)=atomXYZ(ia3)
        end if
        end do !k 
        end do !ia1
        end if !ir2 .lt. nres
        end if !addAnchorAtF
c
        end do! ilp 
        end if !nLoop .ge. 1
c
        if(CONTROL)then
        write(kanalp,*)'atomInLoopList:nAtomInLoop=',nAtomInLoop
        do i=1,nAtomInLoop
        ia = atomInLoopList(i)
        write(kanalp,'(i4,1x,i4,2x,a4,1x,a4,1x,i4)')
     &  i,ia,atomName(ia),resName(ia),resNumb(ia)
        end do !i
c
         write(kanalp,*)'anchorAtomList:nAnchorAtom=',nAnchorAtom
        do i=1,nAnchorAtom
        ia = anchorAtomList(i)
        write(kanalp,'(i4,1x,i4,2x,a4,1x,a4,1x,i4,3x,3(f8.3))')
     &  i,ia,atomName(ia),resName(ia),resNumb(ia),
     &  (refPosAnchorAtom(3*i-3+k),k=1,3)
        end do !i
        end if !CONTROL
c
        write(kanalp,*)'inLoopAtomList : FINISH '
        write(kanalp,*)
c
	return
	end
