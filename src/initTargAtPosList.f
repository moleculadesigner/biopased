c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c                                                                   *
c  Yury Vorobjev 2005                                               *
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c defines:
c
	subroutine initTargAtPosList(natomTarg,atomNumbTarg,atomNameTarg,
     &             resNameTarg,resNumbTarg,atomXYZTarg,
     &             nRestr1atom,restr1atomList,refAtomPos,
     &             restr1AtConst,restr1AtConstList)   
c
        include 'xyzPDBsize.h'
        include 'xyzPDBinfo.h'
        include 'xyzPDBcrd.h'
c
        include "output.h"
        include "kanalUse.h"
        include "statusMessg_mDyn.h"
c
        integer natomTarg
        character*4 atomNameTarg(*)
        character*4 resNameTarg(*)
        integer resNumbTarg(*),atomNumbTarg(*)
        real atomXYZTarg(*)
c
        integer restr1atomList(*)
        real refAtomPos(*)
        integer nRestr1atom
        real restr1AtConst
        real restr1AtConstList(*)
clocal
        integer kanalp
        integer i,ia,k
        integer i3,ir3
        integer resN
        logical includeInA1List
        logical CONTROL,CONTROL1
c        
c initialize:
        kanalp =  kanalRunOut
        CONTROL = .true.
        CONTROL1 = .false.
c
        if(CONTROL)then
        write(kanalp,*)
     &  'initTargAtPosList: define initRestr1AtList ...',
     &  ' natomTarg = ',natomTarg
        end if
c
c define refAtomPos()
        do i = 1,natomTarg
c
        resN = resNumbTarg(i)
c
        if(CONTROL1)then
        write(kanalp,*)'initTargAtPosList: resNumbTarg(i):',resN
        write(kanalp,7071)
     &  'ATOM  ',i,atomNameTarg(i),resNameTarg(i),
     &  resNumbTarg(i) 
        write(kanalp,*)'startAtInRes(resN),stopAtInRes(resN):',
     &  startAtInRes(resN),stopAtInRes(resN)
        end if
c
        do ia = startAtInRes(resN),stopAtInRes(resN)
c
        if( (atomNameTarg(i) .eq. atomName(ia)) .and.
     &      (resNameTarg(i) .eq. resName(ia)) )then
c
         nRestr1atom = nRestr1atom + 1
         restr1atomList(nRestr1atom) = ia
        i3 = 3*i - 3
        ir3 = 3*nRestr1atom - 3
        refAtomPos(ir3+1) = atomXYZTarg(i3+1)
        refAtomPos(ir3+2) = atomXYZTarg(i3+2)
        refAtomPos(ir3+3) = atomXYZTarg(i3+3)
        restr1AtConstList(i)=restr1AtConst    
        end if ! atom=atomTarg
        end do !ia
        end do !i
c
        if(CONTROL)then
        write(kanalp,*)'initTargAtPosList:'
        write(kanalp,*)'nRestrTarg1atom :',nRestr1atom
        do i = 1,nRestr1atom
        ia=restr1atomList(i)
        ir3=i*3-3
        write(kanalp,'(a6,i5,2x,a4,a4,1x,i4,4x,3f8.3,2x,1x,i4,f6.3)')
     &  'ATOMtg',ia,atomName(ia),resName(ia),resNumb(ia),
     &  (refAtomPos(ir3+k),k=1,3),
     &  i,restr1AtConstList(i)
        end do !i
        end if !CONTROL
c
        write(kanalp,*)'initTargAtPosList: Finish:'
c
	return
7071    format(a6,i5,2x,a4,a4,1x,i4,4x,3f8.3) !orig PDB 
        end
