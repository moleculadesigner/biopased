c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c*  refAtPosRestrain (harmonic) energy/forces                                *
c*                                                                           *
c*     Yury Vorobjev 2002                                                    *
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
        subroutine refAtPosRestr1EF(atomXYZ,
     &           natom,nRestr1atom,restr1atomList,
     &           refAtomPos,restrConstList,restr1Eng,restr1AtForce)
c
c harmonic restaraint1 (type1) 
c erestr1 = K*(ri - riref)**2
c
c calcuted for restr1atomList(ir)=ia, ir=1,nRestr1atom
c              ri(ir) = atomXYZ(ia)
c              refAtomPos(ir), ir=1,nRestr1atom
c
	implicit none
        real atomXYZ(*) 
        integer natom
        integer nRestr1atom
        integer restr1atomList(*)
        real refAtomPos(*)
        real restrConstList(*)
        real restr1Eng
        real restr1AtForce(*)
c
        include "charStringSiz.h"
        include "output.h" 
        real restrConst
        real ri0(3),enr   
        integer i,i3,ia,ia3,k
        integer kanalp
        logical CONTROL
c
        kanalp = kanalRunOut
        CONTROL = .false.
c initialize
        restr1Eng = 0.0
     	do i=1,3*natom
        restr1AtForce(i) = 0.0
        end do 
c
        if(CONTROL)then
        write(kanalp,*)' refAtPosRestr1EF:'
        write(kanalp,*)
     &  '  i  ia eng: atomXYZ(3) refPos(3)  restr1AtForce:'
        end if
c
        if(nRestr1atom .ge. 1) then
c
	do i = 1,nRestr1atom
        i3 = 3*i - 3
        ia = restr1atomList(i)
        ia3 = 3*ia - 3
        restrConst=restrConstList(i)
c
        ri0(1) = atomXYZ(ia3+1) - refAtomPos(i3+1)
        ri0(2) = atomXYZ(ia3+2) - refAtomPos(i3+2)
        ri0(3) = atomXYZ(ia3+3) - refAtomPos(i3+3)
c
        enr = restrConst*(ri0(1)**2 + ri0(2)**2 + ri0(3)**2)
        restr1Eng = restr1Eng + enr
        restr1AtForce(ia3+1) = - 2.0*restrConst*ri0(1)
        restr1AtForce(ia3+2) = - 2.0*restrConst*ri0(2)
        restr1AtForce(ia3+3) = - 2.0*restrConst*ri0(3)
c
        if(CONTROL)then
        write(kanalp,'(2i5,1x,f7.2,2x,3(3f7.2,2x))')i,ia,enr,
     &  (atomXYZ(ia3+k),k=1,3),(refAtomPos(i3+k),k=1,3),
     &  (restr1AtForce(ia3+k),k=1,3)
        end if !CONTROL
c
        end do !i
c
        end if ! .ge. 1

	return
	end
c
