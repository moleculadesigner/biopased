c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c*  refMolPosRestrain1HW (harmonic well) energy/forces                       *
c*  harmonic well for rCM move of 1 Molecule=restr1HWatomList(*)             *
c*     Yury Vorobjev 2003                                                    *
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
        subroutine refMolPosRestr1HWEF(atomXYZ,natom,
     &           nRestr1HWatom,restr1HWatomList,atMs1HWatomList,
     &           refMol1HWPos,restr1HWConst,size1HW,
     &           restr1HWEng,restr1HWAtForce)
c
c wall harmonic restaraint1w (type1HW) for the set of atoms=1 Molecule
c erestr1HW = K*[(rCM - riref)**2 - size1HW]**2
c
c calcuted for restr1HWatomList(ir)=ia, ir=1,nRestr1atom
c              ri(ir) = atomXYZ(ia)
c              refAtomPos(ir), ir=1,nRestr1atom
c
	implicit none
        real atomXYZ(*) 
        integer natom
        integer nRestr1HWatom
        integer restr1HWatomList(*)
        real atMs1HWatomList(*)     
        real refMol1HWPos(3)
        real restr1HWConst
        real size1HW
        real restr1HWEng
        real restr1HWAtForce(*)
c
        include "charStringSiz.h"
        include "output.h"
c local
        real ri0(3),enr  
        real rcm(3)
        integer i,i3,ia,ia3,k
        real sz1HW2,en1,d1,dsh2,ms
        integer kanalp
        logical CONTROL
c
        kanalp = kanalRunOut
c        CONTROL = .true.  
         CONTROL = .false.
c initialize
        restr1HWEng = 0.0
        rcm(1) = 0.0
        rcm(2) = 0.0
        rcm(3) = 0.0
     	do i=1,3*natom
        restr1HWAtForce(i) = 0.0
        end do 
c
        sz1HW2 = size1HW**2 
c
        if(CONTROL)then
        write(kanalp,*)' refPosRestr1MHWenForce:'
        write(kanalp,*)'nRestr1HWatom,restr1HWConst,size1HW:',
     &  nRestr1HWatom,restr1HWConst,size1HW 
        write(kanalp,*)(atMs1HWatomList(i),i=1,nRestr1HWatom)
c
        write(kanalp,*)
     &  '  i  ia eng: atomXYZ(3) refPos1MHW(3)  restr1HWAtForce Ms:'
        end if
c
        if(nRestr1HWatom .ge. 1) then
c rcm()
        do i = 1,nRestr1HWatom
        i3 = 3*i - 3
        ia = restr1HWatomList(i)
        ia3 = 3*ia - 3
        rcm(1) = rcm(1) + atomXYZ(ia3+1)*atMs1HWatomList(i)
        rcm(2) = rcm(2) + atomXYZ(ia3+2)*atMs1HWatomList(i)
        rcm(3) = rcm(3) + atomXYZ(ia3+3)*atMs1HWatomList(i)
        end do !i
c        
        ri0(1) = rcm(1) - refMol1HWPos(1)
        ri0(2) = rcm(2) - refMol1HWPos(2)
        ri0(3) = rcm(3) - refMol1HWPos(3)
        dsh2 = ri0(1)**2 + ri0(2)**2 + ri0(3)**2
c
        enr = 0.0
        en1 = 0.0
        if(dsh2 .gt. sz1HW2)then
        d1 = dsh2 - sz1HW2
        en1 = restr1HWConst*d1 
        enr = en1*d1
        end if
        restr1HWEng =  enr
c
        do i = 1,nRestr1HWatom
        ia = restr1HWatomList(i)
        ia3 = 3*ia - 3
c
        ms = -4.0*atMs1HWatomList(i)
        restr1HWAtForce(ia3+1) =  ms*en1*ri0(1)
        restr1HWAtForce(ia3+2) =  ms*en1*ri0(2)
        restr1HWAtForce(ia3+3) =  ms*en1*ri0(3)
c
        if(CONTROL)then
        write(kanalp,'(2i5,1x,f10.2,2x,3(3f7.2,2x),f7.3)')
     &  i,ia,restr1HWEng,
     &  (atomXYZ(ia3+k),k=1,3),(refMol1HWPos(k),k=1,3),
     &  (restr1HWAtForce(ia3+k),k=1,3), -0.25*ms
        end if !CONTROL
c
        end do !i
c
        end if ! .ge. 1
c
	return
	end
c
