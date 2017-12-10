c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c*  refMolPosRestrain1RHW (harmonic well) energy/forces                      *
c*  harmonic well for ResidueCMass move                                      *
c*     Yury Vorobjev 2005                                                    *
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
        subroutine refResPosRestr1RHWEF(atomXYZ,natom,
     &           atomMass,startAtInRes,stopAtInRes,
     &           nRestr1RHWAtom,restr1RHWatomList,
     &           nRestr1ResidHW,restr1ResidHWList,
     &           restr1RefResidPosHW,restr1ResidMass,
     &           restr1RHWConst,sizeRestr1RHW,
     &           restr1RHWEng,restr1RHWAtForce)
c
c wall harmonic restaraint1w (type1HW) 
c erestr1HW = K*[(rCM - riref)**2 - size1HW]**2
c
c calcuted for restr1HWatomList(ir)=ia, ir=1,nRestr1atom
c              ri(ir) = atomXYZ(ia)
c              refAtomPos(ir), ir=1,nRestr1atom
c
	implicit none
        real atomXYZ(*) 
        integer natom
        real atomMass(*)
        integer startAtInRes(*),stopAtInRes(*)
        integer nRestr1RHWAtom
        integer restr1RHWatomList(*)
        integer nRestr1ResidHW
        integer restr1ResidHWList(*)
        real restr1RefResidPosHW(*)
        real restr1ResidMass(*)
        real restr1RHWConst
        real sizeRestr1RHW 
        real restr1RHWEng
        real restr1RHWAtForce(*)
c
        include "charStringSiz.h"
        include "output.h"
c local
        real ri0(3),enr  
        real rcm(3)
        integer i,i3,ia,ia3
        integer ir,k,iL,iL3
        real sz1HW2,en1,d1,dsh2,ms
        integer kanalp
        logical CONTROL
c
        kanalp =  kanalRunOut
cx       CONTROL = .true.  
        CONTROL = .false.
c initialize
        restr1RHWEng = 0.0
     	do i=1,3*natom
        restr1RHWAtForce(i)= 0.0
        end do 
c
        sz1HW2 = sizeRestr1RHW**2 
c
        if(CONTROL)then
        write(kanalp,*)'refPosRestr1RHWenForce:'
        write(kanalp,*)
     &  'nRestr1ResidHW, restr1RHWConst,sizeRestr1RHW:',
     &  nRestr1ResidHW, restr1RHWConst,sizeRestr1RHW
        write(kanalp,*)(restr1ResidMass(i),i=1,nRestr1ResidHW)
        end if!C
c
        if(nRestr1ResidHW .ge. 1) then
c rcm()
        do iL = 1,nRestr1ResidHW
        ir = restr1ResidHWList(iL)
        iL3 = 3*iL - 3
        rcm(1) = 0.0
        rcm(2) = 0.0
        rcm(3) = 0.0
        do ia = startAtInRes(ir),stopAtInRes(ir)
        ia3 = 3*ia - 3
        ms = atomMass(ia)/restr1ResidMass(iL)
        rcm(1) = rcm(1) + atomXYZ(ia3+1)*ms            
        rcm(2) = rcm(2) + atomXYZ(ia3+2)*ms                
        rcm(3) = rcm(3) + atomXYZ(ia3+3)*ms                
        end do !ia
c        
        ri0(1) = rcm(1) - restr1RefResidPosHW(iL3+1)
        ri0(2) = rcm(2) - restr1RefResidPosHW(iL3+2)
        ri0(3) = rcm(3) - restr1RefResidPosHW(iL3+3)
        dsh2 = ri0(1)**2 + ri0(2)**2 + ri0(3)**2
c
        enr = 0.0
        en1 = 0.0
        if(dsh2 .gt. sz1HW2)then
        d1 = dsh2 - sz1HW2
        en1 = restr1RHWConst*d1 
        enr = en1*d1
        end if
        restr1RHWEng =  enr
c
        do ia = startAtInRes(ir),stopAtInRes(ir)
        ia3 = 3*ia - 3
c
        ms = -4.0*en1*atomMass(ia)/restr1ResidMass(iL)
        restr1RHWAtForce(ia3+1) =  ms*ri0(1)
        restr1RHWAtForce(ia3+2) =  ms*ri0(2)
        restr1RHWAtForce(ia3+3) =  ms*ri0(3)
c
        if(CONTROL)then
        write(kanalp,'(2i5,1x,f10.2,2x,3(3f7.2,2x),f7.3)')
     &  iL,ia,restr1RHWEng,
     &  (rcm(k),k=1,3),(restr1RefResidPosHW(iL3+k),k=1,3),
     &  (restr1RHWAtForce(ia3+k),k=1,3)
        end if !CONTROL
c
        end do !ia
        end do !iL
c
        end if ! .ge. 1
c
	return
	end
c
