c READ in PDBXYZ 
c atmXYZ() real*8 to USE in add_h  subroutine
c
        subroutine readPDB05_nh
c
c copy xyzPDBinfo.h to addHxyz.h                                    
c
c        implicit none
        include 'xyzPDBsize.h'
        include 'xyzPDBinfo.h'
        include 'xyzPDBcrd.h'
        include 'addHxyz.h'
c
        include "output.h"
        include "kanalUse.h"
        include "statusMessg_mDyn.h"
c
        integer ia,i3,ir  
        integer k,k3,j
        integer kanalp
        logical CONTROL
c
        kanalp=kanalRunOut
        CONTROL=.true.
c
        natom_nh=natom
        do ia =1,natom 
        atomNumb_nh(ia)=atomNumb(ia)
        atomName_nh(ia)=atomName(ia)
        resName_nh(ia) =resName(ia)
        chName_nh(ia) = chName(ia)
        resNumb_nh(ia) = resNumb(ia)
        atomNameEx_nh(ia)=atomNameEx(ia)
        resNameEx_nh(ia)=resNameEx(ia)
        i3=3*ia-2
        atomXYZ_nh(1,ia)=atomXYZ(i3)
        atomXYZ_nh(2,ia)=atomXYZ(i3+1)
        atomXYZ_nh(3,ia)=atomXYZ(i3+2)
        end do !ia
c
        nres_nh=nres
        do ir=1,nres
        resNameRes_nh(ir)=resNameRes(ir)
        chNameRes_nh(ir)=chNameRes(ir)
        startAtInRes_nh(ir)=startAtInRes(ir)
        stopAtInRes_nh(ir)=stopAtInRes(ir)
        realAtomFlag_nh(ir)=realAtomFlag(ir)
        atHvyNbInRes_nh(ir)=atHvyNbInRes(ir)
        end do !ir
c
        if(CONTROL)then
           write(kanalp,*)'PDBread05_nh:pdbXYZInPut: Natom:',natom
           do k=1,natom
           k3=3*k-3
c
           write(kanalp,7072)
     &     k, 'ATOM  ',atomNumb_nh(k),atomName_nh(k),
     &     resName_nh(k),chName_nh(k),
     &     resNumb_nh(k), (atomXYZ_nh(j,k),j=1,3)
     &     ,atomNameEx_nh(k),resNameEx_nh(k)
           end do !k
           write(kanalp,*)
     &  'pdbRead05_nh:strtAtInRes(k) stpAtInRes(k) resTyp atHvyNbInRes'
           do k=1,nres
           if(startAtInRes(k) .gt. 0)then
           write(kanalp,'(i6,a6,a4,i5,1x,i5,1x,i5)')
     &     k,resNameRes_nh(k),chNameRes_nh(k),
     &     startAtInRes_nh(k),stopAtInRes_nh(k),atHvyNbInRes_nh(k)
           else
           write(kanalp,'(i6,a6,a4,i5,1x,i5)')
     &     k,'      ','    ',startAtInRes(k),stopAtInRes(k)
           end if!startAtInRes(k) .gt. 0 !nonLoop
           end do!k
           end if !C
        return
7071    format(a6,i5,2x,a4,a4,a1,i4,4x,3f8.3) !orig PDB
7072    format(i5,1x,a6,1x,i4,2x,a4,a4,a1,i4,4x,3f8.3,2x,a8,1x,a8)
        end
