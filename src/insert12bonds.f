c insert12bonds.f
c
c Yuri Vorobjev 2005
c
	subroutine insert12bond
     &          (ia,ja,natom,nPairL12,startPairL12,pair12List)
c
c insert new bond 12 ia-ja into nPairL12(ia),startPairL12(ia),pair12List(*)
c                         data structure
c natom - number of atoms
c nPairL12(ia) - number of 12 bonds for ia atom
c startPairL12(ia) - pointer to the first 12 bond for atom=ia 
c in the list of 12 bonds pair12List(*) 
c pair12List(*) - list of 12 bonds symmetrical, i.e. i-j and j-i
c
        implicit none
c
        include "charStringSiz.h"
        include "output.h"
        include "kanalUse.h"
        include "statusMessg_mDyn.h"
c
        integer ia,ja,natom
        integer nPairL12(*)
        integer startPairL12(*)
        integer pair12List(*)
c
        integer i,ib,ia1,ia2,shift,k
        logical CONTROL
        integer kanalp
c
        CONTROL = .false.
        kanalp =  kanalRunOut 
c 
        if(CONTROL)then
        write(kanalp,*)'start  insert12bonds: '
        do i=1,natom
        write(kanalp,'(a17,3i5,6i5)')'ia,nP12,st,pL12:',
     &  i,nPairL12(i),startPairL12(i),
     &  (pair12List(k),k=startPairL12(i),startPairL12(i)+nPairL12(i)-1)
        end do!i
        end if!C
c
c update nPairL12
        nPairL12(ia) = nPairL12(ia) + 1
        nPairL12(ja) = nPairL12(ja) + 1
c
        ia1=ia
        ia2=ja
        if(ia .gt. ja)then
        ia1=ja
        ia2=ia
        end if
c
        do i = ia1+1,ia2 
        startPairL12(i) = startPairL12(i) + 1
        end do !i
c
        do i = ia2+1,natom
        startPairL12(i) = startPairL12(i) + 2
        end do!i
c
c make a shift
        do i = natom,ia1+1,-1
        if(i .gt. ia2)then
        shift=2
        else 
        shift = 1
        end if
c
        do ib = startPairL12(i)+nPairL12(i)-1,startPairL12(i),-1
        pair12List(ib) =  pair12List(ib-shift)
        end do!ib 
        end do !i
c insert new
        pair12List(startPairL12(ia1+1)-1) = ia2
        pair12List(startPairL12(ia2+1)-1) = ia1
c
        if(CONTROL)then
        write(kanalp,*)'done  insert12bonds: '
        do i=1,natom
        write(kanalp,'(a17,3i5,6i5)')'ia,nP12,st,pL12:',
     &  i,nPairL12(i),startPairL12(i),
     &  (pair12List(k),k=startPairL12(i),startPairL12(i)+nPairL12(i)-1)
        end do!i 
        end if!C
c
	return
	end
