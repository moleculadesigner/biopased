c get atomList by specified atom names
c
c Vorobjev YN 2005
c
	subroutine getAtomListByName(iv,natom,atomName,atResName,
     &    nSpecAtType,specAtNameList,natomExList,natomExListMX,
     &    atomNumbExList)
c
c iv=0, extract by atom name comparision
c collect atomNumbExList(*) by specified specAtNameList(*) 
c                from full atomName(*) list for all atoms
c
	implicit none
        integer iv
        character*4 atomName(*)
        character*4 atResName(*)
        integer natom,nSpecAtType
        character*4 specAtNameList(*)
        integer natomExList,natomExListMX
        integer atomNumbExList(*)
c
        include "charStringSiz.h"
        include "output.h"
        include "kanalUse.h"
        include "statusMessg_mDyn.h"
c local
        integer ia,is
        logical inList
        logical CONTROL
        integer kanalp
c
        CONTROL=.false.
        kanalp = kanalRunOut
c
c iv = 0, then extraction by atomName only
c
        natomExList=0
        do ia = 1,natom
        inList=.false.
        do is = 1,nSpecAtType
        if(atomName(ia) .eq. specAtNameList(is))then
        inList=.true.
        goto 101
        end if 
        end do!is
101     continue
        if(inList)then
        natomExList=natomExList+1
c
        if(natomExListMX .lt. natomExList)then
         iStatus = iStatus + 1          
        write(kanalPStat,*)mStatus,iStatus,messageStatus
     &  ,' error extracAtList too long ! natomExListMX:',natomExListMX
        write(kanalp,*)
     &  'getAtomListByName: ERROR!: extracAtList too long:',
     &  ' natomExListMX:', natomExListMX
        stop
        end if!error
c
        atomNumbExList(natomExList)=ia
        end if !inList 
        end do!ia
c
        if(CONTROL)then
        write(kanalp,*)'getAtomListByName:natomExList=',natomExList
        do is=1,natomExList
        ia=atomNumbExList(is)
        write(kanalp,'(a16,i5,i5,1x,a4)')
     &                'is,ia,atName:',is,ia,atomName(ia)
        end do!is
        end if!C
c
	return
	end
