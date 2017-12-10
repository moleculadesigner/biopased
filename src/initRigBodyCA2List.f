c initialization of CA-CA restAtomA2List to keep RigBody FixGeom
c
       subroutine initRigBodyCA2List(caRigBodyHK,nRigBody,
     &                          startAtInRGB,stopAtInRGB,
     &                          atomName,atomXYZ,
     &                          nRestrDistA2MX,nRestrDistA2,
     &                 restrDistA2List,restrDistA2DistHK)   
c
c
        implicit none
        real caRigBodyHK
        integer nRigBody
        integer startAtInRGB(*),stopAtInRGB(*) 
        character*4 atomName(*)
        real atomXYZ(*)
        integer nRestrDistA2MX,nRestrDistA2
        integer restrDistA2List(*)
        real restrDistA2DistHK(*)
c
        include "charStringSiz.h"
        include "output.h"
        include "kanalUse.h"
        include "statusMessg_mDyn.h"
c
        real d
        integer i,i2,ia1,ia2,ia13,ia23 
        integer nDistA2in,nDistA2out
        integer ir,nA1,nA2,nDistA2forAt
        integer kanalp
        logical CONTROL
c initialize:
        kanalp = kanalRunOut
        CONTROL = .true.
        nDistA2forAt = 10
c
         if(nRigBody .ge. 1)then 
         do ir = 1,nRigBody
c make 
         nA1=startAtInRGB(ir)
         nA2=stopAtInRGB(ir)
c
         nDistA2in = nRestrDistA2
         call makeAt2DistList(nDistA2forAt,nA1,nA2,
     &        nRestrDistA2,nRestrDistA2MX,restrDistA2List)
         nDistA2out = nRestrDistA2
c 
         if(CONTROL)then
         write(kanalp,*)'initRigBodyCA2Lis: iRigB:',ir,
     &   ' nA1,na2:',nA1,nA2,' nRestrDistA2;',nRestrDistA2
         end if !C
c define restrDistA2DistHK(*)
         do i=nDistA2in+1,nDistA2out
         i2 = i*2-1
         ia1=restrDistA2List(i2)
         ia2=restrDistA2List(i2+1)
         ia13=ia1*3-2
         ia23=ia2*3-2
         d=sqrt((atomXYZ(ia13)-atomXYZ(ia23))**2 + 
     &     (atomXYZ(ia13+1)-atomXYZ(ia23+1))**2 +
     &     (atomXYZ(ia13+2)-atomXYZ(ia23+2))**2)
c          
         restrDistA2DistHK(i2) = d
         restrDistA2DistHK(i2+1) = caRigBodyHK
c
        if(CONTROL)then
        write(kanalp,'(a20,2x,a4,1x,i5,2x,a4,1x,i5,2x,f7.2,f7.2)')
     &  'atA2List(): i1,i2:',
     &  atomName(ia1),restrDistA2List(i2),
     &  atomName(ia2),restrDistA2List(i2+1),  
     &  restrDistA2DistHK(i2),restrDistA2DistHK(i2+1)
        end if !COntrol
c
         end do !i
         end do !irigB
         end if ! nRigBody .ge. 1
c
        return
        end
c
        subroutine makeAt2DistList(nDist,nA1,nA2,
     &     natA2List,natA2ListMX,atA2List)
c
c IN:
c     nDist,nA1,nA2
c     natA2List,natA2ListMX
c     atA2List(*)
c
c OUT: atA2List(*)
c      natA2List
c
cX        implicit none
        include 'xyzPDBsize.h'
c
        integer nDist    ! nDistPerAtom
        integer nA1,nA2  ! distances for atom  i>na1, i<na2
        integer natA2List,natA2ListMX
        integer atA2List(*)
c
        include 'xyzPDBinfo.h'
c 
        include "output.h"
        include "kanalUse.h"
        include "statusMessg_mDyn.h"
c
        integer nAtSet
        integer nAtSetMX
        integer nNameSetMX
        parameter (nNameSetMX = 3)
        parameter (nAtSetMX = nresMAX*nNameSetMX)
        character*4 atNameSet(nNameSetMX)
        integer atInSetList(nAtSetMX)
cX        data atNameSet / 'H   ','CA  ','O   '/
c
        integer nDist1OPT
        integer aAtSet
        integer i,j,k
        integer i1,i2,j2,k1,ia1,ia2
        integer n2,nBLk
        integer nA2L
        logical CONTROL
        integer kanalp
c
        kanalp = kanalRunOut 
c
        CONTROL = .true. 
c       define natNameSet()
        atNameSet(2) = 'CA  '
        atNameSet(1) = 'H   '
        atNameSet(3) = 'O   '
c
        nA2L=0
        nDist1OPT = 3   ! n shortDist
        nAtSet = 0
        do i = nA1,nA2
        do j=1,nNameSetMX
c
        if(atomName(i) .eq. atNameSet(j))then
        nAtSet = nAtSet + 1
        if(nAtSet .gt. nAtSetMX) then
        write(kanalp,*) 'makeAt2DistList: ERROR! nAtSet.gt.nAtSetMX',
     &  ' nAtSet:',nAtSet
        stop
        end if
c
        atInSetList(nAtSet) = i 
        end if !atomName(i) .eq. atNameSet(j)
        end do!j         
        end do !i
c
        if(CONTROL)then
        write(kanalp,*)'initRigBodyCA2List: atInSetList(*):'
        do i=1,nAtSet 
        write(kanalp,*)'i,atInSetList,atName:',
     &  i,atInSetList(i),atomName(atInSetList(i))
        end do !i
        end if!C
c
        if(nAtSet .gt. 1) then
        if(nDist1OPT .gt. nAtSet-1)nDist1OPT = nAtSet-1 
        nBLk = nAtSet/nDist   ! atomInBlockSize per One A-A dist
        do i1 = 1,nAtSet
        ia1 = atInSetList(i1)
        do i2 = 1,nDist
        if(i2 .eq. 1)then
        do k1=1, nDist1OPT               !do shortDist
        j2=i1 + k1
        if(j2 .gt. nAtSet) j2=j2 - nAtSet
        ia2 = atInSetList(j2)
        if(ia1 .ne. ia2)then
        nA2L=nA2L + 1
        natA2List = natA2List + 1
c
        if(natA2List .gt. natA2ListMX)then
        write(kanalp,*)
     &  ' initRigBodyCA2List: ERROR! restrDistA2.h: ', 
     &  ' natA2List .gt. natA2ListMX !'
         write(kanalp,*)'natA2List,natA2ListMX:',natA2List,natA2ListMX
c
         write(kanalPStat,*)mError,
     &   ' initRigBodyCA2List: ERROR! restrDistA2.h: ',
     &   ' natA2ListMAX is too small '
c
        stop
        end if !natA2List .gt. natA2ListMX
c
        n2= natA2List*2-1
        atA2List(n2) = ia1
        atA2List(n2+1) = ia2
        end if !ia1 .ne. ia2
        end do !k1
        end if !i2 .eq. 1       
c
       j2 = i1 + (i2-1)*nBLk
       if(j2 .gt. nAtSet)j2 = j2 - nAtSet 
       ia2 = atInSetList(j2)
       if(ia1 .ne. ia2)then
        nA2L=nA2L + 1
        natA2List = natA2List + 1
c
        if(natA2List .gt. natA2ListMX)then
        write(kanalp,*)
     &  'initRigBodyCA2List: ERROR! natA2List .gt. natA2ListMX !!'
         write(kanalp,*)'natA2List,natA2ListMX:',natA2List,natA2ListMX 
        stop
        end if !natA2List .gt. natA2ListMX
c
        n2 = natA2List*2-1
        atA2List(n2) = ia1
        atA2List(n2+1) = ia2
        end if !ia1 .ne. ia2
        end do !i2
        end do !i1 
c
        end if !nAtSet .gt. 1
c
        if(CONTROL)then
        write(kanalp,*)'initRigBodyCA2List: nA2L:',nA2L,
     &  ' atA2List(*):ia1,ia2:'
 
        do i=natA2List-nA2L+1,natA2List
        ia1=atA2List(2*i-1)
        ia2=atA2List(2*i)  
        write(kanalp,'(a20,2x,a4,1x,i5,2x,a4,1x,i5)')
     &  'atA2List(): i1,i2:', atomName(ia1),ia1,atomName(ia2),ia2
        end do !i 
        end if !COntrol
c
        return
        end
c
