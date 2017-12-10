c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c* get Hb12-8 H-bond triplet X-H...Y for pair H...Y = i,j atoms(glob numb)   *
c* get Hb128angle                                                            *
c* from Hb128 Table                                                          *
c*                                        hbond128.h                         *
c*     Yury Vorobjev 2004                                                    *
c* In: i,j   - global atom numbers to check them on Hb128 formation          *
c* OUT:                                                                      *
c*      tijk(3) = Y,H,X  of X-H...Y triplet, 
c*               tijk(1)=iY, tijk(2)=jH, tijk(3)=kX                          *
c*      hb128status = .true./.false. hb128Exists/NonExist                    *
c*      rm,em,A12,B12 - hBvdw parameters                                     *
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
	subroutine getHb128ijkTripPar
     &             (i,j,tijk,rm,em,A12,B12,hb128status,hB128Type)
c
        include 'xyzPDBsize.h'
c
        include "output.h"
        include "kanalUse.h"
        include "statusMessg_mDyn.h"
c
        include 'hbond128.h'
        include 'xyzPDBinfo.h'
        include 'pair1234array.h'
c implicit none
	integer i,ib,j
        integer tijk(3)
        real rm,em,A12,B12
        logical hb128status
        integer hB128Type
c
        integer getPointHB128ParTable
c local
        character*2 iFF,jFF
        integer acceptor,donor
        integer h,y
        integer iHxN,iYyN,kShift
c back-backbone Hb Only
        logical OPT_backBackHb128
        logical bbX,bbY
        integer nAtProtBBHb128max
        parameter (nAtProtBBHb128max=2)
        character*4 protBBoneHb128(nAtProtBBHb128max)
        integer nAtProtBBHb128
        integer kanalp
        logical CONTROL
c
        data protBBoneHb128/'O   ','N   '/
c                            xxxx
        kanalp = kanalRunOut
c        CONTROL = .true. 
        CONTROL = .false.
        nAtProtBBHb128=nAtProtBBHb128max
c
        OPT_backBackHb128 = .true.
        iFF = ffAtomName(i) 
        jFF = ffAtomName(j)
c
        hb128status = .false.
        hB128Type = 1
c
c define ATHx_N     
        iHxN = 0
        iYyN = 0
        rm = 0.0
        em = 0.0
        A12 = 0.0
        B12 = 0.0
c
        iHxN = hB128atomType(i)
        iYyN = hB128atomType(j)
        if(iHxN .eq. 0 .or. iYyN .eq. 0 ) goto 101
        if(iHxN .le. nHXhb128 .and. iYyN .le. nHXhb128)goto 101  
        if(iHxN .gt. nHXhb128 .and. iYyN .gt. nHXhb128)goto 101
c
        if(iHxN .gt. nHXhb128) then
        iYyN = hB128atomType(i) - nHXhb128
        iHxN = hB128atomType(j)
        tijk(1) = i
        tijk(2) = j
        else
        iYyN = hB128atomType(j) - nHXhb128
        tijk(1) = j
        tijk(2) = i
        end if
c atom H = tijk(2) : X=atom bonded to H
        tijk(3) = pair12List(startPairL12(tijk(2)))
c
         hb128status = .true. 
c
c protBBone Filter:
         if(OPT_backBackHb128)then
         bbY = .false.
         bbX = .false.
         do ib=1,nAtProtBBHb128
         if(atomName(tijk(3)) .eq. protBBoneHb128(ib))bbX = .true.
         if(atomName(tijk(1)) .eq. protBBoneHb128(ib))bbY = .true.
         end do !ib
         if(bbX .and. bbY) then
         hB128Type = 1           ! back-back type
         else
         hB128Type = 2           ! back-side, or side-side
         end if 
c
         end if ! OPT_backBackHb128
c
        kShift =  getPointHB128ParTable(iHxN,iYyN,nYYhb128)
c
        rm = hB128parRE(kShift+1)
        em = hB128parRE(kShift+2)
        A12 = hb128parAB(kShift+1)
        B12 = hb128parAB(kShift+2)
c
c 101     continue
c
        if(CONTROL)then
        write(kanalp,*)
     &  'getHb128ijTripPar: i,j,iFF,jFF,rm,em,A12,B12,hb128status:',
     &  ' X-H..Y:'
        write(kanalp,*) 
     &  '            ',i,j,iFF,jFF,rm,em,A12,B12,hb128status,
     &  tijk(3),tijk(2),tijk(1)
        end if !C
c
101     continue
c
	return
	end
