c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c* get Hb12-8 H-bond potential parameters for i,j atoms(glob numb)           *
c* from Hb128 Table                                                          *
c*                                        hbond128.h                         *
c*     Yury Vorobjev 2004                                                    *
c* In: i,j   - global aotm numbers                                           *
c* OUT:                                                                      *
c*      hb128status = .true./.false. hb128Exists/NonExist                    *
c*      rm,em,A12,B12 - hBvdw parameters                                     *
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
	subroutine getHb128ijPar(i,j,rm,em,A12,B12,hb128status)
c
        include 'xyzPDBsize.h'
        include 'hbond128.h'
        include 'xyzPDBinfo.h'
c
        include "output.h"
        include "kanalUse.h"
        include "statusMessg_mDyn.h"
c implicit none
	integer i,j
        real rm,em,A12,B12
        logical hb128status
c
        integer getPointHB128ParTable
c local
        character*2 iFF,jFF
        integer acceptor,donor
        integer h,y
        integer iHxN,iYyN,kShift
        integer kanalp
        logical CONTROL
c
        kanalp = kanalRunOut
c        CONTROL = .true. 
        CONTROL = .false.
c
        iFF = ffAtomName(i) 
        jFF = ffAtomName(j)
c
        hb128status = .false.
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
        else
        iYyN = hB128atomType(j) - nHXhb128
        end if
c
         hb128status = .true. 
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
     &  'getHb128ijPar: i,j,iFF,jFF,rm,em,A12,B12,hb128status:'
        write(kanalp,*) 
     &  '            ',i,j,iFF,jFF,rm,em,A12,B12,hb128status  
        end if !C
c
101     continue
c
	return
	end
