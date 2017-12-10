c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c* get getHb128AtType : create file hB128atomType(*)  in  hbond128.h         *
c* from Hb128 Table                                                          *
c*                                                                           *
c*     Yury Vorobjev 2004                                                    *
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
	subroutine getHb128AtType
c 
c defines hB128atomType(ia) = iH128  = hB128atomType(Code)  = H=ATHX_L;
c                           Acceptor =(ATYY_L+ATHX_N)
c                           H=ATHX_L, ATYY_L : parameters in h128ParamFile
c
        include 'xyzPDBsize.h'
c
        include "output.h"
        include "kanalUse.h"
        include "statusMessg_mDyn.h"
c
        include 'hbond128.h'
        include 'xyzPDBinfo.h'
c implicit none
	integer ia
c
        character*2 iFFname
        integer h,y
        integer iH128    
        integer kanalp
        logical CONTROL
c
        kanalp = kanalRunOut
        CONTROL = .true.
c
c loop over all atoms
c
        if(CONTROL)then
        write(kanalp,*)'getHb128atomType:'
        end if
c
        do ia = 1,natom
        iFFname = ffAtomName(ia) 
c define ATHx_N     
        iH128 = 0
c
        do h=1,nHXhb128
        if(hb128HxList(h) .eq. iFFname) iH128 = h
        end do !ih
c
        if(iH128 .eq. 0)then
        do y=1,nYYhb128 
        if(hb128YYList(y) .eq. iFFname) iH128 = y + nHXhb128
        end do !y
        end if ! iH128 .eq. 0
c
        hB128atomType(ia) = iH128
c
        if(CONTROL)then
        write(kanalp,*) 'getHb128atomType:',
     &  ia,atomName(ia),resName(ia),ffAtomName(ia),hB128atomType(ia)
        end if !C
c
        end do !ia
c
	return
	end
