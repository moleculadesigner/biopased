c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c*  compactization to GeoCentr energy/forces                                 *
c*                                                                           *
c*     Yury Vorobjev 2004                                                    *
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
        subroutine compactGeoCntEF(atomXYZ,
     &           natom,nmoveatom,moveAtomList,
     &           compactGeoFpar,compactGeoEn,compactGeoForce)
c
c SLOW Forces: USed in : mdRun.f
c compactGeoForceParam = compactGeoFpar(): alfa,d1,beta,d2
c
c calcuted for moveAtomList(im)=ia, im=1,nmoveatom  
c              ri(ir) = atomXYZ(ia)
c
	implicit none
        real atomXYZ(*) 
        integer natom
        integer nmoveatom  
        integer moveAtomList(*)
        real compactGeoFpar(*)
        real compactGeoEn
        real compactGeoForce(*)
c
        include "charStringSiz.h"
        include "output.h"
        include "kanalUse.h"
        include "statusMessg_mDyn.h"
c
        real rgc(3)
        real ri0(3),enr   
        real a,aa,bb,rM0,rMX
        real e0,f0,cf,rri
        real fX,eX,scale
        real fcsum(3)
        integer i,i3,ia,ia3,k
        integer kanalp
        logical CONTROL
c
cx        kanalp = 6
        kanalp = kanalRunOut
c
c        CONTROL = .true.  
        CONTROL = .false.
c
        scale = 1.0/(compactGeoFpar(3)+compactGeoFpar(4))
        scale = scale - 1.0/(compactGeoFpar(3)+compactGeoFpar(2))
        if(scale .gt. 0.0)then
        scale = 1.0/scale
        else 
         scale = 1.0
        end if
c
c initialize
        compactGeoEn=0.0
     	do i=1,3*natom
        compactGeoForce(i) = 0.0
        end do 
c calculate geoCenter
        do i=1,3
        rgc(i)=0.0
        fcsum(i)=0.0
        end do!i
c
        do i=1,natom
        i3 = 3*i-3
        rgc(1)=rgc(1) + atomXYZ(i3+1) 
        rgc(2)=rgc(2) + atomXYZ(i3+2) 
        rgc(3)=rgc(3) + atomXYZ(i3+3) 
        end do!i
c
        rgc(1)=rgc(1)/natom
        rgc(2)=rgc(2)/natom
        rgc(3)=rgc(3)/natom
c max dist from center
        rMX = -1.0e10
        do i=1,natom
        i3 = 3*i-3
        a =   (atomXYZ(i3+1) - rgc(1))**2
     &      + (atomXYZ(i3+2) - rgc(2))**2
     &      + (atomXYZ(i3+3) - rgc(3))**2
        if( a .gt. rMX ) rMX = a
        end do !i
c
        rMX = sqrt(rMX) + compactGeoFpar(4) 
        rM0 = rMX - compactGeoFpar(2)
        if(rM0 .lt. 0.0 ) rM0=0.0
c
        aa = compactGeoFpar(1)
        bb = compactGeoFpar(3)
c e0, f0 - fenergy and force at r=rM0 
        e0 = aa/(bb+compactGeoFpar(2))
        f0 = aa/(bb+compactGeoFpar(2))**2
        eX = aa/bb
        fX = eX/bb
c
        if(CONTROL)then
        write(kanalp,*)' getCompactGeoEF: geoCxyz:',rgc
        write(kanalp,*)' getCompactGeoEF: rMX, rM0:',rMX, rM0
        write(kanalp,*)'compactGeoFpar:aa,d1,bb,d2:',
     &  (compactGeoFpar(k),k=1,4)
c
        write(kanalp,*)
     &  '  i  ia eng: atomXYZ(3)  compactGeoForce:'
        end if
c
        if(nmoveatom .ge. 1) then
c
	do i = 1,nmoveatom
        i3 = 3*i - 3
        ia = moveAtomList(i)
        ia3 = 3*ia - 3
c
        ri0(1) = atomXYZ(ia3+1) - rgc(1)
        ri0(2) = atomXYZ(ia3+2) - rgc(2)
        ri0(3) = atomXYZ(ia3+3) - rgc(3)
        rri = sqrt(ri0(1)**2+ri0(2)**2+ri0(3)**2)
c
        enr = 0.0
        cf = 0.0
c
        if ( rri .ge. rM0 .and. rri .le. rMX)then
        enr = aa/(bb+rMX-rri) - e0 - f0*(rri-rM0)
        cf = aa/(bb+rMX-rri)**2 - f0
        end if ! rri < rM0
c
        if ( rri .ge. rMX)then
        enr = eX - e0 + (fX - f0)*(rri-rMX)
        cf = fX - f0
        end if ! rri > rMX
c
        enr = enr*scale
        compactGeoEn = compactGeoEn + enr
c
        cf = cf*scale
        compactGeoForce(ia3+1) = - cf*ri0(1)
        compactGeoForce(ia3+2) = - cf*ri0(2)
        compactGeoForce(ia3+3) = - cf*ri0(3)
c
        fcsum(1) = fcsum(1) + compactGeoForce(ia3+1)
        fcsum(2) = fcsum(2) + compactGeoForce(ia3+2)
        fcsum(3) = fcsum(3) + compactGeoForce(ia3+3)
c
        if(CONTROL)then
        write(kanalp,'(2i5,1x,f5.2,1x,f6.1,2x,3(3f7.2,2x))')
     &  i,ia,enr,rri,
     &  (atomXYZ(ia3+k),k=1,3),(compactGeoForce(ia3+k),k=1,3)
        end if !CONTROL
c
        end do !i
c
        if(CONTROL)then 
        write(kanalp,'(a10,3f8.3)')'fcsum : ',fcsum
        end if !C
c correct forces
        fcsum(1)=fcsum(1)/nmoveatom
        fcsum(2)=fcsum(2)/nmoveatom
        fcsum(3)=fcsum(3)/nmoveatom
c

        if(CONTROL)then
        write(kanalp,'(a33)')'Corrected fCompact:'
        end if !C
c
        do i = 1,nmoveatom
        i3 = 3*i - 3
        ia = moveAtomList(i)
        ia3 = 3*ia - 3
        compactGeoForce(ia3+1) = compactGeoForce(ia3+1)-fcsum(1) 
        compactGeoForce(ia3+2) = compactGeoForce(ia3+2)-fcsum(2) 
        compactGeoForce(ia3+3) = compactGeoForce(ia3+3)-fcsum(3) 
c
        if(CONTROL)then
        write(kanalp,'(2i5,1x,f5.2,1x,f6.1,2x,3(3f7.2,2x))')
     &  i,ia,enr,rri,
     &  (atomXYZ(ia3+k),k=1,3),(compactGeoForce(ia3+k),k=1,3)
        end if !CONTROL
c
        end do ! i
c
        end if ! nmoveatom .ge. 1
c
        if(CONTROL) STOP
c
	return
	end
