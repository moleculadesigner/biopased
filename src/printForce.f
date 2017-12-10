c utility subr
c control forces
c print all forces
c
	subroutine printForces
c
        include 'xyzPDBsize.h'
        include 'xyzPDBinfo.h'
        include 'enForce.h'
        include 'mdAtomXYZvel.h'
c
        include "output.h"
        include "kanalUse.h"
        include "statusMessg_mDyn.h"
c
	integer i,k
        real a
        integer kanalp
        logical CONTROLP,CONTROLV
c
        kanalp = kanalRunOut
        CONTROLP = .false.
        CONTROLV = .false.
c
        write(kanalp,*)'Print-ControlForces: '
c
        write(kanalp,*)'vbdefForce:CONTROL:'
        do i=1,3*natom
        a=vbdefForce(i)
        CONTROLV = .false.
        if(a .ge. 0.0 .or. a .lt. 0.0)CONTROLV = .true.
        end do!i
c
        if(.not.CONTROLV)then
        write(kanalp,*)'Print Forces:'
        write(kanalp,*)'vbdefForce:'
        do i=1,natom
        write(kanalp,*)'ATOM  ',i,(vbdefForce(3*i-3+k),k=1,3)
        end do !i
        end if
c
        write(kanalp,*)'vAngdefForce: CONTROL:'
        do i=1,3*natom
        a=vAngdefForce(i)
        CONTROLV = .false.
        if(a .ge. 0.0 .or. a .lt. 0.0)CONTROLV = .true.
        end do!i
c
        if(.not.CONTROLV)then
        do i=1,natom 
        write(kanalp,*)'ATOM  ',i,(vAngdefForce(3*i-3+k),k=1,3) 
        end do !i
        end if
c
        write(kanalp,*)'impDefForce: CONTROL:'
        do i=1,3*natom
        a=impDefForce(i)
        CONTROLV = .false.
        if(a .ge. 0.0 .or. a .lt. 0.0)CONTROLV = .true.
        end do!i
        if(.not.CONTROLV)then
        write(kanalp,*)'impDefForce:'
        do i=1,natom
        write(kanalp,*)'ATOM  ',i,(impDefForce(3*i-3+k),k=1,3)
        end do !i
        end if
c
        write(kanalp,*)'torsAngForce: CONTROL:'
        do i=1,3*natom
        a=torsAngForce(i)
        CONTROLV = .false.
        if(a .ge. 0.0 .or. a .lt. 0.0)CONTROLV = .true.
        end do!i
        if(.not.CONTROLV)then 
        write(kanalp,*)'torsAngForce:'
        do i=1,natom
        write(kanalp,*)'ATOM  ',i,(torsAngForce(3*i-3+k),k=1,3)
        end do !i
        end if
c
        write(kanalp,*)'vdwForceR1: CONTROL:'
        do i=1,3*natom
        a=vdwForceR1(i)
        CONTROLV = .false.
        if(a .ge. 0.0 .or. a .lt. 0.0)CONTROLV = .true.
        end do!i
        if(.not.CONTROLV)then
        write(kanalp,*)'vdwForceR1:'
        do i=1,natom
        write(kanalp,*)'ATOM  ',i,(vdwForceR1(3*i-3+k),k=1,3)
        end do !i
        end if
c
        write(kanalp,*)'coulForceR1: CONTROL:'
        do i=1,3*natom
        a=coulForceR1(i)
        CONTROLV = .false.
        if(a .ge. 0.0 .or. a .lt. 0.0)CONTROLV = .true.
        end do!i
        if(.not.CONTROLV)then
        write(kanalp,*)'coulForceR1:'
        do i=1,natom
        write(kanalp,*)'ATOM  ',i,(coulForceR1(3*i-3+k),k=1,3)
        end do!i
        end if
c
        write(kanalp,*)'coulForceR2: CONTROL:'
        do i=1,3*natom
        a=coulForceR2(i)
        CONTROLV = .false.
        if(a .ge. 0.0 .or. a .lt. 0.0)CONTROLV = .true.
        end do!i
        if(.not.CONTROLV)then 
        write(kanalp,*)'coulForceR2:'
        do i=1,natom
        write(kanalp,*)'ATOM  ',i,(coulForceR2(3*i-3+k),k=1,3)
        end do !i
        end if
c
        write(kanalp,*)'restr1AtForce: CONTROL:'
        do i=1,3*natom
        a=restr1AtForce(i)
        CONTROLV = .false.
        if(a .ge. 0.0 .or. a .lt. 0.0)CONTROLV = .true.
        end do!i
        if(.not.CONTROLV)then 
        write(kanalp,*)'restr1AtForce:'
        do i=1,natom
        write(kanalp,*)'ATOM  ',i,(restr1AtForce(3*i-3+k),k=1,3)
        end do !i
        end if
c        
         write(kanalp,*)'atomSolFr: CONTROL:'
        do i=1,3*natom
        a=atomSolFr(i)
        CONTROLV = .false.
        if(a .ge. 0.0 .or. a .lt. 0.0)CONTROLV = .true.
        end do!i
        if(.not.CONTROLV)then
        write(kanalp,*)'atomSolFr:'
        do i=1,natom
        write(kanalp,*)'ATOM  ',i,(atomSolFr(3*i-3+k),k=1,3)
        end do !i
        end if
c
        if(.not.CONTROLV)then
        write(kanalp,*)'atomVelp: atomVelm:'
        do i=1,natom
        write(kanalp,*)'ATOM  ',i,
     &  (atomVelp(3*i-3+k),k=1,3),(atomVelm(3*i-3+k),k=1,3)
        end do !i
        end if
c
        return
        end
c
        subroutine printHB128EngForces
c
        include 'xyzPDBsize.h'
        include 'xyzPDBinfo.h'
        include 'enForce.h'
        include 'hbond128.h'
        include 'mdAtomXYZvel.h'
        include 'optionPar.h'
c
        include "output.h"
        include "kanalUse.h"
        include "statusMessg_mDyn.h"
        include "filedat.h"
c
        integer i,k,i4
        integer ix,ih,iy
        integer ix3,ih3,iy3
        real a,s1,s2,s3
        real xh(3),yh(3),xy(3)
c
        real hbAngPar(6)
        real eHbAng
        real x3(3),h2(3),y1(3)
        real f1(3),f2(3),f3(3)
        real dxh,dyh,dxy
        real ch
        integer kanalp
        logical CONTROLHB128
        character*(charLenMAX) hb128EngFile
        logical OPT_wrhb128File 
        integer kanalHb128
c
        kanalp = kanalRunOut
c
        CONTROLHB128 = .true.
        OPT_wrhb128File = .true.
        hb128EngFile = './mdHb128Efin.res'
        kanalHb128 = kanalpHB128
c
        if(nMolNameLett .ge. 1)
     &  hb128EngFile = molNameArgLine(1:nMolNameLett)//
     &                 ".mdHb128Efin.res"
c
        if(OPT_resFileDir)then
        hb128EngFile=resultFileDir(1:nResFileDirLett)//'mdHb128Efin.res'
        if(nMolNameLett .ge. 1)then
        hb128EngFile=resultFileDir(1:nResFileDirLett)//
     &  molNameArgLine(1:nMolNameLett)//'.mdHb128Efin.res'
        end if
        end if
c
        if(OPT_wrhb128File)then
        open(unit=kanalHb128, file=hb128EngFile,form='formatted',
     &       status='unknown')
c
c         rewind kanalHb128res
        end if !OPT_wrhb128File
c
        if(CONTROLHB128)then
        a = 0.0
        write(kanalp,*)'hB128 energy List:'
        write(kanalp,'(a58)')
     &  ' ihB:   Y    Hx   X   eHB128   dyh     dxy    aXhY     dxh'
        if(OPT_wrhb128File)then
        write(kanalHb128,*)'hB128 energy List:'
        write(kanalHb128,'(a58)')
     &  ' ihB:   Y    Hx   X   eHB128   dyh     dxy    aXhY     dxh'
        end if !OPT_wrhb128File
c
        do i = 1,nHBbondHxY
        iy=hB128List(3*i-2)
        ih=hB128List(3*i-1)
        ix=hB128List(3*i)
        iy3=iy*3-3
        ih3=ih*3-3
        ix3=ix*3-3
c..
        dxy=0.0
        dyh=0.0
        dxh=0.0
        do k =1,3
        x3(k)=atomXYZ0(ix3+k)
        h2(k)=atomXYZ0(ih3+k)
        y1(k)=atomXYZ0(iy3+k)
        xh(k)=atomXYZ0(ix3+k)-atomXYZ0(ih3+k)
        yh(k)=atomXYZ0(iy3+k)-atomXYZ0(ih3+k)
        xy(k)=atomXYZ0(ix3+k)-atomXYZ0(iy3+k)
        dxy=dxy + xy(k)**2
        dyh=dyh + yh(k)**2
        dxh=dxh + xh(k)**2
        end do !k
c
        dxy = sqrt(dxy)
        dyh = sqrt(dyh)
        dxh = sqrt(dxh)
c 
        call vectorNrm1(xh)
        call vectorNrm1(yh)
        call scalarp(xh,yh,ch)
c
       ch = acos(ch)*(180.0/3.1415927)
cx       ch = 180.0 - ch
c
c  X...Y dist and X-H..Y cosA
c        
        write(kanalp,'(4(i4,1x),5(f7.2,1x))')
     &  i,(hB128List(3*i-3+k),k=1,3), hBHxYeng128List(i), 
     &  dyh,dxy,ch,dxh
c
        if(OPT_wrhb128File)then 
            write(kanalHb128,'(4(i4,1x),5(f7.2,1x))')
     &  i,(hB128List(3*i-3+k),k=1,3), hBHxYeng128List(i),
     &  dyh,dxy,ch,dxh
        end if ! OPT_wr
c
         a = a + hBHxYeng128List(i)
c
c check the hB128 energy/forces
c assign hB128 Param
        i4=i*4-4
        do k=1,4
        hbAngPar(k)=hB128PotParList(i4+k)
        end do !k
c
        hbAngPar(5) = hB128AgLpar(1)
        hbAngPar(6) = hB128AgLpar(2)
c
        call getHb128angEforce
     &       (y1,h2,x3,
     &                  hbAngPar,eHbAng,f1,f2,f3)
c
        end do !i hBonds
c
        write(kanalp,*)'hB128EngTot:    ',a
        if(OPT_wrhb128File)then 
        write(kanalHb128,*)'hB128EngTot:    ',a
c
        close (kanalHb128)
        end if !OPT_wr 
c write allHb128 forces
         write(kanalp,*)' hb128 forces:'
         s1 = 0.0
         s2 = 0.0
         s3 = 0.0
         write(kanalp,*)' atom     fx       fy      fz',
     &   '     fsx   fsy    fsz'
         do i = 1,natom
         s1 = s1 + hBHxY128force(3*i-2)
         s2 = s2 + hBHxY128force(3*i-1)
         s3 = s3 + hBHxY128force(3*i)
         write(kanalp,'(a6,i5,1x,6f8.2)')
     &   'ATOM  ',i,(hBHxY128force(3*i-3+k),k=1,3),s1,s2,s3
         end do !i
c
        end if !CONTROLHB128
c
        return
        end 
