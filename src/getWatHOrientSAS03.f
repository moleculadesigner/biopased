c calculate waterBridge Hydrogen orientation
c for solvateMoL03  solvation model ElectrFieldWaterShell
c
	subroutine getWatHOrientSAS03(atomXYZs,isd,sdXYZ,wDip,
     &    wMolXYZ,wMolAtQ,wMolAtName,wMolResName,wMolEng)
c
        include 'xyzPDBsize.h'
c
        include "output.h"
        include "kanalUse.h"
        include "statusMessg_mDyn.h"
c
        include 'xyzPDBinfo.h'
        include 'dataSASdr.h'
        include 'nbPairSAS03.h'
        include 'solvWBrg01.h'
c
        real atomXYZs(*)
        integer isd
        real sdXYZ(*)
        real wDip(*)           !unit wDip vector
        real wMolXYZ(*)
        real wMolAtQ(*)
        real wMolEng
        character*4 wMolAtName(*)
        character*4 wMolResName(*)
c local
        real rcut(3)
        real Trot(9)   
        real csfi,snfi,qe
        integer ngridMAX
        parameter(ngridMAX=17)
        integer ngridMX,ishOr
        integer nOrMX,ior,iorMax
        integer ia,ia3,ja,ja3,ip,k,jn,id3,i2
        integer ig,ig3,ig1,ig2,ig13,ig23
        real v10(3),v1(3),v20(3),v2(3),v30(3),v3(3)
        real s2,s3,zscale
        real grid0(3*ngridMAX)
        real grid1(3*ngridMAX)
        real gridePot(ngridMAX),grideField(3*ngridMAX)
        real dgridePot,dgrideField(3)
        real gridePotMiN,gridePotS2,zWH12
        real cReps
        real qHWbgMol
        integer epotVar_OPT
        integer igridePotMin
c
	data grid0/0.0 , 0.0,  0.0,
     &             0.757, 0.0,  0.00,-0.757,  0.0,  0.00,
     &             0.699,0.289, 0.00,-0.699,-0.289, 0.00,
     &             0.534,0.534,0.00,-0.534,-0.534,0.00,
     &             0.289,0.699,0.00,-0.289,-0.699,0.00,
     &             0.0,  0.757, 0.00, 0.0,  -0.757, 0.00,
     &             -0.289,0.699,0.00,0.289,-0.699,0.00,
     &             -0.534, 0.534,0.00,0.534,-0.534,0.00,
     &             -0.699,0.289, 0.00,0.699,-0.289, 0.00/  
c
        logical CONTROL,CONTROL0
        integer kanalp
c
        ngridMX = ngridMAX
        nOrMX = ngridMAX/2
        ishOr = ngridMAX - 2*nOrMX
        gridePotMiN = 1.0e10
        igridePotMiN = 1
c 
        CONTROL=.false. 
        CONTROL0=.false. 
        kanalp = kanalRunOut
c
        qHWbgMol = 0.417
cx        epotVar_OPT=1
        epotVar_OPT=2  ! cReps
        cReps = 3.3
        zWH12 = 0.577 !  dipMWEFS_OPT/0.834
c define loc xyz system
        do k=1,3
        v10(k)=0.0
        v20(k)=0.0
        v30(k)=0.0
        end do
        v10(1)=1.0
        v20(2)=1.0
        v30(3)=1.0
        do k=1,3
        v3(k)=wDip(k)
        end do!k
        call vectorNrm1(v3)
        v2(1)=0.0
        v2(2)=0.0
        v2(3)=1.0
        call scalarp(v3,v2,s2)
        if(s2 .ne. 1.0)then
        do k=1,3
        v2(k)=v2(k)-v3(k)*s2
        end do
        call vectorNrm1(v2)
        else
        v2(1)=0.0
        v2(2)=1.0
        v2(3)=0.0
        end if
        call vectorp(v2,v3,v1)
c
        call getRotMatrix02(v10,v20,v30,v1,v2,v3,Trot)
c 
        rcut(1) = RcutEFSMAX
        rcut(2) = rcut(1)**2
        rcut(3) = rcut(2)*rcut(1)
c
        ip=1
        ia = dot_IATNUM(isd,ip)    ! SAS atom for iProbe=1
        ia3=3*ia-3
        id3=3*isd-3
c
        do ig=1,ngridMAX
        ig3=3*ig-3
c
        call vectMatrx3Prod(grid0(ig3+1),Trot,grid1(ig3+1))
c
        if(CONTROL)then
        write(kanalp,*)'Rotation:'
        write(kanalp,*)'wDip:            ',(wDip(k),k=1,3)
        write(kanalp,*)'ig: grid0(ig3+k):',(grid0(ig3+k),k=1,3)
        write(kanalp,*)'ig: grid1(ig3+k):',(grid1(ig3+k),k=1,3)
        end if
        zscale = 0.5
        if(ig .eq. 1)zscale = -0.5
        do k=1,3
        grid1(ig3+k)=grid1(ig3+k)+sdXYZ(k)+zscale*wDip(k)*zWH12  !in global XYZ
        end do!k
c
        if(CONTROL)then
        write(kanalp,*)'Final=Rotation:+Transl:'
        write(kanalp,*)'ig: grid1(ig3+k):',(grid1(ig3+k),k=1,3)
        end if
c
c calculate gridePot in grid1(ig3+k) point from atom ia
c
        qe=eLectScaleEFSMod*atomQ(ia)
        if(epotVar_OPT .eq. 1)
     &  call ePotEfieldSAS(atomXYZs(ia3+1),qe,grid1(ig3+1),rcut,
     &                   gridePot(ig),grideField(ig3+1))
        if(epotVar_OPT .eq. 2)
     &   call ePotEfieldSAS03(atomXYZs(ia3+1),qe,grid1(ig3+1),rcut,
     &                   cReps,gridePot(ig),grideField(ig3+1))
c
c field from neighbours of ia:
c
        if(nnbPairLSAS03(ia) .ge. 1) then
        do jn = startnbPairLSAS03(ia),
     &          startnbPairLSAS03(ia) + nnbPairLSAS03(ia)-1
           ja = nbpairListSAS03(jn)   ! atom j
           ja3 = ja*3-3
           qe=eLectScaleEFSMod*atomQ(ja)
c
           if(atomQ(ja) .ne. 0.0 )then
           qe = eLectScaleEFSMod*atomQ(ja)
c
           if(epotVar_OPT .eq. 1)
     &     call ePotEfieldSAS(atomXYZs(ja3+1),qe,grid1(ig3+1),rcut,
     &                   dgridePot,dgrideField)
          if(epotVar_OPT .eq. 2)
     &    call ePotEfieldSAS03(atomXYZs(ja3+1),qe,grid1(ig3+1),rcut,
     &                  cReps,dgridePot,dgrideField)
c
           gridePot(ig)=gridePot(ig)+dgridePot
c
           do k=1,3
           grideField(ig3+k)=grideField(ig3+k)+dgrideField(k)
           end do !k
c
           end if !atomQ(ja) .ne. 0.0
c
           end do !jn
           end if !nnbPairLSAS03(ia) .ge. 1
c
        end do !ig
c find Optimal Orient due to MIN Electrostatic energy gridePot()
        do ior=1,nOrMX 
        i2=ior*2-1+ishOr
        gridePotS2 = gridePot(i2)+gridePot(i2+1) - 2.0*gridePot(1)
        gridePotS2 =  gridePotS2*qHWbgMol ! electrEng
        if(CONTROL)then
        write(kanalp,*)
     &  'getWatHOrientSAS03:ior,gridePotS2:',ior,gridePotS2,
     &  ' gridePotO:',gridePot(1)	
       
        end if !C
        if(gridePotS2 .lt. gridePotMiN)then
        gridePotMiN = gridePotS2 
        iorMaX=ior
        end if
        end do!ior
c
c
        if(CONTROL)then
        write(kanalp,*)'getWatHOrientSAS03:iorMiN:',iorMaX,
     &  ' gridePotS2(eCol):',gridePotMin
        end if        
c
        wMolEng = gridePotMin
        ig1=2*iorMaX-1 +ishOr
        ig2=ig1+1
        ig13=3*ig1-3
        ig23=3*ig2-3
        do k=1,3
cx        wMolXYZ(k)=sdXYZ(k)
        wMolXYZ(k)=grid1(k)
        wMolXYZ(k+3)=grid1(ig13+k)
        wMolXYZ(k+6)=grid1(ig23+k)
        end do!
        wMolAtQ(1)=-0.834
        wMolAtQ(2)= 0.417
        wMolAtQ(3)= 0.417
        wMolAtName(1) = 'Owb'
        wMolAtName(2) = 'H1wb'
        wMolAtName(3) = 'H2wb'
        wMolResName(1)='HOHb'   
c
        if(CONTROL0)then
        write(kanalp,*)'getWatHOrientSAS03:wMolXYZ:'
               write(kanalp,7002)
     &   'ATOMwb',isd       ,'O   ',ia,                        
     &   (wMolXYZ(k),k=1,3)
         write(kanalp,7002)
     &   'ATOMwb',isd       ,'H1  ',ia,
     &   (wMolXYZ(k+3),k=1,3)
         write(kanalp,7002)
     &   'ATOMwb',isd       ,'H2  ',ia,
     &   (wMolXYZ(k+6),k=1,3)
        end if        
c        
	return
7002   format(a4,1x,i6,2x,a4,5x,i4,4x,3f8.3,1x,f6.2,3f6.2,f5.2) ! PDB 
        end
