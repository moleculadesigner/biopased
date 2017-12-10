c calculw BornRad for molecule atom
c
	subroutine initBornRadSAS(atomXYZ)
c
cx        implicit none
        include 'xyzPDBsize.h'
        real atomXYZ(*)
        include 'xyzPDBinfo.h'
        include 'dataSASdr.h'
        include 'solvGBorn.h'
        include "output.h"  
c
        real  probeD_bornRad_OPT
        integer ip
        integer ndotMX
        integer kanalp
        logical CONTROL
c
        kanalp = kanalRunOut
        CONTROL = .false.
c        
        ndotMX = ndotMAX
        probeD_bornRad_OPT = solMolRad_OPT ! 2.8
c atomic radii
        call getAtomRadSAS03(natom,atomName,atomRad)
c
        ip = 1
c
cx        dotden = 4.0   ! daet tochnost' ~ 0.1 A, rascheta BornRad
c                        ! in inputMDpar.f 
        dprobe(ip) = probeD_bornRad_OPT  ! 
c
        call  surf_SAS04(atomXYZ,atomRad,
     &                 natom,dprobe(ip),dotden,
     &                 dotXYZ(1,ip),dotnrm(1,ip),dotarea(1,ip),
     &                 dot_IATNUM(1,ip),ndot(ip),ndotMX,
     &                 nsurfAt(ip),nsurfAtList(1,ip),
     &                 atSurfAr(1,ip),atSurfNrm(1,ip),
     &                 atSurfXYZ(1,ip),atSurfProj(1,ip),
     &                 bindSiteAt01(1,ip),
     &                 dot_jATNUM(1,ip),
     &                 head_dotNum,linkListDotAt,
     &                 atDistToSAS,atSASexp,
     &                 dot_eField,dot_ePot)
c
        nSAScall = nSAScall + 1 
c
        write(kanalp,*)'initBornRadSAS : SAS call nSAScall=',nSAScall
c
c
        call getBornRadSAS(atomXYZ,atomRad,
     &              natom,
     &              dotXYZ(1,ip),dotnrm(1,ip),dotarea(1,ip),
     &              ndot(ip),
     &              nsurfAt(ip),nsurfAtList(1,ip),
     &              atSurfAr(1,ip),atSurfNrm(1,ip),atSurfXYZ(1,ip),
     &              head_dotNum,linkListDotAt,
     &              atBornRad,atBornRadDr)
c
         return
         end
