c initialization distant dependent dielectric model of Lazaridis 
c is NOT USED in bioPASED
c Yurii Vorobjev  2005
c
	subroutine initSASdielModLaZ(atomXYZ_s) 

        include 'xyzPDBsize.h'
        real atomXYZ_s(*)
        include 'xyzPDBinfo.h'
        include 'coulEnPar.h'
        include 'dataSASdr.h'
        include 'controlCall.h'
c
        include "output.h"
        include "kanalUse.h"
        include "statusMessg_mDyn.h"
c
        integer ip 
        integer ndotMX 
        integer kanalp
        logical CONTROL
c
       kanalp =  kanalRunOut
       CONTROL = .true.
c define param for COntrolPrint : controlCall.h
       ip = 1
       ncallSubTot(ip)=0
       nPrintNcall(ip) = 31 ! 99 
c
       if(coulVar_OPT .ne. 3)return
c define sasDielModelParameters Lazaridis:
c
c r05=p(3)*LN[p(4)-p(5)*atDistToSAS+p(6)*atSASexp] - tochka peregiba-inflection point of D(r)
c D(r) = p(1)+p(2)*(r/r05)**p(7)/[1. + (r/r05)**p(7)]
c
cx        epsSol_OPT = 80.0
cx        epsMol_OPT = 6.0 ! t4:2.0  ! t3:4.0  ! 12. md res for EPS inside proteins
cx        solMolRad_OPT = 2.8 ! diameter W mol
c
        parDrLZmodel(1) = epsMol_OPT
        parDrLZmodel(2) = epsSol_OPT - parDrLZmodel(1)
        parDrLZmodel(3) = 7.6
        parDrLZmodel(4) = 1.50 !1.25= ploho! t3:1.5
        parDrLZmodel(5) = 4.4   !  5.8
        parDrLZmodel(6) = 2.8   !  4.0
        parDrLZmodel(7) = 2.6    ! average
        parDrLZmodel(8) = 1.0   !  
c
c calculate SAS and SAS dielectric model parameters
        ndotMX = ndotMAX
        natomMX = natomMAX
c
        ip = 1
c
cx        dotden = 4.0  ! inputMDSApar.f
c
        dprobe(ip) = solMolRad_OPT  ! water
c
cx         call  surf_SAS02(atomXYZ_s,atomRad,
cx     &                 natom,dprobe(ip),dotden,
cx     &                 dotXYZ(1,ip),dotnrm(1,ip),dotarea(1,ip),
cx     &                 dot_IATNUM(1,ip),ndot(ip),ndotMX,
cx     &                 nsurfAt(ip),nsurfAtList(1,ip),atSurfAr(1,ip),
cx     &                 atSurfNrm(1,ip),bindSiteAt01(1,ip),
cx     &                 atDistToSAS,atSASexp)
c
               call  surf_SAS04(atomXYZ_s,atomRad,
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
         write(kanalp,*)
     & 'initSASdielModLaZ: doneInitialization dielSASLaZaridis model:'
c
         return
         end
