c calculate forces from hydrophobic cavity
c solvateMoL03  
c
	subroutine getHPolvForceSAS03(atomXYZs,
     &                   hpSoLEnergy,hpSoLForceAt)
c
c calculates watBrgEnergySoL,watBrgForceAt(*) : enForce.h
c
        include 'xyzPDBsize.h'
c
        include "output.h"
        include "kanalUse.h"
        include "statusMessg_mDyn.h"
c
        include 'xyzPDBinfo.h'
        include 'dataSASdr.h'
        include 'solvWBrg01.h'
	include 'solvate01.h'
	include 'optionPar.h'
c
        real atomXYZs(*)
        real hpSoLForceAt(*)
        real hpSoLEnergy
c local
        integer natomSAS,ndotMX
        real atSurfArTot
        real hpSoLForceSUM(3)
        real hpSolForceSumP(3)
        real hpSolForceSumM(3)
        real scalP(3),scalM(3)
        real gamaHP,pi,snrm
        integer k,ia,ia3,ip,is
c
        integer iTypeCorr
        logical CONTROL,CONTROL1
        integer kanalp
c
	CONTROL = .false.  ! control print & stop execution
        CONTROL1 = .false.
        kanalp = kanalRunOut
c
cx        dotden = 4.0   ! daet tochnost' ~ 0.1 A, rascheta BornRad
c                        ! in inputMDpar.f
        ip = 1
        dprobe(ip) = solMolRad_OPT  !
c
        ndotMX = ndotMAX
        natomSAS = natom
        call getAtomRadSAS03(natomSAS,atomName,atomRad)
c calculate SAS for soluteMolecule
        if(OPT_SolvateExWat)natomSAS = natomSMOL
c
        call  surf_SAS04(atomXYZs,atomRad,
     &                 natomSAS,dprobe(ip),dotden,
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
        gamaHP = 0.038  ! reproduce Wat-Wat PMF in water! 0.070=expWater!cal/A^2
	if(OPT_SolvateExWat)gamaHP = 0.069      !vdwP-W is explicit term
c
        pi=3.1415927
        iTypeCorr = 1    ! 0=equal correction, 1=proportiona correction
c
        do k=1,3*natomMAX
        hpSoLForceAt(k) = 0.0
        end do
        hpSoLEnergy = 0.0 
        do k=1,3
        hpSoLForceSUM(k) = 0.0
        hpSoLForceSumP(k) = 0.0
        hpSoLForceSumM(k) = 0.0 
        end do!k
c
        if(CONTROL)then
        write(kanalp,*)'getHPSAS03Force: start!'
        end if !
c
        ip = 1
c
        atSurfArTot = 0.0
c
        do is=1,nsurfAt(ip)
        ia = nsurfAtList(is,ip)
        ia3=3*ia-3
        snrm = sqrt(atSurfNrm(ia3+1,ip)**2+atSurfNrm(ia3+2,ip)**2
     &         +atSurfNrm(ia3+3,ip)**2)
c
        atSurfArTot=atSurfArTot + atSurfAr(ia,ip)
        hpSoLEnergy = hpSoLEnergy + atSurfAr(ia,ip)*gamaHP
c
        if(atSurfAr(ia,ip) .gt. 0.0)then
        do k=1,3
        hpSoLForceAt(ia3+k)= -2.0*gamaHP*sqrt(pi*snrm*atSurfAr(ia,ip))*
     &  atSurfNrm(ia3+k,ip)/snrm 
c
        hpSoLForceSUM(k) = hpSoLForceSUM(k) + hpSoLForceAt(ia3+k)
        if( hpSoLForceAt(ia3+k) .gt. 0.0 )
     &     hpSoLForceSumP(k) = hpSoLForceSumP(k) + hpSoLForceAt(ia3+k)
         if( hpSoLForceAt(ia3+k) .lt. 0.0 )
     &     hpSoLForceSumM(k) = hpSoLForceSumM(k) + hpSoLForceAt(ia3+k)
        end do !k
        end if !atSurfAr(ia) .gt. 0.0
c
        if(CONTROL)then
         write(kanalp,7003)
     &  'ATOM  ',ia,atomName(ia),resName(ia),resNumb(ia),
     &  (atomXYZs(ia3+k),k=1,3),(hpSoLForceAt(ia3+k),k=1,3),
     &  atSurfAr(ia,ip)
        end if !C
        end do !is
c
        if(CONTROL)then
        write(kanalp,*)'hpSoLEnergy:',hpSoLEnergy
        write(kanalp,7003)
     &  'SUMhpF',ia,atomName(ia),resName(ia),resNumb(ia),
     &  (hpSoLForceSUM(k),k=1,3),
     &   atSurfArTot
         end if
c
c correct forces by condition of zero sum F
c
        do k=1,3
        scalP(k) = 0.0
        scalM(k) = 0.0
        if(hpSoLForceSumP(k) .gt. 0.0)
     &        scalP(k)=hpSoLForceSUM(k)*0.5/hpSoLForceSumP(k)
        if(hpSoLForceSumM(k) .lt. 0.0)
     &        scalM(k)=hpSoLForceSUM(k)*0.5/hpSoLForceSumM(k)
        hpSoLForceSUM(k) = hpSoLForceSUM(k)/nsurfAt(ip)
        end do !k
c
        do is=1,nsurfAt(ip)
        ia = nsurfAtList(is,ip)
        ia3=3*ia-3
        do k=1,3
        if(iTypeCorr .eq. 1)then
        if(hpSoLForceAt(ia3+k) .gt. 0.0)
     &  hpSoLForceAt(ia3+k) = hpSoLForceAt(ia3+k)*(1.0 - scalP(k))
        if(hpSoLForceAt(ia3+k) .lt. 0.0)
     &  hpSoLForceAt(ia3+k) = hpSoLForceAt(ia3+k)*(1.0 - scalM(k))
c
        else        
        hpSoLForceAt(ia3+k)=hpSoLForceAt(ia3+k) - hpSoLForceSUM(k) 
        end if !iTypeCorr
        end do!k
        end do !is
c
        if(CONTROL)then
        do k=1,3
        hpSoLForceSUM(k)=0.0
        do is=1,nsurfAt(ip)
        ia = nsurfAtList(is,ip)
        ia3=3*ia-3
        hpSoLForceSUM(k)=hpSoLForceSUM(k)+hpSoLForceAt(ia3+k)
        end do !is
        end do!k 
        write(kanalp,*)'scalP(k):',scalP
        write(kanalp,*)'scalM(k):',scalM 
        write(kanalp,*)'hpSoLEForceSUM: corrected:',
     &  (hpSoLForceSUM(k),k=1,3),
     &   ' atSurfTot:',atSurfArTot
        end if!C1
c
        if(CONTROL)STOP
c
	return
c
7003   format(a4,2x,i5,1x,a4,a4,1x,i4,4x,3f8.3,1x,3f8.3,1x,f8.3) ! PDB 
        end
