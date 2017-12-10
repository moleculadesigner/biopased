c init_ExWatShellRestrPot
c calculate restrain Potential for ExpWaterShell molecules      
c for water molecules
c YN Vorobjev   2009
c
	subroutine initExWatShRestrSASPot
     &  	(ic_exWatSASPot,natomSOLV,atomXYZsolv,
     &                          	eng_exWatSASP,ff_exWatSASP)
c
c INPut: ic_exWatSASPot - counter of the routine call
c        natomSOLV,atomXYZsolv(*) - LIG XYZ
c OUT: eng_exWatSASP - enrergy of expuLsionPotential
c      ff_exWatSASP(iw) - forces on WatSAS SHELL atoms, iw=1,..,3*nWatSheLLMolec
c       
	include 'xyzPDBsize.h'
c
        include "output.h"
        include "kanalUse.h"
        include "statusMessg_mDyn.h"
c
        include "xyzPDBinfo.h"
        include "solvate01.h"
        include 'dataSASdr.h'
        include 'nbPairSAS03.h'
	include "exWatShRestrSASPot.h"
	include "enForce_PW.h"
c
        integer ic_exWatSASPot
	integer natomSOLV
	real atomXYZsolv(*)
        real eng_exWatSASP
	real ff_exWatSASP(*)
c local
        real ff_exWatLoc(3*nAtomWatMAX)
c
        logical CONTROL,CONTROL1
	integer kanalp
        integer nnbpLMX 
	integer nnbpLT 
        real rcut
	integer is,is3,ia,ia3,k,i
	real KexpL, ddPentMX
c
        kanalp = kanalRunOut
	CONTROL = . false. !.true.
	CONTROL1 = .false. !.true. 
c 
        ic_exWatSASPot=ic_exWatSASPot+1
	areaProtSAS04=0.0
	natomExWatSh = natomWSolv
c
        if(CONTROL1)then
        write(*,*)'initExWatShRestrSASPot Start!ic_exWatSASPot:',
     &	ic_exWatSASPot
        end if
c define expLSASPot.h data structure from dataSASdr.h
        if(CONTROL)
     &  write(kanalp,*)
     &	'initExWatShRestrSASPot: init expLSASPot.h data structure'
c
        nsurfSASP01=nsurfAt(1) 
        do is=1,nsurfSASP01
	ia=nsurfAtList(is,1)     
	atSurfArP01(is)=atSurfAr(ia,1)
        areaProtSAS04=areaProtSAS04 + atSurfArP01(is) 
c
	ia3=3*ia-3
	is3=3*is-3
	do k=1,3
        atSurfNrmP01(is3+k)=atSurfNrm(ia3+k,1) 
	atSurfXYZP01(is3+k)=atSurfXYZ(ia3+k,1)
	end do!k
	end do!is
c define initial dSOLVSHELLVol
        if(ic_exWatSASPot.eq.1)
     &        dSOLVSHELLVolFix=areaProtSAS04*dSOLVSheLL
c redefine dSOLVSheLL
        dSOLVSheLLvar=dSOLVSHELLVolFix/areaProtSAS04
c
c use nbPairSAS03.h data structure for nbpairList,startnbPairL,nnbPairL,nbpairListD2
c
        if(CONTROL1)then
        write(kanalp,*)
     &	'initExWatShRestrSASPot: w2: start expLPotLgSASneigList!'
        write(*,*)
     &	'initExWatShRestrSASPot: w2: start expLPotLgSASneigList!'
	end if !C
c
        nnbpLMX = nnbpLSAS03MAX
	nnbpLT = 0
	rcut_exWSASPot = 12.0 
	rcut = rcut_exWSASPot
c
c neighbor List of proteinSAS atoms from Ligand atoms: exWatShRestrSASPot.h
	call exWatShSASpAtNeigList(ic_exWatSASPot,
     &    nsurfSASP01,atSurfXYZP01,rcut,
     &    natomSOLV,atomXYZsolv,
     &    nbpairListSAS03,startnbPairLSAS03,nnbPairLSAS03,
     &    nbpairListSAS03D2,indxSASatNeighExWatSh,
     &    d2SASatNeighExWatSh,nnbpLT,nnbpLMX)
c
        if(CONTROL1)then
        write(kanalp,*)
     &  'initExWatShRestrSASPot: w3: finish exWatShSASpAtNeigList!'
        write(*,*)
     &  'initExWatShRestrSASPot: w3: finish exWatShSASpAtNeigList!'
        end if !C
c expulsion potential parameters: set:
        KhSASexWatSh = 3.0   ! 1.2  !kcal/mol/A^2
	ddSASexWatShMX = 2.5 ! 1.5 A
c calculate eng_exWatSASP and forces ff_exWatSASP() for WatSolvSHeLL :
        call get_expLSASPot02(natomSOLV,atomXYZsolv,
     &                          eng_exWatSASP,ff_exWatSASP)
c
       if(CONTROL1)then
       write(kanalp,*)'initExWatShRestrSASPot: ff_exWatSASP:'
       do i=1,natomSolutMol+natomWSolv
       write(kanalp,*)"ATOM  ", 
     &    i,atomName(i),resName(i),(ff_exWatSASP(3*i-3+k),k=1,3) 
       end do!i
       write(*,*)'initExWatShRestrSASPot FINish!'
       STOP
	end if !C
c
	return
	end
