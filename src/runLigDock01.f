c run ligand docking method 01:
c     finds optimal ligand position
c     in the vicinity of the point: LigGcPosXYZin
c
c     mdSA start from uniform over 4PI orientation 
c     nOrient = nXZ*nZ
c     nXZ : vertor direction on the unut sphere
c     nZ  : rotation around vector direction
c
c fineGraneRotGrid: nXZ*nZ= nV2 = 18*8 = 144
c coarseGraneRotGrid:       nV1 = 6 *4 = 24
c
c Y.N. Vorobjev,  2003
c
	subroutine runLigDock01
     &     (ligGcPosXYZin,nLigDockPosFin,engLigDockPosFin,
     &                                     atLigDockPosXYZ)
c
c In: ligGcPosXYZin(3) - initial ligGeoCentrPos
c 
c OUT: nLigDockPosFin,
c      engLigDockPosFin
c      atLigDockPosXYZ
c
c        implicit none
        include 'xyzPDBsize.h'
c
        include "output.h"
        include "kanalUse.h"
        include "statusMessg_mDyn.h"
c
        include 'xyzPDBinfo.h'
        include 'xyzPDBcrd.h'
        include 'engXYZtra.h'
        include 'enForce.h'
        include 'mdRunPar.h'
        include 'optionPar.h'
        include 'simAnneal.h'
c
        real ligGcPosXYZin(*)
        real atLigDockPosXYZ(*)
        real engLigDockPosFin(*)
        integer nLigDockPosFin
c
        include 'ligInfo.h'
c
c data to construct a set of rotation matrix
c
        integer nsZ2MAX
        integer nsXZMAX
cv2        parameter (nsZ2MAX=8)
cv2        parameter (nsXZMAX=18)
        parameter (nsZ2MAX=4)
        parameter (nsXZMAX=6)
        real angXZ(2*nsXZMAX)
        real angZ2(nsZ2MAX)
c local
        character*(charLenMAX) fileLigEng,fileLigEngPdb
        integer iz2,ixz
        integer i,ia,i2,i3,j,ia3
        integer nOrient
        real PI
        logical CONTROL
        logical OPT_doMDloc
        logical OPT_MDSAloc
        integer kanalLigEng
        integer kanalp
c
        real tz1(3,3),tz2(3,3),tx(3,3)
        real tet(3)     ! tet(1)=tetx, tet(3) = tetz, tet(2)=tety
c rotation angles in PI rad unit
c   angXZ = (X,z1)
cv2        data angXZ /0.0,0.0, 0.0,0.25, 0.5,0.25, 1.0,0.25, 1.5,0.25, 
cv2     &              0.0,0.5, 0.25,0.5, 0.5,0.5, 0.75,0.5,  
cv2     &              1.0,0.50, 1.25,0.5, 1.5,0.50, 1.75,0.5, 
cv2     &              0.0,0.75, 0.5,0.75, 1.0,0.75, 1.5,0.75, 0.0,1.0/ 
cv1: 
         data angXZ /0.0,0.0,
     &               0.0,0.5, 0.5,0.5, 1.0,0.5, 1.5,0.5,
     &               0.0,1.0/
c 
c angZ2 = z2
cv2        data angZ2 /0.0,0.25,0.50,0.75,1.0,1.25,1.50,1.75/
cv1:
         data angZ2 /0.0,0.50,1.0,1.50/ 
c  T = Tz2*Tx*Tz1,   xout = T*xin
c
        kanalp = kanalRunOut
        CONTROL=.true.
        OPT_doMDloc = .false.
        OPT_MDSAloc = .true.
c
        PI = 3.1415927
c
        fileLigEngPdb = 'mdSALigOr.pdb'
        fileLigEng = 'mdSALigOrEng.res'
c
         kanalLigEng = 33
         open(unit=kanalLigEng, file=fileLigEng, form='formatted',
     &       status='unknown')
c
c make and run Lig orientations
        nOrient = 0
        do ixz = 1,nsXZMAX
c make initial orientation atomLigXYZ0()
        i2 = 2*ixz-1
        tet(1)=angXZ(i2+1)*PI
        tet(2)=angXZ(i2)*PI
        do iz2 = 1,nsZ2MAX
        tet(3) = angZ2(iz2)*PI
c
c rotate LigXYZ = atomLigXYZ0c (*) = atLigXYZ relative ligCMass
c
        nOrient = nOrient + 1
        call rotateMolZXZ(nAtomInLig,atomLigXYZ0c,tet,
     &                        ligGcPosXYZin,atomLigXYZ) 
c               
c rotated Lig start ligandXYZ = atomLigXYZ(*) [ligInfo.h]
c
        do i=1,nAtomInLig
        i3 = 3*i-3
        ia = atomInLigList(i)
        ia3 = ia*3-3
c
        atomXYZ(ia3+1)=atomLigXYZ(i3+1)
        atomXYZ(ia3+2)=atomLigXYZ(i3+2)
        atomXYZ(ia3+3)=atomLigXYZ(i3+3)
c
        end do ! i
c
        if(CONTROL)then
        write(kanalp,*)'runLigDock01: ligStartXYZ: nOrient:',nOrient
        write(kanalp,*)'              ixz,iz2:',ixz,iz2,
     &                 ' tet(1,2,3):',tet
        do i=1,nAtomInLig
        i3 = 3*i-3
        ia = atomInLigList(i)
        write(kanalp,7071)"ATOM  ",
     &     atomNumb(ia),atomName(ia),resName(ia),chName(ia),resNumb(ia),
     &     (atomLigXYZ(i3+j),j=1,3), atomQ(ia)
        end do
        end if !CONTROL
c
c do md + mdSA run, keep final LigXYZ and engPot
c       
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
        if(OPT_doMDloc)then
        call initMDStart(tempT0)
c 
        optra = 1     !open file : './job/xyzMd.tra'
        wtra = 1
        cltra = 0
        atype = 1         ! 0- free MD, 1- NTV
c
        call mdRun(ntimeMX,ntime0,ntime,ntimeR1,ntimeR2,
     &             ntimeF1,ntimeF2,ntimeF3,deltat,
     &             tempTg,tempTg,tauTRF,atype,optra,wtra,nwtra,cltra)
c
        end if ! OPT_doMDloc
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
        if(OPT_MDSAloc) then
c
        if(nSAstep .gt. nSAstepMAX)then
        write(kanalp,*)'ERROR!!: nSAstep .gt. nSAstepMAX!!'
        stop
        end if!nSAstep
c
        call simAnnealing(nSAstep,nSAparMX,SAProtcol)
c
        end if !MDSAloc
c
c result for the orientation:
c
        fileLigEngPdb = 'mdSALigOr.pdb'
        call wrLigEngXyz(fileLigEngPdb,nOrient)
c
        write(kanalLigEng,'(a22,i4)')'$ENELIG: nOrient:',nOrient
        write(kanalLigEng,'(a12,f14.5)')'eVbondDefLG:',eVbondDefLg
        write(kanalLigEng,'(a12,f14.5)')'eVangDefLG :',eVangDefLg
        write(kanalLigEng,'(a12,f14.5)')'eImpDefLG  :',eImpDefLg
        write(kanalLigEng,'(a12,f14.5)')'eTorsDefLG :',eTorsDefLg
        write(kanalLigEng,'(a12,f14.5)')'engVDWR1LG :',engVDWR1Lg
        write(kanalLigEng,'(a12,f14.5)')'engCOULR1LG:',engCOULR1Lg
        write(kanalLigEng,'(a12,f14.5)')'engCOULR2LG:',engCOULR2Lg
        write(kanalLigEng,'(a12,f14.5)')'restr1EngLG:',restr1EngLg
        write(kanalLigEng,'(a12,f14.5)')'eRstHW1MLLG:',restr1MHWEngLg
        write(kanalLigEng,'(a12,f14.5)')'eGeoDefLG  :',eGeoDefLg
        write(kanalLigEng,'(a12,f14.5)')'engCOULLG  :',engCOULLg
        write(kanalLigEng,'(a12,f14.5)')'engSolvLG  :',molSolEnLg
        write(kanalLigEng,'(a12,f14.5)')'engPOTENTLG:',engPOTENTLg
        write(kanalLigEng,'(a7)')'ENDLIG '
c        
        end do !iz2
        end do !ixz
c
	return
7071    format(a6,1x,i4,2x,a4,a4,a1,i4,4x,3f8.3,7x,f8.5) ! PDB rq
        end

