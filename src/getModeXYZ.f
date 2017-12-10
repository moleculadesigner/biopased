c
c Yury Vorobjev, 2005
c
c calculate XYZ coords of main Essential modes
c
	subroutine getModeXYZ(essModeTypeOPT,nMode)
c
c nMode - number of modes to be calculated
c
cx        implicit none
        include 'xyzPDBsize.h'
        include 'xyzPDBinfo.h' 
        include 'engXYZtra.h'
        include 'modeXYZtra.h'
        include 'filedat.h'
c
        include "output.h"
        include "kanalUse.h"
        include "statusMessg_mDyn.h"
c
        integer essModeTypeOPT,nMode
c local
        integer ipmMX_OPT
	integer i,im,iXX,iYY
        integer ia,ia3,ib,ib3
        integer iv,iv3,k,iegV,ipm
        real as
        character*2 iXXc2,iYYc2
        character*(charLenMAX) fileEssMod           
c
        real xyz(3)
        logical CONTROL
        integer kanalp
c
        kanalp = kanalRunOut
        CONTROL = .true.
c
        if(CONTROL)then
        write(kanalp,*)'getModeXYZ start: '
        write(kanalp,*)'write XYZ of 1,..,nMode:',nMode
        end if
c init
cx        nProtBBatRes = 4  ! number of protein BackBone atoms per residue
        ipmMX_OPT = 2 ! =4 number of files per one Mode !MAX=nModeSnapMAX=4
c
        do i=1,nModeXYZmAX*nModeSnapMAX*3*nAtomMAX
         atEssModeXYZ(i)=0.0
        end do!i
c
cx        kanalModXYZ = 37
c
c calculate essential Modes 1,..,nMode
        do im = 1,nMode
        iegV = im
        iXX = im 
        iXXc2 = '00'
c define fileName:
         call conv_numb_to_char(iXX,iXXc2)
         if(iXX .le. 9) then
         iXXc2 = '0'//iXXc2(1:1)
          else
           if(iXX .le. 99 )then
            iXXc2 = iXXc2(1:2)
                 else
                 write(kanalp,*)
     &           'getModeXYZ: ERROR!:number of iXX > 99 '
                 stop
                end if
               end if
c
c make number of files per one Mode
        do ipm = 1,ipmMX_OPT
c define amplitude scale
        if(ipm .le. ipmMX_OPT/2)then
        as = (ipm-1)*2.0/ipmMX_OPT -1.0 
        else
        as = ipm*2.0/ipmMX_OPT -1.0
        end if
c
c write file mdEssMode.XX.YY.pdb     : XX=nMode, YY=snap Per one mode
c define file name 
        iYY = ipm
        iYYc2 = '00'
         call conv_numb_to_char(iYY,iYYc2)
c
         if(iYY .le. 9) then
         iYYc2 = '0'//iYYc2(1:1)
          else
           if(iYY .le. 99 )then
            iYYc2 = iYYc2(1:2)
                 else
                 write(kanalp,*)
     &           'getModeXYZ: ERROR!:number of iYY > 99 '
                 stop
               end if
             end if
c
        fileEssMod = 'mdEssModeXYZ.'//iXXc2//'.'//iYYc2//'.pdb'
        if(nMolNameLett .ge. 1)then
        fileEssMod = molNameArgLine(1:nMolNameLett)//
     &                 '.mdEssModeXYZ.'//iXXc2//'.'//iYYc2//'.pdb'
        end if
c
        open(unit=kanalModXYZ, file=fileEssMod, form='formatted',
     &       status='unknown')
c
        write(kanalModXYZ,*)fileEssMod,' iegV=',iegV,
     &  ' lambda(A^2): ',lambdaMPF(iegV),' as=',as,
     &  ' xyz = <xyz> + as*eigVect(*,iegV)'
c
        do ib = 1,nAtomProtBB
        ib3=3*ib-3
c define eigVector:
        if(essModeTypeOPT .eq. 2) iv=ib  ! all protBB atom in MPF list
        if(essModeTypeOPT .eq. 1) iv = int((ib-1)/nProtBBatRes) + 1 ! only CA atoms in  MPF list
        iv3=3*iv-3 
c
        ia=atomNumbProtBBList(ib)
        ia3=3*ia-3
        do k=1,3
        xyz(k) = atomXYZtraAve(ia3+k) + 
     &  as*sqrt(lambdaMPF(iegV))*eigVectMPF(iv3+k,iegV)
        end do!k   
c
         write(kanalModXYZ,7071)"ATOM  ",
     &   atomNumb(ia),atomName(ia),resName(ia),chName(ia),resNumb(ia),
     &   (xyz(k),k=1,3), atomQ0(ia)
c
        end do !ib
c
          iStatus = iStatus + 1
        write(kanalPStat,*)mStatus,iStatus,messageStatus
     &  ,' file:',fileEssMod, ' is written ...'
c
        close(kanalModXYZ)
c
        end do !ipm
        end do !im
c
	return
7071    format(a6,1x,i4,2x,a4,a4,a1,i4,4x,3f8.3,7x,f8.5) ! PDB rq  
        end
