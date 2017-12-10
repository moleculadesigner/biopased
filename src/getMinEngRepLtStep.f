c select snaps of MinEng for ALL Replics
c
c USING data of replicaEx.h
c
	subroutine getMinEngRepLtStep(nMdCycle)
c
        include 'xyzPDBsize.h' 
        include 'xyzPDBinfo.h'
        include 'ligInfo.h'
        include 'optionPar.h' 
        include 'replicaEx.h'
        integer nMdCycle
c
        include "output.h"
        include "kanalUse.h"
        include "statusMessg_mDyn.h"
c local
        integer nLfilePdb0,nLfileMinPdb
        integer nRpL,i,is,j,k,k3,ich
        integer indx0,indx
        integer engTermMX  
        integer rL0eShift
        integer xShift,xmShift
        integer na3,isMin
        integer kanalPDBmin
        integer nEtermRef
        character*4 nMdCycle4C,nRpL4C
        character*(charLenMAX) pdbMinEngFile
        logical doWRpdb
        real big
        logical OPT_minEngSnap
        integer kanalp
        logical CONTROL
c
        kanalp=kanalRunOut
        CONTROL = .true.
c
        OPT_minEngSnap = .true.
        nEtermRef = 2       ! ePot' = ePot-eGeo
        engTermMX=engTermMAX
	big = 1.0e10
c
        call conv_numb_to_char00n(nMdCycle,nMdCycle4C)
c 
        nRpL=0
c
100     nRpL = nRpL+1
        if(nRpL .gt. nRepLic)goto 101
c
        rL0eShift=engTermMAX*(nRpL-1)         
        indx = rL0eShift

c write XYZFileMinEng
         nLfilePdb0=60
          call clean_spacelf
     &       (mMdRepLpdbf(nRpL),nLfilePdb0,pdbMinEngFile,nLfileMinPdb)
c
c          call conv_numb_to_char00n(nRpL,nRpL4C)
c
          pdbMinEngFile = 
     &    pdbMinEngFile(1:nLfileMinPdb)//'.eMin.'//nMdCycle4C//'.pdb'
c
           if(CONTROL)then
           write(kanalp,*)'getMinEngRepLTra: wrFile:',pdbMinEngFile
           end if
c
          kanalPDBmin = 61
          open(unit=kanalPDBmin, file=pdbMinEngFile, form='formatted',
     &       status='unknown')
c        
           write(kanalPDBmin,'(a33,i6,a10,i6)')
     &   '$REMARK: MdMinEngTstep xyzRepL:   ',nRpL,' nMDcylce:',nMdCycle
          write(kanalPDBmin,'(a8,2x,a6)')'$ENERGY: ',' :Kcal'
          write(kanalPDBmin,'(a23,2x,i4)')
     &    '$nRepEngTern: engPOTENT1:',nEtermRef
          write(kanalPDBmin,'(a10,f14.5)')
     &    'eGeoDef  :',mdEngRepL0eMinRun(indx+3) 
          write(kanalPDBmin,'(a10,f14.5)')
     &    'engVDWR1 :',mdEngRepL0eMinRun(indx+4)
          write(kanalPDBmin,'(a12,f14.5)')
     &    'hBHxYeng128:',mdEngRepL0eMinRun(indx+5)
          write(kanalPDBmin,'(a10,f14.5)')
     &    'engCOULR1:',mdEngRepL0eMinRun(indx+6)
          write(kanalPDBmin,'(a10,f14.5)')
     &    'restrEng :',mdEngRepL0eMinRun(indx+8)
          write(kanalPDBmin,'(a10,f14.5)')
     &    'engSolv  :',mdEngRepL0eMinRun(indx+7)
          write(kanalPDBmin,'(a10,f14.5)')
     &    'engPOTENT:',mdEngRepL0eMinRun(indx+1)
          write(kanalPDBmin,'(a10,f14.5)')
     &    'engPOTENT1:',mdEngRepL0eMinRun(indx+2)
          write(kanalPDBmin,'(a10,f14.5)')
     &    'Temp(K)  :',mdEngRepL0eMinRun(indx+9)
c
          write(kanalPDBmin,'(a12)')'REMARK: PDB: '
c
c update mdEngRepL0MixMin(eTerm,nRpL)
           do k=1,engTermMAX
           mdEngRepL0MixMin(k,nRpL)=mdEngRepL0eMinRun(indx+k)
           end do!k
c
           na3 = 3*natom
           xmShift= na3*(nRpL-1)                                   
c update replicaXYZV
           xShift = na3*(nRpL - 1)
           do k = 1,na3
           atomXYZRepL(xShift+k) = atomXYZRepL0eMinRun(xmShift+k)
           atomVelRepL(xShift+k) = atomVelRepL0eMinRun(xmShift+k)
           end do !k
c writeToFile
           do ich=1,nChain
           do k = startAtInCha(ich),stopAtInCha(ich)
           k3 = 3*k-3
c
           doWRpdb=.true.
           if(OPT_doLigDock .eq. 1)then
           doWRpdb=.false.
           if(resName(k) .eq. resName(atomInLigList(1))) doWRpdb=.true.
           end if !OPT_doLigDock .eq. 1
c
           if(doWRpdb)then
           write(kanalPDBmin,7071)"ATOM  ",
     &     atomNumb(k),atomName(k),resName(k),chName(k),resNumb(k),
     &     (atomXYZRepL0MixRun(xmShift+k3+j),j=1,3), atomQ0(k)
           end if !doWRpdb
           end do !k
c
c TERminal chainLine
           if( (resName(stopAtInCha(ich))(1:3) .ne. 'HOH' .and.
     &        resName(stopAtInCha(ich))(1:3) .ne. 'ION') .or.
     &        ich .eq. nChain)
     &     write(kanalPDBmin,7071)"TER   ", stopAtInCha(ich)
           end do !ich
c END      line
           write(kanalPDBmin,7071)"END   "
c
           close(kanalPDBmin)           
c 
        goto 100
101     continue
c        
        return
7071    format(a6,1x,i4,2x,a4,a4,a1,i4,4x,3f8.3,7x,f8.5) ! PDB rq 
        end
