c writeMolXYZfile
c 
c Yuri Vorobjev  2005                           
c
	subroutine writeMOLXYZfile(iv,fileMOLX,nfile,wAtomXYZ)
c
c write PDB in file=fileMOLX
c iv=0/1 
c iv=1 write snaps sequence file name = fileMOLX.nfile.pdb
c
        include 'xyzPDBsize.h'
c
        include "output.h"
        include "kanalUse.h"
        include "statusMessg_mDyn.h"
c
        include 'xyzPDBinfo.h'
        include 'enForce.h'
        include 'filedat.h'  
        include 'mdAtomXYZvel.h'
        include 'optionPar.h'
        include 'engXYZtra.h'
        include 'compGeoEF.h'
c
        real wAtomXYZ(*)
        character*60 fileMOLX
        integer iv,nfile
c  local
        integer nrecMoL
        character*4 nRecMoLC4,nRecMoLC
        character*60 fileMoL,fileMOLX0 
        integer nLfileMoL0,nLfileMoL
        integer kanalMoLXYZ
        integer kanalp
        integer ich,k,k3,j
        logical doWRpdb
        logical CONTROL
c
        kanalp = kanalRunOut
        CONTROL = .true.     
c
        if(CONTROL)then
        write(kanalp,*)'writeMoLEngXyz start: '
        end if
c
        nrecMoL = nfile
c
        nLfileMoL0 = 60   ! n letters in file name
        nRecMoLC = '0000'
c
        kanalMoLXYZ = 32
c
        fileMOLX0 = fileMOLX
        call clean_spacelf(fileMOLX0,nLfileMoL0,fileMoL,nLfileMoL)
c
c initial fileName
        fileMOLX0 = fileMoL(1:nLfileMoL)
c define extended fileName:       
         call conv_numb_to_char(nRecMoL,nRecMoLC)
c
         if(iv .ne. 0)then 
         if(nRecMoL .le. 9) then
         nRecMoLC4 = '000'//nRecMoLC(1:1)
          else 
           if(nRecMoL .le. 99 )then
            nRecMoLC4 = '00'//nRecMoLC(1:2)
             else
              if(nRecMoL .le. 999 )then
               nRecMoLC4 = '0'//nRecMoLC(1:3)
               else 
                if(nRecMoL .le. 9999 )then
                 nRecMoLC4 = nRecMoLC(1:4)
                else 
                 write(kanalp,*)
     &           'mdRun: ERROR!:number of nRecMoL > 9999 '
                 stop
                end if
               end if
              end if
             end if
c
        fileMOLX0 = fileMoL(1:nLfileMoL-4)//nRecMoLC4//'.pdb'
c
        end if !iv
c
        kanalMoLXYZ = 32
        open(unit=kanalMoLXYZ, file=fileMOLX0, form='formatted',
     &       status='unknown')
c
          write(kanalMoLXYZ,'(a8)')'########'
          write(kanalMoLXYZ,'(a18,i4)')'$ENELIG: nOrient: ',nfile
          write(kanalMoLXYZ,'(a12,f14.5)')'eVbondDef  :',eVbondDef  
          write(kanalMoLXYZ,'(a12,f14.5)')'eVangDef   :',eVangDef  
          write(kanalMoLXYZ,'(a12,f14.5)')'eImpDef    :',eImpDef  
          write(kanalMoLXYZ,'(a12,f14.5)')'eTorsDef   :',eTorsDef  
          write(kanalMoLXYZ,'(a12,f14.5)')'engVDWR1   :',engVDWR1  
          write(kanalMoLXYZ,'(a12,f14.5)')'hBHxYeng128:',hBHxYeng128  
          write(kanalMoLXYZ,'(a12,f14.5)')'engCOULR1  :',engCOULR1  
          write(kanalMoLXYZ,'(a12,f14.5)')'engCOULR2  :',engCOULR2  
          write(kanalMoLXYZ,'(a12,f14.5)')'restr1Eng  :',restr1Eng  
          write(kanalMoLXYZ,'(a12,f14.5)')'eRstHW1ML  :',restr1MHWEng 
          write(kanalMoLXYZ,'(a10,f14.5)')'eRstA2Eng:',restrDistA2Eng 
          write(kanalMoLXYZ,'(a10,f14.5)')'eCompcEng:',compactGeoEn 
          write(kanalMoLXYZ,'(a12,f14.5)')'eGeoDef    :',eGeoDef  
          write(kanalMoLXYZ,'(a12,f14.5)')'engCOUL    :',engCOUL  
          write(kanalMoLXYZ,'(a12,f14.5)')'engSolv    :',molSolEn  
          write(kanalMoLXYZ,'(a10,f14.5)')'engSolHP :',hpSoLEnergy
          write(kanalMoLXYZ,'(a10,f14.5)')'bornPolz :',bornPolzEng
          write(kanalMoLXYZ,'(a10,f14.5)')'engSolEFS:',watBrgEnergySoL
          write(kanalMoLXYZ,'(a10,f14.5)')'engSolTot:',solvMolEngF3
          write(kanalMoLXYZ,'(a12,f14.5)')'engPOTENT  :',engPOTENT 
          write(kanalMoLXYZ,'(a10,f14.5)')'kinEng   :',kinEng
          write(kanalMoLXYZ,'(a10,f14.5)')'Temp(K)  :',tempT0
          write(kanalMoLXYZ,'(a10,f14.5)')'engTotal :',engTotal 
          write(kanalMoLXYZ,'(a7)')'$PDBXYZ '
c
          write(kanalMoLXYZ,'(a12)')'REMARK: PDB: '
c
           do ich=1,nChain
           do k = startAtInCha(ich),stopAtInCha(ich) 
           k3 = 3*k-3
c
           doWRpdb=.true.
c
           if(doWRpdb)then
           write(kanalMoLXYZ,7071)"ATOM  ",
     &     atomNumb(k),atomName(k),resName(k),chName(k),resNumb(k),
     &     (wAtomXYZ(k3+j),j=1,3), atomQ0(k)
           end if !doWRpdb
           end do !k
c TERminal chainLine
           write(kanalMoLXYZ,7071)"TER   ", stopAtInCha(ich)
           end do !ich
c END      line
           write(kanalMoLXYZ,7071)"END   " 
c
           close(kanalMoLXYZ)
c
         return
7071    format(a6,1x,i4,2x,a4,a4,a1,i4,4x,3f8.3,7x,f8.5) ! PDB rq
	 end
