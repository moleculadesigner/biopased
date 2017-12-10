c writeLigEngXyz
c 
c Yuri Vorobjev  2003                           
c
	subroutine wrLigEngXyz(fileLigN,nfile)
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
        include 'ligInfo.h'
        character*(charLenMAX) fileLigN
        integer nfile
c
c  local
        integer nrecLig
        character*4 nRecLigC4,nRecLigC
        character*(charLenMAX) fileLig 
        integer nLfileLig0,nLfileLig
        integer kanalLigEng
        integer kanalp
        integer ich,k,k3,j
        logical doWRpdb
        logical CONTROL
c
        kanalp = kanalRunOut
        CONTROL = .true.     
c
        if(CONTROL)then
        write(kanalp,*)'writeLigEngXyz start: '
        end if
c
        nrecLig = nfile
c
        nLfileLig0 = charLenMAX   ! n letters in file name
        nRecLigC = '0000'
c
        kanalLigEng = 32
c
        call clean_spacelf(fileLigN,nLfileLig0,fileLig,nLfileLig)
c
c define fileName:       
         call conv_numb_to_char(nRecLig,nRecLigC)
c
         if(nRecLig .le. 9) then
         nRecLigC4 = '000'//nRecLigC(1:1)
          else 
           if(nRecLig .le. 99 )then
            nRecLigC4 = '00'//nRecLigC(1:2)
             else
              if(nRecLig .le. 999 )then
               nRecLigC4 = '0'//nRecLigC(1:3)
               else 
                if(nRecLig .le. 9999 )then
                 nRecLigC4 = nRecLigC(1:4)
                else 
                 write(kanalp,*)
     &           'mdRun: ERROR!:number of nRecLig > 9999 '
                 stop
                end if
               end if
              end if
             end if
c
        fileLigN = fileLig(1:nLfileLig-4)//nRecLigC4//'.pdb'
c
        kanalLigEng = 32
        open(unit=kanalLigEng, file=fileLigN, form='formatted',
     &       status='unknown')
c
          write(kanalLigEng,'(a8)')'########'
          write(kanalLigEng,'(a18,i4)')'$ENELIG: nOrient: ',nfile
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
          write(kanalLigEng,'(a12)')'REMARK: PDB: '
c
           do ich=1,nChain
           do k = startAtInCha(ich),stopAtInCha(ich) 
           k3 = 3*k-3
c
           doWRpdb=.true.
           doWRpdb=.false.
           if(resName(k) .eq. resName(atomInLigList(1))) doWRpdb=.true.
c
           if(doWRpdb)then
           write(kanalLigEng,7071)"ATOM  ",
     &     atomNumb(k),atomName(k),resName(k),chName(k),resNumb(k),
     &     (atomXYZ0(k3+j),j=1,3), atomQ0(k)
           end if !doWRpdb
           end do !k
c TERminal chainLine
           write(kanalLigEng,7071)"TER   ", stopAtInCha(ich)
           end do !ich
c END      line
           write(kanalLigEng,7071)"END   " 
c
           close(kanalLigEng)
c
         return
7071    format(a6,1x,i4,2x,a4,a4,a1,i4,4x,3f8.3,7x,f8.5) ! PDB rq
	 end
