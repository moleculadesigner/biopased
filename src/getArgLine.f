      subroutine getArgLine
C
      implicit none
c
c
      include "charStringSiz.h"
      include "filedat.h"
      include "output.h"
      include "kanalUse.h" 
      include "statusMessg_mDyn.h"
c
      character*(charLenMAX) arg
      integer i,n,kanalp
      logical opErr,opOut,userErr,userOut
c
c        ... temp for each of the whitespace delimited command line words
      integer iarg
c        ... arg pointer, final number of arguments
      integer indx
      integer iargc
c
      opErr = .false.
      opOut = .false.
      userOut = .false.
      userErr = .false.
c
c     --- default file names ---
c
      mdynParFile = './MdynPar.inp'
      xyzInpFile = './molec.pdb'
      loopFile = './moveRes.inp'
      saProtFile = './SAprotocol.inp'
      restrA1type = './restrAt1.inp'
      restrA2type = './distRestrA2.inp'
      molNameArgLine = "" 
      outputFile = "./mDynSB.out"
      outErrFile  = "./mDynSB.err"
c
      kanalp=kanalRunOut
c
        nMolNameLett = 0
c
c     --- default status of output: new
c
      iarg = 0
      indx = iargc()
c
      if (indx .eq. 0) goto 20
c
   10 continue
           iarg = iarg + 1
           call getarg(iarg,arg)
c
           if (arg .eq. '-i') then
                iarg = iarg + 1
                call getarg(iarg,mdynParFile)
           elseif (arg .eq. '-c') then
                iarg = iarg + 1
                call getarg(iarg,xyzInpFile)
           elseif (arg .eq. '-mv') then
                iarg = iarg + 1
                call getarg(iarg,loopFile)
           elseif (arg .eq. '-sa') then
                iarg = iarg + 1
                call getarg(iarg,saProtFile)
           elseif (arg .eq. '-r') then
                iarg = iarg + 1
                call getarg(iarg,restrA1type)
c
          elseif (arg .eq. '-d') then
                iarg = iarg + 1
                call getarg(iarg,restrA2type)
c
           elseif (arg .eq. '-o') then
                iarg = iarg + 1
                call getarg(iarg,outputFile)
                userOut=.true.
c
           elseif (arg .eq. '-er') then
                iarg = iarg + 1
                call getarg(iarg,outErrFile)
c
           elseif (arg .eq. '-mn') then
                iarg = iarg + 1
                call getarg(iarg,molNameArgLine)
c
          nMolNameLett = 0
          n = 80
          do i = 1,n
          if(molNameArgLine(i:i).ne.' ')then
          nMolNameLett = nMolNameLett + 1
          else
          goto 101
          end if
          end do !i
c
101       continue
c
           else
         write(kanalPstat,'(/,5x,a,a)') 'getArgLine: unknown flag: ',arg
           endif ! 1
c
       if (iarg .lt. indx) goto 10
c
20     continue
c
         if( .not. opOut)then
         if(.not. userOut)then
         if(nMolNameLett .ge. 1)
     &    outputFile = molNameArgLine(1:nMolNameLett)//".mDynSB.out"
         end if
         if(kanalp .ne. kanalSTD6)
     &    open(unit=kanalp,file=outputFile,form='formatted',
     &      status='unknown')
          opOut = .true.
          end if
c
          if(.not. opErr)then
          if(.not. userErr)then
          if(nMolNameLett .ge. 1)
     &   outErrFile = molNameArgLine(1:nMolNameLett)//".mDynSB.err"
          end if
         if(kanalPStat .ne. kanalSTD6)
     &   open(unit=kanalPStat,file=outErrFile,form='formatted',
     &      status='unknown')
          opErr = .true.
         end if
c
      if(indx .eq. 0)then
      write(kanalp,*)'ArgLine: ALL defalt fileNames are taken:'
      else 
      write(kanalp,*)'User defined fileNames are taken:'
      write(kanalp,9000)mdynParFile,xyzInpFile,loopFile,
     &                  restrA1type,restrA2type,saProtFile,
     &            molNameArgLine,outputFile,outErrFile
       end if
c
      return

 9000 format(/,5x, ' job file name are taken:',
     +   /,5x, 'usage: mdynSB -i inProtcol  -c inPDB -mv moveRes ', 
     +   ' -r inRestrainA1 -d inRestrainA2 -sa saProtocol ',
     +   ' -mn molName -o runOutFile -er errorFile',
     +   /,5x, 'inProtcol      =',a80,
     +   /,5x, 'inPDB          =',a80,
     +   /,5x, 'moveRes        =',a80,
     +   /,5x, 'inRestrainA1   =',a80,
     +   /,5x, 'inRestrainA2   =',a80,
     +   /,5x, 'saProtocol     =',a80,
     +   /,5x, 'molName        =',a80,
     +   /,5x, 'runOutFile     =',a80,
     +   /,5x, 'errorFile      =',a80, 
     +   /, '* * * * * * * * * * * * * * * * * * * * * * * * * ' )
c
      end
