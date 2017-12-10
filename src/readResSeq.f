c readIn molResSeq Info file
c
c Yuri Vorobjev, 2002
c
	subroutine readResSeq(mtinfile,nres,resList,
     &                    nLoop,resLoopFlag,startEndResLoop)

c mtinfile - resSeq Input file name
cEXAMPle:
c$MOLEC
c$SEQ A 22
cALAN VAR  LYS  ARG  TYR* ALA* PHE* LEU* HIS* PRO  MET  ILE  ASP
cGLU  SER  THR  TYR  CYS  ASN  GLN  TRP  GLUN
cENDSEQ
c  TYR*  (or TYR * ) - Ltag=loopResTag=* in pos 5 
c nres - nresidues
c resList(k), k=1,..,nres resList(k)=extendedResName=name(1:4)//Ltag//chTag
c nLoop , number of loops with unknown XYZ
c resLoopFlag(ir), =0,1,2,3,.., belong to core,Loop1,Loop2,..
c startEndResLoop(2*nLoop) - nResSatrt,End of Loop
c 
	implicit none
        character*(*) mtinfile 
        integer nres
        character*6 resList(*)     
        integer nLoop,resLoopFlag(*)
        integer startEndResLoop(*)
c
        include "charStringSiz.h"
        include "output.h"
        include "kanalUse.h"
        include "statusMessg_mDyn.h"
c local
        character*80 linein
        character*(charLenMAX) molName
        integer i,j,k
        character*1 chTag
        integer nin,nresch,nlseq,nost
        integer nlsr,nfirst,nlast,nmol
        integer nresinline,nls,nrch
        integer kanalp
        logical CONTROL,CONTROL0
        logical molBlock,seqBlock
        character*1 LOOPmark
        logical loopRes,loopResP
       
         CONTROL = .false. 
         CONTROL0 = .true.
         kanalp = kanalRunOut
c
         write(kanalp,*)'Start readResSeq:'
c 
c open files
        open(unit=11,file=mtinfile, form='formatted', status='old')
c
c initialize
       loopRes = .false.
       LOOPmark = '*'
    	nres = 0
        nresinline = 13  ! numbers of residues in 1 line of mtinfile
c
c read mtinfile
        nLoop = 0
        nmol = 1
        nfirst = 1
        molBlock = .false.
        nin = 0
100     read(11,'(a80)', end=101)linein
c
        if(CONTROL)then
        write(kanalp,'(a3,a80)')'W1:',linein
        end if !control
c
        if(linein(1:1) .eq. '!' .or. linein(1:1) .eq. '#')goto 100
        nin = nin + 1  ! number of informative line
c
        if(linein(1:6) .eq. '$MOLEC' ) then
        molBlock = .true.
        read(linein,'(7x,a20)') molName
        goto 100
        end if
c
c start read seqBlock
        if(molBlock)then
        if(linein(1:4) .eq. '$SEQ' )then
        seqBlock = .true.
        read (linein, '(5x,a1,1x,i5)') chTag, nresch    
c
        if(CONTROL)then
        write(kanalp,*) 'W2:chTag, nresch:',chTag, nresch
        end if!CONTROL
c
        end if ! $SEQ  
        end if !molBlock
c
c read SEQ of nresch ReS
c
        if(seqBlock)then
        nlseq=int(nresch/nresinline)
        nost = nresch - nlseq*nresinline 
        nlsr = nlseq
        if(nost .gt. 0 ) nlsr = nlseq + 1
c init counters
        nls = 0
        nrch = 0
c read next nlsr lines of seqBlock
200     read(11,'(a80)', end=101)linein
        if(linein(1:1) .eq. '!' .or. linein(1:1) .eq. '#')goto 200

        if(CONTROL) then
        write(kanalp,*)'W3:',linein
        end if !C
c
        nls = nls + 1
c
        if(nls .le.  nlseq) then        
        read(linein, '(13(a5))')
     &  (resList(nres+nrch+k), k=1,nresinline)
        nrch = nrch + nresinline
        goto 200
        end if !nls .le.  nlseq
c
        if(CONTROL)then
        write(kanalp,*)'W3b: nrch:',nrch
        end if !C
c
c read last line
        if(nlsr .eq. nlseq + 1)then
        read(linein, '(13(a5))')
     &  (resList(nres+nrch+k), k=1,nost)
        nrch = nrch + nost      

        if(CONTROL)then
        write(kanalp,*)'W3c: nrch:',nrch
        end if !C
c seqBlock ending
        seqBlock = .false.
c
        goto 100
c
        end if !nlsr .eq. nlseq + 1
c
        end if ! seqBlock
c
        if(linein(1:7) .eq. '$ENDSEQ' )then
c forced seqBlock ending
        seqBlock = .false.
c accumulate        
        nres = nres + nresch
        nlast = nres
c
        if(CONTROL)then
        write(kanalp,*)'W4:Mol,nrch:',nmol,nrch
	call printSEQ(nfirst,nlast,nresinline,resList)
c
        if(nrch .ne. nresch ) then
        write(kanalp,*) 'moltop.f:ERROR!:readRESseq:',
     &  ' nrch, nresch: ',  nrch, nresch
        write(kanalp,*)'W6:Mol:',nmol
	call printSEQ(nfirst,nlast,nresinline,resList)
        stop
        end if
        end if !CONTROL
c
c make extended resName
        do k = 0,nresch-1
        resList(nfirst+k) = resList(nfirst+k)(1:5)//chTag
        end do
c
        if(CONTROL0)then
        write(kanalp,*)'readResSeq: nmol,nresch:',nmol,nresch
	call printSEQ(nfirst,nlast,nresinline,resList)
        end if
c
        end if ! ENDSEQ
c
        if(linein(1:7) .eq. '$ENDMOL' )then
c forced molBlock ending
        molBlock = .false.
        end if !$ENDMOL
c
        goto 100

101     continue

c find loopRes
       if(CONTROL)then 
       write(kanalp,*)'readResSeq: W9: find Loops'
       end if !Contr
c
       nLoop = 0
       resLoopFlag(1) = 0
       if(resList(1)(5:5) .eq. LOOPmark .or.
     &    resList(1)(4:4) .eq. LOOPmark) loopRes=.true.

       if(loopRes)then
       nLoop = nLoop + 1
       resLoopFlag(1) = nLoop
       end if
c
       loopResP = loopRes
c
       do k = 2,nres
       resLoopFlag(k) = 0
       loopRes = .false.
       if(resList(k)(5:5) .eq. LOOPmark .or.
     &    resList(k)(4:4) .eq. LOOPmark) loopRes=.true.
       
       if(loopRes)then
       if(.not. loopResP)then
       nLoop = nLoop + 1
       startEndResLoop(2*nLoop-1) = k       !startLoop
       end if
       resLoopFlag(k) = nLoop
       else
       if(loopResP) startEndResLoop(2*nLoop) = k-1   ! endLoop
       end if
c
       loopResP = loopRes 
       end do !k
c
       if(CONTROL0)then
       write(kanalp,*)'readResSeq: nLoop=',nLoop
       write(kanalp,*)'resLoopFlag():'
       write(kanalp,'(13(i3,4x))')(resLoopFlag(k),k=1,nres)
       write(kanalp,*)'startEndResLoop:'
       write(kanalp,*) (startEndResLoop(k),k=1,2*nLoop)
       end if
c
        close(unit=11)
c
        return
        end         
c
	subroutine printSEQ(nfirst,nlast,ninline,resList)
c print resList
        implicit none
        include "charStringSiz.h"
        include "output.h" 
        integer nfirst,nlast,ninline
        character*6 resList(*)
c
        integer nl,kanalp
        integer nost,k 
        integer nseq,np,l

        kanalp = kanalRunOut 
        nseq=nlast-nfirst+1
c
        np = nfirst-1
        nl = int(nseq/ninline)
        nost = nseq - nl*ninline
c        
        write(kanalp,'(a26,2i6)')
     &  'SEQUENCE:nlast, nfirst:', nfirst, nlast
        do l=1,nl
        write(kanalp,'(13(a6,1x))')
     &  (resList(np+k), k=1,ninline)
        np = np+ninline
        end do
c
        if( nost .gt. 0 )then
        write(kanalp,'(13(a6,1x))')
     &  (resList(np+k), k=1,nost)
        end if
c        
        return
        end
c
