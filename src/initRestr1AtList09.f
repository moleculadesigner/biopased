c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c                                                                   *
c  Yury Vorobjev 2003, 2004, 2009                                   *
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c defines:
c restr1atomList(ia) for a given list of nRestr1Seg,resEndRestr1(2* nRestr1Seg)
c
	subroutine initRestr1AtList09(nRestr1Seg,resEndRestr1,
     &                 resResSegKhConst,resRest1AtLine,
     &                 nRestr1atom,restr1atomList,refAtomPos,
     &                 restr1AtConst,restr1AtConstList)   
c
        include 'xyzPDBsize.h'
        include 'xyzPDBinfo.h'
        include 'xyzPDBcrd.h'
        include 'filedat.h'
c
        include "output.h"
        include "kanalUse.h"
        include "statusMessg_mDyn.h"
c
        integer nRestr1Seg,nRestr1atom
        integer resEndRestr1(*)
        character*40 resRest1AtLine(*)
        integer restr1atomList(*)
        real refAtomPos(*)
        real restr1AtConst
	real resResSegKhConst(*)
        real restr1AtConstList(*)
clocal
        integer nWordMAX
        parameter (nWordMAX=8)
        integer nWord,wLengthinLine(nWordMAX)
        character*4 wordInLine(nWordMAX)
c
        integer nPBBatMAX
        parameter (nPBBatMAX=4)
        character*2 pBBatNameList(nPBBatMAX)
c
        integer nNASP3atMAX
        parameter (nNASP3atMAX=3)
        integer nNASP1atMAX
        parameter (nNASP1atMAX=1)
        character*1 naPSnameList3(nNASP3atMAX)
        character*1 naPSnameList1(nNASP1atMAX)
        character*1 anc1,anc3
c        
        data pBBatNameList/'CA','C ','O ','N '/
c atomCharacter3 or 1 for atoms in the RibPh backbone
        data naPSnameList3/ "'","T","P"/
        data naPSnameList1/"P"/
       
c
        character*4 atNameia
        integer kanalp
        integer i,ia
        integer nLtInLineMX
        integer ip,ip2,ips1,ips2
        integer i3,ia3
        integer ivar,ipba
        integer nltA,nlt,nlt4
        logical includeInA1List
        logical OPT_HatomFree
        logical CONTROL
c        
c initialize:
        kanalp =  kanalRunOut
        CONTROL = .true.
c
        OPT_HatomFree = .true.    ! h-atoms is not harmonically restricted
        nLtInLineMX = 40   ! max letters in line to analysed
        nlt4=4
c
c all realAtom(ia) [protCORE + LOOP] are moving 
         nRestr1atom=0    
c
         if(CONTROL)then
         write(kanalp,*)'define initRestr1AtList ...'
         write(kanalp,*)'pBBatNameList:',pBBatNameList
         end if
c
         if(nRestr1Seg .ge. 1)then
         do ip = 1,nRestr1Seg
         ip2 = ip*2-1
         ips1 = resEndRestr1(ip2)
         ips2 = resEndRestr1(ip2+1)
c
         if(CONTROL)then
         write(kanalp,*)'define initRestr1AtList ...'
         write(kanalp,*)
     &   'ip:',ip,' resRest1AtLine(ip):',resRest1AtLine(ip)
c
         if(ips1 .lt. 0 .or. ips1 .gt. nres 
     &      .or. ips2 .lt. 0 .or. ips2 .gt. nres )then
         write(kanalp,*)
     &  'ERROR!initRestr1AtList: WRONG! resNumb of restSegm:',ip,
     &  ' resNst, resNfin:', ips1,ips2
         write(kanalp,*)' correct file: restrA1.inp ',
     &  ' are out of 1-nREs range !!! nres = ',nres
          write(kanalPStat,*)mError,
     &  ' restraint.inp error file -r restrA1.inp '
         stop
          end if !
          end if !Control
c
c read subroutine char_line_word_analys(n0,line0,nw,word,wlength) 
c read resRest1AtLine() 
c
         call char_line_word_analys(nLtInLineMX,resRest1AtLine(ip),
     &        nWord,wordInLine,wLengthinLine)         
c
         if (nWord .eq. 0) ivar=0
         if (nWord .eq. 1) then
          if(wordInLine(1)(1:wLengthinLine(1)) .eq. 'ALL') ivar = 1
          if(wordInLine(1)(1:wLengthinLine(1)) .eq. 'PBB') ivar = 2      !Prot BackB
          if(wordInLine(1)(1:wLengthinLine(1)) .eq. 'NAB') ivar = 4      !NA Base 
         end if
c
        if (nWord .gt. 1) then 
         if(wordInLine(1)(1:wLengthinLine(1)) .eq. 'PBB') ivar = 3
        end if
c
         do ia = startAtInRes(ips1),stopAtInRes(ips2)
         includeInA1List = .false.
         if(ivar .eq. 0 .or. ivar .eq. 1)includeInA1List = .true.
c
         if( ivar .eq. 2 )then   !ALL ProtBackBoneAtoms
         do ipba = 1,nPBBatMAX 
         nlt = LEN(pBBatNameList(ipba))
         if(pBBatNameList(ipba)(1:nlt) .eq. atomName(ia)(1:nlt))
     &   then 
         includeInA1List = .true.
         goto 2001
         end if !pBBatNameList
         end do!ipba
c
2001     continue
         end if !iv=2
c
         if( ivar .eq. 3 )then   !name selected backBoneAtoms
         do ipba = 2,nWord
         nlt = wLengthinLine(ipba)
         call clean_spacelf(atomName(ia),nlt4,atNameia,nltA)
c
cx         if(CONTROL)then
cx         write(kanalp,*)'nlt,nltA:',nlt,nltA, ' at:',atomName(ia)
cx         end if !C
c
         if(nlt .eq. nltA .and.
     &   (wordInLine(ipba)(1:nlt) .eq. atomName(ia)(1:nlt)))
     &   then
         includeInA1List = .true.
         goto 2002
         end if !(wordInLine(ipba)
         end do!ipba
c
2002     continue
         end if ! iv=3
c
c NAcid analysis: include in RESTRAINList NA BASE atoms
         if(ivar .eq. 4) then
         anc1 = atomName(ia)(1:1)
         anc3 = atomName(ia)(3:3)
         includeInA1List = .true.
c
         do i = 1,nNASP3atMAX
         if(anc3 .eq. naPSnameList3(i))includeInA1List = .false.
         end do !i3
         do i = 1,nNASP1atMAX
         if(anc1 .eq. naPSnameList1(i))includeInA1List = .false.
         end do !i1 
         end if ! ivar=4         
c
         if(includeInA1List) then
         nRestr1atom = nRestr1atom + 1
         restr1atomList(nRestr1atom) = ia    
cc define refAtomPos()
         ia3 = 3*ia - 3
	 i3=3*nRestr1atom-3
         refAtomPos(i3+1) = atomXYZ(ia3+1)
         refAtomPos(i3+2) = atomXYZ(ia3+2)
         refAtomPos(i3+3) = atomXYZ(ia3+3)
         restr1AtConstList(nRestr1atom)=resResSegKhConst(ip)
c
         if(OPT_HatomFree)then
         if(atomName(ia)(1:1) .eq. 'H')
     &   restr1AtConstList(nRestr1atom)=0.0
         end if
c
         end if !includeInA1List
c        
         end do !ia
         end do !ip
c
        if(CONTROL)then
        write(kanalp,*)'initRestr1AtList: nRestr1Seg:',nRestr1Seg
        write(kanalp,*)'nRestr1atom :',nRestr1atom
        write(kanalp,*)'i    restr1atomList: restr1AtConst: '
        do i = 1,nRestr1atom
        ia=restr1atomList(i)
        write(kanalp,'(2i6,1x,a4,1x,a4,i4,1x,f6.3)') 
     &  i,ia,atomName(ia),resName(ia),resNumb(ia),restr1AtConstList(i)
        end do !i
c
        end if !CONTROL
c
        write(kanalp,*)'initRestr1AtList: Finish:'
c
        end if !nRestr1Seg .ge. 1
c
	return
        end
c
        subroutine char_line_word_analys(n0,line0,nw,word,wlength)   
c extract words from line 
c line = aaa  bb ccc 
c
        implicit none
        integer n0
        character*(*) line0
        integer nw
        character*4 word(*)
        integer wlength(*)
c
        character*1 space
        integer i,rw,rwm,nlt
        integer kanalp
        logical CONTROL
c
        kanalp = 6
        CONTROL = .false.
        space =' '
c
        if(CONTROL)then
        write(kanalp,*)
     &  'char_line_word_analys: n0,line0:',n0,line0 
        do i=1,n0
        write(kanalp,*)'i:',i,' line0(i:i):',line0(i:i)
        end do !
        end if
c
        nw = 0
        if(n0 .eq. 0) return
c
        rw = 0        
        rwm = 0     
        nlt = 0
        nw = 0    
c
        do i = 1,n0
        rw = 0
        if( line0(i:i) .ne. space )rw = 1
c start new word
        if( rwm .eq. 0 .and. rw .eq. 1)then
        nw = nw+1
        nlt=1
        word(nw)(nlt:nlt) = line0(i:i)
        end if!
c continue word
        if( rwm .eq. 1 .and. rw .eq. 1)then
        nlt=nlt+1
        word(nw)(nlt:nlt) = line0(i:i)
        end if
c finish the word
        if( rwm .eq. 1 .and. rw .eq. 0)then
        wlength(nw) = nlt
        end if !c
        rwm = rw
c 
        if(CONTROL)then
        write(kanalp,*)'i:',i, ' line0(i:i):',line0(i:i)
        write(kanalp,*)
     &  'rw,rwm:',rw,rwm,' nw:',nw,' word(nw):',word(nw)(1:nlt)
        end if!C
c
        end do !i   
c
        if(CONTROL)then
        write(kanalp,*)'char_line_word_analys: nWord:',nw
        if (nw .ge. 1)then
        do i=1,nw 
        write(kanalp,*)'iw,wlength,word:',i, wlength(i),word(i)
        end do 
        end if !nw
        end if !C
c
        return
        end
