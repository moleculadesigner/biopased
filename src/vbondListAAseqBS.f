c defines vbonded 123, 1234 neighbours
c from ZMATRIX pair12List
c Yuri Vorobjev 2002
c 
        subroutine vbondListAAseq(ivar,atomXYZ,realAtomFlag,
     &           natom,atomNumb,atomName,resName,chName,resNumb,
     &           nres,resNameRes,atomBlockName,
     &           atomNameEx,startAtInRes,
     &           nmoveatom,moveAtomList,moveFlag,
     &           pair12List,startPairL12,nPairL12,np12MAX,
     &           pair13List,startPairL13,nPairL13,np13MAX,
     &           pair14List,startPairL14,nPairL14,np14MAX,
     &           pair14VList,startPairL14V,nPairL14V,
     &           bond12List,nbond12,nbond12noH,
     &           trip123List,nTrip123,np123MAX,
     &           quar1234List,nQuar1234,np1234MAX,
     &           quarImp1234L,nImp1234,nImp1234MAX)
c
c INPUT:
c ivar - option to coollect bonds12(), triplets(), quartets(*)
c  ivar = 0  bonds etc. between all realAtomFlag()          - fullProtMd            
c  ivar = 1  bonds etc. between all i,j with equal realAtomFlag(i) = realAtomFlag(j) 
c                       in pairs,triplets etc. - for addLoop
c
c atomXYZ(),realAtomFlag(),
c natom,atomNumb,atomName,resName,chName,resNumb,
c nres,resNameRes,atomBlockName,atomNameEx,startAtInRes,
c nmoveatom,moveAtomList,
c pair12List,startPairL12,nPairL12,np12MAX
c
c pair12List(i2) - list of atoms 2 of 12 pair, symmetrical
c startPairL12(i1) - point to the first position in pair12List() for 
c                    the sequence of atoms2 for atom1=i1
c nPairL12(i1) - length of sequence for atoms2 
c
c correction:
c 1)all bond12List(),trip123List(),quar1234List(),quarImp1234L() are calculated 
c     between atoms with equal realAtomFlag()
c OUT:
c bond12List(), nbond12 - total pairs12,
c nbond12noH  - pairs12 heavy-heavyAt
c pair13List() - list of atoms 3 of 13 pair
c pair14List() - list of atoms 4 of 14 pair 
c pair14VList() - list of atoms 4 of 14 pair for VDW/COUL interactions
c trip123List() - list of atoms 1,2,3 in 123 triplets
c quar1234List() - list of atoms 1,2,3,4 in  1234 quarets
c quarImp1234L() - list of atoms 1,2,3,4 in 1234 improper quartets
c 
        implicit none
        integer ivar
        real atomXYZ(*)
        integer realAtomFlag(*)
        integer natom
        character*4 atomName(*)
        character*4 resName(*)
        character*1 chName(*)
        integer atomNumb(*),resnumb(*)
        integer nres
        character*4 resNameRes(*)
        character*4 atomBlockName(*) 
        character*8 atomNameEx(*)
        integer startAtInRes(*)
        integer nmoveatom,moveAtomList(*)
        integer moveFlag(*)
c
        include "charStringSiz.h"
        include "output.h"
        include "kanalUse.h"
        include "statusMessg_mDyn.h"
cRESULT:
        integer np12MAX,np13MAX,np14MAX
        integer np123MAX,np1234MAX
        integer pair12List(*),startPairL12(*),nPairL12(*)
        integer pair13List(*),startPairL13(*),nPairL13(*)
        integer pair14List(*),startPairL14(*),nPairL14(*)
        integer pair14VList(*),startPairL14V(*),nPairL14V(*)
        integer bond12List(*),nbond12,nbond12noH
        integer trip123List(*),nTrip123
        integer quar1234List(*),nQuar1234
        integer nImp1234,nImp1234MAX
        integer quarImp1234L(*)
c local variables
        integer i,j,k,k3,i2,i4
        integer i3,j3,ir,jr,jrn
        integer ia,ja,ka,ma,iam
        integer n12H,i12,i2nH,i2H
        integer sjr,fjr
        integer np12s
        integer nb12s,nb12s2
        integer np13,np13s
        integer np14,np14s
        integer np14V,np14Vs
        integer np123,np1234
        integer np123s,np1234s
        integer np123s2,np1234s3
        integer p12i,p12is,p12if
        integer p12j,p12js,p12jf
        integer p12k,p12ks,p12kf
        integer p14i
        real radi,radj,dij2 
        integer kanalp,tri
        logical addNew 
        logical start_ia
        logical find
        logical OPT_removePairInHCycBase
        logical CONTROL,CONTROL1
c
        kanalp=kanalRunOut
        CONTROL = .true.   !control print
        CONTROL1= .false.
        OPT_removePairInHCycBase = .true.  ! remove 14LV HCycle system from the List
        
        if(CONTROL)then 
        write(kanalp,*)'Start vbondListAAseq:'
        end if
c
c initialize
        np13s = 0
        np14s = 0
        np14Vs = 0
        nb12s = 0
        nb12s2 = 0
        np123s = 0
        np1234s = 0
        np123s2 = 0 
        np1234s3 = 0
        addNew = .false.
        find = .false.
c numbers of record in Pair12List()
        np12s = startPairL12(natom) + nPairL12(natom) - 1
c
        do i=1,np12MAX
        bond12List(i)=0
        end do !i
        do i=1,np13MAX
        pair13List(i) = 0
        end do !i
        do i=1,np14MAX
        pair14List(i) = 0
        pair14VList(i) = 0
        end do !i
        do i=1,np123MAX
        trip123List(i) = 0
        end do !i
        do i=1,np1234MAX
        quar1234List(i) = 0
        end do !i
        do i = 1,natom
        startPairL13(i)=0
        startPairL14(i)=0
        startPairL14V(i)=0
        nPairL13(i)=0
        nPairL14(i)=0
        nPairL14V(i)=0
        end do!i
        do i=1,nImp1234MAX
        quarImp1234L(i)=0
        end do
c
c define vbondPair
         do iam = 1,nmoveatom
         ia = moveAtomList(iam)
c         if(realAtomFlag(ia) .eq. 1 )then
         if(realAtomFlag(ia) .ge. 1 )then
         if(nPairL12(ia) .ge. 1)then
         do p12i = startPairL12(ia),startPairL12(ia)+nPairL12(ia)-1
c take atom ja
         ja = pair12List(p12i)
c define addNew value
          addNew = .false.
          if(ivar .eq. 0 ) then
          if(realAtomFlag(ja) .ge. 1 ) addNew = .true.
          end if
c
          if(ivar .eq. 1 ) then
          if(realAtomFlag(ja) .eq. realAtomFlag(ia) )addNew = .true.
          end if
c
         if(addNew)then  !correct addNew
         addNew = .false.
         if(moveFlag(ja) .eq. 0) then
         addNew = .true.
         else
         if(ja .gt. ia )addNew = .true. 
         end if !
         end if ! correct addNew 
c
        if(CONTROL1)then
        write(kanalp,*)'vbondListAAseq:ia,ja,realFlg(ja),moveFlg(ja):',
     &  ia,ja,realAtomFlag(ja),moveFlag(ja)
        write(kanalp,*)'vbondListAAseq: addNew ja:',addNew
        end if !C
c
         if(addNew)then       !addNew bond12()
cx         if(ja .gt. ia )then
         nb12s = nb12s+1
         nb12s2 = 2*nb12s
c control
         if(nb12s .gt. np12MAX)then
         write(kanalp,*)'ERROR:vbondList:np12MAX is low'
         stop
c
         write(kanalPStat,*)mError,
     &    ' vbondList:np12MAX is low' 
c
         end if!control nb12s
         bond12List(nb12s2-1)=ia
         bond12List(nb12s2)=ja
cx        end if !ja .gt. ia
         end if !ja real
         end do !p12i
         end if ! nPairL12(ia) .ge. 1
         end if ! addNew 
         end do ! iam
         nbond12 = nb12s
c
c pair12 is counted twice: ia-ja and ja-ia
c bond12 is counted once
c
        if(CONTROL1)then
        write(kanalp,*)'vbondListAAseq:INput: startPairL12:'
        call print10(natom,startPairL12)
        write(kanalp,*)'vbondListAAseq:INput:nPairL12:'
        call print10(natom,nPairL12)
        write(kanalp,*)'vbondListAAseq:INput:pair12List:np12s:',np12s
        call print10(np12s,pair12List)
c
        write(kanalp,*)'vbondListAAseq:OUt:bond12List12:nbond12:',
     &                  nbond12
        call print10(nb12s2,bond12List)
        end if !Control
c
c for Shake: the sorted  bond12List
        nbond12noH = 0
        do i12 = 1,nbond12
        i2 = i12*2-1
c copy temporary to trip123List
         trip123List(i2) = bond12List(i2)
         trip123List(i2+1) = bond12List(i2+1)
         end do !i12
c
         n12H = 0
         do i12 = 1,nbond12
         i2 = i12*2 - 1
         ia = trip123List(i2)
         ja = trip123List(i2+1)
         if(atomName(ia)(1:1) .ne. 'H' .and.
     &     atomName(ja)(1:1) .ne. 'H') then
         nbond12noH = nbond12noH + 1 
         i2nH = 2*nbond12noH - 1 
         bond12List(i2nH) = ia
         bond12List(i2nH+1) = ja         
         else
c put H (containing) pairs on the end of List
         n12H = n12H +1
         i2H = (nbond12+1-n12H)*2 -1
         bond12List(i2H) = trip123List(i2)
         bond12List(i2H+1) = trip123List(i2+1)
         end if
c 12 pairs to be Shaked are nbond12noH+1 -> nbond12 (iShake=1, H-HeavyAt)
c                                      1 -> nbond12 (iShake=2, all bonds)
c
         trip123List(i2)=0
         trip123List(i2+1)=0  
c
         end do !i12
c end of sorting  bond12List()
        if(CONTROL1)then
        write(kanalp,*)'vbondList: sorted:bond12List12:nbond12nH:',
     &                  nbond12noH
        call print10(nbond12,bond12List)
        end if !Control
c        
c pair13List,startPairL13,nPairL13
c
        do iam=1,nmoveatom
        ia = moveAtomList(iam)
cx        if(realAtomFlag(ia) .eq. 1 )then !realAtom
        if(realAtomFlag(ia) .ge. 1 )then !realAtom
c
        p12is = startPairL12(ia)
        if ( nPairL12(ia) .eq. 0 )then       
         if(CONTROL1)then
         write(kanalp,*)' zero 12pairList: ia:',ia
         end if!C
        else                            !nonZerop12List
c
        p12if = p12is+nPairL12(ia)-1

c        if(CONTROL1)then
c        write(kanalp,*)'p13: ia:',ia,'p12is,p12if:',p12is,p12if
c        end if        

c counters for ia block 
        np13=0
        np14=0
        np123=0
        np1234=0
c loop over ja atoms
        do p12i = p12is,p12if
c take atom ja
        ja = pair12List(p12i)
cx        if(realAtomFlag(ja) .eq. 1 )then !realAtom
        if(realAtomFlag(ja) .ge. 1 )then !realAtom
c
c       if(CONTROL1)then
c       write(kanalp,*)'p13: ja:',ja,' p12i:',p12i
c       end if        
c
c loop over ka atoms bonded to ja
        p12js = startPairL12(ja)
        if ( nPairL12(ja) .gt. 0 )then       !nonzero List
        p12jf = p12js+nPairL12(ja)-1

c       if(CONTROL1)then
c       write(kanalp,*)'p13: ja:',ja,'p12js,p12jf:',p12js,p12jf
c       end if        
c
        do p12j = p12js,p12jf
        
        ka =  pair12List(p12j)  ! triplet  1-2-3=ia-ja-ka

        if(realAtomFlag(ka) .ge. 1 )then !realAtom

c       if(CONTROL1)then
c       write(kanalp,*)'p13: ka:',ka,' p12j:',p12j
c       end if        

        if( ka .ne. ia)then
        np13 = np13+1
        np13s = np13s+1
        if(np13s .gt. np13MAX)then
        write(kanalp,*)'ERROR:vbondList:np13MAX is low'
c
         write(kanalPStat,*)mError,
     &  ' vbondList:np13MAX is low ..'
c
        stop
        end if !control np13s       
c symmetrical 13 pairList
        pair13List(np13s) =  ka
c 
c collect triplets
c define  addNew value
          addNew = .false.
          if(ivar .eq. 0 )then
          if(moveFlag(ka) .eq. 1 .and. ka .gt. ia ) addNew = .true. 
          if(moveFlag(ka) .eq. 0 ) addNew = .true. 
          end if !ivar 0
c
          if(ivar .eq. 1 )then
          if(realAtomFlag(ja) .eq. realAtomFlag(ia) .and. 
     &       realAtomFlag(ka) .eq. realAtomFlag(ia) )then
          if(moveFlag(ka) .eq. 1 .and. ka .gt. ia ) addNew = .true.
          if(moveFlag(ka) .eq. 0 ) addNew = .true.
          end if
          end if !ivar 1
c
        if(addNew)then    !addNew triplet
        np123 = np123 + 1
        np123s = np123s + 1
        np123s2=np123s*3
c control
        if(CONTROL1)then
        write(kanalp,*)'ia,ja,ka:',ia,ja,ka
        end if
c control        
        if(np123s2 .ge. np123MAX)then
        write(kanalp,*)'ERROR:vbondList: low:np123MAX:',np123MAX
c
            write(kanalPStat,*)mError,
     &  ' vbondList:np123MAX is low ..'
c
        stop
        end if !control np13s
c
        trip123List(np123s2-2) =  ia
        trip123List(np123s2-1) =  ja
        trip123List(np123s2) =    ka
cx	end if !ka .gt. ia !collect triplet
        end if !addNew triplet
c
c 14 pair
c   loop over m atoms bonded to ka
        p12ks = startPairL12(ka)
        if ( nPairL12(ka) .gt. 0)then       !nonzero List
        p12kf = p12ks+nPairL12(ka) - 1
c        
        do p12k = p12ks,p12kf
        ma = pair12List(p12k)

cx        if(realAtomFlag(ma) .eq. 1 )then !realAtom ma
        if(realAtomFlag(ma) .ge. 1 )then !realAtom ma
c
        if( ma .ne. ja )then 
        np14=np14 + 1
        np14s = np14s + 1
c control
        if(np14s .gt. np14MAX)then
        write(kanalp,*)'ERROR:vbondList:np14MAX is low:np14MAX='
     &  ,np14MAX, 'Edit pair1234array.h :npL14MAX: '
c
        write(kanalPStat,*)mError,
     &  ' vbondList:np14MAX is low !'
c
        stop
        end if !control np14s
c
c collect quartet : ia-ja-ka-ma
c define addNew 
        addNew = .false.
        if(ivar .eq. 0 )then
        if(moveFlag(ma) .eq. 1 .and. ma .gt. ia )addNew = .true. 
        if(moveFlag(ma) .eq. 0 )addNew = .true. 
        end if !ivar 0
        
        if(ivar .eq. 1 .and. 
     &    realAtomFlag(ja) .eq. realAtomFlag(ia) .and. 
     &    realAtomFlag(ka) .eq. realAtomFlag(ia) .and. 
     &    realAtomFlag(ma) .eq. realAtomFlag(ia)) then
        if(moveFlag(ma) .eq. 1 .and. ma .gt. ia )addNew = .true. 
        if(moveFlag(ma) .eq. 0 )addNew = .true. 
        end if !ivar 1  

        if(addNew) then   !addNew quartet        
        np1234=np1234 + 1
        np1234s=np1234s + 1
        np1234s3=4*np1234s
c control
        if(np1234s3 .ge. np1234MAX)then
        write(kanalp,*)'ERROR:vbondList:np1234MAX is low'
c
        write(kanalPStat,*)mError,
     &    ' vbondList:np1234MAX is low '
c
        stop
        end if !control np14s
c control
        if(CONTROL1)then
        write(kanalp,*)'ia,ja,ka,ma:',ia,ja,ka,ma
        end if
c
        quar1234List(np1234s3-3)=ia
        quar1234List(np1234s3-2)=ja
        quar1234List(np1234s3-1)=ka
        quar1234List(np1234s3)  =ma
c
cx        end if !ma .gt. ja !collect quartet
          end if ! addNew quartet
c
c symmetrical 14 pairList
        pair14List(np14s) = ma 
        end if ! ma .ne. ja       
        end if ! realAtom ma
        end do! p12k 
        end if ! nonzero List (ka)
c
        end if ! ka.ne.ia
        end if ! realAtom ka
        end do ! p12j
        end if ! nPairL12(ja)>0
c
        end if ! realAtom ja
        end do !p12i
c
        if(np13 .gt. 0 )then
        nPairL13(ia) = np13               !number of records in the List for ia
        startPairL13(ia) = np13s-np13+1  !position of the first rec in the List
        if(np14 .gt. 0 )then
        nPairL14(ia) = np14               !number of records in the List for ia
        startPairL14(ia) = np14s-np14+1  !position of the first rec in the List
        end if ! np14 >0
        end if ! np13 >0
c 
        end if ! p12si#0 nonZero p12List 
        end if ! realAtom ia
        end do !iam
c
c trip123 : ia-ja-ka , ka>ia       
        nTrip123 = np123s
        nQuar1234 = np1234s
c       
        if(CONTROL1)then
        np123s2=np123s*3
        write(kanalp,*)'trip123List:nTrip123:',np123s
        call print10(np123s2,trip123List)
c pair13 is counted twice i.e. as 3-1 ; ia-ja-ka and ka-ja-ia
        write(kanalp,*)'startPairL13:'
        call print10(natom,startPairL13)
        write(kanalp,*)'nPairL13:'
        call print10(natom,nPairL13)
        write(kanalp,*)'pair13List:np13s:',np13s
        call print10(np12s,pair13List)
c
        write(kanalp,*)'quar1234List:nQuar1234:',np1234s
c quartets are counted once: i-j-k-m
        call print10(np1234s3,quar1234List)
c
c pair14 is counted twice as 1-4 and 4-1; ia-ja-ka-ma and ma-ja-ka-ia
        write(kanalp,*)'startPairL14:'
        call print10(natom,startPairL14)
        write(kanalp,*)'nPairL14:'
        call print10(natom,nPairL14)
        write(kanalp,*)'pair14List:np14s:',np14s
        call print10(np14s,pair14List)
c
        end if !Control
c
c define IMPROPER quaret in plane  j2-ia-j1
c                                     !
c                                     j3
c zero torsion angle for atoms j1-ia-j2-j3
        nImp1234=0
        do iam=1,nmoveatom
        ia = moveAtomList(iam)
c        if(realAtomFlag(ia) .eq. 1 )then !realAtom ia
        if(realAtomFlag(ia) .ge. 1 )then !realAtom ia
c
        tri = 3
        if ( nPairL12(ia) .eq. tri )then   !planar center
        p12is = startPairL12(ia) 
c
c define addNew
        addNew = .false.
        if(ivar .eq. 0 .and.  
     &    (realAtomFlag(pair12List(p12is)) .ge. 1 .and.
     &     realAtomFlag(pair12List(p12is+1)) .ge. 1 .and.
     &     realAtomFlag(pair12List(p12is+2)) .ge. 1 )) addNew = .true. 
c
        if(ivar .eq. 1 .and.  
     &    (realAtomFlag(pair12List(p12is)) .eq. realAtomFlag(ia) .and.
     &     realAtomFlag(pair12List(p12is+1)) .eq. realAtomFlag(ia) .and.
     &     realAtomFlag(pair12List(p12is+2)) .eq. realAtomFlag(ia) ))
     &    addNew = .true.
c        
        if(addNew)then
        nImp1234=nImp1234+1
c
        if(nImp1234*4 .ge. nImp1234MAX)then
        write(kanalp,*)'vbondList2:ERROR!:nImp1234MAX is low',nImp1234*4
c
        write(kanalPStat,*)mError,
     &    ' nImp1234MAX is low '
c
        stop
        end if !
c
cx        p12is = startPairL12(ia)
        i4=4*nImp1234-4
        quarImp1234L(i4+1)=pair12List(p12is)     
        quarImp1234L(i4+2)= ia     
        quarImp1234L(i4+3)=pair12List(p12is+1)     
        quarImp1234L(i4+4)=pair12List(p12is+2)     
c
cx       end if ! realAtoms j,k,l
        end if ! addNew
        end if! ! planar center
        end if !realAtom ia
        end do!iam
c 
        if(CONTROL1)then
        write(kanalp,*)'nImp1234:',nImp1234
        i4=4*nImp1234
         write(kanalp,*)'quarImp1234L:'
        call print10(i4,quarImp1234L)
        end if !C1
c
c define pair14ListV() : list of 14 atoms which interect via VDW/COUL
c                      : remove at4 from list if at4 is equal to 2 or 3 neigb
c                      : for another path, i.e. 3-,4-,5- atomic rings
        if(CONTROL1)then 
        write(kanalp,*)'vbondListSB: start pair14ListV'
        end if
c
         do iam = 1,nmoveatom
         ia = moveAtomList(iam)
         if(realAtomFlag(ia) .ge. 1 )then
c
         start_ia = .false.                   ! flag defined that 14V neigh for ia exists
         do p14i=startPairL14(ia),startPairL14(ia)+nPairL14(ia)-1        
         ja = pair14List(p14i) 
c
c BlockUnit
c
cx        write(kanalp,*)'vbondListAA: ia,ja,resia,resja,blia,blja:',
cx   &  ia,ja,resNumb(ia),resNumb(ja),atomBlockName(ia),
cx   &  atomBlockName(ja)
c
         if(OPT_removePairInHCycBase)then
        if((resNumb(ia) .eq. resNumb(ja)) .and.
     &      (atomBlockName(ia) .eq. atomBlockName(ja)))then
c
cx        write(kanalp,*)'vbondListAA: ia,ja,resia,resja,blia,blja:',
cx     &  ia,ja,resNumb(ia),resNumb(ja),atomBlockName(ia),
cx     &  atomBlockName(ja),
cx     &  ' skip Ja'
c
          goto 1401            ! skip atom ja
         end if ! atomBlockName
        end if !OPT_removePairInHCycBase 
c
         call compAtomPairList(ja,ia,
     &             pair12List,startPairL12,npairL12,find)
         if(find)goto 1401   !skip ja
c
	 call compAtomPairList(ja,ia,
     &             pair13List,startPairL13,npairL13,find)
c
        if(find)goto 1401   ! skip ja
c
        if(.not. start_ia)then
        start_ia = .true.
        startPairL14V(ia) = np14Vs+1
        end if! start_ia
c
	np14Vs = np14Vs+1
        pair14VList(np14Vs) = ja 
        npairL14V(ia) = npairL14V(ia) + 1
c
1401    continue
        end do !p14i
c
         end if! realAtomFlag(ia) .ge. 1
         end do !iam
c
        if(CONTROL1)then 
        write(kanalp,*)'startPairL14:'
        call print10(natom,startPairL14)
        write(kanalp,*)'nPairL14:'
        call print10(natom,nPairL14)
        write(kanalp,*)'pair14List:np14s:',np14s
        call print10(np14s,pair14List)
c
       write(kanalp,*)'startPairL14V:'
        call print10(natom,startPairL14V)
        write(kanalp,*)'nPairL14V:'
        call print10(natom,nPairL14V)
        write(kanalp,*)'pair14VList:np14Vs:',np14Vs
        call print10(np14Vs,pair14VList)
        end if!CONTROL
c
        if(CONTROL)then 
        write(kanalp,*)'end of vbondList:' 
        end if
c
            return
            end
