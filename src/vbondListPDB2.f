c defines valence bonded 12, 123, 1234 neighbours
c Y Vorobjev 2002
c calculate connectivity from at-at distances and atomic COValent rad
c 
        subroutine vbondListPDB2(atomXYZ,
     &           natom,atomNumb,atomName,resName,chName,resNumb,
     &           nres,resNameRes,chNameRes,
     &           atomNameEx,startAtInRes,
     &           nmoveatom,moveAtomList, 
     &           pair12List,startPairL12,nPairL12,np12MAX,
     &           pair13List,startPairL13,nPairL13,np13MAX,
     &           pair14List,startPairL14,nPairL14,np14MAX,
     &           bond12List,nbond12,
     &           trip123List,nTrip123,np123MAX,
     &           quar1234List,nQuar1234,np1234MAX,
     &           quarImp1234L,nImp1234,nImp1234MAX)
c
c pair12List(i2) - list of atoms 2 of 12 pair
c startPairL12(i1) - poit to first position in pair12List() for 
c                    the sequence of atoms2 for atom1=i1
c nPairL12(i1) - length of sequence for atoms2 
c pair13List() - list of atoms 3 of 13 pair
c pair14List() - list of atoms 4 of 14 pair 
c trip123List() - list of atoms 2,3 in 123 triplets
c quar1234List() - list of atoms 2,3,4 in  1234 quarets
c quarImp1234L() - list of atoms 1,2,3,4 in 1234 improper quartets
c 
        implicit none
c
        include "charStringSiz.h"
        include "output.h"
        include "kanalUse.h"
        include "statusMessg_mDyn.h"
c
        real atomXYZ(*)
        integer natom
        character*4 atomName(*)
        character*4 resName(*)
        character*1 chName(*)
        integer atomNumb(*),resnumb(*)
        integer nres
        character*4 resNameRes(*)
        character*1 chNameRes(*) 
        character*8 atomNameEx(*)
        integer startAtInRes(*)
        integer nmoveatom,moveAtomList(*)
cRESULT:
        integer np12MAX,np13MAX,np14MAX
        integer np123MAX,np1234MAX
        integer pair12List(*),startPairL12(*),nPairL12(*)
        integer pair13List(*),startPairL13(*),nPairL13(*)
        integer pair14List(*),startPairL14(*),nPairL14(*)
        integer bond12List(*),nbond12
        integer trip123List(*),nTrip123
        integer quar1234List(*),nQuar1234
        integer nImp1234,nImp1234MAX
        integer quarImp1234L(*)
c local CYS numbers
        integer nCYSresMAX
        parameter (nCYSresMAX = 51)
        integer nCYSres
        integer cysResList(nCYSresMAX) !defined here
        integer nRESneigh
        integer resNeighList(nCYSresMAX+3)
c local variables
        integer i,j,k,k3,i2,i4
        integer i3,j3,ir,jr,jrn
        integer ia,ja,ka,ma,iam
        integer sjr,fjr
        integer np12,np12s
        integer nb12s,nb12s2
        integer np13,np13s
        integer np14,np14s
        integer np123,np1234
        integer np123s,np1234s
        integer np123s2,np1234s3
        integer p12i,p12is,p12if
        integer p12j,p12js,p12jf
        integer p12k,p12ks,p12kf
        real radi,radj,dij2 
        real xyzi(3)
        integer kanalp,tri
        integer case12
        logical doPair
        logical vBondij
        logical CONTROL,CONTROL1
c
        kanalp=kanalRunOut
        CONTROL = .true.   !control print
        CONTROL1= .false.
        case12 = 0    !rule to define 12 neighbours. 
c                     !0 - symmetrical all pairs, 1 - nonsymmetrical i1 < j2
        doPAir = .false.
        vBondij = .false.
        
        if(CONTROL)then 
        write(kanalp,*)'Start vbondList2SB:'
        end if
c
c calculate connectivity from at-at distances and atomic COValent rad
c initialize
        np12s = 0
        np13s = 0
        np14s = 0
        nb12s = 0
        nb12s2 = 0
        np123s = 0
        np1234s = 0
        np123s2 = 0 
        np1234s3 = 0
c
        do i=1,np12MAX
        pair12List(i) = 0
        bond12List(i)=0
        end do !i
        do i=1,np13MAX
        pair13List(i) = 0
        end do !i
        do i=1,np14MAX
        pair14List(i) = 0
        end do !i
        do i=1,np123MAX
        trip123List(i) = 0
        end do !i
        do i=1,np1234MAX
        quar1234List(i) = 0
        end do !i
        do i = 1,natom
        startPairL12(i)=0
        startPairL13(i)=0
        startPairL14(i)=0
        nPairL12(i)=0
        nPairL13(i)=0
        nPairL14(i)=0
        end do!i
        do i=1,nImp1234MAX
        quarImp1234L(i)=0
        end do
c calculate CYSresList()
c CYS list is needed to find S-S bonds
c CYX - amber name for CYS involved into S-S
        nCYSres = 0
        do i = 1,nres
        if(resNameRes(i)(1:2) .eq. 'CY')then
        nCYSres = nCYSres + 1
        cysResList(nCYSres) = i
        if(nCYSres .ge. nCYSresMAX)then
        write(kanalp,*)'ERROR!:vbondListPDB2:nCYSresMAX is low'
c
        write(kanalPStat,*)mError,
     &    ' vbondListPDB2:nCYSresMAX is low !'
c
        stop
        end if
        end if
        end do ! i 
c
c calculates 12 bonds if 1 is a moving atom
        do iam=1,nmoveatom
        ia = moveAtomList(iam)
        ir = resNumb(ia)
c loop over residues and atoms in the res
c define range for jres
         if(case12 .eq. 0)then     ! ja>#<ia
         sjr = ir - 1             
        if(sjr .lt. 1)sjr=1
        if(ir .gt. 1 .and. 
     &      chNameRes(ir) .ne. chNameRes(ir-1))sjr=ir
        else          
        sjr = ir                  !case ja>ia 
        end if !case12

        fjr = ir+1
        if(fjr .gt. nres)fjr=nres
        if(chNameRes(ir) .ne. chNameRes(ir+1))fjr=ir

c take atoms in resi 
         np12 = 0                                 !atom counter
        call getCovAtRad(atomName(ia)(1:1),radi)

        i3=3*ia-3
        xyzi(1)=atomXYZ(i3+1)
        xyzi(2)=atomXYZ(i3+2)
        xyzi(3)=atomXYZ(i3+3)
c
c make residue List to loop over 
        nRESneigh = 0
        do jr = sjr,fjr
        nRESneigh = nRESneigh + 1
        resNeighList(nRESneigh) = jr 
        end do !jr
c add CYS list
        if(resName(ia)(1:2) .eq. 'CY' 
     &     .and. atomName(ia)(1:2) .eq. 'SG')then
        do jr = 1,nCYSres
        if(resNumb(ia) .ne. cysResList(jr))then
        nRESneigh = nRESneigh + 1
        resNeighList(nRESneigh) = cysResList(jr)
        end if !exclude ia res
        end do !jr
        end if ! CYS
c
c loop over neihbour residues to find vbond12
        do jrn = 1,nRESneigh
        jr = resNeighList(jrn)

c take atom in resJ
        do ja=startAtInRes(jr),startAtInRes(jr+1)-1
         doPair=.false.
c
        if( chName(ia) .eq. chName(ja) )then  ! vBonds inside chain
        if(case12 .eq. 0 .and. ja .ne. ia) doPair=.true. !count twice
        end if ! chName(ia) .eq. chName(ja)
c
        vBondij = .false.
c       
        if( doPair )then 
        call getCovAtRad(atomName(ja),radj) 
        j3=3*ja-3
        dij2=       (atomXYZ(j3+1)-xyzi(1))**2
        dij2=dij2 + (atomXYZ(j3+2)-xyzi(2))**2
        dij2=dij2 + (atomXYZ(j3+3)-xyzi(3))**2
c
c       if(CONTROL1)then
c       write(kanalp,*)
c    &  'ia:',ia,' ja: ',ja,' dij2:',dij2,'radi,radj:',radi,radj
c       end if
c
c general bond
        if(resNumb(ia) .eq. resNumb(ja) ) then
        if(dij2 .le. (radi+radj)**2) vBondij = .true.
        else
c
c make peptide C-N bond
        if( (resNumb(ja) - resNumb(ia)) .eq. 1 ) then
        if( (atomName(ia)(1:3) .eq. 'C  ' .and. 
     &       atomName(ja)(1:3) .eq. 'N  ') ) vBondij = .true.
c
        if(vBondij .and. CONTROL)then
        write(kanalp,*)
     & 'C-N vBondij: ia:atNam,resNumb:',ia,atomName(ia),resNumb(ia), 
     & ' ja:atNam,resNumb:',ja,atomName(ja),resNumb(ja),
     & ' dij2:',dij2,' np12:',np12
       end if ! vBondij .and. CONTROL
c
       end if ! resNumb(ja) - resNumb(ia)) .eq. 1
       end if ! resNumb(ia) .eq. resNumb(ja) 
c
c  add S-S bond
        if( resName(ia)(1:2) .eq. 'CY' 
     &     .and. atomName(ia)(1:2) .eq. 'SG'
     &     .and. atomName(ja)(1:2) .eq. 'CY' 
     &     .and. atomName(ja)(1:2) .eq. 'SG' )then
        if(dij2 .le. (radi+radj)**2) vBondij = .true. 
        end if
c
        if(vBondij)then
        np12 = np12+1
        np12s = np12s+1
        
       if(CONTROL1)then
       write(kanalp,*)
     & 'vBondij: ia:atNam,resNumb:',ia,atomName(ia),resNumb(ia), 
     & ' ja:atNam,resNumb:',ja,atomName(ja),resNumb(ja),
     & ' dij2:',dij2,' np12:',np12
       end if
c control        
        if(np12s .gt. np12MAX)then
        write(kanalp,*)'ERROR:vbondList:np12MAX is low'
c
         write(kanalPStat,*)mError,
     &    ' vbondList:np12MAX is low'
c
        stop
        end if!control np12s
c
        pair12List(np12s) = ja
c
c collect bonds12
        if(ja .gt. ia)then
        nb12s = nb12s+1
        nb12s2 = 2*nb12s
c control
        if(nb12s .gt. np12MAX)then
        write(kanalp,*)'ERROR:vbondList:np12MAX is low'
c
         write(kanalPStat,*)mError,
     &    ' vbondList:np12MAX is low'
c
        stop
        end if!control np12s
        bond12List(nb12s2-1)=ia
        bond12List(nb12s2)=ja
        end if !ja .gt. ia
c
        end if !vBondij
        end if !doPair 
        end do !ja
        end do !jrn
c
        if(np12 .gt. 0 )then
        nPairL12(ia) = np12             !number of records in the List for ia
        startPairL12(ia) = np12s-np12+1 !position of the first rec in the List
        end if !np12>0
c        
        end do !iam
        nbond12 = nb12s
c
c pair12 is counted twice: ia-ja and ja-ia
c bond12 is counted once
        if(CONTROL)then
        write(kanalp,*)'startPairL12:'
        call print10(natom,startPairL12)
        write(kanalp,*)'nPairL12:'
        call print10(natom,nPairL12)
        write(kanalp,*)'pair12List:np12s:',np12s
        call print10(np12s,pair12List)
c
        write(kanalp,*)'bond12List12:nbond12:',nbond12
        call print10(nb12s2,bond12List)
        end if !Control
c
c pair13List,startPairL13,nPairL13
c
        do iam=1,nmoveatom
        ia = moveAtomList(iam)
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

c       if(CONTROL1)then
c       write(kanalp,*)'p13: ja:',ja,' p12i:',p12i
c       end if        

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
     &    ' vbondList:np13MAX is low' 
c
        stop
        end if !control np13s       
c symmetrical 13 pairList
        pair13List(np13s) =  ka
c 
c collect triplets
        if( ka .gt. ia)then !uporjadochennii spisok, sortedList of triplets
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
     &    ' vbondList: low:np123MAX:'
c
        stop
        end if !control np13s
c
        trip123List(np123s2-2) =  ia
        trip123List(np123s2-1) =  ja
        trip123List(np123s2) =    ka
	end if !ka .gt. ia !collect triplet
c
c 14 pair
c   loop over m atoms bonded to ka
        p12ks = startPairL12(ka)
        if ( nPairL12(ka) .gt. 0)then       !nonzero List
        p12kf = p12ks+nPairL12(ka) - 1
c        
        do p12k = p12ks,p12kf
        ma = pair12List(p12k)
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
     &    ' np14MAX is low !'
c
        stop
        end if !control np14s
c
c collect quartet : ia-ja-ka-ma
        if( ma .gt. ia )then  !uporjadochenni spisok, sortedList
        np1234=np1234 + 1
        np1234s=np1234s + 1
        np1234s3=4*np1234s
c control
        if(np1234s3 .ge. np1234MAX)then
        write(kanalp,*)'ERROR:vbondList:np1234MAX is low'
c
        write(kanalPStat,*)mError,
     &    ' np1234MAX is low'
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
        end if !ma .gt. ja !collect quaret
c
c symmetrical 14 pairList
        pair14List(np14s) = ma 
        end if ! ma .ne. ja       
        end do! p12k 
        end if ! nonzero List (ka)
c
        end if ! ka.ne.ia
        end do ! p12j
        end if ! nPairL12(ja)>0
c
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
        end do !iam
c trip123 : ia-ja-ka , ka>ia       
        nTrip123 = np123s
        nQuar1234 = np1234s
c       
        if(CONTROL)then
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
c
        tri = 3
        if ( nPairL12(ia) .eq. tri )then   !planar center
        nImp1234=nImp1234+1
        if(nImp1234*4 .ge. nImp1234MAX)then
        write(kanalp,*)'vbondList2:ERROR!:nImp1234MAX is low',nImp1234*4
c
        write(kanalPStat,*)mError,
     &    ' nImp1234MAX is low'
c
        stop
        end if !
        p12is = startPairL12(ia)
        i4=4*nImp1234-4
        quarImp1234L(i4+1)=pair12List(p12is)     
        quarImp1234L(i4+2)= ia     
        quarImp1234L(i4+3)=pair12List(p12is+1)     
        quarImp1234L(i4+4)=pair12List(p12is+2)     
        end if! ! planar center
        end do!iam
c 
        if(CONTROL)then
        write(kanalp,*)'nImp1234:',nImp1234
        i4=4*nImp1234
         write(kanalp,*)'quarImp1234L:'
        call print10(i4,quarImp1234L)
        end if
c
        if(CONTROL)then 
        write(kanalp,*)'end of vbondListSB:' 
        end if
c
7071    format(a6,1x,i4,2x,a4,a4,a1,i4,4x,3f8.3) !orig PDB
c
            return
            end
c
c calculate atomic radii
	subroutine  getCovAtRad(atomName,rad)

        implicit none
        character*1 atomName
        real rad
c
        integer natmTypMX
        parameter (natmTypMX = 6)
        real radV(natmTypMX)
        character*1 radN(natmTypMX)
        data radV /0.9, 0.8, 0.8, 0.40, 1.2, 1.2/
        data radN /'C','N','O','H','S','X'/
        integer i
        logical find
c
        find = .false.
        do i = 1,natmTypMX
        if(atomName .eq. radN(i))then
        rad = radV(i)
        find = .true.
        goto 101     
        end if
        end do !i
        if(.not.find)rad = radV(natmTypMX)
c  
101     return
        end
