c getSSbonds
c
c YN Vorobjev, 2005
c
	subroutine getSSbonds
c
c calculates: nSSbonds,ssBondAt12List(2*nSS)[iS,jS]:
c             ssBondAtFlag(ia) = 0/1
c             ssBondDist0(nSS)
c
        include 'xyzPDBsize.h'
c
        include "output.h"
        include "kanalUse.h"
        include "statusMessg_mDyn.h"
c
        include 'xyzPDBinfo.h'
        include 'xyzPDBcrd.h'
        include 'pair1234array.h' 
        include 'ssbond.h'
c
        integer i,ia,ia3,k
        integer j,ja,ja3
        integer ic,ic2   
        real dss2
        logical CONTROL
        integer kanalp
c
        kanalp = kanalRunOut
        CONTROL = .true.
c
c c calculate CYSresList()
c CYS list is needed to find S-S bonds
c CYX - amber name for CYS involved into S-S
c ssCYSname = "CYX "  ssAtname="SG  "
c
        nCYSres = 0
        nSSbondAtom = 0
        nSSbonds = 0
        if(OPT_SSBonds)then
        do i = 1,nres
        if(resNameRes(i)(1:4) .eq. ssCYSname)then
        nCYSres = nCYSres + 1
        cysResList(nCYSres) = i
c
        do ia = startAtInRes(i),stopAtInRes(i)
        if (atomName(ia)(1:4) .eq. ssAtname )then
        nSSbondAtom = nSSbondAtom +1
        ssAtomList(nSSbondAtom) = ia
        ssBondAtFlag(nSSbondAtom) = 0
        end if !.eq. ssAtname
        end do !ia
c
        if(nCYSres .ge. nCYSresMAX)then
        write(kanalp,*)'ERROR!:initMolecTopSeq01_SB nCYSresMAX is low',
     &  ' increase nCYSresMAX in ssbond.h and recompile program !!'
        write(kanalPStat,*)mError,
     &    ' initMolecTopSeq01 nCYSresMAX is low !'
        stop
        end if !.ge. nCYSresMAX
c
        end if !.eq. ssCYSname
        end do !i:nres
c
        if(CONTROL)then
        write(kanalp,*)'SS residue list:'
        do i=1,nCYSres
        write(kanalp,'(a40,i4,1x,i4,1x,a4)')
     &  'initMolecTopSeq01: cysResList:',
     &  i,cysResList(i),resNameRes(cysResList(i))
        end do!i
c
        write(kanalp,*)'SS atoms S-S bond candidates:'
        do i=1,nSSbondAtom
        write(kanalp,'(a40,i4,1x,i4,1x,a4,1x,a4,i4,i4)')
     &  'initMolecTopSeq01: ssAtomList:',
     &  i,ssAtomList(i),atomName(ssAtomList(i)),
     &  resName(ssAtomList(i)),resNumb(ssAtomList(i)),
     &  ssBondAtFlag(i)
        end do!i
        end if!C
c
c define S-S bonds atom pairs
c
         if(nSSbondAtom .gt. 1)then
         do i=1,nSSbondAtom-1
         ia=ssAtomList(i)
c
cx         write(kanalp,*)'initMolecTopSeq01:SS:i,ia=',i,ia
c
         if(ssBondAtFlag(i) .eq. 0)then
         ia3=3*ia-3
         do j=i+1,nSSbondAtom
         if(ssBondAtFlag(j) .eq. 0)then
         ja=ssAtomList(j)
c
cx         write(kanalp,*)'initMolecTopSeq01:SS:j,ja=',j,ja
c
         ja3=3*ja-3
          if( ia .ne. ja)then
          dss2 = (atomXYZ(ia3+1)-atomXYZ(ja3+1))**2 +
     &               (atomXYZ(ia3+2)-atomXYZ(ja3+2))**2 +
     &               (atomXYZ(ia3+3)-atomXYZ(ja3+3))**2
c
cx         write(kanalp,*)'ia:XYZ:',ia,
cx     &   atomXYZ(ia3+1),atomXYZ(ia3+2),atomXYZ(ia3+3)
cx         write(kanalp,*)'ja:XYZ:',ja,
cx     &   atomXYZ(ja3+1),atomXYZ(ja3+2),atomXYZ(ja3+3)
cx         write(kanalp,*)'dist2: ia-ja:',dss2,' dss2MAX_OPT:',dss2MAX_OPT
c
          if(dss2 .lt. dss2MAX_OPT)then
          nSSbonds = nSSbonds +1
          ssBondAt12List(2*nSSbonds-1)=ia
          ssBondAt12List(2*nSSbonds)= ja
          ssBondAtFlag(i)=1
          ssBondAtFlag(j)=1
          ssBondDist0(nSSbonds)=sqrt(dss2)
c
cx          write(kanalp,*)'nSSbonds,ia,ja:',nSSbonds,ia,ja
c
          end if !dss2.lt.
          end if ! ia .ne. ja
          end if ! jFlag
          end do!j
          end if ! iFlag
          end do!i
          end if ! nSSbondAtom .gt. 1
c
          if(CONTROL)then
          write(kanalp,*)
     &    'initMolecTopSeq01: SS bonds list:nSSbonds=',nSSbonds
          if(nSSbonds .gt. 0)then
          write(kanalp,*)
     &    'nSS  ia1  Name  Res  ia2  Name  Res ssBondDist0'
          do i=1,nSSbonds
          ia=ssBondAt12List(2*i-1)
          ja=ssBondAt12List(2*i)
          write(kanalp,'(i3,1x,2(i5,1x,a4,1x,a4,1x,i4),f8.3)')
     &                i,ia,atomName(ia),resName(ia),resNumb(ia),
     &                ja,atomName(ja),resName(ja),resNumb(ja),
     &                ssBondDist0(i)
          end do!i
          end if !nSSbonds .gt. 0
          end if !C
c
c correct Pair12ListREseq() and nPairL12REseq(ia)
c insert SS bonds into nPairL12(*),startPairL12(*),pair12List(*) data structure
c
          do ic = 1,nSSbonds
          ic2 = 2*ic-1
          ia=ssBondAt12List(ic2)
          ja=ssBondAt12List(ic2+1) 
          call insert12bond
     &          (ia,ja,natom,nPairL12,startPairL12,pair12List)
          end do !ic
c 
         end if !OPT_SSBonds
c
	 return
	 end
