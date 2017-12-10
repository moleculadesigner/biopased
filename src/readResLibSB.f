c c Yury Vorobjev 2002
c
c readResMTLibrary : unitResLib
c find the given residue by name resName6
c assemble molecular ZMatrix by adding the ZMatrix for the Residue
c to the existing ZMatrix
c
c Single AA Chain version
c
	subroutine readResLibSB(unitResLib,resName6,
     &                  dumFlag,natom, natinres,
     &                  molZmatrx1,molZmatrx2,molZmatrx3,
     &                  atomQ,ffatomName,
     &               natPairCycMX,natPairCyc,atPairCycList)
c
c IN:
c unitResLib=unit file
c resName6 (char*6)
c dumFlag - 0/1 NO/Yes to remove DUM atoms of the Res from final zMatrix
c
c IN/OUT:
c natom - number of atoms before EXecution, and is a total n atoms 
c         natom = natom + natinres   at the final RETURN  
c natinres - number of atom in this residue
c molZmatrx1():atname,atChemName,BackSidename :3*natom
c molZmatrx2():4,3,2,1 - atoms 4-3-2-1  :4*natom
c molZmatrx3():d34,ang234,fi1234    :3*natom
c startMolZmatrx(iMol) = first position in molZmatrx1,2,3
c atomQ(ia)   - charge on atom
c natPairCyc - number of Cycle in MolecTop
c atPairCycList(2*natPairCyc)=(i1,j1;i2,j2;...) - list of atNumb(global)
c                             wich make vbond to close Cycle
c
       implicit none
       integer unitResLib
       integer natinres,natom
       character*6 resName6
       character*4 molZmatrx1(*)
       integer     molZmatrx2(*)
       real        molZmatrx3(*)
       real        atomQ(*)
       character*2 ffatomName(*)
       integer dumFlag
       integer nAtPairCyc,nAtPairCycMX
       integer atPairCycList(*)
c
        include "charStringSiz.h"
        include "output.h"
        include "kanalUse.h"
        include "statusMessg_mDyn.h"
c
c local
       integer nl,nin
       integer natresMXloc
       character*80 linein
       character*7 rType   !INT  AA, internal AminoAcid
       integer ns,nsN,nsCA, nsC,nsT
       integer naM,naCA,naC,naO,naN
       integer naH,natl
       integer nCycLine
       integer i,i2,i3,i4,i5,i6,k
       integer n1,nrimp,ni
       real qres
       logical OPT_removeDUMM
       integer kanalp
       logical CONTROL2,CONTROL1,CONTROL0
c
       CONTROL2 = .false.
       CONTROL1 = .false.
       CONTROL0 = .true. 
       kanalp = kanalRunOut
c
       OPT_removeDUMM = .false. 
       if(dumFlag .eq. 1) OPT_removeDUMM = .true.
c
       if(CONTROL0)then
       write(kanalp,*)'start readResLib:ResName:',resName6(1:4)
       end if
c
       rewind unitResLib
c init
       natresMXloc=50  ! max atoms per RES 
       nl = 0
       nin = 0
       ns = 0        !natm to remove from resid  when assemble Zmatrx
       nsN = 0       !nN(i-1)    atom Number for aaBack
       nsCA = 0      !nCA(i-1)
       nsC = 0       !nC(i-1)
       nsT = natom   !initial natom at the call of the SUBRoutine
       natinres = 0
       naM = 0       !nN(i+1)
       naCA = 0      !nCA(i)
       naC = 0       !nC(i)
       naO = 0       !nO(i)
       naN = 0       !nN(i)
       naH = 0       !nH(i)
c
c read ResLib file
100   nl = nl + 1
      read(unitResLib,'(a80)', end = 999 )linein
c
      if(CONTROL2)then
      write(kanalp,*)'readResLib:W1:nl:',nl
      write(kanalp,*)linein
      end if
c
      if(linein(1:1) .eq. '!' .or. linein(1:1) .eq. '#')goto 100 
      nin = nin + 1  ! number of informative line
c
      if(linein(1:6) .eq. "$MTRES" ) then  !start $MTRES
c start ResBlock 
101   nl = nl + 1
      read(unitResLib,'(a80)', end = 999 )linein
c
      if(linein(1:1) .eq. '!' .or. linein(1:1) .eq. '#')goto 101 
      nin = nin + 1  ! number of informative line
c
      if(linein(1:4) .eq. resName6(1:4) )then    ! startMyRes
c read RESidue lines
c
      if(CONTROL1)then
      write(kanalp,*)'W2:MyREs',linein
      end if
c
      read(linein, '(5x,a7,1x,i5,f8.4)')rType,natl,qres
c
      if(CONTROL0)then
      write(kanalp,*)'readResLib: rType,natl,qres:',rType,natl,qres
      end if
c control
      if(natl .gt. natresMXloc)then
      write(kanalp,*)'readResLib:ERROR! natresMXloc is low!!'
      stop
      end if

c define SHIFTs for atom numbering 
c      nsT = natom
c 
c depends on rType(residueType): should be defined for all rType's
cx      if( rType .eq. 'INT  AA') then  !if AimoAcid
        if( rType(6:7) .eq. 'AA') then  !if AimoAcid
c
       if(OPT_removeDUMM)then
       ns = 3      !remove DUMM atoms
       else
       ns = 0      !collect all(DUMM) atoms
       end if !OPT_removeDUMM
c
       naN = 0       !nN(i)
       naH = 0       !nH(i)
c
c find Shifts
c find numbers in preceeding residue: nN,nCA,nC,ntot
      if(natom .gt. 0)then
      if(natom .gt. natresMXloc)then
      n1 = natom-natresMXloc+1
      else
      n1 = 1
      end if !natom .gt. natresMXloc
c
      do i=n1,natom
      i3 = 3*i-3
      if (molZmatrx1(i3+1) .eq. 'N   ')nsN=i
      if (molZmatrx1(i3+1) .eq. 'CA  ')nsCA=i
      if (molZmatrx1(i3+1) .eq. 'C   ')nsC=i
      end do!i
c
      if(CONTROL1)then
      write(kanalp,*)'Shift: nsN, nsCA, nsC:',nsN, nsCA, nsC
      end if !C
      end if! natom .gt.0
c END find Shifts
c
       naN = natom+1       !nN(i)
       naH = naN+1         !nH(i)
       naCA = naN + 2      !nCA(i)
c
       if(CONTROL1)then
       write(kanalp,*)'nsT,ns,naN,naH,naCA :',
     & nsT,ns,naN,naH,naCA
       end if !C
c
      end if ! INT AA   residue
c
c NucleicAcid:
      if( rType .eq. 'INT  NA') then  !if NuclAcid 
c NOT DEFINED yet !!
      continue  
c expand for NA
      end if ! NA
c
c  read next natl lines as Zmatrix  
      i = 0
103   read(unitResLib,'(a80)')linein
      if(linein(1:1) .eq. '!' .or. linein(1:1) .eq. '#')goto 103
      i = i + 1    
c
      if (i .gt. ns )then ! readResLibAtom
      natom = natom+1
      natinres = natinres + 1
      i3 = 3*natom-3
      i4 = 4*natom-4
c
      if(CONTROL1)then
      write(kanalp,*)'W4:Zmtx',linein
      end if
c
c assemble molecZmatrix      
      read(linein, '(6x, a4,2x,a4,2x,a4)')
     &   molZmatrx1(i3+1),molZmatrx1(i3+2),
     &   molZmatrx1(i3+3)
      read(linein,'(i4,16x,i4,i4,i4)')
     &  molZmatrx2(i4+1),molZmatrx2(i4+2),
     &  molZmatrx2(i4+3),molZmatrx2(i4+4)
c
        if(CONTROL1)then
        write(kanalp,*)'before Sh:nat=',natom,' zm1,2,3,4:',
     &  molZmatrx2(i4+1),molZmatrx2(i4+2),
     &  molZmatrx2(i4+3),molZmatrx2(i4+4)
        end if

c make a shift in atomNumbering to assembe protein or isolatedMolec
        if(nsT .ge. 1)then !SHift
       molZmatrx2(i4+1)=molZmatrx2(i4+1)+nsT - ns
c
       if(molZmatrx2(i4+2) .eq. 1)then
        molZmatrx2(i4+2) = nsN 
        goto 202
        end if !1

       if(molZmatrx2(i4+2) .eq. 2)then
        molZmatrx2(i4+2) = nsCA  
        goto 202
        end if !2

       if(molZmatrx2(i4+2) .eq. 3)then
        molZmatrx2(i4+2) = nsC
        goto 202
        end if !3

       if(molZmatrx2(i4+2) .ge. 4) then
        molZmatrx2(i4+2)=molZmatrx2(i4+2)+nsT - ns  
        goto 202
        end if !4

202     continue
c
       if(molZmatrx2(i4+3) .eq. 1) then
        molZmatrx2(i4+3) = nsN     
        goto 203
        end if !1

       if(molZmatrx2(i4+3) .eq. 2) then
        molZmatrx2(i4+3) = nsCA      
        goto 203
        end if !2

       if(molZmatrx2(i4+3) .eq. 3) then
        molZmatrx2(i4+3) = nsC
        goto 203
        end if !3

       if(molZmatrx2(i4+3) .ge. 4)
     &  molZmatrx2(i4+3)=molZmatrx2(i4+3)+nsT-ns    
203     continue   
c
        if(molZmatrx2(i4+4) .eq. 1) then
        molZmatrx2(i4+4) = nsN     
        goto 204
        end if

        if(molZmatrx2(i4+4) .eq. 2) then
        molZmatrx2(i4+4) = nsCA      
        goto 204
        end if

        if(molZmatrx2(i4+4) .eq. 3) then
        molZmatrx2(i4+4) = nsC
        goto 204
        end if 

        if(molZmatrx2(i4+4) .ge. 4)
     &  molZmatrx2(i4+4)=molZmatrx2(i4+4)+nsT - ns
204     continue
c     
        if(CONTROL1)then
        write(kanalp,*)'Shift: nsN, nsCA, nsC:',nsN, nsCA, nsC
        end if !C
c
        if(CONTROL1)then
        write(kanalp,*)'after Sh:zm1,2,3,4:',
     &  molZmatrx2(i4+1),molZmatrx2(i4+2),
     &  molZmatrx2(i4+3),molZmatrx2(i4+4)
        end if
c
        else            ! if first RESidue
         if(OPT_removeDUMM)then
          do k=1,4
          molZmatrx2(i4+k) = molZmatrx2(i4+k)-ns
          end do!k
         end if! _removeDUMM
c
      end if ! SHift
c
       read(linein, '(36x,4f9.4)')
     & molZmatrx3(i3+1),molZmatrx3(i3+2),
     & molZmatrx3(i3+3),atomQ(natom)
c
       ffatomName(natom) = molZmatrx1(i3+2)(1:2)
       molZmatrx1(i3+3) = resName6(1:4)
c
      end if ! readResLibAtom
c
      if( i .lt. natl)goto 103
c 
c define atomNumbers depends on rType !!
c should be defined for all rType's in MTresLib !!
      if( rType .eq. 'INT  AA') then  !if AimoAcid
      naM = natom+1  ! atmN for N (of next res) - atom +M of backbone
      naC = natom-1     !nC(i)
      naO = natom       !nO(i)
      end if ! 'INT  AA'
c
c read LOOP lines
c
301   nl = nl + 1
      read(unitResLib,'(a80)', end = 999 )linein 
c
      if(linein(1:1) .eq. '!' .or. linein(1:1) .eq. '#')goto 301
      nin = nin + 1  ! number of informative line
c
      if(linein(1:4) .eq. "LOOP" )then    ! startLOOP 
c read LOOP lines   
c
c      read(linein, '(8x,i5)')nCycLine
      read(linein, '(4x,i4)')nCycLine
c read number of LOOP lines=nCycLine
c
      if(CONTROL0)then
      write(kanalp,*)'readResLib: rType,nCyclLine:',rType,nCycLine
      end if
c
c read next nCycLine lines
       i = 0
303    read(unitResLib,'(a80)')linein
       if(linein(1:1) .eq. '!' .or. linein(1:1) .eq. '#')goto 303
       i = i + 1
c read atPairCycList()
       nAtPairCyc = nAtPairCyc + 1
c
       if(nAtPairCyc .gt. nAtPairCycMX)then
       write(kanalp,*)'ERROR! nAtPairCycMX is LOW:',nAtPairCyc
       stop
       end if
c
       i2 = 2*nAtPairCyc - 2
       read(linein,'(2x,2i4)')
     &      atPairCycList(i2+1),atPairCycList(i2+2)
c correct atNumeration to global (from local in RESidue)
       atPairCycList(i2+1) = atPairCycList(i2+1) + nsT - ns
       atPairCycList(i2+2) = atPairCycList(i2+2) + nsT - ns
c
       if(CONTROL0)then
       write(kanalp,*)'nAtPairCyc : ',nAtPairCyc
       write(kanalp,*)
     & 'atPairCycList:i-j:',
     &  atPairCycList(i2+1),atPairCycList(i2+2) 
       end if !C
c
       if( i .lt. nCycLine ) goto 303
       end if !startLOOP
c       
c ENDLOOP block
c
      end if !MyResidue
      goto 100
      end if !$MTRES
      goto 100
c 
999     continue    ! EOFile
c
       if(CONTROL1)then
c print molZmatrx1,molZmatrx2,molZmatrx3
       write(kanalp,*)'readResLibSB: result: natom:',natom
       write(kanalp,*)'molZmatrx1,molZmatrx2,molZmatrx3:'
       do i = 1,natom
       i3 = 3*i-3
       i4 = 4*i-4 
       write(kanalp,'(i3,1x,3(a4,2x),4i4,4(f9.4,1x),a2)') i, 
     & (molZmatrx1(i3+k),k=1,3),(molZmatrx2(i4+k),k=1,4),    
     & ( molZmatrx3(i3+k),k=1,3),atomQ(i),ffatomName(i)
       end do !i
c
       write(kanalp, *) 'nAtPairCyc: ', nAtPairCyc 
       write(kanalp, '(12i6)')	
     & (atPairCycList(i), i = 1,2*nAtPairCyc)
c
       write(kanalp,*)'readResLib: final:'
       end if
c
	return
	end
