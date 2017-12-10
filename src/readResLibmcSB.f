c
c Yury Vorobjev 2002-2004
c
c               Okt. 2006: arbitrary type of res-res linking
c                          topo-defined
c
c readResMTLibrary : unitResLib
c find the given residue by name resName6
c assemble molecular ZMatrix by adding the ZMatrix for the Residue
c to the existing ZMatrix
c MULtyChain/MULtyMolec version
c
	subroutine readResLibBSmc(unitResLib,iread,
     &              resName6,resTypeIn,resFound,
     &              dumFlag,atom3ToLink,natom,natinres,nHvyAtinres,
     &              molZmatrx1,molZmatrx2,molZmatrx3,
     &              atomQ,ffatomName,
     &              nAtPairCycMX,nAtPairCyc,atPairCycList)
c
c IN:
c unitResLib=unit file, iread=resNumbtoRead
c resName6 (char*6)
c resTypeIn,
c dumFlag - 0/1 NO/Yes to remove DUM atoms of the Res from final zMatrix
c
c IN/OUT:
c nres - number of RES before EXecution, and is total n Residues
c        after final RETURN
c natom - number of atoms before EXecution, and is a total n atoms 
c         natom = natom + natinres   at the final RETURN  
c natinres - number of atom in this residue
c molZmatrx1():atname,atFFName,BackSidename :3*natom
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
       integer iread           
       integer natinres,natom
       integer nHvyAtinres
       character*6 resName6   
       character*9 resTypeIn(*)    
       character*12 atom3ToLink(*)   !three atoms to link with 3 DUMM atom of the NEXT residue
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
        logical resFound 
c local
       integer qColMAX
       parameter (qColMAX = 3)
       integer natresMAXloc
       parameter (natresMAXloc = 100) 
c dummAtomNames: for AA and NA seq
       integer iPolymLinkTypeMAX
       parameter (iPolymLinkTypeMAX=4) ! P-P,NA-NA,NA-P,TOPO
       character*4 duAtName(3,iPolymLinkTypeMAX)
       integer polyType
       character*3 resTypeInL
c
       integer LCharMAX
       parameter (LCharMAX=88)
       integer nl,nin
       integer natresMXloc
       character*(LCharMAX) linein
       character*9 rType   !INT  AA, internal AminoAcid
       integer ns,nsN,nsCA, nsC,nsT
       integer naM,naCA,naC,naO,naN
       integer naH,natl
       integer nCycLine
       integer i,i2,i3,i4,i5,i6,k
       integer n1,nrimp,ni
       real qres
       real qColmn(qColMAX)
       integer qCol
cx       logical resFound
       logical OPT_removeDUMM
       integer kanalp
       logical CONTROL2,CONTROL1,CONTROL0
c to link with 3 DUMM 1,2,3 atom of the NEXT residue
       data duAtName/ "N   ","CA  ","C   ",     ! AAcid->AA
     &                "C4' ","C3' ","O3' ",     ! NAcid->NA
     &                "C3' ","O3' ","P3' ",     ! NA->-AA 
     &                "TOPO","TOPO","TOPO"/     ! TOPO defined general LINK
c
       CONTROL2 = .false.   !
       CONTROL1 = .false.  !regular
cx2009       CONTROL0 = .false.   !
       CONTROL0 = .true.
       if(OPT_OUTfull)CONTROL0 = .true.
       resFound = .false.
       kanalp =  kanalRunOut
       OPT_removeDUMM = .false. 
       if(dumFlag .eq. 1) OPT_removeDUMM = .true.
c
       if(CONTROL0)then
       write(kanalp,*)'readResLibmc:w01:START!readResLib:ResName:',
     & resName6(1:4)
       end if
c
       rewind unitResLib
c init
       nHvyAtinres = 0
       nCycLine = 0
       natresMXloc=natresMAXloc ! max atoms per RES 
       polyType = 0   !  0-LG, 1-AA, 2-NA            
       nl = 0
       nin = 0
       ns = 0        !natm to remove from resid  when assemble Zmatrx
       nsN = 0       !nN(i-1)    atom Numbers to assemble polymer chain
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
      read(unitResLib,'(a88)', end = 999 )linein       ! 88=LCharMAX
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
c
      read(unitResLib,'(a80)', end = 999 )linein
c
      if(linein(1:1) .eq. '!' .or. linein(1:1) .eq. '#')goto 101 
      nin = nin + 1  ! number of informative line
c
      if(linein(1:4) .eq. resName6(1:4) )then    ! startMyRes
c startMyRes
      resFound = .true.
c read MyRESidue lines
c
      if(CONTROL1)then
       write(kanalp,*)'readResLibmc:W2:MyREs',linein,' resName6=',
     & resName6(1:4)
       write(kanalp,*)'resFound =',resFound
      end if
c
       qCol = 0
       read(linein, '(5x,a9,1x,i3,f8.4,1x,i2)')rType,natl,qres,qCol
c
       resTypeInL = rType(1:3)     ! defineResType: INT,BEG,FIN
       resTypeIn(iread)=rType      ! extResType:  INT,BEG,FIN + AA,NA,NAP
c
c define atom3ToLink(iread) :
       atom3ToLink(iread) = "            "
       if(rType(6:9) .eq. "AA  ") atom3ToLink(iread)="N   CA  C   "
       if(rType(6:9) .eq. "NA  ") atom3ToLink(iread)="C4' C3' O3' "
       if(rType(6:9) .eq. "NAP ") atom3ToLink(iread)="C3' O3' P3' "
       if(rType(6:9) .eq. "TOPO")
     &   read(linein, '(31x,a12)')atom3ToLink(iread)
c
       if(CONTROL0)then
       write(kanalp,*)
     & 'readResLib: rType,natl,qres,qCol:',
     & rType,natl,qres,qCol
c
       write(kanalp,*)'readResLibmc:W3:atom3ToLink(iread): iread=',
     & iread,' ', atom3ToLink(iread)
       end if !C0 
c control
      if(natl .gt. natresMXloc)then
      write(kanalp,*)'readResLib:ERROR! natresMXloc is low!!'
      natresMXloc = natl + 6
c      stop
      end if
c
c define atom ReNumberingRulles
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
c find atomNumbers Shifts to make atom reNumbering
c find numbers in preceeding AA (i-1) residue: nN,nCA,nC,ntot
c                            NA (i-1)        : nC4',C3',O3'
       if(iread .gt. 1)then
       if(resTypeIn(iread-1)(6:7) .eq. 'AA') polyType = 1
       if(resTypeIn(iread-1)(6:7) .eq. 'NA') polyType = 2  !if NuclAcid 
       if(resTypeIn(iread-1)(6:9) .eq. 'NAP ')polyType = 3 !after NAPLinker
       if(resTypeIn(iread-1)(6:9) .eq. 'TOPO')then
       polyType = 4 ! after TOPO defined linker
       duAtName(1,polyType) = atom3ToLink(iread-1)(1:4)
       duAtName(2,polyType) = atom3ToLink(iread-1)(5:8)
       duAtName(3,polyType) = atom3ToLink(iread-1)(9:12)
       end if!TOPO defined LinkerAtoms
       end if!iread .gt. 1
c
      if(natom .gt. 0)then
      if(natom .gt. natresMXloc)then
      n1 = natom-natresMXloc+1
      else
      n1 = 1
      end if !natom .gt. natresMXloc
c
      do i=n1,natom
      i3 = 3*i-3
      if (molZmatrx1(i3+1) .eq. duAtName(1,polyType))nsN=i
      if (molZmatrx1(i3+1) .eq. duAtName(2,polyType))nsCA=i
      if (molZmatrx1(i3+1) .eq. duAtName(3,polyType))nsC=i
      end do!i
c
      if(CONTROL1)then
      write(kanalp,*)'C1:Shift: nsN, nsCA, nsC:',nsN, nsCA, nsC
      write(kanalp,*)'duAtName(1,polyType): ',  
     & (duAtName(k,polyType),k=1,3), ' polyType=',polyType 
      end if !C1
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
c  read next natl lines as Zmatrix  
      i = 0
103   read(unitResLib,'(a88)')linein
c
      if(CONTROL1)then
      write(kanalp,*)'W3:MyREs',linein
      end if
c
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
      read(linein, '(6x, a4,2x,a4,2x,a3)')
     &   molZmatrx1(i3+1),molZmatrx1(i3+2),
     &   molZmatrx1(i3+3)
c     
       if(molZmatrx1(i3+1)(1:2) .ne. 'DU' 
     &    .and. molZmatrx1(i3+1)(1:1) .ne. 'H')then
       nHvyAtinres = nHvyAtinres + 1
       end if ! notH
c 
      read(linein,'(i4,16x,i4,i4,i4)')
     &  molZmatrx2(i4+1),molZmatrx2(i4+2),
     &  molZmatrx2(i4+3),molZmatrx2(i4+4)
c
        if(CONTROL1)then
        write(kanalp,*)'before Sh:nat=',natom,' zm1,2,3,4:',
     &  molZmatrx2(i4+1),molZmatrx2(i4+2),
     &  molZmatrx2(i4+3),molZmatrx2(i4+4)
        end if
c
c make a shift in atomNumbering to add RES to PolymerChain 
c                               to start NEW PolymerChain
c                               to start NEW IsolatedMolec
       if(resTypeInL .eq. 'INT' .or. 
     &    resTypeInL .eq. 'FIN') then ! RES is added to the Polymer CHain
c                                    ! make the SHift in the atom numbering
c       
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
        else            ! if start NEW CHain with first RESidue
c
cx         if(OPT_removeDUMM)then
          do k=1,4
          molZmatrx2(i4+k) = molZmatrx2(i4+k)-ns
          if(molZmatrx2(i4+k) .ge. 1)then
          molZmatrx2(i4+k) = molZmatrx2(i4+k) + nsT
          end if
          end do!k
cx         end if! _removeDUMM
c
      end if ! SHift
c
       read(linein, '(36x,4f9.4)')
     & molZmatrx3(i3+1),molZmatrx3(i3+2),
     & molZmatrx3(i3+3),
     & atomQ(natom)
c
c read qAtom from defined qColomn
       if(qCol .gt. 0)then
       if(qCol .gt. qColMAX)then
       write(kanalp,*)'readResLibBS ERROR!: qCol>qColMAX:',qCol,qColMAX
       stop
       end if
c read 3 atQ colomns:
       read(linein, '(62x,3f9.4)')qColmn(1),qColmn(2),qColmn(3)
       atomQ(natom) = qColmn(qCol)
       end if
c
       ffatomName(natom) = molZmatrx1(i3+2)(1:2)
cx       molZmatrx1(i3+3) = resName6(1:4)        !jan 16.03
c
      end if ! readResLibAtom
c
      if( i .lt. natl)goto 103
c 
c read LOOP (closure) lines
c
301   nl = nl + 1
      read(unitResLib,'(a80)', end = 999 )linein 
c
      if(CONTROL1)then
      write(kanalp,*)'W2:MyREs:LOOPblock: nl:',nl
      write(kanalp,*)linein
      end if
c
       if(linein(1:1) .eq. '!' .or. linein(1:1) .eq. '#')goto 301
       nin = nin + 1  ! number of informative line
c
       if(linein(1:4) .eq. "LOOP" )then    ! startLOOP 
cstartLOOPblock
       read(linein, '(4x,i4)')nCycLine
c read number of LOOP lines=nCycLine
c
       if(CONTROL0)then
       write(kanalp,*)'readResLibmc:W4:readResLib: rType,nCyclLine:',
     & rType,nCycLine
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
       write(kanalp,*)'.increase parameter nAtPairCycMAX in zmatrix.h'
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
       if(CONTROL1)then
       write(kanalp,*)'nAtPairCyc : ',nAtPairCyc
       write(kanalp,*)
     & 'atPairCycList:i-j:',
     &  atPairCycList(i2+1),atPairCycList(i2+2) 
       end if !C
c
       if( i .lt. nCycLine ) goto 303
       end if !startLOOP
c ENDLOOP bloc
c
       if(linein(1:4) .ne. "DONE")goto 301
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
       write(kanalp,*)'readResLibBS: result: natom:',natom
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
       if(CONTROL0)then
c
       if(resFound)
     &  write(kanalp,*)
     & 'readResLibmc: Finish: rType,nCyclLine,nHvyAtinres:',
     & rType,nCycLine,nHvyAtinres
c
       if(.not.resFound)then
       write(kanalp,*)
     & 'WARNING! res=',resName6(1:4),
     & ' does NOT exist in topo.dat file!:unit=', unitResLib
       end if

       end if !CONTROL0
c
	return
	end
