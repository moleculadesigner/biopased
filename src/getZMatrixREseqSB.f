c generate ZMatrix for REsidue (AA+NA)  sequence 
c         multiCHain, multyMolecular
c         PDB names etc. 
c         pair12List
c
c  Yuri Vorobjev, 2002, 2006, 2007
c
        subroutine getZMatrixREseq(nresREseq,resNameREseq,      
     &          resTypeInREseq,atom3ToLinkInREseq,atomBlockName,
     &              natREseq,atNameREseq,ffatNameREseq,
     &              atResNameREseq,atResNumbREseq,atNameREseqEx,
     &              startAtResREseq,stopAtResREseq,
     &              nChain,startAtInCha,stopAtInCha,
     &              startResInCha,stopResInCha,
     &              nMolec,startAtInMol,stopAtInMol,
     &              startResInMol,stopResInMol,
     &              Pair12ListREseq,startPairL12REseq,
     &              nPairL12REseq,atQREseq)         
c
c IN:
c nresREseq - number Res in REsidue sequence
c resNameREseq(ir), ir=1,nresREseq
c resTypeInREseq(ir) : BEG,ISO,INT,FIN residue type in sequence  
c nChain,startResInCha(ic),stopResInCha(ic) - start/stop Res in Chain ic,
c nMolec,startResInMol(im),stopResInMol(im) - start/stop Res in Molec im
c
c OUT:
c natREseq - total NAtoms 
c atNameREseq(ia) - atom names PDB
c atomBlockName(ia) - blockName(Group) for atom ia
c ffatNameREseq(ia) ff atom names
c atResNameREseq(ia)
c atResNumbREseq(ia)
c startAtResREseq(ir) - starting atom of residue ir
c stopAtResREseq(ir) stop atom in Res ir
c startAtInCha(ic),stopAtInCha(ic) - start/stop At  in Chain ic
c startAtInMol(im),stopAtInMol(im) - start/stop At  in Molec im 
c 
c Pair12ListREseq(iL) - list of pair 12, symmetrical
c startPairL12REseq(ia) - start in the Pair12ListREseq() for atom ia
c nPairL12REseq(ia) - number of Vbond neighbours in the 12List
c 
        include 'xyzPDBsize.h'
        integer nresREseq
        character*6 resNameREseq(*)
        character*4 atNameREseq(*)
        character*4 atResNameREseq(*)
        character*8 atNameREseqEx(*)
        integer natREseq
        character*2 ffatNameREseq(*)
        integer startAtResREseq(*)
        integer stopAtResREseq(*)
        integer atResNumbREseq(*)
c
        integer nChain
        integer startAtInCha(*),stopAtInCha(*)
        integer startResInCha(*),stopResInCha(*) 
        integer nMolec
        integer startAtInMol(*),stopAtInMol(*)
        integer startResInMol(*),stopResInMol(*)
c
        integer nPairL12REseq(*)
        integer Pair12ListREseq(*)
        integer startPairL12REseq(*)
        character*9 resTypeInREseq(*)
        character*12 atom3ToLinkInREseq
        character*4 atomBlockName(*)
        real atQREseq(*)
c
        include "output.h"
        include "kanalUse.h"
        include "statusMessg_mDyn.h"
c        
        include 'zmatrix.h'
cx	include "keyWordOpt.h"
        include 'filedat.h'
        include 'ssbond.h'
c
c local
        integer nAtPairCycMX
        integer nres,np12L
        integer ir,ia,irs,irf
        integer jz,jr,jz4
        integer k,i,i2,i3,i4
        integer ic,ic2
        integer kanalp
        integer unitResLib,unitLib
        integer unitResLibB
        integer unitResLibE
        integer unitResLibLIG
        character*80 line
        integer natinres
        integer dumFlag
        logical resFound             ! 2007 YV
        logical CONTROL,CONTROL1
c
        kanalp = kanalRunOut
        CONTROL = .true.
        CONTROL1 = .false.
c : Linux:
       zmfileINT = HOMEBS_dir(1:ndirLett)//'/dat/bs_int_all94.dat'
       zmfileBEG = HOMEBS_dir(1:ndirLett)//'/dat/bs_one_all94.dat'
       zmfileFIN = HOMEBS_dir(1:ndirLett)//'/dat/bs_fin_all94.dat'
       zmfileLIG = HOMEBS_dir(1:ndirLett)//'/dat/bs_lig_all94.dat'
c : DOS
cx yv21.11.08
cx       zmfileINT = HOMEBS_dir(1:ndirLett)//'bs_int_all94.dat'
cx       zmfileBEG = HOMEBS_dir(1:ndirLett)//'bs_one_all94.dat'
cx       zmfileFIN = HOMEBS_dir(1:ndirLett)//'bs_fin_all94.dat'
cx       zmfileLIG = HOMEBS_dir(1:ndirLett)//'bs_lig_all94.dat'
c
       write(kanalp,*)'ZMfileINT:',zmfileINT
       write(kanalp,*)'ZMfileBEG:',zmfileBEG
       write(kanalp,*)'ZMfileFIN:',zmfileFIN
       write(kanalp,*)'ZMfileLIG:',ZMfileLIG
       write(kanalp,*)'defFlag_zmfileLIG=',defFlag_zmfileLIG
c
c keep zmfileLIG taken from -Lt command line
       if(defFlag_zmfileLIG .eq. 1) then 
c! take LigTopoData from -tL command line
        zmfileLIG = zmfileLIG                              
	iStatus = iStatus + 1                
        write(kanalPStat,*)mStatus,iStatus,messageStatus
     &  ,' LigTopoFile s taken from command line = ...', zmfileLIG
       end if !
c
       write(kanalp,*)'ZMfileLIG:ligandTopoFile=',zmfileLIG
c
       unitResLib = kanal_ResTop01
       unitResLibB = kanal_ResTop02
       unitResLibE = kanal_ResTop03
       unitResLibLIG = kanal_ResTop04
c
       open(unit=unitResLib,file=zmfileINT,form='formatted',
     &      status='old')
c
       open(unit=unitResLibB,file=zmfileBEG,form='formatted',
     &      status='old')
c
       open(unit=unitResLibE,file=zmfileFIN,form='formatted',
     &      status='old')
c
cx       open(unit=unitResLibLIG,file=zmfileLIG,form='formatted',
cx     &      status='old')
c
c read RESidues and assembe zMatrix for the AASequence
c
       nres = nresREseq
       natREseq= 0
       natPairCyc = 0
       natPairCycMX = natPairCycMAX
c
       dumFlag = 1     !remove all DUMM atoms from Zmatrix
       do ir = 1,nres
c
       startAtResREseq(ir) = natREseq+1
c
       atHvyNbInResZm(ir) = 0   !init
c
c define REsLib if start/stop REsidue in Chain
        if(resTypeInREseq(ir)(1:3) .eq. 'BEG' .or. 
     &    resTypeInREseq(ir)(1:3) .eq. 'ISO') unitLib=unitResLibB
c
        if(resTypeInREseq(ir)(1:3) .eq. 'INT') unitLib=unitResLib
        if(resTypeInREseq(ir)(1:3) .eq. 'FIN') unitLib=unitResLibE        
c           
       call readResLibBSmc(unitLib,ir,
     &              resNameREseq(ir),resTypeInREseq,resFound,
     &              dumFlag,atom3ToLinkInREseq,
     &              natREseq,natinres,atHvyNbInResZm(ir),
     &              molZmatrx1,molZmatrx2,molZmatrx3,
     &              atQREseq,ffatNameREseq,
     &              natPairCycMX, natPairCyc,atPairCycList)
c
c YV2009 read LIGmainLib topo 
         if(.not.resFound ) then
         if(resTypeInREseq(ir)(1:3) .eq. 'BEG' .or.
     &    resTypeInREseq(ir)(1:3) .eq. 'ISO') then   !read LIG topo
          unitLib=unitResLibLIG
c
          open(unit=unitResLibLIG,file=zmfileLIG,form='formatted',
     &      status='old')
c
         call readResLibBSmc(unitResLibLig,ir,
     &              resNameREseq(ir),resTypeInREseq,resFound,
     &              dumFlag,atom3ToLinkInREseq,
     &              natREseq,natinres,atHvyNbInResZm(ir),
     &              molZmatrx1,molZmatrx2,molZmatrx3,
     &              atQREseq,ffatNameREseq,
     &              natPairCycMX, natPairCyc,atPairCycList)
c
          close(unit=unitResLibLig)
          end if !unitResLibLig
          end if !read zmfileLIG
c
c YV2009 read -tL LIGtopo.dat file
         if(.not.resFound ) then
         write(kanalPStat,*)'NO RESidue TopoDAT !!'
         stop
         end if
cx       if(.not.resFound ) then
cx         if(resTypeInREseq(ir)(1:3) .eq. 'BEG' .or.
cx     &    resTypeInREseq(ir)(1:3) .eq. 'ISO') then   !read -tL LIG topo
c keep zmfileLIG taken from -Lt command line
cx NOTE! no DOCKING in the bioPASED
cx
cx       if(defFlag_zmfileLIG .eq. 1) then
c! take LigTopoData from -tL command line
cx        iStatus = iStatus + 1
cx        write(kanalPStat,*)mStatus,iStatus,messageStatus
cx     &  ,' LigTopoFile taken from command line = ', ligMolTopDatFile
cx       end if !
c
cx        open(unit=unitResLibLIG,file=ligMolTopDatFile,form='formatted',
cx     &   status='old')
c
cx         call readResLibBSmc(unitResLibLig,ir,
cx     &              resNameREseq(ir),resTypeInREseq,resFound,
cx     &              dumFlag,atom3ToLinkInREseq,
cx     &              natREseq,natinres,atHvyNbInResZm(ir),
cx     &              molZmatrx1,molZmatrx2,molZmatrx3,
cx     &              atQREseq,ffatNameREseq,
cx     &              natPairCycMX, natPairCyc,atPairCycList)
c
cx          close(unit=unitResLibLig)
cx          end if !unitResLibLig
cx          end if !read zmfileLIG
c
c
          if(.not.resFound)then
         write(kanalp,*)'ERROR!getZMatrixREseq: res:',resNameREseq(ir),
     &   ' is NOT FOUND in TOPO Lib: zmfileLIG = ',zmfileLIG
         stop
          end if!.not.resFound
c
       stopAtResREseq(ir) = startAtResREseq(ir)+natinres-1
c
       end do! ir
c
       startAtResREseq(nres+1) = stopAtResREseq(nres) + 1
c
       close(unitResLib) 
       close(unitResLibE) 
       close(unitResLibB) 
c
c define startAtInCha(),stopAtInCha(),startAtInMol(),stopAtInMol()
       do i = 1,nChain
       startAtInCha(i) = startAtResREseq(startResInCha(i))
       stopAtInCha(i) = stopAtResREseq(stopResInCha(i))
       end do!i
c
       do i = 1,nMolec
       startAtInMol(i) = startAtResREseq(startResInMol(i))
       stopAtInMol(i) = stopAtResREseq(stopResInMol(i))
       end do!i
c
c define atomPDB names, etc
       do ia = 1,natREseq
       atNameREseq(ia) = molZmatrx1(ia*3-2)
       atomBlockName(ia) = molZmatrx1(ia*3)
       end do !ia
c
       if(CONTROL .and. OPT_OUTfull)then
       write(kanalp,*)'getZMatrixREseq:'
       do ir=1,nresREseq
       write(kanalp,*)
     & 'ir,startAtResREseq(ir),stopAtResREseq(ir):',
     &  ir,startAtResREseq(ir),stopAtResREseq(ir) 
       end do!ir
c
       do ir=1,nChain
              write(kanalp,*)
     & 'iCh,startAtinCha(ir),stopAtInCha(ir):',
     &  ir,startAtInCha(ir),stopAtInCha(ir)
       end do
c
       do ir=1,nMolec
              write(kanalp,*)
     & 'iMol,startAtinMol(ir),stopAtInMol(ir):',
     &  ir,startAtInMol(ir),stopAtInMol(ir)
       end do
c
       end if!C
c
c define 12PairList
c
       np12L = 0
c
c bigLoop over RES and atoms in RES
       do ir = 1,nresREseq
       do ia = startAtResREseq(ir),stopAtResREseq(ir) 
       atResNameREseq(ia) = resNameREseq(ir)
       atResNumbREseq(ia) = ir
       startPairL12REseq(ia) = np12L+1
       nPairL12REseq(ia) = 0
c
       atNameREseqEx(ia) = atNameREseq(ia)//'    '
c
       if(resTypeInREseq(ir)(1:3) .eq. 'BEG')then
       atNameREseqEx(ia) = atNameREseq(ia)//'NT  '  
       end if
       if(resTypeInREseq(ir)(1:3) .eq. 'FIN')then
       atNameREseqEx(ia) = atNameREseq(ia)//'CT  ' 
       end if
c add blockName to atNameREseqEx(ia):
cx       atNameREseqEx(ia)(8:8) = molZmatrx1(ia*3)(4:4)
c
c define irs,irf: start/final REs to search vbonds
        if(resTypeInREseq(ir)(1:3) .eq. 'BEG')irs = ir
        if(resTypeInREseq(ir)(1:3) .eq. 'ISO')then
        irs = ir
        irf = ir
        end if ! ISO
        if(resTypeInREseq(ir)(1:3) .eq. 'INT')irs = ir-1
        if(resTypeInREseq(ir)(1:3) .eq. 'FIN')irs = ir-1
c
        if(resTypeInREseq(ir)(1:3) .eq. 'BEG')irf = ir+1
        if(resTypeInREseq(ir)(1:3) .eq. 'INT')irf = ir+1 
        if(resTypeInREseq(ir)(1:3) .eq. 'FIN')irf = ir
c
c control:
        if(irs .lt. 1 )then
        write(kanalp,*)'getZMatrixREseq: ERROR! in irs:', irs
        stop
        end if
        if(irf .gt. nresREseq)then
c       irf = nresREseq
        write(kanalp,*)'getZMatrixREseq: ERROR! in irf:', irf
        write(kanalp,*)
     &  'getZMatrixREseq:ERROR! no TERCHA/ENDMOL lines in pdb!?'
c
        write(kanalPStat,*)mError,
     &  'getZMatrixREseq:ERROR! no TER/END lines in (-c in.pdb) !?'
c
        stop
        end if
c
       do jr = irs,irf
       do jz = startAtResREseq(jr),stopAtResREseq(jr)
       jz4=4*jz-4
       if(ia .eq. molZmatrx2(jz4+1)) then
       if(molZmatrx2(jz4+2) .ge. 1)then
       np12L = np12L+1
       Pair12ListREseq(np12L) = molZmatrx2(jz4+2)
       nPairL12REseq(ia) = nPairL12REseq(ia)+1
       end if
       end if
c       
       if(ia .eq. molZmatrx2(jz4+2) )then
       np12L = np12L+1
       Pair12ListREseq(np12L) = molZmatrx2(jz4+1)
       nPairL12REseq(ia) = nPairL12REseq(ia)+1
       end if
c
       end do !jz
       end do !jr
c
c add cycleClosure pairs in 12PairList
       do ic = 1,natPairCyc
       ic2 = 2*ic-1 
       if(atPairCycList(ic2) .eq. ia)then
       np12L = np12L+1
       Pair12ListREseq(np12L) = atPairCycList(ic2+1)
       nPairL12REseq(ia) = nPairL12REseq(ia)+1
       end if !
c
       if(atPairCycList(ic2+1) .eq. ia)then
       np12L = np12L+1
       Pair12ListREseq(np12L) = atPairCycList(ic2)
       nPairL12REseq(ia) = nPairL12REseq(ia)+1
       end if !
       end do !ic
c
c add S-S bonds
       do ic = 1,nSSbonds  
       ic2 = 2*ic-1
       if(ssBondAt12List(ic2) .eq. ia)then
       np12L = np12L+1
       Pair12ListREseq(np12L) = ssBondAt12List(ic2+1)
       nPairL12REseq(ia) = nPairL12REseq(ia)+1
       end if ! ia
c
       if(ssBondAt12List(ic2+1) .eq. ia)then
       np12L = np12L+1
       Pair12ListREseq(np12L) = atPairCycList(ic2)
       nPairL12REseq(ia) = nPairL12REseq(ia)+1
       end if !ia
       end do !ic 
c endS-S
c
       end do !ia bigLoop
       end do !ir bigLoop
c
       if(CONTROL)then
       write(kanalp,*)'getZMatrixREseqe: done: natREseq:',natREseq
       end if!C
       if(CONTROL .and. OPT_OUTfull)then
c print molZmatrx1,molZmatrx2,molZmatrx3
cx      write(kanalp,*)'getZMatrixREseqe: natREseq:',natREseq
       write(kanalp,*)'molZmatrx1,molZmatrx2,molZmatrx3,atomBlockName:'
       do i = 1,natREseq
       i3 = 3*i-3
       i4 = 4*i-4
       write(kanalp,'(i4,1x,3(a4,2x),4i4,4(f9.4,1x),a2,1x,a4)') i,
     & (molZmatrx1(i3+k),k=1,3),(molZmatrx2(i4+k),k=1,4),
     & ( molZmatrx3(i3+k),k=1,3),atQREseq(i),ffatNameREseq(i),
     & atomBlockName(i)
       end do !i
       end if !Control
c
       if(CONTROL1)then
       write(kanalp,*)'getZMatrixREseq: aaSeq Atom List:'
       do ia = 1,natREseq
       write(kanalp,'(a4,2x,i6,1x,a4,a4,1x,i6)')
     & 'ATOM',ia,atNameREseq(ia),atResNameREseq(ia),atResNumbREseq(ia)
       end do !ia
       end if !Control
c
        if(CONTROL1)then
        write(kanalp,*)'getZMatrixREseq:Result:'
        write(kanalp,*)'startPairL12:'
        call print10(natREseq,startPairL12REseq)
        write(kanalp,*)'nPairL12:'
        call print10(natREseq,nPairL12REseq)
        write(kanalp,*)'pair12List:np12L:',np12L
        call print10(np12L,pair12ListREseq)
        end if !Control
c 
       return
       end
