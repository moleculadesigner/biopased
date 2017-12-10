C* ** * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C*     SolventAtmParDef()                                                     *
C* Defines  the solvent atomic parameters                                     *
C* for Lazaridis & Karplus Gaussian shell                                     *
C* solvation model PROTEINS 1999, 35, 133-152                                 *
C*                                                                            *
C*     Yury Vorobjev 2002                                                     *
C* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	subroutine solvAtmParDef
     &             (solvAtmTypeFile,solvAtmDatFile,
     &              natom,atmNameEx,atmResName,nres,
     &              resName,startAtInRes,
     &              chNameRes,atmCharge,
     &              atmChargeNeutr,atmSolPar)
c
c INput:
C* solvAtmTypeFile - fullPath/fileName solvAtmType.dat file (aminAcids)
C* solvAtmDatFile - fullPath/fileName for solvGSPar.dat file 
C* natom
C* atmNameEx(ia) - atomNames Extended
C* atmResName(ia) - resName for atom ia
c* nres - number of residues in system
c* resName(ir) - list of residues names 1,..,nres
c* startAtInRes(ir) - Numb of the first atm in residue ir
C* chNameRes(ir) - chain name for residue list
C* atmCharge(ia) - atomic partial charges
c* OUT:result:
c* atmChargeNeutr(ia) - atomic partial charges Neutralized
C* atmSolPar(6*natom): ia: delGref,delV,lamb,alfa(*),R,qmod
C*                            qmod-neutralized set of atmic Q
        implicit none
	character*(*) solvAtmTypeFile,solvAtmDatFile
        integer natom,nres
        character*8 atmNameEx(*)     !extendedName atm+termTag
        character*4 atmResName(*)
        character*4 resName(*)
        character*1 chNameRes(*)
        integer startAtInRes(*)
        real atmCharge(*)
        real atmChargeNeutr(*)
        real atmSolPar(*)
c
        include "charStringSiz.h"
        include "output.h"
        include "kanalUse.h"
        include "statusMessg_mDyn.h"
c*      
        integer ntbmax                  !solvParTabsiz.h
        parameter (ntbmax=2047)         !MaxNumb of table lines for 
                                        !solvAtmTypeFile
        character*4 atm_tb(ntbmax)      !aacid atm - typeIn solAtmParDat()
        character*4 term_tb(ntbmax)
        character*4 res_tb(ntbmax)
        character*4 rnum_tb(ntbmax)
        character*1 chn_tb(ntbmax)
        integer partype(ntbmax)
        integer nline
        character*4 rnumb
        character*1 chname
        character*80 line
        integer nsolvAtmTypeMAX            !number of solvationAtomTypes
        parameter (nsolvAtmTypeMAX=30)
        character*4 solvGSAtmName(nsolvAtmTypeMAX)
        real solvGSParDat(5*nsolvAtmTypeMAX)     !atomType - solvParData() table
c
        integer atmresmax
        parameter (atmresmax=100)  !max atom in residue
        real atmqin(atmresmax),atmqout(atmresmax) 
        integer nat
c
        real qtRes
        real lambq
        logical find
        character*4 atmN,termN,aresN
        character*8 aresNtermTg
        integer ia,ja,i5,i6,j5,k
        integer ias,iaf,jas,jaf,ir
        integer contr
        integer kanalp
        logical CONTROL
        logical CONTROL0
        logical OPT_neutrQ
        integer natypetab
        integer linematch
        integer kanal11,kanal12

c CONTROL printOut
        kanal11 = kanal_solvAtType
        kanal12 = kanal_solvGSdat
        kanalp = kanalRunOut     ! print on standartOUT
        CONTROL0 = .true.
        CONTROL = .false. ! .true.
c
        OPT_neutrQ = .false.  !NO neutralization
c
        if(CONTROL0)then
        write(kanalp,*)'In solvAtmParDef:'
        write(kanalp,*)'solvAtmTypeFile:',solvAtmTypeFile
        end if
c* read in file  solvAtmTypeFile : solvGSPar_all_amb.dat
        open(unit=kanal11, file=solvAtmTypeFile, form='formatted',
     &       status= 'old')

        nline = 0

100     read(kanal11,'( a80 )', end=101 ) line 
        if(line(1:3) .eq. 'end' .or. 
     &      line(1:3) .eq. 'END')goto 101
c skip if comment characters ! or #
        if( line(1:1) .eq. '!' .or. line(1:1) .eq. '#')then
        goto 100 
        end if
c
        nline=nline+1
        if(nline .gt. ntbmax) then
        write(kanalp, *)'ERROR!in solvGSPar_all_amb.dat:',
     &  ' ntbmax (param) is SMALL :solvAtmParDef.f:increase ntbmax !'
c
        write(kanalPStat,*)mError, 
     &  '  ntbmax (param) is SMALL:solvAtmParDef.f:increase ntbmax !'
        stop
c
        end if

        read(line, '(a4,2x,a4,a4,1x,a1,3x,i2)') 
     &  atm_tb(nline),term_tb(nline),res_tb(nline),chn_tb(nline),
     &  partype(nline) 
        rnum_tb(nline)='    ' 

        if(CONTROL)then
        write(kanalp,'(a4,2x,a4,a4,1x,a1,3x,i2)') 
     &  atm_tb(nline),term_tb(nline),res_tb(nline),chn_tb(nline),
     &  partype(nline)
        end if     
    
        goto 100
       
101     continue
        close(unit=kanal11)
        natypetab=nline

c initialize
       chname=' '
       rnumb='    '

c control print
        if(CONTROL)then
        write(kanalp,*)'ia,atmNameEx(ia),atmResName(ia)'
        do ia = 1,natom
        write(kanalp,'(i6,2x,a8,a4)')
     &  ia,atmNameEx(ia),atmResName(ia)
        end do
        end if !CONTROL

c* read in file  solvAtmDatFile : solvGSPar.dat
        open(unit=kanal12, file=solvAtmDatFile, form='formatted',
     &       status= 'old')
c       
        if(CONTROL0)then
        write(kanalp,*)'solvAtmDatFile:',solvAtmDatFile
        end if
c        
        nline=0
200     read(kanal12,'(a80)', end=201)line
        if(line(1:3) .eq. 'end' .or. 
     &      line(1:3) .eq. 'END')goto 201

c        write(kanalp,'(a80)')line

        if(line(1:1) .eq. '!' .or. line(1:1) .eq. '#')then 
        goto 200
        end if

        nline=nline+1
        i5=(nline-1)*5
c
        if(nline .gt. nsolvAtmTypeMAX)then
        write(kanalp,*)'ERROR:solvAtmParDef:nsolvAtmTypeMAXMX is low!'
        stop
        end if
        read(line,'(a4,6x,5f8.3)')solvGSAtmName(nline),
     &                           (solvGSParDat(i5+k),k=1,5)

        if(CONTROL0)then
        write(kanalp,'(a4,6x,5f8.3)')
     &        solvGSAtmName(nline),(solvGSParDat(i5+k),k=1,5)         
        end if
c
        goto 200
c
201     continue
        close(unit=kanal12)
c
c do loop over all atoms and define solvAtmPar
       contr = 0
       do ia = 1,natom
c find atom record in table
       i6=(ia-1)*6
       atmN = atmNameEx(ia)(1:4)
       termN = atmNameEx(ia)(5:8)
       aresN = atmResName(ia) 
       aresNtermTg = aresN//termN
       atmSolPar(i6+1) = 0.0          
       atmSolPar(i6+2) = 0.0          
       atmSolPar(i6+3) = 0.0                         
       atmSolPar(i6+4) = 0.0                                       
       atmSolPar(i6+5) = 0.0                 
       atmSolPar(i6+6) = 0.0        
      
       if( atmN(1:1) .eq. 'H' )then
       find = .false.
       linematch = 1
       else
       call find_match_word(atmN,aresNtermTg,rnumb,chname,
     &  natypetab,atm_tb,term_tb,res_tb,rnum_tb,chn_tb,
     &  linematch,contr,find) 
       end if
                                           
       if(find) then
c assign parameters

       if(CONTROL)then
       write(kanalp,*)
     & 'partype(): linematch:',partype(linematch),linematch
       end if
c
       j5=(partype(linematch)-1)*5
       atmSolPar(i6+1) = solvGSParDat(j5+1)          
       atmSolPar(i6+2) = solvGSParDat(j5+2)          
       atmSolPar(i6+3) = solvGSParDat(j5+4)          
       atmSolPar(i6+4) = 0.0898 
     &                   *solvGSParDat(j5+3)/solvGSParDat(j5+4)
       atmSolPar(i6+5) = solvGSParDat(j5+5) 
       atmSolPar(i6+6) = 0.0      
c
        if(CONTROL)then 
        write(kanalp,'(i6,1x,a8,a4,6f8.3)')
     &  ia,atmNameEx(ia),atmResName(ia),
     &  (atmSolPar(i6+k), k=1,6)
        end if
c
       else
c print ERROR if not Hydrogen
       if(atmNameEx(ia)(1:1) .ne. 'H')then 
       write(kanalp,*)'ERROR:! atom ia:',ia, atmNameEx(ia),
     &    ' is not found in solvAtmDatFile:', solvAtmDatFile  
       end if
       end if
       end do !ia       

c neutralize charged residues
c modify lambda for charged residues
        if(CONTROL)then 
        write(kanalp,*)'* * * * * *  * * * * * * * * * * * * * '
        if(OPT_neutrQ)then
        write(kanalp,*)'solvAtmParDef: Neutralize Ionized Res: '
        else 
        write(kanalp,*)
     &  'solvAtmParDef: NO neutralizition of Ionized Res: '
        end if
        end if
c
        lambq = 6.0     ! Gauss shell for ionized res
        do ir = 1,nres
        ias = startAtInRes(ir)
        iaf = startAtInRes(ir+1) - 1
        nat = iaf-ias+1
c
        if(CONTROL)then 
        write(kanalp,'(a4,2x,i5,a10,2i5)')
     &        resName(ir),ir,' ias,iaf:',ias,iaf
        end if
c        
c collect atom charges for the RESir
        qtRes = 0.0
        do ia = 1,nat
        atmqin(ia) = atmCharge(ias+ia-1)
        atmqout(ia) = atmqin(ia)
        qtRes=qtRes + atmqin(ia)
        end do !ia
c
        if(OPT_neutrQ)
     &   call  neutralizeQ(nat,atmqin,qtRes,atmqout) 
        
        do ia = 1,nat
        i6 = 6*(ias+ia-1)-6
        if(OPT_neutrQ)then
        atmChargeNeutr(ias+ia-1) = atmqout(ia)
        atmSolPar(i6+6) = atmqout(ia) 
        else
        atmSolPar(i6+6) = atmqin(ia)
        atmChargeNeutr(ias+ia-1) = atmqin(ia) 
        end if
c
c modify parameter lambda for charged residues
        if(abs(qtRes) .ge. 0.999999) then
         atmSolPar(i6+3) = lambq
        end if 
        end do !ia
        end do !ir 

c control print
        if(CONTROL0)then
       write(kanalp,*)
       write(kanalp,*)'SolvationGSmodel Parameters : ' 
       write(kanalp,'(a6,1x,a8,a4,6(a5,4x))')
     & 'atomNb','atomName','res','   dV','dGref ',
     & ' lamb ',' alfa', 'Ra  ','qat '
       do ia = 1,natom
        i6 = 6*ia-6
        write(kanalp,'(i6,1x,a8,a4,6f8.3)')
     &  ia,atmNameEx(ia),atmResName(ia),
     &  (atmSolPar(i6+k), k=1,6)
        end do
c
       write(kanalp,*)'End of solvAtmParDef: '
       end if

	return
        end
