c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c*  getVDWatMass()                                                            *
c* Defines  the ForceField  parameters                                        *
c*                                                                            *
c*     Yury Vorobjev 2002                                                     *
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
        subroutine getVDWatMass(ffParFile,
     &              natom,atomNameEx,ResName,chName,ffAtomName,
     &             nVDWtypeMX,nVDWtype,atomVDWtype,
     &                                     atomVDW12ab,atomMass)
c
c InPut:
c       ffParFile - ffParameters file 
c       ffAtomName(ia) - FFatomName to search table
c PDBinfo: natom,atomNameEx,ResName,chName,ffAtomName
c
c RESULT: nVDWtype,atomVDWtype(ia)- vdwType of atom ia,
c         atomVDW12ab(t12) = A,B vdw coeff for all pair of vdwType 
c         atomMass(iA) = mass of atom ia
c
        implicit none
        character*(*) ffParFile     
        integer natom
        character*8 atomNameEx(*)
        character*4 ResName(*)
        character*1 chName(*)
        character*2 ffAtomName(*)
        integer atomVDWtype(*),nVDWtype
        integer nVDWtypeMX
        real atomVDW12ab(*)
        real atomMass(*)
c
        include "charStringSiz.h"
        include "output.h"
        include "kanalUse.h"
        include "statusMessg_mDyn.h"
clocal
        integer nVDWatomParMX
        integer ntbmax                  !convertionTablsiz.h
        parameter (ntbmax=2047)         !MaxNumb of table lines
c
        integer tlen
        parameter (tlen = 2)        
        character*(tlen) atName_tb(ntbmax)
        real    rm_tb(ntbmax)
        real    ev_tb(ntbmax)
        real    ams_tb(ntbmax)
        integer atmType_tb(ntbmax)
        character*20 line_tb(ntbmax) 
c hashTable
        integer nhashmx
        parameter (nhashmx = ntbmax )
        integer ntabl_ihash(nhashmx)
        integer hashtblink(nhashmx)    
        integer nhash
c
        character*80 line
        character*(tlen) aabb,bbaa    
        character*2 aa,bb
        character*9 recBlock
        character*4 endBlock 
        integer i,i1,i2,len,i4
        integer nrb,t1,t2,t12,p4
        integer ntbrec,contr
        integer nstr1
        integer kanalp
        integer linematch
        real Avdw12,Bvdw12,rm12,ev12
        logical find
        logical myblock
        logical CONTROL,CONTROL1
        logical OPT_hash
c
c initialize
        kanalp = kanalRunOut
        CONTROL  = .true.  
        CONTROL1 = .false.
        OPT_hash = .true.

c read in ffParfile
        open(unit=11, file=ffParFile, form='formatted',
     &       status= 'old')
c
        if(CONTROL)then
          write(kanalp,*)'In  getVDWatMass :'
          write(kanalp,*)'ffParFile : ',ffParFile
        end if
c
c READ ATOMMASS:
        recBlock = '$ATOMMASS'
        nrb = 9
        endBlock = '$END'
        myblock = .false.
        ntbrec = 0
c
100     read(11,'( a80 )', end=101 ) line
c
        if(CONTROL1)then
        write(kanalp,*)line
        write(kanalp,*) 'myblock :', myblock 
        end if
c
        if(myblock .and. (line(1:4) .eq. endBlock) ) goto 101
        if(line(1:nrb) .eq. recBlock )then
         myblock = .true.  
         goto 100
         end if 
c
        if(.not.myblock)goto 100
c skip if comment characters ! or #
        if( line(1:1) .eq. '!' .or. line(1:1) .eq. '#')then
        goto 100
        end if
c
        if(myblock)then
        ntbrec=ntbrec+1
        if(ntbrec .gt. ntbmax) then
        write(kanalp, *)'ERROR!in getVDWatMass:',
     &  ' ntbmax (param) is SMALL '
        stop
        end if
c 
        read(line, '(a2,1x,f7.3)')
     &  atName_tb(ntbrec),ams_tb(ntbrec)
c
        if(CONTROL1)then
	write(kanalp,'(a2,1x,f7.3,i6)')
     &  atName_tb(ntbrec),ams_tb(ntbrec),ntbrec
        end if
c
        end if !myblock
c
        goto 100
c
101     continue
c end of readIn ffParFile
c
c calculate hashTable for dictionary table
       len = tlen
       nstr1 = len
       contr = 0     !control print
c init hashtable
       do i = 1,nhashmx
       hashtblink(i) = 0
       ntabl_ihash(i) = 0
       end do!i
c 
       do i = 1,ntbrec
       line_tb(i) = atName_tb(i)   ! string length
c calculate hash number from string strtabl1(1:nstr1) simbols
       call hashCalculator(atName_tb(i),nstr1,nhashmx,nhash,contr)
 
c hashtblink() - linkedList to resolve degenerate nhash
c ntabl_ihash() - head() of linked list for given nhash
 
       hashtblink(i) = ntabl_ihash(nhash)
       ntabl_ihash(nhash) = i
 
       if(CONTROL1)then
       write(kanalp,'(a27,i6,i6,a5,a5)')
     & 'hash_dict_tabl:itab,nhash:',i,nhash,' str:',
     &      atName_tb(i)
       end if

       end do !i
c
c find match for atom
        if(CONTROL1)then
        write(kanalp,*)'getVDWatMass: result: atomMass: '
        write(kanalp,*)'N name  Mass  ms_tb :'
        end if
c
        do i = 1, natom
        aa = ffAtomName(i)
c
        contr = 0  !
c
        call find_match_hashtb0(aa,len,
     &         ntbrec,line_tb,
     &         ntabl_ihash,hashtblink,nhashmx,linematch,contr,find)
c
        atomMass(i) = 99.9
        if(find)then
        atomMass(i) = ams_tb(linematch)
c
        if(CONTROL1)then
        write(kanalp,'(i6,1x,a2,1x,2f8.3)')
     &   i, aa, atomMass(i), ams_tb(linematch)
        end if !control
c
        else 
        if(CONTROL)then
        write(kanalp,*)'WARNING!:getMass: atomMASS does NOT FOUND!!'
        end if
        end if !find

        end do!i
c END read AtomMass
c
c READ atom VDW parameters        
        rewind 11
        recBlock = '$VDWPAR  '
        nrb = 7
        endBlock = '$END'
        myblock = .false.
        ntbrec = 0
        nVDWtype = 0
c
200     read(11,'( a80 )', end=201 ) line
c
        if(CONTROL1)then
        write(kanalp,*)line
        end if
        if(myblock .and. (line(1:4) .eq. endBlock) ) goto 201
        if(line(1:nrb) .eq. recBlock(1:nrb) )then
         myblock = .true.  
         goto 200
         end if 
c
        if(.not.myblock)goto 200
c skip if comment characters ! or #
        if( line(1:1) .eq. '!' .or. line(1:1) .eq. '#')then
        goto 200
        end if
c
        if(myblock)then
        ntbrec=ntbrec+1
        if(ntbrec .gt. ntbmax) then
        write(kanalp, *)'ERROR!in getVDWatMass: ',
     &  ' ntbmax (param) is SMALL '
        stop
        end if
c 
        read(line, '(i2,2x,a2,10x,2f8.4)')
     &  atmType_tb(ntbrec),atName_tb(ntbrec),
     &  rm_tb(ntbrec),ev_tb(ntbrec)

c define MAX number of atomVWDtypes
        if(nVDWtype .le. atmType_tb(ntbrec))then
        nVDWtype = atmType_tb(ntbrec)
        if(nVDWtype .gt. nVDWtypeMX)then
        write(kanalp,*)
     & 'ERROR!:getVDWatMass: nVDVtypeMX is small',nVDWtypeMX
        stop
        end if
        end if
c
        if(CONTROL1)then
	write(kanalp,'(i2,2x,a2,1x,2f8.4,i6)')
     &  atmType_tb(ntbrec),atName_tb(ntbrec),
     &  rm_tb(ntbrec),ev_tb(ntbrec),ntbrec
        end if
c
        end if !myblock
c
        goto 200
c
201     continue
        close(unit=11)
c end of readIn ffParFile
c
c calculate hashTable for dictionary table
       len = tlen
       nstr1 = len
       contr = 0     !control print
c init hashtable
       do i = 1,nhashmx
       hashtblink(i) = 0
       ntabl_ihash(i) = 0
       end do!i
c 
       do i = 1,ntbrec
       line_tb(i) = atName_tb(i)   ! string length
c calculate hash number from string strtabl1(1:nstr1) simbols
       call hashCalculator(atName_tb(i),nstr1,nhashmx,nhash,contr)
 
c hashtblink() - linkedList to resolve degenerate nhash
c ntabl_ihash() - head() of linked list for given nhash
 
       hashtblink(i) = ntabl_ihash(nhash)
       ntabl_ihash(nhash) = i
 
       if(CONTROL1)then
       write(kanalp,'(a27,i6,i6,a5,a5)')
     & 'hash_dict_tabl:itab,nhash:',i,nhash,' str:',
     &      atName_tb(i)
       end if

       end do !i
c
        if(CONTROL1)then
        write(kanalp,*)'getVDWatMass: result:nVDWtype:',nVDWtype
        write(kanalp,*)'ia name VDWtype '
        end if
c
c assign VDW atomType to all atoms
        do i = 1, natom
        i2 = 2*i
        aa = ffAtomName(i)
c
        contr = 0  !
c
        call find_match_hashtb0(aa,len,
     &         ntbrec,line_tb,
     &         ntabl_ihash,hashtblink,nhashmx,linematch,contr,find)
c
        atomVDWtype(i) = 0
        if(find)then
        atomVDWtype(i) = atmType_tb(linematch)
c
        if(CONTROL1)then
        write(kanalp,'(i6,1x,a2,1x,i2)')
     &  i,aa,atomVDWtype(i)
        end if !control
c
        else
        if(CONTROL)then
        write(kanalp,*)'WARNING:getVDW: VWD param does NOT FOND!!'
        end if
c
        end if !find
c
        end do!i
c
c make table of VDW parameters for all pair of vdwTypes
        if(CONTROL)then
        write(kanalp,*)
     &    '   i1  i2  t1 t2 t12   rm12  ev12        Avdw12     Bvdw12:'
        end if
c
        do i1 = 1,ntbrec
        t1 = atmType_tb(i1)
        do i2 = i1,ntbrec
        t2 = atmType_tb(i2)
        call vdw12TablePos(nVDWtype,t1,t2,t12)

c read: RM/2  Emin from ffdata table
        call vdw12Par(rm_tb(i1),ev_tb(i1),rm_tb(i2),ev_tb(i2),
     &       rm12,ev12,Avdw12,Bvdw12)	
c nVDWatomParMX - number of VDW parameters stored per atomPair = 4
        nVDWatomParMX = 4        ! nVDWatomParMAX in vdw12Par.h  ! 
        p4 = t12*nVDWatomParMX
        atomVDW12ab(p4-3) = Avdw12
        atomVDW12ab(p4-2) = Bvdw12
        atomVDW12ab(p4-1) = rm12  
        atomVDW12ab(p4)   = ev12  
c atomVDW12ab() keeps pair A and B : vdwE=A/r*12 - B/R*6
c
        if(CONTROL)then
        write(kanalp,'(5i4,2f7.4,f14.3,f10.3)')
     &   i1,i2,t1,t2,t12,rm12,ev12,Avdw12,Bvdw12
        end if
	end do! i2
        end do !i1
c        
        return
	end
c
	subroutine vdw12TablePos(n,t1,t2,t12)
c calculate position in the atom12DWab table for atoms types t1,t2
        implicit none
        integer n,t1,t2,t12,i1,i2
        if(t2 .gt. t1)then
        i1=t1-1
        i2=t2
        else
        i1 = t2-1
        i2 = t1
        end if
        t12 = n*i1-i1*(i1-1)/2+i2-i1
c
        return
        end 
c
        subroutine  vdw12Par(rm1,ev1,rm2,ev2,
     &               rm12,ev12,A12,B12) 
c calculate VDW parameters for mixed pairs to form
c  A/r*12 - B/r*6
c
        implicit none
        real rm1,rm2,ev1,ev2
        real rm12,ev12,A12,B12,a
c
        rm12 = rm1+rm2
        ev12 = sqrt(ev1*ev2)
        a = rm12**6
        B12 = ev12*a
        A12 = B12*a           
        B12 = 2.0*B12   
c
        return
        end
