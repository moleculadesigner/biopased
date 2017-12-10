c assign h-tag to atoms(PDBfile) to add H atoms if required: data h_add.dat file 
c
c Y.N.Vorobjev     2002
c	
        subroutine assign_htg(htgfile,
     &              natm,atmnam,resnam,resnumb,chnam)
c
c       note that resName[char*8]type: resName(1:4)=resNameRegular,
c                                      resName(5:8)=terminalTag(N or C)
c..................................................................
       implicit none
       include "assign_htg.h"
c
        include "charStringSiz.h"
        include "output.h"
        include "kanalUse.h"
        include "statusMessg_mDyn.h"
c
       character*(*) htgfile
       integer natm
       integer resnumb(*)
       character*4 atmnam(*)
       character*8 resnam(*)
       character*1 chnam(*)
c.local..................
c hashTable.h
	integer nhashmx
	parameter (nhashmx = 2047 )
	integer ntabl_ihash(nhashmx)
	integer hashtblink(nhashmx)
c
        logical CONTROL
        integer i,j,k,iun11
        integer ntbrec
	integer iiprint0,iiprint1
        character*80 liner
	character*80 comment
        character*4 atn 
	character*4 rnumch4
	integer kanalpl
	logical find
	integer fvar,fcontr
	integer linematch
	logical OPT_hash
	logical fext
c
c initialize
cx        kanalpl = 6	
        kanalpl = kanalRunOut
c
cx        CONTROL = .true. !local print
        CONTROL = .false.
	OPT_hash = .true. 
c
	if(CONTROL)then
	write(kanalpl,*)'Start:assign_htg:  '
        do i = 1,natm
        write(kanalpl,'(a6,i5,1x,a4,2x,a8,1x,a1,i6)')
     &  'Atom  ',i,atmnam(i),resnam(i),chnam(i),resnumb(i)
 	end do !i
	end if !C
c read atmhtag.dat file 
c  check file
        inquire(file=htgfile, exist=fext)
	if(.not. fext) then
	write(kanalpl,*) 'ERROR:Dat file:',htgfile,' does not exist'
c
         write(kanalPStat,*)mError,
     &   ' file ', htgfile,' does not exist'
c
	stop
	end if

	iun11 =  kanal_htg
	open(unit=iun11,file=htgfile,status='old',err=901)
	
	write(kanalpl,*)'reading hatoms info from file'
	write(kanalpl,*)htgfile
c
c skip/print comments (start with # !) and one column header line 
c......read file ........................................

        ntbrec = 0
100     continue
	 read(iun11,204, end=300)liner
       
        if(CONTROL) write(kanalpl,204)liner

         if(liner(1:3).eq."   " .or. 
     &	    liner(1:1).eq."#" .or. liner(1:1).eq."!")then 
            goto 100 
	  end if !
           if(liner(1:3) .eq."end" )then
	     goto 300
	   end if

	  ntbrec = ntbrec + 1
	  if(ntbrec.gt.ntbmax) then
	    write(kanalpl,*)'assingHtag:maximum # of records exceeded'
	    write(kanalpl,*)' increase nhtmax in assign_htg.h'
	    stop 
	  end if

	  read(liner,200)atnam_tb(ntbrec),term_tb(ntbrec),
     &	  rnam_tb(ntbrec),
     &	  rnumb_tb(ntbrec),chnam_tb(ntbrec),htag_tb(ntbrec),
     &    hname_tb(1,ntbrec),hname_tb(2,ntbrec),hname_tb(3,ntbrec),
     &    atmJKL(1,ntbrec),atmJKL(2,ntbrec),atmJKL(3,ntbrec)

c test
	if(CONTROL)then
	write(kanalpl,*)'CONTROL htgfile: line:',ntbrec
	write(kanalpl,200)atnam_tb(ntbrec),term_tb(ntbrec),
     &	  rnam_tb(ntbrec),
     &	  rnumb_tb(ntbrec),chnam_tb(ntbrec),htag_tb(ntbrec),
     &    hname_tb(1,ntbrec),hname_tb(2,ntbrec),hname_tb(3,ntbrec),
     &    atmJKL(1,ntbrec),atmJKL(2,ntbrec),atmJKL(3,ntbrec)
         end if

	 goto 100
c..................................................................
300   continue
	close(iun11)
  	write(kanalpl,*)'# htgfile:of parameter records:',ntbrec
c
c calculate hashTable for the dictionaryRecords
         fcontr = 0   ! control
c
	 call hash_dict_tabl
     &    (ntbrec,atnam_tb,term_tb,rnam_tb,rnumb_tb,chnam_tb,
     &     ntabl_ihash,hashtblink,nhashmx,fcontr)

	if(CONTROL)then
	do i = 1,nhashmx
	write(kanalpl,*)'i:',i,'  ntabl_ihash:',ntabl_ihash(i),
     &                  ' hashtblink :',hashtblink(i)
	end do
	end if

c  find record number in the table(htgfile) for the atom i
	fcontr = 0      !control
c
	do i = 1,natm

	call conv_numb_to_char(resnumb(i),rnumch4)

	if(.not.OPT_hash)then
        call find_match_word(atmnam(i),resnam(i),rnumch4,chnam(i),
     &        ntbrec,atnam_tb,term_tb,rnam_tb,rnumb_tb,chnam_tb,
     &        linematch,fcontr,find)
	else
	call find_match_hashtb2(atmnam(i),resnam(i),rnumch4,chnam(i),
     &        ntbrec,atnam_tb,term_tb,rnam_tb,rnumb_tb,chnam_tb,
     &	      ntabl_ihash,hashtblink,nhashmx,linematch,fcontr,find)
	end if !OPT_hash

	if(find)then
        atmrecn(i) = linematch 
	else
c
	if(CONTROL)then
	write(kanalpl,*)'WARNING!: at rec does not exist
     &  in the atmhtag.dat for atN:',i
	write(kanalpl,*)atmnam(i),resnam(i),resnumb(i),chnam(i)
        atmrecn(i) = 0	
	end if!CONTROL
c	
	end if !find
c
        if(CONTROL)then 
        write(kanalpl,'(a6,i5,1x,a4,2x,a8,1x,a1,2x,a14,i6,1x,L2)')
     &  'Atom  ',i,atmnam(i),resnam(i),chnam(i),
     &  ' nlineHtag_tb:',atmrecn(i),find
        end if !C
c
        end do !i
c 
	if(CONTROL)then
	write(kanalpl,*)'Result:assign_htg:  Htag assignment:'
        do i = 1,natm
        write(kanalpl,'(a6,i5,1x,a4,2x,a8,1x,a1,2x,a14,i6,1x,L2)')
     &  'Atom  ',i,atmnam(i),resnam(i),chnam(i),
     &  ' nlineHtag_tb:',atmrecn(i)
 	end do !i
	end if
c        
	write(kanalpl,*)'END assign_htg:'
c
200   format(A4,2x,A4,A4,A4,A1,i2,4x,3(A4,2x),3(A4,2x))
204   format(a80)
c
        return
c error messages
901   write(kanalpl,*)'ERRORS while readining file:', htgfile
c
	end
c ENDfile
