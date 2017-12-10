c find_match_hashtb
       subroutine find_match_hashtb2(atm,res,rnum,chn,
     &         ntbrec,atnam_tb,term_tb,rnam_tb,rnumb_tb,chnam_tb,
     &         ntabl_ihash,hashtblink,nhashmx,linematch,contr,find)
c
c	find match of atom extendedRecord: atmName,resName,rnum,chn
c       note that resName[char*8]type: resName(1:4)=resNameRegular,
c                                      resName(5:8)=terminalTag(N or C)
c       with entry in the table: atm_tb,term_tb,res_tb,rnum_tb,chn_tb 
c       on the WORD by WORD basis
c
c       ntb - number of lines in table
c       find = .true. if match is found
c       linematch - line number in the table which match with 
c       match is found if 
c         tableline(1:lenght_tbstring)=atmstring(1:length_tbstring)
        
	implicit none
	integer ntb
	integer contr
	character*4 atm
	character*8 res
	character*1 chn
	character*4 rnum
	integer ntbrec
	character*4 atnam_tb(*)
	character*4 term_tb(*)
	character*4 rnam_tb(*)
	character*4 rnumb_tb(*)
	character*1 chnam_tb(*)
	integer ntabl_ihash(*)
	integer hashtblink(*)
	integer nhashmx
	integer linematch
	logical find
c local
	character*24 ch24,atm_st0,atm_st
	character*24 atm_st_set(6)
	integer nats_set(6)
	integer nsetmax,nset,nats 
	character*24 tb_st0, tb_st
	integer kanalpl,ns,i
	logical CONTROL

	kanalpl = 6
	CONTROL = .false.
 	if(contr.eq.1)CONTROL = .true.
c ----------------------------------------------------------------
c atom extended-name string:
	nsetmax = 6
	ch24 = '                       '
	atm_st0 = ch24
	tb_st0  = ch24
	do i = 1,nsetmax
	atm_st_set(i) = ch24
	end do!i
c make extName to match TabNames
	atm_st_set(nsetmax) = atm
	atm_st_set(5) = atm//res(5:8)
	atm_st_set(4) = atm//'    '//res(1:4)
	atm_st_set(3) = atm//res(5:8)//res(1:4)
	atm_st_set(2) = atm//'    '//res(1:4)//rnum
	atm_st_set(1) = atm//'    '//res(1:4)//rnum//chn
	nats_set(1) = 17
	nats_set(2) = 16
	nats_set(3) = 12
	nats_set(4) = 12
	nats_set(5) = 8 
	nats_set(6) = 4

c notice ! residue extended name is split in two parts and their order is chabged
c the res(5:8) part presents N,C-TERMINALtag for residue 

c initialize
	find = .false.
	linematch = 0
	nset = 0

c compare substring of words (nset) in two strings
100     nset = nset + 1
        
	atm_st0 = atm_st_set(nset)
	call clean_spacelf(atm_st0,nats_set(nset),atm_st,nats)
	call upcase(atm_st,nats)
ctest:
	if(CONTROL)then
	write(kanalpl,*)'atm_st0:',atm_st0,' nset:',nset
	write(kanalpl,*)'atm_st :',atm_st, 'nats:', nats
        end if

	call find_match_hashtb1(atm_st,nats,
     &         ntbrec,atnam_tb,term_tb,rnam_tb,rnumb_tb,chnam_tb,
     &         ntabl_ihash,hashtblink,nhashmx,linematch,contr,find)

	
	if(.not.find.and.nset.lt.nsetmax)then
	if(CONTROL)then
	write(kanalpl,*)'GOTO 100'
	end if
	goto 100
	end if
     
	return
	end
c
	subroutine find_match_hashtb1(text,len,
     &         ntbrec,atnam_tb,term_tb,rnam_tb,rnumb_tb,chnam_tb,
     &         ntabl_ihash,hashtblink,nhashmx,linematch,contr,find)
c
c find linematch in dictionaryTable which match with text(string)
c
	implicit none
	character*(*) text
	integer len
	integer ntbrec
	character*4 atnam_tb(*)
	character*4 term_tb(*)
	character*4 rnam_tb(*)
	character*4 rnumb_tb(*)
	character*1 chnam_tb(*)
	integer nhashmx
	integer linematch
	integer contr
	integer ntabl_ihash(*)
	integer hashtblink(*)
	logical find
c local
	character*20 strtabl0,strtabl1
	integer nhash
	integer nstr0,nstr1
	integer n,i,kanalpl
	logical CONTROL

	CONTROL = .false.
	if(contr.eq.1)CONTROL = .true.

	kanalpl = 6
	nstr0 = 17
	find = .false.
	linematch = 0

         call hashCalculator(text,len,nhashmx,nhash,contr)
	 if(CONTROL)then
	 write(kanalpl,*)
	 write(kanalpl,*)'find_match_hashtb1: inText:',text(1:len)
	 end if

c to resolve hashCalculator degeneracy(virozhdennost')
c make EXPLICIT comparision of text with all tableLine strings 
c linked by linked list hashtblink()

	 n = ntabl_ihash(nhash)

c check equality of TblString(n) with text()
100     continue        
	 
	 if(CONTROL)then
	 write(kanalpl,*)'find_match_hashtb1:Check datTbLine:',n
	 end if

	 if(n.eq.0)then
c	 write(kanalpl,*)'ERROR!:find_match_hashtb1:noMatchwith Table :'
c	 write(kanalpl,*)'For atomString: ',text
	 return
	 end if

	strtabl0 = atnam_tb(n)//term_tb(n)//rnam_tb(n)
     &             //rnumb_tb(n)//chnam_tb(n)
	call clean_spacelf(strtabl0,nstr0,strtabl1,nstr1)
	call upcase(strtabl1,nstr1)

	if(CONTROL)then
        write(kanalpl,*)'find_match_hashtb1:CheckTabLine=',n,
     &  ' tabLineTxt:',strtabl1(1:nstr1),' nstr1=',nstr1
	end if

	if(text(1:len) .eq. strtabl1(1:nstr1))then
	 linematch = n
	 find = .true.
	 if(CONTROL)write(kanalpl,*)'find_match_hashtb1:FIND=True'
	 return

	   else
	   n = hashtblink(n) 
	   goto 100
	end if

	return
	end
c
	subroutine find_match_hashtb0(text,len,
     &         ntbrec,line_tb,
     &         ntabl_ihash,hashtblink,nhashmx,linematch,contr,find)
c
c find linematch in dictionaryTable which match with text(string)
c
	implicit none
	character*(*) text
	integer len
	integer ntbrec
	character*20 line_tb(*)
	integer nhashmx
	integer linematch
	integer contr
	integer ntabl_ihash(*)
	integer hashtblink(*)
	logical find
c local
	character*20 strtabl1
	integer nhash
	integer nstr1
	integer n,i,kanalpl
	logical CONTROL

	CONTROL = .false.
	if(contr.eq.1)CONTROL = .true.

	kanalpl = 6
	find = .false.
	linematch = 0
        nstr1 = len
        if (len .gt. 20 )then
         nstr1 = 20 !max length of textline
         len = nstr1
         end if
c
         call hashCalculator(text,len,nhashmx,nhash,contr)
	 if(CONTROL)then
	 write(kanalpl,*)
	 write(kanalpl,*)'find_match_hashtb0:W1: inText:',text(1:len)
     &   ,' nhash: ',nhash
	 end if

c to resolve hashCalculator degeneracy(virozhdennost')
c make EXPLICIT comparision of text with all tableLine strings 
c linked by linked list hashtblink()

	 n = ntabl_ihash(nhash)

c check equality of TblString(n) with text()
100     continue        
	 
	 if(CONTROL)then
	 write(kanalpl,*)'find_match_hashtb0:W2:Check datTbLine:',n
	 end if
c
        if(n .eq. 0 ) return
c
	strtabl1(1:nstr1) = line_tb(n)(1:nstr1)
	call upcase(strtabl1,nstr1)

	if(CONTROL)then
        write(kanalpl,*)'find_match_hashtb0:W3:CheckTabLine:',n,
     &  'tabLineTxt:',strtabl1(1:nstr1),' nstr1:',nstr1
	end if

	if(text(1:len) .eq. strtabl1(1:nstr1))then
	 linematch = n
	 find = .true.
c
	 if(CONTROL)then
         write(kanalpl,*)
     &    'find_match_hashtb0:W4:FIND=True,linematch:',linematch
         end if !CONTROL
c
	 return
c
	   else
	   n = hashtblink(n) 
	   goto 100
	end if

	return
	end
c
	subroutine find_match_hashtb01(text,len,
     &         ntbrec,line_tb,
     &         ntabl_ihash,hashtblink,nhashmx,
     &         nlinematch,nlinematchMX,linematch,contr,find)
c
c find all nlinematch lines in dictionaryTable which match with text(string)
c           linematch(n),n=1,nlinematch get lineNumbes which match with text
c
	implicit none
	character*(*) text
	integer len
	integer ntbrec
	character*20 line_tb(*)
	integer nhashmx
	integer linematch(*)
        integer nlinematchMX,nlinematch
	integer contr
	integer ntabl_ihash(*)
	integer hashtblink(*)
	logical find
c local
	character*20 strtabl1
	integer nhash
	integer nstr1
	integer n,i,kanalpl
	logical CONTROL

	CONTROL = .false.
	if(contr.eq.1)CONTROL = .true.

	kanalpl = 6
	find = .false.
        nstr1 = len
        if (len .gt. 20 )then
         nstr1 = 20 !max length of textline
         len = nstr1
         end if
c
         call hashCalculator(text,len,nhashmx,nhash,contr)
	 if(CONTROL)then
	 write(kanalpl,*)
	 write(kanalpl,*)'find_match_hashtb01:W1:inText:',text(1:len)        
     &   ,' nhash: ',nhash
	 end if

c to resolve hashCalculator degeneracy(virozhdennost')
c make EXPLICIT comparision of text with all tableLine strings 
c linked by linked list hashtblink()

	 n = ntabl_ihash(nhash)

c check equality of TblString(n) with text()
100     continue        
	 
	 if(CONTROL)then
	 write(kanalpl,*)'find_match_hashtb01:W2:Check datTbLine:',n
	 end if
c
        if(n .eq. 0 ) return
c
	strtabl1(1:nstr1) = line_tb(n)(1:nstr1)
	call upcase(strtabl1,nstr1)

	if(CONTROL)then
        write(kanalpl,*)'find_match_hashtb01:W3:CheckTabLine:',n,
     &  'tabLineTxt:',strtabl1(1:nstr1),' nstr1:',nstr1
	end if !C
c
	if(text(1:len) .eq. strtabl1(1:nstr1))then
        nlinematch = nlinematch + 1
c
        if(nlinematch .gt. nlinematchMX)then
        write(kanalpl,*)'ERROR:find_match_hashtb01:nlinematch too large'        
        stop
        end if ! nlinematch .gt. nlinematchMX
c
	linematch(nlinematch) = n
	find = .true.
c
	 if(CONTROL)then
         write(kanalpl,*)
     &    'find_match_hashtb01:W4:FIND=True,linematch:',n
         end if !CONTROL
c
         end if !text(1:len) .eq. strtabl1(1:nstr1
c
         n = hashtblink(n) 
c
        goto 100
cx	end if !text(1:len) .eq. strtabl1(1:nstr1
c
	return
	end
c
c hashing of distionary table
        subroutine hash_dict_tabl
     &    (ntbrec,atnam_tb,term_tb,rnam_tb,rnumb_tb,chnam_tb,
     &     ntabl_ihash,hashtblink,nhashmx,contr)
c
c ntbrec : numbers of record in dictionry tabl
c atnam_tb(),term_tb(),rnam_tb(),rnumb_tb(),chnam_tb() : table records
c nhashmx : max size of table = 2047  - prostoe chislo
	implicit none
	integer ntbrec
	character*4 atnam_tb(*)
	character*4 term_tb(*)
	character*4 rnam_tb(*)
	character*4 rnumb_tb(*)
	character*1 chnam_tb(*)
	integer nhashmx
	integer contr
	integer ntabl_ihash(*)
	integer hashtblink(*)
c local
	integer nhash
	character*20 strtabl0,strtabl1
	integer nstr0,nstr1
	integer i,kanalpl
	logical CONTROL

	CONTROL = .false.
	if(contr.eq.1)CONTROL = .true.

	kanalpl = 6
c
	if(CONTROL)write(kanalpl,*)'start hash_dict_tabl: contr=',contr
c initiate ntabl_ihash()
	do i =1,nhashmx
	ntabl_ihash(i) = 0
	hashtblink(i) = 0
	end do

	do i = 1, ntbrec
        strtabl0 = atnam_tb(i)//term_tb(i)//rnam_tb(i)
     &             //rnumb_tb(i)//chnam_tb(i)
	nstr0 = 17


       if(CONTROL)then
       write(kanalpl,'(a60)')
     &   'atnam_tb(i),term_tb(i),rnam_tb(i),rnumb_tb(i),chnam_tb(i)'
       write(kanalpl, '(a4,a4,a4,a4,a1)')
     &   atnam_tb(i),term_tb(i),rnam_tb(i),rnumb_tb(i),chnam_tb(i)

       write(kanalpl,'(a27,i6,i6,a7,a21)')
     & 'hash_dict_tabl:itab,nhash:',i,77777,' strIn:',
     &      strtabl0
       end if

	call clean_spacelf(strtabl0,nstr0,strtabl1,nstr1)
	call upcase(strtabl1,nstr1)

c calculate hash number from string strtabl1(1:nstr1) simbols
       call hashCalculator(strtabl1,nstr1,nhashmx,nhash,contr) 
      
c hashtblink() - linkedList to resolve degenerate nhash
c ntabl_ihash() - head() of linked list for given nhash

       hashtblink(i) = ntabl_ihash(nhash)
       ntabl_ihash(nhash) = i
      
       if(CONTROL)then
       write(kanalpl,'(a27,i6,i6,a5,a21,a21)')
     & 'hash_dict_tabl:itab,nhash:',i,nhash,' str:',
     &      strtabl1,strtabl0
       end if

       end do !i
 
       return
       end
c
      subroutine hashCalculator(txt,len,nhashmx,nhash,contr)
c
c calculates hash N for string of len characters 
c
	implicit none
        character*(*) txt
	integer len,contr
	integer nhash,nhashmx
c local
	logical CONTROL
	integer kanalpl
	integer i,j,n
	integer ibase 
	character*38 string
	data string /' *0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ'/

	ibase = 3
        CONTROL = .false.
	if(contr.eq.1)CONTROL = .true.
	kanalpl = 6

	call upcase(txt,len)

        n = 1
        do i = 1,len

	j = index(string,txt(i:i))
	n = n*ibase + j
	n = abs(n)

c	if(CONTROL)then
c	write(kanalpl,*)'in:hashCalc: len:',len,'txt:',txt(1:len),
c     &  'i:',i,'txt(i:i):', txt(i:i), ' j:',j,' n:', n
c	end if
        end do

        nhash = mod(n,nhashmx) + 1

	if(CONTROL)then
	write(kanalpl,*)'in:hashCalcul: len:',len,' txt:',txt(1:len),
     &                   ' nhash:', nhash
	end if

	return
	end
c
       subroutine find_match_word(atm,res,rnum,chn,
     &  ntb,atm_tb,term_tb,res_tb,rnum_tb,chn_tb,linematch,contr,find)
c
c	find match of atom extendedRecord: atmName,resName,rnum,chn
c       note that resName[char*8]type: resName(1:4)=resNameRegular,
c                                      resName(5:8)=terminalTag(N or C)
c       with entry in the table: atm_tb,term_tb,res_tb,rnum_tb,chn_tb 
c       on the WORD by WORD basis
c
c       ntb - number of lines in table
c       find = .true. if match is found
c       linematch - line number in the table which match with 
c       match is found if 
c         tableline(1:lenght_tbstring)=atmstring(1:length_tbstring)
        
	implicit none
	integer ntb
	integer contr
	character*4 atm,atm_tb(ntb)
	character*4 term_tb(ntb)
	character*4 res_tb(ntb)
c!! char*4->char*8
	character*8 res
	character*1 chn, chn_tb(ntb)
c!!	integer rnum
	character*4 rnum
	character*4 rnum_tb(ntb)
	character*4 rnumch4
        logical find
	integer linematch
c local
	integer ntbmaxl
	parameter (ntbmaxl=2047)
        character*24 tb_stlist(ntbmaxl)
	integer ntb_lenght(ntbmaxl)
	integer ntb_lenght_list(ntbmaxl)
	integer numb_tbline(ntbmaxl)
	integer numb_tbline_list(ntbmaxl)
	integer i,nats,ntbs,ntbl,n0,im,ntb0
	integer nlist,nlist1, npos
	integer nmatch
	character*24 ch24,atm_st0,atm_st
	character*24 atm_st_set(6)
	integer nsetmax,nset,lword
	character*24 tb_st0, tb_st
	integer kanalpl,ns
	logical CONTROL

	kanalpl = 6
	CONTROL = .false.
	if(contr.eq.1)CONTROL = .true.
c ----------------------------------------------------------------
c atom extended-name string:
	nsetmax = 6
	lword = 4
	rnumch4 = rnum
	ch24 = '                       '
	atm_st0 = ch24
	tb_st0  = ch24
	do i = 1,nsetmax
	atm_st_set(i) = ch24
	end do!i

	atm_st_set(nsetmax) = atm
	atm_st_set(5) = atm//res(5:8)
	atm_st_set(4) = atm//'    '//res(1:4)
	atm_st_set(3) = atm//res(5:8)//res(1:4)
	atm_st_set(2) = atm//'    '//res(1:4)//rnumch4
	atm_st_set(1) = atm//'    '//res(1:4)//rnumch4//chn

c notice ! residue extended name is split in two part and changed the order
c the res(5:8) part presents N,C-TERMINALtag for residue 

c initialize
	find = .false.
	linematch = 0
	nlist=ntb
	nset = 0
        nmatch = 0

c compare substring of words (nset) in two strings
100     nset = nset + 1
        
	atm_st0 = atm_st_set(nset)
	n0 = lword*nsetmax
	call clean_spacelf(atm_st0,n0,atm_st,nats)
ctest:
	if(CONTROL)then
	write(kanalpl,*)'atm_st0:',atm_st0,' nset:',nset
	write(kanalpl,*)'atm_st :',atm_st, 'nats:', nats
        end if

c compare 
	nlist1 = 0
	do i = 1,nlist

c table string:
        tb_st0 = atm_tb(i)//term_tb(i)//res_tb(i)
     &                  //rnum_tb(i)//chn_tb(i)

	ntb0 = lword*nsetmax
	call clean_spacelf(tb_st0,ntb0,tb_st,ntbs)

ctest:
	if(CONTROL)then
	write(kanalpl,*)'atom String:',atm_st(1:nats),' nats:',nats,
     &          	' nset:',nset
	write(kanalpl,*)'tableString:',tb_st(1:ntbs),' ntbs:',ntbs
        end if

	if(atm_st(1:nats) .eq. tb_st(1:ntbs) ) then
	find = .true.
	 nlist1 = nlist1 + 1
	 numb_tbline(nlist1) = i
	 goto 999
        end if

	end do ! i
ctest
	if(CONTROL)then
	write(kanalpl,*)'END loop i: nlist1, nset:',nlist1, nset
	end if

C vibrat' vse tb_strings kotorie javlayutsja substring dlja atm_str ,
C iz nih vzjat' tb_string s MAX dlinnoi
	
	if(.not.find.and.nset.lt.nsetmax)then
	if(CONTROL)then
	write(kanalpl,*)'GOTO 100'
	end if
	goto 100
	end if
     
999     continue

	if(find)then
	linematch = numb_tbline(nlist1)  
        end if

c	CONTROL = .true.
	if(CONTROL)then
	atm_st_set(1) = atm//res(5:8)//res(1:4)//rnumch4//chn
	write(kanalpl,*)'END of find: atom:',atm_st_set(1),
     &	' find, linematch:',find, linematch
	if(nlist1.gt.1)then
	write(kanalpl,*)'in find:WARNING!: multiple match!:',nlist1
	end if
	end if

        return
        end
c------------------------------------------------------------------
c
	 subroutine find_match_lett(fvar,atm,res,rnum,chn,
     &        ntb,atm_tb,term_tb,res_tb,rnum_tb,chn_tb,find,linematch)
c
c	find match of atom extendedRecord: atmName,resName,rnum,chn
c       with entry in the table: atm_tb,res_tb,rnum_tb,chn_tb 
c       on the LETTER by LETTER basis
c       ntb - number of lines in table
c       find = .true. if match is found
c       linematch - line number in the table which match with 
c       match is found if 
c         tableline(1:lenght_tbstring)=atmstring(1:length_tbstring)
        
	implicit none
	integer fvar,ntb
	character*4 atm,atm_tb(ntb)
	character*4 term_tb(ntb)
	character*4 res_tb(ntb)
c!! char*4->char*8
	character*8 res
	character*1 chn, chn_tb(ntb)
c!!	integer rnum
	character*4 rnum
	character*4 rnum_tb(ntb)
	character*4 rnumch4
        logical find
	integer linematch
c local
	integer ntbmaxl
	parameter (ntbmaxl=1000)
        character*13 tb_stlist(ntbmaxl)
	integer ntb_lenght(ntbmaxl)
	integer ntb_slist(ntbmaxl)
	integer numb_tbline(ntbmaxl)
	integer numb_tbline_list(ntbmaxl)
	integer i,nats,ntbs,ntbl,n0,im
	integer nlist,nlist1, npos
	integer nmatch
	character*20 atm_st0,atm_st
	character*20 tb_st0, tb_st
	integer kanalpl,ns
	logical CONTROL

	kanalpl = 6
	CONTROL = .false.
c ----------------------------------------------------------------
c atom extended-name string:
	rnumch4 = rnum
	if(fvar.eq.0)then
	atm_st0 = atm//res(5:8)//res(1:4)//rnumch4//chn
	else
	atm_st0 = atm//'    '//res(1:4)//rnumch4//chn
	end if
c notice ! residue extended name is split in two part and changed the order
c the res(5:8) part presents N,C-TERMINALtag for residue 
	n0 = 17
	call clean_spacelf(atm_st0,n0,atm_st,nats)
ctest:
	if(CONTROL)then
	write(kanalpl,*)'atm_st0:',atm_st0,' n0:',n0
	write(kanalpl,*)'atm_st :',atm_st, 'nats:', nats
        end if

c initialize
	find = .false.
	linematch = 0
	nlist=ntb
	npos = 0
        nmatch = 0
c compare position pos in two strings
100     npos = npos + 1

c compare 
	nlist1 = 0
	do i = 1,nlist

	if(npos.eq.1)then
c table string:
        tb_st0 = atm_tb(i)//term_tb(i)//res_tb(i)
     &          	//rnum_tb(i)//chn_tb(i)

	call clean_spacelf(tb_st0,n0,tb_st,ntbs)
ctest:
	if(CONTROL)then
	write(kanalpl,*)'tableString:'
	write(kanalpl,*)'tb_st0:',tb_st0, ' n0:',n0
	write(kanalpl,*)'tb_st :',tb_st, ' ntbs:',ntbs
        end if

	tb_stlist(i)   = tb_st
	ntb_lenght(i)  = ntbs
	numb_tbline(i) = i
	ntbl = i

	else
	tb_st = tb_stlist(i)
	ntbs = ntb_lenght(i)
	ntbl = numb_tbline(i) 
	end if

c fill in array of partial_match
	if(atm_st(npos:npos) .eq. tb_st(npos:npos)) then
	 nlist1 = nlist1 + 1
	 tb_stlist(nlist1) = tb_st
	 ntb_lenght(nlist1) = ntbs
	 numb_tbline(nlist1) = ntbl
        end if

c test:
         if(CONTROL)then
         write(kanalpl,*)'in find: npos:', npos
         write(kanalpl,'(a42,2i6,2x,a15,i6,i6)')
     &  'in find:i, nlist1,tb_st,ntbs,ntbl:',
     &   i,nlist1,tb_st,ntbs,ntbl
         end if

	 if(nlist1.ge.1)then
         if( npos .eq. ntb_lenght(nlist1) ) then
	   nmatch = nmatch + 1
	   numb_tbline_list(nmatch) = numb_tbline(nlist1)
	   ntb_slist(nmatch) = numb_tbline(nlist1)
	   find = .true.

           if(CONTROL)then
	   write(kanalpl,*)'in find: npos:', npos
	   write(kanalpl,*)'nmatch, ntb_slist(nmatch):',
     &     nmatch, ntb_slist(nmatch)
           end if
	 end if
	 end if

	end do ! i
ctest
	if(CONTROL)then
	write(kanalpl,*)'END loop i: nlist1, npos:',nlist1, npos
	end if

C vibrat' vse tb_strings kotorie javlayutsja substring dlja atm_str ,
C iz nih vzjat' tb_string s MAX dlinnoi
	
	if(nlist1 .ge. 1 ) then
	nlist = nlist1
	if(CONTROL)then
	write(kanalpl,*)'GOTO 100'
	end if
	goto 100
	end if
     
999     continue

c find line whith MAX simbols
	ns = 0
	if(nmatch.ge.1)then
	do i = 1,nmatch
        if(ntb_slist(i).gt.ns)then
	ns=ntb_slist(i)
	im = i
	end if
	end do
	linematch = numb_tbline_list(im)
	end if

c	CONTROL = .true.
	if(CONTROL)then
	write(kanalpl,*)"fvar:",fvar," nmatch: ", nmatch
	write(kanalpl,*)'END of find: atom:',atm_st0,
     &	' find, linematch:',find, linematch
	end if

        return
        end
c ENDfile
