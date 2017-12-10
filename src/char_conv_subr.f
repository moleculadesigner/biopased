c set of SUBR for character convertion
c
c n1
      subroutine upcase(txt,len)
c
c convert character string txt to upper case
c
	implicit none
        character*(*) txt
	integer len
c
	integer match,i
	integer kanalpl
	logical CONTROL
        character*80 save
  	character*26 ualpha,lalpha
        data ualpha /'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
        data lalpha /'abcdefghijklmnopqrstuvwxyz'/

	CONTROL = .false.
clinux  CONTROL = .true.
	kanalpl = 6

	  do i=1,len
	  save(i:i) = txt(i:i)
	  if((txt(i:i).ge.'a').and.(txt(i:i).le.'z')) then
	    match = index(lalpha,txt(i:i))
            save(i:i) = ualpha(match:match)
	  end if
	  end do

	  if(CONTROL)then
	  write(kanalpl,*) 'upcase: len:',len,
     &	  ' txt:',txt, ' UPsave:',save(1:len)
	  end if

          txt = save(1:len)

	  return
	  end
c
c  delete and compres spaces from string0
	subroutine clean_space(strn0,n0,strnc,n)
c n is number of final imbols
c
	character*(*) strn0,strnc
	integer n0,n

	integer i
	integer pos(80)
	character*1 ch1

c initialize
	strnc = "                    "                            
	do i=1,n0
	pos(i) = 0
	end do
c new positions
	n = 0
	do i = 1,n0
        if(.not. (strn0(i:i).eq." " .or. strn0(i:i).eq. "*"))
     &  then	
	n = n + 1
	pos(i)=n
	end if
	end do
c 
	if(n.gt.0)then
	do i = 1,n0
	j = pos(i)
	if( j .ne. 0 )then
	strnc(j:j) = strn0(i:i)
	end if
	end do
	end if

	return
	end
c ..........................................................................
c  delete right spaces from string0
	subroutine clean_spacelf(strn0,n0,strnc,n)
clinux
	character*(*) strn0,strnc
	integer n0,n
	logical CONTROL
	integer i,nspace
	integer kanalpl
	kanalpl = 6
c initialize 
	CONTROL = .false.
	strnc = "                        "  

	if(CONTROL)then
	write(kanalpl,'(a26,a24)')'in clean_spacelf1: strIn=',strn0
	end if

	nspace = 0
	do i = n0,1,-1
	if( strn0(i:i) .eq.' ' .or. strn0(i:i) .eq.'*')then
	nspace =  nspace + 1
	else
	goto 101
	end if
	end do

101     n = n0 - nspace
	strnc(1:n) = strn0(1:n)

        if(CONTROL)then
	write(kanalpl,'(a26,i6,a24,a24)')
     &  'in clean_spacelf2: strIn=',n,strn0(1:n),strnc(1:n)
	end if

	return
	end
c .........................................................................
	subroutine conv_numb_to_char(numb,rnumbc)

	implicit none
	integer numb
	character*(*) rnumbc

	integer a,b,i,j,k,n
	character*10 cnumb
	data cnumb /'0123456789'/
c initialize
	rnumbc = '    '
        a = 10
	n = 0
	b = 1
100     b = b*a
        n = n + 1
	rnumbc(n:n) = ' '
	if(numb.ge.b)goto 100

	b = numb
	do i = 1,n
	j = mod(b,a) + 1
	k = n-i+1
	rnumbc(k:k) = cnumb(j:j)
	b = b/a
        end do

	return
	end
c .........................................................................
	subroutine conv_numb_to_char00n(numb,rnumbc4)
c
	implicit none
	integer numb
	character*(*) rnumbc4
c loc
        character*10 rnumbc
	integer a,b,i,j,k,n
        integer kanalp
	character*10 cnumb
	data cnumb /'0123456789'/
c initialize
        kanalp = 6
	rnumbc = '    '
        a = 10
	n = 0
	b = 1
100     b = b*a
        n = n + 1
	rnumbc(n:n) = ' '
	if(numb.ge.b)goto 100

	b = numb
	do i = 1,n
	j = mod(b,a) + 1
	k = n-i+1
	rnumbc(k:k) = cnumb(j:j)
	b = b/a
        end do
c
        if(numb .le. 9)then
          rnumbc4 = '000'//rnumbc(1:1)
           else
            if(numb .le. 99)then
             rnumbc4 = '00'//rnumbc(1:2) 
              else
               if(numb .le. 999)then
                rnumbc4 = '0'//rnumbc(1:3)
                 else
                  if(numb .le. 9999)then
                   rnumbc4 = rnumbc(1:4)
                    else
         write(kanalp,*)
     &    'conv_numb_to_char00n: ERROR!:number of nRecPdbGL > 9999 '
                 stop
                end if
               end if
              end if
             end if
	return
	end
c
c  delete right spaces from string0
        subroutine clean_space_right(strn0,n0,strnc,n)
clinux
        character*(*) strn0,strnc
c local
        include "charStringSiz.h"
        character*(charLenMAX) strnLoc
        integer n0,n
        logical CONTROL
        integer i,nspace
        integer kanalpl
        kanalpl = 6
c initialize
        CONTROL = .false.
        strnc = " "
c 
        if(CONTROL)then
        write(kanalpl,'(a26,a24)')'in clean_spacel_right: strIn=',strn0
        end if
c
c define number of space from right to 1 character 
        nspace = 0
        do i = n0,1,-1
        if( strn0(i:i) .eq.' ' .or. strn0(i:i) .eq.'"')then
        nspace =  nspace + 1
        else
        goto 101
        end if
        end do
 
101     n = n0 - nspace
        strnLoc(1:n) = strn0(1:n)
 
        if(CONTROL)then
        write(kanalpl,'(a26,i6,a24,a24)')
     &  'in clean_space_right1: strIn=',n,strn0(1:n),strnLoc(1:n)
        end if
c
         nspace =0
        do i = 1,n
         if(strnLoc(i:i).eq.' ' .or.
     &      strnLoc(i:i).eq.'"' )then
        nspace = nspace +1
        else
        goto 102
        end if
        end do !i
c
102     strnc(1:n-nspace)=strnLoc(nspace+1:n)
        n = n - nspace
c 
        if(CONTROL)then
        write(kanalpl,'(a26,i6,a24,a24)')
     &  'in clean_space_right2: strIn=',n,strn0(1:n0),strnc(1:n)
        end if
c
        return
        end
c
c ENDfile
