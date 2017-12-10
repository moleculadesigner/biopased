c distribute eng P-P, P-Wat, W-W
c
	subroutine engPWatDistribute(iPWtag1,iPWtag2,ed,eTotPW)
c
c iPW1,iPW2 = 1, 2 = P(W) atomTag
c ed -delta Eng, eTot
	integer iPWtag1,iPWtag2
	real ed,eTotPW(3)
c
        if(iPWtag1.eq.1 .and. iPWtag2.eq.1)
     &    eTotPW(1)=eTotPW(1) + ed
        if((iPWtag1.eq.1 .and. iPWtag2.eq.2) .or.
     &       (iPWtag1.eq.2 .and. iPWtag2.eq.1))
     &    eTotPW(2)=eTotPW(2) + ed
        if(iPWtag1.eq.2 .and. iPWtag2.eq.2)
     &    eTotPW(3)=eTotPW(3) + ed
c
         return
	 end
c
        subroutine writePWatSh(kanalPDBres,engName,engPW)
c Write energy term for OPT_SolvateExWat=.true.
c
        integer kanalPDBres
        character*(*) engName
	real engPW(*)
c local
        write(kanalPDBres,'(a9,a5,f12.3,2(a6,f12.3))')
     &   engName,':P-P:',engPW(1),
     &   ' P-W:',engPW(2),' W-W:',engPW(3)
c
         return
	 end
c
   	subroutine scalePWatShEng(dsh,scale)
c
c empirical scaling coeff to increas Esolvation for thin hydration shell
c
	implicit none
	real dsh,scale
c local
        real scaleTab(4)
	real dshTab(4)
	data scaleTab/1.42, 1.08, 1.01, 1.00/
	data dshTab/3.5, 4.5, 6.5, 8.5/
	integer i,ii,nMX,nMN
c
        scale = 1.0
	return
c -------------------------------------------
        nMX=4
	nMN=1
        if(dsh .le. dshTab(1))then
	scale=scaleTab(1)
	return
	end if
	if(dsh .ge. dshTab(nMX))then
	scale=scaleTab(nMX)
	return
	end if
c
        if(dsh .lt. dshTab(nMX) .and. dsh .gt. dshTab(1))then
        do i=1,nMX-1
        if(dsh .le. dshTab(i+1) .and. dsh .ge. dshTab(i))then
        ii = i
        goto 101
	end if!dsh .lt. dshTab(i+1)
	end do!i
101	continue
c
       scale = ((scaleTab(ii+1)-scaleTab(ii))/
     &  (dshTab(ii+1)-dshTab(ii)))*(dsh - dshTab(ii))+scaleTab(ii)
c
        end if!dsh .lt. dshTab(nMX)
c
        return
	end
c
