c correct forces to zero sum
c
	subroutine corr_force_to_zeroSum(ctype,nnat,ff)
c
c ctype = 0/1 - equal, proportional method
c 
	implicit none
cx        include'par-bam-sims.h'        
c
	integer ctype
	integer nnat
	real*4 ff(*)
c local
	real*8 sf(3),df(3),ssf(3),dfPM(3)
	real*8 scalP(3),scalM(3)
	real*8 sfP(3),sfM(3),fk,mf
	real*8 zero,one
	integer i,k,nnZf,i3,k1
	integer kanalpl
	logical CONTROL
c
        kanalpl= 6 !  kanalp
	CONTROL = .false.
c
        one=1.0d00
        zero=0.0d00
c
        do k=1,3
	sf(k) = zero
	df(k) = zero
        scalP(k) = zero
	scalM(k) = zero
	ssf(k) = zero
	dfPM(k) = zero
	end do!k
c
        do i=1,nnat
	i3=3*i-3
c
        mf = zero
	do k=1,3
	fk = ff(i3+k)
	mf=mf + fk**2
	sf(k) = sf(k) + fk
        if(fk .gt. zero)sfP(k)=sfP(k)+fk
        if(fk .lt. zero)sfM(k)=sfM(k)+fk
	end do!k
	if(mf .gt. zero)nnZf=nnZf+1  !nonZero forces
c
        end do!i
c
c correct forces: proportional method ctype=1
c
        do k=1,3
	 ssf(k) = zero
	 scalP(k) = zero
	 scalM(k) = zero
	 df(k) = zero
	 if(sfP(k) .gt. zero)scalP(k)=sf(k)*0.5/sfP(k)
         if(sfM(k) .lt. zero)scalM(k)=sf(k)*0.5/sfM(k)
         if(nnZf .ge. 1)df(k) = sf(k)/nnZf        
	end do !k
c
        if(CONTROL)then
	write(kanalpl,*)'corr_force_to_zeroSum:ia,ffC,df,fSum:'
	end if
c
        do i=1,nnat
	i3=3*i-3
	do k=1,3
c
        dfPM(k) = zero
	if(ctype .eq. 1)then
        if(ff(i3+k) .gt. zero)then
        ff(i3+k) = ff(i3+k)*(one - scalP(k))
	dfPM(k) = -ff(i3+k)*scalP(k)
	end if !P
c
        if(ff(i3+k) .lt. zero)then
        ff(i3+k) = ff(i3+k)*(one - scalM(k))                      
	dfPM(k) = -ff(i3+k)*scalM(k)
	end if !M
c
        end if ! ctype .eq. 1
c
	if(ctype .eq. 0)ff(i3+k) = ff(i3+k) - df(k)
c
        end do!k
c
        if(CONTROL)then
        do k=1,3
        ssf(k) = ssf(k) + ff(i3+k)
	end do !k
c
	write(kanalpl,'(a4,i5,1x,3f8.3,3f8.3,3f8.3)')
     &  'ssf:', i,(ff(i3+k1),k1=1,3),(dfPM(k1), k1=1,3),
     &    (ssf(k1),k1=1,3)
	end if!C
c
	end do!i
c
	return
	end
