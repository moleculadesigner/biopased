c 
c Yuri Vorobjev, 2003
c
c SHake bonds
	subroutine shakeBonds(iShake,nbond12,nbond12noH,nbondShake,
     &                     bond12List,bond12ParL,shitMX,shtol,
     &              moveFlag,atomMass,atomXYZ0,atomXYZ1,
     &              shExitFlag)
c
c INput: iShake - flag for ShakeType =0,1,2 (noShake,shakeH,shakeAll)
c        nbond12,nbond12noH
c bond12List(k)= k1,k2; 
c bond12ParL(k)=(k0,b0k)
c moveFlag() - moving atom flag      
c atomMass() - atomic mass
c atomXYZ0() - atomXYZ at t=t
c atomXYZ1() - (unconstrained step t+dt) coordinates
c shitMX - MAX shake iterations
c shtol - accuracy to shake bonds
c OUTput:
c atomXYZ1() - (result in constrained step t+dt) coordinates
c shExitFlag - exitFlag =1 if succesfull; =0 if failure
c
c        implicit none
        include 'xyzPDBsize.h'
c
	integer iShake
        integer nbond12,nbond12noH 
        integer nbondShake
        integer bond12List(*)
        real bond12ParL(*)
        integer shitMX
        real shtol
        integer moveFlag(*)
        real atomMass(*)
        real atomXYZ0(*)
        real atomXYZ1(*)
        integer shExitFlag
c
        include "output.h"
        include "kanalUse.h"
        include "statusMessg_mDyn.h"
c local 
        integer nShBondMX
        parameter (nShBondMX = natomMAX)
        real shBondVec(5*nShBondMX) ! r12, r12^2, d120 for each bond
c
        real dx1(3),dx2(3)
        real tol2
        integer it,i,i2,i5,ib
        integer k,k1,k2,k31,k32
        integer nss
        integer itol
        integer kanalp
c        logical CONTROL
c
        kanalp = kanalRunOut
c        CONTROL = .false.
c
c        if(CONTROL)then
c        write(kanalp,*)'shakeBond Start:iShake:',iShake
c        end if
c
        nbondShake = 0
        shExitFlag = 1
        if(iShake .eq. 0 )return

c define bond list to shake
        tol2 = shtol**2
c
        if(iShake .eq. 1) nss = nbond12noH  ! shake H-anyAtom
        if(iShake .ge. 2) nss = 0           ! shake all bonds
        nbondShake = nbond12 - nss
c control size
        if( nbondShake .gt. nShBondMX)then
        write(kanalp,*) 'shakeBond : ERROR low nShBondMX param:',
     &  nShBondMX
        stop
        end if !
c
        it = 0
100     it = it + 1
        if( it .gt. shitMX ) then
        shExitFlag = 0
        write(kanalp,*)'Shake exeeds maxNumb Iterations:',it
        return
        end if
c
c        if(CONTROL)then
c        write(kanalp,*)'shakeBond: iteration:',it ,'*it*it*it*it*'
c        end if
c 
        itol = 1
c
        do i = nss+1,nbond12 
        i2 = i*2-1
        k1=bond12List(i2)            ! atom k1-k2 bond
        k2=bond12List(i2+1)
        k31 = k1*3-2
        k32 = k2*3-2
c
c        if(CONTROL)then
c        write(kanalp,*)'shakeBond: i,k1,k2:',i,k1,k2
c         write(kanalp,*)'shakeBond: kb,b0: ',
c     &                bond12ParL(i2),bond12ParL(i2+1)
c        end if
c 
        ib = i - nss
        i5 = 5*ib-4
c keep all bond vectors ia time=t to optimize shake calculations
        if(it .eq. 1) then
        shBondVec(i5) = atomXYZ0(k31) - atomXYZ0(k32)  
        shBondVec(i5+1) = atomXYZ0(k31+1) - atomXYZ0(k32+1)  
        shBondVec(i5+2) = atomXYZ0(k31+2) - atomXYZ0(k32+2)  
        shBondVec(i5+3) = shBondVec(i5)**2 + shBondVec(i5+1)**2 + 
     &                    shBondVec(i5+2)**2
c d12 = bond12ParL(i2+1)
        shBondVec(i5+4) =  bond12ParL(i2+1)**2
        end if !it = 1
c
	call shakeOneBond(it,moveFlag(k1),moveFlag(k2),
     &                    atomMass(k1),atomMass(k2), 
     &                    atomXYZ1(k31),atomXYZ1(k32),
     &                    shBondVec(i5),dx1,dx2,shExitFlag)
c
c check accuracy:
        do k=1,3
        if( dx1(k)**2 .gt. tol2) itol=0
        if( dx2(k)**2 .gt. tol2) itol=0 
        end do !k
c
         end do !i =ibond
c
        if (itol .eq. 0) goto 100     ! next iteration
c
c        if(CONTROL)then
c        write(kanalp,*)'SHAKE succesful,iter:',it,' toler:', shtol
c        end if
c
 	return
	end
c
        subroutine shakeOneBond(it,mf1,mf2,m1,m2,
     &                        x1uc,x2uc,b12v,dx1,dx2,iex)
c
	implicit none
        include "charStringSiz.h"
        include "output.h"
        include "kanalUse.h"
        include "statusMessg_mDyn.h"
c
        integer it          ! flag Number of call
        integer mf1,mf2     ! moveFlag
        integer iex         ! exit code; failure=0,Ok=1!
        real m1,m2
        real x1uc(3),x2uc(3)
        real b12v(5)        ! r12(i),r12^2, d12 bondVect
        real dx1(3),dx2(3)  ! corrections
c local
        real r(3)
        real m12,gt
        real rr0,r02,rr02,rr2
        real mu1,mu2,diffd2r2,shDet
        integer i
        integer kanalp
c        logical CONTROL
c        logical CORRECTION
c
        kanalp = kanalRunOut
c        CONTROL = .false.
c        CORRECTION = .true. 
c
c        if(CONTROL)then
c        write(kanalp,*)'shakeOneBond: start'
c        write(kanalp,*)'it,mf1,mf2: ',it,mf1,mf2
c        write(kanalp,*)'m1,m2:',m1,m2 
c        write(kanalp,*)'r12(i):',(b12v(i),i=1,3) 
c        write(kanalp,*)'r12^2,d12:b12v(5):',b12v(4),b12v(5)
c        write(kanalp,*)'init:x1uc:',x1uc
c        write(kanalp,*)'init:x2uc:',x2uc
c        end if
c
        r(1) = x1uc(1)-x2uc(1)
        r(2) = x1uc(2)-x2uc(2)
        r(3) = x1uc(3)-x2uc(3)
c
c        r0(i) = b12v(i)
c        r02 = b12v(4)=r0(1)*r0(1)+r0(2)*r0(2)+r0(3)*r0(3)
c        d12 = b12v(5)
c
        rr0 = r(1)*b12v(1)+r(2)*b12v(2)+r(3)*b12v(3)
        rr2 = r(1)**2+r(2)**2+r(3)**2   
c
c        if(CONTROL)then
c        write(kanalp,*)'ruc2:',rr2,' dbErr:', (rr2-b12v(5))
c        end if
c
        diffd2r2 = b12v(5) - rr2
c
        if(it .eq. 1) then
c check existance of the SHAKE solution: max shakeShift for iter=1
        shDet = rr0**2/b12v(4) + diffd2r2
c
        if(shDet .lt. 0.0 ) then
        write(kanalp,*) 'SHAKE failure ... MAXshift exeeded ...'
        write(kanalp,*) 'FATAL ERROR ',
     &   'equilibrate, decrease Temp, decrease integr dt time..!'
c
         write(kanalPStat,*)mError,
     &   'SHAKE failure ... MAXshift exeeded ...'
        stop
c       
        iex = 0
        return
        end if !
        end if !it=1
c
c do shake
        iex = 1
        m12 = m1+m2
        gt = 0.5*diffd2r2/rr0/m12
c shake corrections
        mu1 = gt*m2
        mu2 = -gt*m1
c
        do i=1,3
        dx1(i) = mu1*b12v(i) 
        dx2(i) = mu2*b12v(i) 
        end do!i
c
c do CORRECTION
        do i=1,3
        if(mf1 .ne. 0 .and. mf2 .ne. 0)then
        x1uc(i) = x1uc(i) + dx1(i) 
        x2uc(i) = x2uc(i) + dx2(i) 
        end if
c       
        if(mf1 .ne. 0 .and. mf2 .eq. 0)then
        x1uc(i) = x1uc(i) + dx1(i) - dx2(i)
        end if
c
        if(mf2 .ne. 0 .and. mf1 .eq. 0)then
        x2uc(i) = x2uc(i) + dx2(i) - dx1(i)
        end if
        end do !i
c
c        if(CONTROL)then
c        write(kanalp,*)'shakeRes: dx1:',dx1
c        write(kanalp,*)'shakeRes: dx2:',dx2
c        write(kanalp,*)'shakeRes:x1uc:',x1uc
c        write(kanalp,*)'shakeRes:x2uc:',x2uc 
c
c        r(1) = x1uc(1)-x2uc(1)
c        r(2) = x1uc(2)-x2uc(2)
c        r(3) = x1uc(3)-x2uc(3)
c        rr2 = r(1)**2+r(2)**2+r(3)**2   
c        write(kanalp,*)'ruc2:corrected:',rr2,
c     &  ' dbErr:',(rr2-b12v(5))
c        end if
c
c        if(CONTROL)then
c        write(kanalp,*)'shakeOneBond: finish'
c        end if
c
        return
        end 
