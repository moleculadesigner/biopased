* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C*     neutralizeQ(natom,inatmq,neutratmq)                                    *
C* Defines  the neutralized atomic charges of the set                         *
C* for Lazaridis & Karplus Gaussian shell                                     *
C* solvation model PROTEINS 1999, 35, 133-152                                 *
C*                                                                            *
C*     Yury Vorobjev 2002                                                     *
C* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C neutralization of total charge Q for the set of atoms
	subroutine neutralizeQ(natom,inatmq,qt,neutratmq)

c natom - number of atonms ia=1,..,natom
c inatmq(ia) - Initial atomic q
c qt - total initial charge
c neutratmq(ia) - neutralized atomic q

	implicit none
c
        include "charStringSiz.h"
        include "output.h" 
        integer natom
        real inatmq(*)
	real neutratmq(*)
        real qt
c
        integer i
        real kp,km  
        real sqp,sqm
        real sqpc,sqmc
        integer kanalp
        logical CONTROL
c
        CONTROL = .false.
        kanalp= kanalRunOut  
c
        sqp=0.0
        sqm=0.0
        do i = 1,natom
        if(inatmq(i) .gt. 0.0 )then
        sqp = sqp + inatmq(i)
        else
        sqm = sqm + inatmq(i)
        end if
        end do !i
        qt = sqp + sqm
c
        kp=0.0
        km=0.0
        if(sqp .ne. 0.0 )then
        kp = - qt/(2.0*sqp)
        else
        kp = -1.0
        end if
        if(sqm .ne. 0.0 )then
        km = - qt/(2.0*sqm)
        else 
        km = -1.0
        end if
c control
        sqpc = 0.0
        sqmc = 0.0
        do i = 1,natom
        if(inatmq(i) .gt. 0.0 )then
        neutratmq(i) = inatmq(i)*(1.0 + kp)
        sqpc = sqpc + neutratmq(i)
        else
        neutratmq(i) = inatmq(i)*(1.0 + km) 
        sqmc = sqmc + neutratmq(i)
        end if
        end do !i
c
        if(CONTROL)then
        write(kanalp, '(a22,f6.3,a8,2f8.4,a17,2f8.4)')
     &  'neutralizeQ: qtotRes:',qt,' kp,km: ',
     &  (1.0+kp),(1.0+km),' Neutr: sQp,sQm:',sqpc,sqmc
        end if !CONTR
c
        return
        end
