c initialize :
c               bond12LigFlag()
c               trip123LigFlag()
c               qImp1234LigFlag()
c               quar1234LigFlag()
c   in data structure: pair1234array.h
c
c Yuri V. 2003
c
	subroutine initLigFlag1234
c
        include 'xyzPDBsize.h'
c
        include "output.h"
        include "kanalUse.h"
        include "statusMessg_mDyn.h"
c
        include 'pair1234array.h'
        include 'ligInfo.h'
c local
        integer i,k,i2,i3,i4
        integer j1,j2,j3,j4,iL
        integer kanalp
        logical CONTROL
c
        CONTROL = .true.
        kanalp = kanalRunOut
c
        if(CONTROL)then
        write(kanalp,*)'initLigFlag1234 :'
        end if
c
        do i=1,npL12MAX
        bond12LigFlag(i) = 0
        trip123LigFlag(i) = 0
        quar1234LigFlag(i) = 0
        end do !i
c
        do i=1,nImp1234MAX
        qImp1234LigFlag(i) = 0
        end do
c
c bond12LigFlag:
        do i = 1,nbond12
        i2 = 2*i
        j1 = bond12List(i2-1)
        j2 = bond12List(i2)
        do k=1,nAtomInLig
        iL = atomInLigList(k)
        if(iL .eq. j1 .or. iL .eq.j2) then
        bond12LigFlag(i)=1
        goto 101
        end if
        end do!k
101     continue
        if(CONTROL)then
        write(kanalp,*)'iBond12:',i,
     &  'bond12LigFlag(i):',bond12LigFlag(i),' j12:',j1,j2
        end if !C
        end do !ib 
c         
        do i = 1,nTrip123
        i3 = 3*i
        j1 = trip123List(i3-2)
        j2 = trip123List(i3-1)
        j3 = trip123List(i3)
c
        do k=1,nAtomInLig
        iL = atomInLigList(k)
        if(iL .eq. j1 .or. iL .eq.j2 .or. iL .eq. j3) then
        trip123LigFlag(i)=1
        goto 102
        end if
        end do!k
102     continue
        if(CONTROL)then
        write(kanalp,*)'iTrip123:',i,
     &  'trip123LigFlag(i):',trip123LigFlag(i),' j123:',j1,j2,j3
        end if !C
        end do !ib
c
        do i = 1,nQuar1234
        i4 = 4*i
        j1 = quar1234List(i4-3)
        j2 = quar1234List(i4-2)
        j3 = quar1234List(i4-1)
        j4 = quar1234List(i4)
c
        do k=1,nAtomInLig
        iL = atomInLigList(k)
        if(iL .eq. j1 .or. iL .eq.j2  
     &     .or. iL .eq. j3 .or. iL .eq. j4) then
        quar1234LigFlag(i) = 1
        goto 103
        end if
        end do!k
103     continue
        if(CONTROL)then
        write(kanalp,*)'iQuart:',i,
     &  'quar1234LigFlag(i):',quar1234LigFlag(i),' j1234:',j1,j2,j3,j4
        end if !C
        end do !ib
c 	
        do i = 1,nImp1234 
        i4 = 4*i
        j1 = quarImp1234L(i4-3)
        j2 = quarImp1234L(i4-2)
        j3 = quarImp1234L(i4-1)
        j4 = quarImp1234L(i4)
c
        do k=1,nAtomInLig
        iL = atomInLigList(k)
        if(iL .eq. j1 .or. iL .eq.j2
     &     .or. iL .eq. j3 .or. iL .eq. j4) then
        qImp1234LigFlag(i) = 1
        goto 104
        end if
        end do!k
104     continue
        if(CONTROL)then
        write(kanalp,*)'iImpQuart:',i,
     &  'qImp1234LigFlag(i):',qImp1234LigFlag(i),' j1234:',j1,j2,j3,j4
        end if !C
        end do !ib
c
        if(CONTROL)then
        write(kanalp,*)'END initLigFlag1234 !'
        end if
c 
	return
        end
