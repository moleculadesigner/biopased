c
        subroutine print9r(kanalp,na,arr)
c
        implicit none
        integer fmt
        parameter (fmt=9)
        integer nn,na
        integer arr(*)
        integer i,k,na10,nars,i10
        integer kanalp

cx       kanalp = 6

        nn=fmt
        na10=na/nn
        nars=na - na10*nn
        if(na10 .ge. 1)then
        do i = 1,na10
        i10 = i*nn - nn
        write(kanalp,'(9f8.3)')
     &   (arr(i10+k),k=1,nn)

        end do !i
        end if
c
        if( nars .gt. 0)then
        i10 = na10*nn
         write(kanalp,'(9f8.3)')
     &   (arr(i10+k),k=1,nars)
        end if
c
        return
        end
c
        subroutine print10(na,arr)
c
        implicit none
        include "charStringSiz.h"
        include "output.h"
        integer nn,na
        integer arr(*)
        integer i,k,na10,nars,i10
        integer kanalp

        kanalp = kanalRunOut

        nn=12
        na10=na/nn
        nars=na - na10*nn
        if(na10 .ge. 1)then
        do i = 1,na10
        i10 = i*nn - nn
        write(kanalp,'(12(i6,1x))')
     &   (arr(i10+k),k=1,nn)
        end do !i
        end if
c
        if( nars .gt. 0)then
        i10 = na10*nn
         write(kanalp,'(12(i6,1x))')
     &   (arr(i10+k),k=1,nars)
        end if
c
        return
        end
c
c print 10 real in line
c
        subroutine print10r(kanalp,na,arr)
c
        implicit none
        integer fmt
        parameter (fmt=10)
        integer nn,na
        real arr(*)
        integer i,k,na10,nars,i10
        integer kanalp
 
cx       kanalp = 6
 
        nn=fmt
        na10=na/nn
        nars=na - na10*nn
        if(na10 .ge. 1)then
        do i = 1,na10
        i10 = i*nn - nn
        write(kanalp,'(10f7.1)')
     &   (arr(i10+k),k=1,nn)
 
        end do !i
        end if
c
        if( nars .gt. 0)then
        i10 = na10*nn
         write(kanalp,'(10f7.1)')
     &   (arr(i10+k),k=1,nars)
        end if
c
        return
        end
