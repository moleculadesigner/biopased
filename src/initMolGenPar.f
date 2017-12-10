c initMolGenPar
c
c Yurii Vorobjev 2003
c
c initialize general molecular parameters: molGenPar.h
c gyration rad : molRGyr_c 
c
        subroutine    initMolGenPar(n,atxyz,rGyr,rC)
c
        implicit none
        integer n
        real atxyz(*)
        real rC(3)
        real rGyr
c
        integer i,i3
c
        rC(1) = 0.0
        rC(2) = 0.0
        rC(3) = 0.0
        rGyr = 0.0
        do i =1,n
        i3=3*i-3
        rC(1) = rC(1) + atxyz(i3+1)
        rC(2) = rC(2) + atxyz(i3+2)
        rC(3) = rC(3) + atxyz(i3+3)
        end do
c
        rC(1) = rC(1)/n
        rC(2) = rC(2)/n
        rC(3) = rC(3)/n
c
        do i =1,n
        i3=3*i-3
        rGyr = rGyr + (atxyz(i3+1)-rC(1))**2+
     &                (atxyz(i3+2)-rC(2))**2+
     &                (atxyz(i3+3)-rC(3))**2
        end do!i
        rGyr = sqrt(1.6667*rGyr/n)
c
c        write(*,*)'initMolGenPar:nat, rGyr :',n,rGyr,' rC:',rC
c
        return
        end 
