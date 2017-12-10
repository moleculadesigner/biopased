c size of data structure for IonAccessibleSurface calculation
c
        integer ndotIasMAX
        parameter (ndotIasMAX = natomMAX*200)
c
        integer ndotAtIasMAX
c        parameter (ndotAtIasMAX = 1000)  ! dots ~ 1A^2 for atom=9 A radius
        parameter (ndotAtIasMAX = 1200)   !   ~ dotdens = 8
c
