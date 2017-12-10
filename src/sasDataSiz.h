c size of data structure for SolventAccessibleSurface and IonAS calculation
c
        integer nprobeMAX
        parameter (nprobeMAX=1)
c
        integer ndotMAX
        parameter (ndotMAX = natomMAX*200)
c
        integer ndotAtMAX
c        parameter (ndotAtMAX = 1000)  ! dots ~ 1A^2 for atom=9 A radius
        parameter (ndotAtMAX = 1200)   !  ~ dotdens = 8
c
