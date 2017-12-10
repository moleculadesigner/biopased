c xyzPDBsize.h  ! should be first in main program
c
        implicit none
c parameters for siz
        integer natomMAX
cx        parameter (natomMAX=10000)  !OK
        parameter (natomMAX=20000)
cx        parameter (natomMAX=40000)    !to big
        integer nresMAX
        parameter (nresMAX = natomMAX/5)
        integer ncovbMAX
        parameter (ncovbMAX = 6 )
c
        include "charStringSiz.h"           !defines size of fileName strings
c end
