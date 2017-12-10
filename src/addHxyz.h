c pdb arrays noH and addedH
c noHydrogens
        character*4 atomName_nh(natomMAX)
        character*8 atomNameEx_nh(natomMAX) 
        character*8 resNameEx_nh(natomMAX)
        character*4 resName_nh(natomMAX)
        character*4 resNameRes_nh(nresMAX)
        character*1 chName_nh(natomMAX)
        character*1 chNameRes_nh(nresMAX) 
        character*6 head_nh(natomMAX)
        common/addH01/atomName_nh,atomNameEx_nh,resNameEx_nh,
     &                resName_nh,resNameRes_nh,chName_nh,
     &                head_nh,chNameRes_nh
        integer resNumb_nh(natomMAX)
        integer atomNumb_nh(natomMAX)
        integer natom_nh,nres_nh                
        integer startAtInRes_nh(nresMAX+1)
        integer stopAtInRes_nh(nresMAX+1)
        integer realAtomFlag_nh(natomMAX)
        integer atHvyNbInRes_nh(nresMAX) 
        common/addH02/resNumb_nh,atomNumb_nh,natom_nh,nres_nh,
     &                startAtInRes_nh,stopAtInRes_nh,
     &                realAtomFlag_nh,atHvyNbInRes_nh
c
        real*8  atomXYZ_nh(3,natomMAX)        
        common/addH03/ atomXYZ_nh
c
c with addedHydrogens
        character*4 atomName_h(natomMAX)
        character*4 resName_h(natomMAX)
        character*4 resNameRes_h(nresMAX)
        character*8 atomNameEx_h(natomMAX) 
        character*8 resNameEx_h(natomMAX)
        character*1 chName_h(natomMAX)
        character*1 chNameRes_h(nresMAX) 
        character*6 head_h(natomMAX)
        integer resNumb_h(natomMAX)
        integer atomNumb_h(natomMAX)
        integer natom_h, nres_h
        integer startAtInRes_h(nresMAX+1)
        integer stopAtInRes_h(nresMAX+1)                
        integer realAtomFlag_h(natomMAX)
        real*8  atomXYZ_h(3,natomMAX)        
c
