c z-matrix (temporary) arrays for missingAtoms/loopXYZ calculation
c         
        character*4 molZmatrx1Tmp(3*natomMAX)
c molZmatrx1():atname,atFFName,blockName
        integer     molZmatrx2Tmp(4*natomMAX) 
c  molZmatrx2():4,3,2,1 - atoms 4-3-2-1 
        real        molZmatrx3Tmp(3*natomMAX)
c molZmatrx3():d34,ang234,fi1234
        common/zmatrix01Tmp/molZmatrx1Tmp
        common/zmatrxi02Tmp/molZmatrx2Tmp
        common/zmatrix03Tmp/molZmatrx3Tmp
c number of heavy atoms in res in ZMatrix
        integer atHvyNbInResZmTmp(nresMAX)
        common/zmatrix05Tmp/atHvyNbInResZmTmp
c CycleClosure 
        integer nAtPairCycMAXTmp
        parameter (nAtPairCycMAXTmp = nresMAX)
        integer nAtPairCycTmp
        integer atPairCycListTmp(2*nAtPairCycMAXTmp)
        common/zmatrix04Tmp/nAtPairCycTmp,atPairCycListTmp
c
