c z-matrix arrays
        character*(charLenMAX) zmfile
        character*(charLenMAX) zmfileINT   ! Zmatrix Lib file INT res
        character*(charLenMAX) zmfileBEG   ! Zmatrix Lib file BEG res
        character*(charLenMAX) zmfileFIN   ! Zmatrix Lib file FIN res
        character*(charLenMAX) zmfileLIG   ! Zmatrix Lib file LIG res 
c
        integer defFlag_zmfileLIG
        common/zmatrix00/zmfileINT,zmfileBEG,zmfileFIN,
     &  zmfileLIG,defFlag_zmfileLIG
c         
        character*4 molZmatrx1(3*natomMAX)
c molZmatrx1():atname,atFFName,blockName
        integer     molZmatrx2(4*natomMAX) 
c  molZmatrx2():4,3,2,1 - atoms 4-3-2-1 
        real        molZmatrx3(3*natomMAX)
c molZmatrx3():d34,ang234,fi1234
        common/zmatrix01/molZmatrx1
        common/zmatrxi02/molZmatrx2
        common/zmatrix03/molZmatrx3
c
c number of heavy atoms of res in ZMatrix Lib
        integer atHvyNbInResZm(nresMAX)
        common/zmatrix05/atHvyNbInResZm
c CycleClosure 
        integer nAtPairCycMAX
        parameter (nAtPairCycMAX = nresMAX)
        integer nAtPairCyc
        integer atPairCycList(2*nAtPairCycMAX)
        common/zmatrix04/nAtPairCyc,atPairCycList
c
