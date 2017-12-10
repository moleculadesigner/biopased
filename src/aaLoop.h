c aacidLoop info for ONE Loop
c
	integer nresLoopMAX
        parameter (nresLoopMAX = 12)
        integer natLoopMAX
        parameter (natLoopMAX = 300)
        integer nresLoop
        character*4 resNameLoop(nresLoopMAX+2)
        real FiPsiLoop(2*(nresLoopMAX+2))
        integer natLoop
        real atLoopXYZ(3*natLoopMAX) 
        real atQLoop(natLoopMAX)
        character*4 atNameLoop(natLoopMAX)
        character*2 ffatNameLoop(natLoopMAX)
        integer startAtResLoop(nresLoopMAX)
c
        integer nAnchAtLoop1
        integer anchAtListLoop1(natomMAX)
c
        integer atLoopCMatrx(5*natLoopMAX)   !ConnectivityMatrix 
c
c z-matrix for the loop
c
        character*4 LoopZmatrx1(3*natomMAX)
c LoopZmatrx1():atname,atChemName,BSname
        integer     LoopZmatrx2(4*natomMAX)
c  LoopZmatrx2():4,3,2,1 - atoms 4-3-2-1
        real        LoopZmatrx3(3*natomMAX)
c LoopZmatrx3():d34,ang234,fi1234
        common/zmatrix01L/LoopZmatrx1
        common/zmatrxi02L/LoopZmatrx2
        common/zmatrix03L/LoopZmatrx3
c                                                j2
c  fixed 5position per atom: i,j1,j2,j3,j4 : j1- i -j3 : j1,..,j4 - vbonded atoms
c                                                j4
        integer pair12ListLoop(5*natLoopMAX) !standart sequential 12 pairList
        integer start12PairLLoop(natLoopMAX) 
        integer nPair12LLoop(natLoopMAX)
c
cEND
