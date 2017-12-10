c neigbours 12,13,14
c bonds12,
c triplets123,
c quartets1234, quartes1234Imp
        integer npL12MAX
        parameter (npL12MAX = natomMAX*3)
        integer npL13MAX
        parameter (npL13MAX = natomMAX*4)
        integer npL14MAX
        parameter (npL14MAX = natomMAX*6)
        integer npL123MAX
        parameter (npL123MAX = npL12MAX*3)
        integer npL1234MAX
        parameter (npL1234MAX = npL12MAX*4)
        integer nImp1234MAX
        parameter (nImp1234MAX = natomMAX)
c
        integer npL12MX,npL13MX,npL14MX
        integer npL123MX,npL1234MX,nImp1234MX
c
        integer startPairL12(natomMAX),nPairL12(natomMAX)
        integer startPairL13(natomMAX),nPairL13(natomMAX)
        integer startPairL14(natomMAX),nPairL14(natomMAX)
        integer startPairL14V(natomMAX),nPairL14V(natomMAX)
        integer pair12List(npL12MAX)
        integer pair13List(npL13MAX)
        integer pair14List(npL14MAX)
        integer pair14VList(npL14MAX)  ! 14list for VWD/COUL interactions
        integer bond12List(npL12MAX*2) ! i1-i2 numbers, sorted List: h-heavy,hv-H pairs
        integer nbond12
        integer nbond12noH             ! number of 12bonds heavy-heavyAtom
        real bond12ParL(npL12MAX*2)          ! Kb,b0 bondDef parameters
        integer trip123List(npL123MAX)
        integer nTrip123
        real ang123ParL(npL123MAX)
c
        integer quar1234List(npL1234MAX)
        integer nQuar1234
        integer nTorsHarmMAX,nTorsParAng
        parameter (nTorsHarmMAX = 4)
        parameter (nTorsParAng = 4)
        real quar1234ParL(npL1234MAX*nTorsHarmMAX*nTorsParAng)
        integer quar1234nPar(npL1234MAX)    !number of torsHarmonics
c
        integer quarImp1234L(nImp1234MAX*4)
        integer nImp1234
        real impAng1234ParL(nImp1234MAX*2)
c
	common/pair1234/npL12MX,startPairL12,nPairL12,pair12List,
     &                  npL13MX,startPairL13,nPairL13,pair13List,
     &                  npL14MX,startPairL14,nPairL14,pair14List,
     &                  startPairL14V,nPairL14V,pair14VList
c
	common/defGeom12/nbond12,nbond12noH,bond12List,bond12ParL
c
	common/defGeom123/npL123MX,nTrip123,trip123List,ang123ParL
c
        common/defGeom1234/npL1234MX,nQuar1234,quar1234List,
     &                     quar1234ParL,quar1234nPar
c
	common/defGeomImp1234/nImp1234MX,nImp1234,quarImp1234L,
     &                     impAng1234ParL
c
c flag defines bond,triplet,quartets belongs to LIGand=(group of atoms)
c  = 0/1  
        integer vLigFlag                !0/1 
 	integer bond12LigFlag(npL12MAX)
        integer trip123LigFlag(npL12MAX) 
        integer qImp1234LigFlag(nImp1234MAX)
        integer quar1234LigFlag(npL12MAX)
c
	common/ligFlag01/vLigFlag,bond12LigFlag,trip123LigFlag,
     &                   qImp1234LigFlag,quar1234LigFlag
c
c end
