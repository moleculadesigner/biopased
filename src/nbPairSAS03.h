c pairList data structure for SAS03
c full PairList = allSoluteMolecule Atoms included
        integer nnbpLSAS03MAX
cx     parameter (nnbpLSAS03MAX = natomMAX*400)
        parameter (nnbpLSAS03MAX = natomMAX*800)    !
        integer nnbPairLSAS03(natomMAX)
        integer startnbPairLSAS03(natomMAX)
        integer nbpairListSAS03(nnbpLSAS03MAX)
        real nbpairListSAS03D2(nnbpLSAS03MAX)
        common/nbPSAS0301/nnbPairLSAS03,startnbPairLSAS03
        common/nbPSAS0302/nbpairListSAS03
cEND
