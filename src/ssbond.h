c S-S bonds  data structure
c
        logical OPT_SSbondAuto      !all CYS in the inPDB will be replaced by CYX (involved in SS)
c                                   ! default: OPT_SSbondAuto=.true.
        logical OPT_SSBonds         !constructs S-S bonds on the set of CYX/SG residues
        integer nCYSresMAX
        parameter (nCYSresMAX = 50)    !max Numb CYS residues
        integer nCYSres,nSSbondAtom
        integer ssBondAt12List(nCYSresMAX) ! (i,j),(k,L), ... -atoms forming SS bonds
        integer ssBondAtFlag(nCYSresMAX) ! 0/1 =1 if involved
        integer nSSbonds                !total SS bonds 
        integer cysResList(nCYSresMAX) ! S-S bond cysteine residue List
        integer ssAtomList(nCYSresMAX) ! atomNumb can be involved in S-S bonds
        character*4 ssCYSname          ! CYX default in amber top Lib
        character*4 ssAtname           ! SG default
        real dss2MAX_OPT               ! max SS distant
        real ssBondDist0(nCYSresMAX)
c
        common/ssbond01/dss2MAX_OPT,ssBondDist0,
     &                  nSSbonds,nCYSres,nSSbondAtom,
     &                  ssBondAt12List,ssBondAtFlag,
     &                  cysResList,ssAtomList
        common/ssbond02/ssCYSname,ssAtname
c
        common/ssbond03/OPT_SSBonds,OPT_SSbondAuto
c  
