c arays to keep restrainDistA2 Info
        character*(charLenMAX) restrA2File
	integer nRestrDistA2MAX
        parameter (nrestrDistA2MAX = natomMAX*4)
        integer nRestrDistA2MX
c
c restrains of type2 = twoAtom restarins
        integer nRestrDistA2  ! number of restrained distancesA2    
        integer restrDistA2List(2*nRestrDistA2MAX)
        real restrDistA2DistHK(2*nRestrDistA2MAX)
c
        common/restrDistA20/restrA2File
	common/restrDistA21/nRestrDistA2,restrDistA2List,
     &                      restrDistA2DistHK  
c
cend
