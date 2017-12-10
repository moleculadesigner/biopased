c global variables
c _PWat : P-P, P-W, W-W interactions
c energy 
	real eVbondDef_PWat(3)    !eng 1,2,3 = P-P, P-W, W-Wat interactions
	real eVangDef_PWat(3)
	real engVDWR1_PWat(3)
	real engCOULR1_PWat(3)
	real hBHxYeng128_PWat(3)
c 
c in range R1 - R2 mediumForce F2
	real engCOULR2_PWat(3)
	real molSolEn_PWat(3)
	real eGeoDef_PWat(3)
	real engCOUL_PWat(3)
	real engPOTENT_PWat(3)
c
        common/enforce_PW01/engVDWR1_PWat,engCOULR1_PWat,
     &  	engCOULR2_PWat,eVbondDef_PWat,eVangDef_PWat,
     &          hBHxYeng128_PWat,eGeoDef_PWat,molSolEn_PWat,
     &          engPOTENT_PWat,engCOUL_PWat
c
        integer atomPWatTag(natomMAX)
	integer natomSolutMol,natomWSolv
	common/enforce_PW02/natomSolutMol,natomWSolv,atomPWatTag
c
