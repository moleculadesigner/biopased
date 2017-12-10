c# Hbond 12-8 HX-Y
c Use with xyzPDBsize.h
	integer nHXhb128MAX,nYYhb128MAX
        parameter (nHXhb128MAX = 4)
        parameter (nYYhb128MAX = 7)
        integer nHXhb128MX,nYYhb128MX
        integer nHXhb128,nYYhb128
        integer nHbond128Types
        common /hB128_01/nHXhb128,nYYhb128,nHbond128Types,
     &  hB128pariHjY
c Hx,Y atoms FFieldname List
        character*2 hb128HxList(nHXhb128MAX)
        character*2 hb128YYList(nYYhb128MAX)
        common /hB128_01b/ hb128HxList,hb128YYList
c
c Y-HX   Hn  Ym Rmin(A)   Em(kcal/m) : vdw 12-8  parameters
c (Hx,YY): (1,1), (1,2), (1,3),
        integer hB128pariHjY(nYYhb128MAX*nHXhb128MAX*2)
        real hB128parRE(nYYhb128MAX*nHXhb128MAX*2)               ! Rmin, Emin
        real hb128parAB(nYYhb128MAX*nHXhb128MAX*2)               ! A/r^12-B/r^8
        common /hB128_02/hB128parRE,hb128parAB
c hB128 Angle Parameters:
        real hB128AgLpar(2)            ! cosTet0, sigma**2
        common/hm128_02b/hB128AgLpar
c 
c Hbond List
        integer nHBbondHxYMAX
cx        parameter (nHBbondHxYMAX = 4*nresMAX)
	parameter (nHBbondHxYMAX = 8*nresMAX)
c List of triplets hB128List(): (iatY,jatHx,katX)  
c                                H-bond Atoms in the structure
        integer nHBbondHxYMX
        integer nHBbondHxY
        integer hB128List(3*nHBbondHxYMAX)
        integer hB128TypeList(nHBbondHxYMAX)  != 1(backb-backbone),=2(b-S,S-SideCh)
        real hB128WfList(nHBbondHxYMAX)       ! wieghtFactor for the HbondEng/Force
c                                               for the Hbond. Depends on hB128TypeList(*)
        real hB128PotParList(4*nHBbondHxYMAX) ! A12,B12,rm,em
        common/hB128_03/nHBbondHxY,nHBbondHxYMX,hB128List 
        common/hB128_03b/ hB128PotParList
c hB128 engList:
        real hB128WfScaleALL
        real hBHxYeng128List(nHBbondHxYMAX)  
        common/hB128_04/hBHxYeng128List,hB128TypeList,
     &                   hB128WfList,hB128WfScaleALL
c
c hB128force --> enForce.h
cx        real hbHxYeng128 
cx        real hBHxY128force(3*natomMAX)
cx        common/hB128_05/hBHxYeng128,hBHxY128force
c
        integer hB128atomType(natomMAX)           ! =1,2,3 if Hx(donor);=4,5,6,7,8 if acceptorY
        common/hB128_06/hB128atomType
c
        real rcutHb128                            ! rcutof hB128 =  rMax(Hx-Y)
        common/hB128_07/ rcutHb128
c
