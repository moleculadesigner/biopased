c pdbAtXYZ.h       
        real    atomXYZ(3*natomMAX)        ! initialMD PDB xyz  
        real    atomXYZin(3*natomMAX)      ! initial PDBin xyz
        real    atomXYZfin(3*natomMAX)     ! final PDB xyz
	real    atomXYZtmp(3*natomMAX)     ! kabshRot
        common /pdbInfo4/ atomXYZ
        common/pdbInfo5/ atomXYZin,atomXYZfin
c
        integer iDefineMOLec               ! =0 initial
c                                          ! =1 defined XYZ of main Molecule(Prot/NA)
c                                          ! >1 SolvExWatShell
        common/pdbInfo6/iDefineMOLec
c
c end
