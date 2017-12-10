c option Parameters:
        logical OPT_resSeqInp
        logical OPT_PDBresSeq
        logical OPT_AddLoop
        logical OPT_fullProtMD
        logical OPT_LoopMD 
        logical OPT_residMD   
        logical OPT_doMD
        logical OPT_engOptim
        logical OPT_engCalc
        logical OPT_SolvGS
        logical OPT_MDSA
        logical OPT_Hread
        logical OPT_harmAt1PosRst
        logical OPT_HW1MolPosRst
        logical OPT_HW1ResPosRst 
        logical OPT_restrDistA2
        logical OPT_LigRes
        logical OPT_CompactForce
        logical OPT_WRmolAtIntLig
        logical OPT_mdRestart
        logical OPT_WRxyzvtra
        logical OPT_HBond128
        logical OPT_MDREx           !replicaExchange MD
        logical OPT_RigidBodyMd
        logical OPT_RigidBodyCADist
        logical OPT_SolvateExWat
        logical OPT_SolvEFS         !ElectrostaticFieldSolvationModel
        logical OPT_zeroROT
        logical OPT_SolvGBorn
        logical OPT_SolvSASHP
        logical OPT_SolvWbr
        logical OPT_TargAtPos
        logical OPT_EssModeAnalys
        logical OPT_coIonShell
        logical OPT_freeEnergy
        logical OPT_resFileDir
c
        common/option_par/OPT_LoopMD,OPT_residMD,OPT_SolvGS,
     &                    OPT_MDSA,OPT_AddLoop,
     &                    OPT_fullProtMD,OPT_Hread,OPT_engOptim,
     &                    OPT_resSeqInp,OPT_PDBresSeq,
     &                    OPT_harmAt1PosRst,OPT_HW1MolPosRst,
     &                    OPT_HW1ResPosRst,
     &                    OPT_engCalc,OPT_doMD,OPT_LigRes,
     &                    OPT_restrDistA2,OPT_CompactForce,
     &                    OPT_WRmolAtIntLig,OPT_mdRestart,
     &                    OPT_WRxyzvtra,OPT_HBond128,
     &                    OPT_MDREx,OPT_RigidBodyMd,
     &                    OPT_RigidBodyCADist,OPT_SolvateExWat,
     &                    OPT_coIonShell,
     &                    OPT_SolvEFS,OPT_zeroROT,
     &                    OPT_SolvGBorn,OPT_SolvSASHP,OPT_SolvWbr,
     &                    OPT_TargAtPos,OPT_EssModeAnalys,
     &                    OPT_freeEnergy,OPT_resFileDir
c
         integer iSolvateExWat
         integer OPT_iHread
         integer OPT_Shake
         integer OPT_doLigDock 
         common/option_par2/OPT_Shake,OPT_doLigDock
         common/option_par3/OPT_iHread,iSolvateExWat
c
