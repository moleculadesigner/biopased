c arrays solvationGSmodel
        character*(charLenMAX) solvAtmTypeFile,solvAtmDatFile
c
        real atmChargeNeutr(natomMAX)
        real atomSolPar(6*natomMAX)    ! dV,dGref,lamb,alfa,Ra,qat
c
        common/solvatGS01/atomSolPar,atmChargeNeutr
c
        common/solvatGS02/solvAtmTypeFile,solvAtmDatFile
c
        real atomHydPhSolv
        integer ihydrPhSolv
        common /solvatGS03/ atomHydPhSolv,ihydrPhSolv
c
