c replica Exchange MD protocol
c 2006 in development
c replica Exchange MD optimization take LARGE RAM
c
        integer nMCexCycleMAX
        integer nMDsimCycleMAX
        parameter (nMCexCycleMAX = 100 )         ! =1000 regular
        parameter (nMDsimCycleMAX = nMCexCycleMAX)
 	integer nRepLicMAX
        parameter (nRepLicMAX = 8 )   !!y
c
        integer nSnapMixMdRunMAX !MAX number of snaps per one freeMdCycle per RepLica
        parameter (nSnapMixMdRunMAX=10)   !=100 regular
c
        integer nMDsnapRepLMAX   !MAX number of all snaps per RepLica
c                                  accumulated in all nMDrunCycle
        parameter (nMDsnapRepLMAX = 5000)    ! =50000 regular
c
        integer nMDsnapRepL      !current number of snaps DONE per RepLica
        integer nSnapMixMdRun    !current number of snaps per one freeMdCycle per RepLica
        integer nwtraRpLsnap     !store each nwtraRpLsnap to find eMin XYZ
c
        integer engTermMAX
        parameter (engTermMAX = 9)  ! engTermMX in mdRunRepLica.f !
c
        character*(charLenMAX) fileMDRExProtocol
c
c eMdRepL00x.tra files
        character*(charLenMAX) mdEngRepLtraf(nRepLicMAX)  ! files replicasEngTra
        integer kanalmdEngRepLtraf(nRepLicMAX)  ! kanal for file = 51,52,...,60,...
c mMdRepL00x.000n.pdb
        character*(charLenMAX) mMdRepLPDBf(nRepLicMAX)  ! files to write replXYZ.pdb
        common/repL000/mdEngRepLtraf,mMdRepLPDBf,kanalmdEngRepLtraf
c
c INPUT pameters: file: ./MdRExProtocol.inp
c
        integer nRepLic            ! the current total number of XYZ replics
        integer ntimeREMDinit      ! ntimeREMD initial termalization 
        integer ntimeREMDrun1      ! ntimeREMDrun1 free MD run of each replica
c                                  ! in one MD cycle between MC exchange
        integer nMCexCycleMX         ! maximalNumb of MC cycles in the RUN
c ccc
        integer nMDsimCycle          ! current number of MD cycles have DONE
        integer nMCexCycle           ! current number of MC cycles have DONE
        integer nMDcyclePerMCex      ! number of MDcycle per one MCexCycle
c
        common/repL002/ntimeREMDinit,ntimeREMDrun1,
     &         nMCexCycleMX,nMCexCycle,nMDsimCycle,
     &         nMDcyclePerMCex,nSnapMixMdRun,
     &         nMDsnapRepL,nwtraRpLsnap 
c
        integer nStepTraRepL(nRepLicMAX)    ! number of traSnaps for Replics
        integer nRecPdbRepL(nRepLicMAX)     ! number of XYZ.pdb for Replics
        integer nTimeStepRepL(nRepLicMAX)   ! nTimeStep done for replica nRep
        common /repL010/nRecPdbRepL,nStepTraRepL,
     &                 nTimeStepRepL
c replicaTemperature
        real    temp0RepL,deltaT0RepL
c INPUT parameters::
c temp0RepL - lowTemp(~100-200K), 
c deltaT0RepL - tempIncrement (~ 39 K)
c
        real tempRepLic(nRepLicMAX)
        common/repL012/nRepLic,temp0RepL,deltaT0RepL,
     &                tempRepLic
c
        real atomXYZRepL(3*natomMAX*nRepLicMAX) ! {xyz,vxyz;...} for current RepLics
c                                               ! arranged from lowT to highT
        real atomVelRepL(3*natomMAX*nRepLicMAX)
        common/repL020/atomXYZRepL,atomVelRepL
c
        real atomXYZRepLnew(3*natomMAX*nRepLicMAX) ! {xyz,vxyz;...} for RepLics
        real atomVelRepLnew(3*natomMAX*nRepLicMAX)
        common/repL021/atomXYZRepLnew,atomVelRepLnew
c
        real repLijExProbabilityTra(2*nRepLicMAX*nMCexCycleMAX)  !
c probabilityOfExchange between replics i->,j=i+1
        integer repLijTranzitionTra(2*nRepLicMAX*nMCexCycleMAX)  !
c tranzition diagonal(i+1,->i) matrix Traj of replica exchange
        common/repL030/repLijExProbabilityTra,
     &                  repLijTranzitionTra
c engTraj for all repLicas
        real mdEngRepLtra(engTermMAX*nRepLicMAX*nMDsnapRepLMAX)
        real mdEngAvRepLtra(engTermMAX*nRepLicMAX*nMCexCycleMAX)
c 
c engTerms: eP,eP',eG,eV,eHb,eCl,eSol,eRest,Temp 
        common/repL031/mdEngRepLtra,mdEngAvRepLtra
c        
        real mdEngAvRepL1tra(engTermMAX*nMCexCycleMAX) ! engSnap  for LowTemp Replica
        real mdEngRepL1tra(engTermMAX*nMDsnapRepLMAX) !
        common/repL040/mdEngAvRepL1tra,mdEngRepL1tra
c
        real atomXYZRepL0MixRun(3*natomMAX*nRepLicMAX*nSnapMixMdRunMAX) !XYZ 
c                                            allReplics in lastFreeMDMixCycle
        real atomVelRepL0MixRun(3*natomMAX*nRepLicMAX*nSnapMixMdRunMAX) !Vxyz
        real mdEngRepL0MixRun(engTermMAX*nRepLicMAX*nSnapMixMdRunMAX)   !eng of
c                                            allReplics in lastFreeMDMixCycle 
        integer  indMdEngRepL0MixRunMin(nRepLicMAX)
        common/repL042/atomXYZRepL0MixRun,atomVelRepL0MixRun,
     &                 mdEngRepL0MixRun,indMdEngRepL0MixRunMin
c 
c minimal EngTermValue() and minEng snapNumber() for the LAST! lastFreeMDMixCycle
        real mdEngRepL0MixMin(engTermMAX,nRepLicMAX)  !min EngTermValue
        integer iSnapMinRepLEngEterm(engTermMAX,nRepLicMAX) !minEng snapNumber
c                                             for given iRepL,ieTerm
        common/repL043/mdEngRepL0MixMin,iSnapMinRepLEngEterm
c
        real atomXYZRepL0eMinRun(3*natomMAX*nRepLicMAX)     !XYZ for eMin 
c                           over allTimeSteps allReplics in lastFreeMDMixCycle
        real atomVelRepL0eMinRun(3*natomMAX*nRepLicMAX)     !Vxyz
c                           over allTimeSteps allReplics in lastFreeMDMixCycle 
        real mdEngRepL0eMinRun(engTermMAX*nRepLicMAX) 
c
        common/repL044/atomXYZRepL0eMinRun,atomVelRepL0eMinRun,
     &                  mdEngRepL0eMinRun
c
c startPosition in sequential files for the given Replica in the nMCexCycle        
        integer startPosmdEngRepLtra(nRepLicMAX,nMCexCycleMAX)   !startPos in 
c                                                         file mdEngRepLtra( )
c                                                         ! for nRep in the nMCexCycle
        integer startPosmdEngAvRepLtra(nRepLicMAX,nMCexCycleMAX)
        integer startPosmdEngAvRepL1tra(nMCexCycleMAX)
        integer startPosmdEngRepL1tra(nMCexCycleMAX)
        integer startPosatomXYZRepL1tra(nMCexCycleMAX)
        integer startPosatomXYZRepL0MixRun(nRepLicMAX)
        integer startPosEngZRepL0MixRun(nRepLicMAX)
c
        common/repL050/startPosmdEngRepLtra,startPosmdEngAvRepLtra,
     &                  startPosmdEngRepL1tra,startPosmdEngAvRepL1tra,
     &                  startPosatomXYZRepL1tra,
     &                  startPosatomXYZRepL0MixRun,
     &                  startPosEngZRepL0MixRun
c
