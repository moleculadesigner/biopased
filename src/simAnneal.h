c simulated annealing protocol
 	integer nSAstepMAX
        parameter ( nSAstepMAX = 50 )
        integer nSAstep
        integer delTemp
        integer nSAparMAX,nSAparMX
        parameter (nSAparMAX = 5)   ! number of VariableFF SA param for md run
c
c SAProtcol:  nStep  T  a(SC)  wf(hb128BB)  wf(hb128SB,SS) 
c
c nSApar: nStep,T, sc(vdw),wf(hB128),wf(coulEn)
c
        real SAProtcol(nSAparMAX*nSAstepMAX)  ! nstepMX,Temp,aSoftCore,wfHb128,wfCoul
        character*(charLenMAX) fileSAProtocol
        common/SAProt/nSAstep,SAProtcol,fileSAProtocol,nSAparMX
c
