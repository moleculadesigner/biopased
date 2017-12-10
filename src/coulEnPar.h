c
c OPTion parameters to calculate  dielectric const and electrostatic energy
c 
c Y.N.Vorobjev, 2002
c
	integer coulVar_OPT  !0-clasCol,1-col+nBackgr,2-col+nBackgr+eps=c*R
c                            ! 3-sasDielModelLaZaridis, 4-GenBornModel
        real cReps_OPT
        real epsMol_OPT,epsSol_OPT
cx        real solMolRad_OPT
c
        common /coulPar01/cReps_OPT,coulVar_OPT,
     &                    epsMol_OPT,epsSol_OPT
cx   &                    ,solMolRad_OPT
c 
c eta dModel ne Rabotaet ! error
cx        real AdcPar_OPT           ! param of DistDepDiel function
cx        integer ndcPar_OPT        !          DistDepDiel function
cx        real atomDconst(natomMAX)
cx        common/coulPar02/AdcPar_OPT,ndcPar_OPT,atomDconst   ! param of DistDepDiel function
c
        real qPhosphSc_OPT        !scaling coeff for PhsphGroupQ
        integer iqPhosphSc_OPT    !flag to scale charges of phosphate Gr
        common/coulPar03/qPhosphSc_OPT,iqPhosphSc_OPT
c
c       SASdielModel LaZaridis : dataSASdr.h
        
