c engOptimization
	real engEps     !Eng tolerance in engMinimizationStop
        real engPoTMin  !result: minimal potEng
        integer nOptIter  !max number of iterations
c
	common/engOpt1/ engEps,engPoTMin,nOptIter
c   
        integer nOptXYZ
        real atomXYZmin(3*natomMAX)
        common/engOpt2/ atomXYZmin,nOptXYZ      !optimalXYZ
c
