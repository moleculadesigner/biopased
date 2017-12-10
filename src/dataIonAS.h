c
c dataIonAS.h data structure for IonAS 
c
        include 'ionASdataSiz.h'
c
        real dprobeIon
        real ionMolRad_OPT
        real dotIonden
        integer nIonAScall
        common/sasdat00/dotIonden,dprobeIOn,nIonAScall,
     &  ionMolRad_OPT
c
        real dotIonXYZ(3*ndotIasMAX)
        real dotIonnrm(3*ndotIasMAX)
        real dotIonArea(ndotIasMAX)
        common/sasdat01/ dotIonXYZ,dotIonnrm,dotIonArea
c
        integer dotIon_IATNUM(ndotIasMAX)
        integer dotIon_JATNUM(ndotIasMAX)
        integer ndotIon
        integer nsurfIonAt                 !number of atoms on IonAS
        integer nsurfIonAtList(natomMAX)   !list of atoms on IonAS
        real atSurfIonAr(natomMAX)
        real atSurfIonNrm(3*natomMAX)
        real atSurfIonXYZ(3*natomMAX)
        real atSurfIonProj(natomMAX)
        common/sasdat02/dotIon_IATNUM,ndotIon,dotIon_JATNUM,
     &        nsurfIonAt,nsurfIOnAtList,atSurfIonAr,
     &        atSurfIonNrm,atSurfIonXYZ,atSurfIonProj
        integer head_dotIonNum(natomMAX)
        integer linkListDotIonAt(ndotIasMAX)
        common/sasdat02b/head_dotIonNum,linkListDotIonAt
c
c -- bindSite---------------------------------------
        real bindSiteIonAt01(3*natomMAX)
        common/ionASdat03/bindSiteIonAt01
c
c -------------------------------------------------
        real atDistToIonAS(natomMAX)
        real atIonASexp(natomMAX)
        common/ionASDat04/atDistToIonAS,atIonASexp
c  -----------------------------------------------
c dotIon_ePotential(*), dotIon_eField(*)
        integer dotIon_ePot_occFlag(ndotIasMAX)
        real dotIon_eField(3*ndotIasMAX)
        real dotIon_ePot(ndotIasMAX)
        common/ionASDat06/dotIon_ePot,dotIon_eField
        common/ionASDat07/dotIon_ePot_occFlag
c
c-------------------------------------------------
c
