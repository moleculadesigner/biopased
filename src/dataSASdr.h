c
c dataSASdr.h data structure for SolventAccessibleSurfece calculation
c
        include 'sasDataSiz.h'
c
        integer nProbe
c
        real dprobe(nprobeMAX)
        real solMolRad_OPT
        real dotden
        integer nSAScall
        common/sasdat00/nProbe,dotden,dprobe,nSAScall,
     &  solMolRad_OPT
c
        real dotXYZ(3*ndotMAX,nprobeMAX)
        real dotnrm(3*ndotMAX,nprobeMAX)
        real dotarea(ndotMAX,nprobeMAX)
        common/sasdat01/ dotXYZ,dotnrm,dotarea
c
        integer dot_IATNUM(ndotMAX,nprobeMAX)
        integer dot_JATNUM(ndotMAX,nprobeMAX)
        integer ndot(nprobeMAX)
        integer nsurfAt(nprobeMAX)
        integer nsurfAtList(natomMAX,nprobeMAX)
        real atSurfAr(natomMAX,nprobeMAX)
        real atSurfNrm(3*natomMAX,nprobeMAX)
        real atSurfXYZ(3*natomMAX,nprobeMAX)
        real atSurfProj(natomMAX,nprobeMAX)
        common/sasdat02/dot_IATNUM,ndot,dot_JATNUM,
     &        nsurfAt,nsurfAtList,atSurfAr,
     &        atSurfNrm,atSurfXYZ,atSurfProj
        integer head_dotNum(natomMAX)
        integer linkListDotAt(ndotMAX)
        common/sasdat03/head_dotNum,linkListDotAt
c
c -- bindSite---------------------------------------
        real bindSiteAt01(3*natomMAX,nprobeMAX)
        common/sasdat04/bindSiteAt01
c
c -------------------------------------------------
c data for D(r) model of Lazaridis
        real atDistToSAS(natomMAX)
        real atSASexp(natomMAX)
        common/sasDat04/atDistToSAS,atSASexp
        real parDrLZmodel(8)
        common/sasDat05/parDrLZmodel
c  -----------------------------------------------
c dot_ePotential(*), dot_eField(*)
        integer dot_eField_occFlag(ndotMAX)
        real dot_eField(3*ndotMAX)
        real dot_ePot(ndotMAX)
        common/sasDat06/dot_ePot,dot_eField,dot_eField_occFlag
c-------------------------------------------------
        character*(charLenMAX) fileSASdotXYZ
        logical OPT_writeSASdotXYZ
        common/sasDat07/fileSASdotXYZ,OPT_writeSASdotXYZ
c
