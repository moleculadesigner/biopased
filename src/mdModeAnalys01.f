c mdModeAnalys01
c
c Yury Vorobjev, 2005
c
c calculation of positionalFluctuationMatrix
c
	subroutine mdModeAnalys01 
c
cx        implicit none
        include 'xyzPDBsize.h'
        include 'xyzPDBinfo.h' 
        include 'engXYZtra.h'
        include 'modeXYZtra.h'
        include 'filedat.h'
c
        include "output.h"
        include "kanalUse.h"
        include "statusMessg_mDyn.h"
c
	integer i,j,ia,ia3,ja,ja3
        integer i3,j3,is,iSH,k
        integer na3,na3MAX,iv
        integer nfile
        character*(charLenMAX) fileMOLX,fileMOLX0
        character*(charLenMAX) fileEssModeEigv
c
        real xyz(3)
        integer nAtomMPFL3
        integer nAtomMPFMX3
        integer iL,jL,ki,kj
        integer iegV,iL3,jL3
        integer nbb
        integer essModeTypeOPT,nModeWrite
        real massAtAv
        real lambdaMIN,big,essModeRMSDperAtomMIN
        real ant,segV,segV2
        real scaleCA,scalePBB
        real scaleEntropy
        logical CONTROL,CONTROL1
        integer kanalp
c
        data atomTypeMPFList /"CA  ","N   ","C   ","O   "/
c
cx        OPT_EssModeType = 2 ! =1 CA of prot backbone for MatrixPositionaFluct analysis
cx                            ! =2 "CA  ","N   ","C   ","O : protBBone atoms
        lambdaMIN = 1.0e-10
        essModeRMSDperAtomMIN = 0.001
        big = 1.0e+10
        scaleCA = 0.320      !empirical to scale Entropy [CA]
        scalePBB = 1.0
c
        massAtAv = 6.0  !  atomWeight CA,N,C,O 
        nAtomTypeMPF=1             ! default: only CA of protein backbone
        if(OPT_EssModeType .eq. 1)then
        nAtomTypeMPF=1
        scaleEntropy=scaleCA
        end if
        if(OPT_EssModeType .eq. 2)then
         nAtomTypeMPF=4 ! all protein backbone atoms
        scaleEntropy=scalePBB  
        end if
c
        nAtomMPFMX = natomMPFMAX 
        nAtomMPFMX3=3*nAtomMPFMX
c
        kanalp = kanalRunOut
        CONTROL = .true.
        CONTROL1 = .false.
c
        if(CONTROL)then
        write(kanalp,*)'mdModeAnalys01 start: '
        write(kanalp,*)'nSnap:',nStepTra
        end if
c
c define proteinBackBone atomList
        iv = 0    ! extract protein backBond atoms by name
        nbb = 4
        nProtBBatRes=4
c
        call getAtomListByName(iv,natom,atomName,ResName,
     &    nbb,atomTypeMPFList,nAtomProtBB,natomMX,
     &    atomNumbProtBBList)
c
        if(CONTROL)then
        write(kanalp,*)'mdModeAnalys01:atomNumbProtBBList:'
        do is=1,nAtomProtBB 
        write(kanalp,'(a8,i5,i5,1x,a4)')'ibb,ia:',
     &  is,atomNumbProtBBList(is),atomName(atomNumbProtBBList(is))
        end do!is
        end if !C
c
        na3 = natom*3
        na3MAX=3*natomMAX
        do i=1,na3MAX
        atomXYZtraAve(i)=0.0
        end do!i
c
	do is = 1,nStepTra
        iSH = na3*(is-1) 
        do ia3=1,na3
        atomXYZtraAve(ia3)=atomXYZtraAve(ia3) +
     &                       atomXYZtra(ia3+iSH)
        end do!ia
        end do!is
c
        ant=1.0/nStepTra
        do ia3=1,na3    
        atomXYZtraAve(ia3)=atomXYZtraAve(ia3)*ant
        end do!ia3
c
        iv=0
        nfile = 0
        fileMOLX='molXYZtraAve.pdb'
        if(nMolNameLett .ge. 1)then
        fileMOLX = molNameArgLine(1:nMolNameLett)//
     &                 '.molXYZtraAve.pdb'
        end if
c
        call writeMOLXYZfile(iv,fileMOLX,nfile,atomXYZtraAve)
c
c atomicRMSD for mdTra
c
        do is=1,nStepTra
        iSH=na3*(is-1)
        atomXYZRmsdTra(is)=0.0
        do i=1,na3
        atomXYZRmsdTra(is)=atomXYZRmsdTra(is)+
     &     (atomXYZtra(i+iSH)-atomXYZtraAve(i))**2
        end do!i
        atomXYZRmsdTra(is)=sqrt(atomXYZRmsdTra(is)/natom)
        end do !is
c
        if(CONTROL)then
        write(kanalp,*)'mdModeAnalys01: iSnap atomXYZRmsdTra:'
        do is=1,nStepTra 
        write(kanalp,*)'atomXYZRmsdTra:',is,atomXYZRmsdTra(is)
        end do!is
        end if
c
c matrix of positional fluctuations
c define atomList included in List
c
        iv = 0    ! extract protein backBond atoms by name [CA,N,C,O]
c
cx      nAtomTypeMPF=4(1) ! full BBlist(CA atom only)
c
        call getAtomListByName(iv,natom,atomName,ResName,
     &    nAtomTypeMPF,atomTypeMPFList,nAtomMPFList,nAtomMPFMX,
     &    atomNumbMPFList)
c
        if(CONTROL)then
        if(nAtomTypeMPF .eq. 1)then
        write(kanalp,*)'mdModeAnalys01: atoms In EssMode analysis :CA'
        elseif (nAtomTypeMPF .eq. 4) then
        write(kanalp,*)
     &  'mdModeAnalys01: atoms In EssMode analysis :CA N C O'
        end if
c
        write(kanalp,*)
     &  'mdModeAnalys01: atoms In EssMode analysis: atomNumbMPFList '
        do is=1,nAtomMPFList
        write(kanalp,'(a12,i5,i5,1x,a4)')'iEs  iaTom:',
     &  is,atomNumbMPFList(is),atomName(atomNumbMPFList(is))
        end do!is
        end if !C
c initialization:
        do i=1,3*nAtomMPFList
        do j=1,3*nAtomMPFList
        atomXiXjAve(i,j)=0.0
        end do!j
        end do!i
c PFM matrix calculation:
        do is=1,nStepTra
        iSH=na3*(is-1)
        do iL=1,nAtomMPFList
        iL3=3*iL-3
        ia=atomNumbMPFList(iL)
        ia3=3*ia-3
        do jL=1,nAtomMPFList
        jL3=jL*3-3
        ja=atomNumbMPFList(jL)
        ja3=3*ja-3
        do ki=1,3
        do kj=1,3
        atomXiXjAve(iL3+ki,jL3+kj)=atomXiXjAve(iL3+ki,jL3+kj) + 
     &   (atomXYZtra(ia3+iSH+ki)-atomXYZtraAve(ia3+ki))*
     &   (atomXYZtra(ja3+iSH+kj)-atomXYZtraAve(ja3+kj))
        end do!jk
        end do!ik
        end do!jL
        end do!iL
        end do !is
c
         iStatus = iStatus + 1              
        write(kanalPStat,*)mStatus,iStatus,messageStatus
     &  ,' positional fluct matrix is done ...'
        write(kanalPStat,*)'matrix order is ', 3*nAtomMPFList
     &  , '...'
c
        if(CONTROL)then
        write(kanalp,*)'mdModeAnalys01:Matrix PositionFluct: is done'
        end if
c
        nAtomMPFL3=3*nAtomMPFList
        do i=1,nAtomMPFL3    
        do j=1,nAtomMPFL3    
        atomXiXjAve(i,j)=atomXiXjAve(i,j)*ant
        end do !j
        if(CONTROL1)then
        write(kanalp,'(a31,i5,1x,2f8.4)')
     &  'atomXiXjAve(i,i):im,ia,m:',i,atomXiXjAve(i,i)
     &   ,atomXiXjAve(i,i+5)
        end if!C
        end do !i
c
c diagonalization
c
        iStatus = iStatus + 1
        write(kanalPStat,*)mStatus,iStatus,messageStatus
     &  ,' start matrix diagonalization ...'
c
        call jacobiGen(atomXiXjAve,nAtomMPFL3,nAtomMPFMX3,
     &       lambdaMPF,eigVectMPF,nrotJac,bJacobi,zJacobi)
c sortr eigenValues by decreasing
        call eigsrt(lambdaMPF,eigVectMPF,nAtomMPFL3,nAtomMPFMX3)         
c
        write(kanalp,*)
     &  'mdModeAnalys01: Essential Mode eigenValues(A^2):'
c
c write eigenValues to file
        fileEssModeEigv = 'molEssModeFreq'
        if(nMolNameLett .ge. 1)then
        fileEssModeEigv = molNameArgLine(1:nMolNameLett)//
     &                 '.molEssModeFreq'//'.dat'
        end if
c
        open(unit=kanalModFreq,file=fileEssModeEigv,form='formatted',
     &       status='unknown')        

c 
        write(kanalModFreq,*)'Essential mode amplitude/Frequency: '
        write(kanalModFreq,*) fileEssModeEigv
c
        segV=0.0
        essModeEntropy = 0.0
        iegV=1
        do i = 1,nAtomMPFL3
        segV = segV + eigVectMPF(i,iegV)**2 
c
        if(lambdaMPF(i) .lt. 0.0)lambdaMPF(i) = 0.0
        essModeRMSDperAtom(i) = sqrt(lambdaMPF(i)/nAtomMPFList)
c
        if(essModeRMSDperAtom(i) .gt.  essModeRMSDperAtomMIN)then
        essModeFreq(i) = 1.15*       ! in units = hw/kT
     &  sqrt(300.0/(lambdaMPF(i)*tempT0traAv*massAtAv))
c
        write(kanalp,'(a5,i5,a22,1x,f6.3,a26,f8.4)')
     &  'mode ',i,
     &  ' essModeRMSD/atom:', sqrt(lambdaMPF(i)/nAtomMPFList),
     &  ' essModeFreq(hw/kT)= ', essModeFreq(i)
c
        write(kanalModFreq,'(a5,i5,a22,1x,f6.3,a26,f8.4)')
     &  'mode ',i,
     &  ' essModeRMSD/atom:', sqrt(lambdaMPF(i)/nAtomMPFList),
     &  ' essModeFreq(hw/kT)= ', essModeFreq(i)
c
        essModeEntropy = essModeEntropy 
     &  - ALOG(1-exp(-essModeFreq(i)))
     &  + essModeFreq(i)/(exp(essModeFreq(i))-1.0)
c
        else
        essModeFreq(i) = big
        essModeRMSDperAtom(i) = 0.0
        end if
c
        end do!i
c 
c Essential Mode Entropy:
        essModeEntropy = essModeEntropy*0.6*tempT0traAv/300.0  ! =T*S in kcal/mol 
c scale by max total mode
        essModeEntropy = essModeEntropy*(nAtom/nAtomMPFList)
     &                    *scaleEntropy
c
        write(kanalp,*)
        write(kanalp,*)'essModeEntropy*T (kcal/mol): ',essModeEntropy,
     &  ' tempT0Av: ',tempT0traAv
c
        write(kanalModFreq,*)
        write(kanalModFreq,*)
     &  '#essModeEntropy*T (kcal/mol): ',essModeEntropy,
     &  ' tempT0Av: ',tempT0traAv
c 
        close(kanalModFreq)
c
        if(CONTROL1)then
        do i=1,nAtomMPFList
        i3=3*i-3
        ia=atomNumbMPFList(i)
        ia3=3*ia-3
        do k=1,3
        xyz(k) = atomXYZtraAve(ia3+k) + 
     &  sqrt(lambdaMPF(iegV))*eigVectMPF(i3+k,iegV)
        end do!k   
c
         write(kanalp,7071)"ATOM  ",
     &   atomNumb(ia),atomName(ia),resName(ia),chName(ia),resNumb(ia),
     &   (xyz(k),k=1,3), atomQ0(ia)
        end do !i
        end if!C1
c
        essModeTypeOPT = OPT_EssModeType
        nModeWrite = 3                   ! write first 3 modes
c
        write(kanalp,*)
     &  'write filesXYZ for:',nModeWrite, ' EssModes ...'
c
        call getModeXYZ(essModeTypeOPT,nModeWrite)
c
	return
7071    format(a6,1x,i4,2x,a4,a4,a1,i4,4x,3f8.3,7x,f8.5) ! PDB rq  
        end
