      PROGRAM PRETLUS
C
C ****************** Version PRETLUS205 ********************
C
C
C     PREPROCESSOR FOR TLUTY - LISTING THE DIMENSIONS OF ARRAYS
C
C *********************************************************
C
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'ITERAT.FOR'
      INCLUDE 'ALIPAR.FOR'
C
      INIT=1
      ITER=0
      CALL START
      STOP
      END

C
C
C     ****************************************************************
C
C

      BLOCK DATA
C     ==========
C
C     Hydrogenic oscillator strentghs
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'ATOMIC.FOR'
C
      DATA ((OSH(I,J),I=1,20),J=1,16)/20*0.,
     * 0.4162,19*0.,7.910D-2,0.6407,18*0.,2.899D-2,0.1193,
     * 0.8421,17*0.,1.394D-2,4.467D-2,0.1506,1.038,16*0.,7.799D-3,
     * 2.209D-2,5.584D-2,0.1793,1.231,15*0.,4.814D-3,1.270D-2,2.768D-2,
     * 6.549D-2,0.2069,1.424,14*0.,3.183D-3,8.036D-3,1.604D-2,3.23D-2,
     * 7.448D-2,0.234,1.616,13*0.,2.216D-3,5.429D-3,1.023D-2,1.87D-2,
     * 3.645D-2,8.315D-2,0.2609,1.807,12*0.,1.605D-3,3.851D-3,6.98D-3,
     * 1.196D-2,2.104D-2,4.038D-2,9.163D-2,0.2876,1.999,11*0.,1.201D-3,
     * 2.835D-3,4.996D-3,8.187D-3,1.344D-2,2.32D-2,4.416D-2,0.1,0.3143,
     * 2.19,10*0.,9.214D-4,2.151D-3,3.711D-3,5.886D-3,9.209D-3,1.479D-2,
     * 2.525D-2,4.787D-2,0.1083,0.3408,2.381,9*0.,7.227D-4,1.672D-3,
     * 2.839D-3,4.393D-3,6.631D-3,1.012D-2,1.605D-2,2.724D-2,5.152D-2,
     * 0.1166,0.3673,2.572,8*0.,5.744D-4,1.326D-3,2.224D-3,3.375D-3,
     * 4.959D-3,7.289D-3,1.097D-2,1.726D-2,2.918D-2,5.513D-2,0.1248,
     * 0.3938,2.763,7*0.,4.686D-4,1.07D-3,1.776D-3,2.656D-3,3.821D-3,
     * 5.455D-3,7.891D-3,1.177D-2,1.843D-2,3.109D-2,5.872D-2,0.133,
     * 0.4202,2.954,6*0.,3.856D-4,8.764D-4,1.443D-3,2.131D-3,3.014D-3,
     * 4.207D-3,5.905D-3,8.456D-3,1.254D-2,1.958D-2,3.298D-2,6.228D-2,
     * 0.1412,0.4467,3.145,5*0./
      DATA ((OSH(I,J),I=1,20),J=17,20)/3.211D-4,
     * 7.270D-4,1.188D-3,1.739D-3,
     * 2.425D-3,3.324D-3,4.556D-3,6.323D-3,8.995D-3,.01328,.0207,.03486,
     * .06584,.1494,0.4731,3.336,4*0.,2.702D-4,6.099D-4,9.916D-4,
     * 1.439D-3,1.984D-3,2.679D-3,3.602D-3,4.877D-3,6.719D-3,9.515D-3,
     * 0.01402,.02182,.03672,.06938,.1575,.4995,3.527,3*0.,2.296D-4,
     * 5.167D-4,8.361D-4,1.204D-3,1.646D-3,2.196D-3,2.905D-3,3.856D-3,
     * 5.180D-3,7.099D-3,.01002,.01474,.02292,.03858,.07292,.1657,.5259,
     * 3.718,2*0.,1.967D-4,4.416D-4,7.118D-4,1.019D-3,1.382D-3,1.825D-3,
     * 2.383D-3,3.112D-3,4.094D-3,5.468D-3,7.468D-3,.01052,.01545,
     * .02402,.04043,.07644,0.1738,.5523,3.909,0./
      END
C
C
C     ****************************************************************
C
C

      SUBROUTINE START
C     ================
C
C     General input and initialization procedure
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      common/hediff/ hcmass,radstr
C
      read(1,*,end=10,err=10) idisk
   10 continue
      call initia
      if(hcmass.gt.0.) call hedif
      nn0=nn
      CALL COMSET
      call prdini
      return
      end
C
C
C     ****************************************************************
C
C
      SUBROUTINE INITIA
C     =================
C
C     driver for input and initializations - "new" routine
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'ATOMIC.FOR'
      INCLUDE 'MODELQ.FOR'
      INCLUDE 'ITERAT.FOR'
      INCLUDE 'ODFPAR.FOR'
      INCLUDE 'ALIPAR.FOR'
      PARAMETER (T15=1.D-15)
      parameter (xcon=8.0935d-21)
      CHARACTER*20 FINSTD
      CHARACTER*40 FILEI
      CHARACTER*4 TYPIOI
      DIMENSION IGLE(18),IGMN(25),IGFE(26),IGNI(28)
      COMMON/STRPAR/IMER,ITR,IC,IL,IP,NLASTE,NHOD,LASV
      COMMON/INUNIT/IUNIT
      common/freqcl/frmin,frmax,nfrecl
      common/nflmax/nflmx,kijmx
      common/ncfe00/ncfe,nlife,nkulev
      common/nfitma/nfitmx,nmer,nvoigt
      DATA IGLE/2,1,2,1,6,9,4,9,6,1,2,1,6,9,4,9,6,1/
      DATA IGMN/2,1,2,1,6,9,4,9,6,1,2,1,6,9,4,9,6,1,
     *          10,21,28,25,6,7,6/
      DATA IGFE/2,1,2,1,6,9,4,9,6,1,2,1,6,9,4,9,6,1,
     *          10,21,28,25,6,25,30,25/
      DATA IGNI/2,1,2,1,6,9,4,9,6,1,2,1,6,9,4,9,6,1,
     *          10,21,28,25,6,25,28,21,10,21/
C
C ----------------------
C Basic input parameters
C ----------------------
C
      CALL READBF(0)
C
C     for a stellar atmosphere
C
      IF(IDISK.EQ.0) THEN
         READ(IBUFF,*) TEFF,GRAV
c        WRITE(6,601) TEFF,GRAV
         GRAV=EXP(2.3025851*GRAV)
C
C        for an accretion disk
C
       ELSE
         READ(IBUFF,*) XMSTAR,XMDOT,RSTAR,RELDST
      END IF
c
c     the rest of input is the same for atmospheres and disks
C
      READ(IBUFF,*) LTE,LTGREY
      READ(IBUFF,*) FINSTD
      CALL NSTPAR(FINSTD)
C
      IF(IDISK.EQ.1) CALL INPDIS
C
C ----------------------------
C Frequency points and weights
C ----------------------------
C
C     NFREAD     - number of "continuum" frequency points which are
C                  read. May be less than NFREQ - the actual number of
C                  frequencies, because the frequency points in lines
C                  can be calculated.
C     FREQ(IJ)   - frequency (in sec**-1) of the IJ-th freq. point
C     W(IJ)      - corresponding frequency quadrature weight
C
      READ(IBUFF,*) NFREAD
      NJREAD=NFREAD
      if(njread.gt.mfreq) CALL QUIT('njread.gt.mfreq',njread,mfreq)
C
      IF(NJREAD.LT.0) THEN
         frmin = FRCMIN
         frmax = FRCMAX
         NFREQ=-NJREAD
         nfreqc=nfreq
         if(nfreq.eq.1) then
            freq(1)=frmin
            w(1)=1.
            ijali(1)=1
         else
         do ij=1,nfreq
            fr=log(frmin)+(log(frmax)-log(frmin))*(ij-1)/(nfreq-1)
            freq(nfreq-ij+1)=exp(fr)
            ijali(ij)=1
         end do 
         end if
C
C        setting simple quadrature weights - trapezoidal integration
C  
         if(nfreq.gt.1) then           
         w(1)=0.5*(freq(1)-freq(2))
         w(nfreq)=0.5*(freq(nfreq-1)-freq(nfreq)) 
         do ij=2,nfreq-1
            w(ij)=0.5*(freq(ij-1)-freq(ij+1))
         end do
         end if
       ELSE
         NFREQC=NJREAD
      END IF
c
c     Initialize opacity table, if there is one
c
      if(ioptab.lt.0) call tabini
c
      if(nfrecl.gt.nfreq) nfrecl=nfreq
      do ij=1,nfrecl
         ijali(ij)=0
      end do
C
c  30 CONTINUE
      IF(NFFIX.EQ.2) THEN
         DO IJ=1,NFREQC
            IJALI(IJ)=1
         END DO
      END IF
C
C ----------------------------------------------------
C     Initialize 1/i*i, 1/i*i*i; and turbulent velocities
C ----------------------------------------------------
C
      DO I=1,NLMX
         X=I
         XI2(I)=UN/(X*X)
         XI3(I)=XI2(I)/X
      END DO
C
      IF(ABS(VTB).LT.1.E3) VTB=VTB*1.E5
      DO ID=1,ND
         IF(VTB.GT.0.) VTURB(ID)=VTB
         IF(IPTURB.EQ.0) VTURB(ID)=0.
         VTURBS(ID)=ABS(VTB)
      END DO
C
C ----------------------------------------------------
C Input parameters for explicit and non-explicit atoms
C ----------------------------------------------------
C
C     Input parameters are read by procedure STATE
C     (see description there)
C
      CALL STATE(0,ID,X1,X2)
      ID=1
      WRITE(10,702) YTOT(ID),WMY(ID),WMM(ID)
      DO I=1,MLEVEL
         ILK(I)=0
         DO J=1,MLEVEL
            ITRA(J,I)=0
         END DO
      END DO
C
C --------------------------------------------------------------
C Input of parameters for explicit ions, levels, and transitions
C --------------------------------------------------------------
C
      ILEV=0
      IATLST=0
      ION=0
      IA=0
      IUNIT=20
      NATOM=0
      NLEVEL=0
      nftmax=0
      WRITE(10,703)
   10 CONTINUE
      READ(IBUFF,*,END=20,ERR=20) IATII,IZII,NLEVSI,ILASTI,ILVLIN,
     *              NONSTD,TYPIOI,FILEI
      IF(ILASTI.EQ.0) THEN
         ION=ION+1
         if(ion.gt.mion) CALL QUIT('ion.gt.mion',ion,mion)
         IATI(ION)=IATII
         IZI(ION)=IZII
         NLEVS(ION)=NLEVSI
         TYPION(ION)=TYPIOI
         FIDATA(ION)=FILEI
         NLLIM(ION)=ILVLIN
         INODF1(ION)=0
         INODF2(ION)=0
         IUPSUM(ION)=0
         ICUP(ION)=16
         FF(ION)=0.
         MODEFF=1
         NFF=0
         IF(IATI(ION).EQ.1.AND.IZI(ION).EQ.0) THEN 
            IUPSUM(ION)=-100
            MODEFF=2
         END IF
         IF(IATI(ION).EQ.1.AND.IZI(ION).EQ.-1) THEN 
            IUPSUM(ION)=0
            ICUP(ION)=0
            MODEFF=3
         END IF
         IF(IATI(ION).EQ.2.AND.IZI(ION).EQ.1) THEN
            MODEFF=2
            ICUP(ION)=32
         END IF
         IF(NONSTD.GT.0) THEN
            READ(IBUFF,*) IUPSUM(ION),ICUP(ION),MODEFF,NFF
          ELSE IF(NONSTD.LT.0) THEN
            READ(IBUFF,*) INODF1(ION),INODF2(ION),FIODF1(ION),
     *                    FIODF2(ION),FIBFCS(ION)
            IKOBS(ION)=IABS(NONSTD)
            if(ispodf.ge.1) then
            IF(INODF1(ION).EQ.0) THEN
               IUNIT=IUNIT+1
               INODF1(ION)=IUNIT
            END IF
            IF(INODF2(ION).EQ.0) THEN
               IUNIT=IUNIT+1
               INODF2(ION)=IUNIT
            END IF
            end if
            IF(FIBFCS(ION).NE.' ') THEN
               IUNIT=IUNIT+1
               INBFCS(ION)=IUNIT
            END IF
            IUPSUM(ION)=1
            ICUP(ION)=0
         END IF
C
         IF(IATI(ION).EQ.IATLST) THEN
            NFIRST(ION)=ILEV
          ELSE
            NFIRST(ION)=ILEV+1
            IATLST=IATI(ION)
            IA=IATEX(IATLST)
            N0A(IA)=NFIRST(ION)
            NATOM=MAX(NATOM,IA)
         END IF
         NLAST(ION)=NFIRST(ION)+NLEVS(ION)-1
         NNEXT(ION)=NLAST(ION)+1
         ILEV=NNEXT(ION)
         if(ilev.gt.mlevel) CALL QUIT('ilev.gt.mlevel',ilev,mlevel)
         IZ(ION)=IZI(ION)+1 
         CHARG2(ION)=IZ(ION)*IZ(ION)
         IF(NFF.GT.0) FF(ION)=EH/H*IZ(ION)*IZ(ION)/NFF/NFF
C
         N0I=NFIRST(ION)
         N1I=NLAST(ION)
         NKI=NNEXT(ION)
         ITRA(NKI,NKI)=MODEFF
         DO II=N0I,N1I
            IEL(II)=ION
            IATM(II)=IA
         END DO
         ILK(NKI)=ION
         IATM(NKI)=IA
C
         IF(NUMAT(IA).EQ.1) THEN
            IATH=IA
            IF(IZ(ION).EQ.1) IELH=ION
            IF(IZ(ION).EQ.0) IELHM=ION
         END IF
         IF(NUMAT(IA).EQ.2) THEN
            IATHE=IA
            IF(IZ(ION).EQ.1) IELHE1=ION
            IF(IZ(ION).EQ.2) IELHE2=ION
         END IF
C         
         WRITE(10,704) ION,TYPION(ION),N0I,N1I,NKI,IZ(ION),
     *                IUPSUM(ION),ICUP(ION),FF(ION)
C   
       ELSE IF(ILASTI.GT.0) THEN
         ENION(ILEV)=0.
         G(ILEV)=ILASTI
         NQUANT(ILEV)=1
         TYPLEV(ILEV)=TYPIOI
         IMODL(ILEV)=0
         IFWOP(ILEV)=0
         IEL(ILEV)=ION
         NKA(IA)=NNEXT(ION)
         if(modref.ge.0) nref(ia)=nka(ia)
         IF(ILASTI.EQ.1.AND.IATII.GT.IZII) THEN
            IF(IATII.LT.25) THEN
               G(ILEV)=IGLE(IATII-IZII)
             ELSE IF(IATII.EQ.25) THEN
               G(ILEV)=IGMN(IATII-IZII)
             ELSE IF(IATII.EQ.26) THEN
               G(ILEV)=IGFE(IATII-IZII)
             ELSE IF(IATII.EQ.28) THEN
               G(ILEV)=IGNI(IATII-IZII)
            ENDIF
         ENDIF
       ELSE
         GO TO 20
      END IF 
      GO TO 10    
   20 CONTINUE
      NION=ION
      NLEVEL=NKI
C     
      IMER=0
      ITR=0
      ITRX=0
      IC=0
      IL=0
      IP=0
      IF(NFREAD.LT.0) THEN
         NLASTE=NFREQC
         NFREQ=NFREQC
       ELSE
         NLASTE=0
         NFREQ=0
      END IF
      LBPFX=.TRUE.
C
      lasv=.false.
      DO ION=1,NION
         CALL RDATA(ION) 
         IF(IFMOFF.gt.0) THEN
            NFF=NQUANT(NLAST(ION))+1
            IF(NFF.GT.0) FF(ION)=EH/H*IZ(ION)*IZ(ION)/NFF/NFF
         END IF
         if(nfitmx.gt.nftmax) nftmax=nfitmx
      END DO
C
      NTRANS=ITR
      NFREQ=NLASTE
      NCON=IC
      if(ntrans.gt.mtrans) CALL QUIT('ntrans.gt.mtrans',ntrans,mtrans)
C
      CALL LEVSET
C
      WRITE(10,705)
      DO I=1,NLEVEL
         WRITE(10,706) I,TYPLEV(I),TYPION(IEL(I)),ENION(I),G(I),
     *                NQUANT(I),IEL(I),ILK(I),IATM(I),
     *                IMODL(I),ILTLEV(I),IIEXP(I),IIFOR(I)
      END DO
C
C  -----------------------------------------------------------
C     useful quantities for later use
C  -----------------------------------------------------------
C
      DO 125 IT=1,NTRANS
         FRQMX(IT)=0.
         LINEXP(IT)=(.NOT.LINE(IT).OR.INDEXP(IT).EQ.0)
         IF(LINEXP(IT)) GO TO 125
         FFMX=0.
         DO IJ=IFR0(IT),IFR1(IT)
            FFMX=MAX(FFMX,FREQ(IJ))
         END DO
         FRQMX(IT)=FFMX
  125 CONTINUE
C
C ------------------------------------------------------
C     setup continuum frequencies (if they are not read)
C ------------------------------------------------------
C
      if(nfreq.gt.mfreq) CALL QUIT('nfreq.gt.mfreq',nfreq,mfreq)
      if(nfreqc.gt.mfreqc) CALL QUIT('nfreqc.gt.mfreqc',nfreqc,mfreqc)
      IF(NFREAD.GT.0 .AND. ISPODF.EQ.0) CALL INIFRC(0)
      IF(ISPODF.GE.1) CALL INIFRS
      IF(FRLMAX.LE.0.) FRLMAX=FREQ(1)
C
C ---------------------------------------------------------
C     setup depth-independent line profiles (sampling mode)
C ---------------------------------------------------------
C
      IF(ISPODF.GE.1) THEN
         TSTD=0.75*TEFF
         DO 80 IT=1,NTRANS
            IF(LINEXP(IT)) GOTO 80
            INDXPA=IABS(INDEXP(IT))
            IF(INDXPA.GE.2 .AND. INDXPA.LE.4) GO TO 80
            CALL DOPGAM(IT,1,TSTD,DOP,AGAM)
            CALL LINSPL(IT,DOP,AGAM)
   80    CONTINUE
      END IF
C
      call rdatax(0,ic,0)
C
      if(icompt.ne.0) ndeptc=nd
c
C  -----------------------------------------------------------
C     read the input model
C  -----------------------------------------------------------
C
c     IF(.NOT.LTGREY) THEN
         if(ispodf.ne.0) CALL INPMOD(ifinpm)
c        IF(ICHANG.NE.0) CALL CHANGE
         if(ifinpm.eq.0) then
         do id=1,nd
            dm(id)=10.**(-7.+0.1*float(id))
            temp(id)=teff*0.75
            elec(id)=1.e12
            dens(id)=1.e-12
            do ii=1,nlevel
               popul(ii,id)=1.e10
            end do
         end do
         temp(nd)=4.*teff
         end if
c     END IF
C
      CALL LINSET(0,0,0,0,0.d0,0.d0,0.d0)
C
C     -----------------------------------------------------------
C     There is no additional input for ODF transitions (i.e. ABS(MODE)=3)
C     from standard input; all information is taken from corresponding 
C     input files (subroutine ODFSET)
C
C     Input parameters for ODF lines or for iron line sampling
C  -----------------------------------------------------------
C
      DOPSTD=SQRT(TWO*BOLK*TEFF/HMASS+VTB*VTB)
      IF(NHOD.GT.0) CALL ODFHYS(DOPSTD)
      IF(ISPODF.EQ.0) THEN
         CALL ODFSET
       ELSE IF(ISPODF.GE.1) THEN
         CALL IROSET
      END IF
C
C  -----------------------------------------------------------
C     select explicit continuum frequencies
C  -----------------------------------------------------------
C
      IF(NFREAD.GT.0 .AND. ISPODF.EQ.0) CALL INIFRC(1)
C
      CALL TRAINI
c
      nncdw=0
      do ibft=1,ntranc
          itr=itrbf(ibft)
          icdw=mcdw(itr)
          if(icdw.ge.1) nncdw=nncdw+1
      end do
      if(nncdw.gt.mmcdw) 
     &  call quit(' Too many pseudo-continua',nncdw,mmcdw)
c
C
C  -----------------------------------------------------------
C     sorting frequencies & new weights
C  -----------------------------------------------------------
C
      DO IJ=NFREQC+1,NFREQ
         IJX(IJ)=1
      END DO
      if(ioptab.ge.0) then
         CALL SRTFRQ
       else
         do ij=1,nfreq
            ijfr(ij)=ij
            jik(ij)=ij
         end do
      end if
C
C  -----------------------------------------------------------
C     other frequency-dependent quantities
C  -----------------------------------------------------------
C
      if(ifryb.gt.0) then
         do ij=1,nfreq
            ijali(ij)=0
         end do
      end if
c
      CALL CORRWM
C
      do ij=1,nfreq
         if(icompt.eq.0) then
            sigec(ij)=sige
         endif
      end do
c
      CALL RTEANG
C
C     -----------------------------------------------------------
C     print out some important transition parameters
C     -----------------------------------------------------------
C
      if(iprini.gt.1) then
c     WRITE(6,607) NTRANS
      DO IT=1,NTRANS
         IF(FR0(IT).GT.0.) THEN
           ALAM=2.997925D18/FR0(IT)
         ELSE
           ALAM=0.
         END IF
c        WRITE(6,608) IT,ILOW(IT),IUP(IT),INDEXP(IT),ICOL(IT),
c    *                IFR0(IT),IFR1(IT),OSC0(IT),FR0(IT),ALAM
      END DO
      end if
C
C     --------------------------------------------------------
C     Evaluation of the photoionization cross-sections in both
C     explicit and fixed-option continuum transitions
C
C     array   CROSS(transition,frequency)
C     --------------------------------------------------------
C
      NFREQB=NFREQ
      IF(IBFINT.GT.0) NFREQB=NFREQC
      if(lasv) call sigave
      DO 220 ITR=1,NTRANS
         IF(LINE(ITR).OR.INDEXP(ITR).EQ.0) GO TO 220
         IC=ITRCON(ITR)
         if(ibf(ic).gt.49.and.ibf(ic).lt.100) go to 220
         ISIK=0
         MODW=IABS(INDEXP(ITR))
         IF(MODW.EQ.5 .OR. MODW.EQ.15) ISIK=1
         DO 210 IJ=1,NFREQB
            FR=FREQ(IJ)
            IF(ISPODF.GE.1) FR=FREQ(IFREQB(IJ))
            CS=SIGK(FR,ITR,ISIK)
            BFCS(IC,IJ)=real(CS)
            IF(FR.LT.FR0PC(ITR)) BFCS(IC,IJ)=0.
            IF(IFWOP(ILOW(ITR)).LT.0) BFCS(IC,IJ)=1.E-30
  210    CONTINUE
  220 CONTINUE
C
      call rdatax(-1,ic,0)
c
      IF(IOPHL1.NE.0 .OR. IOPHL2.NE.0) then
         if(ielh.le.0)
     *   CALL QUIT('ielh.le.0 for iophl1.gt.0',ielh,iophl1)
         CALL OPAHST
      end if
      if(iophe1.gt.0.and.ielhe1.le.0)
     *   CALL QUIT('iophe1.gt.0.and.ielhe1.le.0',iophe1,ielhe1)
      if(iophe2.gt.0.and.ielhe2.le.0)
     *   CALL QUIT('iophe2.gt.0.and.ielhe2.le.0',iophe2,ielhe2)
      if(iphe2c.gt.0.and.ielhe2.le.0) 
     *   CALL QUIT('iophe2c.gt.0.and.ielhe2.le.0',iophe2c,ielhe2)
C
C --------------------------------------------------------------------
C Parameters specifying the overall organization and structure of the
C global matrices of complete linearization
C --------------------------------------------------------------------
C
      if(hmix0.gt.0..and.iconv.eq.0) iconv=1
      IF(HMIX0.GT.0.and.inpc.gt.0.and.inse.gt.0) THEN
         INDL  = 3
         INPC  = 4
         INSE  = 5
         INZD  = 0
      END IF
c
      if(ifryb.gt.0) then
         inhe=0
         inre=1
         indl=0
         inpc=0
         inmp=0
         inse=0
         inzd=0
         do ij=1,nfreq
            ijali(ij)=0
         end do
      end if
c
      nn=nfreqe
      if(inhe.gt.0) nn=nn+1
      if(inre.gt.0) nn=nn+1
      if(inpc.gt.0) nn=nn+1
      if(inmp.gt.0) nn=nn+1
      if(indl.gt.0) nn=nn+1
      if(inzd.gt.0) nn=nn+1
      nngg=nn
      if(inse.gt.0) nn=nn+NLVEXP
      LCHMAT=.FALSE.
c

C
      if(nd.gt.mdepth) CALL QUIT('nd.gt.mdepth',nd,mdepth)
      if(nfreq.gt.mfreq) CALL QUIT('nfreq.gt.mfreq',nfreq,mfreq)
      if(nfreqe.gt.mfrex) CALL QUIT('nfreqe.gt.mfrex',nfreqe,mfrex)
      if(nlevel.gt.mlevel) CALL QUIT('nlevel.gt.mlevel',nlevel,mlevel)
      if(ntrans.gt.mtrans) CALL QUIT('ntrans.gt.mtrans',ntrans,mtrans)
      if(natom.gt.matom) CALL QUIT('natom.gt.matom',natom,matom)
      if(nion.gt.mion) CALL QUIT('nion.gt.mion',nion,mion)
      if(nn.gt.mtot) CALL QUIT('nn.gt.mtot',nn,mtot)
      if(nlambd.gt.mlambd) CALL QUIT('nlambd.gt.mlambd',nlambd,mlambd)
      if(nhod.gt.mhod) CALL QUIT('nhod.gt.mhod',nhod,mhod)
      if(ntrans.gt.32767)
     * call quit(' Too many transitions to define ITRLIN as INTEGER*2',
     *            ntrans,ntrans)
c
      if (icompt.ne.0) ndeptc=nd
      NLIMAX=0
      DO IJ=1,NFREQ
         NLINES(IJ)=0
         DO 50 IT=1,NTRANS
            IF(LINEXP(IT)) GO TO 50
            IF(KIJ(IJ).LT.KFR0(IT)) GO TO 50
            IF(KIJ(IJ).GT.KFR1(IT)) GO TO 50
            IF(IJLIN(IJ).EQ.IT) GO TO 50
            NLINES(IJ)=NLINES(IJ)+1
   50    CONTINUE
         IF(NLINES(IJ).GT.NLIMAX) NLIMAX=NLINES(IJ)
      END DO
c
      nato=99
      WRITE(6,601) NATO, NION, NLEVEL,NLVEXP, NTRANS, ND
  601 FORMAT(' MATOM  = ',I7/
     *       ' MION   = ',I7/
     *       ' MLEVEL = ',I7/
     *       ' MLVEXP = ',I7/
     *       ' MTRANS = ',I7/
     *       ' MDEPTH = ',I7)
c
      WRITE(6,602) NFREQ,KIJMX,NFREQC,NFREQE,NFLMX,NN,NMU
  602 FORMAT(' MFREQ  = ',I7/
     *       ' MFREQP = ',I7/
     *       ' MFREQC = ',I7/
     *       ' MFREX  = ',I7/
     *       ' MFREQL = ',I7/
     *       ' MTOT   = ',I7/
     *       ' MMU    = ',I7)
c
      WRITE(6,603) NFTMAX, NLIMAX 
  603 FORMAT(' MFIT   = ',I7/
     *       ' MITJ   = ',I7)
c
      nzz=6
      nlmxx=80
      nsmxn=1
      nfrq1=1
      if(ispli.eq.5) nfrq1=nfreq
      write(6,605) nncdw,nmer,nvoigt,nzz,nlmxx,nsmxn,nfrq1
  605 format(' MMCDW  = ',I7/
     *       ' MMER   = ',I7/
     *       ' MVOIGT = ',I7/
     *       ' MZZ    = ',I7/
     *       ' NLMX   = ',I7/
     *       ' MSMX   = ',I7/
     *       ' MFREQ1 = ',I7)
c
      write(6,606) numfreq,numtemp,numrho
  606 format(' MFRTAB = ',I7/
     *       ' MTABT  = ',I7/
     *       ' MTABR  = ',I7) 
c
      ncrs=nlevel+5
      nbf=nlevel
      write(6,607) ncrs,nbf,ndeptc,nmuc
  607 format(' MCROSS = ',I7/
     *       ' MBF    = ',I7/
     *       ' MDEPTC = ',I7/
     *       ' MMUC   = ',I7)

      write(6,611) jidn,nkulev,nlife,ncfe
  611 format(/'parameters in ODFPAR.FOR:'//
     *       ' MDODF  = ',I10/
     *       ' MKULEV = ',I10/
     *       ' MLINE  = ',I10/
     *       ' MCFE   = ',I10)
C
c 701 FORMAT(1H1,'*******************************'//
c    *           ' M O D E L   A T M O S P H E R E'//
c    *           ' *******************************'//
c    * ' TEFF   =',F10.0/' LOG G  =',F10.2/)
c  602 FORMAT(1H0//' YTOT   =',F11.5/
c     * ' WMY    =',1PD15.5/' WMM    =',D15.5/)
  702 FORMAT(' YTOT WMY WMM ',F11.5,1P2D15.5)
  703 FORMAT(1H0//' EXPLICIT IONS INCLUDED'/
     *            ' ----------------------'//
     *  ' NO.  ION     N0    N1    NK    IZ  IUPSUM ICUP       FF'/)
  704 FORMAT(1H ,I3,2X,A4,6I6,1PD15.3)
  705 FORMAT(1H0//' EXPLICIT ENERGY LEVELS INCLUDED'/
     *            ' -------------------------------'//
     * ' NO.   LEVEL    ION   ION.EN.(ERG)        G   NQ',
     * '  IEL ILK IAT IMOD ILT IIE IIF'/)
  706 FORMAT(1H ,I4,1X,A10,A4,1PD15.7,0PF10.2,8I5)
c 707 FORMAT(1H0//' TRANSITION PARAMETERS (TOTAL OF',I6,' )'/
c    *            ' ---------------------'//
c    * '    ITR ILOW IUP INDEXP ICOL   IFR0   IFR1     OSC',
c    * '        FR0',8X,' LAMBDA'/)
c 708 FORMAT(1H ,I6,4I5,2I7,1P2D12.3,0PF12.3)
c 709 FORMAT(1H0//' FREQUENCY POINTS AND WEIGHTS - EXPLICIT'/
c    *            ' ---------------------------------------'//
c    * '    IJ  ',7X,'FREQ',13X,'WEIGHT',11X,'PROF'/)
c 710 FORMAT(1H ,I7,1P2D17.8,D15.5,I5,D17.8)
c 712 FORMAT(1H0,
c    *'IOPADD,IRSCT,IOPHMI,IOPH2P,IOPHE1,IOPHE2,IOPHL1,IOPHL2,IPHE2C:'
c    */9I5)
C
      RETURN
      END
C     
C     
C ************************************************************************    
C    
C     
C  
      SUBROUTINE RDATA(ION)
C     =====================
C 
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'ATOMIC.FOR'
      INCLUDE 'MODELQ.FOR'
      INCLUDE 'ITERAT.FOR'
      INCLUDE 'ODFPAR.FOR'
      INCLUDE 'ALIPAR.FOR'
      COMMON/STRPAR/IMER,ITR,IC,IL,IP,NLASTE,NHOD,LASV
      COMMON/INUNIT/IUNIT
      common/nfitma/nfitmx,nmer,nvoigt
      PARAMETER (T15=1.D-15)
      CHARACTER*1 A
c
      ihydp0=0
      IUNIT=21
      OPEN(IUNIT,FILE=FIDATA(ION),STATUS='OLD')
C
C     read the first record - a label for the energy level input
C  
      READ(IUNIT,501) A
  501 FORMAT(A1)
C
C   -----------------------------------------------------
C   input parameters for explicit energy levels
C   -----------------------------------------------------
C  
      DO IL=1,NLEVS(ION)
         I=IL+NFIRST(ION)-1
         IE=IEL(I)
         N0I=NFIRST(IE)
         NKI=NNEXT(IE)
         IQ=I-N0I+1
         X=IQ*IQ
         READ(IUNIT,*) 
     *   ENION(I),G(I),NQUANT(I),TYPLEV(I),ifwop(i),FRODF(I),IMODL(I)
         if(ifwop(i).lt.0.and.i.ne.nlast(ie)) 
     *   CALL QUIT('ifwop(i).lt.0.and.i.ne.nlast(ie)',i,nlast(ie))
         E=ABS(ENION(I))
         E0=E
         IF(E.EQ.0.) E0=EH*IZ(IE)*IZ(IE)/X
         IF(E.GT.1.D-7.AND.E.LT.100.) E0=1.6018D-12*E
         IF(E.GT.100..AND.E.LT.1.D7) E0=1.9857D-16*E
         IF(E.GT.1.D7) E0=H*E
         IF(ENION(I).GE.0.) THEN
            ENION(I)=E0
          ELSE
            ENION(I)=-E0
         END IF
         IF(G(I).EQ.0.) G(I)=2.D0*X
         IF(NQUANT(I).EQ.0) NQUANT(I)=IQ
         iltlev(i)=0
         if(nquant(i).lt.0) then
            iltlev(i)=1
            nquant(i)=iabs(nquant(i))
         end if
         if(ispodf.eq.0 .and. ifwop(i).ge.2) ifwop(i)=0
         if(ifwop(i).lt.0) then
            enion(i)=0.
            IMER=IMER+1
            IMRG(I)=IMER
            IIMER(IMER)=I
         endif
         LBPFX = LBPFX .AND. IMODL(I).EQ.0 .AND. IIFIX(IATM(I)).EQ.0
         if(imodl(i).gt.100) then
            iguide(i)=imodl(i)-100+NFIRST(ION)-1
            imodl(i)=6
         end if
      END DO
      nmer=imer
c
C ----------------------------------------------------------------------
C
C   skip lines if more levels than needed, and skip the continuum 
C   transition label
C
    5 READ(IUNIT,501) A
      IF(A.NE.'*') GO TO 5
      II0=NFIRST(ION)-1
      ILLIM=NLLIM(ION)+II0
      JCORR=0
      JJCX=0
C
C   -----------------------------------------------------
C   input parameters for continuum transitions
C   -----------------------------------------------------
C
   10 CONTINUE
      READ(IUNIT,*,END=19,ERR=15) II,JJ,MODE,IFANCY,ICOLIS,
     *                            IFRQ0,IFRQ1,OSC,CPARAM
      IFFR0=0
      IF(IABS(MODE).GT.100) THEN
         IF(MODE.GT.0) MODE=MODE-100
         IF(MODE.LT.0) MODE=MODE+100
         IFFR0=1
         READ(IUNIT,*) FR0INP
         if(ion.eq.ielh) then
            if(ii.eq.1.and.cutlym.ne.0) fr0inp=cutlym
            if(ii.eq.2.and.cutbal.ne.0) fr0inp=cutbal
         end if
         IF(FR0INP.LT.1.E10) FR0INP=2.997925D18/FR0INP
      END IF
      ITR=ITR+1
      IF(IABS(MODE).EQ.5 .OR. IABS(MODE).EQ.15) THEN
         READ(IUNIT,*) FR0PC(ITR)
      END IF
      IF(JJ.LT.1000) THEN
         IF(II.EQ.1 .OR. JJCX.EQ.1) JCORR=NLEVS(ION)+1-JJ
         II=II+II0
         JJ=JJ+II0+JCORR
       ELSE
         II=II+II0
         JJCX=1
      END IF
      ILOW(ITR)=II
      IUP(ITR)=JJ
      INDEXP(ITR)=MODE
      OSC0(ITR)=OSC
      ICOL(ITR)=ICOLIS
      CPAR(ITR)=CPARAM
      LALI(ITR)=.FALSE.
      LEXP(ITR)=.FALSE.
      IC=IC+1
      NTRANC=IC
      if(ifancy.eq.15) then
         ibf(ic)=ifancy
         itrbf(ic)=itr
         call rdatax(itr,ic,IUNIT)
         go to 10
      end if
      IF(ITRA(II,JJ).EQ.0) THEN
         ITRA(II,JJ)=ITR
       ELSE
         ICOL(ITR)=99
      END IF
      IE=IEL(II)
      N0I=NFIRST(IE)
      NKI=NNEXT(IE)
      LINE(ITR)=.FALSE.
      IFC0(ITR)=IFRQ0
      IFC1(ITR)=IFRQ1
      FR0(ITR)=(ENION(II)-ENION(JJ)+ENION(NKI))/H
      IF(IFFR0.EQ.1) FR0(ITR)=FR0INP
C
C   -----------------------------------------------------
C   Additional input parameters for continuum transitions
C   -----------------------------------------------------
C
C     Only for IFANCY = 2, 3, or 4
C     S0, ALF, BET, GAM  - parameters for evaluation the
C         photoionization cross-section
C
      NTRANC=IC
      IF(IFANCY.GE.2.AND.IFANCY.LE.4) THEN
         READ(IUNIT,*) S0CS(IC),ALFCS(IC),BETCS(IC),GAMCS(IC)
       ELSE IF(IFANCY.EQ.6) THEN
         READ(IUNIT,*) (CTOP(IFIT,IC),IFIT=1,6)
      END IF
C
C   -----------------------------------------------------
C   Additional input parameters for continuum transitions -TOPBASE DATA
C   -----------------------------------------------------
C
C     Only for IFANCY > 100 there are IFANCY-100 fit points
C
C     XTOP(MFIT,MCROSS) -  x = alog10(nu/nu0) of a fit point
C     CTOP(MFIT,MCROSS) -  sigma = alog10(sigma/10^-18) of a fit point
C
C     there are IFANCY-100 fit points
C
      IF(IFANCY.GT.100) THEN
         NFIT=IFANCY-100
         if(nfit.gt.nfitmx) nfitmx=nfit
         IF(NFIT.GT.MFIT) CALL QUIT('nfit.gt.mfit',nfit,mfit)
         READ(IUNIT,*) (XTOP(IFIT,IC),IFIT=1,NFIT)
         READ(IUNIT,*) (CTOP(IFIT,IC),IFIT=1,NFIT)
      END IF
      IBF(IC)=IFANCY
      if(ifancy.gt.49.and.ifancy.lt.100) lasv=.true.
      ITRA(JJ,II)=IC
      ITRBF(IC)=ITR
      ITRCON(ITR)=IC
      IF(FR0(ITR).GT.0.) THEN
         ALAM=2.997925D18/FR0(ITR)
       ELSE
         ALAM=0.
      END IF
C
      if(icolis.eq.0.and.ifancy.le.1.and.osc.eq.0.) then
         zz=iz(iel(ii))
         xq=nquant(ii)
         if(fr0(itr).gt.0..and.xq.gt.0.)
     *   sig0=2.815d-20*zz*zz/(fr0(itr)*1.d-16)**3/xq**5
         if(zz.gt.1.9) sig0=sig0*2.
         if(zz.gt.2.9) sig0=sig0*1.5
         osc0(itr)=sig0
      end if
C
      IF(II.LT.NLAST(ION)) GO TO 10
   15 READ(IUNIT,501,end=19,err=19) A
      IF(A.NE.'*') GO TO 15
C
C  -----------------------------------------------------------
C  Input parameters for line transitions
C  -----------------------------------------------------------
C
   19 CONTINUE
      IIP=0
      JJP=0
   20 CONTINUE
      READ(IUNIT,*,END=30,ERR=30) II,JJ,MODE,IFANCY,ICOLIS,
     *                              IFRQ0,IFRQ1,OSC,CPARAM
      IF(ISPODF.GE.1) THEN
         INDXPA=IABS(MODE)
         LIJP=II.EQ.IIP .AND. JJ.EQ.JJP
         LOD34=INDXPA.EQ.3 .OR. INDXPA.EQ.4
         IF(LIJP .AND. LOD34) GO TO 20
         IIP=II
         JJP=JJ
      END IF
      IFFR0=0
      IF(IABS(MODE).GT.100) THEN
         IF(MODE.GT.0) MODE=MODE-100
         IF(MODE.LT.0) MODE=MODE+100
         IFFR0=1
         READ(IUNIT,*) FR0INP
         IF(FR0INP.LT.1.E10) FR0INP=2.997925D18/FR0INP
      END IF
      IF(JJ.GT.NLEVS(ION)) THEN
         IF(IABS(MODE).EQ.2) THEN
            READ(IUNIT,*) K1,K2,K3,X1,X2,X3,K4
            GO TO 20
         END IF
         IF(IABS(MODE).EQ.1) READ(IUNIT,*) LCMP
         IF(IABS(IFANCY).EQ.1) READ(IUNIT,*) GAMRL
         GO TO 20
      END IF
      ITR=ITR+1
      II=II+II0
      JJ=JJ+II0
      ILOW(ITR)=II
      IUP(ITR)=JJ
      INDEXP(ITR)=MODE
      OSC0(ITR)=OSC
      ICOL(ITR)=ICOLIS
      CPAR(ITR)=CPARAM
      LALI(ITR)=.FALSE.
      LEXP(ITR)=.FALSE.
      IFC0(ITR)=IFRQ0
      IFC1(ITR)=IFRQ1
      ITRCON(ITR)=0
      IF(ITRA(II,JJ).EQ.0) THEN
         ITRA(II,JJ)=ITR
       ELSE
         ICOL(ITR)=99
      END IF
      IE=IEL(II)
      N0I=NFIRST(IE)
      NKI=NNEXT(IE)
      LINE(ITR)=.TRUE.
      IFR0(ITR)=IFRQ0
      IFR1(ITR)=IFRQ1
      FR0(ITR)=(ENION(II)-ENION(JJ))/H
      IF(IFFR0.EQ.1) FR0(ITR)=FR0INP
      IF(OSC.EQ.0..and.nquant(ii).le.20.and.nquant(jj).le.20) THEN
         IF(MODE.NE.3.AND.MODE.NE.4) THEN
            GH=2.D0*NQUANT(II)*NQUANT(II)
            OSC0(ITR)=OSH(NQUANT(II),NQUANT(JJ))*G(II)/GH
         ENDIF
         IF(ifwop(jj).lt.0) THEN
            OSC0(ITR)=0.
            JJ0=NQUANT(JJ-1)
            J20=MIN(20,NLMX)
            IF(J20.GE.JJ0) THEN
               DO JTR=JJ0,J20
                  OSC0(ITR)=OSC0(ITR)+OSH(NQUANT(II),JTR)
               END DO
            END IF
            IF(NLMX.GT.20) THEN
               XII=NQUANT(II)*NQUANT(II)
               SUF=0.
               DO JTR=21,NLMX
                  XJ=JTR
                  XJJ=XJ*XJ
                  XJTR=XJ/(XJJ-XII)
                  SUF=SUF+XJTR*XJTR*XJTR
               END DO
               XITR=(400.-XII)/20.
               OSC0(ITR)=OSC0(ITR)+OSH(NQUANT(II),20)*SUF*XITR*XITR*XITR
            END IF
         END IF
      END IF
C
      IF(FR0(ITR).GT.0.) THEN
         ALAM=2.997925D18/FR0(ITR)
       ELSE
         ALAM=0.
      END IF
      IF(MODE.EQ.0) GO TO 20
      IF(IABS(MODE).EQ.2) GO TO 25
C
C     change the status of treatment of superlines, 
C     if specified in the "offset-setting" record, or
C     because of too loo or too high frequency
C
      IF(IABS(MODE).EQ.3.OR.IABS(MODE).EQ.4) THEN
         IF(II.LT.ILLIM.OR.FR0(ITR).LE.FRLMIN.OR.
     *      FR0(ITR).GE.FRLMAX) INDEXP(ITR)=0
         GO TO 20
      END IF
C
C  -----------------------------------------------------------
C  Additional input parameters for "clasical" line transitions
C   (i.e. those not represented by ODF's - ie ABS(MODE)=1)
C  -----------------------------------------------------------
C
C     LCOMP   -  mode of considering absorption profile:
C             =  false - depth-independent profile
C             =  true  - depth-dependent profile
C     INTMOD  -  mode of evaluating frequency points in the line
C             =  0  - means that frequency points and weights have
C                     already been read amongst the NJREAD or NFREAD
C                     frequencies;
C             ne 0  - frequency points and weights will be evaluated:
C             the meaning of the individual values:
C             =  1  - equidistant frequencies, trapezoidal integration
C             =  2  - equidistant frequencies, Simpson integration
C             =  3  - modified Simpson integration
C                     3-point Simpson integrations with each subsequent
C                     integration interval doubled, until the whole
C                     integration area is covered
C             =  4  - frequencies (in units of standard x) and weights
C                     (for integration over x) are read;
C
C      NF     -  number of frequency points in the line
C                (has the meaning only for INTMOD ne 0)
C
C      XMAX   =  0  - program sets up default XMAX=4.55
C             >  0  - means that the line is assumed symmetric around the
C                     center, frequency points are set up between x=0 and
C                     x=XMAX, where x is frequency difference from the
C                     line center in units of the standard Doppler width
C                     (standard Doppler width is the Doppler width
C                     calculated with standard temperature TSTD and
C                     standard turbulent velocity VTB);
C             <  0  - frequency points are set between x=XMAX and x=-XMAX
C
C  IMPORTANT NOTE: all lines are
C  set by default to the full-profile mode. Thefore, if XMAX was
C  positive, it is reset to XMAX --> -XMAX, and NF --> 2*NF - 1
C
C     TSTD    >  0  - standard temperature (see above);
C             =  0  - the program assigns for the standard temperature the
C                     default value TSTD = 0.75*TEFF
C
      READ(IUNIT,*) LCOMP(ITR),INTMOD(ITR),NF,XMAX,TSTD
      IPROF(ITR)=IFANCY
      IF(IABS(IPROF(ITR)).EQ.1) THEN
         IP=IP+1
         if(ip.gt.mvoigt) CALL QUIT('ip.gt.mvoigt',ip,mvoigt)
         ITRA(JJ,II)=IP
      END IF
      if(ion.eq.ielh) then
         if(ihydpr.eq.1) then
            if(nquant(ii).le.4.and.nquant(jj).lt.20) iprof(itr)=3
          else if(ihydpr.eq.2) then
            if(nquant(ii).le.2.and.nquant(jj).le.10) iprof(itr)=4
         end if
         if(iprof(itr).eq.3) ihydp0=1
         if(iprof(itr).eq.4) ihydp0=2
      end if
         
C
C     if Voigt profile is assumed (ie. if IPROF = 1), an additional
C     input record is required which specifies an evaluation of the
C     relevant damping parameter - see procedure DOPGAM)
C
C      GAMR             - for > 0  - has the meaning of natural damping
C                                    parameter (=Einstein coefficient for
C                                    spontaneous emission)
C                             = 0  - classical natural damping assumed
C                             < 0  - damping is given by a non-standard,
C                                    user supplied procedure GAMSP
C      STARK1           -     = 0  - Stark broadening neglected
C                             < 0  - scaled classical expression
C                                    (ie gam = -STARK1 * classical Stark)
C                             > 0  - Stark broadening given by
C                                    n(el)*[STARK1*T**STARK2 + STARK3]
C      STARK2, STARK3   -   see above
C      VDWH             -   .le.0  - Van der Waals broadening neglected
C                             > 0  - scaled classical expression
C
      IF(IABS(IPROF(ITR)).EQ.1) THEN
         READ(IUNIT,*) GAMAR(IP),STARK1(IP),STARK2(IP),STARK3(IP),
     *                 VDWH(IP) 
      END IF
C
C     change the status of treatment of lines, in virtue of:
C     a) too low frequency (FRLMIN); or
C     b) if specified in the "offset-setting" record
C
      IF(II.LT.ILLIM.OR.FR0(ITR).LE.FRLMIN.OR.FR0(ITR).GE.FRLMAX) THEN
         INDEXP(ITR)=0
         GO TO 20
      END IF
C      
      IF(TSTD.EQ.0.) TSTD=0.75*TEFF
      IF(XMAX.EQ.0.) XMAX=4.55
      if(xmax.gt.0.) then
         xmax=-xmax
         nf=2*nf-1
      end if
      CALL DOPGAM(ITR,1,TSTD,DOP,AGAM)
      INTM=INTMOD(ITR)
      IF(INTM.NE.0) IFR0(ITR)=NLASTE+1
      IF(INTM.NE.0) IFR1(ITR)=NLASTE+NF
      IF(INTM.NE.0) NLASTE=IFR1(ITR)
      CALL LINSET(ITR,IUNIT,IFRQ0,IFRQ1,XMAX,DOP,AGAM)
      GO TO 20
C
C  -----------------------------------------------------------
C  Additional input parameters for a "merged superline" transition
C   (i.e. transition to a merged level, treated by means of an ODF
C     - i.e. for ABS(MODE)=2)
C  -----------------------------------------------------------
C
C     KDO(1),KDO(2),KDO(3),XDO(1),XDO(2),XDO(3),KDO(4) - 
C            have the following meaning:
C            The superline is represented by four frequency intervals.
C            Going away from the peak of ODF, the first interval is 
c            represented by a KDO(1)-point Simpson integration, with a 
C            distance XDO(1) fiducial Doppler widths between the points.
C            The same for the second and third interval.
C            The rest (the interval between the last point and the 
C            coresponding edge) is represented by a KDO(4)-point
C            Simpson integration.
C            The fiducial Doppler width is that corresponding to the
C            effective temperature and the standard turbulent velocity.
C
   25 NHOD=NHOD+1
      if(ifwop(jj).ge.0) CALL QUIT('ifwop(jj).ge.0',ifwop(jj),jj)
      READ(IUNIT,*) (KDO(IHI,NHOD),IHI=1,3),(XDO(IHI,NHOD),IHI=1,3),
     *           KDO(4,NHOD) 
C
C     again, change status of treatment of lines if required
C     
      IF(II.LT.ILLIM.OR.FR0(ITR).LE.FRLMIN.OR.FR0(ITR).GE.FRLMAX) THEN
         INDEXP(ITR)=0
         GO TO 20
      END IF
      write(10,*) 'itr,indexp',itr,indexp(itr)
      JNDODF(ITR)=NHOD
      IF(II.LT.NLAST(ION)) GO TO 20
   30 CONTINUE
C
C  -----------------------------------------------------------
C  Additional input parameters for iron-peak superlevels: 
C   Energy bands  (only read in sampling mode!)
C  -----------------------------------------------------------
C
      IF(ISPODF.EQ.1 .AND. INODF1(ION).GT.0) THEN
         READ(IUNIT,*,ERR=50,END=50) NEVKU(ION)
         DO I=1,NEVKU(ION)
            READ(IUNIT,*,ERR=50,END=50) XEV(I,ION)
         END DO
         READ(IUNIT,*,ERR=50,END=50) NODKU(ION)
         DO I=1,NODKU(ION)
            READ(IUNIT,*,ERR=50,END=50) XOD(I,ION)
         END DO
         IF(NEVKU(ION).GT.MLEVEL) 
     *   CALL QUIT('NEVKU(ION).GT.MLEVEL',NEVKU(ION),MLEVEL)
         IF(NODKU(ION).GT.MLEVEL) 
     *   CALL QUIT('NODKU(ION).GT.MLEVEL',NODKU(ION),MLEVEL)
c  50    CALL QUIT(' Missing data for iron superlevels',0,0)
   50    CONTINUE
      END IF
C
      CLOSE(IUNIT)
      nvoigt=ip
      if(itr.gt.mtrans) CALL QUIT('ntrans.gt.mtrans',itr,mtrans)
c
c     initialization of Lemke hydrogen line broadening tables, if needed
c
      if(ion.eq.ielh) then
         if(ihydp0.eq.1) ihydpr=21
         if(ihydpr.eq.2) ihydpr=22
         if(ihydp0.gt.0) then
            call lemini
c           call xenini
         end if
      end if
c
      RETURN
      END
c
c
C    *****************************************************************
c
C
      SUBROUTINE NSTPAR(FINSTD)
C     ==========================
C
C     setting up the default values of various input flags, and
C     input of non-standard values of various input flags and parameters
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'ATOMIC.FOR'
      INCLUDE 'MODELQ.FOR'
      INCLUDE 'ITERAT.FOR'
      INCLUDE 'ODFPAR.FOR'
      INCLUDE 'ALIPAR.FOR'
      COMMON/FLXAUX/T4,PGAS,PRAD,PGM,PRADM,ITGMAX,ITGMX0
      common/freqcl/frmin,frmax,nfrecl
      common/quasun/tqmprf,iquasi,nunalp,nunbet,nungam,nunbal
      common/hediff/ hcmass,radstr
      common/irwint/iirwin
      common/deridt/dert
      common/icnrsp/iconrs
      common/imucnn/imucon
      common/ichndm/ichanm
      common/ipricr/iprcrs,nprcrs
C 
      PARAMETER(MVAR=213)
      PARAMETER(INPFI=4)
      CHARACTER*(*) FINSTD
      CHARACTER*6 PVALUE(MVAR)
      CHARACTER*80 TEXT
      CHARACTER VARNAM(MVAR)*6
      CHARACTER*20 BLNK
      CHARACTER*6 BLNK6
C
      DATA VARNAM /'ISPLIN','IRTE  ','IBC   ','ILMCOR','ILPSCT',
     *             'ILASCT','DJMAX ','NTRALI','IPSLTE','IOPTAB',
     *             'IFMOL ','IFENTR','NFRECL','IFRYB ','IFRAYL',
     *             'HCMASS','RADSTR','BERGFC','IHYDPR','IIRWIN',
     *             'ICOMPT','IZSCAL','IBCHE ','IVISC ','ALPHAV',
     *             'ZETA0 ','ZETA1 ','FRACTV','DMVISC','REYNUM',
     *             'IFZ0  ','IHESO6','ICOLHN',
     *             'IFALI ','IFPOPR','JALI  ','IFRALI','IFALIH',
     *             'IFPREC','IELCOR','ICHC  ','IRSPLT','IATREF',
     *             'MODREF','IACPP ','IACDP ','IFLEV ','IDLTE ',
     *             'POPZER','POPZR2','POPZCH','NITZER','RADZER',
     *             'IFDIEL','IFCHTR','SHFAC ',
     *             'QTLAS ','ITLUCY','IACLT ','IACLDT','IFMOFF',
     *             'IOVER ','ITLAS ','NITER ','NLAMBD','IFRSET',
     *             'ND    ','JIDS  ','IDMFIX','ITNDRE',
     *             'NMU   ','IOSCOR',
     *             'NELSC ','IHECOR','IBFINT','IRDER ',
     *             'CHMAX ','ILDER ','IBPOPE','CHMAXT','NLAMT ',
     *             'INTRPL','ICHANG','IFIXMO',
     *             'INHE  ','INRE  ','INPC  ','INZD  ','INSE  ',
     *             'INMP  ','INDL  ','NDRE  ','TAUDIV','IDLST ',
     *             'NRETC ','ICONV ','IPRESS','ITEMP ','IPRINT',
     *             'ITMCOR','ICONRE','IDEEPC','NDCGAP','CRFLIM',
     *             'IOPADD','IRSCT ','IOPHMI','IOPH2P','IOPHE1',
     *             'IOPHE2','IFCIA ',
     *             'IQUASI','NUNALP','NUNBET','TQMPRF',
     *             'IACC  ','IACD  ','KSNG  ','ITEK  ','ORELAX',
     *             'IWINBL','ICOMGR',
     *             'ICRSW ','SWPFAC','SWPLIM','SWPINC',
     *             'TAUFIR','TAULAS','ABROS0','TSURF ','ALBAVE',
     *             'DION0 ','NDGREY','IDGREY','NCONIT','IPRING',
     *             'DM1   ','ABPLA0','ABPMIN','ITGMAX','NNEWD ',
     *             'IHM   ','IH2   ','IH2P  ','IFTENE',
     *             'TRAD  ','WDIL  ',
     *             'TDISK ',
     *             'HMIX0 ','MLTYPE','VTB   ','IPTURB','ILGDER',
     *             'XGRAD ','STRL1 ','STRL2 ','STRLX ',
     *             'FRCMAX','FRCMIN','FRLMAX','FRLMIN','CFRMAX',
     *             'DFTAIL','NFTAIL','TSNU  ','VTNU  ','DDNU  ',
     *             'IELNU ','CNU1  ','CNU2  ','ISPODF',
     *             'DPSILG','DPSILT','DPSILN','DPSILD',
     *             'ICOMST','ICOMDE','ICOMBC','ICOMVE','ICOMRT',
     *             'ICMDRA','KNISH ','FRLCOM','ICHCOO',
     *             'NCFOR1','NCFOR2','NCCOUP','NCITOT','NCFULL',
     *             'IFPRD ','XPDIV ',
     *             'IPRINI','IBINOP','IDCONZ',
     *             'ICOOLP','IPRIND','IPRINP','ICHCKP','IPOPAC',
     *             'ILBC  ','IUBC  ','DERT  ','ICONRS','IMUCON',
     *             'IFPRAD','ICHANM','CUTLYM','CUTBAL','IHXENB',
     *             'IHGOM ','HGLIM ','IPRCRS','NPRCRS'/
C
      DATA PVALUE /'     0','     0','     3','     3','     1',
     *             '     0',' 1.D-3','     3','     0','     0',
     *             '     0','     0','     0','     0','     0',
     *             '    0.','    0.','    1.','     0','     0',
     *             '     0','     0','     1','     0','   0.1',
     *             '   0.0','   0.0','    -1','  0.01','    0.',
     *             '     9','     0','     0',
     *             '     5','     4','     1','     0','     0',
     *             '     1','   100','     0','     1','     1',
     *             '     1','     7','     4','     0','  1000',
     *             '1.D-20','1.D-20','1.D-15','     1','1.D-20',
     *             '     0','     0','    0.',
     *             ' 1.D30','     0','     7','     4','     0',
     *             '     1','   100','    30','     2','     0',
     *             '    70','     0','     1','     1',
     *             '     3','     0',
     *             '     0','     0','     1','     3',
     *             ' 1.D-3','     0','     1','  0.01','     1',
     *             '     0','     0','     0',
     *             '     1','     2','     3','     0','     4',
     *             '     0','     0','     0','   0.5','     5',
     *             '     0','     0','     0','     0','     0',
     *             '     0','     1','     2','     2','   0.7',
     *             '     0','     0','     0','     0','     0',
     *             '     0','     0',
     *             '     0','     3','     0','    0.',
     *             '     7','     4','     0','     4','  1.D0',
     *             '    -1','     0',
     *             '     0',' 1.D-1',' 1.D-3','  3.D0',
     *             ' 1.D-7',' 316.0','   0.4','    0.','    0.',
     *             '    1.','     0','     0','     0','     0',
     *             ' 1.D-3',' 3.D-1',' 1.D-5','    10','     0',
     *             '     0','     0','     0','     0',
     *             '    0.','    0.',
     *             '    0.',
     *             '   -1.','     1','    0.','     1','     0',
     *             '    0.',' 0.001','  0.02','1.D-10',
     *             '    0.',' 1.D12','    0.',' 1.D13','    0.',
     *             '  0.25','    21','    0.','    0.','  0.75',
     *             '     0','   4.5','    3.','     0',
     *             '   10.','  1.25','   10.','  1.25',
     *             '     1','     1','     1','     0','     0',
     *             '     0','     0','8.2D14','     1',
     *             '     0','     1','     0','     1','     1',
     *             '     0','  3.D0',
     *             '     0','     1','    31',
     *             '     0','     0','     1','     0','     0',
     *             '     0','     0','  0.01','    10','    10',
     *             '     1','     1','    0,','    0.','     0',
     *             '     0',' 1.D18','     0','     0'/
C
      DATA BLNK/'                    '/,BLNK6/'      '/
C
      IF(FINSTD.NE.BLNK) 
     *   OPEN(UNIT=INPFI,FILE=FINSTD,STATUS='UNKNOWN')
C
      DO ID=1,MDEPTH
         CRSW(ID)=UN
      END DO
C
      INDV=-1
C
C    go through the input file line by line
c
   10 CONTINUE
      K0=1
      READ(INPFI,500,END=70,ERR=70) TEXT
  500 FORMAT(A)
   20 CONTINUE
      CALL GETWRD(TEXT,K0,K1,K2)
      IF(K1.EQ.0) GO TO 60
      K0=K2+2
      IF(TEXT(K1:K2).EQ.'=') GO TO 20
      INDV=-INDV
      IF(INDV.EQ.1) THEN
         DO 40 I=1,MVAR
            IF(TEXT(K1:K2).EQ.VARNAM(I)(1:K2-K1+1)) GO TO 50
   40    CONTINUE
         CALL GETWRD(TEXT,K0,K1,K2)
         IF(K1.EQ.0) THEN
            K0=1
   45       READ(INPFI,500,END=70) TEXT
            CALL GETWRD(TEXT,K0,K1,K2)
            IF(K1.EQ.0) GO TO 45
         END IF
         K0=K2+2
         INDV=-INDV
         GO TO 20
   50    CONTINUE
         IVAR=I
       ELSE
         PVALUE(IVAR)=BLNK6
         PVALUE(IVAR)(6-K2+K1:6)=TEXT(K1:K2)
      END IF
      GO TO 20
   60 CONTINUE
      GO TO 10
   70 CONTINUE
C
      DO I=1,MVAR
         WRITE(84,684) PVALUE(I)
  684    FORMAT(1X,A)
      END DO
C
      CLOSE(UNIT=84)
      REWIND(84)
      READ(84,*) 
     *             ISPLIN,IRTE  ,IBC   ,ILMCOR,ILPSCT,
     *             ILASCT,DJMAX ,NTRALI,IPSLTE,IOPTAB,
     *             IFMOL ,IFENTR,NFRECL,IFRYB ,IFRAYL,
     *             HCMASS,RADSTR,BERGFC,IHYDPR,IIRWIN,
     *             ICOMPT,IZSCAL,IBCHE ,IVISC ,ALPHAV,
     *             ZETA0 ,ZETA1 ,FRACTV,DMVISC,REYNUM,
     *             IFZ0  ,IHESO6,ICOLHN,
     *             IFALI ,IFPOPR,JALI  ,IFRALI,IFALIH,
     *             IFPREC,IELCOR,ICHC  ,IRSPLT,IATREF,
     *             MODREF,IACPP ,IACDP ,IFLEV ,IDLTE ,
     *             POPZER,POPZR2,POPZCH,NITZER,RADZER,
     *             IFDIEL,IFCHTR,SHFAC ,
     *             QTLAS ,ITLUCY,IACLT ,IACLDT,IFMOFF,
     *             IOVER ,ITLAS ,NITER ,NLAMBD,IFRSET,
     *             ND    ,JIDS  ,IDMFIX,ITNDRE,
     *             NMU   ,IOSCOR,
     *             NELSC ,IHECOR,IBFINT,IRDER ,
     *             CHMAX ,ILDER ,IBPOPE,CHMAXT,NLAMT ,
     *             INTRPL,ICHANG,IFIXMO,
     *             INHE  ,INRE  ,INPC  ,INZD  ,INSE  ,INMP  ,
     *             INDL  ,NDRE  ,TAUDIV,IDLST ,NRETC ,
     *             ICONV ,IPRESS,ITEMP ,IPRINT,
     *             ITMCOR,ICONRE,IDEEPC,NDCGAP,CRFLIM,
     *             IOPADD,IRSCT ,IOPHMI,IOPH2P,IOPHE1,
     *             IOPHE2,IFCIA ,
     *             IQUASI,NUNALP,NUNBET,TQMPRF,
     *             IACC  ,IACD  ,KSNG  ,ITEK  ,ORELAX,
     *             IWINBL,ICOMGR,
     *             ICRSW ,SWPFAC,SWPLIM,SWPINC,
     *             TAUFIR,TAULAS,ABROS0,TSURF ,ALBAVE,
     *             DION0 ,NDGREY,IDGREY,NCONIT,IPRING,
     *             DM1   ,ABPLA0,ABPMIN,ITGMAX,NNEWD ,
     *             IHM   ,IH2   ,IH2P  ,IFTENE,
     *             TRAD  ,WDIL  ,
     *             TDISK ,
     *             HMIX0 ,MLTYPE,VTB   ,IPTURB,ILGDER,
     *             XGRAD ,STRL1 ,STRL2 ,STRLX ,
     *             FRCMAX,FRCMIN,FRLMAX,FRLMIN,CFRMAX,
     *             DFTAIL,NFTAIL,TSNU  ,VTNU  ,DDNU  ,
     *             IELNU ,CNU1  ,CNU2  ,ISPODF,
     *             DPSILG,DPSILT,DPSILN,DPSILD,
     *             ICOMST,ICOMDE,ICOMBC,ICOMVE,ICOMRT,
     *             ICMDRA,KNISH ,FRLCOM,ICHCOO,
     *             NCFOR1,NCFOR2,NCCOUP,NCITOT,NCFULL,
     *             IFPRD ,XPDIV ,
     *             IPRINI,IBINOP,IDCONZ,
     *             ICOOLP,IPRIND,IPRINP,ICHCKP,IPOPAC,
     *             ILBC  ,IUBC  ,DERT  ,ICONRS,IMUCON,
     *             IFPRAD,ICHANM,CUTLYM,CUTBAL,IHXENB,
     *             IHGOM ,HGLIM ,IPRCRS,NPRCRS
C      
      IF(LTGREY) ISPODF=0
      IF(LTE) IFLEV=1
      LCHC=.FALSE.
      IF(ICHC.EQ.1) LCHC=.TRUE.
      NFFIX=IFRALI
      if(frcmax.lt.1.e6) frcmax=frcmax*1.e15
      IF(FRLMAX.EQ.0.) FRLMAX=1.D11*CNU1*TEFF
      if(idisk.eq.0.and.cfrmax.eq.0.) cfrmax=2.
      if(trad.ne.0.) iwinbl=-1
      if(nitzer.gt.itek) nitzer=itek
      if(nitzer.gt.iacc-iacd) nitzer=iacc-iacd
      if(irsct.ne.0.or.iophmi.ne.0.or.ioph2p.ne.0.or.
     *   iophe1.ne.0.or.iophe2.ne.0) iopadd=1
c
      if(ioptab.lt.0.or.ifmol.gt.0) ielcor=-1
c
      RRDIL=un
      IF(IDISK.EQ.0) IFZ0=-1
      ITGMX0=ITGMAX
      DO ITL=1,NITER+1
      NITLAM(ITL)=0
      END DO
      IF(NLAMBD.LT.0) THEN
         NLAMBD=-NLAMBD
         IF(LTE) NLAMBD=1
         DO ITL=1,12
            NITLAM(ITL)=NLAMBD
         END DO
         DO ITL=13,NITER+1
            NITLAM(ITL)=2
         END DO
       ELSE IF(NLAMBD.GT.0) THEN
         IF(LTE) NLAMBD=1
          DO ITL=1,NITER+1
             NITLAM(ITL)=NLAMBD
         END DO
      END IF
      IF(ILMCOR.GE.3) ILPSCT=1
C
      IF(IDISK.EQ.1.AND.INZD.EQ.0.AND.IZSCAL.EQ.0.AND.IVISC.LE.1) THEN
         if(ifryb.eq.0) then
            INZD=4
            INSE=5
         end if
      END IF
C
      IF(IFIXMO.GT.0) THEN
         INHE=0
         INRE=0
         INPC=0
         INZD=0
         INSE=1
      END IF
c
      if(iprcrs.gt.0) then
         niter=0
         nlambd=1
      end if
C
c     initialize the convection parameters
c     
      aconml=1./8.
      bconml=half
      cconml=16.
      if(mltype.eq.2) then
         aconml=1.
         bconml=2.
         cconml=16.
      end if
      nungam=0
      nunbal=0
      if(iquasi.gt.0) call getlal
c      
      if(nd.gt.mdepth) CALL QUIT('nd.gt.mdepth',nd,mdepth)
      if(ndgrey.gt.mdepth) CALL QUIT('ndgrey.gt.mdepth',ndgrey,mdepth)
      if(nlambd.gt.mlambd) CALL QUIT('nlambd.gt.mlambd',nlambd,mlambd)
      if(iacc.le.2) CALL QUIT('Ng too early',iacc,iacc)
      if(nmu.gt.mmu) CALL QUIT('nmu.gt.mmu',nmu,mmu)
      RETURN
      END
C
C
C    ***************************************************************
C
C
      SUBROUTINE NSTOUT
C     =================
C
C     Diagnostic print of the input flags and parameters
c
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'ATOMIC.FOR'
      INCLUDE 'MODELQ.FOR'
      INCLUDE 'ITERAT.FOR'
      INCLUDE 'ODFPAR.FOR'
      INCLUDE 'ALIPAR.FOR'
C
      if(lchc) ichc=1
      write(6,602)
     *             ISPLIN,IRTE  ,IBC   ,ILMCOR,ILPSCT,
     *             ILASCT,DJMAX ,NTRALI,IPSLTE,
     *             ICOMPT,IZSCAL,IBCHE ,IVISC ,
     *             IFALI ,IFPOPR,JALI  ,IFRALI,
     *             IFPREC,IELCOR,ICHC  ,IRSPLT,IATREF,
     *             MODREF,IACPP ,IACPD ,IFLEV, IDLTE,
     *             POPZER,POPZR2,RADZER,NITZER,IFDIEL,
     *             IOVER ,ITLAS ,NITER ,NLAMBD,ND,
     *             JIDS  ,IDMFIX,NMU   ,
     *             NELSC ,IHECOR,IBFINT,IRDER ,
     *             CHMAX ,ILDER ,IBPOPE,CHMAXT,NLAMT
      write(6,603) INTRPL,ICHANG,
     *             INHE  ,INRE  ,INPC  ,INSE  ,INMP  ,
     *             INDL  ,NDRE  ,TAUDIV,IDLST ,NRETC ,
     *             ICONV ,IPRESS,ITEMP ,IPRINT,
     *             IOPADD,IRSCT ,IOPHMI,IOPH2P,IOPHE1,IOPHE2,
     *             IACC  ,IACD  ,KSNG  ,ITEK  ,ORELAX,
     *             IWINBL,
     *             ICRSW ,SWPFAC,SWPLIM,SWPINC,
     *             IFPRD ,XPDIV
      write(6,604) TRAD  ,WDIL  ,
     *             HMIX0 ,VTB/1.e5   ,IPTURB,
     *             XGRAD ,STRL1 ,STRL2 ,STRLX ,
     *             FRCMAX,FRCMIN,FRLMAX,FRLMIN,CFRMAX,
     *             DFTAIL,NFTAIL,TSNU  ,VTNU/1.e5  ,DDNU  ,
     *             IELNU ,CNU1  ,CNU2  ,ISPODF,
     *             DPSILG,DPSILT,DPSILN,DPSILD
      write(6,605) ICOMST,ICOMDE,ICOMBC,
     *             ICMDRA,KNISH,
     *             NCFOR1,NCFOR2,NCCOUP,NCITOT,NCFULL
      IF(LTGREY)   WRITE(6,606)
     *             TAUFIR,TAULAS,ABROS0,TSURF ,ALBAVE,
     *             ABPLA0,ABPMIN,
     *             DION0 ,NDGREY,IDGREY,NCONIT,IPRING,
     *             IHM   ,IH2   ,IH2P
c
  602 FORMAT('ISPLIN=',I6,2X,'IRTE  =',I6,2X,'IBC   =',I6,2X,
     *       'ILMCOR=',I6,2X,'ILPSCT=',I6,2X/
     *       'ILASCT=',I6,2X,'DJMAX =',F6.3,2X,'NTRALI=',I6,2X,
     *       'IPSLTE=',I6,2X/
     *       'ICOMPT=',I6/
     *       'IZSCAL=',I6,2X,'IBCHE =',I6,2X,'IVISC =',I6,2X/
     *       'IFALI =',I6,2X,'IFPOPR=',I6,2X,'JALI  =',I6,2X,
     *       'IFRALI=',I6,2X/
     *       'IFPREC=',I6,2X,'IELCOR=',I6,2X,'ICHC  =',I6,2X,
     *       'IRSPLT=',I6,2X,'IATREF=',I6,2X/
     *       'MODREF=',I6,2X,'IACPP =',I6,2X,'IACPD =',I6,2X,
     *       'IFLEV =',I6,2X,'IDLTE =',I6,2X/
     *       'POPZER=',1PE6.0,2X,'POPZR2=',1PE6.0,2X,
     *       'RADZER=',1PE6.0,2X,'NITZER=',I6,2X,'IFDIEL=',I6,2X/
     *       'IOVER =',I6,2X,'ITLAS =',I6,2X/
     *       'NITER =',I6,2X,'NLAMBD=',I6,2X,'ND    =',I6/
     *       'JIDS  =',I6,2X,'IDMFIX=',I6/
     *       'NMU   =',I6/
     *       'NELSC =',I6,2X,'IHECOR=',I6,2X,'IBFINT=',I6,2X,
     *       'IRDER =',I6,2X,'CHMAX =',0PF6.3,2X/
     *       'ILDER =',I6,2X,'IBPOPE=',I6,2X,
     *       'CHMAXT=',F6.3,2X,'NLAMT =',I6,2X)
  603 FORMAT('INTRPL=',I6,2X,'ICHANG=',I6,2X/
     *       'INHE  =',I6,2X,'INRE  =',I6,2X,'INPC  =',I6,2X,
     *       'INSE  =',I6,2X,'INMP  =',I6,2X/
     *       'INDL  =',I6,2X,'NDRE  =',I6,2X,'TAUDIV=',F6.3,2X,
     *       'IDLST =',I6,2X,'NRETC =',I6,2X/
     *       'ICONV =',I6,2X,'IPRESS=',I6,2X,'ITEMP =',I6,2X,
     *       'IPRINT=',I6,2X/
     *       'IOPADD=',I6,2X,'IRSCT =',I6,2X,'IOPHMI=',I6,2X,
     *       'IOPH2P=',I6,2X,'IOPHE1=',I6,2X/'IOPHE2=',I6,2X/
     *       'IACC  =',I6,2X,'IACD  =',I6,2X,'KSNG  =',I6,2X,
     *       'ITEK  =',I6,2X,'ORELAX=',F6.3,2X/
     *       'IWINBL=',I6,2X/
     *       'ICRSW =',I6,2X,'SWPFAC=',F6.3,2X,'SWPLIM=',F6.3,2X,
     *       'SWPINC=',F6.3,2X/
     *       'IFPRD =',I6,2X,'XPDIV =',F6.1/)
  604 FORMAT('TRAD  =',F6.0,2X,'WDIL  =',F5.3/
     *       'HMIX0 =',F6.1,2X,'VTB   =',F6.0,2X,I6/
     *       'XGRAD =',F6.2,2X,'STRL1 =',1PE6.0,2X,'STRL2 =',E6.0/
     *       'STRLX =',1PE6.0/
     *       'FRCMAX=',1PE6.0,2X,'FRCMIN=',1PE6.0,2X,
     *       'FRLMAX=',1PE6.0,2X,'FRLMIN=',1PE6.0,/,
     *       'CFRMAX=',0PF6.2,/
     *       'DFTAIL=',0PF6.3,2X,'NFTAIL=',I6,/
     *       'TSNU  =',F6.0,2X,'VTNU  =',F6.2,2X,'DDNU  =',F6.3,/
     *       'IELNU =',I6,2X,'CNU1  =',F6.2,2X,'CNU2  =',F6.2,/
     *       'ISPODF=',I6,/
     *       'DPSILG=',F6.2,2X,'DPSILT=',F6.2,2X,'DPSILN=',F6.2,2X,
     *       'DPSILD=',F6.2,/)
  605 FORMAT('ICOMST=',I6,2X,'ICOMDE=',I6,2X,'ICOMBC=',I6/
     *       'ICMDRA=',I6,2X,'KNISH =',I6/
     *       'NCFOR1=',I6,2X,'NCFOR2=',I6,2X,'NCCOUP=',I6,2X,
     *       'NCITOT=',I6,2X,'NCFULL=',I6)
  606 FORMAT('TAUFIR=',1PE6.0,2X,'TAULAS=',0PF6.1,2X,'ABROS0=',F6.3,2X,
     *       'TSURF =',F6.3,2X,'ALBAVE=',F6.3/
     *       'ABPLAO=',F6.3,2X,'ABPMIN=',1PE6.0/
     *       'DION0 =',0PF6.3,2X,'NDGREY=',I6,2X,'IDGREY=',I6,2X,
     *       'NCONIT=',I6,2X,'IPRING=',I6,2X/
     *       'IHM   =',I6,2X,'IH2   =',I6,2X,'IH2P  =',I6,2X/)
C
C     Outdated options (or options not yet implemented 
C                       for distributed processing!)
C
      IF(ISPODF.GE.1) THEN
         IF(IFPREC.EQ.0) 
     *   CALL QUIT('inconsistent ispodf and ipfrec',ispodf,ifprec) 
      END IF
C
      RETURN
      END
C
C
C     ****************************************************************
C
C
      SUBROUTINE GETWRD(TEXT,K0,K1,K2)
C
C  FINDS NEXT WORD IN TEXT FROM INDEX K0. NEXT WORD IS TEXT(K1:K2)
C  THE NEXT WORD STARTS AT THE FIRST ALPHANUMERIC CHARACTER AT K0
C  OR AFTER. IT ENDS WITH THE LAST ALPHANUMERIC CHARACTER IN A ROW
C  FROM THE START
C
C  TAKEN FROM MULTI - M. CARLSSON (1976)
C
C     INCLUDE 'IMPLIC.FOR'
      PARAMETER (MSEPAR=7)
      CHARACTER*(*) TEXT
      CHARACTER SEPAR(MSEPAR)
      DATA SEPAR/' ','(',')','=','*','/',','/
C
      K1=0
      DO 400 I=K0,LEN(TEXT)
        IF(K1.EQ.0) THEN
          DO 100 J=1,MSEPAR
            IF(TEXT(I:I).EQ.SEPAR(J)) GOTO 200
  100     CONTINUE
          K1=I
C
C  NOT START OF WORD
C
  200     CONTINUE
        ELSE
          DO 300 J=1,MSEPAR
            IF(TEXT(I:I).EQ.SEPAR(J)) GOTO 500
  300     CONTINUE
        ENDIF
  400 CONTINUE
C
C  NO NEW WORD. RETURN K1=K2=0
C
      K1=0
      K2=0
      GOTO 999
C
C  NEW WORD IN TEXT(K1:I-1)
C
  500 CONTINUE
      K2=I-1
C
  999 CONTINUE
      RETURN
      END
C
C
C     ****************************************************************
C
C
      SUBROUTINE STATE(MODE,ID,T,ANE)
C     ===============================
C
C     For MODE=0 - initialization of the basic parameters for the
C                  chemical species
C     for MODE=2 - solves the set of LTE Saha equations to determine
C                  the total charge due to ionization of non-explicit
C                  chemical species
C     for MODE=3 - similar as MODE=2, but also with derivatives wrt
C                  temperature and electron density (called from BPOP,
C                  ie. from the linearization step)
C     for MODE=1 - similar as MODE=2, but the total charge is evaluated
C                  summing the contributions of both non-explicit and
C                  explicit chemical species (called from LTEGR)
C
C     Input for MODE > 0:
C
C        T   - temperature
C        ANE - electron density
C
C     Output for MODE > 0 (through COMMON/STATEP)
C
C
C        Q    - total charge, relative to the reference atom
C        QM   - charge of H- (in LTE), evaluated only if H- is not
C               explicit ion, and if desired (if IHM=1)
C        DQT  - derivative of Q wrt temperature
C        DQN  - derivative of Q wrt electron density
C        DQM  - derivative of QM wrt temperature
C        ENER - internal energy (of all species in all considered
C               ionization stages (needed only for evaluating
C               thermodynamic derivatoves if convection is considerd)
C        QREF - total charge due to the reference species
C        DQTR - derivative of QREF wrt temperature
C        DQNR - derivative of QREF wrt electron density
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'ATOMIC.FOR'
      INCLUDE 'MODELQ.FOR'
      CHARACTER*4 TYPAT,DYP
      COMMON/PFSTDS/PFSTD(99,30),MODPF(matom)
      DIMENSION TYPAT(matom),ABND(matom),D(3,matom),
     *          abnref(mdepth),DYP(matom),
     *          xio(8,matom)
      dimension xio2(9,22), xio3(9,13)
      dimension ffi(matom),pfstu(matom),pfstt(matom),pfstn(matom),
     *          entot(matom)
      dimension idat(30),uu(30,17)
      dimension u6(6),u7(7),u8(8),u10(10),u11(11)
      dimension u12(12),u13(13),u14(14),u16(16),u18(18),u20(20)
      dimension u24(24),u25(25),u26(26),u28(28)
      PARAMETER (TH0=5.0404D3,XMX0=2.154D4,THL0=2.3025851,FI0=3.6113D1,
     *           TRHA=1.5D0,c1qm=1.0353d-16,c2qm=8762.9,
     *           ev2erg=1.6018d-12)
C
      equivalence (u6(1),uu(1,3)),(u7(1),uu(1,4)),(u8(1),uu(1,5))
      equivalence (u10(1),uu(1,6)),(u11(1),uu(1,7)),(u12(1),uu(1,8))
      equivalence (u13(1),uu(1,9)),(u14(1),uu(1,10)),(u16(1),uu(1,11))
      equivalence (u18(1),uu(1,12)),(u20(1),uu(1,13)),(u24(1),uu(1,14))
      equivalence (u25(1),uu(1,15)),(u26(1),uu(1,16)),(u28(1),uu(1,17))
C
      DATA DYP/' H  ',' He ',' Li ',' Be ',' B  ',' C  ',
     *         ' N  ',' O  ',' F  ',' Ne ',' Ne ',' Mg ',
     *         ' Al ',' Si ',' P  ',' S  ',' Cl ',' Ar ',
     *         ' K  ',' Ca ',' Sc ',' Ti ',' V  ',' Cr ',
     *         ' Mn ',' Fe ',' Co ',' Ni ',' Cu ',' Zn ',
     *         ' Ga ',' Ge ',' As ',' Se ',' Br ',' Kr ',
     *         ' Rb ',' Sr ',' Y  ',' Zr ',' Nb ',' Mo ',
     *         ' Tc ',' Ru ',' Rh ',' Pd ',' Ag ',' Cd ', 
     *         ' In ',' Sn ',' Sb ',' Te ',' I  ',' Xe ',
     *         ' Cs ',' Ba ',' La ',' Ce ',' Pr ',' Nd ', 
     *         ' Pm ',' Sm ',' Eu ',' Gd ',' Tb ',' Dy ', 
     *         ' Ho ',' Er ',' Tm ',' Yb ',' Lu ',' Hf ', 
     *         ' Ta ',' W  ',' Re ',' Os ',' Ir ',' Pt ', 
     *         ' Au ',' Hg ',' Tl ',' Pb ',' Bi ',' Po ', 
     *         ' At ',' Rn ',' Fr ',' Ra ',' Ac ',' Th ', 
     *         ' Pa ',' U  ',' Np ',' Pu ',' Am ',' Cm ', 
     *         ' Bk ',' Cf ',' Es '/
C
C    Standard atomic constants for first 99 species 
C      Abundances for the first 30 from Grevesse & Sauval,
C         (1998, Space Sci. Rev. 85, 161)
C
C            Element Atomic  Solar    Std.
C                    weight abundance highest 
C 
C                                     ionization stage 
      DATA D/ 1.008, 1.0D0, 2.,
     *        4.003, 1.00D-1, 3.,
     *        6.941, 1.26D-11, 3.,
     *        9.012, 2.51D-11, 3.,
     *       10.810, 5.0D-10, 4.,
     *       12.011, 3.31D-4, 5.,
     *       14.007, 8.32D-5, 5.,
     *       16.000, 6.76D-4, 5.,
     *       18.918, 3.16D-8, 4.,
     *       20.179, 1.20D-4, 4.,
     *       22.990, 2.14D-6, 4.,
     *       24.305, 3.80D-5, 4.,
     *       26.982, 2.95D-6, 4.,
     *       28.086, 3.55D-5, 5.,
     *       30.974, 2.82D-7, 5.,
     *       32.060, 2.14D-5, 5.,
     *       35.453, 3.16D-7, 5.,
     *       39.948, 2.52D-6, 5.,
     *       39.098, 1.32D-7, 5.,
     *       40.080, 2.29D-6, 5.,
     *       44.956, 1.48D-9, 5.,
     *       47.900, 1.05D-7, 5.,
     *       50.941, 1.00D-8, 5.,
     *       51.996, 4.68D-7, 5.,
     *       54.938, 2.45D-7, 5.,
     *       55.847, 3.16D-5, 5.,
     *       58.933, 8.32D-8, 5.,
     *       58.700, 1.78D-6, 5.,
     *       63.546, 1.62D-8, 5.,
     *       65.380, 3.98D-8, 5.,
     *       69.72 ,   1.34896324e-09  ,  3.,  
     *       72.60 ,   4.26579633e-09  ,  3.,  
     *       74.92 ,   2.34422821e-10  ,  3.,  
     *       78.96 ,   2.23872066e-09  ,  3.,  
     *       79.91 ,   4.26579633e-10  ,  3.,  
     *       83.80 ,   1.69824373e-09  ,  3.,  
     *       85.48 ,   2.51188699e-10  ,  3.,  
     *       87.63 ,   8.51138173e-10  ,  3.,  
     *       88.91 ,   1.65958702e-10  ,  3.,  
     *       91.22 ,   4.07380181e-10  ,  3.,  
     *       92.91 ,   2.51188630e-11  ,  3.,   
     *       95.95 ,   9.12010923e-11  ,  3.,   
     *       99.00 ,   1.00000000e-24  ,  3.,   
     *       101.1 ,   6.60693531e-11  ,  3.,   
     *       102.9 ,   1.23026887e-11  ,  3.,   
     *       106.4 ,   5.01187291e-11  ,  3.,   
     *       107.9 ,   1.73780087e-11  ,  3.,   
     *       112.4 ,   5.75439927e-11  ,  3.,   
     *       114.8 ,   6.60693440e-12  ,  3.,   
     *       118.7 ,   1.38038460e-10  ,  3.,   
     *       121.8 ,   1.09647810e-11  ,  3.,   
     *       127.6 ,   1.73780087e-10  ,  3.,   
     *       126.9 ,   3.23593651e-11  ,  3.,   
     *       131.3 ,   1.69824373e-10  ,  3.,   
     *       132.9 ,   1.31825676e-11  ,  3.,   
     *       137.4 ,   1.62181025e-10  ,  3.,   
     *       138.9 ,   1.58489337e-11  ,  3.,   
     *       140.1 ,   4.07380293e-11  ,  3.,   
     *       140.9 ,   6.02559549e-12  ,  3.,   
     *       144.3 ,   2.95120943e-11  ,  3.,   
     *       147.0 ,   1.00000000e-24  ,  3.,   
     *       150.4 ,   9.33254366e-12  ,  3.,   
     *       152.0 ,   3.46736869e-12  ,  3.,   
     *       157.3 ,   1.17489770e-11  ,  3.,   
     *       158.9 ,   2.13796216e-12  ,  3.,   
     *       162.5 ,   1.41253747e-11  ,  3.,   
     *       164.9 ,   3.16227767e-12  ,  3.,   
     *       167.3 ,   8.91250917e-12  ,  3.,   
     *       168.9 ,   1.34896287e-12  ,  3.,   
     *       173.0 ,   8.91250917e-12  ,  3.,   
     *       175.0 ,   1.31825674e-12  ,  3.,   
     *       178.5 ,   5.37031822e-12  ,  3.,   
     *       181.0 ,   1.34896287e-12  ,  3.,   
     *       183.9 ,   4.78630102e-12  ,  3.,   
     *       186.3 ,   1.86208719e-12  ,  3.,   
     *       190.2 ,   2.39883290e-11  ,  3.,   
     *       192.2 ,   2.34422885e-11  ,  3.,   
     *       195.1 ,   4.78630036e-11  ,  3.,   
     *       197.0 ,   6.76082952e-12  ,  3.,   
     *       200.6 ,   1.23026887e-11  ,  3.,   
     *       204.4 ,   6.60693440e-12  ,  3.,   
     *       207.2 ,   1.12201834e-10  ,  3.,   
     *       209.0 ,   5.12861361e-12  ,  3.,   
     *       210.0 ,   1.00000000e-24  ,  3.,   
     *       211.0 ,   1.00000000e-24  ,  3.,   
     *       222.0 ,   1.00000000e-24  ,  3.,   
     *       223.0 ,   1.00000000e-24  ,  3.,   
     *       226.1 ,   1.00000000e-24  ,  3.,   
     *       227.1 ,   1.00000000e-24  ,  3.,   
     *       232.0 ,   1.20226443e-12  ,  3.,   
     *       231.0 ,   1.00000000e-24  ,  3.,  
     *       238.0 ,   3.23593651e-13  ,  3.,  
     *       237.0 ,   1.00000000e-24  ,  3.,  
     *       244.0 ,   1.00000000e-24  ,  3.,  
     *       243.0 ,   1.00000000e-24  ,  3.,  
     *       247.0 ,   1.00000000e-24  ,  3.,  
     *       247.0 ,   1.00000000e-24  ,  3.,  
     *       251.0 ,   1.00000000e-24  ,  3.,  
     *       254.0 ,   1.00000000e-24  ,  3./
C
C
C     Ionization potentials for first 99 species:
      DATA XIo/
C
C     Element Ionization potentials (eV) 
C              I     II      III     IV       V     VI     VII    VIII
C
     *       13.595,  0.   ,  0.   ,  0.   ,  0.  ,  0.  ,  0.  ,  0.  ,
     *       24.580, 54.400,  0.   ,  0.   ,  0.  ,  0.  ,  0.  ,  0.  ,
     *        5.392, 75.619,122.451,  0.   ,  0.  ,  0.  ,  0.  ,  0.  ,
     *        9.322, 18.206,153.850,217.713,  0.  ,  0.  ,  0.  ,  0.  ,
     *        8.296, 25.149, 37.920,259.298,340.22,  0.  ,  0.  ,  0.  ,
     *       11.264, 24.376, 47.864, 64.476,391.99,489.98,  0.  ,  0.  ,
     *       14.530, 29.593, 47.426, 77.450, 97.86,551.93,667.03,  0.  ,
     *       13.614, 35.108, 54.886, 77.394,113.87,138.08,739.11,871.39,
     *       17.418, 34.980, 62.646, 87.140,114.21,157.12,185.14,953.6 ,
     *       21.559, 41.070, 63.500, 97.020,126.30,157.91,207.21,239.0 ,
     *        5.138, 47.290, 71.650, 98.880,138.37,172.09,208.44,264.16,
     *        7.664, 15.030, 80.120,102.290,141.23,186.49,224.9 ,265.96, 
     *        5.984, 18.823, 28.440,119.960,153.77,190.42,241.38,284.53, 
     *        8.151, 16.350, 33.460, 45.140,166.73,205.11,246.41,303.07, 
     *       10.484, 19.720, 30.156, 51.354, 65.01,220.41,263.31,309.26,
     *       10.357, 23.400, 35.000, 47.290, 72.50, 88.03,280.99,328.8 ,
     *       12.970, 23.800, 39.900, 53.500, 67.80, 96.7 ,114.27,348.3 ,
     *       15.755, 27.620, 40.900, 59.790, 75.00, 91.3 ,124.0 ,143.46,
     *        4.339, 31.810, 46.000, 60.900, 82.6 , 99.7 ,118.0 ,155.0 ,
     *        6.111, 11.870, 51.210, 67.700, 84.39,109.0 ,128.0 ,147.0 ,
     *        6.560, 12.890, 24.750, 73.900, 92.0 ,111.1 ,138.0 ,158.7 ,
     *        6.830, 13.630, 28.140, 43.240, 99.8 ,120.0 ,140.8 ,168.5 ,
     *        6.740, 14.200, 29.700, 48.000, 65.2 ,128.9 ,151.0 ,173.7 ,
     *        6.763, 16.490, 30.950, 49.600, 73.0 , 90.6 ,161.1 ,184.7 ,
     *        7.432, 15.640, 33.690, 53.000, 76.0 , 97.0 ,119.24,196.46,
     *        7.870, 16.183, 30.652, 54.800, 75.0 , 99.1 ,125.0 ,151.06,
     *        7.860, 17.060, 33.490, 51.300, 79.5 ,102.0 ,129.0 ,157.0 ,
     *        7.635, 18.168, 35.170, 54.900, 75.5 ,108.0 ,133.0 ,162.0 ,
     *        7.726, 20.292, 36.830, 55.200, 79.9 ,103.0 ,139.0 ,166.0 ,
     *        9.394, 17.964, 39.722, 59.400, 82.6 ,108.0 ,134.0 ,174.0 ,
     *        6.000,  20.509,   30.700, 99.99,99.99,99.99,99.99,99.99,  
     *        7.89944,15.93462, 34.058, 45.715,99.99,99.99,99.99,99.99,    
     *        9.7887, 18.5892,  28.351, 99.99,99.99,99.99,99.99,99.99,    
     *        9.750,21.500, 32.000, 99.99,99.99,99.99,99.99,99.99,    
     *       11.839,21.600, 35.900, 99.99,99.99,99.99,99.99,99.99,    
     *       13.995,24.559, 36.900, 99.99,99.99,99.99,99.99,99.99,    
     *        4.175,27.500, 40.000, 99.99,99.99,99.99,99.99,99.99,    
     *        5.692,11.026, 43.000, 99.99,99.99,99.99,99.99,99.99,    
     *        6.2171,12.2236, 20.5244,60.607,99.99,99.99,99.99,99.99,    
     *        6.63390,13.13,23.17,34.418,80.348,99.99,99.99,99.99,    
     *        6.879,14.319, 25.039, 99.99,99.99,99.99,99.99,99.99,   
     *        7.099,16.149, 27.149, 99.99,99.99,99.99,99.99,99.99,   
     *        7.280,15.259, 30.000, 99.99,99.99,99.99,99.99,99.99,   
     *        7.364,16.759, 28.460, 99.99,99.99,99.99,99.99,99.99,   
     *        7.460,18.070, 31.049, 99.99,99.99,99.99,99.99,99.99,   
     *        8.329,19.419, 32.920, 99.99,99.99,99.99,99.99,99.99,   
     *        7.574,21.480, 34.819, 99.99,99.99,99.99,99.99,99.99,   
     *        8.990,16.903, 37.470, 99.99,99.99,99.99,99.99,99.99,   
     *        5.784,18.860, 28.029, 99.99,99.99,99.99,99.99,99.99,   
     *        7.342,14.627, 30.490,72.3,99.99,99.99,99.99,99.99,   
     *        8.639,16.500, 25.299,44.2,55.7,99.99,99.99,99.99,   
     *        9.0096,18.600, 27.96, 37.4,58.7,99.99,99.99,99.99,   
     *       10.454,19.090, 32.000, 99.99,99.99,99.99,99.99,99.99,   
     *       12.12984,20.975,31.05,45.,54.14,99.99,99.99,99.99,   
     *        3.893,25.100, 35.000, 99.99,99.99,99.99,99.99,99.99,   
     *        5.210,10.000, 37.000, 99.99,99.99,99.99,99.99,99.99,   
     *        5.580,11.060, 19.169, 99.99,99.99,99.99,99.99,99.99,   
     *        5.650,10.850, 20.080, 99.99,99.99,99.99,99.99,99.99,   
     *        5.419,10.550, 23.200, 99.99,99.99,99.99,99.99,99.99,   
     *        5.490,10.730, 20.000, 99.99,99.99,99.99,99.99,99.99,   
     *        5.550,10.899, 20.000, 99.99,99.99,99.99,99.99,99.99,   
     *        5.629,11.069, 20.000, 99.99,99.99,99.99,99.99,99.99,   
     *        5.680,11.250, 20.000, 99.99,99.99,99.99,99.99,99.99,   
     *        6.159,12.100, 20.000, 99.99,99.99,99.99,99.99,99.99,   
     *        5.849,11.519, 20.000, 99.99,99.99,99.99,99.99,99.99,   
     *        5.930,11.670, 20.000, 99.99,99.99,99.99,99.99,99.99,   
     *        6.020,11.800, 20.000, 99.99,99.99,99.99,99.99,99.99,   
     *        6.099,11.930, 20.000, 99.99,99.99,99.99,99.99,99.99,   
     *        6.180,12.050, 23.700, 99.99,99.99,99.99,99.99,99.99,   
     *        6.250,12.170, 20.000, 99.99,99.99,99.99,99.99,99.99,   
     *        6.099,13.899, 19.000, 99.99,99.99,99.99,99.99,99.99,   
     *        7.000,14.899, 23.299, 99.99,99.99,99.99,99.99,99.99,   
     *        7.879,16.200, 24.000, 99.99,99.99,99.99,99.99,99.99,   
     *        7.86404,17.700, 25.000, 99.99,99.99,99.99,99.99,99.99,   
     *        7.870,16.600, 26.000, 99.99,99.99,99.99,99.99,99.99,   
     *        8.500,17.000, 27.000, 99.99,99.99,99.99,99.99,99.99,   
     *        9.100,20.000, 28.000, 99.99,99.99,99.99,99.99,99.99,   
     *        8.95868,18.563,33.227, 99.99,99.99,99.99,99.99,99.99,   
     *        9.220,20.500, 30.000, 99.99,99.99,99.99,99.99,99.99,   
     *       10.430,18.750, 34.200, 99.99,99.99,99.99,99.99,99.99,   
     *        6.10829,20.4283,29.852,50.72,99.99,99.99,99.99,99.99,   
     *        7.416684,15.0325,31.9373,42.33,69.,99.99,99.99,99.99,   
     *        7.285519,16.679, 25.563,45.32,56.0,88.,99.99,99.99,   
     *        8.430,19.000, 27.000, 99.99,99.99,99.99,99.99,99.99,   
     *        9.300,20.000, 29.000, 99.99,99.99,99.99,99.99,99.99,   
     *       10.745,20.000, 30.000, 99.99,99.99,99.99,99.99,99.99,   
     *        4.000,22.000, 33.000, 99.99,99.99,99.99,99.99,99.99,   
     *        5.276,10.144, 34.000, 99.99,99.99,99.99,99.99,99.99,   
     *        6.900,12.100, 20.000, 99.99,99.99,99.99,99.99,99.99,   
     *        6.000,12.000, 20.000, 99.99,99.99,99.99,99.99,99.99,   
     *        6.000,12.000, 20.000, 99.99,99.99,99.99,99.99,99.99,    
     *        6.000,12.000, 20.000, 99.99,99.99,99.99,99.99,99.99,    
     *        6.000,12.000, 20.000, 99.99,99.99,99.99,99.99,99.99,    
     *        6.000,12.000, 20.000, 99.99,99.99,99.99,99.99,99.99,    
     *        6.000,12.000, 20.000, 99.99,99.99,99.99,99.99,99.99,    
     *        6.000,12.000, 20.000, 99.99,99.99,99.99,99.99,99.99,    
     *        6.000,12.000, 20.000, 99.99,99.99,99.99,99.99,99.99,    
     *        6.000,12.000, 20.000, 99.99,99.99,99.99,99.99,99.99,     
     *        6.000,12.000, 20.000, 99.99,99.99,99.99,99.99,99.99/
C
c     additional ionization potentials
c            ix      x     xi   xii   xiii  xiv   xv   xvi    xvii
c
c    energies added for Sc, Ti, V, Cr, Co,Ni, Cy,Zn
c    source G. C. Rodrigues et al. ,
c   Systematic calculation of total atomic energies of ground state configurations,
c   Atomic Data and Nuclear Data Tables, Volume 86, Issue 2, March 2004, Pages 117-233.

      data xio2/
     *      1103.,    0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,      ! F
     *      1196., 1362.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,      ! Ne
     *       300., 1465.,1649.,   0.,   0.,   0.,   0.,   0.,   0.,      ! Na
     *       328.,  367.,1762.,1963.,   0.,   0.,   0.,   0.,   0.,      ! Mg
     *       330.,  398., 442.,2085.,2304.,   0.,   0.,   0.,   0.,      ! Al
     *       351.,  401., 476., 523.,2438.,2673.,   0.,   0.,   0.,      ! Si
     *       372.,  424., 480., 560., 612.,2816.,3069.,   0.,   0.,      ! P
     *       379.,  447., 505., 565., 652., 707.,3223.,3494.,   0.,      ! S
     *       400.,  456., 529., 592., 657., 750., 809.,3658.,3946.,      ! Cl
     *       422.,  479., 539., 618., 686., 756., 854., 918.,4121.,      ! Ar
     *       176.,  503., 564., 629., 714., 787., 862., 968.,1034.,      ! K
     *       188.,  211., 591., 656., 726., 817., 895., 974.,1087.,      ! Ca
     *       180.,  225., 250., 686., 755., 830., 926.,1010.,1094.,      ! Sc
     *       193.,  216., 265., 291., 787., 861., 940.,1042.,1132.,      ! Ti
     *       206.,  230., 255., 308., 336., 896., 974.,1060.,1165.,      ! V
     *       209.,  244., 271., 298., 355., 384.,1010.,1095.,1185.,      ! Cr
     *       222.,  248., 286., 314., 344., 404., 435.,1136.,1222.,      ! Mn
     *       235.,  262., 290., 331., 361., 392., 457., 490.,1266.,      ! Fe
     *       186.,  276., 305., 336., 379., 411., 444., 512., 547.,      ! Co
     *       193.,  224., 321., 352., 384., 430., 464., 499., 571.,      ! Ni
     *       199.,  232., 266., 369., 401., 435., 484., 520., 557.,      ! Cu
     *       203.,  238., 274., 311., 420., 454., 490., 542., 579./      ! Zn
c
c     even higher ionization potentials
c            18     19     20    21    22    23    24    25     26
c
      data xio3/
     *      4426.,    0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,      ! Ar
     *      4611., 4934.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,      ! K
     *      1158., 5129.,5470.,   0.,   0.,   0.,   0.,   0.,   0.,      ! Ca
     *      1206., 1288.,5675.,6034.,   0.,   0.,   0.,   0.,   0.,      ! Sc
     *      1222., 1346.,1425.,6249.,   0.,   0.,   0.,   0.,   0.,      ! Ti
     *      1261., 1356.,1480.,1569.,   0.,   0.,   0.,   0.,   0.,      ! V
     *      1294., 1397.,1497.,1627.,   0.,   0.,   0.,   0.,   0.,      ! Cr
     *      1318., 1431.,1540.,1645.,1782.,   0.,   0.,   0.,   0.,      ! Mn
     *      1357., 1459.,1574.,1689.,1799.,1958.,2346.,8828.,9278.,      ! Fe
     *      1396., 1496.,1603.,1723.,1847.,1963.,2112.,   0.,   0.,      ! Co
     *       606., 1538.,1643.,1756.,1880.,2011.,2133.,2288.,   0.,      ! Ni
     *       628.,  670.,1688.,1797.,1915.,2043.,2183.,2310.,2472.,      ! Cu
     *       616.,  693., 737.,1844.,1958.,2082.,2214.,2362.,2494./      ! Zn
C
C
C     data for additional ionization potentials for the Opacity
C     project species (IDAT sets the internal OP indexing)
C
      data idat   / 1, 2, 0, 0, 0, 3, 4, 5, 0, 6,
     *              7, 8, 9,10, 0,11, 0,12, 0,13,
     *              0, 0, 0,14,15,16, 0,17,0,0/ 
      data u10/173.93,330.391,511.8,783.3,1018.,1273.8,1671.792,
     *         1928.462,9645.005,10986.876/
      data u11/41.449,381.395,577.8,797.8,1116.2,1388.5,1681.5,2130.8,
     *         2418.7,11817.061,13297.676/
      data u12/61.671,121.268,646.41,881.1,1139.4,1504.3,1814.3,2144.7,
     *         2645.2,2964.4,14210.261,15829.951/
      data u13/48.278,151.86,229.446,967.8,1239.8,1536.3,1947.3,2295.4,
     *         2663.4,3214.8,3565.6,16825.022,18584.138/
      data u14/65.748,131.838,270.139,364.093,1345.1,1653.9,1988.4,
     *         2445.3,2831.9,3237.8,3839.8,4222.4,19661.693,21560.63/
      data u16/83.558,188.2,280.9,381.541,586.2,710.184,2265.9,2647.4,
     *         3057.7,3606.1,4071.4,4554.3,5255.9,5703.6,26002.663,
     *         28182.535/
      data u18/127.11,222.848,328.6,482.4,605.1,734.04,1002.73,1157.08,
     *         3407.3,3860.9,4347.,4986.6,5533.8,6095.5,6894.2,7404.4,
     *         33237.173,35699.936/
      data u20/49.306,95.752,410.642,542.6,681.6,877.4,1026.,1187.6,
     *         1520.64,1704.047,4774.,5301.,5861.,6595.,7215.,7860.,
     *         8770.,9338.,41366.,44177.41/
      data u24/54.576,132.966,249.7,396.5,560.2,731.02,1291.9,1490.,
     *         1688.,1971.,2184.,2404.,2862.,3098.52,8151.,8850.,
     *         9560.,10480.,11260.,12070.,13180.,13882.,60344.,63675.9/
      data u25/59.959,126.145,271.55,413.,584.,771.1,961.44,1569.,
     *         1789.,2003.,2307.,2536.,2771.,3250.,3509.82,9152.,
     *         9872.,10620.,11590.,12410.,13260.,14420.,15162.,
     *         65660.,69137.4/
      data u26/63.737,130.563,247.22,442.,605.,799.,1008.,1218.38,
     *         1884.,2114.,2341.,2668.,2912.,3163.,3686.,3946.82,
     *         10180.,10985.,11850.,12708.,13620.,14510.,15797.,
     *         16500.,71203.,74829.6/
      data u28/61.6,146.542,283.8,443.,613.5,870.,1070.,1310.,1560.,
     *         1812.,2589.,2840.,3100.,3470.,3740.,4020.,4606.,
     *         4896.2,12430.,13290.,14160.,15280.,16220.,17190.,
     *         18510.,19351.,82984.,86909.4/
C     
      IF(MODE.NE.0) GO TO 50
C
C     For MODE=0, STATE serves as an auxiliary procedure for START
C     Input of basic parameters for individual chemical species
C
C     An element (hydrogen through zinc) can be considered in one of
C     the three following options:
C     1. explicitly - some of energy levels of some of its ionization
C                     states are considered explicitly, ie. their
C                     populations are determined by solving statistical
C                     equilibrium
C     2. implicitly - the atom is assumed not to contribute to
C                     opacity; but is allowed to contribute to the
C                     total number of particles and to the total charge;
C                     the latter is evaluated assuming LTE ionization
C                     balance, ie. by solving a set of Saha equations
C     3. not considered at all
C
C     Input:
C     NATOMS -  the highest atomic number of an element that is
C               considered (explicitly or non-explicitly)
C            <  0 - then NATOMS=-NATOMS, and all Opacity Project
C                   species are treated as with MODPF>0, i.e.
C                   their partition functions are evaluated from
C                   the Opacity Project ionization fractions
C
C     For each element from 1 (hydrogen) to NATOMS, the following
C     parameters:
C
C     MA     =  0  - if the element is not considered (option 3)
C            =  1  - if the element is non-explicit (option 2)
C            =  2  - if the element is explicit (option 1)
C     ABN    -  if ABN=0, solar abundance is assumed (given above;
C                         abundance here is assumed as relative
C                         to hydrogen by number
C               if ABN>0, non-solar abundance ABN is assumed; in an
C                         arbitrary scale
C               if ABN<0, non-solar abundance ABN is assumed;
C                         (-ABN times the solar value)
C               if ABN>1.e6, depth-dpendent abundance (additional input)
C
      READ(IBUFF,*) NATOMS
      if(natoms.eq.0) then
         do id=1,nd
            ytot(id)=1.11
            wmy(id)=1.41
            wmm(id)=2.17e-24
         end do
         iatref=0
         return
      end if
c
c     WRITE(6,600)
      IAT=0
      IREFA=0
      IFOPPF=0
      IF(NATOMS.LT.0) THEN
         NATOMS=-NATOMS
         IFOPPF=2
      END IF
C
      natms=natoms
      if(ifmol.gt.0) natms=30
      DO 20 I=1,NATMS
         LGR(I)=.TRUE.
         LRM(I)=.TRUE.
         IATEX(I)=-1
         if(i.le.natoms) then
            READ(IBUFF,*) MA,ABN,MODPF0
          else
            ma=1
            abn=0.
            modpf0=0
         end if
         IF(MA.EQ.0) GO TO 20
         TYPAT(I)=DYP(I)
         AMAS(I)=D(1,I)
         ABND(I)=D(2,I)
         IONIZ(I)=int(D(3,I))
         MODPF(I)=MODPF0
C
C        increase the standard highest ionization for Teff larger
C        than 50000 K for N, O, Ne, and Fe
C
         IF(TEFF.GT.5.D4) THEN
            IF(I.EQ.7) IONIZ(I)=6
            IF(I.EQ.8) IONIZ(I)=7
            IF(I.EQ.10) IONIZ(I)=9
            IF(I.EQ.26) IONIZ(I)=9
         END IF
C
         DO J=1,8
            ENEV(I,J)=xio(J,I)
         END DO
c
         if(i.ge.9) then
         do j=9,17
            enev(i,j)=xio2(j-8,i-8)
         end do
         endif
         if(i.ge.18) then
         do j=18,30
            enev(i,j)=xio3(j-17,i-17)
         end do
         end if
c
         LGR(I)=.FALSE.
         IF(MODPF(I).GT.0) IFOPPF=1
         IF(IFOPPF.EQ.2.and.idat(i).gt.0) MODPF(I)=1
         if(modpf(i).gt.0) then
            if(idat(i).eq.0) modpf(i)=0
            if(modpf(i).gt.0) then
               ioniz(i)=i+1
               if(i.ge.10) then
                  do j=9,i
                     enev(i,j)=uu(j,idat(i))*0.1239529d0
                  end do
               end if
            end if
         end if
         IF(ABN.GT.0) ABND(I)=ABN
         IF(ABN.LT.0) ABND(I)=ABS(ABN)*D(2,I)
         IF(ION.NE.0) IONIZ(I)=ION
         IF(ABN.GT.1.E6) THEN
            READ(IBUFF,*) (ABNDD(I,ID),ID=1,ND)
          ELSE
            DO ID=1,ND
               ABNDD(I,ID)=ABND(I)
            END DO
         END IF
         IF(MA.EQ.1) THEN
            LRM(I)=.FALSE.
            IATEX(I)=0
          ELSE
            IAT=IAT+1
            IIFIX(IAT)=0
            IF(MA.EQ.-2) IIFIX(IAT)=1
            IATEX(I)=IAT
            IF(IAT.EQ.IATREF) THEN
               IREFA=I
               DO ID=1,ND
                  ABNREF(ID)=ABNDD(I,ID)
               END DO
            END IF
C
C           store parameters for explicit atoms
C
            AMASS(IAT)=AMAS(I)*HMASS
            NUMAT(IAT)=I
         END IF
   20 CONTINUE
C
C     renormalize abundances to have the standard element abundance
C     equal to unity
C
      if(abnref(1).eq.0.) then
         do id=1,nd
            abnref(id)=1.
         end do
      end if
c
      DO 30 I=1,NATOMS
         IAT=IATEX(I)
         IF(IAT.LT.0) GO TO 30
         DO ID=1,ND
            ABNDD(I,ID)=ABNDD(I,ID)/ABNREF(ID)
            YTOT(ID)=YTOT(ID)+ABNDD(I,ID)
            WMY(ID)=WMY(ID)+ABNDD(I,ID)*AMAS(I)
         END DO
         ABNR=ABND(I)/D(2,I)
         IF(IAT.EQ.0) THEN
c           WRITE(6,601) I,TYPAT(I),ABND(I),ABNR
          ELSE
            DO ID=1,ND
               ABUND(IAT,ID)=ABNDD(I,ID)
            END DO 
c           WRITE(6,602) I,TYPAT(I),ABND(I),ABNR,IAT
         END IF
   30 CONTINUE
      DO ID=1,ND
         WMM(ID)=WMY(ID)*HMASS/YTOT(ID)
      END DO
c
c     initialization of the Opacity Project ionization fractions
c     (if required)
c
      if(ifoppf.gt.0) call opfrac(0,0,t,ane,pf,fra)
c
c 600 FORMAT(1H0//' CHEMICAL ELEMENTS INCLUDED'/
c    *            ' --------------------------'//
c    * ' NUMBER  ELEMENT           ABUNDANCE'/1H ,16X,
c    * 'A=N(ELEM)/N(H)  A/A(SOLAR)'/)
c 601 FORMAT(1H ,I4,3X,A5,1P2E14.2)
c 602 FORMAT(1H ,I4,3X,A5,1P2E14.2,3X,
c    *       'EXPLICIT: IAT=',I3)
      RETURN
C
   50 TLN=LOG(T)*TRHA
      TK=BOLK*T
      THET=TH0/T
      THL=THL0*THET
      XMX=XMX0*SQRT(SQRT(T/ANE))
      DCH=EH/(XMX*XMX*TK)
      Q=0.
      QM=0.
      DQT=0.
      DQN=0.
      DQM=0.
      ENER=0.
      DO 70 I=1,NATOMS
         IF(MODE.GT.1.AND.LRM(I).OR.MODE.EQ.1.AND.LGR(I)) GO TO 70
         ION=IONIZ(I)
         DRQT=0.
         DRQN=0.
         DRST=0.
         DRSN=0.
         DFT=0.
         DFN=0.
         ENTOT(1)=0.
         RS=UN
         CALL PARTF(I,1,T,ANE,XMX,UM,DUTM,DUNM)
         if(i.eq.1) pfhyd=um
         pfstu(1)=um
         pfstt(1)=dutm
         pfstn(1)=dunm
         JMAX=1
         DO J=2,ION
            J1=J-1
            DCHT=DCH*J1
            TE=ENEV(I,J1)*THL
            ENTOT(J)=ENTOT(J1)+TE
            dcht=0.
            FI=FI0+TLN-TE+DCHT
            X=J
            XMAX=XMX*SQRT(X)
            CALL PARTF(I,J,T,ANE,XMAX,U,DUT,DUN)
            pfstu(j)=u
            pfstt(j)=dut
            pfstn(j)=dun
            FFI(J)=0.
            IF(FI.GT.-20.) FFI(J)=EXP(FI)*U/UM/ANE
            IF(FFI(J).GT.UN) JMAX=J
            UM=U
         END DO
         RQ=JMAX-1
         RI=ENTOT(JMAX)
         RE=PFSTT(JMAX)/PFSTU(JMAX)*T
         if(jmax.lt.ion) then
            R=UN
            DO J=JMAX+1,ION
               J1=J-1
               DCHT=DCH*J1
               TE=ENEV(I,J1)*THL
               R=R*FFI(J)
               RR(I,J)=R/pfstu(j)
               RR(I,J)=R
               RS=RS+R
               RQ=RQ+J1*R
               RI=RI+R*ENTOT(J)
               RE=RE+R*PFSTT(J)/PFSTU(J)*T
               DFIT=pfstt(j)/pfstu(j)-pfstt(j1)/pfstu(j1)
     .              +(TRHA+TE-TRHA*DCHT)/T
               DFIN=pfstn(j)/pfstu(j)-pfstn(j1)/pfstu(j1)
     .              +(HALF*DCHT-UN)/ANE
               DFT=DFT+DFIT
               DFN=DFN+DFIN
               DFIT=DFT*R
               DFIN=DFN*R
               DRST=DRST+DFIT
               DRSN=DRSN+DFIN
               DRQT=DRQT+J1*DFIT
               DRQN=DRQN+J1*DFIN
            END DO
         end if
         if(jmax.gt.1) then
            R=UN
            DFT=0.
            DFN=0.
            jmin=min(4,jmax-1)
            DO JJ=1,JMIN
               J=JMAX-JJ
               J1=J-1
               JP1=J+1
               DCHT=DCH*J
               TE=ENEV(I,J)*THL
               R=R/FFI(JP1)
C              RR(I,J)=R/pfstu(j)
               RR(I,J)=R
               RS=RS+R
               RQ=RQ+J1*R
               RI=RI+R*ENTOT(J)
               RE=RE+R*PFSTT(J)/PFSTU(J)*T
               DFIT=pfstt(jp1)/pfstu(jp1)-pfstt(j)/pfstu(j)
     *              +(TRHA+TE-TRHA*DCHT)/T
               DFIN=pfstn(jp1)/pfstu(jp1)-pfstn(j)/pfstu(j)
     *              +(HALF*DCHT-UN)/ANE
               DFT=DFT-DFIT
               DFN=DFN-DFIN
               DFIT=DFT*R
               DFIN=DFN*R
               DRST=DRST+DFIT
               DRSN=DRSN+DFIN
               DRQT=DRQT+J1*DFIT
               DRQN=DRQN+J1*DFIN
            END DO
         endif
         X=RQ/RS
         ABND(I)=ABNDD(I,ID)
         X1=ABND(I)/RS
         RR(I,JMAX)=X1
         DO J=1,ION
            IF(J.NE.JMAX) RR(I,J)=RR(I,J)*X1
         END DO
C         RR(I,JMAX)=RR(I,JMAX)/PFSTU(JMAX)
c         if(i.ge.2) then
c            if(iftene.eq.1) then 
c               ener=ener+ri*x1
c             else
c               ENER=ENER+(RI+RE)*X1
c            end if
c         end if
c
c        internal energy
c
         do j=2,ion
            ener=ener+enev(i,j-1)*rr(i,j)*ev2erg
         end do
         if(iftene.gt.1) then
            do j=1,ion
               ener=ener+rr(i,j)*pfstt(j)/pfstu(j)*t*tk
            end do
         end if
c   
         aref=dens(id)/wmm(id)/ytot(id)
c
         IF(I.EQ.IREFA) THEN
            QREF=X*ABND(I)
            DQTR=(DRQT-X*DRST)*X1
            DQNR=(DRQN-X*DRSN)*X1
          ELSE
            Q=X*ABND(I)+Q
            DQT=DQT+(DRQT-X*DRST)*X1
            DQN=DQN+(DRQN-X*DRSN)*X1
         END IF
   70 CONTINUE
C
C     Negative hydrogen ion
C
      IF(IHM.EQ.1) THEN
         tinv=un/t
         QM=C1QM*Tinv/SQRT(T)*EXP(C2QM*Tinv)
         DQM=-QM*Tinv*(TRHA+C2QM*Tinv)
      END IF
      RETURN
      END
C
C
C     ****************************************************************
C
      SUBROUTINE INPMOD(ifinpm)
C     =========================
C
C     Read an initial model atmosphere from unit 8
C     File 8 contains:
C      1. NDPTH -  number of depth points in which the initial model is
C                  given (if not equal to ND, routine interpolates
C                  automatically to the set DM by linear interpolation
C                  in log(DM)
C         NUMPAR - number of input model parameters in each depth
C                  = 3 for LTE model - ie. N, T, N(electron);
C                  > 3 for NLTE model)
C      2. DEPTH(ID),ID=1,NDPTH - mass-depth points for the input model
C      3. for each depth:
C                 T   - temperature
C                 ANE - electron density
C                 RHO - mass density
C                 level populations - only for NLTE input model
C                       Number of input level populations need not be
C                       equal to NLEVEL; in that case the procedure
C                       CHANGE is called from START to calculate the
C                       remaining level populations
C
C     Note: The output file 7, which is created by this program
C           (procedure OUTPUT) has the same structure as file 8
C           and may thus be used as input to another run of the
C           program
C
C     INTRPL - switch indicating whether (and, if so, how) interpolate
C              the initial model if the depth scales for the input model
C              and the present depth scale are different
C            = 0  -  no interpolation, scale DEPTH replaces DM;
C                    i.e. DEPTH will be used as DM, regardless of which
C                    values of DM were read in subroutine START
C            > 0  -  polynomial interpolation of the (INTRPL-1)th order
C
C            < 0  -  reads different initial models (eg. -1 : Kurucz model)
C
C     If INTRPL > 0, there is an additional input from unit 8, namely
C     new depth scale DM, the one which will be used in the present run
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'ATOMIC.FOR'
      INCLUDE 'MODELQ.FOR'
      PARAMETER (MINPUT=MLEVEL+5,
     *           MDEPTI=MDEPTH)
      COMMON POPUL0(MLEVEL,MDEPTI),ESEMAT(MLEVEL,MLEVEL),BESE(MLEVEL),
     *       TEMP0(MDEPTI),ELEC0(MDEPTI),DENS0(MDEPTI),PPL0(MDEPTI),
     *       DEPTH(MDEPTI),PPL(MDEPTH),POPLTE(MLEVEL),X(MINPUT),
     *       ZD0(MDEPTH)
c     dimension a(mlevel,mlevel),b(mlevel),iifor0(mlevel)
C
      LCHC0=LCHC
      LCHC=.TRUE.
      LTE0=LTE
      LTE=.TRUE.
      NUMLT=3
      IF(IDISK.EQ.1) NUMLT=4
      RRDIL=1.
      TEMPBD=0.
      IF(INTRPL.GE.0) THEN
         READ(8,*,err=50,end=50) NDPTH,NUMPAR
         ifinpm=1
         IF(NDPTH.LE.0) CALL QUIT('NDPTH.LE.0 in Unit 8',ndpth,0)
         IF(NDPTH.GT.MDEPTI) 
     *      CALL QUIT('NDPTH.GT.MDEPTI in Unit 8',ndpth,mdepti)
         READ(8,*) (DEPTH(I),I=1,NDPTH)
         IF(NDPTH.ne.nd.and.intrpl.eq.0) 
     *      CALL QUIT('ndpth.ne.nd in Unit 8',ndpth,nd)
         DO 30 ID=1,NDPTH
            READ(8,*) (X(I),I=1,NUMPAR)
            do i=1,numpar
               if(x(i).lt.0.) x(i)=0.
            end do
            TEMP0(ID)=X(1)
            ELEC0(ID)=X(2)
            DENS0(ID)=X(3)
            IF(IDISK.EQ.1) ZD0(ID)=X(4)
            IF(NUMPAR.GT.NUMLT.AND..NOT.LTE0.AND.ICHANG.NE.-2) THEN
               NLEV0=NUMPAR-NUMLT
               DO 10 I=1,NLEV0
   10            POPUL0(I,ID)=X(NUMLT+I)
             ELSE
               NLEV0=NLEVEL
               TEMP(ID)=X(1)
               ELEC(ID)=X(2)
               DENS(ID)=X(3)
               IF(IDISK.EQ.1) ZD(ID)=X(4)
               CALL WNSTOR(ID)
               CALL SABOLF(ID)
c              DO I=1,NLEV0
c                 IIFOR0(I)=I
c              END DO
c              CALL RATMAT(ID,IIFOR0,-1,A,B)
c              CALL LEVSOL(A,B,POPLTE,IIFOR0,NLEV0,1) 
               DO I=1,NLEV0
                  POPUL0(I,ID)=POPLTE(I)
               END DO
            END IF
   30    CONTINUE
         READ(8,*,END=31,ERR=31) INTRPL
       ELSE IF(INTRPL.GT.-10) THEN
         CALL KURUCZ(NDPTH)
         NUMPAR=3
         IF(ND.NE.NDPTH .AND. INTRPL.EQ.0)
     *      CALL QUIT('ND.NE.NDPTH in KURUCZ',nd,ndpth)
       ELSE
         CALL INCLDY(NDPTH)
         NUMPAR=3
         IF(ND.NE.NDPTH .AND. INTRPL.EQ.0)
     *      CALL QUIT('ND.NE.NDPTH in INCLDY',nd,ndpth)
      END IF
   31 LCHC=LCHC0
      LTE=LTE0
C
C     if INTRPL = 0  -  no interpolation
C        i.e. scale DEPTH replaces DM
C
      IF(INTRPL.EQ.0) THEN
         IDSTD=0
         DO ID=1,ND
            DM(ID)=DEPTH(ID)
            TEMP(ID)=TEMP0(ID)
            ELEC(ID)=ELEC0(ID)
            DENS(ID)=DENS0(ID)
            ANMA(ID)=DENS(ID)/WMM(ID)
            ANTO(ID)=ANMA(ID)+ELEC(ID)
            IF(IDISK.EQ.1) ZD(ID)=ZD0(ID)
            CALL WNSTOR(ID)
            IF(TEMP(ID).LT.TEFF) IDSTD=ID
            DO I=1,NLEVEL
               POPUL(I,ID)=POPUL0(I,ID)
            END DO
         END DO
         ELSTD=ELEC(IDSTD)
         IF(IDISK.EQ.1) THEN
            ZND=ZD(ND)
            IF(ZND.GT.0.) IFZ0=-1
         END IF
         ELSTD=ELEC(IDSTD)
         GO TO 100
       ELSE
C
C    if  INTRPL > 0  - Interpolation from log(DEPTH) to log(DM)
C        INTRPL = 1 or 2 - linear interpolation,
C               = 3  -  parabolic interpolation,
C               = 4  -  cubic interpolation, etc.
C
C        DM(ID)     - mass-depth variable (ie. column mass in g*cm**-2)
C                    at the depth point ID
C
         READ(8,*) (DM(ID),ID=1,ND)
         NTRPL=INTRPL
         CALL INTERP(DEPTH,TEMP0,DM,TEMP,NDPTH,ND,NTRPL,1,0)
         CALL INTERP(DEPTH,ELEC0,DM,ELEC,NDPTH,ND,NTRPL,1,1)
         CALL INTERP(DEPTH,DENS0,DM,DENS,NDPTH,ND,NTRPL,1,1)
         IF(IDISK.EQ.1) 
     *   CALL INTERP(DEPTH,ZD0,DM,ZD,NDPTH,ND,INTRPL,1,0)
         DO I=1,NLEV0
            DO ID=1,NDPTH
               PPL0(ID)=POPUL0(I,ID)
            END DO
            CALL INTERP(DEPTH,PPL0,DM,PPL,NDPTH,ND,NTRPL,1,1)
            DO ID=1,ND
               ANMA(ID)=DENS(ID)/WMM(ID)
               ANTO(ID)=ANMA(ID)+ELEC(ID)
               CALL WNSTOR(ID)
               POPUL(I,ID)=PPL(ID)
            END DO
         END DO
c
c !!!!! attention - temporary fix
c
         if(idisk.eq.1) zd(nd)=0.
C      
         IF(IDISK.EQ.1) THEN
            ZND=ZD(ND)
            IF(ZND.GT.0.) IFZ0=-1
         END IF
      END IF
c
  100 if(ioptab.ge.0) return
      do id=1,nd
         ptotal(id)=dm(id)*grav
         an=ptotal(id)/(bolk*temp(id))
         elec(id)=1.e-16*an
         wmm(id)=dens(id)/an
      end do
      RETURN
c
   50 ifinpm=0
      return
      END
C
C
C     ****************************************************************
C
C
      SUBROUTINE KURUCZ(NDPTH)
C     ========================
C
C     Read an initial model atmosphere from unit 8
C      in Kurucz ATLAS' format
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'ATOMIC.FOR'
      INCLUDE 'MODELQ.FOR'
      PARAMETER (MINPUT=7)
      CHARACTER KUR*15
      DIMENSION KKFIX0(MLEVEL)
      COMMON POPUL0(MLEVEL,MDEPTH),ESEMAT(MLEVEL,MLEVEL),BESE(MLEVEL),
     *       TEMP0(MDEPTH),ELEC0(MDEPTH),DENS0(MDEPTH),PPL0(MDEPTH),
     *       DEPTH(MDEPTH),PPL(MDEPTH),POPLTE(MLEVEL),X(MINPUT)
c     dimension a(mlevel,mlevel),b(mlevel),iifor0(mlevel)
C
      do iat=1,natom
         kkfix0(iat)=iifix(iat)
         iifix(iat)=0
      end do
C
      READ(8,801) KUR,GRAVK
      READ(KUR,802) TEFFK
  801 FORMAT(A15,6X,F8.5)
  802 FORMAT(4X,F8.0)
C
      IF(KUR(1:4).NE.'TEFF')
     *  CALL QUIT(' Unit 8 is NOT a Kurucz model as expected',0,0)
      IF(ABS(TEFFK-TEFF).GT.50.) then
        ieff=int(teff)
        ieffk=int(teffk)
        CALL QUIT(' Teff not corresponding to Kurucz model',ieff,ieffk)
      END IF
      IF(ABS(GRAVK-LOG10(GRAV)).GT.0.02) then
        irav=int(log10(grav)+0.001)
        iravk=int(gravk)
        CALL QUIT(' Gravity not corresponding to Kurucz model',
     *             irav,iravk)
      END IF
C
      DO WHILE(KUR(1:9).NE.'READ DECK')
         READ(8,'(A15)') KUR
      END DO
      READ(KUR,803) NDPTH
  803 FORMAT(10X,I3)
      NDPTH=NDPTH-1
      NUMPAR=3
      NLEV0=NLEVEL
      READ(8,*) TTT
      IF(NDPTH.gt.mdepth.and.intrpl.eq.0) 
     *   CALL QUIT('ndpth.gt.mdepth in KURUCZ',ndpth,mdepth)
      DO ID=1,NDPTH
         READ(8,*) (X(I),I=1,MINPUT)
         DEPTH(ID)=X(1)
         TEMP0(ID)=X(2)
         ELEC0(ID)=X(4)
         AN=X(3)/BOLK/TEMP0(ID)
         DENS0(ID)=WMM(ID)*(AN-ELEC0(ID))
         TEMP(ID)=TEMP0(ID)
         ELEC(ID)=ELEC0(ID)
         DENS(ID)=DENS0(ID)
         ANMA(ID)=DENS(ID)/WMM(ID)
         ANTO(ID)=ANMA(ID)+ELEC(ID)
         CALL WNSTOR(ID)
         CALL SABOLF(ID)
c        DO I=1,NLEV0
c           IIFOR0(I)=I
c        END DO
c        CALL RATMAT(ID,IIFOR0,-1,A,B)
c        CALL LEVSOL(A,B,POPLTE,IIFOR0,NLEV0,1) 
         DO I=1,NLEV0
            POPUL0(I,ID)=POPLTE(I)
         END DO
      END DO
C
      INTRPL=0
      DO WHILE(KUR(1:5).NE.'BEGIN')
         READ(8,'(A15)',END=100,ERR=100) KUR
      END DO
      READ(8,*,END=100,ERR=100) INTRPL
C
  100 do iat=1,natom
         iifix(iat)=kkfix0(iat)
      end do
C
      RETURN
      END
C
C
C     ****************************************************************
C
C
      SUBROUTINE INCLDY(NDPTH)
C     ========================
C
C     Read an initial model atmosphere from unit 8
C      in Cloudy's format as provided by Katya Verner
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'ATOMIC.FOR'
      INCLUDE 'MODELQ.FOR'
      PARAMETER (MINPUT=6)
      DIMENSION RS(MDEPTH),KKFIX0(MLEVEL)
      COMMON POPUL0(MLEVEL,MDEPTH),ESEMAT(MLEVEL,MLEVEL),BESE(MLEVEL),
     *       TEMP0(MDEPTH),ELEC0(MDEPTH),DENS0(MDEPTH),PPL0(MDEPTH),
     *       DEPTH(MDEPTH),PPL(MDEPTH),POPLTE(MLEVEL),X(MINPUT)
c     dimension a(mlevel,mlevel),b(mlevel),iifor0(mlevel)
C
      do iat=1,natom
         kkfix0(iat)=iifix(iat)
         iifix(iat)=0
      end do
C
      READ(8,*) NDPTH
      NUMPAR=3
      IF(NDPTH.gt.mdepth.and.intrpl.eq.0) 
     *   CALL QUIT('ndpth.gt.mdepth in INCLDY',ndpth,mdepth)
      NLEV0=NLEVEL
      READ(8,*) TSTARY,RSTARY
      IF(RSTARY.LT.1.E6) RSTARY=RSTARY*6.96E10
C
      DO ID=NDPTH,1,-1
         READ(8,*) (X(I),I=1,MINPUT)
         RS(ID)=X(1)
         TEMP0(ID)=X(2)
         ELEC0(ID)=X(4)
         AN=X(3)/BOLK/TEMP0(ID)
         DENS0(ID)=WMM(ID)*X(3)
         TEMP(ID)=TEMP0(ID)
         ELEC(ID)=ELEC0(ID)
         DENS(ID)=DENS0(ID)
         ANMA(ID)=DENS(ID)/WMM(ID)
         ANTO(ID)=ANMA(ID)+ELEC(ID)
      END DO
      DM(1)=0.
      DEPTH(1)=0.
      DO ID=2,NDPTH
         DDDM=(DENS(ID-1)+DENS(ID))*(RS(ID-1)-RS(ID))
      DM(ID)=DM(ID-1)+0.5*DDDM
      DEPTH(ID)=DM(ID)
      END DO
C      
      RRDIL=(RSTARY/RS(NDPTH))*(RSTARY/RS(NDPTH))
      TEMPBD=TSTARY
C
      DO ID=1,NDPTH
         CALL WNSTOR(ID)
         CALL SABOLF(ID)
c        DO I=1,NLEV0
c           IIFOR0(I)=I
c        END DO
c        CALL RATMAT(ID,IIFOR0,-1,A,B)
c        CALL LEVSOL(A,B,POPLTE,IIFOR0,NLEV0,1) 
         DO I=1,NLEV0
            POPUL0(I,ID)=POPLTE(I)
         END DO
      END DO
C
      INTRPL=0
      do iat=1,natom
         iifix(iat)=kkfix0(iat)
      end do
C
      RETURN
      END
C
C
C     ****************************************************************
C
C

      SUBROUTINE CHANGE
C     =================
C
C     This procedure controls an evaluation of initial level
C     populations in case where the system of explicit levels
C     (ie. the choice of explicit level, their numbering, or their
C     total number) is not consistent with that for the input level
C     populations read by procedure INPMOD.
C     Obviously, this procedure need be used only for NLTE input models.
C
C     ICHANG < 0 - general change of populations as described below
C            > 0 - a simplified change; original data for the input
C                    model are required to assign the input NLTE populations
C                    to the levels in the new models; all additional
C                    levels are assumed having LTE populations.
C            ICHANG is the unit number for the data file of old model.
C
C     Case ICHANG < 0:
C
C     Input from unit 5:
C     For each explicit level, II=1,NLEVEL, the following parameters:
C      IOLD   -  NE.0 - means that population of this level is
C                       contained in the set of input populations;
C                       IOLD is then its index in the "old" (i.e. input)
C                       numbering.
C                       All the subsequent parameters have no meaning
C                       in this case.
C             -  EQ.0 - means that this level has no equivalent in the
C                       set of "old" levels. Population of this level
C                       has thus to be evaluated.
C      MODE   -    indicates how the population is evaluated:
C             = 0  - population is equal to the population of the "old"
C                    level with index ISIOLD, multiplied by REL;
C             = 1  - population assumed to be LTE, with respect to the
C                    first state of the next ionization degree whose
C                    population must be contained in the set of "old"
C                    (ie. input) populations, with index NXTOLD in the
C                    "old" numbering.
C                    The population determined of this way may further
C                    be multiplied by REL.
C             = 2  - population determined assuming that the b-factor
C                    (defined as the ratio between the NLTE and
C                    LTE population) is the same as the b-factor of
C                    the level ISINEW (in the present numbering). The
C                    level ISINEW must have the equivalent in the "old"
C                    set; its index in the "old" set is ISIOLD, and the
C                    index of the first state of the next ionization
C                    degree, in the "old" numbering, is NXTSIO.
C                    The population determined of this way may further
C                    be multiplied by REL.
C             = 3  - level corresponds to an ion or atom which was not
C                    explicit in the old system; population is assumed
C                    to be LTE.
C      NXTOLD  -  see above
C      ISINEW  -  see above
C      ISIOLD  -  see above
C      NXTSIO  -  see above
C      REL     -  population multiplier - see above
C                 if REL=0, the program sets up REL=1
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'ATOMIC.FOR'
      INCLUDE 'MODELQ.FOR'
      character*20 fnstd
      dimension n0old(30,30),n1old(30,30)
      dimension katold(2,30),vtbold(mdepth)
      COMMON POPUL0(MLEVEL,MDEPTH),POPULL(MLEVEL,MDEPTH),
     *       ESEMAT(MLEVEL,MLEVEL),BESE(MLEVEL),POPL(MLEVEL)
C
      PARAMETER (S = 2.0706D-16)
      IF(ICHANG.LT.0) THEN
      IFESE=0
      DO 100 II=1,NLEVEL
         READ(IBUFF,*) IOLD,MODE,NXTOLD,ISINEW,ISIOLD,NXTSIO,REL
         IF(REL.EQ.0.) REL=1.
         IF(MODE.GE.3) IFESE=IFESE+1
         DO 90 ID=1,ND
            IF(IOLD.EQ.0) GO TO 10
            POPUL0(II,ID)=POPUL(IOLD,ID)
            GO TO 90
   10       IF(MODE.NE.0) GO TO 20
            POPUL0(II,ID)=POPUL(ISIOLD,ID)*REL
            GO TO 90
   20       T=TEMP(ID)
            ANE=ELEC(ID)
            IF(MODE.GE.3) GO TO 40
            NXTNEW=NNEXT(IEL(II))
            SB=S/T/SQRT(T)*G(II)/G(NXTNEW)*EXP(ENION(II)/T/BOLK)
            IF(MODE.GT.1) GO TO 30
            POPUL0(II,ID)=SB*ANE*POPUL(NXTOLD,ID)*REL
            GO TO 90
   30       KK=ISINEW
            KNEXT=NNEXT(IEL(KK))
            SBK=S/T/SQRT(T)*G(KK)/G(KNEXT)*EXP(ENION(KK)/T/BOLK)
            POPUL0(II,ID)=SB/SBK*POPUL(NXTOLD,ID)/POPUL(NXTSIO,ID)*
     *                    POPUL(ISIOLD,ID)*REL
            GO TO 90
   40       IF(IFESE.EQ.1) THEN
               LTE0=LTE
               LTE=.TRUE.
               do iii=1,nlevel
                  if(wop(iii,id).eq.0.) wop(iii,id)=1.
               end do
c              CALL STEQEQ(ID,POPL,0)
               DO III=1,NLEVEL
                  POPULL(III,ID)=POPL(III)
               END DO
               LTE=LTE0
            END IF
            POPUL0(II,ID)=POPULL(II,ID)
   90    CONTINUE
  100 CONTINUE
      DO 110 I=1,NLEVEL
         DO 110 ID=1,ND
            POPUL(I,ID)=POPUL0(I,ID)
  110 CONTINUE
C
C     simplified change - no additional input (the case ICHANG > 0)
C
      ELSE
       LTE0=LTE
       LTE=.TRUE.
       DO 150 ID=1,ND
         do ii=1,nlevel
         if(wop(ii,id).eq.0.) wop(ii,id)=1.
         end do
c      CALL STEQEQ(ID,POPL,0)
         DO 120 II=1,NLEVEL
            POPUL0(II,ID)=POPL(II)
  120    CONTINUE
  150  CONTINUE
C
       IF(ICHANG.EQ.1) THEN
         DO II=NLEV0+1,NLEVEL
            DO ID=1,ND
               POPUL(II,ID)=POPUL0(II,ID)
            END DO
         END DO
C
       ELSE
         modr=0
         rewind 1
         read(1,*,err=200,end=200) modr
  200    continue
         call readbf(ichang)
         if(modr.eq.0) then
            read(95,*) tfold,grold
            read(95,*) ltd1,ltd2
            read(95,*) fnstd
            read(95,*) nfrd
            if(nfrd.lt.0) then
               nfrd=-nfrd
               do ij=1,nfrd
                  read(95,*) frold
               end do
            endif
            read(95,*) natold
            if(natold.lt.0) natold=-natold
            do ia=1,natold
               read(95,*) iao,abnold
               if(abnold.gt.1.e6) read(95,*) (vtbold(i),i=1,ndold)
            end do
            nlold=0
            read(95,*) iato,izo,nlvo,ilasti,ilvi,instd
            if(instd.ne.0) read(95,*) idui
            do while (ilasti.ge.0)
               n0old(iato,izo+1)=nlold+1
               n1old(iato,izo+1)=nlold+nlvo
               nlold=nlold+nlvo
               read(95,*) iato,izo,nlvo,ilasti,ilvi,instd
               if(instd.ne.0) read(95,*) idui
            end do
          else
            read(95,*) tfold,grold,hmold
            read(95,*) ltd1,ltd2,lcold,ispold,chmold
            if(ispold.lt.0) read(95,*,err=203) iol1,iol2,iol3,iol4,
     .                                       iol5,iol6,iol7
            if(iol6.ge.2) read(95,*) djmold
  203       read(95,*,err=204) nitold,ndold,natold,niold,nlvold,
     .                       iol1,iol2,iarold,iol4
  204       continue
            if(iol1.gt.10) then
               read(95,*,err=205) iol1,iol2,iol3,iol4,iol5
               read(95,*,err=205) iol1,iol2
               read(95,*,err=205) iol1,iol2
            end if
  205       continue
            if(niold.lt.0) then
               niold=-niold
               read(95,*,err=206) iol1,iol2,iol3
            end if
  206       continue
            if(iarold.le.-100 .and. iarold.gt.-200) then
               iarold=-iarold-100
               read(95,*) iol1
            endif
            read(95,*) nfrd
            if(nfrd.gt.0) then
               nfrd=-nfrd
               do ij=1,nfrd
                  read(95,*) frold
               end do
             else
               nfrd=-nfrd
               read(95,*) frold
            end if
            read(95,*,err=211) iol1,iol2,iol3
            if(iol3.lt.0) read(95,*) pzold
  211       continue
            read(95,*) iol1,vtbol
            if(iol1.ne.0) read(95,*) (vtbold(i),i=1,ndold)
            read(95,*) natsold
            if(natsold.lt.0) natsold=-natsold
            iat=0
            do ia=1,natsold
               read(95,*) iol1,iol2,iol3,iol4,iol5,abnold
               if(abnold.gt.1.e6) read(95,*) (vtbold(i),i=1,ndold)
               if(iol1.eq.2) then
                  iat=iat+1
                  katold(1,iat)=iol2
                  katold(2,iat)=iol3
               end if
            end do
            do ii=1,niold
               read(95,*) k0old,k1old,k2old,izo
               if(k0old.lt.0) then
                  k0old=-k0old
                  read(95,*) iol1
               end if
               do ia=1,iat
                  if(k0old.ge.katold(1,ia) .and. k1old.ge.katold(2,ia))
     .            iaol=ia
               end do
               n0old(iaol,izo)=k0old
               n1old(iaol,izo)=k1old
               n0old(iaol,izo+1)=k2old
               n1old(iaol,izo+1)=k2old
            end do
        end if
C
        WRITE(6,600)
  600   FORMAT(' Levels: OLD model -> NEW model',/
     .         ' ------------------------------')
        DO 300 II=1,NION
           N0NEW=NFIRST(II)
           N1NEW=NLAST(II)
           IANEW=NUMAT(IATM(N0NEW))
           IZNEW=IZ(IEL(N0NEW))
           IF(N1OLD(IANEW,IZNEW).EQ.0) GO TO 300
           KOLD=N1OLD(IANEW,IZNEW)-N0OLD(IANEW,IZNEW)
           KNEW=NLAST(II)-NFIRST(II)
           IF(KOLD.LT.KNEW) KNEW=KOLD
           JL=N0OLD(IANEW,IZNEW)-1
           DO IL=NFIRST(II),NFIRST(II)+KNEW
              JL=JL+1
              WRITE(6,601) JL,IL
  601         FORMAT(10X,I8,5X,I8)
              DO ID=1,ND
                 POPUL0(IL,ID)=POPUL(JL,ID)
              END DO
           END DO
  300   CONTINUE
        DO 310 II=1,NATOM
           N0NEW=NKA(II)
           IANEW=NUMAT(IATM(N0NEW))
           IZNEW=IZ(IEL(N0NEW))+1
           IF(N0OLD(IANEW,IZNEW).EQ.0) GO TO 310
           WRITE(6,601) N0OLD(IANEW,IZNEW),N0NEW
           DO ID=1,ND
              POPUL0(N0NEW,ID)=POPUL(N0OLD(IANEW,IZNEW),ID)
           END DO
  310   CONTINUE
C
        DO II=1,NLEVEL
           DO ID=1,ND
              POPUL(II,ID)=POPUL0(II,ID)
           END DO
        END DO
      END IF
      LTE=LTE0
C
      END IF
      RETURN
      END
C
C
C     ****************************************************************
C
C

      SUBROUTINE ROSSTD(IJ)
C     =====================
C
C     Rosseland mean opacity
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'ATOMIC.FOR'
      INCLUDE 'MODELQ.FOR'
      INCLUDE 'ITERAT.FOR'
      INCLUDE 'ALIPAR.FOR'
      dimension pld(mdepth),abpld(mdepth)

C
C     for IJ > 0 - contribution from the frequency IJ to the Rosseland
C                  opacity integrals
C
      IF(IJ.GT.0) THEN
         if(ij.eq.1) then
         do id=1,nd
            pld(id)=0.
            abpld(id)=0.
         end do
         end if
         IF(IJX(IJ).GE.0) THEN
           DO ID=1,ND
             PLAN=XKFB(ID)/XKF1(ID)*W(IJ)
             DPLAN=PLAN/XKF1(ID)*FREQ(IJ)*HKT21(ID)
             ABROSD(ID)=ABROSD(ID)+DPLAN/ABSO1(ID)
             SUMDPL(ID)=SUMDPL(ID)+DPLAN
             pld(id)=pld(id)+plan
             abpld(id)=abpld(id)+(abso1(id)-scat1(id))*pld(id)
           END DO
         END IF
       ELSE   
C
C     for IJ=0 - evaluation of the Rosseland optical depth and
C                the radiative equilibrium division point
C
         ID=1
       IDR=0
         TAURS(ID)=HALF*DEDM1*ABROSD(ID)*DENS(ID)
         DO ID=2,ND         
            DTAUR=DELDM(ID-1)*(ABROSD(ID)+ABROSD(ID-1))
            TAURS(ID)=TAURS(ID-1)+DTAUR
          IF(TAURS(ID).LE.TAUDIV) IDR=ID
         END DO
C
C        in the last iteration, output of The Rosseland opacity and 
C        optical depth; skip the rest
C
         IF(LFIN) THEN
            DO ID=1,ND
               TROSS(ID)=TAURS(ID)
               abpl=abpld(id)/pld(id)/dens(id)
               pll=sig4p*4.*temp(id)**4
               WRITE(11,611) ID,DM(ID),TAURS(ID),ABROSD(ID),
     *                       TEMP(ID),ELEC(ID),DENS(ID),
     *                       pld(id),pll,abpl

            END DO
  611       FORMAT(I4,2X,1P6E12.4,2x,3e12.4)
            RETURN
         END IF
C
C        determination of the radiative equilibrium parameters
C        REINT and REDIF;
C        REINT=1 for all ID .le. ND-idlst
C        REDIF=1 for Rosseland optical depth > taudiv  (taur.gt.taudiv)
C                typically 0.1 - 0.5;
C              - this value has been empirically shown to yield
C                the best convergence rate
C
         IF(ITER.GT.ITNDRE) RETURN
         if(ndre.gt.1) then
           write(6,600)
           do id=1,nd
              if(id.lt.ndre) then
                 reint(id)=1.
                 redif(id)=0.
               else
                 reint(id)=0.
                 redif(id)=1.
               end if
            write(6,602) id,redif(id),reint(id)
            end do
          idr=ndre
          else if(ndre.eq.-1) then
           write(6,600)
            do id=1,nd
               reint(id)=1.
               redif(id)=0.
            end do
            redif(1)=1.
            do id=1,nd
            write(6,602) id,redif(id),reint(id)
            end do
          else
c
       DO 10 ID=1,ND
            REINT(ID)=UN
            REDIF(ID)=UN
            IF(ID.GT.ND-IDLST) REINT(ID)=0.
            IF(ID.LT.IDR) THEN
               REDIF(ID)=0.
               IF(MOD(NDRE,10).EQ.-1) THEN
                  REDIF(ID)=TAURS(ID)
                ELSE IF(MOD(NDRE,10).EQ.-2) THEN
                  REDIF(ID)=TAURS(ID)*TAURS(ID)
               END IF
            END IF
            IF(NDRE.LE.-10) THEN
               REDIF(ID)=REDIF(ID)/SIG4P/TEFF**4
               REINT(ID)=REINT(ID)/ABPLAD(ID)*DENS(ID)
            END IF
   10    CONTINUE
         ID=1
c         REDIF(ID)=0.
         IF(MOD(NDRE,10).EQ.-5) REDIF(ID)=UN
         IF(NDRE.LE.-10) THEN 
            REDIF(ID)=REDIF(ID)/SIG4P/TEFF**4
         END IF
         NDRE=1
         if(iter.eq.1) then
            write(6,600)
            do id=1,nd
               write(6,602) id,redif(id),reint(id)
            end do
          WRITE(6,601) IDR,IDR+1,ND-idlst
          WRITE(10,601) IDR,IDR+1,ND-idlst
       endif
       endif
          WRITE(6,601) IDR,IDR+1,ND-idlst
  600    FORMAT(/'  id        redif          reint'/)
  601    FORMAT(/' SCHEME OF RADIATIVE EQUIL. DETERMINED IN RESOLV'/
     *          ' ONLY INTEGRAL EQUATION FOR ID          <= ',I3/
     *          ' BOTH FOR                  ',I5,' <= ID <= ',I3/)
  602    format(i4,1p2e15.3)
      END IF
      RETURN
      END
C
C
C     ****************************************************************
C
C
      SUBROUTINE NEWPOP(ID,POP1)
C     ==========================
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'ATOMIC.FOR'
      INCLUDE 'MODELQ.FOR'
      INCLUDE 'ITERAT.FOR'
      DIMENSION POP1(MLEVEL),DPOP(MLEVEL),DPMAX(MDEPTH),
     *          SBW(MLEVEL)
c
      if(ioptab.lt.0) return
C
      DPMAX(ID)=0.
      DO 10 I=1,NLEVEL
         IF(POPUL(I,ID).GT.0.) 
     *      DPOP(I)=(POP1(I)-POPUL(I,ID))/POPUL(I,ID)
         IF(ABS(DPOP(I)).GT.DPMAX(ID)) THEN
            DPMAX(ID)=ABS(DPOP(I))
            IMAX=I
         END IF
         POPUL(I,ID)=POP1(I) 
   10 CONTINUE   
      WRITE(18,601) ITER,ILAM,ID,DPMAX(ID),IMAX
  601 FORMAT(3I5,1PE10.2,I6)
C
C     array of b-factors
C
      DO 110 I=1,NLEVEL
         BFAC(I,ID)=UN      
         SBW(I)=ELEC(ID)*SBF(I)*WOP(I,ID)
  110 CONTINUE      
      IF(.NOT.LTE.AND.IPSLTE.EQ.0) THEN        
      DO 120 ION=1,NION
         DO 120 I=NFIRST(ION),NLAST(ION)
             IF(POPUL(NNEXT(ION),ID).GT.0.)
     *       BFAC(I,ID)=POPUL(I,ID)/(POPUL(NNEXT(ION),ID)*SBW(I))
  120 CONTINUE 
      END IF     
      RETURN
      END
C
C
C     ****************************************************************
C
C

      SUBROUTINE SABOLF(ID)
C     =====================
C
C     Saha-Boltzmann factors (SBF)
C     and "upper sums" - sum of Saha-Boltzmann factors for upper, LTE,
C     levels which are not included explicitly (USUM), and derivatives
C     wrt. temperature (T) and electron density (DUSUMN)
C
C     Input: ID  - depth index
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'ATOMIC.FOR'
      INCLUDE 'MODELQ.FOR'
      PARAMETER (UH=1.5)
      PARAMETER (CMAX=2.154D4,CCON=2.0706D-16)
C
C     DCHI - approximate lowering of ionization potential for neutrals
C      Actual lowering is DCHI*effective charge, and is considered only
C      if IUPSUM(ION).GT.0
C
c
      if(ioptab.lt.0) return
c
      T=TEMP(ID)
      SQT=SQRT(T)
      ANE=ELEC(ID)
      STANE=SQRT(T/ANE)
      XMAX=CMAX*SQRT(STANE)
      TK=BOLK*T
      CON=CCON/T/SQT
C
C     Saha-Boltzmann factors
C
      DO 50 ION=1,NION
         QZ=IZ(ION)
         CFN=CON/G(NNEXT(ION))
         DCH=0.
         IUPS=IUPSUM(ION)
         SSBF=0.
         DSSBFT=0.
         USUM(ION)=0.
         DUSUMT(ION)=0.
         DUSUMN(ION)=0.
         nlst=nlast(ion)
         if(ifwop(nlst).ge.0) then
             nl1up=nquant(nlst)+1
          else
             nl1up=nquant(nlst)
         end if
         DO 10 II=NFIRST(ION),NLAST(ION)
            if(ifwop(ii).lt.0) then
               E=EH*QZ*QZ/TK
               SUM=0.
               DO 5 J=nl1up,NLMX
                  XJ=J
                  XI=J*J
                  X=E/XI
                  FI=XI*EXP(X)*WNHINT(J,ID)
                  SUM=SUM+FI
    5          CONTINUE
               g(ii)=sum*two
               gmer(imrg(ii),id)=g(ii)
            end if
          X=ENION(II)/TK
            if(x.gt.110.) x=110.
          SB=CFN*G(II)*EXP(X)
            SBF(II)=SB
            SSBF=SSBF+SB
            DSBF(II)=-(UH+ENION(II)/TK)/T
            DSSBFT=DSSBFT+SB*DSBF(II)
   10    CONTINUE
C
C     Upper sums
C
         if(ifwop(nlst).lt.0) go to 50
         IF(ION.EQ.IELHM) THEN
           USUM(ION)=0.
           DUSUMT(ION)=0.
           DUSUMN(ION)=0.
         GO TO 50
       END IF
C        if(ion.ge.0) go to 50
         if(iups.eq.0) then
C
C     1. More exact approach - using (exact) partition functions
C
         IAT=NUMAT(IATM(NFIRST(ION)))
         XMX=XMAX*SQRT(QZ)
         CALL PARTF(IAT,IZ(ION),T,ANE,XMX,U,DUT,DUN)
         EE=ENION(NFIRST(ION))/TK
       if(ee.gt.110.) ee=110.
         CFE=CFN*EXP(EE)
         USUM(ION)=CFE*U-SSBF
         DUSUMT(ION)=CFE*(DUT-U/T*(UH+EE))-DSSBFT
         DUSUMN(ION)=CFE*DUN
C        dusumt(ion)=0.
C        dusumn(ion)=0.
       xx=(ssbf-sbf(nfirst(ion)))/sbf(nfirst(ion))
c         IF(USUM(ION).LT.0.) THEN 
c         IF(USUM(ION).LT.0.or.ee.ge.109.OR.IAT.GT.2) THEN 
         IF(USUM(ION).LT.0.or.ee.ge.109.or.xx.lt.1.e-7) THEN 
            USUM(ION)=0.
            DUSUMT(ION)=0.
            DUSUMN(ION)=0.
         END IF
C
C     2. Approximate approach - summation over fixed number of upper
C        levels, assumed hydrogenic (ie. their ionization energy and
C        statistical weight hydrogenic)
C
         else if(iups.gt.0) then
         SUM=0.
         DSUM=0.
         E=EH*QZ*QZ/TK
         DO 30 J=NQUANT(NLAST(ION))+1,IUPS
            XI=J*J
            X=E/XI
            FI=XI*EXP(X)
            SUM=SUM+FI
            DSUM=DSUM-FI*(UH+X)/T
   30    CONTINUE
         USUM(ION)=SUM*CON*TWO
         DUSUMT(ION)=DSUM*CON*TWO
         DUSUMN(ION)=0.
C
c        3. occupation probability form
c
         else 
         SUM=0.
         DSUM=0.
         E=EH*QZ*QZ/TK
         DO 40 J=NQUANT(NLAST(ION))+1,NLMX
            XJ=J
            XI=J*J
            X=E/XI
            FI=XI*EXP(X)*WNHINT(J,ID)
            SUM=SUM+FI
            DSUM=DSUM-FI*(UH+X)/T
   40    CONTINUE
         USUM(ION)=SUM*CON*TWO
         DUSUMT(ION)=DSUM*CON*TWO
         DUSUMN(ION)=0.
         end if
   50 CONTINUE
      RETURN
      END
C
C
C     ****************************************************************
C
C
C

      SUBROUTINE OPACF1(IJ)
C     =====================
C
C     Absorption, emission, and scattering coefficients
C     at frequency IJ and for all depths
C
C     Input: IJ   opacity and emissivity is calculated for the
C                 frequency points with index IJ
C     Output: ABSO1 -  array of absorption coefficient
C             EMIS1 -  array of emission coefficient
C             SCAT1 -  array of scattering coefficient (all scattering
C                       mechanisms except electron scattering)
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'ATOMIC.FOR'
      INCLUDE 'MODELQ.FOR'
      INCLUDE 'ODFPAR.FOR'
      INCLUDE 'ALIPAR.FOR'
      common/hmolab/anh2(mdepth),anhm(mdepth)
      common/ipricr/iprcrs,nprcrs
      PARAMETER (C14=2.99793D14, c10=c14*1.d-4,CFF1=1.3727D-25)
      dimension anm(mdepth),opatab(mdepth),etatab(mdepth)
      dimension pold(mlevel),abtrh(mtrans)
c
      if(ioptab.lt.0) then
         call opact1(ij)
         return
      end if
C
C     initialize
c
c      IF(ICOMPT.GT.0) THEN
         DO ID=1,ND
            ELSCAT(ID)=ELEC(ID)*SIGEC(IJ)
         END DO
c      END IF
C
      DO ID=1,ND
         ABSO1(ID)=ELSCAT(ID)
         EMIS1(ID)=0.
         SCAT1(ID)=ELSCAT(ID)
      END DO
C
C     basic frequency- and depth-dependent quantities
C
      FR=FREQ(IJ)
      FRINV=UN/FR
      FR3INV=FRINV*FRINV*FRINV
      DO 10 ID=1,ND
         XKF(ID)=EXP(-HKT1(ID)*FR)
         XKF1(ID)=UN-XKF(ID)
         XKFB(ID)=XKF(ID)*BNUE(IJ)
   10 CONTINUE
c 
      if(iprcrs.gt.0) then
         abso1(iprcrs)=0.
         do ii=nfirst(ielh),nlast(ielh)
            if(ii.ne.nprcrs+nfirst(ielh)-1) then
               pold(ii)=popul(ii,iprcrs)
               popul(ii,iprcrs)=0.
               do jj=ii+1,nnext(ielh)
                  itrh=itra(ii,jj)
                  abtrh(itrh)=abtra(itrh,iprcrs)
                  abtra(itrh,iprcrs)=0.
               end do
            end if
         end do
      end if
C
C ********  1a. bound-free contribution - without dielectronic rec.
C
      if(ifdiel.eq.0) then
      DO IBFT=1,NTRANC
         ITR=ITRBF(IBFT)
         SG=CROSS(IBFT,IJ)
         IF(SG.GT.0.) THEN
         II=ILOW(ITR)
         JJ=IUP(ITR)
         IZZ=IZ(IEL(II))
         IMER=IMRG(II)
         DO ID=1,ND
            SGD=SG
            IF(MCDW(ITR).GT.0) THEN
               CALL DWNFR1(FR,FR0(ITR),ID,IZZ,DW1)
               DWF1(MCDW(ITR),ID)=DW1
               SGD=SG*DW1
            END IF
            IF(IFWOP(II).LT.0) THEN
               CALL SGMER1(FRINV,FR3INV,IMER,ID,SGME1)
               SGMG(IMER,ID)=SGME1
               SGD=SGME1
            END IF
            EMISBF=SGD*EMTRA(ITR,ID)
            ABSO1(ID)=ABSO1(ID)+SGD*ABTRA(ITR,ID)
            EMIS1(ID)=EMIS1(ID)+EMISBF
       END DO
       END IF
      END DO
      else
C
C ********  1b. bound-free contribution - with dielectronic rec.
C
      DO IBFT=1,NTRANC
         ITR=ITRBF(IBFT)
         SG=CROSS(IBFT,IJ)
         IF(SG.GT.0.) THEN
         II=ILOW(ITR)
         JJ=IUP(ITR)
         IZZ=IZ(IEL(II))
         IMER=IMRG(II)
         DO ID=1,ND
            SG=CROSSD(IBFT,IJ,ID)
            IF(SG.GT.0.) THEN
            SGD=SG
            IF(MCDW(ITR).GT.0) THEN
               CALL DWNFR1(FR,FR0(ITR),ID,IZZ,DW1)
               DWF1(MCDW(ITR),ID)=DW1
               SGD=SG*DW1
            END IF
            IF(IFWOP(II).LT.0) THEN
               CALL SGMER1(FRINV,FR3INV,IMER,ID,SGME1)
               SGMG(IMER,ID)=SGME1
               SGD=SGME1
            END IF
            EMISBF=SGD*EMTRA(ITR,ID)
            ABSO1(ID)=ABSO1(ID)+SGD*ABTRA(ITR,ID)
            EMIS1(ID)=EMIS1(ID)+EMISBF
          END IF
       END DO
       END IF
      END DO
      end if
C
C ******** 2. free-free contribution
C
      DO 40 ION=1,NION
         IT=ITRA(NNEXT(ION),NNEXT(ION))
C
C        hydrogenic with Gaunt factor = 1
C
         IF(IT.EQ.1) THEN
            DO ID=1,ND
               SF1=SFF3(ION,ID)*FR3INV
               SF2=SFF2(ION,ID)
               IF(FR.LT.FF(ION)) SF2=UN/XKF(ID)
               ABSOFF=SF1*SF2
               ABSO1(ID)=ABSO1(ID)+ABSOFF
               EMIS1(ID)=EMIS1(ID)+ABSOFF
            END DO
C
C         hydrogenic with exact Gaunt factor 
C
          ELSE IF(IT.EQ.2) THEN
            DO ID=1,ND
               SF1=SFF3(ION,ID)*FR3INV
               SF2=SFF2(ION,ID)
               IF(FR.LT.FF(ION)) SF2=UN/XKF(ID)
               X=C14*CHARG2(ION)/FR
               SF2=SF2-UN+GFREE1(ID,X)
               ABSOFF=SF1*SF2
               ABSO1(ID)=ABSO1(ID)+ABSOFF
               EMIS1(ID)=EMIS1(ID)+ABSOFF
            END DO
C
C         H minus free-free opacity
C
          ELSE IF(IT.EQ.3) THEN
            DO ID=1,ND
               ABSOFF=(CFF1+CFFT(ID)*FRINV)*CFFN(ID)*FRINV
               ABSO1(ID)=ABSO1(ID)+ABSOFF
               EMIS1(ID)=EMIS1(ID)+ABSOFF
            END DO
C
C         special evaluation of the cross-section
C
          ELSE IF(IT.LT.0) THEN
            DO ID=1,ND
               ABSOFF=FFCROS(ION,IT,TEMP(ID),FR)*
     *                POPUL(NNEXT(ION),ID)*ELEC(ID)
               ABSO1(ID)=ABSO1(ID)+ABSOFF
               EMIS1(ID)=EMIS1(ID)+ABSOFF
            END DO
         END IF
   40 CONTINUE
C
C     ********  3. - additional continuum opacity (OPADD)
C
      IF(IOPADD.NE.0) THEN
         ICALL=1
         DO ID=1,ND
            CALL OPADD(0,ICALL,IJ,ID)
            ABSO1(ID)=ABSO1(ID)+ABAD
            EMIS1(ID)=EMIS1(ID)+EMAD
            SCAT1(ID)=SCAT1(ID)+SCAD
         END DO
      END IF
C
C ********  4. - opacity and emissivity in lines
C
      IF(ISPODF.EQ.0) THEN
      IF(IJLIN(IJ).GT.0) THEN
C
C     the "primary" line at the given frequency
C
         ITR=IJLIN(IJ)
         DO 50 ID=1,ND
            SG=PRFLIN(ID,IJ)
            ABSO1(ID)=ABSO1(ID)+SG*ABTRA(ITR,ID)
            EMIS1(ID)=EMIS1(ID)+SG*EMTRA(ITR,ID)
   50    CONTINUE
      ENDIF
      IF(NLINES(IJ).LE.0) GO TO 200
C
C     the "overlapping" lines at the given frequency
C
      DO 100 ILINT=1,NLINES(IJ)
         ITR=ITRLIN(ILINT,IJ)
       if(linexp(itr)) goto 100
         IJ0=IFR0(ITR)
         DO 60 IJT=IJ0,IFR1(ITR)
            IF(FREQ(IJT).LE.FR) THEN
               IJ0=IJT
               GO TO 70
            END IF
   60    CONTINUE
   70    IJ1=IJ0-1
         A1=(FR-FREQ(IJ0))/(FREQ(IJ1)-FREQ(IJ0))
         A2=UN-A1
         DO 80 ID=1,ND
            SG=A1*PRFLIN(ID,IJ1)+A2*PRFLIN(ID,IJ0)
            ABSO1(ID)=ABSO1(ID)+SG*ABTRA(ITR,ID)
            EMIS1(ID)=EMIS1(ID)+SG*EMTRA(ITR,ID)
   80    CONTINUE
  100 CONTINUE
  200 CONTINUE
C
C     Opacity sampling option
C
      ELSE
      IF(NLINES(IJ).LE.0) GO TO 400
      DO 300 ILINT=1,NLINES(IJ)
        ITR=ITRLIN(ILINT,IJ)
      KJ=IJ-IFR0(ITR)+KFR0(ITR)
      INDXPA=IABS(INDEXP(ITR))
      IF(INDXPA.NE.3 .AND. INDXPA.NE.4) THEN
        DO ID=1,ND
            SG=PRFLIN(ID,KJ)
            ABSO1(ID)=ABSO1(ID)+SG*ABTRA(ITR,ID)
            EMIS1(ID)=EMIS1(ID)+SG*EMTRA(ITR,ID)
        END DO
      ELSE
        DO ID=1,ND
            KJD=JIDI(ID)
c     SG=EXP(XJID(ID)*SIGFE(KJD,KJ)+(UN-XJID(ID))*SIGFE(KJD+1,KJ))
            ABSO1(ID)=ABSO1(ID)+SG*ABTRA(ITR,ID)
            EMIS1(ID)=EMIS1(ID)+SG*EMTRA(ITR,ID)
        END DO
      ENDIF
  300 CONTINUE
  400 CONTINUE
      ENDIF
C
c     Lyman alpha and beta quasimolecular opacity 
c
      call quasim(ij)
c
C     ----------------------------
C     total opacity and emissivity
C     ----------------------------
C
      DO ID=1,ND
         ABSO1(ID)=ABSO1(ID)-EMIS1(ID)*XKF(ID)
         EMIS1(ID)=EMIS1(ID)*XKFB(ID)
         absot(id)=abso1(id)
      END DO
c
c     ---------------------------------
c     hydrogen pacity from Gomez tables
c     ---------------------------------
c
c     call ghydop(ij)
c
c     --------------------------------------------------------
c     contribution from precalculated background opacity table
c     --------------------------------------------------------`
c
      if(ioptab.gt.0) then
         if(ioptab.eq.1) then
            do id=1,nd
               anm(id)=elec(id)
            end do
          else
            do id=1,nd
               anm(id)=dens(id)/wmm(id)
            end do
         end if
         xlam=2.997925e18/fr
         call opa_ot(nd,xlam,temp,anm,opatab,etatab)
         abso1(id)=abso1(id)+opatab(id)
         emis1(id)=emis1(id)+etatab(id)
         absot(id)=abso1(id)
      end if
c
c     H2-H2 CIA opacity
c
      if(ifcia.gt.0) then
      do id=1,nd
c        call cia_sub(temp(id),anh2(id),fr,oph2)
         abso1(id)=abso1(id)+oph2
         emis1(id)=emis1(id)+oph2*xkfb(id)/xkf1(id)
      end do
      end if
C
C     if needed, evaluate the opacity per gram
C
      if(izscal.eq.0) then
         do id=1,nd
            absot(id)=abso1(id)*dens1(id)
         end do
      end if
c
      if(ifprd.gt.0) call prd(ij)
c
      if(iprcrs.gt.0) then
         ih=nfirst(ielh)+nprcrs-1
         crs=abso1(iprcrs)/(popul(ih,iprcrs)*g(ih)*
     *       0.0265*4.1347e-15)
         write(37,637),ij,fr,abso1(iprcrs),crs
         do ii=nfirst(ielh),nlast(ielh)
            if(ii.ne.ih) then
               popul(ii,iprcrs)=pold(ii)
               do jj=ii+1,nnext(ielh)
                  itrh=itra(ii,jj)
                  abtra(itrh,iprcrs)=abtrh(itrh)
               end do
            end if
         end do
      end if
  637 format(i5,1p3e15.6)

      RETURN
      END
C
C
C     ****************************************************************
C
C
      SUBROUTINE DOPGAM(ITR,ID,T,DOP,AGAM)
C     ====================================
C
C     Doppler width and the Voigt damping parameter for the line ITR
C
C     Input:
C      ITR   - index of transition
C      ID    - depth index
C      T     - temperature
C     Output:
C      DOP   - Doppler width
C      AGAM  - total damping parameter (in units of Doppler widths;
C              ie. = gam/4pi/DOP, where gam is the "physical" damping
C              parameter expresed in circular frequencies)
C
C     Damping parameter is calculated only for transitions with
C     |IPROF(ITR)| = 1 (ie. those for which either Voigt or some non-
C     standard profile is assumed)
C     Determination of AGAM:
C     is controlled by input parameters transmitted by COMMON/VOIPAR:
C
C      GAMAR(IP)        - for > 0  - has the meaning of natural damping
C                                    parameter (=Einstein coefficient for
C                                    spontaneous emission)
C                             = 0  - classical natural damping assumed
C                             < 0  - damping is given by a non-standard,
C                                    user supplied procedure GAMSP
C      STARK1(IP)       -     = 0  - Stark broadening neglected
C                             < 0  - scaled classical expression
C                                    (ie gam = -STARK1 * classical Stark)
C                             > 0  - Stark broadening given by
C                                    n(el)*[STARK1*T**STARK2 + STARK3]
C      STRAK2, STARK3   -   see above
C      VDWH(IP)         -   .le.0  - Van der Waals broadening neglected
C                             > 0  - scaled classical expression
C
C      the corresponding index IP is given by ITRA(IUP(ITR),ILOW(ITR))
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'ATOMIC.FOR'
      INCLUDE 'MODELQ.FOR'
      PARAMETER (BOL2=2.76108D-16, CIN=UN/2.997925D10)
C
      J=IUP(ITR)
      IAT=IATM(J)
      FR=FR0(ITR)
      IE=IEL(J)
      AM=BOL2/AMASS(IAT)*T
      AGAM=0.
C
C     Doppler width
C
      DOP=FR*CIN*SQRT(AM+VTURBS(ID)*VTURBS(ID))
C
C     -----------------
C     damping parameter - only for IPROF = 1
C
      IF(IABS(IPROF(ITR)).NE.1) RETURN
      IP=ITRA(J,ILOW(ITR))
      ANE=ELEC(ID)
      If(GAMAR(IP).GT.0.) THEN
C
C        Natural (radiation) broadening
C
         AGAM=GAMAR(IP)
       ELSE IF (GAMAR(IP).EQ.0.) THEN
         AGAM=2.47342D-22*FR*FR
       ELSE
C
C        Non-standard expression - for the total damping parameter,
C        not only for radiation damping
C
         CALL GAMSP(ITR,T,ANE,AGAM)
      END IF
C
C     Stark broadening
C
      IF(STARK1(IP).LT.0.) THEN
         ANFF=IZ(IE)*IZ(IE)*EH/ENION(J)
         AGAM=AGAM-1.D-8*ANFF**2.5*ANE*STARK1(IP)
       ELSE IF (STARK1(IP).GT.0.) THEN
         AGAM=AGAM+ANE*(STARK1(IP)*T**STARK2(IP)+STARK3(IP))
      END IF
C
C     Van der Waals broadening
C
      IF(VDWH(IP).GT.0.) THEN
         AH=DENS(ID)/WMM(ID)/YTOT(ID)
         AHE=AH*0.1
         IF(IELHE1.NE.0) AHE=AHE*10.*ABUND(IATHE,ID)
         AGAM=AGAM+T**0.3*(AH+0.42*AHE)*VDWH(IP)
      END IF
C
C     Total damping parameter in units of Doppler widths
C
      AGAM=AGAM/DOP/12.566370614
      RETURN
      END
C
C
C     ****************************************************************
C
C
      SUBROUTINE GAMSP(ITR,T,ANE,AGAM)
C     ================================
C
C     Non-standard expression for the damping parameter  -
C     a user-supplied procedure
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      AGAM=0.
      if(itr.le.0) return
      t1=t
      ane1=ane
      RETURN
      END
C
C
C     ****************************************************************
C
C
      FUNCTION PROFIL(FR,A,DOP,ITR,IP,ID)
C     ===================================
C
C     Standard absorption profile - normalized to unity
C
C     Input:
C     FR  - frequency
C     A   - Voigt damping parameter
C     DOP - Doppler width
C     ITR - transition index
C     ID  - depth index
C
C     Profile is evaluated differently for different IP=IPROF(ITR):
C     IP = 0   -  Doppler profile
C     IP = 1   -  Voigt profile
C     IP = 2   -  approximate Stark (+ Doppler) profile for hydrogen lines;
C                 however, the routine is called with IP=2 only from
C                 START, i.e. for the initialization
C     IP > 9   -  non-standard profile, given by a user-supplied
C                 procedure PROFSP
C     IP = 12  -  approximate Stark profile for hydrogen lines
C                  (Klaus Wener's routines)
C
C     V - frequency displacement from the line center in units of
C         Doppler width
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'ATOMIC.FOR'
      INCLUDE 'MODELQ.FOR'
      common/quasun/tqmprf,iquasi,nunalp,nunbet,nungam,nunbal
      PARAMETER (PISQ=1.77245385090551D0,PISQ1=UN/PISQ)
C
      PROFIL=0.
      V=(FR-FR0(ITR))/DOP
      IPA=IABS(IP)
      IF(IPA.EQ.0) THEN
         IF(V.LE.13.) PROFIL=EXP(-V*V)*PISQ1
       ELSE IF(IPA.EQ.1) THEN
         PROFIL=VOIGT(V,A)*PISQ1
       ELSE IF(IPA.EQ.2) THEN
         IF(ID.GT.0) THEN
            ANE=ELEC(ID)
          ELSE
            ANE=1.e9*grav
         END IF
         if(ane.le.0.) ane=1.e14
         F000=EXP(0.666666667*LOG(ANE))
         II=NQUANT(ILOW(ITR))
         JJ=NQUANT(IUP(ITR))
         IZZ=IZ(IEL(ILOW(ITR)))
         FAC=TWO
         if(iquasi.gt.0.and.ii.eq.1) then
            if(jj.eq.2) fac=un
            if(jj.eq.3.and.iquasi.gt.1) fac=un
         end if
         F00=1.25D-9*F000
         IF(IZZ.EQ.2) THEN
            FAC=UN
            F00=3.906D-11*F000
         END IF
         CALL STARK0(II,JJ,IZZ,XKIJ,WL0,FIJ)
         FXK=F00*XKIJ
         DBETA=WL0*WL0/2.997925D18/FXK
         BETAD=DOP*DBETA
         CALL DIVSTR(IZZ)
         BETA=DBETA*ABS(FR-FR0(ITR))
         SG=STARKA(BETA,fac)*BETAD
         PROFIL=SG
       ELSE IF(IPA.GT.10) THEN
         PROFIL=PROFSP(FR,DOP,ITR,ID)
      END IF
      RETURN
      END
C
C
C     ****************************************************************
C
C

      FUNCTION VOIGT(V,AGAM)
C     ======================
C
C     Voigt function
C     Procedure after Matta and Reichel, 1971, Math.Comp. 25, 339.
C
      INCLUDE 'IMPLIC.FOR'
      DIMENSION HN(12),EN(12)
      PARAMETER (PI=3.141592653589793D0, M=12, HH=0.5D0, UN=1.D0)
      PARAMETER (PISQ=1.77245385090551D0,PISQ1=UN/PISQ)
      DATA ICOMP /0/
      SAVE EN,HN,PH,HP,ICOMP
C
C     Initialization of auxiliary quantities
C
      IF(ICOMP.EQ.0) THEN
         HP=HH*PISQ1
         PH=PI/HH
         DO I=1,M
            XI=I
            U=XI*XI*HH*HH
            EN(I)=EXP(-U)
            HN(I)=4.D0*U
         END DO
         ICOMP=1
      END IF
C
C     Main term
C
      AGAM1=UN/AGAM
      X=V*AGAM1
      T=0.25D0*AGAM1*AGAM1
      X2=X*X
      X4=4.D0*X2
      S1=UN+X2
      S2=UN-X2
      U0=0.
      DO I=1,M
         S0=HN(I)*T
         U=EN(I)/((S2+S0)*(S2+S0)+X4)
         U0=U0+U*(S1+S0)
      END DO
      S2=UN/S1
      U0=HP*(S2+2.D0*U0)
C
C     Correction term
C
      IF(T.LT.0.25D0/PH/PH) GO TO 50
      U=X/2.D0/T
      A=COS(U)
      B=SIN(U)
      TSQ1=UN/SQRT(T)
      S1=PH*TSQ1
      S2=S1*X
      C=EXP(-S1)-COS(S2)
      D=SIN(S2)
      T4=0.25D0/T
      U=EXP(-X2*T4-S1+T4)*PISQ*TSQ1/(C*C+D*D)
      U0=U0+U*(A*C-B*D)
   50 VOIGT=U0*AGAM1*PISQ1
      RETURN
      END
C
C
C     ****************************************************************
C
C
      function voigte(vs,a)
c     =====================
c
c  computes a voigt function  h = h(a,v)
c  a=gamma/(4*pi*dnud)   and  v=(nu-nu0)/dnud.  this  is  done after
c  traving (landolt-b\rnstein, p. 449).
c
      INCLUDE 'IMPLIC.FOR'
      real*4 vs,a
      PARAMETER (UN=1., TWO=2.)
      dimension ak(19),a1(5)
      data ak      /-1.12470432, -0.15516677,  3.28867591, -2.34357915,
     ,  0.42139162, -4.48480194,  9.39456063, -6.61487486,  1.98919585,
     , -0.22041650, 0.554153432, 0.278711796,-0.188325687, 0.042991293,
     ,-0.003278278, 0.979895023,-0.962846325, 0.532770573,-0.122727278/
      data sqp/1.772453851/,sq2/1.414213562/
c
      v = abs(vs)
      u = a + v
      v2 = v*v
      if (a.eq.0.0) go to 140
      if (a.gt.0.2) go to 120
      if (v.ge.5.0) go to 121
c
      ex=0.
      if(v2.lt.100.) ex = exp(-v2)
      k = 1
c
  100 quo = un
      if (v.lt.2.4) go to 101
      quo = un/(v2 - 1.5)
      m = 11
      go to 102
c
  101 m = 6
      if (v.lt.1.3) m = 1
  102 do 103 i=1,5
         a1(i) = ak(m)
         m = m + 1
  103 continue
      h1 = quo*(a1(1) + v*(a1(2) + v*(a1(3) + v*(a1(4) + v*a1(5)))))
      if (k.gt.1) go to 110
c
c a le 0.2  and v lt 5.
c
      h = h1*a + ex*(un + a*a*(un - two*v2))
      voigte=h
      return
c
  110 pqs = two/sqp
      h1p = h1 + pqs*ex
      h2p = pqs*h1p - two*v2*ex
      h3p = (pqs*(un - ex*(un - two*v2)) - two*v2*h1p)/3. + pqs*h2p
      h4p = (two*v2*v2*ex - pqs*h1p)/3. + pqs*h3p
      psi = ak(16) + a*(ak(17) + a*(ak(18) + a*ak(19)))
c
c 0.2 lt a le 1.4  and  a + v le 3.2
c
      h = psi*(ex + a*(h1p + a*(h2p + a*(h3p + a*h4p))))
      voigte=h
      return
c
  120 if (a.gt.1.4.or.u.gt.3.2) go to 130
      ex=0.
      if(v2.lt.100.)ex = exp(-v2)
      k = 2
      go to 100
c
c a le 0.2  and  v ge 5.
c
  121 h = a*(15. + 6.*v2 + 4.*v2*v2)/(4.*v2*v2*v2*sqp)
      voigte=h
      return
c
  130 a2 = a*a
      u = sq2*(a2 + v2)
      u2 = un/(u*u)
c
c a gt 1.4  or  a + v gt 3.2
c
      h = sq2/sqp*a/u*(1. + u2*(3.*v2 - a2) +
     ,        u2*u2*(15.*v2*v2 - 30.*v2*a2 + 3.*a2*a2))
      voigte=h
      return
c
c a eq 0.
c
  140 h=0.
      if(v2.lt.100.) h=exp(-v2)
      voigte=h
      return
      end
C
C
C     ****************************************************************
C
C
      FUNCTION PROFSP(FR,DOP,ITR,ID)
C     ================================
C
C     Non-standard absorption profile - normalized to unity;
C     a user supplied procedure
C
C     Input:
C     FR  - frequency
C     A   - Voigt damping parameter
C     DOP - Doppler width
C     ITR - transition index
C     ID  - depth index
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'ATOMIC.FOR'
      INCLUDE 'MODELQ.FOR'
C
      PROFSP=0.
      IP=IPROF(ITR)
C
C     Klaus Werner's Voigt+Stark wing profile (formula A.3.4)
C
      IF(ABS(IP).NE.12) RETURN
      
C     1- Stark wings

      II=NQUANT(ILOW(ITR))
      JJ=NQUANT(IUP(ITR))
      SIJ=JJ*(JJ-1)+II*(II-1)
C
C     Micro-field
C     
      ZMIKRO=0.
      DO IAT=1,NATOM
         N0I=N0A(IAT)
         NKI=NKA(IAT)
         DO I=N0I,NKI-1
            IE=IEL(I)
            CH=IZ(IE)-1
            CH32=CH*SQRT(CH)
            ZMIKRO=ZMIKRO+CH32*POPUL(I,ID)
         END DO
         CH=CH+UN
         CH32=CH*SQRT(CH)
         ZMIKRO=ZMIKRO+CH32*POPUL(NKI,ID)
      END DO
      CALL SABOLF(ID)
      DO ION=1,NION
         CH=IZ(ION)-1
         CH32=CH*SQRT(CH)
         ZMIKRO=ZMIKRO+CH32*USUM(ION)
      END DO
      ZMIKRO=ZMIKRO**0.6666667
C     
      IAT=IATM(ILOW(ITR))
      IE=IEL(ILOW(ITR))
      CH=IZ(IE)
      DBETA=1.385*CH/SIJ/ZMIKRO
C
C     empirical correction in PRO2
C
      CORRE=UN
      IF(IE.EQ.IELHE2 .AND.ILOW(ITR).GT.(NFIRST(IE)+1))
     *   CORRE=HALF
      IF(IAT.NE.IATH .AND. IAT.NE.IATHE) THEN
         CORRE=UN/(CH-UN)
      END IF
      DBETA=DBETA/CORRE
C     
      BETAD=DOP*DBETA
      BETA=DBETA*ABS(FR-FR0(ITR))
      SIGST=UBETA(BETA)*BETAD
      
C     2- Voigt profile

      AGAMS=5.E-5*ELEC(ID)/SQRT(TEMP(ID))*JJ*JJ/CH/CH
      AGAM=2.47342D-22*FR*FR+AGAMS
      AA=AGAM/12.56637/DOP
      V=(FR-FR0(ITR))/DOP
      SIGVT=VOIGT(V,AA)/1.77245385090551D0
      SGA=SIGVT
      IF(SIGST.GT.SIGVT) SGA=SIGST
      PROFSP=SGA
      RETURN
      END
C
C
C     ****************************************************************
C
C
      FUNCTION UBETA(BETA)                                             
C     =====================                                                       
     
            
C
C **********************************************************************
C ***                                                                   
C *** CALLED BY: SUBROUTINE PROFSP                                        
C *** CALLS:     SUBROUTINE LAGRAN                                      
C ***                                                                   
C *** Interpolation of function U(beta) in table from Dien (ApJ 109,452)
C ***                                                                   
C **********************************************************************
      INCLUDE 'IMPLIC.FOR'
C                                                                        
      DIMENSION B0(46),U0(46)                                           
C
      DATA (B0(I),U0(I),I= 1,46) /                                      
     1 0.0,.287,0.1,.286,0.2,.283,0.3,.278,0.4,.271,0.5,.262,0.6,.252,  
     1 0.7,.240,0.8,.228,0.9,.215,1.0,.202,1.1,.188,1.2,.174,1.3,.161,  
     1 1.4,.148,1.5,.135,1.6,.124,1.7,.113,1.8,.1024,1.9,.0928,         
     1 2.0,.0839,2.1,.0758,2.2,.0684,2.3,.0617,2.4,.0557,2.5,.0502,     
     1 2.6,.0454,2.7,.0411,2.8,.0373,2.9,.0338,3.0,.0310,3.2,.0260,     
     1 3.4,.0220,3.6,.0187,3.8,.0160,4.0,.0238,4.2,.0120,4.4,.0104,     
     1 4.6,.0091,4.8,.0080,5.0,.0071,6.0,.0041,7.0,.0027,8.0,.0018,     
     1 9.0,.0014,10.0,.0011/                                            
C                                                                        
C *** asymptotic value                                                  
C
      IF(BETA.GT.10.) THEN                                             
         UBETA=0.2992*BETA**(-2.5)  
       ELSE
C *** Interpolation                                                              
C
         DO I=3,46                                                       
            IF (BETA.LT.B0(I)) GO TO 2
         END DO                                        
   2     CALL LAGRAN (B0(I-2),B0(I-1),B0(I),U0(I-2),U0(I-1),
     *                U0(I), BETA,UBETA)        
      END IF                                          
      RETURN                                                            
      END                                                               
C
C
C     ****************************************************************
C
C
      SUBROUTINE LAGRAN(X0,X1,X2,Y0,Y1,Y2,X,Y)                         
C     =========================================                               
C
C ***                                                                   
C *** Lagrange interpolation for three points                           
C ***
      INCLUDE 'IMPLIC.FOR'
C                                                                        
      XL0=(X-X1)*(X-X2)/(X0-X1)/(X0-X2)                                 
      XL1=(X-X0)*(X-X2)/(X1-X0)/(X1-X2)                                 
      XL2=(X-X0)*(X-X1)/(X2-X0)/(X2-X1)                                 
                                                                        
      Y=Y0*XL0+Y1*XL1+Y2*XL2                                            
                                                                        
      RETURN                                                            
      END                                                               

C
C
C     ****************************************************************
C
C
      SUBROUTINE LINSET(ITR,IUNIT,IFRQ0,IFRQ1,XMAX,DOP,AGAM)
C     ======================================================
C
C     Set up frequency points and weights for a line
C     Auxiliary procedure for START
C
C     Input:
C      ITR    -  index of the transition
C      INMOD  -  mode of evaluating frequency points in the line
C             =  0  - means that frequency points and weights have
C                     already been read;
C             ne 0  - frequency points and weights will be evaluated:
C             the meaning of the individual values:
C             =  1  - equidistant frequencies, trapezoidal integration
C             =  2  - equidistant frequencies, Simpson integration
C             =  3  - modified Simpson integration, i.e. a series of
C                     3-point Simpson integrations with each subsequent
C                     integration interval doubled, until the whole
C                     integartion area is covered
C             =  4  - frequencies (in units of standard x) and weights
C                     (for integration over x) are read;
C
C      XMAX   >  0  - means that the line is assumed symmetric around the
C                     center, frequency points are set up between x=0 and
C                     x=XMAX, where x is frequency difference from the
C                     line center in units of the standard Doppler width
C             <  0  - frequency points are set between x=XMAX and x=-XMAX
C      DOP    -  Doppler width
C      AGAM   -  damping parameter (for lines with Voigt or non-standard
C                profile only)
C
C     Output (to COMMON/FRQEXP)
C      FREQ   -  array of frequencies
C                Note that LINSET calculates values of frequencies only
C                for frequency points between IFR0(ITR) and IFR1(ITR).
C      W      -  corresponding integration weights
C      PROF   -  corresponding values of absorption profile
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'ATOMIC.FOR'
      INCLUDE 'MODELQ.FOR'
      PARAMETER (BOL2=2.76108D-16, CIN=UN/2.997925D10, OS0=0.02654,
     *           F0C1=1.25D-9, TTW=2./3., PISQ1=UN/1.77245385090551D0,
     *           C18IN=UN/CAS,F0C2=3.906D-11)
      DIMENSION X(MFREQL),W0(MFREQL)
C
      IF(ITR.EQ.0) GO TO 200
      INMOD=INTMOD(ITR)
      INMOD0=MOD(IABS(INMOD),10)
      IJ0=IFR0(ITR)
      IJ1=IFR1(ITR)
      N=IJ1-IJ0+1

      IF(INMOD.EQ.0) THEN
         S=OSC0(ITR)*0.02654D0
         IP0=IPROF(ITR)
         IP=IABS(IP0)
         DO I=1,N
            PROF(I+IJ0-1)=PROFIL(FREQ(I+IJ0-1),AGAM,DOP,ITR,IP,0)*
     *                    S/DOP
            IJLIN(I+IJ0-1)=ITR
         END DO
         IF(IP0.GE.0) THEN
            IF(XMAX.LT.0.) THEN
               PROF(IJ0)=0.
               PROF(IJ1)=0.
             ELSE
               PROF(IJ1)=0.
            END IF
         END IF
         IF(IPROF(ITR).GE.0) THEN
            IF(XMAX.LT.0.) THEN
               INTMOD(ITR)=-2
             ELSE
               INTMOD(ITR)=-1
            END IF
          ELSE
            IF(XMAX.LT.0.) THEN
               INTMOD(ITR)=2
             ELSE
               INTMOD(ITR)=1
            END IF
         END IF
         RETURN
      END IF

      X0=0.
      IF(XMAX.LT.0) X0=XMAX
      M=(N-1)/2
      X(1)=0.
      W0(1)=UN
      IF(N.LE.1) GO TO 100
C
      IF(INMOD0.LE.2) THEN
C
C        Trapezoidal integration
C
         HH=ABS(X0+XMAX)/(N-1)
         DO I=1,N
            X(I)=X0+(I-1)*HH
         END DO
         IF(INMOD0.EQ.2) GO TO 40
         DO I=1,N
            W0(I)=HH
         END DO
         W0(1)=HALF *HH
         W0(N)=HALF *HH
         GO TO 100
C
C        Ordinary Simpson integration
C
   40    HH=HH/3.D0
         IF(MOD(N,2).NE.1) 
     *   CALL QUIT('even number of points in Simpson - LINSET',n,n)
         DO I=1,M
            W0(2*I)=4.D0*HH
            W0(2*I+1)=2.D0*HH
         END DO
         W0(1)=HH
         W0(N)=HH
       ELSE IF(INMOD0.EQ.3) THEN
C
C      Modified Simpson integration - a set of 3-point Simpson
C      integrations with continuosly increasing distance between the
C      integration points (schematically, distances between points are
C      h,h,2h,2h,4h,4h, etc.)
C
         TWI=UN
         MM=M
         IF(MOD(N,2).NE.1) 
     *   CALL QUIT('even number of points in MSimpson - LINSET',n,n)
         IF(XMAX.LT.0) MM=M/2
         DO I=1,MM
            TWI=TWI*2.D0
            X(2*I+1)= TWI-UN
            X(2*I)=TWI-UN-TWI/4.D0
            W0(2*I)=2.D0*TWI
            W0(2*I+1)=1.5D0*TWI
         END DO
         TWN=TWI-UN
         TWA=ABS(XMAX)/TWN
         HH=TWA/6.D0
         DO I=1,MM
            X(2*I+1)=X(2*I+1)*TWA
            X(2*I)=X(2*I)*TWA
            W0(2*I)=W0(2*I)*HH
            W0(2*I+1)=W0(2*I+1)*HH
         END DO
         W0(1)=HH
         W0(N)=TWI*HH/2.D0
         X(1)=0.
         IF(M.EQ.MM) GO TO 100
         IF(MOD(N,4).NE.1) 
     *   CALL QUIT('conflict in MSimpson - LINSET',n,n)
         DO I=1,M
            X(M+1+I)=X(I+1)
            W0(M+1+I)=W0(I+1)
         END DO
         M2=2*(M+1)
         DO I=1,M
            X(I)=-X(M2-I)
            W0(I)=W0(M2-I)
         END DO
         X(M+1)=0.
         W0(M+1)=2.D0*HH
C
C        frequencies (in units of standard x) and weights
C        (for integration over x) are read;
C
       ELSE IF(INMOD0.EQ.4) THEN
         READ(IUNIT,*) (X(I),I=1,N),(W0(I),I=1,N)
      END IF
C
C     For all types of integration:
C     set up arrays  FR,WW,PRF
C
  100 S=OSC0(ITR)*0.02654D0
      IP0=IPROF(ITR)
      IP=IABS(IP0)
      DO I=1,N
         FREQ(I+IJ0-1)=FR0(ITR)-DOP*X(I)
         W(I+IJ0-1)=DOP*W0(I)
         PROF(I+IJ0-1)=PROFIL(FREQ(I+IJ0-1),AGAM,DOP,ITR,IP,0)*S/DOP
      END DO
C
C     for IPROF(ITR) ge 0  -  endpoint(s) of the line profile are forced to
C                             have zero cross-section
C
      IF(IP0.GE.0) THEN
         IF(XMAX.LT.0.) THEN
            PROF(IJ0)=0.
            PROF(IJ1)=0.
          ELSE
            PROF(IJ1)=0.
         END IF
      END IF
C
C     Recalculation of quadrature weights in order to enforce exact
C     normalization of the integral (absorption profile * weights)
C
      SUM=0.
      DO I=1,N
         SUM=SUM+PROF(I+IJ0-1)*W(I+IJ0-1)
      END DO
      SUM=S/SUM
      DO I=1,N
         W(I+IJ0-1)=W(I+IJ0-1)*SUM
      END DO
C
C     reset the switch INTMOD
C
      IF(IPROF(ITR).GE.0) THEN
         IF(XMAX.LT.0.) THEN
            INTMOD(ITR)=-2
          ELSE
            INTMOD(ITR)=-1
         END IF
       ELSE
         IF(XMAX.LT.0.) THEN
            INTMOD(ITR)=2
          ELSE
            INTMOD(ITR)=1
         END IF
      END IF
      IF(ABS(INMOD).GE.10) INTMOD(ITR)=0
C
      IF(INDEXP(ITR).NE.0) THEN
         CALL IJALIS(ITR,IFRQ0,IFRQ1)
      END IF
      RETURN
C
  200 CONTINUE
      IF(ISPODF.GT.0) RETURN
      DO 220 IT=1,NTRANS
         IF(.NOT.LINE(IT)) GO TO 220
         IP=IABS(IPROF(IT))
         IF(IP.NE.2) GO TO 220
         IAT=IATM(ILOW(IT))
         AM=BOL2/AMASS(IAT)*TEFF
         DOPP=FR0(IT)*CIN*SQRT(AM+VTB*VTB)
         DOP1=UN/DOPP
         ANE=ELSTD
         IF(ANE.LE.0.) ANE=1.E14
         F000=EXP(TTW*LOG(ANE))
         II=NQUANT(ILOW(IT))
         JJ=NQUANT(IUP(IT))
         IZZ=IZ(IEL(ILOW(IT)))
         FAC=TWO
         F00=F0C1*F000
         IF(IZZ.EQ.2) THEN
            FAC=UN
            F00=F0C2*F000
         END IF
         CALL STARK0(II,JJ,IZZ,XKIJ,WL0,FIJ)
         FXK=F00*XKIJ
         DBETA=WL0*WL0*C18IN/FXK
         BETAD=DOPP*DBETA
         FID=OS0*FIJ*DBETA
         FID0=OS0*(OSC0(IT)-FIJ)*DOP1*PISQ1
         CALL DIVSTR(IZZ)
         DO IJ=IFR0(IT),IFR1(IT)
            BETA=DBETA*ABS(FREQ(IJ)-FR0(IT))
            SG=STARKA(BETA,fac)*FID
            SG0=0.
            V=(FREQ(IJ)-FR0(IT))*DOP1
            IF(ABS(V).LE.13.) SG0=EXP(-V*V)*FID0
            PROF(IJ)=SG+SG0
         END DO
  220 CONTINUE
      RETURN
      END
C
C
C     ****************************************************************
C
C
      SUBROUTINE LINSPL(ITR,DOP,AGAM)
C     ===============================
C
C     Set up depth-independent profile for a line
C     Analog to LINSET used in sampling mode
C
C     Input:
C      ITR    -  index of the transition
C      DOP    -  Doppler width
C      AGAM   -  damping parameter (for lines with Voigt or non-standard
C                profile only)
C
C     Output (to COMMON/FRQEXP)
C      PROF   -  values of absorption profile
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'ATOMIC.FOR'
      INCLUDE 'MODELQ.FOR'
      PARAMETER (OS0=0.02654)
C
      IJ0=IFR0(ITR)
      IJ1=IFR1(ITR)
      N=IJ1-IJ0+1
      KJ0=KFR0(ITR)
      KJ1=KFR1(ITR)
C
C     For all types of integration:
C
      S=OSC0(ITR)*OS0
      IP0=IPROF(ITR)
      IP=IABS(IP0)
      DO I=1,N
         PROF(KJ0+I-1)=PROFIL(FREQ(IJ0+I-1),AGAM,DOP,ITR,IP,0)*S/DOP
      END DO
C
C     for IPROF(ITR) ge 0  -  endpoint(s) of the line profile are forced to
C                             have zero cross-section
C
      IF(IP0.GE.0) THEN
         PROF(KJ0)=0.
         PROF(KJ1)=0.
      END IF
      RETURN
      END
C
C
C     ****************************************************************
C
C

      SUBROUTINE LINPRO(ITR,ID,PRF)
C     =============================
C
C     Line profile coefficient for the line ITR;
C     for "classical" lines (i.e. those which are not represented by ODF's).
C     It is either calculated for each depth (if LCOMP=true), or is
C     just taken as constant, depth-independent quantity, already
C     calculated in START and stored in array PROF
C
C     Input:  ITR    - index of transition
C             ID     - depth index
C     Output: PRF    - array of absorption profile
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'ATOMIC.FOR'
      INCLUDE 'MODELQ.FOR'
      INCLUDE 'ODFPAR.FOR'
      common/quasun/tqmprf,iquasi,nunalp,nunbet,nungam,nunbal
      PARAMETER (BOL2=2.76108D-16, CIN=UN/2.997925D10, OS0=0.02654,
     *           F0C1=1.25D-9, TTW=2./3., PISQ1=UN/1.77245385090551D0,
     *           C18IN=UN/CAS,F0C2=3.906D-11,cca=2.997925e18,
     *           CINV=UN/cca,AL10=2.3025851)
      DIMENSION PRF(MFREQL),PRF0(MHWL),prfb(mfreql),prfr(mfreql)
C
      MODE=IABS(INDEXP(ITR))
      IJ0=IFR0(ITR)
      IJ1=IFR1(ITR)
      INTM0=INTMOD(ITR)
      INTM=IABS(INTM0)
      IP=IABS(IPROF(ITR))
C
C     Doppler width
C
      IF(LCOMP(ITR).OR.IP.EQ.2) THEN
         IAT=IATM(ILOW(ITR))
         AM=BOL2/AMASS(IAT)*TEMP(ID)
         DOP=FR0(ITR)*CIN*SQRT(AM+VTURBS(ID)*VTURBS(ID))
         DOP1=UN/DOP
      END IF
      IF(LCOMP(ITR).AND.IP.EQ.1.OR.IP.LT.0) THEN
         CALL DOPGAM(ITR,ID,TEMP(ID),DOP,AGAM)
         DOP1=UN/DOP
      END IF
      S=OSC0(ITR)*OS0
      XNORM=PISQ1*S*DOP1
C
C     Depth-independent profile
C
      IF(ISPODF.EQ.0) THEN
         DO IJ=IJ0,IJ1
            PRF(IJ-IJ0+1)=PROF(IJ)
         END DO
       ELSE
         DO IJ=KFR0(ITR),KFR1(ITR)
            PRF(IJ-KFR0(ITR)+1)=PROF(IJ)
         END DO
      END IF
      IF(.NOT.LCOMP(ITR)) RETURN
c
c     evaluation of the depth-dependent profile
C
      IF(IP.EQ.0) THEN
         DO IJ=IJ0,IJ1
            V=(FREQ(IJ)-FR0(ITR))*DOP1
            IF(ABS(V).LE.13.) PRF(IJ-IJ0+1)=EXP(-V*V)*XNORM
         END DO
       ELSE IF(IP.EQ.1) THEN
         DO IJ=IJ0,IJ1
            V=(FREQ(IJ)-FR0(ITR))*DOP1
            PRF(IJ-IJ0+1)=VOIGT(V,AGAM)*XNORM
         END DO
       ELSE IF(IP.EQ.2) THEN
C
C      for IP=2 - approximate Stark profile for hydrogen lines, 
C                 evaluate necessary frequency-independent quantities 
C
         ANE=ELEC(ID)
         F000=EXP(TTW*LOG(ANE))
         II=NQUANT(ILOW(ITR))
         JJ=NQUANT(IUP(ITR))
         IZZ=IZ(IEL(ILOW(ITR)))
         FAC=TWO
         F00=F000*F0C1
         if(iquasi.gt.0.and.ii.eq.1) then
            if(jj.eq.2) fac=un
            if(jj.eq.3.and.iquasi.gt.1) fac=un
         end if
         IF(IZZ.EQ.2) THEN
           FAC=UN
           F00=F000*F0C2
         END IF
         CALL STARK0(II,JJ,IZZ,XKIJ0,WL00,FIJ0)
         FXK=F00*XKIJ0
         DBETA=WL00*WL00*C18IN/FXK
         BETAD=DOP*DBETA
         FID=OS0*FIJ0*DBETA
         FID0=OS0*(OSC0(ITR)-FIJ0)*DOP1*PISQ1
         CALL DIVSTR(IZZ)
C
C        loop over frequencies
C
         DO IJ=IJ0,IJ1
            BETA=DBETA*ABS(FREQ(IJ)-FR0(ITR))
            SG=STARKA(BETA,fac)*FID
            SG0=0.
            V=(FREQ(IJ)-FR0(ITR))*DOP1
            IF(ABS(V).LE.13.) SG0=EXP(-V*V)*FID0
            PRF(IJ-IJ0+1)=SG+SG0
         END DO
       ELSE IF(IP.EQ.3.or.ip.eq.4) THEN
C
C     for IP=3 - exact Stark profile for hydrogen lines (Lemke tables)
C
         II=NQUANT(ILOW(ITR))
         JJ=NQUANT(IUP(ITR))
         ANE=ELEC(ID)
         F00=F0C1*EXP(TTW*LOG(ANE))
         FID=OS0*OSC0(ITR)
         anel=log10(ane)
         tl=log10(temp(id))
c
c        switch to either original Lemke/Tremblay of Xenomorph
c
         if(ilxen(ii,jj).eq.0.or.anel.lt.xnemin) then
c
c        original Lemke/Tremblay
c
            WLI0=cca/fr0(itr)
            fac=un
            if(iquasi.gt.0.and.ii.eq.1) then
               if(jj.eq.2) fac=half
               if(jj.eq.3.and.iquasi.gt.1) fac=half
            end if
            iline=0
            if(ip.eq.3) then
               IF(II.LE.4.AND.JJ.LE.22) iline=ilinh(ii,jj)
             else
               IF(II.LE.2.AND.JJ.LE.10) iline=ilinh(ii,jj)
            end if
c                
            if(iline.gt.0) then
            CALL INTLEM(PRF0,WLI0,ILINE,ID)
            NWL=NWLHYD(ILINE)
            DO IJ=IJ0,IJ1
               al=abs(wli0-cca/freq(ij))
               IF(AL.LT.1.E-4) AL=1.E-4
               AL=AL/F00
               AL=LOG10(AL)
               DO IWL=1,NWL-1
                  IW0=IWL
                  IF(AL.LE.WLHYD(ILINE,IWL+1)) GO TO 40
               END DO
   40          IW1=IW0+1
               PRFF=(PRF0(IW0)*(WLHYD(ILINE,IW1)-AL)+PRF0(IW1)*
     *              (AL-WLHYD(ILINE,IW0)))/
     *              (WLHYD(ILINE,IW1)-WLHYD(ILINE,IW0))
               SG=EXP(PRFF*AL10)*FID
               SG=SG*WLI0**2*CINV/F00
              PRF(IJ-IJ0+1)=SG*fac
            END DO
            end if
c
c         XENOMORPH data for selected lines
c
          else
            ixn=ilxen(ii,jj)
            nwl=nwlxen(ixn)
            do iwl=1,nwl
c              call intxen(prfb0,prfr0,tl,anel,iwl,ixn,id)
               prfb(iwl)=prfb0
               prfr(iwl)=prfb0
            end do
            do ij=ij0,ij1
               al=(freq(ij)-fr0(itr))/f00
               if(abs(al).lt.1.e-4) al=1.e-4
               all=log10(abs(al))
               do iwl=1,nwl-1
                  iw0=iwl
                  if(all.le.alxen(ixn,iwl+1)) go to 50
               end do
   50          iw1=iw0+1
               if(al.gt.0.) then
                  prff=(prfb(iw0)*(alxen(ixn,iw1)-all)+
     *                  prfb(iw1)*(all-alxen(ixn,iw0)))/
     *                  (alxen(ixn,iw1)-alxen(ixn,iw0))
                else
                  prff=(prfr(iw0)*(alxen(ixn,iw1)-all)+
     *                  prfr(iw1)*(all-alxen(ixn,iw0)))/
     *                  (alxen(ixn,iw1)-alxen(ixn,iw0))
               end if
               sg=exp(prff*al10)*fid/f00
               prf(ij-ij0+1)=sg
            end do
         END IF 
c
c     for IP ge 10 - special line profile
c
       ELSE IF(IP.GE.10) THEN
         DO IJ=IJ0,IJ1
            PRF(IJ-IJ0+1)=PROFSP(FREQ(IJ),DOP,ITR,ID)
         END DO
      END IF
C
C     if required, force the profile at the endpoints to zero
C
      IF(INTM0.EQ.-1) THEN
         PRF(IJ1-IJ0+1)=0.
       ELSE IF(INTM0.EQ.-2) THEN
         PRF(1)=0.
         PRF(IJ1-IJ0+1)=0.
      END IF
      RETURN
      END
C
C
C
C     ****************************************************************
C
C
 
      FUNCTION SIGK(FR,ITR,MODE)
C     ==========================
C
C     Photoionization cross-sections
C
C     Input: FR  -  frequency
C            ITR -  index of the transition
C            MODE-  =0 - cross-section equal to zero longward of edge
C                   >0 - cross-section non-zero (extrapolated) longward
C                        of edge
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'ATOMIC.FOR'
      PARAMETER (SIH0=2.815D29, E10=2.3025851)
      DIMENSION XFIT(MFIT) , ! local array containing x     for OP data
     +          SFIT(MFIT)   ! local array containing sigma for OP data
C
      PEACH(X,S,A,B)  =A*X**S*(B+X*(1.-B))*1.D-18
      HENRY(X,S,A,B,C)=A*X**S*(C+X*(B-2.*C+X*(1.+C-B)))*1.D-18
C
      SIGK=0.
      IF(INDEXP(ITR).EQ.0) RETURN
      IF(MODE.EQ.0.AND.FR.LT.FR0(ITR)) RETURN
C
C     IC is the index of the given transition in the special numbering
C        of continua (given by ITRCON(ITR)
C     IBF(IC) is the switch controlling the mode of evaluation of the
C        cross-section:
C      = 0  hydrogenic cross-section, with Gaunt factor set to 1
C      = 1  hydrogenic cross-section with exact Gaunt factor
C      = 2  Peach-type expression (see function PEACH)
C      = 3  Henry-type expression (see function HENRY)
C      = 4  Butler fit formula (polynomial fits to the OP results);
C      = 5  Verner fit formula (OP results & HDS calculations at high
C              energies); ONLY for GROUND states;
C      = 6  DETAIL's fit formula from Klaus Werner
C              (similar to Butler's fit but up to 5th order)
C      = 7  hydrogenic cross-section with Gaunt factor after Werner
C      = 9  Opacity project fits (routine TOPBAS - interpolations)
C      > 100 - cross-sections extracted form TOPBASE, for several points
C           In this case, IBF-100 is the number of points
C      < 0  non-standard, user supplied expression (user should update
C           subroutine SPSIGK)
C
C      for H- : for any IBF > 0  - standard expression
C
C      for He I:
C       for IBF = 11 or = 13  -  Opacity Project cross section
C                Seaton-Ferney's cubic fits, Hummer's procedure (HEPHOT)
C           IBF = 11  means that the multiplicity S=1 (singlet)
C           IBF = 13  means that the multiplicity S=3 (triplet)
C       for IBF = 10  - cross section, based on Opacity Project, but
C                       appropriately averaged for an averaged level
C                       (see explanation in SBFHE1)
C
C       for IBF = 21 or 23  Koester's fit (A&A 149, 423) 
C
C           IBF = 21  means that the multiplicity S=1 (singlet)
C           IBF = 23  means that the multiplicity S=3 (triplet)
C
C
c      IC=ITRA(IUP(ITR),ILOW(ITR))
      IC=ITRCON(ITR)
      IB=IBF(IC)
      II=ILOW(ITR)
      IQ=NQUANT(II)
      IE=IEL(ILOW(ITR))
      IF(IB.LT.0) GO TO 60
      IF(IE.EQ.IELHM) GO TO 40
      IF(IE.EQ.IELHE1.AND.IB.GE.10.AND.IB.LE.23) GO TO 50
      CH=IZ(IE)*IZ(IE)
      IQ5=IQ*IQ*IQ*IQ*IQ
C
      IF(IB.EQ.0) THEN
C
C        hydrogenic expression (for IBF = 0)
C
         SIGK=SIH0/FR/FR/FR*CH*CH/IQ5
C
C        exact hydrogenic - with Gaunt factor (for IBF=1)
C
       ELSE IF(IB.EQ.1) THEN
         SIGK=SIH0/FR/FR/FR*CH*CH/IQ5
         SIGK=SIGK*GAUNT(IQ,FR/CH)
       ELSE IF(IB.EQ.2) THEN
C
C        Peach-type formula (for IBF=2)
C   
         FREL=FR0(ITR)/FR
         IF(GAMCS(IC).GT.0) THEN
            IF(GAMCS(IC).LT.1.E6) THEN
              FR00=2.997925E18/GAMCS(IC)
            ELSE
              FR00=GAMCS(IC)
            END IF
            IF(FR.LT.FR00) RETURN        
         FREL=FR00/FR
         END IF 
         if(frel.gt.0.) 
     *   SIGK=PEACH (FREL,S0CS(IC),ALFCS(IC),BETCS(IC))
       ELSE IF(IB.EQ.3) THEN
C
C        Henry-type formula (for IBF=3)
C
         FREL=FR0(ITR)/FR
         if(frel.gt.0.) 
     *   SIGK=HENRY(FREL,S0CS(IC),ALFCS(IC),BETCS(IC),GAMCS(IC))
C
C        Butler expression
C
       ELSE IF(IB.EQ.4) THEN
         FREL=FR0(ITR)/FR
         XL=LOG(FREL)
         SL=S0CS(IC)+XL*(ALFCS(IC)+XL*BETCS(IC))
         SIGK=EXP(SL)
C
C        Verner expression
C
       ELSE IF(IB.EQ.5) THEN
         SIGK=VERNER(FR,ITR)
C
C        DETAIL expression
C
       ELSE IF(IB.EQ.6) THEN
         FREL=FR0(ITR)/FR
         XL=LOG(FREL)
       XL2=XL*XL
       XL3=XL2*XL
         SL=CTOP(1,IC)+XL*CTOP(2,IC)+XL2*CTOP(3,IC)+XL3*CTOP(4,IC)
       SL=SL+XL2*XL2*CTOP(5,IC)+XL3*XL2*CTOP(6,IC)
         SIGK=EXP(SL)
C
C     exact hydrogenic - with Gaunt factor from K Werner (for IBF=7)
C
       ELSE IF(IB.EQ.7) THEN
         IQ5=IQ*IQ*IQ*IQ*IQ
         SIGK=SIH0/(FR*FR*FR)*CH*CH/IQ5*GNTK(IQ,FR/CH)
C
C     selected Opacity Project data (for IBF=9)
C     (c.-s. evaluated by routine TOPBAS which needs an input file RBF.DAT)
C
       ELSE IF(IB.EQ.9) THEN
         SIGK=TOPBAS(FR,FR0(ITR),TYPLEV(II))
C
C     other Opacity Project data (for IBF>100)
C     (c.-s. evaluated by interpolating from direct input data)
C
       ELSE IF(IB.GT.100) THEN
         NFIT=IB-100
         X = LOG10(FR/FR0(ITR))
         SIGM=0.
         IF(X.GE.XTOP(1,IC)) THEN
           DO IFIT = 1,NFIT
              XFIT(IFIT) = XTOP(IFIT,IC)
              SFIT(IFIT) = CTOP(IFIT,IC)
           END DO
           SIGM  = YLINTP (X,XFIT,SFIT,NFIT,MFIT)
           SIGM  = 1.D-18*EXP(E10*SIGM)
         ENDIF
         SIGK=SIGM
      END IF
      RETURN
C
C     special expression for H-
C
   40 SIGK=SBFHMI(FR)
      RETURN
C
C     He I cross-sections
C
   50 SIGK=SBFHE1(II,IB,FR,G(II))
      RETURN
C
C     non-standard, user supplied form of cross-section (for IBF < 0)
C
   60 CALL SPSIGK(IB,FR,SIGSP)
      SIGK=SIGSP
      RETURN
      END
C
C
C     ****************************************************************
C
C
      FUNCTION VERNER(FR,ITR)
C     =======================
C
C     Photoionization cross-sections for ground states
C       of atoms and ions. Analytical fits from:
C           Verner D.A. et al. 1996, ApJ 465
C           Verner & Yakovlev 1995, A&AS 109, 125
C
C     10-July-1996: Version for H to Si, S, Ar, Ca and Fe.
C          No test on threshold energy as given by Verner et al.
C          to avoid inconsistencies with limits FR0 as read by Tlusty
C
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'ATOMIC.FOR'
C
      PARAMETER (HHEV=H/1.6022D-12,T18=1.D-18,MVER=105)
      DIMENSION S0(MVER),E0(MVER),EMX(MVER)
      DIMENSION Y0(MVER),Y1(MVER),YW(MVER),YA(MVER),PV(MVER)
      DIMENSION S95(MVER),E95(MVER),Y95(MVER),YW95(MVER),P95(MVER)
      DIMENSION IV0(14)
C
      DATA IV0/0,1,3,6,10,15,21,28,36,45,55,66,78,91/
      DATA S0/5.475D4,9.492D2,1.369D4,6.245D1,3.201D2,6.083D3,2.932D5,
     +        2.678D2,5.458D2,3.422D3,5.466D0,1.859D4,5.393D1,2.846D2,
     +        2.190D3,5.027D2,8.709D0,1.539D4,1.068D2,2.344D2,1.521D3,
     +        8.235D2,1.944D0,9.375D-1,1.690D4,8.376D1,1.519D2,1.117D3,
     +        1.745D3,5.967D1,6.753D2,8.659D-1,1.642D4,6.864D1,1.329D2,
     +        8.554D2,3.803D3,8.013D1,1.541D2,3.165D0,4.690D-1,1.157D4,
     +        6.930D1,1.039D2,6.759D2,4.287D3,1.583D3,5.708D0,1.685D3,
     +        2.430D0,9.854D-1,1.198D4,5.631D1,6.695D1,5.475D2,1.601D0,
     +        1.04D3,1.885D3,2.33D0,2.346D1,7.101D1,1.609D0,7.215D3,
     +        4.729D1,3.995D1,4.525D2,1.372D8,3.278D0,5.377D2,1.394D3,
     +        1.728D0,2.185D0,3.104D0,6.344D-2,9.008D2,4.427D1,6.14D1,
     +        3.802D2,7.195D0,6.948D-2,4.915D0,1.513D1,2.925D2,1.962D-2,
     +        1.889D2,2.35D-1,4.982D-1,1.774D4,3.388D1,4.036D1,3.239D2,
     +        2.506D1,4.14D0,5.79D-4,6.083D0,8.863D-1,7.293D1,6.68D-2,
     +        3.477D-1,1.465D-1,1.95D-1,1.992D4,2.539D1,4.754D1,2.793D2/
      DATA E0/0.4298,13.61,1.72,3.107,20.06,3.871,9.539,1.181,17.6,
     +        6.879,0.5213,2.869,1.041,33.36,10.75,2.144,0.4058,4.614,
     +        3.506,46.24,15.48,4.034,0.06128,0.242,5.494,4.471,69.43,
     +        21.08,1.24,1.386,0.1723,0.2044,2.854,7.824,87.09,27.54,
     +        12.97,1.763,2.542,0.7744,0.7286,4.008,2.563,113.1,34.85,
     +        4.87,12.47,0.7753,5.566,1.248,1.499,4.888,10.03,158.6,
     +        43.04,6.139,8.203,10.69,0.669,5.408,48.46,2.096,153.5,
     +        13.91,226.8,52.11,11.97,8.139,10.86,29.12,0.9762,1.711,
     +        3.57,0.4884,34.82,14.52,204.2,62.03,13.81,0.2048,10.27,
     +        3.13,24.14,0.3483,2.636,0.4866,1.842,8.044,23.55,273.8,
     +        72.81,23.17,2.556,0.1659,1.288,0.7761,63.05,0.3277,0.7655,
     +        0.3343,0.8787,12.05,35.6,275.2,84.47/
      DATA EMX/5.D4,5.D4,5.D4,64.39,5.D4,5.D4,11.93,12.99,5.D4,5.D4,
     +         194.0,209.8,227.4,5.D4,5.D4,291.0,307.6,328.9,352.2,
     +         5.D4,5.D4,404.8,423.6,447.3,475.3,504.3,5.D4,5.D4,
     +         538.0,558.1,584.0,614.4,649.1,683.7,5.D4,5.D4,694.0,
     +         712.2,739.2,770.9,809.1,850.2,890.5,5.D4,5.D4,870.1,
     +         883.1,913.1,948.0,987.3,1031.,1078.,1125.,5.D4,5.D4,
     +         38.14,1074.,1110.,1143.,1185.,1230.,1281.,1335.,1386.,
     +         5.D4,5.D4,54.9,65.69,1317.,1356.,1400.,1449.,1501.,
     +         1558.,1618.,1675.,5.D4,5.D4,80.4,89.97,102.6,1588.,
     +         1634.,1688.,1739.,1799.,1862.,1929.,1992.,5.D4,5.D4,
     +         106.,118.6,131.1,146.6,1887.,1946.,2001.,2058.,2125.,
     +         2194.,2268.,2336.,5.D4,5.D4/
      DATA Y0/0.,0.4434,0.,0.,0.,0.,8.278D-4,0.,0.,0.,1.319D1,4.96D-3,
     +        0.,0.,0.,1.133,4.929D1,4.378D-3,0.,0.,0.,8.598D-1,4.28D2,
     +        1.877D2,6.415D-3,0.,0.,0.,8.698,2.131D1,3.839D-3,3.328D2,
     +        3.036D-2,0.,0.,0.,1.701D-4,1.715D1,1.641D1,9.531D1,
     +        1.506D2,2.071D-2,0.,0.,0.,4.236D-2,1.52,7.654D1,5.149,
     +        9.169D1,1.042D2,2.536D-2,0.,0.,0.,0.,3.375,3.725,1.383D2,
     +        2.204D1,9.603D-4,9.94D1,0.1667,0.,0.,0.,0.,0.,4.86,
     +        9.402D-1,1.276D2,1.007D2,5.452D1,5.348D2,5.444,0.,0.,0.,
     +        2.041D-1,9.149D1,0.,3.994D1,3.495,5.675D-2,3.552D1,
     +        5.704D2,1.719D2,2.538D-2,0.,0.,0.,1.672D-5,6.634,9.613D1,
     +        0.,2.009D2,1.115,1.149D-2,3.85D2,1.036D3,4.528D2,1.99D-2,
     +        0.,0.,0./
      DATA Y1/0.,2.136,0.,0.,0.,0.,1.269D-2,0.,0.,0.,4.556,3.4D-2,0.,
     +        0.,0.,1.607,3.234,2.528D-2,0.,0.,0.,2.325,2.03D1,3.999,
     +        1.937D-2,0.,0.,0.,1.271D-1,1.503D-2,4.569D-1,4.285D1,
     +        5.554D-2,0.,0.,0.,1.345D-2,7.724D-1,5.124,9.781,2.574D-1,
     +        3.998D-2,0.,0.,0.,5.873,1.084D-1,2.023,6.687,3.702D-1,
     +        1.435,4.417D-2,0.,0.,0.,0.,4.01,0.2279,4.26,0.7577,
     +        6.378D-3,3.278,1.766D-2,0.,0.,0.,0.,0.,3.722,0.1135,
     +        3.979,1.729,2.078,3.997D-3,7.918D-2,0.,0.,0.,0.4753,
     +        0.6565,0.,4.803,0.2701,0.2768,8.223D-3,0.158,6.595,
     +        1.203D-2,0.,0.,0.,0.4207,0.1272,0.6442,0.,4.537,8.051D-2,
     +        0.6396,8.999D-2,0.2936,1.015,1.007D-2,0.,0.,0./
      DATA YW/0.,2.039,0.,0.,0.,0.,3.655D-1,0.,0.,0.,1.887D1,3.503,0.,
     +        0.,0.,9.157D-2,2.093,5.922,0.,0.,0.,9.097D-2,1.043D1,
     +        1.85,7.904,0.,0.,0.,7.589D-2,1.934D-2,1.191D-1,3.143,
     +        2.836D1,0.,0.,0.,2.17D-3,5.103D-1,1.115,6.812,2.57D-1,
     +        2.411D1,0.,0.,0.,2.434D-1,6.558D-2,4.633D-1,8.29D-3,
     +        6.855D-1,1.656,2.811D1,0.,0.,0.,0.,2.328,8.579D-2,0.7365,
     +        0.9275,1.285D-2,1.895,0.9121,0.,0.,0.,0.2805,0.,2.604,
     +        4.326D-2,0.809,0.6325,1.422,0.6666,2.751,0.,0.,0.,0.3166,
     +        0.4615,0.,5.342,0.1,8.839,1.836,0.2773,0.6945,29.53,0.,
     +        0.,0.,0.2837,1.57,0.8626,0.,1.303,2.989D-3,3.28,1.476D-3,
     +        1.646,0.4489,2.392D1,0.,0.,0./
      DATA YA/3.288D1,1.469D0,3.288D1,1.501D1,7.391D0,3.288D1,4.301D-1,
     +        5.645D0,1.719D1,3.288D1,8.618D0,1.783D0,1.767D1,2.163D1,
     +        3.288D1,6.216D1,1.261D2,1.737D0,1.436D1,2.183D1,3.288D1,
     +        8.033D1,8.163D2,2.788D2,1.714D0,3.297D1,2.627D1,3.288D1,
     +        3.784D0,3.175D1,3.852D2,4.931D2,1.792D0,3.210D1,2.535D1,
     +        3.288D1,2.587D0,1.667D1,5.742D1,1.099D2,1.400D2,1.848D0,
     +        7.547D1,2.657D1,3.288D1,5.798D0,3.935D0,6.725D1,6.409D2,
     +        1.066D2,1.350D2,1.788D0,3.628D1,3.352D1,3.288D1,6.148D3,
     +        8.259D0,3.613D0,1.205D2,2.913D1,3.945D1,2.473D2,3.886D-1,
     +        3.889D1,5.315D1,3.288D1,2.228D-1,4.341D7,9.779D0,2.895D0,
     +        9.184D1,9.350D1,6.060D1,5.085D2,1.823D0,3.826D1,2.778D1,
     +        3.288D1,1.621D3,5.675D2,1.990D6,1.674D1,6.973D0,1.856D1,
     +        1.338D2,7.216D2,2.568D2,1.653D0,3.432D1,3.567D1,3.288D1,
     +        2.057D1,1.337D1,1.474D2,1.356D6,1.541D2,1.558D2,4.132D1,
     +        3.733D2,1.404D3,7.461D2,1.582D0,3.307D1,2.848D1,3.288D1/
      DATA PV/2.963,31.88,2.963,4.895,2.916,2.963,10.52,11.7,3.157,
     +        2.963,17.28,16.18,9.54,2.624,2.963,5.101,8.578,15.93,
     +        7.457,2.581,2.963,3.928,8.773,9.156,17.06,6.003,2.315,
     +        2.963,17.64,8.943,6.822,8.785,26.47,5.495,2.336,2.963,
     +        7.275,10.5,6.614,9.203,9.718,24.46,6.448,2.255,2.963,
     +        8.355,7.81,10.05,3.056,8.999,8.836,25.5,5.585,2.002,
     +        2.963,3.839,7.362,9.803,9.714,8.26,2.832,7.681,8.476,
     +        5.265,1.678,2.963,15.74,3.61,7.117,6.487,10.06,9.202,
     +        8.857,9.385,14.44,5.46,2.161,2.963,3.642,9.049,3.477,
     +        11.8,6.724,20.84,6.204,8.659,8.406,26.55,5.085,1.915,
     +        2.963,3.546,11.91,13.36,3.353,9.98,2.4,16.06,8.986,
     +        8.503,8.302,24.25,4.728,2.135,2.963/
      DATA S95/5.475D4,4.47D3,1.369D4,1.564D2,3.201D2,6.083D3,1.306D2,
     +         9.796D1,5.458D2,3.422D3,9.698D1,1.037D2,8.605D1,2.846D2,
     +         2.19D3,7.421D1,6.649D1,8.067D1,8.111D1,2.344D2,1.521D3,
     +         4.748D1,5.002D1,5.235D1,7.046D1,7.304D1,1.519D2,1.117D3,
     +         3.237D1,3.584D1,3.939D1,4.123D1,5.735D1,6.029D1,1.329D2,
     +         8.554D2,2.295D1,4.798D1,3.144D1,5.302D1,3.680D1,4.668D1,
     +         4.890D1,1.039D2,6.759D2,1.664D1,2.783D1,2.943D1,3.027D1,
     +         3.097D1,3.232D1,5.011D1,3.719D1,6.695D1,5.475D2,2.486D2,
     +         1.889D1,3.466D1,2.551D1,3.753D1,2.654D1,2.745D1,3.613D1,
     +         3.850D1,3.995D1,4.525D2,2.023D2,2.049D2,1.877D1,2.212D1,
     +         2.350D1,3.278D1,3.407D1,2.721D1,3.699D1,3.290D1,6.140D1,
     +         3.802D2,1.735D2,1.842D2,1.990D2,2.494D1,2.572D1,2.786D1,
     +         2.653D1,2.989D1,2.857D1,3.223D1,3.272D1,4.036D1,3.239D2,
     +         1.532D2,1.832D2,2.197D2,2.326D2,3.100D1,2.398D1,2.474D1,
     +         3.135D1,2.641D1,3.343D1,2.832D1,3.880D1,4.754D1,2.793D2/
      DATA E95/0.4298,5.996,1.72,27.4,20.06,3.871,40.93,47.59,17.6,
     +         6.879,61.55,59.84,65.92,33.36,10.75,86.55,91.13,83.7,
     +         84.12,46.24,15.48,127.,124.2,122.,107.,106.,69.43,21.08,
     +         177.4,169.,162.,159.3,137.7,135.4,87.09,27.54,239.,166.,
     +         205.5,160.3,192.5,173.7,171.1,113.1,34.85,314.4,245.2,
     +         238.7,236.,234.3,230.8,188.6,219.3,158.6,43.04,36.55,
     +         326.6,242.6,284.1,236.,280.6,277.7,246.4,240.6,226.8,
     +         52.11,49.37,49.4,360.7,334.5,325.9,277.7,274.2,308.2,
     +         267.1,285.8,204.2,62.03,64.45,61.54,60.16,348.2,343.2,
     +         326.8,339.6,319.3,330.5,311.7,313.7,273.8,72.81,78.08,
     +         71.54,66.52,64.82,346.,379.9,375.9,344.2,367.9,335.6,
     +         359.9,315.,275.2,84.47/
      DATA Y95/3.288D1,2.199D0,3.288D1,3.382D1,7.391D0,3.288D1,1.212D2,
     +         1.166D3,1.719D1,3.288D1,7.354D1,7.915D1,1.906D2,2.163D1,
     +         3.288D1,5.498D1,9.609D1,7.471D1,7.459D1,2.183D1,3.288D1,
     +         1.380D2,9.100D1,9.428D1,5.342D1,5.547D1,2.627D1,3.288D1,
     +         3.812D2,1.894D2,1.104D2,1.141D2,5.486D1,5.682D1,2.535D1,
     +         3.288D1,1.257D3,5.000D1,1.230D2,5.000D1,7.933D1,5.876D1,
     +         6.137D1,2.657D1,3.288D1,2.042D5,6.075D1,6.525D1,6.734D1,
     +         6.907D1,7.167D1,5.000D1,8.181D1,3.352D1,3.288D1,3.222D2,
     +         2.527D2,5.000D1,6.977D1,5.000D1,7.237D1,7.512D1,4.968D1,
     +         5.198D1,5.315D1,3.288D1,1.079D4,4.112D4,1.401D2,7.147D1,
     +         6.173D1,5.000D1,5.000D1,5.108D1,5.000D1,5.384D1,2.778D1,
     +         3.288D1,1.131D4,2.404D3,5.082D2,3.260D1,3.376D1,5.000D1,
     +         3.478D1,5.000D1,3.660D1,5.000D1,4.295D1,3.567D1,3.288D1,
     +         5.765D6,3.537D2,1.169D2,1.188D2,1.979D1,5.000D1,5.000D1,
     +         2.260D1,5.000D1,2.487D1,5.000D1,3.034D1,2.848D1,3.288D1/
      DATA YW95/55*0.,0.1465,10*0.,1.463D-2,2.223D-5,10*0.,2.337D-2,
     +          7.839D-3,2.016D-2,10*0.,2.774D-4,2.87D-4,8.658D-4,
     +          8.417D-4,10*0./
      DATA P95/2.963,6.098,2.963,1.49,2.916,2.963,1.348,1.022,3.157,
     +         2.963,1.438,1.436,1.21,2.624,2.963,1.503,1.338,1.442,
     +         1.428,2.581,2.963,1.252,1.335,1.335,1.552,1.538,2.315,
     +         2.963,1.083,1.185,1.289,1.287,1.540,1.533,2.336,2.963,
     +         0.9638,1.65,1.263,1.65,1.377,1.511,1.501,2.255,2.963,
     +         0.845,1.42,1.424,1.422,1.416,1.411,1.65,1.396,2.002,
     +         2.963,3.57,1.12,1.65,1.414,1.65,1.405,1.397,1.579,
     +         1.575,1.678,2.963,2.96,2.995,1.238,1.413,1.464,1.65,
     +         1.65,1.539,1.65,1.56,2.161,2.963,2.762,2.92,3.11,1.732,
     +         1.73,1.65,1.721,1.65,1.712,1.65,1.676,1.915,2.963,2.639,
     +         3.133,3.529,3.547,2.094,1.65,1.65,2.02,1.65,1.982,1.65,
     +         1.925,2.135,2.963/
C
      VERNER=0.
      E=HHEV*FR
C
      II=ILOW(ITR)
      N1=NFIRST(IEL(II))
      IF(II.NE.N1) 
     +   CALL QUIT('Verner fits only for ground states',ii,n1)
      IAT=NUMAT(IATM(II))
      IZZ=IZ(IEL(II))
      IF(IAT.GT.14) GO TO 10
      IVER=IV0(IAT)+IZZ
C
C     1996 Expression
C
      IF(E.LT.EMX(IVER)) THEN
        XX=E/E0(IVER)-Y0(IVER)
        YY=SQRT(XX*XX+Y1(IVER)*Y1(IVER))
        AA=(XX-UN)*(XX-UN)+YW(IVER)*YW(IVER)
        BB=YY**(HALF*PV(IVER)-5.5)
        CC=(UN+SQRT(YY/YA(IVER)))**PV(IVER)
        FY=AA*BB/CC
        VERNER=S0(IVER)*T18*FY
      ELSE
C
C     1995 Expression for high energies
C          (ionization of inner shell electron)
C
        YY=E/E95(IVER)
        XL=0.
        IF((IAT-IZZ).GE.10) XL=UN
        Q=HALF*P95(IVER)-5.5-XL
        AA=(YY-UN)*(YY-UN)+YW95(IVER)*YW95(IVER)
        BB=YY**Q
        CC=(UN+SQRT(YY/Y95(IVER)))**P95(IVER)
        FY=AA*BB/CC
        VERNER=S95(IVER)*T18*FY
      END IF
      RETURN
C
C     Heavier elements
C
   10 IF(IAT.EQ.26) THEN
        VERNER=VERN26(E,IZZ)
      ELSE IF(IAT.EQ.16) THEN
        VERNER=VERN16(E,IZZ)
      ELSE IF(IAT.EQ.18) THEN
        VERNER=VERN18(E,IZZ)
      ELSE IF(IAT.EQ.20) THEN
        VERNER=VERN20(E,IZZ)
      ELSE
        CALL QUIT('VERNER - No data for this element',iat,izz)
      ENDIF
C
      RETURN
      END
C
C
C     ****************************************************************
C
C
      FUNCTION VERN26(E,IZZ)
C     ======================
C
C     Photoionization cross-sections for ground states
C       of all Fe ions.
C           Verner D.A. et al. 1996, ApJ 465
C           Verner & Yakovlev 1995, A&AS 109, 125
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      PARAMETER (T18=1.D-18,MVER=26)
      DIMENSION S0(MVER),E0(MVER),EMX(MVER)
      DIMENSION Y0(MVER),Y1(MVER),YW(MVER),YA(MVER),PV(MVER)
      DIMENSION S95(MVER),E95(MVER),Y95(MVER),YW95(MVER),P95(MVER)
C
      DATA S0/3.062D-1,4.365D3,6.107D0,3.653D2,1.523D-3,5.259D-1,2.42D4,
     +        1.979D1,2.687D1,6.470D1,3.281D0,1.738D0,2.791D-3,1.454D-1,
     +        2.108D2,1.207D1,1.452D0,2.388D0,6.066D-5,4.455D-1,1.098D1,
     +        7.204D-2,2.580D4,1.276D1,1.195D1,8.099D1/
      DATA E0/0.05461,0.1761,0.1698,25.44,0.7256,2.656,5.059,0.07098,
     +        6.741,68.86,8.284,6.295,0.1317,0.8509,0.05555,28.73,
     +        0.3444,31.9,7.519D-4,20.11,9.243,9.713,45.75,73.26,
     +        1057.,293.2/
      DATA EMX/66.,76.17,87.05,106.7,128.8,152.7,178.3,205.5,921.1,
     +         959.,998.3,1039.,1081.,1125.,1181.,1216.,7651.,7769.,
     +         7918.,8041.,8184.,8350.,8484.,8638.,5.D4,5.D4/
      DATA Y0/1.382D2,9.272D1,1.76D2,0.,8.871D1,3.361D1,0.4546,2.542D3,
     +        2.494D1,1.19D-5,2.971D1,4.671D1,2.17D3,4.505D2,2.706D-4,
     +        0.,2.891D1,3.805D1,1.915D6,6.847D1,4.446D1,1.702D2,
     +        3.582D-2,0.,0.,0./
      DATA Y1/0.2481,1.075D2,1.847D1,0.,5.28D-2,3.743D-3,2.683D1,
     +        4.672D2,8.251,6.57D-3,0.522,0.1425,6.852D-3,2.504,
     +        1.628,0.,3.404,0.4805,3.14D1,3.989,3.512,4.263,8.712D-3,
     +        0.,0.,0./
      DATA YW/2.069D1,1.141D1,8.698,0.5602,5.064D1,1.558D1,2.516D-3,
     +        2.158D2,2.387D-4,2.778D-4,0.3279,0.3096,0.6938,0.4937,
     +        1.885D-3,0.,1.264,2.902D-2,4.398,2.757,1.748,9.551D-3,
     +        2.723D1,0.,0.,0./
      DATA YA/2.671D7,6.298D3,1.555D3,8.913D0,3.736D1,1.450D1,4.850D4,
     +        1.745D4,1.807D2,2.062D1,5.360D1,1.130D2,2.487D3,1.239D3,
     +        2.045D4,5.150D2,3.960D2,2.186D1,1.606D6,4.236D1,7.637D1,
     +        1.853D2,1.358D0,4.914D1,5.769D1,3.288D1/
      DATA PV/7.923,5.204,8.055,6.538,17.67,16.32,2.374,6.75,6.29,
     +        4.111,8.571,8.037,9.791,8.066,6.033,3.846,10.13,9.589,
     +        8.813,9.724,7.962,8.843,26.04,4.941,1.718,2.963/
      DATA S95/6.298D1,4.624D1,4.422D1,4.81D1,5.143D1,5.246D1,5.21D1,
     +         5.336D1,2.205D2,2.392D2,2.449D2,3.325D2,3.316D2,3.367D2,
     +         1.496D2,3.383D2,1.15D1,8.327D0,1.155D1,8.619D0,8.773D0,
     +         1.181D1,9.098D1,1.157D1,1.195D1,8.099D1/
      DATA E95/76.3,77.5,77.77,76.25,72.73,72.6,74.33,75.56,171.5,
     +         164.7,163.2,138.3,139.2,138.7,213.6,139.6,1067.,1249.,
     +         1068.,1235.,1228.,1066.,1215.,1087.,1057.,293.2/
      DATA Y95/1.479D1,2.155D1,2.336D1,2.286D1,2.428D1,2.751D1,3.306D1,
     +         3.855D1,5.298D1,5.276D1,5.452D1,5.09D1,5.237D1,5.279D1,
     +         7.000D1,5.459D1,3.412D1,5.000D1,3.578D1,5.000D1,5.000D1,
     +         4.116D1,5.000D1,5.086D1,5.769D1,3.288D1/
      DATA YW95/0.2646,0.2599,0.2557,0.2449,0.1365,0.02105,0.02404,
     +          0.02667,1.508D-5,1.574D-5,1.594D-4,1.114D-5,1.107D-5,
     +          1.111D-5,0.1,1.179D-5,10*0./ 
      DATA P95/7.672,7.138,7.017,7.043,7.028,6.823,6.509,6.265,4.154,
     +         4.204,4.187,4.446,4.41,4.407,3.7,4.366,1.922,1.65,1.895,
     +         1.65,1.65,1.827,1.65,1.722,1.718,2.963/
C
      VERN26=0.
      IVER=IZZ
C
C     1996 Expression
C
      IF(E.LT.EMX(IVER)) THEN
        XX=E/E0(IVER)-Y0(IVER)
        YY=SQRT(XX*XX+Y1(IVER)*Y1(IVER))
        AA=(XX-UN)*(XX-UN)+YW(IVER)*YW(IVER)
        BB=YY**(HALF*PV(IVER)-5.5)
        CC=(UN+SQRT(YY/YA(IVER)))**PV(IVER)
        FY=AA*BB/CC
        VERN26=S0(IVER)*T18*FY
      ELSE
C
C     1995 Expression for high energies
C          (ionization of inner shell electron)
C
        YY=E/E95(IVER)
        XL=0.
        IF(IZZ.LE.16) XL=UN
        Q=HALF*P95(IVER)-5.5-XL
        AA=(YY-UN)*(YY-UN)+YW95(IVER)*YW95(IVER)
        BB=YY**Q
        CC=(UN+SQRT(YY/Y95(IVER)))**P95(IVER)
        FY=AA*BB/CC
        VERN26=S95(IVER)*T18*FY
      END IF
C
      RETURN
      END
C
C
C     ****************************************************************
C
C
      FUNCTION VERN16(E,IZZ)
C     ======================
C
C     Photoionization cross-sections for ground states
C       of all Sulfur ions.
C           Verner D.A. et al. 1996, ApJ 465
C           Verner & Yakovlev 1995, A&AS 109, 125
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      PARAMETER (T18=1.D-18,MVER=16)
      DIMENSION S0(MVER),E0(MVER),EMX(MVER)
      DIMENSION Y0(MVER),Y1(MVER),YW(MVER),YA(MVER),PV(MVER)
      DIMENSION S95(MVER),E95(MVER),Y95(MVER),YW95(MVER),P95(MVER)
C
      DATA S0/4.564D4,3.136D2,6.666D0,2.606D0,5.072D-4,9.139D0,5.703D-1,
     +        3.161D1,9.646D3,5.364D1,1.275D1,3.49D-1,2.294D4,2.555D1,
     +        2.453D1,2.139D2/
      DATA E0/18.08,8.787,2.027,2.173,0.1713,14.13,0.3757,14.62,0.1526,
     +        10.4,6.485,2.443,14.74,33.1,439.,110.4/
      DATA EMX/170.,184.6,199.5,216.4,235.,255.7,2569.,2641.,2705.,
     +         2782.,2859.,2941.,3029.,3107.,5.D4,5.D4/
      DATA Y0/0.9935,2.782,15.68,19.75,94.24,0.,222.2,18.69,1.615D-3,
     +        17.75,34.26,227.9,2.203D-2,0.,0.,0./
      DATA Y1/0.2486,0.1788,9.421,3.361,0.6265,0.,4.606,0.3037,0.4049,
     +        1.663,0.137,1.172,1.073D-2,0.,0.,0./
      DATA YW/0.6385,0.7354,4.109,1.863,0.788,0.,1.503,1.153D-3,1.492,
     +        2.31,1.678,0.7033,27.38,0.,0.,0./
      DATA YA/1.,3.442,54.54,66.41,198.6,1656.,146.,16.11,1438.,36.41,
     +        65.83,541.1,1.529,38.21,44.05,32.88/
      DATA PV/13.61,12.81,8.611,8.655,13.07,3.626,11.35,8.642,5.977,
     +        7.09,7.692,7.769,25.68,5.037,1.765,2.963/
      DATA S95/1.883D2,1.896D2,1.780D2,2.037D2,2.919D2,4.712D2,1.916D1,
     +         1.931D1,1.946D1,2.041D1,2.101D1,2.087D1,2.233D1,2.293D1,
     +         2.453D1,2.139D2/
      DATA E95/91.52,90.58,92.46,87.44,74.11,57.47,495.2,489.1,493.7,
     +         480.2,475.8,482.8,466.9,466.7,439.,110.4/
      DATA Y95/71.93,75.38,149.8,93.1,48.64,36.1,35.55,50.,35.68,50.,
     +         50.,37.42,50.,44.59,44.05,32.88/
      DATA YW95/0.2485,0.2934,0.02142,9.497D-3,0.02785,0.0248,10*0./
      DATA P95/3.633,3.635,3.319,3.565,4.142,4.742,1.742,1.65,1.737,
     +         1.65,1.65,1.72,1.65,1.668,1.765,2.963/
C
      VERN16=0.
      IVER=IZZ
C
C     1996 Expression
C
      IF(E.LT.EMX(IVER)) THEN
        XX=E/E0(IVER)-Y0(IVER)
        YY=SQRT(XX*XX+Y1(IVER)*Y1(IVER))
        AA=(XX-UN)*(XX-UN)+YW(IVER)*YW(IVER)
        BB=YY**(HALF*PV(IVER)-5.5)
        CC=(UN+SQRT(YY/YA(IVER)))**PV(IVER)
        FY=AA*BB/CC
        VERN16=S0(IVER)*T18*FY
      ELSE
C
C     1995 Expression for high energies
C          (ionization of inner shell electron)
C
        YY=E/E95(IVER)
        XL=0.
        IF(IZZ.LE.6) XL=UN
        Q=HALF*P95(IVER)-5.5-XL
        AA=(YY-UN)*(YY-UN)+YW95(IVER)*YW95(IVER)
        BB=YY**Q
        CC=(UN+SQRT(YY/Y95(IVER)))**P95(IVER)
        FY=AA*BB/CC
        VERN16=S95(IVER)*T18*FY
      END IF
C
      RETURN
      END
C
C
C     ****************************************************************
C
C
      FUNCTION VERN18(E,IZZ)
C     ======================
C
C     Photoionization cross-sections for ground states
C       of all Argon ions.
C           Verner D.A. et al. 1996, ApJ 465
C           Verner & Yakovlev 1995, A&AS 109, 125
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      PARAMETER (T18=1.D-18,MVER=18)
      DIMENSION S0(MVER),E0(MVER),EMX(MVER)
      DIMENSION Y0(MVER),Y1(MVER),YW(MVER),YA(MVER),PV(MVER)
      DIMENSION S95(MVER),E95(MVER),Y95(MVER),YW95(MVER),P95(MVER)
C
      DATA S0/2.106D1,2.503D1,3.58D1,2.035D1,9.946D0,1.080D0,3.693D0,
     +        3.295D1,8.279D-1,8.204D0,1.76D3,7.018D-1,2.459D-2,
     +        4.997D-2,2.571D4,2.135D1,3.108D1,1.69D2/
      DATA E0/17.09,24.94,14.17,6.953,10.31,0.544,0.02966,3.844,
     +        0.1926,10.4,0.1257,5.31,0.3209,1.557,18.88,41.54,
     +        446.8,139.9/
      DATA EMX/249.2,266.2,280.1,298.7,320.,342.6,366.7,392.5,3361.,
     +         3446.,3523.,3613.,3702.,3798.,3898.,3988.,5.D4,5.D4/
      DATA Y0/1.688,0.9299,2.384,7.501,6.406,1.7D2,4.383D-4,0.,38.14,
     +        38.04,3.286D-3,1.099D2,2.068D3,4.552D2,2.445D-2,0.,0.,0./
      DATA Y1/0.8943,0.7195,1.794,0.1806,3.659D-3,15.87,2.513,0.,4.649,
     +        0.639,0.3226,0.2202,21.13,6.459,1.054D-2,0.,0.,0./
      DATA YW/0.4185,0.5108,0.6316,0.8842,0.4885,11.07,1.363D-2,0.,
     +        1.434,9.203D-4,1.975,0.4987,0.6692,0.2938,29.09,0.,0.,0./
      DATA YA/2.645D2,1.272D2,3.776D1,1.4D1,7.444D1,9.419D2,9.951D3,
     +        7.082D2,2.392D2,1.495D1,1.579D3,1.001D2,2.285D3,5.031D2,
     +        1.475D0,4.118D1,3.039D1,3.288D1/
      DATA PV/4.796,4.288,5.742,9.595,6.261,7.582,7.313,4.645,11.21,
     +        11.15,6.714,8.939,8.81,8.966,26.34,4.945,2.092,2.963/
      DATA S95/8.372D1,1.937D2,2.281D2,2.007D2,2.474D2,2.786D2,3.204D2,
     +         4.198D2,2.931D1,1.585D1,2.796D1,1.666D1,2.888D1,2.874D1,
     +         2.883D1,3.003D1,3.108D1,1.69D2/
      DATA E95/164.7,108.5,102.5,107.3,98.35,92.33,85.63,73.68,467.9,
     +         612.6,478.9,602.8,473.1,474.9,475.6,468.,466.8,139.9/
      DATA Y95/54.52,70.,43.8,70.,42.84,42.2,42.3,44.19,17.44,50.,
     +         19.17,50.,20.42,22.35,26.15,28.54,30.39,32.88/
      DATA YW95/0.627,0.1,7.167D-3,0.1,7.283D-3,7.408D-3,7.258D-3,
     +          7.712D-3,10*0./
      DATA P95/3.328,3.7,4.046,3.7,4.125,4.227,4.329,4.492,2.362,
     +         1.65,2.271,1.65,2.234,2.171,2.074,2.037,2.092,2.963/
C
      VERN18=0.
      IVER=IZZ
C
C     1996 Expression
C
      IF(E.LT.EMX(IVER)) THEN
        XX=E/E0(IVER)-Y0(IVER)
        YY=SQRT(XX*XX+Y1(IVER)*Y1(IVER))
        AA=(XX-UN)*(XX-UN)+YW(IVER)*YW(IVER)
        BB=YY**(HALF*PV(IVER)-5.5)
        CC=(UN+SQRT(YY/YA(IVER)))**PV(IVER)
        FY=AA*BB/CC
        VERN18=S0(IVER)*T18*FY
      ELSE
C
C     1995 Expression for high energies
C          (ionization of inner shell electron)
C
        YY=E/E95(IVER)
        XL=0.
        IF(IZZ.LE.8) XL=UN
        Q=HALF*P95(IVER)-5.5-XL
        AA=(YY-UN)*(YY-UN)+YW95(IVER)*YW95(IVER)
        BB=YY**Q
        CC=(UN+SQRT(YY/Y95(IVER)))**P95(IVER)
        FY=AA*BB/CC
        VERN18=S95(IVER)*T18*FY
      END IF
C
      RETURN
      END
C
C
C     ****************************************************************
C
C
      FUNCTION VERN20(E,IZZ)
C     ======================
C
C     Photoionization cross-sections for ground states
C       of all Calcium ions.
C           Verner D.A. et al. 1996, ApJ 465
C           Verner & Yakovlev 1995, A&AS 109, 125
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      PARAMETER (T18=1.D-18,MVER=20)
      DIMENSION S0(MVER),E0(MVER),EMX(MVER)
      DIMENSION Y0(MVER),Y1(MVER),YW(MVER),YA(MVER),PV(MVER)
      DIMENSION S95(MVER),E95(MVER),Y95(MVER),YW95(MVER),P95(MVER)
C
      DATA S0/5.37D5,1.064D7,3.815D1,7.736D0,1.523D-1,7.642D1,4.76D-1,
     +        6.641D-1,2.076D2,1.437D1,9.384D-1,1.227D1,1.849D3,
     +        1.116D0,5.513D1,1.293D0,2.028D4,1.105D1,1.936D1,1.369D2/
      DATA E0/12.78,15.53,24.36,4.255,0.6882,9.515,0.808,1.366,0.0552,
     +        16.05,0.2288,23.45,10.08,9.98,130.9,4.293,26.18,94.72,
     +        629.7,172.9/
      DATA EMX/34.43,40.9,373.1,394.4,417.5,442.3,468.7,496.7,527.,
     +         556.9,4265.,4362.,4453.,4555.,4659.,4767.,4880.,4982.,
     +         5.D4,5.D4/
      DATA Y0/1.012D-3,2.161D-3,1.802,14.67,121.,4.829,148.7,103.9,
     +        2.826D-4,0.,24.78,24.17,6.138D-3,71.04,1.833D-2,0.9363,
     +        2.402D-2,0.,0.,0./
      DATA Y1/1.851D-2,6.706D-2,1.233,3.298D-2,3.876,5.824,1.283,3.329,
     +        1.657,0.,3.1,0.5469,69.31,5.311,0.9359,4.589D-2,9.323D-3,
     +        0.,0.,0./
      DATA YW/0.4477,0.6453,0.3126,1.369,8.277,2.471,0.572,0.2806,
     +        1.843D-3,0.,1.39,6.842D-4,241.,3.879,9.084D-2,3.461D-5,
     +        28.03,0.,0.,0./
      DATA YA/0.3162,0.779,293.1,13.55,150.2,89.73,368.2,318.8,1.79D4,
     +        698.9,254.9,13.12,1.792D4,59.18,382.8,16.91,1.456,38.18,
     +        39.21,32.88/
      DATA PV/12.42,21.3,3.944,12.36,10.61,5.141,8.634,8.138,5.893,
     +        3.857,11.03,9.771,2.868,9.005,2.023,14.38,25.6,4.192,
     +        1.862,2.963/
      DATA S95/9.017D1,7.314D1,1.945D2,1.542D2,1.622D2,1.855D2,2.181D2,
     +         2.788D2,1.934D2,6.616D2,1.547D1,1.324D1,1.57D1,1.384D1,
     +         1.417D1,1.665D1,1.486D1,1.82D1,1.936D1,1.369D2/
      DATA E95/44.87,44.98,126.,141.3,138.4,130.3,120.8,107.,129.3,
     +         65.11,701.,750.3,698.9,739.6,734.2,686.2,723.5,664.,
     +         629.7,172.9/
      DATA Y95/14.65,18.98,68.19,99.06,88.11,69.93,58.16,47.68,70.,
     +         43.71,31.97,50.,32.18,50.,50.,34.43,50.,39.79,39.21,
     +         32.88/
      DATA YW95/0.2754,0.2735,4.791D-4,1.107D-3,4.384D-4,1.4D-5,
     +          4.346D-6,4.591D-6,0.1,7.881D-6,10*0./
      DATA P95/7.498,7.152,3.77,3.446,3.521,3.707,3.907,4.2,3.7,4.937,
     +         1.858,1.65,1.851,1.65,1.65,1.823,1.65,1.777,1.862,
     +         2.963/
C
      VERN20=0.
      IVER=IZZ
C
C     1996 Expression
C
      IF(E.LT.EMX(IVER)) THEN
        XX=E/E0(IVER)-Y0(IVER)
        YY=SQRT(XX*XX+Y1(IVER)*Y1(IVER))
        AA=(XX-UN)*(XX-UN)+YW(IVER)*YW(IVER)
        BB=YY**(HALF*PV(IVER)-5.5)
        CC=(UN+SQRT(YY/YA(IVER)))**PV(IVER)
        FY=AA*BB/CC
        VERN20=S0(IVER)*T18*FY
      ELSE
C
C     1995 Expression for high energies
C          (ionization of inner shell electron)
C
        YY=E/E95(IVER)
        XL=0.
        IF(IZZ.LE.10) XL=UN
        Q=HALF*P95(IVER)-5.5-XL
        AA=(YY-UN)*(YY-UN)+YW95(IVER)*YW95(IVER)
        BB=YY**Q
        CC=(UN+SQRT(YY/Y95(IVER)))**P95(IVER)
        FY=AA*BB/CC
        VERN20=S95(IVER)*T18*FY
      END IF
C
      RETURN
      END
C
C
C     ****************************************************************
C
C
 
      FUNCTION GAUNT(I,FR)
C     ====================
C
C     Hydrogenic bound-free Gaunt factor for the principal quantum
C     number I and frequency FR
C
      INCLUDE 'IMPLIC.FOR'
      PARAMETER (UN=1.)
      DIMENSION CGT(7,10),X(7),FRKW(10)
      DATA CGT/0.,12.803223,-5.5759888,1.2302628,-2.9094219D-3,
     *  7.3993579D-6,-8.7356966D-9,-2.0244141,2.1325684,-1.2709045,
     *  1.1595421,-2.0735860D-3,2.7033384D-6,0.,-0.23387146,
     *  0.52471924,-0.55936432,1.1450949,-1.9366592D-3,2.3572356D-6,
     *  0.,-5.4418565D-2,0.19683564,-0.31190730,1.1306695,
     *  -1.3482273D-3,-4.6949424D-6,2.3548636D-8,-8.9182854D-3,
     *  5.5545091D-2,-0.16051018,1.1190904,-1.0401085D-3,
     *  -6.9943488D-6,2.8496742D-8,-5.5303574D-3,4.1921183D-2,
     *  -0.13075417,1.1168376,-8.9466573D-4,-8.8393133D-6,
     *  3.4696768D-8,-2.2752881D-3,2.3350812D-2,-9.5441161D-2,
     *  1.1128632,-7.4833260D-4,-1.0244504D-5,3.8595771D-8,
     *  -9.7200274D-4,1.3298411D-2,-7.1010560D-2,1.1093137,
     *  -6.2619148D-4,-1.1342068D-5,4.1477731D-8,-4.9576163D-4,
     *  8.5139736D-3,-5.6046560D-2,1.1078717,-5.4837392D-4,
     *  -1.2157943D-5,4.3796716D-8,-2.9467046D-4,6.1516856D-3,
     *  -4.7326370D-2,1.1052734,-4.4341570D-4,-1.3235905D-5,
     *  4.7003140D-8/
      DATA FRKW/6.6D15,9*3.3D15/

      IF(I.LE.10) THEN
C      IF(FR.LT.FRKW(I)) THEN
          X(5)=FR/2.99793D14
          X(6)=X(5)*X(5)
          X(7)=X(6)*X(5)
          X(4)=UN
          X(3)=UN/X(5)
          X(2)=X(3)*X(3)
          X(1)=X(2)*X(3)
          GAUNT=0.
          DO 10 J=1,7
            GAUNT=GAUNT+CGT(J,I)*X(J)
   10     CONTINUE
c        ELSE
C        GAUNT=GNTK(I,FR)
C      END IF
      ELSE
        GAUNT=UN
      ENDIF
      RETURN
      END
C
C
C     ****************************************************************
C
C
 
      FUNCTION GNTK(I,FR)
C     ===================
C
C     Hydrogenic bound-free Gaunt factor for the principal quantum
C     number I and frequency FR (from Klaus Werner)
C
      INCLUDE 'IMPLIC.FOR'
      GNTK=1.
      IF(I.GT.3) GO TO 16
      Y=1./FR      
      GO TO (1,2,3),I
    1 GNTK=0.9916+Y*(2.71852D13-Y*2.26846D30)
      GO TO 16
    2 GNTK=1.1050-Y*(2.37490D14-Y*4.07677D28)
      GO TO 16
    3 GNTK=1.1010-Y*(0.98632D14-Y*1.03540D28)
   16 RETURN
      END
C
C
C     ****************************************************************
C
C
      SUBROUTINE SPSIGK(IB,FR,SIGSP)
C     ==============================
C
C     Non-standard evaluation of the photoionization cross-sections
C     Basically user-suppled procedure; here are some examples
C
      INCLUDE 'IMPLIC.FOR'
      SIGSP=0.
C
C     Special formula for the He I ground state
C
      IF(IB.EQ.-201) SIGSP=7.3D-18*EXP(1.373-2.311E-16*FR)
C
C     Special formula for the averaged <n=2> level of He I
C
      IF(IB.EQ.-202) SIGSP=SGHE12(FR)
C
C     Carbon ground configuration levels 2p2 1D and 1S
C
      IF(IB.EQ.-602.OR.IB.EQ.-603) THEN
         CALL CARBON(IB,FR,SG)
         SIGSP=SG
      END IF
C
C     Hidalgo (Ap.J. 153, 981, 1968) photoionization data
C
      IF(IB.LE.-101.AND.IB.GE.-137) SIGSP=HIDALG(IB,FR)
C
C     Reilman and Manson (Ap.J. Suppl. 40, 815, 1979) photoionization data
C
      IF(IB.LE.-301.AND.IB.GE.-337) SIGSP=REIMAN(IB,FR)
      RETURN
      END
C
C
C     *********************************************************************
C
C

      FUNCTION HIDALG(IB,FR)
C     ======================
C
C     Read table of wavelengths and photo-ionization cross-sections
C     from Hidalgo (1968, Ap. J., 153, 981) for the species indicated by IB
C     (Hidalgo's number = INDEX = -IB-100).
C     Compute linearly interpolated value of the cross-section 
C     at the frequency FR.
C
C     (At the moment, only a few transitions are considered)
C
      INCLUDE 'IMPLIC.FOR'
      DIMENSION WL1(20),WL2(20),WL(20),SIG0(20,24),SIGS(20)
C
      DATA WL1 /
     *  39.1, 80.9, 97.6,100.1,104.3,107.2,108.7,111.9,113.6,115.4,
     * 117.1,119.0,124.8,126.9,129.1,131.3,133.6,136.0,138.5,141.1/
      DATA WL2 /
     *  68.5, 80.9,100.1,120.9,158.8,165.7,177.3,190.6,200.7,206.2,
     * 211.9,218.0,224.5,231.3,246.3,5*0./
      DATA SIG0 /
     * 120*0.,
     *.0460,.2400,.3500,.3700,.4000,.4300,.4400,.4600,.4700,.4900, 
     *.5000,.5200,.5700,.6200, 6*0.,
     * 80*0.,
     *.0092,.1000,.1900,.2100,.2300,.2500,.2600,.2900,.3000,.3200,
     *.3400,.3500,.4100,.4300,.4500,.4800,.5000,.5300,.5600,.5900,
     * 20*0.,
     *.3400,.4600,.6300,.7700,.9100,1.080, 14*0.,
     * 20*0.,
     *.0064,.1100,.2200,.4100,.9400,1.000,1.300,1.600, 12*0.,
     * 80*0.,
     *.0370,.0650,.1300,.2400,.5500,.6300,.7700,.9500,1.100,1.250,
     * 10*0.,
     * 20*0.,
     * 20*0.,
     *.0220,.0390,.0800,.1500,.3500,.4000,.4900,.6200,.7200,.7800,
     *.8500,.9300,1.020,
     * 7*0./
      SAVE WL1,WL2,SIG0
C

      INDEX=-IB-100
      NUM=20  
      IF(INDEX.GE.13.AND.INDEX.LE.27) NUM=15
      DO 10 I=1,NUM
         IF(INDEX.LT.13) WL(I)=WL1(I)
         IF(INDEX.GE.13) WL(I)=WL2(I)
         SIGS(I)=SIG0(I,INDEX)
   10 CONTINUE
C 
      WLAM=2.997925D18/FR
      IL=1
      IR=NUM
      DO 50 I=1,NUM-1
        IF(WLAM.GE.WL(I).AND.WLAM.LE.WL(I+1)) THEN
          IL=I
          IR=I+1
          GO TO 60
        ENDIF
 50   CONTINUE
C
C     LINEAR INTERPOLATION:
C
 60   SIGM=(SIGS(IR)-SIGS(IL))*(WLAM-WL(IL))/(WL(IR)-WL(IL))
     *      + SIGS(IL)
C
C     IF OUTSIDE WAVELENGTH RANGE SET TO FIRST(LAST) VALUE:
C
       IF(WLAM.LE.WL(1)) SIGM=SIGS(1)
       IF(WLAM.GE.WL(NUM)) SIGM=SIGS(NUM)
C
C     IF LAST NON-ZERO SIG VALUES, NO INTERPOLATION:
C
c       IF(SIGS(IR).EQ.0.) SIGM=SIGS(IL)   
C
      HIDALG=SIGM*1.D-18
      RETURN
      END
C
C
C     *********************************************************************
C
C

      FUNCTION REIMAN(IB,FR)
C     ======================
C
C     Read table of photon energies and photo-ionization cross-sections
C     from Reilman & Manson (1979, Ap. J. Suppl., 40, 815) for the species 
C     indicated by IB
C    
C     Compute linearly interpolated value of the cross-section 
C     at the frequency FR.
C
C     (At the moment, only a few transitions are considered)
C
      INCLUDE 'IMPLIC.FOR'
      DIMENSION HEV(30),F0(30),SIG0(30,2),SIGS(30)
C
      DATA HEV /
     * 130.,160.,190.,210.,240.,270.,300.,330.,360.,390.,
     * 420.,450.,480.,510.,540.,570.,600.,630.,660.,690.,
     * 720.,750.,780.,810.,840.,870.,900.,930.,960.,990./
      DATA SIG0 /
     * 3*0.,     4.422D-1, 3.478D-1, 
     * 2.794D-1, 2.286D-1, 1.899D-1, 1.598D-1, 1.360D-1,
     * 1.169D-1, 1.013D-1, 8.845D-2, 7.776D-2, 6.877D-2,
     * 6.114D-2, 5.463D-2, 4.904D-2, 4.419D-2, 3.998D-2,
     * 3.629D-2, 3.305D-2, 3.019D-2, 2.766D-2, 2.540D-2,
     * 2.339D-2, 2.158D-2, 1.996D-2, 1.850D-2, 1.718D-2,
     * 4*0.,     1.981D-1, 1.584D-1, 
     * 1.290D-1, 1.066D-1, 8.932D-2, 7.567D-2, 6.475D-2,
     * 5.589D-2, 4.862D-2, 4.259D-2, 3.754D-2, 3.329D-2,
     * 2.966D-2, 2.656D-2, 2.388D-2, 2.157D-2, 1.954D-2,
     * 1.777D-2, 1.621D-2, 1.484D-2, 1.362D-2, 1.253D-2,
     * 1.155D-2, 1.067D-2, 9.888D-3, 9.179D-3/
      SAVE HEV,SIG0
C
      INDEX=-IB-300
      NUM=30  
      DO 10 I=1,NUM
         F0(I)=HEV(I)*2.418573D14
         SIGS(I)=SIG0(I,INDEX)
   10 CONTINUE
C 
      IL=1
      IR=NUM
      DO 50 I=1,NUM-1
        IF(FR.GE.F0(I).AND.FR.LE.F0(I+1)) THEN
          IL=I
          IR=I+1
          GO TO 60
        ENDIF
 50   CONTINUE
C
C     LINEAR INTERPOLATION:
C
 60   SIGM=(SIGS(IR)-SIGS(IL))*(FR-F0(IL))/(F0(IR)-F0(IL))
     *      + SIGS(IL)
C
C     IF OUTSIDE WAVELENGTH RANGE SET TO FIRST(LAST) VALUE:
C
       IF(FR.LE.F0(1)) SIGM=SIGS(1)
       IF(FR.GE.F0(NUM)) SIGM=SIGS(NUM)
C
C     IF LAST NON-ZERO SIG VALUES, NO INTERPOLATION:
C
c       IF(SIGS(IR).EQ.0.) SIGM=SIGS(IL)   
C
      REIMAN=SIGM*1.D-18
      RETURN
      END
C
C
C     ****************************************************************
C
C
 
      SUBROUTINE CARBON(IB,FR,SG)
C     ===========================
C
C     Photoionization cross-section for neutral carbon 2p1D and 2p1S
C     levels (G.B.Taylor - private communication)
C
      INCLUDE 'IMPLIC.FOR'
      DIMENSION FR2(34),SG2(34),FR3(45),SG3(45)
      PARAMETER (FR0=3.28805D15, NC2=34, NC3=45)
      DATA FR2/ 0.74, 0.75, 0.76, 0.77, 0.78, 0.79, 0.80, 0.81, 0.82,
     *                0.83,       0.85, 0.86, 0.87, 0.88, 0.89, 0.90,
     *          0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99,
     *          1.00, 1.10, 1.20, 1.30, 1.45, 1.50, 1.60, 1.80, 2./
      DATA SG2/ 12.04, 12.03, 12.09, 12.26, 12.60, 13.24, 14.36, 16.24,
     *          19.28, 23.94, 37.41, 42.88, 44.76, 43.41, 40.46, 37.19,
     *          34.26, 31.82, 29.96, 28.57, 27.68, 27.37, 27.84, 29.69,
     *          34.45, 46.35, 13.80, 11.54, 10.40,  8.96,  8.54,  7.47,
     *           6.53,  5.66/
      DATA FR3/ 0.66, 0.68, 0.70, 0.72, 0.74, 0.76, 0.78, 0.80, 0.82,
     *          0.84, 0.86, 0.864,0.866,0.868,0.87, 0.874,0.876,0.88,
     *          0.882,0.884,0.886,0.888,0.89 ,0.894,0.896,0.898,0.90,
     *          0.904,0.908,0.910,0.920,0.94, 0.98, 1.00, 1.10, 1.20,
     *          1.26, 1.34, 1.36, 1.40, 1.46, 1.60, 1.70, 1.80, 2./
      DATA SG3/ 13.94, 13.29, 12.56, 11.73, 10.82, 10.18,  8.62,  7.27,
     *           5.74,  4.14,  4.61,  5.92,  6.94,  8.34, 10.21, 16.12,
     *          20.64, 34.56, 44.82, 57.71, 73.09, 89.99,106.38,127.08,
     *         128.38,124.44,117.17, 99.32, 82.95, 76.05, 52.65, 33.23,
     *          21.29, 18.69, 12.62, 11.44,  9.77,  7.53, 10.47,  9.65,
     *          10.19,  7.28,  6.70,  6.11,  4.96/
      SAVE FR2,SG2,FR3,SG3
C
      F=FR/FR0
      IF(IB.NE.-602) GO TO 25
      J=2
      IF(F.LE.FR2(1)) GO TO 20
      DO 10 I=2,NC2
         J=I
         IF(F.GT.FR2(I-1).AND.F.LE.FR2(I)) GO TO 20
   10 CONTINUE
   20 SG=(F-FR2(J-1))/(FR2(J)-FR2(J-1))*(SG2(J)-SG2(J-1))+SG2(J-1)
      SG=SG*1.D-18
   25 IF(IB.NE.-603) GO TO 50
      J=2
      IF(F.LE.FR3(1)) GO TO 40
      DO 30 I=2,NC3
         J=I
         IF(F.GT.FR3(I-1).AND.F.LE.FR3(I)) GO TO 40
   30 CONTINUE
   40 SG=(F-FR3(J-1))/(FR3(J)-FR3(J-1))*(SG3(J)-SG3(J-1))+SG3(J-1)
      SG=SG*1.D-18
   50 CONTINUE
      RETURN
      END
C
C
C     ****************************************************************
C

      FUNCTION SGHE12(FR)
C     ===================
C
C     Special formula for the photoionization cross-section from the
C     averaged <n=2> level of He I
C
      INCLUDE 'IMPLIC.FOR'
      PARAMETER (C1=3.D0,C2=9.D0,C3=1.6D1,T15=1.D-15,
     * A1=6.45105D-18,A2=3.02D-19,A3=9.9847D-18,A4=1.1763673D-17,
     * A5=3.63662D-19,A6=-2.783D2,A7=1.488D1,A8=-2.311D-1,
     * E1=3.5D0,E2=3.6D0,E3=1.91D0,E4=2.9D0,E5=3.3D0)
C
      X=FR*T15
      XX=LOG(FR)
      SGHE12=(C1*(A1/X**E1+A2/X**E2)+A3/X**E3+C2*(A4/X**E4+A5/X**E5)+
     *       C1*EXP(A6+XX*(A7+XX*A8)))/C3
      RETURN
      END
C
C
C     ****************************************************************
C
C

      FUNCTION SBFHE1(II,IB,FR,GG)
C     ============================
C
C     Calculates photoionization cross sections of neutral helium
C     from states with n = 1, 2, 3, 4.
C
C     The levels are either non-averaged (l,s) states, or some
C     averaged levels.
C     The program allows only two standard possibilities of
C     constructing averaged levels:
C     i)  all states within given principal quantum number n (>1) are
C         lumped together
C     ii) all siglet states for given n, and all triplet states for
C         given n are lumped together separately (there are thus two
C         explicit levels for a given n)
C
C     The cross sections are calculated using appropriate averages
C     of the Opacity Project cross sections, calculated by procedure
C     HEPHOT
C
C     For IB=21 or IB=23, Koester (1985, AA 149, 423) fits
C
C     Input parameters:
C      II    - index of the lower level (in the numbering of explicit
C              levels)
C      IB    - photoionization switch IBF for the given transition
C            = 10        -  means that the given transition is from an
C                           averaged level
C            = 11 or 21  -  the given transition is from non-averaged
C                           singlet state
C            = 13 or 23  -  the given transition is from non-averaged
C                           triplet state
C      FR    - frequency
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'ATOMIC.FOR'
C
      NI=NQUANT(II)
      IGI=INT(G(II)+0.01)
      IS=IB-10
      IF(IB.GT.20) IS=IB-20
      sbfhe1=0.
C
C     ----------------------------------------------------------------
C     IB=11 or 13  - photoionization from an non-averaged (l,s) level
C     ----------------------------------------------------------------
C
      IF(IS.EQ.1.OR.IS.EQ.3) THEN
         IL=(IGI/IS-1)/2
       IF(IB.LT.20) THEN
           SBFHE1=HEPHOT(IS,IL,NI,FR)
       ELSE
         SBFHE1=CKOEST(IS,IL,NI,FR,GG)
       ENDIF
      END IF
C
C     ----------------------------------------------------------------
C     IS=10 - photoionization from an averaged level
C     ----------------------------------------------------------------
C
      IF(IS.EQ.0) THEN
         IF(NI.EQ.2) THEN
C
C ********    photoionization from an averaged level with n=2
C
            IF(IGI.EQ.4) THEN
C
C      a) lower level is an averaged singlet state
C
               SBFHE1=(HEPHOT(1,0,2,FR)+3.D0*HEPHOT(1,1,2,FR))/9.D0
            ELSE IF(IGI.EQ.12) THEN
C
C      b) lower level is an averaged triplet state
C
               SBFHE1=(HEPHOT(3,0,2,FR)+3.D0*HEPHOT(3,1,2,FR))/9.D0
            ELSE IF(IGI.EQ.16) THEN
C
C      c) lower level is an average of both singlet and triplet states
C
               SBFHE1=(HEPHOT(1,0,2,FR)+3.D0*(HEPHOT(1,1,2,FR)+
     *                HEPHOT(3,0,2,FR))+9.D0*HEPHOT(3,1,2,FR))/1.6D1
            ELSE
               GO TO 10
            END IF
C
C
C ********    photoionization from an averaged level with n=3
C
         ELSE IF(NI.EQ.3) THEN
            IF(IGI.EQ.9) THEN
C
C      a) lower level is an averaged singlet state
C
               SBFHE1=(HEPHOT(1,0,3,FR)+3.D0*HEPHOT(1,1,3,FR)+
     *                5.D0*HEPHOT(1,2,3,FR))/9.D0
            ELSE IF(IGI.EQ.27) THEN
C
C      b) lower level is an averaged triplet state
C
               SBFHE1=(HEPHOT(3,0,3,FR)+3.D0*HEPHOT(3,1,3,FR)+
     *                5.D0*HEPHOT(3,2,3,FR))/9.D0
            ELSE IF(IGI.EQ.36) THEN
C
C      c) lower level is an average of both singlet and triplet states
C
               SBFHE1=(HEPHOT(1,0,3,FR)+3.D0*HEPHOT(1,1,3,FR)+
     *                5.D0*HEPHOT(1,2,3,FR)+
     *                3.D0*HEPHOT(3,0,3,FR)+9.D0*HEPHOT(3,1,3,FR)+
     *                15.D0*HEPHOT(3,2,3,FR))/3.6D0
            ELSE
               GO TO 10
            END IF
         ELSE IF(NI.EQ.4) THEN
C
C ********    photoionization from an averaged level with n=4
C
            IF(IGI.EQ.16) THEN
C
C      a) lower level is an averaged singlet state
C
               SBFHE1=(HEPHOT(1,0,4,FR)+3.D0*HEPHOT(1,1,4,FR)+
     *                 5.D0*HEPHOT(1,2,4,FR)+
     *                 7.D0*HEPHOT(1,3,4,FR))/1.6D1
            ELSE IF(IGI.EQ.48) THEN
C
C      b) lower level is an averaged triplet state
C
               SBFHE1=(HEPHOT(3,0,4,FR)+3.D0*HEPHOT(3,1,4,FR)+
     *                 5.D0*HEPHOT(3,2,4,FR)+
     *                 7.D0*HEPHOT(3,3,4,FR))/1.6D1
            ELSE IF(IGI.EQ.64) THEN
C
C      c) lower level is an average of both singlet and triplet states
C
               SBFHE1=(HEPHOT(1,0,4,FR)+3.D0*HEPHOT(1,1,4,FR)+
     *                 5.D0*HEPHOT(1,2,4,FR)+
     *                 7.D0*HEPHOT(1,3,4,FR)+
     *                 3.D0*HEPHOT(3,0,4,FR)+
     *                 9.D0*HEPHOT(3,1,4,FR)+
     *                 15.D0*HEPHOT(3,2,4,FR)+
     *                 21.D0*HEPHOT(3,3,4,FR))/6.4D1
            ELSE
               GO TO 10
            END IF
         ELSE
            GO TO 10
         END IF
      END IF
      RETURN
   10 WRITE(6,601) NI,IGI,IS
      WRITE(10,601) NI,IGI,IS
  601 FORMAT(1H0/' INCONSISTENT INPUT TO PROCEDURE SBFHE1'/
     * ' QUANTUM NUMBER =',I3,'  STATISTICAL WEIGHT',I4,'  S=',I3)
      call quit(' ',ni,igi)
      
      END
C
C
C     ****************************************************************
C
C

      FUNCTION HEPHOT(S,L,N,FREQ)
C     ===========================
C
C
C
C
C           EVALUATES HE I PHOTOIONIZATION CROSS SECTION USING SEATON
C           FERNLEY'S CUBIC FITS TO THE OPACITY PROJECT CROSS SECTIONS
C           UP TO SOME ENERGY "EFITM" IN THE RESONANCE-FREE ZONE.  BEYOND
C           THIS ENERGY LINEAR FITS TO LOG SIGMA IN LOG (E/E0) ARE USED.
C           THIS EXTRAPOLATION SHOULD BE USED UP TO THE BEGINNING OF THE
C           RESONANCE ZONE "XMAX", BUT AT PRESENT IT IS USED THROUGH IT.
C           BY CHANGING A FEW LINES THAT ARE PRESENTLY COMMENTED OUT,
C           FOR ENERGIES IN THE RESONANCE ZONE A VALUE OF 1/100 OF THE
C           THRESHOLD CROSS SECTION IS USED -- THIS IS PURELY AD HOC AND
C           ONLY A TEMPORARY MEASURE.  OBVIOUSLY ANY OTHER VALUE OR FUNCTIONAL
C           FORM CAN BE INSERTED HERE.
C
C           CALLING SEQUENCE INCLUDES:
C                S = MULTIPLICITY, EITHER 1 OR 3
C                L = ANGULAR MOMENTUM, 0, 1, OR 2;
C                    for L > 2 - hydrogenic expresion
C                FREQ = FREQUENCY
C
C           DGH JUNE 1988 JILA, slightly modified by I.H.
C
      INCLUDE 'IMPLIC.FOR'
      INTEGER S,L,SS,LL
      PARAMETER (TENM18=1.D-18, FRH=3.28805D15, TWO=2.D0,
     *           TENLG=2.302585093, PHOT0=2.815D29)
      DIMENSION COEF(4,53),IST(3,2),N0(3,2),
     *          FL0(53),A(53),B(53),XFITM(53)
C
      DATA IST/1,36,20,11,45,28/
      DATA N0/1,2,3,2,2,3/
C
      DATA FL0/
     . 2.521D-01,-5.381D-01,-9.139D-01,-1.175D+00,-1.375D+00,-1.537D+00,
     .-1.674D+00,-1.792D+00,-1.896D+00,-1.989D+00,-4.555D-01,-8.622D-01,
     .-1.137D+00,-1.345D+00,-1.512D+00,-1.653D+00,-1.774D+00,-1.880D+00,
     .-1.974D+00,-9.538D-01,-1.204D+00,-1.398D+00,-1.556D+00,-1.690D+00,
     .-1.806D+00,-1.909D+00,-2.000D+00,-9.537D-01,-1.204D+00,-1.398D+00,
     .-1.556D+00,-1.690D+00,-1.806D+00,-1.909D+00,-2.000D+00,-6.065D-01,
     .-9.578D-01,-1.207D+00,-1.400D+00,-1.558D+00,-1.692D+00,-1.808D+00,
     .-1.910D+00,-2.002D+00,-5.749D-01,-9.352D-01,-1.190D+00,-1.386D+00,
     .-1.547D+00,-1.682D+00,-1.799D+00,-1.902D+00,-1.995D+00/
C
      DATA XFITM/
     . 3.262D-01, 6.135D-01, 9.233D-01, 8.438D-01, 1.020D+00, 1.169D+00,
     . 1.298D+00, 1.411D+00, 1.512D+00, 1.602D+00, 7.228D-01, 1.076D+00,
     . 1.206D+00, 1.404D+00, 1.481D+00, 1.464D+00, 1.581D+00, 1.685D+00,
     . 1.777D+00, 9.586D-01, 1.187D+00, 1.371D+00, 1.524D+00, 1.740D+00,
     . 1.854D+00, 1.955D+00, 2.046D+00, 9.585D-01, 1.041D+00, 1.371D+00,
     . 1.608D+00, 1.739D+00, 1.768D+00, 1.869D+00, 1.803D+00, 7.360D-01,
     . 1.041D+00, 1.272D+00, 1.457D+00, 1.611D+00, 1.741D+00, 1.855D+00,
     . 1.870D+00, 1.804D+00, 9.302D-01, 1.144D+00, 1.028D+00, 1.210D+00,
     . 1.362D+00, 1.646D+00, 1.761D+00, 1.863D+00, 1.954D+00/
C
      DATA A/
     . 6.95319D-01, 1.13101D+00, 1.36313D+00, 1.51684D+00, 1.64767D+00,
     . 1.75643D+00, 1.84458D+00, 1.87243D+00, 1.85628D+00, 1.90889D+00,
     . 9.01802D-01, 1.25389D+00, 1.39033D+00, 1.55226D+00, 1.60658D+00,
     . 1.65930D+00, 1.68855D+00, 1.62477D+00, 1.66726D+00, 1.83599D+00,
     . 2.50403D+00, 3.08564D+00, 3.56545D+00, 4.25922D+00, 4.61346D+00,
     . 4.91417D+00, 5.19211D+00, 1.74181D+00, 2.25756D+00, 2.95625D+00,
     . 3.65899D+00, 4.04397D+00, 4.13410D+00, 4.43538D+00, 4.19583D+00,
     . 1.79027D+00, 2.23543D+00, 2.63942D+00, 3.02461D+00, 3.35018D+00,
     . 3.62067D+00, 3.85218D+00, 3.76689D+00, 3.49318D+00, 1.16294D+00,
     . 1.86467D+00, 2.02110D+00, 2.24231D+00, 2.44240D+00, 2.76594D+00,
     . 2.93230D+00, 3.08109D+00, 3.21069D+00/
C
      DATA B/
     .-1.29000D+00,-2.15771D+00,-2.13263D+00,-2.10272D+00,-2.10861D+00,
     .-2.11507D+00,-2.11710D+00,-2.08531D+00,-2.03296D+00,-2.03441D+00,
     .-1.85905D+00,-2.04057D+00,-2.02189D+00,-2.05930D+00,-2.03403D+00,
     .-2.02071D+00,-1.99956D+00,-1.92851D+00,-1.92905D+00,-4.58608D+00,
     .-4.40022D+00,-4.39154D+00,-4.39676D+00,-4.57631D+00,-4.57120D+00,
     .-4.56188D+00,-4.55915D+00,-4.41218D+00,-4.12940D+00,-4.24401D+00,
     .-4.40783D+00,-4.39930D+00,-4.25981D+00,-4.26804D+00,-4.00419D+00,
     .-4.47251D+00,-3.87960D+00,-3.71668D+00,-3.68461D+00,-3.67173D+00,
     .-3.65991D+00,-3.64968D+00,-3.48666D+00,-3.23985D+00,-2.95758D+00,
     .-3.07110D+00,-2.87157D+00,-2.83137D+00,-2.82132D+00,-2.91084D+00,
     .-2.91159D+00,-2.91336D+00,-2.91296D+00/
C
      DATA ((COEF(I,J),I=1,4),J=1,10)/
     . 8.734D-01,-1.545D+00,-1.093D+00, 5.918D-01, 9.771D-01,-1.567D+00,
     .-4.739D-01,-1.302D-01, 1.174D+00,-1.638D+00,-2.831D-01,-3.281D-02,
     . 1.324D+00,-1.692D+00,-2.916D-01, 9.027D-02, 1.445D+00,-1.761D+00,
     .-1.902D-01, 4.401D-02, 1.546D+00,-1.817D+00,-1.278D-01, 2.293D-02,
     . 1.635D+00,-1.864D+00,-8.252D-02, 9.854D-03, 1.712D+00,-1.903D+00,
     .-5.206D-02, 2.892D-03, 1.782D+00,-1.936D+00,-2.952D-02,-1.405D-03,
     . 1.845D+00,-1.964D+00,-1.152D-02,-4.487D-03/
      DATA ((COEF(I,J),I=1,4),J=11,19)/
     . 7.377D-01,-9.327D-01,-1.466D+00, 6.891D-01, 9.031D-01,-1.157D+00,
     .-7.151D-01, 1.832D-01, 1.031D+00,-1.313D+00,-4.517D-01, 9.207D-02,
     . 1.135D+00,-1.441D+00,-2.724D-01, 3.105D-02, 1.225D+00,-1.536D+00,
     .-1.725D-01, 7.191D-03, 1.302D+00,-1.602D+00,-1.300D-01, 7.345D-03,
     . 1.372D+00,-1.664D+00,-8.204D-02,-1.643D-03, 1.434D+00,-1.715D+00,
     .-4.646D-02,-7.456D-03, 1.491D+00,-1.760D+00,-1.838D-02,-1.152D-02/
      DATA ((COEF(I,J),I=1,4),J=20,27)/
     . 1.258D+00,-3.442D+00,-4.731D-01,-9.522D-02, 1.553D+00,-2.781D+00,
     .-6.841D-01,-4.083D-03, 1.727D+00,-2.494D+00,-5.785D-01,-6.015D-02,
     . 1.853D+00,-2.347D+00,-4.611D-01,-9.615D-02, 1.955D+00,-2.273D+00,
     .-3.457D-01,-1.245D-01, 2.041D+00,-2.226D+00,-2.669D-01,-1.344D-01,
     . 2.115D+00,-2.200D+00,-1.999D-01,-1.410D-01, 2.182D+00,-2.188D+00,
     .-1.405D-01,-1.460D-01/
      DATA ((COEF(I,J),I=1,4),J=28,35)/
     . 1.267D+00,-3.417D+00,-5.038D-01,-1.797D-02, 1.565D+00,-2.781D+00,
     .-6.497D-01,-5.979D-03, 1.741D+00,-2.479D+00,-6.099D-01,-2.227D-02,
     . 1.870D+00,-2.336D+00,-4.899D-01,-6.616D-02, 1.973D+00,-2.253D+00,
     .-3.972D-01,-8.729D-02, 2.061D+00,-2.212D+00,-3.072D-01,-1.060D-01,
     . 2.137D+00,-2.189D+00,-2.352D-01,-1.171D-01, 2.205D+00,-2.186D+00,
     .-1.621D-01,-1.296D-01/
      DATA ((COEF(I,J),I=1,4),J=36,44)/
     . 1.129D+00,-3.149D+00,-1.910D-01,-5.244D-01, 1.431D+00,-2.511D+00,
     .-3.710D-01,-1.933D-01, 1.620D+00,-2.303D+00,-3.045D-01,-1.391D-01,
     . 1.763D+00,-2.235D+00,-1.829D-01,-1.491D-01, 1.879D+00,-2.215D+00,
     .-9.003D-02,-1.537D-01, 1.978D+00,-2.213D+00,-2.066D-02,-1.541D-01,
     . 2.064D+00,-2.220D+00, 3.258D-02,-1.527D-01, 2.140D+00,-2.225D+00,
     . 6.311D-02,-1.455D-01, 2.208D+00,-2.229D+00, 7.977D-02,-1.357D-01/
      DATA ((COEF(I,J),I=1,4),J=45,53)/
     . 1.204D+00,-2.809D+00,-3.094D-01, 1.100D-01, 1.455D+00,-2.254D+00,
     .-4.795D-01, 6.872D-02, 1.619D+00,-2.109D+00,-3.357D-01,-2.532D-02,
     . 1.747D+00,-2.065D+00,-2.317D-01,-5.224D-02, 1.853D+00,-2.058D+00,
     .-1.517D-01,-6.647D-02, 1.943D+00,-2.055D+00,-1.158D-01,-6.081D-02,
     . 2.023D+00,-2.070D+00,-6.470D-02,-6.800D-02, 2.095D+00,-2.088D+00,
     .-2.357D-02,-7.250D-02, 2.160D+00,-2.107D+00, 1.065D-02,-7.542D-02/
C
      SAVE COEF,IST,N0,FL0,A,B,XFITM
C
      IF(L.GT.2) GO TO 20
C
C          SELECT BEGINNING AND END OF COEFFICIENTS
C
      SS=(S+1)/2
      LL=L+1
      NSL0=N0(LL,SS)
      I=IST(LL,SS)+N-NSL0
C
C          EVALUATE CROSS SECTION
C
      FL=LOG10(FREQ/FRH)
      X=FL-FL0(I)
      IF(X.GE.-0.001D0) THEN
         IF(X.LT.XFITM(I)) THEN
            P=COEF(4,I)
            DO 10 K=1,3
               P=X*P+COEF(4-K,I)
   10       CONTINUE
            HEPHOT=TENM18*EXP(TENLG*P)
          ELSE
            HEPHOT=TENM18*EXP(TENLG*(A(I)+B(I)*X))
          END IF
       ELSE
          HEPHOT=0.
      END IF
      RETURN
C
C     Hydrogenic expression for L > 2
C      [multiplied by relative population of state (s,l,n), ie.
C       by  stat.weight(s,l)/stat.weight(n)]
C
   20 GN=TWO*N*N
      HEPHOT=PHOT0/FREQ/FREQ/FREQ/N**5*(2*L+1)*S/GN
      RETURN
      END
C
C
C     ****************************************************************
C
C

      FUNCTION CKOEST(S,L,N,FREQ,GG)
C     ==============================
C
C           EVALUATES HE I PHOTOIONIZATION CROSS SECTION USING
C           KOESTER'S FITS (1985, AA 149, 423)
C
C           CALLING SEQUENCE INCLUDES:
C                S = MULTIPLICITY, EITHER 1 OR 3
C                L = ANGULAR MOMENTUM, 0, 1, OR 2;
C                    for L > 2 - hydrogenic expresion
C                N = PRINCIPAL QUANTUM NUMBER
C                FREQ = FREQUENCY
C                GG  = STATISTICAL WEIGHT
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INTEGER S,L,SS,LL
      PARAMETER (PHOT0=2.815D29)
      DIMENSION COEF(3,11),IST(3)
C
      DATA IST/1,2,6/
C
      DATA COEF/
     . -58.229,  4.3965, -0.22134 ,
     . -68.438,  5.7453, -0.26277 ,
     . -67.310,  6.1831, -0.32244 ,
     . -92.020, 10.313 , -0.45090 ,
     . -68.936,  5.2666, -0.15812 ,
     . -63.408,  3.8797, -0.12479 ,
     . -63.778,  4.5102, -0.18213 ,
     . -76.903,  6.3639, -0.21565 ,
     . -61.027,  3.1833, -0.043675,
     . -83.287,  7.1751, -0.20821 ,
     . -83.287,  7.1751, -0.20821 /
C
      SAVE COEF,IST
C
      IF(L.GT.2) GO TO 20
C
C          SELECT BEGINNING AND END OF COEFFICIENTS
C
      SS=(S-1)/2
      LL=2*L
      NSL=IST(N)+LL+SS
C
C          EVALUATE CROSS SECTION
C
      X=LOG(CAS/FREQ)
      CKOEST=EXP(COEF(1,NSL)+X*(COEF(2,NSL)+X*COEF(3,NSL)))/GG
      RETURN
C
C     Hydrogenic expression for L > 2
C      [multiplied by relative population of state (s,l,n), ie.
C       by  stat.weight(s,l)/stat.weight(n)]
C
   20 GN=TWO*N*N
      CKOEST=PHOT0/FREQ/FREQ/FREQ/N**5*(2*L+1)*S/GN
      RETURN
      END
C
C
C     ******************************************************************
C
C
      FUNCTION TOPBAS(FREQ,FREQ0,TYPLV)
C     ==================================
C
C     Procedure calculates the photo-ionisation cross section SIGMA in 
C     [cm^2] at frequency FREQ. FREQ0 is the threshold frequency from
C     level I of ion KI. Threshold cross-sections will be of the order
C     of the numerical value of 10^-18.
C     Opacity-Project (OP) interpolation fit formula
C
      INCLUDE 'IMPLIC.FOR'
      PARAMETER    (E10=2.3025851)
      PARAMETER    (MMAXOP = 200,! maximum number of levels in OP data
     +              MOP    =  15 )! maximum number of fit points per level
      CHARACTER*10  IDLVOP(MMAXOP) ! level identifyer Opacity-Project data
      CHARACTER*10  TYPLV
      LOGICAL*2     LOPREA   
      COMMON /TOPB/ SOP(MOP,MMAXOP) ,! sigma = alog10(sigma/10^-18) of fit point
     +              XOP(MOP,MMAXOP) ,! x = alog10(nu/nu0) of fit point 
     +              NOP(MMAXOP)     ,! number of fit points for current level
     +              NTOTOP          ,! total number of levels in OP data
     +              IDLVOP         ,! level identifyer Opacity-Project data
     +              LOPREA   ! .T. OP data read in; .F. OP data not yer read in
      DIMENSION XFIT(MOP) ,! local array containing x     for OP data
     +          SFIT(MOP)  ! local array containing sigma for OP data
C
      TOPBAS=0.
C
C     Read OP data if not yet done
C
      IF (.NOT.LOPREA) CALL OPDATA
      X = LOG10(FREQ/FREQ0)
      DO IOP = 1,NTOTOP
         IF (IDLVOP(IOP).EQ.TYPLV) THEN
C           level has been detected in OP-data file
            IF (NOP(IOP).LE.0) GO TO 20
            DO IFIT = 1,NOP(IOP)
               XFIT(IFIT) = XOP(IFIT,IOP)
               SFIT(IFIT) = SOP(IFIT,IOP)
            END DO
            SIGM  = YLINTP (X,XFIT,SFIT,NOP(IOP),MOP)
            SIGM  = 1.D-18*EXP(E10*SIGM)
            TOPBAS=SIGM
            GO TO 10
         END IF
      END DO
   10 RETURN
C     Level is not found ,or no data for this level, in RBF.DAT
   20 WRITE (61,100) TYPLV    
  100 FORMAT ('SIGMA.......: OP DATA NOT AVAILABLE FOR LEVEL ',A10)
      RETURN
      END
C

C     ******************************************************************
C
C
      SUBROUTINE OPDATA
C     =================
C
C     Procedure reads photo-ionization cross sections fit coefficients
C     based on Opacity-Project (OP) data from file RBF.DAT
C     Data, as stored, requires linear interpolation.
C
C     Meaning of global variables:
C        NTOTOP    = total number of levels in Opacity Project data
C        IDLVOP() = level identifyer of current level
C        NOP()     = number of fit points for current level
C        XOP(,)    = x     = alog10(nu/nu0)       of fit point
C        SOP(,)    = sigma = alog10(sigma/10^-18) of fit point 
C      
      INCLUDE 'IMPLIC.FOR'
      PARAMETER    (MMAXOP = 200,! maximum number of levels in OP data
     +              MOP    =  15 )! maximum number of fit points per level
      CHARACTER*10  IDLVOP(MMAXOP) ! level identifyer Opacity-Project data
      LOGICAL*2     LOPREA   
      COMMON /TOPB/ SOP(MOP,MMAXOP) ,! sigma = alog10(sigma/10^-18) of fit point
     +              XOP(MOP,MMAXOP) ,! x = alog10(nu/nu0) of fit point 
     +              NOP(MMAXOP)     ,! number of fit points for current level
     +              NTOTOP          ,! total number of levels in OP data
     +              IDLVOP         ,! level identifyer Opacity-Project data
     +              LOPREA   ! .T. OP data read in; .F. OP data not yer read in
      CHARACTER*4 IONID
C
      OPEN (UNIT=40,FILE='RBF.DAT',STATUS='OLD')
C     Skip header
      DO IREAD = 1, 21
         READ (40,*)
      END DO
      IOP = 0
C         = initialize sequential level index op Opacity Project data 
C     Read number of elements in file
      READ (40,*) NEOP
      DO IEOP = 1, NEOP
C        Skip element name header
         DO IREAD = 1, 3
            READ (40,*)
         END DO
C        Read number of ionization stages of current element in  file
         READ (40,*) NIOP
         DO IIOP = 1, NIOP
C           Read ion identifyer, atomic & electron number, # of levels 
C           for current ion
            READ (40,*) IONID, IATOM_OP, IELEC_OP, NLEVEL_OP
            DO ILOP = 1, NLEVEL_OP
C              Increase sequential level index of Opacity Project data
               IOP = IOP+1
C              Read level identifyer and number of sigma fit points
               READ (40,*) IDLVOP(IOP), NOP(IOP)
C              Read normalized log10 frequency and log10 cross section values
               DO IS = 1, NOP(IOP)
                  READ (40,*) INDEX, XOP(IS,IOP), SOP(IS,IOP)
               END DO
            END DO
         END DO
      END DO
      NTOTOP  = IOP
C             = total number of levels in Opacity Project data
      LOPREA  = .TRUE.
C             = set flag as data has been read in
C
      RETURN
      END
C
C
C
C     ******************************************************************
C
C
      FUNCTION YLINTP (XINT,X,Y,N,NTOT)
C     =================================
C
C     linear interpolation routine. Determines YINT = Y(XINT) from 
C     grid Y(X) with N points and dimension NTOT. 
C
      INCLUDE 'IMPLIC.FOR'
      DIMENSION X(NTOT),Y(NTOT)
C
C     bisection (see Numerical Recipes par 3.4 page 90)
      JL = 0
      JU = N+1
10    IF (JU-JL.GT.1) THEN
         JM = (JU+JL)/2
         IF ((X(N).GT.X(1)).EQV.(XINT.GT.X(JM))) THEN
            JL = JM
         ELSE
            JU = JM
         END IF
         GO TO 10
      END IF
      J = JL
      IF (J.EQ.N) J = J-1
      IF (J.EQ.0) J = J+1
CM42 version: enable extrapolation
C
      RC         = (Y(J+1)-Y(J))/(X(J+1)-X(J))
      YLINTP = RC*(XINT-X(J))+Y(J)
C
      RETURN
      END
C
C
C     ****************************************************************
C
C

      SUBROUTINE GFREE0(ID)
C     =====================
C
C     depth-dependent quantities for the hydrogenic free-free Gaunt factor
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'MODELQ.FOR'
      PARAMETER (THET0=5.0404D3,
     *           THET1=UN/THET0,
     *           A0=1.0823d0,
     *           B0=2.98D-2,
     *           C0=6.70D-3,
     *           D0=1.12D-2,
     *           A1=3.9999187d-3,
     *           B1=-7.8622889d-5,
     *           C1=1.070192d0,
     *           A2=6.4628601d-2,
     *           B2=-6.1953813d-4,
     *           C2=2.6061249d-1,
     *           A3=3.7542343d-2,
     *           B3=1.3983474d-5,
     *           C3=5.7917786d-1,
     *           A4=3.4169006d-1,
     *           B4=1.1852264d-2,
     *           XMIN=0.2D0,
     *           XMINI=UN/XMIN,
     *           THMIN=4.0D-2)
C
      T=TEMP(ID)
      THET=UN/MAX(THET0/T,THMIN)
      GF0(ID)=(A0+B0*THET)
      GF1(ID)=(A1+B1*THET)*THET+C1
      GF2(ID)=(A2+B2*THET)*THET+C2
      GF3(ID)=(A3+B3*THET)*THET+C3
      GF4(ID)=A4+B4*THET
      GF5(ID)=C0+D0*THET
      GF6(ID)=GF0(ID)+GF5(ID)*XMINI
C
C     auxiliary quantities for derivatives
C
      THT=THET0/T
      IF(THT.GE.THMIN) THEN 
         THET=THT
         GF0D(ID)=B0*THET1
         GF1D(ID)=(A1+B1*THET*TWO)*THET1
         GF2D(ID)=(A2+B2*THET*TWO)*THET1
         GF3D(ID)=(A3+B3*THET*TWO)*THET1
         GF4D(ID)=B4*THET1
         GF5D(ID)=D0*THET1
         GF6D(ID)=GF0D(ID)+GF5D(ID)*XMINI
       ELSE
         GF0D(ID)=0.
         GF1D(ID)=0.
         GF2D(ID)=0.
         GF3D(ID)=0.
         GF4D(ID)=0.
         GF5D(ID)=0.
         GF6D(ID)=0.
      END IF
      RETURN
      END
C
C
C     ****************************************************************
C
C

      FUNCTION GFREE1(ID,X)
C     =====================
C
C     Hydrogenic free-free Gaunt factor, for depth ID,
C     frequency FR, and charge CH
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'MODELQ.FOR'
      PARAMETER (XMIN=0.2D0,
     *           XMINI=UN/XMIN)
C
      IF(X.LT.UN) THEN
         GFREE1=((GF4(ID)*X-GF3(ID))*X+GF2(ID))*X+GF1(ID)
       ELSE IF(X.LT.XMINI) THEN
         GFREE1=GF0(ID)+GF5(ID)*X
       ELSE
         GFREE1=GF6(ID)
      END IF
      RETURN
      END
C
C
C     ****************************************************************
C
C

      SUBROUTINE GFREED(ID,FR,CH,GFR,GFRD)
C     ====================================
C
C     Hydrogenic free-free Gaunt factor (GFR) and its derivative (GFRD)
C     for depth ID, frequency FR, and charge CH
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'MODELQ.FOR'
      PARAMETER (C14=2.997925D14,
     *           XMIN=0.2D0,
     *           XMINI=UN/XMIN)
C
      X=C14*CH/FR
      IF(X.LT.UN) THEN
         GFR=((GF4(ID)*X-GF3(ID))*X+GF2(ID))*X+GF1(ID)
         GFRD=((GF4D(ID)*X-GF3D(ID))*X+GF2D(ID))*X+GF1D(ID)
       ELSE IF(X.LT.XMINI) THEN
         GFR=GF0(ID)+GF5(ID)*X
         GFRD=GF0D(ID)+GF5D(ID)*X
       ELSE
         GFR=GF6(ID)
         GFRD=GF6D(ID)
      END IF
      RETURN
      END
C
C
C     ****************************************************************
C
C

      FUNCTION FFCROS(IEL,IFR,T,FR)
C     ==================================
C
C     Non-standard evaluation of free-free cross section;
C     a user supplied procedure
C
      INCLUDE 'IMPLIC.FOR'
      FFCROS=0.
      if(iel.eq.0.or.ifr.eq.0) return
      t1=t
      fr1=fr
      RETURN
      END
C
C
C     ****************************************************************
C
C

      FUNCTION SBFHMI_old(FR)
C     ===================
C
C     Bound-free cross-section for H- (negative hydrogen ion)
C
      INCLUDE 'IMPLIC.FOR'
      PARAMETER (UN=1.D0)
      SBFHMI_old=0.
      FR0=1.8259D14
      IF(FR.LT.FR0) RETURN
      IF(FR.LT.2.111D14) GO TO 10
      X=2.997925D15/FR
      SBFHMI=(6.80133D-3+X*(1.78708D-1+X*(1.6479D-1+X*(-2.04842D-2+X*
     *        5.95244D-4))))*1.D-17
      SBFHMI_old=sbfhmi
      RETURN
   10 X=2.997925D15*(UN/FR0-UN/FR)
      SBFHMI=(2.69818D-1+X*(2.2019D-1+X*(-4.11288D-2+X*2.73236D-3)))
     1       *X*1.D-17
      SBFHMI_old=sbfhmi
      RETURN
      END
C
C
C ********************************************************************
C
C
      FUNCTION SBFHMI(FR)
C     ===================
C    
C     Bound-free cross-section for H- (negative hydrogen ion)
C     Taken from Kurucz ATLAS9
C
C     FROM MATHISEN (1984), AFTER WISHART(1979) AND BROAD AND REINHARDT (1976)
C
      INCLUDE 'IMPLIC.FOR'
      DIMENSION WBF(85),BF(85)
      DATA WBF/  18.00,  19.60,  21.40,  23.60,  26.40,  29.80,  34.30,
     1   40.40,  49.10,  62.60, 111.30, 112.10, 112.67, 112.95, 113.05,
     2  113.10, 113.20, 113.23, 113.50, 114.40, 121.00, 139.00, 164.00,
     3  175.00, 200.00, 225.00, 250.00, 275.00, 300.00, 325.00, 350.00,
     4  375.00, 400.00, 425.00, 450.00, 475.00, 500.00, 525.00, 550.00,
     5  575.00, 600.00, 625.00, 650.00, 675.00, 700.00, 725.00, 750.00,
     6  775.00, 800.00, 825.00, 850.00, 875.00, 900.00, 925.00, 950.00,
     7  975.00,1000.00,1025.00,1050.00,1075.00,1100.00,1125.00,1150.00,
     8 1175.00,1200.00,1225.00,1250.00,1275.00,1300.00,1325.00,1350.00,
     9 1375.00,1400.00,1425.00,1450.00,1475.00,1500.00,1525.00,1550.00,
     A 1575.00,1600.00,1610.00,1620.00,1630.00,1643.91/
      DATA BF/   0.067,  0.088,  0.117,  0.155,  0.206,  0.283,  0.414,
     1   0.703,   1.24,   2.33,  11.60,  13.90,  24.30,  66.70,  95.00,
     2   56.60,  20.00,  14.60,   8.50,   7.10,   5.43,   5.91,   7.29,
     3   7.918,  9.453,  11.08,  12.75,  14.46,  16.19,  17.92,  19.65,
     4   21.35,  23.02,  24.65,  26.24,  27.77,  29.23,  30.62,  31.94,
     5   33.17,  34.32,  35.37,  36.32,  37.17,  37.91,  38.54,  39.07,
     6   39.48,  39.77,  39.95,  40.01,  39.95,  39.77,  39.48,  39.06,
     7   38.53,  37.89,  37.13,  36.25,  35.28,  34.19,  33.01,  31.72,
     8   30.34,  28.87,  27.33,  25.71,  24.02,  22.26,  20.46,  18.62,
     9   16.74,  14.85,  12.95,  11.07,  9.211,  7.407,  5.677,  4.052,
     A   2.575,  1.302, 0.8697, 0.4974, 0.1989,    0. /
C    Bell and Berrington J.Phys.B,vol. 20, 801-806,1987.
c
      HMINBF=0.
      IF(FR.GT.1.82365E14) THEN
         WAVE=2.99792458E17/FR
         HMINBF=YLINTP(WAVE,WBF,BF,85,85)*1.E-18
      END IF
      SBFHMI=HMINBF
      RETURN
      END
C
C
C     ****************************************************************
C
C

      SUBROUTINE COLIS(ID,T,COL)
C     ==========================
C
C     Driving procedure for evaluation of the collisional rates
C
C     Input:  T   - temperature
C     Output: COL - array of quantities proportional to the
C                   collisional rates in all transitions
C                   for a given temperature (ie. at a given depth)
C                   Precisely,  COL(IT) * n(el)  is the upward
C                   collisional rate in the transition IT
C
C     Procedure COLIS calls procedures COLH and COLHE for evaluating
C     the collisional rates in hydrogen and helium,
C     and itself evaluates collisional rates for other species
C
C     Evaluation is controlled by input parameter ICOL(ITR).
C     Meaning of ICOL for all species, excluding hydrogen and helium:
C
C      a) for ionization
C         ICOL = 0  - Seatons formula ; here the value of the photo-
C                     ionization cross section at the threshold is
C                     transmitted in array OSC0
C              = 1  - Allen's formula; again, OSC0 has the meaning of
C                     the necessary multiplicative parameter
C              = 2  - the so-called SIMPLE1 mode - see below
C              = 3  - the so-called SIMPLE2 mode - see below
C      b) for excitation
C         ICOL = 0  - Van Regemorter formula, with standard g(bar)=0.25
C              = 1  - Van Regemorter formula, with "exact" g(bar)
C              = 2  - the so-called SIMPLE1 mode - see below
C              = 3  - the so-called SIMPLE2 mode - see below
C              = 4  - Eissner-Seaton formula - see below
C
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'ATOMIC.FOR'
      INCLUDE 'MODELQ.FOR'
      INCLUDE 'ODFPAR.FOR'
      PARAMETER (EXPIA1=-0.57721566,EXPIA2=0.99999193,
     *           EXPIA3=-0.24991055,EXPIA4=0.05519968,
     *           EXPIA5=-0.00976004,EXPIA6=0.00107857,
     *           EXPIB1=0.2677734343,EXPIB2=8.6347608925,
     *           EXPIB3=18.059016973,EXPIB4=8.5733287401,
     *           EXPIC1=3.9584969228,EXPIC2=21.0996530827,
     *           EXPIC3=25.6329561486,EXPIC4=9.5733223454)
      DIMENSION COL(MTRANS)
C
      COMMON/CTRTEMP/ te
C

      CREGER(X,U,A,GG)=19.7363*X*EXP(-U)/U*GG*A
      CSEATN(X,U,A)=1.55D13*X/ABS(U)*EXP(-U)*A
      CALLEN(X,U,A)=X*A*EXP(-U)/U/U
      CSMPL1(X,U,A)=5.465D-11*X*EXP(-U)*A
      CSMPL2(X,U,A)=5.465D-11*X*EXP(-U)*A*(1.+U)
      CUPSX(X,U,A)=8.631D-6/X*EXP(-U)*A
C
      DO 10 I=1,NTRANS
   10    COL(I)=0.
      HKT=HK/T
      SRT=SQRT(T)
      T32=UN/T/SRT
      TK=HKT/H
      CSTD=0.25
      T0=TEFF
      if(t0.gt.0.) then
         TT0=UN-T/T0
         SRT0=SQRT(T0)/SRT
      end if
C
C     Call procedures COLH and COLHE for hydrogen and helium
C
      IF(IATH.NE.0) CALL COLH(ID,T,COL)
      IF(IATHE.NE.0) CALL COLHE(T,COL)
C
C     Loop over all explicit species, excluding hydrogen and helium
C
      DO 100 IAT=1,NATOM
      IF(IAT.EQ.IATH.OR.IAT.EQ.IATHE) GO TO 100
      N0I=N0A(IAT)
      NKI=NKA(IAT)
      DO 50 I=N0I,NKI-1
         IE=IEL(I)
         DO 40 J=I+1,NKI
            IT=ITRA(I,J)
            IF(IT.EQ.0) GO TO 40
            IC=ICOL(IT)
            C1=OSC0(IT)
            C2=CPAR(IT)
            U0=FR0(IT)*HKT
            IF(LINE(IT)) GO TO 30

C
C ********* Charge-transfer reactions
C
             IF(IC.GE.10) THEN
C               radiative charge transfer ionization of neutrals in the ground
C               state with protons
                te=T
                COL(IT)=HCTIon(1,NUMAT(IAT))/
     *                ELEC(ID)*POPUL(NNEXT(IATH),ID)
                IC=IC-10

           ENDIF
             

C
C ********* Collisional ionization
C
            IF(IC.EQ.0) THEN
               COL(IT)=CSEATN(UN/SRT,U0,C1)
             ELSE IF(IC.EQ.1) THEN
               COL(IT)=CALLEN(T32,U0,C1)
             ELSE IF(IC.EQ.2) THEN
               COL(IT)=CSMPL1(SRT,U0,C1)
             ELSE IF(IC.EQ.3) THEN
               COL(IT)=CSMPL2(SRT,U0,C1)
             ELSE IF(IC.EQ.4) THEN
               ia=numat(iatm(i))
               COL(IT)=cion(ia,iz(iel(i)),enion(i)*6.24298e11,t)
             ELSE IF(IC.EQ.5) THEN
               ia=numat(iatm(i))
               izc=iz(ie)
               rno=16.
               ii=i-nfirst(ie)+1
               call irc(ii,t,izc,rno,cs)
               col(it)=cs
             ELSE IF(IC.LT.0) THEN
               CALL CSPEC(I,J,IC,C1,C2,U0,T,COL(IT))
            END IF
            if(ic.eq.4) go to 40
C
C      collisional excitations from level I to higher, non-explicit
C      levels are lumped into the collisional ionization rate
C      (the so-called modified collision ionization rate)
C
            N0Q=NQUANT(NLAST(IE))+1
            N1Q=ICUP(IE)
            IF(N1Q.EQ.0) GO TO 40
            IQ=NQUANT(I)
            REL=G(I)/2./IQ/IQ
            DO 20 JQ=N0Q,N1Q
               XJ=JQ
               U0=(ENION(I)-EH/XJ/XJ)*TK
               IF(JQ.LE.20) CC1=OSH(IQ,JQ)*REL
               IF(JQ.GT.20) CC1=OSH(IQ,20)*(20./XJ)**3*REL
               GG=CSTD
               if(u0.gt.35.) go to 20
               IF(U0.LE.UN) THEN
                  EXPIU0=-LOG(U0)+EXPIA1+U0*(EXPIA2+
     *            U0*(EXPIA3+U0*(EXPIA4+
     *            U0*(EXPIA5+U0*EXPIA6))))
                ELSE
                  EXPIU0=EXP(-U0)*((EXPIB1+U0*(EXPIB2+
     *            U0*(EXPIB3+
     *            U0*(EXPIB4+U0))))/(EXPIC1+U0*(EXPIC2+
     *            U0*(EXPIC3+U0*(EXPIC4+U0)))))/U0
               END IF
               GG0=0.276*EXP(U0)*EXPIU0
               IF(GG0.GT.GG) GG=GG0
               COL(IT)=COL(IT)+CREGER(T32,U0,CC1,GG)
   20       CONTINUE
            GO TO 40
C
C ********* Collisional excitation
C
   30       CONTINUE
            IF(IC.LE.1.AND.IC.GT.0) THEN
               GG=CSTD
               IF(IC.EQ.1) GG=C2
               IF(U0.LE.UN) THEN
                  EXPIU0=-LOG(U0)+EXPIA1+
     *            U0*(EXPIA2+U0*(EXPIA3+U0*(EXPIA4+
     *            U0*(EXPIA5+U0*EXPIA6))))
                ELSE
                  EXPIU0=EXP(-U0)*((EXPIB1+U0*
     *            (EXPIB2+U0*(EXPIB3+
     *            U0*(EXPIB4+U0))))/(EXPIC1+U0*(EXPIC2+
     *            U0*(EXPIC3+U0*(EXPIC4+U0)))))/U0
               END IF
               GG0=0.276*EXP(U0)*EXPIU0
               IF(GG0.GT.GG) GG=GG0
               COL(IT)=CREGER(T32,U0,C1,GG)
             ELSE IF(IC.EQ.2) THEN
               COL(IT)=CSMPL1(SRT,U0,C1*C2)
             ELSE IF(IC.EQ.3) THEN
               COL(IT)=CSMPL2(SRT,U0,C2)
             ELSE IF(IC.EQ.4) THEN
               COL(IT)=CUPSX(SRT,U0,C2/G(I))
             ELSE IF(IC.EQ.9) THEN
               COL(IT)=OMECOL(I,J)*SRT0*EXP(-U0*TT0)
             ELSE IF(IC.LT.0) THEN
               CALL CSPEC(I,J,IC,C1,C2,U0,T,COL(IT))
            END IF
   40    CONTINUE
   50 CONTINUE
  100 CONTINUE
      RETURN
      END

C
C     ****************************************************************
C
C
*     Code by Jim Kingdon, in collaboration with G.J. Ferland
      function HCTRecom(ion,nelem)
*     ion is stage of ionization, 2 for the ion going to the atom
*     nelem is atomic number of element, 2 up to 30
*     Example:  O+ + H => O + H+ is HCTRecom(2,8)
*      integer ion , nelem
      INCLUDE 'IMPLIC.FOR'
      common/CTRTEMP/ te
      common/CTRecomb/ CTRecomb(6,4,30)
*
*     local variables
*      real tused 
*      integer ipIon
*
      ipIon = ion - 1
*
      if( ipIon.gt.4 ) then
*       use statistical charge transfer for ion > 4
        HCTRecom = 1.92e-9 * ipIon
        return
      endif
*
*     Make sure te is between temp. boundaries; set constant outside of range
      tused = max( te,CTRecomb(5,ipIon,nelem) )
      tused = min( tused , CTRecomb(6,ipIon,nelem) )
      tused = tused * 1e-4
*
*     the interpolation equation
      HCTRecom = CTRecomb(1,ipIon,nelem)* 1e-9 * 
     1 (tused**CTRecomb(2,ipIon,nelem)) *
     2 (1. + 
     3 CTRecomb(3,ipIon,nelem) * exp(CTRecomb(4,ipIon,nelem)*tused) )
*
      end
******************************************************************************
      function HCTIon(ion,nelem)
*     ion is stage of ionization, 1 for atom
*     nelem is atomic number of element, 2 up to 30
*     Example:  O + H+ => O+ + H is HCTIon(1,8)
*      integer ion , nelem
      INCLUDE 'IMPLIC.FOR'
      common/CTRTEMP/ te
      common/CTIon/ CTIon(7,4,30)
*
*     local variables
*      real tused 
*      integer ipIon
*
      ipIon = ion
*
*     Make sure te is between temp. boundaries; set constant outside of range
      tused = max( te,CTIon(5,ipIon,nelem) )
      tused = min( tused , CTIon(6,ipIon,nelem) )
      tused = tused * 1e-4
*
*     the interpolation equation
      HCTIon = CTIon(1,ipIon,nelem)* 1e-9 * 
     1 (tused**CTIon(2,ipIon,nelem)) *
     2 (1. + 
     3 CTIon(3,ipIon,nelem) * exp(CTIon(4,ipIon,nelem)*tused) ) *
     4 exp(-CTIon(7,ipIon,nelem)/tused)
*
      end

********************************************************************************
      block data ctdata
*
*      real CTIon
*     second dimension is ionization stage,
*     1=+0 for parent, etc
*     third dimension is atomic number of atom
      INCLUDE 'IMPLIC.FOR'
      common/CTIon/ CTIon(7,4,30)
*      real CTRecomb
*     second dimension is ionization stage,
*     1=+1 for parent, etc
*     third dimension is atomic number of atom
      common/CTRecomb/ CTRecomb(6,4,30)
*
*     local variables
*      integer i
*
*     digital form of the fits to the charge transfer
*     ionization rate coefficients 
*
*     Note: First parameter is in units of 1e-9!
*     Note: Seventh parameter is in units of 1e4 K
*     ionization
      data (CTIon(i,1,3),i=1,7)/2.84e-3,1.99,375.54,-54.07,1e2,1e4,0.0/
      data (CTIon(i,2,3),i=1,7)/7*0./
      data (CTIon(i,3,3),i=1,7)/7*0./
      data (CTIon(i,1,4),i=1,7)/7*0./
      data (CTIon(i,2,4),i=1,7)/7*0./
      data (CTIon(i,3,4),i=1,7)/7*0./
      data (CTIon(i,1,5),i=1,7)/7*0./
      data (CTIon(i,2,5),i=1,7)/7*0./
      data (CTIon(i,3,5),i=1,7)/7*0./
      data (CTIon(i,1,6),i=1,7)/1.07e-6,3.15,176.43,-4.29,1e3,1e5,0.0/
      data (CTIon(i,2,6),i=1,7)/7*0./
      data (CTIon(i,3,6),i=1,7)/7*0./
      data (CTIon(i,1,7),i=1,7)/4.55e-3,-0.29,-0.92,-8.38,1e2,5e4,1.086/
      data (CTIon(i,2,7),i=1,7)/7*0./
      data (CTIon(i,3,7),i=1,7)/7*0./
      data (CTIon(i,1,8),i=1,7)/7.40e-2,0.47,24.37,-0.74,1e1,1e4,0.023/
      data (CTIon(i,2,8),i=1,7)/7*0./
      data (CTIon(i,3,8),i=1,7)/7*0./
      data (CTIon(i,1,9),i=1,7)/7*0./
      data (CTIon(i,2,9),i=1,7)/7*0./
      data (CTIon(i,3,9),i=1,7)/7*0./
      data (CTIon(i,1,10),i=1,7)/7*0./
      data (CTIon(i,2,10),i=1,7)/7*0./
      data (CTIon(i,3,10),i=1,7)/7*0./
      data (CTIon(i,1,11),i=1,7)/3.34e-6,9.31,2632.31,-3.04,1e3,2e4,0.0/
      data (CTIon(i,2,11),i=1,7)/7*0./
      data (CTIon(i,3,11),i=1,7)/7*0./
      data (CTIon(i,1,12),i=1,7)/9.76e-3,3.14,55.54,-1.12,5e3,3e4,0.0/
      data (CTIon(i,2,12),i=1,7)/7.60e-5,0.00,-1.97,-4.32,1e4,3e5,1.670/
      data (CTIon(i,3,12),i=1,7)/7*0./
      data (CTIon(i,1,13),i=1,7)/7*0./
      data (CTIon(i,2,13),i=1,7)/7*0./
      data (CTIon(i,3,13),i=1,7)/7*0./
      data (CTIon(i,1,14),i=1,7)/0.92,1.15,0.80,-0.24,1e3,2e5,0.0/
      data (CTIon(i,2,14),i=1,7)/2.26,7.36e-2,-0.43,-0.11,2e3,1e5,
     1 3.031/
      data (CTIon(i,3,14),i=1,7)/7*0./
      data (CTIon(i,1,15),i=1,7)/7*0./
      data (CTIon(i,2,15),i=1,7)/7*0./
      data (CTIon(i,3,15),i=1,7)/7*0./
      data (CTIon(i,1,16),i=1,7)/1.00e-5,0.00,0.00,0.00,1e3,1e4,0.0/
      data (CTIon(i,2,16),i=1,7)/7*0./
      data (CTIon(i,3,16),i=1,7)/7*0./
      data (CTIon(i,1,17),i=1,7)/7*0./
      data (CTIon(i,2,17),i=1,7)/7*0./
      data (CTIon(i,3,17),i=1,7)/7*0./
      data (CTIon(i,1,18),i=1,7)/7*0./
      data (CTIon(i,2,18),i=1,7)/7*0./
      data (CTIon(i,3,18),i=1,7)/7*0./
      data (CTIon(i,1,19),i=1,7)/7*0./
      data (CTIon(i,2,19),i=1,7)/7*0./
      data (CTIon(i,3,19),i=1,7)/7*0./
      data (CTIon(i,1,20),i=1,7)/7*0./
      data (CTIon(i,2,20),i=1,7)/7*0./
      data (CTIon(i,3,20),i=1,7)/7*0./
      data (CTIon(i,1,21),i=1,7)/7*0./
      data (CTIon(i,2,21),i=1,7)/7*0./
      data (CTIon(i,3,21),i=1,7)/7*0./
      data (CTIon(i,1,22),i=1,7)/7*0./
      data (CTIon(i,2,22),i=1,7)/7*0./
      data (CTIon(i,3,22),i=1,7)/7*0./
      data (CTIon(i,1,23),i=1,7)/7*0./
      data (CTIon(i,2,23),i=1,7)/7*0./
      data (CTIon(i,3,23),i=1,7)/7*0./
      data (CTIon(i,1,24),i=1,7)/7*0./
      data (CTIon(i,2,24),i=1,7)/4.39,0.61,-0.89,-3.56,1e3,3e4,3.349/
      data (CTIon(i,3,24),i=1,7)/7*0./
      data (CTIon(i,1,25),i=1,7)/7*0./
      data (CTIon(i,2,25),i=1,7)/2.83e-1,6.80e-3,6.44e-2,-9.70,1e3,3e4,
     1 2.368/
      data (CTIon(i,3,25),i=1,7)/7*0./
      data (CTIon(i,1,26),i=1,7)/7*0./
      data (CTIon(i,2,26),i=1,7)/2.10,7.72e-2,-0.41,-7.31,1e4,1e5,3.005/
      data (CTIon(i,3,26),i=1,7)/7*0./
      data (CTIon(i,1,27),i=1,7)/7*0./
      data (CTIon(i,2,27),i=1,7)/1.20e-2,3.49,24.41,-1.26,1e3,3e4,4.044/
      data (CTIon(i,3,27),i=1,7)/7*0./
      data (CTIon(i,1,28),i=1,7)/7*0./
      data (CTIon(i,2,28),i=1,7)/7*0./
      data (CTIon(i,3,28),i=1,7)/7*0./
      data (CTIon(i,1,29),i=1,7)/7*0./
      data (CTIon(i,2,29),i=1,7)/7*0./
      data (CTIon(i,3,29),i=1,7)/7*0./
      data (CTIon(i,1,30),i=1,7)/7*0./
      data (CTIon(i,2,30),i=1,7)/7*0./
      data (CTIon(i,3,30),i=1,7)/7*0./
*
*     digital form of the fits to the charge transfer
*     recombination rate coefficients (total)
*
*     Note: First parameter is in units of 1e-9!
*     recombination
      data (CTRecomb(i,1,2),i=1,6)/7.47e-6,2.06,9.93,-3.89,6e3,1e5/
      data (CTRecomb(i,2,2),i=1,6)/1.00e-5,0.,0.,0.,1e3,1e7/
      data (CTRecomb(i,1,3),i=1,6)/6*0./
      data (CTRecomb(i,2,3),i=1,6)/1.26,0.96,3.02,-0.65,1e3,3e4/
      data (CTRecomb(i,3,3),i=1,6)/1.00e-5,0.,0.,0.,2e3,5e4/
      data (CTRecomb(i,1,4),i=1,6)/6*0./
      data (CTRecomb(i,2,4),i=1,6)/1.00e-5,0.,0.,0.,2e3,5e4/
      data (CTRecomb(i,3,4),i=1,6)/1.00e-5,0.,0.,0.,2e3,5e4/
      data (CTRecomb(i,4,4),i=1,6)/5.17,0.82,-0.69,-1.12,2e3,5e4/
      data (CTRecomb(i,1,5),i=1,6)/6*0./
      data (CTRecomb(i,2,5),i=1,6)/2.00e-2,0.,0.,0.,1e3,1e9/
      data (CTRecomb(i,3,5),i=1,6)/1.00e-5,0.,0.,0.,2e3,5e4/
      data (CTRecomb(i,4,5),i=1,6)/2.74,0.93,-0.61,-1.13,2e3,5e4/
      data (CTRecomb(i,1,6),i=1,6)/4.88e-7,3.25,-1.12,-0.21,5.5e3,1e5/
      data (CTRecomb(i,2,6),i=1,6)/1.67e-4,2.79,304.72,-4.07,5e3,5e4/
      data (CTRecomb(i,3,6),i=1,6)/3.25,0.21,0.19,-3.29,1e3,1e5/
      data (CTRecomb(i,4,6),i=1,6)/332.46,-0.11,-9.95e-1,-1.58e-3,1e1,
     1 1e5/
      data (CTRecomb(i,1,7),i=1,6)/1.01e-3,-0.29,-0.92,-8.38,1e2,5e4/
      data (CTRecomb(i,2,7),i=1,6)/3.05e-1,0.60,2.65,-0.93,1e3,1e5/
      data (CTRecomb(i,3,7),i=1,6)/4.54,0.57,-0.65,-0.89,1e1,1e5/
      data (CTRecomb(i,4,7),i=1,6)/2.95,0.55,-0.39,-1.07,1e3,1e6/
      data (CTRecomb(i,1,8),i=1,6)/1.04,3.15e-2,-0.61,-9.73,1e1,1e4/
      data (CTRecomb(i,2,8),i=1,6)/1.04,0.27,2.02,-5.92,1e2,1e5/
      data (CTRecomb(i,3,8),i=1,6)/3.98,0.26,0.56,-2.62,1e3,5e4/
      data (CTRecomb(i,4,8),i=1,6)/2.52e-1,0.63,2.08,-4.16,1e3,3e4/
      data (CTRecomb(i,1,9),i=1,6)/6*0./
      data (CTRecomb(i,2,9),i=1,6)/1.00e-5,0.,0.,0.,2e3,5e4/
      data (CTRecomb(i,3,9),i=1,6)/9.86,0.29,-0.21,-1.15,2e3,5e4/
      data (CTRecomb(i,4,9),i=1,6)/7.15e-1,1.21,-0.70,-0.85,2e3,5e4/
      data (CTRecomb(i,1,10),i=1,6)/6*0./
      data (CTRecomb(i,2,10),i=1,6)/1.00e-5,0.,0.,0.,5e3,5e4/
      data (CTRecomb(i,3,10),i=1,6)/14.73,4.52e-2,-0.84,-0.31,5e3,5e4/
      data (CTRecomb(i,4,10),i=1,6)/6.47,0.54,3.59,-5.22,1e3,3e4/
      data (CTRecomb(i,1,11),i=1,6)/6*0./
      data (CTRecomb(i,2,11),i=1,6)/1.00e-5,0.,0.,0.,2e3,5e4/
      data (CTRecomb(i,3,11),i=1,6)/1.33,1.15,1.20,-0.32,2e3,5e4/
      data (CTRecomb(i,4,11),i=1,6)/1.01e-1,1.34,10.05,-6.41,2e3,5e4/
      data (CTRecomb(i,1,12),i=1,6)/6*0./
      data (CTRecomb(i,2,12),i=1,6)/8.58e-5,2.49e-3,2.93e-2,-4.33,1e3,
     1 3e4/
      data (CTRecomb(i,3,12),i=1,6)/6.49,0.53,2.82,-7.63,1e3,3e4/
      data (CTRecomb(i,4,12),i=1,6)/6.36,0.55,3.86,-5.19,1e3,3e4/
      data (CTRecomb(i,1,13),i=1,6)/6*0./
      data (CTRecomb(i,2,13),i=1,6)/1.00e-5,0.,0.,0.,1e3,3e4/
      data (CTRecomb(i,3,13),i=1,6)/7.11e-5,4.12,1.72e4,-22.24,1e3,3e4/
      data (CTRecomb(i,4,13),i=1,6)/7.52e-1,0.77,6.24,-5.67,1e3,3e4/
      data (CTRecomb(i,1,14),i=1,6)/6*0./
      data (CTRecomb(i,2,14),i=1,6)/6.77,7.36e-2,-0.43,-0.11,5e2,1e5/
      data (CTRecomb(i,3,14),i=1,6)/4.90e-1,-8.74e-2,-0.36,-0.79,1e3,
     1 3e4/
      data (CTRecomb(i,4,14),i=1,6)/7.58,0.37,1.06,-4.09,1e3,5e4/
      data (CTRecomb(i,1,15),i=1,6)/6*0./
      data (CTRecomb(i,2,15),i=1,6)/1.74e-4,3.84,36.06,-0.97,1e3,3e4/
      data (CTRecomb(i,3,15),i=1,6)/9.46e-2,-5.58e-2,0.77,-6.43,1e3,3e4/
      data (CTRecomb(i,4,15),i=1,6)/5.37,0.47,2.21,-8.52,1e3,3e4/
      data (CTRecomb(i,1,16),i=1,6)/3.82e-7,11.10,2.57e4,-8.22,1e3,1e4/
      data (CTRecomb(i,2,16),i=1,6)/1.00e-5,0.,0.,0.,1e3,3e4/
      data (CTRecomb(i,3,16),i=1,6)/2.29,4.02e-2,1.59,-6.06,1e3,3e4/
      data (CTRecomb(i,4,16),i=1,6)/6.44,0.13,2.69,-5.69,1e3,3e4/
      data (CTRecomb(i,1,17),i=1,6)/6*0./
      data (CTRecomb(i,2,17),i=1,6)/1.00e-5,0.,0.,0.,1e3,3e4/
      data (CTRecomb(i,3,17),i=1,6)/1.88,0.32,1.77,-5.70,1e3,3e4/
      data (CTRecomb(i,4,17),i=1,6)/7.27,0.29,1.04,-10.14,1e3,3e4/
      data (CTRecomb(i,1,18),i=1,6)/6*0./
      data (CTRecomb(i,2,18),i=1,6)/1.00e-5,0.,0.,0.,1e3,3e4/
      data (CTRecomb(i,3,18),i=1,6)/4.57,0.27,-0.18,-1.57,1e3,3e4/
      data (CTRecomb(i,4,18),i=1,6)/6.37,0.85,10.21,-6.22,1e3,3e4/
      data (CTRecomb(i,1,19),i=1,6)/6*0./
      data (CTRecomb(i,2,19),i=1,6)/1.00e-5,0.,0.,0.,1e3,3e4/
      data (CTRecomb(i,3,19),i=1,6)/4.76,0.44,-0.56,-0.88,1e3,3e4/
      data (CTRecomb(i,4,19),i=1,6)/1.00e-5,0.,0.,0.,1e3,3e4/
      data (CTRecomb(i,1,20),i=1,6)/6*0./
      data (CTRecomb(i,2,20),i=1,6)/0.,0.,0.,0.,1e1,1e9/
      data (CTRecomb(i,3,20),i=1,6)/3.17e-2,2.12,12.06,-0.40,1e3,3e4/
      data (CTRecomb(i,4,20),i=1,6)/2.68,0.69,-0.68,-4.47,1e3,3e4/
      data (CTRecomb(i,1,21),i=1,6)/6*0./
      data (CTRecomb(i,2,21),i=1,6)/0.,0.,0.,0.,1e1,1e9/
      data (CTRecomb(i,3,21),i=1,6)/7.22e-3,2.34,411.50,-13.24,1e3,3e4/
      data (CTRecomb(i,4,21),i=1,6)/1.20e-1,1.48,4.00,-9.33,1e3,3e4/
      data (CTRecomb(i,1,22),i=1,6)/6*0./
      data (CTRecomb(i,2,22),i=1,6)/0.,0.,0.,0.,1e1,1e9/
      data (CTRecomb(i,3,22),i=1,6)/6.34e-1,6.87e-3,0.18,-8.04,1e3,3e4/
      data (CTRecomb(i,4,22),i=1,6)/4.37e-3,1.25,40.02,-8.05,1e3,3e4/
      data (CTRecomb(i,1,23),i=1,6)/6*0./
      data (CTRecomb(i,2,23),i=1,6)/1.00e-5,0.,0.,0.,1e3,3e4/
      data (CTRecomb(i,3,23),i=1,6)/5.12,-2.18e-2,-0.24,-0.83,1e3,3e4/
      data (CTRecomb(i,4,23),i=1,6)/1.96e-1,-8.53e-3,0.28,-6.46,1e3,3e4/
      data (CTRecomb(i,1,24),i=1,6)/6*0./
      data (CTRecomb(i,2,24),i=1,6)/5.27e-1,0.61,-0.89,-3.56,1e3,3e4/
      data (CTRecomb(i,3,24),i=1,6)/10.90,0.24,0.26,-11.94,1e3,3e4/
      data (CTRecomb(i,4,24),i=1,6)/1.18,0.20,0.77,-7.09,1e3,3e4/
      data (CTRecomb(i,1,25),i=1,6)/6*0./
      data (CTRecomb(i,2,25),i=1,6)/1.65e-1,6.80e-3,6.44e-2,-9.70,1e3,
     1 3e4/
      data (CTRecomb(i,3,25),i=1,6)/14.20,0.34,-0.41,-1.19,1e3,3e4/
      data (CTRecomb(i,4,25),i=1,6)/4.43e-1,0.91,10.76,-7.49,1e3,3e4/
      data (CTRecomb(i,1,26),i=1,6)/6*0./
      data (CTRecomb(i,2,26),i=1,6)/1.26,7.72e-2,-0.41,-7.31,1e3,1e5/
      data (CTRecomb(i,3,26),i=1,6)/3.42,0.51,-2.06,-8.99,1e3,1e5/
      data (CTRecomb(i,4,26),i=1,6)/14.60,3.57e-2,-0.92,-0.37,1e3,3e4/
      data (CTRecomb(i,1,27),i=1,6)/6*0./
      data (CTRecomb(i,2,27),i=1,6)/5.30,0.24,-0.91,-0.47,1e3,3e4/
      data (CTRecomb(i,3,27),i=1,6)/3.26,0.87,2.85,-9.23,1e3,3e4/
      data (CTRecomb(i,4,27),i=1,6)/1.03,0.58,-0.89,-0.66,1e3,3e4/
      data (CTRecomb(i,1,28),i=1,6)/6*0./
      data (CTRecomb(i,2,28),i=1,6)/1.05,1.28,6.54,-1.81,1e3,1e5/
      data (CTRecomb(i,3,28),i=1,6)/9.73,0.35,0.90,-5.33,1e3,3e4/
      data (CTRecomb(i,4,28),i=1,6)/6.14,0.25,-0.91,-0.42,1e3,3e4/
      data (CTRecomb(i,1,29),i=1,6)/6*0./
      data (CTRecomb(i,2,29),i=1,6)/1.47e-3,3.51,23.91,-0.93,1e3,3e4/
      data (CTRecomb(i,3,29),i=1,6)/9.26,0.37,0.40,-10.73,1e3,3e4/
      data (CTRecomb(i,4,29),i=1,6)/11.59,0.20,0.80,-6.62,1e3,3e4/
      data (CTRecomb(i,1,30),i=1,6)/6*0./
      data (CTRecomb(i,2,30),i=1,6)/1.00e-5,0.,0.,0.,1e3,3e4/
      data (CTRecomb(i,3,30),i=1,6)/6.96e-4,4.24,26.06,-1.24,1e3,3e4/
      data (CTRecomb(i,4,30),i=1,6)/1.33e-2,1.56,-0.92,-1.20,1e3,3e4/
*
      end
                
C
C
C     ****************************************************************
C
C


      SUBROUTINE COLH(ID,T,COL)
C     =========================
C
C     Hydrogen collision rates
C
C     All standard expressions  are taken from Mihalas, Heasley, and
C     Auer, NCAR-TN/STR-104 (1975)
C
C     New expressions (also from Mihalas) for collisional ionization
C     for first 10 levels taken from Klaus Werner.
C
C     New standard expressions from Giovanardi et al. (1987, AAS, 70, 29)
C      for collisional excitation (valid from 3000K to 500000K)
C
C     Meaning of ICOL:
C      a) for ionization - .ge.0  - standard expression
C                           <  0  - non-standard, user suplied formula
C      b) for 1 - 2 transition =0 - standard theoretical formula
C                             = 1 - experimental fit (formula quted in
C                                   Mihalas et al.)
C                             = 2 - formula by Crandall et al (procedure
C                                   CEH12)
C      c) for all other line transitions
C                           .ge.0 - standard expression
C                            <  0 - non-standard, user supplie formula
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'ATOMIC.FOR'
      INCLUDE 'MODELQ.FOR'
      PARAMETER (CC0   =   5.465D-11,
     *           CEX1  = -30.20581,
     *           CEX2  =   3.8608704,
     *           CEX3  = 305.63574,
     *           CI1   =   0.3,
     *           CI2   =   0.435,
     *           CA1   =   5.444416D7,
     *           CA2   =  -2.8185937D4,
     *           CA3   =  19.987261,
     *           CA4   =  -5.8906298D-5,
     *           CB1   =   1.3935312D3,
     *           CB2   =  -1.6805859D2,
     *           CB3   =  -2.539D3,
     *           CC1   =   2.0684609D3,
     *           CC2   =  -3.341582D2,
     *           CC3   =  -7.6440625D3,
     *           CD1   =   3.2174844D3,
     *           CD2   =  -5.5882422D2,
     *           CD3   =  -6.86325D3,
     *           CE1   =   5.759125D3,
     *           CE2   =  81.75,
     *           CE3   =  -1.5163D3,
     *           CF1   =   1.461475D4,
     *           CF2   = 393.4,
     *           CF3   =  -4.8284D3,
     *           ALF0  =   1.8,
     *           ALF1  =   0.4,
     *           BET0  =   3.0,
     *           BET1  =   1.2,
     *           O148  =   0.148,
     *           CHMI  =   5.59D-15)
      PARAMETER (EXPIA1=-0.57721566,EXPIA2=0.99999193,
     *           EXPIA3=-0.24991055,EXPIA4=0.05519968,
     *           EXPIA5=-0.00976004,EXPIA6=0.00107857,
     *           EXPIB1=0.2677734343,EXPIB2=8.6347608925,
     *           EXPIB3=18.059016973,EXPIB4=8.5733287401,
     *           EXPIC1=3.9584969228,EXPIC2=21.0996530827,
     *           EXPIC3=25.6329561486,EXPIC4=9.5733223454)
      DIMENSION COL(MTRANS),A(6,10)
      DIMENSION CCOOL(4,14,15),CHOT(4,14,15),XTT(4)
      DATA ((A(I,J),J=1,10),I=1,6) /
     * -86.7633398, 2632.8369 , 7478.9556 ,-4202.8442 ,-47995.930    ,
     * -120942.89 ,-202300.81 ,-261373.03 ,-266337.91 ,-192293.20    ,
     * 100.919188 ,-2738.7485 ,-8495.4590 , 1937.3763 , 45825.371    ,
     *  122209.39 , 211928.67 , 285044.75 , 309455.47 , 258802.22    ,
     * -45.7813807, 1121.3976 , 3794.6826 , 340.35764 ,-16617.055    ,
     * -47390.313 ,-84973.688 ,-117833.95 ,-133243.61 ,-120363.95    ,
     * 10.1978559 ,-224.30670 ,-822.83636 ,-290.10489 , 2905.7393    ,
     *  8944.6025 , 16556.992 , 23544.543 , 27419.742 , 26002.143    ,
     * -1.11223557, 21.923729 , 86.619110 , 48.840523 ,-246.99014    ,
     * -828.41028 ,-1581.2722 ,-2297.9321 ,-2738.1743 ,-2686.4087    ,
     *  .0474198818,-.83974838,-3.5534720 ,-2.6097214 , 8.1972208    ,
     *  30.267115 , 59.521984 , 88.178680 , 107.05288 , 107.73775    /
C
      DATA ((CCOOL(I, 1, K),I=1,4),K=1,15)/   4*0.,
     &        5.742D-01, 1.818D-05,-1.093D-10, 8.687D-16,
     &        1.934D-01,-4.698D-07, 8.352D-11,-5.576D-16,
     &        6.323D-03, 2.237D-06,-1.620D-11, 8.955D-17,
     &        2.035D-02, 6.076D-07,-2.175D-13,-2.495D-18,
     &        1.136D-02, 3.428D-07,-1.467D-13,-1.300D-18,
     &        6.999D-03, 2.126D-07,-9.963D-14,-7.672D-19,
     &        4.624D-03, 1.410D-07,-6.969D-14,-4.927D-19,
     &        3.217D-03, 9.836D-08,-5.031D-14,-3.361D-19,
     &        2.329D-03, 7.135D-08,-3.737D-14,-2.400D-19,
     &        1.741D-03, 5.342D-08,-2.845D-14,-1.775D-19,
     &        1.336D-03, 4.103D-08,-2.213D-14,-1.351D-19,
     &        1.048D-03, 3.220D-08,-1.754D-14,-1.053D-19,
     &        8.369D-04, 2.574D-08,-1.413D-14,-8.368D-20,
     &        6.791D-04, 2.090D-08,-1.154D-14,-6.763D-20/
      DATA ((CCOOL(I, 2, K),I=1,4),K=1,15)/   8*0.,
     &        2.253D+01, 9.350D-04, 1.215D-08,-9.969D-14,
     &        7.816D-01, 5.414D-04,-1.827D-09, 5.140D-17,
     &        1.459D+00, 2.858D-04,-2.207D-09, 9.028D-15,
     &        7.172D-01, 1.440D-04,-1.139D-09, 4.755D-15,
     &        4.107D-01, 8.360D-05,-6.699D-10, 2.823D-15,
     &        2.591D-01, 5.319D-05,-4.293D-10, 1.819D-15,
     &        1.747D-01, 3.608D-05,-2.925D-10, 1.243D-15,
     &        1.237D-01, 2.567D-05,-2.087D-10, 8.891D-16,
     &        9.097D-02, 1.893D-05,-1.539D-10, 6.585D-16,
     &        6.896D-02, 1.438D-05,-1.174D-10, 5.017D-16,
     &        5.356D-02, 1.119D-05,-9.150D-11, 3.913D-16,
     &        4.247D-02, 8.887D-06,-7.272D-11, 3.112D-16,
     &        3.425D-02, 7.176D-06,-5.877D-11, 2.516D-16/
      DATA ((CCOOL(I, 3, K),I=1,4),K=1,15)/  12*0.,
     &       -1.290D+01, 2.059D-02, 5.461D-08,-9.082D-13,
     &        3.562D+02, 7.337D-03,-9.622D-08, 5.596D-13,
     &        5.744D+00, 3.570D-03,-3.259D-08, 1.452D-13,
     &        2.968D+00, 1.813D-03,-1.703D-08, 7.744D-14,
     &        1.756D+00, 1.065D-03,-1.016D-08, 4.667D-14,
     &        1.135D+00, 6.865D-04,-6.601D-09, 3.053D-14,
     &        7.802D-01, 4.713D-04,-4.558D-09, 2.116D-14,
     &        5.615D-01, 3.390D-04,-3.292D-09, 1.532D-14,
     &        4.189D-01, 2.528D-04,-2.461D-09, 1.148D-14,
     &        3.213D-01, 1.939D-04,-1.891D-09, 8.833D-15,
     &        2.523D-01, 1.522D-04,-1.487D-09, 6.953D-15,
     &        2.018D-01, 1.218D-04,-1.192D-09, 5.576D-15/
      DATA ((CCOOL(I, 4, K),I=1,4),K=1,15)/  16*0.,
     &        4.139D+03, 4.645D-01,-7.097D-06, 4.388D-11,
     &        1.794D+03, 4.443D-02,-6.484D-07, 3.936D-12,
     &        1.536D+01, 2.042D-02,-2.065D-07, 9.734D-13,
     &        8.730D+00, 1.033D-02,-1.074D-07, 5.161D-13,
     &        5.434D+00, 6.084D-03,-6.423D-08, 3.116D-13,
     &        3.628D+00, 3.938D-03,-4.196D-08, 2.048D-13,
     &        2.554D+00, 2.718D-03,-2.914D-08, 1.428D-13,
     &        1.873D+00, 1.967D-03,-2.119D-08, 1.041D-13,
     &        1.418D+00, 1.476D-03,-1.594D-08, 7.843D-14,
     &        1.102D+00, 1.138D-03,-1.232D-08, 6.075D-14,
     &        8.744D-01, 8.987D-04,-9.746D-09, 4.809D-14/
      DATA ((CCOOL(I, 5, K),I=1,4),K=1,15)/  20*0.,
     &       -9.122D+02, 1.260D+00,-1.070D-05, 4.290D-11,
     &        3.959D+01, 2.108D-01,-2.162D-06, 1.020D-11,
     &        3.691D+01, 7.806D-02,-8.485D-07, 4.166D-12,
     &        2.352D+01, 3.911D-02,-4.365D-07, 2.179D-12,
     &        1.542D+01, 2.296D-02,-2.601D-07, 1.310D-12,
     &        1.062D+01, 1.487D-02,-1.699D-07, 8.608D-13,
     &        7.642D+00, 1.029D-02,-1.183D-07, 6.014D-13,
     &        5.695D+00, 7.464D-03,-8.621D-08, 4.394D-13,
     &        4.368D+00, 5.617D-03,-6.508D-08, 3.323D-13,
     &        3.430D+00, 4.348D-03,-5.051D-08, 2.583D-13/
      DATA ((CCOOL(I, 6, K),I=1,4),K=1,15)/  24*0.,
     &       -3.431D+03, 4.116D+00,-3.853D-05, 1.679D-10,
     &        4.397D+01, 6.434D-01,-7.008D-06, 3.431D-11,
     &        8.927D+01, 2.325D-01,-2.667D-06, 1.350D-11,
     &        6.153D+01, 1.152D-01,-1.354D-06, 6.957D-12,
     &        4.165D+01, 6.729D-02,-8.024D-07, 4.156D-12,
     &        2.923D+01, 4.349D-02,-5.232D-07, 2.724D-12,
     &        2.130D+01, 3.008D-02,-3.641D-07, 1.902D-12,
     &        1.603D+01, 2.185D-02,-2.656D-07, 1.391D-12,
     &        1.239D+01, 1.647D-02,-2.008D-07, 1.054D-12/
      DATA ((CCOOL(I, 7, K),I=1,4),K=1,15)/  28*0.,
     &       -9.280D+03, 1.116D+01,-1.122D-04, 5.167D-10,
     &        6.658D+01, 1.651D+00,-1.884D-05, 9.487D-11,
     &        2.172D+02, 5.833D-01,-6.977D-06, 3.615D-11,
     &        1.535D+02, 2.858D-01,-3.499D-06, 1.838D-11,
     &        1.049D+02, 1.660D-01,-2.060D-06, 1.090D-11,
     &        7.412D+01, 1.070D-01,-1.339D-06, 7.118D-12,
     &        5.428D+01, 7.389D-02,-9.304D-07, 4.963D-12,
     &        4.103D+01, 5.366D-02,-6.786D-07, 3.629D-12/
      DATA ((CCOOL(I, 8, K),I=1,4),K=1,15)/  32*0.,
     &       -2.069D+04, 2.637D+01,-2.802D-04, 1.342D-09,
     &        2.055D+02, 3.731D+00,-4.420D-05, 2.276D-10,
     &        5.123D+02, 1.292D+00,-1.599D-05, 8.442D-11,
     &        3.578D+02, 6.265D-01,-7.922D-06, 4.235D-11,
     &        2.438D+02, 3.616D-01,-4.633D-06, 2.494D-11,
     &        1.721D+02, 2.322D-01,-3.000D-06, 1.622D-11,
     &        1.260D+02, 1.601D-01,-2.081D-06, 1.129D-11/
      DATA ((CCOOL(I, 9, K),I=1,4),K=1,15)/  36*0.,
     &       -4.032D+04, 5.614D+01,-6.231D-04, 3.073D-09,
     &        6.989D+02, 7.655D+00,-9.352D-05, 4.903D-10,
     &        1.141D+03, 2.605D+00,-3.313D-05, 1.777D-10,
     &        7.755D+02, 1.250D+00,-1.624D-05, 8.808D-11,
     &        5.234D+02, 7.175D-01,-9.437D-06, 5.153D-11,
     &        3.677D+02, 4.590D-01,-6.087D-06, 3.338D-11/
      DATA ((CCOOL(I,10, K),I=1,4),K=1,15)/  40*0.,
     &       -7.097D+04, 1.101D+02,-1.266D-03, 6.390D-09,
     &        2.018D+03, 1.455D+01,-1.824D-04, 9.708D-10,
     &        2.383D+03, 4.875D+00,-6.348D-05, 3.449D-10,
     &        1.569D+03, 2.319D+00,-3.081D-05, 1.691D-10,
     &        1.046D+03, 1.323D+00,-1.779D-05, 9.830D-11/
      DATA ((CCOOL(I,11, K),I=1,4),K=1,15)/  44*0.,
     &       -1.150D+05, 2.020D+02,-2.392D-03, 1.231D-08,
     &        4.988D+03, 2.601D+01,-3.334D-04, 1.797D-09,
     &        4.675D+03, 8.595D+00,-1.142D-04, 6.273D-10,
     &        2.986D+03, 4.054D+00,-5.491D-05, 3.046D-10/
      DATA ((CCOOL(I,12, K),I=1,4),K=1,15)/  48*0.,
     &       -1.737D+05, 3.511D+02,-4.263D-03, 2.227D-08,
     &        1.094D+04, 4.419D+01,-5.774D-04, 3.146D-09,
     &        8.673D+03, 1.442D+01,-1.950D-04, 1.082D-09/
      DATA ((CCOOL(I,13, K),I=1,4),K=1,15)/  52*0.,
     &       -2.459D+05, 5.829D+02,-7.233D-03, 3.830D-08,
     &        2.191D+04, 7.194D+01,-9.561D-04, 5.259D-09/
      DATA ((CCOOL(I,14, K),I=1,4),K=1,15)/  56*0.,
     &       -3.273D+05, 9.312D+02,-1.178D-02, 6.306D-08/
      DATA ((CHOT(I, 1, K),I=1,4),K=1,15)/   4*0.,
     &        5.856D-01, 1.551D-05,-9.669D-12, 5.716D-19,
     &        1.537D-01, 3.548D-06,-3.224D-12, 7.626D-19,
     &        2.400D-02, 1.419D-06,-2.008D-12, 1.356D-18,
     &        2.002D-02, 6.325D-07,-7.070D-13, 4.096D-19,
     &        1.123D-02, 3.549D-07,-3.998D-13, 2.331D-19,
     &        6.940D-03, 2.194D-07,-2.483D-13, 1.453D-19,
     &        4.593D-03, 1.453D-07,-1.648D-13, 9.667D-20,
     &        3.199D-03, 1.012D-07,-1.150D-13, 6.758D-20,
     &        2.318D-03, 7.334D-08,-8.349D-14, 4.910D-20,
     &        1.727D-03, 5.493D-08,-6.270D-14, 3.695D-20,
     &        1.326D-03, 4.218D-08,-4.821D-14, 2.844D-20,
     &        1.040D-03, 3.310D-08,-3.786D-14, 2.236D-20,
     &        8.305D-04, 2.645D-08,-3.028D-14, 1.790D-20,
     &        6.740D-04, 2.147D-08,-2.460D-14, 1.455D-20/
      DATA ((CHOT(I, 2, K),I=1,4),K=1,15)/   8*0.,
     &        1.710D+01, 1.530D-03,-2.553D-09, 1.924D-15,
     &        8.237D+00, 3.554D-04,-7.566D-10, 6.420D-16,
     &        5.932D+00, 1.301D-04,-2.912D-10, 2.535D-16,
     &        2.987D+00, 6.419D-05,-1.444D-10, 1.260D-16,
     &        1.733D+00, 3.689D-05,-8.324D-11, 7.267D-17,
     &        1.102D+00, 2.334D-05,-5.273D-11, 4.605D-17,
     &        7.472D-01, 1.576D-05,-3.567D-11, 3.116D-17,
     &        5.312D-01, 1.118D-05,-2.532D-11, 2.212D-17,
     &        3.919D-01, 8.232D-06,-1.865D-11, 1.630D-17,
     &        2.977D-01, 6.245D-06,-1.416D-11, 1.237D-17,
     &        2.315D-01, 4.855D-06,-1.101D-11, 9.622D-18,
     &        1.838D-01, 3.851D-06,-8.734D-12, 7.635D-18,
     &        1.484D-01, 3.108D-06,-7.050D-12, 6.164D-18/
      DATA ((CHOT(I, 3, K),I=1,4),K=1,15)/  12*0.,
     &        1.940D+02, 1.949D-02,-3.832D-08, 3.137D-14,
     &        4.729D+02, 1.927D-03,-4.171D-09, 3.628D-15,
     &        6.741D+01, 1.315D-03,-3.145D-09, 2.814D-15,
     &        3.444D+01, 6.477D-04,-1.560D-09, 1.399D-15,
     &        2.031D+01, 3.744D-04,-9.054D-10, 8.130D-16,
     &        1.311D+01, 2.388D-04,-5.789D-10, 5.203D-16,
     &        9.007D+00, 1.629D-04,-3.955D-10, 3.556D-16,
     &        6.484D+00, 1.166D-04,-2.835D-10, 2.550D-16,
     &        4.837D+00, 8.666D-05,-2.108D-10, 1.896D-16,
     &        3.711D+00, 6.631D-05,-1.614D-10, 1.452D-16,
     &        2.914D+00, 5.194D-05,-1.265D-10, 1.138D-16,
     &        2.332D+00, 4.150D-05,-1.011D-10, 9.100D-17/
      DATA ((CHOT(I, 4, K),I=1,4),K=1,15)/  16*0.,
     &        7.204D+03, 1.627D-01,-5.181D-07, 5.605D-13,
     &        2.507D+03, 9.370D-03,-2.091D-08, 1.842D-14,
     &        3.823D+02, 6.480D-03,-1.600D-08, 1.452D-14,
     &        1.950D+02, 3.161D-03,-7.869D-09, 7.157D-15,
     &        1.154D+02, 1.823D-03,-4.561D-09, 4.154D-15,
     &        7.486D+01, 1.165D-03,-2.924D-09, 2.665D-15,
     &        5.178D+01, 7.977D-04,-2.006D-09, 1.830D-15,
     &        3.752D+01, 5.737D-04,-1.444D-09, 1.318D-15,
     &        2.816D+01, 4.283D-04,-1.080D-09, 9.858D-16,
     &        2.174D+01, 3.293D-04,-8.307D-10, 7.587D-16,
     &        1.717D+01, 2.592D-04,-6.544D-10, 5.978D-16/
      DATA ((CHOT(I, 5, K),I=1,4),K=1,15)/  20*0.,
     &        2.166D+04, 4.690D-01,-1.122D-06, 1.008D-12,
     &        3.874D+03, 6.443D-02,-1.596D-07, 1.452D-13,
     &        1.465D+03, 2.207D-02,-5.556D-08, 5.082D-14,
     &        7.410D+02, 1.062D-02,-2.698D-08, 2.476D-14,
     &        4.374D+02, 6.096D-03,-1.556D-08, 1.430D-14,
     &        2.841D+02, 3.889D-03,-9.962D-09, 9.167D-15,
     &        1.969D+02, 2.663D-03,-6.838D-09, 6.297D-15,
     &        1.431D+02, 1.918D-03,-4.935D-09, 4.547D-15,
     &        1.078D+02, 1.436D-03,-3.698D-09, 3.409D-15,
     &        8.353D+01, 1.107D-03,-2.854D-09, 2.632D-15/
      DATA ((CHOT(I, 6, K),I=1,4),K=1,15)/  24*0.,
     &        7.146D+04, 1.379D+00,-3.346D-06, 3.023D-12,
     &        1.187D+04, 1.794D-01,-4.501D-07, 4.118D-13,
     &        4.380D+03, 5.990D-02,-1.527D-07, 1.405D-13,
     &        2.192D+03, 2.846D-02,-7.324D-08, 6.759D-14,
     &        1.288D+03, 1.621D-02,-4.197D-08, 3.881D-14,
     &        8.351D+02, 1.031D-02,-2.678D-08, 2.480D-14,
     &        5.790D+02, 7.050D-03,-1.837D-08, 1.702D-14,
     &        4.213D+02, 5.079D-03,-1.326D-08, 1.229D-14,
     &        3.179D+02, 3.804D-03,-9.944D-09, 9.226D-15/
      DATA ((CHOT(I, 7, K),I=1,4),K=1,15)/  28*0.,
     &        1.954D+05, 3.426D+00,-8.397D-06, 7.624D-12,
     &        3.057D+04, 4.266D-01,-1.080D-06, 9.917D-13,
     &        1.103D+04, 1.392D-01,-3.582D-07, 3.309D-13,
     &        5.458D+03, 6.530D-02,-1.696D-07, 1.572D-13,
     &        3.189D+03, 3.693D-02,-9.653D-08, 8.966D-14,
     &        2.062D+03, 2.338D-02,-6.136D-08, 5.707D-14,
     &        1.429D+03, 1.595D-02,-4.200D-08, 3.910D-14,
     &        1.039D+03, 1.148D-02,-3.029D-08, 2.282D-14/
      DATA ((CHOT(I, 8, K),I=1,4),K=1,15)/  32*0.,
     &        4.651D+05, 7.527D+00,-1.859D-05, 1.694D-11,
     &        6.930D+04, 9.038D-01,-2.302D-06, 2.121D-12,
     &        2.450D+04, 2.891D-01,-7.487D-07, 6.939D-13,
     &        1.200D+04, 1.340D-01,-3.505D-07, 3.260D-13,
     &        6.970D+03, 7.523D-02,-1.981D-07, 1.846D-13,
     &        4.493D+03, 4.741D-02,-1.254D-07, 1.170D-13,
     &        3.106D+03, 3.226D-02,-8.559D-08, 7.997D-14/
      DATA ((CHOT(I, 9, K),I=1,4),K=1,15)/  36*0.,
     &        9.956D+05, 1.506D+01,-3.741D-05, 3.418D-11,
     &        1.425D+05, 1.754D+00,-4.489D-06, 4.146D-12,
     &        4.949D+04, 5.510D-01,-1.435D-06, 1.333D-12,
     &        2.401D+04, 2.526D-01,-6.645D-07, 6.196D-13,
     &        1.386D+04, 1.408D-01,-3.729D-07, 3.485D-13,
     &        8.904D+03, 8.835D-02,-2.350D-07, 2.200D-13/
      DATA ((CHOT(I,10, K),I=1,4),K=1,15)/  40*0.,
     &        1.961D+06, 2.798D+01,-6.982D-05, 6.394D-11,
     &        2.715D+05, 3.175D+00,-8.158D-06, 7.551D-12,
     &        9.279D+04, 9.821D-01,-2.567D-06, 2.390D-12,
     &        4.460D+04, 4.458D-01,-1.178D-06, 1.100D-12,
     &        2.561D+04, 2.468D-01,-6.566D-07, 6.150D-13/
      DATA ((CHOT(I,11, K),I=1,4),K=1,15)/  44*0.,
     &        3.613D+06, 4.898D+01,-1.227D-04, 1.125D-10,
     &        4.861D+05, 5.434D+00,-1.401D-05, 1.299D-11,
     &        1.638D+05, 1.658D+00,-4.348D-06, 4.054D-12,
     &        7.810D+04, 7.456D-01,-1.976D-06, 1.850D-12/
      DATA ((CHOT(I,12, K),I=1,4),K=1,15)/  48*0.,
     &        6.300D+06, 8.163D+01,-2.051D-04, 1.884D-10,
     &        8.271D+05, 8.881D+00,-2.296D-05, 2.131D-11,
     &        2.753D+05, 2.676D+00,-7.037D-06, 6.571D-12/
      DATA ((CHOT(I,13, K),I=1,4),K=1,15)/  52*0.,
     &        1.049D+07, 1.305D+02,-3.288D-04, 3.025D-10,
     &        1.348D+06, 1.396D+01,-3.617D-05, 3.361D-11/
      DATA ((CHOT(I,14, K),I=1,4),K=1,15)/  56*0.,
     &        1.680D+07, 2.016D+02,-5.089D-04, 4.687D-10/
C
C     experimental fit formula for Lyman-alpha
C
      CEXP12(CC,X,U)=CC*EXP(-U)*(CEX1+CEX2*X+CEX3/X/X)
C
      HKT=HK/T
      CT=CC0*SQRT(T)
      TK=HKT/H
      t0=t
      X=LOG10(T)
      X2=X*X
      X3=X*X2
      X4=X2*X2
      X5=X3*X2
      XTT(1)=1.
      XTT(2)=T
      XTT(3)=T*T
      XTT(4)=T*T*T
      SQT=SQRT(T)
      N0HN=NFIRST(IELH)
      N1H=NLAST(IELH)
      NKH=NKA(IATH)
      N0Q=NQUANT(N1H)+1
      N1Q=ICUP(IELH)
      DO 200 II=N0HN,N1H
         I=II-N0HN+1
         IT=ITRA(II,NKH)
         IF(IT.EQ.0) GO TO 100
C
C *************** Collisional ionization
C
c        for high temperature, use XSTAR formulae
C
         if(t0.gt.1.e6) then
            rno=16.
            izc=1
            call irc(i,t0,izc,rno,cs)
            col(it)=cs
          go to 100
         end if
C
         IC=ICOL(IT)
         U0=FR0(IT)*HKT
         IF(IC.LT.0) GO TO 90
         if(ifwop(ii).lt.0) go to 95
         GAM=I*I*I
         IF(I.GT.10) GO TO 80
         GAM=A(1,I)+A(2,I)*X+A(3,I)*X2+A(4,I)*X3
     *      +A(5,I)*X4+A(6,I)*X5
   80    COL(IT)=CT*EXP(-U0)*GAM
         GO TO 100
C
C        non-standard (user supplied) formula
C
   90    CALL CSPEC(II,NKH,IC,OSC0(IT),CPAR(IT),U0,T,COL(IT))
         go to 100
c
c      ionization from the merged state
c
   95    sum1=0.
         sum2=0.
         ehk=eh/tk
         n00q=nquant(n1h-1)+1
         n11q=nlmx
         do 96 img=n00q,n11q
            xi=img
            xii=xi*xi
            sum1=sum1+xii*xii*xi*wnhint(img,id)
            sum2=sum2+xii*wnhint(img,id)*exp(ehk/xii)
   96    continue
         col(it)=ct*sum1/sum2
         go to 200
C
C ***************** Collisional excitation
C
  100    I1=I+1
         XI=I
         VI=XI*XI
         ALF=ALF0-ALF1/VI
         BET=BET0-BET1/XI
         NHL=N1H-N0HN+1
         IF(N1Q.GT.0) NHL=N1Q
         N1HC=N1H
         IF(IFWOP(N1H).LT.0) THEN
            NHL=NLMX
            N1HC=N1H-1
         END IF
         CSUM=0.
         IF(I1.GT.NHL) GO TO 200
         CSCA=8.63D-6/2./VI/SQT
         DO 190 J=I1,NHL
            XJ=J
            VJ=XJ*XJ
            IC=0
            JJ=J+N0HN-1
            IF(JJ.GT.N1HC) GO TO 150
            ICT=ITRA(II,JJ)
            IF(ICT.EQ.0) GO TO 190
            IC=ICOL(ICT)
            U0=FR0(ICT)*HKT
            E=U0/EH/TK
            C1=OSC0(ICT)
            IF(IC) 110,160,120
  110       CALL CSPEC(II,JJ,IC,C1,CPAR(ICT),U0,T,COL(ICT))
            GO TO 190
  120       IF(IC-1) 130,130,140
  130       COL(ICT)=CEXP12(CT,X,U0)
            GO TO 190
  140       COL(ICT)=CEH12(T)
            GO TO 190
C
C      collisional excitations from level I to higher, non-explicit
C      levels are lumped into the collisional ionization rate
C      (the so-called modified collision ionization rate)
C
  150       CONTINUE
            E=UN/VI-UN/VJ
            U0=EH*E*TK
            IF(J.LE.20) C1=OSH(I,J)
            IF(J.GT.20) THEN
               C1=OSH(I,20)*((400.-VI)/20.*XJ/(VJ-VI))**3
            end if
  160       CONTINUE
C
            IF(ICOLHN.EQ.1.AND.J.LE.7) GO TO 250
            IF(ICOLHN.EQ.2.AND.J.LE.15) GO TO 260
C
C      Old standard formula for the collisional excitation rate - used for
C      rates in explicit transitions as well as for evaluation of the
C      modified collisional rate
C
            IF(ICOLHN.EQ.1.AND.J.LE.7) GO TO 250
            IF(ICOLHN.EQ.2.AND.J.LE.15) GO TO 260
            CS=4.*CT*C1/E/E
            EX=EXP(-U0)
            IF(U0.LE.UN) THEN
              E1=-LOG(U0)+EXPIA1+U0*(EXPIA2+U0*(EXPIA3+U0*(EXPIA4+
     *            U0*(EXPIA5+U0*EXPIA6))))
            ELSE
              E1=EXP(-U0)*((EXPIB1+U0*(EXPIB2+U0*(EXPIB3+
     *           U0*(EXPIB4+U0))))/(EXPIC1+U0*(EXPIC2+
     *           U0*(EXPIC3+U0*(EXPIC4+U0)))))/U0
            END IF
            E5=E1
            DO 170 IX=1,4
  170          E5=(EX-U0*E5)/IX
            CS=CS*U0*(E1+O148*U0*E5)
            IF(J-I.NE.1) CS=CS*(BET+TWO*(ALF-BET)/(XJ-XI))
          GO TO 180
C      End of the old standard formula (Mihalas et al 1975)
C
c      Butler new calculations
c
  250       call butler(i,j,t,u0,cs,ierr)
            go to 180
c
C      Giovanardi et al. 1987, AAS, 70, 269
C        Cool: T<=60000K ; Hot: T>60000K
C
  260       IF(T.GT.60000.) GOTO 270
            CS=CCOOL(1,I,J)
          DO ICA=2,4
             CS=CS+CCOOL(ICA,I,J)*XTT(ICA)
          END DO
            GOTO 280
  270       CS=CHOT(1,I,J)
          DO ICA=2,4
             CS=CS+CHOT(ICA,I,J)*XTT(ICA)
          END DO
  280       CS=CSCA*CS*EXP(-U0)
C
  180       IF(JJ.GT.N1HC) THEN
               CSUM=CSUM+CS
             ELSE
               COL(ICT)=CS
            END IF
  190    CONTINUE
         IF(IT.NE.0.AND.N1Q.GT.0) COL(IT)=COL(IT)+CSUM
         ITH=ITRA(II,N1H)
         IF(IFWOP(N1H).LT.0.AND.ITH.GT.0) COL(ITH)=CSUM
  200 CONTINUE
C
C     special standard formula for collisional ionization of H-
C
      IF(IELHM.EQ.0) RETURN
      IT=ITRA(NFIRST(IELHM),N0HN)
      IF(IT.EQ.0) RETURN
      IC=ICOL(IT)
      IF(IC) 220,210,210
  210 COL(IT)=CHMI*T*SQRT(T)
      RETURN
C
C     if desired, non-standard, user supplied, formula for H-
C
  220 U0=ENION(NFIRST(IELHM))*TK
      CALL CSPEC(NFIRST(IELHM),N0HN,IC,OSC0(IT),CPAR(IT),U0,T,COL(IT))
      RETURN
      END
C
C
C     ****************************************************************
C
C

      SUBROUTINE BUTLER (NI,NJ,T,U0,COL,IERR)
C     =======================================
C
C Rate coefficients for collisional excitation of hydrogen
C by electrons. Interpolates in Table 3 of Przybilla & Butler
C (2004, ApJ).
C
C
C Input:
C NI Principal quantum number lower level
C NJ "" upper level
C T Temperature
C U0 =h*nu/K/T
C Output:
C COL collisional rate (cm3 s-1)
C IERR error flat (0=ok, 1=T exceeds table range,
C 2=NI higher than 6 or lower than 1
C NJ higher than 7 or lower than 2)
C
      INCLUDE 'IMPLIC.FOR'

      DIMENSION COLSTR(16,21),TREF(16)

      DATA (TREF(I), I=1,16) /
     * 2.5d3, 5d3, 7.5d3, 1d4, 1.5d4, 2d4, 2.5d4, 3d4, 4d4, 5d4, 6d4,
     * 8d4, 1d5, 1.5d5, 2d5, 2.5d5 /

      DATA ((COLSTR(I,J),J=1,21),I=1,16) /
C J=1,21 corresponds to (NI,NJ)={(1,2),(1,3),...,(1,NL),(2,3),...}
C where NL=7 (higher n covered in Table)
C I=1,16 corresponds to T={2.5e3,5e3,7.5e3,1e4,1.5e4,2e4,2.5e4,3e4,
C 4e4,5e4,6e4,8e4,1e5,1.5e5,2e5,2.5e5}
     * 6.40d-1, 2.20d-1, 9.93d-2, 4.92d-2, 2.97d-2, 5.03d-2, 2.35d+1,
     * 1.07d+1, 5.22d+0, 2.91d+0, 5.25d+0, 1.50d+2, 7.89d+1, 4.13d+1,
     * 7.60d+1, 5.90d+2, 2.94d+2, 4.79d+2, 1.93d+3, 1.95d+3, 6.81d+3,

     * 6.98d-1, 2.40d-1, 1.02d-1, 5.84d-2, 4.66d-2, 6.72d-2, 2.78d+1,
     * 1.15d+1, 5.90d+0, 4.53d+0, 7.26d+0, 1.90d+2, 9.01d+1, 6.11d+1,
     * 1.07d+2, 8.17d+2, 4.21d+2, 7.06d+2, 2.91d+3, 3.24d+3, 1.17d+4,

     * 7.57d-1, 2.50d-1, 1.10d-1, 7.17d-2, 6.28d-2, 7.86d-2, 3.09d+1,
     * 1.23d+1, 6.96d+0, 6.06d+0, 8.47d+0, 2.28d+2, 1.07d+2, 8.21d+1,
     * 1.25d+2, 1.07d+3, 5.78d+2, 8.56d+2, 4.00d+3, 4.20d+3, 1.50d+4,

     * 8.09d-1, 2.61d-1, 1.22d-1, 8.58d-2, 7.68d-2, 8.74d-2, 3.38d+1,
     * 1.34d+1, 8.15d+0, 7.32d+0, 9.27d+0, 2.70d+2, 1.26d+2, 1.01d+2,
     * 1.37d+2, 1.35d+3, 7.36d+2, 9.66d+2, 5.04d+3, 4.95d+3, 1.73d+4,

     * 8.97d-1, 2.88d-1, 1.51d-1, 1.12d-1, 9.82d-2, 1.00d-1, 4.01d+1,
     * 1.62d+1, 1.04d+1, 9.17d+0, 1.03d+1, 3.64d+2, 1.66d+2, 1.31d+2,
     * 1.52d+2, 1.93d+3, 1.02d+3, 1.11d+3, 6.81d+3, 6.02d+3, 2.03d+4,

     * 9.78d-1, 3.22d-1, 1.80d-1, 1.33d-1, 1.14d-1, 1.10d-1, 4.71d+1,
     * 1.90d+1, 1.23d+1, 1.05d+1, 1.08d+1, 4.66d+2, 2.03d+2, 1.54d+2,
     * 1.61d+2, 2.47d+3, 1.26d+3, 1.21d+3, 8.20d+3, 6.76d+3, 2.21d+4,

     * 1.06d+0, 3.59d-1, 2.06d-1, 1.50d-1, 1.25d-1, 1.16d-1, 5.45d+1,
     * 2.18d+1, 1.39d+1, 1.14d+1, 1.12d+1, 5.70d+2, 2.37d+2, 1.72d+2,
     * 1.68d+2, 2.96d+3, 1.46d+3, 1.29d+3, 9.29d+3, 7.29d+3, 2.33d+4,

     * 1.15d+0, 3.96d-1, 2.28d-1, 1.64d-1, 1.33d-1, 1.21d-1, 6.20d+1,
     * 2.44d+1, 1.52d+1, 1.21d+1, 1.14d+1, 6.72d+2, 2.68d+2, 1.86d+2,
     * 1.72d+2, 3.40d+3, 1.64d+3, 1.34d+3, 1.02d+4, 7.70d+3, 2.41d+4,

     * 1.32d+0, 4.64d-1, 2.66d-1, 1.85d-1, 1.45d-1, 1.27d-1, 7.71d+1,
     * 2.89d+1, 1.74d+1, 1.31d+1, 1.17d+1, 8.66d+2, 3.19d+2, 2.08d+2,
     * 1.78d+2, 4.14d+3, 1.92d+3, 1.41d+3, 1.15d+4, 8.26d+3, 2.52d+4,

     * 1.51d+0, 5.26d-1, 2.95d-1, 2.01d-1, 1.53d-1, 1.31d-1, 9.14d+1,
     * 3.27d+1, 1.90d+1, 1.38d+1, 1.18d+1, 1.04d+3, 3.62d+2, 2.24d+2,
     * 1.81d+2, 4.75d+3, 2.15d+3, 1.46d+3, 1.26d+4, 8.63d+3, 2.60d+4,

     * 1.68d+0, 5.79d-1, 3.18d-1, 2.12d-1, 1.58d-1, 1.34d-1, 1.05d+2,
     * 3.60d+1, 2.03d+1, 1.44d+1, 1.19d+1, 1.19d+3, 3.98d+2, 2.36d+2,
     * 1.83d+2, 5.25d+3, 2.33d+3, 1.50d+3, 1.34d+4, 8.88d+3, 2.69d+4,

     * 2.02d+0, 6.70d-1, 3.55d-1, 2.29d-1, 1.65d-1, 1.35d-1, 1.29d+2,
     * 4.14d+1, 2.23d+1, 1.51d+1, 1.19d+1, 1.46d+3, 4.53d+2, 2.53d+2,
     * 1.85d+2, 6.08d+3, 2.61d+3, 1.55d+3, 1.49d+4, 9.21d+3, 2.90d+4,

     * 2.33d+0, 7.43d-1, 3.83d-1, 2.39d-1, 1.70d-1, 1.37d-1, 1.51d+2,
     * 4.56d+1, 2.37d+1, 1.56d+1, 1.20d+1, 1.67d+3, 4.95d+2, 2.65d+2,
     * 1.86d+2, 6.76d+3, 2.81d+3, 1.57d+3, 1.63d+4, 9.43d+3, 3.17d+4,

     * 2.97d+0, 8.80d-1, 4.30d-1, 2.59d-1, 1.77d-1, 1.39d-1, 1.93d+2,
     * 5.31d+1, 2.61d+1, 1.63d+1, 1.19d+1, 2.08d+3, 5.68d+2, 2.83d+2,
     * 1.87d+2, 8.08d+3, 3.15d+3, 1.61d+3, 1.97d+4, 9.78d+3, 3.94d+4,

     * 3.50d+0, 9.79d-1, 4.63d-1, 2.71d-1, 1.82d-1, 1.39d-1, 2.26d+2,
     * 5.83d+1, 2.78d+1, 1.68d+1, 1.19d+1, 2.39d+3, 6.16d+2, 2.94d+2,
     * 1.86d+2, 9.13d+3, 3.36d+3, 1.62d+3, 2.27d+4, 1.00d+4, 4.73d+4,

     * 3.95d+0, 1.06d+0, 4.88d-1, 2.81d-1, 1.85d-1, 1.40d-1, 2.52d+2,
     * 6.23d+1, 2.89d+1, 1.71d+1, 1.19d+1, 2.62d+3, 6.51d+2, 3.02d+2,
     * 1.87d+2, 1.00d+4, 3.51d+3, 1.63d+3, 2.54d+4, 1.02d+4, 5.50d+4 /

C do i=1,21
C write(*,'(1X,16(E10.4,1X))')(colstr(J,i), J=1,16)
C enddo

      NL=7
      IERR=0
      COL=0.0d0

      IF (T.LT.2.5d3.OR.T.GE.2.5d5) IERR=1
      IF (NI.LT.1.OR.NI.GT.NL-1.OR.NJ.LT.2.OR.NJ.GT.NL) IERR=2

      IF (IERR.EQ.0) THEN
      J=0
      DO I=1,NI-1
      J=J+(NL-I)
      ENDDO
      DO K=I+1,NJ
      J=J+1
      ENDDO
C write(*,*)NI,NJ,J

C find out nearest points in TREF

      ILOW=1
      DO WHILE (T.GE.TREF(ILOW+1))
      ILOW=ILOW+1
      ENDDO

      IHIG=16
      DO WHILE (T.LT.TREF(IHIG-1))
      IHIG=IHIG-1
      ENDDO

      IF (IHIG.EQ.ILOW) IHIG=IHIG+1

C write(*,*)'ilow=',ILOW,' ihig=',IHIG

C interpolate linearly (log-log) the collision strength

      SL=LOG10(COLSTR(IHIG,J))-LOG10(COLSTR(ILOW,J))
      SL=SL/(LOG10(TREF(IHIG))-LOG10(TREF(ILOW)))
      OR=LOG10(COLSTR(IHIG,J))-SL*LOG10(TREF(IHIG))
      COL=LOG10(T)*SL+OR
      COL=10.**COL

C write(*,*)col,colstr(ilow,j),colstr(ihig,j)

C derive the rate

      COL=8.631d-6/(2.d0*NI**2)/SQRT(T)*EXP(-U0)*COL

      ENDIF

      RETURN
      END


C
C
C     ****************************************************************
C
C
      SUBROUTINE COLHE(T,COL)
C     =======================
C
C     Helium (both neutral and ionized) collision rates
C
C     Meaning of ICOL: for all kinds of transitions and for both HeI
C     and HeII:
C      ICOL =  0  - approximate expressions taken from Mihalas, Heasley
C                   and Auer, NCAR-TN/STR-104 (1975) - for He I and II
C
C     New expression for He II collisional ionization from Klaus Werner
C       (ICOL >= 1)  (also originally from Mihalas)
C
C     For He I bound-bound transitions, the following standard
C     possibilities are also available:
C
C      ICOL =  1, 2, or 3  - much more accurate Storey's rates,
C                   subroutine written by D.G.Hummer (COLLHE).
C                   This procedure can be used only for transitions
C                   between states with n = 1, 2, 3, 4.
C      ICOL =  1  - means that a given transition is a transition
C                   between non-averaged (l,s) states. In this case,
C                   labeling of the He I energy levels must agree
C                   with that given in subroutine COLLHE, ie. states
C                   have to be labeled sequentially in order of
C                   increasing frequency.
C      ICOL =  2  - means that a given transition is a transition between
C                   a non-averaged (l,s) lower state and averaged upper
C                   state.
C      ICOL =  3  - means that a given transition is a transition between
C                   two averaged states.
C      Note:
C     The program allows only two standard possibilities of
C     constructing averaged levels of He I:
C     i)  all states within given principal quantum number n (>1) are
C         lumped together
C     ii) all siglet states for given n, and all triplet states for
C         given n are lumped together separately (there are thus two
C         explicit levels for a given n).
C     If the user wants to use another averaging, he had to take care
C     of appropriate averaged collisional rates himself (by updating
C     subroutine CSPEC)
C
C
C      ICOL <  0  - non-standard, user supplied formula
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'ATOMIC.FOR'
      DIMENSION COL(MTRANS)
      DIMENSION FHE1(16),G0(3),G1(3),G2(3),G3(3),A(6,10)
      PARAMETER (EXPIA1=-0.57721566,EXPIA2=0.99999193,
     *           EXPIA3=-0.24991055,EXPIA4=0.05519968,
     *           EXPIA5=-0.00976004,EXPIA6=0.00107857,
     *           EXPIB1=0.2677734343,EXPIB2=8.6347608925,
     *           EXPIB3=18.059016973,EXPIB4=8.5733287401,
     *           EXPIC1=3.9584969228,EXPIC2=21.0996530827,
     *           EXPIC3=25.6329561486,EXPIC4=9.5733223454)
      DATA FHE1/0.,2.75D-1,7.29D-2,2.96D-2,1.48  D-2,8.5D-3,5.3D-3,
     *          3.5D-3,2.5D-3,1.8D-3,1.5D-3,1.2D-3,9.4D-4,7.5D-4,
     *          6.1D-4,5.3D-4/
      DATA G0/ 7.3399521D-2,  1.7252867,      8.6335087    /,
     *     G1/-1.4592763D-7,  2.0944117D-6,   2.7575544D-5 /,
     *     G2/ 7.6621299D5,   5.4254879D6,    6.6395519D6  /,
     *     G3/ 2.3775439D2,   2.2177891D3,    5.20725D3    /
      DATA ((A(I,J),J=1,10),I=1,6) /
     * -8.5931587    , 85.014091    , 923.64099, 2018.6470, 1551.5061  ,
     * -2327.4819    ,-10701.481    ,-27619.789,-41099.602,-61599.023  ,
     *  9.3868790    ,-78.834488    ,-969.18451,-2243.1768,-2059.9768  ,
     *  1546.7107    , 9834.3447    , 27067.436, 41421.254, 63594.133  ,
     * -4.0027571    , 28.360615    , 401.23965, 983.83374, 1051.4103  ,
     * -204.82320    ,-3335.4211    ,-10100.119,-15863.257,-24949.125  ,
     * 0.83941799    ,-4.7963457    ,-81.122566,-209.86169,-251.30855  ,
     * -43.175175    , 530.37292    , 1826.1049, 2941.6460, 4740.8364  ,
     * -.86396709E-01,0.37385577    , 8.0078983, 21.757591, 28.375637  ,
     *  11.890312    ,-39.536087    ,-161.52513,-266.86011,-440.88257  ,
     * 0.34853835E-02,-.10401310E-01,-.30957383,-.87988985,-1.2254572  ,
     * -.72724497    , 1.0879648    , 5.6239786, 9.5323009, 16.150818  /
      SAVE FHE1,G0,G1,G2,G3
C
      HKT=HK/T
      TK=HKT/H
      SRT=SQRT(T)
      t0=t
      CT=5.465D-11*SRT
      CT1=5.4499487/T/SRT
C
C     --------------
C     Neutral helium
C     --------------
C
      IF(IELHE1.EQ.0) GO TO 60
      ICALL=0
      N0I=NFIRST(IELHE1)
      N1I=NLAST(IELHE1)
      NKI=NNEXT(IELHE1)
      N0Q=NQUANT(NLAST(IELHE1))+1
      N1Q=ICUP(IELHE1)
      DO 50 II=N0I,N1I
         IT=ITRA(II,NKI)
         IF(IT.EQ.0) GO TO 10
C
C ******** Collisional ionization
C
         IC=ICOL(IT)
         C1=OSC0(IT)
         C2=CPAR(IT)
         U0=ENION(II)*TK
         IF(IC.GE.0) THEN
            U1=U0+0.27
            U2=(U0+3.43)/(U0+1.43)**3
            IF(U0.LE.UN) THEN
         EXPIU0=-LOG(U0)+EXPIA1+U0*(EXPIA2+U0*(EXPIA3+U0*(EXPIA4+
     *            U0*(EXPIA5+U0*EXPIA6))))
            ELSE
         EXPIU0=EXP(-U0)*((EXPIB1+U0*(EXPIB2+U0*(EXPIB3+
     *           U0*(EXPIB4+U0))))/(EXPIC1+U0*(EXPIC2+
     *           U0*(EXPIC3+U0*(EXPIC4+U0)))))/U0
            END IF
            IF(U1.LE.UN) THEN
         EXPIU1=-LOG(U1)+EXPIA1+U1*(EXPIA2+U1*(EXPIA3+U1*(EXPIA4+
     *            U1*(EXPIA5+U1*EXPIA6))))
            ELSE
         EXPIU1=EXP(-U1)*((EXPIB1+U1*(EXPIB2+U1*(EXPIB3+
     *           U1*(EXPIB4+U1))))/(EXPIC1+U1*(EXPIC2+
     *           U1*(EXPIC3+U1*(EXPIC4+U1)))))/U1
            END IF
            COL(IT)=CT*C1*U0*(EXPIU0-U0*(0.728*EXPIU1/U1+
     *              0.189*EXP(-U0)*U2))
          ELSE
            CALL CSPEC(II,NKI,IC,C1,C2,U0,T,COL(IT))
         END IF
   10    IF(II.GE.N1I) GO TO 30
C
C ********* Collisional excitation
C
         DO 20 JJ=II+1,N1I
            ICT=ITRA(II,JJ)
            IF(ICT.EQ.0) GO TO 20
            IC=ICOL(ICT)
            C1=OSC0(ICT)
            C2=CPAR(ICT)
            U0=FR0(ICT)*HKT
            IF(IC.EQ.0) THEN
C
C   *** ICOL = 0    Formula used by Mihalas, Heasley, and Auer
C
               IF(U0.LE.UN) THEN
            EX=-LOG(U0)+EXPIA1+U0*(EXPIA2+U0*(EXPIA3+U0*(EXPIA4+
     *            U0*(EXPIA5+U0*EXPIA6))))
               ELSE
            EX=EXP(-U0)*((EXPIB1+U0*(EXPIB2+U0*(EXPIB3+
     *           U0*(EXPIB4+U0))))/(EXPIC1+U0*(EXPIC2+
     *           U0*(EXPIC3+U0*(EXPIC4+U0)))))/U0
               END IF
               IF(II.EQ.N0I) THEN
C
C          excitation from the ground state
C
                  COL(ICT)=CT1*EX/U0*C1
                ELSE
C
C          transitions between excited states
C
                  U1=U0+0.2
                  IF(U1.LE.UN) THEN
          EXPIU1=-LOG(U1)+EXPIA1+U1*(EXPIA2+U1*(EXPIA3+U1*(EXPIA4+
     *            U1*(EXPIA5+U1*EXPIA6))))
                  ELSE
          EXPIU1=EXP(-U1)*((EXPIB1+U1*(EXPIB2+U1*(EXPIB3+
     *           U1*(EXPIB4+U1))))/(EXPIC1+U1*(EXPIC2+
     *           U1*(EXPIC3+U1*(EXPIC4+U1)))))/U1
                  END IF
                  COL(ICT)=CT1/U0*(EX-U0/U1*0.81873*EXPIU1)*C1
               END IF
             ELSE IF(IC.EQ.1) THEN
C
C   *** ICOL = 1   Storey - Hummer collisional rates between
C                  non-averaged states
C                  (Note: procedure COLLHE, which calculates all rates,
C                  is called only once)
C
               IF(ICALL.EQ.0) CALL COLLHE(T,COLHE1)
               ICALL=1
               COL(ICT)=COLHE1(II-N0I+1,JJ-N0I+1)
             ELSE IF(IC.EQ.2.OR.IC.EQ.3) THEN
C
C   *** ICOL = 2 or 3  Storey - Hummer collisional rates between
C                      averaged states
C
               IF(ICALL.EQ.0) CALL COLLHE(T,COLHE1)
               ICALL=1
               COL(ICT)=CHEAV(II,JJ,IC)
             ELSE IF(IC.LT.0) THEN
C
C       Non-standard, user supplied formula
C
               CALL CSPEC(II,JJ,IC,C1,CPAR(ICT),U0,T,COL(ICT))
            END IF
   20    CONTINUE
C
C      collisional excitations from level II to higher, non-explicit
C      levels are lumped into the collisional ionization rate
C      (the so-called modified collision ionization rate);
C      the individual rates are calculated by expressions used by
C      Mihalas, Heasley, and Auer
C
   30    IF(N1Q.EQ.0.OR.IT.EQ.0) GO TO 50
         I=NQUANT(II)
         REL=G(II)/2./I/I
         DO 40 J=N0Q,N1Q
            XJ=J
            U0=(ENION(II)-EH/XJ/XJ)*TK
            IF(I.EQ.1) THEN
               GAM=0.
               C1=FHE1(J)
             ELSE
               C1=OSH(I,J)*REL
               U1=U0+0.2
               IF(U1.LE.UN) THEN
          EXPIU1=-LOG(U1)+EXPIA1+U1*(EXPIA2+U1*(EXPIA3+U1*(EXPIA4+
     *            U1*(EXPIA5+U1*EXPIA6))))
               ELSE
          EXPIU1=EXP(-U1)*((EXPIB1+U1*(EXPIB2+U1*(EXPIB3+
     *           U1*(EXPIB4+U1))))/(EXPIC1+U1*(EXPIC2+
     *           U1*(EXPIC3+U1*(EXPIC4+U1)))))/U1
               END IF
               GAM=U0/U1*0.81873*EXPIU1
            END IF
            IF(U0.LE.UN) THEN
          EXPIU0=-LOG(U0)+EXPIA1+U0*(EXPIA2+U0*(EXPIA3+U0*(EXPIA4+
     *            U0*(EXPIA5+U0*EXPIA6))))
            ELSE
          EXPIU0=EXP(-U0)*((EXPIB1+U0*(EXPIB2+U0*(EXPIB3+
     *           U0*(EXPIB4+U0))))/(EXPIC1+U0*(EXPIC2+
     *           U0*(EXPIC3+U0*(EXPIC4+U0)))))/U0
            END IF
            COL(IT)=COL(IT)+CT1/U0*C1*(EXPIU0-GAM)
   40    CONTINUE
   50 CONTINUE
C
C     --------------
C     Ionized helium
C     --------------
C
   60 IF(IELHE2.EQ.0) RETURN
      N0I=NFIRST(IELHE2)
      N1I=NLAST(IELHE2)
      NKI=NNEXT(IELHE2)
      N0Q=NQUANT(NLAST(IELHE2))+1
      N1Q=ICUP(IELHE2)
      X=LOG10(T)
      X2=X*X
      X3=X2*X
      X4=X3*X
      X5=X4*X
      CT2=3.7036489/T/SRT
C
      DO 200 II=N0I,N1I
         I=II-N0I+1
         IT=ITRA(II,NKI)
         IF(IT.EQ.0) GO TO 100
C
C ********* Collisional ionization
C
c        for high temperature, use XSTAR formulae
C
         if(t0.gt.1.e5) then
            rno=16.
            izc=2
            call irc(i,t0,izc,rno,cs)
          col(it)=cs
          go to 100
       end if
C
         IC=ICOL(IT)
       U0=FR0(IT)*HKT
       IF(IC.EQ.0) THEN
            IF(I.LE.3) THEN
               GAM=G0(I)-G1(I)*T+(G2(I)/T-G3(I))/T
             ELSE IF(I.EQ.4) THEN
               GAM=-95.23828+(62.656249-8.1454078*X)*X
           ELSE IF(I.EQ.5) THEN
               GAM=472.99219-74.144287*X-1869.6562/X2
           ELSE IF(I.EQ.6) THEN
               GAM=825.17186-134.23096*X-2739.4375/X2
           ELSE IF(I.EQ.7) THEN
               GAM=1181.3516-200.71191*X-2810.7812/X2
           ELSE IF(I.EQ.8) THEN
               GAM=1440.1016-259.75781*X-1283.5625/X2
           ELSE IF(I.EQ.9) THEN
               GAM=2492.1250-624.84375*X+30.101562*X2
           ELSE IF(I.EQ.10) THEN
               GAM=4663.3129-1390.1250*X+97.671874*X2
             ELSE
               GAM=I*I*I
            END IF
            COL(IT)=CT*EXP(-U0)*GAM
        ELSE IF(IC.GE.1) THEN
          GAM=I*I*I
          IF(I.LE.10) GAM=A(1,I)+A(2,I)*X+A(3,I)*X2+
     *            A(4,I)*X3+A(5,I)*X4+A(6,I)*X5
            COL(IT)=CT*EXP(-U0)*GAM
          ELSE
            CALL CSPEC(II,NKI,IC,OSC0(IT),CPAR(IT),U0,T,COL(IT))
         END IF
C
  100    I1=I+1
         XI=I
         VI=XI*XI
         NHL=N1I-N0I+1
         IF(N1Q.GT.0) NHL=N1Q
         IF(I1.GT.NHL) GO TO 200
C
C ********** collisional excitation
C
C     both explicit transitions as well as contributions to the
C     modified collisional ionization rate
C
         DO 150 J=I1,NHL
            JJ=J+N0I-1
            IC=0
            IF(JJ.GT.N1I) GO TO 110
            ICT=ITRA(II,JJ)
            IF(ICT.EQ.0) GO TO 150
            IC=ICOL(ICT)
  110       XJ=J
            VJ=XJ*XJ
            U0=ENION(N0I)*(1./VI-1./VJ)*TK
            IF(J.LE.20) C1=OSH(I,J)
            IF(J.GT.20) C1=OSH(I,20)*(20./XJ)**3
            IF(IC.LT.0) GO TO 120
            GAM=XI-(XI-1.)/(XJ-XI)
            IF(GAM.GT.XJ-XI) GAM=XJ-XI
            IF(I.GT.1) GAM=GAM*1.1
            IF(U0.LE.UN) THEN
          EXPIU0=-LOG(U0)+EXPIA1+U0*(EXPIA2+U0*(EXPIA3+U0*(EXPIA4+
     *            U0*(EXPIA5+U0*EXPIA6))))
            ELSE
          EXPIU0=EXP(-U0)*((EXPIB1+U0*(EXPIB2+U0*(EXPIB3+
     *           U0*(EXPIB4+U0))))/(EXPIC1+U0*(EXPIC2+
     *           U0*(EXPIC3+U0*(EXPIC4+U0)))))/U0
            END IF
            CS=CT2/U0*C1*(0.693*EXP(-U0)+EXPIU0)*GAM
            GO TO 130
  120       CALL CSPEC(II,JJ,IC,C1,CPAR(ICT),U0,T,COL(ICT))
            GO TO 150
  130       IF(JJ.GT.N1I) GO TO 140
            COL(ICT)=CS
            GO TO 150
  140       IF(IT.NE.0) COL(IT)=COL(IT)+CS
  150    CONTINUE
  200 CONTINUE
      RETURN
      END
C
C
C     ****************************************************************
C
C

      FUNCTION CEH12(T)
C     =================
C
C     Special formula for collisional rate in hydrogen Lyman-alpha
C     transition
C     After Crandall et al. Ap.J. 191, 789 (1974)
C
      INCLUDE 'IMPLIC.FOR'
      DIMENSION A(6),B(8)
      PARAMETER (C=-118353.41)
      DATA A/ 2.579997D-10, -1.629166D-10,  7.713069D-11,
     *       -2.668768D-11,  6.642513D-12, -9.422885D-13/
      SAVE A
c
      DO 10 I=1,8
   10    B(I)=0.
      X=LOG10(T)-4.
      DO 20 I=1,6
         J=7-I
         B(J)=2.*X*B(J+1)-B(J+2)+A(J)
   20 CONTINUE
      CEH12=2.4*SQRT(T)*(B(1)-B(3))*EXP(C/T)
      RETURN
      END
C
C
C     ****************************************************************
C
C

      SUBROUTINE CSPEC(I,J,IC,OS,CP,U0,T,CS)
C     ======================================
C
C     Non-standard evaluation of collision rates
C     Basically user-supplied procedure; here is an example
C
C      Van Regemorter's formula following the recommendations of
C         Mihalas (1978, Stellar Atmospheres, 2nd edition)
C       IC=-1 for neutrals
C        IC=-2 for ions
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'ATOMIC.FOR'
      DIMENSION CHE1FB(3,4)
      DATA CHE1FB/ 9.63675,-2.22941,-17.30103,
     *            10.85578,-2.40931,-27.00903,
     *             8.38043,-2.04791,-7.36621,
     *             6.95825,-2.01967,-5.98779/
      PARAMETER (EXPIA1=-0.57721566,EXPIA2=0.99999193,
     *           EXPIA3=-0.24991055,EXPIA4=0.05519968,
     *           EXPIA5=-0.00976004,EXPIA6=0.00107857,
     *           EXPIB1=0.2677734343,EXPIB2=8.6347608925,
     *           EXPIB3=18.059016973,EXPIB4=8.5733287401,
     *           EXPIC1=3.9584969228,EXPIC2=21.0996530827,
     *           EXPIC3=25.6329561486,EXPIC4=9.5733223454)
 
      CS=0.

      IF(IC.GT.-10) THEN
 
               IF(U0.LE.UN) THEN
         EXPIU0=-LOG(U0)+EXPIA1+U0*(EXPIA2+U0*(EXPIA3+U0*(EXPIA4+
     *            U0*(EXPIA5+U0*EXPIA6))))
               ELSE
         EXPIU0=EXP(-U0)*((EXPIB1+U0*(EXPIB2+U0*(EXPIB3+
     *           U0*(EXPIB4+U0))))/(EXPIC1+U0*(EXPIC2+
     *           U0*(EXPIC3+U0*(EXPIC4+U0)))))/U0
               END IF

CCCCCC Neutrals (See Auer & Mihalas 1973)
        IF(IC.EQ.-1) THEN

           IF(U0.LE.14.) THEN
               GG=0.276*EXP(U0)*EXPIU0
           ELSE 
               GG=0.066*(1.+1.5/U0)/SQRT(U0)
           ENDIF

CCCCCC Ions (See Mihalas 1972)
        ELSE IF(IC.EQ.-2) THEN

           GG0=0.276*EXP(U0)*EXPIU0
           GG=CP
           IF(GG0.GT.CP) GG=GG0

        END IF
        T32=T**(-1.5)
        CS=19.7363*T32*EXP(-U0)/U0*GG*OS
 
        RETURN
      END IF
C
      IF(IC.EQ.-11) THEN
        XR=-1.68D0
      CS=2.16*U0**XR/T/SQRT(T)*EXP(-U0)*OS
C
C     Forbidden transitions between n=2 He I sublevels
C     (from Klaus Werner)
C
      ELSE IF(IC.EQ.-12) THEN
        N0I=NFIRST(IELHE1)
      I=I-N0I+1
      J=J-N0I+1
      IFORB=0
      IF(I.EQ.2 .AND. J.EQ.3) IFORB=1
      IF(I.EQ.2 .AND. J.EQ.5) IFORB=2
      IF(I.EQ.3 .AND. J.EQ.4) IFORB=3
      IF(I.EQ.4 .AND. J.EQ.5) IFORB=4
      IF(IFORB.EQ.0) CALL QUIT(' Inconsistent ICOL - CSPEC',iforb,0)
      XT=LOG10(T)
      GAM=CHE1FB(1,IFORB)+CHE1FB(2,IFORB)*XT+CHE1FB(3,IFORB)/XT/XT
      GAM=EXP(2.30258509299405*GAM)
      CS=5.465D-11*SQRT(T)*EXP(-U0)*GAM
 
      END IF
 
      RETURN
      END
C
C
C     ****************************************************************
C
C

      SUBROUTINE COLLHE(TEMP,COLHE1)
C     ==============================
C
C          GENERATES COLLISIONAL RATE COEFFICIENTS AMONG THE 19 STATES OF
C          HELIUM WITH N = 1, 2, 3, AND 4, USING RATES EVALUATED FROM THE
C          CROSS SECTIONS CALCULATED BY BERRINGTON AND KINGSTON (J. PHYS.B. 20,
C          6631(1987)).  COLLISIONAL RATE COEFFICIENTS HAVE BEEN EVALUATED
C          NUMERICALLY BY P.J.STOREY FROM THE UNPUBLISHED COMPUTER OUTPUT
C          FILES OF BERRINGTON AND KINGSTON.
C
C          THE STATES INCLUDED IN THE CALCULATION ARE LABELLED SEQUENTIALLY
C          IN ORDER OF INCREASING ENERGY:
C                           1          1 SING S
C                           2          2 TRIP S
C                           3          2 SING S
C                           4          2 TRIP P
C                           5          2 SING P
C                           6          3 TRIP S
C                           7          3 SING S
C                           8          3 TRIP P
C                           9          3 TRIP D
C                          10          3 SING D
C                          11          3 SING P
C                          12          4 TRIP S
C                          13          4 SING S
C                          14          4 TRIP P
C                          15          4 TRIP D
C                          16          4 SING D
C                          17          4 TRIP F
C                          18          4 SING F
C                          19          4 SING P
C
C          THIS ORDERING DIFFERS SLIGHTLY FROM THAT OF BERRINGTON AND KINGSTON,
C          IN WHICH 15 AND 16, AND 17 AND 18, WERE INTERCHANGED.
C
C          THE INTRINSIC ACCURACY OF TRANSITIONS AMONG STATES WITH N = 1, 2,
C          AND 3 IS EXPECTED TO BE CONSIDERABLY BETTER THAN THOSE WITH N = 4
C          AS DISCUSSED BY BERRINGTON AND KINGSTON.  THE FITTING ACCURACY IS
C          EVERYWHERE BETTER THAN 2%.  THE ENERGIES OF THE LEVELS ARE TAKEN
C          FROM W. C. MARTIN, PHYS. CHEM. REF. DATA., VOL.2, 257 (1973) AND
C          ARE GIVEN IN ELECTRON VOLTS.  (BOLTZMANN'S CONSTANT = 8.62E-5)
C
C          FIRST REVISED VERSION: D.G.HUMMER, MAY 1988, JILA
C          slightly modified by I.H.,  July 1988
C
      INCLUDE 'IMPLIC.FOR'
      PARAMETER (UN=1.D0,
     *           C1=3.849485D0,
     *           C2=8.49485002D-1,
     *           N=19)
      DIMENSION ENER(19),COLHE1(19,19),
     *          NSTART(172),A(929),B(10),STWT(19)
C
      DATA ENER/ 0.0D0,19.8198D0,20.6160D0,20.96432D0,21.2182D0,
     .22.7187D0,22.9206D0,23.00731D0,23.0739D0,23.0743D0,23.0873D0,
     .23.5942D0,23.6738D0,23.7081D0,23.7363D0,23.7366D0,23.7373D0,
     .23.7373D0,23.7423D0/
C
      DATA STWT/1.0D0,3.0D0,1.0D0,9.0D0,3.0D0,3.0D0,1.0D0,9.0D0,
     .1.5D1,5.0D0,3.0D0,3.0D0,1.0D0,9.0D0,1.5D1,5.0D0,2.1D1,
     .7.0D0,3.0D0/
C
      DATA NSTART/
     .  1,  6, 11, 16, 20, 28, 32, 40, 44, 52, 57, 62, 67, 72, 77, 82,
     . 88, 92, 98,104,110,114,120,125,129,135,139,147,151,157,164,170,
     .177,183,190,195,202,208,213,220,225,232,236,243,247,251,260,266,
     .273,278,285,290,300,304,309,316,320,324,329,333,338,343,347,352,
     .357,362,367,372,376,382,386,391,395,401,405,410,414,421,425,431,
     .435,440,445,449,454,459,465,470,475,480,487,491,497,503,508,515,
     .520,525,530,536,542,547,552,559,564,571,576,581,587,592,598,603,
     .608,613,617,623,630,635,642,646,650,655,660,666,671,677,683,689,
     .695,702,707,713,718,723,728,732,737,741,745,750,754,759,765,771,
     .777,782,789,796,801,805,810,815,819,824,831,837,844,850,856,861,
     .868,873,877,882,890,895,905,909,913,920,925,930/
C
      DATA (A(I),I=1,95)/
     . 1.7339D-07, 2.7997D-08,-1.3812D-08, 2.6639D-09, 1.7776D-09,
     . 2.9820D-07, 7.5210D-08,-3.5975D-09, 3.2270D-09, 1.5245D-09,
     . 1.5601D-05, 1.5340D-06,-2.2122D-06,-1.1073D-07, 1.9249D-07,
     . 2.3682D-08, 1.0638D-08, 2.0959D-09, 2.8381D-10, 3.0497D-05,
     . 1.9252D-05, 6.3109D-06, 6.9098D-07,-2.8039D-07,-2.1128D-07,
     .-1.2192D-07,-4.4417D-08, 1.3896D-06, 1.5715D-07,-8.4358D-08,
     .-2.8800D-08, 5.9599D-08, 3.4756D-08, 1.2183D-08, 3.7999D-09,
     . 8.5500D-10,-3.9428D-10,-5.3999D-10,-2.2962D-10, 2.2510D-06,
     . 6.1436D-07,-1.2437D-07,-7.1718D-08, 5.9026D-05, 3.8150D-05,
     . 1.1426D-05, 9.2886D-07,-6.4827D-07,-4.4270D-07,-1.8611D-07,
     .-5.6403D-08, 5.8752D-06, 2.5167D-06, 2.0787D-07,-2.3353D-07,
     .-8.9900D-08, 5.6334D-08, 2.9313D-09, 1.7775D-09, 5.5494D-10,
     .-5.8914D-10, 8.1939D-06, 2.4014D-07, 3.8681D-07, 2.8446D-07,
     .-6.1936D-08, 1.8173D-06,-4.7530D-07,-6.6432D-08, 5.4898D-08,
     .-1.2377D-08, 2.0732D-05, 6.5991D-07, 1.5840D-06, 4.9920D-07,
     .-1.8332D-07, 3.2273D-06,-4.9880D-07,-2.4929D-07, 1.1964D-07,
     .-2.1996D-08, 9.7096D-08, 1.3557D-08, 7.4404D-09, 1.3858D-09,
     .-1.3778D-09,-8.4885D-10, 3.5068D-06,-5.0675D-07,-1.2252D-07,
     . 6.0514D-08, 7.0524D-06, 1.4454D-06, 8.3966D-07, 2.7203D-07/
      DATA (A(I),I=96,190)/
     .-2.3854D-08,-8.6693D-08, 7.1193D-06, 8.3111D-09,-9.3916D-07,
     .-2.7944D-08, 1.7803D-07,-3.9216D-08, 9.2760D-06, 2.4761D-06,
     . 1.0095D-06, 3.4039D-07,-8.6900D-08,-1.1156D-07, 2.5036D-05,
     .-5.7791D-06,-1.8197D-06, 1.0630D-06, 9.8746D-09, 3.1048D-09,
     . 1.2172D-09, 1.8411D-10,-1.5835D-10,-1.4443D-10, 1.9240D-06,
     . 1.5012D-07, 7.9741D-08, 2.9323D-08,-1.2796D-08, 5.0457D-07,
     .-5.5146D-08,-2.7506D-08, 1.4341D-09, 1.3321D-05, 2.6960D-06,
     . 2.8059D-07, 1.7548D-08,-2.3623D-07,-1.4619D-07, 1.5294D-06,
     .-9.4925D-08,-1.1347D-07, 8.1980D-09, 1.3324D-04, 7.8068D-05,
     . 2.9238D-05, 4.2718D-06,-1.6556D-06,-1.4529D-06,-8.6046D-07,
     .-2.1062D-07, 2.8510D-06,-4.7936D-07,-2.4042D-07, 6.2333D-08,
     . 1.3493D-09, 3.5113D-10,-4.7269D-11, 2.4872D-11,-9.7136D-12,
     .-1.4217D-11, 1.5680D-06, 8.9272D-07, 2.8313D-07, 5.2456D-08,
     .-2.3404D-08,-2.7304D-08,-8.4539D-09, 1.7967D-07, 6.2133D-08,
     .-3.2814D-09,-3.7299D-09,-1.9587D-09,-1.6685D-09, 1.2707D-05,
     . 7.5029D-06, 2.8330D-06, 6.7478D-07,-1.5012D-07,-2.3916D-07,
     .-8.9594D-08, 7.6160D-07, 2.0422D-07,-2.5881D-08,-1.7624D-08,
     .-8.3518D-09,-4.1743D-09, 4.6044D-05, 2.2425D-05, 3.3079D-06,
     . 3.6752D-07,-8.4476D-08,-4.8455D-07,-2.5391D-07, 7.0824D-07/
      DATA (A(I),I=191,285)/
     . 1.4917D-07,-1.2005D-07,-1.6983D-08, 1.1545D-08, 7.2866D-04,
     . 3.8907D-04, 9.6365D-05, 2.3153D-05, 1.8462D-07,-9.0627D-06,
     .-5.4332D-06, 1.0458D-08, 1.6430D-09, 5.6453D-10, 2.0276D-10,
     .-1.6501D-10,-1.0656D-10, 5.2198D-07, 4.2898D-08,-1.0332D-08,
     .-5.6754D-09,-7.1960D-09, 3.0390D-06, 1.5385D-06, 6.2110D-07,
     . 1.4111D-07,-3.7415D-08,-4.8395D-08,-1.7219D-08, 2.7227D-06,
     . 1.4381D-07,-5.5030D-08,-2.2192D-10,-2.8583D-08, 1.8417D-05,
     . 9.3860D-06, 3.9999D-06, 1.0424D-06,-2.1337D-07,-3.3271D-07,
     .-1.2684D-07, 4.6063D-06,-6.7185D-07,-2.5830D-07, 2.9976D-08,
     . 7.8012D-05, 3.4921D-05, 8.3012D-06, 1.5643D-06,-4.3892D-07,
     .-9.5318D-07,-4.6534D-07, 1.7584D-05,-2.8817D-06,-1.0081D-06,
     . 1.3259D-07, 1.8561D-05, 2.6602D-06,-1.7243D-06,-3.5533D-07,
     . 2.3315D-08, 1.6239D-08, 7.3739D-09, 1.9570D-09,-2.5667D-10,
     .-5.8101D-10,-2.3947D-10, 4.5821D-12, 4.5024D-11, 3.8796D-07,
     . 1.1239D-07,-2.9033D-08,-5.7237D-09, 2.4186D-09,-2.9670D-09,
     . 1.0887D-06, 4.3737D-07, 8.3753D-08, 3.6771D-08, 1.6545D-09,
     .-1.3849D-08,-5.8855D-09, 2.1028D-06, 4.3725D-07,-1.7621D-07,
     .-3.0743D-08, 1.7119D-08, 8.4698D-06, 4.5395D-06, 1.5774D-06,
     . 4.1647D-07,-7.9865D-08,-1.6003D-07,-6.0171D-08, 3.1692D-06/
      DATA (A(I),I=286,380)/
     . 3.6965D-07,-4.9333D-07,-2.8856D-08, 4.2819D-08, 2.0492D-04,
     . 1.3705D-04, 4.0972D-05, 2.9565D-06,-1.4061D-06,-1.8775D-06,
     .-1.6306D-06,-3.6380D-07, 3.1213D-07, 2.0912D-07, 1.1965D-05,
     . 9.9088D-07,-1.3565D-06,-2.3128D-07, 1.1379D-05, 2.6632D-06,
     .-1.6877D-06,-4.0541D-07, 9.8921D-08, 5.9262D-03, 3.2332D-03,
     . 9.2954D-04, 1.6597D-04,-3.7972D-05,-7.1646D-05,-4.0073D-05,
     . 2.5789D-08, 5.2124D-09,-2.2517D-10,-7.9400D-10, 3.3181D-06,
     . 1.9117D-07, 1.0299D-07,-7.7443D-08, 2.2953D-06,-1.0879D-06,
     . 3.4368D-07,-1.1230D-07, 1.6826D-08, 6.9774D-06, 9.2721D-07,
     .-1.8983D-08,-1.6068D-07, 1.5843D-06,-3.5616D-08,-1.1262D-07,
     .-4.3051D-08, 8.8140D-09, 3.4244D-05, 8.0163D-06, 2.9703D-06,
     .-1.6480D-07,-6.3282D-07, 2.3791D-06,-4.9305D-07,-1.7122D-07,
     . 1.2147D-07, 7.2131D-05, 1.7699D-05, 8.4281D-06,-9.3966D-07,
     .-6.8048D-07, 5.3785D-05, 6.9745D-06,-2.3554D-06,-7.2141D-07,
     . 3.8501D-07, 4.4794D-06,-6.2435D-07,-2.7340D-07, 2.3127D-08,
     . 4.2314D-08, 3.3784D-06,-5.1356D-07,-2.4321D-07,-6.2698D-09,
     . 3.0812D-08, 5.6419D-08, 1.6889D-08, 2.7037D-09,-1.7433D-09,
     .-9.4507D-10, 1.9714D-06, 4.2743D-08,-7.4114D-08,-2.9304D-08,
     . 2.8973D-06, 1.0079D-06, 2.9561D-07,-6.2977D-09,-4.5782D-08/
      DATA (A(I),I=381,475)/
     .-2.4331D-08, 4.3284D-06, 2.8398D-07,-3.2705D-07,-1.5079D-07,
     . 5.2686D-06, 1.6539D-06, 3.6079D-07,-1.1589D-07,-5.4904D-08,
     . 7.0939D-06,-1.9534D-06, 3.2391D-08, 5.5702D-08, 3.9072D-05,
     . 1.7299D-05, 5.1530D-06,-6.0911D-07,-1.2652D-06,-4.6507D-07,
     . 1.1506D-05,-2.0422D-06,-6.1195D-07, 5.8641D-08, 1.0150D-05,
     .-4.6977D-07,-6.9446D-07, 2.2516D-08, 1.1109D-07, 3.6715D-05,
     . 5.6186D-06,-3.2309D-06,-1.6403D-06, 4.1051D-05, 2.3335D-05,
     . 7.8106D-06, 2.0279D-07,-8.1139D-07,-3.7619D-07,-1.4982D-07,
     . 2.5621D-05,-1.0798D-05, 1.4607D-06, 3.2421D-07, 6.5478D-09,
     . 2.7233D-09, 3.8646D-10,-1.9143D-10,-1.0483D-10,-4.8664D-11,
     . 1.0698D-06, 2.7752D-07,-2.6636D-08,-3.7583D-08, 2.6989D-07,
     . 3.0325D-08,-2.4613D-08,-1.0828D-08, 2.2522D-09, 8.0777D-06,
     . 2.2558D-06, 7.4760D-08,-2.4140D-07,-5.4292D-08, 8.8362D-07,
     . 9.5045D-08,-5.0543D-08,-1.9134D-08, 1.1284D-05, 4.0659D-06,
     . 7.6246D-07,-9.0394D-08,-8.7408D-08, 1.2192D-06,-2.0362D-07,
     .-1.0700D-07, 2.1990D-08, 1.2519D-08, 5.2758D-05, 2.0175D-05,
     . 3.6690D-06,-4.9014D-07,-7.0474D-07,-3.9931D-07, 4.7638D-05,
     . 1.5546D-05, 6.2294D-07,-1.3428D-06,-1.8177D-07, 4.1203D-06,
     .-4.9340D-07,-3.0198D-07,-1.5902D-08, 2.8567D-08, 2.2397D-06/
      DATA (A(I),I=476,570)/
     .-2.1306D-08,-1.7897D-07, 2.8178D-08, 4.1722D-08, 2.0594D-04,
     . 1.2838D-04, 5.8219D-05, 1.1735D-05,-6.9778D-06,-6.2486D-06,
     .-1.3800D-06, 2.9170D-06,-7.0833D-07,-1.2673D-07, 4.5069D-08,
     . 1.0974D-09, 4.7908D-10, 3.0294D-11,-7.0815D-11,-1.2039D-11,
     . 5.6357D-12, 8.1274D-07, 4.5158D-07, 1.1655D-07,-2.1481D-08,
     .-2.2493D-08,-6.5884D-09, 9.8696D-08, 3.6499D-08, 1.4768D-09,
     .-8.0141D-09,-2.8147D-09, 5.8287D-06, 3.2469D-06, 9.4927D-07,
     .-7.2550D-08,-1.4231D-07,-5.8408D-08,-1.4468D-08, 4.6070D-07,
     . 1.5956D-07,-3.2311D-09,-3.3487D-08,-7.9013D-09, 4.0092D-06,
     . 1.2195D-06, 1.1192D-07,-1.3870D-07,-6.2194D-08, 5.4918D-07,
     .-5.1133D-08,-5.4844D-08,-3.1054D-09, 6.1049D-09, 2.4833D-05,
     . 1.1280D-05, 2.2680D-06,-4.4637D-07,-2.8942D-07,-1.2382D-07,
     . 6.8223D-05, 3.6505D-05, 7.8792D-06,-1.7946D-06,-1.4925D-06,
     .-4.6018D-07, 2.5677D-06,-7.2440D-08,-2.2441D-07,-4.1769D-08,
     . 2.3833D-08, 1.1885D-06, 4.0457D-09,-1.2273D-07,-2.3271D-08,
     . 1.5149D-08, 9.1912D-05, 5.6621D-05, 1.4338D-05,-2.6072D-06,
     .-2.3688D-06,-6.2490D-07,-1.4386D-07, 1.0821D-06,-1.6304D-07,
     .-7.9756D-08,-4.9881D-09, 9.9527D-09, 1.5288D-03, 9.8517D-04,
     . 3.4966D-04, 3.4238D-05,-4.3610D-05,-3.2734D-05,-8.3676D-06/
      DATA (A(I),I=571,665)/
     . 9.6630D-09, 3.2976D-09, 2.1545D-10,-3.3938D-10,-1.4397D-10,
     . 3.4140D-07, 7.6536D-08,-2.2485D-08,-1.8244D-08,-2.0611D-09,
     . 1.3128D-06, 6.4165D-07, 1.4728D-07,-2.7600D-08,-3.0125D-08,
     .-1.0320D-08, 1.6295D-06, 3.3697D-07,-5.7547D-08,-7.3134D-08,
     .-1.5301D-08, 8.0375D-06, 4.0695D-06, 1.1441D-06,-5.8423D-08,
     .-1.6488D-07,-7.5920D-08, 2.1173D-06,-1.9184D-07,-2.4986D-07,
     . 7.2656D-09, 2.7560D-08, 7.9040D-06, 3.3963D-06, 4.8082D-07,
     .-3.1619D-07,-1.6501D-07, 5.8017D-06,-5.5173D-07,-6.3613D-07,
     . 1.9633D-08, 7.8681D-08, 9.4568D-06, 9.8227D-08,-7.8983D-07,
     .-2.0343D-07, 7.0351D-05, 3.6017D-05, 8.0504D-06,-1.7276D-06,
     .-1.5329D-06,-4.7799D-07, 3.5263D-05, 1.8639D-05, 4.2070D-06,
     .-4.2643D-07,-5.2726D-07,-2.5481D-07,-1.0144D-07, 4.4613D-06,
     .-8.7802D-07,-3.9633D-07, 6.7953D-08, 4.1539D-08, 1.8780D-04,
     . 1.0547D-04, 2.0983D-05,-3.7076D-06,-3.7090D-06,-1.5934D-06,
     .-2.6146D-07, 1.6384D-05,-3.3765D-06,-7.8390D-07, 1.2382D-07,
     . 1.5748D-05,-9.4352D-07,-4.8936D-07,-4.9967D-07, 3.8177D-10,
     . 9.6781D-11,-1.9622D-11,-1.7859D-11,-4.4781D-12, 3.0310D-07,
     . 1.1193D-07,-5.4059D-09,-1.4249D-08,-1.9549D-09, 4.7664D-08,
     . 8.5684D-09,-5.0948D-09,-3.1313D-09, 4.9557D-10, 4.2702D-10/
      DATA (A(I),I=666,760)/
     . 2.3433D-06, 8.3824D-07, 2.6585D-08,-7.5944D-08,-1.5093D-08,
     . 1.8775D-07, 2.7614D-08,-2.1193D-08,-1.0264D-08, 1.9660D-09,
     . 1.1637D-09, 7.7890D-06, 3.2901D-06,-2.9753D-07,-5.3153D-07,
     .-5.2098D-08, 3.3947D-08, 5.5937D-07, 5.1363D-08,-8.5390D-08,
     .-2.2580D-08, 1.0857D-08, 3.5157D-09, 4.1303D-05, 2.1947D-05,
     . 3.9480D-06,-1.4234D-06,-9.0220D-07,-2.4485D-07, 1.1031D-04,
     . 6.6726D-05, 1.9581D-05,-4.3637D-07,-2.8911D-06,-1.4332D-06,
     .-3.2332D-07, 3.0922D-06, 2.0624D-07,-4.0771D-07,-1.0402D-07,
     . 4.8807D-08, 1.2491D-06, 2.1187D-07,-1.5951D-07,-7.2537D-08,
     . 1.4035D-08, 9.6514D-09, 2.5146D-05,-9.4490D-08,-3.4103D-06,
     . 2.4020D-07, 2.9564D-07, 6.6566D-07,-3.0772D-08,-1.0055D-07,
     .-3.1082D-09, 1.4897D-08, 7.7189D-04, 3.3710D-04, 3.9258D-05,
     .-1.5271D-05,-6.0652D-06, 2.1986D-04,-3.3135D-05,-8.5682D-06,
     .-9.2988D-07, 3.8085D-06,-1.3259D-07,-3.8517D-07,-5.7276D-08,
     . 3.7226D-08, 2.0825D-09, 1.6707D-10,-1.2443D-10,-3.9117D-11,
     . 1.0068D-07, 6.2278D-09,-5.9169D-09,-3.9535D-09, 5.3334D-07,
     . 2.2022D-07, 1.4522D-08,-2.8440D-08,-8.1406D-09, 6.0187D-07,
     . 5.8318D-08,-2.7772D-08,-2.6639D-08, 3.0929D-06, 1.0648D-06,
     . 1.2317D-07,-9.9233D-08,-4.0796D-08, 1.5217D-06, 8.3095D-08/
      DATA (A(I),I=761,855)/
     .-1.9819D-07,-6.0417D-08, 2.9553D-08, 9.0905D-09, 8.6651D-06,
     . 3.9125D-06,-2.2871D-07,-6.4651D-07,-7.4440D-08, 4.3543D-08,
     . 4.4505D-06, 2.7350D-07,-4.6428D-07,-2.2293D-07, 5.5275D-08,
     . 3.6505D-08, 8.3471D-06, 5.0215D-07,-1.1189D-06,-2.5404D-07,
     . 1.3175D-07, 1.1687D-04, 7.0674D-05, 2.0219D-05,-1.2136D-06,
     .-3.2299D-06,-1.3981D-06,-2.8271D-07, 3.7298D-05, 2.2054D-05,
     . 5.5510D-06,-9.1515D-07,-1.0832D-06,-3.6657D-07,-4.8581D-08,
     . 2.2506D-06,-2.9469D-07,-2.1641D-07,-3.9630D-08, 5.2879D-08,
     . 2.7894D-05,-3.9312D-06,-2.6413D-06, 5.5848D-07, 8.7829D-06,
     .-1.4800D-06,-7.0982D-07,-2.1645D-08, 1.3258D-07, 1.2248D-05,
     .-1.4244D-06,-7.6696D-07,-1.5029D-07, 6.8467D-08, 1.9392D-04,
     .-2.6676D-05,-8.0536D-06,-6.3609D-07, 2.3074D-05, 7.2856D-07,
     .-2.3751D-06,-5.8464D-07, 1.3685D-07, 1.9646D-08, 1.2918D-08,
     . 4.6924D-09, 5.1593D-10,-4.4715D-10,-3.4091D-10,-9.7254D-11,
     . 3.4117D-07, 1.2731D-07,-2.2300D-08,-2.5297D-08,-1.1877D-10,
     . 2.3436D-09, 9.4463D-07, 4.9464D-07, 8.9138D-08,-2.2158D-08,
     .-1.2008D-08,-4.3715D-09,-2.8719D-09, 1.4091D-06, 4.8184D-07,
     .-9.6135D-08,-1.0945D-07,-9.3174D-09, 9.0166D-09, 4.3747D-06,
     . 2.1053D-06, 1.7702D-07,-3.0956D-07,-1.2902D-07,-8.9418D-09/
      DATA (A(I),I=856,929)/
     . 1.4297D-06, 7.7612D-08,-1.7188D-07, 6.4348D-10, 1.1572D-08,
     . 7.4828D-06, 4.0177D-06, 1.0890D-06, 8.6393D-08,-5.4701D-08,
     .-5.7656D-08,-3.3030D-08, 5.0379D-06, 1.2861D-07,-7.1313D-07,
     .-3.3703D-09, 7.9367D-08, 6.2235D-06, 9.8196D-07,-4.6416D-07,
     .-2.3278D-07, 2.9519D-05, 1.5275D-05, 3.1448D-06,-6.4571D-07,
     .-3.4445D-07, 3.6508D-05, 2.0524D-05, 3.1203D-06,-2.3606D-06,
     .-1.5012D-06,-2.5807D-07, 1.3807D-07, 9.9656D-08, 3.1628D-06,
     .-5.6361D-08,-3.8861D-07,-1.0385D-09, 2.9930D-08, 3.1351D-04,
     . 2.3774D-04, 1.0342D-04, 1.2429D-05,-1.4107D-05,-9.2597D-06,
     .-1.7338D-06, 8.4903D-07, 7.3795D-07, 2.2065D-07, 1.2134D-05,
     . 4.6651D-07,-9.3659D-07,-2.7485D-07, 1.2462D-05, 6.2559D-07,
     .-6.1634D-07,-5.1080D-07, 1.2493D-02, 8.2057D-03, 2.9562D-03,
     . 3.3090D-04,-3.7349D-04,-2.8886D-04,-7.4606D-05, 1.2726D-05,
     .-1.6480D-07,-1.5006D-06,-1.1181D-07, 1.4846D-07, 8.7160D-03,
     . 4.2652D-03, 5.2455D-04,-2.6363D-04,-7.9836D-05/
      SAVE ENER,NSTART,A,STWT
C
C          SET ALL ELEMENTS OF COLHE1 = 0
C
           DO 10 I=1,N
              DO 10 J=1,N
   10            COLHE1(I,J)=0.
C
C        EVALUATE GAMMAS FOR REQUESTED LEVELS
C
      XXX=2.D0*(LOG10(TEMP)-C1)/C2
      TFAC=UN/SQRT(TEMP)
C
C          LOOPS OVER LEVELS
C
      DO 30 IL=1,N-1
           DO 30 IU=IL+1,N
                J=((IU*IU-3*IU+4)/2)+IL-1
                N1=NSTART(J)
                NF=NSTART(J+1)-1
                NT=NF-N1+1
                NTM2=NT-2
C
C          CLENSHAW SUMMATION
C
                B(NT)=A(NF)
                B(NT-1)=XXX*B(NT)+A(NF-1)
                IR=NTM2
                JJ=NF-2
                DO 20 J=1,NTM2
                     B(IR)=XXX*B(IR+1)-B(IR+2)+A(JJ)
                     IR=IR-1
                     JJ=JJ-1
   20           CONTINUE
                COLHE1(IU,IL)=(B(1)-B(3))*TFAC
                X=(ENER(IL)-ENER(IU))/8.62D-5/TEMP
                COLHE1(IL,IU)=COLHE1(IU,IL)*STWT(IU)/STWT(IL)*EXP(X)
   30  CONTINUE
       RETURN
       END
C
C
C     ****************************************************************
C
C
      FUNCTION CHEAV(II,JJ,IC)
C     ========================
C
C     Calculates collisional excitation rates of neutral helium
C     between states with n= 1, 2, 3, 4; with either the upper state
C     alone, or both upper and lower states are some averaged states
C     The program allows only two standard possibilities of
C     constructing averaged levels:
C     i)  all states within given principal quantum number n (>1) are
C         lumped together
C     ii) all siglet states for given n, and all triplet states for
C         given n are lumped together separately (there are thus two
C         explicit levels for a given n)
C
C     The rates are calculated using appropriate summations and/or
C     averages of the Storey-Hummer rates (calculated by procedure
C     COLLHE and stored in array COLHE1)
C
C     Input parameters:
C      II,JJ - indices of the lower and the upper level (in the
C              numbering of the explicit levels)
C      IC    - collisional switch ICOL for the given transition
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'ATOMIC.FOR'
C
      CHEAV=0.
      NI=NQUANT(II)
      NJ=NQUANT(JJ)
      IGI=INT(G(II)+0.01)
      IGJ=INT(G(JJ)+0.01)
C
C     ----------------------------------------------------------------
C     IC=2 - transition from an (l,s) lower level to an averaged upper
C            level
C     ----------------------------------------------------------------
C
      IF(IC.EQ.2) THEN
         I=II-NFIRST(IELHE1)+1
         CHEAV=CHEAVJ(I,NJ,IGJ)
      END IF
C
C     ----------------------------------------------------------------
C     IC=3 - transition from an averaged lower level to an averaged
C            upper level
C     ----------------------------------------------------------------
C
      IF(IC.EQ.3) THEN
         IF(NI.EQ.2) THEN
C
C ********    transitions from an averaged level with n=2
C
            IF(IGI.EQ.4) THEN
C
C      a) lower level is an averaged singlet state
C
               CHEAV=(CHEAVJ(3,NJ,IGJ)+3.D0*CHEAVJ(5,NJ,IGJ))/4.D0
            ELSE IF(IGI.EQ.12) THEN
C
C      b) lower level is an averaged triplet state
C
               CHEAV=(CHEAVJ(2,NJ,IGJ)+3.D0*CHEAVJ(4,NJ,IGJ))/4.D0
            ELSE IF(IGI.EQ.16) THEN
C
C      c) lower level is an average of both singlet and triplet states
C
               CHEAV=(CHEAVJ(3,NJ,IGJ)+3.D0*(CHEAVJ(5,NJ,IGJ)+
     *                CHEAVJ(2,NJ,IGJ))+9.D0*CHEAVJ(4,NJ,IGJ))/1.6D1
            ELSE
               GO TO 10
            END IF
C
C
C ********    transitions from an averaged level with n=3
C
         ELSE IF(NI.EQ.3) THEN
            IF(IGI.EQ.9) THEN
C
C      a) lower level is an averaged singlet state
C
               CHEAV=(CHEAVJ(7,NJ,IGJ)+3.D0*CHEAVJ(11,NJ,IGJ)+
     *                5.D0*CHEAVJ(10,NJ,IGJ))/9.D0
            ELSE IF(IGI.EQ.27) THEN
C
C      b) lower level is an averaged triplet state
C
               CHEAV=(CHEAVJ(6,NJ,IGJ)+3.D0*CHEAVJ(8,NJ,IGJ)+
     *                5.D0*CHEAVJ(9,NJ,IGJ))/9.D0
            ELSE IF(IGI.EQ.36) THEN
C
C      c) lower level is an average of both singlet and triplet states
C
               CHEAV=(CHEAVJ(7,NJ,IGJ)+3.D0*CHEAVJ(11,NJ,IGJ)+
     *                5.D0*CHEAVJ(10,NJ,IGJ)+
     *                3.D0*CHEAVJ(6,NJ,IGJ)+9.D0*CHEAVJ(8,NJ,IGJ)+
     *                1.5D1*CHEAVJ(9,NJ,IGJ))/3.6D1
            ELSE
               GO TO 10
            END IF
C
C ********    transitions from an averaged level with n=4
C
         ELSE IF(NI.EQ.4) THEN
            IF(IGI.EQ.16) THEN
C
C      a) lower level is an averaged singlet state
C
               CHEAV=(CHEAVJ(13,NJ,IGJ)+
     *                3.D0*CHEAVJ(19,NJ,IGJ)+
     *                5.D0*CHEAVJ(16,NJ,IGJ)+
     *                7.D0*CHEAVJ(18,NJ,IGJ))/1.6D1
            ELSE IF(IGI.EQ.48) THEN
C
C      b) lower level is an averaged triplet state
C
               CHEAV=(CHEAVJ(12,NJ,IGJ)+
     *                3.D0*CHEAVJ(14,NJ,IGJ)+
     *                5.D0*CHEAVJ(15,NJ,IGJ)+
     *                7.D0*CHEAVJ(17,NJ,IGJ))/1.6D1
            ELSE IF(IGI.EQ.64) THEN
C
C      c) lower level is an average of both singlet and triplet states
C
               CHEAV=(CHEAVJ(13,NJ,IGJ)+
     *                3.D0*CHEAVJ(19,NJ,IGJ)+
     *                5.D0*CHEAVJ(16,NJ,IGJ)+
     *                7.D0*CHEAVJ(18,NJ,IGJ)+
     *                3.D0*CHEAVJ(12,NJ,IGJ)+
     *                9.D0*CHEAVJ(14,NJ,IGJ)+
     *                15.D0*CHEAVJ(15,NJ,IGJ)+
     *                21.D0*CHEAVJ(17,NJ,IGJ))/6.4D1
            ELSE
               GO TO 10
            END IF
         ELSE
            GO TO 10
         END IF
      END IF
      RETURN
      
   10 WRITE(6,601) NI,NJ,IGI,IGJ
      write(10,601) NI,NJ,IGI,IGJ
  601 FORMAT(1H0/' INCONSISTENT INPUT TO PROCEDURE CHEAV'/
     * ' QUANTUM NUMBERS =',2I3,'    STATISTICAL WEIGHTS',2I4)
      call quit(' ',ni,nj)
      
      END
C
C
C     ****************************************************************
C
C

      FUNCTION CHEAVJ(I,NJ,IGJ)
C     =========================
C
C     Calculates collisional excitation rates from a non-averaged (l,s)
C     state of He I, with n=1, 2, 3,  to some averaged state
C     with n = 2, 3, 4.
C
C     The rates are calculated using appropriate summations of the
C     Storey-Hummer rates (calculated by procedure COLLHE, and stored
C     in array COLHE1)
C
C     Input:
C     I   -  index of the lower state, using the ordering defined in
C            COLLHE, ie. I=1 for 1 sing S, I=2 for 2 trip S, etc.
C     NJ  - principal quantum number of the (averaged) upper level
C     IGJ - statistical weight of the upper level
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'ATOMIC.FOR'
      
      CHEAVJ=0.
C
C -----------------------------------------------------
C ********    transitions to an averaged level with n=2
C -----------------------------------------------------
C
      IF(NJ.EQ.2) THEN
         IF(IGJ.EQ.4) THEN
C
C      a) upper level is an averaged singlet state
C
            CHEAVJ=COLHE1(1,3)+COLHE1(1,5)
         ELSE IF(IGJ.EQ.12) THEN
C
C      b) upper level is an averaged triplet state
C
            CHEAVJ=COLHE1(1,2)+COLHE1(1,4)
         ELSE IF(IGJ.EQ.16) THEN
C
C      c) upper level is an average of both siglet and triplet states
C
            CHEAVJ=COLHE1(1,3)+COLHE1(1,5)+COLHE1(1,2)+COLHE1(1,4)
         ELSE
            GO TO 10
         END IF
C
C -----------------------------------------------------
C ********    transitions to an averaged level with n=3
C -----------------------------------------------------
C
       ELSE IF(NJ.EQ.3) THEN
         IF(IGJ.EQ.9) THEN
C
C      a) upper level is an averaged singlet state
C
            CHEAVJ=COLHE1(I,7)+COLHE1(I,11)+COLHE1(I,10)
         ELSE IF(IGJ.EQ.27) THEN
C
C      b) upper level is an averaged triplet state
C
            CHEAVJ=COLHE1(I,6)+COLHE1(I,8)+COLHE1(I,9)
         ELSE IF(IGJ.EQ.36) THEN
C
C      c) upper level is an average of both siglet and triplet states
C
            CHEAVJ=COLHE1(I,7)+COLHE1(I,11)+COLHE1(I,10)+
     *             COLHE1(I,6)+COLHE1(I,8)+COLHE1(I,9)
         ELSE
            GO TO 10
         END IF
C
C -----------------------------------------------------
C ********    transitions to an averaged level with n=4
C -----------------------------------------------------
C
       ELSE IF(NJ.EQ.4) THEN
         IF(IGJ.EQ.16) THEN
C
C      a) upper level is an averaged singlet state
C
            CHEAVJ=COLHE1(I,13)+COLHE1(I,19)+COLHE1(I,16)+
     *             COLHE1(I,18)
         ELSE IF(IGJ.EQ.48) THEN
C
C      b) upper level is an averaged triplet state
C
            CHEAVJ=COLHE1(I,12)+COLHE1(I,14)+COLHE1(I,15)+
     *             COLHE1(I,17)
         ELSE IF(IGJ.EQ.64) THEN
C
C      c) upper level is an average of both siglet and triplet states
C
            CHEAVJ=COLHE1(I,13)+COLHE1(I,19)+COLHE1(I,16)+
     *             COLHE1(I,18)+COLHE1(I,12)+COLHE1(I,14)+
     *             COLHE1(I,15)+COLHE1(I,17)
         ELSE
            GO TO 10
         END IF
       ELSE
         GO TO 10
      END IF
      RETURN
      
   10 WRITE(6,601) NJ,IGJ
      WRITE(10,601) NJ,IGJ
  601 FORMAT(1H0/' INCONSISTENT INPUT TO PROCEDURE CHEAVJ'/
     * ' QUANTUM NUMBER =',I3,'    STATISTICAL WEIGHT',2I4)
      call quit(' ',nj,igj)
      
      END
C
C
C     ****************************************************************
C
C

      SUBROUTINE MATINV(A,N,NR)
C     =========================
C
C     Matrix inversion
C      by LU decomposition
C
C      A  -  matrix of actual size (N x N) and maximum size (NR x NR)
C            to be inverted;
C      Inversion is accomplished in place and the original matrix is
C      replaced by its inverse
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      DIMENSION A(NR,NR)
C
C  --------- Cray version ------------------
C
c      parameter(mt2=2*mtot)
c      dimension scr(mt2)
c      call minv(a,n,nr,scr,det,1e-135,0,1)
c      return
c      end
c  _________________________________________
C
C  For Cray - all what is needed is to comment out the rest of subroutine
C  Note:
C  Obviously, the normal version works for Cray as well, but MINV is a
C  routine devised specifically for Cray and is therefore optimalized.
C  Numerical experience shows that the execution time when using MINV
C  is about 75 - 80 % of the time using ordinary MATINV.
C
C     Normal version follows:
C
      DO 50 I=2,N
         IM1=I-1
         DO 20 J=1,IM1
            JM1=J-1
            DIV=A(J,J)
            SUM=0.
            IF(JM1.LT.1) GO TO 20
            DO 10 K=1,JM1
   10          SUM=SUM+A(I,K)*A(K,J)
   20       A(I,J)=(A(I,J)-SUM)/DIV
         DO 40 J=I,N
            SUM=0.
            DO 30 K=1,IM1
   30          SUM=SUM+A(I,K)*A(K,J)
   40       A(I,J)=A(I,J)-SUM
   50 CONTINUE
      DO 80 II=2,N
         I=N+2-II
         IM1=I-1
         IF(IM1.LT.1) GO TO 80
         DO 70 JJ=1,IM1
            J=I-JJ
            JP1=J+1
            SUM=0.
            IF(JP1.GT.IM1) GO TO 70
        DO 60 K=JP1,IM1
   60          SUM=SUM+A(I,K)*A(K,J)
   70    A(I,J)=-A(I,J)-SUM
   80 CONTINUE
      DO 110 II=1,N
         I=N+1-II
         DIV=A(I,I)
         IP1=I+1
         IF(IP1.GT.N) GO TO 110
         DO 100 JJ=IP1,N
            J=N+IP1-JJ
            SUM=0.
            DO 90 K=IP1,J
   90          SUM=SUM+A(I,K)*A(K,J)
            A(I,J)=-SUM/DIV
  100    CONTINUE
  110 A(I,I)=1.0D0/A(I,I)
C 
      DO 240 I=1,N
         DO 230 J=1,N
            K0=I
            IF(J.GE.I) GO TO 220
            SUM=0.
  200       DO 210 K=K0,N
  210          SUM=SUM+A(I,K)*A(K,J)
            GO TO 230
  220       K0=J
            SUM=A(I,K0)
            IF(K0.EQ.N) GO TO 230
            K0=K0+1
            GO TO 200
  230    A(I,J)=SUM
  240 CONTINUE
      RETURN
      END
C
C
C     ****************************************************************
C
C

      SUBROUTINE MINV3(A)
C     ===================
C
C     Special routine for an invresion of a 3 x 3 matrix
C
      INCLUDE 'IMPLIC.FOR'
      PARAMETER (UN=1.D0)
      DIMENSION A(3,3)
C
      A(2,1)=A(2,1)/A(1,1)
      A(2,2)=A(2,2)-A(2,1)*A(1,2)
      A(2,3)=A(2,3)-A(2,1)*A(1,3)
      A(3,1)=A(3,1)/A(1,1)
      A(3,2)=(A(3,2)-A(3,1)*A(1,2))/A(2,2)
      A(3,3)=A(3,3)-A(3,1)*A(1,3)-A(3,2)*A(2,3)
C
      A(3,2)=-A(3,2)
      A(3,1)=-A(3,1)-A(3,2)*A(2,1)
      A(2,1)=-A(2,1)
C
      A(3,3)=UN/A(3,3)
      A(2,3)=-A(2,3)*A(3,3)/A(2,2)
      A(2,2)=UN/A(2,2)
      A(1,3)=-(A(1,2)*A(2,3)+A(1,3)*A(3,3))/A(1,1)
      A(1,2)=-A(1,2)*A(2,2)/A(1,1)
      A(1,1)=UN/A(1,1)
C
      A(1,1)=A(1,1)+A(1,2)*A(2,1)+A(1,3)*A(3,1)
      A(1,2)=A(1,2)+A(1,3)*A(3,2)
      A(2,1)=A(2,2)*A(2,1)+A(2,3)*A(3,1)
      A(2,2)=A(2,2)+A(2,3)*A(3,2)
      A(3,1)=A(3,3)*A(3,1)
      A(3,2)=A(3,3)*A(3,2)
C
      RETURN
      END
C
C
C     ******************************************************************
C
C

      SUBROUTINE LINEQS(A,B,X,N,NR)
C     =============================
C
C     Solution of the linear system A*X=B
C     by Gaussian elimination with partial pivoting
C
C     Input: A  - matrix of the linear system, with actual size (N x N),
C                and maximum size (NR x NR)
C            B  - the rhs vector
C     Output: X - solution vector
C     Note that matrix A and vector B are destroyed here
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      DIMENSION A(NR,NR),B(NR),X(NR),D(MLEVEL),IP(MLEVEL)
c
      if(n.eq.2) then
        a11=a(1,1)
        a12=a(1,2)
        a21=a(2,1)
        a22=a(2,2)
        x(1)=(a(2,2)*b(1)-a(1,2)*b(2))/
     *       (a(1,1)*a(2,2)-a(1,2)*a(2,1))
        x(2)=(b(2)-a(2,1)*x(1))/a(2,2)
        return
      end if
c
      DO 70 I=1,N
         DO 10 J=1,N
   10       D(J)=A(J,I)
         IM1=I-1
         IF(IM1.LT.1) GO TO 40
         DO 30 J=1,IM1
            IT=IP(J)
            A(J,I)=D(IT)
            D(IT)=D(J)
            JP1=J+1
            DO 20 K=JP1,N
   20          D(K)=D(K)-A(K,J)*A(J,I)
   30    CONTINUE
   40    AM=ABS(D(I))
         IP(I)=I
         DO 50 K=I,N
            IF(AM.GE.ABS(D(K))) GO TO 50
            IP(I)=K
            AM=ABS(D(K))
   50    CONTINUE
         IT=IP(I)
         A(I,I)=D(IT)
         D(IT)=D(I)
         IP1=I+1
         IF(IP1.GT.N) GO TO 80
         DO 60 K=IP1,N
   60       A(K,I)=D(K)/A(I,I)
   70 CONTINUE
   80 DO 100 I=1,N
         IT=IP(I)
         X(I)=B(IT)
         B(IT)=B(I)
         IP1=I+1
         IF(IP1.GT.N) GO TO 110
         DO 90 J=IP1,N
   90       B(J)=B(J)-A(J,I)*X(I)
  100 CONTINUE
  110 DO 140 I=1,N
         K=N-I+1
         SUM=0.
         KP1=K+1
         IF(KP1.GT.N) GO TO 130
         DO 120 J=KP1,N
  120       SUM=SUM+A(K,J)*X(J)
  130    X(K)=(X(K)-SUM)/A(K,K)
  140 CONTINUE
      RETURN
      END
C
C
C     ****************************************************************
C
C

      FUNCTION EXPINT(X)
C     ==================
C
C     First exponential integral function E1(X)
C
      INCLUDE 'IMPLIC.FOR'
      PARAMETER (A1  =  -0.57721566,
     *           A2  =   0.99999193,
     *           A3  =  -0.24991055,
     *           A4  =   0.05519968,
     *           A5  =  -0.00976004,
     *           A6  =   0.00107857,
     *           B1  =   0.2677734343,
     *           B2  =   8.6347608925,
     *           B3  =  18.059016973,
     *           B4  =   8.5733287401,
     *           C1  =   3.9584969228,
     *           C2  =  21.0996530827,
     *           C3  =  25.6329561486,
     *           C4  =   9.5733223454,
     *           UN  =   1.0)
C
      IF(X.LE.UN) THEN
         EXPINT=-LOG(X)+A1+X*(A2+X*(A3+X*(A4+X*(A5+X*A6))))
       ELSE
         EXPINT=EXP(-X)*((B1+X*(B2+X*(B3+X*(B4+X))))/
     *                   (C1+X*(C2+X*(C3+X*(C4+X)))))/X
      END IF
      RETURN
      END
C
C
C     ****************************************************************
C
C

      SUBROUTINE INTERP(X,Y,XX,YY,NX,NXX,NPOL,ILOGX,ILOGY)
C     ====================================================
C
C     General interpolation procedure of the (NPOL-1)-th order
C
C     for  ILOGX = 1  logarithmic interpolation in X
C     for  ILOGY = 1  logarithmic interpolation in Y
C
C     Input:
C      X    - array of original x-coordinates
C      Y    - array of corresponding functional values Y=y(X)
C      NX   - number of elements in arrays X or Y
C      XX   - array of new x-coordinates (to which is to be
C             interpolated
C      NXX  - number of elements in array XX
C     Output:
C      YY   - interpolated functional values YY=y(XX)
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      DIMENSION X(MDEPTH),Y(MDEPTH),XX(MDEPTH),YY(MDEPTH)
C
C     no interpolation for NPOL.LE.0 or NX.le.0
C
      IF(NPOL.LE.0.OR.NX.LE.0) THEN
         N=NX
         IF(NXX.GE.NX) N=NXX
         DO I=1,N
            XX(I)=X(I)
            YY(I)=Y(I)
         END DO
       RETURN
      END IF
C
C     interpolation
C
C     if required, compute logarithms to be interpolated
C
      IF(ILOGX.GT.0) THEN
      DO I=1,NX
         X(I)=LOG10(X(I))
      END DO
      DO I=1,NXX
         XX(I)=LOG10(XX(I))
      END DO
      END IF 
      IF(ILOGY.GT.0) THEN
      DO I=1,NX
         Y(I)=LOG10(Y(I))
      END DO
      END IF 
C
      NM=(NPOL+1)/2
      NM1=NM+1
      NUP=NX+NM1-NPOL
      DO 100 ID=1,NXX
         XXX=XX(ID)
         DO 60 I=NM1,NUP
            IF(XXX.LE.X(I)) GO TO 70
   60    CONTINUE
         I=NUP
   70    J=I-NM
         JJ=J+NPOL-1
         YYY=0.
         DO 90 K=J,JJ
            T=1.
            DO 80 M=J,JJ
               IF(K.EQ.M) GO TO 80
               T=T*(XXX-X(M))/(X(K)-X(M))
   80       CONTINUE
   90    YYY=Y(K)*T+YYY
         YY(ID)=YYY
  100 CONTINUE
C
      IF(ILOGX.GT.0) THEN
      DO I=1,NX
         X(I)=EXP(X(I)*2.30258509299405)
      END DO
      DO I=1,NXX
         XX(I)=EXP(XX(I)*2.30258509299405)
      END DO
      END IF
      IF(ILOGY.GT.0) THEN
      DO I=1,NX
         Y(I)=EXP(Y(I)*2.30258509299405)
      END DO
      DO I=1,NXX
         YY(I)=EXP(YY(I)*2.30258509299405)
      END DO
      END IF
C
      RETURN
      END
C
C
C     ****************************************************************
C
C
      SUBROUTINE OUTPUT
C     =================
C
C     Output of computed model atmosphere on file 7
C     This file may be used as input file 8 (initial model atmosphere)
C     for a subsequent run of the program
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'MODELQ.FOR'

C
      NUMLT=3
      IF(IDISK.EQ.1) NUMLT=4
      NUMPAR=NLEVEL+NUMLT
      IF(LTE.AND.IPRINP.EQ.0) NUMPAR=NUMLT
C
C     NUMPAR  - number of model parameters in each depth
C             = NUMLT  for LTE model, ie. TEMP - temperature
C                                     ELEC - electron density
C                                     DENS - density
C             = NUMLT+NLEVEL  for NLTE model, ie. the above + populations
C ---------------------------------------------------------------------
C      2. DM(ID),ID=1,ND - mass-depth points for the input model
C ---------------------------------------------------------------------
C      3. for each depth:
C                 T   - temperature
C                 ANE - electron density
C                 RHO - mass density
C                 level populations - only for NLTE input model
C
      REWIND 7
      WRITE(7,501) ND,NUMPAR
      WRITE(7,502) (DM(ID),ID=1,ND)
      IF(IDISK.EQ.0) THEN
      DO ID=1,ND
         IF(LTE.AND.IPRINP.EQ.0) THEN
            WRITE(7,503) TEMP(ID),ELEC(ID),DENS(ID)
          ELSE
            WRITE(7,503) TEMP(ID),ELEC(ID),DENS(ID),
     *                   (POPUL(J,ID),J=1,NLEVEL)
         END IF
      END DO
      ELSE
      DO ID=1,ND
         IF(LTE.AND.IPRINP.EQ.0) THEN
            WRITE(7,503) TEMP(ID),ELEC(ID),
     *                                   DENS(ID),ZD(ID)
          ELSE
            WRITE(7,503) TEMP(ID),ELEC(ID),DENS(ID),ZD(ID),
     *                   (POPUL(J,ID),J=1,NLEVEL)
         END IF
      END DO
      END IF
      CLOSE(7)

      IF(IPRIND.GT.0) THEN
        WRITE(17,501) ND,NUMPAR
        WRITE(17,502) (DM(I),I=1,ND)
        IF(IDISK.EQ.0) THEN
        IF(LTE) THEN
        DO ID=1,ND
          WRITE(17,503) TEMP(ID),ELEC(ID),DENS(ID)
        END DO
        ELSE
        WRITE(20,501) ND,NUMPAR
          WRITE(20,502) (DM(ID),ID=1,ND)
        DO ID=1,ND
           WRITE(17,503) TEMP(ID),ELEC(ID),DENS(ID),
     *                  (POPUL(J,ID),J=1,NLEVEL)
           WRITE(20,503) TEMP(ID),ELEC(ID),DENS(ID),
     *                   (BFAC(J,ID),J=1,NLEVEL)
        END DO
      ENDIF
        ELSE
        IF(LTE) THEN
        DO ID=1,ND
          WRITE(17,503) TEMP(ID),ELEC(ID),DENS(ID),ZD(ID)
        END DO
        ELSE
        WRITE(20,501) ND,NUMPAR
          WRITE(20,502) (DM(ID),ID=1,ND)
        DO ID=1,ND
           WRITE(17,503) TEMP(ID),ELEC(ID),DENS(ID),ZD(ID),
     *                  (POPUL(J,ID),J=1,NLEVEL)
           WRITE(20,503) TEMP(ID),ELEC(ID),DENS(ID),ZD(ID),
     *                   (BFAC(J,ID),J=1,NLEVEL)
        END DO
      ENDIF
        ENDIF
      ENDIF
  501 FORMAT(2I5)
  502 FORMAT(1P6D13.6)
  503 FORMAT(1P5D15.6)
      RETURN
      END
C
C
C     ****************************************************************
C
C

      SUBROUTINE OPADD(MODE,ICALL,IJ,ID)
C     ==================================
C
C     Additional opacities
C     This is basically user-supplied procedure; here are some more
C     important non-standard opacity sources, namely
C     Rayleigh scattering, H- opacity, H2+ opacity, and additional
C     opacity of He I and He II.
C     Inclusion of these opacities is contolled by switches transmitted
C     by COMMON/OPCKEY - see description in START.
C
C     Input parameters:
C     MODE  - controls the nature and the amount of calculations
C           = -1 - (OPADD called from START) evaluation of relevant
C                  depth-dependent quantities (usually photoionization
C                  cross-sections, but also possibly other), which are
C                  stored in array CROS
C           = 0  - evaluation of an additional opacity, emissivity, and
C                  scattering - for procedure OPACF1
C           = 1  - the same, but also derivatives wrt temperature,
C                  electron density, and level populations
C     IJ    - frequency index
C     ID    - depth index
C
C     Output: transmitted by COMMON/OPACAD:
C
C     ABAD  - absorption coefficient (at frequency point IJ and depth ID)
C     EMAD  - emission coefficient (at frequency point IJ and depth ID)
C     SCAD  - scattering coefficient (at frequency point IJ and depth ID)
C     Dxy   - derivatives of x (=A for Absorption, =E for Emission,
C             =S for scattering) coefficient wrt y (=T for temperature,
C             =N for electron density)
C     DDN(I)  - quantity proportional to a derivative of absorption (and
C               emission) coefficient wrt population of I-th level; ie.
C               d(abs)/d(pop) = DDN * [1-exp(-h*nu/kT)]
C               d(em)/d(pop) = DDN * (2h*nu**3/c**2)*exp(-h*nu/kT)
C
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'ATOMIC.FOR'
      INCLUDE 'MODELQ.FOR'
      PARAMETER (FRRAY =  2.463D15,
     *           C18   =  2.997925D18,
     *           CR0   =  5.799D-13,
     *           CR1   =  1.422D-6,
     *           CR2   =  2.784D0,
     *           TENM4 =  1.0D-4,
     *           THM0  =  8.7629D3,
     *           SBHM  =  1.0353D-16,
     *           TRHA  =  1.5D0,
     *           SFF0  =  1.3727D-25,
     *           SFF1  =  4.3748D-10,
     *           SFFM2 =  -2.5993D-7,
     *           F0HE1 =  3.29D15,
     *           F0HE2 =  1.316D16,
     *           SBH0  =  4.1412D-16,
     *           SG01  =  2.815D-16,
     *           SG02  =  4.504D-15)
      SAVE T,DELTAT,ANE,HKT,T32,XHM,POPI,SB00
C
      AB0=0.
      AB1=0.
      DAB0=0.
      DAB1=0.
      ABAD=0.
      EMAD=0.
      SCAD=0.
      DBF=0.
      DFF=0.
      DAT=0.
      DET=0.
      DAN=0.
      DEN=0.
      DST=0.
      DSN=0.
      DO 10 I=1,NLEVEL
   10    DDN(I)=0.
C
      FR=FREQ(IJ)
      if(ielh.gt.0) then
      N0HN=NFIRST(IELH)
      NKH=NKA(IATH)
      end if
C
      IF(ICALL.GT.0) THEN
         T=TEMP(ID)
         DELTAT=TENM4*T
         ANE=ELEC(ID)
         HKT=HK/T
         T32=UN/T/SQRT(T)
         XHM=THM0/T
         if(ielh.gt.0) POPI=POPUL(N0HN,ID)
         SB00=SBHM*T32*EXP(XHM)*POPI*ANE
      END IF
C
      IT=NCON
C
C   -----------------------
C   HI  Rayleigh scattering
C   -----------------------
C
      IF(IRSCT.NE.0) THEN
      IT=IT+1
      IF(MODE.LT.0) THEN
         FRM=MIN(FR,FRRAY)
         X=(C18/FRM)**2
         if(it.gt.mcross) 
     *    CALL QUIT('it.gt.mcross in opadd',it,mcross)
         if(ielh.le.0) 
     *    CALL QUIT('conflict in irstc - no hydrogen',irsct,ielh)
         BFCS(IT,IJ)=real((CR0+(CR1+CR2/X)/X)/X/X)
       ELSE
         ABAD=POPI*CROSS(IT,IJ)
         SCAD=ABAD
      END IF
      END IF
      IF(IOPHMI.NE.0) THEN
C
C   ----------------------------
C   H-  bound-free and free-free
C   ----------------------------
C     Note: IOPHMI must not by taken non-zero if H- is considered
C           explicitly, because H- opacity would be taken twice
C
      IT=IT+1
      IF(MODE.LT.0) THEN
         if(it.gt.mcross) 
     *    CALL QUIT('it.gt.mcross in opadd',it,mcross)
         if(ielh.le.0) 
     *    CALL QUIT('conflict in iophmi - no hydrogen',iophmi,ielh)
         BFCS(IT,IJ)=real(SBFHMI(FR))
       ELSE
          SB=SB00*CROSS(IT,IJ)
          SF=SFF0+(SFF1+SFFM2/T)/FR
          PFR=POPI/FR*ANE
          AB0=SB+SF*PFR  
          IF(POPI.GT.0.) DDN(N0HN)=DDN(N0HN)+(SB+SF*PFR)/POPI
          IF(MODE.EQ.1) THEN 
            DSF=SFFM2/T/FR
            DAB0=DAB0+SB*(TRHA+XHM)+DSF*PFR
          END IF
      END IF
      END IF

C
      IF(IOPH2P.GT.0) THEN
C
C   -----------------------------
C   H2+  bound-free and free-free
C   -----------------------------
C
      IT=IT+1
      IF(MODE.LT.0) THEN
         X=FR*1.D-15
         BFCS(IT,IJ)=real((-7.342D-3+(-2.409+(1.028+(-4.23D-1+
     *        (1.224D-1-1.351D-2*X)*X)*X)*X)*X)*1.602D-12/BOLK)
         IT=IT+1
         if(it.gt.mcross) 
     *    CALL QUIT('it.gt.mcross in opadd',it,mcross)
         if(ielh.le.0) 
     *    CALL QUIT('conflict in ioph2p - no hydrogen',ioph2p,ielh)
         X=LOG(FR)
         CS0=-3.0233D3+(3.7797D2+(-1.82496D1+(3.9207D-1-
     *        3.1672D-3*X)*X)*X)*X
         BFCS(IT,IJ)=real(cs0)
       ELSE
         X2=-CROSS(IT,IJ)/T+CROSS(IT+1,IJ)
         IT=IT+1
         SB=0.
         IF(X2.GT.-150..and.fr.lt.3.28e15)
     *      SB=POPUL(N0HN,ID)*EXP(X2)*POPUL(NKH,ID)
         AB1=SB
         DAB1=SB*CROSS(IT-1,IJ)/T/T
         DDN(N0HN)=DDN(N0HN)+SB/POPI
         IF(POPUL(NKH,ID).GT.0) DDN(NKH)=DDN(NKH)+SB/POPUL(NKH,ID)
      END IF
      END IF
C
C   ----------------------
C   neutral helium opacity
C   ----------------------
C     approximate, hydrogenic, opacity of neutral helium
C     given as a sum of bound-free transitions from averaged
C     levels with principal quantum numbers between that next
C     to the highest level considered explicitly and IOPHE1
C
      IF(IOPHE1.GT.0) THEN
      NHE2=NNEXT(IELHE1)
      I1=NQUANT(NLAST(IELHE1))+1
      IF(MODE.LT.0) THEN
         SG0=SG01/(FR*1.D-15)**3
         DO 20 I=I1,IOPHE1
            IT=IT+1
         if(it.gt.mcross) 
     *    CALL QUIT('it.gt.mcross in opadd',it,mcross)
            FRI=F0HE1/I/I
            SG=0.
            IF(FR.GE.FRI) SG=SG0*GAUNT(I,FR)/I**5
            BFCS(IT,IJ)=real(SG)
   20    CONTINUE
       ELSE
         SB=SBH0*T32
         DO 30 I=I1,IOPHE1
            IT=IT+1
            FRI=F0HE1/I/I
            IF(FR.LT.FRI) GO TO 30
            XI=HKT*FRI
            SG=SB*EXP(XI)*I*I*POPUL(NHE2,ID)*ANE*CROSS(IT,IJ)
            AB0=AB0+SG
            DAB0=DAB0+(TRHA+XI)*SG
            IF(POPUL(NHE2,ID).GT.0) DDN(NHE2)=DDN(NHE2)+
     *                             SG/POPUL(NHE2,ID)
   30    CONTINUE
      END IF
      END IF
C
C   --------------
C   ionized helium (analogously as for neutral helium)
C   --------------
C
      IF(IOPHE2.GT.0) THEN
      NHE3=NNEXT(IELHE2)
      I1=NQUANT(NLAST(IELHE2))+1
      IF(MODE.LT.0) THEN
         SG0=SG02/(FR*1.D-15)**3
         DO 40 I=I1,IOPHE2
            IT=IT+1
         if(it.gt.mcross) 
     *    CALL QUIT('it.gt.mcross in opadd',it,mcross)
            FRI=F0HE2/I/I
            SG=0.
            IF(FR.GE.FRI) SG=SG0*GAUNT(I,FR/4.D0)/I**5
            BFCS(IT,IJ)=real(SG)
   40    CONTINUE
       ELSE
         SB=SBH0*T32
         DO 50 I=I1,IOPHE2
            IT=IT+1
            FRI=F0HE2/I/I
            IF(FR.LT.FRI) GO TO 50
            XI=HKT*FRI
            SG=SB*EXP(XI)*I*I*POPUL(NHE3,ID)*ANE*CROSS(IT,IJ)
            AB0=AB0+SG
            DAB0=DAB0+(TRHA+XI)*SG
            IF(POPUL(NHE3,ID).GT.0) DDN(NHE3)=DDN(NHE3)+
     *                              SG/POPUL(NHE3,ID)
   50    CONTINUE
      END IF
      END IF
C
C     ----------------------------------------------
C     The user may supply more opacity sources here:
C     ----------------------------------------------
C
C     Finally, actual absorption and emission coefficients and
C     their derivatives

      IF(MODE.LT.0) RETURN
      X=EXP(-HKT*FR)
      X1=UN-X
      FR15=FR*1.D-15
      BNX=BN*FR15*FR15*FR15*X
      ABAD=ABAD+AB0+AB1
      EMAD=EMAD+AB0+AB1
      IF(MODE.EQ.1) THEN
         HKFT=HKT*FR/T
         DB=HKFT*AB0
         DAT=DAB1
         DET=DAB1
         DAN=AB0/ANE
         DEN=AB0/ANE
      END IF
      RETURN
      END
C
C
C     ****************************************************************
C
C
C

      SUBROUTINE PARTF(IAT,IZI,T,ANE,XMAX,U,DUT,DUN)
C     ==============================================
C
C     Partition functions 
C     The standard evaluation is for hydrogen through zinc, for
C     neutrals and first four ionization degrees.
C     Basically after Traving, Baschek, and Holweger, Abhand. Hamburg.
C     Sternwarte. Band VIII, Nr. 1 (1966)
C
C     For higher atomic numbers  modified Kurucz routine PFSAHA,
C     called PFHEAV here is used. The routine was provided by
C     Charles Proffitt.
C
C     The routine calls special procedures for Fe and Ni; or
C     the values based on the tabulated Opacity Project ionization
C     fractions
C
C     Input:
C      IAT   - atomic number
C      IZI   - ionic charge (=1 for neutrals, =2 for once ionized, etc)
C      T     - temperature
C      ANE   - electron density
C      XMAX  - principal quantum number of the last bound level
C
C     Output:
C      U     - partition function
C      DUT   - derivative dU/dT
C      DUN   - derivative dU/d(ANE)
C
C     Quantities in COMMON/PFSTDS:
C      for MODPF(IAT) < 0  - non-standard, user supplied procedure
C                            for evaluating partition functions PFSPEC
C                     = 0  - standard expressions
C                            After Traving, Baschek, and Holweger, Abhand.
C                            Hamburg. Sternwarte, Band VIII, Nr. 1 (1966)
C                     > 0  - partition functions evaluated from the
C                            Opacity Project ionization fractions
C                            (by routine OPFRAC)
C      PFSTD(IAT,IZI) - see above
C
      INCLUDE 'IMPLIC.FOR'
      COMMON/PFSTDS/PFSTD(99,30),MODPF(99)
      common/irwint/iirwin
      PARAMETER (NIONS=123, NSS=222)
      PARAMETER (UN=1.D0, HALF=0.5D0, TWO=2.D0, TRHA=1.5D0,
     *           THIRD=UN/3.D0, SIXTH=UN/6.D0)
      REAL*4 AHH( 6),  ALB(12),  AB (11),  AC (19),  AN (30),  AO (49),
     *       AF (34),  ANN(23),  ANA(19),  AMG(15),  AAL(17),  ASI(23),
     *       AP (19),  AS (29),  ACL(28),  AAR(25),  AK (30),  ACA(17),
     *       ASC(24),  ATI(33),  AV (33),  ACR(29),  AMN(28),  AFE(35),
     *       ACO(29),  ANI(23),  ACU(20),  AZN(18)
      REAL*4 GHH( 6),  GLB(12),  GB (11),  GC (19),  GN (30),  GO (49),
     *       GF (34),  GNN(23),  GNA(19),  GMG(15),  GAL(17),  GSI(23),
     *       GP (19),  GS (29),  GCL(28),  GAR(25),  GK (30),  GCA(17),
     *       GSC(24),  GTI(33),  GV (33),  GCR(29),  GMN(28),  GFE(35),
     *       GCO(29),  GNI(23),  GCU(20),  GZN(18)
      REAL*4 XL1(99), XL2(123),  XL(222),
     *       CH1(66),  CH2(72),  CH3(55),  CH4(29),  CHION(222)
      REAL*4 ALF(678), GAM(678)
C     INTEGER*2 II1(5,15),II2(5,15),INDEX0(5,30),
      INTEGER   II1(5,15),II2(5,15),INDEX0(5,30),
     *          IS1(53),IS2(70),IS(123),INDEXS(123),
     *          IM1(99),IM2(123),IM(222),INDEXM(222),
     *          IGP1(99),IGP2(123),IGPR(222),
     *          IG01(53),IG02(70),IG0(123)
      DIMENSION IGLE(28)
C
      EQUIVALENCE   ( AHH(1), ALF(  1)),( ALB(1), ALF(  7)),
     *              ( AB (1), ALF( 19)),
     *              ( AC (1), ALF( 30)),( AN (1), ALF( 49)),
     *              ( AO (1), ALF( 79)),( AF (1), ALF(128)),
     *              ( ANN(1), ALF(162)),( ANA(1), ALF(185)),
     *              ( AMG(1), ALF(204)),( AAL(1), ALF(219)),
     *              ( ASI(1), ALF(236)),( AP (1), ALF(259)),
     *              ( AS (1), ALF(278)),( ACL(1), ALF(307)),
     *              ( AAR(1), ALF(335)),( AK (1), ALF(360)),
     *              ( ACA(1), ALF(390)),( ASC(1), ALF(407)),
     *              ( ATI(1), ALF(431)),( AV (1), ALF(464)),
     *              ( ACR(1), ALF(497)),( AMN(1), ALF(526)),
     *              ( AFE(1), ALF(554)),( ACO(1), ALF(589)),
     *              ( ANI(1), ALF(618)),( ACU(1), ALF(641)),
     *              ( AZN(1), ALF(661))
      EQUIVALENCE   ( GHH(1), GAM(  1)),( GLB(1), GAM(  7)),
     *              ( GB (1), GAM( 19)),
     *              ( GC (1), GAM( 30)),( GN (1), GAM( 49)),
     *              ( GO (1), GAM( 79)),( GF (1), GAM(128)),
     *              ( GNN(1), GAM(162)),( GNA(1), GAM(185)),
     *              ( GMG(1), GAM(204)),( GAL(1), GAM(219)),
     *              ( GSI(1), GAM(236)),( GP (1), GAM(259)),
     *              ( GS (1), GAM(278)),( GCL(1), GAM(307)),
     *              ( GAR(1), GAM(335)),( GK (1), GAM(360)),
     *              ( GCA(1), GAM(390)),( GSC(1), GAM(407)),
     *              ( GTI(1), GAM(431)),( GV (1), GAM(464)),
     *              ( GCR(1), GAM(497)),( GMN(1), GAM(526)),
     *              ( GFE(1), GAM(554)),( GCO(1), GAM(589)),
     *              ( GNI(1), GAM(618)),( GCU(1), GAM(641)),
     *              ( GZN(1), GAM(661))
      EQUIVALENCE   ( CH1(1), CHION(  1)),
     *              ( CH2(1), CHION( 67)),
     *              ( CH3(1), CHION(139)),
     *              ( CH4(1), CHION(194)),
     *              ( XL1(1),    XL(  1)),
     *              ( XL2(1),    XL(100))
      EQUIVALENCE   ( IS1(1),  IS(1)),   ( IS2(1),  IS( 54)),
     *              ( IM1(1),  IM(1)),   ( IM2(1),  IM(100)),
     *              (IGP1(1),IGPR(1)),   (IGP2(1),IGPR(100)),
     *              (IG01(1), IG0(1)),   (IG02(1), IG0( 54)),
     *              (II1(1,1),INDEX0(1,1)),(II2(1,1),INDEX0(1,16))
C
      DATA IGLE/2,1,2,1,6,9,4,9,6,1,2,1,6,9,4,9,6,1,
     *          10,21,28,25,6,25,28,21,10,21/
C
C    Internal numbering of ions; each line corresponds to one species,
C    starting from hydrogen.
C    0 means that the ion does not exist
C    negative number means that partition function is assumed constant
C    (ie even for MODPF=0) and equal to ABS(that value)
C
      DATA II1      /   1,  -1,   0,   0,   0,
     *                  2,   3,  -1,   0,   0,
     *                  4,   5,  -2,  -1,   0,
     *                  6,   7,  -1,  -2,  -1,
     *                  8,   9,  10,  -1,  -2,
     *                 11,  12,  13,  14,  -1,
     *                 15,  16,  17,  18,  19,
     *                 20,  21,  22,  23,  24,
     *                 25,  26,  27,  28,  -6,
     *                 29,  30,  31,  32,  -9,
     *                 33,  34,  35,  36,  -4,
     *                 37,  38,  39,  40,  -9,
     *                 41,  42,  43,  44,  -6,
     *                 45,  46,  47,  48,  -1,
     *                 49,  50,  51,  52,  53                         /
      DATA II2      /  54,  55,  56,  57,  58,
     *                 59,  60,  61,  62,  63,
     *                 64,  65,  66,  67,  68,
     *                 69,  70,  71,  72,  73,
     *                 74,  75,  76,  77,  -9,
     *                 78,  76,  80,  81,  82,
     *                 83,  84,  85,  86,  87,
     *                 88,  89,  90,  91,  92,
     *                 93,  94,  95,  96,  97,
     *                 98,  99, 100, 101, 102,
     *                103, 104, 105, 106, 107,
     *                108, 109, 110, 111, -25,
     *                112, 113, 114, 115,  -1,
     *                116, 117, 118, 119,  -1,
     *                120, 121, 122, 123,  -1                         /
C
C     Statistical weights of the ground states
C
      DATA IG01     /   2,
     *                  1,   2,
     *                  2,   1,
     *                  1,   2,
     *                  2,   1,   2,
     *                  1,   2,   1,   2,
     *                  4,   1,   2,   1,   2,
     *                  5,   4,   1,   2,   1,
     *                  4,   5,   4,   1,
     *                  1,   4,   5,   4,
     *                  2,   1,   4,   5,
     *                  1,   2,   1,   4,
     *                  2,   1,   2,   1,
     *                  1,   2,   1,   2,
     *                  4,   1,   2,   1,   2                         /
      DATA  IG02    /   5,   4,   1,   2,   1,
     *                  4,   5,   4,   1,   2,
     *                  1,   4,   5,   4,   1,
     *                  2,   1,   4,   5,   4,
     *                  1,   2,   1,   4,
     *                  4,   3,   4,   1,   4,
     *                  5,   4,   5,   4,   1,
     *                  4,   1,   4,   5,   4,
     *                  7,   6,   1,   4,   5,
     *                  6,   7,   6,   1,   4,
     *                  9,  10,   9,   6,   1,
     *                 10,   9,  10,  20,
     *                  9,   6,   9,  28,
     *                  2,   1,   6,  21,
     *                  1,   2,   1,  10                              /
C
      DATA IS1      /   1,
     *                  1,   1,
     *                  1,   1,
     *                  2,   1,
     *                  1,   2,   1,
     *                  1,   2,   2,   1,
     *                  2,   2,   3,   2,   1,
     *                  3,   4,   3,   5,   2,
     *                  2,   3,   4,   3,
     *                  2,   2,   3,   2,
     *                  1,   2,   2,   3,
     *                  1,   1,   2,   2,
     *                  2,   2,   1,   2,
     *                  1,   2,   2,   1,
     *                  2,   1,   1,   1,   1                         /
      DATA  IS2     /   3,   2,   1,   2,   2,
     *                  2,   3,   2,   1,   1,
     *                  2,   2,   3,   1,   1,
     *                  1,   2,   3,   3,   2,
     *                  2,   1,   2,   2,
     *                  3,   1,   1,   1,   1,
     *                  3,   2,   1,   1,   1,
     *                  2,   3,   1,   1,   1,
     *                  3,   2,   1,   1,   1,
     *                  3,   2,   1,   1,   1,
     *                  3,   2,   2,   1,   1,
     *                  4,   2,   1,   1,
     *                  2,   2,   1,   1,
     *                  3,   2,   1,   1,
     *                  3,   3,   1,   1                              /
C
      DATA IM1      /   2,
     *                  2,   2,
     *                  2,   2,
     *                  3,   2,   3,
     *                  3,   3,   2,   3,
     *                  4,   3,   3,   3,   3,   3,
     *                  3,   3,   4,   3,   3,   4,   2,   3,   2,   3,
     *                  4,   2,   2,   4,   2,   3,   3,   4,   4,   2,
     *                  3,   4,   2,   2,   2,   3,   3,
     *                  3,   3,   4,   2,   2,
     *                  4,   2,   3,   2,   5,   2,   2,
     *                  2,   2,   3,   2,   4,   2,   2,   4,   2,
     *                  2,   2,   2,   3,   2,   4,   2,   2,
     *                  3,   3,   2,   2,   3,   2,
     *                  3,   2,   3,   2,   3,   2,   2,
     *                  5,   4,   4,   4,   3,   3,
     *                  3,   2,   4,   4,   3,   3                    /
      DATA  IM2     /   4,   2,   2,   4,   2,   5,   4,   2,   3,   1,
     *                  3,   2,   5,   2,   2,   4,   2,   4,   4,
     *                  2,   2,   3,   2,   4,   2,   2,   4,   4,
     *                  3,   2,   3,   3,   2,   3,
     *                  4,   2,   2,   4,   2,
     *                  3,   2,   3,   2,   2,   3,   2,
     *                  4,   3,   3,   5,   4,   2,   3,
     *                  6,   4,   3,   6,   3,   5,   4,   2,
     *                  5,   3,   5,   4,   4,   4,   4,   4,
     *                  3,   3,   3,   4,   4,   4,   4,   4,
     *                  3,   2,   3,   4,   4,   4,   4,   4,
     *                  4,   4,   3,   5,   3,   4,   4,   4,   4,
     *                  5,   3,   3,   3,   5,   4,   5,   1,
     *                  6,   3,   5,   3,   5,   1,
     *                  2,   3,   3,   4,   3,   4,   1,
     *                  2,   2,   2,   3,   3,   2,   3,   1          /
C
      DATA IGP1     /   2,
     *                  4,   2,
     *                  2,   4,
     *                  4,  12,   2,
     *                  2,   4,  12,   2,
     *                 12,   2,  18,   4,  12,   2,
     *                 18,  10,  12,  24,   2,  18,   6,   4,  12,   2,
     *                  8,  20,  12,  18,  10,   2,  10,  12,  24,  20,
     *                  2,  18,   6,  18,  10,   4,  12,
     *                 18,  10,   8,  20,  12,
     *                 18,  10,   2,  10,  12,  24,  20,
     *                  8,   4,  18,  10,   8,  20,  12,  18,  10,
     *                  2,   8,   4,  18,  10,   8,  20,  12,
     *                  4,   2,   8,   4,  18,  10,
     *                  2,  18,   4,  12,   2,   8,   4,
     *                 12,   2,  18,   4,  12,   2,
     *                 18,  10,  12,   2,   4,   2                    /
      DATA  IGP2    /   8,  20,  12,  18,  10,  12,   2,  18,   4,  12,
     *                 18,  10,   8,  20,  12,  18,  10,  12,   2,
     *                  8,   4,  18,  10,   8,  20,  12,  18,  12,
     *                  2,   8,   4,  18,  10,   2,
     *                  8,  20,  12,  18,  10,
     *                  4,  20,   2,   8,   4,  18,  10,
     *                 30,  42,  18,  20,   2,  12,  18,
     *                 56,  56,  28,  42,  10,  20,   2,  12,
     *                 50,  70,  56,  72,  64,  42,  20,   2,
     *                 12,  60,  40,  50,  18,  56,  42,  20,
     *                 14,  10,  50,  12,  72,  50,  56,  42,
     *                 60,  56,  40,  50,  18,  12,  72,  50,  56,
     *                 42,  70,  42,  18,  56,  24,  50,  12,
     *                 20,  56,  42,  18,  56,  50,
     *                  2,  30,  10,  20,  56,  42,  56,
     *                  4,   8,  12,   2,  30,  10,  20,  42          /
C
      DATA XL1      /11.0,
     *                8.0,12.0,
     *                6.0, 6.0,
     *                6.0, 4.0, 8.0,
     *                9.0, 6.0, 4.0, 6.0,
     *                6.0, 6.0, 5.0, 6.1, 5.0, 6.0,
     *                6.1, 4.0, 5.0, 3.9, 6.0, 5.0, 4.0, 6.0, 6.3, 6.0,
     *                8.0, 6.0, 3.4, 6.0, 5.0, 3.9, 3.9, 6.0, 4.9, 4.0,
     *                5.9, 5.0, 4.9, 4.0, 4.0, 6.0, 6.0,
     *                4.0, 4.0, 5.0, 4.0, 4.0,
     *                5.0, 4.0, 3.9, 4.0, 5.0, 5.0, 4.0,
     *                6.0, 6.0, 5.0, 4.0, 3.9, 4.0, 4.0, 5.0, 5.0,
     *                7.0, 4.0, 4.0, 4.0, 4.0, 5.0, 5.0, 5.0,
     *                7.0, 7.0, 5.0, 5.0, 5.0, 5.0,
     *                7.0, 4.0, 7.0, 4.0, 7.0, 5.0, 5.0,
     *                6.1, 5.9, 5.0, 5.0, 5.0, 7.0,
     *                5.0, 5.0, 5.0, 7.0, 8.6, 8.0                    /
      DATA  XL2     / 6.0, 5.0, 5.0, 5.0, 5.0, 3.5, 5.0,14.4, 5.0, 4.0,
     *                6.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.2,
     *                6.0, 6.0, 5.1, 5.0, 5.0, 5.0, 5.0, 5.0, 4.0,
     *                7.0, 5.0, 5.0, 6.0, 6.0, 5.0,
     *                6.0, 5.0, 5.0, 3.6, 4.0,
     *                5.9, 6.0, 7.0, 5.0, 4.9, 5.0, 4.3,
     *                4.9, 4.9, 5.0, 5.0, 6.0, 4.6, 3.8,
     *                5.0, 4.7, 5.0, 5.0, 5.0, 5.0, 6.0, 4.8,
     *                5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0,11.2,
     *                5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.2,
     *                6.0, 5.0, 6.0, 7.0, 5.0, 5.0, 5.0, 5.0,
     *                5.0, 5.0, 5.0, 5.0, 5.0, 6.0, 5.0, 3.6, 3.8,
     *                5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 3.0,
     *                5.4, 5.0, 9.0, 5.0, 5.0, 3.0,
     *                8.0, 6.0, 5.0, 7.0, 5.0, 5.0, 2.9,
     *                8.0, 5.0, 5.0, 8.0, 5.0, 5.0, 5.0, 2.8          /
C
C
      DATA CH1      /  13.595 ,
     *                 24.580 ,  54.403 ,
     *                  5.390 , 75.619 ,
     *                  9.320 ,  13.278 ,  18.206 ,
     *                  8.296 ,  25.149 ,  31.146 ,  37.920 ,
     *                 11.256 ,  24.376 ,  30.868 ,  47.871 ,  55.873 ,
     *                 64.476 ,
     *                 14.529 ,  16.428 ,  29.593 ,  36.693 ,
     *                 47.426 ,  55.765 ,  63.626 ,  77.450 ,  87.445 ,
     *                 97.863 ,
     *                 13.614 ,  16.938 ,  18.630 ,
     *                 35.108 ,  37.621 ,  40.461 ,  42.584 ,
     *                 54.886 ,  63.733 ,  70.556 ,
     *                 77.394 ,  87.609 ,  97.077 , 103.911 , 106.116 ,
     *                113.873 , 125.863 ,
     *                 17.418 ,  20.009 ,  34.977 ,  39.204 ,  41.368 ,
     *                 62.646 ,  65.774 ,  69.282 ,  71.882 ,
     *                 87.139 ,  97.852 , 106.089 ,
     *                 21.559 ,  21.656 ,  41.071 ,  44.274 ,
     *                 63.729 ,  68.806 ,  71.434 ,  97.162 , 100.917 /
      DATA CH2      /   5.138 ,  47.290 ,  47.459 ,  71.647 ,  75.504 ,
     *                 98.880 , 104.778 , 107.864 ,
     *                  7.644 ,  15.031 ,  80.117 ,  80.393 ,
     *                109.294 , 113.799 ,
     *                  5.984 ,  10.634 ,  18.823 ,  25.496 ,
     *                 28.441 , 119.957 , 120.383 ,
     *                  8.149 ,  16.339 ,  22.894 ,
     *                 33.459 ,  42.333 ,  45.130 ,
     *                 10.474 ,  11.585 ,  19.720 ,
     *                 30.156 ,  51.354 ,  65.007 ,
     *                 10.357 ,  12.200 ,  13.401 ,  23.405 ,  24.807 ,
     *                 35.047 ,  47.292 ,  57.681 ,  72.474 ,  85.701 ,
     *                 13.014 ,  14.458 ,  23.798 ,  26.041 ,  27.501 ,
     *                 39.904 ,  41.610 ,  53.450 ,  67.801 ,
     *                 15.755 ,  15.933 ,  27.619 ,  29.355 ,
     *                 40.899 ,  42.407 ,  45.234 ,  59.793 ,  75.002 ,
     *                  4.339 ,  31.810 ,  32.079 ,
     *                 45.738 ,  47.768 ,  50.515 ,
     *                 60.897 ,  63.890 ,  65.849 ,  82.799 ,  85.150 /
      DATA CH3      /   6.111 ,   7.808 ,  11.868 ,
     *                 51.207 ,  51.596 ,  67.181 ,  69.536 ,
     *                  6.538 ,   7.147 ,   8.042 ,
     *                 12.891 ,  24.752 ,  74.090 ,  91.847 ,
     *                  6.818 ,  6.953 ,  7.411 ,
     *                 13.635 ,  14.685 ,  28.137 ,  43.236 , 100.083 ,
     *                  6.738 ,   7.101 ,  14.205 ,  15.670 ,  16.277 ,
     *                 29.748 ,  48.464 ,  65.198 ,
     *                  6.763 ,   8.285 ,   9.221 ,
     *                 16.493 ,  18.662 ,  30.950 ,  49.580 ,  73.093 ,
     *                  7.432 ,   8.606 ,   9.240 ,  15.636 ,  18.963 ,
     *                 33.690 ,  53.001 ,  76.006 ,
     *                  7.896 ,   8.195 ,   8.927 ,  16.178 ,  18.662 ,
     *                 30.640 ,  34.607 ,  56.001 ,  79.001           /
      DATA CH4      /   7.863 ,   8.378 ,   9.160 ,   9.519 ,
     *                 17.052 ,  18.958 ,  33.491 ,  53.001 ,
     *                  7.633 ,   8.793 ,  18.147 ,  20.233 ,  35.161 ,
     *                 56.025 ,
     *                  7.724 ,  10.532 ,  10.980 ,
     *                 20.286 ,  27.985 ,  36.826 ,  61.975 ,
     *                  9.391 ,  17.503 ,  17.166 ,
     *                 17.959 ,  27.757 ,  28.310 ,  39.701 ,  65.074 /
C
C     data for hydrogen and helium
C
      DATA AHH      /  20.4976, 747.5023,
     *                 28.1703, 527.8296,  22.2809, 987.7189          /
      DATA GHH      /  10.853 ,  13.342 ,
     *                 21.170 ,  24.125 ,  43.708 ,  53.542           /
C
C     data for lithium and beryllium
C
      DATA ALB      /   8.4915,  97.5015,  23.3299, 192.6701,
     *                  9.1849,  32.9263, 183.8887,  19.9563,  88.0437,
     *                  6.0478,  35.9723, 233.9798                    /
      DATA GLB      /   2.022 ,   4.604 ,  62.032 ,  72.624 ,
     *                  2.735 ,   6.774 ,   8.569 ,  10.750 ,  11.672 ,
     *                  3.967 ,  12.758 ,  16.692                     /
C
C     data for boron
C
      DATA AB       /   4.0086,  19.6741, 402.3110,
     *                  9.7257,  30.9262, 186.3466,  44.1629,  60.8371,
     *                  6.0084,  23.5767,  76.4149                    /
      DATA GB       /   0.002 ,   3.971 ,   7.882 ,
     *                  4.720 ,  13.477 ,  22.103 ,  23.056 ,  24.734 ,
     *                  6.000 ,  24.540 ,  32.300                     /
C
C     data for carbon
C
      DATA AC       /   8.0158,   5.8833,  33.7521, 595.3432,
     *                  4.0003,  17.0841,  82.9154,
     *                 15.9808,  48.2044, 435.8093,
     *                 10.0281,  15.7574, 186.2109,
     *                 15.4127,  55.9559, 243.6311,
     *                  6.0057,  23.5757,  76.4185                    /
      DATA GC       /   0.004 ,   1.359 ,   6.454 ,  10.376 ,
     *                  0.008 ,  16.546 ,  21.614 ,
     *                  5.688 ,  15.801 ,  26.269 ,
     *                  6.691 ,  25.034 ,  40.975 ,
     *                 17.604 ,  36.180 ,  47.133 ,
     *                  8.005 ,  40.804 ,  54.492                     /
C
C     data for nitrogen
C
      DATA AN       /  14.0499,  30.8008, 883.1443,
     *                 10.0000,  16.0000,  64.0000,
     *                  8.0462,   6.2669,  17.8696, 282.8084,
     *                  7.3751,  33.1390, 215.4829,
     *                  4.0003,  19.3533,  80.6462,
     *                 13.0998,  19.6425,  94.3035, 370.9539,
     *                 16.0000,  38.0000,
     *                 10.3289,  14.5021, 187.1624, 108.1615, 191.8383,
     *                  6.0044,  23.5612,  76.4344                    /
      DATA GN       /   2.554 ,   9.169 ,  13.651 ,
     *                 12.353 ,  13.784 ,  14.874 ,
     *                  0.014 ,   2.131 ,  15.745 ,  24.949 ,
     *                  6.376 ,  14.246 ,  29.465 ,
     *                  0.022 ,  31.259 ,  41.428 ,
     *                  7.212 ,  15.228 ,  34.387 ,  46.708 ,
     *                 46.475 ,  49.468 ,
     *                  8.693 ,  37.650 ,  65.479 ,  61.155 ,  79.196 ,
     *                  9.999 ,  60.991 ,  82.262                     /
C
C     data for oxygen
C
      DATA AO       /   4.0029,   5.3656,  36.2853,1044.3447,
     *                131.0217, 868.9779,  14.8533,  93.1466,
     *                 12.7843,   5.6828,  98.0919, 829.4396,
     *                 50.9878, 199.0120,   2.0000,   6.0000,  10.0000,
     *                 10.0000,  30.0000,  50.0000,
     *                  8.0703,   5.7144,  84.1156, 529.0927,
     *                  5.6609,  28.9355, 111.3620, 494.0413,
     *                 45.5249, 134.4751,
     *                  4.0003,  21.2937,  78.7058,
     *                 12.8293,  16.2730, 123.6578, 327.2396,
     *                 48.7883, 102.2117,  20.0060, 161.9903,
     *                 28.4184,  61.5816,
     *                 10.5563,  13.2950, 188.1390,
     *                 14.6560, 129.4922, 470.8512                    /
      DATA GO       /   0.022 ,   2.019 ,   9.812 ,  13.087 ,
     *                 13.804 ,  16.061 ,  14.293 ,  16.114 ,
     *                  3.472 ,   7.437 ,  22.579 ,  32.035 ,
     *                 27.774 ,  33.678 ,  28.118 ,  31.019 ,  34.204 ,
     *                 30.892 ,  33.189 ,  36.181 ,
     *                  0.032 ,   2.760 ,  35.328 ,  48.277 ,
     *                  7.662 ,  16.786 ,  42.657 ,  54.522 ,
     *                 50.204 ,  56.044 ,
     *                  0.048 ,  50.089 ,  66.604 ,
     *                  8.954 ,  18.031 ,  57.755 ,  72.594 ,
     *                 68.388 ,  82.397 ,  31.960 ,  76.876 ,
     *                 75.686 ,  80.388 ,
     *                 10.747 ,  52.323 ,  94.976 ,
     *                 27.405 ,  86.350 , 109.917                     /
C
C     data for fluor
C
      DATA AF       /   2.0001,  39.9012, 122.0986,
     *                 10.0000,  30.0000,  50.0000,
     *                  4.0199,   5.5741,  22.1839, 190.2179,
     *                 53.0383, 126.9616,  31.6894,  75.3105,
     *                 13.5014,   7.9936,  55.7981, 298.7039,
     *                 26.2496,  63.7503,   2.0000,   6.0000,  10.0000,
     *                 28.7150,  71.2850,
     *                  8.0153,   6.1931,  21.7287,  48.7780, 278.2782,
     *                178.5560, 421.4435,  51.7632,  95.2368          /
      DATA GF       /   0.050 ,  13.317 ,  15.692 ,
     *                 15.361 ,  17.128 ,  18.498 ,
     *                  0.048 ,   2.735 ,  20.079 ,  30.277 ,
     *                 27.548 ,  32.532 ,  30.391 ,  34.707 ,
     *                  4.479 ,  12.072 ,  31.662 ,  51.432 ,
     *                 44.283 ,  50.964 ,  46.193 ,  50.436 ,  54.880 ,
     *                 50.816 ,  57.479 ,
     *                  0.058 ,   3.434 ,  14.892 ,  37.472 ,  69.883 ,
     *                 67.810 ,  83.105 ,  72.435 ,  79.747           /
C
C     data for neon
C
      DATA ANN      /  34.5080, 365.4919,  16.5768, 183.4231,
     *                  2.0007,  89.5607, 380.4381,  26.4473,  63.5527,
     *                  4.0342,   5.6162,  11.5176,  72.8273,
     *                 48.5684, 131.4315,  31.1710,  76.8290,
     *                 14.0482,  13.3077,  52.7897, 467.8487,
     *                 54.2196, 195.7800                              /
      DATA GNN      /  17.796 ,  20.730 ,  17.879 ,  20.855 ,
     *                  0.097 ,  29.878 ,  37.221 ,  31.913 ,  37.551 ,
     *                  0.092 ,   3.424 ,  24.806 ,  46.616 ,
     *                 45.643 ,  54.147 ,  48.359 ,  57.420 ,
     *                  5.453 ,  18.560 ,  46.583 ,  80.101 ,
     *                 70.337 ,  85.789                               /
C
C     data for sodium
C
      DATA ANA      /  11.6348, 158.3593,
     *                 21.0453,  50.9546,  10.1389,  25.8611,
     *                  2.0019,  38.0569, 137.9398,  28.3106,  61.6893,
     *                  4.0334,   5.8560,  18.1786, 208.9142,
     *                 93.6895, 406.3095,  60.4276, 239.5719          /
      DATA GNA      /   2.400 ,   4.552 ,
     *                 34.367 ,  40.566 ,  34.676 ,  40.764 ,
     *                  0.170 ,  44.554 ,  57.142 ,  51.689 ,  60.576 ,
     *                  0.152 ,   4.260 ,  36.635 ,  83.254 ,
     *                 72.561 ,  89.475 ,  75.839 ,  92.582           /
C
C     data for magnesium
C
      DATA AMG      /  10.7445, 291.5057,  53.7488,
     *                  6.2270,  31.1291, 132.6438,
     *                 40.4379, 159.5618,  20.3845,  79.6154,
     *                  2.0007, 106.8977, 343.1010,  10.1326, 237.8581/
      DATA GMG      /   2.805 ,   6.777 ,   9.254 ,
     *                  4.459 ,   9.789 ,  13.137 ,
     *                 57.413 ,  71.252 ,  58.010 ,  71.660 ,
     *                  0.276 ,  74.440 ,  94.447 ,  54.472 ,  95.858 /
C
C     data for aluminium
C
      DATA AAL      /   4.0009,  11.7804, 142.2179,  13.6585,  96.3371,
     *                 10.0807,  49.5843, 285.3343,  14.6872,  59.3122,
     *                  6.3277,  29.5086, 134.1634,
     *                 46.3164, 153.6833,  22.9896,  77.0103          /
      DATA GAL      /   0.014 ,   3.841 ,   5.420 ,   3.727 ,   8.833 ,
     *                  4.749 ,  11.902 ,  16.719 ,  11.310 ,  18.268 ,
     *                  6.751 ,  16.681 ,  24.151 ,
     *                 83.551 , 104.787 ,  84.293 , 105.171           /
C
C     data for silicon
C
      DATA ASI      /   7.9658,   4.6762,   1.3512, 123.2267, 443.7797,
     *                  4.0000,   7.4186,  24.1754, 60.4060,
     *                 14.4695,  11.9721,  26.5062, 269.0521,
     *                  9.1793,   4.8766,  29.1442,  52.7998,
     *                 13.2674,  36.0417, 180.6910,
     *                  6.4839,  27.6851, 135.8301                    /
      DATA GSI      /   0.020 ,   0.752 ,   1.614 ,   5.831 ,   7.431 ,
     *                  0.036 ,   8.795 ,  11.208 ,  13.835 ,
     *                  5.418 ,   7.825 ,  14.440 ,  19.412 ,
     *                  6.572 ,  11.449 ,  18.424 ,  25.457 ,
     *                 15.682 ,  27.010 ,  34.599 ,
     *                  9.042 ,  24.101 ,  37.445                     /
C
C     data for phosphorus
C
      DATA AP       /  13.5211,  22.2130, 353.2583,  10.0000, 150.0000,
     *                  8.0241,   5.8085,  51.7542, 252.4002,
     *                  4.0021,  20.7985,  62.4194, 200.7786,
     *                 11.7414,  63.5124, 179.7420,
     *                  6.8835,  32.7777, 228.3366                    /
      DATA GP       /   1.514 ,   5.575 ,   9.247 ,   8.076 ,  10.735 ,
     *                  0.043 ,   1.212 ,   8.545 ,  15.525 ,
     *                  0.074 ,   7.674 ,  16.639 ,  25.118 ,
     *                  8.992 ,  24.473 ,  40.704 ,
     *                 11.464 ,  33.732 ,  55.455                     /
C
C     data for sulphur
C
      DATA AS       /   3.9615,   5.0780,  15.0944, 362.8588,
     *                 51.5995, 268.4002,  12.0000, 276.0000,
     *                 11.4377,   5.5126, 141.0009, 254.0478,
     *                 33.0518, 126.9479,
     *                  4.0707,   4.0637,   5.7245, 144.6376, 106.4909,
     *                  4.0011,  19.2813,  27.5990,  35.1179,
     *                 94.7454, 283.2486,
     *                 10.5474,  28.7137,  65.7378,  24.0000          /
      DATA GS       /   0.053 ,   1.121 ,   5.812 ,   9.425 ,
     *                  8.936 ,  11.277 ,   9.600 ,  12.551 ,
     *                  1.892 ,   3.646 ,  13.550 ,  19.376 ,
     *                 16.253 ,  21.062 ,
     *                  0.043 ,   0.123 ,   1.590 ,  13.712 ,  22.050 ,
     *                  0.118 ,   9.545 ,  18.179 ,  31.441 ,
     *                 30.664 ,  56.150 ,
     *                 10.704 ,  27.075 ,  50.599 ,  43.034           /
C
C     data for chlorine
C
      DATA ACL      /   2.0007,  62.5048, 669.4942,  29.0259, 130.9740,
     *                  3.9064,   0.3993,   5.3570,  60.3424, 119.9913,
     *                138.1567, 278.8418, 102.3681, 158.6314,
     *                 12.6089,   5.9527, 110.5635, 262.8715,
     *                 69.2035, 100.7960,
     *                  7.3458,   5.6638,  44.1256, 202.7846,
     *                  4.0037,  21.8663,  40.5363,  57.5919          /
      DATA GCL      /   0.110 ,   9.919 ,  12.280 ,  11.017 ,  13.532 ,
     *                  0.092 ,   0.581 ,   1.620 ,  13.121 ,  19.787 ,
     *                 16.365 ,  21.988 ,  18.065 ,  23.594 ,
     *                  2.358 ,   5.708 ,  19.084 ,  30.683 ,
     *                 24.880 ,  33.229 ,
     *                  0.102 ,   1.391 ,  14.709 ,  36.968 ,
     *                  0.185 ,  11.783 ,  25.653 ,  44.698           /
C
C     data for argon
C
      DATA AAR      /  43.6623, 324.3375,  20.8298, 163.1701,
     *                  2.0026, 137.4515, 258.5445,  62.8129, 149.1867,
     *                  4.0495,  14.4466,  46.8234, 124.6651,
     *                151.9828, 268.0157, 101.1302, 150.8691,
     *                 13.3718,   8.6528,  60.4614, 285.5072,
     *                  6.7655,   4.7684,  12.8631,  54.5260          /
      DATA GAR      /  12.638 ,  14.958 ,  12.833 ,  15.139 ,
     *                  0.178 ,  17.522 ,  23.584 ,  20.464 ,  25.150 ,
     *                  0.151 ,   1.561 ,  17.399 ,  30.871 ,
     *                 24.684 ,  33.978 ,  27.091 ,  36.481 ,
     *                  2.810 ,   8.877 ,  24.351 ,  44.489 ,
     *                  0.144 ,   1.160 ,  10.210 ,  27.178           /
C
C     data for potassium
C
      DATA AK       /  12.9782, 148.6673,   6.3493,
     *                 66.3444, 101.6553,   4.0001,  13.4465,  46.5534,
     *                  2.0171, 116.4767, 713.4965,  63.5907, 396.4079,
     *                  2.0000,  10.0000,  30.0000,
     *                  4.0702,   5.7791,  52.6795, 327.4539,
     *                 62.8604, 357.1331,  55.9337, 196.0646,
     *                 10.9275,   5.5398,  43.2761,  76.2560,
     *                 42.0000,  18.0000                              /
      DATA GK       /   1.871 ,   3.713 ,  18.172 ,
     *                 21.185 ,  27.705 ,   2.059 ,  23.709 ,  28.542 ,
     *                  0.273 ,  26.709 ,  39.640 ,  31.220 ,  41.865 ,
     *                 29.955 ,  37.557 ,  42.862 ,
     *                  0.228 ,   2.274 ,  21.703 ,  50.191 ,
     *                 32.145 ,  49.262 ,  34.155 ,  51.718 ,
     *                  3.043 ,   5.479 ,  20.547 ,  30.680 ,
     *                 36.275 ,  47.345                               /
C
C     data for calcium
C
      DATA ACA      /  18.2366,  27.5012, 149.2617,  94.5242, 705.4711,
     *                 11.8706,  14.0710, 106.0547,
     *                 57.2414, 110.7567,  29.8121,  54.1874,
     *                  2.0184,  97.5784, 282.3939, 209.1871, 252.8129/
      DATA GCA      /   2.050 ,   3.349 ,   5.321 ,   4.873 ,   7.017 ,
     *                  1.769 ,   5.109 ,   9.524 ,
     *                 27.271 ,  41.561 ,  29.172 ,  42.140 ,
     *                  0.394 ,  28.930 ,  52.618 ,  38.593 ,  49.646 /
C
C     data for scandium
C
      DATA ASC      /   6.0014,  83.1958,  67.3666, 329.4354,
     *                 44.0793, 169.9969, 533.9195,
     *                 34.1642, 124.8475, 228.9879,
     *                 11.9979,  16.9280,  28.4778,  82.0418, 234.5360,
     *                  6.0042,   2.7101,  13.9801,  65.3039,
     *                 12.0000,  12.0000,
     *                  2.0051,   2.9621,  29.0306                    /
      DATA GSC      /   0.021 ,   2.056 ,   3.551 ,   5.465 ,
     *                  1.535 ,   3.797 ,   6.203 ,
     *                  2.389 ,   4.858 ,   7.141 ,
     *                  0.011 ,   0.430 ,   1.156 ,   3.711 ,   8.863 ,
     *                  0.025 ,   3.499 ,  10.463 ,  18.606 ,
     *                 41.779 ,  57.217 ,
     *                  0.539 ,  24.442 ,  51.079                     /
C
C     data for titanium
C
      DATA ATI      /   7.0887,   8.9186,  17.5633, 206.6832, 438.5735,
     *                654.1721,
     *                 38.0462,  69.6271, 364.2845, 832.0408,
     *                 98.8562,  57.9934, 442.1498,
     *                 19.7843,  32.0637,  37.0895, 110.6682, 288.4946,
     *                521.8837,
     *                 10.0000,  34.0000, 120.0000,
     *                 16.1691,  22.3550,  24.1646,  83.5128, 222.7963,
     *                  6.0020,   4.6177,  25.2636,  52.1162,
     *                 12.0000,   8.0000                              /
      DATA GTI      /   0.021 ,   0.048 ,   1.029 ,   2.183 ,   4.109 ,
     *                  5.785 ,
     *                  0.846 ,   1.792 ,   3.836 ,   5.787 ,
     *                  2.561 ,   4.869 ,   6.340 ,
     *                  0.023 ,   0.124 ,   0.774 ,   1.810 ,   4.980 ,
     *                  9.585 ,
     *                  1.082 ,   4.928 ,  11.279 ,
     *                  0.041 ,   1.375 ,   4.768 ,  10.985 ,  19.769 ,
     *                  0.048 ,  11.577 ,  24.531 ,  36.489 ,
     *                 54.436 ,  75.373                               /
C
C     data for vanadium
C
      DATA AV       /  15.2627,  23.9869,  51.3053, 570.3384,1650.9417,
     *                162.2829, 298.8303, 908.8852,
     *                 23.6736,  37.1624,  86.8011, 300.7440, 864.5880,
     *                 57.8961,  79.4605, 214.9007, 864.7425,
     *                 61.8508,  64.0845, 192.8298, 718.2349,
     *                 23.8116,  68.2495, 135.0613, 536.7632,
     *                 15.9543,  22.5542,  71.4921, 248.9544,
     *                  6.0006,   5.8785,  50.5077,  97.6129          /
      DATA GV       /   0.026 ,   0.145 ,   0.718 ,   2.586 ,   5.458 ,
     *                  2.171 ,   4.153 ,   6.097 ,
     *                  0.009 ,   0.366 ,   1.504 ,   5.294 ,  10.126 ,
     *                  1.796 ,   2.353 ,   6.068 ,  12.269 ,
     *                  2.560 ,   3.674 ,   6.593 ,  12.880 ,
     *                  0.045 ,   1.684 ,   8.162 ,  21.262 ,
     *                  0.065 ,   1.746 ,  15.158 ,  33.141 ,
     *                  0.077 ,  21.229 ,  44.134 ,  60.203           /
C
C     data for chromium
C
      DATA ACR      /  30.1842,  79.2847, 149.5293,
     *                215.3696, 119.1974, 741.4321,
     *                184.9946,1352.5038, 784.4937,
     *                 46.6191, 160.1361, 488.0449, 657.1928,
     *                 47.1742, 267.0275, 441.1324, 150.6650,
     *                 24.3768, 122.8359, 285.5092, 794.1654,
     *                 24.2296,  75.0258, 172.9452, 543.6511,
     *                 15.9819,  17.6800,  95.2003, 225.0947          /
      DATA GCR      /   0.993 ,   3.070 ,   5.673 ,
     *                  3.339 ,   4.801 ,   7.198 ,
     *                  2.829 ,   4.990 ,   7.643 ,
     *                  1.645 ,   3.727 ,   7.181 ,  12.299 ,
     *                  2.902 ,   4.273 ,   8.569 ,  14.912 ,
     *                  0.047 ,   2.566 ,   9.441 ,  21.198 ,
     *                  0.078 ,   2.242 ,  15.638 ,  32.725 ,
     *                  0.103 ,   2.146 ,  26.153 ,  49.381           /
C
C     data for manganese
C
      DATA AMN      /  53.9107,  81.3931, 546.6945 ,
     *                144.1893, 407.8029,  45.6177, 298.4423,2410.9335,
     *                 22.6382,  93.8419, 183.9367, 907.5765,
     *                137.0409, 168.6783, 329.0287, 773.2513,
     *                 70.1925,  72.3372, 213.9512, 539.5165,
     *                 24.2373,  93.5415, 456.6167, 506.5484,
     *                 24.7687,  66.9896, 264.1853, 484.0161          /
      DATA GMN      /   2.527 ,   4.204 ,   6.602 ,
     *                  4.155 ,   7.321 ,   2.285 ,   5.631 ,   8.448 ,
     *                  1.496 ,   3.839 ,   7.751 ,  13.484 ,
     *                  3.681 ,   6.054 ,   9.934 ,  14.936 ,
     *                  3.531 ,   6.967 ,  15.222 ,  25.069 ,
     *                  0.071 ,   2.896 ,  20.725 ,  37.383 ,
     *                  0.126 ,   2.660 ,  28.528 ,  53.413           /
C
C     data for iron
C
      DATA AFE      /  14.4102,   2.7050, 421.6612, 940.1484,
     *                 36.2187,  22.8883, 239.5997, 825.2919,
     *                110.0242, 992.3040, 640.6715,
     *                 17.0494,  32.3783,  34.3184, 420.9626,1067.2064,
     *                154.0059, 462.1117, 329.8618,
     *                 15.7906,  47.1186, 279.9292, 692.1005,
     *                 91.0206, 206.3082, 706.9927, 836.6689,
     *                 40.0790,  27.6965,  28.2243,  18.0001,
     *                 24.0899,  89.6340,  51.5756, 241.6980          /
      DATA GFE      /   0.066 ,   0.339 ,   2.897 ,   6.585 ,
     *                  0.923 ,   1.679 ,   4.620 ,   7.053 ,
     *                  4.249 ,   5.875 ,   7.781 ,
     *                  0.062 ,   0.283 ,   1.504 ,   5.430 ,  11.210 ,
     *                  2.792 ,   7.627 ,  13.623 ,
     *                  0.077 ,   3.723 ,  12.137 ,  23.700 ,
     *                  2.688 ,   7.595 ,  15.444 ,  25.587 ,
     *                  3.982 ,   4.677 ,   6.453 ,  23.561 ,
     *                  0.102 ,   3.354 ,  22.954 ,  33.796           /
C
C     data for cobalt
C
      DATA ACO      /  11.9120,  20.4424,  28.3863, 132.5038, 600.7461,
     *                 33.3092, 237.4331, 977.2502,
     *                 55.5396, 318.8169, 619.6366,
     *                 32.6900,  83.8694, 107.4378,
     *                 11.2593,  38.2239,  22.9964, 261.3486, 637.1485,
     *                 23.0233,  41.6599, 264.6460, 181.6699,
     *                 16.0356,   7.8633,  70.3158, 423.3512, 742.3553,
     *                  0.                                            /
      DATA GCO      /   0.112 ,   0.341 ,   0.809 ,   3.808 ,   6.723 ,
     *                  2.057 ,   3.484 ,   7.210 ,
     *                  2.405 ,   5.133 ,   8.097 ,
     *                  2.084 ,   5.291 ,   8.426 ,
     *                  0.135 ,   0.517 ,   1.606 ,   6.772 ,  12.622 ,
     *                  2.512 ,   4.348 ,   8.253 ,  15.377 ,
     *                  0.132 ,   0.863 ,   3.086 ,  11.789 ,  23.263 ,
     *                  0.                                            /
C
C     data for nickel
C
      DATA ANI      /   7.1268,  12.4486,  11.9953,  10.0546, 114.1658,
     *                391.2064,
     *                 26.3908, 213.8081, 938.7927,
     *                  4.1421,  37.3781,  25.9712, 333.3397, 311.1633,
     *                 33.1031, 184.1854, 136.7072,
     *                 11.1915,   5.4174,  53.6793, 460.6781, 380.0056,
     *                  0.                                            /
      DATA GNI      /   0.026 ,   0.137 ,   0.315 ,   1.778 ,   4.029 ,
     *                  6.621 ,
     *                  2.249 ,   4.042 ,   7.621 ,
     *                  0.191 ,   1.235 ,   3.358 ,   8.429 ,  17.096 ,
     *                  3.472 ,   9.065 ,  16.556 ,
     *                  0.194 ,   1.305 ,   5.813 ,  14.172 ,  26.169 ,
     *                  0.                                            /
C
C     data for copper
C
      DATA ACU      /  11.0549, 238.9423,  10.3077, 126.2990,1073.3876,
     *                 30.0000,  50.0000,  60.0000,
     *                 19.2984,  50.5974, 240.2021,1216.9016,
     *                 48.3048, 583.2011, 320.4931,
     *                  4.0155,  70.3264, 313.1213, 536.5331,
     *                  0.                                            /
      DATA GCU      /   4.212 ,   7.227 ,   1.493 ,   5.859 ,   9.709 ,
     *                  7.081 ,   9.362 ,  10.130 ,
     *                  2.865 ,   8.260 ,  14.431 ,  18.292 ,
     *                  9.650 ,  14.640 ,  24.320 ,
     *                  0.337 ,   8.520 ,  16.925 ,  28.342 ,
     *                  0.                                            /
C
C     data for zinc
C
      DATA AZN      /  15.9880, 484.0042,  18.5863, 123.4134,
     *                  3.0000, 189.0000,
     *                  6.1902,  38.9317, 204.8780,
     *                 10.2588,  89.3771, 370.3640,  30.0000, 128.0000,
     *                 24.6904, 106.7491, 439.5586,
     *                  0.                                            /
      DATA GZN      /   4.546 ,   8.840 ,  10.247 ,  16.620 ,
     *                 11.175 ,  16.321 ,
     *                  6.113 ,  12.964 ,  16.444 ,
     *                  7.926 ,  13.633 ,  24.353 ,  16.286 ,  24.910 ,
     *                 10.291 ,  20.689 ,  32.077 ,
     *                  0.                                            /
C
C
      DATA ICOMP /0/
C
c     SAVE ALF,GAM,XL,CHION,INDEX0,INDEXS,INDEXM,IGPR,IG0,ICOMP
C
C     Initialization of auxiliary arrays (executed only once)
C
      IF(ICOMP.NE.0) GO TO 30
      IND=1
      DO 10 K=1,NIONS
         INDEXS(K)=IND
         IND=IND+IS(K)
   10 CONTINUE
      IND=1
      DO 20 K=1,NSS
         INDEXM(K)=IND
         IND=IND+IM(K)
   20 CONTINUE
      ICOMP=1
C
   30 CONTINUE
      MODE=MODPF(IAT)
      IF(MODE.GT.0) GO TO 80

      IF(IAT.EQ.26 .AND. IZI.GE.4 .AND. IZI.LE.9) GO TO 170
      IF(IAT.EQ.28 .AND. IZI.GE.4 .AND. IZI.LE.9) GO TO 171
      IF(IAT.GT.30.AND.IZI.LE.3) GO TO 90
      IF(IZI.LE.0.OR.IAT.LE.0.OR.IAT.GT.30) GO TO 70

      IF(IZI.GT.5) THEN
         IF(IAT.LT.IZI) THEN
            U=1.
            DUT=0.
            DUN=0.
            RETURN
         END IF
         IF(IAT.GT.8) THEN
            U=IGLE(IAT-IZI+1)
            DUT=0.
            DUN=0.
            RETURN
         END IF
         CALL PFCNO(IAT,IZI,T,ANE,U)
         RETURN  
      END IF
c
c     Irwin partition functions by default
c
      if(iirwin.gt.0.and.t.lt.16000.) then
      if(izi.le.2) then
         call mpartf(iat,izi,0,t,u0,du0)
         u=u0
         return
      end if
      end if
c
      IF(MODE.LT.0) GO TO 70
      I0=INDEX0(IZI,IAT)
      IF(I0.LE.0) GO TO 60
C
C     Traving, Baschek, Holweger formula
C
      QZ=IZI
      THET=5.0404D3/T
      A=31.321*QZ*QZ*THET
      XMAX2=XMAX*XMAX
      QQ=XMAX/4.*(XMAX2+XMAX+SIXTH+A+A*A*HALF/XMAX2)
      QAS1=XMAX*THIRD*(XMAX2+TRHA*XMAX+HALF)
      IS0=INDEXS(I0)
      ISS=IS0+IS(I0)-1
      SU1=0.
      SU2=0.
      SQA=0.
      SQQ=0.
      SQT=0.
      SQ2=0.
      DO 50 K=IS0,ISS
         XXL=XL(K)
         GPR=IGPR(K)
         X=CHION(K)*THET
         EX=0.
         IF(X.LT.30) EX=EXP(-X*2.30258029299405)
         SQQ=SQQ+GPR*EX
         QAS=(QAS1-XXL*THIRD*(XXL*XXL+TRHA*XXL+HALF)+(XMAX-XXL)*
     *       (UN+A*HALF/XXL/XMAX)*A)*GPR*EX
         SQA=SQA+QAS
         SQ2=SQ2+QAS*CHION(K)
         SQT=SQT+GPR*(XMAX-XXL)*(UN+A/XMAX/XXL)*EX
         M0=INDEXM(K)
         M1=M0+IM(K)-1
         AL1=0.
         AL2=0.
         DO 40 M=M0,M1
            XG=GAM(M)*THET
            IF(XG.GT.20.) GO TO 40
            XM=EXP(-XG*2.30258029299405)*ALF(M)
            AL1=AL1+XM
            AL2=AL2+GAM(M)*XM
   40    CONTINUE
         SU1=SU1+AL1
         SU2=SU2+AL2
   50 CONTINUE
      U=IG0(I0)
      U=U+SU1+SQA
      IF(U.LT.0.) U=IG0(I0)
      DUT=(2.302580293*THET*(SU2+SQ2)+QQ*SQQ-A*SQT)/T
      DUN=-QQ*SQQ/ANE
      RETURN
C
C     constant value of partition function for some ions (even if
C     MODPF = 0)
C
   60 U=-I0
      DUT=0.
      DUN=0.
      RETURN
C
C     non-standard, user supplied formula
C
   70 CALL PFSPEC(IAT,IZI,U,DUT,DUN)
      RETURN
C
C     Partition functions for Iron  (From Sparks and Fischel)
C
  170 CALL PFFE(IZI,T,ANE,U,DUT,DUN)
      RETURN
C
C     Partition functions for Nickel (from Kurucz predicted levels)
C
  171 CALL PFNI(IZI,T,U,DUT,DUN)
      RETURN
C
C     Opacity Project value
C
   80 call opfrac(iat,izi,t,ane,u,opfra)
      DUT=0.
      DUN=0.
      RETURN
C
C     Modified Kurucz partition functions for IAT > 30
C
   90 CALL PFHEAV(IAT,IZI,3,T,ANE,U)
      DUT=0.
      DUN=0.
      RETURN
      END

C
C ********************************************************************
C

      SUBROUTINE PFCNO(IAT,IZI,T,ANE,PF)
c     ===================================

c     Partition functions for the high ions of CNO:

c       (a) H-like ions:  C VI, N VII, O VIII
c       (b) He-like ions: N VI, O VII
c       (c) O VI from Sparks & Fischel, 1971, NASA SP-3066

c     Output:  PF   partition function

      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'

      PARAMETER (P1=0.1402,P2=0.1285,P3=1.,P4=3.15,P5=4.)
      DIMENSION TT(35),PN(10)
      DIMENSION P6A(24),P6B(10,11) 

      DATA TT /
     * 18.,19.,20.,21.,22.,23.,24.,25.,26.,27.,28.,29.,30.,
     * 32.,34.,36.,38.,40.,42.,44.,46.,48.,
     * 50.,55.,60.,65.,70.,75.,80.,85.,90.,95.,100.,125.,150./

      DATA PN /-2.,-1.,0.,1.,2.,3.,4.,5.,6.,7./
c
      DATA P6A /
     * 0.302, 0.302, 0.302, 0.303, 0.303, 0.304, 0.305, 0.306,
     * 0.307, 0.308, 0.310, 0.312, 0.313, 0.318, 0.322, 0.327,
     * 0.333, 0.339, 0.346, 0.353, 0.360, 0.367, 0.375, 0.394/
c
      DATA P6B /
     * 0.414,9*0.413,
     * 0.436,0.433,0.433,7*0.432,
     * 0.472,0.458,0.453,7*0.451,
     * 0.560,0.499,0.478,0.471,0.469,5*0.468,
     * 0.762,0.593,0.522,0.497,0.489,0.486,4*0.485,
     * 1.090,0.782,0.611,0.539,0.513,0.505,0.502,3*0.501,
     * 1.478,1.070,0.775,0.615,0.550,0.527,0.519,0.517,0.516,0.516,
     * 1.867,1.408,1.018,0.749,0.612,0.557,0.539,0.533,0.531,0.530,
     * 2.233,1.752,1.306,0.944,0.713,0.604,0.564,0.550,0.545,0.544,
     * 3.665,3.166,2.668,2.176,1.700,1.269,0.934,0.735,0.648,0.616,
     * 4.633,4.133,3.633,3.134,2.636,2.146,1.674,1.254,0.942,0.763/     
     
     
      IF (IAT.LT.6 .OR. IAT.GT.8 .OR. IZI.LE.5) THEN
         PF=0.
         PRINT *,' Routine PFCNO does not provide U for this ion '
     *         ,IAT,IZI
         STOP
      END IF
      IF (IZI.GT.IAT) THEN
         PF=1.
         RETURN
      END IF
     
      TK = BOLK*T
      IZIT = IAT-IZI
     
C    1. H-like case

      IF(IZIT.EQ.0) THEN
        ANEL = LOG(ANE)
        AA = 0.09*EXP(ANEL/6.)/SQRT(T)
        ANEL23 = EXP(-2./3.*ANEL)
        Z  = FLOAT(IZI)
        Z2 = Z*Z
        CBZ = 2.*8.59D14*Z*Z*Z
        E0KT = EH*Z2/TK
        
        PF=0.
        DO II=1,NLMX
          XN  = FLOAT(II)
          XN2 = XN*XN
          
          IF(XN.LE.3.01) THEN
            XKN = 1.
          ELSE
            XN1 = 1./(XN+1.)
            XKN = 16./3.*XN*XN1*XN1
          END IF
          BETA = CBZ*XKN/XN2/XN2*ANEL23
          X=EXP(P4*LOG(1.+P3*AA))
          C1=P1*(X+P5*(Z-1.)*AA*AA*AA)
          C2=P2*X
          F=(C1*BETA*BETA*BETA)/(1.+C2*BETA*SQRT(BETA))
          WI=F/(1.+F)
          
          EE = EXP(-E0KT/XN2)
          PF = PF+XN2*WI*EE
        END DO
        PF = 2.*PF
      END IF

C    2. He-like case

      IF(IZIT.EQ.1) THEN
        PF=1.
      END IF

C    3. O VI

      IF(IZIT.EQ.2) THEN
     
        IF(T.LT.18000.) THEN
           PF=2.
        ELSE
          NA = 24
          NB = 11
          PNE = LOG10(ANE*TK)
          T0  = 0.001*T
          J = 1
          IF(PNE.LT.PN(1)) GO TO 15
          IF(PNE.GT.PN(10)) THEN
            J1 = 10
            J2 = 10
            GO TO 16
          END IF
          DO J=1,9
            IF(PNE.GE.PN(J) .AND. PNE.LT.PN(J+1)) GO TO 15
          END DO
   15     J1 = J
          J2 = J+1
          IF(PNE.LT.PN(1)) J2 = 1
   16     DO I=1,34
            IF(T0.GE.TT(I) .AND. T0.LT.TT(I+1)) GO TO 25
          END DO
   25     I1 = I
          I2 = I+1
          IF(T0.GT.TT(35)) THEN
            I1 = 35
            I2 = 35
          END IF
          IF(I2.LE.24) THEN
            PX1=P6A(I1)
            PX2=P6A(I1)
            PY1=P6A(I2)
            PY2=P6A(I2)
          ELSE IF (I1.EQ.24) THEN
            PX1=P6A(I1)
            PX2=P6A(I1)
            PY1=P6B(J1,I2-24)
            PY2=P6B(J2,I2-24)
          ELSE
            PX1=P6B(J1,I1-24)
            PX2=P6B(J2,I1-24)
            PY1=P6B(J1,I2-24)
            PY2=P6B(J2,I2-24)
          END IF
          DLGUNX=PX2-PX1
          PX=PX1+(PNE-PN(J1))*DLGUNX
          DLGUNY=PY2-PY1
          PY=PY1+(PNE-PN(J1))*DLGUNY
          DELT=TT(I2)-TT(I1)
          IF(DELT.NE.0.) THEN
            DLGUT=(PY-PX)/DELT
            PF=PX+(T0-TT(I1))*DLGUT
          ELSE
            PF=PX
          ENDIF
          PF=EXP(2.302585093*PF)
        END IF
      END IF

      RETURN
      END

C
C
C *******************************************************************
C
C
      SUBROUTINE PFSPEC(IAT,IZI,U,DUT,DUN)
C     ====================================
 
C     Non-standard evaluation of the partition function
C     user supplied procedure
C
C     Input:
C      IAT   - atomic number
C      IZI   - ionic charge (=1 for neutrals, =1 for once ionized, etc)
C      T     - temperature
C      ANE   - electron density
C      XMAX  - principal quantum number of the last bound level
C
C     Output:
C      U     - partition function
C      DUT   - derivative dU/dT
C      DUN   - derivative dU/d(ANE)
C
      INCLUDE 'IMPLIC.FOR'
C
      U=0.
      DUT=0.
      DUN=0.
      IF(IAT.EQ.10) THEN
         IF(IZI.EQ.5) U=9.
         IF(IZI.EQ.6) U=6.
         IF(IZI.EQ.7) U=12.
         IF(IZI.EQ.8) U=8.
         IF(IZI.EQ.9) U=1.
         RETURN
      END IF
      IF(IAT.EQ.16) THEN
         IF(IZI.EQ.5) U=1.
         IF(IZI.EQ.6) U=2.
         IF(IZI.EQ.7) U=1.
         IF(IZI.EQ.8) U=6.
         IF(IZI.EQ.9) U=9.
         RETURN
      END IF

      RETURN
        END
C
C
C ********************************************************************
C
C
C
      subroutine pffe(ion,t,ane,pf,dut,dun)
c     =====================================
c
c     partition functions for Fe IV to Fe IX
c     after Fischel and Sparks, 1971, NASA SP-3066
c
c     Output:  PF   partition function
c              DUT  d(PF)/dT
c              DUN  d(PF)/d(ANE)
c
      include 'IMPLIC.FOR'
c
      dimension tt(50),pn(10),nca(9)
      dimension p4a(22),p4b(10,28), 
     *          p5a(30),p5b(10,20),
     *          p6a(37),p6b(10,13),
     *          p7a(40),p7b(10,10),
     *          p8a(41),p8b(10,9),
     *          p9a(45),p9b(10,5)
c
      parameter (xen=2.302585093,xmil=0.001,xmilen=xmil*xen)
      parameter (xbtz=1.38054d-16)
c
      data nca /3*0,22,30,37,40,41,45/
     *     nne /10/
c
      data tt /
     * 3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.,15.,16.,17.,18.,19.,
     * 20.,21.,22.,23.,24.,25.,26.,27.,28.,29.,30.,
     * 32.,34.,36.,38.,40.,42.,44.,46.,48.,
     * 50.,55.,60.,65.,70.,75.,80.,85.,90.,95.,100.,125.,150./
c
      data pn /-2.,-1.,0.,1.,2.,3.,4.,5.,6.,7./
c
      data p4a /
     * 0.778, 0.778, 0.778, 0.779, 0.783, 0.789, 0.801, 0.818,
     * 0.842, 0.871, 0.906, 0.945, 0.987, 1.030, 1.074, 1.117,
     * 1.160, 1.201, 1.242, 1.280, 1.317, 1.353/
c
      data p4b /
     * 1.406,1.393,1.389,7*1.387,
     * 1.464,1.434,1.424,1.421,1.420,5*1.419,
     * 1.546,1.483,1.461,1.454,1.451,1.451,4*1.450,
     * 1.665,1.547,1.503,1.488,1.482,1.481,4*1.480,
     * 1.826,1.636,1.553,1.524,1.514,1.510,4*1.509,
     * 2.024,1.755,1.618,1.564,1.546,1.540,1.538,3*1.537,
     * 2.480,2.087,1.814,1.674,1.619,1.599,1.593,1.591,1.590,1.590,
     * 2.945,2.489,2.105,1.846,1.717,1.667,1.649,1.643,1.641,1.640,
     * 3.379,2.897,2.452,2.089,1.859,1.751,1.710,1.696,1.691,1.689,
     * 3.774,3.283,2.808,2.381,2.054,1.864,1.782,1.751,1.741,1.738,
     * 4.133,3.637,3.150,2.688,2.292,2.015,1.871,1.814,1.793,1.786,
     * 4.460,3.962,3.468,2.989,2.549,2.199,1.984,1.886,1.848,1.835,
     * 4.757,4.258,3.762,3.274,2.809,2.406,2.121,1.972,1.908,1.886,
     * 5.029,4.530,4.032,3.539,3.061,2.624,2.279,2.073,1.976,1.939,
     * 5.279,4.780,4.281,3.785,3.299,2.840,2.450,2.189,2.051,1.996,
     * 5.510,5.010,4.511,4.013,3.522,3.050,2.628,2.318,2.136,2.057,
     * 6.014,5.514,5.014,4.515,4.018,3.530,3.065,2.666,2.381,2.228,
     * 6.435,5.935,5.435,4.936,4.437,3.943,3.460,3.022,2.658,2.422,
     * 6.794,6.294,5.794,5.294,4.794,4.297,3.807,3.343,2.939,2.631,
     * 7.102,6.602,6.102,5.602,5.102,4.604,4.110,3.638,3.194,2.845,
     * 7.370,6.870,6.370,5.870,5.370,4.871,4.375,3.892,3.439,3.052,
     * 7.606,7.106,6.606,6.106,5.605,5.106,4.608,4.125,3.661,3.249,
     * 7.815,7.315,6.814,6.314,5.814,5.314,4.816,4.333,3.851,3.418,
     * 8.001,7.501,7.001,6.500,6.000,5.500,5.001,4.511,4.032,3.586,
     * 8.168,7.668,7.168,6.668,6.168,5.667,5.168,4.680,4.197,3.741,
     * 8.319,7.819,7.319,6.819,6.319,5.818,5.319,4.832,4.347,3.884,
     * 8.900,8.399,7.899,7.399,6.899,6.398,5.898,5.405,4.917,4.431,
     * 9.294,8.794,8.294,7.793,7.293,6.793,6.292,5.799,5.306,4.824/
c
      data p5a /
     * 1.235, 1.276, 1.301, 1.321, 1.339, 1.359, 1.381, 1.405,
     * 1.432, 1.460, 1.489, 1.518, 1.546, 1.574, 1.601, 1.627,
     * 1.652, 1.675, 1.697, 1.718, 1.738, 1.757, 1.775, 1.792,
     * 1.808, 1.823, 1.838, 1.851, 1.877, 1.900/
c
      data p5b /
     * 1.943,1.928,1.923,7*1.921,
     * 2.011,1.964,1.947,1.942,1.941,5*1.940,
     * 2.144,2.025,1.980,1.965,1.960,1.958,4*1.957,
     * 2.361,2.137,2.032,1.993,1.980,1.976,1.975,3*1.974,
     * 2.646,2.315,2.121,2.035,2.004,1.994,1.991,1.990,1.989,1.989,
     * 2.960,2.553,2.260,2.102,2.037,2.015,2.007,2.005,2.004,2.004,
     * 3.274,2.823,2.450,2.205,2.086,2.040,2.025,2.020,2.018,2.018,
     * 3.575,3.101,2.674,2.348,2.158,2.075,2.045,2.036,2.032,2.031,
     * 4.251,3.757,3.275,2.829,2.466,2.234,2.124,2.083,2.069,2.064,
     * 4.822,4.324,3.829,3.346,2.895,2.522,2.278,2.161,2.116,2.100,
     * 5.308,4.808,4.310,3.816,3.334,2.888,2.525,2.297,2.187,2.145,
     * 5.725,5.225,4.726,4.228,3.736,3.260,2.828,2.496,2.294,2.206,
     * 6.088,5.589,5.089,4.590,4.093,3.604,3.139,2.733,2.447,2.291,
     * 6.407,5.907,5.407,4.908,4.409,3.915,3.433,2.988,2.629,2.399,
     * 6.689,6.189,5.689,5.189,4.690,4.193,3.704,3.236,2.832,2.535,
     * 6.940,6.440,5.940,5.440,4.941,4.443,3.949,3.469,3.038,2.687,
     * 7.166,6.666,6.166,5.666,5.166,4.667,4.171,3.684,3.237,2.847,
     * 7.370,6.870,6.369,5.869,5.369,4.870,4.373,3.882,3.417,3.008,
     * 8.150,7.649,7.149,6.649,6.149,5.649,5.149,4.651,4.167,3.700,
     * 8.677,8.177,7.676,7.176,6.676,6.176,5.676,5.176,4.687,4.203/
c
      data p6a /
     * 1.218, 1.273, 1.309, 1.335, 1.358, 1.379, 1.400, 1.421,
     * 1.442, 1.463, 1.484, 1.504, 1.523, 1.542, 1.560, 1.577,
     * 1.594, 1.609, 1.624, 1.638, 1.652, 1.664, 1.677, 1.688,
     * 1.699, 1.709, 1.719, 1.729, 1.746, 1.762, 1.777, 1.790,
     * 1.803, 1.814, 1.825, 1.834, 1.843/
c
      data p6b /
     * 1.862,1.855,1.853,7*1.852,
     * 1.958,1.900,1.880,1.874,1.872,5*1.871,
     * 2.264,2.045,1.944,1.906,1.894,1.890,4*1.888,
     * 2.776,2.386,2.119,1.984,1.930,1.912,1.906,1.904,2*1.903,
     * 3.321,2.856,2.453,2.165,2.012,1.949,1.927,1.920,1.918,1.917,
     * 3.821,3.333,2.868,2.465,2.178,2.025,1.963,1.941,1.934,1.932,
     * 4.266,3.771,3.285,2.825,2.434,2.164,2.027,1.972,1.953,1.947,
     * 4.662,4.164,3.670,3.187,2.739,2.372,2.135,2.022,1.980,1.965,
     * 5.015,4.516,4.019,3.527,3.052,2.624,2.295,2.102,2.019,1.988,
     * 5.332,4.832,4.344,3.838,3.351,2.889,2.493,2.217,2.075,2.017,
     * 5.618,5.118,4.619,4.121,3.628,3.149,2.711,2.364,2.155,2.058,
     * 6.710,6.210,5.710,5.210,4.711,4.213,3.719,3.241,2.807,2.462,
     * 7.446,6.946,6.446,5.946,5.446,4.946,4.447,3.952,3.474,3.022/
c
      data p7a /
     * 1.074,1.130,1.167,1.194,1.215,1.234,1.250,1.266,1.280,1.293,
     * 1.306,1.318,1.329,1.340,1.350,1.360,1.369,1.378,1.386,1.394,
     * 1.401,1.408,1.415,1.421,1.427,1.433,1.439,1.444,1.454,1.463,
     * 1.471,1.479,1.486,1.492,1.498,1.504,1.509,1.514,1.525,1.534/
c
      data p7b /
     * 1.555,1.546,1.544,1.543,6*1.542,
     * 1.617,1.572,1.557,1.552,1.550,1.550,4*1.549,
     * 1.798,1.648,1.587,1.566,1.559,1.557,4*1.556,
     * 2.134,1.832,1.666,1.597,1.573,1.565,1.563,1.562,2*1.561,
     * 2.550,2.138,1.836,1.671,1.602,1.578,1.570,1.568,2*1.567,
     * 2.968,2.504,2.102,1.816,1.665,1.603,1.582,1.575,2*1.572,
     * 3.359,2.875,2.419,2.037,1.779,1.651,1.601,1.584,1.579,1.577,
     * 3.718,3.224,2.745,2.305,1.953,1.736,1.636,1.599,1.586,1.582,
     * 5.097,4.598,4.098,3.601,3.110,2.638,2.217,1.899,1.719,1.643,
     * 6.026,5.526,5.026,4.527,4.028,3.531,3.042,2.576,2.170,1.885/
c
      data p8a /
     * 0.809,0.849,0.875,0.894,0.908,0.918,0.927,0.934,0.939,0.944,
     * 0.948,0.952,0.955,0.958,0.960,0.962,0.964,0.966,0.967,0.969,
     * 0.970,0.971,0.973,0.974,0.975,0.975,0.976,0.977,0.978,0.980,
     * 0.981,0.982,0.983,0.984,0.984,0.985,0.986,0.986,0.987,0.988,
     * 0.989/
c
      data p8b /
     * 0.992,0.991,8*0.990,
     * 1.000,0.994,0.992,7*0.991,
     * 1.032,1.005,0.996,0.993,0.992,5*0.991,
     * 1.129,1.040,1.008,0.997,0.993,5*0.992,
     * 1.335,1.132,1.042,1.009,0.998,0.994,0.993,0.993,2*0.992,
     * 1.640,1.312,1.121,1.038,1.007,0.998,0.994,3*0.993,
     * 1.987,1.573,1.269,1.101,1.030,1.005,0.997,2*0.994,0.993,
     * 3.514,3.017,2.526,2.053,1.628,1.305,1.119,1.039,1.010,1.000,
     * 4.569,4.069,3.569,3.072,2.580,2.103,1.671,1.336,1.136,1.048/
c
      data p9a /39*0.000,0.001,0.002,0.005,0.008,0.014,0.021/
c
      data p9b /
     * 2*0.032,8*0.031,
     * 0.048,0.045,8*0.044,
     * 0.076,0.065,0.061,0.060,6*0.059,
     * 1.128,0.722,0.429,0.271,0.207,0.184,0.177,0.174,2*0.173,
     * 2.696,2.200,1.712,1.249,0.848,0.564,0.415,0.354,0.333,0.327/
c
      na=nca(ion)
      nb=50-na
      pne=log10(ane*xbtz*t)
      t0=xmil*t
      j=1
      if(pne.lt.pn(1)) go to 15
      if(pne.gt.pn(nne)) then
        j1=nne
        j2=nne
        goto 16
      endif
      do 10 j=1,nne-1      
        if(pne.ge.pn(j).and.pne.lt.pn(j+1)) go to 15
   10 continue
   15 j1=j
      j2=j1+1
      if(pne.lt.pn(1)) j2=1
   16 do 20 i=1,49
         if(t0.ge.tt(i).and.t0.lt.tt(i+1)) go to 25
   20 continue
   25 i1=i
      i2=i+1
      if(t0.gt.tt(50)) then
        i1=50
        i2=50
      endif
      if(i2.le.na) then
         if(ion.eq.4) then
           px1=p4a(i1)
           px2=p4a(i1)
           py1=p4a(i2)
           py2=p4a(i2)
         else if(ion.eq.5) then
           px1=p5a(i1)
           px2=p5a(i1)
           py1=p5a(i2)
           py2=p5a(i2)
         else if(ion.eq.6) then
           px1=p6a(i1)
           px2=p6a(i1)
           py1=p6a(i2)
           py2=p6a(i2)
         else if(ion.eq.7) then
           px1=p7a(i1)
           px2=p7a(i1)
           py1=p7a(i2)
           py2=p7a(i2)
         else if(ion.eq.8) then
           px1=p8a(i1)
           px2=p8a(i1)
           py1=p8a(i2)
           py2=p8a(i2)
         else if(ion.eq.9) then
           px1=p9a(i1)
           px2=p9a(i1)
           py1=p9a(i2)
           py2=p9a(i2)
         endif
       else if(i1.eq.na) then
         if(ion.eq.4) then
           px1=p4a(i1)
           px2=p4a(i1)
           py1=p4b(j1,i2-na)
           py2=p4b(j2,i2-na)
         else if(ion.eq.5) then
           px1=p5a(i1)
           px2=p5a(i1)
           py1=p5b(j1,i2-na)
           py2=p5b(j2,i2-na)
         else if(ion.eq.6) then
           px1=p6a(i1)
           px2=p6a(i1)
           py1=p6b(j1,i2-na)
           py2=p6b(j2,i2-na)
         else if(ion.eq.7) then
           px1=p7a(i1)
           px2=p7a(i1)
           py1=p7b(j1,i2-na)
           py2=p7b(j2,i2-na)
         else if(ion.eq.8) then
           px1=p8a(i1)
           px2=p8a(i1)
           py1=p8b(j1,i2-na)
           py2=p8b(j2,i2-na)
         else if(ion.eq.9) then
           px1=p9a(i1)
           px2=p9a(i1)
           py1=p9b(j1,i2-na)
           py2=p9b(j2,i2-na)
         endif
      else
         if(ion.eq.4) then
           px1=p4b(j1,i1-na)
           px2=p4b(j2,i1-na)
           py1=p4b(j1,i2-na)
           py2=p4b(j2,i2-na)
         else if(ion.eq.5) then
           px1=p5b(j1,i1-na)
           px2=p5b(j2,i1-na)
           py1=p5b(j1,i2-na)
           py2=p5b(j2,i2-na)
         else if(ion.eq.6) then
           px1=p6b(j1,i1-na)
           px2=p6b(j2,i1-na)
           py1=p6b(j1,i2-na)
           py2=p6b(j2,i2-na)
         else if(ion.eq.7) then
           px1=p7b(j1,i1-na)
           px2=p7b(j2,i1-na)
           py1=p7b(j1,i2-na)
           py2=p7b(j2,i2-na)
         else if(ion.eq.8) then
           px1=p8b(j1,i1-na)
           px2=p8b(j2,i1-na)
           py1=p8b(j1,i2-na)
           py2=p8b(j2,i2-na)
         else if(ion.eq.9) then
           px1=p9b(j1,i1-na)
           px2=p9b(j2,i1-na)
           py1=p9b(j1,i2-na)
           py2=p9b(j2,i2-na)
         endif
      end if
      dlgunx=px2-px1
      px=px1+(pne-pn(j1))*dlgunx
      dlguny=py2-py1
      py=py1+(pne-pn(j1))*dlguny
      delt=tt(i2)-tt(i1)
      if(delt.ne.0.) then
        dlgut=(py-px)/delt
        pf=px+(t0-tt(i1))*dlgut
        dlgun=dlgunx+(t0-tt(i1))/delt*(dlguny-dlgunx)
      else
        dlgut=0.
        pf=px
        dlgun=dlgunx
      endif
      pf=exp(xen*pf)
      dut=xmilen*pf*dlgut
      dun=(dlgun*pf-t*dut)/ane
      return
      end
C
C
C ********************************************************************
C
C
      subroutine pfni(ion,t,pf,dut,dun)
c     =================================
c
c     partition functions for Ni IV to Ni IX
c
c     this routine interpolates within a grid
c     calculated from all levels predicted by
c     Kurucz (1992), i.e. over 12,000 levels per ion.
c     the partition functions depend only on T !
c     (i.e. no level dissolution with increasing density)
c     TL  27-DEC-1994, 23-JAN-1995
c
c     Output:  PF   partition function
c              DUT  d(PF)/dT
c              DUN  d(PF)/d(ANE)
c
      include 'IMPLIC.FOR'
c
      dimension g0(6)
      dimension p4a(190),p4b(170)
      dimension p5a(190),p5b(170)
      dimension p6a(190),p6b(170)
      dimension p7a(190),p7b(170)
      dimension p8a(190),p8b(170)
      dimension p9a(190),p9b(170)
      parameter (xen=2.302585093,xmil=0.001)
c
      data g0/28.,25.,6.,25.,28.,21./
c
      data p4a/
     .    1.447,1.464,1.482,1.501,1.518,1.535,1.551,1.567,1.582,1.596,
     .    1.610,1.623,1.636,1.648,1.659,1.671,1.681,1.692,1.702,1.711,
     .    1.721,1.730,1.739,1.748,1.757,1.765,1.774,1.782,1.791,1.799,
     .    1.808,1.816,1.824,1.833,1.841,1.850,1.859,1.868,1.877,1.886,
     .    1.895,1.905,1.914,1.924,1.934,1.945,1.955,1.966,1.977,1.989,
     .    2.000,2.012,2.025,2.037,2.050,2.063,2.077,2.091,2.105,2.119,
     .    2.134,2.149,2.164,2.179,2.195,2.211,2.227,2.243,2.260,2.276,
     .    2.293,2.310,2.327,2.344,2.362,2.379,2.397,2.414,2.432,2.449,
     .    2.467,2.484,2.502,2.519,2.537,2.554,2.571,2.588,2.606,2.623,
     .    2.640,2.657,2.674,2.690,2.707,2.723,2.740,2.756,2.772,2.788,
     .    2.804,2.819,2.835,2.850,2.866,2.881,2.896,2.911,2.925,2.940,
     .    2.954,2.969,2.983,2.997,3.010,3.024,3.038,3.051,3.064,3.077,
     .    3.090,3.103,3.116,3.128,3.141,3.153,3.165,3.177,3.189,3.201,
     .    3.213,3.224,3.235,3.247,3.258,3.269,3.280,3.291,3.301,3.312,
     .    3.322,3.332,3.343,3.353,3.363,3.373,3.382,3.392,3.402,3.411,
     .    3.421,3.430,3.439,3.448,3.457,3.466,3.475,3.484,3.492,3.501,
     .    3.509,3.518,3.526,3.534,3.542,3.550,3.558,3.566,3.574,3.582,
     .    3.589,3.597,3.604,3.612,3.619,3.626,3.634,3.641,3.648,3.655,
     .    3.662,3.669,3.676,3.682,3.689,3.696,3.702,3.709,3.715,3.722/
      data p4b/
     .    3.589,3.597,3.604,3.612,3.619,3.626,3.634,3.641,3.648,3.655,
     .    3.662,3.669,3.676,3.682,3.689,3.696,3.702,3.709,3.715,3.722,
     .    3.728,3.734,3.740,3.747,3.753,3.759,3.765,3.771,3.777,3.782,
     .    3.788,3.794,3.800,3.805,3.811,3.816,3.822,3.827,3.833,3.838,
     .    3.843,3.849,3.854,3.859,3.864,3.869,3.874,3.879,3.884,3.889,
     .    3.894,3.899,3.904,3.909,3.913,3.918,3.923,3.927,3.932,3.936,
     .    3.941,3.945,3.950,3.954,3.959,3.963,3.967,3.972,3.976,3.980,
     .    3.984,3.988,3.993,3.997,4.001,4.005,4.009,4.013,4.017,4.021,
     .    4.024,4.028,4.032,4.036,4.040,4.043,4.047,4.051,4.055,4.058,
     .    4.062,4.065,4.069,4.072,4.076,4.079,4.083,4.086,4.090,4.093,
     .    4.097,4.100,4.103,4.107,4.110,4.113,4.116,4.120,4.123,4.126,
     .    4.129,4.132,4.135,4.138,4.141,4.144,4.148,4.151,4.154,4.157,
     .    4.159,4.162,4.165,4.168,4.171,4.174,4.177,4.180,4.182,4.185,
     .    4.188,4.191,4.193,4.196,4.199,4.202,4.204,4.207,4.210,4.212,
     .    4.215,4.217,4.220,4.223,4.225,4.228,4.230,4.233,4.235,4.238,
     .    4.240,4.243,4.245,4.247,4.250,4.252,4.255,4.257,4.259,4.262,
     .    4.264,4.266,4.268,4.271,4.273,4.275,4.278,4.280,4.282,4.284/
      data p5a/
     .    1.398,1.408,1.427,1.446,1.466,1.486,1.506,1.526,1.545,1.564,
     .    1.583,1.601,1.619,1.636,1.652,1.668,1.683,1.698,1.712,1.725,
     .    1.738,1.751,1.763,1.775,1.786,1.797,1.808,1.818,1.828,1.837,
     .    1.846,1.855,1.864,1.873,1.881,1.889,1.897,1.904,1.912,1.919,
     .    1.926,1.933,1.940,1.946,1.953,1.960,1.966,1.972,1.979,1.985,
     .    1.991,1.997,2.003,2.009,2.016,2.022,2.028,2.034,2.040,2.046,
     .    2.052,2.058,2.065,2.071,2.077,2.084,2.090,2.097,2.103,2.110,
     .    2.117,2.124,2.131,2.138,2.145,2.152,2.160,2.167,2.175,2.183,
     .    2.191,2.199,2.207,2.216,2.224,2.233,2.241,2.250,2.259,2.268,
     .    2.278,2.287,2.297,2.306,2.316,2.326,2.336,2.346,2.356,2.367,
     .    2.377,2.387,2.398,2.409,2.419,2.430,2.441,2.452,2.463,2.474,
     .    2.485,2.497,2.508,2.519,2.530,2.542,2.553,2.564,2.576,2.587,
     .    2.599,2.610,2.621,2.633,2.644,2.655,2.667,2.678,2.689,2.701,
     .    2.712,2.723,2.734,2.745,2.757,2.768,2.779,2.790,2.801,2.812,
     .    2.822,2.833,2.844,2.855,2.865,2.876,2.886,2.897,2.907,2.918,
     .    2.928,2.938,2.948,2.958,2.968,2.978,2.988,2.998,3.008,3.018,
     .    3.027,3.037,3.046,3.056,3.065,3.075,3.084,3.093,3.102,3.111,
     .    3.120,3.129,3.138,3.147,3.156,3.164,3.173,3.182,3.190,3.198,
     .    3.207,3.215,3.223,3.232,3.240,3.248,3.256,3.264,3.272,3.279/
      data p5b/
     .    3.120,3.129,3.138,3.147,3.156,3.164,3.173,3.182,3.190,3.198,
     .    3.207,3.215,3.223,3.232,3.240,3.248,3.256,3.264,3.272,3.279,
     .    3.287,3.295,3.303,3.310,3.318,3.325,3.333,3.340,3.347,3.355,
     .    3.362,3.369,3.376,3.383,3.390,3.397,3.404,3.411,3.417,3.424,
     .    3.431,3.438,3.444,3.451,3.457,3.464,3.470,3.476,3.483,3.489,
     .    3.495,3.501,3.507,3.514,3.520,3.526,3.531,3.537,3.543,3.549,
     .    3.555,3.561,3.566,3.572,3.578,3.583,3.589,3.594,3.600,3.605,
     .    3.610,3.616,3.621,3.626,3.632,3.637,3.642,3.647,3.652,3.657,
     .    3.662,3.667,3.672,3.677,3.682,3.687,3.692,3.697,3.701,3.706,
     .    3.711,3.716,3.720,3.725,3.729,3.734,3.738,3.743,3.747,3.752,
     .    3.756,3.761,3.765,3.769,3.774,3.778,3.782,3.786,3.790,3.795,
     .    3.799,3.803,3.807,3.811,3.815,3.819,3.823,3.827,3.831,3.835,
     .    3.839,3.843,3.846,3.850,3.854,3.858,3.862,3.865,3.869,3.873,
     .    3.876,3.880,3.884,3.887,3.891,3.894,3.898,3.901,3.905,3.908,
     .    3.912,3.915,3.918,3.922,3.925,3.929,3.932,3.935,3.939,3.942,
     .    3.945,3.948,3.951,3.955,3.958,3.961,3.964,3.967,3.970,3.974,
     .    3.977,3.980,3.983,3.986,3.989,3.992,3.995,3.998,4.001,4.004/
      data p6a/
     .    0.778,0.804,0.817,0.834,0.854,0.876,0.901,0.928,0.957,0.987,
     .    1.017,1.048,1.079,1.109,1.139,1.169,1.197,1.225,1.253,1.279,
     .    1.304,1.329,1.353,1.376,1.398,1.419,1.440,1.459,1.478,1.497,
     .    1.515,1.532,1.548,1.564,1.580,1.594,1.609,1.623,1.636,1.649,
     .    1.662,1.674,1.686,1.698,1.709,1.720,1.730,1.740,1.750,1.760,
     .    1.769,1.779,1.788,1.796,1.805,1.813,1.821,1.829,1.837,1.845,
     .    1.852,1.860,1.867,1.874,1.881,1.888,1.894,1.901,1.907,1.914,
     .    1.920,1.926,1.932,1.938,1.944,1.950,1.956,1.962,1.968,1.974,
     .    1.979,1.985,1.991,1.996,2.002,2.007,2.013,2.018,2.024,2.029,
     .    2.035,2.041,2.046,2.052,2.057,2.063,2.068,2.074,2.080,2.086,
     .    2.091,2.097,2.103,2.109,2.115,2.121,2.127,2.133,2.139,2.145,
     .    2.152,2.158,2.164,2.171,2.177,2.184,2.190,2.197,2.204,2.211,
     .    2.218,2.225,2.232,2.239,2.246,2.253,2.261,2.268,2.276,2.283,
     .    2.291,2.298,2.306,2.314,2.322,2.330,2.338,2.346,2.354,2.362,
     .    2.370,2.379,2.387,2.395,2.404,2.412,2.420,2.429,2.438,2.446,
     .    2.455,2.463,2.472,2.481,2.489,2.498,2.507,2.516,2.524,2.533,
     .    2.542,2.551,2.560,2.569,2.577,2.586,2.595,2.604,2.613,2.622,
     .    2.631,2.639,2.648,2.657,2.666,2.675,2.683,2.692,2.701,2.710,
     .    2.718,2.727,2.736,2.744,2.753,2.761,2.770,2.779,2.787,2.796/
      data p6b/
     .    2.631,2.639,2.648,2.657,2.666,2.675,2.683,2.692,2.701,2.710,
     .    2.718,2.727,2.736,2.744,2.753,2.761,2.770,2.779,2.787,2.796,
     .    2.804,2.812,2.821,2.829,2.838,2.846,2.854,2.862,2.871,2.879,
     .    2.887,2.895,2.903,2.911,2.919,2.927,2.935,2.943,2.951,2.958,
     .    2.966,2.974,2.982,2.989,2.997,3.005,3.012,3.020,3.027,3.035,
     .    3.042,3.049,3.057,3.064,3.071,3.078,3.086,3.093,3.100,3.107,
     .    3.114,3.121,3.128,3.135,3.141,3.148,3.155,3.162,3.169,3.175,
     .    3.182,3.188,3.195,3.202,3.208,3.214,3.221,3.227,3.234,3.240,
     .    3.246,3.252,3.259,3.265,3.271,3.277,3.283,3.289,3.295,3.301,
     .    3.307,3.313,3.319,3.325,3.330,3.336,3.342,3.348,3.353,3.359,
     .    3.364,3.370,3.376,3.381,3.386,3.392,3.397,3.403,3.408,3.413,
     .    3.419,3.424,3.429,3.434,3.440,3.445,3.450,3.455,3.460,3.465,
     .    3.470,3.475,3.480,3.485,3.490,3.495,3.499,3.504,3.509,3.514,
     .    3.518,3.523,3.528,3.533,3.537,3.542,3.546,3.551,3.555,3.560,
     .    3.564,3.569,3.573,3.578,3.582,3.586,3.591,3.595,3.599,3.604,
     .    3.608,3.612,3.616,3.621,3.625,3.629,3.633,3.637,3.641,3.645,
     .    3.649,3.653,3.657,3.661,3.665,3.669,3.673,3.677,3.681,3.685/
      data p7a/
     .    1.398,1.398,1.398,1.398,1.406,1.425,1.443,1.461,1.480,1.498,
     .    1.516,1.534,1.551,1.568,1.585,1.601,1.616,1.631,1.646,1.660,
     .    1.674,1.687,1.700,1.712,1.724,1.736,1.747,1.758,1.768,1.778,
     .    1.788,1.797,1.806,1.815,1.824,1.832,1.840,1.848,1.855,1.863,
     .    1.870,1.877,1.883,1.890,1.896,1.902,1.908,1.914,1.920,1.925,
     .    1.931,1.936,1.941,1.946,1.951,1.956,1.960,1.965,1.969,1.974,
     .    1.978,1.982,1.986,1.990,1.994,1.998,2.001,2.005,2.009,2.012,
     .    2.016,2.019,2.022,2.026,2.029,2.032,2.035,2.038,2.041,2.044,
     .    2.047,2.050,2.053,2.056,2.059,2.061,2.064,2.067,2.069,2.072,
     .    2.075,2.077,2.080,2.082,2.085,2.088,2.090,2.093,2.095,2.098,
     .    2.100,2.103,2.105,2.107,2.110,2.112,2.115,2.117,2.120,2.122,
     .    2.125,2.127,2.130,2.132,2.135,2.137,2.140,2.142,2.145,2.148,
     .    2.150,2.153,2.155,2.158,2.161,2.163,2.166,2.169,2.172,2.175,
     .    2.178,2.180,2.183,2.186,2.189,2.192,2.195,2.198,2.202,2.205,
     .    2.208,2.211,2.215,2.218,2.221,2.225,2.228,2.232,2.235,2.239,
     .    2.243,2.246,2.250,2.254,2.258,2.261,2.265,2.269,2.273,2.277,
     .    2.282,2.286,2.290,2.294,2.299,2.303,2.307,2.312,2.316,2.321,
     .    2.325,2.330,2.335,2.339,2.344,2.349,2.354,2.359,2.364,2.369,
     .    2.374,2.379,2.384,2.389,2.394,2.399,2.405,2.410,2.415,2.420/
      data p7b/
     .    2.325,2.330,2.335,2.339,2.344,2.349,2.354,2.359,2.364,2.369,
     .    2.374,2.379,2.384,2.389,2.394,2.399,2.405,2.410,2.415,2.420,
     .    2.426,2.431,2.437,2.442,2.448,2.453,2.459,2.464,2.470,2.476,
     .    2.481,2.487,2.493,2.498,2.504,2.510,2.516,2.521,2.527,2.533,
     .    2.539,2.545,2.551,2.556,2.562,2.568,2.574,2.580,2.586,2.592,
     .    2.598,2.604,2.610,2.616,2.622,2.628,2.634,2.640,2.646,2.652,
     .    2.658,2.664,2.670,2.676,2.682,2.687,2.693,2.699,2.705,2.711,
     .    2.717,2.723,2.729,2.735,2.741,2.747,2.753,2.759,2.764,2.770,
     .    2.776,2.782,2.788,2.794,2.799,2.805,2.811,2.817,2.823,2.828,
     .    2.834,2.840,2.846,2.851,2.857,2.863,2.868,2.874,2.879,2.885,
     .    2.891,2.896,2.902,2.907,2.913,2.918,2.924,2.929,2.935,2.940,
     .    2.945,2.951,2.956,2.962,2.967,2.972,2.978,2.983,2.988,2.993,
     .    2.999,3.004,3.009,3.014,3.019,3.025,3.030,3.035,3.040,3.045,
     .    3.050,3.055,3.060,3.065,3.070,3.075,3.080,3.085,3.090,3.095,
     .    3.099,3.104,3.109,3.114,3.119,3.123,3.128,3.133,3.138,3.142,
     .    3.147,3.152,3.156,3.161,3.165,3.170,3.175,3.179,3.184,3.188,
     .    3.193,3.197,3.202,3.206,3.210,3.215,3.219,3.224,3.228,3.232/
      data p8a/
     .    1.447,1.447,1.447,1.447,1.447,1.447,1.459,1.475,1.489,1.504,
     .    1.518,1.531,1.544,1.556,1.568,1.580,1.591,1.602,1.612,1.622,
     .    1.631,1.640,1.649,1.658,1.666,1.674,1.682,1.689,1.696,1.703,
     .    1.710,1.716,1.722,1.728,1.734,1.740,1.745,1.751,1.756,1.761,
     .    1.766,1.770,1.775,1.779,1.784,1.788,1.792,1.796,1.800,1.804,
     .    1.807,1.811,1.814,1.818,1.821,1.824,1.827,1.831,1.834,1.836,
     .    1.839,1.842,1.845,1.848,1.850,1.853,1.855,1.858,1.860,1.863,
     .    1.865,1.867,1.870,1.872,1.874,1.876,1.878,1.880,1.882,1.884,
     .    1.886,1.888,1.890,1.892,1.894,1.896,1.898,1.900,1.902,1.903,
     .    1.905,1.907,1.909,1.911,1.912,1.914,1.916,1.917,1.919,1.921,
     .    1.923,1.924,1.926,1.928,1.929,1.931,1.933,1.934,1.936,1.938,
     .    1.939,1.941,1.943,1.945,1.946,1.948,1.950,1.951,1.953,1.955,
     .    1.957,1.959,1.960,1.962,1.964,1.966,1.968,1.970,1.971,1.973,
     .    1.975,1.977,1.979,1.981,1.983,1.985,1.987,1.989,1.991,1.993,
     .    1.995,1.998,2.000,2.002,2.004,2.006,2.009,2.011,2.013,2.015,
     .    2.018,2.020,2.023,2.025,2.027,2.030,2.032,2.035,2.037,2.040,
     .    2.043,2.045,2.048,2.051,2.053,2.056,2.059,2.062,2.064,2.067,
     .    2.070,2.073,2.076,2.079,2.082,2.085,2.088,2.091,2.094,2.097,
     .    2.100,2.103,2.107,2.110,2.113,2.116,2.120,2.123,2.126,2.130/
      data p8b/
     .    2.070,2.073,2.076,2.079,2.082,2.085,2.088,2.091,2.094,2.097,
     .    2.100,2.103,2.107,2.110,2.113,2.116,2.120,2.123,2.126,2.130,
     .    2.133,2.137,2.140,2.143,2.147,2.151,2.154,2.158,2.161,2.165,
     .    2.168,2.172,2.176,2.180,2.183,2.187,2.191,2.195,2.198,2.202,
     .    2.206,2.210,2.214,2.218,2.222,2.226,2.230,2.233,2.237,2.241,
     .    2.245,2.250,2.254,2.258,2.262,2.266,2.270,2.274,2.278,2.282,
     .    2.286,2.291,2.295,2.299,2.303,2.307,2.312,2.316,2.320,2.324,
     .    2.329,2.333,2.337,2.341,2.346,2.350,2.354,2.359,2.363,2.367,
     .    2.371,2.376,2.380,2.384,2.389,2.393,2.397,2.402,2.406,2.410,
     .    2.415,2.419,2.423,2.428,2.432,2.436,2.441,2.445,2.449,2.454,
     .    2.458,2.462,2.467,2.471,2.475,2.480,2.484,2.488,2.493,2.497,
     .    2.501,2.506,2.510,2.514,2.519,2.523,2.527,2.531,2.536,2.540,
     .    2.544,2.548,2.553,2.557,2.561,2.565,2.570,2.574,2.578,2.582,
     .    2.586,2.591,2.595,2.599,2.603,2.607,2.611,2.616,2.620,2.624,
     .    2.628,2.632,2.636,2.640,2.644,2.648,2.652,2.656,2.661,2.665,
     .    2.669,2.673,2.677,2.681,2.685,2.689,2.693,2.696,2.700,2.704,
     .    2.708,2.712,2.716,2.720,2.724,2.728,2.732,2.736,2.739,2.743/
      data p9a/
     .    1.322,1.322,1.322,1.322,1.322,1.322,1.322,1.322,1.322,1.325,
     .    1.334,1.342,1.351,1.358,1.366,1.373,1.380,1.386,1.392,1.398,
     .    1.404,1.409,1.415,1.420,1.425,1.429,1.434,1.438,1.442,1.446,
     .    1.450,1.454,1.457,1.461,1.464,1.467,1.470,1.473,1.476,1.479,
     .    1.482,1.485,1.487,1.490,1.492,1.495,1.497,1.499,1.501,1.503,
     .    1.505,1.507,1.509,1.511,1.513,1.515,1.517,1.519,1.520,1.522,
     .    1.524,1.525,1.527,1.528,1.530,1.531,1.533,1.534,1.535,1.537,
     .    1.538,1.539,1.541,1.542,1.543,1.545,1.546,1.547,1.548,1.549,
     .    1.551,1.552,1.553,1.554,1.555,1.556,1.558,1.559,1.560,1.561,
     .    1.562,1.563,1.565,1.566,1.567,1.568,1.569,1.570,1.571,1.573,
     .    1.574,1.575,1.576,1.577,1.579,1.580,1.581,1.582,1.584,1.585,
     .    1.586,1.588,1.589,1.590,1.592,1.593,1.594,1.596,1.597,1.599,
     .    1.600,1.602,1.603,1.605,1.606,1.608,1.609,1.611,1.612,1.614,
     .    1.616,1.617,1.619,1.621,1.622,1.624,1.626,1.628,1.630,1.631,
     .    1.633,1.635,1.637,1.639,1.641,1.643,1.645,1.647,1.649,1.651,
     .    1.653,1.655,1.657,1.659,1.661,1.664,1.666,1.668,1.670,1.673,
     .    1.675,1.677,1.679,1.682,1.684,1.686,1.689,1.691,1.694,1.696,
     .    1.699,1.701,1.704,1.706,1.709,1.711,1.714,1.716,1.719,1.722,
     .    1.724,1.727,1.729,1.732,1.735,1.738,1.740,1.743,1.746,1.749/
      data p9b/
     .    1.699,1.701,1.704,1.706,1.709,1.711,1.714,1.716,1.719,1.722,
     .    1.724,1.727,1.729,1.732,1.735,1.738,1.740,1.743,1.746,1.749,
     .    1.751,1.754,1.757,1.760,1.763,1.765,1.768,1.771,1.774,1.777,
     .    1.780,1.783,1.786,1.789,1.792,1.795,1.798,1.801,1.804,1.807,
     .    1.810,1.813,1.816,1.819,1.822,1.825,1.828,1.831,1.834,1.837,
     .    1.840,1.843,1.847,1.850,1.853,1.856,1.859,1.862,1.865,1.869,
     .    1.872,1.875,1.878,1.881,1.884,1.888,1.891,1.894,1.897,1.901,
     .    1.904,1.907,1.910,1.913,1.917,1.920,1.923,1.926,1.930,1.933,
     .    1.936,1.939,1.943,1.946,1.949,1.952,1.956,1.959,1.962,1.965,
     .    1.969,1.972,1.975,1.978,1.982,1.985,1.988,1.992,1.995,1.998,
     .    2.001,2.005,2.008,2.011,2.014,2.018,2.021,2.024,2.027,2.031,
     .    2.034,2.037,2.040,2.044,2.047,2.050,2.053,2.057,2.060,2.063,
     .    2.066,2.070,2.073,2.076,2.079,2.083,2.086,2.089,2.092,2.095,
     .    2.099,2.102,2.105,2.108,2.111,2.115,2.118,2.121,2.124,2.127,
     .    2.131,2.134,2.137,2.140,2.143,2.146,2.149,2.153,2.156,2.159,
     .    2.162,2.165,2.168,2.171,2.175,2.178,2.181,2.184,2.187,2.190,
     .    2.193,2.196,2.199,2.202,2.205,2.208,2.212,2.215,2.218,2.221/
c
      if(t.lt.12000.) then
        pf=g0(ion-3)
        dut=0.
        dun=0.
        return
      endif
c
      it=int(t/1000.)
      if(it.ge.350) it=349
      t1=1000.*it
      t2=t1+1000.
      if(ion.eq.4) then
        if(t.le.200000.) then
          xu1=p4a(it-10)
          xu2=p4a(it-9)
        else
          xu1=p4b(it-180)
          xu2=p4b(it-179)
        endif
      else if(ion.eq.5) then
        if(t.le.200000.) then
          xu1=p5a(it-10)
          xu2=p5a(it-9)
        else
          xu1=p5b(it-180)
          xu2=p5b(it-179)
        endif
      else if(ion.eq.6) then
        if(t.le.200000.) then
          xu1=p6a(it-10)
          xu2=p6a(it-9)
        else
          xu1=p6b(it-180)
          xu2=p6b(it-179)
        endif
      else if(ion.eq.7) then
        if(t.le.200000.) then
          xu1=p7a(it-10)
          xu2=p7a(it-9)
        else
          xu1=p7b(it-180)
          xu2=p7b(it-179)
        endif
      else if(ion.eq.8) then
        if(t.le.200000.) then
          xu1=p8a(it-10)
          xu2=p8a(it-9)
        else
          xu1=p8b(it-180)
          xu2=p8b(it-179)
        endif
      else if(ion.eq.9) then
        if(t.le.200000.) then
          xu1=p9a(it-10)
          xu2=p9a(it-9)
        else
          xu1=p9b(it-180)
          xu2=p9b(it-179)
        endif
      endif
c
      dxt=xmil*(xu2-xu1)
      xu=xu1+(t-t1)*dxt
      pf=exp(xen*xu)
      dut=xen*pf*dxt
      dun=0.
      return
      end

c 
c 
c ******************************************************************
c 
c
      subroutine opfrac(iat,ion,t,ane,pf,fra)
c     =======================================
c
      include 'IMPLIC.FOR'
      parameter (mtemp =100,
     *           melec = 60,
     *           mion1 = 30,
     *           mion2 = 32,
     *           mdat  = 17,
     *           mstag = 258)
      parameter (inp=71)
      dimension frac0(mion2),ioo(mion2),idat(mion1)
      dimension gg(mion1,mdat),g0(mion1),z0(mion2)
      dimension uu(mion1,mdat),u0(mion1)
      dimension indxat(mion1,mdat),indxa(mion1)
      dimension u6(6),u7(7),u8(8),u10(10),u11(11)
      dimension u12(12),u13(13),u14(14),u16(16),u18(18),u20(20)
      dimension u24(24),u25(25),u26(26),u28(28)
      common/pfoptb/pfop(mtemp,melec,mstag),pfophm(mtemp,melec),
     *             frac(mtemp,melec,mstag),
     *             frop(mtemp,melec,mstag),itemp(mtemp)
      equivalence (u6(1),uu(1,3)),(u7(1),uu(1,4)),(u8(1),uu(1,5))
      equivalence (u10(1),uu(1,6)),(u11(1),uu(1,7)),(u12(1),uu(1,8))
      equivalence (u13(1),uu(1,9)),(u14(1),uu(1,10)),(u16(1),uu(1,11))
      equivalence (u18(1),uu(1,12)),(u20(1),uu(1,13)),(u24(1),uu(1,14))
      equivalence (u25(1),uu(1,15)),(u26(1),uu(1,16)),(u28(1),uu(1,17))
      data idat   / 1, 2, 0, 0, 0, 3, 4, 5, 0, 6,
     *              7, 8, 9,10, 0,11, 0,12, 0,13,
     *              0, 0, 0,14,15,16, 0,17, 0, 0/ 
      data indxat /1,2,28*0,
     *            3,4,5,27*0,
     *            6,7,8,9,10,11,12,23*0,
     *            13,14,15,16,17,18,19,20,22*0,
     *            21,22,23,24,25,26,27,28,29,21*0,
     *            30,31,32,33,34,35,36,37,38,39,40,19*0,
     *            41,42,43,44,45,46,47,48,49,50,51,52,18*0,
     *            53,54,55,56,57,58,59,60,61,62,63,64,65,17*0,
     *            66,67,68,69,70,71,72,73,74,75,76,77,78,79,16*0,
     *            80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,15*0,
     *            95,96,97,98,99,100,101,102,103,104,105,106,107,
     *            108,109,110,111,13*0,
     *            112,113,114,115,116,117,118,119,120,121,122,123,
     *            124,125,126,127,128,129,130,11*0,
     *            131,132,133,134,135,136,137,138,139,140,141,142,
     *            143,144,145,146,147,148,149,150,151,9*0,
     *            152,153,154,155,156,157,158,159,160,161,162,163,
     *            164,165,166,167,168,169,170,171,172,173,174,175,
     *            176,5*0,
     *            177,178,179,180,181,182,183,184,185,186,187,188,
     *            189,190,191,192,193,194,195,196,197,198,199,200,
     *            201,202,4*0,
     *            203,204,205,206,207,208,209,210,211,212,213,214,
     *            215,216,217,218,219,220,221,222,223,224,225,226,
     *            227,228,229,3*0,
     *            230,231,232,233,234,235,236,237,238,239,240,241,
     *            242,243,244,245,246,247,248,249,250,251,252,253,
     *            254,255,256,257,258,0/
      data gg/2.,29*0.,
     *        2.,1.,28*0.,
     *        2.,1.,2.,1.,6.,9.,24*0.,
     *        2.,1.,2.,1.,6.,9.,4.,23*0.,
     *        2.,1.,2.,1.,6.,9.,4.,9.,22*0.,
     *        2.,1.,2.,1.,6.,9.,4.,9.,6.,1.,20*0.,
     *        2.,1.,2.,1.,6.,9.,4.,9.,6.,1.,2.,19*0.,
     *        2.,1.,2.,1.,6.,9.,4.,9.,6.,1.,2.,1.,18*0.,
     *        2.,1.,2.,1.,6.,9.,4.,9.,6.,1.,2.,1.,6.,17*0.,
     *        2.,1.,2.,1.,6.,9.,4.,9.,6.,1.,2.,1.,6.,9.,16*0.,
     *        2.,1.,2.,1.,6.,9.,4.,9.,6.,1.,2.,1.,6.,9.,4.,9.,14*0.,
     *        2.,1.,2.,1.,6.,9.,4.,9.,6.,1.,2.,1.,6.,9.,4.,9.,6.,1.,
     *        12*0.,2.,1.,2.,1.,6.,9.,4.,9.,6.,1.,2.,1.,6.,9.,4.,9.,
     *        6.,1.,2.,1.,10*0.,2.,1.,2.,1.,6.,9.,4.,9.,6.,1.,2.,1.,
     *        6.,9.,4.,9.,6.,1.,10.,21.,28.,25.,6.,7.,6*0.,
     *        2.,1.,2.,1.,6.,9.,4.,9.,6.,1.,2.,1.,6.,9.,4.,9.,
     *           6.,1.,10.,21.,28.,25.,6.,7.,6.,5*0.,
     *        2.,1.,2.,1.,6.,9.,4.,9.,6.,1.,2.,1.,6.,9.,4.,9.,
     *           6.,1.,10.,21.,28.,25.,6.,25.,30.,25.,4*0.,
     *        2.,1.,2.,1.,6.,9.,4.,9.,6.,1.,2.,1.,6.,9.,4.,9.,
     *           6.,1.,10.,21.,28.,25.,6.,25.,28.,21.,10.,21.,0.,0./
      data uu(1,1)/109.6787/
      data uu(1,2)/198.3108/
      data uu(2,2)/438.9089/
      data u6/90.82,196.665,386.241,520.178,3162.395,3952.061/
      data u7/117.225,238.751,382.704,624.866,789.537,4452.758,5380.089/
      data u8/109.837,283.24,443.086,624.384,918.657,1114.008,5963.135,
     *        7028.393/
      data u10/173.93,330.391,511.8,783.3,1018.,1273.8,1671.792,
     *         1928.462,9645.005,10986.876/
      data u11/41.449,381.395,577.8,797.8,1116.2,1388.5,1681.5,2130.8,
     *         2418.7,11817.061,13297.676/
      data u12/61.671,121.268,646.41,881.1,1139.4,1504.3,1814.3,2144.7,
     *         2645.2,2964.4,14210.261,15829.951/
      data u13/48.278,151.86,229.446,967.8,1239.8,1536.3,1947.3,2295.4,
     *         2663.4,3214.8,3565.6,16825.022,18584.138/
      data u14/65.748,131.838,270.139,364.093,1345.1,1653.9,1988.4,
     *         2445.3,2831.9,3237.8,3839.8,4222.4,19661.693,21560.63/
      data u16/83.558,188.2,280.9,381.541,586.2,710.184,2265.9,2647.4,
     *         3057.7,3606.1,4071.4,4554.3,5255.9,5703.6,26002.663,
     *         28182.535/
      data u18/127.11,222.848,328.6,482.4,605.1,734.04,1002.73,1157.08,
     *         3407.3,3860.9,4347.,4986.6,5533.8,6095.5,6894.2,7404.4,
     *         33237.173,35699.936/
      data u20/49.306,95.752,410.642,542.6,681.6,877.4,1026.,1187.6,
     *         1520.64,1704.047,4774.,5301.,5861.,6595.,7215.,7860.,
     *         8770.,9338.,41366.,44177.41/
      data u24/54.576,132.966,249.7,396.5,560.2,731.02,1291.9,1490.,
     *         1688.,1971.,2184.,2404.,2862.,3098.52,8151.,8850.,
     *         9560.,10480.,11260.,12070.,13180.,13882.,60344.,63675.9/
      data u25/59.959,126.145,271.55,413.,584.,771.1,961.44,1569.,
     *         1789.,2003.,2307.,2536.,2771.,3250.,3509.82,9152.,
     *         9872.,10620.,11590.,12410.,13260.,14420.,15162.,
     *         65660.,69137.4/
      data u26/63.737,130.563,247.22,442.,605.,799.,1008.,1218.38,
     *         1884.,2114.,2341.,2668.,2912.,3163.,3686.,3946.82,
     *         10180.,10985.,11850.,12708.,13620.,14510.,15797.,
     *         16500.,71203.,74829.6/
      data u28/61.6,146.542,283.8,443.,613.5,870.,1070.,1310.,1560.,
     *         1812.,2589.,2840.,3100.,3470.,3740.,4020.,4606.,
     *         4896.2,12430.,13290.,14160.,15280.,16220.,17190.,
     *         18510.,19351.,82984.,86909.4/
c
c     If the routine is called with IAT=0, initialization:
c     read the ionization fractions form the Opacity Project tables,
c     and evaluate partition functions (assuming that the partition
c     function for the highest ion is eqaul to the statistical weight of
c     the ground state.
c     The table of partition functions, in the OP temperatures and
c     electron densities, are stored in the array PFOP
c
C     Read is done in a loop over the OP species
c
      fra=1.
      if(iat.gt.0) go to 50
      do 40 iatnum=1,28
      if(idat(iatnum).eq.0) go to 40
c
      g0(iatnum+1)=1.
      indxa(iatnum+1)=indxat(iatnum+1,idat(iatnum))
      do i=1,iatnum
        ig0=iatnum-i+1
        g0(ig0)=gg(i,idat(iatnum))
        indxa(i)=indxat(i,idat(iatnum))
        u0(i)=uu(i,idat(iatnum))*1000.
      enddo
c
c     initializion of partition functions by the statistical weights of
c     the ground state
c
      do i=1,iatnum+1
         indx=indxa(i)
         pf0=g0(i)
         do it=1,mtemp
            do ieind=1,melec
               pfop(it,ieind,indx)=pf0
            enddo
         enddo
      enddo
c
      if(iatnum.eq.1) open(inp,file='ioniz.dat',status='old')
      read(inp,*)
      read(inp,*) it0,it1,itstp
      ntt=(it1-it0)/itstp+1
c
      do 30 it=1,ntt
         read(inp,*) itt,ie0,ie1,iestp
         itemp(it)=itt
         net=(ie1-ie0)/iestp+1
         t=exp(2.3025851*0.025*itt)
         safac0=sqrt(t)*t/2.07d-16
         tkcm=0.69496*t
         do 20 ie=1,net
            read(inp,601) iee,ion0,ion1,
     *                    (ioo(i),frac0(i),i=ion0+2,min(ion1,ion0+3)+2)
            ane=exp(2.3025851*0.25*iee)
            safac=safac0/ane
            nio=ion1-ion0
            if(nio.ge.3) then
               nlin=nio/4
               do ilin=1,nlin
                  read(inp,602) (ioo(i),frac0(i),
     *                 i=ion0+4*ilin+2,min(ion1,ion0+4*ilin+3)+2)
               end do
            end if
            ieind=iee/2
            ion0=ion0+2
            ion1=ion1+2
            do 10 ionn=ion0,ion1
              iation=iatnum+2-ionn
              if(ionn.lt.iatnum+2) then
                 if(ionn.eq.ion0) then
                    z0(ionn)=g0(iation)
                  else
                    z0(ionn)=frac0(ionn)/frac0(ionn-1)*safac*z0(ionn-1)
                    z0(ionn)=z0(ionn)*exp(-u0(iation)/tkcm)
                 endif
                 pfop(it,ieind,indxa(iation))=z0(ionn)
                 frac(it,ieind,indxa(iation))=frac0(ionn)
              else
                 u0hm=6090.5
                 z0hm=frac0(ionn)/frac0(ionn-1)*safac
                 z0hm=z0hm*exp(-u0hm/tkcm)
                 pfophm(it,ieind)=z0hm
              end if
   10       continue
   20    continue
   30 continue
   40 continue
  601 format(3i4,2x,4(i4,1x,e9.3))
  602 format(14x,4(i4,1x,e9.3))
      return
C
C     ----------------------------------------------------------
C
C     If the routine is called with IAT>0, evaluate the partition
C     function of atom IAT, ion ION,
C     for temperature T and electron density ANE
C     Evaluation is done by interpolation from previously computed
C     Opacity Project values
C
   50 continue
c
      xt=log10(t)
      kt0=2*int(20.*xt)
      xne=log10(ane)
      kn0=int(2.*xne)
c
      iatnum=iat
      if(idat(iatnum).eq.0) then
         write(6,600) iatnum
         iatnum=-1
  600    format(' data for element no. ',i3,' do not exist')
         return
      end if
      indx=indxat(ion,idat(iatnum))
      if(kt0.lt.itemp(1)) then
         kt1=1
         write(6,611) t
  611    format(' (FRACOP) Extrapol. in T (low)',f7.0)
         go to 120
      endif
      if(kt0.ge.itemp(ntt)) then
         kt1=ntt-1
         write(6,612) t
  612    format(' (FRACOP) Extrapol. in T (high)',f12.0)
         go to 120
      endif
      do 110 it=1,ntt
         if(kt0.eq.itemp(it)) then
            kt1=it
            go to 120
         endif
  110 continue
  120 continue
      if(kn0.lt.1) then
         kn1=1
         go to 130
      endif
      if(kn0.ge.60) then
         kn1=59
         write(6,614) xne
  614    format(' (FRACOP) Extrapol. in Ne (high)',f9.4)
         go to 130
      endif
      kn1=kn0
  130 continue
      xt1=0.025*itemp(kt1)
      dxt=0.05
      at1=(xt-xt1)/dxt
      xn1=0.5*kn1
      dxn=0.5
      an1=(xne-xn1)/dxn
      x11=pfop(kt1,kn1,indx)
      x21=pfop(kt1+1,kn1,indx)
      x12=pfop(kt1,kn1+1,indx)
      x22=pfop(kt1+1,kn1+1,indx)
      x1221=x11*x21*x12*x22
      if(x1221.eq.0.) then
         xx1=x11+at1*(x21-x11)
         xx2=x12+at1*(x22-x12)
         rrx=xx1+an1*(xx2-xx1)
       else
         x11=log10(x11)
         x21=log10(x21)
         x12=log10(x12)
         x22=log10(x22)
         xx1=x11+at1*(x21-x11)
         xx2=x12+at1*(x22-x12)
         rrx=xx1+an1*(xx2-xx1)
         rrx=exp(2.3025851*rrx)
      endif
      pf=rrx
c
      return
      end
c 
c ******************************************************************
c
C
      SUBROUTINE PFHEAV(IIZ,JNION,MODE,t,ane,u)
C     =========================================
C
c     subset of kurucz's pfsaha for Z>28.
c     removed code for Z<28; crp- 28 aug, 1995
C     EDITED 27 JULY 1994 BY GMW - REPLACED PT III PF COEFF. AND IP             
C     MODE 3 RETURNS PARTITION FUNCTION   
C                                      
      INCLUDE 'IMPLIC.FOR' 
      REAL*8 IP
      PARAMETER (DEBCON=1./2.8965E-18,
     *           TVCON=8.6171E-5,
     *           HIONEV=13.595,
     *           ONE=1.,
     *           HALF=0.5,
     *           THIRD=1./3.,
     *           X18=1./18.,
     *           X120=1./120.,
     *           T211=2000./11.)
c                                                 
C     DIMENSION F(6),
      DIMENSION IP(6),PART(6),POTLO(6)
C     DIMENSION FSAVE(6)
      DIMENSION SCALE(4)                                                        
      DIMENSION NNN(6*218)                                                      
      DIMENSION NNN16(54),NNN17(54),NNN18(54),NNN19(54),NNN20(54)               
      DIMENSION NNN21(54),NNN22(54),NNN23(54),NNN24(54),NNN25(54)               
      DIMENSION NNN26(54),NNN27(54),NNN28(54),NNN29(54),NNN30(54)               
      DIMENSION NNN31(54),NNN32(54),NNN33(54),NNN34(54),NNN35(54)               
      DIMENSION NNN36(54),NNN37(54),NNN38(54),NNN39(54),NNN40(12)               
      EQUIVALENCE (NNN( 811-810),NNN16(1))                     
      EQUIVALENCE (NNN( 865-810),NNN17(1)),(NNN( 919-810),NNN18(1))
      EQUIVALENCE (NNN( 973-810),NNN19(1)),(NNN(1027-810),NNN20(1))
      EQUIVALENCE (NNN(1081-810),NNN21(1)),(NNN(1135-810),NNN22(1))
      EQUIVALENCE (NNN(1189-810),NNN23(1)),(NNN(1243-810),NNN24(1))
      EQUIVALENCE (NNN(1297-810),NNN25(1)),(NNN(1351-810),NNN26(1)) 
      EQUIVALENCE (NNN(1405-810),NNN27(1)),(NNN(1459-810),NNN28(1)) 
      EQUIVALENCE (NNN(1513-810),NNN29(1)),(NNN(1567-810),NNN30(1))
      EQUIVALENCE (NNN(1621-810),NNN31(1)),(NNN(1675-810),NNN32(1)) 
      EQUIVALENCE (NNN(1729-810),NNN33(1)),(NNN(1783-810),NNN34(1)) 
      EQUIVALENCE (NNN(1837-810),NNN35(1)),(NNN(1891-810),NNN36(1)) 
      EQUIVALENCE (NNN(1945-810),NNN37(1)),(NNN(1999-810),NNN38(1)) 
      EQUIVALENCE (NNN(2053-810),NNN39(1)),(NNN(2107-810),NNN40(1)) 
C      ( 1)( 2)   ( 3)( 4)   ( 5)( 6)   ( 7)( 8)   ( 9)(10)   ( IP ) G  REF     
      DATA NNN16/                                                               
     1 227027622, 306233052, 356839222, 446052912, 652382292,   763314, 
     2 108416342, 222428472, 353944332, 577378932, 110314303,  1814900, 
     3 198724282, 293236452, 468362702,  86511123, 136016073,  3516000, 
     4 279836622, 461857562, 720693022, 124915873, 192522633,  5600000, 
     5 262136422, 501167232,  87911303, 138916483, 190721673,  7900000, 
     6 201620781, 231026761, 314737361, 450555381, 692386911,   772301, 
     7 109415761, 247938311,  58910042, 190937022,  68311693,  2028903, 
     8 897195961, 107212972, 165021182, 260230862, 356940532,  3682900, 
     9 100010001, 100410231, 108712611, 167124841, 388460411,   939102/ 
      DATA NNN17/                                                               
     1 200020021, 201620761, 223726341, 351352061,  80812472,  1796001, 
     2 100610471, 122617301, 300566361, 149924112, 332342352,  3970000, 
     3 403245601, 493151431, 529654331, 559358091, 611065171,   600000, 
     4  99710051, 104511541, 135016501, 208226431, 321837921,  2050900, 
     5 199820071, 204521391, 229124761, 266028451, 302932131,  3070000, 
     6 502665261, 755183501, 901496201, 102410942, 117912812,   787900, 
     7 422848161, 512153401, 557458941, 636270361, 794489061,  1593000, 
     8 100010261, 114613921, 175221251, 249828711, 324436181,  3421000, 
     9 403143241, 491856701, 649173781, 840396751, 113013392,   981000/ 
      DATA NNN18/                                                               
     1 593676641, 884697521, 105911572, 129515012, 180322212,  1858700, 
     2 484470541,  91510972, 125614082, 157017612, 199722912,  2829900, 
     3 630172361, 799686381, 919797221, 102810942, 117712832,   975000, 
     4 438055511, 691582151,  94510732, 121413672, 152016732,  2150000, 
     5 651982921,  94610382, 113212492, 139515462, 169718482,  3200000, 
     6 437347431, 498951671, 538559501,  74710812, 169126672,  1183910, 
     7 705183611,  93510092, 111614162, 222932532, 427652992,  2160000, 
     8 510869921,  87410312, 123116552, 236530712, 377744832,  3590000, 
     9 100010001, 100010051, 105012781, 198535971,  65911422,  1399507/ 
      DATA NNN19/                                                               
     1 461049811, 522254261, 609088131, 168935052,  68612253,  2455908, 
     2 759990901, 101911142, 129017782, 302856642,  99414333,  3690000, 
     3 200020011, 200720361, 211523021, 269434141, 459163351,   417502, 
     4 100010001, 100110321, 129524961,  61014202, 291753192,  2750004, 
     5 473650891, 533156051,  66810932, 232950852,  99915303,  4000000, 
     6 100110041, 104111741, 146019721, 281941411, 607785251,   569202, 
     7 202621931, 255331271, 384347931, 624085761, 122417632,  1102600, 
     8 100010001, 100110321, 129524961,  61014202, 291753192,  4300000, 
     9 791587851, 100012192, 155119942, 254031782, 389946932,   637900/ 
      DATA NNN20/                                                               
     1 118217102, 220827002, 319036792, 416646512, 513256072,  1223000, 
     2  92510012, 104710862, 112311612, 120212472, 132814282,  2050000, 
     3 141320802, 291439702, 531170262,  92712273, 162521053,   684000, 
     4 354454352, 724689652, 107212643, 148517093, 193321573,  1312900, 
     5 209727032, 324537052, 415446282, 510255752, 604965222,  2298000, 
     6 256636022, 465759302, 749693962, 116514243, 171520333,   687900, 
     7 335157222,  84511463, 147718363, 221826083, 299933893,  1431900, 
     8 223725352, 280830972, 340937362, 406844002, 473150632,  2503900, 
     9 703972941,  82610822, 154822682, 327244912, 571469372,   709900/ 
      DATA NNN21/                                                               
     1  75714552, 274347322, 718897632, 123414913, 174920063,  1614900, 
     2 267645462, 669890262, 115514323, 173620673, 242528083,  2714900, 
     3  90613732, 184823562, 291735332, 419949102, 565764332,   728000, 
     4 131318312, 227126932, 311735452, 397644072, 483852692,  1525900, 
     5 204721673, 234725733, 284031463, 348738613, 426546943,  3000000, 
     6 176824122, 318941082, 515263202, 761790472, 106112303,   736400, 
     7 221934642, 501968372,  88911173, 136316243, 189221613,  1675900, 
     8 210622722, 241025422, 267928262, 297731272, 327834282,  2846000, 
     9 148520202, 255230902, 364942462, 489656082, 638872352,   746000/ 
      DATA NNN22/                                                               
     1 153421292, 288137912, 484660322, 720187062, 101011483,  1807000, 
     2 254537212, 492362292, 770592182, 107312243, 137615273,  3104900, 
     3 115919651, 320746011, 607576761,  95011642, 141817172,   832900, 
     4 755087211, 105913442, 173122222, 282034722, 412247732,  1941900, 
     5 180223462, 289735212, 414247632, 538460052, 662672472,  3292000, 
     6 200020001, 200220141, 206422141, 257633021, 455164681,   757403, 
     7 100810581, 125817401, 260641031,  66210072, 135316982,  2148000, 
     8 795887491,  97711762, 156620252, 248329422, 340038582,  3481900, 
     9 100010001, 100410241, 109212891, 176827421, 444268771,   899003/ 
      DATA NNN23/                                                               
     1 200020021, 201720921, 233329881, 451475371, 127520782,  1690301, 
     2 100310281, 114815371, 246138311, 519265531, 791492761,  3747000, 
     3 252431921, 368440461, 433746521, 512259221, 723389021,   578400, 
     4 100110071, 104611651, 146118581, 225426511, 304734431,  1886000, 
     5 200120111, 205021611, 243628031, 317035371, 390442701,  2802900, 
     6 232637101, 488058571, 669074381, 816189091,  97210632,   734200, 
     7 286335941, 408144471, 479351961, 571862901, 686274341,  1462700, 
     8 100010251, 114013811, 175321601, 256829751, 338337901,  3049000, 
     9 404043481, 494656811, 646772781, 813490751, 101411372,   863900/ 
      DATA NNN24/                                                               
     1 303147981, 618472951, 827392621, 103711702, 131214532,  1650000, 
     2 313037601, 429347901, 536260591, 689477591, 862494881,  2529900, 
     3 526258801, 657372351, 784284071, 897095741, 102711082,   900900, 
     4 440855541, 686481251,  93810792, 125414792, 176321132,  1860000, 
     5 349054751, 699883081,  96611302, 134216202, 197724212,  2800000, 
     6 405342041, 438645621, 475751071, 587974491, 102214572,  1045404, 
     7 568567471, 773485861,  94510362, 112712182, 130914002,  1909000, 
     8 514269581,  86910562, 130716652, 215327742, 351843662,  3200000, 
     9 100010001, 100010091, 109515351, 291060661, 119621482,  1212716/ 
      DATA NNN25/                                                               
     1 414844131, 465649111, 538464651,  87112232, 158019362,  2120000, 
     2 615475101, 867797531, 112213462, 157618062, 203622662,  3209900, 
     3 200020001, 201020501, 215623871, 283536181, 462756261,   389300, 
     4 100010001, 100310371, 119016501, 269146361,  77912412,  2510000, 
     5 424445601, 481750061, 516953311, 549356551, 581759791,  3500000, 
     6 101210791, 135119351, 282340571, 574580391, 111015062,   521002, 
     7 262638611, 504160621, 698579371,  91010692, 129115952,  1000000, 
     8 100010001, 100310351, 118416321, 264945521,  76512182,  3700000, 
     9  71111992, 172323592, 312540402, 510763182, 765791012,   558000/ 
      DATA NNN26/                                                               
     1 204529582, 383647882, 582469262, 807992692, 104911723,  1106000, 
     2  94712552, 148416582, 179819212, 203621522, 227424042,  1916900, 
     3 295959132, 103515693, 215527593, 335939413, 449650223,   565000, 
     4  79718153, 289639443, 495159253, 686877533, 863794813,  1085000, 
     5 298640242, 475053692, 596965912, 725379692, 872094692,  2008000, 
     6 460693672, 158523823, 327242303, 519661563, 709379783,   541900, 
     7 455480232, 114014653, 178521013, 240927073, 299232633,  1055000, 
     8  46410533, 183826893, 354443773, 518459633, 674375243,  2320000, 
     9 139623042, 364860002,  96114603, 209828633, 373446973,   549000/ 
      DATA NNN27/                                                               
     1 460493692, 158523823, 327142303, 519661563, 709279783,  1073000, 
     2 455480232, 114014653, 178521013, 240927073, 299232633,  2000000, 
     3 131720482, 280535692, 441254492, 676583972, 103412583,   555000, 
     4 139623042, 364860002,  96114603, 209828633, 373446973,  1089900, 
     5 460493682, 158523823, 327142303, 519661563, 709279783,  2000000, 
     6  92915672, 222431062, 444763802,  89612173, 159520253,   562900, 
     7 315059662,  97114563, 204627093, 342541693, 490556383,  1106900, 
     8 269037812, 520270372,  91111273, 133915483, 172719093,  2000000, 
     9 800080571, 851699301, 127617362, 240433032, 444958442,   568000/ 
      DATA NNN28/                                                               
     1 125416052, 211828182, 375549622, 644381732, 101112213,  1125000, 
     2 800080571, 851699301, 127617362, 240433032, 444958442,  2000000, 
     3 240432982, 427555202, 708489962, 112613853, 167319843,   615900, 
     4 534793262, 139219123, 247730843, 371043333, 495055893,  1210000, 
     5 364145232, 514756362, 604864112, 673870372, 732276072,  2000000, 
     6 480767202,  89011393, 144118243, 230028753, 354142883,   584900, 
     7 480767192,  89011393, 144118243, 230028753, 354142883,  1151900, 
     8 480767202,  89011393, 144118243, 230028753, 354142883,  2000000, 
     9 343147532, 645887152, 115314793, 183322063, 257729373,   593000/ 
      DATA NNN29/                                                               
     1 343147532, 645887142, 115314793, 183322063, 257729373,  1167000, 
     2 343147532, 645887142, 115314793, 183322063, 257729373,  2000000, 
     3 222635002, 542276772, 100312353, 145716713, 187020703,   602000, 
     4 222635002, 542276772, 100312353, 145716713, 187020703,  1180000, 
     5 222635002, 542276772, 100312353, 145716713, 187020703,  2000000, 
     6 133715382, 209130152, 429859382,  79410293, 129815983,   609900, 
     7 265934782, 497877532, 120517733, 245032063, 400448073,  1193000, 
     8 265934782, 497877532, 120517733, 245032063, 400448073,  2000000, 
     9 800381111,  87510702, 147621462, 310343462, 585475982,   618000/ 
      DATA NNN30/                                                               
     1 156718872, 279244452, 678196342, 128316243, 197823443,  1205000, 
     2  93517192, 364666132, 103414613, 192624193, 293334613,  2370000, 
     3 100010011, 101310651, 118613951, 169120661, 250629971,   625000, 
     4 200120901, 270345231,  81714042, 223533112, 461959862,  1217000, 
     5 100312561, 250851931,  91914182, 198626022, 323638692,  2000000, 
     6 514664441, 759086851,  99211442, 133315612, 182721252,   609900, 
     7 125924831, 438667801,  98714112, 199727872, 380850742,  1389900, 
     8 323948621, 661297271, 158626482, 426865032,  93712843,  1900000, 
     9 659294081, 128016962, 222528952, 372047062, 585171462,   700000/ 
      DATA NNN31/                                                               
     1  99117882, 274638812, 520867322,  84410313, 123314453,  1489900, 
     2 187427702, 343739872, 448049452, 539358282, 625266642,  2329900, 
     3  65210892, 171325762, 373552252, 705192012, 116414343,   787900, 
     4 192837842, 600784802, 111113823, 165419233, 218524383,  1620000, 
     5  99117872, 274638812, 520867312,  84410313, 123314453,  2400000, 
     6 398981651, 130019172, 273438022, 516168382,  88411163,   797900, 
     7 131429482, 523279952, 111414623, 183422233, 262130233,  1770000, 
     8 192837842, 600784792, 111113823, 165419233, 218524383,  2500000, 
     9 600963001,  75910412, 150121572, 301940972, 539168952,   787000/ 
      DATA NNN32/                                                               
     1  73710852, 190731262, 464964142,  83810503, 127315053,  1660000, 
     2 131429482, 523279952, 111414623, 183422233, 262130233,  2600000, 
     3 110815502, 216829732, 398752322, 672484682, 104612673,   850000, 
     4 168225972, 362046562, 566766422, 757484612,  93010103,  1700000, 
     5  73710852, 190731262, 464964142,  83810503, 127315053,  2700000, 
     6 129117892, 239430882, 388748292, 596173252,  89510843,   910000, 
     7 110815502, 216829732, 398752322, 672484682, 104612673,  2000000, 
     8 168225972, 362046562, 566766422, 757484612,  93010103,  2800000, 
     9 158918512, 207523002, 254328242, 316335762, 407246582,   900000/ 
      DATA NNN33/                                                               
     1  98115462, 224930742, 401150612, 623475412,  89910583,  1855900, 
C    2 110815502, 216829732, 398752322, 672484682, 104612673,  2900000, 
     2 146323292, 354651802,  74810923, 161723953, 348749363,  3322700, 
     3 203222611, 265731251, 364042301, 494958601, 702084731,   922000, 
     4 120521331, 357753801,  75310062, 130516572, 206925452,  2050000, 
     5 651780821, 108814772, 195925252, 316338622, 460853882,  3000000, 
     6 100010001, 100110111, 105211851, 152122101, 341552811,  1043002, 
     7 200320211, 210023021, 268834231, 480472341, 111416912,  1875000, 
     8 104012871, 186129471, 458664151,  82410072, 119013732,  3420000, 
     9 200420711, 222424271, 265429161, 325637371, 442853911,   610500/ 
      DATA NNN34/                                                               
     1 100010021, 101910801, 121414641, 189525811, 358949721,  2041900, 
     2 200020311, 216624611, 296337451, 489064791,  85711212,  2979900, 
     3 103411711, 147819101, 244331781, 434862751,  93113762,   741404, 
     4 204122231, 248227841, 311535621, 429153941, 651976431,  1502800, 
     5 100210131, 106812201, 154522671, 381665951,  95512512,  3192900, 
     6 400140351, 416944121, 474851591, 564362181, 690477231,   728700, 
     7 106814451, 204427341, 350744811, 586879131, 108314772,  1667900, 
     8 205523051, 264830231, 345439921, 469156001, 675281671,  2555900, 
     9 500950661, 518153561, 559058941, 628968071, 748483501,   843000/ 
      DATA NNN35/                                                               
     1 443756241, 696282451,  95411012, 128615262, 182922012,  1900000, 
     2 336953201, 682481011,  93810882, 127915272, 184622442,  2700000, 
     3 402841621, 431544771, 463148311, 520059491, 734896851,   930000, 
     4 576168741, 788387631,  96910642, 116012552, 135014462,  2000000, 
     5 490265341, 812797201, 116614322, 179622692, 285035302,  2900000, 
     6 100010001, 100010031, 102311051, 133018071, 264539391,  1074500, 
     7 402841621, 431544771, 463148311, 520059491, 734996851,  2000000, 
     8 576168741, 788387631,  96910642, 116012552, 135014462,  3000000, 
     9 200020011, 201220591, 218124481, 296538611, 488859141,   400000/ 
      DATA NNN36/                                                               
     1 100010001, 100010031, 102311051, 133018071, 264539401,  2200000, 
     2 421645151, 477449611, 511852711, 542455761, 572958821,  3300000, 
     3 100010041, 105212131, 153220271, 270435641, 460258111,   527600, 
     4 201221791, 258131471, 381645781, 546365131, 777592781,  1014400, 
     5 100010001, 100010031, 102311051, 133018071, 264539391,  3400000, 
     6 510064491,  82710872, 142718412, 232328712, 348341572,   690000, 
     7 228951571,  88513232, 183324132, 305537492, 448152402,  1210000, 
     8 723989131, 103511752, 130814352, 155416652, 177018682,  2000000, 
     9 620099241, 162725772, 391457072,  80110833, 141818023,   600000/ 
      DATA NNN37/                                                               
     1 620099241, 162725772, 391457072,  80110833, 141818023,  1200000, 
     2 620099251, 162725772, 391457072,  80110833, 141818023,  2000000, 
     3 347877992, 129318323, 240730533, 380546863, 570368573,   600000, 
     4 347877992, 129318323, 240730533, 380546863, 570368573,  1200000, 
     5 347777992, 129318323, 240730533, 380546863, 570368573,  2000000, 
     6 209530092, 450866762,  96613623, 186524763, 318839893,   600000, 
     7 209530092, 450866762,  96613623, 186524763, 318839893,  1200000, 
     8 209530092, 450866762,  96613623, 186524763, 318839893,  2000000, 
     9 209530092, 450866762,  96613623, 186524763, 318839893,   600000/ 
      DATA NNN38/                                                               
     1 209530092, 450866762,  96613623, 186524763, 318839893,  1200000, 
     2 209530092, 450866762,  96613623, 186524763, 318839893,  2000000, 
     3 209530092, 450866762,  96613623, 186524763, 318839893,   600000, 
     4 209530092, 450866762,  96613623, 186524763, 318839893,  1200000, 
     5 209530092, 450866762,  96613623, 186524763, 318839893,  2000000, 
     6 209530092, 450866762,  96613623, 186524763, 318839893,   600000, 
     7 209530092, 450866762,  96613623, 186524763, 318839893,  1200000, 
     8 209530092, 450866762,  96613623, 186524763, 318839893,  2000000, 
     9 209530092, 450866762,  96613623, 186524763, 318839893,   600000/ 
      DATA NNN39/                                                               
     1 209530092, 450866762,  96613623, 186524763, 318839893,  1200000, 
     2 209530092, 450866762,  96613623, 186524763, 318839893,  2000000, 
     3 209530092, 450866762,  96613623, 186524763, 318839893,   600000, 
     4 209530092, 450866762,  96613623, 186524763, 318839893,  1200000, 
     5 209530092, 450866762,  96613623, 186524763, 318839893,  2000000, 
     6 209530092, 450866762,  96613623, 186524763, 318839893,   600000, 
     7 209530092, 450866762,  96613623, 186524763, 318839893,  1200000, 
     8 209530092, 450866762,  96613623, 186524763, 318839893,  2000000, 
     9 209530092, 450866762,  96613623, 186524763, 318839893,   600000/ 
      DATA NNN40/                                                               
     1 209530092, 450866762,  96613623, 186524763, 318839893,  1200000, 
     2 209530092, 450866762,  96613623, 186524763, 318839893,  2000000/ 
      DATA SCALE/.001,.01,.1,1./    
C
      if(mode.lt.0) return
      tk=1.38054d-16*t
      tv=8.6171d-5*t
C     LOWERING OF THE IONIZATION POTENTIAL IN VOLTS FOR UNIT ZEFF               
      CHARGE=ANE*2.                                                         
      DEBYE=SQRT(TK*DEBCON/CHARGE)                                       
C     DEBYE=SQRT(TK/12.5664/4.801E-10**2/CHARGE)                             
      POTLOW=MIN(1.D0,1.44E-7/DEBYE)                                            
      IF(IIZ.LE.28)then
         write(6,*) 'Error, routine PFHEAV for Z.GE.28 only'
         stop23
       endif
c removed elements with z<28
      if(iiz.eq.28) n=1                                                    
      IF(IIZ.GT.28) N=3*IIZ+54-135                                                     
      IF(IIZ.eq.28) NIONS=4
      IF(IIZ.GT.28) NIONS=3                                                       
      NION2=MIN0(JNION+2,NIONS)                                                  
      N=N-1                                                                     
C                                                                               
      DO 18 ION=1,NION2                                                         
      Z=ION                                                                     
      POTLO(ION)=POTLOW*Z                                                       
      N=N+1 
      nnn6n=nnn(6+6*(N-1))  
c     nnn6n=nnn(6,n)                                                                  
      NNN100=NNN6N/100 
      XN1= NNN100                                                     
      IP(ION)=XN1*1.e-3                                               
      IG=NNN6N-NNN100*100  
      GGG=IG                                                   
      T2000=IP(ION)*T211                                                  
      IT=MAX0(1,MIN0(9, INT(T/T2000-HALF))) 
      XIT=IT                                   
      DT=T/T2000-XIT-HALF                                                
      PMIN=ONE                                                                   
      I=(IT+1)/2                                                                
      nnnin=nnn(i+6*(N-1))  
c     nnnin=nnn(i,n)                                                                  
      K1=NNNIN/100000                                                        
      K2=NNNIN-K1*100000                                                     
      K3=K2/10                                                                  
      xk1=k1                                                
      xk3=k3
      KSCALE=K2-K3*10                                                           
      IF(MOD(IT,2).EQ.0)GO TO 12                                                
      P1=XK1*SCALE(KSCALE)                                                
      P2=XK3*SCALE(KSCALE)                                                
      IF(DT.GE.0.)GO TO 13                                                      
      IF(KSCALE.GT.1)GO TO 13                                                   
      KP1=int(P1)                                                                   
      IF(KP1.NE. INT(P2+.5))GO TO 13                                            
      PMIN=KP1                                                                  
      GO TO 13                                                                  
   12 continue
      xk3=k3
      P1=XK3*SCALE(KSCALE)                                                
      nnni1n=nnn(i+1+6*(N-1))  
c     nnni1n=nnn(i+1,n)                                                                  
      K1=NNNI1N/100000                                                      
      KSCALE=MOD(NNNI1N,10) 
      xk1=k1                                                
      P2=XK1*SCALE(KSCALE)                                                
   13 PART(ION)= MAX (PMIN,P1+(P2-P1)*DT)                                       
      IF(GGG.EQ.0..OR.POTLO(ION).LT..1.OR.T.LT.T2000*4.)GO TO 18               
      IF(T.GT.(T2000*11.)) TV=(T2000*11.)*TVCON                          
      D1=.1/TV                                                                  
      D2=POTLO(ION)/TV  
      DX=SQRT(HIONEV*Z*Z/TV/D2)**3                                                       
      PART(ION)=PART(ION)+GGG*EXP(-IP(ION)/TV)*       
     *          (DX*(THIRD+(ONE-(HALF+(X18+D2*X120)*D2)*D2)*D2)-                              
     *           DX*(THIRD+(ONE-(HALF+(X18+D1*X120)*D1)*D1)*D1))                              
   18 CONTINUE                                                                  
      u=part(jnion)
      RETURN                                                                    
      END                                                                       
C
C
C ********************************************************************
C
C
      FUNCTION XK2DOP(TAU)
C     ====================
C
C     KERNEL FUNCTION K2 
C     AFTER  HUMMER,  1981, J.Q.S.R.T. 26, 187
C
      INCLUDE 'IMPLIC.FOR'
      PARAMETER (PI2SQ=2.506628275D0, PISQ=1.772453851D0, UN=1.D0,
     *   A1= -1.117897000D-1,   A2= -1.249099917D-1, 
     *   A3= -9.136358767D-3,   A4= -3.370280896D-4,
     *   B1=  1.566124168D-1,   B2=  9.013261660D-3,   
     *   B3=  1.908481163D-4,   B4= -1.547417750D-7, 
     *   B5= -6.657439727D-9,
     *   C1=  1.915049608D01,   C2=  1.007986843D02,   
     *   C3=  1.295307533D02,   C4= -3.143372468D01,
     *   D1=  1.968910391D01,   D2=  1.102576321D02,   
     *   D3=  1.694911399D02,   D4= -1.669969409D01,  
     *   D5= -3.666448000D01)
      XK2DOP=UN
      IF(TAU.LE.0.) RETURN
      IF(TAU.GT.11.) GO TO 10
      P=UN+TAU*(A1+TAU*(A2+TAU*(A3+TAU*A4)))
      Q=UN+TAU*(B1+TAU*(B2+TAU*(B3+TAU*(B4+TAU*B5))))
      XK2DOP=TAU/PI2SQ*LOG(TAU/PISQ)+P/Q
      RETURN
   10 X=UN/LOG(TAU/PISQ)
      P=UN+X*(C1+X*(C2+X*(C3+X*C4)))
      Q=UN+X*(D1+X*(D2+X*(D3+X*(D4+X*D5))))
      XK2DOP=P/Q/2.D0/TAU/SQRT(LOG(TAU/PISQ))
      RETURN
      END
C
C
C ********************************************************************
C
C
 
      SUBROUTINE LTEGR
C     ================
C
C     Driving procedure for computing the initial LTE-grey model
C     atmosphere
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'ATOMIC.FOR'
      INCLUDE 'MODELQ.FOR'
C
      COMMON ESEMAT(MLEVEL,MLEVEL),BESE(MLEVEL),
     *       DEPTH(MDEPTH),DEPTH0(MDEPTH),TAU(MDEPTH),TAU0(MDEPTH),
     *       TEMP0(MDEPTH),ELEC0(MDEPTH),DENS0(MDEPTH),DM0(MDEPTH)
C
C ----------------
C Input parameters
C ----------------
C
C     NDEPTH  - number of depth points for evaluating LTE-grey model
C             = 0 - NDEPTH is taken to be ND-1
C     TAUFIR  - Rosseland optical depth in the first depth point
C     TAULAS  - Rosseland optical depth in the last depth point
C     ABROS0  - an estimate of the Rosseland opacity (per gram) at the
C               first depth point
C     TSURF   = 0  - surface temperature and the Hopf function are
C                    evaluated exactly
C             > 0  - value of surface temperature is set to TSURF, and
C                    the Hopf function is assumed to be constant,
C                    corresponding to TSURF
C     ALBAVE  = 0  - wind blanketing is not considered
C             > 0  - wind blanketing is considered; the averaged value
C                    of albedo [precisely the quantity (1+rho)/(1-rho)
C                    in the notation of Hummer (Ap.J. 257, 724, 1982)
C                    see his Eq. (3.1)] is ALBAVE
C     DION0        - initial estimate of the degree of ionization at
C                    the first depth point (=1 for completely ionized;
C                    =1/2 for completely neutral)
C
C --------------------------------------------------------------------
C     IDEPTH  -  mode of determining the mass-depth scale to be used
C                in linearization
C             =  0 - depth scale DM (in g*cm**-2) is evaluated as mass
C                    corresponding to Rosseland optical depths which
C                    are equidistantly spaced in logarithms between
C                    the first point TAUFIR and the last point TAULAS
C                    the last-but-one point is, however, set to
C                    TAULAS-1.
C             =  2 - depth scale DM is evaluated as that corresponding
C                    to input values of Rosseland optical depth -
C                    array TAU0(ID), ID=1,ND
C             =  1 - similar, but now DM is evaluated as mass
C                    corresponding to Rosseland optical depths which
C                    are equidistantly spaced in logarithms between
C                    the first point TAU1 and the last-but-one point
C                    TAU2; the last point is TAUL
C                    (i.e. similar to option 0, now TAU1 and TAUL may
C                    be different from TAUFIR nad TAULAS)
C             =  3 - depth scale DM has already been read in START
C     NCONIT       - number of internal iterations for calculating the
C                    gray model with convection
C             =  0 - and HMIX0>0, then NCONIT is set to 10
C     IPRING       - controls diagnostic output of the LTE-gray model
C                    calculations
C              = 0 - no output
C              = 1 - only final LTE-gray model
C              = 2 - results of all internal iterations
C     IHM     >  0 - negative hydrogen ion considered in particle and
C                    charge conservation in ELDENS
C     IH2     >  0 - hydrogen molecule considered in particle
C                    conservation in ELDENS
C     IH2P    >  0 - ionized hydrogen molecule considered in particle
C                    and charge conservation in ELDENS
C
      IF(NDGREY.EQ.0) THEN
         NDEPTH=ND
       ELSE
         NDEPTH=NDGREY
      END IF
      IF(NDEPTH.GT.MDEPTH) 
     *   CALL QUIT('ndepth.gt.mdepth in LTEGR',ndepth,mdepth)
      IDEPTH=IDGREY
      IF(ALBAVE.GT.0.AND.TSURF.EQ.0.) TSURF=(0.433*ALBAVE)**0.25
      HOPF0=0.
      IF(TSURF.NE.0.) HOPF0=4.D0*TSURF**4/3.D0
      T4=TEFF**4
      ANEREL=(DION0-HALF )/DION0
      IF(NCONIT.EQ.0.AND.HMIX0.GT.0.) NCONIT=10
      IF(ANEREL.LT.1.D-3) ANEREL=1.D-3
      LCHC0=LCHC
      LCHC=.TRUE.
      LTE0=LTE
      LTE=.TRUE.
      IRSPL0=IRSPLT
      IRSPLT=0
      IF(NDEPTH.EQ.0) NDEPTH=ND-1
      ND0=ND
      ND=NDEPTH
      DO 1 I=1,ND0
    1    DM0(I)=DM(I)
C
C     tau(ross) scale - logarithmically equidistant points between
C     input TAUFIR and TAULAS
C
      DML0=LOG(TAUFIR)
      DLGM=(LOG(TAULAS)-DML0)/(NDEPTH-1)
      DO 10 I=1,NDEPTH
         TAU0(I)=DML0+(I-1)*DLGM
         TAU(I)=EXP(TAU0(I))
         TAUROS(I)=TAU(I)
   10 CONTINUE
C
      DPRAD=1.891204931D-15*T4
      if(ifprad.eq.0) dprad=0.
C
      PRAD0=DPRAD/1.732D0
      ABROS=ABROS0
      PLOG1=0.
      PLOG2=0.
      PLOG3=0.
      DPLOG1=0.
      DPLOG2=0.
      IF(IPRING.GT.0) WRITE(6,601)
C
C -------------------------------------------------------------------
C
C     1.part
C     Integration of the hydrostatic equilibrium equation on the
C     tau(ross) scale;
C     basically by a predictor-corrector method
C     (similar to Kurucz's ATLAS code)
C
      DO 50 I=1,NDEPTH
         J=0
         TAUR=TAU(I)
C
C        predictor step
C
c         IF(I.EQ.1) PLOG=LOG(GRAV/ABROS*TAUR)
         IF(I.EQ.1) PLOG=LOG(GRAV/ABROS*TAUR+prad0)
         IF(I.GT.1.AND.I.LE.4) PLOG=PLOG1+DPLOG1
         IF(I.GT.4) PLOG=(3.*PLOG4+8.*DPLOG1-4.*DPLOG2+8.*DPLOG3)/3.
         ERROR=1.
         GO TO 40
C
C        corrector step
C ------ iterate between hydrostatic equilibrium (which determines an
C           increment of the total pressure) and state equations (which
C           determine relevant number densities and then the Rosseland
C           opacity)
C
c   30    IF(I.EQ.1) PNEW=LOG(GRAV/ABROS*TAUR)
   30    IF(I.EQ.1) PNEW=LOG(GRAV/ABROS*TAUR+prad0)
         IF(I.GT.1.AND.I.LE.4)  PNEW=(PLOG+2.*PLOG1+DPLOG+DPLOG1)/3.
         IF(I.GT.4) PNEW=(126.*PLOG1-14.*PLOG3+9.*PLOG4+42.*DPLOG+
     *                    108.*DPLOG1-54.*DPLOG2+24.*DPLOG3)/121.
         ERROR=ABS(PNEW-PLOG)
         PLOG=PNEW
   40    PTOT=EXP(PLOG)
c         P=PTOT-TAUR*DPRAD
         P=PTOT-TAUR*DPRAD-prad0
         J=J+1
         CALL ROSSOP(I,P,TAUR,HOPF0,T4,T,ANE,ABROS)
         DPLOG=GRAV/ABROS*TAUR/PTOT*DLGM
         IF(ERROR.GT.1.D-4.AND.J.LT.10) GO TO 30
C
C ------ end of the iteration loop;
C        set up necessary quantities for the next depth step of the
C        hydrostatic equilibrium
C
         PLOG4=PLOG3
         PLOG3=PLOG2
         PLOG2=PLOG1
         PLOG1=PLOG
         DPLOG3=DPLOG2
         DPLOG2=DPLOG1
         DPLOG1=DPLOG
         TEMP0(I)=T
         ELEC0(I)=ANE
         AN=P/T/BOLK
c         DEPTH(I)=PTOT/GRAV
         DEPTH(I)=(PTOT-prad0)/GRAV
         DM(I)=DEPTH(I)
         DENS0(I)=WMM(I)*(AN-ANE)
         IF(IPRING.GT.0) WRITE(6,602) I,TAU(I),DEPTH(I),T,
     *                                AN,ANE,P,ABROS
c         PTOTAL(I)=PTOT+PRAD0
         PTOTAL(I)=PTOT
         PGS(I)=P
   50 CONTINUE
C
C -------------------------------------------------------------------
C
C     2. Second part - taking into account convection
C
      IF(HMIX0.GT.0.) THEN
         CALL CONTMP
       GO TO 110
      END IF
C
C -------------------------------------------------------------------
C
C     3. Third part
C     Interpolation of the computed model to the depth scale which is
C     going to be used in the subsequent - complete-linearization -
C     part of the model atmosphere construction
C
      ND=ND0
C
C     First option - logarithmically equidistant Rosseland opt.depths
C                    the same first and last depth as in the first part
C
      IF(IDEPTH.EQ.0) THEN
         TAU1=TAUFIR
         TAUL=TAULAS
         TAU2=TAULAS-1.
C
C     Second option - logarithmically equidistant Rosseland opt.depths
C                     the first, last-but-one, and last depths are read
C
        ELSE IF(IDEPTH.EQ.1) THEN
          READ(IBUFF,*) TAU1,TAU2,TAUL
      END IF
C
      IF(IDEPTH.LE.1) THEN
         DML0=LOG(TAU1)
         IF(TAUL.GT.0.) THEN
            DLGM=(LOG(TAU2)-DML0)/(ND-2)
            DO 80 I=1,ND-1
   80          TAU0(I)=DML0+(I-1)*DLGM
            TAU0(ND)=LOG(TAUL)
          ELSE
            DLGM=(LOG(TAU2)-DML0)/(ND-1)
            DO 81 I=1,ND
   81          TAU0(I)=DML0+(I-1)*DLGM
         END IF
       ELSE IF(IDEPTH.EQ.2) THEN
C
C     Third option - prescribed set of Rosseland optical depths
C
         READ(IBUFF,*) (TAU0(I),I=1,ND)
         DO 70 I=1,ND
   70       TAU0(I)=LOG(TAU0(I))
       ELSE if(idepth.eq.3) then
C
C     Fourth option - interpolation to the prescribed mass scale DM
C
         DO 60 I=1,ND
            DM(I)=DM0(I)
   60       DM0(I)=LOG(DM(I))
         CALL INTERP(DEPTH0,TEMP0,DM0,TEMP,NDEPTH,ND,2,0,0)
         CALL INTERP(DEPTH0,ELEC0,DM0,ELEC,NDEPTH,ND,2,0,1)
         CALL INTERP(DEPTH0,DENS0,DM0,DENS,NDEPTH,ND,2,0,1)
      END IF
C
C     in the first three options - interpolation from the previous
C     Rosseland opacity scale to the new scale and from the previous
C     mass depth scale to the new one
C
      IF(IDEPTH.LE.2) THEN
         DO 85 I=1,NDEPTH
            TEMP0(I)=TEMP(I)
            ELEC0(I)=ELEC(I)
            DENS0(I)=DENS(I)
            TAU(I)=LOG(TAUROS(I))
            DEPTH0(I)=LOG(DM(I))
   85    CONTINUE
         CALL INTERP(TAU,DEPTH0,TAU0,DM0,NDEPTH,ND,3,0,0)
         CALL INTERP(TAU,TEMP0,TAU0,TEMP,NDEPTH,ND,3,0,0)
         CALL INTERP(TAU,ELEC0,TAU0,ELEC,NDEPTH,ND,3,0,1)
         CALL INTERP(TAU,DENS0,TAU0,DENS,NDEPTH,ND,3,0,1)
         DO 90 I=1,ND
            DM(I)=EXP(DM0(I))
            PTOTAL(I)=DM(I)*GRAV+PRAD0
            PGS(I)=(DENS(I)/WMM(I)+ELEC(I))*BOLK*TEMP(I)
   90    CONTINUE
      END IF
C
C     Recalculation of the populations
C
c     DO 100 ID=1,ND
c        CALL WNSTOR(ID)
c        CALL STEQEQ(ID,POP,1)
c 100 CONTINUE
  110 CONTINUE
      IF(HMIX0.GE.0.) THEN
         WRITE(6,600)
         CALL CONOUT(2,IPRING)
      END IF
      LCHC=LCHC0
      LTE=LTE0
      IRSPLT=IRSPL0
  600 FORMAT(1H1,' CONVECTIVE FLUX: AT THE END OF LTEGR'/)
  601 FORMAT(1H1, 'COMPUTED LTE-GREY MODEL'//'    ID    TAU',7X,
     * 'MASS',5X,'TEMP',7X,'N',10X,'NE',9X,'P',9X,'ROSS.OP'/)
  602 FORMAT(1H ,I4,1P2D11.3,0PF8.0,1P4D11.3)
      RETURN
      END
C
C
C ********************************************************************
C
C

      SUBROUTINE ROSSOP(ID,P,TAUR,HOPF,T4,T,ANE,ABROSS)
C     =================================================
C
C     Auxiliary procedure for LTEGR
C     Evaluation of temperature, electron density, and Rosseland
C     opacity for a given TAUR (Rosseland optical depth) and P (total
C     pressure)
C
C     Input parameters:
C     ID     - depth index
C     P      - total pressure
C     TAUR   - Rosseland optical depth
C     HOPF   - mode of evaluating Hopf function;
C            = 0  -  exact Hopf function
C            > 0  -  constant Hopf function to HOPF
C     T4     = effective temperature ** 4
C
C     Output:
C     T      - temperature
C     ANE    - electron density
C     ABROSS - Rosseland opacity (per gram)
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'ATOMIC.FOR'
      INCLUDE 'MODELQ.FOR'
      INCLUDE 'ALIPAR.FOR'
      DIMENSION A(5)
C
      DATA A/0.71044609D0,-0.2830385D0,0.57975839D0,-0.75751038D0,
     *       0.45026781D0/
      SAVE A
C
C     Hopf function
C
      X=HOPF
      IF(X.GT.0.) GO TO 20
      X=A(1)
      IF(TAUR.GT.160.) GO TO 20
      EX=EXP(-TAUR)
      E1=EXPINT(TAUR)
      E=E1
      DO 10 I=1,4
         E=(EX-TAUR*E)/I
   10    X=X+E*A(I+1)
C
C     Temperature
C
   20 T=(0.75*T4*(TAUR+X)+EXTOT)**0.25
C
C     Determination of electron density from the total pressure
C
      if(ioptab.ge.-1) then
         AN=P/T/BOLK
         CALL ELDENS(ID,T,AN,ANE,ENRG)
         RHO=WMM(ID)*(AN-ANE)
         DENS(ID)=RHO
C
C     temperature and electron density are transmitted to SABOLF
C     and RATMAT through arrays TEMP and ELEC
C
         TEMP(ID)=T
         ELEC(ID)=ANE
C
C     Corresponding LTE populations
C
         if(ioptab.ge.0) then
c           CALL WNSTOR(ID)
c           CALL STEQEQ(ID,POP,1)
C
C     Finally, evaluation of the Rosseland opacity for the new values
C     of temperature, electron density, and populations
C     (ROSS - Rosseland opacity per 1 cm**3)
C
            CALL OPACF0(ID,NFREQ)
            CALL MEANOP(T,ABSO,SCAT,OPROS,OPPLA)
            ABROSS=OPROS/RHO
            ABROSD(ID)=ABROSS
            ABPLAD(ID)=OPPLA/RHO
           else
            call meanopt(t,id,rho,opros,oppla)
            abrosd(id)=opros
            abplad(id)=oppla
            abross=opros
         end if
       else
         temp(id)=t
         rho=rhoeos(t,p)
         dens(id)=rho
         call meanopt(t,id,rho,opros,oppla)
         abrosd(id)=opros
         abplad(id)=oppla
         abross=opros
      end if
      RETURN
      END
C
C
C     ****************************************************************
C
C
 
      SUBROUTINE CONTMP
C     =================
C
C     Auxiliary procedure for LTEGR
C     Determination of temperature in convectively unstable layers
C     This is done by solving the energy balance equation
C     F(rad)+F(conv)=F(mech), which yields a cubic equation for
C     the logarithmic temperature gradient
C
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'ATOMIC.FOR'
      INCLUDE 'MODELQ.FOR'
      INCLUDE 'ALIPAR.FOR'
      COMMON ESEMAT(MLEVEL,MLEVEL),BESE(MLEVEL),
     *       DEPTH(MDEPTH),DEPTH0(MDEPTH),TAU(MDEPTH),TAU0(MDEPTH),
     *       TEMP0(MDEPTH),ELEC0(MDEPTH),DENS0(MDEPTH),DM0(MDEPTH)
      DIMENSION DELTR(MDEPTH),TEMPR(MDEPTH),ICON0(MDEPTH)
      COMMON/CUBCON/A,B,DEL,GRDADB,DELMDE,RHO,FLXTOT,GRAVD
      common/ichndm/ichanm
      PARAMETER  (ERRT=1.D-3)
C
C     First, store the temperature(rad) and gradient Delta(rad) -
C     quantities for the  purely raditive equilibrium model
C
      T4=TEFF**4
      FLXTO0=SIG4P*T4
      DPRAD=1.891204931D-15*T4
      if(ifprad.eq.0) dprad=0.
C
      PRAD0=DPRAD/1.732D0
      DO 10 ID=1,ND
         TEMPR(ID)=TEMP(ID)
         IF(ID.EQ.1) THEN
            DELTR(ID)=0.
          ELSE
            if(ilgder.eq.0) then
            DELTR(ID)=
     *      (TEMP(ID)-TEMP(ID-1))/(PTOTAL(ID)-PTOTAL(ID-1))*
     *      (PTOTAL(ID)+PTOTAL(ID-1))/(TEMP(ID)+TEMP(ID-1))
            else
            DELTR(ID)=
     *      LOG(TEMP(ID)/TEMP(ID-1))/LOG(PTOTAL(ID)/PTOTAL(ID-1))
            end if
         END IF
   10 CONTINUE
      ICONIT=0
C
C     ------------------------------------------------------
C     Global iteration loop for calculating convective model
C     ------------------------------------------------------
C
   20 ICONIT=ICONIT+1
      ICONBE=0
      CHANTM=0.
      DO 50 ID=1,ND
         T=TEMP(ID)
         PTOT=PTOTAL(ID)
         PGAS=PGS(ID)
         PTURB=HALF*DENS(ID)*VTURB(ID)*VTURB(ID)
         PRAD=PTOT-PGAS-PTURB
         PRADT(ID)=PRAD
         FLXTOT=FLXTO0
         IF(IDISK.EQ.1) THEN
c           FLXTOT=FLXTO0*(UN-THETAV(ID))
            FLXTOT=FLXTO0*(UN-THETA(ID))
            GRAVD=ZD(ID)*QGRAV
         END IF
         ICON0(ID)=0
C
         IF(ID.EQ.1) GO TO 40
         J=0
         IF(ICONIT.EQ.1) T=T-TEMPR(ID-1)+TEMP(ID-1)
         TM=TEMP(ID-1)
         IF(T.LT.0.) T=TM
         PTOTM=PTOTAL(ID-1)
         if(ilgder.eq.0) then
         PT0=HALF*(PTOT+PTOTM)
         else
         pt0=sqrt(ptot*ptotm)
         end if
         DELR=DELTR(ID)
C
C        Inner iteration loop for determining temperature in the
C        conectively unstable layers
C
   30    J=J+1
         TOLD=T
         if(ilgder.eq.0) then
         T0=HALF*(T+TM)
         PG0=HALF*(PGAS+PGM)
         PR0=HALF*(PRAD+PRADM)
         AB0=HALF*(ABROSD(ID)+ABROSD(ID-1))
         else
         t0=sqrt(t*tm)
         pg0=sqrt(pgas*pgm)
         pr0=sqrt(prad*pradm)
         ab0=sqrt(abrosd(id)*abrosd(id-1))
         end if
         IF(ID.GE.ND-2.AND.ICONBE.EQ.0) GO TO 40
         CALL CONVEC(ID,T0,PT0,PG0,PR0,AB0,DELR,FLXCNV,VCON)
         IF(FLXCNV.EQ.0..or.id.lt.idconz) GO TO 40
         ICON0(ID)=1
         ICONBE=1
         CALL CUBIC(DELTA0)
         REFF=DELTA0/DELR
         PRAD=PRADM+(TAUROS(ID)-TAUROS(ID-1))*DPRAD*REFF
         PRADT(ID)=PRAD
         PGAS=PTOT-PRAD-PTURB
         if(ilgder.eq.0) then
         IF(REFF.GT.UN) REFF=UN
         IF(REFF.LT.0.) REFF=0.
         FAC=DELTA0*(PTOT-PTOTM)/(PTOT+PTOTM)
         T=TM*(UN+FAC)/(UN-FAC)
         IF(T.LT.TM) T=TM
         else
         T=TM*(PTOT/PTOTM)**DELTA0
         IF(T.LT.TM) T=TM*1.0001
         end if
         if(iconit.eq.1) then
c           flxcnv=flco0*delmde**1.5
            t00=half*(t+tm)
            CALL CONVC1(ID,T00,PT0,PG0,PR0,AB0,DELR,FLXCNV,FC0)
            flxcn0=fc0*delmde**1.5
            write(6,613) id,j,t0,delta0,delr,flxcn0/flxto0,
     *                   flxcnv/flxto0,fac,temp(id),t
         end if
c  613    format(2i4,f8.1,1p4e13.5,0p2f8.1)
  613    format(2i3,f8.1,2f8.4,2f9.3,f10.6,2f8.1)
         IF(ABS(UN-T/TOLD).GT.ERRT.AND.J.LT.10) GO TO 30
C
C     Store the final quantitites
C
   40    IF(ID.GT.1.AND.ICON0(ID).EQ.0.AND.ICON0(ID-1).EQ.1)
     *      DELTC=DELT0
         DELT0=TEMP(ID)-T
         IF(TEMP(ID).NE.0.) CHANT0=ABS((T-TEMP(ID))/TEMP(ID))
         IF(CHANT0.GT.CHANTM) CHANTM=CHANT0
         TEMP(ID)=T
         IF(ICONIT.GT.1.AND.ICON0(ID).EQ.0.AND.ICONBE.EQ.1)
     *      TEMP(ID)=T-DELTC
         PGM=PGAS
         PRADM=PRAD
         PGS(ID)=PGAS
c         WRITE(110,603) ID,ICONIT,J,t,CHANT0,pradt(id)
   50 CONTINUE
c 603 format(3i4,f12.1,1p2e11.3)
C
C     Diagnostic outprint
C
      IF(IPRING.EQ.2) THEN
         WRITE(6,600) ICONIT
         CALL CONOUT(1,IPRING)
      END IF
  600 FORMAT(1H1,' CONVECTIVE FLUX: AT CONTMP, ITER=',I2/)
C
C     2. New values of electron density, density, sound spped,
C        and mean opacities and optical depths
C
c
      ANEREL=ELEC(1)/(DENS(1)/WMM(1)+ELEC(1))
      DO 70 ID=1,ND
         T=TEMP(ID)
         P=PTOTAL(ID)
         ITINT=0
   60    ITINT=ITINT+1
         if(ioptab.ge.-1) then
            AN=PGS(ID)/T/BOLK
            CALL ELDENS(ID,T,AN,ANE,ENRG)
            ELEC(ID)=ANE
            DENS(ID)=WMM(ID)*(AN-ANE)
            PHMOL(ID)=AHMOL
C
C     Corresponding LTE populations
C
            if(ioptab.ge.0) then
c              CALL WNSTOR(ID)
c              CALL STEQEQ(ID,POP,1)
C
C        Evaluation of the Rosseland and Planck mean opacities
C
               CALL OPACF0(ID,NFREQ)
               CALL MEANOP(T,ABSO,SCAT,OPROS,OPPLA)
               ABROS=OPROS/DENS(ID)
               ABPLA=OPPLA/DENS(ID)
             else
               rho=dens(id)
               call meanopt(t,id,rho,abros,abpla)
               abrosd(id)=abros
               abplad(id)=abpla
            end if
          else
            rho=rhoeos(t,p)
            dens(id)=rho
            call meanopt(t,id,rho,abros,abpla)
            abrosd(id)=abros
            abplad(id)=abpla
         end if
C
C        New values of the the column mass
C
         PTOLD=PTOTAL(ID)
         if(idisk.eq.0.and.ichanm.gt.0) then
         IF(ID.EQ.1) THEN
            DM(ID)=TAUROS(ID)/ABROS
            PTOTAL(ID)=DM(ID)*GRAV+PRAD0
          ELSE
            DM(ID)=DM(ID-1)+(TAUROS(ID)-TAUROS(ID-1))/
     *             HALF/(ABROSD(ID-1)+ABROS)
            PTOTAL(ID)=DM(ID)*GRAV+PRAD0
c            DELTR(ID)=
c     *      (TEMPR(ID)-TEMPR(ID-1))/(PTOTAL(ID)-PTOTAL(ID-1))*
c     *      (PTOTAL(ID)+PTOTAL(ID-1))/(TEMPR(ID)+TEMPR(ID-1))
         END IF
C
C        Store the final quantitites
C
         PTURB=HALF*DENS(ID)*VTURB(ID)*VTURB(ID)
         PGS(ID)=PTOTAL(ID)-PRADT(ID)-PTURB
         end if
         ABROSD(ID)=ABROS
         ABPLAD(ID)=ABPLA
         IF((PTOTAL(ID)-PTOLD)/PTOLD.LT.1.D-3) GO TO 70
         IF(ITINT.GT.5) THEN
            WRITE(6,601) ID,PTOLD,PTOTAL(ID)
            GO TO 70
          ELSE
            GO TO 60
         END IF
   70 CONTINUE
C *** TEMPORARY
C
C     IF(CHANTM.GT.ERRT.AND.ICONIT.LT.NCONIT) GO TO 20
      IF(ICONIT.LT.NCONIT) GO TO 20
  601 FORMAT(1H0,'SLOW CONVERGENCE OF INTERNAL ITERATIONS IN',
     *       ' CONTMP: ID, PTOT(OLD), PTOT(NEW) ='/I3,1P2D10.2/)
      RETURN
      END
C
C
C     ****************************************************************
C
C
 
      SUBROUTINE CONTMD
C     =================
C
C     Auxiliary procedure for LTEGRD
C     Determination of temperature in convectively unstable layers
C     for disks.
C     This is done by solving the energy balance equation
C     F(rad)+F(conv)=F(mech), which yields a cubic equation for
C     the logarithmic temperature gradient
C
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'ATOMIC.FOR'
      INCLUDE 'MODELQ.FOR'
      INCLUDE 'ALIPAR.FOR'
      COMMON ESEMAT(MLEVEL,MLEVEL),BESE(MLEVEL),
     *       DEPTH(MDEPTH),DEPTH0(MDEPTH),TAU(MDEPTH),TAU0(MDEPTH),
     *       TEMP0(MDEPTH),ELEC0(MDEPTH),DENS0(MDEPTH),DM0(MDEPTH)
      DIMENSION DELTR(MDEPTH),TEMPR(MDEPTH),ICON0(MDEPTH)
      COMMON/CUBCON/A,B,DEL,GRDADB,DELMDE,RHO,FLXTOT,GRAVD
      COMMON/PRSAUX/VSND2(MDEPTH),HG1,HR1,RR1
      PARAMETER  (ERRT=1.D-3)
C
C     First, store the temperature(rad) and gradient Delta(rad) -
C     quantities for the  purely raditive equilibrium model
C
      T4=TEFF**4
      FLXTO0=SIG4P*T4
      DPRAD=1.891204931D-15*T4
      if(ifprad.eq.0) dprad=0.
      PRAD0=DPRAD/1.732D0
      DO 10 ID=1,ND
         TEMPR(ID)=TEMP(ID)
         IF(ID.EQ.1) THEN
            DELTR(ID)=0.
          ELSE
            DELTR(ID)=
     *      (TEMP(ID)-TEMP(ID-1))/(PTOTAL(ID)-PTOTAL(ID-1))*
     *      (PTOTAL(ID)+PTOTAL(ID-1))/(TEMP(ID)+TEMP(ID-1))
         END IF
   10 CONTINUE
      ICONIT=0
C
C     ------------------------------------------------------
C     Global iteration loop for calculating convective model
C     ------------------------------------------------------
C
   20 ICONIT=ICONIT+1
      ICONBE=0
      HR1=FLXTO0*PCK*ABROSD(1)/QGRAV
      CHANTM=0.
      DO 50 ID=1,ND
         T=TEMP(ID)
         PTOT=PTOTAL(ID)
         PGAS=PGS(ID)
         PTURB=HALF*DENS(ID)*VTURB(ID)*VTURB(ID)
         PRAD=PRADT(ID)
c        FLXTOT=FLXTO0*(UN-THETAV(ID))
         FLXTOT=FLXTO0*(UN-THETA(ID))
         GRAVD=ZD(ID)*QGRAV
         ICON0(ID)=0
C
         IF(ID.EQ.1) GO TO 40
         J=0
         IF(ICONIT.EQ.1) T=T-TEMPR(ID-1)+TEMP(ID-1)
         TM=TEMP(ID-1)
         IF(T.LT.0.) T=TM
         PGM=PGS(ID-1)
         PTOTM=PTOTAL(ID-1)
         PT0=HALF*(PTOT+PTOTM)
         DELR=DELTR(ID)
C
C        Inner iteration loop for determining temperature in the
C        conectively unstable layers
C
   30    J=J+1
         TOLD=T
         T0=HALF*(T+TM)
         PG0=HALF*(PGAS+PGM)
         PR0=HALF*(PRAD+PRADM)
         AB0=HALF*(ABROSD(ID)+ABROSD(ID-1))
         IF(ID.GE.ND-2.AND.ICONBE.EQ.0) GO TO 40
         CALL CONVEC(ID,T0,PT0,PG0,PR0,AB0,DELR,FLXCNV,VCON)
         IF(FLXCNV.EQ.0.) GO TO 40
         ICON0(ID)=1
         ICONBE=1
         if(id.eq.nd) then
            pip=(ptot+ptotm)/(ptot-ptotm)
            t=tm*(pip+delr)/(pip-delr)
            go to 40
         end if
         CALL CUBIC(DELTA0)
         FAC=DELTA0*(PTOT-PTOTM)/(PTOT+PTOTM)
         T=TM*(UN+FAC)/(UN-FAC)
         IF(T.LT.TM) T=TM
         IF(ABS(UN-T/TOLD).GT.ERRT.AND.J.LT.10) GO TO 30
C
C     Store the final quantitites
C
   40    IF(ID.GT.1.AND.ICON0(ID).EQ.0.AND.ICON0(ID-1).EQ.1)
     *      DELTC=DELT0
         if(id.eq.nd) then
            pip=(ptot+ptotm)/(ptot-ptotm)
            t=tm*(pip+delr)/(pip-delr)
         end if
         DELT0=TEMP(ID)-T
         PRADT(ID)=PRADT(ID)*(T/TEMP(ID))**4
         DENS(ID)=DENS(ID)*(TEMP(ID)/T)
         IF(TEMP(ID).NE.0.) CHANT0=ABS((T-TEMP(ID))/TEMP(ID))
         IF(CHANT0.GT.CHANTM) CHANTM=CHANT0
         TEMP(ID)=T
         IF(ICONIT.GT.1.AND.ICON0(ID).EQ.0.AND.ICONBE.EQ.1)
     *      TEMP(ID)=T-DELTC
         PRADM=PRADT(ID)
c         WRITE(110,603) ID,ICONIT,J,t,CHANT0,pradt(id)
   50 CONTINUE
c 603 format(3i4,f12.1,1p2e11.3)
C
C     Diagnostic outprint
C
      IF(IPRING.EQ.2) THEN
         WRITE(6,600) ICONIT
         CALL CONOUT(1,IPRING)
      END IF
  600 FORMAT(1H1,' CONVECTIVE FLUX: AT CONTMD, ITER=',I2/)
C
C     2. New values of electron density and density
C
c     CALL HESOL6
C
C     Evaluation of the Rosseland and Planck mean opacities
C
      DO ID=1,ND
         T=TEMP(ID)
c        CALL WNSTOR(ID)
c        CALL STEQEQ(ID,POP,1)
         CALL OPACF0(ID,NFREQ)
         CALL MEANOP(T,ABSO,SCAT,OPROS,OPPLA)
         ABROS=OPROS/DENS(ID)
         ABPLA=OPPLA/DENS(ID)
         ABROSD(ID)=ABROS
         ABPLAD(ID)=ABPLA
      END DO
      IF(CHANTM.GT.ERRT.AND.ICONIT.LT.NCONIT) GO TO 20
      RETURN
      END
C
C
C     ****************************************************************
C
C
 
      SUBROUTINE ELDENS(ID,T,AN,ANE,ENRG)
C     ===================================
C
C     Evaluation of the electron density and the total hydrogen
C     number density for a given total particle number density
C     and temperature;
C     by solving the set of Saha equations, charge conservation and
C     particle conservation equations (by a Newton-Raphson method)
C
C     Input parameters:
C     T    - temperature
C     AN   - total particle number density
C
C     Output:
C     ANE   - electron density
C     ANP   - proton number density
C     AHTOT - total hydrogen number density
C     AHMOL - relativer number of hydrogen molecules with respect to the
C             total number of hydrogens
C     ENERG - part of the internal energy: excitation and ionization
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'MODELQ.FOR'
      INCLUDE 'ATOMIC.FOR'
      DIMENSION R(3,3),S(3),P(3)
      if(ioptab.lt.-1) return
C
      if(ifmol.gt.0) then
         aein=an*anerel
         call moleq(id,t,an,aein,ane,entt,enrg)
         if(ifentr.gt.0) enrg=entt
         anerel=ane/an
         return
      end if
c
      QM=0.
      Q2=0.
      QP=0.
      Q=0.
      DQN=0.
      TK=BOLK*T
      THET=5.0404D3/T
C
C     Coefficients entering ionization (dissociation) balance of:
C     atomic hydrogen          - QH;
C     negative hydrogen ion    - QM   (considered only if IHM>0);
C     hydrogen molecule        - QP   (considered only if IH2>0);
C     ion of hydrogen molecule - Q2   (considered only if IH2P>0).
C
      IF(IATREF.EQ.IATH) THEN
      IF(IHM.GT.0)  QM=1.0353D-16/T/SQRT(T)*EXP(8762.9/T)
      IF(IH2P.GT.0) QP=TK*EXP((-11.206998+THET*(2.7942767+THET*
     *             (0.079196803-0.024790744*THET)))*2.30258509299405)
      IF(IH2.GT.0) Q2=TK*EXP((-12.533505+THET*(4.9251644+THET*
     *            (-0.056191273+0.0032687661*THET)))*2.30258509299405)
      QH0=EXP((15.38287+1.5*LOG10(T)-13.595*THET)*2.30258509299405)*two
c     QH=EXP((15.38287+1.5*LOG10(T)-13.595*THET)*2.30258509299405)
      END IF
C
C     Initial estimate of the electron density
C
      if(anerel.le.0.) then
         if(t.gt.1.e4) then
            anerel=0.5
          else
            if(elec(id).gt.0..and.dens(id).gt.0.) then
               anerel=elec(id)/(elec(id)+dens(id)/wmm(id))
             else
               anerel=0.1
            end if
         end if
      end if
c
      ANE=AN*ANEREL
      IT=0
C
C     Basic Newton-Raphson loop - solution of the non-linear set
C     for the unknown vector P, consistiong of AH, ANH (neutral
C     hydrogen number density) and ANE.
C
   10 IT=IT+1
C
C     procedure STATE determines Q (and DQN) - the total charge (and its
C     derivative wrt temperature) due to ionization of all atoms which
C     are considered (both explicit and non-explicit), by solving the set
C     of Saha equations for the current values of T and ANE
C
      CALL STATE(1,ID,T,ANE)
C
C     Auxiliary parameters for evaluating the elements of matrix of
C     linearized equations.
C     Note that complexity of the matrix depends on whether the hydrogen
C     molecule is taken into account
C     Treatment of hydrogen ionization-dissociation is based on
C     Mihalas, in Methods in Comput. Phys. 7, p.10 (1967)
C
      IF(IATREF.EQ.IATH) THEN
      qh=qh0/pfhyd
      G2=QH/ANE
      G3=0.
      G4=0.
      G5=0.
      D=0.
      E=0.
      G3=QM*ANE
      A=UN+G2+G3
      D=G2-G3
      IF(IT.GT.1) GO TO 60
      IF(IH2.EQ.0.AND.IH2P.EQ.0) GO TO 40
      IF(IH2.EQ.0) GO TO 20
      E=G2*QP/Q2
      B=TWO*(UN+E)
      GG=ANE*Q2
      GO TO 30
   20 B=TWO
      E=UN
      GG=G2*ANE*QP
   30 C1=B*(GG*B+A*D)-E*A*A
      C2=A*(TWO*E+B*Q)-D*B
      C3=-E-B*Q
      F1=(SQRT(C2*C2-4.*C1*C3)-C2)*HALF/C1
      FE=F1*D+E*(UN-A*F1)/B+Q
      GO TO 50
   40 F1=UN/A
      FE=D/A+Q
   50 AH=ANE/FE
      ANH=AH*F1
   60 AE=ANH/ANE
      GG=AE*QP
      E=ANH*Q2
      B=ANH*QM
C
C     Matrix of the linearized system R, and the rhs vector S
C
      R(1,1)=YTOT(ID)
      R(1,2)=0.
      R(1,3)=UN
      R(2,1)=-Q
      R(2,2)=-D-TWO*GG
      R(2,3)=UN+B+AE*(G2+GG)-DQN*AH
      R(3,1)=-UN
      R(3,2)=A+4.*(E+GG)
      R(3,3)=B-AE*(G2+TWO*GG)
      S(1)=AN-ANE-YTOT(ID)*AH
      S(2)=ANH*(D+GG)+Q*AH-ANE
      S(3)=AH-ANH*(A+TWO*(E+GG))
C
C     Solution of the linearized equations for the correction vector P
C
      CALL LINEQS(R,S,P,3,3)
C
C     New values of AH, ANH, and ANE
C
      AH=AH+P(1)
      ANH=ANH+P(2)
      DELNE=P(3)
      ANE=ANE+DELNE
C
C     hydrogen is not the reference atom
C
      ELSE
C
C     Matrix of the linearized system R, and the rhs vector S
C
      IF(IT.EQ.1) THEN
         ANE=AN*HALF
         AH=ANE/YTOT(ID)
      END IF
      R(1,1)=YTOT(ID)
      R(1,2)=UN
      R(2,1)=-Q-QREF
      R(2,2)=UN-(DQN+DQNR)*AH
      S(1)=AN-ANE-YTOT(ID)*AH
      S(2)=(Q+QREF)*AH-ANE
C
C     Solution of the linearized equations for the correction vector P
C
      CALL LINEQS(R,S,P,2,3)
      AH=AH+P(1)
      DELNE=P(2)
      ANE=ANE+DELNE
      END IF
C
C     Convergence criterion
C
      IF(ANE.LE.0.) ANE=1.D-3*AN
      IF(ABS(DELNE/ANE).GT.1.D-3.AND.IT.LE.10) GO TO 10
C
C     ANEREL is the exact ratio betwen electron density and total
C     particle density, which is going to be used in the subseguent
C     call of ELDENS
C
      ANEREL=ANE/AN
      AHTOT=AH
      IF(IATREF.EQ.IATH) THEN
      AHMOL=TWO*ANH*(ANH*Q2+ANH/ANE*QP)/AH
      ANP=ANH/ANE*QH
      END IF
C
      RETURN
      END
C
C
C ***********************************************************************
C
C
 
      SUBROUTINE MEANOP(T,ABSO,SCAT,OPROS,OPPLA)
C     ==========================================
C
C     Rosseland and Planck mean opacities
C
C     Input parameters:
C      T    - temperature
C      ABSO - array of absorption coefficients in all explicit
C             frequency points
C      SCAT - array of scttering coefficients
C     Output:
C      OPROS - Rosseland opacity (per 1 cm**3)
C      OPPLA - Planck mean opacity (per 1 cm**3)
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'MODELQ.FOR'
      INCLUDE 'ATOMIC.FOR'
      DIMENSION ABSO(MFREQ),SCAT(MFREQ)
C        
      ABR=0.
      SUMDB=0.
      ABP=0.
      SUMB=0.
      HKT=HK/T
C          
      DO IJ=1,NFREQC
         FR=FREQ(IJ)
       X=HKT*FR
       IF(X.GT.150.) X=150.
         EX=EXP(X)
         E1=UN/(EX-UN)
         PLAN=BNUE(IJ)*E1*W(IJ)
         DPLAN=PLAN*HKT*FR*EX*E1
         ABR=ABR+DPLAN/ABSO(IJ)
         ABP=ABP+PLAN*(ABSO(IJ)-SCAT(IJ))
         SUMDB=SUMDB+DPLAN
         SUMB=SUMB+PLAN
      END DO
      OPROS=SUMDB/ABR
      OPPLA=ABP/SUMB
      RETURN
      END
 
C
C
C     ****************************************************************
C
C
      SUBROUTINE MEANOPT(T,ID,RHO,OPROS,OPPLA)
C     ========================================
C
C     Rosseland and Planck mean opacities
C
C     Input parameters:
C      T    - temperature
C      RHO  - density
C     Output:
C      OPROS - Rosseland opacity (per gram)
C      OPPLA - Planck mean opacity (per gram)
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'MODELQ.FOR'
C        
      ABR=0.
      SUMDB=0.
      ABP=0.
      SUMB=0.
      HKT=HK/T
C          
      DO 10 IJ=1,NFREQ
         FR=FREQ(IJ)
         EX=EXP(HKT*FR)
         E1=UN/(EX-UN)
         PLAN=BNUE(IJ)*E1*W(IJ)
         DPLAN=PLAN*HKT*FR*EX*E1
         CALL OPCTAB(FR,IJ,ID,T,RHO,AB,SC,SCT,1)
         ABR=ABR+DPLAN/(AB+SCT)
         ABP=ABP+PLAN*AB
         SUMDB=SUMDB+DPLAN
         SUMB=SUMB+PLAN
   10 CONTINUE
      OPROS=SUMDB/ABR
      OPPLA=ABP/SUMB
      RETURN
      END
 
C
C
C     ****************************************************************
C
C
 
      SUBROUTINE CONVEC(ID,T,PTOT,PG,PRAD,ABROS,DELTA,FLXCNV,VCONV)
C     =============================================================
C
C     Determination of the mixing-lengths convective flux
C
C     Input:  T     - temperature
C             PTOT  - total pressure
C             PG    - gas pressure
C             PRAD  - radiation pressure
C             ABROS - Rosseland opacity (per gram)
C             DELTA - corresponding temperature gradient
C     Output: FLXCNV - convective flux (expressed as H, ie F/4/pi)
C             VCONV  - convective velocity
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      COMMON/CUBCON/A,B,DDEL,GRDADB,DLT,RHO,FLXTOT,GRAVD
C
      VCONV=0.
      FLXCNV=0.
      DLT=0.
      IF(HMIX0.LT.0.) RETURN
C 
C     Thermodynamic derivatives
C
      if(ioptab.ge.-1) then
         CALL TRMDER(ID,T,PG,PRAD,TAURS(ID),HEATCP,DLRDLT,GRDADB,RHO)
       else
         call trmdrt(id,t,ptot,heatcp,dlrdlt,grdadb,rho)
      end if
      DDEL=DELTA-GRDADB
C
C     Convective instability criterion
C
      IF(DDEL.LT.0.) GO TO 100
      if(idisk.eq.0) then
         HSCALE=PTOT/RHO/GRAV
       else
         if(gravd.eq.0.) go to 100
         hscale=ptot/rho/gravd
      end if
      HMIX=HMIX0
      if(hmix0.eq.0.) hmix=1.
      VCO=HMIX*SQRT(ABS(aconml*PTOT/RHO*DLRDLT))
      FLCO=bconml*RHO*HEATCP*T*HMIX/12.5664
      TAUE=HMIX*ABROS*RHO*HSCALE
      FAC=TAUE/(UN+HALF *TAUE*TAUE)
C
C     Set up parameters A and B (see Mihalas, Eq. 7-76, 7-79, etc)
C
      B=5.67d-5*T**3/(rho*heatcp*VCO)*FAC*cconml*half
      IF(FLXTOT.GT.0.) A=FLCO*VCO/FLXTOT*DELTA
C
C     Determination of   Delta - Delta(E)
C
      D=B*B/2.D0
      DLT=D+DDEL-B*SQRT(D/2.D0+DDEL)
      IF(DLT.LT.0.) DLT=0.
C
C     Resulting convective velocity VCONV and flux FLXCNV
C
      VCONV=VCO*SQRT(DLT)
      FLXCNV=FLCO*VCONV*DLT
  100 CONTINUE
      RETURN
      END
C
C
C     ****************************************************************
C
C
 
      SUBROUTINE CONVC1(ID,T,PTOT,PG,PRAD,ABROS,DELTA,FLXCNV,FC0)
C     ===========================================================
C
C     Determination of the mixing-lengths convective flux
C
C     Input:  T     - temperature
C             PTOT  - total pressure
C             PG    - gas pressure
C             PRAD  - radiation pressure
C             ABROS - Rosseland opacity (per gram)
C             DELTA - corresponding temperature gradient
C     Output: FLXCNV - convective flux (expressed as H, ie F/4/pi)
C             VCONV  - convective velocity
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      COMMON/CUBCON/A,B,DDEL,GRDADB,DLT,RHO,FLXTOT,GRAVD
C
      VCONV=0.
      FLXCNV=0.
      DLT=0.
      IF(HMIX0.LT.0.) RETURN
C 
C     Thermodynamic derivatives
C
      if(ioptab.ge.-1) then
         CALL TRMDER(ID,T,PG,PRAD,TAURS(ID),HEATCP,DLRDLT,GRDADB,RHO)
       else
         call trmdrt(id,t,ptot,heatcp,dlrdlt,grdadb,rho)
      end if
      DDEL=DELTA-GRDADB
C
C     Convective instability criterion
C
      if(idisk.eq.0) then
         HSCALE=PTOT/RHO/GRAV
       else
         if(gravd.eq.0.) go to 100
         hscale=ptot/rho/gravd
      end if
      HMIX=HMIX0
      if(hmix0.eq.0.) hmix=1.
      VCO=HMIX*SQRT(ABS(aconml*PTOT/RHO*DLRDLT))
      FLCO=bconml*RHO*HEATCP*T*HMIX/12.5664
      FC0=FLCO*VCO
      IF(DDEL.LT.0.) GO TO 100
c
      TAUE=HMIX*ABROS*RHO*HSCALE
      FAC=TAUE/(UN+HALF *TAUE*TAUE)
C
C     Set up parameters A and B (see Mihalas, Eq. 7-76, 7-79, etc)
C
      B=5.67d-5*T**3/(rho*heatcp*VCO)*FAC*cconml*half
      IF(FLXTOT.GT.0.) A=FLCO*VCO/FLXTOT*DELTA
C
C     Determination of   Delta - Delta(E)
C
      D=B*B/2.D0
      DLT=D+DDEL-B*SQRT(D/2.D0+DDEL)
      IF(DLT.LT.0.) DLT=0.
C
C     Resulting convective velocity VCONV and flux FLXCNV
C
      VCONV=VCO*SQRT(DLT)
      FLXCNV=FLCO*VCONV*DLT
  100 CONTINUE
      RETURN
      END
C
C
C     ****************************************************************
C
C
 
      SUBROUTINE TRMDER(ID,T,PG,PRAD,TAU,HEATCP,DLRDLT,GRDADB,RHO)
C     ============================================================
C
C     Thermodynamic derivatives
C     Evaluation similar as in Kurucz's ATLAS
C
C     Input:  T     -  temperature
C             PG    -  gas pressure
C             PRAD  -  radiation pressure
C     Output: DEDT  -  d(energy)/d(T)
C             DRDT  -  d(rho)/d(T)
C             DEDPG -  d(energy)/d(PG)
C             DRDPG -  d(rho)/d(PG)
C             RHO   -  density
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      DIMENSION T1(3),P1(3),RHON(3),ENERGY(3),entrp(3)
C
C   NUMERICAL EVALUATION OF THERMODYNAMIC DERIVATIVES
C
      T1(1)=T*1.001D0
      T1(2)=T
      T1(3)=T
      P1(1)=PG
      P1(2)=PG
      P1(3)=PG*1.001D0
      DO I=1,3
         TT=T1(I)
         TKN=TT*1.38054D-16
         ANT=P1(I)/TKN
         CALL ELDENS(ID,TT,ANT,ANE,ENRG)
         RHON(I)=WMM(ID)*(ANT-ANE)
         ENERGY(I)=(1.5D0*P1(I)+ENRG+3.D0*PRAD*(TT/T)**4)/RHON(I)
         entrp(i)=enrg
      END DO
c
      DRDT=(RHON(1)-RHON(2))/T*1.D3
      DRDPG=(RHON(3)-RHON(2))/PG*1.D3
      RHO=RHON(2)
      DPDPG=1.
      dpdt=0.
      if(tau.lt.50.) DPDT=4.*PRAD/T*(un-exp(-tau))
      DLRDLT=T/RHO*(DRDT-DRDPG*DPDT/DPDPG)
      ptot=pg+prad
c
      if(ifentr.le.0) then
         DEDT=(ENERGY(1)-ENERGY(2))/T*1.D3
         DEDPG=(ENERGY(3)-ENERGY(2))/PG*1.D3
         HEATCV=DEDT-DEDPG*DRDT/DRDPG
         HEATCP=DEDT-DEDPG*DPDT/DPDPG-
     *          PTOT/RHO/RHO*(DRDT-DRDPG*DPDT/DPDPG)
         GRDADB=-PTOT/RHO/T*DLRDLT/HEATCP
c
       else
         dsdt=(entrp(1)-entrp(2))/t*1.e3
         dsdp=(entrp(3)-entrp(2))/pg*1.e3
         grdadb=-dsdp/dsdt*pg/t
         heatcp=-PTOT/RHO/T*DLRDLT/grdadb
      end if
c
      RETURN
      END
C
C
C     ***************************************************************
C
C
 
      SUBROUTINE CUBIC(DELTA)
C     =======================
C
C     Solution of the cubic equation for determination of
C     the true gradient DELTA
C
C     Input:   A,B   - coefficients; transmitted by COMMON/CUBCON
C              DEL   - DELTA(RAD) - DELTA(ADIAB); also transm. by CUBCON
C     Output:  DELTA - true gradient
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      COMMON/CUBCON/A,B,DEL,GRDADB,DELMDE,RHO,FLXTOT,GRAVD
      PARAMETER (THIRD = 0.333333333333333D0)
C
C     first solve the cubic equation
C         A*X**3 + X**2 + B*X = DEL
C     where X = (DELTA - DELTA(ELEM))**(1/2)
C
      AA=THIRD/A
      BB=B/A
      CC=-DEL/A
      P=BB*THIRD-AA*AA
      Q=AA**3-(BB*AA-CC)/2.D0
      D=Q*Q+P*P*P
      IF(D.GT.0.) THEN
         D=SQRT(D)
c        IF(D.EQ.ABS(Q)) THEN
         if(d-abs(q).lt.1.e-14*d) then
            SOL=(2.D0*D)**THIRD-AA
          ELSE
            D1=ABS(D-Q)
            D2=ABS(D+Q)
            SOL=D1/(D-Q)*D1**THIRD-D2/(D+Q)*D2**THIRD-AA
          END IF
       ELSE
         COSF=-Q/SQRT(ABS(P*P*P))
         TANF=SQRT(UN-COSF*COSF)/COSF
         FI=ATAN(TANF)*THIRD
         SOL=2.D0*SQRT(ABS(P))*COS(FI)-AA
      END IF
C
C     if the previous formalism gives an unphysical solution
C     x > DEL, then find the physical solution in the range (0, DEL)
C     by a Newton-Raphson solution of the cubic equation
C
      DELDA=SOL*(B+SOL)
      IF(DELDA.GT.DEL.OR.DELDA.LT.0.) THEN
         X0=sol
         J=0
   10    DELX=(DEL-X0*(B+X0+A*X0*X0))/(3.D0*A*X0*X0+2.D0*X0+B)
         X0=X0+DELX
         J=J+1
         IF(ABS(DELX/X0).GT.1D-6.AND.J.LT.50) GO TO 10
         SOL=X0
      END IF
C
C     finally, the actual gradient delta
C
      DELTA=GRDADB+B*SOL+SOL*SOL
      RETURN
      END
C
C
C     ***************************************************************
C
C
 
      SUBROUTINE CONOUT(IMOD,IPRIN)
C     =============================
C
C     Diagnostic outprint of temperature gradients, convective flux,
C     and their derivatives
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'MODELQ.FOR'
      INCLUDE 'ALIPAR.FOR'
      COMMON/CUBCON/A,B,DEL,GRDADB,DELMDE,RHO,FLXTOT,GRAVD
C
      IF(IPRIN.GT.0) WRITE(6,600)
      ANEREL=ELEC(1)/(DENS(1)/WMM(1)+ELEC(1))
      ICBEG=0
      FLXTO0=SIG4P*TEFF**4
      DO ID=1,ND
         T=TEMP(ID)
         PTOT=PTOTAL(ID)
         PG=PGS(ID)
         PRAD=PTOT-PG-HALF*DENS(ID)*VTURB(ID)**2
         IF(IMOD.EQ.2) THEN
            if(ioptab.ge.0) then
            CALL OPACF0(ID,NFREQ)
            CALL MEANOP(T,ABSO,SCAT,OPROS,OPPLA)
            ABROSD(ID)=OPROS/DENS(ID)
            else
            call meanopt(t,id,dens(id),opros,oppla)
            if(hmix0.gt.0.) abrosd(id)=opros
            end if
         END IF
         FLXTOT=flxto0
         if(idisk.eq.1) then
c           flxtot=flxto0*(un-thetav(id))
            flxtot=flxto0*(un-theta(id))
            gravd=zd(id)*qgrav
            prad=pradt(id)
         end if
         IF(ID.EQ.1) THEN
            TAU=DM(ID)*ABROSD(ID)
            FLXCR=0.
            GRDADB=0.
            DELTA(ID)=0.
            FLXC(ID)=0.
            FLXCNV=0.
          ELSE
            TM=TEMP(ID-1)
            TAU=TAUM+HALF*(DM(ID)-DM(ID-1))*(ABROSD(ID)+ABROSD(ID-1))
            PTOTM=PTOTAL(ID-1)
            PGM=PGS(ID-1)
            PRADM=PTOTM-PGM-HALF*DENS(ID-1)*VTURB(ID-1)**2
            if(idisk.eq.1) pradm=pradt(id-1)
            if(ilgder.eq.0) then
               T0=HALF*(T+TM)
               PT0=HALF*(PTOT+PTOTM)
               PG0=HALF*(PG+PGM)
               PR0=HALF*(PRAD+PRADM)
               AB0=HALF*(ABROSD(ID)+ABROSD(ID-1))
               DLT=(T-TM)/(PTOT-PTOTM)*PT0/T0
             else
               T0=SQRT(T*TM)
               PT0=SQRT(PTOT*PTOTM)
               PG0=SQRT(PG*PGM)
               PR0=SQRT(PRAD*PRADM)
               AB0=SQRT(ABROSD(ID)*ABROSD(ID-1))
               DLT=LOG(T/TM)/LOG(PTOT/PTOTM)  
            end if       
c           if(hmix0.gt.0.) DELTA(ID)=DLT
            DELTA(ID)=DLT
            flxcnv=0.
            if(idisk.ne.1.or.id.lt.nd)
     *      CALL CONVEC(ID,T0,PT0,PG0,PR0,AB0,DLT,FLXCNV,VCON)
            if(hmix0.gt.0.) FLXC(ID)=FLXCNV
            IF(ICBEG.EQ.0.AND.FLXC(ID).GT.0..AND.FLXC(ID-1).EQ.0..
     *         AND.ID.GT.25) ICBEG=ID
            if(icbeg.gt.0.and.flxc(id).gt.0.) icend=id
         END IF
       PRADR=PRAD/PTOT
         conrel=0.
         radrel=1.
         if(flxtot.gt.0.) then
           conrel=FLXCNV/FLXTOT
           radrel=flrd(id)/flxtot
         end if
c         IF(IPRIN.GT.0) WRITE(6,601) ID,TAU,T,PTOT,PRADR,
c     *                   PGS(ID),DELTA(ID),GRDADB,conrel
         IF(IPRIN.GT.0) WRITE(6,601) ID,TAU,T,DELTA(ID),
     *                  GRDADB,conrel,radrel,conrel+radrel

         TAUM=TAU
      END DO
      if(iprin.gt.0) write(6,603) icbeg,icend
  603 format(/' convective zone between depths (inclusive)  ',2i4/)
C
c     if ICONV=3, then:
c     in the convective zone the radiative (+ convective) equilibrium
c     is taken obligatorily in the differential form, i.e. 
c     NDRE is modified to have the value just below the beginning of
c     convection zone, and
c     rediff(id) has to be set to unity, and reint(id) to 0 for id => ndre
c
      IF(ICBEG.GT.3.AND.ICONV.EQ.3) THEN 
         NDRE=ICBEG-1
         DO ID=1,ND
            IF(ID.GE.NDRE) THEN
               REINT(ID)=0.
               REDIF(ID)=1.
             else
               reint(id)=1.
               redif(id)=0.
            END IF
         END DO
         WRITE(6,602) NDRE
      END IF
c
      IF(ICBEG.GT.3.AND.ICONV.EQ.2) THEN 
         NDRE=ICBEG-1
         DO ID=1,ND
            IF(ID.GE.NDRE) THEN
               REDIF(ID)=1.
            END IF
         END DO
         WRITE(6,602) NDRE
      END IF
c
c  600 FORMAT(1H0,'  ID',4X,'TAUR',5X,'TEMP',5X,'PTOT',3X,'PR/PTOT',
c     *       3X,'HYD.MOL',4X,'DELTA',2X,'DELTA(AD)',2X,'CON/TOT'//)
c  601 FORMAT(1H ,I4,1PD9.2,0PF9.1,1P3D9.2,2X,2D9.2,1X,D9.2)
  600 FORMAT(//'  ID',3X,'TAUR',7X,'TEMP',6X,
     *       'DELTA',2X,'DELTA(AD)',2X,'CON/TOT  RAD/TOT  (C+R)/TOT'//)
  601 FORMAT(1H ,I4,1PD9.2,0PF9.1,1P5D10.2)
  602 FORMAT(//' NDRE IS RESET IN CONOUT DUE TO THE EXISTENCE OF'
     *       ,' CONVECTIVE ZONE'/'  NDRE= ',I3/)
      RETURN
      END
C
C
C     ***************************************************************
C
C
 
      FUNCTION ERFCX(X)
C     ================
C
C     complementary error function
C     expression from Abramowitz and Stegun, p.299, Eq. 7.1.26
C
      INCLUDE 'IMPLIC.FOR'
      PARAMETER (P  = 0.3275911D0,
     *           A1 = 0.254829592D0,
     *           A2 =-0.284496736D0,
     *           A3 = 1.421413741D0,
     *           A4 =-1.453152027D0,
     *           A5 = 1.061405429D0,
     *           UN = 1.D0)
      T=UN/(UN+P*X)
      ERFCX=0.
      IF(X.GT.13.) RETURN
      ERFCX=T*(A1+T*(A2+T*(A3+T*(A4+T*A5))))*EXP(-X*X)
      RETURN
      END
 
C
C
C     ******************************************************************
C
C
C
      SUBROUTINE LEMINI
C     =================
C
C     Initializes necessary arrays for evaluating hydrogen line profiles
C     from the Lemke tables
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'MODELQ.FOR'
C      
      DO I=1,4
         DO J=1,22
            ILINH(I,J)=0
         END DO
      END DO
      if(ihydpr.eq.21) then
         open(unit=ihydpr,file='lemke.dat',status='old')
         write(6,641) ihydpr
       else if(ihydpr.eq.22) then
         open(unit=ihydpr,file='tremblay.dat',status='old')
         write(6,642) ihydpr
      end if
  641 format(' -----------'/
     *       ' reading Lemke tables; ihydpr =',i3,/
     *       ' -----------')
  642 format(' -----------'/
     *       ' reading Tremblay tables; ihydpr =',i3,/
     *       ' -----------')
C
C ---------------------------------
C     read Lemke or Tremblay tables
C ---------------------------------
C
      ILINE=0
      READ(IHYDPR,*) NTAB
      DO ITAB=1,NTAB
      ILINEB=ILINE
      READ(IHYDPR,*) NLLY
      DO ILI=1,NLLY
         ILINE=ILINE+1
         READ(IHYDPR,*) I,J,ALMIN,ANEMIN,TMIN,DLA,DLE,DLT,
     *                  NWL,NE,NT
         write(6,643) ntab,nlly,iline,i,j
  643    format(' ntab,nlly,iline,i,j ',5i4)
         ILINH(I,J)=ILINE
         NWLH(ILINE)=NWL
         NWLHYD(ILINE)=NWL
         NTH(ILINE)=NT
         NEH(ILINE)=NE
         DO IWL=1,NWL
            WLH(IWL,ILINE)=ALMIN+(IWL-1)*DLA
            WLHYD(ILINE,IWL)=WLH(IWL,ILINE)
            WLH(IWL,ILINE)=EXP(2.3025851*WLH(IWL,ILINE))
         END DO
         DO INE=1,NE
            XNELEM(INE,ILINE)=ANEMIN+(INE-1)*DLE
         END DO
         DO IT=1,NT
            XTLEM(IT,ILINE)=TMIN+(IT-1)*DLT
         END DO 
      END DO
c
      DO ILI=1,NLLY         
         ILNE=ILINEB+ILI
         NWL=NWLH(ILNE)
         READ(IHYDPR,500)
         DO INE=1,NEH(ILNE)
            DO IT=1,NTH(ILNE)
               READ(IHYDPR,*) QLT,(PRFHYD(ILNE,IWL,IT,INE),IWL=1,NWL)
            END DO
         END DO
C
C        coefficient for the asymptotic profile is determined from
C        the input data
C
         XCLOG=PRFHYD(ILNE,NWL,1,1)+2.5*WLHYD(ILNE,NWL)-0.477121
         XKLOG=0.6666667*XCLOG
         XK0(ILNE)=EXP(XKLOG*2.3025851)
      END DO
      END DO
      CLOSE(IHYDPR)
  500 FORMAT(1X)
C
      RETURN
      END
C
C
C ********************************************************************
C
C
      SUBROUTINE INTLEM(PRFH,WL0,ILINE,ID)
C     ====================================
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'MODELQ.FOR'
      PARAMETER (FOC1=1.25e-9,TTW=2./3.,VTBC=6.06e-9)
      DIMENSION PRFH(MHWL)
C
C     temperature is modified in order to account for the
C     effect of turbulent velocity on the Doppler width
C
      T=TEMP(ID)+VTBC*VTURBS(ID)*VTURBS(ID)
      ANE=ELEC(ID)
      TL=LOG10(T)
      ANEL=LOG10(ANE)
      F00=FOC1*EXP(TTW*LOG(ANE))
      XK=XK0(ILINE)
      FXK=F00*XK
      DOP=1.E8/WL0*SQRT(1.65E8*T)
      DBETA=WL0*WL0/2.997925E18/FXK
      BETAD=DBETA*DOP
C
C     interpolation to the actual values of temperature and electron
C     density
C
      NWL=NWLHYD(ILINE)
      DO IWL=1,NWL
         CALL INTHYD(PRFH0,TL,ANEL,IWL,ILINE)
         PRFH(IWL)=PRFH0
      END DO
      RETURN
      END
C
C
C ********************************************************************
C
C
      SUBROUTINE INTHYD(W0,X0,Z0,IWL,ILINE)
C     =====================================
C
C     Interpolation in temperature and electron density from the
C     Lemke tables for hydrogen lines to the actual valus of
C     temperature and electron density
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'MODELQ.FOR'
      DIMENSION ZZ(3),XX(3),WX(3),WZ(3)
C
      NX=2
      NZ=2
      NT=NTH(ILINE)
      NE=NEH(ILINE)
      BETA=WLH(IWL,ILINE)/XK
      IZH=1
C
C     for values lower than the lowest grid value of electron density
C     the profiles are determined by the approximate expression
C     (see STARKA); not by an extrapolation in the HYD tables which may
C     be very inaccurate
C
c      IF(Z0.LT.XNELEM(1,ILINE)*0.99.OR.Z0.GT.XNELEM(NE,ILINE)*1.01) THEN
      IF(Z0.LT.XNELEM(1,ILINE)*0.99) THEN
         CALL DIVSTR(IZH)
         W0=STARKA(BETA,TWO)*DBETA
         W0=LOG10(W0)
         GO TO 500
      END IF
C
C     Otherwise, one interpolates (or extrapolates for higher than the
C     highes grid value of electron density) in the HYD tables
C
      DO 10 IZZ=1,NE-1
         IPZ=IZZ
         IF(Z0.LE.XNELEM(IZZ+1,ILINE)) GO TO 20
   10 CONTINUE
   20 N0Z=IPZ-NZ/2+1
      IF(N0Z.LT.1) N0Z=1
      IF(N0Z.GT.NE-NZ+1) N0Z=NE-NZ+1
      N1Z=N0Z+NZ-1
C
      DO 300 IZZ=N0Z,N1Z
         I0Z=IZZ-N0Z+1
         ZZ(I0Z)=XNELEM(IZZ,ILINE)
C
C     Likewise, the approximate expression instead of extrapolation
C     is used for higher that the highest grid value of temperature,
C     if the Doppler width expressed in beta units (BETAD) is
C     sufficiently large (> 10)
C
         IF(X0.GT.1.01*XTLEM(NT,ILINE).AND.BETAD.GT.10.) THEN
            CALL DIVSTR(IZH)
            W0=STARKA(BETA,TWO)*DBETA
            W0=LOG10(W0)
            GO TO 500
         END IF
C
C     Otherwise, normal inter- or extrapolation
C
C     Both interpolations (in T as well as in electron density) are
C     by default the quadratic interpolations in logarithms
C
         DO 30 IX=1,NT-1
            IPX=IX
            IF(X0.LE.XTLEM(IX+1,ILINE)) GO TO 40
   30    CONTINUE
   40    N0X=IPX-NX/2+1
         IF(N0X.LT.1) N0X=1
         IF(N0X.GT.NT-NX+1) N0X=NT-NX+1
         N1X=N0X+NX-1
         DO 200 IX=N0X,N1X
            I0=IX-N0X+1
            XX(I0)=XTLEM(IX,ILINE)
            WX(I0)=PRFHYD(ILINE,IWL,IX,IZZ)
  200    CONTINUE
         IF(WX(1).LT.-99..OR.WX(2).LT.-99..OR.WX(3).LT.-99.) THEN
            CALL DIVSTR(IZH)
            W0=STARKA(BETA,TWO)*DBETA
            W0=LOG10(W0)
            GO TO 500
          ELSE
            WZ(I0Z)=YINT(XX,WX,X0)
         END IF
  300 CONTINUE
      W0=YINT(ZZ,WZ,Z0)
  500 CONTINUE
      RETURN
      END
C
C
C ********************************************************************
C
C
      FUNCTION YINT(XL,YL,XL0)
C     ========================
C
C     Quadratic interpolation routine
C
C     Input:  XL - array of x
C             YL - array of f(x)
C             XL0 - the point x(0) to which one interpolates
C
      INCLUDE 'IMPLIC.FOR'
      DIMENSION XL(3),YL(3)
      A0=(XL(2)-XL(1))*(XL(3)-XL(2))*(XL(3)-XL(1))
      A1=(XL0-XL(2))*(XL0-XL(3))*(XL(3)-XL(2))
      A2=(XL0-XL(1))*(XL(3)-XL0)*(XL(3)-XL(1))
      A3=(XL0-XL(1))*(XL0-XL(2))*(XL(2)-XL(1))
      YINT=(YL(1)*A1+YL(2)*A2+YL(3)*A3)/A0
      RETURN
      END
C
C
C     ******************************************************************
C
C 
      SUBROUTINE STARK0(I,J,IZZ,XKIJ,WL0,FIJ)
C     =======================================
C
C     Auxiliary procedure for evaluating the approximate Stark profile
C     of hydrogen lines - sets up necessary frequency independent
C     parameters
C
C     Input:  I     - principal quantum number of the lower level
C             J     - principal quantum number of the upper level
C             IZZ   - ionic charge (IZZ=1 for hydrogen, etc.)
C     Output: XKIJ  - coefficients K(i,j) for the Hotzmark profile;
C                     exact up to j=6, asymptotic for higher j
C             WL0   - wavelength of the line i-j
C             FIJ   - Stark f-value for the line i-j
C
      INCLUDE 'IMPLIC.FOR'
      PARAMETER (RYD1=911.763811,RYD2=911.495745/4.,CXKIJ=5.5E-5)
      PARAMETER (WI1=911.753578, WI2=227.837832)
      PARAMETER (UN=1.,TEN=10.,TWEN=20.,HUND=100.)
      DIMENSION FSTARK(10,4),XKIJT(5,4)
      DATA XKIJT/3.56D-4,5.23D-4,1.09D-3,1.49D-3,2.25D-3,.0125,.0177,
     * .028,.0348,.0493,.124,.171,.223,.261,.342,.683,.866,1.02,1.19,
     * 1.46/
      DATA FSTARK/.1387,.0791,.02126,.01394,.00642,4.814D-3,2.779D-3,
     * 2.216D-3,1.443D-3,1.201D-3,.3921,.1193,.03766,.02209,.01139,
     * 8.036D-3,5.007D-3,3.85D-3,2.658D-3,2.151D-3,.6103,.1506,.04931,
     * .02768,.01485,.01023,6.588D-3,4.996D-3,3.524D-3,2.838D-3,.8163,
     * .1788,.05985,.03189,.01762,.01196,7.825D-3,5.882D-3,4.233D-3,
     * 3.375D-3/
      SAVE XKIJT,FSTARK
C
      II=I*I
      JJ=J*J
      JMIN=J-I
      IF(JMIN.LE.5) THEN
         XKIJ=XKIJT(JMIN,I)
      ELSE
         XKIJ=CXKIJ*(II*JJ)*(II*JJ)/(JJ-II)
      END IF
      IF(JMIN.LE.10) THEN
         FIJ=FSTARK(JMIN,I)
      ELSE
         CFIJ=((TWEN*I+HUND)*J/(I+TEN)/(JJ-II))
         FIJ=FSTARK(10,I)*CFIJ*CFIJ*CFIJ
      END IF

      WL0=WI1
      IF(IZZ.EQ.2) WL0=WI2
      WL0=WL0/(UN/II-UN/JJ)
      RETURN
      END
C
C
C     ****************************************************************
C
C
 
      FUNCTION STARKA(BETA,FAC)
C     =========================
C
C     Approximate expressions for the hydrogen Stark profile
C
C     Input: BETA  - delta lambda in beta units,
C            BETAD - Doppler width in beta units
C            A     - auxiliary parameter
C                    A=1.5*LOG(BETAD)-1.761
C            DIV   - only for A > 1; division point between Doppler
C                    and asymptotic Stark wing, expressed in units
C                    of betad.
C                    DIV = solution of equation
C                    exp(-(beta/betad)**2)/betad/sqrt(pi)=
C                     = 1.5*FAC*beta**-5/2
C                    (ie. the point where Doppler profile is equal to
C                     the asymptotic Holtsmark)
C                    In order to save computer time, the division point
C                    DIV is calculated in advance by routine DIVSTR.
C            FAC   - Multiplicative factor (2. for H I; 1. for He II)
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'MODELQ.FOR'
      PARAMETER (F0=-0.5758228,F1=0.4796232,F2=0.07209481,AL=1.26)
      PARAMETER (SD=0.5641895,SLO=-2.5,THRA=1.5,BL1=1.14,BL2=11.4)
      PARAMETER (SAC=0.08, PISQ1=UN/1.77245385090551D0)
C
C     for a > 1 Doppler core + asymptotic Holtzmark wing with division
C               point DIV
C
      BETAD1=UN/BETAD
      IF(ADH.GT.AL) THEN
         XD=BETA*BETAD1
         IF(XD.LE.DIVH) THEN
            STARKA=EXP(-XD*XD)*BETAD1*PISQ1
          ELSE
            STARKA=THRA*FAC*EXP(SLO*LOG(BETA))
         END IF
       ELSE
C
C     empirical formula for a < 1
C
         IF(BETA.LE.BL1) THEN
            STARKA=SAC
          ELSE IF(BETA.LT.BL2) THEN
            XL=LOG(BETA)
            FL=(F0*XL+F1)*XL
            STARKA=F2*EXP(FL)
          ELSE
            STARKA=THRA*FAC*EXP(SLO*LOG(BETA))
         END IF
      END IF
      RETURN
      END
C
C
C *******************************************************************
C
C
 
      SUBROUTINE DIVSTR(IAH)
C     ======================
C
C     Auxiliary procedure for STARKA - determination of the division
C     point between Doppler and asymptotic Stark profiles
C
C     Input:  BETAD - Doppler width in beta units
C     Output: A     - auxiliary parameter
C                     A=1.5*LOG(BETAD)-1.671
C             DIV   - only for A > 1; division point between Doppler
C                     and asymptotic Stark wing, expressed in units
C                     of betad.
C                     DIV = solution of equation
C                     exp(-(beta/betad)**2)/betad/sqrt(pi)=3*beta**-5/2
C
C     He II: different definition of parameter ADH !
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'MODELQ.FOR'
      PARAMETER (UNQ=1.25,UNH=1.5,TWH=2.5,FO=4.,FI=5.)
      PARAMETER (CA=1.671,BL=5.821,AL=1.26,CX=0.28,DX=0.0001)
      PARAMETER (CA2=0.978,XA2=0.69314718)
C
      ADH=UNH*LOG(BETAD)-CA
      IF(IAH.EQ.2) ADH=ADH+XA2
      IF(BETAD.LT.BL) RETURN
      IF(ADH.GE.AL) THEN
         X=SQRT(ADH)*(UN+UNQ*LOG(ADH)/(FO*ADH-FI))
      ELSE
         X=SQRT(CX+ADH)
      ENDIF
      DO 10 I=1,5
         X2=X*X
         XN=X*(UN-(X2-TWH*LOG(X)-ADH)/(TWO*X2-TWH))
         IF(ABS(XN-X).LE.DX) GO TO 20
         X=XN
   10 CONTINUE
   20 DIVH=X
      RETURN
      END
C
C
C ********************************************************************
C
C
      SUBROUTINE OPAHST
C     =================
C
C     Auxiliary routine for START
C     sets up necessary parameters for routines OPAHYL and OPHYL1, i.e.
C     for opacity and emissivity in higher hydrogen lines
C     Also sets up Stark parameters for OPAHYL
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'ODFPAR.FOR'
C
      ALLIM1=1450.
      ABLIM1=6650.
      ABLIM2=5000.
      ABLIM3=6500.
C
C     Lyman lines
C
      ILOW=1
      IF(IABS(IOPHL1).EQ.1) IOPHL2=IOPHL2*2
      DO I=1,4
         M1FILE(I,ILOW)=MAX(I,IABS(IOPHL1))
         M2FILE(I,ILOW)=I+1
      END DO
      DO I=5,NLMX
         M1FILE(I,ILOW)=MAX(I-1,IABS(IOPHL1))
         M2FILE(I,ILOW)=MIN(I+3,NLMX)
      END DO
      M1FILE(NLMX,ILOW)=NLMX+1
      M2FILE(NLMX,ILOW)=NLMX
C
      IF(IABS(IOPHL1).GT.100) THEN
         IOPHL1=MOD(IOPHL1,100)
         ISET=0
   40    CONTINUE
         READ(IBUFF,*,ERR=90) IL1,IU1,IM1,IP1
         ISET=ISET+1
         IF(IL1.LE.0.AND.ISET.EQ.1) THEN
            IL1=1
            IU1=4
            IM1=0
            IP1=1
         END IF
         IF(IL1.LE.0.AND.ISET.EQ.2) THEN
            IL1=5
            IU1=NLMX
            IM1=1
            IP1=3
         END IF
         IUP1=MIN(IU1,NLMX)
         DO I=IL1,IUP1
            M1FILE(I,ILOW)=MAX(I-IM1,IABS(IOPHL1))
            M2FILE(I,ILOW)=MIN(I+IP1,NLMX)
         END DO
         IF(IU1.LT.NLMX) GO TO 40
   90 CONTINUE   
      READ(IBUFF,*,ERR=100) ALLIM1
         IF(ALLIM1.LE.0) ALLIM1=1450.
      END IF
      M1FILE(NLMX,ILOW)=NLMX+1
      M2FILE(NLMX,ILOW)=NLMX
C
C     Balmer lines
C
  100 ILOW=2
      IF(IABS(IOPHL2).EQ.1) IOPHL2=IOPHL2*3
      IF(IABS(IOPHL2).EQ.2) IOPHL2=IOPHL2*3/2
      DO I=1,6
         M1FILE(I,ILOW)=MAX(I,IABS(IOPHL2))
         M2FILE(I,ILOW)=I+1
      END DO
      DO I=7,NLMX
         M1FILE(I,ILOW)=MAX(I-1,IABS(IOPHL2))
         M2FILE(I,ILOW)=MIN(I+3,NLMX)
      END DO
      IF(IABS(IOPHL2).GT.100) THEN
         IOPHL2=MOD(IOPHL2,100)
         ISET=0
  140    CONTINUE
         READ(IBUFF,*,ERR=190) IL1,IU1,IM1,IP1
         ISET=ISET+1
         IF(IL1.LE.0.AND.ISET.EQ.1) THEN
            IL1=1
            IU1=6
            IM1=0
            IP1=1
         END IF
         IF(IL1.LE.0.AND.ISET.EQ.2) THEN
            IL1=7
            IU1=NLMX
            IM1=1
            IP1=3
         END IF
         IUP1=MIN(IU1,NLMX)
         DO I=IL1,IUP1
            M1FILE(I,ILOW)=MAX(I-IM1,IABS(IOPHL2))
            M2FILE(I,ILOW)=MIN(I+IP1,NLMX)
         END DO
         IF(IU1.LT.NLMX) GO TO 140
  190    CONTINUE
         READ(IBUFF,*,ERR=200,END=200) ABLIM1,ABLIM2,ABLIM3
         IF(ABLIM1.LE.0) ABLIM1=6650.
         IF(ABLIM2.LE.0) ABLIM2=5000.
         IF(ABLIM3.LE.0) ABLIM3=6500.
      END IF
  200 CONTINUE
      M1FILE(NLMX,ILOW)=NLMX+1
      M2FILE(NLMX,ILOW)=NLMX
c
C    -------------------
C    Stark paramereters
C    -------------------
C
      izzh=1
      IF(IOPHL1.NE.0) THEN
         I=1
         I1=MAX(2,IABS(IOPHL1))
         DO J=I1,NLMX
            CALL STARK0(I,J,izzh,XKIJ(I,J),WL0(I,J),FIJ(I,J))
         END DO
      END IF
      IF(IOPHL2.NE.0) THEN
         I=2
         I2=MAX(3,IABS(IOPHL2))
         DO J=I2,NLMX
            CALL STARK0(I,J,izzh,XKIJ(I,J),WL0(I,J),FIJ(I,J))
         END DO
      END IF
      RETURN
      END
C
C
C ********************************************************************
C
C
      SUBROUTINE WNSTOR(ID)
C     =====================
C
C     Stores occupation probabilities for hydrogen levels
C     in common WNCOM for further use
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'ATOMIC.FOR'
      INCLUDE 'MODELQ.FOR'
      PARAMETER (SIXTH=UN/6.,CCOR=0.09)
      parameter (p1=0.1402,p2=0.1285,p3=un,p4=3.15,p5=4.)
      parameter (tkn=3.01,ckn=5.33333333,cb0=8.59d14,f23=-2./3.)
      if(ioptab.lt.0) return
C
      cb=cb0*bergfc
      ANE=ELEC(ID)
      A=CCOR*EXP(SIXTH*LOG(ANE))/SQRT(TEMP(ID))
      z=un
      x=exp(p4*log(un+p3*a))
      c1=p1*(x+p5*(z-un)*a*a*a)
      c2=p2*x
      beta0=cb*z*z*z*exp(f23*log(ane))
      DO 20 I=1,NLMX
         XN=I
         if(xn.le.tkn) then
            xkn=un
         else
            xn1=un/(xn+un)
            xkn=ckn*xn*xn1*xn1
         end if
         beta=beta0*xkn*xi2(i)*xi2(i)
         f=(c1*beta*beta*beta)/(un+c2*beta*sqrt(beta))
         WNHINT(I,ID)=f/(un+f)
   20 CONTINUE
C
C     array WOP - occupation probabilities for explicit levels
C       (if ifwop>1, occ. probabilities have been initialized
C          for iron-peak elements and are not updated)
C
      do ii=1,nlevel
      
        if(ifwop(ii).le.0) then
          wop(ii,id)=un
        
        else if(ifwop(ii).eq.1) then
          ie=iel(ii)
          nq=nquant(ii)
          if(iz(ie).eq.1) then
            wop(ii,id)=wnhint(nq,id)
          else
            z=iz(ie)
            xn=nq
            wop(ii,id)=wn(xn,a,ane,z)
          end if
        
        end if
      
        if(ifwop(ii).gt.1.and.lte) wop(ii,id)=un
      
      end do
      
      RETURN
      END
C
C
C ********************************************************************
C
C

      function wn(xn,a,ane,z)
c
c     evaluation of the occupation probablities for a hydrogenic ion
c     using eqs (4.26), and (4.39) of Hummer,Mihalas Ap.J. 331, 794, 1988.
c     approximate evaluation of Q(beta) - Hummer
c
c     Input: xn  - real number corresponding to quantum number n
c            a   - correlation parameter
c            ane - electron density
c            z   - ionic charge
c
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      parameter (p1=0.1402,p2=0.1285,p3=un,p4=3.15,p5=4.)
      parameter (tkn=3.01,ckn=5.33333333,cb0=8.59d14)
      parameter (f23=-2./3.)
c
      cb=cb0*bergfc
c
c     evaluation of k(n)
c
      if(xn.le.tkn) then
         xkn=un
       else
         xn1=un/(xn+un)
         xkn=ckn*xn*xn1*xn1
      end if
c
c     evaluation of beta
c
      beta=cb*z*z*z*xkn/(xn*xn*xn*xn)*exp(f23*log(ane))
c
c     approximate expression for Q(beta)
c
      x=exp(p4*log(un+p3*a))
      c1=p1*(x+p5*(z-un)*a*a*a)
      c2=p2*x
      f=(c1*beta*beta*beta)/(un+c2*beta*sqrt(beta))
      wn=f/(un+f)
      return
      end
C
C
C ********************************************************************
C
C

      SUBROUTINE DWNFR(MODE,N,FRE,A,ANE,Z,FR,DW)
C     ==========================================
C
C     Auxiliary routine to compute set of dissolved fractions
C     for all frequencies
C      MODE=0  ->  DW=1
C      MODE>0  ->  DW=1-w
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'MODELQ.FOR'
      parameter (p1=0.1402,p2=0.1285,p3=un,p4=3.15,p5=4.)
      parameter (tkn=3.01,ckn=5.33333333,cb0=8.59d14,f23=-2./3.)
      PARAMETER (FRH=3.28805D15,SQFRH=5.734152D7)
      DIMENSION FR(N),DW(N)
C
      cb=cb0*berfc
c
      IF(MODE.EQ.0) THEN
       DO 20 IJ=1,N
         DW(IJ)=UN
   20  CONTINUE
      ELSE
       DO 10 IJ=1,N
         IF(FR(IJ).LT.FRE) THEN
            XN=SQFRH*Z/SQRT(FRE-FR(IJ))
            if(xn.le.tkn) then
               xkn=un
            else
               xn1=un/(xn+un)
               xkn=ckn*xn*xn1*xn1
            end if
            beta=cb*z*z*z*xkn/(xn*xn*xn*xn)*exp(f23*log(ane))
            x=exp(p4*log(un+p3*a))
            c1=p1*(x+p5*(z-un)*a*a*a)
            c2=p2*x
            f=(c1*beta*beta*beta)/(un+c2*beta*sqrt(beta))
            DW(IJ)=UN-f/(un+f)
         ELSE
            DW(IJ)=UN
         END IF
   10  CONTINUE
      END IF
      RETURN
      END
C
C
C ********************************************************************
C
C

      SUBROUTINE ODF1(IMODE,IL,IU,ID,ODF)
C     ===================================
C
C     opacity distribution function for overlapping lines near the series limit
C
C     The lines converge to the edge of the (continuum) transition IL - IU,
C     IL - index of the lower level
C     IU - index of the upper level (usually the ground state of teh next ion)
C     ID - depth index
C
C     Output: ODF - opacity distribution function interpolated to the set of
C                   explicit frequencies
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'ATOMIC.FOR'
      INCLUDE 'MODELQ.FOR'
      INCLUDE 'ODFPAR.FOR'
      PARAMETER (FRH=3.28805D15,CQT=1.284523D12)
      PARAMETER (CCOR=0.09,C00=1.25D-9,CID=0.02654,SIXTH=UN/6.)
      DIMENSION FRO(MFRO),ODF0(MFRO),ABSO0(MFRO),ODF(MFREQ),SGT(MFRO),
     *          ALAM(MFRO),
     *          FROD(MFRO),SGFR(MFRO),IODF(MFRO),IODR(MFRO),DWF(MFRO)
      SAVE FRO,SGFR,ODF0,IODF
C
      I=NQUANT(IL)
      KL=INDODF(IL)
      IELO=IEL(IL)
      N1H=NLAST(IELO)
      NQ1=NQLODF(IL)
      FRE=ENION(IL)/H
      T=TEMP(ID)
      SQT=SQRT(T)
      ANE=ELEC(ID)
      ANES=EXP(SIXTH*LOG(ANE))
c      ACOR=CCOR*ANES/SQT 
      F00=C00*ANES*ANES*ANES*ANES
      DOP0=CQT*SQT
      QZ=IZ(IELO)
C
C     pseudocontinuum opacity (non-zero in all frequencies);
C     formulated through the dissolved fraction
C
      ITR=ITRA(IL,IU)
      NFR0=NFRODF(KL)
      IF(IMODE.EQ.0) THEN
         DO IJ=1,NFR0
            FRO(IJ)=FROS(IJ,KL)
            SGFR(IJ)=SIGK(FRO(IJ),ITR,1)
          ALAM(IJ)=CAS/FRO(IJ)
         END DO
      END IF
C        
C        function D(nu) - dissolved fraction
C        
c      CALL DWNFR(1,NFR0,FRE,ACOR,ANE,QZ,FRO,DWF)
      DO 20 IJ=1,NFR0
         ABSO0(IJ)=SGFR(IJ)*DWF(IJ)
   20 CONTINUE
C
C     summation over individual lines
C
      DO 40 J=NQ1,NLMX
         XJ=J
         FXK=F00*XKIJ(KL,J)
         DOP=DOP0/WL0(KL,J)
         DBETA=WL0(KL,J)*WL0(KL,J)/CAS/FXK
         BETAD=DOP*DBETA
         FID=CID*FIJ(KL,J)*DBETA
         CALL DIVSTR(1)
         WPROB=WNHINT(J,ID)
         CALL ODFHST(NFR0,FXK,FID,WPROB,WL0(KL,J),ALAM,SGT)
         DO 30 IJ=1,NFR0
   30       ABSO0(IJ)=ABSO0(IJ)+SGT(IJ)
   40 CONTINUE
C
C     opacity distribution function in the internal set of frequencies
C
      IF(IMODE.EQ.0) THEN
         ODF0(1)=ABSO0(1)
         IODF(1)=1
         DO 60 IJ=2,NFR0
            ODF0(IJ)=ABSO0(IJ)
            IODF(IJ)=IJ
            IF(ODF0(IJ).LT.ODF0(IJ-1)) THEN
               AB=ODF0(IJ)
               IJODF=IODF(IJ)
               DO 70 IJ0=1,IJ-1
                  IJ1=IJ-IJ0+1
                  IF(ODF0(IJ1).GE.ODF0(IJ1-1)) GO TO 71
                  ODF0(IJ1)=ODF0(IJ1-1)
                  ODF0(IJ1-1)=AB
                  IODF(IJ1)=IODF(IJ1-1)
                  IODF(IJ1-1)=IJODF
   70          CONTINUE
   71          CONTINUE
            END IF
         if(odf0(ij).gt.0.001) write(6,603) ij,id,odf0(ij)
  603    format(' ij,id,odf0',2i5,1pd10.3)
   60    CONTINUE
      ELSE
         ODF0(1)=ABSO0(IODF(1))
         IODR(1)=IODF(1)
         DO 80 IJ=2,NFR0
            ODF0(IJ)=ABSO0(IODF(IJ))
            IODR(IJ)=IODF(IJ)
            IF(ODF0(IJ).LT.ODF0(IJ-1)) THEN
               AB=ODF0(IJ)
               IJODF=IODR(IJ)
               DO 85 IJ0=1,IJ-1
                  IJ1=IJ-IJ0+1
                  IF(ODF0(IJ1).GE.ODF0(IJ1-1)) GO TO 86
                  ODF0(IJ1)=ODF0(IJ1-1)
                  ODF0(IJ1-1)=AB
                  IODR(IJ1)=IODR(IJ1-1)
                  IODR(IJ1-1)=IJODF
   85          CONTINUE
   86          CONTINUE
            END IF
         if(odf0(ij).gt.0.001) write(6,603) ij,id,odf0(ij)
   80    CONTINUE
         DO 90 IJ=1,NFR0
            IODF(IJ)=IODR(IJ)
   90    CONTINUE
      ENDIF
C
C     Reinitialization of the internal frequencies set
C
      FROD(1)=FRO(1)
      IW=IODF(1)
      IF(IW.GT.1 .AND. IW.LT.NFR0) THEN
         W1=FRO(IW-1)-FRO(IW+1)
      ELSE IF (IW.EQ.1) THEN
         W1=FRO(1)-FRO(2)
      ELSE
         W1=FRO(NFR0-1)-FRO(NFR0)
      ENDIF
      DO 200 IJ=2,NFR0-1
         IW=IODF(IJ)
         IF(IW.GT.1 .AND. IW.LT.NFR0) THEN
            W2=HALF*(FRO(IW-1)-FRO(IW+1))
         ELSE IF (IW.EQ.1) THEN
            W2=HALF*(FRO(1)-FRO(2))
         ELSE
            W2=HALF*(FRO(NFR0-1)-FRO(NFR0))
         ENDIF
         FROD(IJ)=FROD(IJ-1)-HALF*(W1+W2)
         W1=W2
  200 CONTINUE
      IJ=NFR0
      IW=IODF(IJ)
      IF(IW.GT.1 .AND. IW.LT.NFR0) THEN
         W2=FRO(IW-1)-FRO(IW+1)
      ELSE IF (IW.EQ.1) THEN
         W2=FRO(1)-FRO(2)
      ELSE
         W2=FRO(NFR0-1)-FRO(NFR0)
      ENDIF
      FROD(IJ)=FROD(IJ-1)-HALF*(W1+W2)
C
C     Interpolated opacity distribution function to explicit frequencies
C
      DO 150 IJ=2,NFREQ
         IF(FREQ(IJ).GT.FREQ(IJ-1)) RETURN
         ODF(IJ)=0.
         IF(FREQ(IJ).GT.FROD(1).OR.FREQ(IJ).LT.FROD(NFR0)) GO TO 150
         IF(ID.EQ.1) THEN
            IF(FREQ(IJ-1).GT.FROD(1)) I1ODF(IL)=IJ
            I2ODF(IL)=IJ
         END IF
         DO 110 IJ1=2,NFR0
            IJ0=IJ1
            IF(FREQ(IJ).GE.FROD(IJ1)) GO TO 120
  110    CONTINUE
  120    ODF(IJ)=ODF0(IJ0-1)+(ODF0(IJ0)-ODF0(IJ0-1))/
     *           (FROD(IJ0)-FROD(IJ0-1))*(FREQ(IJ)-FROD(IJ0-1))
  150 CONTINUE
      RETURN
      END
C
C
C ********************************************************************
C
C
      SUBROUTINE ODFHST(N,FXK,FID,WP,WL,ALAM,SG)
C     ==========================================
C
C     Auxiliary routine for ODF1 (replaces multiple calls to STARKA)
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'MODELQ.FOR'
      INCLUDE 'ODFPAR.FOR'
      PARAMETER (F0=-0.5758228,F1=0.4796232,F2=0.07209481,AL=1.26)
      PARAMETER (SD=0.5641895,SLO=-2.5,THRA=1.5,BL1=1.14,BL2=11.4)
      PARAMETER (SAC=0.08, THR=THRA*TWO)
      DIMENSION ALAM(MFRO),SG(MFRO)
C
      BETAD1=UN/BETAD
      FXK1=UN/FXK
      FIDWP=FID*WP
C
C     for a > 1 Doppler core + asymptotic Holtzmark wing with division
C               point DIV
C
      IF(ADH.GT.AL) THEN
         DO IJ=1,N
            BETA=ABS(ALAM(IJ)-WL)*FXK1
            XD=BETA*BETAD1
            IF(XD.LE.DIVH) THEN
               ST=SD*EXP(-XD*XD)*BETAD1
             ELSE
               ST=THR*EXP(SLO*LOG(BETA))
            END IF
            SG(IJ)=ST*FIDWP
         END DO
       ELSE
C
C      empirical formula for a < 1
C
         DO IJ=1,N
            BETA=ABS(ALAM(IJ)-WL)*FXK1
            XD=BETA*BETAD1
            IF(BETA.LE.BL1) THEN
               ST=SAC
             ELSE IF(BETA.LT.BL2) THEN
               XL=LOG(BETA)
               FL=(F0*XL+F1)*XL
               ST=F2*EXP(FL)
             ELSE
               ST=THR*EXP(SLO*LOG(BETA))
            END IF
            SG(IJ)=ST*FIDWP
         END DO
      END IF
C
      RETURN
      END
C
C
C ********************************************************************
C
C
      SUBROUTINE ODFFR(IL,IU)
C     =======================
C
C     Internal frequencies set for opacity distribution function
C     for overlapping lines near the series limit
C
C     The lines converge to the edge of the (continuum) transition IL - IU,
C     IL - index of the lower level
C     IU - index of the upper level (usually the ground state of the next ion)
C          or a mean level in the line ODF formalism 
C
C     Output: FROS - set of internal frequencies
C                    in common ODFFRQ
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'ATOMIC.FOR'
      INCLUDE 'MODELQ.FOR'
      INCLUDE 'ODFPAR.FOR'
      PARAMETER (FRH=3.28805D15,CDOP=2.84511D-7,CDOM=14.)
      PARAMETER (SIX=6.,SEPT=7.)
      DIMENSION FFRO(MFRO)
C
      CH=IZ(IEL(IL))*IZ(IEL(IL))
      FRION=CH*FRH
      FRE=ENION(IL)/H
      NF=1
        NQ1=NQUANT(IU)
        XL2=UN/(NQUANT(IL)*NQUANT(IL))
        XU1=UN/((NQUANT(IU)-1)*(NQUANT(IU)-1))
        XU2=UN/(NQUANT(IU)*NQUANT(IU))
        FRC=FRION*(XL2-XU2)
        FFRO(NF)=HALF*(FRC+FRION*(XL2-XU1))
        KT=ITRA(IL,IU)
        KL=JNDODF(KT)
        DOPO=CDOP*SQRT(TEFF)*FRC
      DOPM=CDOM*DOPO
      FR1=FFRO(1)
C
      DO 10 I=NQ1,NLMX
         II=I*I
         FR2=FRE-FRION/II
         DF=FR2-FR1
         IF(DF.GT.DOPM) THEN
            DO 20 J=1,7
               NF=NF+1
               FFRO(NF)=FR1+J*DOPO
   20       CONTINUE
            DF=FR2-SEPT*DOPO-FFRO(NF)
            NI=int(DF/SIX/DOPO)
            DDF=DF/(NI+1)
            DO 21 J=1,NI
               NF=NF+1
               FFRO(NF)=FR1+SEPT*DOPO+J*DDF
   21       CONTINUE
            DO 22 J=7,0,-1
               NF=NF+1
               FFRO(NF)=FR2-J*DOPO
   22       CONTINUE
            FR1=FR2
         ELSE
            NI=int(DF/DOPO)
            DDF=DF/(NI+1)
            DO 30 J=1,NI
               NF=NF+1
             if(nf.gt.mfro-3) then
                nf=nf-1
                go to 15
             end if
               FFRO(NF)=FR1+J*DDF
   30       CONTINUE
            NF=NF+1
            FFRO(NF)=FR2
            FR1=FR2
         ENDIF
   10 CONTINUE
   15 CONTINUE
      NF=NF+1
      FFRO(NF)=FRE*0.999999999
      NFRODF(KL)=NF
      if(nf.gt.mfro) 
     *   CALL QUIT('too many points for hydrogen ODF - nf.gt.mfro',
     *             nf,mfro)
      DO 40 I=1,NF
         FROS(I,KL)=FFRO(NF-I+1)
   40 CONTINUE
C
C   Associated weights
C
      WNUS(1,KL)=HALF*(FROS(1,KL)-FROS(2,KL))
      WNUS(NF,KL)=HALF*(FROS(NF-1,KL)-FROS(NF,KL))
      DO 50 I=2,NF-1
        WNUS(I,KL)=HALF*(FROS(I-1,KL)-FROS(I+1,KL))
   50 CONTINUE
C
      RETURN
      END
C
C
C ********************************************************************
C
C
      SUBROUTINE TIMING(MOD,ITER)
C     ===========================
C
C     Timing procedure (call machine dependent routine!!)
C
      CHARACTER ROUT*20
      dimension dummy(2)
      DATA T0/0./
      SAVE T0
C
      TIME=etime(dummy)
      DT=TIME-T0
      T0=TIME
      IF(MOD.EQ.1) THEN
         IP=ITER-1
         ROUT='     FORMAL SOLUTION'
      ELSE IF(MOD.EQ.2) THEN
         IP=ITER
         ROUT='       LINEARIZATION'
      ENDIF
      WRITE(69,600) IP,MOD,TIME,DT,ROUT
  600 FORMAT(2I4,2F11.2,2X,A20)
      RETURN
      END
C
C
C ********************************************************************
C
C
      subroutine quit(text,i1,i2)
c
c     stops the program and writes a text
c
      character*(*) text
c     write(6,10) text,i1,i2
      write(10,10) text,i1,i2
   10 format(1x,a,2x,2i10)
c     stop
      end
C
C
C ********************************************************************
C
C
      SUBROUTINE ODFSET
C     =================
C
C     Initialization of line ODF's
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'ATOMIC.FOR'
      INCLUDE 'MODELQ.FOR'
      INCLUDE 'ODFPAR.FOR'
      COMMON/STFCR/OFR(MFODF),OW(MFODF),OWSUB(MFODF),
     *          ODFL0(MDODF,MFODF),ODF2(MDEPTH),IFTRA(MTRANS),
     *          IDODF(MDODF),NDODF
      DIMENSION DML(MDEPTH)
C
      IDSTD=ND*2/3
      NLASTE=NFREQ
      ITR0=0
      DO ID=1,ND
         IF(DM(ID).GT.0) THEN
            DML(ID)=LOG(DM(ID))
          ELSE
            DML(ID)=ID
         END IF
      END DO
C
      DO 500 ION=1,NION
         IND=INODF1(ION)
         IF(IND.LE.0) GO TO 500
         IND2=INODF2(ION)
         IF(FIODF1(ION).NE.' ') OPEN(IND,FILE=FIODF1(ION),STATUS='OLD')
         IF(FIODF2(ION).NE.' ') OPEN(IND2,FILE=FIODF2(ION),STATUS='OLD')
         READ(IND,*,END=500) NDODF
         IF(NDODF.GT.MDODF) 
     *   CALL QUIT('too many depths for an ODF - ndodf.gt.mdodf',
     *             ndodf,mdodf)
         READ(IND,*) (IDODF(ID),ID=1,NDODF)
         IREC=0
   10    CONTINUE
         READ(IND,*,END=500) II,JJ,FR,NFRO,FAV
         IF(NFRO.GT.MFODF) 
     *   CALL QUIT('too many frequencies for an ODF - nfro.gt.mfodf',
     *             nfro,mfodf)
         DO IJ=1,NFRO
            READ(IND,*) OFR(IJ),OW(IJ),OWSUB(IJ)
         END DO
         IND2=INODF2(ION)
         READ(IND2,*) ((ODFL0(ID,IF),ID=1,NDODF),IF=1,NFRO)
         IREC=IREC+1
C
         N0=NFIRST(ION)-1
         I=II+N0
         J=JJ+N0
         IF(J.GT.NLAST(ION)) GO TO 10
         IF(I.GE.NLAST(ION)) GO TO 500
         ITR=ITRA(I,J)
         IF(ITR.EQ.ITR0) THEN
            ITR1=0
            IF(IF1.EQ.1) THEN
               IFIJ=0
               DO 30 IT=1,NTRANS
                  IF(ILOW(IT).NE.I.OR.IUP(IT).NE.J) GO TO 30
                  IF(IT.EQ.ITR) GO TO 30
                  IFIJ=IFIJ+1
                  IFTRA(IT)=IFIJ
   30          CONTINUE
               IF1=0
            END IF
            DO IT=1,NTRANS
               IF(IFTRA(IT).GT.0) THEN
                  ITR1=IT
                  GO TO 50
               END IF
            END DO
   50       CONTINUE
            IF(ITR1.EQ.0) THEN
               WRITE(6,601) ITR,N0,II,JJ
               STOP
  601          FORMAT(' CONFLICT IN ODF INPUT; ITR=',4I5)
            END IF
            ITR=ITR1
            IFTRA(ITR)=0
            OSC0(ITR)=FAV
          ELSE
            ITR0=ITR
            IF1=1
            OSC0(ITR)=FAV
         END IF
C
         MODE=IABS(INDEXP(ITR))
         IF(MODE.EQ.3.OR.MODE.EQ.4) THEN
            LCOMP(ITR)=.FALSE.
            INTMOD(ITR)=5
         END IF
         IFRQ0=IFR0(ITR)
         IFRQ1=IFR1(ITR)
         IF(OFR(1).GE.OFR(NFRO)) THEN
            IF(MODE.EQ.3) THEN
               IFR0(ITR)=NLASTE+1
               IFR1(ITR)=NLASTE+NFRO
               DO IJ=1,NFRO
                  FREQ(IJ+NLASTE)=OFR(IJ)
                  W(IJ+NLASTE)=OW(IJ)
               END DO
               IF(NDODF.EQ.1) THEN
                  DO ID=1,ND
                     DO IJ=1,NFRO
                        PRFLIN(ID,IJ+NLASTE)=real(ODFL0(1,IJ))
                     END DO
                  END DO
                ELSE
                  DO ID=1,ND
                     ID1=1
                     DO IDO=1,NDODF-1
                        IF(ID.GE.IDODF(IDO).AND.ID.LE.IDODF(IDO+1)) THEN
                           ID1=IDO
                          ID2=IDO+1
                          GO TO 140
                        END IF
                     END DO                           
  140                CONTINUE 
                     IF(ID2.GT.NDODF) ID2=NDODF
                     IF(ID1.EQ.ID2) THEN
                        A1=1.
                        A2=0.
                      ELSE
                        X=DML(IDODF(ID2))-DML(IDODF(ID1))
                        A1=(DML(IDODF(ID2))-DML(ID))/X                          
                        A2=UN-A1
                     END IF
                     DO IJ=1,NFRO
                        IF(ODFL0(ID1,IJ).LE.0.OR.
     *                     ODFL0(ID2,IJ).LE.0) THEN
                           PRFLIN(ID,IJ+NLASTE)=0.
                         ELSE
                           X=EXP(A1*LOG(ODFL0(ID1,IJ))+
     *                       A2*LOG(ODFL0(ID2,IJ)))
                           PRFLIN(ID,IJ+NLASTE)=real(X)
                        END IF
                     END DO                         
                  END DO 
               END IF                        
C
               IF(IPROF(ITR).EQ.0) THEN
                  DO ID=1,ND
                     PRFLIN(ID,IFR1(ITR))=0.
                  END DO
               END IF
               DO IJ=1,NFRO
                  PROF(IJ+NLASTE)=PRFLIN(IDSTD,IJ+NLASTE)
               END DO 
               NLASTE=IFR1(ITR)
            END IF
          ELSE
            IF(MODE.EQ.3) THEN
               IFR0(ITR)=NLASTE+1
               IFR1(ITR)=NLASTE+NFRO
               DO IJ=1,NFRO
                  FREQ(IJ+NLASTE)=OFR(NFRO-IJ+1)
                  W(IJ+NLASTE)=OW(NFRO-IJ+1)
               END DO
               IF(NDODF.EQ.1) THEN
                  DO ID=1,ND
                     DO IJ=1,NFRO
                        PRFLIN(ID,IJ+NLASTE)=real(ODFL0(1,NFRO-IJ+1))
                     END DO
                  END DO
                ELSE
                  DO ID=1,ND
                     ID1=1
                     DO IDO=1,NDODF-1
                        IF(ID.GE.IDODF(IDO).AND.ID.LE.IDODF(IDO+1)) THEN
                           ID1=IDO
                           ID2=IDO+1
                           GO TO 240
                        END IF
                     END DO                           
  240                CONTINUE 
                     IF(ID2.GT.NDODF) ID2=NDODF
                     IF(ID1.EQ.ID2) THEN
                        A1=1.
                        A2=0.
                      ELSE
                        X=DML(IDODF(ID2))-DML(IDODF(ID1))
                        A1=(DML(IDODF(ID2))-DML(ID))/X                          
                        A2=UN-A1
                     END IF
                     DO IJ=1,NFRO
                        IJ0=NFRO-IJ+1
                        IF(ODFL0(ID1,IJ0).LE.0.OR.ODFL0(ID2,IJ0).LE.0) 
     *                     THEN
                           PRFLIN(ID,IJ+NLASTE)=0.
                         ELSE
                           X=EXP(A1*LOG(ODFL0(ID1,IJ0))+
     *                     A2*LOG(ODFL0(ID2,IJ0)))
                           PRFLIN(ID,IJ+NLASTE)=REAL(X)
                        END IF
                     END DO                         
                  END DO 
               END IF                        
C
               IF(IPROF(ITR).EQ.0) THEN
                  DO ID=1,ND
                     PRFLIN(ID,IFR0(ITR))=0.
                  END DO
               END IF
               DO IJ=1,NFRO
                  PROF(IJ+NLASTE)=PRFLIN(IDSTD,IJ+NLASTE)
               END DO 
               NLASTE=IFR1(ITR)
            END IF
         END IF
         IF(NLASTE.GT.MFREQ)
     *   CALL QUIT(' too many frequencies in ODFSET - nlaste.gt.mfreq',
     *              nlaste,mfreq)
         IF(INDEXP(ITR).NE.0) THEN
            CALL IJALIS(ITR,IFRQ0,IFRQ1)
         END IF
         GO TO 10
  500 CONTINUE
C
      NFREQ=NLASTE
      RETURN
      END
C
C
C ********************************************************************
C
C
      SUBROUTINE ODFHYS(DOPO)
C     =======================
C
C     Initialization of line ODF's for hydrogen
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'ATOMIC.FOR'
      INCLUDE 'MODELQ.FOR'
      INCLUDE 'ODFPAR.FOR'
      PARAMETER (CCM=UN/2.997925D10,THIRD=UN/3.,FRH=3.28805D15)
      DIMENSION FFRO(MFRO)
C
      izzh=1
      IF(ISPODF.GE.1) THEN
         DO 200 ITR=1,NTRANS
            JND=JNDODF(ITR)
            MODE=IABS(INDEXP(ITR))
            IF(JND.LE.0 .OR. MODE.NE.2) GO TO 200
            LINEXP(ITR)=.FALSE.
            LCOMP(ITR)=.FALSE.
            INTMOD(ITR)=6
            I=ILOW(ITR)
            J=IUP(ITR)
            NQLODF(I)=IABS(IPROF(ITR))
            IF(NQLODF(I).EQ.0) NQLODF(I)=NQUANT(J)
            OSC0(ITR)=0.
            IS=NQUANT(I)
            DO K=NQUANT(J),NLMX
               CALL STARK0(IS,K,izzh,XKIJ(JND,K),WL0(JND,K),
     *                    FIJ(JND,K))
               OSC0(ITR)=OSC0(ITR)+FIJ(JND,K)
            END DO
  200    CONTINUE
         RETURN
      END IF
C
      NLASTE=NFREQ
      DO 100 ITR=1,NTRANS
         JND=JNDODF(ITR)
         MODE=IABS(INDEXP(ITR))
         IF(JND.LE.0 .OR. MODE.NE.2) GO TO 100
         LCOMP(ITR)=.FALSE.
         INTMOD(ITR)=6
         I=ILOW(ITR)
         J=IUP(ITR)
         NQLODF(I)=IABS(IPROF(ITR))
         IF(NQLODF(I).EQ.0) NQLODF(I)=NQUANT(J)
         XJ2A=HALF*(XI2(NQUANT(J))+XI2(NQUANT(J)-1))
C
C        set up explicit frequencies & weights
C
         NFRO=0
         DO IFQ=1,4
            NFRO=NFRO+KDO(IFQ,JND)
         END DO    
         NFRO=NFRO-2
         FRION=FRH*IZ(IEL(I))*IZ(IEL(I))
         FRA=FRION*(XI2(NQUANT(I))-XJ2A)
         DOPI=DOPO*FRA*CCM
         FRB=0.99999999*FRION*XI2(NQUANT(I))
         IFRQ0=IFR0(ITR)
         IFRQ1=IFR1(ITR)
         IFR0(ITR)=NLASTE+1
         IFR1(ITR)=NLASTE+NFRO
         I1ODF(I)=IFR0(ITR)
         I2ODF(I)=IFR1(ITR)-1
         FFRO(1)=0.99999999*FRA
         FFRO(2)=FRA
         IJ00=1
         DO IK=1,3
            DO IJ=2,KDO(IK,JND)
               IJQ=IJ00+IJ
               FFRO(IJQ)=FFRO(IJQ-1)+XDO(IK,JND)*DOPI
            END DO
            IJ00=IJ00+KDO(IK,JND)-1
         END DO
         do ij=1,ij00
            if(ffro(ij).lt.frb) nfrb=ij
         end do
         if(nfrb.eq.ij00) then
            IJ00=IJ00+1
            FFRO(NFRO)=0.99999999*FRION*XI2(NQUANT(I))
            do while (ffro(ij00).ge.ffro(nfro))
               xdo(3,jnd)=0.75*xdo(3,jnd)
               ij00=ij00-kdo(3,jnd)
               do ij=2,kdo(3,jnd)
                  ijq=ij00+ij
                  ffro(ijq)=ffro(ijq-1)+xdo(3,jnd)*dopi
               end do
              ij00=ij00+kdo(3,jnd)
            enddo
            TIDO=(FFRO(NFRO)-FFRO(IJ00))/FLOAT(KDO(4,JND)-1)
            DO IJ=1,KDO(4,JND)-2
               IJQ=NFRO-IJ
               FFRO(IJQ)=FFRO(NFRO)-FLOAT(IJ)*TIDO
            END DO
          else
            TIDO=(FRB-FFRO(nfrb))*third
            ffro(nfrb+1)=FFRO(nfrb)+tido
            ffro(nfrb+2)=frb-tido
            ffro(nfrb+3)=frb
            nfro=nfrb+3
            IFR1(ITR)=NLASTE+NFRO
            I2ODF(I)=IFR1(ITR)-1
         endif
         DO IJ=1,NFRO
            FREQ(NLASTE+IJ)=FFRO(NFRO-IJ+1)
         END DO
         W(NLASTE+NFRO)=HALF*(FREQ(NLASTE+NFRO-1)-FREQ(NLASTE+NFRO))
         W(NLASTE+NFRO-1)=W(NLASTE+NFRO)
         DO IJ=2,NFRO-2,2
            TIDO=(FREQ(NLASTE+IJ)-FREQ(NLASTE+IJ+1))*THIRD
            W(NLASTE+IJ-1)=W(NLASTE+IJ-1)+TIDO
            W(NLASTE+IJ)=W(NLASTE+IJ)+4.*TIDO
            W(NLASTE+IJ+1)=W(NLASTE+IJ+1)+TIDO
         END DO
         NLASTE=IFR1(ITR)
C
C        set up internal frequencies & Stark parameters
C
         CALL ODFFR(I,J)
         OSC0(ITR)=0.
         IS=NQUANT(I)
         DO K=NQUANT(J),NLMX
            CALL STARK0(IS,K,izzh,XKIJ(JND,K),WL0(JND,K),FIJ(JND,K))
            OSC0(ITR)=OSC0(ITR)+FIJ(JND,K)
         END DO
         IF(INDEXP(ITR).NE.0) THEN
            CALL IJALIS(ITR,IFRQ0,IFRQ1)
         END IF
  100 CONTINUE
C
      NFREQ=NLASTE
      RETURN
      END
C
C
C ********************************************************************
C
C
      SUBROUTINE ODFMER
C     =================
C
C     Opacity distribution function for superlines to merged states
C     (calculated only if DT/T>=CHTL, kept constant after).
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'ATOMIC.FOR'
      INCLUDE 'MODELQ.FOR'
      INCLUDE 'ODFPAR.FOR'
      PARAMETER (CHTL=1.D-3)
C
      DO 10 ITR=1,NTRANS
         IF(.NOT.LINE(ITR).OR.IABS(INDEXP(ITR)).NE.2) GO TO 10
         DO ID=1,ND
            IF(INIT.EQ.1 .OR. ABS(CHANT(ID)).GE.CHTL)
     *         CALL ODFHYD(ID,ITR)
         END DO
   10 CONTINUE
      RETURN
      END
C
C
C
C ********************************************************************
C
C
      SUBROUTINE ODFHYD(ID,ITR)
C     =========================
C
C     Line ODF's for hydrogen line series
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'ATOMIC.FOR'
      INCLUDE 'MODELQ.FOR'
      INCLUDE 'ODFPAR.FOR'
      PARAMETER (CDOP=TWO*BOLK/HMASS)
      PARAMETER (CA=2.997925D18,CCM=CA/1.D8,FRH=3.28805D15)
      PARAMETER (RYDEL=911.764,TTW=2./3.)
      PARAMETER (C00=1.25D-9,CID=0.02654)
      DIMENSION SIG(MFRO),SGT(MFRO),ODF(MFRO),IODF(MFRO)
      DIMENSION YNUS(MFRO),ALAM(MFRO)
C
      JO=JNDODF(ITR)
      IF(ISPODF.EQ.0) THEN
        NF=NFRODF(JO)
        DO IJ=1,NF
          IODF(IJ)=0
          SIG(IJ)=0.
          ODF(IJ)=0.
        YNUS(IJ)=FROS(IJ,JO)
        ALAM(IJ)=CAS/YNUS(IJ)
        END DO
      ELSE
        NF=IFR1(ITR)-IFR0(ITR)+1
      DO IJ=1,NF
          SIG(IJ)=0.
        YNUS(IJ)=FREQ(IFR0(ITR)+IJ-1)
        ALAM(IJ)=CAS/YNUS(IJ)
        END DO
      END IF
C
      II=ILOW(ITR)
      JJ=IUP(ITR)
      ANES=EXP(TTW*LOG(ELEC(ID)))
      F00=C00*ANES
      FRA=FRH*(XI2(NQUANT(II))-XI2(NQUANT(JJ)))
      DOPO=FRA/CCM*SQRT(CDOP*TEMP(ID)+VTB*VTB)
      DO 30 J=NQLODF(II),NLMX
         WL=RYDEL/(XI2(NQUANT(II))-XI2(J))
         FXK=F00*XKIJ(JO,J)
         DBETA=WL*WL/CA/FXK
         BETAD=DBETA*DOPO
         FID=CID*FIJ(JO,J)*DBETA
         CALL DIVSTR(1)
         WPROB=WNHINT(J,ID)
         CALL ODFHST(NF,FXK,FID,WPROB,WL,ALAM,SGT)
         DO IJ=1,NF
           SIG(IJ)=SIG(IJ)+SGT(IJ)
         END DO
   30 CONTINUE
C
      IF(ISPODF.EQ.0) THEN
        CALL INDEXX(NF,SIG,IODF)
        DO 40 IJ=1,NF
          ODF(IJ)=SIG(IODF(IJ))
   40   CONTINUE
        I0=IFR0(ITR)
        I1=IFR1(ITR)
        IF(IABS(INDEXP(ITR)).EQ.2) YNUS(1)=FREQ(I0)
        IW1=IODF(1)
        DO 50 IJ=2,NF
          IW2=IODF(IJ)
          IF(IJ.GT.2 .AND. IJ.LT.NF) THEN
            YNUS(IJ)=YNUS(IJ-1)-HALF*(WNUS(IW1,JO)+WNUS(IW2,JO))
          ELSE IF (IJ.EQ.2) THEN
            YNUS(IJ)=YNUS(IJ-1)-HALF*(TWO*WNUS(IW1,JO)+WNUS(IW2,JO))
          ELSE IF (IJ.EQ.NF) THEN
            YNUS(IJ)=YNUS(IJ-1)-HALF*(WNUS(IW1,JO)+TWO*WNUS(IW2,JO))
          ENDIF
          IW1=IW2
   50   CONTINUE
      END IF
C
      IF(ISPODF.EQ.0) THEN
        PRFLIN(ID,I1)=1.E-35
        DO 80 IJQ=I0,I1-1
           DO 60 IJ=2,NF
              JI=IJ
              IF(YNUS(IJ).LE.FREQ(IJQ)) GO TO 70
   60      CONTINUE
   70      PRFLN=ODF(JI-1)+(ODF(JI)-ODF(JI-1))*
     *       (FREQ(IJQ)-YNUS(JI-1))/(YNUS(JI)-YNUS(JI-1))
           PRFLIN(ID,IJ0)=real(PRFLN)
   80   CONTINUE
      ELSE
        DO IJ=1,NF
        PRFLIN(ID,KFR0(ITR)+IJ-1)=real(SIG(IJ))
      END DO
      END IF
      RETURN
      END
C
C
C ********************************************************************
C
C
      SUBROUTINE INDEXX(N,ARRIN,INDX)
C     ===============================
C
C     Sorting routine
C
      INCLUDE 'IMPLIC.FOR'
      DIMENSION ARRIN(N),INDX(N)
      DO 11 J=1,N
        INDX(J)=J
11    CONTINUE
      M=N/2+1
      IR=N
10    CONTINUE
        IF(M.GT.1)THEN
          M=M-1
          INDXT=INDX(M)
          Q=ARRIN(INDXT)
        ELSE
          INDXT=INDX(IR)
          Q=ARRIN(INDXT)
          INDX(IR)=INDX(1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            INDX(1)=INDXT
            RETURN
          ENDIF
        ENDIF
        I=M
        J=M+M
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF(ARRIN(INDX(J)).LT.ARRIN(INDX(J+1)))J=J+1
          ENDIF
          IF(Q.LT.ARRIN(INDX(J)))THEN
            INDX(I)=INDX(J)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        INDX(I)=INDXT
      GO TO 10
      END
C
C
C ********************************************************************
C
C
      SUBROUTINE SIGAVE
C     =================
C
C     Read bound-free cross-sections for averaged levels
C     from the unit INSA (given by IFANCY), with increasing frequencies
C     It assumes that all continuum transitions for a given ion are
C     given in a successive order in the data
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'ATOMIC.FOR'
      INCLUDE 'MODELQ.FOR'
      INCLUDE 'ODFPAR.FOR'
      PARAMETER (HCCM=H*2.997925D10)
      PARAMETER (TX=2.30258509299405,BAM=1.D-18)
      DIMENSION XIFE(8)
      DIMENSION FRINSG(MFREQ),CRIN(MFREQ),JKF(MFREQ)
C
      DATA XIFE/63480.,130563.,247220.,442000.,605000.,799000.,
     *          1008000.,1218380./
C
      NFREQB=NFREQ
      IF(IBFINT.GT.0) NFREQB=NFREQC
      ITR=0
   10 ITR=ITR+1
      IF(ITR.GT.NTRANS) RETURN
      IC=ITRA(IUP(ITR),ILOW(ITR))
      INSA=IBF(IC)
      IF(INSA.LT.50 .OR. INSA.GT.100) GO TO 10
      IE=IEL(ILOW(ITR))
      ITR=ITR-1
      NL1=NFIRST(IE)
      NL2=NLAST(IE)
      IF(FIBFCS(IE).NE.' ') THEN
         INSA=INBFCS(IE)
         OPEN(INSA,FILE=FIBFCS(IE),STATUS='OLD')
      END IF
      READ(INSA,*,END=500,ERR=500) IERR,IZRR,NLRR
      DO 100 I=NL1,NL2
        ITR=ITR+1
        IF(INDEXP(ITR).EQ.0) GO TO 100
        IC=ITRA(IUP(ITR),ILOW(ITR))
        READ(INSA,*) INL,ECMR,GDUM,NFIS
        IF(IERR.NE.26) GO TO 20
        ECMR=XIFE(IZRR)-ECMR
        DE=ABS((ENION(ILOW(ITR))-HCCM*ECMR)/ENION(ILOW(ITR)))
        IF(DE.GT.2.D-2)  THEN
         PRINT *,INSA,IE,ITR,I
c        CALL QUIT('Incorrect energy level in SIGAVE; itr,i:',itr,i)
        ENDIF
   20   DO 30 IJ=1,NFIS
          JI=NFIS-IJ+1
          READ(INSA,*,END=500,ERR=500) FRINSG(JI),CRIN(JI)
   30   CONTINUE
        DO 50 IJ=1,NFREQB
          JK=0
        FR=FREQ(IJ)
        IF(ISPODF.GE.1) FR=FREQ(IFREQB(IJ))
          DO 35 IK=1,NFIS
              IF(FR.GT.FRINSG(IK)) THEN
                 JK=IK
                 GO TO 40
              ENDIF
   35     CONTINUE
          JK=NFIS
   40     IF(JK.EQ.1) JK=2
          JKF(IJ)=JK
   50   CONTINUE
        DO 60 IJ=1,NFREQB
          JK=JKF(IJ)
        FR=FREQ(IJ)
        IF(ISPODF.GE.1) FR=FREQ(IFREQB(IJ))
        IF(CRIN(JK-1).EQ.0. .OR. CRIN(JK).EQ.0.) THEN
            BFCS (IC,IJ)=real(CRIN(JK)+(FR-FRINSG(JK))/
     *        (FRINSG(JK-1)-FRINSG(JK))*(CRIN(JK-1)-CRIN(JK)))
          ELSE
            XF1=LOG10(FRINSG(JK-1))
            XF2=LOG10(FRINSG(JK))
            YS1=LOG10(CRIN(JK-1))
            YS2=LOG10(CRIN(JK))
            XXF=LOG10(FR)
            YYF=(XXF-XF2)/(XF1-XF2)*(YS1-YS2)+YS2
            EXTX=EXP(TX*YYF)
            BFCS (IC,IJ)=real(EXTX)
          ENDIF
          BFCS(IC,IJ)=real(BAM)*BFCS(IC,IJ)
   60   CONTINUE 
  100 CONTINUE
      GO TO 10
  500 CALL QUIT
     *('error in data for bf-cs of averaged levels - itr,ie:',
     *          itr,ie)
      RETURN
      END
C
C
C *********************************************************************
C
c
      subroutine readbf(ius)
c     ======================
c
c     auxiliary subroutine for enabling reading of input data with
c     comments
c
c     lines beginning with ! or * are understood as comments
c
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      character*80 buff
c     
      iur=ius 
      if(iur.eq.0) iur=5
      iuw=ibuff
      if(iur.ne.5) iuw=95
c      
   10 continue
      read(iur,501,end=20) buff
      if(buff(1:1).eq.'!'.or.buff(1:1).eq.'*') go to 10
      write(iuw,501) buff
      go to 10
c      
   20 continue
      rewind iuw
      return 
  501 format(a)
      end
C
C
C     *******************************************************************
C
C
   
      SUBROUTINE CORRWM
C     =================
C
C     The routine for management of various flags for treating 
C     frequency points; in particular those connected to the so-called 
C     "subtraction weights" (in the non-overlapping mode only)
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'ATOMIC.FOR'
      INCLUDE 'MODELQ.FOR'
      PARAMETER (T15=1.D-15)
C
      NFREQE=0
      DO 10 IJ=1,NFREQ
         IJEX(IJ)=0
         DO ID=1,ND
            LSKIP(ID,IJ)=.FALSE.
         END DO
         if(ifprad.eq.0) then
            do id=1,nd
               lskip(id,ij)=.true.
            end do
         end if
         IF(IJALI(IJ).NE.0) GO TO 10
         NFREQE=NFREQE+1
         IJEX(IJ)=NFREQE
         IJFR(NFREQE)=IJ
   10 CONTINUE
c
      if(ifryb.ne.0) then
         nfreqe=0
         do ij=1,nfreq
            ijex(ij)=0
         end do
      end if
C
      IF(IBFINT.LE.0) THEN
         DO IJ=1,NFREQ
            IJBF(IJ)=IJ
            AIJBF(IJ)=UN
         END DO
       ELSE
         IF(ISPODF.EQ.0) THEN
            DO IJ=1,NFREQC
               IJBF(IJ)=IJ
               AIJBF(IJ)=UN
            END DO
            IF(NFREQ.GT.NFREQC) THEN
               DO IJ=NFREQC+1,NFREQ
                  FR=FREQ(IJ)
                  IJ0=1
                  DO IJT=1,NFREQC
                     IF(FREQ(IJT).LE.FR) THEN
                     IJ0=IJT
                     GO TO 12
                     END IF
                  END DO
   12             IJ1=IJ0-1
                  A1=(FR-FREQ(IJ0))/(FREQ(IJ1)-FREQ(IJ0))
                  IJBF(IJ)=IJ1
                  AIJBF(IJ)=A1
               END DO
            END IF
          ELSE
            DO IJ=1,NFREQC-1
               IJ0=IFREQB(IJ)
               IJ1=IFREQB(IJ+1)
               DO KJ=IJ0,IJ1-1
                  IJBF(KJ)=IJ
                  AIJBF(KJ)=(FREQ(KJ)-FREQ(IJ1))/(FREQ(IJ0)-FREQ(IJ1))
               END DO
            END DO
            IJ0=IFREQB(NFREQC)
            IJBF(IJ0)=NFREQC
            AIJBF(IJ0)=UN
         END IF
      END IF
C
      if(nfreqe.gt.mfrex) CALL QUIT('nfreqe.gt.mfrex',nfreqe,mfrex)
C
      DO 100 ITR=1,NTRANS
         IF(.NOT.LINE(ITR)) GO TO 100
C
C        first set up array LSKIP(ID,IJ), which has values
C        TRUE - if the radiation at frequency point IJ does not contribute
C               radiation pressure (ie. this point belongs to a transition
C               for which the user required the radiation pressure to be
C               skipped - IABS(INDEXP) chosen as 9 or 19)
C        FALSE - normal calculation of radiation pressure
C
         INX=IABS(INDEXP(ITR))
         IF(INX.EQ.9.OR.INX.GE.19) THEN
            DO IJ=IFR0(ITR),IFR1(ITR)
               DO ID=1,ND
                  LSKIP(ID,IJ)=.TRUE.
               END DO
            END DO
         END IF
  100 CONTINUE
C   
c     WRITE(6,609)
      DO 110 IJ=1,NFREQ
         FR15=FREQ(IJ)*T15
         W0E(IJ)=W(IJ)*PI4H/FREQ(IJ)
         BNUE(IJ)=BN*FR15*FR15*FR15
         IF(IJALI(IJ).NE.0.or.ifryb.gt.0) GO TO 110
c        if(ispodf.eq.0) then
c           WRITE(6,610) IJ,FREQ(IJ),W(IJ),PROF(IJ)
c         else
c           WRITE(6,610) IJ,FREQ(IJ),W(IJ)
c        end if
  110 CONTINUE
C
      DO IJ=1,NFREQ
         WC(IJ)=W(IJ)
         IF(IJALI(IJ).LE.0) WC(IJ)=0.
      END DO
C
c 609 FORMAT(1H0//' FREQUENCY POINTS AND WEIGHTS - EXPLICIT'/
c    *            ' ---------------------------------------'//
c    * '       IJ',7X,'FREQ',13X,'WEIGHT',11X,'PROF'/)
c 610 FORMAT(1H ,I8,1P2D17.8,D15.5,D17.8)
      RETURN 
      END
C
C
C     *******************************************************************
C
C
   
      SUBROUTINE IJALIS(ITR,IFRQ0,IFRQ1)
C     ==================================
C
C     auxiliary routine - sets up the necessary flags for ALI treatment
C      of individual transitions (in the fully hybrid CL/ALI scheme)
C
C     Output:
C
C     IJALI(IJ) = 0 - frequency point IJ is an explicit point
C               = 1 - frequency point IJ is an ALI point
C
C     LEXP(ITR) = T - at least one point within transition ITR is explicit
C     LEXP(ITR) = F - no point within transition ITR is explicit
C     LALI(ITR) = T - at least one point within transition ITR is ALI
C     LALI(ITR) = F - no point within transition ITR is ALI
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'ATOMIC.FOR'
      INCLUDE 'MODELQ.FOR'
C
      INDXP=INDEXP(ITR)
      I0=IFR0(ITR)
      I1=IFR1(ITR)
      NF=I1-I0+1
      DO IJ=I0,I1
         IF(INDXP.GT.0) THEN
            IJALI(IJ)=0
          ELSE IF(INDXP.LT.0) THEN
            IJALI(IJ)=1
         END IF
      END DO
C
C     primarily explicit transitions
C
      IF(INDXP.GT.0) THEN
         LEXP(ITR)=.TRUE.
         LALI(ITR)=.FALSE.
         IF(IFRQ0.GT.0) THEN
            LALI(ITR)=.TRUE.
            IF(IFRQ1.EQ.0.OR.IFRQ1.GT.NF) IFRQ1=NF
            DO I=IFRQ0,IFRQ1
                IJALI(I0+I-1)=1
            END DO
          ELSE IF(IFRQ0.LT.0) THEN
            LALI(ITR)=.TRUE.
            READ(57,*) (IJALI(IJ),IJ=I0,I1)
         END IF
         IF(IFRQ0.EQ.1.AND.IFRQ1.EQ.NF) LEXP(ITR)=.FALSE.
       ELSE IF(INDXP.LT.0) THEN
C
C      primarily ALI transitions
C
         LALI(ITR)=.TRUE.
         LEXP(ITR)=.FALSE.
         IF(IFRQ0.GT.0) THEN
            LEXP(ITR)=.TRUE.
            IF(IFRQ1.EQ.0.OR.IFRQ1.GT.NF) IFRQ1=NF
            DO I=IFRQ0,IFRQ1
                IJALI(I0+I-1)=0
            END DO
          ELSE IF(IFRQ0.LT.0) THEN
            LEXP(ITR)=.TRUE.
            READ(57,*) (IJALI(IJ),IJ=I0,I1)
         END IF
         IF(IFRQ0.EQ.1.AND.IFRQ1.EQ.NF) LALI(ITR)=.FALSE.
      END IF
      IF(NFFIX.GT.0) THEN
         DO IJ=I0,I1
            IJALI(IJ)=1
         END DO
         LALI(ITR)=.TRUE.
         LEXP(ITR)=.FALSE.
      END IF
      RETURN
      END
C
C
C     *******************************************************************
C
C
   
      SUBROUTINE IJALI2
C     =================
C
C     auxiliary routine - sets up the necessary flags for ALI treatment
C     of individual transitions (in the fully hybrid CL/ALI scheme)
C
C     Version for opacity sampling mode
C
C     Output:
C
C     IJALI(IJ) = 0 - frequency point IJ is an explicit point
C               = 1 - frequency point IJ is an ALI point
C
C     LEXP(ITR) = T - at least one point within transition ITR is explicit
C     LEXP(ITR) = F - no point within transition ITR is explicit
C     LALI(ITR) = T - at least one point within transition ITR is ALI
C     LALI(ITR) = F - no point within transition ITR is ALI
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'ATOMIC.FOR'
      INCLUDE 'MODELQ.FOR'
      INCLUDE 'ODFPAR.FOR'
C
      DO IJ=1,NFREQ
         IJALI(IJ)=1
         IJX(IJ)=1
         NLINES(IJ)=0
      END DO
      NLITOT=0
      NLIMAX=0
C
C     Overlapping lines at frequency IJ
C
      DO 10 IT=1,NTRANS
         IF(LINEXP(IT)) GO TO 10
         DO IJ=IFR0(IT),IFR1(IT)
            NLINES(IJ)=NLINES(IJ)+1
            ITRLIN(NLINES(IJ),IJ)=int2(IT)
         END DO
   10 CONTINUE
      DO IJ=1,NFREQ
         NLITOT=NLITOT+NLINES(IJ)
         IF(NLINES(IJ).GT.MITJ)
     *   CALL QUIT('Too many overlappins-nlines(ij).gt.mitj',
     *   nlines(ij),mitj)
         IF(NLINES(IJ).GT.NLIMAX) NLIMAX=NLINES(IJ)
      END DO
      WRITE(10,*) ' Max. number of line overlaps:       ',NLIMAX
      WRITE(10,*) ' Total number of line overlaps:      ',NLITOT
C
C     Switches for ALI and explicit transitions
C
      IF(NFFIX.EQ.2) THEN
         DO ITR=1,NTRANS
            LEXP(ITR)=.FALSE.
            LALI(ITR)=.TRUE.
         END DO
         RETURN
      END IF
C
      XFRMA=DLOG10(FRS1)
      DO 100 ITR=1,NTRANS
         INDXP=INDEXP(ITR)
         I0=IFR0(ITR)
         I1=IFR1(ITR)
         NF=I1-I0+1
         IF(FR0(ITR).GT.FRS1) GO TO 100
         IJL=IJTC(ITR)
C
C        primarily explicit line transitions
C
         IF(LINE(ITR)) THEN
            IF(INDXP.GT.0) THEN
               LEXP(ITR)=.TRUE.
               LALI(ITR)=.FALSE.
               IF(IFC0(ITR).EQ.0) THEN
                  DO IJ=I0,I1
                     IJALI(IJ)=0
                  END DO
                ELSE
                   LALI(ITR)=.TRUE.
                   NFC=IABS(IFC1(ITR)-IFC0(ITR)+1)
                   IF(NFC.EQ.NF) THEN
                      LEXP(ITR)=.FALSE.
                    ELSE
                      NFC=NFC/2
                      DO IJ=I0,IJL-NFC
                         IJALI(IJ)=0
                      END DO
                      DO IJ=IJL+NFC,I1
                         IJALI(IJ)=0
                      END DO
                   END IF
               END IF
             ELSE IF(INDXP.LT.0) THEN
C
C            primarily ALI line transitions
C
               LEXP(ITR)=.FALSE.
               LALI(ITR)=.TRUE.
               IF(IFC0(ITR).NE.0) THEN
                  LEXP(ITR)=.TRUE.
                  NFC=IABS(IFC1(ITR)-IFC0(ITR)+1)
                  IF(NFC.EQ.NF) THEN
                     LALI(ITR)=.FALSE.
                     DO IJ=I0,I1
                        IJALI(IJ)=0
                     END DO
                   ELSE
                     NFC=NFC/2
                     DO IJ=IJL-NFC,IJL+NFC
                        IJALI(IJ)=0
                     END DO
                 END IF
               END IF
            END IF
C
C         continuum transitions
C
          ELSE
            IF(IFC0(ITR).GT.0) THEN
               DO IJ=1,IFC1(ITR)-IFC0(ITR)+1
                  IJALI(IJL-IJ+1)=0
               END DO
            END IF
         END IF
  100 CONTINUE
C
      RETURN
      END
C
C
C     *******************************************************************
C
C
   
      SUBROUTINE LEVSET
C     =================
C
C     sets up level parameters IIEXP and IIFOR which control the
C     treatment of levels
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'ATOMIC.FOR'
      INCLUDE 'MODELQ.FOR'
      if(ioptab.lt.0) return
C
C     1. case - treatment of levels determined by IMODL
C
      IF(IFLEV.EQ.0) THEN
         DO I=1,NLEVEL
            IIEXP(I)=0
            IIFOR(I)=0
         END DO
         IIE=0
         IIF=0
         DO 20 IAT=1,NATOM
            IGRP=0
            IF(IIFIX(IAT).EQ.1) GO TO 20
            DO I=N0A(IAT),NKA(IAT)
               INEW=1
               IF(IMODL(I).EQ.0) THEN
                  IIE=IIE+1
                  IIF=IIF+1
                  IIEXP(I)=IIE
                  IIFOR(I)=IIF
                  INDLEV(IIEXP(I))=I
                ELSE IF(IMODL(I).GT.0) THEN
                  IIF=IIF+1
                  IIFOR(I)=IIF
                  IF(ILTLEV(I).GE.1) THEN
                     IIE=IIE+1
                     IIEXP(I)=IIE
                  END IF
                  IF(I.EQ.NFIRST(IEL(I)).OR.I.EQ.NNEXT(IEL(I))) THEN
                     IIE=IIE+1
                     IIEXP(I)=IIE
                  END IF
                ELSE IF(IMODL(I).LT.-100) THEN
                  IF(I.GT.1) THEN
                     IF(IMODL(I).EQ.IMODL(I-1)) INEW=0
                  END IF
                  IIEXP(I)=-IIE
                  IF(INEW.EQ.1) THEN
                     IIE=IIE+1
                     IIEXP(I)=-IIE
                     IM=NFIRST(IEL(I))
                     LML=.TRUE.
                     DO WHILE (IM.LT.I-1 .AND. LML)
                        IF(IMODL(I).EQ.IMODL(IM)) THEN
                           IIEXP(I)=IIEXP(IM)
                           IIE=IIE-1
                           LML=.FALSE.
                        END IF
                        IM=IM+1
                     END DO
                  END IF
                  IGRP=1
                  IIF=IIF+1
                  IIFOR(I)=IIF
                ELSE IF(IMODL(I).LT.-200) THEN
                  IF(I.GT.1) THEN
                     IF(IMODL(I).EQ.IMODL(I-1)) INEW=0
                  END IF
                  IF(INEW.EQ.1) IIE=IIE+1
                  IF(INEW.EQ.1) IIF=IIF+1
                  IIEXP(I)=-IIE
                  IIFOR(I)=-IIF
                  IGRP=1
               END IF
            END DO
            IF(IGRP.EQ.1) THEN
               DO I=N0A(IAT),NKA(IAT)
                  IF(IIEXP(I).GT.0) IIEXP(I)=-IIEXP(I)
                  IF(IMODL(I).EQ.0) IMODL(I)=7
               END DO
            END IF
   20    CONTINUE
         NLVEXP=IABS(IIE)
         if(nlvexp.gt.mlvexp) 
     *   CALL QUIT('nlvexp.gt.mlvexp',nlvexp,mlvexp)
         NLVFOR=IABS(IIF)
         DO 30 I=1,NLEVEL
            IF(IMODL(I).EQ.1.OR.IMODL(I).EQ.3) THEN
               IIEXP(I)=0
             ELSE IF(IMODL(I).EQ.4.OR.IMODL(I).EQ.5) THEN
               IIEXP(I)=0
             ELSE IF(IMODL(I).EQ.-1.OR.IMODL(I).EQ.-3) THEN
               IIEXP(I)=0
               IIFOR(I)=0
             ELSE IF(IMODL(I).EQ.6) THEN
               IIEXP(I)=0
             ELSE IF(IMODL(I).EQ.-5.OR.IMODL(I).EQ.-6) THEN
               IIEXP(I)=0
               IIFOR(I)=0
             ELSE IF(IMODL(I).LT.-100) THEN
               IMODL(I)=7
             ELSE IF(IMODL(I).LT.-200) THEN
                IMODL(I)=-7
            END IF
            DO ID=1,ND
               ILTREF(I,ID)=NNEXT(IEL(I))
            END DO      
   30    CONTINUE
         IF(IGRP.EQ.1) THEN
            DO I=N0A(IAT),NKA(IAT)
               IMODL(I)=7
            END DO
         END IF
C
C     2. case - treatment of levels automatic - all levels with ILK=0
C               in updated LTE mode
C
       ELSE
         IIF=0
         DO 110 I=1,NLEVEL
            IF(IIFIX(IATM(I)).EQ.1) GO TO 110
            IMODL(I)=5
            IF(I.EQ.NFIRST(IEL(I)).OR.I.EQ.NNEXT(IEL(I))) THEN
               IIF=IIF+1
               IIFOR(I)=IIF
            END IF
  110    CONTINUE
         NLVEXP=IIF
         if(nlvexp.gt.mlvexp) 
     *   CALL QUIT('nlvexp.gt.mlvexp',nlvexp,mlvexp)
         NLVFOR=IIF
         DO I=1,NLEVEL
            IF(I.NE.NFIRST(IEL(I)).AND.I.NE.NNEXT(IEL(I)))
     *      IIFOR(I)=0
         END DO
         DO I=1,NLEVEL
            IIEXP(I)=IIFOR(I)
            IF(IIEXP(I).GT.0) INDLEV(IIEXP(I))=I
            DO ID=1,ND
               ILTREF(I,ID)=NNEXT(IEL(I))
            END DO
         END DO
         if(.not.lte) then
            do i=1,nlevel
               iifor(i)=i
            end do
           nlvfor=nlevel
         end if
c      
      END IF
C     
C     initialize b-factors
C
      DO I=1,NLEVEL
         DO ID=1,ND
            BFAC(I,ID)=UN
         END DO
      END DO
C     
      do ii=1,nlvexp
         indlev(ii)=0
         do id=1,nd
            igzero(ii,id)=0
         end do
      end do
      do i=1,nlevel
         do id=1,nd
            ipzero(i,id)=0
         end do
         if(iabs(imodl(i)).le.6) then
            IF(IIEXP(I).GT.0) INDLEV(IIEXP(I))=I
         end if
      end do
C 
      RETURN
      END
C
C
C     *******************************************************************
C
C

      SUBROUTINE DWNFR0(ID)
C     =====================
C
C     Auxiliary quantities for dissolved fractions
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'MODELQ.FOR'
      PARAMETER (SIXTH=UN/6.,CCOR=0.09)
      parameter (p1=0.1402,p2=0.1285,p3=un,p4=3.15,p5=4.)
      parameter (f23=-2./3.)
C
      ANE=ELEC(ID)
      ELEC23(ID)=EXP(F23*LOG(ANE))
      ANES=EXP(SIXTH*LOG(ANE))
      ACOR(ID)=CCOR*ANES/SQT1(ID)
      X=EXP(P4*LOG(UN+P3*ACOR(ID)))
      DWC2(ID)=P2*X
      A3=ACOR(ID)*ACOR(ID)*ACOR(ID)
      DO 10 IZZ=1,MZZ
         Z3(IZZ)=IZZ*IZZ*IZZ
         DWC1(IZZ,ID)=P1*(X+P5*(IZZ-UN)*A3)
   10 CONTINUE
      RETURN
      END
C
C
C ********************************************************************
C
C
      SUBROUTINE DWNFR1(FR,FR0,ID,IZZ,DW1)
C     ====================================
C
C     dissolved fraction for frequency FR
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'MODELQ.FOR'
      PARAMETER (TKN=3.01,CKN=5.33333333,CB0=8.59d14)
      PARAMETER (SQFRH=5.734152D7)
C
      cb=cb0*bergfc
c
      IF(FR.LT.FR0) THEN
         XN=SQFRH*IZZ/SQRT(FR0-FR)
         if(xn.le.tkn) then
            xkn=un
          else
            xn1=un/(xn+un)
            xkn=ckn*xn*xn1*xn1
         end if
         BETA=CB*Z3(IZZ)*XKN/(XN*XN*XN*XN)*ELEC23(ID)
         BETA3=BETA*BETA*BETA
         BETA32=SQRT(BETA3)
         F=(DWC1(IZZ,ID)*BETA3)/(UN+DWC2(ID)*BETA32)
         DW1=UN-F/(UN+F)
       ELSE
         DW1=UN
      END IF
      RETURN
      END
C
C
C     ****************************************************************
C
C
      SUBROUTINE SGMER0
C     =================
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'ATOMIC.FOR'
      INCLUDE 'MODELQ.FOR'
      PARAMETER (FRH=3.28805E15, PH2=2.815D29*2., EHB=157802.77355)
      DIMENSION FREDG(NLMX),S(NLMX),SUM(NLMX),SUD(NLMX)
C
      IMER=0
      DO 100 II=1,NLEVEL
         IF(IFWOP(II).GE.0) GO TO 100
         IMER=IMER+1
         IMRG(II)=IMER
         IIMER(IMER)=II
         IE=IEL(II)
         CH=IZ(IE)*IZ(IE)
         FRCH(IMER)=FRH*CH
         SGM0(IMER)=PH2*CH*CH
         II0=NQUANT(II-1)+1
         DO 90 ID=1,ND
            EX=EHB*CH*TEMP1(ID)
            DO 10 I=II0,NLMX
               FREDG(I)=FRCH(IMER)*XI2(I)
               EXI=EXP(EX*XI2(I))
               S(I)=EXI*WNHINT(I,ID)*XI3(I)
               SUM(I)=0.
   10       CONTINUE
            SUM(NLMX)=S(NLMX)
            SUD(NLMX)=S(NLMX)*XI2(NLMX)
            DO 20 I=NLMX-1,II0,-1
               SUM(I)=SUM(I+1)+S(I)
   20       CONTINUE
            DO 30 I=1,II0-1
               SUM(I)=SUM(II0)
   30       CONTINUE
            SGEM=SGM0(IMER)/GMER(IMER,ID)
            DO 60 I=1,NLMX
               SGMSUM(I,IMER,ID)=SUM(I)*SGEM
   60       CONTINUE
   90    CONTINUE
  100 CONTINUE
      RETURN
      END
C
C
C     ****************************************************************
C
C
      SUBROUTINE SGMER1(FRINV,FR3INV,IMER,ID,SGME1)
C     =============================================
C
C     photoionization cross-section for a merged level
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'ATOMIC.FOR'
      INCLUDE 'MODELQ.FOR'
C
      ISU=INT(SQRT(FRCH(IMER)*FRINV))+1
      SGME1=SGMSUM(ISU,IMER,ID)*FR3INV
      RETURN
      END
C
C
C     ****************************************************************
C
C
      SUBROUTINE SGMERD(FRINV,FR3INV,IMER,ID,SGME1,DSGME1)
C     ====================================================
C
C     photoionization cross-section for a merged level
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'ATOMIC.FOR'
      INCLUDE 'MODELQ.FOR'
C
      ISU=INT(SQRT(FRCH(IMER)*FRINV))+1
      SGME1=SGMSUM(ISU,IMER,ID)*FR3INV
      DSGME1=-SGMSUD(ISU,IMER,ID)*FR3INV
      RETURN
      END
C
C
C     ****************************************************************
C
C
      SUBROUTINE TDPINI
C     =================
C
C     initialization of only temperature dependent quantities
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'ATOMIC.FOR'
      INCLUDE 'MODELQ.FOR'
      INCLUDE 'ODFPAR.FOR'
      INCLUDE 'ALIPAR.FOR'
      PARAMETER (CFF1=1.3727D-25,CFF2=4.3748D-10,CFF3=2.5993D-7)
      PARAMETER (SGFF0 = 3.694D8)
C
C     temperature-dependent quantities
C
      DO 10 ID=1,ND
         T=TEMP(ID)
         T1=UN/T
         HKT1(ID)=HK*T1
         HKT21(ID)=HKT1(ID)*T1
         TK1(ID)=HKT1(ID)/H
         SQT1(ID)=SQRT(T)
         TEMP1(ID)=T1
         CALL GFREE0(ID)
         EMEL1(ID)=UN
   10 CONTINUE
C
C     delta m (for evaluation of optical depths)
C
      DO 20 ID=1,ND-1
         DELDM(ID)=HALF*(DM(ID+1)-DM(ID))
         deldmz(id)=deldm(id)
         if(izscal.eq.1) deldmz(id)=half*(zd(id)-zd(id+1))
   20 CONTINUE
      DEDM1=DM(1)/DENS(1)
      RETURN
      END
C
C
C     ****************************************************************
C
C
C
      SUBROUTINE OPAINI(IMOD)
C     =======================
C
C     initialization of only depth-dependent quantities
C     for evaluation of opacities
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'ATOMIC.FOR'
      INCLUDE 'MODELQ.FOR'
      INCLUDE 'ODFPAR.FOR'
      INCLUDE 'ALIPAR.FOR'
      DIMENSION PRF(MFREQL)
      PARAMETER (CFF1=1.3727D-25,CFF2=4.3748D-10,CFF3=2.5993D-7)
      PARAMETER (SIXTH=UN/6.,CCOR=0.09,T32=1.5D0)
      PARAMETER (SGFF0 = 3.694D8)
      DATA ICOMP /0/
C
      imod1=imod
      DO 100 ID=1,ND
         WMT=WMM(ID)*YTOT(ID)
         T=TEMP(ID)
         ANE=ELEC(ID)
         ELEC1(ID)=UN/ANE
         DENS1(ID)=UN/DENS(ID)
         DENSI(ID)=DENS1(ID)
         DENSIM(ID)=DENSI(ID)*WMM(ID)
         ELSCAT(ID)=ANE*SIGE
         CALL DWNFR0(ID)
         CALL WNSTOR(ID)
c        CALL SABOLF(ID)
c        CALL REFLEV(ID,IMOD)
c        CALL LEVGRP(ID,IIEXP,0,POPP)
         DO II=1,NLEVEL
            POPINV(II,ID)=0.
            IF(POPUL(II,ID).NE.0.) POPINV(II,ID)=UN/POPUL(II,ID)
         END DO
         DO II=1,NLEVEL
            IIE=IIEXP(II)
            IF(IIE.EQ.0) THEN
              IE=ILTREF(II,ID)
               PP(II,ID)=POPUL(II,ID)*POPINV(IE,ID)
c              IF(IABS(IMODL(II)).LE.5) THEN
c                 PT(II,ID)=POPUL(II,ID)*DSBPST(II,ID)
c                PN(II,ID)=POPUL(II,ID)*DSBPSN(II,ID)
c              END IF
c            ELSE IF(IIE.LT.0) THEN
c              PP(II,ID)=SBPSI(II,ID)
            END IF
         END DO
         DO ION=1,NION
            USUMS(ION,ID)=USUM(ION)
            DUSMT(ION,ID)=DUSUMT(ION)
            DUSMN(ION,ID)=DUSUMN(ION)
         ENDDO
c
c        quantities for the bound-free opacity
c
         DO 10 IBFT=1,NTRANC
            ITR=ITRBF(IBFT)
            II=ILOW(ITR)
            JJ=IUP(ITR)
            IT=ITRA(JJ,II)
            IE=IEL(II)
            NKE=NNEXT(IE)
            CORR=UN
            IF(NKE.NE.JJ) CORR=G(NKE)/G(JJ)*
     *         EXP((ENION(NKE)-ENION(JJ))*TK1(ID))
            ABTRA(ITR,ID)=POPUL(II,ID)
            EMTRA(ITR,ID)=POPUL(JJ,ID)*ANE*SBF(II)*WOP(II,ID)*CORR
            DEMLT(ITR,ID)=-(T32+FR0(ITR)*HKT1(ID))/TEMP(ID)
   10    CONTINUE
c
c        quantities for the free-free opacity
c
         IF(IELHM.GT.0) THEN
            CFFN(ID)=POPUL(NFIRST(IELH),ID)*ANE
            CFFT(ID)=CFF2-CFF3/T
         END IF
         SGFF=SGFF0/SQT1(ID)*ANE
         DO 20 ION=1,NION
            SFF2(ION,ID)=EXP(FF(ION)*HKT1(ID))
            SFF3(ION,ID)=POPUL(NNEXT(ION),ID)*CHARG2(ION)*SGFF
            DSFF(ION,ID)=(FF(ION)*HKT1(ID)+HALF)*TEMP1(ID)
   20    CONTINUE
  100 CONTINUE
      if(izscal.eq.1) then
         do id=1,nd
            densi(id)=un
            densim(id)=0.
         end do
      end if
      CALL SGMER0
C
C     initialization of the line opacity
C
      LASER=ITER.GT.ITLAS
      DO 200 ITR=1,NTRANS
         INDXA=IABS(INDEXP(ITR))
         IF(.NOT.LINE(ITR)) GO TO 200
         II=ILOW(ITR)
         JJ=IUP(ITR)
         IF(INDXA.NE.0.AND.INTMOD(ITR).NE.0 .AND. ICOMP.EQ.0) THEN
           IJL0=IFR0(ITR)
           IJL1=IFR1(ITR)
         IF(ISPODF.GE.1) THEN
             IJL0=KFR0(ITR)
             IJL1=KFR1(ITR)
         END IF
           IF(INDXA.LT.2.OR.INDXA.GT.4) THEN
             DO 150 ID=1,ND
               CALL LINPRO(ITR,ID,PRF)
               DO 110 IJ=IJL0,IJL1
                  PRFLIN(ID,IJ)=real(PRF(IJ-IJL0+1))
  110          CONTINUE
  150        CONTINUE
           ENDIF
         END IF
         GG=G(II)/G(JJ)
         DO 160 ID=1,ND
            IF(IFWOP(JJ).GE.0) THEN
               PI=POPUL(II,ID)*WOP(JJ,ID)
               PJ=GG*POPUL(JJ,ID)*WOP(II,ID)
             ELSE
               PI=POPUL(II,ID)
               PJ=G(II)/GMER(IMRG(JJ),ID)*POPUL(JJ,ID)*WOP(II,ID)
            END IF
            ABTRA(ITR,ID)=PI
            EMTRA(ITR,ID)=PJ*EXP(FR0(ITR)*HKT1(ID))
            DEMLT(ITR,ID)=-FR0(ITR)*HKT21(ID)
            IF(LASER) THEN
            qtt=0.
            if(pi.ne.pj) QTT=PJ/(PI-PJ)*(EXP(FR0(ITR)*HKT1(ID))-UN)
            IF(QTT.LT.0. .OR. QTT.GT.QTLAS) THEN
               ABTRA(ITR,ID)=0.
               EMTRA(ITR,ID)=0.
               DEMLT(ITR,ID)=0.
            END IF
            END IF
c
c           set up ABTRA and EMTRA to 0 in the range where
c           the hydrogen opacity is calculated from Gomez tables
c
            if(ihgom.gt.0.and.elec(id).gt.hglim) then
               if(ii.ge.n0hn.and.ii.le.n0hn-1+ihgom) then
                  abtra(itr,id)=0.
                  emtra(itr,id)=0.
                  demlt(itr,id)=0.
               end if
            end if
c
  160    CONTINUE
  200 CONTINUE
      ICOMP=1
      RETURN
      END
C
C
C     ****************************************************************
C
C

      SUBROUTINE TRAINI
C     =================
C
C     initialization of depth-independent quantities
C     for evaluation of opacities
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'ATOMIC.FOR'
      INCLUDE 'MODELQ.FOR'
      INCLUDE 'ODFPAR.FOR'
c 
      do itr=1,ntrans
         idiel(itr)=0
      end do
C
C     bound-free transitions
C
      NCDW=0
      DO 10 IBFT=1,NTRANC
         ITR=ITRBF(IBFT)
         ii=ilow(itr)
         if(ilk(iup(itr)).ne.0.and.nfirst(iel(ii)).eq.ii.
     *      and.IFDIEL.NE.0) idiel(itr)=1
         MODW=IABS(INDEXP(ITR))
         IF(MODW.NE.5.AND.MODW.NE.15) GO TO 10
         NCDW=NCDW+1
         MCDW(ITR)=NCDW
         ITRCDW(NCDW)=ITR
   10 CONTINUE
      IF(ISPODF.GE.1) RETURN
C
C     bound-bound transitions
C
      DO IJ=1,NFREQ
         NLINES(IJ)=0
      END DO
C
      DO 100 ITR=1,NTRANS
         IF(LINEXP(ITR)) GOTO 100
         DO IJ=IFR0(ITR),IFR1(ITR)
            IJLIN(IJ)=ITR
         END DO
  100 CONTINUE
      RETURN
      END

C
C

C     ****************************************************************
C
C
      SUBROUTINE RTEDF1(IJ)
C     =====================
C
C     Solution of the radiative transfer equation - for one frequency 
C     for the known source function.
C     Determination of the radiation field and variable Eddington
C     factors.
C
C     The numerical method used:
c     Discontinuous Finite Element method
c     Castor, Dykema, Klein, 1992, ApJ 387, 561.
C
C     different formulation of the boundary conditions
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'MODELQ.FOR'
      INCLUDE 'ALIPAR.FOR'
      COMMON/OPTDPT/DT(MDEPTH)
      PARAMETER (SIXTH=UN/6.D0,
     *           THIRD=UN/3.D0,
     *           TWOTHR=TWO/3.D0)
      DIMENSION RDK(MDEPTH),FKK(MDEPTH),
     *          ST0(MDEPTH),SA0(MDEPTH),SS0(MDEPTH),
     *          dtau(mdepth),rip(mdepth),rim(mdepth),
     *          riin(mdepth),riup(mdepth),u(mdepth),
     *          aip(mdepth),aim(mdepth),al0(mdepth),
     *          aiin(mdepth),aiup(mdepth),
     *          ali0(mdepth),ss0c(mdepth),
     *          AAA(MDEPTH),BBB(MDEPTH),CCC(MDEPTH),EEE(MDEPTH),
     *          ZZZ(MDEPTH),ALRH(MDEPTH),ALRM(MDEPTH),ALRP(MDEPTH),
     *          DDD(MDEPTH),AANU(MDEPTH)
C
      FR=FREQ(IJ)
C
C     optical depth scale
C
      DO ID=1,ND-1
         DT(ID)=DELDMZ(ID)*(ABSOT(ID+1)+ABSOT(ID))
         SA0(ID)=EMIS1(ID)/ABSO1(ID)
         SS0(ID)=-SCAT1(ID)/ABSO1(ID)
      END DO
      SA0(ND)=EMIS1(ND)/ABSO1(ND)
      SS0(ND)=-SCAT1(ND)/ABSO1(ND)
C
      TAUMIN=ABSO1(1)*DEDM1
C
C     Allowance for wind blanketing
C
      ALB1=0.
      IF(IWINBL.GT.0) ALB1=TWO*ALBE(IJ)/(UN+ALBE(IJ))
C
C     quantities for the lower boundary condition
C
      FR15=FR*1.D-15
      BNU=BN*FR15*FR15*FR15
      PLAND=BNU/(EXP(HK*FR/TEMP(ND))-UN)*RRDIL
      DPLAN=BNU/(EXP(HK*FR/TEMP(ND-1))-UN)*RRDIL
      IF(TEMPBD.GT.0.) THEN
        PLAND=BNU/(EXP(HK*FR/TEMPBD)-UN)*RRDIL
        DPLAN=BNU/(EXP(HK*FR/TEMPBD)-UN)*RRDIL
      ENDIF
      DPLAN=(PLAND-DPLAN)/DT(ND-1)
C
C     global ALI loop for treating electron scattering
C
      itrali=0
   10 itrali=itrali+1
C
C     total source function
C  
      DO 20 ID=1,ND
         ST0(ID)=SA0(ID)-SS0(ID)*RAD(IJ,ID)
         RAD1(ID)=0.
         RDK(ID)=0.
         ALI1(ID)=0.
   20 CONTINUE
      AH=0.
      ahout=0.
      ahd=0.
      U0=0.
      QQ0=0.
      US0=0.
c
c     loop over angle poits
c
      DO 100 I=1,NMU
         AMU2=AMU(I)*AMU(I)*WTMU(I)
         do id=1,nd-1
            dtau(id)=dt(id)/amu(i)
         enddo
c               
c        incoming intensity 
c   
         ID=1
         P0=0.
         EX=UN
C
C        allowance for non-zero optical depth at the first depth point
C
c        rim(id)=EXTRAD(IJ)
         rim(id)=EXTINT(IJ,I)
c         IF(IWINBL.EQ.0) THEN
c            TAMM=TAUMIN/AMU(I)
c            EX=EXP(-TAMM)
c            P0=UN-EX
c            QQ0=QQ0+P0*AMU(I)*WTMU(I)
c            U0=U0+EX*WTMU(I)
c            US0=US0+P0/TAMM*WTMU(I)
c            rim(id)=st0(id)*p0
c         END IF
c               
c          incoming intensity 
c   
            aim(id)=0.
            do id=1,nd-1
               dt0=dtau(id)
               dtaup1=dt0+un
               dtau2=dt0*dt0
               bb=two*dtaup1
               cc=dt0*dtaup1
               aa=un/(dtau2+bb)
               rim(id+1)=(two*rim(id)+dt0*st0(id)+cc*st0(id+1))*aa
               rip(id)=(bb*rim(id)+cc*st0(id)-dt0*st0(id+1))*aa
               aim(id+1)=cc*aa
               aip(id)=(cc+bb*aim(id))*aa
            enddo
            do id=2,nd-1
               dtt=un/(dtau(id-1)+dtau(id))
               riin(id)=(rim(id)*dtau(id)+rip(id)*dtau(id-1))*dtt
               aiin(id)=(aim(id)*dtau(id)+aip(id)*dtau(id-1))*dtt
            enddo
            riin(1)=rim(1)
            riin(nd)=rim(nd)
            aiin(1)=aim(1)
            aiin(nd)=aim(nd)
C
c           outgoing intensity
c 
            if(idisk.eq.0) rim(nd)=PLAND+AMU(I)*DPLAN
            do id=nd-1,1,-1
               dt0=dtau(id)
               dtaup1=dt0+un
               dtau2=dt0*dt0
               bb=two*dtaup1
               cc=dt0*dtaup1
               aa=un/(dtau2+bb)
               rim(id)=(two*rim(id+1)+dt0*st0(id+1)+cc*st0(id))*aa
               rip(id+1)=(bb*rim(id+1)+cc*st0(id+1)-dt0*st0(id))*aa
               aim(id)=cc*aa
               aip(id+1)=(cc+bb*aim(id+1))*aa
            enddo
            do id=2,nd-1
               dtt=un/(dtau(id-1)+dtau(id))
               riup(id)=(rim(id)*dtau(id-1)+rip(id)*dtau(id))*dtt
               aiup(id)=(aim(id)*dtau(id-1)+aip(id)*dtau(id))*dtt
            enddo
            riup(1)=rim(1)
            riup(nd)=rim(nd)
            aiup(1)=aim(1)
            aiup(nd)=aim(nd)
c
c           final symmetrized (Feautrier) intensity -- (riin+riup)/2
c
            do id=1,nd
               u(id)=(riin(id)+riup(id))*half
               al0(id)=(aiin(id)+aiup(id))*half
            enddo
c
         DO 70 ID=1,ND
            RAD1(ID)=RAD1(ID)+WTMU(I)*U(ID)
            RDK(ID)=RDK(ID)+AMU2*U(ID)
            ALI1(ID)=ALI1(ID)+WTMU(I)*AL0(ID)
   70    CONTINUE
         AH=AH+AMU(I)*WTMU(I)*U(1)
         ahd=ahd+amu(i)*wtmu(i)*u(nd)
         ahout=ahout+amu(i)*wtmu(i)*riup(1)
c
c     end of the loop over angle points
c
  100 CONTINUE
C
C     solution of the transfer equation
C     Variables:
C     RAD1    - mean intensity
C     FAK1    - Eddington factor f(K) = K/J
C     FH      - the "surface" Eddington factor
C     ALI1    - diagonal element of the lambda operator
C
      IF(IBC.EQ.0) THEN
         ALI1(ND)=RAD1(ND)/ST0(ND)
         ALI1(ND-1)=RAD1(ND-1)/ST0(ND-1)
      END IF
C
      DJTOT=0.
      DO 110 ID=1,ND
         DELTAJ=(RAD1(ID)-RAD(IJ,ID))/(UN+SS0(ID)*ALI1(ID))
         RAD(IJ,ID)=RAD(IJ,ID)+DELTAJ
         DJTOT=MAX(DJTOT,ABS(DELTAJ/RAD(IJ,ID)))
  110 CONTINUE
      IF(DJTOT.GT.DJMAX.AND.ITRALI.LE.NTRALI) GO TO 10
C
C     end of ALI loop for electron scattering
C
      DO 120 ID=1,ND
         RAD1(ID)=RAD(IJ,ID)
         FAK1(ID)=RDK(ID)/RAD(IJ,ID)
         FKK(ID)=FAK1(ID)
  120 CONTINUE
      FLUX(IJ)=AHout*half
      FH(IJ)=AH/RAD1(1)-HALF*ALB1
      FH0=FH(IJ)
      fhd(ij)=ahd/rad1(nd)
C
C ********************
C
C     Again solution of the transfer equation, now with Eddington
C     FKK and FH determined above, to insure strict consistency of the
C     radiation field and Eddington factors
C
C     Upper boundary condition
C
      if(ilmcor.eq.3) then
         do id=1,nd
            sa0(id)=st0(id)
            ss0c(id)=ss0(id)
            ss0(id)=0.
         end do
      end if
      ID=1
      DTP1=DT(ID)
      IF(MOD(ISPLIN,3).EQ.0) THEN
         B=DTP1*HALF
         C=0.
       ELSE
         B=DTP1*THIRD
         C=B*HALF
      END IF
      BQ=UN/(B+QQ0)
      CQ=C*BQ
      BBB(ID)=(FKK(ID)/DTP1+FH0+B)*BQ+SS0(ID)
      CCC(ID)=(FKK(ID+1)/DTP1)*BQ-CQ*(UN+SS0(ID+1))
      ZZZ(ID)=UN/BBB(ID)
      VLL=SA0(ID)+CQ*SA0(ID+1)
      IF(IWINBL.LT.0) VLL=VLL+HEXTRD(IJ)*BQ
      AANU(ID)=VLL*ZZZ(ID)
      DDD(ID)=CCC(ID)*ZZZ(ID)
      IF(ISPLIN.GT.2) FFF=BBB(ID)/CCC(ID)-UN
C
C     Normal depth point
C
      DO 280 ID=2,ND-1
         DTM1=DTP1
         DTP1=DT(ID)
         DT0=TWO/(DTP1+DTM1)
         AL=UN/DTM1*DT0
         GA=UN/DTP1*DT0
         IF(MOD(ISPLIN,3).EQ.0) THEN
            A=0.
            C=0.
          ELSE IF(ISPLIN.EQ.1) THEN
            A=DTM1*DT0*SIXTH
            C=DTP1*DT0*SIXTH
          ELSE
            A=(UN-HALF*DTP1*DTP1*AL)*SIXTH
            C=(UN-HALF*DTM1*DTM1*GA)*SIXTH
         END IF
         AAA(ID)=AL*FKK(ID-1)-A*(UN+SS0(ID-1))
         CCC(ID)=GA*FKK(ID+1)-C*(UN+SS0(ID+1))
         BBB(ID)=(AL+GA)*FKK(ID)+(UN-A-C)*(UN+SS0(ID))
         VLL=A*SA0(ID-1)+C*SA0(ID+1)+(UN-A-C)*SA0(ID)
         AANU(ID)=VLL+AAA(ID)*AANU(ID-1)
         IF(ISPLIN.LE.2) THEN
            ZZZ(ID)=UN/(BBB(ID)-AAA(ID)*DDD(ID-1))
            DDD(ID)=CCC(ID)*ZZZ(ID)
            AANU(ID)=AANU(ID)*ZZZ(ID)
         ELSE
            SUM=-AAA(ID)+BBB(ID)-CCC(ID)
            FFF=(SUM+AAA(ID)*FFF*DDD(ID-1))/CCC(ID)
            DDD(ID)=UN/(UN+FFF)
            AANU(ID)=AANU(ID)*DDD(ID)/CCC(ID)
         ENDIF
  280 CONTINUE
C
C     Lower boundary condition
C
      ID=ND
      IF(IBC.EQ.0) THEN
         BBB(ID)=FKK(ID)/DTP1+HALF
         AAA(ID)=FKK(ID-1)/DTP1
         VLL=HALF*PLAND+THIRD*DPLAN
       ELSE IF(IBC.LT.4) THEN
         B=UN/DTP1
         A=TWO*B*B
c        BBB(ID)=UN+SS0(ID)+B*TWO*FHD(IJ)+A*FKK(ID)
         BBB(ID)=UN+SS0(ID)+B+A*FKK(ID)
         AAA(ID)=A*FKK(ID-1)
         VLL=SA0(ID)+B*(PLAND+TWOTHR*DPLAN)
       ELSE 
         B=UN/DTP1
         A=TWO*B*B
         BBB(ID)=B+A*FKK(ID)
         AAA(ID)=A*FKK(ID-1)
         VLL=B*(PLAND+TWOTHR*DPLAN)
      END IF
      EEE(ND)=AAA(ID)/BBB(ID)
      ZZZ(ID)=UN/(BBB(ID)-AAA(ID)*DDD(ID-1))
      RAD1(ID)=(VLL+AAA(ID)*AANU(ID-1))*ZZZ(ID)
      FAK1(ID)=FKK(ND)
      ALRH(ID)=ZZZ(ID)
C
C     Backsolution
C
      DO 290 ID=ND-1,1,-1
         EEE(ID)=AAA(ID)/(BBB(ID)-CCC(ID)*EEE(ID+1))
         RAD1(ID)=AANU(ID)+DDD(ID)*RAD1(ID+1)
         FAK1(ID)=FKK(ID)
         ALRH(ID)=ZZZ(ID)/(UN-DDD(ID)*EEE(ID+1))
         ALRM(ID)=0
         ALRP(ID)=0
  290 CONTINUE
c     FLUX(IJ)=FH(IJ)*RAD1(1)-half*half*extrad(ij)
c
C        evaluate approximate Lambda operator
C
C        a) Rybicki-Hummer Lambda^star operator (diagonal)
C           (for JALI = 1)
C
         DO ID=1,ND
            ALIM1(ID)=0.
            ALIP1(ID)=0.
         END DO
         IF(JALI.EQ.1) THEN
           DO ID=1,ND 
             ALI1(ID)=ALRH(ID)
           END DO
c
           IF(IBC.EQ.0) THEN
             ali1(nd-1)=rad1(nd-1)/sa0(nd-1)
             ali1(nd)=rad1(nd)/sa0(nd)
           END IF
C
C        for IFALI>5:
C        tridiagonal Rybicki-Hummer operator (off-diagonal terms)
C
           IF(IFALI.GE.6) THEN
             ALIP1(1)=ALRH(2)*DDD(1)
             DO ID=2,ND-1 
               ALIM1(ID)=ALRH(ID-1)*EEE(ID)
               ALIP1(ID)=ALRH(ID+1)*DDD(ID)
             END DO
             ALIM1(ND)=ALRH(ND-1)*EEE(ND)
             IF(IBC.EQ.0) THEN
               ALIM1(nd)=0.
               ALIM1(nd-1)=0.
               ALIP1(nd)=0.
               ALIP1(nd-1)=0.
             END IF
           END IF
c
C        b) diagonal Olson-Kunasz Lambda^star operator, 
C           (for JALI = 2)
C
         ELSE IF(JALI.EQ.2) THEN
           DO ID=1,ND-1
             ALI0(ID)=0.
             DO I=1,NMU
               DIV=DT(ID)/AMU(I)
               ALI0(ID)=ALI0(ID)+(UN-EXP(-DIV))/DIV*WTMU(I)
             END DO
           END DO
           DO ID=2,ND-1
             ALI1(ID)=UN-HALF*(ALI0(ID)+ALI0(ID-1))
           END DO
           ALI1(1)=UN-HALF*(ALI0(1)+US0)
           ALI1(ND)=UN-ALI0(ND-1)
           ali1(nd-1)=rad1(nd-1)/sa0(nd-1)
           ali1(nd)=rad1(nd)/sa0(nd)
       END IF
C          
C        correction of Lambda^star for scattering
C  
         IF(ILMCOR.EQ.1) THEN     
           DO ID=1,ND
             ALI1(ID)=ALI1(ID)*(UN+SS0(ID))
             ALIM1(ID)=ALIM1(ID)*(UN+SS0(ID))
             ALIP1(ID)=ALIP1(ID)*(UN+SS0(ID))
           END DO 
           IF(IBC.EQ.4) THEN
             ALI1(ND)= ALI1(ND)/(UN+SS0(ND))
             ALIM1(ND)= ALIM1(ND)/(UN+SS0(ND))
             ALIP1(ND)= ALIP1(ND)/(UN+SS0(ND))
           END IF
         END IF
C
      RETURN
      END


C
C
C     ****************************************************************
C
C
      SUBROUTINE RTEDF2(IJ)
C     =====================
C
C     Solution of the radiative transfer equation - for one frequency 
C     for the known source function.
C     Determination of the radiation field and variable Eddington
C     factors.
C     Analogous to RTEDF1, but using opacity and emissivity
c     instead of the source function
C
C     The numerical method used:
c     Discontinuous Finite Element method
c     Castor, Dykema, Klein, 1992, ApJ 387, 561.
C
C     different formulation of the boundary conditions
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'MODELQ.FOR'
      INCLUDE 'ALIPAR.FOR'
      PARAMETER (SIXTH=UN/6.D0,
     *           THIRD=UN/3.D0,
     *           TWOTHR=TWO/3.D0,
     *           three=3.d0, 
     *           quart=0.25d0)
      DIMENSION DT(MDEPTH),RDK(MDEPTH),FKK(MDEPTH),
     *          ST0(MDEPTH),SA0(MDEPTH),SS0(MDEPTH),
     *          dtau(mdepth),rip(mdepth),rim(mdepth),
     *          riin(mdepth),riup(mdepth),u(mdepth),
     *          aip(mdepth),aim(mdepth),al0(mdepth),
     *          aiin(mdepth),aiup(mdepth),
     *          chip0(mdepth),chim0(mdepth),ddm0(mdepth),
     *          chip(mdepth),chim(mdepth),
     *          etap(mdepth),etam(mdepth)
C
      FR=FREQ(IJ)
C
C     optical depth scale
C
      DO ID=1,ND-1
         DT(ID)=DELDMZ(ID)*(ABSOT(ID+1)+ABSOT(ID))
         SA0(ID)=EMIS1(ID)/ABSO1(ID)
         SS0(ID)=-SCAT1(ID)/ABSO1(ID)
         ddm0(id)=(dm(id+1)-dm(id))
         chip0(id)=(abso1(id)*dens1(id)+three*
     *              abso1(id+1)*dens1(id+1))*quart*ddm0(id)
         chim0(id)=(abso1(id)*dens1(id)*three+
     *              abso1(id+1)*dens1(id+1))*quart*ddm0(id)
      END DO
      SA0(ND)=EMIS1(ND)/ABSO1(ND)
      SS0(ND)=-SCAT1(ND)/ABSO1(ND)
C
      TAUMIN=ABSO1(1)*DEDM1
C
C     Allowance for wind blanketing
C
      ALB1=0.
      IF(IWINBL.GT.0) ALB1=TWO*ALBE(IJ)/(UN+ALBE(IJ))
C
C     quantities for the lower boundary condition
C
      FR15=FR*1.D-15
      BNU=BN*FR15*FR15*FR15
      PLAND=BNU/(EXP(HK*FR/TEMP(ND))-UN)*RRDIL
      DPLAN=BNU/(EXP(HK*FR/TEMP(ND-1))-UN)*RRDIL
      IF(TEMPBD.GT.0.) THEN
        PLAND=BNU/(EXP(HK*FR/TEMPBD)-UN)*RRDIL
        DPLAN=BNU/(EXP(HK*FR/TEMPBD)-UN)*RRDIL
      ENDIF
      DPLAN=(PLAND-DPLAN)/DT(ND-1)
C
C     global ALI loop for treating electron scattering
C
      itrali=0
   10 itrali=itrali+1
C
C     total source function
C  
      DO 20 ID=1,ND
         ST0(ID)=SA0(ID)-SS0(ID)*RAD(IJ,ID)
         RAD1(ID)=0.
         RDK(ID)=0.
         ALI1(ID)=0.
   20 CONTINUE
      AH=0.
      U0=0.
      QQ0=0.
      US0=0.
c
c     loop over angle poits
c
      DO 100 I=1,NMU
         AMU2=AMU(I)*AMU(I)*WTMU(I)
         do id=1,nd-1
            dtau(id)=dt(id)/amu(i)
            chip(id)=un+chip0(id)/amu(i)
            chim(id)=un+chim0(id)/amu(i)
            etap(id)=emis1(id+1)*dens1(id+1)/amu(i)*ddm0(id)
            etam(id)=emis1(id)*dens1(id)/amu(i)*ddm0(id)
         enddo
c               
c        incoming intensity 
c   
         ID=1
         P0=0.
         EX=UN
C
C        allowance for non-zero optical depth at the first depth point
C
c        rim(id)=EXTRAD(IJ)
         rim(id)=EXTINT(IJ,I)
         IF(IWINBL.EQ.0) THEN
            TAMM=TAUMIN/AMU(I)
            EX=EXP(-TAMM)
            P0=UN-EX
            QQ0=QQ0+P0*AMU(I)*WTMU(I)
            U0=U0+EX*WTMU(I)
            US0=US0+P0/TAMM*WTMU(I)
            rim(id)=st0(id)*p0
         END IF
c               
c          incoming intensity 
c   
            aim(id)=0.
            do id=1,nd-1
               dt0=dtau(id)
               dtaup1=dt0+un
               dtau2=dt0*dt0
               bb=two*dtaup1
               cc=dt0*dtaup1
               aa=un/(dtau2+bb)
               aam=un/(un+chim(id)*chip(id))
               rim(id+1)=(two*rim(id)+etap(id)*chim(id)+etam(id))*aam
               rip(id)=(two*rim(id)*chim(id)+etam(id)*chip(id)-
     *                  etap(id))*aam
               aim(id+1)=bb*aa
               aip(id)=(cc+bb*aim(id))*aa
            enddo
            do id=2,nd-1
c               riin(id)=(rim(id)*dtau(id)+rip(id)*dtau(id-1))*dtt
               dtt=un/(dtau(id-1)+dtau(id))
               dtm=un/(ddm0(id-1)+ddm0(id))
               riin(id)=(rim(id)*ddm0(id)+rip(id)*ddm0(id-1))*dtm
               aiin(id)=(aim(id)*dtau(id)+aip(id)*dtau(id-1))*dtt
            enddo
            riin(1)=rim(1)
            riin(nd)=rim(nd)
            aiin(1)=aim(1)
            aiin(nd)=aim(nd)
C              
c           outgoing intensity
c               
            rim(nd)=PLAND+AMU(I)*DPLAN
            do id=nd-1,1,-1
               dt0=dtau(id)
               dtaup1=dt0+un
               dtau2=dt0*dt0
               bb=two*dtaup1
               cc=dt0*dtaup1
               aa=un/(dtau2+bb)
               aam=un/(un+chim(id)*chip(id))
               rim(id)=(two*rim(id+1)+etam(id)*chip(id)+etap(id))*aam
               rip(id+1)=(two*rim(id+1)*chip(id)+etap(id)*chim(id)-
     *                    etam(id))*aam
               aim(id)=cc*aa
               aip(id+1)=(cc+bb*aim(id+1))*aa
            enddo
            do id=2,nd-1
c               riup(id)=(rim(id)*dtau(id-1)+rip(id)*dtau(id))*dtt
               dtt=un/(dtau(id-1)+dtau(id))
               dtm=un/(ddm0(id-1)+ddm0(id))
               riup(id)=(rim(id)*ddm0(id-1)+rip(id)*ddm0(id))*dtm
               aiup(id)=(aim(id)*dtau(id-1)+aip(id)*dtau(id))*dtt
            enddo
            riup(1)=rim(1)
            riup(nd)=rim(nd)
            aiup(1)=aim(1)
            aiup(nd)=aim(nd)
c                
c           final symmetrized (Feautrier) intensity -- (riin+riup)/2
c                
            do id=1,nd
               u(id)=(riin(id)+riup(id))*half
               al0(id)=(aiin(id)+aiup(id))*half
            enddo
c             
         DO 70 ID=1,ND
            RAD1(ID)=RAD1(ID)+WTMU(I)*U(ID)
            RDK(ID)=RDK(ID)+AMU2*U(ID)
            ALI1(ID)=ALI1(ID)+WTMU(I)*AL0(ID)
   70    CONTINUE
         AH=AH+AMU(I)*WTMU(I)*U(1)
c
c     end of the loop over angle points
c
  100 CONTINUE
C
C     solution of the transfer equation
C     Variables:
C     RAD1    - mean intensity
C     FAK1    - Eddington factor f(K) = K/J
C     FH      - the "surface" Eddington factor
C     ALI1    - diagonal element of the lambda operator
C
      IF(IBC.EQ.0) THEN
         ALI1(ND)=RAD1(ND)/ST0(ND)
         ALI1(ND-1)=RAD1(ND-1)/ST0(ND-1)
      END IF
C
      DJTOT=0.
      DO 110 ID=1,ND
         DELTAJ=(RAD1(ID)-RAD(IJ,ID))/(UN+SS0(ID)*ALI1(ID))
         RAD(IJ,ID)=RAD(IJ,ID)+DELTAJ
         DJTOT=MAX(DJTOT,ABS(DELTAJ/RAD(IJ,ID)))
  110 CONTINUE
      IF(DJTOT.GT.DJMAX.AND.ITRALI.LE.NTRALI) GO TO 10
C
C     end of ALI loop for electron scattering
C
      DO 120 ID=1,ND
         RAD1(ID)=RAD(IJ,ID)
         FAK1(ID)=RDK(ID)/RAD(IJ,ID)
         FKK(ID)=FAK1(ID)
  120 CONTINUE
      FLUX(IJ)=AH
      FH(IJ)=AH/RAD1(1)-HALF*ALB1
      FH0=FH(IJ)
C          
C        correction of Lambda^star for scattering
C  
         IF(ILMCOR.EQ.1) THEN     
         DO ID=1,ND
             ALI1(ID)=ALI1(ID)*(UN+SS0(ID))
         END DO   
         IF(IBC.EQ.4) THEN
             ALI1(ND)= ALI1(ND)/(UN+SS0(ND))
         END IF
         ELSE IF(ILMCOR.EQ.3) THEN     
         DO ID=1,ND
             ALI1(ID)=ALI1(ID)/(UN+SS0(ID)*ALI1(ID))
         END DO   
         IF(IBC.EQ.4) THEN
             ALI1(ND)= ALI1(ND)*(UN+SS0(ND)*ALI1(ID))
         END IF
         END IF
C
      RETURN
      END
C
C
C     ****************************************************************
C
C

      SUBROUTINE RTEFR1(IJ)
C     =====================
C
C     Solution of the radiative transfer equation - for one frequency
C     - for the known source function.
C     Determination of the radiation field and variable Eddington
C     factors.
C
C     The numerical method used:
C      for ISPLIN = 0  -  the ordinary Feautrier scheme
C                 = 1  -  the spline collocation method
C                 = 2  -  Hermitian fourth-order method
C                 = 3  -  improved Feautrier scheme
C                         (Rybicki & Hummer 1991, A&A 245, 171.)
C
C     In all cases, the overall matrix system is solved by the standard
C     Gaussian elimination, analogous to that described in SOLVE
C     (auxiliary matrix D is called ALF in SOLVE; auxiliary vector ANU
C     is called BET in SOLVE)
C
C      U0     - derivative of Q0 wrt taumin;
C               for "fixed" frequencies, U0 has the meaning of
C               absorption coefficient * second moment H
C               ( a quantity needed for lower boundary condition of the
C               hydrostatic equilibrium equation, specifically for
C               accounting for an effect of fixed-option transitions)
C
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'MODELQ.FOR'
      INCLUDE 'ALIPAR.FOR'
      INCLUDE 'ITERAT.FOR'
      PARAMETER (SIXTH=UN/6.D0,
     *           THIRD=UN/3.D0,
     *           TWOTHR=TWO/3.D0)
      DIMENSION AANU(MDEPTH),DDD(MDEPTH),FKK(MDEPTH),
     *          RDD(MDEPTH),ST0(MDEPTH),SS0(MDEPTH),AB0(MDEPTH),
     *          AA(MMU,MMU),BB(MMU,MMU),CC(MMU,MMU),VL(MMU),
     *          FFD(MMU,MMU),FF0D(MMU,MMU),
     *          FFPD(MMU,MMU),ali0(mdepth),ss0c(mdepth),
     *          AAA(MDEPTH),BBB(MDEPTH),CCC(MDEPTH),EEE(MDEPTH),
     *          ZZZ(MDEPTH),ALRH(MDEPTH),ALRM(MDEPTH),ALRP(MDEPTH),
     *          D(MMU,MMU,MDEPTH),ANU(MMU,MDEPTH),scor(mdepth)
      DIMENSION rmmu(2*MMU),wmmu(2*MMU),rwmu(2*MMU),
     *          dtau(mdepth),ri(mdepth),ali(mdepth),alij1(mdepth)
      COMMON/OPTDPT/DT(MDEPTH)
C
      WW=W(IJ)
      ISPL=ISPLIN
      IF(ISPLIN.GE.5) THEN
         ISPLIN=ISPL-5
         IF(IJALI(IJ).GT.0) THEN
            IF(IRTE.EQ.0) THEN
               CALL RTEDF1(IJ)
             ELSE 
               CALL RTEDF2(IJ)
            END IF
            ISPLIN=ISPL
            if(ifprad.eq.0) return
            DO ID=1,ND
               if(.not.lskip(ID,IJ))
     *         PRADT(ID)=PRADT(ID)+RAD1(ID)*FAK1(ID)*W(ij)
            END DO
            if(.not.lskip(1,IJ))
     *         PRD0=PRD0+ABSO1(1)*W(IJ)*(RAD1(1)*FH(IJ)-HEXTRD(IJ))
            DO ID=1,ND
               PRADA(ID)=PRADA(ID)+RAD1(ID)*FAK1(ID)*WW
            END DO
            RETURN
         END IF
      END IF
C
      if(icompt.gt.0.and.(iter.gt.1.or.ilam.gt.0)) then
         call rtecf1(ij)
       return
      end if
c
C
      FR=FREQ(IJ)
C
C     total source function
C  
      AH=0.
      DO 20 ID=1,ND
         AB0(ID)=ABSO1(ID)
         ST0(ID)=EMIS1(ID)/AB0(ID)
         RAD1(ID)=0.
         ALI1(ID)=0.
   20 CONTINUE
C
C     non-coherent electron scattering by lambda iteration
C  
      IF(NELSC.LE.0) THEN
         DO 22 ID=1,ND
            SS0(ID)=-SCAT1(ID)/AB0(ID)
   22    CONTINUE
       ELSE
         DO 24 ID=1,ND
            ST0(ID)=ST0(ID)+SCAT1(ID)*EMEL1(ID)*RAD1(ID)/AB0(ID)
            SS0(ID)=0.
   24    CONTINUE
      END IF
C
C     optical depth scale
C
      DO ID=1,ND-1
         DT(ID)=DELDMZ(ID)*(ABSOT(ID+1)+ABSOT(ID))
      END DO
c
      U0=0.
      QQ0=0.
      US0=0.
      TAUMIN=ABSO1(1)*DEDM1
C
      ALB1=0.
C
C     Allowance for wind blanketing
C
      IF(IWINBL.GT.0) ALB1=TWO*ALBE(IJ)/(UN+ALBE(IJ))
C
C ************** Forward elimination
C
C     Upper boundary condition
C
      ID=1
      DTP1=DT(1)
      P0=0.
      EX=UN
      IF(MOD(ISPLIN,3).EQ.0) THEN
         B=DTP1*HALF
         C=0.
       ELSE
         B=DTP1*THIRD
         C=B*HALF
      END IF
      QQ0=0.
      US0=0.
      DO 60 I=1,NMU
c         IF(IWINBL.EQ.0) THEN
C
C        allowance for non-zero optical depth at the first depth point
C
c            TAMM=TAUMIN/AMU(I)
c            EX=EXP(-TAMM)
c            P0=UN-EX
c            QQ0=QQ0+P0*AMU(I)*WTMU(I)
c            U0=U0+EX*WTMU(I)
c            US0=US0+P0/TAMM*WTMU(I)
c         END IF
C
         BI=B/AMU(I)
         CI=C/AMU(I)
         VL(I)=(BI+P0)*ST0(ID)+CI*ST0(ID+1)
c        IF(IWINBL.LT.0) VL(I)=VL(I)+EXTRAD(IJ)
         IF(IWINBL.LT.0) VL(I)=VL(I)+EXTINT(IJ,I)
         DO 40 J=1,NMU
            BB(I,J)=SS0(ID)*WTMU(J)*(BI+P0)-ALB1*WTMU(J)
   40       CC(I,J)=-CI*SS0(ID+1)*WTMU(J)
         BB(I,I)=BB(I,I)+AMU(I)/DTP1+UN+BI
         CC(I,I)=CC(I,I)+AMU(I)/DTP1-CI
         ANU(I,ID)=0.
   60 CONTINUE
C
      IF(ISPLIN.LE.2) THEN
         CALL MATINV(BB,NMU,MMU)
         DO 70 I=1,NMU
            DO 70 J=1,NMU
               D(I,J,ID)=0.
               DO 71 K=1,NMU
   71             D(I,J,ID)=D(I,J,ID)+BB(I,K)*CC(K,J)
               ANU(I,1)=ANU(I,1)+BB(I,J)*VL(J)
   70    CONTINUE
      ELSE
         DO 80 I=1,NMU
            DO 81 J=1,NMU
               FF0D(I,J)=BB(I,J)/CC(I,I)
   81       CONTINUE
            FF0D(I,I)=FF0D(I,I)-UN
   80    CONTINUE
C
c        CALL MINV3(BB)
         CALL MATINV(BB,NMU,MMU)
         DO 82 I=1,NMU
            ANU(I,ID)=0.
            DO 82 J=1,NMU
               D(I,J,ID)=BB(I,J)*CC(J,J)
               ANU(I,ID)=ANU(I,ID)+BB(I,J)*VL(J)
   82    CONTINUE
      ENDIF
C
C     Normal depth points  1 < ID < ND
C
      DO 200 ID=2,ND-1
         DTM1=DTP1
         DTP1=DT(ID)
         DT0=TWO/(DTM1+DTP1)
         AL=UN/DTM1*DT0
         GA=UN/DTP1*DT0
         BE=AL+GA
         IF(MOD(ISPLIN,3).EQ.0) THEN
            A=0.
            C=0.
          ELSE IF(ISPLIN.EQ.1) THEN
            A=DTM1*DT0*SIXTH
            C=DTP1*DT0*SIXTH
          ELSE
            A=(UN-HALF*AL*DTP1*DTP1)*SIXTH
            C=(UN-HALF*GA*DTM1*DTM1)*SIXTH
         END IF
         B=UN-A-C
         VL0=A*ST0(ID-1)+B*ST0(ID)+C*ST0(ID+1)
         DO 110 I=1,NMU
            DO 110 J=1,NMU
               AA(I,J)=-A*SS0(ID-1)*WTMU(J)
               CC(I,J)=-C*SS0(ID+1)*WTMU(J)
               BB(I,J)=B*SS0(ID)*WTMU(J)
  110    CONTINUE
         DO 120 I=1,NMU
            VL(I)=VL0
            DIV=AMU(I)*AMU(I)
            AA(I,I)=AA(I,I)+DIV*AL-A
            CC(I,I)=CC(I,I)+DIV*GA-C
            BB(I,I)=BB(I,I)+DIV*BE+B
  120    CONTINUE
         DO 130 I=1,NMU
            DO 130 J=1,NMU
               VL(I)=VL(I)+AA(I,J)*ANU(J,ID-1)
  130    CONTINUE
         IF(ISPLIN.LE.2) THEN
            DO 140 I=1,NMU
               DO 140 J=1,NMU
                  S=0.
                  DO 141 K=1,NMU
  141                S=S+AA(I,K)*D(K,J,ID-1)
                  BB(I,J)=BB(I,J)-S
  140       CONTINUE
C
            CALL MATINV(BB,NMU,MMU)
            DO 150 I=1,NMU
               DO 150 J=1,NMU
                  D(I,J,ID)=0.
                  DO 151 K=1,NMU
  151                D(I,J,ID)=D(I,J,ID)+BB(I,K)*CC(K,J)
  150       CONTINUE
         ELSE
            DO 160 I=1,NMU
               BB(I,I)=-AA(I,I)+BB(I,I)-CC(I,I)
               DO 160 J=1,NMU
                  FFPD(I,J)=AA(I,I)*FF0D(I,J)
  160       CONTINUE
            DO 162 I=1,NMU
               DO 162 J=1,NMU
                  S=0.
                  DO 163 K=1,NMU
  163                S=S+FFPD(I,K)*D(K,J,ID-1)
                  FFD(I,J)=(BB(I,J)+S)/CC(I,I)
  162       CONTINUE
            DO 164 I=1,NMU
               DO 165 J=1,NMU
  165             FF0D(I,J)=FFD(I,J)
               FFD(I,I)=FFD(I,I)+UN
  164       CONTINUE
C
            CALL MATINV(FFD,NMU,MMU)
            DO 170 I=1,NMU
               DO 170 J=1,NMU
                  D(I,J,ID)=FFD(I,J)
                  BB(I,J)=FFD(I,J)/CC(J,J)
  170       CONTINUE
         ENDIF
         DO 180 I=1,NMU
            ANU(I,ID)=0.
            DO 180 J=1,NMU
               ANU(I,ID)=ANU(I,ID)+BB(I,J)*VL(J)
  180    CONTINUE
  200 CONTINUE
C
C     Lower boundary condition
C
      ID=ND
C
C     First option: 
C     b.c. is different from stellar atmospheres; expresses symmetry
C     at the central plane   I(taumax,-mu,nu)=I(taumax,+mu,nu)
C
      IF(IFZ0.GE.0.AND.IDISK.EQ.1) THEN
         B=DTP1*HALF
         A=0.
         DO 204 I=1,NMU
            BI=B/AMU(I)
            AI=A/AMU(I)
            VL(I)=ST0(ID)*BI+ST0(ID-1)*AI
            DO 202 J=1,NMU
               AA(I,J)=-AI*SS0(ID-1)*WTMU(J)
  202          BB(I,J)=BI*SS0(ID)*WTMU(J)
            AA(I,I)=AA(I,I)+AMU(I)/DTP1-AI
            BB(I,I)=BB(I,I)+AMU(I)/DTP1+BI
  204    CONTINUE
         DO 208 I=1,NMU
            S1=0.
            DO 207 J=1,NMU
               S=0.
               S1=S1+AA(I,J)*ANU(J,ID-1)
               DO 206 K=1,NMU
  206             S=S+AA(I,K)*D(K,J,ID-1)
  207          BB(I,J)=BB(I,J)-S
            VL(I)=VL(I)+S1
  208    CONTINUE
C
C     Second option:
C     b.c. is the same as in stellar atmospheres - the last depth point
C     is not at the central plane
C
      ELSE
      FR15=FR*1.D-15
      BNU=BN*FR15*FR15*FR15
      PLAND=BNU/(EXP(HK*FR/TEMP(ND ))-UN)*RRDIL
      DPLAN=BNU/(EXP(HK*FR/TEMP(ND-1))-UN)*RRDIL
      IF(TEMPBD.GT.0.) THEN
        PLAND=BNU/(EXP(HK*FR/TEMPBD)-UN)*RRDIL
        DPLAN=BNU/(EXP(HK*FR/TEMPBD)-UN)*RRDIL
      ENDIF
      DPLAN=(PLAND-DPLAN)/DT(ND-1)
      IF(IBC.EQ.0.OR.IBC.EQ.4) THEN
         DO 220 I=1,NMU
            AA(I,I)=AMU(I)/DTP1
            VL(I)=PLAND+AMU(I)*DPLAN+AA(I,I)*ANU(I,ID-1)
            DO 210 J=1,NMU
  210          BB(I,J)=-AA(I,I)*D(I,J,ID-1)
            BB(I,I)=BB(I,I)+AA(I,I)+UN
  220    CONTINUE
       ELSE
         DO 224 I=1,NMU
            A=AMU(I)/DTP1
            B=HALF/A
            AA(I,I)=A
            VL(I)=B*ST0(ID)+PLAND+AMU(I)*DPLAN+AA(I,I)*ANU(I,ID-1)
            DO 222 J=1,NMU
  222          BB(I,J)=B*SS0(ID)*WTMU(J)-AA(I,I)*D(I,J,ID-1)
            BB(I,I)=BB(I,I)+A+B+UN
  224    CONTINUE
      END IF
      END IF
C
         CALL MATINV(BB,NMU,MMU)
C
      DO 230 I=1,NMU
         ANU(I,ID)=0.
         DO 230 J=1,NMU
            D(I,J,ID)=0.
            ANU(I,ID)=ANU(I,ID)+BB(I,J)*VL(J)
  230 CONTINUE
C
C ***************** Backsolution
C
      FKK(ND)=THIRD
      AJ=0.
      AH=0.
      AK=0.
      DO 235 I=1,NMU
         RMU=WTMU(I)*ANU(I,ID)
         AJ=AJ+RMU
         AH=AH+RMU*AMU(I)
         AK=AK+RMU*AMU(I)*AMU(I)
  235 CONTINUE
      RDD(ID)=AJ
      IF(IBC.EQ.0) THEN
         FKK(ND)=THIRD
       ELSE
         FKK(ID)=AK/AJ
         FHD(IJ)=AH/AJ
      END IF
C
      DO 260 ID=ND-1,1,-1
         DO 240 I=1,NMU
            DO 240 J=1,NMU
               ANU(I,ID)=ANU(I,ID)+D(I,J,ID)*ANU(J,ID+1)
  240    CONTINUE
         AJ=0.
         AK=0.
         DO 250 I=1,NMU
            RMU=WTMU(I)*ANU(I,ID)
            AJ=AJ+RMU
            AK=AK+RMU*AMU(I)*AMU(I)
  250    CONTINUE
C
C         solution of the transfer equation
C         Variables:
C          ANU     - Feautrier intensity
C          RDD     - mean intensity
C          FKK     - Eddington factor f(K) = K/J
C
         FKK(ID)=AK/AJ
         RDD(ID)=AJ
  260 CONTINUE
C
      if(idisk.eq.1) then
        do id=1,nd
          fak(ij,id)=fkk(id)
        end do
      end if
C
C     the "surface" Eddington factor fH
C
      AH=0.
      DO 270 I=1,NMU
  270    AH=AH+WTMU(I)*AMU(I)*ANU(I,1)
      FH0=AH/AJ-HALF*ALB1
      FH(IJ)=FH0
C
C ********************
C
C     Again solution of the transfer equation, now with Eddington
C     FKK and FH determined above, to insure strict consistency of the
C     radiation field and Eddington factors
C
      if(ilmcor.eq.2) then
         do id=1,nd
            scor(id)=un/(un+ss0(id))
         end do
      else if(ilmcor.eq.3) then
         do id=1,nd
            ss0c(id)=ss0(id)
            st0(id)=st0(id)-ss0(id)*rdd(id)
            ss0(id)=0.
         end do
      end if
C
C     Upper boundary condition
C
      ID=1
      DTP1=DT(ID)
      IF(MOD(ISPLIN,3).EQ.0) THEN
         B=DTP1*HALF
         C=0.
       ELSE
         B=DTP1*THIRD
         C=B*HALF
      END IF
      BQ=UN/(B+QQ0)
      CQ=C*BQ
      BBB(ID)=(FKK(ID)/DTP1+FH0+B)*BQ+SS0(ID)
      CCC(ID)=(FKK(ID+1)/DTP1)*BQ-CQ*(UN+SS0(ID+1))
      VLL=ST0(ID)+CQ*ST0(ID+1)
      IF(IWINBL.LT.0) VLL=VLL+HEXTRD(IJ)*BQ
      if(ilmcor.eq.2) then
         bbb(id)=bbb(id)*scor(id)
         ccc(id)=ccc(id)*scor(id)
         vll=vll*scor(id)
      end if
      ZZZ(ID)=UN/BBB(ID)
      AANU(ID)=VLL*ZZZ(ID)
      DDD(ID)=CCC(ID)*ZZZ(ID)
      IF(ISPLIN.GT.2) FFF=BBB(ID)/CCC(ID)-UN
C
C     Normal depth point
C
      DO 280 ID=2,ND-1
         DTM1=DTP1
         DTP1=DT(ID)
         DT0=TWO/(DTP1+DTM1)
         AL=UN/DTM1*DT0
         GA=UN/DTP1*DT0
         IF(MOD(ISPLIN,3).EQ.0) THEN
            A=0.
            C=0.
          ELSE IF(ISPLIN.EQ.1) THEN
            A=DTM1*DT0*SIXTH
            C=DTP1*DT0*SIXTH
          ELSE
            A=(UN-HALF*DTP1*DTP1*AL)*SIXTH
            C=(UN-HALF*DTM1*DTM1*GA)*SIXTH
         END IF
         AAA(ID)=AL*FKK(ID-1)-A*(UN+SS0(ID-1))
         CCC(ID)=GA*FKK(ID+1)-C*(UN+SS0(ID+1))
         BBB(ID)=(AL+GA)*FKK(ID)+(UN-A-C)*(UN+SS0(ID))
         VLL=A*ST0(ID-1)+C*ST0(ID+1)+(UN-A-C)*ST0(ID)
         if(ilmcor.eq.2) then
            aaa(id)=aaa(id)*scor(id)
            bbb(id)=bbb(id)*scor(id)
            ccc(id)=ccc(id)*scor(id)
            vll=vll*scor(id)
         end if
         AANU(ID)=VLL+AAA(ID)*AANU(ID-1)
         IF(ISPLIN.LE.2) THEN
            ZZZ(ID)=UN/(BBB(ID)-AAA(ID)*DDD(ID-1))
            DDD(ID)=CCC(ID)*ZZZ(ID)
            AANU(ID)=AANU(ID)*ZZZ(ID)
         ELSE
            SUM=-AAA(ID)+BBB(ID)-CCC(ID)
            FFF=(SUM+AAA(ID)*FFF*DDD(ID-1))/CCC(ID)
            DDD(ID)=UN/(UN+FFF)
            AANU(ID)=AANU(ID)*DDD(ID)/CCC(ID)
         ENDIF
  280 CONTINUE
C
C     Lower boundary condition
C
      ID=ND
c
c     stellar atmospheric
c
      IF(IDISK.EQ.0.OR.IFZ0.LT.0) then
      IF(IBC.EQ.0) THEN
         BBB(ID)=FKK(ID)/DTP1+HALF
         AAA(ID)=FKK(ID-1)/DTP1
         VLL=HALF*PLAND+THIRD*DPLAN
       ELSE IF(IBC.LT.4) THEN
         B=UN/DTP1
         A=TWO*B*B
         BBB(ID)=UN+SS0(ID)+B*TWO*FHD(IJ)+A*FKK(ID)
         AAA(ID)=A*FKK(ID-1)
         VLL=ST0(ID)+B*(PLAND+TWOTHR*DPLAN)
       ELSE 
         B=UN/DTP1
         A=TWO*B*B
         BBB(ID)=B+A*FKK(ID)
         AAA(ID)=A*FKK(ID-1)
         VLL=B*(PLAND+TWOTHR*DPLAN)
      END IF
c
c     accretion disk - symmetric boundary
c
      ELSE
         B=TWO/DTP1
         BBB(ID)=FKK(ID)/DTP1*B+UN+SS0(ND)
         AAA(ID)=FKK(ID-1)/DTP1*B
         VLL=ST0(ID)
      END IF
      if(ilmcor.eq.2) then
         aaa(id)=aaa(id)*scor(id)
         bbb(id)=bbb(id)*scor(id)
         vll=vll*scor(id)
      end if
      EEE(ND)=AAA(ID)/BBB(ID)
      ZZZ(ID)=UN/(BBB(ID)-AAA(ID)*DDD(ID-1))
      RAD1(ID)=(VLL+AAA(ID)*AANU(ID-1))*ZZZ(ID)
      FAK1(ID)=FKK(ND)
      ALRH(ID)=ZZZ(ID)
C
C     Backsolution
C
      DO 290 ID=ND-1,1,-1
         EEE(ID)=AAA(ID)/(BBB(ID)-CCC(ID)*EEE(ID+1))
         RAD1(ID)=AANU(ID)+DDD(ID)*RAD1(ID+1)
         FAK1(ID)=FKK(ID)
         ALRH(ID)=ZZZ(ID)/(UN-DDD(ID)*EEE(ID+1))
         ALRM(ID)=0
         ALRP(ID)=0
  290 CONTINUE
c     FLUX(IJ)=FH(IJ)*RAD1(1)-half*half*extrad(ij)
      FLUX(IJ)=FH(IJ)*RAD1(1)-half*hextrd(ij)
c
C        evaluate approximate Lambda operator
C
C        a) Rybicki-Hummer Lambda^star operator (diagonal)
C           (for JALI = 1)
C
         DO ID=1,ND
            ALIM1(ID)=0.
            ALIP1(ID)=0.
         END DO
         IF(JALI.EQ.1) THEN
           DO ID=1,ND 
             ALI1(ID)=ALRH(ID)
           END DO
c
           IF(IBC.EQ.0) THEN
             ali1(nd-1)=rad1(nd-1)/st0(nd-1)
             ali1(nd)=rad1(nd)/st0(nd)
           END IF
C
C        for IFALI>5:
C        tridiagonal Rybicki-Hummer operator (off-diagonal terms)
C
           IF(IFALI.GE.6) THEN
             ALIP1(1)=ALRH(2)*DDD(1)
             DO ID=2,ND-1 
              ALIM1(ID)=ALRH(ID-1)*EEE(ID)
              ALIP1(ID)=ALRH(ID+1)*DDD(ID)
             END DO
             ALIM1(ND)=ALRH(ND-1)*EEE(ND)
             IF(IBC.EQ.0) THEN
              ALIM1(nd)=0.
              ALIM1(nd-1)=0.
              ALIP1(nd)=0.
              ALIP1(nd-1)=0.
             END IF
           END IF
c
C        b) diagonal Olson-Kunasz Lambda^star operator, 
C           (for JALI = 2)
C
         ELSE IF(JALI.EQ.2) THEN
           DO ID=1,ND-1
             ALI0(ID)=0.
             DO I=1,NMU
               DIV=DT(ID)/AMU(I)
               ALI0(ID)=ALI0(ID)+(UN-EXP(-DIV))/DIV*WTMU(I)
             END DO
           END DO
           DO ID=2,ND-1
             ALI1(ID)=UN-HALF*(ALI0(ID)+ALI0(ID-1))
           END DO
           ALI1(1)=UN-HALF*(ALI0(1)+US0)
           ALI1(ND)=UN-ALI0(ND-1)
           ali1(nd-1)=rad1(nd-1)/st0(nd-1)
           ali1(nd)=rad1(nd)/st0(nd)
       END IF
C          
C        correction of Lambda^star for scattering
C  
         IF(ILMCOR.EQ.1) THEN     
           DO ID=1,ND
             ALI1(ID)=ALI1(ID)*(UN+SS0(ID))
             ALIM1(ID)=ALIM1(ID)*(UN+SS0(ID))
             ALIP1(ID)=ALIP1(ID)*(UN+SS0(ID))
           END DO   
           IF(IBC.EQ.4) THEN
             ALI1(ND)= ALI1(ND)/(UN+SS0(ND))
             ALIM1(ND)= ALIM1(ND)/(UN+SS0(ND))
             ALIP1(ND)= ALIP1(ND)/(UN+SS0(ND))
           END IF
         END IF
C
      if(ifalih.gt.0) then
c
c     solution for the individual angles - to get Lambda^star_H
C
      do id=1,nd
         alih1(id)=0.
         alij1(id)=0.
      end do
      nw=nmu
      do i=1,nw
        rmmu(i)=-amu(nw-i+1)
        rmmu(i+nw)=amu(i)
      wmmu(i)=wtmu(nw-i+1)
      wmmu(i+nw)=wtmu(i)
      end do
      do i=1,2*nw
         rwmu(i)=rmmu(i)*wmmu(i)*half
      end do
C
c     --------------------- loop over angles
c
      do i=1,2*nw
         do id=1,nd-1
            dtau(id)=dt(id)/abs(rmmu(i))
         end do
c
c        boundary conditions
c
c        rup=extrad(ij)
         rup=extint(ij,i)
C
C        diffusion approximation for semi-infinite atmospheres
C
         rdown=pland+rmmu(i)*dplan
c
c        solution of the transfer equation
c
         call rtesol(dtau,st0,rup,rdown,rmmu(i),ri,ali)
c
       DO ID=1,ND
            alih1(id)=alih1(id)+rwmu(i)*ali(id)
            alij1(id)=alij1(id)+wmmu(i)*ali(id)*half
         END DO
      end do
c     if(ij.eq.10.or.ij.eq.20.or.ij.eq.50.or.ij.eq.100) then
c        do id=1,nd
c        write(109,609) ij,id,fr,ali1(id),alij1(id),alih1(id)
c        end do
c     end if
c 609 format(2i5,1p4e15.6)
      end if
C
c     --------------------- end of loop over angles
c
c 390 CONTINUE
      ISPLIN=ISPL
c
      if(idisk.ne.0) then
        iji=nfreq-kij(ij)+1
        DO ID=1,ND
        rad(iji,id)=rad1(id)
        END DO
      endif
C
C     radiation pressure
C
      if(ifprad.gt.0) then
      if(.not.lskip(1,IJ))
     *   PRD0=PRD0+ABSO1(1)*WW*(RAD1(1)*FH(IJ)-HEXTRD(IJ))
      DO ID=1,ND
         if(.not.lskip(ID,IJ))
     *      PRADT(ID)=PRADT(ID)+RAD1(ID)*FAK1(ID)*WW
         PRADA(ID)=PRADA(ID)+RAD1(ID)*FAK1(ID)*WW
      END DO
      end if
c
      if(chmax.ge.1.91e-3.and.chmax.le.2.03e-3) then
      tauij=taumin
      do id=1,nd
      if(ilmcor.eq.3) ss0(id)=ss0c(id)
      if(id.gt.1) tauij=tauij+dt(id-1)
      write(97,697) ij,id,tauij,rad1(id),st0(id)/(un+ss0(id)),st0(id),
     *              un+ss0(id),ali1(id)
      end do
  697 format(2i4,1p6e12.4)
      end if
C
      RETURN
      END
C
C
C     ****************************************************************
C
C
      SUBROUTINE OPACF0(ID,NFRQ)
C     ==========================
C
C     Absorption, emission, and scattering coefficients
C     at depth ID
C
C     Input: ID   opacity and emissivity is calculated for the
C                 depth point ID
C     Output: ABSO -  array of absorption coefficient
C             EMIS -  array of emission coefficient
C             SCAT -  array of scattering coefficient 
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'ATOMIC.FOR'
      INCLUDE 'MODELQ.FOR'
      INCLUDE 'ODFPAR.FOR'
      INCLUDE 'ALIPAR.FOR'
      PARAMETER (FRH=3.28805E15, PH2=2.815D29*2., EHB=157802.77355)
      PARAMETER (CFF1=1.3727D-25,CFF2=4.3748D-10,CFF3=2.5993D-7)
      PARAMETER (C14=2.99793D14)
      PARAMETER (SGFF0 = 3.694D8)
      common/hmolab/anh2(mdepth),anhm(mdepth)
      DIMENSION FREDG(NLMX),S(NLMX),SUM(NLMX),PRF(MFREQL)
C
C     initialize
c
C        this part is analogous to TDPINI - for one depth only
C
         T=TEMP(ID)
         T1=UN/T
         HKT1(ID)=HK*T1
         HKT21(ID)=HKT1(ID)*T1
         TK1(ID)=HKT1(ID)/H
         SQT1(ID)=SQRT(T)
         TEMP1(ID)=T1
         CALL GFREE0(ID)
         EMEL1(ID)=UN
         LASER=ITER.GT.ITLAS
c
C        this part is analogous to OPAINI - for one depth only
C
         ANE=ELEC(ID)
         ELEC1(ID)=UN/ANE
         DENS1(ID)=UN/DENS(ID)
         DENSI(ID)=DENS1(ID)
         DENSIM(ID)=DENSI(ID)*WMM(ID)
         ELSCAT(ID)=ANE*SIGE
         if(izscal.eq.1) then
            densi(id)=un
            densim(id)=0.
         end if
         CALL DWNFR0(ID)
         CALL WNSTOR(ID)
         CALL SABOLF(ID)
c
c        quantities for the bound-free opacity
c
         DO 10 IBFT=1,NTRANC
            ITR=ITRBF(IBFT)
          if(indexp(itr).eq.0) go to 10
            II=ILOW(ITR)
            JJ=IUP(ITR)
            IT=ITRA(JJ,II)
            IE=IEL(II)
            NKE=NNEXT(IE)
            CORR=UN
            IF(NKE.NE.JJ) CORR=G(NKE)/G(JJ)*
     *         EXP((ENION(NKE)-ENION(JJ))*TK1(ID))
            ABTRA(ITR,ID)=POPUL(II,ID)
            EMTRA(ITR,ID)=POPUL(JJ,ID)*ANE*SBF(II)*WOP(II,ID)*CORR
   10    CONTINUE
c
c        quantities for the free-free opacity
c
         IF(IELHM.GT.0) THEN
            CFFN(ID)=POPUL(NFIRST(IELH),ID)*ANE
            CFFT(ID)=CFF2-CFF3/T
         END IF
         SGFF=SGFF0/SQT1(ID)*ANE
         DO 20 ION=1,NION
            SFF2(ION,ID)=EXP(FF(ION)*HKT1(ID))
            SFF3(ION,ID)=POPUL(NNEXT(ION),ID)*CHARG2(ION)*SGFF
   20    CONTINUE
c
C        this part is analogous to SGMER0 - for one depth only
C
      IMER=0
      DO 25 II=1,NLEVEL
         IF(IFWOP(II).GE.0) GO TO 25
         IMER=IMER+1
         IMRG(II)=IMER
         IIMER(IMER)=II
         IE=IEL(II)
         CH=IZ(IE)*IZ(IE)
         FRCH(IMER)=FRH*CH
         SGM0(IMER)=PH2*CH*CH
         II0=NQUANT(II-1)+1
            EX=EHB*CH*TEMP1(ID)
            DO 21 I=II0,NLMX
               FREDG(I)=FRCH(IMER)*XI2(I)
               EXI=EXP(EX*XI2(I))
               S(I)=EXI*WNHINT(I,ID)*XI3(I)
               SUM(I)=0.
   21       CONTINUE
            SUM(NLMX)=S(NLMX)
            DO 22 I=NLMX-1,II0,-1
               SUM(I)=SUM(I+1)+S(I)
   22       CONTINUE
            DO 23 I=1,II0-1
               SUM(I)=SUM(II0)
   23       CONTINUE
            SGEM=SGM0(IMER)/GMER(IMER,ID)
            DO 24 I=1,NLMX
               SGMSUM(I,IMER,ID)=SUM(I)*SGEM
   24       CONTINUE
   25 CONTINUE
C
C     initialization of the line opacity
C
      IF(NFRQ.GT.NFREQC) THEN
      DO 29 ITR=1,NTRANS
         IF(.NOT.LINE(ITR)) GO TO 29
         IF(INTMOD(ITR).EQ.0) GO TO 29
         INDXA=IABS(INDEXP(ITR))
         IJL0=IFR0(ITR)
         IJL1=IFR1(ITR)
       IF(ISPODF.GE.1) THEN
           IJL0=KFR0(ITR)
           IJL1=KFR1(ITR)
       END IF
         II=ILOW(ITR)
         JJ=IUP(ITR)
         IF(INDXA.LT.2.OR.INDXA.GT.4) THEN
         CALL LINPRO(ITR,ID,PRF)
         DO 28 IJ=IJL0,IJL1
            PRFLIN(ID,IJ)=real(PRF(IJ-IJL0+1))
   28    CONTINUE
         ENDIF
   29 CONTINUE
      GG=G(II)/G(JJ)
      IF(IFWOP(JJ).GE.0) THEN
         PI=POPUL(II,ID)*WOP(JJ,ID)
         PJ=POPUL(JJ,ID)*WOP(II,ID)*GG
       ELSE
         PI=POPUL(II,ID)
         PJ=POPUL(JJ,ID)*WOP(II,ID)*G(II)/GMER(IMRG(JJ),ID)
      END IF
      ABTRA(ITR,ID)=PI
      EMTRA(ITR,ID)=PJ*EXP(FR0(ITR)*HKT1(ID))
      IF(LASER) THEN
        qtt=0.
      if(pi.ne.pj) QTT=PJ/(PI-PJ)*(EXP(FR0(ITR)*HKT1(ID))-UN)
      IF(QTT.LT.0. .OR. QTT.GT.QTLAS) THEN
          ABTRA(ITR,ID)=0.
          EMTRA(ITR,ID)=0.
        END IF
      END IF
      END IF
C
C     ---------------------------------------------------------
C     
C     loop over frequency points
C
      ICALL=1
      DO 200 IJ=1,NFRQ
         if(icompt.gt.0) ELSCAT(ID)=ELEC(ID)*SIGEC(IJ)
         ABSO(IJ)=ELSCAT(ID)
         EMIS(IJ)=0.
         SCAT(IJ)=ELSCAT(ID)
C
C     basic frequency- and depth-dependent quantities
C
         FR=FREQ(IJ)
         FRINV=UN/FR
         FR3INV=FRINV*FRINV*FRINV
         XKF(ID)=EXP(-HKT1(ID)*FR)
         XKF1(ID)=UN-XKF(ID)
         XKFB(ID)=XKF(ID)*BNUE(IJ)
C
C ********  1. bound-free contribution
C
         DO 30 IBFT=1,NTRANC
            ITR=ITRBF(IBFT)
            II=ILOW(ITR)
            JJ=IUP(ITR)
          if(ifdiel.eq.0) then
               SG=CROSS(IBFT,IJ)
           else
               SG=CROSSD(IBFT,IJ,ID)
          endif
            IF(IFWOP(II).LT.0) THEN
               IMER=IMRG(II)
               CALL SGMER1(FRINV,FR3INV,IMER,ID,SGME1)
               SGMG(IMER,ID)=SGME1
               SG=SGME1
            END IF
            IF(SG.LE.0.) GO TO 30
            IF(MCDW(ITR).GT.0) THEN
               IZZ=IZ(IEL(II))
               CALL DWNFR1(FR,FR0(ITR),ID,IZZ,DW1)
               DWF1(MCDW(ITR),ID)=DW1
               SG=SG*DW1
            END IF
            EMISBF=SG*EMTRA(ITR,ID)
            ABSO(IJ)=ABSO(IJ)+SG*ABTRA(ITR,ID)
            EMIS(IJ)=EMIS(IJ)+EMISBF
   30    CONTINUE
C
C ******** 2. free-free contribution
C
         DO 40 ION=1,NION
            IT=ITRA(NNEXT(ION),NNEXT(ION))
C
C           hydrogenic with Gaunt factor = 1
C
            IF(IT.EQ.1) THEN
               SF1=SFF3(ION,ID)*FR3INV
               SF2=SFF2(ION,ID)
               IF(FR.LT.FF(ION)) SF2=UN/XKF(ID)
               ABSOFF=SF1*SF2
               ABSO(IJ)=ABSO(IJ)+ABSOFF
               EMIS(IJ)=EMIS(IJ)+ABSOFF
C
C            hydrogenic with exact Gaunt factor 
C
             ELSE IF(IT.EQ.2) THEN
               SF1=SFF3(ION,ID)*FR3INV
               SF2=SFF2(ION,ID)
               IF(FR.LT.FF(ION)) SF2=UN/XKF(ID)
               X=C14*CHARG2(ION)/FR
               SF2=SF2-UN+GFREE1(ID,X)
               ABSOFF=SF1*SF2
               ABSO(IJ)=ABSO(IJ)+ABSOFF
               EMIS(IJ)=EMIS(IJ)+ABSOFF
C
C            H minus free-free opacity
C
             ELSE IF(IT.EQ.3) THEN
               ABSOFF=(CFF1+CFFT(ID)*FRINV)*CFFN(ID)*FRINV
               ABSO(IJ)=ABSO(IJ)+ABSOFF
               EMIS(IJ)=EMIS(IJ)+ABSOFF
C
C            special evaluation of the cross-section
C
             ELSE IF(IT.LT.0) THEN
               ABSOFF=FFCROS(ION,IT,TEMP(ID),FR)*
     *                POPUL(NNEXT(ION),ID)*ELEC(ID)
               ABSO(IJ)=ABSO(IJ)+ABSOFF
               EMIS(IJ)=EMIS(IJ)+ABSOFF
            END IF
   40    CONTINUE
C
C     ********  3. - additional continuum opacity (OPADD)
C
         IF(IOPADD.NE.0) THEN
            CALL OPADD(0,ICALL,IJ,ID)
            ABSO(IJ)=ABSO(IJ)+ABAD
            EMIS(IJ)=EMIS(IJ)+EMAD
            SCAT(IJ)=SCAT(IJ)+SCAD
         END IF
C
C ********  4. - opacity and emissivity in lines
C      
         IF(ISPODF.EQ.0) THEN
         IF(IJLIN(IJ).GT.0) THEN
C
C        the "primary" line at the given frequency
C
            ITR=IJLIN(IJ)
            SG=PRFLIN(ID,IJ)
            ABSO(IJ)=ABSO(IJ)+SG*ABTRA(ITR,ID)
            EMIS(IJ)=EMIS(IJ)+SG*EMTRA(ITR,ID)
         ENDIF
         IF(NLINES(IJ).LE.0) GO TO 110
C
C        the "overlapping" lines at the given frequency
C
         DO 100 ILINT=1,NLINES(IJ)
            ITR=ITRLIN(ILINT,IJ)
          if(linexp(itr)) goto 100
            IJ0=IFR0(ITR)
            DO 60 IJT=IJ0,IFR1(ITR)
               IF(FREQ(IJT).LE.FR) THEN
                  IJ0=IJT
                  GO TO 70
               END IF
   60       CONTINUE
   70       IJ1=IJ0-1
            X=UN/(FREQ(IJ1)-FREQ(IJ0))
            A1=(FR-FREQ(IJ0))*X
            X=UN/(FREQ(IJ1)-FREQ(IJ0))
            A2=(FREQ(IJ1)-FR)*X
            SG=A1*PRFLIN(ID,IJ1)+A2*PRFLIN(ID,IJ0)
            ABSO(IJ)=ABSO(IJ)+SG*ABTRA(ITR,ID)
            EMIS(IJ)=EMIS(IJ)+SG*EMTRA(ITR,ID)
  100    CONTINUE
  110    CONTINUE
C
C     Opacity sampling option
C
      ELSE
      IF(NLINES(IJ).LE.0) GO TO 400
      DO 300 ILINT=1,NLINES(IJ)
        ITR=ITRLIN(ILINT,IJ)
c       IF(LINEXP(ITR)) GO TO 300
      KJ=IJ-IFR0(ITR)+KFR0(ITR)
      INDXPA=IABS(INDEXP(ITR))
      IF(INDXPA.NE.3 .AND. INDXPA.NE.4) THEN
        DO ID=1,ND
            SG=PRFLIN(ID,KJ)
            ABSO(IJ)=ABSO(IJ)+SG*ABTRA(ITR,ID)
            EMIS(IJ)=EMIS(IJ)+SG*EMTRA(ITR,ID)
        END DO
      ELSE
        DO ID=1,ND
          KJD=JIDI(ID)
c      SG=EXP(XJID(ID)*SIGFE(KJD,KJ)+(UN-XJID(ID))*SIGFE(KJD+1,KJ))
            ABSO(IJ)=ABSO(IJ)+SG*ABTRA(ITR,ID)
            EMIS(IJ)=EMIS(IJ)+SG*EMTRA(ITR,ID)
        END DO
      ENDIF
  300 CONTINUE
  400 CONTINUE
      ENDIF
C
C        ----------------------------
C        total opacity and emissivity
C        ----------------------------
C
         ABSO(IJ)=ABSO(IJ)-EMIS(IJ)*XKF(ID)
         EMIS(IJ)=EMIS(IJ)*XKFB(ID)
c
c        contribution from precalculated opacity table
c
c
c        H2-H2 CIA opacity
c  
         if(ifcia.gt.0) then
c           call cia_sub(t,anh2(id),freq(ij),oph2)
            abso(ij)=abso(ij)+oph2
            emis(ij)=emis(ij)+oph2*xkfb(id)/xkf1(id)
         end if
         
  200 CONTINUE
      RETURN
      END
C
C
C     ****************************************************************
C
C
      SUBROUTINE SRTFRQ
C     =================
C
C     Sort the frequency sets, and assign to each frequency
C     a list of contributing transitions
C     Select final frequency set.
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'ATOMIC.FOR'
      INCLUDE 'MODELQ.FOR'
      PARAMETER (SIXTH=UN/6.,FTH=4./3.,V0X=4.D-4,VCX=10.*V0X)
      DIMENSION SX(500)
C     DIMENSION SGZ(MTRANS)
C
      if(ioptab.lt.0) return
      if(ispodf.ge.1) return
C
C     Sort frequencies and assign primary line
C
      CALL INDEXX(NFREQ,FREQ,NLINES)
      DO 10 IJ=1,NFREQ
         KIJ(NLINES(IJ))=NFREQ-IJ+1
   10 CONTINUE
      DO 20 IT=1,NTRANS
         IF(LINEXP(IT)) GOTO 20
         KFR0(IT)=KIJ(IFR0(IT))
         KFR1(IT)=KIJ(IFR1(IT))
         DO 25 IJ=IFR0(IT),IFR1(IT)
            IJLIN(IJ)=IT
   25    CONTINUE
   20 CONTINUE
      DO 30 IJ=1,NFREQ
         JIK(KIJ(IJ))=IJ
   30 CONTINUE
      JK1=JIK(1)
      IF(IJLIN(JK1).NE.0)
     *  CALL QUIT(' Largest freq. is a line freq. - (SRTFRQ)',
     *             JK1,IJLIN(JK1))
      JK1=JIK(NFREQ)
      IF(IJLIN(JK1).NE.0)
     *  CALL QUIT(' Smallest freq. is a line freq. - (SRTFRQ)',
     *             JK1,IJLIN(JK1))
C
C     lines or ODFs associated with each frequency
C
      NLIMAX=0
      DO 40 IJ=1,NFREQ
         NLINES(IJ)=0
         DO 50 IT=1,NTRANS
            IF(LINEXP(IT)) GOTO 50
            IF(KIJ(IJ).LT.KFR0(IT)) GOTO 50
            IF(KIJ(IJ).GT.KFR1(IT)) GOTO 50
            IF(IJLIN(IJ).EQ.IT) GOTO 50
            NLINES(IJ)=NLINES(IJ)+1
            IF(NLINES(IJ).GT.MITJ)
     *      CALL QUIT('Too many overlappins-nlines(ij).gt.mitj',
     *      nlines(ij),mitj)
            ITRLIN(NLINES(IJ),IJ)=int2(IT)
   50    CONTINUE
         IF(NLINES(IJ).GT.NLIMAX) NLIMAX=NLINES(IJ)
   40 CONTINUE
c     WRITE(6,611) NLIMAX
c 611 FORMAT(' MAXIMUM NUMBER OF OVERLAPPING TRANSITIONS:  ',I3/)
C
C     Select final set of frequencies:
C        IJX = 1 : included frequency
C        IJX =-1 : rejected frequency
C        IJX = 0 : used for rates, but no contribution of primary
C                   transition to opacity
C
      NPPX=NFREQ-NFREQC
      DO 310 IT=1,NTRANS
         IF(LINEXP(IT)) GOTO 310
         IF(ABS(INDEXP(IT)).NE.3) GOTO 310
         IF(PROF(IFR0(IT)+1).GT.PROF(IFR1(IT)-1)) THEN
           DO 315 IJ=IFR0(IT)+5,IFR1(IT)-1
              IJX(IJ)=-1
              NPPX=NPPX-1
  315      CONTINUE
         ELSE
           DO 316 IJ=IFR0(IT)+1,IFR1(IT)-5
              IJX(IJ)=-1
              NPPX=NPPX-1
  316      CONTINUE
         ENDIF
  310 CONTINUE
      ISX=0
      DO 320 IJ=1,NFREQ
        ISX=ISX-1
        IF(ISX.GT.0) GOTO 320
        IJP=JIK(IJ)
        DX0=0.
        IF(IJX(IJP).EQ.1) GOTO 320
        IF(PROF(IJP).EQ.0.) GOTO 320
        DX0=V0X*FREQ(IJP)
        DNUX=ABS(FREQ(JIK(IJ-1))-FREQ(IJP))
        IF(DNUX.GT.DX0) THEN
          IJX(IJP)=1
          NPPX=NPPX+1
        ELSE
          NPX=0
          DO WHILE (DNUX.LT.DX0 .AND. IJX(JIK(IJ+NPX)).EQ.-1)
            ITRX=IJLIN(JIK(IJ+NPX))
            PSX0=PROF(IFR0(ITRX+1))
            IF(PSX0.GT.0.) THEN
              SX0=PROF(JIK(IJ+NPX))/PSX0
              SX(NPX+1)=PROF(JIK(IJ+NPX))/PROF(IJP)*SX0
            ELSE
              SX(NPX+1)=0.
            ENDIF
            NPX=NPX+1
            DNUX=ABS(FREQ(JIK(IJ-1))-FREQ(JIK(IJ+NPX)))
          END DO
          IF(NPX.EQ.1) THEN
            IJX(IJP)=1
            NPPX=NPPX+1
          ELSE
            SXX=-1.
            DO 325 IPX=1,NPX
              IF(SX(IPX).GT.SXX) THEN
                SXX=SX(IPX)
                ISX=IPX
              ENDIF
  325       CONTINUE
            IJX(JIK(IJ+ISX))=1
            NPPX=NPPX+1
          ENDIF
        ENDIF
  320 CONTINUE
      DO 330 IJ=1,NFREQ
        IJP=JIK(IJ)
        IF(IJP.GT.NFREQC) GOTO 330
        IF(IJX(IJP).EQ.1) THEN
          NPPX=NPPX+1
          GOTO 330
        ENDIF
        DX0=VCX*FREQ(IJP)
        NIXA=0
        DO WHILE (IJX(JIK(IJ-NIXA)).NE.1)
           NIXA=NIXA+1
        ENDDO
        NIXB=0
        DO WHILE (IJX(JIK(IJ+NIXB)).NE.1)
           NIXB=NIXB+1
        ENDDO
        DNUXA=ABS(FREQ(JIK(IJ-NIXA))-FREQ(IJP))
        DNUXB=ABS(FREQ(JIK(IJ+NIXB))-FREQ(IJP))
        IF(DNUXA.GT.DX0 .AND. DNUXB.GT.DX0) THEN
          IJX(IJP)=1
          NPPX=NPPX+1
        ELSE
          IJX(IJP)=-1
        ENDIF
  330 CONTINUE
c
c     correction
c
      if(icompt.eq.0) then
      do ij=1,nfreqc
         ijx(ij)=1
      end do
      do ije=1,nfreqe
         ijx(ijfr(ije))=1
      end do
      end if
C
C     weights
C
      DO 100 IJ=1,NFREQ
         W(IJ)=0.
         KJ0=KIJ(IJ)
         IF(IJX(JIK(KJ0)).EQ.-1) GOTO 100
         IF(KJ0.GE.2 .AND. KJ0.LT.NFREQ) THEN
           IK1=KJ0-1
           DO WHILE (IJX(JIK(IK1)).EQ.-1)
             IK1=IK1-1
           ENDDO
           IK2=KJ0+1
           DO WHILE (IJX(JIK(IK2)).EQ.-1)
             IK2=IK2+1
           ENDDO
           W(IJ)=HALF*ABS(FREQ(JIK(IK1))-FREQ(JIK(IK2)))
         ELSE IF(KJ0.EQ.1) THEN
           W(IJ)=HALF*ABS(FREQ(JIK(KJ0))-FREQ(JIK(KJ0+1)))
         ELSE IF(KJ0.EQ.NFREQ) THEN
           W(IJ)=HALF*ABS(FREQ(JIK(KJ0-1))-FREQ(JIK(KJ0)))
         ENDIF
  100 CONTINUE
C
C     Correction for Simpson weights
C
      JK1=JIK(1)
      DO 120 IJ=2,NFREQ,2
         JK2=JIK(IJ)
         JK3=JIK(IJ+1)
         IF(IJLIN(JK2).NE.0 .OR. IJLIN(JK3).NE.0) GOTO 130
         IF(WCH(JK2).NE.0.) GOTO 130
         W(JK1)=W(JK1)-SIXTH*W(JK2)
         W(JK3)=W(JK3)-SIXTH*W(JK2)
         W(JK2)=W(JK2)*FTH
         JK1=JK3
  120 CONTINUE
  130 JK1=JIK(NFREQ)
      DO 140 IJ=NFREQ-1,1,-2
         JK2=JIK(IJ)
         JK3=JIK(IJ-1)
         IF(IJLIN(JK2).NE.0 .OR. IJLIN(JK3).NE.0) GOTO 150
         IF(WCH(JK2).NE.0.) GOTO 150
         W(JK1)=W(JK1)-SIXTH*W(JK2)
         W(JK3)=W(JK3)-SIXTH*W(JK2)
         W(JK2)=W(JK2)*FTH
         JK1=JK3
  140 CONTINUE
  150 CONTINUE
C
C     check accuracy of weights for integration
C
c 190 Z0=0.
      Z0=0.
      Z1=0.
      Z2=0.
      ZH=0.
      T1=TEFF
      T2=TWO*TEFF
      T3=HALF*TEFF
      X1=HK/T1
      X2=HK/T2
      X3=HK/T3
      DO 200 IJ=1,NFREQ
         Z0=Z0+W(IJ)
         X15=FREQ(IJ)*1.D-15
         BNZ=BN*X15*X15*X15
         FX1=FREQ(IJ)*X1
         IF(FX1.GT.100.) GOTO 200
         Z1=Z1+W(IJ)*BNZ/(EXP(FREQ(IJ)*X1)-1)
         Z2=Z2+W(IJ)*BNZ/(EXP(FREQ(IJ)*X2)-1)
         ZH=ZH+W(IJ)*BNZ/(EXP(FREQ(IJ)*X3)-1)
  200 CONTINUE
      T1S=SQRT(SQRT(0.25*Z1/SIG4P))
      T1ER=T1S/T1-UN
      T2S=SQRT(SQRT(0.25*Z2/SIG4P))
      T2ER=T2S/T2-UN
      T3S=SQRT(SQRT(0.25*ZH/SIG4P))
      T3ER=T3S/T3-UN
      JK1=JIK(1)
      JK2=JIK(NFREQ)
      Z00=FREQ(JK1)-FREQ(JK2)
c     WRITE(6,601) FREQ(JK1),FREQ(JK2),Z00,Z0,T3,T3ER,T1,T1ER,T2,T2ER
c 601 FORMAT(/' ACCURACY OF INTEGRATIONS:',/,
c    * ' Interval:',1p4e16.8,/,
c    * 15x,' Planck functions:',9x,0pf12.0,4x,1pe12.4,/,
c    * 42x,0pf12.0,4x,1pe12.4,/,42x,0pf12.0,4x,1pe12.4,/)
c     WRITE(6,602) NFREQ,NPPX
c 602 FORMAT(' TOTAL NUMBER OF FREQUENCIES:',I8,/,
c    *       ' SELECTED FREQUENCIES:       ',I8)
C
      RETURN
      END
C
C
C     ****************************************************************
C
C
      SUBROUTINE INIFRC(IALIEX)
C     =========================
C
C     Setup continuum frequencies, including frequencies around
C       ionization limits
C          IALIEX=0 :  setup frequencies, all ALI
C          IALIEX=1 :  change IJALI for explicit frequencies
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'ATOMIC.FOR'
      INCLUDE 'MODELQ.FOR'
      INCLUDE 'ODFPAR.FOR'
      PARAMETER (THIRD=UN/3.,FTH=4./3.,TENLG = 2.302585093)
c      parameter (dfedg=0.000001)
      DIMENSION FRLEV(MLEVEL),IENS(MLEVEL),IJFL(MLEVEL),IJXCO(MFREQC),
     *          FREQCO(MFREQC),WCO(MFREQC),WCHCO(MFREQC)
C
      SAVE IJFL
      if(ioptab.lt.0) return
C
      dfedg=0.000001
      if(icompt.gt.0) dfedg=0.01
      DFEDG=0.01
      IF(IALIEX.EQ.1) THEN
C
        if(nffix.lt.0) then
           do ij=1,nfreqc
              ijali(ij)=0
           end do
           return
        end if
c
        IF(NFFIX.EQ.2) RETURN
        DO 10 IT=1,NTRANS
          IF(LINE(IT)) GOTO 10
          IF(IFC1(IT).EQ.0) GOTO 10
          IJFL0=IJFL(ILOW(IT))+1
          IF(IFC1(IT).LT.100) THEN
            DO 11 IJ=IFC0(IT),IFC1(IT)
              IJFLS=IJFL0-IJ
              IF(IJFLS.GE.1) THEN
                IJALI(IJFLS)=0
                IJX(IJFLS)=1
              ENDIF
   11       CONTINUE
          ELSE
            DO 12 IJ=IJFL0,1,-1
              IJALI(IJ)=0
              IJX(IJ)=1
   12       CONTINUE
          ENDIF
   10   CONTINUE
C
        IF(ICOMPT.GT.0.AND.FRLCOM.GT.0) THEN
      DO IJ=1,NFREQ
         IF(FREQ(IJ).GT.FRLCOM) THEN
              IJALI(IJ)=0
              IJX(IJ)=1
         END IF
      END DO
      END IF
        RETURN
      ENDIF
C
      NEND=NFTAIL
      DIVEND=DFTAIL
      NJC=NFREQC/5
      DNX=UN-UN/FLOAT(NJC)
      CALL INDEXX(NLEVEL,ENION,IENS)
      DO 50 IL=1,NLEVEL
         ILS=IENS(NLEVEL-IL+1)
         FRLEV(IL)=ENION(ILS)/H
         IJFL(IL)=0
   50 CONTINUE
c      NFREQC=NEND+1
      IF(FRCMAX.LE.0.) THEN
         FREQCO(1)=8.D11*TEFF
      ELSE
         FREQCO(1)=FRCMAX
      ENDIF
      IL0=1
      IF(FREQCO(1).LT.CFRMAX*FRLEV(IL0) .AND. CFRMAX.GT.UN) THEN
        FREQCO(1)=CFRMAX*FRLEV(IL0)
      GO TO 60
      ENDIF
      DO WHILE (FREQCO(1).LT.FRLEV(IL0) .AND. IL0.LT.NLEVEL)
        ILS=IENS(NLEVEL-IL0+1)
      ILN=NNEXT(IEL(ILS))
      ITR0=ITRA(ILS,ILN)
      INDEXP(ITR0)=0
      IFR0(ITR0)=0
      IFR1(ITR0)=0
        IF(FRLEV(IL0).GT.0.) WRITE(10,159) IL0,FRLEV(IL0),ILS,ILN,ITR0
  159   FORMAT(' Edge at frequency larger than FRCMAX',I5,E12.4,3I7)
        IL0=IL0+1
      END DO
   60 CONTINUE
C 
      if(nftail.le.0) then 
      if(nftail.eq.-1.or.dftail.eq.0.) then
      NEND1=NFREQC
      FREQCO(NEND1)=FRCMIN
      IJXCO(1)=1
      IJXCO(NEND1)=1
      XEND=UN/FLOAT(NEND1-1)
      D121=(FREQCO(1)/FREQCO(NEND1))**XEND
      DO IJ=2,NEND1-1
         FREQCO(IJ)=FREQCO(IJ-1)/D121
         IJXCO(IJ)=2
      END DO
      D121=THIRD*(FREQCO(1)-FREQCO(2))
      DO IJ=2,NEND1-1,2
         WCO(IJ)=4.*D121
         WCO(IJ-1)=WCO(IJ-1)+D121
         WCO(IJ+1)=WCO(IJ+1)+D121
         WCHCO(IJ)=0.
      END DO
      else
c      
      NEND1=-NFTAIL
      FREQCO(NEND1)=DFTAIL*FREQCO(1)
      IJXCO(1)=1
      IJXCO(NEND1)=1
      XEND=UN/FLOAT(NEND1-1)
      D121=(FREQCO(1)/FREQCO(NEND1))**XEND
      DO IJ=2,NEND1-1
         FREQCO(IJ)=FREQCO(IJ-1)/D121
         IJXCO(IJ)=2
      END DO
      D121=THIRD*(FREQCO(1)-FREQCO(2))
      DO IJ=2,NEND1-1,2
         WCO(IJ)=4.*D121
         WCO(IJ-1)=WCO(IJ-1)+D121
         WCO(IJ+1)=WCO(IJ+1)+D121
         WCHCO(IJ)=0.
      END DO
C
      NEND2=NFREQC
c      FREQCO(NEND2)=FRCMIN
      IJXCO(NEND1)=1
      IJXCO(NEND1+1)=1
      XEND=UN/FLOAT(NEND2-NEND1)
      d121=(freqco(nend1)/frcmin)**xend
      DO IJ=NEND1+1,NEND2
         freqco(ij)=freqco(ij-1)/d121
         IJXCO(IJ)=2
      END DO
      D121=THIRD*(FREQCO(NEND1)-FREQCO(NEND1+1))
      DO IJ=NEND1+1,NEND2-1,2
         WCO(IJ)=4.*D121
         WCO(IJ-1)=WCO(IJ-1)+D121
         WCO(IJ+1)=WCO(IJ+1)+D121
         WCHCO(IJ)=0.
      END DO
      end if
C
      else
C
      NEND=NFTAIL
      DIVEND=DFTAIL
      NJC=NFREQC/5
      DNX=UN-UN/FLOAT(NJC)
      NFREQC=NEND+1
   
      FREQCO(NEND)=(un+dfedg)*FRLEV(il0)
      FREQCO(NEND+1)=(un-dfedg)*FRLEV(il0)
      NEND1=NEND/2+1
c      if(icompt.gt.0) nend1=nend
      XEND=UN/FLOAT(NEND1-1)
      FREQCO(NEND1)=FREQCO(1)-(UN-DIVEND)*(FREQCO(1)-FREQCO(NEND))
      IJXCO(NEND+1)=1
C
C     high-frequency tail of the spectrum - a two-part Simpson integration
C
C     1st part - from the highest frequency FRCMAX to a division freq.
C      
      IJXCO(1)=1
      IJXCO(NEND1)=1
      D121=XEND*(FREQCO(1)-FREQCO(NEND1))
      if(icompt.gt.0) d121=(freqco(1)/freqco(nend1))**xend
      DO IJ=2,NEND1-1
         FREQCO(IJ)=FREQCO(IJ-1)-D121
         if(icompt.gt.0) freqco(ij)=freqco(ij-1)/d121
         IJXCO(IJ)=2
      END DO
      D121=THIRD*(FREQCO(1)-FREQCO(2))
      DO IJ=2,NEND1-1,2
         WCO(IJ)=4.*D121
         WCO(IJ-1)=WCO(IJ-1)+D121
         WCO(IJ+1)=WCO(IJ+1)+D121
         WCHCO(IJ)=0.
      END DO
C
C     2nd part - from the division freq to the first discontinuity
C      
      if(nend1.lt.nend) then     
      IJXCO(NEND)=1
      IJXCO(NEND+1)=1
      D121=XEND*(FREQCO(NEND1)-FREQCO(NEND))
      if(icompt.gt.0) d121=(freqco(nend1)/freqco(nend))**xend
      DO IJ=NEND1+1,NEND-1
         FREQCO(IJ)=FREQCO(IJ-1)-D121
         if(icompt.gt.0) freqco(ij)=freqco(ij-1)/d121
         IJXCO(IJ)=2
      END DO
      D121=THIRD*(FREQCO(NEND1)-FREQCO(NEND1+1))
      DO IJ=NEND1+1,NEND-1,2
         WCO(IJ)=4.*D121
         WCO(IJ-1)=WCO(IJ-1)+D121
         WCO(IJ+1)=WCO(IJ+1)+D121
         WCHCO(IJ)=0.
      END DO
      end if
C 
C     the 1st discontinuity - the one with the highest frequency 
C     
      HAEND=HALF*(FREQCO(NEND)-FREQCO(NEND+1))
      XEND=UN/FLOAT(NEND-1)
      WCO(NEND)=WCO(NEND)+HAEND
      WCO(NEND+1)=WCO(NEND+1)+HAEND
      WCHCO(1)=0.
      WCHCO(NEND)=HAEND
      ILS=IENS(NLEVEL)
      IJFL(ILS)=NEND
      IL0=NLEVEL
      IF(FRCMIN.LE.0.) FRCMIN=1.D12
      FRCLST=FRLEV(IL0)
      DO WHILE(FRCLST.LT.FRCMIN)
         IF(FRLEV(IL0).GT.0.) WRITE(10,160) IL0,FRLEV(IL0)
  160    FORMAT(' Edge at frequency smaller than 1.d12',i5,e12.4)
         IL0=IL0-1
         FRCLST=FRLEV(IL0)
      END DO
      IL0=2
C
  100 FRC0=DNX*FREQCO(NFREQC)
      IF(FRC0.LT.FRCLST) THEN
        NFREQC=NFREQC+2
        FREQCO(NFREQC-1)=(un+dfedg)*FRCLST
        FREQCO(NFREQC)=(un-dfedg)*FRCLST
        IJXCO(NFREQC-1)=1
        IJXCO(NFREQC)=1
        WCO(NFREQC)=WCO(NFREQC)+HALF*(FREQCO(NFREQC-1)-FREQCO(NFREQC))
        WCO(NFREQC-1)=WCO(NFREQC-1)+
     *                HALF*(FREQCO(NFREQC-2)-FREQCO(NFREQC))
        WCO(NFREQC-2)=WCO(NFREQC-2)+
     *                HALF*(FREQCO(NFREQC-2)-FREQCO(NFREQC-1))
        WCHCO(NFREQC-1)=HALF*(FREQCO(NFREQC-1)-FREQCO(NFREQC))
        WCHCO(NFREQC-2)=HALF*(FREQCO(NFREQC-2)-FREQCO(NFREQC-1))
        ILS=IENS(NLEVEL-IL0+1)
        IJFL(ILS)=NFREQC-1
        IF(IL0.LT.NLEVEL) THEN
          DO 110 IL=IL0+1,NLEVEL
            ILS=IENS(NLEVEL-IL+1)
            IJFL(ILS)=NFREQC-1
            IF(FRLEV(IL).LT.FRCMIN) IJFL(ILS)=0
  110     CONTINUE 
        ENDIF
        D121=XEND*(FREQCO(NFREQC)-FRCMIN)
        DO 120 IJ=NFREQC+1,NFREQC+NEND-1
          FREQCO(IJ)=FREQCO(IJ-1)-D121
          IJXCO(IJ)=2
          WCHCO(IJ)=0.
  120   CONTINUE
        IJXCO(NFREQC+NEND-1)=1
        DO 130 IJ=NFREQC+1,NFREQC+NEND-2,2
           WCO(IJ)=FTH*D121
           WCO(IJ-1)=WCO(IJ-1)+THIRD*D121
           WCO(IJ+1)=WCO(IJ+1)+THIRD*D121
  130   CONTINUE
        WCHCO(NFREQC)=THIRD*D121
        NFREQC=NFREQC+NEND-1
        GO TO 200
      ENDIF
      DF0=FRLEV(IL0)+0.1*(FREQCO(NFREQC)-FRC0)
      FRTL=(un+dfedg)*FRLEV(IL0)
      IF(FRC0.GT.DF0) THEN
        NFREQC=NFREQC+1
        FREQCO(NFREQC)=FRC0
        IJXCO(NFREQC)=2
        WCO(NFREQC)=WCO(NFREQC)+HALF*(FREQCO(NFREQC-1)-FREQCO(NFREQC))
        WCO(NFREQC-1)=WCO(NFREQC-1)+
     *                HALF*(FREQCO(NFREQC-1)-FREQCO(NFREQC))
        WCHCO(NFREQC-1)=HALF*(FREQCO(NFREQC-1)-FREQCO(NFREQC))
      ELSE IF(FRTL.LT.FREQCO(NFREQC)) THEN
        NFREQC=NFREQC+2
        FREQCO(NFREQC-1)=FRTL
        FREQCO(NFREQC)=(un-dfedg)*FRLEV(IL0)
        IJXCO(NFREQC-1)=1
        IJXCO(NFREQC)=1
        WCO(NFREQC)=WCO(NFREQC)+HALF*(FREQCO(NFREQC-1)-FREQCO(NFREQC))
        WCO(NFREQC-1)=WCO(NFREQC-1)+
     *                HALF*(FREQCO(NFREQC-2)-FREQCO(NFREQC))
        WCO(NFREQC-2)=WCO(NFREQC-2)+
     *                HALF*(FREQCO(NFREQC-2)-FREQCO(NFREQC-1))
        WCHCO(NFREQC-1)=HALF*(FREQCO(NFREQC-1)-FREQCO(NFREQC))
        WCHCO(NFREQC-2)=HALF*(FREQCO(NFREQC-2)-FREQCO(NFREQC-1))
        ILS=IENS(NLEVEL-IL0+1)
        IJFL(ILS)=NFREQC-1
        IL0=IL0+1
      ELSE
        ILS=IENS(NLEVEL-IL0+1)
        IJFL(ILS)=NFREQC-1
        IL0=IL0+1
      ENDIF
      GO TO 100
      end if
C       
  200 SUMWC=0.
C
      IF(NFREAD.GT.0) THEN
      DO IJ=NFREQ,1,-1 
         FREQ(IJ+NFREQC)=FREQ(IJ)
         W(IJ+NFREQC)=W(IJ)
         PROF(IJ+NFREQC)=PROF(IJ)    
         IJALI(IJ+NFREQC)=IJALI(IJ)
      END DO
      DO IJ=1,NFREQC
         FREQ(IJ)=FREQCO(IJ)
         W(IJ)=WCO(IJ)
         WCH(IJ)=WCHCO(IJ)
         IJALI(IJ)=1
         IJX(IJ)=IJXCO(IJ)
         PROF(IJ)=0.
      END DO
      DO 320 ITR=1,NTRANS
         IF(.NOT.LINE(ITR).OR.INDEXP(ITR).EQ.0) GO TO 320
         IF(IFR0(ITR).GT.0) IFR0(ITR)=IFR0(ITR)+NFREQC
         IF(IFR1(ITR).GT.0) IFR1(ITR)=IFR1(ITR)+NFREQC
  320 CONTINUE
      NFREQ=NFREQ+NFREQC
      END IF
C
C     determination of the first (IFR0) and the last (IFR1) frequency
C     for each explicit continuum - only if they were not already read
C
      DO 340 ITR=1,NTRANS
         IF(LINE(ITR)) GO TO 340
            IF(IFR0(ITR).LE.0) IFR0(ITR)=1
            IF(IFR1(ITR).GT.0) GO TO 340
            IF1=0
            FR01=FR0(ITR)
            MODE=INDEXP(ITR)
            IF(IABS(MODE).EQ.5.OR.IABS(MODE).EQ.15) FR01=FR0PC(ITR)
            DO 330 IJ=1,NFREQC
               IF(FREQ(IJ).GE.FR01) IF1=IJ
  330       CONTINUE
            IFR1(ITR)=IF1
  340 CONTINUE
      RETURN
      END
C
C
C     ****************************************************************
C
C
      SUBROUTINE INIFRS
C     =================
C
C     Setup frequencies in opacity sampling mode
C       
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'ATOMIC.FOR'
      INCLUDE 'MODELQ.FOR'
      INCLUDE 'ODFPAR.FOR'
      common/nflmax/nflmx,kijmx
      PARAMETER (THIRD=UN/3.,FTH=4./3.)
      DIMENSION FRLEV(MLEVEL),IENS(MLEVEL),ITRL(MLEVEL)
      DIMENSION FRL0(MTRANS),FRL1(MTRANS),FRLC(5*MTRANS)
      DIMENSION IKC(5*MTRANS),ITKC(5*MTRANS),ITJNU(5*MTRANS)
      DIMENSION FLNU(2*MATOM+3),DLNU(2*MATOM+3),ILNU(2*MATOM+3)
      DIMENSION XMASS(30)
      DATA XMASS/ 1.008, 4.003, 6.941, 9.012,10.810,12.011,14.007,
     &           16.000,18.918,20.179,22.990,24.305,26.982,28.086,
     &           30.974,32.060,35.453,39.948,39.098,40.080,44.956,
     &           47.900,50.941,51.996,54.938,55.847,58.933,58.700,
     &           63.546,65.380/

C
      IF(TSNU.EQ.0.) TSNU=TEFF
      IF(VTNU.EQ.0.) VTNU=VTB
      IF(VTNU.LT.1.D4) VTNU=VTNU*1.D5
      FRS1=CNU1*1.D11*TSNU
      FRS2=3.28805D15/CNU2/CNU2
C
      DO IAT=1,NATOM
         CDOP=TWO*BOLK/AMASS(IAT)
         DLNU(IAT)=0.375/2.997925D10*SQRT(CDOP*TSNU+VTNU*VTNU)
         DLNU(IAT+NATOM)=20.*DLNU(IAT)
         FLNU(IAT)=DLOG(FRS1)
         FLNU(IAT+NATOM)=DLOG(FRS1)
      END DO
C
      XPNU=24.
      CDOP=TWO*BOLK/XMASS(1)/HMASS
      DLNU(2*NATOM+1)=50./2.997925D10*SQRT(CDOP*TSNU+VTNU*VTNU)
      DLNU(2*NATOM+2)=5.*DLNU(2*NATOM+1)
      FLNU(2*NATOM+1)=DLOG(FRS2)
      FLNU(2*NATOM+2)=DLOG(FRCMIN)
      NNU=2*NATOM+3
      IF(ISPODF.EQ.1 .AND. DDNU.GT.0.) THEN 
         CDOP=TWO*BOLK/AMASS(NATOM)
         IF(IELNU.GT.0) CDOP=TWO*BOLK/XMASS(IELNU)/HMASS
         DLNU(NNU)=DDNU/2.997925D10*SQRT(CDOP*TSNU+VTNU*VTNU)
         FLNU(NNU)=DLOG(FRS2)
       ELSE     
         DLNU(NNU)=DLNU(2*NATOM+2)
         FLNU(NNU)=DLOG(FRS1)
      END IF
      CALL INDEXX(NNU,DLNU,ILNU)
C
C     Store line and continua frequencies
C
      NLIC=0
      DO 10 ITR=1,NTRANS
         IJTC(ITR)=0
         INDXPA=IABS(INDEXP(ITR))
         IF(INDXPA.EQ.0) GO TO 10
         IF(INDXPA.EQ.3 .OR. INDXPA.EQ.4) GO TO 10
         IF(FR0(ITR).EQ.0.) GO TO 10
         IF(LINE(ITR)) THEN
            ILV0=ILOW(ITR)
            IAT=IATM(ILV0)
            ITC=ITRA(ILV0,NNEXT(IEL(ILV0)))
            IF(ITC.EQ.0) ITC=ITRA(ILV0,NNEXT(IEL(ILV0))+1)
            IF(INDXPA.NE.2) THEN
               FRLC(NLIC+1)=FR0(ITR)
               ITKC(NLIC+1)=ITR
               IJTC(ITR)=NLIC+1
               ITJNU(NLIC+1)=IAT
               FRLC(NLIC+2)=FREQ(IFR0(ITR))
               FRL0(ITR)=FRLC(NLIC+2)
               ITKC(NLIC+2)=ITR
               IF(ITC.GT.0 .AND. FREQ(IFR0(ITR)).GT.FR0(ITC) .AND.
     *         FR0(ITC).GT.FR0(ITR)) FRLC(NLIC+2)=0.999999*FR0(ITC)
               ITJNU(NLIC+2)=IAT
               FRLC(NLIC+3)=FREQ(IFR1(ITR))
               FRL1(ITR)=FRLC(NLIC+3)
               ITKC(NLIC+3)=ITR
               ITJNU(NLIC+3)=2*NATOM+1
               NLIC=NLIC+3
               D0=DLOG(FRL0(ITR))-DLOG(FRL1(ITR))
               IF(D0.GT.XPNU*DLNU(IAT)) THEN
                  ITJNU(NLIC-1)=IAT+NATOM
                  FRLC(NLIC+1)=EXP(DLOG(FR0(ITR))+XPNU*DLNU(IAT))
                  ITKC(NLIC+1)=ITR
                  ITJNU(NLIC+1)=IAT
                  FRLC(NLIC+2)=EXP(DLOG(FR0(ITR))-XPNU*DLNU(IAT))
                  ITKC(NLIC+2)=ITR
                  ITJNU(NLIC+2)=IAT+NATOM
                  NLIC=NLIC+2
               END IF
             ELSE
               FRLC(NLIC+1)=0.999999*FR0(ITC)
               FRL0(ITR)=FRLC(NLIC+1)
               ITKC(NLIC+1)=ITR
               ITJNU(NLIC+1)=IAT
               FRLC(NLIC+2)=FREQ(IFR1(ITR-1))
               FRL1(ITR)=FRLC(NLIC+2)
               ITKC(NLIC+2)=ITR
               ITJNU(NLIC+2)=2*NATOM+1
               NLIC=NLIC+2
            END IF
          ELSE
            NLIC=NLIC+1
            FRLC(NLIC)=FR0(ITR)
            ITKC(NLIC)=ITR
            IJTC(ITR)=NLIC
            ITJNU(NLIC)=0
         END IF
   10 CONTINUE
C
      IKC(1)=1
      IF(NLIC.GT.1) CALL INDEXX(NLIC,FRLC,IKC)
      DO IJ=1,MFREQ
         FREQ(IJ)=0.
         W(IJ)=0.
         WCH(IJ)=0.
         NLINES(IJ)=0
      END DO
C
C     Sort continuum limits
C
      CALL INDEXX(NLEVEL,ENION,IENS)
      DO IL=1,NLEVEL
         ILS=IENS(NLEVEL-IL+1)
         FRLEV(IL)=ENION(ILS)/H
         ITRL(IL)=ITRA(ILS,NNEXT(IEL(ILS)))
      END DO
      IF(FRCMAX.LT.1.01*FRLEV(1) .AND. FRCMAX.GT.0.) THEN
         ILS=IENS(NLEVEL)
         ILN=NNEXT(IEL(ILS))
         ITR0=ITRA(ILS,ILN)
         WRITE(10,640) FRLEV(1),ILS,ILN,ITR0
  640    FORMAT(1PE12.4,3I7)
         CALL QUIT(' Edge at frequency larger than FRCMAX; ii,itr:',
     *            ils,itr0)
      END IF
C
C     Highest frequency tail
C
      IF(FRCMAX.LE.0.) FRCMAX=FRLEV(1)*CFRMAX
      IF(FRS1.GT.FRCMAX) THEN
         FRCMAX=FRS1
         NFTAIL=0
      END IF
      IF(NFTAIL.GT.0) THEN
         NFTA1=NFTAIL/2+1
         FREQ(1)=FRCMAX
         NEND=0
         IL=1
         KJ=1
         DO WHILE(FRLEV(IL).GT.FRS1)
            NEND1=NEND+NFTA1
            NEND=NEND+NFTAIL
            ITR=ITRL(IL)
            IFR0(ITR)=1
            IFR1(ITR)=NEND
            FREQ(NEND)=1.000001*FRLEV(IL)
            FREQ(NEND+1)=0.999999*FRLEV(IL)
            FREQ(NEND1)=FREQ(KJ)-(UN-DFTAIL)*(FREQ(KJ)-FREQ(NEND))
            XEND=UN/FLOAT(NFTA1-1)
            D121=XEND*(FREQ(KJ)-FREQ(NEND1))
            DO IJ=KJ+1,NEND1-1
               FREQ(IJ)=FREQ(IJ-1)-D121
            END DO
            D121=THIRD*(FREQ(KJ)-FREQ(KJ+1))
            DO IJ=KJ+1,NEND1-1,2
               W(IJ)=4.*D121
               W(IJ-1)=W(IJ-1)+D121
               W(IJ+1)=W(IJ+1)+D121
            END DO
            D121=XEND*(FREQ(NEND1)-FREQ(NEND))
            DO IJ=NEND1+1,NEND-1
               FREQ(IJ)=FREQ(IJ-1)-D121
            END DO
            D121=THIRD*(FREQ(NEND1)-FREQ(NEND1+1))
            DO IJ=NEND1+1,NEND-1,2
               W(IJ)=4.*D121
               W(IJ-1)=W(IJ-1)+D121
               W(IJ+1)=W(IJ+1)+D121
            END DO
            D121=HALF*(FREQ(NEND)-FREQ(NEND+1))
            W(NEND)=W(NEND)+D121
            W(NEND+1)=W(NEND+1)+D121
            IL=IL+1
            KJ=NEND+1
         END DO
         NEND1=NEND+NFTA1
         NEND=NEND+NFTAIL
         FREQ(NEND)=FRS1
         FREQ(NEND1)=FREQ(KJ)-(UN-DFTAIL)*(FREQ(KJ)-FREQ(NEND))
         XEND=UN/FLOAT(NFTA1-1)
         D121=XEND*(FREQ(KJ)-FREQ(NEND1))
         DO IJ=KJ+1,NEND1-1
            FREQ(IJ)=FREQ(IJ-1)-D121
         END DO
         D121=THIRD*(FREQ(KJ)-FREQ(KJ+1))
         DO IJ=KJ+1,NEND1-1,2
            W(IJ)=4.*D121
            W(IJ-1)=W(IJ-1)+D121
            W(IJ+1)=W(IJ+1)+D121
         END DO
         D121=XEND*(FREQ(NEND1)-FREQ(NEND))
         DO IJ=NEND1+1,NEND-1
            FREQ(IJ)=FREQ(IJ-1)-D121
         END DO
         D121=THIRD*(FREQ(NEND1)-FREQ(NEND1+1))
         DO IJ=NEND1+1,NEND-1,2
            W(IJ)=4.*D121
            W(IJ-1)=W(IJ-1)+D121
            W(IJ+1)=W(IJ+1)+D121
         END DO
       ELSE
         FREQ(1)=FRS1
         NEND=1
      END IF
      NFREQC=NEND
      DO IJ=1,NFREQC
         IFREQB(IJ)=IJ
      END DO
      NFRS1=NFREQC
C
C     Setup frequency points
C
      DO IT=1,NTRANS
         IFR0(IT)=0
         IFR1(IT)=0
      END DO
C
      IL=NLIC
      DO WHILE(FRLC(IKC(IL)).GT.FRS1)
         IL=IL-1
      END DO
      NFREQ=NEND
      XFRA=DLOG(FRS1)
      DO WHILE(IL.GT.0)
         ITR=ITKC(IKC(IL))
         NFS=0
         XFRB=DLOG(FRLC(IKC(IL)))
         IF(XFRA.GT.XFRB) THEN
            IKNU=ITJNU(IKC(IL))
            IDN=1
            DO WHILE(FLNU(ILNU(IDN)).GE.XFRA .AND. IDN.LT.NNU)
               IDN=IDN+1
            END DO
            DXNU=DLNU(ILNU(IDN))
            IF(IKNU.EQ.0) XFRB=DLOG(1.000001*FRLC(IKC(IL)))
            NFS=INT((XFRA-XFRB)/DXNU)+1
            XFS0=(XFRA-XFRB)/FLOAT(NFS)
            DO IJ=NFREQ+1,NFREQ+NFS
               XFR=DLOG(FREQ(IJ-1))-XFS0
               FREQ(IJ)=EXP(XFR)
            END DO
            NFREQ=NFREQ+NFS
            IF(DLOG(FR0(ITR)).EQ.XFRB) IJTC(ITR)=NFREQ
            IF(IKNU.EQ.0) THEN
               IFR0(ITR)=1
               IFR1(ITR)=NFREQ
               NFREQ=NFREQ+1
               FREQ(NFREQ)=0.999999*FRLC(IKC(IL))
               XFRB=DLOG(FREQ(NFREQ))
             ELSE IF(IKNU.LE.2*NATOM) THEN
               IF(IFR0(ITR).EQ.0) THEN
                  IFR0(ITR)=NFREQ
                  D0=DLOG(FR0(ITR))
                  IF(FLNU(IKNU).GT.D0) FLNU(IKNU)=D0
                  IF(IABS(INDEXP(ITR)).EQ.2) THEN
                     D0=DLOG(FR0(ITR-1))
                     IF(FLNU(IKNU).GT.D0) FLNU(IKNU)=D0
                  END IF
                ELSE
                  IFR1(ITR)=NFREQ
                  IF(IKNU.LE.NATOM) THEN
                     IAT=IATM(ILOW(ITR))
                     D0=DLOG(FR0(ITR))-XPNU*DLNU(IAT)
                     D1=DLOG(FRL1(ITR))
                     IF(D1.GT.D0) D0=D1
                     IF(FLNU(IKNU).GT.D0) FLNU(IKNU)=D0
                   ELSE
                     D0=DLOG(FRL1(ITR))
                     IF(FLNU(IKNU).GT.D0) FLNU(IKNU)=D0 
                  END IF
               END IF
             ELSE
               IF(IFR0(ITR).EQ.0) THEN
                  IFR0(ITR)=NFREQ
                ELSE
                  IFR1(ITR)=NFREQ
               END IF
            END IF
          ELSE IF(XFRA.EQ.XFRB) THEN
            IKNU=ITJNU(IKC(IL))
            IF(IKNU.EQ.0) THEN
               FREQ(NFREQ)=1.000001*FRLC(IKC(IL))
               FREQ(NFREQ+1)=0.999999*FRLC(IKC(IL))
               IFR0(ITR)=1
               IFR1(ITR)=NFREQ
               NFREQ=NFREQ+1
               XFRB=DLOG(FREQ(NFREQ))
             ELSE IF(IKNU.LE.2*NATOM) THEN
               IF(IFR0(ITR).EQ.0) THEN
                  IFR0(ITR)=NFREQ
                  D0=DLOG(FR0(ITR))
                  IF(FLNU(IKNU).GT.D0) FLNU(IKNU)=D0
                  IF(IABS(INDEXP(ITR)).EQ.2) THEN
                     D0=DLOG(FR0(ITR-1))
                     IF(FLNU(IKNU).GT.D0) FLNU(IKNU)=D0
                  END IF
                ELSE
                  IFR1(ITR)=NFREQ
                  IF(IKNU.LE.NATOM) THEN
                     IAT=IATM(ILOW(ITR))
                     D0=DLOG(FR0(ITR))-XPNU*DLNU(IAT)
                     D1=DLOG(FRL1(ITR))
                     IF(D1.GT.D0) D0=D1
                     IF(FLNU(IKNU).GT.D0) FLNU(IKNU)=D0
                   ELSE
                     D0=DLOG(FRL1(ITR))
                     IF(FLNU(IKNU).GT.D0) FLNU(IKNU)=D0 
                  END IF
                  IF(DLOG(FR0(ITR)).EQ.XFRB) IJTC(ITR)=NFREQ
               END IF
             ELSE
               IF(IFR0(ITR).EQ.0) THEN
                  IFR0(ITR)=NFREQ
                ELSE
                  IFR1(ITR)=NFREQ
                  IF(DLOG(FR0(ITR)).EQ.XFRB) IJTC(ITR)=NFREQ
               END IF
            END IF
         END IF
         IL=IL-1
         XFRA=XFRB
         IF(XPNU.EQ.24. .AND. FREQ(NFREQ).LT.FRS2) THEN
            XPNU=HALF*XPNU
            DO IAT=1,NATOM
               DLNU(IAT)=TWO*DLNU(IAT)
               DLNU(IAT+NATOM)=TWO*DLNU(IAT+NATOM)
            END DO
         END IF
      END DO
C
      XFRB=DLOG(FRCMIN)
      IF(XFRA.GT.XFRB) THEN
         DXNU=DLNU(NNU-1)
         NFS=INT((XFRA-XFRB)/DXNU)+1
         XFS0=(XFRA-XFRB)/FLOAT(NFS)
         DO IJ=NFREQ+1,NFREQ+NFS
            XFR=DLOG(FREQ(IJ-1))-XFS0
            FREQ(IJ)=EXP(XFR)
         END DO
         NFREQ=NFREQ+NFS
      END IF
      FREQ(NFREQ)=FRCMIN
      DO 20 ITR=1,NTRANS
         IF(LINEXP(ITR)) GO TO 20
         DO IJ=IFR0(ITR),IFR1(ITR)
            NLINES(IJ)=NLINES(IJ)+1
         END DO
   20 CONTINUE
C
C     Choose continuum frequency points in the global set
C
      FRLEV(NLEVEL+1)=FRCMIN
      IL=1
      DO WHILE(FRLEV(IL).GT.FRS1)
         IL=IL+1
      END DO
      IB0=NFRS1
      NUB=2*NATOM+1
      XFRA=DLOG(FRS1)
      DO WHILE(IL.LE.NLEVEL+1)
         IF(FRLEV(IL).LT.FRCMIN) GO TO 490
         IF(IL.GT.1 .AND. IL.LE.NLEVEL) THEN
            IF(FRLEV(IL).GE.FRLEV(IL-1)) GO TO 490
         END IF
         IF(IL.LE.NLEVEL) ITR=ITRL(IL)
         FRLV0=FRLEV(IL)
         IB1=IB0
         DO WHILE(FREQ(IB1).GT.FRLV0)
            IB1=IB1+1
            XFRB=DLOG(FREQ(IB1))
            IF(IFREQB(NFREQC).LT.IB1) THEN
               IF(NLINES(IB1).EQ.0 .AND. ISPODF.GT.1) THEN
                  NFREQC=NFREQC+1
                  IFREQB(NFREQC)=IB1
                  XFRA=XFRB
                ELSE IF((XFRA-XFRB).GT.DLNU(NUB)) THEN
                  NFREQC=NFREQC+1
                  IFREQB(NFREQC)=IB1
                  XFRA=XFRB
               END IF
            END IF
         END DO
         IF(IL.LE.NLEVEL) THEN
            IFR0(ITR)=1
            IFR1(ITR)=IB1-1
            IJTC(ITR)=IFR1(ITR)
         END IF
         IF(IFREQB(NFREQC).LT.(IB1-1)) THEN
            NFREQC=NFREQC+1 
           IFREQB(NFREQC)=IB1-1
         END IF  
         IF(IFREQB(NFREQC).LT.IB1) THEN
            NFREQC=NFREQC+1
            IFREQB(NFREQC)=IB1
         END IF
         XFRA=DLOG(FREQ(IB1))
         IB0=IB1
 490     IL=IL+1
         IF(FRLEV(IL).LT.FRS2) NUB=2*NATOM+2
      END DO
C
      IF(IFREQB(NFREQC).LT.NFREQ) THEN
         NFREQC=NFREQC+1
         IFREQB(NFREQC)=NFREQ
      END IF
C
      NFREQL=0
      XBL=DLOG(FRS1)
      NFLX=0
      nflmx=0
      kijmx=0
      DO 410 IT=1,NTRANS
         IF(LINEXP(IT)) GO TO 410
         IF(FR0(IT).LT.FRCMIN) GO TO 410
         INDXPA=ABS(INDEXP(IT))
         IF(INDXPA.GT.2 .AND. INDXPA.LE.4) GO TO 410
         IL0=ILOW(IT)
         ITC=ITRA(IL0,NNEXT(IEL(IL0)))
         IF(ITC.EQ.0) ITC=ITRA(IL0,NNEXT(IEL(IL0))+1)
         IF(IFR1(IT).LE.IFR1(ITC)) GO TO 411
         IF(IFR0(IT).LE.IFR1(ITC) .AND. ITC.GT.0) IFR0(IT)=IFR1(ITC)+1
  411    NF=IFR1(IT)-IFR0(IT)+1
         KFR0(IT)=NFREQL+1
         KFR1(IT)=NFREQL+NF
         NFREQL=NFREQL+NF
         IF(INDXPA.EQ.2) THEN
            FR02H=HALF*(FREQ(IFR0(IT))+FREQ(IFR1(IT)))
            IJTC(IT)=IFR1(IT)
            DO WHILE(FR02H.GT.FREQ(IJTC(IT)) .AND. IJTC(IT).GT.1)
               IJTC(IT)=IJTC(IT)-1
            END DO
         END IF
         if(nf.gt.nflmx) nflmx=nf
         if(kfr1(it).gt.kijmx) kijmx=kfr1(it)
         IF(NF.GT.MFREQL) THEN
            WRITE(10,*) IL0,IT,NF
            CALL QUIT('Too many frequencies in a line - nf.gt.mfreql',
     *               nf,mfreql)
         END IF
         IF(NF.GT.NFLX) NFLX=NF
         IF(KFR1(IT).GT.MFREQP)
     *      CALL QUIT('Too many cross-sections to store in PRFLIN',
     *                kfr1(it),mfreqp)
  410 CONTINUE
c
c     write(6,504) nflmx,kijmx
c 504 format(' MFREQL = ',I7/
c    *       ' MFREQP = ',I7) 
C
      DO 350 IT=1,NTRANS
         MODW=IABS(INDEXP(IT))
         IF(MODW.NE.5 .AND. MODW.NE.15) GO TO 350
         IJTC(IT)=IFR1(IT)
         FRLV0=FR0PC(IT)
         IB1=NFRS1
         DO WHILE(FREQ(IB1).GT.FRLV0)
            IB1=IB1+1
         END DO
         IFR0(IT)=1
         IFR1(IT)=IB1-1
  350 CONTINUE
C
C     Weights
C
      DO IJ=NFRS1+1,NFREQ
         D121=HALF*(FREQ(IJ-1)-FREQ(IJ))
         W(IJ-1)=W(IJ-1)+D121
         W(IJ)=W(IJ)+D121
      END DO
C
      DO IJ=1,NFREQ
         IJALI(IJ)=1
         IJX(IJ)=1
         JIK(IJ)=IJ
      END DO
      NPPX=NFREQ
C
      write(10,*) nfreq,nfreqc,nfreql,nflx
      IF(NFREQ.GT.MFREQ) THEN
         WRITE(10,1000) NFREQ
 1000    FORMAT(' Number of frequencies:',I10)
         CALL QUIT('nfreq.gt.mfreq',nfreq,mfreq)
      END IF
C
      RETURN
      END
C
C
C     ****************************************************************
C
C
      FUNCTION CROSS(IBFT,IJ)
C     =======================
C
C     Evaluation of the photoionization cross-section 
C     IBF - index ot the b-f transition
C     IJ  - frequency index
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'ATOMIC.FOR'
      INCLUDE 'MODELQ.FOR'
C
      IJ0=IJBF(IJ)
      A1=AIJBF(IJ)
      CROSS=A1*BFCS(IBFT,IJ0)+(UN-A1)*BFCS(IBFT,IJ0+1)
c
      RETURN
      END
C
C
C     ****************************************************************
C
C
      FUNCTION CROSSD(IBFT,IJ,ID)
C     ===========================
C
C     Evaluation of the photoionization cross-section 
C     IBF - index ot the b-f transition
C     IJ  - frequency index
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'ATOMIC.FOR'
      INCLUDE 'MODELQ.FOR'
C
      IJ0=IJBF(IJ)
      A1=AIJBF(IJ)
      CROSSD=A1*BFCS(IBFT,IJ0)+(UN-A1)*BFCS(IBFT,IJ0+1)
c
c     contribution from dielectronic recombination
c
      if(ifdiel.eq.0) return
      ITR=ITRBF(IBFT)
      if(idiel(itr).gt.0.and.id.gt.0) then
         i=ilow(itr)
         ion=iel(i)
         if(i.eq.nfirst(ion).and.iup(itr).eq.nnext(ion)) then
            if(freq(ij).ge.fr0(itr).and.freq(ij).le.fr0(itr)*1.1) 
     *      crossd=crossd+diesig(ion,id)
         end if
      end if
c
      RETURN
      END
C
C
C     ****************************************************************
C
C
      SUBROUTINE DIETOT
C     =================
C
C     modification of the photoionization cross-section
C     for taking into account dielectronic recombination
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'ATOMIC.FOR'
      INCLUDE 'MODELQ.FOR'
C
      do ion=1,nion
         i=nfirst(ion)
         ia=numat(iatm(i))
         io=iz(ion)
         do id=1,nd
            t=temp(id)
            xpx=dens(id)/wmm(id)/ytot(id)
            call dielrc(ia,io,t,xpx,dirt,sig0)
            diesig(ion,id)=sig0
            if(id.eq.1.or.id.eq.35.or.id.eq.nd) then
            write(99,699) ion,ia,io,id,i,nnext(ion),dirt,sig0
            end if
         end do
      end do
  699 format(6i5,1p2e12.4)
      return
      end 
C
C
C     ****************************************************************
C
C
      SUBROUTINE LINSEL
C     =================
C
C     Exclude weakest lines
C     (i.e. set them to detailed radiative balance, and
C     exclude their frequencies)
C
C     Selection of lines based on ratio of line core flux to
C     continuum flux
C
C     Non-standard parameters: STRL1 (default value 0.001)
C                              STRL2 ( 0.02)
C
C     STRL2 allows to reduce the number of frequency points in
C       "intermediate-strength" lines
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'ATOMIC.FOR'
      INCLUDE 'MODELQ.FOR'
      INCLUDE 'ODFPAR.FOR'
      INCLUDE 'ALIPAR.FOR'
      DIMENSION PRFTMP(MDEPTH)
C
      parameter (SIXTH=UN/6.,FTH=4./3.)
C
      NLSW=0
      NLSI=0
      NLSS=0
      NLSTO=0
C
      IF(ISPODF.GE.1) GO TO 190
C
      CALL OPAINI(1)
C
C     Normal lines
C
      DO 10 ITR=1,NTRANS
         IF(LINEXP(ITR)) GO TO 10
         IF(LEXP(ITR) .OR. IEL(ILOW(ITR)).EQ.IELH) THEN
            NLSTO=NLSTO+1
            NLSS=NLSS+1
            GO TO 10
         ENDIF
         MODE=IABS(INDEXP(ITR))
         IF(MODE.GE.2 .AND. MODE.LE.4) GO TO 10
         IKA=IFR0(ITR)
         IF(IKA.EQ.0) GO TO 10
         IKB=IFR1(ITR)
         CALL OPACF1(IKA)
         CALL RTEFR1(IKA)
         FLUXA=FH(IKA)*RAD1(1)
         CALL OPACF1(IKB)
         CALL RTEFR1(IKB)
         FLUXB=FH(IKB)*RAD1(1)
         IK0=(IKA+IKB)/2
         IF(MODE.EQ.2) IK0=IFR1(ITR)-1
         CALL OPACF1(IK0)
         CALL RTEFR1(IK0)
         FLUX0=FH(IK0)*RAD1(1)
         RHAB=UN-TWO*FLUX0/(FLUXA+FLUXB)
         NFK0=IKB-IKA+1
         NLSTO=NLSTO+1
         IF(ABS(RHAB).LT.STRL1) THEN
            NLSW=NLSW+1
            LINEXP(ITR)=.TRUE.
            INDEXP(ITR)=0
            DO IJ=IFR0(ITR),IFR1(ITR)
               IJX(IJ)=-1
               IJLIN(IJ)=0
            END DO
            IFR0(ITR)=0
            IFR1(ITR)=0
            KFR0(ITR)=0
            KFR1(ITR)=0
            LALI(ITR)=.FALSE.
          ELSE IF(ABS(RHAB).LT.STRL2) THEN
            NLSI=NLSI+1
            IFR0(ITR)=IKA+3
            IFR1(ITR)=IKB-3
            KFR0(ITR)=KIJ(IFR0(ITR))
            KFR1(ITR)=KIJ(IFR1(ITR))
            DO IJ=IKA,IFR0(ITR)
               IJX(IJ)=-1
               IJLIN(IJ)=0
               DO ID=1,ND
                  PRFLIN(ID,IJ)=0.
               END DO
            END DO
            IJX(IFR0(ITR))=1
            IJLIN(IFR0(ITR))=ITR
            DO IJ=IFR1(ITR),IKB
               IJX(IJ)=-1
               IJLIN(IJ)=0
               DO ID=1,ND
                  PRFLIN(ID,IJ)=0.
               END DO
            END DO
            IJX(IFR1(ITR))=1
            IJLIN(IFR1(ITR))=ITR
          ELSE
            NLSS=NLSS+1
         ENDIF
         NFK1=IFR1(ITR)-IFR0(ITR)
         IF(NFK1.GT.0) NFK1=NFK1+1
         WRITE(82,601) ITR,ILOW(ITR),IUP(ITR),
     *                 NFK0,NFK1,RHAB
  601    FORMAT(I6,2I5,2I7,1PE12.4)
   10 CONTINUE
C
C     superlines
C
      DO 20 ITR=1,NTRANS
         IF(LINEXP(ITR)) GO TO 20
         MODE=IABS(INDEXP(ITR))
         IF(MODE.NE.3 .AND. MODE.NE.4) GO TO 20
         IF(LEXP(ITR)) THEN
            NLSTO=NLSTO+1
            NLSS=NLSS+1
            GO TO 20
         ENDIF
         IKA=IFR0(ITR)
         IF(IKA.EQ.0) GO TO 20
         IKB=IFR1(ITR)
         NFK0=IFR1(ITR)-IFR0(ITR)+1
         NFK1=NFK0
         RHAB=0.
         RHABMX=0.
         PRFA=PRFLIN(1,IKA+1)
         PRFB=PRFLIN(1,IKB-1)
         IF(PRFA.GT.PRFB) THEN
            IK2=IKB-1
            DO WHILE (IK2.GT.IKA .AND. RHAB.LT.STRL1)
               CALL OPACF1(IK2)
               CALL RTEFR1(IK2)
               FLUX2=FH(IK2)*RAD1(1)
               DO ID=1,ND
                  PRFTMP(ID)=PRFLIN(ID,IK2)
                  PRFLIN(ID,IK2)=0.
               END DO
               CALL OPACF1(IK2)
               CALL RTEFR1(IK2)
               FLUX1=FH(IK2)*RAD1(1)
               RHAB=ABS(UN-FLUX2/FLUX1)
               IF(RHAB.GT.RHABMX) RHABMX=RHAB
               IK20=IK2
               IF(RHAB.LT.STRL1) THEN
                  IJX(IK2+1)=-1
                  IJLIN(IK2+1)=0
                  IK2=IK2-1
                ELSE
                  IK2=IK2+1
                  IFR1(ITR)=IK2
                  KFR1(ITR)=KIJ(IK2)
                  DO ID=1,ND
                     PRFLIN(ID,IK20)=real(PRFTMP(ID))
                  END DO
               END IF
            END DO
            NFK1=IFR1(ITR)-IFR0(ITR)+1
            IF(IK2.EQ.IKA) NFK1=0
          ELSE
            IK2=IKA+1
            DO WHILE (IK2.LT.IKB .AND. RHAB.LT.STRL1)
               CALL OPACF1(IK2)
               CALL RTEFR1(IK2)
               FLUX2=FH(IK2)*RAD1(1)
               DO ID=1,ND
                  PRFTMP(ID)=PRFLIN(ID,IK2)
                  PRFLIN(ID,IK2)=0.
               END DO
               CALL OPACF1(IK2)
               CALL RTEFR1(IK2)
               FLUX1=FH(IK2)*RAD1(1)
               RHAB=ABS(UN-FLUX2/FLUX1)
               IF(RHAB.GT.RHABMX) RHABMX=RHAB
               IK20=IK2
               IF(RHAB.LT.STRL1) THEN
                  IJX(IK2-1)=-1
                  IJLIN(IK2-1)=0
                  IK2=IK2+1
                ELSE
                  IK2=IK2-1
                  IFR0(ITR)=IK2
                  KFR0(ITR)=KIJ(IK2)
                  DO ID=1,ND
                     PRFLIN(ID,IK20)=real(PRFTMP(ID))
                  END DO
               END IF
            END DO
            NFK1=IFR1(ITR)-IFR0(ITR)+1
            IF(IK2.EQ.IKB) NFK1=0
         END IF
         IF(NFK1.EQ.0) THEN
            NLSW=NLSW+1
            LINEXP(ITR)=.TRUE.
            INDEXP(ITR)=0
            IFR0(ITR)=0
            IFR1(ITR)=0
            KFR0(ITR)=0
            KFR1(ITR)=0
            DO IJ=IKA,IKB
               IJX(IJ)=-1
               IJLIN(IJ)=0
            END DO
            LALI(ITR)=.FALSE.
          ELSE IF(NFK1.EQ.NFK0) THEN
            NLSS=NLSS+1
          ELSE
            NLSI=NLSI+1
         END IF
         NLSTO=NLSTO+1
         WRITE(82,601) ITR,ILOW(ITR),IUP(ITR),
     *                 NFK0,NFK1,RHABMX
   20 CONTINUE
C
      WRITE(6,602) NLSTO,NLSW,NLSI,NLSS
  602 FORMAT(' Total number of lines :',i8,/,
     *          ' Number of weak lines  :',i8,/,
     *          ' Intermediate lines    :',i8,/,
     *          ' Number of strong lines:',i8/)
C
C     lines or ODFs associated with each frequency
C
      NLIMAX=0
      DO 40 IJ=1,NFREQ
         NLINES(IJ)=0
         DO 50 IT=1,NTRANS
            IF(LINEXP(IT)) GO TO 50
            IF(KIJ(IJ).LT.KFR0(IT)) GO TO 50
            IF(KIJ(IJ).GT.KFR1(IT)) GO TO 50
            IF(IJLIN(IJ).EQ.IT) GO TO 50
            NLINES(IJ)=NLINES(IJ)+1
            IF(NLINES(IJ).GT.MITJ)
     *      CALL QUIT('Too many overlappins-nlines(ij).gt.mitj',
     *      nlines(ij),mitj)
            ITRLIN(NLINES(IJ),IJ)=int2(IT)
   50    CONTINUE
         IF(NLINES(IJ).GT.NLIMAX) NLIMAX=NLINES(IJ)
   40 CONTINUE
      WRITE(6,603) NLIMAX
  603 FORMAT(' MAXIMUM NUMBER OF OVERLAPPING TRANSITIONS:  ',I3/)
C
C     recalculate weights for frequency integration
C     after the exclusion of some frequencies
C
      NPPX=0
      DO 100 IJ=1,NFREQ
         IF(IJX(IJ).GT.0) NPPX=NPPX+1
         W(IJ)=0.
         KJ0=KIJ(IJ)
         IF(IJX(JIK(KJ0)).EQ.-1) GO TO 100
         IF(KJ0.GE.2 .AND. KJ0.LT.NFREQ) THEN
            IK1=KJ0-1
            DO WHILE (IJX(JIK(IK1)).EQ.-1)
               IK1=IK1-1
            END DO
            IK2=KJ0+1
            DO WHILE (IJX(JIK(IK2)).EQ.-1)
               IK2=IK2+1
            END DO
            W(IJ)=HALF*ABS(FREQ(JIK(IK1))-FREQ(JIK(IK2)))
          ELSE IF(KJ0.EQ.1) THEN
            W(IJ)=HALF*ABS(FREQ(JIK(KJ0))-FREQ(JIK(KJ0+1)))
          ELSE IF(KJ0.EQ.NFREQ) THEN
            W(IJ)=HALF*ABS(FREQ(JIK(KJ0-1))-FREQ(JIK(KJ0)))
         END IF
  100 CONTINUE
C
C     Correction for Simpson weights
C
      JK1=JIK(1)
      DO 120 IJ=2,NFREQ,2
         JK2=JIK(IJ)
         JK3=JIK(IJ+1)
         IF(IJLIN(JK2).NE.0 .OR. IJLIN(JK3).NE.0) GO TO 130
         IF(WCH(JK2).NE.0.) GO TO 130
         W(JK1)=W(JK1)-SIXTH*W(JK2)
         W(JK3)=W(JK3)-SIXTH*W(JK2)
         W(JK2)=W(JK2)*FTH
         JK1=JK3
  120 CONTINUE
  130 JK1=JIK(NFREQ)
      DO 140 IJ=NFREQ-1,1,-2
         JK2=JIK(IJ)
         JK3=JIK(IJ-1)
         IF(IJLIN(JK2).NE.0 .OR. IJLIN(JK3).NE.0) GO TO 150
         IF(WCH(JK2).NE.0.) GO TO 150
         W(JK1)=W(JK1)-SIXTH*W(JK2)
         W(JK3)=W(JK3)-SIXTH*W(JK2)
         W(JK2)=W(JK2)*FTH
         JK1=JK3
  140 CONTINUE
  150 CONTINUE
C
      DO IJ=1,NFREQ
         W0E(IJ)=W(IJ)*PI4H/FREQ(IJ)
         IF(IJALI(IJ).GT.0) WC(IJ)=W(IJ)
      END DO
C
C     check accuracy of weights for integration
C
  190 Z0=0.
      Z1=0.
      Z2=0.
      ZH=0.
      T1=TEFF
      T2=TWO*TEFF
      T3=HALF*TEFF
      X1=HK/T1
      X2=HK/T2
      X3=HK/T3
      DO 200 IJ=1,NFREQ
         Z0=Z0+W(IJ)
         X15=FREQ(IJ)*1.D-15
         BNZ=BN*X15*X15*X15
         FX1=FREQ(IJ)*X1
         IF(FX1.GT.100.) GO TO 200
         Z1=Z1+W(IJ)*BNZ/(EXP(FREQ(IJ)*X1)-1)
         Z2=Z2+W(IJ)*BNZ/(EXP(FREQ(IJ)*X2)-1)
         ZH=ZH+W(IJ)*BNZ/(EXP(FREQ(IJ)*X3)-1)
  200 CONTINUE
      T1S=SQRT(SQRT(0.25*Z1/SIG4P))
      T1ER=T1S/T1-UN
      T2S=SQRT(SQRT(0.25*Z2/SIG4P))
      T2ER=T2S/T2-UN
      T3S=SQRT(SQRT(0.25*ZH/SIG4P))
      T3ER=T3S/T3-UN
      JK1=JIK(1)
      JK2=JIK(NFREQ)
      Z00=FREQ(JK1)-FREQ(JK2)
      WRITE(6,701) FREQ(JK1),FREQ(JK2),Z00,Z0,T3,T3ER,T1,T1ER,T2,T2ER
  701 FORMAT(/' ACCURACY OF INTEGRATIONS:',/,
     *   ' Interval:',1p4e16.8,/,
     *   15x,' Planck functions:',9x,0pf12.0,4x,1pe12.4,/,
     *   42x,0pf12.0,4x,1pe12.4,/,42x,0pf12.0,4x,1pe12.4,/)
        IF(ISPODF.GT.0) NPPX=NFREQ
      WRITE(6,702) NFREQ,NPPX
  702 FORMAT(' TOTAL NUMBER OF FREQUENCIES:',I8,/,
     *         ' SELECTED FREQUENCIES:       ',I8)
C 
      RETURN
      END
C
C
C ********************************************************************
C
C
 
      SUBROUTINE IROSET
C     =================
C
C     Initialization of opacity sampling for iron-peak lines
C
C     IOBS = 2 : ALL lines except lines to autoionized levels
C     IOBS = 1 : ALL lines
C     IOBS = 0 : only lines between observed levels
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'ATOMIC.FOR'
      INCLUDE 'MODELQ.FOR'
      INCLUDE 'ODFPAR.FOR'
      common/ncfe00/ncfe,nlife,nkulev
c     REAL*4 VDOP,AGAM,SIG0,SIGT
      PARAMETER (CSIG=0.0149736)
      COMMON/LINED/WAVE(MLINE),
c    *             VDOP(MLINE,MDODF),
c    &             AGAM(MLINE,MDODF),SIG0(MLINE,MDODF),
     &             JTR(MLINE,2)
c     DIMENSION SIGT(MDODF,MFREQ) 
      DIMENSION DML(MDEPTH)
C
      JIDR(1)=1
      IF(JIDS.EQ.0) THEN
         JIDR(2)=INT(0.7*ND)
         JIDR(3)=ND
         JIDN=3
         JIDC=2
       ELSE
         I=1
         DO WHILE(JIDR(I).LT.ND)
            I=I+1
            JIDR(I)=JIDR(I-1)+JIDS
            IF(JIDR(I).LE.INT(0.7*ND)) JIDC=JIDR(I)
         END DO
         JIDN=I
         IF(JIDR(I).GT.ND) THEN
            IF(JIDR(I-1).GE.ND-5) JIDN=JIDN-1
            JIDR(JIDN)=ND
         END IF
         IF(JIDN.GT.MDODF)
     &   CALL QUIT(' Too many depths for Fe x-sections',JIDN,MDODF)
      END IF
c
      DO ID=1,ND
         DML(ID)=LOG(DM(ID))
      END DO
      DO I=1,JIDN-1
         DXI=DML(JIDR(I+1))-DML(JIDR(I))
         DO ID=JIDR(I)+1,JIDR(I+1)
            JIDI(ID)=I
            XJID(ID)=(DML(JIDR(I+1))-DML(ID))/DXI
         END DO
      END DO
      JIDI(1)=1
      XJID(1)=1.

      XFRMA=DLOG(FRS1)
      IJD=INT(9./DDNU)
      IF(IJD.LT.2) IJD=2
      NFTT=0
      NFTMX=0
      nkulev=0
      nlife=0
c     write(6,401) jidn
c     write(6,401) (jidr(i),i=1,jidn)
c 401 format(3i5)
C
      DO 500 ION=1,NION
         IND=INODF1(ION)
         IF(IND.LE.0) GO TO 500
         IF(NLLIM(ION).GE.NLEVS(ION)) THEN
            DO ID=1,ND
               DO I=NFIRST(ION),NLAST(ION)
                  WOP(I,ID)=UN
               END DO
            END DO
            GO TO 500
         END IF
C
C       Set up superlevels and read line data
C     
        IOBS=IKOBS(ION)
        CALL LEVCD(ION,IOBS)
        CALL INKUL(ION,IOBS)
C
        if(keve.gt.nkulev) nkulev=keve
        if(kodd.gt.nkulev) nkulev=kodd
        if(nlinku.gt.nlife) nlife=nlinku
C
C       Assign line to supertransition and compute cross-section
C
        N1=NFIRST(ION)
        NLII=NLAST(ION)-N1+1
        DO IL=1,NLII-1
           KEVL=0
           KODL=0
           IF(JEN(IL).LE.NEVKU(ION)) THEN
              KEVL=JEN(IL)
            ELSE
              KODL=JEN(IL)-NEVKU(ION)
           END IF
           ILOK=N1+IL-1
           DO IU=IL+1,NLII
              IUPK=N1+IU-1
              ITR=ITRA(ILOK,IUPK)
              INDXPA=ABS(INDEXP(ITR))
              W1=0.
              W2=0.
              IFRKU=0
              NFT=0
              NLT=0
              KEVU=0
              KODU=0
              IF(JEN(IU).LE.NEVKU(ION)) THEN
                 KEVU=JEN(IU)
               ELSE
                 KODU=JEN(IU)-NEVKU(ION)
              END IF
              IF(KEVL.NE.0) THEN
                 KEV=KEVL
                 KOD=KODU
                 IEO=0
                 GSUPER=YMKU(JEN(IL),1)
               ELSE
                 KEV=KEVU
                 KOD=KODL
                 IEO=1
                 GSUPER=YMKU(JEN(IL)-NEVKU(ION),2)
              END IF
c
       write(10,211) ion,il,kevl,kodl,ilok,iu,kevu,kodu,iupk,kev,kod
  211  format(i4,4i6,2x,4i6,3x,2i7)
c             DO IJ=1,MFREQ
c                DO I=1,MDODF
c                   SIGT(I,IJ)=0.
c                END DO
c             END DO
              FCOL=0.
              DO 10 K=1,NLINKU
                 IF(KSEV(JTR(K,1)).NE.KEV) GO TO 10
                 IF(KSOD(JTR(K,2)).NE.KOD) GO TO 10
                 NLT=NLT+1
                 FRL=CAS/WAVE(K)
                 IJL=NINT((XFRMA-DLOG(FRL))/DXNU)+NFRS1
                 D0=ABS((FREQ(IJL)-FRL)/(FREQ(IJL)-FREQ(IJL+1)))
                 IF(D0.GT.HALF) THEN
                    DO WHILE(FRL.GT.FREQ(IJL))
                       IJL=IJL-1
                    END DO
                    DO WHILE(FRL.LT.FREQ(IJL))
                       IJL=IJL+1
                    END DO
                    D1=FRL-FREQ(IJL)
                    D2=FREQ(IJL-1)-FRL
                    IF(D2.LT.D1) IJL=IJL-1
                 END IF  
                 IJ0=IJL-IJD
                 IJ1=IJL+IJD
                 IF(IJ0.LT.1) IJ0=1
                 IF(IJ1.GT.NFREQ) IJ1=NFREQ
                 IF(IFRKU.EQ.0) IFRKU=IJ0
                 NFT=IJ1-IFRKU+1
c
c                write(10,210) ion,il,iu,k,ijl,ij0,ij1,ifrki,nft
c 210 format(' ion,il,iu,k,ijl,ij0,ij1,ifrku,nft',9i7)
                 DO IJ=IJ0,IJ1
                    DNU=FREQ(IJ)-FRL
c                   DO I=1,JIDN
c                      VV=DNU*dble(VDOP(K,I))
c                      PRFK=VOIGTE(real(VV),AGAM(K,I))/GSUPER
c                      SIGT(I,IJ)=SIGT(I,IJ)+SIG0(K,I)*real(PRFK)
c                   END DO
                 END DO
c                FCOL=FCOL+SIG0(K,JIDC)/VDOP(K,JIDC)
   10         CONTINUE
              OSC0(ITR)=0.
              IF(INDXPA.EQ.3 .OR. INDXPA.EQ.4) OSC0(ITR)=STRLX
              IF(FCOL.GT.0.) OSC0(ITR)=FCOL/GSUPER/CSIG
              IF(NLT.GT.0) THEN
                 W1=CAS/FREQ(IFRKU)
                 W2=CAS/FREQ(IFRKU+NFT-1)
              END IF
              IF(NFT.GT.0) THEN
                 IFR0(ITR)=IFRKU
                 IFR1(ITR)=IFRKU+NFT-1
                 KFR0(ITR)=NFTT+1
                 KFR1(ITR)=NFTT+NFT
                 NFTT=NFTT+NFT
c                if(iprint.ge.4)
c    *           write(6,611) itr,ifr0(itr),ifr1(itr),
c    *                        kfr0(itr),kfr1(itr)
c 611            format(i7,2x,2i9,2x,2i9)
                 DO IJ=IFR0(ITR),IFR1(ITR)
                    KJ=IJ-IFR0(ITR)+KFR0(ITR)
c                   DO I=1,JIDN
c                      sxx=log(dble(sigt(i,ij))+1.d-40)
c                      SIGFE(I,KJ)=real(sxx)
c                   END DO
                 END DO
              END IF
              IF(KJ.GT.MCFE) 
     &           CALL QUIT(' Too many Fe cross-sect. to store',
     &                    KJ,MCFE)
c                WRITE(41,313) IL,IU,W1,W2,IFRKU,NFT,NLT,OSC0(ITR)
C             END IF
              IF(NFT.GT.NFTMX) NFTMX=NFT
           END DO
        END DO
  500 CONTINUE
      ncfe=kj
      WRITE(10,*) ' Max. number of freq. per transition:',NFTMX
      WRITE(10,*) ' Number of superline cross-sections: ',NFTT
C     
      CALL IJALI2   
C     
      RETURN
C     
c 310 FORMAT(F6.2,2I5)
c 311 FORMAT(1P2E17.8,I10,F9.0,I5)
c 312 FORMAT(F9.0,1PE14.5,0P,F9.2)
c 313 FORMAT(2I4,2F12.3,3I10,1PE12.3)
c 314 FORMAT(1P10E12.5)
C     
      END
C
C
C     ***************************************************************
C
C
      SUBROUTINE LEVCD(ION,IOBS)
C     ==========================
C
C     Mean energy and statistical weights of superlevels.
C     Read atomic data from Kurucz CD-ROM file (gf*.gam)
C
C     Setup collisional strengths between superlevels, using
C      Eissner-Seaton formula for each possible individual
C      transition. Assumes Gamma(T)=0.05, and T=Teff
C      Contributions from allowed transitions will be superseded
C      in routine INKUL.
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'ATOMIC.FOR'
      INCLUDE 'MODELQ.FOR'
      INCLUDE 'ODFPAR.FOR'
C
      PARAMETER (BOLCM=1.D8/HK/CAS,CCOR=0.09,SIXTH=UN/6.,GES=0.05)
      DIMENSION GWE(MDEPTH,MLEVEL,2),GWB(MDEPTH,MLEVEL,2),AA(MDEPTH)
      DIMENSION E0FE(10),E0NI(10),E0CR(10)
      COMMON/COLKUR/OMES(100,100),EKU(15000),GKU(15000),GST,
     &              KKU(15000)
      
      DATA E0FE/63480.,130563.,247220.,442000.,605000.,799000.,
     &            1008000.,1218380.,1884000.,2114000./
      DATA E0NI/61590.,146560.,283700.,443000.,613500.,871000.,
     &            1070000.,1310000.,1560000.,1812000./
      DATA E0CR/54576.,132966.,249700.,396500.,560200.,731020.,
     &            1291900.,1490000.,1688000.,1971000./
C
C     Initialization
C
      IF(IOBS.NE.1 .AND. IOBS.NE.2) IOBS=0
      DO I=1,NEVKU(ION)
         YMKU(I,1)=0.
         EMKU(I,1)=0.
         DO ID=1,ND
            GWE(ID,I,1)=0.
            GWB(ID,I,1)=0.
         END DO
      END DO
      DO I=1,NODKU(ION)
         YMKU(I,2)=0.
         EMKU(I,2)=0.
         DO ID=1,ND
            GWE(ID,I,2)=0.
            GWB(ID,I,2)=0.
         END DO
      END DO
      NEVOD=NEVKU(ION)+NODKU(ION)
      IF(NEVOD.GT.100) 
     &   CALL QUIT(' Too many superlevels in a single Fe ion',
     &             NEVKU(ION),NODKU(ION))
      DO I=1,NEVOD
         DO J=1,NEVOD
            OMES(I,J)=0.
         END DO
      END DO
      IWSUP=IFWOP(NFIRST(ION))
      IF(IWSUP.GE.2) THEN
         DO ID=1,ND
            TEMP1(ID)=UN/TEMP(ID)
            AA(ID)=CCOR*EXP(SIXTH*LOG(ELEC(ID)))/SQRT(TEMP(ID))
         END DO
         ZZ=IZ(ION)
         IF(IZ(ION).GT.10) 
     *   CALL QUIT(' Too high Fe, Ni or Cr ion: ion,iz',ion,iz(ion))
         IAT=IATM(NFIRST(ION))
         E0=E0FE(IZ(ION))
         IF(NUMAT(IAT).EQ.28) E0=E0NI(IZ(ION)) 
         IF(NUMAT(IAT).EQ.24) E0=E0CR(IZ(ION)) 
      END IF
C
C     CD-ROM format: Read energy levels 
C
      IUN1=31
      OPEN(IUN1,FILE=FIODF1(ION),STATUS='OLD')
      READ(IUN1,170) NLINKU,KEVE,KODD
      IF(KEVE+KODD.GT.15000)
     &   CALL QUIT(' Too many levels in Kurucz file',keve,kodd)
C
C     Even parity
C
      DO K=1,KEVE
         KSL=0
         READ(IUN1,171) YJ,E,AR,SR,WR
         GEV=(TWO*YJ+UN)
         IF(E.LT.0.) THEN
            E=-E
            IF(IOBS.EQ.0) GO TO 10
         END IF
         IF(E.LE.XEV(1,ION)) KSL=1
         DO I=2,NEVKU(ION)
            IF(E.LE.XEV(I,ION) .AND. E.GT.XEV(I-1,ION)) KSL=I
         END DO
         IF(KSL.EQ.0) WRITE(10,*) ' Error with even levels',E,YJ
         KKU(K)=KSL
         GKU(K)=GEV
         EKU(K)=E
         YMKU(KSL,1)=YMKU(KSL,1)+GEV
         EMKU(KSL,1)=EMKU(KSL,1)+GEV*E
c        IF(IWSUP.EQ.2) THEN
c           EBCM=E/BOLCM
c           DO ID=1,ND
c              GWX=GEV*EXP(-EBCM*TEMP1(ID))
c              GWB(ID,KSL,1)=GWB(ID,KSL,1)+GWX
c              GWE(ID,KSL,1)=GWE(ID,KSL,1)+GWX*E
c           END DO
c         ELSE IF(IWSUP.EQ.3) THEN
c           EBCM=E/BOLCM
c           IF(E.LT.E0) THEN
c              XN=SQRT(E0/(E0-E))
c              DO ID=1,ND
c                 WID=WN(XN,AA(ID),ELEC(ID),ZZ)
c                 GWX=GEV*WID*EXP(-EBCM*TEMP1(ID))
c                 GWB(ID,KSL,1)=GWB(ID,KSL,1)+GWX
c                 GWE(ID,KSL,1)=GWE(ID,KSL,1)+GWX*E
c              END DO
c            ELSE
c              DO ID=1,ND
c                 WID=UN
c                 GWX=GEV*WID*EXP(-EBCM*TEMP1(ID))
c                 GWB(ID,KSL,1)=GWB(ID,KSL,1)+GWX
c                 GWE(ID,KSL,1)=GWE(ID,KSL,1)+GWX*E
c              END DO
c           END IF
c        END IF
   10    EEV(K)=E
         AEV(K)=AR
         SEV(K)=SR
         WEV(K)=WR
         KSEV(K)=KSL
      END DO
      DO I=1,NEVKU(ION)
        IF(YMKU(I,1).EQ.0.)
     *     call quit(' No levels in even superlevel',i,i)
      EMKU(I,1)=EMKU(I,1)/YMKU(I,1)
      END DO
C
C     Odd parity
C
      DO K=1,KODD
         KSL=0
         READ(IUN1,171) YJ,E,AR,SR,WR
         GOD=(TWO*YJ+UN)
         IF(E.LT.0.) THEN
            E=-E
            IF(IOBS.EQ.0) GO TO 20
         END IF
         IF(E.LE.XOD(1,ION)) KSL=1
         DO I=2,NODKU(ION)
            IF(E.LE.XOD(I,ION) .AND. E.GT.XOD(I-1,ION)) KSL=I
         END DO
         IF(KSL.EQ.0) WRITE(10,*) ' Error with odd levels',E,YJ
         KKU(K+KEVE)=KSL+NEVKU(ION)
         GKU(K+KEVE)=GOD
         EKU(K+KEVE)=E
         YMKU(KSL,2)=YMKU(KSL,2)+GOD
         EMKU(KSL,2)=EMKU(KSL,2)+GOD*E
c        IF(IWSUP.EQ.2) THEN
c           EBCM=E/BOLCM
c           DO ID=1,ND
c              GWX=GOD*EXP(-EBCM*TEMP1(ID))
c              GWB(ID,KSL,2)=GWB(ID,KSL,2)+GWX
c              GWE(ID,KSL,2)=GWE(ID,KSL,2)+GWX*E
c           END DO
c         ELSE IF(IWSUP.EQ.3) THEN
c           EBCM=E/BOLCM
c           IF(E.LT.E0) THEN
c              XN=SQRT(E0/(E0-E))
c              DO ID=1,ND
c                 WID=WN(XN,AA(ID),ELEC(ID),ZZ)
c                 GWX=GOD*WID*EXP(-EBCM*TEMP1(ID))
c                 GWB(ID,KSL,2)=GWB(ID,KSL,2)+GWX
c                 GWE(ID,KSL,2)=GWE(ID,KSL,2)+GWX*E
c              END DO
c            ELSE
c              DO ID=1,ND
c                 WID=UN
c                 GWX=GOD*WID*EXP(-EBCM*TEMP1(ID))
c                 GWB(ID,KSL,2)=GWB(ID,KSL,2)+GWX
c                 GWE(ID,KSL,2)=GWE(ID,KSL,2)+GWX*E
c              END DO
c           END IF
c        END IF
   20    EOD(K)=E
         AOD(K)=AR
         SOD(K)=SR
         WOD(K)=WR
         KSOD(K)=KSL
      END DO
      DO I=1,NODKU(ION)
         IF(YMKU(I,2).EQ.0.)
     *   call quit(' No levels in odd superlevel',I,I)
         EMKU(I,2)=EMKU(I,2)/YMKU(I,2)
      END DO
      CLOSE(IUN1)
C
C     Collisional strengths of transitions between super-levels
C      Eissner-Seaton formula (Gamma=0.05)
C
c     GST=8.63E-6*GES/SQRT(TEFF)
c     TK0=UN/BOLCM/TEFF
c     DO I=1,KEVE+KODD-1
c        KI=KKU(I)
c        DO J=I+1,KEVE+KODD
c           KJ=KKU(J)
c           U0=ABS(EKU(I)-EKU(J))*TK0
c           OMES(KI,KJ)=OMES(KI,KJ)+GST*EXP(-U0)
c           OMES(KJ,KI)=OMES(KI,KJ)
c        END DO
c     END DO
C
C     Sort superlevel energies
C
      NLEVKU=NEVKU(ION)+NODKU(ION)
      DO I=1,NEVKU(ION)
         EU(I)=EMKU(I,1)
      END DO
      DO I=1,NODKU(ION)
         EU(I+NEVKU(ION))=EMKU(I,2)
      END DO
      CALL INDEXX(NLEVKU,EU,JEN)
c
C     Superlevel generalized occupation probabilities
C      
c     IF(IWSUP.GE.2) THEN
c        DO I=1,NLEVKU
c           II=NFIRST(ION)+I-1
c           JJ=JEN(I)
c           JK=1
c           IF(JJ.GT.NEVKU(ION)) THEN
c              JJ=JEN(I)-NEVKU(ION)
c              JK=2
c           END IF
c           DO ID=1,ND
c              ESUP=GWE(ID,JJ,JK)/GWB(ID,JJ,JK)
c              WSUP=EXP(ESUP/BOLCM*TEMP1(ID))/YMKU(JJ,JK)
c              WOP(II,ID)=WSUP*GWB(ID,JJ,JK)
c           END DO
c        END DO
c     END IF
C
  170 FORMAT(I7,13X,I6,12X,I6)
  171 FORMAT(8X,F4.1,4X,F13.3,18X,3E9.2)
C
      RETURN
      END
C
C
C     ***************************************************************
C
C
      SUBROUTINE INKUL(ION,IOBS)
C     ==========================
C
C     Read line list from Kurucz CD-ROM files  (gf*.lin)
C
C     INPUT:  Unit 18
C               WMIN : Min. wavelength (lines at smaller wave are NOT selected)
C               WMAX : Max. wavelength (lines at larger wave are NOT selected)
C               IOBS : Type of selected lines
C
C     OUTPUT:  fill up common/lined/
C     *******
C           -  WAVE     : wavelength in ANGSTROMS
C           -  SIG0     : 0.02654/sqrt(pi)*gf/VDOP (divided by g(super) later)
C           -  VDOP     : Doppler width
C           -  AGAM     : Damping parameter
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'ATOMIC.FOR'
      INCLUDE 'MODELQ.FOR'
      INCLUDE 'ODFPAR.FOR'
c     REAL*4 VDOP,AGAM,SIG0
C
      PARAMETER (TEN=10.,TENLG=2.302585093,GES=0.05)
      PARAMETER (BOL2=2.76108D-16,BOLCM=1.D8/HK/CAS)
      PARAMETER (CSTK=3.54,PSTK=2./3.,TSTK=UN/6.)
      PARAMETER (CVDW=3.74,PVDW=0.4,TVDW=0.3)
      PARAMETER (PI4V=0.25/3.141592654,CSIG=0.0149736)
      PARAMETER (EXPIA1=-0.57721566,EXPIA2=0.99999193,
     *           EXPIA3=-0.24991055,EXPIA4=0.05519968,
     *           EXPIA5=-0.00976004,EXPIA6=0.00107857,
     *           EXPIB1=0.2677734343,EXPIB2=8.6347608925,
     *           EXPIB3=18.059016973,EXPIB4=8.5733287401,
     *           EXPIC1=3.9584969228,EXPIC2=21.0996530827,
     *           EXPIC3=25.6329561486,EXPIC4=9.5733223454)
      DIMENSION E0FE(10),E0NI(10),E0CR(10)
      COMMON/COLKUR/OMES(100,100),EKU(15000),GKU(15000),GST,
     &              KKU(15000)
      COMMON/LINED/WAVE(MLINE),
c    *             VDOP(MLINE,MDODF),
c    *             AGAM(MLINE,MDODF),SIG0(MLINE,MDODF),
     &             JTR(MLINE,2)
      
      DATA E0FE/63480.,130563.,247220.,442000.,605000.,799000.,
     &            1008000.,1218380.,1884000.,2114000./
      DATA E0NI/61590.,146560.,283700.,443000.,613500.,871000.,
     &            1070000.,1310000.,1560000.,1812000./
      DATA E0CR/54576.,132966.,249700.,396500.,560200.,731020.,
     &            1291900.,1490000.,1688000.,1971000./
C
      IAT=IATM(NFIRST(ION))
      E0=E0FE(IZ(ION))
      IF(NUMAT(IAT).EQ.28) E0=E0NI(IZ(ION)) 
      IF(NUMAT(IAT).EQ.24) E0=E0CR(IZ(ION)) 
      CDOP=BOL2/AMASS(IAT)
      TK35=UN/BOLCM/TEMP(JIDR(2))
      CVR=19.7363/TEMP(JIDR(2))/SQRT(TEMP(JIDR(2)))
      if(jidn.gt.3) then
         TK35=UN/BOLCM/TEFF
         CVR=19.7363/TEFF/SQRT(TEFF)
      end if
      
      XION=0.
      DO K=1,JIDN
         XIONI=0.
         XIATI=0.
         KJJ=JIDR(K)
         DO I=NFIRST(ION),NLAST(ION)
            XIONI=XIONI+POPUL(I,KJJ)
         END DO
         DO I=N0A(IAT),NKA(IAT)
             XIATI=XIATI+POPUL(I,KJJ)
         END DO
         XIONK=XIONI/XIATI
         IF(XIONK.GT.XION) XION=XIONK
c        write(6,800) k,ion,nfirst(ion),popul(nfirst(ion),k),
c    *   xioni,xiati
c 800 format(3i6,1p3e10.2)
      END DO
      if(xion.le.1.) xion=1.
c     write(6,801) ion,xion
c 801 format(' ion, xion ',i4,1pe10.2)
  
      NLINKU=0
      WMIN=CAS/FRS1/TEN
      WMAX=CAS/FRS2/TEN
      IUN2=32
      OPEN(IUN2,FILE=FIODF2(ION),STATUS='OLD')
c     iii=0
   10 READ(IUN2,180,ERR=20,END=20) WA,GFR,JEVR,JODR,IFPLI
      GF=EXP(TENLG*GFR)
      IF(WA.GT.WMAX) GO TO 11
      IF(WA.LT.WMIN) GO TO 11
      IF(IOBS.EQ.0 .AND. IFPLI.EQ.1) GO TO 10
      IF(IOBS.EQ.2 .AND. EOD(JODR).GT.E0) GO TO 10
      IF(IOBS.EQ.2 .AND. EEV(JEVR).GT.E0) GO TO 10
      E00=EEV(JEVR)
      IF(EOD(JODR).LT.EEV(JEVR)) E00=EOD(JODR)
      XLSTR=XION*GF*EXP(-E00*TK35)
      if(jidn.gt.3) XLSTR=XION*GF*EXP(-E00*8./E0)
c      iii=iii+1
c      if(iii.le.20) write(6,802) gf,e00,xlstr,strlx
c 802  format(1p4e10.2)
      IF(XLSTR.LT.STRLX) GO TO 10
      NLINKU=NLINKU+1
      WAVE(NLINKU)=WA*TEN
      JTR(NLINKU,1)=JEVR
      JTR(NLINKU,2)=JODR
c     GR=AEV(JEVR)+AOD(JODR)
c     C4=SEV(JEVR)
c     C4P=SOD(JODR)
c     SMX=DMAX1(ABS(C4P-C4),DMIN1(ABS(C4),ABS(C4P)))
c     GSLOG0=CSTK+PSTK*DLOG10(SMX)
c     DO I=1,JIDN
c        ID=JIDR(I)
c        VT2=VTURBS(ID)*VTURBS(ID)
c        VDOP(NLINKU,I)=real(WAVE(NLINKU)/1D+8/SQRT(CDOP*TEMP(ID)+VT2))
c        GS=EXP(TENLG*(GSLOG0+TSTK*DLOG10(TEMP(ID))))
c        AGAM(NLINKU,I)=real((GR+GS*ELEC(ID))*PI4V*VDOP(NLINKU,I))
c        SIG0(NLINKU,I)=real(CSIG*GF*VDOP(NLINKU,I))
c     END DO
 
 11   KA=KKU(JEVR)
c     KB=KKU(JODR+KEVE)
c     IF(KA.LE.KB) THEN
c        K1=KA
c        K2=KB
c      ELSE
c        K1=KB
c        K2=KA
c     END IF
c     IF(K1.EQ.K2) GO TO 10
c     U0=1.E7/WA*TK35
c     IF(U0.LE.UN) THEN
c        EXPIU0=-LOG(U0)+EXPIA1+U0*(EXPIA2+U0*(EXPIA3+U0*(EXPIA4+
c    *           U0*(EXPIA5+U0*EXPIA6))))
c      ELSE
c        EXPIU0=EXP(-U0)*((EXPIB1+U0*(EXPIB2+U0*(EXPIB3+
c    *          U0*(EXPIB4+U0))))/(EXPIC1+U0*(EXPIC2+
c    *          U0*(EXPIC3+U0*(EXPIC4+U0)))))/U0
c     END IF
c     GB=0.276*EXP(U0)*EXPIU0
c     IF(GB.LT.0.25) GB=0.25
c     OMES(K1,K2)=OMES(K1,K2)+(CVR/U0*GF*GB-GST)*EXP(-U0)

      GO TO 10
   20 CLOSE(IUN2)
      WRITE(10,600) NUMAT(IAT),IZ(ION),NLINKU
C
C     Store collisional excitation strengths
C
      NLEVKU=NEVKU(ION)+NODKU(ION)
      write(10,610) ion,nlevku
  610 format('  ion,nlevku',2i7)
c     DO I=1,NLEVKU-1
c        II=NFIRST(ION)+I-1
c        I1=JEN(I)
c        IF(I1.LE.NEVKU(ION)) THEN
c           GSUP=YMKU(I1,1)
c         ELSE
c           GSUP=YMKU(I1-NEVKU(ION),2)
c        END IF
c        DO J=I+1,NLEVKU
c           JJ=NFIRST(ION)+J-1
c           J1=JEN(J)
c           IT=ITRA(II,JJ)
c           C2=CPAR(IT)
c           OMECOL(II,JJ)=OMES(I1,J1)/GSUP*C2/GES
c           OMECOL(JJ,II)=OMECOL(II,JJ)
c        END DO
c     END DO
c      
      RETURN
C       
  180 FORMAT(F11.4,F7.3,2I4,I1)
  600 FORMAT(' Ion',2I3,' : ',I9,' Lines included')
      END
C
C
C     ***************************************************************
C
C
C

      SUBROUTINE OPACFL(IJ)
C     =====================
C
C     Absorption, emission, and scattering coefficients
C     at frequency IJ and for all depths
C
C     Input: IJ   opacity and emissivity is calculated for the
C                 frequency points with index IJ
C     Output: ABSO1 -  array of absorption coefficient
C             EMIS1 -  array of emission coefficient
C             SCAT1 -  array of scattering coefficient (all scattering
C                       mechanisms except electron scattering)
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'ATOMIC.FOR'
      INCLUDE 'MODELQ.FOR'
      INCLUDE 'ODFPAR.FOR'
      INCLUDE 'ALIPAR.FOR'
      PARAMETER (C14=2.99793D14, CFF1=1.3727D-25)
C
C     initialize
c
      DO ID=1,ND
         ABSO1(ID)=ELSCAT(ID)
         EMIS1(ID)=0.
         SCAT1(ID)=ELSCAT(ID)
         ABSO1L(ID)=0.
         EMIS1L(ID)=0.
      END DO
C
C     basic frequency- and depth-dependent quantities
C
      FR=FREQ(IJ)
      FRINV=UN/FR
      FR3INV=FRINV*FRINV*FRINV
      DO 10 ID=1,ND
         XKF(ID)=EXP(-HKT1(ID)*FR)
         XKF1(ID)=UN-XKF(ID)
         XKFB(ID)=XKF(ID)*BNUE(IJ)
   10 CONTINUE
C
C ********  1a. bound-free contribution - without dielectronic rec.
C
      if(ifdiel.eq.0) then
      DO IBFT=1,NTRANC
         ITR=ITRBF(IBFT)
         SG=CROSS(IBFT,IJ)
         IF(SG.GT.0.) THEN
         II=ILOW(ITR)
         JJ=IUP(ITR)
         IZZ=IZ(IEL(II))
         IMER=IMRG(II)
         DO ID=1,ND
            SGD=SG
            IF(MCDW(ITR).GT.0) THEN
               CALL DWNFR1(FR,FR0(ITR),ID,IZZ,DW1)
               DWF1(MCDW(ITR),ID)=DW1
               SGD=SG*DW1
            END IF
            IF(IFWOP(II).LT.0) THEN
               CALL SGMER1(FRINV,FR3INV,IMER,ID,SGME1)
               SGMG(IMER,ID)=SGME1
               SGD=SGME1
            END IF
            EMISBF=SGD*EMTRA(ITR,ID)
            ABSO1(ID)=ABSO1(ID)+SGD*ABTRA(ITR,ID)
            EMIS1(ID)=EMIS1(ID)+EMISBF
            if(iluctr(itr).gt.0) then
               ABSO1L(ID)=ABSO1L(ID)+SGD*ABTRA(ITR,ID)
               EMIS1L(ID)=EMIS1L(ID)+EMISBF
            end if
         END DO
         END IF
      END DO
      else
C
C ********  1. bound-free contribution - with dielectronic rec.
C
      DO IBFT=1,NTRANC
         ITR=ITRBF(IBFT)
         SG=CROSS(IBFT,IJ)
         IF(SG.GT.0.) THEN
         II=ILOW(ITR)
         JJ=IUP(ITR)
         IZZ=IZ(IEL(II))
         IMER=IMRG(II)
         DO ID=1,ND
            SG=CROSSD(IBFT,IJ,ID)
            IF(SG.GT.0.) THEN
            SGD=SG
            IF(MCDW(ITR).GT.0) THEN
               CALL DWNFR1(FR,FR0(ITR),ID,IZZ,DW1)
               DWF1(MCDW(ITR),ID)=DW1
               SGD=SG*DW1
            END IF
            IF(IFWOP(II).LT.0) THEN
               CALL SGMER1(FRINV,FR3INV,IMER,ID,SGME1)
               SGMG(IMER,ID)=SGME1
               SGD=SGME1
            END IF
            EMISBF=SGD*EMTRA(ITR,ID)
            ABSO1(ID)=ABSO1(ID)+SGD*ABTRA(ITR,ID)
            EMIS1(ID)=EMIS1(ID)+EMISBF
            if(iluctr(itr).eq.0) then
              ABSO1L(ID)=ABSO1L(ID)+SGD*ABTRA(ITR,ID)
              EMIS1L(ID)=EMIS1L(ID)+EMISBF
            endif
            END IF
         END DO
         END IF
      END DO
      end if
C
C ******** 2. free-free contribution
C
      DO 40 ION=1,NION
         IT=ITRA(NNEXT(ION),NNEXT(ION))
C
C        hydrogenic with Gaunt factor = 1
C
         IF(IT.EQ.1) THEN
            DO ID=1,ND
               SF1=SFF3(ION,ID)*FR3INV
               SF2=SFF2(ION,ID)
               IF(FR.LT.FF(ION)) SF2=UN/XKF(ID)
               ABSOFF=SF1*SF2
               ABSO1(ID)=ABSO1(ID)+ABSOFF
               EMIS1(ID)=EMIS1(ID)+ABSOFF
            END DO
C
C         hydrogenic with exact Gaunt factor 
C
          ELSE IF(IT.EQ.2) THEN
            DO ID=1,ND
               SF1=SFF3(ION,ID)*FR3INV
               SF2=SFF2(ION,ID)
               IF(FR.LT.FF(ION)) SF2=UN/XKF(ID)
               X=C14*CHARG2(ION)/FR
               SF2=SF2-UN+GFREE1(ID,X)
               ABSOFF=SF1*SF2
               ABSO1(ID)=ABSO1(ID)+ABSOFF
               EMIS1(ID)=EMIS1(ID)+ABSOFF
            END DO
C
C         H minus free-free opacity
C
          ELSE IF(IT.EQ.3) THEN
            DO ID=1,ND
               ABSOFF=(CFF1+CFFT(ID)*FRINV)*CFFN(ID)*FRINV
               ABSO1(ID)=ABSO1(ID)+ABSOFF
               EMIS1(ID)=EMIS1(ID)+ABSOFF
            END DO
C
C         special evaluation of the cross-section
C
          ELSE IF(IT.LT.0) THEN
            DO ID=1,ND
               ABSOFF=FFCROS(ION,IT,TEMP(ID),FR)*
     *                POPUL(NNEXT(ION),ID)*ELEC(ID)
               ABSO1(ID)=ABSO1(ID)+ABSOFF
               EMIS1(ID)=EMIS1(ID)+ABSOFF
            END DO
         END IF
   40 CONTINUE
C
C     ********  3. - additional continuum opacity (OPADD)
C
      IF(IOPADD.NE.0) THEN
         ICALL=1
         DO ID=1,ND
            CALL OPADD(0,ICALL,IJ,ID)
            ABSO1(ID)=ABSO1(ID)+ABAD
            EMIS1(ID)=EMIS1(ID)+EMAD
            SCAT1(ID)=SCAT1(ID)+SCAD
         END DO
      END IF
C
C ********  4. - opacity and emissivity in lines
C
      IF(ISPODF.EQ.0) THEN
      IF(IJLIN(IJ).GT.0) THEN
C
C     the "primary" line at the given frequency
C
         ITR=IJLIN(IJ)
         DO ID=1,ND
            SG=PRFLIN(ID,IJ)
            ABSO1(ID)=ABSO1(ID)+SG*ABTRA(ITR,ID)
            EMIS1(ID)=EMIS1(ID)+SG*EMTRA(ITR,ID)
         END DO
         if(iluctr(itr).gt.0) then
            DO ID=1,ND
               SG=PRFLIN(ID,IJ)
               ABSO1L(ID)=ABSO1L(ID)+SG*ABTRA(ITR,ID)
               EMIS1L(ID)=EMIS1L(ID)+SG*EMTRA(ITR,ID)
            END DO
         end if
      END IF
      IF(NLINES(IJ).LE.0) GO TO 200
C
C     the "overlapping" lines at the given frequency
C
      DO 100 ILINT=1,NLINES(IJ)
         ITR=ITRLIN(ILINT,IJ)
         if(linexp(itr)) go to 100
         IJ0=IFR0(ITR)
         DO IJT=IJ0,IFR1(ITR)
            IF(FREQ(IJT).LE.FR) THEN
               IJ0=IJT
               GO TO 70
            END IF
         END DO
   70    IJ1=IJ0-1
         A1=(FR-FREQ(IJ0))/(FREQ(IJ1)-FREQ(IJ0))
         A2=UN-A1
         DO ID=1,ND
            SG=A1*PRFLIN(ID,IJ1)+A2*PRFLIN(ID,IJ0)
            ABSO1(ID)=ABSO1(ID)+SG*ABTRA(ITR,ID)
            EMIS1(ID)=EMIS1(ID)+SG*EMTRA(ITR,ID)
         END DO
  100 CONTINUE
  200 CONTINUE
C
C     Opacity sampling option
C
      ELSE
      IF(NLINES(IJ).LE.0) GO TO 400
      DO 300 ILINT=1,NLINES(IJ)
         ITR=ITRLIN(ILINT,IJ)
c        IF(LINEXP(ITR)) GO TO 300
C        IF(IJ.LT.IFR0(ITR) .OR. IJ.GT.IFR1(ITR)) GO TO 300
         KJ=IJ-IFR0(ITR)+KFR0(ITR)
         INDXPA=IABS(INDEXP(ITR))
         IF(INDXPA.NE.3 .AND. INDXPA.NE.4) THEN
         DO ID=1,ND
            SG=PRFLIN(ID,KJ)
            ABSO1(ID)=ABSO1(ID)+SG*ABTRA(ITR,ID)
            EMIS1(ID)=EMIS1(ID)+SG*EMTRA(ITR,ID)
         END DO
      ELSE
        DO ID=1,ND
           KJD=JIDI(ID)
c          SG=EXP(XJID(ID)*SIGFE(KJD,KJ)+
c    *        (UN-XJID(ID))*SIGFE(KJD+1,KJ))
           ABSO1(ID)=ABSO1(ID)+SG*ABTRA(ITR,ID)
           EMIS1(ID)=EMIS1(ID)+SG*EMTRA(ITR,ID)
        END DO
      END IF
  300 CONTINUE
  400 CONTINUE
      END IF
C
C     ----------------------------
C     total opacity and emissivity
C     ----------------------------
C
      DO ID=1,ND
         ABSO1(ID)=ABSO1(ID)-EMIS1(ID)*XKF(ID)
         EMIS1(ID)=EMIS1(ID)*XKFB(ID)
         ABSO1L(ID)=ABSO1L(ID)-EMIS1L(ID)*XKF(ID)
         EMIS1L(ID)=EMIS1L(ID)*XKFB(ID)
         ABSO1L(ID)=ABSO1(ID)-ABSO1L(ID)
         EMIS1L(ID)=EMIS1(ID)-EMIS1L(ID)
      END DO
      RETURN
      END
C
C
C     ***************************************************************
C
C
      subroutine rdatax(itr,ic,iunit)
c     ===============================
c
c     for itr, itrx ne 0 - read input data for an individual transition
c     and prepare necessary auxiliary arrays
c
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'ATOMIC.FOR'
      INCLUDE 'MODELQ.FOR'
c
      parameter (mtrx=1000)
      dimension iex(mtrx),itrind(mtrx),izx0(mtrx),izx1(mtrx),
     *          nmaxx(mtrx),izx(mtrx),nshx(mtrx),nax(mtrx),icx(mtrx),
     *          etx(mtrx),ssx(mtrx),dx(mtrx),
     *          aphx(11,5,mtrx),bphx(5,mtrx),a(11,5),b(5)
c
      if(itr.gt.0) then
         itrx=itrx+1
         ntrx=itrx
         ii=ilow(itr)
         ie=iel(ii)
         jj=iup(itr)
         iex(itrx)=ie
         icx(itrx)=ic
         itrind(itrx)=itr
         izx0(itrx)=iz(ie)
         izx1(itrx)=jj-1000
c
c        read inner-shell photoionization data from Omer's tables
c
         read(iunit,*) etx(itrx)
         read(iunit,*) nmaxx(itrx),izx(itrx),nshx(itrx)
         read(iunit,*) ssx(itrx)
         read(iunit,*) nax(itrx)
         read(iunit,*) dx(itrx)
         do i=1,nax(itrx)
            read(iunit,*) bphx(i,itrx)
         end do
         do j=1,nax(itrx)
            do i=1,11
               read(iunit,*) aphx(i,j,itrx)
            end do
         end do
       else if(itr.eq.0) then
c
c      for itr=0 - set up array BFCS with actual cross-sections
c
         if(ntrx.gt.0) then
            do itx=1,ntrx
               ie=iex(itx)
               it=itrind(itx)
               ii=ilow(it)
               ia=iatm(ii)
               iz1=izx1(itx)
               ic=icx(itx)
               jj=0
               do i=1,nlevel
                  if(iatm(i).eq.ia) then
                     if(iz(iel(i)).eq.iz1) then
                        jj=i
                        go to 10
                     end if
                  end if
               end do
   10          continue
               if(jj.eq.0) then
                  if(iz1.eq.iz(iel(nka(ia)-1))+1) jj=nka(ia)
               end if
               if(jj.eq.0) indexp(it)=0
c
               if(indexp(it).ne.0) then
                  iup(it)=jj
                  fr0(it)=etx(itx)/4.1357e-15
                  line(it)=.false.
                  itrcon(it)=ic
                  if(icol(it).ne.99) then
                     itra(ii,jj)=it
                     itra(jj,ii)=ic
                  end if
c  
                  write(6,601) itx,it,ii,jj,ia,
     *                         itra(ii,jj),itra(jj,ii),
     *                         icol(it),fr0(it),etx(itx)
  601             format(8i4,1pe12.4,0pf10.2)
               end if
            end do
         end if
c
       else
         nfreqb=nfreq
         if(ibfint.gt.0) nfreqb=nfreqc
         if(ntrx.gt.0) then
            do itx=1,ntrx
               it=itrind(itx)
               ic=icx(itx)
               if(indexp(it).ne.0) then
                  na=nax(itx)
                  do i=1,na
                     b(i)=bphx(i,itx)
                  end do
                  do j=1,na
                     do i=1,11
                        a(i,j)=aphx(i,j,itx)
                     end do
                  end do
c         
                  do ij=1,nfreqb
                     call bkhsgo(freq(ij),etx(itx),dx(itx),b,na,a,
     *               ssx(itx),nmaxx(itx),izx(itx),nshx(itx),sg)
                     bfcs(ic,ij)=real(sg)
                  end do
                  write(97,681) it,ic,ilow(it),iup(it),bfcs(ic,1)
  681             format(4i5,1p1e15.5)
               end if
            end do
         end if
      end if
      return
      end
C     
C     
C ************************************************************************    
C     
C     
      subroutine bkhsgo(freq,et,d,b,na,a,ss,nmax,iz,nsh,sg)
c     =====================================================
c
c
c subroutine to calculate K and L shell photoionization cross-sections
c       -based on Tim Kallman's bkhsgo subroutine from XSTAR, modified
c        by Omer Blaes 5-7-98.
c       -na.ne.2 bug corrected on 2/24/00 by Omer Blaes
c
c freq is photon frequency in Hz (note that this subroutine immediately
c converts it into eV)
c
c et is threshold energy in eV
c
c iz is the ionization stage of the species being photoionized (1=neutral etc.)
c
c ss is the iz'th element of the array ss(nmax) in Tim's original version
c
c sg is returned as the contribution to the photoionization cross-section
c in cm^2 due to whatever process is being considered.
c
c     this routine does the work in computing cross sections by the
c     method of barfield, et. al.
c
c
c
      INCLUDE 'IMPLIC.FOR'
      dimension b(5),a(11,5)
c
c      data sigth/6.65e-25/
      data sigth/1.e-34/
c
c
      tmp1 = 0.
      jj = 1
      epii = 4.1357e-15*freq
      sg=0.
      if ( epii.gt.et ) then
        xx = epii*(1.e-3) - d
        if ( xx.le.0. ) goto 100
        do nna=1,na
          if ( xx.ge.b(nna) ) jj = jj + 1
        enddo
        if ( jj.le.na ) then
        if(xx.lt.0.) xx=0.
          yy = log10(xx)
          tmp = 0.
          do 5 kkk = 1,11
            kk = 12 - kkk
            tmp = a(kk,jj) + yy*tmp
 5        continue
          if(tmp.lt.-50.) tmp=-50.
          if(tmp.gt.24.) tmp=24.
          sgtmp = 10.**(tmp-24.)
        nelec=nmax+1-iz
        if(nelec.gt.nsh) nelec=nsh
          enelec = float(nelec)
          tmp1o = tmp1
          tmp1=sgtmp*ss
        if(tmp1.lt.sigth*enelec) tmp1=sigth*enelec
          if ( epii.ge.5.e+4 ) then
           if(tmp1.gt.tmp1o) tmp1 = tmp1o
          end if
          sg = sg + tmp1
        endif
      endif
 100  continue
      return
      end
C     
C     
C ************************************************************************    
C     
C     
c
      function cion(n,j,e,t)
c
c     collisional ionization rate from Raymond
c     rate is returned in units cm^3/s
c     inputs are: n=nuclear charge,
c                 j=ion stage
c                 e=valence shell ionization threshold (in eV)
c                 t=temperature in K.
c
c Tim Kallman's routine for calculating collisional ionization rates.
c    Note that this routine only accounts for valence shell ionization.
c    It should be called only once for each ion stage and the valence
c    shell ionization threshold should be the lowest one (e.g. the 2p
c    ionization potential for CI).
c   
c
      INCLUDE 'IMPLIC.FOR'
c
c  sm younger jqsrt 26, 329; 27, 541; 29, 61   with moores for undone
c  a0 for b-like ion has twice 2s plus one 2p  as in summers et al
c  chi = kt / i
c
      dimension a0(30),a1(30),a2(30),a3(30),b0(30),b1(30),
     &          b2(30),b3(30),c0(30),c1(30),c2(30),c3(30),
     &          d0(30),d1(30),d2(30),d3(30)
c
      data a0/13.5,27.0,9.07,11.80,20.2,28.6,37.0,45.4,
     &     53.8,62.2,11.7,38.8,37.27,46.7,57.4,67.0,
     &     77.8,90.1,106.,120.8,135.6,150.4,165.2,180.0,
     &     194.8,209.6,224.4,239.2,154.0,268.8/
      data a1/ - 14.2,-60.1,4.30,27*0./
      data a2/40.6,140.,7.69,27*0./
      data a3/ - 17.1,-89.8,-7.53,27*0./
c
      data b0/ - 4.81,-9.62,-2.47,-3.28,-5.96,-8.64,-11.32,
     &     -14.00,-16.68,-19.36,-4.29,-16.7,-14.58,-16.95,
     &     -19.93,-23.05,-26.00,-29.45,-34.25,-38.92,
     &     -43.59,-48.26,-52.93,-57.60,-62.27,-66.94,
     &     -71.62,-76.29,-80.96,-85.63/
      data b1/9.77,33.1,-3.78,27*0./
      data b2/ - 28.3,-82.5,-3.59,27*0./
      data b3/11.4,54.6,3.34,27*0./
c
      data c0/1.85,3.69,1.34,1.64,2.31,2.984,3.656,4.328,
     &     5.00,5.672,1.061,1.87,3.26,5.07,6.67,8.10,
     &     9.92,11.79,7.953,8.408,8.863,9.318,9.773,
     &     10.228,10.683,11.138,11.593,12.048,12.505,12.96/
      data c1/0.,4.32,.343,27*0./
      data c2/0.,-2.527,-2.46,27*0./
      data c3/0.,.262,1.38,27*0./
c
      data d0/ - 10.9,-21.7,-5.37,-7.58,-12.66,-17.74,
     &     -22.82,-27.9,-32.98,-38.06,-7.34,-28.8,-24.87,
     &     -30.5,-37.9,-45.3,-53.8,-64.6,-54.54,-61.70,
     &     -68.86,-76.02,-83.18,-90.34,-97.50,-104.66,
     &     -111.82,-118.98,-126.14,-133.32/
      data d1/8.90,42.5,-12.4,27*0./
      data d2/ - 35.7,-131.,-8.09,27*0./
      data d3/16.5,87.4,1.23,27*0./
c
      cion = 0.
      chir = t/(11590.*e)
      if ( chir.le..0115 ) return
      chi=chir
      if(chi.lt.0.1) ch=0.1
      ch2 = chi*chi
      ch3 = ch2*chi
      alpha = (.001193+.9764*chi+.6604*ch2+.02590*ch3)
     &        /(1.0+1.488*chi+.2972*ch2+.004925*ch3)
      beta = (-.0005725+.01345*chi+.8691*ch2+.03404*ch3)
     &       /(1.0+2.197*chi+.2457*ch2+.002503*ch3)
      j2 = j*j
      j3 = j2*j
      iso = n - j + 1
c
      a = a0(iso) + a1(iso)/j + a2(iso)/j2 + a3(iso)/j3
      b = b0(iso) + b1(iso)/j + b2(iso)/j2 + b3(iso)/j3
      c = c0(iso) + c1(iso)/j + c2(iso)/j2 + c3(iso)/j3
      d = d0(iso) + d1(iso)/j + d2(iso)/j2 + d3(iso)/j3
c
c  fe ii experimental ionization montague et al: d. neufeld fit
      if ( n.eq.26 .and. j.eq.2 ) then
         a = -13.825
         b = -11.0395
         c = 21.07262
         d = 0.
      endif
c
      ch = 1./chi
      fchi = 0.3*ch*(a+b*(1.+ch)+(c-(a+b*(2.+ch))*ch)*alpha+d*beta*ch)
c      x=-1./chir
c      if(x.lt.-30) x=-30.
c      if(x.gt.30) x=30.
c      expo=exp(x)
c     cion = 2.2e-6*sqrt(chir)*fchi*expo/(e*sqrt(e))
      cion = 2.2e-6*sqrt(chir)*fchi*expo(-1./chir)/(e*sqrt(e))
c     cion = 2.2e-6*sqrt(chir)*fchi*exp(-1./chir)/(e*sqrt(e))
      return
      end
C     
C     
C ************************************************************************    
C     
C     
c
c
      subroutine dielrc(iatom,iont,temp,xpx,dirt,sig0)
c     ================================================
c
      INCLUDE 'IMPLIC.FOR'      
c
c
c Modification of Tim Kallman's XSTAR routine rrrec to calculate dielectronic
c recombination rates (only) to individual ionic species (modified by Omer
c Blaes 5-8-98)
c
c Here temp=temperature in K, xpx is the number density of atomic nuclei in
c cm^{-3} (hydrogen density=xpx/1.1, helium density=xpx*0.1/1.1),
c dirt  =the dielectronic rate in cm^3/s, 
c sig0  = the value of sigma_0, the corresponding pseudo-cross-section
c iatom - atomic number (1=H, 2=He, 6=C, etc.)
c iont  - ionization stage (1 for neutrals, 2 for once ionized, etc.)
c
c     this routine computes radiative recombination rates, both rr 
c     and dr.  rates are output in units of cm**3/s for each ion 
c     stage, where ions are numbered from 1-168:
c       1=HI, 2=HeI, 3=HeII, 4=C I, 5=C II, ..., 9=C VI, 10=N I, ...,
c       16=N VII, 17=O I, ..., 24=O VIII, 25=Ne I, ..., 34=Ne X, 
c       35=Mg I, ..., 46=Mg XII, 47=Si I, ..., 60=Si XIV, 61=S I, ...,
c       76=S XVI, 77=Ar I, ..., 94=Ar XVIII, 95=Ca I, ... 114=Ca XX,
c       115=Fe I, ..., 140=Fe XXVI, 141=Ni I, ... 168=Ni XXVIII.
c     inputs are rate coefficients from Aldrovandi and Pequignot, Storey,
c     and from Arnaud and Raymond for iron
c
      parameter (nni=168)
      parameter (cons=0.1239529*3.28805e15/13.595)
      dimension inid(28,28),uu(28,28)
      dimension adi(nni),bdi(nni),t0(nni),t1(nni),cdd(nni)
      dimension dcfe(26,4),defe(26,4)
      dimension gli(20),gfe(26),gni(28)
      dimension istorey(13),rstorey(5,13)
c
c
c Each non-indented line in the following data statements corresponds
c to each of the elements H, He, C, N, O, Ne, Mg, Si, S, Ar, Ca, Fe,
c and Ni.
c
      data inid/1, 27*0,
     *          2, 3, 26*0,
     *         28*0,
     *         28*0,
     *         28*0,
     *          4, 5, 6, 7, 8, 9, 22*0,
     *         10,11,12,13,14,15,16, 21*0,
     *         17,18,19,20,21,22,23,24, 20*0,
     *         28*0,
     *         25,26,27,28,29,30,31,32,33,34, 18*0,
     *         28*0,
     *         35,36,37,38,39,40,41,42,43,44,45,46, 16*0,
     *         28*0,
     *         47,48,49,50,51,52,53,54,55,56,57,58,59,60,14*0,
     *         28*0,
     *         61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,12*0,
     *         28*0,
     *         77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,
     *         10*0,
     *         28*0,
     *         95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,
     *         110,111,112,113,114,8*0,
     *         28*0,
     *         28*0,
     *         28*0,
     *         28*0,
     *         28*0,
     *         115,116,117,118,119,120,121,122,123,124,125,126,127,128,
     *         129,130,131,132,133,134,135,136,137,138,139,140,2*0,     
     *         28*0,
     *         141,142,143,144,145,146,147,148,149,150,151,152,153,154,
     *         155,156,157,158,159,160,161,162,163,164,165,166,167,168/
      data adi/0.,
     %         1.9E-03,0.,
     %         6.9E-04,7.0E-03,3.8E-03,4.8E-02,4.8E-02,0.,
     %         5.2E-04,1.7E-03,1.2E-02,5.5E-03,7.6E-02,6.6E-02,0.,
     %         1.4E-03,1.4E-03,2.8E-03,1.7E-02,7.1E-03,1.1E-01,8.6E-02,
     %          0.0,
     %         1.3E-03,3.1E-03,7.5E-03,5.7E-03,1.0E-02,4.0E-02,1.1E-02,
     %          1.8E-01,1.3E-01,0.,
     %         1.7E-3,3.5E-3,3.9E-3,9.3E-3,1.5E-2,1.2E-2,
     %          1.4E-2,3.8E-2,1.4E-2,2.6E-1,1.7E-1,0.,
     %         6.2E-03,1.4E-02,1.1E-02,1.4E-02,7.8E-03,1.6E-02,2.3E-02,
     %          1.1E-02,1.1E-02,4.8E-02,1.8E-02,3.4E-01,2.1E-01,0.,
     %         7.3E-05,4.9E-03,9.1E-03,4.3E-02,2.5E-02,3.1E-02,1.3E-02,
     %          2.1E-02,3.5E-02,3.0E-02,3.1E-02,6.3E-02,2.3E-02,
     %          4.2E-01,2.5E-01,0.,
     %         .0001,.011,.034,.0685,.090,.0635,.0260,.017,
     %          .0210,.0350,.0540,.0713,.0960,.0850,.0170,
     %          .476,.297,0.,
     %         3.28E-4,5.84E-02,1.12E-01,1.32E-01,1.33E-01,1.26E-01,
     %          1.39E-01,9.55E-02,4.02E-01,4.19E-02,2.57E-02,4.45E-02,
     %          5.48E-02,7.13E-02,9.03E-02,1.10E-01,2.05E-02,5.49E-01,
     %          3.55E-01,0.,
     %         1.8E-3,3.6E-2,7.8E-2,2.2E-1,1.4E-1,1.4E-1,
     %          1.1E-1,6.3E-1,5.5E-1,3.6E-1,2.6E-1,1.6E-1,
     %          6.6E-2,2.5E-1,1.2E-1,5.0E+0,3.7E-2,6.3E-2,
     %          7.0E-2,1.1E-1,1.0E-1,1.1E-1,3.6E-2,7.5E-1,
     %          5.2E-1,0.,
     %         1.41E-03, 5.20E-03, 1.38E-02, 2.30E-02, 4.19E-02,
     %          6.83E-02, 1.22E-01, 3.00E-01, 1.50E-01, 6.97E-01,
     %          7.09E-01, 6.44E-01, 5.25E-01, 4.46E-01, 3.63E-01,
     %          3.02E-01, 1.02E-01, 2.70E-01, 4.67E-02, 8.35E-02,
     %          9.96E-02, 1.99E-01, 2.40E-01, 1.15E-01, 3.16E-02,
     %          8.03E-01, 5.75E-01, 0./
      data bdi/0.,
     %         0.3,0.,
     %         3.0,    0.5,    2.0,    0.2,    0.2,    0.,
     %         3.8,    4.1,    1.4,    3.0,    0.2,    0.2,    0.,
     %         2.5,    3.3,    6.0,    2.0,    3.2,    0.2,    0.2,
     %          0.,
     %         1.9,    0.6,    0.7,    4.3,    4.8,    1.6,    5.0,
     %          0.2,    0.2,    0.,
     %         0.,   0.,  3., 3.2,  3.2, 6.7,
     %          4.4,  3.5, 10., 0.2, 0.2, 0.,
     %         0.,     0.,     0.,     0.,     10.,    4.,     8.,
     %          6.3,    6.,     5.,     10.5,   0.2,    0.2,    0.,
     %         0.,     2.5,    6.0,    0.,     0.,     0.,     22.,
     %          6.4,    13.,    6.8,    6.3,    4.1,    12.,    0.2,
     %          0.2,    0.,
     %         .005, .045, .057, .087, .0769, .140, .120, .100, 1.92,
     %          1.66, 1.67, 1.40, 1.31, 1.02, .245, .294, .277, 0.,
     %         0.0907,.110,.0174,.132,.114,.162,.0878,.263,.0627,
     %          .0616,2.77,2.23,2.00,1.82,.424,.243,.185,.292,.275,
     %          0.,
     %         6*0.,1.3,4*0.4,  0.8,     2.7,     0.1,     1.9, 0.1,
     %          26.,    23.,    17.,    8.,     11.7,   15.4,   29.,
     %          0.3,    0.3,    0.,
     %         .469, .357, .281, .128, .0417, .0558, .0346, 0.,
     %          1.90, .277, .135, .134, .192, .332, .337, .121,
     %          .0514, .183, 7.56, 4.55, 4.87, 2.19, 1.15, 1.23,
     %          .132, .289, .286, 0./
      data t0/0.,
     %        47.,0.,
     %        11.,    15.,    9.1,    340.,   410.,   0.,
     %        13.,    14.,    18.,    11.,    470.,   540.,   0.,
     %        17.,    17.,    18.,    22.,    13.,    620.,   700.,
     %         0.,
     %        31.,    29.,    26.,    24.,    24.,    29.,    17.,
     %         980.,   1100.,  0.,
     %        5.1, 61., 44., 39., 34., 31.,
     %         31., 36., 21., 1400., 1500., 0.,
     %        11.,     12.,     10.,     120.,    55.,     49.,
     %         42.,     38.,     37.,     42.,     25.,     1900.,
     %         2000.,   0.,
     %        11.,     12.,     13.,     18.,     15.,     190.,
     %         67.,     59.,     55.,     47.,     42.,     50.,
     %         30.,     2400.,   2500.,   0.,
     %        32., 29., 23.9, 25.6, 25.0, 21.0, 18., 270., 83.,
     %         69.5, 60.5, 66.8, 65.0, 53.0, 35.5, 3010.,3130.,0.,
     %        3.46,38.5,40.8,38.2,35.3,31.9,32.2,24.7,22.9,373.,92.6,
     %         79.6,69.0,67.0,47.2,56.7,42.1,3650.,3780.,0.,
     %        5.8,    13.,     28.,     37.,     49.,     63.,
     %         68.,    77.,     73.,     71.,     68.,     61.,  59.,
     %         43.,    35.,    770.,   100.,   87.,    62.,    69.,
     %         68.,    67.,    41.,    5800.,  5900.,  0.,
     %        9.82, 20.1, 30.5, 42.0, 55.6, 67.2, 79.3, 90.0, 100.,
     %         78.1, 76.4, 74.4, 66.5, 59.7, 52.4, 49.6, 44.6, 849.,
     %         136., 123., 106., 125., 123., 33.2, 64.5, 6650.,
     %         6810., 0./
      data t1/0.,
     %        9.4,0.,
     %        4.9,    23.,    37.,    51.,    76.,    0.,
     %        4.8,    6.8,    38.,    59.,    72.,    98.,    0.,
     %        13.,    5.8,    9.1,    59.,    80.,    95.,    130.,
     %         0.,
     %        15.,    17.,    45.,    17.,    35.,    110.,   130.,
     %         140.,   260.,   0.,
     %        0.,  0.,  41., 87., 100., 54.,
     %         36., 160.,210.,240.,350., 0.,
     %        0.,     0.,     0.,     0.,     100.,   130.,   170.,
     %         60.,    110.,   250.,   280.,   310.,   440.,   0.,
     %        0.,      8.8,     15.,     0.,      0.,      0.,
     %         180.,    200.,    230.,    120.,    130.,    340.,
     %         360.,    460.,    550.,    0.,
     %        31., 55., 60., 38.1, 33., 21.5, 21.5, 330., 350.,
     %         360., 380., 290., 360., 280., 110., 605., 654., 0.,
     %        1.64,24.5,42.7,69.2,87.8,74.3,69.9,44.3,28.1,584.,
     %         489.,462.,452.,332.,137.,441.,227.,725.,768.,0.,
     %        6*0.,36.,63.,    85.,     89.,     100.,    120., 190.,
     %         190.,   250.,   90.,    630.,   770.,   620.,   510.,
     %         870.,   990.,   1000.,  980.,   1200.,  0.,
     %        10.1, 19.1, 23.2, 31.8, 45.5, 55.1, 52.8, 0.00, 55.0,
     %         88.7, 180.,125.,  189., 88.4, 129., 62.4, 159., 801.,
     %         932., 945., 945., 801., 757., 264., 193., 1190.,
     %         908., 0./
      data gli /2.,1.,2.,1.,6.,9.,4.,9.,6.,1.,2.,1.,6.,9.,4.,9.,6.,1.,
     *          2.,1./
      data gfe /2.,1.,2.,1.,6.,9.,4.,9.,6.,1.,2.,1.,6.,9.,4.,9.,
     *          6.,1.,10.,21.,28.,25.,6.,25.,30.,25./
      data gni /2.,1.,2.,1.,6.,9.,4.,9.,6.,1.,2.,1.,6.,9.,4.,9.,
     *          6.,1.,10.,21.,28.,25.,6.,25.,28.,21.,10.,21./
c
c parameters for calculating density dependent correction ap from Raymond
c
      DATA cdd/.24E-02,.1430E-01,.9094E-03,.3500E-01,.3050E-01,
     $  .9043E-02,.1077E-01,.2585E-03,.1953E-03,.8000E-01,
     $  .8715E-02,.1346E-01,.4753E-02,.6304E-02,.1601E-03,
     $  .1574E-03,.5600E-01,.1610E-01,.4081E-02,.7718E-02,
     $  .2910E-02,.4070E-02,.1059E-03,.1306E-03,.3370E-01,
     $  .1023E-01,.4726E-02,.3410E-02,.1660E-02,.3649E-02,
     $  .1412E-02,.2040E-02,.5280E-04,12*0.,.9555E-04,.8100E-01,
     $  .4168E-01,.2792E-01,.2585E-01,.2137E-02,.7325E-03,
     $  .8059E-03,.7821E-03,.6306E-03,.1501E-02,.5546E-03,
     $  .7711E-03,.1760E-04,.5965E-04,.6600E-01,.2842E-01,
     $  .1740E-01,.1579E-01,.1355E-01,.1221E-01,.1272E-02,
     $  .3673E-03,.4921E-03,.4976E-03,.4592E-03,.1108E-02,
     $  .3973E-03,.5326E-03,.1098E-04,.4948E-04,18*0.,20*0.,9*0.,
     $  .2030E-02,.2299E-02,.2313E-02,.2233E-02,.2734E-02,
     $  .2934E-02,.2319E-02,.3406E-03,.5245E-04,.1246E-03,
     $  .1320E-03,.1711E-03,.4206E-03,.1339E-03,.1461E-03,
     $  .1015E-05,.2508E-04,28*0./
c
      data dcfe/2.2e-4,2.3e-3,1.5e-2,3.8e-2,8.0e-2,9.2e-2,
     &     0.16,0.18,0.14,0.1,0.225,0.24,0.26,0.19,
     &     0.12,1.23,2.53e-3,5.67e-3,1.6e-2,1.85e-2,9.2e-4,
     &     0.131,1.1e-2,0.256,0.43,0.,1.e-4,2.7e-3,
     &     4.7e-3,1.6e-2,2.4e-2,4.1e-2,3.6e-2,0.07,0.26,
     &     0.28,0.231,0.17,0.16,0.09,0.12,0.,3.36e-2,
     &     7.82e-2,7.17e-2,9.53e-2,0.129,8.49e-2,4.88e-2,
     &     0.452,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     &     0.,0.,0.,0.,0.,0.,0.6,0.,0.181,3.18e-2,
     &     9.06e-2,7.9e-2,0.192,0.613,8.01e-2,0.,0.,0.,
     &     0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     &     0.,0.,0.,0.,1.92,1.26,0.739,1.23,0.912,0.,
     &     0.529,0.,0.,0./
      data defe/5.12,16.7,28.6,37.3,54.2,45.5,66.7,66.1,
     &     21.6,22.2,59.6,75.,36.3,39.4,24.6,560.,22.5,
     &     16.2,23.7,13.2,39.1,73.2,0.1,4.625e3,5.3e3,
     &     0.,12.9,31.4,52.1,67.4,100.,360.,123.,129.,
     &     136.,144.,362.,205.,193.,198.,248.,0.,117.,
     &     96.,85.1,66.6,80.3,316.,36.2,6.e3,0.,0.,
     &     0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     &     0.,0.,560.,0.,341.,330.,329.,297.,392.,
     &     877.,306.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     &     0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,683.,
     &     729.,787.,714.,919.,0.,928.,0.,0.,0./
c
      data istorey/5,6,7,11,12,13,14,18,19,20,21,22,0/
      data rstorey/ 0.0108,-0.1075, 0.2810,-0.0193,-0.1127,
     $              1.8267, 4.1012, 4.8443, 0.2261, 0.5960,
     $              2.3196,10.7328, 6.8830,-0.1824, 0.4101,
     $              0.0000, 0.6310, 0.1990,-0.0197, 0.4398,
     $              0.0320,-0.6624, 4.3191, 0.0003, 0.5946,
     $             -0.8806,11.2406,30.7066,-1.1721, 0.6127,
     $              0.4134,-4.6319,25.9172,-2.2290,-0.2360,
     $              0.0000, 0.0238, 0.0659, 0.0349, 0.5334,
     $             -0.0036, 0.7519, 1.5252,-0.0838, 0.2769,
     $              0.0000,21.8790,16.2730,-0.7020, 1.1899,
     $              0.0061, 0.2269,32.1419, 1.9939,-0.0646,
     $             -2.8425, 0.2283,40.4072,-3.4956, 1.7558,
     $               5*0./
c
      data uu/109.6787, 27*0.,
     *        198.3108,438.9089,26*0.,
     *        28*0.,
     *        28*0.,
     *        28*0.,
     *        90.82,196.665,386.241,520.178,3162.395,3952.06,
     *        22*0.,
     *        117.225,238.751,382.704,624.866,789.537,4452.758,
     *        5380.089,21*0.,
     *        109.837,283.24,443.086,624.384,918.657,1114.008,
     *        5963.135,7028.393,20*0.,
     *        28*0.,
     *        173.93,330.391,511.8,783.3,1018.,1273.8,1671.792,
     *        1928.462,9645.005,10986.876,18*0.,
     *        28*0.,
     *        61.671,121.268,646.41,881.1,1139.4,1504.3,1814.3,2144.7,
     *        2645.2,2964.4,14210.261,15829.951,16*0.,
     *        28*0.,
     *        65.748,131.838,270.139,364.093,1345.1,1653.9,1988.4,
     *        2445.3,2831.9,3237.8,3839.8,4222.4,19661.693,21560.63,
     *        14*0.,
     *        28*0.,
     *        83.558,188.2,280.9,381.541,586.2,710.184,2265.9,2647.4,
     *        3057.7,3606.1,4071.4,4554.3,5255.9,5703.6,26002.663,
     *        28182.535,12*0.,
     *        28*0.,
     *        127.11,222.848,328.6,482.4,605.1,734.04,1002.73,1157.08,
     *        3407.3,3860.9,4347.,4986.6,5533.8,6095.5,6894.2,7404.4,
     *        33237.173,35699.936,10*0.,
     *        28*0.,
     *        49.306,95.752,410.642,542.6,681.6,877.4,1026.,1187.6,
     *        1520.64,1704.047,4774.,5301.,5861.,6595.,7215.,7860.,
     *        8770.,9338.,41366.,44177.4,8*0.,
     *        28*0.,
     *        28*0.,
     *        28*0.,
     *        28*0.,
     *        28*0.,
     *        63.737,130.563,247.22,442.,605.,799.,1008.,1218.38,
     *        1884.,2114.,2341.,2668.,2912.,3163.,3686.,3946.82,
     *        10180.,10985.,11850.,12708.,13620.,14510.,15797.,
     *        16500.,71203.,74829.,2*0.,
     *        28*0.,
     *        61.6,146.542,283.8,443.,613.5,870.,1070.,1310.,1560.,
     *        1812.,2589.,2840.,3100.,3470.,3740.,4020.,4606.,
     *        4896.2,12430.,13290.,14160.,15280.,16220.,17190.,
     *        18510.,19351.,82984.,86909.4/
c     data hfrac/0.75/
      data hfrac/1.0/
      data ergsev/1.602192e-12/
      data cc1/1.e-06/
c
c Tim Kallman works with temperatures in units of 10^4 K
c
      t=temp/1.e4
c
      dirt=0.
      ini=inid(iont,iatom)
      if(ini.le.0) return
c
      ekt = t*(0.861707)
      xst = sqrt(t)
      hconst = hfrac*ekt*ergsev
      t3s2=1./(t*xst)
      tmr = 1.e-6*t3s2
      alogt = log10(t)
      kk = 0
      ist=1
      j=ini
         dirt = 0.
c        dr for iron from Arnaud and Raymond
         if ( j.lt.115 .or. j.gt.139 ) go to 2901
            kk = kk + 1
            do 20 n = 1,4
               dirt = dirt + dcfe(kk,n)*expo(-defe(kk,n)/ekt)
 20         continue
            dirt = dirt*tmr
            go to 101
 2901       continue
c           aldrovandi and Pequignot rates
c  The reference is Aldrovandi, S. M. V. and P\'equignot, D. (1973)
c  A&A, 25, 137
c
c     ap is the density dependent correction to dr from Raymond
c
            enn = xpx**(0.2)
            ap = 1./(1.+cdd(j)*enn)
            dirt = adi(j)*ap*cc1*expo(-t0(j)/t)
     &                *(1.+bdi(j)*expo(-t1(j)/t))
     $                /(t*sqrt(t))
            dirtemp=0.
c           storey dr rates
            if ((j.ne.(istorey(ist)-1)).or.(t.gt.6.).or.(ist.gt.12)) 
     $        go to 101
              dirtemp=
     $          (1.e-12)*(rstorey(1,ist)/t+rstorey(2,ist)
     $          +t*(rstorey(3,ist)+t*rstorey(4,ist)))*t3s2
     $          *expo(-rstorey(5,ist)/t)
              dirt=dirt+dirtemp
              ist=ist+1
 101          continue
c
c     pseudo cross-section
c
      if(iatom.le.20) then
         gp=gli(iont+1)
       if(gp.le.0) gp=1.
         gg=gp/gli(iont)
       else if(iatom.le.26) then
         gp=gfe(iont+1)
       if(gp.le.0) gp=1.
         gg=gp/gfe(iont)
       else if(iatom.le.28) then
         gp=gni(iont+1)
       if(gp.le.0) gp=1.
         gg=gp/gni(iont)
      end if
      frq0=cons*uu(iont,iatom)
      frq1=1.1*frq0
      delfr=frq1-frq0
      fra=0.5*(frq0+frq1)
      x=1.-expo(-4.79928e-11*delfr/temp)
      sig0=dirt*8.47272e24*gg*sqrt(temp)/fra**2/x
      return
      end
c
c
C     
C     
C ************************************************************************    
C     
c
      function expo(x)
c     ================
c     
      INCLUDE 'IMPLIC.FOR'
      crit=80.  
      if(x.lt.-crit) x=-crit
      if(x.gt.crit) x=crit
      expo=exp(x)    
      return
      end
c
C     
C     
C ************************************************************************    
C     
c
      SUBROUTINE IRC(N,T,IC,RNO,SE)
c     =============================
c
C IRC CALCULATES THE EXCITATION RATE, SE [cm**3/s], FOR IONIZATION
C OF HYDROGEN ATOMS FROM STATE N DUE TO ELECTRON COLLISIONS, ASSUMING
C THE CONTINUUM STARTS AT LEVEL RNO.  
C REF. JOHNSON (1972)
C
c     a modification of Tim Kallman's XSTAR routine
c
      INCLUDE 'IMPLIC.FOR'
      IF(IC.NE.1) THEN                          ! MAB
         CALL SZIRC(N,T,IC,rno,SE)
         RETURN
      END IF
c
      XO=1.-N*N/RNO/RNO
      YN=XO*157803./(T*N*N)
      IF(N-2) 100,200,300
 100   AN=1.9603*N*(1.133/3./XO**3-0.4059/4./XO**4+0.07014/5./XO**5)
      BN=2./3.*N*N/XO*(3.+2./XO-0.603/XO/XO)
      RN=0.45
      GO TO 400
C
 200   AN=1.9603*N*(1.0785/3./XO**3-0.2319/4./XO**4+0.02947/5./XO**5)
      BN=(4.-18.63/N+36.24/(N*N)-28.09/(N*N*N))/N
      BN=2./3.*N*N/XO*(3.+2./XO+BN/XO/XO)
      RN=0.653
      GO TO 400
C
 300   G0=(0.9935+0.2328/N-0.1296/(N*N))/3./XO**3
      G1=-(0.6282-0.5598/N+0.5299/(N*N))/(N*4.)/XO**4
      G2=(0.3887-1.181/N+1.470/(N*N))/(N*N*5.)/XO**5
      AN=1.9603*N*(G0+G1+G2)
      BN=(4.-18.63/N+36.24/(N*N)-28.09/(N*N*N))/N
      BN=(3.+2./XO+BN/XO/XO)*2.*N*N/3./XO
      RN=1.94*N**(-1.57)
C
 400   CONTINUE
      RN=RN*XO
      ZN=RN+YN
      CALL EXPINX(YN,EY)
      CALL EXPINX(ZN,EZ)
      SE=AN*(EY/YN/YN-EXP(-RN)*EZ/ZN/ZN)
      EY=1.+1./YN-EY*(2./YN+1.)
      EZ=EXP(-RN)*(1.+1./ZN-EZ*(2./ZN+1.))
      SE=SE+(BN-AN*LOG(2.*N*N/XO))*(EY-EZ)
      SE=SE*SQRT(T)*YN*YN*N*N*1.095E-10/XO
      RETURN
      END
c
C     
C     
C ************************************************************************    
C     
c
      subroutine szirc(nn,T,ic,rno,cii)
c     =================================
c
c     calculates electron impact ionizition rates from semiempirical
c     formula (eq.35) from Sampson & Zhang (1988, ApJ 335, 516)
c
c     a modification of Tim Kallman's XSTAR routine
c
      INCLUDE 'IMPLIC.FOR'
       real abethe(11), hbethe(11), rbethe(11)
       DATA(abethe(i),i=1,11)/ 1.134, 0.603, 0.412, 0.313, 0.252,
     1       0.211, 0.181, 0.159, 0.142, 0.128, 1.307 /
       DATA(hbethe(i),i=1,11)/ 1.48, 3.64, 5.93, 8.32, 10.75, 12.90,
     1       15.05, 17.20, 19.35, 21.50, 2.15 /
       DATA(rbethe(i),i=1,11)/ 2.20, 1.90, 1.73, 1.65, 1.60, 1.56,
     1       1.54, 1.52, 1.52, 1.52, 1.52 /
       rz=ic
       Boltz=1.38066e-16
       Eion=2.179874e-11
       const=4.6513e-3
C
       rc=float(int(rno))
       if (nn.lt.11) then
         an=abethe(nn)
         hn=hbethe(nn)
         rrn=rbethe(nn)
       else
         an=abethe(11)/float(nn)
         hn=hbethe(11)*float(nn)
         rrn=rbethe(11)
       endif
       tt= T*Boltz
       rn=float(nn)
       yy=rz*rz*Eion/tt*(1./rn/rn-1./rc/rc-.25*(1./(rc-1.)**2-
     c    1./rc/rc))
       call eint(yy,e1,e2,e3)
       cii=const*sqrt(tt)*(rn**5)/(rz**4)*an*yy* (
     1   e1/rn-(exp(-yy)-yy*e3)/(3.*rn)+(yy*e2-2.*yy*e1+exp(-yy))*
     2   3.*hn/rn/(3.-rrn)+(e1-e2)*3.36*yy)
       return
       end
C     
C     
C ************************************************************************    
C     
c
      subroutine expinx(x,em1)
c     ========================
c
c expinx is a subroutine to calculate the value of e1, the exponential
c integral or em1=x*expo(x)*e1 at the point x.  the polynomial
c expressions that are used come from abromowitz and stegen
c
c     a modification of Tim Kallman's XSTAR routine
c
c
      INCLUDE 'IMPLIC.FOR'
      if(x.le.1.) go to 100
c
      b1=9.5733223454
      b2=25.6329561486
      b3=21.0996530827
      b4=3.9584969228
      c1=8.5733287401
      c2=18.0590169730
      c3=8.6347608925
      c4=0.2677737343
      em1=x**4+c1*x**3+c2*x*x+c3*x+c4
      em1=em1/(x**4+b1*x*x*x+b2*x*x+b3*x+b4)
      go to 200
c
 100  continue
      a0=-0.57721566
      a1=0.99999193
      a2=-0.24991055
      a3=0.05519968
      a4=-0.00976004
      a5=0.00107857
      if (x.gt.0)then
      e1= a0+a1*x+a2*x*x+a3*x**3+a4*x**4+a5*x**5-log(x)
      else
      e1=-a0+a1*x+a2*x*x+a3*x**3+a4*x**4+a5*x**5-log(-x)
      endif
      em1=e1*x*expo(x)
c
 200  continue
      return
      end
C     
C     
C ************************************************************************    
C     
       subroutine eint(t,e1,e2,e3)
c     ============================
c
c  returns the values of the exponential integral function of order
c  1, 2, and 3
c
c     a modification of Tim Kallman's XSTAR routine
c
      INCLUDE 'IMPLIC.FOR'
       e1=0.
       e2=0.
       e3=0.
       call expinx(t,ss)
c       ss=expint(x)
       e1=ss/t/expo(t)
       e2=exp(-t)-t*e1
       e3=0.5*(expo(-t)-t*e2)
       return
       end
c
c
C     
C     
C ************************************************************************    
C     
C     
      SUBROUTINE COMSET
C     =================
C
C     sets up necessary parameters for treating the Compton scattering
c
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'MODELQ.FOR'
      dimension freqi(mfreq)
      parameter (xcon=8.0935d-21,YCON=1.68638E-10)
      parameter (t15=1.d-15)
      common/auxcbc/cden1m(mdepth),cden10(mdepth),
     *              cden2m(mdepth),cden20(mdepth)
      common/comgfs/gfm(mfreq,mdeptc),gfp(mfreq,mdeptc)
      DIMENSION PL(MDEPTH),PLM(MDEPTH)
C
      if(icompt.le.0) go to 100
      nmuc=3
      nsti=0
      nedd=3
C
C     frequency-dependent universal parameters
C
      do ij=1,nfreq
         cder10(ij)=0.
         cder1p(ij)=0.
         cder1m(ij)=0.
         cder20(ij)=0.
         cder2p(ij)=0.
         cder2m(ij)=0.
       iji=nfreq-kij(ij)+1
       ijorig(iji)=ij
       freqi(iji)=freq(ij)
       fr=freqi(iji)
       bnus(iji)=two*xcon*fr/(bn*(fr*1.d-15)**3)
      end do
C
      ij=1
      dlnfr(ij)=log(freqi(ij+1)/freqi(ij))
      do ij=2,nfreq-1
       dlnfr(ij)=log(freqi(ij+1)/freqi(ij))
         delp=dlnfr(ij)
       delm=dlnfr(ij-1)
       del0=delp+delm
       cd0=two/del0
       cder2m(ij)=cd0/delm
       cder2p(ij)=cd0/delp
       cder20(ij)=-cder2m(ij)-cder2p(ij)
      end do
c
      do ij=1,nfreq-1
       frj0=freqi(ij)
       frjp=freqi(ij+1)
       frz=sqrt(frj0*frjp)
       do id=1,nd
C           to avoid over/underflow problems:
            IF(HK*FRJ0/TEMP(ID).LT.200.) THEN
            fjb0=un/(exp(hk*frj0/temp(id))-un)
            ELSE
              fjb0=0.
            ENDIF
            IF(HK*FRJP/TEMP(ID).LT.200.) THEN
            fjbp=un/(exp(hk*frjp/temp(id))-un)
            ELSE
              fjbp=0.
            ENDIF
          fjz0=fjb0*(bn*(frj0*t15)**3)
          fjzp=fjbp*(bn*(frjp*t15)**3)
          if(ichcoo.eq.0) then
             zj0=hk*frz/temp(id)
             dfjz=fjz0-fjzp
             dfjb=fjb0-fjbp
             fzz=un+fjbp-3./zj0
             aa=dfjz*dfjb
             bb=dfjz*fzz+fjzp*dfjb
             cc=fjzp*fzz-dfjz/dlnfr(ij)/zj0
           else
               e2=ycon*temp(id)
               zxxp=xcon*frjp*(un+fjbp)-3.*e2
               zxx0=xcon*frj0*(un+fjb0)-3.*e2
             dzxx=zxx0-zxxp
             dfjb=fjb0-fjbp
               dfjz=fjz0-fjzp
             aa=dfjz*dzxx
             bb=dfjz*zxxp+fjzp*dzxx
             cc=fjzp*zxxp-e2*dfjz/dlnfr(ij)
          end if
CXXX        to avoid division by zero:
          if(abs(aa).eq.0.and.abs(bb).eq.0.) then
              xx1=0.
            elseif(abs(aa).lt.1.e-7*abs(bb)) then
             xx1=-cc/bb
          else
             dd=bb*bb-4.*aa*cc
             if(dd.lt.0.) dd=0.
             dd=sqrt(dd)
             xx1=(dd-bb)*half/aa
             if(ichcoo.gt.0) then
                  xx2=-(dd+bb)*half/aa
                  dxx1=abs(xx1-half)
                  dxx2=abs(xx2-half)
                  if(dxx2.lt.dxx1) xx1=xx2
                  if((xx1.gt.1.).or.(xx1.lt.0.)) xx1=half
             end if
            end if
          delj(ij,id)=xx1
         end do
      end do
c
C     angle-dependent universal parameters
C
      call angset
c
C     frequency-dependent universal parameters
C
  100 continue
      do ij=1,nfreq
c
c  first-order expression
c
         if(knish.eq.0) then
            SIGEC(IJ)=SIGE*(un-two*freq(ij)*xcon)
c
C  Use full Klein-Nishina cross section (Rybicki & Lightman 1975):
c
          else
            xf=xcon*freq(ij)
            if(xf.lt.1.d-1) then
               SIGEC(IJ)=SIGE*(1.-xf*(2.-xf*(26./5.-xf*(13.3
     *                  -xf*(1144./35.-xf*(544./7.-xf*(3784./21.
     *                  -xf*(6148./15.-xf*(151552./165.
     *                  -xf*111872./55.)))))))))
              else if(xf.gt.1.d3) then
                 SIGEC(IJ)=SIGE*3./8./xf*(log(2.*xf)+0.5)
              else
                 SIGEC(IJ)=SIGE*0.75*((1.+xf)/xf**3*(2.*xf*(1.+xf)/
     *           (1.+2.*xf)-log(1.+2.*xf))+0.5*log(1.+2.*xf)/xf
     *           -(1.+3.*xf)/(1.+2.*xf)**2)
            endif
         endif
      end do
c
      if(icompt.le.0) return
      IJ=1
      IJO=ijorig(ij)
      DO ID=1,ND
         PLM(ID)=BNUE(IJO)/(EXP(HK/Temp(ID)*FREQ(IJO))-UN)
C        GFM(IJ,ID)=0.
C        GFP(IJO,ID)=0.
      END DO
C      
      DO IJ=2,NFREQ
       IJO=ijorig(ij)
         DO ID=1,ND
C           to avoid over/underflow problems:
            IF(HK/TEMP(ID)*FREQ(IJO).LT.200.) THEN
               PL(ID)=BNUE(IJO)/(EXP(HK/temp(ID)*FREQ(IJO))-UN)
            ELSE
               PL(ID)=PLM(ID)
            ENDIF
c            GFM(IJ,ID)=PLM(ID)/PL(ID)
c            GFP(IJ-1,ID)=PL(ID)/PLM(ID)
            PLM(ID)=PL(ID)
         END DO
      END DO
C      
      return
      end
C
C
C     ******************************************************************
C
c
c
      subroutine angset
c     =================
c
c     sets up angles points and angle-dependent quantities for treating 
c     the Compton scattering
c
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      parameter(three=3.d0, five=5.d0, zero=0.d0, tr16=3.d0/16.d0)
      dimension amu0(mmuc),wtmu0(mmuc)
c    
c     amu=cos(angle between line of sight and normal to slab) grid and
c     gauss-legendre integration weights for the interval mu=[0,1]
c
      call gauleg(zero,un,amu0,wtmu0,nmuc,mmuc)
c
      do i=1,nmuc
         amuc(i)=-amu0(nmuc-i+1)
         amuc(i+nmuc)=amu0(i)
         wtmuc(i)=wtmu0(nmuc-i+1)
         wtmuc(i+nmuc)=wtmu0(i)
      end do
      nmuc=2*nmuc
c
      do i=1,nmuc
         amuc1(i)=amuc(i)*wtmuc(i)
         amuc2(i)=amuc(i)*amuc(i)*wtmuc(i)
         amuc3(i)=amuc(i)*amuc(i)*amuc(i)*wtmuc(i)
         a1=amuc(i)
         a2=a1*a1
         a3=a1*a2
         do i1=1,nmuc
            b1=amuc(i1)
            b2=b1*b1
            b3=b1*b2
            trw=tr16*wtmuc(i1)
            calph(i,i1)=(three*a2*b2-a2-b2+three)*trw
            cbeta(i,i1)=(five*(a1*b1+a3*b3)-three*(a3*b1+a1*b3))*trw
            cgamm(i,i1)=a1*b1*trw
         end do
      end do
c  
      return
      end

C
C
C *********************************************************************
C
C

      SUBROUTINE GAULEG(X1,X2,X,W,N,M)
C     ================================
C
C     set up angle points
C
      INCLUDE 'IMPLIC.FOR'
      DIMENSION X(M),W(M)
      PARAMETER (EPS=3.D-14,half=0.5d0,pi=3.141592654d0,quart=0.25,
     *            un=1.d0, two=2.d0)
c
      N2=(N+1)/2
      XM=HALF*(X2+X1)
      XL=HALF*(X2-X1)
      DO 12 I=1,N2
         Z=COS(PI*(I-quart)/(N+half))
1        CONTINUE
         P1=1.D0
         P2=0.D0
         DO 11 J=1,N
            P3=P2
            P2=P1
            P1=((TWO*J-UN)*Z*P2-(J-UN)*P3)/J
11        CONTINUE
          PP=N*(Z*P1-P2)/(Z*Z-un)
          Z1=Z
          Z=Z1-P1/PP
        IF(ABS(Z-Z1).GT.EPS) GO TO 1
        X(I)=XM-XL*Z
        X(N+1-I)=XM+XL*Z
        W(I)=TWO*XL/((UN-Z*Z)*PP*PP)
        W(N+1-I)=W(I)
12    CONTINUE
      RETURN
      END
C
C
C
C     ****************************************************************
C
C
      subroutine rte_sc(dtau,st0,rup,rdown,amu0,ri,ali)
C     ================================================
C
C     formal solver of the radiative transfer equation 
C     for one frequency, angle, and for completely known source function;
c     using short characteristics
C
c
c     input:    dtau - optical depth increments Delta tau
c               st0  - total source function
c               rup  - intensity at the upper boundary (id=1)
c               rdown- intensity at the lower boundary (id=nd)
c               amu0 - cosine of angle of propagation (wrt. the normal)
c     output:   ri   - radiation intensity
c
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      dimension dtau(mdepth),st0(mdepth),ri(mdepth),ali(mdepth),
     *          dtx1(mdepth),dtx2(mdepth),dtx0(mdepth)
C  
      do id=1,nd-1
       dtx1(id)=exp(-dtau(id))
       dtx2(id)=(un-dtx1(id))/dtau(id)
         dtx0(id)=un-dtx2(id)
      enddo
c               
c     incoming intensity 
c 
      if(amu0.lt.0) then
c  
         ID=1
         ri(id)=rup
         do id=1,nd-1
          ri(id+1)=ri(id)*dtx1(id)+st0(id)*(dtx2(id)-dtx1(id))+
     *               st0(id+1)*dtx0(id)
            ali(id+1)=dtx0(id)
         enddo
       ali(1)=0.  
C
c     outgoing intensity
c
       else
c
         ri(nd)=rdown
         do id=nd-1,1,-1
          ri(id)=ri(id+1)*dtx1(id)+st0(id)*dtx0(id)+
     *             st0(id+1)*(dtx2(id)-dtx1(id))
            ali(id)=dtx0(id)
         enddo
       ali(nd)=0.
      end if
      return
      end
C
C
C
C     ****************************************************************
C
C
      subroutine rtesol(dtau,st0,rup,rdown,amu0,ri,ali)
C     ================================================
C
C     formal solver of the radiative transfer equation 
C     for one frequency, angle, and for completely known source function;
c     by the Discontinuous Finite Element method
c     Castor, Dykema, Klein, 1992, ApJ 387, 561.
c
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      parameter (one=1.d0)
      dimension dtau(mdepth),st0(mdepth),ri(mdepth),ali(mdepth),
     *          rim(mdepth),rip(mdepth),aim(mdepth),aip(mdepth)
c   
c     incoming intensity 
c   
      if(amu0.lt.0) then
c   
            id=1
            rip(id)=rup
            dt0=dtau(id)
            dtaup1=dt0+one
            dtau2=dt0*dt0
            bb=two*dtaup1
            cc=dt0*dtaup1
            aa=dtau2+bb
            rim(id)=(aa*rip(id)-cc*st0(id)+dt0*st0(id+1))/bb
            do id=1,nd-1
               dt0=dtau(id)
               dtaup1=dt0+one
               dtau2=dt0*dt0
               bb=two*dtaup1
               cc=dt0*dtaup1
               aa=dtau2+bb
               rim(id+1)=(two*rim(id)+dt0*st0(id)+cc*st0(id+1))/aa
               rip(id)=(bb*rim(id)+cc*st0(id)-dt0*st0(id+1))/aa
               aim(id+1)=cc/aa
               aip(id)=(cc+bb*aim(id))/aa
            enddo
            do id=2,nd-1
               dtt=un/(dtau(id-1)+dtau(id))
               ri(id)=(rim(id)*dtau(id)+rip(id)*dtau(id-1))*dtt
               ali(id)=(aim(id)*dtau(id)+aip(id)*dtau(id-1))*dtt
            enddo
            ri(1)=rip(1)
            ri(nd)=rim(nd)
            ali(1)=aim(1)
            ali(nd)=aim(nd)
C
c     outgoing intensity
c
       else
c
            rip(nd)=rdown
            id=nd-1
            dt0=dtau(id)
            dtaup1=dt0+one
            dtau2=dt0*dt0
            bb=two*dtaup1
            cc=dt0*dtaup1
            aa=dtau2+bb
            rim(id+1)=(aa*rip(id+1)-cc*st0(id+1)+dt0*st0(id))/bb
            do id=nd-1,1,-1
               dt0=dtau(id)
               dtaup1=dt0+one
               dtau2=dt0*dt0
               bb=two*dtaup1
               cc=dt0*dtaup1
               aa=dtau2+bb
               rim(id)=(two*rim(id+1)+dt0*st0(id+1)+cc*st0(id))/aa
               rip(id+1)=(bb*rim(id+1)+cc*st0(id+1)-dt0*st0(id))/aa
               aim(id)=cc/aa
               aip(id+1)=(cc+bb*aim(id+1))/aa
            enddo
            do id=2,nd-1
               dtt=un/(dtau(id-1)+dtau(id))
               ri(id)=(rim(id)*dtau(id-1)+rip(id)*dtau(id))*dtt
               ali(id)=(aim(id)*dtau(id-1)+aip(id)*dtau(id))*dtt
            enddo
            ri(1)=rim(1)
            ri(nd)=rip(nd)
            ali(1)=aim(1)
            ali(nd)=aim(nd)
      end if
c
      return
      end
C
C
C     ****************************************************************
C
C
      subroutine rtefe2(dtau,s,rup,rdown,ri)
C     ======================================
C
C     formal solver of the radiative transfer equation 
C     for one frequency, angle, and for completely known source function;
c     original Feautrier (second-order) scheme
c
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      parameter (one=1.d0)
      dimension dtau(mdepth),s(mdepth),ri(mdepth)
      dimension a(mdepth),b(mdepth),c(mdepth),d(mdepth),
     *          f(mdepth),v(mdepth),z(mdepth)
c
c     set up the global tridiagonal matrix
c               
c     upper boundary condition 
c                 
      id=1
      cc=two/dtau(id)
      c(id)=cc/dtau(id)
      b(id)=one+cc+c(id)
      a(id)=0.
      v(id)=s(id)+cc*rup
c               
c     normal depth points 
c              
      do id=2,nd-1
         dtinv=two/(dtau(id-1)+dtau(id))  
         a(id)=dtinv/dtau(id-1)
         c(id)=dtinv/dtau(id)  
         b(id)=one+a(id)+c(id)
         v(id)=s(id)
      enddo
c               
c     lower boundary condition 
c              
      id=nd
      aa=two/dtau(id-1)
      a(id)=aa/dtau(id-1)
      b(id)=one+aa+a(id)
      if(rdown.eq.0.) b(id)=one+a(id)
      c(id)=0.
      v(id)=s(id)+aa*rdown
c               
c     ---------------------------------------------------
c     solution by elimination
c     1. forward sweep
c               
      f(1)=(b(1)-c(1))/c(1)
      d(1)=one/(one+f(1))
      z(1)=v(1)/b(1)
c 
c     ii) normal depth points
c
      do id=2,nd-1
         f(id)=(b(id)-a(id)-c(id)+a(id)*f(id-1)*d(id-1))/c(id)
         d(id)=one/(one+f(id))
         z(id)=(v(id)+a(id)*z(id-1))*d(id)/c(id)
      enddo
c
c     iii) upper boundary
c
      id=nd
      z(id)=(v(id)+a(id)*z(id-1))/(b(id)-a(id)*d(id-1))
c
c     2. backward elimination 
c   
      ri(nd)=z(nd)
      do id=nd-1,1,-1
         ri(id)=ri(id+1)*d(id)+z(id)
      enddo
c
      return
      end
C
C
C
C     ****************************************************************
C
C

      SUBROUTINE RTECF0(IJ)
C     =====================
C
C     Setup of the individual matrix elements of matrices A,B,C, E,U,V,
C     and alpha, beta, gamma, for solving the coupled transfer equation
C     with Compton scattering
C     Evaluation for a given frequency point IJ.
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'MODELQ.FOR'
      INCLUDE 'ALIPAR.FOR'
      INCLUDE 'ITERAT.FOR'
      COMMON/OPTDPT/DT(MDEPTH)
      PARAMETER (XCON=8.0935D-21,YCON=1.68638E-10)
      COMMON/AUXRTE/
     *       COMA(MDEPTH),COMB(MDEPTH),COMC(MDEPTH),VL(MDEPTH),
     *       COME(MDEPTH),U(MDEPTH),V(MDEPTH),BS(MDEPTH),
     *       AL(MDEPTH),BE(MDEPTH),GA(MDEPTH)
      common/auxcbc/cden1m(mdepth),cden10(mdepth),
     *       cden2m(mdepth),cden20(mdepth)
c
      IJI=NFREQ-KIJ(IJ)+1
      FR=FREQ(IJ)
      frp=freq(ijorig(iji+1))
      frm=freq(ijorig(iji-1))
      xcomp=fr*xcon
C
C     optical depth scale 
C   
      do id=1,nd-1
         dt(id)=deldmz(id)*(absot(id+1)+absot(id))
      end do
C
C     depth discretization matrices
C 
c     1. upper boundary
c  
      id=1
      dtp1=dt(id)
      bb0=un/dtp1
      bb1=two*bb0*bb0
      be(id)=bb0*two*fh(ij)+bb1*fak(ij,id)
      ga(id)=bb1*fak(ij,id+1)
      sext=two*bb0*hextrd(ij)
c
c     2. normal depth point
c  
      do id=2,nd-1
         dtm1=dtp1
         dtp1=dt(id)
         dt0=two/(dtm1+dtp1)
         al(id)=fak(ij,id-1)/dtm1*dt0
         ga(id)=fak(ij,id+1)/dtp1*dt0
         be(id)=fak(ij,id)*dt0*(un/dtm1+un/dtp1)
      end do
c
c     3. lower boundary
c  
      id=nd
c
c     stellar atmospheric
c
      IF(IDISK.EQ.0.OR.IFZ0.LT.0) then
         IF(IBC.EQ.0) THEN
            be(ID)=fak(ij,id)/DTP1+HALF
            al(ID)=fak(ij,id-1)/DTP1
          ELSE IF(IBC.LT.4) THEN
            B=UN/DTP1
            A=TWO*B*B
            be(id)=B*TWO*fhd(ij)+a*fak(ij,id)
            al(id)=a*fak(ij,id-1)
          ELSE 
            B=UN/DTP1
            A=TWO*B*B
            be(id)=b+a*fak(ij,id)
            al(id)=a*fak(ij,id-1)
         END IF
c
c     accretion disk - symmetric boundary
c
       ELSE
         bb0=un/dtp1
         bb1=two*bb0*bb0
         B=TWO/DTP1
         be(id)=bb1*fak(ij,id)
         al(id)=bb1*fak(ij,id-1)
      END IF
C
C     scattering matrices
C   
      do id=1,nd
         scat0=elec(id)*sige
         sa0=emis1(id)/abso1(id)
         ss0=scat0/abso1(id)
         epsnu=(abso1(id)-scat1(id))/abso1(id)
         x0=ss0
         e2=ycon*temp(id)+0.7*xcomp*xcomp
         e1=xcomp-3.*e2-0.7*xcomp*xcomp
         e0=1.-xcomp-4.2*xcomp*xcomp
         coma(id)=0.
         comc(id)=0.
         u(id)=0.
         v(id)=0.
         vl(id)=sa0
         if(id.eq.1.and.iwinbl.lt.0) vl(id)=sa0+sext
         bs(id)=0.
         if(iji.eq.1) then
            comc(id)=0.
            comb(id)=x0*(un-2.*xcomp)
          else if(iji.lt.nfreq) then
            del0=two/(dlnfr(iji)+dlnfr(iji-1))
            cder1p(iji)=(un-delj(iji,id))*del0
            cder1m(iji)=-delj(iji-1,id)*del0
            if(ichcoo.eq.0) then
               cder10(iji)=-cder1m(iji)-cder1p(iji)
               coma(id)=x0*(e1*cder1m(iji)+e2*cder2m(iji))
               comb(id)=x0*(e0+e1*cder10(iji)+e2*cder20(iji))
               comc(id)=x0*(e1*cder1p(iji)+e2*cder2p(iji))
               x0=ss0*bnus(iji)
               IF(ICOMST.EQ.0) X0=0.
               come(id)=x0*(cder10(iji)-un)
               u(id)=x0*cder1m(iji)
               v(id)=x0*cder1p(iji)
               bs(id)=come(id)*rad(iji,id)+
     *                u(id)*rad(iji-1,id)+v(id)*rad(iji+1,id)
             else
               cder10(iji)=-del0*(un-delj(iji-1,id)-delj(iji,id))
               zxxp=xcon*frp+0.5*bnus(iji+1)*rad(iji+1,id)-3.*e2
               zxx0=xcomp+0.5*bnus(iji)*rad(iji,id)-3.*e2
               zxxm=xcon*frm+0.5*bnus(iji-1)*rad(iji-1,id)-3.*e2
               zxxp12=((un-delj(iji,id))*zxxp+delj(iji,id)*zxx0)*del0
               zxxm12=((un-delj(iji-1,id))*zxx0+delj(iji-1,id)*zxxm)*
     *                del0
               coma(id)=x0*(-delj(iji-1,id)*zxxm12+e2*cder2m(iji))
               comc(id)=x0*((un-delj(iji,id))*zxxp12+e2*cder2p(iji))
               comb(id)=x0*(delj(iji,id)*zxxp12-(un-delj(iji-1,id))*
     *                  zxxm12+e2*cder20(iji))-epsnu+1
            end if
          else
            dlt=delj(iji-1,id)
            zj1=exp(-hk*freq(ij)/temp(id))
            zj2=exp(-hk*freq(ij+1)/temp(id))
            if(ichcoo.eq.0) then
               zj0=un/(hk*sqrt(freq(ij)*freq(ij+1))/temp(id))
               zxx=un-3.*zj0+(un-dlt)*zj1+dlt*zj2
               comb(id)=zj0/dlnfr(iji-1)+(un-dlt)*zxx
               coma(id)=-zj0/dlnfr(iji-1)+dlt*zxx
             else
               zxx0=xcomp*(un+zj1)-3.*e2
               zxxm=xcon*frm*(un+zj2)-3.*e2
               zxx=(un-dlt)*zxx0+dlt*zxxm
               comb(id)=e2/dlnfr(iji-1)+(un-dlt)*zxx
               coma(id)=-e2/dlnfr(iji-1)+dlt*zxx
            end if
            vl(id)=0.
            if(icomde.ne.0) then
               al(id)=0.
               be(id)=-un
               ga(id)=0.
            end if
         end if
         if(icomde.eq.0) then
            coma(id)=0.
            comc(id)=0.
            comb(id)=x0*(un-2.*xcomp)
         end if
      end do
C     
      return
      end
C
C
C     ****************************************************************
C
C
      SUBROUTINE INICOM
C     =================
C
C     Auxiliary procedure for INILAM
C     initialization of g-factors for the Compton scattering
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'ATOMIC.FOR'
      INCLUDE 'MODELQ.FOR'
      common/comgfs/gfm(mfreq,mdeptc),gfp(mfreq,mdeptc)
      DIMENSION PL(MDEPTH),PLM(MDEPTH)
C
      IJ=1
      IJO=ijorig(ij)
      DO ID=1,ND
         PLM(ID)=BNUE(IJO)/(EXP(HKT1(ID)*FREQ(IJO))-UN)
      END DO
C      
      DO IJ=2,NFREQ
         IJO=ijorig(ij)
         DO ID=1,ND
            PL(ID)=BNUE(IJO)/(EXP(HKT1(ID)*FREQ(IJO))-UN)
            PLM(ID)=PL(ID)
         END DO
      END DO
C      
      RETURN
      END
C
C
C     ****************************************************************
C
C

      SUBROUTINE RTECOM
C     =================
C
C     Solution of the radiative transfer equation with Compton scattering
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'MODELQ.FOR'
      INCLUDE 'ALIPAR.FOR'
      INCLUDE 'ITERAT.FOR'
      COMMON/OPTDPT/DT(MDEPTH)
      COMMON/AUXRTE/
     *          COMA(MDEPTH),COMB(MDEPTH),COMC(MDEPTH),VL(MDEPTH),
     *          COME(MDEPTH),U(MDEPTH),V(MDEPTH),BS(MDEPTH),
     *          AL(MDEPTH),BE(MDEPTH),GA(MDEPTH)
      common/comgfs/gfm(mfreq,mdeptc),gfp(mfreq,mdeptc)
      dimension aa(mdepth),bb(mdepth),cc(mdepth),
     *          d(mdepth),f(mdepth),z(mdepth),rd(mdepth)
c           
      PRD0=0.
      DO ID=1,ND
         PRADT(ID)=0.
      END DO
c           
c     ---------------------------------------------
c     "1st formal solution" to update Eddington factors 
c     ---------------------------------------------
c 
      if(ncfor1.gt.0) then
      do iform=1,ncfor1
         ij0=1
         if(icombc.gt.0) ij0=2
         DO IJ=ij0,NFREQ
            CALL OPACF1(IJ)
            CALL RTECF1(IJ)
         END DO
         if(icombc.gt.0) then
            ij=1
            iji=nfreq
            call rtecf0(ij)
            do id=1,nd
               rad(iji,id)=-rad(iji-1,id)*coma(id)/(comb(id)+bs(id))
            end do
         end if
      end do
      end if
c           
c     -----------------------------------------------
c     coupled solution for the frequency derivatives terms 
c     -----------------------------------------------
c
c     fully coupled treatment - traditional formulation
c
      if(ncfull.gt.0) then
      do icfull=1,ncfull
         call rtecmc 
c
c        iterative treatment of the derivative terms
c
         if(ncitot.gt.0) then
         do ictot=1,ncitot
c 
            if(nccoup.gt.0) then
            do iccoup=1,nccoup
               do ij=1,nfreq
                  ijo=ijorig(ij)
                  fr=freq(ijo)
                  call opacf1(ijo)
                  call rtecf0(ijo)
                  do id=1,nd
                     comb(id)=comb(id)+bs(id)
                     bb(id)=be(id)+un-comb(id)
                     aa(id)=al(id)
                     cc(id)=ga(id)
                     vl(id)=vl(id)+
     *                      (coma(id)*gfm(ij,id)+
     *                      comc(id)*gfp(ij,id))*rad(ij,id)
                  end do
c
c                 ----------------
c                 forward sweep
c                 ----------------
c
c                 i) upper boundary
c
                  f(1)=(bb(1)-cc(1))/cc(1)
                  d(1)=un/(un+f(1))
                  z(1)=vl(1)/bb(1)
c 
c                 ii) normal depth points
c
                  do id=2,nd-1
                     f(id)=(bb(id)-aa(id)-cc(id)+
     *                     aa(id)*f(id-1)*d(id-1))/cc(id)
                     d(id)=un/(un+f(id))
                     z(id)=(vl(id)+aa(id)*z(id-1))*d(id)/cc(id)
                  end do
c
c                 iii) lower boundary
c
                  id=nd
                  z(id)=(vl(id)+aa(id)*z(id-1))/(bb(id)-aa(id)*d(id-1))
c
c                 --------------------
c                 backward elimination 
c                 --------------------
c   
                  rd(nd)=z(nd)
                  do id=nd-1,1,-1
                     rd(id)=rd(id+1)*d(id)+z(id)
                  end do
c   
                  do id=1,nd
                     rad(ij,id)=rd(id)
                  end do
               end do
            end do
            end if
c           
c           ---------------------------------------------
c           "2nd formal solution" to update Eddington factors 
c           ---------------------------------------------
c 
            if(ncfor2.gt.0) then
            do iform=1,ncfor2
               ij0=nfreq
               if(icombc.gt.0) ij0=nfreq-1
               PRD0=0.
               DO ID=1,ND
                  PRADT(ID)=0.
               END DO
               DO IJ=1,ij0
                  ijo=ijorig(ij)
                  CALL OPACF1(IJo)
                  CALL RTECF1(IJo)
               END DO
               PRD0=PRD0*PCK
               DO ID=1,ND
                  PRADT(ID)=PRADT(ID)*PCK
               END DO
               if(icombc.gt.0) then
                  ij=1
                  iji=nfreq
                  call rtecf0(ij)
c                 do id=1,nd
c                    radcm(iji,id)=-radcm(iji-1,id)*coma(id)/
c    *               (comb(id)+bs(id))
c                 end do
c                 flux(1)=radcm(iji,1)*fh(2)
               end if
c   
c              do id=1,nd
c                 DO IJ=1,NFREQ
c                    rad(ij,id)=radcm(ij,id)
c                 END DO
c              end do
c
            end do
            end if
c 
         end do
         end if
c 
c     ---------------------------------------------
c     end of formal solutions
c     ---------------------------------------------
c
      end do
      end if
      return
      end
c

C
C
C     ****************************************************************
C
C
C

      SUBROUTINE RTECF1(IJ)
C     =====================
C
C     Solution of the radiative transfer equation with Compton scattering
C     for one frequency (assuming the radiation intensity in i
C     other frequencies is given
C     solution is done for individual angles, and new Eddington factors 
C     are determined
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'MODELQ.FOR'
      INCLUDE 'ALIPAR.FOR'
      INCLUDE 'ITERAT.FOR'
      PARAMETER (SIXTH=UN/6.D0,
     *           THIRD=UN/3.D0,
     *           TWOTHR=TWO/3.D0)
      COMMON/OPTDPT/DT(MDEPTH)
      COMMON/SURFEX/EXTJ(MFREQ),EXTH(MFREQ)
      COMMON/EXTINT/WANGLE,EXTIN(MFREQ)
      COMMON/AUXRTE/
     *          COMA(MDEPTH),COMB(MDEPTH),COMC(MDEPTH),VL(MDEPTH),
     *          COME(MDEPTH),U(MDEPTH),V(MDEPTH),BS(MDEPTH),
     *          AL(MDEPTH),BE(MDEPTH),GA(MDEPTH)
      common/comgfs/gfm(mfreq,mdeptc),gfp(mfreq,mdeptc)
      DIMENSION RI(MDEPTH),RDH(MDEPTH),RDK(MDEPTH),RDN(MDEPTH),
     *          DTAU(MDEPTH),ST0(MDEPTH),RDWN(MMUC),
     *          ali(mdepth)
      DIMENSION AANU(MDEPTH),DDD(MDEPTH),FKK(MDEPTH),ali0(mdepth),
     *          SS0(MDEPTH),
     *          AAA(MDEPTH),BBB(MDEPTH),CCC(MDEPTH),EEE(MDEPTH),
     *          ZZZ(MDEPTH),ALRH(MDEPTH),ALRM(MDEPTH),ALRP(MDEPTH),
     *          ss0c(mdepth)
C
      IF(IJ.EQ.1) THEN
         if(icompt.gt.0.and.icombc.gt.0) then
            IJE=IJEX(IJ)
            DO ID=1,ND
               rad1(id)=rad(nfreq,id)
               fak1(id)=0.333333
               ali1(id)=0.
            END DO
            return
         end if
      END IF
C
      WW=W(IJ)
      IJI=NFREQ-KIJ(IJ)+1
      FR=FREQ(IJ)
      CALL RTECF0(IJ)
c
      do id=1,nd
         rad1(id)=0.
         ali1(id)=0.
         rdh(id)=0.
         rdk(id)=0.
         rdn(id)=0.
         st0(id)=vl(id)+(comb(id)+bs(id))*rad(iji,id)
         ss0(id)=0.
      end do
      rdh1=0.
      rdhd=0.
c
      if(iji.gt.1) then
         do id=1,nd
            st0(id)=st0(id)+coma(id)*rad(iji-1,id)
         end do
      end if
      if(iji.lt.nfreq) then
         do id=1,nd
            st0(id)=st0(id)+comc(id)*rad(iji+1,id)
         end do
      end if
c
      if(idisk.eq.0.or.ifz0.lt.0) then
         FR15=FR*1.D-15
         BNU=BN*FR15*FR15*FR15
         PLAND=BNU/(EXP(HK*FR/TEMP(ND))-UN)*RRDIL
         DPLAN=BNU/(EXP(HK*FR/TEMP(ND-1))-UN)*RRDIL
         IF(TEMPBD.GT.0.) THEN
           PLAND=BNU/(EXP(HK*FR/TEMPBD)-UN)*RRDIL
           DPLAN=BNU/(EXP(HK*FR/TEMPBD)-UN)*RRDIL
         ENDIF
         DPLAN=(PLAND-DPLAN)/DT(ND-1)
      end if
c
      if(icomrt.eq.0) then
c
c     ========================================================
c     Formal angle-dependent solution done by Feautrier scheme
c     ========================================================
c
c     loop over angles points
c
      do i=1,nmu
         do id=1,nd-1
            dtau(id)=dt(id)/amu(i)
         end do
c
c        boundary conditions
c
         rup=0.
         rdown=0.
         rup=extint(ij,i)
         if(idisk.eq.0.or.ifz0.lt.0) rdown=pland+amu(i)*dplan
c
c        solution of the transfer equation
c
         call rtefe2(dtau,st0,rup,rdown,ri)
         ttau=0.
         do id=1,nd
            riid=wtmu(i)*ri(id)
            rad1(id)=rad1(id)+riid
            rdk(id)=rdk(id)+amu(i)*amu(i)*riid
         end do
         rdh1=rdh1+amu(i)*wtmu(i)*ri(1)
         rdhd=rdhd+amu(i)*wtmu(i)*ri(nd)
      end do
      rdh1=rdh1-half*hextrd(ij)
c
c     ----------------------
c     end of the loop over angle points
c
c     ===========================================
c     Formal angle-dependent solution done by DFE
c     ===========================================
c
      else
c
c     loop over angle points
c     ----------------------
c
      do i=1,nmuc
         do id=1,nd-1
            dtau(id)=dt(id)/abs(amuc(i))
         end do
c
c        boundary conditions
c
         rup=0.
         rdown=0.
         if(amuc(i).lt.0.) rup=extint(ij,i)
C
C        diffusion approximation for semi-infinite atmospheres
C
         if(idisk.eq.0.or.ifz0.lt.0) rdown=pland+amuc(i)*dplan
c
c        the case of finite slab - irradiation of the back side
c
         if(amuc(i).gt.0.) rdown=rdwn(nmuc-i+1)
c
c        solution of the transfer equation
c
         call rtesol(dtau,st0,rup,rdown,amuc(i),ri,ali)
         ttau=0.
         do id=1,nd
            riid=ri(id)*half
            rad1(id)=rad1(id)+wtmuc(i)*riid
            ali1(id)=ali1(id)+wtmuc(i)*ali(id)
            rdh(id)=rdh(id)+amuc1(i)*riid
            rdk(id)=rdk(id)+amuc2(i)*riid
            rdn(id)=rdn(id)+amuc3(i)*riid
         end do
         rdwn(i)=ri(nd)
         if(amuc(i).gt.0.) rdh1=rdh1+amuc1(i)*ri(1)*half
         rdhd=rdhd+abs(amuc1(i))*ri(nd)*half
      end do
c
c     ----------------------
c     end of the loop over angle points
c
      end if
c
      do id=1,nd
         fak1(id)=fak(ij,id)
         radk(ij,id)=rdk(id)
         if(icomve .gt. 0) then
            fkk(id)=rdk(id)/rad1(id)
          else
            fkk(id)=fak(ij,id)
         endif
         ss0(id)=0.
      end do
      if(icomve.gt.0) then
         do id=1,nd
            fak(ij,id)=rdk(id)/rad1(id)
            fak1(id)=fak(ij,id)
            fkk(id)=fak(ij,id)
         end do
      end if
      if(rad1(1).gt.0.) then
         flux(ij)=rdh1
         fhd(ij)=rdhd/rad1(nd)
      end if
c
      ah=rdh1
      if(iwinbl.lt.0) ah=ah+half*hextrd(ij)
      aj=rad1(1)
      fh(ij)=ah/aj
C
C     ********************
C
C     Again solution of the transfer equation, now with Eddington
C     FKK and FH determined above, to insure strict consistency of the
C     radiation field and Eddington factors
C
C     Upper boundary condition
C
      U0=0.
      QQ0=0.
      US0=0.
      TAUMIN=ABSO1(1)*DEDM1
      NMU=3
      DO I=1,NMU
         IF(IWINBL.EQ.0.AND.WANGLE.EQ.0.) THEN
C
C           allowance for non-zero optical depth at the first depth point
C
            TAMM=TAUMIN/AMU(I)
            EX=EXP(-TAMM)
            P0=UN-EX
            QQ0=QQ0+P0*AMU(I)*WTMU(I)
            U0=U0+EX*WTMU(I)
            if(tamm.gt.0.) US0=US0+P0/TAMM*WTMU(I)
         END IF
      END DO
      ID=1
      DTP1=DT(ID)
      IF(MOD(ISPLIN,3).EQ.0) THEN
         B=DTP1*HALF
         C=0.
       ELSE
         B=DTP1*THIRD
         C=B*HALF
      END IF
      BQ=UN/(B+QQ0)
      CQ=C*BQ
      BBB(ID)=(FKK(ID)/DTP1+FH(IJ)+B)*BQ+SS0(ID)
      CCC(ID)=(FKK(ID+1)/DTP1)*BQ-CQ*(UN+SS0(ID+1))
      ZZZ(ID)=UN/BBB(ID)
      VLL=ST0(ID)+CQ*ST0(ID+1)
c     IF(IWINBL.LT.0) VLL=VLL+TWO*HEXTRD(IJ)/DTP1
      AANU(ID)=VLL*ZZZ(ID)
      DDD(ID)=CCC(ID)*ZZZ(ID)
      IF(ISPLIN.GT.2) FFF=BBB(ID)/CCC(ID)-UN
C
C     Normal depth point
C
      DO 280 ID=2,ND-1
         DTM1=DTP1
         DTP1=DT(ID)
         DT0=TWO/(DTP1+DTM1)
         ALP=UN/DTM1*DT0
         GAM=UN/DTP1*DT0
         IF(MOD(ISPLIN,3).EQ.0) THEN
            A=0.
            C=0.
          ELSE IF(ISPLIN.EQ.1) THEN
            A=DTM1*DT0*SIXTH
            C=DTP1*DT0*SIXTH
          ELSE
            A=(UN-HALF*DTP1*DTP1*ALP)*SIXTH
            C=(UN-HALF*DTM1*DTM1*GAM)*SIXTH
         END IF
         AAA(ID)=ALP*FKK(ID-1)-A*(UN+SS0(ID-1))
         CCC(ID)=GAM*FKK(ID+1)-C*(UN+SS0(ID+1))
         BBB(ID)=(ALP+GAM)*FKK(ID)+(UN-A-C)*(UN+SS0(ID))
         VLL=A*ST0(ID-1)+C*ST0(ID+1)+(UN-A-C)*ST0(ID)
         AANU(ID)=VLL+AAA(ID)*AANU(ID-1)
         IF(ISPLIN.LE.2) THEN
            ZZZ(ID)=UN/(BBB(ID)-AAA(ID)*DDD(ID-1))
            DDD(ID)=CCC(ID)*ZZZ(ID)
            AANU(ID)=AANU(ID)*ZZZ(ID)
          ELSE
            SUM=-AAA(ID)+BBB(ID)-CCC(ID)
            FFF=(SUM+AAA(ID)*FFF*DDD(ID-1))/CCC(ID)
            DDD(ID)=UN/(UN+FFF)
            AANU(ID)=AANU(ID)*DDD(ID)/CCC(ID)
         ENDIF
  280 CONTINUE
C
C     Lower boundary condition
C
      ID=ND
c
c     stellar atmospheric
c
      IF(IDISK.EQ.0.OR.IFZ0.LT.0) then
      IF(IBC.EQ.0) THEN
         BBB(ID)=FKK(ID)/DTP1+HALF
         AAA(ID)=FKK(ID-1)/DTP1
         VLL=HALF*PLAND+THIRD*DPLAN
       ELSE IF(IBC.LT.4) THEN
         B=UN/DTP1
         A=TWO*B*B
         BBB(ID)=UN+SS0(ID)+B*TWO*FHD(IJ)+A*FKK(ID)
         AAA(ID)=A*FKK(ID-1)
         VLL=ST0(ID)+B*(PLAND+TWOTHR*DPLAN)
       ELSE 
         B=UN/DTP1
         A=TWO*B*B
         BBB(ID)=B+A*FKK(ID)
         AAA(ID)=A*FKK(ID-1)
         VLL=B*(PLAND+TWOTHR*DPLAN)
      END IF
c
c     accretion disk - symmetric boundary
c
      ELSE
         B=TWO/DTP1
         BBB(ID)=FKK(ID)/DTP1*B+UN+SS0(ND)
         AAA(ID)=FKK(ID-1)/DTP1*B
         VLL=ST0(ID)
      END IF
C
      EEE(ND)=AAA(ID)/BBB(ID)
      ZZZ(ID)=UN/(BBB(ID)-AAA(ID)*DDD(ID-1))
      RAD1(ID)=(VLL+AAA(ID)*AANU(ID-1))*ZZZ(ID)
      FAK1(ID)=FKK(ND)
      ALRH(ID)=ZZZ(ID)
      frd=bbb(nd)*rad1(nd)-aaa(nd)*rad1(nd-1)
      frd1=(bbb(nd)-un)*rad1(nd)-aaa(nd)*rad1(nd-1)
C
C     Backsolution
C
      DO 290 ID=ND-1,1,-1
         EEE(ID)=AAA(ID)/(BBB(ID)-CCC(ID)*EEE(ID+1))
         RAD1(ID)=AANU(ID)+DDD(ID)*RAD1(ID+1)
         FAK1(ID)=FKK(ID)
C        write(42,642),ij,id,rad1(id),st0(id),fak1(id)
         ALRH(ID)=ZZZ(ID)/(UN-DDD(ID)*EEE(ID+1))
         ALRM(ID)=0
         ALRP(ID)=0
  290 CONTINUE
c     flux(ij)=fh(ij)*rad1(1)-half*hextrd(ij)
c
C        evaluate approximate Lambda operator
C
C        a) Rybicki-Hummer Lambda^star operator (diagonal)
C           (for JALI = 1)
C
         DO 301 ID=1,ND
            ALIM1(ID)=0.
            ALIP1(ID)=0.
  301    CONTINUE
         IF(JALI.EQ.1) THEN
         DO 310 ID=1,ND 
            ALI1(ID)=ALRH(ID)
  310    CONTINUE
c
         IF(IBC.EQ.0) THEN
            ali1(nd-1)=rad1(nd-1)/st0(nd-1)
            ali1(nd)=rad1(nd)/st0(nd)
         END IF
C
C        for IFALI>5:
C        tridiagonal Rybicki-Hummer operator (off-diagonal terms)
C
         IF(IFALI.GE.6) THEN
         ALIP1(1)=ALRH(2)*DDD(1)
         DO 340 ID=2,ND-1 
            ALIM1(ID)=ALRH(ID-1)*EEE(ID)
            ALIP1(ID)=ALRH(ID+1)*DDD(ID)
  340    CONTINUE
         ALIM1(ND)=ALRH(ND-1)*EEE(ND)
         IF(IBC.EQ.0) THEN
            ALIM1(nd)=0.
            ALIM1(nd-1)=0.
            ALIP1(nd)=0.
            ALIP1(nd-1)=0.
         END IF
         END IF
c
C        b) diagonal Olson-Kunasz Lambda^star operator, 
C           (for JALI = 2)
C
         ELSE IF(JALI.EQ.2) THEN
         DO 360 ID=1,ND-1
            ALI0(ID)=0.
            DO 350 I=1,NMU
               DIV=DT(ID)/AMU(I)
               ALI0(ID)=ALI0(ID)+(UN-EXP(-DIV))/DIV*WTMU(I)
  350       CONTINUE
  360    CONTINUE
         DO 370 ID=2,ND-1
            ALI1(ID)=UN-HALF*(ALI0(ID)+ALI0(ID-1))
  370    CONTINUE
         ALI1(1)=UN-HALF*(ALI0(1)+US0)
         ALI1(ND)=UN-ALI0(ND-1)
         ali1(nd-1)=rad1(nd-1)/st0(nd-1)
         ali1(nd)=rad1(nd)/st0(nd)
       END IF
C          
C        correction of Lambda^star for scattering
C  
         IF(ILMCOR.EQ.1) THEN     
         DO ID=1,ND
             ALI1(ID)=ALI1(ID)*(UN+SS0(ID))
             ALIM1(ID)=ALIM1(ID)*(UN+SS0(ID))
             ALIP1(ID)=ALIP1(ID)*(UN+SS0(ID))
         END DO   
         ELSE IF(ILMCOR.EQ.3) THEN     
         DO ID=1,ND
             ALI1(ID)=ALI1(ID)/(UN+SS0C(ID)*ALI1(ID))
             ALIM1(ID)=ALIM1(ID)/(UN+SS0C(ID)*ALIM1(ID))
             ALIP1(ID)=ALIP1(ID)/(UN+SS0C(ID)*ALIP1(ID))
         END DO   
         END IF
C
c     DO ID=1,ND
c        radcm(iji,id)=rad1(id)
c     END DO
C
C     radiation pressure
C
      if(.not.lskip(1,IJ))
     *PRD0=PRD0+ABSO1(1)*WW*(RAD1(1)*FH(IJ)-HEXTRD(IJ))
      DO ID=1,ND
         if(.not.lskip(ID,IJ))
     *   PRADT(ID)=PRADT(ID)+RAD1(ID)*FAK1(ID)*WW
         PRADA(ID)=PRADA(ID)+RAD1(ID)*FAK1(ID)*WW
      END DO
c
      if(chmax.ge.1.91e-3.and.chmax.le.2.03e-3) then
         tauij=taumin
         do id=1,nd
         if(id.gt.1) tauij=tauij+dt(id-1)
            write(97,697) ij,id,tauij,rad1(id),st0(id)/(un+ss0(id)),
     *      st0(id),un+ss0(id),ali1(id)
         end do
  697    format(2i4,1p6e12.4)
      end if
c
      do id=1,nd
         fak(ij,id)=fak1(id)
      end do
C
      RETURN
      END
C
C
C     ****************************************************************
C
C

      SUBROUTINE RTECMC
C     =================
C
C     Solution of the radiative transfer equation with Compton scattering
C
      INCLUDE 'IMPLIC.FOR'                                                      
      INCLUDE 'BASICS.FOR'
      INCLUDE 'MODELQ.FOR'
      INCLUDE 'ALIPAR.FOR'
      INCLUDE 'ITERAT.FOR'
      COMMON/AUXRTE/
     *          COMA(MDEPTH),COMB(MDEPTH),COMC(MDEPTH),VL(MDEPTH),
     *          COME(MDEPTH),U(MDEPTH),V(MDEPTH),BS(MDEPTH),
     *          AL(MDEPTH),BE(MDEPTH),GA(MDEPTH)
      common/comgfs/gfm(mfreq,mdeptc),gfp(mfreq,mdeptc)
      DIMENSION BB(MDEPTC,MDEPTC),AA(MDEPTC),CC(MDEPTC),
     *          Z(MFREQ,MDEPTC),D(MFREQ,MDEPTC,MDEPTC),
     *          FF(MDEPTC,MDEPTC),ZZ(MDEPTC),
     *          drad(mfreq,mdeptc)
c
      nsti=1
      if(icomst.gt.1) nsti=icomst
      do isti=1,nsti
      DO IJ=1,NFREQ
         IJO=ijorig(ij)
         FR=FREQ(IJO)
         CALL OPACF1(IJO)
         CALL RTECF0(IJO)
         do id=1,nd
            do id1=1,nd
               bb(id,id1)=0.
            end do
         end do
         id=1
         bb(id,id)=be(id)
         bb(id,id+1)=-ga(id)
         do id=2,nd-1
            bb(id,id)=be(id)
            bb(id,id-1)=-al(id)
            bb(id,id+1)=-ga(id)
         end do
         id=nd
         bb(id,id)=be(id)
         bb(id,id-1)=-al(id)
         do id=1,nd
            if(ichcoo.eq.0) then
               bb(id,id)=bb(id,id)+un-comb(id)-bs(id)
             else
               bb(id,id)=bb(id,id)+un-comb(id)
            end if
            aa(id)=coma(id)
            cc(id)=comc(id)
         end do
c
c        linearization matrices for stimulated emission
c
         if(isti.gt.1) then
            do id=1,nd
               vl(id)=vl(id)-bb(id,id)*rad(ij,id)
               bb(id,id)=bb(id,id)-come(id)*rad(ij,id)
               aa(id)=aa(id)+u(id)*rad(ij,id) 
               cc(id)=cc(id)+v(id)*rad(ij,id) 
            end do
            id=1
            vl(id)=vl(id)-bb(id,id+1)*rad(ij,id+1)
            do id=2,nd-1
               vl(id)=vl(id)-bb(id,id-1)*rad(ij,id-1)-
     *                bb(id,id+1)*rad(ij,id+1)
            end do
            id=nd
            vl(id)=vl(id)-bb(id,id-1)*rad(ij,id-1)
            if(ij.gt.1) then
               do id=1,nd
                  vl(id)=vl(id)+aa(id)*rad(ij-1,id)
               end do
            end if
            if(ij.lt.nfreq) then
               do id=1,nd
                  vl(id)=vl(id)+cc(id)*rad(ij+1,id)
               end do
            end if
         end if
c
c        forward sweep of the grand matrix
c
         if(ij.eq.1) then
            call matinv(bb,nd,mdepth)
            do id=1,nd
               sum=0.
               do id1=1,nd
                  d(ij,id,id1)=bb(id,id1)*cc(id1)
                  sum=sum+bb(id,id1)*vl(id1)
               end do
               z(ij,id)=sum
            end do
c
          else 
            do id=1,nd
               do id1=1,nd
                  ff(id,id1)=bb(id,id1)-aa(id)*d(ij-1,id,id1)
               end do
            end do
            call matinv(ff,nd,mdepth)
            do id=1,nd
               do id1=1,nd
                  d(ij,id,id1)=ff(id,id1)*cc(id1)
               end do
            end do
            do id=1,nd
               zz(id)=vl(id)+aa(id)*z(ij-1,id)
            end do
            do id=1,nd
               sum=0.
               do id1=1,nd
                  sum=sum+ff(id,id1)*zz(id1)
               end do
               z(ij,id)=sum
            end do
         end if
      END DO
c
c     ----------------------------------
c     backward sweep of the grand matrix
c     ----------------------------------
c
      if(isti.eq.1) then
         ij=nfreq
         do id=1,nd
            rad(ij,id)=z(ij,id)
         end do
c
         DO IJ=NFREQ-1,1,-1
            do id=1,nd
                sum=0.
                do id1=1,nd
                   sum=sum+d(ij,id,id1)*rad(ij+1,id1)
                end do
                rad(ij,id)=z(ij,id)+sum
            end do
         END DO
      end if
c   
      if(isti.gt.1) then
         ij=nfreq
         do id=1,nd
            drad(ij,id)=z(ij,id)
         end do
c
         DO IJ=NFREQ-1,1,-1
            do id=1,nd
               sum=0.
               do id1=1,nd
                  sum=sum+d(ij,id,id1)*drad(ij+1,id1)
               end do
               drad(ij,id)=z(ij,id)+sum
            end do
         END DO
c
         chmax=0.
         DO IJ=1,NFREQ
            dri=0.
            do id=1,nd
               if(rad(ij,id).gt.0.) dr=drad(ij,id)/rad(ij,id)
               if(abs(dr).gt.chmax) chmax=abs(dr)
               if(abs(dr).gt.dri) dri=abs(dr)
               if(dr.gt.9.) dr=9.
               if(dr.lt.-0.999) dr=-0.999
               rad(ij,id)=rad(ij,id)*(un+dr)
            end do
         END DO
      end if
c
      if(isti.gt.1.and.chmax.lt.1.e-3) go to 100
      end do
c 
  100 continue  
      return
      end
          
C     
C     
C ************************************************************************    
C     
C     

      SUBROUTINE COMPT0(IJ,ID,ab,compa,compb,compc,compe,comps,compd)
C     ===============================================================
C
c     auxiliary quantities for the Compton scattering source function
c
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'MODELQ.FOR'
      INCLUDE 'ALIPAR.FOR'
      INCLUDE 'ITERAT.FOR'
      PARAMETER (XCON=8.0935D-21,YCON=1.68638E-10)
      common/auxcbc/cden1m(mdepth),cden10(mdepth),
     *              cden2m(mdepth),cden20(mdepth)
c
      IJI=NFREQ-KIJ(IJ)+1
      if(iji.eq.1) then
         compa=0.
         compb=0.
         compc=0.
         compd=0.
         compe=0.
         comps=0.
         return
      end if
c
      FR=FREQ(IJ)
      frp=freq(ijorig(iji+1))
      frm=freq(ijorig(iji-1))
      xcomp=fr*xcon
      e2=ycon*temp(id)
      e1=xcomp-3.*e2
c
      del0=two/(dlnfr(iji)+dlnfr(iji-1))
      cder1p(iji)=(un-delj(iji,id))*del0
      cder1m(iji)=-delj(iji-1,id)*del0
      cder10(iji)=-del0*(un-delj(iji-1,id)-delj(iji,id))
      ss0=elec(id)*sige/ab
      if(ichcoo.eq.0) then
         cder10(iji)=-cder1m(iji)-cder1p(iji)
         compa=ss0*(e1*cder1m(iji)+e2*cder2m(iji))
         compb=ss0*(un-xcomp-sigec(ij)/sige+e1*cder10(iji)+
     *         e2*cder20(iji))
         compc=ss0*(e1*cder1p(iji)+e2*cder2p(iji))
       else
         epsnu=(ab-elec(id)*sigec(ij))/ab
         zxxp=xcon*frp+0.5*bnus(iji+1)*rad(iji+1,id)-3.*e2
         zxx0=xcomp+0.5*bnus(iji)*rad(iji,id)-3.*e2
         zxxm=xcon*frm+0.5*bnus(iji-1)*rad(iji-1,id)-3.*e2
         zxxp12=((un-delj(iji,id))*zxxp+delj(iji,id)*zxx0)*del0
         zxxm12=((un-delj(iji-1,id))*zxx0+delj(iji-1,id)*zxxm)*del0
         compa=ss0*(-delj(iji-1,id)*zxxm12+e2*cder2m(iji))
         compc=ss0*((un-delj(iji,id))*zxxp12+e2*cder2p(iji))
         compb=ss0*(delj(iji,id)*zxxp12-(un-delj(iji-1,id))*zxxm12+
     *         e2*cder20(iji)-sigec(ij)/sige)-epsnu+1.
         compe=0.
      end if
      compd=(-3.*cder10(iji)+cder20(iji))*rad(iji,id)
c
      IF(ICOMDE.EQ.0) THEN
         COMPA=0.
         COMPC=0.
         COMPB=0.
      END IF
c
      x0=ss0*bnus(iji)
      if(icomst.eq.0) x0=0.
      if(ichcoo.eq.0) then
         compe=x0*(cder10(iji)-un)
         compu=x0*cder1m(iji)
         compv=x0*cder1p(iji)
         cbs=compe*rad(iji,id)
         compe=cbs
      end if
      comps=compb*rad(iji,id)
      if(iji.gt.1) then 
         if(ichcoo.eq.0) cbs=cbs+compu*rad(iji-1,id)
         comps=comps+compa*rad(iji-1,id)
         compd=compd+(-3.*cder1m(iji)+cder2m(iji))*rad(iji-1,id)
      end if
      if(iji.lt.nfreq) then
         if(ichcoo.eq.0) cbs=cbs+compv*rad(iji+1,id)
         comps=comps+compc*rad(iji+1,id)
         compd=compd+(-3.*cder1p(iji)+cder2p(iji))*rad(iji+1,id)
      end if
      if(ichcoo.eq.0) then
         compb=compb+cbs
         compa=compa+compu*rad(iji,id)
         compc=compc+compv*rad(iji,id)
         comps=comps+cbs*rad(iji,id)
      end if
      compd=compd*ss0*ycon
      IF(ICOMDE.EQ.0) COMPD=0.
c
c     a variant with ICOMPT=2 - no off-diagonal terms in intensity
c
      if(icompt.eq.2) then
         if(iji.gt.1) compb=compb+compa*rad(iji-1,id)
         if(iji.lt.nfreq) compb=compb+compc*rad(iji+1,id)
         compa=0.
         compc=0.
       else if(icompt.eq.3) then
         compa=0.
         compb=0.
         compc=0.
      end if
      return
      end
C
C
C     ****************************************************************
C
C

      SUBROUTINE TAUFR1(IJ)
C     =====================
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'MODELQ.FOR'
      INCLUDE 'ALIPAR.FOR'
      INCLUDE 'ITERAT.FOR'
      COMMON/OPTDPT/DT(MDEPTH)
      dimension ST0(MDEPTH),SS0(MDEPTH),ab0(mdepth),
     *          tau(mdepth),taus(mdepth)
c      PARAMETER (TAUREF = 0.6666666666667)
      PARAMETER (TAUREF = 1.)
      PARAMETER (XCON=8.0935D-21,YCON=1.68638E-10)
C
      FR=FREQ(IJ)
      DO ID=1,ND
         AB0(ID)=ABSO1(ID)
C        put in a floor to avoid division by zero:
         IF(AB0(ID)-SCAT1(ID).GT.1.D-100) THEN
           ST0(ID)=EMIS1(ID)/(AB0(ID)-scat1(id))
          ELSE
           ST0(ID)=0.
         ENDIF
         if(st0(id).eq.0) st0(id)=1.d-20*scat1(id)
         SS0(ID)=-SCAT1(ID)/AB0(ID)
      END DO
C
      id=1
      TAUMIN=ABSO1(1)*DEDM1
      tau(1)=taumin
C     to avoid a negative square root:
      taus(1)=sqrt(3.*ab0(id)*max(ab0(id)-scat1(id),0.d0))*DEDM1
C  
      IREF=1
      IREFs=1
      DO ID=1,ND-1
         DT(ID)=DELDMZ(ID)*(ABSOT(ID+1)+ABSOT(ID))
         tau(id+1)=tau(id)+dt(id)
C        to avoid negative square root:
         eps0=sqrt(ab0(id)*max(ab0(id)-scat1(id),0.d0))
         eps1=sqrt(ab0(id+1)*max(ab0(id+1)-scat1(id+1),0.d0))
         dts=deldm(id)*(eps0*dens1(id)+eps1*dens1(id+1))*sqrt(3.)
         taus(id+1)=taus(id)+dts
         IF(TAU(Id).LE.TAUREF.AND.TAU(Id+1).GT.TAUREF) IREF=Id
         IF(TAUs(Id).LE.TAUREF.AND.TAUs(Id+1).GT.TAUREF) IREFs=Id
      END DO
      if(iref.eq.1.and.tau(nd).le.tauref) iref=nd
      if(irefs.eq.1.and.taus(nd).le.tauref) irefs=nd
C 
      t0=1. 
      iref0=iref
      iref=irefs
      if(irefs.lt.nd) then 
         T0=LOG(TAUs(IREF+1)/TAUs(IREF))
         X0=LOG(TAUs(IREF+1)/TAUREF)/T0
         X1=LOG(TAUREF/TAUs(IREF))/T0
         DMREF=EXP(LOG(DM(IREF))*X0+LOG(DM(IREF+1))*X1)
         TREF=EXP(LOG(TEMP(IREF))*X0+LOG(TEMP(IREF+1))*X1)
         abREF=EXP(LOG(ab0(IREF))*X0+LOG(ab0(IREF+1))*X1)
         scREF=EXP(LOG(scat1(IREF))*X0+LOG(scat1(IREF+1))*X1)
         STREF=EXP(LOG(ST0(IREF))*X0+LOG(ST0(IREF+1))*X1)
         tauef=EXP(LOG(TAU(IREF))*X0+LOG(TAU(IREF+1))*X1)
       else
         x0=1.
         x1=0.
         dmref=dm(nd)
         tref=temp(nd)
         abref=ab0(nd)
         scref=scat1(nd)
         stref=st0(nd)
         tauef=tau(nd)
      end if
      epref=(abref-scref)/abref
CX    add if statement to avoid overflow:
      IF(hk*fr/tref.lt.200.) then
         bref=1.4743e-2*(fr*1.e-15)**3/(exp(hk*fr/tref)-1.)
       ELSE
         bref=0.
      END IF
      taur=tauef*tauef
      if(tauef.gt.taur) taur=tauef
      yref=4.*ycon*tref*taur
      ALM=2.997925E18/FREQ(IJ)
      r1=rad(nfreq-kij(ij)+1,1)
      if(epref.ge.0) rb1=sqrt(epref)*bref
      if(epref.ge.0) rs1=sqrt(epref)*stref
C  
      return
      end
C
C
C     ****************************************************************
C
C
      SUBROUTINE RTECMU
C     =================
C
C     Solution of the radiative transfer equation with Compton scattering
C     for one frequency at a time (assuming the radiation intensity in other 
C     frequencies given), for a number of specific intensities (Gaussian)
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'MODELQ.FOR'
      INCLUDE 'ALIPAR.FOR'
      INCLUDE 'ITERAT.FOR'
      PARAMETER (XCON=8.0935D-21,YCON=1.68638E-10)
      COMMON/OPTDPT/DT(MDEPTH)
      COMMON/AUXRTE/
     *          COMA(MDEPTH),COMB(MDEPTH),COMC(MDEPTH),VL(MDEPTH),
     *          COME(MDEPTH),U(MDEPTH),V(MDEPTH),BS(MDEPTH),
     *          AL(MDEPTH),BE(MDEPTH),GA(MDEPTH)
      DIMENSION RI(MDEPTH),DTAU(MDEPTH),ST0(MDEPTH),ali(mdepth)
      PARAMETER (MW=10, nw=10, zero=0.)
      DIMENSION rmu(MW),b(MW),rintmu(mw),rintpo(mw),rdwn(2*mw),
     *          rmmu(2*MW),wmmu(2*mw)
      DIMENSION RJTOT(MDEPTH),RJNUT(MDEPTH),RDJ1(MDEPTH),
     *          ABSCAD(MDEPTH),ABRAD(MDEPTH),PLTOT(MDEPTH),
     *          retot(mdepth),re1(mdepth),re2(mdepth),
     *          recm(mdepth),recm0(mdepth),scom(mdepth)
c
c     set-up Gaussian angle points
c
      call gauleg(zero,un,rmu,b,nw,mw)
      do i=1,nw
         rmmu(i)=-rmu(nw-i+1)
         rmmu(i+nw)=rmu(i)
         wmmu(i)=b(nw-i+1)
         wmmu(i+nw)=b(i)
      end do
C
      SUMW=0.
      DO ID=1,ND
         RJTOT(ID)=0.
         RJNUT(ID)=0.
         ABSCAD(ID)=0.
         ABPLAD(ID)=0.
         ABRAD(ID)=0.
         pltot(id)=0.
         retot(id)=0.
         re1(id)=0.
         re2(id)=0.
         recm(id)=0.
         recm0(id)=0.
      END DO
C
c     --------------------- loop over frequencies
c
      do 500 ij=1,nfreq
      IJI=NFREQ-KIJ(IJ)+1
      FR=FREQ(IJ)
      xcomp=xcon*fr
      CALL OPACF1(IJ)
      if(icompt.gt.0) then
         CALL RTECF0(IJ)
       else
         do id=1,nd
            if(id.lt.nd) dt(id)=deldmz(id)*(absot(id+1)+absot(id))
            comb(id)=elec(id)*sige/abso1(id)
            vl(id)=emis1(id)/abso1(id)
            coma(id)=0.
            comc(id)=0.
           bs(id)=0.
         end do
      end if
c
      SUMW=SUMW+W(IJ)
      do id=1,nd
         x0=elec(id)*sige/abso1(id)
         vl(id)=emis1(id)/abso1(id)
         st0(id)=vl(id)+(comb(id)+bs(id))*rad(iji,id)
         RDJ1(ID)=0.
         ABSCAD(ID)=ABSCAD(ID)+SCAT1(ID)*W(IJ)
         scom(id)=(comb(id)-x0*(un-two*xcomp)+bs(id))*rad(iji,id)
      end do
      rdh1=0.
      rdhd=0.
c
      if(iji.gt.1) then
         do id=1,nd
            st0(id)=st0(id)+coma(id)*rad(iji-1,id)
            scom(id)=scom(id)+coma(id)*rad(iji-1,id)
         end do
      end if
      if(iji.lt.nfreq) then
         do id=1,nd
            st0(id)=st0(id)+comc(id)*rad(iji+1,id)
            scom(id)=scom(id)+comc(id)*rad(iji+1,id)
         end do
      end if
c
      if(ifz0.lt.0) then
         fr15=fr*1.d-15
         bnu=bn*fr15*fr15*fr15
         pland=0.
         x=hk*fr/temp(nd)
         ex=exp(-x)
         pland=bnu*ex/(un-ex)
         dplan=0.
         x=hk*fr/temp(nd-1)
         ex=exp(-x)
         dplan=bnu*ex/(un-ex)
         dplan=(pland-dplan)/dt(nd-1)
      end if
C
c     --------------------- loop over angles
c
      do i=1,2*nw
         do id=1,nd-1
            dtau(id)=dt(id)/abs(rmmu(i))
         end do
c
c        boundary conditions
c
         rup=extint(ij,i)
C
C        diffusion approximation for semi-infinite atmospheres
C
         if(ifz0.lt.0) rdown=pland+rmmu(i)*dplan
c
c        the case of finite slab - irradiation of the back side
c
         if(rmmu(i).gt.0.) rdown=rdwn(nw-i+1)
c
c        solution of the transfer equation
c
         call rtesol(dtau,st0,rup,rdown,rmmu(i),ri,ali)
c
         if(rmmu(i).gt.0.) then
            if(ri(1).lt.1.e-35) ri(1)=1.e-35
            rintmu(i-nw)=ri(1)
            rintpo(i-nw)=0.
            rdh1=rdh1+rmmu(i)*wmmu(i)*ri(1)*half
         end if
         rdwn(i)=ri(nd)
         rdhd=rdhd+abs(rmmu(i)*wmmu(i))*ri(nd)*half
         DO ID=1,ND
            RDJ1(ID)=RDJ1(ID)+WMMU(I)*RI(ID)*HALF
         END DO
      end do
C
c     --------------------- end of loop over angles
c
      BBN=1.4743E-2*(FR*1.E-15)**3
      DO ID=1,ND
         pla=0.
         x=hk*fr/temp(nd)
         ex=exp(-x)
         pla=bbn*ex/(un-ex)*w(ij)
         RJTOT(ID)=RJTOT(ID)+RDJ1(ID)*W(IJ)
         RJNUT(ID)=RJNUT(ID)+RDJ1(ID)*FREQ(IJ)*W(IJ)
         ABRAD(ID)=ABRAD(ID)+RDJ1(ID)*W(IJ)*(ABSO1(ID)-SCAT1(ID))
         ABPLAD(ID)=ABPLAD(ID)+PLA*(ABSO1(ID)-SCAT1(ID))
         PLTOT(ID)=PLTOT(ID)+PLA
         retot(id)=retot(id)+abso1(id)*(st0(id)-rdj1(id))*w(ij)
         re1(id)=re1(id)+(abso1(id)-scat1(id))*rdj1(id)*w(ij)
         re2(id)=re2(id)+emis1(id)*w(ij)
         recm(id)=recm(id)+(st0(id)-vl(id)-
     *            scat1(id)/abso1(id)*rdj1(id))*w(ij)
         recm0(id)=recm0(id)+scom(id)*w(ij)
      END DO
c
      wll=2.997925e18/fr
      WRITE(14,641) wll,rdh1,(RINTMU(I),RINTPO(I),I=1,NW)
  500 CONTINUE
C
c     --------------------- end of loop over frequencies
c
  641 FORMAT(1H ,f15.3,1pe15.5/(1P5E15.5))
c
      tautot=dm(nd)*elec(nd)*sige/dens(nd)
      DO ID=1,ND
         ABSCAD(ID)=elec(ID)*sige/DENS(ID)
         ABRAD(ID)=ABRAD(ID)/DENS(ID)/RJTOT(ID)
         ABPLAD(ID)=ABPLAD(ID)/DENS(ID)/PLTOT(ID)
         XNU=RJNUT(ID)/RJTOT(ID)
         re1(id)=re1(id)/dens(id)
         re2(id)=re2(id)/dens(id)
         retot(id)=retot(id)/dens(id)
         taurr=dm(id)*abscad(id)
         xl=abplad(id)*(temp(id)/teff)**4
         xr1=0.75*(1./sqrt(3.)+taurr*(1.-0.5*taurr/tautot))
         xr3a=4.*temp(id)*ycon
         xr3b=xnu*xcon
         xr3=xr3a-xr3b
         xr4=abscad(id)*xr3
         xx1=xr1*(abrad(id)-xr4)
         xx2=0.25/dm(nd)
         xr=xx1+xx2
         xtj=sig4p*4.*teff**4*xr1
         XH1=ABPLAD(ID)*PLTOT(ID)
         XH2=ABRAD(ID)*RJTOT(ID)
         XH12=XH1-XH2
         XH3=XR4*RJTOT(ID)
         XH123=XH12+XH3
         XHR=SIG4P*TEFF**4/DM(ND)
      END DO
      RETURN
      END
C
C
C     ****************************************************************
C
C
      SUBROUTINE RTEANG
C     =================
C
C     initialization of the angle quadrature points for the radiative
C     transfer equation
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'MODELQ.FOR'
      INCLUDE 'ALIPAR.FOR'
      PARAMETER (NMU3=3, NMU5=5, ZERO=0.D0)
      COMMON/EXTINT/WANGLE,EXTIN(MFREQ)
      COMMON/SURFEX/EXTJ(MFREQ),EXTH(MFREQ)
      DIMENSION AMU0(MMU),WTMU0(MMU)
C
C     If irradiation is neglected, the angular quadrature is a standard
C     NMU-point Gaussian quadrature
C
      X=WANGLE*HALF
      XJ=0.
      XH=0.
      IF(X.LE.0.) THEN
         call gauleg(zero,un,amu0,wtmu0,nmu,mmu)
         do i=1,nmu
            amu(i)=amu0(i)
            wtmu(i)=wtmu0(i)
            fmu(i)=0.
         end do
       ELSE
C
C     Here, allowance is made for irradiation by central star.
C     First, establish angular integration that takes into account
C     angles with mu < 0; instead of the standard 3-point integration
C     over angles, we have now a more general NMU5-point integration
C
         X0=HALF-X
         X1=HALF+X
         call gauleg(-un,un,amu0,wtmu0,nmu3,mmu)
         DO I=1,NMU3
            AMU(I)=X0*AMU0(I)+X1
            WTMU(I)=X0*WTMU0(I)
            FMU(I)=0.
         END DO
         NMU=NMU5
         i4=nmu3+1
         i5=nmu3+2
         AMU(i4)=X*(UN+0.577350269189626D0)
         AMU(i5)=X*(UN-0.577350269189626D0)
         DO I=NMU3+1,NMU5
            WTMU(I)=X
            FMU(I)=ASIN(SQRT((WANGLE**2-AMU(I)**2)/(UN-AMU(I)**2)))/
     *             3.141592653589793D0
            XJ=XJ+WTMU(I)*FMU(I)
            XH=XH+WTMU(I)*AMU(I)*FMU(I)
         END DO
      END IF
C
      DO IJ=1,NFREQ
         EXTJ(IJ)=XJ*EXTIN(IJ)*HALF
         EXTH(IJ)=XH*EXTIN(IJ)*HALF
      END DO
C
      RETURN
      END
C
C
C
C     ****************************************************************
C
C
      SUBROUTINE PRD(IJ)
c     ==================
c     
c     modification of the line emission coefficient 
c     and the scattering coefficient in the case of PRD
c
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'ATOMIC.FOR'
      INCLUDE 'MODELQ.FOR'
      INCLUDE 'ITERAT.FOR'
      parameter(a21=4.699e8,pi2=6.28318531,gr=2.*4.8e-8)
c
      if(ij.gt.0) then
c     if(ilam.le.1) return
      FR=FREQ(IJ)
      IF(ISPODF.EQ.0) THEN
      IF(IJLIN(IJ).GT.0) THEN
C
C     the "primary" line at the given frequency
C
         ITR=IJLIN(IJ)
         ITRPRD=IPRD(ITR)
         IF(ITRPRD.GT.0) THEN
         DFR=ABS(FREQ(IJ)-FR0(ITR))
         if(ilow(itr).eq.nfirst(ielh)) then
            omeg=dfr*pi2
            gra=a21+gr*popul(nfirst(ielh),id)
            do id=1,nd
               coher(itrprd,id)=a21/
     *             (gra+gami(2,'elec',omeg,temp(id),elec(id)))
            end do
         end if
         DO ID=1,ND
            SG=PRFLIN(ID,IJ)
            IF(DFR/DOPTR(ITRPRD,ID).LE.XPDIV) SG=0.
            SCALIN=SG*ABTRA(ITR,ID)*COHER(ITRPRD,ID)
            SCAT1(ID)=SCAT1(ID)+SCALIN
            scem=sg*emtra(itr,id)*coher(itrprd,id)*xkfb(id)
c           EMIS1(ID)=EMIS1(ID)-SCALIN*RJBAR(ITRPRD,ID)
            EMIS1(ID)=EMIS1(ID)-SCEM
         END DO 
         END IF
      END IF
      IF(NLINES(IJ).GT.0) THEN
C
C     the "overlapping" lines at the given frequency
C
      DO 100 ILINT=1,NLINES(IJ)
         ITR=ITRLIN(ILINT,IJ)
         ITRPRD=IPRD(ITR)
         IF(ITRPRD.EQ.0) GO TO 100
         IJ0=IFR0(ITR)
         DO 60 IJT=IJ0,IFR1(ITR)
            IF(FREQ(IJT).LE.FR) THEN
               IJ0=IJT
               GO TO 70
            END IF
   60    CONTINUE
   70    IJ1=IJ0-1
         A1=(FR-FREQ(IJ0))/(FREQ(IJ1)-FREQ(IJ0))
         A2=UN-A1
         DFR=ABS(FREQ(IJ)-FR0(ITR))
         if(ilow(itr).eq.nfirst(ielh)) then
            omeg=dfr*pi2
            gra=a21+gr*popul(nfirst(ielh),id)
            do id=1,nd
               coher(itrprd,id)=a21/
     *              (gra+gami(2,'elec',omeg,temp(id),elec(id)))
            end do
         end if
         DO ID=1,ND
            SG=A1*PRFLIN(ID,IJ1)+A2*PRFLIN(ID,IJ0)
            IF(DFR/DOPTR(ITRPRD,ID).LE.XPDIV) SG=0.
            SCALIN=SG*ABTRA(ITR,ID)*COHER(ITRPRD,ID)
            scem=sg*emtra(itr,id)*coher(itrprd,id)*xkfb(id)
            SCAT1(ID)=SCAT1(ID)+SCALIN
c           EMIS1(ID)=EMIS1(ID)-SCALIN*RJBAR(ITRPRD,ID)
            EMIS1(ID)=EMIS1(ID)-SCEM
         END DO 
  100 CONTINUE
      END IF
C
C     Opacity sampling option
C
      ELSE
      IF(NLINES(IJ).GT.0) THEN
      DO 300 ILINT=1,NLINES(IJ)
         ITR=ITRLIN(ILINT,IJ)
         ITRPRD=IPRD(ITR)
         IF(ITRPRD.EQ.0) GO TO 300
       KJ=IJ-IFR0(ITR)+KFR0(ITR)
       INDXPA=IABS(INDEXP(ITR))
       IF(INDXPA.NE.3 .AND. INDXPA.NE.4) THEN
         DFR=ABS(FREQ(IJ)-FR0(ITR))
         if(ilow(itr).eq.nfirst(ielh)) then
            omeg=dfr*pi2
            gra=a21+gr*popul(nfirst(ielh),id)
            do id=1,nd
               coher(itrprd,id)=a21/
     *              (gra+gami(2,'elec',omeg,temp(id),elec(id)))
            end do
         end if
       DO ID=1,ND
            SG=PRFLIN(ID,KJ)
            IF(DFR/DOPTR(ITRPRD,ID).LE.XPDIV) SG=0.
            SCALIN=SG*ABTRA(ITR,ID)*COHER(ITRPRD,ID)
            SCAT1(ID)=SCAT1(ID)+SCALIN
c           EMIS1(ID)=EMIS1(ID)-SCALIN*PJBAR(ITRPRD,ID)
            EMIS1(ID)=EMIS1(ID)-SCEM
       END DO
         END IF
  300 CONTINUE
      END IF
      END IF
      RETURN
c


      end if
c
      do itrp=1,ntrprd
         itr=itrtot(itrp)
         aji=osc0(itr)*g(ilow(itr))/g(iup(itr))*7.42163e-22*
     *       fr0(itr)**2
         omeg=0.
         do id=1,nd
            t=temp(id)
            ane=elec(id)
            call dopgam(itr,id,t,dop,agam) 
            doptr(itrp,id)=dop
            coher(itrp,id)=0.99
            if(agam.gt.0.) coher(itrp,id)=aji/(12.5664*dop*agam)
            if(coher(itrp,id).gt.0.999) coher(itrp,id)=0.999
c            
c           special expression for Lyman-alpha
c
            coher(itrp,id)=aji/(aji+9.8e-8*popul(nfirst(ielh),id)+
     *                     0.667*(gami(2,'iont',omeg,t,ane)+
     *                     gami(2,'elec',omeg,t,ane)))
            rjbar(itrp,id)=pjbar(itrp,id)
c            IF(ID.EQ.1.OR.ID.EQ.10.OR.ID.EQ.20.OR.ID.EQ.40)
c     *      WRITE(107,607) ID,ITR,ITRP,AJI,DOPTR(ITRP,ID),
c     *                     COHER(ITRP,ID),RJBAR(ITRP,ID)
c  607 FORMAT(3I5,1P4E12.4)
         end do
      end do
      return
      END

C
C
C
C     ****************************************************************
C
C
      SUBROUTINE PRDINI
c     =================
c     
c     initialization of PRD
c
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'ATOMIC.FOR'
      INCLUDE 'MODELQ.FOR'
c
      ntrprd=0
      do itr=1,ntrans
         iprd(itr)=0
         if(ifprd.gt.0.and.line(itr).and.indexp(itr).ne.0) then
         ii=ilow(itr)
         jj=iup(itr)
         iat=iatm(ii)
c
c        select Lyman alpha for PRD
c
         if(iat.eq.iath.and.ii.eq.nfirst(ielh).and.
     *      fr0(itr).lt.2.5e15) then
            ntrprd=ntrprd+1
            iprd(itr)=ntrprd
            itrtot(ntrprd)=itr
         end if
c
c        select Mg I resonance line for PRD
c
         if(numat(iat).eq.12.and.iz(iel(ii)).eq.1.and.
     *      ii.eq.nfirst(iel(ii)).and.fr0(itr).lt.1.06e15) then
            ntrprd=ntrprd+1
            iprd(itr)=ntrprd
            itrtot(ntrprd)=itr
         end if
c
c        select Mg II resonance lines for PRD
c
         if(numat(iat).eq.12.and.iz(iel(ii)).eq.2.and.
     *      ii.eq.nfirst(iel(ii)).and.fr0(itr).lt.1.08e15) then
            ntrprd=ntrprd+1
            iprd(itr)=ntrprd
            itrtot(ntrprd)=itr
         end if
         end if
c         if(iprd(itr).gt.0) write(107,607) 
c     * itr,iprd(itr),ii,jj,iat,numat(iat),iel(ii),iz(iel(ii)),
c     * fr0(itr)
c 607 format(8i5,1pe13.5)
      end do
c
      do itrp=1,ntrprd
         do id=1,nd
            pjbar(itrp,id)=0.
         end do
      end do 
      return
      end
C
C
C     ****************************************************************
C
C
 
      function gami(j,aper,omeg,t,ane)
c     ================================
c
c     function i(j) defined by eqs. (4.5)-(4.9) of 
c            cooper, ballagh, and hubeny (1989), ap.j. 344, 949.
c     j    = principal quantum number
c     aper = either 'iont', or 'elec', whether one calculates the
c            ion or electron contribution
c     omeg = delta omega (circular frequency)
c     t    = temperature
c     ane  = electron density (assumed equal to proton density)
c
      INCLUDE 'IMPLIC.FOR'
      character*4 aper
      dimension xx(3)
      data xx/0., 50.6205, 68.6112/
c
      if(omeg.gt.0.) then
         gami=xx(j)*ane/sqrt(omeg)
         return
      end if
c
      x=j*j
      omegp=5.64e4*sqrt(ane)
      amu=1.
      if(aper.eq.'iont') then
         amu=30.2
         omegp=omegp/42.85
      end if
      omegc=1.7455e11*t/amu/amu/j
      corr=0.27-log(8.356e-13*x*amu*amu*ane/t/t)
      gami=3.885e-5*amu*x*ane/sqrt(t)*corr
      if(omeg.lt.omegp) return
      gamp=gami
      gam0=22.58*x**0.75*ane
      gamc=gam0/sqrt(omegc)
      if(omeg.lt.omegc) then
         gami=log(omeg/omegp)/log(omegc/omegp)*log(gamc/gamp)+
     *        log(gamp)
         gami=exp(gami)
       else
         gami=gam0/sqrt(omeg)
      end if
      return
      end
c
C
C
C ********************************************************************
C
C
      SUBROUTINE INPDIS
C     =================
C
C     driver for input specific for disks
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'ATOMIC.FOR'
      INCLUDE 'MODELQ.FOR'
      INCLUDE 'ITERAT.FOR'
      INCLUDE 'ODFPAR.FOR'
      INCLUDE 'ALIPAR.FOR'
      PARAMETER (VELC=2.997925E10, 
     *           pi4=12.5663706d0)
      PARAMETER (GRCON = 6.668D-8)
      common/relcor/arh,brh,crh,drh
C
C ----------------------
C Basic input parameters
C ----------------------
C
C     The user may choose one of the two
C     following possibilities to input the basic physical parameters:
C
C     XMSTAR     - M(star), either in M(Sun), or in grams;
C     XMDOT      - M(dot), either in M(Sun)/year; or in g/s
C     RSTAR      - R(star), either in R(Sun), or in cm
C     RELDST     - R/R(star)
C
C     or to directly give parameters with which the program works, ie
C
C     TEFF       - effective temperature
C     QGRAV      - coefficient in the hydrostatic equilibrium equation
C                     QGRAV=G*M(star)/R**3
C     DMTOT      - total column mass at the midplane
C
      WRITE(6,660)
      IF(FRACTV.LT.0.) THEN
         AMUV0=DMVISC**(ZETA0+UN)   
         FRACTV=UN/(UN+(ZETA0+UN)/(ZETA1+UN)*AMUV0/(UN-AMUV0))  
      END IF
      IF(DMVISC.LT.0.) DMVISC=(UN/(UN+(ZETA0+UN)/(ZETA1+UN)*
     *                 FRACTV/(UN-FRACTV)))**(UN/(ZETA0+UN))
      alpha0=alphav
      WRITE(6,600) XMSTAR,XMDOT,RSTAR,RELDST,ALPHAV
C
c     if XMSTAR<0, turn on general relativistic corrections
c     RSTAR now has the input meaning of dimensionless angular
c                momentum of the Kerr black hole, but we will call this
c                value AA and RSTAR will take on the meaning of 1 radius-
c                equivalent of mass of the black hole
c           RELDST now is expressed in multiples of 1*G*XMSTAR/c^2
c                (note that for Schwarzschild black hole, the
c                horizon is at RELDST=2 and the smallest radius of
c                stable circular orbit is at RELDST=6)
C
      IF(XMSTAR.NE.0.) THEN
         IF(XMSTAR.LT.0) THEN
            AA=RSTAR
            RSTAR=-XMSTAR*1.989D33*GRCON/VELC/VELC
          ELSE
            AA=0.
         END IF
         IF(abs(XMSTAR).GT.1.D16) XMSTAR=XMSTAR/1.989D33
         IF(XMDOT.GT.1.D3) XMDOT=XMDOT/6.3029D25
         IF(RSTAR.GT.1.D3) RSTAR=RSTAR/6.9598D10
         R=RSTAR*ABS(RELDST)
         QGRAV=5.9D0*GRCON*abs(XMSTAR)/R**3
         OMEG32=SQRT(QGRAV)*1.5
c
c        apply general relativistic corrections to
c        QGRAV and TEFF;  keep MSTAR<0 for future use
c
         RR0=RELDST
         CALL GRCOR(AA,RR0,XMSTAR,QCOR,TCOR,ARH,BRH,CRH,DRH)
         TEFF0=(1.79049311D-1*QGRAV*3.34379D24*XMDOT/SIG4P)**0.25
         TEFF=TEFF0*TCOR
         QGRAV=QGRAV*QCOR
         OMEG32=OMEG32*ARH/BRH
         rr0=abs(rr0)
         xmas9=abs(Xmstar)*1.d-9
         xmdt=Xmdot/xmas9/2.22
         XMD=XMDOT*6.3029D25
         alpav=abs(alphav)
c
c        following is the old procedure
c
         if(alphav.le.0.) then
C        -------------------------
            chih=0.39
            if(reynum.le.0.) then
               reynum=(rr0/xmdt)**2/alpav*arh*crh/drh/drh
c              REYNUM=(R/XMDOT*1.10422E-15*12.5663*VELC/CHIH)**2*
c    *                2./ALPAV*ARH*CRH/DRH/DRH 
             else
               alpav=(rr0/xmdt)**2/reynum*arh*crh/drh/drh
c              ALPAV=(R/XMDOT*1.10422E-15*12.5663*VELC/CHIH)**2*
c    *               2./REYNUM*ARH*CRH/DRH/DRH 
            endif
            VISC=1.176565D22*SQRT(GRCON*abs(XMSTAR)*R)/REYNUM
            DMTOT=3.34379D24*XMDOT/VISC*BRH*DRH/ARH/ARH
C
            if(alphav.lt.0.) then
c
C ****************************************************************
C Compute the Keplerian rotation frequency (omega=c/r_g/x^1.5):
C
               RE=ABS(RR0)
               OMEGA=VELC/RSTAR/6.9698D10/RE**1.5D0
C
C Compute Relativistic factors, using Krolik (1998) notation:
C
               RELT=DRH/ARH
               RELR=DRH/BRH
               RL2=RE*(1.D0-2.D0*AA/RE**1.5D0+AA*AA/RE/RE)**2/BRH
               EINF=(1.D0-2.D0/RE+AA/RE**1.5D0)/SQRT(BRH)
               RELZ=(RL2-AA*AA*(EINF-1.D0))/RE
C
C Compute the surface mass density (assuming pure electron scattering, 
C pure hydrogen composition, and that T_rphi =:w \alpha P_tot):
C
               DMTOT=0.5D0*SIGMAR(ALPAV,XMD,TEFF,OMEGA,RELR,RELT,RELZ)
            end if
c           ----------------------------
c
c           new procedure
c
          else
            call column
         end if
c        
         EDISC=SIG4P*TEFF**4/DMTOT
         WBARM=XMDOT*6.3029D25/6./3.1415926*BRH*DRH/(ARH*ARH)
         reynum=dmtot/wbarm*sqrt(xmstar*r)*3.03818e18
         WRITE(6,601) TEFF,QGRAV,DMTOT,
     *                ZETA0,ZETA1,FRACTV,DMVISC,TSTAR
         write(6,321) tcor,qcor,arh,brh,crh,drh,
     *                   reynum,alpav,
     *                   dmtot,edisc,wbarm
  321 FORMAT(
     * ' tcor      =',1PD10.3/
     * ' qcor      =',D10.3/
     * ' A(RH)     =',D10.3/
     * ' B(RB)     =',D10.3/
     * ' C(RH)     =',D10.3/
     * ' D(RH)     =',D10.3/
     * ' Re        =',D10.3/
     * ' alpha     =',D10.3//
     * ' DMTOT     =',D10.3/
     * ' EDISC     =',D10.3/
     * ' WBARM     =',D10.3//)
C
       ELSE
         TEFF=XMDOT
         QGRAV=RSTAR
         DMTOT=RELDST
         EDISC=SIG4P*TEFF**4/DMTOT
         OMEG32=SQRT(QGRAV)*1.5
         WRITE(6,601) TEFF,QGRAV,DMTOT,
     *                ZETA0,ZETA1,FRACTV,DMVISC,TSTAR
      END IF
C
C     set up the maximum frequency
C
c     if(idgrey.le.2) then
c     IF(FRCMAX.EQ.0.) FRCMAX=2.83e11*(dmtot*0.39)**0.25*teff
c     end if
      IF(FRLMAX.EQ.0.) FRLMAX=1.D11*CNU1*TEFF
C
  660 FORMAT(1H1,'***************************************'//
     *           ' M O D E L   O F   A   D I S K   R I N G'//
     *           ' ***************************************'//)
  600 FORMAT(
     * ' M(STAR)   =',1PD10.3/
     * ' M(DOT)    =',D10.3/
     * ' R(STAR)   =',D10.3/
     * ' R/R(STAR) =',D10.3/
     * ' ALPHA     =',D10.3//)
  601 FORMAT(
     * ' TEFF      =',F10.0/
     * ' QGRAV     =',1PD10.3/
     * ' DMTOT     =',D10.3/
     * ' ZETA0     =',D10.3/
     * ' ZETA1     =',D10.3/
     * ' FRACTV    =',D10.3/
     * ' DMVISC    =',D10.3/
     * ' TSTAR     =',0PF10.0//)
      RETURN
      END
C
C
C     ****************************************************************
C
C
      subroutine column
c     =================
c
c     approximate determination of the total disk column
c     mass, DMTOT
c
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'MODELQ.FOR'
      common/relcor/arh,brh,crh,drh
c
      parameter (xmdsun  = 6.3029e25,
     *           xmsun   = 1.989e33,
     *           rsun    = 6.9598e10,
     *           grcon   = 6.668e-8,
     *           velc    = 2.997925e10,
     *           rgas    = 1.3e8,
     *           xkram0  = 7.e25,
     *           chiel   = 0.39,
     *           pi      = 3.14159265e0,
     *           pi4     = 4.*pi)
c
      alpha= abs(alphav)
      corkra=10.
      xkap0=xkram0/corkra
      r = rstar*abs(reldst)
      ga=xmdot*xmdsun/pi4*sqrt(5.9*grcon*abs(xmstar)/r**3)*drh/arh
c
      be=0.8*rgas*xkap0**0.125*(two*qgrav/pi/rgas)**0.0625*sqrt(teff)
      al=(sig4p*pi4*teff**4*chiel/velc)**2/(3.*qgrav)
c
      dm00=(ga/alpha/be)**0.8
      write(6,640) ga,al,be,dm00
  640 format(/' new procedure to determine M_tot'/
     *      ' gam, al, be, dm0 ',1p4e11.3/
     *      ' iter  M  delta(M)/M   p, jac'/) 
         itdm=0
   10    itdm=itdm+1
         p0=alpha*dm00*(al+be*dm00**0.25)-ga
         ppr=alpha*(al+1.25*be*dm00**0.25)
         ddm0=-p0/ppr
         write(6,641) itdm,dm00,ddm0/dm00,p0,ppr
  641    format(i4,1p4e11.3)
         dm00=dm00+ddm0
         if(abs(ddm0/dm00).gt.1.e-2.and.itdm.lt.20) go to 10
      dmtot=dm00
      visc=3.34379D24*XMDOT/dmtot*BRH*DRH/ARH/ARH
c
      return
      end   
C
C
C     ****************************************************************
C
        SUBROUTINE GRCOR(AA,RR,XMSTAR,QCOR,TCOR,ARH,BRH,CRH,DRH)
C      =======================================================
C
C       Procedure for computing general-relativistic correction
C       factors to gravitational factor (QGRAV) and effective 
C       temperature (TEFF)
C       Also calculates all frour quantities in the Riffer-Herlod (RH)
C       notation - ARH, BRH, CRH, DRH
C
C       Input:
C             AA   - angular momentum (0.98 maximum)
C             RR   - R/R_g = r/(GM/c^2)
C       Outout:
C             QCOR - g-correction  = C/B   in RH notation
C             TCOR - T-correction  = (D/B)^(1/4)   in RH notation
C             ARH  - A  in RH notation
C             BRH  - B  in RH notation
C             CRH  - C  in RH notation
C             DRH  - D  in RH notation
C
      INCLUDE 'IMPLIC.FOR'
      PARAMETER (THIRD=1.D0/3.D0, PI3=1.0471976)
C
C  ----------------
C  Imput parameters
C  ----------------
C
C       AA      - specific angular momentum/mass
C                  of the Kerr black hole
C       RR      - distance/mass of the Kerr black hole
C
C  -----------------------------------
C  Classical case  - no GR  corrections
C  ------------------------------------
C
        if(Xmstar.gt.0.) then
           arh=1.
           brh=1.
           crh=1.
           drh=1.-sqrt(1./rr)
           qcor=1.
           tcor=drh**0.25
           return
        end if
c
C  ---------------------------------
C  Set correcion factors A through G  (see Novikov & Thorne,'73, eq.5.4.1a-g)
C  ---------------------------------
C
        rror=rr
        rr=abs(rr)
      AA2=AA*AA
      RR1=1/RR
      RR12=SQRT(RR1)
      RR2=RR1*RR1
      A2R2=AA2*RR2
      A4R4=A2R2*A2R2
      A2R3=AA2*RR2*RR1
      AR32=SQRT(A2R3)
C
      A = 1 +   A2R2 + 2*A2R3
      B = 1 +   AR32
      C = 1 - 3*RR1  + 2*AR32
      D = 1 - 2*RR1  +   A2R2
      E = 1 + 4*A2R2 - 4*A2R3 + 3*A4R4
C
C  -------------------------------
C  Set correction factor for QGRAV  (see Novikov & Thorne,'73, eq.5.7.2)
C  -------------------------------
C
      if(rror.lt.0) QCOR = B*B*D*E/(A*A*C)
c
c  correction - after Riffert and Harold 
c
        if(rror.gt.0) QCOR = (1. - 4.*AR32 + 3.*A2R2)/C
C
C  -----------------------
C  Set correction factor Q  (see Page & Thorne,'73, eq.35)
C  -----------------------
C
C      Minimum radius for last stable circular orbit per unit mass, X0
C
      Z1 = 1 + (1-AA2)**THIRD * ((1+AA)**THIRD + (1-AA)**THIRD)
      Z2 = SQRT(3*AA2 + Z1*Z1)
      X0 = SQRT(3 + Z2 - SQRT((3-Z1)*(3+Z1+2*Z2)))
C
C      Roots of x^3 - 3x + 2a = 0
C
      CA3 = THIRD * ACOS(AA)
      X1 = 2*COS(CA3-PI3)
      X2 = 2*COS(CA3+PI3)
      X3 = -2*COS(CA3)
C
C       FB = '[]' term in eq. (35) of Page&Thorne '73
C
      X = SQRT(RR)
      C1 = 3*(X1-AA)*(X1-AA)/(X1*(X1-X2)*(X1-X3))
      C2 = 3*(X2-AA)*(X2-AA)/(X2*(X2-X1)*(X2-X3))
      C3 = 3*(X3-AA)*(X3-AA)/(X3*(X3-X1)*(X3-X2))
      AL0 = 1.5*AA*log(X/X0)
      AL1 = log((X-X1)/(X0-X1))
      AL2 = log((X-X2)/(X0-X2))
      AL3 = log((X-X3)/(X0-X3))
      FB = (X-X0 - AL0 - C1*AL1 - C2*AL2 - C3*AL3)
      Q = FB*(1+AR32)*RR12/SQRT(1-3*RR1+2*AR32)

C  ------------------------------
C  Set correction factor for TEFF  (see Novikov & Thorne,'73, eq.5.5.14b)
C  ------------------------------
C
      TCOR = (Q/B/SQRT(C))**0.25
C
C  ------------------------------
C  RH quantities
C  ------------------------------
C 
        ARH = D
        BRH = C
        CRH = 1. - 4.*AR32 + 3.*A2R2
        DRH = Q/B*SQRT(C)
C 
      RETURN
      END
C
C
C     ****************************************************************
C
C
      SUBROUTINE DMDER
C     ================
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'ATOMIC.FOR'
      INCLUDE 'MODELQ.FOR'
      COMMON/DEPTDR/DDM(MDEPTH),DDP(MDEPTH),DD0(MDEPTH),
     *              DDMIN(MDEPTH),DDPLU(MDEPTH),DDA(MDEPTH),
     *              DDC(MDEPTH),DDB(MDEPTH)
C
      DO ID=2,ND-1
         DDM(ID)=DM(ID)-DM(ID-1)
         DDP(ID)=DM(ID+1)-DM(ID)
         DD0(ID)=DM(ID+1)-DM(ID-1)
         DDMIN(ID)=DDP(ID)/DD0(ID)
         DDPLU(ID)=DDM(ID)/DD0(ID)
         DDA(ID)=DDMIN(ID)/DDM(ID)
         DDC(ID)=DDPLU(ID)/DDP(ID)
      END DO
C
      DDM(1)=0.
      DDM(ND)=DM(ND)-DM(ND-1)
      DDP(1)=DM(2)-DM(1)
      DDP(ND)=0.
      DDMIN(1)=0.
      DDMIN(ND)=1.
      DDPLU(1)=1.
      DDPLU(ND)=0.
      DDA(1)=0.
      DDA(ND)=UN/DDM(ND)
      DDC(1)=UN/DDP(1)
      DDC(ND)=0.
      DO ID=1,ND
         DDB(ID)=DDA(ID)-DDC(ID)
      END DO
C
      RETURN
      END
C
C
C     ****************************************************************
C
C
      FUNCTION SIGMAR(ALPHA,XMDT,TEF,OMEGA,RELR,RELT,RELZ)
C     =====================================================
c
C--------------------------------------------------------------------
C The following function takes as inpute various disk parameters computed
C at a certain radius in cgs units, and outputs the surface mass density, 
C assuming that the opacity is electron-scattering dominated (or that kappa 
C is independent of the density).  The equations were derived assuming a 1-zone
C model (i.e. rho, T_g, mu are constant with height), assuming t_r,phi =
C -alpha P, and that the dissipation is constant per unit optical depth. 
C See chapter 7 of Krolik for information on the notation used here 
C (especially relativistic factors).
C     
C--------------------------------------------------------------------
C Uses Numerical Recipes subroutine LAGUER
C--------------------------------------------------------------------
C ALPHA - viscosity parameter
C XMDT  - accretion rate (in g/s)
C TEF   - temperature in Kelvins
C OMEGA - Keplerian frequency (in Hz)
C RELR,RELT, RELZ  - relativistic factors 
C KAPPA - opacity (cm^2/g)
C MU    - mean atomic mass (g) = rho/N
C--------------------------------------------------------------------
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      REAL*8 KAPPA,MU
      COMPLEX*16 COEFF(11),XGUESS
C
C We should check that the physical constants used here agree with those in
C disk195g:
C
      parameter (ZERO=0.D0,ONE=1.D0,TRES=3.D0,
     *           FOUR=4.D0,THIRD=ONE/TRES,FOURTH=0.25D0,EPS=1.D-5)
      parameter (C=2.9979D10,SIGMAB=5.6703D-5)
      parameter (BK=1.3807D-16)
      PI=ACOS(-ONE)
C
C We'll assume fully ionized, pure hydrogen:
C
      KAPPA=0.39D0
      MU=0.5D0*1.6726D-24
      FAC1=RELZ*(HALF*TRES*C*OMEGA/ALPHA/KAPPA/SIGMAB/TEF**4)**2
      FAC2=(HALF*KAPPA)**FOURTH*BK*TEF/MU
      FAC3=XMDT*OMEGA*RELT/PI
C
C Coefficients of the equation for x^4=Sigma:
C
      COEFF(1)=DCMPLX(FAC1*(HALF*FAC3)**2,ZERO)
      COEFF(2)=ZERO
      COEFF(3)=ZERO
      COEFF(4)=ZERO
      COEFF(5)=DCMPLX(-(TRES*FAC3)/(8.D0*ALPHA),ZERO)
      COEFF(6)=DCMPLX(-FAC1*FAC3*ALPHA*FAC2,ZERO)
      COEFF(7)=ZERO
      COEFF(8)=ZERO
      COEFF(9)=ZERO
      COEFF(10)=DCMPLX(FOURTH*FAC2,ZERO)
      COEFF(11)=DCMPLX(FAC1*(ALPHA*FAC2)**2,ZERO)
C
C At small radii, we'll approximate P_rad >> P_gas
C First, compute sigma assuming radiation pressure dominates:
C
      SIGRAD=FOUR*OMEGA*C*C*RELT*RELZ/ALPHA/KAPPA**2/SIGMAB/TEF**4/RELR
C
C Next, compute sigma assuming gas pressure dominates:
C
      SIGGAS=((MU*XMDT*OMEGA*RELT/PI/ALPHA/BK/TEF)**4/8./KAPPA)**0.2D0
C
C Use a starting guess which has the correct value for P_gas >> P_rad
C or P_rad >> P_gas.
C
      XGUESS=DCMPLX(ONE/(ONE/SIGRAD+ONE/SIGGAS)**FOURTH,ZERO)
C
C Look for root of the 10th order equation for x:
C
      CALL LAGUER(COEFF,10,XGUESS,ITS)
C
C Make sure that we haven't landed a wrong root:
C
      IF(ABS(DIMAG(XGUESS)).LT.EPS.AND.DBLE(XGUESS).GT.ZERO) THEN
        SIGMAR=DBLE(XGUESS)**4
      ELSE
        SIGMAR=ONE/(ONE/SIGRAD+ONE/SIGGAS)
        WRITE(6,*) 'Surface density approximated'
      ENDIF
      WRITE(6,2000) TEF,SIGRAD,SIGGAS,SIGMAR
      RETURN
 2000 FORMAT(20(2x,1pe12.5))
      END
C
C
C     ****************************************************************
C
C
      SUBROUTINE laguer(a,m,x,its)
C     ============================
C
C     Routine from Numerical Recipees
C
      INCLUDE 'IMPLIC.FOR'
      COMPLEX*16 a(m+1),x
      PARAMETER (EPSS=2.e-7,MR=8,MT=10,MAXIT=MT*MR)
      REAL frac(MR)
      COMPLEX*16 dx,x1,b,d,f,g,h,sq,gp,gm,g2
      SAVE frac
      DATA frac /.5,.25,.75,.13,.38,.62,.88,1./
C
      do 12 iter=1,MAXIT
        its=iter
        b=a(m+1)
        err=abs(b)
        d=dcmplx(0.d0,0.d0)
        f=dcmplx(0.d0,0.d0)
        abx=abs(x)
        do 11 j=m,1,-1
          f=x*f+d
          d=x*d+b
          b=x*b+a(j)
          err=abs(b)+abx*err
11      continue
        err=EPSS*err
        if(abs(b).le.err) then
          return
        else
          g=d/b
          g2=g*g
          h=g2-2.d0*f/b
          sq=sqrt(dble(m-1)*(dble(m)*h-g2))
          gp=g+sq
          gm=g-sq
          abp=abs(gp)
          abm=abs(gm)
          if(abp.lt.abm) gp=gm
          if (max(abp,abm).gt.0.D0) then
            dx=dble(m)/gp
          else
            dx=exp(dcmplx(log(1.d0+abx),dble(iter)))
          endif
        endif
        x1=x-dx
        if(x.eq.x1)return
        if (mod(iter,MT).ne.0) then
          x=x1
        else
          x=x-dx*frac(iter/MT)
        endif
12    continue
c
      write(6,601) x,x1
  601 format(' too many iterations in laguer, x,x1 ',1p2e9.1)
      return
      END
C
C
C     ****************************************************************
C
C
      SUBROUTINE LTEGRD
C     =================
C
C     Driving procedure for computing the initial LTE-grey disk model
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'MODELQ.FOR'
C
      PARAMETER  (ERRT=1.D-3, THIRD=UN/3.D0, FOUR=4.D0)
      DIMENSION TAU0(MDEPTH),TEMP0(MDEPTH),ELEC0(MDEPTH),
     *       DENS0(MDEPTH),ZD0(MDEPTH),DM0(MDEPTH)
      COMMON/PRSAUX/VSND2(MDEPTH),HG1,HR1,RR1
      COMMON/FACTRS/GAMJ(MDEPTH),GAMH,FAK0
      COMMON/TOTJHK/TOTJ(MDEPTH),TOTH(MDEPTH),TOTK(MDEPTH),
     *              RDOPAC(MDEPTH),FLOPAC(MDEPTH)
      COMMON/FLXAUX/T4,PGAS,PRAD,PGM,PRADM,ITGMAX,ITGMX0
      COMMON/CUBCON/A,B,DEL,GRDADB,DELMDE,RHO,FLXTOT,GRAVD
C
C ----------------
C Input parameters
C ----------------
C
C     NDEPTH  -  number of depth point for evaluating LTE-grey model
C             <  0 - the program computes an isothermal structure;
C                    the temperature is specified by an extra record
C                    (see below)
C     DM1     -  mass at the first depth point
C     ABROS0  -  initial estimate of the Rosseland opacity (per gram)
C                at the first depth point
C     ABPLA0  -  initial estimate of the Planck mean opacity (per
C                gram) at the first depth point
C     DION0   -  initial estimate of the degree of ionization at the
C                first depth point (=1 for completely ionized; =1/2 for
C                completely neutral)
C     ITGMAX  -  maximum number of global iterations of the procedure
C                for determining the Eddington and opacity factors
C        -  number of iterations in recalculating new depth scale
C                in order to get a better coverage of optical depths
C             >  0 - subroutine NEWDM (blind addition of points)
C             <  0 - subroutine NEWDMT - new depths determined as the
C                    equal-segment endpoints of the curve y=tauros(m)
C
C --------------------------------------------------------------------
C
C     IDEPTH  -  mode of determining the mass-depth scale to be used
C                in linearization
C             =  0 - depth scale DM is set up by the program
C             =  1 - interpolation to the depth scale DM (in g*cm**-2),
C                    which has been read in START
C             =  2 - DM is evaluated as mass corresponding to Rosseland
C                    optical depths which are equidistantly spaced in
C                    logarithms between the first point TAU1 and the
C                    last-but-one point TAU2=0.99*TAUMAX; the last point
C                    is TAUMAX, where TAUMAX is the Ross. optical depth
C                    corresponding to the last depth point (midplane),
C                    set up by the program
C     NCONIT       - number of internal iterations for calculating the
C                    LTE-gray model with convection
C             =  0 - if HMIX0>0, program sets NCONIT=10.
C     IPRING       - switch for determining an amount of output from
C                    the calculation of LTE-gray model
C             =  0 - only final LTE-gray model is printed
C             =  1 - more tables are printed
C             =  2 - complete output
C     IHM     >  0 - negative hydrogen ion considered in particle and
C                    charge conservation in ELDENS
C     IH2     >  0 - hydrogen molecule considered in particle
C                    conservation in ELDENS
C     IH2P    >  0 - ionized hydrogen molecule considered in particle
C                    and charge conservation in ELDENS
C
C --------------------------------------------------------------------
C
      IF(NDGREY.EQ.0) THEN
         NDEPTH=ND
       ELSE
         NDEPTH=NDGREY
      END IF
      IF(NDEPTH.GT.MDEPTH) call quit(' NDEPTH too large in LTEGR',
     *                     ndept,mdepth)
      IDEPTH=IDGREY
      ITGMAX=ITGMX0
      IF(HMIX0.GT.0..AND.NCONIT.EQ.0) NCONIT=10
      if(dion0.lt.0) then
         abpmin=-dion0
         dion0=1.
      end if
      T4=TEFF**4
      TOTF=SIG4P*T4
      ABFL0=SIGE/WMM(1)
      if(idmfix.eq.1) then
         t0=teff
         DMTOT=TOTF/EDISC
       else
         call greyd
         t0=temp(1)
         EDISC=TOTF/DMTOT
      end if
c
      ZND=0.
      VSND20=2.76D-16*T0/WMM(1)*DION0+VTB*VTB
      HSCALG=SQRT(TWO*VSND20/QGRAV)
      HSCALR=4.19168946D-10*TOTF*ABFL0/QGRAV
      R=HSCALR/HSCALG
      WRITE(6,615) HSCALG,HSCALR,R
  615 FORMAT(/' GAS PRESSURE SCALE HEIGHT  = ',1PD10.3/
     *        ' RAD.PRESSURE SCALE HEIGHT  = ',1PD10.3/
     *        ' RATIO                      = ',1PD10.3/)
      GAMH=UN
      FAK0=THIRD
      ANEREL=(DION0-HALF)/DION0
      IF(ANEREL.LT.ERRT) ANEREL=ERRT
      IF(NDEPTH.EQ.0) NDEPTH=ND
      LCHC0=LCHC
      LCHC=.TRUE.
      ND0=ND
      ND=NDEPTH
      DO 10 ID=1,ND0
   10    DM0(ID)=DM(ID)
C
C     mass-depth scale
C     Initial estimate of the density, geometrical distance z, and
C     pressure
C
      CALL ZMRHO(R,HSCALG)
C
      if(ipring.eq.2) then
      xdm=dm(1)
      DO ID=1,ND
         if(id.gt.1) xdm=xdm-half*(dens(id)+dens(id-1))*
     *                           (zd(id)-zd(id-1))
         WRITE(6,602) ID,DM(ID),TAUROS(ID),
     *                TEMP(ID),ELEC(ID),PTOTAL(ID),
     *                ZD(ID),ABROSD(ID),ABPLAD(ID)
     *                ,dens(id),xdm
      end do
      end if
c
      ITGREY=-1
      AMUV0=DMVISC**(ZETA0+UN)
      AMUV1=UN-AMUV0
      DO 30 ID=1,ND
         PGS(ID)=DENS(ID)*VSND20
         IF(DM(ID).LE.DMVISC*DM(ND)) THEN
            VISCD(ID)=(UN-FRACTV)*(ZETA1+UN)/
     *                DMVISC**(ZETA1+UN)*(DM(ID)/DM(ND))**ZETA1
            THETA(ID)=(UN-FRACTV)*(DM(ID)/DMVISC/DM(ND))**(ZETA1+UN)
          ELSE
            VISCD(ID)=FRACTV*(ZETA0+UN)/AMUV1*
     *                (DM(ID)/DM(ND))**ZETA0
            THETA(ID)=(UN-FRACTV)+FRACTV*((DM(ID)/DM(ND))**(ZETA0+UN)-
     *                 AMUV0)/AMUV1
         END IF
         GAMJ(ID)=UN
C
C     First estimates of the values of the Rosseland opacity
C     and function tauthe
C
         IF(ID.EQ.1) THEN
            TAUR=DM(ID)*ABROS0
            TAUTHE(ID)=TAUR*THETA(ID)/(ZETA1+TWO)
            ABROSD(ID)=ABROS0
            ABPLAD(ID)=ABPLA0
         ELSE
            DDM=DM(ID)-DM(ID-1)
            TAUR=TAUROS(ID-1)+DDM*ABROSD(ID-1)
            TAUTHE(ID)=TAUTHE(ID-1)+DDM*ABROSD(ID-1)*THETA(ID)
            ABROSD(ID)=ABROSD(ID-1)
            ABPLAD(ID)=ABPLAD(ID-1)
         END IF
C   
         do 20 ii=1,nlevel
            wop(ii,id)=un
   20    continue
C
C     Determination of temperature and mean opacities
C
         CALL TEMPER(ID,TAUR,ITGREY)
   30 CONTINUE
C
      IF(IPRING.GE.2) THEN
      WRITE(6,601)
      xdm=dm(1)
      DO 40 ID=1,ND
         if(id.gt.1) xdm=xdm-half*(dens(id)+dens(id-1))*
     *                           (zd(id)-zd(id-1))
         WRITE(6,602) ID,DM(ID),TAUROS(ID),
     *                TEMP(ID),ELEC(ID),PTOTAL(ID),
     *                ZD(ID),ABROSD(ID),ABPLAD(ID)
     *                ,dens(id),xdm
   40 CONTINUE
      END IF
C
C
C     Simultaneous solution of the hydrostatic equilibrium and the
C     z-m relation, assuming sound speed fixed
C
      if(nconit.ge.0) CALL HESOLV
C
      if(nconit.lt.-2) then
         do id=2,nd
            dm(id)=dm(id-1)-half*(dens(id)+dens(id-1))*
     *                           (zd(id)-zd(id-1))
         end do
      end if
      IF(IPRING.GE.2) THEN
      xdm=dm(1)
      WRITE(6,601)
      DO 50 ID=1,ND
         if(id.gt.1) xdm=xdm-half*(dens(id)+dens(id-1))*
     *                           (zd(id)-zd(id-1))
         WRITE(6,602) ID,DM(ID),TAUROS(ID),
     *                TEMP(ID),ELEC(ID),PTOTAL(ID),
     *                ZD(ID),ABROSD(ID),ABPLAD(ID)
     *                ,dens(id),xdm
   50 CONTINUE
      END IF
C
C -------------------------------------------------------------------
C
C     Outer iteration loop for the pseudo-grey model;
C     basically generalized Unsold-Lucy procedure
C
C     1.iteration - assumes that
C           Rosseland opacity = flux mean opacity
C           Planck mean opac  = absorption mean opacity
C           Eddington factors = 1/3 and 1/sgrt(3)
C
C     next iterations
C           improvement of mean opacities and Eddington factors;
C           corrections of the temperature (generalized Unsold-Lucy)
C
C
C     1.part
C
  100 ITGREY=ITGREY+1
C
      ANEREL=ELEC(1)/(DENS(1)/WMM(1)+ELEC(1))
      DO 110 ID=1,ND
         TAUR=TAUROS(ID)
         IF(ITGREY.GT.1) TAUR=TAUFLX(ID)
         CALL TEMPER(ID,TAUR,ITGREY)
  110 CONTINUE
C
C     Again simultaneous solution of the hydrostatic equilibrium
C     and the z-m relation, assuming sound speed fixed
C
      if(nconit.ge.0) CALL HESOLV
C
      IF(IPRING.GE.1) THEN
      WRITE(6,601)
      xdm=dm(1)
      DO 120 ID=1,ND
         WRITE(6,602) ID,DM(ID),TAUROS(ID),
     *                TEMP(ID),ELEC(ID),PTOTAL(ID),
     *                ZD(ID),ABROSD(ID),ABPLAD(ID)
     *                ,dens(id),xdm
  120 CONTINUE
      END IF
C
  601 FORMAT(1H1,' ID    DM     TAUROSS   TEMP        NE       P',
     * 8X,'ZD      ROSS.MEAN  PLANCK','   dens     dm0'/)
  602 FORMAT(1H ,I3,1P2D9.2,0PF11.0,1P3D9.2,2X,2D9.2,2x,2d9.2)
C
C     If required, modification of the depth scale (logarithmically
C     equidistant not in m, but in Tau(ross)
C
C     *****************************************************
C
      IF(NNEWD.GT.0.AND.ITGREY.EQ.0.AND.TAUROS(ND).GT.10.) THEN
         DO 125 III=1,NNEWD
            CALL NEWDM
  125    CONTINUE
      END IF
C
C     another modification of the depth scale 
C
      IF(NNEWD.LT.0) THEN
         DO 126 III=1,-NNEWD
            CALL NEWDMT
  126    CONTINUE
      END IF
C
      IF(HMIX0.GT.0.) THEN
         CALL CONTMD
         GO TO 170
      END IF
C
C     *****************************************************
C
C     If ITGMAX = 0  -  no iterative improvement of the pseudo-grey
C     model is required
C
      IF(ITGMAX.EQ.0) GO TO 170
      IF(ITGREY.EQ.0) ITGREY=1
C
C     Opacities and mean intensities in all frequency points ;
C     evaluation of appropriate integrals over frequency
C
      CALL RADTOT
C
C     Interpolation of TOTH and FLOPAC, which are determined by RTE
C     at the intermediate depth ponts DM(ID+1/2) to the grid DM
C
      DO 140 ID=2,ND-1
       A1=DM(ID+1)-DM(ID-1)
       A0=(DM(ID)-DM(ID-1))/A1
       A1=(DM(ID+1)-DM(ID))/A1
       TOTH(ID)=A0*TOTH(ID+1)+A1*TOTH(ID)
       FLOPAC(ID)=A0*FLOPAC(ID+1)+A1*FLOPAC(ID)
  140 CONTINUE
      TOTH(ND)=0.
      FLOPAC(ND)=FLOPAC(ND-1)
C
C     Determination of new temperature
C
      IF(IPRING.GE.1) WRITE(6,613) ITGREY
      DO 160 ID=1,ND
       HMECH=TOTF*(UN-THETA(ID))
       DFLUX=TOTH(ID)-HMECH
       FAKK=TOTK(ID)/TOTJ(ID)
       ABRAD=RDOPAC(ID)/DENS(ID)/TOTJ(ID)
       GAMJ(ID)=ABRAD/ABPLAD(ID)/FAKK*THIRD
       IF(ID.NE.ND) THEN
          ABFLX=FLOPAC(ID)/TOTH(ID)
        ELSE
          ABFLX=ABFLXM
         END IF
         abflx=abrosd(id)
       IF(ID.EQ.1) THEN
          FHH=TOTH(ID)/TOTJ(ID)
          GAMH=FAKK/FHH/5.7753D-1
            TAUFLX(ID)=ABFLX*DM(ID)
            TAUTHE(ID)=TAUFLX(ID)*THETA(ID)/(ZETA1+TWO)
          DFINT=TAUFLX(ID)*DFLUX
          DB0=FAKK/FHH*DFLUX
          ELSE
            IF(DM(ID).LE.DMVISC*DM(ND)) THEN
               ZETAD=ZETA1
             ELSE
               ZETAD=ZETA0
            END IF
            DDM=DM(ID)-DM(ID-1)
          A0=(ABFLXM*DM(ID)-ABFLX*DM(ID-1))/DDM/(ZETAD+TWO)
          A1=(ABFLX-ABFLXM)/DDM/(ZETAD+3.D0)
            TAUFLX(ID)=TAUFLX(ID-1)+DDM*HALF*(ABFLXM+ABFLX)
          TAUTHE(ID)=TAUTHE(ID-1)+
     *         A0*(THETA(ID)*DM(ID)-THETA(ID-1)*DM(ID-1))+
     *         A1*(THETA(ID)*DM(ID)**2-THETA(ID-1)*DM(ID-1)**2)
          DFINT=DFINT+DDM*HALF*(ABFLXM*DFLUXM+ABFLX*DFLUX)
         END IF
       ABFLXM=ABFLX
       DFLUXM=DFLUX
       IF(ITGMAX.LT.0) GO TO 150
C
C        generalized Unsold-Lucy procedure
C
       B0=FOUR*SIG4P*TEMP(ID)**4
       DIS=TOTF*VISCD(ID)/ABPLAD(ID)/DM(ND)
       DB1=ABRAD/ABPLAD(ID)*TOTJ(ID)-B0+DIS
       DB=DB1-3.D0*GAMJ(ID)*(DB0+DFINT)
       BNEW=FOUR*SIG4P*TEMP(ID)**4+DB
       TEMP(ID)=SQRT(SQRT(BNEW/4.D0/SIG4P))
       BREL=DB/B0
C
C        diagnostic output of iterative improvement
C
  150    R2=ABFLX/ABROSD(ID)
       R3=ABRAD/ABPLAD(ID)
         IF(IPRING.GE.1) WRITE(6,614) ID,FAKK,TAUROS(ID),
     *                ABROSD(ID),ABFLX,R2,
     *                ABRAD,ABPLAD(ID),R3,TOTH(ID),HMECH,BREL
  160 CONTINUE
C
  613 FORMAT(1H1,'ITERATIVE IMPROVEMENT,   ITGREY =',I2//
     * ' ID    FK    TAUROS     ABROS',4X,
     * 'ABFLUX    RATIO     ABRAD   ABPLA     RATIO',10X,'FLUX',
     *  7X,'MECH',4X,'DELTA(B)/B'/)
  614 FORMAT(1H ,I3,1P2D9.2,1X,3D9.2,1X,3D9.2,3X,2D13.5,3X,D10.2)
C
      IF(ITGREY.LE.IABS(ITGMAX)) GO TO 100
C
C     End of iteration loop for the pseudo-grey model
C -------------------------------------------------------------------
C
C
C     2. The final part
C     Interpolation of the computed model to the depth scale which is
C     going to be used in the subsequent - complete-linearization -
C     step of the model construction
C
C
C     First option - no interpolation
C
C     Second option - interpolation to the prescribed mass scale DM
C
  170 CONTINUE
      IF(IDEPTH.EQ.1) THEN
         CALL INTERP(DM,TEMP,DM0,TEMP0,ND,ND0,2,1,0)
         CALL INTERP(DM,ELEC,DM0,ELEC0,ND,ND0,2,1,1)
         CALL INTERP(DM,DENS,DM0,DENS0,ND,ND0,2,1,1)
         CALL INTERP(DM,ZD,DM0,ZD0,ND,ND0,2,1,0)
      END IF
C
C     Third option - logarithmically equidistant Rosseland opt.depths
C
      IF(IDEPTH.EQ.2) THEN
         READ(IBUFF,*) TAU1
         TAU0(ND0)=TAUROS(ND)
         TAU2=TAU0(ND0)*0.99
         DML0=LOG(TAU1)
         DLGM=(LOG(TAU2)-DML0)/(ND0-2)
         DO 180 I=1,ND0-1
  180       TAU0(I)=EXP(DML0+(I-1)*DLGM)
         CALL INTERP(TAUROS,DM,TAU0,DM0,ND,ND0,2,1,0)
         CALL INTERP(TAUROS,TEMP,TAU0,TEMP0,ND,ND0,2,1,0)
         CALL INTERP(TAUROS,ELEC,TAU0,ELEC0,ND,ND0,2,1,1)
         CALL INTERP(TAUROS,DENS,TAU0,DENS0,ND,ND0,2,1,1)
         CALL INTERP(TAUROS,ZD,TAU0,ZD0,ND,ND0,2,1,0)
      END IF
C
C     Fourth option - truncation of the disk and computing only
C                     a disk atmosphere
C
      IF(IDEPTH.GE.3) THEN
         TAU1=TAUROS(1)
         IF(IDEPTH.EQ.3) THEN
            READ(IBUFF,*) TDIV
          ELSE IF(IDEPTH.EQ.4) THEN
            READ(IBUFF,*) TAU0(ND)
          ELSE IF(IDEPTH.EQ.5) THEN
            READ(IBUFF,*) TDIV
          ELSE IF(IDEPTH.EQ.6) THEN
            READ(IBUFF,*) TAU0(1),TAU0(ND)
         END IF
         IF(IDEPTH.EQ.3.OR.IDEPTH.EQ.5) THEN
            DO 182 ID=1,ND
               IF(TAUROS(ID).LE.TDIV.AND.TAUROS(ID+1).GT.TDIV)
     *         ID1=ID
  182       CONTINUE
         END IF
         if(tauros(nd).le.tdiv) ID1=ND
         IF(IDEPTH.EQ.3) THEN
            ND0=ID1
            DO 184 ID=1,ND0
               DM0(ID)=DM(ID)
               TEMP0(ID)=TEMP(ID)
               ELEC0(ID)=ELEC(ID)
               DENS0(ID)=DENS(ID)
               ZD0(ID)=ZD(ID)
  184       CONTINUE
            nd=nd0
          ELSE IF(IDEPTH.EQ.5) THEN
            TAU0(ND0)=TAUROS(ID1)
         END IF
         IF(IDEPTH.GE.4) THEN
            TAU2=TAU0(ND0)*0.99
            DML0=LOG(TAU1)
            DLGM=(LOG(TAU2)-DML0)/(ND0-2)
            DO 186 I=1,ND0-1
  186          TAU0(I)=EXP(DML0+(I-1)*DLGM)
            CALL INTERP(TAUROS,DM,TAU0,DM0,ND,ND0,2,1,0)
            CALL INTERP(TAUROS,TEMP,TAU0,TEMP0,ND,ND0,2,1,0)
            CALL INTERP(TAUROS,ELEC,TAU0,ELEC0,ND,ND0,2,1,1)
            CALL INTERP(TAUROS,DENS,TAU0,DENS0,ND,ND0,2,1,1)
            CALL INTERP(TAUROS,ZD,TAU0,ZD0,ND,ND0,2,1,0)
         END IF
c        IFZ0=-1
         ZND=ZD0(ND0)
         IF(INZD.GT.0) THEN
            INZD=0
            IF(INSE.GT.0) INSE=INSE-1
         END IF
      END IF
C
C     in the two last options - interpolation from the previous
C     Rosseland opacity scale to the new scale and from the previous
C     mass depth scale to the new one
C
      IF(IDEPTH.GT.0) THEN
         ND=ND0
         DO 190 I=1,ND
            DM(I)=DM0(I)
            TEMP(I)=TEMP0(I)
            ELEC(I)=ELEC0(I)
            DENS(I)=DENS0(I)
            ZD(I)=ZD0(I)
  190    CONTINUE
      END IF
C
C     Recalculation of the populations
C
      DO 200 ID=1,ND
       AN=DENS(ID)/WMM(ID)+ELEC(ID)
       ANEREL=ELEC(ID)/AN
       CALL ELDENS(ID,TEMP(ID),AN,ANE,ENRG)
         ELEC(ID)=ANE
         DENS(ID)=WMM(ID)*(AN-ANE)
         PGS(ID)=AN*BOLK*TEMP(ID)
       PHMOL(ID)=AHMOL
c        CALL WNSTOR(ID)
c        CALL STEQEQ(ID,POP,1)
  200 CONTINUE
      if(nconit.lt.0) CALL PSOLVE
      IF(HMIX0.GE.0.AND.IPRING.GT.0) CALL CONOUT(2,IPRING)
      LCHC=LCHC0
      RETURN
      END
C
C
C     ****************************************************************
C
C
 
      SUBROUTINE TEMPER(ID,TAUF,ITGR)
C     ===============================
C
C     Auxiliary procedure for LTEGR
C     Evaluation of temperature, electron density, Rosseland opacity
C     and Planck mean opacity for at a given depth point
C
C     Input parameters:
C     ID     - depth index
C     TAUF   - Rosseland optical depth (if ITGR = -1, 0 or 1
C            - flux mean opacity (if ITGR > 1)
C     ITGR   = -1, 0, 1 -  means that TEMPER is called in the first
C                          iteration of the pseudo-grey model;
C                          temperature is evaluated;
C            > 1      -    next iterations; temperature is given,
C                          only evaluation of electron density and
C                          populations
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'MODELQ.FOR'
      INCLUDE 'ALIPAR.FOR'
      COMMON/PRSAUX/VSND2(MDEPTH),HG1,HR1,RR1
      COMMON/FLXAUX/T4,PGAS,PRAD,PGM,PRADM,ITGMAX,ITGMX0
      COMMON/FACTRS/GAMJ(MDEPTH),GAMH,FAK0
      PARAMETER  (ERRT=1.D-3)
C
      IT=0
      IF(ITGR.GT.1.AND.ITGMAX.GT.0) THEN
       T=TEMP(ID)
       GO TO 10
      END IF
C
      IF(ID.EQ.1) THEN
         DDM=HALF*DM(ID)
       ELSE
         DDM=DM(ID)-DM(ID-1)
      END IF
C
C
C     Initial estimate of temperature for current values of Rosseland
C     and Planck mean opacities
C
      call tlocal(id,tauf,t)
C
C ********** Iteration loop for determining temperature at depth ID
C            for a given total pressure
C
   10 IT=IT+1
      TEMP(ID)=T
C
C     Estimate of the gas pressure
C
      PRAD=1.8912D-15*T4*(GAMH*5.7735D-1+TAUF-TAUTHE(ID))
      PTURB=HALF*DENS(ID)*VTURB(ID)*VTURB(ID)
      PGAS=PGS(ID)
      PTOT=PGAS+PRAD+PTURB
      PTOTAL(ID)=PTOT
      PRADT(ID)=PRAD
      IF(PGAS.LE.0.) WRITE(6,603) ID,IT,PGAS,PTOT,PTURB,PRAD
  603 format(' negative gas pressure!! id,it,pgas,p,pturb,prad =',
     *       2i3,1p4d9.2)
C
C     Determination of electron density from the known temperature
C     and total pressure
C
      AN=PGAS/T/BOLK
      CALL ELDENS(ID,T,AN,ANE,ENRG)
      ELEC(ID)=ANE
      DENS(ID)=WMM(ID)*(AN-ANE)
      PHMOL(ID)=AHMOL
      VSND2(ID)=PTOTAL(ID)/DENS(ID)
c     CALL WNSTOR(ID)
c     CALL STEQEQ(ID,POP,1)
      IF(IT.GT.1.AND.REL.LT.ERRT) GO TO 30
C
C     For itgr.gt.1 - only new electron density and populations
C
      IF(ITGR.GT.1) RETURN
C
C     Evaluation of the Rosseland and Planck mean opacities
C     for the new values of temperature, electron density, and
C     populations (OPROS - Rosseland opacity per 1 cm**3; OPPLA - Planck
C     mean opacity per 1 cm**3)
C
      CALL OPACF0(ID,NFREQ)
      CALL MEANOP(T,ABSO,SCAT,OPROS,OPPLA)
      ABROS=OPROS/DENS(ID)
      ABPLA=OPPLA/DENS(ID)
      if(abpla.lt.abpmin) abpla=abpmin
      ABFLX=ABROS
      ABPLAD(ID)=ABPLA
      ABROSD(ID)=ABROS
C
C     New values of the Rosseland opacity and function tauthe
C
      IF(ID.EQ.1) THEN
         TAUR=DM(ID)*ABROS
         TAUTHE(ID)=DM(ID)*ABFLX*THETA(ID)/(ZETA1+TWO)
       ELSE
         TAUR=TAUROS(ID-1)+DDM*HALF*(ABROSD(ID-1)+ABROS)
       ABFLXM=ABROSD(ID-1)
       ZETAD=ZETA0
       IF(DM(ID).LE.DMVISC*DM(ND)) ZETAD=ZETA1
       A0=(ABFLXM*DM(ID)-ABFLX*DM(ID-1))/DDM/(ZETAD+TWO)
       A1=(ABFLX-ABFLXM)/DDM/(ZETAD+3.D0)
         TAUTHE(ID)=TAUTHE(ID-1)+
     *              A0*(THETA(ID)*DM(ID)-THETA(ID-1)*DM(ID-1))+
     *              A1*(THETA(ID)*DM(ID)**2-THETA(ID-1)*DM(ID-1)**2)
      END IF
      TAUF=TAUR
C
C     New value of temperature
C
      call tlocal(id,tauf,t)
C
C     Convergence criterion for temperature
C     (if REL < 1e-3, temperature is not recalculated again, but for
C     consistency the electron density and pressures are still
C     calculated consistently with the last temperature)
C
      REL=ABS(T-TEMP(ID))/TEMP(ID)
      IF(IT.LE.5) GO TO 10
C
C     Store the final quantitites
C
   30 TEMP(ID)=T
      PGS(ID)=PGAS
      VSND2(ID)=PTOTAL(ID)/DENS(ID)
      ABROSD(ID)=ABROS
      ABPLAD(ID)=ABPLA
      TAUROS(ID)=TAUR
      TAUFLX(ID)=TAUF
      IF(ID.NE.1) RETURN
      DPRAD=1.8912D-15*T4*(TAUF-TAUTHE(ID))
      HG1=SQRT(2.D0*PGS(1)/DENS(1)/QGRAV)
      HR1=DPRAD/DM(1)/QGRAV
      RR1=HR1/HG1
      RETURN
      END
C
C
C     ****************************************************************
C
C
      subroutine tlocal(id,tauf,t)
c     ============================
c
c     local temperature as a a function of optical depth
c     for a grey model
c
c     input;   ID   - depth index
c              TAUF - current estimate of the flux-mean opacity
c     output:  T    - local temperature
c
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'MODELQ.FOR'
      parameter (c1=0.8112,c2=3.966e14,c3=6.745e-10,c4=0.96,
     *           c23=c2*c3,c34=c3*c4)
c
      COMMON/FLXAUX/T4,PGAS,PRAD,PGM,PRADM,ITGMAX,ITGMX0
      COMMON/FACTRS/GAMJ(MDEPTH),GAMH,FAK0
C
      if(tdisk.gt.0.) then
         t=tdisk
         return
      end if
c
      vis=viscd(id)/(3.*dm(nd))
      extra=4.*fak0*wdil*(tstar/teff)**4
      gj=gamj(id)
      gh=gamh*5.7735e-1
      gg=gh+tauf-tauthe(id)+extra
c
      if(icompt.eq.0.or.icomgr.eq.0) then
         t=(0.75*t4*(gj*gg+vis/abplad(id)))**0.25
         return
      end if
c
      epsbar=abplad(id)/abrosd(id)
      tfor=c1*teff*epsbar**(-0.125)
      tf0=tfor
      if(tauf.gt.un.and.tfor.lt.temp(id).or.tauf.ge.100.) then
         tfor=0.
         b=gg*(c3-c34) 
         b=0.
       else
         b=gg*c3
      end if
      a=epsbar/(0.75*t4)
      c=gg*(epsbar*gj+c34*tfor)+vis/abrosd(id)
      call quartc(a,b,c,t1)
      t=t1
      write(81,691) id,tauf,tauthe(id),
     *              abrosd(id),abplad(id),gg,epsbar,tf0,b,t
  691 format(i4,1p10e10.2)
c
      return
      end

C
C
C     ***************************************************************
C
C
 
      subroutine quartc(a,b,c,x)      
C     ==========================
c    
c     solver for the algebraic equation of the fourth order
c     a*x**4 + b*x =c
c
c     solution done by the Newton-Raphson method
c
c     Initial estimate
c
      INCLUDE 'IMPLIC.FOR'
c
      if(a.gt.b) then
         x=(c/a)**0.25
        else
        x=c/b
      end if
c          
      it=0
   10 continue
      it=it+1
      ax=a*x**3
      v=c-b*x-x*ax
      d=4.*ax+b
      if(d.ne.0.) dx=v/d
      x=x+dx
      if(abs(dx/x).gt.1.e-3) then
         if(it.lt.20) go to 10
       else
         if(it.ge.20) write(6,601) a,b,c,dx,x 
      end if
  601 format(' slow convergence of quartic solver'/
     *       ' a,b,c,dx,x = ',1p5e13.4)
      return
      end
      
C
C
C     ******************************************************************
C
C

      SUBROUTINE NEWDM
C     ================
C
C     New m-scale, calculated as that corresponding to a new
C     tau(Ross)-scale, which is logarithmically equidistant, with
C     a generally different step in six different regions:
C
C     1. region between the original tau(1) and T0 (taken as 0.01);
C     2. region between T0 and TC0 (taken as 0.01 and 0.1) - denser
C        mesh (with number of points N0)
C     3. region between TC0 and TC1 (taken as 0.1 and 10.) - the
C        central region with densest mesh (N1 points)
C     4. region between TC1 and T1 (taken as 10. and 100.) - as dens
C        as the second region (N0 points)
C     5. the remaining region between T1 and the original last tau
C        (more precisely, the last-but-one point TAUROS(ND-1).
C
C     If T1 is greater than the original last tau, the new tau-scale
C     is equidistant between T0 and the last tau.
C
C     The procedure also calulates all the necessary state parameters
C     for the new depth scale (density, z, pressure, opacities, and
C     temperature)
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'MODELQ.FOR'
      PARAMETER (TEN=1.D1)
      DIMENSION TAU(MDEPTH),TAUL(MDEPTH),DM0(MDEPTH),DENS0(MDEPTH),
     *          ABRS0(MDEPTH),ABPL0(MDEPTH)
      COMMON/PRSAUX/VSND2(MDEPTH),HG1,HR1,RR1
      COMMON/FACTRS/GAMJ(MDEPTH),GAMH,FAK0
      COMMON/FLXAUX/T4,PGAS,PRAD,PGM,PRADM,ITGMAX,ITGMX0
C
      DATA N0,NC0 /8,24/
      DATA T0,TC0,TC1,T1 /-2.D0,-1.D0,1.D0,2.D0/
C
      DO 10 ID=1,ND
         DM0(ID)=DM(ID)
         DENS0(ID)=DENS(ID)
       ABRS0(ID)=ABROSD(ID)
         ABPL0(ID)=ABPLAD(ID)
         TAUL(ID)=LOG10(TAUROS(ID))
         IF(TAUL(ID).LT.T0) IMIN=ID
         IF(TAUL(ID).LT.T1) IMAX=ID
   10 CONTINUE
      ND1=ND-1
      NC=2*N0+NC0
      NB=ND1-NC
      IF(IMAX.GE.ND1) THEN
         IC=0
       ELSE
         X=(TAUL(IMIN)-TAUL(1))/(TAUL(ND1)-TAUL(IMAX))
         X1=FLOAT(NB)/(X+UN)
         IC=int(X1)
      END IF
      NB0=NB-IC
C
C     New tau-scale
C
C     First, logarithmically equidistant tau-points between tau(1)
C     (which is the previous TAUROS(1)), and log tau = T0.
C     Their number is NB0.
c
      DT=(T0-TAUL(1))/FLOAT(NB0-1)
      TAU(1)=TAUL(1)
      DO 20 ID=2,NB0
         TAU(ID)=TAU(ID-1)+DT
   20 CONTINUE
C
      IF(IC.GT.0) THEN
C
C     2.region - between log tau = T0 and TC0
C
      DT=(TC0-T0)/FLOAT(N0)
      DO 40 I=1,N0
         TAU(NB0+I)=TAU(NB0+I-1)+DT
   40 CONTINUE
C
C     3. The most dense region between TC0 and TC1 (central part)
C
      NB1=NB0+N0
      DT=(TC1-TC0)/FLOAT(NC0)
      DO 50 I=1,NC0
         TAU(NB1+I)=TAU(NB1+I-1)+DT
   50 CONTINUE
C
C     4. The part similar to that between T0 and TC0, this time
C        betwen TC1 and T1
C
      NB2=NB1+NC0
      DT=(T1-TC1)/FLOAT(N0)
      DO 60 I=1,N0
         TAU(NB2+I)=TAU(NB2+I-1)+DT
   60 CONTINUE
C
C     5. The remainig part between T1 and the last-but-one tau
C
      NB3=NB2+N0
      DT=(TAUL(ND1)-T1)/FLOAT(IC)
      DO 70 I=1,IC
         TAU(NB3+I)=TAU(NB3+I-1)+DT
   70 CONTINUE
      TAU(ND)=TAUL(ND)
C
C     The case where the last tau is smaller than T1; in this case
C     the points are logarithmically equidistant between T0 and the
C     last-but-one tau
C
      ELSE
      DT=(TAUL(ND1)-T0)/FLOAT(NC)
      DO 80 I=1,NC
         TAU(NB0+I)=TAU(NB0+I-1)+DT
   80 CONTINUE
      TAU(ND)=TAUL(ND)
C
      END IF
C
C     ---------------------------------------
C
C     Final new Rosseland optical depth scale
C
      DO 100 ID=1,ND
         TAU(ID)=TEN**TAU(ID)
  100 CONTINUE
C
C     Interpolation from the old to the new tau(Ross) scale to
C     get the new m-scale, density and Planck mean opacity
C
      CALL INTERP(TAUROS,DM0,TAU,DM,ND,ND,2,1,1)
      CALL INTERP(DM0,DENS0,DM,DENS,ND,ND,2,1,1)
      CALL INTERP(DM0,ABRS0,DM,ABROSD,ND,ND,2,1,1)
      CALL INTERP(DM0,ABPL0,DM,ABPLAD,ND,ND,2,1,1)
C
C     New Rosseland opacity and functions theta and tauthe
C
      AMUV0=DMVISC**(ZETA0+UN)
      AMUV1=UN-AMUV0
      DO 110 ID=1,ND
         IF(DM(ID).LE.DMVISC*DM(ND)) THEN
            VISCD(ID)=(UN-FRACTV)*(ZETA1+UN)/
     *                DMVISC**(ZETA1+UN)*(DM(ID)/DM(ND))**ZETA1
            THETA(ID)=(UN-FRACTV)*(DM(ID)/DMVISC/DM(ND))**(ZETA1+UN)
          ELSE
            VISCD(ID)=FRACTV*(ZETA0+UN)/AMUV1*
     *                (DM(ID)/DM(ND))**ZETA0
            THETA(ID)=(UN-FRACTV)+FRACTV*((DM(ID)/DM(ND))**(ZETA0+UN)-
     *                 AMUV0)/AMUV1
         END IF
         GAMJ(ID)=UN
         IF(ID.EQ.1) THEN
            TAUROS(ID)=DM(ID)*ABROSD(ID)
            TAUTHE(ID)=TAUROS(ID)*THETA(ID)/(ZETA1+TWO)
          ANEREL=ELEC(ID)/(DENS(ID)/WMM(ID)+ELEC(ID))
         ELSE
            DDM=DM(ID)-DM(ID-1)
            TAUROS(ID)=TAUROS(ID-1)+DDM*HALF*(ABROSD(ID-1)+ABROSD(ID))
            ZETAD=ZETA0
            IF(DM(ID).LE.DMVISC*DM(ND)) ZETAD=ZETA1
          A0=(ABROSD(ID-1)*DM(ID)-ABROSD(ID)*DM(ID-1))/DDM/
     *         (ZETAD+TWO)
          A1=(ABROSD(ID)-ABROSD(ID-1))/DDM/(ZETAD+3.D0)
            TAUTHE(ID)=TAUTHE(ID-1)+
     *              A0*(THETA(ID)*DM(ID)-THETA(ID-1)*DM(ID-1))+
     *              A1*(THETA(ID)*DM(ID)**2-THETA(ID-1)*DM(ID-1)**2)
         END IF
         TAUR=TAUROS(ID)
         CALL TEMPER(ID,TAUR,1)
  110 CONTINUE
C
C     Next step - simultaneous solution of the hydrostatic
C     equilibrium and the z-m relation
C
      if(nconit.ge.0) CALL HESOLV
C
C     New temperature and mean opacities for the current density
C     and pressure
C
      DO 140 ID=1,ND
         TAUR=TAUROS(ID)
         CALL TEMPER(ID,TAUR,1)
  140 CONTINUE
C
C     Once again - simultaneous solution of the hydrostatic
C     equilibrium and the z-m relation
C
      if(nconit.ge.0) CALL HESOLV
C
      IF(IPRING.GE.1) THEN
      WRITE(6,601)
      DO 230 ID=1,ND
         WRITE(6,602) ID,DM(ID),TAUROS(ID),
     *                TEMP(ID),ELEC(ID),PTOTAL(ID),
     *                ZD(ID),ABROSD(ID),ABPLAD(ID)
  230 CONTINUE
      END IF
C
  601 FORMAT(1H1,' NEW DEPTH GRID ESTABLISHED, NEW MODEL:'/
     *           ' --------------------------------------'/
     * ' ID    DM     TAUROSS   TEMP     NE       P',
     * 8X,'ZD     ROSS.MEAN  PLANCK'/)
  602 FORMAT(1H ,I3,1P2D9.2,0PF8.0,1P3D9.2,2X,2D9.2)
C
      RETURN
      END
C
C
C     ****************************************************************
C
C
 
      SUBROUTINE NEWDMT
C     =================
C
C     New m-scale, calculated as that corresponding to the new
C     grid better representing temperature variations
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'MODELQ.FOR'
      DIMENSION DM0(MDEPTH),DM11(MDEPTH),DENS0(MDEPTH),ZD0(MDEPTH),
     *          T0(MDEPTH),T1(MDEPTH),ELEC0(MDEPTH),PT0(MDEPTH),
     *          ABRS0(MDEPTH),ABPL0(MDEPTH)
      COMMON/PRSAUX/VSND2(MDEPTH),HG1,HR1,RR1
      COMMON/FACTRS/GAMJ(MDEPTH),GAMH,FAK0
      COMMON/FLXAUX/T4,PGAS,PRAD,PGM,PRADM,ITGMAX,ITGMX0
C
      DO 10 ID=1,ND
         DM0(ID)=LOG10(DM(ID))
c        T0(ID)=LOG10(TEMP(ID))
         T0(ID)=LOG10(TAUROS(ID))
         ELEC0(ID)=ELEC(ID)
         DENS0(ID)=DENS(ID)
         PT0(ID)=PTOTAL(ID)
         ZD0(ID)=ZD(ID)
       ABRS0(ID)=ABROSD(ID)
         ABPL0(ID)=ABPLAD(ID)
   10 CONTINUE
      ND1=ND-1
      CALL GRIDP(DM0,T0,DM11,T1,ND1)
      DM11(ND)=DM0(ND)
      T1(ND)=T0(ND)
      DO 20 ID=1,ND
         DM(ID)=EXP(2.3025851*DM11(ID))
C        TEMP(ID)=EXP(2.3025851*T1(ID))
         TAUROS(ID)=EXP(2.3025851*T1(ID))
   20 CONTINUE
      CALL INTERP(DM0,ELEC0,DM11,ELEC,ND,ND,2,0,1)
      CALL INTERP(DM0,DENS0,DM11,DENS,ND,ND,2,0,1)
      CALL INTERP(DM0,PT0,DM11,PTOTAL,ND,ND,2,0,1)
      CALL INTERP(DM0,ZD0,DM11,ZD,ND,ND,2,0,0)
      CALL INTERP(DM0,ABRS0,DM11,ABROSD,ND,ND,2,0,1)
      CALL INTERP(DM0,ABPL0,DM11,ABPLAD,ND,ND,2,0,1)
      DO 30 ID=1,ND
         VSND2(ID)=PTOTAL(ID)/DENS(ID)
   30 CONTINUE
C
C     New Rosseland opacity and functions theta and tauthe
C
      AMUV0=DMVISC**(ZETA0+UN)
      AMUV1=UN-AMUV0
      DO 40 ID=1,ND
         IF(DM(ID).LE.DMVISC*DM(ND)) THEN
            VISCD(ID)=(UN-FRACTV)*(ZETA1+UN)/
     *                DMVISC**(ZETA1+UN)*(DM(ID)/DM(ND))**ZETA1
            THETA(ID)=(UN-FRACTV)*(DM(ID)/DMVISC/DM(ND))**(ZETA1+UN)
          ELSE
            VISCD(ID)=FRACTV*(ZETA0+UN)/AMUV1*
     *                (DM(ID)/DM(ND))**ZETA0
            THETA(ID)=(UN-FRACTV)+FRACTV*((DM(ID)/DM(ND))**(ZETA0+UN)-
     *                 AMUV0)/AMUV1
         END IF
         GAMJ(ID)=UN
         IF(ID.EQ.1) THEN
            TAUROS(ID)=DM(ID)*ABROSD(ID)
            TAUTHE(ID)=TAUROS(ID)*THETA(ID)/(ZETA1+TWO)
          ANEREL=ELEC(ID)/(DENS(ID)/WMM(ID)+ELEC(ID))
         ELSE
            DDM=DM(ID)-DM(ID-1)
            TAUROS(ID)=TAUROS(ID-1)+DDM*HALF*(ABROSD(ID-1)+ABROSD(ID))
            ZETAD=ZETA0
            IF(DM(ID).LE.DMVISC*DM(ND)) ZETAD=ZETA1
          A0=(ABROSD(ID-1)*DM(ID)-ABROSD(ID)*DM(ID-1))/DDM/
     *         (ZETAD+TWO)
          A1=(ABROSD(ID)-ABROSD(ID-1))/DDM/(ZETAD+3.D0)
            TAUTHE(ID)=TAUTHE(ID-1)+
     *              A0*(THETA(ID)*DM(ID)-THETA(ID-1)*DM(ID-1))+
     *              A1*(THETA(ID)*DM(ID)**2-THETA(ID-1)*DM(ID-1)**2)
         END IF
         TAUR=TAUROS(ID)
         CALL TEMPER(ID,TAUR,1)
   40 CONTINUE
C
C     Next step - simultaneous solution of the hydrostatic
C     equilibrium and the z-m relation
C
      if(nconit.ge.0) CALL HESOLV
C
C     New temperature and mean opacities for the current density
C     and pressure
C
      DO 50 ID=1,ND
         TAUR=TAUROS(ID)
         CALL TEMPER(ID,TAUR,1)
   50 CONTINUE
C
C     Once again - simultaneous solution of the hydrostatic
C     equilibrium and the z-m relation
C
      if(nconit.ge.0) CALL HESOLV
C
      IF(IPRING.GE.1) THEN
      WRITE(6,601)
      DO 60 ID=1,ND
         WRITE(6,602) ID,DM(ID),TAUROS(ID),
     *                TEMP(ID),ELEC(ID),PTOTAL(ID),
     *                ZD(ID),ABROSD(ID),ABPLAD(ID)
   60 CONTINUE
      END IF
C
  601 FORMAT(1H1,' NEW DEPTH GRID ESTABLISHED, NEW MODEL:'/
     *           ' --------------------------------------'/
     * ' ID    DM     TAUROSS   TEMP     NE       P',
     * 8X,'ZD     ROSS.MEAN  PLANCK'/)
  602 FORMAT(1H ,I3,1P2D9.2,0PF8.0,1P3D9.2,2X,2D9.2)
C
      RETURN
      END
C
C
C     ****************************************************************
C
C
 
      SUBROUTINE GRIDP(X,Y,XNEW,YNEW,N)
c     =================================
c
c     evaluation of new grid points for a function; grid points
c     determined by dividing the curve Y=f(x) into n-1 equal segments;
c     the x-coordinates of the endpoints of the individual segments 
c     define new grid points
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      DIMENSION X(N),Y(N),XNEW(N),YNEW(N),Z(MDEPTH)
C
C     original segments - lengths (Z), and directional cosines (CD);
C     ZTOT - total length of the curve;
C     Z0 - length of a new segment
C
      ZTOT=0.
      DO 10 I=2,N
         Z(I-1)=SQRT((X(I)-X(I-1))**2+(Y(I)-Y(I-1))**2)
         ZTOT=ZTOT+Z(I-1)
   10 CONTINUE
      Z0=ZTOT/(N-1)
C
      ISEG=1
      XLAST=X(ISEG)
      YLAST=Y(ISEG)
      ZREST=Z(ISEG)
      ZREM=Z0
      IP=1
      XNEW(IP)=X(1)
      YNEW(IP)=Y(1)
C
   20 CONTINUE
      IF(ZREM.LT.ZREST) THEN
         ZREST=ZREST-ZREM
         XLAST=XLAST+ZREM*(X(ISEG+1)-X(ISEG))/Z(ISEG)
         YLAST=YLAST+ZREM*(Y(ISEG+1)-Y(ISEG))/Z(ISEG)
         IP=IP+1
         XNEW(IP)=XLAST
         YNEW(IP)=YLAST
         ZREM=Z0
         IF(IP.GE.N-1) GO TO 50
       ELSE
         ZREM=ZREM-ZREST
         ISEG=ISEG+1
         XLAST=X(ISEG)
         YLAST=Y(ISEG)
         ZREST=Z(ISEG)
      END IF
      GO TO 20
   50 XNEW(N)=X(N)
      YNEW(N)=Y(N)
      RETURN
      END
C
C
C     ****************************************************************
C
C
 
      SUBROUTINE HESOLV
C     =================
C
C     Solution of the coupled system of hydrostatic equilibrium equation
C     and the z-m relation; with a given (generally depth-dependent) sound
C     speed, defined as total pressure/density.
C     Numerical solution by a Newton-Raphson method
C
C     Input:  P     - initial total pressure
C             VSND2 - sound speed squared
C             HG1   - gas pressure scale height at the surface
C             RR1   - ratio of radiation and gas pressure scale heights
C                     at the surface
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'MODELQ.FOR'
      COMMON/PRSAUX/VSND2(MDEPTH),HG1,HR1,RR1
      DIMENSION P(MDEPTH),B(2,2),C(2,2),VL(2),
     *          D(2,2,MDEPTH),ANU(2,MDEPTH)
C
      DATA ERROR /1.D-4/
C
C     Density for a given total pressure and sound speed
C
      DO 10 ID=1,ND
       P(ID)=PTOTAL(ID)
c         DENS(ID)=P(ID)/VSND2(ID)
         vsnd2(id)=p(id)/dens(id)
   10 CONTINUE
C
C     Consistent z-values
C
c     ZD(ND)=ZND
c     DO 20 IID=1,ND-1
c        ID=ND-IID
c        ZD(ID)=ZD(ID+1)+HALF*(DM(ID+1)-DM(ID))*(UN/DENS(ID)+
c    *                   UN/DENS(ID+1))
c  20 CONTINUE
C
C     Basic Newton-Raphson iteration loop
C
      ITERH=0
   30 CONTINUE
      ITERH=ITERH+1
C
C     -------------------
C     Forward elimination
C     -------------------
C
C     Upper boundary condition
C
      ID=1
      X=ZD(1)/HG1-RR1
      IF(X.LT.3.) THEN
       IF(X.LT.0.) X=0.
         F1=1.772453851D0*EXP(X*X)*ERFCX(X)
       ELSE
       F1=(UN-HALF/X/X)/X
      END IF
      BET0=HALF/DENS(ID)/P(ID)
      BETP=HALF/DENS(ID+1)/P(ID+1)
      GAMA=UN/(DM(ID+1)-DM(ID))
      B(1,1)=F1
      B(1,2)=TWO*(X*F1-UN)*P(1)/HG1
      B(2,1)=BET0
      B(2,2)=GAMA
      C(1,1)=0.
      C(1,2)=0.
      C(2,1)=-BETP
      C(2,2)=GAMA
      VL(1)=DM(ID)*2.D0*VSND2(ID)/HG1-P(ID)*F1
      VL(2)=BET0*P(ID)+BETP*P(ID+1)-GAMA*(ZD(ID)-ZD(ID+1))
      ANU(1,ID)=0.
      ANU(2,ID)=0.
      CALL MATINV(B,2,2)
      DO 120 I=1,2
         DO 120 J=1,2
            S=0.
            DO 110 K=1,2
  110          S=S+B(I,K)*C(K,J)
            D(I,J,ID)=S
            ANU(I,ID)=ANU(I,ID)+B(I,J)*VL(J)
  120 CONTINUE
C
C     Normal depth points  1 < ID < ND
C
      DO 150 ID=2,ND-1
         BET0=BETP
         BETP=HALF/DENS(ID+1)/P(ID+1)
         GAMA=UN/(DM(ID+1)-DM(ID))
         DMD=HALF*(DM(ID+1)-DM(ID-1))
         AA=UN/(DM(ID)-DM(ID-1))/DMD
         CC=GAMA/DMD
         BB=AA+CC
         BQ=QGRAV/P(ID)/DENS(ID)
         B(1,1)=BB+BQ-AA*D(1,1,ID-1)
         B(1,2)=-AA*D(1,2,ID-1)
         B(2,1)=BET0
         B(2,2)=GAMA
         C(1,1)=CC
         C(1,2)=0.
         C(2,1)=-BETP
         C(2,2)=GAMA
         VL(1)=AA*P(ID-1)+CC*P(ID+1)-(BB-BQ)*P(ID)+AA*ANU(1,ID-1)
         VL(2)=BET0*P(ID)+BETP*P(ID+1)-GAMA*(ZD(ID)-ZD(ID+1))
         CALL MATINV(B,2,2)
         ANU(1,ID)=0.
         ANU(2,ID)=0.
         DO 140 I=1,2
            DO 140 J=1,2
               S=0.
               DO 130 K=1,2
  130             S=S+B(I,K)*C(K,J)
               D(I,J,ID)=S
               ANU(I,ID)=ANU(I,ID)+B(I,J)*VL(J)
  140    CONTINUE
  150 CONTINUE
C
C     Lower boundary condition
C
      ID=ND
      AA=TWO/(DM(ID)-DM(ID-1))**2
      BQ=QGRAV/P(ID)/DENS(ID)
      B(1,1)=AA+BQ-AA*D(1,1,ID-1)
      B(1,2)=-AA*D(1,2,ID-1)
      B(2,1)=0.
      B(2,2)=UN
      VL(1)=QGRAV/DENS(ID)+AA*(P(ID-1)-P(ID)+ANU(1,ID-1))
      VL(2)=0.
      CALL MATINV(B,2,2)
      ANU(1,ID)=0.
      ANU(2,ID)=0.
      DO 170 I=1,2
         DO 170 J=1,2
            S=0.
            DO 160 K=1,2
  160          S=S+B(I,K)*C(K,J)
            D(I,J,ID)=S
            ANU(I,ID)=ANU(I,ID)+B(I,J)*VL(J)
  170 CONTINUE
C
C     ------------
C     Backsolution
C     ------------
C
      P(ID)=P(ID)+ANU(1,ID)
      ZD(ID)=ZND
      CHMAXX=ABS(ANU(1,ID)/P(ID))
      DO 200 IID=1,ND-1
         ID=ND-IID
         DO 180 I=1,2
            DO 180 J=1,2
  180          ANU(I,ID)=ANU(I,ID)+D(I,J,ID)*ANU(J,ID+1)
       CH1=ANU(1,ID)/P(ID)
       CH2=ANU(2,ID)/ZD(ID)
       IF(ABS(CH1).GT.CHMAXX) CHMAXX=ABS(CH1)
       IF(ABS(CH2).GT.CHMAXX) CHMAXX=ABS(CH2)
       IF(CH1.LT.-0.9D0) CH1=-0.9D0
       IF(CH1.GT.9.D0) CH1=9.D0
       P(ID)=P(ID)*(UN+CH1)
  200 CONTINUE
C
C     Recalculate density for the new total pressure
C
      DO 210 ID=1,ND
         DENS(ID)=P(ID)/VSND2(ID)
  210 CONTINUE
C
C     New  z-values
C
      ZD(ND)=ZND
      DO 220 IID=1,ND-1
         ID=ND-IID
         ZD(ID)=ZD(ID+1)+HALF*(DM(ID+1)-DM(ID))*(UN/DENS(ID)+
     *                   UN/DENS(ID+1))
  220 CONTINUE
C
C     Convergence criterion for the Newton-Raphson method
C
      IF(IPRING.GE.1) WRITE(6,601) ITERH,CHMAXX
  601 FORMAT(/' solution of hydrostatic eq. + z-m relation:',
     *       'iter = ',I3,' max.rel.chan. =',1PD10.2)
      IF(CHMAXX.GT.ERROR.AND.ITERH.LT.10) GO TO 30
C
      DO ID=1,ND
       X=PGS(ID)/PTOTAL(ID)
       PTOTAL(ID)=P(ID)
       PGS(ID)=X*P(ID)
      END DO
C
C     Recalculation of the populations
C
      if(ih2p.ge.0) then
      ID=1
      ANEREL=ELEC(ID)/(DENS(ID)/WMM(ID)+ELEC(ID))
      DO ID=1,ND
         CALL RHONEN(ID,TEMP(ID),DENS(ID),AN,ANE)
         ELEC(ID)=ANE
c         PGS(ID)=AN*BOLK*TEMP(ID)
c        CALL WNSTOR(ID)
c        CALL STEQEQ(ID,POP,1)
      END DO
      end if
      RETURN
      END
C
C
C     ****************************************************************
C
C

      SUBROUTINE PSOLVE
C     =================
C
C     "formal" solution of the second-order equation for the total
C     pressure  -  d**2 P/d m**2 = Q/DENS;
C     with known density
C     the resulting tridiagonal system is solved by the standard
C     elimination
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'MODELQ.FOR'
      DIMENSION D(MDEPTH),ANU(MDEPTH)
C
C     forward elimination
C
      ID=1
      B=1.D0
      VL=PTOTAL(ID)
      D(ID)=0.
      ANU(ID)=VL/B
      DO 10 ID=2,ND-1
         DMD=HALF*(DM(ID+1)-DM(ID-1))
         A=UN/DMD/(DM(ID)-DM(ID-1))
         C=UN/DMD/(DM(ID+1)-DM(ID))
         B=A+C-A*D(ID-1)
         VL=QGRAV/DENS(ID)
         D(ID)=C/B
         ANU(ID)=(VL+A*ANU(ID-1))/B
   10 CONTINUE
      ID=ND
      A=TWO/(DM(ID)-DM(ID-1))**2
      B=A-A*D(ID-1)
      VL=QGRAV/DENS(ID)            
      ANU(ID)=(VL+A*ANU(ID-1))/B
      PTOTAL(ND)=ANU(ND)
C
C     backsubstitution
C
      DO 20 IID=1,ND-1
         ID=ND-IID
         PTOTAL(ID)=ANU(ID)+D(ID)*PTOTAL(ID+1)
   20 CONTINUE
      RETURN
      END
C
C
C     ****************************************************************
C
C
 
      SUBROUTINE ZMRHO(R,HG)
C     ============================
C
C     Initial estimate of DM, DENS, and ZD
C     by an approximate solution of the hydrostatic equilibrium
C     equation
C     Both gas and radiation pressure contribute;
C
C     Input: R     - ration of radiation to gas pressure scale heights
C            HG    - gas pressure scale height
C            DM1   - mass at the first depth point
C            DMTOT - mass at the last depth point (central plane)
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'MODELQ.FOR'
      PARAMETER (PISQ=1.77245385090551D0,pisq2=pisq*half)
C
C     Mass-depth grid - logarithmicaly equidistant
C
      if(nd.gt.mdepth) nd=mdepth
C
      if(dm1.gt.0) then
         DM(ND)=DMTOT
         DM(ND-1)=0.99D0*DMTOT
         DML=LOG(DM(ND-1)/DM1)/(ND-2)
         DML1=LOG(DM1)
         DO ID=1,ND-1
            DM(ID)=EXP(DML1+(ID-1)*DML)
         END DO
       else if(dm1.lt.-1.e-20.and.dm1.gt.-1.e-10) then
         DM(ND)=DMTOT
         DM1=ABS(DM1*1.e20)
         DML=LOG(DM(ND)/DM1)/(ND-1)
         DML1=LOG(DM1)
         DO ID=1,ND-1
            DM(ID)=EXP(DML1+(ID-1)*DML)
         end do
       else if(dm1.gt.-1.e-10) then
         if(mod(nd,2).eq.0) nd=nd-1
         dmha=dmtot*half
         dm1=abs(dm1*1.e10)
         ndha=nd/2
         dml=log(dmha/dm1)/ndha
         dml1=log(dm1)
         dm(nd)=dmtot
         do id=1,ndha
            dm(id)=exp(dml1+(id-1)*dml)
            dm(nd-id)=dm(nd)-exp(dml1+id*dml)
         end do
C
C        Determination of the total pressure scale height - function BETAH
C        HH - total pressure scale height
C
      end if
      HH=BETAH(R)*R
      DMH=PISQ/2.D0*EXP(-R*(HH-R))*ERFCX(HH-R)/HH
      RHO0=DM(ND)/HH/HG
C
C     Approximate solution of the hydrostatic equilibrium
C
      DO ID=1,ND
       DMREL=DM(ID)/DM(ND)
         IF(DMREL.LE.DMH) THEN
          X=R+ERFCIN(DMREL*2.D0/PISQ*HH*EXP(R*(HH-R)))
            RHO=EXP(-(X-R)*(X-R)-(HH-R)*R)
          ELSE IF(DMREL.LT.UN) THEN
          HSQ=SQRT(HH*(HH-R))
          X=ERFCIN(2.D0/PISQ*HSQ*(DMREL-DMH)+ERFCX(HSQ))*HH/HSQ
          RHO=EXP(-X*X*(UN-R/HH))
          ELSE
          X=0.
          RHO=UN
         END IF
         DENS(ID)=RHO0*RHO
         ZD(ID)=X*HG
      END DO
c
      if(dm1.lt.-1.e-10) then
c        dmha=dmtot*half
c        dm1=abs(dm1)
c        if(mod(nd,2).eq.0) nd=nd-1
c        ndha=nd/2
c        dml=log(dmha/dm1)/ndha
c        dml1=log(dm1)
c        dmdi=dmha/ndha
c        dm(nd)=dmtot
c        do id=1,ndha
c           dm(id)=exp(dml1+(id-1)*dml)
c           dm(nd-id)=dm(nd)-id*dmdi
c        end do
C
         DM(ND)=DMTOT
         DM(ND-1)=0.99D0*DMTOT
         dm1=abs(dm1)
         DML=LOG(DM(ND-1)/DM1)/(ND-2)
         DML1=LOG(DM1)
         DO ID=1,ND-1
            DM(ID)=EXP(DML1+(ID-1)*DML)
         end do
c
         hr=r*hg
         hg2=hg*pisq2
         hrg=hr+hg2
         dmh=hg2/hrg
         rho0=dmtot/hrg
         do id=1,nd
            dmrel=dm(id)/dm(nd)
            if(dmrel.ge.dmh) then
               zd(id)=hrg*(un-dmrel)
               dens(id)=rho0
             else
               zd(id)=hr+hg*erfcin(dmrel*hrg/hg2)
               x=(zd(id)-hr)/hg
               dens(id)=rho0*exp(-x*x)
            end if
         end do
      end if
      RETURN
      END
C
C
C     ********************************************************************
C
C
 
      FUNCTION BETAH(R)
C     =================
C
C     Determination of the total pressure scale height
C     Solution of the transcendental equation by the Newton-Raphson method
C
      INCLUDE 'IMPLIC.FOR'
      PARAMETER (UN=1.D0,
     *           PISQ=1.77245385090551D0)
      IF(R.LT.0.88) THEN
       BET0=PISQ/2.D0/R
       ELSE
       BET0=UN+UN/3.D0/R/R
      END IF
C
      ITER=0
      BETA=BET0
   10 ITER=ITER+1
      B1=BETA-UN
      RB1=R*B1
      BSQ=SQRT(BETA*B1)
      ERF1=ERFCX(R*BSQ)
      ERF2=ERFCX(RB1)
      RHS=BSQ/B1*(UN-ERF1)+EXP(-R*RB1)*ERF2
      DP=R/PISQ*(2.D0-EXP(-R*BETA*RB1))+(UN-ERF1)/2.D0/B1/BSQ+
     *   R*R*EXP(-R*RB1)*ERF2
      DBETA=(RHS-2.D0/PISQ*BETA*R)/DP
      DEL=DBETA/BETA
      BETA=BETA+DBETA
      IF(ABS(DEL).GT.1.D-5.AND.ITER.LE.10) GO TO 10
      BETAH=BETA
      RETURN
      END
C
C
C     ***********************************************************
C
C
      FUNCTION ERFCIN(X)
C     ==================
C
C     Approximate inverse complementary error function inverfc(x)
C
      INCLUDE 'IMPLIC.FOR'
      parameter (pisq=1.77245385090551D0, pisq2=pisq/2.d0)
      XL=-LOG(X)
      REL=0.88623D0+XL*(7.4871471D-3-XL*1.7726701D-4)
      E=SQRT(-LOG(X*(2.D0-X)))*REL
      iterr=0
   10 continue
      iterr=iterr+1
      dele=(erfcx(e)-x)*pisq2*exp(e*e)
      err=abs(dele/e)
      e=e+dele
      if(err.gt.1.e-6.and.iterr.lt.10) go to 10
      erfcin=e
      RETURN
      END
 
C
C
C     ****************************************************************
C

      SUBROUTINE RADTOT
C     =================
C
C     Evaluation of integrated radiative intensities and moments
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'ATOMIC.FOR'
      INCLUDE 'MODELQ.FOR'
      INCLUDE 'ODFPAR.FOR'
      INCLUDE 'ALIPAR.FOR'
      INCLUDE 'ITERAT.FOR'
      COMMON/OPTDPT/DT(MDEPTH)
      COMMON/SURFEX/EXTJ(MFREQ),EXTH(MFREQ)
      COMMON/TOTJHK/TOTJ(MDEPTH),TOTH(MDEPTH),TOTK(MDEPTH),
     *              RDOPAC(MDEPTH),FLOPAC(MDEPTH)
      DIMENSION SUMPL(MDEPTH)
C
C     zero the quantities
C
      DO ID=1,ND
         ABROSD(ID)=0.
         SUMDPL(ID)=0.
         ABPLAD(ID)=0.
         SUMPL(ID)=0.
         TOTJ(ID)=0.
         TOTH(ID)=0.
         TOTK(ID)=0.
       RDOPAC(ID)=0.
       FLOPAC(ID)=0.
         if(id.lt.nd) THEN
            DELDM(ID)=HALF*(DM(ID+1)-DM(ID))
            deldmz(id)=deldm(id)
            if(izscal.eq.1) deldmz(id)=half*(zd(id)-zd(id+1)) 
         end if
      END DO
      DEDM1=DM(1)/DENS(1)
C
C     loop over frequencies
C
      CALL OPAINI(1)
      DO 100 IJ=1,NFREQ
         FR=FREQ(IJ)
         CALL OPACF1(IJ)
         CALL RTEFR1(IJ)
         WW=W(IJ)
         DO ID=1,ND
            PLAN=XKFB(ID)/XKF1(ID)*WW
            DPLAN=PLAN/XKF1(ID)*FREQ(IJ)*HKT21(ID)
            ABROSD(ID)=ABROSD(ID)+DPLAN/ABSO1(ID)
            ABPLAD(ID)=ABPLAD(ID)+PLAN*(ABSO1(ID)-SCAT1(ID))
            SUMDPL(ID)=SUMDPL(ID)+DPLAN
            SUMPL(ID)=SUMPL(ID)+PLAN
            TOTJ(ID)=TOTJ(ID)+WW*RAD1(ID)
            TOTK(ID)=TOTK(ID)+WW*RAD1(ID)*FAK1(ID)
            RDOPAC(ID)=RDOPAC(ID)+WW*RAD1(ID)*(ABSO1(ID)-SCAT1(ID))
            IF(ID.LT.ND) THEN
               FLUX1=RAD1(ID+1)*FAK1(ID+1)-RAD1(ID)*FAK1(ID)
               TOTH(ID+1)=TOTH(ID+1)+WW*FLUX1/DT(ID)
            END IF
         END DO
         WF=WW*(FH(IJ)*RAD1(1)-HEXTRD(IJ))
c        WF=WW*FH(IJ)*RAD1(1)
         TOTH(1)=TOTH(1)+WF
         FLOPAC(1)=FLOPAC(1)+WF*ABSO1(1)/DENS(1)
c        write(6,601) ij,freq(ij),fh(ij),rad1(1),hextrd(ij),wf,toth(1),
c    *                flopac(1)
  100 CONTINUE
c 601 format(i4,1pe15.5,6e10.2)
c
c     Rosseland and Planck mean opacities
C
      DO ID=1,ND
         ABROSD(ID)=SUMDPL(ID)/ABROSD(ID)
         ABPLAD(ID)=ABPLAD(ID)/SUMPL(ID)
      END DO
C
c     Rosseland optical depth scale; flux mean
c
      ID=1
      TAUROS(ID)=HALF*DEDM1*ABROSD(ID)
      DO ID=2,ND         
         DTAUR=DELDM(ID-1)*(ABROSD(ID)*DENS1(ID)+
     *                      ABROSD(ID-1)*DENS1(ID-1))
       TAUROS(ID)=TAUROS(ID-1)+DTAUR
       FLOPAC(ID)=(TOTK(ID)-TOTK(ID-1))/(DM(ID)-DM(ID-1))
      END DO
c
c     final Rosseland and Planck mean opacities
C
      DO ID=1,ND
         ABROSD(ID)=ABROSD(ID)/DENS(ID)
         ABPLAD(ID)=ABPLAD(ID)/DENS(ID)
c        write(6,602) id,tauros(id),abrosd(id),abplad(id),flopac(id),
c    *                totj(id),toth(id),totk(id)
c 602 format(i4,1p7e10.2)
      END DO
C
      RETURN
      END
C
C
C     ****************************************************************
C
C

      SUBROUTINE OPACFA(IJ)
C     =====================
C
C     Absorption, emission, and scattering coefficients
C     at frequency IJ and for all depths
C
C     Saves additionally contributions per ion (for computing
C      ionic cooling and heating rates, see routine COOLRT)
C
C     Input: IJ   opacity and emissivity is calculated for the
C                 frequency points with index IJ
C     Output: ABSO1 -  array of absorption coefficient
C             EMIS1 -  array of emission coefficient
C             SCAT1 -  array of scattering coefficient (all scattering
C                       mechanisms except electron scattering)
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'ATOMIC.FOR'
      INCLUDE 'MODELQ.FOR'
      INCLUDE 'ODFPAR.FOR'
      INCLUDE 'ALIPAR.FOR'
      COMMON/COOLCO/ABSOTI(MION,MDEPTH),EMISTI(MION,MDEPTH),
     *               ABSOC1(MDEPTH),EMISC1(MDEPTH)
      PARAMETER (C14=2.99793D14, CFF1=1.3727D-25)
C
C     initialize
c
      IF(ICOMPT.GT.0) THEN
         DO ID=1,ND
            ELSCAT(ID)=ELEC(ID)*SIGEC(IJ)
         END DO
      END IF
C
      DO ID=1,ND
         ABSO1(ID)=ELSCAT(ID)
         EMIS1(ID)=0.
         SCAT1(ID)=ELSCAT(ID)
         ABSOC1(ID)=ABSO1(ID)
         EMISC1(ID)=0.
       DO ION=1,NION
          ABSOTI(ION,ID)=0.
          EMISTI(ION,ID)=0.
       END DO
      END DO
C
C     basic frequency- and depth-dependent quantities
C
      FR=FREQ(IJ)
      FRINV=UN/FR
      FR3INV=FRINV*FRINV*FRINV
      DO 10 ID=1,ND
         XKF(ID)=EXP(-HKT1(ID)*FR)
         XKF1(ID)=UN-XKF(ID)
         XKFB(ID)=XKF(ID)*BNUE(IJ)
   10 CONTINUE
C
C ********  1a. bound-free contribution - without dielectronic rec.
C
      if(ifdiel.eq.0) then
      DO IBFT=1,NTRANC
         ITR=ITRBF(IBFT)
         SG=CROSS(IBFT,IJ)
         IF(SG.GT.0.) THEN
         II=ILOW(ITR)
         JJ=IUP(ITR)
         IZZ=IZ(IEL(II))
         IMER=IMRG(II)
         DO ID=1,ND
            SGD=SG
            IF(MCDW(ITR).GT.0) THEN
               CALL DWNFR1(FR,FR0(ITR),ID,IZZ,DW1)
               DWF1(MCDW(ITR),ID)=DW1
               SGD=SG*DW1
            END IF
            IF(IFWOP(II).LT.0) THEN
               CALL SGMER1(FRINV,FR3INV,IMER,ID,SGME1)
               SGMG(IMER,ID)=SGME1
               SGD=SGME1
            END IF
            EMISBF=SGD*EMTRA(ITR,ID)
            ABSO1(ID)=ABSO1(ID)+SGD*ABTRA(ITR,ID)
            EMIS1(ID)=EMIS1(ID)+EMISBF
          ABSOTI(IEL(II),ID)=ABSOTI(IEL(II),ID)+SGD*ABTRA(ITR,ID)
          EMISTI(IEL(II),ID)=EMISTI(IEL(II),ID)+EMISBF
       END DO
       END IF
      END DO
      else
C
C ********  1b. bound-free contribution - with dielectronic rec.
C
      DO IBFT=1,NTRANC
         ITR=ITRBF(IBFT)
         SG=CROSS(IBFT,IJ)
         IF(SG.GT.0.) THEN
         II=ILOW(ITR)
         JJ=IUP(ITR)
         IZZ=IZ(IEL(II))
         IMER=IMRG(II)
         DO ID=1,ND
            SG=CROSSD(IBFT,IJ,ID)
            IF(SG.GT.0.) THEN
            SGD=SG
            IF(MCDW(ITR).GT.0) THEN
               CALL DWNFR1(FR,FR0(ITR),ID,IZZ,DW1)
               DWF1(MCDW(ITR),ID)=DW1
               SGD=SG*DW1
            END IF
            IF(IFWOP(II).LT.0) THEN
               CALL SGMER1(FRINV,FR3INV,IMER,ID,SGME1)
               SGMG(IMER,ID)=SGME1
               SGD=SGME1
            END IF
            EMISBF=SGD*EMTRA(ITR,ID)
            ABSO1(ID)=ABSO1(ID)+SGD*ABTRA(ITR,ID)
            EMIS1(ID)=EMIS1(ID)+EMISBF
          ABSOTI(IEL(II),ID)=ABSOTI(IEL(II),ID)+SGD*ABTRA(ITR,ID)
          EMISTI(IEL(II),ID)=EMISTI(IEL(II),ID)+EMISBF
          END IF
       END DO
       END IF
      END DO
      end if
C
C ******** 2. free-free contribution
C
      DO 40 ION=1,NION
         IT=ITRA(NNEXT(ION),NNEXT(ION))
C
C        hydrogenic with Gaunt factor = 1
C
         IF(IT.EQ.1) THEN
            DO ID=1,ND
               SF1=SFF3(ION,ID)*FR3INV
               SF2=SFF2(ION,ID)
               IF(FR.LT.FF(ION)) SF2=UN/XKF(ID)
               ABSOFF=SF1*SF2
               ABSO1(ID)=ABSO1(ID)+ABSOFF
               EMIS1(ID)=EMIS1(ID)+ABSOFF
             ABSOTI(ION,ID)=ABSOTI(ION,ID)+ABSOFF
             EMISTI(ION,ID)=EMISTI(ION,ID)+ABSOFF
            END DO
C
C         hydrogenic with exact Gaunt factor 
C
          ELSE IF(IT.EQ.2) THEN
            DO ID=1,ND
               SF1=SFF3(ION,ID)*FR3INV
               SF2=SFF2(ION,ID)
               IF(FR.LT.FF(ION)) SF2=UN/XKF(ID)
               X=C14*CHARG2(ION)/FR
               SF2=SF2-UN+GFREE1(ID,X)
               ABSOFF=SF1*SF2
               ABSO1(ID)=ABSO1(ID)+ABSOFF
               EMIS1(ID)=EMIS1(ID)+ABSOFF
             ABSOTI(ION,ID)=ABSOTI(ION,ID)+ABSOFF
             EMISTI(ION,ID)=EMISTI(ION,ID)+ABSOFF
            END DO
C
C         H minus free-free opacity
C
          ELSE IF(IT.EQ.3) THEN
            DO ID=1,ND
               ABSOFF=(CFF1+CFFT(ID)*FRINV)*CFFN(ID)*FRINV
               ABSO1(ID)=ABSO1(ID)+ABSOFF
               EMIS1(ID)=EMIS1(ID)+ABSOFF
             ABSOTI(ION,ID)=ABSOTI(ION,ID)+ABSOFF
             EMISTI(ION,ID)=EMISTI(ION,ID)+ABSOFF
            END DO
C
C         special evaluation of the cross-section
C
          ELSE IF(IT.LT.0) THEN
            DO ID=1,ND
               ABSOFF=FFCROS(ION,IT,TEMP(ID),FR)*
     *                POPUL(NNEXT(ION),ID)*ELEC(ID)
               ABSO1(ID)=ABSO1(ID)+ABSOFF
               EMIS1(ID)=EMIS1(ID)+ABSOFF
             ABSOTI(ION,ID)=ABSOTI(ION,ID)+ABSOFF
             EMISTI(ION,ID)=EMISTI(ION,ID)+ABSOFF
            END DO
         END IF
   40 CONTINUE
C
C     ********  3. - additional continuum opacity (OPADD)
C
      IF(IOPADD.NE.0) THEN
         ICALL=1
         DO ID=1,ND
            CALL OPADD(0,ICALL,IJ,ID)
            ABSO1(ID)=ABSO1(ID)+ABAD
            EMIS1(ID)=EMIS1(ID)+EMAD
            SCAT1(ID)=SCAT1(ID)+SCAD
          ABSOTI(IELH,ID)=ABSOTI(IELH,ID)+ABAD
          EMISTI(IELH,ID)=EMISTI(IELH,ID)+EMAD
         END DO
      END IF
C
      DO ID=1,ND
         ABSOC1(ID)=ABSO1(ID)
         EMISC1(ID)=EMIS1(ID)
      END DO
      IF(ICOOLP.EQ.0) GO TO 500
C
C ********  4. - opacity and emissivity in lines
C
      IF(ISPODF.EQ.0) THEN
      IF(IJLIN(IJ).GT.0) THEN
C
C     the "primary" line at the given frequency
C
         ITR=IJLIN(IJ)
       ION=IEL(ILOW(ITR))
         DO 50 ID=1,ND
            SG=PRFLIN(ID,IJ)
            ABSO1(ID)=ABSO1(ID)+SG*ABTRA(ITR,ID)
            EMIS1(ID)=EMIS1(ID)+SG*EMTRA(ITR,ID)
          ABSOTI(ION,ID)=ABSOTI(ION,ID)+SG*ABTRA(ITR,ID)
          EMISTI(ION,ID)=EMISTI(ION,ID)+SG*EMTRA(ITR,ID)
   50    CONTINUE
      ENDIF
      IF(NLINES(IJ).LE.0) GO TO 200
C
C     the "overlapping" lines at the given frequency
C
      DO 100 ILINT=1,NLINES(IJ)
         ITR=ITRLIN(ILINT,IJ)
       if(linexp(itr)) goto 100
         IJ0=IFR0(ITR)
         DO 60 IJT=IJ0,IFR1(ITR)
            IF(FREQ(IJT).LE.FR) THEN
               IJ0=IJT
               GO TO 70
            END IF
   60    CONTINUE
   70    IJ1=IJ0-1
         A1=(FR-FREQ(IJ0))/(FREQ(IJ1)-FREQ(IJ0))
         A2=UN-A1
       ION=IEL(ILOW(ITR))
         DO 80 ID=1,ND
            SG=A1*PRFLIN(ID,IJ1)+A2*PRFLIN(ID,IJ0)
            ABSO1(ID)=ABSO1(ID)+SG*ABTRA(ITR,ID)
            EMIS1(ID)=EMIS1(ID)+SG*EMTRA(ITR,ID)
          ABSOTI(ION,ID)=ABSOTI(ION,ID)+SG*ABTRA(ITR,ID)
          EMISTI(ION,ID)=EMISTI(ION,ID)+SG*EMTRA(ITR,ID)
   80    CONTINUE
  100 CONTINUE
  200 CONTINUE
C
C     Opacity sampling option
C
      ELSE
      IF(NLINES(IJ).LE.0) GO TO 400
      DO 300 ILINT=1,NLINES(IJ)
        ITR=ITRLIN(ILINT,IJ)
      ION=IEL(ILOW(ITR))
      KJ=IJ-IFR0(ITR)+KFR0(ITR)
      INDXPA=IABS(INDEXP(ITR))
      IF(INDXPA.NE.3 .AND. INDXPA.NE.4) THEN
        DO ID=1,ND
            SG=PRFLIN(ID,KJ)
            ABSO1(ID)=ABSO1(ID)+SG*ABTRA(ITR,ID)
            EMIS1(ID)=EMIS1(ID)+SG*EMTRA(ITR,ID)
          ABSOTI(ION,ID)=ABSOTI(ION,ID)+SG*ABTRA(ITR,ID)
          EMISTI(ION,ID)=EMISTI(ION,ID)+SG*EMTRA(ITR,ID)
        END DO
      ELSE
        DO ID=1,ND
          KJD=JIDI(ID)
c     SG=EXP(XJID(ID)*SIGFE(KJD,KJ)+(UN-XJID(ID))*SIGFE(KJD+1,KJ))
            ABSO1(ID)=ABSO1(ID)+SG*ABTRA(ITR,ID)
            EMIS1(ID)=EMIS1(ID)+SG*EMTRA(ITR,ID)
          ABSOTI(ION,ID)=ABSOTI(ION,ID)+SG*ABTRA(ITR,ID)
          EMISTI(ION,ID)=EMISTI(ION,ID)+SG*EMTRA(ITR,ID)
        END DO
      ENDIF
  300 CONTINUE
  400 CONTINUE
      ENDIF
  500 CONTINUE
C
C     ----------------------------
C     total opacity and emissivity
C     ----------------------------
C
      DO ID=1,ND
         ABSO1(ID)=ABSO1(ID)-EMIS1(ID)*XKF(ID)
         ABSOC1(ID)=ABSOC1(ID)-EMISC1(ID)*XKF(ID)
       DO ION=1,NION
         ABSOTI(ION,ID)=ABSOTI(ION,ID)-EMISTI(ION,ID)*XKF(ID)
       END DO
         EMIS1(ID)=EMIS1(ID)*XKFB(ID)
         EMISC1(ID)=EMISC1(ID)*XKFB(ID)
       DO ION=1,NION
         EMISTI(ION,ID)=EMISTI(ION,ID)*XKFB(ID)
       END DO
         absot(id)=abso1(id)
      END DO
      if(izscal.eq.0) then
         do id=1,nd
            absot(id)=abso1(id)*dens1(id)
         end do
      end if
c
      if(ifprd.gt.0) call prd(ij)
c
      RETURN
      END
C
C
C     ****************************************************************
C
C
      
      SUBROUTINE VISINI
C     =================
C
C     Auxiliary procedure for RESOLV
C     initialization of necessary quantities for treating the viscosity
C     in disks
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'ATOMIC.FOR'
      INCLUDE 'MODELQ.FOR'
      INCLUDE 'ITERAT.FOR'
C
      AMUV0=DMVISC**(ZETA0+UN)
      AMUV1=UN-AMUV0
      GP=0.
      GN=UN
      IF(INMP.GT.0) THEN
         GP=UN
         GN=0.
      END IF
C
      IF(IVISC.LE.1) THEN 
         x=0.
         DO ID=1,ND
            DMD=DM(1)
            IF(ID.GT.1) DMD=(DM(ID)+DM(ID-1))*HALF
            IF(DM(ID).LE.DMVISC*DM(ND)) THEN
               VISCD(ID)=(UN-FRACTV)*(ZETA1+UN)/
     *                   DMVISC**(ZETA1+UN)*(DM(ID)/DM(ND))**ZETA1
               THETAV(ID)=(UN-FRACTV)*(DMD/DMVISC/DM(ND))**(ZETA1+UN)
             ELSE
               VISCD(ID)=FRACTV*(ZETA0+UN)/AMUV1*
     *                   (DM(ID)/DM(ND))**ZETA0
               THETAV(ID)=(UN-FRACTV)+FRACTV*((DMD/DM(ND))**(ZETA0+UN)-
     *                    AMUV0)/AMUV1
             END IF
             TVISC(ID)=EDISC*VISCD(ID)*DENS(ID)
             DTVIST(ID)=0.
             DTVISR(ID)=EDISC*VISCD(ID)*WMM(ID)
             DTVISN(ID)=0.
             if(id.gt.1) then
c            X=X+HALF*(TVISC(ID)+TVISC(ID-1))*(DM(ID)-DM(ID-1))
             X=X+HALF*(TVISC(ID)/DENS(ID)+TVISC(ID-1)/DENS(ID-1))*
     *         (DM(ID)-DM(ID-1))
             end if
          END DO
          vtot=x
          x=0.
          if(iter.eq.niter) write(6,600)
          do id=1,nd
             AN=DENS(ID)/WMM(ID)+ELEC(ID)
             PGS(ID)=BOLK*TEMP(ID)*AN
             if(id.gt.1) then
c            X=X+HALF*(TVISC(ID)+TVISC(ID-1))*(DM(ID)-DM(ID-1))
             X=X+HALF*(TVISC(ID)/DENS(ID)+TVISC(ID-1)/DENS(ID-1))*
     *         (DM(ID)-DM(ID-1))
             end if
             alp=tvisc(id)/omeg32/pgs(id)*12.5664
             if(iter.eq.niter)
     *       write(6,601) id,dm(id),tvisc(id),thetav(id),x/vtot,
     *                    viscd(id),alp
             write(96,601) id,dm(id),tvisc(id),thetav(id),x/vtot,
     *                    viscd(id),alp
          if(id.eq.nd) 
     *    write(96,601) id,edisc,viscd(id),dens(id),pgs(id),omeg32
          end do
  601     format(i5,1p6e12.4)
       ELSE IF(IVISC.EQ.2) THEN
          X=0.
          THETAV(1)=0.
          DO ID=1,ND
             AN=DENS(ID)/WMM(ID)+ELEC(ID)
             PGS(ID)=BOLK*TEMP(ID)*AN
             TVISC(ID)=OMEG32*ALPHAV*PGS(ID)/12.5664
             DTVIST(ID)=TVISC(ID)/TEMP(ID)
             DTVISN(ID)=TVISC(ID)/AN
             DTVISR(ID)=0.
             if(id.gt.1) then
c            X=X+HALF*(TVISC(ID)+TVISC(ID-1))*(DM(ID)-DM(ID-1))
             X=X+HALF*(TVISC(ID)/DENS(ID)+TVISC(ID-1)/DENS(ID-1))*
     *         (DM(ID)-DM(ID-1))
             end if
          END DO
          VTOT=X
          X=0.
          write(6,602)
  600 format(/'    ID    DM          TVISC   THETAV(orig)   THETAV',
     *        '      VISCD       ALPHA'/)
  602 format(/'    ID    DM          TVISC     THETAV        PGAS',
     *        '       VISCD       ALPHA'/)
          DO ID=1,ND
             if(id.gt.1) then
c            X=X+HALF*(TVISC(ID)+TVISC(ID-1))*(DM(ID)-DM(ID-1))
             X=X+HALF*(TVISC(ID)/DENS(ID)+TVISC(ID-1)/DENS(ID-1))*
     *         (DM(ID)-DM(ID-1))
             end if
             THETAV(ID)=X/VTOT
             viscd(id)=tvisc(id)/dens(id)/edisc
             write(6,601) id,dm(id),tvisc(id),thetav(id),pgs(id),
     *                    viscd(id),alphav
             write(96,601) id,dm(id),tvisc(id),thetav(id),pgs(id),
     *                    viscd(id),alphav
             if(id.eq.nd) 
     *       write(96,601) id,edisc,viscd(id),dens(id),pgs(id),omeg32
          END DO
      END IF
      RETURN
      END
C
C
C     ****************************************************************
C
C

      SUBROUTINE GREYD
C     ================
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'MODELQ.FOR'
      INCLUDE 'ATOMIC.FOR'
      INCLUDE 'ALIPAR.FOR'
      PARAMETER (ERRM0=1.E-3, NTRM=50)
C
      CHI0=20.
      IF(TEFF.GT.10000.) THEN
         XION=2.
       ELSE IF(TEFF.LT.6000.) THEN
         XION=1.
       ELSE
         XION=1.+(TEFF-6000.)/4000.
      END IF
      ANEREL=1.-1./XION
C
      DMP=0.
      ID=1
      C1=BOLK*XION/WMM(ID)
      C2=3.1415926/2./QGRAV
      C3=SQRT(C1*C2)
      C4=WBARM*OMEG32/ALPHAV
      C5=C4/C1
      C6=C1/WBARM/OMEG32
C
      ITRM=0
   10 CONTINUE
      ITRM=ITRM+1
      C1=BOLK*XION/WMM(ID)
      C3=SQRT(C1*C2)
      C5=C4/C1
      T=(0.375*TEFF**4*CHI0*C5)**0.2
      DM0=C5/T
      RHO=DM0/SQRT(T)/C3
      TEMP(ID)=T
      DENS(ID)=RHO
      CALL RHONEN(ID,T,RHO,AN,ANE)
      ELEC(ID)=ANE
      XION=AN/(AN-ANE)
c     ANEREL=1.-1./XION
c     CALL WNSTOR(ID)
c     CALL STEQEQ(ID,POP,1)
C
C     evaluation of the Rosseland and Planck mean opacities
C     for the new values of temperature, electron density, and
C     populations (OPROS - Rosseland opacity per 1 cm**3; OPPLA - Planck
C     mean opacity per 1 cm**3)
C
      CALL OPACF0(ID,NFREQ)
      CALL MEANOP(T,ABSO,SCAT,OPROS,OPPLA)
      ABROS=OPROS/DENS(ID)
      ABPLA=OPPLA/DENS(ID)
      if(abpla.lt.abpmin) abpla=abpmin
      CHI0=(ABROS+chi0)/2.
      WRITE(6,601) ITRM,T,DM0,RHO,CHI0,abros,CHI0*DM0,XION,ANE
  601 FORMAT(i5,1p8e9.2)
      ERRM=ABS(DM0-DMP)/DM0
      DMP=DM0
      IF(ERRM.GT.ERRM0.AND.ITRM.LT.NTRM) GO TO 10
C
      DMTOT=DM0
      RETURN
      END
C
C
C     ****************************************************************
C
C
 
      SUBROUTINE RHONEN(ID,T,RHO,AN,ANE)
C     ====================================
C
C     Evaluation of the electron density and the total hydrogen
C     number density for a given total particle number density
C     and temperature;
C     by solving the set of Saha equations, charge conservation and
C     particle conservation equations (by a Newton-Raphson method)
C
C     Input parameters:
C     T    - temperature
C     RHO  - density                        
C
C     Output:
C     AN    - total particle number density
C     ANE   - electron density
C     ANP   - proton number density
C     AHTOT - total hydrogen number density
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'MODELQ.FOR'
      INCLUDE 'ATOMIC.FOR'
      DIMENSION R(3,3),S(3),P(3)
C
      QM=0.
      Q2=0.
      QP=0.
      Q=0.
      DQN=0.
      TK=BOLK*T
      THET=5.0404D3/T
C
C     Coefficients entering ionization (dissociation) balance of:
C     atomic hydrogen          - QH;
C     negative hydrogen ion    - QM   (considered only if IHM>0);
C     hydrogen molecule        - QP   (considered only if IH2>0);
C     ion of hydrogen molecule - Q2   (considered only if IH2P>0).
C
      IF(IATREF.EQ.IATH) THEN
      IF(IHM.GT.0)  QM=1.0353D-16/T/SQRT(T)*EXP(8762.9/T)
      IF(IH2P.GT.0) QP=TK*EXP((-11.206998+THET*(2.7942767+THET*
     *             (0.079196803-0.024790744*THET)))*2.30258509299405)
      IF(IH2.GT.0) Q2=TK*EXP((-12.533505+THET*(4.9251644+THET*
     *            (-0.056191273+0.0032687661*THET)))*2.30258509299405)
      QH=EXP((15.38287+1.5*LOG10(T)-13.595*THET)*2.30258509299405)
      END IF
C
C     Initial estimate of the electron density
C
      AN=RHO/WMM(ID)/(1.-ANEREL)
      ANE=AN*ANEREL
      IT=0
C
C     Basic Newton-Raphson loop - solution of the non-linear set
C     for the unknown vector P, consistiong of AH, ANH (neutral
C     hydrogen number density) and ANE.
C
   10 IT=IT+1
C
C     procedure STATE determines Q (and DQN) - the total charge (and its
C     derivative wrt temperature) due to ionization of all atoms which
C     are considered (both explicit and non-explicit), by solving the set
C     of Saha equations for the current values of T and ANE
C
      CALL STATE(1,ID,T,ANE)
C
C     Auxiliary parameters for evaluating the elements of matrix of
C     linearized equations.
C     Note that complexity of the matrix depends on whether the hydrogen
C     molecule is taken into account
C     Treatment of hydrogen ionization-dissociation is based on
C     Mihalas, in Methods in Comput. Phys. 7, p.10 (1967)
C
      IF(IATREF.EQ.IATH) THEN
      G2=QH/ANE
      G3=0.
      G4=0.
      G5=0.
      D=0.
      E=0.
      G3=QM*ANE
      A=UN+G2+G3
      D=G2-G3
      IF(IT.GT.1) GO TO 60
      IF(IH2.EQ.0.AND.IH2P.EQ.0) GO TO 40
      IF(IH2.EQ.0) GO TO 20
      E=G2*QP/Q2
      B=TWO*(UN+E)
      GG=ANE*Q2
      GO TO 30
   20 B=TWO
      E=UN
      GG=G2*ANE*QP
   30 C1=B*(GG*B+A*D)-E*A*A
      C2=A*(TWO*E+B*Q)-D*B
      C3=-E-B*Q
      F1=(SQRT(C2*C2-4.*C1*C3)-C2)*HALF/C1
      FE=F1*D+E*(UN-A*F1)/B+Q
      GO TO 50
   40 F1=UN/A
      FE=D/A+Q
   50 AH=ANE/FE
      ANH=AH*F1
   60 AE=ANH/ANE
      GG=AE*QP
      E=ANH*Q2
      B=ANH*QM
C
C     Matrix of the linearized system R, and the rhs vector S
C
      R(1,1)=YTOT(ID)
      R(1,2)=0.
      R(1,3)=0.
      R(2,1)=-Q
      R(2,2)=-D-TWO*GG
      R(2,3)=UN+B+AE*(G2+GG)-DQN*AH
      R(3,1)=-UN
      R(3,2)=A+4.*(E+GG)
      R(3,3)=B-AE*(G2+TWO*GG)
      S(1)=RHO/WMM(ID)-YTOT(ID)*AH
      S(2)=ANH*(D+GG)+Q*AH-ANE
      S(3)=AH-ANH*(A+TWO*(E+GG))
C
C     Solution of the linearized equations for the correction vector P
C
      CALL LINEQS(R,S,P,3,3)
C
C     New values of AH, ANH, and ANE
C
      AH=AH+P(1)
      ANH=ANH+P(2)
      DELNE=P(3)
      ANE=ANE+DELNE
      AN=YTOT(ID)*AH+ANE
C
C     hydrogen is not the reference atom
C
      ELSE
C
C     Matrix of the linearized system R, and the rhs vector S
C
      IF(IT.EQ.1) THEN
         ANE=AN*HALF
         AH=ANE/YTOT(ID)
      END IF
      R(1,1)=YTOT(ID)
      R(1,2)=0.
      R(2,1)=-Q-QREF
      R(2,2)=UN-(DQN+DQNR)*AH
      S(1)=RHO/WMM(ID)-YTOT(ID)*AH
      S(2)=(Q+QREF)*AH-ANE
C
C     Solution of the linearized equations for the correction vector P
C
      CALL LINEQS(R,S,P,2,3)
      AH=AH+P(1)
      DELNE=P(2)
      ANE=ANE+DELNE
      AN=YTOT(ID)*AH+ANE
      END IF
C
C     Convergence criterion
C
      IF(ANE.LE.0.) ANE=1.D-3*AN
      IF(ABS(DELNE/ANE).GT.1.D-3.AND.IT.LE.10) GO TO 10
C
C     ANEREL is the exact ratio betwen electron density and total
C     particle density, which is going to be used in the subseguent
C     call of ELDENS
C
      ANEREL=ANE/AN
      AHTOT=AH
      IF(IATREF.EQ.IATH) THEN
      ANP=ANH/ANE*QH
      END IF
C
      RETURN
      END

C
C
C     **************************************************************
C
C
      subroutine quasim(ij)
c     =====================
c
c     quasi-molecular opacity for Lyman alpha and beta
c
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'ATOMIC.FOR'
      INCLUDE 'MODELQ.FOR'
      common/quasun/tqmprf,iquasi,nunalp,nunbet,nungam,nunbal
      dimension sgd(mdepth)
c
      if(iquasi.le.0) return
      fr=freq(ij)
      wlam=2.997925e18/fr
      if(wlam.lt.911..or.wlam.gt.1727.) return
      ii=nfirst(ielh)
c
      do jup=2,iquasi+1
       jj=ii+1
       itr=itra(ii,jj)
       do id=1,nd
            anp=popul(nnext(ielh),id)
            t=temp(id)
            if(tqmprf.gt.0.) t=tqmprf
            call allard(wlam,t,popul(ii,id),anp,sg,1,jup)
          sgd(id)=sg
       end do
       if(ijlin(ij).eq.itr) then
          do id=1,nd
C             sgd(id)=sgd(id)-half*prflin(id,ij)
          end do
        else
          sg0=0.
          do ilint=1,nlines(ij)
             itt=itrlin(ilint,ij)
             if(itt.eq.itr) then
                  IJ0=IFR0(ITR)
                  DO 10 IJT=IJ0,IFR1(ITR)
                     IF(FREQ(IJT).LE.FR) THEN
                     IJ0=IJT
                     GO TO 20
                     END IF
   10             CONTINUE
   20             IJ1=IJ0-1
                  A1=(FR-FREQ(IJ0))/(FREQ(IJ1)-FREQ(IJ0))
                  A2=UN-A1
                do id=1,nd
                     SG0=A1*PRFLIN(ID,IJ1)+A2*PRFLIN(ID,IJ0)
C                   sgd(id)=sgd(id)-half*sg0
                end do
             end if
          end do
       end if
       do id=1,nd
          abso1(id)=abso1(id)+sgd(id)*abtra(itr,id)
            emis1(id)=emis1(id)+sgd(id)*emtra(itr,id)
       end do
      end do
c
      return
      end
C
C
C ********************************************************************
C
      subroutine getlal
c     =================
c
c     getlal reads in the profile functions for Lyman alpha, beta, gamma,
c     and Balmer alpha, including the quasi-molecular satellites;
c     valid for first and second order in neutral and ionized H density
c     modified routine provided originally by D. Koester
c           
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      parameter (NXMAX=1400,NNMAX=5,NTAMAX=6)
      common/quasun/tqmprf,iquasi,nunalp,nunbet,nungam,nunbal
      common /callarda/xlalp(NXMAX),plalp(NXMAX,NNMAX),stnnea,stncha,
     *     vneua,vchaa,nxalp,iwarna
      common /callardb/xlbet(NXMAX),plbet(NXMAX,NNMAX),stnneb,stnchb,
     *     vneub,vchab,nxbet,iwarnb
      common /callardg/xlgam(NXMAX),plgam(NXMAX,NNMAX),stnneg,stnchg,
     *     vneug,vchag,nxgam,iwarng
      common /callardc/xlbal(NXMAX),plbal(NXMAX,NNMAX),stnnec,stnchc,
     *     vneuc,vchac,nxbal,iwarnc
      common /calphatd/xlalpd(NXMAX,NTAMAX),plalpd(NXMAX,NNMAX,NTAMAX),
     *       stnead(ntamax),stnchd(ntamax),
     *       vneuad(ntamax),vchaad(ntamax),
     *       talpd(ntamax),nxalpd(ntamax),ntalpd
c
c     Lyman alpha
c
      nxalp=0
      if(nunalp.gt.0) then
         read(nunalp,*) nxalp,stnnea,stncha,vneua,vchaa
         do i=1,nxalp
            read(nunalp,*) xlalp(i),(plalp(i,j),j=1,NNMAX)
         end do
         close(nunalp)
         stnnea=10.0**stnnea
         stncha=10.0**stncha
         iwarna=0
       else if(nunalp.lt.0) then
c
c      input of temperature-dependent profile
c
         nualp=-nunalp
         read(nualp,*) ntalpd
         do it=1,ntalpd
            read(nualp,*) talpd(it)
            read(nualp,*) nxalpd(it),stnead(it),stnchd(it),vneuad(it),
     *                    vchaad(it)
            do i=1,nxalpd(it)
               read(nualp,*) xlalpd(i,it),(plalpd(i,j,it),j=1,NNMAX)
            end do
            stnead(it)=10.0**stnead(it)
            stnchd(it)=10.0**stnchd(it)
         end do
         close(nualp)
      end if
c
c     Lyman beta
c
      nxbet=0
      if(nunbet.gt.0) then
      read(nunbet,*) nxbet,stnneb,stnchb,vneub,vchab
      do i=1,nxbet
         read(nunbet,*) xlbet(i),(plbet(i,j),j=1,NNMAX)
      end do
      close(nunbet)
      stnneb=10.0**stnneb
      stnchb=10.0**stnchb
      iwarnb=0
      end if
c
c     Lyman gamma
c
      nxgam=0
      if(nungam.gt.0) then
      read(nungam,*) nxgam,stnneg,stnchg,vneug,vchag
      do i=1,nxgam
         read(nungam,*) xlgam(i),(plgam(i,j),j=1,NNMAX)
      end do
      close(nungam)
      stnneg=10.0**stnneg
      stnchg=10.0**stnchg
      iwarng=0
      end if
c
c     Balmer alpha
c
      nxbal=0
      if(nunbal.gt.0) then
      read(nunbal,*) nxbal,stnnec,stnchc,vneuc,vchac
      do i=1,nxbal
         read(nunbal,*) xlbal(i),(plbal(i,j),j=1,NNMAX)
      end do
      close(nunbal)
      stnnec=10.0**stnnec
      stnchc=10.0**stnchc
      iwarnc=0
      end if
      return
      end
c
C
C
C ********************************************************************
C
C
      subroutine allard(xl,t,hneutr,hcharg,prof,iq,jq)
c     ================================================
c
c     quasi-molecular opacity for Lyman alpha, beta, and Balmer alpha
c     modified routine provided originally by D. Koester
c
c     Input:  xl:  wavelength in [A]
c             hneutr:  neutral H particle density [cm-3]
c             hcharg: ionized H particle density [cm-3]
c             iq:   quantum number of the lower level
c             jq:   quantum number of the upper level;
c                   =2  -  Lyman alpha
c                   =3  -  Lyman beta
c     Output: prof:  Lyman alpha line profile, normalized to 1.0e8
c             if integrated over A;
c             It then renormalized by multiplying by
c             8.853e-29*lambda_0^2*f_ij
c
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      parameter (NXMAX=1400,NNMAX=5,NTAMAX=6)
      parameter (xnorma=8.8528e-29*1215.6*1215.6*0.41618,
     *           xnormb=8.8528e-29*1025.73*1025.7*0.0791,
     *           xnormg=8.8528e-29*972.53*972.53*0.0290,
     *           xnormc=8.8528e-29*6562.*6562.*0.6407)
      common /callarda/xlalp(NXMAX),plalp(NXMAX,NNMAX),stnnea,stncha,
     *     vneua,vchaa,nxalp,iwarna
      common /callardb/xlbet(NXMAX),plbet(NXMAX,NNMAX),stnneb,stnchb,
     *     vneub,vchab,nxbet,iwarnb
      common /callardg/xlgam(NXMAX),plgam(NXMAX,NNMAX),stnneg,stnchg,
     *     vneug,vchag,nxgam,iwarng
      common /callardc/xlbal(NXMAX),plbal(NXMAX,NNMAX),stnnec,stnchc,
     *     vneuc,vchac,nxbal,iwarnc
      common /calphatd/xlalpd(NXMAX,NTAMAX),plalpd(NXMAX,NNMAX,NTAMAX),
     *       stnead(ntamax),stnchd(ntamax),
     *       vneuad(ntamax),vchaad(ntamax),
     *       talpd(ntamax),nxalpd(ntamax),ntalpd
      common/quasun/tqmprf,iquasi,nunalp,nunbet,nungam,nunbal
c
      prof=0.
c
c     Lyman alpha
c
      if(iq.eq.1.and.jq.eq.2) then    
      if(nunalp.lt.0) then
      call allardt(xl,t,hneutr,hcharg,prof)
      else
      if(xl.lt.xlalp(1).or.xl.gt.xlalp(nxalp)) return
      vn1=hneutr/stnnea
      vn2=hcharg/stncha
      vns=vn1*vneua+vn2*vchaa
      if(iwarna.eq.0) then
         if(vn1*vneua.gt.0.3.or.vn2*vchaa.gt.0.3) then
            write(*,*) '          warning: density too high for',
     *           ' Lyman alpha expansion'
            iwarna=1
         endif
      endif
      vn11=vn1*vn1
      vn22=vn2*vn2
      vn12=vn1*vn2
      xnorm=1.0/(1.0+vns+0.5*vns*vns)
c
      jl=0
      ju=nxalp+1
 10   if(ju-jl.gt.1) then
         jm=(ju+jl)/2
         if((xlalp(nxalp).gt.xlalp(1)).eqv.(xl.gt.xlalp(jm))) then
            jl=jm
         else
            ju=jm
         endif
         go to 10
      endif
      j=jl
c
      if(j.eq.0) j=1
      if(j.eq.nxalp) j=j-1
      a1=(xl-xlalp(j))/(xlalp(j+1)-xlalp(j))
      p1=  vn1*((1.0-a1)*plalp(j,1)+a1*plalp(j+1,1))
      p11=vn11*((1.0-a1)*plalp(j,2)+a1*plalp(j+1,2))
      p2=  vn2*((1.0-a1)*plalp(j,3)+a1*plalp(j+1,3))
      p22=vn22*((1.0-a1)*plalp(j,4)+a1*plalp(j+1,4))
      p12=vn12*((1.0-a1)*plalp(j,5)+a1*plalp(j+1,5))
      prof=(p1+p2+p11+p22+p12)*xnorm*xnorma
      return
      end if
      end if
c
c     Lyman beta
c
      if(iq.eq.1.and.jq.eq.3) then
      if(nxbet.eq.0) return
      if(xl.lt.xlbet(1).or.xl.gt.xlbet(nxbet)) return
      vn1=hneutr/stnneb
      vn2=hcharg/stnchb
      vns=vn1*vneub+vn2*vchab
      if(iwarnb.eq.0) then
         if(vn1*vneub.gt.0.3.or.vn2*vchab.gt.0.3) then
            write(*,*) '          warning: density too high for',
     *           ' Lyman beta expansion'
            iwarnb=1
         endif
      endif
      vn11=vn1*vn1
      vn22=vn2*vn2
      vn12=vn1*vn2
      xnorm=1.0/(1.0+vns+0.5*vns*vns)
c
      jl=0
      ju=nxbet+1
 20   if(ju-jl.gt.1) then
         jm=(ju+jl)/2
         if((xlbet(nxbet).gt.xlbet(1)).eqv.(xl.gt.xlbet(jm))) then
            jl=jm
         else
            ju=jm
         endif
         go to 20
      endif
      j=jl
c
      if(j.eq.0) j=1
      if(j.eq.nxbet) j=j-1
      a1=(xl-xlbet(j))/(xlbet(j+1)-xlbet(j))
      p1=  vn1*((1.0-a1)*plbet(j,1)+a1*plbet(j+1,1))
      p11=vn11*((1.0-a1)*plbet(j,2)+a1*plbet(j+1,2))
      p2=  vn2*((1.0-a1)*plbet(j,3)+a1*plbet(j+1,3))
      p22=vn22*((1.0-a1)*plbet(j,4)+a1*plbet(j+1,4))
      p12=vn12*((1.0-a1)*plbet(j,5)+a1*plbet(j+1,5))
      prof=(p1+p2+p11+p22+p12)*xnorm*xnormb
      return
      end if
c
c     Lyman gamma
c
      if(iq.eq.1.and.jq.eq.4) then
      if(nxgam.eq.0) return
      if(xl.lt.xlgam(1).or.xl.gt.xlgam(nxgam)) return
      vn1=hneutr/stnneg
      vn2=hcharg/stnchg
      vns=vn1*vneug+vn2*vchag
      if(iwarng.eq.0) then
         if(vn1*vneug.gt.0.3.or.vn2*vchag.gt.0.3) then
            write(*,*) '          warning: density too high for',
     *           ' Lyman gamma expansion'
            iwarng=1
         endif
      endif
      vn11=vn1*vn1
      vn22=vn2*vn2
      vn12=vn1*vn2
      xnorm=1.0/(1.0+vns+0.5*vns*vns)
c
      jl=0
      ju=nxgam+1
 30   if(ju-jl.gt.1) then
         jm=(ju+jl)/2
         if((xlgam(nxgam).gt.xlgam(1)).eqv.(xl.gt.xlgam(jm))) then
            jl=jm
         else
            ju=jm
         endif
         go to 30
      endif
      j=jl
c
      if(j.eq.0) j=1
      if(j.eq.nxgam) j=j-1
      a1=(xl-xlgam(j))/(xlgam(j+1)-xlgam(j))
      p1=  vn1*((1.0-a1)*plgam(j,1)+a1*plgam(j+1,1))
      p11=vn11*((1.0-a1)*plgam(j,2)+a1*plgam(j+1,2))
      p2=  vn2*((1.0-a1)*plgam(j,3)+a1*plgam(j+1,3))
      p22=vn22*((1.0-a1)*plgam(j,4)+a1*plgam(j+1,4))
      p12=vn12*((1.0-a1)*plgam(j,5)+a1*plgam(j+1,5))
      prof=(p1+p2+p11+p22+p12)*xnorm*xnormg
      return
      end if
c
c     Balmer alpha
c
      if(iq.eq.2.and.jq.eq.3) then    
      if(xl.lt.xlbal(1).or.xl.gt.xlbal(nxbal)) return
c      vn1=hneutr/stnnec
      vn1=0.
      vn2=hcharg/stnchc
      vns=vn1*vneuc+vn2*vchac
      vn11=vn1*vn1
      vn22=vn2*vn2
      vn12=vn1*vn2
      xnorm=1.0/(1.0+vns+0.5*vns*vns)
c
      jl=0
      ju=nxbal+1
 40   if(ju-jl.gt.1) then
         jm=(ju+jl)/2
         if((xlbal(nxbal).gt.xlbal(1)).eqv.(xl.gt.xlbal(jm))) then
            jl=jm
         else
            ju=jm
         endif
         go to 40
      endif
      j=jl
c
      if(j.eq.0) j=1
      if(j.eq.nxbal) j=j-1
      a1=(xl-xlbal(j))/(xlbal(j+1)-xlbal(j))
      p1=  vn1*((1.0-a1)*plbal(j,1)+a1*plbal(j+1,1))
      p11=vn11*((1.0-a1)*plbal(j,2)+a1*plbal(j+1,2))
      p2=  vn2*((1.0-a1)*plbal(j,3)+a1*plbal(j+1,3))
      p22=vn22*((1.0-a1)*plbal(j,4)+a1*plbal(j+1,4))
      p12=vn12*((1.0-a1)*plbal(j,5)+a1*plbal(j+1,5))
      prof=(p1+p2+p11+p22+p12)*xnorm*xnormc
      end if
c     
      return
      end
C
C
C ********************************************************************
C
C
      subroutine allardt(xl,t,hneutr,hcharg,prof)
c     ===========================================
c
c     quasi-molecular opacity for Lyman alpha, with T-dependent
c     profile
c
c     Input:  xl:  wavelength in [A]
c             hneutr:  neutral H particle density [cm-3]
c             hcharg: ionized H particle density [cm-3]
c     Output: prof:  Lyman alpha line profile, normalized to 1.0e8
c             if integrated over A;
c             It then renormalized by multiplying by
c             8.853e-29*lambda_0^2*f_ij
c
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      parameter (NXMAX=1400,NNMAX=5,NTAMAX=6)
      parameter (xnorma=8.8528e-29*1215.6*1215.6*0.41618)
      common /calphatd/xlalpd(NXMAX,NTAMAX),plalpd(NXMAX,NNMAX,NTAMAX),
     *       stnead(ntamax),stnchd(ntamax),
     *       vneuad(ntamax),vchaad(ntamax),
     *       talpd(ntamax),nxalpd(ntamax),ntalpd
c
      prof=0.
c
c     find the two partial tables close to actual T
c
      it0=0
      do it=1,ntalpd
         it0=it
         if(t.lt.talpd(it)) then
            it0=it-1
            go to 10
         end if
      end do
   10 continue
      if(it0.eq.0) then
         it0=1
         go to 20
      end if
      if(it0.ge.ntalpd) then
         it0=ntalpd
         go to 20
      end if
      go to 30
   20 continue
c
      if(xl.lt.xlalpd(1,it0).or.xl.gt.xlalpd(nxalpd(it0),it0)) return
      vn1=hneutr/stnead(it0)
      vn2=hcharg/stnchd(it0)
      vns=vn1*vneuad(it0)+vn2*vchaad(it0)
      vn11=vn1*vn1
      vn22=vn2*vn2
      vn12=vn1*vn2
      xnorm=1.0/(1.0+vns+0.5*vns*vns)
c
      jl=0
      ju=nxalpd(it0)+1
 110  if(ju-jl.gt.1) then
         jm=(ju+jl)/2
         if(xl.gt.xlalpd(jm,it0)) then
            jl=jm
         else
            ju=jm
         endif
         go to 110
      endif
      j=jl
c
      if(j.eq.0) j=1
      if(j.eq.nxalpd(it0)) j=j-1
      a1=(xl-xlalpd(j,it0))/(xlalpd(j+1,it0)-xlalpd(j,it0))
      p1=  vn1*((1.0-a1)*plalpd(j,1,it0)+a1*plalpd(j+1,1,it0))
      p11=vn11*((1.0-a1)*plalpd(j,2,it0)+a1*plalpd(j+1,2,it0))
      p2=  vn2*((1.0-a1)*plalpd(j,3,it0)+a1*plalpd(j+1,3,it0))
      p22=vn22*((1.0-a1)*plalpd(j,4,it0)+a1*plalpd(j+1,4,it0))
      p12=vn12*((1.0-a1)*plalpd(j,5,it0)+a1*plalpd(j+1,5,it0))
      prof=(p1+p2+p11+p22+p12)*xnorm*xnorma
      return
c
   30 continue
c
c     interpolate in the tables for different T
c
c     the lower T
c
      if(xl.lt.xlalpd(1,it0).or.xl.gt.xlalpd(nxalpd(it0),it0)) return
      vn1=hneutr/stnead(it0)
      vn2=hcharg/stnchd(it0)
      vns=vn1*vneuad(it0)+vn2*vchaad(it0)
      vn11=vn1*vn1
      vn22=vn2*vn2
      vn12=vn1*vn2
      xnorm=1.0/(1.0+vns+0.5*vns*vns)
      jl=0
      ju=nxalpd(it0)+1
 120  if(ju-jl.gt.1) then
         jm=(ju+jl)/2
         if(xl.gt.xlalpd(jm,it0)) then
            jl=jm
         else
            ju=jm
         endif
         go to 120
      endif
      j=jl
c
      if(j.eq.0) j=1
      if(j.eq.nxalpd(it0)) j=j-1
      a1=(xl-xlalpd(j,it0))/(xlalpd(j+1,it0)-xlalpd(j,it0))
      p1=  vn1*((1.0-a1)*plalpd(j,1,it0)+a1*plalpd(j+1,1,it0))
      p11=vn11*((1.0-a1)*plalpd(j,2,it0)+a1*plalpd(j+1,2,it0))
      p2=  vn2*((1.0-a1)*plalpd(j,3,it0)+a1*plalpd(j+1,3,it0))
      p22=vn22*((1.0-a1)*plalpd(j,4,it0)+a1*plalpd(j+1,4,it0))
      p12=vn12*((1.0-a1)*plalpd(j,5,it0)+a1*plalpd(j+1,5,it0))
      prof0=(p1+p2+p11+p22+p12)*xnorm*xnorma
c
c     the higher T
c
      it0=it0+1
      if(xl.lt.xlalpd(1,it0).or.xl.gt.xlalpd(nxalpd(it0),it0)) return
      vn1=hneutr/stnead(it0)
      vn2=hcharg/stnchd(it0)
      vns=vn1*vneuad(it0)+vn2*vchaad(it0)
      vn11=vn1*vn1
      vn22=vn2*vn2
      vn12=vn1*vn2
      xnorm=1.0/(1.0+vns+0.5*vns*vns)
      jl=0
      ju=nxalpd(it0)+1
 130  if(ju-jl.gt.1) then
         jm=(ju+jl)/2
         if(xl.gt.xlalpd(jm,it0)) then
            jl=jm
         else
            ju=jm
         endif
         go to 130
      endif
      j=jl
c
      if(j.eq.0) j=1
      if(j.eq.nxalpd(it0)) j=j-1
      a1=(xl-xlalpd(j,it0))/(xlalpd(j+1,it0)-xlalpd(j,it0))
      p1=  vn1*((1.0-a1)*plalpd(j,1,it0)+a1*plalpd(j+1,1,it0))
      p11=vn11*((1.0-a1)*plalpd(j,2,it0)+a1*plalpd(j+1,2,it0))
      p2=  vn2*((1.0-a1)*plalpd(j,3,it0)+a1*plalpd(j+1,3,it0))
      p22=vn22*((1.0-a1)*plalpd(j,4,it0)+a1*plalpd(j+1,4,it0))
      p12=vn12*((1.0-a1)*plalpd(j,5,it0)+a1*plalpd(j+1,5,it0))
      prof1=(p1+p2+p11+p22+p12)*xnorm*xnorma
c
c     final profile coefficient
c
      dt=talpd(it0)-talpd(it0-1)
      prof=(prof0*(talpd(it0)-t)+prof1*(t-talpd(it0-1)))/dt
c
      return
      end
C
C
C     **************************************************************
C
C
      subroutine hedif
c       ================
c
c subroutine to calculate the depth dependent abundance profile for
c a layered H+He atmosphere.
c
      INCLUDE 'IMPLIC.FOR'                                                      
      INCLUDE 'BASICS.FOR'
      INCLUDE 'MODELQ.FOR'
      INCLUDE 'ATOMIC.FOR'
      common/hediff/ hcmass,radstr
c     real depth(mdepth+1),qs(mdepth+1),
      dimension depth(mdepth+1),qs(mdepth+1),
     *     ps(mdepth+1),gams(mdepth+1),
     *     abunds(mdepth+1),hms(mdepth+1)
c
      data smas,srad /1.9891e33,6.9599e10/
      data z1,z2,a1,a2 /1.,2.,1.,4./
      data bigg,pi  / 6.6732e-8,3.141592654/
c
c Set up starting values
c
        do id=1,nd
           depth(id+1)=dm(id)
        end do
        if(radstr.lt.1.e3) radius=radstr*srad
        if(hcmass.gt.1.e-10) hcmass=hcmass*1.e-13
c
      gam=1.e-30
      gams(1)=1.e-30
c
999    continue
      depth(1)=1.e-10
        q1=(depth(1)*4*pi*radius**2/smas)    
        p1=(q1*grav**2/(4*pi*bigg))
      ps(1)=p1
      qs(1)=q1
      hms(1)=0
      abunds(1)=0.0
      dpsl=-6
      hm=0.0
      do i=2,nd+1,1
            q2=(depth(i)*4*pi*radius**2/smas)    
            p2=(q2*grav**2/(4*pi*bigg))
          dp=p2-p1
          dlp=log(p2)-log(p1)
C          if(q2.ge.1.0e-19) gam=gam+raph(gam,z1,z2,a1,a2)*dlp
          gam=gam+raph(gam,z1,z2,a1,a2)*dlp
          abun0=gam
          hm=hm+(q2-q1)/((1+gam*a2/a1))
          p1=p2
          ps(i)=p2
          qs(i)=q2
          gams(i)=gam
          abunds(i)=abun0
          hms(i)=hm
        end do
c
      dh1=(log10(hcmass)-log10(hms(nd+1)))
      dh=hcmass/hms(nd+1)
      if(dh.ge.0.99)go to 99
      gam=gams(1)*1.1
      gams(1)=gam
      hm=0.0
      go to 999
c
c Now work backwards to get the full profiles
c
99      q1=(depth(nd+1)*4*pi*radius**2/smas)    
        p1=(q1*grav**2/(4*pi*bigg))
c
c       store new helium abundance and corresponding new YTOT, MMY, WMM
c
        write(6,600)
        do id=1,nd
           aheold=0.
           ahenew=abunds(id+1)
           if(iathe.gt.0) then
              aheold=abund(iathe,id)
              abund(iathe,id)=ahenew
           end if
           ytot(id)=ytot(id)-aheold+ahenew
           wmy(id)=wmy(id)+(ahenew-aheold)*4.003
           wmm(id)=wmy(id)*hmass/ytot(id)
           write(6,601) id,aheold,ahenew,ytot(id),wmy(id),wmm(id)
        end do
  600   format(' stratified helium abundance'/
     *'  id     He(old)     He(new)    ytot     wmm     wmy'/)
  601   format(i4,1p5e11.3)
c
        return
      end
C
C
C     ****************************************************************
C
C


      function raph(gam,z1,z2,a1,a2)
c     ==============================
c
c     auxiliary function for subroutine hedif
c
      INCLUDE 'IMPLIC.FOR'
c
      b=1+gam
      c=z1+z2*gam
      d=a1+a2*gam
      e=(1+z1)+gam*(1+z2)
      den=(c*d/(gam*b))+(d*(z1-z2)**2/(b*e))
      dnum=e*(a2*z1-a1*z2)+d*(z2-z1)
      dgam=dnum/den
      raph=dgam
      return
      end
C
C
C     ****************************************************************
C
C
      subroutine opa_ot(nd0,wl,te,de,op,et)
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      dimension te(mdepth),de(mdepth),op(mdepth),et(mdepth)
      call opa_ot_mock(nd0,wl,iopt,te,de,op,et)
      return
      end
C
c     ****************************************************************
c

      subroutine opa_ot_mock(nd0,wl,iopt,te,de,op,et)
c
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      dimension te(mdepth),de(mdepth),op(mdepth),et(mdepth)
c
      if(iopt.eq.1) return
      fr=2.997925e18/wl
      fr15=fr*1.e-15
      bn0=1.4743e-2*fr15**3
c
      do id=1,nd0
         den=2.*de(id)*2.7e-24
         op(id)=0.3*den
         et(id)=op(id)*bn0/(exp(hk*fr/te(id))-1.d0)
      end do
      return
      end
C
C
C ********************************************************************
C
C
      SUBROUTINE TABINI
C     =================
C
C     Initialization and reading of the opacity table for thermal processe
C     and Rayleigh scattering
c     raytab: scattering opacities in cm^2/gm at 5.0872638d14 Hz (sodium D)
c     (NOTE: Quantities in rayleigh.tab are in log_e)
C
c     tempvec: array of temperatures
c     rhovec: array of densities (gm/cm^3)
c     nu:     array of frequencies
c     table:  absorptive opacities in cm^2/gm
c     (NOTE:  Quantities in absorption.tab are in log_e)
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'MODELQ.FOR'
      common/intcff/yint(mfreq),jint(mfreq)
c
c     dimension frlt(mfrtab)
      character*(80) optable
c 
      if(ibinop.eq.0) then
         read(15,*,err=10,end=10) optable
         go to 20
   10    optable='absopac.dat'
   20    open(53,file=optable,status='old')
      end if
c
      if(ibinop.eq.0) then
         read(53,*) numfreq,numtemp,numrho
c        read(53,*)
c        read(53,*) (tempvec(i),i=1,numtemp)
c        read(53,*)
c        read(53,*) (rhovec(j),j=1,numrho)
       else      
         read(53) numfreq,numtemp,numrho
c        read(53) (tempvec(i),i=1,numtemp)
c        read(53) (rhovec(j),j=1,numrho)
      end if
      close(53)
c
      RETURN
      END
C
C
C ********************************************************************
C
C
      subroutine rayset
c     ===================
c  
c     set up a table of Rayleigh scattering opacity
c
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'MODELQ.FOR'
c
      do id=1,nd
       t=temp(id)
       rho=dens(id)
       RL = LOG(RHO) 
       DELTAR=(RL-RTAB1)/(RTAB2-RTAB1)*FLOAT(numrho-1)
       JR = 1 + IDINT(DELTAR)
       IF(JR.LT.1) JR = 1
       IF(JR.GT.(numrho-1)) JR = numrho-1
       r1i=rhovec(jr)
       r2i=rhovec(jr+1)
       dri=(RL-R1i)/(R2i-R1i)
       if(JR .eq. 1) dri = 0.d0
C
       TL=LOG(T)
       DELTAT=(TL-TTAB1)/(TTAB2-TTAB1)*FLOAT(numtemp-1)
       JT = 1 + IDINT(DELTAT)
       IF(JT.LT.1) JT = 1
       IF(JT.GT.numtemp-1) JT = numtemp-1
       t1i=tempvec(jt)
       t2i=tempvec(jt+1)
       dti=(TL-T1i)/(T2i-T1i)
       if(JT .eq. 1) dti = 0.d0
C
       opr1=raytab(jt,jr)+dti*
     *      (raytab(jt+1,jr)-raytab(jt,jr))
       opr2=raytab(jt,jr+1)+dti*
     *      (raytab(jt+1,jr+1)-raytab(jt,jr+1))
       opac=opr1+dri*(opr2-opr1)
       raysc(id) = exp(opac)
      end do
      return
      end
C
C
C ***********************************************************************
C
C
 
      SUBROUTINE RAYLEIGH(MODE,IJ,ID,SCR)
C     ===================================
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'ATOMIC.FOR'
      INCLUDE 'MODELQ.FOR'
      PARAMETER (FRRAY =  2.463D15,
     *           C18   =  2.997925D18,
     *           CR0   =  5.799D-13,
     *           CR1   =  1.422D-6,
     *           CR2   =  2.784D0)
      COMMON/RAYSCT/RCS(MFREQ)
C
      IF(MODE.EQ.0) THEN
         DO IK=1,NFREQ
            FRM=MIN(FREQ(IK),FRRAY)
            X=(C18/FRM)**2
            RCS(IK)=(CR0+(CR1+CR2/X)/X)/X/X
         END DO
       ELSE
         SCR=RCS(IJ)*DENS(ID)/(WMM(ID)*YTOT(ID))
      END IF
      RETURN
      END

C
C
C ********************************************************************
C
C
      SUBROUTINE OPCTAB(FR,IJ,ID,T,RHO,AB,SC,SCT,IGRAM)
C     =================================================
C
C     opacity for a given temperature and density computed 
C     by an interpolation of the precalculated opacity table
C
C     This is a simplified routine with all interpolations linear
C
C     Input:   FR  - frequency   (Hz)
C              T   - temperature (K)
C              RHO - density     (g cm^-3)
C      Outout: AB  - absorptive opacity  (per gram)
C              SC  - scattering opacity  (per gram) 
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'MODELQ.FOR'
      parameter (nmiean = 20,
     *           frray0 = 5.0872638d14)
C
      data inirea /0/
C
C      initialize opacity tables
C
c       if(inirea.eq.0) then
c          call tabini
c          if(ifrayl.gt.0) call rayset
c          inirea=1
c       end if
c
       jf=ij
       frij=fr
c
       RL = LOG(RHO) 
       DELTAR=(RL-RTAB1)/(RTAB2-RTAB1)*FLOAT(numrho-1)
       JR = 1 + IDINT(DELTAR)
       IF(JR.LT.1) JR = 1
       IF(JR.GT.(numrho-1)) JR = numrho-1
       r1i=rhovec(jr)
       r2i=rhovec(jr+1)
       dri=(RL-R1i)/(R2i-R1i)
       if(JR .eq. 1) dri = 0.d0
C
       TL=LOG(T)
       DELTAT=(TL-TTAB1)/(TTAB2-TTAB1)*FLOAT(numtemp-1)
       JT = 1 + IDINT(DELTAT)
       IF(JT.LT.1) JT = 1
       IF(JT.GT.numtemp-1) JT = numtemp-1
       t1i=tempvec(jt)
       t2i=tempvec(jt+1)
       dti=(TL-T1i)/(T2i-T1i)
       if(JT .eq. 1) dti = 0.d0
C
       opr1=absopac(jt,jr,jf)+dti*
     *      (absopac(jt+1,jr,jf)-absopac(jt,jr,jf))
       opr2=absopac(jt,jr+1,jf)+dti*
     *      (absopac(jt+1,jr+1,jf)-absopac(jt,jr+1,jf))
       opac=opr1+dri*(opr2-opr1)
       opac = exp(opac)
C
       AB=opac
C
C     ************************************************************
C      scattering
C     ************************************************************
C      1. Rayleighh scattering
C
       sc=0.
       if(ifrayl.lt.0) then
          sc=raysc(id)*(freq(jf)/frray0)**4
        else if(ifrayl.gt.0) then
          call rayleigh(1,ij,id,scr)
          sc=scr
       end if
       sct=sc
C
c      2. cloud scattering (not yet implemented)
c
       if(iter.le.0) return
c
c       ab=ab+abscld(id,jf)
c       sct=sct+scacld(id,jf,1)
c
       if(igram.eq.0) then
       ab = ab*dens(id)
       sc = sc*dens(id)
       sct=sct*dens(id)
       end if
c
       RETURN
       END
C
C
C ********************************************************************
C
C
      SUBROUTINE OPACT1(IJ)
C     =====================
C
C     Absorption, emission, and scattering coefficients
C     at frequency IJ and for all depths
C
C     Input: IJ   opacity and emissivity is calculated for the
C                 frequency points with index IJ
C     Output: ABSO1 -  array of absorption coefficient
C             EMIS1 -  array of emission coefficient
C             SCAT1 -  array of scattering coefficient
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'MODELQ.FOR'
      INCLUDE 'ALIPAR.FOR'
      common/hmolab/anh2(mdepth),anhm(mdepth)
C
      FR=FREQ(IJ)
      DO ID=1,ND
         T=TEMP(ID)
         RHO=DENS(ID)
         XKF(ID)=EXP(-HKT1(ID)*FR)
         XKF1(ID)=UN-XKF(ID)
         XKFB(ID)=XKF(ID)*BNUE(IJ)
         PLAN=XKFB(ID)/XKF1(ID)
         CALL OPCTAB(FR,IJ,ID,T,RHO,AB,SC,SCT,0)
         ABSO1(ID)=AB+sct
         SCAT1(ID)=SCT
         EMIS1(ID)=AB*PLAN
         ABSOT(ID)=ABSO1(ID)/DENS(ID)
      END DO
C
      RETURN
      END
C
C
C     ****************************************************************
C
C
       SUBROUTINE SETTRM
C      =================
C    
C      reads the equation-of-state tables for the pressure (P)
C      and entropy (S), as a function of T and rho;
C
C      stores P(rho,T) and S(rho,t) in arrays PL and SL
C
       INCLUDE 'IMPLIC.FOR'
       COMMON/THERM/SL(330,100),PL(330,100)
       COMMON/TABLTD/R1,R2,T1,T2,T12,T22,INDEX
       common/tdedge/redge,pedge(100),sedge(100),cvedge(100),
     &               cpedge(100),gammaedge(100),tedge(100)
       common/tdflag/JON
       parameter (RCON=8.31434E7) 
C
       open(58,file='stab.dat',status='old')
       open(59,file='ptab.dat',status='old')
C
       READ(58,*) YHEA,INDEX,R1,R2,T1,T2,T12,T22
       
       DO 200 JR = 1,INDEX
       DO 150 JQS=1,10
       JL = 1 + (JQS-1)*10
       JU = JL + 9
       READ(58,130) (SL(JR,JQ),JQ=JL,JU)
130    FORMAT(10F8.5)
150    CONTINUE
200    CONTINUE
C
       READ(59,*) YHEA,INDEX,R1,R2,T1,T2,T12,T22
       DO 300 JR=1,INDEX
       DO 160 JQP=1,10
       JL = 1 + (JQP-1)*10
       JU = JL + 9
       READ(59,130) (PL(JR,JQ),JQ=JL,JU)
160    CONTINUE
300    CONTINUE
C
       CLOSE(58)
       CLOSE(59)
c
c  Edge arrays
c
       r = 1.5d0*10.d0**r1
       tmin = 1.5d0*10.d0**t1
       tmax = 0.9d0*10.d0**t2
       redge = r
       do i = 1, 100
          t = t1 + (t2-t1)*dfloat(i-1)/dfloat(99)
          t = 10.d0**t
          t = min(tmax,max(t,tmin))
          tedge(i) = t
        rho=r
          CALL PRSENT(RHO*1.1,T,P1,S1)
          CALL PRSENT(RHO,T*1.1,P2,S2)
          CALL PRSENT(RHO,T,P0,S0)
          S1=RCON*S1
          S2=RCON*S2
          S0=RCON*S0
          DPDR=(P1-P0)/(.1*RHO)
          DPDT=(P2-P0)/(.1*T)
          DSDT=(S2-S0)/(.1*T)
          DSDR=(S1-S0)/(.1*RHO)
          DEN=DPDR*DSDT-DPDT*DSDR
          P=P0
          S=S0/RCON
          CV=T*DSDT
          CP=T*DEN/DPDR
          DQ=DSDT*P/(DEN*RHO)
          GAMMA=1.d0/DQ
c
          pedge(i) = p
          sedge(i) = s
          cvedge(i) = cv
          cpedge(i) = cp
          gammaedge(i) = gamma
          write(44,45) i,tedge(i),cvedge(i),cpedge(i),sedge(i),
     &                 gammaedge(i)
45        format(i4,5e14.5)
       enddo
c
       RETURN
       END
c
C
C
C     ***************************************************************
C
C
      FUNCTION RHOEOS(T,P)
C     ====================
C
C     equation of state - determining density from
C     temperature and pressure 
C
C     Input:   T - temperature (K)
C              P - total pressure (cgs)
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'MODELQ.FOR'
      parameter(wmol0=1.67333E-24/2.3)
      data inirea /0/
C
C      initialize EOS tables
C
       if(inirea.eq.0) then
          call settrm
          inirea=1
       end if
C
      AN=P/BOLK/T
      RHO=AN*wmol0
C
      niteos=0
10    niteos=niteos+1
      CALL PRSENT(RHO,T,P0,S0)
      CALL PRSENT(RHO*1.01,T,P1,S1)
      DPDR=(P1-P0)/(.01*RHO)
      DRXX=(P-P0)/DPDR/rho
c     if(drxx.lt.-0.9) drxx=-0.9
c     rho=rho*(un+drxx)
c      write(62,601) niteos,t,p,rho,drxx,dpdr
c  601 format(i5,f10.1,1p3e12.4)
      IF(ABS(DRXX).GT.1.d-5.and.niteos.lt.20) GO TO 10
C
      rhoeos=rho
      return
      end
C
C
C  ********************************************************************
C
C
 
      SUBROUTINE SETDRT
C     =================
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'MODELQ.FOR'
      COMMON/RHODER/DRHODT(MDEPTH)
      PARAMETER(DDTMIN=0., DDTPLU=0.001)
C
      DO ID=1,ND
         T=TEMP(ID)
         P=PTOTAL(ID)
         RHO1=RHOEOS(T*(UN-DDTMIN),P)
         RHO2=RHOEOS(T*(UN+DDTPLU),P)
         DRHODT(ID)=(RHO2-RHO1)/T/(DDTMIN+DDTPLU)
      END DO
      RETURN
      END
C
C
C  ********************************************************************
C
C
 
      SUBROUTINE TRMDRT(ID,T,P,HEATCP,DLRDLT,GRDADB,RHO)
C     ==================================================
C
C     Thermodynamic derivatives - based on statew equation and entropy 
C     tables
C
C     Input:  T     -  temperature
C             P     -  gas pressure
C
C     Output: HEATCP - specific heat at constant pressure
C             DLRDLT - d(ln rho)/d(ln T)
C             GRDADB - adiabatic gradient d(ln T)/d(ln P)_ad
C             RHO   -  density
C
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      COMMON/CC/DPDR,DPDT,DSDT,DSDR,CV,S,GAMMA
      COMMON/CONVOUT/CFLX(MDEPTH),VELCON(MDEPTH),GRADAD(MDEPTH),
     &               ENT(MDEPTH)
C
      parameter (RCON=8.31434E7) 
      parameter(wmol0=1.67333E-24/2.3)
      common/tdedge/redge,pedge(100),sedge(100),cvedge(100),
     &              cpedge(100),gammaedge(100),tedge(100)
      common/tdflag/JON
C
C     numerical evaluation of thermodynamic derivatives
C
      rho=rhoeos(t,p)
      drho=0.01*rho
      dt=0.01*t
      call prsent(rho,t,p0,s0)
      call prsent(rho+drho,t,p1,s1)
      call prsent(rho-drho,t,p2,s2)
      call prsent(rho,t+dt,p3,s3)
      call prsent(rho,t-dt,p4,s4)
      dpdr=(p1-p2)/(2.*drho)
      dpdt=(p3-p4)/(2.*dt)
      dsdr=(s1-s2)/(2.*drho)*rcon
      dsdt=(s3-s4)/(2.*dt)*rcon
      DEN=DPDR*DSDT-DPDT*DSDR
c
      if(jon .eq. 0) then
         HEATCV=T*DSDT
         HEATCP=T*DEN/DPDR
         DQ=DSDT*P/(DEN*RHO)
         GAMMA=1.d0/DQ
         DLRDLT = -RHO*DPDR/(T*DPDT)
         DLRDLT = 1.D0/DLRDLT
         GRDADB = -P/(HEATCP*RHO*T)*DLRDLT
         TDPT=T*DPDT
       else if(jon .ne. 0) then
         HEATCV = cvedge(JON) 
         HEATCP = cpedge(JON)
         DLRDLT = -1.d0
         GRDADB = -P/(HEATCP*RHO*T)*DLRDLT  ! 0.4d0
         GAMMA = gammaedge(JON)
      endif
C
      grdadb=p/t*(dsdr/(dsdr*dpdt-dsdt*dpdr))
      GRADAD(ID) = GRDADB
      ENT(ID)    = S0
C
c      write(61,601) jon,t,p,rho,s,dpdr,dpdt,heatcp,dlrdlt,grdadb
c  601 format(i4,1p9e10.2)
      RETURN
      END

C
C
C     ***************************************************************
C
C
       SUBROUTINE PRSENT(R,T,FP,FS)
C      ============================
C
C      interpolates pressure and entropy from tables
C
       INCLUDE 'IMPLIC.FOR'
       COMMON/THERM/SL(330,100),PL(330,100)
       COMMON/TABLTD/R1,R2,T1,T2,T12,T22,INDEX
       common/tdedge/redge,pedge(100),sedge(100),cvedge(100),
     &               cpedge(100),gammaedge(100),tedge(100)
       common/tdflag/JON
C
       JON=0
       RL = DLOG10(R) 
       ALPHA=T1+(RL-R1)/(R2-R1)*(T12-T1)
       BETA=T2-T1+((T22-T12)-(T2-T1))*(RL-R1)/(R2-R1)
       QL = (DLOG10(T) - ALPHA)/BETA
       DELTA=(RL-R1)/(R2-R1)*FLOAT(INDEX-1)
       JR = 1 + IDINT(DELTA)
       JQ = 1 + IDINT(99.*QL)
       IF(JR.LT.2) GO TO 300
       IF(JR.GT.(INDEX-1)) GO TO 300
       IF(JQ.LT.2) GO TO 300
       IF(JQ.GT.99) GO TO 300
       P = DELTA - (JR-1)
       Q = 99.*QL - (JQ-1)
C        interpolate:
       FS = 0.5D0*Q*(Q-1.D0)*SL(JR,JQ-1)
     1  + 0.5D0*P*(P-1.D0)*SL(JR-1,JQ)
     2  + (1.D0+P*Q-P*P-Q*Q)*SL(JR,JQ)
     3  + 0.5D0*P*(P-2.D0*Q+1.D0)*SL(JR+1,JQ)
     4  + 0.5D0*Q*(Q-2.D0*P+1.D0)*SL(JR,JQ+1)
     5  + P*Q*SL(JR+1,JQ+1)
        FS = 10.D0**FS
C
       FP = 0.5D0*Q*(Q-1.D0)*PL(JR,JQ-1)
     1  + 0.5D0*P*(P-1.D0)*PL(JR-1,JQ)
     2  + (1.D0+P*Q-P*P-Q*Q)*PL(JR,JQ)
     3  + 0.5D0*P*(P-2.D0*Q+1.D0)*PL(JR+1,JQ)
     4  + 0.5D0*Q*(Q-2.D0*P+1.D0)*PL(JR,JQ+1)
     5  + P*Q*PL(JR+1,JQ+1)
        FP = 10.D0**FP
      RETURN
C
C       off the table
C
300    CONTINUE
C
       write(60,*) ' Off the table!'
C
       JQ = min(98,max(JQ, 2))
       JON = JQ
       FP = pedge(JQ)*R*T/(redge*tedge(JQ))
       FS = sedge(JQ) + 1.d0/(gammaedge(JQ)-1.d0)*
     &      dlog(FP/pedge(JQ)*(redge/R)**gammaedge(JQ))
       write(60,*) JQ, R, T, FP, FS
C
       RETURN
       END
C
C
C     ***************************************************************
C
C
      SUBROUTINE TABIN2
C     =================
C
C     Initialization and reading of the opacity table for thermal processes
C     and Rayleigh scattering
c     raytab: scattering opacities in cm^2/gm AT 5.0872638d14 Hz (sodium D)
c     (NOTE: Quantities in rayleigh.tab are in log_e)
C
c     tempvec: array of temperatures
c     rhovec: array of densities (gm/cm^3)
c     nu:     array of frequencies
c     table:  absorptive opacities in cm^2/gm
c     (NOTE:  Quantities in absorption.tab are in log_e)
C
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'MODELQ.FOR'
c
      dimension absort(mfrtab)
      character*(80) optable
c 
      if(ibinop.eq.0) then
        read(15,*,err=10,end=10) optable
        go to 20
   10   optable='absopac.dat'
   20   open(53,file=optable,status='old')
      end if
c
      if(ibinop.eq.0) then
      read(53,*) numfreq,numtemp,numrho
      read(53,*)
      read(53,*) (tempvec(i),i=1,numtemp)
      read(53,*)
      read(53,*) (rhovec(j),j=1,numrho)
c 30  format(10f12.5)
       else      
        read(53) numfreq,numtemp,numrho
        read(53) (tempvec(i),i=1,numtemp)
        read(53) (rhovec(j),j=1,numrho)
      end if
c
      INDEX = numtemp
      R1 = rhovec(1)
      R2 = rhovec(numrho)
      T1 = tempvec(1)
      T2 = tempvec(numtemp)
      T12 = T1
      T22 = T2
c
      if(ibinop.eq.0) then
      do k = 1, numfreq
       read(53,*)
       read(53,*)
       read(53,*) frtab(k)
       do j = 1, numrho
            read(53,*) (absopac(i,j,k),i=1,numtemp)
       end do
      end do
      close(53)
       else
        do k = 1, numfreq
           read(53) frtab(k)
           do j = 1, numrho
              read(53) (absopac(i,j,k),i=1,numtemp)
           end do
        end do
        close(53)
      end if
c

      FR1 = dlog10(frtab(1))
      FR2 = dlog10(frtab(numfreq))
c
c     interpolate ABSOPAC from tabular to actual frequencies
c
      do i=1,numtemp
         do j=1,numrho
            do k=1,numfreq
               absort(k)=absopac(i,j,k)
            end do
            do k=1,nfreq
               fr=freq(k)
               FRL = LOG10(fr)
               DELTAF=(FRL-FR1)/(FR2-FR1)*FLOAT(numfreq-1)
               JF = 1 + IDINT(DELTAF)
               IF(JF.LT.1) JF = 1
               IF(JF.GT.numfreq-1) JF = numfreq-1
               xn1=log10(frtab(JF))
               xn2=log10(frtab(JF+1))
               djf=(FRL-xn1)/(xn2-xn1)
               if(JF .eq. 1) djf = 0.d0
               opac=absort(jf)+djf*(absort(jf+1)-absort(jf))
               absopac(i,j,k)=opac
            end do
         end do
      end do
c
      numfreq2 = numfreq
c
c     read Rayleigh scattering opacity table
c
      open(52,file='rayleigh.tab',status='old')
      read(52,*) numfreq2,numtemp,numrho
      read(52,*)
      read(52,*) (tempvec(i),i=1,numtemp)
      read(52,*)
      read(52,*) (rhovec(j),j=1,numrho)
      read(52,*)
      do j = 1, numrho
         read(52,*) (raytab(i,j),i=1,numtemp)
      end do
      close(52)
c
      RETURN
      END
C
C
C ********************************************************************
C
C
      subroutine moleq(id,tt,an,aein,ane,entt,energ)
c     ==============================================
c
c     calculation of the equilibrium state of atoms and molecules
c
c     Input:  id    - depth point
c             tt    - temperature [K]
c             an    - number density
c             aein  - initial estimate of the electron pressure
c
c     Output: ane    - electron density
c             entt   - entropy
c             energ  - internal energy
c
C     Output through common
c
c     Input data given in file tsuji.molec
c 
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'MODELQ.FOR'
      INCLUDE 'ATOMIC.FOR'
c
      character*128 MOLEC
      COMMON/COMFH1/C(600,5),PPMOL(600),APMLOG(600),
     *              XIP(100),CCOMP(100),UIIDUI(100),P(100),
     *              FP(100),XKP(100),EPS,SWITER,
     *              NELEM(5,600),NATO(5,600),MMAX(600),
     *              NELEMX(50),NMETAL,NMOLEC,NIMAX
      common/eospar/pfmol(500,mdepth),anmol(500,mdepth),
     *              pfato(100,mdepth),anato(100,mdepth),
     *              pfion(100,mdepth),anion(100,mdepth)
      common/entrop/entato(100,mdepth),ention(100,mdepth),
     *              entmol(500,mdepth),enttot(mdepth)
      common/hmolab/anh2(mdepth),anhm(mdepth)
      DIMENSION NATOMM(5),NELEMM(5),
     *          emass(100),uelem(100),ull(100),anden(500),
     *          aelem(100),epadd(12),abadd(12),emadd(12),
     *          iatnx(38),ammol(500),cmol(500)
      dimension denso(mdepth),eleco(mdepth),wmmo(mdepth)
c
      data nmetal/38/
      data iatnx/ 1, 2, 3, 4, 5, 6, 7, 8, 9,
     *            11,12,13,14,15,16,17,19,20,
     *            21,22,23,24,25,26,28,29,32,
     *            35,37,38,39,40,41,53,56,57,58,60/
      data epadd/7.899,11.840,4.176,5.692,6.380,6.840,
     *           6.880,10.454,5.210,5.610,5.470,5.490/
      data abadd/-8.59,-9.50,-9.40,-9.10,-9.76,-9.40,
     *           -10.6,-10.5,-9.87,-10.8,-10.4,-10.5/
      data emadd/72.59, 79.90, 85.47, 87.62, 88.91, 91.22,
     *           92.91,126.90,137.34,138.91,140.12,144.24/

      data iread/1/
c
c      if(tt.gt.1.e4.or.ifmol.eq.0) return  
      if(ifmol.eq.0) return  
      molec ='tsuji.molec'
c
        ECONST=4.342945E-1
        AVO=0.602217E+24
        SPA=0.196E-01
        GRA=0.275423E+05
        AHE=0.1
        tk=1./(tt*1.38054e-16)
        pgas=an/tk
        sahcon=1.87840e20*tt*sqrt(tt)
        nimax=3000
        eps=0.001
        switer=1
C
      if(iread.eq.1) then
c
      do i=1,nmetal
         ia=iatnx(i)
         nelemx(i)=ia
         if(ia.le.30) then
            ccomp(ia)=abndd(ia,id)
            xip(ia)=enev(ia,1)
            emass(ia)=amas(ia)
          else
            ccomp(ia)=exp(abadd(i-26)/econst)
            xip(ia)=epadd(i-26)
            emass(ia)=emadd(i-26)
         end if
      end do
c
c---- read molecular data from tsuji.molec ----------------------
c
        J=0
        OPEN(UNIT=26,FILE=MOLEC,STATUS='OLD')
  10    J=J+1
        READ (26,511) CMOL(J),(C(J,K),K=1,5),MMAX(J),
     *               (NELEMM(M),NATOMM(M),M=1,4)
 511    FORMAT (A8,E11.5,4E12.5,I1,(I2,I3),3(I2,I2))
c        
        MMAXJ=MMAX(J)
        IF(MMAXJ.EQ.0) GO TO 20
        DO M=1,MMAXJ
           NELEM(M,J)=NELEMM(M)
           NATO(M,J)=NATOMM(M)
        END DO
        GO TO 10
  20    CONTINUE
        NMOLEC=J-1
        close(26)
c        
        DO I=1,NMETAL
           NELEMI=NELEMX(I)
           P(NELEMI)=1.0E-20
        END DO
        iread=0
      endif
c
c---- end of reading atomic and molecular data  ----------------------
c
        p(99)= aein/tk
        pesave=p(99)
        p(99)=pesave
c
        THET=5040./tt
        TEM=tt
        PGLOG=log10(Pgas)
        PG=Pgas
        tkln25=-log(tk)*2.5
        tkln15=log(bolk*tt)*1.5
        tkev=5040./tt
c
        CALL RUSSEL(TEM,PG)
c
        PE=P(99)
        ane=pe*tk
        PELOG=log10(PE)
        emass(99)=5.486e-4
        uelem(99)=2.
        aelem(99)=pe*tk/(2.*sahcon*emass(nelemi)**1.5)
        ull(99)=log10(aelem(99))
c
        idstd=nd*2/3
c        entcon=83.23
        entcon=103.97
        ann=an-ane
c
c----atoms-----------------------------------------------------------------
c
        entt=0
        antt=0.
c        energ=1.5*an
        energ=0.
        tmass=0.
        DO I=1,NMETAL
           NELEMI=NELEMX(I)
           FPLOG=log10(FP(NELEMI))
           anden(i)=(p(nelemi)+1.e-20)*tk
           nelemi=nelemx(i)
         tmass=tmass+anden(i)*emass(nelemi)
           call mpartf(nelemi,1,0,tt,u0,dulog)
           uelem(nelemi)=u0
           aelem(nelemi)=anden(i)/(u0*sahcon*emass(nelemi)**1.5)
           ull(nelemi)=log10(aelem(nelemi))
           anato(nelemi,id)=anden(i) 
           pfato(nelemi,id)=u0
c           dulog=0.
           entato(nelemi,id)=tkln15-log(anden(i))+log(u0)+
     *                       1.5*log(emass(nelemi))+entcon+dulog*tt
           anx=anden(i)/ann
           antt=antt+anx
           entt=entt+entato(nelemi,id)*anx
           energ=energ+dulog/tk*anden(i)
           if(iprint.gt.3.and.mod(id,10).eq.1.and.nelemi.le.26)
     *     write(80,680) id,nelemi,anden(i),anx,antt,uelem(nelemi),
     *           entato(nelemi,id),entt,dulog,dulog/tk*anden(i),
     *           energ
        END DO
        an1=anden(1)
c
c---- positive ions ---------------------------------------------------------
c
        DO I=1,NMETAL
           NELEMI=NELEMX(I)
           PLOG= log10(P(NELEMI)+1.0D-30)
           XKPLOG=log10(XKP(NELEMI)+1.0D-30)
           PIONL=PLOG+XKPLOG-PELOG
           anden(i+nmetal)=exp(pionl/econst)*tk
         tmass=tmass+anden(i+nmetal)*emass(nelemi)
           call mpartf(nelemi,2,0,tt,u1,dulog)
c           dulog=0.
           anion(nelemi,id)=anden(i+nmetal)
           pfion(nelemi,id)=u1
           ent0=tkln25-pionl/econst+log(u1)+
     *                       1.5*log(emass(nelemi))+entcon+dulog/tk
           ention(nelemi,id)=tkln15-log(anden(i+nmetal))+log(u1)+
     *                       1.5*log(emass(nelemi))+entcon+dulog*tt
           anx=anion(nelemi,id)/ann
           antt=antt+anx
           entt=entt+ention(nelemi,id)*anx
           energ=energ+(xip(nelemi)*1.602e-12+dulog/tk)*anion(nelemi,id)
           enp=(xip(nelemi)*1.602e-12+dulog/tk)*anion(nelemi,id)
           ia=iatex(nelemi)
           abold=abund(ia,id)
           if(ia.gt.0) 
     *        abund(ia,id)=(anato(nelemi,id)+anion(nelemi,id))/
     *                     (anato(1,id)+anion(1,id))
           IF(iprint.gt.3.and.ID.EQ.IDSTD) THEN
              WRITE(71,671) NELEMI,ANATO(NELEMI,ID)/an,
     *                      ANION(NELEMI,ID)/an,
     *                      PFATO(NELEMI,ID),PFION(NELEMI,ID),
     *                      entato(nelemi,id)*anato(nelemi,id)/an,
     *                      entt
  671      FORMAT(I4,1P4E11.3,1p4e11.3)
           END IF
           if(iprint.gt.3.and.mod(id,10).eq.1.and.nelemi.le.26)
     *     write(80,680) id,nelemi,anden(i+nmetal),anx,antt,
     *           pfion(nelemi,id),
     *           ention(nelemi,id),entt,dulog,enp,energ
  680  format(2i4,1p11e11.3)
         END DO
c
c     H- (incorrect in Tsuji's table)
c
        j=1
        anmol(1,id)=1.0353e-16/tt/sqrt(tt)*exp(8762.9/tt)*
     *              anato(1,id)*ane
        pfmol(1,id)=1.
        entmol(1,id)=tkln15-log(anmol(1,id))+1.5*log(emass(1))+
     *               entcon
        anx=anmol(1,id)/ann
        antt=antt+anx
        entt=entt+entmol(j,id)*anx
        IF(iprint.gt.3.and.ID.EQ.IDSTD) 
     *  WRITE(72,672) J,CMOL(J),ANMOL(J,ID)/an,PFMOL(J,ID),
     *                   entmol(j,id)*anx,entt
c
c---- molecules-------------------------------------------------------------
c
      DO J=2,NMOLEC
         jm=j+2*nmetal
         PMOLL=log10(PPMOL(J)+1.0D-20)
           anden(jm)=exp(pmoll/econst)*tk
           umoll=1.
           if(pmoll.gt.-20.) then
              umoll=log10(anden(jm))+c(j,2)*thet
              amasm=0.
              do jjj=1,mmax(j)
                 i=nelem(jjj,j)
                 amasm=amasm+NATO(jjj,j)*emass(i)
                 umoll=umoll-NATO(jjj,j)*ull(i)
              end do
              ammol(j)=amasm
            tmass=tmass+anden(jm)*amasm
              umoll=exp(umoll/econst)/(sahcon*amasm**1.5)
c
c     replace with Irwin data whenever available
c
             call mpartf(0,0,j,tt,um,dulog)
           if(um.gt.0.) umoll=um
c
c     replace the TiO partition function by Kurucz-Schwenke data
c
            IF(J.eq.29) CALL TIOPF(TT,UMOLL)
           end if
           anmol(j,id)=anden(jm)
           pfmol(j,id)=umoll
c           dulog=0.
           entmol(j,id)=tkln15-log(anden(jm))+log(umoll)+
     *                  1.5*log(amasm)+entcon+dulog*tt
           anx=anden(jm)/ann
           antt=antt+anx
           entt=entt+entmol(j,id)*anx
           energ=energ+dulog/tk*anden(jm)
           IF(iprint.gt.3.and.ID.EQ.IDSTD) THEN
           WRITE(72,672) J,CMOL(J),ANMOL(J,ID)/an,PFMOL(J,ID),
     *                   entmol(j,id)*anden(jm)/an,dulog*tt*anden(jm)
  672      FORMAT(I4,2X,A8,1P2E12.3,1p2e11.3)
           END IF
           if(iprint.gt.3.and.mod(id,10).eq.1.and.j.le.10)
     *     write(80,680) id,j,anmol(j,id),anx,antt,pfmol(nelemi,id),
     *          entmol(j,id),entt,dulog,dulog/tk*anden(jm),energ
        END DO
c
c       electrons
c
        entel=tkln15-log(ane)+1.5*log(emass(99))+entcon
        entt=entt+entel*ane/ann
c
c       approximate internal energy (atomic + H_2)
c
        tkvn=8.6171e-5*tt
        eh2=-4.476/tkvn+(1.8031e-3+(-2.*5.0739e-7+(3.*8.1424e-11-
     *       4.*5.0501e-15*tt)*tt)*tt)*tt
c        energ=energ+eh2*anmol(2,id)
c
           j=2
           if(iprint.gt.3.and.mod(id,10).eq.1)
     *     write(80,680) id,j,anmol(j,id),anx,antt,pfmol(nelemi,id),
     *                   entmol(j,id),entt,eh2,energ
        if(iprint.gt.3) then
        write(73,673) id,tt,an,ane,anato(1,id),anion(1,id),
     *   entato(1,id),ention(1,id),entel,entt
        write(75,673) id,tt,an,ane,enp,eh2*anmol(2,id),energ
  673   format(i4,0pf8.1,1p9e11.3)
        end if
c
        jm=2*nmetal
        anhm(id)=anden(1+jm)
        anh2(id)=anden(2+jm)
C
        denso(id)=dens(id)
      eleco(id)=elec(id)
      wmmo(id)=wmm(id)
c       dens(id)=tmass*hmass
c       elec(id)=pe*tk
C       wmm(id)=dens(id)/(an-elec(id))
c       ytot(id)=(an-elec(id))/(anato(1,id)+anion(1,id))
c       wmy(id)=wmm(id)*ytot(id)/hmass
C
        RETURN
        END
C
C
C *************************************************************************
C
C     
      SUBROUTINE RUSSEL(TEM,PG)
c     =========================
c
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'MODELQ.FOR'
      COMMON/COMFH1/C(600,5),PPMOL(600),APMLOG(600),
     *              XIP(100),CCOMP(100),UIIDUI(100),P(100),
     *              FP(100),XKP(100),EPS,SWITER,
     *              NELEM(5,600),NATO(5,600),MMAX(600),
     *              NELEMX(50),NMETAL,NMOLEC,NIMAX
        DIMENSION FX(100),DFX(100),Z(100),PREV(100),WA(50)
C
        ECONST=4.342945E-1
        EPSDIE=5.0E-3
        T=5040.0/TEM
        PGLOG=log10(PG)
C
C    HEH=helium/hydrogen ratio by number
C
        HEH=CCOMP(2)/CCOMP(1)
C
C    evaluation of log XKP(MOL)
C
        DO J=1,NMOLEC
           APLOGJ=C(J,5)
           DO K=1,4
              KM5=5-K
              APLOGJ=APLOGJ*T + C(J,KM5)
           END DO
           APMLOG(J)=APLOGJ
        END DO
        DHH=(((0.1196952E-02*T-0.2125713E-01)*T+0.1545253E+00)*T
     *     -0.5161452E+01)*T+0.1277356E+02
        DHH=EXP(DHH/ECONST)
C
C  evaluation of the ionization constants
C
        TEM25=TEM**2*SQRT(TEM)
        DO I=1,NMETAL
           NELEMI = NELEMX(I)
*
* calculation of the partition functions following Irwin (1981)
C
           call mpartf(nelemi,1,0,tem,g0,dulog)
           call mpartf(nelemi,2,0,tem,g1,dulog)
           uiidui(nelemi)=g1/g0*0.6665
c        
           XKP(NELEMI)=UIIDUI(NELEMI)*TEM25*
     *                 EXP(-XIP(NELEMI)*T/ECONST)
        END DO
C
C   preliminary value of PH at high temperatures
C
        HKP=XKP(1)
        IF(T.LT.0.6) THEN
           PPH=SQRT(HKP*(PG/(1.0+HEH)+HKP))-HKP
           PH=PPH**2/HKP
         ELSE
           IF(PG/DHH.LE.0.1) THEN
              PH=PG/(1.0+HEH)
            ELSE
              PH=0.5 * (SQRT(DHH*(DHH+4.0 *PG/(1.0+HEH)))-DHH)
           END IF
        END IF
C
C  evaluation of the fictitious pressures of hydrogen
C     PG=PH+PHH+2.0*PPH+HEH*(PH+2.0*PHH+PPH)
C
        U=(1.0+2.0*HEH)/DHH
        Q=1.0+HEH
        R=(2.0+HEH)*SQRT(HKP)
        S=-1.0*PG
        X=SQRT(PH)
C
C       Russell iterations
C
        ITERAT=0
   10   CONTINUE
        F=((U*X**2+Q)*X+R)*X+S
        DF=2.0*(2.0*U*X**2+Q)*X+R
        XR=X-F/DF
C
        IF(ABS((X-XR)/XR).GT.EPSDIE) THEN
           ITERAT=ITERAT+1
           IF(ITERAT.GT.50) THEN
              WRITE(6,710) TEM,PG,X,XR,PH
  710    FORMAT(1H1, ' NOT CONVERGE IN RUSSEL '/// 'TEM=',F9.2,5X,'PG=',
     *        E12.5,5X,'X1=',E12.5,5X,'X2=',E12.5,5X,'PH=',E12.5/////)
            ELSE
              X=XR
              GO TO 10
           END IF 
        END IF
        PH=XR**2
        PHH=PH**2/DHH
        PPH=SQRT(HKP*PH)
        FPH=PH+2.0*PHH+PPH
        P(100)=PPH
C
C   evaluation of the fictitious pressure of each element
C
        DO I=1,NMETAL
           NELEMI=NELEMX(I)
           FP(NELEMI)=CCOMP(NELEMI)*FPH
        END DO
C
C   check of initialization
C
        PE=P(99)
        IF(PH.GT.P(1)) THEN
           DO I=1,NMETAL
              NELEMI=NELEMX(I)
              P(NELEMI)=FP(NELEMI)*EXP(-5.0*T/ECONST)
              IF (P(NELEMI).LT.1.E-20) P(NELEMI)=1.E-20
           END DO
           P(1)=PH
        END IF
C
C    Russell equations
C
        NITERR = 0
   20   CONTINUE
        DO I=1,NMETAL
           NELEMI=NELEMX(I)
           FX(NELEMI)=-FP(NELEMI)+P(NELEMI)*(1.0+XKP(NELEMI)/PE)
           DFX(NELEMI)=1.0+XKP(NELEMI)/PE
        END DO
C
        SPNION=0.0
        DO J=1,NMOLEC
           MMAXJ=MMAX(J)
           PMOLJL=-APMLOG(J)
           DO M=1,MMAXJ
              NELEMJ=NELEM(M,J)
              NATOMJ=NATO(M,J)
              PMOLJL=PMOLJL+DFLOAT(NATOMJ)*log10(P(NELEMJ))
           END DO
C
           IF(PMOLJL.GT.(PGLOG+1.0)) THEN
              DO M=1,MMAXJ
                 NELEMJ=NELEM(M,J)
                 NATOMJ=NATO(M,J)
                 P(NELEMJ)=1.0E-2*P(NELEMJ)
                 PMOLJL=PMOLJL+DFLOAT(NATOMJ)*(-2.0)
              END DO
           END IF
           PMOLJ=EXP(PMOLJL/ECONST)
           DO M=1,MMAXJ
              NELEMJ=NELEM(M,J)
              NATOMJ=NATO(M,J)
              ATOMJ=DFLOAT(NATOMJ)
              IF(NELEMJ.EQ.99) SPNION=SPNION+PMOLJ
              DO I=1,NMETAL
                 NELEMI=NELEMX(I)
                 IF(NELEMJ.EQ.NELEMI) THEN
                    FX(NELEMI)=FX(NELEMI)+ATOMJ*PMOLJ
                    DFX(NELEMI)=DFX(NELEMI)+ATOMJ**2*
     *                          PMOLJ/P(NELEMI)
                 END IF
              END DO
           END DO
           PPMOL(J)=PMOLJ
        END DO
C
C   solution of the Russell equations by Newton-Raphson method
C
        DO I=1,NMETAL
           NELEMI=NELEMX(I)
           WA(I)=log10(P(NELEMI)+1.0D-20)
        END DO
        IMAXP1=NMETAL+1
        WA(IMAXP1)=log10(PE+1.0D-20)                  
        DELTRS = 0.0 
        DO I=1,NMETAL
           NELEMI=NELEMX(I)
           PREV(NELEMI)=P(NELEMI)-FX(NELEMI)/DFX(NELEMI)
           PREV(NELEMI)=ABS(PREV(NELEMI))
           IF(PREV(NELEMI).LT.1.0E-20) PREV(NELEMI)=1.0E-20
           Z(NELEMI)=PREV(NELEMI)/P(NELEMI)
           DELTRS=DELTRS+ABS(Z(NELEMI)-1.0)
           IF(SWITER.GT.0.0) THEN
              P(NELEMI)=(PREV(NELEMI)+P(NELEMI))*0.5
            ELSE
              P(NELEMI)=PREV(NELEMI)
           END IF
        END DO
C
C   ionization equilibrium
C
        PEREV =0.0
        DO I=1,NMETAL
           NELEMI = NELEMX(I)
           PEREV=PEREV+XKP(NELEMI)*P(NELEMI)
        END DO
C
        PEREV=SQRT(PEREV/(1.0+SPNION/PE))
        DELTRS=DELTRS+ABS((PE-PEREV)/PE)
        PE=(PEREV+PE)*0.5
        P(99)=PE
        IF(DELTRS.GT.EPS) THEN
           NITERR=NITERR+1
           IF(NITERR.LE.NIMAX) THEN
              GO TO 20
            ELSE
              WRITE(6,605) NIMAX
           END IF
        END IF
  605   FORMAT(1H0,'*DOES NOT CONVERGE AFTER ',I4,' ITERATIONS')
C
        RETURN
        END

C
C
C ********************************************************************
C
C
c 
      subroutine mpartf(jatom,ion,indmol,t,u,dulog)
c     =============================================
c
c     yields partition functions with polynomial data from
c     ref. Irwin, A.W., 1981, ApJ Suppl. 45, 621.
c     ln u(temp)=sum(a(i)*(ln(temp))**(i-1)) 1<=a<=6
c
c     Input:
c       jatom = element number in periodic table
c       ion   = 1 for neutral, 2 for once ionized and 3 for twice ionized
c       indmol= index of a molecular specie (Tsuji index)
c       temp  = temperature
c     Output:
c       u     = partf.(linear scale) for iat,ion, or indmol, and temperature t
c       dulog = d ln(u)/d ln(T)
c
c
      implicit real*8 (a-h,o-z)
      real*8 a(6,3,92),aa(6),am(6,300)
      dimension indtsu(66),irw(300),igle(28)
      save iread,a,am
      DATA IGLE/2,1,2,1,6,9,4,9,6,1,2,1,6,9,4,9,6,1,
     *          10,21,28,25,6,25,28,21,10,21/
C
      data indtsu / 2,  5, 12, 4, 8, 7, 6,
     *              9, 11, 10, 29, 50, 59, 46, 132, 52, 19,
     *             13, 42, 38, 39, 37, 44, 36, 14, 118, 33,
     *              3, 16, 57, 32, 49, 60, 54, 41, 107,  0,
     *            148, 152, 153, 155, 0, 17, 24, 25, 28, 51,
     *            112, 119,   0,   0,21, 15, 43, 56,  0, 64,
     *             47,  65,   0,  61, 0, 62,118, 40, 66/
      data iread /0/
c
c     read data if first call:
c
      if(iread.ne.1) then
        open(67,file= 'irwin.dat',status='old')
        read(67,*)
        read(67,*)
        do j=1,92
          do i=1,3
            if(j.eq.1.and.i.eq.3) goto 10
            sp=float(j)+float(i-1)/100.
            read(67,*) spec,aa
            do k=1,6
                a(k,i,j)=aa(k)
            end do
   10    continue
         end do
       end do 
c
       read(67,*)
       read(67,*)
       read(67,*)
       do i=1,500
          irw(i)=0
       end do
       do i=1,66 
          read(67,*) spec,aa
        indm=indtsu(i)
          if(indm.gt.0) then
           irw(indm)=i
           do j=1,6
              am(j,indm)=aa(j)
           end do
          end if
      end do
        close(67)
        iread=1
      endif
c
c     evaluation of the partition function
c     stop if T is out of limits of Irwin's tables
c
        u=1.
        dulog=0.
        if(t.lt.1000.) then
          stop 'partf; temp<1000 K'
        else if(t.gt.16000.) then
c          stop 'partf; temp>16000 K'
c          write(6,601) t
c  601     format(' warning! T = ',f12.1,  'larger than 16000.'/)
           if(indmol.eq.0) then
              if(jatom.le.28.and.ion.le.jatom) u=igle(jatom-ion+1)
           end if
           return
        endif
        tl=log(t)
        u=0.
c
c     atomic species
c
      if(jatom.gt.0.and.ion.gt.0) then
        ulog=    a(1,ion,jatom)+
     *       tl*(a(2,ion,jatom)+
     *       tl*(a(3,ion,jatom)+
     *       tl*(a(4,ion,jatom)+
     *       tl*(a(5,ion,jatom)+
     *       tl*(a(6,ion,jatom))))))
        u=exp(ulog)
        dulog=   a(2,ion,jatom)+
     *       tl*(a(3,ion,jatom)*2.+
     *       tl*(a(4,ion,jatom)*3.+
     *       tl*(a(5,ion,jatom)*4.+
     *       tl*(a(6,ion,jatom)*5.))))
      end if
c
c     molecular species
c
      if(indmol.gt.0) then
        indm=indmol
        if(irw(indm).gt.0) then
           ulog=    am(1,indm)+
     *          tl*(am(2,indm)+
     *          tl*(am(3,indm)+
     *          tl*(am(4,indm)+
     *          tl*(am(5,indm)+
     *          tl*(am(6,indm))))))
           u=exp(ulog)
        dulog=   am(2,indm)+
     *       tl*(am(3,indm)*2.+
     *       tl*(am(4,indm)*3.+
     *       tl*(am(5,indm)*4.+
     *       tl*(am(6,indm)*5.))))
        end if
      end if
      return
      end
C      
C
C *************************************************************************
C
C
      subroutine tiopf(t,pf)
c     ======================
c
c     TiO partition function (data from Kurucz web site)
c
      INCLUDE 'IMPLIC.FOR'
      dimension pf0(800)
      data pf0/       
     *    29.107,    55.425,    82.417,   111.190,   142.564,   176.916,
     *   214.340,   254.774,   298.065,   344.021,   392.431,   443.089,
     *   495.795,   550.365,   606.632,   664.449,   723.686,   784.230,
     *   845.981,   908.862,   972.800,  1037.739,  1103.636,  1170.451,
     *  1238.155,  1306.723,  1376.144,  1446.403,  1517.492,  1589.409,
     *  1662.152,  1735.724,  1810.122,  1885.352,  1961.428,  2038.351,
     *  2116.119,  2194.758,  2274.260,  2354.633,  2435.907,  2518.063,
     *  2601.125,  2685.096,  2769.992,  2855.809,  2942.560,  3030.257,
     *  3118.897,  3208.496,  3299.067,  3390.598,  3483.106,  3576.598,
     *  3671.095,  3766.569,  3863.048,  3960.522,  4059.035,  4158.545,
     *  4259.074,  4360.642,  4463.259,  4566.905,  4671.582,  4777.321,
     *  4884.105,  4991.937,  5100.852,  5210.813,  5321.838,  5433.972,
     *  5547.154,  5661.417,  5776.789,  5893.211,  6010.774,  6129.422,
     *  6249.173,  6370.026,  6491.973,  6615.042,  6739.240,  6864.542,
     *  6990.959,  7118.533,  7247.214,  7377.053,  7508.012,  7640.121,
     *  7773.370,  7907.764,  8043.309,  8180.032,  8317.835,  8456.861,
     *  8597.055,  8738.396,  8880.926,  9024.672,  9169.570,  9315.610,
     *  9462.927,  9611.339,  9760.963,  9911.798, 10063.900, 10217.148,
     * 10371.572, 10527.253, 10684.109, 10842.173, 11001.469, 11161.970,
     * 11323.751, 11486.758, 11650.978, 11816.415, 11983.159, 12151.134,
     * 12320.243, 12490.668, 12662.333, 12835.234, 13009.470, 13184.926,
     * 13361.601, 13539.660, 13718.891, 13899.456, 14081.252, 14264.326,
     * 14448.643, 14634.341, 14821.225, 15009.476, 15199.021, 15389.829,
     * 15581.955, 15775.377, 15970.188, 16166.239, 16363.513, 16562.006,
     * 16761.930, 16963.301, 17165.906, 17369.881, 17575.236, 17781.814,
     * 17989.816, 18198.996, 18409.707, 18621.680, 18835.068, 19049.715,
     * 19265.768, 19483.375, 19702.006, 19922.209, 20143.668, 20366.555,
     * 20590.742, 20816.402, 21043.338, 21271.672, 21501.369, 21732.563,
     * 21965.119, 22199.068, 22434.432, 22671.266, 22909.307, 23148.898,
     * 23389.893, 23632.322, 23875.969, 24121.160, 24367.707, 24615.848,
     * 24865.471, 25116.320, 25368.604, 25622.342, 25877.512, 26134.055,
     * 26392.404, 26651.764, 26912.826, 27175.250, 27439.197, 27704.539,
     * 27971.287, 28239.572, 28509.373, 28780.707, 29053.516, 29327.602,
     * 29603.338, 29880.539, 30159.105, 30439.322, 30721.055, 31004.254,
     * 31288.818, 31575.061, 31862.693, 32151.781, 32442.586, 32734.619,
     * 33027.777, 33323.023, 33619.535, 33917.707, 34217.711, 34518.996,
     * 34821.676, 35126.195, 35432.141, 35739.602, 36048.926, 36359.488,
     * 36672.023, 36985.633, 37300.863, 37617.965, 37936.469, 38256.309,
     * 38578.074, 38901.668, 39226.461, 39552.969, 39880.852, 40210.785,
     * 40541.852, 40874.691, 41209.359, 41545.535, 41883.602, 42222.715,
     * 42563.895, 42906.508, 43250.656, 43596.902, 43944.355, 44293.695,
     * 44644.504, 44997.621, 45351.590, 45707.242, 46065.008, 46424.367,
     * 46785.605, 47148.023, 47512.496, 47878.418, 48246.426, 48615.895,
     * 48987.336, 49360.082, 49734.758, 50111.004, 50489.383, 50868.996,
     * 51250.250, 51633.691, 52018.945, 52405.715, 52794.090, 53184.340,
     * 53576.375, 53970.605, 54366.176, 54763.148, 55162.430, 55563.215,
     * 55966.391, 56371.000, 56777.176, 57185.570, 57596.074, 58007.617,
     * 58421.418, 58837.172, 59254.539, 59673.418, 60094.066, 60517.410,
     * 60941.844, 61368.660, 61797.395, 62227.590, 62659.789, 63094.238,
     * 63529.695, 63967.488, 64407.887, 64849.496, 65292.867, 65735.922,
     * 66182.000, 66631.266, 67082.055, 67534.391, 67988.992, 68446.117,
     * 68904.789, 69365.180, 69827.914, 70292.781, 70759.352, 71228.500,
     * 71699.375, 72171.672, 72647.086, 73123.984, 73603.023, 74083.516,
     * 74566.359, 75050.555, 75537.758, 76027.258, 76518.125, 77012.008,
     * 77507.063, 78003.813, 78503.977, 79006.125, 79509.320, 80015.375,
     * 80522.461, 81031.938, 81544.164, 82058.313, 82574.352, 83093.914,
     * 83614.367, 84136.820, 84662.211, 85188.867, 85719.375, 86249.977,
     * 86783.781, 87319.219, 87857.180, 88396.797, 88939.805, 89484.266,
     * 90032.023, 90580.930, 91132.563, 91686.148, 92242.742, 92799.406,
     * 93360.016, 93923.453, 94488.313, 95055.211, 95625.297, 96197.477,
     * 96771.531, 97348.156, 97926.922, 98507.453, 99091.563, 99677.938,
     *100267.234,100856.438,101449.828,102045.750,102643.094,103244.117,
     *103846.969,104450.313,105057.641,105667.188,106279.516,106894.937,
     *107512.789,108133.117,108754.758,109377.687,110005.039,110634.602,
     *111266.141,111902.133,112537.984,113178.891,113819.766,114464.312,
     *115110.969,115760.687,116412.469,117068.055,117724.547,118384.383,
     *119047.469,119712.469,120380.187,121051.336,121724.102,122399.250,
     *123076.266,123756.977,124441.195,125126.406,125816.453,126506.766,
     *127202.367,127899.086,128598.266,129299.969,130004.969,130712.016,
     *131409.266,132117.719,132828.969,133544.016,134262.750,134986.344,
     *135712.891,136439.937,137170.969,137905.562,138641.578,139380.266,
     *140122.937,140868.641,141615.484,142366.703,143123.078,143880.000,
     *144638.484,145401.594,146168.125,146935.359,147707.484,148482.641,
     *149256.578,150037.281,150821.953,151606.750,152396.094,153188.766,
     *153983.391,154782.141,155582.203,156387.234,157192.719,158003.156,
     *158815.125,159632.437,160450.766,161274.750,162098.172,162926.000,
     *163756.609,164593.141,165430.859,166270.937,167114.750,167960.797,
     *168811.562,169663.906,170517.203,171376.531,172239.469,173105.891,
     *173975.250,174847.203,175721.453,176597.250,177480.984,178366.094,
     *179253.828,180145.734,181038.000,181936.031,182837.969,183739.922,
     *184645.937,185558.281,186470.844,187387.422,188307.234,189232.281,
     *190156.000,191088.234,192022.062,192957.250,193899.328,194842.984,
     *195788.391,196736.156,197687.828,198645.719,199603.422,200569.234,
     *201536.437,202508.641,203481.000,204459.016,205438.750,206424.312,
     *207409.953,208398.734,209393.234,210391.047,211390.984,212395.516,
     *213401.547,214420.141,215431.812,216453.453,217476.734,218501.266,
     *219530.219,220560.719,221597.891,222637.875,223677.750,224725.500,
     *225777.406,226829.297,227893.125,228954.547,230020.969,231086.453,
     *232157.469,233233.047,234315.406,235395.625,236480.953,237572.125,
     *238666.484,239765.125,240863.281,241969.750,243079.250,244191.719,
     *245304.812,246427.937,247548.234,248673.562,249804.984,250942.781,
     *252078.953,253222.812,254369.641,255519.359,256671.406,257827.906,
     *258988.859,260154.734,261322.281,262458.781,263606.437,264770.625,
     *265947.750,267125.156,268314.125,269507.687,270702.344,271905.156,
     *273110.156,274318.937,275531.687,276751.344,277970.781,279198.531,
     *280425.750,281663.250,282897.469,284138.906,285383.594,286637.031,
     *287891.156,289147.625,290413.312,291678.719,292946.031,294225.875,
     *295501.344,296782.656,298070.094,299363.875,300652.250,301953.750,
     *303260.062,304563.781,305874.375,307191.437,308517.031,309835.750,
     *311159.375,312490.937,313827.469,315166.781,316511.031,317860.406,
     *319214.969,320565.875,321929.344,323296.906,324660.219,326035.687,
     *327413.844,328794.406,330173.156,331566.156,332953.469,334356.187,
     *335757.625,337165.562,338566.094,339984.750,341402.937,342828.125,
     *344257.562,345686.750,347123.125,348564.250,350008.906,351453.219,
     *352908.062,354361.469,355828.000,357292.500,358765.719,360233.687,
     *361713.562,363200.187,364685.656,366174.500,367673.594,369174.906,
     *370678.969,372191.125,373708.937,375225.281,376743.719,378270.406,
     *379804.500,381334.250,382879.125,384420.812,385969.531,387519.812,
     *389078.937,390639.781,392213.875,393782.437,395359.156,396943.625,
     *398527.625,400110.937,401711.750,403310.344,404908.937,406513.875,
     *408125.781,409741.906,411356.875,412979.500,414613.125,416245.500,
     *417889.094,419530.000,421179.906,422831.531,424484.344,426153.187,
     *427816.406,429489.094,431161.312,432840.656,434517.000,436215.281,
     *437896.000,439602.594,441300.625,443016.156,444722.906,446445.437,
     *448164.812,449885.937,451615.094,453351.594,455090.125,456833.281,
     *458582.719,460335.344,462094.844,463857.094,465629.906,467402.781,
     *469178.406,470963.750,472745.906,474539.594,476333.312,478131.125,
     *479934.000,481740.750,483557.844,485376.625,487202.937,489033.562,
     *490868.031,492709.281,494547.375,496401.094,498249.594,500110.250,
     *501966.594,503836.062,505704.437,507580.687,509469.187,511349.781,
     *513239.000,515137.187,517038.812,518942.906,520858.156,522767.094,
     *524610.625,526433.812,528331.062,530253.437,532185.500,534127.875,
     *536073.937,538028.312,539983.375,541954.687,543916.312,545902.500,
     *547874.812,549857.125,551850.937,553839.937,555836.625,557838.500,
     *559849.937,561859.375,563880.625,565889.875,567916.000,569953.625,
     *571990.375,574034.937,576085.062,578127.375,580188.937,582251.000,
     *584328.812,586385.562,588464.062,590551.875,592644.625,594722.250,
     *596829.937,598931.375,601029.687,603142.812,605262.812,607384.625,
     *609513.125,611644.000,613775.875,615930.375,618073.750,620218.437,
     *622381.937,624524.312,626697.500,628869.000,631040.937,633223.562,
     *635409.187,637597.562,639800.187,642002.125,644212.562,646416.250,
     *648633.562,650864.187,653083.687,655315.312,657549.687,659795.500,
     *662032.250,664292.875,666542.312,668806.250,671071.312,673340.937,
     *675626.938,677898.750/
c
      it=int(t/10.)
      if(it.gt.800) it=800
      pf=pf0(it)
      return
      end

C
C
C ********************************************************************
C
C
