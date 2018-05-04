C
C     Parameters that specify dimensions of arrays
C
      PARAMETER (MATOM  =   99, ! max.num. of explicit atoms
     *           MION   =  170, ! max.num. of explicit ions
     *           MLEVEL = 1134, ! max.num. of explicit levels
     *           MLVEXP =  233, ! max.num. of explicit linearized levels
     *           MTRANS =21000, ! max.num. of all transitions
     *           MDEPTH =   70, ! max.num. of depth points
     *           MFREQ =125000, ! max.num. of frequency points
     *           MFREQP=220000, ! working arrays of frequency
     *           MFREQC =15000, ! max.num. of freq.points in continuum
     *           MFREX  =   34, ! max.num. of linearized frequencies
     *           MFREQL =15798, ! max.num. of frequencies per line
     *           MTOT   =  270, ! max.num. of linearized parameters
     *           MMU    =    3, ! max.num. of angle points
     *           MFIT   =  157, ! max.num. of fit points (OP b-f c.s)
     *           MITJ   =  280, ! max.num. of overlapping transitions
     *           MMCDW  =    6, ! max.num. of levels with pseudocont.
     *           MMER   =    2, ! max.num. of merged levels
     *           MVOIGT =   80, ! max.num. of lines with Voigt profile
     *           MZZ    =   10, ! maximum charge for occup.prob. ions
     *           NLMX   =   80, ! highest hydrogenic level considered
     *           MSMX   =    1, ! size of matrix kept in memory in SOLVE
     *           MFREQ1 =    1, ! =1 for ISPLIN<5; =MFREQ otherwise
     *           MFRTAB =30000, ! max.num. of freqeuncies in opac.table
c    *           MFRTAB =    1, ! max.num. of freqeuncies in opac.table
     *           MTABT  =   50, ! max.number of temps in opac.table
     *           MTABR  =   50, ! max.number of densities in opac.table
     *           MCROSS=MLEVEL+5,! max.num. of b-f cross.secs.
     *           MBF =  MLEVEL, ! max.num. of b-f transitions
     *           MDEPTC =    1, ! max.num. of depth points (Compton)
     *           MMUC   =    2, ! max.num. of angle points (Compton)
     *           MLEVE3 =    1, ! =1 for diag.prec; =MLEVEL for trid.prec.
     *           MLVEX3 =    1, ! =1 for diag.oper.; =MLVEXP for tridiag.   .
     *           MTRAN3 =    1) ! =1 for diag.oper; =MTRANS for tridiag.
C
C     Basic physical constants
C
      PARAMETER (H     = 6.6256D-27,      ! Planck constant     h 
     *           BOLK  = 1.38054D-16,     ! Boltzmann constant  k
     *           HK    = 4.79928144D-11,  ! h/k
     *           CAS   = 2.997925D18,     ! light speed c (A/s)
     *           EH    = 2.17853041D-11,  ! ionizaton energy of hydrogen
     *           BN    = 1.4743D-2,       ! 2*h/c**3, c -light speed
     *           SIGE  = 6.6516D-25,      ! Thomson scattering c-s
     *           SIG4P = 4.5114062D-6,    ! Stefan-Boltzmann const/4pi
     *           PI4H  = 1.8966D27,       ! 4pi/h
     *           PCK   = 4.19168946D-10,  ! 4pi/c
     *           HMASS = 1.67333D-24)     ! mass of hydrogen atom
C
C     Basic mathematical constants
C
      PARAMETER (UN    = 1.0D0,
     *           HALF  = 0.5D0,
     *           TWO   = 2.0D0)
C
C     Unit number
C
      PARAMETER (IBUFF = 95)
C
C     Basic parameters
C
      COMMON/BASNUM/NATOM,NION,NLEVEL,NTRANS,ND,NFREQ,NFREQC,NFREQE,
     *              IOPTAB,IDISK,IZSCAL,IDMFIX,IHESO6,IFMOL,IFENTR,
     *              NFREQL,NLEV0,ICOLHN,IOSCOR,ILGDER,IFRYB,IFRSET,
     *              NFREAD,NELSC,NTRANC,IOVER,JALI,IBC,IUBC,
     *              IRDER,ILMCOR,IFDIEL,IFALIH,IFTENE,ITNDRE,
     *              ILPSCT,ILASCT,IRTE,IDLTE,IBFINT,INTRPL,ICHANG,
     *              NATOMS,IPSLTE,ISPODF,ITLUCY,NRETC,IFRAYL,IFPRAD
      COMMON/INPPAR/TEFF,GRAV,
     *              YTOT(MDEPTH),WMM(MDEPTH),WMY(MDEPTH),
     *              xmstar,xmdot,rstar,alpha0,reynum,
     *              QGRAV,EDISC,DZETA,RELDST,
     *              visc,zeta0,zeta1,dmvisc,fractv,
     *              omeg32,wbarm,wbar,alphav,pgas0,
     *              bergfc,cutlym,cutbal,
     *              ISPLIN,IRSPLT,ivisc,ibche,LTE,LTGREY,LCHC,LRESC
      COMMON/MATKEY/NN,NN0,INHE,INRE,INPC,INSE,INZD,INMP,NDRE,insel
      COMMON/FIXDEN/IFIXDE
      COMMON/INVINT/XI2(NLMX),XI3(NLMX)
      COMMON/RUNKEY/CHMAX,ITER,NITER,NITZER,INIT,LAC2,LFIN
      COMMON/CONKEY/HMIX0,crflim,
     *              NCONIT,ICONV,INDL,IPRESS,ITEMP,ICBEG,
     *              itmcor,iconre,ideepc,ndcgap,IDCONZ
      COMMON/OPCKEY/NCON,IRSCT,IOPHMI,IOPH2P,IOPHE1,IOPHE2,
     *              IOPHL1,IOPHL2,IPHE2C,IFMOFF
      COMMON/PRINTS/IPRINT,IPRING,IPRIND,IPRINP,ICOOLP,ICHCKP,
     *              IPOPAC,IPRINI
      COMMON/PSILIM/DPSILG,DPSILT,DPSILN,DPSILD
      COMMON/CENTRL/ZND,IFZ0
c
      COMMON/ANGLES/AMU(MMU),WTMU(MMU),FMU(MMU),NMU
c
      common/comptn/amuc(mmuc),wtmuc(mmuc),
     *              amuc1(mmuc),amuc2(mmuc),amuc3(mmuc),
     *              amuj(mmuc),amuk(mmuc),amuh(mmuc),amun(mmuc),
     *              calph(mmuc,mmuc),cbeta(mmuc,mmuc),
     *              cgamm(mmuc,mmuc),RADZER,FRLCOM,
     *              SIGEC(MFREQ),ijorig(mfreq)
      common/angnum/nmuc
      common/compti/nedd,nsti,islab,ilbc,icompt,icomst,icomde,
     *              icombc,icmdra,knish,itcomp,icomve,icomrt,
     *              ichcoo,icomgr
      common/comite/ncfor1,ncfor2,nccoup,ncitot,ncfull
      common/mlcons/aconml,bconml,cconml
      common/taursl/taurs(mdepth)
