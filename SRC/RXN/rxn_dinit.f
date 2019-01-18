*= This file contains the following subroutines, related to reading/loading
*= the product (cross section) x (quantum yield) for photo-reactions of
*= organic dinitrates in MCM/GECKO-A, which where not yet present in TUV5.2:
*=
*=     dn01 through dn07

*=============================================================================*

      SUBROUTINE dn01(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! isopropylene dinitrate

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for                =*
*=  CH3CH(NO3)CH2NO3 photolysis:                                             =*
*=        CH3CH(NO3)CH2NO3 + hv -> products                                  =*
*=                                                                           =*
*=  Cross section:  Barnes et al. (1993)                                     =*
*=  Quantum yield:  assumed unity with equal branching                       =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE '../../params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)

      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER(lcl) jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:

      INTEGER j

* data arrays

      INTEGER kdata
      PARAMETER(kdata=580)

      INTEGER i, n
      REAL x1(kdata)
      REAL y1(kdata)

* local

      REAL yg(kw)
      REAL qy1, qy2
      REAL sig
      INTEGER iw

      INTEGER mabs, myld, xsvers(3), qyvers(3)
      CHARACTER(llog) :: logmsg(nlog)
      CHARACTER(lspc) :: spc


***************** isopropylene dinitrate photolysis *************************

* Setting photolysis index and cross section/quantum yield options
      j = j+1
      jlabel(j) = 'CH3CH(NO3)CH2NO3 -> CH3CH(NO3)CH2O + NO2'
      j = j+1
      jlabel(j) = 'CH3CH(NO3)CH2NO3 -> CH3CH(O)CH2NO3 + NO2'

* cross section options
      spc = "CH3CH(NO3)CH2NO3"
      xsvers = (/1, 2, 2/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "from least square fit by Barnes et al"
      logmsg(2) = "from data by Barnes et al"
      CALL set_option(2,spc,"abs",logmsg,xsvers,mabs)

* quantum yield options
      qyvers = (/1, 2, 1/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "estimated unity with equal branching"
      logmsg(2) = "scaled externally"
      CALL set_option(2,spc,"yld",logmsg,qyvers,myld)


* cross sections

      IF(mabs==2) THEN
        OPEN(UNIT=kin,FILE='DATAJ1/MCMext/DINIT/12PropNit.abs',
     $       STATUS='old')
        do i = 1, 6
          read(kin,*)
        enddo

        n  = 20
        DO i = 1, n
          READ(kin,*) x1(i), y1(i)
        ENDDO
        CLOSE(kin)

        CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)
      ENDIF


* quantum yields
      IF(myld==1) THEN
        qy1 = 0.5
        qy2 = 0.5
      ELSEIF(myld==2) THEN
        qy1 = 1.0
        qy2 = 1.0
      ENDIF


* combine:

      DO iw = 1, nw - 1
        sig = exp(-5.99e-4*wc(iw)**2+0.2915*wc(iw)-79.24)
        IF(mabs==2 .AND. wc(iw)>=245. .AND. wc(iw)<=340.)THEN
          sig = yg(iw)
        ENDIF
        DO i = 1, nz
          sq(j-1,i,iw) = sig * qy2
          sq(j  ,i,iw) = sig * qy1
        ENDDO
      ENDDO

      END

* ============================================================================*

      SUBROUTINE dn02(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! 1,2-butylene dinitrate

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for                =*
*=  CH3CH2CH(NO3)CH2NO3 photolysis:                                          =*
*=        CH3CH2CH(NO3)CH2NO3 + hv -> products                               =*
*=                                                                           =*
*=  Cross section:  Barnes et al. (1993)                                     =*
*=  Quantum yield:  assumed unity with equal branching                       =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE '../../params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)

      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER(lcl) jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:

      INTEGER j

* data arrays

      INTEGER kdata
      PARAMETER(kdata=580)

      INTEGER i, n
      REAL x1(kdata)
      REAL y1(kdata)

* local

      REAL yg(kw)
      REAL qy1, qy2
      REAL sig
      INTEGER iw

      INTEGER mabs, myld, xsvers(3), qyvers(3)
      CHARACTER(llog) :: logmsg(nlog)
      CHARACTER(lspc) :: spc


***************** 1,2-butylene dinitrate photolysis *************************

* Setting photolysis index and cross section/quantum yield options
      j = j+1
      jlabel(j) = 'C2H5CH(NO3)CH2NO3 -> C2H5CH(NO3)CH2O + NO2'
      j = j+1
      jlabel(j) = 'C2H5CH(NO3)CH2NO3 -> C2H5CH(O)CH2NO3 + NO2'

* cross section options
      spc = "C2H5CH(NO3)CH2NO3"
      xsvers = (/1, 2, 2/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "from least square fit by Barnes et al"
      logmsg(2) = "from data by Barnes et al"
      CALL set_option(2,spc,"abs",logmsg,xsvers,mabs)

* quantum yield options
      qyvers = (/1, 2, 1/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "estimated unity with equal branching"
      logmsg(2) = "scaled externally"
      CALL set_option(2,spc,"yld",logmsg,qyvers,myld)


* cross sections

      IF(mabs==2) THEN
        OPEN(UNIT=kin,FILE='DATAJ1/MCMext/DINIT/12ButNit.abs',
     $       STATUS='old')
        do i = 1, 6
          read(kin,*)
        enddo

        n  = 20
        DO i = 1, n
          READ(kin,*) x1(i), y1(i)
        ENDDO
        CLOSE(kin)

        CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)
      ENDIF


* quantum yields
      IF(myld==1) THEN
        qy1 = 0.5
        qy2 = 0.5
      ELSEIF(myld==2) THEN
        qy1 = 1.0
        qy2 = 1.0
      ENDIF


* combine:

      DO iw = 1, nw - 1
        sig = exp(-6.217e-4*wc(iw)**2+0.3025*wc(iw)-80.41)
        IF(mabs==2 .AND. wc(iw)>=245. .AND. wc(iw)<=340.)THEN
          sig = yg(iw)
        ENDIF
        DO i = 1, nz
          sq(j-1,i,iw) = sig * qy2
          sq(j  ,i,iw) = sig * qy1
        ENDDO
      ENDDO

      END

* ============================================================================*

      SUBROUTINE dn03(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! 2,3-butylene dinitrate

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for                =*
*=  CH3CH(NO3)CH(NO3)CH3 photolysis:                                         =*
*=        CH3CH(NO3)CH(NO3)CH3 + hv -> CH3CH(NO3)CH(O.)CH3 + NO2             =*
*=                                                                           =*
*=  Cross section:  Barnes et al. (1993)                                     =*
*=  Quantum yield:  assumed unity                                            =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE '../../params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)

      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER(lcl) jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:

      INTEGER j

* data arrays

      INTEGER kdata
      PARAMETER(kdata=580)

      INTEGER i, n
      REAL x1(kdata)
      REAL y1(kdata)

* local

      REAL yg(kw)
!     REAL qy
      REAL sig
      INTEGER iw

      INTEGER mabs, xsvers(3)
      CHARACTER(llog) :: logmsg(nlog)
      CHARACTER(lspc) :: spc


***************** 2,3-butylene dinitrate photolysis *************************

* Setting photolysis index and cross section options
      j = j+1
      jlabel(j) = 'CH3CH(NO3)CH(NO3)CH3 -> RO. + NO2'

* cross section options
      spc = "CH3CH(NO3)CH(NO3)CH3"
      xsvers = (/1, 2, 2/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "from least square fit by Barnes et al"
      logmsg(2) = "from data by Barnes et al"
      CALL set_option(2,spc,"abs",logmsg,xsvers,mabs)


* cross sections

      IF(mabs==2) THEN
        OPEN(UNIT=kin,FILE='DATAJ1/MCMext/DINIT/23ButNit.abs',
     $       STATUS='old')
        do i = 1, 6
          read(kin,*)
        enddo

        n  = 20
        DO i = 1, n
          READ(kin,*) x1(i), y1(i)
        ENDDO
        CLOSE(kin)

        CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)
      ENDIF


* quantum yields
!      qy = 1.0


* combine:

      DO iw = 1, nw - 1
        sig = exp(-5.74e-4*wc(iw)**2+0.2771*wc(iw)-77.47)
        IF(mabs==2 .AND. wc(iw)>=245. .AND. wc(iw)<=340.)THEN
          sig = yg(iw)
        ENDIF
        DO i = 1, nz
          sq(j  ,i,iw) = sig! * qy
        ENDDO
      ENDDO

      END

* ============================================================================*

      SUBROUTINE dn04(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! 1,4-dinitrooxy-2-butene

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for                =*
*=  CH3CH(NO3)CH(NO3)CH3 photolysis:                                         =*
*=        CH2(NO3)CH=CHCH2NO3 + hv -> CH2(NO3)CH=CHCH2O + NO2                =*
*=                                                                           =*
*=  Cross section:  Barnes et al. (1993)                                     =*
*=  Quantum yield:  assumed unity                                            =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE '../../params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)

      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER(lcl) jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:

      INTEGER j

* data arrays

      INTEGER kdata
      PARAMETER(kdata=580)

      INTEGER i, n
      REAL x1(kdata)
      REAL y1(kdata)

* local

      REAL yg(kw)
!     REAL qy
      REAL sig
      INTEGER iw

      INTEGER mabs, xsvers(3)
      CHARACTER(llog) :: logmsg(nlog)
      CHARACTER(lspc) :: spc


***************** 1,4-dinitrooxy-2-butene photolysis *************************

* Setting photolysis index and cross section options
      j = j+1
      jlabel(j) = 'CH2(NO3)CH=CHCH2NO3 -> RO. + NO2'

* cross section options
      spc = "CH3CH(NO3)CH(NO3)CH3"
      xsvers = (/1, 2, 2/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "from least square fit by Barnes et al"
      logmsg(2) = "from data by Barnes et al"
      CALL set_option(2,spc,"abs",logmsg,xsvers,mabs)


* cross sections

      IF(mabs==2) THEN
        OPEN(UNIT=kin,FILE='DATAJ1/MCMext/DINIT/14Nit2Butene.abs',
     $       STATUS='old')
        do i = 1, 6
          read(kin,*)
        enddo

        n  = 18
        DO i = 1, n
          READ(kin,*) x1(i), y1(i)
        ENDDO
        CLOSE(kin)

        CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)
      ENDIF


* quantum yields

!      qy = 1.0


* combine:

      DO iw = 1, nw - 1
        sig = exp(-5.432e-4*wc(iw)**2+0.2631*wc(iw)-75.92)
        IF(mabs==2 .AND. wc(iw)>=245. .AND. wc(iw)<=330.)THEN
          sig = yg(iw)
        ENDIF
        DO i = 1, nz
          sq(j  ,i,iw) = sig! * qy
        ENDDO
      ENDDO

      END

* ============================================================================*

      SUBROUTINE dn05(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! 3,4-dinitrooxy-1-butene

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for                =*
*=  CH3CH(NO3)CH(NO3)CH3 photolysis:                                         =*
*=        CH2=CHCH(NO3)CH2NO3 + hv -> CH2=CHCH(NO3)CH2O + NO2                =*
*=                                                                           =*
*=  Cross section:  Barnes et al. (1993)                                     =*
*=  Quantum yield:  assumed unity with equal branching                       =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE '../../params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)

      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER(lcl) jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:

      INTEGER j

* data arrays

      INTEGER kdata
      PARAMETER(kdata=580)

      INTEGER i, n
      REAL x1(kdata)
      REAL y1(kdata)

* local

      REAL yg(kw)
      REAL qy1, qy2
      REAL sig
      INTEGER iw

      INTEGER mabs, myld, xsvers(3), qyvers(3)
      CHARACTER(llog) :: logmsg(nlog)
      CHARACTER(lspc) :: spc


***************** 3,4-dinitrooxy-1-butene photolysis *************************

* Setting photolysis index and cross section/quantum yield options
      j = j+1
      jlabel(j) = 'CH2=CHCH(NO3)CH2NO3 -> CH2=CHCH(NO3)CH2O + NO2'
      j = j+1
      jlabel(j) = 'CH2=CHCH(NO3)CH2NO3 -> CH2=CHCH(O)CH2NO3 + NO2'

* cross section options
      spc = "CH2=CHCH(NO3)CH2NO3"
      xsvers = (/1, 2, 2/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "from least square fit by Barnes et al"
      logmsg(2) = "from data by Barnes et al"
      CALL set_option(2,spc,"abs",logmsg,xsvers,mabs)

* quantum yield options
      qyvers = (/1, 2, 1/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "estimated unity with equal branching"
      logmsg(2) = "scaled externally"
      CALL set_option(2,spc,"yld",logmsg,qyvers,myld)


* cross sections

      IF(mabs==2) THEN
        OPEN(UNIT=kin,FILE='DATAJ1/MCMext/DINIT/34Nit1Butene.abs',
     $       STATUS='old')
        do i = 1, 6
          read(kin,*)
        enddo

        n  = 20
        DO i = 1, n
          READ(kin,*) x1(i), y1(i)
        ENDDO
        CLOSE(kin)

        CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)
      ENDIF


* quantum yields
      IF(myld==1) THEN
        qy1 = 0.5
        qy2 = 0.5
      ELSEIF(myld==2) THEN
        qy1 = 1.0
        qy2 = 1.0
      ENDIF


* combine:

      DO iw = 1, nw - 1
        sig = exp(-6.217e-4*wc(iw)**2+0.3025*wc(iw)-80.41)
        IF(mabs==2 .AND. wc(iw)>=245. .AND. wc(iw)<=340.)THEN
          sig = yg(iw)
        ENDIF
        DO i = 1, nz
          sq(j-1,i,iw) = sig * qy2
          sq(j  ,i,iw) = sig * qy1
        ENDDO
      ENDDO

      END

* ============================================================================*

      SUBROUTINE dn06(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! 1-methyl-cyclohexyl-1,2-dinitrate

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for                =*
*=  1-methyl-cyclohexyl-1,2-dinitrate photolysis:                            =*
*=        C6H9-1-CH3-1,2-NO3 + hv -> products                                =*
*=                                                                           =*
*=  Cross section:  Wängberg et al. (1996)                                   =*
*=  Quantum yield:  assumed unity with equal branching                       =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE '../../params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)

      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER(lcl) jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:

      INTEGER j

* data arrays

      INTEGER kdata
      PARAMETER(kdata=580)

      INTEGER i, n
      REAL x1(kdata)
      REAL y1(kdata)

* local

      REAL yg(kw), dum
      REAL qy1, qy2
      REAL sig
      INTEGER iw

      INTEGER mabs, myld, xsvers(3), qyvers(3)
      CHARACTER(llog) :: logmsg(nlog)
      CHARACTER(lspc) :: spc


***************** 1-methyl-cyclohexyl-1,2-dinitrate photolysis *************************

* Setting photolysis index and cross section/quantum yield options
      j = j+1
      jlabel(j) = 'CH3C1(NO3)[CH2]4C1H(NO3) -> R1O + NO2'
      j = j+1
      jlabel(j) = 'CH3C1(NO3)[CH2]4C1H(NO3) -> R2O + NO2'

* cross section options
      spc = "CH3C1(NO3)[CH2]4C1H(NO3)"
      xsvers = (/1, 2, 2/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "from Wangberg et al"
      logmsg(2) = "as plotted in the Calvert et al. 2011 book"
      CALL set_option(2,spc,"abs",logmsg,xsvers,mabs)

* quantum yield options
      qyvers = (/1, 2, 1/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "estimated unity with equal branching"
      logmsg(2) = "scaled externally"
      CALL set_option(2,spc,"yld",logmsg,qyvers,myld)


* cross sections

      IF(mabs==1) THEN
        OPEN(UNIT=kin,FILE='DATAJ1/MCMext/POLY/Wangberg96.abs',
     $       STATUS='old')
        do i = 1, 9
          read(kin,*)
        enddo

        n  = 14
        DO i = 1, n
          READ(kin,*) x1(i), dum, dum, y1(i)
          y1(i) = y1(i)*1.e-20
        ENDDO
        CLOSE(kin)

       ELSEIF(mabs==2) THEN

        OPEN(UNIT=kin,
     $       FILE  ='DATAJ1/MCMext/DINIT/1MecHex12DiNit_calv.abs',
     $       STATUS='old')
        do i = 1, 8
          read(kin,*)
        enddo

        n  = 111
        DO i = 1, n
          READ(kin,*) x1(i), y1(i)
          y1(i) = y1(i)*1.e-20
        ENDDO
        CLOSE(kin)
      ENDIF

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)


* quantum yields
      IF(myld==1) THEN
        qy1 = 0.5
        qy2 = 0.5
      ELSEIF(myld==2) THEN
        qy1 = 1.0
        qy2 = 1.0
      ENDIF


* combine:

      DO iw = 1, nw - 1
        sig = yg(iw)
        DO i = 1, nz
          sq(j-1,i,iw) = sig * qy2
          sq(j  ,i,iw) = sig * qy1
        ENDDO
      ENDDO

      END

* ============================================================================*

      SUBROUTINE dn07(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! uDINIT

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for                =*
*=  generic unsaturated dinitrate photolysis:                                =*
*=        RNO3 + hv -> RO + NO2                                              =*
*=                                                                           =*
*=  Cross section:  Average of 1,4-dinitroxy-2-butene and                    =*
*=                  3,4-dinitroxy-1-butene (Barnes et al., 1993)             =*
*=  Quantum yield:  assumed unity                                            =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE '../../params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)

      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER(lcl) jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:

      INTEGER j

* data arrays

      INTEGER kdata
      PARAMETER(kdata=580)

      INTEGER i, n
      REAL x1(kdata)
      REAL y1(kdata)

* local

      REAL yg(kw)
!      REAL qy
      REAL sig
      INTEGER iw


***************** uDINIT photolysis *************************

* Setting photolysis index
      j = j+1
      jlabel(j) = 'uDINIT -> RO. + NO2'


      OPEN(UNIT=kin,FILE='DATAJ1/MCMext/DINIT/uDINIT.abs',
     $     STATUS='old')
      do i = 1, 7
        read(kin,*)
      enddo

      n  = 20
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
      ENDDO
      CLOSE(kin)

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)


* quantum yields
!      qy = 1.0


* combine:

      DO iw = 1, nw - 1
        sig = yg(iw)
        DO i = 1, nz
          sq(j,i,iw) = sig! * qy
        ENDDO
      ENDDO

      END

* ============================================================================*
