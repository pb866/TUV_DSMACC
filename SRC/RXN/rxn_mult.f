*= This file contains the following subroutines, related to reading/loading
*= the product (cross section) x (quantum yield) for photo-reactions of
*= organic compounds with multiple chromphores in MCM-GECKO,
*= which where not yet present in TUV5.2:
*=
*=     mm01 through mm05
*=
*=
*= Further routines concern the calculation of generic j values, where
*= quantum yields of carbonyl compounds (including ketenes) have been
*= set to 1 and previous cross sections (MCM preferences) are used.
*= These j values are used in the MCM to synthesise total j values,
*= where absorption cross sections are added up and the total
*= quantum yield is set to 1:
*=
*=    mm06 through mm12

*=============================================================================*

      SUBROUTINE mm01(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! NMEK, nitroxymethyl ethyl ketone

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for                =*
*=  nitroxymethyl ethyl ketone photolysis:                                   =*
*=        C2H5COCH2NO3 + hv -> C2H5COCH2O + NO2                              =*
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
!      REAL qy
      REAL sig
      INTEGER iw

      INTEGER mabs, xsvers(3)
      CHARACTER(llog) :: logmsg(nlog)
      CHARACTER(lspc) :: spc


***************** C2H5COCH2NO3 photolysis *************************

* Setting photolysis index and cross section options
      j = j+1
      jlabel(j) = 'C2H5COCH2NO3 -> C2H5COCH2O + NO2'

* cross section options
      spc = "C2H5COCH2NO3"
      xsvers = (/2, 2, 2/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "from least square fit by Barnes et al"
      logmsg(2) = "from data by Barnes et al"
      CALL set_option(2,spc,"abs",logmsg,xsvers,mabs)


* cross sections

      IF(mabs==2) THEN
        OPEN(UNIT=kin,FILE='DATAJ1/MCMext/MULT/NMEK_Bar93.abs',
     $       STATUS='old')
        do i = 1, 4
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
!        qy = 1.0


* combine:

      DO iw = 1, nw - 1
        IF(mabs==1) THEN
          sig = exp(-1.011e-3*wc(iw)**2+0.5437*wc(iw)-116.9)
         ELSEIF(mabs==2) THEN
          sig = yg(iw)
        ENDIF

        DO i = 1, nz
          sq(j,i,iw) = sig! * qy
        ENDDO

      ENDDO

      END

* ============================================================================*

      SUBROUTINE mm02(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! M1NEK, methyl 1-nitroxyethyl ketone

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for                =*
*=  methyl 1-nitroxyethyl ketone photolysis:                                 =*
*=        CH3COCH(NO3)CH3 + hv -> CH3COCH(O.)CH3 + NO2                       =*
*=                                                                           =*
*=  Cross section:  Barnes et al. (1993)                                     =*
*=  Quantum yield:  see options below                                        =*
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
      REAL qy
      REAL sig
      INTEGER iw

      INTEGER mabs, myld, xsvers(3), qyvers(3)
      CHARACTER(llog) :: logmsg(nlog)
      CHARACTER(lspc) :: spc


***************** CH3COCH(NO3)CH3 photolysis *************************


* Setting photolysis index and cross section/quantum yield options
      j = j+1
      jlabel(j) = 'CH3COCH(NO3)CH3 + hv -> CH3COCH(O.)CH3 + NO2'

* cross section options
      spc = "CH3COCH(NO3)CH3"
      xsvers = (/2, 2, 2/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "from least square fit by Barnes et al"
      logmsg(2) = "from data by Barnes et al"
      CALL set_option(2,spc,"abs",logmsg,xsvers,mabs)

* quantum yield options
      qyvers = (/2, 2, 2/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "estimated unity (Calvert et al. 2011)"
      logmsg(2) = "estimated 0.75 (Muller et al. 2014)"
      CALL set_option(2,spc,"yld",logmsg,qyvers,myld)


* cross sections

      IF(mabs==2) THEN
        OPEN(UNIT=kin,FILE='DATAJ1/MCMext/MULT/M1NEK_Bar93.abs',
     $       STATUS='old')
        do i = 1, 4
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
        qy = 1.0
       ELSEIF(myld==2) THEN
        qy = 0.75
      ENDIF


* combine:

      DO iw = 1, nw - 1
        IF(mabs==1) THEN
          sig = exp(-1.044e-3*wc(iw)**2+0.578*wc(iw)-123.5)
         ELSEIF(mabs==2) THEN
          sig = yg(iw)
        ENDIF

        DO i = 1, nz
          sq(j,i,iw) = sig * qy
        ENDDO

      ENDDO

      END

* ============================================================================*

      SUBROUTINE mm03(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! 2-oxo-cyclohexyl nitrate

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for                =*
*=  2-oxo-cyclohexyl nitrate photolysis:                                     =*
*=        C6H9-2-=O-1-NO3 + hv -> RO. + NO2                                  =*
*=                                                                           =*
*=  Cross section:  Wängberg et al. (1996)                                   =*
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

      REAL yg(kw), dum
!     REAL qy
      REAL sig
      INTEGER iw


      INTEGER mabs, xsvers(3)
      CHARACTER(llog) :: logmsg(nlog)
      CHARACTER(lspc) :: spc


***************** 2-oxo-cyclohexyl nitrate photolysis *************************


* Setting photolysis index and cross section options

      j = j+1
      jlabel(j) = 'O=C1H[CH2]4C1HNO3 -> RO. + NO2'

* cross section options
      spc = "O=C1H[CH2]4C1HNO3"
      xsvers = (/1, 1, 2/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "from Wangberg et al"
      logmsg(2) = "as plotted in Calvert et al. 2011 book"
      CALL set_option(2,spc,"abs",logmsg,xsvers,mabs)


* cross sections

      IF(mabs==1) THEN
        OPEN(UNIT=kin,FILE='DATAJ1/MCMext/MULT/Wangberg96.abs',
     $       STATUS='old')
        do i = 1, 11
          read(kin,*)
        enddo

        n  = 14
        DO i = 1, n
          READ(kin,*) x1(i), dum, y1(i)
          y1(i) = y1(i)*1.e-20
        ENDDO
        CLOSE(kin)

       ELSEIF(mabs==2) THEN

        OPEN(UNIT=kin,
     $       FILE  ='DATAJ1/MCMext/MULT/2oxo-cHexNit_calv.abs',
     $       STATUS='old')
        do i = 1, 8
          read(kin,*)
        enddo

        n  = 94
        DO i = 1, n
          READ(kin,*) x1(i), y1(i)
          y1(i) = y1(i)*1.e-20
        ENDDO
        CLOSE(kin)
      ENDIF

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)


* quantum yields
!     qy = 1.0


* combine:

      DO iw = 1, nw - 1
        sig = yg(iw)
        DO i = 1, nz
          sq(j,i,iw) = sig! * qy
        ENDDO
      ENDDO

      END

* ============================================================================*

      SUBROUTINE mm04(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! 2-hexanone-5-hydroperoxide

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for                =*
*=  2-hexanone-5-hydroperoxide photolysis:                                   =*
*=        CH3COCH2CH2CH(OOH)CH3 + hv -> RO. + OH                             =*
*=                                                                           =*
*=  Cross section:  Jorand et al. (2000)                                     =*
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
      PARAMETER(kdata=100)

      INTEGER i, n
      REAL x1(kdata)
      REAL y1(kdata)

* local

      REAL yg(kw)
!     REAL qy
      REAL sig
      INTEGER iw


***************** 2-hexanone-5-hydroperoxide photolysis *************************

      j = j+1
      jlabel(j) = 'CH3COCH2CH2CH(OOH)CH3 + hv -> RO. + OH'


*cross sections

      OPEN(UNIT=kin,FILE='DATAJ1/MCMext/MULT/2oxo5oohHex.abs',
     $     STATUS='old')
      do i = 1, 6
        read(kin,*)
      enddo

      n  = 81
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
      ENDDO
      CLOSE(kin)

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)


* quantum yields
!      qy = 1.0


* combine

      DO iw = 1, nw - 1
        sig = yg(iw)
        DO i = 1, nz
          sq(j,i,iw) = sig! * qy
        ENDDO
      ENDDO

      END

* ============================================================================*

      SUBROUTINE mm05(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! oxohexyl-hydroperoxide

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for                =*
*=  oxohexyl-hydroperoxide photolysis:                                       =*
*=        oxohexyl-hydroperoxide + hv -> RO. + OH                            =*
*=                                                                           =*
*=  Cross section:  Jorand et al. (2000)                                     =*
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
      PARAMETER(kdata=100)

      INTEGER i, n
      REAL x1(kdata)
      REAL y1(kdata)

* local

      REAL yg(kw)
!     REAL qy
      REAL sig
      INTEGER iw


***************** oxohexyl-hydroperoxide photolysis *************************

      j = j+1
      jlabel(j) = 'oxohexyl-hydroperoxide + hv -> RO. + OH'


      OPEN(UNIT=kin,FILE='DATAJ1/MCMext/MULT/oxo+oohHex.abs',
     $     STATUS='old')
      do i = 1, 7
        read(kin,*)
      enddo

      n  = 81
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

      SUBROUTINE mm06(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! C2ALDqy1

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Generation of generic product (cross section) x (quantum yield) for      =*
*=  acetaldehyde with single channel with QY = 1. MCM preference used.       =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*lcl, string identifier for each photolysis         (O)=*
*=           reaction defined                                                =*
*=  lcl    - INTEGER, length of character for labels                         =*
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
      PARAMETER(kdata=150)

      INTEGER i, n
      REAL x1(kdata)
      REAL y1(kdata)

* local

      REAL yg(kw)
      INTEGER  iw


***************** C2ALDqy1 photolysis *************************

      j = j+1
      jlabel(j) = 'C2ALDqy1'

      OPEN(UNIT=kin,FILE='DATAJ1/CH3CHO/CH3CHO_iup.abs',STATUS='old')
        do i = 1, 4
          read(kin,*)
        enddo
      n = 106
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
        y1(i) = y1(i) * 1.e-20
      ENDDO
      CLOSE(kin)

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)


* Map onto grid:
      DO iw = 1, nw - 1
         DO i = 1, nz
            sq(j,i,iw) = yg(iw)
         ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE mm07(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! C3ALDqy1

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Generation of generic product (cross section) x (quantum yield) for      =*
*=  propionaldehyde with single channel with QY = 1. MCM preference used.    =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*lcl, string identifier for each photolysis         (O)=*
*=           reaction defined                                                =*
*=  lcl    - INTEGER, length of character for labels                         =*
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
      PARAMETER(kdata=150)

      INTEGER i, n, n1
      REAL x1(kdata), x1oh(kdata)
      REAL y1(kdata), y1oh(kdata)

* local

      REAL yg(kw), ygoh(kw)
      INTEGER iw


***************** C3ALDqy1 photolysis *************************

      j = j+1
      jlabel(j) = 'C3ALDqy1'
      j = j+1
      jlabel(j) = 'C3ALDqy1oh'

      OPEN(UNIT=kin,FILE='DATAJ1/C2H5CHO/C2H5CHO_iup.abs',
     $     STATUS='old')
      do i = 1, 4
        read(kin,*)
      enddo
      n = 106
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
        y1(i) = y1(i) * 1.e-20
      ENDDO
      CLOSE(kin)

      x1oh(:) = 0.
      x1oh(:) = x1(:) - 10.
      y1oh(:) = y1(:)
      n1 = n
      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j-1),0.,0.,yg)
      n = n1
      CALL interpol(x1oh,y1oh,kdata,n,nw,wl,jlabel(j),0.,0.,ygoh)

* Map onto grid:
      DO iw = 1, nw - 1
         DO i = 1, nz
           sq(j-1,i,iw) = yg(iw)
           sq(j  ,i,iw) = ygoh(iw)
         ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE mm08(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! C4ALDqy1

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Generation of generic product (cross section) x (quantum yield) for      =*
*=  butyraldehyde with single channel with QY = 1. MCM preference used.      =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*lcl, string identifier for each photolysis         (O)=*
*=           reaction defined                                                =*
*=  lcl    - INTEGER, length of character for labels                         =*
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

      INTEGER i, n, n1
      REAL x1(kdata), x1oh(kdata)
      REAL y1(kdata), y1oh(kdata)

* local

      REAL yg(kw),ygoh(kw)
      REAL sig,sigoh
      INTEGER iw


***************** C4ALDqy1 photolysis *************************

      j = j+1
      jlabel(j) = 'C4ALDqy1'
      j = j+1
      jlabel(j) = 'C4ALDqy1oh'


      OPEN(UNIT=kin,FILE='DATAJ1/MCMext/ALD/nC3H7CHO_iup.abs',
     $     STATUS='old')
      do i = 1, 5
         read(kin,*)
      enddo

      n = 106
      DO i = 1, n
         READ(kin,*) x1(i), y1(i)
         y1(i) = y1(i)*1.e-20
      ENDDO
      CLOSE(kin)

      x1oh(:) = 0.
      x1oh(:) = x1(:) - 10.
      y1oh(:) = y1(:)
      n1 = n
      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j-1),0.,0.,yg)
      n = n1
      CALL interpol(x1oh,y1oh,kdata,n,nw,wl,jlabel(j),0.,0.,ygoh)


* quantum yields: no pressure dependence

* Map onto grid:
      DO iw = 1, nw - 1
         DO i = 1, nz
           sq(j-1,i,iw) = yg(iw)
           sq(j  ,i,iw) = ygoh(iw)
         ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE mm09(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! C5ALDqy1

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Generation of generic product (cross section) x (quantum yield) for      =*
*=  valeraldehyde with single channel with QY = 1. MCM preference used.      =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*lcl, string identifier for each photolysis         (O)=*
*=           reaction defined                                                =*
*=  lcl    - INTEGER, length of character for labels                         =*
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

      INTEGER i, n, n1
      REAL x1(kdata), x1oh(kdata)
      REAL y1(kdata), y1oh(kdata)


* local

      REAL yg(kw), ygoh(kw), dum
      INTEGER iw


***************** C5ALDqy1 photolysis *************************

      j = j+1
      jlabel(j) = 'C5ALDqy1'
      j = j+1
      jlabel(j) = 'C5ALDqy1oh'


      OPEN(UNIT=kin,FILE='DATAJ1/MCMext/ALD/Tadic.abs',
     $     STATUS='old')

      do i = 1, 5
        read(kin,*)
      enddo

      n = 121
      DO i = 1, n
        READ(kin,*) x1(i), dum, y1(i)
      ENDDO
      CLOSE(kin)

      x1oh(:) = 0.
      x1oh(:) = x1(:) - 10.
      y1oh(:) = y1(:)
      n1 = n
      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j-1),0.,0.,yg)
      n = n1
      CALL interpol(x1oh,y1oh,kdata,n,nw,wl,jlabel(j),0.,0.,ygoh)


* quantum yields: no pressure dependence

* Map onto grid:
      DO iw = 1, nw - 1
        DO i = 1, nz
          sq(j-1,i,iw) = yg(iw)
          sq(j  ,i,iw) = ygoh(iw)
        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE mm10(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! C6ALDqy1

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Generation of generic product (cross section) x (quantum yield) for      =*
*=  capronaldehyde with single channel with QY = 1. MCM preference used.     =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*lcl, string identifier for each photolysis         (O)=*
*=           reaction defined                                                =*
*=  lcl    - INTEGER, length of character for labels                         =*
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

      INTEGER i, n, n1
      REAL x1(kdata), x1oh(kdata)
      REAL y1(kdata), y1oh(kdata)

* local

      REAL yg(kw), ygoh(kw), dum
      INTEGER iw


***************** C6ALDqy1 photolysis *************************

      j = j+1
      jlabel(j) = 'C6ALDqy1'
      j = j+1
      jlabel(j) = 'C6ALDqy1oh'

      OPEN(UNIT=kin,FILE='DATAJ1/MCMext/ALD/hexanal_rad',
     $     STATUS='old')
      do i = 1, 2
        read(kin,*)
      enddo

      n = 89
      DO i = 1, n
        READ(kin,*) x1(i), y1(i), dum
      ENDDO
      CLOSE(kin)

      x1oh(:) = 0.
      x1oh(:) = x1(:) - 10.
      y1oh(:) = y1(:)
      n1 = n
      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j-1),0.,0.,yg)
      n = n1
      CALL interpol(x1oh,y1oh,kdata,n,nw,wl,jlabel(j),0.,0.,ygoh)


* quantum yields: no pressure dependence

* Map onto grid:
      DO iw = 1, nw - 1
        DO i = 1, nz
          sq(j-1,i,iw) = yg(iw)
          sq(j  ,i,iw) = ygoh(iw)
        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE mm11(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! C7ALDqy1

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Generation of generic product (cross section) x (quantum yield) for      =*
*=  heptyl aldehyde with single channel with QY = 1. MCM preference used.    =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*lcl, string identifier for each photolysis         (O)=*
*=           reaction defined                                                =*
*=  lcl    - INTEGER, length of character for labels                         =*
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

      INTEGER i, n, n1
      REAL x1(kdata), x1oh(kdata)
      REAL y1(kdata), y1oh(kdata)

* local

      REAL yg(kw), ygoh(kw)
      INTEGER iw


***************** C7ALDqy1 photolysis *************************

      j = j+1
      jlabel(j) = 'C7ALDqy1'
      j = j+1
      jlabel(j) = 'C7ALDqy1oh'

      OPEN(UNIT=kin,FILE='DATAJ1/MCMext/ALD/nHeptanal.abs',
     $     STATUS='old')
      do i = 1, 5
        read(kin,*)
      enddo

      n = 11
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
      ENDDO
      CLOSE(kin)

      x1oh(:) = 0.
      x1oh(:) = x1(:) - 10.
      y1oh(:) = y1(:)
      n1 = n
      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j-1),0.,0.,yg)
      n = n1
      CALL interpol(x1oh,y1oh,kdata,n,nw,wl,jlabel(j),0.,0.,ygoh)


* quantum yields: no pressure dependence

* Map onto grid:
      DO iw = 1, nw - 1
        DO i = 1, nz
          sq(j-1,i,iw) = yg(iw)
          sq(j  ,i,iw) = ygoh(iw)
        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE mm12(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! C8ALDqy1

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Generation of generic product (cross section) x (quantum yield) for      =*
*=  octyl aldehyde with single channel with QY = 1. MCM preference used.     =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*lcl, string identifier for each photolysis         (O)=*
*=           reaction defined                                                =*
*=  lcl    - INTEGER, length of character for labels                         =*
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

      INTEGER i, n, n1
      REAL x1(kdata), x1oh(kdata), y1(kdata), y1oh(kdata)

* local

      REAL yg(kw), ygoh(kw)
      REAL sig,sigoh
      INTEGER iw


***************** C8ALDqy1 photolysis *************************

      j = j+1
      jlabel(j) = 'C8ALDqy1'
      j = j+1
      jlabel(j) = 'C8ALDqy1oh'


* cross sections

      OPEN(UNIT=kin,FILE='DATAJ1/MCMext/ALD/nAldehydes.abs',
     $     STATUS='old')

      do i = 1, 9
        read(kin,*)
      enddo

      n = 121
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
        y1(i) = y1(i)*1.e-20
      ENDDO
      CLOSE(kin)

      x1oh(:) = 0.
      x1oh(:) = x1(:) - 10.
      y1oh(:) = y1(:)
      n1 = n

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j-7),0.,0.,yg)
      n = n1
      CALL interpol(x1oh,y1oh,kdata,n,nw,wl,jlabel(j-4),0.,0.,ygoh)


* quantum yields: no pressure dependence

* Map onto grid:
      DO iw = 1, nw - 1
         DO i = 1, nz
             sq(j-1,i,iw) = sig
             sq(j  ,i,iw) = sigOH
         ENDDO
      ENDDO

      END

*=============================================================================*
