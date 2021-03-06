*= This file contains the following subroutines, related to reading/loading
*= the product (cross section) x (quantum yield) for photo-reactions of
*= organic nitrates in MCM/GECKO-A, which were not yet present in TUV5.2:
*=
*=     mn01 through mn14

*=============================================================================*

      SUBROUTINE mn01(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! 1-C5H11ONO2

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for n-C5H11ONO2    =*
*=  1-pentyl nitrate photolysis:                                             =*
*=           n-C5H11ONO2 + hv -> n-C5H11O + NO2                              =*
*=                                                                           =*
*=  Absorption cross sections: options for the following:                    =*
*=  (from MPI-Mainz Spectral Atlas)                                          =*
*=                                                                           =*
*=  K.C. Clemitshaw, J. Williams, O.V. Rattigan, D.E. Shallcross, K.S. Law,  =*
*=  and R.A. Cox, "Gas-phase ultraviolet absorption cross-sections and       =*
*=  atmospheric lifetimes of several C2-C5 alkyl nitrates,"                  =*
*=  J. Photochem. Photobiol. A: Chem. 102, 117-126 (1996).                   =*
*=                                                                           =*
*=  L. Zhu and D. Kellis, "Temperature dependence of the UV absorption cross =*
*=  sections and photodissociation products of C3 - C5 alkyl nitrates,"      =*
*=  Chem. Phys. Lett. 278, 41-48 (1997).                                     =*
*=                                                                           =*
*=  Quantum yield assumed unity                                              =*
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
      PARAMETER(kdata=200)

      INTEGER i, n
      REAL x1(kdata), x2(kdata), xs2(kdata)
      REAL y1(kdata), y2(kdata), ys1(kdata), ys2(kdata)

* local

      REAL yg1(kw), yg2(kw)
      REAL sig ! qy
      INTEGER iw
      INTEGER mabs, xsvers(3)
      CHARACTER(llog) :: logmsg(nlog)
      CHARACTER(lspc) :: spc

***************** n-C5H11ONO2 photolysis *************************

      j = j+1
      jlabel(j) = 'n-C5H11ONO2 -> n-C5H11O + NO2'

* cross section options
      spc = "n-C5H11ONO2"
      xsvers = (/2, 3, 2/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "from Clemitshaw et al. 1996"
      logmsg(2) = "from Zhu and Kellis 1997"
      logmsg(3) = "mean of Zhu and Kellis/Clemitshaw"
      CALL set_option(3,spc,"abs",logmsg,xsvers,mabs)


* cross sections

      IF( mabs.EQ.1 .OR. mabs.EQ.3) THEN
         OPEN(UNIT=kin,FILE='DATAJ1/MCMext/NIT/1PentNit_Cle96.abs',
     $        STATUS='old')
         DO i = 1, 6
            READ(kin,*)
         ENDDO
         n = 15
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
         ENDDO
         CLOSE(kin)
         IF(mabs==3) THEN
           ys1(:) = y1(:)
          ELSE
           CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg1)
         ENDIF
      ENDIF

      IF( mabs.EQ.2 .OR. mabs.EQ.3) THEN
         OPEN(UNIT=kin,FILE='DATAJ1/MCMext/NIT/1PentNit_Z+K97.abs',
     $        STATUS='old')
         DO i = 1, 6
            READ(kin,*)
         ENDDO
         n = 16
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
         ENDDO
         CLOSE(kin)
         IF(mabs==3) THEN
           xs2(:) = x1(:)
           ys2(:) = y1(:)
          ELSE
           CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg1)
         ENDIF
      ENDIF

* Calculate Average of above cross sections
      IF(mabs==3) THEN
        n = 16

        x1(1) = xs2(1)
        y1(1) = ys2(1)
        DO i = 2, n
          x1(i) = xs2(i)
          y1(i) = (ys1(i-1)+ys2(i))/2.
        ENDDO

        CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg1)
      ENDIF

! Temperature dependency (assume for all cases)
      OPEN(UNIT=kin,FILE='DATAJ1/MCMext/NIT/1PentNit_Z+K97_Tdep.abs',
     $     STATUS='old')
      DO i = 1, 9
        READ(kin,*)
      ENDDO
      n = 16
      DO i = 1, n
        READ(kin,*) x2(i), y2(i)
        y2(i) = y2(i) * 1.E-3
      ENDDO
      CLOSE(kin)

      CALL interpol(x2,y2,kdata,n,nw,wl,jlabel(j),0.,y2(n),yg2)


* quantum yield = 1

* combine
      DO iw = 1, nw - 1
         DO i = 1, nz

!        qy = 1.
!        assume T-dep. from Zhu and Kellis 1997 for all cases
         sig = yg1(iw)*exp(yg2(iw)*(tlev(i)-298.))

         sq(j,i,iw) = sig! * qy

         ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE mn02(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! 2-C5H11ONO2

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for 2-C5H11ONO2        =*
*=  photolysis:                                                              =*
*=          2-C5H11ONO2 + hv -> 2-C5H11O + NO2                               =*
*=                                                                           =*
*=  Cross section: Roberts and Fajer (1989)                                  =*
*=  Quantum yield: Assumed to be unity                                       =*
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

      INTEGER i
      INTEGER iw

* local

!      REAL qy
      REAL sig

***************** 2-C5H11ONO2 photolysis *************************

      j = j + 1
      jlabel(j) = '2-C5H11ONO2 -> 2-C5H11O + NO2'

* cross sections calculated with least square fit from Roberts and Fajer (1989)
* see also DATAJ1/MCMext/NIT/2PentNit.abs
* quantum yield = 1

!      qy = 1.

      DO iw = 1, nw - 1
         sig = exp(-1.231e-3*wc(iw)**2+0.6454*wc(iw)-129.1)
         DO i = 1, nz
            sq(j,i,iw) = sig ! * qy
         ENDDO
      ENDDO

      RETURN
      END

* ============================================================================*

      SUBROUTINE mn03(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! 3-C5H11ONO2

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for 3-C5H11ONO2        =*
*=  photolysis:                                                              =*
*=          3-C5H11ONO2 + hv -> 2-C5H11O + NO2                               =*
*=                                                                           =*
*=  Cross section: Roberts and Fajer (1989)                                  =*
*=  Quantum yield: Assumed to be unity                                       =*
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

      INTEGER i
      INTEGER iw

* local

!      REAL qy
      REAL sig

***************** 3-C5H11ONO2 photolysis *************************

      j = j + 1
      jlabel(j) = '3-C5H11ONO2 -> 3-C5H11O + NO2'

* cross sections calculated with least square fit from Roberts and Fajer (1989)
* see also DATAJ1/MCMext/NIT/3PentNit.abs
* quantum yield = 1

!      qy = 1.

      DO iw = 1, nw - 1
         sig = exp(-1.446e-3*wc(iw)**2+0.7712*wc(iw)-147.4)
         DO i = 1, nz
            sq(j,i,iw) = sig ! * qy
         ENDDO
      ENDDO

      RETURN
      END

* ============================================================================*

      SUBROUTINE mn04(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! c-C5H11ONO2

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for c-C5H11ONO2        =*
*=  photolysis:                                                              =*
*=          c-C5H11ONO2 + hv -> c-C5H11O + NO2                               =*
*=                                                                           =*
*=  Cross section: Roberts and Fajer (1989)                                  =*
*=  Quantum yield: Assumed to be unity                                       =*
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

      INTEGER i
      INTEGER iw

* local

!      REAL qy
      REAL sig

***************** c-C5H11ONO2 photolysis *************************

      j = j + 1
      jlabel(j) = 'c-C5H11ONO2 -> c-C5H11O + NO2'

* cross sections calculated with least square fit from Roberts and Fajer (1989)
* see also DATAJ1/MCMext/NIT/3PentNit.abs
* quantum yield = 1

!      qy = 1.

      DO iw = 1, nw - 1
        sig = exp(-1.884e-3*wc(iw)**2+1.0109*wc(iw)-180.5)
         DO i = 1, nz
            sq(j,i,iw) = sig ! * qy
         ENDDO
      ENDDO

      RETURN
      END

* ============================================================================*

      SUBROUTINE mn05(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! isobutyl nitrate

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for                =*
*=  isobutyl nitrate photolysis:                                             =*
*=                    (CH3)2CHCH2ONO2 + hv -> i-C4H9O + NO2                  =*
*=                                                                           =*
*=  Cross section: Clemitshaw et al. (1989)                                  =*
*=  Quantum yield: Assumed to be unity                                       =*
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
      PARAMETER(kdata=200)

      INTEGER i, n
      REAL x1(kdata), y1(kdata)

* local

      REAL yg1(kw)
      REAL sig ! qy
      INTEGER iw

***************** i-C4H9ONO2 photolysis *************************

      j = j+1
      jlabel(j) = 'i-C4H9ONO2 -> i-C4H9O + NO2'

* cross sections
      OPEN(UNIT=kin,FILE='DATAJ1/MCMext/NIT/isobutylNit.abs',
     $     STATUS='old')
      DO i = 1, 6
        READ(kin,*)
      ENDDO
      n = 25
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
      ENDDO
      CLOSE(kin)

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg1)

* quantum yield = 1

* combine
      DO iw = 1, nw - 1
         DO i = 1, nz
         sig = yg1(iw)
         sq(j,i,iw) = sig! * qy
         ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE mn06(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! isopentyl nitrate

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for                =*
*=  isopentyl nitrate photolysis:                                            =*
*=                    (CH3)2CHCH2CH2ONO2 + hv -> i-C5H11O + NO2              =*
*=                                                                           =*
*=  Cross section: Turberg et al. (1990)                                     =*
*=  Quantum yield: Assumed to be unity                                       =*
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
      PARAMETER(kdata=200)

      INTEGER i, n
      REAL x1(kdata), y1(kdata)

* local

      REAL yg1(kw)
      REAL sig ! qy
      INTEGER iw

      INTEGER mabs, xsvers(3)
      CHARACTER(llog) :: logmsg(nlog)
      CHARACTER(lspc) :: spc

***************** i-C5H11ONO2 photolysis *************************

      j = j+1
      jlabel(j) = 'i-C5H11ONO2 -> i-C5H11O + NO2'

* cross section options
      spc = "i-C5H11ONO2"
      xsvers = (/1, 2, 2/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "from experimental data by Zhy and Kronin, 1997"
      logmsg(2) = "extended by scaled n-C5H11ONO2 cross sections."
      CALL set_option(2,spc,"abs",logmsg,xsvers,mabs)

* cross sections
      IF(mabs==1) THEN
        OPEN(UNIT=kin,FILE='DATAJ1/MCMext/NIT/isopentylNit.abs',
     $       STATUS='old')
        DO i = 1, 6
          READ(kin,*)
        ENDDO
        n = 25
        DO i = 1, n
          READ(kin,*) x1(i), y1(i)
        ENDDO
        CLOSE(kin)
       ELSEIF(mabs==2) THEN
        OPEN(UNIT=kin,FILE='DATAJ1/MCMext/NIT/iNIText.abs',
     $       STATUS='old')
        DO i = 1, 7
          READ(kin,*)
        ENDDO
        n = 33
        DO i = 1, n
          READ(kin,*) x1(i), y1(i)
          y1(i) = 1e-20*y1(i)
        ENDDO
        CLOSE(kin)
      ENDIF

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg1)


* quantum yield = 1

* combine
      DO iw = 1, nw - 1
         DO i = 1, nz
           sig = yg1(iw)
           sq(j,i,iw) = sig! * qy
         ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE mn07(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! C1 OH-subst. alkyl nitrates
*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for                    =*
*=  C1 OH-subst. alkyl nitrates photolysis:                                  =*
*=          C(X)2(OH)ONO2 + hv -> RO + NO2                                   =*
*=                                                                           =*
*=  Cross section: Choices as given below in the source code                 =*
*=  Quantum yield: Assumed to be unity                                       =*
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
      PARAMETER (kdata = 2000)

      INTEGER i, n
      INTEGER iw
      REAL x1(kdata)
      REAL y1(kdata)

* local

      REAL yg(kw), yg1(kw), yg2(kw), yg3(kw)
!      REAL qy
      REAL sig, fOH


***************** C1 OH-subst. alkyl nitrates photolysis *************************

      j = j + 1
      jlabel(j) = 'C1(OH)NO3 -> C1(OH)O + NO2'

* Only MCM/GECKO-A choices considered for generic species:
* IUPAC (2002) recommendation

* cross sections

      OPEN(UNIT=kin,FILE='DATAJ1/RONO2/CH3NO3_iup.abs',
     $    STATUS='old')
      DO i = 1, 5
        READ(kin,*)
      ENDDO
      n = 21
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
      ENDDO
      CLOSE (kin)

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)

      OPEN(UNIT=kin,FILE='DATAJ1/RONO2/CH3NO3_Tdep_iup.abs',
     $     STATUS='old')
      DO i = 1, 8
        READ(kin,*)
      ENDDO
      n = 19
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
        y1(i) = y1(i) * 1.e-3
      ENDDO
      CLOSE (kin)

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg1)


* OH scaling factor

! 2-hydroxyethyl nitrate
      OPEN(UNIT=kin,FILE='DATAJ1/MCMext/NIT/NitEtOH.abs',
     $    STATUS='old')
      DO i = 1, 6
        READ(kin,*)
      ENDDO
      n = 23
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
      ENDDO
      CLOSE (kin)

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg2)

! ethyl nitrate
      OPEN(UNIT=kin,FILE='DATAJ1/RONO2/C2H5NO3_iup.abs',
     $    STATUS='old')
      DO i = 1, 5
        READ(kin,*)
      ENDDO
      n = 32
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
        y1(i) = y1(i) * 1.e-20
      ENDDO
      CLOSE (kin)

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),y1(1),0.,yg3)


* quantum yield = 1
!      qy = 1.

      DO iw = 1, nw - 1
        DO i = 1, nz
          IF(wc(iw)<270.) THEN
            fOH = 0.577
           ELSEIF(yg3(iw)>0.) THEN
            fOH = yg2(iw)/yg3(iw)
           ELSE
            fOH = 0.
          ENDIF
          sig = yg(iw) * exp (yg1(iw) * (tlev(i)-298.)) * fOH
          sq(j,i,iw) = sig! * qy
        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE mn08(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! C≥3 OH-subst. alkyl nitrates
*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for                    =*
*=  C≥3 OH-subst. alkyl nitrates photolysis:                                 =*
*=          R(OH)ONO2 + hv -> RO + NO2                                       =*
*=                                                                           =*
*=  Cross section: Choices as given below in the source code                 =*
*=  Quantum yield: Assumed to be unity                                       =*
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
      PARAMETER (kdata = 2000)

      INTEGER i, n
      INTEGER iw
      REAL x1(kdata)
      REAL y1(kdata)

* local

      REAL yg1(kw), yg2(kw), yg3(kw)
!      REAL qy
      REAL sig, fOH


***************** C≥3 OH-subst. alkyl nitrates photolysis *************************

      j = j + 1
      jlabel(j) = 'R(OH)NO3 -> R(OH)O + NO2'

* Only MCM/GECKO-A choices considered for generic species:
* IUPAC (2002) recommendation

* cross sections

      OPEN(UNIT=kin,FILE='DATAJ1/RONO2/nC3H7ONO2_iup2006.abs',
     $     STATUS='old')
      DO i = 1, 3
         READ(kin,*)
      ENDDO
      n = 32
      DO i = 1, n
         READ(kin,*) x1(i), y1(i)
         y1(i) = y1(i) * 1.e-20
      ENDDO
      CLOSE (kin)

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg1)


* OH scaling factor

! 2-hydroxyethyl nitrate
      OPEN(UNIT=kin,FILE='DATAJ1/MCMext/NIT/NitEtOH.abs',
     $    STATUS='old')
      DO i = 1, 6
        READ(kin,*)
      ENDDO
      n = 23
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
      ENDDO
      CLOSE (kin)

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),y1(1),0.,yg2)

! ethyl nitrate
      OPEN(UNIT=kin,FILE='DATAJ1/RONO2/C2H5NO3_iup.abs',
     $    STATUS='old')
      DO i = 1, 5
        READ(kin,*)
      ENDDO
      n = 32
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
        y1(i) = y1(i) * 1.e-20
      ENDDO
      CLOSE (kin)

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),y1(1),0.,yg3)


* quantum yield = 1
!      qy = 1.

      DO iw = 1, nw - 1
        IF(wc(iw)<270.) THEN
          fOH = 0.577
         ELSEIF(yg3(iw)>0.) THEN
          fOH = yg2(iw)/yg3(iw)
         ELSE
          fOH = 0.
        ENDIF
        sig = yg1(iw) * fOH
        DO i = 1, nz
          sq(j,i,iw) = sig! * qy
        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE mn09(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! OH-subst. iNIT
*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for                    =*
*=  iso-branched OH-subst. alkyl nitrates photolysis:                        =*
*=          iR(OH)ONO2 + hv -> iRO + NO2                                     =*
*=                                                                           =*
*=  Cross section: Choices as given below in the source code                 =*
*=  Quantum yield: Assumed to be unity                                       =*
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
      PARAMETER (kdata = 2000)

      INTEGER i, n
      INTEGER iw
      REAL x1(kdata)
      REAL y1(kdata)

* local

      REAL yg1(kw), yg2(kw), yg3(kw)
!      REAL qy
      REAL sig, fOH


***************** OH-subst. iNIT photolysis *************************

      j = j + 1
      jlabel(j) = 'iR(OH)NO3 -> iR(OH)O + NO2'

* Only MCM/GECKO-A choices considered for generic species:
* IUPAC (2002) recommendation

* cross sections
      OPEN(UNIT=kin,FILE='DATAJ1/MCMext/NIT/isopentylNit.abs',
     $     STATUS='old')
      DO i = 1, 6
        READ(kin,*)
      ENDDO
      n = 25
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
      ENDDO
      CLOSE(kin)

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg1)

* OH scaling factor

! 2-hydroxyethyl nitrate
      OPEN(UNIT=kin,FILE='DATAJ1/MCMext/NIT/NitEtOH.abs',
     $    STATUS='old')
      DO i = 1, 6
        READ(kin,*)
      ENDDO
      n = 23
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
      ENDDO
      CLOSE (kin)

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),y1(1),0.,yg2)

! ethyl nitrate
      OPEN(UNIT=kin,FILE='DATAJ1/RONO2/C2H5NO3_iup.abs',
     $    STATUS='old')
      DO i = 1, 5
        READ(kin,*)
      ENDDO
      n = 32
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
        y1(i) = y1(i) * 1.e-20
      ENDDO
      CLOSE (kin)

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),y1(1),0.,yg3)


* quantum yield = 1
!      qy = 1.

      DO iw = 1, nw - 1
        IF(wc(iw)<270.) THEN
          fOH = 0.577
         ELSEIF(yg3(iw)>0.) THEN
          fOH = yg2(iw)/yg3(iw)
         ELSE
          fOH = 0.
        ENDIF
        sig = yg1(iw) * fOH
        DO i = 1, nz
          sq(j,i,iw) = sig! * qy
        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE mn10(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! OH-subst. tNIT
*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for                    =*
*=  tert-branched OH-subst. alkyl nitrates photolysis:                       =*
*=          iR(OH)ONO2 + hv -> iRO + NO2                                     =*
*=                                                                           =*
*=  Cross section: Choices as given below in the source code                 =*
*=  Quantum yield: Assumed to be unity                                       =*
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
      PARAMETER (kdata = 2000)

      INTEGER i, n
      INTEGER iw
      REAL x1(kdata)
      REAL y1(kdata)

* local

      REAL yg2(kw), yg3(kw)
!      REAL qy
      REAL sig, fOH, a, b, c


***************** OH-subst. tNIT photolysis *************************

      j = j + 1
      jlabel(j) = 'tR(OH)NO3 -> tR(OH)O + NO2'

* Only MCM/GECKO-A choices considered for generic species:
* IUPAC (2002) recommendation

* cross sections
* coefficients from Roberts and Fajer 1989, over 270-330 nm

      a = -0.993E-3
      b = 0.5307
      c = -115.5


* OH scaling factor

! 2-hydroxyethyl nitrate
      OPEN(UNIT=kin,FILE='DATAJ1/MCMext/NIT/NitEtOH.abs',
     $    STATUS='old')
      DO i = 1, 6
        READ(kin,*)
      ENDDO
      n = 23
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
      ENDDO
      CLOSE (kin)

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),y1(1),0.,yg2)

! ethyl nitrate
      OPEN(UNIT=kin,FILE='DATAJ1/RONO2/C2H5NO3_iup.abs',
     $    STATUS='old')
      DO i = 1, 5
        READ(kin,*)
      ENDDO
      n = 32
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
        y1(i) = y1(i) * 1.e-20
      ENDDO
      CLOSE (kin)

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),y1(1),0.,yg3)


* quantum yield = 1
!      qy = 1.


      DO iw = 1, nw - 1
        IF(wc(iw)<270.) THEN
          fOH = 0.577
         ELSEIF(yg3(iw)>0.) THEN
          fOH = yg2(iw)/yg3(iw)
         ELSE
          fOH = 0.
        ENDIF
        sig = EXP(a*wc(iw)**2 + b*wc(iw) + c) * fOH
        DO i = 1, nz
          sq(j,i,iw) = sig! * qy
        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE mn11(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! PANs
*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for                    =*
*=  photolysis of PANs:                                                      =*
*=          iR(OH)ONO2 + hv -> iRO + NO2                                     =*
*=                                                                           =*
*=  Cross section: average of PAN and PPN                                    =*
*=  Quantum yield: estimate based on HNO4 and other PANs                     =*
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
      PARAMETER (kdata = 100)

      INTEGER i, n, n2
      INTEGER iw
      REAL x1(kdata), x2(kdata)
      REAL y1(kdata), y2(kdata)

* local

      REAL yg(kw), yg1(kw)
      REAL qyNO2, qyNO3
      REAL sig


***************** PANs photolysis *************************

      j = j+1
      jlabel(j) = 'PAN -> RCO(OO) + NO2'
      j = j+1
      jlabel(j) = 'PAN -> RCO(O) + NO3'

* Only MCM/GECKO-A choices considered for generic species:

* PAN cross sections

      OPEN(UNIT=kin,FILE='DATAJ1/RONO2/meanPAN+PPN.abs',STATUS='OLD')
      DO i = 1, 9
         READ(kin,*)
      ENDDO
      n = 77
      n2 = n
      DO i = 1, n
         READ(kin,*) x1(i), y1(i), y2(i)
         y1(i) = y1(i) * 1.E-20
         y2(i) = y2(i) * 1E-3
         x2(i) = x1(i)
      ENDDO
      CLOSE(kin)

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)
      CALL interpol(x2,y2,kdata,n2,nw,wl,jlabel(j),0.,0.,yg1)


* quantum yields
      qyNO2 = 0.61
      qyNO3 = 0.39


      DO iw = 1, nw - 1
        DO i = 1, nz
          sig = yg(iw) * EXP(yg1(iw)*(tlev(i)-298.))
          sq(j-1,i,iw) = sig * qyNO2
          sq(j  ,i,iw) = sig * qyNO3
        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE mn12(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! 2/3-C5H11ONO2

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for average of 2-      =*
*=  and 3-C5H11ONO2 photolysis:                                              =*
*=              C5H11ONO2 + hv -> C5H11O + NO2                               =*
*=                                                                           =*
*=  Cross section: Roberts and Fajer (1989)                                  =*
*=  Quantum yield: Assumed to be unity                                       =*
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

      INTEGER i
      INTEGER iw

* local

!      REAL qy
      REAL sig


***************** 3-C5H11ONO2 photolysis *************************

      j = j + 1
      jlabel(j) = 'C5H11ONO2 -> C5H11O + NO2'

* cross sections calculated with least square fit
* to data of 2- adn 3-pentyl nitrate from Roberts and Fajer (1989)
* see also DATAJ1/MCMext/NIT/avrge23PentNit.abs
* quantum yield = 1

!      qy = 1.

      DO iw = 1, nw - 1
         sig = exp(-1.339e-3*wc(iw)**2+0.7086*wc(iw)-138.3)
         DO i = 1, nz
            sq(j,i,iw) = sig ! * qy
         ENDDO
      ENDDO

      RETURN
      END

* ============================================================================*

      SUBROUTINE mn13(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! RNO2

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for RNO2           =*
*=  RNO2 photolysis:                                                         =*
*=           RNO2 + hv -> R + NO2                                            =*
*=                                                                           =*
*=  Absorption cross sections: Taylor et al. (1980) for nitromethane         =*
*=  (from MPI-Mainz Spectral Atlas)                                          =*
*=  Quantum yields:            Estimates based on Calvert et al. 2011 book   =*
*=                                                                           =*
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
      PARAMETER(kdata=200)

      INTEGER i, n
      REAL x1(kdata)
      REAL y1(kdata)

* local

      REAL yg1(kw)
      REAL sig, qy1, qy2
      INTEGER iw

      INTEGER myld, qyvers(3)
      CHARACTER(llog) :: logmsg(nlog)
      CHARACTER(lspc) :: spc


***************** RNO2 photolysis *************************

* Setting photolysis index and cross section/quantum yield options
      j = j+1
      jlabel(j) = 'CH3NO2 -> CH3 + NO2'
      j = j+1
      jlabel(j) = 'RNO2 -> alkene + HONO'

* quantum yield options
      spc = "RNO2"
      qyvers = (/1, 2, 1/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "with HONO channel of 0.08"
      logmsg(2) = "with HONO channel scaled externally"
      CALL set_option(2,spc,"yld",logmsg,qyvers,myld)


       OPEN(UNIT=kin,FILE='DATAJ1/MCMext/NIT/CH3NO2_Taylor80.abs',
     $      STATUS='old')
       DO i = 1, 8
          READ(kin,*)
       ENDDO
       n = 35
       DO i = 1, n
          READ(kin,*) x1(i), y1(i)
       ENDDO
       CLOSE(kin)
       CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg1)


* quantum yields
* channel 1 – RNO2 + hv -> R + NO2:
* @238nm: 0.2
* @254nm: 0.7
* @296nm: 0.7
* @337nm: 0.3
* below 254nm: lin. interpol. to qy1=0
* above 296nm: lin. interpol. to qy1=0
* in-between: qy1 = 0.7

* channel 2 – RNO2 + hv -> 1-alkene + HONO:
      IF(myld==1) THEN
        qy2 = 0.08
      ELSEIF(myld==2) THEN
        qy2 = 1.0
      ENDIF


* combine
      DO iw = 1, nw - 1
         sig = yg1(iw)
         IF(wc(iw)<254.) THEN
           qy1 = MAX(0.,0.2 + 0.5/16.*(wc(iw)-238.))
          ELSEIF(wc(iw)>296.) THEN
           qy1 = MAX(0.,0.7 - 0.4/41.*(wc(iw)-296.))
          ELSE
           qy1 = 0.7
         ENDIF

         DO i = 1, nz
           sq(j-1,i,iw) = sig * qy1
           sq(j  ,i,iw) = sig * qy2
         ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE mn14(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! C2H5NO2

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for C2H5NO2        =*
*=  photolysis:                                                              =*
*=           C2H5NO2 + hv -> C2H5 + NO2                                      =*
*=                                                                           =*
*=  Absorption cross sections: see options below                             =*
*=  Quantum yields:            estimated same as CH3NO2                      =*
*=                                                                           =*
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
      PARAMETER(kdata=200)

      INTEGER i, n
      REAL x1(kdata), x2(kdata)
      REAL y1(kdata), y2(kdata)

* local

      REAL yg1(kw), yg2(kw), yg(kw)
      REAL sig, qy1, qy2
      INTEGER iw

      INTEGER mabs, myld, xsvers(3), qyvers(3)
      CHARACTER(llog) :: logmsg(nlog)
      CHARACTER(lspc) :: spc


***************** C2H5NO2 photolysis *************************


* Setting photolysis index and cross section/quantum yield options
      j = j+1
      jlabel(j) = 'C2H5NO2 -> C2H5 + NO2'
      j = j+1
      jlabel(j) = 'C2H5NO2 -> C2H4 + HONO'

* cross section options
      spc = "C2H5NO2"
      xsvers = (/2, 2, 2/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "from Wallace-Goodeve 1943"
      logmsg(2) = "from McMillan 1966, pers. comm"
      logmsg(3) ="mean of option Wallace-Goodeve 1943 and McMillan 1966"
      CALL set_option(3,spc,"abs",logmsg,xsvers,mabs)

* quantum yield options
      qyvers = (/1, 2, 1/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "with HONO channel of 0.08"
      logmsg(2) = "with HONO channel scaled externally"
      CALL set_option(2,spc,"yld",logmsg,qyvers,myld)

*** Cross sections
* Wallace-Goodeve 1934
      IF(mabs==1 .or. mabs==3) THEN
         OPEN(UNIT=kin,
     $        FILE='DATAJ1/MCMext/NIT/C2H5NO2_Wallace-Goodeve34.abs',
     $        STATUS='old')
         DO i = 1, 7
           READ(kin,*)
         ENDDO
         n = 55
         DO i = 1, n
           READ(kin,*) x1(i), y1(i)
         ENDDO
         CLOSE(kin)
         CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg1)
       ENDIF

* McMillan 1966
       IF(mabs==2 .or. mabs==3) THEN
         OPEN(UNIT=kin,FILE='DATAJ1/MCMext/NIT/C2H5NO2_McMillan66.abs',
     $        STATUS='old')
         DO i = 1, 7
           READ(kin,*)
         ENDDO
         n = 27
         DO i = 1, n
           READ(kin,*) x2(i), y2(i)
         ENDDO
         CLOSE(kin)
         CALL interpol(x2,y2,kdata,n,nw,wl,jlabel(j),0.,0.,yg2)
      ENDIF

* Average of the above
      IF(mabs==1) THEN
        yg = yg1
      ELSEIF(mabs==2) THEN
        yg = yg2
      ELSEIF(mabs==3) THEN
        DO i = 1,kw
          yg(i) = (yg1(i)+yg2(i))/2.
        ENDDO
      ENDIF


* quantum yields
* channel 1 – RNO2 + hv -> R + NO2:
* @238nm: 0.2
* @254nm: 0.7
* @296nm: 0.7
* @337nm: 0.3
* below 254nm: lin. interpol. to qy1=0
* above 296nm: lin. interpol. to qy1=0
* in-between: qy1 = 0.7

* channel 2 – RNO2 + hv -> 1-alkene + HONO:
      IF(myld==1) THEN
        qy2 = 0.08
      ELSEIF(myld==2) THEN
        qy2 = 1.0
      ENDIF


* combine
      DO iw = 1, nw - 1
         sig = yg(iw)
         IF(wc(iw)<254.) THEN
           qy1 = MAX(0.,0.2 + 0.5/16.*(wc(iw)-238.))
          ELSEIF(wc(iw)>296.) THEN
           qy1 = MAX(0.,0.7 - 0.4/41.*(wc(iw)-296.))
          ELSE
           qy1 = 0.7
         ENDIF

         DO i = 1, nz
           sq(j-1,i,iw) = sig * qy1
           sq(j  ,i,iw) = sig * qy2
         ENDDO
      ENDDO

      END

*=============================================================================*
