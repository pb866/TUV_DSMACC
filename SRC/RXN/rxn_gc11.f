*= This file contains the following subroutines, related to reading/loading
*= the product (cross section) x (quantum yield) for photo-reactions of
*= from fast-JX and the GEOS-CHEM mechanism, which were not yet present in TUV5.2:
*=
*=     gc01 through gc05

*=============================================================================*

      SUBROUTINE gc01(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! NO

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for NO photolysis: =*
*=        NO + hv -> N + O(3P)                                               =*
*=                                                                           =*
*=  Cross section:  MMPI-Mainz UV/VIS Spectral Atlas                         =*
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
      PARAMETER(kdata=180)

      INTEGER i, n
      REAL x1(kdata)
      REAL y1(kdata)

* local

      REAL yg(kw)
!      REAL qy
!      REAL sig
      INTEGER iw

      j = j+1
      jlabel(j) = 'NO -> N + O(3P)'

*----- NO photolysis

* cross sections
      OPEN(UNIT=kin,FILE='DATAJ1/GC11/NO.abs',STATUS='old')
      do i = 1, 17
        read(kin,*)
      enddo

      n  = 164
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
      ENDDO
      CLOSE(kin)

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),y1(1),0.,yg)

* quantum yields
!      qy = 1.

* combine:
      DO iw = 1, nw - 1
        DO i = 1, nz
          sq(j,i,iw) = yg(iw)! xs * qy
        ENDDO
      ENDDO

      END

* ==================================================================== *

      SUBROUTINE gc02(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! OCS

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for OCS photolysis:=*
*=        OCS + hv -> CO + S                                                 =*
*=                                                                           =*
*=  Cross section:  JPL No. 18                                               =*
*=  Quantum yield:  JPL No. 18 (for 220 - 250nm)                             =*
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
      PARAMETER(kdata=180)

      INTEGER i, n, n2
      REAL x1(kdata), x2(kdata)
      REAL y1(kdata), y2(kdata)

* local

      REAL yg(kw), yg2(kw)
!      REAL qy
      REAL sig
      INTEGER iw

      j = j+1
      jlabel(j) = 'OCS -> CO + SO2'

*----- OCS photolysis

* cross sections
      OPEN(UNIT=kin,FILE='DATAJ1/GC11/OCS_jpl15.abs',STATUS='old')
      do i = 1, 11
        read(kin,*)
      enddo

      n  = 41
      n2 = n
      DO i = 1, n
        READ(kin,*) x1(i), y1(i), y2(i)
        x2(i) = x1(i)
      ENDDO
      CLOSE(kin)

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),y1(1),0.,yg)

      CALL interpol(x2,y2,kdata,n2,nw,wl,jlabel(j),0.,0.,yg2)

* quantum yields
!      qy = 1.

* combine:
      DO iw = 1, nw - 1
        DO i = 1, nz
          IF(tlev(i)<=225.) THEN
            sig = yg2(iw)
          ELSEIF(tlev(i)>=295.) THEN
            sig = yg(iw)
          ELSE
            sig = yg(iw)+(yg(iw)-yg2(iw))/70.*(tlev(i)-295.)
          ENDIF
          sq(j,i,iw) = sig! xs * qy
        ENDDO
      ENDDO

      END

* ==================================================================== *

      SUBROUTINE gc03(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! H2SO4

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for H2SO4          =*
*=  photolysis:                                                              =*
*=        H2SO4 + hv -> 2 OH + SO2                                           =*
*=                                                                           =*
*=  Cross section:  MMPI-Mainz UV/VIS Spectral Atlas                         =*
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
      PARAMETER(kdata=180)

      INTEGER i, n
      REAL x1(kdata)
      REAL y1(kdata)

* local

      REAL yg(kw)
!      REAL qy
!      REAL sig
      INTEGER ierr
      INTEGER iw

      j = j+1
      jlabel(j) = 'H2SO4 -> 2 OH + SO2'

*----- H2SO4 photolysis

* QY = 1
* cross sections
      OPEN(UNIT=kin,FILE='DATAJ1/GC11/NO.abs',STATUS='old')
      do i = 1, 17
        read(kin,*)
      enddo

      n  = 28
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
      ENDDO
      CLOSE(kin)

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)

* quantum yields
!      qy = 1.

* combine:
      DO iw = 1, nw - 1
        DO i = 1, nz
          sq(j,i,iw) = yg(iw)! xs * qy
        ENDDO
      ENDDO

      END

* ==================================================================== *

      SUBROUTINE gc04(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! MACR

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for CH2=C(CH3)CHO      =*
*=  (methacrolein) overall photolysis reaction:                              =*
*=       CH2=C(CH3)CHO + hv -> Products                                      =*
*=                                                                           =*
*=  Cross section: from JPL 2006 (originally from Gierczak et al.)           =*
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

      INTEGER iw
      INTEGER i, n
      REAL x1(kdata)
      REAL y1(kdata)

* local

      REAL yg(kw), yg1(kw), dum
      real qy
      REAL sig

**************** methacrolein photodissociation

* same channels estimated as for acrolein
      j = j+1
      jlabel(j) = 'MACR-> products'

* cross section from
*      JPL 2006 (originally from Gierczak et al.), also IUPAC recommended

      OPEN(UNIT=kin,FILE='DATAJ1/ABS/Methacrolein_jpl11.abs',
     $    STATUS='OLD')
      DO i = 1, 7
         READ(kin,*)
      ENDDO
      n = 146
      DO i = 1, n
         READ(kin,*) x1(i), y1(i)
         y1(i) = y1(i) * 1.E-20
      ENDDO
      CLOSE(kin)

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)


      OPEN(UNIT=kin,FILE='DATAJ1/ABS/MACR_calv.dat',
     $     STATUS='OLD')
      DO i = 1, 35
        READ(kin,*)
      ENDDO
      n = 90
      DO i = 1, n
        READ(kin,*) x1(i), dum, y1(i)
      ENDDO
      CLOSE(kin)

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg1)


      DO iw = 1, nw-1
        sig = yg(iw)
* quantum yield:
        qy = yg1(iw)
        DO i = 1, nz
          sq(j,i,iw) = qy * sig
        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE gc05(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! CH2Br2

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for CH2Br2         =*
*=  photolysis:                                                              =*
*=        CH2Br2 + hv -> CH2Br + Br                                          =*
*=                                                                           =*
*=  Cross section:  IUPAC(2008)                                              =*
*=  Quantum yield:  IUPAC(2008)                                              =*
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
      PARAMETER(kdata=180)

      INTEGER i, n
      REAL x1(kdata)
      REAL y1(kdata)

* local

      REAL yg(kw), yg2(kw)
!      REAL qy
      REAL sig
      INTEGER iw

*----- CH2Br2 photolysis
      j = j+1
      jlabel(j) = 'CH2Br2 -> CH2Br + Br'


* cross sections
      OPEN(UNIT=kin,FILE='DATAJ1/GC11/CH2Br2_295K_iup08.abs',
     &     STATUS='old')
      do i = 1, 13
        read(kin,*)
      enddo

      n  = 46
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
      ENDDO
      CLOSE(kin)
      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),y1(1),0.,yg)

      OPEN(UNIT=kin,FILE='DATAJ1/GC11/CH2Br2_210K_iup08.abs',
     &     STATUS='old')
      do i = 1, 14
        read(kin,*)
      enddo

      n  = 33
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
      ENDDO
      CLOSE(kin)
      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg2)

* quantum yields
!      qy = 1.

* combine:
      DO iw = 1, nw - 1
        DO i = 1, nz
          IF(tlev(i)<=210.) THEN
            sig = yg2(iw)
          ELSEIF(tlev(i)>=295.) THEN
            sig = yg(iw)
          ELSE
            sig = yg(iw)+(yg(iw)-yg2(iw))/85.*(tlev(i)-295.)
          ENDIF
          sq(j,i,iw) = sig! xs * qy
        ENDDO
      ENDDO

      END

* ==================================================================== *
