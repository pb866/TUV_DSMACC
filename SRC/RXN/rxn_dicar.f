*= This file contains the following subroutines, related to reading/loading
*= the product (cross section) x (quantum yield) for photo-reactions of
*= dicarbonyls in MCM-GECKO, which were not yet present in TUV5.2:
*=
*=     mb01 through mb06

*=============================================================================*

      SUBROUTINE mb01(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! butenedial

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for                =*
*=  CHOCH=CHCHO photolysis:                                                  =*
*=        HC(O)CH=CHCHO + hv -> 3H-furan-2-one                               =*
*=                                                                           =*
*=  Cross section:  trans-butenedial (see options below)                     =*
*=  Quantum yield:  trans-butenedial (see options below)                     =*
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

      REAL yg(kw), yg1(kw), dum
      REAL qy
      REAL sig
      INTEGER iw

      INTEGER mabs, myld, xsvers(3), qyvers(3)
      CHARACTER(llog) :: logmsg(nlog)
      CHARACTER(lspc) :: spc


***************** butenedial photolysis *************************

* Setting photolysis index and cross section/quantum yield options
      j = j+1
      jlabel(j) = 'CHOCH=CHCHO -> 3H-furan-2-one'
      j = j+1
      jlabel(j) = 'uDICARaa(poly)'

* cross section options
      spc = "CHOCH=CHCHO"
      xsvers = (/1, 2, 2/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "from SAPRC-99"
      logmsg(2) = "from IUPAC"
      CALL set_option(2,spc,"abs",logmsg,xsvers,mabs)

* quantum yield options
      qyvers = (/1, 2, 2/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      ! logmsg same as for xs
      CALL set_option(2,spc,"yld",logmsg,qyvers,myld)


* cross sections

      IF(mabs==1) THEN
        OPEN(UNIT=kin,FILE='DATAJ1/MCMext/DICAR/butenedial1',
     $       STATUS='old')
        do i = 1, 6
          read(kin,*)
        enddo

        n  = 56
        DO i = 1, n
          READ(kin,*) x1(i), y1(i)
        ENDDO
        CLOSE(kin)

       ELSEIF(mabs==2) THEN

        OPEN(UNIT=kin,FILE='DATAJ1/MCMext/DICAR/ButDial_iup.abs',
     $       STATUS='old')
        do i = 1, 5
          read(kin,*)
        enddo

        n  = 33
        DO i = 1, n
          READ(kin,*) x1(i), y1(i)
          y1(i) = y1(i)*1.e-20
        ENDDO
        CLOSE(kin)
      ENDIF

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)

* quantum yields

      IF(myld==1) THEN
        OPEN(UNIT=kin,FILE='DATAJ1/MCMext/DICAR/butenedial1',
     $       STATUS='old')
        do i = 1, 6
          read(kin,*)
        enddo

        n  = 56
        DO i = 1, n
          READ(kin,*) x1(i), dum, y1(i)
        ENDDO
        CLOSE(kin)

        CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg1)
      ENDIF


* combine:

      DO iw = 1, nw - 1
        sig = yg(iw)
        IF(myld==1) THEN
          qy = yg1(iw)
         ELSEIF(myld==2) THEN
          !cis/trans conversion channel available from IUPAC (<0.4)
          IF(wc(iw)>248.) THEN
            qy = 0.012
           ELSE
            qy = 0.028
          ENDIF
        ENDIF
        DO i = 1, nz
          sq(j-1,i,iw) = sig * qy
          sq(j  ,i,iw) = sig
        ENDDO
      ENDDO

      END

* =============================================================================

      SUBROUTINE mb02(nw,wl,wc,nz,tlev,airlev,j,sq,jlabel) ! 4-oxo-2-pentenal

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for                =*
*=  CH3C(O)CH=CHCHO photolysis:                                              =*
*=        CH3C(O)CH=CHCHO + hv -> 5methyl-3H-furan-2-one                     =*
*=                                                                           =*
*=  Cross section:  4 oxo2pentenal Bierbach 94                               =*
*=  Quantum yield:  4 oxo2pentenal Bierbach 94                               =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE '../../params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)

      INTEGER nz

      REAL tlev(kz)
      REAL airlev(kz)

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

      REAL yg(kw), yg1(kw), dum
      REAL qy1,qy2,qy3,qy4
      REAL sig
      INTEGER iw

      INTEGER mabs, myld, xsvers(3), qyvers(3)
      CHARACTER(llog) :: logmsg(nlog)
      CHARACTER(lspc) :: spc


***************** 4-oxo-2-pentenal photolysis *************************

* Setting photolysis index and cross section/quantum yield options
      j = j+1
      jlabel(j) = 'CH3COCH=CHCHO -> 5Me-3H-2-furanone'
      j = j+1
      jlabel(j) = 'CH3COCH=CHCHO -> CH3 + CHOCH=CHCO'
      j = j+1
      jlabel(j) = 'CH3COCH=CHCHO -> CH3COCH=CH2 + CO'
      j = j+1
      jlabel(j) = 'CH3COCH=CHCHO -> maleic anhydride + HO2. + R.'
      j = j+1
      jlabel(j) = 'uDICARak(poly)'
* further channel possible: cis/trans conversion with qy < 0.2 (see IUPAC)

* cross section options
      spc = "CH3COCH=CHCHO"
      xsvers = (/1, 2, 2/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "from SAPRC-99"
      logmsg(2) = "from IUPAC"
      CALL set_option(2,spc,"abs",logmsg,xsvers,mabs)

* quantum yield options
      qyvers = (/1, 2, 2/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      ! logmsg same as for xs
      CALL set_option(2,spc,"yld",logmsg,qyvers,myld)


* cross sections
* currently not distinguished between cis/trans

      IF(mabs == 1) THEN
        OPEN(UNIT=kin,FILE='DATAJ1/MCMext/DICAR/4oxo2pentenal1',
     $       STATUS='old')
        do i = 1, 6
          read(kin,*)
        enddo

        n = 56
        DO i = 1, n
          READ(kin,*) x1(i), y1(i)
        ENDDO
        CLOSE(kin)

       ELSEIF(mabs == 2) THEN

        OPEN(UNIT=kin,FILE='DATAJ1/MCMext/DICAR/4oPentenal_iup.abs',
     $       STATUS='old')
        do i = 1, 6
          read(kin,*)
        enddo

        n = 21
        DO i = 1, n
          READ(kin,*) x1(i), dum, dum, y1(i)
          y1(i) = y1(i) * 1E-20
        ENDDO
        CLOSE(kin)
      ENDIF
      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j-3),0.,0.,yg)


* quantum yields from data file

      IF(myld == 1) THEN
        OPEN(UNIT=kin,FILE='DATAJ1/MCMext/DICAR/4oxo2pentenal1',
     $       STATUS='old')
        do i = 1, 6
          read(kin,*)
        enddo

        n = 56
        DO i = 1, n
          READ(kin,*) x1(i), dum, y1(i)
        ENDDO
        CLOSE(kin)
        CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j-3),0.,0.,yg1)
      ENDIF


* combine xs and qy:

      DO iw = 1, nw - 1
* cross section:
        sig = yg(iw)
* quantum yields:
        IF(myld == 1) THEN
          qy1 = yg1(iw)
          qy2 = 0.
          qy3 = 0.
          qy4 = qy1
         ELSEIF(myld==2) THEN
          qy1 = 0.05
          IF(wl(iw)<=308.) THEN
            qy2 = 0.3
            qy3 = 0.4
           ELSEIF(wl(iw)>=351.) THEN
            qy2 = 0.23
            qy3 = 0.33
           ELSE
            qy2 = 0.3 + (0.23-0.3)*(wl(iw)-308.)/(351.-308.)
            qy3 = 0.4 + (0.33-0.4)*(wl(iw)-308.)/(351.-308.)
          ENDIF
          qy4 = 0.
        ENDIF


         DO i = 1, nz
            sq(j-4,i,iw) = sig * qy1
            sq(j-3,i,iw) = sig * qy2
            sq(j-2,i,iw) = sig * qy3
            sq(j-1,i,iw) = sig * qy4
            sq(j  ,i,iw) = sig
         ENDDO
      ENDDO

      END

* =============================================================================

      SUBROUTINE mb03(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! hexadienedial

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for                =*
*=  CHOCH=CHCH=CHCHO photolysis:                                             =*
*=        HC(O)CH=CHCH=CHCHO + hv -> Z-3,4-Diformyl-cyclobutene              =*
*=                                                                           =*
*=  Cross section:  see options below                                        =*
*=  Quantum yield:  0.1 (estimate)                                           =*
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


***************** hexadienedial photolysis *************************

* Setting photolysis index and cross section options
      j = j+1
      jlabel(j) = 'CHOCH=CHCH=CHCHO -> diformyl cyclobutene'

* cross section options
      spc = "CHOCH=CHCH=CHCHO"
      xsvers = (/1, 2, 2/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "from Klotz et al. 1995"
      logmsg(2) = "from Xiang and Zhu 2007"
      CALL set_option(2,spc,"abs",logmsg,xsvers,mabs)

* quantum yield options
      qyvers = (/1, 2, 1/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "estimated 0.1"
      logmsg(2) = "scaled externally"
      CALL set_option(2,spc,"yld",logmsg,qyvers,myld)


* cross sections

      IF(mabs==1) THEN
        OPEN(UNIT=kin,FILE='DATAJ1/MCMext/DICAR/EEhexadienedial1.prn',
     $       STATUS='old')
        do i = 1, 5
          read(kin,*)
        enddo

        n  = 13
        DO i = 1, n
          READ(kin,*) x1(i), y1(i)
        ENDDO
        CLOSE(kin)

       ELSEIF(mabs==2) THEN

        OPEN(UNIT=kin,
     $       FILE='DATAJ1/MCMext/DICAR/hexadienedial_X+Z07.abs',
     $       STATUS='old')
        do i = 1, 6
          read(kin,*)
        enddo

        n  = 29
        DO i = 1, n
          READ(kin,*) x1(i), y1(i)
        ENDDO
        CLOSE(kin)
      ENDIF

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)

* quantum yields
      IF(myld==1) THEN
        qy = 0.1
      ELSEIF(myld==2) THEN
        qy = 1.0
      ENDIF


* combine:

      DO iw = 1, nw - 1
        sig = yg(iw)
        DO i = 1, nz
          sq(j,i,iw) = sig*qy
        ENDDO
      ENDDO

      END

* =============================================================================

      SUBROUTINE mb04(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! 3-hexene-2,5-dione

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for                =*
*=  CH3COCH=CHCOCH3 photolysis:                                              =*
*=        CH3COCH=CHCOCH3 + hv -> CH3CO + CH=CHCOCH3                         =*
*=                                                                           =*
*=  Cross section:  Liu et al. 1999 in acetonitrile                          =*
*=  Quantum yield:  estimated 0.1                                            =*
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

      INTEGER   myld, qyvers(3)
      CHARACTER(llog) :: logmsg(nlog)
      CHARACTER(lspc) :: spc


***************** CH3COCH=CHCOCH3 photolysis *************************

* Setting photolysis index and quantum yield options
      j = j+1
      jlabel(j) = 'CH3COCH=CHCOCH3 -> CH3CO + CH=CHCOCH3'


* quantum yield options
      spc = "CH3COCH=CHCOCH3"
      qyvers = (/1, 2, 1/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "estimated 0.1"
      logmsg(2) = "scaled externally"
      CALL set_option(2,spc,"yld",logmsg,qyvers,myld)


* cross sections



* main fraction >~80% undergoes cis/trans-isomerisation

      OPEN(UNIT=kin,FILE='DATAJ1/MCMext/DICAR/hexenedione.abs',
     $     STATUS='old')
      do i = 1, 9
        read(kin,*)
      enddo

      n  = 160
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
        y1(i) = y1(i)*1.E-20
      ENDDO
      CLOSE(kin)

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)


* quantum yields
      IF(myld==1) THEN
        qy = 0.1
      ELSEIF(myld==2) THEN
        qy = 1.0
      ENDIF


* combine:

      DO iw = 1, nw - 1
        sig = yg(iw)
        DO i = 1, nz
          sq(j  ,i,iw) = sig*qy
        ENDDO
      ENDDO

      END

* =============================================================================

      SUBROUTINE mb05(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! pinonaldehyde

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section)x(quantum yield) for                      =*
*=  pinonaldehyde photolysis                                                 =*
*=          pinonaldehyde + hv -> R + CO + HO2                               =*
*=                                                                           =*
*=  Cross section from IUPAC                                                 =*
*=  Quantum yield see options below                                          =*
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

      INTEGER i, n
      REAL x1(kdata)
      REAL y1(kdata)

* local

      REAL yg(kw)
      REAL qy
      REAL sig
      INTEGER iw

      INTEGER myld, qyvers(3)
      CHARACTER(llog) :: logmsg(nlog)
      CHARACTER(lspc) :: spc


***************** pinonaldehyde photolysis *************************

* Setting photolysis index and quantum yield options
      j = j+1
      jlabel(j) = 'pinonaldehyde -> R + CO + HO2'

      spc = "pinonaldehyde"
      qyvers = (/1, 4, 3/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = 'from EUPHORE chamber experiments' ! (RADICAL 2002)
      logmsg(2) = 'from Jaoui and Kamens 2003'
      logmsg(3) = 'as approximate average of opt. 1 and 2'
      logmsg(4) = 'scaled externally'
      CALL set_option(4,spc,"yld",logmsg,qyvers,myld)


      OPEN(UNIT=kin,FILE='DATAJ1/MCMext/DICAR/Pinonaldehyde_iup.abs',
     $     STATUS='old')
      do i = 1, 5
         read(kin,*)
      enddo

      n = 14
      DO i = 1, n
         READ(kin,*) x1(i), y1(i)
         y1(i) = y1(i) * 1E-20
      ENDDO
      CLOSE(kin)

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j-1),0.,0.,yg)


* quantum yields
      IF(myld == 1) THEN
        qy = 0.14
       ELSEIF(myld == 2) THEN
        qy = 0.4
       ELSEIF(myld == 3) THEN
        qy = 0.25
       ELSEIF(myld == 4) THEN
        qy = 1.0
      ENDIF


* combine:

      DO iw = 1, nw - 1
         sig = yg(iw)
         DO i = 1, nz
            sq(j,i,iw) = sig * qy
         ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE mb06(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! caronaldehyde

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section)x(quantum yield) for                      =*
*=  caronaldehyde photolysis                                                 =*
*=          caronaldehyde + hv -> R + CO + HO2                               =*
*=                                                                           =*
*=  Cross section from Hallquist et al. 1997                                 =*
*=  Quantum yield estimated same as pinonaldehyde                            =*
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

      INTEGER i, n
      REAL x1(kdata)
      REAL y1(kdata)

* local

      REAL yg(kw)
      REAL qy
      REAL sig
      INTEGER iw

      INTEGER myld, qyvers(3)
      CHARACTER(llog) :: logmsg(nlog)
      CHARACTER(lspc) :: spc


***************** Caronaldehyde photolysis *************************


* Setting photolysis index and quantum yield options
      j = j+1
      jlabel(j) = 'caronaldehyde -> R + CO + HO2'

* quantum yield options
* Criegee channel replaced by CH3 fission
      spc = "caronaldehyde"
      qyvers = (/1, 4, 3/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "from pinonaldehyde EUPHORE chamber experiments"
      logmsg(2) = "same as pinonaldehyde from Jaoui and Kamens 2003"
      logmsg(3) = "as approximate average of opt. 1 and 2"
      logmsg(4) = 'scaled externally'
      CALL set_option(4,spc,"yld",logmsg,qyvers,myld)


      OPEN(UNIT=kin,FILE='DATAJ1/MCMext/DICAR/caronaldehyde.abs',
     $     STATUS='old')
      do i = 1, 6
         read(kin,*)
      enddo

      n = 11
      DO i = 1, n
         READ(kin,*) x1(i), y1(i)
      ENDDO
      CLOSE(kin)

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j-1),0.,0.,yg)


* quantum yields
      IF(myld == 1) THEN
        qy = 0.14
       ELSEIF(myld == 2) THEN
        qy = 0.4
       ELSEIF(myld == 3) THEN
        qy = 0.25
       ELSEIF(myld == 4) THEN
        qy = 1.0
      ENDIF


* combine:

      DO iw = 1, nw - 1
         sig = yg(iw)
         DO i = 1, nz
            sq(j,i,iw) = sig * qy
         ENDDO
      ENDDO

      END

*=============================================================================*
