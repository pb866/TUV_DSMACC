*= This file contains the following subroutines, related to reading/loading
*= the product (cross section) x (quantum yield) for photo-reactions of
*= organic nitrates in MCM-GECKO, which were not yet present in TUV5.2:
*=
*=     mn01 through mn11

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
      INCLUDE 'params'

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
      INTEGER ierr,idum
      INTEGER iw
      INTEGER mabs

***************** n-C5H11ONO2 photolysis *************************

      j = j+1
      jlabel(j) = 'n-C5H11ONO2 -> n-C5H11O + NO2'


      IF(vers==1)THEN
        mabs = 2 ! GECKO-A database
       ELSEIF(vers==2)THEN
        mabs = 3
       ELSEIF(vers==0) THEN
        mabs = 2
       ELSE
        STOP "'vers' not set. Choose value between 0 and 2 in 'params'."
      ENDIF

      IF(vers==1 .OR. vers==2) THEN
        CONTiNUE
       ELSEIF(mabs.EQ.1) THEN
        WRITE(kout,'(A)')
     &       ' n-C5H11ONO2 cross sections from Clemitshaw et al. 1996.'
       ELSEIF(mabs.EQ.2) THEN
        WRITE(kout,'(A)')
     &       ' n-C5H11ONO2 cross sections from Zhu and Kellis 1997.'
       ELSEIF(mabs.EQ.3) THEN
        WRITE(kout,'(2A)')
     &       ' n-C5H11ONO2 cross sections average of',
     &       ' Zhu and Kellis/Clemitshaw.'
       ELSE
        STOP "'mabs' not defined for n-C5H11ONO2 photolysis."
      ENDIF

      IF( mabs.EQ.1 .OR. mabs.EQ.3) THEN
         OPEN(UNIT=kin,FILE='DATAJ1/MCMext/NIT/1PentNit_Cle96.abs',
     $        STATUS='old')
         DO i = 1, 6
            READ(kin,*)
         ENDDO
         n = 15
         DO i = 1, n
            READ(kin,*) idum, y1(i)
            x1(i) = float(idum)
         ENDDO
         CLOSE(kin)
         IF(mabs==3) THEN
           ys1(:) = y1(:)
          ELSE
           CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
           CALL addpnt(x1,y1,kdata,n,               0.,0.)
           CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
           CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
           CALL inter2(nw,wl,yg1,n,x1,y1,ierr)
           IF (ierr .NE. 0) THEN
             WRITE(*,*) ierr, jlabel(j)
             STOP
           ENDIF
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
            READ(kin,*) idum, y1(i)
            x1(i) = float(idum)
         ENDDO
         CLOSE(kin)
         IF(mabs==3) THEN
           xs2(:) = x1(:)
           ys2(:) = y1(:)
          ELSE
           CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
           CALL addpnt(x1,y1,kdata,n,               0.,0.)
           CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
           CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
           CALL inter2(nw,wl,yg1,n,x1,y1,ierr)
           IF (ierr .NE. 0) THEN
             WRITE(*,*) ierr, jlabel(j)
             STOP
           ENDIF
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

        CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
        CALL addpnt(x1,y1,kdata,n,               0.,0.)
        CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
        CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
        CALL inter2(nw,wl,yg1,n,x1,y1,ierr)
        IF (ierr .NE. 0) THEN
          WRITE(*,*) ierr, jlabel(j)
          STOP
        ENDIF
      ENDIF

! Temperature dependency (assume for all cases)
      OPEN(UNIT=kin,FILE='DATAJ1/MCMext/NIT/1PentNit_Z+K97_Tdep.abs',
     $     STATUS='old')
      DO i = 1, 9
        READ(kin,*)
      ENDDO
      n = 16
      DO i = 1, n
        READ(kin,*) idum, y2(i)
        x2(i) = float(idum)
        y2(i) = y2(i) * 1.E-3
      ENDDO
      CLOSE(kin)

      CALL addpnt(x2,y2,kdata,n,x2(1)*(1.-deltax),0.)
      CALL addpnt(x2,y2,kdata,n,               0.,0.)
      CALL addpnt(x2,y2,kdata,n,x2(n)*(1.+deltax),y2(n))
      CALL addpnt(x2,y2,kdata,n,           1.e+38,y2(n))
      CALL inter2(nw,wl,yg2,n,x2,y2,ierr)
      IF (ierr .NE. 0) THEN
        WRITE(*,*) ierr, jlabel(j)
        STOP
      ENDIF


* quantum yield  = 1

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
      INCLUDE 'params'

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

**************** 2-C5H11ONO2 photodissociation

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
      INCLUDE 'params'

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

**************** 3-C5H11ONO2 photodissociation

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
      INCLUDE 'params'

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

**************** c-C5H11ONO2 photodissociation

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
      INCLUDE 'params'

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
      INTEGER ierr, idum
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
        READ(kin,*) idum, y1(i)
        x1(i) = float(idum)
      ENDDO
      CLOSE(kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,               0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
      CALL inter2(nw,wl,yg1,n,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
        WRITE(*,*) ierr, jlabel(j)
        STOP
      ENDIF


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
      INCLUDE 'params'

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
      INTEGER ierr, idum
      INTEGER iw

***************** i-C5H11ONO2 photolysis *************************

      j = j+1
      jlabel(j) = 'i-C5H11ONO2 -> i-C5H11O + NO2'

* cross sections
      OPEN(UNIT=kin,FILE='DATAJ1/MCMext/NIT/isopentylNit.abs',
     $     STATUS='old')
      DO i = 1, 6
        READ(kin,*)
      ENDDO
      n = 25
      DO i = 1, n
        READ(kin,*) idum, y1(i)
        x1(i) = float(idum)
      ENDDO
      CLOSE(kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,               0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
      CALL inter2(nw,wl,yg1,n,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
        WRITE(*,*) ierr, jlabel(j)
        STOP
      ENDIF


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
      INCLUDE 'params'

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
      INTEGER n1, n2
      REAL x1(kdata), x2(kdata)
      REAL y1(kdata), y2(kdata)

* local

      REAL yg(kw), yg1(kw), yg2(kw), yg3(kw)
!      REAL qy
      REAL sig, fOH
      INTEGER ierr

      INTEGER mabs

**************** RNO3 photodissociation

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

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,               0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
      CALL inter2(nw,wl,yg,n,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
        WRITE(*,*) ierr, jlabel(j)
        STOP
      ENDIF

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

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,               0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
      CALL inter2(nw,wl,yg1,n,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
        WRITE(*,*) ierr, jlabel(j)
        STOP
      ENDIF


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

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),y1(1))
      CALL addpnt(x1,y1,kdata,n,               0.,y1(1))
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
      CALL inter2(nw,wl,yg2,n,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
        WRITE(*,*) ierr, jlabel(j)
        STOP
      ENDIF

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

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),y1(1))
      CALL addpnt(x1,y1,kdata,n,               0.,y1(1))
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
      CALL inter2(nw,wl,yg3,n,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
        WRITE(*,*) ierr, jlabel(j)
        STOP
      ENDIF

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
      INCLUDE 'params'

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
      INTEGER n1, n2
      REAL x1(kdata), x2(kdata)
      REAL y1(kdata), y2(kdata)

* local

      REAL yg1(kw), yg2(kw), yg3(kw)
!      REAL qy
      REAL sig, fOH
      INTEGER ierr

      INTEGER mabs

**************** RNO3 photodissociation

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

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,               0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
      CALL inter2(nw,wl,yg1,n,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF


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

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),y1(1))
      CALL addpnt(x1,y1,kdata,n,               0.,y1(1))
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
      CALL inter2(nw,wl,yg2,n,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
        WRITE(*,*) ierr, jlabel(j)
        STOP
      ENDIF

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

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),y1(1))
      CALL addpnt(x1,y1,kdata,n,               0.,y1(1))
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
      CALL inter2(nw,wl,yg3,n,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
        WRITE(*,*) ierr, jlabel(j)
        STOP
      ENDIF

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
      INCLUDE 'params'

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
      INTEGER n1, n2
      REAL x1(kdata), x2(kdata)
      REAL y1(kdata), y2(kdata)

* local

      REAL yg1(kw), yg2(kw), yg3(kw)
!      REAL qy
      REAL sig, fOH
      INTEGER ierr

      INTEGER mabs

**************** RNO3 photodissociation

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

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,               0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
      CALL inter2(nw,wl,yg1,n,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
        WRITE(*,*) ierr, jlabel(j)
        STOP
      ENDIF


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

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),y1(1))
      CALL addpnt(x1,y1,kdata,n,               0.,y1(1))
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
      CALL inter2(nw,wl,yg2,n,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
        WRITE(*,*) ierr, jlabel(j)
        STOP
      ENDIF

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

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),y1(1))
      CALL addpnt(x1,y1,kdata,n,               0.,y1(1))
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
      CALL inter2(nw,wl,yg3,n,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
        WRITE(*,*) ierr, jlabel(j)
        STOP
      ENDIF

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
      INCLUDE 'params'

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
      INTEGER n1, n2
      REAL x1(kdata), x2(kdata)
      REAL y1(kdata), y2(kdata)

* local

      REAL yg2(kw), yg3(kw)
!      REAL qy
      REAL sig, fOH, a, b, c
      INTEGER ierr

      INTEGER mabs

**************** RNO3 photodissociation

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

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),y1(1))
      CALL addpnt(x1,y1,kdata,n,               0.,y1(1))
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
      CALL inter2(nw,wl,yg2,n,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
        WRITE(*,*) ierr, jlabel(j)
        STOP
      ENDIF

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

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),y1(1))
      CALL addpnt(x1,y1,kdata,n,               0.,y1(1))
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
      CALL inter2(nw,wl,yg3,n,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
        WRITE(*,*) ierr, jlabel(j)
        STOP
      ENDIF

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
      INCLUDE 'params'

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
      INTEGER n1, n2
      REAL x1(kdata), x2(kdata)
      REAL y1(kdata), y2(kdata)

* local

      REAL yg(kw), yg1(kw), yg2(kw), yg3(kw)
      REAL qyNO2, qyNO3
      REAL sig, sig1, sig2
      INTEGER ierr

      INTEGER mabs

**************** RNO3 photodissociation

      j = j+1
      jlabel(j) = 'PAN -> RCO(OO) + NO2'
      j = j+1
      jlabel(j) = 'PAN -> RCO(O) + NO3'

* Only MCM/GECKO-A choices considered for generic species:

* PAN cross sections

      OPEN(UNIT=kin,FILE='DATAJ1/RONO2/PAN_talukdar.abs',STATUS='OLD')
      DO i = 1, 14
         READ(kin,*)
      ENDDO
      n = 78
      DO i = 1, n
         READ(kin,*) x1(i), y1(i), y2(i)
         y1(i) = y1(i) * 1.E-20
         y2(i) = y2(i) * 1E-3
         x2(i) = x1(i)
      ENDDO
      n2 = n
      CLOSE(kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,               0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
      CALL inter2(nw,wl,yg,n,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j-1)
         STOP
      ENDIF

      CALL addpnt(x2,y2,kdata,n2,x2(1)*(1.-deltax) ,0.)
      CALL addpnt(x2,y2,kdata,n2,                0.,0.)
      CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),0.)
      CALL addpnt(x2,y2,kdata,n2,            1.e+38,0.)
      CALL inter2(nw,wl,yg1,n2,x2,y2,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF

* PPN cross sections

      OPEN(UNIT=kin,FILE='DATAJ1/ABS/PPN_Harwood.txt',STATUS='OLD')
      DO i = 1, 10
         READ(kin,*)
      ENDDO
      n = 66
      DO i = 1, n
         READ(kin,*) x1(i), y1(i), y2(i)
         y1(i) = y1(i) * 1.E-20
         y2(i) = y2(i) * 1E-3
         x2(i) = x1(i)
      ENDDO
      n2 = n
      CLOSE(kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,               0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
      CALL inter2(nw,wl,yg2,n,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF

      CALL addpnt(x2,y2,kdata,n2,x2(1)*(1.-deltax),0.)
      CALL addpnt(x2,y2,kdata,n2,          0.,0.)
      CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),0.)
      CALL addpnt(x2,y2,kdata,n2,      1.e+38,0.)
      CALL inter2(nw,wl,yg3,n2,x2,y2,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF

* quantum yields
      qyNO2 = 0.61
      qyNO3 = 0.39


      DO iw = 1, nw - 1
        DO i = 1, nz
          sig1 = yg(iw) * EXP(yg1(iw)*(tlev(i)-298.))
          sig2 = yg2(iw) * EXP(yg3(iw)*(tlev(i)-298.))
          sig = (sig1+sig2)/2
          sq(j,i,iw) = sig * qyNO2
          sq(j,i,iw) = sig * qyNO3
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
      INCLUDE 'params'

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

**************** 3-C5H11ONO2 photodissociation

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
