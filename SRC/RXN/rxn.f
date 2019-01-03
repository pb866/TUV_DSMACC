* This file contains the following subroutines, related to reading/loading
* the product (cross section) x (quantum yield) for photo-reactions:
*     r01 through r47
*     r101 through r105
*     r107 through r115
*     r118 through r145
*     pxCH2O
*=============================================================================*

      SUBROUTINE r01(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! O3

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product of (cross section) x (quantum yield) for the two     =*
*=  O3 photolysis reactions:                                                 =*
*=             (a) O3 + hv -> O2 + O(1D)                                     =*
*=             (b) O3 + hv -> O2 + O(3P)                                     =*
*=  Cross section:  Combined data from WMO 85 Ozone Assessment (use 273K     =*
*=                  value from 175.439-847.5 nm) and data from Molina and    =*
*=                  Molina (use in Hartley and Huggins bans (240.5-350 nm)   =*
*=                  IUPAC recommendations added (pb)                         =*
*=  Quantum yield:  Choice between                                           =*
*=                   (1) data from Michelsen et al, 1994                     =*
*=                   (2) JPL 87 recommendation                               =*
*=                   (3) JPL 90/92 recommendation (no "tail")                =*
*=                   (4) data from Shetter et al., 1996                      =*
*=                   (5) JPL 97 recommendation                               =*
*=                   (6) JPL 00 recommendation                               =*
*=                   (7) Matsumi et al., 2002                                =*
*=                   (8) IUPAC recommendation                                =*
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

      INTEGER n1, n2, n3, n4
      INTEGER kdata
      PARAMETER (kdata = 500)
      REAL x1(kdata), x2(kdata), x3(kdata), x4(kdata)
      REAL y1(kdata), y2(kdata), y3(kdata), y4(kdata)

* local

      REAL xs(kz,kw)

      REAL yg1(kw), yg2(kw), yg3(kw), yg4(kw)
      REAL qy1d, qy3p
      REAL tau, tau2, tau3
      REAL a, b, c
      REAL a0, a1, a2, a3, a4, a5, a6
      REAL xl, xl0
      INTEGER kmich, kjpl87, kjpl92, kshet, kjpl97, kjpl00, kmats, kiup
      INTEGER i, iw, n, idum

      REAL fo3qy, fo3qy2
      EXTERNAL fo3qy, fo3qy2

      INTEGER mabs, myld, xsvers(3), qyvers(3)
      CHARACTER(llog) :: logmsg(nlog)
      CHARACTER(lspc) :: spc


***************** ozone photolysis *************************

* Setting photolysis index and cross section/quantum yield options
      j = j + 1
      jlabel(j) = 'O3 -> O2 + O(1D)'
      j = j + 1
      jlabel(j) = 'O3 -> O2 + O(3P)'

* cross section options
* call cross section read/interpolate routine
* cross sections from WMO 1985 Ozone Assessment
* from 175.439 to 847.500 nm. Using value at 273 K.
* Values are over-written in Hartly and Huggins bands, using different
* options depending on value of mabs:

*     mabs = 1: mostly Reims grp (Malicet, Brion) (TUV standard)
*     mabs = 2: JPL 2006
*     mabs = 3: exclusively Molina and Molina
*     mabs = 4: exclusively Bass et al.
*     mabs = 5: IUPAC recommendation, low res. data with interpolation between 362.5 and 410nm
*     mabs = 6: IUPAC recommendation, low res. data forced xs to zero between 362.5 and 410nm
*     mabs = 7: IUPAC recommendation, high res. data by Malicet et al. (1995)
*                                     and xs forced to 0 between 362.5 and 410nm
*                                     (MCM standard)
      spc = "O3"
      xsvers = (/1, 7, 6/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "from Reims group"
      logmsg(2) = "from JPL 2006 evaluation"
      logmsg(3) = "exclusively from Molina and Molina"
      logmsg(4) = "exclusively from Bass et al"
      logmsg(5) =
     & "from IUPAC low res. data interpolated between 362.5 and 410nm"
      logmsg(6) =
     & "from IUPAC low res. data forced to zero between 362.5 and 410nm"
      logmsg(7) = "from IUPAC high res. data by Malicet et al. (1995)"
      CALL set_option(7,spc,"abs",logmsg,xsvers,mabs)

* quantum yield options
* choose quantum yield recommendation:
*    kjpl87:  JPL recommendation 1987                - JPL 87, 90, 92 do not "tail"
*    kjpl92:  JPL recommendations 1990/92 (identical) - still with no "tail"
*    kjpl97:  JPL recommendation 1997, includes tail, similar to Shetter et al.
*    kmich :  Michelsen et al., 1994
*    kshet :  Shetter et al., 1996
*    kjpl00:  JPL 2000
*    kmats:   Matsumi et al., 2002 (TUV standard)
*    kiup:    IUPAC recommendation (last update 2001) (MCM standard)
      kmich = 1
      kjpl87 = 2
      kjpl92 = 3
      kshet = 4
      kjpl97 = 5
      kjpl00 = 6
      kmats = 7
      kiup = 8
      qyvers = (/kmats, kiup, kjpl00/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "from Michelsen et al. (1994)"
      logmsg(2) = "from JPL 1987 evaluation"
      logmsg(3) = "from JPL 1992 evaluation"
      logmsg(4) = "from Shetter et al. (1996)"
      logmsg(5) = "from JPL 1997 evaluation"
      logmsg(6) = "from JPL 2000 evaluation"
      logmsg(7) =
     &      "from Matsumi et al. (2002) (same as JPL 2006 evaluation)"
      logmsg(8) = "from IUPAC recommendation"
      CALL set_option(8,spc,"yld",logmsg,qyvers,myld)


* cross sections
      CALL rdo3xs(mabs, nz,tlev,nw,wl, xs)

******* quantum yield:

* read parameters from JPL'97

      IF (myld .EQ. kjpl97) THEN
        OPEN(UNIT=kin,FILE='DATAJ1/YLD/O3.param_jpl97.yld',STATUS='old')
        READ(kin,*)
        READ(kin,*)
        READ(kin,*)
        n1 = 21
        n2 = n1
        DO i = 1, n1
           READ(kin,*) x1(i), y1(i), y2(i)
           x2(i) = x1(i)
        ENDDO
        CLOSE(kin)

        CALL interpol(x1,y1,kdata,n1,nw,wl,jlabel(j),y1(1),y1(n1),yg1)
        CALL interpol(x2,y2,kdata,n2,nw,wl,jlabel(j),y2(1),y2(n2),yg2)

      ENDIF

* read parameters from Michelsen, H. A., R.J. Salawitch, P. O. Wennber,
* and J. G. Anderson, Geophys. Res. Lett., 21, 2227-2230, 1994.

      IF (myld .EQ. kmich) THEN
        OPEN(UNIT=kin,FILE='DATAJ1/YLD/O3.param.yld',STATUS='old')
        READ(kin,*)
        READ(kin,*)
        READ(kin,*)
        n1 = 21
        n2 = n1
        DO i = 1, n1
           READ(kin,*) x1(i), y1(i), y2(i)
           x2(i) = x1(i)
        ENDDO
        CLOSE(kin)

        CALL interpol(x1,y1,kdata,n1,nw,wl,jlabel(j),y1(1),y1(n1),yg1)
        CALL interpol(x2,y2,kdata,n2,nw,wl,jlabel(j),y2(1),y2(n2),yg2)

      ENDIF

* quantum yield data from
* Shetter et al, J.Geophys.Res., v 101 (D9), pg. 14,631-14,641, June 20, 1996

      IF (myld .EQ. kshet) THEN
        OPEN(UNIT=kin,FILE='DATAJ1/YLD/O3_shetter.yld',STATUS='OLD')
        READ(kin,*) idum, n
        DO i = 1, idum-2
          READ(kin,*)
        ENDDO
        n = n-2
        DO i = 1, n
          READ(kin,*) x1(i),y3(i),y4(i),y1(i),y2(i)
          x2(i) = x1(i)
          x3(i) = x1(i)
          x4(i) = x1(i)
        ENDDO
        DO i = n+1, n+2
           READ(kin,*) x3(i),y3(i),y4(i)
           x4(i) = x3(i)
        ENDDO
        CLOSE(kin)

        n1 = n
        n2 = n
        n3 = n+2
        n4 = n+2

* coefficients for exponential fit:
        CALL interpol(x1,y1,kdata,n1,nw,wl,jlabel(j),y1(1),0.,yg1)
        CALL interpol(x2,y2,kdata,n2,nw,wl,jlabel(j),y2(1),0.,yg2)


* phi data at 298 and 230 K
        CALL interpol(x3,y3,kdata,n3,nw,wl,jlabel(j),y3(1),0.,yg3)
        CALL interpol(x4,y4,kdata,n4,nw,wl,jlabel(j),y4(1),0.,yg4)
      ENDIF


* compute cross sections and yields at different wavelengths, altitudes:

      DO 10 iw = 1, nw-1

         DO 20 i = 1, nz

* quantum yields
* coefficients from jpl 87:

             IF (myld .EQ. kjpl87) THEN
               tau = tlev(i) - 230.
               tau2 = tau*tau
               tau3 = tau2*tau
               xl = wc(iw)
               xl0 = 308.2 + 4.4871e-2*tau + 6.938e-5*tau2 -
     >               2.5452e-6*tau3
               a = 0.9*(0.369 + 2.85e-4*tau + 1.28e-5*tau2 +
     >                  2.57e-8*tau3)
               b     = -0.575 + 5.59e-3*tau - 1.439e-5*tau2 -
     >                  3.27e-8*tau3
               c = 0.9*(0.518 + 9.87e-4*tau - 3.94e-5*tau2 +
     >                  3.91e-7*tau3)
               qy1d = a*atan(b*(xl-xl0)) + c
               qy1d = amax1(0.,qy1d)
               qy1d = amin1(0.9,qy1d)
             ENDIF

* from jpl90, jpl92:
* (caution: error in JPL92 for first term of a3)

             IF (myld .EQ. kjpl92) THEN
               tau = 298. - tlev(i)
               tau2 = tau*tau
               xl0 = wc(iw) - 305.
               a0 =   .94932   - 1.7039e-4*tau + 1.4072E-6*tau2
               a1 = -2.4052e-2 + 1.0479e-3*tau - 1.0655e-5*tau2
               a2 =  1.8771e-2 - 3.6401e-4*tau - 1.8587e-5*tau2
               a3 = -1.4540e-2 - 4.7787e-5*tau + 8.1277e-6*tau2
               a4 =  2.3287e-3 + 1.9891e-5*tau - 1.1801e-6*tau2
               a5 = -1.4471e-4 - 1.7188e-6*tau + 7.2661e-8*tau2
               a6 =  3.1830e-6 + 4.6209e-8*tau - 1.6266e-9*tau2
               qy1d = a0 + a1*xl0 + a2*(xl0)**2 + a3*(xl0)**3 +
     >                a4*(xl0)**4 + a5*(xl0)**5 + a6*(xl0)**6
               IF (wc(iw) .LT. 305.) qy1d = 0.95
               IF (wc(iw) .GT. 320.) qy1d = 0.
               IF (qy1d .LT. 0.02) qy1d = 0.
             ENDIF

* from JPL'97

           IF (myld .EQ. kjpl97) THEN
             IF (wc(iw) .LT. 271.) THEN
                qy1d = 0.87
             ELSE IF (wc(iw) .GE. 271. .AND. wc(iw) .LT. 290.) THEN
                qy1d = 0.87 + (wc(iw)-271.)*(.95-.87)/(290.-271.)
             ELSE IF (wc(iw) .GE. 290. .AND. wc(iw) .LT. 305.) THEN
                qy1d = 0.95
             ELSE IF (wc(iw) .GE. 305. .AND. wc(iw) .LE. 325.) THEN
                qy1d = yg1(iw) * EXP ( -yg2(iw) /tlev(i) )
             ELSE
                qy1d = 0.
             ENDIF
           ENDIF

* from Michelsen, H. A., R.J. Salawitch, P. O. Wennber, and J. G. Anderson
* Geophys. Res. Lett., 21, 2227-2230, 1994.

           IF (myld .EQ. kmich) THEN
             IF (wc(iw) .LT. 271.) THEN
                qy1d = 0.87
             ELSE IF (wc(iw) .GE. 271. .AND. wc(iw) .LT. 305.) THEN
                qy1d = 1.98 - 301./wc(iw)
             ELSE IF (wc(iw) .GE. 305. .AND. wc(iw) .LE. 325.) THEN
                qy1d = yg1(iw) * EXP (-yg2(iw) /(0.6951*tlev(i)))
             ELSE
                qy1d = 0.
             ENDIF
           ENDIF

* Shetter et al.:
* phi = A * exp(-B/T), A and B are based on meas. at 298 and 230 K
* do linear interpolation between phi(298) and phi(230) for wavelengths > 321
* as phi(230)=0. for those wavelengths, so there are no A and B factors

           IF (myld .EQ. kshet) THEN
             IF (wl(iw+1) .LE. 321.) THEN
               qy1d = yg1(iw) * EXP(-1. * yg2(iw)/tlev(i))
             ELSE
               qy1d = (yg3(iw) - yg4(iw))/(298.-230.) * (tlev(i)-230.) +
     >                 yg4(iw)
             ENDIF
           ENDIF

* JPL 2000:

           IF (myld == kjpl00) THEN
              qy1d = fo3qy(wc(iw),tlev(i))
           ENDIF

* Matsumi et al. / IUPAC:

           IF (myld == kmats .or. myld == kiup) THEN
              qy1d = fo3qy2(wc(iw),tlev(i))
           ENDIF

* compute product

           sq(j-1,i,iw) = qy1d*xs(i,iw)
           qy3p = 1.0 - qy1d
           sq(j,i,iw) = qy3p*xs(i,iw)

 20     CONTINUE
 10   CONTINUE

      RETURN
      END

*=============================================================================*

      SUBROUTINE r02(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! NO2

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for NO2            =*
*=  photolysis:                                                              =*
*=         NO2 + hv -> NO + O(3P)                                            =*
*=  Cross section from IUPAC recommendations, last update 16.7.2001          =*
*=  Quantum yield from IUPAC recommendations, last update 16.7.2001          =*
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

      INTEGER nw, iw
      REAL wl(kw), wc(kw)

      INTEGER nz, iz

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

      INTEGER n2
      REAL x1(kdata), x2(kdata)
      REAL y1(kdata), y2(kdata)

* local

      REAL yg1(kw), yg2(kw)
      REAL no2xs(kz,kw), qy(kz,kw)
      INTEGER i, n

      INTEGER mabs, myld, xsvers(3), qyvers(3)
      CHARACTER(llog) :: logmsg(nlog)
      CHARACTER(lspc) :: spc


***************** NO2 photolysis *************************

* Setting photolysis index and cross section/quantum yield options
      j = j + 1
      jlabel(j) = 'NO2 -> NO + O(3P)'

* cross section options

      spc = "NO2"
      xsvers = (/4, 6, 6/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "from Davidson et al. (1988), indepedent of T"
      logmsg(2) = "from JPL 1994 - 2002 evaluations"
      logmsg(3) = "from Harder et al"
      logmsg(4) =
     &          "from JPL 2006 interpolating between midpoints of bins"
      logmsg(5) = "from JPL 2006, bin-to-bin interpolation"
      logmsg(6) = "from IUPAC recommendation"
      CALL set_option(6,spc,"abs",logmsg,xsvers,mabs)

* quantum yield options
      qyvers = (/2, 3, 3/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "from JPL 2002 evaluation"
      logmsg(2) = "from JPL 2006 evaluation"
      logmsg(3) = "from IUPAC recommendation"
      CALL set_option(3,spc,"yld",logmsg,qyvers,myld)


* cross sections
      IF (mabs .EQ. 1) CALL no2xs_d(nz,tlev,nw,wl, no2xs)
      IF (mabs .EQ. 2) CALL no2xs_jpl94(nz,tlev,nw,wl, no2xs)
      IF (mabs .EQ. 3) CALL no2xs_har(nz,tlev,nw,wl, no2xs)
      IF (mabs .EQ. 4) CALL no2xs_jpl06a(nz,tlev,nw,wl, no2xs)
      IF (mabs .EQ. 5) CALL no2xs_jpl06b(nz,tlev,nw,wl, no2xs)
      IF (mabs .EQ. 6) CALL no2xs_iup(nz,tlev,nw,wl, no2xs)

* quantum yields
*     myld = 1   NO2_calvert.yld  (same as JPL2002)
*     myld = 2   NO2_jpl11.yld (same as jpl2006) (TUV standard)
*     myld = 3   NO2_IUPAC.yld (IUPAC recommendation) (MCM standard)

      IF (myld .EQ. 1) THEN

* from Gardiner, Sperry, and Calvert

         OPEN(UNIT=kin,FILE='DATAJ1/YLD/NO2_calvert.yld',STATUS='old')
         DO i = 1, 8
            READ(kin,*)
         ENDDO
         n = 66
         DO i = 1, n
            READ(kin,*) x1(i),y1(i)
         ENDDO
         CLOSE(kin)
        CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),y1(1),0.,yg1)


         do iw = 1, nw - 1
            do iz = 1, nz
               qy(iz,iw) = yg1(iw)
            enddo
         enddo

      ELSE IF(myld. EQ. 2) THEN

* from jpl 2011

         OPEN(UNIT=kin,FILE='DATAJ1/YLD/NO2_jpl11.yld',STATUS='old')
         DO i = 1, 2
            READ(kin,*)
         ENDDO
         n = 25
         n2 = n
         DO i = 1, n
            READ(kin,*) x1(i),y1(i),y2(i)
            x2(i) = x1(i)
         ENDDO
         CLOSE(kin)

         CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),y1(1),0.,yg1)
         CALL interpol(x2,y2,kdata,n2,nw,wl,jlabel(j),y2(1),0.,yg2)

         DO iw = 1, nw - 1
            DO iz = 1, nz
               qy(iz,iw) = yg1(iw) +
     $              (yg1(iw)-yg2(iw)) * (tlev(iz)-298)/50.
               qy(iz,iw) = amax1(qy(iz,iw),0.)
            ENDDO
         ENDDO

      ELSE IF(myld. EQ. 3) THEN

* from IUPAC recommendation

         OPEN(UNIT=kin,FILE='DATAJ1/YLD/NO2_IUPAC.yld',STATUS='old')
         DO i = 1, 5
            READ(kin,*)
         ENDDO
         n = 20
         n2 = n
         DO i = 1, n
            READ(kin,*) x1(i),y1(i),y2(i)
            x2(i) = x1(i)
         ENDDO
         CLOSE(kin)

         CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),y1(1),0.,yg1)
         CALL interpol(x2,y2,kdata,n2,nw,wl,jlabel(j),y2(1),0.,yg2)

         DO iw = 1, nw - 1
            DO iz = 1, nz
               qy(iz,iw) = yg1(iw) +
     $              (yg1(iw)-yg2(iw)) * (tlev(iz)-298.)/50.
               qy(iz,iw) = amax1(qy(iz,iw),0.)
            ENDDO
         ENDDO

      ENDIF

* combine

      DO iw = 1, nw - 1
         DO iz = 1, nz
            sq(j,iz,iw) = no2xs(iz,iw)*qy(iz,iw)
         ENDDO
      ENDDO

      RETURN
      END

*=============================================================================*

      SUBROUTINE r03(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! NO3

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (absorptioon cross section) x (quantum yield) for    =*
*=  both channels of NO3 photolysis:                                         =*
*=          (a) NO3 + hv -> NO + O2                                          =*
*=          (b) NO3 + hv -> NO2 + O(3P)                                      =*
*=  Cross section from IUPAC recommendation                                  =*
*=  Quantum yield from Madronich (1988)                                      =*
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
      PARAMETER(kdata=350)

      REAL x(kdata), x1(kdata), x2(kdata)
      REAL y1(kdata), y2(kdata)
      real q1_298(kdata), q1_230(kdata), q1_190(kdata)
      real q2_298(kdata), q2_230(kdata), q2_190(kdata)


* local

      REAL yg(kw), yg1(kw), yg2(kw)
      REAL qy, qy1, qy2
      real yg1_298(kw), yg1_230(kw), yg1_190(kw)
      real yg2_298(kw), yg2_230(kw), yg2_190(kw)
      REAL xs(kz,kw)

      INTEGER irow, icol
      INTEGER i, iw, iz, n, n1, n2, idum

      INTEGER mabs, myld, xsvers(3), qyvers(3)
      CHARACTER(llog) :: logmsg(nlog)
      CHARACTER(lspc) :: spc


***************** NO3 photolysis *************************

* Setting photolysis index and cross section/quantum yield options
      j = j + 1
      jlabel(j) = 'NO3 -> NO + O2'
      j = j + 1
      jlabel(j) = 'NO3 -> NO2 + O(3P)'

* cross section options
      spc = "NO3"
      xsvers = (/3, 4, 2/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "from Graham and Johnston (1978)"
      logmsg(2) = "from JPL 1994 evaluation"
      logmsg(3) = "from JPL 2011 evaluation"
      logmsg(4) = "from IUPAC recommendation"
      CALL set_option(4,spc,"abs",logmsg,xsvers,mabs)

* quantum yield options
      qyvers = (/2, 3, 1/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "from Madronich (1988)"
      logmsg(2) = "from JPL 2011 evaluation"
      logmsg(3) = "from IUPAC & JPL recommendation"
      CALL set_option(3,spc,"yld",logmsg,qyvers,myld)
* myld = 1  from Madronich (1988) see CEC NO3 book.
* myld = 2  from JPL-2011 (TUV standard)
* myld = 3  below 587nm Phi(2) = 1 (according to IUPAC), above JPL recomm.
*           (MCM standard)


* cross sections

      IF(mabs. eq. 1) then

*     measurements of Graham and Johnston 1978

         OPEN(UNIT=kin,FILE='DATAJ1/ABS/NO3_gj78.abs',STATUS='old')
         DO i = 1, 9
            READ(kin,*)
         ENDDO
         n = 305
         DO irow = 1, 30
            READ(kin,*) ( y1(10*(irow-1) + icol), icol =  1, 10 )
         ENDDO
         READ(kin,*) ( y1(300 + icol), icol = 1, 5 )
         CLOSE (kin)
         DO i = 1, n
            y1(i) =  y1(i) * 1.E-19
            x1(i) = 400. + 1.*FLOAT(i-1)
         ENDDO

         CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)

         DO iw = 1, nw-1
          DO iz = 1, nz
            xs(iz,iw) = yg(iw)
          ENDDO
         ENDDO

      ELSEIF(mabs .EQ. 2) THEN

*     cross section from JPL94:

         OPEN(UNIT=kin,FILE='DATAJ1/ABS/NO3_jpl94.abs',STATUS='old')
         READ(kin,*) idum, n
         DO i = 1, idum-2
            READ(kin,*)
         ENDDO
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
            y1(i) = y1(i)*1E-20
         ENDDO
         CLOSE (kin)
         CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg1)

* (use JPL94 for wavelengths longer than 600 nm)

         DO iw = 1, nw-1
          DO iz = 1, nz
            xs(iz,iw) = yg1(iw)
          ENDDO
         ENDDO

      ELSEIF(MABS .EQ. 3) THEN

* cross sections from JPL2011

         OPEN(UNIT=kin,FILE='DATAJ1/ABS/NO3_jpl11.abs',STATUS='old')
         DO i = 1, 6
            READ(kin,*)
         ENDDO
         DO i = 1, 289
            READ(kin,*) x1(i), y1(i)
            y1(i) = y1(i)*1E-20
         ENDDO
         CLOSE (kin)

         n = 289
         CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg1)
         DO iw = 1, nw-1
          DO iz = 1, nz
            xs(iz,iw) = yg1(iw)
          ENDDO
         ENDDO

      ELSEIF(mabs .EQ. 4) THEN

* IUPAC recommendation

        OPEN(UNIT=kin,FILE='DATAJ1/ABS/NO3_IUPAC.abs',STATUS='old')
        DO i = 1, 5
         read(kin,*)
        ENDDO
        n  = 291
        n1 = n
        n2 = n
        DO i = 1, n
           READ(kin,*) x1(i), y2(i), y1(i)
           x2(i) = x1(i)
           y1(i) = y1(i) * 1.e-19
           y2(i) = y2(i) * 1.e-19
        ENDDO
        CLOSE (kin)


        CALL interpol(x1,y1,kdata,n1,nw,wl,jlabel(j),0.,0.,yg1)
        CALL interpol(x2,y2,kdata,n2,nw,wl,jlabel(j),0.,0.,yg2)
        DO iw = 1, nw-1
          DO iz = 1, nz

            xs(iz,iw) = yg1(iw) +
     $           (yg2(iw) - yg1(iw))*(tlev(iz) - 230.)/(298. - 230.)
            IF (tlev(iz) .LE. 230.) xs(iz,iw) = yg1(iw)
            IF (tlev(iz) .GE. 298.) xs(iz,iw) = yg2(iw)

          ENDDO
         ENDDO
      ENDIF

* quantum yield:

      if (myld .eq. 1) then

* for   NO3 ->NO+O2

         DO iw = 1, nw - 1
            IF (wc(iw).LT.584.) THEN
               qy = 0.
            ELSEIF (wc(iw).GE.640.) THEN
               qy = 0.
            ELSEIF (wc(iw).GE.595.) THEN
               qy = 0.35*(1.-(wc(iw)-595.)/45.)
            ELSE
               qy = 0.35*(wc(iw)-584.)/11.
            ENDIF
            DO i = 1, nz
               sq(j-1,i,iw) = xs(i,iw)*qy
            ENDDO
         ENDDO

* for  NO3 ->NO2+O

         DO iw = 1, nw - 1
            IF (wc(iw).LT.584.) THEN
               qy = 1.
            ELSEIF (wc(iw).GT.640.) THEN
               qy = 0.
            ELSEIF (wc(iw).GT.595.) THEN
               qy = 0.65*(1-(wc(iw)-595.)/45.)
            ELSE
               qy = 1.-0.35*(wc(iw)-584.)/11.
            ENDIF
            DO i = 1, nz
               sq(j,i,iw) = xs(i,iw)*qy
            ENDDO
         ENDDO

* yields from JPL2011:

      ELSEIF(myld .EQ. 2 .OR. myld .EQ. 3) THEN

         open(unit=kin,file='DATAJ1/YLD/NO3_jpl2011.qy',status='old')
         do i = 1, 5
            read(kin,*)
         enddo
         do i = 1, 56
            read(kin,*) x(i), q1_298(i), q1_230(i), q1_190(i),
     $           q2_298(i), q2_230(i), q2_190(i)

            q1_298(i) = q1_298(i)/1000.
            q1_230(i) = q1_230(i)/1000.
            q1_190(i) = q1_190(i)/1000.
            q2_298(i) = q2_298(i)/1000.
            q2_230(i) = q2_230(i)/1000.
            q2_190(i) = q2_190(i)/1000.

         enddo
         close(kin)

         n = 56
         do i = 1, n
            x1(i) = x(i)
         enddo
         CALL interpol(x1,q1_298,kdata,n,nw,wl,jlabel(j),0.,0.,yg1_298)

         n = 56
         do i = 1, n
            x1(i) = x(i)
         enddo
         CALL interpol(x1,q1_230,kdata,n,nw,wl,jlabel(j),0.,0.,yg1_230)
         n = 56
         do i = 1, n
            x1(i) = x(i)
         enddo
         CALL interpol(x1,q1_190,kdata,n,nw,wl,jlabel(j),0.,0.,yg1_190)

         n = 56
         do i = 1, n
            x1(i) = x(i)
         enddo
         CALL interpol(x1,q2_298,kdata,n,nw,wl,jlabel(j),1.,0.,yg2_298)

         n = 56
         do i = 1, n
            x1(i) = x(i)
         enddo
         CALL interpol(x1,q2_230,kdata,n,nw,wl,jlabel(j),1.,0.,yg2_230)

         n = 56
         do i = 1, n
            x1(i) = x(i)
         enddo
         CALL interpol(x1,q2_190,kdata,n,nw,wl,jlabel(j),1.,0.,yg2_190)

* compute T-dependent quantum yields

         DO iw = 1, nw-1
            DO iz = 1, nz

               if(tlev(iz) .le. 190.) then

                  qy1 = yg1_190(iw)
                  qy2 = yg2_190(iw)

               elseif(tlev(iz) .gt. 190. .and. tlev(iz) .le. 230.) then

                  qy1 = yg1_190(iw) + (yg1_230(iw) - yg1_190(iw))*
     $                 (tlev(iz) - 190.)/(230.-190.)
                  qy2 = yg2_190(iw) + (yg2_230(iw) - yg2_190(iw))*
     $                 (tlev(iz) - 190.)/(230.-190.)

               elseif(tlev(iz) .gt. 230. .and. tlev(iz) .le. 298.) then

                  qy1 = yg1_230(iw) + (yg1_298(iw) - yg1_230(iw))*
     $                 (tlev(iz) - 230.)/(298.-230.)

                  qy2 = yg2_230(iw) + (yg2_298(iw) - yg2_230(iw))*
     $                 (tlev(iz) - 230.)/(298.-230.)

               elseif(tlev(iz) .gt. 298.) then

                  qy1 = yg1_298(iw)
                  qy2 = yg2_298(iw)

               endif

               IF(myld .EQ. 3 .AND. wl(iw) .LT. 587.) THEN
* Use qy2 = 1 for wavelength below 587nm (as recommended in IUPAC),
* JPL recommendation else
                 qy1 = 0
                 qy2 = 1
               ENDIF

               sq(j-1, iz, iw) = qy1 * xs(iz,iw)
               sq(j,   iz, iw) = qy2 * xs(iz,iw)

            ENDDO
         ENDDO

      ENDIF

      RETURN
      END

*=============================================================================*

      SUBROUTINE r04(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! N2O5

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product of (cross section) x (quantum yiels) for N2O5 photolysis =*
*=  reactions:                                                               =*
*=       (a) N2O5 + hv -> NO3 + NO + O(3P)                                   =*
*=       (b) N2O5 + hv -> NO3 + NO2                                          =*
*=  Cross section from JPL2011: use tabulated values for 300K, correct for   =*
*=  temperature.
*=  Quantum yield: Analysis of data in JPL94 (->DATAJ1/YLD/N2O5.qy)          =*
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

      REAL x1(kdata), x2(kdata)
      REAL y1(kdata), A(kdata), B(kdata)
      INTEGER n1, n2

* local

      REAL yg1(kw), yg2(kw)
      INTEGER i, iz, iw
      REAL t, xs, dum

***************** N2O5 photolysis *************************

      j = j + 1
      jlabel(j) = 'N2O5 -> NO3 + NO + O(3P)'
      j = j + 1
      jlabel(j) = 'N2O5 -> NO3 + NO2'


* cross section from jpl2011, at 300 K

      OPEN(UNIT=kin,FILE='DATAJ1/ABS/N2O5_jpl11.abs',STATUS='old')
      DO i = 1, 4
         READ(kin,*)
      ENDDO
      n1 = 103
      DO i = 1, n1
         READ(kin,*) x1(i), y1(i)
         y1(i) = y1(i) * 1.E-20
      ENDDO

* read temperature dependence coefficients:

      DO i = 1, 4
         READ(kin,*)
      ENDDO
      n2 = 8
      DO i = 1, n2
         READ(kin,*) x2(i), A(i), B(i)
      ENDDO
      CLOSE(kin)

      CALL interpol(x1,y1,kdata,n1,nw,wl,jlabel(j),0.,0.,yg1)
      CALL interpol(x2,B,kdata,n2,nw,wl,jlabel(j),0.,0.,yg2)


* combine with quantum yields
      DO iw = 1, nw - 1
         DO iz = 1, nz

* temperature dependence only valid for 233 - 295 K.  Extend to 300.

            t = MAX(233.,MIN(tlev(iz),300.))

* Apply temperature correction to 300K values. Do not use A-coefficients
* because they are inconsistent with the values at 300K.

            dum = 1000.*yg2(iw)*((1./t) - (1./300.))
            xs = yg1(iw) * 10.**(dum)

* quantum yield = 1 for NO2 + NO3, zero for other channels

            sq(j-1, iz, iw) = 0. * xs
            sq(j  , iz, iw) = 1. * xs

         ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r05(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! HNO2

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for HNO2 photolysis=*
*=     HNO2 + hv -> NO + OH                                                  =*
*=  Cross section:  from IUPAC recommendation, last update: 16.7.2001        =*
*=  Quantum yield:  unity (from IUPAC recommendation, last update: 16.7.2001)=*
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
*=  EDIT HISTORY:                                                            =*
*=  05/98  Original, adapted from former JSPEC1 subroutine                   =*
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

      REAL x1(kdata)
      REAL y1(kdata)

* local

      REAL yg(kw)
      INTEGER i, iw, n

      INTEGER mabs, xsvers(3)
      CHARACTER(llog) :: logmsg(nlog)
      CHARACTER(lspc) :: spc


***************** HNO2 photolysis *************************

* Setting photolysis index and cross section options
      j = j + 1
      jlabel(j) = 'HNO2 -> OH + NO'

      spc = "HNO2"
      xsvers = (/2, 3, 4/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "from JPL 2006 or older recommendation"
      logmsg(2) = "from JPL 2011 recommendation"
      logmsg(3) = "from IUPAC forced to 0 between 280 and 295nm"
      logmsg(4) = "from IUPAC TUV interpolated between 280 and 295nm"
      CALL set_option(4,spc,"abs",logmsg,xsvers,mabs)


* cross sections

      IF(mabs .eq. 1) then

         OPEN(UNIT=kin,FILE='DATAJ1/ABS/HONO_jpl92.abs',STATUS='old')
         DO i = 1, 13
            READ(kin,*)
         ENDDO
         n = 91
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
            y1(i) = y1(i) * 1.E-20
         ENDDO
         CLOSE (kin)

         CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)

       ELSEIF(mabs .eq. 2) then

         OPEN(UNIT=kin,FILE='DATAJ1/ABS/HONO_jpl11.abs',STATUS='old')
         DO i = 1, 3
            READ(kin,*)
         ENDDO
         n = 192
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
            y1(i) = y1(i) * 1.E-20
         ENDDO
         CLOSE (kin)

         CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)

       ELSEIF(mabs .eq. 3) then

* IUPAC recommendation, forced to 0 between 280 and 295nm

         OPEN(UNIT=kin,FILE='DATAJ1/ABS/HONO_IUPACa.abs',STATUS='old')
         DO i = 1, 6
            READ(kin,*)
         ENDDO
         n = 42
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
            y1(i) = y1(i) * 1.E-20
         ENDDO
         CLOSE (kin)

         CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)

       ELSEIF(mabs .eq. 4) then

* IUPAC recommendation, TUV interpolated between 280 and 295nm

         OPEN(UNIT=kin,FILE='DATAJ1/ABS/HONO_IUPACb.abs',STATUS='old')
         DO i = 1, 7
            READ(kin,*)
         ENDDO
         n = 38
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
            y1(i) = y1(i) * 1.E-20
         ENDDO
         CLOSE (kin)

         CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)

      ENDIF

* quantum yield = 1

!      qy = 1.
      DO iw = 1, nw - 1
         DO i = 1, nz
            sq(j,i,iw) = yg(iw)!*qy
         ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r06(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! HNO3

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product of (cross section) x (quantum yield) for HNO3 photolysis =*
*=        HNO3 + hv -> OH + NO2                                              =*
*=  Cross section: IUPAC recommendation (Burkholder et al., 1993)            =*
*=  Quantum yield: IUPAC rcommendationn or unity                             =*
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
      PARAMETER(kdata=100)

      INTEGER n, n1, n2
      REAL x1(kdata), x2(kdata)
      REAL y1(kdata), y2(kdata)

* local

      REAL yg1(kw), yg2(kw), yg(kz,kw)
      REAL qy
      INTEGER i, iw

      INTEGER mabs, myld, xsvers(3), qyvers(3)
      CHARACTER(llog) :: logmsg(nlog)
      CHARACTER(lspc) :: spc


***************** HNO3 photolysis *************************

* Setting photolysis index and cross section/quantum yield options
      j = j + 1
      jlabel(j) = 'HNO3 -> OH + NO2'
      j = j + 1
      jlabel(j) = 'HNO3 -> OH + NO + O(3P)'
      j = j + 1
      jlabel(j) = 'HNO3 -> OH + NO + O(1D)'
      j = j + 1
      jlabel(j) = 'HNO3 -> HONO + O(3P)'
      j = j + 1
      jlabel(j) = 'HNO3 -> HONO + O(1D)'

* cross section options
      spc = "HNO3"
      xsvers = (/2, 3, 1/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "from JPL 1985 recommendation"
      logmsg(2) = "from TUV provided Burkholder data"
      logmsg(3) = "from IUPAC recommended Burkholder data"
      CALL set_option(3,spc,"abs",logmsg,xsvers,mabs)

* quantum yield options
      qyvers = (/1, 2, 1/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "set to unity"
      logmsg(2) = "from IUPAC recommendation"
      CALL set_option(2,spc,"yld",logmsg,qyvers,myld)


* cross sections

      IF(mabs==1) THEN
* cross section from JPL85

        OPEN(UNIT=kin,FILE='DATAJ1/ABS/HNO3.abs',STATUS='old')
        DO i = 1, 9
          READ(kin,*)
        ENDDO
        n = 29
        DO i = 1, n
          READ(kin,*) x1(i), y1(i)
          y1(i) = y1(i) * 1.E-20
        ENDDO
        CLOSE (kin)
        CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg1)

       ELSEIF(mabs==2) THEN

* HNO3 cross section parameters from Burkholder et al. 1993

        OPEN(UNIT=kin,FILE='DATAJ1/ABS/HNO3_burk.abs',STATUS='old')
        DO i = 1, 6
          READ(kin,*)
        END DO
        n1 =  83
        n2 = n1
        DO i = 1, n1
          READ(kin,*) y1(i), y2(i)
          x1(i) = 184. + i*2.
          x2(i) = x1(i)
        END DO
        CLOSE (kin)
        CALL interpol(x1,y1,kdata,n1,nw,wl,jlabel(j),0.,0.,yg1)
        CALL interpol(x2,y2,kdata,n2,nw,wl,jlabel(j),y2(1),y2(n2),yg2)

       ELSEIF(mabs==3) THEN

* HNO3 IUPAC recommended cross section parameters from Burkholder et al. 1993

        OPEN(UNIT=kin,FILE='DATAJ1/ABS/HNO3_IUPAC.abs',STATUS='old')
        DO i = 1, 7
          READ(kin,*)
        END DO
        n1 = 33
        n2 = n1
        DO i = 1, n1
          READ(kin,*) x1(i), y1(i), y2(i)
          x2(i) = x1(i)
        END DO
        CLOSE (kin)

        CALL interpol(x1,y1,kdata,n1,nw,wl,jlabel(j),0.,0.,yg1)
        CALL interpol(x2,y2,kdata,n2,nw,wl,jlabel(j),y2(1),y2(n2),yg2)
      ENDIF

      IF(mabs==1) THEN
        DO iw = 1, nw - 1
         DO i = 1, nz
           sq(j-4,i,iw) = yg1(iw)
         ENDDO
        ENDDO
       ELSEIF(mabs==2 .OR. mabs==3) THEN
* Introduce Temperature dependencies to Burkholder data and correct for units
        DO iw = 1, nw - 1
         DO i = 1, nz
           yg(i,iw) = yg1(iw) * 1.E-20
     &                  * exp( yg2(iw)/1.e3*(tlev(i)-298.))
         ENDDO
        ENDDO
      ENDIF

* myld = 1: unity (TUV standard)
* myld = 2: IUPAC (assuming 5 possible channels, MCM standard):
!                 • OH + NO2
!                 • OH + NO + O(3P)
!                 • OH + NO + O(1D)
!                 • HONO + O(3P)
!                 • HONO + O(1D)



      IF(myld==1) THEN
* quantum yield = 1
        DO iw = 1, nw - 1
         DO i = 1, nz
           sq(j-4,i,iw) = yg(i,iw)
         ENDDO
        ENDDO
         sq(j-3,:,:) = 0.
         sq(j-2,:,:) = 0.
         sq(j-1,:,:) = 0.
         sq(j,:,:)   = 0.


       ELSEIF(myld==2) THEN
* for   HNO3 -> OH + NO2

         DO iw = 1, nw - 1
            IF (wc(iw).GE.248.) THEN
               qy = 0.97
            ELSEIF (wc(iw).LT.248. .AND. wc(iw).GE.222.) THEN
               qy = 0.8+(0.97-0.8)*((wc(iw)-222.)/(248.-222.))
            ELSEIF (wc(iw).LT.222. .AND. wc(iw).GE.193.) THEN
               qy = 0.19+(0.8-0.19)*((wc(iw)-193.)/(222.-193.))
            ELSE
               qy = 0.19
            ENDIF
            DO i = 1, nz
               sq(j-4,i,iw) = yg(i,iw)*qy
            ENDDO
         ENDDO


* for   HNO3 -> OH + NO + O(3P)

         DO iw = 1, nw - 1
            IF (wc(iw).GE.248.) THEN
               qy = 0.027
            ELSEIF (wc(iw).LT.248. .AND. wc(iw).GE.222.) THEN
               qy = 0.06+(0.027-0.06)*((wc(iw)-222.)/(248.-222.))
            ELSEIF (wc(iw).LT.222. .AND. wc(iw).GE.193.) THEN
               qy = 0.438+(0.06-0.438)*((wc(iw)-193.)/(222.-193.))
            ELSE
               qy = 0.438
            ENDIF
            DO i = 1, nz
               sq(j-3,i,iw) = yg(i,iw)*qy
            ENDDO
         ENDDO


* for   HNO3 -> OH + NO + O(1D)

         DO iw = 1, nw - 1
            IF (wc(iw).GE.248.) THEN
               qy = 0.003
            ELSEIF (wc(iw).LT.248. .AND. wc(iw).GE.222.) THEN
               qy = 0.04+(0.003-0.04)*((wc(iw)-222.)/(248.-222.))
            ELSEIF (wc(iw).LT.222. .AND. wc(iw).GE.193.) THEN
               qy = 0.232+(0.04-0.232)*((wc(iw)-193.)/(222.-193.))
            ELSE
               qy = 0.232
            ENDIF
            DO i = 1, nz
               sq(j-2,i,iw) = yg(i,iw)*qy
            ENDDO
         ENDDO


* for   HNO3 -> HONO + O(3P)

         DO iw = 1, nw - 1
            IF (wc(iw).GE.248.) THEN
               qy = 0.0
            ELSEIF (wc(iw).LT.248. .AND. wc(iw).GE.222.) THEN
               qy = 0.06+(0.0-0.06)*((wc(iw)-222.)/(248.-222.))
            ELSEIF (wc(iw).LT.222. .AND. wc(iw).GE.193.) THEN
               qy = 0.092+(0.06-0.092)*((wc(iw)-193.)/(222.-193.))
            ELSE
               qy = 0.092
            ENDIF
            DO i = 1, nz
               sq(j-1,i,iw) = yg(i,iw)*qy
            ENDDO
         ENDDO


* for   HNO3 -> HONO + O(1D)

         DO iw = 1, nw - 1
            IF (wc(iw).GE.248.) THEN
               qy = 0.0
            ELSEIF (wc(iw).LT.248. .AND. wc(iw).GE.222.) THEN
               qy = 0.04+(0.0-0.04)*((wc(iw)-222.)/(248.-222.))
            ELSEIF (wc(iw).LT.222. .AND. wc(iw).GE.193.) THEN
               qy = 0.048+(0.04-0.048)*((wc(iw)-193.)/(222.-193.))
            ELSE
               qy = 0.048
            ENDIF
            DO i = 1, nz
               sq(j,i,iw) = yg(i,iw)*qy
            ENDDO
         ENDDO

      ENDIF

      END

*=============================================================================*

      SUBROUTINE r07(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! HNO4

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product of (cross section) x (quantum yield) for HNO4 photolysis =*
*=       HNO4 + hv -> HO2 + NO2                                              =*
*=  Cross section:  from IUPAC                                               =*
*=  Quantum yield:  from IUPAC                                               =*
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
      PARAMETER(kdata=100)

      REAL x1(kdata)
      REAL y1(kdata)

C* local

      REAL yg(kw)
      REAL qy1, qy2
      INTEGER i, iw, n

      INTEGER mabs, myld, xsvers(3), qyvers(3)
      CHARACTER(llog) :: logmsg(nlog)
      CHARACTER(lspc) :: spc


***************** HNO4 photolysis *************************

* Setting photolysis index and cross section/quantum yield options
      j = j + 1
      jlabel(j) = 'HNO4 -> HO2 + NO2'
      j = j + 1
      jlabel(j) = 'HNO4 -> OH + NO3'

* cross section options
      spc = "HNO4"
      xsvers = (/1, 2, 2/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "from JPL 2011"
      logmsg(2) = "ffrom IUPAC"
      CALL set_option(2,spc,"abs",logmsg,xsvers,mabs)

* quantum yield options
      qyvers = (/1, 2, 2/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "set to unity"
      logmsg(2) = "from IUPAC"
      CALL set_option(2,spc,"yld",logmsg,qyvers,myld)


* cross sections

      IF(mabs==1) THEN

*--------------- TUV standard version (start) --------------------------
      OPEN(UNIT=kin,FILE='DATAJ1/ABS/HNO4_jpl11.abs',STATUS='old')
      DO i = 1, 2
         READ(kin,*)
      ENDDO
      n = 54
      DO i = 1, n
         READ(kin,*) x1(i), y1(i)
         y1(i) = y1(i) * 1.E-20
      ENDDO
      CLOSE (kin)
*--------------- TUV standard version (end) ----------------------------

      ELSEIF(mabs==2) THEN

*--------------- MCM standard version (start) --------------------------
      OPEN(UNIT=kin,FILE='DATAJ1/ABS/HNO4_IUPAC.abs',STATUS='old')
      DO i = 1, 9
         READ(kin,*)
      ENDDO
      n = 29
      DO i = 1, n
         READ(kin,*) x1(i), y1(i)
      ENDDO
      CLOSE (kin)
*--------------- MCM standard version (end) ----------------------------
      ENDIF

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)


* quantum yields

      IF(myld == 1) THEN
        qy1 = 1.
        qy2 = 0.
       ELSEIF(myld==2) THEN
        qy1 = 0.59
        qy2 = 0.41
      ENDIF
      DO iw = 1, nw - 1
         DO i = 1, nz
            sq(j-1,i,iw) = yg(iw)*qy1
            sq(j,i,iw)   = yg(iw)*qy2
         ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r08(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! H2O2

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product of (cross section) x (quantum yield) for H2O2 photolysis =*
*=         H2O2 + hv -> 2 OH                                                 =*
*=  Cross section:  From IUPAC, tabulated values between 195 and 350nm       =*
*=  Quantum yield:  IUPAC, Phi = 1 for lambda > 220nm, Phi = 0.85 else       =*
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
      REAL sq(kj,kz,kw), absq(kz,kw)

* input/output:
      INTEGER j

* data arrays

      INTEGER kdata
      PARAMETER(kdata=600)

C     INTEGER n1, n2, n3, n4, n5
      REAL x1(kdata)
      REAL y1(kdata)

* local

      REAL yg(kw)
      REAL qy
      REAL a0, a1, a2, a3, a4, a5, a6, a7
      REAL b0, b1, b2, b3, b4
      REAL xs
      REAL t
      INTEGER i, iw, n, idum
      REAL lambda
      REAL sumA, sumB, chi

      INTEGER mabs, myld, xsvers(3), qyvers(3)
      CHARACTER(llog) :: logmsg(nlog)
      CHARACTER(lspc) :: spc


***************** H2O2 photolysis *************************

* Setting photolysis index and cross section/quantum yield options
      j = j + 1
      jlabel(j) = 'H2O2 -> 2 OH'

* cross section options
      spc = "H2O2"
      xsvers = (/2, 3, 3/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "from old TUV data"
      logmsg(2) = "from JPL94"
      logmsg(3) = "from JPL 2011 evaluation"
      logmsg(4) = "from IUPAC"
      CALL set_option(4,spc,"abs",logmsg,xsvers,mabs)

* quantum yield options
      qyvers = (/1, 2, 2/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "set to unity"
      logmsg(2) = "IUPAC recommended"
      CALL set_option(2,spc,"yld",logmsg,qyvers,myld)


* cross sections from Lin et al. 1978

      IF(mabs==1) THEN
*--------------- TUV old version (start) -------------------------------
      OPEN(UNIT=kin,FILE='DATAJ1/ABS/H2O2_lin.abs',STATUS='old')
      DO i = 1, 7
         READ(kin,*)
      ENDDO
      n = 32
      DO i = 1, n
         READ(kin,*) x1(i), y1(i)
         y1(i) = y1(i) * 1.E-20
      ENDDO
      CLOSE (kin)

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)

*--------------- TUV old version (end) ---------------------------------

      ELSEIF(mabs==2) THEN

*--------------- TUV standard version (start) --------------------------
* cross section from JPL94 (identical to JPL97)
* tabulated data up to 260 nm

      OPEN(UNIT=kin,FILE='DATAJ1/ABS/H2O2_jpl94.abs',STATUS='old')
      READ(kin,*) idum,n
      DO i = 1, idum-2
         READ(kin,*)
      ENDDO
      DO i = 1, n
         READ(kin,*) x1(i), y1(i)
         y1(i) = y1(i) * 1.E-20
      ENDDO
      CLOSE (kin)


      OPEN(UNIT=kin,FILE='DATAJ1/ABS/H2O2_Kahan.abs',STATUS='old')
      DO i = 1, 494
         n = n + 1
         READ(kin,*) x1(n), y1(n)
      ENDDO
      CLOSE (kin)

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)

* JPL XS parameterisation values
      A0 = 6.4761E+04
      A1 = -9.2170972E+02
      A2 = 4.535649
      A3 = -4.4589016E-03
      A4 = -4.035101E-05
      A5 = 1.6878206E-07
      A6 = -2.652014E-10
      A7 = 1.5534675E-13

      B0 = 6.8123E+03
      B1 = -5.1351E+01
      B2 = 1.1522E-01
      B3 = -3.0493E-05
      B4 = -1.0924E-07

      DO iw = 1, nw - 1

* Parameterization (JPL94)
* Range 260-350 nm; 200-400 K

         IF ((wl(iw) .GE. 260.) .AND. (wl(iw) .LT. 350.)) THEN

           lambda = wc(iw)
           sumA = ((((((A7*lambda + A6)*lambda + A5)*lambda +
     >                  A4)*lambda +A3)*lambda + A2)*lambda +
     >                  A1)*lambda + A0
           sumB = (((B4*lambda + B3)*lambda + B2)*lambda +
     >               B1)*lambda + B0

           DO i = 1, nz
              t = MIN(MAX(tlev(i),200.),400.)
              chi = 1./(1.+EXP(-1265./t))
              xs = (chi * sumA + (1.-chi)*sumB)*1E-21
              absq(i,iw) = xs
           ENDDO
         ELSE
           DO i = 1, nz
              absq(i,iw) = yg(iw)
           ENDDO
         ENDIF

      ENDDO



*--------------- TUV standard version (end) ----------------------------

      ELSEIF(mabs==3) THEN

*--------------- MCM standard version (start) --------------------------
* cross section from IUPAC recommendation, last update 2.10.2001
* tabulated data between 195 and 350nm

      OPEN(UNIT=kin,FILE='DATAJ1/ABS/H2O2_IUPAC.abs',STATUS='old')
      DO i = 1, 5
         READ(kin,*)
      ENDDO

      n = 32
      DO i = 1, n
         READ(kin,*) x1(i), y1(i)
         y1(i) = y1(i) * 1.E-20
      ENDDO
      CLOSE (kin)

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)

*--------------- MCM standard version (end) ----------------------------
      ENDIF

*--- Cross section for cases 1 and 3:
      IF(mabs==1 .OR. mabs==3) THEN
      DO iw = 1, nw - 1
        DO i = 1, nz
          absq(i,iw) = yg(iw)
        ENDDO
      ENDDO
      ENDIF


* quantum yields

      DO iw = 1, nw - 1
        IF(myld == 1) THEN
          qy = 1
         ELSEIF(myld == 2) THEN
          IF(wl(iw) .GE. 220.) THEN
            qy = 1.
           ELSE
            qy = 0.85
          ENDIF
        ENDIF
        DO i = 1, nz
          sq(j,i,iw) = absq(i,iw)*qy
        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r09(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! CHBr3

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product of (cross section) x (quantum yield) for CHBr3 photolysis=*
*=          CHBr3 + hv -> Products                                           =*
*=  Cross section: Choice of data from Atlas (?Talukdar???) or JPL97         =*
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

      INTEGER n1, n2, n3, n4, n5
      REAL x1(kdata), x2(kdata), x3(kdata), x4(kdata), x5(kdata)
      REAL y1(kdata), y2(kdata), y3(kdata), y4(kdata), y5(kdata)

* local

      REAL yg(kw), yg1(kw), yg2(kw), yg3(kw), yg4(kw), yg5(kw)

      real t
C      real qy

      INTEGER i, iw
      INTEGER iz

      INTEGER kopt, opt(3)
      CHARACTER(llog) :: logmsg(nlog)
      CHARACTER(lspc) :: spc

*_______________________________________________________________________

* calculate centre wavelengths
      DO 5, iw = 1, nw - 1
         wc(iw) = (wl(iw) + wl(iw+1))/2.
 5    CONTINUE


***************** CHBr3 photolysis *************************

* Setting photolysis index and cross section/quantum yield options
      j = j + 1
      jlabel(j) = 'CHBr3 -> Products'

* cross section options
      spc = "H2O2"
      opt = (/2, 2, 1/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "from Elliot Atlas, 1997"
      logmsg(2) = "JPL 1997"
      CALL set_option(2,spc,"abs",logmsg,opt,kopt)


* cross sections

      if (kopt .eq. 1) then

      OPEN(UNIT=kin,FILE='DATAJ1/ABS/CHBr3.abs',STATUS='old')
      DO i = 1, 5
         READ(kin,*)
      ENDDO

      n5 = 25
      n4 = 27
      n3 = 29
      n2 = 31
      n1 = 39
      DO i = 1, n5
         READ(kin,*) x1(i), y1(i), y2(i), y3(i), y4(i), y5(i)
      ENDDO
      do i = n5 + 1, n4
         READ(kin,*) x1(i), y1(i), y2(i), y3(i), y4(i)
      enddo
      do i = n4 + 1, n3
         READ(kin,*) x1(i), y1(i), y2(i), y3(i)
      enddo
      do i = n3 + 1, n2
         READ(kin,*) x1(i), y1(i), y2(i)
      enddo
      do i = n2 + 1, n1
         READ(kin,*) x1(i), y1(i)
      enddo
      CLOSE (kin)

      do i = 1, n1
         y1(i) = y1(i) * 1.e-23
      enddo
      do i = 1, n2
         x2(i) = x1(i)
         y2(i) = y2(i) * 1.e-23
      enddo
      do i = 1, n3
         x3(i) = x1(i)
         y3(i) = y3(i) * 1.e-23
      enddo
      do i = 1, n4
         x4(i) = x1(i)
         y4(i) = y4(i) * 1.e-23
      enddo
      do i = 1, n5
         x5(i) = x1(i)
         y5(i) = y5(i) * 1.e-23
      enddo

      CALL interpol(x1,y1,kdata,n1,nw,wl,jlabel(j),y1(1),0.,yg1)
      CALL interpol(x2,y2,kdata,n2,nw,wl,jlabel(j),y2(1),0.,yg2)
      CALL interpol(x3,y3,kdata,n3,nw,wl,jlabel(j),y3(1),0.,yg3)
      CALL interpol(x4,y4,kdata,n4,nw,wl,jlabel(j),y4(1),0.,yg4)
      CALL interpol(x5,y5,kdata,n5,nw,wl,jlabel(j),y5(1),0.,yg5)


* quantum yield = 1

C      qy = 1.
      DO iw = 1, nw - 1
         DO iz = 1, nz

            t = tlev(iz)

            if (t .ge. 296.) then
               yg(iw) = yg1(iw)

            else if(t .ge. 286.) then
               yg(iw) = yg1(iw) + (t-286.)*(yg2(iw)-yg1(iw))/10.

            else if(t .ge. 276.) then
               yg(iw) = yg2(iw) + (t-276.)*(yg3(iw)-yg2(iw))/10.

            else if(t .ge. 266.) then
               yg(iw) = yg3(iw) + (t-266.)*(yg4(iw)-yg3(iw))/10.

            else if(t .ge. 256.) then
               yg(iw) = yg4(iw) + (t-256.)*(yg5(iw)-yg4(iw))/10.

            else if(t .lt. 256.) then
               yg(iw) = yg5(iw)

            endif

            sq(j,iz,iw) = yg(iw)!*qy

         ENDDO
      ENDDO

* jpl97, with temperature dependence formula,
*w = 290 nm to 340 nm,
*T = 210K to 300 K
*sigma, cm2 = exp((0.06183-0.000241*w)*(273.-T)-(2.376+0.14757*w))

      ELSEIF (kopt .EQ. 2) THEN

      OPEN(UNIT=kin,FILE='DATAJ1/ABS/CHBr3.jpl97',STATUS='old')
      DO i = 1, 6
         READ(kin,*)
      ENDDO
      n1 = 87
      DO i = 1, n1
         READ(kin,*) x1(i), y1(i)
         y1(i) = y1(i) * 1.e-20
      ENDDO
      CLOSE(kin)

      CALL interpol(x1,y1,kdata,n1,nw,wl,jlabel(j),y1(1),0.,yg1)


* quantum yield = 1

C      qy = 1.
      DO iw = 1, nw - 1
         DO iz = 1, nz

            t = tlev(iz)
            yg(iw) = yg1(iw)

            IF (wc(iw) .GT. 290. .AND. wc(iw) .LT. 340.
     $           .AND. t .GT. 210 .AND. t .LT. 300) THEN
               yg(iw) = EXP((0.06183-0.000241*wc(iw))*(273.-T)-
     $              (2.376+0.14757*wc(iw)))
            ENDIF

            sq(j,iz,iw) = yg(iw)!*qy
         ENDDO
      ENDDO

      ENDIF

      END

*=============================================================================*

      SUBROUTINE pxCH2O(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! HCHO (new SR)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  JPL 2011 or IUPAC 2013 recommendation.                                   =*
*=  Provide product of (cross section) x (quantum yield) for CH2O photolysis =*
*=        (a) CH2O + hv -> H + HCO                                           =*
*=        (b) CH2O + hv -> H2 + CO                                           =*
*=  written by s. madronich march 2013                                       =*
*=  extended by P. Bräuer, Nov. 2015                                         =*
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

      INTEGER kdata
      PARAMETER(kdata=200)

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

      INTEGER i, j, iz, iw

* data arrays

      INTEGER n, n1, n2
      REAL x1(kdata), x2(kdata)
      REAL y298(kdata), tcoef(kdata)
      REAL qr(kdata), qm(kdata), x

* local

      REAL yg1(kw), yg2(kw), yg3(kw), yg4(kw)
      REAL ak300, akt, sig
      real qyr300, qym300, qymt, fact
      REAL dum

      INTEGER mabs, myld, xsvers(3), qyvers(3)
      CHARACTER(llog) :: logmsg(nlog)
      CHARACTER(lspc) :: spc


***************** CH2O photolysis *************************

* working grid arrays:
* yg  = T-indep. high res. IUPAC cross sections
* yg1 = cross section at 298K
* yg2 = temperature correction coefficient for cross section
* yg3 = quantum yields for radical channel, H + HCO
* yg4 = quantum yields for molecular channel, H2 + CO.

* Setting photolysis index and cross section/quantum yield options
      j = j+1
      jlabel(j) = 'CH2O -> H + HCO'
      j = j+1
      jlabel(j) = 'CH2O -> H2 + CO'

* cross section options
      spc = "CH2O"
      xsvers = (/1, 4, 1/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "from JPL 2011"
      logmsg(2) = "from IUPAC high res. data"
      logmsg(3) = "from IUPAC T-dep. low res. data"
      logmsg(4) = "from IUPAC mixed high res. data with T-dep"
      CALL set_option(4,spc,"abs",logmsg,xsvers,mabs)

* quantum yield options
      qyvers = (/1, 2, 1/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "from JPL 2011"
      logmsg(2) = "from IUPAC 2013"
      CALL set_option(2,spc,"yld",logmsg,qyvers,myld)




*_______________________________________________________________________

* calculate centre wavelengths
      DO 5, iw = 1, nw - 1
         wc(iw) = (wl(iw) + wl(iw+1))/2.
 5    CONTINUE


* cross sections

      IF(mabs.EQ.1) THEN
* read JPL2011 cross section data:

        OPEN(UNIT=kin,FILE='DATAJ1/CH2O/CH2O_jpl11.abs'
     $      ,STATUS='old')
        do i = 1, 4
           read(kin,*)
        enddo
        n = 150
        n1 = n
        n2 = n
        DO i = 1, n
           READ(kin,*) x1(i), y298(i), tcoef(i)
           x2(i) = x1(i)
           y298(i) = y298(i) * 1.e-20
           tcoef(i) = tcoef(i) * 1.e-24
        ENDDO
        CLOSE(kin)

*     terminate endpoints and interpolate to working grid

        CALL interpol(x1,y298,kdata,n1,nw,wl,jlabel(j),0.,0.,yg1)
        CALL interpol(x2,tcoef,kdata,n2,nw,wl,jlabel(j),0.,0.,yg2)

      ENDIF

      IF(mabs .EQ. 2 .OR. mabs .EQ. 4) THEN

* read IUPAC 2013 high res. data

         OPEN(UNIT=kin,FILE='DATAJ1/CH2O/HCHO_IUPAChr.abs'
     $        ,STATUS='old')
         DO i = 1, 6
            READ(kin,*)
         ENDDO
         n = 150
         DO i = 1, n
            READ(kin,*) x1(i), y298(i)
         ENDDO
         CALL interpol(x1,y298,kdata,n,nw,wl,jlabel(j),0.,0.,yg1)
      ENDIF

      IF(mabs .EQ. 3 .OR. mabs .EQ. 4) THEN

* read IUPAC 2013 T-dep. low res. data

            OPEN(UNIT=kin,FILE='DATAJ1/CH2O/HCHO_IUPAClr.abs',
     &           STATUS='old')
            DO i = 1, 8
               READ(kin,*)
            ENDDO
            n = 37
            n1 = n
            n2 = n
            DO i = 1, n
               READ(kin,*) x1(i), y298(i), tcoef(i)
               y298(i)  = y298(i) * 1.e-20
               tcoef(i) = tcoef(i) * 1.e-24
               x2(i)    = x1(i)
            ENDDO
            CLOSE(kin)

*     terminate endpoints and interpolate to working grid
        IF(mabs==3) THEN
          CALL interpol(x1,y298,kdata,n1,nw,wl,jlabel(j),0.,0.,yg1)
        ENDIF

        CALL interpol(x2,tcoef,kdata,n2,nw,wl,jlabel(j),0.,0.,yg2)

      ENDIF


* quantum yields: Read, terminate, interpolate:

      IF(myld.EQ.1) THEN
* read JPL2011 quantum yield data:

        OPEN(UNIT=kin,FILE='DATAJ1/CH2O/CH2O_jpl11.yld',STATUS='old')
           DO i = 1, 4
              READ(kin,*)
           ENDDO
           n = 112
           n1 = n
           n2 = n
           DO i = 1, n
              READ(kin,*) x1(i), qr(i), qm(i)
              x2(i) = x1(i)
           ENDDO
           CLOSE(kin)

           CALL interpol(x1,qr,kdata,n1,nw,wl,jlabel(j),qr(1),0.,yg3)
           CALL interpol(x2,qm,kdata,n2,nw,wl,jlabel(j),qm(1),0.,yg4)

      ELSE IF(myld .EQ. 2) THEN
* read IUPAC 2013 quantum yield data:

         OPEN(UNIT=kin,FILE='DATAJ1/CH2O/HCHO_IUPAC.yld',STATUS='old')
         DO i = 1, 9
            READ(kin,*)
         ENDDO
         n = 86
         n1 = n
         n2 = n
         DO i = 1, n
            READ(kin,*) x1(i), qr(i), qm(i), dum
            x2(i) = x1(i)
         ENDDO
         CLOSE(kin)

         CALL interpol(x1,qr,kdata,n1,nw,wl,jlabel(j),qr(1),0.,yg3)
         CALL interpol(x2,qm,kdata,n2,nw,wl,jlabel(j),qm(1),0.,yg4)
      ENDIF

* combine gridded quantities:
* yg1 = cross section at 298K
* yg2 = temperature correction coefficient for cross section
* yg3 = quantum yields for radical channel, H + HCO
* yg4 = quantum yields for molecular channel, H2 + CO.

      DO iz = 1, nz

         DO iw = 1, nw - 1

* correct cross section for temperature dependence:
           sig = yg1(iw)
           IF(mabs==1 .OR. mabs>=3 .AND.
     &       (wl(iw) .GE. 251.7 .AND. wl(iw) .LE. 350.)) THEN
               sig = yg1(iw) + yg2(iw) * (tlev(iz) - 298.)
           ENDIF

* assign room temperature quantum yields for radical and molecular channels

           qyr300 = yg3(iw)
           qym300 = yg4(iw)
           qymt = qym300

* between 330 ande 360 nm, molecular channel is pressure and temperature dependent.

* JPL quantum yield data
           IF (wc(iw) .ge. 330. .and. wc(iw) .lt. 360. .and.
     $         qym300 .gt. 0.) then

             IF(myld==1) THEN
               ak300 = 1./qym300  - 1./(1. - qyr300)
               ak300 = ak300/2.45e19
               akt   = ak300 * (1. + 0.05 * (wc(iw) - 329.)
     $                 * (300. - tlev(iz))/80.)

               qymt = 1./(1./(1.-qyr300) + akt*airden(iz))

              ELSE
* IUPAC quantum yield data
               CALL intermat(wc(iw),tlev(iz),329.,353.,220.,300.,
     &                       1.12,0.26,2.47,0.39,fact)
               fact = fact * 1e-19
               x = (1.-qyr300)*fact*2.465E19
               qymt=MAX(0.,MIN((1.-qyr300),
     &              qym300*(1.+x)/(1.+x*airden(iz)/2.465E19)))
           ENDIF
         ENDIF

         sq(j-1,iz,iw) = sig * qyr300
         sq(j  ,iz,iw) = sig * qymt

         ENDDO
      ENDDO

      RETURN
      END

*=============================================================================*

      SUBROUTINE r10(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! obsolete HCHO

*-----------------------------------------------------------------------------*
*=  C U R R E N T L Y   O B S O L E T E ! ! !                                =*
*=  SUBROUTINE pxCH2O used instead.                                          =*
*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product of (cross section) x (quantum yield) for CH2O photolysis =*
*=        (a) CH2O + hv -> 2 HO2 + CO (- 2 O2)                               =*
*=        (b) CH2O + hv -> H2 + CO                                           =*
*=  Cross section: Choice between                                            =*
*=                 1) Bass et al., 1980 (resolution: 0.025 nm)               =*
*=                 2) Moortgat and Schneider (resolution: 1 nm)              =*
*=                 3) Cantrell et al. (orig res.) for > 301 nm,              =*
*=                    IUPAC 92, 97 elsewhere                                 =*
*=                 4) Cantrell et al. (2.5 nm res.) for > 301 nm,            =*
*=                    IUPAC 92, 97 elsewhere                                 =*
*=                 5) Rogers et al., 1990                                    =*
*=                 6) new NCAR recommendation, based on averages of          =*
*=                    Cantrell et al., Moortgat and Schneider, and Rogers    =*
*=                    et al.                                                 =*
*=  Quantum yield: Choice between                                            =*
*=                 1) Evaluation by Madronich 1991 (unpublished)             =*
*=                 2) IUPAC 89, 92, 97                                       =*
*=                 3) Madronich, based on 1), updated 1998.                  =*
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

      INTEGER kdata
      PARAMETER(kdata=16000)

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

      INTEGER j, iw

* data arrays

      INTEGER n
      real x(kdata), y(kdata)
      INTEGER n1, n2, n3
      REAL x1(kdata), x2(kdata), x3(kdata)
      REAL y1(kdata), y2(kdata), y3(kdata)

* local

      REAL yg1(kw), yg2(kw), yg3(kw), yg4(kw), yg5(kw)
      REAL phi1, phi2, phi20, ak300, akt
      REAL qy1, qy2

      REAL sig, slope
      REAL t
      REAL dum
      INTEGER idum

      INTEGER i
      INTEGER irev

      INTEGER mabs, myld, xsvers(3), qyvers(3)
      CHARACTER(llog) :: logmsg(nlog)
      CHARACTER(lspc) :: spc


***************** CH2O photolysis *************************

* working grid arrays:
*     yg1 = cross section at a specific temperature
*     yg2, yg3 = cross sections at different temp or slope, for calculating
*                temperature depedence
*     yg4 = quantum yield data for radical channel
*     yg5 = quantum yield data for molecular channel

* Setting photolysis index and cross section/quantum yield options
      j = j+1
      jlabel(j) = 'CH2O -> H + HCO'
      j = j+1
      jlabel(j) = 'CH2O -> H2 + CO'

* cross section options
      spc = "CH2O"
      xsvers = (/6, 6, 2/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "from Bass et al."
      logmsg(2) = "from Moortgat and Schneider (IUPAC)"
      logmsg(3) =
     &  "in high res. from Cantrell et al. 1990 completed by IUPAC data"
      WRITE(logmsg(4), "(2A)") "in low res. with T-dep. ",
     & "from Cantrell et al. 1990 completed by IUPAC data"
      logmsg(5) = "from Rogers et al. 1990"
      logmsg(6) = "NCAR recommendation"
      CALL set_option(6,spc,"abs",logmsg,xsvers,mabs)

* quantum yield options
      qyvers = (/1, 2, 1/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "from Madronich, 1991 (unpublished)"
      logmsg(2) = "from IUPAC recommendation"
      logmsg(3) = "from JPL recommendation"
      CALL set_option(3,spc,"yld",logmsg,qyvers,myld)


*_______________________________________________________________________

* calculate centre wavelengths
      DO 5, iw = 1, nw - 1
         wc(iw) = (wl(iw) + wl(iw+1))/2.
 5    CONTINUE

* Input data options:
* mabs for absorption:
* 1:  DATAJ1/CH2O/CH2O_nbs.abs'
*     from Bass et al., Planet. Space. Sci. 28, 675, 1980.
*     over 258.750-359.525 in 0.025 nm steps
* 2:  DATAJ1/CH2O_iupac1.abs
*     Moortgat and Schneider, personal communication as reported in IUPAC 89, 92, 97
*     at 285K.  Over 240-360 nm in 1 nm bins (note that IUPAC 89,92,97 incorectly
*     claims 0.5 nm intervals in footnote)
* 3:  DATAJ1/CH2O/ch2o_can_hr.abs for wc > 301 nm, temperature dependent
*     DATAJ1/CH2O/ch2o_iupac1.abs elsewhere
*     from Cantrell et al. 1990 for wc > 301 nm.  Original data from Cantrell,
*     at high resolution
* 4:  DATAJ1/CH2O/CH2O_can_lr.abs for wc > 301 nm, temperature dependent
*     DATAJ1/CH2O/CH2O_iupac1.abs elsewhere
*     from Cantrell et al. 1990 for wc > 301 nm.  Data from Cantrell et al., as
*     reported by IUPAC'92,'97.  On 2.5 nm intervals.
* 5:  DATAJ1/CH2O/CH2O_rog.abs'
*     from Rogers et al., J. Phys. Chem. 94, 4011, 1990.
* 6:  DATAJ2/CH2O_ncar.abs
*     new NCAR recommendation, based on averages of Moortgat and Schneider, Cantrell et al.,
*     and Rogers.
* mopt2 for quantum yields:
* 1:  DATAJ1/CH2O/CH2O_i_mad.yld and
*     DATAJ1/CH2O/CH2O_ii_mad.yld
*     evaluated by Madronich, 1991, unpublished
* 2:  DATAJ1/CH2O/CH2O_iupac.yld
*     from IUPAC'89, '92, '97
* 3:  DATAJ1/CH2O/CH2O_jpl97.abs'
*     based on Madronich 1991 unpublished evaluation, updated Jan 1998.


      IF (mabs .EQ. 1) THEN

* read NBS/Bass data

         OPEN(UNIT=kin,FILE='DATAJ1/CH2O/CH2O_nbs.abs'
     $        ,STATUS='old')
         n = 4032
         DO i = 1, n
            READ(kin,*) x(i), y(i)
         ENDDO
         CALL interpol(x,y,kdata,n,nw,wl,jlabel(j),0.,0.,yg1)

      ELSEIF (mabs .EQ. 2 .OR. mabs .EQ. 3 .OR. mabs .EQ. 4) THEN

         OPEN(UNIT=kin,FILE='DATAJ1/CH2O_iupac1.abs',STATUS='old')
         DO i = 1, 4
            READ(kin,*)
         ENDDO
         n = 121
         DO i = 1, n
            READ(kin,*) x(i), y(i)
            y(i) = y(i) * 1.e-20
         ENDDO
         CLOSE(kin)
         CALL interpol(x,y,kdata,n,nw,wl,jlabel(j-1),0.,0.,yg1)

         IF(mabs .EQ. 3) THEN

* data are on wavenumber grid (cm-1), so convert to wavelength in nm:
* grid was on increasing wavenumbers, so need to reverse to get increasing
* wavelengths
* cross section assumed to be zero for wavelengths longer than 360 nm
* if y1 < 0, then make = 0 (some negative cross sections, actually 273 K intercepts
* are in the original data,  Here, make equal to zero)

         OPEN(kin,FILE='DATAJ1/CH2O/CH2O_can_hr.abs',STATUS='old')
         READ(kin,*) idum, n
         DO i = 1, idum-2
            READ(kin,*)
         ENDDO
         DO i = 1, n
            READ(kin,*) x1(i), y1(i), y2(i)
            x1(i) = 1./x1(i) * 1E7
            IF (x1(i) .GT. 360.) THEN
               y1(i) = 0.
               y2(i) = 0.
            ENDIF
         ENDDO
         CLOSE(kin)

         DO i = 1, n/2
            irev = n+1-i
            dum = x1(i)
            x1(i) = x1(irev)
            x1(irev) = dum
            dum = y1(i)
            y1(i) = y1(irev)
            y1(irev) = dum
            dum = y2(i)
            y2(i) = y2(irev)
            y2(irev) = dum
         ENDDO
         DO i = 1, n
            x2(i) = x1(i)
            y1(i) = max(y1(i),0.)
         ENDDO
         n1 = n
         n2 = n

         CALL interpol(x1,y1,kdata,n1,nw,wl,jlabel(j),0.,0.,yg2)
         CALL interpol(x2,y2,kdata,n2,nw,wl,jlabel(j),0.,0.,yg3)

      ELSEIF(mabs .eq. 4) THEN

            OPEN(UNIT=kin,FILE='DATAJ1/CH2O/CH2O_can_lr.abs',
     $        STATUS='old')
            DO i = 1, 4
               READ(kin,*)
            ENDDO
            n = 23
            DO i = 1, n
               READ(kin,*) x2(i), y2(i), y3(i)
               x3(i) = x2(i)
            ENDDO
            CLOSE(kin)
            n2 = n
            n3 = n

            CALL interpol(x2,y2,kdata,n2,nw,wl,jlabel(j),0.,0.,yg2)
            CALL interpol(x3,y3,kdata,n3,nw,wl,jlabel(j),0.,0.,yg3)
         ENDIF

      ELSEIF (mabs .EQ. 5) THEN

* read Rodgers data

         OPEN(UNIT=kin,FILE='DATAJ1/CH2O/CH2O_rog.abs'
     $        ,STATUS='old')
         DO i = 1, 10
            READ(kin,*)
         ENDDO
         n = 261
         DO i = 1, n
            READ(kin,*) x(i), y(i)
            y(i) = y(i) * 1.e-20
         ENDDO
         CALL interpol(x,y,kdata,n,nw,wl,jlabel(j-1),0.,0.,yg1)

      ELSEIF(mabs .EQ. 6) THEN

            OPEN(UNIT=kin,FILE='DATAJ1/CH2O/CH2O_ncar.abs',STATUS='old')
            DO i = 1, 3
               READ(kin,*)
            ENDDO
            n = 126
            DO i = 1, n
               READ(kin,*) x2(i), y2(i), y3(i)
               x3(i) = x2(i)
            ENDDO
            CLOSE(kin)
            n2 = n
            n3 = n

            CALL interpol(x2,y2,kdata,n2,nw,wl,jlabel(j),0.,0.,yg2)
            CALL interpol(x3,y3,kdata,n3,nw,wl,jlabel(j),0.,0.,yg3)

      ENDIF

* quantum yield

      IF (myld .EQ. 1) THEN

         OPEN(UNIT=kin,FILE='DATAJ1/CH2O/CH2O_i_mad.yld',STATUS='old')
         DO i = 1, 11
            READ(kin,*)
         ENDDO
         n = 20
         DO i = 1, n
            READ(kin,*) x(i), y(i)
         ENDDO
         CLOSE(kin)
         CALL interpol(x,y,kdata,n,nw,wl,jlabel(j),y(1),0.,yg4)

         OPEN(UNIT=kin,FILE='DATAJ1/CH2O/CH2O_ii_mad.yld',STATUS='old')
         DO i = 1, 9
            READ(kin,*)
         ENDDO
         n = 33
         DO i = 1, n
            READ(kin,*) x(i), y(i)
         ENDDO
         CLOSE(kin)
         CALL interpol(x,y,kdata,n,nw,wl,jlabel(j),y(1),0.,yg5)

      ELSEIF(myld .EQ. 2) then

         OPEN(UNIT=kin,FILE='DATAJ1/CH2O/CH2O_iupac.yld',STATUS='old')
         DO i = 1, 7
            READ(kin,*)
         ENDDO
         n = 13
         DO i = 1, n
            READ(kin,*) x1(i), y1(i), y2(i)
            x2(i) = x1(i)
         ENDDO
         CLOSE(kin)
         n1 = n
         n2 = n

         CALL interpol(x1,y1,kdata,n1,nw,wl,jlabel(j),y1(1),0.,yg4)
         CALL interpol(x2,y2,kdata,n2,nw,wl,jlabel(j),y2(1),0.,yg5)

* box-filling interpolation.
c         DO i = 1, n
c            READ(kin,*) x1(i), y1(i), y2(i)
c            x1(i) = x1(i) - 5.0
c            x2(i) = x1(i)
c         ENDDO
c         n = n + 1
c         x1(n) = x1(n-1) + 5.0
c         x2(n) = x1(n)
c         CLOSE(kin)
c         DO i = 1, n-1
c            y1(i) = y1(i) * (x1(i+1)-x1(i))
c         ENDDO
c         CALL inter3(nw,wl,yg4,n,x1,y1,0)
c         DO iw = 1, nw-1
c            yg4(iw) = yg4(iw)/(wl(iw+1)-wl(iw))
c         ENDDO
c         DO i = 1, n-1
c            y2(i) = y2(i) * (x2(i+1)-x2(i))
c         ENDDO
c         CALL inter3(nw,wl,yg5,n,x2,y2,0)
c         DO iw = 1, nw-1
c            yg5(iw) = yg5(iw)/(wl(iw+1)-wl(iw))
c         ENDDO

      ELSEIF(myld .EQ. 3) THEN

         OPEN(UNIT=kin,FILE='DATAJ1/CH2O/CH2O_jpl97.abs',STATUS='old')
         DO i = 1, 4
            READ(kin,*)
         ENDDO
         n = 23
         DO i = 1, n
            READ(kin,*) x1(i), dum, dum, dum, dum, y1(i), y2(i)
            x2(i) = x1(i)
         ENDDO
         CLOSE(kin)
         n1 = n
         n2 = n

         CALL interpol(x1,y1,kdata,n1,nw,wl,jlabel(j),y1(1),0.,yg4)
         CALL interpol(x2,y2,kdata,n2,nw,wl,jlabel(j),y2(1),0.,yg5)

      ENDIF

* combine
* y1 = xsect
* y2 = xsect(223), Cantrell et al.
* y3 = xsect(293), Cantrell et al.
* y4 = qy for radical channel
* y5 = qy for molecular channel
* pressure and temperature dependent for w > 330.

      DO iw = 1, nw - 1

         IF (mabs .eq. 6) THEN
            sig = yg2(iw)
         ELSE
            sig = yg1(iw)
         ENDIF

         DO i = 1, nz

* correct cross section for temperature dependence for > 301. nm

            IF (wl(iw) .GE. 301.) THEN
               t = MAX(223.15, MIN(tlev(i), 293.15))
               IF (mabs .EQ. 3 .OR. mabs .EQ. 6) THEN
                  sig = yg2(iw) + yg3(iw) * (t - 273.15)

               ELSEIF (mabs .EQ. 4) THEN
                  slope = (yg3(iw) - yg2(iw)) / (293. - 223.)
                  sig = yg2(iw) + slope * (t - 223.)

               ENDIF

            ENDIF
            sig = MAX(sig, 0.)

* quantum yields:
* temperature and pressure dependence beyond 330 nm

            qy1 = yg4(iw)
            IF ( (wc(iw) .GE. 330.) .AND. (yg5(iw) .GT. 0.) ) THEN
               phi1 = yg4(iw)
               phi2 = yg5(iw)
               phi20 = 1. - phi1
               ak300=((1./phi2)-(1./phi20))/2.54E+19
               akt=ak300*(1.+61.69*(1.-tlev(i)/300.)*(wc(iw)/329.-1.))
               qy2 = 1. / ( (1./phi20) + airden(i)*akt)

            ELSE
               qy2 = yg5(iw)
            ENDIF
            qy2 = MAX(0.,qy2)
            qy2 = MIN(1.,qy2)

            sq(j-1,i,iw) = sig * qy1
            sq(j  ,i,iw) = sig * qy2

         ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r11(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! CH3CHO

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for CH3CHO photolysis: =*
*=      (a)  CH3CHO + hv -> CH3 + HCO                                        =*
*=      (b)  CH3CHO + hv -> CH4 + CO                                         =*
*=      (c)  CH3CHO + hv -> CH3CO + H                                        =*
*=  Cross section:  Choice between                                           =*
*=                   (1) IUPAC 97 data, from Martinez et al.                 =*
*=                   (2) Calvert and Pitts                                   =*
*=                   (3) Martinez et al., Table 1 scanned from paper         =*
*=                   (4) KFA tabulations                                     =*
*=                   (5) JPL 2011                                            =*
*=  Quantum yields: Choice between                                           =*
*=                   (1) IUPAC 97, pressure correction using Horowith and    =*
*=                                 Calvert, 1982                             =*
*=                   (2) NCAR data file, from Moortgat, 1986                 =*
*=                   (3) IUPAC 2013 including pressure correction            =*
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

      INTEGER i, n, n1, n2, n3
      REAL x1(kdata), x2(kdata), x3(kdata)
      REAL y1(kdata), y2(kdata), y3(kdata)

* local

      REAL yg(kw), yg1(kw), yg2(kw), yg3(kw), yg4(kw), yg5(kw)
      REAL qy1, qy2, qy3, qytot
      REAL qy1_n0, qy1_0, x
      REAL sig
      REAL quench
      INTEGER  iw

      INTEGER mabs, myld, mpress, xsvers(3), qyvers(3), pvers(3)
      CHARACTER(llog) :: logmsg(nlog)
      CHARACTER(lspc) :: spc


***************** HNO3 photolysis *************************

* Setting photolysis index and cross section/quantum yield/pressure dep. options
      j = j+1
      jlabel(j) = 'CH3CHO -> CH3 + HCO'
      j = j+1
      jlabel(j) = 'CH3CHO -> CH4 + CO'
      j = j+1
      jlabel(j) = 'CH3CHO -> CH3CO + H'
      j = j+1
      jlabel(j) = 'genCH3CHO(poly)'

* cross section options
      spc = "CH3CHO"
      xsvers = (/5, 1, 1/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "from IUPAC"
      logmsg(2) = "from Calvert and Pitts"
      logmsg(3) = "from Martinez et al"
      logmsg(4) = "from KFA tabulations"
      logmsg(5) = "from JPL 2011"
      CALL set_option(5,spc,"abs",logmsg,xsvers,mabs)

* quantum yield options
      qyvers = (/1, 3, 1/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "from IUPAC 97"
      logmsg(2) = "from NCAR"
      logmsg(3) = "from IUPAC 2013"
      CALL set_option(3,spc,"yld",logmsg,qyvers,myld)

* pressure-dependence options
      pvers = (/1, 2, 1/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "estimated from IUPAC 2006"
      logmsg(2) = "from IUPAC 2013"
      CALL set_option(2,spc,"que",logmsg,pvers,mpress)


* cross sections

      IF (mabs .EQ. 1) THEN

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

      ELSEIF(mabs .EQ. 2) THEN

* cross section from Calvert and  Pitts

         OPEN(UNIT=kin,FILE='DATAJ1/CH3CHO/d021_cp.abs',STATUS='old')
         DO i = 1, 14
            READ(kin,*)
         ENDDO
         n = 54
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
            x1(i) = x1(i)/10.
            y1(i) = y1(i) * 3.82E-21
         ENDDO
         CLOSE (kin)

      ELSEIF(mabs .EQ. 3) THEN

         OPEN(UNIT=kin,FILE='DATAJ1/CH3CHO/CH3CHO_mar.abs',STATUS='old')
         DO i = 1, 3
            READ(kin,*)
         ENDDO
         n = 106
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
         ENDDO
         CLOSE (kin)

      ELSEIF(mabs .EQ. 4) THEN

* cross section from KFA tables
* ch3cho.001 - Calvert and Pitts 1966
* ch3cho.002 - Meyrahn thesis 1984
* ch3cho.003 - Schneider and Moortgat, priv comm. MPI Mainz 1989, 0.012 nm resol.
* ch3cho.004 - Schneider and Moortgat, priv comm. MPI Mainz 1989, 0.08  nm resol.
* ch3cho.005 - IUPAC'92
* ch3cho.006 - Libuda, thesis Wuppertal 1992

c         OPEN(UNIT=kin,FILE='DATAJ2/KFA/ch3cho.001',STATUS='old')
C         n = 217
c         OPEN(UNIT=kin,FILE='DATAJ2/KFA/ch3cho.002',STATUS='old')
c         n = 63
c         OPEN(UNIT=kin,FILE='DATAJ2/KFA/ch3cho.003',STATUS='old')
c         n = 13738
c         OPEN(UNIT=kin,FILE='DATAJ2/KFA/ch3cho.004',STATUS='old')
c         n = 2053
         OPEN(UNIT=kin,FILE='DATAJ2/KFA/ch3cho.005',STATUS='old')
         n = 18
c         OPEN(UNIT=kin,FILE='DATAJ2/KFA/ch3cho.006',STATUS='old')
c         n = 1705

         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
         ENDDO
         CLOSE (kin)


      ELSEIF (mabs .EQ. 5) THEN

         OPEN(UNIT=kin,
     $        FILE='DATAJ1/CH3CHO/CH3CHO_jpl11.abs',STATUS='old')
         do i = 1, 2
            read(kin,*)
         enddo
         n = 101
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
            y1(i) = y1(i) * 1.e-20
         ENDDO
         CLOSE(kin)

      ENDIF

       CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)

* quantum yields

      IF (myld .EQ. 1) THEN

         OPEN(UNIT=kin,FILE='DATAJ1/CH3CHO/CH3CHO_iup.yld',STATUS='old')
         do i = 1, 4
            read(kin,*)
         enddo
         n = 12
         DO i = 1, n
            READ(kin,*) x1(i), y2(i), y1(i)
            x2(i) = x1(i)
         ENDDO
         CLOSE(kin)
         n1 = n
         n2 = n

         CALL interpol(x1,y1,kdata,n1,nw,wl,jlabel(j),0.,0.,yg1)
         CALL interpol(x2,y2,kdata,n2,nw,wl,jlabel(j),0.,0.,yg2)

         DO iw = 1, nw-1
            yg3(iw) = 0.
         ENDDO

      ELSEIF (myld .EQ. 2) THEN

         OPEN(UNIT=kin,FILE='DATAJ1/CH3CHO/d021_i.yld',STATUS='old')
         DO i = 1, 18
            READ(kin,*)
         ENDDO
         n = 10
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
         ENDDO
         CLOSE (kin)

         CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),y1(1),0.,yg1)

         OPEN(UNIT=kin,FILE='DATAJ1/CH3CHO/d021_ii.yld',STATUS='old')
         DO i = 1, 10
            READ(kin,*)
         ENDDO
         n = 9
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
         ENDDO
         CLOSE (kin)

         CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),y1(1),0.,yg2)

         OPEN(UNIT=kin,FILE='DATAJ1/CH3CHO/d021_iii.yld',STATUS='old')
         DO i = 1, 10
            READ(kin,*)
         ENDDO
         n = 9
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
         ENDDO
         CLOSE (kin)

         CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),y1(1),0.,yg3)

       ELSEIF (myld .EQ. 3) THEN

         OPEN(UNIT=kin,FILE='DATAJ1/CH3CHO/CH3CHO_iup13.yld',
     &        STATUS='old')
         do i = 1, 14
            read(kin,*)
         enddo
         n = 17
         DO i = 1, n
            READ(kin,*) x1(i), y2(i), y1(i), y3(i)
            x2(i) = x1(i)
            x3(i) = x1(i)
         ENDDO
         CLOSE(kin)
         n1 = n
         n2 = n
         n3 = n

         CALL interpol(x1,y1,kdata,n1,nw,wl,jlabel(j-2),y1(1),0.,yg1)
         CALL interpol(x2,y2,kdata,n2,nw,wl,jlabel(j-1),y2(1),0.,yg2)
         CALL interpol(x3,y3,kdata,n3,nw,wl,jlabel(j),y3(1),0.,yg3)

      ENDIF


* pressure-dependence

      IF(mpress == 2) THEN

         OPEN(UNIT=kin,FILE='DATAJ1/CH3CHO/CH3CHO_pdep_iup13.yld',
     $     STATUS='old')
         do i = 1, 7
            read(kin,*)
         enddo
         n = 18
         n1= n
         n2= n
         DO i = 1, n
            READ(kin,*) x1(i), y1(i), y2(i)
            x2(i) = x1(i)
            y2(i) = y2(i)*1.E-21
         ENDDO
         CLOSE (kin)

         CALL interpol(x1,y1,kdata,n1,nw,wl,jlabel(j),1.,0.,yg4)
         CALL interpol(x2,y2,kdata,n2,nw,wl,jlabel(j),0.,0.,yg5)

      ENDIF


* combine:

      DO iw = 1, nw - 1

* absorption cross sections:
         sig = yg(iw)

* quantum yields:
         qy1 = yg1(iw)
         qy2 = yg2(iw)
         qy3 = yg3(iw)

* pressure correction for channel 1, CH3 + CHO
* based on Horowitz and Calvert 1982.
* or on overall quantum yield based on IUPAC 2013
         DO i = 1, nz
            IF(mpress==1) THEN

              qy1_n0 = yg1(iw)
              qy1_0  = 1. - qy2 - qy3
              IF (qy1_n0 .GT. 0.) THEN
                x = qy1_0/qy1_n0 - 1.
               ELSE
                x = 0.
              ENDIF

              qy1 = MAX(0.,MIN(1.,qy1_n0 * (1. + x)/(1. + x
     &              *airden(i)/2.465E19)))
!             qy1=MAX(0.,MIN(1.,
!    &            1./(yg5(iw)+yg4(iw)*airden(i)/2.465E19)))
             ELSEIF(mpress==2) THEN
              quench = MAX(0.,MIN(1.,1./(1./yg4(iw)+yg5(iw)*airden(i))))
              qytot = qy1+qy2+qy3
              IF(qytot>0 .AND. yg5(iw)>0.) THEN
                qy1 = qy1/qytot*quench
                qy2 = qy2/qytot*quench
                qy3 = qy3/qytot*quench
              ENDIF
            ENDIF

            sq(j-3,i,iw) = sig * qy1
            sq(j-2,i,iw) = sig * qy2
            sq(j-1,i,iw) = sig * qy3
            sq(j  ,i,iw) = sig

         ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r12(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! C2H5CHO

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for C2H5CHO        =*
*=  photolysis:                                                              =*
*=         C2H5CHO + hv -> C2H5 + HCO                                        =*
*=                                                                           =*
*=  Cross section:  Choice between                                           =*
*=                   (1) IUPAC 97 - 2002 data, from Martinez et al.          =*
*=                   (2) Calvert and Pitts, as tabulated by KFA              =*
*=  Quantum yield:  IUPAC 97 or 2002 recommendation                          =*
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

      REAL yg(kw), yg1(kw), ygoh(kw)
      REAL qy1
      REAL sig, sigoh
      INTEGER iw

      INTEGER  myld, qyvers(3)!, mabs,xsvers(3)
      CHARACTER(llog) :: logmsg(nlog)
      CHARACTER(lspc) :: spc


***************** n-propanal photolysis *************************

* Setting photolysis index and cross section/quantum yield options
      j = j+1
      jlabel(j) = 'C2H5CHO -> C2H5 + HCO'
      j = j+1
      jlabel(j) = 'ALD3OH -> R(OH) + HCO'
      j = j+1
      jlabel(j) = 'genC2H5CHO(poly)'
      j = j+1
      jlabel(j) = 'genC2H5CHO(OHpoly)'

      spc = "C2H5CHO"
* cross section options
C      xsvers = (/1, 1, 2/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
C      logmsg(1) = "from IUPAC '97"
C      logmsg(2) = "from KFA tables"
C      CALL set_option(2,spc,"abs",logmsg,xsvers,mabs)

* quantum yield options
      qyvers = (/1, 2, 2/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "from IUPAC '97"
      logmsg(2) = "from IUPAC '02"
      CALL set_option(2,spc,"yld",logmsg,qyvers,myld)


* cross sections
* second cross sections option removed (pb)
* Data folder/file not available

!      IF (mabs .EQ. 1) THEN

         OPEN(UNIT=kin,FILE='DATAJ1/C2H5CHO/C2H5CHO_iup.abs',
     $        STATUS='old')
         do i = 1, 4
            read(kin,*)
         enddo
         n = 106
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
            y1(i) = y1(i) * 1.e-20
         ENDDO
         CLOSE(kin)

!      ELSEIF(mabs .EQ. 2) THEN

* cross section from KFA tables
* c2h5cho.001 - Calvert and Pitts 1966

!         OPEN(UNIT=kin,FILE='DATAJ2/KFA/c2h5cho.001',STATUS='old')
!         n = 83

!         DO i = 1, n
!            READ(kin,*) x1(i), y1(i)
!         ENDDO
!         CLOSE (kin)

!      ENDIF

      x1oh(:) = 0.
      x1oh(:) = x1(:) - 10.
      y1oh(:) = y1(:)
      n1 = n
      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)
      n = n1
      CALL interpol(x1oh,y1oh,kdata,n,nw,wl,jlabel(j),0.,0.,ygoh)


* quantum yields

      IF (myld .EQ. 1) THEN

         OPEN(UNIT=kin,FILE='DATAJ1/C2H5CHO/C2H5CHO_iup.yld',
     $        STATUS='old')
         do i = 1, 4
            read(kin,*)
         enddo
         n = 5
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
         ENDDO
         CLOSE(kin)

         CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg1)

      ELSEIF (myld .EQ. 2) THEN

         OPEN(UNIT=kin,FILE='DATAJ1/C2H5CHO/C2H5CHO_iup13.yld',
     $        STATUS='old')
         do i = 1, 8
            read(kin,*)
         enddo
         n = 7
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
         ENDDO
         CLOSE(kin)

         CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg1)

      ELSEIF (myld .EQ. 3) THEN

         OPEN(UNIT=kin,FILE='DATAJ1/C2H5CHO/C2H5CHO_calv.yld',
     $        STATUS='old')
         do i = 1, 9
            read(kin,*)
         enddo
         n = 70
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
         ENDDO
         CLOSE(kin)

         CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,1.,yg1)
      ENDIF

* combine:

      DO iw = 1, nw - 1

            sig = yg(iw)
            sigoh = ygoh(iw)

         DO i = 1, nz

* quantum yields:
* use Ster-Volmer pressure dependence:

           IF (yg1(iw) .LT. pzero) THEN
             qy1 = 0.
            ELSE
             qy1 = 1./(1. + (1./yg1(iw) - 1.)*airden(i)/2.45e19)
           ENDIF
           qy1 = MIN(qy1,1.)

           sq(j-3,i,iw) = sig * qy1
           IF(swOH==1) THEN
             sq(j-2,i,iw) = sigoh * qy1
           ELSEIF(swOH==2) THEN
             sq(j-2,i,iw) = sigoh * qyoh
           ENDIF
           sq(j-1,i,iw) = sig
           sq(j  ,i,iw) = sigoh

         ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r13(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! CHOCHO

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for CHOCHO         =*
*=  photolysis:                                                              =*
*=              CHOCHO + hv -> Products                                      =*
*=                                                                           =*
*=  Cross section: Choice between                                            =*
*=                  (1) Plum et al., as tabulated by IUPAC 97                =*
*=                  (2) Plum et al., as tabulated by KFA.                    =*
*=                  (3) Orlando et al.                                       =*
*=                  (4) Horowitz et al., 2001                                =*
*=                  (5) Volkamer et al. 2005 recommended by JPL2011/IUPAC2013=*
*=  Quantum yield:      Salter et al. 2013                                   =*
*=                      (incorrectly referred by IUPAC 2013)                 =*
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
      PARAMETER(kdata=500)

      INTEGER i, n, n1, n2, n3
      REAL x(kdata), x1(kdata)
      REAL y1(kdata), y2(kdata), y3(kdata)

* local

      REAL yg(kw), yg1(kw), yg2(kw), yg3(kw)
      REAL qy, qyI, qyII, qyIII, qy1, qy3
      REAL A1, A2, A3 ! quenching parameters
      REAL sig, dum
      INTEGER iw

      INTEGER mabs, myld, xsvers(3), qyvers(3)
      CHARACTER(llog) :: logmsg(nlog)
      CHARACTER(lspc) :: spc


***************** glyoxal photolysis *************************
* see review by Madronich, Chapter VII in "The Mechansims of
*  Atmospheric Oxidation of the Alkanes, Calvert et al, Oxford U.
*  Press, 2000.
* Four possible channels:
*     I     H2 + 2 CO
*     II    2 HCO
*     III   HCHO + CO
*     IV    HCO + H + CO
*
*  Based on that review, the following quantum yield assignments are made:
*
*     qy_I = 0
*     qy_II = 0.63 for radiation between 280 and 380 nm
*     qy_III = 0.2  for radiation between 280 and 380 nm
*     qy_IV = 0
* The yields for channels II and III were determined by Bauerle et al. (personal
* communication from G. Moortgat, still unpublished as of Dec 2000).
* Bauerle et al. used broad-band irradiation 280-380 nm.
* According to Zhu et al., the energetic threshold (for II) is 417 nm.  Therefore,
* here the quantum yields were set to zero for wc > 417.  Furthermore, the
* qys of Bauerle et al. were reduced to give the same J values when using full solar
* spectrum.  The reduction factor was calculated by comparing the J-values (for
* high sun) using the 380 and 417 cut offs.  The reduction factor is 7.1


* Setting photolysis index and cross section/quantum yield options
      j = j + 1
      jlabel(j) = 'CHOCHO -> 2 HO2 + 2 CO'
      j = j + 1
      jlabel(j) = 'CHOCHO -> H2 + 2 CO'
      j = j + 1
      jlabel(j) = 'CHOCHO -> CH2O + CO'
! generic rate constant for aldehyde groups within di- and polycarbonyls
! (XS = 0.5*XS(GLY), QY = 1)
      j = j + 1
      jlabel(j) = 'Ald (mult)'
      j = j + 1
      jlabel(j) = 'Ald (poly)'

* cross section options
      spc = "CHOCHO"
      xsvers = (/5, 5, 5/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "from IUPAC 97"
      logmsg(2) = "from from KFA tables"
      logmsg(3) = "from Orlando and Tyndall 2001"
      logmsg(4) = "from Horrowitz et al. 2001"
      logmsg(5) = "from from JPL 2011 (also in IUPAC 2013)"
      CALL set_option(5,spc,"abs",logmsg,xsvers,mabs)

* quantum yield options
      qyvers = (/2, 3, 3/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "from IUPAC 97"
      logmsg(2) = "from JPL 2011"
      logmsg(3) = "from Salter et al. 2013"
      CALL set_option(3,spc,"yld",logmsg,qyvers,myld)


* cross sections

      IF (mabs .EQ. 1) THEN

         OPEN(UNIT=kin,FILE='DATAJ1/CHOCHO/CHOCHO_iup.abs',
     $        STATUS='old')
         DO i = 1, 4
            read(kin,*)
         ENDDO
         n = 110
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
            y1(i) = y1(i) * 1.e-20
         ENDDO
         CLOSE(kin)

         CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)


      ELSEIF(mabs .EQ. 2) THEN

* cross section from KFA tables
* chocho.001 - Plum et al. 1983

         OPEN(UNIT=kin,FILE='DATAJ2/KFA/chocho.001',STATUS='old')
         n = 219

         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
         ENDDO
         CLOSE (kin)

         CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)

      ELSEIF(mabs .EQ. 3) THEN

* cross section from Orlando et la.
* Orlando, J. J.; G. S. Tyndall, 2001:  The atmospheric chemistry of the
* HC(O)CO radical. Int. J. Chem. Kinet., 33, 149-156.

         OPEN(UNIT=kin,
     $        FILE='DATAJ1/CHOCHO/glyoxal_orl.abs',STATUS='old')

         do i = 1, 6
            read(kin,*)
         enddo
         n = 481
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
         ENDDO
         CLOSE (kin)

         CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)

      ELSEIF(mabs .EQ. 4) THEN

         OPEN(UNIT=kin,
     $        FILE='DATAJ1/CHOCHO/glyoxal_horowitz.abs',STATUS='old')

         DO i = 1, 8
            read(kin,*)
         ENDDO
         n = 270
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
            y1(i) = y1(i) * 1.e-20
         ENDDO
         CLOSE (kin)

         CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)

      ELSEIF(mabs .eq. 5) then

         open(unit=kin,
     $        FILE='DATAJ1/CHOCHO/glyoxal_jpl11.abs',STATUS='old')

         DO i = 1, 2
            read(kin,*)
         ENDDO
         n = 277
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
            y1(i) = y1(i) * 1.e-20
         ENDDO
         CLOSE (kin)

         CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)

      ENDIF

* quantum yields

      IF(myld .eq. 2) then

         open(unit=kin,
     $        FILE='DATAJ1/CHOCHO/glyoxal_jpl11.qy',STATUS='old')

         DO i = 1, 3
            read(kin,*)
         ENDDO
         n = 40
         DO i = 1, n
            READ(kin,*) x(i), dum, y1(i), y2(i), y3(i)
         ENDDO
         CLOSE (kin)

         n1 = n
         do i = 1, n
            x1(i) = x(i)
         enddo

         CALL interpol(x1,y1,kdata,n1,nw,wl,jlabel(j),y1(1),0.,yg1)
         n2 = n

         do i = 1, n
            x1(i) = x(i)
         enddo
         CALL interpol(x1,y2,kdata,n2,nw,wl,jlabel(j),y2(1),0.,yg2)

         n3 = n
         do i = 1, n
            x1(i) = x(i)
         enddo
         CALL interpol(x1,y3,kdata,n3,nw,wl,jlabel(j),y3(1),0.,yg3)

      ENDIF

* combine:

      DO iw = 1, nw - 1

         sig = yg(iw)

* quantum yields:

         IF(myld .EQ. 1) THEN

* Use values from Bauerle, but corrected to cutoff at 417 rather than 380.
* this correction is a reduction by 7.1.
* so that qyI = 0.63/7.1  and qyII = 0.2/7.1

            qyII = 0.
            if(wc(iw) .lt. 417. ) then
               qyI = 0.089
               qyIII = 0.028
            else
               qyI = 0.
               qyIII = 0.
            endif

          ELSEIF(myld .EQ. 2) THEN

           qyI   = yg1(iw)
           qyII  = yg2(iw)
           qyIII = yg3(iw)

          ELSEIF(myld .EQ. 3) THEN

* calculate wavelength-dependent zero-pressure quantum yields
* from parameterisations given in Salter et al. 2013
* qyI is the sum of their channels P1a, P1b, and P1g

           qy1   = 1.-1./(1.+EXP((wc(iw)-330.)/14.1))
     &           + 0.39-0.39/(1.+EXP((wc(iw)-331.7)/(-14.)))
     &           + 0.28/(1.+EXP((wc(iw)-276.7)/12.))
           qy3   = 1. - qy1
           qyII  = 0.

         ENDIF

         DO i = 1, nz

* calculate pressure dependency of total qy by Salter et al., 2013:
          A1 = 6.48e-19*(tlev(i)/295.)**(-1.83)
     &       * exp((-7.6e-4*(tlev(i)/295.)**(-0.515))
     &       * (1.e7/wc(iw)-2.38e4))
          A2 = 112.8*(tlev(i)/295.)**(-1.53)
     &       * exp((-4.61e-3*(tlev(i)/295.)**0.507)
     &       * (1.e7/wc(iw)-2.38e4))
          A3 = 2.25e-16*(tlev(i)/295.)**(-9.18)
     &       * exp((-7.8e-4*(tlev(i)/295.)**(-7.03))
     &       * (1.e7/wc(iw)-2.38e4))

          qy = min(1.,max(0.,1./
     &          (((1+A1*airden(i)+A2)*(1+A3*airden(i)))/
     &            (1+A2+A3*airden(i)))))
          qyI   = qy*qy1
          qyIII = qy*qy3


          sq(j-4,i,iw) = sig * qyI
          sq(j-3,i,iw) = sig * qyII
          sq(j-2,i,iw) = sig * qyIII
           ! Divide by 2 to use 1 aldehyde group in glyoxal for
           ! compounds with an aldehyde group, an adjacent carbonyl group
           ! and further polyfunctional chromphores
          sq(j-1,i,iw) = sig * qy
          sq(j,  i,iw) = sig

         ENDDO

      ENDDO

      RETURN
      END

*=============================================================================*

      SUBROUTINE r14(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! CH3COCHO

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for CH3COCHO       =*
*=  photolysis:                                                              =*
*=           CH3COCHO + hv -> CH3CO + HCO                                    =*
*=                                                                           =*
*=  Cross section: Choices as given below in the source code                 =*
*=  Quantum yield: Choices as given below in the source code                 =*
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
      PARAMETER(kdata=500)

      INTEGER i, n
      INTEGER n1, n2
      REAL x1(kdata), x2(kdata)
      REAL y1(kdata), y2(kdata)
      real x(kdata), y(kdata)

* local

      REAL yg(kw), yg1(kw), yg2(kw)
      REAL qy
      REAL sig
      INTEGER iw

      REAL phi0, kq

      INTEGER mabs, myld, xsvers(3), qyvers(3)
      CHARACTER(llog) :: logmsg(nlog)
      CHARACTER(lspc) :: spc


***************** methylglyoxal photolysis *************************

* Setting photolysis index and cross section/quantum yield options
      j = j+1
      jlabel(j) = 'CH3COCHO -> CH3CO + HCO'

* cross section options
* 1:  from Meller et al. (1991), as tabulated by IUPAC
*     for wc < 402, use coarse data (5 nm, table 1)
*     for wc > 402, use finer data (2 nm, table 2)
* 2: average at 1nm of  Staffelbach et al. 1995 and Meller et al. 1991
*     Cross section from KFA tables:
* 3: ch3cocho.001 - Plum et al. 1983
* 4: ch3cocho.002 - Meller et al. 1991, 0.033 nm resolution
* 5: ch3cocho.003 - Meller et al. 1991, 1.0   nm resolution
* 6: ch3cocho.004 - Staffelbach et al. 1995
* 7: use synthetic spectrum, average of CHOCHO and CH3COCOCH3:
* 8: cross section from JPL2011
      spc = "CH3COCHO"
      xsvers = (/8, 1, 7/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "from IUPAC"
      logmsg(2) =
     &          "mean of Staffelbach et al. 1995 and Meller et al. 1991"
      logmsg(3) = "from Plum et al. 1983"
      logmsg(4) = "in high res. from Meller et al. 1991"
      logmsg(5) = "from Meller et al. 1991"
      logmsg(6) = "from Staffelbach et al. 1995"
      logmsg(7) = "synthesised from glyoxal and biacetyl cross sections"
      logmsg(8) = "from JPL 2011"
      CALL set_option(8,spc,"abs",logmsg,xsvers,mabs)

* quantum yield options
* 1:  Plum et al., 0.107
* 2:  Plum et al., divided by two = 0.0535
* 3:  Staffelbach et al., 0.45 at wc .le. 300, 0 for wc .gt. 430, linear
*     interpl in between
* 4:  Koch and Moortgat, prv. comm. 1997. - pressure-dependent
* 5:  Chen, Y., W. Wang, and L. Zhu, Wavelength-dependent photolysis of methylglyoxal
*      in the 290-440 nm region, J Phys Chem A, 104, 11126-11131, 2000.
* 6:  Quantum yields from IUPAC (2003) in N2
* 7:  Quantum yields from IUPAC (2003) assumed as air
      qyvers = (/5, 6, 5/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "from Plum et al"
      logmsg(2) = "from Plum et al. divided by 2"
      logmsg(3) = "from Staffelbach et al"
      logmsg(4) = "from from Koch and Moortgat"
      logmsg(5) = "from Chen, Y. et al"
      logmsg(6) = "from IUPAC in N2"
      CALL set_option(6,spc,"yld",logmsg,qyvers,myld)


* cross sections

      IF (mabs .EQ. 1) THEN

         OPEN(UNIT=kin,FILE='DATAJ1/CH3COCHO/CH3COCHO_iup1.abs',
     $        STATUS='old')
         do i = 1, 4
            read(kin,*)
         enddo
         n = 38
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
            y1(i) = y1(i) * 1.e-20
         ENDDO
         CLOSE(kin)

         CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg1)

         OPEN(UNIT=kin,FILE='DATAJ1/CH3COCHO/CH3COCHO_iup2.abs',
     $        STATUS='old')
         do i = 1, 4
            read(kin,*)
         enddo
         n = 75
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
            y1(i) = y1(i) * 1.e-20
         ENDDO
         CLOSE(kin)

         CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),y1(1),0.,yg2)

         DO iw = 1, nw-1
            IF(wc(iw) .LT. 401.) THEN
               yg(iw) = yg1(iw)
            ELSE
               yg(iw) = yg2(iw)
            ENDIF
         ENDDO

      ELSEIF(mabs .EQ. 2) THEN

         OPEN(UNIT=kin,FILE='DATAJ1/CH3COCHO/CH3COCHO_ncar.abs',
     $        STATUS='old')
         n = 271
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
         ENDDO
         CLOSE(kin)

         CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)

      ELSEIF(mabs .GT. 2 .and. mabs .lt. 7) THEN

* cross section from KFA tables
* ch3cocho.001 - Plum et al. 1983
* ch3cocho.002 - Meller et al. 1991, 0.033 nm resolution
* ch3cocho.003 - Meller et al. 1991, 1.0   nm resolution
* ch3cocho.004 - Staffelbach et al. 1995

         IF(mabs .EQ. 3) THEN
            OPEN(UNIT=kin,FILE='DATAJ2/KFA/ch3cocho.001',STATUS='old')
            n = 136
         ELSEIF(mabs .EQ. 4) THEN
            OPEN(UNIT=kin,FILE='DATAJ2/KFA/ch3cocho.002',STATUS='old')
            n = 8251
         ELSEIF(mabs .EQ. 5) THEN
            OPEN(UNIT=kin,FILE='DATAJ2/KFA/ch3cocho.003',STATUS='old')
            n = 275
         ELSEIF(mabs .EQ. 6) THEN
            OPEN(UNIT=kin,FILE='DATAJ2/KFA/ch3cocho.004',STATUS='old')
            n = 162
         ENDIF

         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
         ENDDO
         CLOSE (kin)

         CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)

      ELSEIF(mabs .EQ. 7) THEN

         OPEN(UNIT=kin,FILE='DATAJ1/CH3COCOCH3/biacetyl_plum.abs',
     $        STATUS='old')
         DO i = 1, 7
            READ(kin,*)
         ENDDO
         n = 55
         DO i = 1, n
            READ(kin,*) x(i), y(i)
            y(i) = y(i) * 1.e-20
         ENDDO
         CLOSE(kin)

         CALL interpol(x,y,kdata,n,nw,wl,jlabel(j),0.,0.,yg1)


         OPEN(UNIT=kin,
     $        FILE='DATAJ1/CHOCHO/glyoxal_orl.abs',STATUS='old')
         do i = 1, 6
            read(kin,*)
         enddo
         n = 481
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
         ENDDO
         CLOSE (kin)

         CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg2)

         do iw = 1, nw-1
            yg(iw) = 0.5*(yg1(iw) + yg2(iw))
         enddo

      ELSEIF(mabs .eq. 8) then

         OPEN(UNIT=kin,
     $        FILE='DATAJ1/CH3COCHO/CH3COCHO_jpl11.abs',STATUS='old')
         do i = 1, 2
            read(kin,*)
         enddo
         n = 294
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
            y1(i) = y1(i) * 1.e-20
         ENDDO
         CLOSE (kin)

         CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)


      ENDIF

* quantum yields

         IF(myld .EQ. 4) THEN
            OPEN(UNIT=kin,FILE='DATAJ1/CH3COCHO/CH3COCHO_km.yld',
     $           STATUS='old')
            DO i = 1, 5
               READ(kin,*)
            ENDDO
            n = 5
            DO i = 1, n
               READ(kin,*) x1(i), y1(i), y2(i)
               x2(i) = x1(i)
            ENDDO
            CLOSE (kin)
            n1 = n
            n2 = n

            CALL interpol(x1,y1,kdata,n1,nw,wl,jlabel(j),1.,0.,yg1)
            CALL interpol(x2,y2,kdata,n2,nw,wl,jlabel(j),1.,0.,yg2)

         ENDIF


* combine:

      DO iw = 1, nw - 1

         sig = yg(iw)

         DO i = 1, nz

* quantum yields:

            IF (myld .EQ. 1) THEN
               qy = 0.107

            ELSEIF(myld .EQ. 2) THEN
               qy = 0.107/2.

            ELSEIF(myld .EQ. 3) THEN
               IF(wc(iw) .LE. 300.) THEN
                  qy = 0.45
               ELSE IF (wc(iw) .GE. 430.) THEN
                  qy = 0.
               ELSE
                  qy = 0.45 + (0-0.45)*(wc(iw)-300.)/(430.-300.)
               ENDIF

            ELSEIF(myld .EQ. 4) THEN

               IF (yg1(iw) .GT. 0.) THEN

                  qy = yg2(iw)/( 1. + (airden(i)/2.465E19)
     $                 * ( (yg2(iw)/yg1(iw)) - 1.))

               ELSE
                  qy = 0.
               ENDIF

            ELSEIF(myld .EQ. 5) THEN

* zero pressure yield:
* 1.0 for wc < 380 nm
* 0.0 for wc > 440 nm
* linear in between:

               phi0 = 1. - (wc(iw) - 380.)/60.
               phi0 = MIN(phi0,1.)
               phi0 = MAX(phi0,0.)

* Pressure correction: quenching coefficient, torr-1
* in air, Koch and Moortgat:

               kq = 1.36e8 * EXP(-8793/wc(iw))

* in N2, Chen et al:

c               kq = 1.93e4 * EXP(-5639/wc(iw))

               IF(phi0 .GT. 0.) THEN
                  IF (wc(iw) .GE. 380. .AND. wc(iw) .LE. 440.) THEN
                     qy = phi0 / (phi0 + kq * airden(i) * 760./2.456E19)
                  ELSE
                     qy = phi0
                  ENDIF
               ELSE
                  qy = 0.
               ENDIF

            ELSEIF(myld .EQ. 6) THEN

             phi0 = 1.
             IF(wc(iw) .GE. 380. .AND. wc(iw) .LE. 440.) THEN
               phi0 = MAX(0.,MIN(1.,3.63e-7*exp(5693./wc(iw))))
              ELSEIF(wc(iw) .GT. 440.) THEN
               phi0 = 0.
             ENDIF

             kq = 1.93e4*exp(-5639./wc(iw))

             IF(phi0>0.) THEN
               qy = 1. / (1./phi0 + kq * airden(i) * 760./2.456E19)
              ELSE
               qy = 0.
             ENDIF

            ENDIF

            sq(j,i,iw) = sig * qy

         ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r15(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! CH3COCH3

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for CH3COCH3 photolysis=*
*=          CH3COCH3 + hv -> Products                                        =*
*=                                                                           =*
*=  Cross section: Choices as given below in the source code                 =*
*=  Quantum yield: Choices as given below in the source code                 =*
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
      INTEGER n1, n2, n3, n4
      REAL x1(kdata), x2(kdata), x3(kdata), x4(kdata)
      REAL y1(kdata), y2(kdata), y3(kdata), y4(kdata)

* local

      REAL yg(kw), yg1(kw), yg2(kw), yg3(kw), yg4(kw)
      REAL qy1, qy2
      REAL sig
      INTEGER iw

      REAL a, b, t, m, w
      real fco, fac

      INTEGER mabs, myld, xsvers(3), qyvers(3)
      CHARACTER(llog) :: logmsg(nlog)
      CHARACTER(lspc) :: spc


***************** acetone photolysis *************************

* Setting photolysis index and cross section/quantum yield options
      j = j + 1
      jlabel(j) = 'CH3COCH3 -> CH3CO + CH3'
      j = j + 1
      jlabel(j) = 'CH3COCH3 -> CO + 2 CH3'

* cross section options
* 1:  cross section from Calvert and  Pitts
* 2:  Martinez et al. 1991, also in IUPAC'97
* 3:  NOAA 1998, unpublished as of Jan 98.
* 4:  JPL-2011/IUPAC 2013
      spc = "CH3COCH3"
      xsvers = (/4, 4, 4/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "from Calvert and Pitts"
      logmsg(2) = "from IUPAC 1997"
      logmsg(3) = "from NOAA 1998"
      logmsg(4) = "JPL 2011 / IUPAC 2013"
      CALL set_option(4,spc,"abs",logmsg,xsvers,mabs)

* quantum yield options
* 1:  Gardiner et al. 1984
* 2:  IUPAC 97
* 3:  McKeen, S. A., T. Gierczak, J. B. Burkholder, P. O. Wennberg, T. F. Hanisco,
*       E. R. Keim, R.-S. Gao, S. C. Liu, A. R. Ravishankara, and D. W. Fahey,
*       The photochemistry of acetone in the upper troposphere:  a source of
*       odd-hydrogen radicals, Geophys. Res. Lett., 24, 3177-3180, 1997.
* 4:  Blitz, M. A., D. E. Heard, M. J. Pilling, S. R. Arnold, and M. P. Chipperfield
*       (2004), Pressure and temperature-dependent quantum yields for the
*       photodissociation of acetone between 279 and 327.5 nm, Geophys.
*       Res. Lett., 31, L06111, doi:10.1029/2003GL018793.
      qyvers = (/4, 4, 4/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "from Gardiner et al. 1984"
      logmsg(2) = "from IUPAC 97"
      logmsg(3) = "from McKeen et al. 1997"
      logmsg(4) = "from Blitz et al. 2004 (IUPAC 2013)"
      CALL set_option(4,spc,"yld",logmsg,qyvers,myld)


* cross sections

      IF (mabs .EQ. 1) THEN

         OPEN(UNIT=kin,FILE='DATAJ1/CH3COCH3/CH3COCH3_cp.abs',
     $        STATUS='old')
         DO i = 1, 6
            READ(kin,*)
         ENDDO
         n = 35
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
            y1(i) = y1(i) * 3.82E-21
         ENDDO
         CLOSE (kin)

         CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j-1),0.,0.,yg)

      ELSEIF(mabs .EQ. 2) THEN

         OPEN(UNIT=kin,FILE='DATAJ1/CH3COCH3/CH3COCH3_iup.abs',
     $        STATUS='old')
         DO i = 1, 4
            READ(kin,*)
         ENDDO
         n = 96
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
            y1(i) = y1(i) * 1.e-20
         ENDDO
         CLOSE (kin)

         CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j-1),0.,0.,yg)

      ELSEIF(mabs .EQ. 3) THEN

         OPEN(UNIT=kin,FILE='DATAJ1/CH3COCH3/CH3COCH3_noaa.abs',
     $        STATUS='old')
         DO i = 1, 12
            READ(kin,*)
         ENDDO
         n = 135
         DO i = 1, n
            READ(kin,*) x1(i), y1(i), y2(i), y3(i)
            x2(i) = x1(i)
            x3(i) = x1(i)
         ENDDO
         CLOSE (kin)
         n1 = n
         n2 = n
         n3 = n

         CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j-1),0.,0.,yg)
         CALL interpol(x2,y2,kdata,n2,nw,wl,jlabel(j-1),0.,0.,yg2)
         CALL interpol(x3,y3,kdata,n3,nw,wl,jlabel(j-1),0.,0.,yg3)

      ELSEIF(mabs.eq.4) then
         OPEN(UNIT=kin,FILE='DATAJ1/CH3COCH3/CH3COCH3_jpl11.abs',
     $        STATUS='old')
         DO i = 1, 5
            READ(kin,*)
         ENDDO
         n = 135
         DO i = 1, n
            READ(kin,*) x1(i), y1(i), y2(i), y3(i), y4(i)
            x2(i) = x1(i)
            x3(i) = x1(i)
            x4(i) = x1(i)
            y1(i) = y1(i) * 1.e-20
            y2(i) = y2(i) / 1.e3
            y3(i) = y3(i) / 1.e5
            y4(i) = y4(i) / 1.e8
         ENDDO
         CLOSE (kin)
         n1 = n
         n2 = n
         n3 = n
         n4 = n

         CALL interpol(x1,y1,kdata,n1,nw,wl,jlabel(j-1),0.,0.,yg)
         CALL interpol(x2,y2,kdata,n2,nw,wl,jlabel(j-1),0.,0.,yg2)
         CALL interpol(x3,y3,kdata,n3,nw,wl,jlabel(j-1),0.,0.,yg3)
         CALL interpol(x4,y4,kdata,n4,nw,wl,jlabel(j-1),0.,0.,yg4)
      ENDIF


* quantum yields

      IF (myld .EQ. 2) THEN

         OPEN(UNIT=kin,FILE='DATAJ1/CH3COCH3/CH3COCH3_iup.yld',
     $        STATUS='old')
         DO i = 1, 4
            READ(kin,*)
         ENDDO
         n = 9
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
         ENDDO
         CLOSE (kin)

         CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j-1),0.,0.,yg1)

      ENDIF


* combine

      DO iw = 1, nw - 1

         DO i = 1, nz

            sig = yg(iw)

            IF(mabs .EQ. 3) THEN
*!!! this definition of t is not consistent with JPL2011

               t = 298. - tlev(i)
               t = MIN(t, 298.-235.)
               t = MAX(t, 0.)

               sig = yg(iw)*(1. + yg2(iw)*t + yg3(iw)*t*t)

            ELSEIF (mabs .eq. 4)THEN

               t = MIN(MAX(tlev(i), 235.),298.)
               sig = yg(iw)*(1. + yg2(iw)*t + yg3(iw)*t*t)

            ENDIF

            IF (myld .EQ. 1) THEN
               qy1 = 0.0766 + 0.09415*EXP(-airden(i)/3.222e18)
               qy2 = 0.

            ELSEIF (myld .EQ. 2) THEN
               qy1 = yg1(iw)
               qy2 = 0.

            ELSEIF (myld .EQ. 3) THEN
               IF (wc(iw) .LE. 292.) THEN
                  qy1 = 1.
                  qy2 = 0.
               ELSEIF (wc(iw) .GE. 292.  .AND. wc(iw) .LT. 308. ) THEN
                  a = -15.696 + 0.05707*wc(iw)
                  b = EXP(-88.81+0.15161*wc(iw))
                  qy1 = 1./(a + b*airden(i))
                  qy2 = 0.
               ELSEIF (wc(iw) .GE. 308.  .AND. wc(iw) .LT. 337. ) THEN
                  a = -130.2 + 0.42884*wc(iw)
                  b = EXP(-55.947+0.044913*wc(iw))
                  qy1 = 1./(a + b*airden(i))
                  qy2 = 0.
               ELSEIF (wc(iw) .GE. 337.) THEN
                  qy1 = 0.
                  qy2 = 0.
               ENDIF

               qy1 = max(0., qy1)
               qy1 = min(1., qy1)

            ELSEIF (myld .eq. 4) then
                  w = wc(iw)
                  t = tlev(i)
                  m = airden(i)
                  CALL qyacet(w, t, m, fco, fac)
                  qy1 = min(1., max(0.,fac))
                  qy2 = min(1., max(0.,fco))

            ENDIF

            sq(j-1,i,iw) = sig*qy1
            sq(j,i,iw)   = sig*qy2

         ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r16(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! CH3OOH

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for CH3OOH photolysis: =*
*=         CH3OOH + hv -> CH3O + OH                                          =*
*=                                                                           =*
*=  Cross section: Choices as given below in the source code                 =*
*=  Quantum yield: unity (IUPAC)                                             =*
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
      PARAMETER(kdata=100)

      INTEGER i, n, idum
      REAL x1(kdata)
      REAL y1(kdata)

* local

      REAL yg(kw)
!      REAL qy
      INTEGER iw

      INTEGER mabs, xsvers(3)
      CHARACTER(llog) :: logmsg(nlog)
      CHARACTER(lspc) :: spc


***************** CH3OOH photolysis *************************

* Setting photolysis index and cross section options
* 1:  JPL data base (1985,92,94,97). 1997 is from  Vaghjiani and Ravishankara (1989), at
*     10 nm resolution
* 2:  IUPAC97 (from  Vaghjiani and Ravishankara (1989) at 5 nm resolution).
* 3:  Cox and Tyndall (1978), only for wavelengths < 280 nm
* 4:  Molina and Arguello (1979).  According to Vaghjiani and Ravishankara (1989),
*     Molina and Arguello had a problem measuring CH3OOH, cross sections 40% too high.
* 5:  JPL2011
         j = j + 1
         jlabel(j) = 'CH3OOH -> CH3O + OH'

* cross section options
      spc = "CH3OOH"
      xsvers = (/5, 2, 2/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "from JPL 85 - 02"
      logmsg(2) = "from IUPAC 1997"
      logmsg(3) = "from Cox and Tyndall 1978"
      logmsg(4) = "from Molina and Arguello (1979)"
      logmsg(5) = "from JPL 2011"
      CALL set_option(5,spc,"abs",logmsg,xsvers,mabs)


* cross sections

      IF (mabs .EQ. 1) THEN

c         OPEN(UNIT=kin,FILE='DATAJ1/CH3OOH/CH3OOH_jpl85.abs',
c     $        STATUS='old')
c         OPEN(UNIT=kin,FILE='DATAJ1/CH3OOH/CH3OOH_jpl92.abs',
c     $        STATUS='old')
c         OPEN(UNIT=kin,FILE='DATAJ1/CH3OOH/CH3OOH_jpl94.abs',
c     $        STATUS='old')
         OPEN(UNIT=kin,FILE='DATAJ1/CH3OOH/CH3OOH_jpl94.abs',
     $        STATUS='old')
         READ(kin,*) idum, n
         DO i = 1, idum-2
            READ(kin,*)
         ENDDO
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
            y1(i) = y1(i) * 1.E-20
         ENDDO
         CLOSE (kin)

         CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)

      ELSEIF (mabs .EQ. 2) THEN

         OPEN(UNIT=kin,FILE='DATAJ1/CH3OOH/CH3OOH_iup.abs',
     $        STATUS='old')
         DO i = 1, 4
            READ(kin,*)
         ENDDO
         n = 32
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
            y1(i) = y1(i) * 1.E-20
         ENDDO
         CLOSE (kin)

         CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)

      ELSEIF (mabs .EQ. 3) THEN

         OPEN(UNIT=kin,FILE='DATAJ1/CH3OOH/CH3OOH_ct.abs',
     $        STATUS='old')
         DO i = 1, 4
            READ(kin,*)
         ENDDO
         n = 12
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
         ENDDO
         CLOSE (kin)

         CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)

      ELSEIF (mabs .EQ. 4) THEN

         OPEN(UNIT=kin,FILE='DATAJ1/CH3OOH/CH3OOH_ma.abs',
     $        STATUS='old')
         DO i = 1, 4
            READ(kin,*)
         ENDDO
         n = 15
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
         ENDDO
         CLOSE (kin)

         CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)

      ELSEIF (mabs .EQ. 5) THEN

         OPEN(UNIT=kin,FILE='DATAJ1/CH3OOH/CH3OOH_jpl11.abs',
     $        STATUS='old')
         DO i = 1, 2
            READ(kin,*)
         ENDDO
         n = 40
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
            y1(i) = y1(i) * 1.e-20
         ENDDO
         CLOSE (kin)

         CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)

      ENDIF


* quantum yield = 1
!      qy = 1.

* Map onto grid
      DO iw = 1, nw - 1

         DO i = 1, nz
            sq(j,i,iw) = yg(iw)!*qy
         ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r17(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! CH3NO3

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for CH3ONO2            =*
*=  photolysis:                                                              =*
*=          CH3ONO2 + hv -> CH3O + NO2                                       =*
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
      INTEGER n1, n2
      REAL x1(kdata), x2(kdata)
      REAL y1(kdata), y2(kdata)

* local

      REAL yg(kw), yg1(kw)
!      REAL qy
      REAL sig

      INTEGER mabs, xsvers(3)
      CHARACTER(llog) :: logmsg(nlog)
      CHARACTER(lspc) :: spc


***************** CH3ONO2 photolysis *************************

* Setting photolysis index and cross section options
* 1:  Calvert and  Pitts 1966
* 2:  Talukdar, Burkholder, Hunter, Gilles, Roberts, Ravishankara, 1997.
* 3:  IUPAC-97, table of values for 298K.
* 4:  IUPAC-97, temperature-dependent equation
* 5:  Taylor et al. 1980
* 6:  fit from Roberts and Fajer, 1989
* 7:  Rattigan et al. 1992
* 8:  Libuda and Zabel 1995
* 9:  JPL2011 including T-dependence
*10:  IUPAC (2002) recommendation
      j = j + 1
      jlabel(j) = 'CH3ONO2 -> CH3O + NO2'

* cross section options
      spc = "CH3ONO2"
      xsvers = (/9, 10, 9/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg( 1) = "from Calvert and Pitts (1966)"
      logmsg( 2) = "from Talukdar et al. (1997)"
      logmsg( 3) = "at 298K from IUPAC 1997"
      logmsg( 4) = "from IUAPC 1997 with T-dep"
      logmsg( 5) = "from Taylor et al. (1980)"
      logmsg( 6) = "from fitted Roberts and Fajer 1989 data"
      logmsg( 7) = "from Rattigan et al. (1992)"
      logmsg( 8) = "from Libuda and Zabel (1995)"
      logmsg( 9) = "from JPL 2011"
      logmsg(10) = "from IUPAC 2002"
      CALL set_option(10,spc,"abs",logmsg,xsvers,mabs)


* cross sections

      IF (mabs .EQ. 1) THEN

         OPEN(UNIT=kin,FILE='DATAJ1/RONO2/CH3ONO2_cp.abs',STATUS='old')
         DO i = 1, 3
            READ(kin,*)
         ENDDO
         n = 15
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
         ENDDO
         CLOSE (kin)

         n1 = n
         CALL interpol(x1,y1,kdata,n1,nw,wl,jlabel(j),0.,0.,yg)

      ELSEIF(mabs .EQ. 2) THEN

*        sigma(T,lambda) = sigma(298,lambda) * exp(B * (T-298))

         OPEN(UNIT=kin,FILE='DATAJ1/RONO2/CH3ONO2_tal.abs',STATUS='old')
         DO i = 1, 4
            READ(kin,*)
         ENDDO
         n = 55
         DO i = 1, n
            READ(kin,*) x1(i), y1(i), y2(i)
            x2(i) = x1(i)
            y1(i) = y1(i) * 1.e-20
         ENDDO
         CLOSE (kin)

         n1 = n
         n2 = n
         CALL interpol(x1,y1,kdata,n1,nw,wl,jlabel(j),0.,0.,yg)
         CALL interpol(x2,y2,kdata,n2,nw,wl,jlabel(j),y2(1),y2(n2),yg1)

      ELSEIF (mabs .EQ. 3) THEN

         OPEN(UNIT=kin,FILE='DATAJ1/RONO2/CH3ONO2_iup1.abs',
     $        STATUS='old')
         DO i = 1, 4
            READ(kin,*)
         ENDDO
         n = 13
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
            y1(i) = y1(i)*1e-20
         ENDDO
         CLOSE (kin)

         n1 = n
         CALL interpol(x1,y1,kdata,n1,nw,wl,jlabel(j),0.,0.,yg)

      ELSEIF(mabs .EQ. 4) THEN

*        sigma(T,lambda) = sigma(298,lambda) * 10**(B * T)

         OPEN(UNIT=kin,FILE='DATAJ1/RONO2/CH3ONO2_iup2.abs',
     $        STATUS='old')
         DO i = 1, 4
            READ(kin,*)
         ENDDO
         n = 7
         DO i = 1, n
            READ(kin,*) x1(i), y1(i), y2(i)
            x2(i) = x1(i)
            y1(i) = y1(i) * 1.e-21
            y2(i) = y2(i) * 1.e-3
         ENDDO
         CLOSE (kin)

         n1 = n
         n2 = n
         CALL interpol(x1,y1,kdata,n1,nw,wl,jlabel(j),-36.,-36.,yg)
         CALL interpol(x2,y2,kdata,n2,nw,wl,jlabel(j),y2(1),y2(n2),yg1)

      ELSEIF (mabs .EQ. 5) THEN

         OPEN(UNIT=kin,FILE='DATAJ1/RONO2/CH3ONO2_tay.abs',
     $        STATUS='old')
         DO i = 1, 4
            READ(kin,*)
         ENDDO
         n = 13
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
         ENDDO
         CLOSE (kin)

         n1 = n
         CALL interpol(x1,y1,kdata,n1,nw,wl,jlabel(j),0.,0.,yg)

      ELSEIF (mabs .EQ. 6) THEN

         DO iw = 1, nw-1
            IF(wc(iw) .GT. 284.) THEN
               yg(iw) = EXP(-1.044e-3*wc(iw)*wc(iw) +
     $              0.5309*wc(iw) - 112.4)
            ELSE
               yg(iw) = 0.
            ENDIF
         ENDDO

      ELSEIF (mabs .EQ. 7) THEN

         OPEN(UNIT=kin,FILE='DATAJ1/RONO2/CH3ONO2_rat.abs',
     $        STATUS='old')
         DO i = 1, 4
            READ(kin,*)
         ENDDO
         n = 24
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
         ENDDO
         CLOSE (kin)

         n1 = n
         CALL interpol(x1,y1,kdata,n1,nw,wl,jlabel(j),0.,0.,yg)

      ELSEIF (mabs .EQ. 8) THEN

         OPEN(UNIT=kin,FILE='DATAJ1/RONO2/CH3ONO2_lib.abs',
     $        STATUS='old')
         DO i = 1, 4
            READ(kin,*)
         ENDDO
         n = 1638
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
         ENDDO
         CLOSE (kin)

         n1 = n
         CALL interpol(x1,y1,kdata,n1,nw,wl,jlabel(j),0.,0.,yg)

      ELSEIF (mabs. eq. 9) THEN

         OPEN(UNIT=kin,FILE='DATAJ1/RONO2/CH3ONO2_jpl11.abs',
     $        STATUS='old')
         DO i = 1, 2
            READ(kin,*)
         ENDDO
         n = 65
         DO i = 1, n
            READ(kin,*) x1(i), y1(i), y2(i)
            y1(i) = y1(i) * 1.e-20
            x2(i) = x1(i)
            y2(i) = y2(i) * 1.e-3
         ENDDO
         CLOSE (kin)

         n1 = n
         CALL interpol(x1,y1,kdata,n1,nw,wl,jlabel(j),0.,0.,yg)
         n2 = n
         CALL interpol(x2,y2,kdata,n2,nw,wl,jlabel(j),0.,0.,yg1)

      ELSEIF (mabs .eq. 10) THEN

         OPEN(UNIT=kin,FILE='DATAJ1/RONO2/CH3NO3_iup.abs',
     $        STATUS='old')
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
     $        STATUS='old')
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

      ENDIF

* quantum yield = 1
!      qy = 1.

* combine
      DO iw = 1, nw - 1
         sig = yg(iw)

         DO i = 1, nz

            IF(mabs .EQ. 2 .OR. mabs .EQ. 9 .OR. mabs .EQ. 10) THEN
               sig = yg(iw) * exp (yg1(iw) * (tlev(i)-298.))
            ELSEIF (mabs .EQ. 4) THEN
               sig = yg(iw)*10.**(yg1(iw)*tlev(i))
            ENDIF

            sq(j,i,iw) = sig! * qy

         ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r18(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! PAN

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for PAN photolysis:    =*
*=       PAN + hv -> Products                                                =*
*=                                                                           =*
*=  Cross section: from Talukdar et al., 1995                                =*
*=  Quantum yield: see below in source code                                  =*
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
      PARAMETER(kdata=100)

      INTEGER iw
      INTEGER i, n, n1, n2
      REAL x1(kdata), x2(kdata)
      REAL y1(kdata), y2(kdata)

* local

      REAL yg(kw), yg1(kw), yg2(kw), yg3(kw)
      real qyNO2, qyNO3
      REAL sig

      INTEGER myld, qyvers(3)
      CHARACTER(llog) :: logmsg(nlog)
      CHARACTER(lspc) :: spc


***************** PAN photolysis *************************

* Setting photolysis index and quantum yield options
* 1: from JPL 2011 values for >300 nm.
* 2: from IUPAC 2005 (linearly interpolated)
* 3: from Calvert et al. 2011 book (ammended by values at 248nm from IUPAC)
      j = j+1
      jlabel(j) = 'CH3CO(OONO2) -> CH3CO(OO) + NO2'
      j = j+1
      jlabel(j) = 'CH3CO(OONO2) -> CH3CO(O) + NO3'

* quantum yield options
      spc = "PAN"
      qyvers = (/1, 2, 2/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "from JPL"
      logmsg(2) = "from IUPAC"
      logmsg(3) = "from the Calvert et al. 2011 book"
      CALL set_option(3,spc,"yld",logmsg,qyvers,myld)


* cross sections

* cross section from Senum et al., 1984, J.Phys.Chem. 88/7, 1269-1270

C     OPEN(UNIT=kin,FILE='DATAJ1/RONO2/PAN_senum.abs',STATUS='OLD')
C     DO i = 1, 14
C        READ(kin,*)
C     ENDDO
C     n = 21
C     DO i = 1, n
C        READ(kin,*) x1(i), y1(i)
C        y1(i) = y1(i) * 1.E-20
C     ENDDO
C     CLOSE(kin)

C      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
C      CALL addpnt(x1,y1,kdata,n,               0.,0.)
C      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
C      CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
C      CALL inter2(nw,wl,yg,n,x1,y1,ierr)
C      IF (ierr .NE. 0) THEN
C         WRITE(*,*) ierr, jlabel(j)
C         STOP
C      ENDIF

* cross section from
*      Talukdar et al., 1995, J.Geophys.Res. 100/D7, 14163-14174

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

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j-1),0.,0.,yg)
      CALL interpol(x2,y2,kdata,n2,nw,wl,jlabel(j),0.,0.,yg1)


* quantum yields

      IF(myld==3) THEN

        OPEN(UNIT=kin,FILE='DATAJ1/RONO2/PAN_calv.yld',STATUS='OLD')
        DO i = 1, 8
          READ(kin,*)
        ENDDO
        n = 74
        n1 = n
        n2 = n
        DO i = 1, n
          READ(kin,*) x1(i), y1(i), y2(i)
          x2(i) = x1(i)
        ENDDO
        CLOSE(kin)

        CALL interpol(x1,y1,kdata,n1,nw,wl,jlabel(j),y1(1),y1(n1),yg2)
        CALL interpol(x2,y2,kdata,n2,nw,wl,jlabel(j),y2(1),y2(n2),yg3)

      ENDIF


      DO iw = 1, nw-1
        IF (myld == 1) THEN
          qyNO2 = 0.7
          qyNO3 = 0.3
         ELSEIF (myld == 2) THEN
          IF(wc(iw).LE.248.) THEN
            qyNO2 = 0.76
            qyNO3 = 0.24
           ELSEIF(wc(iw).GE.308.) THEN
            qyNO2 = 0.59
            qyNO3 = 0.41
           ELSE
            qyNO2 = 0.76-(wc(iw)-248.)/(308.-248.)*(0.76-0.59)
            qyNO3 = 1 - qyNO2
          ENDIF
         ELSEIF(myld == 3) THEN
          qyNO2 = yg2(iw)
          qyNO3 = yg3(iw)
        ENDIF


* combine

        DO i = 1, nz

          sig = yg(iw) * EXP(yg1(iw)*(tlev(i)-298.))

          sq(j-1,i,iw) = qyNO2 * sig
          sq(j,i,iw)   = qyNO3 * sig

        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r19(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! CCl2O

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for CCl2O photolysis:  =*
*=        CCl2O + hv -> Products                                             =*
*=                                                                           =*
*=  Cross section: JPL 94 recommendation                                     =*
*=  Quantum yield: Unity (Calvert and Pitts)                                 =*
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
      PARAMETER(kdata=100)

      REAL x1(kdata)
      REAL y1(kdata)

* local

      REAL yg(kw)
      REAL qy
      INTEGER i, iw, n, idum
      INTEGER iz

***************** CCl2O photolysis *************************

      j = j+1
      jlabel(j) = 'CCl2O -> Products'

*** cross sections from JPL94 recommendation

      OPEN(kin,FILE='DATAJ1/ABS/CCl2O_jpl94.abs',STATUS='OLD')
      READ(kin,*) idum, n
      DO i = 1, idum-2
        READ(kin,*)
      ENDDO
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
        y1(i) = y1(i) * 1E-20
      ENDDO
      CLOSE(kin)

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)

*** quantum yield unity (Calvert and Pitts)
      qy = 1.
      DO iw = 1, nw-1
        DO iz = 1, nz
           sq(j,iz,iw) = qy * yg(iw)
        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r20(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! CCl4

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for CCl4 photolysis:   =*
*=      CCl4 + hv -> Products                                                =*
*=  Cross section: from JPL 97 recommendation                                =*
*=  Quantum yield: assumed to be unity                                       =*
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
      PARAMETER(kdata=100)

      REAL x1(kdata)
      REAL y1(kdata)

* local

      REAL yg(kw)
      REAL qy
      INTEGER i, iw, n, idum
      INTEGER iz
      REAL b0, b1, b2, b3, b4, tcoeff, sig
      REAL w1, w2, w3, w4, temp

      INTEGER mabs, xsvers(3)
      CHARACTER(llog) :: logmsg(nlog)
      CHARACTER(lspc) :: spc


***************** CCl4 photolysis *************************

* Setting photolysis index and cross section options
      j = j+1
      jlabel(j) = 'CCl4 -> Products'

* cross section options
      spc = "CCl4"
      xsvers = (/2, 2, 1/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg( 1) = "from JPL 1997 recommendation"
      logmsg( 2) = "from JPL 2011 with T dependence"
      CALL set_option(2,spc,"abs",logmsg,xsvers,mabs)


*** cross sections from JPL97 recommendation (identical to 94 data)

      IF(mabs .EQ. 1) THEN

         OPEN(kin,FILE='DATAJ1/ABS/CCl4_jpl94.abs',STATUS='OLD')
         READ(kin,*) idum, n
         DO i = 1, idum-2
            READ(kin,*)
         ENDDO
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
            y1(i) = y1(i) * 1E-20
         ENDDO
         CLOSE(kin)

         CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)

      ELSEIF(mabs .EQ. 2) THEN

         OPEN(kin,FILE='DATAJ1/ABS/CCl4_jpl11.abs',STATUS='OLD')
         DO i = 1, 5
            READ(kin,*)
         ENDDO
         n = 44
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
            y1(i) = y1(i) * 1E-20
         ENDDO
         CLOSE(kin)

         CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)

      ENDIF

* compute temperature correction factors:

      b0 = 1.0739
      b1 = -1.6275e-2
      b2 = 8.8141e-5
      b3 = -1.9811e-7
      b4 = 1.5022e-10

*** quantum yield assumed to be unity

      qy = 1.
      DO iw = 1, nw-1

* compute temperature correction coefficients:

         tcoeff = 0.
         IF(wc(iw) .GT. 194. .AND. wc(iw) .LT. 250.) THEN
            w1 = wc(iw)
            w2 = w1**2
            w3 = w1**3
            w4 = w1**4
            tcoeff = b0 + b1*w1 + b2*w2 + b3*w3 + b4*w4
         ENDIF

         DO iz = 1, nz

            IF(mabs .EQ. 1) THEN
               sig = yg(iw)
            ELSEIF (mabs .EQ. 2) THEN

               temp = tlev(iz)
               temp = min(max(temp,210.),300.)

               sig = yg(iw) * 10.**(tcoeff*(temp-295.))
            ENDIF

            sq(j,iz,iw) = qy * sig

         ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r21(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! CClFO

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for CClFO photolysis:  =*
*=         CClFO + hv -> Products                                            =*
*=  Cross section: from JPL 97                                               =*
*=  Quantum yield: assumed to be unity                                       =*
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
      PARAMETER(kdata=100)

      REAL x1(kdata)
      REAL y1(kdata)

* local

      REAL yg(kw)
C      REAL qy
      INTEGER i, iw, n, idum
      INTEGER iz


***************** CClFO photolysis *************************

      j = j+1
      jlabel(j) = 'CClFO -> Products'

*** cross sections from JPL97 recommendation (identical to 94 recommendation)

      OPEN(kin,FILE='DATAJ1/ABS/CClFO_jpl94.abs',STATUS='OLD')
      READ(kin,*) idum, n
      DO i = 1, idum-2
        READ(kin,*)
      ENDDO
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
        y1(i) = y1(i) * 1E-20
      ENDDO
      CLOSE(kin)

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)

*** quantum yield unity
C      qy = 1.
      DO iw = 1, nw-1
        DO iz = 1, nz
          sq(j,iz,iw) = yg(iw)! * qy
        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r22(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! CF2O

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for CF2O photolysis:   =*
*=        CF2O + hv -> Products                                              =*
*=  Cross section:  from JPL 97 recommendation                               =*
*=  Quantum yield:  unity (Nolle et al.)                                     =*
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
      PARAMETER(kdata=100)

      REAL x1(kdata)
      REAL y1(kdata)

* local

      REAL yg(kw)
C      REAL qy
      INTEGER i, iw, n
      INTEGER iz


***************** CF2O photolysis *************************

      j = j+1
      jlabel(j) = 'CF2O -> Products'

**** cross sections from JPL97 recommendation (identical to 94 recommendation)

      OPEN(kin,FILE='DATAJ1/ABS/CF2O_jpl11.abs',STATUS='OLD')
      DO i = 1, 5
        READ(kin,*)
      ENDDO
      n = 21
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
        y1(i) = y1(i) * 1E-20
      ENDDO
      CLOSE(kin)

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)

*** quantum yield unity (Nolle et al.)
C      qy = 1.
      DO iw = 1, nw-1
        DO iz = 1, nz
           sq(j,iz,iw) = yg(iw)! * qy
        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r23(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! CF2ClCFCl2

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for CFC-113 photolysis:=*
*=          CF2ClCFCl2 + hv -> Products                                      =*
*=  Cross section:  from JPL 97 recommendation, linear interp. between       =*
*=                  values at 210 and 295K                                   =*
*=  Quantum yield:  assumed to be unity                                      =*
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
      PARAMETER(kdata=100)

      INTEGER n1, n2
      REAL x1(kdata), x2(kdata)
      REAL y1(kdata), y2(kdata)

* local

      REAL yg1(kw), yg2(kw)
C      REAL qy
      REAL t
      INTEGER i, iw, n, idum
      INTEGER iz
      REAL slope


***************** CF2ClCFCl2 (CFC-113) photolysis *************************

      j = j+1
      jlabel(j) = 'CF2ClCFCl2 (CFC-113) -> Products'

*** cross sections from JPL97 recommendation (identical to 94 recommendation)

      OPEN(kin,FILE='DATAJ1/ABS/CFC-113_jpl94.abs',STATUS='OLD')
      READ(kin,*) idum, n
      DO i = 1, idum-2
        READ(kin,*)
      ENDDO
      DO i = 1, n
        READ(kin,*) x1(i), y1(i), y2(i)
        y1(i) = y1(i) * 1E-20
        y2(i) = y2(i) * 1E-20
        x2(i) = x1(i)
      ENDDO
      CLOSE(kin)

      n1 = n
      n2 = n

** sigma @ 295 K
      CALL interpol(x1,y1,kdata,n1,nw,wl,jlabel(j),0.,0.,yg1)

* sigma @ 210 K
      CALL interpol(x2,y2,kdata,n2,nw,wl,jlabel(j),0.,0.,yg2)

*** quantum yield assumed to be unity
C      qy = 1.

      DO iz = 1, nz
        t = MAX(210.,MIN(tlev(iz),295.))
        slope = (t-210.)/(295.-210.)
        DO iw = 1, nw-1
            sq(j,iz,iw) = (yg2(iw) + slope*(yg1(iw)-yg2(iw)))! * qy
        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r24(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! CF2ClCF2Cl

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for CFC-144 photolysis:=*
*=              CF2ClCF2Cl + hv -> Products                                  =*
*=  Cross section: from JPL 97 recommendation, linear interp. between values =*
*=                 at 210 and 295K                                           =*
*=  Quantum yield: assumed to be unity                                       =*
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
      PARAMETER(kdata=100)

      INTEGER n1, n2
      REAL x1(kdata), x2(kdata)
      REAL y1(kdata), y2(kdata)

* local

      REAL yg1(kw), yg2(kw)
C      REAL qy
      REAL t
      INTEGER i, iw, n, idum
      INTEGER iz
      REAL slope


***************** CF2ClCF2Cl (CFC-114) photolysis *************************

      j = j+1
      jlabel(j) = 'CF2ClCF2Cl (CFC-114) -> Products'

**** cross sections from JPL97 recommendation (identical to 94 recommendation)

      OPEN(kin,FILE='DATAJ1/ABS/CFC-114_jpl94.abs',STATUS='OLD')
      READ(kin,*) idum, n
      DO i = 1, idum-2
        READ(kin,*)
      ENDDO
      DO i = 1, n
        READ(kin,*) x1(i), y1(i), y2(i)
        y1(i) = y1(i) * 1E-20
        y2(i) = y2(i) * 1E-20
        x2(i) = x1(i)
      ENDDO
      CLOSE(kin)

      n1 = n
      n2 = n

** sigma @ 295 K
      CALL interpol(x1,y1,kdata,n1,nw,wl,jlabel(j),0.,0.,yg1)

* sigma @ 210 K
      CALL interpol(x2,y2,kdata,n2,nw,wl,jlabel(j),0.,0.,yg2)

*** quantum yield assumed to be unity
C      qy = 1.

      DO iz = 1, nz
        t = MAX(210.,MIN(tlev(iz),295.))
        slope = (t-210.)/(295.-210.)
        DO iw = 1, nw-1
            sq(j,iz,iw) = (yg2(iw) + slope*(yg1(iw)-yg2(iw)))! * qy
        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r25(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! CF3CF2Cl

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for CFC-115 photolysis =*
*=             CF3CF2Cl + hv -> Products                                     =*
*=  Cross section: from JPL 97 recommendation                                =*
*=  Quantum yield: assumed to be unity                                       =*
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
      PARAMETER(kdata=100)

      REAL x1(kdata)
      REAL y1(kdata)

* local

      REAL yg(kw)
C      REAL qy
      INTEGER i, iw, n, idum
      INTEGER iz


***************** CF3CF2Cl (CFC-115) photolysis *************************

      j = j+1
      jlabel(j) = 'CF3CF2Cl (CFC-115) -> Products'

**** cross sections from JPL97 recommendation (identical to 94 recommendation)

      OPEN(kin,FILE='DATAJ1/ABS/CFC-115_jpl94.abs',STATUS='OLD')
      READ(kin,*) idum, n
      DO i = 1, idum-2
        READ(kin,*)
      ENDDO
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
        y1(i) = y1(i) * 1E-20
      ENDDO
      CLOSE(kin)

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)

**** quantum yield assumed to be unity
C      qy = 1.

      DO iw = 1, nw-1
        DO iz = 1, nz
          sq(j,iz,iw) = yg(iw)! * qy
        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r26(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! CCl3F

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for CFC-11  photolysis =*
*=          CCl3F + hv -> Products                                           =*
*=  Cross section: from JPL 97 recommendation                                =*
*=  Quantum yield: assumed to be unity                                       =*
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
      PARAMETER(kdata=100)

      REAL x1(kdata)
      REAL y1(kdata)

* local

      REAL yg(kw)
C      REAL qy
      REAL t
      INTEGER i, iw, n, idum
      INTEGER iz


***************** Cl3F (CFC-11) photolysis *************************

      j = j+1
      jlabel(j) = 'CCl3F (CFC-11) -> Products'

**** cross sections from JPL97 recommendation (identical to 94 recommendation)

      OPEN(kin,FILE='DATAJ1/ABS/CFC-11_jpl94.abs',STATUS='OLD')
      READ(kin,*) idum, n
      DO i = 1, idum-2
        READ(kin,*)
      ENDDO
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
        y1(i) = y1(i) * 1E-20
      ENDDO
      CLOSE(kin)

** sigma @ 298 K
      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)

**** quantum yield assumed to be unity
C      qy = 1.

      DO iz = 1, nz
        t = 1E-04 * (tlev(iz)-298.)
        DO iw = 1, nw-1
          sq(j,iz,iw) = yg(iw) * EXP((wc(iw)-184.9) * t)! * qy
        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r27(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! CCl2F2

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for CFC-12  photolysis:=*
*=         CCl2F2 + hv -> Products                                           =*
*=  Cross section: from JPL 97 recommendation                                =*
*=  Quantum yield: assumed to be unity                                       =*
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
      PARAMETER(kdata=100)

      REAL x1(kdata)
      REAL y1(kdata)

* local

      REAL yg(kw)
C      REAL qy
      REAL t
      INTEGER i, iw, n, idum
      INTEGER iz


***************** CCl2F2 (CFC-12)  photolysis *************************

      j = j+1
      jlabel(j) = 'CCl2F2 (CFC-12) -> Products'

**** cross sections from JPL97 recommendation (identical to 94 recommendation)

      OPEN(kin,FILE='DATAJ1/ABS/CFC-12_jpl94.abs',STATUS='OLD')
      READ(kin,*) idum, n
      DO i = 1, idum-2
        READ(kin,*)
      ENDDO
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
        y1(i) = y1(i) * 1E-20
      ENDDO
      CLOSE(kin)

** sigma @ 298 K
      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)

**** quantum yield assumed to be unity
C      qy = 1.

      DO iz = 1, nz
        t = 1E-04 * (tlev(iz)-298.)
        DO iw = 1, nw-1
          sq(j,iz,iw) = yg(iw) * EXP((wc(iw)-184.9) * t)! * qy
        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r28(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! CH3Br

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for CH3Br photolysis:  =*
*=         CH3Br + hv -> Products                                            =*
*=  Cross section: from JPL 97 recommendation                                =*
*=  Quantum yield: assumed to be unity                                       =*
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
      PARAMETER(kdata=100)

      REAL x1(kdata)
      REAL y1(kdata)

* local

      REAL yg(kw)
C      REAL qy
      INTEGER i, iw, n, idum
      INTEGER iz


***************** CH3Br photolysis *************************

* data from JPL97 (identical to 94 recommendation)

      j = j+1
      jlabel(j) = 'CH3Br -> Products'

      OPEN(kin,FILE='DATAJ1/ABS/CH3Br_jpl94.abs',STATUS='OLD')
      READ(kin,*) idum, n
      DO i = 1, idum-2
        READ(kin,*)
      ENDDO
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
        y1(i) = y1(i) * 1E-20
      ENDDO
      CLOSE(kin)
      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)

**** quantum yield assumed to be unity
C      qy = 1.

      DO iw = 1, nw-1
        DO iz = 1, nz
          sq(j,iz,iw) = yg(iw)! * qy
        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r29(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! CH3CCl3

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for CH3CCl3 photolysis =*
*=           CH3CCl3 + hv -> Products                                        =*
*=  Cross section: from JPL 97 recommendation, piecewise linear interp.      =*
*=                 of data at 210, 250, and 295K                             =*
*=  Quantum yield: assumed to be unity                                       =*
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
      PARAMETER(kdata=100)

      INTEGER n1, n2, n3
      REAL x1(kdata), x2(kdata), x3(kdata)
      REAL y1(kdata), y2(kdata), y3(kdata)

* local

      REAL yg1(kw), yg2(kw), yg3(kw)
C      REAL qy
      REAL t
      INTEGER i, iw, n, idum
      INTEGER iz
      REAL slope


***************** CH3CCl3 photolysis *************************

      j = j+1
      jlabel(j) = 'CH3CCl3 -> Products'

**** cross sections from JPL97 recommendation (identical to 94 recommendation)

      OPEN(kin,FILE='DATAJ1/ABS/CH3CCl3_jpl94.abs',STATUS='OLD')
      READ(kin,*) idum, n
      DO i = 1, idum-2
        READ(kin,*)
      ENDDO
      DO i = 1, n
        READ(kin,*) x1(i), y1(i), y2(i), y3(i)
        y1(i) = y1(i) * 1E-20
        y2(i) = y2(i) * 1E-20
        y3(i) = y3(i) * 1E-20
        x2(i) = x1(i)
        x3(i) = x1(i)
      ENDDO
      CLOSE(kin)

      n1 = n
      n2 = n
      n3 = n

** sigma @ 295 K
      CALL interpol(x1,y1,kdata,n1,nw,wl,jlabel(j),0.,0.,yg1)

** sigma @ 250 K
      CALL interpol(x2,y2,kdata,n2,nw,wl,jlabel(j),0.,0.,yg2)

** sigma @ 210 K
      CALL interpol(x3,y3,kdata,n3,nw,wl,jlabel(j),0.,0.,yg3)

**** quantum yield assumed to be unity
C      qy = 1.

      DO iz = 1, nz
        t = MIN(295.,MAX(tlev(iz),210.))
        IF (t .LE. 250.) THEN
          slope = (t-210.)/(250.-210.)
          DO iw = 1, nw-1
            sq(j,iz,iw) = (yg3(iw) + slope*(yg2(iw)-yg3(iw)))! * qy
          ENDDO
        ELSE
          slope = (t-250.)/(295.-250.)
          DO iw = 1, nw-1
            sq(j,iz,iw) = (yg2(iw) + slope*(yg1(iw)-yg2(iw)))! * qy
          ENDDO
        ENDIF
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r30(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! CH3Cl

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for CH3Cl photolysis:  =*
*=            CH3Cl + hv -> Products                                         =*
*=  Cross section: from JPL 97 recommendation, piecewise linear interp.      =*
*=                 from values at 255, 279, and 296K                         =*
*=  Quantum yield: assumed to be unity                                       =*
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
      PARAMETER(kdata=100)

      INTEGER n1, n2, n3
      REAL x1(kdata), x2(kdata), x3(kdata)
      REAL y1(kdata), y2(kdata), y3(kdata)

* local

      REAL yg1(kw), yg2(kw), yg3(kw)
C      REAL qy
      REAL t
      INTEGER i, iw, n, idum
      INTEGER iz
      REAL slope


***************** CH3Cl photolysis *************************

      j = j+1
      jlabel(j) = 'CH3Cl -> Products'

**** cross sections from JPL97 recommendation (identical to 94 recommendation)

      OPEN(kin,FILE='DATAJ1/ABS/CH3Cl_jpl94.abs',STATUS='OLD')
      READ(kin,*) idum, n
      DO i = 1, idum-2
        READ(kin,*)
      ENDDO
      DO i = 1, n
        READ(kin,*) x1(i), y1(i), y2(i), y3(i)
        y1(i) = y1(i) * 1E-20
        y2(i) = y2(i) * 1E-20
        y3(i) = y3(i) * 1E-20
        x2(i) = x1(i)
        x3(i) = x1(i)
      ENDDO
      CLOSE(kin)

      n1 = n
      n2 = n
      n3 = n

** sigma @ 296 K
      CALL interpol(x1,y1,kdata,n1,nw,wl,jlabel(j),0.,0.,yg1)

** sigma @ 279 K
      CALL interpol(x2,y2,kdata,n2,nw,wl,jlabel(j),0.,0.,yg2)

** sigma @ 255 K
      CALL interpol(x3,y3,kdata,n3,nw,wl,jlabel(j),0.,0.,yg3)

**** quantum yield assumed to be unity
C      qy = 1.

      DO iz = 1, nz
        t = MAX(255.,MIN(tlev(i),296.))
        IF (t .LE. 279.) THEN
          slope = (t-255.)/(279.-255.)
          DO iw = 1, nw-1
            sq(j,iz,iw) = (yg3(iw)+slope*(yg2(iw)-yg3(iw)))! * qy
          ENDDO
        ELSE
          slope = (t-279.)/(296.-279.)
          DO iw = 1, nw-1
            sq(j,iz,iw) = (yg2(iw)+slope*(yg1(iw)-yg2(iw)))! * qy
          ENDDO
        ENDIF
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r31(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! ClOO

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for ClOO photolysis:   =*
*=          ClOO + hv -> Products                                            =*
*=  Cross section: from JPL 97 recommendation                                =*
*=  Quantum yield: assumed to be unity                                       =*
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
      PARAMETER(kdata=100)

C     INTEGER n1, n2, n3, n4, n5
      REAL x1(kdata)
      REAL y1(kdata)

* local

      REAL yg(kw)
C      REAL qy
      INTEGER i, iw, n, idum
      INTEGER iz


***************** ClOO photolysis *************************

      j = j+1
      jlabel(j) = 'ClOO -> Products'

**** cross sections from JPL97 recommendation (identical to 94 recommendation)
** also identical to JPL2011 recommendation

      OPEN(kin,FILE='DATAJ1/ABS/ClOO_jpl94.abs',STATUS='OLD')
      READ(kin,*) idum, n
      DO i = 1, idum-2
        READ(kin,*)
      ENDDO
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
        y1(i) = y1(i) * 1E-20
      ENDDO
      CLOSE(kin)
      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)

**** quantum yield assumed to be unity
C      qy = 1.

      DO iw = 1, nw-1
        DO iz = 1, nz
          sq(j,iz,iw) = yg(iw)! * qy
        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r32(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! CF3CHCl2

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for HCFC-123 photolysis=*
*=       CF3CHCl2 + hv -> Products                                           =*
*=  Cross section: from Orlando et al., 1991                                 =*
*=  Quantum yield: assumed to be unity                                       =*
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

* local

C      REAL qy
      REAL t
      INTEGER i, iw, idum
      INTEGER iz, k
      REAL lambda, sum
      CHARACTER*120 inline

      REAL coeff(4,3), TBar, LBar


***************** CF3CHCl2 (HCFC-123) photolysis *************************

      j = j+1
      jlabel(j) = 'CF3CHCl2 (HCFC-123) -> Products'

**** cross sections from JPL94 recommendation

C     OPEN(kin,FILE='DATAJ1/ABS/HCFC-123_jpl94.abs',STATUS='OLD')
C     READ(kin,*) idum, n
C     DO i = 1, idum-2
C       READ(kin,*)
C     ENDDO
C     DO i = 1, n
C       READ(kin,*) x1(i), y1(i)
C       y1(i) = y1(i) * 1E-20
C     ENDDO
C     CLOSE(kin)

C     CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
C     CALL addpnt(x1,y1,kdata,n,          0.,0.)
C     CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
C     CALL addpnt(x1,y1,kdata,n,        1E38,0.)

C     CALL inter2(nw,wl,yg,n,x1,y1,ierr)

C     IF (ierr .NE. 0) THEN
C        WRITE(*,*) ierr, jlabel(j)
C        STOP
C     ENDIF

**** quantum yield assumed to be unity
C     qy = 1.

C     DO iw = 1, nw-1
C       DO iz = 1, nz
C         sq(j,iz,iw) = qy * yg(iw)
C       ENDDO
C     ENDDO


**** cross section from Orlando et al., 1991

      OPEN(kin,FILE='DATAJ1/ABS/HCFCs_orl.abs',STATUS='OLD')
      READ(kin,*) idum
      DO i = 1, idum-2
        READ(kin,*)
      ENDDO
      READ(kin,100) inline
 100  FORMAT(A120)
      READ(inline(6:),*) TBar,i,(coeff(i,k),k=1,3)
      READ(kin,*)           i,(coeff(i,k),k=1,3)
      READ(kin,*)           i,(coeff(i,k),k=1,3)
      READ(kin,*)           i,(coeff(i,k),k=1,3)
      CLOSE(kin)

      LBar = 206.214

**** quantum yield assumed to be unity

C      qy = 1.

      DO iw = 1, nw-1

        lambda = wc(iw)

C use parameterization only up to 220 nm, as the error bars associated with
C the measurements beyond 220 nm are very large (Orlando, priv.comm.)

        IF (lambda .GE. 190. .AND. lambda .LE. 220.) THEN
          DO iz = 1, nz
             t = MIN(295.,MAX(tlev(i),203.))-TBar
             sum = 0.
             DO i = 1, 4
                sum = (coeff(i,1)+t*(coeff(i,2)+t*coeff(i,3))) *
     >                (lambda-LBar)**(i-1) + sum
             ENDDO
             sq(j,iz,iw) = EXP(sum)! * qy
          ENDDO
        ELSE
          DO iz = 1, nz
            sq(j,iz,iw) = 0.
          ENDDO
        ENDIF
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r33(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! CF3CHFCl

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for HCFC-124 photolysis=*
*=        CF3CHFCl + hv -> Products                                          =*
*=  Cross section: from Orlando et al., 1991                                 =*
*=  Quantum yield: assumed to be unity                                       =*
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


* local

C      REAL qy
      REAL t
      INTEGER i, iw, idum
      INTEGER iz, k
      REAL lambda, sum
      CHARACTER*120 inline

      REAL coeff(4,3), TBar, LBar


***************** CF3CHFCl (HCFC-124) photolysis *************************

      j = j+1
      jlabel(j) = 'CF3CHFCl (HCFC-124) -> Products'

**** cross sections from JPL94 recommendation

C     OPEN(kin,FILE='DATAJ1/ABS/HCFC-124_jpl94.abs',STATUS='OLD')
C     READ(kin,*) idum, n
C     DO i = 1, idum-2
C       READ(kin,*)
C     ENDDO
C     DO i = 1, n
C       READ(kin,*) x1(i), y1(i)
C       y1(i) = y1(i) * 1E-20
C     ENDDO
C     CLOSE(kin)

C     CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
C     CALL addpnt(x1,y1,kdata,n,          0.,0.)
C     CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
C     CALL addpnt(x1,y1,kdata,n,        1E38,0.)

C     CALL inter2(nw,wl,yg,n,x1,y1,ierr)

C     IF (ierr .NE. 0) THEN
C       WRITE(*,*) ierr, jlabel(j)
C       STOP
C     ENDIF

**** quantum yield assumed to be unity
C     qy = 1.

C     DO iw = 1, nw-1
C       DO iz = 1, nz
C         sq(j,iz,iw) = qy * yg(iw)
C       ENDDO
C     ENDDO

**** cross section from Orlando et al., 1991

      OPEN(kin,FILE='DATAJ1/ABS/HCFCs_orl.abs',STATUS='OLD')
      READ(kin,*) idum
      idum = idum+5
      DO i = 1, idum-2
        READ(kin,*)
      ENDDO
      READ(kin,100) inline
 100  FORMAT(A120)
      READ(inline(6:),*) TBar,i,(coeff(i,k),k=1,3)
      READ(kin,*)           i,(coeff(i,k),k=1,3)
      READ(kin,*)           i,(coeff(i,k),k=1,3)
      READ(kin,*)           i,(coeff(i,k),k=1,3)
      CLOSE(kin)

      LBar = 206.214

**** quantum yield assumed to be unity

C      qy = 1.

      DO iw = 1, nw-1
        lambda = wc(iw)
        IF (lambda .GE. 190. .AND. lambda .LE. 230.) THEN
          DO iz = 1, nz
             t = MIN(295.,MAX(tlev(i),203.))-TBar
             sum = 0.
             DO i = 1, 4
                sum = (coeff(i,1)+t*(coeff(i,2)+t*coeff(i,3))) *
     >                (lambda-LBar)**(i-1) + sum
             ENDDO
             sq(j,iz,iw) = EXP(sum)! * qy
          ENDDO
        ELSE
          DO iz = 1, nz
            sq(j,iz,iw) = 0.
          ENDDO
        ENDIF
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r34(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! CH3CFCl2

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for HCFC-141b          =*
*=  photolysis:                                                              =*
*=         CH3CFCl2 + hv -> Products                                         =*
*=  Cross section: from JPL97 recommendation                                 =*
*=  Quantum yield: assumed to be unity                                       =*
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
      PARAMETER(kdata=100)

      REAL x1(kdata)
      REAL y1(kdata)

* local

      REAL yg(kw)
C      REAL qy
      INTEGER i, iw, n, idum
      INTEGER iz


***************** CH3CFCl2 (HCFC-141b) photolysis *************************

      j = j+1
      jlabel(j) = 'CH3CFCl2 (HCFC-141b) -> Products'

**** cross sections from JPL97 recommendation (identical to 94 recommendation)

      OPEN(kin,FILE='DATAJ1/ABS/HCFC-141b_jpl94.abs',STATUS='OLD')
      READ(kin,*) idum, n
      DO i = 1, idum-2
        READ(kin,*)
      ENDDO
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
        y1(i) = y1(i) * 1E-20
      ENDDO
      CLOSE(kin)

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)

**** quantum yield assumed to be unity
C      qy = 1.

      DO iw = 1, nw-1
        DO iz = 1, nz
          sq(j,iz,iw) = yg(iw)! * qy
        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r35(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! HCFC-142b

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for HCFC-142b          =*
*=  photolysis:                                                              =*
*=          CH3CF2Cl + hv -> Products                                        =*
*=  Cross section: from Orlando et al., 1991                                 =*
*=  Quantum yield: assumed to be unity                                       =*
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

* local

C      REAL qy
      REAL t
      INTEGER i, iw, idum
      INTEGER iz, k
      REAL lambda, sum
      CHARACTER*80 inline

      REAL coeff(4,3), TBar, LBar


***************** CH3CF2Cl (HCFC-142b) photolysis *************************

      j = j+1
      jlabel(j) = 'CH3CF2Cl (HCFC-142b) -> Products'

**** cross sections from JPL94 recommendation

C     OPEN(kin,FILE='DATAJ1/ABS/HCFC-142b_jpl94.abs',STATUS='OLD')
C     READ(kin,*) idum, n
C     DO i = 1, idum-2
C       READ(kin,*)
C     ENDDO
C     DO i = 1, n
C       READ(kin,*) x1(i), y1(i)
C       y1(i) = y1(i) * 1E-20
C     ENDDO
C     CLOSE(kin)

C     CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
C     CALL addpnt(x1,y1,kdata,n,          0.,0.)
C     CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
C     CALL addpnt(x1,y1,kdata,n,        1E38,0.)

C     CALL inter2(nw,wl,yg,n,x1,y1,ierr)

C     IF (ierr .NE. 0) THEN
C       WRITE(*,*) ierr, jlabel(j)
C       STOP
C     ENDIF

**** quantum yield assumed to be unity
C     qy = 1.

C     DO iw = 1, nw-1
C       DO iz = 1, nz
C         sq(j,iz,iw) = qy * yg(iw)
C       ENDDO
C     ENDDO

**** cross section from Orlando et al., 1991

      OPEN(kin,FILE='DATAJ1/ABS/HCFCs_orl.abs',STATUS='OLD')
      READ(kin,*) idum
      idum = idum+10
      DO i = 1, idum-2
        READ(kin,*)
      ENDDO
      READ(kin,100) inline
 100  FORMAT(A80)
      READ(inline(6:),*) TBar,i,(coeff(i,k),k=1,3)
      READ(kin,*)           i,(coeff(i,k),k=1,3)
      READ(kin,*)           i,(coeff(i,k),k=1,3)
      READ(kin,*)           i,(coeff(i,k),k=1,3)
      CLOSE(kin)

      LBar = 206.214

**** quantum yield assumed to be unity

C      qy = 1.

      DO iw = 1, nw-1

        lambda = wc(iw)
        IF (lambda .GE. 190. .AND. lambda .LE. 230.) THEN

          DO iz = 1, nz
             t = MIN(295.,MAX(tlev(i),203.))-TBar
             sum = 0.
             DO i = 1, 4
                sum = (coeff(i,1)+t*(coeff(i,2)+t*coeff(i,3))) *
     >                (lambda-LBar)**(i-1) + sum
             ENDDO

* offeset exponent by 40 (exp(-40.) = 4.248e-18) to prevent exp. underflow errors
* on some machines.

c             sq(j,iz,iw) = qy * EXP(sum)
             sq(j,iz,iw) = 4.248e-18 * EXP(sum + 40.)! * qy

          ENDDO

        ELSE
          DO iz = 1, nz
            sq(j,iz,iw) = 0.
          ENDDO
        ENDIF

      ENDDO


      END

*=============================================================================*

      SUBROUTINE r36(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! CF3CF2CHCl2

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for HCFC-225ca         =*
*=  photolysis:                                                              =*
*=           CF3CF2CHCl2 + hv -> Products                                    =*
*=  Cross section: from JPL 97 recommendation                                =*
*=  Quantum yield: assumed to be unity                                       =*
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
      PARAMETER(kdata=100)

      REAL x1(kdata)
      REAL y1(kdata)

* local

      REAL yg(kw)
C      REAL qy
      INTEGER i, iw, n, idum
      INTEGER iz


***************** CF3CF2CHCl2 (HCFC-225ca) photolysis *************************

      j = j+1
      jlabel(j) = 'CF3CF2CHCl2 (HCFC-225ca) -> Products'

**** cross sections from JPL97 recommendation (identical to 94 recommendation)

      OPEN(kin,FILE='DATAJ1/ABS/HCFC-225ca_jpl94.abs',STATUS='OLD')
      READ(kin,*) idum, n
      DO i = 1, idum-2
        READ(kin,*)
      ENDDO
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
        y1(i) = y1(i) * 1E-20
      ENDDO
      CLOSE(kin)

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)

**** quantum yield assumed to be unity
C      qy = 1.

      DO iw = 1, nw-1
        DO iz = 1, nz
          sq(j,iz,iw) = yg(iw)! * qy
        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r37(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! CF2ClCF2CHFCl

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for HCFC-225cb         =*
*=  photolysis:                                                              =*
*=          CF2ClCF2CHFCl + hv -> Products                                   =*
*=  Cross section: from JPL 97 recommendation                                =*
*=  Quantum yield: assumed to be unity                                       =*
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
      PARAMETER(kdata=100)

      REAL x1(kdata)
      REAL y1(kdata)

* local

      REAL yg(kw)
C      REAL qy
      INTEGER i, iw, n, idum
      INTEGER iz


***************** CF2ClCF2CHFCl (HCFC-225cb) photolysis *************************

      j = j+1
      jlabel(j) = 'CF2ClCF2CHFCl (HCFC-225cb) -> Products'

**** cross sections from JPL97 recommendation (identical to 94 recommendation)

      OPEN(kin,FILE='DATAJ1/ABS/HCFC-225cb_jpl94.abs',STATUS='OLD')
      READ(kin,*) idum, n
      DO i = 1, idum-2
        READ(kin,*)
      ENDDO
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
        y1(i) = y1(i) * 1E-20
      ENDDO
      CLOSE(kin)

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)

**** quantum yield assumed to be unity
C      qy = 1.

      DO iw = 1, nw-1
        DO iz = 1, nz
          sq(j,iz,iw) = yg(iw)! * qy
        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r38(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! CHClF2

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for HCFC-22 photolysis =*
*=          CHClF2 + hv -> Products                                          =*
*=  Cross section: from JPL 97 recommendation, piecewise linear interp.      =*
*=                 from values at 210, 230, 250, 279, and 295 K              =*
*=  Quantum yield: assumed to be unity                                       =*
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
      PARAMETER(kdata=100)

      INTEGER n1, n2, n3, n4, n5
      REAL x1(kdata), x2(kdata), x3(kdata), x4(kdata), x5(kdata)
      REAL y1(kdata), y2(kdata), y3(kdata), y4(kdata), y5(kdata)

* local

      REAL yg1(kw), yg2(kw), yg3(kw), yg4(kw), yg5(kw)
C      REAL qy
      REAL t
      INTEGER i, iw, n, idum
      INTEGER iz
      REAL slope


***************** CHClF2 (HCFC-22) photolysis *************************

      j = j+1
      jlabel(j) = 'CHClF2 (HCFC-22) -> Products'

**** cross sections from JPL97 recommendation (identical to 94 recommendation)

      OPEN(kin,FILE='DATAJ1/ABS/HCFC-22_jpl94.abs',STATUS='OLD')
      READ(kin,*) idum, n
      DO i = 1, idum-2
        READ(kin,*)
      ENDDO
      DO i = 1, n
        READ(kin,*) x1(i), y1(i), y2(i), y3(i), y4(i), y5(i)
        y1(i) = y1(i) * 1E-20
        y2(i) = y2(i) * 1E-20
        y3(i) = y3(i) * 1E-20
        y4(i) = y4(i) * 1E-20
        y5(i) = y5(i) * 1E-20
        x2(i) = x1(i)
        x3(i) = x1(i)
        x4(i) = x1(i)
        x5(i) = x1(i)
      ENDDO
      CLOSE(kin)

      n1 = n
      n2 = n
      n3 = n
      n4 = n
      n5 = n

** sigma @ 295 K
      CALL interpol(x1,y1,kdata,n1,nw,wl,jlabel(j),0.,0.,yg1)

** sigma @ 270 K
      CALL interpol(x2,y2,kdata,n2,nw,wl,jlabel(j),0.,0.,yg2)

** sigma @ 250 K
      CALL interpol(x3,y3,kdata,n3,nw,wl,jlabel(j),0.,0.,yg3)

** sigma @ 230 K
      CALL interpol(x4,y4,kdata,n4,nw,wl,jlabel(j),0.,0.,yg4)

** sigma @ 210 K
      CALL interpol(x5,y5,kdata,n5,nw,wl,jlabel(j),0.,0.,yg5)

**** quantum yield assumed to be unity
C      qy = 1.

      DO iz = 1, nz
         t = MIN(295.,MAX(tlev(iz),210.))
         IF (t .LE. 230.) THEN
            slope = (t-210.)/(230.-210.)
            DO iw = 1, nw-1
              sq(j,iz,iw) = (yg5(iw)+slope*(yg4(iw)-yg5(iw)))! * qy
            ENDDO
         ELSEIF (t .LE. 250.) THEN
            slope = (t-230.)/(250.-230.)
            DO iw = 1, nw-1
              sq(j,iz,iw) = (yg4(iw)+slope*(yg3(iw)-yg4(iw)))! * qy
            ENDDO
         ELSEIF (t .LE. 270.) THEN
            slope = (t-250.)/(270.-250.)
            DO iw = 1, nw-1
              sq(j,iz,iw) = (yg3(iw)+slope*(yg2(iw)-yg3(iw)))! * qy
            ENDDO
         ELSE
            slope = (t-270.)/(295.-270.)
            DO iw = 1, nw-1
              sq(j,iz,iw) = (yg2(iw)+slope*(yg1(iw)-yg2(iw)))! * qy
            ENDDO
         ENDIF
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r39(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! HO2

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for HO2 photolysis:    =*
*=          HO2 + hv -> OH + O                                               =*
*=  Cross section: from JPL 97 recommendation                                =*
*=  Quantum yield: assumed shape based on work by Lee, 1982; normalized      =*
*=                 to unity at 248 nm                                        =*
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
      PARAMETER(kdata=100)

      REAL x1(kdata)
      REAL y1(kdata)

* local

      REAL yg(kw)
      REAL qy
      INTEGER i, iw, n
      INTEGER iz


***************** HO2 photolysis *************************

      j = j+1
      jlabel(j) = 'HO2 -> OH + O'

**** cross sections from JPL11 recommendation

      OPEN(kin,FILE='DATAJ1/ABS/HO2_jpl11.abs',STATUS='OLD')
      DO i = 1, 10
        READ(kin,*)
      ENDDO
      n = 15
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
        y1(i) = y1(i) * 1E-20
      ENDDO
      CLOSE(kin)
      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)

**** quantum yield:  absolute quantum yield has not been reported yet, but
****                 Lee measured a quantum yield for O(1D) production at 248
****                 nm that was 15 time larger than at 193 nm
**** here:  a quantum yield of unity is assumed at 248 nm and beyond, for
****        shorter wavelengths a linear decrease with lambda is assumed

      DO iw = 1, nw-1
         IF (wc(iw) .GE. 248.) THEN
            qy = 1.
         ELSE
            qy = 1./15. + (wc(iw)-193.)*(14./15.)/(248.-193.)
            qy = MAX(qy,0.)
         ENDIF
         DO iz = 1, nz
           sq(j,iz,iw) = qy * yg(iw)
         ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r40(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! CF2Br2

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) Halon-1202 photolysis: =*
*=         CF2Br2 + hv -> Products                                           =*
*=  Cross section: from JPL 97 recommendation                                =*
*=  Quantum yield: unity (Molina and Molina)                                 =*
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
      PARAMETER(kdata=100)

      REAL x1(kdata)
      REAL y1(kdata)

* local

      REAL yg(kw)
C      REAL qy
      INTEGER i, iw, n, idum
      INTEGER iz


***************** CF2Br2 (Halon-1202) photolysis *************************

      j = j+1
      jlabel(j) = 'CF2Br2 (Halon-1202) -> Products'

**** cross sections from JPL97 recommendation (identical to 94 recommendation)

      OPEN(kin,FILE='DATAJ1/ABS/Halon-1202_jpl97.abs',STATUS='OLD')
      READ(kin,*) idum, n
      DO i = 1, idum-2
        READ(kin,*)
      ENDDO
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
        y1(i) = y1(i) * 1E-20
      ENDDO
      CLOSE(kin)

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)

**** quantum yield unity (Molina and Molina)
C      qy = 1.

      DO iw = 1, nw-1
        DO iz = 1, nz
           sq(j,iz,iw) = yg(iw)! * qy
        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r41(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! CF2ClBr

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for Halon-1211         =*
*=  photolysis:                                                              =*
*=           CF2ClBr + hv -> Products                                        =*
*=  Cross section: from JPL 97 recommendation                                =*
*=  Quantum yield: assumed to be unity                                       =*
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
      PARAMETER(kdata=100)

      REAL x1(kdata)
      REAL y1(kdata)

* local

      REAL yg(kw)
C      REAL qy
      INTEGER i, iw, n, idum
      INTEGER iz


***************** CF2BrCl (Halon-1211) photolysis *************************

      j = j+1
      jlabel(j) = 'CF2BrCl (Halon-1211) -> Products'

**** cross sections from JPL97 recommendation (identical to 94 recommendation)

      OPEN(kin,FILE='DATAJ1/ABS/Halon-1211_jpl97.abs',STATUS='OLD')
      READ(kin,*) idum, n
      DO i = 1, idum-2
        READ(kin,*)
      ENDDO
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
        y1(i) = y1(i) * 1E-20
      ENDDO
      CLOSE(kin)

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)

**** quantum yield assumed to be unity
C      qy = 1.

      DO iw = 1, nw-1
        DO iz = 1, nz
           sq(j,iz,iw) = yg(iw)! * qy
        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r42(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! CF3Br

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for Halon-1301         =*
*=  photolysis:                                                              =*
*=         CF3Br + hv -> Products                                            =*
*=  Cross section: from JPL 97 recommendation                                =*
*=  Quantum yield: assumed to be unity                                       =*
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
      PARAMETER(kdata=100)

      REAL x1(kdata)
      REAL y1(kdata)

* local

      REAL yg(kw)
C      REAL qy
      INTEGER i, iw, n, idum
      INTEGER iz


***************** CF3Br (Halon-1301) photolysis *************************

      j = j+1
      jlabel(j) = 'CF3Br (Halon-1301) -> Products'

**** cross sections from JPL97 recommendation (identical to 94 recommendation)

      OPEN(kin,FILE='DATAJ1/ABS/Halon-1301_jpl97.abs',STATUS='OLD')
      READ(kin,*) idum, n
      DO i = 1, idum-2
        READ(kin,*)
      ENDDO
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
        y1(i) = y1(i) * 1E-20
      ENDDO
      CLOSE(kin)
      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)

**** quantum yield assumed to be unity
C      qy = 1.

      DO iw = 1, nw-1
        DO iz = 1, nz
           sq(j,iz,iw) = yg(iw)! * qy
        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r43(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! CF2BrCF2Br

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for Halon-2402         =*
*=  photolysis:                                                              =*
*=           CF2BrCF2Br + hv -> Products                                     =*
*=  Cross section: from JPL 97 recommendation                                =*
*=  Quantum yield: assumed to be unity                                       =*
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
      PARAMETER(kdata=100)

      REAL x1(kdata)
      REAL y1(kdata)

* local

      REAL yg(kw)
C      REAL qy
      INTEGER i, iw, n, idum
      INTEGER iz


***************** CF2BrCF2Br (Halon-2402) photolysis *************************

      j = j+1
      jlabel(j) = 'CF2BrCF2Br (Halon-2402) -> Products'

**** cross sections from JPL97 recommendation (identical to 94 recommendation)

      OPEN(kin,FILE='DATAJ1/ABS/Halon-2402_jpl97.abs',STATUS='OLD')
      READ(kin,*) idum, n
      DO i = 1, idum-2
        READ(kin,*)
      ENDDO
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
        y1(i) = y1(i) * 1E-20
      ENDDO
      CLOSE(kin)
      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)

**** quantum yield assumed to be unity
C      qy = 1.

      DO iw = 1, nw-1
        DO iz = 1, nz
           sq(j,iz,iw) = yg(iw)! * qy
        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r44(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! N2O

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for N2O photolysis:    =*
*=              N2O + hv -> N2 + O(1D)                                       =*
*=  Cross section: from JPL 97 recommendation                                =*
*=  Quantum yield: assumed to be unity, based on Greenblatt and Ravishankara =*
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

* local

C      REAL qy
      REAL a, b
      REAL a0, a1, a2, a3, a4
      REAL b0, b1, b2, b3
      REAL t
      INTEGER iw, iz
      REAL lambda


***************** N2O photolysis *************************

      j = j+1
      jlabel(j) = 'N2O -> N2 + O(1D)'

**** cross sections according to JPL97 recommendation (identical to 94 rec.)
**** see file DATAJ1/ABS/N2O_jpl94.abs for detail

      A0 = 68.21023
      A1 = -4.071805
      A2 = 4.301146E-02
      A3 = -1.777846E-04
      A4 = 2.520672E-07

      B0 = 123.4014
      B1 = -2.116255
      B2 = 1.111572E-02
      B3 = -1.881058E-05

**** quantum yield of N(4s) and NO(2Pi) is less than 1% (Greenblatt and
**** Ravishankara), so quantum yield of O(1D) is assumed to be unity
C      qy = 1.

      DO iw = 1, nw-1
         lambda = wc(iw)
         IF (lambda .GE. 173. .AND. lambda .LE. 240.) THEN
           DO iz = 1, nz
             t = MAX(194.,MIN(tlev(iz),320.))
             A = (((A4*lambda+A3)*lambda+A2)*lambda+A1)*lambda+A0
             B = (((B3*lambda+B2)*lambda+B1)*lambda+B0)
             B = (t-300.)*EXP(B)
             sq(j,iz,iw) = EXP(A+B)! * qy
           ENDDO
         ELSE
           DO iz = 1, nz
             sq(j,iz,iw) = 0.
           ENDDO
         ENDIF
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r45(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! ClONO2

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for ClONO2 photolysis: =*
*=        ClONO2 + hv -> Products                                            =*
*=                                                                           =*
*=  Cross section: JPL 97 recommendation                                     =*
*=  Quantum yield: JPL 97 recommendation                                     =*
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

      REAL x1(kdata),x2(kdata),x3(kdata)
      REAL y1(kdata),y2(kdata),y3(kdata)
      INTEGER n1, n2, n3

* local

      REAL yg1(kw), yg2(kw), yg3(kw)
      REAL qy1, qy2
      REAL xs
      INTEGER i, iw, n
      INTEGER iz


***************** ClONO2 photolysis *************************

      j = j+1
      jlabel(j) = 'ClONO2 -> Cl + NO3'
      j = j+1
      jlabel(j) = 'ClONO2 -> ClO + NO2'


*** cross sections from JPL97 recommendation.  Same in JPL-2011.

      OPEN(kin,FILE='DATAJ1/ABS/ClONO2_jpl97.abs',STATUS='OLD')
      DO i = 1, 2
         READ(kin,*)
      ENDDO
      n = 119
      DO i = 1, n
        READ(kin,*) x1(i), y1(i), y2(i), y3(i)
        y1(i) = y1(i) * 1E-20
        x2(i) = x1(i)
        x3(i) = x1(i)
      ENDDO
      CLOSE(kin)

      n1 = n
      CALL interpol(x1,y1,kdata,n1,nw,wl,jlabel(j),0.,0.,yg1)
      n2 = n
      CALL interpol(x2,y2,kdata,n2,nw,wl,jlabel(j),0.,0.,yg2)
      n3 = n
      CALL interpol(x3,y3,kdata,n3,nw,wl,jlabel(j),0.,0.,yg3)

      DO iw = 1, nw-1

*** quantum yields (from jpl97, same in jpl2011)

         IF( wc(iw) .LT. 308.) THEN
            qy1 = 0.6
         ELSEIF( (wc(iw) .GE. 308) .AND. (wc(iw) .LE. 364.) ) THEN
            qy1 = 7.143e-3 * wc(iw) - 1.6
         ELSEIF( wc(iw) .GT. 364. ) THEN
            qy1 = 1.0
         ENDIF
         qy2 = 1. - qy1

* compute T-dependent cross section

         DO iz = 1, nz
            xs = yg1(iw)*( 1. +
     $           yg2(iw)*(tlev(iz)-296) +
     $           yg3(iw)*(tlev(iz)-296)*(tlev(iz)-296))
            sq(j-1,iz,iw) = qy1 * xs
            sq(j,iz,iw) = qy2 * xs

         ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r46(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! BrONO2

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for BrONO2 photolysis: =*
*=        BrONO2 + hv -> Products                                            =*
*=                                                                           =*
*=  Cross section: JPL 03 recommendation                                     =*
*=  Quantum yield: JPL 03 recommendation                                     =*
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
      PARAMETER(kdata=100)

      REAL x1(kdata)
      REAL y1(kdata)

* local

      REAL yg1(kw)
      REAL qyNO2, qyNO3
      INTEGER i, iw, n
      INTEGER iz


***************** BrONO2 photolysis *************************

      j = j+1
      jlabel(j) = 'BrONO2 -> BrO + NO2'
      j = j+1
      jlabel(j) = 'BrONO2 -> Br + NO3'


*** cross sections from JPL03 recommendation

      OPEN(kin,FILE='DATAJ1/ABS/BrONO2_jpl03.abs',STATUS='OLD')
      DO i = 1, 13
         READ(kin,*)
      ENDDO
      n = 61
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
        y1(i) = y1(i) * 1E-20
      ENDDO
      CLOSE(kin)

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg1)

*** quantum yields (from jpl97)

      qyNO2 = 0.15
      qyNO3 = 0.85
      DO iw = 1, nw-1
         DO iz = 1, nz
            sq(j-1,iz,iw) = qyNO2 * yg1(iw)
            sq(j,iz,iw) = qyNO3 * yg1(iw)
         ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r47(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! Cl2

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for Cl2 photolysis:    =*
*=        Cl2 + hv -> 2 Cl                                                   =*
*=                                                                           =*
*=  Cross section: JPL 97 recommendation                                     =*
*=  Quantum yield: 1     (Calvert and Pitts, 1966)                           =*
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
C      REAL qy
      INTEGER iz, iw

      real aa, bb, ex1, ex2, sig, alpha(kz)

      INTEGER mabs, xsvers(3)
      CHARACTER(llog) :: logmsg(nlog)
      CHARACTER(lspc) :: spc


***************** Cl2 photolysis *************************

* Setting photolysis index and cross section options
      j = j+1
      jlabel(j) = 'Cl2 -> Cl + Cl'

* cross section options
      spc = "Cl2"
      xsvers = (/2, 2, 1/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg( 1) = "from Finlayson-Pitts and Pitts"
      logmsg( 2) = "from JPL 2011 formula"
      CALL set_option(2,spc,"abs",logmsg,xsvers,mabs)


      IF (mabs .EQ. 1) THEN

*** cross sections from JPL97 recommendation (as tab by Finlayson-Pitts
* and Pitts, 1999.

         OPEN(kin,FILE='DATAJ1/ABS/CL2_fpp.abs',STATUS='OLD')
         do i = 1, 5
            read(kin,*)
         enddo
         n = 22
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
            y1(i) = y1(i) * 1E-20
         ENDDO
         CLOSE(kin)

         CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)

      ELSEIF(mabs .EQ. 2) THEN

         DO iz = 1, nz
            aa = 402.7/tlev(iz)
            bb = exp(aa)
            alpha(iz) = (bb - 1./bb) / (bb + 1./bb)
         ENDDO

      ENDIF


*** quantum yield = 1 (Calvert and Pitts, 1966)
C      qy = 1.

      DO iw = 1, nw-1

         if(mabs .eq. 1) sig = yg(iw)

         DO iz = 1, nz

            if (mabs .eq. 2) then

         ex1 = 27.3  * exp(-99.0 * alpha(iz) * (log(329.5/wc(iw)))**2)
         ex2 = 0.932 * exp(-91.5 * alpha(iz) * (log(406.5/wc(iw)))**2)
         sig = 1e-20 * alpha(iz)**0.5 * (ex1 + ex2)

            ENDIF

            sq(j,iz,iw) = sig! * qy

         ENDDO
      ENDDO

      RETURN
      END

*=============================================================================*

      SUBROUTINE r101(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! CH2(OH)CHO

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for CH2(OH)CHO     =*
*=  (glycolaldehye, hydroxy acetaldehyde) photolysis:                        =*
*=           CH2(OH)CHO + hv -> Products                                     =*
*=                                                                           =*
*=  Quantum yield about 50%                                                  =*
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

      CHARACTER(lcl) jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:

      INTEGER j

* data arrays

      INTEGER kdata
      PARAMETER(kdata=300)

      INTEGER i, n
      REAL x(kdata), y(kdata)

* local

      REAL yg(kw)
      INTEGER iw
      REAL qy1, qy2, qy3

      INTEGER mabs, myld, xsvers(3), qyvers(3)
      CHARACTER(llog) :: logmsg(nlog)
      CHARACTER(lspc) :: spc


***************** glycolaldehyde photolysis *************************

* Setting photolysis index and cross section/quantum yield options
      j = j+1
      jlabel(j) = 'HOCH2CHO -> CH2OH + HCO'
      j = j+1
      jlabel(j) = 'HOCH2CHO -> CH3OH + CO'
      j = j+1
      jlabel(j) = 'HOCH2CHO -> CH2CHO + OH'

* cross section options
      spc = "CH2(OH)CHO"
      xsvers = (/2, 1, 1/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "from IUPAC"
      logmsg(2) = "from JPL"
      CALL set_option(2,spc,"abs",logmsg,xsvers,mabs)

* quantum yield options
      qyvers = (/1, 3, 2/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "from JPL"
      logmsg(2) = "from IUPAC"
      logmsg(3) = "scaled externally"
      CALL set_option(3,spc,"yld",logmsg,qyvers,myld)


      IF(mabs .EQ. 1) THEN

*=  Cross section from                                                       =*
*= The Atmospheric Chemistry of Glycolaldehyde, C. Bacher, G. S. Tyndall     =*
*= and J. J. Orlando, J. Atmos. Chem., 39 (2001) 171-189.                    =*
*= (JPL vs. IUPAC recommendation)                                            =*


         OPEN(UNIT=kin,FILE='DATAJ1/CH2OHCHO/glycolaldehyde.abs',
     $        STATUS='old')
         DO i = 1, 15
            READ(kin,*)
         ENDDO
         n = 131
         DO i = 1, n
            READ(kin,*) x(i), y(i)
         ENDDO
         CLOSE(kin)

         CALL interpol(x,y,kdata,n,nw,wl,jlabel(j),0.,0.,yg)

      ELSEIF(mabs .EQ. 2) THEN

         OPEN(UNIT=kin,FILE='DATAJ1/CH2OHCHO/glycolaldehyde_jpl11.abs',
     $        STATUS='old')
         DO i = 1, 2
            READ(kin,*)
         ENDDO
         n = 63
         DO i = 1, n
            READ(kin,*) x(i), y(i)
            y(i) = y(i) * 1.e-20
         ENDDO
         CLOSE(kin)

         CALL interpol(x,y,kdata,n,nw,wl,jlabel(j),0.,0.,yg)

      ENDIF

* Quantum yields now scaled externally
      IF (myld == 1) THEN
        qy1 = 0.83
        qy2 = 0.1
        qy3 = 0.07
       ELSEIF (myld == 2) THEN
        qy1 = 0.75
        qy2 = 0.
        qy3 = 0.
       ELSEIF(myld == 3) THEN
        qy1 = 1.0
        qy2 = 1.0
        qy3 = 1.0
      ENDIF


* combine:
      DO iw = 1, nw - 1
         DO i = 1, nz
           sq(j-2,i,iw) = yg(iw) * qy1
           sq(j-1,i,iw) = yg(iw) * qy2
           sq(j  ,i,iw) = yg(iw) * qy3
         ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r102(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! CH3COCOCH3

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for CH3COCOCH3     =*
*=  (biacetyl) photolysis:                                                   =*
*=           CH3COCOCH3 + hv -> Products                                     =*
*=                                                                           =*
*=  Cross section from either                                                =*
*= 1.  Plum et al., Environ. Sci. Technol., Vol. 17, No. 8, 1983, p.480      =*
*= 2.  Horowitz et al., J. Photochem Photobio A, 146, 19-27, 2001.           =*
*=                                                                           =*
*=  Quantum yield =0.158                                                     =*
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
      PARAMETER(kdata=300)

      INTEGER i, n
      REAL x(kdata), y(kdata)

* local

      REAL yg(kw)
      REAL qy
      INTEGER iw

      INTEGER mabs, myld, xsvers(3), qyvers(3)
      CHARACTER(llog) :: logmsg(nlog)
      CHARACTER(lspc) :: spc


***************** biacetyl photolysis *************************

* Setting photolysis index and cross section options
      j = j+1
      jlabel(j) = 'CH3COCOCH3 -> Products'

* cross section options
      spc = "CH3COCOCH3"
      xsvers = (/2, 2, 1/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg( 1) = "from Plum et al"
      logmsg( 2) = "from Horowitz et al"
      CALL set_option(2,spc,"abs",logmsg,xsvers,mabs)

      qyvers = (/1, 2, 1/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg( 2) = "scaled externally" ! opt 1, same as xs option
      CALL set_option(2,spc,"yld",logmsg,qyvers,myld)


* cross sections

      IF( mabs. EQ. 1) THEN
         OPEN(UNIT=kin,FILE='DATAJ1/CH3COCOCH3/biacetyl_plum.abs',
     $        STATUS='old')
         DO i = 1, 7
            READ(kin,*)
         ENDDO
         n = 55
         DO i = 1, n
            READ(kin,*) x(i), y(i)
            y(i) = y(i) * 1.e-20
         ENDDO
         CLOSE(kin)

      ELSEIF(mabs. EQ. 2) THEN
         OPEN(UNIT=kin,FILE='DATAJ1/CH3COCOCH3/biacetyl_horowitz.abs',
     $        STATUS='old')
         DO i = 1, 8
            READ(kin,*)
         ENDDO
         n = 287
         DO i = 1, n
            READ(kin,*) x(i), y(i)
            y(i) = y(i) * 1.e-20
         ENDDO
         CLOSE(kin)

      ENDIF

      CALL interpol(x,y,kdata,n,nw,wl,jlabel(j),0.,0.,yg)

* quantum yield from Plum et al.
      IF(myld==1) THEN
        qy = 0.158
       ELSEIF(myld==2) THEN
        qy = 1.
      ENDIF

      DO iw = 1, nw - 1
         DO i = 1, nz
            sq(j,i,iw) = yg(iw) * qy
         ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r103(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! MVK

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for CH3COCHCH2     =*
*=  Methyl vinyl ketone photolysis:                                          =*
*=           CH3COCH=CH2 + hv -> Products                                    =*
*=                                                                           =*
*=  Cross section from                                                       =*
*= W. Schneider and G. K. Moorgat, priv. comm, MPI Mainz 1989 as reported by =*
*= Roeth, E.-P., R. Ruhnke, G. Moortgat, R. Meller, and W. Schneider,        =*
*= UV/VIS-Absorption Cross Sections and QUantum Yields for Use in            =*
*= Photochemistry and Atmospheric Modeling, Part 2: Organic Substances,      =*
*= Forschungszentrum Julich, Report Jul-3341, 1997.                          =*
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
      PARAMETER(kdata=20000)

      INTEGER i, n
      REAL x(kdata), y(kdata)

* local

      REAL yg(kw)
      REAL qy, qy1, qy2, qy3
      INTEGER iw

      INTEGER mabs, myld, xsvers(3), qyvers(3)
      CHARACTER(llog) :: logmsg(nlog)
      CHARACTER(lspc) :: spc


***************** MVK photolysis *************************

* Setting photolysis index and cross section options
      j = j+1
      jlabel(j) = 'CH3COCH=CH2 -> CH3 + C2H3CO'
      j = j+1
      jlabel(j) = 'CH3COCH=CH2 -> C2H3 + CH3CO'
      j = j+1
      jlabel(j) = 'CH3COCH=CH2 -> C3H6 + CO'

* cross section options
      spc = "MVK"
      xsvers = (/2, 2, 1/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg( 1) = "from Schneider and Moortgat"
      logmsg( 2) = "from JPL/IUPAC"
      CALL set_option(2,spc,"abs",logmsg,xsvers,mabs)

* quantum yield options
      qyvers = (/1, 2, 2/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = 'from IUPAC'
      logmsg(2) = 'from IUPAC with branching scaled externally'
      CALL set_option(2,spc,"yld",logmsg,qyvers,myld)


* cross sections

      IF(mabs .EQ. 1) THEN

         OPEN(UNIT=kin,FILE='DATAJ1/ABS/MVK_schneider.abs',
     $        STATUS='old')
         DO i = 1, 9
            READ(kin,*)
         ENDDO
         n = 19682
         DO i = 1, n
            READ(kin,*) x(i), y(i)
         ENDDO
         CLOSE(kin)

         CALL interpol(x,y,kdata,n,nw,wl,jlabel(j),0.,0.,yg)

      ELSEIF(mabs .EQ. 2) THEN

         OPEN(UNIT=kin,FILE='DATAJ1/ABS/MVK_jpl11.abs',
     $        STATUS='old')
         DO i = 1, 2
            READ(kin,*)
         ENDDO
         n = 146
         DO i = 1, n
            READ(kin,*) x(i), y(i)
            y(i) = y(i) * 1.e-20
         ENDDO
         CLOSE(kin)

         CALL interpol(x,y,kdata,n,nw,wl,jlabel(j),0.,0.,yg)
      ENDIF

* quantum yield from
* Gierczak, T., J. B. Burkholder, R. K. Talukdar, A. Mellouki, S. B. Barone,
* and A. R. Ravishankara, Atmospheric fate of methyl vinyl ketone and methacrolein,
* J. Photochem. Photobiol A: Chemistry, 110 1-10, 1997.
* depends on pressure and wavelength, set upper limit to 1.0
* also recommended by IUPAC

* branching of 0.6:0.4 for CH3 and vinyl bond fission estimated

      DO iw = 1, nw - 1
         DO i = 1, nz
            qy  = exp(-0.055*(wc(iw)-308.)) /
     $           (5.5 + 9.2e-19*airden(i))
            qy  = min(qy, 1.)
            IF(myld==1) THEN
              qy1 = 0.2 * qy
              qy2 = 0.2 * qy
              qy3 = 0.6 * qy
            ELSEIF(myld==2) THEN
              qy1 = qy
              qy2 = qy
              qy3 = qy
            ENDIF
            sq(j-2,i,iw) = yg(iw) * qy1
            sq(j-1,i,iw) = yg(iw) * qy2
            sq(j  ,i,iw) = yg(iw) * qy3
         ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r104(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! MACR

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for CH2=C(CH3)CHO      =*
*=  (methacrolein) photolysis:                                               =*
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
      real qy, qy1, qy2, qy3, qy4!, qyHPALD
      REAL sig

      INTEGER myld, mprod, qyvers(3), prvers(3)
      CHARACTER(llog) :: logmsg(nlog)
      CHARACTER(lspc) :: spc


***************** methacrolein photolysis *************************

* Setting photolysis index and quantum yield/product distribution options
* same channels estimated as for acrolein
      j = j+1
      jlabel(j) = 'CH2=C(CH3)CHO -> CH2=CCH3 + CHO'
      j = j+1
      jlabel(j) = 'CH2=C(CH3)CHO -> CH3CH=CH2 + CO'
      j = j+1
      jlabel(j) = 'CH2=C(CH3)CHO -> CH3C.(OO.)CH3 + CO'
      j = j+1
      jlabel(j) = 'CH2=C(CH3)CHO -> CH2=C(CH3)CO + H'
      j = j+1
      jlabel(j) = 'HPALD -> Products'

* quantum yield options
      spc = "MACR"
      qyvers = (/1, 2, 2/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "estimated with 0.1 (upper limit)"
      logmsg(2) = "from Calvert et al. 2011 book"
      CALL set_option(2,spc,"yld",logmsg,qyvers,myld)

* product distribution options
      prvers = (/1, 2, 2/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "with estimated Criegie formation in main channel"
      logmsg(2) = "with estimated propylene formation in main channel"
      CALL set_option(2,spc,"prd",logmsg,prvers,mprod)


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

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j-3),0.,0.,yg)


* product distribution
* 1: assuming CH3CH: + O2 -> CH3CHOO for main channel
* 2: assuming CH3CH: + M -> C2H4 + M for main channel

      IF(myld==2) THEN
        OPEN(UNIT=kin,FILE='DATAJ1/ABS/MACR_calv.dat',
     $       STATUS='OLD')
        DO i = 1, 35
          READ(kin,*)
        ENDDO
        n = 90
        DO i = 1, n
          READ(kin,*) x1(i), dum, y1(i)
        ENDDO
        CLOSE(kin)

        CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j-4),0.,0.,yg1)

      ENDIF

! QY for HPALD:
C      qyHPALD = 1.0

      DO iw = 1, nw-1

        sig = yg(iw)
* quantum yields
        IF(myld==1) THEN
          qy = 0.01
         ELSEIF(myld==2) THEN
          qy = yg1(iw)
        ENDIF

        DO i = 1, nz
           IF(mprod == 1) THEN
             qy1 = qy*(-0.0173*(airden(i)/1.E19)**2
     &                +0.083*airden(i)/1.E19+0.0492)
             qy2 = qy*(-0.0111*(airden(i)/1.E19)**2
     &                +0.1014*airden(i)/1.E19+0.0532)
             qy3 = qy*(0.0518*(airden(i)/1.E19)**2
     &                -0.2675*airden(i)/1.E19+0.7953)
             qy4 = qy*(-0.0217*(airden(i)/1.E19)**2
     &                +0.0788*airden(i)/1.E19+0.1029)
            ELSEIF(mprod == 2) THEN
             qy1 = qy*(-0.0173*(airden(i)/1.E19)**2
     &                +0.083*airden(i)/1.E19+0.0492)
             qy2 = qy*(0.0407*(airden(i)/1.E19)**2
     &                -0.1661*airden(i)/1.E19+0.8485)
             qy3 = 0.0
             qy4 = qy*(-0.0217*(airden(i)/1.E19)**2
     &                +0.0788*airden(i)/1.E19+0.1029)
           ENDIF

* combine
           sq(j-4,i,iw) = qy1 * sig
           sq(j-3,i,iw) = qy2 * sig
           sq(j-2,i,iw) = qy3 * sig
           sq(j-1,i,iw) = qy4 * sig
           sq(j  ,i,iw) = sig!* qyHPALD

        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r105(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! CH3COCO(OH)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for CH3COCO(OH)    =*
*=  pyruvic acid        photolysis:                                          =*
*=           CH3COCO(OH) + hv -> Products                                    =*
*=                                                                           =*
*=  Cross section from                                                       =*
*= Horowitz, A., R. Meller, and G. K. Moortgat, The UV-VIS absorption cross  =*
*= section of the a-dicarbonyl compounds: pyruvic acid, biacetyl, and        =*
*= glyoxal. J. Photochem. Photobiol. A:Chemistry, v.146, pp.19-27, 2001.     =*
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
      PARAMETER(kdata=20000)

      INTEGER i, n
      REAL x(kdata), y(kdata)

* local

      REAL yg(kw)
      REAL qy1, qy2, qy3, qy4
      INTEGER iw

      INTEGER mabs, myld, xsvers(3), qyvers(3)
      CHARACTER(llog) :: logmsg(nlog)
      CHARACTER(lspc) :: spc


***************** pyruvic acid photolysis *************************

* Setting photolysis index and cross section/quantum yield options
      j = j+1
      jlabel(j) = 'CH3COCOOH -> CH3CHO + CO2'
      j = j+1
      jlabel(j) = 'CH3COCOOH -> CH3CO + COOH'
      j = j+1
      jlabel(j) = 'CH3COCOOH -> CH3COOH + CO'
      j = j+1
      jlabel(j) = 'CH3COCOOH -> CH3CO + CO + OH'

* cross section options
      spc = "CH3COCOOH"
      xsvers = (/1, 1, 2/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "from Horowitz et al"
      logmsg(2) = "from JPL 2011"
      CALL set_option(2,spc,"abs",logmsg,xsvers,mabs)

* quantum yield options
      qyvers = (/1, 1, 1/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "from TUV with no products identified" !also used for external scaling
      logmsg(2) = "from JPL 2011"
      CALL set_option(2,spc,"yld",logmsg,qyvers,myld)


* cross sections

      IF (mabs .EQ. 1) THEN

         OPEN(UNIT=kin,FILE='DATAJ1/CH3COCOOH/pyruvic_horowitz.abs',
     $        STATUS='old')
         DO i = 1, 8
            READ(kin,*)
         ENDDO
         n = 148
         DO i = 1, n
            READ(kin,*) x(i), y(i)
            y(i) = y(i) * 1.e-20
         ENDDO
         CLOSE(kin)

         CALL interpol(x,y,kdata,n,nw,wl,jlabel(j),0.,0.,yg)

      ELSEIF (mabs .eq. 2) then

         OPEN(UNIT=kin,FILE='DATAJ1/CH3COCOOH/pyruvic_jpl11.abs',
     $        STATUS='old')
         DO i = 1, 2
            READ(kin,*)
         ENDDO
         n = 139
         DO i = 1, n
            READ(kin,*) x(i), y(i)
            y(i) = y(i) * 1.e-20
         ENDDO
         CLOSE(kin)

         CALL interpol(x,y,kdata,n,nw,wl,jlabel(j),0.,0.,yg)

      ENDIF


* quantum yields

      IF(myld==1) THEN
        qy1 = 1.
        qy2 = 0.
        qy3 = 0.
        qy4 = 0.
       ELSEIF(myld==2) THEN
        qy1 = 0.48
        qy2 = 0.39
        qy3 = 0.08
        qy4 = 0.05
      ENDIF

* combine
      DO iw = 1, nw - 1
         DO i = 1, nz
            sq(j-3,i,iw) = yg(iw) * qy1
            sq(j-2,i,iw) = yg(iw) * qy2
            sq(j-1,i,iw) = yg(iw) * qy3
            sq(j  ,i,iw) = yg(iw) * qy4
         ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r107(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! CH3CHONO2CH3

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for CH3CHONO2CH3   =*
*=  isopropyl nitrate   photolysis:                                          =*
*=           CH3CHONO2CH3 + hv -> CH3CHOCH3 + NO2                            =*
*=                                                                           =*
*= Absorption cross sections as given below in source code                   =*
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

      INTEGER i, n, n1, n2
      REAL x1(kdata), y1(kdata)
      REAL x2(kdata), y2(kdata)

* local

      REAL dum
      REAL yg1(kw), yg2(kw)
      REAL sig!, qy
      INTEGER iw

      INTEGER mabs, xsvers(3)
      CHARACTER(llog) :: logmsg(nlog)
      CHARACTER(lspc) :: spc


***************** CH3CHONO2CH3 photolysis *************************

* Setting photolysis index and cross section options
      j = j+1
      jlabel(j) = 'CH3CHONO2CH3 -> CH3CHOCH3 + NO2'

* cross section options
      spc = "CH3COCOOH"
      xsvers = (/1, 2, 2/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "from Talukdar et al"
      logmsg(2) = "from IUPAC"
      CALL set_option(2,spc,"abs",logmsg,xsvers,mabs)


      IF(mabs==1) THEN
        OPEN(UNIT=kin,FILE='DATAJ1/RONO2/RONO2_talukdar.abs',
     $       STATUS='old')
        DO i = 1, 10
           READ(kin,*)
        ENDDO
        n1 = 0
        n2 = 0
        DO i = 1, 63
           READ(kin,*) x1(i), dum, dum, dum, dum, y1(i), y2(i)
           if (y1(i) .gt. 0.) n1 = n1 + 1
           if (y2(i) .gt. 0.) n2 = n2 + 1
           x2(i) = x1(i)
           y1(i) = y1(i) * 1.e-20
           y2(i) = y2(i) * 1.e-3
        ENDDO
        CLOSE(kin)

       ELSEIF(mabs==2) THEN

        OPEN(UNIT=kin,FILE='DATAJ1/RONO2/iC3H7NO3_iup.abs',
     $       STATUS='old')
        DO i = 1, 5
           READ(kin,*)
        ENDDO
        n = 37
        n1 = n
        n2 = n
        DO i = 1, n
           READ(kin,*) x1(i), y1(i), y2(i)
           x2(i) = x1(i)
           y1(i) = y1(i) * 1.e-20
           y2(i) = y2(i) * 1.e-3
        ENDDO
        CLOSE(kin)
      ENDIF

      CALL interpol(x1,y1,kdata,n1,nw,wl,jlabel(j),0.,0.,yg1)
      CALL interpol(x2,y2,kdata,n2,nw,wl,jlabel(j),y2(1),y2(n2),yg2)


* quantum yield  = 1
!      qy = 1.

* Map onto grid
      DO iw = 1, nw - 1
         DO i = 1, nz
            sig = yg1(iw)*exp(yg2(iw)*(tlev(i)-298.))
            sq(j,i,iw) = sig ! * qy
         ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r108(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! CH2(OH)CH2(ONO2)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for                =*
*=   nitroxy ethanol CH2(OH)CH2(ONO2) + hv -> CH2(OH)CH2(O.) + NO2           =*
*=                                                                           =*
*=  Cross section from Roberts, J. R. and R. W. Fajer, UV absorption cross   =*
*=    sections of organic nitrates of potential atmospheric importance and   =*
*=    estimation of atmospheric lifetimes, Env. Sci. Tech., 23, 945-951,     =*
*=    1989.
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

* local

      REAL sig!, qy
      INTEGER iw, i
      REAL a, b, c


***************** CH2(OH)CH2(ONO2) photolysis *************************

      j = j+1
      jlabel(j) = 'CH2(OH)CH2(ONO2) -> CH2(OH)CH2(O.) + NO2'


* coefficients from Roberts and Fajer 1989, over 270-306 nm

      a = -2.359E-3
      b = 1.2478
      c = -210.4

* quantum yield  = 1
!      qy = 1.

      DO iw = 1, nw - 1
         sig = EXP(a*wc(iw)*wc(iw) + b*wc(iw) + c)
         DO i = 1, nz
            sq(j,i,iw) = sig! * qy
         ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r109(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! CH3COCH2(ONO2)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for                =*
*=   nitroxy acetone CH3COCH2(ONO2) + hv -> CH3COCH2(O.) + NO2               =*
*=                                                                           =*
*=  Cross section:  see options below                                        =*
*=  Quantum yield:  see options below                                        =*
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
      REAL sig, sig1, sig2
      INTEGER iw

      INTEGER mabs, myld, xsvers(3), qyvers(3)
      CHARACTER(llog) :: logmsg(nlog)
      CHARACTER(lspc) :: spc


***************** nitroxy acetone photolysis *************************

* Setting photolysis index and cross section/quantum yield options
      j = j+1
      jlabel(j) = 'CH3COCH2(ONO2) -> CH3COCH2(O.) + NO2'

* cross section options
      spc = "CH3COCH2(ONO2)"
      xsvers = (/2, 6, 2/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "from least square fit by Barnes et al"
      logmsg(2) = "from data by Barnes et al"
      logmsg(3) = "from least square fit by Roberts and Fajer"
      logmsg(4) = "from data by Roberts and Fajer"
      logmsg(5) =
     &   "mean of least square fits by Barnes et al./Roberts and Fajer"
      logmsg(6) = "mean of data by Barnes et al./Roberts and Fajer"
      CALL set_option(6,spc,"abs",logmsg,xsvers,mabs)

* quantum yield options
      qyvers = (/2, 2, 1/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "estimated unity (Calvert et al.)"
      logmsg(2) = "estimated 0.9 (Muller et al.)"
      CALL set_option(2,spc,"yld",logmsg,qyvers,myld)


* cross sections

      IF(mabs==2) THEN
        OPEN(UNIT=kin,FILE='DATAJ1/MCMext/POLY/NOA_Bar93.abs',
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

      ELSEIF(mabs==4) THEN
        OPEN(UNIT=kin,FILE='DATAJ1/MCMext/POLY/NOA_R+F89.abs',
     $       STATUS='old')
        do i = 1, 6
          read(kin,*)
        enddo

        n  = 17
        DO i = 1, n
          READ(kin,*) x1(i), y1(i)
        ENDDO
        CLOSE(kin)
        CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)
      ELSEIF(mabs==6) THEN
        OPEN(UNIT=kin,FILE='DATAJ1/MCMext/POLY/NOA_mean.abs',
     $       STATUS='old')
        do i = 1, 9
          read(kin,*)
        enddo

        n  = 20
        DO i = 1, n
          READ(kin,*) x1(i), y1(i)
          y1(i) = y1(i) * 1.e-20
        ENDDO
        CLOSE(kin)
        CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)

      ENDIF


* quantum yields

      IF(myld==1) THEN
        qy = 1.0
       ELSEIF(myld==2) THEN
        qy = 0.9
      ENDIF


* combine:

      DO iw = 1, nw - 1
        sig1 = exp(-9.804e-4*wc(iw)**2+0.5404*wc(iw)-118.3)
        sig2 = exp(-1.365e-3*wc(iw)**2+0.7834*wc(iw)-156.8)
        IF(mabs==1)THEN
          sig = sig1
         ELSEIF(mabs==3)THEN
          sig = sig2
         ELSEIF(mabs==5)THEN
          sig = (sig1+sig2)/2.
         ELSEIF(mabs==2 .or. mabs==4 .or. mabs==6)THEN
          sig = yg(iw)
        ENDIF
        DO i = 1, nz
          sq(j,i,iw) = sig * qy
        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r110(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! C(CH3)3(ONO2)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for                =*
*=  t-butyl nitrate C(CH3)3(ONO2) + hv -> C(CH3)(O.) + NO2                   =*
*=                                                                           =*
*=  Cross section from Roberts, J. R. and R. W. Fajer, UV absorption cross   =*
*=    sections of organic nitrates of potential atmospheric importance and   =*
*=    estimation of atmospheric lifetimes, Env. Sci. Tech., 23, 945-951,     =*
*=    1989.                                                                  =*
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

* local

      REAL sig!, qy
      INTEGER iw, i
      REAL a, b, c


***************** C(CH3)3(ONO2) photolysis *************************

      j = j+1
      jlabel(j) = 'C(CH3)3(ONO2) -> C(CH3)3(O.) + NO2'


* coefficients from Roberts and Fajer 1989, over 270-330 nm

      a = -0.993E-3
      b = 0.5307
      c = -115.5

* quantum yield  = 1
!      qy = 1.

      DO iw = 1, nw - 1
         sig = EXP(a*wc(iw)**2 + b*wc(iw) + c)
         DO i = 1, nz
            sq(j,i,iw) = sig! * qy
         ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r111(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! ClOOCl

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for ClOOCl         =*
*=  ClO dimer           photolysis:                                          =*
*=           ClOOCl + hv -> Cl + ClOO                                        =*
*=                                                                           =*
*=  Cross section from  JPL2011                                              =*
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
      PARAMETER(kdata=20000)

      INTEGER i, n
      REAL x(kdata), y(kdata)

* local

      REAL yg(kw)
C      REAL qy
      INTEGER iw


***************** ClOOCl photolysis *************************

* from JPL-2011

      j = j+1
      jlabel(j) = 'ClOOCl -> Cl + ClOO'

      OPEN(UNIT=kin,FILE='DATAJ1/ABS/ClOOCl_jpl11.abs',
     $     STATUS='old')
      DO i = 1, 3
         READ(kin,*)
      ENDDO
      n = 111
      DO i = 1, n
         READ(kin,*) x(i), y(i)
         y(i) = y(i) * 1.e-20
      ENDDO
      CLOSE(kin)

      CALL interpol(x,y,kdata,n,nw,wl,jlabel(j),0.,0.,yg)

* quantum yield  = 1
!      qy = 1.

      DO iw = 1, nw - 1
         DO i = 1, nz
            sq(j,i,iw) = yg(iw)! * qy
         ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r112(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! CH2(OH)COCH3

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for hydroxyacetone =*
*=  CH2(OH)COCH3 photolysis:                                                 =*
*=           CH2(OH)COCH3 + hv -> CH3CO + CH2OH                              =*
*=                             -> CH2(OH)CO + CH3                            =*
*=                                                                           =*
*=  Cross section and quantum yields as given below in the source code       =*
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
      REAL x(kdata), y(kdata)

* local

      REAL yg(kw)
      REAL qy1, qy2
      INTEGER iw

      INTEGER mabs, myld, xsvers(3), qyvers(3)
      CHARACTER(llog) :: logmsg(nlog)
      CHARACTER(lspc) :: spc


***************** hydroxyacetone photolysis *************************

* Setting photolysis index and cross section/quantum yield options
      j = j+1
      jlabel(j) = 'CH2(OH)COCH3 -> CH3CO + CH2(OH)'
      j = j+1
      jlabel(j) = 'CH2(OH)COCH3 -> CH2(OH)CO + CH3'

* cross section options
      spc = "CH2(OH)COCH3"
      xsvers = (/2, 1, 2/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "from IUPAC"
      logmsg(2) = "from JPL"
      logmsg(3) = "from Messaadia et al"
      CALL set_option(3,spc,"abs",logmsg,xsvers,mabs)

* quantum yield options
      qyvers = (/1, 3, 2/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "from Orlando et al. 2006"
      logmsg(2) = "from IUPAC"
      logmsg(3) = "scaled externally"
      CALL set_option(3,spc,"yld",logmsg,qyvers,myld)


* cross sections

      if (mabs.eq.1) then
         OPEN(UNIT=kin,FILE='DATAJ1/ABS/Hydroxyacetone.abs',
     $        STATUS='old')
         DO i = 1, 8
            READ(kin,*)
         ENDDO
         n = 101
         DO i = 1, n
            READ(kin,*) x(i), y(i)
         ENDDO
         CLOSE(kin)

      ELSEIF(mabs .eq. 2) then
         OPEN(UNIT=kin,FILE='DATAJ1/ABS/Hydroxyacetone_jpl11.abs',
     $        STATUS='old')
         DO i = 1, 2
            READ(kin,*)
         ENDDO
         n = 96
         DO i = 1, n
            READ(kin,*) x(i), y(i)
            y(i) = y(i) * 1.e-20
         ENDDO
         CLOSE(kin)

      ELSEIF(mabs .eq. 3) then
         OPEN(UNIT=kin,FILE='DATAJ1/ABS/Hydroxyacetone_Mes12.abs',
     $        STATUS='old')
         DO i = 1, 4
            READ(kin,*)
         ENDDO
         n = 161
         DO i = 1, n
            READ(kin,*) x(i), y(i)
         ENDDO
         CLOSE(kin)

      ENDIF

      CALL interpol(x,y,kdata,n,nw,wl,jlabel(j),0.,0.,yg)

* Total quantum yield  = 0.65, fromm Orlando et al. 2006
* Assume equal branching for each of the two channels

      IF(myld==1) THEN
        qy1 = 0.325
        qy2 = 0.325
       ELSEIF(myld==2) THEN
        qy1 = 0.31
        qy2 = 0.29
       ELSEIF(myld==2) THEN
        qy1 = 1.0
        qy2 = 1.0
      ENDIF

      DO iw = 1, nw - 1
         DO i = 1, nz
            sq(j-1,i,iw) = yg(iw) * qy1
            sq(j,i,iw) = yg(iw) * qy2
         ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r113(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! HOBr

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for HOBr           =*
*=  HOBr -> OH + Br                                                          =*
*=  Cross section from JPL 2003                                              =*
*=  Quantum yield assumed unity as in JPL2003                                =*
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

      INTEGER j, i

* data arrays

* local

      REAL    sig!, qy
      INTEGER iw


***************** HOBr photolysis *************************

* from JPL2003

      j = j+1
      jlabel(j) = 'HOBr -> OH + Br'

C      qy = 1.
      DO iw = 1, nw - 1
         sig = 24.77 * EXP( -109.80*(LOG(284.01/wc(iw)))**2 ) +
     $         12.22 * exp(  -93.63*(LOG(350.57/wc(iw)))**2 ) +
     $         2.283 * exp(- 242.40*(LOG(457.38/wc(iw)))**2 )
         sig = sig * 1.e-20
         IF(wc(iw) .LT. 250. .OR. wc(iw) .GT. 550.) sig = 0.

         DO i = 1, nz
            sq(j,i,iw) = sig! * qy
         ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r114(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! BrO

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for BrO            =*
*=  BrO -> Br + O                                                            =*
*=  Cross section from JPL 2003                                              =*
*=  Quantum yield assumed unity as in JPL2003                                =*
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

* output weighting functions

      INTEGER j
      CHARACTER(lcl) jlabel(kj)
      REAL sq(kj,kz,kw)

* data arrays

      INTEGER n, i
      REAL x(20), y(20)

* local

      INTEGER iw
      REAL yg(kw), dum!, qy


***************** BrO photolysis *************************

* from JPL2003

      j = j+1
      jlabel(j) = 'BrO -> Br + O'

      OPEN(UNIT=kin,FILE='DATAJ1/ABS/BrO.jpl03',
     $     STATUS='old')
      DO i = 1, 14
         READ(kin,*)
      ENDDO
      n = 15
      DO i = 1, n
         READ(kin,*) x(i), dum, y(i)
         y(i) = y(i) * 1.e-20
      ENDDO
      n = n + 1
      x(n) = dum
      CLOSE(kin)

* use bin-to-bin interpolation

      CALL inter4(nw,wl,yg,n,x,y,1)

C      qy = 1.
      DO iw = 1, nw - 1
         DO i = 1, nz
            sq(j,i,iw) = yg(iw)! * qy
         ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r115(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! Br2

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for BrO            =*
*=  Br2 -> Br + Br                                                           =*
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

* output weighting functions

      INTEGER j
      CHARACTER(lcl) jlabel(kj)
      REAL sq(kj,kz,kw)

* data arrays

      INTEGER kdata
      PARAMETER(kdata=50)
      INTEGER n, i
      REAL x(kdata), y(kdata)

* local

      INTEGER iw
      REAL    yg(kw)!, qy


***************** Br2 photolysis *************************

      j = j + 1
      jlabel(j) = 'Br2 -> Br + Br'

* Absorption cross section from:
* Seery, D.J. and D. Britton, The continuous absorption spectra of chlorine,
* bromine, bromine chloride, iodine chloride, and iodine bromide, J. Phys.
* Chem. 68, p. 2263 (1964).

      OPEN(UNIT=kin,FILE='DATAJ1/ABS/Br2.abs',
     $     STATUS='old')

      DO i = 1, 6
         READ(kin,*)
      ENDDO
      n = 29
      DO i = 1, n
         READ(kin,*) x(i),  y(i)
      ENDDO
      CLOSE(kin)

      CALL interpol(x,y,kdata,n,nw,wl,jlabel(j),0.,0.,yg)

C      qy = 1.
      DO iw = 1, nw - 1
         DO i = 1, nz
            sq(j,i,iw) = yg(iw)! * qy
         ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r118(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! NO3-(aq)
*-----------------------------------------------------------------------------*
*= NO3-(aq) photolysis for snow simulations                                  =*
*=        a) NO3-(aq) + hv -> NO2 + O-                                       =*
*=        b) NO3-(aq) + hv -> NO2- + O(3P)                                   =*
*=  Cross section:                                                           =*
*=  Burley & Johnston, Geophys. Res. Lett., 19, 1359-1362 (1992)             =*
*=  Chu & Anastasio, J. Phys. Chem. A, 107, 9594-9602 (2003)                 =*
*=  Quantum yield:                                                           =*
*=  Warneck & Wurzinger, J. Phys. Chem., 92, 6278-6283 (1988)                =*
*=  Chu & Anastasio, J. Phys. Chem. A, 107, 9594-9602 (2003)                 =*
*-----------------------------------------------------------------------------*
*= NOTE: user may have to manually add these reactions to the end of the     =*
*= reaction list in file usrinp to include these reactions for a snow run:   =*
*= T 74 NO3-(aq) -> NO2 + O-                                                 =*
*= T 75 NO3-(aq) -> NO2- + O(3P)                                             =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE '../../params'

* input
      INTEGER nw,nz
      REAL wl(kw), wc(kw), tlev(kz), airden(kz)

* weighting functions
      CHARACTER(lcl) jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:
      INTEGER j

* data arrays
      INTEGER kdata
      PARAMETER(kdata=50)
      REAL x1(kdata),x2(kdata)
      REAL y1(kdata),y2(kdata)     ! y1 = 20'C, y2 = -20'C

* local
      REAL yg1(kw),yg2(kw)
      REAL qy1(kz), qy2
      INTEGER i, iw, n, n1, n2, iz

      INTEGER mabs, xsvers(3)
      CHARACTER(llog) :: logmsg(nlog)
      CHARACTER(lspc) :: spc


***************** nitrate photolysis *************************

*** NO3-(aq) quantum yields
* O- (OH and NO2) production
      j = j + 1
      jlabel(j) = 'NO3-(aq) -> NO2(aq) + O-'
      DO iz = 1, nz
*        qy1(iz) = 9.3e-3  ! Warneck & Wurzinger 1988
        qy1(iz) = exp(-2400./tlev(iz) + 3.6) ! Chu & Anastasio, 2003
      ENDDO

* O(3P) (NO2-(aq) ....> NO) production
      j = j + 1
      jlabel(j) = 'NO3-(aq) -> NO2-(aq) + O(3P)'
      qy2 = 1.1e-3  ! Warneck & Wurzinger '88


* cross section options
      spc = "Nitrate"
      xsvers = (/2, 2, 2/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "from Burley and Johnston"
      logmsg(2) = "from Chu and Anastasio 2003"
      CALL set_option(2,spc,"abs",logmsg,xsvers,mabs)


* cross sections

      if (mabs .eq. 1) then
*** NO3-(aq) cross sections from Burley & Johnston (header lines = 24,
* data lines = 19)
         OPEN(kin,FILE='DATAJ1/ABS/NO3-_BJ92.abs',STATUS='OLD')

         n = 24
         DO i = 1, n
            READ(kin,*)
         ENDDO
         n = 19
         DO i = 1, n
            READ(kin,*) x1(i), y1(i), y2(i)
            x2(i) = x1(i)
            y1(i)=y1(i)*1e-20
            y2(i)=y2(i)*1e-20
         ENDDO
         CLOSE(kin)
         n1 = n
         CALL interpol(x1,y1,kdata,n1,nw,wl,jlabel(j),0.,0.,yg1)

         n2 = n
         CALL interpol(x2,y2,kdata,n2,nw,wl,jlabel(j),0.,0.,yg2)

      elseif (mabs .eq. 2) then

*** NO3-(aq) cross sections from Chu and Anastasio 2003:
* convert from molar abs log10 to cm2 per molec

         OPEN(kin,FILE='DATAJ1/ABS/NO3-_CA03.abs',STATUS='OLD')
         n = 7
         do i = 1, n
            read(kin,*)
         enddo
         n = 43
         DO i = 1, n
            read(kin,*) x1(i), y1(i)
            y1(i) = y1(i) * 3.82e-21
         enddo
         n1 = n
         CALL interpol(x1,y1,kdata,n1,nw,wl,jlabel(j),0.,0.,yg2)

      endif


* combine

      DO iw = 1, nw-1
         DO iz = 1, nz
! Use yg1 and mabs 1 for 20'C
! comment out qy1/qy2 for qy = 1
            sq(j-1,iz,iw) = qy1(iz)*yg2(iw)
            sq(j,iz,iw) = qy2*yg2(iw)
         ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r119(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! MEK

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for                =*
*=    methylethylketone                                                      =*
*=  CH3COCH2CH3 photolysis:                                                  =*
*=           CH3COCH2CH3  -> CH3CO + CH2CH3                                  =*
*=                                                                           =*
*=  Cross section from Martinez et al. (1992)                                =*
*=                                                                           =*
*=  Quantum yield assumed 0.325 for each channel (J. Orlando, priv.comm.2003)=*
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
      PARAMETER(kdata=20000)

      INTEGER i, n
      REAL x(kdata), y(kdata)

* local

      REAL yg(kw), dum
      REAL qy, eta, ptorr
      INTEGER iw

      INTEGER myld, qyvers(3)
      CHARACTER(llog) :: logmsg(nlog)
      CHARACTER(lspc) :: spc


***************** MEK photolysis *************************

* Setting photolysis index and quantum yield options
      j = j+1
      jlabel(j) = 'CH3COCH2CH3 -> CH3CO + CH2CH3'

* quantum yield options
* myld 1:
* Raber, W.H. (1992) PhD Thesis, Johannes Gutenberg-Universitaet, Mainz, Germany.
* other channels assumed negligible (less than 10%).
* Total quantum yield  = 0.38 at 760 Torr.
* Stern-Volmer form given:  1/phi = 0.96 + 2.22e-3*P(torr)

* myld 2:
* IUPAC recommendation of 0.34 with quenching from above


      spc = "MEK"
      qyvers = (/1, 2, 2/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "from Raber, W.H. (1992)"
      logmsg(2) = "from IUPAC"
      CALL set_option(2,spc,"yld",logmsg,qyvers,myld)


* cross sections

      OPEN(UNIT=kin,FILE='DATAJ1/ABS/Martinez.abs',
     $     STATUS='old')
      DO i = 1, 4
         READ(kin,*)
      ENDDO
      n = 96
      DO i = 1, n
         READ(kin,*) x(i), dum, y(i)
         y(i) = y(i) * 1.e-20
      ENDDO
      CLOSE(kin)

      CALL interpol(x,y,kdata,n,nw,wl,jlabel(j),0.,0.,yg)

* quantum yields

*     compute local pressure in torr
      DO i = 1, nz
        IF(myld==1) THEN
          ptorr = 760.*airden(i)/2.69e19
          qy = 1./(0.96 + 2.22E-3*ptorr)
          qy = MIN(qy, 1.0)
         ELSEIF(myld==2) THEN
          eta = 760.*2.22E-3
          qy = 0.34 * (1.+eta)/(1.+eta*airden(i)/2.465e19)
        ENDIF
        DO iw = 1, nw-1
           sq(j,i,iw) = yg(iw) * qy
        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r120(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! PPN

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for PPN photolysis:    =*
*=       PPN + hv -> Products                                                =*
*=                                                                           =*
*=  Cross section: from JPL 2006 (originally from Harwood et al. 2003)       =*
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
      PARAMETER(kdata=100)

      INTEGER iw
      INTEGER i, n, n1, n2
      REAL x1(kdata), x2(kdata)
      REAL y1(kdata), y2(kdata)

* local

      REAL yg(kw), yg1(kw), yg2(kw), yg3(kw)
      real qyNO2, qyNO3
      REAL sig

      INTEGER myld, qyvers(3)
      CHARACTER(llog) :: logmsg(nlog)
      CHARACTER(lspc) :: spc


***************** PPN photolysis *************************

* Setting photolysis index and quantum yield options
      j = j+1
      jlabel(j) = 'CH3CH2CO(OONO2) -> CH3CH2CO(OO) + NO2'
      j = j+1
      jlabel(j) = 'CH3CH2CO(OONO2) -> CH3CH2CO(O) + NO3'

* quantum yield options
      spc = "PPN"
      qyvers = (/1, 2, 2/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "from Harwood et al (JPL)"
      logmsg(2) = "from Calvert et al. 2011 book"
      CALL set_option(2,spc,"yld",logmsg,qyvers,myld)


* cross sections
*      from JPL 2011 (originally from Harwood et al. 2003)

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

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)
      CALL interpol(x2,y2,kdata,n2,nw,wl,jlabel(j),0.,0.,yg1)


* quantum yields

      IF(myld == 1) THEN
        qyNO2 = 0.61
        qyNO3 = 0.39
       ELSEIF(myld == 2) THEN

        OPEN(UNIT=kin,FILE='DATAJ1/RONO2/PPN_calv.yld',STATUS='OLD')
        DO i = 1, 8
          READ(kin,*)
        ENDDO
        n = 63
        DO i = 1, n
          READ(kin,*) x1(i), y1(i), y2(i)
          x2(i) = x1(i)
        ENDDO
        n1 = n
        n2 = n
        CLOSE(kin)

        CALL interpol(x1,y1,kdata,n1,nw,wl,jlabel(j),x1(1),x1(n1),yg2)
        CALL interpol(x2,y2,kdata,n2,nw,wl,jlabel(j),x2(1),x2(n2),yg3)

      ENDIF

      DO iw = 1, nw-1
        IF(myld==2) THEN
          qyNO2 = yg2(iw)
          qyNO3 = yg3(iw)
        ENDIF
        DO i = 1, nz

          sig = yg(iw) * EXP(yg1(iw)*(tlev(i)-298.))

          sq(j-1,i,iw) = qyNO2 * sig
          sq(j,i,iw)   = qyNO3 * sig

        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r121(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! CH2(OH)(OOH)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for CH2(OH)(OOH)       =*
*=  (hydroxy methyl hydroperoxide) photolysis:                               =*
*=       CH2(OH)(OOH) + hv -> CH2(OH)(O.) + OH                               =*
*=                                                                           =*
*=  Cross section: from JPL 2006 (originally from Bauerle and Moortgat 1999  =*
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
      PARAMETER(kdata=100)

      INTEGER iw
      INTEGER i, n
      REAL x1(kdata)
      REAL y1(kdata)

* local

      REAL yg(kw)
!      real qy
      REAL sig

      INTEGER mabs, xsvers(3)
      CHARACTER(llog) :: logmsg(nlog)
      CHARACTER(lspc) :: spc


***************** HOCH2OOH photolysis *************************

* Setting photolysis index and cross section options
      j = j+1
      jlabel(j) = 'HOCH2OOH -> HOCH2O. + OH'

      spc = "HOCH2OOH"
      xsvers = (/1, 2, 1/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "from JPL 2011"
      logmsg(2) = "from  IUPAC 2002"
      CALL set_option(2,spc,"abs",logmsg,xsvers,mabs)


* cross sections

      IF(mabs==1) THEN

      OPEN(UNIT=kin,FILE='DATAJ1/ABS/HOCH2OOH_jpl11.abs',STATUS='OLD')
      DO i = 1, 3
         READ(kin,*)
      ENDDO
      n = 32
      DO i = 1, n
         READ(kin,*) x1(i), y1(i)
         y1(i) = y1(i) * 1.E-20
      ENDDO
      CLOSE(kin)

      ELSEIF(mabs==2)THEN

      OPEN(UNIT=kin,FILE='DATAJ1/ABS/HOCH2OOH_iup.abs',STATUS='OLD')
      DO i = 1, 5
         READ(kin,*)
      ENDDO
      n = 33
      DO i = 1, n
         READ(kin,*) x1(i), y1(i)
         y1(i) = y1(i) * 1.E-20
      ENDDO
      CLOSE(kin)

      ENDIF

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)

* quantum yields assumed unity (also IUPAC recommended for > 290 nm
!      qy = 1.

      DO iw = 1, nw-1
        sig = yg(iw)
        DO i = 1, nz
          sq(j,i,iw)   = sig ! * qy
        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r122(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! CH2=CHCHO

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for CH2=CHCHO          =*
*=  (acrolein) photolysis:                                                   =*
*=       CH2=CHCHO + hv -> Products                                          =*
*=                                                                           =*
*=  Cross section: from JPL 2006 (originally from Magneron et al.            =*
*=  Quantum yield: P-dependent, JPL 2006 orig. from Gardner et al.           =*
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

      INTEGER iw
      INTEGER i, n
      REAL x1(kdata)
      REAL y1(kdata)

* local

      REAL yg(kw)
      real qy, qym1, qy1, qy2, qy3, qy4
      REAL sig

      INTEGER mabs, myld, mprod, xsvers(3), qyvers(3), prvers(3)
      CHARACTER(llog) :: logmsg(nlog)
      CHARACTER(lspc) :: spc


***************** acrolein photolysis *************************

* Setting photolysis index and cross section/quantum yield/pressure dep. options
      j = j+1
      jlabel(j) = 'CH2=CHCHO -> CH2=CH + CHO'
      j = j+1
      jlabel(j) = 'CH2=CHCHO -> CH2=CH2 + CO'
      j = j+1
      jlabel(j) = 'CH2=CHCHO -> CH3CHOO + CO'
      j = j+1
      jlabel(j) = 'CH2=CHCHO -> CH2=CHCO + H'

* cross section options
      spc = "Acrolein"
      xsvers = (/1, 1, 2/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "from JPL"
      logmsg(2) = "from Calvert et al. 2011 book"
      CALL set_option(2,spc,"abs",logmsg,xsvers,mabs)

* quantum yield options
      qyvers = (/1, 1, 2/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      ! log messages same as for cross sections
      CALL set_option(2,spc,"yld",logmsg,qyvers,myld)

* pressure-dependence options
      prvers = (/1, 2, 1/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "with estimated Criegie formation in main channel"
      logmsg(2) = "with estimated ethylene formation in main channel"
      CALL set_option(2,spc,"prd",logmsg,prvers,mprod)


* cross sections

      IF(mabs==1) THEN
        OPEN(UNIT=kin,FILE='DATAJ1/ABS/Acrolein.txt',STATUS='OLD')
        DO i = 1, 6
          READ(kin,*)
        ENDDO
        n = 55
        DO i = 1, n
          READ(kin,*) x1(i), y1(i)
          y1(i) = y1(i) * 1.E-20
        ENDDO
        CLOSE(kin)

       ELSEIF(mabs==2) THEN

        OPEN(UNIT=kin,FILE='DATAJ1/ABS/Acrolein_calv.dat',STATUS='OLD')
        DO i = 1, 6
          READ(kin,*)
        ENDDO
        n = 149
        DO i = 1, n
          READ(kin,*) x1(i), y1(i)
          y1(i) = y1(i) * 1.E-20
        ENDDO
        CLOSE(kin)
      ENDIF

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j-2),0.,0.,yg)

* quantum yields are pressure dependent between air number densities
* of 8e17 and 2.6e19, Gardner et al.:
* Product distribution from curve fits of data given in Calvert et al. 2011:

* product information taken from Calvert et al. 2011
* (fits to data of fractions of each channel)


      DO iw = 1, nw-1
        DO i = 1, nz
           if(airden(i) .gt. 2.6e19) then
              qy = 0.004
           elseif(airden(i) .gt. 8.e17 .and. airden(i) .lt. 2.6e19) then
             if(myld==1) then
               qym1 = 0.086 + 1.613e-17 * airden(i)
               qy = 0.004 + 1./qym1
              elseif(myld==2) then
               qym1 = -0.836+1.159e-17*airden(i)-2.166e-37*airden(i)**2
               qy = 1./qym1
             endif
           elseif(airden(i) .lt. 8.e17) then
             if(myld==1) then
              qym1 = 0.086 + 1.613e-17 * 8.e17
              qy = 0.004 + 1./qym1
              elseif(myld==2) then
               qym1 = -0.836+1.159e-17*airden(i)-2.166e-37*airden(i)**2
               qy = 1./qym1
             endif
           endif
* product distribution from Calvert et al. 2011:
* fits of direct qy given behind do not match 0.0065 at 1 atm
           IF(mprod == 1) THEN
             qy1 = qy*(-0.0173*(airden(i)/1.E19)**2
     &                +0.083*airden(i)/1.E19+0.0492)
             qy2 = qy*(-0.0111*(airden(i)/1.E19)**2
     &                +0.1014*airden(i)/1.E19+0.0532)
             qy3 = qy*(0.0518*(airden(i)/1.E19)**2
     &                -0.2675*airden(i)/1.E19+0.7953)
             qy4 = qy*(-0.0217*(airden(i)/1.E19)**2
     &                +0.0788*airden(i)/1.E19+0.1029)
            ELSEIF(mprod == 2) THEN
             qy1 = qy*(-0.0173*(airden(i)/1.E19)**2
     &                +0.083*airden(i)/1.E19+0.0492)
             qy2 = qy*(0.0407*(airden(i)/1.E19)**2
     &                -0.1661*airden(i)/1.E19+0.8485)
             qy3 = 0.0
             qy4 = qy*(-0.0217*(airden(i)/1.E19)**2
     &                +0.0788*airden(i)/1.E19+0.1029)
           ENDIF
           sig = yg(iw)
           sq(j-3,i,iw) = qy1 * sig
           sq(j-2,i,iw) = qy2 * sig
           sq(j-1,i,iw) = qy3 * sig
           sq(j  ,i,iw) = qy4 * sig
        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r123(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! CH3CO(OOH) + hv

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for peracetic acid     =*
*=      photolysis:                                                          =*
*=       CH3CO(OOH) + hv -> Products                                         =*
*=                                                                           =*
*=  Cross section: from JPL 2006 (originally from Orlando and Tyndall 2003   =*
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
      PARAMETER(kdata=100)

      INTEGER iw
      INTEGER i, n
      REAL x1(kdata)
      REAL y1(kdata)

* local

      REAL yg(kw)
C      real qy
      REAL sig


***************** peracetic acid photolysis *************************

      j = j+1
      jlabel(j) = 'CH3CO(OOH) -> Products'

* cross section from
*      JPL 2006 (originally from Orlando and Tyndall 2003)

      OPEN(UNIT=kin,FILE='DATAJ1/ABS/Peracetic_acid.txt',STATUS='OLD')
      DO i = 1, 6
         READ(kin,*)
      ENDDO
      n = 66
      DO i = 1, n
         READ(kin,*) x1(i), y1(i)
         y1(i) = y1(i) * 1.E-20
      ENDDO
      CLOSE(kin)

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)

* quantum yields assumed unity
!      qy = 1.

      DO iw = 1, nw-1
        DO i = 1, nz
          sig = yg(iw)
          sq(j,i,iw)   = sig! * qy
        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r124(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! (CH3)2NNO

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for dimethyl nitroso   =*
*=  amine photolysis:                                                        =*
*=       (CH3)2NNO + hv -> Products                                          =*
*=                                                                           =*
*=  Cross section: from Lindley 1978 (cited by Calvert et al. 2009)
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
      PARAMETER(kdata=150)

      INTEGER iw
      INTEGER i, n
      REAL x1(kdata)
      REAL y1(kdata)

* local

      REAL yg(kw)
C      REAL qy


***************** DMNA photolysis *************************

      j = j+1
      jlabel(j) = '(CH3)2NNO -> Products'

* cross section from
* Lindley (1978, PhD Thesis Ohio State U., Jack Calvert advisor), cited by Calvert et al. (2009).

      OPEN(UNIT=kin,FILE='DATAJ1/ABS/dmna.abs',STATUS='OLD')
      DO i = 1, 5
         READ(kin,*)
      ENDDO
      n = 132
      DO i = 1, n
         READ(kin,*) x1(i), y1(i)
         y1(i) = y1(i) * 1.E-19
      ENDDO
      CLOSE(kin)

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)

* quantum yields assumed unity
!      qy = 1.

      DO iw = 1, nw-1
        DO i = 1, nz
          sq(j,i,iw)   = yg(iw)! qy
        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r125(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! ClO

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for ClO photolysis     =*
*=       ClO + hv -> Cl + O                                                  =*
*=                                                                           =*
*=  Cross section: from Maric and Burrows 1999                               =*
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
      PARAMETER(kdata=500)

      INTEGER iw
      INTEGER i, n
      REAL x1(kdata)
      REAL y1(kdata)

* local

      REAL yg(kw)
      REAL qy1, qy2

      real tmp(12), x(kdata), y(kdata,12), tx
      real xdum
      integer m, nn, ii
      real ygt(kw, 12), yy
      INTEGER m1, m2


***************** ClO photolysis *************************

      j = j+1
      jlabel(j) = 'ClO -> Cl + O(1D)'
      j = j+1
      jlabel(j) = 'ClO -> Cl + O(3P)'

* cross section from
* Maric D. and J.P. Burrows, J. Quantitative Spectroscopy and
* Radiative Transfer 62, 345-369, 1999.  Data was downloaded from
* their web site on 15 September 2009.


      OPEN(UNIT=kin,FILE='DATAJ1/ABS/ClO_spectrum.prn',STATUS='OLD')
      DO i = 1, 2
         READ(kin,*)
      ENDDO
      nn = 453
      DO ii = 1, nn
         i = nn - ii + 1
         READ(kin,*) xdum, x(i), xdum, (y(i,m), m = 1, 12)
      ENDDO
      CLOSE(kin)

      DO m = 1, 12
         tmp(m) = 190. + 10.*FLOAT(m-1)
         IF(m .EQ. 1) tmp(m) = 180.

         DO i = 1, nn
            x1(i) = x(i)
            y1(i) = y(i,m)
         ENDDO
         n = nn
         CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)

         DO iw = 1, nw-1
            ygt(iw,m) = yg(iw)
         ENDDO

      ENDDO

      DO i = 1, nz

         tx = tlev(i)

* locate temperature indices for interpolation:
         m1 = 1 + INT((tx - 190.)/10.)
         m1 = MAX(1 ,m1)
         m1 = MIN(11,m1)
         m2 = m1 + 1

         DO iw = 1, nw-1

            yy = ygt(iw,m1) + (ygt(iw,m2)-ygt(iw,m1))
     $           *(tx-tmp(m1))/(tmp(m2)-tmp(m1))

* threshold for O(1D) productionis 263.4 nm:

            if(wc(iw) .lt. 263.4) then
               qy1 = 1.
            else
               qy1 = 0.
            endif
            qy2 = 1. - qy1

            sq(j-1,i,iw) = qy1 * yy
            sq(j,i,iw)   = qy2 * yy

         ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r126(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! ClNO2

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for nitryl chloride    =*
*=       ClNO2 -> Cl + NO2                                                   =*
*=                                                                           =*
*=  Cross section: from JPL 2006                                             =*
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
      PARAMETER(kdata=150)

      INTEGER iw
      INTEGER i, n
      REAL x1(kdata)
      REAL y1(kdata)

* local

      REAL yg(kw)
C      REAL qy

      INTEGER mabs, xsvers(3)
      CHARACTER(llog) :: logmsg(nlog)
      CHARACTER(lspc) :: spc

***************** ClNO2 photolysis *************************

* Setting photolysis index and cross section options
      j = j+1
      jlabel(j) = 'ClNO2 -> Cl + NO2'

* cross section options
* mabs = 1:   JPL 2006, same as JPL-2011
* mabs = 2:   IUPAC 2007
      spc = "ClNO2"
      xsvers = (/1, 2, 2/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "from JPL"
      logmsg(2) = "from IUPAC"
      CALL set_option(2,spc,"abs",logmsg,xsvers,mabs)


* cross sections

      mabs = 1
      if(mabs.eq.1) then

         OPEN(UNIT=kin,FILE='DATAJ1/ABS/ClNO2.abs',STATUS='OLD')
         DO i = 1, 2
            READ(kin,*)
         ENDDO
         n = 26
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
            y1(i) = y1(i) * 1.E-20
         ENDDO
         CLOSE(kin)

      elseif (mabs .eq. 2) then

         OPEN(UNIT=kin,FILE='DATAJ1/ABS/ClNO2_iupac.abs',STATUS='OLD')
         DO i = 1, 6
            READ(kin,*)
         ENDDO
         n = 17
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
         ENDDO
         CLOSE(kin)

      endif

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)


* quantum yields assumed unity
!      qy = 1.

      DO iw = 1, nw-1
        DO i = 1, nz
          sq(j,i,iw)   = yg(iw)! * qy
        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r127(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! BrNO

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for nitrosyl bromide   =*
*=       BrNO -> Br + NO                                                   =*
*=                                                                           =*
*=  Cross section: from JPL 2006                                             =*
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
      PARAMETER(kdata=150)

      INTEGER iw
      INTEGER i, n
      REAL x1(kdata)
      REAL y1(kdata)

* local

      REAL yg(kw)
C      REAL qy


***************** BrNO photolysis *************************

      j = j+1
      jlabel(j) = 'BrNO -> Br + NO'

* cross section from
* JPL 2006

      OPEN(UNIT=kin,FILE='DATAJ1/ABS/BrNO.abs',STATUS='OLD')
      DO i = 1, 3
         READ(kin,*)
      ENDDO
      n = 27
      DO i = 1, n
         READ(kin,*) x1(i), y1(i)
      ENDDO
      CLOSE(kin)

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)

* quantum yields assumed unity
!      qy = 1.

      DO iw = 1, nw-1
        DO i = 1, nz
          sq(j,i,iw)   = yg(iw)! * qy
        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r128(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! BrNO2

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for bromine nitrite    =*
*=       BrNO2 -> Br + NO2                                                   =*
*=                                                                           =*
*=  Cross section: from JPL 2006                                             =*
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
      PARAMETER(kdata=150)

      INTEGER iw
      INTEGER i, n
      REAL x1(kdata)
      REAL y1(kdata)

* local

      REAL yg(kw)
C      REAL qy


***************** BrNO2 photolysis *************************

      j = j+1
      jlabel(j) = 'BrNO2 -> Br + NO2'

* cross section from
* IUPAC (vol III) 2007

      OPEN(UNIT=kin,FILE='DATAJ1/ABS/BrNO2.abs',STATUS='OLD')
      DO i = 1, 6
         READ(kin,*)
      ENDDO
      n = 54
      DO i = 1, n
         READ(kin,*) x1(i), y1(i)
      ENDDO
      CLOSE(kin)

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)

* quantum yields assumed unity
!      qy = 1.

      DO iw = 1, nw-1
        DO i = 1, nz
          sq(j,i,iw)   = yg(iw)! * qy
        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r129(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! BrONO

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for bromine nitrite    =*
*=       BrONO -> Br + NO2                                                   =*
*=       BrONO -> BrO + NO                                                   =*
*=                                                                           =*
*=  Cross section: from IUPAC (vol.3)                                        =*
*=  Quantum yield: Assumed to be 0.5 for each channel                        =*
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

      REAL yg(kw)
      REAL qy1, qy2


***************** BrONO photolysis *************************

      j = j+1
      jlabel(j) = 'BrONO -> Br + NO2'
      j = j+1
      jlabel(j) = 'BrONO -> BrO + NO'

* cross section from
* IUPAC (vol III) 2007

      OPEN(UNIT=kin,FILE='DATAJ1/ABS/BrONO.abs',STATUS='OLD')
      DO i = 1, 8
         READ(kin,*)
      ENDDO
      n = 32
      DO i = 1, n
         READ(kin,*) x1(i), y1(i)
      ENDDO
      CLOSE(kin)

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)

* quantum yields assumed unity

      qy1 = 0.5
      qy2 = 0.5

      DO iw = 1, nw-1
        DO i = 1, nz
          sq(j-1,i,iw)   = qy1 * yg(iw)
          sq(j,i,iw)     = qy2 * yg(iw)
        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r130(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! HOCl

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for
*=       HOCl -> HO + Cl                                                     =*
*=  Cross section: from IUPAC (vol.3)                                        =*
*=  Quantum yield: Assumed to be 1                                           =*
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

      REAL yg(kw)
C      REAL qy

******************** HOCl photodissociation

      j = j + 1
      jlabel(j) = 'HOCl -> HO + Cl'

* cross section from
* IUPAC (vol III) 2007

      OPEN(UNIT=kin,FILE='DATAJ1/ABS/HOCl.abs',STATUS='OLD')
      DO i = 1, 7
         READ(kin,*)
      ENDDO
      n = 111
      DO i = 1, n
         READ(kin,*) x1(i), y1(i)
      ENDDO
      CLOSE(kin)

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)

* quantum yields assumed unity
!      qy = 1

      DO iw = 1, nw-1
        DO i = 1, nz
          sq(j,i,iw) = yg(iw)! * qy
        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r131(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! NOCl

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for
*=       NOCl -> NO + Cl                                                     =*
*=  Cross section: from IUPAC (vol.3)                                        =*
*=  Quantum yield: Assumed to be 1                                           =*
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
      INTEGER i, n, ii
      REAL x1(kdata), y1(kdata)
      integer nn
      REAL x223(kdata),x243(kdata),x263(kdata),x298(kdata),
     $     x323(kdata), x343(kdata)
      REAL y223(kdata),y243(kdata),y263(kdata),y298(kdata),
     $     y323(kdata), y343(kdata)

* local

      REAL yg223(kw),yg243(kw),yg263(kw),yg298(kw),
     $     yg323(kw), yg343(kw)
      REAL sig!, qy


***************** NOCl photolysis *************************

      j = j + 1
      jlabel(j) = 'NOCl -> NO + Cl'

* cross section from
* IUPAC (vol III) 2007

      OPEN(UNIT=kin,FILE='DATAJ1/ABS/NOCl.abs',STATUS='OLD')
      DO i = 1, 7
         READ(kin,*)
      ENDDO
      n = 80
      DO i = 1, n
         READ(kin,*) x1(i), y1(i)
         y223(i) = y1(i)
         y243(i) = y1(i)
         y263(i) = y1(i)
         y298(i) = y1(i)
         y323(i) = y1(i)
         y343(i) = y1(i)

         x223(i) = x1(i)
         x243(i) = x1(i)
         x263(i) = x1(i)
         x298(i) = x1(i)
         x323(i) = x1(i)
         x343(i) = x1(i)

      ENDDO
      READ(kin,*)
      n = 61
      do i = 1, n
         ii = i + 80
         read(kin,*) x1(ii), y223(ii), y243(ii), y263(ii),
     $        y298(ii), y323(ii), y343(ii)

         x223(ii) = x1(ii)
         x243(ii) = x1(ii)
         x263(ii) = x1(ii)
         x298(ii) = x1(ii)
         x323(ii) = x1(ii)
         x343(ii) = x1(ii)

      enddo
      n = ii
      CLOSE(kin)

      nn = n
      CALL interpol(x223,y223,kdata,nn,nw,wl,jlabel(j),0.,0.,yg223)
      nn = n
      CALL interpol(x243,y243,kdata,nn,nw,wl,jlabel(j),0.,0.,yg243)
      nn = n
      CALL interpol(x263,y263,kdata,nn,nw,wl,jlabel(j),0.,0.,yg263)
      nn = n
      CALL interpol(x298,y298,kdata,nn,nw,wl,jlabel(j),0.,0.,yg298)
      nn = n
      CALL interpol(x323,y323,kdata,nn,nw,wl,jlabel(j),0.,0.,yg323)
      nn = n
      CALL interpol(x343,y343,kdata,nn,nw,wl,jlabel(j),0.,0.,yg343)


* quantum yields assumed unity
!      qy = 1

      DO iw = 1, nw-1
        DO i = 1, nz

           if(tlev(i) .le. 223.) then
              sig = yg223(iw)

           elseif (tlev(i) .gt. 223. .and. tlev(i) .le. 243.) then
              sig = yg223(iw) +
     $             (yg243(iw) - yg223(iw))*(tlev(i) - 223.)/20.

           elseif (tlev(i) .gt. 243. .and. tlev(i) .le. 263.) then
              sig = yg243(iw) +
     $             (yg263(iw) - yg243(iw))*(tlev(i) - 243.)/20.

           elseif (tlev(i) .gt. 263. .and. tlev(i) .le. 298.) then
              sig = yg263(iw) +
     $             (yg298(iw) - yg263(iw))*(tlev(i) - 263.)/35.

           elseif (tlev(i) .gt. 298. .and. tlev(i) .le. 323.) then
              sig = yg298(iw) +
     $             (yg323(iw) - yg298(iw))*(tlev(i) - 298.)/25.

           elseif (tlev(i) .gt. 323. .and. tlev(i) .le. 343.) then
              sig = yg323(iw) +
     $             (yg343(iw) - yg323(iw))*(tlev(i) - 323.)/20.

           endif

           sq(j,i,iw) = sig! * qy

        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r132(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! OClO

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for
*=       OClO -> Products                                                    =*
*=  Cross section: from Wahner et al., J. Phys. Chem. 91, 2734, 1987         =*
*=  Quantum yield: Assumed to be 1                                           =*
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
      PARAMETER(kdata=2000)

      INTEGER iw
      INTEGER i
      integer nn, n204, n296, n378
      REAL x204(kdata),x296(kdata),x378(kdata)
      REAL y204(kdata),y296(kdata),y378(kdata)


* local

      REAL yg204(kw),yg296(kw),yg378(kw)
      REAL sig!, qy


***************** OClO photolysis *************************

      j = j + 1
      jlabel(j) = 'OClO -> Products'

* cross section from
*A. Wahner, G.S. tyndall, A.R. Ravishankara, J. Phys. Chem., 91, 2734, (1987).
*Supplementary Data, as quoted at:
*http://www.atmosphere.mpg.de/enid/26b4b5172008b02407b2e47f08de2fa1,0/Spectra/Introduction_1rr.html

      OPEN(UNIT=kin,FILE='DATAJ1/ABS/OClO.abs',STATUS='OLD')
      DO i = 1, 6
         READ(kin,*)
      ENDDO
      n204 = 1074-6
      DO i = 1, n204
         READ(kin,*) x204(i), y204(i)
      ENDDO

      READ(kin,*)
      n296 = 1067
      do i = 1, n296
         read(kin,*) x296(i), y296(i)
      enddo

      read(kin,*)
      n378 = 1068
      do i = 1, n378
         read(kin,*) x378(i), y378(i)
      enddo

      CLOSE(kin)

      nn = n204
      CALL interpol(x204,y204,kdata,nn,nw,wl,jlabel(j),0.,0.,yg204)
      nn = n296
      CALL interpol(x296,y296,kdata,nn,nw,wl,jlabel(j),0.,0.,yg296)
      nn = n378
      CALL interpol(x378,y378,kdata,nn,nw,wl,jlabel(j),0.,0.,yg378)


* quantum yields assumed unity
!      qy = 1

      DO iw = 1, nw-1
        DO i = 1, nz

           if(tlev(i) .le. 204.) then
              sig = yg204(iw)

           elseif (tlev(i) .gt. 204. .and. tlev(i) .le. 296.) then
              sig = yg204(iw) +
     $             (yg296(iw) - yg204(iw))*(tlev(i) - 204.)/92.

           elseif (tlev(i) .gt. 296. .and. tlev(i) .le. 378.) then
              sig = yg296(iw) +
     $             (yg378(iw) - yg296(iw))*(tlev(i) - 296.)/82.

           elseif (tlev(i) .gt. 378.) then
              sig = yg378(iw)
           endif

          sq(j,i,iw) = sig! * qy

        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r133(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! BrCl

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for                    =*
*=       BrCl -> Br + Cl                                                     =*
*=  Cross section: from Maric et al., J. Phtoochem Photobiol. A: Chem        =*
*=   83, 179-192, 1994.                                                      =*
*=  Quantum yield: Assumed to be 1                                           =*
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

      INTEGER iw
      INTEGER i, n
      REAL x(kdata), y(kdata)

* local

      REAL yg(kw)
C      REAL qy


***************** BrCl photolysis *************************

      j = j + 1
      jlabel(j) = 'BrCl -> Br + Cl'

* cross section from
* D. Maric, J.P. Burrows, and G.K. Moortgat, "A study of the UV-visible
* absorption spectra of Br2 and BrCl," J. Photochem. Photobiol. A: Chem.
* 83, 179-192 (1994).

      OPEN(UNIT=kin,FILE='DATAJ1/ABS/BrCl.abs',STATUS='OLD')
      DO i = 1, 9
         READ(kin,*)
      ENDDO
      n = 81
      DO i = 1, n
         READ(kin,*) x(i), y(i)
      ENDDO
      CLOSE(kin)

      CALL interpol(x,y,kdata,n,nw,wl,jlabel(j),0.,0.,yg)

* quantum yields assumed unity
!      qy = 1.

      DO iw = 1, nw-1
        DO i = 1, nz

          sq(j,i,iw) = yg(iw)! * qy

        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r134(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! CH3NO4

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for                    =*
*=       CH3(OONO2) -> CH3(OO) + NO2                                         =*
*=  Cross section: from
*= I. Bridier, R. Lesclaux, and B. Veyret, "Flash photolysis kinetic study
*= of the equilibrium CH3O2 + NO2 « CH3O2NO2," Chemical Physics Letters
*= 191, 259-263 (1992).
*=  Quantum yield: Assumed to be 1                                           =*
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

      INTEGER iw
      INTEGER i, n
      REAL x(kdata), y(kdata)

* local

      REAL yg(kw)
      REAL qy1, qy2

      INTEGER myld, qyvers(3)
      CHARACTER(llog) :: logmsg(nlog)
      CHARACTER(lspc) :: spc


***************** CH3(OONO2) photolysis *************************

* Setting photolysis index and quantum yield options
      j = j + 1
      jlabel(j) = 'CH3(OONO2) -> CH3(OO) + NO2'
      j = j + 1
      jlabel(j) = 'CH3(OONO2) -> CH3(O) + NO3'

* quantum yield options
      spc = "CH3(OONO2)"
      qyvers = (/1, 1, 2/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "assumed to unity" !also used for scaling externally
      logmsg(2) = "estimated same as HNO4"
      CALL set_option(2,spc,"yld",logmsg,qyvers,myld)


* cross sections

*= from
*= I. Bridier, R. Lesclaux, and B. Veyret, "Flash photolysis kinetic study
*= of the equilibrium CH3O2 + NO2 « CH3O2NO2," Chemical Physics Letters
*= 191, 259-263 (1992).

      OPEN(UNIT=kin,FILE='DATAJ1/ABS/CH3OONO2.abs',STATUS='OLD')
      DO i = 1, 9
         READ(kin,*)
      ENDDO
      n = 26
      DO i = 1, n
         READ(kin,*) x(i), y(i)
      ENDDO
      CLOSE(kin)

      CALL interpol(x,y,kdata,n,nw,wl,jlabel(j-1),0.,0.,yg)


* quantum yields assumed unity

      IF(myld == 1) THEN
        qy1 = 1.
        qy2 = 0.
       ELSEIF(myld==2) THEN
        qy1 = 0.59
        qy2 = 0.41
      ENDIF

* combine
      DO iw = 1, nw-1
        DO i = 1, nz

          sq(j-1,i,iw) = qy1 * yg(iw)
          sq(j  ,i,iw) = qy2 * yg(iw)

        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r135(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! C(CH3)3(ONO) + hv

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for t-butyl nitrite    =*
*=       C(CH3)3(ONO) -> C(CH3)3(O) + NO                                     =*
*=  Cross section: from                                                      =*
*=  V. McMillan, 1966, private communication to J.G. Calvert, J.N.Pitts, Jr.,=*
*=  Photochemistry, London, 1966, p. 455.                                    =*
*=  Quantum yield: Assumed to be 1                                           =*
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

      INTEGER iw
      INTEGER i, n
      REAL x(kdata), y(kdata)

* local

      REAL yg(kw)
C      REAL qy


***************** tert-butyl nitrite photolysis *************************

      j = j + 1
      jlabel(j) = 'C(CH3)3(ONO) -> C(CH3)3(O) + NO'

* cross section from
*=  V. McMillan, 1966, private communication to J.G. Calvert, J.N.Pitts, Jr.,
*=  Photochemistry, London, 1966, p. 455.

      OPEN(UNIT=kin,FILE='DATAJ1/ABS/t-butyl-nitrite.abs',STATUS='OLD')
      DO i = 1, 4
         READ(kin,*)
      ENDDO
      n = 96
      DO i = 1, n
         READ(kin,*) x(i), y(i)
      ENDDO
      CLOSE(kin)

      CALL interpol(x,y,kdata,n,nw,wl,jlabel(j),0.,0.,yg)


* quantum yields assumed unity
!      qy = 1.

      DO iw = 1, nw-1
        DO i = 1, nz

          sq(j,i,iw) = yg(iw)! * qy

        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r136(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! ClONO

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for ClONO              =*
*=        ClONO -> Cl + NO2                                                  =*
*=  cross section from IPUAC, orig from Molina and Molina (1977)             =*
*=  Quantum yield: Assumed to be 1                                           =*
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

      INTEGER iw
      INTEGER i, n
      REAL x(kdata), y(kdata)

* local

      REAL yg(kw)
C      REAL qy


***************** ClONO photolysis *************************

      j = j + 1
      jlabel(j) = 'ClONO -> Cl + NO2'

* cross section from JPL-2011
* Also published (with some minor differences) as:
* R. Atkinson, D.L. Baulch, R.A. Cox, J.N. Crowley, R.F. Hampson, R.G. Hynes, M.E. Jenkin, M.J. Rossi,
* and J. Troe, "Evaluated kinetic and photochemical data for atmospheric chemistry: Volume III - gas
* phase reactions of inorganic halogens", Atmos. Chem. Phys. 7, 981-1191 (2007).Comments:
* IUPAC (2005, 2007) recommendation:
* The preferred values of the absorption cross-sections at 231 K are the values reported by
* L.T. Molina and M.J. Molina, "Ultraviolet absorption spectrum of chlorine nitrite, ClONO,"
* Geophys. Res. Lett. 4, 83-86 (1977).

      OPEN(UNIT=kin,FILE='DATAJ1/ABS/ClONO_jpl11.abs',STATUS='OLD')
      DO i = 1, 3
         READ(kin,*)
      ENDDO
      n = 34
      DO i = 1, n
         READ(kin,*) x(i), y(i)
      ENDDO
      CLOSE(kin)

      CALL interpol(x,y,kdata,n,nw,wl,jlabel(j),0.,0.,yg)

* quantum yields assumed unity
!      qy = 1.

      DO iw = 1, nw-1
        DO i = 1, nz

          sq(j,i,iw) = yg(iw)! * qy

        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r137(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! HCl

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for HCl                =*
*=        HCl -> H + Cl                                                      =*
*=  cross section from JPL2011                                               =*
*=  Quantum yield: Assumed to be 1                                           =*
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
      PARAMETER(kdata=100)

      INTEGER iw
      INTEGER i, n
      REAL x(kdata), y(kdata)

* local

      REAL yg(kw)
C      REAL qy


***************** HCl photolysis *************************

      j = j + 1
      jlabel(j) = 'HCl -> H + Cl'

* cross section from JPL2011

      OPEN(UNIT=kin,FILE='DATAJ1/ABS/HCl_jpl11.abs',STATUS='OLD')
      DO i = 1, 3
         READ(kin,*)
      ENDDO
      n = 31
      DO i = 1, n
         READ(kin,*) x(i), y(i)
         y(i) = y(i) * 1.e-20
      ENDDO
      CLOSE(kin)

      CALL interpol(x,y,kdata,n,nw,wl,jlabel(j),0.,0.,yg)

* quantum yields assumed unity
!      qy = 1.

      DO iw = 1, nw-1
        DO i = 1, nz

          sq(j,i,iw) = yg(iw)! * qy

        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r138(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! CH3COOH

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for acetic acid        =*
*=        CH3COOH -> CH3 + COOH                                              =*
*=  cross section from JPL2011                                               =*
*=  Quantum yield: Assumed to be 0.55                                        =*
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
      PARAMETER(kdata=100)

      INTEGER iw
      INTEGER i, n
      REAL x(kdata), y(kdata)

* local

      REAL yg(kw)
      REAL qy


***************** acetic acid photolysis *************************

      j = j + 1
      jlabel(j) = 'CH3COOH -> CH3 + COOH'

* cross section from JPL2011

      OPEN(UNIT=kin,FILE='DATAJ1/ABS/CH3COOH_jpl11.abs',STATUS='OLD')
      DO i = 1, 2
         READ(kin,*)
      ENDDO
      n = 18
      DO i = 1, n
         READ(kin,*) x(i), y(i)
         y(i) = y(i) * 1.e-20
      ENDDO
      CLOSE(kin)

      CALL interpol(x,y,kdata,n,nw,wl,jlabel(j),0.,0.,yg)

* quantum yields assumed unity

      qy = 0.55


* combine
      DO iw = 1, nw-1
        DO i = 1, nz

          sq(j,i,iw) = qy * yg(iw)

        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r139(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! CH3OCl

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x(quantum yield) for methyl hypochlorite =*
*=        CH3OCl -> CH3O + Cl                                                =*
*=  cross section from JPL2011                                               =*
*=  Quantum yield: Assumed to be 1                                           =*
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
      PARAMETER(kdata=100)

      INTEGER iw
      INTEGER i, n
      REAL x(kdata), y(kdata)

* local

      REAL yg(kw)
C      REAL qy


***************** methyl hypochlorite photolysis *************************

      j = j + 1
      jlabel(j) = 'CH3OCl -> CH3O + Cl'


* cross section from JPL2011

      OPEN(UNIT=kin,FILE='DATAJ1/ABS/CH3OCl_jpl11.abs',STATUS='OLD')
      DO i = 1, 3
         READ(kin,*)
      ENDDO
      n = 83
      DO i = 1, n
         READ(kin,*) x(i), y(i)
         y(i) = y(i) * 1.e-20
      ENDDO
      CLOSE(kin)

      CALL interpol(x,y,kdata,n,nw,wl,jlabel(j),0.,0.,yg)

* quantum yields assumed unity
!      qy = 1.

      DO iw = 1, nw-1
        DO i = 1, nz

          sq(j,i,iw) = yg(iw)! * qy

        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r140(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! CHCl3

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for CHCl3 photolysis:  =*
*=      CHCL3 + hv -> Products                                               =*
*=  Cross section: from JPL 2011 recommendation                              =*
*=  Quantum yield: assumed to be unity                                       =*
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
      PARAMETER(kdata=100)

      REAL x1(kdata)
      REAL y1(kdata)

* local

      REAL yg(kw)
C      REAL qy
      INTEGER i, iw, n
      INTEGER iz
      REAL b0, b1, b2, b3, b4, tcoeff, sig
      REAL w1, w2, w3, w4, temp


***************** CHCl3 photolysis *************************

      j = j+1
      jlabel(j) = 'CHCl3 -> Products'


* cross sections

      OPEN(kin,FILE='DATAJ1/ABS/CHCl3_jpl11.abs',STATUS='OLD')
      DO i = 1, 3
         READ(kin,*)
      ENDDO
      n = 39
      DO i = 1, n
         READ(kin,*) x1(i), y1(i)
         y1(i) = y1(i) * 1E-20
      ENDDO
      CLOSE(kin)

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)

* compute temperature correction factors:

      b0 = 3.7973
      b1 = -7.0913e-2
      b2 = 4.9397e-4
      b3 = -1.5226e-6
      b4 = 1.7555e-9

*** quantum yield assumed to be unity
!      qy = 1.


      DO iw = 1, nw-1

* compute temperature correction coefficients:

         tcoeff = 0.
         IF(wc(iw) .GT. 190. .AND. wc(iw) .LT. 240.) THEN
            w1 = wc(iw)
            w2 = w1**2
            w3 = w1**3
            w4 = w1**4
            tcoeff = b0 + b1*w1 + b2*w2 + b3*w3 + b4*w4
         ENDIF

         DO iz = 1, nz
            temp = tlev(iz)
            temp = min(max(temp,210.),300.)
            sig = yg(iw) * 10.**(tcoeff*(temp-295.))
            sq(j,iz,iw) = sig! * qy
         ENDDO

      ENDDO

      END

*=============================================================================*

      SUBROUTINE r141(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! C2H5ONO2

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for CH3CH2ONO2     =*
*=  ethyl nitrate       photolysis:                                          =*
*=           CH3CH2ONO2 + hv -> CH3CH2O + NO2                                =*
*=                                                                           =*
*=  Absorption cross sections as given below in source code                  =*
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

      INTEGER i, n, n1, n2
      REAL x1(kdata), y1(kdata)
      REAL x2(kdata), y2(kdata)

* local

      REAL dum
      REAL yg1(kw), yg2(kw)
!     REAL sig, qy
      INTEGER iw

      INTEGER mabs, xsvers(3)
      CHARACTER(llog) :: logmsg(nlog)
      CHARACTER(lspc) :: spc


***************** ethyl nitrate photolysis *************************

* Setting photolysis index and cross section options
      j = j+1
      jlabel(j) = 'C2H5ONO2 -> C2H5O + NO2'

      spc = "C2H5ONO2"
      xsvers = (/1, 2, 1/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "from from Talukdar et al"
      logmsg(2) = "from IUPAC "
      CALL set_option(2,spc,"abs",logmsg,xsvers,mabs)


* cross sections

      IF(mabs==1) THEN

      OPEN(UNIT=kin,FILE='DATAJ1/RONO2/RONO2_talukdar.abs',
     $     STATUS='old')
      DO i = 1, 10
         READ(kin,*)
      ENDDO
      n1 = 0
      n2 = 0
      DO i = 1, 63
         READ(kin,*) x1(i), dum, dum, y1(i), y2(i)
         if (y1(i) .gt. 0.) n1 = n1 + 1
         if (y2(i) .gt. 0.) n2 = n2 + 1
         x2(i) = x1(i)
         y1(i) = y1(i) * 1.e-20
         y2(i) = y2(i) * 1.e-3
      ENDDO
      CLOSE(kin)

      ELSEIF(mabs==2) THEN

      OPEN(UNIT=kin,FILE='DATAJ1/RONO2/C2H5NO3_iup.abs',
     $     STATUS='old')
      DO i = 1, 5
         READ(kin,*)
      ENDDO
      n = 32
      n1 = n
      n2 = n
      DO i = 1, n
         READ(kin,*) x1(i), y1(i), y2(i)
         x2(i) = x1(i)
         y1(i) = y1(i) * 1.e-20
         y2(i) = y2(i) * 1.e-3
      ENDDO
      CLOSE(kin)

      ENDIF

      CALL interpol(x1,y1,kdata,n1,nw,wl,jlabel(j),0.,0.,yg1)
      CALL interpol(x2,y2,kdata,n2,nw,wl,jlabel(j),0.,y2(n2),yg2)

* quantum yield  = 1

      DO iw = 1, nw - 1
         DO i = 1, nz
            sq(j,i,iw) = yg1(iw)*exp(yg2(iw)*(tlev(i)-298.)) ! sig * qy
         ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r142(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! n-C3H7ONO2

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for n-C3H7ONO2         =*
*=  photolysis:                                                              =*
*=          n-C3H7ONO2 + hv -> C3H7O + NO2                                     =*
*=                                                                           =*
*=  Cross section:  IUPAC 2006 (Atkinson et al., ACP, 6, 3625-4055, 2006)    =*
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
      PARAMETER (kdata = 200)

      INTEGER i, n
      INTEGER iw
      REAL x1(kdata)
      REAL y1(kdata)

* local

      REAL yg1(kw)
!     REAL qy


***************** n-propyl nitrate photolysis *************************

      j = j + 1
      jlabel(j) = 'n-C3H7ONO2 -> C3H7O + NO2'

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


* quantum yield = 1

      DO iw = 1, nw - 1
         DO i = 1, nz
            sq(j,i,iw) = yg1(iw) ! * qy
         ENDDO
      ENDDO

      RETURN
      END

*=============================================================================*

      SUBROUTINE r143(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! 1-C4H9ONO2

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for 1-C4H9ONO2         =*
*=  photolysis:                                                              =*
*=          1-C4H9ONO2 + hv -> 1-C4H9O + NO2                                 =*
*=                                                                           =*
*=  Cross section:  IUPAC 2006 (Atkinson et al., ACP, 6, 3625-4055, 2006)    =*
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
      PARAMETER (kdata = 200)

      INTEGER i, n
      INTEGER iw
      REAL x1(kdata)
      REAL y1(kdata)

* local

      REAL yg1(kw)
!     REAL qy


***************** n-butyl nitrate photolysis *************************

      j = j + 1
      jlabel(j) = '1-C4H9ONO2 -> 1-C4H9O + NO2'

      OPEN(UNIT=kin,FILE='DATAJ1/RONO2/1C4H9ONO2_iup2006.abs',
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

* quantum yield = 1

      DO iw = 1, nw - 1
         DO i = 1, nz
            sq(j,i,iw) = yg1(iw) ! * qy
         ENDDO
      ENDDO

      RETURN
      END

*=============================================================================*

      SUBROUTINE r144(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! 2-C4H9ONO2

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for 2-C4H9ONO2         =*
*=  photolysis:                                                              =*
*=          2-C4H9ONO2 + hv -> 2-C4H9O + NO2                                 =*
*=                                                                           =*
*=  Cross section:  IUPAC 2006 (Atkinson et al., ACP, 6, 3625-4055, 2006)    =*
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
      PARAMETER (kdata = 200)

      INTEGER i, n
      INTEGER iw
      REAL x1(kdata)
      REAL y1(kdata)

* local

      REAL yg1(kw)
!     REAL qy


***************** sec-butyl nitrate photolysis *************************

      j = j + 1
      jlabel(j) = '2-C4H9ONO2 -> 2-C4H9O + NO2'

* 1:  IUPAC 2006

      OPEN(UNIT=kin,FILE='DATAJ1/RONO2/2C4H9ONO2_iup2006.abs',
     $     STATUS='old')
      DO i = 1, 3
         READ(kin,*)
      ENDDO
      n = 15
      DO i = 1, n
         READ(kin,*) x1(i), y1(i)
         y1(i) = y1(i) * 1.e-20
      ENDDO
      CLOSE (kin)

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg1)

* quantum yield = 1

      DO iw = 1, nw - 1
         DO i = 1, nz
            sq(j,i,iw) = yg1(iw) ! * qy
         ENDDO
      ENDDO

      RETURN
      END

*=============================================================================*

      SUBROUTINE r145(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! perfluoro n-iodo propane (H24)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for                    =*
*=     perfluoro n-iodo propane (H24)                                        =*
*=  cross section from JPL2011                                               =*
*=  Quantum yield: Assumed to be 0.55                                        =*
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
      PARAMETER(kdata=100)

      INTEGER iw
      INTEGER i, n
      REAL x(kdata), y(kdata)

* local

      REAL yg(kw)
C      REAL qy


***************** perfluoro 1-iodopropane photolysis *************************

      j = j + 1
      jlabel(j) = 'perfluoro 1-iodopropane -> products'


* cross section from JPL2011

      OPEN(UNIT=kin,FILE='DATAJ1/ABS/PF-n-iodopropane.abs',STATUS='OLD')
      DO i = 1, 2
         READ(kin,*)
      ENDDO
      n = 16
      DO i = 1, n
         READ(kin,*) x(i), y(i)
         y(i) = y(i) * 1.e-20
      ENDDO
      CLOSE(kin)

      CALL interpol(x,y,kdata,n,nw,wl,jlabel(j),0.,0.,yg)


* quantum yields assumed unity
!      qy = 1.

      DO iw = 1, nw-1
        DO i = 1, nz

          sq(j,i,iw) = yg(iw)! * qy

        ENDDO
      ENDDO

      END

*=============================================================================*
