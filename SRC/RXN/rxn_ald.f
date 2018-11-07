*= This file contains the following subroutines, related to reading/loading
*= the product (cross section) x (quantum yield) for photo-reactions of
*= aldehydes in MCM-GECKO, which were not yet present in TUV5.2:
*=
*=     ma01 through ma19

*=============================================================================*

      SUBROUTINE ma01(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! n-C3H7CHO

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section)x(quantum yield) for n-C3H7CHO photolysis =*
*=          n-C3H7CHO + hv -> products (Norish type I + II)                  =*
*=                                                                           =*
*=  Cross section and quantum yield from IUPAC                               =*
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

      REAL yg(kw),yg1(kw),ygoh(kw)
      REAL qy1, qy2, qy0, eta
      REAL sig,sigoh
      INTEGER iw

      INTEGER myld, qyvers(3)
      CHARACTER(llog) :: logmsg(nlog)
      CHARACTER(lspc) :: spc


***************** n-C3H7CHO photolysis *************************

* Setting photolysis index and quantum yield options
* Norish type I:
      j = j+1
      jlabel(j) = 'n-C3H7CHO -> n-C3H7 + CHO'
* Norish type II:
      j = j+1
      jlabel(j) = 'n-C3H7CHO -> C2H4 + CH2CHOH'
* test OH substituted butyraldehyde:
      j = j+1
      jlabel(j) = 'ALD4OH -> NI products'
      j = j+1
      jlabel(j) = 'ALD4OH -> NII products'
* estimates in compounds with polyfunctional chromophores
      j = j+1
      jlabel(j) = 'C4nALDpoly'
      j = j+1
      jlabel(j) = 'C4nALDOHpoly'

      spc = "n-C3H7CHO"
      qyvers = (/1, 2, 2/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "from Moortgat 1999"
      logmsg(2) = "from IUPAC"
      CALL set_option(2,spc,"yld",logmsg,qyvers,myld)


* absorption cross sections

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
      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)

      n = n1
      CALL interpol(x1oh,y1oh,kdata,n,nw,wl,jlabel(j),0.,0.,ygoh)


* quantum yields

      IF(myld==1) THEN
        qy1 = 0.19
        qy2 = 0.06
       ELSEIF(myld==2) THEN
        qy0 = 0.21
        qy2 = 0.1
      ENDIF


* pressure dependence

      OPEN(UNIT=kin,FILE='DATAJ1/MCMext/ALD/nC3H7CHO_Chen02.yld',
     $     STATUS='old')
      do i = 1, 4
         read(kin,*)
      enddo

      n = 11
      DO i = 1, n
         READ(kin,*) x1(i), y1(i)
      ENDDO
      CLOSE(kin)

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j-1),qy0,qy0,yg1)


* combine:

      DO iw = 1, nw - 1
         sig = yg(iw)
         eta = MAX(0.,yg1(iw)/qy0-1.)
         sigoh = ygoh(iw)
         DO i = 1, nz
           qy1 = qy0*(1.+eta)/(1.+eta*airden(i)/2.465e19)
           sq(j-5,i,iw) = sig * qy1
           sq(j-4,i,iw) = sig * qy2
           IF(swOH==1) THEN
             sq(j-3,i,iw) = sigoh * qy1
             sq(j-2,i,iw) = sigoh * qy2
           ELSEIF(swOH==2) THEN
             IF(qy1+qy2 > 0.) THEN
               sq(j-3,i,iw) = sigoh * qyoh*qy1/(qy1+qy2)
               sq(j-2,i,iw) = sigoh * qyoh*qy2/(qy1+qy2)
              ELSE
               sq(j-3,i,iw) = 0.
               sq(j-2,i,iw) = 0.
             ENDIF
           ENDIF
           sq(j-1,i,iw) = sig
           sq(j,  i,iw) = sigoh
         ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE ma02(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! i-C3H7CHO

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section)x(quantum yield) for i-C3H7CHO photolysis =*
*=          i-C3H7CHO + hv -> i-C3H7 + CHO                                   =*
*=                                                                           =*
*=  Cross section and quantum yield as given below in source code            =*
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
      REAL x1(kdata), x2(kdata)
      REAL y1(kdata), y2(kdata)

* local

      REAL yg(kw), yg1(kw), dum
      REAL qy
      REAL sig
      INTEGER iw

      INTEGER mabs, myld, xsvers(3), qyvers(3)
      CHARACTER(llog) :: logmsg(nlog)
      CHARACTER(lspc) :: spc


***************** i-C3H7CHO photolysis *************************

* Setting photolysis index and cross section/quantum yield options
* only Norish type I for i-C3H7CHO!
      j = j+1
      jlabel(j) = 'i-C3H7CHO -> i-C3H7 + CHO'

* cross section options
      spc = "i-C3H7CHO"
      xsvers = (/1, 2, 2/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "from Martinez et al. 1992"
      logmsg(2) = "from IUPAC"
      CALL set_option(2,spc,"abs",logmsg,xsvers,mabs)

* quantum yield options
      qyvers = (/1, 4, 4/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "from Desai et al. 1986"
      logmsg(2) = "from Calvert book, opt. 1"
      logmsg(3) = "from Calvert book, opt. 2"
      logmsg(4) = "from IUPAC"
      CALL set_option(4,spc,"yld",logmsg,qyvers,myld)


* cross sections

      IF(mabs==1) THEN

        OPEN(UNIT=kin,FILE='DATAJ1/MCMext/ALD/i_butyraldehyde_R.prn',
     $       STATUS='old')
        do i = 1, 3
           read(kin,*)
        enddo

        n = 101
        DO i = 1, n
           READ(kin,*) x1(i), y1(i), dum
        ENDDO
        CLOSE(kin)

       ELSEIF(mabs==2) THEN

        OPEN(UNIT=kin,FILE='DATAJ1/MCMext/ALD/iC3H7CHO_iup.abs',
     $       STATUS='old')
        do i = 1, 5
           read(kin,*)
        enddo

        n = 121
        DO i = 1, n
           READ(kin,*) x1(i), y1(i)
           y1(i) = y1(i) * 1E-20
        ENDDO
        CLOSE(kin)
      ENDIF

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)


* quantum yields
! myld 3 needs adjustments of zero-pressure qy!

      IF(myld==1) THEN

        OPEN(UNIT=kin,FILE='DATAJ1/MCMext/ALD/i_butyraldehyde_R.prn',
     $       STATUS='old')
        do i = 1, 3
           read(kin,*)
        enddo

        n = 101
        DO i = 1, n
           READ(kin,*) x1(i), dum, y1(i)
        ENDDO
        CLOSE(kin)

       ELSEIF(myld==2 .OR. myld==3) THEN

        OPEN(UNIT=kin,FILE='DATAJ1/MCMext/ALD/iC3H7CHO_calv.yld',
     $       STATUS='old')
        do i = 1, 6
           read(kin,*)
        enddo

        n = 73
        DO i = 1, n
           READ(kin,*) x1(i), dum, y1(i),y2(i)
           x2(i) = x1(i)
        ENDDO
        CLOSE(kin)

       ELSEIF(myld==4) THEN

        OPEN(UNIT=kin,FILE='DATAJ1/MCMext/ALD/iC3H7CHO_iup.yld',
     $       STATUS='old')
        do i = 1, 7
           read(kin,*)
        enddo

        n = 11
        DO i = 1, n
           READ(kin,*) x1(i), y1(i)
        ENDDO
        CLOSE(kin)

      ENDIF

      IF(myld==3)THEN
        CALL interpol(x2,y2,kdata,n,nw,wl,jlabel(j),0.,0.,yg1)
       ELSE
        CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg1)
      ENDIF

*combine
      DO iw = 1, nw - 1
         sig = yg (iw)
         qy  = yg1(iw)
         DO i = 1, nz
            sq(j  ,i,iw) = sig * qy
         ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE ma03(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! n-C4H9CHO

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for n-C4H9CHO      =*
*=  photolysis:                                                              =*
*=         n-C4H9CHO + hv -> Norish type I + II products                     =*
*=                                                                           =*
*=  Cross section:  see options below in the source code                     =*
*=  Quantum yield:  see options below in the source code                     =*
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

      REAL yg(kw), yg1(kw), yg2(kw), ygoh(kw), dum
      REAL qy1, qy2, qy3, qy0, eta, ptorr
      REAL sig, sigoh
      INTEGER iw

      INTEGER mabs, myld, xsvers(3), qyvers(3)
      CHARACTER(llog) :: logmsg(nlog)
      CHARACTER(lspc) :: spc


***************** n-C4H9CHO photolysis *************************


* Setting photolysis index and cross section/quantum yield options


      j = j+1
      jlabel(j) = 'n-C4H9CHO -> C4H9 +  CHO'
      j = j+1
      jlabel(j) = 'n-C4H9CHO -> CH3CH=CH2 +  CH2=CHOH'
      j = j+1
      jlabel(j) = 'n-C4H9CHO -> 2-methylcyclobutanol'
      j = j+1
      jlabel(j) = 'C5nALDOH -> NI products'
      j = j+1
      jlabel(j) = 'C5nALDOH -> NII products'
      j = j+1
      jlabel(j) = 'C5nALDOH -> cycl. products'
* estimates in compounds with polyfunctional chromophores
      j = j+1
      jlabel(j) = 'C5nALDpoly'
      j = j+1
      jlabel(j) = 'C5nALDOHpoly'

* cross section options
      spc = "n-C4H9CHO"
      xsvers = (/1, 2, 2/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "from Zhu et al. 99"
      logmsg(2) = "from Tadic et al. 2001"
      CALL set_option(2,spc,"abs",logmsg,xsvers,mabs)

* quantum yield options
      qyvers = (/1, 4, 3/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "from Zhu et al. 99"
      logmsg(2) = "from Calvert et al. 2011 book"
      logmsg(3) = "based on Tadic et al. 2001"
      logmsg(4) =
     &        "based on Tadic et al. 2001 (branching scaled externally)"
      CALL set_option(4,spc,"yld",logmsg,qyvers,myld)


* cross sections

      IF(mabs==1) THEN
        OPEN(UNIT=kin,FILE='DATAJ1/MCMext/ALD/n_pentanal_rad',
     $       STATUS='old')

        do i = 1, 3
          read(kin,*)
        enddo

        n = 23
        DO i = 1, n
          READ(kin,*) x1(i), y1(i), dum
        ENDDO
        CLOSE(kin)

       ELSEIF(mabs==2) THEN

        OPEN(UNIT=kin,FILE='DATAJ1/MCMext/ALD/Tadic.abs',
     $       STATUS='old')

        do i = 1, 5
          read(kin,*)
        enddo

        n = 121
        DO i = 1, n
          READ(kin,*) x1(i), dum, y1(i)
        ENDDO
        CLOSE(kin)
      ENDIF

      x1oh(:) = 0.
      x1oh(:) = x1(:) - 10.
      y1oh(:) = y1(:)
      n1 = n
      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j-8),0.,0.,yg)
      n = n1
      CALL interpol(x1oh,y1oh,kdata,n,nw,wl,jlabel(j-2),0.,0.,ygoh)


* quantum yields

      IF(myld==1) THEN
        OPEN(UNIT=kin,FILE='DATAJ1/MCMext/ALD/n_pentanal_rad',
     $       STATUS='old')

        do i = 1, 3
          read(kin,*)
        enddo

        n = 23
        DO i = 1, n
          READ(kin,*) x1(i), dum, y1(i)
        ENDDO
        CLOSE(kin)

       ELSEIF(myld==2) THEN

        OPEN(UNIT=kin,FILE='DATAJ1/MCMext/ALD/pentanal_calv.yld',
     $       STATUS='old')

        do i = 1, 4
          read(kin,*)
        enddo

        n = 66
        DO i = 1, n
          READ(kin,*) x1(i), y1(i)
        ENDDO
        CLOSE(kin)
      ENDIF

      IF(myld<=2) THEN
        CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg1)
      ENDIF


* pressure dependency

      IF(myld<=2) THEN
        OPEN(UNIT=kin,FILE='DATAJ1/MCMext/ALD/valeraldehyde_C+Z98.yld',
     $       STATUS='old')

        do i = 1, 4
          read(kin,*)
        enddo

        n = 11
        DO i = 1, n
          READ(kin,*) x1(i), y1(i)
        ENDDO
        CLOSE(kin)

        CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j-3),0.,0.,yg2)
      ENDIF


* combine:

      DO iw = 1, nw - 1
        sig = yg(iw)
        sigoh = ygoh(iw)
        IF(myld<=2) THEN
          qy0 = yg1(iw)
          qy2 = 0.19
          qy3 = 0.00
          eta = MAX(0.,yg2(iw)/qy0-1.)
        ENDIF
        DO i = 1, nz
          IF(myld<=2) THEN
            qy1 = qy0*(1.+eta)/(1.+eta*airden(i)/2.465e19)
          ELSEIF(myld>=3) THEN
            ptorr = 760.*airden(i)/2.55e19
            qy3 = 1.0/(2.44 + 7.771E-4*ptorr)
            IF(myld==3) THEN
              qy1 = 0.13*qy3
              qy2 = 0.70*qy3
              qy3 = 0.17*qy3
            ELSE
              qy1 = qy3
              qy2 = qy3
            ENDIF
          ENDIF
          sq(j-7,i,iw) = sig * qy1
          sq(j-6,i,iw) = sig * qy2
          sq(j-5,i,iw) = sig * qy3
          IF(swOH==1) THEN
            sq(j-4,i,iw) = sigoh * qy1
            sq(j-3,i,iw) = sigoh * qy2
            sq(j-2,i,iw) = sigoh * qy3
          ELSEIF(swOH==2) THEN
            IF(qy1+qy2+qy3 > 0.) THEN
              sq(j-4,i,iw) = sigoh * qyoh*qy1/(qy1+qy2+qy3)
              sq(j-3,i,iw) = sigoh * qyoh*qy2/(qy1+qy2+qy3)
              sq(j-2,i,iw) = sigoh * qyoh*qy3/(qy1+qy2+qy3)
             ELSE
              sq(j-4,i,iw) = 0.
              sq(j-3,i,iw) = 0.
              sq(j-2,i,iw) = 0.
            ENDIF
          ENDIF
          sq(j-1,i,iw) = sig
          sq(j,  i,iw) = sigoh
        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE ma04(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! i-C4H9CHO

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for i-C4H9CHO      =*
*=  photolysis:                                                              =*
*=         i-C4H9CHO + hv -> Norish type I + II products                     =*
*=                                                                           =*
*=  Cross section:  see options below in the source code                     =*
*=  Quantum yield:  see options below in the source code                     =*
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

      INTEGER i, n, n1, n2
      REAL x1(kdata), x2(kdata)
      REAL y1(kdata), y2(kdata)

* local

      REAL yg(kw), yg1(kw), yg2(kw), dum
      REAL qy1, qy2, qy3, qy0, eta, ptorr
      REAL sig
      INTEGER iw

      INTEGER mabs, myld, xsvers(3), qyvers(3)
      CHARACTER(llog) :: logmsg(nlog)
      CHARACTER(lspc) :: spc


***************** i-C4H9CHO photolysis *************************


* Setting photolysis index and cross section/quantum yield options
      j = j+1
      jlabel(j) = 'i-C4H9CHO -> C4H9 + CHO'
      j = j+1
      jlabel(j) = 'i-C4H9CHO -> CH3CH=CH2 + CH2=CHOH'
      j = j+1
      jlabel(j) = 'i-C4H9CHO -> 3-methyl-1-cyclobutanol'

* cross section options
      spc = "i-C4H9CHO"
      xsvers = (/1, 5, 1/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      WRITE(logmsg(1),"(2A)") "from Zhu et al. 99 extented to 390nm ",
     &       "as given in GECKO-A database at 0 pressure"
      WRITE(logmsg(2),"(2A)") "from Zhu et al. 99 extented to 390nm ",
     &       "as given by Calvert et al. 2008 at 0 pressure"
      logmsg(3) = "from Lanza et al. 2008 extented to 390nm"
      WRITE(logmsg(4),"(2A)") "mean of Lanza et al. 2008 and ",
     &       "Calvert et al. 2011 data extented to 390nm"
      WRITE(logmsg(5),"(2A)") "from Calvert et al. 2011 data ",
     &       "extented with a tail between 330 and 350nm"
      CALL set_option(5,spc,"abs",logmsg,xsvers,mabs)

* quantum yield options
      qyvers = (/1, 2, 3/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      WRITE(logmsg(1),"(2A)") "from Zhu et al. 99 extented to 390nm ",
     &       "as given in GECKO-A database at 0 pressure"
      WRITE(logmsg(2),"(2A)") "with estimated quenching ",
     &       "above 310nm according to data in Calvert et al. (2008)"
      logmsg(3) = "from Tadic et al. (2012)"
      CALL set_option(3,spc,"yld",logmsg,qyvers,myld)


      IF(mabs==1) THEN
        OPEN(UNIT=kin,FILE='DATAJ1/MCMext/ALD/i_pentanal_rad',
     $       STATUS='old')
        do i = 1, 2
          read(kin,*)
        enddo

        n = 23
        DO i = 1, n
          READ(kin,*) x1(i), y1(i), dum
        ENDDO
        CLOSE(kin)

       ELSEIF(mabs==2) THEN

        OPEN(UNIT=kin,
     $       FILE='DATAJ1/MCMext/ALD/isovaleraldehyde_calv.dat',
     $       STATUS='old')
        do i = 1, 8
          read(kin,*)
        enddo

        n = 61
        DO i = 1, n
          READ(kin,*) x1(i), y1(i), dum, dum
          y1(i) = y1(i) * 1.E-20
        ENDDO
        CLOSE(kin)

       ELSEIF(mabs==3) THEN

        OPEN(UNIT=kin,
     $       FILE='DATAJ1/MCMext/ALD/isovaleraldehyde_Lan08ext.abs',
     $       STATUS='old')
        do i = 1, 4
          read(kin,*)
        enddo

        n = 88
        DO i = 1, n
          READ(kin,*) x1(i), y1(i)
        ENDDO
        CLOSE(kin)

       ELSEIF(mabs==4) THEN

        OPEN(UNIT=kin,
     $       FILE='DATAJ1/MCMext/ALD/isovaleraldehyde_avrgext.abs',
     $       STATUS='old')
        do i = 1, 9
          read(kin,*)
        enddo

        n = 63
        DO i = 1, n
          READ(kin,*) x1(i), y1(i)
          y1(i) = y1(i) * 1.E-20
        ENDDO
        CLOSE(kin)

       ELSEIF(mabs==5) THEN

        OPEN(UNIT=kin,
     $       FILE='DATAJ1/MCMext/ALD/isovaleraldehyde_calvext.abs',
     $       STATUS='old')
        do i = 1, 6
          read(kin,*)
        enddo

        n = 71
        DO i = 1, n
          READ(kin,*) x1(i), y1(i)
          y1(i) = y1(i) * 1.E-20
        ENDDO
        CLOSE(kin)
      ENDIF

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j-1),0.,0.,yg)


* quantum yields

      IF(myld==1) THEN
        OPEN(UNIT=kin,FILE='DATAJ1/MCMext/ALD/i_pentanal_rad',
     $       STATUS='old')
        do i = 1, 2
          read(kin,*)
        enddo

        n = 23
        DO i = 1, n
          READ(kin,*) x1(i), dum, y1(i)
        ENDDO
        CLOSE(kin)

        CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,y1(n),yg1)

       ELSEIF(myld==2) THEN

! yg1: qy1 at 1 atm; yg2: qy1 at 0 atm
        OPEN(UNIT=kin,
     $       FILE='DATAJ1/MCMext/ALD/isovaleraldehyde_qu.dat',
     $       STATUS='old')
        do i = 1, 11
          read(kin,*)
        enddo

        n  = 61
        n1 = n
        n2 = n
        DO i = 1, n
          READ(kin,*) x1(i), dum, y2(i), y1(i)
          x2(i) = x1(i)
        ENDDO
        CLOSE(kin)

        CALL interpol(x1,y1,kdata,n1,nw,wl,jlabel(j-1),0.,0.,yg1)
        CALL interpol(x2,y2,kdata,n2,nw,wl,jlabel(j-1),0.,y2(n),yg2)
      ENDIF

* set minor pressure-independent pathway
      IF(myld < 3) THEN
        qy2 = 0.13
        qy3 = 0.
      ENDIF

* combine:

      DO iw = 1, nw - 1
         sig = yg(iw)
         IF(myld==1) THEN
           qy1 = yg1(iw)
          ELSEIF(myld==2) THEN
           qy0 = yg1(iw)
           IF(qy0>=pzero) THEN
             eta = MAX(0.,yg2(iw)/qy0-1.)
            ELSE
             eta = 0.
           ENDIF
          ELSEIF(myld==3) THEN
           ptorr = 760.*airden(i)/2.55e19
           qy0 = 1.0/(1.96 + 1.798E-3*ptorr)
           qy1 = 0.4*qy0
           qy2 = 0.38*qy0
           qy3 = 0.22*qy0
         ENDIF
         DO i = 1, nz
           IF(myld==2) qy1 = qy0*(1.+eta)/(1.+eta*airden(i)/2.465e19)
           sq(j-2,i,iw) = sig * qy1
           sq(j-1,i,iw) = sig * qy2
           sq(j,  i,iw) = sig * qy3
         ENDDO
      ENDDO

      END

* ============================================================================*

      SUBROUTINE ma05(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! sec-C4H9CHO

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for sec-C4H9CHO    =*
*=  photolysis:                                                              =*
*=         sec-C4H9CHO + hv -> Norish type I + II products                   =*
*=                                                                           =*
*=  Cross section:  estimated same as i-C4H9CHO                              =*
*=  Quantum yield:  Calvert et al. (2011) book                               =*
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
      REAL qy0, qy1, qy2, ptorr
      REAL sig
      INTEGER iw

      INTEGER mabs, myld, xsvers(3), qyvers(3)
      CHARACTER(llog) :: logmsg(nlog)
      CHARACTER(lspc) :: spc


***************** sec-C4H9CHO photolysis *************************


* Setting photolysis index and cross section/quantum yield options
      j = j+1
      jlabel(j) = 'sec-C4H9CHO -> C4H9 + CHO'
      j = j+1
      jlabel(j) = 'sec-C4H9CHO -> CH3CH=CHOH + CH2=CH2'

* cross section options
      spc = "sec-C4H9CHO"
      xsvers = (/1, 2, 2/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "estimated same as i-C4H9CHO from GECKO-A database"
      WRITE(logmsg(2),"(2A)") "estimated same as i-C4H9CHO ",
     &   "based on Calvert et al. (2011)"
      CALL set_option(2,spc,"abs",logmsg,xsvers,mabs)

* quantum yield options
      qyvers = (/2, 2, 1/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "from Calvert et al. (2011) book"
      logmsg(2) = "from Calvert et al. (2011) with estimated quenching"
      CALL set_option(2,spc,"yld",logmsg,qyvers,myld)

* Absorption cross sections
* estimated same as i-C4H9CHO

      IF(mabs==1) THEN
        OPEN(UNIT=kin,FILE='DATAJ1/MCMext/ALD/i_pentanal_rad',
     $       STATUS='old')
        do i = 1, 2
          read(kin,*)
        enddo

        n = 23
        DO i = 1, n
          READ(kin,*) x1(i), y1(i), dum
        ENDDO
        CLOSE(kin)

       ELSEIF(mabs==2) THEN

        OPEN(UNIT=kin,
     $       FILE='DATAJ1/MCMext/ALD/isovaleraldehyde_avrgext.abs',
     $       STATUS='old')
        do i = 1, 9
          read(kin,*)
        enddo

        n = 63
        DO i = 1, n
          READ(kin,*) x1(i), y1(i)
          y1(i) = y1(i) * 1.E-20
        ENDDO
        CLOSE(kin)
      ENDIF

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j-1),0.,0.,yg)

* quantum yields
* from Calvert et al. (2011) book
* with possible option for estimation of pressure dependence

* zero pressure qy:
      qy0 = 0.55
      qy2 = 0.15


* combine:
      DO iw = 1, nw - 1
         sig = yg(iw)
         DO i = 1, nz
           IF(myld==2) THEN
             ptorr = 760.*airden(i)/2.55e19
             qy1 = MAX(0.,MIN(1.,1./(1./qy0+1.061e-3*ptorr)))
            ELSE
             qy1 = qy0
           ENDIF
           sq(j-1,i,iw) = sig * qy1
           sq(j  ,i,iw) = sig * qy2
         ENDDO
      ENDDO

      END

* ============================================================================*

      SUBROUTINE ma06(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! t-C4H9CHO/tALD

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for t-C4H9CHO      =*
*=  photolysis:                                                              =*
*=         t-C4H9CHO + hv -> HCO. +  t-C4H9.                                 =*
*=  Also provides values for generalised t-aldehydes with possible OH-subst. =*
*=                                                                           =*
*=  Cross section:  see options below in the source code                     =*
*=  Quantum yield:  see options below in the source code                     =*
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

      INTEGER i, n, n1, n2
      REAL x1(kdata), x2(kdata)
      REAL y1(kdata), y2(kdata)

* local

      REAL yg(kw), yg1(kw), yg2(kw), dum
      REAL qy, qy0, eta
      REAL sig
      INTEGER iw

      INTEGER mabs, myld, xsvers(3), qyvers(3)
      CHARACTER(llog) :: logmsg(nlog)
      CHARACTER(lspc) :: spc


***************** t-C4H9CHO photolysis *************************


* Setting photolysis index and cross section/quantum yield options
      j = j+1
      jlabel(j) = 't-C4H9CHO -> C4H9 + CHO'
      j = j+1
      jlabel(j) = 'tALD -> products'

* cross section options
      spc = "t-C4H9CHO"
      xsvers = (/1, 1, 2/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      WRITE(logmsg(1),"(2A)") "from Zhu et al. 99 extented to 390nm ",
     &                "as given in GECKO-A database at 0 pressure"
      WRITE(logmsg(2),"(2A)") "from Zhu et al. 99 extented to 335nm ",
     &                "as given by Calvert et al. 2008 at 0 pressure"
      CALL set_option(2,spc,"abs",logmsg,xsvers,mabs)

* quantum yield options
      qyvers = (/1, 4, 4/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      WRITE(logmsg(1),"(2A)") "from Zhu et al. 99 extented to 390nm ",
     &       "as given in GECKO-A database at 0 pressure"
      WRITE(logmsg(2),"(2A)") "from Zhu et al. 99 extented to 340nm ",
     &       "as given by Calvert et al. (2008) at 0 pressure"
      logmsg(3) = "from Zhu et al. 99 with estimated quenching"
      logmsg(4) =
     &         "from Zhu et al. 99 with estimated quenching above 310nm"
      CALL set_option(4,spc,"yld",logmsg,qyvers,myld)


      IF(mabs==1) THEN
        OPEN(UNIT=kin,FILE='DATAJ1/MCMext/ALD/t_pentanal_rad',
     $       STATUS='old')
        do i = 1, 3
          read(kin,*)
        enddo

        n = 23
        DO i = 1, n
          READ(kin,*) x1(i), y1(i), dum
        ENDDO
        CLOSE(kin)

       ELSEIF(mabs==2) THEN

        OPEN(UNIT=kin,
     $       FILE='DATAJ1/MCMext/ALD/pivaldehyde_calv.dat',
     $       STATUS='old')
        do i = 1, 8
          read(kin,*)
        enddo

        n = 61
        DO i = 1, n
          READ(kin,*) x1(i), y1(i), dum, dum
          y1(i) = y1(i) * 1.E-20
        ENDDO
        CLOSE(kin)
      ENDIF

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j-2),0.,0.,yg)


* quantum yields

      IF(myld==1) THEN
        OPEN(UNIT=kin,FILE='DATAJ1/MCMext/ALD/t_pentanal_rad',
     $       STATUS='old')
        do i = 1, 3
          read(kin,*)
        enddo

        n = 23
        DO i = 1, n
          READ(kin,*) x1(i), dum, y1(i)
        ENDDO
        CLOSE(kin)

          CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j-1),0.,y1(n),yg1)

       ELSE

        OPEN(UNIT=kin,
     $       FILE='DATAJ1/MCMext/ALD/pivaldehyde_calv.dat',
     $       STATUS='old')
        do i = 1, 8
          read(kin,*)
        enddo

        n = 61
        n1 = n
        n2 = n
        DO i = 1, n
          READ(kin,*) x1(i), dum, y2(i), y1(i)
          x2(i) = x1(i)
        ENDDO
        CLOSE(kin)

        CALL interpol(x1,y1,kdata,n1,nw,wl,jlabel(j),0.,0.,yg1)
        CALL interpol(x2,y2,kdata,n2,nw,wl,jlabel(j),0.,0.,yg2)

      ENDIF

! General quantum yields scaled externally
C      qyoh = 0.75
C      qyg1 = 0.4
C      qyg2 = 0.13

* combine:

      DO iw = 1, nw - 1
        sig = yg(iw)
        IF(myld==1) THEN
          qy = yg1(iw)
         ELSEIF(myld==2) THEN
          qy = yg2(iw)
         ELSEIF(myld>=3) THEN
          qy0 = yg1(iw)
          IF(qy0>=pzero) THEN
            eta = MAX(0.,yg2(iw)/qy0-1.)
           ELSE
            eta = 0.
          ENDIF
        ENDIF
        DO i = 1, nz
          IF(myld==3 .OR. (myld==4 .AND. wc(iw)>=310.)) THEN
            qy = qy0*(1.+eta)/(1.+eta*airden(i)/2.465e19)
          ENDIF
          sq(j-1,i,iw) = sig * qy
          sq(j  ,i,iw) = sig
        ENDDO
      ENDDO

      END

* ============================================================================*

      SUBROUTINE ma07(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! n-C5H11CHO

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for n-C5H11CHO     =*
*=  photolysis:                                                              =*
*=         n-C5H11CHO + hv -> Norish type I + II products                    =*
*=                                                                           =*
*=  Cross section:  see options below in the source code                     =*
*=  Quantum yield:  see options below in the source code                     =*
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

      REAL yg(kw), yg1(kw), yg2(kw), ygoh(kw), dum
      REAL qy1, qy2, qy3, qy0, eta, ptorr
      REAL sig, sigoh
      INTEGER iw

      INTEGER mabs, myld, xsvers(3), qyvers(3)
      CHARACTER(llog) :: logmsg(nlog)
      CHARACTER(lspc) :: spc


***************** n-C5H11CHO photolysis *************************


* Setting photolysis index and cross section/quantum yield options
      j = j+1
      jlabel(j) = 'n-C5H11CHO -> C5H11 +  CHO'
      j = j+1
      jlabel(j) = 'n-C5H11CHO -> C2H5CH=CH2 +  CH2=CHOH'
      j = j+1
      jlabel(j) = 'n-C5H11CHO -> 2-ethylcyclobutanol'
      j = j+1
      jlabel(j) = 'C6nALDOH -> NI products'
      j = j+1
      jlabel(j) = 'C6nALDOH -> NII products'
      j = j+1
      jlabel(j) = 'C6nALDOH -> cycl. product'
* estimates in compounds with polyfunctional chromophores
      j = j+1
      jlabel(j) = 'C6nALDpoly'
      j = j+1
      jlabel(j) = 'C6nALDOHpoly'

* cross section options
      spc = "n-C5H11CHO"
      xsvers = (/1, 1, 1/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "from Plagens et al. 1998"
      logmsg(2) = "from O'Connor et al. 2006"
      logmsg(3) = "from Jimenez et al. 2007"
      CALL set_option(4,spc,"abs",logmsg,xsvers,mabs)

* quantum yield options
      qyvers = (/1, 4, 3/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "from GECKO-A database"
      logmsg(2) = "from Calvert et al. 2011 book"
      logmsg(3) = "from Tadic et al. 2001"
      logmsg(4) = "from Tadic et al. 2001 (branching scaled externally)"
      CALL set_option(4,spc,"yld",logmsg,qyvers,myld)


      IF(mabs==1) THEN
        OPEN(UNIT=kin,FILE='DATAJ1/MCMext/ALD/hexanal_rad',
     $       STATUS='old')
        do i = 1, 2
          read(kin,*)
        enddo

        n = 89
        DO i = 1, n
          READ(kin,*) x1(i), y1(i), dum
        ENDDO
        CLOSE(kin)

       ELSEIF(mabs==2) THEN

        OPEN(UNIT=kin,FILE='DATAJ1/MCMext/ALD/nHexanal_OCon06.abs',
     $       STATUS='old')

        do i = 1, 5
          read(kin,*)
        enddo

        n = 64
        DO i = 1, n
          READ(kin,*) x1(i), y1(i)
        ENDDO
        CLOSE(kin)

       ELSEIF(mabs==3) THEN

        OPEN(UNIT=kin,FILE='DATAJ1/MCMext/ALD/nHexanal_Jim07.abs',
     $       STATUS='old')

        do i = 1, 5
          read(kin,*)
        enddo

        n = 73
        DO i = 1, n
          READ(kin,*) x1(i), y1(i), dum
        ENDDO
        CLOSE(kin)
      ENDIF

      x1oh(:) = 0.
      x1oh(:) = x1(:) - 10.
      y1oh(:) = y1(:)
      n1 = n
      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)
      n = n1
      CALL interpol(x1oh,y1oh,kdata,n,nw,wl,jlabel(j),0.,0.,ygoh)


* quantum yields


      IF(myld==2) THEN

        OPEN(UNIT=kin,FILE='DATAJ1/MCMext/ALD/nHexanal_calv11.yld',
     $       STATUS='old')

        do i = 1, 5
          read(kin,*)
        enddo

        n = 68
        DO i = 1, n
          READ(kin,*) x1(i), y1(i)
        ENDDO
        CLOSE(kin)

        CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg1)
      ENDIF


* pressure dependency

      OPEN(UNIT=kin,FILE='DATAJ1/MCMext/ALD/hexanal_T+Z04.yld',
     $     STATUS='old')

      do i = 1, 4
        read(kin,*)
      enddo

      n = 11
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
      ENDDO
      CLOSE(kin)

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg2)


* combine:

      DO iw = 1, nw - 1
        sig = yg(iw)
        sigoh = ygoh(iw)
        IF(myld==1) THEN
          qy0 = 0.075
          qy2 = 0.175
          qy3 = 0.
         ELSEIF(myld==2) THEN
          qy0 = yg1(iw)
          qy2 = 0.29
          qy3 = 0.
        ENDIF
        IF(myld<=2) eta = MAX(0.,yg2(iw)/qy0-1.)
        DO i = 1, nz
          IF(myld<=2) THEN
            qy1 = qy0*(1.+eta)/(1.+eta*airden(i)/2.465e19)
          ELSEIF(myld>=3) THEN
            ptorr = 760.*airden(i)/2.55e19
            qy3 = 1.0/(2.26 + 4.758E-4*ptorr)
            IF(myld==3) THEN
              qy1 = 0.21*qy3
              qy2 = 0.61*qy3
              qy3 = 0.18*qy3
            ELSE
              qy1 = qy3
              qy2 = qy3
            ENDIF
          ENDIF
          sq(j-7,i,iw) = sig * qy1
          sq(j-6,i,iw) = sig * qy2
          sq(j-5,i,iw) = sig * qy3
          IF(swOH==1) THEN
            sq(j-4,i,iw) = sigoh * qy1
            sq(j-3,i,iw) = sigoh * qy2
            sq(j-2,i,iw) = sigoh * qy3
          ELSEIF(swOH==2) THEN
            IF(qy1+qy2+qy3 > 0.) THEN
              sq(j-4,i,iw) = sigoh * qyoh*qy1/(qy1+qy2+qy3)
              sq(j-3,i,iw) = sigoh * qyoh*qy2/(qy1+qy2+qy3)
              sq(j-2,i,iw) = sigoh * qyoh*qy3/(qy1+qy2+qy3)
             ELSE
              sq(j-4,i,iw) = 0.
              sq(j-3,i,iw) = 0.
              sq(j-2,i,iw) = 0.
            ENDIF
          ENDIF
          sq(j-1,i,iw) = sig
          sq(j,  i,iw) = sigoh
        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE ma08(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! n-C6H13CHO

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for n-C6H13CHO     =*
*=  photolysis:                                                              =*
*=         n-C6H13CHO + hv -> Norish type I + II products                    =*
*=                                                                           =*
*=  Cross section:  see options below in the source code                     =*
*=  Quantum yield:  see options below in the source code                     =*
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

      REAL yg(kw), yg1(kw), yg2(kw), ygoh(kw)
      REAL qy1, qy2, qy3, qy0, eta, ptorr
      REAL sig, sigoh
      INTEGER iw

      INTEGER myld, qyvers(3)
      CHARACTER(llog) :: logmsg(nlog)
      CHARACTER(lspc) :: spc

***************** n-C6H13CHO photolysis *************************


* Setting photolysis index and quantum yield options
      j = j+1
      jlabel(j) = 'n-C6H13CHO -> C6H13 + CHO'
      j = j+1
      jlabel(j) = 'n-C6H13CHO -> C3H7CH=CH2 + CH2=CHOH'
      j = j+1
      jlabel(j) = 'n-C6H13CHO -> 2-propylcyclobutanol'
      j = j+1
      jlabel(j) = 'C7nALDOH -> NI products'
      j = j+1
      jlabel(j) = 'C7nALDOH -> NII products'
      j = j+1
      jlabel(j) = 'C7nALDOH -> cycl. product'
* estimates in compounds with polyfunctional chromophores
      j = j+1
      jlabel(j) = 'C7nALDpoly'
      j = j+1
      jlabel(j) = 'C7nALDOHpoly'


      spc = "n-C6H13CHO"
      qyvers = (/1, 3, 2/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = 'from Calvert et al. 2011 book'
      logmsg(2) = 'from Tadic et al. 2001'
      logmsg(3) = 'from Tadic et al. 2001 (branching scaled externally)'
      CALL set_option(3,spc,"yld",logmsg,qyvers,myld)


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

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j-7),0.,0.,yg)
      n = n1
      CALL interpol(x1oh,y1oh,kdata,n,nw,wl,jlabel(j-6),0.,0.,ygoh)


* quantum yields
      IF(myld==1) THEN
        OPEN(UNIT=kin,FILE='DATAJ1/MCMext/ALD/nHeptanal_calv.yld',
     $       STATUS='old')

        do i = 1, 7
          read(kin,*)
        enddo

        n = 71
        DO i = 1, n
          READ(kin,*) x1(i), y1(i)
        ENDDO
        CLOSE(kin)

        CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j-8),0.,0.,yg1)
      ENDIF


* pressure dependency

        OPEN(UNIT=kin,FILE='DATAJ1/MCMext/ALD/heptanal_T+Z04.yld',
     $       STATUS='old')

        do i = 1, 4
          read(kin,*)
        enddo

        n = 11
        DO i = 1, n
          READ(kin,*) x1(i), y1(i)
        ENDDO
        CLOSE(kin)

        CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j-7),0.,0.,yg2)


* combine:

      DO iw = 1, nw - 1
        sig = yg(iw)
        sigoh = ygoh(iw)
        IF(myld==1) THEN
          qy0 = yg1(iw)
          qy2 = 0.12
          qy3 = 0.00
          eta = MAX(0.,yg2(iw)/qy0-1.)
        ENDIF
        DO i = 1, nz
          IF(myld==1) THEN
            qy1 = qy0*(1.+eta)/(1.+eta*airden(i)/2.465e19)
          ELSEIF(myld>=2) THEN
            ptorr = 760.*airden(i)/2.55e19
            qy3 = 1.0/(2.408 + 1.169E-3*ptorr)
            IF(myld==2) THEN
              qy1 = 0.10*qy3
              qy2 = 0.72*qy3
              qy3 = 0.18*qy3
            ELSE
              qy1 = qy3
              qy2 = qy3
            ENDIF
          ENDIF
          sq(j-7,i,iw) = sig * qy1
          sq(j-6,i,iw) = sig * qy2
          sq(j-5,i,iw) = sig * qy3
          IF(swOH==1) THEN
            sq(j-4,i,iw) = sigoh * qy1
            sq(j-3,i,iw) = sigoh * qy2
            sq(j-2,i,iw) = sigoh * qy3
          ELSEIF(swOH==2) THEN
            IF(qy1+qy2+qy3 > 0.) THEN
              sq(j-4,i,iw) = sigoh * qyoh*qy1/(qy1+qy2+qy3)
              sq(j-3,i,iw) = sigoh * qyoh*qy2/(qy1+qy2+qy3)
              sq(j-2,i,iw) = sigoh * qyoh*qy3/(qy1+qy2+qy3)
             ELSE
              sq(j-4,i,iw) = 0.
              sq(j-3,i,iw) = 0.
              sq(j-2,i,iw) = 0.
            ENDIF
          ENDIF
          sq(j-1,i,iw) = sig
          sq(j,  i,iw) = sigoh
        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE ma09(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! Glycidaldehyde

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for glycidaldehyde =*
*=  photolysis:                                                              =*
*=         Glycidaldehyde + hv ->  products                                  =*
*=                                                                           =*
*=  Cross section and quantum yield:  taken from Calvert et al. 2011         =*
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
      INTEGER ::         myld, qyvers(3)
      CHARACTER(llog) :: logmsg(nlog)
      CHARACTER(lspc) :: spc

***************** Glycidaldehyde photolysis *************************


* Ratio for both channels estimated from ratio of channel I + II
* from glycolaldehyde photolysis with overall quantum yield = 0.6
* (see Calvert et al. 2011 for total qy)
! 89%:
      j = j+1
      jlabel(j) = 'Glycidaldehyde -> oxyranyl radical + CHO'
! 11%:
      j = j+1
      jlabel(j) = 'Glycidaldehyde -> oxyrane + CO'


      spc = "Glycidaldehyde"
      qyvers = (/1, 2, 1/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      WRITE(logmsg(1),"(2A)") 'from Calvert et al. 2011 book ',
     &                        'with branching ratio from glycolaldehyde'
      logmsg(2) = 'scaled externally'
      CALL set_option(2,spc,"yld",logmsg,qyvers,myld)


      OPEN(UNIT=kin,FILE='DATAJ1/MCMext/ALD/Glycidaldehyde.abs',
     $     STATUS='old')

      do i = 1, 5
        read(kin,*)
      enddo

      n = 84
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
        y1(i) = y1(i)*1.e-20
      ENDDO
      CLOSE(kin)

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j-1),0.,0.,yg)

* quantum yields
      IF(myld==1) THEN
        qy1 = 0.6*0.89
        qy2 = 0.6*0.11
       ELSE
        qy1 = 1.0
        qy2 = 1.0
      ENDIF

      DO iw = 1, nw - 1
        sig = yg(iw)
        DO i = 1, nz
          sq(j-1,i,iw) = sig * qy1
          sq(j  ,i,iw) = sig * qy2
        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE ma10(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! CH3CH=CHCHO

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for CH3CH=CHCHO        =*
*=  (crotonaldehyde) photolysis:                                             =*
*=       CH3CH=CHCHO + hv -> Products                                        =*
*=                                                                           =*
*=  Cross section: from Magneron et al.                                      =*
*=  Quantum yield: estimated 10 times that of acrolein                       =*
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
      PARAMETER(kdata=4500)

      INTEGER iw
      INTEGER i, n
      REAL x1(kdata)
      REAL y1(kdata)

* local

      REAL yg1(kw), yg2(kw)
      real qy, qym1, qy1, qy2, qy3
      REAL sig

      INTEGER mabs, myld, xsvers(3), qyvers(3)
      CHARACTER(llog) :: logmsg(nlog)
      CHARACTER(lspc) :: spc


***************** crotonaldehyde photolysis *************************


* Setting photolysis index and cross section/quantum yield options

      j = j+1
      jlabel(j) = 'CH3CH=CHCHO -> CH3CH=CH + CHO'
      j = j+1
      jlabel(j) = 'CH3CH=CHCHO -> CH3CH=CH2 + CO'
      j = j+1
      jlabel(j) = 'CH3CH=CHCHO -> CH3CH=CHCO + H'

* cross section options
* UV-C: Lee et al. 2007
* UV/VIS:
      spc = "CH3CH=CHCHO"
      xsvers = (/2, 2, 2/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "from Magneron et al. 2002 low res. data"
      logmsg(2) = "from Magneron et al. 1999 unpublished high res. data"
      CALL set_option(2,spc,"abs",logmsg,xsvers,mabs)

* quantum yield options
* Criegee channel replaced by CH3 fission
      qyvers = (/1, 1, 2/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "estimated 10x acrolein with JPL 2006 data"
      logmsg(2) = "estimated 10x acrolein with Calvert et al. 2011 data"
      CALL set_option(2,spc,"yld",logmsg,qyvers,myld)


      OPEN(UNIT=kin,FILE='DATAJ1/MCMext/ALD/crotonaldehyde_Lee07.abs',
     &     STATUS='OLD')
      DO i = 1, 6
        READ(kin,*)
      ENDDO
      n = 4368
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
      ENDDO
      CLOSE(kin)

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j-1),0.,0.,yg1)


      IF(mabs==1) THEN
        OPEN(UNIT=kin,FILE='DATAJ1/MCMext/ALD/crotonaldehyde_Mag02.abs',
     &       STATUS='OLD')
        DO i = 1, 6
          READ(kin,*)
        ENDDO
        n = 68
        DO i = 1, n
          READ(kin,*) x1(i), y1(i)
        ENDDO
        CLOSE(kin)

       ELSEIF(mabs==2) THEN

        OPEN(UNIT=kin,FILE='DATAJ1/MCMext/ALD/crotonaldehyde_Mag99.abs',
     &       STATUS='OLD')
        DO i = 1, 6
          READ(kin,*)
        ENDDO
        n = 3202
        DO i = 1, n
          READ(kin,*) x1(i), y1(i)
        ENDDO
        CLOSE(kin)
      ENDIF

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j-1),0.,0.,yg2)

* combine xs with qy:

      DO iw = 1, nw-1

* cross sections: combine UV-C and UV/VIS-data:
           IF(mabs==1) THEN
              IF(wc(iw)<=250.) THEN
                sig = max(0.,yg1(iw))
               ELSE
                sig = max(0.,yg2(iw))
              ENDIF
            ELSEIF(mabs==2) THEN
              IF(wc(iw)<=225.) THEN
                sig = max(0.,yg1(iw))
               ELSE
                sig = max(0.,yg2(iw))
              ENDIF
           ENDIF

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
               qym1 = -0.836+1.159e-17*8.e17-2.166e-37*8.e17**2
               qy = 1./qym1
             endif
           endif
* product distribution estimated from Calvert et al. 2011:
           qy  = MIN(10. * qy,1.) ! qy estimated 10 times higher than acrolein
           qy1 = qy*(-0.0173*(airden(i)/1.E19)**2
     &               +0.083*airden(i)/1.E19+0.0492)
           qy2 = qy*(0.0407*(airden(i)/1.E19)**2
     &               -0.1661*airden(i)/1.E19+0.8485)
           qy3 = qy*(-0.0217*(airden(i)/1.E19)**2
     &               +0.0788*airden(i)/1.E19+0.1029)
           sq(j-2,i,iw) = sig * qy1
           sq(j-1,i,iw) = sig * qy2
           sq(j  ,i,iw) = sig * qy3
        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE ma11(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! 2-hexenal

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for 2-hexeanal         =*
*=  photolysis:                                                              =*
*=       2-hexenal + hv -> Products                                          =*
*=                                                                           =*
*=  Cross section: see options below                                         =*
*=  Quantum yield: scaled to give j(total) ~ 3e-5 s-1                        =*
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
      PARAMETER(kdata=120)

      INTEGER iw
      INTEGER i, n
      REAL x1(kdata)
      REAL y1(kdata)

* local

      REAL yg(kw)
      real qy, qym1, qy1, qy2, qy3
      REAL sig

      INTEGER mabs, myld, xsvers(3), qyvers(3)
      CHARACTER(llog) :: logmsg(nlog)
      CHARACTER(lspc) :: spc


***************** 2-hexenal photolysis *************************


* Setting photolysis index and cross section/quantum yield options
      j = j+1
      jlabel(j) = '2-hexenal -> 1-pentenyl radical + CHO'
      j = j+1
      jlabel(j) = '2-hexenal -> 1-pentene + CO'
      j = j+1
      jlabel(j) = '2-hexenal -> C3H7CH=CHCO + H'

* cross section options
      spc = "2-hexenal"
      xsvers = (/1, 2, 1/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "from O'Connor et al. 2006"
      logmsg(2) = "from Jiminez et al. 2007"
      CALL set_option(2,spc,"abs",logmsg,xsvers,mabs)

* quantum yield options
* Criegee channel replaced by CH3 fission
      qyvers = (/1, 1, 2/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "estimated 10x acrolein with JPL 2006 data"
      logmsg(2) = "estimated 10x acrolein with Calvert et al. 2011 data"
      CALL set_option(2,spc,"yld",logmsg,qyvers,myld)


* cross sections
      IF(mabs==1) THEN
        OPEN(UNIT=kin,FILE='DATAJ1/MCMext/ALD/hexenal_OCon06.abs',
     &       STATUS='OLD')
        DO i = 1, 6
          READ(kin,*)
        ENDDO

        n = 111
        DO i = 1, n
          READ(kin,*) x1(i), y1(i)
        ENDDO
        CLOSE(kin)

       ELSEIF(mabs==2) THEN

        OPEN(UNIT=kin,FILE='DATAJ1/MCMext/ALD/hexenal_Jim07.abs',
     &       STATUS='OLD')
        DO i = 1, 6
          READ(kin,*)
        ENDDO

        n = 81
        DO i = 1, n
          READ(kin,*) x1(i), y1(i)
        ENDDO
        CLOSE(kin)
      ENDIF

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j-1),0.,0.,yg)

* combine xs with qy:

      DO iw = 1, nw-1

        sig = yg(iw)

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
               qym1 = -0.836+1.159e-17*8.e17-2.166e-37*8.e17**2
               qy = 1./qym1
             endif
           endif
* product distribution estimated from Calvert et al. 2011:
           qy = MIN(12. * qy,1.) ! qy estimated 12 times higher than acrolein
           qy1 = qy*(-0.0173*(airden(i)/1.E19)**2
     &               +0.083*airden(i)/1.E19+0.0492)
           qy2 = qy*(0.0407*(airden(i)/1.E19)**2
     &               -0.1661*airden(i)/1.E19+0.8485)
           qy3 = qy*(-0.0217*(airden(i)/1.E19)**2
     &               +0.0788*airden(i)/1.E19+0.1029)
           sq(j-2,i,iw) = sig * qy1
           sq(j-1,i,iw) = sig * qy2
           sq(j  ,i,iw) = sig * qy3
        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE ma12(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! C4H9CH(C2H5)CHO / alkyl-subst. ald.

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for                =*
*=  2-ethyl hexanal photolysis:                                              =*
*=         C4H9CH(C2H5)CHO + hv -> products                                  =*
*=                                                                           =*
*=  Cross section:  Fraire et al. (2011)                                     =*
*=  Quantum yield:  Fraire et al. (2011) @254nm                              =*
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
      real qy1, qy2
      REAL sig, sigoh
      INTEGER iw
      INTEGER myld, qyvers(3)
      CHARACTER(llog) :: logmsg(nlog)
      CHARACTER(lspc) :: spc


***************** C4H9CH(C2H5)CHO photolysis *************************

* Setting photolysis index and quantum yield options

      j = j+1
      jlabel(j) = "C4H9CH(C2H5)CHO -> C7H15 + CHO"
      j = j+1
      jlabel(j) = "C4H9CH(C2H5)CHO -> C7H16 + CO"
      j = j+1
      jlabel(j) = "intAldOH -> R + CHO"
      j = j+1
      jlabel(j) = "intAldOH -> R' + CO"

* quantum yield options
      spc = "C4H9CH(C2H5)CHO"
      qyvers = (/1, 2, 1/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "by Fraire et al. (2011) @254nm"
      logmsg(2) = "scaled externally"
      CALL set_option(2,spc,"yld",logmsg,qyvers,myld)


* Absorption cross sections

      OPEN(UNIT=kin,FILE='DATAJ1/MCMext/ALD/2EthylHexanal.abs',
     $     STATUS='old')
      do i = 1, 5
        read(kin,*)
      enddo

      n = 33
      DO i = 1, n
        READ(kin,*) x1(i), y1(i), dum
      ENDDO
      CLOSE(kin)

      x1oh(:) = 0.
      x1oh(:) = x1(:) - 10.
      y1oh(:) = y1(:)
      n1 = n
      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j-3),0.,0.,yg)
      n = n1
      CALL interpol(x1oh,y1oh,kdata,n,nw,wl,jlabel(j-1),0.,0.,ygoh)

* quantum yields
* measurements at 254nm extended to all other wavelength
* branching of qy1/qy2 of 0.8:0.2 assumed (total qy = 0.51)

* zero pressure qy:
      IF(myld==1) THEN
        qy1  = 0.51*0.8
        qy2  = 0.51*0.2
      ELSEIF(myld==2) THEN
        qy1 = 1.0
        qy2 = 1.0
      ENDIF


* combine:

      DO iw = 1, nw - 1
         sig = yg(iw)
         sigoh = ygoh(iw)
         DO i = 1, nz
           sq(j-3,i,iw) = sig * qy1
           sq(j-2,i,iw) = sig * qy2
           sq(j-1,i,iw) = sigoh * qy1
           sq(j  ,i,iw) = sigoh * qy2
         ENDDO
      ENDDO

      END

* ============================================================================*

      SUBROUTINE ma13(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! CH3CH=C(CH3)CHO

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for CH3CH=C(CH3)CHO    =*
*=  photolysis:                                                              =*
*=       CH3CH=C(CH3)CHO + hv -> Products                                    =*
*=                                                                           =*
*=  Cross section: Lanza et al. (2008)                                       =*
*=  Quantum yield: estimated 10 times that of acrolein                       =*
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
      INTEGER i, n, n1
      REAL x1(kdata), x1oh(kdata)
      REAL y1(kdata), y1oh(kdata)

* local

      REAL yg1(kw), ygoh(kw)
      real qy, qym1, qy1, qy2, qy3
      REAL sig, sigoh

      INTEGER myld, qyvers(3)
      CHARACTER(llog) :: logmsg(nlog)
      CHARACTER(lspc) :: spc

***************** 2-methyl crotonaldehyde photolysis *************************


* Setting photolysis index and quantum yield options

      j = j+1
      jlabel(j) = 'CH3CH=C(CH3)CHO -> CH3CH=CCH3 + CHO'
      j = j+1
      jlabel(j) = 'CH3CH=C(CH3)CHO -> CH3CH=CHCH3 + CO'
      j = j+1
      jlabel(j) = 'CH3CH=C(CH3)CHO -> CH3CH=C(CH3)CO + H'

      j = j+1
      jlabel(j) = 'aMeC4uALDOH -> NI products'
      j = j+1
      jlabel(j) = 'aMeC4uALDOH -> alkene + CO'
      j = j+1
      jlabel(j) = 'aMeC4uALDOH -> acyl + H'

      j = j+1
      jlabel(j) = 'genaMeuAld(poly)'
      j = j+1
      jlabel(j) = 'genaMeuAldOH(poly)'

* quantum yield options
* Criegee channel replaced by CH3 fission
      spc = "CH3CH=C(CH3)CHO"
      qyvers = (/1, 1, 2/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "estimated 10x acrolein with JPL 2006 data"
      logmsg(2) = "estimated 10x acrolein with Calvert et al. 2011 data"
      CALL set_option(2,spc,"yld",logmsg,qyvers,myld)


* cross section

      OPEN(UNIT=kin,FILE='DATAJ1/MCMext/ALD/2Me2butenal.abs',
     &     STATUS='OLD')
      DO i = 1, 4
        READ(kin,*)
      ENDDO
      n = 84
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
        y1(i) = y1(i)*1.E-20
      ENDDO
      CLOSE(kin)

      x1oh(:) = 0.
      x1oh(:) = x1(:) - 10.
      y1oh(:) = y1(:)
      n1 = n
      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg1)
      n = n1
      CALL interpol(x1oh,y1oh,kdata,n,nw,wl,jlabel(j),0.,0.,ygoh)


* combine xs with qy:

      DO iw = 1, nw-1
        sig = yg1(iw)
        sigoh = ygoh(iw)
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
               qym1 = -0.836+1.159e-17*8.e17-2.166e-37*8.e17**2
               qy = 1./qym1
             endif
           endif
* product distribution as estimated 10 x acrolein from Calvert et al. 2011:
* fits of direct qy given behind do not match 0.0065 at 1 atm
           qy = MIN(10. * qy,1.)
           qy1 = qy*(-0.0173*(airden(i)/1.E19)**2
     &              +0.083*airden(i)/1.E19+0.0492)
           qy2 = qy*(0.0407*(airden(i)/1.E19)**2
     &              -0.1661*airden(i)/1.E19+0.8485)
           qy3 = qy*(-0.0217*(airden(i)/1.E19)**2
     &              +0.0788*airden(i)/1.E19+0.1029)
           sq(j-7,i,iw) = sig * qy1
           sq(j-6,i,iw) = sig * qy2
           sq(j-5,i,iw) = sig * qy3
           IF(swOH==1) THEN
             sq(j-4,i,iw) = sigoh * qy1
             sq(j-3,i,iw) = sigoh * qy2
             sq(j-2,i,iw) = sigoh * qy3
           ELSEIF(swOH==2) THEN
             IF(qy1+qy2+qy3 > 0.) THEN
               sq(j-4,i,iw) = sigoh * qyoh*qy1/(qy1+qy2+qy3)
               sq(j-3,i,iw) = sigoh * qyoh*qy2/(qy1+qy2+qy3)
               sq(j-2,i,iw) = sigoh * qyoh*qy3/(qy1+qy2+qy3)
              ELSE
               sq(j-4,i,iw) = 0.
               sq(j-3,i,iw) = 0.
               sq(j-2,i,iw) = 0.
             ENDIF
           ENDIF
           sq(j-1,i,iw) = sig
           sq(j  ,i,iw) = sigoh
        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE ma14(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! CH3C(CH3)=CHCHO

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for CH3C(CH3)=CHCHO    =*
*=  photolysis:                                                              =*
*=       CH3C(CH3)=CHCHO + hv -> Products                                    =*
*=                                                                           =*
*=  Cross section: Lanza et al. (2008)                                       =*
*=  Quantum yield: estimated 10 times that of acrolein                       =*
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
      INTEGER i, n, n1
      REAL x1(kdata), x1oh(kdata)
      REAL y1(kdata), y1oh(kdata)

* local

      REAL yg1(kw), ygoh(kw)
      real qy, qym1, qy1, qy2, qy3
      REAL sig, sigoh

      INTEGER myld, qyvers(3)
      CHARACTER(llog) :: logmsg(nlog)
      CHARACTER(lspc) :: spc

***************** 3-methyl crotonaldehyde photolysis *************************


* Setting photolysis index and quantum yield options

      j = j+1
      jlabel(j) = 'CH3C(CH3)=CHCHO -> (CH3)2C=CH + CHO'
      j = j+1
      jlabel(j) = 'CH3C(CH3)=CHCHO -> (CH3)2C=CH2 + CO'
      j = j+1
      jlabel(j) = 'CH3C(CH3)=CHCHO -> (CH3)2C=CHCO + H'

      j = j+1
      jlabel(j) = 'bMeC4uALDOH -> NI products'
      j = j+1
      jlabel(j) = 'bMeC4uALDOH -> alkene + CO'
      j = j+1
      jlabel(j) = 'bMeC4uALDOH -> acyl + H'

      j = j+1
      jlabel(j) = 'genbMeuAld(poly)'
      j = j+1
      jlabel(j) = 'genbMeuAldOH(poly)'


* quantum yield options
* Criegee channel replaced by CH3 fission
      spc = "CH3C(CH3)=CHCHO"
      qyvers = (/1, 1, 2/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "estimated 10x acrolein with JPL 2006 data"
      logmsg(2) = "estimated 10x acrolein with Calvert et al. 2011 data"
      CALL set_option(2,spc,"yld",logmsg,qyvers,myld)


* cross section

      OPEN(UNIT=kin,FILE='DATAJ1/MCMext/ALD/3Me2butenal.abs',
     &     STATUS='OLD')
      DO i = 1, 4
        READ(kin,*)
      ENDDO
      n = 84
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
      ENDDO
      CLOSE(kin)

      x1oh(:) = 0.
      x1oh(:) = x1(:) - 10.
      y1oh(:) = y1(:)
      n1 = n
      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg1)
      n = n1
      CALL interpol(x1oh,y1oh,kdata,n,nw,wl,jlabel(j),0.,0.,ygoh)


* combine xs with qy:

      DO iw = 1, nw-1
        sig = yg1(iw)
        sigoh = ygoh(iw)
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
               qym1 = -0.836+1.159e-17*8.e17-2.166e-37*8.e17**2
               qy = 1./qym1
             endif
           endif
* product distribution as estimated 10 x acrolein from Calvert et al. 2011:
* fits of direct qy given behind do not match 0.0065 at 1 atm
           qy = MIN(10. * qy,1.)
           qy1 = qy*(-0.0173*(airden(i)/1.E19)**2
     &              +0.083*airden(i)/1.E19+0.0492)
           qy2 = qy*(0.0407*(airden(i)/1.E19)**2
     &              -0.1661*airden(i)/1.E19+0.8485)
           qy3 = qy*(-0.0217*(airden(i)/1.E19)**2
     &                +0.0788*airden(i)/1.E19+0.1029)
           sq(j-7,i,iw) = sig * qy1
           sq(j-6,i,iw) = sig * qy2
           sq(j-5,i,iw) = sig * qy3
           IF(swOH==1) THEN
             sq(j-4,i,iw) = sigoh * qy1
             sq(j-3,i,iw) = sigoh * qy2
             sq(j-2,i,iw) = sigoh * qy3
           ELSEIF(swOH==2) THEN
             IF(qy1+qy2+qy3 > 0.) THEN
               sq(j-4,i,iw) = sigoh * qyoh*qy1/(qy1+qy2+qy3)
               sq(j-3,i,iw) = sigoh * qyoh*qy2/(qy1+qy2+qy3)
               sq(j-2,i,iw) = sigoh * qyoh*qy3/(qy1+qy2+qy3)
              ELSE
               sq(j-4,i,iw) = 0.
               sq(j-3,i,iw) = 0.
               sq(j-2,i,iw) = 0.
             ENDIF
           ENDIF
           sq(j-1,i,iw) = sig
           sq(j  ,i,iw) = sigoh
        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE ma15(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! 2,4-hexadienal

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for 2,4-hexadienal     =*
*=  photolysis:                                                              =*
*=       2,4-hexadienal + hv -> Products                                     =*
*=                                                                           =*
*=  Cross section: O'Connor et al. (2006)                                    =*
*=  Quantum yield: estimated 10*acrolein                                     =*
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
      PARAMETER(kdata=120)

      INTEGER iw
      INTEGER i, n, n1
      REAL x1(kdata), x1oh(kdata)
      REAL y1(kdata), y1oh(kdata)

* local

      REAL yg(kw), ygoh(kw)
      real qy, qym1, qy1, qy2, qy3
      REAL sig, sigoh

***************** 2-hexenal photolysis *************************

      j = j+1
      jlabel(j) = 'hexadienal -> 1-pentenyl radical + CHO'
      j = j+1
      jlabel(j) = 'hexadienal -> 1,3-pentadiene + CO'
      j = j+1
      jlabel(j) = 'hexadienal -> CH3CH=CHCH=CHCO + H'

      j = j+1
      jlabel(j) = 'uuALDOH -> NI products'
      j = j+1
      jlabel(j) = 'uuALDOH -> diene + CO'
      j = j+1
      jlabel(j) = 'uuALDOH -> acyl + H'

      j = j+1
      jlabel(j) = 'genluuALD(poly)'
      j = j+1
      jlabel(j) = 'genluuALD(OHpoly)'


* cross sections
      OPEN(UNIT=kin,FILE='DATAJ1/MCMext/ALD/hexadienal.abs',
     &     STATUS='OLD')
      DO i = 1, 5
        READ(kin,*)
      ENDDO

      n = 111
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
      ENDDO
      CLOSE(kin)

      x1oh(:) = 0.
      x1oh(:) = x1(:) - 10.
      y1oh(:) = y1(:)
      n1 = n
      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)
      n = n1
      CALL interpol(x1oh,y1oh,kdata,n,nw,wl,jlabel(j),0.,0.,ygoh)

* combine xs with qy:

      DO iw = 1, nw-1

        sig = yg(iw)
        sigoh = ygoh(iw)

        DO i = 1, nz

           if(airden(i) .gt. 2.6e19) then
             qy = 0.004
            elseif(airden(i) .gt. 8.e17 .and. airden(i) .lt. 2.6e19)
     &        then
             qym1 = 0.086 + 1.613e-17 * airden(i)
             qy = 0.004 + 1./qym1
            elseif(airden(i) .lt. 8.e17) then
             qym1 = 0.086 + 1.613e-17 * 8.e17
             qy = 0.004 + 1./qym1
           endif
* product distribution estimated from Calvert et al. 2011:
           qy = MIN(10. * qy,1.) ! qy estimated 10 times higher than acrolein
           qy1 = qy*(-0.0173*(airden(i)/1.E19)**2
     &               +0.083*airden(i)/1.E19+0.0492)
           qy2 = qy*(0.0407*(airden(i)/1.E19)**2
     &               -0.1661*airden(i)/1.E19+0.8485)
           qy3 = qy*(-0.0217*(airden(i)/1.E19)**2
     &               +0.0788*airden(i)/1.E19+0.1029)
           sq(j-7,i,iw) = sig * qy1
           sq(j-6,i,iw) = sig * qy2
           sq(j-5,i,iw) = sig * qy3
           IF(swOH==1) THEN
             sq(j-4,i,iw) = sigoh * qy1
             sq(j-3,i,iw) = sigoh * qy2
             sq(j-2,i,iw) = sigoh * qy3
           ELSEIF(swOH==2) THEN
             IF(qy1+qy2+qy3 > 0.) THEN
               sq(j-4,i,iw) = sigoh * qyoh*qy1/(qy1+qy2+qy3)
               sq(j-3,i,iw) = sigoh * qyoh*qy2/(qy1+qy2+qy3)
               sq(j-2,i,iw) = sigoh * qyoh*qy3/(qy1+qy2+qy3)
              ELSE
               sq(j-4,i,iw) = 0.
               sq(j-3,i,iw) = 0.
               sq(j-2,i,iw) = 0.
             ENDIF
           ENDIF
           sq(j-1,i,iw) = sig
           sq(j  ,i,iw) = sigoh
        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE ma16(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! (CH3)2C(OH)CHO

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section)x(quantum yield) for (CH3)2C(OH)CHO       =*
*=  photolysis                                                               =*
*=          (CH3)2C(OH)CHO + hv -> (CH3)2COH + CHO                           =*
*=                                                                           =*
*=  Cross section: Chakir et al. (2004)                                      =*
*=  Quantum yield: estimated same as i-C3H7CHO                               =*
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

      REAL yg(kw), yg1(kw), dum
      REAL qy
      REAL sig
      INTEGER iw

      INTEGER myld, qyvers(3)
      CHARACTER(llog) :: logmsg(nlog)
      CHARACTER(lspc) :: spc

***************** (CH3)2C(OH)CHO photolysis *************************


* Setting photolysis index and quantum yield options

* only Norish type I based on i-C3H7CHO
      j = j+1
      jlabel(j) = 'CH32C(OH)CHO -> (CH3)2COH + CHO'

      spc = "CH32C(OH)CHO"
      qyvers = (/1, 4, 4/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = 'from Desai et al. 1986 for i-C3H7CHO'
      logmsg(2) = 'opt. 1 from Calvert book for i-C3H7CHO'
      logmsg(3) = 'opt. 2 from Calvert book for i-C3H7CHO'
      logmsg(4) = 'from IUPAC for i-C3H7CHO'
      CALL set_option(4,spc,"yld",logmsg,qyvers,myld)


* cross sections

      OPEN(UNIT=kin,FILE='DATAJ1/MCMext/ALD/bOHisobutyraldehyde.abs',
     $     STATUS='old')
      do i = 1, 4
         read(kin,*)
      enddo

      n = 12
      DO i = 1, n
         READ(kin,*) x1(i), y1(i)
      ENDDO
      CLOSE(kin)

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)


* quantum yields

      IF(myld==1) THEN

        OPEN(UNIT=kin,FILE='DATAJ1/MCMext/ALD/i_butyraldehyde_R.prn',
     $       STATUS='old')
        do i = 1, 3
           read(kin,*)
        enddo

        n = 101
        DO i = 1, n
           READ(kin,*) x1(i), dum, y1(i)
        ENDDO
        CLOSE(kin)

       ELSEIF(myld==2 .OR. myld==3) THEN

        OPEN(UNIT=kin,FILE='DATAJ1/MCMext/ALD/iC3H7CHO_calv.yld',
     $       STATUS='old')
        do i = 1, 6
           read(kin,*)
        enddo

        n = 73
        DO i = 1, n
           READ(kin,*) x1(i), dum, y1(i),y2(i)
           x2(i) = x1(i)
        ENDDO
        CLOSE(kin)

       ELSEIF(myld==4) THEN

        OPEN(UNIT=kin,FILE='DATAJ1/MCMext/ALD/iC3H7CHO_iup.yld',
     $       STATUS='old')
        do i = 1, 7
           read(kin,*)
        enddo

        n = 11
        DO i = 1, n
           READ(kin,*) x1(i), y1(i)
        ENDDO
        CLOSE(kin)

      ENDIF

      IF(myld==3)THEN
        CALL interpol(x2,y2,kdata,n,nw,wl,jlabel(j),0.,0.,yg1)
       ELSE
        CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg1)
      ENDIF

*combine
      DO iw = 1, nw - 1
         sig = yg (iw)
         qy  = yg1(iw)
         DO i = 1, nz
            sq(j  ,i,iw) = sig * qy
         ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE ma17(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) !C≥8 n-aldehyde

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section)x(quantum yield) for n-octanal photolysis =*
*=          n-C7H15CHO + hv -> products                                      =*
*=  Additional estimates for larger n-aldehydes                              =*
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
      PARAMETER(kdata=580)

      INTEGER i, n, n1
      REAL x1(kdata), x1oh(kdata)
      REAL y1(kdata), y1oh(kdata)

* local

      REAL yg(kw), ygoh(kw)
      REAL yld1, yld2, yld3
      REAL qy1, qy2, qy3, qyN, qy, ptorr, eta
      REAL sig, sigoh
      INTEGER iw

      INTEGER myld, qyvers(3)
      CHARACTER(llog) :: logmsg(nlog)
      CHARACTER(lspc) :: spc

***************** n-C7H15CHO photolysis *************************


* Setting photolysis index and quantum yield options
      j = j+1
      jlabel(j) = 'n-C7H15CHO -> C7H15 + CHO'
      j = j+1
      jlabel(j) = 'n-C7H15CHO -> C4H7CH=CH2 + CH2=CHOH'
      j = j+1
      jlabel(j) = 'n-C7H15CHO -> 2-butylcyclobutanol'

      j = j+1
      jlabel(j) = 'C8nALDOH -> NI products'
      j = j+1
      jlabel(j) = 'C8nALDOH -> NII products'
      j = j+1
      jlabel(j) = 'C8nALDOH -> cycl. product'
* estimates in compounds with polyfunctional chromophores
      j = j+1
      jlabel(j) = 'nALDpoly(C>7)'
      j = j+1
      jlabel(j) = 'nALDOHpoly(C>7)'

      spc = "n-C7H15CHO"
      qyvers = (/1, 3, 1/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = 'from Tadic et al. 2011'
      logmsg(2) =
     &  'estimated 7% Norrish I, 74.4% Norrish II and 18.6% cyclisation'
      logmsg(3) = "scaled externally"
      CALL set_option(3,spc,"yld",logmsg,qyvers,myld)


* cross sections

      OPEN(UNIT=kin,FILE='DATAJ1/MCMext/ALD/nAldehydes.abs',
     $     STATUS='old')

      do i = 1, 12
        read(kin,*)
      enddo

      n = 128
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


* quantum yields

      IF(myld==1) THEN
        yld1 = 0.07
        yld2 = 0.34
        yld3 = 1. - yld1 - yld2
        qyN  = 0.32
      ELSEIF(myld==2) THEN
        yld1 = 0.07
        yld2 = 0.74
        yld3 = 0.19
        qyN  = 0.32
      ELSEIF(myld==3) THEN
        yld1 = 1.0
        yld2 = 1.0
        yld3 = 1.0
        qyN  = 1.0
      ENDIF

* pressure dependence
      eta = 700. * 1.061e-3 / 2.368


* combine:

      DO iw = 1, nw - 1
        sig = yg(iw)
        sigoh = ygoh(iw)
        DO i = 1, nz
          ptorr = 760.*airden(i)/2.55e19 !2.46?
          qy = qyN * (1. + eta) / (1 + eta * airden(i)/2.27e19) !2.27e19 = 700Torr
          qy1 = yld1*qy
          qy2 = yld2*qy
          qy3 = yld3*qy
          sq(j-7,i,iw) = sig * qy1
          sq(j-6,i,iw) = sig * qy2
          sq(j-5,i,iw) = sig * qy3
          IF(swOH==1) THEN
            sq(j-4,i,iw) = sigoh * qy1
            sq(j-3,i,iw) = sigoh * qy2
            sq(j-2,i,iw) = sigoh * qy3
          ELSEIF(swOH==2) THEN
            IF(qy1+qy2+qy3 > 0.) THEN
              sq(j-4,i,iw) = sigoh * qyoh*qy1/(qy1+qy2+qy3)
              sq(j-3,i,iw) = sigoh * qyoh*qy2/(qy1+qy2+qy3)
              sq(j-2,i,iw) = sigoh * qyoh*qy3/(qy1+qy2+qy3)
             ELSE
              sq(j-4,i,iw) = 0.
              sq(j-3,i,iw) = 0.
              sq(j-2,i,iw) = 0.
            ENDIF
            sq(j-1,i,iw) = sig
            sq(j  ,i,iw) = sigoh
          ENDIF
        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE ma18(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! linear alpha,beta-unsaturated aldehydes

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section)x(quantum yield) for                      =*
*=  linear alpha,beta-unsaturated aldehyde photolysis:                       =*
*=          aldehydes + hv ->  products                                      =*
*=  Cross section:  average of crotonaldehyde and 2-hexenal                  =*
*=  Quantum yield:  10 * acrolein                                            =*
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
      PARAMETER(kdata=7000)

      INTEGER iw
      INTEGER i, n, n1
      REAL x1(kdata), x1oh(kdata)
      REAL y1(kdata), y1oh(kdata)

* local

      REAL yg(kw), ygoh(kw)
      real qy, qym1, qy1, qy2, qy3
      REAL sig, sigoh

**************** crotonaldehyde photodissociation

      j = j+1
      jlabel(j) = 'luALD -> NI products'
      j = j+1
      jlabel(j) = 'luALD -> alkene + CO'
      j = j+1
      jlabel(j) = 'luALD -> acyl + H'

      j = j+1
      jlabel(j) = 'luALDOH -> NI products'
      j = j+1
      jlabel(j) = 'luALDOH -> alkene + CO'
      j = j+1
      jlabel(j) = 'luALDOH -> acyl + H'

      j = j+1
      jlabel(j) = 'genluALD(poly)'
      j = j+1
      jlabel(j) = 'genluALD(OHpoly)'


* cross section from
* UV-C: Lee et al. 2007
* UV/VIS: unpublished high res data by Magneron et al. 1999
* quantum yields estimated 10x acrolein

      OPEN(UNIT=kin,FILE='DATAJ1/MCMext/ALD/uAld.abs',
     &     STATUS='OLD')
      DO i = 1, 11
        READ(kin,*)
      ENDDO
      n = 6848
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

* combine xs with qy:

      DO iw = 1, nw-1

        sig = yg(iw)
        sigoh = ygoh(iw)

        DO i = 1, nz

           if(airden(i) .gt. 2.6e19) then
             qy = 0.004
            elseif(airden(i) .gt. 8.e17 .and. airden(i) .lt. 2.6e19)
     &        then
             qym1 = 0.086 + 1.613e-17 * airden(i)
             qy = 0.004 + 1./qym1
            elseif(airden(i) .lt. 8.e17) then
             qym1 = 0.086 + 1.613e-17 * 8.e17
             qy = 0.004 + 1./qym1
           endif
* product distribution estimated from Calvert et al. 2011:
           qy  = MIN(10. * qy,1.) ! qy estimated 10 times higher than acrolein
           qy1 = qy*(-0.0173*(airden(i)/1.E19)**2
     &               +0.083*airden(i)/1.E19+0.0492)
           qy2 = qy*(0.0407*(airden(i)/1.E19)**2
     &               -0.1661*airden(i)/1.E19+0.8485)
           qy3 = qy*(-0.0217*(airden(i)/1.E19)**2
     &               +0.0788*airden(i)/1.E19+0.1029)
           sq(j-7,i,iw) = sig * qy1
           sq(j-6,i,iw) = sig * qy2
           sq(j-5,i,iw) = sig * qy3
           IF(swOH==1) THEN
             sq(j-4,i,iw) = sigoh * qy1
             sq(j-3,i,iw) = sigoh * qy2
             sq(j-2,i,iw) = sigoh * qy3
           ELSEIF(swOH==2) THEN
             IF(qy1+qy2+qy3 > 0.) THEN
               sq(j-4,i,iw) = sigoh * qyoh*qy1/(qy1+qy2+qy3)
               sq(j-3,i,iw) = sigoh * qyoh*qy2/(qy1+qy2+qy3)
               sq(j-2,i,iw) = sigoh * qyoh*qy3/(qy1+qy2+qy3)
              ELSE
               sq(j-4,i,iw) = 0.
               sq(j-3,i,iw) = 0.
               sq(j-2,i,iw) = 0.
             ENDIF
           ENDIF
           sq(j-1,i,iw) = sig
           sq(j  ,i,iw) = sigoh
        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE ma19(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! neopentylaldehyde

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for                =*
*=  neopentyl aldehyde photolysis:                                           =*
*=         neoC5H11CHO + hv -> Norish type I + II products                   =*
*=                                                                           =*
*=  Cross section:  estimated same as isopentyl aldehyde                     =*
*=  Quantum yield:  Tadic et al. 2012                                        =*
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
      REAL yld1, yld2, yld3, qy1, qy2, qy3, qy0, ptorr
      REAL sig
      INTEGER iw

      INTEGER myld, qyvers(3)
      CHARACTER(llog) :: logmsg(nlog)
      CHARACTER(lspc) :: spc


***************** neopenty laldehyde photolysis *************************


* Setting photolysis index and cross section/quantum yield options
      j = j+1
      jlabel(j) = 'neoC5H11CHO -> neoC5H11 + CHO'
      j = j+1
      jlabel(j) = 'neoC5H11CHO -> CH32C=CH2 + CH2=CHOH'
      j = j+1
      jlabel(j) = 'neoC5H11CHO -> 3,3-dimethyl-1-cyclobutanol'

* quantum yield options
      spc = "neoC5H11CHO"
      qyvers = (/1, 2, 2/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "from Tadic et al. 2012"
      logmsg(2) =
     & "from Tadic et al. 2012 with cyclisation estimated as Norrish-II"
      CALL set_option(2,spc,"yld",logmsg,qyvers,myld)


* cross sections
      OPEN(UNIT=kin,
     $     FILE='DATAJ1/MCMext/ALD/isovaleraldehyde_calvext.abs',
     $     STATUS='old')
      do i = 1, 6
        read(kin,*)
      enddo

      n = 71
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
        y1(i) = y1(i) * 1.E-20
      ENDDO
      CLOSE(kin)

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j-1),0.,0.,yg)


* product yields
      IF(myld==1) THEN
        yld1 = 0.2
        yld2 = 0.55
        yld3 = 0.25
       ELSEIF(myld==2) THEN
         yld1 = 0.2
         yld2 = 0.8
         yld3 = 0.0
      ENDIF

* combine:

      DO iw = 1, nw - 1
         sig = yg(iw)
         ptorr = 760.*airden(i)/2.55e19
         qy0 = 1.0/(1.85 + 1.381E-3*ptorr)
         qy1 = yld1*qy0
         qy2 = yld2*qy0
         qy3 = yld3*qy0
         DO i = 1, nz
           sq(j-2,i,iw) = sig * qy1
           sq(j-1,i,iw) = sig * qy2
           sq(j  ,i,iw) = sig * qy3
         ENDDO
      ENDDO

      END

* ============================================================================*
