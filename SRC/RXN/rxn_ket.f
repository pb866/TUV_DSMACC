*= This file contains the following subroutines, related to reading/loading
*= the product (cross section) x (quantum yield) for photo-reactions of
*= ketones in MCM-GECKO, which were not yet present in TUV5.2:
*=
*=     mk01 through mk36

*=============================================================================*

      SUBROUTINE mk01(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! diethyl ketone

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for                =*
*=  C2H5COC2H5 photolysis:                                                   =*
*=        C2H5COC2H5 -> products                                        =*
*=                                                                           =*
*=  Cross section:  see options below                                        =*
*=  Quantum yield:  estimated (0.85:0.15)                                    =*
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
      REAL A, lc, w

      INTEGER mabs, myld, xsvers(3), qyvers(3)
      CHARACTER(llog) :: logmsg(nlog)
      CHARACTER(lspc) :: spc

***************** C2H5COC2H5 (DEK) photolysis *************************

* Setting photolysis index and cross section options
      j = j+1
      jlabel(j) = 'C2H5COC2H5 -> 2 C2H5 + CO'
      j = j+1
      jlabel(j) = 'C2H5COC2H5 -> C2H5CO + C2H5'

      spc = "C2H5COC2H5"
      xsvers = (/1, 4, 1/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "from Martinez et al. 1992"
      logmsg(2) = "from Horowitz 1999"
      logmsg(3) = "from Koch et al. 2008"
      logmsg(4) = "from T-dep. gaussian approximation"
      CALL set_option(4,spc,"abs",logmsg,xsvers,mabs)

      qyvers = (/1, 2, 1/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = 'estimated based on Calvert et al. (2011)' ! (RADICAL 2002)
      logmsg(2) = 'scaled externally'
      CALL set_option(2,spc,"yld",logmsg,qyvers,myld)


* cross sections

      IF(mabs==1) THEN
        OPEN(UNIT=kin,FILE='DATAJ1/ABS/Martinez.abs',
     $       STATUS='old')
        do i = 1, 4
          read(kin,*)
        enddo

        n  = 96
        DO i = 1, n
          READ(kin,*) x1(i), dum, dum, dum, y1(i)
          y1(i) = y1(i)*1.E-20
        ENDDO
        CLOSE(kin)

       ELSEIF(mabs==2) THEN

        OPEN(UNIT=kin,FILE='DATAJ1/MCMext/KET/DEK_Hor99.abs',
     $       STATUS='old')
        do i = 1, 7
          read(kin,*)
        enddo

        n  = 2196
        DO i = 1, n
          READ(kin,*) x1(i), y1(i)
        ENDDO
        CLOSE(kin)

       ELSEIF(mabs>=3) THEN

        OPEN(UNIT=kin,FILE='DATAJ1/MCMext/KET/DEK_Koch08.abs',
     $       STATUS='old')
        do i = 1, 6
          read(kin,*)
        enddo

        n  = 159
        DO i = 1, n
          READ(kin,*) x1(i), y1(i)
        ENDDO
        CLOSE(kin)
      ENDIF

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)


* quantum yields
      IF(myld==1) THEN
        qy1 = 0.85
        qy2 = 0.15
       ELSE
        qy1 = 1.0
        qy2 = 1.0
      ENDIF


* combine:
      DO iw = 1, nw - 1
        sig = yg(iw)
        DO i = 1, nz
          IF(mabs==4 .AND. wc(iw)>=230. .AND. wc(iw)<=330.) THEN
            A   = 4.77e-20 + 4.87e-23*tlev(i)
            lc  = 273.7    + 0.0186  *tlev(i)
            w   =  25.1    + 0.0101  *tlev(i)
            sig = A * exp(-((wc(iw)-lc)/w)**2)
          ENDIF
          sq(j-1,i,iw) = sig * qy1
          sq(j  ,i,iw) = sig * qy2
        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE mk02(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! methyl n-propyl ketone

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for                =*
*=  C3H7COCH3 photolysis:                                                    =*
*=        C3H7COCH3 -> Norish type I + II products                      =*
*=                                                                           =*
*=  Cross section:  see options below                                        =*
*=  Quantum yield:  estimates based on Calvert et al. 2011 book              =*
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

      INTEGER i, n, n1, n2, n3
      REAL x1(kdata),x2(kdata),x3(kdata)
      REAL y1(kdata),y2(kdata),y3(kdata)

* local

      REAL yg(kw), yg1(kw), yg2(kw), yg3(kw), dum
      REAL qy1, qy2, qy3, qy4
      REAL sig
      INTEGER iw
      REAL A, lc, w

      INTEGER mabs, xsvers(3)
      CHARACTER(llog) :: logmsg(nlog)
      CHARACTER(lspc) :: spc

***************** C3H7COCH3 photolysis *************************


* Setting photolysis index and cross section options
      j = j+1
      jlabel(j) = 'C3H7COCH3 -> CH3CO + C3H7'
      j = j+1
      jlabel(j) = 'C3H7COCH3 -> C3H7CO + CH3'
      j = j+1
      jlabel(j) = 'C3H7COCH3 -> C3H7 + CO + CH3'
      j = j+1
      jlabel(j) = 'C3H7COCH3 -> CH3C(OH)=CH2 + CH2=CH2'

      spc = "C3H7COCH3"
      xsvers = (/1, 2, 1/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "from Martinez et al. 1992"
      logmsg(2) = "from T-dep. gaussian approximation"
      CALL set_option(2,spc,"abs",logmsg,xsvers,mabs)


* cross sections

      OPEN(UNIT=kin,FILE='DATAJ1/ABS/Martinez.abs',
     $     STATUS='old')
      do i = 1,4
         read(kin,*)
      enddo

      n = 96
      DO i = 1, n
         READ(kin,*) x1(i),dum,dum,y1(i),dum
         y1(i) = y1(i)*1.E-20
      ENDDO
      CLOSE(kin)

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)


* quantum yields

      OPEN(UNIT=kin,FILE='DATAJ1/MCMext/KET/2pentanone.yld',
     $     STATUS='old')
      do i = 1,6
         read(kin,*)
      enddo

      n  = 4
      n1 = n
      n2 = n
      n3 = n
      DO i = 1, n
         READ(kin,*) x1(i), y1(i), y2(i), y3(i)
         x2(i) = x1(i)
         x3(i) = x1(i)
      ENDDO
      CLOSE(kin)

      CALL interpol(x1,y1,kdata,n1,nw,wl,jlabel(j-1),0.058,0.2,yg1)
      CALL interpol(x2,y2,kdata,n2,nw,wl,jlabel(j-1),0.007,0.025,yg2)
      CALL interpol(x3,y3,kdata,n3,nw,wl,jlabel(j-1),0.05,0.05,yg3)


* combine xs and qy:
      DO iw = 1, nw - 1
        sig = yg(iw)
        qy1 = yg1(iw)
        qy2 = yg2(iw)
        qy3 = yg3(iw)
        IF(wc(iw)<=254.) THEN
          qy4 = 0.39
         ELSEIF(wc(iw)>=313.) THEN
          qy4 = 0.28
         ELSE
          qy4 = 0.39 - (wc(iw)-254.) * (0.11/59.)
        ENDIF
        DO i = 1, nz
          IF(mabs==2 .AND. wc(iw)>=240. .AND. wc(iw)<=330.) THEN
            A   = 4.72e-20 + 4.87e-23*tlev(i)
            lc  = 274.16   + 0.0186  *tlev(i)
            w   =  25.1    + 0.0101  *tlev(i)
            sig = A * exp(-((wc(iw)-lc)/w)**2)
          ENDIF
          sq(j-3,i,iw) = sig * qy1
          sq(j-2,i,iw) = sig * qy2
          sq(j-1,i,iw) = sig * qy3
          sq(j  ,i,iw) = sig * qy4
        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE mk03(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! methyl n-butyl ketone

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for                =*
*=  C4H9COCH3 photolysis:                                                    =*
*=        C4H9COCH3 -> CH3CH=CH2 + CH2=C(OH)CH3                         =*
*=                                                                           =*
*=  Cross section:  McMillan (1966)                                          =*
*=  Quantum yield:  estimates based on Calvert et al. 2011 book              =*
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

      REAL yg(kw), yg0(kw), yg1(kw)
      REAL sig, qy, eta
      INTEGER iw

      INTEGER myld, qyvers(3)
      CHARACTER(llog) :: logmsg(nlog)
      CHARACTER(lspc) :: spc

***************** C4H9COCH3 photolysis *************************


* Setting photolysis index and quantum yield options

      j = j+1
      jlabel(j) = 'C4H9COCH3 -> CH3CH=CH2 + CH2=C(OH)CH3'

      spc = "C4H9COCH3"
      qyvers = (/1, 2, 1/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = 'without quenching from Calvert et al. 2011 book'
      logmsg(2) = 'with quenching from Calvert et al. 2011 book'
      CALL set_option(2,spc,"yld",logmsg,qyvers,myld)

* cross sections

      OPEN(UNIT=kin,FILE='DATAJ1/MCMext/KET/2hexanone.abs',
     $     STATUS='old')
      do i = 1,6
         read(kin,*)
      enddo

      n = 68
      DO i = 1, n
         READ(kin,*) x1(i), y1(i)
      ENDDO
      CLOSE(kin)

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)


      OPEN(UNIT=kin,FILE='DATAJ1/MCMext/KET/2hexanone.yld',
     $     STATUS='old')
      do i = 1,7
         read(kin,*)
      enddo

      n  = 72
      n1 = n
      n2 = n
      DO i = 1, n
        READ(kin,*) x1(i), y1(i), y2(i)
        x2(i) = x1(i)
      ENDDO
      CLOSE(kin)

      CALL interpol(x1,y1,kdata,n1,nw,wl,jlabel(j-1),0.,0.,yg0)
      CALL interpol(x2,y2,kdata,n2,nw,wl,jlabel(j-1),0.,0.,yg1)


* combine xs and qy:
      DO iw = 1, nw - 1
        sig = yg(iw)
        IF(myld == 1) THEN
          qy = yg0(iw)
         ELSEIF(myld == 2) THEN
          eta = MAX(0.,yg0(iw)/yg1(iw) - 1.)
        ENDIF
        DO i = 1, nz
          IF(yg1(iw)>0.) THEN
            qy = yg1(iw)*(1.+eta)/(1.+eta*airden(i)/2.465e19)
           ELSE
            qy = 0
          ENDIF
          sq(j  ,i,iw) = sig * qy
        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE mk04(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! ethyl n-propyl ketone

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for                =*
*=  C3H7COC2H5 photolysis:                                                   =*
*=        C3H7COC2H5 -> products                                        =*
*=                                                                           =*
*=  Cross section:  Horrowitz (1999) unpublished data                        =*
*=  Quantum yield:  estimates based on Calvert et al. 2011 book              =*
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
      PARAMETER(kdata=2500)

      INTEGER i, n
      REAL x1(kdata)
      REAL y1(kdata)

* local

      REAL yg(kw)
      REAL qyI, qy1, qy2, qy3, qy4
      REAL sig
      INTEGER iw

***************** C3H7COC2H5 photolysis *************************

* Setting photolysis index
      j = j+1
      jlabel(j) = 'C3H7COC2H5 -> C2H5CO + C3H7'
      j = j+1
      jlabel(j) = 'C3H7COC2H5 -> C3H7CO + C2H5'
      j = j+1
      jlabel(j) = 'C3H7COC2H5 -> C3H7 + CO + C2H5'
      j = j+1
      jlabel(j) = 'C3H7COC2H5 -> C2H5C(OH)=CH2 + CH2=CH2'



* cross sections

      OPEN(UNIT=kin,FILE='DATAJ1/MCMext/KET/3hexanone.abs',
     $     STATUS='old')
      do i = 1,7
         read(kin,*)
      enddo

      n = 2318
      DO i = 1, n
         READ(kin,*) x1(i), y1(i)
      ENDDO
      CLOSE(kin)

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)


* quantum yields
      qyI = 0.61
      qy4 = 0.21

* combine xs and qy:
      DO iw = 1, nw - 1
        sig = yg(iw)
        IF(wc(iw)<=290.) THEN
          qy3 = 0.45*qyI
          qy1 = (qyI-qy3)/13*8
          qy2 = (qyI-qy3)/13*5
         ELSE
          qy3 = 0.25*qyI
          qy1 = (qyI-qy3)/13*8
          qy2 = (qyI-qy3)/13*5
        ENDIF
        DO i = 1, nz
          sq(j-3,i,iw) = sig * qy1
          sq(j-2,i,iw) = sig * qy2
          sq(j-1,i,iw) = sig * qy3
          sq(j  ,i,iw) = sig * qy4
        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE mk05(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! MIPK

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for                =*
*=  MIPK photolysis:                                                         =*
*=        MIPK -> products                                              =*
*=                                                                           =*
*=  Cross section:  estimated same as MIBK                                   =*
*=  Quantum yield:  estimates based on Calvert et al. 2011 book              =*
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
      REAL qyI, qy1, qy2, qy3, qy4, qydum
      REAL sig
      INTEGER iw

      INTEGER mabs, xsvers(3)
      CHARACTER(llog) :: logmsg(nlog)
      CHARACTER(lspc) :: spc

***************** MIPK photolysis *************************


* Setting photolysis index and cross section options
      j = j+1
      jlabel(j) = 'MIPK -> CH3CO + i-C3H7'
      j = j+1
      jlabel(j) = 'MIPK -> i-C3H7CO + CH3'
      j = j+1
      jlabel(j) = 'MIPK -> i-C3H7 + CO + CH3'
      j = j+1
      jlabel(j) = 'MIPK -> CH2=CHOH + CH3CH=CH2'

      spc = "MIPK"
      xsvers = (/1, 2, 1/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "same as MIBK"
      logmsg(2) = "same as DIPK"
      CALL set_option(2,spc,"abs",logmsg,xsvers,mabs)


* cross sections

      IF(mabs == 1) THEN
        OPEN(UNIT=kin,FILE='DATAJ1/MCMext/KET/MIBK.abs',
     $       STATUS='old')
        do i = 1,6
          read(kin,*)
        enddo

        n = 111
        DO i = 1, n
           READ(kin,*) x1(i), y1(i)
        ENDDO
        CLOSE(kin)

      ELSEIF(mabs == 2) THEN
        OPEN(UNIT=kin,FILE='DATAJ1/MCMext/KET/di-isopropylket.abs',
     $       STATUS='old')
        do i = 1,6
          read(kin,*)
        enddo

        n = 111
        DO i = 1, n
          READ(kin,*) x1(i), y1(i)
        ENDDO
        CLOSE(kin)
      ENDIF

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)


* combine xs and qy:
      DO iw = 1, nw - 1
        sig = yg(iw)
        DO i = 1, nz
          CALL qyacet(wc(iw),tlev(i),airden(i),qy3,qydum)
          qy3 = min(1., max(0.,qy3))
          IF(wc(iw)<=254.) THEN
            qy4 = 0.36
           ELSEIF(wc(iw)>313.) THEN
            qy4 = 0.0
           ELSE
            qy4 = 0.36 - (wc(iw)-254.) * (0.35/59.)
          ENDIF
          qyI = max(0.,0.75 - qy3 - qy4)
          qy1 = qyI*6/7
          qy2 = qyI/7
          sq(j-3,i,iw) = sig * qy1
          sq(j-2,i,iw) = sig * qy2
          sq(j-1,i,iw) = sig * qy3
          sq(j  ,i,iw) = sig * qy4
        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE mk06(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! MIBK

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for                =*
*=  MIBK photolysis:                                                         =*
*=        MIBK -> products                                              =*
*=                                                                           =*
*=  Cross section:  Yujing and Mellouki 2000                                 =*
*=  Quantum yield:  estimates based on Calvert et al. 2011 book              =*
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
      REAL qy1, qy2, qy3, qy4
      REAL sig
      INTEGER iw

      INTEGER myld, qyvers(3)
      CHARACTER(llog) :: logmsg(nlog)
      CHARACTER(lspc) :: spc

***************** MIBK photolysis *************************


* Setting photolysis index and quantum yield options
      j = j+1
      jlabel(j) = 'MIBK -> CH3CO + i-C4H9'
      j = j+1
      jlabel(j) = 'MIBK -> i-C4H9CO + CH3'
      j = j+1
      jlabel(j) = 'MIBK -> i-C4H9 + CO + CH3'
      j = j+1
      jlabel(j) = 'MIBK -> CH3C(OH)=CH2 + CH3CH=CH2'

      spc = "MIBK"
      qyvers = (/1, 4, 2/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = 'estimates (0.34*0.3/0.7) from GECKO-A database'
      WRITE(logmsg(2),"(2A)") 'estimates (0.15/0.35) based on ',
     &                        'Calvert et al. 2011 book'
      WRITE(logmsg(3),"(2A)") 'estimates (0.15/0.21) based on ',
     &                        'Calvert et al. 2011 book'
      logmsg(4) = 'scaled externally'
      CALL set_option(4,spc,"yld",logmsg,qyvers,myld)


* cross sections

      OPEN(UNIT=kin,FILE='DATAJ1/MCMext/KET/MIBK.abs',
     $     STATUS='old')
      do i = 1,6
         read(kin,*)
      enddo

      n = 111
      DO i = 1, n
         READ(kin,*) x1(i), y1(i)
      ENDDO
      CLOSE(kin)

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)


* quantum yields

      IF(myld==1) THEN
        qy1 = 0.34*0.3
        qy2 = 0.
        qy3 = 0.
        qy4 = 0.34*0.7
       ELSEIF(myld==2) THEN
        qy1 = 0.035
        qy2 = 0.035
        qy3 = 0.08
        qy4 = 0.35
       ELSEIF(myld==3) THEN
        qy1 = 0.035
        qy2 = 0.035
        qy3 = 0.08
        qy4 = 0.21
       ELSEIF(myld==4) THEN
        qy1 = 1.
        qy2 = 1.
        qy3 = 1.
        qy4 = 1.
      ENDIF


* combine xs and qy:
      DO iw = 1, nw - 1
        sig = yg(iw)
        DO i = 1, nz
          sq(j-3,i,iw) = sig * qy1
          sq(j-2,i,iw) = sig * qy2
          sq(j-1,i,iw) = sig * qy3
          sq(j  ,i,iw) = sig * qy4
        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE mk07(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! 4-Me-2-hexanone

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for                =*
*=  4-Me-2-hexanone photolysis:                                              =*
*=        4-Me-2-hexanone -> Norish type II products                    =*
*=                                                                           =*
*=  Cross section:  estimated same as 5-Me-2-hexanone                        =*
*=  Quantum yield:  estimates based on Calvert et al. 2011 book and MIPK     =*
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
      REAL qy, qy1, qy2, qyrat
      REAL sig
      INTEGER iw

      INTEGER mabs, myld, xsvers(3), qyvers(3)
      CHARACTER(llog) :: logmsg(nlog)
      CHARACTER(lspc) :: spc

***************** 4-Me-2-hexanone photolysis *************************


* Setting photolysis index and cross section/quantum yield options

! channels I - III ignored: qy(I-III)/qy(IVa/b) = 0.02 @ 42˚C
!     j = j+1
!     jlabel(j) = '4-Me-2-hexanone -> CH3CO + CH2CH(CH3)CH2CH3'
!     j = j+1
!     jlabel(j) = '4-Me-2-hexanone -> CH3CH2CH(CH3)CH2CO + CH3'
!     j = j+1
!     jlabel(j) = '4-Me-2-hexanone -> CH2CH(CH3)CH2CH3 + CO + CH3'
! explicit channels:
      j = j+1
      jlabel(j) = '4-Me-2-hexanone -> CH3C(OH)=CH2 + 2-butene'
      j = j+1
      jlabel(j) = '4-Me-2-hexanone -> CH3C(OH)=CH2 + 1-butene'

      spc = "4-Me-2-hexanone"
      xsvers = (/2, 1, 2/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "same as MIBK"
      logmsg(2) = "same as 5-Me-2-hexanone"
      CALL set_option(2,spc,"abs",logmsg,xsvers,mabs)

      qyvers = (/1, 2, 2/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = 'estimates (0.34*0.7) from GECKO-A database'
      WRITE(logmsg(2),"(A)") 'qy = 0.35 (estimate)'
      WRITE(logmsg(3),"(A)") 'qy = 0.21 (estimate)'
      CALL set_option(3,spc,"yld",logmsg,qyvers,myld)


* cross sections

      IF(mabs == 1) THEN
        OPEN(UNIT=kin,FILE='DATAJ1/MCMext/KET/MIBK.abs',
     $       STATUS='old')
        do i = 1,6
          read(kin,*)
        enddo

        n = 111
        DO i = 1, n
           READ(kin,*) x1(i), y1(i)
        ENDDO
        CLOSE(kin)

      ELSEIF(mabs == 2) THEN
        OPEN(UNIT=kin,FILE='DATAJ1/MCMext/KET/5Me2hexanone.abs',
     $       STATUS='old')
        do i = 1,6
           read(kin,*)
        enddo

        n = 111
        DO i = 1, n
           READ(kin,*) x1(i), y1(i)
        ENDDO
        CLOSE(kin)
      ENDIF

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)


* quantum yields

      IF(myld==1) THEN
        qy = 0.34*0.7
       ELSEIF(myld==2) THEN
        qy = 0.35
       ELSEIF(myld==3) THEN
        qy = 0.21
      ENDIF


* combine xs and qy:
      DO iw = 1, nw - 1
        sig = yg(iw)
        DO i = 1, nz
          qyrat = (-2.11e-3*wc(iw) + 0.7745)
     &          * 10**(171.92 * (1/317. - 1/tlev(i)))
          IF (wl(iw) <= 360.) THEN
            qy1 = qy / (1. + qyrat)
            qy2 = qy / (1. + 1./qyrat)
           ELSE
            qy1 = 0.
            qy2 = 0.
          ENDIF
          sq(j-1,i,iw) = sig * qy1
          sq(j  ,i,iw) = sig * qy2
        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE mk08(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! 5-Me-2-hexanone

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for                =*
*=  5-Me-2-hexanone photolysis:                                              =*
*=        5-Me-2-hexanone -> products                                   =*
*=                                                                           =*
*=  Cross section:  Yujing and Mellouki 2000                                 =*
*=  Quantum yield:  estimates based on Calvert et al. 2011 book              =*
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
      REAL qyI, qy1, qy2, qy3, qy4, qydum
      REAL sig
      INTEGER iw

      INTEGER myld, qyvers(3)
      CHARACTER(llog) :: logmsg(nlog)
      CHARACTER(lspc) :: spc

***************** 5-Me-2-hexanone photolysis *************************


* Setting photolysis index and quantum yield options

      j = j+1
      jlabel(j) = '5-Me-2-hexanone -> CH3CO + CH2CH2CH(CH3)2'
      j = j+1
      jlabel(j) = '5-Me-2-hexanone -> (CH3)2CHCH2CH2CO + CH3'
      j = j+1
      jlabel(j) = '5-Me-2-hexanone -> CH2CH2CH(CH3)2 + CO + CH3'
      j = j+1
      jlabel(j) = '5-Me-2-hexanone -> CH3C(OH)=CH2 + CH2=C(CH3)2'

      spc = "5-Me-2-hexanone"
      qyvers = (/1, 5, 3/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = 'estimates (0.34*0.7/0.3) from GECKO-A database'
      logmsg(2) = 'estimated same as MIPK'
      logmsg(3) = 'estimated same as MIBK'
      logmsg(4) = 'estimated same as 4-Me-2-hexanone'
      logmsg(5) = 'scaled externally'
      CALL set_option(5,spc,"yld",logmsg,qyvers,myld)


* cross sections

      OPEN(UNIT=kin,FILE='DATAJ1/MCMext/KET/5Me2hexanone.abs',
     $     STATUS='old')
      do i = 1,6
        read(kin,*)
      enddo

      n = 111
      DO i = 1, n
         READ(kin,*) x1(i), y1(i)
      ENDDO
      CLOSE(kin)

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)

* wavelength-independent quantum yields
      ylc: IF(myld==1) THEN
        qy1 = 0.34*0.3
        qy2 = 0.
        qy3 = 0.
        qy4 = 0.34*0.7
      ELSEIF(myld==3) THEN ylc
        qy1 = 0.035
        qy2 = 0.035
        qy3 = 0.08
        qy4 = 0.35
      ELSEIF(myld==5) THEN ylc
        qy1 = 1.
        qy2 = 1.
        qy3 = 1.
        qy4 = 1.
      ENDIF ylc


* combine xs and qy:
      DO iw = 1, nw - 1
        sig = yg(iw)
        DO i = 1, nz
          yl: IF(myld==2 .OR. myld==4) THEN
            IF(wc(iw)<=254.) THEN
              qy4 = 0.36
             ELSEIF(wc(iw)>313.) THEN
              qy4 = 0.0
             ELSE
              qy4 = 0.36 - (wc(iw)-254.) * (0.35/59.)
            ENDIF
            IF(myld==2) THEN
              CALL qyacet(wc(iw),tlev(i),airden(i),qy3,qydum)
              qy3 = min(1., max(0.,qy3))
              qyI = max(0.,0.75 - qy3 - qy4)
              qy2 = qyI/7
              qy1 = qyI*6/7
             ELSEIF(myld==4) THEN
              qy3 = 0.
              qy2 = 0.
              qy1 = 0.
            ENDIF
          ENDIF yl

          sq(j-3,i,iw) = sig * qy1
          sq(j-2,i,iw) = sig * qy2
          sq(j-1,i,iw) = sig * qy3
          sq(j  ,i,iw) = sig * qy4

        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE mk09(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! di-isopropyl ketone

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for                =*
*=  CH3CH(CH3)COCH(CH3)2 photolysis:                                         =*
*=        CH3CH(CH3)COCH(CH3)2 -> Norish type I products                =*
*=                                                                           =*
*=  Cross section:  Yujing and Mellouki (2000)                               =*
*=  Quantum yield:  estimates based on Calvert et al. 2011 book              =*
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

      INTEGER myld, qyvers(3)
      CHARACTER(llog) :: logmsg(nlog)
      CHARACTER(lspc) :: spc

***************** CH3CH(CH3)COCH(CH3)2 photolysis *************************


* Setting photolysis index and quantum yield options

      j = j+1
      jlabel(j) = 'CH3CH(CH3)COCH(CH3)2 -> i-C3H7CO + i-C3H7'
      j = j+1
      jlabel(j) = 'CH3CH(CH3)COCH(CH3)2 -> 2 i-C3H7 + CO'

      spc = "CH3CH(CH3)COCH(CH3)2"
      qyvers = (/1, 2, 1/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = 'estimate by Calvert et al. 2011'
      logmsg(2) = 'scaled externally'
      CALL set_option(2,spc,"yld",logmsg,qyvers,myld)


* cross sections

      OPEN(UNIT=kin,FILE='DATAJ1/MCMext/KET/di-isopropylket.abs',
     $     STATUS='old')
      do i = 1,6
         read(kin,*)
      enddo

      n = 111
      DO i = 1, n
         READ(kin,*) x1(i), y1(i)
      ENDDO
      CLOSE(kin)

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)


* quantum yields

      IF(myld==1) THEN
        qy1 = 0.5
        qy2 = 0.5
      ELSEIF(myld==2) THEN
        qy1 = 1.0
        qy2 = 1.0
      ENDIF


* combine xs and qy:
      DO iw = 1, nw - 1
        sig = yg(iw)
        DO i = 1, nz
          sq(j-1,i,iw) = sig * qy1
          sq(j  ,i,iw) = sig * qy2
        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE mk10(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! c-C4H6O

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for                =*
*=  cyclobutanone photolysis:                                                =*
*=        c-C4H6O ->  products                                          =*
*=                                                                           =*
*=  Cross section:  Calvert et al. 2008/2011 book                            =*
*=  Quantum yield:  estimates based on Calvert et al. 2011 book              =*
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

      INTEGER i, n, n1, n2, n3
      REAL x1(kdata),x2(kdata),x3(kdata)
      REAL y1(kdata),y2(kdata),y3(kdata)

* local

      REAL yg(kw), yg1(kw), yg2(kw)
      REAL qy1, qy2, qy3, qyII, qyrat
      REAL sig
      INTEGER iw

***************** c-C4H6O photolysis *************************


      j = j+1
      jlabel(j) = 'c-C4H6O -> C2H4 + CH2=C=O'
      j = j+1
      jlabel(j) = 'c-C4H6O -> C3H6 + CO'
      j = j+1
      jlabel(j) = 'c-C4H6O -> c-C3H6 + CO'
      j = j+1
      jlabel(j) = 'genC4cKet'


* cross sections

      OPEN(UNIT=kin,FILE='DATAJ1/MCMext/KET/cC4H6O.dat',
     $     STATUS='old')
      do i = 1,9
         read(kin,*)
      enddo

      n  = 63
      n1 = n
      n2 = n
      n3 = n
      DO i = 1, n
         READ(kin,*) x1(i), y1(i), y2(i), y3(i)
         y1(i) = y1(i)*1e-20
         x2(i) = x1(i)
         x3(i) = x1(i)
      ENDDO
      CLOSE(kin)

      CALL interpol(x1,y1,kdata,n1,nw,wl,jlabel(j),0.,0.,yg)


* quantum yields

      CALL interpol(x2,y2,kdata,n2,nw,wl,jlabel(j-2),y2(1),y2(n2),yg1)
      CALL interpol(x3,y3,kdata,n3,nw,wl,jlabel(j-1),y3(1),y3(n3),yg2)


* combine xs and qy:
      DO iw   = 1, nw - 1
        sig   = yg(iw)
        qy1   = yg1(iw)
        qyII  = yg2(iw)
        qyrat = 3.77e4*0.96**wc(iw)
        qy2   = qyII/(1.+1./qyrat)
        qy3   = qyII/(1.+ qyrat)
        DO i = 1, nz
          sq(j-3,i,iw) = sig * qy1
          sq(j-2,i,iw) = sig * qy2
          sq(j-1,i,iw) = sig * qy3
          sq(j  ,i,iw) = sig
        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE mk11(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! c-C6H10O

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for                =*
*=  cyclohexanone photolysis:                                                =*
*=        c-C6H10O ->  products                                         =*
*=                                                                           =*
*=  Cross section:  Iwasaki et al. 2008                                      =*
*=  Quantum yield:  estimates based on Calvert et al. 2011 book              =*
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

      REAL yg(kw), yg1(kw), yg2(kw), yg3(kw)
      REAL qy1, qy2, qy3
      REAL sig
      INTEGER iw

***************** c-C6H10O photolysis *************************


      j = j+1
      jlabel(j) = 'c-C6H10O -> 5-hexenal'
      j = j+1
      jlabel(j) = 'c-C6H10O -> cyclopentane + CO'
      j = j+1
      jlabel(j) = 'c-C6H10O -> 1-pentene + CO'
      j = j+1
      jlabel(j) = 'genC6cKet'


* cross sections

      OPEN(UNIT=kin,FILE='DATAJ1/MCMext/KET/cC6H10O.abs',
     $     STATUS='old')
      do i = 1,6
         read(kin,*)
      enddo

      n  = 23
      DO i = 1, n
         READ(kin,*) x1(i), y1(i)
      ENDDO
      CLOSE(kin)

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)


* quantum yields

      OPEN(UNIT=kin,FILE='DATAJ1/MCMext/KET/cC6H10O_1.yld',
     $     STATUS='old')
      do i = 1, 9
         read(kin,*)
      enddo

      n  = 59
      DO i = 1, n
         READ(kin,*) x1(i), y1(i)
      ENDDO
      CLOSE(kin)
      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j-3),0.,y1(n),yg1)

      OPEN(UNIT=kin,FILE='DATAJ1/MCMext/KET/cC6H10O_2.yld',
     $     STATUS='old')
      do i = 1, 9
         read(kin,*)
      enddo

      n  = 48
      DO i = 1, n
         READ(kin,*) x1(i), y1(i)
      ENDDO
      CLOSE(kin)
      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j-2),y1(1),0.,yg2)

      OPEN(UNIT=kin,FILE='DATAJ1/MCMext/KET/cC6H10O_3.yld',
     $     STATUS='old')
      do i = 1, 9
         read(kin,*)
      enddo

      n  = 67
      DO i = 1, n
         READ(kin,*) x1(i), y1(i)
      ENDDO
      CLOSE(kin)
      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j-1),0.,0.,yg3)


* combine xs and qy:
      DO iw = 1, nw - 1
        sig = yg(iw)
        qy1 = yg1(iw)
        qy2 = yg2(iw)
        qy3 = yg3(iw)
        IF(wc(iw)<=245.) THEN
          qy3 = max(0.,9.77e-3*wc(iw)-1.9646)
        ENDIF
        DO i = 1, nz
          sq(j-3,i,iw) = sig * qy1
          sq(j-2,i,iw) = sig * qy2
          sq(j-1,i,iw) = sig * qy3
          sq(j  ,i,iw) = sig
        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE mk12(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! c-C5H8O

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for                =*
*=  cyclopentanone photolysis:                                               =*
*=        c-C5H8O ->  products                                          =*
*=                                                                           =*
*=  Cross section:  Calvert et al. 2008/2011 book                            =*
*=  Quantum yield:  estimates based on Calvert et al. 2011 book              =*
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
      REAL qy1, qy2, qy3, qyI, qyrat
      REAL sig
      INTEGER iw

***************** c-C5H8O photolysis *************************


      j = j+1
      jlabel(j) = 'c-C5H8O -> 2 C2H4 + CO'
      j = j+1
      jlabel(j) = 'c-C5H8O -> c-C4H8 + CO'
      j = j+1
      jlabel(j) = 'c-C5H8O -> CH2=CHCH2CH2CHO'
      j = j+1
      jlabel(j) = 'genC5cKet'


* cross sections

      OPEN(UNIT=kin,FILE='DATAJ1/MCMext/KET/cC5H8O.abs',
     $     STATUS='old')
      do i = 1,8
         read(kin,*)
      enddo

      n  = 255
      DO i = 1, n
         READ(kin,*) x1(i), y1(i)
         y1(i) = y1(i)*1e-20
      ENDDO
      CLOSE(kin)

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j-1),0.,0.,yg)


* combine xs and qy:
      DO iw   = 1, nw - 1
        sig   = yg(iw)
        IF(wc(iw)<=280.) THEN
          qyI = 0.9
          qy3 = 0.
         ELSEIF(wc(iw)>=330.) THEN
          qyI = 0.05
          qy3 = 0.85
         ELSE
          qyI = 0.9 - (wc(iw)-290.)*0.85/36.
          qy3 = (wc(iw)-290.)*0.85/36.
        ENDIF
        IF(wc(iw)>=240.) THEN
          qyrat = 0.64
         ELSE
          qyrat = 0.56
        ENDIF
        qy1   = qyI/(1.+ qyrat)
        qy2   = qyI/(1.+1./qyrat)
        DO i = 1, nz
          sq(j-3,i,iw) = sig * qy1
          sq(j-2,i,iw) = sig * qy2
          sq(j-1,i,iw) = sig * qy3
          sq(j  ,i,iw) = sig
        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE mk13(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! c-C3H4O

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for                =*
*=  cyclopropanone photolysis:                                               =*
*=        c-C3H4O ->  products                                          =*
*=                                                                           =*
*=  Cross section:  Thomas and Rodiguez (1971)                               =*
*=  Quantum yield:  estimates based on Calvert et al. 2011 book              =*
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

      REAL yg(kw), yg1(kw)
      REAL qy1, qy2, qy
      REAL sig
      INTEGER iw

***************** c-C3H4O photolysis *************************


      j = j+1
      jlabel(j) = 'c-C3H4O -> C2H4 + CO'
      j = j+1
      jlabel(j) = 'c-C3H4O -> further products'
      j = j+1
      jlabel(j) = 'genC3cKet'


* cross sections

      OPEN(UNIT=kin,FILE='DATAJ1/MCMext/KET/cC3H4O_calv.dat',
     $     STATUS='old')
      do i = 1,9
         read(kin,*)
      enddo

      n  = 122
      n1 = n
      n2 = n
      DO i = 1, n
         READ(kin,*) x1(i), y1(i), y2(i)
         y1(i) = y1(i)*1e-20
         x2(i) = x1(i)
      ENDDO
      CLOSE(kin)

      CALL interpol(x1,y1,kdata,n1,nw,wl,jlabel(j),0.,0.,yg)


* quantum yields

      CALL interpol(x2,y2,kdata,n2,nw,wl,jlabel(j),1.,0.,yg1)


* combine xs and qy:
      DO iw = 1, nw - 1
        sig = yg(iw)
        IF(wc(iw)<=320.) THEN
          qy = 1.
         ELSEIF(wc(iw)>=370.) THEN
          qy = 0.9
         ELSE
          qy = 1. - (wc(iw)-320.)*0.1/50.
        ENDIF
        qy1 = yg1(iw)
        qy2 = qy - qy1
        DO i = 1, nz
          sq(j-2,i,iw) = sig * qy1
          sq(j-1,i,iw) = sig * qy2
          sq(j  ,i,iw) = sig
        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE mk14(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! EVK

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for CH3CH2COCH=CH  =*
*=  Ethyl vinyl ketone photolysis:                                           =*
*=           CH3CH2COCH=CH2 -> Products                                 =*
*=                                                                           =*
*=  Cross section same as MVK from                                           =*
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
      PARAMETER(kdata=50)

      INTEGER i, n
      REAL x(kdata), y(kdata)

* local

      REAL yg(kw)
      REAL qy, qy1, qy2, qy3
      INTEGER iw

      INTEGER myld, qyvers(3)
      CHARACTER(llog) :: logmsg(nlog)
      CHARACTER(lspc) :: spc

***************** ethyl vinyl ketone photolysis *************************


      j = j+1
      jlabel(j) = 'CH3CH2COCH=CH2 -> C2H5 + C2H3CO'
      j = j+1
      jlabel(j) = 'CH3CH2COCH=CH2 -> C2H3 + C2H5CO'
      j = j+1
      jlabel(j) = 'CH3CH2COCH=CH2 -> 1-C4H8 + CO'
      j = j+1
      jlabel(j) = 'genuKet(poly)'


* quantum yield options
      spc = "EVK"
      qyvers = (/1, 2, 2/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = 'same as MVK'
      logmsg(2) = 'same as MVK with branching scaled externally'
      CALL set_option(2,spc,"yld",logmsg,qyvers,myld)

* cross sections
      OPEN(UNIT=kin,FILE='DATAJ1/MCMext/KET/EVK.abs',
     $     STATUS='old')
      DO i = 1, 4
        READ(kin,*)
      ENDDO
      n = 24
      DO i = 1, n
        READ(kin,*) x(i), y(i)
      ENDDO
      CLOSE(kin)

      CALL interpol(x,y,kdata,n,nw,wl,jlabel(j),0.,0.,yg)

* quantum yield estimated same as MVK from
* Gierczak, T., J. B. Burkholder, R. K. Talukdar, A. Mellouki, S. B. Barone,
* and A. R. Ravishankara, Atmospheric fate of methyl vinyl ketone and methacrolein,
* J. Photochem. Photobiol A: Chemistry, 110 1-10, 1997.
* depends on pressure and wavelength, set upper limit to 1.0
* also recommended by IUPAC

* branching from IUPAC

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
            sq(j-3,i,iw) = yg(iw) * qy1
            sq(j-2,i,iw) = yg(iw) * qy2
            sq(j-1,i,iw) = yg(iw) * qy3
            sq(j  ,i,iw) = yg(iw)
         ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE mk15(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! methyl 2-hydroxyethyl ketone

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for                =*
*=  methyl 2-hydroxethyl ketone (CH3COC2H4OH) photolysis:                    =*
*=           CH3COC2H4OH  -> 20% Norrish type I + 80% Norrish type II   =*
*=                                                                           =*
*=  Cross section:  Messaadia et al. (2012)                                  =*
*=  Quantum yield:  0.08 Bouzidi et al., ES&T, 2015 (0.11 general estimate   =*
*=                  for ketones with beta-OH on secondary C)                 =*
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

      INTEGER myld, qyvers(3)
      CHARACTER(llog) :: logmsg(nlog)
      CHARACTER(lspc) :: spc


***************** CH3COC2H4OH photolysis *************************


* Setting photolysis index and quantum yield options

      j = j+1
      jlabel(j) = 'CH3COC2H4OH -> CH3 + COCH2CH2OH'
      j = j+1
      jlabel(j) = 'CH3COC2H4OH -> CH3COCH3 + HCHO'

      spc = "CH3COC2H4OH"
      qyvers = (/1, 2, 1/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = 'by Bouzidi et al. 2015'
      logmsg(2) = 'scaled externally'
      CALL set_option(2,spc,"yld",logmsg,qyvers,myld)


* cross sections
      OPEN(UNIT=kin,FILE='DATAJ1/MCMext/KET/4OH2butanone.abs',
     $     STATUS='old')
      DO i = 1, 4
        READ(kin,*)
      ENDDO
      n = 146
      DO i = 1, n
        READ(kin,*) x(i), y(i)
      ENDDO
      CLOSE(kin)

      CALL interpol(x,y,kdata,n,nw,wl,jlabel(j),0.,0.,yg)


* quantum yields

      IF(myld == 1) THEN
        qy1 = 0.2*0.08
        qy2 = 0.8*0.08
      ELSEIF(myld == 2) THEN
        qy1 = 1.0
        qy2 = 1.0
      ENDIF


* combine
      DO iw = 1, nw - 1
         DO i = 1, nz
            sq(j-1,i,iw) = yg(iw) * qy1
            sq(j  ,i,iw) = yg(iw) * qy2
         ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE mk16(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! acetoin

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for                =*
*=  acetoin (CH3COCH(OH)CH3) photolysis:                                      =*
*=           CH3COCH(OH)CH3  -> CH3CO + CH3CHOH                         =*
*=                                                                           =*
*=  Cross section:  Messaadia et al. (2012)                                  =*
*=  Quantum yield:  estimates based on hydroxyacetone, MEK, and 2-pentanone  =*
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
      INTEGER myld, qyvers(3)
      CHARACTER(llog) :: logmsg(nlog)
      CHARACTER(lspc) :: spc


***************** CH3COCH(OH)CH3 photolysis *************************


* Setting photolysis index and quantum yield options

      j = j+1
      jlabel(j) = 'CH3COCH(OH)CH3 -> CH3CO + CH3CHOH'
      j = j+1
      jlabel(j) = 'CH3COCH(OH)CH3 -> CH3CH(OH)CO + CH3'

      spc = "CH3COCH(OH)CH3"
      qyvers = (/1, 4, 2/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "same as hydroxyacetone"
      logmsg(2) = "same as MEK"
      WRITE(logmsg(3),"(2A)") "mean of methylacetoin, ",
     &     "4-hydroxy-2-butanone, and 4-hydroxy-4-methyl-2-pentanone"
      logmsg(4) = "scaled externally"
      CALL set_option(4,spc,"yld",logmsg,qyvers,myld)


* cross sections
      OPEN(UNIT=kin,FILE='DATAJ1/MCMext/KET/3OH2butanone.abs',
     $     STATUS='old')
      DO i = 1, 4
        READ(kin,*)
      ENDDO
      n = 133
      DO i = 1, n
        READ(kin,*) x(i), y(i)
      ENDDO
      CLOSE(kin)

      CALL interpol(x,y,kdata,n,nw,wl,jlabel(j),0.,0.,yg)


* quantum yields

      IF(myld==1) THEN
        qy1 = 0.31
        qy2 = 0.29
       ELSEIF(myld==2) THEN
        qy1 = 0.34
        qy2 = 0.00
       ELSEIF(myld==3) THEN
        qy1 = 0.00
        qy2 = 0.11
       ELSEIF(myld==4) THEN
        qy1 = 1.
        qy2 = 1.
      ENDIF

* combine
      DO iw = 1, nw - 1
         DO i = 1, nz
            sq(j-1,i,iw) = yg(iw) * qy1
            sq(j  ,i,iw) = yg(iw) * qy2
         ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE mk17(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! CH3COC(CH3)2OH -> CH3CO + (CH3)2COH

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for                =*
*=  methylacetoin (CH3COC(CH3)2OH) photolysis:                               =*
*=           CH3COC(CH3)2(OH)  -> CH3CO + (CH3)2COH                     =*
*=                                                                           =*
*=  Cross section:  Messaadia et al. (2012)                                  =*
*=  Quantum yield:  0.1 Bouzidi et al., Atmos Environ, 2014 (0.11 general    =*
*=                  estimate for ketones with alpha-OH on tert-C)            =*
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
      REAL qy
      INTEGER iw

      INTEGER myld, qyvers(3)
      CHARACTER(llog) :: logmsg(nlog)
      CHARACTER(lspc) :: spc


***************** CH3COC(CH3)2OH photolysis *************************


* Setting photolysis index and quantum yield options

      j = j+1
      jlabel(j) = 'CH3COC(CH3)2OH -> CH3 + (CH3)2C(OH)CO'

      spc = "CH3COC(CH3)2OH"
      qyvers = (/1, 2, 1/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "by Bouzidi et al. 2014"
      logmsg(2) = "scaled externally"
      CALL set_option(2,spc,"yld",logmsg,qyvers,myld)


* cross sections
      OPEN(UNIT=kin,FILE='DATAJ1/MCMext/KET/3OH3Me2butanone.abs',
     $     STATUS='old')
      DO i = 1, 4
        READ(kin,*)
      ENDDO
      n = 151
      DO i = 1, n
        READ(kin,*) x(i), y(i)
      ENDDO
      CLOSE(kin)

      CALL interpol(x,y,kdata,n,nw,wl,jlabel(j),0.,0.,yg)


* quantum yields
      IF(myld == 1) THEN
        qy = 0.1
      ELSEIF(myld == 2) THEN
        qy = 1.0
      ENDIF

* combine
      DO iw = 1, nw - 1
         DO i = 1, nz
            sq(j,i,iw) = yg(iw) * qy
         ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE mk18(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! 4-hydroxy-4-methyl-2-pentanone

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for                =*
*=  4-hydroxy-4-methyl-2-pentanone (CH3COCH2C(CH3)2OH) photolysis:           =*
*=     CH3COCH2C(CH3)2OH  -> 60% Norrish type I + 40% Norrish type II   =*
*=                                                                           =*
*=  Cross section:  Aslan et al., Atmos Environ, 2017                        =*
*=                  or Magneron et al., ES&T, 2003                           =*
*=  Quantum yield:  0.15 Aslan et al., Atmos Environ, 2015 (0.11 general     =*
*=                  estimate for ketones with beta-OH on tert-C)             =*
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


***************** CH3COCH2C(CH3)2OH photolysis *************************

* Setting photolysis index and cross section/quantum yield options

      j = j+1
      jlabel(j) = 'CH3COCH2C(CH3)2OH -> CH3COCH2 + CH3C(OH)CH3'
      j = j+1
      jlabel(j) = 'CH3COCH2C(CH3)2OH -> 2 CH3COCH3'

* cross section options
      spc = "CH3COCH2C(CH3)2OH"
      xsvers = (/2, 2, 1/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "from Azlan et al. 2017"
      logmsg(2) = "from Magneron et al. 2003"
      CALL set_option(2,spc,"abs",logmsg,xsvers,mabs)

* quantum yield options
      qyvers = (/1, 2, 1/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = "by Azlan et al. 2017"
      logmsg(2) = "scaled externally"
      CALL set_option(2,spc,"yld",logmsg,qyvers,myld)


* cross sections

      IF(mabs == 1) THEN
        OPEN(UNIT=kin,
     &       FILE='DATAJ1/MCMext/KET/4OH4Me2Pentanone_Asl17.abs',
     &       STATUS='old')
        DO i = 1, 4
          READ(kin,*)
        ENDDO
        n = 139
        DO i = 1, n
          READ(kin,*) x(i), y(i)
          y(i) = y(i)*1.e-20
        ENDDO
        CLOSE(kin)
        CALL interpol(x,y,kdata,n,nw,wl,jlabel(j),0.,0.,yg)

      ELSEIF(mabs == 2) THEN

        OPEN(UNIT=kin,
     &       FILE='DATAJ1/MCMext/KET/4OH4Me2Pentanone_Mag03.abs',
     &       STATUS='old')
        DO i = 1, 4
          READ(kin,*)
        ENDDO
        n = 49
        DO i = 1, n
          READ(kin,*) x(i), y(i)
          y(i) = y(i)*1.e-20
        ENDDO
        CLOSE(kin)

        CALL interpol(x,y,kdata,n,nw,wl,jlabel(j),0.,0.,yg)
      ENDIF


* quantum yields

      IF(myld == 1) THEN
        qy1 = 0.6*0.15
        qy2 = 0.4*0.15
      ELSEIF(myld == 2) THEN
        qy1 = 1.0
        qy2 = 1.0
      ENDIF


* combine
      DO iw = 1, nw - 1
         DO i = 1, nz
            sq(j-1,i,iw) = yg(iw) * qy1
            sq(j  ,i,iw) = yg(iw) * qy2
         ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE mk19(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! Ketene

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for                =*
*=  Ketene (CH2=C=O) photolysis:                                             =*
*=           CH2=C=O -> CO2 + CO + H2                                   =*
*=                                                                           =*
*=  Cross section:  Laufer and Keller (1971)                                 =*
*=  Quantum yield:  Calvert et al. (2011)                                    =*
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

      INTEGER i, n1, n2
      REAL x1(kdata),x2(kdata)
      REAL y1(kdata),y2(kdata)

* local

      REAL yg(kw),yg1(kw)
      INTEGER iw


***************** Ketene photolysis *************************

      j = j+1
      jlabel(j) = 'CH2=C=O -> CO2 + CO + H2'


* read data

      OPEN(UNIT=kin,FILE='DATAJ1/MCMext/KET/Ketene.dat',
     $     STATUS='old')
      DO i = 1, 9
        READ(kin,*)
      ENDDO
      n1 = 151
      n2 = n1
      DO i = 1, n1
        READ(kin,*) x1(i), y1(i), y2(i)
        x2(i) = x1(i)
        y1(i) = y1(i)*1.E-20
      ENDDO
      CLOSE(kin)


* cross sections
      CALL interpol(x1,y1,kdata,n1,nw,wl,jlabel(j),0.,0.,yg)

* quantum yields
      CALL interpol(x2,y2,kdata,n2,nw,wl,jlabel(j),y2(1),0.,yg1)


* combine
      DO iw = 1, nw - 1
        DO i = 1, nz
          sq(j,i,iw) = yg(iw) * yg1(iw)
        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE mk20(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! Methylketene

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for                =*
*=  Methylketene (CH3CH=C=O) photolysis:                                     =*
*=           CH3CH=C=O -> C2H4 + CO                                     =*
*=                                                                           =*
*=  Cross section:  Chong and Kistiakowsky (1964)                            =*
*=  Quantum yield:  Calvert et al. (2011)                                    =*
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

      INTEGER i, n1, n2
      REAL x1(kdata),x2(kdata)
      REAL y1(kdata),y2(kdata)

* local

      REAL yg(kw),yg1(kw),dum
      INTEGER iw

***************** Methylketene photolysis *************************

      j = j+1
      jlabel(j) = 'CH3CH=C=O -> C2H4 + CO'
      j = j+1
      jlabel(j) = 'genKete(poly)'


* cross sections

      OPEN(UNIT=kin,FILE='DATAJ1/MCMext/KET/methylketene.abs',
     $     STATUS='old')
      DO i = 1, 8
        READ(kin,*)
      ENDDO
      n1 = 136
      DO i = 1, n1
        READ(kin,*) x1(i), y1(i)
        y1(i) = y1(i)*1.E-20
      ENDDO
      CLOSE(kin)

      CALL interpol(x1,y1,kdata,n1,nw,wl,jlabel(j),0.,0.,yg)


* quantum yields

      OPEN(UNIT=kin,FILE='DATAJ1/MCMext/KET/Ketene.dat',
     $     STATUS='old')
      DO i = 1, 9
        READ(kin,*)
      ENDDO
      n2 = 151
      DO i = 1, n2
        READ(kin,*) x2(i), dum, y2(i)
      ENDDO
      CLOSE(kin)

      CALL interpol(x2,y2,kdata,n2,nw,wl,jlabel(j),y2(1),0.,yg1)


* combine
      DO iw = 1, nw - 1
         DO i = 1, nz
            sq(j-1,i,iw) = yg(iw) * yg1(iw)
            sq(j  ,i,iw) = yg(iw)
         ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE mk21(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! Generic unbranched or ext. cyclic ketones

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for                =*
*=  RCOR'  photolysis:                                                       =*
*=           RCOR'+ hv -> R' + RCO                                           =*
*=  and alpha-hydroxy substituted equivalents.                               =*
*=  Routine calculates j values for unbranched ketones with or without OH    =*
*=  substitutions or for ketones with "external" ring structures.            =*
*=                                                                           =*
*=  Cross section:  Parameterisations sig = A*exp(-((lambda-lambda_c)/w)**2) =*
*=  Quantum yield:  estimate (0.75) with equal branching for NorI and NorII  =*
*=                  for unbranched ketones and estimate (0.14) for ketones   =*
*=                  with external rings.                                     =*
*=                                                                           =*
*-----------------------------------------------------------------------------*
*=                                                                           =*
*= Only C3 compounds are computed explicitly. Higher unbranched ketones are  =*
*= fitted with fit functions:                                                =*
*= j(unsubst)  = j(C3)*(1.25*exp(CN/11.39)-0.55)                             =*
*= j(alpha-OH) = j(C3)*(0.73*exp(CN/9.74)-0.09)                              =*
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
*-----------------------------------------------------------------------------*'

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

      REAL    sig,sigOH,sigtOH,sig2OH !,qy,qy1,qy2,qyc
      REAL    A,lc,w,AOH,lcOH,wOH,AtOH,lctOH,wtOH,A2OH,lc2OH
      INTEGER iz,iw

***************** Generic unbranched or ext. cyclic ketones photolysis ***

      j = j+1
      jlabel(j) = 'lcKet -> products'
      j = j+1
      jlabel(j) = 'lcKet(OH) -> products'
      j = j+1
      jlabel(j) = 'lcKet(2OH) -> products'
      j = j+1
      jlabel(j) = 'cKet(t-OH) -> products'


* cross sections are parameterised with Gaussian fit
*
* sig (CN, wl, T) = A(CN,T)*exp{-[(wl - lc(CN,T)) / w(CN,T)]**2}
* variables:  CN = carbon number; wl wavelength; T = temperature
* parameters: A = amplitude; lc = centre wavelength; w = curve width
*
* with (A(CN, T) / 1E-20 cm2) = 4.845e-3 * (T/K) + 2.364*ln(CN) + 1.035
*      (lc(CN, T) / nm)       = 1.71e-2  * (T/K) + 1.718 * CN   + 265.344
*       w(CN, T)              = 9.8e-3   * (T/K) - 0.728 * CN   + 29.367
*
* for alpha-OH subst. ketones, parameters are adjusted as follows:
*
* (A(OH) / 1E-20 cm2) = A + 1.75
* (lc(OH) / nm)       = lc - 8.0
* w(OH)               = w
*
* if OH is on quaternary carbon atom, parameters change to:
*
* (A(t-OH) / 1E-20 cm2) = A + 1.75
* (lc(t-OH) / nm)       = lc - 8.0
* w(t-OH)               = w


* quantum yields
! reaction pathways (NI and NII channels) are now represented by
! scaling the final j value in the mechanism
C      qy  = 0.75*0.5
C      qy1 = 0.34
C      qy2 = 0.28
C      qyc = 0.14


* combine
      DO iw = 1, nw - 1
         DO iz = 1, nz
             A     = 4.845e-3 * tlev(iz) + 2.364*log(3.) + 1.035
             A     = A*1.e-20
             lc    = 1.71e-2   * tlev(iz) + 1.718 * 3.   + 265.344
             w     = 9.8e-3    * tlev(iz) - 0.728 * 3.   + 29.367
             AOH   = A + 1.75e-20
             A2OH   = A + 2*1.75e-20
             lcOH  = lc - 8.
             lc2OH  = lc - 16.
             wOH   = w
             AtOH   = A*0.86
             lctOH  = lc - 8.
             wtOH   = w
             sig   = A*exp(-((wc(iw) - lc) / w)**2)
             sigOH = AOH*exp(-((wc(iw) - lcOH) / wOH)**2)
             sig2OH = A2OH*exp(-((wc(iw) - lc2OH) / wOH)**2)
             sigtOH = AtOH*exp(-((wc(iw) - lctOH) / wtOH)**2)
             sq(j-3,iz,iw) = sig
             sq(j-2,iz,iw) = sigOH
             sq(j-1,iz,iw) = sig2OH
             sq(j  ,iz,iw) = sigtOH
         ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE mk22(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! alpha-branched ketones (with Norrish II and OH-subst.)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for                =*
*=  alpha-branched ketone photolysis:                                        =*
*=                                                                           =*
*=  Cross section:  parameterised DIPK cross sections (scaled for OH-subst.) =*
*=  Quantum yield:  estimate                                                 =*
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

      INTEGER i
      REAL    AOH,lcOH,wOH,AtOH,lctOH,wtOH
      REAL    sigOH, sigtOH
C      REAL    qy1, qy2, qy1OH, qy2OH
      INTEGER iw

***************** alpha-branched ketones (with Norrish II and OH-subst.) ***

      j = j+1
      jlabel(j) = 'a-br. Ket(OH) -> products'

      j = j+1
      jlabel(j) = 'a-br. Ket(t-OH) -> products'


* cross sections are parameterised as given below

* quantum yields
! Now scaling j values externally with effective quantum yields
C      qy1   = 0.15
C      qy2   = 0.35
C      qy1OH = 0.34
C      qy2OH = 0.28



* combine xs and qy:
      DO iw = 1, nw - 1
        AOH   = 9.65e-20
        lcOH  = 280.1
        wOH   = 25.4
        AtOH  = 6.79e-20
        lctOH = 280.1
        wtOH  = 25.4
        sigOH = AOH*exp(-(wc(iw) - lcOH)**2 / wOH**2)
        sigtOH = AtOH*exp(-(wc(iw) - lctOH)**2 / wtOH**2)
        DO i = 1, nz
          sq(j-1,i,iw) = sigOH
          sq(j  ,i,iw) = sigtOH
        ENDDO
      ENDDO

      END

* ============================================================================*

      SUBROUTINE mk23(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! OH-subst. bKet

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for                =*
*=  beta-branched ketone photolysis:                                         =*
*=                                                                           =*
*=  Cross section:  parameterised MIBK cross sections (scaled for OH-subst.) =*
*=  Quantum yield:  estimate                                                 =*
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

      INTEGER i
      REAL    A,lc,w,AOH,lcOH,wOH
      REAL    sigOH
!     REAL    qy, qy1, qy2, qy3
      INTEGER iw


***************** OH-subst. bKet photolysis *************************

      j = j+1
      jlabel(j) = 'b-br. Ket(OH) -> products'
!     j = j+1
!     jlabel(j) = 'b-br. Ket(OH) -> NII products'


* cross sections are parameterised as given below

* quantum yields
! Now externally scaled
!     qy1 = 0.34
!     qy2 = 0.28

* combine xs and qy:
      DO iw = 1, nw - 1
        A     = 6.38e-20
        lc    = 281.22
        w     = 29.16
        AOH   = A + 1.75e-20
        lcOH  = lc - 8.
        wOH   = w
        sigOH = AOH*exp(-((wc(iw) - lcOH) / wOH)**2)
        DO i = 1, nz
!         sq(j-1,i,iw) = sigOH * qy1
          sq(j  ,i,iw) = sigOH!* qy2
        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE mk24(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! OH-subst. uKet

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for                =*
*=  beta-branched ketone photolysis:                                         =*
*=                                                                           =*
*=  Cross section:  parameterised MIBK cross sections (scaled for OH-subst.) =*
*=  Quantum yield:  estimate                                                 =*
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

      INTEGER i
      REAL    A,lc,w,AOH,lcOH,wOH,AtOH,lctOH,wtOH
      REAL    sigOH,sigtOH
      REAL    qy, qy1, qy2, qy3
      INTEGER iw


***************** OH-subst. uKet photolysis *************************

      j = j+1
      jlabel(j) = 'uKet(OH) -> RdCO + ROH'
      j = j+1
      jlabel(j) = 'uKet(OH) -> ROHCO + Rd'
      j = j+1
      jlabel(j) = 'uKet(OH) -> alkene + CO'
      j = j+1
      jlabel(j) = 'uKet(t-OH) -> RdCO + ROH'
      j = j+1
      jlabel(j) = 'uKet(t-OH) -> ROHCO + Rd'
      j = j+1
      jlabel(j) = 'uKet(t-OH) -> alkene + CO'
      j = j+1
      jlabel(j) = 'genuKet(OHpoly)'
      j = j+1
      jlabel(j) = 'genuKet(t-OHpoly)'


* cross sections are parameterised as given below

* quantum yields


* combine xs and qy:
      DO iw = 1, nw - 1
        A     = 6.4e-20
        lc    = 324.8
        w     = 36.4
        AOH   = A + 1.75e-20
        lcOH  = lc - 8.
        wOH   = w
        AtOH  = A*0.86
        lctOH = lc - 8.
        wtOH  = w
        sigOH = AOH*exp(-((wc(iw) - lcOH) / wOH)**2)
        sigtOH = AtOH*exp(-((wc(iw) - lctOH) / wtOH)**2)
        DO i = 1, nz
          qy  = exp(-0.055*(wc(iw)-308.)) /
     $         (5.5 + 9.2e-19*airden(i))
          qy  = min(qy, 1.)
          qy1 = 0.2 * qy
          qy2 = 0.2 * qy
          qy3 = 0.6 * qy
          sq(j-7,i,iw) = sigOH * qy1
          sq(j-6,i,iw) = sigOH * qy2
          sq(j-5,i,iw) = sigOH * qy3
          sq(j-4,i,iw) = sigtOH * qy1
          sq(j-3,i,iw) = sigtOH * qy2
          sq(j-2,i,iw) = sigtOH * qy3
          sq(j-1,i,iw) = sigOH
          sq(j  ,i,iw) = sigtOH
        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE mk25(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! C5-Ketone

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for                =*
*=  C5 ketone photolysis:                                                    =*
*=        C5 ketone -> products                                         =*
*=                                                                           =*
*=  Cross section:  average of Martinez et al. 1992 for 2-pentanone          =*
*=                  and Koch et al. 2008 for 3-pentanone                     =*
*=  Quantum yield:  scaled externally                                        =*
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
      REAL sig
      INTEGER iw


***************** C5-Ketone photolysis *************************

      j = j+1
      jlabel(j) = 'lKET5 -> products'



* cross sections
      OPEN(UNIT=kin,FILE='DATAJ1/MCMext/KET/Ket5.abs',
     $     STATUS='old')
      do i = 1,9
        read(kin,*)
      enddo

      n  = 255
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
        y1(i) = y1(i) * 1.e-20
      ENDDO
      CLOSE(kin)

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)

* combine:
      DO iw = 1, nw - 1
        sig = yg(iw)
        DO i = 1, nz
          sq(j,i,iw) = sig
        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE mk26(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! c-C7H12O

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for                =*
*=  cycloheptanone photolysis:                                               =*
*=        c-C7H12O ->  products                                              =*
*=                                                                           =*
*=  Cross section:  based on Hamer and Huber (1978)                          =*
*=  Quantum yield:  based on Hamer and Huber (1978)                          =*
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

      INTEGER i, n, ns
      REAL x1(kdata), x2(kdata), x3(kdata)
      REAL y1(kdata), y2(kdata), y3(kdata)

* local

      REAL yg(kw), yg1(kw), yg2(kw), yg3(kw)
      REAL qy1, qy2, qy3
      REAL sig
      INTEGER iw

***************** c-C7H12O photolysis *************************

* Setting photolysis index and reaction names
      j = j+1
      jlabel(j) = 'c-C7H12O -> 6-heptenal'
      j = j+1
      jlabel(j) = 'c-C7H12O -> cyclohexane + CO'
      j = j+1
      jlabel(j) = 'c-C7H12O -> 1-hexene + CO'
      j = j+1
      jlabel(j) = 'genC7cKet'

* cross sections

      OPEN(UNIT=kin,FILE='DATAJ1/MCMext/KET/cHeptanone.abs',
     $     STATUS='old')
      do i = 1,8
         read(kin,*)
      enddo

      n  = 90
      DO i = 1, n
         READ(kin,*) x1(i), y1(i)
         y1(i) = y1(i) * 1.e-20
      ENDDO
      CLOSE(kin)

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)


* quantum yields

      OPEN(UNIT=kin,FILE='DATAJ1/MCMext/KET/cHeptanone.yld',
     $     STATUS='old')
      do i = 1, 7
         read(kin,*)
      enddo

      n  = 9
      ns = n
      DO i = 1, n
         READ(kin,*) x1(i), y1(i), y2(i), y3(i)
         x2(i) = x1(i)
         x3(i) = x1(i)
      ENDDO
      CLOSE(kin)
      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j-3),0.,y1(n),yg1)

      n  = ns
      CALL interpol(x2,y2,kdata,n,nw,wl,jlabel(j-2),0.,y2(n),yg2)

      n  = ns
      CALL interpol(x3,y3,kdata,n,nw,wl,jlabel(j-1),0.,y3(n),yg3)


* combine xs and qy:
      DO iw = 1, nw - 1
        sig = yg(iw)
        qy1 = yg1(iw)
        qy2 = yg2(iw)
        qy3 = yg3(iw)
        DO i = 1, nz
          sq(j-3,i,iw) = sig * qy1
          sq(j-2,i,iw) = sig * qy2
          sq(j-1,i,iw) = sig * qy3
          sq(j  ,i,iw) = sig
        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE mk27(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! 4-heptanone

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for                =*
*=  4-heptanone photolysis:                                                  =*
*=        4-C7H14O ->  products                                              =*
*=                                                                           =*
*=  Cross section:  C5 ketone average                                        =*
*=  Quantum yield:  based on Calvert et al. (2011 )                          =*
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

      INTEGER myld, qyvers(3)
      CHARACTER(llog) :: logmsg(nlog)
      CHARACTER(lspc) :: spc

***************** 4-heptanone photolysis *************************


* Setting photolysis index and quantum yield options
      j = j+1
      jlabel(j) = '4-heptanone -> NI products'
      j = j+1
      jlabel(j) = '4-heptanone -> NII products'

      spc = "4-heptanone"
      qyvers = (/1, 2, 1/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = 'estimates from Calvert et al. 2011 book'
      logmsg(2) = 'scaled externally'
      CALL set_option(2,spc,"yld",logmsg,qyvers,myld)


* cross sections
      OPEN(UNIT=kin,FILE='DATAJ1/MCMext/KET/Ket5.abs',
     $     STATUS='old')
      do i = 1,9
        read(kin,*)
      enddo

      n  = 255
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
        y1(i) = y1(i) * 1.e-20
      ENDDO
      CLOSE(kin)

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)


* quantum yields
      IF(myld==1) THEN
        qy1 = 0.34
        qy2 = 0.21
      ELSEIF(myld==2) THEN
        qy1 = 1.0
        qy2 = 1.0
      ENDIF


* combine xs and qy
      DO iw = 1, nw - 1
        sig = yg(iw)
        DO i = 1, nz
          sq(j-1,i,iw) = sig * qy1
          sq(j,  i,iw) = sig * qy2
        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE mk28(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! 4-octanone

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for                =*
*=  4-octanone photolysis:                                                   =*
*=        4-C8H16O ->  products                                              =*
*=                                                                           =*
*=  Cross section:  C5 ketone average                                        =*
*=  Quantum yield:  based on Calvert et al. (2011 )                          =*
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
      REAL qy1, qy2, qy3
      REAL sig
      INTEGER iw

      INTEGER myld, qyvers(3)
      CHARACTER(llog) :: logmsg(nlog)
      CHARACTER(lspc) :: spc

***************** 4-octanone photolysis *************************


* Setting photolysis index and quantum yield options
      j = j+1
      jlabel(j) = '4-octanone -> NI products'
      j = j+1
      jlabel(j) = '4-octanone -> n-C4H9COCH3 + C2H4'
      j = j+1
      jlabel(j) = '4-octanone -> n-C3H7COCH3 + C3H6'

      spc = "4-octanone"
      qyvers = (/1, 2, 1/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = 'estimates from Calvert et al. 2011 book'
      logmsg(2) = 'scaled externally'
      CALL set_option(2,spc,"yld",logmsg,qyvers,myld)


* cross sections
      OPEN(UNIT=kin,FILE='DATAJ1/MCMext/KET/Ket5.abs',
     $     STATUS='old')
      do i = 1,9
        read(kin,*)
      enddo

      n  = 255
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
        y1(i) = y1(i) * 1.e-20
      ENDDO
      CLOSE(kin)

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)


* quantum yields
      IF(myld==1) THEN
        qy1 = 0.078
        qy2 = 0.03
        qy3 = 0.054
      ELSEIF(myld==2) THEN
        qy1 = 1.0
        qy2 = 1.0
        qy3 = 1.0
      ENDIF


* combine xs and qy
      DO iw = 1, nw - 1
        sig = yg(iw)
        DO i = 1, nz
          sq(j-2,i,iw) = sig * qy1
          sq(j-1,i,iw) = sig * qy2
          sq(j,  i,iw) = sig * qy3
        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE mk29(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! n-propyl isopropyl ketone

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for                =*
*=  n-propyl isopropyl ketone photolysis:                                    =*
*=        n-C3H7COCH(CH3)2 ->  products                                      =*
*=                                                                           =*
*=  Cross section:  same as DIPK                                             =*
*=  Quantum yield:  based on Calvert et al. (2011 )                          =*
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

      INTEGER myld, qyvers(3)
      CHARACTER(llog) :: logmsg(nlog)
      CHARACTER(lspc) :: spc

***************** n-propyl isopropyl ketone photolysis *************************


* Setting photolysis index and quantum yield options
      j = j+1
      jlabel(j) = 'n-C3H7COCH(CH3)2 -> NI products'
      j = j+1
      jlabel(j) = 'n-C3H7COCH(CH3)2 -> i-C3H7COCH3 + C2H4'

      spc = "n-C3H7COCH(CH3)2"
      qyvers = (/1, 2, 1/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = 'estimates from Calvert et al. 2011 book'
      logmsg(2) = 'scaled externally'
      CALL set_option(2,spc,"yld",logmsg,qyvers,myld)


* cross sections

      OPEN(UNIT=kin,FILE='DATAJ1/MCMext/KET/di-isopropylket.abs',
     $     STATUS='old')
      do i = 1,6
         read(kin,*)
      enddo

      n = 111
      DO i = 1, n
         READ(kin,*) x1(i), y1(i)
      ENDDO
      CLOSE(kin)

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)


* quantum yields
      IF(myld==1) THEN
        qy1 = 0.69
        qy2 = 0.12
      ELSEIF(myld==2) THEN
        qy1 = 1.0
        qy2 = 1.0
      ENDIF


* combine xs and qy
      DO iw = 1, nw - 1
        sig = yg(iw)
        DO i = 1, nz
          sq(j-1,i,iw) = sig * qy1
          sq(j,  i,iw) = sig * qy2
        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE mk30(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! methyl neopentyl ketone

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for                =*
*=  methyl neopentyl ketone photolysis:                                      =*
*=        CH3COCH2C(CH3)3 ->  products                                       =*
*=                                                                           =*
*=  Cross section:  same as MIBK                                             =*
*=  Quantum yield:  based on Calvert et al. (2011 )                          =*
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

      INTEGER myld, qyvers(3)
      CHARACTER(llog) :: logmsg(nlog)
      CHARACTER(lspc) :: spc

***************** methyl neopentyl ketone photolysis *************************


* Setting photolysis index and quantum yield options
      j = j+1
      jlabel(j) = 'CH3COCH2C(CH3)3 -> NI products'
      j = j+1
      jlabel(j) = 'CH3COCH2C(CH3)3 -> CH3COCH3 + i-C4H8'

      spc = "CH3COCH2C(CH3)3"
      qyvers = (/1, 2, 1/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = 'estimates from Calvert et al. 2011 book'
      logmsg(2) = 'scaled externally'
      CALL set_option(2,spc,"yld",logmsg,qyvers,myld)


* cross sections

      OPEN(UNIT=kin,FILE='DATAJ1/MCMext/KET/MIBK.abs',
     $     STATUS='old')
      do i = 1,6
         read(kin,*)
      enddo

      n = 111
      DO i = 1, n
         READ(kin,*) x1(i), y1(i)
      ENDDO
      CLOSE(kin)

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)


* quantum yields
      IF(myld==1) THEN
        qy1 = 0.042
        qy2 = 0.27
      ELSEIF(myld==2) THEN
        qy1 = 1.0
        qy2 = 1.0
      ENDIF


* combine xs and qy
      DO iw = 1, nw - 1
        sig = yg(iw)
        DO i = 1, nz
          sq(j-1,i,iw) = sig * qy1
          sq(j,  i,iw) = sig * qy2
        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE mk31(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! n-propyl isobutyl ketone

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for                =*
*=  n-propyl isobutyl ketone photolysis:                                     =*
*=        n-propyl isobutyl ketone ->  products                              =*
*=                                                                           =*
*=  Cross section:  same as MIBK                                             =*
*=  Quantum yield:  based on Calvert et al. (2011 )                          =*
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
      REAL qy1, qy2, qy3
      REAL sig
      INTEGER iw

      INTEGER myld, qyvers(3)
      CHARACTER(llog) :: logmsg(nlog)
      CHARACTER(lspc) :: spc

***************** n-propyl isobutyl ketone photolysis *************************


* Setting photolysis index and quantum yield options
      j = j+1
      jlabel(j) = '2-Me-4-heptanone -> NI products'
      j = j+1
      jlabel(j) = '2-Me-4-heptanone -> i-C4H9COCH3 + C2H4'
      j = j+1
      jlabel(j) = '2-Me-4-heptanone -> n-C3H7COCH3 + C3H6'

      spc = "2-Me-4-heptanone"
      qyvers = (/1, 2, 1/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = 'estimates from Calvert et al. 2011 book'
      logmsg(2) = 'scaled externally'
      CALL set_option(2,spc,"yld",logmsg,qyvers,myld)


* cross sections

      OPEN(UNIT=kin,FILE='DATAJ1/MCMext/KET/MIBK.abs',
     $     STATUS='old')
      do i = 1,6
         read(kin,*)
      enddo

      n = 111
      DO i = 1, n
         READ(kin,*) x1(i), y1(i)
      ENDDO
      CLOSE(kin)

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)


* quantum yields
      IF(myld==1) THEN
        qy1 = 0.21
        qy2 = 0.055
        qy3 = 0.11
      ELSEIF(myld==2) THEN
        qy1 = 1.0
        qy2 = 1.0
        qy3 = 1.0
      ENDIF


* combine xs and qy
      DO iw = 1, nw - 1
        sig = yg(iw)
        DO i = 1, nz
          sq(j-2,i,iw) = sig * qy1
          sq(j-1,i,iw) = sig * qy2
          sq(j,  i,iw) = sig * qy3
        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE mk32(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! 3-Me-4-heptanone

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for                =*
*=  3-Me-4-heptanone photolysis:                                             =*
*=        3-Me-4-heptanone ->  products                                      =*
*=                                                                           =*
*=  Cross section:  same as DIPK                                             =*
*=  Quantum yield:  based on Calvert et al. (2011 )                          =*
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

      INTEGER myld, qyvers(3)
      CHARACTER(llog) :: logmsg(nlog)
      CHARACTER(lspc) :: spc

***************** 3-Me-4-heptanone photolysis *************************


* Setting photolysis index and quantum yield options
      j = j+1
      jlabel(j) = '3-Me-4-heptanone -> NI products'
      j = j+1
      jlabel(j) = '3-Me-4-heptanone -> NII products'

      spc = "3-Me-4-heptanone"
      qyvers = (/1, 2, 1/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = 'estimates from Calvert et al. 2011 book'
      logmsg(2) = 'scaled externally'
      CALL set_option(2,spc,"yld",logmsg,qyvers,myld)


* cross sections

      OPEN(UNIT=kin,FILE='DATAJ1/MCMext/KET/di-isopropylket.abs',
     $     STATUS='old')
      do i = 1,6
         read(kin,*)
      enddo

      n = 111
      DO i = 1, n
         READ(kin,*) x1(i), y1(i)
      ENDDO
      CLOSE(kin)

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)


* quantum yields
      IF(myld==1) THEN
        qy1 = 0.39
        qy2 = 0.1
      ELSEIF(myld==2) THEN
        qy1 = 1.0
        qy2 = 1.0
      ENDIF


* combine xs and qy
      DO iw = 1, nw - 1
        sig = yg(iw)
        DO i = 1, nz
          sq(j-1,i,iw) = sig * qy1
          sq(j,  i,iw) = sig * qy2
        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE mk33(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! 2,2-Me-3-hexanone

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for                =*
*=  2,2-Me-3-hexanone photolysis:                                            =*
*=        2,2-Me-3-hexanone ->  products                                     =*
*=                                                                           =*
*=  Cross section:  same as DIPK                                             =*
*=  Quantum yield:  based on Calvert et al. (2011 )                          =*
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

      INTEGER myld, qyvers(3)
      CHARACTER(llog) :: logmsg(nlog)
      CHARACTER(lspc) :: spc

***************** 2,2-Me-3-hexanone photolysis *************************


* Setting photolysis index and quantum yield options
      j = j+1
      jlabel(j) = '2,2-Me-3-hexanone -> NI products'
      j = j+1
      jlabel(j) = '2,2-Me-3-hexanone -> t-C4H9COCH3 + C2H4'

      spc = "2,2-Me-3-hexanone"
      qyvers = (/1, 2, 1/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = 'estimates from Calvert et al. 2011 book'
      logmsg(2) = 'scaled externally'
      CALL set_option(2,spc,"yld",logmsg,qyvers,myld)


* cross sections

      OPEN(UNIT=kin,FILE='DATAJ1/MCMext/KET/di-isopropylket.abs',
     $     STATUS='old')
      do i = 1,6
         read(kin,*)
      enddo

      n = 111
      DO i = 1, n
         READ(kin,*) x1(i), y1(i)
      ENDDO
      CLOSE(kin)

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)


* quantum yields
      IF(myld==1) THEN
        qy1 = 0.57
        qy2 = 0.097
      ELSEIF(myld==2) THEN
        qy1 = 1.0
        qy2 = 1.0
      ENDIF


* combine xs and qy
      DO iw = 1, nw - 1
        sig = yg(iw)
        DO i = 1, nz
          sq(j-1,i,iw) = sig * qy1
          sq(j,  i,iw) = sig * qy2
        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE mk34(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! DIBK

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for                =*
*=  DIBK photolysis:                                                         =*
*=        (CH3)3CCOC(CH3)3 ->  products                                      =*
*=                                                                           =*
*=  Cross section:  same as MIBK                                             =*
*=  Quantum yield:  based on Calvert et al. (2011 )                          =*
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

      INTEGER myld, qyvers(3)
      CHARACTER(llog) :: logmsg(nlog)
      CHARACTER(lspc) :: spc

***************** DIBK photolysis *************************


* Setting photolysis index and quantum yield options
      j = j+1
      jlabel(j) = 'DIBK -> NI products'
      j = j+1
      jlabel(j) = 'DIBK -> i-C4H9COCH3 + C3H6'

      spc = "DIBK"
      qyvers = (/1, 2, 1/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = 'estimates from Calvert et al. 2011 book'
      logmsg(2) = 'scaled externally'
      CALL set_option(2,spc,"yld",logmsg,qyvers,myld)


* cross sections

      OPEN(UNIT=kin,FILE='DATAJ1/MCMext/KET/MIBK.abs',
     $     STATUS='old')
      do i = 1,6
         read(kin,*)
      enddo

      n = 111
      DO i = 1, n
         READ(kin,*) x1(i), y1(i)
      ENDDO
      CLOSE(kin)

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)


* quantum yields
      IF(myld==1) THEN
        qy1 = 0.18
        qy2 = 0.36
      ELSEIF(myld==2) THEN
        qy1 = 1.0
        qy2 = 1.0
      ENDIF


* combine xs and qy
      DO iw = 1, nw - 1
        sig = yg(iw)
        DO i = 1, nz
          sq(j-1,i,iw) = sig * qy1
          sq(j,  i,iw) = sig * qy2
        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE mk35(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! di-sec-butyl ketone

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for                =*
*=  di-sec-butyl ketone photolysis:                                          =*
*=        di-sec-butyl ketone ->  products                                   =*
*=                                                                           =*
*=  Cross section:  same as DIPK                                             =*
*=  Quantum yield:  based on Calvert et al. (2011 )                          =*
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

      INTEGER myld, qyvers(3)
      CHARACTER(llog) :: logmsg(nlog)
      CHARACTER(lspc) :: spc

***************** di-sec-butyl ketone photolysis *************************


* Setting photolysis index and quantum yield options
      j = j+1
      jlabel(j) = 'di-sec-butyl ketone -> NI products'
      j = j+1
      jlabel(j) = 'di-sec-butyl ketone -> sec-C4H9COCH2CH3 + C2H4'

      spc = "di-sec-butyl ketone"
      qyvers = (/1, 2, 1/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = 'estimates from Calvert et al. 2011 book'
      logmsg(2) = 'scaled externally'
      CALL set_option(2,spc,"yld",logmsg,qyvers,myld)


* cross sections

      OPEN(UNIT=kin,FILE='DATAJ1/MCMext/KET/di-isopropylket.abs',
     $     STATUS='old')
      do i = 1,6
         read(kin,*)
      enddo

      n = 111
      DO i = 1, n
         READ(kin,*) x1(i), y1(i)
      ENDDO
      CLOSE(kin)

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)


* quantum yields
      IF(myld==1) THEN
        qy1 = 0.3
        qy2 = 0.29
      ELSEIF(myld==2) THEN
        qy1 = 1.0
        qy2 = 1.0
      ENDIF


* combine xs and qy
      DO iw = 1, nw - 1
        sig = yg(iw)
        DO i = 1, nz
          sq(j-1,i,iw) = sig * qy1
          sq(j,  i,iw) = sig * qy2
        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE mk36(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! di-tert-butyl ketone

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for                =*
*=  di-tert-butyl ketone photolysis:                                         =*
*=        di-tert-butyl ketone ->  products                                  =*
*=                                                                           =*
*=  Cross section:  same as DIPK                                             =*
*=  Quantum yield:  based on Calvert et al. (2011 )                          =*
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

***************** di-tert-butyl ketone photolysis *************************


* Setting photolysis index and quantum yield options
      j = j+1
      jlabel(j) = 'di-t-butyl ketone -> NI products'

      spc = "di-t-butyl ketone"
      qyvers = (/1, 2, 1/) ! options for vers = 1/2/0 (TUV/MCM&GECKO-A/individual)
      logmsg(1) = 'estimates from Calvert et al. 2011 book'
      logmsg(2) = 'scaled externally'
      CALL set_option(2,spc,"yld",logmsg,qyvers,myld)


* cross sections

      OPEN(UNIT=kin,FILE='DATAJ1/MCMext/KET/di-isopropylket.abs',
     $     STATUS='old')
      do i = 1,6
         read(kin,*)
      enddo

      n = 111
      DO i = 1, n
         READ(kin,*) x1(i), y1(i)
      ENDDO
      CLOSE(kin)

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)


* quantum yields
      IF(myld==1) THEN
        qy = 0.84
      ELSEIF(myld==2) THEN
        qy = 1.0
      ENDIF


* combine xs and qy
      DO iw = 1, nw - 1
        sig = yg(iw)
        DO i = 1, nz
          sq(j,  i,iw) = sig * qy
        ENDDO
      ENDDO

      END

*=============================================================================*
