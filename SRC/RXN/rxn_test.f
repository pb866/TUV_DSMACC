*= This file contains the following subroutines, related to testing
*= the product (cross section) x (quantum yield) for photo-reactions of
*= certain species in MCM/GECKO-A, which were not yet present in TUV5.2
*= and derive possible parameterisations:
*=
*=     t01 through t02

*=============================================================================*

      SUBROUTINE t01(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! t-OH cKet

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

      REAL    qy,sig,sigOH,sigtOH,qy1,qy2,qyc
      REAL    A,lc,w,AOH,lcOH,wOH,AtOH,lctOH,wtOH
      INTEGER iz,iw,CN

************************* CH3COCH(OH)CH3 photolysis


      j = j+1
      jlabel(j) = 'lKet3 -> products'
      j = j+1
      jlabel(j) = 'lKet3(OH) -> NI products'
      j = j+1
      jlabel(j) = 'lKet3(OH) -> NII products'
      j = j+1
      jlabel(j) = 'cKet3 -> products'
      j = j+1
      jlabel(j) = 'cKet3(OH) -> products'
      j = j+1
      jlabel(j) = 'cKet3(t-OH) -> products'
      j = j+1
      jlabel(j) = 'lKet4 -> products'
      j = j+1
      jlabel(j) = 'lKet4(OH) -> NI products'
      j = j+1
      jlabel(j) = 'lKet4(OH) -> NII products'
      j = j+1
      jlabel(j) = 'cKet4 -> products'
      j = j+1
      jlabel(j) = 'cKet4(OH) -> products'
      j = j+1
      jlabel(j) = 'cKet4(t-OH) -> products'
      j = j+1
      jlabel(j) = 'lKet5 -> products'
      j = j+1
      jlabel(j) = 'lKet5(OH) -> NI products'
      j = j+1
      jlabel(j) = 'lKet5(OH) -> NII products'
      j = j+1
      jlabel(j) = 'cKet5 -> products'
      j = j+1
      jlabel(j) = 'cKet5(OH) -> products'
      j = j+1
      jlabel(j) = 'cKet5(t-OH) -> products'
      j = j+1
      jlabel(j) = 'lKet6 -> products'
      j = j+1
      jlabel(j) = 'lKet6(OH) -> NI products'
      j = j+1
      jlabel(j) = 'lKet6(OH) -> NII products'
      j = j+1
      jlabel(j) = 'cKet6 -> products'
      j = j+1
      jlabel(j) = 'cKet6(OH) -> products'
      j = j+1
      jlabel(j) = 'cKet6(t-OH) -> products'
      j = j+1
      jlabel(j) = 'lKet7 -> products'
      j = j+1
      jlabel(j) = 'lKet7(OH) -> NI products'
      j = j+1
      jlabel(j) = 'lKet7(OH) -> NII products'
      j = j+1
      jlabel(j) = 'cKet7 -> products'
      j = j+1
      jlabel(j) = 'cKet7(OH) -> products'
      j = j+1
      jlabel(j) = 'cKet7(t-OH) -> products'
      j = j+1
      jlabel(j) = 'lKet8 -> products'
      j = j+1
      jlabel(j) = 'lKet8(OH) -> NI products'
      j = j+1
      jlabel(j) = 'lKet8(OH) -> NII products'
      j = j+1
      jlabel(j) = 'cKet8 -> products'
      j = j+1
      jlabel(j) = 'cKet8(OH) -> products'
      j = j+1
      jlabel(j) = 'cKet8(t-OH) -> products'
      j = j+1
      jlabel(j) = 'lKet9 -> products'
      j = j+1
      jlabel(j) = 'lKet9(OH) -> NI products'
      j = j+1
      jlabel(j) = 'lKet9(OH) -> NII products'
      j = j+1
      jlabel(j) = 'cKet9 -> products'
      j = j+1
      jlabel(j) = 'cKet9(OH) -> products'
      j = j+1
      jlabel(j) = 'cKet9(t-OH) -> products'
      j = j+1
      jlabel(j) = 'lKet10 -> products'
      j = j+1
      jlabel(j) = 'lKet10(OH) -> NI products'
      j = j+1
      jlabel(j) = 'lKet10(OH) -> NII products'
      j = j+1
      jlabel(j) = 'cKet10 -> products'
      j = j+1
      jlabel(j) = 'cKet10(OH) -> products'
      j = j+1
      jlabel(j) = 'cKet10(t-OH) -> products'
      j = j+1
      jlabel(j) = 'lKet11 -> products'
      j = j+1
      jlabel(j) = 'lKet11(OH) -> NI products'
      j = j+1
      jlabel(j) = 'lKet11(OH) -> NII products'
      j = j+1
      jlabel(j) = 'cKet11 -> products'
      j = j+1
      jlabel(j) = 'cKet11(OH) -> products'
      j = j+1
      jlabel(j) = 'cKet11(t-OH) -> products'
      j = j+1
      jlabel(j) = 'lKet12 -> products'
      j = j+1
      jlabel(j) = 'lKet12(OH1) -> NI products'
      j = j+1
      jlabel(j) = 'lKet12(OH1) -> NII products'
      j = j+1
      jlabel(j) = 'cKet12 -> products'
      j = j+1
      jlabel(j) = 'cKet12(OH) -> products'
      j = j+1
      jlabel(j) = 'cKet12(t-OH) -> products'
      j = j+1
      jlabel(j) = 'lKet13 -> products'
      j = j+1
      jlabel(j) = 'lKet13(OH) -> NI products'
      j = j+1
      jlabel(j) = 'lKet13(OH) -> NII products'
      j = j+1
      jlabel(j) = 'cKet13 -> products'
      j = j+1
      jlabel(j) = 'cKet13(OH) -> products'
      j = j+1
      jlabel(j) = 'cKet13(t-OH) -> products'
      j = j+1
      jlabel(j) = 'lKet14 -> products'
      j = j+1
      jlabel(j) = 'lKet14(OH) -> NI products'
      j = j+1
      jlabel(j) = 'lKet14(OH) -> NII products'
      j = j+1
      jlabel(j) = 'cKet14 -> products'
      j = j+1
      jlabel(j) = 'cKet14(OH) -> products'
      j = j+1
      jlabel(j) = 'cKet14(t-OH) -> products'
      j = j+1
      jlabel(j) = 'lKet15 -> products'
      j = j+1
      jlabel(j) = 'lKet15(OH) -> NI products'
      j = j+1
      jlabel(j) = 'lKet15(OH) -> NII products'
      j = j+1
      jlabel(j) = 'cKet15 -> products'
      j = j+1
      jlabel(j) = 'cKet15(OH) -> products'
      j = j+1
      jlabel(j) = 'cKet15(t-OH) -> products'
      j = j+1
      jlabel(j) = 'lKet16 -> products'
      j = j+1
      jlabel(j) = 'lKet16(OH) -> NI products'
      j = j+1
      jlabel(j) = 'lKet16(OH) -> NII products'
      j = j+1
      jlabel(j) = 'cKet16 -> products'
      j = j+1
      jlabel(j) = 'cKet16(OH) -> products'
      j = j+1
      jlabel(j) = 'cKet16(t-OH) -> products'
      j = j+1
      jlabel(j) = 'lKet17 -> products'
      j = j+1
      jlabel(j) = 'lKet17(OH) -> NI products'
      j = j+1
      jlabel(j) = 'lKet17(OH) -> NII products'
      j = j+1
      jlabel(j) = 'cKet17 -> products'
      j = j+1
      jlabel(j) = 'cKet17(OH) -> products'
      j = j+1
      jlabel(j) = 'cKet17(t-OH) -> products'
      j = j+1
      jlabel(j) = 'lKet18 -> products'
      j = j+1
      jlabel(j) = 'lKet18(OH) -> NI products'
      j = j+1
      jlabel(j) = 'lKet18(OH) -> NII products'
      j = j+1
      jlabel(j) = 'cKet18 -> products'
      j = j+1
      jlabel(j) = 'cKet18(OH) -> products'
      j = j+1
      jlabel(j) = 'cKet18(t-OH) -> products'
      j = j+1
      jlabel(j) = 'lKet19 -> products'
      j = j+1
      jlabel(j) = 'lKet19(OH) -> NI products'
      j = j+1
      jlabel(j) = 'lKet19(OH) -> NII products'
      j = j+1
      jlabel(j) = 'cKet19 -> products'
      j = j+1
      jlabel(j) = 'cKet19(OH) -> products'
      j = j+1
      jlabel(j) = 'cKet19(t-OH) -> products'
      j = j+1
      jlabel(j) = 'lKet20 -> products'
      j = j+1
      jlabel(j) = 'lKet20(OH) -> NI products'
      j = j+1
      jlabel(j) = 'lKet20(OH) -> NII products'
      j = j+1
      jlabel(j) = 'cKet20 -> products'
      j = j+1
      jlabel(j) = 'cKet20(OH) -> products'
      j = j+1
      jlabel(j) = 'cKet20(t-OH) -> products'


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
* (A1(OH) / 1E-20 cm2) = A + 1.75
* (A2(OH) / 1E-20 cm2) = A * 1.33
* (lc(OH) / nm)       = lc - 8.0
* w(OH)               = w
*
* if OH is on quaternary carbon atom, parameters change to:
*
* (A1(t-OH) / 1E-20 cm2) = A *1.2
* (A2(t-OH) / 1E-20 cm2) = A
* (lc(t-OH) / nm)       = lc - 9.6
* w(t-OH)               = w


* quantum yields
! factor of 0.5 is used to represent the equal branching between NI and NII reactions
! Output has to be applied to every Norrish channel separately
      qy  = 0.75*0.5
      qy1 = 0.34
      qy2 = 0.28
      qyc = 0.14


* combine
      DO iw = 1, nw - 1
         DO iz = 1, nz
           DO CN = 3,20
             A     = 4.845e-3 * tlev(iz) + 2.364*log(FLOAT(CN)) + 1.035
             A     = A*1.e-20
             lc    = 1.71e-2   * tlev(iz) + 1.718 * FLOAT(CN) + 265.344
             w     = 9.8e-3    * tlev(iz) - 0.728 * FLOAT(CN) + 29.367
             AOH  = A + 1.75e-20 ! alternatively: *1.33
             lcOH  = lc - 8.
             wOH   = w
             AtOH = A*0.86
             lctOH = lc - 8.
             wtOH  = w
             sig    = A*exp(-((wc(iw) - lc) / w)**2)
             sigOH  = AOH*exp(-((wc(iw) - lcOH) / wOH)**2)
             sigtOH = AtOH*exp(-((wc(iw) - lctOH) / wtOH)**2)
             sq(j-125+CN*6,iz,iw) = sig * qy
             sq(j-124+CN*6,iz,iw) = sigOH * qy1
             sq(j-123+CN*6,iz,iw) = sigOH * qy2
             sq(j-122+CN*6,iz,iw) = sig * qyc
             sq(j-121+CN*6,iz,iw) = sigOH * qyc
             sq(j-120+CN*6,iz,iw) = sigtOH * qyc
           ENDDO
         ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE t02(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) !C10 – C80 aldehydes

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

      INTEGER i, n, cn
      REAL x1(kdata)
      REAL y1(kdata)

* local

      REAL yg(kw)
      REAL qy
      REAL sig
      INTEGER iw

***************** n-aldehyde photolysis *************************


* Setting photolysis index and quantum yield options
      DO cn = 10, 80
        j = j + 1
        WRITE(jlabel(j),"('n-C',I2,'Ald')") cn
      ENDDO


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

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j-7),0.,0.,yg)


* quantum yields


* combine with quantum yields and map on grid:

      DO cn = 10, 80
        DO iw = 1, nw - 1
          sig = yg(iw)
          qy = 1. / (3. + 1e-3*exp(float(cn)/2*wc(iw)/979.8))
          DO i = 1, nz
            sq(j+cn-80,i,iw) = sig * qy
          ENDDO
        ENDDO
      ENDDO

      END

*=============================================================================*
