      SUBROUTINE mcmext(nw,wl,nz,tlev,airden,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Load various "weighting functions" (products of cross section and        =*
*=  quantum yield at each altitude and each wavelength).  The altitude       =*
*=  dependence is necessary to ensure the consideration of pressure and      =*
*=  temperature dependence of the cross sections or quantum yields.          =*
*=  The actual reading, evaluation and interpolation is done in separate     =*
*=  subroutines for ease of management and manipulation.  Please refer to    =*
*=  the inline documentation of the specific subroutines for detail          =*
*=  information.                                                             =*
*=  In this subroutine an addition to the reactions in swchem is treated     =*
*=  as needed for MCM-GECKO.                                                 =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section * quantum yield (cm^2) for each          (O)=*
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
      REAL wl(kw)

      INTEGER nz
      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER(lcl) jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:
      INTEGER j

* local:
      REAL wc(kw)
      INTEGER iw
*_______________________________________________________________________

* complete wavelength grid

      DO 5, iw = 1, nw - 1
         wc(iw) = (wl(iw) + wl(iw+1))/2.
 5    CONTINUE

*_______________________________________________________________________

******** Aldehyde photochemistry
*ma01.  n-C3H7CHO + hv -> products (Norrish type I + II)
      CALL ma01(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*ma02.  i-C3H7CHO + hv -> products
      CALL ma02(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*ma03.  n-C4H9CHO + hv -> products (Norrish type I + II)
      CALL ma03(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*ma04.  i-C4H9CHO + hv -> products (Norrish type I + II)
      CALL ma04(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*ma05.  sec-C4H9CHO + hv -> products (Norrish type I + II)
      CALL ma05(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*ma06.  t-C4H9CHO + hv -> C4H9 + CHO
      CALL ma06(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*ma19.  neoC5H11CHO + hv -> Norrish type I + II products
      CALL ma19(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*ma07.  n-C5H11CHO + hv -> products (Norrish type I + II)
      CALL ma07(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*ma08.  n-C6H13CHO + hv -> products (Norrish type I + II)
      CALL ma08(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*ma09.  Glycidaldehyde + hv -> products (both channels)
      CALL ma09(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*ma10.  Crotonaldehyde + hv -> products (all 3 channels)
      CALL ma10(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*ma11.  2-hexenal + hv -> products (all 3 channels)
      CALL ma11(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*ma12.  C4H9C(C2H5)CHO + hv -> products
      CALL ma12(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*ma13.  CH3CH=C(CH3)CHO + hv -> Products (all 3 channels)
      CALL ma13(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*ma14.  CH3C(CH3)=CHCHO + hv -> Products (all 3 channels)
      CALL ma14(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*ma15.  2,4-Hexadienal + hv -> Products (all 3 channels)
      CALL ma15(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*ma16.  (CH3)2C(OH)CHO + hv -> (CH3)2COH + CHO
      CALL ma16(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*ma17.  octanal + hv -> products and larger n-aldehydes
      CALL ma17(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*ma18.  linear alpha,beta-unsaturated aldehydes + hv -> products
      CALL ma18(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)


******** Keto photochemistry

*mk01.  diethyl ketone + hv -> products (both channels)
      CALL mk01(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*mk02.  methyl n-propyl ketone + hv -> products (all 4 channels)
      CALL mk02(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*mk03.  C4H9COCH3 + hv -> CH3CH=CH2 + CH2=C(OH)CH3
      CALL mk03(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*mk04.  C3H7COC2H5 + hv -> products
      CALL mk04(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*mk27.  4-heptanone + hv -> products
      CALL mk27(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*mk28.  4-octanone + hv -> products
      CALL mk28(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*mk29.  n-propyl isopropyl ketone + hv -> products
      CALL mk29(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*mk30.  methyl neopentyl ketone + hv -> products
      CALL mk30(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*mk31.  n-propyl isobutyl ketone + hv -> products
      CALL mk31(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*mk32.  3-Me-4-heptanone + hv -> products
      CALL mk32(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*mk33.  2,2-Me-3-hexanone + hv -> products
      CALL mk33(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*mk34.  DIBK + hv -> products
      CALL mk34(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*mk35.  di-sec-butyl ketone + hv -> products
      CALL mk35(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*mk36.  di-tert-butyl ketone + hv -> products
      CALL mk36(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*mk05.  MIPK + hv -> products
      CALL mk05(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*mk06.  MIBK + hv -> products
      CALL mk06(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*mk07.  4-Me-2-hexanone + hv -> Norrish type II products
      CALL mk07(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*mk08.  5-Me-2-hexanone + hv -> products
      CALL mk08(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*mk09.  di-isopropyl ketone + hv -> Norrish type I products
      CALL mk09(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*mk13.  cyclopropanone + hv -> products
      CALL mk13(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*mk10.  cyclobutanone + hv -> products
      CALL mk10(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*mk12.  cyclopentanone + hv -> products
      CALL mk12(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*mk11.  cyclohexanone + hv -> products
      CALL mk11(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*mk26.  cycloheptanone + hv -> products
      CALL mk26(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*mk14.  ethyl vinyl ketone + hv -> Norrish type I products
      CALL mk14(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*mk15.  CH3COC2H4OH  + hv -> CH3CO + CH2CH2OH
      CALL mk15(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*mk16.  CH3COCH(OH)CH3  + hv -> CH3CO + CH3CHOH
      CALL mk16(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*mk17.  methylacetoin
      CALL mk17(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*mk18.  4-hydroxy-4-methyl-2-pentanone
      CALL mk18(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*mk21.  RCOR' + hv -> products
*       RCOR(OH) + hv -> products
      CALL mk21(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)


******** Ketene photochemistry

*mk19.  ketene + hv -> CO2 + CO + H2
      CALL mk19(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*mk20.  methylketene + hv -> C2H4 + CO
      CALL mk20(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)


******** Dicarbonyl photochemistry

*dc01.  butenedial + hv -> 3H-furan-2-one
      CALL dc01(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*dc02.  2-oxo pentenedial + hv -> products (all three channels)
      CALL dc02(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*dc03.  E,E-2,4-hexadienedial + hv -> Z-3,4-Diformyl-cyclobutene
      CALL dc03(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*dc04.  CH3COCH=CHCOCH3 + hv -> CH3CO + CH=CHCOCH3
      CALL dc04(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*dc05.  pinonaldehyde + hv -> R + CO + HO2
      CALL dc05(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*dc06.  caronaldehyde + hv -> R + CO + HO2
      CALL dc06(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)


******** Nitrate photochemistry

*mn01.  1-C5H11ONO2 + hv -> 1-C5H11O + NO2
      CALL mn01(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*mn02.  2-C5H11ONO2 + hv -> 2-C5H11O + NO2
      CALL mn02(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*mn03.  3-C5H11ONO2 + hv -> 3-C5H11O + NO2
      CALL mn03(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*mn12.  internal linear mono-nitrates
      CALL mn12(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*mn04.  c-C5H11ONO2 + hv -> c-C5H11O + NO2
      CALL mn04(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*mn05.  (CH3)2CHCH2ONO2 + hv -> i-C4H9O + NO2
      CALL mn05(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*mn06.  (CH3)2CHCH2CH2ONO2 + hv -> i-C5H11O + NO2
      CALL mn06(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*mn07.  C1 OH-subst. alkyl nitrates
      CALL mn07(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*mn08.  C≥3 OH-subst. alkyl nitrates
      CALL mn08(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*mn09.  OH-subst. iNIT
      CALL mn09(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*mn10.  OH-subst. tNIT
      CALL mn10(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)


******** PAN photochemistry

*mn01.  PANs + hv -> products
      CALL mn11(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)


******** Dinitrate photochemistry

*dn01.  CH3CH(NO3)CH2NO3 + hv -> products (both channels)
      CALL dn01(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*dn02.  CH3CH2CH(NO3)CH2NO3 + hv -> products (both channels)
      CALL dn02(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*dn03.  CH3CH(NO3)CH(NO3)CH3 -> CH3CH(NO3)CH(O.)CH3 + NO2
      CALL dn03(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*dn04.  CH2(NO3)CH=CHCH2NO3 + hv -> CH2(NO3)CH=CHCH2O + NO2
      CALL dn04(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*dn05.  CH2=CHCH(NO3)CH2NO3 + hv -> products (both channels)
      CALL dn05(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*dn07.  generic unsaturated dinitrates
      CALL dn07(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*dn06.  1-methyl-cyclohexyl-1,2-dinitrate + hv -> products (both channels)
      CALL dn06(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)


******** Hydroperoxide photochemistry

*mh01.  (CH3)3COOH + hv -> (CH3)3CO + OH
      CALL mh01(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)


******** Multifunctional chromophore photochemistry

*mf01.  CH3CH2COCH2NO3 -> CH3CH2COCH2O + NO2
      CALL mf01(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*mf02.  CH3COCH(NO3)CH3 + hv -> CH3COCH(O.)CH3 + NO2
      CALL mf02(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*mf03.  2-oxo-cyclohexyl nitrate + hv -> RO. + NO2
      CALL mf03(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*mf04.  CH3COCH2CH2CH(OOH)CH3 + hv -> RO. + OH
      CALL mf04(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*mf05.  oxohexyl-hydroperoxide + hv -> RO. + OH
      CALL mf05(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)


******** Criegee radical photochemistry

*ci01.  CH2OO + hv -> HCHO + O(3P)
      CALL ci01(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*ci02.  CH3CHOO + hv -> CH3CHO + O(3P)
      CALL ci02(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*ci03.  C2H5CHOO + hv -> C2H5CHO + O(3P)
      CALL ci03(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*ci04.  (CH3)2COO + hv -> CH3COCH3 + O(3P)
      CALL ci04(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)


******** Nitrocompound photochemistry

*mn13.  CH3NO2 + hv -> CH3 + NO2 (+ HONO channel for RNO2)
      CALL mn13(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)
*mn14.  C2H5NO2 + hv -> products
      CALL mn14(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)


******** Generic j values for multifunctional chromophore photochemistry

*mk25.  lKET5 --> products
      CALL mk25(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)


******** Halogen photochemistry

*h01.  IO + hv -> I + O(3P)
      CALL h01(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)
*h02.  HOI + hv -> I + OH
      CALL h02(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)
*h03.  OIO + hv -> I + O2
      CALL h03(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*mf06.  C2ALDqy1
!     CALL mf06(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*mf07.  C3ALDqy1
 !    CALL mf07(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*mf08.  C4ALDqy1
!     CALL mf08(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*mf09.  C5ALDqy1
!     CALL mf09(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*mf10.  C6ALDqy1
!     CALL mf10(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*mf11.  C7ALDqy1
 !    CALL mf11(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*mf12.  C8ALDqy1
!     CALL mf12(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)



******** Sensitivity studies / Tests

*t01.  n-Aldehyde parameterisations
C      CALL t02(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

****************************************************************

      IF (j .GT. kj) STOP 'Limit of rxns reached. Increase kj.'
      RETURN
      END
