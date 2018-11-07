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

******** Aldehyde Photochemistry
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


******** Keto Photochemistry

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


******** Ketene Photochemistry

*mk19.  ketene + hv -> CO2 + CO + H2
      CALL mk19(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*mk20.  methylketene + hv -> C2H4 + CO
      CALL mk20(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)


******** Dicarbonyl Photochemistry

*mb01.  butenedial + hv -> 3H-furan-2-one
      CALL mb01(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*mb02.  2-oxo pentenedial + hv -> products (all three channels)
      CALL mb02(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*mb03.  E,E-2,4-hexadienedial + hv -> Z-3,4-Diformyl-cyclobutene
      CALL mb03(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*mb04.  CH3COCH=CHCOCH3 + hv -> CH3CO + CH=CHCOCH3
      CALL mb04(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*ma03.  pinonaldehyde + hv -> R + CO + HO2
      CALL mb05(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*ma13.  caronaldehyde + hv -> R + CO + HO2
      CALL mb06(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)


******** Nitrate Photochemistry

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


******** PAN Photochemistry

*mn01.  PANs + hv -> products
      CALL mn11(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)


******** Dinitrate Photochemistry

*md01.  CH3CH(NO3)CH2NO3 + hv -> products (both channels)
      CALL md01(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*md02.  CH3CH2CH(NO3)CH2NO3 + hv -> products (both channels)
      CALL md02(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*md03.  CH3CH(NO3)CH(NO3)CH3 -> CH3CH(NO3)CH(O.)CH3 + NO2
      CALL md03(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*md04.  CH2(NO3)CH=CHCH2NO3 + hv -> CH2(NO3)CH=CHCH2O + NO2
      CALL md04(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*md05.  CH2=CHCH(NO3)CH2NO3 + hv -> products (both channels)
      CALL md05(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*md07.  generic unsaturated dinitrates
      CALL md07(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*md06.  1-methyl-cyclohexyl-1,2-dinitrate + hv -> products (both channels)
      CALL md06(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)


******** Hydroperoxide Photochemistry

*md01.  (CH3)3COOH + hv -> (CH3)3CO + OH
      CALL mh01(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)


******** Multifunctional chromophore Photochemistry

*mm01.  CH3CH2COCH2NO3 -> CH3CH2COCH2O + NO2
      CALL mm01(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*mm02.  CH3COCH(NO3)CH3 + hv -> CH3COCH(O.)CH3 + NO2
      CALL mm02(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*mm03.  2-oxo-cyclohexyl nitrate + hv -> RO. + NO2
      CALL mm03(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*mm04.  CH3COCH2CH2CH(OOH)CH3 + hv -> RO. + OH
      CALL mm04(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*mm05.  oxohexyl-hydroperoxide + hv -> RO. + OH
      CALL mm05(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)


******** Criegee Radical Photochemistry

*mr01.  CH2OO + hv -> HCHO + O(3P)
      CALL mr01(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*mr02.  CH3CHOO + hv -> CH3CHO + O(3P)
      CALL mr02(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*mr03.  C2H5CHOO + hv -> C2H5CHO + O(3P)
      CALL mr03(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*mr04.  (CH3)2COO + hv -> CH3COCH3 + O(3P)
      CALL mr04(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)


******** nitrocompound Photochemistry

*mn13.  CH3NO2 + hv -> CH3 + NO2 (+ HONO channel for RNO2)
      CALL mn13(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)
*mn14.  C2H5NO2 + hv -> products
      CALL mn14(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)


******** Generic j values for multifunctional chromophore photochemistry

*mk25.  lKET5 --> products
      CALL mk25(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*mm06.  C2ALDqy1
!     CALL mm06(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*mm07.  C3ALDqy1
 !    CALL mm07(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*mm08.  C4ALDqy1
!     CALL mm08(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*mm09.  C5ALDqy1
!     CALL mm09(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*mm10.  C6ALDqy1
!     CALL mm10(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*mm11.  C7ALDqy1
 !    CALL mm11(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*mm12.  C8ALDqy1
!     CALL mm12(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)



******** Sensitivity studies / Tests

*mr01.  n-Aldehyde parameterisations
C      CALL t02(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

****************************************************************

      IF (j .GT. kj) STOP 'Limit of rxns reached. Increase kj.'
      RETURN
      END
