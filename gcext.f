      SUBROUTINE gcext(nw,wl,nz,tlev,airden,j,sq,jlabel)

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
*=  In this subroutine is an addition to the reactions in swchem             =*
*=  and hosts additional photo-reactions of the GEOS-CHEM mechanisms.        =*
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
      INCLUDE 'params'

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
      INTEGER  j

* local:
      REAL wc(kw)
      INTEGER iw
*_______________________________________________________________________

* complete wavelength grid

      DO 5, iw = 1, nw - 1
         wc(iw) = (wl(iw) + wl(iw+1))/2.
5     CONTINUE

*_______________________________________________________________________

******** Inorganic GEOS-CHEM Photochemistry
*gc01.  NO + hv -> N + O(3P)
      CALL gc01(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*gc01.  OCS + hv -> CO + SO2
      CALL gc02(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*gc03.  H2SO4 + hv -> 2 OH + SO2
      CALL gc03(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*gc04.  MACR + hv -> products
      CALL gc04(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*gc05.  CH2Br2 + hv -> CH2Br + Br
      CALL gc05(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)



****************************************************************

      IF (j .GT. kj) STOP 'Limit of rxns reached. Increase kj.'
      RETURN
      END
