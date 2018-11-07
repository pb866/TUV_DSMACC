*= This file contains the following subroutines, related to reading/loading
*= the product (cross section) x (quantum yield) for photo-reactions of
*= organic hydroperoxides in MCM-GECKO, which were not yet present in TUV5.2:
*=
*=     md01

*=============================================================================*

      SUBROUTINE mh01(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! tert-butyl hydroperoxide

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for                =*
*=  tert-butyl hydroperoxide photolysis:                                     =*
*=        (CH3)3COOH + hv -> (CH3)3CO + OH                                   =*
*=                                                                           =*
*=  Cross section:  Baasandorj et al. (2010)                                 =*
*=  Quantum yield:  Baasandorj et al. (2010)                                 =*
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

      INTEGER i, n
      REAL x1(kdata)
      REAL y1(kdata)

* local

      REAL yg(kw)
!      REAL qy
      REAL sig
      INTEGER iw


***************** tert-butyl hydroperoxide photolysis *************************

      j = j+1
      jlabel(j) = '(CH3)3COOH + hv -> (CH3)3CO + OH'


      OPEN(UNIT=kin,FILE='DATAJ1/MCMext/ROOH/tButHP.abs',
     $     STATUS='old')
      do i = 1, 6
        read(kin,*)
      enddo

      n  = 46
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
      ENDDO
      CLOSE(kin)

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)


* quantum yields
!      qy = 1.0


* combine:

      DO iw = 1, nw - 1
        sig = yg(iw)
        DO i = 1, nz
          sq(j,i,iw) = sig! * qy
        ENDDO
      ENDDO

      END

* ============================================================================*
