*= This file contains the following subroutines, related to reading/loading
*= the product (cross section) x (quantum yield) for photo-reactions of
*= halogen compounds in MCM/GECKO-A, which were not yet present in TUV5.2:
*=
*=     h01 through h03

*=============================================================================*

      SUBROUTINE h01(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! IO

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for                =*
*=  iodine oxide photolysis:                                                 =*
*=        IO + hv -> I + O(3P)                                               =*
*=                                                                           =*
*=  Cross section:  IUPAC                                                    =*
*=  Quantum yield:  IUPAC                                                    =*
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
!      REAL qy
!      REAL sig
      INTEGER iw

      j = j+1
      jlabel(j) = 'IO -> I + O(3P)'

* cross sections
      OPEN(UNIT=kin,FILE='DATAJ1/MCMext/HALO/IO_IUPAC2007.abs',
     $     STATUS='old')
      do i = 1, 2
        read(kin,*)
      enddo

      n  = 25
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
      ENDDO
      CLOSE(kin)
      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)

* quantum yields
!      qy = 1.

* combine:
      DO iw = 1, nw - 1
        DO i = 1, nz
          sq(j  ,i,iw) = yg(iw)! xs * qy
        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE h02(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! HOI

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for                =*
*=  hypoiodous acid photolysis:                                              =*
*=        HOI + hv -> I + OH                                                 =*
*=                                                                           =*
*=  Cross section:  IUPAC                                                    =*
*=  Quantum yield:  IUPAC                                                    =*
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
!      REAL qy
!      REAL sig
      INTEGER iw

      j = j+1
      jlabel(j) = 'HOI -> I + OH'

* cross sections
      OPEN(UNIT=kin,FILE='DATAJ1/MCMext/HALO/HOI_iup07.abs',
     $     STATUS='old')
      do i = 1, 11
        read(kin,*)
      enddo

      n  = 43
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
        y1(i) = y1(i)*1E-20
      ENDDO
      CLOSE(kin)

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)

* quantum yields
!      qy = 1.

* combine:
      DO iw = 1, nw - 1
        DO i = 1, nz
          sq(j  ,i,iw) = yg(iw)! xs * qy
        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE h03(nw,wl,wc,nz,tlev,airden,j,sq,jlabel) ! OIO

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for                =*
*=  iodine dioxide photolysis:                                               =*
*=        OIO + hv -> I + O2                                                 =*
*=                                                                           =*
*=  Cross section:  JPL2006                                                  =*
*=  Quantum yield:  JPL2006                                                  =*
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
!      REAL qy
!      REAL sig
      INTEGER iw

      j = j+1
      jlabel(j) = 'OIO -> I + O2'

* cross sections
      OPEN(UNIT=kin,FILE='DATAJ1/MCMext/HALO/OIO_jpl06.abs',
     $     STATUS='old')
      do i = 1, 11
        read(kin,*)
      enddo

      n  = 57
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
        y1(i) = y1(i)*1E-20
      ENDDO
      CLOSE(kin)

      CALL interpol(x1,y1,kdata,n,nw,wl,jlabel(j),0.,0.,yg)

* quantum yields
!      qy = 1.

* combine:
      DO iw = 1, nw - 1
        DO i = 1, nz
          sq(j  ,i,iw) = yg(iw)! xs * qy
        ENDDO
      ENDDO

      END
