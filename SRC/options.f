      SUBROUTINE set_option(nopt, spc, settype, logmsg, opt, mopt)
*======================================================================*
*=                                                                    =*
*=  PURPOSE:                                                          =*
*=  Set cross section, quantum yield or pressure dependence options   =*
*=  to choose from different input data and make log entry.           =*
*=                                                                    =*
*=--------------------------------------------------------------------=*
*=  PARAMETERS:                                                       =*
*=  - nopt:       Number of different log entries (data options)      =*
*=  - spc:        Species name or formula for log entries             =*
*=  - settype:    Set to "abs" or "yld" for cross sections/yields     =*
*=  - logmsg:     Log message (preceeded by "<spc> <option type>")    =*
*=  - opt:        Array of cross section or yield options for         =*
*=                  i) TUV                                            =*
*=                 ii) MCM/GECKO-A                                    =*
*=                iii) individually chosen values                     =*
*=  - mopt:       Return value with chosen cross section/yield option =*
*=                                                                    =*
*======================================================================*

      IMPLICIT NONE
      INCLUDE '../params'

!     I/O variables
      CHARACTER(llog), INTENT(IN) ::  logmsg(nlog)
      CHARACTER(lspc), INTENT(IN) ::  spc
      CHARACTER(3), INTENT(IN) ::     settype
      INTEGER, INTENT(IN) ::          opt(3), nopt
      INTEGER, INTENT(OUT) ::         mopt
!     Internal variables
      INTEGER ::                      i
      CHARACTER(35) ::                typestr

!     Bounds check for options
      IF(ANY(opt>nopt)) THEN
        WRITE(*,"(3A)") "Options out of bounds for species ", trim(spc),
     &                  " in set_option!"
        WRITE(*,"(2A)") "Adjust the array in the rxn routine calling",
     &                  " set_option."
        STOP
      ENDIF

!     Check for standard or individual settings in params
!     and set options for cross sections or quantum yields respectively
      IF(vers==1)THEN
        mopt = opt(1)
       ELSEIF(vers==2)THEN
        mopt = opt(2)
       ELSEIF(vers==0) THEN
        mopt = opt(3)
       ELSE
        STOP "'vers' not set. Choose value between 0 and 2 in 'params'."
      ENDIF

!     Define type of options
      IF(settype == "abs") THEN
        typestr = " cross sections"
      ELSEIF(settype == "yld") THEN
        typestr = " quantum yields"
      ELSEIF(settype == "que") THEN
        typestr = " quenching"
      ELSEIF(settype == "prd") THEN
        typestr = " product distribution"
      ELSE
        WRITE(*,"(2A)") "'settype' in SR set_option has to be set to ",
     &                  "either 'abs', 'yld', 'que' or 'prd'."
        WRITE(*,"(2A)") "Currently set to: ", settype
        STOP
      ENDIF

!     Print log entry for individual choices
      IF(vers==1 .or. vers==2) THEN
        CONTINUE
      ELSE
        DO i = 1, nopt
          IF(mopt==i) THEN
            WRITE(kout,101) trim(spc), trim(typestr), trim(logmsg(i))
            EXIT
          ELSEIF(i==nopt) THEN
            WRITE(*,102) trim(settype), trim(spc)
            STOP "TUV"
          ENDIF
        ENDDO
      ENDIF

! Print formats
  101 FORMAT(X,2A,X,A,".")
  102 FORMAT("'m",A,"' not defined for ",A," photolysis.")

      END SUBROUTINE set_option
