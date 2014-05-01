!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: GC_error_module.f90
!
! !DESCRIPTION: This module contains error routines
!\\
!\\
! !INTERFACE: 
!
MODULE GC_error_module
! 
!
! !USES
!
  USE GC_parameters_module, ONLY: errunit, max_ch_len
!
  IMPLICIT NONE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: error_exit
  PUBLIC :: write_err_message
  PUBLIC :: write_log_message
  PUBLIC :: open_process_logfile
!
! !PUBLIC DATA MEMBERS:
!
  ! ----------
  ! Error file
  ! ----------
  CHARACTER(LEN=12), PARAMETER :: errfname = 'GC_error.log'
!
! !REVISION HISTORY:
!  April 2013 - G. Gonzalez Abad - Initial Version
!
!EOP
!------------------------------------------------------------------------------

  CONTAINS

!BOP
!
! !IROUTINE:  error_exit
!
! !DESCRIPTION: This routine provides a message with the overall success
!               of the run.
!\\
!\\
! !INTERFACE: error_extit (yn_fail)
!
  SUBROUTINE error_exit ( yn_fail )

    IMPLICIT NONE
!
! !INPUT/OUTPUT PARAMETERS: 
!
    ! ---------------
    ! Input variables
    ! ---------------
    LOGICAL, INTENT (IN) :: yn_fail
!
! !REVISION HISTORY: 
!  Apri 2013 - G. Gonzalez Abad - Initial Version
!
!EOP
!BOC    
    IF ( yn_fail ) THEN
       CALL write_err_message ( &
            .TRUE., "ERRORs encoutered! Look in file "//errfname//" for messages" )
       STOP 1
    ELSE
       CALL write_log_message ( .TRUE., "AMF program executed normally." )
       STOP 0
    END IF
    RETURN

  END SUBROUTINE error_exit
!EOC

!BOP
!
! !IROUTINE: open_process_logfile
!
! !DESCRIPTION: This routine opens the error file
!\\
!\\
! !INTERFACE: open_process_logfile(yn_fail)
!
  SUBROUTINE open_process_logfile ( yn_fail )

    IMPLICIT NONE
!
! !OUTPUT PARAMETERS: 
!
  ! ---------------
  ! Output variable
  ! ---------------
    LOGICAL, INTENT (OUT) :: yn_fail
!
! !REVISION HISTORY: 
!  Apri 2013 - G. Gonzalez Abad - Initial Version
!
!EOP
!BOC    
    ! --------------
    ! Local variable
    ! --------------
    INTEGER :: ios

    yn_fail   = .FALSE.

    OPEN (UNIT=errunit, FILE=errfname, STATUS='UNKNOWN', ACTION='WRITE', IOSTAT=ios)
    IF ( ios /= 0 ) THEN
       yn_fail = .TRUE.
       CALL write_err_message ( &
            .TRUE., 'ERROR: Unable to open log-file '//errfname//' for writing.' )
       RETURN
    END IF

    WRITE (UNIT=errunit, FMT='(A,/,A,/,A,/)') &
         "! ############################################################ !", &
         "!        Process and Error Log File for GC tool                !", &
         "! ############################################################ !"

    RETURN
    
  END SUBROUTINE open_process_logfile
!EOC

!BOP
!
! !IROUTINE: write_err_message
!
! !DESCRIPTION: This routine outputs error message to error file and screen if
!               yn_verb set .TRUE.
!\\
!\\
! !INTERFACE: write_err_message(yn_verb, msg)
!  
  SUBROUTINE write_err_message ( yn_verb, msg )

    IMPLICIT NONE
!
! !INPUT PARAMETERS: 
!
    ! ---------------
    ! Input variables
    ! ---------------
    LOGICAL,           INTENT (IN) :: yn_verb
    CHARACTER (LEN=*), INTENT (IN) :: msg
!
! !REVISION HISTORY: 
!  Apri 2013 - G. Gonzalez Abad - Initial Version
!
!EOP
!BOC
    ! --------------
    ! Local variable
    ! --------------
    CHARACTER (LEN=max_ch_len) :: lmsg

    ! ---------------------
    ! Compose total message
    ! ---------------------
    lmsg = '<err> --- '//TRIM(ADJUSTL(msg))

    ! ---------------------------------
    ! Write message to process log file
    ! ---------------------------------
    WRITE (UNIT=errunit, FMT='(A)') TRIM(ADJUSTL(lmsg))

    ! ---------------------------------------
    ! Write message to screen if echo actived
    ! ---------------------------------------
    IF ( yn_verb) WRITE (*,'(A)')  TRIM(ADJUSTL(lmsg))

    RETURN
  END SUBROUTINE write_err_message
!EOC

!BOP
!
! !IROUTINE: write_log_message
!
! !DESCRIPTION: This routine writes a log message to the error file and
!               screen if yn_verb set .TRUE.
!\\
!\\
! !INTERFACE: write_log_message(yn_verb, msg)
!

  SUBROUTINE write_log_message ( yn_verb, msg )

    IMPLICIT NONE

    ! ---------------------------------------------------------
    ! Almost the same as WRITE_ERR_MESSAGE but with a slightly
    ! different formatting. No energy to have an All-In-One.
    ! ---------------------------------------------------------
!
! !INPUT PARAMETERS: 
!
    ! ---------------
    ! Input variables
    ! ---------------
    LOGICAL,           INTENT (IN) :: yn_verb
    CHARACTER (LEN=*), INTENT (IN) :: msg
!
! !REVISION HISTORY: 
!  Apri 2013 - G. Gonzalez Abad - Initial Version
!
!EOP
!BOC
    ! --------------
    ! Local variable
    ! --------------
    CHARACTER (LEN=max_ch_len) :: lmsg
    
    ! ---------------------
    ! Compose total message
    ! ---------------------
    lmsg = '<log> --- '//TRIM(ADJUSTL(msg))
    
    ! ---------------------------------
    ! Write message to process log file
    ! ---------------------------------
    WRITE (UNIT=errunit, FMT='(A)') TRIM(ADJUSTL(lmsg))
    
    ! ------------------------------------------
    ! If in verbose mode, echo message to screen
    ! ------------------------------------------
    IF ( yn_verb) WRITE (*,'(A)')  TRIM(ADJUSTL(lmsg))

    RETURN
  END SUBROUTINE write_log_message
!EOC

END MODULE GC_error_module
