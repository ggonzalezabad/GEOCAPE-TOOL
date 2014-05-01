!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: saort_Module.f90
!
! !DESCRIPTION: This module contains the data type to declare a Something
!   object and the methods to work with the Something object. 
!\\
!\\
! !INTERFACE: 
!
MODULE saort_Module.f90
! 
! !USES:
!
  USE ESMF_Mod
!
  IMPLICIT NONE
!
! !PUBLIC TYPES:
!
  TYPE t_saortSomething 
     !... declare stuff here
  END TYPE t_saortSomething
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: saort_SomethingRoutine1
  PUBLIC :: saort_SomethingFunction1
!
! !PUBLIC DATA MEMBERS:
!
  INTEGER(ESMF_KIND_I4), PUBLIC :: myPublicVariable  ! public data variable
!
! !REVISION HISTORY:
!  March 2013 - G. Gonzalez Abad - Initial Version
!
! !REMARKS:
! Protex is great!
!
!EOP
!------------------------------------------------------------------------------

CONTAINS

!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  saort_SomethingRoutine1
!
! !DESCRIPTION: This routine does something to the input variable and returns
! the result in the output variable.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE saort_SomethingRoutine1( input, inpout, output, status )
!
! !INPUT PARAMETERS: 
!
    INTEGER(ESMF_KIND_I4), INTENT(IN) :: input   ! Input variable 
!
! !INPUT/OUTPUT PARAMETERS: 
!
    INTEGER(ESMF_KIND_I4), INTENT(IN) :: inpout  ! In/out variable
!
! !OUTPUT PARAMETERS:
!
    INTEGER(ESMF_KIND_I4), INTENT(IN) :: output  ! Output variable
!
! !REVISION HISTORY: 
!  March 2013 - G. Gonzalez Abad - Initial Version
!
! !REMARKS:
! Protex is great!
!
!EOP
!------------------------------------------------------------------------------
!BOC

    !%%% Your code goes here %%%

  END SUBROUTINE saort_SomethingRoutine1
!EOC

!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  saort_SomethingFunction1
!
! !DESCRIPTION: This function does something to the input variable and returns
! the result in the value variable.
!\\
!\\
! !INTERFACE:
!
  FUNCTION saort_SomethingFunction1( input ) RESULT( value )
!
! !INPUT PARAMETERS: 
!
    INTEGER(ESMF_KIND_I4), INTENT(IN) :: input   ! Input variable 
!
! !OUTPUT PARAMETERS:
!
    INTEGER(ESMF_KIND_I4), INTENT(IN) :: value   ! Output variable
!
! !REVISION HISTORY: 
!  March 2013 - G. Gonzalez Abad - Initial Version
!
! !REMARKS:
! Protex is great!
!
!EOP
!------------------------------------------------------------------------------
!BOC
    
    !%%% Your code goes here! %%%

  END FUNCTION GC_SomethingFunction1
!EOC

END MODULE saort_SomethingMod
