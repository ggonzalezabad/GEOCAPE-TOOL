!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: GC_SomethingMod.F90
!
! !DESCRIPTION: This module contains the data type to declare a Something
!   object and the methods to work with the Something object. 
!\\
!\\
! !INTERFACE: 
!
MODULE GC_SomethingMod
! 
! !USES:
!
  USE ESMF_Mod
!
  IMPLICIT NONE
!
! !PUBLIC TYPES:
!
  TYPE t_GeosChemSomething 
     !... declare stuff here
  END TYPE t_GeosChemSomething
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: GC_SomethingRoutine1
  PUBLIC :: GC_SomethingFunction1
!
! !PUBLIC DATA MEMBERS:
!
  INTEGER(ESMF_KIND_I4), PUBLIC :: myPublicVariable  ! public data variable
!
! !REVISION HISTORY:
!  21 May 2008 - R. Yantosca - Initial Version
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
! !IROUTINE:  GC_SomethingRoutine1
!
! !DESCRIPTION: This routine does something to the input variable and returns
! the result in the output variable.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GC_SomethingRoutine1( input, inpout, output, status )
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
!  21 May 2008 - R. Yantosca - Initial Version
!
! !REMARKS:
! Protex is great!
!
!EOP
!------------------------------------------------------------------------------
!BOC

    !%%% Your code goes here %%%

  END SUBROUTINE GC_SomethingRoutine1
!EOC

!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  GC_SomethingFunction1
!
! !DESCRIPTION: This function does something to the input variable and returns
! the result in the value variable.
!\\
!\\
! !INTERFACE:
!
  FUNCTION GC_SomethingFunction1( input ) RESULT( value )
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
!  21 May 2008 - R. Yantosca - Initial Version
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

END MODULE GC_SomethingMod
