!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: saort_Routine.f90
!
! !DESCRIPTION: This routine takes in an input variable, does something to it,
!  and then sends out an output variable.
!\\
!\\
! !INTERFACE:
!
SUBROUTINE saort_Routine( input, output )
!
! !USES:
!
  USE saort_SomethingMod
!
! !INPUT PARAMETERS: 
!
  REAL(ESMF_KIND_R8), INTENT(IN) :: input    ! input variable
!
! !OUTPUT PARAMETERS:
!
  REAL(ESMF_KIND_R8), INTENT(IN) :: output   ! output variable
!
! !BUGS:  
! None known at this time
!
! !SEE ALSO: 
! saort_SomethingMod.f90
!
! !SYSTEM ROUTINES: 
! None
!
! !FILES USED:  
! saort_SomethingMod.F90
!
! !REVISION HISTORY: 
! March 2013 - G. Gonzalez Abad - Initial version
!
! !REMARKS:
! Protex is great!
! 
!EOP
!------------------------------------------------------------------------------
!BOC

  !%%% Your code goes here! %%%

END SUBROUTINE saort_Routine
!EOC
