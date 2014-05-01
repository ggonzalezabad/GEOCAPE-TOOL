!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: GC_Routine.F90
!
! !DESCRIPTION: This routine takes in an input variable, does something to it,
!  and then sends out an output variable.
!\\
!\\
! !INTERFACE:
!
SUBROUTINE GC_Routine( input, output )
!
! !USES:
!
  USE GC_SomethingMod
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
! GC_SomethingMod.F90
!
! !SYSTEM ROUTINES: 
! None
!
! !FILES USED:  
! GC_SomethingMod.F90
!
! !REVISION HISTORY: 
! 21 May 2008 - R. Yantosca - Initial version
!
! !REMARKS:
! Protex is great!
! 
!EOP
!------------------------------------------------------------------------------
!BOC

  !%%% Your code goes here! %%%

END SUBROUTINE GC_Routine
!EOC
