!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: GC_SomethingIncludeFile.h
!
! !DESCRIPTION: This include file contains the various parameters that will
!   allow the module and routine to do stuff to various things in various
!   routines in various places.
!\\
!\\
! !PUBLIC TYPES:
!
  TYPE t_GeosChemSomething
     !%%% declare stuff here %%%
  END TYPE t_GeosChemSomething
!
! !PUBLIC MEMBER FUNCTIONS:
! None
!
! !PUBLIC DATA MEMBERS:
!
  INTEGER(ESMF_KIND_I8), PUBLIC, PARAMETER :: myIntParam   ! INTEGER value
  REAL(ESMF_KIND_I8),    PUBLIC, PARAMETER :: myRealParam  ! REAL*8 value
!
! !REVISION HISTORY: 
!  21 May 2008 - R. Yantosca - Initial Version
!
! !REMARKS:
!
!EOP
!------------------------------------------------------------------------------
