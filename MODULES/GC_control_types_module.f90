MODULE GC_control_types_module

  USE GC_parameters_module, ONLY: max_ch_len

  IMPLICIT NONE

  PUBLIC :: aerosol_control

  TYPE :: aerosol_control
     LOGICAL :: do_aerosols = .false. ! Do aerosols
     LOGICAL :: use_aerprof = .false. ! Use profile from atmospheric model file
     CHARACTER(LEN=max_ch_len) :: aerfile ! Aerosol control file
     LOGICAL :: do_aod_jacobians = .false. ! Aerosol optical depth jacobians
     LOGICAL :: do_assa_jacobians = .false. ! Aerosol single scattering albedo jacobians
     LOGICAL :: do_aerph_jacobians = .false. ! Aerosol peak heigh jacobians
     LOGICAL :: do_aerhw_jacobians = .false. ! Aerosol half width jacobians
     LOGICAL :: do_aerre_jacobians = .false. ! Aerosol exponential relaxation jacobians
     LOGICAL :: do_aer_columnwf = .false. ! Aerosol column weighting function
     REAL(KIND=8) :: reflambda
     INTEGER :: naer ! Number of aerosol types
     CHARACTER(LEN=2), ALLOCATABLE, DIMENSION(:) :: types ! SU, DU, BC, OC, SC, SF
     CHARACTER(LEN=3), ALLOCATABLE, DIMENSION(:) :: profiles ! BOX EXP GDF
     REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: tau0s, & ! Aerosol total loading
     z_upperlimit, & ! Aerosol layer upper limit
     z_lowerlimit, & ! Aerosol layer lower limit
     z_peakheight, & ! Aerosol layer peak height (only GDF)
     half_width, & ! Aerosol layer half width (only GDF)
     relaxation ! Aerosol layer exponential relaxation (only exponential)
  END TYPE aerosol_control

END MODULE GC_control_types_module
