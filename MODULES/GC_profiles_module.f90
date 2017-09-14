MODULE GC_profiles_module

  USE GC_parameters_module, ONLY: pmaxl, pmaxc, maxaer
  USE GC_variables_module,  ONLY: GC_nlayers, profile_data_filename, &
                                  profile_data, footprint_data,      &
                                  profids, profile_is_level, pnc,    &
                                  messages, nmessages,               &
                                  use_footprint_info, year, month,   &
                                  day, utc, longitude, latitude,     &
                                  GC_n_sun_positions,                &
                                  GC_sun_positions, GC_n_view_angles,&
                                  GC_view_angles, GC_n_azimuths,     &
                                  GC_azimuths, wind_speed, is_ocean, &
                                  use_lambertian, cfrac, cld_tops,   &
                                  ngases, which_gases, heights,      &
                                  temperatures, pressures,           &
                                  aircolumns, daircolumns_dT,        &
                                  gas_partialcolumns,                &
                                  gas_totalcolumns, &
                                  naer0, &
                                  aer_profile0, aer_ctr
  USE GC_error_module

  IMPLICIT NONE
  
CONTAINS
  
  SUBROUTINE prepare_profiles(error)
    
    IMPLICIT NONE
    
    LOGICAL, INTENT(INOUT) :: error
    
    ! ----------------
    ! Code starts here
    ! ----------------
    ! ------------------------
    ! Reading the profile file
    ! ------------------------
    CALL geocape_profile_reader_1                 &
         ( profile_data_filename, pmaxl, pmaxc,   & ! input 
         profile_data, footprint_data, profids, & ! Output
         profile_is_level,GC_nlayers, pnc,      & ! Output
         messages(nmessages+1), error )           ! Exception handling
    
    ! ----------------------------------------------
    ! Override footprint input information if needed
    ! ----------------------------------------------
    IF (use_footprint_info) THEN
       
       year      = INT(footprint_data(1))
       month     = INT(footprint_data(2))
       day       = INT(footprint_data(3))
       utc       =     footprint_data(4)
       longitude = footprint_data(5)
       
       IF (longitude .GT. 180.d0) longitude = longitude-360.d0  
       latitude  = footprint_data(6)
       GC_n_sun_positions  = 1
       GC_sun_positions(1) = footprint_data(7)
       GC_n_view_angles    = 1
       GC_view_angles(1)   = footprint_data(8)
       GC_n_azimuths       = 1
       GC_azimuths(1)      = footprint_data(9)
       wind_speed          = footprint_data(13)
       
       IF (INT(footprint_data(12)) == 1) THEN
          is_ocean = .TRUE.
       ELSE
          is_ocean = .FALSE.
       ENDIF
       
       ! --------------------------------------------------------------------
       ! ************ Use lambertian surface albedo over land ***************
       ! Only use BRDF for ocean surface when use_lambertian is .false.
       ! --------------------------------------------------------------------
!!$       IF (.NOT. use_lambertian .AND. .NOT. is_ocean) use_lambertian = .TRUE. 
       
       cfrac        = footprint_data(14)
       cld_tops(1)  = footprint_data(15)
    ENDIF
    
    ! ------------------
    ! Exception handling
    ! ------------------
    IF (error) CALL write_err_message (.FALSE., messages(nmessages+1))
    
    ! ----------------
    ! Profile settings
    ! ----------------
    !  Given this data set:
    !    - the LEVEL height grid [km]
    !    - the LAYER temperature grid [K]
    !    - the LEVEL pressure grid [hPa]
    !    - the LAYER partial columns of Air         [mol/cm^2]
    !    - the LAYER partial columns of trace gases [mol/cm^2]
    !    - the TOTAL columns of the trace gases     [mol/cm^2]
    ! --------------------------------------------------------
    
    CALL geocape_profile_setter_1                             &
         ( pmaxl, pmaxc, GC_nlayers, pnc, ngases, profile_data,  &  ! Input
         profids, profile_is_level, which_gases(1:ngases),       &  ! Input
         heights(0:GC_nlayers), temperatures(0:GC_nlayers),      &  ! Output
         pressures(0:GC_nlayers),                                &  ! Output
         aircolumns(1:GC_nlayers), daircolumns_dT(1:GC_nlayers), &  ! Output
         gas_partialcolumns(1:GC_nlayers,1:ngases),              &  ! Output
         gas_totalcolumns(1:ngases),                             &  ! Output
         error, messages(nmessages+1) )                             ! Exception handling
    
    ! ------------------
    ! Exception handling
    ! ------------------
    IF (error) CALL write_err_message (.FALSE., messages(nmessages+1))
    
    ! ------------------------------------------------------------------
    ! Read aerosol profile before modifying the altitude grids by clouds
    ! ------------------------------------------------------------------
    IF (aer_ctr%do_aerosols .AND. aer_ctr%use_aerprof) THEN
       CALL get_aerprof_from_atmosprof(pmaxl, pmaxc, GC_nlayers, pnc, &
                                       profids, profile_data, maxaer, &
                                       naer0, aer_ctr%aer_types, aer_profile0)
    ELSE
       naer0 = 0; aer_profile0 = 0.0d0
    ENDIF
    
  END SUBROUTINE prepare_profiles
  
END MODULE GC_profiles_module
