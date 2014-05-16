! Things to do as of March 20, 2013
! cfrac_jacobians is unnormalized
! Update to latest version of VLIDORT (2.6): (separate SC and MS diff layering)
! Ocean BRDF might not be working (inconsistent arguments), Add land BRDF capability
! Add/improve ctp, surface pressure, aerosol plume height, aerosol plume width weighting function
! Implement option whether to keep output at original grid when inserting a cloudy layer
! Merge Solar/thermal infrared depends on spectral region
! Need to get a better high resolution solar reference (better than 0.01 nm, ACE for 
!      mid/thermal infrared, Kitt peak for UV/Visible/Near infrared)  
! Need to find new reflectance spectra on globe (combine ASTER and OMI/MODIS)
! Temperature weighting functions only consider air density and o3 effects 
!     (temperature jacobians from hitran as well)
! Break the main program to several sub-programs
! Check phase functions legendre moments (definition and usage)
! Compress NETCDF data (need to use netcdf4 library)

program geocape_tools_v2p6

!  This is the GEOCAPE simulation tool for UV/Vis (solar)
!  Stand-alone software, apart from the VLIDORT dependencies

!  Based on the tool created by Rob Spurr
!  Versions 1.1-1.4 for the UV/Vis simulation of Stokes-vector
!  Versions 1.5-1.6 for the IR simulation of Stokes-vector
!  Versions 1.1-1.3
!  Radiances and Jacobians in a clear-sky Rayleigh-scattering atmosphere

!  ******************* Important note **************************************
!  All VLIDORT inputs are either HARD_WIRED or DERIVED from Geocape variables

!  Version 1.1: R. Spurr, RT SOLUTIONS Inc., 10-11 June 2009
!  Version 1.2: R. Spurr, RT SOLUTIONS Inc., 22 June 2009
!  Version 1.3: X. Liu, UMBC/Gest.,          23 June 2009
!  Version 1.4: R. Spurr, complete revision of UV, including aerosol
!  Version 1.5: R. Spurr, IR tool development
!  Version 1.6: V. Natraj, Caltech, Vis tool for GEOCAPE for GSFC and TES profiles, 30 Oct 2009
!  Version 1.7: X. Liu, Add interface to convolve cross sections with slit functions using
!          two methods (high-resolution calculation or effective cross section), March 30, 2011
!  Version 1.8: Modify user input atmospheric profiles
!               Output upwelling/downwelling radiances at any altitudes
!               Add ocean CoxMunk BRDF
!               Use wavelength or wavenumber
!               Improve cross section preparation (more flexible) and add more gases
!               Remove dependency on ozone for weighting functions
!               Add AMF and scattering weightings
!               Input surface albedo sepctra
!               Add satellite viewing geometry calculation
!               Implement Lambertian cloud/scattering clouds
!               Generate aerosol plume
!               Implement multiple types of aerosol/clouds, either user input or 
!                   generated profiles, compute aerosol/cloud weighting functions
!  Version 2.0: Modularization, upgraded to Vlidort 2.6, April 2013 gga

!  ##################################################################
!  ##################################################################
!          A R G U M E N T S    and    I N C L U D E  F I L E S
!  ##################################################################
!  ##################################################################

!  ------------------
!  Modules to be used
!  ------------------
  USE VLIDORT_INPUTS
  USE VLIDORT_AUX
  USE VBRDF_LINSUP_MASTERS_M
  USE VLIDORT_SUP_ACCESSORIES

!  USE GC_parameters_module
  USE GC_variables_module
  USE GC_error_module,      ONLY: error_exit
  USE GC_read_input_module, ONLY: read_control_file, check_input, &
                                  open_process_logfile
  USE GC_profiles_module,   ONLY: prepare_profiles
  USE GC_aerosols_module,   ONLY: read_aer_control_file, aerosol_profiles
  USE GC_clouds_module,     ONLY: read_cld_control_file, cloud_profiles
  USE GC_solar_module,      ONLY: prepare_solar_spectrum
  USE GC_xsections_module,  ONLY: read_xsec_filenames, prepare_xsec
  USE GC_surface_module,    ONLY: prepare_surface
  USE GC_Vlidort_module,    ONLY: Vlidort_GC_config, save_results, &
                                  Vlidort_set_optical,             &
                                  Vlidort_cloud_and_calculation,   &
                                  VBRDF_TO_VLIDORT
  USE GC_convolution_module,ONLY: convolve_slit
  USE GC_netcdf_module,     ONLY: netcdf_output

  IMPLICIT NONE

!  ##################################################################
!  ##################################################################
!                  S T A R T    of    C O D E
!  ##################################################################
!  ##################################################################
   ! ---------------
   ! Open error file
   ! ---------------
   CALL open_process_logfile(yn_error)
   IF (yn_error) CALL error_exit (yn_error)

   ! ------------------------------
   ! Save Vlidort control file name
   ! ------------------------------
   CALL GETARG(2, vlidort_control_file)
   vlidort_control_file = TRIM(ADJUSTL(vlidort_control_file))

   ! -----------------------------------
   ! Save Vlidort BRDF control file name
   ! -----------------------------------
   CALL GETARG(2, vlidort_vbrdf_control_file)
   vlidort_vbrdf_control_file = TRIM(ADJUSTL(vlidort_vbrdf_control_file))
   
   ! -----------------
   ! Read control file
   ! -----------------
   CALL read_control_file(yn_error)
   IF (yn_error) CALL error_exit (yn_error)

   ! ------------------------------
   ! Check for input incosistencies
   ! ------------------------------
   CALL check_input(yn_error)
   IF (yn_error) CALL error_exit (yn_error)

!  To be removed gga
!  Exception handling initialization
   nmessages = 0
   messages  = ' '
   fail      = .false.

   du = 66

!  ##################################################################
!  ##################################################################
!                  D A T A   E X T R A C T I O N
!  ##################################################################
!  ##################################################################

   ! ===============
   ! PROFILE SECTION
   ! ===============
   CALL prepare_profiles(yn_error)
   IF (yn_error) CALL error_exit (yn_error)

   ! ==================================================
   ! Cloud settings: read cloud control file and set up
   ! ==================================================
   ! --------------------
   ! Initialize variables
   ! --------------------
   ncld               = 0
   cld_flags          = .FALSE.
   do_cld_columnwf    = .FALSE.
   do_cod_Jacobians   = .FALSE.
   do_cssa_Jacobians  = .FALSE.
   do_cfrac_Jacobians = .FALSE.
   do_ctp_Jacobians   = .FALSE. 
   cld_profile = 0.0d0; tcld_profile = 0.0d0; tcldtau0 = 0.0d0

   ! -----------------------
   ! Read cloud control file
   ! -----------------------
   IF (do_clouds) THEN
      CALL read_cld_control_file(yn_error)
      IF (yn_error) CALL error_exit (yn_error)
   END IF

   ! ---------------------
   ! Sep up cloud profiles
   ! ---------------------
   CALL cloud_profiles(yn_error)
   IF (yn_error) CALL error_exit (yn_error)

   ! ----------------------
   ! Not implemented yet!!!
   ! ----------------------
   do_ctp_jacobians = .FALSE. 
   
   ! -------------------------------
   ! Re-work levels altitudes due to
   ! cloud top pressure
   ! -------------------------------
   mid_pressures(1:GC_nlayers)    = EXP((LOG(pressures(0:GC_nlayers-1)) + &
                                    LOG(pressures(1:GC_nlayers)))/2.0)
   mid_temperatures(1:GC_nlayers) = (temperatures(0:GC_nlayers-1) +       &
                                    temperatures(1:GC_nlayers))/2.0
   mid_heights(1:GC_nlayers)      = (heights(0:GC_nlayers-1) +            &
                                    heights(1:GC_nlayers))/2.0

   ! -------------------------
   ! Debug code, to be removed
   ! -------------------------

   IF ( do_debug_geocape_tool ) THEN
      WRITE(du,'(/a)')'Debug from the Profile setter--------'
      WRITE(du,'(2i5)')GC_nlayers, ngases
      WRITE(du,'(20d14.6)')(gas_totalcolumns(g),g=1,ngases)
      DO n = 1, GC_nlayers
         WRITE(du,5),n,heights(n), pressures(n),mid_temperatures(n),&
         aircolumns(n),daircolumns_dT(n), (gas_partialcolumns(n,g),g=1,ngases)
      ENDDO
5     FORMAT(i3,3f10.4,1p20e15.6)
   END IF

   ! =============================================
   ! Read aerosol control file and set up aerosols
   ! =============================================
   ! --------------------
   ! Initialize variables
   ! --------------------
   naer               = 0
   aer_flags          = .FALSE.
   do_aer_columnwf    = .FALSE.
   do_aod_Jacobians   = .FALSE.
   do_assa_Jacobians  = .FALSE.
   do_aerph_Jacobians = .FALSE.
   do_aerhw_Jacobians = .FALSE.
   aer_profile = 0.0d0; taer_profile = 0.0d0; taertau0 = 0.0d0

   IF (do_aerosols) THEN
      CALL read_aer_control_file(yn_error)
      IF (yn_error) CALL error_exit (yn_error)
    
      ! ---------------------------------
      ! Read or generate aerosol profiles
      ! ---------------------------------
      ! call get_aerprof_from_atmosprof(pmaxl, pmaxc, GC_nlayers, pnc, profids, profile_data, &
      !     maxaer, naer, aer_types, aer_profile)
      ! Already read the profile and interpolate to the cloud-modified grid
      ! ---------------------------------
      CALL aerosol_profiles(yn_error)
      IF (yn_error) CALL error_exit(yn_error)
   END IF
   taertau0 = SUM(taer_profile(1:GC_nlayers)) 

   ! -------------------
   ! Not implemented yet
   ! -------------------
   do_aerph_Jacobians  = .FALSE.
   do_aerhw_Jacobians  = .FALSE.

   IF (.NOT. do_Jacobians) THEN
      do_QU_Jacobians     = .FALSE.
      do_T_Jacobians      = .FALSE.
      do_sfcprs_Jacobians = .FALSE.
      do_aod_Jacobians    = .FALSE.
      do_assa_Jacobians   = .FALSE.
      do_aerph_Jacobians  = .FALSE.
      do_aerhw_Jacobians  = .FALSE.
      do_cod_Jacobians    = .FALSE.
      do_cssa_Jacobians   = .FALSE.
      do_ctp_Jacobians    = .FALSE.
      do_cfrac_Jacobians  = .FALSE.
      do_AMF_calculation  = .FALSE.
   END IF

   ! ================
   ! GEOMETRY SECTION
   ! ================
   ! --------------------------------------------------------
   ! Either read from input file or calculate based on inputs
   ! --------------------------------------------------------
   hour   = INT(utc)
   minute = INT((utc - hour) * 60.d0)
   second = (utc - hour - minute / 60.d0) * 3600.d0
   IF (do_sat_viewcalc) THEN

      tmzn          = 0.0
      geometry_data = satellite_angles(year, month, day, hour, minute,    &
                                       second, tmzn, latitude, longitude, &
                                       satlat, satlon, satalt)
      
      GC_sun_positions(1) = geometry_data(4)
      GC_view_angles(1)   = geometry_data(2)
      GC_azimuths(1)      = geometry_data(1)-(geometry_data(3) + 180.0)

      IF (GC_azimuths(1) .LT.   0.d0) GC_azimuths(1) = GC_azimuths(1)+360.d0
      IF (GC_azimuths(1) .GE. 360.d0) GC_azimuths(1) = GC_azimuths(1)-360.d0
      
   ENDIF

   ! ===========
   ! Solar stuff
   ! ===========
   CALL prepare_solar_spectrum(yn_error)
   IF (yn_error) CALL error_exit (yn_error)

   ! ====================
   ! CROSS SECTIONS STUFF
   ! ====================
   ! ---------------------------------------------------------------
   ! At the end of the control file, xsec filenames must be provided
   ! Files must be in database_dir/SAO_crosssections.
   ! ---------------------------------------------------------------
   CALL read_xsec_filenames(yn_error)
   IF (yn_error) CALL error_exit (yn_error)

   ! ----------------------------------------
   ! Read x-sections and convolve if necessary
   ! ----------------------------------------
   CALL prepare_xsec(yn_error)
   IF (yn_error) CALL error_exit (yn_error)

   !  Debug output
   if ( do_debug_geocape_tool ) then
      write(du,'(/a)')'Debug from the Xsecs setter--------'
      do n = 1, nlambdas
         write(du,'(i5,f10.5,1p150e12.4)') n,lambdas(n), gas_xsecs(n,1,1:ngases), &
              Rayleigh_xsecs(n), Rayleigh_depols(n)
      enddo 
   endif

   ! ===============
   ! SURFACE SECTION
   ! ===============
   CALL prepare_surface(yn_error)
   IF (yn_error) CALL error_exit (yn_error)
   !  Debug output
    
   IF ( do_debug_geocape_tool .AND. (use_lambertian .OR. do_lambertian_cld) ) THEN
      WRITE(du,'(/a)')'Debug from the Surface setter--------'
      DO n = 1, nlambdas
         WRITE(du,'(I5,2f10.4)') n,lambdas(n), ground_ler(n)
      ENDDO
   ENDIF

!  ##################################################################
!  ##################################################################
!      V L I D O R T    I N I T I A L I Z A T I O N   S E C T I O N
!  ##################################################################
!  ##################################################################

   ! -------------------------------------------------------
   ! Read Vlidort VBRDF supplement input file and initialize
   ! control arrays defined in GC_variables module
   ! -------------------------------------------------------
   CALL VBRDF_LIN_INPUTMASTER ( vlidort_vbrdf_control_file, & ! Input
                                VBRDF_Sup_In,               & ! Outputs
                                VBRDF_LinSup_In,            & ! Outputs
                                VBRDF_Sup_InputStatus )       ! Outputs

   ! -----------------------------------------------
   ! Exception handling. Use Type structure directly
   ! -----------------------------------------------
   IF ( VBRDF_Sup_InputStatus%BS_STATUS_INPUTREAD.ne.vlidort_success ) then
      open(1,file = 'VLIDORT_VBRDF_ReadInput.log', status = 'unknown')
      if ( VBRDF_Sup_InputStatus%BS_STATUS_INPUTREAD.ne.VLIDORT_WARNING ) THEN
         WRITE(1,*)' FATAL ERRORS:   Wrong input from VBRDF input file-read'
      else
         WRITE(1,*)' WARNINGS    :   Wrong input from VBRDF input file-read'
      endif
      WRITE(1,*)'  ------ Here are the messages and actions '
      write(1,'(A,I3)')'    ** Number of messages = ',VBRDF_Sup_InputStatus%BS_NINPUTMESSAGES
      DO N = 1, VBRDF_Sup_InputStatus%BS_NINPUTMESSAGES
         write(1,'(A,I3,A,A)')'Message # ',N,' : ',Adjustl(trim(VBRDF_Sup_InputStatus%BS_INPUTMESSAGES(N)))
         write(1,'(A,I3,A,A)')'Action  # ',N,' : ',Adjustl(trim(VBRDF_Sup_InputStatus%BS_INPUTACTIONS(N)))
      ENDDO
      close(1)
      if ( VBRDF_Sup_InputStatus%BS_STATUS_INPUTREAD.ne.VLIDORT_WARNING ) THEN
         STOP'Fatal Read-input fail: Look at file VLIDORT_VBRDF_ReadInput.log'
      else
         Write(*,*)'Going on....but warning from Read-input VBRDF: Look at file VLIDORT_VBRDF_ReadInput.log'
      endif
   ENDIF

   ! -----------------------------
   ! Do not want debug restoration
   ! -----------------------------
   DO_DEBUG_RESTORATION = .false.

   ! ---------------------------------
   ! A normal calculation will require
   ! ---------------------------------
   BS_NMOMENTS_INPUT = 2 * VBRDF_Sup_In%BS_NSTREAMS - 1

   ! -----------------
   ! VLIDORT BRDF call
   ! -----------------
   CALL VBRDF_LIN_MAINMASTER ( &
        DO_DEBUG_RESTORATION,    & ! Inputs
        BS_NMOMENTS_INPUT,       & ! Inputs
        VBRDF_Sup_In,            & ! Inputs
        VBRDF_LinSup_In,         & ! Inputs
        VBRDF_Sup_Out,           & ! Outputs
        VBRDF_LinSup_Out,        & ! Outputs
        VBRDF_Sup_OutputStatus )   ! Output Status

   !  Exception handling
   IF ( VBRDF_Sup_OutputStatus%BS_STATUS_OUTPUT .EQ. VLIDORT_SERIOUS ) THEN
      write(*,*)'geocape_tool_v2p6: program failed, VBRDF calculation aborted'
      write(*,*)'Here are the error messages from the VBRDF supplement : - '
      write(*,*)' - Number of error messages = ',VBRDF_Sup_OutputStatus%BS_NOUTPUTMESSAGES
      do i = 1, VBRDF_Sup_OutputStatus%BS_NOUTPUTMESSAGES
         write(*,*) '  * ',adjustl(Trim(VBRDF_Sup_OutputStatus%BS_OUTPUTMESSAGES(I)))
      enddo
   ENDIF

   ! -------------------------------------------
   ! Copy VBRDF outputs to VLIDORT's BRDF inputs
   ! -------------------------------------------
   CALL VBRDF_TO_VLIDORT(yn_error)

   ! -------------------------------------------------------
   ! Read Vlidort control file and initialize control arrays
   ! defined in GC_variables module
   ! -------------------------------------------------------
   CALL VLIDORT_INPUT_MASTER ( vlidort_control_file, & ! Input
                               VLIDORT_FixIn,        & ! Outputs
                               VLIDORT_ModIn,        & ! Outputs
                               VLIDORT_InputStatus)    ! Outputs
   
   ! ------------------
   ! Exception handling
   ! ------------------
   IF ( VLIDORT_InputStatus%TS_status_inputread .NE. vlidort_success ) THEN
      OPEN(1,file = 'V2p6_VLIDORT_ReadInput.log',status = 'unknown')
      WRITE(1,*)' FATAL:   Wrong input from VLIDORT input file-read'
      WRITE(1,*)'  ------ Here are the messages and actions '
      WRITE(1,'(A,I3)')'    ** Number of messages = ',&
           VLIDORT_InputStatus%TS_NINPUTMESSAGES
      DO N = 1, VLIDORT_InputStatus%TS_NINPUTMESSAGES
         NF = LEN_STRING(VLIDORT_InputStatus%TS_INPUTMESSAGES(N))
         NA = LEN_STRING(VLIDORT_InputStatus%TS_INPUTACTIONS(N))
         WRITE(1,'(A,I3,A,A)')'Message # ',N,' : ',&
              VLIDORT_InputStatus%TS_INPUTMESSAGES(N)(1:NF)
         WRITE(1,'(A,I3,A,A)')'Action  # ',N,' : ',&
              VLIDORT_InputStatus%TS_INPUTACTIONS(N)(1:NA)
      ENDDO
      CLOSE(1)
      STOP'Read-input fail: Look at file V2p6_VLIDORT_ReadInput.log'
   ENDIF


   ! -----------------------------------
   ! Hard wired and GC control dependent
   ! Vlidort options.
   ! -----------------------------------
   CALL Vlidort_GC_config (yn_error)

   ! ------------------------------------------------------
   ! This is the checking routine for compatibility between
   ! VLIDORT and VBRDF setup
   ! ------------------------------------------------------
   CALL VLIDORT_VBRDF_INPUT_CHECKER ( &
        VBRDF_Sup_In,             & ! Inputs
        VLIDORT_FixIn,            & ! Inputs
        VLIDORT_ModIn,            & ! Inputs
        VLIDORT_VBRDFCheck_Status ) ! Outputs

   ! ------------------
   ! Exception handling
   ! ------------------
   IF ( VLIDORT_VBRDFCheck_Status%TS_STATUS_INPUTCHECK .ne. &
        vlidort_success ) then
      WRITE(*,*)' FATAL: Main and BRDFSup inputs are incompatible'
      WRITE(*,*)'  ------ Here are the messages and actions '
      write(*,'(A,I3)')'    ** Number of messages = ',&
           VLIDORT_VBRDFCheck_Status%TS_NCHECKMESSAGES
      DO N = 1, VLIDORT_VBRDFCheck_Status%TS_NCHECKMESSAGES
         NF = LEN_STRING(VLIDORT_VBRDFCheck_Status%TS_CHECKMESSAGES(N))
         NA = LEN_STRING(VLIDORT_VBRDFCheck_Status%TS_ACTIONS(N))
         write(*,'(A,I3,A,A)')'Message # ',N,': ',&
              VLIDORT_VBRDFCheck_Status%TS_CHECKMESSAGES(N)(1:NF)
         write(*,'(A,I3,A,A)')'Action  # ',N,': ',&
              VLIDORT_VBRDFCheck_Status%TS_ACTIONS(N)(1:NA)
      ENDDO
      STOP'Checking fail: Look at file V2p6_VLIDORT_BRDFcheck.log'
   ENDIF
   

!  ##################################################################
!  ##################################################################
!  M A I N    W A V E L E N G T H - L O O P    C A L C U L A T I O N
!  ##################################################################
!  ##################################################################

   ! ---------------------
   ! =====================
   ! Start wavelength loop
   ! =====================
   ! ---------------------
   DO w = 1, nlambdas

      ! ----------------------------------------------
      ! Set optical properties for Vlidort calculation
      ! ----------------------------------------------
      CALL Vlidort_set_optical (yn_error)

      IF (yn_error) CALL error_exit(yn_error)

      ! ---------------------------------------------------------------------
      ! Cloud loop, clear and cloudy part; Vlidort calculation takes place
      ! inside the loop. Intizialize output: stokes, profile & surface jac
      ! obians at each lambda. Set cloud albedo, reset Rayleigh flag, Bulk
      ! properties at each wavelength: gas absorptions, Rayleigh scattering,
      ! optical depth, total scattering, cloudy parts with scattering clouds,
      ! Vlidort Linearized inputs, save profilewf_sum, surfacewf_sum.
      ! ---------------------------------------------------------------------
      CALL Vlidort_cloud_and_calculation (yn_error)

      ! ============
      ! save results
      ! ============
      CALL save_results (yn_error)

          ! -------------------
   END DO ! End wavelength loop
          ! -------------------

   ! ----------------
   ! Close debug file
   ! ----------------
   IF ( do_debug_geocape_tool ) THEN
      CLOSE(du)
   ENDIF

   ! ----------------------------------------------
   ! Need to convolve radiances with slit functions
   ! ----------------------------------------------
   IF (.NOT. do_effcrs .AND. lambda_resolution /= 0.0d0 ) THEN

      CALL convolve_slit(yn_error)
      
   END IF

!  ##################################################################
!  ##################################################################
!                 W R I T E    R E S U L T S
!  ##################################################################
!  ##################################################################

   ! -----------------
   ! NETCDF file write
   ! -----------------
   CALL netcdf_output (yn_error)

   ! ##################
   ! ------------------
   ! Final call to exit
   ! ------------------
   ! ##################
   CALL error_exit(yn_error)

 end program geocape_tools_v2p6
