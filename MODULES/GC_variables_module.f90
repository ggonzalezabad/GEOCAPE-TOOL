MODULE GC_variables_module

  USE VLIDORT_PARS
  USE VLIDORT_IO_DEFS
  USE VLIDORT_LININPUTS_DEF
  USE VLIDORT_LINSUP_INOUT_DEF
  USE VLIDORT_LINOUTPUTS_DEF
  USE VLIDORT_OUTPUTS_DEF

  USE VBRDF_SUP_INPUTS_DEF
  USE VBRDF_LINSUP_INPUTS_DEF
  USE VBRDF_SUP_OUTPUTS_DEF
  USE VBRDF_LINSUP_OUTPUTS_DEF

  USE GC_parameters_module

  IMPLICIT NONE
! ----------------------------
! Number of greek coefficients
! ----------------------------
  INTEGER           :: ngksec, nactgkmatc

! -------------------
! Temporal wavelength
! -------------------
  REAL(KIND=8), DIMENSION(maxflambdas) :: tmpwaves, tmpcwaves

! --------------------------------------------------------------------------------
! Variables to store partial clear/cloudy calculations directly from VLIDORT
! Dimensions, (view_geometries,stokes_components,direction(up/down_welling),cloud)
! --------------------------------------------------------------------------------
  REAL(KIND=8), DIMENSION(GC_maxuserlevels,max_geometries, maxstokes, 2, 2)                       :: stokes_clrcld
  REAL(KIND=8), DIMENSION(GC_maxuserlevels,max_szangles, maxstokes, 2, 2)                         :: stokes_flux
  REAL(KIND=8), DIMENSION(GC_maxuserlevels,max_szangles, maxstokes, 2)                            :: stokes_direct_flux
  REAL(KIND=8), DIMENSION(max_atmoswfs, maxlayers, GC_maxuserlevels,max_geometries, maxstokes, 2) :: profilewf_sum
  REAL(KIND=8), DIMENSION(max_surfacewfs, GC_maxuserlevels,max_geometries, maxstokes, 2, 2)       :: surfacewf_clrcld

  REAL(KIND=8), DIMENSION(0:max_atmoswfs, 0:maxmoms, 1:maxgksec)                 :: l_phasmoms_total_input



! =======================
! ***********************
! =======================
! Control Input variables
! =======================
! ***********************
! =======================

! ------------------------------------------
! Paths to the input/output data directories
! ------------------------------------------
   CHARACTER(LEN=max_ch_len) ::  database_dir, results_dir

! -------------------
! File names for data
! -------------------
   CHARACTER(LEN=max_ch_len) :: profile_data_filename, debug_filename, &
                                albspectra_fname, aerfile, cldfile

!  ----------------------------------------------------------
!  Aerosol/clouds flag (true if you want to include aerosols)
! -----------------------------------------------------------
   LOGICAL :: do_aerosols, do_clouds, do_lambertian_cld
   LOGICAL :: use_aerprof, use_cldprof
  
! --------------------------------------------------------
! Vector calculation flag (true if you want polarization)
! --------------------------------------------------------
   LOGICAL :: do_vector_calculation

! ----------------------------------------------------------
! Optional output of Q and U values and their linearizations
!    ( only if the vector calculation flag has been set )
! ----------------------------------------------------------
   LOGICAL :: do_StokesQU_output

! ----------------------------------------------------------
! Jacobian flags (trace gas, temperatures, surface pressure)
! ----------------------------------------------------------
   LOGICAL :: do_Jacobians, do_QU_Jacobians
   LOGICAL :: do_AMF_calculation
   LOGICAL :: do_T_Jacobians
   LOGICAL :: do_sfcprs_Jacobians

! -----------------------------------------
! Flag for using Normalized Jacobian output
! -----------------------------------------
   LOGICAL :: do_normalized_WFoutput ! T: dI/dx*x, F; dI/dx
   LOGICAL :: do_normalized_radiance ! T: I/F, F: I

! --------------------
! Stream control input
! --------------------
   INTEGER :: nstreams_choice

! -------------------------------
! Upwelling & downwelling control
! -------------------------------
   INTEGER               :: ndir ! Number of directions (1 or 2)
   INTEGER               :: idir ! Directions loop index
   INTEGER, DIMENSION(2) :: idix ! Directions
   INTEGER :: didx ! = 1, upwelling; 2. downwelling

! -------------
! GC geometries
! -------------
   INTEGER      :: GC_n_sun_positions
   REAL(KIND=8) :: GC_sun_positions(GC_maxgeometries)
   INTEGER      :: GC_n_view_angles
   REAL(KIND=8) :: GC_view_angles(GC_maxgeometries)
   INTEGER      :: GC_n_azimuths
   REAL(KIND=8) :: GC_azimuths(GC_maxgeometries)
   INTEGER      :: GC_n_user_levels
   REAL(KIND=8) :: GC_user_levels(GC_maxuserlevels)
   INTEGER      :: GC_n_user_altitudes
   REAL(KIND=8) :: GC_user_altitudes(GC_maxuserlevels)
   LOGICAL      :: GC_do_user_altitudes

! ---------------------
! Output levels control
! ---------------------
   INTEGER :: ilev ! Level index

! ---------------
! Spectral inputs
! ---------------
   REAL(KIND=8) :: lambda_start
   REAL(KIND=8) :: lambda_finish
   LOGICAL      :: use_wavelength        ! T: input in wavelength (nm) F: cm^-1
   LOGICAL      :: do_effcrs             ! Use effective cross sections
   REAL(KIND=8) :: lambda_resolution     ! Spectral resolution at FWHM
   REAL(KIND=8) :: lambda_dw, lambda_dfw ! Spectral interval for output and fine radiances
   
! ---------------------------
! Derived spectral quantities
! ---------------------------
   INTEGER                                 :: nlambdas
   REAL(KIND=8), DIMENSION(maxflambdas)    :: lambdas
   REAL(KIND=8), DIMENSION(maxflambdas)    :: wavnums
  
! -----------------------------------------------------
! Add wavelength grids for fine/coarse wavelength grids
! -----------------------------------------------------
   INTEGER                              :: nflambdas, nclambdas
   REAL(KIND=8)                         :: edge_dw
   REAL(KIND=8), DIMENSION(maxflambdas) :: flambdas, fwavnums
   REAL(KIND=8), DIMENSION(maxlambdas)  :: clambdas, cwavnums

! ------------------------------------------------------
! Use footprint information (lon,lat,yy,mm,dd,utc,
! sza,vza,aza,Ts,Ps,fc,ws,ctp) from atmospheric profiles
! ------------------------------------------------------
   LOGICAL :: use_footprint_info

! ---------------------
! Outputs for debugging 
! ---------------------
   LOGICAL :: do_debug_geocape_tool

! -------------------------
! Use solar photons logical
! -------------------------
   LOGICAL :: use_solar_photons ! use photons/cm^2/s /nm	

! --------------------------------------------------------------------
! Possible gases are labeled as follows: O3, NO2, HCHO, SO2, H2O, GLYX
! , BRO, OCLO, IO, CO, CO2, N2O, CH4, O2, NO, HNO3, OCS
! --------------------------------------------------------------------
   INTEGER                                   :: ngases
   CHARACTER(LEN=4), DIMENSION ( maxgases )  :: which_gases

! ------------------------
! Surface albedo variables
! ------------------------
   LOGICAL                             :: use_lambertian, do_brdf_surface
   LOGICAL                             :: use_fixedalbedo, use_albspectra
   REAL(KIND=8)                        :: fixed_albedo
   LOGICAL                             :: is_ocean
   REAL(KIND=8), DIMENSION(maxlambdas) :: wavlens_um  ! wavelength in microns
   REAL(KIND=8)                        :: wind_speed
   REAL(KIND=8), DIMENSION(maxlambdas) :: water_rn ! real part of refr index
   REAL(KIND=8), DIMENSION(maxlambdas) :: water_cn ! imag part of refr index
  
! -------------------------------------------------
! Geographical, time, satellite positions, geometry
! -------------------------------------------------
   LOGICAL                    :: do_sat_viewcalc
   INTEGER(KIND=4)            :: year, month, day
   REAL(KIND=8)               :: latitude, longitude, satlon, satlat, &
                                 satalt
   REAL(KIND=8)               :: utc
   REAL(KIND=8), DIMENSION(4) :: geometry_data

! ===========================
! ***************************
! ===========================
! End control Input variables
! ===========================
! ***************************
! ===========================

! =========================
! *************************
! =========================
! Aerosol related varaibles
! =========================
! *************************
! =========================

! ----------------------------------
! Logicals for aerosols calculations
! ----------------------------------
   LOGICAL :: do_aod_Jacobians, do_assa_Jacobians, do_aerph_jacobians, &
             do_aerhw_jacobians, do_aer_columnwf

! -------------
! Aerosol Stuff 
! -------------
   INTEGER                         :: naer, naer0
   REAL(KIND=8)                    :: aer_reflambda, taertau0
   REAL(KIND=8), DIMENSION(maxaer) :: aer_tau0s, aer_z_upperlimit,       &
                                      aer_z_lowerlimit, aer_z_peakheight,&
                                      aer_half_width
   REAL(KIND=8), DIMENSION(maxaer, GC_maxlayers) :: aer_profile,         &
                                                    aer_d_profile_dtau,  &
                                                    aer_d_profile_dpkh,  &
                                                    aer_d_profile_dhfw,  &
                                                    aer_profile0
   REAL(KIND=8), DIMENSION(GC_maxlayers)         :: taer_profile
   CHARACTER(LEN=2), DIMENSION(maxaer)           :: aer_types
  
   LOGICAL,      DIMENSION(GC_maxlayers)             :: aer_flags
   REAL(KIND=8), DIMENSION(maxlambdas,GC_maxlayers)  :: aer_opdeps
   REAL(KIND=8), DIMENSION(maxlambdas,GC_maxlayers)  :: aer_ssalbs
  
! --------------------------------------------------------
! Extinction coefficients relative to specified wavelength
! --------------------------------------------------------
   REAL(KIND=8), DIMENSION(maxlambdas,GC_maxlayers)           :: aer_relqext 
  
! --------------------------------------------------------
! To save space, get aerosol properties at each wavelength
! --------------------------------------------------------
   REAL(KIND=8), DIMENSION(GC_maxlayers, 0:maxmoms, maxgksec) :: aer_phfcn  

! =============================
! *****************************
! =============================
! End aerosol related variables
! =============================
! *****************************
! =============================

! ========================
! ************************
! ========================
! Clouds related varaibles
! ========================
! ************************
! ========================

! ---------------------------------------------------------------------
! Cloud Stuff (Either Lambertian or scattering clouds, multiple clouds)
! ---------------------------------------------------------------------
   LOGICAL :: do_cod_Jacobians, do_cssa_Jacobians, do_ctp_Jacobians,     &
              do_cld_columnwf, do_cfrac_Jacobians

   LOGICAL                             :: use_zcldspec
   REAL(KIND=8)                        :: lambertian_cldalb, cfrac, ipafrac
   INTEGER                             :: ncld
   CHARACTER(len=2), DIMENSION(maxcld) :: cld_types

! -------------------------------  
! Assume homogeneuos distribution
! -------------------------------
   REAL(KIND=8), DIMENSION(maxcld)                  :: cld_bots,         &
                                                       cld_tops,         &
                                                       cld_taus  
   LOGICAL,      DIMENSION(GC_maxlayers)            :: cld_flags
   INTEGER,      DIMENSION(maxcld)                  :: cld_lowers,       &
                                                       cld_uppers
   REAL(KIND=8), DIMENSION(maxlambdas,GC_maxlayers) :: cld_opdeps
   REAL(KIND=8), DIMENSION(maxlambdas,GC_maxlayers) :: cld_ssalbs
   REAL(KIND=8), DIMENSION(maxcld, GC_maxlayers)    :: cld_profile
   REAL(KIND=8), DIMENSION(GC_maxlayers)            :: tcld_profile

! --------------------------------------------------------
! Extinction coefficients relative to specified wavelength
! --------------------------------------------------------
   REAL(KIND=8)                                      :: cld_reflambda,   &
                                                        tcldtau0
   REAL(KIND=8), DIMENSION(maxlambdas,GC_maxlayers)  :: cld_relqext 

! -------------------------------------------------------
! To save space, get clouds properties at each wavelength
! -------------------------------------------------------
   REAL(KIND=8), DIMENSION(GC_maxlayers, 0:maxmoms, maxgksec) :: cld_phfcn

! ============================
! ****************************
! ============================
! End clouds related variables
! ============================
! ****************************
! ============================

! =====================
! *********************
! =====================
! Sun related varaibles
! =====================
! *********************
! =====================

! --------------
! Time variables
! --------------
   REAL(KIND=8)           :: jday, jday0, rdis
   REAL(KIND=8), EXTERNAL :: julday

! --------------
! Solar spectrum
! --------------
   REAL(KIND=8), DIMENSION(maxflambdas)  :: solar_spec_data
   REAL(KIND=8), DIMENSION(maxlambdas)   :: solar_cspec_data

! -----------------------
! Solar spectrum filename
! -----------------------
   CHARACTER(LEN=max_ch_len) :: solar_spec_filename

! =========================
! *************************
! =========================
! End Sun related variables
! =========================
! *************************
! =========================

! =======================
! ***********************
! =======================
! X-sec related varaibles
! =======================
! ***********************
! =======================
! -------------------
! Cross-section stuff 
! --------------------------------------------------------
! cross sections types: 
! 1: does not denpendent on P, T (same for all layers)
! 2: 2nd parameterized T-dependent coefficients (e.g., O3)
! 3: depends on P, and T (e.g., those read from HITRAN)
! --------------------------------------------------------
   CHARACTER(LEN=max_ch_len) :: xsec_data_filenames(maxgases)
   CHARACTER(LEN=max_ch_len) :: xsec_data_path

  
   ! Number of absorption coefficents (1, 3, and n for the above example)
   INTEGER, DIMENSION(maxgases) :: gas_nxsecs
   INTEGER, DIMENSION(maxgases) :: gas_xsecs_type

   REAL(KIND=8), DIMENSION(maxflambdas, GC_maxlayers, maxgases) :: gas_xsecs
   REAL(KIND=8), DIMENSION(maxflambdas)                         :: o3c1_xsecs
   REAL(KIND=8), DIMENSION(maxflambdas)                         :: o3c2_xsecs
   REAL(KIND=8), DIMENSION(maxflambdas)                         :: Rayleigh_xsecs
   REAL(KIND=8), DIMENSION(maxflambdas)                         :: Rayleigh_depols

! ===========================
! ***************************
! ===========================
! End X-sec related variables
! ===========================
! ***************************
! ===========================

! =========================
! *************************
! =========================
! Surface related varaibles
! =========================
! *************************
! =========================  
  REAL(KIND=8), dimension(maxlambdas) :: ground_ler
  REAL(KIND=8), dimension(MAXSTOKES_SQ, MAX_USER_STREAMS, &
     MAX_USER_RELAZMS, MAXBEAMS ) :: Total_brdf
  REAL(KIND=8) :: WSA_CALCULATED, BSA_CALCULATED
  LOGICAL :: GCM, OUTPUT_WSABSA
  INTEGER :: NSTOKESSQ

! =============================
! *****************************
! =============================
! End surface related variables
! =============================
! *****************************
! =============================

! =========================
! *************************
! =========================
! Profile related varaibles
! =========================
! *************************
! =========================  
! ------------------------------------------------------
! Profile stuff (Output from "geocape_profile_setter_1")
! ------------------------------------------------------
! Profile data
! ------------
  REAL(kind=8), DIMENSION(pmaxl,pmaxc) :: profile_data
  REAL(kind=8), DIMENSION(15):: footprint_data
  CHARACTER(len=256)         :: profids
  logical                    :: profile_is_level
  INTEGER                    :: pnc   ! # of data fields in atmospheric profiles

! ---------------------------------------------------------------------   
! Number of layers has 'GC_' prefix (distinguish with VLIDORT variable)
! ---------------------------------------------------------------------
  INTEGER                    :: GC_nlayers

! ------------------------
! Derived quantities (PTH)
! ------------------------
  REAL(kind=8), DIMENSION ( 0:GC_maxlayers ) :: heights
  REAL(kind=8), DIMENSION (   GC_maxlayers ) :: mid_temperatures, mid_pressures, mid_heights
  REAL(kind=8), DIMENSION ( 0:GC_maxlayers ) :: pressures, temperatures

! ---------------------------------
! Partial columns and T-derivatives
! ---------------------------------
  REAL(kind=8), DIMENSION ( GC_maxlayers )         :: aircolumns
  REAL(kind=8), DIMENSION ( GC_maxlayers )         :: daircolumns_dT
  REAL(kind=8), DIMENSION ( GC_maxlayers,maxgases ):: gas_partialcolumns
  REAL(kind=8), DIMENSION ( maxgases )             :: gas_totalcolumns


! =============================
! *****************************
! =============================
! End profile related variables
! =============================
! *****************************
! =============================

   interface
      function satellite_angles(yr,mo,day,hr,min,sec,tz,latp,lonp,latss,lonss,hgtss)
        INTEGER(kind=4), intent(in) :: yr
        INTEGER(kind=4), intent(in) :: mo
        INTEGER(kind=4), intent(in) :: day
        INTEGER(kind=4), intent(in) :: hr
        INTEGER(kind=4), intent(in) :: min
        REAL(kind=8),    intent(in) :: sec
        REAL(kind=8),    intent(in) :: tz
        REAL(kind=8),    intent(in) :: latp
        REAL(kind=8),    intent(in) :: lonp
        REAL(kind=8),    intent(in) :: latss
        REAL(kind=8),    intent(in) :: lonss
        REAL(kind=8),    intent(in) :: hgtss
        REAL(kind=8)                :: satellite_angles(4)
      end function satellite_angles
   endinterface

   REAL(kind=8), DIMENSION(maxlambdas,GC_maxlayers)  :: opdeps
   REAL(kind=8), DIMENSION(maxlambdas,GC_maxlayers)  :: ssalbs


! ----------------------------
! Exception handling variables
! ----------------------------
!  Tool status

   logical            :: fail
   INTEGER            :: nmessages
   CHARACTER(Len=200) :: messages (maxmessages)

! -------------------
! Geocape Tool output 
! -------------------

!  I-component of Stokes vector + Jacobians

   REAL(kind=8), DIMENSION(maxlambdas,GC_maxuserlevels,GC_maxgeometries, 2) :: GC_Radiances
   REAL(kind=8), DIMENSION(maxlambdas,GC_maxuserlevels,GC_maxgeometries, 2) :: GC_flux
   REAL(kind=8), DIMENSION(maxlambdas,GC_maxuserlevels,GC_maxgeometries, 2) :: GC_direct_flux


   REAL(kind=8), DIMENSION(maxlambdas,GC_maxuserlevels,GC_maxgeometries,maxgases, 2)              :: GC_AMFs

   REAL(kind=8), DIMENSION(maxlambdas,GC_maxlayers,GC_maxuserlevels,GC_maxgeometries, 2)          :: GC_Temperature_Jacobians
   REAL(kind=8), DIMENSION(maxlambdas,GC_maxlayers,GC_maxuserlevels,GC_maxgeometries,maxgases, 2) :: GC_Tracegas_Jacobians
   REAL(kind=8), DIMENSION(maxlambdas,GC_maxlayers,GC_maxuserlevels,GC_maxgeometries, 2)          :: GC_Scattering_Weights
   REAL(kind=8), DIMENSION(maxlambdas,GC_maxlayers,GC_maxuserlevels,GC_maxgeometries, 2)          :: GC_aod_Jacobians
   REAL(kind=8), DIMENSION(maxlambdas,GC_maxlayers,GC_maxuserlevels,GC_maxgeometries, 2)          :: GC_assa_Jacobians
   REAL(kind=8), DIMENSION(maxlambdas,GC_maxlayers,GC_maxuserlevels,GC_maxgeometries, 2)          :: GC_cod_Jacobians
   REAL(kind=8), DIMENSION(maxlambdas,GC_maxlayers,GC_maxuserlevels,GC_maxgeometries, 2)          :: GC_cssa_Jacobians

   REAL(kind=8), DIMENSION(maxlambdas,GC_maxuserlevels,GC_maxgeometries, 2)                       :: GC_cfrac_Jacobians
   REAL(kind=8), DIMENSION(maxlambdas,GC_maxuserlevels,GC_maxgeometries, 2)                       :: GC_Surfalbedo_Jacobians
   REAL(kind=8), DIMENSION(maxlambdas,GC_maxuserlevels,GC_maxgeometries, 2)                       :: GC_Windspeed_Jacobians
   REAL(kind=8), DIMENSION(maxlambdas,GC_maxuserlevels,GC_maxgeometries, 2)                       :: GC_sfcprs_Jacobians

!  Suggested additional code for Q and U. etc ......................
!    ---Only bothering with Q/U Jacobians for the trace gases.............. GRONK !

   REAL(kind=8), DIMENSION(maxlambdas,GC_maxuserlevels,GC_maxgeometries, 2)                       :: GC_Qvalues
   REAL(kind=8), DIMENSION(maxlambdas,GC_maxuserlevels,GC_maxgeometries, 2)                       :: GC_Uvalues
   REAL(kind=8), DIMENSION(maxlambdas,GC_maxuserlevels,GC_maxgeometries, 2)                       :: GC_Qflux
   REAL(kind=8), DIMENSION(maxlambdas,GC_maxuserlevels,GC_maxgeometries, 2)                       :: GC_Qdirect_flux
   REAL(kind=8), DIMENSION(maxlambdas,GC_maxuserlevels,GC_maxgeometries, 2)                       :: GC_Uflux
   REAL(kind=8), DIMENSION(maxlambdas,GC_maxuserlevels,GC_maxgeometries, 2)                       :: GC_Udirect_flux

   REAL(kind=8), DIMENSION(maxlambdas,GC_maxlayers,GC_maxuserlevels,GC_maxgeometries,maxgases, 2) :: GC_Tracegas_QJacobians
   REAL(kind=8), DIMENSION(maxlambdas,GC_maxlayers,GC_maxuserlevels,GC_maxgeometries,maxgases, 2) :: GC_Tracegas_UJacobians
   REAL(kind=8), DIMENSION(maxlambdas,GC_maxlayers,GC_maxuserlevels,GC_maxgeometries, 2)          :: GC_aod_QJacobians
   REAL(kind=8), DIMENSION(maxlambdas,GC_maxlayers,GC_maxuserlevels,GC_maxgeometries, 2)          :: GC_aod_UJacobians
   REAL(kind=8), DIMENSION(maxlambdas,GC_maxlayers,GC_maxuserlevels,GC_maxgeometries, 2)          :: GC_assa_QJacobians
   REAL(kind=8), DIMENSION(maxlambdas,GC_maxlayers,GC_maxuserlevels,GC_maxgeometries, 2)          :: GC_assa_UJacobians
   REAL(kind=8), DIMENSION(maxlambdas,GC_maxlayers,GC_maxuserlevels,GC_maxgeometries, 2)          :: GC_cod_QJacobians
   REAL(kind=8), DIMENSION(maxlambdas,GC_maxlayers,GC_maxuserlevels,GC_maxgeometries, 2)          :: GC_cod_UJacobians
   REAL(kind=8), DIMENSION(maxlambdas,GC_maxlayers,GC_maxuserlevels,GC_maxgeometries, 2)          :: GC_cssa_QJacobians
   REAL(kind=8), DIMENSION(maxlambdas,GC_maxlayers,GC_maxuserlevels,GC_maxgeometries, 2)          :: GC_cssa_UJacobians

   REAL(kind=8), DIMENSION(maxlambdas,GC_maxuserlevels,GC_maxgeometries, 2)                       :: GC_cfrac_QJacobians
   REAL(kind=8), DIMENSION(maxlambdas,GC_maxuserlevels,GC_maxgeometries, 2)                       :: GC_cfrac_UJacobians
   REAL(kind=8), DIMENSION(maxlambdas,GC_maxuserlevels,GC_maxgeometries, 2)                       :: GC_Surfalbedo_QJacobians
   REAL(kind=8), DIMENSION(maxlambdas,GC_maxuserlevels,GC_maxgeometries, 2)                       :: GC_Surfalbedo_UJacobians
   REAL(kind=8), DIMENSION(maxlambdas,GC_maxuserlevels,GC_maxgeometries, 2)                       :: GC_Windspeed_QJacobians
   REAL(kind=8), DIMENSION(maxlambdas,GC_maxuserlevels,GC_maxgeometries, 2)                       :: GC_Windspeed_UJacobians
   REAL(kind=8), DIMENSION(maxlambdas,GC_maxuserlevels,GC_maxgeometries, 2)                       :: GC_sfcprs_QJacobians
   REAL(kind=8), DIMENSION(maxlambdas,GC_maxuserlevels,GC_maxgeometries, 2)                       :: GC_sfcprs_UJacobians
   
! --------------------------------------
! Local variables for optical properties
! --------------------------------------
   INTEGER                             :: nscatter
   REAL(kind=8), DIMENSION(maxscatter) :: scaco_input
   REAL(kind=8), DIMENSION(0:maxmoms, 1:maxgksec, maxscatter)     :: phasmoms_input
   REAL(kind=8), DIMENSION(0:maxmoms, 1:maxgksec)                 :: phasmoms_total_input
   REAL(kind=8) :: depol, beta_2, pRay_12, pRay_22, pRay_52, pRay_41
   REAL(kind=8) :: total_sca, omega, total_tau, total_wf, total_Qwf, total_Uwf
   REAL(kind=8) :: total_molabs, total_molsca, total_moltau
   REAL(kind=8) :: total_aersca, total_aertau, total_cldtau, total_cldsca
   REAL(kind=8) :: total_vaerssa, total_vcldssa
   REAL(kind=8) :: total_gasabs(GC_maxlayers), temp, temp_sq
   REAL(kind=8) :: ratio, xsec, dxsec_dT, L_gas, L_air, pvar
   REAL(kind=8) :: gasabs(GC_maxlayers,maxgases), abs_o3(GC_maxlayers)
   REAL(kind=8) :: O3_partialcolumns(GC_maxlayers)
   INTEGER      :: aerscaidx, cldscaidx         ! aerosol and clouds indices in the scatters
   ! Index for weighting functions
   INTEGER      :: gaswfidx, twfidx, aodwfidx, assawfidx,  codwfidx, cssawfidx, sfcprswfidx 
   INTEGER      :: N_TOTALPROFILE_WFS_ncld, N_TOTALPROFILE_WFS_wcld


   ! --------------
   ! Time variables
   ! --------------
   INTEGER(KIND=4) :: hour, minute
   REAL(KIND=8)    :: second, tmzn
 
   ! --------------
   ! Help variables
   ! --------------
   INTEGER       :: w, v, g, n, i, j, k, q, du, nspec, ib, um, ua, il, ic, nw1
   REAL(kind=8)  :: deltp
   REAL(kind=8), DIMENSION(maxlambdas, GC_maxgeometries) :: temp_rad
   REAL(kind=8), DIMENSION(maxlambdas, GC_maxlayers)     :: temp_wf

   CHARACTER(LEN=256) :: surface_data_path

   ! --------------
   ! Error handling
   ! --------------
   LOGICAL :: yn_error

   ! -----------------------
   ! File names for the data
   ! -----------------------   
   CHARACTER(LEN=256) :: fnamep, unitname, fname, varname, ds_fname

   INTEGER            :: unitlen, flen, varlen

   ! --------------------------
   ! VBRDF supplement variables
   ! --------------------------
   LOGICAL :: DO_DEBUG_RESTORATION
   INTEGER :: BS_NMOMENTS_INPUT   

   ! =========================
   ! *************************
   ! =========================
   ! Vlidort related varaibles
   ! =========================
   ! *************************
   ! =========================  
   ! ---------------------------------
   ! Filename for Vlidort control file
   ! ---------------------------------
   CHARACTER(LEN=max_ch_len) :: vlidort_control_file
   CHARACTER(LEN=max_ch_len) :: vlidort_vbrdf_control_file
   LOGICAL                   :: OPENERRORFILEFLAG
   INTEGER                   :: NA,NF

   ! VLIDORT file inputs status structure   
   TYPE(VLIDORT_Input_Exception_Handling) :: VLIDORT_InputStatus

   ! VLIDORT input structures   
   TYPE(VLIDORT_Fixed_Inputs)             :: VLIDORT_FixIn
   TYPE(VLIDORT_Modified_Inputs)          :: VLIDORT_ModIn

   ! VLIDORT linearized input structures
   TYPE(VLIDORT_Fixed_LinInputs)          :: VLIDORT_LinFixIn
   TYPE(VLIDORT_Modified_LinInputs)       :: VLIDORT_LinModIn

   ! VLIDORT output structure
   TYPE(VLIDORT_Outputs)                  :: VLIDORT_Out

   ! VLIDORT supplements i/o structure
   TYPE(VLIDORT_Sup_InOut)                :: VLIDORT_Sup

   ! VLIDORT linearized supplements i/o structure
   TYPE(VLIDORT_LinSup_InOut)             :: VLIDORT_LinSup

   ! VLIDORT linearized output structure
   TYPE(VLIDORT_LinOutputs)               :: VLIDORT_LinOut

   ! VBRDF supplement input structure
   TYPE(VBRDF_Sup_Inputs)                 :: VBRDF_Sup_In

   ! VBRDF supplement linearized input structure
   TYPE(VBRDF_LinSup_Inputs)              :: VBRDF_LinSup_In

   !  VLIDORT VBRDF supplement output structure
   TYPE(VBRDF_Sup_Outputs)                :: VBRDF_Sup_Out
   
   !  VLIDORT VBRDF supplement linearized output structure
   TYPE(VBRDF_LinSup_Outputs)             :: VBRDF_LinSup_Out

   ! VBRDF supplement input status
   TYPE(VBRDF_Input_Exception_Handling)   :: VBRDF_Sup_InputStatus

   !  VLIDORT VBRDF supplement output status structure
   TYPE(VBRDF_Output_Exception_Handling)  :: VBRDF_Sup_OutputStatus
   
   !  VBRDF supplement / VLIDORT VBRDF-related inputs consistency check status
   TYPE(VLIDORT_Exception_Handling)       :: VLIDORT_VBRDFCheck_Status


   ! =============================
   ! *****************************
   ! =============================
   ! End vlidort related varaibles
   ! =============================
   ! *****************************
   ! =============================

END MODULE GC_variables_module
