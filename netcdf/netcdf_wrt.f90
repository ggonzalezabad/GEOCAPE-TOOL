subroutine netcdf_wrt ( fname,        & !Filename
     lat, lon, year, month, day, utc, & !Geolocation
     nsza, nvza, naza, ngeom,         & !Number of SZA, VZA, AZA and number of geometries
     sza, vza, aza,                   & !Viewing geometry
     ulvl,                            & !User levels
     nw, ws,                          & !Number of lambdas and lambdas
     nz, zs,                          & !Number of layers and layers
     ps, ts,                          & !Pressure and temperature profile
     ngas, gases, gascol, aircol,     & !Number of gases, gases, gas partial column and air columns
     cfrac, zcldtop, cldalb,          & !Cloud fraction, cloud top and lambertian cloud albedo
     aods0, cods0,                    & !Aerosol and cloud profile
     ods, ssas,                       & !Optical depth and single scattering albedo
     aods, assas,                     & !Aerosol optical depth and single scattering albedo
     cods, cssas,                     & !Cloud optical depth and single scattering albedo
     surfalb, windspeed,              & !Surface albedo and wind speed
     nstreams, aercld_nmoms,          & !Number of streams and number of aerosol and cloud moments
     lam_dw, lam_dfw, use_wavelength, & !Fine lambda, delta fine lambda and use wavelenght or wavenumber
     aer_reflam, cld_reflam,          & !Aerosol and cloud reference wavelength
     do_vector_calc, do_QU_output,    & !LOGICAL SWITCHES
     do_Jacobian, do_QU_Jacobian,     &
     do_AMF_calc, do_T_Jacobian,      &
     do_sfcprs_Jacobian,              &
     do_aod_Jacobian,                 &
     do_assa_Jacobian,                &
     do_cod_Jacobian,                 &
     do_cssa_Jacobian,                &
     do_cfrac_Jacobian,               &
     do_aer_columnwf,                 &
     do_cld_columnwf,                 &
     do_norm_WFoutput,                &
     do_norm_radiance,                &
     use_lambertian,                  &
     do_lambertian_cld, do_upwelling, &
     do_effcrs, use_solar_photons,    &
     lam_reso, irradiance, radiance,  & !Lambda resolution, irradiance and radiance
     q, u, Gas_Jacobian, ScatWeights, & !Qvalues, Uvalues, GasJacobians, and Scattering Weights
     AMF, Gas_QJacobian,              & !Air mass factors, GasQJacobians
     Gas_UJacobian, TJacobian,        & !GasUJacobians and Temperature Jacobians
     SurfalbJacobian,                 & !Surface albedo Jacobians
     SurfalbQJacobian,                & !Surface albedo QJacobians
     SurfalbUJacobian,                & !Surface albedo UJacobians
     WSJacobian,                      & !Wind speed Jacobians
     WSQJacobian, WSUJacobian,        & !Wind speed QJacobians and wind speed UJacobians
     aodJacobian,                     & !Aerosol optical depth Jacobians
     aodQJacobian, aodUJacobian,      & !Aerosol optical depth QJacobians and UJacobians
     assaJacobian,                    & !Aerosol single scattering Jacobians
     assaQJacobian, assaUJacobian,    & !Aerosol single scattering QJacobians and UJacobians
     codJacobian,                     & !Cloud optical depth Jacobians
     codQJacobian, codUJacobian,      & !Cloud optical depth QandUJacobians
     cssaJacobian,                    & !Cloud single scattering Jacobians
     cssaQJacobian, cssaUJacobian,    & !Cloud single scattering QandUJacobians
     cfracJacobian,                   & !Cloud fraction Jacobians
     cfracQJacobian, cfracUJacobian,  & !Cloud fraction QandUJacobians
     sfcprsJacobian, sfcprsQJacobian, sfcprsUJacobian, & !Surface pressure or cloud pressure Jacobians?
     flux, qflux, uflux,              & !Total flux
     direct_flux, qdirect_flux,       & !Direct flux
     udirect_flux,                    & !Direct flux
     brdf, nsq,                       & !Results from VBRDF supplement
     do_brdf,                         & !Results from VBRDF supplement
     message, fail)   
  
  implicit none
  include '../netcdf/netcdf.inc'

  !=================================  
  !     Input/out variables
  !=================================
  character(len=256), intent(in)                  :: fname                   !Filename
  integer(kind=4),    intent(in)                  :: nw, nz, ngas            !Number of lamdas, number of layers and number of gases
  integer(kind=4),    intent(in)                  :: nsza, nvza, naza, ngeom !Number of solar zenith angles, viewing zenith angles, azimuth angles and geometries
  integer(kind=4),    intent(in)                  :: nsq                     !Number of stokes square
  integer(kind=4),    intent(in)                  :: nstreams, aercld_nmoms  !Number of streams and aerosols/clouds moments
  integer(kind=4),    intent(in)                  :: year, month, day
  character(len=4),   intent(in), dimension(ngas) :: gases                   !Name of gas, possible problem with string length
  real(kind=8), intent(in)                        :: lat, lon, utc           !Latitude, longitude, UTC
  real(kind=8), intent(in), dimension(nsza)       :: sza                     !Solar zenith angle
  real(kind=8), intent(in), dimension(nvza)       :: vza                     !Veiwing zenith angle
  real(kind=8), intent(in), dimension(naza)       :: aza                     !Azimuth zenith angle
  real(kind=8), intent(in)                        :: ulvl                    !User levels still limit to one
  real(kind=8), dimension(0:nz), intent(in)       :: zs, ps, ts              !Heights, pressure and temperature profile
  real(kind=8), dimension(nz), intent(in)         :: aods0, cods0, aircol    !Aerosol, clouds and air column profile
  real(kind=8), dimension(nz, ngas), intent(in)   :: gascol                  !Gas partial columns
  real(kind=8), intent(in)                        :: cfrac, zcldtop, cldalb  !Cloud fraction, cloud top and cloud albedo
  real(kind=8), intent(in)                        :: lam_reso, lam_dw, &     !Lambda resolution, lambda fine 
                                                     lam_dfw                 !lambda delta
  real(kind=8), intent(in)                        :: aer_reflam, cld_reflam  !Aerosol and cloud reference wavelength
  real(kind=8), dimension(nw), intent(in)         :: ws, surfalb             !Wavelengths and surface albedo
  real(kind=8), dimension(nsq,nvza,naza,nsza), &
                               intent(in)         :: brdf                    !BRDF
  real(kind=8), dimension(nw,nz), intent(in)      :: aods, assas, cods, &    !Aerosol optical depth, single scattering albedo, cloud optical depth
                                                     cssas, ods, ssas        !cloud single scattering albedo, optical depth and single scattering albedo
  real(kind=8), intent(in)                        :: windspeed               !Wind speed

  logical, intent(in)  :: do_vector_calc, do_QU_output,        & !Self explained
                          do_Jacobian, do_QU_Jacobian,         &
                          do_AMF_calc, do_T_Jacobian,          &
                          do_sfcprs_Jacobian, do_aod_Jacobian, &
                          do_assa_Jacobian, do_cod_Jacobian,   &
                          do_cssa_Jacobian, do_cfrac_Jacobian, &
                          do_aer_columnwf, do_cld_columnwf,    &
                          do_norm_WFoutput, do_norm_radiance,  &
                          use_lambertian, do_lambertian_cld,   &
                          do_upwelling, do_effcrs,             &
                          use_wavelength, use_solar_photons,   &
                          do_brdf

  real(kind=8), dimension(nw),       intent(in) :: irradiance     !Solar spectrum

  real(kind=8), dimension(nw,ngeom), intent(in) :: radiance, q, u !Radiance, q and u values

  real(kind=8), dimension(nw,ngeom), intent(in) :: flux, qflux, uflux !Flux, q and u values
  real(kind=8), dimension(nw,ngeom), intent(in) :: direct_flux, qdirect_flux, udirect_flux !Direct flux, q and u values

  real(kind=8), dimension(nw,ngeom), intent(in) :: SurfalbJacobian,  & !Surface albedo Jacobian
                                                    SurfalbQJacobian, & !Surface albedo QJacobian
                                                    SurfalbUJacobian, & !Surface albedo UJacobian
                                                    WSJacobian,       & !Wind speed Jacobian
                                                    WSQJacobian,      & !Wind speed QJacobian
                                                    WSUJacobian,      & !Wind speed UJacobian
                                                    cfracJacobian,    & !Cloud fraction Jacobian
                                                    cfracQJacobian,   & !Cloud fraction QJacobian
                                                    cfracUJacobian,   & !Cloud fraction UJacobian
                                                    sfcprsJacobian,   & !Surface pressure Jacobian
                                                    sfcprsQJacobian,  & !Surface pressure QJacobian
                                                    sfcprsUJacobian     !Surface pressure UJacobian

  real(kind=8), dimension(nw, nz, ngeom), intent(in) :: ScatWeights,   & !Scattering weigths
                                                        TJacobian,     & !Temperature Jacobian
                                                        aodJacobian,   & !Aerosol optical depth Jacobian
                                                        aodQJacobian,  & !Aerosol optical depth QJacobian
                                                        aodUJacobian,  & !Aerosol optical depth UJacobian
                                                        assaJacobian,  & !Aerosol single scattering albedo Jacobians
                                                        assaQJacobian, & !Aerosol single scattering albedo QJacobians
                                                        assaUJacobian, & !Aerosol single scattering albedo UJacobians
                                                        codJacobian,   & !Cloud optical depth Jacobians
                                                        codQJacobian,  & !Cloud optical depth QJacobians
                                                        codUJacobian,  & !Cloud optical depth UJacobians
                                                        cssaJacobian,  & !Cloud single scattering albedo Jacobians
                                                        cssaQJacobian, & !Cloud single scattering albedo QJacobians
                                                        cssaUJacobian    !Cloud single scattering albedo UJacobians

  real(kind=8), dimension(nw, ngeom, ngas), intent(in) :: AMF !Air Mass Factors

  real(kind=8), dimension(nw, nz, ngeom, ngas), intent(in) :: Gas_Jacobian,  & ! Trace gas Jacobians
                                                              Gas_QJacobian, & ! Trace gas QJacobians
                                                              Gas_UJacobian    ! Trace gas UJacobians

 character*(*), intent(inout) :: message
 logical, intent(out)         :: fail
  
  !=================================  
  !     Local variables
  !=================================
  integer :: ncid, rcode, nlen, i
  integer :: szadim, vzadim, azadim, gasdim, lvldim, laydim, wavdim, geodim, nsqdim
  integer :: szaid, vzaid, azaid, gasid, lvlid, wavid, psid, tsid, airid, &
       aer0id, cld0id, radid, qid, uid, irradid, sfcid, sfcwfid, wswfid, cfracwfid, &
       sfcprswfid, aodwfid, assawfid, codwfid, cssawfid, sfcqwfid, wsqwfid, cfracqwfid, &
       sfcprsqwfid, aodqwfid, assaqwfid, codqwfid, cssaqwfid, sfcuwfid, wsuwfid, cfracuwfid,&
       sfcprsuwfid, aoduwfid, assauwfid, coduwfid, cssauwfid, gascolid, aodsid, assasid, &
       codsid, cssasid, scatwtid, odsid, ssasid, amfid, gaswfid, gasqwfid, &
       gasuwfid, tempwfid, fluxid, dfluxid, qfluxid, ufluxid, qdfluxid, udfluxid, brdfid
  integer, dimension(2)    :: gascol_dims, wavalt_dims, wavgas_dims, wavgeo_dims
  integer, dimension(3)    :: wavaltgeo_dims, wavgeogas_dims
  integer, dimension(4)    :: wavaltgeogas_dims, brdfdim
  integer, dimension(ngas) :: gas_indices
  integer, dimension(1)    :: ndimstart1, ndimcount1
  integer, dimension(2)    :: ndimstart2, ndimcount2
  integer, dimension(3)    :: ndimstart3, ndimcount3
  integer, dimension(4)    :: ndimstart4, ndimcount4
  character(len=256)       :: gasstr
  character(len=30)        :: irrad_unitc, rad_unitc

  fail = .false.; message = ' '

  ncid = nccre(trim(fname), ncclob, rcode)
  if (rcode .eq. -1) then
     message =  ' error in netcdf_out: nccre failed'
     fail = .true.; return
  endif
  
  ! Create the dimensions of the dataset:
  szadim  = ncddef (ncid, 'nsza',   nsza,  rcode)
  vzadim  = ncddef (ncid, 'nvza',   nvza,  rcode)
  azadim  = ncddef (ncid, 'naza',   naza,  rcode)
  gasdim  = ncddef (ncid, 'ngas',   ngas,  rcode)
  lvldim  = ncddef (ncid, 'nlevel', nz+1,  rcode)
  laydim  = ncddef (ncid, 'nlayer', nz,    rcode)
  wavdim  = ncddef (ncid, 'nw',     nw,    rcode)
  geodim  = ncddef (ncid, 'ngeom',  ngeom, rcode)
  nsqdim  = ncddef (ncid, 'nstokessq', nsq, rcode)
  
  !=============================================================================
  ! Create the coordinate (aka independent) variables:
  !============================================================================= 
  szaid  = ncvdef (ncid, 'Solarzenithangle',     ncfloat, 1, szadim, rcode)
  vzaid  = ncvdef (ncid, 'Viewingzenithangle',   ncfloat, 1, vzadim, rcode)
  azaid  = ncvdef (ncid, 'Relativeazimuthangle', ncfloat, 1, azadim, rcode)
  gasid  = ncvdef (ncid, 'gas',                  nclong,  1, gasdim, rcode)
  lvlid  = ncvdef (ncid, 'zs',                   ncfloat, 1, lvldim, rcode)
  wavid  = ncvdef (ncid, 'Wavelength',           ncfloat, 1, wavdim, rcode)

  !=============================================================================
  ! Create a vector containing the often referred to spacetime coordinates:
  !=============================================================================
  gascol_dims(1) = laydim; gascol_dims(2) = gasdim
  wavalt_dims(1) = wavdim; wavalt_dims(2) = laydim
  wavgas_dims(1) = wavdim; wavgas_dims(2) = gasdim
  wavgeo_dims(1) = wavdim; wavgeo_dims(2) = geodim

  wavaltgeo_dims(1) = wavdim; wavaltgeo_dims(2) = laydim; wavaltgeo_dims(3) = geodim
  wavgeogas_dims(1) = wavdim; wavgeogas_dims(2) = geodim; wavgeogas_dims(3) = gasdim

  wavaltgeogas_dims(1) = wavdim; wavaltgeogas_dims(2) = laydim
  wavaltgeogas_dims(3) = geodim; wavaltgeogas_dims(4) = gasdim

  brdfdim(1) = nsqdim; brdfdim(2) = vzadim; brdfdim(3)=azadim; brdfdim(4) = szadim

  !=============================================================================
  ! Create the dependent variables:
  !=============================================================================
  ! ps, ts
  psid = ncvdef(ncid, 'ps', ncfloat, 1, lvldim, rcode)
  tsid = ncvdef(ncid, 'ts', ncfloat, 1, lvldim, rcode)

  ! aircol, aods0, cods0
  airid  = ncvdef(ncid, 'aircol', ncfloat, 1, laydim, rcode)
  aer0id = ncvdef(ncid, 'aods0',  ncfloat, 1, laydim, rcode)
  cld0id = ncvdef(ncid, 'cods0',  ncfloat, 1, laydim, rcode)

  ! varaibles with 1D, wavdim
  irradid = ncvdef(ncid, 'irradiance', ncfloat, 1, wavdim, rcode)
  sfcid   = ncvdef(ncid, 'surfalb',    ncfloat, 1, wavdim, rcode)

  ! variables with 2D, laydim, gasdim
  gascolid = ncvdef(ncid,  'gascol', ncdouble, 2, gascol_dims, rcode)
  
  ! variables with 2D, wavdim, laydim
  odsid   = ncvdef(ncid, 'ods',   ncfloat, 2, wavalt_dims, rcode)
  ssasid  = ncvdef(ncid, 'ssas',  ncfloat, 2, wavalt_dims, rcode)
  aodsid  = ncvdef(ncid, 'aods',  ncfloat, 2, wavalt_dims, rcode)
  assasid = ncvdef(ncid, 'assas', ncfloat, 2, wavalt_dims, rcode)
  codsid  = ncvdef(ncid, 'cods',  ncfloat, 2, wavalt_dims, rcode)
  cssasid = ncvdef(ncid, 'cssas', ncfloat, 2, wavalt_dims, rcode)

  ! variables with 2D, wavdim, geodim
  radid   = ncvdef(ncid, 'radiance',    ncfloat, 2, wavgeo_dims, rcode)
  fluxid  = ncvdef(ncid, 'flux',        ncfloat, 2, wavgeo_dims, rcode)
  dfluxid = ncvdef(ncid, 'direct_flux', ncfloat, 2, wavgeo_dims, rcode)
  if (do_vector_calc .and. do_QU_output) then
     qid      = ncvdef(ncid, 'q',            ncfloat, 2, wavgeo_dims, rcode)
     uid      = ncvdef(ncid, 'u',            ncfloat, 2, wavgeo_dims, rcode)
     qfluxid  = ncvdef(ncid, 'qflux',        ncfloat, 2, wavgeo_dims, rcode)
     ufluxid  = ncvdef(ncid, 'uflux',        ncfloat, 2, wavgeo_dims, rcode)
     qdfluxid = ncvdef(ncid, 'qdirect_flux', ncfloat, 2, wavgeo_dims, rcode)
     udfluxid = ncvdef(ncid, 'udirect_flux', ncfloat, 2, wavgeo_dims, rcode)
  endif

  if (do_Jacobian) then

     if (use_lambertian) then
        sfcwfid = ncvdef(ncid, 'surfalb_jac',   ncfloat, 2, wavgeo_dims, rcode)
     else 
        wswfid  = ncvdef(ncid, 'windspeed_jac', ncfloat, 2, wavgeo_dims, rcode)
     endif

     if (do_cfrac_Jacobian) &
          cfracwfid  = ncvdef(ncid, 'cfrac_jac',     ncfloat, 2, wavgeo_dims, rcode)
     if (do_sfcprs_Jacobian) &
          sfcprswfid = ncvdef(ncid, 'sfcprs_jac',    ncfloat, 2, wavgeo_dims, rcode)
     if (do_aer_columnwf .and. do_aod_Jacobian) &  !No column jacobians yet
          aodwfid    = ncvdef(ncid, 'aodcolwf_jac',  ncfloat, 2, wavgeo_dims, rcode)
     if (do_aer_columnwf .and. do_assa_Jacobian) & !No column jacobians yet
          assawfid   = ncvdef(ncid, 'assacolwf_jac', ncfloat, 2, wavgeo_dims, rcode)
     if (do_cld_columnwf .and. do_cod_Jacobian) &  !No column jacobians yet
          codwfid    = ncvdef(ncid, 'codcolwf_jac',  ncfloat, 2, wavgeo_dims, rcode)
     if (do_cld_columnwf .and. do_cssa_Jacobian) & !No column jacobians yet
          cssawfid   = ncvdef(ncid, 'cssacolwf_jac', ncfloat, 2, wavgeo_dims, rcode)

  endif

  if (do_QU_Jacobian) then

     if (use_lambertian) then
        sfcqwfid = ncvdef(ncid, 'surfalb_qjac',  ncfloat, 2, wavgeo_dims, rcode)
        sfcuwfid = ncvdef(ncid, 'surfalb_ujac',  ncfloat, 2, wavgeo_dims, rcode)
     else 
        wsqwfid = ncvdef(ncid, 'windspeed_qjac', ncfloat, 2, wavgeo_dims, rcode)
        wsuwfid = ncvdef(ncid, 'windspeed_ujac', ncfloat, 2, wavgeo_dims, rcode)
     endif

     if (do_cfrac_Jacobian) then
          cfracqwfid = ncvdef(ncid, 'cfrac_qjac',  ncfloat, 2, wavgeo_dims, rcode)
          cfracuwfid = ncvdef(ncid, 'cfrac_ujac',  ncfloat, 2, wavgeo_dims, rcode)
     endif
     if (do_sfcprs_Jacobian) then
          sfcprsqwfid = ncvdef(ncid, 'sfcprs_qjac',  ncfloat, 2, wavgeo_dims, rcode)
          sfcprsuwfid = ncvdef(ncid, 'sfcprs_ujac',  ncfloat, 2, wavgeo_dims, rcode)
     endif
     if (do_aer_columnwf .and. do_aod_Jacobian) then  !No column jacobians yet
          aodqwfid = ncvdef(ncid, 'aodcol_qjac',  ncfloat, 2, wavgeo_dims, rcode)
          aoduwfid = ncvdef(ncid, 'aodcol_ujac',  ncfloat, 2, wavgeo_dims, rcode)
     endif
     if (do_aer_columnwf .and. do_assa_Jacobian) then !No column jacobians yet
          assaqwfid = ncvdef(ncid, 'assacol_qjac', ncfloat, 2, wavgeo_dims, rcode)
          assauwfid = ncvdef(ncid, 'assacol_ujac', ncfloat, 2, wavgeo_dims, rcode)
     endif
     if (do_cld_columnwf .and. do_cod_Jacobian) then  !No column jacobians yet
          codqwfid = ncvdef(ncid, 'codcol_qjac',  ncfloat, 2, wavgeo_dims, rcode)
          coduwfid = ncvdef(ncid, 'codcol_ujac',  ncfloat, 2, wavgeo_dims, rcode)
     endif
     if (do_cld_columnwf .and. do_cssa_Jacobian) then !No column jacobians yet
          cssaqwfid = ncvdef(ncid, 'cssacol_qjac', ncfloat, 2, wavgeo_dims, rcode)
          cssauwfid = ncvdef(ncid, 'cssacol_ujac', ncfloat, 2, wavgeo_dims, rcode)
     endif

  endif

  ! variables with 3D, wavdim, laydim, geodim
  if (do_AMF_calc) then
     scatwtid = ncvdef(ncid, 'scatweights', ncfloat, 3, wavaltgeo_dims, rcode)
  endif

  if (do_Jacobian) then
     if (do_T_Jacobian) &
          tempwfid = ncvdef(ncid, 't_jac',    ncfloat, 3, wavaltgeo_dims, rcode)
     if (.not. do_aer_columnwf .and. do_aod_Jacobian) &
          aodwfid  = ncvdef(ncid, 'aod_jac',  ncfloat, 3, wavaltgeo_dims, rcode)  
     if (.not. do_aer_columnwf .and. do_assa_Jacobian) &
          assawfid = ncvdef(ncid, 'assa_jac', ncfloat, 3, wavaltgeo_dims, rcode)   
     if (.not. do_cld_columnwf .and. do_cod_Jacobian) &
          codwfid  = ncvdef(ncid, 'cod_jac',  ncfloat, 3, wavaltgeo_dims, rcode)  
     if (.not. do_cld_columnwf .and. do_cssa_Jacobian) &
          cssawfid = ncvdef(ncid, 'cssa_jac', ncfloat, 3, wavaltgeo_dims, rcode)  
  endif

  if (do_QU_Jacobian) then
     if (.not. do_aer_columnwf .and. do_aod_Jacobian) then 
          aodqwfid = ncvdef(ncid, 'aod_qjac', ncfloat, 3, wavaltgeo_dims, rcode)  
          aoduwfid = ncvdef(ncid, 'aod_ujac', ncfloat, 3, wavaltgeo_dims, rcode)  
     endif
     if (.not. do_aer_columnwf .and. do_assa_Jacobian) then
          assaqwfid = ncvdef(ncid, 'assa_qjac', ncfloat, 3, wavaltgeo_dims, rcode)   
          assauwfid = ncvdef(ncid, 'assa_ujac', ncfloat, 3, wavaltgeo_dims, rcode)   
     endif
     if (.not. do_cld_columnwf .and. do_cod_Jacobian) then
          codqwfid = ncvdef(ncid, 'cod_qjac', ncfloat, 3, wavaltgeo_dims, rcode)  
          coduwfid = ncvdef(ncid, 'cod_ujac', ncfloat, 3, wavaltgeo_dims, rcode)  
     endif
     if (.not. do_cld_columnwf .and. do_cssa_Jacobian) then
          cssaqwfid = ncvdef(ncid, 'cssa_qjac', ncfloat, 3, wavaltgeo_dims, rcode)  
          cssauwfid = ncvdef(ncid, 'cssa_ujac', ncfloat, 3, wavaltgeo_dims, rcode) 
     endif
  endif

  ! variables with 3D, wavdim, geodim, gasdim
  if (do_AMF_calc) then
     amfid = ncvdef(ncid, 'amf', ncfloat, 3, wavgeogas_dims, rcode)
  endif

  ! variables with 4D, wavdim, laydim, geodim, gasdim
  if (do_Jacobian) then
     gaswfid = ncvdef(ncid,  'gas_jac', ncfloat, 4, wavaltgeogas_dims, rcode)
     if (do_QU_Jacobian) then
        gasqwfid = ncvdef(ncid,  'gas_qjac', ncfloat, 4, wavaltgeogas_dims, rcode)
        gasuwfid = ncvdef(ncid,  'gas_ujac', ncfloat, 4, wavaltgeogas_dims, rcode)
     endif
  endif

  ! variables with 4d, nsqdim, nvza, naza, nsza
  if (do_brdf) then
     brdfid = ncvdef(ncid, 'BRDF', ncfloat, 4, brdfdim, rcode)
  endif
  
  !=============================================================================
  ! Assign attributes (meta-data) to the various variables:
  !============================================================================= 
  call ncapt (ncid, ncglobal, 'lon',          ncfloat, 1, real(lon, kind=4),          rcode)
  call ncapt (ncid, ncglobal, 'lat',          ncfloat, 1, real(lat, kind=4),          rcode)
  call ncapt (ncid, ncglobal, 'month',        nclong,  1, month,                      rcode)
  call ncapt (ncid, ncglobal, 'year',         nclong,  1, year,                       rcode)
  call ncapt (ncid, ncglobal, 'day',          nclong,  1, day,                        rcode)
  call ncapt (ncid, ncglobal, 'utc',          ncfloat, 1, real(utc, kind=4),          rcode)
  call ncapt (ncid, ncglobal, 'cfrac',        ncfloat, 1, real(cfrac, kind=4),        rcode)
  call ncapt (ncid, ncglobal, 'zcldtop',      ncfloat, 1, real(zcldtop, kind=4),      rcode)
  call ncapt (ncid, ncglobal, 'cldalb',       ncfloat, 1, real(cldalb, kind=4),       rcode)
  call ncapt (ncid, ncglobal, 'userlvl',      ncfloat, 1, real(ulvl, kind=4),         rcode)
  call ncapt (ncid, ncglobal, 'windspeed',    ncfloat, 1, real(windspeed, kind=4),    rcode)
  call ncapt (ncid, ncglobal, 'lam_reso',     ncfloat, 1, real(lam_reso, kind=4),     rcode)
  call ncapt (ncid, ncglobal, 'aer_reflam',   ncfloat, 1, real(aer_reflam, kind=4),   rcode)
  call ncapt (ncid, ncglobal, 'cld_reflam',   ncfloat, 1, real(cld_reflam, kind=4),   rcode)
  call ncapt (ncid, ncglobal, 'lam_dw',       ncfloat, 1, real(lam_dw, kind=4),       rcode)
  call ncapt (ncid, ncglobal, 'lam_dfw',      ncfloat, 1, real(lam_dfw, kind=4),      rcode)
  call ncapt (ncid, ncglobal, 'nstreams',     nclong,  1, nstreams,                   rcode)
  call ncapt (ncid, ncglobal, 'aercld_nmoms', nclong,  1, aercld_nmoms,               rcode)
  call ncapt (ncid, ncglobal, 'nsza',         nclong,  1, nsza,                       rcode)
  call ncapt (ncid, ncglobal, 'nvza',         nclong,  1, nvza,                       rcode)
  call ncapt (ncid, ncglobal, 'naza',         nclong,  1, naza,                       rcode)
  call ncapt (ncid, ncglobal, 'ngeometries',  nclong,  1, ngeom,                      rcode)
  call ncapt (ncid, ncglobal, 'nwavelengths', nclong,  1, nw,                         rcode)
  call ncapt (ncid, ncglobal, 'nlayers',      nclong,  1, nz,                         rcode)
  call ncapt (ncid, ncglobal, 'ngas',         nclong,  1, ngas,                       rcode)

  if (use_lambertian) then
     call ncapt (ncid, ncglobal, 'user_lambertian', ncbyte, 1, 1, rcode)
  else
     call ncapt (ncid, ncglobal, 'user_lambertian', ncbyte, 1, 0, rcode)
  endif 
  if (do_lambertian_cld) then
     call ncapt (ncid, ncglobal, 'do_lambertian_cld', ncbyte, 1, 1, rcode)
  else
     call ncapt (ncid, ncglobal, 'do_lambertian_cld', ncbyte, 1, 0, rcode)
  endif  
  if (do_effcrs) then
     call ncapt (ncid, ncglobal, 'do_effcrs', ncbyte, 1, 1, rcode)
  else
     call ncapt (ncid, ncglobal, 'do_effcrs', ncbyte, 1, 0, rcode)
  endif
  if (use_wavelength) then
     call ncapt (ncid, ncglobal, 'use_wavelength', ncbyte, 1, 1, rcode)
  else
     call ncapt (ncid, ncglobal, 'use_wavelength', ncbyte, 1, 0, rcode)
  endif
  if (do_upwelling) then
     call ncapt (ncid, ncglobal, 'do_upwelling', ncbyte, 1, 1, rcode)
  else
     call ncapt (ncid, ncglobal, 'do_upwelling', ncbyte, 1, 0, rcode)
  endif   
  if (do_norm_WFoutput) then
     call ncapt (ncid, ncglobal, 'do_norm_WFoutput', ncbyte, 1, 1, rcode)
  else
     call ncapt (ncid, ncglobal, 'do_norm_WFoutput', ncbyte, 1, 0, rcode)
  endif
  if (do_norm_radiance) then
     call ncapt (ncid, ncglobal, 'do_norm_radiance', ncbyte, 1, 1, rcode)
  else
     call ncapt (ncid, ncglobal, 'do_norm_radiance', ncbyte, 1, 0, rcode)
  endif   
  if (use_solar_photons) then
     call ncapt (ncid, ncglobal, 'use_solar_photons', ncbyte, 1, 1, rcode)
  else
     call ncapt (ncid, ncglobal, 'use_solar_photons', ncbyte, 1, 0, rcode)
  endif
  if (do_vector_calc) then
     call ncapt (ncid, ncglobal, 'do_vector_calc', ncbyte, 1, 1, rcode)
  else
     call ncapt (ncid, ncglobal, 'do_vector_calc', ncbyte, 1, 0, rcode)
  endif
  if (do_QU_output) then
     call ncapt (ncid, ncglobal, 'do_QU_output', ncbyte, 1, 1, rcode)
  else
     call ncapt (ncid, ncglobal, 'do_QU_output', ncbyte, 1, 0, rcode)
  endif
  if (do_Jacobian) then
     call ncapt (ncid, ncglobal, 'do_Jacobian', ncbyte, 1, 1, rcode)
  else
     call ncapt (ncid, ncglobal, 'do_Jacobian', ncbyte, 1, 0, rcode)
  endif
  if (do_QU_Jacobian) then
     call ncapt (ncid, ncglobal, 'do_QU_Jacobian', ncbyte, 1, 1, rcode)
  else
     call ncapt (ncid, ncglobal, 'do_QU_Jacobian', ncbyte, 1, 0, rcode)
  endif   
  if (do_AMF_calc) then
     call ncapt (ncid, ncglobal, 'do_AMF_calc', ncbyte, 1, 1, rcode)
  else
     call ncapt (ncid, ncglobal, 'do_AMF_calc', ncbyte, 1, 0, rcode)
  endif   
  if (do_T_Jacobian) then
     call ncapt (ncid, ncglobal, 'do_T_Jacobian', ncbyte, 1, 1, rcode)
  else
     call ncapt (ncid, ncglobal, 'do_T_Jacobian', ncbyte, 1, 0, rcode)
  endif  
  if (do_sfcprs_Jacobian) then
     call ncapt (ncid, ncglobal, 'do_sfcprs_Jacobian', ncbyte, 1, 1, rcode)
  else
     call ncapt (ncid, ncglobal, 'do_sfcprs_Jacobian', ncbyte, 1, 0, rcode)
  endif  
  if (do_aod_Jacobian) then
     call ncapt (ncid, ncglobal, 'do_aod_Jacobian', ncbyte, 1, 1, rcode)
  else
     call ncapt (ncid, ncglobal, 'do_aod_Jacobian', ncbyte, 1, 0, rcode)
  endif   
  if (do_assa_Jacobian) then
     call ncapt (ncid, ncglobal, 'do_assa_Jacobian', ncbyte, 1, 1, rcode)
  else
     call ncapt (ncid, ncglobal, 'do_assa_Jacobian', ncbyte, 1, 0, rcode)
  endif 
  if (do_cod_Jacobian) then
     call ncapt (ncid, ncglobal, 'do_cod_Jacobian', ncbyte, 1, 1, rcode)
  else
     call ncapt (ncid, ncglobal, 'do_cod_Jacobian', ncbyte, 1, 0, rcode)
  endif 
  if (do_cssa_Jacobian) then
     call ncapt (ncid, ncglobal, 'do_cssa_Jacobian', ncbyte, 1, 1, rcode)
  else
     call ncapt (ncid, ncglobal, 'do_cssa_Jacobian', ncbyte, 1, 0, rcode)
  endif        
  if (do_cfrac_Jacobian) then
     call ncapt (ncid, ncglobal, 'do_cfrac_Jacobian', ncbyte, 1, 1, rcode)
  else
     call ncapt (ncid, ncglobal, 'do_cfrac_Jacobian', ncbyte, 1, 0, rcode)
  endif    
  if (do_aer_columnwf) then
     call ncapt (ncid, ncglobal, 'do_aer_columnwf', ncbyte, 1, 1, rcode)
  else
     call ncapt (ncid, ncglobal, 'do_aer_columnwf', ncbyte, 1, 0, rcode)
  endif     
  if (do_cld_columnwf) then
     call ncapt (ncid, ncglobal, 'do_cld_columnwf', ncbyte, 1, 1, rcode)
  else
     call ncapt (ncid, ncglobal, 'do_cld_columnwf', ncbyte, 1, 0, rcode)
  endif                            
  if (do_brdf) then
     call ncapt (ncid, ncglobal, 'do_brdf', ncbyte, 1, 1, rcode)
  else
     call ncapt (ncid, ncglobal, 'do_brdf', ncbyte, 1, 0, rcode)
  endif                            
     
  ! write the list of gases in one string
  write(gasstr, '(10(I2,A1,A4,A1))') (i,':',gases(i),',', i=1, ngas)
  nlen=len_trim(gasstr); gasstr=gasstr(1:nlen-1); nlen=nlen-1
  call ncaptc (ncid, ncglobal, 'gases', ncchar,  nlen, trim(gasstr), rcode)

  ! Units for variables
  call ncaptc (ncid, ncglobal, 'windspeed_units', ncchar, 3, 'm/s', rcode)
  call ncaptc (ncid, ncglobal, 'lonlat_units',    ncchar, 7, 'degrees', rcode)
  call ncaptc (ncid, ncglobal, 'zcldtop_units',   ncchar, 2, 'km', rcode)
  call ncaptc (ncid, szaid,    'units', ncchar, 7, 'degrees', rcode)
  call ncaptc (ncid, vzaid,    'units', ncchar, 7, 'degrees', rcode)
  call ncaptc (ncid, azaid,    'units', ncchar, 7, 'degrees', rcode)
  call ncaptc (ncid, lvlid,    'units', ncchar, 2, 'km', rcode)
  if (use_wavelength) then
     call ncaptc(ncid,wavid,'units', ncchar,2, 'nm',rcode)
  else 
     call ncaptc(ncid,wavid,'units', ncchar,4, 'cm-1',rcode)
  endif
  if (use_solar_photons) then
     irrad_unitc='photons/cm2/nm/s'
  else
     irrad_unitc='W/m2/cm-1'
  endif
  call ncaptc(ncid, irradid, 'units', ncchar, len_trim(irrad_unitc), irrad_unitc, rcode)
  if (do_norm_radiance) then
     rad_unitc='unitless'
  else
     rad_unitc=irrad_unitc//'/sr'
  endif
  call ncaptc(ncid, radid,    'units', ncchar, len_trim(rad_unitc), rad_unitc, rcode) 
  call ncaptc(ncid, gascolid, 'units', ncchar, 13, 'molecules/cm2',            rcode)
  call ncaptc(ncid, airid,    'units', ncchar, 13, 'molecules/cm2',            rcode)
  call ncaptc(ncid, psid,     'units', ncchar, 3,  'hPa',                      rcode)
  call ncaptc(ncid, tsid,     'units', ncchar, 1,  'K',                        rcode)    
          

  !=============================================================================
  ! Get out of 'define' mode, and into 'data' mode
  !=============================================================================
  call ncendf (ncid, rcode)
  
  !=============================================================================
  ! Fill in the values of the non time varying variables.
  ! Remember, we're still in "initialization" here - this is only done the
  ! first time through.
  !=============================================================================
  call ncvpt (ncid, szaid, 1, nsza, real(sza, kind=4), rcode)
  call ncvpt (ncid, vzaid, 1, nvza, real(vza, kind=4), rcode)
  call ncvpt (ncid, azaid, 1, naza, real(aza, kind=4), rcode)
  call ncvpt (ncid, lvlid, 1, nz+1, real(zs, kind=4), rcode)
  call ncvpt (ncid, wavid, 1, nw,   real(ws, kind=4), rcode)

  do i = 1, ngas
     gas_indices(i) = i
  enddo
  call ncvpt (ncid, gasid, 1, ngas, gas_indices, rcode)
    
  !=============================================================================
  ! Define the START and COUNT arrays for each of the array variables.
  ! Fill in the values for other variables
  !=============================================================================
  ! write ps and ts
  ndimstart1 = (/ 1 /)
  ndimcount1 = (/ nz+1 /)  
  call ncvpt (ncid, psid, ndimstart1, ndimcount1, real(ps, kind=4), rcode)
  call ncvpt (ncid, tsid, ndimstart1, ndimcount1, real(ts, kind=4), rcode)

  ! write 1D with laydim aircol, aods0, cods0 
  ndimstart1 = (/ 1 /)
  ndimcount1 = (/ nz /)  
  call ncvpt (ncid, airid,  ndimstart1, ndimcount1, real(aircol, kind=4), rcode)
  call ncvpt (ncid, aer0id, ndimstart1, ndimcount1, real(aods0, kind=4),  rcode)
  call ncvpt (ncid, cld0id, ndimstart1, ndimcount1, real(cods0, kind=4),  rcode)

  ! write 1D variables with wavdim
  ndimstart1 = (/ 1 /)
  ndimcount1 = (/ nw /)  
  call ncvpt (ncid, irradid, ndimstart1, ndimcount1, real(irradiance, kind=4), rcode)
  call ncvpt (ncid, sfcid,   ndimstart1, ndimcount1, real(surfalb, kind=4),    rcode)

  ! write 2D variables with wavdim,geodim
  ndimstart2 = (/ 1, 1 /)
  ndimcount2 = (/ nw, ngeom /)  
  call ncvpt (ncid, radid,   ndimstart2, ndimcount2, real(radiance,    kind=4), rcode)
  call ncvpt (ncid, fluxid,  ndimstart2, ndimcount2, real(flux,        kind=4), rcode)
  call ncvpt (ncid, dfluxid, ndimstart2, ndimcount2, real(direct_flux, kind=4), rcode)
  if (do_vector_calc .and. do_QU_output) then
     call ncvpt (ncid, qid,      ndimstart2, ndimcount2, real(q,            kind=4), rcode)
     call ncvpt (ncid, uid,      ndimstart2, ndimcount2, real(u,            kind=4), rcode)
     call ncvpt (ncid, qfluxid,  ndimstart2, ndimcount2, real(qflux,        kind=4), rcode)
     call ncvpt (ncid, ufluxid,  ndimstart2, ndimcount2, real(uflux,        kind=4), rcode)
     call ncvpt (ncid, qdfluxid, ndimstart2, ndimcount2, real(qdirect_flux, kind=4), rcode)
     call ncvpt (ncid, udfluxid, ndimstart2, ndimcount2, real(udirect_flux, kind=4), rcode)
  endif

  if (do_Jacobian) then
     if (use_lambertian) then
        call ncvpt (ncid, sfcwfid, ndimstart2, ndimcount2, real(SurfalbJacobian, kind=4), rcode)
     else 
        call ncvpt (ncid, wswfid, ndimstart2, ndimcount2, real(WSJacobian, kind=4), rcode)
     endif
     if (do_cfrac_Jacobian) &
          call ncvpt (ncid, cfracwfid,  ndimstart2, ndimcount2, real(cfracJacobian, kind=4),         rcode)
     if (do_sfcprs_Jacobian) &
          call ncvpt (ncid, sfcprswfid, ndimstart2, ndimcount2, real(sfcprsJacobian, kind=4),        rcode)
     if (do_aer_columnwf .and. do_aod_Jacobian) &
          call ncvpt (ncid, aodwfid,    ndimstart2, ndimcount2, real(aodJacobian, kind=4),  rcode)
     if (do_aer_columnwf .and. do_assa_Jacobian) &
          call ncvpt (ncid, assawfid,   ndimstart2, ndimcount2, real(assaJacobian, kind=4), rcode)
     if (do_cld_columnwf .and. do_cod_Jacobian) &
          call ncvpt (ncid, codwfid,    ndimstart2, ndimcount2, real(codJacobian, kind=4),  rcode)
     if (do_cld_columnwf .and. do_cssa_Jacobian) &
          call ncvpt (ncid, cssawfid,   ndimstart2, ndimcount2, real(cssaJacobian, kind=4), rcode)
  endif

  if (do_QU_Jacobian) then
     if (use_lambertian) then
        call ncvpt (ncid, sfcqwfid, ndimstart2, ndimcount2, real(SurfalbQJacobian, kind=4), rcode)
     else 
        call ncvpt (ncid, wsqwfid,  ndimstart2, ndimcount2, real(WSQJacobian, kind=4),      rcode)
     endif
     if (do_cfrac_Jacobian) &
          call ncvpt (ncid, cfracqwfid, ndimstart2, ndimcount2, real(cfracQJacobian, kind=4), rcode)
     if (do_sfcprs_Jacobian) &
          call ncvpt (ncid, sfcprsqwfid, ndimstart2, ndimcount2, real(sfcprsQJacobian, kind=4), rcode)
     if (do_aer_columnwf .and. do_aod_Jacobian) &
          call ncvpt (ncid, aodqwfid, ndimstart2, ndimcount2, real(aodQJacobian, kind=4), rcode)
     if (do_aer_columnwf .and. do_assa_Jacobian) &
          call ncvpt (ncid, assaqwfid, ndimstart2, ndimcount2, real(assaQJacobian, kind=4), rcode)
     if (do_cld_columnwf .and. do_cod_Jacobian) &
          call ncvpt (ncid, codqwfid, ndimstart2, ndimcount2, real(codQJacobian, kind=4), rcode)
     if (do_cld_columnwf .and. do_cssa_Jacobian) &
          call ncvpt (ncid, cssaqwfid, ndimstart2, ndimcount2, real(cssaQJacobian, kind=4), rcode)

     if (use_lambertian) then
        call ncvpt (ncid, sfcuwfid, ndimstart2, ndimcount2, real(SurfalbUJacobian, kind=4), rcode)
     else 
        call ncvpt (ncid, wsuwfid, ndimstart2, ndimcount2, real(WSUJacobian, kind=4), rcode)
     endif
     if (do_cfrac_Jacobian) &
          call ncvpt (ncid, cfracuwfid, ndimstart2, ndimcount2, real(cfracUJacobian, kind=4), rcode)
     if (do_sfcprs_Jacobian) &
          call ncvpt (ncid, sfcprsuwfid, ndimstart2, ndimcount2, real(sfcprsUJacobian, kind=4), rcode)
     if (do_aer_columnwf .and. do_aod_Jacobian) &
          call ncvpt (ncid, aoduwfid, ndimstart2, ndimcount2, real(aodUJacobian, kind=4), rcode)
     if (do_aer_columnwf .and. do_assa_Jacobian) &
          call ncvpt (ncid, assauwfid, ndimstart2, ndimcount2, real(assaUJacobian, kind=4), rcode)
     if (do_cld_columnwf .and. do_cod_Jacobian) &
          call ncvpt (ncid, coduwfid, ndimstart2, ndimcount2, real(codUJacobian, kind=4), rcode)
     if (do_cld_columnwf .and. do_cssa_Jacobian) &
          call ncvpt (ncid, cssauwfid, ndimstart2, ndimcount2, real(cssaUJacobian, kind=4), rcode)

  endif

  ! 2D, variables, laydim, gasdim
  ndimstart2 = (/ 1, 1 /)
  ndimcount2 = (/ nz, ngas /)  
  call ncvpt (ncid, gascolid, ndimstart2, ndimcount2, gascol, rcode)

  ! 2D, variables, wavdim, laydim
  ndimstart2 = (/ 1, 1 /)
  ndimcount2 = (/ nw, nz /) 
  call ncvpt (ncid, odsid,   ndimstart2, ndimcount2, real(ods, kind=4), rcode)
  call ncvpt (ncid, ssasid,  ndimstart2, ndimcount2, real(ssas, kind=4), rcode) 
  call ncvpt (ncid, aodsid,  ndimstart2, ndimcount2, real(aods, kind=4), rcode)
  call ncvpt (ncid, assasid, ndimstart2, ndimcount2, real(assas, kind=4), rcode)
  call ncvpt (ncid, codsid,  ndimstart2, ndimcount2, real(cods, kind=4), rcode)
  call ncvpt (ncid, cssasid, ndimstart2, ndimcount2, real(cssas, kind=4), rcode)

  ! 3D variables, wavdim, altdim, geodim
  ndimstart3 = (/ 1, 1, 1 /)
  ndimcount3 = (/ nw, nz, ngeom /) 
  if (do_AMF_calc) then
     call ncvpt (ncid, scatwtid, ndimstart3, ndimcount3, real(ScatWeights, kind=4), rcode)
  endif

  if (do_Jacobian) then
     if (do_T_Jacobian) &
          call ncvpt (ncid, tempwfid, ndimstart3, ndimcount3, real(TJacobian, kind=4), rcode)
     if (.not. do_aer_columnwf .and. do_aod_Jacobian) &
          call ncvpt (ncid, aodwfid, ndimstart3, ndimcount3, real(aodJacobian, kind=4), rcode)
     if (.not. do_aer_columnwf .and. do_assa_Jacobian) &
          call ncvpt (ncid, assawfid, ndimstart3, ndimcount3, real(assaJacobian, kind=4), rcode)
     if (.not. do_cld_columnwf .and. do_cod_Jacobian) &
          call ncvpt (ncid, codwfid, ndimstart3, ndimcount3, real(codJacobian, kind=4), rcode)
     if (.not. do_cld_columnwf .and. do_cssa_Jacobian) &
          call ncvpt (ncid, cssawfid, ndimstart3, ndimcount3, real(cssaJacobian, kind=4), rcode)   
  endif

  if (do_QU_Jacobian) then
     if (.not. do_aer_columnwf .and. do_aod_Jacobian) &
          call ncvpt (ncid, aodqwfid, ndimstart3, ndimcount3, real(aodQJacobian, kind=4), rcode)
     if (.not. do_aer_columnwf .and. do_assa_Jacobian) &
          call ncvpt (ncid, assaqwfid, ndimstart3, ndimcount3, real(assaQJacobian, kind=4), rcode)
     if (.not. do_cld_columnwf .and. do_cod_Jacobian) &
          call ncvpt (ncid, codqwfid, ndimstart3, ndimcount3, real(codQJacobian, kind=4), rcode)
     if (.not. do_cld_columnwf .and. do_cssa_Jacobian) &
          call ncvpt (ncid, cssaqwfid, ndimstart3, ndimcount3, real(cssaQJacobian, kind=4), rcode) 

     if (.not. do_aer_columnwf .and. do_aod_Jacobian) &
          call ncvpt (ncid, aoduwfid, ndimstart3, ndimcount3, real(aodUJacobian, kind=4), rcode)
     if (.not. do_aer_columnwf .and. do_assa_Jacobian) &
          call ncvpt (ncid, assauwfid, ndimstart3, ndimcount3, real(assaUJacobian, kind=4), rcode)
     if (.not. do_cld_columnwf .and. do_cod_Jacobian) &
          call ncvpt (ncid, coduwfid, ndimstart3, ndimcount3, real(codUJacobian, kind=4), rcode)
     if (.not. do_cld_columnwf .and. do_cssa_Jacobian) &
          call ncvpt (ncid, cssauwfid, ndimstart3, ndimcount3, real(cssaUJacobian, kind=4), rcode)    
  endif

  ! 3D variables, wavdim, geodim, gasdim
  if (do_AMF_calc) then
     ndimstart3 = (/ 1, 1, 1/)
     ndimcount3 = (/ nw, ngeom, ngas /) 
     call ncvpt (ncid, amfid, ndimstart3, ndimcount3, real(AMF, kind=4), rcode)
  endif

  ! 4D variables, wavdim, laydim, geodim, gasdim
  ndimstart4 = (/ 1, 1, 1, 1 /)
  ndimcount4 = (/ nw, nz, ngeom, ngas /) 
  if (do_Jacobian) then
     call ncvpt (ncid, gaswfid, ndimstart4, ndimcount4, real(Gas_Jacobian, kind=4), rcode)  
     if (do_QU_Jacobian) then
        call ncvpt (ncid, gasqwfid, ndimstart4, ndimcount4, real(Gas_QJacobian, kind=4), rcode)  
        call ncvpt (ncid, gasuwfid, ndimstart4, ndimcount4, real(Gas_UJacobian, kind=4), rcode)  
     endif
  endif

  ! 4d variables, nsqdim, vzadim, azadim, szadim
  ndimstart4 = (/ 1, 1, 1, 1 /)
  ndimcount4 = (/ nsq, nvza, naza, nsza /) 
  if (do_brdf) then
     call ncvpt (ncid, brdfid, ndimstart4, ndimcount4, real(brdf, kind=4), rcode)  
  endif

  !==============================================================================
  ! CLOSE the NetCDF file
  !==============================================================================
  
  call ncclos (ncid, rcode)
  

return
end subroutine netcdf_wrt
