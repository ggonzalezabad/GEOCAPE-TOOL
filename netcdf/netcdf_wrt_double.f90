subroutine netcdf_wrt ( fname, lat, lon, year, month, day, utc,                      &
     sza, vza, aza, ulvl, nw, ws, nz, zs, ps, ts, ngas, gases, gascol, aircol,       &
     cfrac, zcldtop, cldalb, aods0, cods0, ods, ssas, aods, assas, cods, cssas, surfalb, &
     windspeed, nstreams, aercld_nmoms, lam_dw, lam_dfw, use_wavelength, aer_reflam, &
     cld_reflam, do_vector_calc, do_QU_output, do_Jacobian, do_QU_Jacobian,          &
     do_AMF_calc, do_T_Jacobian, do_sfcprs_Jacobian, do_aod_Jacobian,         &
     do_assa_Jacobian, do_cod_Jacobian, do_cssa_Jacobian, do_cfrac_Jacobian,         &
     do_aer_columnwf, do_cld_columnwf, do_norm_WFoutput, do_norm_radiance,           &
     use_lambertian, do_lambertian_cld, do_upwelling, do_effcrs, use_solar_photons,  &
     lam_reso, irradiance, radiance, q, u, Gas_Jacobian, ScatWeights, AMF,           &
     Gas_QJacobian, Gas_UJacobian, TJacobian, SurfalbJacobian, SurfalbQJacobian,     &
     SurfalbUJacobian, WSJacobian, WSQJacobian, WSUJacobian, aodJacobian,            &
     aodQJacobian, aodUJacobian, assaJacobian, assaQJacobian, assaUJacobian,         &
     codJacobian, codQJacobian, codUJacobian, cssaJacobian, cssaQJacobian,           &
     cssaUJacobian, cfracJacobian, cfracQJacobian, cfracUJacobian,                   &
     sfcprsJacobian, sfcprsQJacobian, sfcprsUJacobian, message, fail)   
  
  implicit none
  include '../netcdf/netcdf.inc'

  !=================================  
  !     Input/out variables
  !=================================
  character(len=256), intent(in) :: fname
  integer(kind=4), intent(in)    :: nw, nz, year, month, day, ngas, nstreams, aercld_nmoms
  character(len=4), dimension(ngas) :: gases
  real(kind=8), intent(in)       :: lat, lon, utc, sza, vza, aza, ulvl, cfrac, zcldtop,&
       cldalb, windspeed, lam_reso, lam_dw, lam_dfw, aer_reflam, cld_reflam
  logical, intent(in)            :: do_vector_calc, do_QU_output,     &
       do_Jacobian, do_QU_Jacobian, do_AMF_calc, do_T_Jacobian,          &
       do_sfcprs_Jacobian, do_aod_Jacobian, do_assa_Jacobian, do_cod_Jacobian, &
       do_cssa_Jacobian, do_cfrac_Jacobian, do_aer_columnwf, do_cld_columnwf,    &
       do_norm_WFoutput, do_norm_radiance, use_lambertian,             &
       do_lambertian_cld, do_upwelling, do_effcrs, use_wavelength, use_solar_photons
  real(kind=8), dimension(0:nz), intent(in) :: zs, ps, ts
  real(kind=8), dimension(nz), intent(in)   :: aods0, cods0, aircol
  real(kind=8), dimension(nw), intent(in)   :: ws, surfalb, radiance, irradiance, q, u, &
       SurfalbJacobian, SurfalbQJacobian, SurfalbUJacobian, WSJacobian, WSQJacobian,    &
       WSUJacobian, cfracJacobian, cfracQJacobian, cfracUJacobian, sfcprsJacobian,      &
       sfcprsQJacobian, sfcprsUJacobian
  real(kind=8), dimension(nz, ngas), intent(in) :: gascol
  real(kind=8), dimension(nw, nz), intent(in)   :: aods, assas, cods, cssas, ods, ssas
  real(kind=8), dimension(nw, nz), intent(in)   :: ScatWeights, TJacobian, aodJacobian, &
       aodQJacobian, aodUJacobian, assaJacobian, assaQJacobian, assaUJacobian, codJacobian, &
       codQJacobian, codUJacobian, cssaJacobian, cssaQJacobian, cssaUJacobian
  real(kind=8), dimension(nw, ngas), intent(in) :: AMF
  real(kind=8), dimension(nw, nz, ngas), intent(in) :: Gas_Jacobian, Gas_QJacobian, &
       Gas_UJacobian

 character*(*), intent(inout) :: message
 logical, intent(out)         :: fail
  
  !=================================  
  !     Local variables
  !=================================
  integer :: ncid, rcode, nlen, i
  integer :: szadim, vzadim, azadim, gasdim, lvldim, laydim, wavdim
  integer :: szaid, vzaid, azaid, gasid, lvlid, wavid, psid, tsid, airid, &
       aer0id, cld0id, radid, qid, uid, irradid, sfcid, sfcwfid, wswfid, cfracwfid, &
       sfcprswfid, aodwfid, assawfid, codwfid, cssawfid, sfcqwfid, wsqwfid, cfracqwfid, &
       sfcprsqwfid, aodqwfid, assaqwfid, codqwfid, cssaqwfid, sfcuwfid, wsuwfid, cfracuwfid,&
       sfcprsuwfid, aoduwfid, assauwfid, coduwfid, cssauwfid, gascolid, aodsid, assasid, &
       codsid, cssasid, scatwtid, odsid, ssasid, tempwdid, amfid, gaswfid, gasqwfid, &
       gasuwfid, tempwfid
  integer, dimension(2)    :: gascol_dims, wavalt_dims, wavgas_dims
  integer, dimension(3)    :: wavaltgas_dims
  integer, dimension(ngas) :: gas_indices
  integer, dimension(1)    :: ndimstart1, ndimcount1
  integer, dimension(2)    :: ndimstart2, ndimcount2
  integer, dimension(3)    :: ndimstart3, ndimcount3
  character(len=256)       :: gasstr
  character(len=30)        :: irrad_unitc, rad_unitc

  fail = .false.; message = ' '

  ncid = nccre(trim(fname), ncclob, rcode)
  if (rcode .eq. -1) then
     message =  ' error in netcdf_out: nccre failed'
     fail = .true.; return
  endif 
  
  ! Create the dimensions of the dataset:
  szadim  = ncddef (ncid, 'nsza',     1, rcode)
  vzadim  = ncddef (ncid, 'nvza',   1, rcode)
  azadim  = ncddef (ncid, 'naza', 1, rcode)
  gasdim  = ncddef (ncid, 'ngas', ngas, rcode)
  lvldim  = ncddef (ncid, 'nlevel', nz+1, rcode)
  laydim  = ncddef (ncid, 'nlayer', nz, rcode)
  wavdim  = ncddef (ncid, 'nw',  nw,  rcode)
  
  !=============================================================================
  ! Create the coordinate (aka independent) variables:
  !============================================================================= 
  szaid  = ncvdef (ncid, 'Solarzenithangle', ncdouble, 1, szadim,  rcode)
  vzaid  = ncvdef (ncid, 'Viewingzenithangle', ncdouble, 1, vzadim, rcode)
  azaid  = ncvdef (ncid, 'Relativeazimuthangle', ncdouble, 1, azadim, rcode)
  gasid  = ncvdef (ncid, 'gas', nclong, 1, gasdim, rcode)
  lvlid  = ncvdef (ncid, 'zs', ncdouble, 1, lvldim,  rcode)
  wavid  = ncvdef (ncid, 'Wavelength', ncdouble, 1, wavdim, rcode)

  !=============================================================================
  ! Create a vector containing the often referred to spacetime coordinates:
  !=============================================================================
  gascol_dims(1) = laydim; gascol_dims(2) = gasdim
  wavalt_dims(1) = wavdim; wavalt_dims(2) = laydim
  wavgas_dims(1) = wavdim; wavgas_dims(2) = gasdim
  wavaltgas_dims(1) = wavdim; wavaltgas_dims(2) = laydim; wavaltgas_dims(3) = gasdim
  
  !=============================================================================
  ! Create the dependent variables:
  !=============================================================================
  ! ps, ts
  psid = ncvdef(ncid, 'ps', ncdouble, 1, lvldim, rcode)
  tsid = ncvdef(ncid, 'ts', ncdouble, 1, lvldim, rcode)

  ! aircol, aods0, cods0
  airid = ncvdef(ncid,  'aircol', ncdouble, 1, laydim, rcode)
  aer0id = ncvdef(ncid, 'aods0',  ncdouble, 1, laydim, rcode)
  cld0id = ncvdef(ncid, 'cods0',  ncdouble, 1, laydim, rcode)

  ! vraibles with 1D, wavdim
  radid = ncvdef(ncid, 'radiance',  ncdouble, 1, wavdim, rcode)
  if (do_vector_calc .and. do_QU_output) then
     qid = ncvdef(ncid, 'q',  ncdouble, 1, wavdim, rcode)
     uid = ncvdef(ncid, 'u',  ncdouble, 1, wavdim, rcode)
  endif
  irradid = ncvdef(ncid, 'irradiance',  ncdouble, 1, wavdim, rcode)
  sfcid = ncvdef(ncid, 'surfalb',  ncdouble, 1, wavdim, rcode)
  if (do_Jacobian) then
     if (use_lambertian) then
        sfcwfid = ncvdef(ncid, 'surfalb_jac',  ncdouble, 1, wavdim, rcode)
     else 
        wswfid = ncvdef(ncid, 'windspeed_jac',  ncdouble, 1, wavdim, rcode)
     endif
     if (do_cfrac_Jacobian) &
          cfracwfid = ncvdef(ncid, 'cfrac_jac',  ncdouble, 1, wavdim, rcode)
     if (do_sfcprs_Jacobian) &
          sfcprswfid = ncvdef(ncid, 'sfcprs_jac',  ncdouble, 1, wavdim, rcode)
     if (do_aer_columnwf .and. do_aod_Jacobian) &
          aodwfid = ncvdef(ncid, 'aodcolwf_jac',  ncdouble, 1, wavdim, rcode)
     if (do_aer_columnwf .and. do_assa_Jacobian) &
          assawfid = ncvdef(ncid, 'assacolwf_jac',  ncdouble, 1, wavdim, rcode)
     if (do_cld_columnwf .and. do_cod_Jacobian) &
          codwfid = ncvdef(ncid, 'codcolwf_jac',  ncdouble, 1, wavdim, rcode)
     if (do_cld_columnwf .and. do_cssa_Jacobian) &
          cssawfid = ncvdef(ncid, 'cssacolwf_jac',  ncdouble, 1, wavdim, rcode)
  endif

  if (do_QU_Jacobian) then
     if (use_lambertian) then
        sfcqwfid = ncvdef(ncid, 'surfalb_qjac',  ncdouble, 1, wavdim, rcode)
     else 
        wsqwfid = ncvdef(ncid, 'windspeed_qjac',  ncdouble, 1, wavdim, rcode)
     endif
     if (do_cfrac_Jacobian) &
          cfracqwfid = ncvdef(ncid, 'cfrac_qjac',  ncdouble, 1, wavdim, rcode)
     if (do_sfcprs_Jacobian) &
          sfcprsqwfid = ncvdef(ncid, 'sfcprs_qjac',  ncdouble, 1, wavdim, rcode)
     if (do_aer_columnwf .and. do_aod_Jacobian) &
          aodqwfid = ncvdef(ncid, 'aodcol_qjac',  ncdouble, 1, wavdim, rcode)
     if (do_aer_columnwf .and. do_assa_Jacobian) &
          assaqwfid = ncvdef(ncid, 'assacol_qjac',  ncdouble, 1, wavdim, rcode)
     if (do_cld_columnwf .and. do_cod_Jacobian) &
          codqwfid = ncvdef(ncid, 'codcol_qjac',  ncdouble, 1, wavdim, rcode)
     if (do_cld_columnwf .and. do_cssa_Jacobian) &
          cssaqwfid = ncvdef(ncid, 'cssacol_qjac',  ncdouble, 1, wavdim, rcode)

     if (use_lambertian) then
        sfcuwfid = ncvdef(ncid, 'surfalb_ujac',  ncdouble, 1, wavdim, rcode)
     else 
        wsuwfid = ncvdef(ncid, 'windspeed_ujac',  ncdouble, 1, wavdim, rcode)
     endif
     if (do_cfrac_Jacobian) &
          cfracuwfid = ncvdef(ncid, 'cfrac_ujac',  ncdouble, 1, wavdim, rcode)
     if (do_sfcprs_Jacobian) &
          sfcprsuwfid = ncvdef(ncid, 'sfcprs_ujac',  ncdouble, 1, wavdim, rcode)
     if (do_aer_columnwf .and. do_aod_Jacobian) &
          aoduwfid = ncvdef(ncid, 'aodcol_ujac',  ncdouble, 1, wavdim, rcode)
     if (do_aer_columnwf .and. do_assa_Jacobian) &
          assauwfid = ncvdef(ncid, 'assacol_ujac',  ncdouble, 1, wavdim, rcode)
     if (do_cld_columnwf .and. do_cod_Jacobian) &
          coduwfid = ncvdef(ncid, 'codcol_ujac',  ncdouble, 1, wavdim, rcode)
     if (do_cld_columnwf .and. do_cssa_Jacobian) &
          cssauwfid = ncvdef(ncid, 'cssacol_ujac',  ncdouble, 1, wavdim, rcode)
  endif

  ! 2d, variables, laydim, gasdim
  gascolid = ncvdef(ncid,  'gascol', ncdouble, 2, gascol_dims, rcode)

  ! 2d, variables, wavdim, laydim
  odsid = ncvdef(ncid,  'ods',  ncdouble, 2, wavalt_dims, rcode)
  ssasid = ncvdef(ncid, 'ssas', ncdouble, 2, wavalt_dims, rcode)
  aodsid = ncvdef(ncid,  'aods',  ncdouble, 2, wavalt_dims, rcode)
  assasid = ncvdef(ncid, 'assas', ncdouble, 2, wavalt_dims, rcode)
  codsid = ncvdef(ncid,  'cods',  ncdouble, 2, wavalt_dims, rcode)
  cssasid = ncvdef(ncid,  'cssas', ncdouble, 2, wavalt_dims, rcode)
  if (do_AMF_calc) then
     scatwtid = ncvdef(ncid,  'scatweights', ncdouble, 2, wavalt_dims, rcode)
  endif
  if (do_Jacobian) then
     if (do_T_Jacobian) &
          tempwfid = ncvdef(ncid,  't_jac', ncdouble, 2, wavalt_dims, rcode)
     if (.not. do_aer_columnwf .and. do_aod_Jacobian) &
          aodwfid = ncvdef(ncid,  'aod_jac', ncdouble, 2, wavalt_dims, rcode)  
     if (.not. do_aer_columnwf .and. do_assa_Jacobian) &
          assawfid = ncvdef(ncid,  'assa_jac', ncdouble, 2, wavalt_dims, rcode)   
     if (.not. do_cld_columnwf .and. do_cod_Jacobian) &
          codwfid = ncvdef(ncid,  'cod_jac', ncdouble, 2, wavalt_dims, rcode)  
     if (.not. do_cld_columnwf .and. do_cssa_Jacobian) &
          cssawfid = ncvdef(ncid,  'cssa_jac', ncdouble, 2, wavalt_dims, rcode)  
  endif
  if (do_QU_Jacobian) then
     if (.not. do_aer_columnwf .and. do_aod_Jacobian) &
          aodqwfid = ncvdef(ncid,  'aod_qjac', ncdouble, 2, wavalt_dims, rcode)  
     if (.not. do_aer_columnwf .and. do_assa_Jacobian) &
          assaqwfid = ncvdef(ncid,  'assa_qjac', ncdouble, 2, wavalt_dims, rcode)   
     if (.not. do_cld_columnwf .and. do_cod_Jacobian) &
          codqwfid = ncvdef(ncid,  'cod_qjac', ncdouble, 2, wavalt_dims, rcode)  
     if (.not. do_cld_columnwf .and. do_cssa_Jacobian) &
          cssaqwfid = ncvdef(ncid,  'cssa_qjac', ncdouble, 2, wavalt_dims, rcode)  

     if (.not. do_aer_columnwf .and. do_aod_Jacobian) &
          aoduwfid = ncvdef(ncid,  'aod_ujac', ncdouble, 2, wavalt_dims, rcode)  
     if (.not. do_aer_columnwf .and. do_assa_Jacobian) &
          assauwfid = ncvdef(ncid,  'assa_ujac', ncdouble, 2, wavalt_dims, rcode)   
     if (.not. do_cld_columnwf .and. do_cod_Jacobian) &
          coduwfid = ncvdef(ncid,  'cod_ujac', ncdouble, 2, wavalt_dims, rcode)  
     if (.not. do_cld_columnwf .and. do_cssa_Jacobian) &
          cssauwfid = ncvdef(ncid,  'cssa_ujac', ncdouble, 2, wavalt_dims, rcode) 
  endif

  ! 2d, vraibles, wavdim, gasdim
  !if (do_AMF_calc) then
  !   amfid = ncvdef(ncid,  'amf', ncdouble, 2, wavgas_dims, rcode)
  !endif

  ! 3d, vraibles, wavdim, laydim, gasdim
  if (do_Jacobian) then
     gaswfid = ncvdef(ncid,  'gas_jac', ncdouble, 3, wavaltgas_dims, rcode)
     if (do_QU_Jacobian) then
        gasqwfid = ncvdef(ncid,  'gas_qjac', ncdouble, 3, wavaltgas_dims, rcode)
        gasuwfid = ncvdef(ncid,  'gas_ujac', ncdouble, 3, wavaltgas_dims, rcode)
     endif
  endif
  
  !=============================================================================
  ! Assign attributes (meta-data) to the various variables:
  !============================================================================= 
  call ncapt (ncid, ncglobal, 'lon',          ncdouble, 1, lon,          rcode)
  call ncapt (ncid, ncglobal, 'lat',          ncdouble, 1, lat,          rcode)
  call ncapt (ncid, ncglobal, 'month',        nclong, 1, month,        rcode)
  call ncapt (ncid, ncglobal, 'year',         nclong, 1, year,         rcode)
  call ncapt (ncid, ncglobal, 'day',          nclong, 1, day,          rcode)
  call ncapt (ncid, ncglobal, 'utc',          ncdouble, 1, utc,          rcode)
  call ncapt (ncid, ncglobal, 'cfrac',        ncdouble, 1, cfrac,        rcode)
  call ncapt (ncid, ncglobal, 'zcldtop',      ncdouble, 1, zcldtop,      rcode)
  call ncapt (ncid, ncglobal, 'cldalb',       ncdouble, 1, cldalb,       rcode)
  call ncapt (ncid, ncglobal, 'userlvl',      ncdouble, 1, ulvl,         rcode)
  call ncapt (ncid, ncglobal, 'windspeed',    ncdouble, 1, windspeed,    rcode)
  call ncapt (ncid, ncglobal, 'lam_reso',     ncdouble, 1, lam_reso,     rcode)
  call ncapt (ncid, ncglobal, 'aer_reflam',   ncdouble, 1, aer_reflam,   rcode)
  call ncapt (ncid, ncglobal, 'cld_reflam',   ncdouble, 1, cld_reflam,   rcode)
  call ncapt (ncid, ncglobal, 'lam_dw',       ncdouble, 1, lam_dw,       rcode)
  call ncapt (ncid, ncglobal, 'lam_dfw',      ncdouble, 1, lam_dfw,      rcode)
  call ncapt (ncid, ncglobal, 'nstreams',     nclong, 1, nstreams,     rcode)
  call ncapt (ncid, ncglobal, 'aercld_nmoms', nclong, 1, aercld_nmoms, rcode)

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
     
  ! write the list of gases in one string
  write(gasstr, '(10(I2,A1,A4,A1))') ((i,':',gases(i),','), i=1, ngas)
  nlen=len_trim(gasstr); gasstr=gasstr(1:nlen-1); nlen=nlen-1
  call ncaptc (ncid, ncglobal, 'gases', ncchar,  nlen, trim(gasstr), rcode)

  ! Units for variables
  call ncaptc (ncid, ncglobal, 'windspeed_units', ncchar, 3, 'm/s', rcode)
  call ncaptc (ncid, ncglobal, 'lonlat_units', ncchar, 7, 'degrees', rcode)
  call ncaptc (ncid, ncglobal, 'zcldtop_units', ncchar, 2, 'km', rcode)
  call ncaptc (ncid, szaid, 'units', ncchar, 7, 'degrees', rcode)
  call ncaptc (ncid, vzaid, 'units', ncchar, 7, 'degrees', rcode)
  call ncaptc (ncid, azaid, 'units', ncchar, 7, 'degrees', rcode)
  call ncaptc (ncid, lvlid, 'units', ncchar, 2, 'km', rcode)
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
     rad_unitc=irrad_unitc // '/sr'
  endif
  call ncaptc(ncid, radid, 'units', ncchar, len_trim(rad_unitc), rad_unitc, rcode)
  
  call ncaptc(ncid, gascolid, 'units', ncchar, 13, 'molecules/cm2', rcode)
  call ncaptc(ncid, airid, 'units', ncchar, 13, 'molecules/cm2', rcode)
  call ncaptc(ncid, psid, 'units', ncchar, 3, 'hPa', rcode)
  call ncaptc(ncid, tsid, 'units', ncchar, 1, 'K', rcode)                

  !=============================================================================
  ! Get out of 'define' mode, and into 'data' mode
  !=============================================================================
  
  call ncendf (ncid, rcode)
  
  !=============================================================================
  ! Fill in the values of the non time varying variables.
  ! Remember, we're still in "initialization" here - this is only done the
  ! first time through.
  !=============================================================================
  call ncvpt (ncid, szaid, 1, 1, sza, rcode)
  call ncvpt (ncid, vzaid, 1, 1, vza, rcode)
  call ncvpt (ncid, azaid, 1, 1, aza, rcode)
  call ncvpt (ncid, lvlid, 1, nz+1, zs, rcode)
  call ncvpt (ncid, wavid, 1, nw, ws, rcode)

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
  call ncvpt (ncid, psid, ndimstart1, ndimcount1, ps, rcode)
  call ncvpt (ncid, tsid, ndimstart1, ndimcount1, ts, rcode)

  ! write aircol, aods0, cods0
  ndimstart1 = (/ 1 /)
  ndimcount1 = (/ nz /)  
  call ncvpt (ncid, airid, ndimstart1, ndimcount1, aircol, rcode)
  call ncvpt (ncid, aer0id, ndimstart1, ndimcount1, aods0, rcode)
  call ncvpt (ncid, cld0id, ndimstart1, ndimcount1, cods0, rcode)

  ! write 1D variables with wavdim
  ndimstart1 = (/ 1 /)
  ndimcount1 = (/ nw /)  
  call ncvpt (ncid, irradid, ndimstart1, ndimcount1, irradiance, rcode)
  call ncvpt (ncid, radid, ndimstart1, ndimcount1, radiance, rcode)
  if (do_vector_calc .and. do_QU_output) then
     call ncvpt (ncid, qid, ndimstart1, ndimcount1, q, rcode)
     call ncvpt (ncid, uid, ndimstart1, ndimcount1, u, rcode)
  endif
  call ncvpt (ncid, sfcid, ndimstart1, ndimcount1, surfalb, rcode)
  if (do_Jacobian) then
     if (use_lambertian) then
        call ncvpt (ncid, sfcwfid, ndimstart1, ndimcount1, SurfalbJacobian, rcode)
     else 
        call ncvpt (ncid, wswfid, ndimstart1, ndimcount1, WSJacobian, rcode)
     endif
     if (do_cfrac_Jacobian) &
          call ncvpt (ncid, cfracwfid, ndimstart1, ndimcount1, cfracJacobian, rcode)
     if (do_sfcprs_Jacobian) &
          call ncvpt (ncid, sfcprswfid, ndimstart1, ndimcount1, sfcprsJacobian, rcode)
     if (do_aer_columnwf .and. do_aod_Jacobian) &
          call ncvpt (ncid, aodwfid, ndimstart1, ndimcount1, aodJacobian(1:nw, 1), rcode)
     if (do_aer_columnwf .and. do_assa_Jacobian) &
          call ncvpt (ncid, assawfid, ndimstart1, ndimcount1, assaJacobian(1:nw, 1), rcode)
     if (do_cld_columnwf .and. do_cod_Jacobian) &
          call ncvpt (ncid, codwfid, ndimstart1, ndimcount1, codJacobian(1:nw, 1), rcode)
     if (do_cld_columnwf .and. do_cssa_Jacobian) &
          call ncvpt (ncid, cssawfid, ndimstart1, ndimcount1, cssaJacobian(1:nw, 1), rcode)
  endif

  if (do_QU_Jacobian) then
     if (use_lambertian) then
        call ncvpt (ncid, sfcqwfid, ndimstart1, ndimcount1, SurfalbQJacobian, rcode)
     else 
        call ncvpt (ncid, wsqwfid, ndimstart1, ndimcount1, WSQJacobian, rcode)
     endif
     if (do_cfrac_Jacobian) &
          call ncvpt (ncid, cfracqwfid, ndimstart1, ndimcount1, cfracQJacobian, rcode)
     if (do_sfcprs_Jacobian) &
          call ncvpt (ncid, sfcprsqwfid, ndimstart1, ndimcount1, sfcprsQJacobian, rcode)
     if (do_aer_columnwf .and. do_aod_Jacobian) &
          call ncvpt (ncid, aodqwfid, ndimstart1, ndimcount1, aodQJacobian(1:nw, 1), rcode)
     if (do_aer_columnwf .and. do_assa_Jacobian) &
          call ncvpt (ncid, assaqwfid, ndimstart1, ndimcount1, assaQJacobian(1:nw, 1), rcode)
     if (do_cld_columnwf .and. do_cod_Jacobian) &
          call ncvpt (ncid, codqwfid, ndimstart1, ndimcount1, codQJacobian(1:nw, 1), rcode)
     if (do_cld_columnwf .and. do_cssa_Jacobian) &
          call ncvpt (ncid, cssaqwfid, ndimstart1, ndimcount1, cssaQJacobian(1:nw, 1), rcode)

     if (use_lambertian) then
        call ncvpt (ncid, sfcuwfid, ndimstart1, ndimcount1, SurfalbUJacobian, rcode)
     else 
        call ncvpt (ncid, wsuwfid, ndimstart1, ndimcount1, WSUJacobian, rcode)
     endif
     if (do_cfrac_Jacobian) &
          call ncvpt (ncid, cfracuwfid, ndimstart1, ndimcount1, cfracUJacobian, rcode)
     if (do_sfcprs_Jacobian) &
          call ncvpt (ncid, sfcprsuwfid, ndimstart1, ndimcount1, sfcprsUJacobian, rcode)
     if (do_aer_columnwf .and. do_aod_Jacobian) &
          call ncvpt (ncid, aoduwfid, ndimstart1, ndimcount1, aodUJacobian(1:nw, 1), rcode)
     if (do_aer_columnwf .and. do_assa_Jacobian) &
          call ncvpt (ncid, assauwfid, ndimstart1, ndimcount1, assaUJacobian(1:nw, 1), rcode)
     if (do_cld_columnwf .and. do_cod_Jacobian) &
          call ncvpt (ncid, coduwfid, ndimstart1, ndimcount1, codUJacobian(1:nw, 1), rcode)
     if (do_cld_columnwf .and. do_cssa_Jacobian) &
          call ncvpt (ncid, cssauwfid, ndimstart1, ndimcount1, cssaUJacobian(1:nw, 1), rcode)

  endif

  ! 2d, variables, laydim, gasdim
  ndimstart2 = (/ 1, 1 /)
  ndimcount2 = (/ nz, ngas /)  
  call ncvpt (ncid, gascolid, ndimstart2, ndimcount2, gascol, rcode)

  ! 2d, variables, wavdim, laydim
  ndimstart2 = (/ 1, 1 /)
  ndimcount2 = (/ nw, nz /) 
  call ncvpt (ncid, odsid, ndimstart2, ndimcount2, ods, rcode)
  call ncvpt (ncid, ssasid, ndimstart2, ndimcount2, ssas, rcode) 
  call ncvpt (ncid, aodsid, ndimstart2, ndimcount2, aods, rcode)
  call ncvpt (ncid, assasid, ndimstart2, ndimcount2, assas, rcode)
  call ncvpt (ncid, codsid, ndimstart2, ndimcount2, cods, rcode)
  call ncvpt (ncid, cssasid, ndimstart2, ndimcount2, cssas, rcode)

  if (do_AMF_calc) then
     call ncvpt (ncid, scatwtid, ndimstart2, ndimcount2, ScatWeights, rcode)
  endif

  if (do_Jacobian) then
     if (do_T_Jacobian) &
          call ncvpt (ncid, tempwfid, ndimstart2, ndimcount2, TJacobian, rcode)
     if (.not. do_aer_columnwf .and. do_aod_Jacobian) &
          call ncvpt (ncid, aodwfid, ndimstart2, ndimcount2, aodJacobian, rcode)
     if (.not. do_aer_columnwf .and. do_assa_Jacobian) &
          call ncvpt (ncid, assawfid, ndimstart2, ndimcount2, assaJacobian, rcode)
     if (.not. do_cld_columnwf .and. do_cod_Jacobian) &
          call ncvpt (ncid, codwfid, ndimstart2, ndimcount2, codJacobian, rcode)
     if (.not. do_cld_columnwf .and. do_cssa_Jacobian) &
          call ncvpt (ncid, cssawfid, ndimstart2, ndimcount2, cssaJacobian, rcode)   
  endif

  if (do_QU_Jacobian) then
     if (.not. do_aer_columnwf .and. do_aod_Jacobian) &
          call ncvpt (ncid, aodqwfid, ndimstart2, ndimcount2, aodQJacobian, rcode)
     if (.not. do_aer_columnwf .and. do_assa_Jacobian) &
          call ncvpt (ncid, assaqwfid, ndimstart2, ndimcount2, assaQJacobian, rcode)
     if (.not. do_cld_columnwf .and. do_cod_Jacobian) &
          call ncvpt (ncid, codqwfid, ndimstart2, ndimcount2, codQJacobian, rcode)
     if (.not. do_cld_columnwf .and. do_cssa_Jacobian) &
          call ncvpt (ncid, cssaqwfid, ndimstart2, ndimcount2, cssaQJacobian, rcode) 

     if (.not. do_aer_columnwf .and. do_aod_Jacobian) &
          call ncvpt (ncid, aoduwfid, ndimstart2, ndimcount2, aodUJacobian, rcode)
     if (.not. do_aer_columnwf .and. do_assa_Jacobian) &
          call ncvpt (ncid, assauwfid, ndimstart2, ndimcount2, assaUJacobian, rcode)
     if (.not. do_cld_columnwf .and. do_cod_Jacobian) &
          call ncvpt (ncid, coduwfid, ndimstart2, ndimcount2, codUJacobian, rcode)
     if (.not. do_cld_columnwf .and. do_cssa_Jacobian) &
          call ncvpt (ncid, cssauwfid, ndimstart2, ndimcount2, cssaUJacobian, rcode)    
  endif

!  ! 2d, vraibles, wavdim, gasdim
!  if (do_AMF_calc) then
!     ndimstart2 = (/ 1, 1/)
!     ndimcount2 = (/ nw, ngas /) 
!     call ncvpt (ncid, amfid, ndimstart2, ndimcount2, AMF, rcode)
!  endif

  ! 3d, vraibles, wavdim, laydim, gasdim
  ndimstart3 = (/ 1, 1, 1 /)
  ndimcount3 = (/ nw, nz, ngas /) 
  if (do_Jacobian) then
     call ncvpt (ncid, gaswfid, ndimstart3, ndimcount3, Gas_Jacobian, rcode)  
     if (do_QU_Jacobian) then
        call ncvpt (ncid, gasqwfid, ndimstart3, ndimcount3, Gas_QJacobian, rcode)  
        call ncvpt (ncid, gasuwfid, ndimstart3, ndimcount3, Gas_UJacobian, rcode)  
     endif
  endif
  
  !==============================================================================
  ! CLOSE the NetCDF file
  !==============================================================================
  
  call ncclos (ncid, rcode)
  

return
end subroutine netcdf_wrt
