subroutine netcdf_wrt ( fname,        & !Filename
     nsza, ngeom,                     & !Number of SZA and number of geometries
     nz, ngas,                        & !Number of layers, number of gases
     ods, ssas,                       & !Optical depth and single scattering albedo
     aods, assas,                     & !Aerosol optical depth and single scattering albedo
     cods, cssas,                     & !Cloud optical depth and single scattering albedo
     surfalb,                         & !Surface albedo
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
     use_lambertian,                  &
     irradiance, radiance,            & !Lambda resolution, irradiance and radiance
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
     flux, qflux, uflux,                               & !Total flux
     direct_flux, qdirect_flux,                        & !Direct flux
     udirect_flux,                                     & !Direct flux
     lambda,                                           & !Wavelength
     outlev,                                           & !Output level
     message, fail)   
  
  implicit none
  include '../netcdf/netcdf.inc'

  !=================================  
  !     Input/out variables
  !=================================
  character(len=256), intent(in)                  :: fname                    !Filename
  integer(kind=4),    intent(in)                  :: nz, ngas, lambda, outlev !Number of layers, number of gases, wavelength (index) and  output level (index)
  integer(kind=4),    intent(in)                  :: nsza, ngeom !Number of solar zenith angles, viewing zenith angles, 
                                                                             !azimuth angles and geometries
  real(kind=8), dimension(1), intent(in)          :: surfalb, irradiance     !Surface albedo and irradiance
  real(kind=8), dimension(1,nz), intent(in)       :: aods, assas, cods, &    !Aerosol optical depth, single scattering albedo, cloud optical depth
                                                     cssas, ods, ssas        !cloud single scattering albedo, optical depth and single scattering albedo
  logical, intent(in)  :: do_vector_calc, do_QU_output,        & !Self explained
                          do_Jacobian, do_QU_Jacobian,         &
                          do_AMF_calc, do_T_Jacobian,          &
                          do_sfcprs_Jacobian, do_aod_Jacobian, &
                          do_assa_Jacobian, do_cod_Jacobian,   &
                          do_cssa_Jacobian, do_cfrac_Jacobian, &
                          do_aer_columnwf, do_cld_columnwf,    &
                          use_lambertian

  real(kind=8), dimension(ngeom), intent(in) :: radiance, q, u ! Radiance, q and u values

  real(kind=8), dimension(nsza), intent(in) :: qdirect_flux, udirect_flux, qflux, uflux !Flux and direct flux, q and u values
  real(kind=8), dimension(nsza), intent(in) :: direct_flux, flux  !Flux and direct flux

  real(kind=8), dimension(ngeom), intent(in) :: SurfalbJacobian,  & !Surface albedo Jacobian
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

  real(kind=8), dimension(nz, ngeom), intent(in) :: ScatWeights,   & !Scattering weigths
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

  real(kind=8), dimension(ngeom, ngas), intent(in) :: AMF !Air Mass Factors

  real(kind=8), dimension(nz, ngeom, ngas), intent(in) :: Gas_Jacobian,  & ! Trace gas Jacobians
                                                          Gas_QJacobian, & ! Trace gas QJacobians
                                                          Gas_UJacobian    ! Trace gas UJacobians

 character*(*), intent(inout) :: message
 logical, intent(out)         :: fail
  
  !=================================  
  !     Local variables
  !=================================
  integer :: ncid, rcode
  integer :: szadim, vzadim, azadim, gasdim, laydim, wavdim, geodim, nsqdim
  integer :: radid, qid, uid, irradid, sfcid, sfcwfid, wswfid, cfracwfid, &
       sfcprswfid, aodwfid, assawfid, codwfid, cssawfid, sfcqwfid, wsqwfid, cfracqwfid, &
       sfcprsqwfid, aodqwfid, assaqwfid, codqwfid, cssaqwfid, sfcuwfid, wsuwfid, cfracuwfid,&
       sfcprsuwfid, aoduwfid, assauwfid, coduwfid, cssauwfid, aodsid, assasid, &
       codsid, cssasid, scatwtid, odsid, ssasid, amfid, gaswfid, gasqwfid, &
       gasuwfid, tempwfid, fluxid, dfluxid, qfluxid, ufluxid, qdfluxid, udfluxid
  integer, dimension(1)    :: ndimstart1, ndimcount1
  integer, dimension(2)    :: ndimstart2, ndimcount2
  integer, dimension(3)    :: ndimstart3, ndimcount3
  integer, dimension(4)    :: ndimstart4, ndimcount4 
  integer, dimension(5)    :: ndimstart5, ndimcount5 

  fail = .false.; message = ' '

  ncid = ncopn(trim(fname), ncwrite, rcode)
  if (rcode .eq. -1) then
     message =  ' error in netcdf_out: nccre failed'
     fail = .true.; return
  endif

  !=============================================================================
  ! Create the dependent variables:
  !=============================================================================
  ! varaibles with 1D, wavdim
  irradid = ncvid(ncid, 'irradiance', rcode)
  sfcid   = ncvid(ncid, 'surfalb', rcode)

  ! variables with 2D, wavdim, laydim
  odsid   = ncvid(ncid, 'ods', rcode)
  ssasid  = ncvid(ncid, 'ssas', rcode)
  aodsid  = ncvid(ncid, 'aods', rcode)
  assasid = ncvid(ncid, 'assas', rcode)
  codsid  = ncvid(ncid, 'cods', rcode)
  cssasid = ncvid(ncid, 'cssas', rcode)

  ! variables with 3D, wavdim, geodim, olvdim
  radid   = ncvid(ncid, 'radiance', rcode)
  fluxid  = ncvid(ncid, 'flux', rcode)
  dfluxid = ncvid(ncid, 'direct_flux', rcode)
  if (do_vector_calc .and. do_QU_output) then
     qid      = ncvid(ncid, 'q', rcode)
     uid      = ncvid(ncid, 'u', rcode)
     qfluxid  = ncvid(ncid, 'qflux', rcode)
     ufluxid  = ncvid(ncid, 'uflux', rcode)
     qdfluxid = ncvid(ncid, 'qdirect_flux', rcode)
     udfluxid = ncvid(ncid, 'udirect_flux', rcode)
  endif

  if (do_Jacobian) then

     if (use_lambertian) then
        sfcwfid = ncvid(ncid, 'surfalb_jac', rcode)
     else 
        wswfid  = ncvid(ncid, 'windspeed_jac', rcode)
     endif

     if (do_cfrac_Jacobian) &
          cfracwfid  = ncvid(ncid, 'cfrac_jac', rcode)
     if (do_sfcprs_Jacobian) &
          sfcprswfid = ncvid(ncid, 'sfcprs_jac', rcode)

  endif

  if (do_QU_Jacobian) then

     if (use_lambertian) then
        sfcqwfid = ncvid(ncid, 'surfalb_qjac', rcode)
        sfcuwfid = ncvid(ncid, 'surfalb_ujac', rcode)
     else 
        wsqwfid = ncvid(ncid, 'windspeed_qjac', rcode)
        wsuwfid = ncvid(ncid, 'windspeed_ujac', rcode)
     endif

     if (do_cfrac_Jacobian) then
          cfracqwfid = ncvid(ncid, 'cfrac_qjac', rcode)
          cfracuwfid = ncvid(ncid, 'cfrac_ujac', rcode)
     endif
     if (do_sfcprs_Jacobian) then
          sfcprsqwfid = ncvid(ncid, 'sfcprs_qjac', rcode)
          sfcprsuwfid = ncvid(ncid, 'sfcprs_ujac', rcode)
     endif

  endif

  ! variables with 4D, wavdim, laydim, geodim, olvdim
  if (do_AMF_calc) then
     scatwtid = ncvid(ncid, 'scatweights', rcode)
  endif

  if (do_Jacobian) then
     if (do_T_Jacobian) &
          tempwfid = ncvid(ncid, 't_jac', rcode)
     if (.not. do_aer_columnwf .and. do_aod_Jacobian) &
          aodwfid  = ncvid(ncid, 'aod_jac', rcode)  
     if (.not. do_aer_columnwf .and. do_assa_Jacobian) &
          assawfid = ncvid(ncid, 'assa_jac', rcode)   
     if (.not. do_cld_columnwf .and. do_cod_Jacobian) &
          codwfid  = ncvid(ncid, 'cod_jac', rcode)  
     if (.not. do_cld_columnwf .and. do_cssa_Jacobian) &
          cssawfid = ncvid(ncid, 'cssa_jac', rcode)  
  endif

  if (do_QU_Jacobian) then
     if (.not. do_aer_columnwf .and. do_aod_Jacobian) then 
          aodqwfid = ncvid(ncid, 'aod_qjac', rcode)  
          aoduwfid = ncvid(ncid, 'aod_ujac', rcode)  
     endif
     if (.not. do_aer_columnwf .and. do_assa_Jacobian) then
          assaqwfid = ncvid(ncid, 'assa_qjac', rcode)   
          assauwfid = ncvid(ncid, 'assa_ujac', rcode)   
     endif
     if (.not. do_cld_columnwf .and. do_cod_Jacobian) then
          codqwfid = ncvid(ncid, 'cod_qjac', rcode)  
          coduwfid = ncvid(ncid, 'cod_ujac', rcode)  
     endif
     if (.not. do_cld_columnwf .and. do_cssa_Jacobian) then
          cssaqwfid = ncvid(ncid, 'cssa_qjac', rcode)  
          cssauwfid = ncvid(ncid, 'cssa_ujac', rcode) 
     endif
  endif

  ! variables with 4D, wavdim, geodim, gasdim, olvdim
  if (do_AMF_calc) then
     amfid = ncvid(ncid, 'amf', rcode)
  endif

  ! variables with 5D, wavdim, laydim, geodim, gasdim, ovdim
  if (do_Jacobian) then
     gaswfid = ncvid(ncid,  'gas_jac', rcode)
     if (do_QU_Jacobian) then
        gasqwfid = ncvid(ncid,  'gas_qjac', rcode)
        gasuwfid = ncvid(ncid,  'gas_ujac', rcode)
     endif
  endif

  ! write 1D variables with wavdim
  ndimstart1 = (/ lambda /)
  ndimcount1 = (/ 1 /)  
  call ncvpt (ncid, irradid, ndimstart1, ndimcount1, real(irradiance(1), kind=4), rcode)
  call ncvpt (ncid, sfcid,   ndimstart1, ndimcount1, real(surfalb(1), kind=4),    rcode)

  ! write 3D variables with wavdim,geodim,olvdim
  ndimstart3 = (/ lambda, 1 , outlev/)
  ndimcount3 = (/  1, ngeom , 1/)  
  call ncvpt (ncid, radid,   ndimstart3, ndimcount3, real(radiance(1:ngeom), kind=4), rcode)
  ndimstart3 = (/ lambda, 1 , outlev/)
  ndimcount3 = (/   1, nsza , 1/)  
  call ncvpt (ncid, fluxid,  ndimstart3, ndimcount3, real(flux(1:nsza), kind=4), rcode)
  call ncvpt (ncid, dfluxid, ndimstart3, ndimcount3, real(direct_flux(1:nsza), kind=4), rcode)
  if (do_vector_calc .and. do_QU_output) then
     ndimstart3 = (/ lambda, 1 , outlev/)
     ndimcount3 = (/ 1 , ngeom , 1/)
     call ncvpt (ncid, qid, ndimstart3, ndimcount3, real(q(1:ngeom), kind=4), rcode)
     call ncvpt (ncid, uid, ndimstart3, ndimcount3, real(u(1:ngeom), kind=4), rcode)
     ndimstart3 = (/ lambda, 1 , outlev/)
     ndimcount3 = (/   1, nsza , 1/)  
     call ncvpt (ncid, qfluxid,  ndimstart3, ndimcount3, real(qflux(1:nsza), kind=4), rcode)
     call ncvpt (ncid, ufluxid,  ndimstart3, ndimcount3, real(uflux(1:nsza), kind=4), rcode)
     call ncvpt (ncid, qdfluxid, ndimstart3, ndimcount3, real(qdirect_flux(1:nsza), kind=4), rcode)
     call ncvpt (ncid, udfluxid, ndimstart3, ndimcount3, real(udirect_flux(1:nsza), kind=4), rcode)
  endif

  ndimstart3 = (/ lambda, 1 , outlev/)
  ndimcount3 = (/  1, ngeom , 1/)  
  if (do_Jacobian) then
     if (use_lambertian) then
        call ncvpt (ncid, sfcwfid, ndimstart3, ndimcount3, real(SurfalbJacobian(1:ngeom), kind=4), rcode)
     else 
        call ncvpt (ncid, wswfid, ndimstart3, ndimcount3, real(WSJacobian(1:ngeom), kind=4), rcode)
     endif
     if (do_cfrac_Jacobian) &
          call ncvpt (ncid, cfracwfid,  ndimstart3, ndimcount3, real(cfracJacobian(1:ngeom), kind=4), rcode)
     if (do_sfcprs_Jacobian) &
          call ncvpt (ncid, sfcprswfid, ndimstart3, ndimcount3, real(sfcprsJacobian(1:ngeom), kind=4), rcode)
  endif

  if (do_QU_Jacobian) then
     if (use_lambertian) then
        call ncvpt (ncid, sfcqwfid, ndimstart3, ndimcount3, real(SurfalbQJacobian(1:ngeom), kind=4), rcode)
     else 
        call ncvpt (ncid, wsqwfid,  ndimstart3, ndimcount3, real(WSQJacobian(1:ngeom), kind=4), rcode)
     endif
     if (do_cfrac_Jacobian) &
          call ncvpt (ncid, cfracqwfid, ndimstart3, ndimcount3, real(cfracQJacobian(1:ngeom), kind=4), rcode)
     if (do_sfcprs_Jacobian) &
          call ncvpt (ncid, sfcprsqwfid, ndimstart3, ndimcount3, real(sfcprsQJacobian(1:ngeom), kind=4), rcode)

     if (use_lambertian) then
        call ncvpt (ncid, sfcuwfid, ndimstart3, ndimcount3, real(SurfalbUJacobian(1:ngeom), kind=4), rcode)
     else 
        call ncvpt (ncid, wsuwfid, ndimstart3, ndimcount3, real(WSUJacobian(1:ngeom), kind=4), rcode)
     endif
     if (do_cfrac_Jacobian) &
          call ncvpt (ncid, cfracuwfid, ndimstart3, ndimcount3, real(cfracUJacobian(1:ngeom), kind=4), rcode)
     if (do_sfcprs_Jacobian) &
          call ncvpt (ncid, sfcprsuwfid, ndimstart3, ndimcount3, real(sfcprsUJacobian(1:ngeom), kind=4), rcode)
  endif

  ! 2D, variables, wavdim, laydim
  ndimstart2 = (/ lambda,  1 /)
  ndimcount2 = (/      1, nz /) 
  call ncvpt (ncid, odsid,   ndimstart2, ndimcount2, real(ods(1,1:nz), kind=4), rcode)
  call ncvpt (ncid, ssasid,  ndimstart2, ndimcount2, real(ssas(1,1:nz), kind=4), rcode)
  call ncvpt (ncid, aodsid,  ndimstart2, ndimcount2, real(aods(1,1:nz), kind=4), rcode)
  call ncvpt (ncid, assasid, ndimstart2, ndimcount2, real(assas(1,1:nz), kind=4), rcode)
  call ncvpt (ncid, codsid,  ndimstart2, ndimcount2, real(cods(1,1:nz), kind=4), rcode)
  call ncvpt (ncid, cssasid, ndimstart2, ndimcount2, real(cssas(1,1:nz), kind=4), rcode)

  ! 4D variables, wavdim, altdim, geodim, olvdim
  ndimstart4 = (/ lambda, 1, 1 , outlev/)
  ndimcount4 = (/ 1, nz, ngeom , 1/) 
  if (do_AMF_calc) then
     call ncvpt (ncid, scatwtid, ndimstart4, ndimcount4, real(ScatWeights(1:nz,1:ngeom), kind=4), rcode)
  endif

  if (do_Jacobian) then
     if (do_T_Jacobian) &
          call ncvpt (ncid, tempwfid, ndimstart4, ndimcount4, real(TJacobian(1:nz,1:ngeom), kind=4), rcode)
     if (.not. do_aer_columnwf .and. do_aod_Jacobian) &
          call ncvpt (ncid, aodwfid, ndimstart4, ndimcount4, real(aodJacobian(1:nz,1:ngeom), kind=4), rcode)
     if (.not. do_aer_columnwf .and. do_assa_Jacobian) &
          call ncvpt (ncid, assawfid, ndimstart4, ndimcount4, real(assaJacobian(1:nz,1:ngeom), kind=4), rcode)
     if (.not. do_cld_columnwf .and. do_cod_Jacobian) &
          call ncvpt (ncid, codwfid, ndimstart4, ndimcount4, real(codJacobian(1:nz,1:ngeom), kind=4), rcode)
     if (.not. do_cld_columnwf .and. do_cssa_Jacobian) &
          call ncvpt (ncid, cssawfid, ndimstart4, ndimcount4, real(cssaJacobian(1:nz,1:ngeom), kind=4), rcode)   
  endif

  if (do_QU_Jacobian) then
     if (.not. do_aer_columnwf .and. do_aod_Jacobian) &
          call ncvpt (ncid, aodqwfid, ndimstart4, ndimcount4, real(aodQJacobian(1:nz,1:ngeom), kind=4), rcode)
     if (.not. do_aer_columnwf .and. do_assa_Jacobian) &
          call ncvpt (ncid, assaqwfid, ndimstart4, ndimcount4, real(assaQJacobian(1:nz,1:ngeom), kind=4), rcode)
     if (.not. do_cld_columnwf .and. do_cod_Jacobian) &
          call ncvpt (ncid, codqwfid, ndimstart4, ndimcount4, real(codQJacobian(1:nz,1:ngeom), kind=4), rcode)
     if (.not. do_cld_columnwf .and. do_cssa_Jacobian) &
          call ncvpt (ncid, cssaqwfid, ndimstart4, ndimcount4, real(cssaQJacobian(1:nz,1:ngeom), kind=4), rcode) 

     if (.not. do_aer_columnwf .and. do_aod_Jacobian) &
          call ncvpt (ncid, aoduwfid, ndimstart4, ndimcount4, real(aodUJacobian(1:nz,1:ngeom), kind=4), rcode)
     if (.not. do_aer_columnwf .and. do_assa_Jacobian) &
          call ncvpt (ncid, assauwfid, ndimstart4, ndimcount4, real(assaUJacobian(1:nz,1:ngeom), kind=4), rcode)
     if (.not. do_cld_columnwf .and. do_cod_Jacobian) &
          call ncvpt (ncid, coduwfid, ndimstart4, ndimcount4, real(codUJacobian(1:nz,1:ngeom), kind=4), rcode)
     if (.not. do_cld_columnwf .and. do_cssa_Jacobian) &
          call ncvpt (ncid, cssauwfid, ndimstart4, ndimcount4, real(cssaUJacobian(1:nz,1:ngeom), kind=4), rcode)    
  endif

  ! 4D variables, wavdim, geodim, gasdim, olvdim
  if (do_AMF_calc) then
     ndimstart4 = (/ lambda, 1, 1, outlev/)
     ndimcount4 = (/ 1, ngeom, ngas, 1 /) 
     call ncvpt (ncid, amfid, ndimstart4, ndimcount4, real(AMF(1:ngeom,1:ngas), kind=4), rcode)
  endif

  ! 5D variables, wavdim, laydim, geodim, gasdim, olvdim
  ndimstart5 = (/ lambda, 1, 1, 1, outlev /)
  ndimcount5 = (/ 1, nz, ngeom, ngas, 1 /) 
  if (do_Jacobian) then
     call ncvpt (ncid, gaswfid, ndimstart5, ndimcount5, real(Gas_Jacobian(1:nz,1:ngeom,1:ngas), kind=4), rcode)  
     if (do_QU_Jacobian) then
        call ncvpt (ncid, gasqwfid, ndimstart5, ndimcount5, real(Gas_QJacobian(1:nz,1:ngeom,1:ngas), kind=4), rcode)  
        call ncvpt (ncid, gasuwfid, ndimstart5, ndimcount5, real(Gas_UJacobian(1:nz,1:ngeom,1:ngas), kind=4), rcode)  
     endif
  endif

  !==============================================================================
  ! CLOSE the NetCDF file
  !==============================================================================
  
  call ncclos (ncid, rcode)
  

return
end subroutine netcdf_wrt
