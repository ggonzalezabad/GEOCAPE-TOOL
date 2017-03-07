subroutine netcdf_wrt ( fname, message, fail)   

  USE GC_variables_module, ONLY: GC_n_sun_positions, VLIDORT_Out, & !# SZA, # geometries (TS_N_GEOMETRIES)
       GC_nlayers, ngases, nlambdas, nclambdas, & !# layers, gases, wavelengths, covolved wavelengths if needed
       GC_n_user_levels, & !# of output levels
       opdeps, ssalbs, & !optical depth and single scattering albedo
       aer_opdeps, aer_ssalbs, & !Aerosol optical depth and single scattering albedo
       cld_opdeps, cld_ssalbs, & !Cloud optical depth and single scattering albedo
       ground_ler, & !Surface albedo
       do_vector_calculation, do_StokesQU_output, do_Jacobians,      & !LOGICAL SWITCHES
       do_QU_Jacobians, do_AMF_calculation, do_T_Jacobians, do_sfcprs_Jacobians,    & !LOGICAL SWITCHES
       do_aod_Jacobians, do_assa_Jacobians, do_cod_Jacobians, do_cssa_Jacobians,    & !LOGICAL SWITCHES
       do_cfrac_Jacobians, do_aer_columnwf, do_cld_columnwf, use_lambertian,        & !LOGICAL SWITCHES
       do_effcrs, lambda_resolution, & !Convolution switch and lambda resolution
       solar_cspec_data, GC_radiances, & !Irradiance and radiance
       GC_Qvalues, GC_Uvalues, & !Qvalues, Uvalues
       GC_Tracegas_Jacobians, GC_Scattering_Weights, & !Gas Jacobians, Scattering Weights
       GC_AMFs, GC_Tracegas_QJacobians, GC_Tracegas_UJacobians, & !Air mass factor, Gas UJacobian, Gas QJacobian
       GC_Temperature_Jacobians, & !Temperature Jacobians, 
       GC_Surfalbedo_Jacobians, GC_Surfalbedo_QJacobians, GC_Surfalbedo_UJacobians, & !Surface albedo Jacobians
       GC_Windspeed_Jacobians, GC_Windspeed_QJacobians, GC_Windspeed_UJacobians, & !Wind speed Jacobians
       GC_aod_Jacobians, GC_aod_QJacobians, GC_aod_UJacobians, & !Aerosol optical depth Jacobians
       GC_assa_Jacobians, GC_assa_QJacobians, GC_assa_UJacobians, & !Aerosol single sccatering Jacobians
       GC_cod_Jacobians, GC_cod_QJacobians, GC_cod_UJacobians, & !Cloud optical depth Jacobians
       GC_cssa_Jacobians, GC_cssa_QJacobians, GC_cssa_UJacobians, & !Cloud single scattering Jacobians
       GC_cfrac_Jacobians, GC_cfrac_QJacobians, GC_cfrac_UJacobians, & !Cloud fraction Jacobians
       GC_sfcprs_Jacobians, GC_sfcprs_QJacobians, GC_sfcprs_UJacobians, & !Surface or cloud pressure Jacobians
       GC_flux, GC_Qflux, GC_Uflux, & !Totla flux
       GC_direct_flux, GC_Qdirect_flux, GC_Udirect_flux, & !Direct flux
       didx !Direction index
  
  IMPLICIT none
  INCLUDE '../netcdf/netcdf.inc'

  !=================================  
  !     Input/out variables
  !=================================
  character(len=256), intent(in) :: fname                    !Filename
  character*(*), intent(inout)   :: message
  logical, intent(out)           :: fail
  
  !=================================  
  !     Local variables
  !=================================
  integer :: ncid, rcode, nwav, iout, ngeom
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

  rcode = NF_OPEN(trim(fname), NF_WRITE, ncid)
  print*, rcode
  if (rcode .ne. NF_NOERR) then
     message =  ' error in netcdf_out: nccre failed'
     fail = .true.; return
  endif

  !=============================================================================
  ! Create the dependent variables:
  !=============================================================================
  ! varaibles with 1D, wavdim
  rcode = NF_INQ_VARID(ncid, 'irradiance', irradid)
  rcode = NF_INQ_VARID(ncid, 'surfalb', sfcid)

  ! variables with 2D, wavdim, laydim
  rcode = NF_INQ_VARID(ncid, 'ods', odsid)
  rcode = NF_INQ_VARID(ncid, 'ssas', ssasid)
  rcode = NF_INQ_VARID(ncid, 'aods', aodsid)
  rcode = NF_INQ_VARID(ncid, 'assas', assasid)
  rcode = NF_INQ_VARID(ncid, 'cods', codsid)
  rcode = NF_INQ_VARID(ncid, 'cssas', cssasid)

  ! variables with 3D, wavdim, geodim, olvdim
  rcode = NF_INQ_VARID(ncid, 'radiance', radid)
  rcode = NF_INQ_VARID(ncid, 'flux', fluxid)
  rcode = NF_INQ_VARID(ncid, 'direct_flux', dfluxid)
  if (do_vector_calculation .and. do_StokesQU_output) then
     rcode = NF_INQ_VARID(ncid, 'q', qid)
     rcode = NF_INQ_VARID(ncid, 'u', uid)
     rcode = NF_INQ_VARID(ncid, 'qflux', qfluxid)
     rcode = NF_INQ_VARID(ncid, 'uflux', ufluxid)
     rcode = NF_INQ_VARID(ncid, 'qdirect_flux', qdfluxid)
     rcode = NF_INQ_VARID(ncid, 'udirect_flux', udfluxid)
  endif

  if (do_Jacobians) then

     if (use_lambertian) then
        rcode = NF_INQ_VARID(ncid, 'surfalb_jac', sfcwfid)
     else 
        rcode = NF_INQ_VARID(ncid, 'windspeed_jac', wswfid)
     endif

     if (do_cfrac_Jacobians) &
          rcode = NF_INQ_VARID(ncid, 'cfrac_jac', cfracwfid)
     if (do_sfcprs_Jacobians) &
          rcode = NF_INQ_VARID(ncid, 'sfcprs_jac', sfcprswfid)

  endif

  if (do_QU_Jacobians) then

     if (use_lambertian) then
        rcode = NF_INQ_VARID(ncid, 'surfalb_qjac', sfcqwfid)
        rcode = NF_INQ_VARID(ncid, 'surfalb_ujac', sfcuwfid)
     else 
        rcode = NF_INQ_VARID(ncid, 'windspeed_qjac', wsqwfid)
        rcode = NF_INQ_VARID(ncid, 'windspeed_ujac', wsuwfid)
     endif

     if (do_cfrac_Jacobians) then
        rcode = NF_INQ_VARID(ncid, 'cfrac_qjac', cfracqwfid)
        rcode = NF_INQ_VARID(ncid, 'cfrac_ujac', cfracuwfid)
     endif
     if (do_sfcprs_Jacobians) then
        rcode = NF_INQ_VARID(ncid, 'sfcprs_qjac', sfcprsqwfid)
        rcode = NF_INQ_VARID(ncid, 'sfcprs_ujac', sfcprsuwfid)
     endif

  endif

  ! variables with 4D, wavdim, laydim, geodim, olvdim
  if (do_AMF_calculation) then
     rcode = NF_INQ_VARID(ncid, 'scatweights', scatwtid)
  endif

  if (do_Jacobians) then
     if (do_T_Jacobians) &
          rcode = NF_INQ_VARID(ncid, 't_jac', tempwfid)
     if (.not. do_aer_columnwf .and. do_aod_Jacobians) &
          rcode = NF_INQ_VARID(ncid, 'aod_jac', aodwfid)  
     if (.not. do_aer_columnwf .and. do_assa_Jacobians) &
          rcode = NF_INQ_VARID(ncid, 'assa_jac', assawfid)   
     if (.not. do_cld_columnwf .and. do_cod_Jacobians) &
          rcode = NF_INQ_VARID(ncid, 'cod_jac', codwfid)  
     if (.not. do_cld_columnwf .and. do_cssa_Jacobians) &
          rcode = NF_INQ_VARID(ncid, 'cssa_jac', cssawfid)  
  endif

  if (do_QU_Jacobians) then
     if (.not. do_aer_columnwf .and. do_aod_Jacobians) then 
        rcode = NF_INQ_VARID(ncid, 'aod_qjac', aodqwfid)  
        rcode = NF_INQ_VARID(ncid, 'aod_ujac', aoduwfid)  
     endif
     if (.not. do_aer_columnwf .and. do_assa_Jacobians) then
        rcode = NF_INQ_VARID(ncid, 'assa_qjac', assaqwfid)   
        rcode = NF_INQ_VARID(ncid, 'assa_ujac', assauwfid)   
     endif
     if (.not. do_cld_columnwf .and. do_cod_Jacobians) then
        rcode = NF_INQ_VARID(ncid, 'cod_qjac', codqwfid)  
        rcode = NF_INQ_VARID(ncid, 'cod_ujac', coduwfid)  
     endif
     if (.not. do_cld_columnwf .and. do_cssa_Jacobians) then
        rcode = NF_INQ_VARID(ncid, 'cssa_qjac', cssaqwfid)  
        rcode = NF_INQ_VARID(ncid, 'cssa_ujac', cssauwfid) 
     endif
  endif

  ! variables with 4D, wavdim, geodim, gasdim, olvdim
  if (do_AMF_calculation) then
     rcode = NF_INQ_VARID(ncid, 'amf', amfid)
  endif

  ! variables with 5D, wavdim, laydim, geodim, gasdim, ovdim
  if (do_Jacobians) then
     rcode = NF_INQ_VARID(ncid,  'gas_jac', gaswfid)
     if (do_QU_Jacobians) then
        rcode = NF_INQ_VARID(ncid,  'gas_qjac', gasqwfid)
        rcode = NF_INQ_VARID(ncid,  'gas_ujac', gasuwfid)
     endif
  endif

  ! Work out the number of lambdas. It depends if we are doing convolved cross-sections or not
  IF (.NOT. do_effcrs .AND. lambda_resolution /= 0.0d0 ) THEN 
     nwav = nclambdas
  ELSE
     nwav = nlambdas
  END IF

  ! Fill in value for # of geometries
  ngeom = VLIDORT_Out%Main%TS_N_GEOMETRIES

  ! Loop over lambdas
  print*, 'Writing output #', didx
  
  ! write 1D variables with wavdim
  ndimstart1 = (/ 1 /)
  ndimcount1 = (/ nwav /)  
  rcode = NF_PUT_VARA_REAL (ncid, irradid, ndimstart1, ndimcount1, real(solar_cspec_data(1:nwav), kind=4))
  rcode = NF_PUT_VARA_REAL (ncid, sfcid,   ndimstart1, ndimcount1, real(ground_ler(1:nwav), kind=4))
  
  ! 2D, variables, wavdim, laydim
  ndimstart2 = (/    1,  1 /)
  ndimcount2 = (/ nwav, GC_nlayers /) 
  rcode = NF_PUT_VARA_REAL (ncid, odsid,   ndimstart2, ndimcount2, real(opdeps(1:nwav,1:GC_nlayers), kind=4))
  rcode = NF_PUT_VARA_REAL (ncid, ssasid,  ndimstart2, ndimcount2, real(ssalbs(1:nwav,1:GC_nlayers), kind=4))
  rcode = NF_PUT_VARA_REAL (ncid, aodsid,  ndimstart2, ndimcount2, real(aer_opdeps(1:nwav,1:GC_nlayers), kind=4))
  rcode = NF_PUT_VARA_REAL (ncid, assasid, ndimstart2, ndimcount2, real(aer_ssalbs(1:nwav,1:GC_nlayers), kind=4))
  rcode = NF_PUT_VARA_REAL (ncid, codsid,  ndimstart2, ndimcount2, real(cld_opdeps(1:nwav,1:GC_nlayers), kind=4))
  rcode = NF_PUT_VARA_REAL (ncid, cssasid, ndimstart2, ndimcount2, real(cld_ssalbs(1:nwav,1:GC_nlayers), kind=4))
  
  ! write 3D variables with wavdim,geodim,olvdim
  DO iout = 1, GC_n_user_levels
     ndimstart3 = (/ 1, 1 , iout/)
     ndimcount3 = (/ nwav, ngeom , 1/)  
     rcode = NF_PUT_VARA_REAL (ncid, radid,   ndimstart3, ndimcount3, &
          real(GC_radiances(1:nwav,iout,1:ngeom,didx), kind=4))
     ndimcount3 = (/   1, GC_n_sun_positions , 1/)  
     rcode = NF_PUT_VARA_REAL (ncid, fluxid,  ndimstart3, ndimcount3, &
          real(GC_flux(1:nwav,iout,1:GC_n_sun_positions,didx), kind=4))
     rcode = NF_PUT_VARA_REAL (ncid, dfluxid, ndimstart3, ndimcount3, &
          real(GC_direct_flux(1:nwav,iout,1:GC_n_sun_positions,didx), kind=4))
     if (do_vector_calculation .and. do_StokesQU_output) then
        ndimcount3 = (/ 1 , ngeom , 1/)
        rcode = NF_PUT_VARA_REAL (ncid, qid, ndimstart3, ndimcount3, &
             real(GC_Qvalues(1:nwav,iout,1:ngeom,didx), kind=4))
        rcode = NF_PUT_VARA_REAL (ncid, uid, ndimstart3, ndimcount3, &
             real(GC_Uvalues(1:nwav,iout,1:ngeom,didx), kind=4))
        ndimcount3 = (/   1, GC_n_sun_positions , 1/)  
        rcode = NF_PUT_VARA_REAL (ncid, qfluxid,  ndimstart3, ndimcount3, &
             real(GC_Qflux(1:nwav,iout,1:GC_n_sun_positions,didx), kind=4))
        rcode = NF_PUT_VARA_REAL (ncid, ufluxid,  ndimstart3, ndimcount3, &
             real(GC_Uflux(1:nwav,iout,1:GC_n_sun_positions,didx), kind=4))
        rcode = NF_PUT_VARA_REAL (ncid, qdfluxid, ndimstart3, ndimcount3, &
             real(GC_Qdirect_flux(1:nwav,iout,1:GC_n_sun_positions,didx), kind=4))
        rcode = NF_PUT_VARA_REAL (ncid, udfluxid, ndimstart3, ndimcount3, &
             real(GC_Udirect_flux(1:nwav,iout,1:GC_n_sun_positions,didx), kind=4))
     endif
     
     ndimstart3 = (/ 1, 1 , iout/)
     ndimcount3 = (/ nwav, ngeom , 1/)  
     if (do_Jacobians) then
        if (use_lambertian) then
           rcode = NF_PUT_VARA_REAL (ncid, sfcwfid, ndimstart3, ndimcount3, &
                real(GC_Surfalbedo_Jacobians(1:nwav,iout,1:ngeom,didx), kind=4))
        else 
           rcode = NF_PUT_VARA_REAL (ncid, wswfid, ndimstart3, ndimcount3, &
                real(GC_Windspeed_Jacobians(1:nwav,iout,1:ngeom,didx), kind=4))
        endif
        if (do_cfrac_Jacobians) &
             rcode = NF_PUT_VARA_REAL (ncid, cfracwfid,  ndimstart3, ndimcount3, &
             real(GC_cfrac_Jacobians(1:nwav,iout,1:ngeom,didx), kind=4))
        if (do_sfcprs_Jacobians) &
             rcode = NF_PUT_VARA_REAL (ncid, sfcprswfid, ndimstart3, ndimcount3, &
             real(GC_sfcprs_Jacobians(1:nwav,iout,1:ngeom,didx), kind=4))
     endif
     
     if (do_QU_Jacobians) then
        if (use_lambertian) then
           rcode = NF_PUT_VARA_REAL (ncid, sfcqwfid, ndimstart3, ndimcount3, &
                real(GC_Surfalbedo_QJacobians(1:nwav,iout,1:ngeom,didx), kind=4))
        else 
           rcode = NF_PUT_VARA_REAL (ncid, wsqwfid,  ndimstart3, ndimcount3, &
                real(GC_Windspeed_QJacobians(1:nwav,iout,1:ngeom,didx), kind=4))
        endif
        if (do_cfrac_Jacobians) &
             rcode = NF_PUT_VARA_REAL (ncid, cfracqwfid, ndimstart3, ndimcount3, &
             real(GC_cfrac_QJacobians(1:nwav,iout,1:ngeom,didx), kind=4))
        if (do_sfcprs_Jacobians) &
             rcode = NF_PUT_VARA_REAL (ncid, sfcprsqwfid, ndimstart3, ndimcount3, &
             real(GC_sfcprs_QJacobians(1:nwav,iout,1:ngeom,didx), kind=4))
        if (use_lambertian) then
           rcode = NF_PUT_VARA_REAL (ncid, sfcuwfid, ndimstart3, ndimcount3, &
                real(GC_Surfalbedo_UJacobians(1:nwav,iout,1:ngeom,didx), kind=4))
        else 
           rcode = NF_PUT_VARA_REAL (ncid, wsuwfid, ndimstart3, ndimcount3, &
                real(GC_Windspeed_UJacobians(1:nwav,iout,1:ngeom,didx), kind=4))
        endif
        if (do_cfrac_Jacobians) &
             rcode = NF_PUT_VARA_REAL (ncid, cfracuwfid, ndimstart3, ndimcount3, &
             real(GC_cfrac_UJacobians(1:nwav,iout,1:ngeom,didx), kind=4))
        if (do_sfcprs_Jacobians) &
             rcode = NF_PUT_VARA_REAL (ncid, sfcprsuwfid, ndimstart3, ndimcount3, &
             real(GC_sfcprs_UJacobians(1:nwav,iout,1:ngeom,didx), kind=4))
     endif
     
     ! 4D variables, wavdim, altdim, geodim, olvdim
     ndimstart4 = (/    1,          1,     1, iout/)
     ndimcount4 = (/ nwav, GC_nlayers, ngeom,    1/) 
     if (do_AMF_calculation) then
        rcode = NF_PUT_VARA_REAL (ncid, scatwtid, ndimstart4, ndimcount4, &
             real(GC_Scattering_Weights(1:nwav,1:GC_nlayers,iout,1:ngeom,didx), kind=4))
     endif
     
     if (do_Jacobians) then
        if (do_T_Jacobians) &
             rcode = NF_PUT_VARA_REAL (ncid, tempwfid, ndimstart4, ndimcount4, &
             real(GC_Temperature_Jacobians(1:nwav,1:GC_nlayers,iout,1:ngeom,didx), kind=4))
        if (.not. do_aer_columnwf .and. do_aod_Jacobians) &
             rcode = NF_PUT_VARA_REAL (ncid, aodwfid, ndimstart4, ndimcount4, &
             real(GC_aod_Jacobians(1:nwav,1:GC_nlayers,iout,1:ngeom,didx), kind=4))
        if (.not. do_aer_columnwf .and. do_assa_Jacobians) &
             rcode = NF_PUT_VARA_REAL (ncid, assawfid, ndimstart4, ndimcount4, &
             real(GC_assa_Jacobians(1:nwav,1:GC_nlayers,iout,1:ngeom,didx), kind=4))
        if (.not. do_cld_columnwf .and. do_cod_Jacobians) &
             rcode = NF_PUT_VARA_REAL (ncid, codwfid, ndimstart4, ndimcount4, &
             real(GC_cod_Jacobians(1:nwav,1:GC_nlayers,iout,1:ngeom,didx), kind=4))
        if (.not. do_cld_columnwf .and. do_cssa_Jacobians) &
             rcode = NF_PUT_VARA_REAL (ncid, cssawfid, ndimstart4, ndimcount4, &
             real(GC_cssa_Jacobians(1:nwav,1:GC_nlayers,iout,1:ngeom,didx), kind=4))   
     endif
     
     if (do_QU_Jacobians) then
        if (.not. do_aer_columnwf .and. do_aod_Jacobians) &
             rcode = NF_PUT_VARA_REAL (ncid, aodqwfid, ndimstart4, ndimcount4, &
             real(GC_aod_QJacobians(1:nwav,1:GC_nlayers,iout,1:ngeom,didx), kind=4))
        if (.not. do_aer_columnwf .and. do_assa_Jacobians) &
             rcode = NF_PUT_VARA_REAL (ncid, assaqwfid, ndimstart4, ndimcount4, &
             real(GC_assa_QJacobians(1:nwav,1:GC_nlayers,iout,1:ngeom,didx), kind=4))
        if (.not. do_cld_columnwf .and. do_cod_Jacobians) &
             rcode = NF_PUT_VARA_REAL (ncid, codqwfid, ndimstart4, ndimcount4, &
             real(GC_cod_QJacobians(1:nwav,1:GC_nlayers,iout,1:ngeom,didx), kind=4))
        if (.not. do_cld_columnwf .and. do_cssa_Jacobians) &
             rcode = NF_PUT_VARA_REAL (ncid, cssaqwfid, ndimstart4, ndimcount4, &
             real(GC_cssa_QJacobians(1:nwav,1:GC_nlayers,iout,1:ngeom,didx), kind=4)) 
        
        if (.not. do_aer_columnwf .and. do_aod_Jacobians) &
             rcode = NF_PUT_VARA_REAL (ncid, aoduwfid, ndimstart4, ndimcount4, &
             real(GC_aod_UJacobians(1:nwav,1:GC_nlayers,iout,1:ngeom,didx), kind=4))
        if (.not. do_aer_columnwf .and. do_assa_Jacobians) &
             rcode = NF_PUT_VARA_REAL (ncid, assauwfid, ndimstart4, ndimcount4, &
             real(GC_assa_UJacobians(1:nwav,1:GC_nlayers,iout,1:ngeom,didx), kind=4))
        if (.not. do_cld_columnwf .and. do_cod_Jacobians) &
             rcode = NF_PUT_VARA_REAL (ncid, coduwfid, ndimstart4, ndimcount4, &
             real(GC_cod_UJacobians(1:nwav,1:GC_nlayers,iout,1:ngeom,didx), kind=4))
        if (.not. do_cld_columnwf .and. do_cssa_Jacobians) &
             rcode = NF_PUT_VARA_REAL (ncid, cssauwfid, ndimstart4, ndimcount4, &
             real(GC_cssa_UJacobians(1:nwav,1:GC_nlayers,iout,1:ngeom,didx), kind=4))    
     endif
     
     ! 4D variables, wavdim, geodim, gasdim, olvdim
     if (do_AMF_calculation) then
        ndimstart4 = (/ 1, 1, 1, iout/)
        ndimcount4 = (/ nwav, ngeom, ngases, 1 /) 
        rcode = NF_PUT_VARA_REAL (ncid, amfid, ndimstart4, ndimcount4, &
             real(GC_AMFs(1:nwav,iout,1:ngeom,1:ngases,didx), kind=4))
     endif
     
     ! 5D variables, wavdim, laydim, geodim, gasdim, olvdim
     ndimstart5 = (/ 1,          1,     1,      1, iout /)
     ndimcount5 = (/ nwav, GC_nlayers, ngeom, ngases,    1 /) 
     if (do_Jacobians) then
        rcode = NF_PUT_VARA_REAL (ncid, gaswfid, ndimstart5, ndimcount5, &
             real(GC_Tracegas_Jacobians(1:nwav,1:GC_nlayers,iout,1:ngeom,1:ngases,didx), kind=4))  
        if (do_QU_Jacobians) then
           rcode = NF_PUT_VARA_REAL (ncid, gasqwfid, ndimstart5, ndimcount5, &
                real(GC_Tracegas_QJacobians(1:nwav,1:GC_nlayers,iout,1:ngeom,1:ngases,didx), kind=4))  
           rcode = NF_PUT_VARA_REAL (ncid, gasuwfid, ndimstart5, ndimcount5, &
                real(GC_Tracegas_UJacobians(1:nwav,1:GC_nlayers,iout,1:ngeom,1:ngases,didx), kind=4))  
        endif
     endif
  end do ! End output levels loop
  
  !==============================================================================
  ! CLOSE the NetCDF file
  !==============================================================================
  rcode = NF_CLOSE(ncid)
     
  return
end subroutine netcdf_wrt
