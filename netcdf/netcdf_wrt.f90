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
  integer :: ncid, rcode, iwav, nwav, iout, ngeom
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
  if (do_vector_calculation .and. do_StokesQU_output) then
     qid      = ncvid(ncid, 'q', rcode)
     uid      = ncvid(ncid, 'u', rcode)
     qfluxid  = ncvid(ncid, 'qflux', rcode)
     ufluxid  = ncvid(ncid, 'uflux', rcode)
     qdfluxid = ncvid(ncid, 'qdirect_flux', rcode)
     udfluxid = ncvid(ncid, 'udirect_flux', rcode)
  endif

  if (do_Jacobians) then

     if (use_lambertian) then
        sfcwfid = ncvid(ncid, 'surfalb_jac', rcode)
     else 
        wswfid  = ncvid(ncid, 'windspeed_jac', rcode)
     endif

     if (do_cfrac_Jacobians) &
          cfracwfid  = ncvid(ncid, 'cfrac_jac', rcode)
     if (do_sfcprs_Jacobians) &
          sfcprswfid = ncvid(ncid, 'sfcprs_jac', rcode)

  endif

  if (do_QU_Jacobians) then

     if (use_lambertian) then
        sfcqwfid = ncvid(ncid, 'surfalb_qjac', rcode)
        sfcuwfid = ncvid(ncid, 'surfalb_ujac', rcode)
     else 
        wsqwfid = ncvid(ncid, 'windspeed_qjac', rcode)
        wsuwfid = ncvid(ncid, 'windspeed_ujac', rcode)
     endif

     if (do_cfrac_Jacobians) then
          cfracqwfid = ncvid(ncid, 'cfrac_qjac', rcode)
          cfracuwfid = ncvid(ncid, 'cfrac_ujac', rcode)
     endif
     if (do_sfcprs_Jacobians) then
          sfcprsqwfid = ncvid(ncid, 'sfcprs_qjac', rcode)
          sfcprsuwfid = ncvid(ncid, 'sfcprs_ujac', rcode)
     endif

  endif

  ! variables with 4D, wavdim, laydim, geodim, olvdim
  if (do_AMF_calculation) then
     scatwtid = ncvid(ncid, 'scatweights', rcode)
  endif

  if (do_Jacobians) then
     if (do_T_Jacobians) &
          tempwfid = ncvid(ncid, 't_jac', rcode)
     if (.not. do_aer_columnwf .and. do_aod_Jacobians) &
          aodwfid  = ncvid(ncid, 'aod_jac', rcode)  
     if (.not. do_aer_columnwf .and. do_assa_Jacobians) &
          assawfid = ncvid(ncid, 'assa_jac', rcode)   
     if (.not. do_cld_columnwf .and. do_cod_Jacobians) &
          codwfid  = ncvid(ncid, 'cod_jac', rcode)  
     if (.not. do_cld_columnwf .and. do_cssa_Jacobians) &
          cssawfid = ncvid(ncid, 'cssa_jac', rcode)  
  endif

  if (do_QU_Jacobians) then
     if (.not. do_aer_columnwf .and. do_aod_Jacobians) then 
          aodqwfid = ncvid(ncid, 'aod_qjac', rcode)  
          aoduwfid = ncvid(ncid, 'aod_ujac', rcode)  
     endif
     if (.not. do_aer_columnwf .and. do_assa_Jacobians) then
          assaqwfid = ncvid(ncid, 'assa_qjac', rcode)   
          assauwfid = ncvid(ncid, 'assa_ujac', rcode)   
     endif
     if (.not. do_cld_columnwf .and. do_cod_Jacobians) then
          codqwfid = ncvid(ncid, 'cod_qjac', rcode)  
          coduwfid = ncvid(ncid, 'cod_ujac', rcode)  
     endif
     if (.not. do_cld_columnwf .and. do_cssa_Jacobians) then
          cssaqwfid = ncvid(ncid, 'cssa_qjac', rcode)  
          cssauwfid = ncvid(ncid, 'cssa_ujac', rcode) 
     endif
  endif

  ! variables with 4D, wavdim, geodim, gasdim, olvdim
  if (do_AMF_calculation) then
     amfid = ncvid(ncid, 'amf', rcode)
  endif

  ! variables with 5D, wavdim, laydim, geodim, gasdim, ovdim
  if (do_Jacobians) then
     gaswfid = ncvid(ncid,  'gas_jac', rcode)
     if (do_QU_Jacobians) then
        gasqwfid = ncvid(ncid,  'gas_qjac', rcode)
        gasuwfid = ncvid(ncid,  'gas_ujac', rcode)
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
!!$  DO iwav = 1, nwav
!!$
     print*, 'Writing output #', didx

     ! write 1D variables with wavdim
     ndimstart1 = (/ 1 /)
     ndimcount1 = (/ nwav /)  
     call ncvpt (ncid, irradid, ndimstart1, ndimcount1, real(solar_cspec_data(1:nwav), kind=4), rcode)
     call ncvpt (ncid, sfcid,   ndimstart1, ndimcount1, real(ground_ler(1:nwav), kind=4),    rcode)

     ! 2D, variables, wavdim, laydim
     ndimstart2 = (/    1,  1 /)
     ndimcount2 = (/ nwav, GC_nlayers /) 
     call ncvpt (ncid, odsid,   ndimstart2, ndimcount2, real(opdeps(1:nwav,1:GC_nlayers), kind=4), rcode)
     call ncvpt (ncid, ssasid,  ndimstart2, ndimcount2, real(ssalbs(1:nwav,1:GC_nlayers), kind=4), rcode)
     call ncvpt (ncid, aodsid,  ndimstart2, ndimcount2, real(aer_opdeps(1:nwav,1:GC_nlayers), kind=4), rcode)
     call ncvpt (ncid, assasid, ndimstart2, ndimcount2, real(aer_ssalbs(1:nwav,1:GC_nlayers), kind=4), rcode)
     call ncvpt (ncid, codsid,  ndimstart2, ndimcount2, real(cld_opdeps(1:nwav,1:GC_nlayers), kind=4), rcode)
     call ncvpt (ncid, cssasid, ndimstart2, ndimcount2, real(cld_ssalbs(1:nwav,1:GC_nlayers), kind=4), rcode)

     ! write 3D variables with wavdim,geodim,olvdim
     DO iout = 1, GC_n_user_levels
        ndimstart3 = (/ 1, 1 , iout/)
        ndimcount3 = (/ nwav, ngeom , 1/)  
        call ncvpt (ncid, radid,   ndimstart3, ndimcount3, &
             real(GC_radiances(1:nwav,iout,1:ngeom,didx), kind=4), rcode)
        ndimcount3 = (/   1, GC_n_sun_positions , 1/)  
        call ncvpt (ncid, fluxid,  ndimstart3, ndimcount3, &
             real(GC_flux(1:nwav,iout,1:GC_n_sun_positions,didx), kind=4), rcode)
        call ncvpt (ncid, dfluxid, ndimstart3, ndimcount3, &
             real(GC_direct_flux(1:nwav,iout,1:GC_n_sun_positions,didx), kind=4), rcode)
        if (do_vector_calculation .and. do_StokesQU_output) then
           ndimcount3 = (/ 1 , ngeom , 1/)
           call ncvpt (ncid, qid, ndimstart3, ndimcount3, &
                real(GC_Qvalues(1:nwav,iout,1:ngeom,didx), kind=4), rcode)
           call ncvpt (ncid, uid, ndimstart3, ndimcount3, &
                real(GC_Uvalues(1:nwav,iout,1:ngeom,didx), kind=4), rcode)
           ndimcount3 = (/   1, GC_n_sun_positions , 1/)  
           call ncvpt (ncid, qfluxid,  ndimstart3, ndimcount3, &
                real(GC_Qflux(1:nwav,iout,1:GC_n_sun_positions,didx), kind=4), rcode)
           call ncvpt (ncid, ufluxid,  ndimstart3, ndimcount3, &
                real(GC_Uflux(1:nwav,iout,1:GC_n_sun_positions,didx), kind=4), rcode)
           call ncvpt (ncid, qdfluxid, ndimstart3, ndimcount3, &
                real(GC_Qdirect_flux(1:nwav,iout,1:GC_n_sun_positions,didx), kind=4), rcode)
           call ncvpt (ncid, udfluxid, ndimstart3, ndimcount3, &
                real(GC_Udirect_flux(1:nwav,iout,1:GC_n_sun_positions,didx), kind=4), rcode)
        endif
        
        ndimstart3 = (/ 1, 1 , iout/)
        ndimcount3 = (/ nwav, ngeom , 1/)  
        if (do_Jacobians) then
           if (use_lambertian) then
              call ncvpt (ncid, sfcwfid, ndimstart3, ndimcount3, &
                   real(GC_Surfalbedo_Jacobians(1:nwav,iout,1:ngeom,didx), kind=4), rcode)
           else 
              call ncvpt (ncid, wswfid, ndimstart3, ndimcount3, &
                   real(GC_Windspeed_Jacobians(1:nwav,iout,1:ngeom,didx), kind=4), rcode)
           endif
           if (do_cfrac_Jacobians) &
                call ncvpt (ncid, cfracwfid,  ndimstart3, ndimcount3, &
                real(GC_cfrac_Jacobians(1:nwav,iout,1:ngeom,didx), kind=4), rcode)
           if (do_sfcprs_Jacobians) &
                call ncvpt (ncid, sfcprswfid, ndimstart3, ndimcount3, &
                real(GC_sfcprs_Jacobians(1:nwav,iout,1:ngeom,didx), kind=4), rcode)
        endif
        
        if (do_QU_Jacobians) then
           if (use_lambertian) then
              call ncvpt (ncid, sfcqwfid, ndimstart3, ndimcount3, &
                   real(GC_Surfalbedo_QJacobians(1:nwav,iout,1:ngeom,didx), kind=4), rcode)
           else 
              call ncvpt (ncid, wsqwfid,  ndimstart3, ndimcount3, &
                   real(GC_Windspeed_QJacobians(1:nwav,iout,1:ngeom,didx), kind=4), rcode)
           endif
           if (do_cfrac_Jacobians) &
                call ncvpt (ncid, cfracqwfid, ndimstart3, ndimcount3, &
                real(GC_cfrac_QJacobians(1:nwav,iout,1:ngeom,didx), kind=4), rcode)
           if (do_sfcprs_Jacobians) &
                call ncvpt (ncid, sfcprsqwfid, ndimstart3, ndimcount3, &
                real(GC_sfcprs_QJacobians(1:nwav,iout,1:ngeom,didx), kind=4), rcode)
           if (use_lambertian) then
              call ncvpt (ncid, sfcuwfid, ndimstart3, ndimcount3, &
                   real(GC_Surfalbedo_UJacobians(1:nwav,iout,1:ngeom,didx), kind=4), rcode)
           else 
              call ncvpt (ncid, wsuwfid, ndimstart3, ndimcount3, &
                   real(GC_Windspeed_UJacobians(1:nwav,iout,1:ngeom,didx), kind=4), rcode)
           endif
           if (do_cfrac_Jacobians) &
                call ncvpt (ncid, cfracuwfid, ndimstart3, ndimcount3, &
                real(GC_cfrac_UJacobians(1:nwav,iout,1:ngeom,didx), kind=4), rcode)
           if (do_sfcprs_Jacobians) &
                call ncvpt (ncid, sfcprsuwfid, ndimstart3, ndimcount3, &
                real(GC_sfcprs_UJacobians(1:nwav,iout,1:ngeom,didx), kind=4), rcode)
        endif
        
        ! 4D variables, wavdim, altdim, geodim, olvdim
        ndimstart4 = (/    1,          1,     1, iout/)
        ndimcount4 = (/ nwav, GC_nlayers, ngeom,    1/) 
        if (do_AMF_calculation) then
           call ncvpt (ncid, scatwtid, ndimstart4, ndimcount4, &
                real(GC_Scattering_Weights(1:nwav,1:GC_nlayers,iout,1:ngeom,didx), kind=4), rcode)
        endif
        
        if (do_Jacobians) then
           if (do_T_Jacobians) &
                call ncvpt (ncid, tempwfid, ndimstart4, ndimcount4, &
                real(GC_Temperature_Jacobians(1:nwav,1:GC_nlayers,iout,1:ngeom,didx), kind=4), rcode)
           if (.not. do_aer_columnwf .and. do_aod_Jacobians) &
                call ncvpt (ncid, aodwfid, ndimstart4, ndimcount4, &
                real(GC_aod_Jacobians(1:nwav,1:GC_nlayers,iout,1:ngeom,didx), kind=4), rcode)
           if (.not. do_aer_columnwf .and. do_assa_Jacobians) &
                call ncvpt (ncid, assawfid, ndimstart4, ndimcount4, &
                real(GC_assa_Jacobians(1:nwav,1:GC_nlayers,iout,1:ngeom,didx), kind=4), rcode)
           if (.not. do_cld_columnwf .and. do_cod_Jacobians) &
                call ncvpt (ncid, codwfid, ndimstart4, ndimcount4, &
                real(GC_cod_Jacobians(1:nwav,1:GC_nlayers,iout,1:ngeom,didx), kind=4), rcode)
           if (.not. do_cld_columnwf .and. do_cssa_Jacobians) &
                call ncvpt (ncid, cssawfid, ndimstart4, ndimcount4, &
                real(GC_cssa_Jacobians(1:nwav,1:GC_nlayers,iout,1:ngeom,didx), kind=4), rcode)   
        endif
        
        if (do_QU_Jacobians) then
           if (.not. do_aer_columnwf .and. do_aod_Jacobians) &
                call ncvpt (ncid, aodqwfid, ndimstart4, ndimcount4, &
                real(GC_aod_QJacobians(1:nwav,1:GC_nlayers,iout,1:ngeom,didx), kind=4), rcode)
           if (.not. do_aer_columnwf .and. do_assa_Jacobians) &
                call ncvpt (ncid, assaqwfid, ndimstart4, ndimcount4, &
                real(GC_assa_QJacobians(1:nwav,1:GC_nlayers,iout,1:ngeom,didx), kind=4), rcode)
           if (.not. do_cld_columnwf .and. do_cod_Jacobians) &
                call ncvpt (ncid, codqwfid, ndimstart4, ndimcount4, &
                real(GC_cod_QJacobians(1:nwav,1:GC_nlayers,iout,1:ngeom,didx), kind=4), rcode)
           if (.not. do_cld_columnwf .and. do_cssa_Jacobians) &
                call ncvpt (ncid, cssaqwfid, ndimstart4, ndimcount4, &
                real(GC_cssa_QJacobians(1:nwav,1:GC_nlayers,iout,1:ngeom,didx), kind=4), rcode) 
           
           if (.not. do_aer_columnwf .and. do_aod_Jacobians) &
                call ncvpt (ncid, aoduwfid, ndimstart4, ndimcount4, &
                real(GC_aod_UJacobians(1:nwav,1:GC_nlayers,iout,1:ngeom,didx), kind=4), rcode)
           if (.not. do_aer_columnwf .and. do_assa_Jacobians) &
                call ncvpt (ncid, assauwfid, ndimstart4, ndimcount4, &
                real(GC_assa_UJacobians(1:nwav,1:GC_nlayers,iout,1:ngeom,didx), kind=4), rcode)
           if (.not. do_cld_columnwf .and. do_cod_Jacobians) &
                call ncvpt (ncid, coduwfid, ndimstart4, ndimcount4, &
                real(GC_cod_UJacobians(1:nwav,1:GC_nlayers,iout,1:ngeom,didx), kind=4), rcode)
           if (.not. do_cld_columnwf .and. do_cssa_Jacobians) &
                call ncvpt (ncid, cssauwfid, ndimstart4, ndimcount4, &
                real(GC_cssa_UJacobians(1:nwav,1:GC_nlayers,iout,1:ngeom,didx), kind=4), rcode)    
        endif
        
        ! 4D variables, wavdim, geodim, gasdim, olvdim
        if (do_AMF_calculation) then
           ndimstart4 = (/ 1, 1, 1, iout/)
           ndimcount4 = (/ nwav, ngeom, ngases, 1 /) 
           call ncvpt (ncid, amfid, ndimstart4, ndimcount4, &
                real(GC_AMFs(1:nwav,iout,1:ngeom,1:ngases,didx), kind=4), rcode)
        endif
        
        ! 5D variables, wavdim, laydim, geodim, gasdim, olvdim
        ndimstart5 = (/ 1,          1,     1,      1, iout /)
        ndimcount5 = (/ nwav, GC_nlayers, ngeom, ngases,    1 /) 
        if (do_Jacobians) then
           call ncvpt (ncid, gaswfid, ndimstart5, ndimcount5, &
                real(GC_Tracegas_Jacobians(1:nwav,1:GC_nlayers,iout,1:ngeom,1:ngases,didx), kind=4), rcode)  
           if (do_QU_Jacobians) then
              call ncvpt (ncid, gasqwfid, ndimstart5, ndimcount5, &
                   real(GC_Tracegas_QJacobians(1:nwav,1:GC_nlayers,iout,1:ngeom,1:ngases,didx), kind=4), rcode)  
              call ncvpt (ncid, gasuwfid, ndimstart5, ndimcount5, &
                   real(GC_Tracegas_UJacobians(1:nwav,1:GC_nlayers,iout,1:ngeom,1:ngases,didx), kind=4), rcode)  
           endif
        endif
     end do ! End output levels loop
!!$  end do ! End wavelength loop

  !==============================================================================
  ! CLOSE the NetCDF file
  !==============================================================================
  
  call ncclos (ncid, rcode)
     
  return
end subroutine netcdf_wrt
