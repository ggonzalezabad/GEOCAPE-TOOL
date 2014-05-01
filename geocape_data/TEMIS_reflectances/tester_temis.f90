program tester

   integer, parameter :: max_lambdas = 1000
   real(kind=8)       :: lat, lon
   integer            :: month, n
   integer            :: n_lambdas
   real(kind=8)       :: lambdas ( max_lambdas )
   real(kind=8)       :: ground_ler ( max_lambdas )

   lat = 40.1
   lon = -90.7
   month = 7
   n_lambdas  = 20
   do n = 1, n_lambdas
     lambdas(n) = 366.0d0 + 30*(n-1)
   enddo

   call get_temis_lers &
     ( lat, lon, month, max_lambdas, n_lambdas, lambdas, ground_ler )

   stop
end program tester
