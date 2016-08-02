!     file ew_sub.f
!     compile with: f2py -c -m ew_sub ew_sub.f
!     This file contains subroutines for generating mock angles on the sky

      function get_mock_angles (costhetas, bal_flags, threshold,&
                                 &max_angle, fluxes, flux_limit, npts, nfluxes)
      
!     subroutine get_mock_angles(costhetas, bal_flags, npts)
      integer i
      logical undetected 
      real(8) PI, costheta, theta
      parameter (PI = 3.14159265359)

!  	declare the arrays we will actually modify and the 
!     variables supplied by the python frontend        
      real(8), intent(inout), dimension(0:npts-1) :: costhetas
      integer, intent(inout), dimension(0:npts-1) :: bal_flags
      real(8), intent(in), dimension(0:nfluxes-1) :: fluxes
      integer, intent(in) :: npts, nfluxes 
      real(8), intent(in) :: threshold, max_angle, flux_limit 


!     loop over npts, the number of mock quasars in our sample
      do 10 i = 1, npts

        undetected = .TRUE.

!       apply the detection limit
        do while (undetected)

          call random_number(costheta)
          theta = dacos(costheta) * 180.d0 / PI

!          write(*,*) costheta, theta, max_angle

          call is_source_detected(costheta, undetected,&
                                  &fluxes, flux_limit, nfluxes)

          if (theta .gt. max_angle) then
            undetected = .TRUE.
          endif

        enddo

        costhetas(i) = costheta

        if (theta .gt. threshold) then
          bal_flags(i) = 1
        else
          bal_flags(i) = 0
        endif

10    enddo
      end





      subroutine is_source_detected(costheta, undetected, fluxes,&
                                    &flux_limit, nfluxes)
      real(8) costheta, flux_limit, flux, z
      integer nfluxes, i
      real(8), intent(in), dimension(0:nfluxes-1) :: fluxes
      logical undetected

      call random_number(z)
      i = FLOOR(nfluxes*z)
      flux = fluxes(i) * costheta

      if (flux .gt. flux_limit) then
        undetected = .FALSE.
      endif

      return
      end