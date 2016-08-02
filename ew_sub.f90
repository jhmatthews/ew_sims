!     file ew_sub.f
!     compile with: f2py -c -m ew_sub ew_sub.f

! -------------------------------------------------
!     function get_mock_angles
!     This function generates mock angles for quasars. 
!     Arguments:
!       costhetas     array-like, float
!                     the array that holds the cosines of the angles 
!       bal_flags     array-like, int32
!                     1 for a BALQSO, 0 for a quasar
!       thetamin      float. thetamin
!       max_angle     float. thetamax
!       fluxes        array-like, float. fluxes of quasars.
!       flux_limit    float. Single value flux limit for rejection limit.
!       npts, nfluxes ints. Lengths of costhetas and fluxes respectively.
! 
!     note: be careful with placing the integer array lengths, as f2py will 
!     move them to the end eithout telling the user.
! -------------------------------------------------

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

          ! random number in costheta space       
          call random_number(costheta)

          ! convert to angle in degrees
          theta = dacos(costheta) * 180.d0 / PI

          ! does this mock object survive the flux limit?
          call is_source_detected(costheta, undetected,&
                                  &fluxes, flux_limit, nfluxes)

          ! is the object "obscured"?
          if (theta .gt. max_angle) then
            undetected = .TRUE.
          endif

        enddo

        ! if we get here then we've survived the various tests
        ! so copy the value over to the array
        costhetas(i) = costheta

        ! flag if this object is a BALQSO or not
        if (theta .gt. threshold) then
          bal_flags(i) = 1
        else
          bal_flags(i) = 0
        endif

10    enddo
      return
      end

! -------------------------------------------------
!     subroutine is_source_detected
!     This function selects a random flux and imposes a flux rejection limit.
!     It modifies the logical variable undetected to reflect this.
!     Arguments:
!       costheta      float. cosine of angle in question
!       undetected    bool. should be false on entry 
!                     and is modified by this routine
!       fluxes        array-like, float. fluxes of quasars.
!       flux_limit    float. Single value flux limit for rejection limit.
!       nfluxes       int.  Length of fluxes.
! -------------------------------------------------

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