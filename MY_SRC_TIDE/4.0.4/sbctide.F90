MODULE sbctide
   !!======================================================================
   !!                       ***  MODULE  sbctide  ***
   !! Initialization of tidal forcing
   !!======================================================================
   !! History :  9.0  !  2007  (O. Le Galloudec)  Original code
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers variables
   USE dom_oce        ! ocean space and time domain
   USE phycst         ! physical constant
   USE daymod         ! calandar
   USE tideini        ! 
   !
   USE in_out_manager ! I/O units
   USE iom            ! xIOs server
   USE ioipsl         ! NetCDF IPSL library
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   ! NB - to access love number 
   USE bdytides
   ! END NB

   IMPLICIT NONE
   PUBLIC

   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   pot_astro   !

   !!----------------------------------------------------------------------
   !!   tidal potential
   !!----------------------------------------------------------------------
   !!   sbc_tide            : 
   !!   tide_init_potential :
   !!----------------------------------------------------------------------

   REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) ::   amp_pot, phi_pot
   REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) ::   amp_load, phi_load

   ! davbyr: Tide drag
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   tdiss, h2rough,N2mean

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: sbctide.F90 10068 2018-08-28 14:09:04Z nicolasmartin $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE sbc_tide( kt )
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE sbc_tide  ***
      !!----------------------------------------------------------------------      
      INTEGER, INTENT( in ) ::   kt     ! ocean time-step
      INTEGER               ::   jk     ! dummy loop index
      INTEGER               ::   nsec_day_orig     ! Temporary variable
      !!----------------------------------------------------------------------
      
      IF( nsec_day == NINT(0.5_wp * rdt) .OR. kt == nit000 ) THEN      ! start a new day
         !
         IF( kt == nit000 )THEN
            ALLOCATE( amp_pot(jpi,jpj,nb_harmo),                      &
               &      phi_pot(jpi,jpj,nb_harmo), pot_astro(jpi,jpj)   )
            IF( ln_read_load )THEN
               ALLOCATE( amp_load(jpi,jpj,nb_harmo), phi_load(jpi,jpj,nb_harmo) )
               CALL tide_init_load
            ENDIF
            IF( ln_int_wave_drag )THEN !davbyr: Allocation and read tdiss
               ALLOCATE( tdiss(jpi, jpj) ) 
               tdiss(:,:) = 0.0_wp       
               IF ( ln_calc_tdiss ) THEN
                ALLOCATE( h2rough(jpi, jpj) )
                ALLOCATE( N2mean(jpi, jpj) )
                h2rough(:,:)     = 0.0_wp
                N2mean(:, :) = 0.0_wp

               ENDIF
               CALL tide_init_diss
            ENDIF
         ENDIF
         !
         IF( ln_read_load )THEN
            amp_pot(:,:,:) = amp_load(:,:,:)
            phi_pot(:,:,:) = phi_load(:,:,:)
         ELSE 
            amp_pot(:,:,:) = 0._wp
            phi_pot(:,:,:) = 0._wp
         ENDIF
         pot_astro(:,:) = 0._wp
         !
         ! If the run does not start from midnight then need to initialise tides
         ! at the start of the current day (only occurs when kt==nit000)
         ! Temporarily set nsec_day to beginning of day.
         nsec_day_orig = nsec_day
         IF ( nsec_day /= NINT(0.5_wp * rdt) ) THEN 
            kt_tide = kt - (nsec_day - 0.5_wp * rdt)/rdt
            nsec_day = NINT(0.5_wp * rdt)
         ELSE
            kt_tide = kt 
         ENDIF
         CALL tide_harmo( omega_tide, v0tide, utide, ftide, ntide, nb_harmo )
         !
         !
         IF(lwp) THEN
            WRITE(numout,*)
            WRITE(numout,*) 'sbc_tide : Update of the components and (re)Init. the potential at kt=', kt
            WRITE(numout,*) '~~~~~~~~ '
            DO jk = 1, nb_harmo
               WRITE(numout,*) Wave(ntide(jk))%cname_tide, utide(jk), ftide(jk), v0tide(jk), omega_tide(jk)
            END DO
         ENDIF
         !
         IF( ln_tide_pot )   CALL tide_init_potential
         !
         ! Reset nsec_day
         nsec_day = nsec_day_orig 
      ENDIF
      !
   END SUBROUTINE sbc_tide


   SUBROUTINE tide_init_potential
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE tide_init_potential  ***
      !!----------------------------------------------------------------------
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      REAL(wp) ::   zcons, ztmp1, ztmp2, zlat, zlon, ztmp, zamp, zcs   ! local scalar
      !!----------------------------------------------------------------------

      DO jk = 1, nb_harmo
         !--- NB 11/2017
         ! love number now provides in tide namelist
         zcons = dn_love_number * Wave(ntide(jk))%equitide * ftide(jk)
         ! ORIGINAL zcons = 0.7_wp * Wave(ntide(jk))%equitide * ftide(jk)
         !--- END NB
         DO ji = 1, jpi
            DO jj = 1, jpj
               ztmp1 =  ftide(jk) * amp_pot(ji,jj,jk) * COS( phi_pot(ji,jj,jk) + v0tide(jk) + utide(jk) )
               ztmp2 = -ftide(jk) * amp_pot(ji,jj,jk) * SIN( phi_pot(ji,jj,jk) + v0tide(jk) + utide(jk) )
               zlat = gphit(ji,jj)*rad !! latitude en radian
               zlon = glamt(ji,jj)*rad !! longitude en radian
               ztmp = v0tide(jk) + utide(jk) + Wave(ntide(jk))%nutide * zlon
               ! le potentiel est composé des effets des astres:
               IF    ( Wave(ntide(jk))%nutide == 1 )  THEN  ;  zcs = zcons * SIN( 2._wp*zlat )
               ELSEIF( Wave(ntide(jk))%nutide == 2 )  THEN  ;  zcs = zcons * COS( zlat )**2
               !--- NB 11/2017
               ! Add tide potential for long period tides
               ELSEIF( Wave(ntide(jk))%nutide == 0 )  THEN  ;  zcs = zcons * (0.5_wp-1.5_wp*SIN(zlat)**2._wp)
               !--- END NB
               ELSE                                         ;  zcs = 0._wp
               ENDIF
               ztmp1 = ztmp1 + zcs * COS( ztmp )
               ztmp2 = ztmp2 - zcs * SIN( ztmp )
               zamp = SQRT( ztmp1*ztmp1 + ztmp2*ztmp2 )
               amp_pot(ji,jj,jk) = zamp
               phi_pot(ji,jj,jk) = ATAN2( -ztmp2 / MAX( 1.e-10_wp , zamp ) ,   &
                  &                        ztmp1 / MAX( 1.e-10_wp,  zamp )   )
            END DO
         END DO
      END DO
      !
   END SUBROUTINE tide_init_potential

   SUBROUTINE tide_init_load
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE tide_init_load  ***
      !!----------------------------------------------------------------------
      INTEGER :: inum                 ! Logical unit of input file
      INTEGER :: ji, jj, itide        ! dummy loop indices
      REAL(wp), DIMENSION(jpi,jpj) ::   ztr, zti   !: workspace to read in tidal harmonics data 
      !!----------------------------------------------------------------------
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'tide_init_load : Initialization of load potential from file'
         WRITE(numout,*) '~~~~~~~~~~~~~~ '
      ENDIF
      !
      CALL iom_open ( cn_tide_load , inum )
      !
      DO itide = 1, nb_harmo
         CALL iom_get  ( inum, jpdom_data,TRIM(Wave(ntide(itide))%cname_tide)//'_z1', ztr(:,:) )
         CALL iom_get  ( inum, jpdom_data,TRIM(Wave(ntide(itide))%cname_tide)//'_z2', zti(:,:) )
         !
         DO ji=1,jpi
            DO jj=1,jpj
               amp_load(ji,jj,itide) =  SQRT( ztr(ji,jj)**2. + zti(ji,jj)**2. )
               phi_load(ji,jj,itide) = ATAN2(-zti(ji,jj), ztr(ji,jj) )
            END DO
         END DO
         !
      END DO
      CALL iom_close( inum )
      !
   END SUBROUTINE tide_init_load

   SUBROUTINE tide_init_diss !davbyr: subroutine
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE tide_init_diss  ***
      !!----------------------------------------------------------------------
      INTEGER :: inum                 ! Logical unit of input file
      INTEGER :: ji, jj,jidbg,jjdbg    ! dummy loop indices

      !!----------------------------------------------------------------------
      IF (ln_calc_tdiss) THEN
        IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'tide_init_diss : Read roughness (h) from file:',cn_h_rough
         WRITE(numout,*) '~~~~~~~~~~~~~~ '
        ENDIF

        CALL iom_open(cn_h_rough, inum)
        CALL iom_get (inum, jpdom_data, 'h', h2rough(:,:))
        CALL iom_close(inum)
        CALL iom_close( inum )
        ! mask hrough in shllow water
        DO ji=2,jpim1
         DO jj=2,jpjm1
          h2rough(ji,jj) = h2rough(ji,jj) * h2rough(ji,jj) ! read in h, need h^2
          IF ( gdepw_0(ji,jj,mbkt(ji,jj)+1) < tdiss_mindepth ) THEN
           h2rough(ji,jj) = 0.0
          ENDIF
         ENDDO
        ENDDO
        CALL lbc_lnk( 'tide_init_diss', h2rough, 'T', 1._wp )

        !!!!!!!
      ELSE
        IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'tide_init_diss : Read tidal dissipation from file:',cn_int_wave_drag
         WRITE(numout,*) '~~~~~~~~~~~~~~ '
        ENDIF
         
        CALL iom_open(cn_int_wave_drag, inum)
        CALL iom_get (inum, jpdom_data, 'tdiss', tdiss(:,:))
        CALL iom_close(inum)
        CALL iom_close( inum )

      ENDIF

   END SUBROUTINE tide_init_diss

  !!======================================================================
END MODULE sbctide
